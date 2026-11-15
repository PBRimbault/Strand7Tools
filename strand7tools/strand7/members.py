"""
Build design members from the model's beam elements. Read-only.

For each row of the member table:
  - collects its beams and puts them in order, end to end
  - checks the member is a single straight chain (no branches or gaps)
  - checks section/orientation consistency along the chain
  - works out the sign factors that convert each beam's forces into the member's axes

Only straight members of constant section are supported for now.
"""

import ctypes
import math
from dataclasses import dataclass, field

from ..member_table import MemberDef, parse_row, read_rows
from .api import S, check, get_str, total_beams


class MemberError(Exception):
    """A member cannot be built; the message says why."""


@dataclass
class BeamInChain:
    beam: int
    forward: bool        # True if the beam runs in the member direction (N1 -> N2)
    length: float
    x0: float            # chainage at the start of the beam
    x1: float
    prop: int
    axes_note: str
    # Sign factors (+1/-1) relating this beam's axes to the member's reference axes.
    # f3: beam axis 3 vs member direction; f1/f2: beam axis 1/2 vs the first beam's axes.
    # Used to convert Strand7 forces into the member's sign convention
    # (shear x f3*f_axis, moment x f_axis, axial x 1, torque x 1).
    f1: int
    f2: int
    f3: int


@dataclass
class Member:
    definition: MemberDef
    chain: list
    nodes: list
    chainage: list
    length: float
    section_name: str
    direction: tuple
    ref_axis1: tuple     # principal axis 1 of the first beam (global components)
    ref_axis2: tuple     # principal axis 2 of the first beam (global components)
    warnings: list = field(default_factory=list)

    @property
    def member_id(self):
        return self.definition.member_id

    def effective_lengths(self):
        """{'Lcr_y', 'Lcr_z', 'Lw'} in m (None where that check is switched off)."""
        return self.definition.effective_lengths(self.length)


# ----------------------------------------------------------------------
# Querying the model
# ----------------------------------------------------------------------

def load_beam_ids(uid, n_beams):
    ids = {}
    for b in range(1, n_beams + 1):
        v = ctypes.c_long()
        check(S.St7GetBeamID(uid, b, ctypes.byref(v)))
        ids[b] = v.value
    return ids


def beam_nodes(uid, beam):
    conn = (ctypes.c_long * (S.kMaxElementNode + 1))()
    check(S.St7GetElementConnection(uid, S.tyBEAM, beam, conn))
    if conn[0] not in (2, 3):  # Beam2, or Beam3 (3rd node is only the orientation node)
        raise MemberError(f"beam {beam} has {conn[0]} nodes; only Beam2/Beam3 are supported")
    return conn[1], conn[2]


def beam_axes(uid, beam):
    d = (ctypes.c_double * 9)()
    check(S.St7GetBeamAxisSystemInitial(uid, beam, d))
    return tuple(d[0:3]), tuple(d[3:6]), tuple(d[6:9])


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


# ----------------------------------------------------------------------
# Building a member
# ----------------------------------------------------------------------

def order_chain(ends):
    """
    Put beams in order from one end of the member to the other. ends = {beam: (node1, node2)}.
    Returns ([(beam, forward)], start node, end node). Pure logic, no Strand7.
    """
    adj = {}
    for b, (n1, n2) in ends.items():
        adj.setdefault(n1, []).append(b)
        adj.setdefault(n2, []).append(b)
    if any(len(bs) > 2 for bs in adj.values()):
        raise MemberError("the beams branch (3 or more beams meet at a node); "
                          "a member must be a single line of beams")
    terminals = sorted(n for n, bs in adj.items() if len(bs) == 1)
    if len(terminals) != 2:
        raise MemberError(f"the beams do not form one open chain "
                          f"(found {len(terminals)} free ends; there may be a gap or a closed loop)")
    start = terminals[0]
    order, used, node = [], set(), start
    while True:
        nxt = [b for b in adj[node] if b not in used]
        if not nxt:
            break
        b = nxt[0]
        used.add(b)
        n1, n2 = ends[b]
        forward = (n1 == node)
        order.append((b, forward))
        node = n2 if forward else n1
    if len(used) != len(ends):
        raise MemberError("the beams are not all connected to each other")
    return order, start, node


def build_member(uid, d):
    coords = {}

    def xyz(n):
        if n not in coords:
            a = (ctypes.c_double * 3)()
            check(S.St7GetNodeXYZ(uid, n, a))
            coords[n] = (a[0], a[1], a[2])
        return coords[n]

    ends = {b: beam_nodes(uid, b) for b in d.beams}
    order, start, finish = order_chain(ends)

    # Properties and axes
    props = {}
    axes = {}
    for b, _ in order:
        p = ctypes.c_long()
        check(S.St7GetElementProperty(uid, S.ptBEAMPROP, b, ctypes.byref(p)))
        props[b] = p.value
        axes[b] = beam_axes(uid, b)

    warnings = []
    section_names = {}
    for p in set(props.values()):
        section_names[p] = get_str(S.St7GetBeamSectionName, uid, p)
        ptype = ctypes.c_long()
        check(S.St7GetBeamPropertyType(uid, p, ctypes.byref(ptype)))
        if ptype.value != S.btBeam:
            warnings.append(f"property {p} is not a flexural beam type")
    if len(set(props.values())) > 1:
        warnings.append("beams use different properties: " +
                        ", ".join(f"{p} ({n})" for p, n in section_names.items()) +
                        " - only constant sections are supported for now")
    elif len(set(section_names.values())) > 1:
        warnings.append("different section names along the member")

    # Chain geometry
    p_start, p_finish = xyz(start), xyz(finish)
    u = tuple((q - p) / math.dist(p_start, p_finish) for p, q in zip(p_start, p_finish))
    nodes, chainage, chain = [start], [0.0], []
    ref_axis1, ref_axis2 = axes[order[0][0]][0], axes[order[0][0]][1]

    def sgn(x):
        return 1 if x >= 0 else -1

    for b, forward in order:
        n1, n2 = ends[b]
        length = math.dist(xyz(n1), xyz(n2))
        a1 = dot(axes[b][0], ref_axis1)
        a2 = dot(axes[b][1], ref_axis2)
        if abs(a1) >= 0.999 and abs(a2) >= 0.999:
            note = "same section axes"
        else:
            raise MemberError(f"beam {b}: section axes are rotated relative to the first beam "
                              f"(dot products {a1:.3f}, {a2:.3f}); a member needs a constant "
                              f"section orientation - check the beam Angle attribute")
        chain.append(BeamInChain(b, forward, length, chainage[-1], chainage[-1] + length,
                                 props[b], note, sgn(a1), sgn(a2), sgn(dot(axes[b][2], u))))
        nodes.append(n2 if forward else n1)
        chainage.append(chainage[-1] + length)
    total = chainage[-1]
    if total <= 0:
        raise MemberError("member has zero length")

    # Straightness
    p0 = p_start
    max_dev = 0.0
    for n in nodes:
        v = tuple(a - b for a, b in zip(xyz(n), p0))
        cross = (v[1] * u[2] - v[2] * u[1], v[2] * u[0] - v[0] * u[2], v[0] * u[1] - v[1] * u[0])
        max_dev = max(max_dev, math.sqrt(sum(c * c for c in cross)))
    if max_dev > 1e-4 * total:
        kink = 0.0
        for a, b, c in zip(nodes, nodes[1:], nodes[2:]):
            va = tuple(q - p for p, q in zip(xyz(a), xyz(b)))
            vb = tuple(q - p for p, q in zip(xyz(b), xyz(c)))
            cosang = dot(va, vb) / (math.sqrt(dot(va, va)) * math.sqrt(dot(vb, vb)))
            kink = max(kink, math.degrees(math.acos(max(-1.0, min(1.0, cosang)))))
        raise MemberError(f"member is not straight (nodes deviate up to {max_dev:.4g} "
                          f"from the line between its ends; the largest change of direction "
                          f"between adjacent beams is {kink:.1f} degrees). "
                          f"Only straight members are supported")

    return Member(d, chain, nodes, chainage, total,
                  next(iter(section_names.values())), u, ref_axis1, ref_axis2, warnings)


def load_members(uid, member_file):
    """
    Build every member in the member table. Returns (members, errors); errors are
    (member id, message) for rows that could not be built. A beam that is in more than one
    member gets a warning on the later member.
    """
    n_beams = total_beams(uid)
    rows = read_rows(member_file)
    need_ids = any("ID=" in r.get("beams", "").upper() for r in rows)
    beam_ids = load_beam_ids(uid, n_beams) if need_ids else {}
    members, errors, seen, owner = [], [], set(), {}
    for row in rows:
        member_id = row.get("member", "?")
        if member_id in seen:
            errors.append((member_id, "this member number is used more than once in the table; "
                                      "give each member a unique number"))
            continue
        seen.add(member_id)
        try:
            m = build_member(uid, parse_row(row, beam_ids, n_beams))
        except (ValueError, MemberError) as e:
            errors.append((member_id, str(e)))
            continue
        for b in m.definition.beams:
            if b in owner:
                m.warnings.append(f"beam {b} is also in member {owner[b]}")
            else:
                owner[b] = member_id
        members.append(m)
    return members, errors
