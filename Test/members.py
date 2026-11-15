"""
Step 3: Build design members from beam elements.  Read-only.

Terms:
  beam   - a Strand7 beam element
  member - a continuous, straight, constant-section chain of beams over which
           buckling can occur (from one restraint to the next)

Reads the member table (members.csv), then for each member:
  - collects its beams and puts them in order, end to end
  - checks the member is a single straight chain (no branches or gaps)
  - checks section/orientation consistency along the chain
  - works out the effective lengths from the k factors (or typed overrides)

Prints the result so you can confirm the structure before any checking is done.
Only straight members of constant section are supported for now.
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))  # st7_setup changes the working dir

import ctypes  # noqa: E402
import csv  # noqa: E402
import math  # noqa: E402
from dataclasses import dataclass, field  # noqa: E402

from st7_setup import St7API as S, check, model_session, MODEL_FILE  # noqa: E402

MEMBER_FILE = os.path.join(HERE, "members.csv")


class MemberError(Exception):
    """A member cannot be built; the message says why."""


# ----------------------------------------------------------------------
# Reading the member table
# ----------------------------------------------------------------------

YES = {"y", "yes", "true", "1", "x"}
NO = {"n", "no", "false", "0"}


def parse_bool(text, default=True):
    t = (text or "").strip().lower()
    if t == "":
        return default
    if t in YES:
        return True
    if t in NO:
        return False
    raise ValueError(f"cannot read {text!r} as yes/no")


def parse_float(text, default=None):
    t = (text or "").strip()
    return default if t == "" else float(t)


def parse_int_list(text):
    """'1-4;7' -> [1, 2, 3, 4, 7]  (';' or ',' separated)."""
    out = []
    for tok in (text or "").replace(",", ";").split(";"):
        tok = tok.strip()
        if not tok:
            continue
        if "-" in tok:
            a, b = (int(x) for x in tok.split("-", 1))
            if b < a:
                raise ValueError(f"bad range {tok!r}")
            out.extend(range(a, b + 1))
        else:
            out.append(int(tok))
    return out


def resolve_beams(text, beam_ids, n_beams):
    """Beam numbers from text like '1-10', '1;2;5' or 'ID=5' (all beams with that Beam ID)."""
    beams = []
    for tok in (text or "").replace(",", ";").split(";"):
        tok = tok.strip()
        if not tok:
            continue
        if tok.upper().startswith("ID="):
            wanted = int(tok[3:])
            found = [b for b, i in beam_ids.items() if i == wanted]
            if not found:
                raise ValueError(f"no beams have Beam ID {wanted}")
            beams.extend(found)
        else:
            beams.extend(parse_int_list(tok))
    unique = list(dict.fromkeys(beams))
    bad = [b for b in unique if not 1 <= b <= n_beams]
    if bad:
        raise ValueError(f"beam numbers not in the model: {bad}")
    return unique


@dataclass
class MemberDef:
    member_id: str
    beams: list
    buckling_y: bool = True
    kcr_y: float = 1.0
    lcr_y: float = None  # optional override of kcr_y x member length
    buckling_z: bool = True
    kcr_z: float = 1.0
    lcr_z: float = None
    ltb: bool = True
    kz: float = 1.0
    kw: float = 1.0
    lw: float = None     # optional override of kz x member length
    grade: str = ""      # steel grade, e.g. S355; blank = the default grade in settings.toml
    comment: str = ""


def read_rows(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Member table not found: {path}")
    with open(path, newline="", encoding="utf-8-sig") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    rows = []
    for raw in csv.DictReader(lines):
        rows.append({k.strip().lower(): (v or "").strip()
                     for k, v in raw.items() if k is not None})
    return rows


def parse_row(row, beam_ids, n_beams):
    member_id = row.get("member", "")
    if not member_id:
        raise ValueError("missing member number")
    beams = resolve_beams(row.get("beams", ""), beam_ids, n_beams)
    if not beams:
        raise ValueError("no beams listed")
    return MemberDef(
        member_id=member_id,
        beams=beams,
        buckling_y=parse_bool(row.get("buckling_y")),
        kcr_y=parse_float(row.get("kcr_y"), 1.0),
        lcr_y=parse_float(row.get("lcr_y")),
        buckling_z=parse_bool(row.get("buckling_z")),
        kcr_z=parse_float(row.get("kcr_z"), 1.0),
        lcr_z=parse_float(row.get("lcr_z")),
        ltb=parse_bool(row.get("ltb")),
        kz=parse_float(row.get("kz"), 1.0),
        kw=parse_float(row.get("kw"), 1.0),
        lw=parse_float(row.get("lw")),
        grade=row.get("grade", ""),
        comment=row.get("comment", ""),
    )


# ----------------------------------------------------------------------
# Querying the model
# ----------------------------------------------------------------------

def get_str(fn, *args):
    buf = ctypes.create_string_buffer(S.kMaxStrLen)
    check(fn(*args, buf, S.kMaxStrLen))
    return buf.value.decode("cp1252", errors="replace")  # Strand7 strings are Windows ANSI


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
    # (shear x f3*f_axis, moment x f_axis, axial x 1, torque x f3).
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


def order_chain(ends):
    """Put beams in order from one end of the member to the other."""
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


# ----------------------------------------------------------------------
# Effective lengths and reporting
# ----------------------------------------------------------------------

def effective_lengths(m):
    """Return {'Lcr_y', 'Lcr_z', 'Lw'} (None where that check is switched off)."""
    d = m.definition

    def pick(on, k, override):
        if not on:
            return None
        return override if override is not None else k * m.length

    return {
        "Lcr_y": pick(d.buckling_y, d.kcr_y, d.lcr_y),
        "Lcr_z": pick(d.buckling_z, d.kcr_z, d.lcr_z),
        "Lw": pick(d.ltb, d.kz, d.lw),
    }


def print_member(m):
    d = m.definition
    eff = effective_lengths(m)
    print("=" * 72)
    print(f"Member {d.member_id}" + (f"   ({d.comment})" if d.comment else ""))
    print(f"  Section : {m.section_name}")
    print(f"  Length  : {m.length:.4f}   direction ({m.direction[0]:.3f}, "
          f"{m.direction[1]:.3f}, {m.direction[2]:.3f})")
    a1, a2 = m.ref_axis1, m.ref_axis2
    print(f"  Axes    : axis 1 = ({a1[0]:.3f}, {a1[1]:.3f}, {a1[2]:.3f})  "
          f"axis 2 = ({a2[0]:.3f}, {a2[1]:.3f}, {a2[2]:.3f})   [beam {m.chain[0].beam}]")
    print("            major axis (y-y) = bending about axis 1 -> Moment 2, Shear 2, SA2 "
          "(axis 2 = web direction)")
    print("            minor axis (z-z) = bending about axis 2 -> Moment 1, Shear 1, SA1")
    print("  Beams   : " + " ".join(f"{c.beam}{'' if c.forward else '*'}" for c in m.chain)
          + "     (* = beam runs against the member direction)")
    print("  Sign factors vs member axes (+1 = same direction, -1 = opposite):")
    for c in m.chain:
        print(f"            beam {c.beam:>4}: axis3 {c.f3:+d}   axis1 {c.f1:+d}   axis2 {c.f2:+d}")
    print("  Nodes   : " + " ".join(str(n) for n in m.nodes))
    last = len(m.nodes) - 1
    print("            " + "--".join("X" if i in (0, last) else "o" for i in range(len(m.nodes))))
    print()
    print("  Effective lengths:")
    rows = (
        ("Flexural buckling about y", eff["Lcr_y"], d.kcr_y, d.lcr_y, "Lcr,y"),
        ("Flexural buckling about z", eff["Lcr_z"], d.kcr_z, d.lcr_z, "Lcr,z"),
        ("Lateral-torsional buckling", eff["Lw"], d.kz, d.lw, "Lw"),
    )
    for label, value, k, override, symbol in rows:
        if value is None:
            print(f"    {label:<27}: not checked")
        else:
            src = "typed value" if override is not None else f"k={k:g} x {m.length:.4f}"
            print(f"    {label:<27}: {symbol:<6} = {value:.4f}   ({src})")
    if m.warnings:
        print()
        for w in m.warnings:
            print(f"  WARNING: {w}")
    print()


def main():
    with model_session() as uid:
        total = ctypes.c_long()
        check(S.St7GetTotal(uid, S.tyBEAM, ctypes.byref(total)))
        n_beams = total.value

        rows = read_rows(MEMBER_FILE)
        if any(r.get("restraints") for r in rows):
            print("NOTE: the 'restraints' column is no longer used - members now end at "
                  "restraints. Split such members in the table.\n")
        need_ids = any("ID=" in r.get("beams", "").upper() for r in rows)
        beam_ids = load_beam_ids(uid, n_beams) if need_ids else {}

        print(f"Model: {MODEL_FILE}")
        print(f"Beams in model: {n_beams}")
        print(f"Member table: {MEMBER_FILE}\n")

        owner, built, seen_ids = {}, 0, set()
        for row in rows:
            member_id = row.get("member", "?")
            if member_id in seen_ids:
                print(f"Member {member_id}: ERROR - this member number is used more than "
                      f"once in the table; give each member a unique number\n")
                continue
            seen_ids.add(member_id)
            try:
                d = parse_row(row, beam_ids, n_beams)
                m = build_member(uid, d)
            except (ValueError, MemberError) as e:
                print(f"Member {member_id}: ERROR - {e}\n")
                continue
            for b in d.beams:
                if b in owner:
                    m.warnings.append(f"beam {b} is also in member {owner[b]}")
                else:
                    owner[b] = member_id
            print_member(m)
            built += 1

        print(f"{built} of {len(rows)} members built.")
        unassigned = [b for b in range(1, n_beams + 1) if b not in owner]
        if unassigned:
            shown = ", ".join(str(b) for b in unassigned[:30])
            more = " ..." if len(unassigned) > 30 else ""
            print(f"Beams not in any member ({len(unassigned)}): {shown}{more}")


if __name__ == "__main__":
    main()
