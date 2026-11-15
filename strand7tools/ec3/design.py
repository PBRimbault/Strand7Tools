"""
Design case evaluation: pure logic, no Strand7 and no printing.

  combine()          adds factored member force diagrams (linear superposition)
  evaluate_member()  runs the EC3 cross-section checks at every station of a member
  summarise()        governing station, class, utilisation and status for a member / design case
  buckling_check()   member buckling for a member / design case, from its station evaluations
  pick_station()     the station for a worked calculation
  status_rank()      sort key: failures first, then incomplete checks, then highest utilisation

A member force diagram ("by_beam") is the structure made by strand7.forces.assemble():
    [(BeamInChain, [(x, {"N":.., "Vy":.., "Vz":.., "T":.., "My":.., "Mz":..}), ...]), ...]
with forces in kN and kNm, in the member axes, tension positive.
"""

from dataclasses import dataclass, field

from ..quantities import QUANTITIES
from .buckling import check_member_buckling
from .cross_section import check_cross_section

CHECKS = ["N", "My", "Mz", "Vy", "Vz", "T", "N+M"]


class CombineError(Exception):
    """The diagrams being added do not line up."""


def combine(parts):
    """
    Add factored diagrams: parts = [(factor, by_beam), ...]. Every diagram must have the
    same beams and the same stations (use the fixed station grid for this).
    """
    if not parts:
        raise CombineError("nothing to combine")
    first = parts[0][1]
    for _, other in parts[1:]:
        if len(other) != len(first):
            raise CombineError("the diagrams have different numbers of beams")
    out = []
    for i, (chain_beam, pts) in enumerate(first):
        new_pts = []
        for j, (x, _) in enumerate(pts):
            total = dict.fromkeys(QUANTITIES, 0.0)
            for factor, other in parts:
                cb, op = other[i]
                if cb.beam != chain_beam.beam or len(op) != len(pts):
                    raise CombineError(f"beam {chain_beam.beam}: the stations differ between cases")
                xx, vals = op[j]
                if abs(xx - x) > 1e-6 * max(1.0, abs(x)):
                    raise CombineError(f"beam {chain_beam.beam}: station {j} is at x = {xx} in "
                                       f"one case and x = {x} in another")
                for q in QUANTITIES:
                    total[q] += factor * vals[q]
            new_pts.append((x, total))
        out.append((chain_beam, new_pts))
    return out


@dataclass
class StationEval:
    beam: int
    x: float
    forces: dict                 # kN, kNm
    result: object               # cross_section.CrossSectionResult


def evaluate_member(sec, steel, gM0, by_beam):
    """Check every station. Returns [StationEval]."""
    out = []
    for chain_beam, pts in by_beam:
        for x, f in pts:
            forces_si = {q: f[q] * 1e3 for q in QUANTITIES}      # kN, kNm -> N, N*m
            out.append(StationEval(chain_beam.beam, x, f,
                                   check_cross_section(sec, steel, gM0, forces_si)))
    return out


def status_of(util, not_checked):
    """FAILS wins even if part of the check was not possible; OK only if everything was checked."""
    if util is not None and util > 1.0:
        return "FAILS (incomplete)" if not_checked else "FAILS"
    if util is None:
        return "NOT CHECKED"
    if not_checked:
        return "INCOMPLETE"
    return "OK"


def status_rank(status, util):
    """
    Sort key for choosing what to show: a failure first (even if the check was also incomplete),
    then an incomplete check, then the highest utilisation. An incomplete check must never hide
    a failing one, nor look better than a complete one.
    """
    if status.startswith("FAILS"):
        group = 0
    elif status in ("INCOMPLETE", "NOT CHECKED"):
        group = 1
    else:
        group = 2
    return (group, -(util if util is not None else -1.0))


@dataclass
class Summary:
    member: str
    design_case: str
    n_stations: int
    max_class: int
    governing_check: str | None
    governing_util: float | None
    beam: int | None
    x: float | None
    max_by_check: dict = field(default_factory=dict)
    not_checked: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    @property
    def status(self):
        return status_of(self.governing_util, self.not_checked)


def summarise(member_id, design_case, evals):
    """Governing station and utilisation over the stations of one member and design case."""
    best = None
    max_by_check = {}
    not_checked, notes = [], []
    for e in evals:
        r = e.result
        for reason in r.not_checked:
            if reason not in not_checked:
                not_checked.append(reason)
        for note in r.notes:
            if note not in notes:
                notes.append(note)
        for name, value in r.util.items():
            if value > max_by_check.get(name, -1.0):
                max_by_check[name] = value
        g_name, g_val = r.governing
        if g_val is not None and (best is None or g_val > best[0]):
            best = (g_val, g_name, e)
    max_class = max((e.result.cls for e in evals), default=0)
    if best is None:
        return Summary(member_id, design_case, len(evals), max_class, None, None, None, None,
                       max_by_check, not_checked, notes)
    g_val, g_name, e = best
    return Summary(member_id, design_case, len(evals), max_class, g_name, g_val, e.beam, e.x,
                   max_by_check, not_checked, notes)


def buckling_check(sec, steel, gM1, bk, evals, member_length, lengths):
    """
    Member buckling (EC3 6.3.1 and 6.3.2) for one member and design case, from its station
    evaluations. lengths = {"Lcr_y", "Lcr_z", "Lw"} in m (None = switched off).
    """
    stations = [(e.x, e.forces["N"] * 1e3, e.forces["My"] * 1e3) for e in evals]
    mz_ed = max((abs(e.forces["Mz"]) * 1e3 for e in evals), default=0.0)
    cls = max((e.result.cls for e in evals), default=0)
    return check_member_buckling(sec, steel, gM1, bk, stations, member_length, cls,
                                 lengths["Lcr_y"], lengths["Lcr_z"], lengths["Lw"], mz_ed=mz_ed)


def pick_station(evals, x=None):
    """
    The station for a worked calculation: the one nearest x (m), or with x = None the governing
    station (if nothing could be checked, the one with the largest |My|). Returns (StationEval,
    description).
    """
    if x is not None:
        return min(evals, key=lambda ev: abs(ev.x - x)), f"station nearest x = {x:g} m"
    scored = [(e.result.governing[1], e) for e in evals if e.result.governing[1] is not None]
    e = (max(scored, key=lambda t: t[0])[1] if scored
         else max(evals, key=lambda ev: abs(ev.forces["My"])))
    return e, "governing station"
