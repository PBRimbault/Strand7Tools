"""
Member force diagrams from beam results. Read-only.

For a member and a result case:
  - reads the beam forces (Strand7 principal axes, End Force convention)
  - converts every beam into the MEMBER's axes and sign convention
    (using the sign factors worked out in members.py)
  - joins the beams into one diagram along the member (x = distance from the
    member's start node)

Quantities (see quantities.py), tension positive, in kN and kNm (model in m, with N, kN or MN).
"""

import ctypes
from contextlib import contextmanager

from .api import S, check, model_units

MIN_STATIONS = 9            # minimum stations per beam with Strand7's own stations
                            # (Strand7 adds more at maxima and minima)
COMBINATION_CODES = {"none": "kNoCombinations", "existing": "kUseExistingCombinations",
                     "generate": "kGenerateNewCombinations"}


@contextmanager
def result_file(uid, path, combinations="existing"):
    """Open a result file; yields (number of primary cases, number of secondary cases)."""
    code = getattr(S, COMBINATION_CODES[combinations])
    n_primary = ctypes.c_long()
    n_secondary = ctypes.c_long()
    check(S.St7OpenResultFile(uid, path.encode(), b"", code,
                              ctypes.byref(n_primary), ctypes.byref(n_secondary)))
    try:
        yield n_primary.value, n_secondary.value
    finally:
        S.St7CloseResultFile(uid)


FORCE_TO_KN = {"fuNEWTON": 1e-3, "fuKILONEWTON": 1.0, "fuMEGANEWTON": 1e3}


def force_scale(uid):
    """
    Factor from model forces to kN, or None if the model is not in metres with forces in N, kN
    or MN (the checks need that).
    """
    units = model_units(uid)
    if units[S.ipLENGTHU] == S.luMETRE:
        for name, factor in FORCE_TO_KN.items():
            if units[S.ipFORCEU] == getattr(S, name):
                return factor
    return None


def beam_forces(uid, beam, case, grid=None, length=None):
    """
    [(distance from end 1, [SF1, BM1, SF2, BM2, axial, torque]), ...] in model units.

    grid=None : Strand7 chooses the stations (it adds extra ones at maxima and minima).
    grid=[fractions 0..1] : results at exactly these stations (needed to add cases together,
                            because every case must use the same stations). length is the
                            beam length, used to turn the fractions into distances.
    """
    idx = [S.ipBeamSF1, S.ipBeamBM1, S.ipBeamSF2, S.ipBeamBM2, S.ipBeamAxialF, S.ipBeamTorque]
    n_cols = ctypes.c_long()
    pos = (ctypes.c_double * S.kMaxBeamResult)()
    res = (ctypes.c_double * S.kMaxBeamResult)()

    if grid is None:
        n_stations = ctypes.c_long()
        check(S.St7SetBeamResultPosMode(uid, S.bpLength))   # positions as lengths
        check(S.St7GetBeamResultArray(uid, S.rtBeamForce, S.stBeamPrincipal, beam,
                                      MIN_STATIONS, case, ctypes.byref(n_stations),
                                      ctypes.byref(n_cols), pos, res))
        n, positions = n_stations.value, [pos[i] for i in range(n_stations.value)]
    else:
        if length is None:
            raise ValueError("beam_forces: a grid needs the beam length")
        if len(grid) * 8 > S.kMaxBeamResult:
            raise ValueError(f"too many stations per beam ({len(grid)}); "
                             f"the API allows at most {S.kMaxBeamResult // 8}")
        for i, f in enumerate(grid):
            pos[i] = f
        check(S.St7SetBeamResultPosMode(uid, S.bpParam))    # positions as fractions
        check(S.St7GetBeamResultArrayPos(uid, S.rtBeamForce, S.stBeamPrincipal, beam, case,
                                         len(grid), pos, ctypes.byref(n_cols), res))
        n, positions = len(grid), [f * length for f in grid]
    nc = n_cols.value
    return [(positions[i], [res[i * nc + k] for k in idx]) for i in range(n)]


def assemble(uid, m, case, scale, grid=None):
    """
    Per beam, the forces converted to member axes: [(BeamInChain, [(x, {quantity: value})])].
    grid: None = Strand7's own stations; or a list of fractions (see beam_forces).
    """
    by_beam = []
    for c in m.chain:
        pts = []
        for pos, (sf1, bm1, sf2, bm2, axial, torque) in beam_forces(
                uid, c.beam, case, grid, c.length):
            x = c.x0 + pos if c.forward else c.x1 - pos
            pts.append((x, {
                "N": axial * scale,                    # tension positive, direction-independent
                "Vy": c.f3 * c.f1 * sf1 * scale,       # shear: flips with axis 3 and the shear axis
                "Vz": c.f3 * c.f2 * sf2 * scale,
                "T": torque * scale,                   # torque: direction-independent, like N
                                                       # (reversing axis 3 reverses both the twist
                                                       # axis and d/dx; checked in check_portal.py)
                "My": c.f2 * bm2 * scale,              # moment: flips with the shear axis only
                "Mz": c.f1 * bm1 * scale,
            }))
        pts.sort(key=lambda p: p[0])
        by_beam.append((c, pts))
    return by_beam


def continuity_warnings(m, by_beam):
    """Moments and torque must be continuous across a node unless a moment is applied there."""
    # Tolerance per quantity, from that quantity's own peak, so a sign error in a small torque is
    # not hidden by a large bending moment; with a small floor from the largest of the three.
    peaks = {q: max((abs(v[q]) for _, pts in by_beam for _, v in pts), default=0.0)
             for q in ("My", "Mz", "T")}
    floor = 1e-6 * max(max(peaks.values()), 1e-3)
    out = []
    for i in range(1, len(by_beam)):
        (c1, pts1), (c2, pts2) = by_beam[i - 1], by_beam[i]
        left, right = pts1[-1][1], pts2[0][1]
        for q in ("My", "Mz", "T"):
            if abs(left[q] - right[q]) > max(1e-3 * peaks[q], floor):
                out.append(f"{q} jumps at node {m.nodes[i]} between beam {c1.beam} "
                           f"({left[q]:.3f}) and beam {c2.beam} ({right[q]:.3f})")
    return out
