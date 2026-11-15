"""
Step 4: Assemble member force diagrams from beam results.  Read-only.

For each member in members.csv and each chosen result case:
  - reads the beam forces (Strand7 principal axes, End Force convention)
  - converts every beam into the MEMBER's axes and sign convention
    (using the sign factors worked out in members.py)
  - joins the beams into one diagram along the member (x = distance from the
    member's start node)
  - checks the diagram is continuous at the joins between beams
  - prints the minimum and maximum of each force and writes member_forces.csv

Quantities, tension positive, in kN and kNm (if the model is in N and m):
  N   axial force
  Vy  shear along axis 1  (Strand7 SF1; with shear area SA1)
  Vz  shear along axis 2  (Strand7 SF2; with shear area SA2, the web)
  T   torque
  My  moment about axis 1 (Strand7 BM2: MAJOR-axis bending)
  Mz  moment about axis 2 (Strand7 BM1: minor-axis bending)

Settings are at the top of this file.
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))  # st7_setup changes the working dir

import csv  # noqa: E402
import ctypes  # noqa: E402
from contextlib import contextmanager  # noqa: E402

from st7_setup import St7API as S, check, model_session, MODEL_FILE  # noqa: E402
from members import (MEMBER_FILE, MemberError, read_rows, parse_row,  # noqa: E402
                     build_member, load_beam_ids, get_str)

# ---------------- SETTINGS (edit these) ----------------
RESULT_FILE = None          # e.g. r"D:\Scratch\test_model.LSA"; None = next to the model
RESULT_CASES = []           # e.g. [3]; [] = list the cases and stop
COMBINATIONS = "existing"   # "none", "existing" or "generate"
MIN_STATIONS = 9            # minimum stations per beam (Strand7 adds more at maxima)
PRINT_MEMBERS = []          # e.g. ["1"] to print the whole assembled diagram of these members
OUTPUT_CSV = "member_forces.csv"
# -------------------------------------------------------

QUANTITIES = ["N", "Vy", "Vz", "T", "My", "Mz"]


@contextmanager
def result_file(uid, path, combination_code):
    n_primary = ctypes.c_long()
    n_secondary = ctypes.c_long()
    check(S.St7OpenResultFile(uid, path.encode(), b"", combination_code,
                              ctypes.byref(n_primary), ctypes.byref(n_secondary)))
    try:
        yield n_primary.value, n_secondary.value
    finally:
        S.St7CloseResultFile(uid)


def force_scale(uid):
    units = (ctypes.c_long * S.kLastUnit)()
    check(S.St7GetUnits(uid, units))
    if units[S.ipLENGTHU] == S.luMETRE and units[S.ipFORCEU] == S.fuNEWTON:
        return 1e-3, "kN, kNm"
    print("WARNING: model is not in N and m; values are left in model units.")
    return 1.0, "model units"


def load_members(uid):
    total = ctypes.c_long()
    check(S.St7GetTotal(uid, S.tyBEAM, ctypes.byref(total)))
    n_beams = total.value
    rows = read_rows(MEMBER_FILE)
    need_ids = any("ID=" in r.get("beams", "").upper() for r in rows)
    beam_ids = load_beam_ids(uid, n_beams) if need_ids else {}
    members, seen = [], set()
    for row in rows:
        member_id = row.get("member", "?")
        if member_id in seen:
            print(f"Member {member_id}: ERROR - member number used more than once; skipped")
            continue
        seen.add(member_id)
        try:
            members.append(build_member(uid, parse_row(row, beam_ids, n_beams)))
        except (ValueError, MemberError) as e:
            print(f"Member {member_id}: ERROR - {e}; skipped")
    return members


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
                "T": c.f3 * torque * scale,            # torque: flips with axis 3
                "My": c.f2 * bm2 * scale,              # moment: flips with the shear axis only
                "Mz": c.f1 * bm1 * scale,
            }))
        pts.sort(key=lambda p: p[0])
        by_beam.append((c, pts))
    return by_beam


def continuity_warnings(m, by_beam):
    """Moments and torque must be continuous across a node unless a moment is applied there."""
    peak = max((abs(v[q]) for _, pts in by_beam for _, v in pts for q in ("My", "Mz", "T")),
               default=0.0)
    tol = max(1e-3, 1e-3 * peak)
    out = []
    for i in range(1, len(by_beam)):
        (c1, pts1), (c2, pts2) = by_beam[i - 1], by_beam[i]
        left, right = pts1[-1][1], pts2[0][1]
        for q in ("My", "Mz", "T"):
            if abs(left[q] - right[q]) > tol:
                out.append(f"{q} jumps at node {m.nodes[i]} between beam {c1.beam} "
                           f"({left[q]:.3f}) and beam {c2.beam} ({right[q]:.3f})")
    return out


def main():
    result_path = RESULT_FILE or (os.path.splitext(MODEL_FILE)[0] + ".LSA")
    if not os.path.isfile(result_path):
        raise FileNotFoundError(f"Result file not found: {result_path}")
    comb_code = {"none": S.kNoCombinations, "existing": S.kUseExistingCombinations,
                 "generate": S.kGenerateNewCombinations}[COMBINATIONS]

    with model_session() as uid:
        scale, unit_label = force_scale(uid)
        members = load_members(uid)
        if not members:
            print("No members could be built; fix members.csv first (run members.py).")
            return

        with result_file(uid, result_path, comb_code) as (n_prim, n_sec):
            print(f"Model        : {MODEL_FILE}")
            print(f"Result file  : {result_path}")
            print("Available result cases:")
            for c in range(1, n_prim + n_sec + 1):
                kind = "primary" if c <= n_prim else "combination"
                print(f"  {c:>3}  {get_str(S.St7GetResultCaseName, uid, c)}   [{kind}]")
            print()
            if not RESULT_CASES:
                print("RESULT_CASES is empty: set it at the top of this file "
                      "(for example RESULT_CASES = [3]) and run again.")
                return

            check(S.St7SetBeamResultPosMode(uid, S.bpLength))
            csv_rows = []
            for case in RESULT_CASES:
                case_name = get_str(S.St7GetResultCaseName, uid, case)
                print("=" * 76)
                print(f"Result case {case}: {case_name}      ({unit_label}; tension positive)")
                for m in members:
                    d = m.definition
                    by_beam = assemble(uid, m, case, scale)
                    points = [(c.beam, x, v) for c, pts in by_beam for x, v in pts]
                    for c, pts in by_beam:
                        for x, v in pts:
                            csv_rows.append([case, case_name, d.member_id, c.beam, x]
                                            + [v[q] for q in QUANTITIES])

                    print(f"\n  Member {d.member_id}   length {m.length:.3f} m   "
                          f"beams {' '.join(str(c.beam) for c in m.chain)}")
                    print(f"    {'':>3} {'min':>10} {'at x':>8} {'max':>10} {'at x':>8}")
                    for q in QUANTITIES:
                        lo = min(points, key=lambda p: p[2][q])
                        hi = max(points, key=lambda p: p[2][q])
                        print(f"    {q:>3} {lo[2][q]:>10.3f} {lo[1]:>8.3f} "
                              f"{hi[2][q]:>10.3f} {hi[1]:>8.3f}")
                    for w in continuity_warnings(m, by_beam):
                        print(f"    CHECK: {w}")

                    if str(d.member_id) in [str(x) for x in PRINT_MEMBERS]:
                        print(f"\n    Assembled diagram, member {d.member_id} "
                              f"(x from start node {m.nodes[0]}):")
                        print(f"    {'beam':>5} {'x':>8} " + " ".join(f"{q:>10}" for q in QUANTITIES))
                        for c, pts in by_beam:
                            for x, v in pts:
                                print(f"    {c.beam:>5} {x:>8.3f} "
                                      + " ".join(f"{v[q]:>10.3f}" for q in QUANTITIES))

            out_path = os.path.join(HERE, OUTPUT_CSV)
            with open(out_path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["case", "case_name", "member", "beam", "x"] + QUANTITIES)
                w.writerows(csv_rows)
            print(f"\nWrote {len(csv_rows)} rows to {out_path}")
            print("If no CHECK lines appear, moments and torque are continuous at every "
                  "join between beams.")


if __name__ == "__main__":
    main()
