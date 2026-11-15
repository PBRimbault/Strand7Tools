"""
List the result cases in the result file, and print the assembled member force diagrams
(member axes, tension positive, kN and kNm) for chosen result cases. Read-only.

    python tools\\list_forces.py [project file]       (default: Test\\test_model2.toml)

With RESULT_CASES empty it only lists the result cases. Writes member_forces.csv to the
project's output folder. A CHECK line means a moment or torque jumps where two beams meet.
"""

import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from strand7tools.project import load_project  # noqa: E402
from strand7tools.quantities import QUANTITIES  # noqa: E402
from strand7tools.strand7.api import S, get_str, model_session  # noqa: E402
from strand7tools.strand7.forces import (assemble, continuity_warnings, force_scale,  # noqa: E402
                                         result_file)
from strand7tools.strand7.members import load_members  # noqa: E402

DEFAULT_PROJECT = os.path.join(ROOT, "Test", "test_model2.toml")

# ---------------- SETTINGS (edit these) ----------------
RESULT_CASES = []           # e.g. [3]; [] = list the cases and stop
PRINT_MEMBERS = []          # e.g. ["1"] to print the whole assembled diagram of these members
OUTPUT_CSV = "member_forces.csv"
# -------------------------------------------------------


def main(project_file=DEFAULT_PROJECT):
    project = load_project(project_file)
    if not os.path.isfile(project.result_file):
        raise FileNotFoundError(f"Result file not found: {project.result_file}")

    with model_session(project.model) as uid:
        scale = force_scale(uid)
        unit_label = "kN, kNm"
        if scale is None:
            print("WARNING: model is not in m with N, kN or MN; values are left in model units.")
            scale, unit_label = 1.0, "model units"
        members, errors = load_members(uid, project.members)
        for mid, msg in errors:
            print(f"Member {mid}: ERROR - {msg}; skipped")
        if not members:
            print("No members could be built; fix the member table first (tools/list_members.py).")
            return

        with result_file(uid, project.result_file, project.combinations) as (n_prim, n_sec):
            print(f"Model        : {project.model}")
            print(f"Result file  : {project.result_file}")
            print("Available result cases:")
            for c in range(1, n_prim + n_sec + 1):
                kind = "primary" if c <= n_prim else "combination"
                print(f"  {c:>3}  {get_str(S.St7GetResultCaseName, uid, c)}   [{kind}]")
            print()
            if not RESULT_CASES:
                print(f"RESULT_CASES is empty: set it at the top of {__file__} "
                      "(for example RESULT_CASES = [3]) and run again.")
                return

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

    out_path = project.output(OUTPUT_CSV)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "case_name", "member", "beam", "x"] + QUANTITIES)
        w.writerows(csv_rows)
    print(f"\nWrote {len(csv_rows)} rows to {out_path}")
    print("If no CHECK lines appear, moments and torque are continuous at every "
          "join between beams.")


if __name__ == "__main__":
    main(*sys.argv[1:2])
