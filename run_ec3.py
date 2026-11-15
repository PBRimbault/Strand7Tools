"""
EC3 checks (EN 1993-1-1 clause 6.2 cross-sections, 6.3.1 and 6.3.2 member buckling) on the
members of a Strand7 model, for the design cases of a project. Read-only: the model is not
changed.

    python run_ec3.py                       uses PROJECT below
    python run_ec3.py Test\\test_model3.toml

The project file names the model, the member table and the design cases (see
strand7tools/project.py). Design settings are in strand7tools/data/settings.toml.

Prints the governing results for each member, and (optionally) step-by-step worked
calculations so they can be checked by hand; writes every station and every member buckling
check to CSV files in the project's output folder.
"""

import os
import sys

from strand7tools.app import run_checks
from strand7tools.project import load_project
from strand7tools.report import (print_report, print_worked_buckling, print_worked_cross_section,
                                 write_buckling_csv, write_station_csv)

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------- SETTINGS (edit these) ----------------
PROJECT = os.path.join(HERE, "Test", "test_model2.toml")
SHOW_DETAIL = False         # True = also print a table per member and design case
WORKED = []                 # cross-section worked calculations, e.g. [("4", "TEST", None)]
                            #   (member, design case, x in m)  x = None: the governing station
WORKED_BUCKLING = []        # member buckling worked calculations, e.g. [("4", "ULS_S7")]
STATION_CSV = "ec3_results.csv"
BUCKLING_CSV = "ec3_buckling.csv"
# -------------------------------------------------------


def main(project_file=PROJECT):
    project = load_project(project_file)
    run = run_checks(project)
    print_report(run, SHOW_DETAIL)

    path = project.output(STATION_CSV)
    n = write_station_csv(path, run)
    print(f"Wrote {n} station results to {path}\n")
    path = project.output(BUCKLING_CSV)
    write_buckling_csv(path, run)
    print(f"Wrote member buckling results to {path}\n")

    for member_id, case_name, x in WORKED:
        print_worked_cross_section(run, member_id, case_name, x)
    for member_id, case_name in WORKED_BUCKLING:
        print_worked_buckling(run, member_id, case_name)


if __name__ == "__main__":
    main(*sys.argv[1:2])
