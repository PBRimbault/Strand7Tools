"""
Print a project's design inputs as they are read: the settings, a table of fy / fu by
thickness, and the design cases. No Strand7 needed.

    python tools\\show_inputs.py [project file]       (default: Test\\test_model2.toml)
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from strand7tools.design_cases import file_warnings, load_design_cases  # noqa: E402
from strand7tools.fmt import eng  # noqa: E402
from strand7tools.project import load_project  # noqa: E402
from strand7tools.settings import SettingsError, load_grades, load_settings, strength  # noqa: E402

DEFAULT_PROJECT = os.path.join(ROOT, "Test", "test_model2.toml")


def main(project_file=DEFAULT_PROJECT):
    project = load_project(project_file)
    s = load_settings(project.settings)
    grades = load_grades(project.steel_grades)

    print(f"Project       : {project.path}")
    print(f"  model           : {project.model}")
    print(f"  result file     : {project.result_file}")
    print(f"  member table    : {project.members}")
    print(f"  output folder   : {project.output_folder}")
    print()
    print(f"Settings file : {project.settings}")
    print(f"  material source : {s.material_source}")
    print(f"  E (if 'ec3')    : {eng(s.E_ec3)} N/mm2")
    print(f"  G (if 'ec3')    : {eng(s.G_ec3)} N/mm2")
    print(f"  default grade   : {s.default_grade}")
    print(f"  thickness rule  : {s.thickness_rule}")
    print(f"  gamma_M0, M1, M2: {s.gamma_M0:g}, {s.gamma_M1:g}, {s.gamma_M2:g}")
    print(f"  eta (Av,z)      : {s.eta:g}")
    b = s.buckling
    print(f"  buckling        : lambda_0 = {b.lambda_0:g}, lambda_LT0 = {b.lambda_LT0:g}, "
          f"beta_LT = {b.beta_LT:g}, alpha = {b.alpha}, alpha_LT = {b.alpha_LT}")
    if s.default_grade not in grades:
        print(f"  WARNING: default grade {s.default_grade!r} is not in the grade table")
    print()

    thicknesses = [10, 16, 20, 40, 50, 63, 80, 100]
    print(f"Grade table   : {project.steel_grades}")
    print("  fy / fu in N/mm2 by thickness (mm)")
    print(f"  {'grade':<7}" + "".join(f"{t:>10}" for t in thicknesses))
    for grade in grades:
        cells = []
        for t in thicknesses:
            try:
                fy, fu = strength(grades, grade, t)
                cells.append(f"{fy:g}/{fu:g}")
            except SettingsError:
                cells.append("-")
        print(f"  {grade:<7}" + "".join(f"{c:>10}" for c in cells))
    print()

    cases = load_design_cases(project.design_cases)
    print(f"Design case file: {project.design_cases}")
    print(f"{len(cases)} design cases\n")
    for case in cases.values():
        kind = ("direct from Strand7 (single row, factor 1)" if case.is_direct
                else "factored sum of linear results")
        print(f"Design case {case.name}   [{kind}]")
        for r in case.rows:
            where = r.result_file if r.result_file else "default result file"
            note = f"   {r.comment}" if r.comment else ""
            print(f"  {r.factor:>8.3f} x result case {r.result_case:<4} ({where}){note}")
        print()
    for w in file_warnings(cases):
        print(f"WARNING: {w}")


if __name__ == "__main__":
    main(*sys.argv[1:2])
