"""
Step 5: EC3 cross-section checks (clause 6.2) on the members in members.csv,
for the design cases in design_cases.csv.  Read-only: the model is not changed.

For every member and design case it
  - reads the beam results from the .LSA result file(s) and converts them to member axes
    (a design case that adds cases together is evaluated on a fixed set of stations
     along each beam; a case taken directly from Strand7 uses Strand7's own stations,
     which include the maxima and minima)
  - classifies the section and checks axial force, bending, shear and their combinations
    at every station (see ec3_checks.py for exactly what is and is not covered)
  - prints the governing station for each member, and (optionally) a step-by-step
    worked calculation for one station so it can be checked by hand
  - writes every station to ec3_results.csv

Member buckling checks (clause 6.3) are NOT included yet.

Settings are at the top of this file; design settings are in settings.toml.
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))  # st7_setup changes the working dir

import csv  # noqa: E402
import ctypes  # noqa: E402
import math  # noqa: E402

from st7_setup import St7API as S, check, model_session, MODEL_FILE  # noqa: E402
from members import get_str, effective_lengths  # noqa: E402
from member_forces import load_members, assemble, force_scale, result_file  # noqa: E402
from section import read_i_section, SectionError  # noqa: E402
from section_math import short_name  # noqa: E402
from settings import load_settings, load_grades, strength, SettingsError  # noqa: E402
from design_cases import load_design_cases, file_warnings, DesignCaseError  # noqa: E402
from design_eval import (combine, evaluate_member, summarise, CombineError,  # noqa: E402
                         buckling_check, QUANTITIES, CHECKS)
from ec3_checks import ISection, check_station  # noqa: E402

# ---------------- SETTINGS (edit these) ----------------
RESULT_FILE = None          # default result file; None = the model's name with .LSA
COMBINATIONS = "existing"   # "none", "existing" or "generate" (Strand7 combination cases)
DESIGN_CASES = []           # names of design cases to run; [] = all in design_cases.csv
STATIONS_PER_BEAM = 21      # stations per beam when design cases are added together (odd is best)
SHOW_DETAIL = False         # True = also print a table per member and design case
WORKED = []                 # cross-section worked calculations, e.g. [("4", "TEST", None)]
                            #   (member, design case, x in m)  x = None: the governing station
WORKED_BUCKLING = []        # member buckling worked calculations, e.g. [("4", "ULS_S7")]
OUTPUT_CSV = "ec3_results.csv"
OUTPUT_BUCKLING_CSV = "ec3_buckling.csv"
# -------------------------------------------------------


def fmt_util(u):
    if u is None:
        return "-"
    if math.isinf(u):
        return "inf"
    return f"{u:.3f}"


# ----------------------------------------------------------------------
# Sections and materials
# ----------------------------------------------------------------------

def member_section(uid, m, settings, grades):
    """(ISection, grade, fy in Pa, warnings) for a member."""
    s = read_i_section(uid, m.chain[0].prop)
    sec = ISection.from_dims(short_name(s["name"]), s["h"], s["b"], s["tw"], s["tf"], s["r"], settings.eta)
    if settings.material_source == "model":
        sec.E, sec.G = s["E"], s["G"]
        if not (150e9 <= sec.E <= 250e9) or not (0.30 <= sec.G / sec.E <= 0.55):
            raise SectionError(f"the beam material has E = {sec.E:g}, G = {sec.G:g}: these should "
                               f"be about 210e9 and 81e9 Pa (model in N and m). Check the "
                               f"material, or set source = \"ec3\" in settings.toml")
    else:
        sec.E, sec.G = settings.E_ec3 * 1e6, settings.G_ec3 * 1e6
    grade = m.definition.grade or settings.default_grade
    fy, _fu = strength(grades, grade, sec.t_max * 1e3)
    warnings = []
    if not s["r_known"]:
        warnings.append("no root radius in the section data: r = 0 was used, so Wpl and Av "
                        "are slightly low")
    return sec, grade, fy * 1e6, warnings


# ----------------------------------------------------------------------
# Getting the forces
# ----------------------------------------------------------------------

def norm(path):
    return os.path.normcase(os.path.abspath(path))


def plan_extraction(cases, default_path):
    """{file key: (path, {(result case, 'native' or 'grid')})} - what has to be read from where."""
    needs = {}
    for dc in cases.values():
        mode = "native" if dc.is_direct else "grid"
        for row in dc.rows:
            path = row.result_file or default_path
            needs.setdefault(norm(path), (path, set()))[1].add((row.result_case, mode))
    return needs


def extract(uid, members, needs, scale, comb_code, grid):
    """Read every needed result case for every member. Returns {(file, case, mode, member): by_beam}."""
    data = {}
    for key, (path, wanted) in needs.items():
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Result file not found: {path}")
        if os.path.splitext(path)[1].lower() != ".lsa":
            raise DesignCaseError(f"{path}: only linear static (.LSA) result files are supported "
                                  f"so far")
        with result_file(uid, path, comb_code) as (n_prim, n_sec):
            total = n_prim + n_sec
            for case, mode in sorted(wanted):
                if not 1 <= case <= total:
                    raise DesignCaseError(f"result case {case} does not exist in {path} "
                                          f"(it has {total})")
                for m in members:
                    data[(key, case, mode, m.definition.member_id)] = assemble(
                        uid, m, case, scale, grid if mode == "grid" else None)
    return data


def design_case_forces(dc, member_id, data, default_path):
    """The member force diagram for one design case."""
    def key(row):
        return norm(row.result_file or default_path)

    if dc.is_direct:
        row = dc.rows[0]
        return data[(key(row), row.result_case, "native", member_id)]
    return combine([(row.factor, data[(key(row), row.result_case, "grid", member_id)])
                    for row in dc.rows])


# ----------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------

def print_members(members, info):
    print("MEMBERS")
    for m in members:
        mid = m.definition.member_id
        sec, grade, fy, _ = info[mid]
        print(f"  Member {mid}: {sec.name}   {grade}   fy = {fy / 1e6:g} N/mm2 "
              f"(t = {sec.t_max * 1e3:.1f} mm)   length {m.length:.3f} m")
    print()


def print_summary(members, info, best):
    print("MEMBER SUMMARY - cross-section checks (EC3 6.2), governing design case and station")
    print(f"{'Member':<8}{'Section':<28}{'Grade':<7}{'Class':>5}  {'Design case':<14}"
          f"{'x [m]':>8}  {'Check':<6}{'Util':>7}  Status")
    for m in members:
        mid = m.definition.member_id
        sec, grade, fy, _ = info[mid]
        s = best.get(mid)
        short = (sec.name[:24] + "..") if len(sec.name) > 26 else sec.name
        if s is None:
            print(f"{mid:<8}{short:<28}{grade:<7}{'-':>5}  {'-':<14}{'-':>8}  {'-':<6}{'-':>7}  "
                  f"NOT CHECKED (no design case)")
            continue
        x = f"{s.x:.3f}" if s.x is not None else "-"
        print(f"{mid:<8}{short:<28}{grade:<7}{s.max_class:>5}  {s.design_case:<14}"
              f"{x:>8}  {s.governing_check or '-':<6}{fmt_util(s.governing_util):>7}  {s.status}")
    print()


def print_detail(members, summaries):
    print("DETAIL - maximum utilisation by check, per member and design case")
    print(f"{'Member':<8}{'Design case':<14}{'Class':>5}" + "".join(f"{c:>8}" for c in CHECKS)
          + "  Status")
    for m in members:
        mid = m.definition.member_id
        for s in summaries[mid]:
            print(f"{mid:<8}{s.design_case:<14}{s.max_class:>5}"
                  + "".join(f"{fmt_util(s.max_by_check.get(c)):>8}" for c in CHECKS)
                  + f"  {s.status}")
    print()


def print_notes(members, info, best, summaries):
    printed = False
    for m in members:
        mid = m.definition.member_id
        lines = [f"WARNING: {w}" for w in info[mid][3]] + [f"WARNING: {w}" for w in m.warnings]
        seen = set()
        for s in summaries.get(mid, []):
            for reason in s.not_checked:
                if reason not in seen:
                    seen.add(reason)
                    lines.append(f"NOT FULLY CHECKED ({s.design_case}): {reason}")
        for s in summaries.get(mid, []):
            for note in s.notes:
                if note not in seen:
                    seen.add(note)
                    lines.append(f"note ({s.design_case}): {note}")
        if lines:
            printed = True
            print(f"Member {mid}:")
            for ln in lines:
                print("  " + ln)
    if printed:
        print()


def buckling_status(res):
    name, value = res.governing
    if value is not None and value > 1.0:
        return "FAILS (incomplete)" if res.not_checked else "FAILS"
    if value is None:
        return "NOT CHECKED"
    return "INCOMPLETE" if res.not_checked else "OK"


def print_buckling(members, buckling, cases):
    print("MEMBER BUCKLING (EC3 6.3.1 flexural, 6.3.2 LTB) - governing design case")
    print("  N_Ed = largest compression and My,Ed = largest |My| anywhere in the member")
    print(f"{'Member':<8}{'Design case':<14}{'N_Ed kN':>10}{'My,Ed kNm':>11}"
          f"{'Nb,y':>8}{'Nb,z':>8}{'LTB':>8}  Status")
    notes = []
    for m in members:
        mid = m.definition.member_id
        rows = [(name, buckling[(mid, name)]) for name in cases if (mid, name) in buckling]
        if not rows:
            continue

        def rank(item):
            res = item[1]
            st = buckling_status(res)
            group = 0 if st.startswith("FAILS") else 1 if st in ("INCOMPLETE", "NOT CHECKED") else 2
            value = res.governing[1]
            return (group, -(value if value is not None else -1.0))

        name, res = min(rows, key=rank)
        v = res.values
        cells = "".join(f"{fmt_util(res.util.get(k)):>8}" for k in ("Nb,y", "Nb,z", "LTB"))
        print(f"{mid:<8}{name:<14}{v.get('N_Ed', 0) / 1e3:>10.1f}{v.get('My_Ed', 0) / 1e3:>11.2f}"
              f"{cells}  {buckling_status(res)}")
        seen = set()
        for _, r in rows:
            for reason in r.not_checked:
                if reason not in seen:
                    seen.add(reason)
                    notes.append(f"Member {mid}: NOT FULLY CHECKED: {reason}")
            for note in r.notes:
                if note not in seen:
                    seen.add(note)
                    notes.append(f"Member {mid}: note: {note}")
    print()
    for n in notes:
        print(n)
    if notes:
        print()


def write_buckling_csv(path, members, buckling, cases):
    cols = ["member", "design_case", "N_Ed_kN", "My_Ed_kNm", "Lcr_y_m", "Lcr_z_m", "Lw_m",
            "chi_y", "chi_z", "C1", "Mcr_kNm", "chi_LT", "U_Nb_y", "U_Nb_z", "U_LTB", "status"]

    def num(v, scale=1.0, f="{:.4f}"):
        return "" if v is None else f.format(v * scale)

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for m in members:
            mid = m.definition.member_id
            for name in cases:
                r = buckling.get((mid, name))
                if r is None:
                    continue
                v = r.values
                w.writerow([mid, name, num(v.get("N_Ed"), 1e-3), num(v.get("My_Ed"), 1e-3),
                            num(v.get("Lcr_y")), num(v.get("Lcr_z")), num(v.get("Lw")),
                            num(v.get("chi_y")), num(v.get("chi_z")), num(v.get("C1")),
                            num(v.get("Mcr"), 1e-3), num(v.get("chi_LT")),
                            fmt_util(r.util.get("Nb,y")), fmt_util(r.util.get("Nb,z")),
                            fmt_util(r.util.get("LTB")), buckling_status(r)])


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["member", "design_case", "beam", "x_m"]
                   + [f"{q}_{'kN' if q[0] in 'NV' else 'kNm'}" for q in QUANTITIES]
                   + ["class"] + [f"U_{c}" for c in CHECKS]
                   + ["governing", "util", "not_checked"])
        w.writerows(rows)


def worked_calculations(requests, members, info, all_evals):
    by_id = {m.definition.member_id: m for m in members}
    for member_id, case_name, x in requests:
        member_id = str(member_id)
        evals = all_evals.get((member_id, case_name))
        print("=" * 76)
        if member_id not in by_id or evals is None:
            print(f"WORKED CALCULATION: member {member_id}, design case {case_name}: "
                  f"nothing to show (not in this run)")
            continue
        if x is None:
            scored = [(e.result.governing[1], e) for e in evals if e.result.governing[1] is not None]
            e = max(scored, key=lambda t: t[0])[1] if scored else \
                max(evals, key=lambda ev: abs(ev.forces["My"]))
            why = "governing station"
        else:
            e = min(evals, key=lambda ev: abs(ev.x - x))
            why = f"station nearest x = {x:g} m"
        sec, grade, fy, _ = info[member_id]
        forces_si = {q: e.forces[q] * 1e3 for q in QUANTITIES}
        r = check_station(sec, fy, info["gamma_M0"], forces_si, trace=True)
        print(f"WORKED CALCULATION: member {member_id}, design case {case_name}, beam {e.beam}, "
              f"x = {e.x:.3f} m ({why}), grade {grade}")
        print("=" * 76)
        print("\n".join(r.trace))
        print()


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    settings = load_settings()
    grades = load_grades()
    all_cases = load_design_cases()
    if DESIGN_CASES:
        missing = [n for n in DESIGN_CASES if n not in all_cases]
        if missing:
            raise DesignCaseError(f"design case(s) not in design_cases.csv: {', '.join(missing)}")
        cases = {n: all_cases[n] for n in DESIGN_CASES}
    else:
        cases = all_cases
    if STATIONS_PER_BEAM < 3:
        raise ValueError("STATIONS_PER_BEAM must be at least 3")
    grid = [i / (STATIONS_PER_BEAM - 1) for i in range(STATIONS_PER_BEAM)]

    default_path = RESULT_FILE or (os.path.splitext(MODEL_FILE)[0] + ".LSA")
    comb_code = {"none": S.kNoCombinations, "existing": S.kUseExistingCombinations,
                 "generate": S.kGenerateNewCombinations}[COMBINATIONS]

    for w in file_warnings(cases):
        print(f"WARNING: {w}")

    with model_session() as uid:
        scale, unit_label = force_scale(uid)
        if unit_label != "kN, kNm":
            raise RuntimeError("the checks need a model in N and m (forces are converted to "
                               "kN and kNm); change the model units or extend force_scale()")
        members = load_members(uid)
        if not members:
            print("No members could be built; fix members.csv first (run members.py).")
            return

        n_beams = ctypes.c_long()
        check(S.St7GetTotal(uid, S.tyBEAM, ctypes.byref(n_beams)))

        info = {"gamma_M0": settings.gamma_M0}
        usable = []
        for m in members:
            mid = m.definition.member_id
            try:
                info[mid] = member_section(uid, m, settings, grades)
            except (SectionError, SettingsError, RuntimeError) as e:
                print(f"Member {mid}: NOT CHECKED - {e}")
                continue
            usable.append(m)

        needs = plan_extraction(cases, default_path)
        data = extract(uid, usable, needs, scale, comb_code, grid)

    # ---- the Strand7 model is closed from here on; everything below is calculation ----
    print(f"Model          : {MODEL_FILE}")
    print(f"Default results: {default_path}")
    print(f"Design cases   : {', '.join(cases)}")
    print(f"gamma_M0 = {settings.gamma_M0:g}   default grade = {settings.default_grade}   "
          f"eta = {settings.eta:g}")
    print("Cross-section checks (EC3 6.2) and member buckling (6.3.1 flexural, 6.3.2 LTB, "
          "general case).")
    print("NOT INCLUDED YET: the bending and axial force interaction for buckling (6.3.3).")
    print()

    summaries = {m.definition.member_id: [] for m in usable}
    buckling = {}
    all_evals = {}
    csv_rows = []
    for m in usable:
        mid = m.definition.member_id
        sec, grade, fy, _ = info[mid]
        for name, dc in cases.items():
            try:
                by_beam = design_case_forces(dc, mid, data, default_path)
            except CombineError as e:
                print(f"Member {mid}, design case {name}: ERROR - {e}")
                continue
            evals = evaluate_member(sec, fy, settings.gamma_M0, by_beam)
            all_evals[(mid, name)] = evals
            summaries[mid].append(summarise(mid, name, evals))
            buckling[(mid, name)] = buckling_check(
                sec, fy, sec.E, sec.G, settings.gamma_M1, settings.buckling, evals, m.length,
                effective_lengths(m), trace=True)
            for e in evals:
                r = e.result
                g_name, g_val = r.governing
                csv_rows.append([mid, name, e.beam, f"{e.x:.4f}"]
                                + [f"{e.forces[q]:.4f}" for q in QUANTITIES]
                                + [r.cls] + [fmt_util(r.util.get(c)) for c in CHECKS]
                                + [g_name or "", fmt_util(g_val), "; ".join(r.not_checked)])

    # Which design case to show for each member: a failure first (even if the check was also
    # incomplete), then an incomplete check, then the highest utilisation. An incomplete check
    # must never hide a failing one, nor look better than a complete one.
    def rank(s):
        if s.status.startswith("FAILS"):
            group = 0
        elif s.status in ("INCOMPLETE", "NOT CHECKED"):
            group = 1
        else:
            group = 2
        return (group, -(s.governing_util if s.governing_util is not None else -1.0))

    best = {mid: min(lst, key=rank) for mid, lst in summaries.items() if lst}

    print_members(usable, info)
    print_summary(usable, info, best)
    if SHOW_DETAIL:
        print_detail(usable, summaries)
    print_notes(usable, info, best, summaries)
    print_buckling(usable, buckling, cases)

    in_members = {b for m in usable for b in (c.beam for c in m.chain)}
    unassigned = [b for b in range(1, n_beams.value + 1) if b not in in_members]
    if unassigned:
        shown = ", ".join(str(b) for b in unassigned[:30]) + (" ..." if len(unassigned) > 30 else "")
        print(f"Beams not in any member ({len(unassigned)}) have NOT been checked: {shown}\n")

    out_path = os.path.join(HERE, OUTPUT_CSV)
    write_csv(out_path, csv_rows)
    print(f"Wrote {len(csv_rows)} station results to {out_path}\n")

    write_buckling_csv(os.path.join(HERE, OUTPUT_BUCKLING_CSV), usable, buckling, cases)
    print(f"Wrote member buckling results to {os.path.join(HERE, OUTPUT_BUCKLING_CSV)}\n")

    if WORKED:
        worked_calculations(WORKED, usable, info, all_evals)
    for member_id, case_name in WORKED_BUCKLING:
        res = buckling.get((str(member_id), case_name))
        print("=" * 76)
        if res is None:
            print(f"WORKED BUCKLING CALCULATION: member {member_id}, design case {case_name}: "
                  f"nothing to show (not in this run)")
            continue
        print(f"WORKED BUCKLING CALCULATION: member {member_id}, design case {case_name}")
        print("=" * 76)
        print("\n".join(res.trace))
        print()


if __name__ == "__main__":
    main()
