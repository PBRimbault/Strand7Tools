"""
Console report and CSV files for a RunResult (app.py). Presentation only: no checking here.
"""

import csv

from .ec3.design import CHECKS, pick_station
from .fmt import fmt_util
from .quantities import QUANTITIES
from .worked import buckling_lines, cross_section_lines

RULE = "=" * 76


def print_header(run):
    p, s = run.project, run.settings
    print(f"Project        : {p.path}")
    print(f"Model          : {p.model}")
    print(f"Default results: {p.result_file}")
    print(f"Design cases   : {', '.join(run.cases)}")
    print(f"gamma_M0 = {s.gamma_M0:g}   default grade = {s.default_grade}   eta = {s.eta:g}")
    print("Cross-section checks (EC3 6.2) and member buckling (6.3.1 flexural, 6.3.2 LTB, "
          "general case).")
    print("NOT INCLUDED YET: the bending and axial force interaction for buckling (6.3.3).")
    print()


def print_problems(run):
    for w in run.file_warnings:
        print(f"WARNING: {w}")
    for mid, msg in run.errors:
        print(f"Member {mid}: {msg}")
    for mid, name, msg in run.case_errors:
        print(f"Member {mid}, design case {name}: ERROR - {msg}")
    if not run.members:
        print("No members could be checked; fix the member table first "
              "(tools/list_members.py shows what was built).")
    if run.file_warnings or run.errors or run.case_errors or not run.members:
        print()


def print_members(run):
    print("MEMBERS")
    for m in run.members:
        mi = run.info[m.member_id]
        print(f"  Member {m.member_id}: {mi.section.name}   {mi.steel.grade}   "
              f"fy = {mi.steel.fy / 1e6:g} N/mm2 (t = {mi.section.t_max * 1e3:.1f} mm)   "
              f"length {m.length:.3f} m")
    print()
    for w in run.material_warnings:
        print(f"WARNING: {w}")
    if run.material_warnings:
        print()


def print_summary(run):
    print("MEMBER SUMMARY - cross-section checks (EC3 6.2), governing design case and station")
    print(f"{'Member':<8}{'Section':<28}{'Grade':<7}{'Class':>5}  {'Design case':<14}"
          f"{'x [m]':>8}  {'Check':<6}{'Util':>7}  Status")
    for m in run.members:
        mid = m.member_id
        mi = run.info[mid]
        s = run.governing_summary(mid)
        name = mi.section.name
        short = (name[:24] + "..") if len(name) > 26 else name
        grade = mi.steel.grade
        if s is None:
            print(f"{mid:<8}{short:<28}{grade:<7}{'-':>5}  {'-':<14}{'-':>8}  {'-':<6}{'-':>7}  "
                  f"NOT CHECKED (no design case)")
            continue
        x = f"{s.x:.3f}" if s.x is not None else "-"
        print(f"{mid:<8}{short:<28}{grade:<7}{s.max_class:>5}  {s.design_case:<14}"
              f"{x:>8}  {s.governing_check or '-':<6}{fmt_util(s.governing_util):>7}  {s.status}")
    print()


def print_detail(run):
    print("DETAIL - maximum utilisation by check, per member and design case")
    print(f"{'Member':<8}{'Design case':<14}{'Class':>5}" + "".join(f"{c:>8}" for c in CHECKS)
          + "  Status")
    for m in run.members:
        mid = m.member_id
        for s in run.summaries[mid]:
            print(f"{mid:<8}{s.design_case:<14}{s.max_class:>5}"
                  + "".join(f"{fmt_util(s.max_by_check.get(c)):>8}" for c in CHECKS)
                  + f"  {s.status}")
    print()


def print_notes(run):
    printed = False
    for m in run.members:
        mid = m.member_id
        lines = ([f"WARNING: {w}" for w in run.info[mid].warnings]
                 + [f"WARNING: {w}" for w in m.warnings])
        seen = set()
        summaries = run.summaries.get(mid, [])
        for s in summaries:
            for reason in s.not_checked:
                if reason not in seen:
                    seen.add(reason)
                    lines.append(f"NOT FULLY CHECKED ({s.design_case}): {reason}")
        for s in summaries:
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


def print_buckling(run):
    print("MEMBER BUCKLING (EC3 6.3.1 flexural, 6.3.2 LTB) - governing design case")
    print("  N_Ed = largest compression and My,Ed = largest |My| anywhere in the member")
    print(f"{'Member':<8}{'Design case':<14}{'N_Ed kN':>10}{'My,Ed kNm':>11}"
          f"{'Nb,y':>8}{'Nb,z':>8}{'LTB':>8}  Status")
    notes = []
    for m in run.members:
        mid = m.member_id
        gov = run.governing_buckling(mid)
        if gov is None:
            continue
        name, res = gov
        v = res.values
        cells = "".join(f"{fmt_util(res.util.get(k)):>8}" for k in ("Nb,y", "Nb,z", "LTB"))
        print(f"{mid:<8}{name:<14}{v.get('N_Ed', 0) / 1e3:>10.1f}{v.get('My_Ed', 0) / 1e3:>11.2f}"
              f"{cells}  {res.status}")
        seen = set()
        for case in run.cases:
            r = run.buckling.get((mid, case))
            if r is None:
                continue
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


def print_unassigned(run):
    b = run.unassigned_beams
    if b:
        shown = ", ".join(str(x) for x in b[:30]) + (" ..." if len(b) > 30 else "")
        print(f"Beams not in any member ({len(b)}) have NOT been checked: {shown}\n")


def print_report(run, show_detail=False):
    print_problems(run)
    print_header(run)
    print_members(run)
    print_summary(run)
    if show_detail:
        print_detail(run)
    print_notes(run)
    print_buckling(run)
    print_unassigned(run)


# ----------------------------------------------------------------------
# Worked calculations
# ----------------------------------------------------------------------

def print_worked_cross_section(run, member_id, case_name, x=None):
    member_id = str(member_id)
    evals = run.evals.get((member_id, case_name))
    print(RULE)
    if evals is None:
        print(f"WORKED CALCULATION: member {member_id}, design case {case_name}: "
              f"nothing to show (not in this run)")
        return
    e, why = pick_station(evals, x)
    print(f"WORKED CALCULATION: member {member_id}, design case {case_name}, beam {e.beam}, "
          f"x = {e.x:.3f} m ({why}), grade {run.info[member_id].steel.grade}")
    print(RULE)
    print("\n".join(cross_section_lines(e.result)))
    print()


def print_worked_buckling(run, member_id, case_name):
    res = run.buckling.get((str(member_id), case_name))
    print(RULE)
    if res is None:
        print(f"WORKED BUCKLING CALCULATION: member {member_id}, design case {case_name}: "
              f"nothing to show (not in this run)")
        return
    print(f"WORKED BUCKLING CALCULATION: member {member_id}, design case {case_name}")
    print(RULE)
    print("\n".join(buckling_lines(res)))
    print()


# ----------------------------------------------------------------------
# CSV files
# ----------------------------------------------------------------------

def write_station_csv(path, run):
    """Every station of every member and design case. Returns the number of rows."""
    rows = []
    for m in run.members:
        mid = m.member_id
        for name in run.cases:
            for e in run.evals.get((mid, name), []):
                r = e.result
                g_name, g_val = r.governing
                rows.append([mid, name, e.beam, f"{e.x:.4f}"]
                            + [f"{e.forces[q]:.4f}" for q in QUANTITIES]
                            + [r.cls] + [fmt_util(r.util.get(c)) for c in CHECKS]
                            + [g_name or "", fmt_util(g_val), "; ".join(r.not_checked)])
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["member", "design_case", "beam", "x_m"]
                   + [f"{q}_{'kN' if q[0] in 'NV' else 'kNm'}" for q in QUANTITIES]
                   + ["class"] + [f"U_{c}" for c in CHECKS]
                   + ["governing", "util", "not_checked"])
        w.writerows(rows)
    return len(rows)


def write_buckling_csv(path, run):
    cols = ["member", "design_case", "N_Ed_kN", "My_Ed_kNm", "Lcr_y_m", "Lcr_z_m", "Lw_m",
            "chi_y", "chi_z", "C1", "Mcr_kNm", "chi_LT", "U_Nb_y", "U_Nb_z", "U_LTB", "status"]

    def num(v, scale=1.0, f="{:.4f}"):
        return "" if v is None else f.format(v * scale)

    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for m in run.members:
            mid = m.member_id
            for name in run.cases:
                r = run.buckling.get((mid, name))
                if r is None:
                    continue
                v = r.values
                w.writerow([mid, name, num(v.get("N_Ed"), 1e-3), num(v.get("My_Ed"), 1e-3),
                            num(v.get("Lcr_y")), num(v.get("Lcr_z")), num(v.get("Lw")),
                            num(v.get("chi_y")), num(v.get("chi_z")), num(v.get("C1")),
                            num(v.get("Mcr"), 1e-3), num(v.get("chi_LT")),
                            fmt_util(r.util.get("Nb,y")), fmt_util(r.util.get("Nb,z")),
                            fmt_util(r.util.get("LTB")), r.status])
