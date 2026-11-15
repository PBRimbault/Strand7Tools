"""
Run the EC3 checks for a project: read the model and results (read-only), then check every
member for every design case. Returns everything in a RunResult; nothing is printed here
(report.py and worked.py present the results, and a GUI can use the same RunResult).

For every member and design case:
  - the beam results are read from the .LSA result file(s) and converted to member axes
    (a design case that adds cases together is evaluated on a fixed set of stations along
    each beam; a case taken directly from Strand7 uses Strand7's own stations, which include
    the maxima and minima)
  - the cross-section checks (6.2) are made at every station
  - the member buckling checks (6.3.1, 6.3.2) are made once, from the worst values
"""

import os
from dataclasses import dataclass, field

from .design_cases import DesignCaseError, file_warnings, load_design_cases
from .ec3.design import (CombineError, buckling_check, combine, evaluate_member, status_rank,
                         summarise)
from .ec3.sections import Steel
from .settings import SettingsError, load_grades, load_settings, strength
from .strand7.api import model_session, total_beams
from .strand7.forces import assemble, force_scale, result_file
from .strand7.members import load_members
from .strand7.sections import SectionError, make_section, read_section

E_G_TOLERANCE = 0.01        # warn if the model's E or G differs from E_ec3 / G_ec3 by more than this


@dataclass
class MemberInfo:
    section: object             # ec3.sections.Section
    steel: object               # ec3.sections.Steel
    warnings: list = field(default_factory=list)


@dataclass
class RunResult:
    project: object
    settings: object
    cases: dict                                     # {name: DesignCase}, in run order
    members: list = field(default_factory=list)     # members that could be checked
    info: dict = field(default_factory=dict)        # member id -> MemberInfo
    errors: list = field(default_factory=list)      # (member id, message): not built / not checked
    case_errors: list = field(default_factory=list)  # (member id, design case, message)
    file_warnings: list = field(default_factory=list)
    material_warnings: list = field(default_factory=list)
    evals: dict = field(default_factory=dict)       # (member id, case) -> [StationEval]
    summaries: dict = field(default_factory=dict)   # member id -> [Summary], in case order
    buckling: dict = field(default_factory=dict)    # (member id, case) -> BucklingResult
    unassigned_beams: list = field(default_factory=list)

    def governing_summary(self, member_id):
        """The cross-section summary to show for a member (see design.status_rank)."""
        lst = self.summaries.get(member_id)
        if not lst:
            return None
        return min(lst, key=lambda s: status_rank(s.status, s.governing_util))

    def governing_buckling(self, member_id):
        """(design case, BucklingResult) to show for a member, or None."""
        rows = [(name, self.buckling[(member_id, name)]) for name in self.cases
                if (member_id, name) in self.buckling]
        if not rows:
            return None
        return min(rows, key=lambda r: status_rank(r[1].status, r[1].governing[1]))


# ----------------------------------------------------------------------
# Sections and materials
# ----------------------------------------------------------------------

def member_section(uid, m, settings, grades):
    """MemberInfo for a member (section from the model, steel from the grade table)."""
    data = read_section(uid, m.chain[0].prop)
    sec, warnings = make_section(data, settings.eta)
    if settings.material_source == "model":
        E, G = data["E"], data["G"]
        if not (150e9 <= E <= 250e9) or not (0.30 <= G / E <= 0.55):
            raise SectionError(f"the beam material has E = {E:g}, G = {G:g}: these should "
                               f"be about 210e9 and 81e9 Pa. Check the "
                               f"material, or set source = \"ec3\" in settings.toml")
    else:
        E, G = settings.E_ec3 * 1e6, settings.G_ec3 * 1e6
    grade = m.definition.grade or settings.default_grade
    fy, fu = strength(grades, grade, sec.t_max * 1e3)
    return MemberInfo(sec, Steel(grade, fy * 1e6, fu * 1e6, E, G), warnings)


def material_warnings(members, info, settings):
    """
    With material source "model", the Strand7 E and G are used on purpose (the same stiffness
    as the analysis). Say so, once per beam property, whenever they differ from the code values
    (E_ec3, G_ec3 in settings.toml) by more than E_G_TOLERANCE, so a wrong material is noticed.
    """
    if settings.material_source != "model":
        return []
    users = {}
    for m in members:
        users.setdefault(m.chain[0].prop, []).append(m.member_id)
    out = []
    for prop, mids in users.items():
        steel = info[mids[0]].steel
        for sym, value, ref_mpa in (("E", steel.E, settings.E_ec3), ("G", steel.G, settings.G_ec3)):
            ref = ref_mpa * 1e6
            if abs(value / ref - 1) > E_G_TOLERANCE:
                out.append(f"{sym} = {value / 1e9:.1f} GPa in the Strand7 material of property "
                           f"{prop} (members {', '.join(mids)}) differs from {ref / 1e9:.1f} GPa "
                           f"(EN 1993-1-1 3.2.6, {sym}_ec3 in settings.toml) by "
                           f"{(value / ref - 1) * 100:+.1f}%. The model value is used; update the "
                           f"Strand7 material if this is not intended.")
    return out


# ----------------------------------------------------------------------
# Getting the forces
# ----------------------------------------------------------------------

def _key(path):
    return os.path.normcase(os.path.abspath(path))


def plan_extraction(cases, default_path):
    """{file key: (path, {(result case, 'native' or 'grid')})} - what has to be read from where."""
    needs = {}
    for dc in cases.values():
        mode = "native" if dc.is_direct else "grid"
        for row in dc.rows:
            path = row.result_file or default_path
            needs.setdefault(_key(path), (path, set()))[1].add((row.result_case, mode))
    return needs


def extract(uid, members, needs, scale, combinations, grid):
    """Read every needed result case for every member. Returns {(file, case, mode, member): by_beam}."""
    data = {}
    for key, (path, wanted) in needs.items():
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Result file not found: {path}")
        if os.path.splitext(path)[1].lower() != ".lsa":
            raise DesignCaseError(f"{path}: only linear static (.LSA) result files are supported "
                                  f"so far")
        with result_file(uid, path, combinations) as (n_prim, n_sec):
            total = n_prim + n_sec
            for case, mode in sorted(wanted):
                if not 1 <= case <= total:
                    raise DesignCaseError(f"result case {case} does not exist in {path} "
                                          f"(it has {total})")
                for m in members:
                    data[(key, case, mode, m.member_id)] = assemble(
                        uid, m, case, scale, grid if mode == "grid" else None)
    return data


def design_case_forces(dc, member_id, data, default_path):
    """The member force diagram for one design case."""
    def key(row):
        return _key(row.result_file or default_path)

    if dc.is_direct:
        row = dc.rows[0]
        return data[(key(row), row.result_case, "native", member_id)]
    return combine([(row.factor, data[(key(row), row.result_case, "grid", member_id)])
                    for row in dc.rows])


# ----------------------------------------------------------------------
# The run
# ----------------------------------------------------------------------

def select_cases(all_cases, names):
    if not names:
        return dict(all_cases)
    missing = [n for n in names if n not in all_cases]
    if missing:
        raise DesignCaseError(f"design case(s) not in the design case file: {', '.join(missing)}")
    return {n: all_cases[n] for n in names}


def run_checks(project):
    """Everything for one project. Raises on problems that stop the whole run."""
    settings = load_settings(project.settings)
    grades = load_grades(project.steel_grades)
    cases = select_cases(load_design_cases(project.design_cases), project.run_cases)
    n = project.stations_per_beam
    grid = [i / (n - 1) for i in range(n)]
    run = RunResult(project, settings, cases, file_warnings=file_warnings(cases))

    with model_session(project.model) as uid:
        scale = force_scale(uid)
        if scale is None:
            raise RuntimeError("the checks need a model in metres with forces in N, kN or MN "
                               "(they are converted to kN and kNm); change the model units")
        members, errors = load_members(uid, project.members)
        run.errors.extend((mid, f"ERROR - {msg}; skipped") for mid, msg in errors)
        n_beams = total_beams(uid)
        for m in members:
            try:
                run.info[m.member_id] = member_section(uid, m, settings, grades)
            except (SectionError, SettingsError, RuntimeError) as e:
                run.errors.append((m.member_id, f"NOT CHECKED - {e}"))
                continue
            run.members.append(m)
        data = extract(uid, run.members, plan_extraction(cases, project.result_file), scale,
                       project.combinations, grid)
    # ---- the Strand7 model is closed from here on; everything below is calculation ----

    run.material_warnings = material_warnings(run.members, run.info, settings)
    for m in run.members:
        mid = m.member_id
        mi = run.info[mid]
        run.summaries[mid] = []
        for name, dc in cases.items():
            try:
                by_beam = design_case_forces(dc, mid, data, project.result_file)
            except CombineError as e:
                run.case_errors.append((mid, name, str(e)))
                continue
            evals = evaluate_member(mi.section, mi.steel, settings.gamma_M0, by_beam)
            run.evals[(mid, name)] = evals
            run.summaries[mid].append(summarise(mid, name, evals))
            run.buckling[(mid, name)] = buckling_check(
                mi.section, mi.steel, settings.gamma_M1, settings.buckling, evals, m.length,
                m.effective_lengths())

    in_members = {c.beam for m in run.members for c in m.chain}
    run.unassigned_beams = [b for b in range(1, n_beams + 1) if b not in in_members]
    return run
