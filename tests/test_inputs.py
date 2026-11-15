"""Reading the inputs (settings, member table, design cases, project file) and member chains."""

import os

import pytest

from strand7tools import settings as settings_module
from strand7tools.design_cases import DesignCaseError, load_design_cases
from strand7tools.member_table import MemberDef, parse_int_list, resolve_beams
from strand7tools.project import ProjectError, load_project
from strand7tools.strand7.members import MemberError, order_chain


def write(path, text):
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_settings_without_a_buckling_section_use_the_defaults(tmp_path, bk):
    with open(settings_module.SETTINGS_FILE, encoding="utf-8") as f:
        text = f.read()
    kept, skipping = [], False           # remove only the [buckling] and [buckling.*] tables
    for line in text.splitlines(keepends=True):
        header = line.strip()
        if header.startswith("["):
            skipping = header.strip("[]").split(".")[0].strip() == "buckling"
        if not skipping:
            kept.append(line)
    old = settings_module.load_settings(write(tmp_path / "settings.toml", "".join(kept)))
    assert old.buckling.alpha == bk.alpha
    assert old.buckling.flexural_rows == bk.flexural_rows
    assert old.buckling.ltb_rows == bk.ltb_rows


def test_grade_strength_by_thickness():
    grades = settings_module.load_grades()
    fy, fu = settings_module.strength(grades, "S355", 10.2)
    assert fy == 355 and fu > fy
    with pytest.raises(settings_module.SettingsError):
        settings_module.strength(grades, "S999", 10)


def test_member_table_parsing():
    assert parse_int_list("1-4;7") == [1, 2, 3, 4, 7]
    assert resolve_beams("ID=5;3", {1: 5, 2: 0, 3: 1, 4: 5}, 4) == [1, 4, 3]
    with pytest.raises(ValueError):
        resolve_beams("9", {}, 4)


def test_effective_lengths():
    d = MemberDef("1", [1], kcr_y=0.7, buckling_z=False, lw=2.5)
    assert d.effective_lengths(10.0) == {"Lcr_y": pytest.approx(7.0), "Lcr_z": None, "Lw": 2.5}


def test_order_chain_handles_reversed_beams():
    order, start, end = order_chain({6: (10, 11), 5: (12, 11), 3: (12, 13)})
    assert (start, end) == (10, 13)
    assert order == [(6, True), (5, False), (3, True)]


@pytest.mark.parametrize("ends", [
    {1: (1, 2), 2: (2, 3), 3: (2, 4)},       # branch
    {1: (1, 2), 2: (3, 4)},                  # gap
])
def test_order_chain_rejects_branches_and_gaps(ends):
    with pytest.raises(MemberError):
        order_chain(ends)


def test_design_cases(tmp_path):
    path = write(tmp_path / "dc.csv",
                 "# comment\ndesign_case,result_file,result_case,factor,comment\n"
                 "ULS,,1,1.35,G\nULS,,2,1.5,Q\nNL,run.NLA,1,1.0,\n")
    cases = load_design_cases(path)
    assert not cases["ULS"].is_direct and cases["NL"].is_direct
    assert cases["NL"].rows[0].result_file == os.path.join(str(tmp_path), "run.NLA")
    bad = write(tmp_path / "bad.csv", "design_case,result_file,result_case,factor\n"
                                      "NL,run.NLA,1,1.0\nNL,run.NLA,2,1.0\n")
    with pytest.raises(DesignCaseError):
        load_design_cases(bad)              # nonlinear results cannot be added together


def test_project_file(tmp_path):
    path = write(tmp_path / "p.toml", 'model = "m.st7"\nmembers = "mem.csv"\n'
                                      'design_cases = "dc.csv"\nrun_cases = ["ULS"]\n')
    p = load_project(path)
    assert p.model == os.path.join(str(tmp_path), "m.st7")
    assert p.result_file == os.path.join(str(tmp_path), "m.LSA")
    assert p.run_cases == ("ULS",) and p.stations_per_beam == 21
    assert p.settings == settings_module.SETTINGS_FILE
    with pytest.raises(ProjectError):
        load_project(write(tmp_path / "typo.toml", 'model = "m.st7"\nmembrs = "x"\n'))
