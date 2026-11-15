"""The worked calculations are produced from the results, and show the same numbers."""

from strand7tools.ec3.buckling import check_member_buckling
from strand7tools.ec3.cross_section import check_cross_section
from strand7tools.worked import buckling_lines, cross_section_lines


def forces(**kw):
    return {q: kw.get(q, 0.0) for q in ("N", "Vy", "Vz", "T", "My", "Mz")}


def test_cross_section_worked_calculation(sec, steel):
    r = check_cross_section(sec, steel, 1.0, forces(N=600e3, Vz=60e3, My=80e3))
    text = "\n".join(cross_section_lines(r))
    assert "1. CLASSIFICATION" in text and "6. BENDING AND AXIAL FORCE (6.2.9)" in text
    name, value = r.governing
    assert text.splitlines()[-1].startswith(f"GOVERNING: {name} = {value:.3f}   OK")


def test_cross_section_worked_calculation_class_3_and_4(sec, steel):
    text = "\n".join(cross_section_lines(
        check_cross_section(sec, steel, 1.0, forces(N=-500e3, My=100e3))))
    assert "class 3: elastic check" in text
    lines = cross_section_lines(check_cross_section(sec, steel, 1.0, forces(N=-1500e3)))
    assert "NOT CHECKED: class 4" in lines[-1]


def test_buckling_worked_calculation(sec, steel, bk):
    st = [(7.0 * i / 20, -50e3, 49e3 * 4 * (i / 20) * (1 - i / 20)) for i in range(21)]
    r = check_member_buckling(sec, steel, 1.0, bk, st, 7.0, 1, 7.0, 7.0, 7.0)
    lines = buckling_lines(r)
    text = "\n".join(lines)
    assert f"C1 from the moment diagram = {r.values['C1']:.3f}" in text
    assert "NOT CHECKED - the interaction (6.61)/(6.62)" in text
    assert lines[-1].startswith("GOVERNING:") and "INCOMPLETE" in lines[-1]
