"""
Member buckling (EC3 6.3), strand7tools/ec3/buckling.py.

305x165x40 UB in S355, E = 210 GPa, G = 81 GPa, gamma_M1 = 1.0. Expected values come from
plain arithmetic below; the C1 values are the published ones for fork-supported beams
(ENV 1993-1-1 Annex F / NCCI SN003) and are NOT calculated by the code under test.
"""

import math

import pytest

from conftest import A, E, FY, G, IT, IW, IY, IZ, WPL_Y, close
from strand7tools.ec3.buckling import c1_from_diagram, check_member_buckling, mcr_uniform
from strand7tools.ec3.cross_section import check_cross_section
from strand7tools.ec3.design import StationEval, buckling_check
from strand7tools.ec3.sections import ISection

GM1 = 1.0
NPL = A * FY


def stations(L, f_n, f_m, n=41):
    return [(L * i / (n - 1), f_n(i / (n - 1)), f_m(i / (n - 1))) for i in range(n)]


@pytest.fixture
def check(sec, steel, bk):
    def _check(st, length, cls, Lcr_y, Lcr_z, Lw, mz_ed=0.0):
        return check_member_buckling(sec, steel, GM1, bk, st, length, cls, Lcr_y, Lcr_z, Lw,
                                     mz_ed=mz_ed)
    return _check


def test_mcr_uniform_closed_form():
    L = 6.0
    want = (math.pi ** 2 * E * IZ / L ** 2) * math.sqrt(IW / IZ + L ** 2 * G * IT
                                                         / (math.pi ** 2 * E * IZ))
    assert close(mcr_uniform(E, G, IZ, IT, IW, L) / 1e3, want / 1e3)


L6 = 6.0
XS = [L6 * i / 200 for i in range(201)]


def c1(f):
    return c1_from_diagram(XS, [f(x / L6) for x in XS], L6, E, G, IZ, IT, IW)


@pytest.mark.parametrize("label, diagram, want, tol", [
    ("uniform moment", lambda s: 1.0, 1.000, 1e-3),
    ("uniform load, simply supported", lambda s: 4 * s * (1 - s), 1.132, 0.01),
    ("central point load", lambda s: 2 * s if s < .5 else 2 * (1 - s), 1.365, 0.01),
    ("end moments psi = 0.5", lambda s: 1 - 0.5 * s, 1.323, 0.01),
    ("end moments psi = 0 (table 1.879)", lambda s: 1 - s, 1.879, 0.03),
    ("end moments psi = -1", lambda s: 1 - 2 * s, 2.752, 0.02),
])
def test_c1_against_published_values(label, diagram, want, tol):
    assert close(c1(diagram), want, tol), label


def test_c1_sign_of_the_diagram_does_not_matter():
    assert close(c1(lambda s: -4 * s * (1 - s)), c1(lambda s: 4 * s * (1 - s)))


def test_flexural_buckling_z(check):
    Lcr = 4.0
    ncr = math.pi ** 2 * E * IZ / Lcr ** 2
    lam = math.sqrt(A * FY / ncr)
    phi = 0.5 * (1 + 0.34 * (lam - 0.2) + lam ** 2)          # h/b = 1.84 > 1.2, tf < 40: curve b
    chi = 1 / (phi + math.sqrt(phi ** 2 - lam ** 2))
    r = check(stations(6.0, lambda s: -400e3, lambda s: 0.0), 6.0, 1, None, Lcr, None)
    assert close(r.values["lam_z"], lam)
    assert close(r.values["chi_z"], chi)
    assert close(r.util["Nb,z"], 400e3 / (chi * A * FY))
    assert "Nb,y" not in r.util and r.flexural_state["y"] == "off"


def test_flexural_buckling_y_curve_a(check):
    r = check(stations(6.0, lambda s: -400e3, lambda s: 0.0), 6.0, 1, 6.0, 6.0, None)
    lam_y = math.sqrt(A * FY / (math.pi ** 2 * E * IY / 36.0))
    phi_y = 0.5 * (1 + 0.21 * (lam_y - 0.2) + lam_y ** 2)
    assert close(r.values["chi_y"], 1 / (phi_y + math.sqrt(phi_y ** 2 - lam_y ** 2)))


def test_tension_needs_no_buckling_check(check):
    r = check(stations(6.0, lambda s: 300e3, lambda s: 0.0), 6.0, 1, 6.0, 6.0, None)
    assert r.util["Nb,y"] == 0.0 and r.util["Nb,z"] == 0.0


def test_ltb_uniform_moment(check):
    r = check(stations(6.0, lambda s: 0.0, lambda s: 100e3), 6.0, 1, None, None, 6.0)
    mcr = (math.pi ** 2 * E * IZ / 36.0) * math.sqrt(IW / IZ + 36.0 * G * IT / (math.pi ** 2 * E * IZ))
    lam_lt = math.sqrt(WPL_Y * FY / mcr)
    phi_lt = 0.5 * (1 + 0.21 * (lam_lt - 0.2) + lam_lt ** 2)      # h/b = 1.84 <= 2: curve a
    chi_lt = min(1.0, 1 / (phi_lt + math.sqrt(phi_lt ** 2 - lam_lt ** 2)), 1 / lam_lt ** 2)
    assert close(r.values["Mcr"] / 1e3, mcr / 1e3)
    assert close(r.values["lam_LT"], lam_lt)
    assert close(r.values["chi_LT"], chi_lt)
    assert close(r.util["LTB"], 100e3 / (chi_lt * WPL_Y * FY))


def test_ltb_udl_diagram_and_lw_override(check):
    st = stations(7.0, lambda s: 0.0, lambda s: 49e3 * 4 * s * (1 - s))
    r = check(st, 7.0, 1, None, None, 7.0)
    mcr_u = (math.pi ** 2 * E * IZ / 49.0) * math.sqrt(IW / IZ + 49.0 * G * IT / (math.pi ** 2 * E * IZ))
    assert close(r.values["Mcr"] / 1e3, 1.132 * mcr_u / 1e3, tol=0.01)
    r2 = check(st, 7.0, 1, None, None, 3.5)
    assert r2.values["Mcr"] > r.values["Mcr"] and r2.util["LTB"] < r.util["LTB"]
    assert any("differs" in nt for nt in r2.notes)


def test_buckling_curves(sec, bk):
    assert (sec.flexural_curve(bk, 355, "y"), sec.flexural_curve(bk, 355, "z")) == ("a", "b")
    assert (sec.flexural_curve(bk, 460, "y"), sec.flexural_curve(bk, 460, "z")) == ("a0", "a0")
    stocky = ISection("HE-like", 0.200, 0.200, 0.009, 0.015, 0.018)       # h/b = 1.0
    assert (stocky.flexural_curve(bk, 355, "y"), stocky.flexural_curve(bk, 355, "z")) == ("b", "c")
    assert sec.ltb_curve(bk) == "a"
    assert ISection("deep", 0.6, 0.2, 0.01, 0.016, 0.012).ltb_curve(bk) == "b"


def test_class_4_is_not_checked(check):
    r = check(stations(6.0, lambda s: -100e3, lambda s: 0.0), 6.0, 4, 6.0, 6.0, 6.0)
    assert r.not_checked and not r.util
    assert r.status == "NOT CHECKED"


def test_through_design_buckling_check(sec, steel, bk):
    L = 7.0
    pts = [(L * i / 20, {"N": 0.0, "Vy": 0.0, "Vz": 0.0, "T": 0.0,
                         "My": 49.0 * 4 * (i / 20) * (1 - i / 20), "Mz": 0.0}) for i in range(21)]
    evals = [StationEval(1, x, f, check_cross_section(sec, steel, 1.0,
                                                      {q: v * 1e3 for q, v in f.items()}))
             for x, f in pts]
    r = buckling_check(sec, steel, GM1, bk, evals, L, {"Lcr_y": L, "Lcr_z": L, "Lw": L})
    assert close(r.values["My_Ed"] / 1e3, 49.0)


# --- 6.3.3 is not implemented: compression with bending must be reported as not fully checked

def interaction(check, n, my, mz=0.0, lengths=(6.0, 6.0, 6.0)):
    r = check(stations(6.0, lambda s: n, lambda s: my * 4 * s * (1 - s)), 6.0, 1, *lengths,
              mz_ed=mz)
    return any("6.3.3" in reason for reason in r.not_checked), r


def test_compression_with_bending_is_flagged(check):
    flagged, r = interaction(check, -100e3, 50e3)
    assert flagged and r.interaction_state == "flagged"
    assert "Nb,z" in r.util and "LTB" in r.util          # the separate results are still given
    assert r.status == "INCOMPLETE"


@pytest.mark.parametrize("label, n, my, mz, lengths, flagged", [
    ("compression + minor-axis bending only", -100e3, 0.0, 5e3, (6.0, 6.0, 6.0), True),
    ("tension + bending", 100e3, 50e3, 0.0, (6.0, 6.0, 6.0), False),
    ("compression, no bending", -100e3, 0.0, 0.0, (6.0, 6.0, 6.0), False),
    ("bending, no axial force", 0.0, 50e3, 0.0, (6.0, 6.0, 6.0), False),
    ("solver noise: 0.5 kN (< 0.1% Npl) + bending", -0.5e3, 50e3, 0.0, (6.0, 6.0, 6.0), False),
    ("all buckling checks switched off", -100e3, 50e3, 0.0, (None, None, None), False),
])
def test_interaction_flag(check, label, n, my, mz, lengths, flagged):
    assert interaction(check, n, my, mz, lengths)[0] == flagged, label
