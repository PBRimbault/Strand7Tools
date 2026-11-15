"""
Cross-section checks (EC3 6.2), strand7tools/ec3/cross_section.py.

Each expected value is calculated here by plain arithmetic from the section numbers in
conftest.py (S355, gamma_M0 = 1.0), written independently of the code under test.
"""

import math

import pytest

from conftest import A, B, EPS, FY, HW, IT, IY, R, TF, TW, WEL_Y, WPL_Y, WPL_Z, close
from strand7tools.ec3.cross_section import check_cross_section

NPL = A * FY
MPL_Y, MPL_Z = WPL_Y * FY, WPL_Z * FY
AV_Z = A - 2 * B * TF + (TW + 2 * R) * TF
VPL_Z = AV_Z * FY / math.sqrt(3)
TAU_RD = FY / math.sqrt(3)
T_HALF = 0.5 * TAU_RD * IT / TF          # torque giving tau_t = half the shear stress limit


@pytest.fixture
def run(sec, steel):
    def _run(N=0.0, Vy=0.0, Vz=0.0, T=0.0, My=0.0, Mz=0.0):
        return check_cross_section(sec, steel, 1.0,
                                   {"N": N, "Vy": Vy, "Vz": Vz, "T": T, "My": My, "Mz": Mz})
    return _run


def test_pure_major_axis_bending_class_1(run):
    r = run(My=100e3)
    assert r.cls == 1
    assert close(r.util["My"], 100e3 / MPL_Y)
    assert close(r.util["N+M"], 100e3 / MPL_Y)


def test_small_compression_is_neglected_class_2(run):
    """200 kN + My 150 kNm + Vz 60 kN: class 2, N below both 6.2.9.1(4) limits."""
    r = run(N=-200e3, Vz=60e3, My=150e3)
    assert r.cls == 2
    assert close(r.util["N"], 200e3 / NPL)
    assert close(r.util["Vz"], 60e3 / VPL_Z)
    assert 200e3 <= min(0.25 * NPL, 0.5 * HW * TW * FY)
    assert close(r.util["N+M"], 150e3 / MPL_Y)


def test_tension_reduces_mpl_y(run):
    r = run(N=600e3, My=80e3)
    n = 600e3 / NPL
    a = min((A - 2 * B * TF) / A, 0.5)
    assert r.cls == 1
    assert close(r.util["N+M"], 80e3 / (MPL_Y * (1 - n) / (1 - 0.5 * a)))


def test_biaxial_bending(run):
    r = run(N=-100e3, My=100e3, Mz=10e3)
    n = 100e3 / NPL
    ry, rz = 100e3 / MPL_Y, 10e3 / MPL_Z                  # n <= a and N small: no reduction
    assert r.cls == 1
    assert close(r.util["N+M"], max(ry, rz, ry ** 2 + rz ** max(5 * n, 1.0)))


def test_uniaxial_overload_reports_the_plain_ratio(run):
    assert close(run(My=250e3).util["N+M"], 250e3 / MPL_Y)


def test_high_shear_reduces_the_moment_resistance(run):
    r = run(Vz=300e3, My=100e3)
    rho = (2 * 300e3 / VPL_Z - 1) ** 2
    mv = (WPL_Y - rho * (HW * TW) ** 2 / (4 * TW)) * FY
    assert close(r.util["Vz"], 300e3 / VPL_Z)
    assert close(r.util["My"], 100e3 / mv)


def test_class_4_is_not_checked(run):
    r = run(N=-1500e3)
    assert r.cls == 4
    assert r.not_checked and not r.util


def test_torsion_constant(sec):
    assert close(sec.props.It * 1e8, IT * 1e8, tol=0.005)


def test_solver_noise_torque_gives_no_note(run):
    r = run(T=1.0, My=10e3)
    assert not r.notes and not r.not_checked and r.util["T"] < 1e-3


def test_st_venant_torsion(run):
    r = run(T=2e3, My=10e3)
    assert close(r.util["T"], 2e3 * TF / IT / TAU_RD)
    assert any("St Venant" in n for n in r.notes) and not r.not_checked


def test_shear_resistance_reduced_for_torsion(run):
    red = math.sqrt(1 - 0.5 / 1.25)
    r = run(T=T_HALF, Vz=100e3, My=10e3)
    assert close(r.util["T"], 0.5)
    assert close(r.util["Vz"], 100e3 / (VPL_Z * red))
    # above half of the REDUCED Vpl,z: moment reduced with rho from Vpl,T,Rd
    r = run(T=T_HALF, Vz=300e3, My=100e3)
    rho = (2 * 300e3 / (VPL_Z * red) - 1) ** 2
    mv = (WPL_Y - rho * (HW * TW) ** 2 / (4 * TW)) * FY
    assert close(r.util["My"], 100e3 / mv)


def test_torsion_alone_can_use_up_the_shear_resistance(run):
    r = run(T=2.5 * T_HALF, Vz=10e3)                  # tau_t above 1.25 tau_Rd
    assert r.util["T"] > 1 and math.isinf(r.util["Vz"])


def test_class_3_elastic_check(run):
    r = run(N=-500e3, My=100e3)
    sigma = max(abs(-500e3 / A + sy * 100e3 / WEL_Y) for sy in (1, -1))
    assert r.cls == 3
    assert close(r.util["N+M"], sigma / FY)
    assert close(r.util["My"], 100e3 / (WEL_Y * FY))


def test_classification_numbers(run):
    r = run(N=-500e3, My=100e3)
    flange, web = r.classification.part("Flange outstand"), r.classification.part("Web")
    d = HW - 2 * R
    assert close(flange.ratio, (B - TW - 2 * R) / 2 / TF)
    assert close(web.ratio, d / TW)
    assert close(web.alpha, 0.5 * (1 + 500e3 / (d * TW * FY)))
    sig_mean, bend = 500e3 / A, 100e3 * (d / 2) / IY
    psi = (sig_mean - bend) / (sig_mean + bend)
    assert close(web.psi, psi)
    assert close(web.limits[2], 42 * EPS / (0.67 + 0.33 * psi))
