"""
Self-test for ec3_buckling.py. Needs no Strand7: run  python test_buckling.py

305x165x40 UB in S355, E = 210 GPa, G = 81 GPa, gamma_M1 = 1.0. Expected values come from
plain arithmetic below; the C1 values are the published ones for fork-supported beams
(EN 1993-1-1 ENV 1993-1-1 Annex F / NCCI SN003) and are NOT calculated by the code under test.
"""

import math
import os
import shutil
import sys
import tempfile
from types import SimpleNamespace

import settings as settings_module
from design_eval import StationEval, buckling_check
from ec3_buckling import (c1_from_diagram, check_member_buckling, flexural_curve, ltb_curve,
                          mcr_uniform)
from ec3_checks import ISection, check_station

E, G, GM1, FY = 210e9, 81e9, 1.0, 355e6
H, B, TW, TF, R = 0.3034, 0.165, 0.006, 0.0102, 0.009
A, IY, IZ = 5.133531e-03, 8.505546e-05, 7.643650e-06
IT = 14.90e-8
IW = IZ * (H - TF) ** 2 / 4
WPL_Y = 6.232880e-04

sec = ISection.from_dims("UB 305x165x40", H, B, TW, TF, R)
sec.E, sec.G = E, G
bk = settings_module.load_settings().buckling
failures = []


def expect(label, got, want, tol=2e-3):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    print(f"  {'ok  ' if ok else 'FAIL'} {label:<52} got {got:10.4f}   expected {want:10.4f}")
    if not ok:
        failures.append(label)


def expect_true(label, cond):
    print(f"  {'ok  ' if cond else 'FAIL'} {label}")
    if not cond:
        failures.append(label)


def stations(L, f_n, f_m, n=41):
    return [(L * i / (n - 1), f_n(i / (n - 1)), f_m(i / (n - 1))) for i in range(n)]


print("1. Mcr for uniform moment, closed form")
L = 6.0
want = (math.pi ** 2 * E * IZ / L ** 2) * math.sqrt(IW / IZ + L ** 2 * G * IT / (math.pi ** 2 * E * IZ))
expect("Mcr,uniform [kNm]", mcr_uniform(E, G, IZ, IT, IW, L) / 1e3, want / 1e3)

print("2. C1 from the moment diagram against published values (L = 6 m)")
n = 200
xs = [L * i / n for i in range(n + 1)]


def c1(f):
    return c1_from_diagram(xs, [f(x / L) for x in xs], L, E, G, IZ, IT, IW)


expect("uniform moment: 1.000", c1(lambda s: 1.0), 1.000, tol=1e-3)
expect("uniform load, simply supported: 1.132", c1(lambda s: 4 * s * (1 - s)), 1.132, tol=0.01)
expect("central point load: 1.365", c1(lambda s: 2 * s if s < .5 else 2 * (1 - s)), 1.365, tol=0.01)
expect("end moments psi = 0.5: 1.323", c1(lambda s: 1 - 0.5 * s), 1.323, tol=0.01)
expect("end moments psi = 0 (table 1.879)", c1(lambda s: 1 - s), 1.879, tol=0.03)
expect("end moments psi = -1: 2.752", c1(lambda s: 1 - 2 * s), 2.752, tol=0.02)
expect("sign of the diagram does not matter", c1(lambda s: -4 * s * (1 - s)), c1(lambda s: 4 * s * (1 - s)))

print("3. Flexural buckling about z-z: 400 kN, Lcr,z = 4 m")
Lcr = 4.0
ncr = math.pi ** 2 * E * IZ / Lcr ** 2
lam = math.sqrt(A * FY / ncr)
phi = 0.5 * (1 + 0.34 * (lam - 0.2) + lam ** 2)             # h/b = 1.84 > 1.2, tf < 40: z-z curve b
chi = 1 / (phi + math.sqrt(phi ** 2 - lam ** 2))
nb = chi * A * FY
st = stations(6.0, lambda s: -400e3, lambda s: 0.0)
r = check_member_buckling(sec, FY, E, G, GM1, bk, st, 6.0, 1, None, Lcr, None)
expect("lambda_z", r.values["lam_z"], lam)
expect("chi_z", r.values["chi_z"], chi)
expect("U(Nb,z) = 400 / Nb,z,Rd", r.util["Nb,z"], 400e3 / nb)
expect_true("y-y switched off: no Nb,y result", "Nb,y" not in r.util)
r = check_member_buckling(sec, FY, E, G, GM1, bk, st, 6.0, 1, 6.0, 6.0, None)
ncr_y = math.pi ** 2 * E * IY / 36.0
lam_y = math.sqrt(A * FY / ncr_y)
phi_y = 0.5 * (1 + 0.21 * (lam_y - 0.2) + lam_y ** 2)       # y-y curve a
chi_y = 1 / (phi_y + math.sqrt(phi_y ** 2 - lam_y ** 2))
expect("chi_y (curve a)", r.values["chi_y"], chi_y)

print("4. Tension only: no buckling check")
st = stations(6.0, lambda s: 300e3, lambda s: 0.0)
r = check_member_buckling(sec, FY, E, G, GM1, bk, st, 6.0, 1, 6.0, 6.0, None)
expect("U(Nb,y) = U(Nb,z) = 0", r.util["Nb,y"] + r.util["Nb,z"], 0.0)

print("5. LTB, uniform moment 100 kNm, Lw = L = 6 m")
st = stations(6.0, lambda s: 0.0, lambda s: 100e3)
r = check_member_buckling(sec, FY, E, G, GM1, bk, st, 6.0, 1, None, None, 6.0)
mcr = (math.pi ** 2 * E * IZ / 36.0) * math.sqrt(IW / IZ + 36.0 * G * IT / (math.pi ** 2 * E * IZ))
lam_lt = math.sqrt(WPL_Y * FY / mcr)
phi_lt = 0.5 * (1 + 0.21 * (lam_lt - 0.2) + lam_lt ** 2)      # h/b = 1.84 <= 2: curve a
chi_lt = min(1.0, 1 / (phi_lt + math.sqrt(phi_lt ** 2 - lam_lt ** 2)), 1 / lam_lt ** 2)
expect("Mcr [kNm]", r.values["Mcr"] / 1e3, mcr / 1e3, tol=2e-3)
expect("lambda_LT", r.values["lam_LT"], lam_lt)
expect("chi_LT", r.values["chi_LT"], chi_lt)
expect("U(LTB) = 100 / Mb,Rd", r.util["LTB"], 100e3 / (chi_lt * WPL_Y * FY))

print("6. LTB with a UDL moment diagram (C1 = 1.132) and the Lw override")
st = stations(7.0, lambda s: 0.0, lambda s: 49e3 * 4 * s * (1 - s))
r = check_member_buckling(sec, FY, E, G, GM1, bk, st, 7.0, 1, None, None, 7.0)
mcr_u = (math.pi ** 2 * E * IZ / 49.0) * math.sqrt(IW / IZ + 49.0 * G * IT / (math.pi ** 2 * E * IZ))
expect("Mcr = 1.132 x Mcr,uniform [kNm]", r.values["Mcr"] / 1e3, 1.132 * mcr_u / 1e3, tol=0.01)
r2 = check_member_buckling(sec, FY, E, G, GM1, bk, st, 7.0, 1, None, None, 3.5)
expect_true("shorter Lw raises Mcr and lowers the utilisation", r2.values["Mcr"] > r.values["Mcr"]
            and r2.util["LTB"] < r.util["LTB"])
expect_true("a note warns that Lw differs from the member length",
            any("differs" in nt for nt in r2.notes))

print("7. Buckling curve tables")
expect_true("305x165x40, S355: y-y a, z-z b", (flexural_curve(bk, sec, 355, "y"),
                                                flexural_curve(bk, sec, 355, "z")) == ("a", "b"))
expect_true("305x165x40, S460: a0, a0", (flexural_curve(bk, sec, 460, "y"),
                                         flexural_curve(bk, sec, 460, "z")) == ("a0", "a0"))
stocky = ISection.from_dims("HE-like", 0.200, 0.200, 0.009, 0.015, 0.018)       # h/b = 1.0
expect_true("h/b <= 1.2, tf <= 100: b, c", (flexural_curve(bk, stocky, 355, "y"),
                                            flexural_curve(bk, stocky, 355, "z")) == ("b", "c"))
expect_true("LTB: h/b <= 2 curve a; h/b > 2 curve b", ltb_curve(bk, sec) == "a" and
            ltb_curve(bk, ISection.from_dims("deep", 0.6, 0.2, 0.01, 0.016, 0.012)) == "b")

print("8. Class 4 is not checked; through design_eval")
r = check_member_buckling(sec, FY, E, G, GM1, bk, stations(6.0, lambda s: -100e3, lambda s: 0.0),
                          6.0, 4, 6.0, 6.0, 6.0)
expect_true("class 4: not checked", bool(r.not_checked) and not r.util)

L = 7.0
pts = [(L * i / 20, {"N": 0.0, "Vy": 0.0, "Vz": 0.0, "T": 0.0,
                      "My": 49.0 * 4 * (i / 20) * (1 - i / 20), "Mz": 0.0}) for i in range(21)]
evals = [StationEval(1, x, f, check_station(sec, FY, 1.0, {q: v * 1e3 for q, v in f.items()}))
         for x, f in pts]
r = buckling_check(sec, FY, E, G, GM1, bk, evals, L, {"Lcr_y": L, "Lcr_z": L, "Lw": L}, trace=True)
expect("via design_eval: |My|,Ed = 49 kNm", r.values["My_Ed"] / 1e3, 49.0)
expect_true("worked calculation text produced", len(r.trace) > 10)

print("8b. Compression with bending needs 6.3.3: reported as not fully checked")
NPL = A * FY


def interaction_flagged(n, my, mz=0.0, lengths=(6.0, 6.0, 6.0)):
    st = stations(6.0, lambda s: n, lambda s: my * 4 * s * (1 - s))
    r = check_member_buckling(sec, FY, E, G, GM1, bk, st, 6.0, 1, *lengths, mz_ed=mz)
    return any("6.3.3" in reason for reason in r.not_checked), r


flagged, r = interaction_flagged(-100e3, 50e3)
expect_true("100 kN compression + 50 kNm My: flagged", flagged)
expect_true("... and the N and LTB results are still given", "Nb,z" in r.util and "LTB" in r.util)
expect_true("100 kN compression + 5 kNm Mz only: flagged", interaction_flagged(-100e3, 0.0, 5e3)[0])
expect_true("100 kN TENSION + 50 kNm My: not flagged", not interaction_flagged(100e3, 50e3)[0])
expect_true("100 kN compression, no bending: not flagged", not interaction_flagged(-100e3, 0.0)[0])
expect_true("50 kNm My, no axial force: not flagged", not interaction_flagged(0.0, 50e3)[0])
expect_true("solver noise (N = 0.5 kN < 0.1% Npl) + 50 kNm: not flagged",
            0.5e3 < 1e-3 * NPL and not interaction_flagged(-0.5e3, 50e3)[0])
expect_true("all buckling checks switched off: not flagged",
            not interaction_flagged(-100e3, 50e3, lengths=(None, None, None))[0])

print("9. settings.toml without a [buckling] section still loads (defaults)")
tmp = tempfile.mkdtemp()
try:
    with open(settings_module.SETTINGS_FILE, encoding="utf-8") as f:
        text = f.read()
    # Remove only the [buckling] block (and any [buckling.*] tables), wherever it sits.
    kept, skipping = [], False
    for line in text.splitlines(keepends=True):
        header = line.strip()
        if header.startswith("["):
            skipping = header.strip("[]").split(".")[0].strip() == "buckling"
        if not skipping:
            kept.append(line)
    path = os.path.join(tmp, "settings.toml")
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(kept))
    s_old = settings_module.load_settings(path)
    expect_true("defaults equal the shipped values", s_old.buckling.alpha == bk.alpha
                and s_old.buckling.flexural_rows == bk.flexural_rows
                and s_old.buckling.ltb_rows == bk.ltb_rows)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
if failures:
    print(f"{len(failures)} FAILED: " + "; ".join(failures))
    sys.exit(1)
print("All checks passed.")
