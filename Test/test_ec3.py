"""
Self-test for ec3_checks.py. Needs no Strand7: run  python test_ec3.py

Each expected value is calculated here by plain arithmetic from the section numbers
below (305x165x40 UB, BS EN 10365 dimensions, S355, gamma_M0 = 1.0), written independently
of ec3_checks.py. The section numbers come from a brute-force numerical integration of
the real section shape, root fillets included.
"""

import math
import sys

from ec3_checks import ISection, check_station

# --- section, from numerical integration of the true shape (m, m2, m3) -----------------
H, B, TW, TF, R = 0.3034, 0.165, 0.006, 0.0102, 0.009
A = 5.133531e-03
WPL_Y, WPL_Z = 6.232880e-04, 1.417429e-04
WEL_Y, WEL_Z = 8.505546e-05 / (H / 2), 7.643650e-06 / (B / 2)
HW = H - 2 * TF
FY = 355e6
EPS = math.sqrt(235 / 355)

NPL = A * FY                         # N
MPL_Y, MPL_Z = WPL_Y * FY, WPL_Z * FY
AV_Z = A - 2 * B * TF + (TW + 2 * R) * TF
VPL_Z = AV_Z * FY / math.sqrt(3)

sec = ISection.from_dims("UB 305x165x40", H, B, TW, TF, R, eta=1.0)
failures = []


def run(label, N=0.0, Vy=0.0, Vz=0.0, T=0.0, My=0.0, Mz=0.0):
    return check_station(sec, FY, 1.0, {"N": N, "Vy": Vy, "Vz": Vz, "T": T, "My": My, "Mz": Mz})


def expect(label, got, want, tol=2e-3):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    print(f"  {'ok  ' if ok else 'FAIL'} {label:<44} got {got:9.4f}   expected {want:9.4f}")
    if not ok:
        failures.append(label)


def expect_true(label, cond):
    print(f"  {'ok  ' if cond else 'FAIL'} {label}")
    if not cond:
        failures.append(label)


print("1. Pure major-axis bending, class 1")
r = run("", My=100e3)
expect("class", r.cls, 1)
expect("U(My) = My / (Wpl,y fy)", r.util["My"], 100e3 / MPL_Y)
expect("U(N+M)", r.util["N+M"], 100e3 / MPL_Y)

print("2. Compression 200 kN + My 150 kNm + Vz 60 kN (class 2, N effect neglected)")
r = run("", N=-200e3, Vz=60e3, My=150e3)
expect("class", r.cls, 2)
expect("U(N) = N / (A fy)", r.util["N"], 200e3 / NPL)
expect("U(Vz) = Vz / Vpl,z", r.util["Vz"], 60e3 / VPL_Z)
n_limit = min(0.25 * NPL, 0.5 * HW * TW * FY)
expect_true("200 kN is below both limits for neglecting N", 200e3 <= n_limit)
expect("U(N+M) = My / Mpl,y", r.util["N+M"], 150e3 / MPL_Y)

print("3. Tension 600 kN + My 80 kNm (class 1, N reduces Mpl,y)")
r = run("", N=600e3, My=80e3)
n = 600e3 / NPL
a = min((A - 2 * B * TF) / A, 0.5)
mn_y = MPL_Y * (1 - n) / (1 - 0.5 * a)
expect("class", r.cls, 1)
expect("U(N+M) = My / MN,y", r.util["N+M"], 80e3 / mn_y)

print("4. Compression 100 kN + My 100 kNm + Mz 10 kNm (class 1, biaxial)")
r = run("", N=-100e3, My=100e3, Mz=10e3)
n = 100e3 / NPL
ry, rz = 100e3 / MPL_Y, 10e3 / MPL_Z                  # n <= a and N small: no reduction
biaxial = ry ** 2 + rz ** max(5 * n, 1.0)
expect("class", r.cls, 1)
expect("U(N+M) = max(ry, rz, ry^2 + rz^beta)", r.util["N+M"], max(ry, rz, biaxial))

print("4b. Uniaxial overload: My 250 kNm reports the plain ratio, not its square")
r = run("", My=250e3)
expect("U(N+M) = My / Mpl,y", r.util["N+M"], 250e3 / MPL_Y)

print("5. High shear: Vz 300 kN (more than half of Vpl,z,Rd) + My 100 kNm")
r = run("", Vz=300e3, My=100e3)
rho = (2 * 300e3 / VPL_Z - 1) ** 2
mv = (WPL_Y - rho * (HW * TW) ** 2 / (4 * TW)) * FY
expect("U(Vz)", r.util["Vz"], 300e3 / VPL_Z)
expect("U(My) uses the reduced MV,y,Rd", r.util["My"], 100e3 / mv)

print("6. Pure compression 1500 kN: web c/t above 42 eps -> class 4, NOT CHECKED")
r = run("", N=-1500e3)
expect("class", r.cls, 4)
expect_true("reported as not checked", bool(r.not_checked) and not r.util)

print("7. Torsion (St Venant only)")
IT = 14.90e-8       # m4: numerical membrane-analogy solution for this section, root fillets included
TAU_RD = FY / math.sqrt(3)
expect("It of the section [cm4] (formula vs numerical)", sec.props["It"] * 1e8, IT * 1e8, tol=0.005)
r = run("", T=1.0, My=10e3)
expect_true("T = 0.001 kNm (solver noise): no note, tiny utilisation",
            not r.notes and not r.not_checked and r.util["T"] < 1e-3)
r = run("", T=2e3, My=10e3)
expect("U(T) = T t_max / It / (fy / sqrt3)", r.util["T"], 2e3 * TF / IT / TAU_RD)
expect_true("a note says only St Venant torsion is included",
            any("St Venant" in n for n in r.notes) and not r.not_checked)

print("7b. Shear resistance reduced for torsion (6.26)")
T_half = 0.5 * TAU_RD * IT / TF                 # tau_t = half of the shear stress limit
r = run("", T=T_half, Vz=100e3, My=10e3)
red = math.sqrt(1 - 0.5 / 1.25)
expect("U(T)", r.util["T"], 0.5)
expect("U(Vz) = Vz / (Vpl,z,Rd sqrt(1 - tau_t / (1.25 tau_Rd)))", r.util["Vz"],
       100e3 / (VPL_Z * red))
r = run("", T=T_half, Vz=300e3, My=100e3)       # above half of the REDUCED Vpl,z: moment reduced
vred = VPL_Z * red
rho = (2 * 300e3 / vred - 1) ** 2
mv = (WPL_Y - rho * (HW * TW) ** 2 / (4 * TW)) * FY
expect("U(My) uses rho from the reduced Vpl,T,Rd", r.util["My"], 100e3 / mv)
r = run("", T=2.5 * T_half, Vz=10e3)            # tau_t above 1.25 tau_Rd: no shear left
expect_true("torsion alone overloads the section: U(T) > 1 and shear is infinite",
            r.util["T"] > 1 and math.isinf(r.util["Vz"]))

print("8. Class 3 web in combined compression and bending: elastic stress check")
r = run("", N=-500e3, My=100e3)
sigma = max(abs(-500e3 / A + sy * 100e3 / WEL_Y) for sy in (1, -1))
expect("class", r.cls, 3)
expect("U(N+M) = max sigma / fy", r.util["N+M"], sigma / FY)
expect("U(My) uses Wel,y", r.util["My"], 100e3 / (WEL_Y * FY))

print("9. Classification numbers from first principles")
expect("flange c/tf", r.classification["ratio_f"], (B - TW - 2 * R) / 2 / TF)
expect("web c/tw", r.classification["ratio_w"], (HW - 2 * R) / TW)
alpha = 0.5 * (1 + 500e3 / ((HW - 2 * R) * TW * FY))
expect("alpha", r.classification["alpha"], alpha)
sig_mean = 500e3 / A
d = HW - 2 * R
bend = 100e3 * (d / 2) / 8.505546e-05
psi = (sig_mean - bend) / (sig_mean + bend)
expect("psi", r.classification["psi"], psi)
expect("class 3 web limit 42 eps / (0.67 + 0.33 psi)", r.classification["lim_w"][2],
       42 * EPS / (0.67 + 0.33 * psi))

print()
if failures:
    print(f"{len(failures)} FAILED: " + "; ".join(failures))
    sys.exit(1)
print("All checks passed.")
