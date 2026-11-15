"""
EC3 (EN 1993-1-1) member buckling checks, clause 6.3: flexural buckling (6.3.1) and
lateral-torsional buckling, general case (6.3.2.2). Doubly symmetric I-sections.
No Strand7 needed. Internal units are SI: N, m, Pa and N*m.

The checks are made once per member and design case, using the worst values over the member:
    N_Ed  = largest compression anywhere in the member
    My_Ed = largest |My| anywhere in the member
This is safe, and it is what a hand check normally does.

Mcr (lateral-torsional buckling)
    EN 1993-1-1 gives no formula for Mcr. Here Mcr is calculated from the real bending moment
    diagram of the member and design case: the elastic buckling problem
        E Iw phi'''' - G It phi'' - (M(x)^2 / (E Iz)) phi = 0
    is solved for a member with fork supports at both ends (twist and lateral bending
    restrained, free to warp and to rotate in plan), with the load acting at the shear centre.
    The result is written as Mcr = C1 x Mcr,uniform, where Mcr,uniform is the closed-form
    value for uniform moment over the effective length Lw (= kz x member length).

    Limits of this Mcr:
      - fork supports and loads at the shear centre; loads on the top flange (destabilising)
        would reduce Mcr, and bottom-flange loads would increase it
      - for kz different from 1 (for example a cantilever) the fork-support C1 is applied
        to the effective length Lw: treat the result as approximate
      - bending about the major axis only
"""

import math
from dataclasses import dataclass, field

from section_math import eng

PI = math.pi


# ----------------------------------------------------------------------
# Small numerical helpers (pure Python, no numpy needed)
# ----------------------------------------------------------------------

def _jacobi_max_eigenvalue(a):
    """Largest eigenvalue of a small symmetric matrix (cyclic Jacobi). a: list of lists."""
    n = len(a)
    a = [row[:] for row in a]
    for _ in range(60):
        off = sum(a[i][j] ** 2 for i in range(n) for j in range(i + 1, n))
        scale = sum(a[i][i] ** 2 for i in range(n)) or 1.0
        if off < 1e-24 * scale:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-300:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1))
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
    return max(a[i][i] for i in range(n))


def _interp(xs, ys, x):
    """Linear interpolation; xs ascending (repeated values allowed)."""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    if xs[hi] == xs[lo]:
        return ys[hi]
    return ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / (xs[hi] - xs[lo])


# ----------------------------------------------------------------------
# Elastic critical moment
# ----------------------------------------------------------------------

def mcr_uniform(E, G, Iz, It, Iw, L):
    """Elastic critical moment for uniform bending, fork supports, effective length L."""
    return (PI ** 2 * E * Iz / L ** 2) * math.sqrt(Iw / Iz + L ** 2 * G * It / (PI ** 2 * E * Iz))


def c1_from_diagram(xs, ms, length, E, G, Iz, It, Iw, modes=16, samples=800):
    """
    C1 for a moment diagram M(x), x from 0 to length (fork supports): Mcr,max / Mcr,uniform
    with both calculated over the member length. xs ascending, ms the moments (any sign).
    """
    peak = max(abs(m) for m in ms)
    if peak <= 0:
        return 1.0
    dx = length / samples
    xm = [(i + 0.5) * dx for i in range(samples)]
    m2 = [(_interp(xs, ms, x) / peak) ** 2 for x in xm]          # (M / Mmax)^2
    sines = [[math.sin(n * PI * x / length) for x in xm] for n in range(1, modes + 1)]
    # Galerkin with phi = sum a_n sin(n pi x / L): K a = (Mmax^2 / E Iz) B a
    k = [(length / 2) * (E * Iw * (n * PI / length) ** 4 + G * It * (n * PI / length) ** 2)
         for n in range(1, modes + 1)]       # the length/2 is the integral of sin^2
    b = [[sum(m2[i] * sines[p][i] * sines[q][i] for i in range(samples)) * dx / (E * Iz)
          for q in range(modes)] for p in range(modes)]
    c = [[b[p][q] / math.sqrt(k[p] * k[q]) for q in range(modes)] for p in range(modes)]
    mu = _jacobi_max_eigenvalue(c)                 # largest eigenvalue = 1 / Mmax,cr^2
    mcr_max = 1 / math.sqrt(mu)
    return mcr_max / mcr_uniform(E, G, Iz, It, Iw, length)


# ----------------------------------------------------------------------
# Buckling curves
# ----------------------------------------------------------------------

def _row_matches(row, hb, tf_mm):
    return (hb > row.get("hb_over", -math.inf) and hb <= row.get("hb_upto", math.inf)
            and tf_mm <= row.get("tf_upto", math.inf))


def flexural_curve(bk, sec, fy_mpa, axis):
    """Buckling curve name for flexural buckling about 'y' or 'z' (Table 6.2, rolled I)."""
    hb, tf_mm = sec.h / sec.b, sec.tf * 1e3
    for row in bk.flexural_rows:
        if _row_matches(row, hb, tf_mm):
            hs = fy_mpa > bk.high_strength_above
            return row.get(axis + "_hs", row[axis]) if hs else row[axis]
    raise ValueError("no flexural buckling curve row fits this section")


def ltb_curve(bk, sec):
    """Buckling curve name for lateral-torsional buckling (Table 6.4, rolled I)."""
    hb, tf_mm = sec.h / sec.b, sec.tf * 1e3
    for row in bk.ltb_rows:
        if _row_matches(row, hb, tf_mm):
            return row["curve"]
    raise ValueError("no LTB curve row fits this section")


def reduction_factor(lam, alpha, lam0, beta, cap_inverse_square):
    """chi from EN 1993-1-1 (6.49) / (6.56); returns (phi, chi). chi = 1 on the plateau."""
    if lam <= lam0:
        return None, 1.0
    phi = 0.5 * (1 + alpha * (lam - lam0) + beta * lam ** 2)
    chi = min(1.0, 1 / (phi + math.sqrt(max(phi ** 2 - beta * lam ** 2, 0.0))))
    if cap_inverse_square:
        chi = min(chi, 1 / lam ** 2)
    return phi, chi


# ----------------------------------------------------------------------
# Member check
# ----------------------------------------------------------------------

@dataclass
class BucklingResult:
    util: dict = field(default_factory=dict)      # "Nb,y", "Nb,z", "LTB" -> utilisation
    values: dict = field(default_factory=dict)    # numbers behind the checks
    not_checked: list = field(default_factory=list)
    notes: list = field(default_factory=list)
    trace: list = field(default_factory=list)

    @property
    def governing(self):
        if not self.util:
            return None, None
        name = max(self.util, key=lambda k: self.util[k])
        return name, self.util[name]


def kN(x):
    return f"{x / 1e3:.3f}"


def check_member_buckling(sec, fy, E, G, gM1, bk, stations, member_length, cls,
                          Lcr_y, Lcr_z, Lw, trace=False):
    """
    sec      : ISection
    stations : [(x, N, My)] over the member, N and N*m, tension positive (x from member start)
    cls      : section class (the highest over the member); 4 -> not checked
    Lcr_y, Lcr_z, Lw : effective lengths in m, None = that check is switched off
    """
    lines = []

    def tr(text=""):
        if trace:
            lines.append(text)

    res = BucklingResult(trace=lines)
    P = sec.props
    A, Iy, Iz, It, Iw = P["A"], P["Iy"], P["Iz"], P["It"], P["Iw"]
    fy_mpa = fy / 1e6
    n_ed = max((-n for _, n, _ in stations), default=0.0)               # compression, N
    x_n = max(stations, key=lambda s: -s[1])[0] if stations else 0.0
    m_ed = max((abs(m) for _, _, m in stations), default=0.0)
    x_m = max(stations, key=lambda s: abs(s[2]))[0] if stations else 0.0
    res.values.update({"N_Ed": n_ed, "x_N": x_n, "My_Ed": m_ed, "x_M": x_m})

    tr("MEMBER BUCKLING CHECK  (EN 1993-1-1, clause 6.3)")
    tr(f"Section : {sec.name}   fy = {fy_mpa:g} N/mm2   E = {E/1e9:.1f} kN/mm2   "
       f"G = {G/1e9:.1f} kN/mm2   gamma_M1 = {gM1:g}")
    tr(f"  A = {eng(A*1e6)} mm2   Iy = {eng(Iy*1e12)} mm4   Iz = {eng(Iz*1e12)} mm4   "
       f"It = {eng(It*1e12)} mm4   Iw = {eng(Iw*1e18)} mm6")
    tr(f"Member length = {member_length:.3f} m   section class (highest over the member) = {cls}")
    tr(f"Design values over the member:  N_Ed = {kN(n_ed)} kN compression (at x = {x_n:.3f} m)   "
       f"|My|,Ed = {kN(m_ed)} kNm (at x = {x_m:.3f} m)")
    tr()

    if cls >= 4:
        res.not_checked.append("class 4 section: member buckling with effective properties is "
                               "not implemented")
        tr("NOT CHECKED: class 4 section.")
        return res

    # --- 6.3.1 flexural buckling ---------------------------------------------------------
    tr("1. FLEXURAL BUCKLING (6.3.1)")
    npl = A * fy
    for axis, I, Lcr, key in (("y", Iy, Lcr_y, "Nb,y"), ("z", Iz, Lcr_z, "Nb,z")):
        if Lcr is None:
            tr(f"   {axis}-{axis}: switched off in members.csv")
            continue
        if n_ed <= 0:
            res.util[key] = 0.0
            tr(f"   {axis}-{axis}: no compression in the member, no check needed")
            continue
        ncr = PI ** 2 * E * I / Lcr ** 2
        lam = math.sqrt(npl / ncr)
        curve = flexural_curve(bk, sec, fy_mpa, axis)
        alpha = bk.alpha[curve]
        if n_ed / ncr <= 0.04:
            phi, chi = None, 1.0
            plateau = f"N_Ed/Ncr = {n_ed/ncr:.3f} <= 0.04: no buckling check needed (6.3.1.2(4))"
        else:
            phi, chi = reduction_factor(lam, alpha, bk.lambda_0, 1.0, False)
            plateau = None
        nb = chi * npl / gM1
        res.util[key] = n_ed / nb
        res.values[f"Lcr_{axis}"] = Lcr
        res.values[f"Ncr_{axis}"] = ncr
        res.values[f"lam_{axis}"] = lam
        res.values[f"chi_{axis}"] = chi
        res.values[f"Nb_{axis}"] = nb
        tr(f"   {axis}-{axis}: Lcr = {Lcr:.3f} m   Ncr = pi^2 E I{axis} / Lcr^2 = {kN(ncr)} kN")
        tr(f"        lambda = sqrt(A fy / Ncr) = {lam:.3f}   curve {curve} (alpha = {alpha:g})")
        if plateau:
            tr(f"        {plateau}")
        elif phi is None:
            tr(f"        lambda <= {bk.lambda_0:g}: chi = 1")
        else:
            tr(f"        phi = 0.5 (1 + alpha (lambda - {bk.lambda_0:g}) + lambda^2) = {phi:.3f}")
        tr(f"        chi = {chi:.3f}   Nb,{axis},Rd = chi A fy / gM1 = {kN(nb)} kN")
        tr(f"        N_Ed / Nb,{axis},Rd = {kN(n_ed)} / {kN(nb)} = {res.util[key]:.3f}")
    tr()

    # --- 6.3.2 lateral-torsional buckling, general case --------------------------------------
    tr("2. LATERAL-TORSIONAL BUCKLING, general case (6.3.2.2)")
    if Lw is None:
        tr("   switched off in members.csv")
    elif m_ed <= 0 or len(stations) < 2:
        res.util["LTB"] = 0.0
        tr("   no bending about y-y: no check needed")
    else:
        stations = sorted(stations)
        xs = [s[0] for s in stations]
        ms = [s[2] for s in stations]
        c1 = c1_from_diagram(xs, ms, member_length, E, G, Iz, It, Iw)
        mcr_u = mcr_uniform(E, G, Iz, It, Iw, Lw)
        mcr = c1 * mcr_u
        W = P["Wpl_y"] if cls <= 2 else P["Wel_y"]
        lam_lt = math.sqrt(W * fy / mcr)
        curve = ltb_curve(bk, sec)
        alpha_lt = bk.alpha_LT[curve]
        if m_ed / mcr <= bk.lambda_LT0 ** 2:
            phi, chi = None, 1.0
            plateau = f"M_Ed/Mcr = {m_ed/mcr:.3f} <= lambda_LT0^2: no LTB check needed (6.3.2.2(4))"
        else:
            phi, chi = reduction_factor(lam_lt, alpha_lt, bk.lambda_LT0, bk.beta_LT, True)
            plateau = None
        mb = chi * W * fy / gM1
        res.util["LTB"] = m_ed / mb
        res.values.update({"Lw": Lw, "C1": c1, "Mcr": mcr, "lam_LT": lam_lt, "chi_LT": chi,
                           "Mb": mb})
        tr(f"   Lw = {Lw:.3f} m   Mcr,uniform(Lw) = {kN(mcr_u)} kNm")
        tr(f"   C1 from the moment diagram = {c1:.3f}   Mcr = C1 Mcr,uniform = {kN(mcr)} kNm")
        tr(f"   W{'pl' if cls <= 2 else 'el'},y = {eng(W*1e9)} mm3   "
           f"lambda_LT = sqrt(W fy / Mcr) = {lam_lt:.3f}")
        tr(f"   h/b = {sec.h/sec.b:.2f}  curve {curve} (alpha_LT = {alpha_lt:g})   "
           f"lambda_LT,0 = {bk.lambda_LT0:g}   beta = {bk.beta_LT:g}")
        if plateau:
            tr(f"   {plateau}")
        elif phi is None:
            tr(f"   lambda_LT <= {bk.lambda_LT0:g}: chi_LT = 1")
        else:
            tr(f"   phi_LT = {phi:.3f}")
        tr(f"   chi_LT = {chi:.3f}   Mb,Rd = chi_LT W fy / gM1 = {kN(mb)} kNm")
        tr(f"   |My|,Ed / Mb,Rd = {kN(m_ed)} / {kN(mb)} = {res.util['LTB']:.3f}")
        res.notes.append("LTB: Mcr assumes fork supports and loads at the shear centre "
                         "(top-flange loading would give a lower Mcr)")
        if abs(Lw - member_length) > 1e-6 * member_length:
            res.notes.append(f"LTB: Lw = {Lw:.3f} m differs from the member length "
                             f"{member_length:.3f} m: the fork-support C1 is applied to Lw "
                             f"(approximate)")
    tr()

    if res.util:
        name, value = res.governing
        tr(f"GOVERNING: {name} = {value:.3f}   " + ("OK" if value <= 1.0 else "FAILS"))
    return res
