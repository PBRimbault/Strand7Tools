"""
EC3 (EN 1993-1-1) member buckling checks, clause 6.3: flexural buckling (6.3.1) and
lateral-torsional buckling, general case (6.3.2.2). Pure calculation: no Strand7 and no
printing (worked.py turns a result into a worked calculation). SI units: N, m, Pa and N*m.

The checks are made once per member and design case, using the worst values over the member:
    N_Ed  = largest compression anywhere in the member
    My_Ed = largest |My| anywhere in the member
This is safe, and it is what a hand check normally does.

Not implemented yet: the bending and axial compression interaction (6.3.3). A member with both
compression and bending (above INTERACTION_TOL of Npl / Mpl, i.e. not solver noise) is reported
as NOT FULLY CHECKED, so it can never show as OK.

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

PI = math.pi
INTERACTION_TOL = 1e-3      # N_Ed / Npl and M_Ed / Mpl below this count as zero (solver noise)


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
# Reduction factor
# ----------------------------------------------------------------------

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
    section: object             # Section
    steel: object               # Steel
    gM1: float
    bk: object                  # settings.Buckling
    member_length: float
    cls: int
    util: dict = field(default_factory=dict)       # "Nb,y", "Nb,z", "LTB" -> utilisation
    values: dict = field(default_factory=dict)     # numbers behind the checks (SI)
    # what happened to each check: "off" (switched off in the member table), "none" (no
    # compression / no bending, nothing to check), "checked"; None = not reached (class 4)
    flexural_state: dict = field(default_factory=dict)      # "y", "z" -> state
    ltb_state: str = None
    interaction_state: str = None   # "flagged", "off" (all buckling off), "not needed"
    not_checked: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    @property
    def governing(self):
        if not self.util:
            return None, None
        name = max(self.util, key=lambda k: self.util[k])
        return name, self.util[name]

    @property
    def status(self):
        """FAILS wins even if part of the check was not possible; OK only if all was checked."""
        _, value = self.governing
        if value is not None and value > 1.0:
            return "FAILS (incomplete)" if self.not_checked else "FAILS"
        if value is None:
            return "NOT CHECKED"
        return "INCOMPLETE" if self.not_checked else "OK"


def check_member_buckling(sec, steel, gM1, bk, stations, member_length, cls,
                          Lcr_y, Lcr_z, Lw, mz_ed=0.0):
    """
    sec      : Section
    steel    : Steel (fy, E, G)
    bk       : settings.Buckling (curves, imperfection factors, plateau lengths)
    stations : [(x, N, My)] over the member, N and N*m, tension positive (x from member start)
    cls      : section class (the highest over the member); 4 -> not checked
    Lcr_y, Lcr_z, Lw : effective lengths in m, None = that check is switched off
    mz_ed    : largest |Mz| anywhere in the member, N*m (only used to decide whether the
               bending and axial compression interaction, 6.3.3, is needed)
    """
    res = BucklingResult(sec, steel, gM1, bk, member_length, cls)
    v = res.values
    P = sec.props
    fy, E, G = steel.fy, steel.E, steel.G
    fy_mpa = fy / 1e6
    n_ed = max(0.0, max((-n for _, n, _ in stations), default=0.0))     # compression, N (>= 0)
    x_n = max(stations, key=lambda s: -s[1])[0] if stations else 0.0
    m_ed = max((abs(m) for _, _, m in stations), default=0.0)
    x_m = max(stations, key=lambda s: abs(s[2]))[0] if stations else 0.0
    v.update({"N_Ed": n_ed, "x_N": x_n, "My_Ed": m_ed, "x_M": x_m, "Mz_Ed": mz_ed})

    if cls >= 4:
        res.not_checked.append("class 4 section: member buckling with effective properties is "
                               "not implemented")
        return res

    # --- 6.3.1 flexural buckling -----------------------------------------------------------
    npl = P.A * fy
    for axis, I, Lcr, key in (("y", P.Iy, Lcr_y, "Nb,y"), ("z", P.Iz, Lcr_z, "Nb,z")):
        if Lcr is None:
            res.flexural_state[axis] = "off"
            continue
        if n_ed <= 0:
            res.flexural_state[axis] = "none"
            res.util[key] = 0.0
            continue
        res.flexural_state[axis] = "checked"
        ncr = PI ** 2 * E * I / Lcr ** 2
        lam = math.sqrt(npl / ncr)
        curve = sec.flexural_curve(bk, fy_mpa, axis)
        alpha = bk.alpha[curve]
        if n_ed / ncr <= 0.04:                       # 6.3.1.2(4): no check needed
            phi, chi, plateau = None, 1.0, True
        else:
            phi, chi = reduction_factor(lam, alpha, bk.lambda_0, 1.0, False)
            plateau = False
        nb = chi * npl / gM1
        res.util[key] = n_ed / nb
        v.update({f"Lcr_{axis}": Lcr, f"Ncr_{axis}": ncr, f"lam_{axis}": lam,
                  f"curve_{axis}": curve, f"alpha_{axis}": alpha, f"phi_{axis}": phi,
                  f"plateau_{axis}": plateau, f"chi_{axis}": chi, f"Nb_{axis}": nb})

    # --- 6.3.2 lateral-torsional buckling, general case ------------------------------------
    if Lw is None:
        res.ltb_state = "off"
    elif m_ed <= 0 or len(stations) < 2:
        res.ltb_state = "none"
        res.util["LTB"] = 0.0
    else:
        res.ltb_state = "checked"
        stations = sorted(stations)
        xs = [s[0] for s in stations]
        ms = [s[2] for s in stations]
        c1 = c1_from_diagram(xs, ms, member_length, E, G, P.Iz, P.It, P.Iw)
        mcr_u = mcr_uniform(E, G, P.Iz, P.It, P.Iw, Lw)
        mcr = c1 * mcr_u
        W = P.Wpl_y if cls <= 2 else P.Wel_y
        lam_lt = math.sqrt(W * fy / mcr)
        curve = sec.ltb_curve(bk)
        alpha_lt = bk.alpha_LT[curve]
        if m_ed / mcr <= bk.lambda_LT0 ** 2:         # 6.3.2.2(4): no check needed
            phi, chi, plateau = None, 1.0, True
        else:
            phi, chi = reduction_factor(lam_lt, alpha_lt, bk.lambda_LT0, bk.beta_LT, True)
            plateau = False
        mb = chi * W * fy / gM1
        res.util["LTB"] = m_ed / mb
        v.update({"Lw": Lw, "Mcr_uniform": mcr_u, "C1": c1, "Mcr": mcr, "W_y": W,
                  "lam_LT": lam_lt, "curve_LT": curve, "alpha_LT": alpha_lt, "phi_LT": phi,
                  "plateau_LT": plateau, "chi_LT": chi, "Mb": mb})
        res.notes.append("LTB: Mcr assumes fork supports and loads at the shear centre "
                         "(top-flange loading would give a lower Mcr)")
        if abs(Lw - member_length) > 1e-6 * member_length:
            res.notes.append(f"LTB: Lw = {Lw:.3f} m differs from the member length "
                             f"{member_length:.3f} m: the fork-support C1 is applied to Lw "
                             f"(approximate)")

    # --- 6.3.3 bending and axial compression: not implemented yet ---------------------------
    # The checks above treat N and M separately. A member with both compression and bending
    # needs the interaction (6.61)/(6.62), so it must not be reported as fully checked.
    compression = n_ed > INTERACTION_TOL * npl
    bending = (m_ed > INTERACTION_TOL * P.Wpl_y * fy
               or mz_ed > INTERACTION_TOL * P.Wpl_z * fy)
    any_buckling = Lcr_y is not None or Lcr_z is not None or Lw is not None
    if compression and bending and any_buckling:
        res.interaction_state = "flagged"
        res.not_checked.append(
            f"bending and axial compression interaction (6.3.3) is not implemented: "
            f"N_Ed = {n_ed/1e3:.1f} kN with |My| = {m_ed/1e3:.1f} kNm, |Mz| = {mz_ed/1e3:.1f} kNm")
    elif not any_buckling:
        res.interaction_state = "off"
    else:
        res.interaction_state = "not needed"
    return res
