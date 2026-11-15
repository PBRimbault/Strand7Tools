"""
Worked calculations: turn a check result into step-by-step text that can be checked by hand.
The calculations themselves are in ec3/; nothing here changes a number.

  cross_section_lines(CrossSectionResult)  -> [str]    EC3 6.2 at one station
  buckling_lines(BucklingResult)           -> [str]    EC3 6.3 for a member and design case
"""

from .fmt import eng, kN


# ----------------------------------------------------------------------
# Cross-section (6.2)
# ----------------------------------------------------------------------

def _section_dimension_lines(sec):
    if sec.shape == "I":
        return [f"  h = {sec.h*1e3:.1f}  b = {sec.b*1e3:.1f}  tw = {sec.tw*1e3:.1f}  "
                f"tf = {sec.tf*1e3:.1f}  r = {sec.r*1e3:.1f}  hw = {sec.hw*1e3:.1f}  "
                f"d = {sec.c_web*1e3:.1f}  (mm)"]
    return [f"  {sec!r}"]


def _classification_lines(cl):
    out = [f"1. CLASSIFICATION (Table 5.2)   eps = sqrt(235/fy) = {cl.eps:.3f}"]
    for p in cl.parts:
        if not p.compressed:
            out.append(f"   {p.name}: no compression".ljust(73) + "-> class 1")
            continue
        l1, l2, l3 = p.limits
        if p.kind == "outstand":
            a, b, c = p.limit_labels or ("class 1", "class 2", "class 3")
            out.append(f"   {p.name}: {p.symbol} = {p.ratio:.2f}   limits {a} = {l1:.2f}, "
                       f"{b} = {l2:.2f}, {c} = {l3:.2f}   -> class {p.cls}")
        else:
            out.append(f"   {p.name}: {p.symbol} = {p.ratio:.2f}   alpha = {p.alpha:.3f}   "
                       f"psi = {p.psi:.3f}")
            out.append(f"        limits class 1 = {l1:.2f}, class 2 = {l2:.2f}, "
                       f"class 3 = {l3:.2f}   -> class {p.cls}")
    out.append(f"   Section class = {cl.cls}")
    return out


def cross_section_lines(r):
    sec, P, v, u = r.section, r.section.props, r.values, r.util
    fy, gM0 = r.steel.fy, r.gM0
    f = r.forces
    out = ["CROSS-SECTION CHECK  (EN 1993-1-1, clause 6.2)", f"Section : {sec.name}"]
    out += _section_dimension_lines(sec)
    out += [
        f"  A = {eng(P.A*1e6)} mm2   Wpl,y = {eng(P.Wpl_y*1e9)} mm3   "
        f"Wpl,z = {eng(P.Wpl_z*1e9)} mm3",
        f"  Wel,y = {eng(P.Wel_y*1e9)} mm3   Wel,z = {eng(P.Wel_z*1e9)} mm3   "
        f"Av,y = {eng(P.Av_y*1e6)} mm2   Av,z = {eng(P.Av_z*1e6)} mm2",
        f"Steel   : fy = {fy/1e6:g} N/mm2 (for t = {sec.t_max*1e3:.1f} mm)   gamma_M0 = {gM0:g}",
        "",
        "Design forces at this station (member axes, tension positive):",
        f"  N = {kN(f['N'])} kN    Vy = {kN(f['Vy'])} kN    Vz = {kN(f['Vz'])} kN",
        f"  T = {kN(f['T'])} kNm   My = {kN(f['My'])} kNm   Mz = {kN(f['Mz'])} kNm",
        "",
    ]
    out += _classification_lines(r.classification)
    out.append("")
    if r.cls == 4:
        out.append("   NOT CHECKED: class 4 sections need effective properties (not implemented).")
        return out

    W = v["W"]
    out.append("2. RESISTANCES")
    out.append(f"   Npl,Rd  = A fy / gM0 = {kN(v['Npl'])} kN                         (6.2.3, 6.2.4)")
    out.append(f"   Mc,y,Rd = {W},y fy / gM0 = {kN(v['Mc_y'])} kNm                      (6.2.5)")
    out.append(f"   Mc,z,Rd = {W},z fy / gM0 = {kN(v['Mc_z'])} kNm")
    if v["tau_t"] > 0:
        red = v["shear_reduction"]
        out.append(f"   Vpl,y,Rd = Av,y fy / (sqrt3 gM0) x {red:.3f} = {kN(v['Vpl_y'])} kN       "
                   f"(6.2.6, reduced for torsion)")
        out.append(f"   Vpl,z,Rd = Av,z fy / (sqrt3 gM0) x {red:.3f} = {kN(v['Vpl_z'])} kN")
    else:
        out.append(f"   Vpl,y,Rd = Av,y fy / (sqrt3 gM0) = {kN(v['Vpl_y'])} kN                 (6.2.6)")
        out.append(f"   Vpl,z,Rd = Av,z fy / (sqrt3 gM0) = {kN(v['Vpl_z'])} kN")
    out.append("")

    out += [
        "3. TORSION AND SHEAR (6.2.6, 6.2.7)",
        f"   It = {eng(P.It*1e12)} mm4   t_max = {sec.t_max*1e3:.1f} mm",
        f"   tau_t = |T| t_max / It = {v['tau_t']/1e6:.2f} N/mm2",
        f"   tau_t / (fy / (sqrt3 gM0)) = {v['tau_t']/1e6:.2f} / {v['tau_Rd']/1e6:.1f} = {u['T']:.3f}",
        "   (all of T is taken by St Venant torsion; warping is not included)",
        f"   |Vy| / Vpl,y,Rd = {kN(abs(f['Vy']))} / {kN(v['Vpl_y'])} = {u['Vy']:.3f}",
        f"   |Vz| / Vpl,z,Rd = {kN(abs(f['Vz']))} / {kN(v['Vpl_z'])} = {u['Vz']:.3f}",
        "",
        "4. BENDING AND SHEAR (6.2.8)",
    ]
    if v["vz_high"]:
        if "rho" not in v:
            out.append("   Vz > 0.5 Vpl,z,Rd on a class 3 section: NOT CHECKED (not implemented).")
        else:
            out.append(f"   |Vz| > 0.5 Vpl,z,Rd: rho = (2 |Vz| / Vpl,z,Rd - 1)^2 = {v['rho']:.3f}")
            out.append(f"   MV,y,Rd = (Wpl,y - rho Aw^2 / (4 tw)) fy / gM0 = {kN(v['Mv_y'])} kNm")
    else:
        out.append("   |Vz| <= 0.5 Vpl,z,Rd: no reduction of the moment resistance about y-y.")
    if v["vy_high"]:
        out.append("   |Vy| > 0.5 Vpl,y,Rd: NOT CHECKED (reduction of Mz resistance not implemented).")
    out.append("")

    out += [
        "5. AXIAL FORCE AND BENDING ON THEIR OWN",
        f"   n = |N| / Npl,Rd = {kN(abs(f['N']))} / {kN(v['Npl'])} = {v['n']:.3f}",
        f"   |My| / Mc,y,Rd = {kN(abs(f['My']))} / {kN(v['My_res'])} = {u['My']:.3f}",
        f"   |Mz| / Mc,z,Rd = {kN(abs(f['Mz']))} / {kN(v['Mc_z'])} = {u['Mz']:.3f}",
        "",
        "6. BENDING AND AXIAL FORCE (6.2.9)",
    ]
    ax = r.axial
    if ax is not None:
        if ax.neglect_y:
            out.append("   |N| <= 0.25 Npl,Rd and |N| <= 0.5 hw tw fy/gM0: effect of N on My "
                       "neglected (6.2.9.1(4))")
        else:
            out.append(f"   a = (A - 2 b tf) / A = {ax.a:.3f} (max 0.5)")
            out.append(f"   MN,y,Rd = Mpl,y,Rd (1 - n) / (1 - 0.5 a) = {kN(ax.MN_y)} kNm")
        out.append(f"   MN,z,Rd = {kN(ax.MN_z)} kNm"
                   + ("  (n <= a: no reduction)" if v["n"] <= ax.a else ""))
        out.append(f"   |My| / MN,y,Rd = {v['ry']:.3f}     |Mz| / MN,z,Rd = {v['rz']:.3f}")
        if "biaxial" in v:
            out.append(f"   biaxial (6.2.9.1(6)): (My/MN,y)^2 + (Mz/MN,z)^beta with beta = "
                       f"max(5n, 1) = {ax.beta:.3f}  ->  {v['biaxial']:.3f}")
            out.append(f"   utilisation = max of the three = {u['N+M']:.3f}")
        else:
            out.append(f"   bending about one axis only: utilisation = {u['N+M']:.3f}")
    else:
        out.append("   class 3: elastic check, sigma = N/A + My/Wel,y + Mz/Wel,z at the corners "
                   "(6.2.9.2)")
        out.append(f"   max |sigma| = {v['sigma_max']/1e6:.1f} N/mm2   fy/gM0 = {fy/gM0/1e6:.1f} "
                   f"N/mm2   -> {u['N+M']:.3f}")
    out.append("")

    name, value = r.governing
    out.append(f"GOVERNING: {name} = {value:.3f}   " + ("OK" if value <= 1.0 else "FAILS")
               + ("   (NOT FULLY CHECKED: " + "; ".join(r.not_checked) + ")"
                  if r.not_checked else ""))
    return out


# ----------------------------------------------------------------------
# Member buckling (6.3)
# ----------------------------------------------------------------------

def buckling_lines(res):
    sec, P, v, bk = res.section, res.section.props, res.values, res.bk
    fy, E, G = res.steel.fy, res.steel.E, res.steel.G
    n_ed, m_ed = v["N_Ed"], v["My_Ed"]
    out = [
        "MEMBER BUCKLING CHECK  (EN 1993-1-1, clause 6.3)",
        f"Section : {sec.name}   fy = {fy/1e6:g} N/mm2   E = {E/1e9:.1f} kN/mm2   "
        f"G = {G/1e9:.1f} kN/mm2   gamma_M1 = {res.gM1:g}",
        f"  A = {eng(P.A*1e6)} mm2   Iy = {eng(P.Iy*1e12)} mm4   Iz = {eng(P.Iz*1e12)} mm4   "
        f"It = {eng(P.It*1e12)} mm4   Iw = {eng(P.Iw*1e18)} mm6",
        f"Member length = {res.member_length:.3f} m   section class (highest over the member) "
        f"= {res.cls}",
        f"Design values over the member:  N_Ed = {kN(n_ed)} kN compression (at x = "
        f"{v['x_N']:.3f} m)   |My|,Ed = {kN(m_ed)} kNm (at x = {v['x_M']:.3f} m)",
        "",
    ]
    if res.cls >= 4:
        out.append("NOT CHECKED: class 4 section.")
        return out

    out.append("1. FLEXURAL BUCKLING (6.3.1)")
    for axis in ("y", "z"):
        state = res.flexural_state[axis]
        if state == "off":
            out.append(f"   {axis}-{axis}: switched off in members.csv")
            continue
        if state == "none":
            out.append(f"   {axis}-{axis}: no compression in the member, no check needed")
            continue
        ncr, nb, phi = v[f"Ncr_{axis}"], v[f"Nb_{axis}"], v[f"phi_{axis}"]
        out.append(f"   {axis}-{axis}: Lcr = {v[f'Lcr_{axis}']:.3f} m   "
                   f"Ncr = pi^2 E I{axis} / Lcr^2 = {kN(ncr)} kN")
        out.append(f"        lambda = sqrt(A fy / Ncr) = {v[f'lam_{axis}']:.3f}   "
                   f"curve {v[f'curve_{axis}']} (alpha = {v[f'alpha_{axis}']:g})")
        if v[f"plateau_{axis}"]:
            out.append(f"        N_Ed/Ncr = {n_ed/ncr:.3f} <= 0.04: no buckling check needed "
                       f"(6.3.1.2(4))")
        elif phi is None:
            out.append(f"        lambda <= {bk.lambda_0:g}: chi = 1")
        else:
            out.append(f"        phi = 0.5 (1 + alpha (lambda - {bk.lambda_0:g}) + lambda^2) "
                       f"= {phi:.3f}")
        out.append(f"        chi = {v[f'chi_{axis}']:.3f}   Nb,{axis},Rd = chi A fy / gM1 = "
                   f"{kN(nb)} kN")
        out.append(f"        N_Ed / Nb,{axis},Rd = {kN(n_ed)} / {kN(nb)} = "
                   f"{res.util['Nb,' + axis]:.3f}")
    out.append("")

    out.append("2. LATERAL-TORSIONAL BUCKLING, general case (6.3.2.2)")
    if res.ltb_state == "off":
        out.append("   switched off in members.csv")
    elif res.ltb_state == "none":
        out.append("   no bending about y-y: no check needed")
    else:
        mcr, mb, phi = v["Mcr"], v["Mb"], v["phi_LT"]
        out += [
            f"   Lw = {v['Lw']:.3f} m   Mcr,uniform(Lw) = {kN(v['Mcr_uniform'])} kNm",
            f"   C1 from the moment diagram = {v['C1']:.3f}   Mcr = C1 Mcr,uniform = "
            f"{kN(mcr)} kNm",
            f"   W{'pl' if res.cls <= 2 else 'el'},y = {eng(v['W_y']*1e9)} mm3   "
            f"lambda_LT = sqrt(W fy / Mcr) = {v['lam_LT']:.3f}",
        ]
        if sec.shape == "I":
            out.append(f"   h/b = {sec.h/sec.b:.2f}  curve {v['curve_LT']} "
                       f"(alpha_LT = {v['alpha_LT']:g})   lambda_LT,0 = {bk.lambda_LT0:g}   "
                       f"beta = {bk.beta_LT:g}")
        else:
            out.append(f"   curve {v['curve_LT']} (alpha_LT = {v['alpha_LT']:g})   "
                       f"lambda_LT,0 = {bk.lambda_LT0:g}   beta = {bk.beta_LT:g}")
        if v["plateau_LT"]:
            out.append(f"   M_Ed/Mcr = {m_ed/mcr:.3f} <= lambda_LT0^2: no LTB check needed "
                       f"(6.3.2.2(4))")
        elif phi is None:
            out.append(f"   lambda_LT <= {bk.lambda_LT0:g}: chi_LT = 1")
        else:
            out.append(f"   phi_LT = {phi:.3f}")
        out.append(f"   chi_LT = {v['chi_LT']:.3f}   Mb,Rd = chi_LT W fy / gM1 = {kN(mb)} kNm")
        out.append(f"   |My|,Ed / Mb,Rd = {kN(m_ed)} / {kN(mb)} = {res.util['LTB']:.3f}")
    out.append("")

    out.append("3. BENDING AND AXIAL COMPRESSION (6.3.3)")
    if res.interaction_state == "flagged":
        out.append(f"   N_Ed = {kN(n_ed)} kN with |My|,Ed = {kN(m_ed)} kNm and |Mz|,Ed = "
                   f"{kN(v['Mz_Ed'])} kNm:")
        out.append("   NOT CHECKED - the interaction (6.61)/(6.62) is not implemented yet.")
    elif res.interaction_state == "off":
        out.append("   all buckling checks switched off in members.csv")
    else:
        out.append("   not needed: the member is not under both compression and bending")
    out.append("")

    if res.util:
        name, value = res.governing
        verdict = "FAILS" if value > 1.0 else "INCOMPLETE" if res.not_checked else "OK"
        out.append(f"GOVERNING: {name} = {value:.3f}   {verdict}"
                   + ("   (NOT FULLY CHECKED: " + "; ".join(res.not_checked) + ")"
                      if res.not_checked else ""))
    return out
