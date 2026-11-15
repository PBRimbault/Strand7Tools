"""
EC3 (EN 1993-1-1) cross-section checks, clause 6.2, for doubly symmetric I-sections.

Pure calculation: no Strand7 needed. Internal units are SI: N, m, Pa and N*m.
(Printed output uses kN, kNm, mm and N/mm2.)

Forces are in the MEMBER axes, tension positive:
    N   axial force            Vy, Vz  shear along y (flanges) and z (web)
    T   torque                 My, Mz  moments about the major (y-y) and minor (z-z) axes

What is checked, per station:
    - section classification (Table 5.2)                      class 4 -> NOT CHECKED
    - axial resistance (6.2.3, 6.2.4)                         gross section, no bolt holes
    - bending resistance about each axis (6.2.5)
    - shear resistance, plastic (6.2.6)                       web shear buckling -> NOT CHECKED
    - bending and shear (6.2.8)                               major axis, class 1 and 2
    - bending and axial force, uniaxial and biaxial (6.2.9)   class 1/2 plastic, class 3 elastic
    - torsion, St Venant only (6.2.7): tau_t = T t_max / It against fy / (sqrt3 gM0), and the
      shear resistance reduced to Vpl,T,Rd (6.26)

Anything this module does not cover is reported as NOT CHECKED with the reason. It is never
silently passed:
    - class 4 sections
    - web shear buckling (hw/tw > 72 eps / eta)
    - shear above 50% of the plastic shear resistance with a class 3 section, or in the y direction
    - warping torsion and its combination with bending (only St Venant torsion is checked)

Classification simplifications, all on the safe side or neutral:
    - flange outstand: uniform-compression limits (9, 10, 14 eps) whenever any compression is present
    - web: from N and My only (Mz is ignored for the web, where its stresses are small)
"""

import math
from dataclasses import dataclass, field

from section_math import eng, i_section_properties

SQRT3 = math.sqrt(3.0)
TORSION_TOL = 1e-3          # torque above this fraction of Mpl,z,Rd gets a note (solver noise below)


# ----------------------------------------------------------------------
# Section
# ----------------------------------------------------------------------

@dataclass
class ISection:
    name: str
    h: float
    b: float
    tw: float
    tf: float
    r: float
    eta: float
    props: dict = field(default_factory=dict)

    @classmethod
    def from_dims(cls, name, h, b, tw, tf, r, eta=1.0):
        """Dimensions in metres."""
        return cls(name, h, b, tw, tf, r, eta, i_section_properties(h, b, tw, tf, r, eta))

    @property
    def t_max(self):
        return max(self.tf, self.tw)


# ----------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------

class _Trace:
    """Collects the lines of a worked calculation (does nothing when switched off)."""

    def __init__(self, on):
        self.on = on
        self.lines = []

    def __call__(self, text=""):
        if self.on:
            self.lines.append(text)


def kN(x):
    return f"{x / 1e3:.3f}"


def _ratio(load, capacity):
    """load / capacity, with a non-positive capacity meaning infinite utilisation under load."""
    if capacity > 0:
        return load / capacity
    return math.inf if load > 1e-9 else 0.0


# ----------------------------------------------------------------------
# Classification (EN 1993-1-1 Table 5.2)
# ----------------------------------------------------------------------

def classify(sec, fy, N, My, Mz):
    """Section class 1-4 and the numbers behind it. fy in Pa, N tension positive."""
    P = sec.props
    eps = math.sqrt(235e6 / fy)
    A, Iy = P["A"], P["Iy"]

    # --- flange outstand ---------------------------------------------------
    c_f = P["c_flange"]
    ratio_f = c_f / sec.tf
    lim_f = (9 * eps, 10 * eps, 14 * eps)
    sigma_mean = -N / A                                       # compression positive
    max_flange_compression = sigma_mean + abs(My) / P["Wel_y"] + abs(Mz) / P["Wel_z"]
    flange_compressed = max_flange_compression > 0

    def pick(ratio, lims):
        return 1 if ratio <= lims[0] else 2 if ratio <= lims[1] else 3 if ratio <= lims[2] else 4

    cls_f = pick(ratio_f, lim_f) if flange_compressed else 1

    # --- web (internal part, straight portion d) ------------------------------
    d = P["c_web"]
    ratio_w = d / sec.tw
    sig_max = sigma_mean + abs(My) * (d / 2) / Iy
    sig_min = sigma_mean - abs(My) * (d / 2) / Iy
    web_compressed = sig_max > 0
    alpha = min(max(0.5 * (1 - N / (d * sec.tw * fy)), 0.0), 1.0)    # compressed part of the web
    psi = sig_min / sig_max if web_compressed else None

    if alpha > 0.5:
        l1, l2 = 396 * eps / (13 * alpha - 1), 456 * eps / (13 * alpha - 1)
    elif alpha > 0:
        l1, l2 = 36 * eps / alpha, 41.5 * eps / alpha
    else:
        l1 = l2 = math.inf
    if psi is None:
        l3 = math.inf
    elif psi > -1:
        l3 = 42 * eps / (0.67 + 0.33 * psi)
    else:
        l3 = 62 * eps * (1 - psi) * math.sqrt(-psi)
    lim_w = (l1, l2, l3)
    cls_w = pick(ratio_w, lim_w) if web_compressed else 1

    return {"eps": eps, "ratio_f": ratio_f, "lim_f": lim_f, "cls_f": cls_f,
            "flange_compressed": flange_compressed,
            "ratio_w": ratio_w, "lim_w": lim_w, "cls_w": cls_w,
            "web_compressed": web_compressed, "alpha": alpha, "psi": psi,
            "cls": max(cls_f, cls_w)}


# ----------------------------------------------------------------------
# Station check
# ----------------------------------------------------------------------

@dataclass
class StationResult:
    cls: int
    classification: dict
    util: dict                  # check name -> utilisation (empty if not checked)
    not_checked: list           # reasons the section could not be fully checked
    notes: list                 # information that does not stop the check
    trace: list                 # worked calculation lines (only if requested)

    @property
    def governing(self):
        if not self.util:
            return None, None
        name = max(self.util, key=lambda k: self.util[k])
        return name, self.util[name]


def check_station(sec, fy, gM0, forces, trace=False):
    """
    Cross-section checks at one station.

    sec    : ISection
    fy     : yield strength, Pa
    gM0    : partial factor
    forces : dict with N, Vy, Vz, T, My, Mz (N and N*m, member axes, tension positive)
    trace  : True to also return the worked calculation as text lines
    """
    tr = _Trace(trace)
    P = sec.props
    A, hw, tw, tf, b = P["A"], P["hw"], sec.tw, sec.tf, sec.b
    N, Vy, Vz, Tq = forces["N"], forces["Vy"], forces["Vz"], forces["T"]
    My, Mz = forces["My"], forces["Mz"]
    notes, not_checked = [], []

    tr("CROSS-SECTION CHECK  (EN 1993-1-1, clause 6.2)")
    tr(f"Section : {sec.name}")
    tr(f"  h = {sec.h*1e3:.1f}  b = {b*1e3:.1f}  tw = {tw*1e3:.1f}  tf = {tf*1e3:.1f}  "
       f"r = {sec.r*1e3:.1f}  hw = {hw*1e3:.1f}  d = {P['c_web']*1e3:.1f}  (mm)")
    tr(f"  A = {eng(A*1e6)} mm2   Wpl,y = {eng(P['Wpl_y']*1e9)} mm3   "
       f"Wpl,z = {eng(P['Wpl_z']*1e9)} mm3")
    tr(f"  Wel,y = {eng(P['Wel_y']*1e9)} mm3   Wel,z = {eng(P['Wel_z']*1e9)} mm3   "
       f"Av,y = {eng(P['Av_y']*1e6)} mm2   Av,z = {eng(P['Av_z']*1e6)} mm2")
    tr(f"Steel   : fy = {fy/1e6:g} N/mm2 (for t = {sec.t_max*1e3:.1f} mm)   gamma_M0 = {gM0:g}")
    tr()
    tr("Design forces at this station (member axes, tension positive):")
    tr(f"  N = {kN(N)} kN    Vy = {kN(Vy)} kN    Vz = {kN(Vz)} kN")
    tr(f"  T = {kN(Tq)} kNm   My = {kN(My)} kNm   Mz = {kN(Mz)} kNm")
    tr()

    # --- classification --------------------------------------------------------
    cl = classify(sec, fy, N, My, Mz)
    eps = cl["eps"]
    cls = cl["cls"]
    tr(f"1. CLASSIFICATION (Table 5.2)   eps = sqrt(235/fy) = {eps:.3f}")
    if cl["flange_compressed"]:
        f1, f2, f3 = cl["lim_f"]
        tr(f"   Flange outstand: c/tf = {cl['ratio_f']:.2f}   limits 9eps = {f1:.2f}, "
           f"10eps = {f2:.2f}, 14eps = {f3:.2f}   -> class {cl['cls_f']}")
    else:
        tr("   Flange outstand: no compression                                       -> class 1")
    if cl["web_compressed"]:
        w1, w2, w3 = cl["lim_w"]
        psi = cl["psi"]
        tr(f"   Web: c/tw = {cl['ratio_w']:.2f}   alpha = {cl['alpha']:.3f}   psi = {psi:.3f}")
        tr(f"        limits class 1 = {w1:.2f}, class 2 = {w2:.2f}, class 3 = {w3:.2f}   "
           f"-> class {cl['cls_w']}")
    else:
        tr("   Web: no compression                                                   -> class 1")
    tr(f"   Section class = {cls}")
    tr()

    if cls == 4:
        not_checked.append("class 4 section: effective section properties (EN 1993-1-5) "
                           "are not implemented")
        tr("   NOT CHECKED: class 4 sections need effective properties (not implemented).")
        return StationResult(cls, cl, {}, not_checked, notes, tr.lines)

    # --- resistances ---------------------------------------------------------
    Npl = A * fy / gM0
    Mpl_y, Mel_y = P["Wpl_y"] * fy / gM0, P["Wel_y"] * fy / gM0
    Mpl_z, Mel_z = P["Wpl_z"] * fy / gM0, P["Wel_z"] * fy / gM0
    Mc_y = Mpl_y if cls <= 2 else Mel_y
    Mc_z = Mpl_z if cls <= 2 else Mel_z
    Vpl_y = P["Av_y"] * fy / (SQRT3 * gM0)
    Vpl_z = P["Av_z"] * fy / (SQRT3 * gM0)
    # --- torsion, St Venant only (6.2.7) -------------------------------------------
    tau_Rd = fy / (SQRT3 * gM0)
    tau_t = abs(Tq) * sec.t_max / P["It"]
    reduction = math.sqrt(max(1 - tau_t / (1.25 * tau_Rd), 0.0))      # (6.26)
    Vpl_y *= reduction
    Vpl_z *= reduction
    which = "Wpl" if cls <= 2 else "Wel"
    tr("2. RESISTANCES")
    tr(f"   Npl,Rd  = A fy / gM0 = {kN(Npl)} kN                         (6.2.3, 6.2.4)")
    tr(f"   Mc,y,Rd = {which},y fy / gM0 = {kN(Mc_y)} kNm                      (6.2.5)")
    tr(f"   Mc,z,Rd = {which},z fy / gM0 = {kN(Mc_z)} kNm")
    if tau_t > 0:
        tr(f"   Vpl,y,Rd = Av,y fy / (sqrt3 gM0) x {reduction:.3f} = {kN(Vpl_y)} kN       "
           f"(6.2.6, reduced for torsion)")
        tr(f"   Vpl,z,Rd = Av,z fy / (sqrt3 gM0) x {reduction:.3f} = {kN(Vpl_z)} kN")
    else:
        tr(f"   Vpl,y,Rd = Av,y fy / (sqrt3 gM0) = {kN(Vpl_y)} kN                 (6.2.6)")
        tr(f"   Vpl,z,Rd = Av,z fy / (sqrt3 gM0) = {kN(Vpl_z)} kN")
    tr()

    # --- things not covered ----------------------------------------------------
    if hw / tw > 72 * eps / sec.eta:
        not_checked.append(f"web shear buckling: hw/tw = {hw/tw:.1f} exceeds 72 eps/eta = "
                           f"{72*eps/sec.eta:.1f}; EN 1993-1-5 is not implemented")
    util = {}
    if abs(Tq) > TORSION_TOL * Mpl_z:
        notes.append("torsion is checked as St Venant torsion only (6.2.7); warping stresses "
                     "and their combination with bending are not included")

    # --- shear -------------------------------------------------------------------
    util["Vy"] = _ratio(abs(Vy), Vpl_y)
    util["Vz"] = _ratio(abs(Vz), Vpl_z)
    util["T"] = tau_t / tau_Rd
    tr("3. TORSION AND SHEAR (6.2.6, 6.2.7)")
    tr(f"   It = {eng(P['It']*1e12)} mm4   t_max = {sec.t_max*1e3:.1f} mm")
    tr(f"   tau_t = |T| t_max / It = {tau_t/1e6:.2f} N/mm2")
    tr(f"   tau_t / (fy / (sqrt3 gM0)) = {tau_t/1e6:.2f} / {tau_Rd/1e6:.1f} = {util['T']:.3f}")
    tr("   (all of T is taken by St Venant torsion; warping is not included)")
    tr(f"   |Vy| / Vpl,y,Rd = {kN(abs(Vy))} / {kN(Vpl_y)} = {util['Vy']:.3f}")
    tr(f"   |Vz| / Vpl,z,Rd = {kN(abs(Vz))} / {kN(Vpl_z)} = {util['Vz']:.3f}")
    tr()

    # --- bending and shear (6.2.8) ---------------------------------------------------
    My_res = Mc_y
    tr("4. BENDING AND SHEAR (6.2.8)")
    if abs(Vz) > 0.5 * Vpl_z and Vpl_z > 0:
        if cls == 3:
            not_checked.append("class 3 section with Vz above 50% of Vpl,z,Rd: the reduced "
                               "moment resistance is not implemented")
            tr("   Vz > 0.5 Vpl,z,Rd on a class 3 section: NOT CHECKED (not implemented).")
        else:
            rho = min((2 * abs(Vz) / Vpl_z - 1) ** 2, 1.0)
            Aw = hw * tw
            Mv_y = (P["Wpl_y"] - rho * Aw ** 2 / (4 * tw)) * fy / gM0
            My_res = min(Mc_y, Mv_y)
            tr(f"   |Vz| > 0.5 Vpl,z,Rd: rho = (2 |Vz| / Vpl,z,Rd - 1)^2 = {rho:.3f}")
            tr(f"   MV,y,Rd = (Wpl,y - rho Aw^2 / (4 tw)) fy / gM0 = {kN(Mv_y)} kNm")
            notes.append("Vz above 50% of Vpl,z,Rd: moment resistance reduced (6.2.8); the "
                         "interaction with N is then approximate")
    else:
        tr("   |Vz| <= 0.5 Vpl,z,Rd: no reduction of the moment resistance about y-y.")
    if abs(Vy) > 0.5 * Vpl_y:
        not_checked.append("Vy above 50% of Vpl,y,Rd: the reduction of Mz resistance is "
                           "not implemented")
        tr("   |Vy| > 0.5 Vpl,y,Rd: NOT CHECKED (reduction of Mz resistance not implemented).")
    tr()

    # --- axial force alone and bending alone ---------------------------------------------
    n = abs(N) / Npl
    util["N"] = n
    util["My"] = _ratio(abs(My), My_res)
    util["Mz"] = _ratio(abs(Mz), Mc_z)
    tr("5. AXIAL FORCE AND BENDING ON THEIR OWN")
    tr(f"   n = |N| / Npl,Rd = {kN(abs(N))} / {kN(Npl)} = {n:.3f}")
    tr(f"   |My| / Mc,y,Rd = {kN(abs(My))} / {kN(My_res)} = {util['My']:.3f}")
    tr(f"   |Mz| / Mc,z,Rd = {kN(abs(Mz))} / {kN(Mc_z)} = {util['Mz']:.3f}")
    tr()

    # --- bending and axial force (6.2.9) ------------------------------------------------------
    tr("6. BENDING AND AXIAL FORCE (6.2.9)")
    if cls <= 2:
        a = min((A - 2 * b * tf) / A, 0.5)
        if abs(N) <= 0.25 * Npl and abs(N) <= 0.5 * hw * tw * fy / gM0:
            MN_y = My_res
            tr("   |N| <= 0.25 Npl,Rd and |N| <= 0.5 hw tw fy/gM0: effect of N on My neglected "
               "(6.2.9.1(4))")
        else:
            MN_y = min(My_res, My_res * (1 - n) / (1 - 0.5 * a))
            tr(f"   a = (A - 2 b tf) / A = {a:.3f} (max 0.5)")
            tr(f"   MN,y,Rd = Mpl,y,Rd (1 - n) / (1 - 0.5 a) = {kN(MN_y)} kNm")
        if n <= a:
            MN_z = Mpl_z
        else:
            MN_z = Mpl_z * (1 - ((n - a) / (1 - a)) ** 2)
        tr(f"   MN,z,Rd = {kN(MN_z)} kNm" + ("  (n <= a: no reduction)" if n <= a else ""))
        ry, rz = _ratio(abs(My), MN_y), _ratio(abs(Mz), MN_z)
        tr(f"   |My| / MN,y,Rd = {ry:.3f}     |Mz| / MN,z,Rd = {rz:.3f}")
        if ry > 1e-3 and rz > 1e-3:
            beta = max(5 * n, 1.0)
            biaxial = ry ** 2 + rz ** beta
            util["N+M"] = max(ry, rz, biaxial)
            tr(f"   biaxial (6.2.9.1(6)): (My/MN,y)^2 + (Mz/MN,z)^beta with beta = "
               f"max(5n, 1) = {beta:.3f}  ->  {biaxial:.3f}")
            tr(f"   utilisation = max of the three = {util['N+M']:.3f}")
        else:
            util["N+M"] = max(ry, rz)
            tr(f"   bending about one axis only: utilisation = {util['N+M']:.3f}")
    else:
        corners = [N / A + sy * abs(My) / P["Wel_y"] + sz * abs(Mz) / P["Wel_z"]
                   for sy in (1, -1) for sz in (1, -1)]
        sig = max(abs(c) for c in corners)
        util["N+M"] = sig / (fy / gM0)
        tr("   class 3: elastic check, sigma = N/A + My/Wel,y + Mz/Wel,z at the corners (6.2.9.2)")
        tr(f"   max |sigma| = {sig/1e6:.1f} N/mm2   fy/gM0 = {fy/gM0/1e6:.1f} N/mm2   "
           f"-> {util['N+M']:.3f}")
    tr()

    name, value = max(util.items(), key=lambda kv: kv[1])
    tr(f"GOVERNING: {name} = {value:.3f}   " + ("OK" if value <= 1.0 else "FAILS")
       + ("   (NOT FULLY CHECKED: " + "; ".join(not_checked) + ")" if not_checked else ""))
    return StationResult(cls, cl, util, not_checked, notes, tr.lines)


# ----------------------------------------------------------------------
# Worked example (no Strand7 needed)
# ----------------------------------------------------------------------

def example():
    """Two worked calculations for the 305x165x40 UB (BS EN 10365 dimensions) in S355."""
    sec = ISection.from_dims("UB 305x165x40", 0.3034, 0.165, 0.006, 0.0102, 0.009, eta=1.0)
    cases = [
        # tension + bending + shear: class 1, N reduces the plastic moment resistance
        {"N": 600e3, "Vy": 0.0, "Vz": 60e3, "T": 0.0, "My": 80e3, "Mz": 0.0},
        # compression + bending: the web becomes class 3, elastic check
        {"N": -500e3, "Vy": 0.0, "Vz": 40e3, "T": 0.0, "My": 100e3, "Mz": 0.0},
    ]
    for i, forces in enumerate(cases):
        if i:
            print("\n" + "=" * 76 + "\n")
        print("\n".join(check_station(sec, 355e6, 1.0, forces, trace=True).trace))


if __name__ == "__main__":
    example()
