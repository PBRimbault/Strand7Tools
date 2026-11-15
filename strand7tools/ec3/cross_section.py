"""
EC3 (EN 1993-1-1) cross-section checks, clause 6.2, at one station. Pure calculation: no
Strand7 and no printing (worked.py turns a result into a worked calculation).
Internal units are SI: N, m, Pa and N*m.

Forces are in the MEMBER axes, tension positive:
    N   axial force            Vy, Vz  shear along y (flanges) and z (web)
    T   torque                 My, Mz  moments about the major (y-y) and minor (z-z) axes

What is checked:
    - section classification (Table 5.2)                      class 4 -> NOT CHECKED
    - axial resistance (6.2.3, 6.2.4)                         gross section, no bolt holes
    - bending resistance about each axis (6.2.5)
    - shear resistance, plastic (6.2.6)                       web shear buckling -> NOT CHECKED
    - bending and shear (6.2.8)                               major axis, class 1 and 2
    - bending and axial force, uniaxial and biaxial (6.2.9)   class 1/2 plastic, class 3 elastic
    - torsion, St Venant only (6.2.7): tau_t against fy / (sqrt3 gM0), and the shear
      resistance reduced to Vpl,T,Rd (6.26)

Anything not covered is reported as NOT CHECKED with the reason; it is never silently passed:
    - class 4 sections
    - web shear buckling (hw/tw > 72 eps / eta)
    - shear above 50% of the plastic shear resistance with a class 3 section, or in the y direction
    - warping torsion and its combination with bending (only St Venant torsion is checked)

The shape-dependent rules (classification, 6.2.8, 6.2.9.1 ...) come from the Section object
(sections.py), so this module works for any section that implements that interface.
"""

import math
from dataclasses import dataclass, field

SQRT3 = math.sqrt(3.0)
TORSION_TOL = 1e-3          # torque above this fraction of Mpl,z,Rd gets a note (solver noise below)


def _ratio(load, capacity):
    """load / capacity, with a non-positive capacity meaning infinite utilisation under load."""
    if capacity > 0:
        return load / capacity
    return math.inf if load > 1e-9 else 0.0


@dataclass
class CrossSectionResult:
    section: object             # Section
    steel: object               # Steel
    gM0: float
    forces: dict                # N, Vy, Vz, T, My, Mz in N and N*m
    classification: object      # Classification
    util: dict = field(default_factory=dict)         # check name -> utilisation
    values: dict = field(default_factory=dict)       # resistances and intermediate values (SI)
    axial: object = None        # AxialReduction (class 1 and 2)
    not_checked: list = field(default_factory=list)  # reasons the section was not fully checked
    notes: list = field(default_factory=list)        # information that does not stop the check

    @property
    def cls(self):
        return self.classification.cls

    @property
    def governing(self):
        if not self.util:
            return None, None
        name = max(self.util, key=lambda k: self.util[k])
        return name, self.util[name]


def check_cross_section(sec, steel, gM0, forces):
    """
    Cross-section checks at one station.

    sec    : Section
    steel  : Steel (only fy is used here)
    gM0    : partial factor
    forces : dict with N, Vy, Vz, T, My, Mz (N and N*m, member axes, tension positive)
    """
    P = sec.props
    fy = steel.fy
    N, Vy, Vz, Tq = forces["N"], forces["Vy"], forces["Vz"], forces["T"]
    My, Mz = forces["My"], forces["Mz"]

    cl = sec.classify(fy, N, My, Mz)
    res = CrossSectionResult(sec, steel, gM0, dict(forces), cl)
    v, util = res.values, res.util
    cls = cl.cls
    eps = cl.eps

    if cls == 4:
        res.not_checked.append("class 4 section: effective section properties (EN 1993-1-5) "
                               "are not implemented")
        return res

    # --- resistances ---------------------------------------------------------------
    Npl = P.A * fy / gM0
    Mpl_y, Mel_y = P.Wpl_y * fy / gM0, P.Wel_y * fy / gM0
    Mpl_z, Mel_z = P.Wpl_z * fy / gM0, P.Wel_z * fy / gM0
    Mc_y = Mpl_y if cls <= 2 else Mel_y
    Mc_z = Mpl_z if cls <= 2 else Mel_z
    # torsion, St Venant only (6.2.7), and the shear resistance reduced for it (6.26)
    tau_Rd = fy / (SQRT3 * gM0)
    tau_t = sec.torsion_shear_stress(Tq)
    reduction = math.sqrt(max(1 - tau_t / (1.25 * tau_Rd), 0.0))
    Vpl_y = P.Av_y * fy / (SQRT3 * gM0) * reduction
    Vpl_z = P.Av_z * fy / (SQRT3 * gM0) * reduction
    v.update({"Npl": Npl, "Mpl_y": Mpl_y, "Mel_y": Mel_y, "Mpl_z": Mpl_z, "Mel_z": Mel_z,
              "Mc_y": Mc_y, "Mc_z": Mc_z, "W": "Wpl" if cls <= 2 else "Wel",
              "tau_Rd": tau_Rd, "tau_t": tau_t, "shear_reduction": reduction,
              "Vpl_y": Vpl_y, "Vpl_z": Vpl_z})

    # --- things not covered --------------------------------------------------------
    sb = sec.web_shear_buckling(eps)
    if sb is not None and sb[0] > sb[1]:
        res.not_checked.append(f"web shear buckling: hw/tw = {sb[0]:.1f} exceeds 72 eps/eta = "
                               f"{sb[1]:.1f}; EN 1993-1-5 is not implemented")
    if abs(Tq) > TORSION_TOL * Mpl_z:
        res.notes.append("torsion is checked as St Venant torsion only (6.2.7); warping stresses "
                         "and their combination with bending are not included")

    # --- shear and torsion -----------------------------------------------------------
    util["Vy"] = _ratio(abs(Vy), Vpl_y)
    util["Vz"] = _ratio(abs(Vz), Vpl_z)
    util["T"] = tau_t / tau_Rd

    # --- bending and shear (6.2.8) -----------------------------------------------------
    My_res = Mc_y
    v["vz_high"] = abs(Vz) > 0.5 * Vpl_z and Vpl_z > 0
    if v["vz_high"]:
        if cls == 3:
            res.not_checked.append("class 3 section with Vz above 50% of Vpl,z,Rd: the reduced "
                                   "moment resistance is not implemented")
        else:
            rho = min((2 * abs(Vz) / Vpl_z - 1) ** 2, 1.0)
            Mv_y = sec.wpl_y_reduced_for_shear(rho) * fy / gM0
            My_res = min(Mc_y, Mv_y)
            v.update({"rho": rho, "Mv_y": Mv_y})
            res.notes.append("Vz above 50% of Vpl,z,Rd: moment resistance reduced (6.2.8); the "
                             "interaction with N is then approximate")
    v["vy_high"] = abs(Vy) > 0.5 * Vpl_y
    if v["vy_high"]:
        res.not_checked.append("Vy above 50% of Vpl,y,Rd: the reduction of Mz resistance is "
                               "not implemented")
    v["My_res"] = My_res

    # --- axial force alone and bending alone ----------------------------------------------
    n = abs(N) / Npl
    v["n"] = n
    util["N"] = n
    util["My"] = _ratio(abs(My), My_res)
    util["Mz"] = _ratio(abs(Mz), Mc_z)

    # --- bending and axial force (6.2.9) ---------------------------------------------------
    if cls <= 2:
        ax = sec.plastic_moments_with_axial(N, Npl, My_res, Mpl_z, fy, gM0)
        res.axial = ax
        ry, rz = _ratio(abs(My), ax.MN_y), _ratio(abs(Mz), ax.MN_z)
        v.update({"ry": ry, "rz": rz})
        if ry > 1e-3 and rz > 1e-3:
            biaxial = ry ** ax.alpha + rz ** ax.beta
            v["biaxial"] = biaxial
            util["N+M"] = max(ry, rz, biaxial)
        else:
            util["N+M"] = max(ry, rz)
    else:
        corners = [N / P.A + sy * abs(My) / P.Wel_y + sz * abs(Mz) / P.Wel_z
                   for sy in (1, -1) for sz in (1, -1)]
        sig = max(abs(c) for c in corners)
        v["sigma_max"] = sig
        util["N+M"] = sig / (fy / gM0)
    return res
