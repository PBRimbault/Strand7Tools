"""
Steel sections and steel for the EC3 checks. Pure calculation, SI units (m, Pa, N).

THE SECTION INTERFACE
The checks (cross_section.py, buckling.py) use only what the Section base class defines, so a
new profile (channel, angle, RHS, CHS ...) is added by writing one new Section subclass; the
checks themselves do not change. A section provides:

    name, shape, props (SectionProperties: A, Iy, Iz, Wel, Wpl, Av, It, Iw), t_max
    and the EC3 rules that depend on the shape:
        classify()                       Table 5.2
        torsion_shear_stress()           6.2.7 (St Venant)
        web_shear_buckling()             6.2.6(6): hw/tw against 72 eps / eta
        wpl_y_reduced_for_shear()        6.2.8(5)
        plastic_moments_with_axial()     6.2.9.1
        flexural_curve(), ltb_curve()    Tables 6.2 and 6.4

Axes: y-y is the major axis, z-z the minor axis.
"""

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass


# ----------------------------------------------------------------------
# Steel
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class Steel:
    grade: str
    fy: float           # Pa, for the section's thickness
    fu: float           # Pa
    E: float            # Pa
    G: float            # Pa


# ----------------------------------------------------------------------
# What every section provides
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class SectionProperties:
    A: float
    Iy: float
    Iz: float
    Wel_y: float
    Wel_z: float
    Wpl_y: float
    Wpl_z: float
    Av_y: float         # shear area for shear along y (the flanges of an I-section)
    Av_z: float         # shear area for shear along z (the web of an I-section)
    It: float           # St Venant torsion constant
    Iw: float           # warping constant


@dataclass(frozen=True)
class PartClass:
    """Classification of one compression part (Table 5.2)."""
    name: str                   # e.g. "Flange outstand", "Web"
    kind: str                   # "outstand" or "internal"
    symbol: str                 # the ratio, e.g. "c/tf"
    ratio: float
    limits: tuple               # class 1, 2 and 3 limits
    compressed: bool
    cls: int
    limit_labels: tuple = None  # e.g. ("9eps", "10eps", "14eps"); None = "class 1" ...
    alpha: float = None         # internal parts: compressed proportion
    psi: float = None           # internal parts: stress ratio


@dataclass(frozen=True)
class Classification:
    eps: float
    parts: tuple

    @property
    def cls(self):
        return max(p.cls for p in self.parts)

    def part(self, name):
        return next(p for p in self.parts if p.name == name)


@dataclass(frozen=True)
class AxialReduction:
    """Plastic moment resistances reduced for axial force (6.2.9.1)."""
    MN_y: float
    MN_z: float
    neglect_y: bool             # True: the effect of N on My is neglected (6.2.9.1(4))
    a: float
    alpha: float                # exponents of the biaxial criterion (6.41)
    beta: float


def class_from_limits(ratio, limits):
    return 1 if ratio <= limits[0] else 2 if ratio <= limits[1] else 3 if ratio <= limits[2] else 4


class Section(ABC):
    """Base class: the interface the EC3 checks rely on (see the module docstring)."""
    shape = "?"

    name: str
    props: SectionProperties

    @property
    @abstractmethod
    def t_max(self):
        """Thickness used to look up fy and fu (the thickest element), m."""

    @abstractmethod
    def classify(self, fy, N, My, Mz):
        """Classification (Table 5.2) under N (tension positive), My, Mz. fy in Pa."""

    @abstractmethod
    def torsion_shear_stress(self, T):
        """Largest St Venant shear stress for a torque T (N*m), Pa."""

    @abstractmethod
    def web_shear_buckling(self, eps):
        """(slenderness, limit) for shear buckling of the web (6.2.6(6)), or None."""

    @abstractmethod
    def wpl_y_reduced_for_shear(self, rho):
        """Plastic modulus about y-y reduced for high shear along z (6.2.8(5)), m3."""

    @abstractmethod
    def plastic_moments_with_axial(self, N, Npl, My_rd, Mpl_z, fy, gM0):
        """AxialReduction for class 1 and 2 sections (6.2.9.1). My_rd: the y-y resistance
        before the axial reduction (it may already be reduced for shear)."""

    @abstractmethod
    def flexural_curve(self, bk, fy_mpa, axis):
        """Buckling curve name for flexural buckling about 'y' or 'z' (Table 6.2)."""

    @abstractmethod
    def ltb_curve(self, bk):
        """Buckling curve name for lateral-torsional buckling (Table 6.4)."""


# ----------------------------------------------------------------------
# Doubly symmetric I-section
# ----------------------------------------------------------------------

def torsion_constant(h, b, tw, tf, r):
    """
    St Venant torsion constant It of a doubly symmetric I-section, root fillets included
    (the usual rolled-section formula, as in the ArcelorMittal tables). Checked against a
    numerical solution of the membrane analogy for a 305x165x40 UB (within 0.1%).
    """
    d_fillet = ((tf + r) ** 2 + (r + tw / 2) * tw) / (2 * r + tf) if r > 0 else 0.0
    alpha = (-0.042 + 0.2204 * tw / tf + 0.1355 * r / tf
             - 0.0865 * r * tw / tf ** 2 - 0.0725 * tw ** 2 / tf ** 2)
    return (2 / 3 * (b - 0.63 * tf) * tf ** 3 + 1 / 3 * (h - 2 * tf) * tw ** 3
            + 2 * alpha * d_fillet ** 4)


def i_section_properties(h, b, tw, tf, r, eta=1.0):
    """
    Properties of a doubly symmetric I-section, in the units of the inputs (use metres).
    Root fillets are included.

    eta: lower limit factor for Av,z (EN 1993-1-1 6.2.6(3): Av >= eta*hw*tw).
         Recommended value 1.2 in EN 1993-1-5; the National Annex may differ - verify.
    """
    hw = h - 2 * tf                       # web depth between flanges
    a_f = (1 - math.pi / 4) * r ** 2      # area of ONE root fillet
    c = r * (10 - 3 * math.pi) / (3 * (4 - math.pi))      # fillet centroid from web/flange faces
    i_own = r ** 4 * (1 - 5 * math.pi / 16) - a_f * c ** 2  # fillet's own I about its centroid

    area = 2 * b * tf + hw * tw + 4 * a_f
    y_f = hw / 2 - c                      # fillet centroid from the y-y axis
    x_f = tw / 2 + c                      # fillet centroid from the z-z axis

    iy = (2 * (b * tf ** 3 / 12 + b * tf * ((h - tf) / 2) ** 2)
          + tw * hw ** 3 / 12 + 4 * (i_own + a_f * y_f ** 2))
    iz = 2 * tf * b ** 3 / 12 + hw * tw ** 3 / 12 + 4 * (i_own + a_f * x_f ** 2)

    wpl_y = b * tf * (h - tf) + tw * hw ** 2 / 4 + 4 * a_f * y_f
    wpl_z = tf * b ** 2 / 2 + hw * tw ** 2 / 4 + 4 * a_f * x_f

    av_z = max(area - 2 * b * tf + (tw + 2 * r) * tf, eta * hw * tw)   # shear along the web
    av_y = area - hw * tw                                              # shear along the flanges

    return SectionProperties(
        A=area, Iy=iy, Iz=iz, Wel_y=iy / (h / 2), Wel_z=iz / (b / 2), Wpl_y=wpl_y, Wpl_z=wpl_z,
        Av_y=av_y, Av_z=av_z, It=torsion_constant(h, b, tw, tf, r),
        Iw=iz * (h - tf) ** 2 / 4,        # thin-walled approximation
    )


def _row_matches(row, hb, tf_mm):
    return (hb > row.get("hb_over", -math.inf) and hb <= row.get("hb_upto", math.inf)
            and tf_mm <= row.get("tf_upto", math.inf))


class ISection(Section):
    """
    Doubly symmetric rolled I-section (UB, UC, HE, IPE ...). Dimensions in metres.

    Classification simplifications, all on the safe side or neutral:
      - flange outstand: uniform-compression limits (9, 10, 14 eps) whenever any compression
        is present
      - web: from N and My only (Mz is ignored for the web, where its stresses are small)
    """
    shape = "I"

    def __init__(self, name, h, b, tw, tf, r, eta=1.0):
        self.name = name
        self.h, self.b, self.tw, self.tf, self.r, self.eta = h, b, tw, tf, r, eta
        self.props = i_section_properties(h, b, tw, tf, r, eta)

    def __repr__(self):
        return (f"ISection({self.name!r}, h={self.h}, b={self.b}, tw={self.tw}, tf={self.tf}, "
                f"r={self.r}, eta={self.eta})")

    # --- geometry --------------------------------------------------------------
    @property
    def hw(self):
        """Web depth between the flanges."""
        return self.h - 2 * self.tf

    @property
    def c_web(self):
        """Straight part of the web, d (for classification)."""
        return self.hw - 2 * self.r

    @property
    def c_flange(self):
        """Flange outstand c (for classification)."""
        return (self.b - self.tw - 2 * self.r) / 2

    @property
    def t_max(self):
        return max(self.tf, self.tw)

    # --- EC3 rules ---------------------------------------------------------------
    def classify(self, fy, N, My, Mz):
        P = self.props
        eps = math.sqrt(235e6 / fy)
        sigma_mean = -N / P.A                                     # compression positive

        # flange outstand (uniform-compression limits whenever there is compression)
        ratio_f = self.c_flange / self.tf
        lim_f = (9 * eps, 10 * eps, 14 * eps)
        flange_compressed = sigma_mean + abs(My) / P.Wel_y + abs(Mz) / P.Wel_z > 0
        cls_f = class_from_limits(ratio_f, lim_f) if flange_compressed else 1
        flange = PartClass("Flange outstand", "outstand", "c/tf", ratio_f, lim_f,
                           flange_compressed, cls_f, limit_labels=("9eps", "10eps", "14eps"))

        # web (internal part, straight portion d), from N and My
        d = self.c_web
        ratio_w = d / self.tw
        sig_max = sigma_mean + abs(My) * (d / 2) / P.Iy
        sig_min = sigma_mean - abs(My) * (d / 2) / P.Iy
        web_compressed = sig_max > 0
        alpha = min(max(0.5 * (1 - N / (d * self.tw * fy)), 0.0), 1.0)   # compressed part
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
        cls_w = class_from_limits(ratio_w, lim_w) if web_compressed else 1
        web = PartClass("Web", "internal", "c/tw", ratio_w, lim_w, web_compressed, cls_w,
                        alpha=alpha, psi=psi)
        return Classification(eps, (flange, web))

    def torsion_shear_stress(self, T):
        return abs(T) * self.t_max / self.props.It

    def web_shear_buckling(self, eps):
        return self.hw / self.tw, 72 * eps / self.eta

    def wpl_y_reduced_for_shear(self, rho):
        aw = self.hw * self.tw
        return self.props.Wpl_y - rho * aw ** 2 / (4 * self.tw)

    def plastic_moments_with_axial(self, N, Npl, My_rd, Mpl_z, fy, gM0):
        P = self.props
        n = abs(N) / Npl
        a = min((P.A - 2 * self.b * self.tf) / P.A, 0.5)
        neglect_y = abs(N) <= 0.25 * Npl and abs(N) <= 0.5 * self.hw * self.tw * fy / gM0
        MN_y = My_rd if neglect_y else min(My_rd, My_rd * (1 - n) / (1 - 0.5 * a))
        MN_z = Mpl_z if n <= a else Mpl_z * (1 - ((n - a) / (1 - a)) ** 2)
        return AxialReduction(MN_y, MN_z, neglect_y, a, 2.0, max(5 * n, 1.0))

    def flexural_curve(self, bk, fy_mpa, axis):
        """Table 6.2, rolled I-sections (the rows are in settings.toml)."""
        hb, tf_mm = self.h / self.b, self.tf * 1e3
        for row in bk.flexural_rows:
            if _row_matches(row, hb, tf_mm):
                hs = fy_mpa > bk.high_strength_above
                return row.get(axis + "_hs", row[axis]) if hs else row[axis]
        raise ValueError("no flexural buckling curve row fits this section")

    def ltb_curve(self, bk):
        """Table 6.4, rolled I-sections (the rows are in settings.toml)."""
        hb, tf_mm = self.h / self.b, self.tf * 1e3
        for row in bk.ltb_rows:
            if _row_matches(row, hb, tf_mm):
                return row["curve"]
        raise ValueError("no LTB curve row fits this section")
