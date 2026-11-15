"""
Section property maths for doubly symmetric I-sections. Pure calculation: no Strand7 needed.

Used by section.py (which reads dimensions from the model) and by ec3_checks.py.
"""

import math
import re


def eng(x, decimals=3):
    """Engineering notation: mantissa 1..1000, exponent a multiple of 3 (e.g. 85.056e+06)."""
    if x == 0:
        return f"{0:.{decimals}f}e+00"
    exp = int(math.floor(math.log10(abs(x)) / 3) * 3)
    return f"{x / 10 ** exp:.{decimals}f}e{exp:+03d}"


def short_name(name):
    """'BS EN - UB 305 x 165 x 40 - BS EN 10365-2017 BGL' -> 'UB 305x165x40' (else unchanged)."""
    m = re.search(r"([A-Za-z]+)\s*(\d+(?:\.\d+)?(?:\s*[x\u00d7]\s*\d+(?:\.\d+)?)+)", name)
    if not m:
        return name
    dims = re.sub(r"\s*[x\u00d7]\s*", "x", m.group(2))
    return f"{m.group(1)} {dims}"


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
    Properties of a doubly symmetric I-section, in the units of the inputs
    (use metres). Root fillets are included.

    eta: lower limit factor for Av,z (EN 1993-1-1 6.2.6(3): Av >= eta*hw*tw).
         Recommended value 1.2 in EN 1993-1-5; the National Annex may differ -
         verify.
    """
    hw = h - 2 * tf                       # web depth between flanges
    d = hw - 2 * r                        # straight part of the web (for classification)
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

    return {
        "A": area, "Iy": iy, "Iz": iz,
        "Wel_y": iy / (h / 2), "Wel_z": iz / (b / 2),
        "Wpl_y": wpl_y, "Wpl_z": wpl_z,
        "Av_z": av_z, "Av_y": av_y,
        "It": torsion_constant(h, b, tw, tf, r),
        "Iw": iz * (h - tf) ** 2 / 4,     # thin-walled approximation
        "hw": hw, "c_flange": (b - tw - 2 * r) / 2, "c_web": d,
    }
