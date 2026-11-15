"""
Section dimensions and material from a Strand7 beam property. Read-only.

read_section() returns the raw data (dimensions in m, E and G in Pa, Strand7's own section
values); make_section() turns it into an EC3 Section object. Supported so far: doubly symmetric
I-sections, from the Beam Geometry Library (with root radius) or the standard I-section
(no root radius: Wpl and Av come out slightly low, with a warning).

Axes: y-y is the MAJOR axis (Strand7 principal axis 1, bending in plane 2), z-z the minor axis.
"""

import ctypes
import re

from ..ec3.sections import ISection
from .api import S, check, get_str, model_units


class SectionError(Exception):
    """The section cannot be handled; the message says why."""


LENGTH_TO_M = {"luMETRE": 1.0, "luCENTIMETRE": 0.01, "luMILLIMETRE": 0.001,
               "luFOOT": 0.3048, "luINCH": 0.0254}
STRESS_TO_PA = {"suPASCAL": 1.0, "suKILOPASCAL": 1e3, "suMEGAPASCAL": 1e6}


def length_factor(uid):
    units = model_units(uid)
    for name, factor in LENGTH_TO_M.items():
        if getattr(S, name) == units[S.ipLENGTHU]:
            return factor
    raise SectionError("unknown model length unit")


def stress_factor(uid):
    """Factor from the model's stress unit (Strand7 keeps a separate one) to Pa."""
    units = model_units(uid)
    for name, factor in STRESS_TO_PA.items():
        if getattr(S, name) == units[S.ipSTRESSU]:
            return factor
    raise SectionError("model stress unit not supported (only Pa, kPa and MPa)")


def short_name(name):
    """'BS EN - UB 305 x 165 x 40 - BS EN 10365-2017 BGL' -> 'UB 305x165x40' (else unchanged)."""
    m = re.search(r"([A-Za-z]+)\s*(\d+(?:\.\d+)?(?:\s*[x×]\s*\d+(?:\.\d+)?)+)", name)
    if not m:
        return name
    dims = re.sub(r"\s*[x×]\s*", "x", m.group(2))
    return f"{m.group(1)} {dims}"


def read_section(uid, prop):
    """Dimensions (in metres), material (Pa) and Strand7's own values for a beam property."""
    integers = (ctypes.c_long * 5)()
    sec = (ctypes.c_double * S.kNumBeamSectionData)()
    mat = (ctypes.c_double * S.kNumMaterialData)()
    check(S.St7GetBeamPropertyData(uid, prop, integers, sec, mat))
    sec_type = integers[S.ipBeamPropSectionType]
    name = get_str(S.St7GetBeamSectionName, uid, prop)
    f = length_factor(uid)

    if sec_type in (S.bsBGLSection, S.bsBXSSection):
        shape = ctypes.c_long()
        dims = (ctypes.c_double * S.kMaxBGLDimensions)()
        try:
            check(S.St7GetBeamSectionGeometryBGL(uid, prop, ctypes.byref(shape), dims))
        except RuntimeError:
            raise SectionError("not a Beam Geometry Library section")
        if shape.value != S.bgISection:
            raise SectionError("library section is not an I-section")
        d_, b1, b2, tw, tf1, tf2, r1, r2 = dims[0:8]      # D, B1, B2, tw, tf1, tf2, Rr1, Rr2
        kind, r_known = "library I-section (with root radius)", True
    elif sec_type == S.bsISection:
        b1, b2, d_ = sec[S.ipD1], sec[S.ipD2], sec[S.ipD3]
        tf1, tf2, tw = sec[S.ipT1], sec[S.ipT2], sec[S.ipT3]
        r1 = r2 = 0.0
        kind, r_known = "standard I-section (no root radius)", False
    else:
        raise SectionError("section type not supported (only I-sections so far)")

    if (abs(b1 - b2) > 1e-6 * d_ or abs(tf1 - tf2) > 1e-6 * d_ or abs(r1 - r2) > 1e-6 * d_):
        raise SectionError("monosymmetric I-section not supported yet")

    material = (ctypes.c_double * 9)()
    check(S.St7GetBeamMaterialData(uid, prop, material))
    # Strand7 gives the modulus and shear modulus in the model's STRESS unit: convert to Pa.
    stress_to_pa = stress_factor(uid)

    return {
        "shape": "I",
        "E": material[S.ipBeamModulus] * stress_to_pa,
        "G": material[S.ipBeamShear] * stress_to_pa,
        "name": name, "kind": kind, "r_known": r_known,
        "h": d_ * f, "b": b1 * f, "tw": tw * f, "tf": tf1 * f, "r": r1 * f,
        "strand7": {"A": sec[S.ipAREA] * f ** 2, "Iy": sec[S.ipI11] * f ** 4,
                    "Iz": sec[S.ipI22] * f ** 4,
                    "SAy": sec[S.ipSA1] * f ** 2, "SAz": sec[S.ipSA2] * f ** 2},
    }


def make_section(data, eta):
    """(Section, warnings) from read_section() data."""
    if data["shape"] != "I":
        raise SectionError(f"no EC3 section class for shape {data['shape']!r}")
    sec = ISection(short_name(data["name"]), data["h"], data["b"], data["tw"], data["tf"],
                   data["r"], eta)
    warnings = []
    if not data["r_known"]:
        warnings.append("no root radius in the section data: r = 0 was used, so Wpl and Av "
                        "are slightly low")
    return sec, warnings
