"""
Step 5a: EC3 section properties for I-sections.  Read-only.

For every beam property in the model that is a doubly symmetric I-section
(from the Beam Geometry Library, or a standard I-section), this computes the
properties the EC3 checks need, from the section DIMENSIONS, including the
root fillets:

    A, Iy, Iz, Wel,y, Wel,z, Wpl,y, Wpl,z, Av,z, Av,y, Iw, and the
    geometry used for classification (c/t of the flange outstand and web).

It prints them next to Strand7's own A, I11 and I22 as a built-in check, and
the values can be compared with the published section tables.

Standard I-section properties (no root radius in the data) are computed with
r = 0, with a warning: Wpl and Av come out slightly low.

Axes: y-y is the MAJOR axis (Strand7 principal axis 1, bending in plane 2),
      z-z is the minor axis.
"""

import ctypes
import math

from st7_setup import St7API as S, check, model_session
from members import get_str
from section_math import eng, i_section_properties  # noqa: F401


class SectionError(Exception):
    """The section cannot be handled; the message says why."""


LENGTH_TO_M = {"luMETRE": 1.0, "luCENTIMETRE": 0.01, "luMILLIMETRE": 0.001,
               "luFOOT": 0.3048, "luINCH": 0.0254}


def length_factor(uid):
    units = (ctypes.c_long * S.kLastUnit)()
    check(S.St7GetUnits(uid, units))
    for name, factor in LENGTH_TO_M.items():
        if getattr(S, name) == units[S.ipLENGTHU]:
            return factor
    raise SectionError("unknown model length unit")


def read_i_section(uid, prop):
    """Dimensions (in metres) and Strand7's own values for a beam property."""
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

    return {
        "name": name, "kind": kind, "r_known": r_known,
        "h": d_ * f, "b": b1 * f, "tw": tw * f, "tf": tf1 * f, "r": r1 * f,
        "strand7": {"A": sec[S.ipAREA] * f ** 2, "Iy": sec[S.ipI11] * f ** 4,
                    "Iz": sec[S.ipI22] * f ** 4,
                    "SAy": sec[S.ipSA1] * f ** 2, "SAz": sec[S.ipSA2] * f ** 2},
    }


def main():
    # (label, unit, key, factor from SI to the unit, Strand7 key or None)
    # All printed in mm (from SI metres), engineering notation, y before z.
    rows = [("A", "mm2", "A", 1e6, "A"), ("Iy", "mm4", "Iy", 1e12, "Iy"),
            ("Iz", "mm4", "Iz", 1e12, "Iz"), ("Wel,y", "mm3", "Wel_y", 1e9, None),
            ("Wel,z", "mm3", "Wel_z", 1e9, None), ("Wpl,y", "mm3", "Wpl_y", 1e9, None),
            ("Wpl,z", "mm3", "Wpl_z", 1e9, None), ("Av,y", "mm2", "Av_y", 1e6, None),
            ("Av,z", "mm2", "Av_z", 1e6, None), ("Iw", "mm6", "Iw", 1e18, None)]

    with model_session() as uid:
        n_props = (ctypes.c_long * S.kMaxEntityTotals)()
        last_prop = (ctypes.c_long * S.kMaxEntityTotals)()
        check(S.St7GetTotalProperties(uid, n_props, last_prop))

        for prop in range(1, last_prop[S.ipBeamPropTotal] + 1):
            try:
                s = read_i_section(uid, prop)
            except SectionError as e:
                print(f"Property {prop}: skipped - {e}\n")
                continue
            except RuntimeError:
                continue  # no property with this number

            p = i_section_properties(s["h"], s["b"], s["tw"], s["tf"], s["r"])
            print("=" * 70)
            print(f"Property {prop}: {s['name']}")
            print(f"  {s['kind']}")
            print(f"  h = {s['h']*1e3:.1f}  b = {s['b']*1e3:.1f}  tw = {s['tw']*1e3:.1f}  "
                  f"tf = {s['tf']*1e3:.1f}  r = {s['r']*1e3:.1f}   (mm)")
            if not s["r_known"]:
                print("  WARNING: no root radius available, so r = 0 was used; "
                      "Wpl and Av are slightly low.")
            print(f"  {'':8} {'unit':>5} {'this app':>13} {'Strand7':>13} {'diff':>8}")
            for label, unit, key, fac, s7key in rows:
                val = p[key] * fac
                if s7key:
                    ref = s["strand7"][s7key] * fac
                    print(f"  {label:<8} {unit:>5} {eng(val):>13} {eng(ref):>13} "
                          f"{(val / ref - 1) * 100:>7.2f}%")
                else:
                    print(f"  {label:<8} {unit:>5} {eng(val):>13}")
            print("  Av,y and Av,z are the EC3 6.2.6(3) areas for plastic shear RESISTANCE.")
            print(f"  Strand7 elastic shear areas (stiffness only, not used in the checks): "
                  f"SA,y = {eng(s['strand7']['SAy'] * 1e6)} mm2   "
                  f"SA,z = {eng(s['strand7']['SAz'] * 1e6)} mm2")
            print(f"  classification geometry: flange outstand c/tf = "
                  f"{p['c_flange'] / s['tf']:.2f}   web c/tw = {p['c_web'] / s['tw']:.2f}")
            print()


if __name__ == "__main__":
    main()
