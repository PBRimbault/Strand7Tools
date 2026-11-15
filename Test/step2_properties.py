"""
Step 2a: List beam properties and beams from the model. Read-only.

Prints:
  - the model units
  - every beam property: name, section name, section type, dimensions,
    section properties and material data
  - every beam: end nodes, property, length, group

Uses only functions documented in the Strand7 R3 API manual.
"""

import ctypes

from st7_setup import St7API as S, check, model_session

BEAM_TYPE_NAMES = ["btNull", "btSpring", "btCable", "btTruss", "btCutoff",
                   "btContact", "btBeam", "btUser", "btPipe", "btConnection"]
SECTION_NAMES = ["bsCircularHollow", "bsCircularSolid", "bsCruciform",
                 "bsISection", "bsLipChannel", "bsLSection", "bsSquareHollow",
                 "bsSquareSolid", "bsTopHatChannel", "bsTrapezoidHollow",
                 "bsTrapezoidSolid", "bsTriangleHollow", "bsTriangleSolid",
                 "bsTSection", "bsZSection", "bsBXSSection", "bsBGLSection",
                 "bsNullSection"]
LENGTH_UNITS = ["luMETRE", "luCENTIMETRE", "luMILLIMETRE", "luFOOT", "luINCH"]
FORCE_UNITS = ["fuNEWTON", "fuKILONEWTON", "fuMEGANEWTON", "fuKILOFORCE",
               "fuPOUNDFORCE", "fuTONNEFORCE", "fuKIPFORCE"]
STRESS_UNITS = ["suPASCAL", "suKILOPASCAL", "suMEGAPASCAL", "suKSCm",
                "suPSI", "suKSI", "suPSF"]


def name_of(value, candidates):
    """Look up which of the named constants has this value."""
    for n in candidates:
        if hasattr(S, n) and getattr(S, n) == value:
            return n
    return f"unknown({value})"


def get_str(fn, *args):
    """Call an API function whose last two args are (buffer, max length)."""
    buf = ctypes.create_string_buffer(S.kMaxStrLen)
    check(fn(*args, buf, S.kMaxStrLen))
    return buf.value.decode("cp1252", errors="replace")  # Strand7 strings are Windows ANSI


def print_units(uid):
    units = (ctypes.c_long * S.kLastUnit)()
    check(S.St7GetUnits(uid, units))
    print("Model units:",
          name_of(units[S.ipLENGTHU], LENGTH_UNITS),
          name_of(units[S.ipFORCEU], FORCE_UNITS),
          name_of(units[S.ipSTRESSU], STRESS_UNITS))
    print()


def list_beam_properties(uid):
    n_props = (ctypes.c_long * S.kMaxEntityTotals)()
    last_prop = (ctypes.c_long * S.kMaxEntityTotals)()
    check(S.St7GetTotalProperties(uid, n_props, last_prop))
    print(f"Beam properties: {n_props[S.ipBeamPropTotal]} "
          f"(highest number {last_prop[S.ipBeamPropTotal]})")
    print("-" * 70)

    found = []
    for prop in range(1, last_prop[S.ipBeamPropTotal] + 1):
        integers = (ctypes.c_long * 5)()
        sec = (ctypes.c_double * S.kNumBeamSectionData)()
        mat = (ctypes.c_double * S.kNumMaterialData)()
        try:
            check(S.St7GetBeamPropertyData(uid, prop, integers, sec, mat))
        except RuntimeError:
            continue  # no property with this number
        beam_type = integers[S.ipBeamPropBeamType]
        if beam_type == S.btNull:
            continue
        found.append(prop)

        sec_type = integers[S.ipBeamPropSectionType]
        print(f"Property {prop}: {get_str(S.St7GetPropertyName, uid, S.ptBEAMPROP, prop)!r}")
        print(f"  beam type     : {name_of(beam_type, BEAM_TYPE_NAMES)}")
        print(f"  section name  : {get_str(S.St7GetBeamSectionName, uid, prop)!r}")
        print(f"  section type  : {name_of(sec_type, SECTION_NAMES)}")
        print(f"  material name : {get_str(S.St7GetMaterialName, uid, S.ptBEAMPROP, prop)!r}")

        # Dimensions. Standard sections use D1..D3, T1..T3; BGL sections have
        # their own dimension array.
        if sec_type in (S.bsBGLSection, S.bsBXSSection):
            # A library (BGL) section can be reported as a BXS section, so try the BGL call.
            shape = ctypes.c_long()
            dims = (ctypes.c_double * S.kMaxBGLDimensions)()
            try:
                check(S.St7GetBeamSectionGeometryBGL(uid, prop, ctypes.byref(shape), dims))
            except RuntimeError as e:
                print(f"  BGL dimensions: not available ({e})")
            else:
                shape_name = name_of(shape.value, ["bgNullSection", "bgRectangularHollow",
                                                   "bgISection", "bgChannel", "bgTSection",
                                                   "bgAngle", "bgBulbFlat"])
                print(f"  BGL shape     : {shape_name}")
                print("  BGL dimensions (index: value):",
                      {i: round(d, 6) for i, d in enumerate(dims) if d != 0.0})
        else:
            print("  D1,D2,D3      :", [sec[S.ipD1], sec[S.ipD2], sec[S.ipD3]])
            print("  T1,T2,T3      :", [sec[S.ipT1], sec[S.ipT2], sec[S.ipT3]])

        print(f"  A={sec[S.ipAREA]:.6g}  I11={sec[S.ipI11]:.6g}  I22={sec[S.ipI22]:.6g}  J={sec[S.ipJ]:.6g}")
        print(f"  SA1={sec[S.ipSA1]:.6g}  SA2={sec[S.ipSA2]:.6g}")
        print(f"  E={mat[S.ipModulus]:.6g}  nu={mat[S.ipPoisson]:.4g}  "
              f"density={mat[S.ipDensity]:.6g}  G={mat[S.ipShearModulus]:.6g}")
        print()
    return found


def list_beams(uid):
    total = ctypes.c_long()
    check(S.St7GetTotal(uid, S.tyBEAM, ctypes.byref(total)))
    print(f"Beams: {total.value}")
    print("-" * 70)
    print(f"{'Beam':>5} {'N1':>6} {'N2':>6} {'Prop':>5} {'Length':>12}  Group")

    for beam in range(1, total.value + 1):
        conn = (ctypes.c_long * (S.kMaxElementNode + 1))()
        check(S.St7GetElementConnection(uid, S.tyBEAM, beam, conn))

        prop = ctypes.c_long()
        check(S.St7GetElementProperty(uid, S.ptBEAMPROP, beam, ctypes.byref(prop)))

        length = ctypes.c_double()
        check(S.St7GetElementData(uid, S.tyBEAM, beam, 0, ctypes.byref(length)))

        group_id = ctypes.c_long()
        check(S.St7GetEntityGroup(uid, S.tyBEAM, beam, ctypes.byref(group_id)))
        group_name = get_str(S.St7GetGroupIDName, uid, group_id.value)

        print(f"{beam:>5} {conn[1]:>6} {conn[2]:>6} {prop.value:>5} "
              f"{length.value:>12.5g}  {group_name}")


def main():
    with model_session() as uid:
        print_units(uid)
        list_beam_properties(uid)
        list_beams(uid)


if __name__ == "__main__":
    main()