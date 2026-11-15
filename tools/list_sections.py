"""
EC3 section properties for every beam property of the model, next to Strand7's own A, Iy and
Iz as a built-in check (the values can also be compared with the published section tables).
Read-only.

    python tools\\list_sections.py [project file]     (default: Test\\test_model2.toml)

Standard I-section properties (no root radius in the data) are computed with r = 0, with a
warning: Wpl and Av come out slightly low.
"""

import ctypes
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from strand7tools.fmt import eng  # noqa: E402
from strand7tools.project import load_project  # noqa: E402
from strand7tools.settings import load_settings  # noqa: E402
from strand7tools.strand7.api import S, check, model_session  # noqa: E402
from strand7tools.strand7.sections import SectionError, make_section, read_section  # noqa: E402

DEFAULT_PROJECT = os.path.join(ROOT, "Test", "test_model2.toml")

# (label, unit, property, factor from SI to the unit, Strand7 key or None); mm, y before z
ROWS = [("A", "mm2", "A", 1e6, "A"), ("Iy", "mm4", "Iy", 1e12, "Iy"),
        ("Iz", "mm4", "Iz", 1e12, "Iz"), ("Wel,y", "mm3", "Wel_y", 1e9, None),
        ("Wel,z", "mm3", "Wel_z", 1e9, None), ("Wpl,y", "mm3", "Wpl_y", 1e9, None),
        ("Wpl,z", "mm3", "Wpl_z", 1e9, None), ("Av,y", "mm2", "Av_y", 1e6, None),
        ("Av,z", "mm2", "Av_z", 1e6, None), ("It", "mm4", "It", 1e12, None),
        ("Iw", "mm6", "Iw", 1e18, None)]


def main(project_file=DEFAULT_PROJECT):
    project = load_project(project_file)
    eta = load_settings(project.settings).eta
    with model_session(project.model) as uid:
        n_props = (ctypes.c_long * S.kMaxEntityTotals)()
        last_prop = (ctypes.c_long * S.kMaxEntityTotals)()
        check(S.St7GetTotalProperties(uid, n_props, last_prop))

        for prop in range(1, last_prop[S.ipBeamPropTotal] + 1):
            try:
                s = read_section(uid, prop)
                sec, warnings = make_section(s, eta)
            except SectionError as e:
                print(f"Property {prop}: skipped - {e}\n")
                continue
            except RuntimeError:
                continue  # no property with this number

            P = sec.props
            print("=" * 70)
            print(f"Property {prop}: {s['name']}   ->   {sec.name}")
            print(f"  {s['kind']}")
            print(f"  h = {s['h']*1e3:.1f}  b = {s['b']*1e3:.1f}  tw = {s['tw']*1e3:.1f}  "
                  f"tf = {s['tf']*1e3:.1f}  r = {s['r']*1e3:.1f}   (mm)   eta = {eta:g}")
            for w in warnings:
                print(f"  WARNING: {w}")
            print(f"  {'':8} {'unit':>5} {'this app':>13} {'Strand7':>13} {'diff':>8}")
            for label, unit, key, fac, s7key in ROWS:
                val = getattr(P, key) * fac
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
            if sec.shape == "I":
                print(f"  classification geometry: flange outstand c/tf = "
                      f"{sec.c_flange / sec.tf:.2f}   web c/tw = {sec.c_web / sec.tw:.2f}")
            print(f"  material: E = {s['E'] / 1e9:.1f} GPa   G = {s['G'] / 1e9:.1f} GPa")
            print()


if __name__ == "__main__":
    main(*sys.argv[1:2])
