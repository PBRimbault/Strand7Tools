"""
Loads the design settings (settings.toml) and the steel grade table (steel_grades.csv).

Run this file on its own to print what was loaded, with a small table of fy and fu by
thickness, so you can check the values. It does not need Strand7.
"""

import csv
import math
import os
import tomllib
from dataclasses import dataclass

HERE = os.path.dirname(os.path.abspath(__file__))
SETTINGS_FILE = os.path.join(HERE, "settings.toml")
GRADES_FILE = os.path.join(HERE, "steel_grades.csv")


class SettingsError(Exception):
    """A settings or grade file problem; the message says what to fix."""


@dataclass(frozen=True)
class Settings:
    material_source: str      # "model" or "ec3"
    E_ec3: float              # N/mm2
    G_ec3: float              # N/mm2
    default_grade: str
    thickness_rule: str       # "max"
    gamma_M0: float
    gamma_M1: float
    gamma_M2: float
    eta: float                # lower limit factor for Av,z


def _section(data, name, path):
    if name not in data or not isinstance(data[name], dict):
        raise SettingsError(f"{path}: missing section [{name}]")
    return data[name]


def _take(section, section_name, key, kind, path):
    if key not in section:
        raise SettingsError(f"{path}: [{section_name}] is missing '{key}'")
    value = section[key]
    if kind is float and isinstance(value, int) and not isinstance(value, bool):
        value = float(value)
    if not isinstance(value, kind):
        raise SettingsError(f"{path}: [{section_name}] {key} must be a {kind.__name__}, "
                            f"got {value!r}")
    return value


def _no_extras(section, section_name, allowed, path):
    extras = set(section) - set(allowed)
    if extras:
        raise SettingsError(f"{path}: unknown setting(s) in [{section_name}]: "
                            f"{', '.join(sorted(extras))} (typo?)")


def load_settings(path=SETTINGS_FILE):
    if not os.path.isfile(path):
        raise SettingsError(f"Settings file not found: {path}")
    with open(path, "rb") as f:
        try:
            data = tomllib.load(f)
        except tomllib.TOMLDecodeError as e:
            raise SettingsError(f"{path}: not valid TOML ({e})")

    mat = _section(data, "material", path)
    pf = _section(data, "partial_factors", path)
    sh = _section(data, "shear", path)
    _no_extras(mat, "material", ["source", "E_ec3", "G_ec3", "default_grade", "thickness_rule"], path)
    _no_extras(pf, "partial_factors", ["gamma_M0", "gamma_M1", "gamma_M2"], path)
    _no_extras(sh, "shear", ["eta"], path)

    s = Settings(
        material_source=_take(mat, "material", "source", str, path),
        E_ec3=_take(mat, "material", "E_ec3", float, path),
        G_ec3=_take(mat, "material", "G_ec3", float, path),
        default_grade=_take(mat, "material", "default_grade", str, path),
        thickness_rule=_take(mat, "material", "thickness_rule", str, path),
        gamma_M0=_take(pf, "partial_factors", "gamma_M0", float, path),
        gamma_M1=_take(pf, "partial_factors", "gamma_M1", float, path),
        gamma_M2=_take(pf, "partial_factors", "gamma_M2", float, path),
        eta=_take(sh, "shear", "eta", float, path),
    )
    if s.material_source not in ("model", "ec3"):
        raise SettingsError(f"{path}: [material] source must be \"model\" or \"ec3\"")
    if s.thickness_rule != "max":
        raise SettingsError(f"{path}: [material] thickness_rule must be \"max\" (the only rule so far)")
    for name in ("E_ec3", "G_ec3", "gamma_M0", "gamma_M1", "gamma_M2", "eta"):
        if getattr(s, name) <= 0:
            raise SettingsError(f"{path}: {name} must be positive")
    return s


def load_grades(path=GRADES_FILE):
    """{grade: [(t_max_mm, fy, fu), ...]} with each grade's bands in ascending thickness."""
    if not os.path.isfile(path):
        raise SettingsError(f"Steel grade table not found: {path}")
    with open(path, newline="", encoding="utf-8-sig") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    grades = {}
    for n, row in enumerate(csv.DictReader(lines), start=1):
        try:
            grade = row["grade"].strip()
            band = (float(row["t_max_mm"]), float(row["fy"]), float(row["fu"]))
        except (KeyError, ValueError, AttributeError):
            raise SettingsError(f"{path}: row {n} is not 'grade,t_max_mm,fy,fu' with numbers")
        if not grade or band[0] <= 0 or band[1] <= 0 or band[2] <= 0:
            raise SettingsError(f"{path}: row {n} has a blank grade or non-positive numbers")
        grades.setdefault(grade, []).append(band)
    if not grades:
        raise SettingsError(f"{path}: no grades found")
    for grade, bands in grades.items():
        thicknesses = [b[0] for b in bands]
        if thicknesses != sorted(thicknesses) or len(set(thicknesses)) != len(thicknesses):
            raise SettingsError(f"{path}: grade {grade} thickness limits must be ascending "
                                f"and not repeated")
    return grades


def strength(grades, grade, thickness_mm):
    """(fy, fu) in N/mm2 for a grade at a thickness in mm."""
    if grade not in grades:
        raise SettingsError(f"unknown steel grade {grade!r}; the table has: "
                            f"{', '.join(grades)}")
    bands = grades[grade]
    for t_max, fy, fu in bands:
        if thickness_mm <= t_max:
            return fy, fu
    raise SettingsError(f"no strength data for {grade} at {thickness_mm:g} mm "
                        f"(the table goes up to {bands[-1][0]:g} mm)")


def eng(x, decimals=3):
    """Engineering notation, e.g. 210.000e+03."""
    if x == 0:
        return f"{0:.{decimals}f}e+00"
    exp = int(math.floor(math.log10(abs(x)) / 3) * 3)
    return f"{x / 10 ** exp:.{decimals}f}e{exp:+03d}"


def main():
    s = load_settings()
    grades = load_grades()

    print(f"Settings file : {SETTINGS_FILE}")
    print(f"  material source : {s.material_source}")
    print(f"  E (if 'ec3')    : {eng(s.E_ec3)} N/mm2")
    print(f"  G (if 'ec3')    : {eng(s.G_ec3)} N/mm2")
    print(f"  default grade   : {s.default_grade}")
    print(f"  thickness rule  : {s.thickness_rule}")
    print(f"  gamma_M0, M1, M2: {s.gamma_M0:g}, {s.gamma_M1:g}, {s.gamma_M2:g}")
    print(f"  eta (Av,z)      : {s.eta:g}")
    if s.default_grade not in grades:
        print(f"  WARNING: default grade {s.default_grade!r} is not in the grade table")
    print()

    thicknesses = [10, 16, 20, 40, 50, 63, 80, 100]
    print(f"Grade table   : {GRADES_FILE}")
    print("  fy / fu in N/mm2 by thickness (mm)")
    print(f"  {'grade':<7}" + "".join(f"{t:>10}" for t in thicknesses))
    for grade in grades:
        cells = []
        for t in thicknesses:
            try:
                fy, fu = strength(grades, grade, t)
                cells.append(f"{fy:g}/{fu:g}")
            except SettingsError:
                cells.append("-")
        print(f"  {grade:<7}" + "".join(f"{c:>10}" for c in cells))


if __name__ == "__main__":
    main()
