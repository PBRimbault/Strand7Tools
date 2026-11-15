"""
Loads the design settings (settings.toml) and the steel grade table (steel_grades.csv).
No Strand7 needed. The shipped defaults are in strand7tools/data; a project file can point
to its own copies. tools/show_inputs.py prints what was loaded.
"""

import csv
import os
import tomllib
from dataclasses import dataclass

DATA_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
SETTINGS_FILE = os.path.join(DATA_FOLDER, "settings.toml")
GRADES_FILE = os.path.join(DATA_FOLDER, "steel_grades.csv")


class SettingsError(Exception):
    """A settings or grade file problem; the message says what to fix."""


@dataclass(frozen=True)
class Buckling:
    """Member buckling parameters (EN 1993-1-1 6.3). Everything a National Annex may change."""
    alpha: dict               # imperfection factors, flexural buckling (Table 6.1)
    alpha_LT: dict            # imperfection factors, lateral-torsional buckling (Table 6.3)
    lambda_0: float           # plateau length, flexural buckling
    lambda_LT0: float         # plateau length, LTB (general case)
    beta_LT: float            # LTB (general case)
    high_strength_above: float  # fy above this (N/mm2) uses the S460 columns of Table 6.2
    flexural_rows: tuple      # rolled I-sections (Table 6.2): first matching row wins
    ltb_rows: tuple           # rolled I-sections (Table 6.4): first matching row wins


DEFAULT_BUCKLING = {
    "lambda_0": 0.2, "lambda_LT0": 0.2, "beta_LT": 1.0, "high_strength_above": 420.0,
    "alpha": {"a0": 0.13, "a": 0.21, "b": 0.34, "c": 0.49, "d": 0.76},
    "alpha_LT": {"a": 0.21, "b": 0.34, "c": 0.49, "d": 0.76},
    "flexural_curves": [
        {"hb_over": 1.2, "tf_upto": 40, "y": "a", "z": "b", "y_hs": "a0", "z_hs": "a0"},
        {"hb_over": 1.2, "tf_upto": 100, "y": "b", "z": "c", "y_hs": "a", "z_hs": "a"},
        {"hb_upto": 1.2, "tf_upto": 100, "y": "b", "z": "c", "y_hs": "a", "z_hs": "a"},
        {"y": "d", "z": "d", "y_hs": "c", "z_hs": "c"},
    ],
    "ltb_curves": [
        {"hb_upto": 2.0, "curve": "a"},
        {"hb_over": 2.0, "curve": "b"},
    ],
}


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
    buckling: Buckling = None


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


_ROW_KEYS_FLEX = {"hb_over", "hb_upto", "tf_upto", "y", "z", "y_hs", "z_hs"}
_ROW_KEYS_LTB = {"hb_over", "hb_upto", "tf_upto", "curve"}


def _load_buckling(data, path):
    raw = data.get("buckling", {})
    if not isinstance(raw, dict):
        raise SettingsError(f"{path}: [buckling] must be a section")
    allowed = set(DEFAULT_BUCKLING)
    extras = set(raw) - allowed
    if extras:
        raise SettingsError(f"{path}: unknown setting(s) in [buckling]: "
                            f"{', '.join(sorted(extras))} (typo?)")
    cfg = {k: raw.get(k, v) for k, v in DEFAULT_BUCKLING.items()}
    for key in ("lambda_0", "lambda_LT0", "beta_LT", "high_strength_above"):
        if not isinstance(cfg[key], (int, float)) or isinstance(cfg[key], bool) or cfg[key] <= 0:
            raise SettingsError(f"{path}: [buckling] {key} must be a positive number")
    for key in ("alpha", "alpha_LT"):
        table = dict(DEFAULT_BUCKLING[key])
        table.update(cfg[key] if isinstance(cfg[key], dict) else {})
        if any(not isinstance(v, (int, float)) or v < 0 for v in table.values()):
            raise SettingsError(f"{path}: [buckling.{key}] values must be numbers >= 0")
        cfg[key] = table

    def rows(name, allowed_keys, needed):
        out = []
        for i, r in enumerate(cfg[name], start=1):
            if not isinstance(r, dict) or set(r) - allowed_keys or not needed <= set(r):
                raise SettingsError(f"{path}: [[buckling.{name}]] row {i} must use only "
                                    f"{', '.join(sorted(allowed_keys))} and include "
                                    f"{', '.join(sorted(needed))}")
            out.append(dict(r))
        if not out:
            raise SettingsError(f"{path}: [[buckling.{name}]] has no rows")
        return tuple(out)

    flex = rows("flexural_curves", _ROW_KEYS_FLEX, {"y", "z"})
    ltb = rows("ltb_curves", _ROW_KEYS_LTB, {"curve"})
    for r in flex:
        for k in ("y", "z", "y_hs", "z_hs"):
            if k in r and r[k] not in cfg["alpha"]:
                raise SettingsError(f"{path}: curve {r[k]!r} has no imperfection factor in "
                                    f"[buckling.alpha]")
    for r in ltb:
        if r["curve"] not in cfg["alpha_LT"]:
            raise SettingsError(f"{path}: curve {r['curve']!r} has no imperfection factor in "
                                f"[buckling.alpha_LT]")
    return Buckling(cfg["alpha"], cfg["alpha_LT"], float(cfg["lambda_0"]),
                    float(cfg["lambda_LT0"]), float(cfg["beta_LT"]),
                    float(cfg["high_strength_above"]), flex, ltb)


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
        buckling=_load_buckling(data, path),
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
