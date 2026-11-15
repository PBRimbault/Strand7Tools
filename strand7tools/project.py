"""
A project file (TOML) names one Strand7 model and the inputs for checking it. Paths in it are
relative to the project file's folder. Example (Test/test_model2.toml):

    model = "test_model2.st7"          # required
    members = "members.csv"            # required: the member table
    design_cases = "design_cases.csv"  # required: the design cases

    # optional, with their defaults:
    result_file = ""                   # default result file; "" = the model's name with .LSA
    combinations = "existing"          # Strand7 combination cases: "none", "existing", "generate"
    run_cases = []                     # names of design cases to run; [] = all of them
    stations_per_beam = 21             # stations per beam when cases are added together
    settings = ""                      # "" = strand7tools/data/settings.toml
    steel_grades = ""                  # "" = strand7tools/data/steel_grades.csv
    output_folder = "output"           # where the result CSV files are written
"""

import os
import tomllib
from dataclasses import dataclass, field

from .settings import GRADES_FILE, SETTINGS_FILE

COMBINATION_MODES = ("none", "existing", "generate")


class ProjectError(Exception):
    """A problem in the project file; the message says what to fix."""


@dataclass(frozen=True)
class Project:
    path: str                   # the project file
    model: str                  # absolute paths from here on
    members: str
    design_cases: str
    result_file: str
    combinations: str = "existing"
    run_cases: tuple = field(default_factory=tuple)
    stations_per_beam: int = 21
    settings: str = SETTINGS_FILE
    steel_grades: str = GRADES_FILE
    output_folder: str = ""

    @property
    def name(self):
        return os.path.splitext(os.path.basename(self.path))[0]

    @property
    def folder(self):
        return os.path.dirname(self.path)

    def output(self, filename):
        """Absolute path of an output file (the output folder is created if needed)."""
        os.makedirs(self.output_folder, exist_ok=True)
        return os.path.join(self.output_folder, filename)


_KEYS = {"model", "members", "design_cases", "result_file", "combinations", "run_cases",
         "stations_per_beam", "settings", "steel_grades", "output_folder"}


def load_project(path):
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        raise ProjectError(f"Project file not found: {path}")
    with open(path, "rb") as f:
        try:
            data = tomllib.load(f)
        except tomllib.TOMLDecodeError as e:
            raise ProjectError(f"{path}: not valid TOML ({e})")
    extras = set(data) - _KEYS
    if extras:
        raise ProjectError(f"{path}: unknown key(s) {', '.join(sorted(extras))} (typo?)")
    folder = os.path.dirname(path)

    def text(key, default=None):
        value = data.get(key, default)
        if value is None:
            raise ProjectError(f"{path}: '{key}' is required")
        if not isinstance(value, str):
            raise ProjectError(f"{path}: '{key}' must be text in quotes")
        return value

    def file(key, default=None):
        value = text(key, default)
        return os.path.normpath(os.path.join(folder, value)) if value else ""

    model = file("model")
    result_file = file("result_file", "") or os.path.splitext(model)[0] + ".LSA"
    combinations = text("combinations", "existing")
    if combinations not in COMBINATION_MODES:
        raise ProjectError(f"{path}: combinations must be one of {', '.join(COMBINATION_MODES)}")
    run_cases = data.get("run_cases", [])
    if not isinstance(run_cases, list) or not all(isinstance(c, str) for c in run_cases):
        raise ProjectError(f"{path}: run_cases must be a list of names, e.g. [\"ULS_A\"]")
    stations = data.get("stations_per_beam", 21)
    if not isinstance(stations, int) or isinstance(stations, bool) or stations < 3:
        raise ProjectError(f"{path}: stations_per_beam must be a whole number, at least 3")
    return Project(
        path=path, model=model, members=file("members"), design_cases=file("design_cases"),
        result_file=result_file, combinations=combinations, run_cases=tuple(run_cases),
        stations_per_beam=stations,
        settings=file("settings", "") or SETTINGS_FILE,
        steel_grades=file("steel_grades", "") or GRADES_FILE,
        output_folder=file("output_folder", "output"),
    )
