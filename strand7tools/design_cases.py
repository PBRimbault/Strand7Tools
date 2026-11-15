"""
Design cases (load combinations) for the EC3 checks.

design_cases.csv defines each design case as a factored sum of Strand7 result cases.
This module reads and validates that file; it does not need Strand7.
tools/show_inputs.py prints what it found.

Rules:
  - results from linear static (.LSA) files can be added together (superposition);
  - results from any other result file type (nonlinear etc.) cannot: such a design
    case must be a single row with factor 1;
  - Strand7's envelope cases are never used: the interaction checks need N, V and M
    that occur together in one real case or combination.
"""

import csv
import os
from dataclasses import dataclass, field

LINEAR_EXTENSIONS = {".lsa"}


class DesignCaseError(Exception):
    """A problem in design_cases.csv; the message says what to fix."""


@dataclass(frozen=True)
class Row:
    result_file: str | None      # None = the default result file; else an absolute path
    result_case: int
    factor: float
    comment: str = ""

    @property
    def linear(self):
        """Default result file is assumed to be a linear static file."""
        if self.result_file is None:
            return True
        return os.path.splitext(self.result_file)[1].lower() in LINEAR_EXTENSIONS


@dataclass
class DesignCase:
    name: str
    rows: list = field(default_factory=list)

    @property
    def is_direct(self):
        """A single row with factor 1: the Strand7 result is used as it is."""
        return len(self.rows) == 1 and self.rows[0].factor == 1.0

    @property
    def needs_superposition(self):
        return not self.is_direct


def load_design_cases(path):
    """
    Return {name: DesignCase} in file order, after validating every row. A relative
    result_file is taken relative to the folder of the design case file.
    """
    if not os.path.isfile(path):
        raise DesignCaseError(f"Design case file not found: {path}")
    with open(path, newline="", encoding="utf-8-sig") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    headers = {h.strip().lower() for h in (reader.fieldnames or [])}
    missing = {"design_case", "result_case", "factor"} - headers
    if missing:
        raise DesignCaseError(f"{path}: missing column(s): {', '.join(sorted(missing))}")

    cases = {}
    for n, raw in enumerate(reader, start=1):
        row = {k.strip().lower(): (v or "").strip() for k, v in raw.items() if k is not None}
        where = f"{path}: data row {n}"
        name = row.get("design_case", "")
        if not name:
            raise DesignCaseError(f"{where}: design_case is blank")
        try:
            case_number = int(row.get("result_case", ""))
        except ValueError:
            raise DesignCaseError(f"{where} ({name}): result_case must be a whole number, "
                                  f"got {row.get('result_case')!r}")
        if case_number < 1:
            raise DesignCaseError(f"{where} ({name}): result_case must be 1 or more")
        try:
            factor = float(row.get("factor", ""))
        except ValueError:
            raise DesignCaseError(f"{where} ({name}): factor must be a number, "
                                  f"got {row.get('factor')!r}")
        if factor == 0:
            raise DesignCaseError(f"{where} ({name}): a factor of zero is almost certainly a "
                                  f"mistake; delete the row instead")
        result_file = row.get("result_file") or None
        if result_file and not os.path.isabs(result_file):
            result_file = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(path)),
                                                        result_file))
        cases.setdefault(name, DesignCase(name)).rows.append(
            Row(result_file, case_number, factor, row.get("comment", "")))

    for case in cases.values():
        seen = set()
        for r in case.rows:
            key = (r.result_file and r.result_file.lower(), r.result_case)
            if key in seen:
                raise DesignCaseError(f"{path}: design case {case.name} uses result case "
                                      f"{r.result_case} more than once; combine the factors")
            seen.add(key)
        if any(not r.linear for r in case.rows) and not case.is_direct:
            raise DesignCaseError(
                f"{path}: design case {case.name} uses a result file that is not a linear "
                f"static (.LSA) file. Such results cannot be added together: it must be a "
                f"single row with factor 1.")
    return cases


def file_warnings(cases):
    """Result files named in the table that do not exist."""
    out = []
    for case in cases.values():
        for r in case.rows:
            if r.result_file and not os.path.isfile(r.result_file):
                out.append(f"{case.name}: result file not found: {r.result_file}")
    return out
