"""
The member table (members.csv): which beams make up each design member, and its buckling
lengths. Pure parsing, no Strand7 (strand7/members.py builds the members from the model).

Terms:
  beam   - a Strand7 beam element
  member - a continuous, straight, constant-section chain of beams over which
           buckling can occur (from one restraint to the next)
"""

import csv
import os
from dataclasses import dataclass

YES = {"y", "yes", "true", "1", "x"}
NO = {"n", "no", "false", "0"}


def parse_bool(text, default=True):
    t = (text or "").strip().lower()
    if t == "":
        return default
    if t in YES:
        return True
    if t in NO:
        return False
    raise ValueError(f"cannot read {text!r} as yes/no")


def parse_float(text, default=None):
    t = (text or "").strip()
    return default if t == "" else float(t)


def parse_int_list(text):
    """'1-4;7' -> [1, 2, 3, 4, 7]  (';' or ',' separated)."""
    out = []
    for tok in (text or "").replace(",", ";").split(";"):
        tok = tok.strip()
        if not tok:
            continue
        if "-" in tok:
            a, b = (int(x) for x in tok.split("-", 1))
            if b < a:
                raise ValueError(f"bad range {tok!r}")
            out.extend(range(a, b + 1))
        else:
            out.append(int(tok))
    return out


def resolve_beams(text, beam_ids, n_beams):
    """Beam numbers from text like '1-10', '1;2;5' or 'ID=5' (all beams with that Beam ID)."""
    beams = []
    for tok in (text or "").replace(",", ";").split(";"):
        tok = tok.strip()
        if not tok:
            continue
        if tok.upper().startswith("ID="):
            wanted = int(tok[3:])
            found = [b for b, i in beam_ids.items() if i == wanted]
            if not found:
                raise ValueError(f"no beams have Beam ID {wanted}")
            beams.extend(found)
        else:
            beams.extend(parse_int_list(tok))
    unique = list(dict.fromkeys(beams))
    bad = [b for b in unique if not 1 <= b <= n_beams]
    if bad:
        raise ValueError(f"beam numbers not in the model: {bad}")
    return unique


@dataclass
class MemberDef:
    member_id: str
    beams: list
    buckling_y: bool = True
    kcr_y: float = 1.0
    lcr_y: float = None  # optional override of kcr_y x member length
    buckling_z: bool = True
    kcr_z: float = 1.0
    lcr_z: float = None
    ltb: bool = True
    kz: float = 1.0
    kw: float = 1.0
    lw: float = None     # optional override of kz x member length
    grade: str = ""      # steel grade, e.g. S355; blank = the default grade in settings.toml
    comment: str = ""

    def effective_lengths(self, length):
        """{'Lcr_y', 'Lcr_z', 'Lw'} in m for a member of this length (None = switched off)."""
        def pick(on, k, override):
            if not on:
                return None
            return override if override is not None else k * length

        return {
            "Lcr_y": pick(self.buckling_y, self.kcr_y, self.lcr_y),
            "Lcr_z": pick(self.buckling_z, self.kcr_z, self.lcr_z),
            "Lw": pick(self.ltb, self.kz, self.lw),
        }


def read_rows(path):
    """The table's rows as dicts (lower-case keys, stripped values); '#' lines are skipped."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Member table not found: {path}")
    with open(path, newline="", encoding="utf-8-sig") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    rows = []
    for raw in csv.DictReader(lines):
        rows.append({k.strip().lower(): (v or "").strip()
                     for k, v in raw.items() if k is not None})
    return rows


def parse_row(row, beam_ids, n_beams):
    member_id = row.get("member", "")
    if not member_id:
        raise ValueError("missing member number")
    beams = resolve_beams(row.get("beams", ""), beam_ids, n_beams)
    if not beams:
        raise ValueError("no beams listed")
    return MemberDef(
        member_id=member_id,
        beams=beams,
        buckling_y=parse_bool(row.get("buckling_y")),
        kcr_y=parse_float(row.get("kcr_y"), 1.0),
        lcr_y=parse_float(row.get("lcr_y")),
        buckling_z=parse_bool(row.get("buckling_z")),
        kcr_z=parse_float(row.get("kcr_z"), 1.0),
        lcr_z=parse_float(row.get("lcr_z")),
        ltb=parse_bool(row.get("ltb")),
        kz=parse_float(row.get("kz"), 1.0),
        kw=parse_float(row.get("kw"), 1.0),
        lw=parse_float(row.get("lw")),
        grade=row.get("grade", ""),
        comment=row.get("comment", ""),
    )
