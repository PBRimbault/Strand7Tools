"""
Step 2b: Read beam forces from a linear static (.LSA) result file. Read-only.

Settings are at the top of this file:
  RESULT_FILE   - path to the .LSA (None = same name/folder as the model).
                  Set it to use a results file in a scratch folder.
  RESULT_CASES  - result case numbers to read. Leave empty [] to just LIST
                  the available cases and stop.
  COMBINATIONS  - "none", "existing" or "generate": whether Strand7's own load
                  case combinations are opened with the result file.
  N_STATIONS    - number of stations per beam (evenly spaced, end to end).

Writes every value to forces.csv and prints the maximum absolute values for
each beam and case, so you can check them against Strand7.

Forces are Strand7 PRINCIPAL-axis results: SF1, BM1, SF2, BM2, axial, torque.
The sign convention and which bending moment belongs to which axis must be
verified against a hand calculation before we rely on them (see notes).
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))  # st7_setup changes the working dir

import csv  # noqa: E402
import ctypes  # noqa: E402
from contextlib import contextmanager  # noqa: E402

from st7_setup import St7API as S, check, model_session, MODEL_FILE  # noqa: E402

# ---------------- SETTINGS (edit these) ----------------
RESULT_FILE = None          # e.g. r"D:\Scratch\test_model.LSA"
RESULT_CASES = [3]           # e.g. [1, 2, 3]; [] = list cases only
COMBINATIONS = "existing"   # "none", "existing" or "generate"
N_STATIONS = 5
OUTPUT_CSV = "forces.csv"
PRINT_BEAMS = [10]            # e.g. [10] to print the signed value at every station for these beams
# -------------------------------------------------------

COLUMNS = ["SF1", "BM1", "SF2", "BM2", "Axial", "Torque"]


def get_str(fn, *args):
    buf = ctypes.create_string_buffer(S.kMaxStrLen)
    check(fn(*args, buf, S.kMaxStrLen))
    return buf.value.decode(errors="replace")


@contextmanager
def result_file(uid, path, combination_code):
    n_primary = ctypes.c_long()
    n_secondary = ctypes.c_long()
    check(S.St7OpenResultFile(uid, path.encode(), b"", combination_code,
                              ctypes.byref(n_primary), ctypes.byref(n_secondary)))
    try:
        yield n_primary.value, n_secondary.value
    finally:
        S.St7CloseResultFile(uid)


def force_scale(uid):
    """Return (force factor to kN, label) if the model is in N and m."""
    units = (ctypes.c_long * S.kLastUnit)()
    check(S.St7GetUnits(uid, units))
    if units[S.ipLENGTHU] == S.luMETRE and units[S.ipFORCEU] == S.fuNEWTON:
        return 1e-3, "kN, kNm"
    print("WARNING: model is not in N and m; values are left in model units.")
    return 1.0, "model units"


def main():
    result_path = RESULT_FILE or (os.path.splitext(MODEL_FILE)[0] + ".LSA")
    if not os.path.isfile(result_path):
        raise FileNotFoundError(f"Result file not found: {result_path}")

    comb_code = {
        "none": S.kNoCombinations,
        "existing": S.kUseExistingCombinations,
        "generate": S.kGenerateNewCombinations,
    }[COMBINATIONS]

    with model_session() as uid:
        scale, unit_label = force_scale(uid)

        with result_file(uid, result_path, comb_code) as (n_prim, n_sec):
            print(f"Result file : {result_path}")
            print(f"Cases       : {n_prim} primary, {n_sec} secondary (combinations)")
            print()
            print("Available result cases:")
            for c in range(1, n_prim + n_sec + 1):
                kind = "primary" if c <= n_prim else "combination"
                print(f"  {c:>3}  {get_str(S.St7GetResultCaseName, uid, c)}   [{kind}]")
            print()

            if not RESULT_CASES:
                print("RESULT_CASES is empty: edit it at the top of this file "
                      "(for example RESULT_CASES = [1, 2]) and run again.")
                return

            # Stations are given as a fraction of beam length (0..1)
            check(S.St7SetBeamResultPosMode(uid, S.bpParam))
            fractions = [i / (N_STATIONS - 1) for i in range(N_STATIONS)]
            pos = (ctypes.c_double * S.kMaxBeamResult)()
            for i, p in enumerate(fractions):
                pos[i] = p
            res = (ctypes.c_double * S.kMaxBeamResult)()
            ncols = ctypes.c_long()
            idx = [S.ipBeamSF1, S.ipBeamBM1, S.ipBeamSF2,
                   S.ipBeamBM2, S.ipBeamAxialF, S.ipBeamTorque]

            n_beams = ctypes.c_long()
            check(S.St7GetTotal(uid, S.tyBEAM, ctypes.byref(n_beams)))

            rows = []
            summary = []
            for case in RESULT_CASES:
                case_name = get_str(S.St7GetResultCaseName, uid, case)
                for beam in range(1, n_beams.value + 1):
                    prop = ctypes.c_long()
                    check(S.St7GetElementProperty(uid, S.ptBEAMPROP, beam,
                                                  ctypes.byref(prop)))
                    length = ctypes.c_double()
                    check(S.St7GetElementData(uid, S.tyBEAM, beam, 0,
                                              ctypes.byref(length)))
                    check(S.St7GetBeamResultArrayPos(
                        uid, S.rtBeamForce, S.stBeamPrincipal, beam, case,
                        N_STATIONS, pos, ctypes.byref(ncols), res))

                    maxabs = [0.0] * len(COLUMNS)
                    for i, frac in enumerate(fractions):
                        vals = [res[i * ncols.value + k] * scale for k in idx]
                        rows.append([case, case_name, beam, prop.value, frac,
                                     frac * length.value] + vals)
                        maxabs = [max(m, abs(v)) for m, v in zip(maxabs, vals)]
                    summary.append((case, beam, maxabs))

            out_path = os.path.join(HERE, OUTPUT_CSV)
            with open(out_path, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["case", "case_name", "beam", "prop", "position_fraction",
                            "x_m"] + COLUMNS)
                w.writerows(rows)
            print(f"Wrote {len(rows)} rows to {out_path}")
            print()

            print(f"Maximum absolute values per beam ({unit_label}):")
            print(f"{'Case':>5} {'Beam':>5} " + " ".join(f"{c:>10}" for c in COLUMNS))
            for case, beam, maxabs in summary:
                print(f"{case:>5} {beam:>5} " + " ".join(f"{m:>10.3f}" for m in maxabs))

            if PRINT_BEAMS:
                print()
                print(f"Signed values at each station ({unit_label}); x = distance from end 1:")
                print(f"{'Case':>5} {'Beam':>5} {'x':>8} "
                      + " ".join(f"{c:>10}" for c in COLUMNS))
                for r in rows:
                    if r[2] in PRINT_BEAMS:
                        print(f"{r[0]:>5} {r[2]:>5} {r[5]:>8.3f} "
                              + " ".join(f"{v:>10.3f}" for v in r[6:]))


if __name__ == "__main__":
    main()