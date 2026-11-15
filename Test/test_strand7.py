"""
Step 1: Check that Python can talk to Strand7.

What it does: opens a model, counts nodes and beams, closes the model.
It does NOT save or change your model.

Before running, edit the three settings in the SETTINGS block below.
"""

import ctypes
import os
import sys
import tempfile

# ---------------- SETTINGS (edit these) ----------------

# Folder that contains St7API.py (the Strand7 Python wrapper).
# Typical location (check your install, the exact path may differ):
API_FOLDER = r"C:\Program Files\Strand7 R31\API Includes\Python"

# Full path to a SMALL test model (.st7). Use a copy, not your only version.
MODEL_FILE = r"C:\Strand7Tools\Test\test_model.st7"

# Scratch folder Strand7 uses for temporary files. Must already exist.
SCRATCH_FOLDER = tempfile.gettempdir()

# -------------------------------------------------------

# Folder containing St7API.dll (and the other Strand7 DLLs)
DLL_FOLDER = r"C:\Program Files\Strand7 R31\Bin64"   # <-- check this path

sys.path.append(API_FOLDER)

os.add_dll_directory(DLL_FOLDER)
os.chdir(DLL_FOLDER)  # fallback in case dependencies are still not found

try:
    import St7API
except ImportError:
    print("Could not import St7API.")
    print(f"I looked in: {API_FOLDER}")
    print("Find St7API.py on your machine and put its folder in API_FOLDER.")
    sys.exit(1)

UID = 1  # model ID number; any small integer works when one model is open


def check(error_code):
    """Stop with a readable message if a Strand7 call fails."""
    if error_code != St7API.ERR7_NoError:
        buf = ctypes.create_string_buffer(St7API.kMaxStrLen)
        St7API.St7GetAPIErrorString(error_code, buf, St7API.kMaxStrLen)
        raise RuntimeError(f"Strand7 error {error_code}: {buf.value.decode()}")


def main():
    if not os.path.isfile(MODEL_FILE):
        print(f"Model file not found: {MODEL_FILE}")
        return

    check(St7API.St7Init())
    print("Strand7 API started.")

    try:
        check(St7API.St7OpenFile(UID, MODEL_FILE.encode(), SCRATCH_FOLDER.encode()))
        print(f"Opened: {MODEL_FILE}")

        n_nodes = ctypes.c_long()
        n_beams = ctypes.c_long()
        check(St7API.St7GetTotal(UID, St7API.tyNODE, ctypes.byref(n_nodes)))
        check(St7API.St7GetTotal(UID, St7API.tyBEAM, ctypes.byref(n_beams)))

        print(f"Nodes: {n_nodes.value}")
        print(f"Beams: {n_beams.value}")

        check(St7API.St7CloseFile(UID))
    finally:
        St7API.St7Release()
        print("Strand7 API released.")


if __name__ == "__main__":
    main()