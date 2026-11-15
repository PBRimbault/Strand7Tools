"""
Shared Strand7 setup. Every script imports this instead of repeating paths.

Usage in another script:

    import st7_setup
    from st7_setup import St7API, check, model_session

    with model_session() as uid:
        ...  # model is open here; it is closed and released automatically
"""

import ctypes
import os
import sys
import tempfile
from contextlib import contextmanager

# ---------------- SETTINGS (edit these) ----------------

# Folder containing St7API.py
API_FOLDER = r"C:\Program Files\Strand7 R31\API Includes\Python"

# Folder containing St7API.dll
DLL_FOLDER = r"C:\Program Files\Strand7 R31\Bin64"

# The test model (.st7). Update if you saved your 10-beam model elsewhere.
MODEL_FILE = r"C:\Strand7Tools\Test\test_model2.st7"

# Scratch folder for Strand7 temporary files. Must already exist.
SCRATCH_FOLDER = tempfile.gettempdir()

# -------------------------------------------------------

os.add_dll_directory(DLL_FOLDER)
os.chdir(DLL_FOLDER)  # fallback so dependent DLLs are found
sys.path.append(API_FOLDER)

import St7API  # noqa: E402  (must come after the path setup above)

UID = 1  # model ID number used for the open model


def check(error_code):
    """Raise a readable error if a Strand7 API call failed."""
    if error_code != St7API.ERR7_NoError:
        buf = ctypes.create_string_buffer(St7API.kMaxStrLen)
        St7API.St7GetAPIErrorString(error_code, buf, St7API.kMaxStrLen)
        raise RuntimeError(f"Strand7 error {error_code}: {buf.value.decode()}")


@contextmanager
def model_session(model_file=None):
    """Start the API, open the model, and always clean up afterwards."""
    model_file = model_file or MODEL_FILE
    if not os.path.isfile(model_file):
        raise FileNotFoundError(f"Model file not found: {model_file}")

    check(St7API.St7Init())
    try:
        check(St7API.St7OpenFile(UID, model_file.encode(), SCRATCH_FOLDER.encode()))
        try:
            yield UID
        finally:
            St7API.St7CloseFile(UID)
    finally:
        St7API.St7Release()
