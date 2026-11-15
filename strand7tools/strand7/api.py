"""
Access to the Strand7 API. The API (St7API.py and St7API.dll) is loaded the first time it is
used, so the rest of the app (and the tests) can be imported without Strand7.

    from strand7tools.strand7.api import S, check, model_session

    with model_session(r"C:\\path\\model.st7") as uid:
        check(S.St7GetTotal(uid, S.tyBEAM, ...))   # the model is closed again afterwards

Models are only ever opened for reading here: nothing saves them.
"""

import ctypes
import os
import sys
import tempfile
from contextlib import contextmanager

# ---------------- SETTINGS (edit these for your Strand7 installation) ----------------
API_FOLDER = r"C:\Program Files\Strand7 R31\API Includes\Python"     # contains St7API.py
DLL_FOLDER = r"C:\Program Files\Strand7 R31\Bin64"                   # contains St7API.dll
SCRATCH_FOLDER = tempfile.gettempdir()      # Strand7 temporary files; must already exist
# --------------------------------------------------------------------------------------

UID = 1  # model ID number used for the open model

_module = None


def load_api():
    """Import St7API (once). St7API.py loads the DLL by name, so the DLL folder must be the
    working directory at that moment; it stays so, as Strand7 may need it later."""
    global _module
    if _module is None:
        if not os.path.isdir(DLL_FOLDER) or not os.path.isdir(API_FOLDER):
            raise RuntimeError(f"Strand7 API not found: check API_FOLDER and DLL_FOLDER in {__file__}")
        os.add_dll_directory(DLL_FOLDER)
        os.chdir(DLL_FOLDER)
        if API_FOLDER not in sys.path:
            sys.path.append(API_FOLDER)
        import St7API
        _module = St7API
    return _module


class _LazyAPI:
    """Stands in for the St7API module; loads it on first attribute access."""

    def __getattr__(self, name):
        return getattr(load_api(), name)


S = _LazyAPI()


def check(error_code):
    """Raise a readable error if a Strand7 API call failed."""
    if error_code != S.ERR7_NoError:
        buf = ctypes.create_string_buffer(S.kMaxStrLen)
        S.St7GetAPIErrorString(error_code, buf, S.kMaxStrLen)
        raise RuntimeError(f"Strand7 error {error_code}: {buf.value.decode()}")


@contextmanager
def model_session(model_file):
    """Start the API, open the model, and always clean up afterwards."""
    model_file = os.path.abspath(model_file)
    if not os.path.isfile(model_file):
        raise FileNotFoundError(f"Model file not found: {model_file}")
    check(S.St7Init())
    try:
        check(S.St7OpenFile(UID, model_file.encode(), SCRATCH_FOLDER.encode()))
        try:
            yield UID
        finally:
            S.St7CloseFile(UID)
    finally:
        S.St7Release()


def get_str(fn, *args):
    """A string from an API function that fills a buffer."""
    buf = ctypes.create_string_buffer(S.kMaxStrLen)
    check(fn(*args, buf, S.kMaxStrLen))
    return buf.value.decode("cp1252", errors="replace")  # Strand7 strings are Windows ANSI


def model_units(uid):
    units = (ctypes.c_long * S.kLastUnit)()
    check(S.St7GetUnits(uid, units))
    return units


def total_beams(uid):
    n = ctypes.c_long()
    check(S.St7GetTotal(uid, S.tyBEAM, ctypes.byref(n)))
    return n.value
