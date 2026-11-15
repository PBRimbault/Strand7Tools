"""
Shared test data: a 305x165x40 UB (BS EN 10365 dimensions) in S355. No Strand7 needed.

The section numbers below come from a brute-force numerical integration of the real section
shape, root fillets included, and are NOT calculated by the code under test.
"""

import math

import pytest

from strand7tools.ec3.sections import ISection, Steel
from strand7tools.settings import load_settings

H, B, TW, TF, R = 0.3034, 0.165, 0.006, 0.0102, 0.009
A = 5.133531e-03
IY, IZ = 8.505546e-05, 7.643650e-06
WPL_Y, WPL_Z = 6.232880e-04, 1.417429e-04
WEL_Y, WEL_Z = IY / (H / 2), IZ / (B / 2)
IT = 14.90e-8               # membrane-analogy solution, root fillets included
IW = IZ * (H - TF) ** 2 / 4
HW = H - 2 * TF
FY = 355e6
E, G = 210e9, 81e9
EPS = math.sqrt(235 / 355)


def close(got, want, tol=2e-3):
    """Relative tolerance, with an absolute floor of tol for values near zero."""
    return abs(got - want) <= tol * max(1.0, abs(want))


@pytest.fixture(scope="session")
def sec():
    return ISection("UB 305x165x40", H, B, TW, TF, R, eta=1.0)


@pytest.fixture(scope="session")
def steel():
    return Steel("S355", FY, 470e6, E, G)


@pytest.fixture(scope="session")
def bk():
    return load_settings().buckling
