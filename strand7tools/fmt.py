"""Number formatting shared by the reports and worked calculations."""

import math


def eng(x, decimals=3):
    """Engineering notation: mantissa 1..1000, exponent a multiple of 3 (e.g. 85.056e+06)."""
    if x == 0:
        return f"{0:.{decimals}f}e+00"
    exp = int(math.floor(math.log10(abs(x)) / 3) * 3)
    return f"{x / 10 ** exp:.{decimals}f}e{exp:+03d}"


def kN(x):
    """N or N*m -> kN or kNm, 3 decimals."""
    return f"{x / 1e3:.3f}"


def fmt_util(u):
    if u is None:
        return "-"
    if math.isinf(u):
        return "inf"
    return f"{u:.3f}"
