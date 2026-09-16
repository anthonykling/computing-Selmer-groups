import re
from GaussianInt import GaussInt


def parse_gaussian_integer(s):
    s = s.replace(" ", "")

    # a + bi
    match = re.fullmatch(r"([+-]?\d+)([+-]\d*)i", s)

    if match:
        a = int(match.group(1))
        b_str = match.group(2)

        if b_str in ("+", ""):
            b = 1
        elif b_str == "-":
            b = -1
        else:
            b = int(b_str)

        return GaussInt(a, b)

    # ordinary integer
    if re.fullmatch(r"[+-]?\d+", s):
        return GaussInt(int(s), 0)

    raise ValueError(
        "Invalid Gaussian integer. "
        "Examples: 3+2i, -5+i, 7-4i"
    )