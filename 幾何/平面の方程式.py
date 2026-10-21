"""点から平面・直線の方程式を求める。

使用例:
    plane((0, 0, 0), (1, 0, 0), (0, 1, 0))  # (0, 0, 1, 0)
    line((0, 0), (2, 2))  # (1, -1, 0)
"""
from math import gcd


def plane(A, B, C):
    """3 点 A, B, C を通る平面 ax+by+cz+d=0 の係数を返す。

    A, B, C: (x, y, z)
    戻り値: (a, b, c, d)
    """
    ax, ay, az = A
    bx, by, bz = B
    cx, cy, cz = C
    a = (by - ay) * (cz - az) - (cy - ay) * (bz - az)
    b = (bz - az) * (cx - ax) - (cz - az) * (bx - ax)
    c = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay)
    d = -(a * ax + b * ay + c * az)
    return (a, b, c, d)


def line(A, B):
    """2 点 A, B を通る直線 ax+by+c=0 の既約係数を返す。

    A, B: (x, y)
    戻り値: (a, b, c)。a > 0、a=0 なら b > 0 に正規化。
    """
    x1, y1 = A
    x2, y2 = B
    _a = y2 - y1
    _b = -(x2 - x1)
    _c = -x1 * y2 + x2 * y1
    g = gcd(gcd(_a, _b), _c)
    if g:
        _a //= g
        _b //= g
        _c //= g
    if _a < 0 or (_a == 0 and _b < 0):
        return (-_a, -_b, -_c)
    return (_a, _b, _c)
