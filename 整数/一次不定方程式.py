"""一次不定方程式 ax + by = d。

使用例:
    axbyd(3, 5, 1)  # (2, -1, 5, -3) など
"""
from math import gcd


def modinv(a, b):
    """a の mod b における逆元。gcd(a, b)=1 を仮定。"""
    p = b
    x, y, u, v = 1, 0, 0, 1
    while b:
        k = a // b
        x -= k * u
        y -= k * v
        x, u = u, x
        y, v = v, y
        a, b = b, a % b
    x %= p
    if x < 0:
        x += p
    return x


def solve_linear_congruence(A, B, mod):
    """Ax ≡ B (mod M) の解を 1 つ返す。解なしなら -1。"""
    a, b, m = A, B, mod
    d = gcd(a, gcd(b, m))
    while d > 1:
        a //= d
        b //= d
        m //= d
        d = gcd(a, gcd(b, m))
    if gcd(a, m) != 1:
        return -1
    return (modinv(a, m) * b) % m


def axbyd(a, b, d):
    """ax + by = d の一般解 (x0, y0, dx, dy) を返す。解なしなら -1。

    一般解は (x0 + t*dx, y0 + t*dy) (t は整数)。
    計算量: O(log max(a, b))
    """
    g = gcd(a, b)
    if d % g != 0:
        return -1
    a //= g
    b //= g
    d //= g
    if a < 0:
        a *= -1
        b *= -1
        d *= -1
    if b < 0:
        mb = -b
        _x0 = solve_linear_congruence(a, 1, mb)
        if _x0 == -1:
            return -1
        _y0 = (a * _x0 - 1) // mb
        return (_x0 * d, _y0 * d, b, -a)
    _x0 = solve_linear_congruence(a, 1, b)
    if _x0 == -1:
        return -1
    _y0 = (1 - a * _x0) // b
    return (_x0 * d, _y0 * d, b, -a)
