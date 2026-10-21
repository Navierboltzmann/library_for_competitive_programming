"""一次合同方程式 Ax ≡ B (mod M)。

使用例:
    solve_linear_congruence(3, 4, 5)  # 3  (3*3=9 ≡ 4)
"""
from math import gcd


def modinv(a, b):
    """a の mod b における逆元を返す。gcd(a, b)=1 を仮定。

    計算量: O(log b)
    """
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
    """Ax ≡ B (mod M) の解を 1 つ返す。解なしなら -1。

    A, B, mod: 整数
    計算量: O(log M)
    """
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
