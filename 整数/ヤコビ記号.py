"""Legendre / Jacobi 記号。

(a/p) = 0, ±1。平方剰余の相互法則を二進 GCD 風に計算する。

使用例:
    legendre(2, 7)   # 1    3^2 ≡ 2
    jacobi(2, 15)    # 1    ただし 15 は合成数なので平方剰余とは限らない
"""


def jacobi(a, n):
    """Jacobi 記号 (a/n)。n は正の奇数。

    n が奇素数なら Legendre 記号に一致する。
    計算量: O(log a log n)
    """
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer")
    a %= n
    ans = 1
    while a:
        while a % 2 == 0:
            a //= 2
            n8 = n & 7
            if n8 == 3 or n8 == 5:
                ans = -ans
        a, n = n, a
        if a & 3 == 3 and n & 3 == 3:
            ans = -ans
        a %= n
    return ans if n == 1 else 0


def legendre(a, p):
    """Legendre 記号 (a/p)。p は奇素数。

    1: 平方剰余, -1: 平方非剰余, 0: p | a。
    """
    return jacobi(a, p)


def is_quadratic_residue(a, p):
    """F_p で a が平方剰余か。p=2 も扱う。"""
    a %= p
    if a == 0 or p == 2:
        return True
    return pow(a, (p - 1) // 2, p) == 1
