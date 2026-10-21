"""半環上の行列累乗。

f が「加法」、g が「乗法」。単位元は問題に合わせて C の対角を書き換える。
デフォルト例は (f, g) = (xor, and) で、単位元は (1<<30)-1。

使用例:
    A = [[1, 2], [3, 4]]
    matrix_pow(A, 3, f, g)
"""
from functools import reduce


def f(x, y):
    """半環の加法(デフォルトは XOR)。"""
    return x ^ y


def g(x, y):
    """半環の乗法(デフォルトは AND)。"""
    return x & y


def matrix_mul(A, B, f, g, mod=None):
    """半環上の行列積。

    f: 加法
    g: 乗法
    計算量: O(n m k)
    """
    nA = len(A)
    mA = len(A[0])
    mB = len(B[0])
    tmp = [[0] * mB for _ in range(nA)]
    if mod is None:
        for i in range(nA):
            for j in range(mB):
                tmp[i][j] = reduce(f, (g(A[i][k], B[k][j]) for k in range(mA)))
        return tmp
    for i in range(nA):
        for j in range(mB):
            tmp[i][j] = reduce(f, (g(A[i][k], B[k][j]) % mod for k in range(mA))) % mod
    return tmp


def matrix_pow(A, n, f, g, mod=None):
    """半環上の行列累乗。対角の単位元は (1<<30)-1 (AND の単位)。

    計算量: O(N^3 log n)
    """
    N = len(A)
    C = [[0] * N for _ in range(N)]
    for i in range(N):
        C[i][i] = (1 << 30) - 1
    B = [row[:] for row in A]
    while n:
        if n & 1:
            C = matrix_mul(C, B, f, g, mod)
        B = matrix_mul(B, B, f, g, mod)
        n >>= 1
    return C
