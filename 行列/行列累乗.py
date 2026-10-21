"""行列積と行列累乗。

使用例:
    A = [[1, 1], [1, 0]]
    matrix_pow(A, 10)  # Fibonacci 行列
"""


def matrix_mul(A, B, mod=None):
    """行列積 A @ B を返す。

    A: n×m
    B: m×k
    mod: 指定時は各要素を mod する
    計算量: O(n m k)
    """
    nA = len(A)
    mA = len(A[0])
    mB = len(B[0])
    tmp = [[0] * mB for _ in range(nA)]
    if mod is None:
        for i in range(nA):
            Ai = A[i]
            ti = tmp[i]
            for k in range(mA):
                aik = Ai[k]
                Bk = B[k]
                for j in range(mB):
                    ti[j] += aik * Bk[j]
        return tmp
    for i in range(nA):
        Ai = A[i]
        ti = tmp[i]
        for k in range(mA):
            aik = Ai[k]
            if aik:
                Bk = B[k]
                for j in range(mB):
                    ti[j] = (ti[j] + aik * Bk[j]) % mod
    return tmp


def matrix_pow(A, n, mod=None):
    """正方行列 A の n 乗を返す。

    n: 非負整数
    計算量: O(N^3 log n)
    """
    N = len(A)
    C = [[0] * N for _ in range(N)]
    for i in range(N):
        C[i][i] = 1
    B = [row[:] for row in A]
    while n:
        if n & 1:
            C = matrix_mul(C, B, mod)
        B = matrix_mul(B, B, mod)
        n >>= 1
    return C
