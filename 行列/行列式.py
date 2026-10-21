"""素数法での行列式。

使用例:
    determinant([[1, 2], [3, 4]], 998244353)  # -2
"""


def determinant(A, mod, replace=False):
    """det(A) mod を返す。

    A: n×n 行列
    mod: 素数
    replace: True なら A を破壊的に掃き出す
    計算量: O(n^3)
    """
    if not replace:
        A = [a[:] for a in A]
    n = len(A)
    res = 1
    for i, a_i in enumerate(A):
        if a_i[i] == 0:
            for j in range(i + 1, n):
                if A[j][i]:
                    break
            else:
                return 0
            A[i], A[j] = A[j], A[i]
            a_i = A[i]
            res = -res
        inv = pow(a_i[i], mod - 2, mod)
        for j in range(i + 1, n):
            a_j = A[j]
            t = a_j[i] * inv % mod
            if t:
                for k in range(i + 1, n):
                    a_j[k] = (a_j[k] - t * a_i[k]) % mod
    for i in range(n):
        res = res * A[i][i] % mod
    return res
