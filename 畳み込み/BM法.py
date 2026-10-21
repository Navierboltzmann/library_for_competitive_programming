"""
Bostan-Mori 法による線形漸化式の第 N 項 (愚直な多項式乗算版)

線形漸化式
    A_k = -C[1] * A_{k-1} - C[2] * A_{k-2} - ... - C[d] * A_{k-d}
の第 N 項 (1-index) を O(d^2 log N) で求める。
C[0] は 1 を想定 (分母多項式 Q(x) = C[0] + C[1] x + ... + C[d] x^d)。
mod は素数でなくてもよい (逆元を使わない)。

使用例 (フィボナッチ数 F_1 = F_2 = 1, F_k = F_{k-1} + F_{k-2}):
    bostan_mori([1, -1, -1], [1, 1], 10, 10**9 + 7)  # -> 55
"""


def bostan_mori(C, A, N, mod):
    """
    線形漸化式の第 N 項 (1-index) を mod で返す

    C: 漸化式の係数 [1, c1, c2, ..., cd] (A_k = -c1*A_{k-1} - ... - cd*A_{k-d})
    A: 初項 [a1, a2, ..., an] (n >= d が必要)
    N: 求める項の番号 (1-index, N >= 1)
    mod: 法
    戻り値: A_N % mod
    計算量: O(n^2 log N) (n = max(len(A), len(C)))
    """
    N -= 1
    n = len(A)
    # P = A * C mod x^n
    P = [0] * n
    for i in range(n):
        ai = A[i]
        for j in range(min(n - i, len(C))):
            P[i + j] = (P[i + j] + ai * C[j]) % mod
    Q = [c % mod for c in C]
    if len(Q) < n:
        Q += [0] * (n - len(Q))
    m = len(Q)

    while N:
        # Qm = Q(-x)
        Qm = [(mod - q) % mod if j & 1 else q for j, q in enumerate(Q)]
        # P <- (P(x) Q(-x)) の N の偶奇に合う側の係数
        par = N & 1
        lp = len(P)
        U = [0] * (lp + m)
        for i in range(lp):
            pi = P[i]
            if pi:
                for j in range((par ^ i) & 1, m, 2):
                    U[i + j] = (U[i + j] + pi * Qm[j]) % mod
        P = U[par:lp + m - 1:2]
        # Q <- Q(x) Q(-x) の偶数次の係数
        V = [0] * (2 * m)
        for i in range(m):
            qi = Q[i]
            if qi:
                for j in range(i & 1, m, 2):
                    V[i + j] = (V[i + j] + qi * Qm[j]) % mod
        Q = V[0:2 * m:2]
        N >>= 1
    return P[0] % mod


if __name__ == "__main__":
    print(bostan_mori([1, -1, -1], [1, 1], 10, 10**9 + 7))  # 55
