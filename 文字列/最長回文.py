"""Manacher 法で各位置を中心とする最長回文長を求める。

R[2*i]     = S[i] を中心とする奇数長の最大回文長
R[2*i+1]   = S[i:i+2] を中心とする偶数長の最大回文長

使用例:
    manacher('abaa')  # [1, 0, 3, 0, 1, 2, 1]
"""


def manacher(S):
    """各中心での最長回文長の配列 R を返す。

    S: 文字列または列
    計算量: O(N)
    """
    C = []
    for a in S:
        C.append(a)
        C.append(0)
    if C:
        C.pop()
    L = len(C)
    R = [0] * L
    i = j = 0
    while i < L:
        while j <= i < L - j and C[i - j] == C[i + j]:
            j += 1
        R[i] = j
        k = 1
        while j - R[i - k] > k <= i < L - k:
            R[i + k] = R[i - k]
            k += 1
        i += k
        j -= k
    for i in range(L):
        if i & 1 == R[i] & 1:
            R[i] -= 1
    return R
