"""KMP (Knuth-Morris-Pratt) / MP 法。

失敗関数(border)で線形時間マッチングする。

使用例:
    kmp_search("ababcabc", "abc")  # [2, 5]
    mp_table("ababaca")            # [0, 0, 1, 2, 3, 0, 1]
"""


def mp_table(S):
    """S の MP 失敗関数。A[i] = S[:i] の最長 border。

    A[0] = 0, len(A) = |S|+1。
    計算量: O(N)
    """
    n = len(S)
    A = [0] * (n + 1)
    j = 0
    for i in range(1, n):
        while j and S[i] != S[j]:
            j = A[j]
        if S[i] == S[j]:
            j += 1
        A[i + 1] = j
    return A


def kmp_search(text, pattern):
    """pattern の開始位置をすべて返す(重なりあり)。

    計算量: O(|text| + |pattern|)
    """
    if not pattern:
        return list(range(len(text) + 1))
    A = mp_table(pattern)
    res = []
    j = 0
    m = len(pattern)
    for i, c in enumerate(text):
        while j and c != pattern[j]:
            j = A[j]
        if c == pattern[j]:
            j += 1
        if j == m:
            res.append(i - m + 1)
            j = A[j]
    return res
