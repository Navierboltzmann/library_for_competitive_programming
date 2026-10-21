"""Z-algorithm。

Z[i] = S と S[i:] の最長共通接頭辞の長さ。Z[0] = |S|。

使用例:
    z_algorithm("aaabaab")  # [7, 2, 1, 0, 2, 1, 0]
    count_occurrences("aaabaab", "aa")  # 2
"""


def z_algorithm(S):
    """文字列 S の Z 配列を返す。

    計算量: O(N)
    """
    n = len(S)
    Z = [0] * n
    if n == 0:
        return Z
    Z[0] = n
    l = r = 0
    for i in range(1, n):
        if i < r:
            Z[i] = min(r - i, Z[i - l])
        while i + Z[i] < n and S[Z[i]] == S[i + Z[i]]:
            Z[i] += 1
        if i + Z[i] > r:
            l, r = i, i + Z[i]
    return Z


def count_occurrences(text, pattern):
    """pattern の text 内出現回数(重なりあり)。"""
    if not pattern:
        return len(text) + 1
    Z = z_algorithm(list(pattern) + [None] + list(text))
    m = len(pattern)
    return sum(1 for v in Z[m + 1:] if v == m)
