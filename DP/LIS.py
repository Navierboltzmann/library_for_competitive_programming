"""最長増加部分列 (LIS) の長さ。

seq の狭義増加部分列のうち最長のものの長さを O(N log N) で求める。

使用例:
    length_LIS([1, 3, 2, 4])  # 3  ([1, 2, 4] など)
"""
from bisect import bisect_left


def length_LIS(seq):
    """seq の最長狭義増加部分列の長さを返す。

    seq: 数列
    計算量: O(N log N)
    """
    if not seq:
        return 0
    lis = [seq[0]]
    for x in seq:
        if x > lis[-1]:
            lis.append(x)
        else:
            lis[bisect_left(lis, x)] = x
    return len(lis)
