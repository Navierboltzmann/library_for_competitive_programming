"""ソート済み列に対する「a 以上の最小値」。

使用例:
    min_upper_or_equal(3, [1, 2, 4, 5])  # 4
    min_upper_or_equal(6, [1, 2, 4, 5])  # -1
"""
from bisect import bisect_left


def min_upper_or_equal(a, D):
    """ソート済み D のうち a 以上の最小値を返す。存在しなければ -1。

    a: 閾値
    D: 昇順ソート済みの列(空でない想定)
    計算量: O(log N)
    """
    i = bisect_left(D, a)
    if i == len(D):
        return -1
    return D[i]
