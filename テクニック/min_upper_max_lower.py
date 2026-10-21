"""ソート済み列に対する下限・上限の index。

使用例:
    D = [1, 3, 3, 7]
    min_upper_or_equal(3, D)  # 1  (最初の 3 以上)
    max_lower_or_equal(3, D)  # 2  (最後の 3 以下)
"""
from bisect import bisect_left, bisect_right


def min_upper_or_equal(a, D):
    """ソート済み D のうち a 以上の最小 index。無ければ -1。

    a: 閾値
    D: 昇順ソート済みの列(空でない想定)
    計算量: O(log N)
    """
    i = bisect_left(D, a)
    if i == len(D):
        return -1
    return i


def max_lower_or_equal(a, D):
    """ソート済み D のうち a 以下の最大 index。無ければ -1。

    a: 閾値
    D: 昇順ソート済みの列(空でない想定)
    計算量: O(log N)
    """
    i = bisect_right(D, a) - 1
    return i
