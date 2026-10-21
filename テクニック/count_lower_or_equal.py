"""ソート済み列に対する「a 以下の個数」。

使用例:
    cle(3, [1, 2, 2, 4, 5])  # 3
"""
from bisect import bisect_right


def cle(a, D):
    """ソート済み D のうち a 以下の要素数を返す。

    a: 閾値
    D: 昇順ソート済みの列(空でない想定)
    計算量: O(log N)
    """
    return bisect_right(D, a)
