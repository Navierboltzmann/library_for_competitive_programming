"""1-index の座標圧縮。

重複を除いて昇順に並べ、元の各値を 1, 2, 3, ... に写す。

使用例:
    zip_coordinate([10, 3, 10, 7])  # [3, 1, 3, 2]
"""


def zip_coordinate(x):
    """列 x を 1-index で座標圧縮した列を返す。

    x: 圧縮対象の列
    計算量: O(N log N)
    """
    tx = sorted(set(x))
    d = {v: i + 1 for i, v in enumerate(tx)}
    return [d[v] for v in x]
