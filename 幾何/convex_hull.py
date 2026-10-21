"""2次元点集合の凸包(Andrew's monotone chain)。

使用例:
    convex_hull([(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)])
    # [(0, 0), (1, 0), (1, 1), (0, 1)]
"""


def convex_hull(ps):
    """点集合 ps の凸包を反時計回りで返す(共線は落とす)。

    ps: (x, y) のリスト。破壊的に sort する。
    計算量: O(N log N)
    """
    ps.sort()

    def cross(o, a, b):
        """oa × ob の z 成分。"""
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in ps:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(ps):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]
