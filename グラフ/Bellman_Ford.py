"""
ベルマンフォード法 (負辺ありの単一始点最短路・負閉路検出)

使用例:
    e = [[(1, 4), (2, 1)], [], [(1, -2)]]   # e[v] = [(行き先, コスト), ...]
    bf(0, e)                               # [0, -1, 1]
    s から到達可能な負閉路があれば -1 を返す。到達不能な頂点の距離は INF。
"""

INF = 1 << 60


def bf(s, e):
    """
    始点 s から各頂点への最短距離を返す O(NM)

    s: 始点
    e: 隣接リスト。e[v] = [(行き先, コスト), ...]
    戻り値: 距離のリスト (到達不能は INF)。s から到達可能な負閉路があれば -1
    """
    n = len(e)
    d = [INF] * n
    d[s] = 0
    for _ in range(n):
        update = False
        for j in range(n):
            dj = d[j]
            if dj == INF:
                continue
            for k, w in e[j]:
                nd = dj + w
                if nd < d[k]:
                    d[k] = nd
                    update = True
        if not update:
            return d
    return -1


if __name__ == "__main__":
    e = [[(1, 4), (2, 1)], [], [(1, -2)]]
    print(bf(0, e))  # [0, -1, 1]
    print(bf(0, [[(1, 1)], [(0, -2)]]))  # -1
