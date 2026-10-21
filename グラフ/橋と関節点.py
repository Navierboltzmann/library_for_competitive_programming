"""橋と関節点 (lowlink)。

無向グラフの橋(除くと非連結になる辺)と関節点を O(N+M) で列挙する。

使用例:
    g = [[1], [0, 2], [1]]
    bridges(g)                # [(0, 1), (1, 2)]
    articulation_points(g)    # [1]
"""
import sys

sys.setrecursionlimit(10**6)


def _lowlink(graph):
    """ord, low, parent を返す。根の parent は -2。"""
    n = len(graph)
    ord_ = [-1] * n
    low = [0] * n
    parent = [-1] * n
    k = 0

    def dfs(v, p):
        nonlocal k
        ord_[v] = low[v] = k
        k += 1
        for to in graph[v]:
            if to == p:
                continue
            if ord_[to] == -1:
                parent[to] = v
                dfs(to, v)
                if low[v] > low[to]:
                    low[v] = low[to]
            elif low[v] > ord_[to]:
                low[v] = ord_[to]

    for s in range(n):
        if ord_[s] == -1:
            parent[s] = -2
            dfs(s, -1)
    return ord_, low, parent


def bridges(graph):
    """橋を (min, max) のリストで返す。

    計算量: O(N + M)
    """
    ord_, low, parent = _lowlink(graph)
    res = []
    for v, p in enumerate(parent):
        if p >= 0 and low[v] > ord_[p]:
            res.append((v, p) if v < p else (p, v))
    res.sort()
    return res


def articulation_points(graph):
    """関節点(切断点)の頂点番号を昇順で返す。

    計算量: O(N + M)
    """
    n = len(graph)
    ord_, low, parent = _lowlink(graph)
    art = set()
    children = [0] * n
    for v, p in enumerate(parent):
        if p >= 0:
            children[p] += 1
            if parent[p] >= 0 and low[v] >= ord_[p]:
                art.add(p)
    for s in range(n):
        if parent[s] == -2 and children[s] >= 2:
            art.add(s)
    return sorted(art)
