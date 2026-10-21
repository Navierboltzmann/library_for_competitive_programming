"""トポロジカルソート (Kahn)。

有向グラフの入次数を潰して線形順序を得る。閉路があれば空リスト。

使用例:
    g = [[1], [2], []]
    topological_sort(g)  # [0, 1, 2]
"""
from collections import deque


def topological_sort(graph):
    """隣接リスト graph のトポロジカル順。閉路なら []。

    graph[v] = v から出る行き先のリスト。0-index。
    同じ深さは番号の小さい順。
    計算量: O(N + M)
    """
    n = len(graph)
    indeg = [0] * n
    for v in range(n):
        for to in graph[v]:
            indeg[to] += 1
    q = deque(v for v in range(n) if indeg[v] == 0)
    order = []
    while q:
        v = q.popleft()
        order.append(v)
        for to in graph[v]:
            indeg[to] -= 1
            if indeg[to] == 0:
                q.append(to)
    if len(order) != n:
        return []
    return order


def is_dag(graph):
    """閉路が無いか。"""
    return len(topological_sort(graph)) == len(graph)
