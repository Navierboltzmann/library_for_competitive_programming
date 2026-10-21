"""最小全域木 (Kruskal)。

辺を重み順に見て、閉路を作らないものだけ採用する。

使用例:
    edges = [(1, 0, 1), (2, 1, 2), (2, 0, 2)]
    mst_kruskal(3, edges)  # (3, [(0, 1, 1), (1, 2, 2)])
"""


class _DSU:
    """内部用 Union-Find。"""

    def __init__(self, n):
        self.p = [-1] * n

    def leader(self, a):
        while self.p[a] >= 0:
            if self.p[self.p[a]] >= 0:
                self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def merge(self, a, b):
        a, b = self.leader(a), self.leader(b)
        if a == b:
            return False
        if self.p[a] > self.p[b]:
            a, b = b, a
        self.p[a] += self.p[b]
        self.p[b] = a
        return True


def mst_kruskal(n, edges):
    """頂点数 n、辺 (cost, u, v) の MST。

    戻り値 (総重み, 使った辺のリスト)。非連結なら全域にならない。
    計算量: O(M log M)
    """
    dsu = _DSU(n)
    used = []
    total = 0
    for cost, u, v in sorted(edges):
        if dsu.merge(u, v):
            used.append((u, v, cost))
            total += cost
            if len(used) == n - 1:
                break
    return total, used
