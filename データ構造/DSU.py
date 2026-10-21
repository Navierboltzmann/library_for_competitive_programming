"""Disjoint Set Union (Union-Find)。

ACL 互換の API。経路圧縮 + union by size。

使用例:
    dsu = DSU(4)
    dsu.merge(0, 1)
    dsu.same(0, 1)  # True
    dsu.size(0)     # 2
"""


class DSU:
    """素集合データ構造。"""

    def __init__(self, n):
        """n 頂点の DSU を作る。

        n: 頂点数
        """
        self._n = n
        self.parent_or_size = [-1] * n

    def merge(self, a, b):
        """a と b を併合し、新しい leader を返す。

        a, b: 0-index
        計算量: ほぼ O(1)
        """
        x = self.leader(a)
        y = self.leader(b)
        if x == y:
            return x
        parent = self.parent_or_size
        if -parent[x] < -parent[y]:
            x, y = y, x
        parent[x] += parent[y]
        parent[y] = x
        return x

    def same(self, a, b):
        """a と b が同じ集合か。

        計算量: ほぼ O(1)
        """
        return self.leader(a) == self.leader(b)

    def leader(self, a):
        """a の代表元を返す(経路圧縮)。

        計算量: ほぼ O(1)
        """
        parent = self.parent_or_size
        root = a
        while parent[root] >= 0:
            root = parent[root]
        while a != root:
            nxt = parent[a]
            parent[a] = root
            a = nxt
        return root

    def size(self, a):
        """a の属する集合のサイズ。

        計算量: ほぼ O(1)
        """
        return -self.parent_or_size[self.leader(a)]

    def groups(self):
        """各集合の頂点リストを返す。

        計算量: O(N)
        """
        n = self._n
        leader_buf = [self.leader(i) for i in range(n)]
        result = [[] for _ in range(n)]
        for i in range(n):
            result[leader_buf[i]].append(i)
        return [r for r in result if r]
