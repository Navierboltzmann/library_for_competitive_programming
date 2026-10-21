"""子リスト付き DSU。

各集合のメンバを merge 時に連結して保持する。members(a) が O(1) で取れる。

使用例:
    dsu = DSU(4)
    dsu.merge(0, 1)
    dsu.members(0)  # [0, 1]
"""


class DSU:
    """メンバリスト付き Union-Find。"""

    def __init__(self, n):
        """n 頂点の DSU を作る。

        n: 頂点数
        """
        self._n = n
        self.parent_or_size = [-1] * n
        self.member = [[i] for i in range(n)]

    def merge(self, a, b):
        """a と b を併合し、新しい leader を返す。

        a, b: 0-index
        計算量: ならし O(log N) (リスト連結込み)
        """
        x = self.leader(a)
        y = self.leader(b)
        if x == y:
            return x
        parent = self.parent_or_size
        if -parent[x] < -parent[y]:
            x, y = y, x
        parent[x] += parent[y]
        self.member[x].extend(self.member[y])
        parent[y] = x
        return x

    def members(self, a):
        """a の属する集合のメンバリストを返す。

        計算量: ほぼ O(1)
        """
        return self.member[self.leader(a)]

    def same(self, a, b):
        """a と b が同じ集合か。"""
        return self.leader(a) == self.leader(b)

    def leader(self, a):
        """a の代表元を返す(経路圧縮、反復)。"""
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
        """a の属する集合のサイズ。"""
        return -self.parent_or_size[self.leader(a)]

    def groups(self):
        """各集合の頂点リストを返す。"""
        n = self._n
        leader_buf = [self.leader(i) for i in range(n)]
        result = [[] for _ in range(n)]
        for i in range(n):
            result[leader_buf[i]].append(i)
        return [r for r in result if r]
