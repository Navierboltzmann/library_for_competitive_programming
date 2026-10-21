"""Union-Find (parents 配列版)。

使用例:
    uf = UnionFind(4)
    uf.union(0, 1)
    uf.same(0, 1)  # True
    uf.group_count()  # 3
"""
from collections import defaultdict


class UnionFind:
    """素集合データ構造。"""

    def __init__(self, n):
        """n 頂点の Union-Find を作る。

        n: 頂点数
        """
        self.n = n
        self.parents = [-1] * n

    def find(self, x):
        """x の代表元を返す(経路圧縮、反復)。

        計算量: ほぼ O(1)
        """
        parents = self.parents
        root = x
        while parents[root] >= 0:
            root = parents[root]
        while x != root:
            nxt = parents[x]
            parents[x] = root
            x = nxt
        return root

    def union(self, x, y):
        """x と y を併合する。

        計算量: ほぼ O(1)
        """
        x = self.find(x)
        y = self.find(y)
        if x == y:
            return
        parents = self.parents
        if parents[x] > parents[y]:
            x, y = y, x
        parents[x] += parents[y]
        parents[y] = x

    def size(self, x):
        """x の属する集合のサイズ。"""
        return -self.parents[self.find(x)]

    def same(self, x, y):
        """x と y が同じ集合か。"""
        return self.find(x) == self.find(y)

    def members(self, x):
        """x の属する集合のメンバを列挙する。

        計算量: O(N)
        """
        root = self.find(x)
        return [i for i in range(self.n) if self.find(i) == root]

    def roots(self):
        """全代表元のリスト。"""
        return [i for i, x in enumerate(self.parents) if x < 0]

    def group_count(self):
        """集合の個数。"""
        return len(self.roots())

    def all_group_members(self):
        """代表元 → メンバリスト の dict。

        計算量: O(N)
        """
        group_members = defaultdict(list)
        for member in range(self.n):
            group_members[self.find(member)].append(member)
        return group_members
