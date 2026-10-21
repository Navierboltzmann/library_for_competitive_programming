"""
LCA (最小共通祖先)・重み付き木

行きがけ順 + Sparse Table で、前計算 O(N log N)、LCA クエリ O(1)。
DFS は非再帰なので setrecursionlimit は不要。
LCA の判定は辺数の深さで行い、辺の重みは距離計算にだけ使うので、重みは 0 や負でもよい。

使用例:
    N = 5
    e = [[] for _ in range(N)]
    for a, b, c in [(0, 1, 3), (0, 2, 1), (1, 3, 4), (1, 4, 2)]:
        e[a].append((b, c))
        e[b].append((a, c))
    T = LCA_wp(e, root=0)
    T.lca(3, 4)               # -> 1
    T.distance_between(3, 2)  # -> 8
    T.distance_from_root(4)   # -> 5
"""


class LCA_wp:
    """
    重み付き木の LCA・2頂点間距離

    属性 (頂点番号で引く):
        IN[v], OUT[v]: オイラーツアーで v に入った / 出た時刻 (0 ～ 2N-1)。
                       v の部分木は時刻の区間 [IN[v], OUT[v]) に対応する。
        DIST[v]: 根から v までの重み付き距離
        depth[v]: 根から v までの辺数
        parents[v]: v の親 (根は -1)
        order: 行きがけ順に並べた頂点列
        pos[v]: order における v の位置
    """

    def __init__(self, e, root=0):
        """
        前計算

        e: 隣接リスト。e[v] = [(隣接頂点, 辺の重み), ...]
        root: 根
        計算量: O(N log N)
        """
        N = len(e)
        self.N = N
        self.e = e
        self.root = root
        self.IN = [-1] * N
        self.OUT = [-1] * N
        self.DIST = [-1] * N
        self.depth = [-1] * N
        self.parents = [-1] * N
        self.order = []
        self.pos = [-1] * N
        self.time = 0
        self.dfs(root)
        self._build_table()

    def dfs(self, x):
        """
        x を根として非再帰 DFS し、IN, OUT, DIST, depth, parents, order を埋める
        (子を訪れる順番は再帰 DFS と同じ)

        x: 根
        計算量: O(N)
        """
        e = self.e
        IN = self.IN
        OUT = self.OUT
        DIST = self.DIST
        depth = self.depth
        par = self.parents
        order = self.order
        t = self.time
        DIST[x] = 0
        depth[x] = 0
        stack = [x]
        while stack:
            v = stack.pop()
            if v >= 0:
                IN[v] = t
                t += 1
                order.append(v)
                stack.append(~v)
                dv = depth[v] + 1
                wv = DIST[v]
                ev = e[v]
                for i in range(len(ev) - 1, -1, -1):
                    u, c = ev[i]
                    if depth[u] == -1:
                        depth[u] = dv
                        DIST[u] = wv + c
                        par[u] = v
                        stack.append(u)
            else:
                OUT[~v] = t
                t += 1
        self.time = t

    def _build_table(self):
        """
        LCA 用の Sparse Table を構築する

        行きがけ順の i 番目の頂点の親 p を depth[p] * N + p と整数に詰めて並べ、区間最小を取る。
        計算量: O(N log N)
        """
        N = self.N
        order = self.order
        pos = self.pos
        par = self.parents
        depth = self.depth
        for i, v in enumerate(order):
            pos[v] = i
        m = len(order)
        base = [0] * m
        for i in range(1, m):
            p = par[order[i]]
            base[i] = depth[p] * N + p
        table = [base]
        k = 1
        while 2 * k <= m:
            prev = table[-1]
            table.append([a if a < b else b for a, b in zip(prev, prev[k:])])
            k <<= 1
        self.table = table

    def lca(self, l, r):
        """
        頂点 l, r の最小共通祖先を返す

        計算量: O(1)
        """
        if l == r:
            return l
        pos = self.pos
        a = pos[l]
        b = pos[r]
        if a > b:
            a, b = b, a
        k = (b - a).bit_length() - 1
        t = self.table[k]
        x = t[a + 1]
        y = t[b - (1 << k) + 1]
        return (x if x < y else y) % self.N

    def distance_between(self, l, r):
        """
        頂点 l, r 間の距離 (重みの和) を返す

        計算量: O(1)
        """
        DIST = self.DIST
        return DIST[l] + DIST[r] - 2 * DIST[self.lca(l, r)]

    distance = distance_between

    def distance_from_root(self, x):
        """根から頂点 x までの距離を返す"""
        return self.DIST[x]


if __name__ == "__main__":
    N = 5
    e = [[] for _ in range(N)]
    for a, b, c in [(0, 1, 3), (0, 2, 1), (1, 3, 4), (1, 4, 2)]:
        e[a].append((b, c))
        e[b].append((a, c))
    T = LCA_wp(e, root=0)
    print(T.lca(3, 4), T.distance_between(3, 2), T.distance_from_root(4))  # 1 8 5
