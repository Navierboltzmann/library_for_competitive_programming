"""
LCA (最小共通祖先)・重みなし木 [部分木クエリ付き]

行きがけ順 + Sparse Table で、前計算 O(N log N)、LCA クエリ O(1)。
さらに頂点重み w を与えると、オイラーツアー上の Segment Tree で
「頂点重みの一点更新」「部分木の頂点重みの集約」を O(log N) で行える。
DFS は非再帰なので setrecursionlimit は不要。

使用例:
    N = 5
    e = [[] for _ in range(N)]
    for a, b in [(0, 1), (0, 2), (1, 3), (1, 4)]:
        e[a].append(b)
        e[b].append(a)
    w = [1, 2, 3, 4, 5]
    T = LCA(e, w, lambda x, y: x + y, 0, root=0)
    T.lca(3, 4)               # -> 1
    T.distance_between(3, 2)  # -> 3
    T.tree_query(1)           # -> 11 (= 2 + 4 + 5)
    T.weight_update(3, 10)
    T.tree_query(1)           # -> 17
"""


class LCA:
    """
    重みなし木の LCA・2頂点間距離・部分木クエリ

    属性 (頂点番号で引く):
        IN[v], OUT[v]: オイラーツアーで v に入った / 出た時刻 (0 ～ 2N-1)。
                       v の部分木は時刻の区間 [IN[v], OUT[v]) に対応する。
        DIST[v]: 根から v までの辺数 (depth[v] と同じリスト)
        parents[v]: v の親 (根は -1)
        order: 行きがけ順に並べた頂点列
        pos[v]: order における v の位置
        W: 頂点重みの Segment Tree (w を与えたときのみ)。
           位置 IN[v] に v の重み、それ以外の位置は単位元。
    """

    def __init__(self, e, w=None, f=None, ide_ele=None, root=0):
        """
        前計算

        e: 隣接リスト。e[v] = [隣接頂点, ...]
        w: 頂点重みのリスト (省略時は部分木クエリなし)
        f: 部分木クエリの二項演算 (結合的であること。f, ide_ele とも省略時は和)
        ide_ele: f の単位元
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
        self.depth = self.DIST
        self.parents = [-1] * N
        self.order = []
        self.pos = [-1] * N
        self.time = 0
        self.dfs(root)
        self._build_table()
        if w is None:
            return
        if f is None and ide_ele is None:
            f = lambda x, y: x + y
            ide_ele = 0
        init_val_w = [ide_ele] * (2 * N)
        IN = self.IN
        for i in range(N):
            if IN[i] >= 0:
                init_val_w[IN[i]] = w[i]
        self.W = self.SegTree(init_val_w, f, ide_ele)

    def weight_update(self, x, y):
        """
        頂点 x の重みを y に更新する

        計算量: O(log N)
        """
        self.W.update(self.IN[x], y)

    def tree_query(self, x):
        """
        頂点 x の部分木に含まれる頂点の重みを f で集約した値を返す

        計算量: O(log N)
        """
        return self.W.query(self.IN[x], self.OUT[x])

    def dfs(self, x):
        """
        x を根として非再帰 DFS し、IN, OUT, DIST, parents, order を埋める
        (子を訪れる順番は再帰 DFS と同じ)

        x: 根
        計算量: O(N)
        """
        e = self.e
        IN = self.IN
        OUT = self.OUT
        DIST = self.DIST
        par = self.parents
        order = self.order
        t = self.time
        DIST[x] = 0
        stack = [x]
        while stack:
            v = stack.pop()
            if v >= 0:
                IN[v] = t
                t += 1
                order.append(v)
                stack.append(~v)
                dv = DIST[v] + 1
                ev = e[v]
                for i in range(len(ev) - 1, -1, -1):
                    u = ev[i]
                    if DIST[u] == -1:
                        DIST[u] = dv
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
        頂点 l, r 間の距離 (辺数) を返す

        計算量: O(1)
        """
        DIST = self.DIST
        return DIST[l] + DIST[r] - 2 * DIST[self.lca(l, r)]

    distance = distance_between

    def distance_from_root(self, x):
        """根から頂点 x までの距離 (深さ) を返す"""
        return self.DIST[x]

    class SegTree:
        """
        Segment Tree (一点更新・区間取得)

        segfunc は結合的であればよい (可換でなくてよい)。
        """

        def __init__(self, init_val, segfunc=None, ide_ele=None):
            """
            初期化

            init_val: 配列の初期値
            segfunc: 二項演算 (省略時は min)
            ide_ele: 単位元 (segfunc, ide_ele とも省略時は (1 << 60, -1))
            計算量: O(n)
            """
            if segfunc is None:
                def segfunc(x, y):
                    """区間クエリの 2 つの結果を結合する。"""
                    return y if y < x else x
                if ide_ele is None:
                    ide_ele = (1 << 60, -1)
            n = len(init_val)
            self.segfunc = segfunc
            self.ide_ele = ide_ele
            num = 1 << (n - 1).bit_length()
            self.num = num
            tree = [ide_ele] * (2 * num)
            tree[num:num + n] = init_val
            for i in range(num - 1, 0, -1):
                tree[i] = segfunc(tree[2 * i], tree[2 * i + 1])
            self.tree = tree

        def update(self, k, x):
            """
            k番目の値をxに更新

            k: index(0-index)
            x: update value
            計算量: O(log n)
            """
            segfunc = self.segfunc
            tree = self.tree
            k += self.num
            tree[k] = x
            k >>= 1
            while k:
                tree[k] = segfunc(tree[2 * k], tree[2 * k + 1])
                k >>= 1

        def query(self, l, r):
            """
            [l, r)のsegfuncしたものを得る

            l: index(0-index)
            r: index(0-index)
            計算量: O(log n)
            """
            segfunc = self.segfunc
            tree = self.tree
            sml = smr = self.ide_ele
            l += self.num
            r += self.num
            while l < r:
                if l & 1:
                    sml = segfunc(sml, tree[l])
                    l += 1
                if r & 1:
                    r -= 1
                    smr = segfunc(tree[r], smr)
                l >>= 1
                r >>= 1
            return segfunc(sml, smr)


if __name__ == "__main__":
    N = 5
    e = [[] for _ in range(N)]
    for a, b in [(0, 1), (0, 2), (1, 3), (1, 4)]:
        e[a].append(b)
        e[b].append(a)
    T = LCA(e, [1, 2, 3, 4, 5], lambda x, y: x + y, 0, root=0)
    print(T.lca(3, 4), T.distance_between(3, 2), T.tree_query(1))  # 1 3 11
    T.weight_update(3, 10)
    print(T.tree_query(1))  # 17
