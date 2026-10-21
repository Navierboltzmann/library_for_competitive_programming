"""
LCA (最小共通祖先)・重みなし木 [部分木クエリ + 頂点 Wavelet]

行きがけ順 + Sparse Table で、前計算 O(N log N)、LCA クエリ O(1)。
さらに
    - tree_wavelet(x, k): x の部分木に含まれる深さ k の頂点数 (Wavelet Matrix, O(log N))
    - 頂点重み w を与えた場合: 頂点重みの一点更新・部分木の集約 (Segment Tree, O(log N))
ができる。DFS は非再帰なので setrecursionlimit は不要。

使用例:
    N = 5
    e = [[] for _ in range(N)]
    for a, b in [(0, 1), (0, 2), (1, 3), (1, 4)]:
        e[a].append(b)
        e[b].append(a)
    T = LCA(e, root=0)
    T.lca(3, 4)          # -> 1
    T.distance(3, 2)     # -> 3
    T.tree_wavelet(0, 2) # -> 2 (頂点 3, 4)
    T.tree_wavelet(1, 1) # -> 1 (頂点 1)
"""


class LCA:
    """
    重みなし木の LCA・2頂点間距離・部分木クエリ・部分木内の深さ別頂点数

    属性 (頂点番号で引く):
        IN[v], OUT[v]: オイラーツアーで v に入った / 出た時刻 (0 ～ 2N-1)。
                       v の部分木は時刻の区間 [IN[v], OUT[v]) に対応する。
        DIST[v]: 根から v までの辺数 (depth[v] と同じリスト)
        parents[v]: v の親 (根は -1)
        order: 行きがけ順に並べた頂点列
        pos[v]: order における v の位置
        W: 頂点重みの Segment Tree (w を与えたときのみ)。
           位置 IN[v] に v の重み、それ以外の位置は単位元。
        Wavelet: 行きがけ順に深さを並べた Wavelet Matrix。
                 v の部分木は区間 [pos[v], pos[v] + 部分木サイズ) に対応する。
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
        DIST = self.DIST
        init_wavelet = [DIST[v] for v in self.order]
        self.Wavelet = self.WaveletMatrix(init_wavelet, max(1, max(init_wavelet).bit_length()))
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

    def tree_wavelet(self, x, k):
        """
        頂点 x の部分木に含まれる、深さ (根からの辺数) が k の頂点の個数を返す

        計算量: O(log N)
        """
        l = self.pos[x]
        return self.Wavelet.rank(k, l, l + ((self.OUT[x] - self.IN[x] + 1) >> 1))

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
                    """頂点値の最小値を結合する。"""
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


    class WaveletMatrix:
        """
        Wavelet Matrix (静的な非負整数列に対する区間クエリ)

        値は 0 <= v < 2^log であること。
        """

        def __init__(self, array, log=32):
            """
            構築

            array: 非負整数のリスト
            log: 値のビット数
            計算量: O(n log)
            """
            n = len(array)
            self.n = n
            self.log = log
            self.mat = []
            self.zs = []
            array = list(array)
            for d in range(log - 1, -1, -1):
                BV = self.BitVector(n + 1)
                bit = BV.bit
                ls = []
                rs = []
                for i, val in enumerate(array):
                    if val >> d & 1:
                        rs.append(val)
                        bit[i] = 1
                    else:
                        ls.append(val)
                BV.build()
                self.mat.append(BV)
                self.zs.append(len(ls))
                ls.extend(rs)
                array = ls
            self.cnts = [BV.cnt for BV in self.mat]

        def access(self, i):
            """
            i 番目の値を返す

            計算量: O(log)
            """
            mat = self.mat
            zs = self.zs
            res = 0
            for d in range(self.log):
                res <<= 1
                BV = mat[d]
                if BV.bit[i]:
                    res |= 1
                    i = BV.cnt[i] + zs[d]
                else:
                    i -= BV.cnt[i]
            return res

        def rank(self, val, l, r):
            """
            区間 [l, r) に含まれる val の個数を返す

            計算量: O(log)
            """
            log = self.log
            if val < 0 or val >> log:
                return 0
            cnts = self.cnts
            zs = self.zs
            for d in range(log):
                cnt = cnts[d]
                if val >> (log - d - 1) & 1:
                    z = zs[d]
                    l = cnt[l] + z
                    r = cnt[r] + z
                else:
                    l -= cnt[l]
                    r -= cnt[r]
            return r - l

        def quantile(self, l, r, k):
            """
            区間 [l, r) で k 番目 (1-index) に大きい値を返す (1 <= k <= r - l)

            計算量: O(log)
            """
            cnts = self.cnts
            zs = self.zs
            res = 0
            for d in range(self.log):
                res <<= 1
                cnt = cnts[d]
                cl = cnt[l]
                cr = cnt[r]
                if cr - cl >= k:
                    z = zs[d]
                    l = cl + z
                    r = cr + z
                    res |= 1
                else:
                    l -= cl
                    r -= cr
                    k -= cr - cl
            return res

        def kth_smallest(self, l, r, k):
            """
            区間 [l, r) で k 番目 (0-index) に小さい値を返す (0 <= k < r - l)

            計算量: O(log)
            """
            return self.quantile(l, r, r - l - k)

        class BitVector:
            """
            rank (先頭からの 1 の個数) を O(1) で返すビット列

            bit[i] にビットを持ち、build 後は cnt[i] = bit[0:i] の 1 の個数。
            """

            def __init__(self, size):
                """長さ size の 0 埋めビット列を作る"""
                self.size = size
                self.bit = [0] * size
                self.cnt = [0] * (size + 1)

            def set(self, i):
                """i ビット目を 1 にする (build 前に呼ぶ)"""
                self.bit[i] = 1

            def build(self):
                """
                rank 用の累積和を構築する

                計算量: O(size)
                """
                cnt = self.cnt
                s = 0
                for i, b in enumerate(self.bit):
                    s += b
                    cnt[i + 1] = s

            def __getitem__(self, i):
                """i ビット目を返す"""
                return self.bit[i]

            def popcount(self, x):
                """32bit 整数 x の立っているビット数を返す"""
                x = x - ((x >> 1) & 0x55555555)
                x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
                x = (x + (x >> 4)) & 0x0f0f0f0f
                x = x + (x >> 8)
                x = x + (x >> 16)
                return x & 0x0000007f

            def rank1(self, r):
                """[0, r) の 1 の個数を返す"""
                return self.cnt[r]

            def rank0(self, r):
                """[0, r) の 0 の個数を返す"""
                return r - self.cnt[r]


if __name__ == "__main__":
    N = 5
    e = [[] for _ in range(N)]
    for a, b in [(0, 1), (0, 2), (1, 3), (1, 4)]:
        e[a].append(b)
        e[b].append(a)
    T = LCA(e, [1, 2, 3, 4, 5], lambda x, y: x + y, 0, root=0)
    print(T.lca(3, 4), T.distance(3, 2), T.tree_query(1))  # 1 3 11
    print(T.tree_wavelet(0, 2), T.tree_wavelet(1, 1), T.tree_wavelet(2, 2))  # 2 1 0
    T.weight_update(3, 10)
    print(T.tree_query(1))  # 17
