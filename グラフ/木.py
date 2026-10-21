"""
根付き木ライブラリ

頂点 index, index+1, ..., index+N-1 からなる根付き木を扱う.
深さ, k 個上の祖先, LCA, 距離, 直径, パス, 部分木サイズ, 重心,
木 DP (葉から / 根から / 全方位), オイラーツアー などを提供する.

使用例:
    N = 5
    E = [(0, 1), (0, 2), (1, 3), (1, 4)]
    T = Making_Tree(N, E, 0)           # 辺リストから根 0 の木を作る
    T.lowest_common_ancestor(3, 4)     # -> 1
    T.distance(3, 2)                   # -> 3
    T.diameter()                       # -> (3, (3, 2))

    # 手動で作る場合
    T = Tree(N)
    T.root_set(0)
    for u, v in E:
        T.child_set(u, v)
    T.seal()                           # 以降は木を変更できない

    # 全方位木 DP (各頂点から最も遠い頂点までの距離)
    T.rerooting(max, 0, lambda x, v, w: x + 1, lambda x, v: x)
"""


class Tree:
    __slots__ = ("N", "index", "parent", "__mutable",
                 "root", "children", "depth", "tower", "upper_list", "deg", "des_count", "preorder_number",
                 "euler_vertex", "euler_edge", "in_time", "out_time", "bfs_order")

    def __init__(self, N, index=0):
        """
        N 頂点 (index, index+1, ..., N-1+index) の根付き木を生成する.

        N: 頂点数
        index: 頂点番号の開始値 (0 or 1 など)
        """
        self.N = N
        self.index = index
        self.parent = [-1] * (N + index)
        self.__mutable = True

    def vertex_exist(self, x):
        """頂点 x が存在するかどうかを判定する."""
        return self.index <= x < self.index + self.N

    def __after_seal_check(self, *vertexes):
        """木が確定していて, vertexes の頂点が全て存在するかどうかをチェックする."""
        if self.__mutable:
            return False
        a = self.index
        b = a + self.N
        for v in vertexes:
            if not a <= v < b:
                return False
        return True

    def is_mutable(self):
        """木が確定して [いない] かどうかを返す."""
        return self.__mutable

    # ---------------- 設定パート ----------------
    def root_set(self, root):
        """頂点 root を根に設定する."""
        assert self.vertex_exist(root)
        assert self.__mutable
        self.root = root

    def parent_set(self, x, y):
        """頂点 x の親を y に設定する."""
        assert self.vertex_exist(x)
        assert self.vertex_exist(y)
        assert self.__mutable
        self.parent[x] = y

    def child_set(self, x, y):
        """頂点 x の子の一つに y を設定する."""
        assert self.vertex_exist(x)
        assert self.vertex_exist(y)
        assert self.__mutable
        self.parent[y] = x

    def seal(self):
        """
        木の情報を確定させる (children を構築する). 以降は木を変更できない.

        計算量: O(N)
        """
        assert self.__mutable
        assert hasattr(self, "root")

        a = self.index
        b = a + self.N
        root = self.root
        C = [[] for _ in range(b)]
        p = self.parent
        for i in range(a, b):
            if i != root:
                pi = p[i]
                assert a <= pi < b
                C[pi].append(i)

        self.__mutable = False
        self.children = C

    # ---------------- データを求める ----------------
    def depth_search(self, Mode=True):
        """
        各頂点の深さを求める (根の深さは 0).

        Mode: True なら深さのリストを返す. False なら計算するだけ.
        副作用: depth (深さ), tower (tower[d] = 深さ d の頂点のリスト), bfs_order (BFS 順) を設定.
        計算量: O(N)
        """
        assert self.__after_seal_check()

        if hasattr(self, "depth"):
            if Mode:
                return self.depth
            else:
                return

        C = self.children
        root = self.root
        D = [-1] * (self.index + self.N)
        D[root] = 0
        order = [root]
        for x in order:
            d = D[x] + 1
            for y in C[x]:
                D[y] = d
                order.append(y)

        E = [[] for _ in range(self.N)]
        for x in order:
            E[D[x]].append(x)

        self.depth = D
        self.tower = E
        self.bfs_order = order

        if Mode:
            return D

    def vertex_depth(self, x):
        """頂点 x の深さを求める. O(1) (初回のみ O(N))"""
        assert self.__after_seal_check(x)

        if not hasattr(self, "depth"):
            self.depth_search(Mode=False)
        return self.depth[x]

    def __upper_list(self):
        """ダブリング表 upper_list[k][x] = x の 2^k 個上の頂点 (根より上は根) を作る. O(N log N)"""
        assert self.__after_seal_check()

        if hasattr(self, "upper_list"):
            return

        if not hasattr(self, "depth"):
            self.depth_search(False)

        b = max(self.depth).bit_length()
        idx = self.index

        Y = self.parent[:]
        Y[self.root] = self.root
        X = [Y] if b else []
        for _ in range(1, b):
            Y = [Y[y] for y in Y]
            for x in range(idx):
                Y[x] = -1
            X.append(Y)
        self.upper_list = X

    def __upper(self, x, k):
        """x の k 個上の頂点 (k <= depth[x] を仮定, 前計算済みを仮定)."""
        U = self.upper_list
        i = 0
        while k:
            if k & 1:
                x = U[i][x]
            k >>= 1
            i += 1
        return x

    def upper(self, x, k, over=True):
        """
        頂点 x から見て k 個親の頂点を求める.

        over: (頂点 x の深さ) < k のときに True ならば根を返し, False ならばエラーを吐く.
        計算量: O(log N) (初回のみ O(N log N))
        """
        assert self.__after_seal_check(x)
        assert 0 <= k

        if not hasattr(self, "upper_list"):
            self.__upper_list()

        if self.depth[x] < k:
            if over:
                return self.root
            else:
                raise ValueError

        return self.__upper(x, k)

    def lowest_common_ancestor(self, x, y):
        """
        頂点 x, y の最小共通先祖 (x, y に共通する先祖で最も深いもの) を求める.

        計算量: O(log N) (初回のみ O(N log N))
        """
        assert self.__after_seal_check(x, y)

        if not hasattr(self, "upper_list"):
            self.__upper_list()

        D = self.depth
        U = self.upper_list
        dx = D[x]
        dy = D[y]
        if dx > dy:
            x, y = y, x
            dx, dy = dy, dx

        dd = dy - dx
        i = 0
        while dd:
            if dd & 1:
                y = U[i][y]
            dd >>= 1
            i += 1

        if x == y:
            return x

        for k in range(dx.bit_length() - 1, -1, -1):
            Uk = U[k]
            px = Uk[x]
            py = Uk[y]
            if px != py:
                x = px
                y = py
        return U[0][x]

    def __degree_count(self):
        """各頂点の次数を計算する. O(N)"""
        assert self.__after_seal_check()

        if hasattr(self, "deg"):
            return

        root = self.root
        ch = self.children
        deg = [0] * (self.index + self.N)
        for v in range(self.index, self.index + self.N):
            deg[v] = len(ch[v]) + (v != root)
        self.deg = deg

    def degree(self, v):
        """頂点 v の次数を求める."""
        assert self.__after_seal_check(v)

        if not hasattr(self, "deg"):
            self.__degree_count()
        return self.deg[v]

    def diameter(self):
        """
        木の直径を求める.

        戻り値: (直径の長さ (辺数), (端点1, 端点2))
        計算量: O(N)
        """
        assert self.__after_seal_check()

        root = self.root
        pa = self.parent
        ch = self.children
        size = self.index + self.N

        def bfs(start):
            """start から最遠の頂点とその距離を求める。"""
            X = [-1] * size
            X[start] = 0
            Q = [start]
            for x in Q:
                dx = X[x] + 1
                if x != root:
                    p = pa[x]
                    if X[p] == -1:
                        X[p] = dx
                        Q.append(p)
                for y in ch[x]:
                    if X[y] == -1:
                        X[y] = dx
                        Q.append(y)
            m = max(X)
            return X.index(m), m

        y, _ = bfs(root)
        z, d = bfs(y)
        return d, (y, z)

    def path(self, u, v):
        """
        頂点 u, v 間のパス (u, ..., v の頂点リスト) を求める.

        計算量: O(パス長 + log N)
        """
        assert self.__after_seal_check(u, v)

        w = self.lowest_common_ancestor(u, v)
        pa = self.parent

        X = [u]
        while u != w:
            u = pa[u]
            X.append(u)

        Y = [v]
        while v != w:
            v = pa[v]
            Y.append(v)
        return X + Y[-2::-1]

    def is_parent(self, u, v):
        """u は v の親か?"""
        assert self.__after_seal_check(u, v)
        return v != self.root and u == self.parent[v]

    def is_children(self, u, v):
        """u は v の子か?"""
        assert self.__after_seal_check(u, v)
        return self.is_parent(v, u)

    def is_brother(self, u, v):
        """2つの頂点 u, v は兄弟 (親が同じ) か?"""
        assert self.__after_seal_check(u, v)

        if u == self.root or v == self.root:
            return False
        return self.parent[u] == self.parent[v]

    def is_ancestor(self, u, v):
        """頂点 u は頂点 v の先祖か? (u == v も True). O(log N)"""
        assert self.__after_seal_check(u, v)

        if not hasattr(self, "upper_list"):
            self.__upper_list()

        dd = self.depth[v] - self.depth[u]
        if dd < 0:
            return False
        return u == self.__upper(v, dd)

    def is_descendant(self, u, v):
        """頂点 u は頂点 v の子孫か? (u == v も True). O(log N)"""
        assert self.__after_seal_check(u, v)
        return self.is_ancestor(v, u)

    def direction(self, u, v):
        """頂点 u から頂点 v へ向かうパスで, 頂点 u の次に通る頂点を求める. O(log N)"""
        assert self.__after_seal_check(u, v)
        assert u != v

        if self.is_ancestor(u, v):
            D = self.depth
            return self.__upper(v, D[v] - (D[u] + 1))
        else:
            return self.parent[u]

    def is_leaf(self, v):
        """頂点 v は葉 (子を持たない) か?"""
        return not self.children[v]

    def distance(self, u, v):
        """2頂点 u, v 間の距離 (辺数) を求める. O(log N)"""
        assert self.__after_seal_check(u, v)

        w = self.lowest_common_ancestor(u, v)
        D = self.depth
        return D[u] + D[v] - 2 * D[w]

    def __descendant_count(self):
        """各頂点の子孫の数 (自身を含む) を計算する. O(N)"""
        assert self.__after_seal_check()
        if hasattr(self, "des_count"):
            return

        if not hasattr(self, "depth"):
            self.depth_search(False)

        des = [1] * (self.index + self.N)
        pa = self.parent
        order = self.bfs_order
        for i in range(len(order) - 1, 0, -1):
            x = order[i]
            des[pa[x]] += des[x]
        self.des_count = des

    def descendant_count(self, v):
        """頂点 v の子孫の数 (v 自身を含む) を求める."""
        assert self.__after_seal_check(v)
        self.__descendant_count()
        return self.des_count[v]

    def subtree_size(self, v):
        """頂点 v を根とした部分根付き木のサイズ (頂点数) を求める."""
        return self.descendant_count(v)

    def preorder(self, v):
        """
        頂点 v の行きがけ順 (DFS preorder, 1-indexed) を求める.

        子は children[v] の順に訪問する.
        計算量: O(1) (初回のみ O(N))
        """
        assert self.__after_seal_check(v)
        if hasattr(self, "preorder_number"):
            return self.preorder_number[v]

        ch = self.children
        T = [-1] * (self.N + self.index)
        stack = [self.root]
        p = 1
        while stack:
            x = stack.pop()
            T[x] = p
            p += 1
            c = ch[x]
            if c:
                stack.extend(reversed(c))
        self.preorder_number = T
        return T[v]

    def dfs_yielder(self, order=None):
        """
        DFS における頂点の出入りを yield する.

        以下のような関数を (仮想的に) 実行する.

        def dfs(v):
            yield (v, 1)  # 頂点 v に入る
            for w in self.children[v]:
                dfs(w)
            yield (v, 0)  # 頂点 v を出る

        order (1変数関数): for w in self.children[v] の順番を指定する (昇順) (※ 無い場合は任意, 破壊的)
        """
        assert self.__after_seal_check()

        root = self.root
        ch = self.children
        pa = self.parent

        if order is not None:
            for w in range(self.index, self.index + self.N):
                ch[w].sort(key=order)

        S = [0] * (self.index + self.N)

        yield (root, 1)
        v = root
        while True:
            c = ch[v]
            s = S[v]
            if s == len(c):  # 進めないので頂点 v を出る
                yield (v, 0)
                if v == root:
                    break
                v = pa[v]
            else:  # 進める
                S[v] = s + 1
                v = c[s]
                yield (v, 1)

    def top_down(self):
        """木の根から (深さの浅い順に) 頂点を yield する."""
        assert self.__after_seal_check()
        if not hasattr(self, "tower"):
            self.depth_search(False)

        yield from self.bfs_order

    def bottom_up(self):
        """木の葉から (深さの深い順に) 頂点を yield する."""
        assert self.__after_seal_check()
        if not hasattr(self, "tower"):
            self.depth_search(False)

        for E in self.tower[::-1]:
            yield from E

    def tree_dp_from_leaf(self, merge, unit, f, g, Mode=False):
        """
        葉から木 DP を行う.

        [input]
        merge: 可換モノイドを成す2項演算 M x M -> M
        unit: M の単位元
        f: X x V x V -> M: f(x, v, w): v が親, w が子
        g: M x V -> X: g(x, v)
        Mode: False -> 根の値のみ, True -> 全ての値 (リスト)

        [補足]
        頂点 v の子が x, y, z, ... のとき, 更新式は * を merge として
            dp[v] = g(f(dp[x], v, x) * f(dp[y], v, y) * f(dp[z], v, z) * ..., v)
        になる.
        計算量: O(N) 回の merge, f, g 呼び出し
        """
        assert self.__after_seal_check()
        if not hasattr(self, "depth"):
            self.depth_search(False)

        data = [unit] * (self.index + self.N)
        ch = self.children
        order = self.bfs_order

        for i in range(len(order) - 1, -1, -1):
            x = order[i]
            a = unit
            for y in ch[x]:
                a = merge(a, f(data[y], x, y))
            data[x] = g(a, x)

        if Mode:
            return data
        else:
            return data[self.root]

    def tree_dp_from_root(self, f, alpha):
        """
        根から木 DP を行う.

        [input]
        alpha: 根の値
        f: X x V x V -> X: f(x, v, w): v が親, w が子

        [補足]
        頂点 v の親が x のとき, 更新式は
            dp[v] = f(dp[x], x, v) (v != root), alpha (v == root)
        になる.
        戻り値: dp のリスト
        計算量: O(N) 回の f 呼び出し
        """
        assert self.__after_seal_check()
        if not hasattr(self, "depth"):
            self.depth_search(False)

        data = [0] * (self.index + self.N)
        ch = self.children

        data[self.root] = alpha
        for x in self.bfs_order:
            dx = data[x]
            for y in ch[x]:
                data[y] = f(dx, x, y)

        return data

    def rerooting(self, merge, unit, f, g):
        """
        全方位木 DP を行う.

        [input]
        merge: 可換モノイドを成す2項演算 M x M -> M
        unit: M の単位元
        f: X x V x V -> M: f(x, v, w): v が親, w が子
        g: M x V -> X: g(x, v)

        ※ tree_dp_from_leaf と同じ形式

        [補足]
        頂点 v の子が x, y, z, ... のとき, 更新式は
            dp[v] = g(f(dp[x], v, x) * f(dp[y], v, y) * f(dp[z], v, z) * ..., v)
        になる. 各頂点を根としたときの dp[根] のリストを返す.
        計算量: O(N) 回の merge, f, g 呼び出し
        """
        assert self.__after_seal_check()

        size = self.index + self.N
        root = self.root
        ch = self.children
        pa = self.parent

        # 葉から (lower[v]: v を根とする部分木の dp)
        lower = self.tree_dp_from_leaf(merge, unit, f, g, True)

        # 根から (upper[c]: c の親 v を根とし, c 側を除いた部分の dp)
        upper = [unit] * size
        A = [unit] * size
        for v in self.bfs_order:
            cc = ch[v]
            deg = len(cc)
            fl = [f(lower[c], v, c) for c in cc]

            R = [unit] * (deg + 1)
            y = unit
            for i in range(deg - 1, -1, -1):
                y = merge(y, fl[i])
                R[i] = y

            x = f(upper[v], v, pa[v]) if v != root else unit
            for i in range(deg):
                upper[cc[i]] = g(merge(x, R[i + 1]), v)
                x = merge(x, fl[i])
            A[v] = g(x, v)
        return A

    def euler_tour_vertex(self, order=None):
        """
        オイラーツアー (vertex) に関する計算を行う.

        euler_vertex: 長さ 2N-1 の頂点列
        in_time[v], out_time[v]: v が euler_vertex に最初 / 最後に現れる位置
        order: 子の順番を指定する key 関数 (破壊的)
        計算量: O(N)
        """
        assert self.__after_seal_check()
        if hasattr(self, "euler_vertex"):
            return

        size = self.index + self.N
        L = 2 * self.N - 1
        X = [-1] * L

        ch = self.children
        if order is not None:
            for i in range(self.index, size):
                ch[i].sort(key=order)

        pa = self.parent
        S = [0] * size
        in_time = [-1] * size
        out_time = [-1] * size

        v = self.root
        for t in range(L):
            X[t] = v
            if in_time[v] == -1:
                in_time[v] = t
            out_time[v] = t
            c = ch[v]
            s = S[v]
            if s == len(c):
                v = pa[v]
            else:
                S[v] = s + 1
                v = c[s]

        self.euler_vertex = X
        self.in_time = in_time
        self.out_time = out_time

    def euler_tour_edge(self):
        """
        オイラーツアー (edge) に関する計算を行う.

        euler_edge[t] = (u, v, k): u から v へ向かう
        (k=+1 のときは葉へ進む向き, k=-1 のときは根へ進む向き)
        計算量: O(N)
        """
        assert self.__after_seal_check()
        if hasattr(self, "euler_edge"):
            return

        if not hasattr(self, "euler_vertex"):
            self.euler_tour_vertex()

        euler = self.euler_vertex
        pa = self.parent
        L = 2 * (self.N - 1)
        Ed = [0] * L
        for t in range(L):
            u = euler[t]
            v = euler[t + 1]
            Ed[t] = (u, v, 1 if u == pa[v] else -1)
        self.euler_edge = Ed

    def centroid(self, all=False):
        """
        木の重心を求める.

        all: False -> 重心のうちの1頂点. True -> 全ての重心のリスト (頂点番号昇順).
        計算量: O(N)
        """
        assert self.__after_seal_check()

        N = self.N
        M = N // 2

        if not hasattr(self, "des_count"):
            self.__descendant_count()

        G = []
        ch = self.children
        des = self.des_count

        for v in range(self.index, self.index + N):
            if N - des[v] > M:
                continue

            for x in ch[v]:
                if des[x] > M:
                    break
            else:
                if all:
                    G.append(v)
                else:
                    return v
        return G

    def generated_subtree(self, S):
        """
        S を含む最小の部分木の頂点を求める.

        戻り値: 頂点のリスト (昇順)
        計算量: O(N + |S| log N)
        """
        assert self.__after_seal_check(*S)

        S = list(set(S))
        if not S:
            return []

        w = S[0]
        for s in S:
            w = self.lowest_common_ancestor(w, s)

        pa = self.parent
        mark = [False] * (self.index + self.N)
        mark[w] = True
        T = [w]
        for s in S:
            while not mark[s]:
                mark[s] = True
                T.append(s)
                s = pa[s]
        T.sort()
        return T

    def generated_subtree_size(self, S):
        """
        S を含む最小の部分木のサイズを求める.

        ※ 戻り値は [辺の数] (頂点数はこれ +1).
        計算量: O(N + |S| log N)
        """
        assert self.__after_seal_check(*S)

        if not hasattr(self, "in_time"):
            self.euler_tour_vertex()

        in_time = self.in_time
        S = sorted(set(S), key=in_time.__getitem__)
        K = len(S)

        dist = self.distance
        X = dist(S[-1], S[0])
        for i in range(K - 1):
            X += dist(S[i], S[i + 1])
        return X // 2


# =================================================
def Making_Tree(N, E, root, index=0):
    """
    辺リストから根付き木を作る.

    N: 頂点数
    E: 辺のリスト [(u, v), ...] (無向)
    root: 根
    index: 頂点番号の開始値
    戻り値: seal 済みの Tree
    計算量: O(N)
    """
    F = [[] for _ in range(index + N)]
    for u, v in E:
        assert index <= u < index + N
        assert index <= v < index + N
        assert u != v

        F[u].append(v)
        F[v].append(u)

    X = [-1] * (index + N)
    X[root] = root

    Q = [root]
    for x in Q:
        for y in F[x]:
            if X[y] == -1:
                X[y] = x
                Q.append(y)

    T = Tree(N, index)
    T.root_set(root)
    T.parent = X
    T.seal()
    return T
# ==================================================


if __name__ == "__main__":
    T = Making_Tree(5, [(0, 1), (0, 2), (1, 3), (1, 4)], 0)
    print(T.lowest_common_ancestor(3, 4))  # 1
    print(T.distance(3, 2))  # 3
    print(T.diameter())  # (3, (3, 2))
    print(T.rerooting(max, 0, lambda x, v, w: x + 1, lambda x, v: x))  # [2, 2, 3, 3, 3]
