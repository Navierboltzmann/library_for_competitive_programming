"""
HL分解 (Heavy-Light Decomposition) と、それを使った木上のセグ木・遅延セグ木。

- HLD             : 頂点にIDを振りなおし、パス・部分木を ID 上の半開区間の集まりに変換する
- HLD_SegTree     : 木上で パス/部分木 の区間取得 + 一点更新
- HLD_LazySegTree : 木上で パス/部分木 の区間作用 + 区間取得
- SegTree, lazy_segtree : 上の2つが内部で使うセグ木 (単体でも使える)

使用例:
    e = [[1, 2], [0, 3], [0], [1]]      # 隣接リスト (0-indexed, 木)
    A = [5, 3, 7, 1]
    T = HLD_SegTree(e, A, min, 1 << 60)
    T.path_query(3, 2)                  # 3-1-0-2 上の min = 1
    T.update(3, 10)
    T.sub_query(1)                      # 頂点1の部分木 {1, 3} の min = 3

lazy_segtree / HLD_LazySegTree に渡す演算の例:

    # 区間加算・区間最小値
    ope = lambda x, y: min(x, y)
    ide_ele = 1 << 60
    mapping = lambda f, x: x + f
    composition = lambda f, g: f + g    # g → f の順に作用
    id_ = 0

    # 区間加算・区間和 (値は (和, 要素数))
    ope = lambda x, y: (x[0] + y[0], x[1] + y[1])
    ide_ele = (0, 0)
    mapping = lambda f, x: (x[0] + f * x[1], x[1])
    composition = lambda f, g: f + g
    id_ = 0

    # 区間等差数列加算・区間和 (値は (左端, 右端(開), 和)、作用 (a, b) は位置 i に a*i+b を加算)
    def ope(x, y):
        return (min(x[0], y[0]), max(x[1], y[1]), x[2] + y[2])
    ide_ele = (1 << 30, -(1 << 30), 0)
    def mapping(f, x):
        a, b = f
        l, r, s = x
        return (l, r, s + (r - l) * (a * l + a * r - a + 2 * b) // 2)
    def composition(f, g):
        return (f[0] + g[0], f[1] + g[1])
    id_ = (0, 0)
"""


class SegTree:
    """
    Segment Tree (一点更新・区間取得)
    """

    def __init__(self, init_val, segfunc, ide_ele):
        """
        初期化 O(N)

        init_val: 配列の初期値
        segfunc: 二項演算 (結合的)
        ide_ele: segfunc の単位元
        """
        n = len(init_val)
        self.segfunc = segfunc
        self.ide_ele = ide_ele
        self.num = num = 1 << (n - 1).bit_length()
        tree = [ide_ele] * (2 * num)
        tree[num:num + n] = init_val
        for i in range(num - 1, 0, -1):
            tree[i] = segfunc(tree[2 * i], tree[2 * i + 1])
        self.tree = tree

    def update(self, k, x):
        """
        k番目の値をxに更新 O(logN)

        k: index(0-index)
        x: update value
        """
        tree = self.tree
        segfunc = self.segfunc
        k += self.num
        tree[k] = x
        k >>= 1
        while k:
            tree[k] = segfunc(tree[2 * k], tree[2 * k + 1])
            k >>= 1

    def get(self, k):
        """k番目の値を返す O(1)"""
        return self.tree[k + self.num]

    def query(self, l, r):
        """
        [l, r)のsegfuncしたものを得る O(logN)

        l: index(0-index)
        r: index(0-index)
        """
        tree = self.tree
        segfunc = self.segfunc
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


class lazy_segtree:
    """
    遅延評価セグ木 (ACL の lazy_segtree 相当)

    ope(x, y): 値の二項演算, e: ope の単位元,
    mapping(f, x): 作用 f を値 x に適用, composition(f, g): g → f の順に作用させる合成,
    id_: 恒等作用。 演算の例はモジュール docstring を参照。
    """

    def __init__(self, lst, ope, e, mapping, composition, id_):
        """
        初期化 O(N)

        lst: 配列の初期値
        """
        self.n = n = len(lst)
        self.log = (n - 1).bit_length()
        self.size = size = 1 << self.log
        data = [e] * (2 * size)
        data[size:size + n] = lst
        for i in range(size - 1, 0, -1):
            data[i] = ope(data[2 * i], data[2 * i + 1])
        self.data = data
        self.lz = [id_] * size
        self.e = e
        self.op = ope
        self.mapping = mapping
        self.composition = composition
        self.identity = id_

    def update(self, k):
        """節点 k の値を子から再計算する"""
        self.data[k] = self.op(self.data[2 * k], self.data[2 * k + 1])

    def all_apply(self, k, f):
        """節点 k に作用 f を適用する (内部節点なら遅延値にも積む)"""
        self.data[k] = self.mapping(f, self.data[k])
        if k < self.size:
            self.lz[k] = self.composition(f, self.lz[k])

    def push(self, k):
        """節点 k の遅延値を子に伝播する"""
        lz = self.lz
        f = lz[k]
        data = self.data
        mapping = self.mapping
        a = 2 * k
        data[a] = mapping(f, data[a])
        data[a + 1] = mapping(f, data[a + 1])
        if a < self.size:
            composition = self.composition
            lz[a] = composition(f, lz[a])
            lz[a + 1] = composition(f, lz[a + 1])
        lz[k] = self.identity

    def _push_path(self, p):
        """葉 p (size 加算済み) の祖先の遅延値を根から順に伝播する"""
        push = self.push
        for i in range(self.log, 0, -1):
            push(p >> i)

    def _update_path(self, p):
        """葉 p (size 加算済み) の祖先を下から順に再計算する"""
        data = self.data
        op = self.op
        p >>= 1
        while p:
            data[p] = op(data[2 * p], data[2 * p + 1])
            p >>= 1

    def set(self, p, x):
        """
        p番目の値をxに置き換える O(logN)

        p: index(0-index)
        x: 新しい値
        """
        p += self.size
        self._push_path(p)
        self.data[p] = x
        self._update_path(p)

    def get(self, p):
        """p番目の値を返す O(logN)"""
        p += self.size
        self._push_path(p)
        return self.data[p]

    def prod(self, l, r):
        """
        [l, r) を ope でまとめた値を返す O(logN)

        l: index(0-index)
        r: index(0-index)
        """
        if l == r:
            return self.e
        l += self.size
        r += self.size
        push = self.push
        for i in range(self.log, 0, -1):
            if (l >> i) << i != l:
                push(l >> i)
            if (r >> i) << i != r:
                push((r - 1) >> i)
        data = self.data
        op = self.op
        sml = smr = self.e
        while l < r:
            if l & 1:
                sml = op(sml, data[l])
                l += 1
            if r & 1:
                r -= 1
                smr = op(data[r], smr)
            l >>= 1
            r >>= 1
        return op(sml, smr)

    def all_prod(self):
        """全体を ope でまとめた値を返す O(1)"""
        return self.data[1]

    def apply_point(self, p, f):
        """
        p番目の値に作用 f を適用する O(logN)

        p: index(0-index)
        f: 作用
        """
        p += self.size
        self._push_path(p)
        self.data[p] = self.mapping(f, self.data[p])
        self._update_path(p)

    def apply(self, l, r, f):
        """
        [l, r) の各値に作用 f を適用する O(logN)

        l: index(0-index)
        r: index(0-index)
        f: 作用
        """
        if l == r:
            return
        size = self.size
        log = self.log
        l += size
        r += size
        push = self.push
        for i in range(log, 0, -1):
            if (l >> i) << i != l:
                push(l >> i)
            if (r >> i) << i != r:
                push((r - 1) >> i)
        data = self.data
        lz = self.lz
        mapping = self.mapping
        composition = self.composition
        l2, r2 = l, r
        while l < r:
            if l & 1:
                data[l] = mapping(f, data[l])
                if l < size:
                    lz[l] = composition(f, lz[l])
                l += 1
            if r & 1:
                r -= 1
                data[r] = mapping(f, data[r])
                if r < size:
                    lz[r] = composition(f, lz[r])
            l >>= 1
            r >>= 1
        l, r = l2, r2
        op = self.op
        for i in range(1, log + 1):
            if (l >> i) << i != l:
                k = l >> i
                data[k] = op(data[2 * k], data[2 * k + 1])
            if (r >> i) << i != r:
                k = (r - 1) >> i
                data[k] = op(data[2 * k], data[2 * k + 1])

    def max_right(self, l, g):
        """
        g(prod(l, r)) が True となる最大の r を返す O(logN)

        l: index(0-index)
        g: 単調な判定関数 (g(e) は True であること)
        """
        if l == self.n:
            return self.n
        l += self.size
        self._push_path(l)
        data = self.data
        op = self.op
        size = self.size
        sm = self.e
        while True:
            while l % 2 == 0:
                l >>= 1
            if not g(op(sm, data[l])):
                while l < size:
                    self.push(l)
                    l *= 2
                    if g(op(sm, data[l])):
                        sm = op(sm, data[l])
                        l += 1
                return l - size
            sm = op(sm, data[l])
            l += 1
            if l & -l == l:
                break
        return self.n

    def min_left(self, r, g):
        """
        g(prod(l, r)) が True となる最小の l を返す O(logN)

        r: index(0-index)
        g: 単調な判定関数 (g(e) は True であること)
        """
        if r == 0:
            return 0
        r += self.size
        push = self.push
        for i in range(self.log, 0, -1):
            push((r - 1) >> i)
        data = self.data
        op = self.op
        size = self.size
        sm = self.e
        while True:
            r -= 1
            while r > 1 and r % 2 == 1:
                r >>= 1
            if not g(op(data[r], sm)):
                while r < size:
                    push(r)
                    r = 2 * r + 1
                    if g(op(data[r], sm)):
                        sm = op(data[r], sm)
                        r -= 1
                return r + 1 - size
            sm = op(data[r], sm)
            if r & -r == r:
                break
        return 0


class HLD:
    """HL分解。頂点にIDを振りなおし、パス・部分木をID上の半開区間の集まりとして返す。
    構築 O(N)、path_query / lca は O(logN) 個の区間・O(logN) 時間。

    SegTreeに載せる配列は reorder で A[ID[v]] = val[v] と並べ替えたもの。
    辺に値を持たせる場合は、辺(par[v], v)の値を頂点vに置き edge=True でクエリする。

    使い方:
        hld = HLD(e)                          # e: 隣接リスト (0-indexed, 木)
        seg = SegTree(hld.reorder(A), min, INF)
        res = INF
        for l, r in hld.path_query(u, v):     # u-v パス上の min
            res = min(res, seg.query(l, r))
        seg.query(*hld.sub_query(u))          # u の部分木の min
        seg.update(hld.ID[u], x)              # 頂点 u の値を x に

    属性 (すべて元の頂点番号で引く):
        ID[v]    : 振りなおしたID
        par[v]   : 親 (根は -1)
        depth[v] : 根からの深さ
        size[v]  : 部分木の頂点数
        head[v]  : v を含む heavy path の一番上の頂点
        heavy[v] : 部分木が最大の子 (葉は -1)
    """

    def __init__(self, e, root=0):
        """
        構築 O(N)

        e (list of list of int): 木の隣接リスト。
        root (int): 根。
        """
        n = len(e)
        self.N = n
        self.e = e
        self.root = root

        par = [-1] * n
        depth = [0] * n
        order = [root]
        for x in order:
            px = par[x]
            dx = depth[x] + 1
            for y in e[x]:
                if y != px:
                    par[y] = x
                    depth[y] = dx
                    order.append(y)

        size = [1] * n
        heavy = [-1] * n
        for i in range(len(order) - 1, 0, -1):
            x = order[i]
            p = par[x]
            size[p] += size[x]
            h = heavy[p]
            if h == -1 or size[h] < size[x]:
                heavy[p] = x

        # heavy child を先に訪れる pre-order でIDを振ると、heavy path と部分木がそれぞれ連続区間になる
        ID = [0] * n
        head = [root] * n
        stack = [root]
        cnt = 0
        while stack:
            x = stack.pop()
            ID[x] = cnt
            cnt += 1
            px = par[x]
            hx = heavy[x]
            for y in e[x]:
                if y != px and y != hx:
                    head[y] = y
                    stack.append(y)
            if hx != -1:
                head[hx] = head[x]
                stack.append(hx)

        self.par = par
        self.depth = depth
        self.size = size
        self.heavy = heavy
        self.head = head
        self.ID = ID

    def reorder(self, val):
        """頂点番号順の配列 val を ID 順に並べ替えた配列を返す。SegTree の初期値に使う。O(N)"""
        A = [None] * self.N
        for v, i in enumerate(self.ID):
            A[i] = val[v]
        return A

    def path_query(self, u, v, edge=False):
        """u-v パス (両端を含む) を覆う ID 上の区間 [l, r) のリストを返す。O(logN)
        区間の順番はパス上の順番と一致しないので、載せる演算は可換であること。

        u, v (int): パスの両端 (元の頂点番号)。
        edge (bool): True なら LCA を除く (辺に値を持たせる場合)。
        """
        ID, head, par = self.ID, self.head, self.par
        res = []
        while head[u] != head[v]:
            if ID[head[u]] < ID[head[v]]:
                u, v = v, u
            hu = head[u]
            res.append((ID[hu], ID[u] + 1))
            u = par[hu]
        if ID[u] > ID[v]:
            u, v = v, u
        if ID[u] + edge <= ID[v]:
            res.append((ID[u] + edge, ID[v] + 1))
        return res

    def sub_query(self, k, edge=False):
        """k の部分木に対応する ID 上の区間 [l, r) を返す。O(1)
        edge=True なら k 自身 (辺(par[k], k)) を除く。
        """
        return (self.ID[k] + edge, self.ID[k] + self.size[k])

    def lca(self, u, v):
        """u と v の最小共通祖先 O(logN)"""
        ID, head, par = self.ID, self.head, self.par
        while head[u] != head[v]:
            if ID[head[u]] < ID[head[v]]:
                u, v = v, u
            u = par[head[u]]
        return u if self.depth[u] < self.depth[v] else v

    def dist(self, u, v):
        """u-v 間の辺の本数 O(logN)"""
        return self.depth[u] + self.depth[v] - 2 * self.depth[self.lca(u, v)]


class HLD_SegTree:
    """木の上の SegTree。パス・部分木の区間取得と一点更新を O(log^2 N) / O(logN) で行う。
    頂点はすべて元の頂点番号で指定する。segfunc は可換であること。

    使い方:
        T = HLD_SegTree(e, A, lambda x, y: x + y, 0)
        T.path_query(u, v)   # u-v パス上の和
        T.sub_query(u)       # u の部分木の和
        T.update(u, x)       # 頂点 u の値を x に
        T.get(u)             # 頂点 u の値

    辺に値を持たせる場合: A[v] = 辺(par[v], v) の値 (A[root] は単位元) として edge=True。
    """

    def __init__(self, e, init_val, segfunc, ide_ele, root=0, edge=False):
        """
        構築 O(N)

        e (list of list of int): 木の隣接リスト。
        init_val (list): 頂点番号順の初期値。
        segfunc (function): 可換な二項演算。
        ide_ele (any): segfunc の単位元。
        root (int): 根。
        edge (bool): True なら辺に値を持たせる。
        """
        self.hld = HLD(e, root=root)
        self.ID = self.hld.ID
        self.N = len(e)
        self.edge = edge
        self.segfunc = segfunc
        self.ide_ele = ide_ele
        self.seg = SegTree(self.hld.reorder(init_val), segfunc, ide_ele)

    def update(self, x, val):
        """頂点 x の値を val に更新する。O(logN)
        edge=True のときは辺(par[x], x) の値を更新する。
        """
        self.seg.update(self.ID[x], val)

    def get(self, x):
        """頂点 x の値を返す。O(1)
        edge=True のときは辺(par[x], x) の値を返す。
        """
        return self.seg.get(self.ID[x])

    def path_query(self, u, v):
        """u-v パス上の値を segfunc でまとめた結果を返す。O(log^2 N)
        edge=False なら両端を含むパス上の頂点、edge=True ならパス上の辺が対象。
        u == v のとき、edge=True なら単位元を返す。
        """
        segfunc = self.segfunc
        query = self.seg.query
        res = self.ide_ele
        for l, r in self.hld.path_query(u, v, self.edge):
            res = segfunc(res, query(l, r))
        return res

    def sub_query(self, x):
        """x の部分木の値を segfunc でまとめた結果を返す。O(logN)
        edge=False なら x を含む部分木の頂点、edge=True なら部分木内の辺 (辺(par[x], x) は含まない) が対象。
        """
        return self.seg.query(*self.hld.sub_query(x, self.edge))


class HLD_LazySegTree:
    """木の上の遅延セグ木。パス・部分木への区間作用と区間取得を O(log^2 N) / O(logN) で行う。
    頂点はすべて元の頂点番号で指定する。ope は可換であること。
    ope, ide_ele, mapping, composition, id_ は lazy_segtree と同じ (モジュール docstring の例を参照)。

    使い方 (区間加算・区間和):
        T = HLD_LazySegTree(e, [(a, 1) for a in A],
                            lambda x, y: (x[0] + y[0], x[1] + y[1]), (0, 0),
                            lambda f, x: (x[0] + f * x[1], x[1]),
                            lambda f, g: f + g, 0)
        T.path_apply(u, v, f)    # u-v パス上の頂点すべてに f を作用
        T.sub_apply(u, f)        # u の部分木すべてに f を作用
        T.path_query(u, v)[0]    # u-v パス上の和
        T.sub_query(u)[0]        # u の部分木の和
        T.set(u, (x, 1))         # 頂点 u の値を x に
        T.get(u)[0]              # 頂点 u の値

    辺に値を持たせる場合: lst[v] = 辺(par[v], v) の値 (lst[root] は単位元) として edge=True。
    """

    def __init__(self, e, lst, ope, ide_ele, mapping, composition, id_, root=0, edge=False):
        """
        構築 O(N)

        e (list of list of int): 木の隣接リスト。
        lst (list): 頂点番号順の初期値。
        ope (function): 可換な二項演算。
        ide_ele (any): ope の単位元。
        mapping (function): mapping(f, x) = 作用 f を値 x に適用した結果。
        composition (function): composition(f, g) = g の後に f を作用させる合成。
        id_ (any): 恒等作用。
        root (int): 根。
        edge (bool): True なら辺に値を持たせる。
        """
        self.hld = HLD(e, root=root)
        self.ID = self.hld.ID
        self.N = len(e)
        self.edge = edge
        self.segfunc = ope
        self.ide_ele = ide_ele
        self.seg = lazy_segtree(self.hld.reorder(lst), ope, ide_ele, mapping, composition, id_)

    def set(self, x, val):
        """頂点 x の値を val に置き換える。O(logN)
        edge=True のときは辺(par[x], x) の値を置き換える。
        """
        self.seg.set(self.ID[x], val)

    def get(self, x):
        """頂点 x の値 (それまでの作用をすべて反映したもの) を返す。O(logN)
        edge=True のときは辺(par[x], x) の値を返す。
        """
        return self.seg.get(self.ID[x])

    def path_query(self, u, v):
        """u-v パス上の値を ope でまとめた結果を返す。O(log^2 N)
        edge=False なら両端を含むパス上の頂点、edge=True ならパス上の辺が対象。
        u == v のとき、edge=True なら単位元を返す。
        """
        segfunc = self.segfunc
        prod = self.seg.prod
        res = self.ide_ele
        for l, r in self.hld.path_query(u, v, self.edge):
            res = segfunc(res, prod(l, r))
        return res

    def sub_query(self, x):
        """x の部分木の値を ope でまとめた結果を返す。O(logN)
        edge=False なら x を含む部分木の頂点、edge=True なら部分木内の辺 (辺(par[x], x) は含まない) が対象。
        """
        return self.seg.prod(*self.hld.sub_query(x, self.edge))

    def path_apply(self, u, v, f):
        """u-v パス上のすべてに作用 f を適用する。O(log^2 N)
        対象は path_query と同じ (edge=True ならパス上の辺)。
        """
        apply = self.seg.apply
        for l, r in self.hld.path_query(u, v, self.edge):
            apply(l, r, f)

    def sub_apply(self, x, f):
        """x の部分木のすべてに作用 f を適用する。O(logN)
        対象は sub_query と同じ (edge=True なら辺(par[x], x) は含まない)。
        """
        self.seg.apply(*self.hld.sub_query(x, self.edge), f)


if __name__ == "__main__":
    e = [[1, 2], [0, 3], [0], [1]]
    T = HLD_SegTree(e, [5, 3, 7, 1], min, 1 << 60)
    print(T.path_query(3, 2))  # 1
    T.update(3, 10)
    print(T.sub_query(1))  # 3
