"""
遅延評価セグメントツリー (ACL の lazy_segtree の移植)

モノイド (ope, e) と作用 (mapping, composition, id_) を与えて、
区間作用・区間積をそれぞれ O(log N) で処理する。
よく使う設定は下の辞書 (range_add_range_min など) を ** で展開して渡せばよい。

使用例:
    A = [3, 1, 4, 1, 5]
    seg = lazy_segtree(A, **range_add_range_min)
    seg.apply(1, 4, 10)     # A[1:4] に 10 を加算 -> [3, 11, 14, 11, 5]
    seg.prod(0, 5)          # -> 3
    seg.get(2)              # -> 14

    # 区間和系は要素を (値, 区間長) のタプルで持つ
    seg = lazy_segtree([(a, 1) for a in A], **range_add_range_sum)
    seg.apply(0, 3, 2)
    seg.prod(1, 4)[0]       # -> (1+2) + (4+2) + 1 = 10
"""

_INF = 1 << 60

# 区間更新区間最大 (id_ = -(1<<60) は「更新なし」を表す番兵。この値への更新はできない)
range_update_range_max = {
    'ope': lambda x, y: x if x > y else y,
    'e': -_INF,
    'mapping': lambda f, x: x if f == -_INF else f,
    'composition': lambda f, g: f if f != -_INF else g,  # g -> f の順に作用
    'id_': -_INF,
}

# 区間更新区間最小 (id_ = 1<<60 は「更新なし」を表す番兵。この値への更新はできない)
range_update_range_min = {
    'ope': lambda x, y: x if x < y else y,
    'e': _INF,
    'mapping': lambda f, x: x if f == _INF else f,
    'composition': lambda f, g: f if f != _INF else g,  # g -> f の順に作用
    'id_': _INF,
}

# 区間加算区間最大
range_add_range_max = {
    'ope': lambda x, y: x if x > y else y,
    'e': -_INF,
    'mapping': lambda f, x: x + f,
    'composition': lambda f, g: f + g,  # g -> f の順に作用
    'id_': 0,
}

# 区間加算区間最小
range_add_range_min = {
    'ope': lambda x, y: x if x < y else y,
    'e': _INF,
    'mapping': lambda f, x: x + f,
    'composition': lambda f, g: f + g,  # g -> f の順に作用
    'id_': 0,
}

# 区間加算区間和: 要素は (和, 区間長)。初期値は [(a, 1) for a in A]
range_add_range_sum = {
    'ope': lambda x, y: (x[0] + y[0], x[1] + y[1]),
    'e': (0, 0),
    'mapping': lambda f, x: (x[0] + f * x[1], x[1]),
    'composition': lambda f, g: f + g,  # g -> f の順に作用
    'id_': 0,
}

# 区間更新区間和: 要素は (和, 区間長)。初期値は [(a, 1) for a in A]
# (id_ = 1<<60 は「更新なし」を表す番兵。この値への更新はできない)
range_update_range_sum = {
    'ope': lambda x, y: (x[0] + y[0], x[1] + y[1]),
    'e': (0, 0),
    'mapping': lambda f, x: (f * x[1], x[1]) if f != _INF else x,
    'composition': lambda f, g: f if f != _INF else g,  # g -> f の順に作用
    'id_': _INF,
}


def _arith_ope(x, y):
    """隣接する等差数列加算用ノードを結合する。"""
    lx, rx, sx = x
    ly, ry, sy = y
    return (lx if lx < ly else ly, rx if rx > ry else ry, sx + sy)


def _arith_mapping(f, x):
    """線式 a*i+b の加算を区間和ノードに適用する。"""
    a, b = f
    l, r, s = x
    if l >= r:
        return x
    return (l, r, s + (r - l) * (a * (l + r - 1) + 2 * b) // 2)


# 区間等差数列加算区間和
# 要素は (区間左端 l, 区間右端 r, 和)。初期値は [(i, i + 1, a) for i, a in enumerate(A)]
# 作用 (a, b) は各 i ∈ [l, r) について A[i] += a * i + b
# (添字 i そのものを使うので、[L, R) に初項 c・公差 d を足したいときは a = d, b = c - d * L)
range_arith_add_range_sum = {
    'ope': _arith_ope,
    'e': (_INF, -_INF, 0),
    'mapping': _arith_mapping,
    'composition': lambda f, g: (f[0] + g[0], f[1] + g[1]),  # g -> f の順に作用
    'id_': (0, 0),
}


def make_range_affine_range_sum(mod=998244353):
    """
    区間一次関数作用区間和 (mod 付き) の設定辞書を作る

    要素は (和, 区間長)。初期値は [(a, 1) for a in A]
    作用 (a, b) は各要素を x -> a * x + b に変える

    mod: 法
    戻り値: lazy_segtree に ** で渡す辞書
    """
    def ope(x, y):
        """区間和と区間長を結合する。"""
        return ((x[0] + y[0]) % mod, x[1] + y[1])

    def mapping(f, x):
        """一次関数 f を区間和ノード x に適用する。"""
        return ((f[0] * x[0] + f[1] * x[1]) % mod, x[1])

    def composition(f, g):  # g -> f の順に作用
        """一次関数 g の後に f を適用する合成作用を返す。"""
        return (f[0] * g[0] % mod, (f[0] * g[1] + f[1]) % mod)

    return {
        'ope': ope,
        'e': (0, 0),
        'mapping': mapping,
        'composition': composition,
        'id_': (1, 0),
    }


class lazy_segtree:
    """
    遅延評価セグメントツリー

    ope(x, y): モノイドの二項演算
    e: ope の単位元
    mapping(f, x): 作用 f を要素 x に適用した結果
    composition(f, g): 作用 g の後に f を行う合成作用
    id_: 恒等作用
    """

    def __init__(self, lst, ope, e, mapping, composition, id_):
        """
        初期化 O(N)

        lst: 配列の初期値
        """
        n = len(lst)
        self.n = n
        self.log = (n - 1).bit_length()
        size = 1 << self.log
        self.size = size
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
        """ノード k の値を子から再計算する"""
        data = self.data
        data[k] = self.op(data[2 * k], data[2 * k + 1])

    def all_apply(self, k, f):
        """ノード k に作用 f を適用し、内部ノードなら遅延値に積む"""
        self.data[k] = self.mapping(f, self.data[k])
        if k < self.size:
            self.lz[k] = self.composition(f, self.lz[k])

    def push(self, k):
        """ノード k (内部ノード) の遅延値を子に伝播する"""
        lz = self.lz
        f = lz[k]
        if f == self.identity:
            return
        data = self.data
        mapping = self.mapping
        a = 2 * k
        b = a + 1
        data[a] = mapping(f, data[a])
        data[b] = mapping(f, data[b])
        if a < self.size:
            composition = self.composition
            lz[a] = composition(f, lz[a])
            lz[b] = composition(f, lz[b])
        lz[k] = self.identity

    def _push_path(self, p):
        """葉 p (内部添字) の祖先を根から順に push する"""
        push = self.push
        for i in range(self.log, 0, -1):
            push(p >> i)

    def _update_path(self, p):
        """葉 p (内部添字) の祖先を下から順に再計算する"""
        data = self.data
        op = self.op
        p >>= 1
        while p:
            data[p] = op(data[2 * p], data[2 * p + 1])
            p >>= 1

    def set(self, p, x):
        """
        p 番目の値を x に更新 O(log N)

        p: index(0-index)
        x: update value
        """
        p += self.size
        self._push_path(p)
        self.data[p] = x
        self._update_path(p)

    def get(self, p):
        """
        p 番目の値を返す O(log N)

        p: index(0-index)
        """
        p += self.size
        self._push_path(p)
        return self.data[p]

    def prod(self, l, r):
        """
        区間 [l, r) の ope による積を返す O(log N)

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
        """全体の積を返す O(1)"""
        return self.data[1]

    def apply_point(self, p, f):
        """
        p 番目の要素に作用 f を適用 O(log N)

        p: index(0-index)
        f: 作用
        """
        p += self.size
        self._push_path(p)
        self.data[p] = self.mapping(f, self.data[p])
        self._update_path(p)

    def apply(self, l, r, f):
        """
        区間 [l, r) の各要素に作用 f を適用 O(log N)

        l: index(0-index)
        r: index(0-index)
        f: 作用
        """
        if l == r:
            return
        l += self.size
        r += self.size
        log = self.log
        push = self.push
        for i in range(log, 0, -1):
            if (l >> i) << i != l:
                push(l >> i)
            if (r >> i) << i != r:
                push((r - 1) >> i)
        data = self.data
        lz = self.lz
        size = self.size
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
        g(prod(l, r)) が True となる最大の r を返す O(log N)

        g は単調 (g(e) = True で、ある所から先はずっと False) であること
        l: index(0-index)
        g: 判定関数
        """
        n = self.n
        if l == n:
            return n
        size = self.size
        l += size
        self._push_path(l)
        data = self.data
        op = self.op
        push = self.push
        sm = self.e
        while True:
            while l & 1 == 0:
                l >>= 1
            if not g(op(sm, data[l])):
                while l < size:
                    push(l)
                    l *= 2
                    t = op(sm, data[l])
                    if g(t):
                        sm = t
                        l += 1
                return l - size
            sm = op(sm, data[l])
            l += 1
            if l & -l == l:
                break
        return n

    def min_left(self, r, g):
        """
        g(prod(l, r)) が True となる最小の l を返す O(log N)

        g は単調 (g(e) = True で、ある所から先はずっと False) であること
        r: index(0-index)
        g: 判定関数
        """
        if r == 0:
            return 0
        size = self.size
        r += size
        push = self.push
        for i in range(self.log, 0, -1):
            push((r - 1) >> i)
        data = self.data
        op = self.op
        sm = self.e
        while True:
            r -= 1
            while r > 1 and r & 1:
                r >>= 1
            if not g(op(data[r], sm)):
                while r < size:
                    push(r)
                    r = 2 * r + 1
                    t = op(data[r], sm)
                    if g(t):
                        sm = t
                        r -= 1
                return r + 1 - size
            sm = op(data[r], sm)
            if r & -r == r:
                break
        return 0


if __name__ == "__main__":
    A = [3, 1, 4, 1, 5]
    seg = lazy_segtree(A, **range_add_range_min)
    seg.apply(1, 4, 10)
    print(seg.prod(0, 5), seg.get(2))  # 3 14
    seg = lazy_segtree([(a, 1) for a in A], **range_add_range_sum)
    seg.apply(0, 3, 2)
    print(seg.prod(1, 4)[0])  # 10
