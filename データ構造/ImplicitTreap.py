"""
ImplicitTreap: 列に対する挿入・削除・区間反転・区間作用・区間積を期待 O(log N) で行う平衡二分木

遅延セグ木(ACL の lazy_segtree)と同じ (op, e, mapping, composition, id) を渡して使う。
    op(x, y)          : 値のモノイド演算
    e                 : op の単位元
    mapping(f, x)     : 作用素 f を値 x に作用させた結果
    composition(f, g) : g を作用させた後に f を作用させる合成作用素(g → f の順)
    id                : 恒等作用素
reverse を使う場合、op は可換であること(反転しても区間積が変わらない前提)。

参考: https://judge.yosupo.jp/submission/31925

使用例(区間加算・区間最小値):
    tr = ImplicitTreap(min, 1 << 60, lambda f, x: x + f, lambda f, g: f + g, 0)
    tr.build([5, 1, 3])
    tr.insert(1, 4)        # [5, 4, 1, 3]
    tr.apply(0, 2, 10)     # [15, 14, 1, 3]
    tr.prod(0, 2)          # -> 14
    tr.reverse(0, 4)       # [3, 1, 14, 15]
    tr.erase(0)            # [1, 14, 15]
    tr.get(1)              # -> 14
    list(tr.iter())        # -> [1, 14, 15]

よく使う (op, e, mapping, composition, id) の例:
    # 区間加算・区間最小値
    ope = lambda x, y: min(x, y); ide_ele = 1 << 60
    mapping = lambda f, x: x + f; composition = lambda f, g: f + g; id_ = 0

    # 区間加算・区間和(値は (和, 長さ) の組で持つ。初期値は (a, 1))
    ope = lambda x, y: (x[0] + y[0], x[1] + y[1]); ide_ele = (0, 0)
    mapping = lambda f, x: (x[0] + f * x[1], x[1]); composition = lambda f, g: f + g; id_ = 0

    # 区間等差数列加算・区間和(値は (l, r, s) = (左端 index, 右端 index + 1, 和)。初期値は (i, i + 1, a))
    #   作用 (a, b): 位置 i に a * i + b を加算
    def ope(x, y):
        lx, rx, sx = x
        ly, ry, sy = y
        return (min(lx, ly), max(rx, ry), sx + sy)
    INF = 1 << 30
    ide_ele = (INF, -INF, 0)
    def mapping(f, x):
        a, b = f
        l, r, s = x
        s += (r - l) * (a * l + a * r - a + 2 * b) // 2
        return (l, r, s)
    def composition(f, g):
        return (f[0] + g[0], f[1] + g[1])
    id_ = (0, 0)

    # 区間一次関数作用・区間和(作用 (a, b): x -> a * x + b、値は (和, 長さ))
    mod = 10**9 + 7
    ope = lambda x, y: ((x[0] + y[0]) % mod, (x[1] + y[1]) % mod); ide_ele = (0, 0)
    mapping = lambda f, x: ((f[0] * x[0] + f[1] * x[1]) % mod, x[1])
    composition = lambda f, g: ((f[0] * g[0]) % mod, (f[0] * g[1] + f[1]) % mod)  # g → f の順
    id_ = (1, 0)
"""


class ImplicitTreap():
    """
    暗黙のキー(位置)を持つ Treap

    ノード 0 は番兵(空)。ノード t について
        val[t]: 値, pri[t]: 優先度(大きいほど根に近い), cnt[t]: 部分木サイズ,
        acc[t]: 部分木の値の積, laz[t]: 子へ未伝播の作用素,
        ptr[3t], ptr[3t+1]: 左右の子, ptr[3t+2]: 反転フラグ
    """
    __slots__ = ['root', 'op', 'e', 'mapping', 'composition', 'id', 'val', 'pri', 'ptr', 'cnt', 'acc', 'laz', 'rand']

    def __init__(self, op, e, mapping, composition, id):
        """
        op, e, mapping, composition, id: モジュール docstring 参照
        """
        self.root = 0
        self.op = op
        self.e = e
        self.mapping = mapping
        self.composition = composition
        self.id = id
        self.val = [e]
        self.pri = [1 << 32]
        self.ptr = [0, 0, 0]
        self.cnt = [0]
        self.acc = [e]
        self.laz = [id]
        self.rand = self._xor64_()

    def _xor64_(self):
        """優先度用の 32bit 擬似乱数を無限に生成するジェネレータ"""
        x = 88172645463325252
        while True:
            x = x ^ ((x << 7) & 0xffffffff)
            x = x ^ (x >> 9)
            yield x & 0xffffffff

    def _newnode_(self, x):
        """値 x のノードを1つ作り、その番号を返す。O(1)"""
        idx = len(self.val)
        self.val.append(x)
        self.pri.append(next(self.rand))
        self.ptr.extend((0, 0, 0))
        self.cnt.append(1)
        self.acc.append(self.e)
        self.laz.append(self.id)
        return idx

    def _push_(self, t):
        """ノード t の反転フラグ・遅延作用を子へ伝播する(t != 0)。O(1)"""
        ptr, laz = self.ptr, self.laz
        t3 = t * 3
        if ptr[t3 + 2]:
            ptr[t3 + 2] = 0
            l, r = ptr[t3], ptr[t3 + 1] = ptr[t3 + 1], ptr[t3]
            if l:
                ptr[l * 3 + 2] ^= 1
            if r:
                ptr[r * 3 + 2] ^= 1
        f = laz[t]
        if f != self.id:
            l, r = ptr[t3], ptr[t3 + 1]
            acc, mapping, composition = self.acc, self.mapping, self.composition
            if l:
                laz[l] = composition(f, laz[l])
                acc[l] = mapping(f, acc[l])
            if r:
                laz[r] = composition(f, laz[r])
                acc[r] = mapping(f, acc[r])
            self.val[t] = mapping(f, self.val[t])
            laz[t] = self.id

    def _update_(self, t):
        """ノード t の cnt, acc を子から再計算する。O(1)"""
        ptr, acc, op = self.ptr, self.acc, self.op
        l, r = ptr[t * 3], ptr[t * 3 + 1]
        self.cnt[t] = self.cnt[l] + self.cnt[r] + 1
        acc[t] = op(op(acc[l], self.val[t]), acc[r])

    def _split_(self, t, k, update=False):
        """
        木 t を先頭 k 個とそれ以外に分割し (l, r) を返す

        update=False のときは通過したノードの acc を再計算しない(後の _merge_ で再計算される)。
        計算量: 期待 O(log N)
        """
        ptr, cnt = self.ptr, self.cnt
        push = self._push_
        l = r = 0
        while t:
            push(t)
            p = cnt[ptr[t * 3]] + 1
            if k < p:
                v, ptr[t * 3] = ptr[t * 3], r
                r, t = t, v
            else:
                v, ptr[t * 3 + 1] = ptr[t * 3 + 1], l
                l, t = t, v
                k -= p
        upd = self._update_
        s = 0
        while l:
            v, ptr[l * 3 + 1] = ptr[l * 3 + 1], s
            if update:
                upd(l)
            s, l = l, v
        l = s
        s = 0
        while r:
            v, ptr[r * 3] = ptr[r * 3], s
            if update:
                upd(r)
            s, r = r, v
        r = s
        return l, r

    def _merge_(self, l, r, push_lt=False, push_rt=False):
        """
        木 l と木 r を連結した木の根を返す

        push_lt / push_rt: l の右端の経路 / r の左端の経路を辿る際に _push_ するか
        計算量: 期待 O(log N)
        """
        ptr, pri = self.ptr, self.pri
        push, upd = self._push_, self._update_
        s = 0
        while l:
            if push_lt:
                push(l)
            v, ptr[l * 3 + 1] = ptr[l * 3 + 1], s
            s, l = l, v
        l = s
        s = 0
        while r:
            if push_rt:
                push(r)
            v, ptr[r * 3] = ptr[r * 3], s
            s, r = r, v
        r = s
        t = 0
        while l or r:
            if pri[l] < pri[r]:
                v, ptr[l * 3 + 1] = ptr[l * 3 + 1], t
                upd(l)
                t, l = l, v
            else:
                v, ptr[r * 3] = ptr[r * 3], t
                upd(r)
                t, r = r, v
        return t

    def _iter_(self, t):
        """木 t の値を先頭から順に生成する(遅延も解消する)。O(部分木サイズ)"""
        ptr, val = self.ptr, self.val
        push = self._push_
        stack = []
        while stack or t:
            if t:
                push(t)
                stack.append(t)
                t = ptr[t * 3]
            else:
                t = stack.pop()
                yield val[t]
                t = ptr[t * 3 + 1]

    def build(self, arr):
        """
        列 arr で初期化する(既存の内容は破棄)

        計算量: O(N)
        """
        n = len(arr)
        rand = self.rand
        val = self.val = [self.e] + list(arr)
        pri = self.pri = [1 << 32] + [next(rand) for _ in range(n)]
        ptr = self.ptr = [0] * (n * 3 + 3)
        cnt = self.cnt = [0] + [1] * n
        acc = self.acc = [self.e] * (n + 1)
        self.laz = [self.id] * (n + 1)
        op = self.op
        # 優先度について Cartesian tree を構築する
        par = [0] * (n + 1)
        for i in range(2, n + 1):
            p = i - 1
            l = 0
            while p and pri[i] > pri[p]:
                pp = par[p]
                if l:
                    par[l] = p
                par[p] = i
                l, p = p, pp
            par[i] = p
        self.root = 0
        for i in range(1, n + 1):
            p = par[i]
            if not p:
                self.root = i
            elif i < p:
                ptr[p * 3] = i
            else:
                ptr[p * 3 + 1] = i
        if not n:
            return
        stack = [self.root]
        order = []
        while stack:
            v = stack.pop()
            order.append(v)
            l, r = ptr[v * 3], ptr[v * 3 + 1]
            if l:
                stack.append(l)
            if r:
                stack.append(r)
        for v in reversed(order):
            l, r = ptr[v * 3], ptr[v * 3 + 1]
            cnt[v] = cnt[l] + cnt[r] + 1
            acc[v] = op(op(acc[l], val[v]), acc[r])

    def size(self):
        """要素数を返す。O(1)"""
        return self.cnt[self.root]

    def insert(self, pos, x):
        """
        位置 pos(0-index)に x を挿入する(挿入後 x が pos 番目になる)

        計算量: 期待 O(log N)
        """
        l, r = self._split_(self.root, pos)
        self.root = self._merge_(self._merge_(l, self._newnode_(x)), r)

    def erase(self, pos):
        """
        位置 pos(0-index)の要素を削除する

        計算量: 期待 O(log N)
        """
        l, r = self._split_(self.root, pos + 1)
        l, _ = self._split_(l, pos)
        self.root = self._merge_(l, r)

    def get(self, pos):
        """
        位置 pos(0-index)の値を返す

        計算量: 期待 O(log N)
        """
        t1, t2 = self._split_(self.root, pos + 1, True)
        t1, t3 = self._split_(t1, pos, True)
        res = self.val[t3]
        self.root = self._merge_(self._merge_(t1, t3), t2)
        return res

    def reverse(self, l, r):
        """
        区間 [l, r) を反転する

        計算量: 期待 O(log N)
        """
        t2, t3 = self._split_(self.root, r)
        t1, t2 = self._split_(t2, l)
        if t2:
            self.ptr[t2 * 3 + 2] ^= 1
        self.root = self._merge_(self._merge_(t1, t2, False, True), t3, True, False)

    def apply(self, l, r, x):
        """
        区間 [l, r) の各要素に作用素 x を作用させる

        計算量: 期待 O(log N)
        """
        t2, t3 = self._split_(self.root, r)
        t1, t2 = self._split_(t2, l)
        if t2:
            self.laz[t2] = self.composition(x, self.laz[t2])
        self.root = self._merge_(self._merge_(t1, t2, False, True), t3, True, False)

    def prod(self, l, r):
        """
        区間 [l, r) の値の積 op(a[l], ..., a[r-1]) を返す(空なら e)

        計算量: 期待 O(log N)
        """
        t2, t3 = self._split_(self.root, r, True)
        t1, t2 = self._split_(t2, l, True)
        res = self.acc[t2]
        self.root = self._merge_(self._merge_(t1, t2), t3)
        return res

    def iter(self):
        """全要素を先頭から順に生成する。O(N)"""
        yield from self._iter_(self.root)


if __name__ == "__main__":
    tr = ImplicitTreap(min, 1 << 60, lambda f, x: x + f, lambda f, g: f + g, 0)
    tr.build([5, 1, 3])
    tr.insert(1, 4)
    tr.apply(0, 2, 10)
    assert tr.prod(0, 2) == 14
    tr.reverse(0, 4)
    assert list(tr.iter()) == [3, 1, 14, 15]
    tr.erase(0)
    assert tr.get(1) == 14 and tr.size() == 3
