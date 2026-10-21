"""
区間更新・区間取得の遅延セグメントツリー

区間 [l, r) を値 x で一様に上書きする更新と、区間 [l, r) の segfunc による畳み込みを
それぞれ O(log N) で処理する。
ノードの値を x で直接上書きするので、segfunc は f(x, x) = x を満たすもの
(min, max, gcd, and, or など) に限る。区間和には使えない。
また None は「遅延なし」の印として使うので、None への更新はできない。

使用例:
    seg = LazySegTree_RUQ([5, 3, 8, 6], segfunc, ide_ele)   # 区間最小
    seg.update(1, 3, 7)     # -> [5, 7, 7, 6]
    seg.query(0, 4)         # -> 5
"""


def segfunc(x, y):
    """区間最小 (デフォルトの演算)"""
    return x if x < y else y


ide_ele = 1 << 60


class LazySegTree_RUQ:
    """
    区間更新・区間取得の遅延セグメントツリー

    segfunc: 演算 (f(x, x) = x を満たすこと)
    ide_ele: segfunc の単位元
    """

    def __init__(self, init_val, segfunc, ide_ele):
        """
        初期化 O(N)

        init_val: 配列の初期値
        """
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
        self.lazy = [None] * (2 * num)

    def gindex(self, l, r):
        """
        区間 [l, r) の処理で伝播・再計算が必要なノードを下から順に列挙する

        l: index(0-index)
        r: index(0-index)
        """
        l += self.num
        r += self.num
        lm = l >> (l & -l).bit_length()
        rm = r >> (r & -r).bit_length()
        while r > l:
            if l <= lm:
                yield l
            if r <= rm:
                yield r
            r >>= 1
            l >>= 1
        while l:
            yield l
            l >>= 1

    def propagates(self, *ids):
        """ids (下から順) のノードの遅延値を上から順に子へ伝播する"""
        lazy = self.lazy
        tree = self.tree
        for i in reversed(ids):
            v = lazy[i]
            if v is None:
                continue
            lazy[i] = None
            lazy[2 * i] = v
            lazy[2 * i + 1] = v
            tree[2 * i] = v
            tree[2 * i + 1] = v

    def update(self, l, r, x):
        """
        区間 [l, r) の値を x に更新 O(log N)

        l: index(0-index)
        r: index(0-index)
        x: update value
        """
        ids = list(self.gindex(l, r))
        self.propagates(*ids)
        tree = self.tree
        lazy = self.lazy
        l += self.num
        r += self.num
        while l < r:
            if l & 1:
                lazy[l] = x
                tree[l] = x
                l += 1
            if r & 1:
                lazy[r - 1] = x
                tree[r - 1] = x
            r >>= 1
            l >>= 1
        f = self.segfunc
        for i in ids:
            tree[i] = f(tree[2 * i], tree[2 * i + 1])

    def query(self, l, r):
        """
        区間 [l, r) の segfunc による畳み込みを返す O(log N)

        l: index(0-index)
        r: index(0-index)
        """
        self.propagates(*self.gindex(l, r))
        tree = self.tree
        f = self.segfunc
        res = self.ide_ele
        l += self.num
        r += self.num
        while l < r:
            if l & 1:
                res = f(res, tree[l])
                l += 1
            if r & 1:
                res = f(res, tree[r - 1])
            l >>= 1
            r >>= 1
        return res


if __name__ == "__main__":
    seg = LazySegTree_RUQ([5, 3, 8, 6], segfunc, ide_ele)
    seg.update(1, 3, 7)
    print(seg.query(0, 4))  # 5
