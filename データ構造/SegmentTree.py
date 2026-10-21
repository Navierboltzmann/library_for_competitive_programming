"""可換モノイド上のセグメント木。

使用例:
    def segfunc(x, y):
        return min(x, y)

    ide_ele = 10**18
    st = SegTree([3, 1, 4, 1, 5], segfunc, ide_ele)
    st.query(1, 4)  # 1
    st.update(2, 0)
"""


def segfunc(x, y):
    """区間の二項演算(デフォルトは min)。"""
    return min(x, y)


ide_ele = 10**18


class SegTree:
    """可換モノイドのセグメント木。"""

    def __init__(self, init_val, segfunc, ide_ele):
        """配列 init_val から構築する。

        init_val: 初期配列
        segfunc: 結合的かつ可換な二項演算
        ide_ele: 単位元
        計算量: O(N)
        """
        n = len(init_val)
        self.segfunc = segfunc
        self.ide_ele = ide_ele
        self.num = 1 << (n - 1).bit_length()
        tree = [ide_ele] * (2 * self.num)
        tree[self.num:self.num + n] = init_val
        for i in range(self.num - 1, 0, -1):
            tree[i] = segfunc(tree[2 * i], tree[2 * i + 1])
        self.tree = tree

    def update(self, k, x):
        """k 番目の値を x に更新する。

        k: 0-index
        x: 新しい値
        計算量: O(log N)
        """
        k += self.num
        tree = self.tree
        op = self.segfunc
        tree[k] = x
        while k > 1:
            tree[k >> 1] = op(tree[k], tree[k ^ 1])
            k >>= 1

    def query(self, l, r):
        """半開区間 [l, r) の演算結果を返す。

        l, r: 0-index
        計算量: O(log N)
        """
        op = self.segfunc
        tree = self.tree
        res = self.ide_ele
        l += self.num
        r += self.num
        while l < r:
            if l & 1:
                res = op(res, tree[l])
                l += 1
            if r & 1:
                res = op(res, tree[r - 1])
            l >>= 1
            r >>= 1
        return res
