"""非可換モノイド上のセグメント木。

左右の畳み込み順を保つので、行列積など非可換な演算でも使える。

使用例:
    st = SegTree([1, 2, 3], lambda x, y: x + y, 0)
    st.query(0, 3)  # 6
"""


class SegTree:
    """非可換モノイドのセグメント木。"""

    def __init__(self, init_val, segfunc, ide_ele):
        """配列 init_val から構築する。

        init_val: 初期配列
        segfunc: 結合的な二項演算(非可換可)
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
            tk = k >> 1
            tree[tk] = op(tree[tk << 1], tree[(tk << 1) + 1])
            k >>= 1

    def query(self, l, r):
        """半開区間 [l, r) の演算結果を返す(左右の順を保持)。

        l, r: 0-index
        計算量: O(log N)
        """
        op = self.segfunc
        tree = self.tree
        res_l = self.ide_ele
        res_r = self.ide_ele
        l += self.num
        r += self.num
        while l < r:
            if l & 1:
                res_l = op(res_l, tree[l])
                l += 1
            if r & 1:
                res_r = op(tree[r - 1], res_r)
            l >>= 1
            r >>= 1
        return op(res_l, res_r)
