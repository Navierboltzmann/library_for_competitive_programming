"""区間の最大連続部分列和。

各ノードに (区間全体の最大部分列和, 左端からの最大, 右端からの最大, 区間和) を持つ。

使用例:
    ls = LargestSubarray([1, -2, 3, 2, -1])
    ls.query(0, 5)  # 5  (3+2)
    ls.update(1, 10)
"""


class SegTree:
    """非可換モノイドのセグメント木。"""

    def __init__(self, init_val, segfunc, ide_ele):
        """配列 init_val から構築する。

        init_val: 初期配列
        segfunc: 結合的な二項演算
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
        """半開区間 [l, r) の演算結果を返す。

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


class LargestSubarray:
    """区間の最大連続部分列和。"""

    ide_ele = (0, 0, 0, 0)

    def segfunc(self, x, y):
        """左右のノードを結合する。

        x, y: (全体max, 左max, 右max, 区間和)
        """
        return (
            max(x[0], y[0], x[2] + y[1]),
            max(x[1], x[3] + y[1]),
            max(y[2], x[2] + y[3]),
            x[3] + y[3],
        )

    def __init__(self, A):
        """列 A から構築する。

        A: 数値の列
        計算量: O(N)
        """
        self.T = SegTree([(a, a, a, a) for a in A], self.segfunc, self.ide_ele)

    def update(self, k, x):
        """k 番目を x に更新する。

        計算量: O(log N)
        """
        self.T.update(k, (x, x, x, x))

    def query(self, l, r):
        """半開区間 [l, r) の最大連続部分列和を返す。

        空区間を含まない(要素を 1 つ以上取る)場合は、値が全て負なら最小の 1 要素になる。
        計算量: O(log N)
        """
        return self.T.query(l, r)[0]
