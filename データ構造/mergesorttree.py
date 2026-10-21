"""
MergeSortTree: 各ノードに区間内の値をソートして持つセグメント木

区間 [l, r) を O(log N) 個のノードに分解し、各ノードのソート済み配列に対して二分探索などを行う。
query は「区間内で x に最も近い値との差 |x - a| の最小値」を返す実装になっている
(https://yukicoder.me/problems/no/1332)。別の問題に使う場合は _f と query のマージ方法を書き換える。
区間内の x 未満の個数 count_less も用意している。

使用例:
    mst = MergeSortTree([5, 1, 9, 3])
    mst.query(0, 4, 7)       # -> 2  (|7 - 5| or |7 - 9|)
    mst.query(1, 2, 7)       # -> 6  (|7 - 1|)
    mst.count_less(0, 4, 5)  # -> 2  (1, 3)
"""
from bisect import bisect_left, bisect_right


class MergeSortTree:
    """
    各ノードが担当区間の値をソート済みリストで持つセグメント木

    構築: O(N log N), メモリ: O(N log N)
    """

    INF = 1 << 62

    def __init__(self, init_val):
        """
        init_val: 初期配列(list)

        計算量: O(N log N)
        """
        n = len(init_val)
        num = self.num = 1 << (n - 1).bit_length()
        tree = self.tree = [[] for _ in range(2 * num)]
        for i, a in enumerate(init_val):
            tree[num + i] = [a]
        for i in range(num - 1, 0, -1):
            # 2つのソート済み列の連結は Timsort がほぼ線形でマージする
            tree[i] = sorted(tree[2 * i] + tree[2 * i + 1])

    def _f(self, x, A):
        """
        1ノード分の処理: ソート済み A の中で x に最も近い値との差を返す

        x: クエリの値
        A: ノードのソート済みリスト(空でないこと)
        計算量: O(log |A|)
        """
        idx = bisect_right(A, x)
        if idx == 0:
            return A[0] - x
        if idx == len(A):
            return x - A[-1]
        d1 = x - A[idx - 1]
        d2 = A[idx] - x
        return d1 if d1 < d2 else d2

    def query(self, l, r, x):
        """
        区間 [l, r) の値 a について |x - a| の最小値を返す

        l: 左端(0-index, 含む)
        r: 右端(0-index, 含まない)
        x: クエリの値
        戻り値: 最小値(区間が空なら INF)
        計算量: O(log^2 N)
        """
        f = self._f
        tree = self.tree
        res = self.INF
        l += self.num
        r += self.num
        while l < r:
            if l & 1:
                v = f(x, tree[l])
                if v < res:
                    res = v
                l += 1
            if r & 1:
                r -= 1
                v = f(x, tree[r])
                if v < res:
                    res = v
            l >>= 1
            r >>= 1
        return res

    def count_less(self, l, r, x):
        """
        区間 [l, r) の中で x 未満の値の個数を返す

        計算量: O(log^2 N)
        """
        tree = self.tree
        res = 0
        l += self.num
        r += self.num
        while l < r:
            if l & 1:
                res += bisect_left(tree[l], x)
                l += 1
            if r & 1:
                r -= 1
                res += bisect_left(tree[r], x)
            l >>= 1
            r >>= 1
        return res


if __name__ == "__main__":
    mst = MergeSortTree([5, 1, 9, 3])
    assert mst.query(0, 4, 7) == 2
    assert mst.query(1, 2, 7) == 6
    assert mst.count_less(0, 4, 5) == 2
