"""Fenwick Tree (Binary Indexed Tree)。

0-index で点加算・prefix 和・区間和。

使用例:
    ft = FenwickTree(5)
    ft.build([1, 2, 3, 4, 5])
    ft.query(1, 4)  # 9
    ft.add(2, 10)
"""


class FenwickTree:
    """区間和 BIT。"""

    def __init__(self, n):
        """長さ n の 0 初期化 BIT を作る。

        n: サイズ
        """
        self.n = n
        self.data = [0] * n

    def build(self, arr):
        """arr の値で O(N) 構築する。

        arr: 長さ n 以下の列
        """
        data = self.data
        n = self.n
        for i, a in enumerate(arr):
            data[i] = a
        for i in range(1, n + 1):
            j = i + (i & -i)
            if j <= n:
                data[j - 1] += data[i - 1]

    def add(self, p, x):
        """p 番目に x を加算する。

        p: 0-index
        x: 加算値
        計算量: O(log N)
        """
        data = self.data
        n = self.n
        p += 1
        while p <= n:
            data[p - 1] += x
            p += p & -p

    def sum(self, r):
        """半開区間 [0, r) の和を返す。

        r: 0-index の右端
        計算量: O(log N)
        """
        data = self.data
        s = 0
        while r:
            s += data[r - 1]
            r -= r & -r
        return s

    def query(self, l, r):
        """半開区間 [l, r) の和を返す。

        計算量: O(log N)
        """
        return self.sum(r) - self.sum(l)
