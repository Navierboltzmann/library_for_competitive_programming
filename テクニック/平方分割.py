"""区間和の平方分割。

使用例:
    sd = SqrtDecomposition([1, 2, 3, 4, 5])
    sd.query(1, 4)  # 9
"""


class SqrtDecomposition:
    """区間和を O(√N) で取る。点更新は未実装。"""

    def __init__(self, a):
        """列 a からバケット和を構築する。

        a: 数値の列
        計算量: O(N)
        """
        self.N = len(a)
        n = int(self.N ** 0.5)
        if n * n < self.N:
            n += 1
        if n == 0:
            n = 1
        self.n = n
        pad = n * n - self.N
        self.a = a + [0] * pad
        B = [0] * n
        aa = self.a
        for i in range(self.N):
            B[i // n] += aa[i]
        self.B = B

    def query(self, l, r):
        """半開区間 [l, r) の和を返す。

        l, r: 0-index
        計算量: O(√N)
        """
        n = self.n
        a = self.a
        B = self.B
        L = (l + n - 1) // n
        R = (r + n - 1) // n
        if L == R:
            s = 0
            for i in range(l, r):
                s += a[i]
            return s
        s = 0
        for i in range(l, n * L):
            s += a[i]
        for i in range(L, R - 1):
            s += B[i]
        for i in range(n * (R - 1), r):
            s += a[i]
        return s
