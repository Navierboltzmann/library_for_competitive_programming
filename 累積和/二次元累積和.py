"""2次元累積和。

使用例:
    cs = cs_2d([[1, 2], [3, 4]])
    cs.query(0, 2, 0, 2)  # 10
"""


class cs_2d:
    """グリッドの矩形和を O(1) で取る。"""

    def __init__(self, x):
        """2次元配列 x から累積和を構築する。

        x: n×m の 2 次元リスト
        計算量: O(nm)
        """
        n = len(x)
        m = len(x[0])
        self.n = n
        self.m = m
        S = [0] * ((n + 1) * (m + 1))
        for i in range(n):
            row = x[i]
            ri = m * (i + 1)
            pi = m * i
            for j in range(m):
                S[ri + j + 1] = S[ri + j] + S[pi + j + 1] - S[pi + j] + row[j]
        self.S = S

    def query(self, ix, jx, iy, jy):
        """半開区間 [ix, jx) × [iy, jy) の和を返す。

        ix, jx: 行の半開区間
        iy, jy: 列の半開区間
        計算量: O(1)
        """
        m = self.m
        S = self.S
        return S[m * jx + jy] - S[m * jx + iy] - S[m * ix + jy] + S[m * ix + iy]
