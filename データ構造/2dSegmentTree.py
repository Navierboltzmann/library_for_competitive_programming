"""
SegTree2D: 2次元セグメント木(1点更新・矩形区間取得)

N×M の2次元配列に対し、1点更新 O(log N log M)、矩形 [Lx, Rx)×[Ly, Ry) の集約 O(log N log M)。
mode で演算を選ぶか、func と default(単位元)を直接指定する。
メモリは 4NM 程度使うので N*M が大きすぎる場合は注意。

使用例:
    st = SegTree2D([[1, 2, 3], [4, 5, 6]], mode='sum')
    st.query(0, 2, 1, 3)   # -> 2 + 3 + 5 + 6 = 16
    st.update(1, 2, 10)    # (1, 2) の値を 10 に
    st.leafvalue(1, 2)     # -> 10
    st2 = SegTree2D(grid, func=lambda a, b: max(a, b), default=-1)
"""
import math


def _lcm(x, y):
    """x と y の最小公倍を返す。"""
    return x * y // math.gcd(x, y)


class SegTree2D():
    """
    2次元セグメント木

    dat[i][j]: 行方向ノード i × 列方向ノード j が担当する矩形の集約値
    """

    DEFAULT = {
        'min': 1 << 60,
        'max': -(1 << 60),
        'sum': 0,
        'prd': 1,
        'gcd': 0,
        'lmc': 1,
        'lcm': 1,
        '^': 0,
        '&': (1 << 60) - 1,
        '|': 0,
    }

    FUNC = {
        'min': min,
        'max': max,
        'sum': (lambda x, y: x + y),
        'prd': (lambda x, y: x * y),
        'gcd': math.gcd,
        'lmc': _lcm,
        'lcm': _lcm,
        '^': (lambda x, y: x ^ y),
        '&': (lambda x, y: x & y),
        '|': (lambda x, y: x | y),
    }

    def __init__(self, ls2D, mode='min', func=None, default=None):
        """
        ls2D: 初期値の2次元配列(N×M, N, M >= 1)
        mode: 'min', 'max', 'sum', 'prd'(積), 'gcd', 'lcm'('lmc' も可), '^', '&', '|'
        func, default: 指定すれば任意の(結合的な)関数と単位元で計算する
        計算量: O(NM)
        """
        N = len(ls2D)
        M = len(ls2D[0])
        self.default = self.DEFAULT[mode] if default is None else default
        self.func = self.FUNC[mode] if func is None else func
        self.N = N
        self.M = M
        self.KN = (N - 1).bit_length()
        self.KM = (M - 1).bit_length()
        self.N2 = 1 << self.KN
        self.M2 = 1 << self.KM
        d = self.default
        width = 2 * self.M2
        self.dat = [[d] * width for _ in range(2 * self.N2)]
        M2 = self.M2
        for i in range(N):
            row, src = self.dat[self.N2 + i], ls2D[i]
            for j in range(M):
                row[M2 + j] = src[j]
        self.build()

    def build(self):
        """
        葉(dat[N2 + i][M2 + j])の値から全ノードを計算する

        計算量: O(NM)
        """
        f = self.func
        dat = self.dat
        M2 = self.M2
        lo, hi = M2, M2 + self.M
        for i in range(self.N2 - 1, 0, -1):
            row, a, b = dat[i], dat[i << 1], dat[i << 1 | 1]
            for j in range(lo, hi):
                row[j] = f(a[j], b[j])
        for row in dat:
            for j in range(M2 - 1, 0, -1):
                row[j] = f(row[j << 1], row[j << 1 | 1])

    def leafvalue(self, x, y):
        """
        (x, y) の値を返す(0-index)

        計算量: O(1)
        """
        return self.dat[x + self.N2][y + self.M2]

    def update(self, x, y, value):
        """
        (x, y) の値を value に更新する(0-index)

        計算量: O(log N log M)
        """
        f = self.func
        dat = self.dat
        i = x + self.N2
        j0 = y + self.M2
        row = dat[i]
        row[j0] = value
        j = j0 >> 1
        while j:
            row[j] = f(row[j << 1], row[j << 1 | 1])
            j >>= 1
        i >>= 1
        while i:
            row = dat[i]
            row[j0] = f(dat[i << 1][j0], dat[i << 1 | 1][j0])
            j = j0 >> 1
            while j:
                row[j] = f(row[j << 1], row[j << 1 | 1])
                j >>= 1
            i >>= 1

    def _query_row(self, row, Ly, Ry):
        """1行分(行ノード row)について列 [Ly, Ry)(内部 index)の集約値を返す。O(log M)"""
        f = self.func
        vL = vR = self.default
        while Ly < Ry:
            if Ly & 1:
                vL = f(vL, row[Ly])
                Ly += 1
            if Ry & 1:
                Ry -= 1
                vR = f(row[Ry], vR)
            Ly >>= 1
            Ry >>= 1
        return f(vL, vR)

    def query(self, Lx, Rx, Ly, Ry):
        """
        矩形 [Lx, Rx) × [Ly, Ry) の集約値を返す(0-index, 半開区間)

        空なら default を返す。計算量: O(log N log M)
        """
        f = self.func
        dat = self.dat
        qrow = self._query_row
        Lx += self.N2
        Rx += self.N2
        Ly += self.M2
        Ry += self.M2
        vL = vR = self.default
        while Lx < Rx:
            if Lx & 1:
                vL = f(vL, qrow(dat[Lx], Ly, Ry))
                Lx += 1
            if Rx & 1:
                Rx -= 1
                vR = f(qrow(dat[Rx], Ly, Ry), vR)
            Lx >>= 1
            Rx >>= 1
        return f(vL, vR)


if __name__ == "__main__":
    st = SegTree2D([[1, 2, 3], [4, 5, 6]], mode='sum')
    assert st.query(0, 2, 1, 3) == 16
    st.update(1, 2, 10)
    assert st.leafvalue(1, 2) == 10
    assert st.query(0, 2, 0, 3) == 25
