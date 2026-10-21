"""Li Chao Tree による Convex Hull Trick。

取りうる x の集合を先に渡し、直線 y = ax+b の追加と点 x での最小値を O(log N) で処理する。

使用例:
    cht = CHT([0, 1, 2, 5])
    cht.add_line(-1, 10)
    cht.add_line(0, 3)
    cht.query(2)  # min(-1*2+10, 0*2+3) = 3
"""


class CHT:
    """Li Chao Tree (区間最小)。"""

    def __init__(self, X):
        """取りうる x のリスト X で構築する。

        X: クエリで使う x の集合
        計算量: O(N)
        """
        X = sorted(set(X))
        self.inf = 10 ** 18
        self.n = 1 << (len(X) - 1).bit_length()
        self.X = X + [self.inf] * (self.n - len(X))
        self.D = {a: i for i, a in enumerate(X)}
        lmr = [(0, 0, 0)] * self.n + [(x, x, x) for x in self.X]
        for i in range(self.n - 1, 0, -1):
            lmr[i] = (lmr[i * 2][0], lmr[i * 2][2], lmr[i * 2 ^ 1][2])
        self.lmr = lmr
        self.F = [None] * (self.n * 2)

    def calc(self, f, x):
        """直線 f=(a, b) の x での値。"""
        a, b = f
        return a * x + b

    def _update(self, i, f):
        """ノード i 以下に直線 f を挿入する。"""
        lmr = self.lmr
        F = self.F
        calc = self.calc
        while True:
            l, m, r = lmr[i]
            fi = F[i]
            if fi is None:
                F[i] = f
                return
            cl = calc(fi, l) > calc(f, l)
            cr = calc(fi, r) > calc(f, r)
            cm = calc(fi, m) > calc(f, m)
            if cl and cr:
                F[i] = f
                return
            if not cl and not cr:
                return
            if cm:
                F[i], f = f, fi
                cl = not cl
            i = i * 2 if cl else i * 2 + 1

    def query(self, x):
        """点 x での最小値を返す。x は構築時の X に含まれていること。

        計算量: O(log N)
        """
        i = self.D[x] + self.n
        mi = self.inf
        F = self.F
        calc = self.calc
        while i > 0:
            if F[i]:
                v = calc(F[i], x)
                if v < mi:
                    mi = v
            i >>= 1
        return mi

    def add_line(self, a, b):
        """直線 y = ax + b を追加する。

        計算量: O(log N)
        """
        self._update(1, (a, b))
