"""一定間隔ごとの累積和・imos。

x ≡ k (mod Mod) となる位置だけを対象にする。

使用例:
    hc = hopping_cumsum([1, 2, 3, 4, 5], 2)
    hc.query(0, 5, 1)  # 2+4 = 6

    hi = hopping_imos(5, 2)
    hi.add(0, 5, 1, 0)
    hi.build()
    hi.query(2)  # 1
"""


class hopping_cumsum:
    """一定間隔の区間和。"""

    def __init__(self, X, Mod):
        """列 X から間隔 Mod の累積和を構築する。

        X: 数値の列
        Mod: 間隔
        計算量: O(N)
        """
        self.Mod = Mod
        self.X = X
        self.N = len(X)
        T = [0] * Mod + X
        for i in range(self.N):
            T[i + Mod] += T[i]
        self.T = T

    def query(self, l, r, k):
        """[l, r) のうち x ≡ k (mod Mod) の和を返す。

        計算量: O(1)
        """
        L = l + (k - l) % self.Mod
        R = r + (k - r) % self.Mod
        return self.T[R] - self.T[L]


class hopping_imos:
    """一定間隔の区間加算。"""

    def __init__(self, N, Mod):
        """長さ N、間隔 Mod の imos を作る。

        N: 列長
        Mod: 間隔
        """
        self.N = N
        self.Mod = Mod
        self.X = [0] * (N + Mod)

    def add(self, l, r, a, k):
        """[l, r) のうち x ≡ k (mod Mod) に a を加える。

        計算量: O(1)
        """
        L = l + (k - l) % self.Mod
        R = r + (k - r) % self.Mod
        if L < self.N:
            self.X[L] += a
        if R < self.N:
            self.X[R] -= a

    def build(self):
        """差分を展開する。

        計算量: O(N)
        """
        X = self.X
        Mod = self.Mod
        for i in range(self.N):
            X[i + Mod] += X[i]

    def query(self, idx):
        """位置 idx の値を返す。build 後に呼ぶ。"""
        return self.X[idx]
