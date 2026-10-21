"""mod 付き行列クラス。

加算・減算・積・スカラー倍・累乗をサポートする。

使用例:
    A = Matrix(2, 2, [[1, 1], [1, 0]])
    print(A ** 10)
"""


class Matrix:
    """法 MOD 上の行列。"""

    def __init__(self, n, m, mat=None, mod=998244353):
        """n×m 行列を作る。mat があればコピーする。

        n, m: 行数・列数
        mat: 初期値
        mod: 法
        """
        self.n = n
        self.m = m
        self.MOD = mod
        self.mat = [[0] * self.m for _ in range(self.n)]
        if mat:
            assert len(mat) == n and len(mat[0]) == m
            for i in range(self.n):
                self.mat[i] = mat[i][:]

    def is_square(self):
        """正方行列か。"""
        return self.n == self.m

    def __getitem__(self, key):
        """行の取得。"""
        if isinstance(key, slice):
            return self.mat[key]
        assert key >= 0
        return self.mat[key]

    def id(n):
        """n 次単位行列。"""
        res = Matrix(n, n)
        for i in range(n):
            res[i][i] = 1
        return res

    def __len__(self):
        """行数。"""
        return len(self.mat)

    def __str__(self):
        """空白区切りの文字列。"""
        return '\n'.join(' '.join(map(str, self[i])) for i in range(self.n))

    def times(self, k):
        """スカラー倍 kA。"""
        res = [[0] * self.m for _ in range(self.n)]
        MOD = self.MOD
        for i in range(self.n):
            res_i, self_i = res[i], self[i]
            for j in range(self.m):
                res_i[j] = k * self_i[j] % MOD
        return Matrix(self.n, self.m, res, self.MOD)

    def __pos__(self):
        """+A。"""
        return self

    def __neg__(self):
        """-A。"""
        return self.times(-1)

    def __add__(self, other):
        """行列和。"""
        assert self.n == other.n and self.m == other.m
        res = [[0] * self.m for _ in range(self.n)]
        MOD = self.MOD
        for i in range(self.n):
            res_i, self_i, other_i = res[i], self[i], other[i]
            for j in range(self.m):
                res_i[j] = (self_i[j] + other_i[j]) % MOD
        return Matrix(self.n, self.m, res, self.MOD)

    def __sub__(self, other):
        """行列差。"""
        assert self.n == other.n and self.m == other.m
        res = [[0] * self.m for _ in range(self.n)]
        MOD = self.MOD
        for i in range(self.n):
            res_i, self_i, other_i = res[i], self[i], other[i]
            for j in range(self.m):
                res_i[j] = (self_i[j] - other_i[j]) % MOD
        return Matrix(self.n, self.m, res, self.MOD)

    def __mul__(self, other):
        """行列積またはスカラー倍。"""
        if other.__class__ == Matrix:
            assert self.m == other.n
            res = [[0] * other.m for _ in range(self.n)]
            MOD = self.MOD
            for i in range(self.n):
                res_i, self_i = res[i], self[i]
                for k in range(self.m):
                    self_ik = self_i[k]
                    if not self_ik:
                        continue
                    other_k = other[k]
                    for j in range(other.m):
                        res_i[j] = (res_i[j] + self_ik * other_k[j]) % MOD
            return Matrix(self.n, other.m, res, self.MOD)
        return self.times(other)

    def __rmul__(self, other):
        """左からのスカラー倍。"""
        return self.times(other)

    def __pow__(self, k):
        """正方行列の k 乗。

        計算量: O(N^3 log k)
        """
        assert self.is_square()
        tmp = Matrix(self.n, self.n, self.mat, self.MOD)
        res = Matrix.id(self.n)
        res.MOD = self.MOD
        while k:
            if k & 1:
                res *= tmp
            tmp *= tmp
            k >>= 1
        return res
