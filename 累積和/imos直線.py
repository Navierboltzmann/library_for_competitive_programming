"""一次関数を区間加算する imos。

区間 [l, r) に ax+b を足し、最後に全点の値を O(N) で確定する。

使用例:
    im = imos_line(5)
    im.add_val(1, 0, 1, 4)  # [1, 4) に x を加算
    im.query()  # [0, 1, 2, 3, 0]
"""


class imos_line:
    """一次関数の区間加算。"""

    def __init__(self, N):
        """長さ N の列を 0 で初期化する。

        N: 列長
        """
        self.N = N
        self.A = [0] * (N + 1)
        self.B = [0] * (N + 1)
        self.C = [0] * (N + 1)

    def add_val(self, a, b, l, r):
        """半開区間 [l, r) に ax+b を加算する。

        a, b: 一次関数の係数
        l, r: 0-index の半開区間
        計算量: O(1)
        """
        self.B[l] += a * l + b
        self.B[r] -= a * l + b
        self.A[l + 1] += a
        self.A[r] -= a
        self.C[r] -= (r - l - 1) * a

    def query(self):
        """imos を展開し、各位置の値のリストを返す。

        計算量: O(N)
        """
        A = self.A
        B = self.B
        C = self.C
        N = self.N
        for i in range(N):
            B[i + 1] += B[i]
            A[i + 1] += A[i]
        for i in range(N):
            A[i + 1] += A[i] + C[i + 1]
        return [A[i] + B[i] for i in range(N)]
