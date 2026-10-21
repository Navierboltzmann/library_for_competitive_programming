"""直角三角形領域への加算 (imos)。

左上が直角の直角二等辺三角形に値を足し、build 後に各マスを O(1) で取る。
参考: https://imoz.jp/algorithms/imos_method.html

使用例:
    t = triangle_imos(4)
    t.add(1, 0, 0, 3)
    t.build()
    t.query(0, 0)  # 1
"""


class triangle_imos:
    """左上が直角の三角形加算。"""

    def __init__(self, N):
        """N×N グリッドを初期化する。

        N: 一辺の長さ
        """
        self.S = [[0] * (N + 3) for _ in range(N + 3)]
        self.N = N

    def add(self, w, x, y, d):
        """左上 (x, y)、隣辺長 d の直角三角形に w を加算する。

        w: 加える値
        x, y: 左上(0-index)
        d: 隣辺の長さ
        計算量: O(1)
        """
        _x = x + 1
        _y = y + 1
        S = self.S
        S[_x][_y] += w
        S[_x - 1][_y + 1] -= w
        S[_x + d][_y] -= w
        S[_x + d][_y + 1] += w
        S[_x][_y + d + 1] -= w
        S[_x - 1][_y + d + 1] += w

    def build(self):
        """差分を展開して各マスの値を確定する。

        計算量: O(N^2)
        """
        S = self.S
        N = self.N
        for i in range(N + 3):
            for j in range(i):
                S[i - j - 1][j + 1] += S[i - j][j]
        for j in range(N + 3):
            for i in range(N + 2):
                S[i + 1][j] += S[i][j]
        for i in range(N + 3):
            row = S[i]
            for j in range(N + 2):
                row[j + 1] += row[j]

    def query(self, x, y):
        """マス (x, y) の値を返す。build 後に呼ぶ。

        計算量: O(1)
        """
        return self.S[x + 1][y + 1]
