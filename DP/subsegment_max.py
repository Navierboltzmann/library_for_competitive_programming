"""区間 [l, r) の区間最大値を前計算する。

b[i][j] が区間 [i, j) の「素の値」で、それに区間の入れ子 max を載せる。

使用例:
    n = 3
    b = [[-1] * (n + 1) for _ in range(n)]
    b[0][1] = 1; b[1][2] = 2; b[2][3] = 3
    b[0][2] = 0; b[1][3] = 0; b[0][3] = 4
    ss = ss_max(b)
    ss.query(0, 3)  # 4
"""


class ss_max:
    """区間の入れ子 max を O(N^2) 前計算する。"""

    def __init__(self, b):
        """b[i][j] を区間 [i, j) の素の値として max テーブルを構築する。

        b: n × (n+1) の 2 次元リスト
        計算量: O(N^2)
        """
        n = len(b)
        self.n = n
        B = [[-1] * (n + 1) for _ in range(n)]
        for i in range(n):
            B[i][i + 1] = b[i][i + 1]
        for d in range(2, n + 1):
            for i in range(n + 1 - d):
                B[i][i + d] = max(B[i][i + d - 1], B[i + 1][i + d], b[i][i + d])
        self.B = B

    def query(self, l, r):
        """区間 [l, r) の max を返す。

        計算量: O(1)
        """
        return self.B[l][r]
