"""1 以上 M 以下のオイラーの φ 関数を篩で一括計算する。

使用例:
    ep = euler_phi(10)
    ep.phi[10]  # 4
"""


class euler_phi:
    """φ(1), ..., φ(M) を保持する。"""

    def __init__(self, M=10**6):
        """1 以上 M 以下の φ を O(M log log M) で構築する。

        M: 上限
        """
        phi = list(range(M + 1))
        for x in range(2, M + 1):
            if phi[x] == x:
                for y in range(x, M + 1, x):
                    phi[y] = phi[y] // x * (x - 1)
        self.phi = phi
