"""Sparse Table。冪等な区間演算を O(1) で取る。

デフォルトは min。構築 O(N log N)、更新不可。

使用例:
    st = SparseTable([3, 1, 4, 1, 5])
    st.query(1, 4)  # min(1,4,1) = 1
    st.query(0, 5)  # 1
"""


class SparseTable:
    """静的区間クエリ。"""

    def __init__(self, A, op=min):
        """列 A を構築する。

        A: 列
        op: 結合的かつ冪等な二項演算 (min / max / gcd など)
        計算量: O(N log N)
        """
        self.op = op
        self.n = n = len(A)
        if n == 0:
            self.table = []
            self.log = []
            return
        log = [0] * (n + 1)
        for i in range(2, n + 1):
            log[i] = log[i >> 1] + 1
        k = log[n] + 1
        table = [A[:]]
        for j in range(1, k):
            prev = table[j - 1]
            span = 1 << (j - 1)
            cur = [op(prev[i], prev[i + span]) for i in range(n - (1 << j) + 1)]
            table.append(cur)
        self.table = table
        self.log = log

    def query(self, l, r):
        """[l, r) の演算結果。空区間は呼んではいけない。

        計算量: O(1)
        """
        k = self.log[r - l]
        row = self.table[k]
        return self.op(row[l], row[r - (1 << k)])
