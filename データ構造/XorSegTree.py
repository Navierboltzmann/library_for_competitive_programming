"""XOR セグメント木。

区間 [L, R) の各 i について A[i xor x] を f で畳み込む。
N は 2 の冪である必要がある。

使用例:
    st = XorSegmentTree([1, 2, 3, 4], lambda a, b: a + b, 0)
    st.range_fold(0, 4, 1)  # A[i xor 1] の総和
"""


class XorSegmentTree:
    """XOR シフト付き区間畳み込み。"""

    def __init__(self, A, f, E):
        """配列 A から構築する。

        A: 長さ 2^k の列
        f: 結合的な二項演算
        E: 単位元
        計算量: O(N log N)
        """
        self.N = len(A)
        ST = [[] for _ in range(self.N * 2 - 1)]
        self.f = f
        self.E = E
        for i in range(self.N):
            ST[self.N - 1 + i].append(A[i])
        for i in range(self.N - 2, -1, -1):
            left = ST[i * 2 + 1]
            right = ST[i * 2 + 2]
            cnt = len(left)
            node = ST[i]
            node_app = node.append
            for j in range(cnt):
                node_app(f(left[j], right[j]))
            for j in range(cnt):
                node_app(f(right[j], left[j]))
        self.ST = ST

    def _range_fold(self, L, R, x, i, l, r):
        """[L, R) をノード i が表す [l, r) で畳み込む。"""
        if r <= L or R <= l:
            return self.E
        if L <= l and r <= R:
            return self.ST[i][x]
        p = (r - l) // 2
        m = (l + r) // 2
        f = self.f
        if (x & p) == 0:
            return f(
                self._range_fold(L, R, x, i * 2 + 1, l, m),
                self._range_fold(L, R, x, i * 2 + 2, m, r),
            )
        resL = self.E
        if R >= m:
            resL = self._range_fold(max(L, m) - p, R - p, x ^ p, i * 2 + 1, l, m)
        resR = self.E
        if L < m:
            resR = self._range_fold(L + p, min(R, m) + p, x ^ p, i * 2 + 2, m, r)
        return f(resR, resL)

    def range_fold(self, L, R, x):
        """Σ_{i in [L, R)} A[i xor x] を f で畳み込む。

        L, R: 半開区間
        x: XOR する値
        計算量: O(log^2 N)
        """
        return self._range_fold(L, R, x, 0, 0, self.N)
