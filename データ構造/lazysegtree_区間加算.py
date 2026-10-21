"""
区間加算・区間和 (BIT 2本による実装)

区間 [l, r) への一様加算と区間 [l, r) の和を、それぞれ O(log N) で処理する。
mod を指定すると和を mod で割った余りで管理する。

使用例:
    bit = RSQandRAQ(5)
    bit.add(1, 4, 3)    # A = [0, 3, 3, 3, 0]
    bit.add(0, 2, 1)    # A = [1, 4, 3, 3, 0]
    bit.query(1, 3)     # -> 7
"""


class RSQandRAQ:
    """
    区間加算、区間和取得クエリをそれぞれ O(log N) で答える

    add: 区間 [l, r) に val を加える
    query: 区間 [l, r) の和を求める
    l, r は 0-indexed
    """

    def __init__(self, n, mod=None):
        """
        初期化 (全要素 0) O(N)

        n: 要素数
        mod: 法 (None なら mod を取らない)
        """
        self.n = n
        self.bit0 = [0] * (n + 1)
        self.bit1 = [0] * (n + 1)
        self.mod = mod

    def _add(self, bit, i, val):
        """BIT bit の位置 i (0-indexed) に val を加える O(log N)"""
        n = self.n
        mod = self.mod
        i += 1
        if mod is None:
            while i <= n:
                bit[i] += val
                i += i & -i
        else:
            while i <= n:
                bit[i] = (bit[i] + val) % mod
                i += i & -i

    def _get(self, bit, i):
        """BIT bit の先頭 i 項の和を返す O(log N)"""
        s = 0
        mod = self.mod
        if mod is None:
            while i > 0:
                s += bit[i]
                i -= i & -i
        else:
            while i > 0:
                s = (s + bit[i]) % mod
                i -= i & -i
        return s

    def _add2(self, i, v0, v1):
        """位置 i (0-indexed) について bit0 に v0、bit1 に v1 を同時に加える"""
        n = self.n
        mod = self.mod
        bit0 = self.bit0
        bit1 = self.bit1
        i += 1
        if mod is None:
            while i <= n:
                bit0[i] += v0
                bit1[i] += v1
                i += i & -i
        else:
            while i <= n:
                bit0[i] = (bit0[i] + v0) % mod
                bit1[i] = (bit1[i] + v1) % mod
                i += i & -i

    def _prefix(self, i):
        """先頭 i 項 (A[0] + ... + A[i-1]) の和を返す (mod 未適用の場合あり)"""
        bit0 = self.bit0
        bit1 = self.bit1
        s0 = s1 = 0
        j = i
        while j > 0:
            s0 += bit0[j]
            s1 += bit1[j]
            j -= j & -j
        return s0 + i * s1

    def add(self, l, r, val):
        """
        区間 [l, r) に val を加える O(log N)

        l: index(0-index)
        r: index(0-index)
        val: 加える値
        """
        self._add2(l, -val * l, val)
        self._add2(r, val * r, -val)

    def query(self, l, r):
        """
        区間 [l, r) の和を求める O(log N)

        l: index(0-index)
        r: index(0-index)
        """
        res = self._prefix(r) - self._prefix(l)
        if self.mod is None:
            return res
        return res % self.mod


if __name__ == "__main__":
    bit = RSQandRAQ(5)
    bit.add(1, 4, 3)
    bit.add(0, 2, 1)
    print(bit.query(1, 3))  # 7
