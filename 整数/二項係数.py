"""素数 mod 上の二項係数。

前計算 O(N) のあと C(n, r) を O(1) で求める。

使用例:
    C = combination(10**6, 998244353)
    C.cmb(10, 3)  # 120
"""


class combination:
    """素数 p を法とする二項係数テーブル。"""

    def __init__(self, N, p):
        """0! ... N! とその逆元を前計算する。

        N: 前計算する最大次数
        p: 素数法
        計算量: O(N)
        """
        fact = [1] * (N + 1)
        factinv = [1] * (N + 1)
        inv = [0] * (N + 1)
        inv[1] = 1
        for i in range(2, N + 1):
            fact[i] = fact[i - 1] * i % p
            inv[i] = -inv[p % i] * (p // i) % p
            factinv[i] = factinv[i - 1] * inv[i] % p
        self.fact = fact
        self.factinv = factinv
        self.inv = inv
        self.p = p

    def cmb(self, n, r):
        """C(n, r) mod p を返す。範囲外は 0。

        n, r: 整数
        計算量: O(1)
        """
        if r < 0 or n < r:
            return 0
        r = min(r, n - r)
        return self.fact[n] * self.factinv[r] % self.p * self.factinv[n - r] % self.p
