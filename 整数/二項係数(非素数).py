"""任意法での二項係数 (Lucas + Garner)。

mod が素数でなくても C(n, k) mod m を計算できる。

使用例:
    bc = BinomialCoefficient(12)
    bc.choose(10, 3)  # 120 % 12 = 0
"""


def modinv(a, MOD):
    """a の mod MOD における逆元。"""
    r0, r1, s0, s1 = a, MOD, 1, 0
    while r1:
        r0, r1, s0, s1 = r1, r0 % r1, s1, s0 - r0 // r1 * s1
    return s0 % MOD


def factorize(N):
    """N を (p, p^e) のリストに分解する。"""
    factorization = []
    i = 2
    while i * i <= N:
        if N % i == 0:
            c = 0
            pe = 1
            while N % i == 0:
                N //= i
                pe *= i
                c += 1
            factorization.append((i, pe))
        i += 1
    if N != 1:
        factorization.append((N, N))
    return factorization


class BinomialCoefficient:
    """任意法 m の二項係数。"""

    def __init__(self, m):
        """法 m で前計算する。

        m: 法
        """
        self.MOD = m
        self.factorization = factorize(m)
        self.facs = []
        self.invs = []
        self.coeffs = []
        self.pows = []
        for p, pe in self.factorization:
            fac = [1] * pe
            for i in range(1, pe):
                fac[i] = fac[i - 1] * (i if i % p else 1) % pe
            inv = [1] * pe
            inv[-1] = fac[-1]
            for i in range(pe - 1, 0, -1):
                inv[i - 1] = inv[i] * (i if i % p else 1) % pe
            self.facs.append(fac)
            self.invs.append(inv)
            c = modinv(m // pe, pe)
            self.coeffs.append(m // pe * c % m)
            powp = [1]
            while powp[-1] * p != pe:
                powp.append(powp[-1] * p)
            self.pows.append(powp)

    def choose(self, n, k):
        """C(n, k) mod m を返す。

        計算量: O(log n * ω(m))
        """
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1 % self.MOD
        res = 0
        for i, (p, pe) in enumerate(self.factorization):
            res += self._choose_pe(n, k, p, pe, self.facs[i], self.invs[i], self.pows[i]) * self.coeffs[i]
            res %= self.MOD
        return res

    def _E(self, n, k, r, p):
        """p 進付値 v_p(C(n, k)) を計算する途中値。"""
        res = 0
        while n:
            n //= p
            k //= p
            r //= p
            res += n - k - r
        return res

    def _choose_pe(self, n, k, p, pe, fac, inv, powp):
        """C(n, k) mod p^e を返す。"""
        r = n - k
        e0 = self._E(n, k, r, p)
        if e0 >= len(powp):
            return 0
        res = powp[e0]
        if (p != 2 or pe == 4) and self._E(n // (pe // p), k // (pe // p), r // (pe // p), p) % 2:
            res = pe - res
        while n:
            res = res * fac[n % pe] % pe * inv[k % pe] % pe * inv[r % pe] % pe
            n //= p
            k //= p
            r //= p
        return res
