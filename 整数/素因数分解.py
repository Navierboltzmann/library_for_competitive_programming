"""線形篩による素因数分解・約数・メビウス関数。

使用例:
    pf = prime_factorize(10**6)
    pf.factors(12)    # [2, 2, 3]
    pf.divisors(12)   # [1, 2, 3, 4, 6, 12]
    pf.is_prime(13)   # True
    pf.mobius(6)      # 1
"""
from collections import Counter
from itertools import product


class prime_factorize:
    """1 以上 M 以下の篩。"""

    def __init__(self, M=10**6):
        """最小素因数・素数判定・メビウスを O(M log log M) で構築する。

        M: 上限
        """
        sieve = [-1] * (M + 1)
        sieve[1] = 1
        p = [False] * (M + 1)
        mu = [1] * (M + 1)
        for i in range(2, M + 1):
            if sieve[i] == -1:
                p[i] = True
                i2 = i * i
                for j in range(i2, M + 1, i2):
                    mu[j] = 0
                for j in range(i, M + 1, i):
                    sieve[j] = i
                    mu[j] *= -1
        self.sieve = sieve
        self.p = p
        self.mu = mu

    def factors(self, x):
        """x の素因数を小さい順とは限らないリストで返す。

        計算量: O(log x)
        """
        sieve = self.sieve
        tmp = []
        while sieve[x] != x:
            tmp.append(sieve[x])
            x //= sieve[x]
        tmp.append(sieve[x])
        return tmp

    def divisors(self, x):
        """x の約数を昇順で返す。

        計算量: O(d(x))
        """
        C = Counter(self.factors(x))
        tmp = []
        for ps in product(*[[k ** i for i in range(v + 1)] for k, v in C.items()]):
            res = 1
            for pp in ps:
                res *= pp
            tmp.append(res)
        tmp.sort()
        return tmp

    def is_prime(self, x):
        """x が素数か。"""
        return self.p[x]

    def mobius(self, x):
        """μ(x) を返す。"""
        return self.mu[x]
