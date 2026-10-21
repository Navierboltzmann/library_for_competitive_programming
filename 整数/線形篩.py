"""線形篩で乗法的関数を一括前計算する。

素数・最小素因数・φ・μ・d(n)・σ(n) を O(N) で構築する。

使用例:
    sp = linear_sieve(20)
    sp.primes[:5]     # [2, 3, 5, 7, 11]
    sp.phi[10]        # 4
    sp.mu[6]          # 1
    sp.tau[12]        # 6
    sp.sigma[12]      # 28
    sp.lpf[12]        # 2
"""


class linear_sieve:
    """1 ... N の線形篩。"""

    def __init__(self, N):
        """0..N を O(N) で篩う。

        N: 上限
        """
        self.N = N
        lpf = [0] * (N + 1)
        primes = []
        phi = [0] * (N + 1)
        mu = [0] * (N + 1)
        tau = [0] * (N + 1)
        sigma = [0] * (N + 1)
        # p^e までの積を持つ補助
        tau_e = [0] * (N + 1)
        sigma_pe = [0] * (N + 1)
        phi[1] = 1
        mu[1] = 1
        tau[1] = 1
        sigma[1] = 1
        for i in range(2, N + 1):
            if lpf[i] == 0:
                lpf[i] = i
                primes.append(i)
                phi[i] = i - 1
                mu[i] = -1
                tau[i] = 2
                sigma[i] = 1 + i
                tau_e[i] = 1
                sigma_pe[i] = 1 + i
            for p in primes:
                v = i * p
                if v > N:
                    break
                lpf[v] = p
                if i % p == 0:
                    phi[v] = phi[i] * p
                    mu[v] = 0
                    tau_e[v] = tau_e[i] + 1
                    tau[v] = tau[i] // (tau_e[i] + 1) * (tau_e[v] + 1)
                    sigma_pe[v] = sigma_pe[i] * p + 1
                    sigma[v] = sigma[i] // sigma_pe[i] * sigma_pe[v]
                    break
                phi[v] = phi[i] * (p - 1)
                mu[v] = -mu[i]
                tau_e[v] = 1
                tau[v] = tau[i] * 2
                sigma_pe[v] = 1 + p
                sigma[v] = sigma[i] * (1 + p)
        self.lpf = lpf
        self.primes = primes
        self.phi = phi
        self.mu = mu
        self.tau = tau
        self.sigma = sigma

    def is_prime(self, x):
        """x が素数か。1 <= x <= N。"""
        return x >= 2 and self.lpf[x] == x

    def factorize(self, x):
        """x を {素数: 指数} に分解。O(log x)。"""
        res = {}
        lpf = self.lpf
        while x > 1:
            p = lpf[x]
            res[p] = res.get(p, 0) + 1
            x //= p
        return res

    def mobius_prefix(self):
        """P[i] = Σ_{k=1}^{i} μ(k)。"""
        P = [0] * (self.N + 1)
        mu = self.mu
        for i in range(1, self.N + 1):
            P[i] = P[i - 1] + mu[i]
        return P
