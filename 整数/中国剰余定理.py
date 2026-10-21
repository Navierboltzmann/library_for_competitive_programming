"""拡張ユークリッドと中国剰余定理。

法が互いに素でなくてもよい CRT と、互いに素なときの Garner。

使用例:
    ext_gcd(21, 15)          # (3, -2, 3)  21*(-2)+15*3=3
    crt([2, 3, 2], [3, 5, 7])  # (23, 105)
    garner([2, 3, 2], [3, 5, 7])  # 23
"""


def ext_gcd(a, b):
    """a x + b y = gcd(a, b) となる (g, x, y) を返す。

    計算量: O(log min(|a|, |b|))
    """
    x, y, x1, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x, x1 = x1, x - q * x1
        y, y1 = y1, y - q * y1
    return a, x, y


def modinv(a, m):
    """a の mod m における逆元。存在しなければ ValueError。"""
    g, x, _ = ext_gcd(a % m, m)
    if g != 1:
        raise ValueError("inverse does not exist")
    return x % m


def crt(rems, mods):
    """x ≡ rems[i] (mod mods[i]) を解く。

    法は互いに素でなくてよい。
    戻り値は (r, m)。解は x ≡ r (mod m)、m = lcm。
    解なしなら (0, 0)。
    計算量: O(n log M)
    """
    r0, m0 = 0, 1
    for r1, m1 in zip(rems, mods):
        if m1 <= 0:
            raise ValueError("modulus must be positive")
        r1 %= m1
        if m0 < m1:
            r0, r1 = r1, r0
            m0, m1 = m1, m0
        if m0 % m1 == 0:
            if r0 % m1 != r1:
                return 0, 0
            continue
        g, p, _ = ext_gcd(m0, m1)
        if (r1 - r0) % g:
            return 0, 0
        u1 = m1 // g
        x = (r1 - r0) // g % u1 * p % u1
        r0 += x * m0
        m0 *= u1
        r0 %= m0
    return r0, m0


def garner(rems, mods, mod=None):
    """互いに素な法に対する Garner。0 ≤ x < ∏ mods を返す。

    mod を指定すると x mod その値を返す(途中で法を落とせる)。
    計算量: O(n^2)
    """
    n = len(mods)
    invs = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i):
            invs[i][j] = modinv(mods[j], mods[i])
    xs = [0] * n
    for i in range(n):
        v = rems[i] % mods[i]
        for j in range(i):
            v = (v - xs[j]) * invs[i][j] % mods[i]
        xs[i] = v
    ans, mul = 0, 1
    M = 1
    for i in range(n):
        ans += xs[i] * mul
        mul *= mods[i]
        M *= mods[i]
        if mod is not None:
            ans %= mod
            mul %= mod
    if mod is None:
        return ans % M
    return ans % mod
