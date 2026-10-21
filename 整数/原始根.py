"""原始根・乗法的位数。

原始根が存在するのは n = 1, 2, 4, p^k, 2 p^k (p は奇素数) のときだけ。

使用例:
    primitive_root(7)            # 3
    multiplicative_order(3, 7)   # 6
    is_primitive_root(3, 7)      # True
"""
from math import gcd


def _factorize(n):
    """試し割りで {素数: 指数}。"""
    res = {}
    if n <= 1:
        return res
    while n % 2 == 0:
        res[2] = res.get(2, 0) + 1
        n //= 2
    p = 3
    while p * p <= n:
        while n % p == 0:
            res[p] = res.get(p, 0) + 1
            n //= p
        p += 2
    if n > 1:
        res[n] = res.get(n, 0) + 1
    return res


def totient(n):
    """Euler の φ(n)。"""
    ans = n
    for p in _factorize(n):
        ans = ans // p * (p - 1)
    return ans


def carmichael(n):
    """Carmichael の λ(n)。"""
    if n <= 0:
        raise ValueError("n must be positive")
    fac = _factorize(n)
    ans = 1
    for p, e in fac.items():
        if p == 2 and e >= 3:
            pe = 1 << (e - 2)
        else:
            pe = (p - 1) * p ** (e - 1)
        ans = ans // gcd(ans, pe) * pe
    return ans


def multiplicative_order(a, n):
    """a^e ≡ 1 (mod n) となる最小の e > 0。gcd(a, n)=1 が必要。"""
    if n == 1:
        return 1
    a %= n
    if gcd(a, n) != 1:
        raise ValueError("a and n must be coprime")
    e = totient(n)
    for p in _factorize(e):
        while e % p == 0 and pow(a, e // p, n) == 1:
            e //= p
    return e


def is_primitive_root(g, n):
    """g が mod n の原始根か。"""
    if n == 1:
        return g % 1 == 0
    if gcd(g, n) != 1:
        return False
    phi = totient(n)
    return multiplicative_order(g, n) == phi


def has_primitive_root(n):
    """mod n に原始根が存在するか。"""
    if n <= 0:
        return False
    if n == 1 or n == 2 or n == 4:
        return True
    m = n >> 1 if n % 2 == 0 else n
    if n % 4 == 0 or m % 2 == 0:
        return False
    return len(_factorize(m)) == 1


def primitive_root(n):
    """mod n の最小の原始根。存在しなければ None。

    計算量: 素因数分解 + 小さい g を試す
    """
    if n <= 0:
        raise ValueError("n must be positive")
    if n == 1:
        return 0
    if n == 2 or n == 4:
        return n - 1
    if not has_primitive_root(n):
        return None
    phi = totient(n)
    ps = list(_factorize(phi))
    for g in range(2, n):
        if gcd(g, n) != 1:
            continue
        if all(pow(g, phi // p, n) != 1 for p in ps):
            return g
    return None


def discrete_sqrt_exists(a, p):
    """a が F_p の平方剰余か (p は奇素数または 2)。Euler 判定。"""
    a %= p
    if a == 0:
        return True
    if p == 2:
        return True
    return pow(a, (p - 1) // 2, p) == 1
