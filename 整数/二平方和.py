"""二平方和 x^2 + y^2 = n と Cornacchia (x^2 + d y^2 = n)。

素数 p = 2 または p ≡ 1 (mod 4) のときだけ p = x^2 + y^2。
一般の n は素因数分解して合成する。

使用例:
    sum_of_two_squares(50)     # (5, 5) など
    cornacchia(1, 13)          # (3, 2)   3^2+2^2=13
    identities_two_squares(1, 2, 3, 4)  # (11, 2) など
"""
from math import isqrt


def _square_root_mod_p(n, p):
    """Tonelli-Shanks。平方非剰余なら None。p は奇素数。"""
    n %= p
    if n == 0:
        return 0
    if pow(n, (p - 1) >> 1, p) != 1:
        return None
    q = p - 1
    s = 0
    while q & 1 == 0:
        q >>= 1
        s += 1
    if s == 1:
        return pow(n, (p + 1) >> 2, p)
    z = 2
    while pow(z, (p - 1) >> 1, p) != p - 1:
        z += 1
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) >> 1, p)
    while t != 1:
        t2i, i = t, 1
        while i < m:
            t2i = t2i * t2i % p
            if t2i == 1:
                break
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        r = r * b % p
        c = b * b % p
        t = t * c % p
        m = i
    return r


def cornacchia(d, n):
    """x^2 + d y^2 = n の非負解 (x, y) を 1 つ返す。無ければ None。

    gcd(x, y) = 1 のもの(原始解)を優先して探す。
    n > 0, d > 0。
    計算量: Tonelli-Shanks + ユークリッド
    """
    if n <= 0 or d <= 0:
        return None
    if n == 1:
        return 1, 0
    # r^2 ≡ -d (mod n) が必要。n が素数ならこれで十分
    # 一般の n ではまず n が素数のときを扱い、合成数は呼び出し側で合成する
    if n % 2 == 0:
        # 偶は直接試す範囲が狭いので √n まで y を見るのは重い。
        # 小さい n だけ落とす
        pass
    # -d が mod n で平方剰余か。n が奇素数のとき
    r = _square_root_mod_p(-d % n, n) if n % 2 else None
    if r is None:
        # n が合成、または非剰余。√(n/d) まで走査(n が小さいとき用)
        ymax = isqrt(n // d)
        for y in range(ymax + 1):
            x = isqrt(n - d * y * y)
            if x * x + d * y * y == n:
                return x, y
        return None
    # ユークリッド: n, r, r^2, ... が n を下回るまで
    a, b = n, r % n
    bound = isqrt(n)
    while b * b > n:
        a, b = b, a % b
    if b == 0:
        return None
    rem = n - b * b
    if rem % d:
        return None
    y2 = rem // d
    y = isqrt(y2)
    if y * y == y2:
        return b, y
    return None


def identities_two_squares(a, b, c, d):
    """(a^2+b^2)(c^2+d^2) の二平方和の 2 通り。"""
    return (abs(a * c - b * d), abs(a * d + b * c)), (abs(a * c + b * d), abs(a * d - b * c))


def _factorize(n):
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
        res[n] = 1
    return res


def sum_of_two_squares(n):
    """n = x^2 + y^2 の非負解 (x, y) (x >= y) を 1 つ。無ければ None。

    p ≡ 3 (mod 4) の指数が奇数なら不能。
    """
    if n < 0:
        return None
    if n == 0:
        return 0, 0
    fac = _factorize(n)
    x, y = 1, 0
    for p, e in fac.items():
        if p == 2:
            for _ in range(e):
                (x, y), _ = identities_two_squares(x, y, 1, 1)
            continue
        if p & 3 == 3:
            if e & 1:
                return None
            p_e2 = p ** (e // 2)
            x *= p_e2
            y *= p_e2
            continue
        one = cornacchia(1, p)
        if one is None:
            return None
        px, py = one
        for _ in range(e):
            (x, y), _ = identities_two_squares(x, y, px, py)
    x, y = abs(x), abs(y)
    if x < y:
        x, y = y, x
    return x, y


def all_sum_of_two_squares(n):
    """n = x^2 + y^2, x >= y >= 0 のすべての解。"""
    if n < 0:
        return []
    res = []
    lim = isqrt(n)
    for y in range(lim + 1):
        x2 = n - y * y
        x = isqrt(x2)
        if x * x == x2 and x >= y:
            res.append((x, y))
        if y * y > n // 2:
            break
    return res
