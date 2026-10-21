"""二次合同方程式 x^2 ≡ a (mod n)。

F_p は Tonelli-Shanks、p^e は Hensel、一般の n は CRT で合成する。

使用例:
    sqrt_mod(2, 7)          # [3, 4]
    sqrt_mod(2, 13)         # []
    sqrt_mod(8, 16)         # [4, 12] など
"""
from math import gcd


def _tonelli_shanks(n, p):
    """mod 奇素数 p の平方根を 1 つ。非剰余なら None。"""
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


def _sqrt_mod_2e(a, e):
    """x^2 ≡ a (mod 2^e) の 0 <= x < 2^e の解。"""
    a &= (1 << e) - 1
    if e == 0:
        return [0]
    if e == 1:
        return [a & 1]
    if e == 2:
        return [x for x in range(4) if x * x % 4 == a]
    # a ≡ 0, 1 (mod 4) が必要。さらに e>=3 なら a ≡ 1 (mod 8) か 4|a
    sols = _sqrt_mod_2e(a, e - 1)
    pe, pe1 = 1 << (e - 1), 1 << e
    res = set()
    for x in sols:
        for t in (x, x + pe):
            t &= pe1 - 1
            if (t * t - a) % pe1 == 0:
                res.add(t)
    # Hensel がつぶれる (x が偶数) ときは 2 ずつ上げた候補も見る
    if not res:
        for x in sols:
            for k in range(2):
                t = (x + k * pe) & (pe1 - 1)
                if (t * t - a) % pe1 == 0:
                    res.add(t)
    return sorted(res)


def _hensel_odd(a, p, e):
    """x^2 ≡ a (mod p^e)、p 奇素数。"""
    pe_full = p ** e
    a %= pe_full
    a0 = a % p
    r = _tonelli_shanks(a0, p)
    if r is None:
        return []
    if r == 0:
        # p | a。p^e | a なら 0 のみ (e=1) など、付値で場合分け
        v = 0
        aa = a
        pe_full = p ** e
        aa %= pe_full
        if aa == 0:
            step = p ** ((e + 1) // 2)
            return list(range(0, pe_full, step))
        tmp = aa
        while tmp % p == 0:
            tmp //= p
            v += 1
        if v % 2:
            return []
        if v >= e:
            # すべて p^{ceil(e/2)} の倍数
            step = p ** ((e + 1) // 2)
            return list(range(0, pe_full, step))
        # a = p^v a', v even, gcd(a', p)=1
        # x = p^{v/2} x', x'^2 ≡ a' (mod p^{e-v})
        half = v // 2
        a_odd = aa // p ** v
        inner = _hensel_odd(a_odd, p, e - v)
        step = p ** half
        inner_mod = p ** (e - v)
        res = []
        for x in inner:
            for t in range(step):
                res.append((x + t * inner_mod) * step % pe_full)
        return sorted(set(res))
    sols = {r, p - r}
    pe = p
    for k in range(1, e):
        pe_next = pe * p
        nxt = set()
        for x in sols:
            # (x + b p^k)^2 ≡ a (mod p^{k+1})
            inv = pow(2 * x, -1, p)
            t = (a - x * x) // pe
            b = t * inv % p
            nxt.add((x + b * pe) % pe_next)
        sols = nxt
        pe = pe_next
    return sorted(sols)


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


def _crt_pair(r1, m1, r2, m2):
    g = gcd(m1, m2)
    if (r2 - r1) % g:
        return None
    # x = r1 + m1 t, t ≡ (r2-r1)/g * inv(m1/g) (mod m2/g)
    m1g, m2g = m1 // g, m2 // g
    t = (r2 - r1) // g * pow(m1g % m2g, -1, m2g) % m2g
    mod = m1 // g * m2
    return (r1 + m1 * t) % mod, mod


def sqrt_mod(a, n):
    """x^2 ≡ a (mod n) の 0 <= x < n の全解。n > 0。

    計算量: n の素因数分解 + 各 p^e で O(log^2 p + e)
    """
    if n <= 0:
        raise ValueError("n must be positive")
    a %= n
    if n == 1:
        return [0]
    fac = _factorize(n)
    parts = []
    for p, e in fac.items():
        pe = p ** e
        if p == 2:
            sols = _sqrt_mod_2e(a, e)
        else:
            sols = _hensel_odd(a, p, e)
        if not sols:
            return []
        parts.append((sols, pe))
    # 全組合せを CRT
    acc = [(0, 1)]
    for sols, pe in parts:
        nxt = []
        for r0, m0 in acc:
            for r1 in sols:
                cr = _crt_pair(r0, m0, r1, pe)
                if cr is not None:
                    nxt.append(cr)
        acc = nxt
        if not acc:
            return []
    return sorted({r for r, _ in acc})
