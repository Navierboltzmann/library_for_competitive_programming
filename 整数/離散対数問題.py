"""離散対数: x^k ≡ y (mod m) となる最小の k ≥ 0。

Baby-step giant-step。解なしは -1。

使用例:
    solveDiscreteLogarithm(3, 13, 17)  # 4  (3^4=81≡13)
"""
from math import gcd, isqrt


def solveDiscreteLogarithm(x, y, m):
    """x^k ≡ y (mod m) の最小非負 k を返す。無ければ -1。

    x, y, m: 整数 (m > 0)
    計算量: O(√φ(m') + 素因数分解)
    """
    if y >= m or y < 0:
        return -1
    if x == 0:
        if m == 1:
            return 0
        if y == 1:
            return 0
        if y == 0:
            return 1
        return -1
    p = 3
    tmp = x - 1
    cnt = 0
    primes = []
    counts = []
    ps = 0
    while tmp & 1:
        tmp >>= 1
        cnt += 1
    if cnt:
        primes.append(2)
        counts.append(cnt)
        ps += 1
    tmp += 1
    while tmp != 1:
        cnt = 0
        while tmp % p == 0:
            tmp //= p
            cnt += 1
        if cnt:
            primes.append(p)
            counts.append(cnt)
            ps += 1
        p += 2
        if tmp != 1 and p * p > x:
            primes.append(tmp)
            counts.append(1)
            ps += 1
            break
    tail = 0
    mp = m
    for i in range(ps):
        f = 0
        while mp % primes[i] == 0:
            mp //= primes[i]
            f += 1
        t = (f + counts[i] - 1) // counts[i]
        if tail < t:
            tail = t
    z = 1
    for i in range(tail):
        if z == y:
            return i
        z = z * x % m
    if y % gcd(z, m):
        return -1
    p = 3
    u = mp
    tmp = mp - 1
    if tmp & 1:
        u >>= 1
        while tmp & 1:
            tmp >>= 1
    tmp += 1
    while tmp != 1:
        if tmp % p == 0:
            u //= p
            u *= p - 1
            while tmp % p == 0:
                tmp //= p
        p += 2
        if tmp != 1 and p * p > mp:
            u //= tmp
            u *= tmp - 1
            break
    p = 1
    loop = u
    while p * p <= u:
        if u % p == 0:
            if z * pow(x, p, m) % m == z:
                loop = p
                break
            ip = u // p
            if z * pow(x, ip, m) % m == z:
                loop = ip
        p += 1
    sq = isqrt(loop)
    if sq * sq < loop:
        sq += 1
    b = pow(pow(x, loop - 1, m), sq, m)
    d = {}
    f = z
    for i in range(sq):
        d[f] = i
        f = f * x % m
    g = y
    for i in range(sq):
        if g in d:
            return i * sq + d[g] + tail
        g = g * b % m
    return -1
