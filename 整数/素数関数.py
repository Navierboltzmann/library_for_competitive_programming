"""素数計数 π(n) (Lucy-Hedgehog / min_25 の簡易版)。

使用例:
    count_primes(10)       # 4
    count_primes(10**10)   # 455052511
"""


def count_primes(n):
    """n 以下の素数の個数を返す。

    n: 非負整数
    計算量: O(n^{3/4} / log n) 程度
    """
    if n < 2:
        return 0
    v = int(n ** 0.5) + 1
    smalls = [i // 2 for i in range(1, v + 1)]
    smalls[1] = 0
    s = v // 2
    roughs = [2 * i + 1 for i in range(s)]
    larges = [(n // (2 * i + 1) + 1) // 2 for i in range(s)]
    skip = [False] * v
    pc = 0
    for p in range(3, v):
        if smalls[p] <= smalls[p - 1]:
            continue
        q = p * p
        pc += 1
        if q * q > n:
            break
        skip[p] = True
        for i in range(q, v, 2 * p):
            skip[i] = True
        ns = 0
        for k in range(s):
            i = roughs[k]
            if skip[i]:
                continue
            d = i * p
            larges[ns] = larges[k] - (larges[smalls[d] - pc] if d < v else smalls[n // d]) + pc
            roughs[ns] = i
            ns += 1
        s = ns
        for j in range((v - 1) // p, p - 1, -1):
            c = smalls[j] - pc
            e = min((j + 1) * p, v)
            for i in range(j * p, e):
                smalls[i] -= c
    for k in range(1, s):
        m = n // roughs[k]
        s2 = larges[k] - (pc + k - 1)
        for l in range(1, k):
            p = roughs[l]
            if p * p > m:
                break
            s2 -= smalls[m // p] - (pc + l - 1)
        larges[0] -= s2
    return larges[0]
