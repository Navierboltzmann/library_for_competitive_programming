"""有理数・平方根の連分数と収束分数。

使用例:
    continued_fraction_rational(22, 7)   # [3, 7]
    a0, period = continued_fraction_sqrt(13)  # 3, [1, 1, 1, 1, 6]
    list(convergents([1, 2, 2, 2]))      # [(1, 1), (3, 2), (7, 5), ...]
"""
from math import isqrt


def continued_fraction_rational(a, b):
    """有理数 a/b (b>0) の連分数係数 [a0; a1, ...] を返す。

    計算量: O(log max(|a|, b))
    """
    if b <= 0:
        raise ValueError("b must be positive")
    terms = []
    while b:
        q, r = divmod(a, b)
        terms.append(q)
        a, b = b, r
    return terms


def continued_fraction_sqrt(D):
    """√D の連分数 [a0; 周期] を (a0, period) で返す。

    D が平方数なら period は空。
    周期は a = 2 a0 で終わる。
    計算量: 周期長に比例。周期は O(√D) 以下
    """
    if D < 0:
        raise ValueError("D must be non-negative")
    a0 = isqrt(D)
    if a0 * a0 == D:
        return a0, []
    period = []
    m, d, a = 0, 1, a0
    while True:
        m = d * a - m
        d = (D - m * m) // d
        a = (a0 + m) // d
        period.append(a)
        if a == 2 * a0:
            return a0, period


def convergents(terms):
    """連分数係数から収束分数 (h, k) = h/k を順に生成する。"""
    h_prev2, h_prev1 = 0, 1
    k_prev2, k_prev1 = 1, 0
    for a in terms:
        h = a * h_prev1 + h_prev2
        k = a * k_prev1 + k_prev2
        yield h, k
        h_prev2, h_prev1 = h_prev1, h
        k_prev2, k_prev1 = k_prev1, k


def evaluate_continued_fraction(terms):
    """連分数係数を既約分数 (h, k) に戻す。"""
    last = None
    for last in convergents(terms):
        pass
    return last if last is not None else (0, 1)


def best_approximations(terms):
    """収束分数と中間分数(半収束分数)を順に生成する。

    |α - h/k| が分母 k 以下で最良のものすべて。
    """
    terms = list(terms)
    if not terms:
        return
    h_prev2, h_prev1 = 0, 1
    k_prev2, k_prev1 = 1, 0
    for i, a in enumerate(terms):
        start = 1 if i else a
        for t in range(start, a + 1):
            h = t * h_prev1 + h_prev2
            k = t * k_prev1 + k_prev2
            yield h, k
        h_prev2, h_prev1 = h_prev1, a * h_prev1 + h_prev2
        k_prev2, k_prev1 = k_prev1, a * k_prev1 + k_prev2


def sqrt_convergents(D, periods=2):
    """√D の収束分数を periods 周期分生成する。平方数なら (a0, 1) のみ。"""
    a0, period = continued_fraction_sqrt(D)
    yield from convergents([a0] + period * periods)
