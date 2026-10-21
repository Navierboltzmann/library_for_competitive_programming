"""Stern-Brocot 木と Farey 数列。

正の既約分数を二分探索する。最良有理近似・隣接分数に使う。

使用例:
    sb_path(2, 5)                 # ['L', 'L', 'R']
    sb_ancestor(2, 5)             # (1, 3) の次が (2, 5) など
    farey_neighbors(2, 5, 5)      # (1, 3), (2, 5), (1, 2)
    best_rational_approx(355, 113, 10)  # 22/7 など
"""
from math import gcd


def sb_path(a, b):
    """1/1 から既約分数 a/b への L/R パス。

    L: 左の子 (p/q の左隣との中項)、R: 右の子。
    """
    g = gcd(a, b)
    a, b = a // g, b // g
    if a <= 0 or b <= 0:
        raise ValueError("require positive fraction")
    path = []
    la, lb, ra, rb = 0, 1, 1, 0
    while True:
        ma, mb = la + ra, lb + rb
        if (ma, mb) == (a, b):
            return path
        # a/b < ma/mb ?
        if a * mb < ma * b:
            path.append("L")
            ra, rb = ma, mb
        else:
            path.append("R")
            la, lb = ma, mb


def sb_from_path(path):
    """L/R パスが指す分数 (a, b)。"""
    la, lb, ra, rb = 0, 1, 1, 0
    for c in path:
        ma, mb = la + ra, lb + rb
        if c == "L":
            ra, rb = ma, mb
        else:
            la, lb = ma, mb
    return la + ra, lb + rb


def farey_neighbors(a, b, N):
    """分母 ≤ N の Farey 列で a/b の左隣・自身・右隣。

    a/b は既約で 0 <= a/b <= 1、b <= N を仮定。
    戻り値 ((pl, ql), (a, b), (pr, qr))
    """
    g = gcd(a, b)
    a, b = a // g, b // g
    aa, bb = a, b
    x0, y0, x1, y1 = 1, 0, 0, 1
    while bb:
        qot, aa, bb = aa // bb, bb, aa % bb
        x0, x1 = x1, x0 - qot * x1
        y0, y1 = y1, y0 - qot * y1
    # a x0 + b y0 = 1 → a q - b p = 1 は q=x0, p=-y0
    q, p = x0, -y0
    if q <= 0:
        t = (-q) // b + 1
        q += t * b
        p += t * a
    t = (N - q) // b
    q += t * b
    p += t * a
    left = (p, q)
    k = (N + q) // b
    right = (k * a - p, k * b - q)
    return left, (a, b), right


def mediate(a, b, c, d):
    """中項 (a+c)/(b+d)。"""
    return a + c, b + d


def best_rational_approx(p, q, max_den):
    """p/q を分母 ≤ max_den の既約分数で最良近似する。

    連分数の収束分数・中間分数から、|p/q - h/k| が最小のものを返す。
    """
    g = gcd(p, q)
    p, q = p // g, q // g
    if q <= max_den:
        return p, q
    h2, h1 = 0, 1
    k2, k1 = 1, 0
    a, b = p, q
    while b:
        ai = a // b
        a, b = b, a % b
        nxt_k = ai * k1 + k2
        if nxt_k > max_den:
            tmax = (max_den - k2) // k1 if k1 else 0
            cands = [(h1, k1)]
            if tmax > 0:
                cands.append((tmax * h1 + h2, tmax * k1 + k2))
            # |p/q - h/k| = |p k - q h| / (q k) を分母同士で比較
            best = cands[0]
            for cand in cands[1:]:
                h, k = cand
                hb, kb = best
                if abs(p * k - q * h) * kb < abs(p * kb - q * hb) * k:
                    best = cand
            return best
        h2, h1 = h1, ai * h1 + h2
        k2, k1 = k1, ai * k1 + k2
    return h1, k1
