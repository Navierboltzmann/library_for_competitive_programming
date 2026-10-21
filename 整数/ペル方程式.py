"""ペル方程式 x^2 - D y^2 = N。

連分数で x^2 - D y^2 = ±1 の基本解を求め、
一般の N は基本単数による同伴と、基本解の探索で解く。

使用例:
    pell(13)                 # (649, 180)
    negative_pell(13)        # (18, 5)
    pell(61)                 # (1766319049, 226153980)
    solve_pell(2, 7)         # [(3, 1), ...] の基本解
    next_pell(3, 1, 3, 2, 2) # (13, 9)   (3+√2)(3+2√2)
"""
from math import isqrt


def _is_square(n):
    """n が平方数なら ⌊√n⌋、違えば -1。"""
    if n < 0:
        return -1
    x = isqrt(n)
    return x if x * x == n else -1


def compose(x1, y1, x2, y2, D):
    """(x1 + y1 √D)(x2 + y2 √D) の成分。Brahmagupta 合成。"""
    return x1 * x2 + D * y1 * y2, x1 * y2 + y1 * x2


def pell_pow(x, y, D, k):
    """(x + y √D)^k の成分。k < 0 なら共役の正冪 (ノルムで割る)。

    ノルム N = x^2 - D y^2 として、k<0 のときは (x - y√D)^{-k} / N^{-k}。
    |N| = 1 なら整数のまま。
    計算量: O(log |k|)
    """
    if k == 0:
        return 1, 0
    if k < 0:
        nrm = x * x - D * y * y
        if nrm == 0:
            raise ValueError("zero norm")
        xp, yp = pell_pow(x, -y, D, -k)
        return xp // nrm, yp // nrm
    rx, ry = 1, 0
    while k:
        if k & 1:
            rx, ry = compose(rx, ry, x, y, D)
        x, y = compose(x, y, x, y, D)
        k >>= 1
    return rx, ry


def pell(D):
    """x^2 - D y^2 = 1 の最小正整数解 (x, y)。

    D が平方数・非正なら None。
    計算量: 連分数の周期に比例
    """
    if D <= 0:
        return None
    a0 = isqrt(D)
    if a0 * a0 == D:
        return None
    h2, h1 = 0, 1
    k2, k1 = 1, 0
    m, d, a = 0, 1, a0
    while True:
        h2, h1 = h1, a * h1 + h2
        k2, k1 = k1, a * k1 + k2
        if h1 * h1 - D * k1 * k1 == 1:
            return h1, k1
        m = d * a - m
        d = (D - m * m) // d
        a = (a0 + m) // d


def negative_pell(D):
    """x^2 - D y^2 = -1 の最小正整数解 (x, y)。無ければ None。

    連分数周期が奇数のときだけ解を持つ。
    """
    if D <= 0:
        return None
    a0 = isqrt(D)
    if a0 * a0 == D:
        return None
    h2, h1 = 0, 1
    k2, k1 = 1, 0
    m, d, a = 0, 1, a0
    while True:
        h2, h1 = h1, a * h1 + h2
        k2, k1 = k1, a * k1 + k2
        nrm = h1 * h1 - D * k1 * k1
        if nrm == -1:
            return h1, k1
        if nrm == 1:
            return None
        m = d * a - m
        d = (D - m * m) // d
        a = (a0 + m) // d


def next_pell(x, y, t, u, D):
    """解 (x, y) に基本単数 t + u √D を掛けた次の解。"""
    return compose(x, y, t, u, D)


def pell_nth(D, n):
    """x^2 - D y^2 = 1 の n 番目の正整数解 (n=1 が基本解)。"""
    fund = pell(D)
    if fund is None or n <= 0:
        return None
    return pell_pow(fund[0], fund[1], D, n)


def associated_solutions(x, y, D, max_count=None):
    """(x, y) を基本解とする x^2 - D y^2 = N の同伴解を生成する。

    基本単数 ε = t+u√D を右から掛け続ける。
    x, y は非負を仮定。符号は呼び出し側で付ける。
    """
    fund = pell(D)
    if fund is None:
        yield x, y
        return
    t, u = fund
    cnt = 0
    while True:
        yield x, y
        cnt += 1
        if max_count is not None and cnt >= max_count:
            return
        x, y = compose(x, y, t, u, D)
        if x < 0:
            x = -x
        if y < 0:
            y = -y


def _divisors_signed(n):
    """n の正負の約数。"""
    n = abs(n)
    small, large = [], []
    i = 1
    while i * i <= n:
        if n % i == 0:
            small.append(i)
            q = n // i
            if q != i:
                large.append(q)
        i += 1
    pos = small + large[::-1]
    return [-x for x in reversed(pos)] + pos


def _solve_square_D(s, N):
    """x^2 - s^2 y^2 = N、つまり (x-sy)(x+sy) = N。"""
    if s < 0:
        s = -s
    res = []
    seen = set()
    for f in _divisors_signed(N):
        if N % f:
            continue
        g = N // f
        if (g - f) % (2 * s) or (g + f) % 2:
            continue
        y = (g - f) // (2 * s)
        x = (g + f) // 2
        if y < 0:
            continue
        if x < 0:
            x = -x
        key = (x, y)
        if key not in seen and x * x - s * s * y * y == N:
            seen.add(key)
            res.append((x, y))
    res.sort()
    return res


def _cf_candidates(D, N):
    """√D の収束分数のうちノルムが N のもの。|N| < 2√D なら十分。"""
    a0 = isqrt(D)
    h2, h1 = 0, 1
    k2, k1 = 1, 0
    m, d, a = 0, 1, a0
    found = []
    first = True
    while True:
        h2, h1 = h1, a * h1 + h2
        k2, k1 = k1, a * k1 + k2
        if h1 * h1 - D * k1 * k1 == N:
            found.append((abs(h1), abs(k1)))
        if not first and a == 2 * a0:
            break
        first = False
        m = d * a - m
        d = (D - m * m) // d
        a = (a0 + m) // d
    return found


def _search_fundamentals(D, N, t, u, max_y=None):
    """Nagell の限界まで y を走査して各類の基本解を探す。

    0 < x <= sqrt(|N|(t+1)/2) を満たす解が各類に 1 つある。
    max_y を指定するとその上限で打ち切る(不完全になり得る)。
    """
    absN = abs(N)
    # x <= sqrt(|N|(t+1)/2)
    # y^2 = (x^2 - N)/D <= (|N|(t+1)/2 - N)/D
    bound_x = isqrt(absN * (t + 1) // 2) + 2
    if bound_x * bound_x < D:
        y_lim = 1
    else:
        num = bound_x * bound_x - N
        y_lim = isqrt(num // D) + 2 if num > 0 else 1
    if N < 0:
        y_min = isqrt((-N + D - 1) // D)
    else:
        y_min = 0
    if max_y is not None:
        y_lim = min(y_lim, max_y)
    res = []
    for y in range(y_min, y_lim + 1):
        x = _is_square(N + D * y * y)
        if x >= 0:
            res.append((x, y))
    return res


def _reduce_fundamental(x, y, D, t, u, N):
    """同伴類の中で x + y√D が最小の非負解へ落とす。"""
    x, y = abs(x), abs(y)
    while True:
        nx = abs(t * x - D * u * y)
        ny = abs(-u * x + t * y)
        if nx * nx - D * ny * ny != N:
            break
        if nx * nx + D * ny * ny >= x * x + D * y * y:
            break
        x, y = nx, ny
    return x, y


def _normalize_class(sols, D, t, u, N):
    """同伴を除き、各類の最小非負解だけ残す。"""
    seen = set()
    res = []
    for x, y in sols:
        x, y = _reduce_fundamental(x, y, D, t, u, N)
        if (x, y) not in seen and x * x - D * y * y == N:
            seen.add((x, y))
            res.append((x, y))
    res.sort()
    return res


def solve_pell(D, N=1, max_y=2_000_000):
    """x^2 - D y^2 = N の、各同伴類の非負基本解を返す。

    D > 0, N ≠ 0。
    D が平方なら因数分解で完全に解く。
    そうでなければ連分数 + Nagell 限界の走査。
    限界が max_y を超える類は取りこぼす可能性がある。
    計算量: =±1 は周期長、一般は O(min(限界, max_y))
    """
    if D <= 0 or N == 0:
        raise ValueError("require D > 0 and N != 0")
    s = isqrt(D)
    if s * s == D:
        return _solve_square_D(s, N)

    if N == 1:
        fund = pell(D)
        return [fund] if fund else []
    if N == -1:
        fund = negative_pell(D)
        return [fund] if fund else []

    t, u = pell(D)
    cands = _cf_candidates(D, N)
    # Nagell 限界が実用範囲なら完全探索
    bound_x = isqrt(abs(N) * (t + 1) // 2) + 2
    y_lim_est = isqrt(max(0, bound_x * bound_x - N) // D) + 2
    if y_lim_est <= max_y:
        cands.extend(_search_fundamentals(D, N, t, u, max_y))
    else:
        cands.extend(_search_fundamentals(D, N, t, u, max_y))

    seen = set()
    uniq = []
    for x, y in cands:
        if x < 0:
            x = -x
        if y < 0:
            y = -y
        if (x, y) in seen:
            continue
        if x * x - D * y * y != N:
            continue
        seen.add((x, y))
        uniq.append((x, y))
    return _normalize_class(uniq, D, t, u, N)


def all_pell_solutions(D, N=1, max_count=10, max_y=2_000_000):
    """各類について最大 max_count 個までの非負解を小さい順に近い形で返す。"""
    res = []
    for x, y in solve_pell(D, N, max_y=max_y):
        for i, sol in enumerate(associated_solutions(x, y, D, max_count)):
            res.append(sol)
            if i + 1 >= max_count:
                break
    res = sorted(set(res))
    return res
