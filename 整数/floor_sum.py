"""商の和 floor_sum。

sum_{i=0}^{n-1} floor((a i + b) / m) を O(log(n+m+|a|+|b|)) で求める。
ACL と同じ仕様。a, b は負でもよい。

使用例:
    floor_sum(4, 10, 6, 3)  # 0+0+1+2 = 3
"""


def floor_sum(n, m, a, b):
    """Σ_{i=0}^{n-1} ⌊(a i + b) / m⌋ を返す。

    n >= 0, m >= 1。a, b は任意の整数。
    計算量: O(log (n + m + |a| + |b|))
    """
    if n < 0 or m <= 0:
        raise ValueError("require n >= 0 and m >= 1")
    ans = 0
    if a < 0:
        a2 = a % m
        ans -= n * (n - 1) // 2 * ((a2 - a) // m)
        a = a2
    if b < 0:
        b2 = b % m
        ans -= n * ((b2 - b) // m)
        b = b2
    return ans + _floor_sum_unsigned(n, m, a, b)


def _floor_sum_unsigned(n, m, a, b):
    """0 <= a, b の内部実装。"""
    ans = 0
    while True:
        if a >= m:
            ans += n * (n - 1) // 2 * (a // m)
            a %= m
        if b >= m:
            ans += n * (b // m)
            b %= m
        y_max = a * n + b
        if y_max < m:
            break
        n = y_max // m
        b = y_max % m
        m, a = a, m
    return ans


def floor_sum_range(L, R, m, a, b):
    """Σ_{i=L}^{R-1} ⌊(a i + b) / m⌋。"""
    if L >= R:
        return 0
    return floor_sum(R, m, a, b) - floor_sum(L, m, a, b)
