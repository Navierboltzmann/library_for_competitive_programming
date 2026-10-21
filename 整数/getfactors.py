"""素因数分解(試し割り)。

使用例:
    get_factors(12)  # {2: 2, 3: 1}
"""


def get_factors(x):
    """x を試し割りで素因数分解し {素数: 指数} を返す。

    x: 1 以上の整数(1 なら空 dict)
    計算量: O(√x)
    """
    res = {}
    if x == 1:
        return res
    now = x
    i = 2
    while i * i <= now:
        if now % i == 0:
            cnt = 0
            while now % i == 0:
                now //= i
                cnt += 1
            res[i] = cnt
        i += 1
    if now != 1:
        res[now] = 1
    return res
