"""約数列挙。

使用例:
    md(12)  # [1, 2, 3, 4, 6, 12]
"""


def md(n):
    """n の約数を昇順で返す。

    n: 正整数
    計算量: O(√n)
    """
    lower, upper = [], []
    i = 1
    while i * i <= n:
        if n % i == 0:
            lower.append(i)
            q = n // i
            if i != q:
                upper.append(q)
        i += 1
    lower.extend(reversed(upper))
    return lower
