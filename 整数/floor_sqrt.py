"""整数平方根 ⌊√n⌋。

使用例:
    floor_sqrt(10)  # 3
"""


def floor_sqrt(n):
    """n の整数平方根 ⌊√n⌋ を返す。

    n: 非負整数
    計算量: 組み込み int.bit_length を使ったニュートン法相当
    """
    if n <= 1:
        return n
    x = 1 << ((n.bit_length() + 1) >> 1)
    while True:
        y = (x + n // x) >> 1
        if y >= x:
            return x
        x = y
