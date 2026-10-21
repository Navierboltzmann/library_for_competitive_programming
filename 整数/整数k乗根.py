"""整数の k 乗根 ⌊n^{1/k}⌋ と、それが完全 k 乗か。

使用例:
    iroot(10, 2)          # 3
    iroot(27, 3)          # 3
    is_perfect_power(16)  # (2, 4)   2^4
"""


def iroot(n, k):
    """⌊n^{1/k}⌋ を返す。n >= 0, k >= 1。

    ニュートン法。Python の任意精度整数で正確。
    計算量: O(log n · M(bit(n))) 程度
    """
    if k <= 0:
        raise ValueError("k must be positive")
    if n < 0:
        if k % 2 == 0:
            raise ValueError("even root of negative")
        return -iroot(-n, k)
    if n <= 1 or k == 1:
        return n
    # 初期値: 2^{ceil(bitlen / k)}
    x = 1 << ((n.bit_length() + k - 1) // k)
    while True:
        y = ((k - 1) * x + n // x ** (k - 1)) // k
        if y >= x:
            return x
        x = y


def is_kth_power(n, k):
    """n が完全 k 乗か。そうなら根を返す。"""
    if n < 0:
        if k % 2 == 0:
            return None
        r = iroot(-n, k)
        return -r if (-r) ** k == n else None
    r = iroot(n, k)
    return r if r ** k == n else None


def is_perfect_power(n):
    """n が完全冪 a^b (b>=2) なら (a, b) を返す。b は最大。

    1 は (1, 2) とする。n < 0 は奇数指数のみ。
    """
    if n == 0:
        return 0, 2
    if n == 1:
        return 1, 2
    sign = 1
    if n < 0:
        sign = -1
        n = -n
    max_k = n.bit_length()
    for k in range(max_k, 1, -1):
        if sign < 0 and k % 2 == 0:
            continue
        r = is_kth_power(n, k)
        if r is not None:
            return sign * r, k
    return None
