"""Lucas の定理による二項係数の剰余。

素数 p に対し C(m, n) mod p を、m, n を p 進展開して計算する。
p が 2 未満のときは通常の二項係数にフォールバックする。

使用例:
    lucas(10, 3, 7)  # C(10, 3) % 7
"""


def combi(n, k):
    """二項係数 C(n, k) を整数で返す。

    n, k: 非負整数
    計算量: O(min(k, n-k))
    """
    if n < k or k < 0:
        return 0
    k = min(k, n - k)
    if k == 0:
        return 1
    val = 1
    for i in range(k):
        val = val * (n - i) // (i + 1)
    return val


def lucas(m, n, p):
    """Lucas の定理で C(m, n) mod p を返す。

    m, n: 非負整数
    p: 素数(2 未満なら combi にフォールバック)
    計算量: O(log_p max(m, n) * p)
    """
    if p < 2:
        return combi(m, n)
    ans = 1
    while m or n:
        m, dm = divmod(m, p)
        n, dn = divmod(n, p)
        if dm < dn:
            return 0
        ans = ans * (combi(dm, dn) % p) % p
    return ans
