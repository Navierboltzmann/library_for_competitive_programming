"""XOR 畳み込み (Walsh-Hadamard Transform)。

長さが 2 の冪である列 f, g に対し
h[k] = Σ_{i⊕j=k} f[i] g[j] を mod で計算する。
グローバル変数 `mod` を参照する。

使用例:
    mod = 998244353
    Xor_convolution([1, 2, 0, 0], [3, 4, 0, 0])
"""
mod = 998244353


def walsh_hadamard_transform(f):
    """f を Walsh-Hadamard 変換する(破壊的、in-place)。

    f: 長さ 2^k の列
    計算量: O(N log N)
    """
    n = len(f)
    h = 1
    while h < n:
        for i in range(0, n, h * 2):
            for j in range(i, i + h):
                x = f[j]
                y = f[j + h]
                f[j] = (x + y) % mod
                f[j + h] = (x - y) % mod
        h *= 2


def Xor_convolution(f, g):
    """XOR 畳み込みを計算し、結果を f に書き込んで返す。

    f, g: 同じ長さ 2^k の列(破壊的)
    計算量: O(N log N)
    """
    k = (len(f) - 1).bit_length()
    inv = pow(499122177, k, mod)
    walsh_hadamard_transform(f)
    walsh_hadamard_transform(g)
    n = 1 << k
    for i in range(n):
        f[i] = f[i] * g[i] % mod
    walsh_hadamard_transform(f)
    for i in range(n):
        f[i] = f[i] * inv % mod
    return f
