"""AND / OR / XOR 畳み込み (高速ゼータ・メビウス)。

長さ 2^n の列 a, b を破壊的に畳み込む。結果は a に書き込まれる。

使用例:
    a = [1, 2, 3, 4]
    b = [5, 6, 7, 8]
    xor_convolution(a, b[:], 2)
"""
MOD = 998244353


def fztAND(vec, n):
    """AND ゼータ変換(上位集合和)。破壊的。"""
    N = 1 << n
    for k in range(n):
        i = 1 << k
        for j in range(N):
            if not i & j:
                vec[j] = (vec[j] + vec[i | j]) % MOD


def fztOR(vec, n):
    """OR ゼータ変換(下位集合和)。破壊的。"""
    N = 1 << n
    for k in range(n):
        i = 1 << k
        for j in range(N):
            if not i & j:
                vec[i | j] = (vec[i | j] + vec[j]) % MOD


def fztXOR(vec, n):
    """XOR ゼータ (Walsh-Hadamard)。破壊的。"""
    N = 1 << n
    for k in range(n):
        i = 1 << k
        for j in range(N):
            if not i & j:
                x = vec[j]
                y = vec[i | j]
                vec[j] = (x + y) % MOD
                vec[i | j] = (x - y) % MOD


def fmtAND(vec, n):
    """AND メビウス変換。破壊的。"""
    N = 1 << n
    for k in range(n):
        i = 1 << k
        for j in range(N):
            if not i & j:
                vec[j] = (vec[j] - vec[i | j]) % MOD


def fmtOR(vec, n):
    """OR メビウス変換。破壊的。"""
    N = 1 << n
    for k in range(n):
        i = 1 << k
        for j in range(N):
            if not i & j:
                vec[i | j] = (vec[i | j] - vec[j]) % MOD


def and_convolution(a, b, n):
    """AND 畳み込み。結果は a に書き込む。

    計算量: O(n 2^n)
    """
    fztAND(a, n)
    fztAND(b, n)
    N = 1 << n
    for i in range(N):
        a[i] = a[i] * b[i] % MOD
    fmtAND(a, n)


def or_convolution(a, b, n):
    """OR 畳み込み。結果は a に書き込む。

    計算量: O(n 2^n)
    """
    fztOR(a, n)
    fztOR(b, n)
    N = 1 << n
    for i in range(N):
        a[i] = a[i] * b[i] % MOD
    fmtOR(a, n)


def xor_convolution(a, b, n):
    """XOR 畳み込み。結果は a に書き込む。

    計算量: O(n 2^n)
    """
    fztXOR(a, n)
    fztXOR(b, n)
    N = 1 << n
    for i in range(N):
        a[i] = a[i] * b[i] % MOD
    fztXOR(a, n)
    invM = pow(N, MOD - 2, MOD)
    for i in range(N):
        a[i] = a[i] * invM % MOD
