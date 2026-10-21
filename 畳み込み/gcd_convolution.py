"""GCD 畳み込みと乗法的ゼータ / メビウス。

配列は 1-index で、A[0] はダミー。

使用例:
    A = [0, 1, 2, 3]
    B = [0, 4, 5, 6]
    gcd_conv(A, B)  # C[g] = Σ_{gcd(i,j)=g} A[i] B[j]
"""


def mul_zeta(A):
    """倍数ゼータ: B[i] = Σ_{i|j} A[j]。

    計算量: O(N log log N)
    """
    X = A[:]
    N = len(A) - 1
    sieve = [False] * (N + 1)
    for i in range(2, N + 1):
        if sieve[i]:
            continue
        for j in range(N // i, 0, -1):
            sieve[i * j] = True
            X[j] += X[i * j]
    return X


def mul_mobius(A):
    """倍数ゼータの逆変換。

    計算量: O(N log log N)
    """
    X = A[:]
    N = len(A) - 1
    sieve = [False] * (N + 1)
    for i in range(2, N + 1):
        if sieve[i]:
            continue
        for j in range(1, N // i + 1):
            sieve[i * j] = True
            X[j] -= X[i * j]
    return X


def div_zeta(A):
    """約数ゼータ: B[i] = Σ_{j|i} A[j]。

    計算量: O(N log log N)
    """
    X = A[:]
    N = len(A) - 1
    sieve = [False] * (N + 1)
    for i in range(2, N + 1):
        if sieve[i]:
            continue
        for j in range(1, N // i + 1):
            sieve[i * j] = True
            X[i * j] += X[j]
    return X


def gcd_conv(A, B):
    """GCD 畳み込み C[g] = Σ_{gcd(i,j)=g} A[i] B[j]。

    計算量: O(N log log N)
    """
    X = mul_zeta(A)
    Y = mul_zeta(B)
    return mul_mobius([a * b for a, b in zip(X, Y)])
