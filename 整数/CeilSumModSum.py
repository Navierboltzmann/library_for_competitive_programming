"""
商の和・余りの和を O(√N) で求めるライブラリ

- CeilSum(N)  : sum_{i=1}^{N} floor(N / i)
- iCeilSum(N) : sum_{i=1}^{N} floor(N / i) * i
- ModSum(N)   : sum_{i=1}^{N} (N mod i)

(名前は CeilSum だが中身は floor の和)

使用例:
    CeilSum(10)   # -> 27
    iCeilSum(10)  # -> 73
    ModSum(10)    # -> 27
"""


def CeilSum(N):
    """
    sum_{i=1}^{N} N//i を求める

    N: 0以上の整数
    return: 和
    計算量: O(√N)
    """
    ans = 0
    i = 1
    while i <= N:
        q = N // i
        j = N // q
        ans += q * (j - i + 1)
        i = j + 1
    return ans


def iCeilSum(N):
    """
    sum_{i=1}^{N} (N//i)*i を求める

    N: 0以上の整数
    return: 和
    計算量: O(√N)
    """
    ans = 0
    i = 1
    while i <= N:
        q = N // i
        j = N // q
        ans += q * (i + j) * (j - i + 1) // 2
        i = j + 1
    return ans


def ModSum(N):
    """
    sum_{i=1}^{N} N%i を求める

    N: 0以上の整数
    return: 和
    計算量: O(√N)
    """
    return N * N - iCeilSum(N)


if __name__ == "__main__":
    print(CeilSum(10), iCeilSum(10), ModSum(10))
