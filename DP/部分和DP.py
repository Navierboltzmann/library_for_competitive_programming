"""部分和の通り数 DP。

A の部分集合で和がちょうど W になるものの個数を mod で数える。
グローバル変数 `mod` を参照する。

使用例:
    mod = 998244353
    subset_sum([1, 2, 3], 3)  # 2  ({3} と {1,2})
"""


def subset_sum(A, W):
    """A の部分集合で和が W になる通り数を mod で返す。

    A: 正整数の列
    W: 目標和
    計算量: O(N W)
    """
    dp = [0] * (W + 1)
    dp[0] = 1
    for a in A:
        for j in range(W, a - 1, -1):
            dp[j] = (dp[j] + dp[j - a]) % mod
    return dp[W]
