"""転倒数 (inversion) をマージソートで数える。

A を破壊的にソートしつつ、i < j かつ A[i] > A[j] となる組の数を返す。

使用例:
    mergecount([3, 1, 2])  # 2
"""


def mergecount(A):
    """列 A の転倒数を返し、A を昇順にソートする。

    A: 比較可能な列(破壊的)
    計算量: O(N log N)
    """
    n = len(A)
    if n <= 1:
        return 0
    mid = n >> 1
    A1 = A[:mid]
    A2 = A[mid:]
    cnt = mergecount(A1) + mergecount(A2)
    i1 = i2 = 0
    n1 = mid
    n2 = n - mid
    for i in range(n):
        if i2 == n2:
            A[i] = A1[i1]
            i1 += 1
        elif i1 == n1:
            A[i] = A2[i2]
            i2 += 1
        elif A1[i1] <= A2[i2]:
            A[i] = A1[i1]
            i1 += 1
        else:
            A[i] = A2[i2]
            i2 += 1
            cnt += n1 - i1
    return cnt
