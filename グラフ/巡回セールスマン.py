"""巡回セールスマン問題 (bit DP)。

cost[v][u] を与え、頂点 0 から全頂点を訪れて 0 に戻る最短路を求める。

使用例:
    cost = [[0, 1, 5], [1, 0, 2], [5, 2, 0]]
    tsp(cost)  # 8
"""

INF = 10**10


def tsp(cost, start=0):
    """距離行列 cost の TSP 最短路長を返す。到達不能なら INF。

    cost: V×V の隣接行列
    start: 始点(終点でもある)
    計算量: O(V^2 2^V)
    """
    V = len(cost)
    full = (1 << V) - 1
    dp = [[INF] * V for _ in range(1 << V)]
    dp[1 << start][start] = 0
    for S in range(1 << V):
        dps = dp[S]
        for v in range(V):
            if not (S >> v) & 1:
                continue
            cur = dps[v]
            if cur >= INF:
                continue
            cv = cost[v]
            for u in range(V):
                if (S >> u) & 1:
                    continue
                nxt = cur + cv[u]
                nS = S | (1 << u)
                if nxt < dp[nS][u]:
                    dp[nS][u] = nxt
    ans = INF
    cf = cost
    for v in range(V):
        t = dp[full][v] + cf[v][start]
        if t < ans:
            ans = t
    return ans
