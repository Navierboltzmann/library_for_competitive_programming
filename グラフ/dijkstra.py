"""
ダイクストラ法 (単一始点最短路, 非負整数重み)

ヒープには (距離 * N + 頂点) の整数を入れて高速化している。

使用例:
    e = [[(1, 2), (2, 5)], [(2, 1)], []]   # e[v] = [(行き先, コスト), ...]
    dijkstra(0, e)                        # [0, 2, 3]
    到達不能な頂点の距離は INF。
"""
from heapq import heappop, heappush

INF = 1 << 60


def dijkstra(s, e):
    """
    始点 s から各頂点への最短距離を返す O((N + M) logN)

    s: 始点
    e: 隣接リスト。e[v] = [(行き先, コスト), ...]  (コストは非負整数)
    戻り値: 距離のリスト (到達不能は INF)
    """
    N = len(e)
    dist = [INF] * N
    dist[s] = 0
    h = [s]
    while h:
        nw, v = divmod(heappop(h), N)
        if dist[v] != nw:
            continue
        for iv, ic in e[v]:
            nc = ic + nw
            if nc < dist[iv]:
                dist[iv] = nc
                heappush(h, nc * N + iv)
    return dist


if __name__ == "__main__":
    e = [[(1, 2), (2, 5)], [(2, 1)], []]
    print(dijkstra(0, e))  # [0, 2, 3]
