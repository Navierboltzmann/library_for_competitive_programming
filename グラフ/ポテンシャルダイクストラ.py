"""
ポテンシャル付きダイクストラ法 (負辺を含むグラフの単一始点最短路)

各辺 x→y (コスト c) について phi[y] - phi[x] <= c を満たすポテンシャル phi があれば、
c' = c - (phi[y] - phi[x]) >= 0 に変換してダイクストラできる。

使用例:
    e = [[(1, -2)], [(2, 3)], []]       # e[v] = [(行き先, コスト), ...]
    phi = [0, -2, 0]                    # phi[y] - phi[x] <= c を満たす (例えば -最長路 など)
    dijkstra_p(0, e, phi)               # [0, -2, 1]
"""
from heapq import heappop, heappush

INF = 1 << 60


def dijkstra_p(s, e, phi):
    """
    ポテンシャル phi を使って始点 s から各頂点への最短距離を返す O((N + M) logN)

    s: 始点
    e: 隣接リスト。e[v] = [(行き先, コスト), ...]  (コストは整数)
    phi: ポテンシャル。各辺 x→y (コスト c) について phi[y] - phi[x] <= c が成立しなければならない
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
        base = nw + phi[v]
        for iv, ic in e[v]:
            nc = base + ic - phi[iv]
            if nc < dist[iv]:
                dist[iv] = nc
                heappush(h, nc * N + iv)
    ps = phi[s]
    return [d + phi[i] - ps if d != INF else INF for i, d in enumerate(dist)]


if __name__ == "__main__":
    e = [[(1, -2)], [(2, 3)], []]
    print(dijkstra_p(0, e, [0, -2, 0]))  # [0, -2, 1]
