"""
ダイクストラ法 (経路復元つき, 非負整数重み)

使用例:
    graph = [[(1, 2), (2, 5)], [(2, 1)], []]   # graph[v] = [(行き先, コスト), ...]
    D = Dijkstra(graph, 0)
    D.min_cost(2)        # 3
    D.min_cost_path(2)   # [0, 1, 2]
    到達不能な頂点のコストは INF。
"""
from heapq import heappop, heappush

INF = 1 << 63


class Dijkstra:
    """
    始点 start からの最短距離と最短路木を構築時に計算する。O((N + M) logN)

    属性:
        cost[v] : start から v への最短距離 (到達不能は INF)
        prev[v] : 最短路における v の直前の頂点 (start・到達不能は -1)
    """

    def __init__(self, graph, start):
        """
        graph: 隣接リスト。graph[v] = [(行き先, コスト), ...]  (コストは非負整数)
        start: 始点
        """
        N = len(graph)
        cost = [INF] * N
        cost[start] = 0
        prev = [-1] * N
        pq = [start]
        while pq:
            w, v = divmod(heappop(pq), N)
            if w > cost[v]:
                continue
            for nv, nw in graph[v]:
                nc = w + nw
                if nc < cost[nv]:
                    cost[nv] = nc
                    prev[nv] = v
                    heappush(pq, nc * N + nv)
        self.cost = cost
        self.prev = prev

    def min_cost(self, goal: int):
        """start から goal への最短距離を返す O(1)"""
        return self.cost[goal]

    def min_cost_path(self, goal: int):
        """
        start から goal への最短路の頂点列を返す O(経路長)

        goal に到達不能な場合は [goal] を返す。
        """
        prev = self.prev
        nd = goal
        path = []
        while nd >= 0:
            path.append(nd)
            nd = prev[nd]
        return path[::-1]


if __name__ == "__main__":
    graph = [[(1, 2), (2, 5)], [(2, 1)], []]
    D = Dijkstra(graph, 0)
    print(D.min_cost(2), D.min_cost_path(2))  # 3 [0, 1, 2]
