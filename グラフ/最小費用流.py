"""下限付き最小費用フロー (Cost Scaling)。

参考: https://judge.yosupo.jp/submission/26494

使用例:
    flow = MinCostFlow(3)
    flow.add_supply(0, 2)
    flow.add_edge(0, 1, 0, 2, 5)
    flow.add_edge(1, 2, 0, 2, 3)
    ok = flow.solve()
"""
from heapq import heappop, heappush

INF = 1 << 60


class MinCostFlow:
    """下限付き最小費用流。"""

    def __init__(self, n):
        """n 頂点のネットワークを作る。

        n: 頂点数
        """
        self.n = n
        self.edges = []
        self._E = [[] for _ in range(n + 1)]
        self._b = [0] * (n + 1)
        self._u = 1

    def add_supply(self, v, amount):
        """頂点 v に供給量 amount を与える。"""
        self.add_edge(self.n, v, amount, amount, 0, False)

    def add_edge(self, u, v, lower, upper, cost, append=True):
        """u → v に流量 [lower, upper]、コスト cost の辺を張る。

        append: True なら self.edges に記録する
        """
        e = [v, upper, cost, 0]
        rev = [u, -lower, -cost, e]
        e[-1] = rev
        self._E[u].append(e)
        self._E[v].append(rev)
        self._u = max(self._u, upper - lower)
        if not lower <= 0 <= upper:
            self._push(e, lower)
        if append:
            self.edges.append((e, upper))

    def _push(self, e, f):
        """辺 e に流量 f を流す。"""
        rev = e[-1]
        self._b[rev[0]] -= f
        self._b[e[0]] += f
        e[1] -= f
        rev[1] += f

    def _saturate(self, delta):
        """負辺を delta 以上の容量で飽和させる。"""
        p = self._p
        for v in range(self.n + 1):
            pv = p[v]
            for e in self._E[v]:
                nv, cap, cost, _ = e
                if cap >= delta and cost + pv - p[nv] < 0:
                    self._push(e, cap)

    def _shortest_path(self, delta):
        """ポテンシャル付きダイクストラで増加路を探す。"""
        n = self.n + 1
        dist = [INF] * n
        heap = []
        for v, b in enumerate(self._b):
            if b >= delta:
                dist[v] = 0
                heap.append(v)
        prev = [None] * n
        ward = n.bit_length()
        mask = (1 << ward) - 1
        p = self._p
        while heap:
            packed = heappop(heap)
            d, v = packed >> ward, packed & mask
            if dist[v] != d:
                continue
            if self._b[v] <= -delta:
                t = v
                f = -self._b[v]
                while prev[v]:
                    cap = prev[v][1]
                    if cap < f:
                        f = cap
                    v = prev[v][-1][0]
                if self._b[v] < f:
                    f = self._b[v]
                while prev[t]:
                    self._push(prev[t], f)
                    t = prev[t][-1][0]
                for i in range(n):
                    di = dist[i]
                    p[i] += d if d < di else di
                return True
            pv = p[v]
            for e in self._E[v]:
                nv, cap, cost, _ = e
                if cap < delta:
                    continue
                nd = d + cost + pv - p[nv]
                if dist[nv] > nd:
                    prev[nv] = e
                    dist[nv] = nd
                    heappush(heap, nd << ward | nv)
        return False

    def solve(self):
        """フローを流し、実行可能なら True を返す。

        辺の実流量は (upper - e[1]) で取れる。
        """
        self._p = [0] * (self.n + 1)
        delta = 1 << (self._u - 1).bit_length()
        while delta:
            self._saturate(delta)
            while self._shortest_path(delta):
                pass
            delta >>= 1
        return not any(self._b)
