"""二部グラフの最大マッチング (Dinic)。

左側 nL 頂点、右側 nR 頂点。辺は左→右。

使用例:
    bm = BipartiteMatching(2, 2)
    bm.add_edge(0, 0)
    bm.add_edge(0, 1)
    bm.add_edge(1, 1)
    bm.max_matching()        # 2
    bm.matching_edges()      # [(0, 0), (1, 1)] など
"""
from collections import deque


class BipartiteMatching:
    """二部最大マッチング。"""

    def __init__(self, nL, nR):
        """左側 nL、右側 nR の空グラフ。"""
        self.nL = nL
        self.nR = nR
        self.n = nL + nR + 2
        self.s = nL + nR
        self.t = self.s + 1
        self.graph = [[] for _ in range(self.n)]

    def add_edge(self, u, v):
        """左 u から右 v へ辺を張る。"""
        self._add(u, self.nL + v, 1)

    def _add(self, fr, to, cap):
        self.graph[fr].append([to, cap, len(self.graph[to])])
        self.graph[to].append([fr, 0, len(self.graph[fr]) - 1])

    def max_matching(self):
        """最大マッチング数。計算量 O(E √V)。"""
        for u in range(self.nL):
            self._add(self.s, u, 1)
        for v in range(self.nR):
            self._add(self.nL + v, self.t, 1)
        flow = 0
        INF = 10 ** 9
        while True:
            level = self._bfs()
            if level[self.t] < 0:
                break
            it = [0] * self.n
            while True:
                f = self._dfs(self.s, INF, level, it)
                if f == 0:
                    break
                flow += f
        self._flow = flow
        return flow

    def _bfs(self):
        level = [-1] * self.n
        q = deque([self.s])
        level[self.s] = 0
        while q:
            v = q.popleft()
            for to, cap, _ in self.graph[v]:
                if cap > 0 and level[to] < 0:
                    level[to] = level[v] + 1
                    q.append(to)
        return level

    def _dfs(self, v, f, level, it):
        if v == self.t:
            return f
        g = self.graph[v]
        while it[v] < len(g):
            to, cap, rev = g[it[v]]
            if cap > 0 and level[v] < level[to]:
                d = self._dfs(to, min(f, cap), level, it)
                if d:
                    g[it[v]][1] -= d
                    self.graph[to][rev][1] += d
                    return d
            it[v] += 1
        return 0

    def matching_edges(self):
        """マッチした (左, 右) のリスト。max_matching の後に呼ぶ。"""
        res = []
        for u in range(self.nL):
            for to, cap, _ in self.graph[u]:
                if self.nL <= to < self.nL + self.nR and cap == 0:
                    res.append((u, to - self.nL))
        return res
