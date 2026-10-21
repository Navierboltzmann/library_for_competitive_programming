"""Dinic 法による最大流。

使用例:
    dinic = Dinic(4)
    dinic.add_link(0, 1, 10)
    dinic.add_link(1, 3, 5)
    dinic.add_link(0, 2, 4)
    dinic.add_link(2, 3, 10)
    dinic.max_flow(0, 3)  # 9
"""
from collections import deque

INF = 10**18


class Dinic:
    """最大流 (Dinic)。"""

    def __init__(self, n):
        """n 頂点のフローネットワークを作る。

        n: 頂点数
        """
        self.n = n
        self.links = [[] for _ in range(n)]
        self.depth = None
        self.progress = None

    def add_link(self, _from, to, cap):
        """容量 cap の辺 _from → to を張る(逆辺も自動で張る)。

        計算量: O(1)
        """
        self.links[_from].append([cap, to, len(self.links[to])])
        self.links[to].append([0, _from, len(self.links[_from]) - 1])

    def bfs(self, s):
        """残留グラフ上で s からのレベルグラフを作る。"""
        n = self.n
        depth = [-1] * n
        depth[s] = 0
        q = deque([s])
        links = self.links
        while q:
            v = q.popleft()
            dv = depth[v]
            for cap, to, _rev in links[v]:
                if cap > 0 and depth[to] < 0:
                    depth[to] = dv + 1
                    q.append(to)
        self.depth = depth

    def dfs(self, v, t, flow):
        """レベルグラフ上で v → t に flow 以下を流す。

        再帰。深いグラフでは sys.setrecursionlimit か pypyjit を併用する。
        """
        if v == t:
            return flow
        links_v = self.links[v]
        depth = self.depth
        progress = self.progress
        for i in range(progress[v], len(links_v)):
            progress[v] = i
            link = links_v[i]
            cap, to, rev = link
            if cap == 0 or depth[v] >= depth[to]:
                continue
            d = self.dfs(to, t, flow if flow < cap else cap)
            if d == 0:
                continue
            link[0] -= d
            self.links[to][rev][0] += d
            return d
        return 0

    def max_flow(self, s, t):
        """s → t の最大流を返す。

        計算量: O(V^2 E)
        """
        flow = 0
        while True:
            self.bfs(s)
            if self.depth[t] < 0:
                return flow
            self.progress = [0] * self.n
            current_flow = self.dfs(s, t, INF)
            while current_flow > 0:
                flow += current_flow
                current_flow = self.dfs(s, t, INF)
