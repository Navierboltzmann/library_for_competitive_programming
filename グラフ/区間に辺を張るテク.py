"""区間から区間へ辺を張るテクニック。

セグ木 2 本(入次数用・出次数用)を使い、
[l1, r1) の全頂点から [l2, r2) の全頂点へコスト c の辺を
O(log N) 頂点を経由して張る。
頂点 i の本体は n+i (n は N 以上の 2 冪)。

使用例:
    g = IntervalEdges(N=5)
    g.add(0, 3, 1, 4, 10)
    # g.e が隣接リスト
"""


class IntervalEdges:
    """区間→区間の辺を対数個の辺で表現する。"""

    def __init__(self, N, Q=0):
        """頂点数 N、追加エッジ最大 Q 本でグラフを構築する。

        N: 元の頂点数
        Q: add の最大呼び出し回数
        """
        n = 1
        while n < N:
            n *= 2
        self.n = n
        self.N3 = 3 * n + 1
        self.q = 0
        self.e = [[] for _ in range(self.N3 + 2 * Q)]
        e = self.e
        for i in range(n - 1, 0, -1):
            a = 2 * i
            b = 2 * i + 1
            e[i].append((a, 0))
            e[i].append((b, 0))
            if n <= a < 2 * n:
                e[a].append((i + 2 * n, 0))
            else:
                e[a + 2 * n].append((i + 2 * n, 0))
            if n <= b < 2 * n:
                e[b].append((i + 2 * n, 0))
            else:
                e[b + 2 * n].append((i + 2 * n, 0))

    def add(self, l1, r1, l2, r2, c):
        """[l1, r1) から [l2, r2) へコスト c の辺を張る。

        計算量: O(log N)
        """
        n = self.n
        e = self.e
        N3 = self.N3
        q = self.q
        source = []
        l1 += n
        r1 += n
        while l1 < r1:
            if l1 & 1:
                source.append(l1 if n <= l1 < 2 * n else l1 + 2 * n)
                l1 += 1
            if r1 & 1:
                v = r1 - 1
                source.append(v if n <= v < 2 * n else v + 2 * n)
            l1 >>= 1
            r1 >>= 1
        target = []
        l2 += n
        r2 += n
        while l2 < r2:
            if l2 & 1:
                target.append(l2)
                l2 += 1
            if r2 & 1:
                target.append(r2 - 1)
            l2 >>= 1
            r2 >>= 1
        for s in source:
            e[s].append((N3 + q, 0))
        for t in target:
            e[N3 + q + 1].append((t, 0))
        e[N3 + q].append((N3 + q + 1, c))
        self.q = q + 2
