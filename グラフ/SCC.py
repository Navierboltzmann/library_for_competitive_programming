"""
強連結成分分解 (SCC, Kosaraju 法) と縮約グラフの構築。非再帰。

ラベルはトポロジカル順 (ラベル a の成分から b の成分へ辺があれば a < b)。

使用例:
    N = 4
    G = [[1], [2], [0, 3], []]
    RG = [[2], [0], [1], [2]]           # 逆辺。省略すると G から自動で作る
    label, group = scc(N, G, RG)        # label = 2, group = [0, 0, 0, 1]
    G0, GP = construct(N, G, label, group)
    # G0 = [{1}, set()]  (縮約後の辺), GP = [[0, 1, 2], [3]]  (各成分に属する頂点)
"""


def scc(N, G, RG=None):
    """
    強連結成分分解 O(N + M)

    N: 頂点数
    G: 順方向の有向グラフ (隣接リスト)
    RG: 逆方向の有向グラフ (None なら G から作る)
    戻り値: (ラベル数, 各頂点のラベル番号のリスト)
    """
    if RG is None:
        RG = [[] for _ in range(N)]
        for v in range(N):
            for w in G[v]:
                RG[w].append(v)

    order = []
    used = [0] * N
    it = [0] * N
    for i in range(N):
        if used[i]:
            continue
        used[i] = 1
        stack = [i]
        while stack:
            s = stack[-1]
            gs = G[s]
            k = it[s]
            ls = len(gs)
            while k < ls and used[gs[k]]:
                k += 1
            if k < ls:
                t = gs[k]
                it[s] = k + 1
                used[t] = 1
                stack.append(t)
            else:
                it[s] = k
                stack.pop()
                order.append(s)

    group = [None] * N
    label = 0
    for s in reversed(order):
        if group[s] is not None:
            continue
        group[s] = label
        stack = [s]
        while stack:
            v = stack.pop()
            for t in RG[v]:
                if group[t] is None:
                    group[t] = label
                    stack.append(t)
        label += 1
    return label, group


def construct(N, G, label, group):
    """
    SCC の結果から縮約後のグラフを構築する O(N + M)

    N: 頂点数
    G: 順方向の有向グラフ
    label, group: scc の戻り値
    戻り値: (G0, GP)
        G0[a]: 成分 a から辺が出ている成分の集合 (set, 自己ループなし)
        GP[a]: 成分 a に属する頂点のリスト (頂点番号の昇順)
    """
    G0 = [set() for _ in range(label)]
    GP = [[] for _ in range(label)]
    for v in range(N):
        lbs = group[v]
        add = G0[lbs].add
        for w in G[v]:
            lbt = group[w]
            if lbs != lbt:
                add(lbt)
        GP[lbs].append(v)
    return G0, GP


if __name__ == "__main__":
    N = 4
    G = [[1], [2], [0, 3], []]
    label, group = scc(N, G)
    print(label, group)  # 2 [0, 0, 0, 1]
    print(construct(N, G, label, group))  # ([{1}, set()], [[0, 1, 2], [3]])
