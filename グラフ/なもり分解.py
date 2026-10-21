"""
なもり分解 (N 頂点 N 辺の連結グラフ = 閉路が1つのグラフ を、閉路と閉路に生えた木に分解する)

各連結成分がなもりグラフであれば、複数成分 (なもりグラフの森) でもそのまま使える。

使用例:
    e = [[1, 2], [0, 2], [0, 1, 3], [2]]   # 無向グラフの隣接リスト (0-1-2 が閉路、3 は 2 に生えた木)
    cycle, root = namori_decomposition(e)
    # cycle = [True, True, True, False]  (閉路上の頂点か)
    # root  = [0, 1, 2, 2]               (その頂点が属する木の根 = 閉路上の頂点)
"""


def namori_decomposition(e):
    """
    なもりグラフを閉路と木に分解する O(N)

    e: 無向グラフの隣接リスト (辺 u-v は e[u], e[v] の両方に入れる)
    戻り値: (cycle, root)
        cycle[v]: v が閉路上の頂点なら True
        root[v] : v が属する木の根 (= v から閉路に入って最初の閉路上の頂点。閉路上なら v 自身)
    """
    n = len(e)
    deg = [len(ei) for ei in e]
    cycle = [True] * n
    que = [i for i in range(n) if deg[i] == 1]
    for x in que:
        cycle[x] = False
        for ix in e[x]:
            deg[ix] -= 1
            if deg[ix] == 1:
                que.append(ix)

    root = [-1] * n
    que = [i for i in range(n) if cycle[i]]
    for i in que:
        root[i] = i
    for x in que:
        rx = root[x]
        for ix in e[x]:
            if root[ix] == -1:
                root[ix] = rx
                que.append(ix)
    return cycle, root


if __name__ == "__main__":
    e = [[1, 2], [0, 2], [0, 1, 3], [2]]
    print(namori_decomposition(e))  # ([True, True, True, False], [0, 1, 2, 2])
