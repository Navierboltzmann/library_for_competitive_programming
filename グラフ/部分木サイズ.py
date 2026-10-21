"""
部分木サイズ

隣接リストで与えられた木について, 根からの距離と各頂点を根とする部分木のサイズを求める.

使用例:
    e = [[1, 2], [0, 3, 4], [0], [1], [1]]
    sub, dist = subtree_size(e, root=0)
    # sub  = [5, 3, 1, 1, 1]
    # dist = [0, 1, 1, 2, 2]
"""


def subtree_size(e, root=0):
    """
    各頂点の部分木サイズと根からの距離を求める.

    e: 木の隣接リスト (無向, 0-indexed, 連結であること)
    root: 根
    戻り値: (sub, dist)
        sub[v]: 頂点 v を根とする部分木の頂点数
        dist[v]: 根から頂点 v までの距離 (辺数)
    計算量: O(N)
    """
    N = len(e)
    par = [-1] * N
    dist = [-1] * N
    dist[root] = 0
    order = [root]
    for x in order:
        dx = dist[x] + 1
        for y in e[x]:
            if dist[y] == -1:
                dist[y] = dx
                par[y] = x
                order.append(y)

    sub = [1] * N
    for i in range(len(order) - 1, 0, -1):
        x = order[i]
        sub[par[x]] += sub[x]

    return sub, dist


if __name__ == "__main__":
    print(subtree_size([[1, 2], [0, 3, 4], [0], [1], [1]]))  # ([5, 3, 1, 1, 1], [0, 1, 1, 2, 2])
