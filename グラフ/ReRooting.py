"""
全方位木DP (ReRooting)

各頂点を根としたときの木DPの値を、全頂点についてまとめて O(N) で求める。
根 0 から BFS で向き付けし、bottom-up → top-down の2パスで計算する。

設定する関数:
    ide_ele             : merge の単位元
    merge(a, b)         : 子たちの値をまとめる演算 (結合的。可換でなくてもよいが子の順に依存する)
    adj_bu(a, size)     : 子をまとめた値 a に頂点自身を加えて「部分木の値」にする (size: 部分木の頂点数)
    adj_td(a, size, size_p) : 親側をまとめた値 a に親頂点を加えて「親側部分木の値」にする
                          (size: 親側部分木の頂点数, size_p: 親の更に親側部分木の頂点数)
    adj_fin(a)          : 全方向をまとめた値 a に頂点自身を加えて最終的な答えにする
                          (根の答えは adj_bu(ME[root], N) で計算されるので、通常 adj_fin(a) = adj_bu(a, N))

使用例 (各頂点から他の全頂点への距離の和):
    N = 3
    G = [[1], [0, 2], [1]]
    ide_ele = 0
    merge = lambda a, b: a + b
    adj_bu = lambda a, size: a + size - 1
    adj_td = lambda a, size, size_p: a + size - 1
    adj_fin = lambda a: a + N - 1
    RR = ReRooting(G, ide_ele, merge, adj_bu, adj_td, adj_fin)
    RR.ans  # [3, 2, 3]
"""


class ReRooting:
    """
    全方位木DP。構築時に計算を行い、結果は self.ans[v] (v を根としたときの値) に入る。O(N)

    属性:
        N       : 頂点数
        G       : 根 0 から向き付けした子のリスト (入力の G は変更しない)
        P[v]    : 親 (根は -1)
        R       : BFS 順
        SIZE[v] : 根 0 のときの部分木サイズ
        ME[v]   : 子たちの値を merge したもの
        TD[v]   : 親側部分木の値
        ans[v]  : v を根としたときの答え
    """

    def __init__(self, G, ide_ele, merge, adj_bu, adj_td, adj_fin):
        """
        G: 木の隣接リスト (0-indexed, 無向)
        ide_ele, merge, adj_bu, adj_td, adj_fin: モジュール docstring 参照
        """
        N = len(G)
        self.N = N
        self.ide_ele = ide_ele
        self.merge = merge
        self.adj_bu = adj_bu
        self.adj_td = adj_td
        self.adj_fin = adj_fin

        P = [-1] * N
        C = [None] * N
        R = [0]
        for i in R:
            p = P[i]
            ch = [a for a in G[i] if a != p]
            C[i] = ch
            for a in ch:
                P[a] = i
            R.extend(ch)
        self.G = C
        self.P = P
        self.R = R

        SIZE = [1] * N
        for k in range(N - 1, 0, -1):
            i = R[k]
            SIZE[P[i]] += SIZE[i]
        self.SIZE = SIZE

        ME = [ide_ele] * N
        ans = [0] * N
        TD = [ide_ele] * N
        for k in range(N - 1, 0, -1):
            i = R[k]
            a = adj_bu(ME[i], SIZE[i])
            ans[i] = a
            p = P[i]
            ME[p] = merge(ME[p], a)
        r = R[0]
        ans[r] = adj_bu(ME[r], SIZE[r])

        for i in R:
            ch = C[i]
            ac = TD[i]
            for j in ch:
                TD[j] = ac
                ac = merge(ac, ans[j])
            ac = ide_ele
            rest_i = N - SIZE[i]
            for j in reversed(ch):
                aj = ans[j]
                TD[j] = adj_td(merge(TD[j], ac), N - SIZE[j], rest_i)
                ac = merge(ac, aj)
                ans[j] = adj_fin(merge(ME[j], TD[j]))
        self.ME = ME
        self.TD = TD
        self.ans = ans


if __name__ == "__main__":
    N = 5
    G = [[1, 2], [0, 3, 4], [0], [1], [1]]
    RR = ReRooting(G, 0, lambda a, b: a + b,
                   lambda a, size: a + size - 1,
                   lambda a, size, size_p: a + size - 1,
                   lambda a: a + N - 1)
    print(RR.ans)  # [6, 5, 9, 8, 8]
