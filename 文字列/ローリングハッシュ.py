"""複数 mod の Rolling Hash。

S は整数列を想定する。部分列 S[l:r] のハッシュを O(1) で取る。

使用例:
    rh = Rolling_hash([ord(c) for c in 'abcab'])
    rh.query(0, 0, 2) == rh.query(0, 3, 5)  # True ('ab' == 'ab')
"""


class Rolling_hash:
    """複数ハッシュによる部分列比較。"""

    def __init__(self, S, mod_list=None, base_num_list=None):
        """ハッシュテーブルを構築する。

        S: 整数列
        mod_list: 法のリスト
        base_num_list: 基数のリスト
        計算量: O(K N)  (K はハッシュの本数)
        """
        if mod_list is None:
            mod_list = [1000000009, 998244353, 3000012541]
        if base_num_list is None:
            base_num_list = [500009, 500029, 500041]
        n = len(mod_list)
        N = len(S)
        self.n = n
        self.S = S
        self.N = N
        self.mod_list = mod_list
        H = [[0] * (N + 1) for _ in range(n)]
        b_inv = [[1] * N for _ in range(n)]
        for j in range(n):
            base_num = base_num_list[j]
            mod = mod_list[j]
            base_inv = pow(base_num, mod - 2, mod)
            tb = 1
            tH = 0
            Hj = H[j]
            invj = b_inv[j]
            tb_inv = 1
            for i in range(N):
                tH = (tH + S[i] * tb) % mod
                Hj[i + 1] = tH
                if i + 1 < N:
                    tb = tb * base_num % mod
                    tb_inv = tb_inv * base_inv % mod
                    invj[i + 1] = tb_inv
        self.H = H
        self.b_inv = b_inv

    def query(self, j, l, r):
        """j 番目のハッシュで S[l:r] のハッシュを返す。

        j: ハッシュ番号
        l, r: 半開区間
        計算量: O(1)
        """
        return (self.H[j][r] - self.H[j][l]) * self.b_inv[j][l] % self.mod_list[j]
