"""
WaveletMatrix: 静的な非負整数列に対する rank / select / 区間 k 番目などを O(log σ) で答える

値は 0 <= a < 2^log の整数であること(log はコンストラクタで指定、既定 32)。

参考: https://judge.yosupo.jp/submission/33990

使用例:
    wm = WaveletMatrix(log=4)
    wm.build([5, 1, 3, 1, 7])
    wm.access(2)              # -> 3
    wm.rank(1, 4)             # -> 2   (A[0:4] 中の 1 の個数)
    wm.select(1, 1)           # -> 3   (2個目(0-index で1)の 1 の位置)
    wm.quantile(0, 5, 2)      # -> 3   (A[0:5] の中で 0-index 2 番目に小さい値)
    wm.range_freq(0, 5, 2, 6) # -> 2   (A[0:5] の中で 2 <= a < 6 の個数)
"""
from itertools import accumulate


class FullyIndexableDictionary():
    """
    完備辞書(ビット列に対する rank / select)

    set で 1 を立ててから build を呼ぶ。build 後は rank が O(1)。
    内部では位置ごとの累積和 cum(cum[i] = [0, i) にある 1 の個数)を持つ。
    """

    def __init__(self, size):
        """
        size: ビット列の長さ

        計算量: O(size)
        """
        self.size = size
        self.block = (size + 31) >> 5
        self.bit = bytearray(size)
        self.cum = None

    def popcount(self, x):
        """32bit 整数 x の立っているビット数を返す。O(1)"""
        x = x - ((x >> 1) & 0x55555555)
        x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
        x = x + (x >> 4) & 0x0f0f0f0f
        x = x + (x >> 8)
        x = x + (x >> 16)
        return x & 0x0000007f

    def set(self, k):
        """k 番目のビットを 1 にする(build 前に呼ぶ)。O(1)"""
        self.bit[k] = 1

    def build(self):
        """rank 用の累積和を構築する。O(size)"""
        self.cum = list(accumulate(self.bit, initial=0))

    def access(self, k):
        """k 番目のビットを返す。O(1)"""
        return self.bit[k]

    def rank(self, k, v):
        """
        [0, k) にあるビット v (0 or 1) の個数を返す

        計算量: O(1)
        """
        r = self.cum[k]
        return r if v else k - r

    def select(self, k, v):
        """
        k 番目(0-index)のビット v の位置を返す(存在しなければ -1)

        計算量: O(log size)
        """
        cum = self.cum
        if k < 0 or self.rank(self.size, v) <= k:
            return -1
        # rank(m, v) >= k + 1 となる最小の m を探し、その m - 1 が答え
        l, r = 0, self.size
        if v:
            while r - l > 1:
                m = (l + r) >> 1
                if cum[m] >= k + 1:
                    r = m
                else:
                    l = m
        else:
            while r - l > 1:
                m = (l + r) >> 1
                if m - cum[m] >= k + 1:
                    r = m
                else:
                    l = m
        return l


class WaveletMatrix():
    """
    ウェーブレット行列

    mat[lv]: 上から lv ビット目で並べ替えた段のビット列(FullyIndexableDictionary)
    mid[lv]: その段の 0 の個数
    """

    def __init__(self, log=32):
        """
        log: 値のビット数(0 <= a < 2^log)
        """
        self.size = 0
        self.log = log
        self.mat = [None] * log
        self.mid = [None] * log

    def build(self, arr):
        """
        数列 arr から構築する

        計算量: O(N log σ)
        """
        n = self.size = len(arr)
        arr = list(arr)
        for lv in range(self.log - 1, -1, -1):
            fid = self.mat[lv] = FullyIndexableDictionary(n + 1)
            bit = fid.bit
            lt = []
            rt = []
            for i, a in enumerate(arr):
                if (a >> lv) & 1:
                    bit[i] = 1
                    rt.append(a)
                else:
                    lt.append(a)
            self.mid[lv] = len(lt)
            fid.build()
            lt.extend(rt)
            arr = lt

    def access(self, k):
        """
        k 番目(0-index)の値を返す

        計算量: O(log σ)
        """
        mat, mid = self.mat, self.mid
        res = 0
        for lv in range(self.log - 1, -1, -1):
            fid = mat[lv]
            cum = fid.cum
            if fid.bit[k]:
                res |= 1 << lv
                k = cum[k] + mid[lv]
            else:
                k -= cum[k]
        return res

    def rank(self, x, r):
        """
        [0, r) にある x の個数を返す

        計算量: O(log σ)
        """
        mat, mid = self.mat, self.mid
        l = 0
        for lv in range(self.log - 1, -1, -1):
            cum = mat[lv].cum
            if (x >> lv) & 1:
                l = cum[l] + mid[lv]
                r = cum[r] + mid[lv]
            else:
                l -= cum[l]
                r -= cum[r]
        return r - l

    def select(self, x, k):
        """
        k 番目(0-index)の x の位置を返す(存在しなければ -1)

        計算量: O(log σ log N)
        """
        if k < 0 or self.rank(x, self.size) <= k:
            return -1
        mat, mid = self.mat, self.mid
        res = 0
        for lv in range(self.log - 1, -1, -1):
            cum = mat[lv].cum
            if (x >> lv) & 1:
                res = mid[lv] + cum[res]
            else:
                res -= cum[res]
        res += k
        for lv in range(self.log):
            if (x >> lv) & 1:
                res = mat[lv].select(res - mid[lv], 1)
            else:
                res = mat[lv].select(res, 0)
        return res

    def quantile(self, l, r, k):
        """
        区間 [l, r) の中で k 番目(0-index)に小さい値を返す

        0 <= k < r - l であること。計算量: O(log σ)
        """
        mat, mid = self.mat, self.mid
        res = 0
        for lv in range(self.log - 1, -1, -1):
            cum = mat[lv].cum
            cl, cr = cum[l], cum[r]
            cnt = (r - cr) - (l - cl)
            if cnt <= k:
                res |= 1 << lv
                k -= cnt
                l = cl + mid[lv]
                r = cr + mid[lv]
            else:
                l -= cl
                r -= cr
        return res

    def rank_less(self, l, r, x):
        """
        区間 [l, r) の中で x 未満の値の個数を返す

        計算量: O(log σ)
        """
        if x <= 0:
            return 0
        if x >> self.log:
            return r - l
        mat, mid = self.mat, self.mid
        res = 0
        for lv in range(self.log - 1, -1, -1):
            cum = mat[lv].cum
            cl, cr = cum[l], cum[r]
            if (x >> lv) & 1:
                res += (r - cr) - (l - cl)
                l = cl + mid[lv]
                r = cr + mid[lv]
            else:
                l -= cl
                r -= cr
        return res

    def range_freq(self, l, r, lower, upper):
        """
        区間 [l, r) の中で lower <= a < upper を満たす値の個数を返す

        計算量: O(log σ)
        """
        if lower >= upper:
            return 0
        return self.rank_less(l, r, upper) - self.rank_less(l, r, lower)


if __name__ == "__main__":
    wm = WaveletMatrix(log=4)
    wm.build([5, 1, 3, 1, 7])
    assert wm.access(2) == 3
    assert wm.rank(1, 4) == 2
    assert wm.select(1, 1) == 3
    assert wm.select(1, 2) == -1
    assert wm.quantile(0, 5, 2) == 3
    assert wm.range_freq(0, 5, 2, 6) == 2
