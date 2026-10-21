"""
SumCntBIT: 多重集合に対する「k番目の値」「上位/下位k個の和」を O(log N) で求めるデータ構造

個数 BIT と総和 BIT を座標圧縮した上で持つ。
追加・削除される可能性のある値を全て先にコンストラクタへ渡しておく必要がある(座圧用)。
ABC306E, ABC287G, ABC281E, ABC062D などで利用可能。

参考: https://github.com/cozysauna/competitive-programming-python/blob/master/graph/sum_cnt_bit.py

使用例:
    bit = SumCntBIT([5, 1, 3, 3])
    bit.add(5, 1); bit.add(1, 1); bit.add(3, 2)   # {1, 3, 3, 5}
    bit.kth_smallest_value(2)  # -> 3   (k は 1-index)
    bit.kth_biggest_value(1)   # -> 5
    bit.kth_smallest_sum(3)    # -> 1 + 3 + 3 = 7
    bit.kth_biggest_sum(2)     # -> 5 + 3 = 8
    bit.remove(3, 1)           # {1, 3, 5}
"""


class SumCntBIT:
    """
    座圧済みの個数 BIT と総和 BIT で多重集合を管理する

    add_items: add, remove される値のリスト(重複可、座圧用)
    """

    def __init__(self, add_items):
        """
        add_items: add, remove される可能性のある値の列

        計算量: O(N log N)
        """
        add_items = sorted(set(add_items))
        self.N = len(add_items)
        self.compress = {x: i for i, x in enumerate(add_items, 1)}
        # original[i]: 圧縮後の番号 i (1-index) に対応する元の値
        self.original = [0] + add_items
        self.sum_data = [0] * (self.N + 1)
        self.cnt_data = [0] * (self.N + 1)
        self.all_cnt = 0
        self.all_sum = 0

    def _add(self, data, i, x):
        """BIT data の i 番目(1-index)に x を加算する。O(log N)"""
        n = self.N
        while i <= n:
            data[i] += x
            i += i & -i

    def add(self, x, y):
        """
        値 x を y 個追加する

        x: 追加する値(コンストラクタで渡した値のいずれか)
        y: 個数(負なら削除)
        計算量: O(log N)
        """
        s = x * y
        self.all_sum += s
        self.all_cnt += y
        sum_data, cnt_data, n = self.sum_data, self.cnt_data, self.N
        i = self.compress[x]
        while i <= n:
            sum_data[i] += s
            cnt_data[i] += y
            i += i & -i

    def remove(self, x, y):
        """
        値 x を y 個削除する

        計算量: O(log N)
        """
        self.add(x, -y)

    def _lower_bound(self, data, w):
        """
        data の累積和が w 以上になる最小の index (1-index) を返す

        w <= 0 なら 1 を返す。計算量: O(log N)
        """
        if w <= 0:
            return 1
        n = self.N
        i = 0
        size = 1 << n.bit_length()
        while size:
            j = i + size
            if j <= n and data[j] < w:
                w -= data[j]
                i = j
            size >>= 1
        return i + 1

    def _sum(self, data, i):
        """BIT data の [1, i] の和を返す。O(log N)"""
        ret = 0
        while i > 0:
            ret += data[i]
            i -= i & -i
        return ret

    def kth_smallest_value(self, k):
        """
        k 番目(1-index)に小さい値を返す

        1 <= k <= 全要素数 であること。計算量: O(log N)
        """
        return self.original[self._lower_bound(self.cnt_data, k)]

    def kth_biggest_value(self, k):
        """
        k 番目(1-index)に大きい値を返す

        1 <= k <= 全要素数 であること。計算量: O(log N)
        """
        return self.kth_smallest_value(self.all_cnt - k + 1)

    def kth_smallest_sum(self, k):
        """
        小さい方から k 個の値の和を返す

        k が全要素数を超える場合は全要素の和。計算量: O(log N)
        """
        k = min(k, self.all_cnt)
        if k <= 0:
            return 0
        # 累積個数が k 未満となる最大の位置 i まで降りつつ、個数と和を同時に集計する
        cnt_data, sum_data, n = self.cnt_data, self.sum_data, self.N
        i = 0
        cnt = 0
        total = 0
        size = 1 << n.bit_length()
        while size:
            j = i + size
            if j <= n and cnt + cnt_data[j] < k:
                cnt += cnt_data[j]
                total += sum_data[j]
                i = j
            size >>= 1
        return total + self.original[i + 1] * (k - cnt)

    def kth_biggest_sum(self, k):
        """
        大きい方から k 個の値の和を返す

        k が全要素数を超える場合は全要素の和。計算量: O(log N)
        """
        k = min(k, self.all_cnt)
        return self.all_sum - self.kth_smallest_sum(self.all_cnt - k)

    def print(self):
        """全要素を昇順に表示する(デバッグ用)。O(N log N)"""
        print(*[self.kth_smallest_value(i) for i in range(1, self.all_cnt + 1)])


if __name__ == "__main__":
    bit = SumCntBIT([5, 1, 3, 3])
    bit.add(5, 1)
    bit.add(1, 1)
    bit.add(3, 2)
    assert bit.kth_smallest_value(2) == 3
    assert bit.kth_biggest_value(1) == 5
    assert bit.kth_smallest_sum(3) == 7
    assert bit.kth_biggest_sum(2) == 8
    bit.remove(3, 1)
    bit.print()  # 1 3 5
