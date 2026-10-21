"""
BinaryTrie: 非負整数の(多重)集合を2進トライで管理する

追加・削除・個数・k番目・lower_bound(bisect)が O(max_log)。
xor_all で全要素に一括で xor を作用できる(xor 最小値クエリなどに使う)。
値は 0 <= x < 2^max_log であること。

参考:
 - https://atcoder.jp/contests/arc028/submissions/19916627
 - https://judge.yosupo.jp/submission/35057

使用例:
    bt = BinaryTrie(max_log=30)
    for x in [5, 1, 3, 3]:
        bt.add(x)
    len(bt)               # -> 4
    bt.kth_elem(1)        # -> 3   (0-index で 1 番目に小さい値)
    bt.bisect_left(3)     # -> 1
    bt.count(3)           # -> 2
    bt.discard(3)         # 3 を1個削除
    bt.xor_all(6)         # 全要素に xor 6 -> {3, 5, 7}
    bt.minimum()          # -> 3   (= min(a ^ 6))
    list(bt)              # -> [3, 5, 7]
"""


class BinaryTrie:
    """
    2進トライによる非負整数の多重集合

    ノード v の子は edges[2v](ビット0), edges[2v+1](ビット1)。存在しなければ -1。
    size[v]: v の部分木の要素数, is_end[v]: v で終わる要素の個数。
    lazy: 全体に作用させている xor の値(格納値 ^ lazy が実際の値)。
    """

    def __init__(
        self,
        max_log: int = 35,
        allow_multiple_elements: bool = True,
        add_query_limit: int = 3 * 10**5,
    ):
        """
        max_log: 値のビット数(0 <= x < 2^max_log)
        allow_multiple_elements: False なら同じ値は1個までしか持たない(set として使う)
        add_query_limit: 互換性のために残している引数(ノードは必要に応じて確保するので上限はない)
        """
        self.max_log = max_log
        self.x_end = 1 << max_log
        self.v_list = [0] * (max_log + 1)
        self.multiset = allow_multiple_elements
        self.add_query_count = 0
        self.add_query_limit = add_query_limit
        self.edges = [-1, -1]
        self.size = [0]
        self.is_end = [0]
        self.max_v = 0
        self.lazy = 0

    def _find_node(self, x: int) -> int:
        """値 x(lazy 適用前の実際の値)の終端ノードを返す。無ければ -1。O(max_log)"""
        x ^= self.lazy
        edges = self.edges
        v = 0
        for i in range(self.max_log - 1, -1, -1):
            v = edges[2 * v + ((x >> i) & 1)]
            if v == -1:
                return -1
        return v

    def add(self, x: int):
        """
        x を1個追加する(set モードで既にあれば何もしない)

        計算量: O(max_log)
        """
        x ^= self.lazy
        edges, size, is_end = self.edges, self.size, self.is_end
        v_list = self.v_list
        v = 0
        for i in range(self.max_log - 1, -1, -1):
            j = 2 * v + ((x >> i) & 1)
            nv = edges[j]
            if nv == -1:
                nv = self.max_v = self.max_v + 1
                edges[j] = nv
                edges.append(-1)
                edges.append(-1)
                size.append(0)
                is_end.append(0)
            v = nv
            v_list[i] = v
        if self.multiset or is_end[v] == 0:
            is_end[v] += 1
            for u in v_list:
                size[u] += 1
        self.add_query_count += 1

    def discard(self, x: int):
        """
        x を1個削除する(無ければ何もしない)

        計算量: O(max_log)
        """
        if not 0 <= x < self.x_end:
            return
        x ^= self.lazy
        edges, size, is_end = self.edges, self.size, self.is_end
        v_list = self.v_list
        v = 0
        for i in range(self.max_log - 1, -1, -1):
            v = edges[2 * v + ((x >> i) & 1)]
            if v == -1:
                return
            v_list[i] = v
        if is_end[v] > 0:
            is_end[v] -= 1
            for u in v_list:
                size[u] -= 1

    def erase(self, x: int, count: int = -1):
        """
        x を count 個削除する(count = -1 なら全て削除。持っている個数より多ければ全て)

        計算量: O(max_log)
        """
        if not 0 <= x < self.x_end:
            return
        x ^= self.lazy
        edges, size, is_end = self.edges, self.size, self.is_end
        v_list = self.v_list
        v = 0
        for i in range(self.max_log - 1, -1, -1):
            v = edges[2 * v + ((x >> i) & 1)]
            if v == -1:
                return
            v_list[i] = v
        if count == -1 or is_end[v] < count:
            count = is_end[v]
        if is_end[v] > 0:
            is_end[v] -= count
            for u in v_list:
                size[u] -= count

    def count(self, x: int) -> int:
        """
        x の個数を返す

        計算量: O(max_log)
        """
        if not 0 <= x < self.x_end:
            return 0
        v = self._find_node(x)
        return 0 if v == -1 else self.is_end[v]

    def bisect_left(self, x: int) -> int:
        """
        x 未満の要素数を返す

        計算量: O(max_log)
        """
        if x < 0:
            return 0
        if self.x_end <= x:
            return len(self)
        edges, size, lazy = self.edges, self.size, self.lazy
        v = 0
        ret = 0
        for i in range(self.max_log - 1, -1, -1):
            lc = edges[2 * v]
            rc = edges[2 * v + 1]
            if (lazy >> i) & 1:
                lc, rc = rc, lc
            if (x >> i) & 1:
                if lc != -1:
                    ret += size[lc]
                if rc == -1:
                    return ret
                v = rc
            else:
                if lc == -1:
                    return ret
                v = lc
        return ret

    def bisect_right(self, x: int) -> int:
        """
        x 以下の要素数を返す

        計算量: O(max_log)
        """
        return self.bisect_left(x + 1)

    def index(self, x: int) -> int:
        """
        x が何番目(0-index)に小さいかを返す。無ければ ValueError

        計算量: O(max_log)
        """
        if x not in self:
            raise ValueError(f"{x} is not in BinaryTrie")
        return self.bisect_left(x)

    def find(self, x: int) -> int:
        """
        x が何番目(0-index)に小さいかを返す。無ければ -1

        計算量: O(max_log)
        """
        if x not in self:
            return -1
        return self.bisect_left(x)

    def kth_elem(self, k: int) -> int:
        """
        k 番目(0-index)に小さい値を返す。負の k は末尾から(-1 で最大値)

        -len <= k < len であること。計算量: O(max_log)
        """
        edges, size, lazy = self.edges, self.size, self.lazy
        if k < 0:
            k += size[0]
        v = 0
        ret = 0
        for i in range(self.max_log - 1, -1, -1):
            lc = edges[2 * v]
            rc = edges[2 * v + 1]
            if (lazy >> i) & 1:
                lc, rc = rc, lc
            if lc == -1:
                v = rc
                ret |= 1 << i
            elif size[lc] <= k:
                k -= size[lc]
                v = rc
                ret |= 1 << i
            else:
                v = lc
        return ret

    def minimum(self) -> int:
        """最小値を返す(空でないこと)。O(max_log)"""
        return self.kth_elem(0)

    def maximum(self) -> int:
        """最大値を返す(空でないこと)。O(max_log)"""
        return self.kth_elem(-1)

    def xor_all(self, x: int):
        """
        全要素に xor x を作用させる

        計算量: O(1)
        """
        self.lazy ^= x

    def __iter__(self):
        """全要素を昇順に(重複込みで)列挙する。O(ノード数)"""
        edges, lazy = self.edges, self.lazy
        q = [(0, 0)]
        for i in range(self.max_log - 1, -1, -1):
            swap = (lazy >> i) & 1
            nq = []
            for v, x in q:
                lc = edges[2 * v]
                rc = edges[2 * v + 1]
                if swap:
                    lc, rc = rc, lc
                if lc != -1:
                    nq.append((lc, 2 * x))
                if rc != -1:
                    nq.append((rc, 2 * x + 1))
            q = nq
        is_end = self.is_end
        for v, x in q:
            for _ in range(is_end[v]):
                yield x

    def __str__(self):
        """格納要素を昇順に並べた文字列表現を返す。"""
        return "BinaryTrie(" + ", ".join(map(str, self)) + ")"

    def __getitem__(self, k):
        """kth_elem(k) と同じ"""
        return self.kth_elem(k)

    def __contains__(self, x: int) -> bool:
        """x が 1 個以上格納されているかを返す。"""
        return bool(self.count(x))

    def __len__(self):
        """重複を含む格納要素数を返す。"""
        return self.size[0]

    def __bool__(self):
        """要素が 1 個以上あるかを返す。"""
        return bool(len(self))

    def __ixor__(self, x: int):
        """bt ^= x で xor_all(x)"""
        self.xor_all(x)
        return self


if __name__ == "__main__":
    bt = BinaryTrie(max_log=30)
    for x in [5, 1, 3, 3]:
        bt.add(x)
    assert len(bt) == 4 and bt.kth_elem(1) == 3
    assert bt.bisect_left(3) == 1 and bt.count(3) == 2
    bt.discard(3)
    bt.xor_all(6)
    assert bt.minimum() == 3 and list(bt) == [3, 5, 7]
    print(bt)  # BinaryTrie(3, 5, 7)
