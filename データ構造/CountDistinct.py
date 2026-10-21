"""静的な列の半開区間 ``[left, right)`` にある値の種類数を求める。

各位置の「次に同じ値が現れる位置」をマージソート木に載せる。
``sortedcontainers`` などの外部パッケージは不要。

使用例:
    counter = CountDistinctElements([1, 2, 1, 3])
    counter.query(0, 4)  # 3
    counter.query(0, 3)  # 2
"""
from bisect import bisect_left


class CountDistinctElements:
    """前処理 O(N log N)、1 クエリ O(log^2 N) の区間種類数データ構造。"""

    def __init__(self, values):
        """列 values からマージソート木を構築する。"""
        self.n = len(values)
        size = 1 << (self.n - 1).bit_length() if self.n else 1
        self.size = size

        # next_index[i] は values[i] と同じ値が次に現れる位置。
        # 再出現しない場合は n とし、right 以上なら区間内の最後の 1 個と数える。
        next_index = [self.n] * self.n
        latest = {}
        for index in range(self.n - 1, -1, -1):
            value = values[index]
            next_index[index] = latest.get(value, self.n)
            latest[value] = index

        tree = [[] for _ in range(2 * size)]
        for index, nxt in enumerate(next_index):
            tree[size + index] = [nxt]
        for node in range(size - 1, 0, -1):
            tree[node] = self._merge(tree[node << 1], tree[node << 1 | 1])
        self.tree = tree

    @staticmethod
    def _merge(left, right):
        """2 つのソート済み列を線形時間でマージする。"""
        merged = []
        append = merged.append
        i = j = 0
        while i < len(left) and j < len(right):
            if left[i] <= right[j]:
                append(left[i])
                i += 1
            else:
                append(right[j])
                j += 1
        merged.extend(left[i:])
        merged.extend(right[j:])
        return merged

    def query(self, left, right):
        """半開区間 ``[left, right)`` に含まれる異なる値の個数を返す。"""
        if not 0 <= left <= right <= self.n:
            raise ValueError(f"query out of range: ({left}, {right})")

        result = 0
        tree = self.tree
        left += self.size
        right_node = right + self.size
        while left < right_node:
            if left & 1:
                node = tree[left]
                result += len(node) - bisect_left(node, right)
                left += 1
            if right_node & 1:
                right_node -= 1
                node = tree[right_node]
                result += len(node) - bisect_left(node, right)
            left >>= 1
            right_node >>= 1
        return result


if __name__ == "__main__":
    counter = CountDistinctElements([1, 2, 1, 3])
    assert counter.query(0, 4) == 3
    assert counter.query(0, 3) == 2
    assert counter.query(2, 2) == 0
