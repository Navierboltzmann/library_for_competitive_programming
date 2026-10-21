"""Mo's algorithm で区間内の異なる値の個数を求める。

クエリは 0-index の半開区間 ``[left, right)`` で指定する。
値を座標圧縮するため、負数や大きな整数を含んでもメモリを浪費しない。

使用例:
    values = [1, 2, 1, 3, 2]
    mo_distinct(values, [(0, 3), (1, 5), (2, 2)])  # [2, 3, 0]
"""
from math import isqrt


def mo_distinct(values, queries):
    """各半開区間 ``[left, right)`` に含まれる異なる値の個数を返す。

    values: 対象の列
    queries: ``(left, right)`` の列
    戻り値: 入力クエリと同じ順番の答え
    計算量: おおむね O((N + Q) √N)
    """
    n = len(values)
    q = len(queries)
    if q == 0:
        return []

    normalized_queries = []
    for index, (left, right) in enumerate(queries):
        if not 0 <= left <= right <= n:
            raise ValueError(f"query out of range: ({left}, {right})")
        normalized_queries.append((left, right, index))

    # 値域ではなく異なる値の数に比例する配列だけを持つ。
    compressed = {}
    encoded = []
    for value in values:
        value_id = compressed.setdefault(value, len(compressed))
        encoded.append(value_id)
    frequency = [0] * len(compressed)

    # N / sqrt(Q) は N と Q の大きさが偏った場合にも移動量を抑えやすい。
    block_size = max(1, n // max(1, isqrt(q)))

    def query_order(query):
        """左端ブロックと、ブロックごとに反転した右端順を返す。"""
        left, right, _ = query
        block = left // block_size
        return block, right if block & 1 == 0 else -right

    normalized_queries.sort(key=query_order)
    answers = [0] * q
    current_left = current_right = distinct = 0

    def add(position):
        """position の要素を現在区間に追加する。"""
        nonlocal distinct
        value_id = encoded[position]
        if frequency[value_id] == 0:
            distinct += 1
        frequency[value_id] += 1

    def remove(position):
        """position の要素を現在区間から除外する。"""
        nonlocal distinct
        value_id = encoded[position]
        frequency[value_id] -= 1
        if frequency[value_id] == 0:
            distinct -= 1

    for left, right, index in normalized_queries:
        while current_left > left:
            current_left -= 1
            add(current_left)
        while current_right < right:
            add(current_right)
            current_right += 1
        while current_left < left:
            remove(current_left)
            current_left += 1
        while current_right > right:
            current_right -= 1
            remove(current_right)
        answers[index] = distinct

    return answers


if __name__ == "__main__":
    assert mo_distinct([1, 2, 1, 3, 2], [(0, 3), (1, 5), (2, 2)]) == [2, 3, 0]
