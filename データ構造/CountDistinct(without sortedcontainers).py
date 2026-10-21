"""区間種類数。SortedMultiset を内蔵し、sortedcontainers 不要。

使用例:
    cd = CountDistinctElements([1, 2, 1, 3])
    cd.query(0, 4)  # 3
"""
import math
from collections import defaultdict
# https://github.com/tatyam-prime/SortedSet/blob/main/SortedMultiset.py
import math
from bisect import bisect_left, bisect_right
from typing import Generic, Iterable, Iterator, List, Tuple, TypeVar, Optional
T = TypeVar('T')

class SortedMultiset(Generic[T]):
    BUCKET_RATIO = 16
    SPLIT_RATIO = 24
    
    def __init__(self, a: Iterable[T] = []) -> None:
        "Make a new SortedMultiset from iterable. / O(N) if sorted / O(N log N)"
        a = list(a)
        n = self.size = len(a)
        if any(a[i] > a[i + 1] for i in range(n - 1)):
            a.sort()
        num_bucket = int(math.ceil(math.sqrt(n / self.BUCKET_RATIO)))
        self.a = [a[n * i // num_bucket : n * (i + 1) // num_bucket] for i in range(num_bucket)]

    def __iter__(self) -> Iterator[T]:
        """要素を昇順に列挙する。"""
        for i in self.a:
            for j in i: yield j

    def __reversed__(self) -> Iterator[T]:
        """要素を降順に列挙する。"""
        for i in reversed(self.a):
            for j in reversed(i): yield j
    
    def __eq__(self, other) -> bool:
        """昇順の要素列として other と等しいか判定する。"""
        return list(self) == list(other)
    
    def __len__(self) -> int:
        """格納している要素数を返す。"""
        return self.size
    
    def __repr__(self) -> str:
        """デバッグ用の文字列表現を返す。"""
        return "SortedMultiset" + str(self.a)
    
    def __str__(self) -> str:
        """集合風の文字列表現を返す。"""
        s = str(list(self))
        return "{" + s[1 : len(s) - 1] + "}"

    def _position(self, x: T) -> Tuple[List[T], int, int]:
        "return the bucket, index of the bucket and position in which x should be. self must not be empty."
        for i, a in enumerate(self.a):
            if x <= a[-1]: break
        return (a, i, bisect_left(a, x))

    def __contains__(self, x: T) -> bool:
        """x が存在するかを O(√N) で判定する。"""
        if self.size == 0: return False
        a, _, i = self._position(x)
        return i != len(a) and a[i] == x

    def count(self, x: T) -> int:
        "Count the number of x."
        return self.index_right(x) - self.index(x)

    def add(self, x: T) -> None:
        "Add an element. / O(√N)"
        if self.size == 0:
            self.a = [[x]]
            self.size = 1
            return
        a, b, i = self._position(x)
        a.insert(i, x)
        self.size += 1
        if len(a) > len(self.a) * self.SPLIT_RATIO:
            mid = len(a) >> 1
            self.a[b:b+1] = [a[:mid], a[mid:]]
    
    def _pop(self, a: List[T], b: int, i: int) -> T:
        """バケット a の i 番目を削除し、空になったバケットも除去する。"""
        ans = a.pop(i)
        self.size -= 1
        if not a: del self.a[b]
        return ans

    def discard(self, x: T) -> bool:
        "Remove an element and return True if removed. / O(√N)"
        if self.size == 0: return False
        a, b, i = self._position(x)
        if i == len(a) or a[i] != x: return False
        self._pop(a, b, i)
        return True

    def lt(self, x: T) -> Optional[T]:
        "Find the largest element < x, or None if it doesn't exist."
        for a in reversed(self.a):
            if a[0] < x:
                return a[bisect_left(a, x) - 1]

    def le(self, x: T) -> Optional[T]:
        "Find the largest element <= x, or None if it doesn't exist."
        for a in reversed(self.a):
            if a[0] <= x:
                return a[bisect_right(a, x) - 1]

    def gt(self, x: T) -> Optional[T]:
        "Find the smallest element > x, or None if it doesn't exist."
        for a in self.a:
            if a[-1] > x:
                return a[bisect_right(a, x)]

    def ge(self, x: T) -> Optional[T]:
        "Find the smallest element >= x, or None if it doesn't exist."
        for a in self.a:
            if a[-1] >= x:
                return a[bisect_left(a, x)]
    
    def __getitem__(self, i: int) -> T:
        "Return the i-th element."
        if i < 0:
            for a in reversed(self.a):
                i += len(a)
                if i >= 0: return a[i]
        else:
            for a in self.a:
                if i < len(a): return a[i]
                i -= len(a)
        raise IndexError
    
    def pop(self, i: int = -1) -> T:
        "Pop and return the i-th element."
        if i < 0:
            for b, a in enumerate(reversed(self.a)):
                i += len(a)
                if i >= 0: return self._pop(a, ~b, i)
        else:
            for b, a in enumerate(self.a):
                if i < len(a): return self._pop(a, b, i)
                i -= len(a)
        raise IndexError

    def index(self, x: T) -> int:
        "Count the number of elements < x."
        ans = 0
        for a in self.a:
            if a[-1] >= x:
                return ans + bisect_left(a, x)
            ans += len(a)
        return ans

    def index_right(self, x: T) -> int:
        "Count the number of elements <= x."
        ans = 0
        for a in self.a:
            if a[-1] > x:
                return ans + bisect_right(a, x)
            ans += len(a)
        return ans


class CountDistinctElements:

    """静的列の区間種類数。"""
    def __init__(self, _init_val):

        """列 _init_val から構築する。"""
        self.index_list = defaultdict(lambda: SortedMultiset())
        for idx, val in enumerate(_init_val):
            self.index_list[val].add(idx)

        self.N = len(_init_val)
        self.a = _init_val
        self.b = [self.N + 1] * self.N

        for v in self.index_list.values():

            n = len(v)
            for i in range(n - 1):
                idx = v[i]
                jdx = v[i + 1]
                self.b[idx] = jdx
            self.b[v[-1]] = self.N + 1

        self.num = 1 << (self.N - 1).bit_length()
        self.tree = [SortedMultiset() for _ in range(2 * self.num)]

        for i in range(self.N):
            self.tree[self.num + i].add(self.b[i])
            k = i + self.num
            while k > 1:
                self.tree[k >> 1].add(self.b[i])
                k >>= 1

    def query(self, l, r):

        """半開区間 [l, r) の種類数。"""
        _x = 0
        lim = r
        l += self.num
        r += self.num
        while l < r:
            if l & 1:
                _x += len(self.tree[l]) - self.tree[l].index(lim)
                l += 1
            if r & 1:
                _x += len(self.tree[r - 1]) - self.tree[r - 1].index(lim)
            l >>= 1
            r >>= 1

        return _x
