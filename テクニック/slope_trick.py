"""Slope Trick。

区分線形凸関数 f を左右のヒープで管理し、
f(x) += max(x-a, 0) / max(a-x, 0) / 定数 を O(log N) で追加する。

使用例:
    st = slope_trick()
    st.add_minus(3)   # += max(3-x, 0)
    st.add_plus(5)    # += max(x-5, 0)
    st.Min()          # 最小値
    st.arg_Min()      # 最小を与える x
"""
from heapq import heappush, heappushpop, heappop


class slope_trick:
    """区分線形凸関数の最小値管理。"""

    def __init__(self):
        """f = 0 で初期化する。"""
        self.L = [10**20]
        self.R = [10**20]
        self.m = 0

    def add_plus(self, a):
        """f(x) += max(x - a, 0) を行う。

        計算量: O(log N)
        """
        b = -heappushpop(self.L, -a)
        heappush(self.R, b)
        self.update(b, a, 1)

    def add_minus(self, a):
        """f(x) += max(a - x, 0) を行う。

        計算量: O(log N)
        """
        b = heappushpop(self.R, a)
        heappush(self.L, -b)
        self.update(b, a, -1)

    def add_const(self, a):
        """f(x) += a を行う。"""
        self.m += a

    def update(self, x, y, k):
        """傾き追加に伴う最小値の補正。"""
        if k == 1:
            self.m += max(x - y, 0)
        else:
            self.m += max(y - x, 0)

    def Min(self):
        """現在の最小値を返す。

        計算量: O(1)
        """
        return self.m

    def arg_Min(self):
        """最小を与える x のひとつを返す。

        計算量: O(log N)
        """
        x = -self.L[0]
        return x
