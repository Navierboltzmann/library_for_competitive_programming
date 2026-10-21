"""文字集合が小さいときの配列 Trie。

デフォルトは英小文字 26 文字。出現回数と prefix 回数を持つ。

使用例:
    t = Trie()
    t.insert("abc")
    t.insert("ab")
    t.count("ab")         # 1
    t.count_prefix("ab")  # 2
    t.find("abc")         # True
"""


class Trie:
    """配列実装の Trie。"""

    def __init__(self, sigma=26, offset=ord("a")):
        """空の Trie を作る。

        sigma: 文字種
        offset: 文字を 0..sigma-1 に写す基準
        """
        self.sigma = sigma
        self.offset = offset
        self.nxt = [[-1] * sigma]
        self.end = [0]
        self.pcnt = [0]

    def _idx(self, c):
        """文字 c を添字に変換する。"""
        if isinstance(c, str):
            return ord(c) - self.offset
        return c - self.offset

    def insert(self, s):
        """文字列 s を挿入する。計算量 O(|s|)。"""
        v = 0
        self.pcnt[v] += 1
        for ch in s:
            i = self._idx(ch)
            if self.nxt[v][i] == -1:
                self.nxt[v][i] = len(self.nxt)
                self.nxt.append([-1] * self.sigma)
                self.end.append(0)
                self.pcnt.append(0)
            v = self.nxt[v][i]
            self.pcnt[v] += 1
        self.end[v] += 1

    def _walk(self, s):
        """s を辿った先の頂点。無ければ -1。"""
        v = 0
        for ch in s:
            i = self._idx(ch)
            v = self.nxt[v][i]
            if v == -1:
                return -1
        return v

    def find(self, s):
        """s が登録されているか。"""
        v = self._walk(s)
        return v != -1 and self.end[v] > 0

    def count(self, s):
        """s と完全一致する挿入回数。"""
        v = self._walk(s)
        return 0 if v == -1 else self.end[v]

    def count_prefix(self, s):
        """s を接頭辞に持つ挿入回数。"""
        v = self._walk(s)
        return 0 if v == -1 else self.pcnt[v]
