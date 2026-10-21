"""
Auxiliary Tree (虚树 / 圧縮木)

木の頂点の部分集合 S に対し、S と S の任意の2頂点の LCA だけからなる木を作る。
LCA は Euler tour + Sparse Table で O(1)。非再帰なので recursionlimit の変更は不要。

使用例:
    e = [[1, 2], [0, 3, 4], [0], [1], [1]]
    AT = AuxiliaryTree(e)
    AT.lca(3, 4)                       # 1
    AT.build_auxiliary_tree([3, 4, 2]) # {1: [4, 3, 0], 4: [1], ...} のような隣接リスト (dict)
"""


class AuxiliaryTree:
    """
    Auxiliary Tree を構築するためのクラス。構築 O(N logN)。

    属性:
        n                : 頂点数
        adj              : 隣接リスト
        depth[v]         : 根からの深さ
        parent[v]        : 親 (根は -1)
        euler            : Euler tour (長さ 2N-1)
        first_occurrence : euler で v が最初に現れる位置
        sparse_table     : euler 上の Sparse Table。値は (depth << shift) | 頂点番号
    """

    def __init__(self, e, root=0):
        """
        初期化 O(N logN)

        e: 木の隣接リスト (0-indexed)
        root: 根 (省略時 0)
        """
        self.n = len(e)
        self.adj = e
        self.root = root
        self.depth = [0] * self.n
        self.parent = [-1] * self.n
        self.euler = []
        self.first_occurrence = [-1] * self.n
        self.log = [0] * (2 * self.n)
        self.sparse_table = []
        self._shift = max(1, self.n.bit_length())
        self._mask = (1 << self._shift) - 1

        self._initialize_log()
        self._dfs(root)
        self._build_sparse_table()

    def _initialize_log(self):
        """log[i] = floor(log2(i)) のテーブルを作る O(N)"""
        log = self.log
        for i in range(2, len(log)):
            log[i] = log[i >> 1] + 1

    def _dfs(self, root):
        """非再帰 DFS で Euler tour・深さ・親を求める O(N)"""
        adj = self.adj
        depth = self.depth
        parent = self.parent
        first = self.first_occurrence
        euler = self.euler
        it = [0] * self.n
        first[root] = 0
        euler.append(root)
        stack = [root]
        while stack:
            v = stack[-1]
            i = it[v]
            adj_v = adj[v]
            if i < len(adj_v):
                it[v] = i + 1
                u = adj_v[i]
                if u == parent[v]:
                    continue
                parent[u] = v
                depth[u] = depth[v] + 1
                first[u] = len(euler)
                euler.append(u)
                stack.append(u)
            else:
                stack.pop()
                if stack:
                    euler.append(stack[-1])

    def _build_sparse_table(self):
        """Euler tour 上の区間最小 (深さ最小の頂点) を引く Sparse Table を作る O(N logN)"""
        depth = self.depth
        shift = self._shift
        row = [(depth[v] << shift) | v for v in self.euler]
        table = [row]
        m = len(row)
        half = 1
        while 2 * half <= m:
            prev = row
            row = [prev[i] if prev[i] < prev[i + half] else prev[i + half]
                   for i in range(len(prev) - half)]
            table.append(row)
            half *= 2
        self.sparse_table = table

    def lca(self, u, v):
        """
        u と v の最小共通祖先を返す O(1)

        u, v: 頂点番号
        """
        first = self.first_occurrence
        left = first[u]
        right = first[v]
        if left > right:
            left, right = right, left
        j = (right - left + 1).bit_length() - 1
        row = self.sparse_table[j]
        a = row[left]
        b = row[right - (1 << j) + 1]
        return (a if a < b else b) & self._mask

    def build_auxiliary_tree(self, nodes):
        """
        頂点集合 nodes から Auxiliary Tree を構築する O(K logK) (K = len(nodes))

        nodes: 含めたい頂点のリスト (重複は無視)
        戻り値: {頂点: [隣接頂点, ...]} の dict (無向)。nodes が1頂点以下なら空 dict。
        """
        first = self.first_occurrence
        table = self.sparse_table
        mask = self._mask
        nodes = sorted(set(nodes), key=first.__getitem__)
        tree = {}

        def add_edge(a, b):
            """仮想木に無向辺 a-b を追加する。"""
            if a in tree:
                tree[a].append(b)
            else:
                tree[a] = [b]
            if b in tree:
                tree[b].append(a)
            else:
                tree[b] = [a]

        stack = []
        for node in nodes:
            if not stack:
                stack.append(node)
                continue

            left = first[stack[-1]]
            right = first[node]
            j = (right - left + 1).bit_length() - 1
            row = table[j]
            a = row[left]
            b = row[right - (1 << j) + 1]
            lca_node = (a if a < b else b) & mask

            first_lca = first[lca_node]
            while len(stack) > 1 and first[stack[-2]] >= first_lca:
                top = stack.pop()
                add_edge(stack[-1], top)

            if stack[-1] != lca_node:
                top = stack.pop()
                add_edge(lca_node, top)
                stack.append(lca_node)

            stack.append(node)

        while len(stack) > 1:
            top = stack.pop()
            add_edge(stack[-1], top)

        return tree


if __name__ == "__main__":
    e = [[1, 2], [0, 3, 4], [0], [1], [1]]
    AT = AuxiliaryTree(e)
    print(AT.lca(3, 4))  # 1
    print(AT.build_auxiliary_tree([3, 4, 2]))
