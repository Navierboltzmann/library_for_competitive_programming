"""
HL 蛻?隗｣ (Heavy-Light Decomposition)

鬆らせ 0 繧呈?ｹ縺ｨ縺吶ｋ譛ｨ繧? HL 蛻?隗｣縺?, 鬆らせ縺ｫ譁ｰ縺励＞ ID (0..N-1) 繧呈険繧顔峩縺?.
蜷後§ heavy path 荳翫?ｮ鬆らせ縺ｯ ID 縺碁｣邯壹＠, 驛ｨ蛻?譛ｨ繧? ID 縺ｮ騾｣邯壼玄髢薙↓縺ｪ繧九?ｮ縺ｧ,
縲栗D 鬆?縺ｫ荳ｦ縺ｹ譖ｿ縺医◆驟榊?励阪ｒ繧ｻ繧ｰ繝｡繝ｳ繝域惠縺ｪ縺ｩ縺ｫ霈峨○繧後?ｰ繝代せ / 驛ｨ蛻?譛ｨ繧ｯ繧ｨ繝ｪ縺悟?ｦ逅?縺ｧ縺阪ｋ.

菴ｿ逕ｨ萓?:
    e = [[1, 2], [0, 3, 4], [0], [1], [1]]   # 髫｣謗･繝ｪ繧ｹ繝? (辟｡蜷?, 0-indexed)
    hld = HL_Decomposition(e)
    A = [10, 20, 30, 40, 50]                 # 蜷?鬆らせ縺ｮ蛟､
    B = [0] * len(e)
    for v in range(len(e)):
        B[hld.ID[v]] = A[v]                  # ID 鬆?縺ｫ荳ｦ縺ｹ譖ｿ縺医※ SegTree 縺ｫ霈峨○繧?
    for l, r in hld.query(3, 2):             # 繝代せ 3-2 荳翫?ｮ鬆らせ = B 縺ｮ蜊企幕蛹ｺ髢? [l, r) 縺ｮ蜥碁寔蜷?
        ...                                  # seg.prod(l, r) 縺ｪ縺ｩ
    l, r = hld.subtree(1)                    # 鬆らせ 1 縺ｮ驛ｨ蛻?譛ｨ = B[l:r]
    hld.lca(3, 4)                            # -> 1
"""


class HL_Decomposition():

    def __init__(self, e):
        """
        鬆らせ 0 繧呈?ｹ縺ｨ縺励※ HL 蛻?隗｣縺吶ｋ.

        e: 譛ｨ縺ｮ髫｣謗･繝ｪ繧ｹ繝? (辟｡蜷?, 0-indexed)

        螻樊ｧ (PAR, HEAD, DEPTH 縺ｯ譁ｰ縺励＞ ID 縺ｧ豺ｻ蟄嶺ｻ倥￠繝ｻ蛟､繧よ眠縺励＞ ID):
            ID[v]: 鬆らせ v 縺ｮ譁ｰ縺励＞ ID
            PAR[i]: ID i 縺ｮ鬆らせ縺ｮ隕ｪ縺ｮ ID (譬ｹ縺ｯ -1)
            HEAD[i]: ID i 縺ｮ鬆らせ縺悟ｱ槭☆繧? heavy path 縺ｮ蜈磯?ｭ縺ｮ ID
            DEPTH[i]: ID i 縺ｮ鬆らせ縺ｮ light edge 縺ｮ豺ｱ縺? (譬ｹ縺九ｉ霎ｿ繧? light edge 縺ｮ譛ｬ謨ｰ)
            INV[i]: ID i 縺ｮ蜈?縺ｮ鬆らせ逡ｪ蜿ｷ
            SUB[v]: 鬆らせ v 縺ｮ驛ｨ蛻?譛ｨ繧ｵ繧､繧ｺ
        險育ｮ鈴㍼: O(N log N) (蟄舌?ｮ荳ｦ縺ｹ譖ｿ縺?)
        """
        N = len(e)
        self.N = N
        self.e = e

        # BFS 縺ｧ隕ｪ縺ｨ險ｪ蝠城??繧呈ｱゅａ, 驛ｨ蛻?譛ｨ繧ｵ繧､繧ｺ繧定ｨ育ｮ?
        par = [-1] * N
        order = [0]
        par[0] = 0
        for x in order:
            for y in e[x]:
                if par[y] == -1:
                    par[y] = x
                    order.append(y)
        par[0] = -1
        sub = [1] * N
        for i in range(N - 1, 0, -1):
            x = order[i]
            sub[par[x]] += sub[x]

        # DFS (heavy child 繧貞?医↓險ｪ蝠?) 縺ｧ ID 繧呈険繧?
        ID = [-1] * N
        head = [0] * N
        depth = [0] * N
        stack = [0]
        cnt = 0
        while stack:
            x = stack.pop()
            ID[x] = cnt
            cnt += 1
            px = par[x]
            # (sub, 鬆らせ逡ｪ蜿ｷ) 縺ｮ譏?鬆?縺ｫ遨阪? -> 譛螟ｧ縺ｮ繧ゅ?ｮ (heavy child) 縺悟?医↓ pop 縺輔ｌ繧?
            ch = [sub[y] * N + y for y in e[x] if y != px]
            if not ch:
                continue
            ch.sort()
            hx = head[x]
            dx = depth[x] + 1
            for k in range(len(ch) - 1):
                y = ch[k] % N
                head[y] = y
                depth[y] = dx
                stack.append(y)
            y = ch[-1] % N
            head[y] = hx
            depth[y] = dx - 1
            stack.append(y)

        PAR = [-1] * N
        HEAD = [0] * N
        DEPTH = [0] * N
        INV = [0] * N
        for v in range(N):
            i = ID[v]
            INV[i] = v
            if v:
                PAR[i] = ID[par[v]]
            HEAD[i] = ID[head[v]]
            DEPTH[i] = depth[v]

        self.ID = ID
        self.PAR = PAR
        self.HEAD = HEAD
        self.DEPTH = DEPTH
        self.INV = INV
        self.SUB = sub

    def query(self, l, r):
        """
        鬆らせ l, r 髢薙?ｮ繝代せ荳翫?ｮ鬆らせ繧?, 譁ｰ縺励＞ ID 縺ｮ蜊企幕蛹ｺ髢薙?ｮ繝ｪ繧ｹ繝医〒霑斐☆.

        l, r: 蜈?縺ｮ鬆らせ逡ｪ蜿ｷ
        謌ｻ繧雁､: [(a, b), ...] (蜷?蛹ｺ髢? [a, b) 縺ｯ a < b, 鬆?荳榊酔)
        險育ｮ鈴㍼: O(log N)
        """
        HEAD = self.HEAD
        PAR = self.PAR
        DEPTH = self.DEPTH
        L = self.ID[l]
        R = self.ID[r]
        res = []
        if DEPTH[L] < DEPTH[R]:
            L, R = R, L

        dR = DEPTH[R]
        while DEPTH[L] != dR:
            h = HEAD[L]
            res.append((h, L + 1))
            L = PAR[h]

        while HEAD[L] != HEAD[R]:
            h = HEAD[L]
            res.append((h, L + 1))
            L = PAR[h]
            h = HEAD[R]
            res.append((h, R + 1))
            R = PAR[h]

        if L > R:
            L, R = R, L
        res.append((L, R + 1))
        return res

    def subtree(self, v):
        """鬆らせ v 縺ｮ驛ｨ蛻?譛ｨ縺ｫ螻槭☆繧矩?らせ繧?, 譁ｰ縺励＞ ID 縺ｮ蜊企幕蛹ｺ髢? (a, b) 縺ｧ霑斐☆. O(1)"""
        i = self.ID[v]
        return i, i + self.SUB[v]

    def lca(self, u, v):
        """鬆らせ u, v 縺ｮ譛蟆丞?ｱ騾夂･門?? (蜈?縺ｮ鬆らせ逡ｪ蜿ｷ) 繧定ｿ斐☆. O(log N)"""
        HEAD = self.HEAD
        PAR = self.PAR
        U = self.ID[u]
        V = self.ID[v]
        while HEAD[U] != HEAD[V]:
            if HEAD[U] > HEAD[V]:
                U = PAR[HEAD[U]]
            else:
                V = PAR[HEAD[V]]
        return self.INV[U if U < V else V]


if __name__ == "__main__":
    e = [[1, 2], [0, 3, 4], [0], [1], [1]]
    hld = HL_Decomposition(e)
    print(hld.ID)  # [0, 1, 4, 3, 2]
    print(hld.query(3, 2))  # [(3, 4), (4, 5), (0, 2)]
    print(hld.subtree(1), hld.lca(3, 4))  # (1, 4) 1
