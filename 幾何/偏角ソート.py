"""偏角ソートと、原点まわりの偏角比較。

atan2 を使わず、象限と外積だけで比較する。

使用例:
    pts = [(1, 0), (0, 1), (-1, 0), (1, 1)]
    argument_sort(pts)  # 反時計、x 軸正から
"""


def argument_cmp(p, q):
    """原点から見た偏角の比較。p < q なら負、等しければ 0、大きければ正。

    (0, 0) は最も小さい。同じ偏角なら原点に近い方が小さい。
    """
    def half(x, y):
        # 上半平面と負の x 軸を 1、それ以外(原点含む下側)を 0
        return int(y < 0 or (y == 0 and x < 0))

    ox, oy = p
    qx, qy = q
    if ox == 0 and oy == 0:
        return 0 if qx == 0 and qy == 0 else -1
    if qx == 0 and qy == 0:
        return 1
    hp, hq = half(ox, oy), half(qx, qy)
    if hp != hq:
        return hp - hq
    cross = ox * qy - oy * qx
    if cross != 0:
        return -1 if cross > 0 else 1
    return (ox * ox + oy * oy) - (qx * qx + qy * qy)


def argument_sort(points, origin=(0, 0)):
    """点列を origin まわりの偏角で安定に近い形でソートして返す。

    破壊しない。計算量 O(N log N)
    """
    ox, oy = origin
    from functools import cmp_to_key

    def cmp(p, q):
        a = (p[0] - ox, p[1] - oy)
        b = (q[0] - ox, q[1] - oy)
        return argument_cmp(a, b)

    return sorted(points, key=cmp_to_key(cmp))
