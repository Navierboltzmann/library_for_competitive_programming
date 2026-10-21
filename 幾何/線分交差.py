"""線分の交差判定と交点。

端点共有・平行・重なりを区別する。

使用例:
    segments_intersect((0, 0), (2, 2), (0, 2), (2, 0))  # True
    segment_intersection((0, 0), (2, 2), (0, 2), (2, 0))  # (1.0, 1.0)
"""


def _cross(ax, ay, bx, by):
    """2次元外積。"""
    return ax * by - ay * bx


def _dot(ax, ay, bx, by):
    """内積。"""
    return ax * bx + ay * by


def ccw(a, b, c):
    """a->b に対する c の位置。

    1: 左, -1: 右, 2: b の先, -2: a の手前, 0: 辺上。
    """
    ax, ay = a
    bx, by = b[0] - ax, b[1] - ay
    cx, cy = c[0] - ax, c[1] - ay
    cr = _cross(bx, by, cx, cy)
    if cr > 0:
        return 1
    if cr < 0:
        return -1
    if _dot(bx, by, cx, cy) < 0:
        return -2
    if bx * bx + by * by < cx * cx + cy * cy:
        return 2
    return 0


def segments_intersect(p1, p2, p3, p4, endpoint=True):
    """線分 p1p2 と p3p4 が交差するか。

    endpoint=False なら端点だけの接触は False。
    """
    a, b, c, d = ccw(p1, p2, p3), ccw(p1, p2, p4), ccw(p3, p4, p1), ccw(p3, p4, p2)
    if a * b * c * d == 0 and min(a, b, c, d) == 0:
        # 端点が他方の線分上
        return endpoint
    return a * b <= 0 and c * d <= 0


def segment_intersection(p1, p2, p3, p4):
    """交点を 1 つ返す。平行で重なるときは端点の 1 つ。交わらなければ None。

    戻り値は float のタプル。
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4
    den = _cross(x2 - x1, y2 - y1, x4 - x3, y4 - y3)
    if den == 0:
        if not segments_intersect(p1, p2, p3, p4):
            return None
        # 重なり: 共有点を探す
        for p in (p1, p2, p3, p4):
            if ccw(p1, p2, p) == 0 and ccw(p3, p4, p) == 0:
                return (float(p[0]), float(p[1]))
        return None
    t = _cross(x3 - x1, y3 - y1, x4 - x3, y4 - y3) / den
    u = _cross(x3 - x1, y3 - y1, x2 - x1, y2 - y1) / den
    if 0 <= t <= 1 and 0 <= u <= 1:
        return (x1 + t * (x2 - x1), y1 + t * (y2 - y1))
    return None
