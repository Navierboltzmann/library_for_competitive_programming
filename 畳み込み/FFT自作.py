"""
浮動小数点 FFT による畳み込み (numpy 不使用・純 Python)

convolution(f, g)        : 整数列の畳み込み (結果を四捨五入した整数で返す)
convolution_mod(f, g, p) : 任意の法 p での畳み込み (各値を 2^15 で上下に分割して誤差を抑える)

複素数は実部・虚部を別々の float のリストで持ち (PyPy で高速)、
回転因子は長さごとに前計算してキャッシュする。

戻り値の長さは元実装との互換のため 1 << (len(f) + len(g) - 1).bit_length() で、
len(f) + len(g) - 1 より後ろは 0 で埋められている。

使用例:
    convolution([1, 2, 3], [4, 5])                # -> [4, 13, 22, 15]
    convolution_mod([10**9, 2], [10**9, 3], 10**9 + 7)[:3]
"""
import math

_fft_cache = {}


def _fft_table(n):
    """
    長さ n (2冪) の FFT 用テーブル (ビット反転表, 回転因子の実部, 虚部) を返す

    rt[k + j] = exp(i * pi * j / k) (k = 1, 2, 4, ..., n/2, 0 <= j < k)
    計算量: O(n) (初回のみ)
    """
    tab = _fft_cache.get(n)
    if tab is not None:
        return tab
    m = max(n, 2)
    rr = [1.0] * m
    ri = [0.0] * m
    pi = math.pi
    cos = math.cos
    sin = math.sin
    k = 2
    while k < n:
        for j in range(k):
            ang = pi * j / k
            rr[k + j] = cos(ang)
            ri[k + j] = sin(ang)
        k <<= 1
    h = n.bit_length() - 1
    rev = [0] * n
    for i in range(1, n):
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (h - 1))
    tab = (rev, rr, ri)
    _fft_cache[n] = tab
    return tab


def _fft(re, im):
    """
    re + i*im を長さ n (2冪) の FFT で変換する (in-place)

    re, im: 実部, 虚部の float のリスト (長さ n)
    計算量: O(n log n)
    """
    n = len(re)
    rev, rr, ri = _fft_table(n)
    for i in range(n):
        j = rev[i]
        if i < j:
            re[i], re[j] = re[j], re[i]
            im[i], im[j] = im[j], im[i]
    if n >= 2:
        for i in range(0, n, 2):
            ur = re[i]
            ui = im[i]
            vr = re[i + 1]
            vi = im[i + 1]
            re[i] = ur + vr
            im[i] = ui + vi
            re[i + 1] = ur - vr
            im[i + 1] = ui - vi
    k = 2
    while k < n:
        k2 = k << 1
        for s in range(0, n, k2):
            for i in range(s, s + k):
                t = i + k - s
                x = rr[t]
                y = ri[t]
                u = re[i + k]
                v = im[i + k]
                zr = x * u - y * v
                zi = x * v + y * u
                ar = re[i]
                ai = im[i]
                re[i + k] = ar - zr
                im[i + k] = ai - zi
                re[i] = ar + zr
                im[i] = ai + zi
        k = k2


def _out_len(f, g):
    """元実装と同じ戻り値の長さ 1 << (len(f) + len(g) - 1).bit_length() を返す"""
    return 1 << (len(bin(len(f) + len(g) - 1)) - 2)


def convolution(f, g):
    """
    整数列 f, g の畳み込みを浮動小数点 FFT で計算する

    f, g: 数値のリスト
    戻り値: h[k] = sum_{i+j=k} f[i] * g[j] を round した int のリスト
            (長さ 1 << (len(f)+len(g)-1).bit_length(), 末尾は 0 埋め)
    計算量: O(N log N) (N = len(f) + len(g))
    注意: 結果の絶対値が 10^15 程度を超えると誤差で正しくなくなる
    """
    out_len = _out_len(f, g)
    L = len(f) + len(g) - 1
    if not f or not g:
        return [0] * out_len
    n = 1 << (L - 1).bit_length()
    re = [0.0] * n
    im = [0.0] * n
    for i, x in enumerate(f):
        re[i] = float(x)
    for i, x in enumerate(g):
        im[i] = float(x)
    _fft(re, im)
    # (f + i g)^2 = f^2 - g^2 + 2i f g
    sr = [0.0] * n
    si = [0.0] * n
    for i in range(n):
        a = re[i]
        b = im[i]
        sr[i] = a * a - b * b
        si[i] = 2.0 * a * b
    # out[i] = in[-i] - conj(in[i]) を FFT すると虚部が 4n * (f*g)
    mask = n - 1
    outr = [0.0] * n
    outi = [0.0] * n
    for i in range(n):
        j = -i & mask
        outr[i] = sr[j] - sr[i]
        outi[i] = si[j] + si[i]
    _fft(outr, outi)
    inv = 1.0 / (4 * n)
    res = [round(outi[i] * inv) for i in range(L)]
    res += [0] * (out_len - L)
    return res


def convolution_mod(f, g, p):
    """
    f, g の畳み込みを法 p で計算する (p は任意, p < 2^31 程度)

    各値を 2^15 で上位・下位 (符号付き) に分割し、FFT 4 回で計算する。
    f, g: 整数のリスト
    p: 法
    戻り値: h[k] = (sum_{i+j=k} f[i] * g[j]) % p のリスト
            (長さ 1 << (len(f)+len(g)-1).bit_length(), 末尾は 0 埋め)
    計算量: O(N log N) (N = len(f) + len(g))
    注意: 浮動小数点誤差のため N ~ 10^6 程度が目安の上限 (N = 5*10^5, p ~ 2^31 まで検証済み)
    """
    out_len = _out_len(f, g)
    L = len(f) + len(g) - 1
    if not f or not g:
        return [0] * out_len
    n = 1 << (L - 1).bit_length()
    lr = [0.0] * n
    li = [0.0] * n
    rr = [0.0] * n
    ri = [0.0] * n
    # 誤差を減らすため各値を (-p/2, p/2] に寄せ、上位・下位とも符号付きで持つ
    half = p >> 1
    for i, x in enumerate(f):
        x %= p
        if x > half:
            x -= p
        lo = ((x + 16384) & 32767) - 16384
        lr[i] = float((x - lo) >> 15)
        li[i] = float(lo)
    for i, x in enumerate(g):
        x %= p
        if x > half:
            x -= p
        lo = ((x + 16384) & 32767) - 16384
        rr[i] = float((x - lo) >> 15)
        ri[i] = float(lo)
    _fft(lr, li)
    _fft(rr, ri)
    mask = n - 1
    inv = 0.5 / n
    ar = [0.0] * n
    ai = [0.0] * n
    br = [0.0] * n
    bi = [0.0] * n
    for i in range(n):
        j = -i & mask
        xr = lr[i]
        xi = li[i]
        yr = lr[j]
        yi = -li[j]
        gr = rr[i] * inv
        gi = ri[i] * inv
        # F_hi = (L[i] + conj(L[j])) / 2, F_lo = (L[i] - conj(L[j])) / 2i
        pr = xr + yr
        pi_ = xi + yi
        qr = xi - yi
        qi = yr - xr
        ar[j] = pr * gr - pi_ * gi
        ai[j] = pr * gi + pi_ * gr
        br[j] = qr * gr - qi * gi
        bi[j] = qr * gi + qi * gr
    _fft(ar, ai)
    _fft(br, bi)
    res = [0] * out_len
    for i in range(L):
        hh = round(ar[i]) % p
        mid = (round(ai[i]) + round(br[i])) % p
        ll = round(bi[i]) % p
        res[i] = ((hh << 30) + (mid << 15) + ll) % p
    return res


if __name__ == "__main__":
    print(convolution([1, 2, 3], [4, 5]))  # [4, 13, 22, 15]
    print(convolution_mod([10**9, 2], [10**9, 3], 10**9 + 7)[:3])
