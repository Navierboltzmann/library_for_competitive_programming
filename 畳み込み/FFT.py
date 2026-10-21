"""numpy FFT による畳み込み。

numpy が必要。大きな係数は 15 bit 分割して mod 畳み込みする。

使用例:
    convolve([1, 2, 3], [4, 5])
    convolve2([1, 2], [3, 4], 998244353)
    convolve_pow([1, 1], 5, 998244353)
"""
try:
    import numpy as np
except ImportError:  # numpy が無い環境でもテンプレート自体は import できるようにする。
    np = None


def _require_numpy():
    """numpy が利用できることを確認し、無ければ明確な例外を出す。"""
    if np is None:
        raise ImportError("FFT.py requires numpy; use FFT自作.py in a standard-library-only environment")


def convolve(f, g):
    """整数列 f, g の線形畳み込み(丸め誤差あり)。

    計算量: O((N+M) log (N+M))
    """
    _require_numpy()
    if len(f) == 0 or len(g) == 0:
        return np.array([], dtype=np.int64)
    tf = np.array(f, np.int64)
    tg = np.array(g, np.int64)
    fft_len = 1
    while 2 * fft_len < len(tf) + len(tg) - 1:
        fft_len *= 2
    fft_len *= 2
    Ff = np.fft.rfft(tf, fft_len)
    Fg = np.fft.rfft(tg, fft_len)
    h = np.fft.irfft(Ff * Fg, fft_len)
    h = np.rint(h).astype(np.int64)
    return h[:len(f) + len(g) - 1]


def convolve2(f, g, p):
    """15 bit 分割で f * g mod p を計算する。

    計算量: O((N+M) log (N+M))
    """
    _require_numpy()
    if len(f) == 0 or len(g) == 0:
        return np.array([], dtype=np.int64)
    f1, f2 = np.divmod(f, 1 << 15)
    g1, g2 = np.divmod(g, 1 << 15)
    a = convolve(f1, g1) % p
    c = convolve(f2, g2) % p
    b = (convolve(f1 + f2, g1 + g2) - a - c) % p
    h = (a << 30) + (b << 15) + c
    return h % p


def convolve_pow(f, n, p):
    """多項式 f の n 乗を mod p で返す(次数は len(f) で打ち切り)。

    計算量: O(N log N log n)
    """
    _require_numpy()
    N = len(f)
    if N == 0:
        return np.array([], dtype=np.int64)
    C = np.zeros(N, dtype=np.int64)
    C[0] = 1
    B = np.asarray(f, dtype=np.int64)
    while n:
        if n & 1:
            C = convolve2(C, B, p)[:N]
        B = convolve2(B, B, p)[:N]
        n >>= 1
    return C
