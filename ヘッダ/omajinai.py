"""PyPy 向け再帰アンロール設定。

CPython では pypyjit が無いので何もしない。
再帰 DFS を PyPy で使うときにファイル先頭へコピペする。
"""
try:
    import pypyjit
    pypyjit.set_param('max_unroll_recursion=-1')
except ImportError:
    pass
