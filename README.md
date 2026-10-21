# Python テンプレート

提出用にコピペするライブラリ。各ファイルは独立している。

```
Python/
  ヘッダ/        提出の定番 import
  整数/          数論・合同式・篩
  畳み込み/      FFT / FPS / ゼータ
  行列/
  データ構造/
  グラフ/        LCA/ を含む
  文字列/
  幾何/
  DP/
  累積和/
  テクニック/    Mo・座標圧縮など
```

## ヘッダ

| ファイル | 内容 |
|---|---|
| へっだ.py | import・再帰上限・高速入力 |
| omajinai.py | PyPy の再帰アンロール |

## 整数

| ファイル | 内容 |
|---|---|
| ペル方程式.py | `x^2 - D y^2 = N`。±1 の基本解、同伴解、一般 N |
| 連分数.py | 有理数 / √D の連分数と収束分数 |
| 中国剰余定理.py | 拡張 Euclid、CRT（非素でも可）、Garner |
| 一次不定方程式.py | `ax + by = d` |
| 一次合同方程式.py | `Ax ≡ B (mod M)` |
| 二次合同方程式.py | `x^2 ≡ a (mod n)`（Hensel + CRT） |
| 有限体の平方根.py | Tonelli-Shanks（mod 素数） |
| ヤコビ記号.py | Legendre / Jacobi、平方剰余判定 |
| 原始根.py | 原始根・位数・φ・λ |
| 離散対数問題.py | Baby-step giant-step |
| 二平方和.py | `x^2 + y^2 = n`、Cornacchia |
| floor_sqrt.py | `⌊√n⌋` |
| 整数k乗根.py | `⌊n^{1/k}⌋`、完全冪判定 |
| floor_sum.py | `Σ ⌊(ai+b)/m⌋`（ACL） |
| CeilSumModSum.py | `Σ ⌊N/i⌋`、`Σ N mod i` |
| SternBrocot.py | Stern-Brocot / Farey、最良有理近似 |
| 素因数分解.py | 線形篩の最小素因数・μ・約数 |
| 素因数分解v2.py | Miller-Rabin + Pollard rho |
| getfactors.py | 試し割り |
| 線形篩.py | 素数・φ・μ・d(n)・σ(n) を O(N) |
| 約数列挙.py | `O(√n)` の約数 |
| オイラーのφ関数.py | 1..M の φ を篩 |
| 素数関数.py | `π(n)`（Lucy-Hedgehog） |
| 二項係数.py | 素数 mod、前計算 O(N) |
| 二項係数(非素数).py | 任意 mod（Lucas + Garner） |
| Lucasの定理.py | `C(m,n) mod p` |
| 階乗.py | 大きな `n! mod 998244353` |

## 畳み込み

| ファイル | 内容 |
|---|---|
| FFT.py | numpy FFT |
| FFT自作.py | 標準ライブラリのみ |
| FFT任意mod.py | NTT 任意 mod |
| FFT逆数.py | 形式的冪級数の逆元 |
| FPS.py | 形式的冪級数 |
| BM法.py | Berlekamp-Massey |
| 乗法畳み込み.py | 添字の積 |
| 添え字に大小関係がある畳み込み.py | 大小つき畳み込み |
| gcd_convolution.py | GCD 畳み込み・乗法ゼータ |
| XOR_AND_OR_convolution.py | 集合ゼータ (AND/OR/XOR) |
| xor_conv.py | Walsh-Hadamard のみ |

## 行列

| ファイル | 内容 |
|---|---|
| 行列クラス.py | 行列演算 |
| 行列累乗.py | 行列冪 |
| 抽象化行列累乗.py | 半環など |
| 行列式.py | 行列式 |

## データ構造

| ファイル | 内容 |
|---|---|
| DSU.py / union_find.py | Union-Find |
| DSU_子のリスト付き.py | 子リストつき |
| FenwickTree.py | BIT |
| SegmentTree.py | セグ木 |
| SegmentTree(非可換).py | 非可換 |
| lazysegtree.py | 遅延セグ木 |
| lazysegtree_区間加算.py / lazysegtree_区間更新.py | よく使う作用 |
| Beats.py | Segment Tree Beats |
| 2dSegmentTree.py | 2 次元 |
| SparseTable.py | 冪等な区間演算 |
| WaveletMatrix.py | Wavelet |
| mergesorttree.py | マージソート木 |
| BinaryTrie.py | 二進 Trie |
| ImplicitTreap.py | implicit treap |
| CHT.py | Convex Hull Trick |
| XorSegTree.py | XOR セグ木 |
| CountDistinct.py | 区間種類数 |
| sumoftopk.py | 上位 K 和 |
| 最大連続部分列.py | 区間最大部分列和 |

## グラフ

| ファイル | 内容 |
|---|---|
| 木.py | 根付き木（LCA・直径・全方位・オイラーツアー） |
| ダイクストラ.py / dijkstra.py | 最短路 |
| ポテンシャルダイクストラ.py | ポテンシャルつき |
| Bellman_Ford.py | 負辺 |
| SCC.py | 強連結成分 |
| トポロジカルソート.py | Kahn |
| Dinic.py | 最大流 |
| 最小費用流.py | 最小費用流 |
| 二部マッチング.py | 二部最大マッチング |
| 最小全域木.py | Kruskal |
| 橋と関節点.py | lowlink |
| 区間に辺を張るテク.py | 区間辺 |
| 巡回セールスマン.py | TSP |
| ReRooting.py | 全方位木 DP |
| HL分解.py / HL分解v2.py | 重軽分解 |
| AuxiliaryTree.py | 補助木 |
| なもり分解.py | なもり分解 |
| 部分木サイズ.py | 部分木サイズ |
| neighbor.py | 隣接 |
| LCA/ | 各種 LCA・部分木クエリ |

## 文字列

| ファイル | 内容 |
|---|---|
| ローリングハッシュ.py | 複数 mod |
| 最長回文.py | Manacher |
| Zアルゴリズム.py | Z-algorithm |
| KMP.py | KMP / MP |
| Trie.py | 配列 Trie |

## 幾何

| ファイル | 内容 |
|---|---|
| convex_hull.py | 凸包 |
| 偏角ソート.py | 偏角比較（atan2 なし） |
| 線分交差.py | 交差判定・交点 |
| 平面の方程式.py | 平面 |

## DP / 累積和 / テクニック

| ファイル | 内容 |
|---|---|
| LIS.py | 最長増加部分列 |
| 部分和DP.py | 部分和 |
| 転倒数.py | 転倒数 |
| subsegment_max.py | 部分列 max |
| hopping_imos.py | 一定間隔の累積 / imos |
| triangle_imos.py | 三角 imos |
| imos直線.py | 直線 imos |
| 二次元累積和.py | 2 次元累積 |
| Mo.py | Mo's algorithm |
| 平方分割.py | 平方分割 |
| 座標圧縮.py | 座標圧縮 |
| slope_trick.py | slope trick |
| min_upper_*.py / count_lower_*.py | 二分探索ヘルパ |
