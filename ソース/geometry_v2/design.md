# 幾何テンプレート 第二世代 設計書

- 作成日：2026-10-07
- 状態：**Phase 1（体系設計）完了**（2026-10-08 ユーザー確認済み。6-1）。**Phase 2（基盤 00_geo_base）完了**（2026-10-08。2-P2 章）。**Phase 3（★ 7 本）完了**（2026-10-08。2-P3 章）。**Phase 4（☆ 8 本）完了**（2026-10-08。2-P4 章）。**最終照合完了**（2026-10-09〜10。7 章：図形の位置関係の網羅・旧テンプレート・過去の解答・Dynalist・典型 48 問）
- 方法：グラフ・グリッド第二世代（`Desktop\claude\graph_grid_v2\design.md`、引き継ぎメモ `Desktop\claude\graph_grid_template_review.md`）と同じ方法・観点で作り直す
- 次回はこのファイルの「6. 未決事項と次にやること」と「7. 最終照合」から再開する
- **2026-10-10 追記**：数学第二世代（`Desktop\claude\math_v2`）と共通の整数の部分を `00_int_common` に切り出し、16 本すべての同梱部分を更新した（7-9）

---

## 0. 目的・観点・グラフ・グリッドから引き継ぐこと

### 0-1. 目的

幾何のテンプレートを「第二世代」としてゼロから設計し直す。既存テンプレートの整理ではなく、新しい体系を作る。

成果物はテンプレートそのものと、**問題の特徴から使う道具に辿れる選択体系**（Dynalist の索引）。例えば次のように辿れる状態を目指す。

- 点の集合 → 外側を囲む → 凸包（04_convex_hull）
- 点が多角形の内側か → 03_polygon（凸でクエリが多いなら 04_convex_hull の O(log N) 版）
- 円が交わるか・つながるか → 05_circle（整数で厳密に）
- 座標が小数で与えられる → 00_geo_base の parse_fixed で 10^k 倍の整数にする
- 角度の順に並べる・ある向きの範囲の点を数える → 01_arg_sort

### 0-2. 観点（グラフ・グリッドと同じ 3 つ）

| 観点 | 幾何での具体化 |
|---|---|
| アクセスしやすさ | 「入力の形（点・直線・多角形・円・点の集合・マンハッタン距離・長方形）→ 聞かれていること → ファイルと関数」で辿れる索引。ファイル名は分類番号＋名前。依存するファイルは同梱し、1 回のコピーで動く |
| 使いやすさ | 型は整数の `Point` と実数の `PointD` の 2 つだけ。判定は整数で厳密（誤差を考えなくてよい）。呼び方・戻り値の形をそろえる（交点は `vector<PointD>`、位置関係は名前つきの定数、点の組は番号を返す）。入力関数・小数の入力を整数に直す関数つき |
| 連携 | 点・直線・線分・円の型を基盤で統一し、すべてのファイルが同じ型を受け取る。グラフ・グリッド第二世代（46 ファイル）と同時に貼れる（名前が衝突しない）。グラフに渡す形（距離の関数・辺のリスト）をそろえ、組み合わせ方を索引に載せる |

### 0-3. 最重要の方針（グラフ・グリッドの 1-0 を最初から適用する）

第二世代のテンプレートは互いに関係し、旧テンプレートとの互換性もないため、**後から追加・修正しなくて済むこと**が最も重要。

- 追加でほしくなる可能性のある機能はすべて入れる。入れるか迷った機能は入れる
- 速さの改善が要らないくらい速くする
- 正しいこと（愚直解とのランダムテスト、サニタイザ、Main_base＋ACL に貼った状態でのコンパイル確認）
- 使いやすいこと

### 0-4. グラフ・グリッドでの経験から、最初から取り入れること

| グラフ・グリッドで起きたこと | 幾何での対応 |
|---|---|
| 1-0 の方針が Phase 2 の途中で決まった | 最初から適用する（0-3） |
| 最終照合（AtCoder の典型との照合）で、足りないファイルが 7 本見つかった | **典型との照合を Phase 1 で行う**（5-5）。最終照合でも問題文を確かめながらもう一度行う |
| △（型）にしたワーシャルフロイドを、間違えやすいのでテンプレートにした | 「数行で書けるが間違えやすい」ものは最初からテンプレートに入れる（偏角の比較、角の大きさの比較、小数の入力の整数化など） |
| std::priority_queue・lower_bound がデバッグビルドで O(n) の検査をして終わらなくなった | 本体の二分探索・ヒープは自分で書く（凸多角形の O(log N) 判定、偏角の lower_bound など） |
| MSVC が手元になく、GCC 拡張を避けた | **MSVC ではビルドしない・`__int128` は使える**（2026-10-08 ユーザー確認）。`__int128` をそのまま使う。ヘッダは今までどおり使うものを書く |
| 名前の衝突を Main_base＋ACL で確認した | **グラフ・グリッド第二世代の 46 ファイルも一緒に貼って**確認する |
| 各 Phase で Dynalist の項目と索引を作った | 同じ（dynalist_items_phaseN.txt、dynalist_tool_index.txt / .opml） |

---

## 1. 決定事項

### 1-1. 進め方・既存テンプレートの扱い

- 既存テンプレートは「過去の実験結果」として扱う。既存の API を正解として引き継がない
- 参照したのは、新しい一覧と機能を決めた後（Phase 1 の終盤、5 章）だけ
- 原本は変更・削除しない。新しく作るものは `Desktop\claude\geometry_v2\` に出力する（最後に `C:\atcoder_template\典型90\ソース\geometry_v2\` へ置く想定。グラフ・グリッドと同じ）
- 既存の場所：`C:\atcoder_template\典型90\ソース`（GeometryTemplateTest.txt など。5-0）
- Phase 1（体系設計）→ Phase 2（基盤）→ Phase 3（★）→ Phase 4（☆）→ 最終照合 の順に進め、区切りごとにユーザーの確認を取る
- ACL は使わない（幾何では ACL に使えるものもない）。std だけで自作する

### 1-2. 共通規約

グラフ・グリッドと同じもの：

| 項目 | 決定 |
|---|---|
| 優先順位 | 速度 → int / long long と変わらない使い勝手 → 正しさ |
| ファイル | 1テンプレート1ファイル。単体で実行できる .cpp |
| テスト | ランダムテスト関数と、全機能を使う使い方関数を同じファイルに入れ、main から呼ぶ（形は 9/23 版と同じ：`xxx_test::random_test()`（Release は 8 倍の量）・`benchmark()`・`xxx_demo::demo()`・`setup_console()`） |
| 使い方の出力 | 関数名の表示に加え、何を返す・提供するかを日本語でも書く |
| 使うライブラリ | std のみ（ACL 不使用）。bits/stdc++.h に頼らず、使うヘッダを書く |
| 番号 | 点は入力順に 0-indexed。点の組・最も遠い点などは**座標ではなく番号で返す**（座標は `ps[i]` で取れる）。「なし」は -1 |
| 名前 | 型は PascalCase（`Point`, `Circle`）、関数・メンバは snake_case（`convex_hull`, `cross_points`） |
| Main_base への依存 | テンプレートの中では Main_base の別名・マクロを使わない |
| 再帰 | 非再帰を基本にする。深さが log N で済むもの（最近点対の分割統治）だけ再帰でよい |

幾何で新しく決めたもの（★ は特に確認してほしいもの。6-1）：

| 項目 | 決定 |
|---|---|
| ★ 座標の型 | **整数の `Point`（`Coord` = long long）と実数の `PointD`（`Real` = long double）の 2 つ**。どちらも `PointT<T>` の別名。`using` の 1 行で型を変えられる |
| ★ 整数で厳密に | 判定（向き・交差・内外・円の位置関係・偏角の比較・角の大きさの比較）と、面積の 2 倍・距離の 2 乗は**整数で厳密に**求める。誤差（EPS）を使わない |
| ★ 座標の範囲 | 整数版はすべての関数で **\|x\|, \|y\| ≤ 10^9（半径も 10^9）まで正しい**ことを保証する。途中の計算で long long をはみ出すところは 128 ビットで計算する。答えそのものが long long に収まらないことがあるもの（全部の組の距離の和など）は、Int128 で返す版か mod を取る版を用意し、説明に書く |
| ★ 128 ビット | **`__int128` をそのまま使う**（型名 `Int128` = `__int128`。MSVC ではビルドしないため自作の代わりは作らない。2026-10-08 確認）。利用者が Int128 を直接扱う場面はなるべく作らない（比較の結果だけを返す）。出力用に `int128_to_string` を用意する |
| 実数 | 答えが整数にならないもの（交点・射影・接点・円の中心・最小包含円など）だけ `PointD` / `Real` で返す。実数どうしの比較は**相対誤差** `GEO_EPS`（= 1e-9、1 行で変えられる）。Real を double に変えても（速さのために変えることがあるので）テストが通るようにする |
| 角度 | 偏角は x 軸の正の向きから反時計回りに [0, 2π)。0 ベクトル（原点そのもの）は先頭。実数の角はラジアン（度で受け取る・返す関数は名前に `_deg` をつける） |
| 向き | 反時計回りが正（`cross(a, b) > 0` は b が a の左）。直線・線分は向き（a → b）を持つ。左側 = 正 |
| 座標の向き | 幾何の (x, y) は「横・縦（上が正）」。グリッドの (i, j)（行・列）とは回転の向きが逆になる。グリッドのマスを幾何で扱うときの注意を基盤の説明に書く |
| 多角形 | `vector<Point>`。始点を末尾に繰り返さない。**返す多角形は反時計回り**。受け取る多角形は向きを問わない（向きが要る関数は説明に書き、`make_ccw` を用意する） |
| 戻り値 | 交点・接点 → `vector<PointD>`（なしは空）。位置関係 → 名前つきの定数（Windows のマクロ `IN` / `OUT` とぶつからない名前。例：`OUTSIDE` / `ON_BORDER` / `INSIDE`）。点の組 → 番号と値を持つ構造体（pair・tuple で返さない） |
| 長方形 | 半開 [x1, x2) × [y1, y2)（グラフ・グリッドの区間と同じ。面積 = (x2 − x1)(y2 − y1)）。格子点を閉区間で数える関数は名前で区別する |
| 名前の衝突 | Main_base（`sgn`, `INF`, `LINF`, `FIXED`, `PQ`, `chmin` など）・ACL・**グラフ・グリッド第二世代 46 ファイル**・std の関数名（`abs`, `arg`, `norm`, `proj`, `polar`, `rotate`, `distance`, `real`, `imag`, `conj`）・よく使う変数名（N, M, P, Q, R, A, B, C, D, X, Y, PI, EPS）・Windows のマクロ（`IN`, `OUT`, `min`, `max`）と同じ名前をグローバルに置かない。定数は `GEO_EPS`・`GEO_PI` |
| 出力 | テンプレートは cout しない。実数の出力のしかた（`fixed << setprecision(15)`）は使い方関数で示す |

### 1-3. 共通API

#### 型：基盤で 1 回だけ定義し、すべてのファイルが同じ型を受け取る

```cpp
using Coord = long long;          // 整数の座標
using Real  = long double;        // 実数
constexpr Real GEO_EPS = 1e-9;    // 実数の比較（相対誤差）
using Int128 = __int128;         // 途中の計算用（出力は int128_to_string）

template<class T> struct PointT { T x, y; };   // + - * / == < など
using Point  = PointT<Coord>;     // 整数の点（基本はこちら）
using PointD = PointT<Real>;      // 実数の点（交点などの答え）

template<class T> struct LineT    { PointT<T> a, b; };  // 直線 a → b（a ≠ b）
template<class T> struct SegmentT { PointT<T> a, b; };  // 線分 a - b
template<class T> struct CircleT  { PointT<T> c; T r; };// 円（中心と半径）
using Line = LineT<Coord>;        using LineD = LineT<Real>;
using Segment = SegmentT<Coord>;  using SegmentD = SegmentT<Real>;
using Circle = CircleT<Coord>;    using CircleD = CircleT<Real>;
```

- 直線と線分は同じ形だが**別の型**にする。`intersects(s1, s2)`（線分どうし）と `intersects(l, s)`（直線と線分）のように、同じ名前の関数を型で使い分けられる
- 整数の型を渡すと厳密な判定、実数の型を渡すと誤差つきの判定になる（同じ関数名）。整数を受け取って実数を返す関数（交点など）は `PointD` を返す。`to_d(p)` で実数に直せる

#### 基本の呼び方

```cpp
Point a = read_point();                 // 読み込み（read_points(N) で N 個）
vector<Point> ps = read_points(N);
Coord k = parse_fixed("1.2345", 4);    // 小数を 10^4 倍した整数 12345 に（誤差なし）

dot(a, b); cross(a, b); norm2(a);       // 整数のまま
dist2(a, b);                            // 距離の 2 乗（整数）
dist(a, b);                             // 距離（Real）
ccw(a, b, c);                           // +1 左に曲がる / -1 右 / 0 一直線（厳密）
ccw_detail(a, b, c);                    // 一直線のときの位置まで（5 分類。名前つき定数）
rot90(a);                               // 反時計回りに 90 度（整数のまま）

auto ch = convex_hull(ps);              // 反時計回りの凸包
auto where = contains(poly, q);         // OUTSIDE / ON_BORDER / INSIDE
auto cp = cross_points(c1, c2);         // 2 円の交点（vector<PointD>。0〜2 個）
auto [i, j, d2] = closest_pair(ps);     // 最も近い 2 点の番号と距離の 2 乗
```

- 1 回で終わる計算は関数（`convex_hull(ps)`、`closest_pair(ps)`）、前処理してから何度も聞くものは構造体（`ArgSorted S(ps, o)`、`ManhattanFarthest F(ps)`）
- 複数の値を返す関数は、名前付きのメンバを持つ構造体で返す（構造化束縛でも受け取れる）

### 1-4. ファイルの形と共通化

- グラフ・グリッドと同じ：**.cpp 1 ファイル**、コピー範囲の目印（`ここからテンプレート本体（コピー範囲）` / `ここまでテンプレート本体`）、依存するテンプレートは**貼る範囲に同梱**して `#ifndef` ガードで重複を防ぐ
- ガード名は**グラフ・グリッドと同じ付け方**（`GG_` ＋ファイル名の大文字。例：`GG_GEO_BASE`、`GG_CONVEX_HULL`）。こうすると `update_bundles.ps1 / .bat` を**そのまま**使える（geometry_v2 にコピーするだけ）。`GG_` は「第二世代テンプレート共通の接頭辞」と読み替える
- ファイル名は「分類番号 2 桁＋名前」（例：`04_convex_hull.cpp`）。**名前の部分はグラフ・グリッドのファイルと重ならない**ようにする（ガードが重ならないため）
- 内部の名前空間は `geo_internal`（グラフ・グリッドの `gg_internal` とは分ける）
- 文字コードは UTF-8（BOM 付き）、main で Windows のコンソールの設定をする（グラフ・グリッドと同じ）

### 1-5. 正しさ・速さの確認のしかた（幾何で追加するもの）

- **小さい座標で退化ケースを大量に作る**：\|x\| ≤ 3〜10 の点で、一直線・重なり・端点で接する・同じ点・面積 0 の多角形などを意図して作り、別の方法（全探索・有理数での計算）と比べる
- **大きい座標**：\|x\| = 10^9 付近の点と、そこで退化するケース（a + k·d の形で一直線、接する円など）を作り、128 ビットではみ出さないことを確かめる
- **実数の答え**：答えが条件を満たすか（交点が両方の図形の上にある、接点で半径と接線が垂直など）と、別の式で求めた値との比較を相対誤差で
- **Real を double にした版**でも全テストを通す（Real は 1 行で変えられるため）
- **名前の衝突**：Main_base＋ACL に、グラフ・グリッド第二世代 46 ファイルと幾何の全ファイルの本体を 2 通りの順で全部貼り、gnu++2b / c++17・LOCAL あり / なしでコンパイル・実行
- **速さ**：N = 2×10^5 で凸包・偏角ソート・最近点対・マンハッタン・長方形の和集合、N = 2000 で O(N² log N) の数え上げ（一直線の組・直線の数）

---

## 2. 分類とテンプレート一覧（確定。2026-10-08）

記号：★＝優先して整備（Phase 3）、☆＝余裕が出てから（Phase 4。ただしグラフ・グリッドと同じく全部作る）、△＝専用テンプレートは作らず「型」として持つ

関数名は案。Phase 2〜4 の実装で確定し、この表に反映する（グラフ・グリッドと同じ）。

### 0. 基盤（Phase 2）

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `00_geo_base` | **作成済み（Phase 2）**。`Coord`・`Real`・`Int128`（= __int128）・`GEO_EPS`・`GEO_PI`。`PointT`（`Point` / `PointD`）：+ − 単項− ×k k× ÷k += −= *= /= == != < > <= >=（x → y の辞書順）。整数の点を実数から作る・実数を掛ける／割るとコンパイルエラー、Point → PointD は `to_d` か明示の変換だけ。`dot` / `cross`（2 点・3 点の版）/ `norm2` / `dist2` / `len` / `dist`、`dot128` / `cross128` / `norm2_128`、`ccw`（3 分類）/ `ccw_detail`（5 分類。`CcwDetail` の名前つき定数）/ `collinear`、`rot90` / `rot90_cw` / `rot180` / `rot90(p, 中心, k)`、`rotated(p, θ [, 中心])` / `rotated_deg(p, 度 [, 中心])`（90 度の倍数は厳密）、`complex_mul` / `conjugate`、`to_d`（点・直線・線分・円）、`midpoint` / `unit` / `from_polar`、`angle_of` / `angle_of_deg` / `angle_between` / `angle_between_deg` / `deg_to_rad` / `rad_to_deg`、`real_sign` / `real_cmp` / `real_eq`（点どうしも）、`direction_key` / `line_direction_key`、`isqrt` / `isqrt_ceil`、`parse_fixed` / `read_fixed` / `read_point_fixed` / `read_points_fixed`、`read_point` / `read_points` / `read_point_d` / `read_points_d` / `read_segment` / `read_line` / `read_circle`、`int128_to_string`、`PointHash`、`LineT` / `SegmentT` / `CircleT` と別名（`Line` / `LineD` / `Segment` / `SegmentD` / `Circle` / `CircleD`）、`Where`（OUTSIDE / ON_BORDER / INSIDE）、`operator<<`（`(x, y)` の形）/ `operator>>`。**最終照合で追加（7 章）：`rotated` / `rotated_deg` の中心の型を自由に（Point を PointD のまわりに回せる）・結果の −0 を 0 に、`AffineT`（`Affine` / `AffineD`。`then` で合成、rot90 / rot90_cw / rot180 / reflect_x / reflect_y / translate / scale / rotate）** | なし | ABC259 B（回転）、ABC168 C（時計の針）、ABC197 D（中心のまわりの回転）、ABC198 C（isqrt）、ABC226 D（direction_key）、ABC189 E（Affine） |

### 1. 角度・偏角

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `01_arg_sort` | **作成済み（Phase 3）**。`arg_half` / `arg_cmp`（−1 / 0 / 1）/ `arg_less` / `same_direction`、`arg_sort(ps, 原点, far_first)` / `arg_sort_idx`（原点と同じ点が先頭、同じ向きは近い順・遠い順）、`arg_sort_from` / `arg_sort_idx_from`（向きを指定して一周）、`unique_directions` / `unique_directions_idx`（向きごとに近い点・遠い点）。`ArgSorted S(ps, 原点)`：order / vec / zeros / size / lower_bound / upper_bound / count(l, r, include_l, include_r)（**0 度をまたいでよい**）/ exists / groups。`angle_cmp` / `angle_less`（∠(a, b) と ∠(c, d) を 128 ビットで厳密に）。`count_in_open_half(ps, 原点)` | geo_base | ABC225 E（角度の区間）、ABC139 F（ベクトルの和の最大）、典型90 009（最大の角） |

### 2. 直線・線分

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `02_line` | **作成済み（Phase 3）**。判定（整数なら厳密）：`is_parallel` / `is_orthogonal`（直線・線分のどの組でも）/ `on_line` / `on_segment` / `side` / `intersects(s1, s2)`（端点・重なりも含む）/ `intersects_strict` / `intersects(l, s)`・`(s, l)`・`(l1, l2)` / `line_relation`（`LINE_CROSS` / `LINE_PARALLEL` / `LINE_SAME`）。実数：`projection` / `reflection` / `dist`（点・直線・線分のどの組でも、順は逆でもよい）/ `cross_point(l1, l2)` / `cross_point(l, s)`・`(s, l)` / `cross_points(s1, s2)`（重なるときは両端）。厳密：`dist2_cmp(l or s, p, D2)`（距離の 2 乗と整数の比較。**設計の「距離の 2 乗を分数で」から変更**：分数どうしの比較は 256 ビットが要るため、整数との比較にした）/ `cross_point_frac`（`FracPoint{x, y, d}`・`FracPointHash`）。作る：`perpendicular` / `parallel_through`（整数のまま）/ `bisector`。`LineKey`（contains / dir_key / line / `LineKeyHash`）/ `line_key(p, q)`・`(l)` / `bisector_key`。`lattice_points_on` / `cells_crossed`。**最終照合で追加：`segment_relation`（`SEGMENT_SEPARATE` / `TOUCH` / `CROSS` / `OVERLAP`）、半直線 `RayT`（`Ray` / `RayD`）：on_ray / intersects（直線・線分・半直線）/ cross_point / cross_points / dist / dist2_cmp** | geo_base | ABC226 D（向きの正規化）、ABC248 E（直線の正規化）、ARC042 B（点と線分の距離）、ABC016 D（intersects_strict） |

### 3. 多角形・三角形

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `03_polygon` | **作成済み（Phase 3）**。`area2`（符号付き面積の 2 倍。整数で厳密）/ `area` / `is_ccw` / `make_ccw`、`contains(poly, p)`（凸でなくてよい）、`is_convex(poly, strict)`（向きはどちらでもよい。一周だけ回るかも確かめ、星形五角形を除く）、`is_simple`（O(N²)）、`perimeter`、`centroid`、`boundary_lattice_points` / `interior_lattice_points`、三角形：`triangle_area2` / `triangle_kind`（`TRIANGLE_ACUTE` / `TRIANGLE_RIGHT` / `TRIANGLE_OBTUSE` / `TRIANGLE_DEGENERATE`。`RIGHT` などとぶつからない名前）/ `in_triangle`。**最終照合で追加：位置関係（周を含む閉じた領域として）`contains(poly, s)`（線分が入るか。128 ビットの半整数点で厳密）/ `contains(P, Q)` / `intersects(poly, s or l or r)` / `intersects(P, Q)` / `dist(poly, p or s)` / `dist(P, Q)` / `boundary_dist`、`PolygonPrefix`（頂点 i..j の部分の面積の 2 倍・重心を O(1)）** | geo_base, **line（最終照合で追加）** | ABC266 C（凸四角形）、ABC362 B（直角三角形）、典型90 041（ピックの定理）、ABC250 F（PolygonPrefix） |

### 4. 凸包・凸多角形

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `04_convex_hull` | **作成済み（Phase 3）**。`convex_hull(ps, keep_collinear)` / `convex_hull_idx`（反時計回り。(x, y) 最小の点から。同じ点・1〜2 点・全部一直線にも対応）/ `lower_hull` / `upper_hull`、`farthest_pair(ps)`（`FarthestPair{i, j, d2}`）/ `convex_diameter(h)`、`convex_contains(h, p)`（O(log N)。頂点 0〜2 個の退化した形・一直線の頂点があってもよい）、**`lattice_points_in_convex(h)`（追加。凸包の中の格子点の数。退化した形も正しく数える）**。**最終照合で追加：`convex_contains(h, Q)`（多角形が入るか）** | geo_base, **line**, polygon | 典型90 041、ABC296 G（凸多角形の内外判定） |
| ☆ `04_convex_polygon` | **作成済み（Phase 4）**。h は convex_hull の結果（頂点 0〜2 個・一直線の形も可）。`extreme_point(h, dir)`（dot(dir, h[i]) が最大の頂点。O(log N)）/ `tangents(h, p)`（`ConvexTangents{right, left}`：p から見て右・左の接点。一直線に並ぶ頂点は p に近い方。内側・周の上なら -1。O(log N)）/ **`convex_intersects_line(h, l)`（追加。O(log N)）** / `convex_cut(poly, l)`（直線の左側。PointD。点・線分・一直線の形は両端の線分として）/ `convex_intersection(P, Q)`（09_halfplane で。空・線分・1 点も）/ `minkowski_sum(P, Q)`（整数で厳密。convex_hull と同じ形）/ `convex_width(h)`。**最終照合で追加：`convex_intersects(P, Q)`**。使い方関数に凸多角形の尺取り（ABC250 F の型。全探索と 28 万例で一致を確認） | geo_base, line, polygon, convex_hull, **halfplane** | ABC250 F（凸多角形を切る尺取り） |

### 5. 円

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `05_circle` | **作成済み（Phase 3）**。厳密：`contains(c, p)` / `circle_relation`（`CIRCLE_SEPARATE`(4) / `CIRCLE_EXTERNAL`(3) / `CIRCLE_INTERSECT`(2) / `CIRCLE_INTERNAL`(1) / `CIRCLE_CONTAIN`(0) / `CIRCLE_SAME`(5)）/ `has_common_point` / `count_cross(c, l)`・`(c, s)`。実数：`cross_points(c, l)`・`(c, s)`・`(c1, c2)` / `tangent_points` / `common_tangents`（`LineD`。a は c1 との接点）/ `circumcircle` / `incircle`（一直線なら r = −1）/ `circle_with_diameter` / `circles_through(a, b, r)` / `arc_length` / `sector_area`。`lattice_points_in_circle(cx, cy, r, scale, include_border)`。**最終照合で追加：`in_circumcircle`（4 点が同じ円の上か。192 ビットで厳密）、円と半直線 `count_cross` / `cross_points`、円板としての `disk_intersects` / `disk_contains` / `disk_dist` と `circumference_dist`、円周の上の弦 `chords_cross` / `chord_direction_key` / `chord_cross_counts` / `count_crossing_chords` / `ChordSet`** | geo_base, line | ABC259 D（円がつながるか）、ABC191 D（円の中の格子点）、ABC157 F（2 円の交点の候補）、ABC338 E・ABC405 F・ABC424 F・典型90 017・ABC402 D（弦） |
| ☆ `05_circle_cover` | **作成済み（Phase 4）**。`min_enclosing_circle(ps)`（最小包含円。Welzl・期待 O(N)。点がなければ r = -1。円の内側の判定の誤差は「半径 × GEO_EPS」と「座標の丸め」の大きい方）、`max_points_in_circle(ps, r)`（`CircleCover{count, center}`。円周の上も数える。O(N² log N)）。点は Point / PointD のどちらでも | geo_base, line, circle | ABC151 F（最小包含円） |
| ☆ `05_circle_area` | **作成済み（Phase 4）**。`intersection_area(c1, c2)`（2 円の共通部分の面積。整数なら場合分けは厳密）、`intersection_area(c, poly)`（円と多角形（凸でなくてよい・向きは自由）の共通部分の面積。O(N)）、`circular_segment_area(r, θ)`（弓形。角が小さいときは級数）、**`circle_cut_area(c, l)`（追加。円の、直線の左側の面積）**。**最終照合で追加：円と多角形の位置関係 `contains(poly, c)` / `disk_contains(c, poly)` / `disk_intersects(c, poly)` / `disk_dist(c, poly)`** | geo_base, line, polygon, circle | ARC051 A（円と長方形） |

### 6. 点の集合（数え上げ・最近点）

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `06_point_set` | **作成済み（Phase 3）**。`closest_pair`（`ClosestPair{i, j, d2}`）、`pairs_within(ps, D)` / `pairs_within_d2(ps, D2)`、`count_collinear_triples` / `count_triangles` / `max_collinear`（`MaxCollinear{count, i, j}`）/ `lines_through`（`LineCount{key, count}`）/ `count_lines_with_at_least`（無限なら −1）、`count_triangles_by_angle`（`TriangleCounts{acute, right, obtuse}`）/ `max_angle`（`MaxAngle{angle, i, j, k}`）、`count_axis_rectangles` / `count_squares` / `count_parallelograms` / `count_trapezoids` | geo_base, line, arg_sort | ABC181 C、ABC224 C、ABC248 E、ABC033 D、典型90 009、ABC218 D、ABC275 C、ABC304 C |
| ☆ `06_congruence` | **作成済み（Phase 4）**。`congruent(A, B, allow_flip = false)` / `congruence(A, B, allow_flip)`（`Congruence{ok, flipped, to, angle, shift}`：B[to[i]] = R(angle)·F(A[i]) + shift。回転は 90 度の倍数でなくてよい。整数で厳密。O(N log N)）、**`polygon_congruent(P, Q, allow_flip)` / `polygon_congruence`（追加。頂点の順番を保つ。始点・回り方は自由。O(N)）** | geo_base | ABC207 D |

### 7. マンハッタン距離・チェビシェフ距離

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `07_manhattan` | **作成済み（Phase 3）**。`manhattan` / `chebyshev` / `in_diamond`、`manhattan_to_chebyshev(p)` / `chebyshev_to_manhattan(q)`（10_grid_transform の (x, y) 版と同じ名前の Point 版）、`ManhattanFarthest F(ps)`：max_pair（`ManhattanPair{i, j, d}`）/ farthest(q)（`ManhattanFar{i, d}`）、`DynamicManhattan`：add / erase / size / alive / point / max_pair / farthest、`sum_manhattan_pairs` / `sum_chebyshev_pairs`（Int128）、**追加：`sum_manhattan_from_each`（各点から全点への距離の和）/ `manhattan_median`（距離の和が最小の格子点）**。**最終照合で追加：`AbsSum`（Σ w\|x − x_i\| を O(log N)・重み付きの中央値）/ `abs_sum_table(xs, lo, hi)`** | geo_base | 典型90 036、ABC351 E、ABC366 E |
| ☆ `07_manhattan_mst` | **作成済み（Phase 4）**。`manhattan_mst(ps)`（`ManhattanMST{cost, edges}`。edges は重みの小さい順。O(N log N)）、`manhattan_mst_candidates(ps)`（候補の辺 4N 本以下。`ManhattanEdges{us, vs, ws}` は graph_grid_v2 の `graph_from_edges(N, e.us, e.vs, e.ws, UNDIRECTED, 0)` にそのまま渡せる）、**`chebyshev_mst(ps)`（追加）**。**最終照合で追加：`manhattan_nearest_each(ps)`（各点の最近点 `ManhattanNear{j, d}`。8 方向）** | geo_base | ABC283 F |

### 8. 軸に平行な長方形・直方体

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ☆ `08_rect` | **作成済み（Phase 4）**。`Rect(x1, y1, x2, y2)`（半開）：width / height / area / empty / contains(p or r2)、`rect_from_points` / `rect_from_closed(L, D, R, U)` / `lattice_points_closed`、`intersect` / `overlaps` / `bounding`、`union_area(rects)`（走査＋非再帰のセグメント木。O(N log N)）、`BoxK<K>`（`Box3`）：volume（Int128。辺が 2×10^9 なら K ≤ 4）/ empty / contains、`intersect` / `box_intersection_volume`、`quadrant_split` / `sign_grid_split` / `split_by_x_axis` / `split_by_y_axis`、`diagonal_lattice(r, c, ±1)`（`DiagonalLattice{exist, first, last, count}`）。**最終照合で追加：`rect_from_ranges`、閉じた長方形としての `where_closed` / `intersection_dim`（Rect・BoxK。隣り合う・角で接するの区別）/ `boundary_position`、`to_polygon`、`rect_fits`（回して入るか。厳密）** | geo_base | ABC343 E（立方体の重なり）、ABC318 B（シートの重なり）、ARC029 B、ARC076 E、ABC312 E |

### 9. 半平面交差

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ☆ `09_halfplane` | **作成済み（Phase 4）**。`HalfPlane(a, b, c)`（a x + b y ≤ c）/ `HalfPlane(l)`（直線の左側）/ `HalfPlaneD`、`contains(h, p)`、`halfplane_intersection(hs or lines, box)`（凸多角形、反時計回り、PointD。空・線分・1 点も。整数なら判定は 128 ビットで厳密、交点を丸めない。**有界でないときだけ** |x|, |y| ≤ box で切る）、**`halfplane_unbounded(hs or lines)`（追加）** | geo_base, **line（最終照合で 03_polygon が使うようになった）**, **polygon**（arg_sort は使わなかった） | ABC251 G（contains(h, p)） |

### 10. 3 次元

| ファイル | 責務・提供する機能 | 依存 | 代表的な問題 |
|---|---|---|---|
| ☆ `10_geo3d` | **作成済み（Phase 4）**。`Point3` / `Point3D`（2 次元と同じ演算子）：dot / cross / norm2 / dist2 / len / dist / dot128 / norm2_128 / dist2_128 / to_d / midpoint / unit、`volume6`（Int128 で厳密）/ `tetra_volume` / `orient3d` / `coplanar` / `collinear` / `triangle_area` / `in_tetrahedron`、`angle_between` / `elevation_angle`（_deg）、`rot90_x` / `rot90_y` / `rot90_z` / `rotation24` / `rotated` / `rotated_deg`、`Line3` / `Segment3`（on_line / on_segment / is_parallel / is_orthogonal / intersects / projection / reflection / dist）、`Plane3(a, b, c, d)` / `plane_from_points`（side / on_plane / dist / projection / reflection / cross_point / cross_line / angle_between）、`Sphere`（contains / sphere_relation / cross_points / intersection_volume / sphere_cap_volume）、`read_point3` など / `Point3Hash`。**最終照合で追加：`line_relation`（`LINE3_CROSS` / `PARALLEL` / `SAME` / `SKEW`）/ `cross_point(l, m)`・`(l, s)` / `cross_points(s, t)` / `FracPoint3`・`cross_point_frac` / `direction_key`・`line_direction_key`（3 次元）/ `LineKey3`・`line_key`、`line_plane_relation` / `intersects`・`dist`（直線・線分と平面）/ `plane_key` / `plane_relation` / `dist(平面, 平面)`、`sphere_plane_relation`（256 ビットで厳密）/ `Circle3D`・`cross_circle(球, 平面 or 球)`、`ball_intersects` / `ball_contains` / `ball_dist`、`rotated` の型を自由に** | geo_base | 典型90 018、ABC301 G |

### △（型）

| 型 | 置き場所 |
|---|---|
| 実数の二分探索・三分探索で答えを求める（半径を二分探索して円の交点を調べる、など） | 05_circle の使い方関数（ABC157 F の形） |
| 角度の区間を並べて貪欲に選ぶ | 01_arg_sort の使い方関数（ABC225 E の形） |
| 凸多角形の上の尺取り | 04_convex_polygon の使い方関数（ABC250 F の形）。最終照合で 03_polygon の `PolygonPrefix` を追加 |
| 小数の座標を 10^k 倍して整数で解く | 00_geo_base の使い方関数（parse_fixed） |

### 対象外（決定済み）

- 2 次元の集計（累積和・いもす法、長方形の中の点の数を数える 2D BIT）：グラフ・グリッドと同じく対象外（データ構造）
- 座標圧縮：汎用テンプレートとして別に持つ
- 1 次元の区間の分割（旧 RangeSplit の sign_split / cut_by_range）：区間のテンプレートの範囲（5-4）

## 2-P2. Phase 2（基盤）の実装メモ（2026-10-08）

### 作ったファイル（`Desktop\claude\geometry_v2\`）

| ファイル | 内容 | テスト |
|---|---|---|
| `00_geo_base.cpp` | 幾何基盤 | 演算（128 ビットの愚直計算と比較）、ccw / ccw_detail（\|x\| ≤ 6 の全退化ケース、\|x\| = 10^9 で a + k·d の形に作った一直線、1 だけずらした点、実数の点は回転・拡大・平行移動しても向きが変わらないこと）、回転（複素数の積と比較。90 度の倍数は厳密に一致）、角度（std::arg・acos と比較）、direction_key（外積・内積での判定と一致）、isqrt（境目の前後と 9.2×10^18 まで）、parse_fixed（7 通りの書き方：末尾の 0 なし・余分な 0・先頭の 0・".25"・"5."・"+"）、入力関数・Int128 の文字列化・ハッシュ。Release で約 278 万回の照合 |
| `update_bundles.ps1` / `update_bundles.bat` | graph_grid_v2 のものをそのままコピー（ガードが `GG_` ＋ファイル名なので変更不要。`GG_GEO_BASE` を正しく見つけることを確認） | － |
| `dynalist_items_phase2.txt` | Dynalist に貼る項目（説明・メモ・パス） | － |

確認したこと：g++ 13（-O2 -Wall -Wextra -Wshadow -pedantic）で警告なし・全テスト通過、clang 18（-pedantic-errors）で警告なし、C++17 / 20 / 23 でコンパイル・実行、`-D_GLIBCXX_DEBUG -fsanitize=address,undefined` で全テスト通過、**Real を double に変えた版**でも全テスト通過（double の精度で区別できる範囲のケースだけ）、Main_base.cpp（ACL 込み）に**グラフ・グリッド第二世代の 46 本と 00_geo_base の本体**を 2 通りの順で貼り、gnu++2b / c++17・LOCAL あり / なしで g++・clang の両方でコンパイル・実行（名前の衝突なし。`debug(p)`、std::midpoint がある C++20 でも `midpoint(a, b)` が幾何の方を呼ぶこと、10_grid_transform の `rotate_point_cw`・`manhattan_to_chebyshev` と一緒に使えることも確認）。

参考：この貼り付けテストで -Wshadow をつけると、グラフ・グリッドの 05_tree_classic と 08_min_cost_flow のローカル変数 `INF` が Main_base の `INF` を隠すという警告が出る（動作は正しい。幾何のファイルからは警告なし）。

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）

| 処理 | 規模 | 時間 |
|---|---|---|
| ccw / ccw_detail | 2×10^6 回 | 6 / 7 ms |
| rotated（毎回違う角。long double の sin・cos） | 2×10^6 回 | 約 400 ms |
| isqrt / direction_key | 2×10^6 回 | 25 / 220 ms |
| read_points / read_points_fixed(N, 4) | N = 2×10^5 | 30〜50 / 20 ms |

### 実装で決めたこと

- `Int128` は `__extension__ typedef __int128 Int128;`（-pedantic でも警告が出ない書き方）
- **整数の点を実数から作る・整数の点に実数を掛ける／割るとコンパイルエラー**（`Point(1.5, 2)`、`p * 0.5` など。黙って切り捨てないため）。Point → PointD は `to_d(p)` か `PointD(p)`（明示）だけで、勝手には変わらない
- `cout << p` は `(x, y)` の形（グラフ・グリッドの `operator<<` と同じくデバッグ表示用）。答えの出力は `p.x << ' ' << p.y`
- 実数の点の `ccw` は、ベクトルの長さの積 \|u\|\|v\| に合わせた相対誤差（sin θ ≤ 1e-9 なら一直線）
- `rotated_deg` は 90 度の倍数なら `rot90` で厳密に回す（-0.000000 のような出力を防ぐ）。`angle_of` は -0 と 2π に丸められた値を 0 にする
- `parse_fixed` は、小数点以下に k 桁を超える 0 でない桁・数として読めない文字があると、メッセージを出して止める
- `Where`（OUTSIDE / ON_BORDER / INSIDE）と `CcwDetail` は基盤に置いた（03_polygon・04_convex_hull・05_circle で共通に使う）
- 線分・直線・円の型は基盤に置き、アルゴリズムは 02_line・05_circle に置く（設計どおり）。入力用に `read_segment` / `read_line` / `read_circle` を足した
- 設計になかったが足したもの（迷ったら入れる）：`rot180`、`dot(o, a, b)` / `cross(o, a, b)`、`collinear`、`complex_mul` / `conjugate`、`midpoint` / `unit` / `from_polar`、`isqrt_ceil`、`read_point_fixed` / `read_points_fixed`、点どうしの `real_eq`
- テストとデモの形はグラフ・グリッドと同じ（`geob_test::random_test()`（Release は 8 倍の量）・`benchmark()`・`geob_demo::demo()`・`setup_console()`）。デモに ABC259 B・ABC168 C・ABC226 D・ABC198 C・ABC191 D の使い方を入れた

- **Phase 3 での変更**：`geo_internal::sgn_cross` / `geo_internal::sgn_dot`（外積・内積の符号。整数なら厳密、実数なら |u||v| に合わせた誤差）を基盤に移した（02_line・03_polygon などで共通に使うため。利用者向けの関数は変えていない）
- **Phase 4 での修正**：`angle_between(a, b)` で片方が 0 ベクトル、もう片方の成分が両方負のとき π を返していた（atan2(+0, -0) = π）→ 説明どおり 0 を返すように直した（テストを追加。全ファイルの同梱部分も更新済み）

### 次回以降の作業のしかた

- 本体の正は各 .cpp の「ここからテンプレート本体（コピー範囲）」の中。基盤を直したら `update_bundles.bat` を実行する（Phase 3 のファイルができてから意味を持つ）
- Phase 3 のファイルも、同梱する基盤は 00_geo_base.cpp から取り出して入れ、最後に `update_bundles.bat -Check` で一致を確かめる（こちらでは同じ処理を Python で確認する）

## 2-P3. Phase 3（★）の実装メモ（2026-10-08）

### 作ったファイル（`Desktop\claude\geometry_v2\`）

| ファイル | 同梱している本体 | テスト（主な愚直解。Release の照合回数） |
|---|---|---|
| `01_arg_sort.cpp` | geo_base, arg_sort | atan2 と外積・内積で作った愚直な比較・並べ替え・範囲の数（端を含む / 含まない、0 度をまたぐ）、角の比較は回転・拡大・入れ替えで等しい角を作る（85 万回） |
| `02_line.cpp` | geo_base, line | 端点が相手の線分の上か＋真の交差で作った愚直な交差判定、有理数の交点、重なりの両端、距離は場合分けした有理数、マスは 4000 点のサンプリング、実数の型は回転・拡大しても判定が同じこと（130 万回） |
| `03_polygon.cpp` | geo_base, polygon | 角度の和で巻き数を求める内外判定、全探索の単純性・凸性（星形の多角形・とげ・同じ点）、三角形分割の重心、格子点を全部数える（46 万回） |
| `04_convex_hull.cpp` | geo_base, polygon, convex_hull | 「ほかの 3 点の三角形に含まれない点」「支える直線がある点」で作った凸包の頂点の集合、全組の最遠点、O(N) の内外判定、頂点の多い凸多角形（12 万回） |
| `05_circle.cpp` | geo_base, line, circle | 線分と円周の交点の数は平方根を使わない別の式で、位置関係は実数の距離で、交点・接点・接線・外接円などは条件を満たすか、格子点は全部数える。10^9 の座標で外接・内接・接線を作る（82 万回） |
| `06_point_set.cpp` | geo_base, line, arg_sort, point_set | 全部の組・3 点の組・4 点の組の全探索（長方形・正方形・平行四辺形・台形）（9 万回） |
| `07_manhattan.cpp` | geo_base, manhattan | 全部の組の距離、追加・削除は 1 回ごとに全部数え直し、中央値は範囲を全部調べる（30 万回） |

確認したこと（00_geo_base を含む 8 ファイルすべて）：g++ 13（-O2 -Wall -Wextra -Wshadow -pedantic）で警告なし・全テスト通過、clang 18（-pedantic-errors）で警告なし・全テスト通過、C++20 / C++23 でコンパイル・実行、`-D_GLIBCXX_DEBUG -fsanitize=address,undefined` で全体（テスト・速度計測・デモ）が通過、**Real を double に変えた版**でも全テスト通過、Main_base.cpp（ACL 込み）に**グラフ・グリッド第二世代 46 本と幾何 8 本**の本体を 2 通りの順で全部貼って gnu++2b / c++17・LOCAL あり / なしでコンパイル・実行（幾何から警告なし）、幾何の各ファイルを 1 本ずつ Main_base に貼っても動くこと、同梱部分が持ち主のファイルと全部一致（update_bundles -Check と同じ確認）。

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）

| 処理 | 規模 | 時間 |
|---|---|---|
| arg_sort_idx / ArgSorted + count × N / count_in_open_half | N = 2×10^5 | 35 / 130 / 70 ms |
| intersects / cross_point / dist(線分, 点) / line_key | 2×10^6 回 | 25 / 80 / 30 / 320 ms |
| contains（N = 2×10^5 の多角形）× 100 / is_convex / is_simple（N = 2000） | | 62 / 5 / 7 ms |
| convex_hull（ランダム / 全部が頂点）/ farthest_pair / convex_contains × N | N = 2×10^5 | 32 / 31 / 32 / 47 ms |
| circle_relation / cross_points(円, 円) / cross_points(円, 線分) | 10^6 回 | 10 / 32 / 40 ms |
| lattice_points_in_circle（ABC191 D の最大：半径 10^9、scale 10^4） | | 1 ms |
| closest_pair / pairs_within | N = 2×10^5 | 50〜60 / 115 ms |
| count_collinear_triples / count_lines_with_at_least（K = 2）/ count_triangles_by_angle / max_angle | N = 2000 | 116 / 700 / 460 / 340 ms |
| 長方形 + 正方形 + 平行四辺形 + 台形（座標 ±40 に 2000 点） | N = 2000 | 680 ms |
| ManhattanFarthest / DynamicManhattan（追加・削除・質問 × N）/ 距離の和 3 種 | N = 2×10^5 | 1 / 215 / 106 ms |

### 実装で決めたこと

- **組み立て方**：本体ごとに書いた元ファイルから、依存する本体（持ち主のファイルのブロック）を同梱して .cpp を組み立てる。テスト・デモの共通部品（CHECK・乱数・許容誤差・to_s・SHOW など）も同じ形で差し込む（作業用の仕組み。できた .cpp はグラフ・グリッドと同じ形）
- 偏角の比較・凸包・内外判定・数え上げはすべて整数（必要なところは Int128）。atan2 は使わない
- `count_triangles_by_angle` / `max_angle` は、点ごとに偏角の順に並べてから尺取り（二分探索より 2 倍速い）。同じ向きの点が前にあるときは範囲から外す
- `pairs_within` は 1 辺 ⌈√D2⌉ のマスに分け、同じマスと隣の 4 マスだけを見る（同じ組を 2 回見ない）
- `count_lines_with_at_least` は K 点以上の直線だけキーを作る（K = 2 のときだけ全部の直線のキーを作るので遅め）
- `convex_contains` は、凸包に一直線の頂点（keep_collinear）があっても正しく判定する（最初・最後の辺の向きの最も遠い頂点を二分探索）
- `lattice_points_in_convex` を追加（ピックの定理は凸包が点・線分のときに使えないため。実は一直線の点の並びでも同じ式で正しくなる）
- 2 円の交点・接線・接点は、半弦などの 2 乗を Int128 で正確に求めてから平方根をとる（接する直前でも精度が落ちない）
- 実数の円と線分の交点の判定は、桁落ちする C = |w|² − r² の誤差の目安を |w|² と r² の大きい方にした（接するのに 2 点と判定する誤りを直した）
- 名前：内外は `OUTSIDE / ON_BORDER / INSIDE`、三角形は `TRIANGLE_*`、直線の関係は `LINE_*`、円の関係は `CIRCLE_*`（Windows のマクロ・よくある名前とぶつからない）
- 設計から変えたもの：02_line の「距離の 2 乗を分数で」→ `dist2_cmp`（整数との比較）。追加したもの：04 の `lattice_points_in_convex`、07 の `sum_manhattan_from_each` / `manhattan_median`

## 2-P4. Phase 4（☆）の実装メモ（2026-10-08）

### 作ったファイル（`Desktop\claude\geometry_v2\`）

| ファイル | 同梱している本体 | テスト（主な愚直解。Release の照合回数） |
|---|---|---|
| `04_convex_polygon.cpp` | geo_base, line, polygon, convex_hull, halfplane, convex_polygon | 全部の頂点を調べた最大・接線（同じ向きは近い方）、有理数で作った共通部分（候補の点の凸包）、半平面交差で切った形・左右の面積の和、全部の和の凸包（ミンコフスキー和）、各辺から最も遠い頂点（幅）。一直線の頂点・退化した形・keep_collinear の一直線の形を含む（42 万回） |
| `05_circle_cover.cpp` | geo_base, line, circle, circle_cover | 1・2・3 点で決まる円を全部試す（最小包含円）、2 点を通る半径 r の円の中心を全部試す（最大の点の数）、座標 10^9 付近に集まった点（原点の近くへ平行移動した答えと同じ）（10 万回） |
| `05_circle_area.cpp` | geo_base, line, polygon, circle, circle_area | 適応シンプソンの数値積分（被積分関数の折れ目で区切る）、拡大しても比が同じ、接する直前の細い共通部分（25 万回） |
| `06_congruence.cpp` | geo_base, congruence | 2 点の行き先を全部試し、変換を有理数で確かめる愚直解、多角形は始点・向き・裏返しを全部試す、重ね方の精度（座標 10^9 付近・重心のすぐ近くの点・2000 点）（32 万回） |
| `07_manhattan_mst.cpp` | geo_base, manhattan_mst | 全部の組を使うプリム法（同じ点・座標 2 以内で重なりが多いもの・10^9）（4 万回） |
| `08_rect.cpp` | geo_base, rect | マス・格子点を全部数える、大きい座標は座標圧縮したマスごとに数える（44 万回） |
| `09_halfplane.cpp` | geo_base, polygon, halfplane | 全部の交点のうち全部を満たすものの凸包（有理数）、箱を大きくして変わるか（有界か）、凸多角形の辺＋余分な直線、隣り合う向きの差が 10^-19 ラジアン、頂点が 4×10^18 の有界な三角形（30 万回） |
| `10_geo3d.cpp` | geo_base, geo3d | 行列式の余因子展開、重心座標、カハンの式（角）、クラメルの公式（3 次元の交差）、正規方程式（直線どうしの距離）、三分探索（線分どうしの距離）、断面の数値積分（球の共通部分）、ほぼ平行な平面の交線、192 ビットの積の差（373 万回） |

確認したこと（00_geo_base を含む 16 ファイルすべて）：g++ 13（-O2 -Wall -Wextra -Wshadow -pedantic）で警告なし・全テスト通過、clang 18（-pedantic-errors）で警告なし・全テスト通過（clang は引数の評価順が違うので別の乱数の例になる）、C++20 / C++23 でコンパイル、`-D_GLIBCXX_DEBUG -fsanitize=address,undefined` で全体（テスト・速度計測・デモ）が通過、`-fsanitize=float-cast-overflow` を足して Release の量のテストが通過、**Real を double に変えた版**でも全テスト通過、Main_base.cpp（ACL 込み）に**グラフ・グリッド第二世代 46 本と幾何 16 本**の本体を 2 通りの順で全部貼って gnu++2b / c++17・LOCAL あり / なしでコンパイル・実行（幾何から警告なし。graph_grid_v2 の kruskal に manhattan_mst_candidates を渡す使い方も確認）、幾何の各ファイルを 1 本ずつ Main_base に貼っても動く、同梱部分 45 個が持ち主のファイルと全部一致（update_bundles -Check と同じ確認）、索引・Dynalist の項目に書いた関数名がすべてテンプレートの中にある。さらに、**作ったものを見ていない別のエージェントに本体を読ませ、全探索・厳密計算で確かめさせた**（見つかった 4 件と説明の誤り 2 件は下の表のとおり直した）。

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）

| 処理 | 規模 | 時間 |
|---|---|---|
| extreme_point × N / tangents × N | 頂点 138344 の凸多角形、N = 2×10^5 | 48 / 240 ms |
| convex_intersection + minkowski_sum + convex_width + convex_cut | 頂点 138344 | 72 ms |
| min_enclosing_circle / max_points_in_circle | N = 2×10^5 / N = 2000 | 10〜20 / 60 ms |
| intersection_area(円, 円) × 10^6 / intersection_area(円, 多角形 N = 2×10^5) × 10 | | 210 / 180〜230 ms |
| congruence × 2 / polygon_congruence | N = 2×10^5 | 320〜345 / 23 ms |
| manhattan_mst / chebyshev_mst | N = 2×10^5 | 280 / 270〜290 ms |
| union_area | N = 2×10^5（座標 10^9 のランダムな長方形） | 250〜340 ms |
| halfplane_intersection（整数 / 実数） | N = 2×10^5 | 55〜70 / 90〜105 ms |
| orient3d × 2×10^6 / dist(線分, 線分) × 10^6 / intersects × 10^6 / 球の共通部分 + 球と直線 × 10^6 / rotated × 10^6 | | 12 / 70 / 10 / 58 / 204 ms |

### 実装で決めたこと

- **04_convex_polygon**：`tangents` は「p から見て右・左」の接点の番号（見える辺は h[left] → … → h[right]）。`convex_intersection` は 09_halfplane で求める（そのため 09 を同梱）。`convex_cut` は点・線分・一直線（keep_collinear の一直線の形も）を両端の線分として扱う。実数の形を切るときの交点は辺の上で求める（判定の誤差と食い違って空の交点を読むことがないように）。`convex_intersects_line` を追加
- **05_circle_cover**：最小包含円は Welzl（固定の種の乱数で並べ替え、期待 O(N)）。円の内側の判定の誤差は「半径 × GEO_EPS」と「中心を Real で表すときの丸め（座標 × 64ε）」の大きい方（座標の大きさ × GEO_EPS にすると、座標 10^9 で 1 だけ離れた点を同じ点とみなす。半径だけにすると Real が double のとき判定が食い違って答えが壊れる）。`max_points_in_circle` は 1 点を円周の上に置いて中心の偏角の区間を走査
- **05_circle_area**：弓形は角が小さいとき級数。円と多角形は中心と各辺の三角形に分け、辺と円周の交点は桁落ちしない二次方程式の解で。2 円は半弦の 2 乗を Int128 で。`circle_cut_area` を追加
- **06_congruence**：重心を原点に（n 倍して整数のまま）、偏角の順に（長さの 2 乗・次との外積・内積）の列を作って KMP で巡回一致。重ね方の角は、重心から最も遠い点の向きを整数のベクトルで比べて求める。多角形の版 `polygon_congruence` を追加
- **07_manhattan_mst**：候補の辺は 45 度ずつに分けた右半分の 4 方向だけ BIT で求める（並べ替えは 4 回を共有）。`chebyshev_mst` を追加。はじめからつながっている点があるときは graph_grid_v2 の `kruskal(G, uf)`（使えない辺がある問題には候補の辺は使えない）
- **08_rect**：`union_area` は非再帰のセグメント木（覆われた回数・覆われた長さ・節点の長さを 1 つの構造体に）。32 ビットの節点にすると 25% 速いが、座標の範囲を超えたときに黙ってあふれるので long long のままにした。旧テンプレート（RangeSplit・RectangleDiagonalQuery）の 2 次元の分割と傾き ±1 の格子点を取り込んだ
- **09_halfplane**：整数の半平面なら偏角の並べ替え・交点を残すかの判定はすべて 128 ビットで厳密（|a|, |b| ≤ 2×10^9、|c| ≤ 4×10^18）。頂点だけ実数（交点を丸めない。旧テンプレートの不具合「高-3」を解消）。**境界の向きが 180 度以上あくときだけ箱で切る**（有界なら頂点が 10^18 より遠くても正しい）。`halfplane_unbounded` を追加。依存は設計の line・arg_sort ではなく polygon（結果に area などをそのまま使えるように）
- **10_geo3d**：`Point3` / `Point3D` は 2 次元の Point と同じ演算子（整数の点に実数を掛けるとコンパイルエラー）。体積・向き・交差は Int128 で厳密。`dot` / `norm2` / `dist2` は long long で返す（座標の差が 1.7×10^9 を超えるときは `_128` の版。使い勝手を優先）。`Plane3` は a x + b y + c z = d（整数なら d は Int128）。2 平面の交線は向きを 128 ビット、点を 192 ビットで誤差なく求めてから割る（ほぼ平行な平面でも正しい）。回転はすべて右ねじの向き。球の位置関係の値は CIRCLE_* と同じ
- **テストの許容誤差**：実数の判定は相対誤差 10^-9 なので、座標 10^9 では「1 だけ離れた頂点」を 1 つにまとめることがある。実数の版と整数の版の面積を比べるときは（座標の大きさ）² × 10^-8 まで許した。始点が「(x, y) 最小の頂点」の多角形は、同じ x の頂点が丸めで入れ替わることがあるので、始点のずれを許して比べた

### 独立したレビューで見つかって直したもの

| ファイル | 見つかったこと | 直し方（回帰テストを追加済み） |
|---|---|---|
| 06_congruence | 重ね方の角・平行移動の精度が悪い（重心のすぐ近くの点を基準にしていた。n = 2×10^5 でずれ 783） | 重心から最も遠い点を基準にし、角を整数のベクトルから求める（ずれ 3×10^-11） |
| 10_geo3d | 整数の平面の `cross_line` が、ほぼ平行な平面で大きくずれる（向きが 0.88 rad） | 向きは 128 ビット、点は 192 ビットで誤差なく求めてから割る（Real が double でも打ち消し合わないように下位 64 ビットを ±2^63 にそろえる） |
| 05_circle_cover | 座標 10^9 付近で 1 だけ離れた点を同じ点とみなす（半径 0.5 が 0 に） | 誤差の目安を「半径」と「座標の丸め」に |
| 04_convex_polygon | `convex_cut` に keep_collinear の一直線の凸包を渡すと同じ点が 2 回出る | 一直線なら両端の線分として扱う |
| 09_halfplane（説明の誤り） | 有界でも箱で切っていた（頂点が 2×10^18 より遠い有界な三角形が切られる） | 有界でないときだけ箱を足す |
| 08_rect（説明の不足） | `BoxK` の体積が K ≥ 5 で 128 ビットを超えうる | 説明に「辺が 2×10^9 なら K ≤ 4」と書いた |

（レビューの前に自分で見つけて直したもの：基盤の `angle_between` の 0 ベクトル、08 の `union_area` を非再帰に（約 1.4 倍速い）、09 の `halfplane_unbounded` の並べ替えを整数で厳密に（10^-19 ラジアンの差を実数の偏角で取り違えていた）、05_circle_area の速度計測の long long へのはみ出し、07 のデモの説明の誤り）

### 設計から変えたもの・足したもの

- 足したもの（迷ったら入れる）：04 `convex_intersects_line`、05_circle_area `circle_cut_area`、06 `polygon_congruent` / `polygon_congruence`、07 `chebyshev_mst`、09 `halfplane_unbounded`、10 の線分・平面・球・回転一式（設計では「点と直線・平面の距離、角度」だけだった）
- 依存：04_convex_polygon は 09_halfplane を同梱（設計は convex_hull, line）、09_halfplane は polygon を同梱（設計は line, arg_sort）

---

## 3. 依存関係

```
00_int_common（2026-10-10 追加。すべてのファイルが 00_geo_base の前に同梱する。数学第二世代と共通）
00_geo_base（すべてのファイルが同梱する）
├─ 01_arg_sort
├─ 02_line
├─ 03_polygon         ← 02_line（最終照合で追加）
├─ 04_convex_hull     ← 02_line, 03_polygon
├─ 04_convex_polygon  ← 02_line, 03_polygon, 04_convex_hull, 09_halfplane
├─ 05_circle          ← 02_line
├─ 05_circle_cover    ← 02_line, 05_circle
├─ 05_circle_area     ← 02_line, 03_polygon, 05_circle
├─ 06_point_set       ← 02_line, 01_arg_sort
├─ 06_congruence
├─ 07_manhattan
├─ 07_manhattan_mst
├─ 08_rect
├─ 09_halfplane       ← 02_line, 03_polygon
└─ 10_geo3d
```

- ← の右に書いたファイル（と基盤）を同梱する。例：`04_convex_hull` には geo_base と line・polygon が入る（凸包の面積・ピックの定理・多角形の位置関係を 1 回のコピーで使える）
- 最終照合で 03_polygon が 02_line を同梱するようにした（多角形と線分の位置関係に線分の判定を使うため）。03_polygon を同梱するファイル（04_convex_hull・09_halfplane・04_convex_polygon・05_circle_area）にも 02_line が入る
- グラフ・グリッド第二世代には依存しない（同梱しない）。組み合わせ方は索引（4 章「グラフ・グリッドと組み合わせる」）と使い方関数のコメントに書く

### 作る順番

- Phase 2（基盤）：00_geo_base（＋ update_bundles のコピー、Dynalist の項目）
- Phase 3（★）：01_arg_sort、02_line、03_polygon、04_convex_hull、05_circle、06_point_set、07_manhattan（7 本）
- Phase 4（☆）：09_halfplane、04_convex_polygon、05_circle_cover、05_circle_area、06_congruence、07_manhattan_mst、08_rect、10_geo3d（8 本。04_convex_polygon が 09_halfplane を同梱するので 09 を先に作った）
- 最終照合：旧ファイルを 1 機能ずつ、過去の解答、AtCoder の典型（問題文を確かめながら）、Dynalist の項目と突き合わせる

合計 16 ファイル（★ 8・☆ 8）。

---

## 4. 道具選択体系

**16 本すべてを `geometry_v2\dynalist_tool_index.txt` / `.opml` に作った**（最終照合で、2 つの図形の関係を引く「▼ 図形の位置関係（早見）」と、追加した関数・見つけにくかった関数の行を足した。7-7）（Phase 3 で ★、Phase 4 で ☆ を足した。☆ の分は大分類「凸包・凸多角形」「円」「点の集合」「マンハッタン距離」に足し、「軸に平行な長方形・直方体」「半平面交差」「3 次元」を新しく作った。索引に書いた関数名がテンプレートの中にあることを確かめた）。下は Phase 1 の案（全体の形）。

グラフ・グリッドと同じく、Dynalist の「00_典型発想インデックス」の中に 1 項目（▼ 幾何道具選択（第二世代））として置く。辿り方は「入力の形（大分類）」→「聞かれていること（小分類）」→「ファイル名と関数」。各 Phase で関数が確定したら、`geometry_v2\dynalist_tool_index.txt` / `.opml` として完成させる。

```
▼ 幾何道具選択（第二世代）
	▼ 使い方
		1. 入力の形で大分類を選ぶ（点・ベクトル／角度／直線・線分／多角形／凸包／円／点の集合／マンハッタン距離／長方形／その他）
		2. 聞かれていることで小分類を選ぶ
		3. 書いてあるファイルを貼る（依存するファイルは同梱済み。1回のコピーで動く）
		判定は整数で厳密（誤差を考えなくてよい）。交点など答えが整数にならないものだけ PointD（実数）
	▼ 点・ベクトル　基盤：00_geo_base
		読み込み → read_point() ／ read_points(N)（小数なら read_points_d）
		座標が小数（小数点以下 k 桁まで）で、誤差なく判定したい → parse_fixed(s, k) ／ read_fixed(k) で 10^k 倍の整数に
		左に曲がるか・右か・一直線か → ccw(a, b, c)
		一直線のとき、どの順に並んでいるか → ccw_detail(a, b, c)
		内積・外積・距離 → dot ／ cross ／ dist2（整数のまま）／ dist
		90 度回転 → rot90(p) ／ rot90(p, 中心, k)　任意の角で回す → rotated(p, θ) ／ rotated_deg
		向き・傾きが同じものをまとめる → direction_key(v) ／ line_direction_key(v)
		平方根を誤差なく（整数の切り捨て） → isqrt(n)
	▼ 角度　01_arg_sort
		偏角の順に並べる → arg_sort(ps, 原点) ／ arg_sort_idx
		ある向きの範囲（扇形）に入る点の数 → ArgSorted S(ps, 原点);　S.count(l, r)（2π をまたいでもよい）
		2 つの角の大きさを誤差なく比べる → angle_less(a1, b1, a2, b2)
		角度の区間を選ぶ問題 → 使い方デモ（ABC225 E）
		いくつか選んでベクトルの和の長さを最大にする → 使い方デモ（ABC139 F）
	▼ 直線・線分　02_line
		線分が交わるか → intersects(s1, s2)（端点で触れるのも含む）／ intersects_strict
		交点 → cross_point(l1, l2) ／ cross_points(s1, s2)
		点との距離・垂線の足・線対称な点 → dist ／ projection ／ reflection
		同じ直線の上にある点をまとめる・直線の数を数える → line_key(a, b)
		垂直二等分線 → bisector(a, b) ／ bisector_key(a, b)
		線分の上の格子点の数・線分が通るマスの数 → lattice_points_on(s) ／ cells_crossed(s)
	▼ 多角形・三角形　03_polygon
		面積 → area2(poly)（整数の 2 倍）／ area(poly)
		点が内側か（凸でなくてよい） → contains(poly, p)
		凸か（凸四角形の判定） → is_convex(poly)
		内部・周上の格子点の数（ピックの定理） → interior_lattice_points ／ boundary_lattice_points
		直角三角形か・鋭角か・鈍角か → triangle_kind(a, b, c)
	▼ 凸包・凸多角形
		点の集合を囲む凸多角形 → 04_convex_hull　convex_hull(ps)
		最も遠い 2 点 → 04_convex_hull　farthest_pair(ps)
		点が凸多角形の内側か（クエリが多い） → 04_convex_hull　convex_contains(hull, p)（O(log N)）
		a x + b y の最大（クエリが多い） → 04_convex_polygon　extreme_point(hull, {a, b})
		直線で切る・2 つの凸多角形の共通部分・ミンコフスキー和 → 04_convex_polygon
	▼ 円
		点が円の内側か・2 円の位置関係・円どうしが接するか → 05_circle　contains ／ circle_relation ／ has_common_point
		交点・接点・共通接線 → 05_circle　cross_points ／ tangent_points ／ common_tangents
		3 点を通る円・内接円 → 05_circle　circumcircle ／ incircle
		円の中の格子点の数（中心・半径が小数） → 05_circle　lattice_points_in_circle
		全部の点を囲む最小の円 → 05_circle_cover　min_enclosing_circle(ps)
		半径 r の円で最大何点を囲めるか → 05_circle_cover　max_points_in_circle(ps, r)
		2 円・円と多角形の共通部分の面積 → 05_circle_area
	▼ 点の集合（数え上げ・最近点）　06_point_set
		最も近い 2 点 → closest_pair(ps)
		距離 D 以内の点どうしをつなぐ → pairs_within(ps, D) を Union-Find へ
		一直線に並ぶ 3 点の組・三角形の数 → count_collinear_triples ／ count_triangles
		K 点以上を通る直線の数 → count_lines_with_at_least(ps, K)
		鋭角・直角・鈍角三角形の数 → count_triangles_by_angle(ps)
		長方形・正方形・平行四辺形・台形の数 → count_axis_rectangles ／ count_squares ／ count_parallelograms ／ count_trapezoids
		回転・平行移動で重なるか → 06_congruence　congruent(A, B)
	▼ マンハッタン距離・チェビシェフ距離　07_manhattan
		最大のマンハッタン距離・ある点から最も遠い点 → ManhattanFarthest F(ps)
		点の追加・削除がある → DynamicManhattan
		全部の組の距離の和 → sum_manhattan_pairs ／ sum_chebyshev_pairs
		45 度回転（マンハッタン ↔ チェビシェフ） → manhattan_to_chebyshev(p)
		マンハッタン距離の最小全域木 → 07_manhattan_mst
	▼ 軸に平行な長方形・直方体　08_rect
		共通部分・面積 → intersect(a, b) ／ area
		長方形の和集合の面積 → union_area(rects)
		直方体の重なりの体積 → box_intersection_volume(a, b)
		傾き ±1 の直線の上の格子点 → diagonal_lattice(rect, c, ±1)
	▼ その他
		不等式（半平面）をすべて満たす領域 → 09_halfplane
		3 次元 → 10_geo3d
		答えを実数で二分探索・三分探索 → △（05_circle の使い方デモ）
	▼ グラフ・グリッドと組み合わせる
		点どうしの距離を重みにした最短路 → graph_grid_v2 02_dijkstra　dijkstra_dense(N, s, cost)（cost の中で dist を使う）
		点どうしの距離で最小全域木 → graph_grid_v2 03_prim　prim_dense(N, cost)
		距離 D 以内の点をつなぐ → 06_point_set pairs_within ＋ graph_grid_v2 03_union_find
		円どうしがつながっているか → 05_circle has_common_point で辺を張り graph_grid_v2 01_components
		マンハッタン MST の辺から → graph_grid_v2 03_kruskal（graph_from_edges）
		グリッドのマスの回転・グリッドの図形の一致 → graph_grid_v2 10_grid_transform
```

---

## 5. 既存テンプレートとの照合（Phase 1-6）

### 5-0. 対象

Dynalist（dynalist-2026-10-7.opml）の「09_幾何」と、「13_グリッド・座標」「15_道具」に載っている幾何関係の項目が指す旧ファイルと、過去の解答。

| 旧ファイル | Dynalist の項目 |
|---|---|
| GeometryTemplateTest.txt（892 行） | AtCoder 幾何フルテンプレ（整数主導） |
| 259D_2.txt | 同上の「円と円の関係」 |
| 275C.txt | 90度回転（任意中心・左右統合）／範囲判定（動的サイズ対応）／BFS用 近傍ジェネレータ |
| ManhattanTemplateTest.txt | マンハッタン距離 座標変換テンプレ（復元付き） |
| ManhattanMSTemplateTest.txt | マンハッタン距離テンプレ（multiset対応・動的版） |
| RangeSplitTemplateTest.txt | 格子点・整数区間の分割をまとめたテンプレート |
| RectangleDiagonalQueryTemplate.txt | 矩形を斜めに切り取ったときの格子点 |
| 過去の解答 | 266C.txt（凸四角形）、374D.txt（距離の和）、259D_2.txt（円のつながり）、275C.txt（正方形の数） |

旧索引（00_典型発想インデックス）には幾何の ▼ はない。

### 5-1. 新設計にあって、既存にないもの

- 実数の点・直線・線分・円の型と、射影・反射・交点（直線・線分・円）・接線・外接円・内接円
- 多角形（凸でない）の面積・内外判定・ピックの定理・三角形の種類
- 偏角ソートの番号版・始める向きの指定、角の大きさの厳密な比較
- 直線の正規化（LineKey）と、それを使う数え上げ（一直線の組・K 点以上を通る直線・三角形の数）
- 長方形・正方形・平行四辺形の数、合同の判定
- 円の中の格子点、半径 r の円で囲める最大の点数、共通部分の面積
- マンハッタン距離の和・最小全域木、長方形の和集合の面積、直方体の重なり、3 次元
- 小数の入力を整数に直す関数、isqrt、128 ビットでの途中計算（座標 10^9 まで保証）

### 5-2. 既存にあって、新設計になかったもの → 追加した

| 追加先 | 追加した機能 | 参考にした既存 |
|---|---|---|
| 01_arg_sort | 同じ向きの点を 1 つにまとめる（近い点 / 遠い点を残す）`unique_directions` | GeometryTemplateTest.txt（ArgSorter::sort_unique_farthest） |
| 01_arg_sort | 同じ向きの点を原点から近い順（遠い順）に並べる | 同上（ArgSorter の比較） |
| 04_convex_hull | 頂点が 1〜2 個の退化した凸包での内外判定 | 同上（convex_contains） |
| 05_circle | 2 円が「完全に一致」を別の値で返す | 259D_2.txt（circle_relation の 5） |
| 00_geo_base | 整数の点を任意の中心のまわりに 90 度回す（左右） | 275C.txt（rot90） |
| 08_rect | 2 次元の分割（象限・9 分割・x 軸・y 軸）と格子点の数 | RangeSplitTemplateTest.txt |
| 08_rect | 傾き ±1 の直線と長方形の交点の格子点（あるか・最初・最後・個数） | RectangleDiagonalQueryTemplate.txt |

### 5-3. 既存の不具合・使いにくさ（新設計で解消）

不具合の詳細は `Desktop\claude\atcoder_template_review_20260925.md` を参照。

| 既存 | 指摘 | 新設計での解消 |
|---|---|---|
| GeometryTemplateTest.txt（高-3） | 半平面交差の交点を整数に丸めていて結果が壊れる | 09_halfplane は交点を PointD のまま持つ |
| 同上（中-22） | 扇形の判定で 2π をまたぐと点を見落とす | ArgSorted::count は 2π をまたぐ区間をテストで確認する |
| 同上（中-23） | __int128 を使わずオーバーフロー（座標 3×10^4 程度から）。点と線分の距離の 2 乗を整数の割り算で切り捨てて大小が逆転する | すべて \|x\| ≤ 10^9 を保証（128 ビット）。距離は実数か分数で返す |
| 同上 | ArgSorter のメンバ関数が自分の基準点 O を使わない。Point に ==・<・スカラー倍がない | ArgSorted は作るときの原点を使う。PointT に演算子をそろえる |
| 同上 | 最小包含円の誤差が絶対誤差 1e-12（座標 10^9 で小さすぎる） | 相対誤差 |
| 同上 | 最近点対が std::function の再帰で遅い。座標で返す | 番号で返し、std::function を使わない |
| 同上 | ccw が 5 分類（+2 / −2 も正負になる）で、`ccw(...) > 0` と書くと一直線の場合も入る | `ccw` は 3 分類、5 分類は `ccw_detail`（名前つき定数）に分けた |
| 275C.txt | struct 名 Circle・Manhattan が幾何フルテンプレと衝突する。GridRange の x が行で、幾何の x（横）と意味が逆 | 型は基盤で 1 回だけ定義。グリッドの範囲判定はグラフ・グリッド第二世代の 00_grid_base にある。座標の向きを 1-2 に明記 |
| ManhattanTemplateTest / ManhattanMSTemplateTest | 空の集合で max_pair / farthest_point が未初期化の値を返す。消した点をもう一度消すと未定義 | 空なら番号 −1 を返す。erase は 2 回目を無視する |

### 5-4. 体系の外に置くもの（今のまま残す）

- 275C.txt の GridRange・TorusGrid・BFS 用の近傍ジェネレータ：グラフ・グリッド第二世代の 00_grid_base に同等以上の機能がある（グラフ・グリッドの 7-2 で照合済み）
- RangeSplitTemplateTest.txt の 1 次元の分割（sign_split・cut_by_range）：区間の道具なので幾何の体系の外（2 次元の分割は 08_rect に入れる）
- FractionTemplateTest.txt（分数）：数学の道具。幾何で分数が要るところ（交点の厳密な値、距離の 2 乗の比較）は 02_line の中に専用の形で持つ

### 5-5. AtCoder の典型との照合

グラフ・グリッドの最終照合で足りないものが見つかった反省から、Phase 1 の時点で典型と突き合わせた。代表例は最終照合で問題文を確かめ直す。

| 典型 | 代表例 | 担当 |
|---|---|---|
| 点を任意の角度・中心のまわりに回す | ABC259 B、ABC197 D | 00_geo_base rotated / rotated_deg |
| 時計の針・三角関数 | ABC168 C、ABC144 D | 00_geo_base（角の補助・回転） |
| 距離の切り上げを整数で（平方根の誤差を避ける） | ABC198 C | 00_geo_base isqrt・dist2 |
| 小数の座標を整数に直す | ABC191 D | 00_geo_base parse_fixed ＋ 05_circle lattice_points_in_circle |
| 3 点が一直線か・三角形になるか | ABC181 C、ABC224 C | 00_geo_base ccw、06_point_set count_triangles |
| 凸四角形の判定 | ABC266 C | 03_polygon is_convex |
| 直角三角形・鋭角・鈍角 | ABC362 B、ABC033 D | 03_polygon triangle_kind、06_point_set count_triangles_by_angle |
| 向き・傾きを gcd で正規化して数える | ABC226 D | 00_geo_base direction_key |
| K 点以上を通る直線の数 | ABC248 E | 02_line line_key、06_point_set |
| 点と線分の距離 | ARC042 B | 02_line dist(s, p) |
| 偏角ソート・角度の区間 | ABC225 E、ABC139 F、典型90 009 | 01_arg_sort |
| 凸包・ピックの定理 | 典型90 041 | 04_convex_hull ＋ 03_polygon |
| 凸多角形の内外判定（クエリ） | ABC296 G | 04_convex_hull convex_contains |
| 凸多角形を直線で分ける・面積の尺取り | ABC250 F | 04_convex_polygon（使い方デモ） |
| 円がつながるか（2 円の位置関係） | ABC259 D | 05_circle has_common_point ＋ グラフの連結成分 |
| 半径を二分探索して円の交点を候補にする | ABC157 F | 05_circle cross_points・circles_through（△ 型） |
| 最小包含円 | ABC151 F | 05_circle_cover |
| 軸に平行な長方形・正方形の数 | ABC218 D、ABC275 C | 06_point_set |
| 距離 D 以内でつながる | ABC304 C | 06_point_set pairs_within ＋ Union-Find |
| 合同（回転・平行移動で重なるか） | ABC207 D | 06_congruence |
| マンハッタン距離の最大・45 度回転 | 典型90 036 | 07_manhattan |
| 全部の組の距離の和（45 度回転） | ABC351 E | 07_manhattan sum_chebyshev_pairs（偶奇で分ける例を使い方デモに） |
| 直方体の重なり | ABC343 E | 08_rect box_intersection_volume |
| 長方形の重なりの面積 | ABC318 B | 08_rect union_area |
| 3 次元の角度 | 典型90 018 | 10_geo3d |
| 点を回る順の距離の和 | ABC374 D | 00_geo_base dist（全探索は問題ごと） |

### 5-6. 追加しないもの（理由）

| 機能 | 理由 |
|---|---|
| ボロノイ図・ドロネー三角形分割・3 次元凸包 | ABC（D〜F）ではほとんど出ない。必要なら個別に書く |
| 線分の交差の数え上げ（平面走査） | 同上 |
| 円の和集合の面積・多角形の和集合 | 同上 |
| 最小面積の外接長方形・長方形の和集合の周長 | 同上 |
| 長方形の中の点の数（オフライン・2D BIT） | 2 次元の集計なので対象外（データ構造） |

---

## 6. 未決事項と次にやること

### 6-1. ユーザーへの確認（2026-10-08 回答済み）

| 確認したこと | 回答・決定 |
|---|---|
| 一覧（2 章の 16 ファイル）と ★ / ☆ の分け方 | このままでよい。☆ も全部作る |
| MSVC 対策（__int128 が使えない環境への切り替え） | **MSVC ではビルドしない。__int128 は使える** → `Int128 = __int128` をそのまま使い、自作の 128 ビット整数は作らない |
| ガード名を `GG_` のままにして update_bundles をそのまま使う | よい |
| 長方形を半開 [x1, x2) × [y1, y2) にする | よい |
| 座標の型を `Point`（整数）と `PointD`（実数）の 2 つにする／索引を「▼ 幾何道具選択（第二世代）」として置く | 個別には聞いていない。設計書の案のまま進める（変えたくなったら Phase 2 の確認のときに） |

### 6-2. 次にやること

1. ~~6-1 の確認~~（済み）
2. ~~Phase 2：`00_geo_base.cpp`~~（済み。2-P2 章）
3. ~~Phase 3：★ の 7 ファイル~~（済み。2-P3 章）
4. ~~Phase 4：☆ の 8 ファイル~~（済み。2-P4 章）。ユーザーが内容を確認し、Dynalist に `dynalist_items_phase2.txt` / `dynalist_items_phase3.txt` / `dynalist_items_phase4.txt` の項目と索引（`dynalist_tool_index.opml`。16 本ぶん）を入れる
5. ~~最終照合：旧ファイル（5-0 の 7 本）を 1 機能ずつ、過去の解答、典型（5-5 の問題文を確かめながら）、Dynalist の項目~~（済み。7 章）
6. Dynalist への反映はユーザーの環境で行う（MSVC での確認は不要になった）
7. ~~後から機能を足すときは、既存のファイルを直さずに新しいファイル（新しいガード）として足す（0-3 の方針）~~ → 最終照合（2026-10-09）でユーザーの指示により変更：**追加・修正は今までの成果物（既存のファイル・索引・Dynalist の項目・この設計書）に入れる**。新しいファイルは作っていない（7 章）
8. ~~最終照合~~（済み。7 章）

---

## 7. 最終照合（2026-10-09〜10）

ユーザーの指示：「最終照合をお願いします。また図形の位置関係は幾何テンプレートとして特に大事だと考えているので、図形の位置関係がすべて含まれているかという観点での見直しもお願いします。追加修正は今までの成果物に入れてください。」

→ 足したもの・直したものは**すべて今までの 16 ファイル・索引・Dynalist の項目・この設計書に入れた**（新しいファイルは作っていない。6-2 の 7 を変更）。

### 7-1. 図形の位置関係の網羅表

2 つの図形の組ごとに「共有点があるか／入っているか／関係の分類／交点・共通部分／距離」を並べた。**太字は最終照合で足したもの**。整数の型なら判定はすべて厳密（座標 10^9 まで。必要なところは 128・192・256 ビット）。

円は 2 通りに扱う：**円周**（count_cross・circle_relation・has_common_point・cross_points）と**円板**（中身を含む。contains(c, p)・`disk_*`）。名前に `disk_` が付く関数が円板。多角形はすべて**周を含む閉じた領域**（凸でなくてよい・単純多角形）。

#### 2 次元

| 組 | 共有点 | 包含 | 関係の分類 | 交点・共通部分 | 距離 | ファイル |
|---|---|---|---|---|---|---|
| 点・点 | == | － | － | － | dist2 / dist | 00 |
| 点・直線 | on_line | － | side（左右） | projection / reflection | dist / dist2_cmp | 02 |
| 点・**半直線** | **on_ray** | － | － | － | **dist / dist2_cmp** | 02 |
| 点・線分 | on_segment | － | ccw_detail（5 分類） | － | dist / dist2_cmp | 02 |
| 点・円 | contains(c, p)（Where） | 同左 | － | tangent_points | **disk_dist / circumference_dist** | 05 |
| 点・3 点を通る円 | **in_circumcircle**（4 点が同じ円の上 ⇔ ON_BORDER） | 同左 | － | circumcircle | － | 05 |
| 点・三角形／多角形 | in_triangle / contains(poly, p)（Where） | 同左 | － | － | **dist(poly, p) / boundary_dist** | 03 |
| 点・凸多角形 | convex_contains(h, p)（O(log N)） | 同左 | － | tangents | － | 04 |
| 点・長方形 | r.contains（半開）/ **where_closed**（閉） | 同左 | **boundary_position**（周の上の位置） | － | － | 08 |
| 点・半平面／ひし形 | contains(h, p) / in_diamond | 同左 | － | － | manhattan / chebyshev | 09 / 07 |
| 直線・直線 | intersects | － | line_relation（CROSS / PARALLEL / SAME）、is_parallel / is_orthogonal | cross_point / cross_point_frac | dist | 02 |
| 直線・**半直線** | **intersects** | － | is_parallel | **cross_point** | **dist** | 02 |
| 直線・線分 | intersects | － | is_parallel | cross_point（重なるなら両端） | dist | 02 |
| **半直線・半直線** | **intersects** | － | is_parallel | **cross_point** | **dist** | 02 |
| **半直線・線分** | **intersects** | － | is_parallel | **cross_points**（重なるなら両端） | **dist** | 02 |
| 線分・線分 | intersects / intersects_strict | － | **segment_relation**（SEPARATE / TOUCH / CROSS / OVERLAP） | cross_points | dist | 02 |
| 円周・直線／線分／**半直線** | count_cross（厳密） | － | 接する ⇔ count_cross == 1 | cross_points | － | 05 |
| **円板・直線／線分／半直線** | **disk_intersects** | **disk_contains(c, s)** | － | － | **disk_dist** | 05 |
| 円・円 | has_common_point（円周）/ **disk_intersects**（円板） | **disk_contains** | circle_relation（共通接線の数・同じ円） | cross_points / common_tangents / intersection_area | **disk_dist / circumference_dist** | 05 / 05_circle_area |
| 円・多角形（長方形） | **disk_intersects(c, poly)** | **contains(poly, c) / disk_contains(c, poly)** | － | intersection_area | **disk_dist(c, poly)** | 05_circle_area |
| 円・直線の片側 | － | － | － | circle_cut_area | － | 05_circle_area |
| 多角形・直線／線分／半直線 | **intersects(poly, l / s / r)** | **contains(poly, s)** | － | （凸なら convex_cut） | **dist(poly, s)** | 03 / 04 |
| 多角形・多角形 | **intersects(P, Q)** | **contains(P, Q)** | － | （凸なら convex_intersection） | **dist(P, Q)** | 03 |
| 凸多角形・直線 | convex_intersects_line（O(log N)） | － | － | convex_cut | － | 04 |
| 凸多角形・凸多角形 | **convex_intersects** | **convex_contains(h, Q)** | － | convex_intersection / minkowski_sum | （dist(P, Q)） | 04 |
| 長方形・長方形 | overlaps（面積が正）| r.contains(r2) / **rect_fits**（回して入るか） | **intersection_dim**（−1 / 0 角 / 1 辺 / 2 面積） | intersect / union_area / bounding | － | 08 |
| 箱・箱（k 次元） | － | b.contains(p) | **intersection_dim** | intersect / box_intersection_volume | － | 08 |
| 半平面どうし | － | － | halfplane_unbounded | halfplane_intersection | － | 09 |
| **円周の上の弦どうし**（位置の番号で） | **chords_cross** | － | **chord_direction_key**（平行） | **count_crossing_chords / chord_cross_counts / ChordSet** | － | 05 |

長方形と円・多角形の関数をつなぐために **to_polygon(r)** を足した。

#### 3 次元（10_geo3d）

| 組 | 共有点 | 包含 | 関係の分類 | 交点・共通部分 | 距離 |
|---|---|---|---|---|---|
| 点・直線／線分 | on_line / on_segment | － | － | projection / reflection | dist |
| 点・平面 | on_plane | － | side（表裏） | projection / reflection | dist |
| 点・四面体／球 | in_tetrahedron / contains(sp, p) | 同左 | － | － | **ball_dist** |
| 直線・直線 | intersects | － | **line_relation**（CROSS / PARALLEL / SAME / SKEW）、**line_key**（同じ直線のキー） | **cross_point / cross_point_frac** | dist |
| 直線・線分／線分・線分 | intersects | － | is_parallel / is_orthogonal | **cross_point / cross_points** | dist |
| 直線／線分・平面 | **intersects** | － | **line_plane_relation**（CROSS / PARALLEL / ON） | cross_point | **dist** / angle_between |
| 平面・平面 | － | － | **plane_relation**（CROSS / PARALLEL / SAME）、**plane_key** | cross_line | **dist** / angle_between |
| 球面・直線 | － | － | － | cross_points | － |
| **球面・平面** | **sphere_plane_relation**（0 / 1 接する / 2） | － | 同左 | **cross_circle**（Circle3D） | － |
| 球・球 | sphere_relation | **ball_contains** | sphere_relation | **cross_circle** / intersection_volume | **ball_dist** |
| **球体・線分／直線／平面／球** | **ball_intersects** | **ball_contains(sp, s)** | － | － | **ball_dist** |
| 箱・箱 | － | － | **intersection_dim**（08_rect の BoxK） | box_intersection_volume | － |

#### 入れなかった位置関係（理由）

| 組 | 理由・代わり |
|---|---|
| 3 次元の三角形・凸多面体と点・線分（点が三角形の上か、線分と三角形の交差、点が凸多面体の中か） | 広い洗い出し（7-6）の 109 問で 0 問。3 次元の凸包も対象外のまま |
| 球面と線分の交点 | cross_points(sp, Line3) から線分の上のものを選ぶ（ball_intersects で共有点の有無は厳密） |
| 円周と多角形の周の交点の数 | 出題がない。必要なら辺ごとに count_cross(c, s) を足す（頂点で二重に数えない注意） |
| 凸でない多角形を直線で切った形 | 出題がない（凸なら convex_cut） |
| 直線と半平面 | halfplane_intersection に境界を足して扱える |

### 7-2. 旧テンプレート 7 本との照合（1 機能ずつ）

記号：◯＝新しいテンプレートに同じかそれ以上の機能がある、体系外＝グラフ・グリッドなど別の体系にある。

**GeometryTemplateTest.txt**

| 旧 | 新 | 判定 |
|---|---|---|
| Point（+ − / dot / cross / norm2） | Point（00。== < スカラー倍も） | ◯ |
| ccw（5 分類：+1 / −1 / +2 c-a-b / −2 a-b-c / 0） | ccw_detail（同じ値。名前つき定数）・ccw（3 分類） | ◯ |
| intersect_ss | intersects(s1, s2) | ◯ |
| dist2_ps（整数の割り算で切り捨てる不具合） | dist2_cmp(s, p, D2)（誤差なし）/ dist(s, p) | ◯（不具合を解消） |
| dist2_ss | dist(s1, s2) | ◯ |
| Circle / point_circle（−1 / 0 / 1） | Circle / contains(c, p)（Where） | ◯ |
| circle_line_relation（−1 / 0 / 1） | count_cross(c, l)（0 / 1 / 2）/ **disk_intersects(c, l)** | ◯ |
| circle_segment_relation（中心と線分の距離と r。円板としての関係） | **disk_intersects(c, s)**（−dist2_cmp(s, c.c, r²) と同じ値）。Phase 4 までは count_cross（円周との交点）しかなく意味が違った → 円板の関数を足した | ◯（最終照合で追加） |
| convex_hull | convex_hull | ◯ |
| area2（絶対値） | area2（符号付き）/ area | ◯ |
| convex_diameter（2 点の組） | convex_diameter（番号）/ farthest_pair | ◯ |
| convex_contains（−1 / 0 / 1。頂点 1〜2 個も） | convex_contains（Where。頂点 0〜2 個も）/ **convex_contains(h, Q)** | ◯ |
| closest_pair（std::function の再帰・座標で返す） | closest_pair（番号・d2） | ◯ |
| CircleD / in_circle / circle_from（2 点・3 点） | CircleD / contains / circle_with_diameter / circumcircle / **in_circumcircle** | ◯ |
| minimum_enclosing_circle（絶対誤差 1e-12） | min_enclosing_circle（相対誤差） | ◯ |
| HalfPlane / inside / intersection_ll（交点を丸める不具合）/ halfplane_intersection | HalfPlane / contains(h, p) / halfplane_intersection（交点は PointD） | ◯ |
| ArgSorter::operator() / sort_points（原点と同じ点を除く） | arg_less / arg_sort（原点と同じ点は先頭）/ unique_directions | ◯ |
| ArgSorter::sort_unique_farthest | unique_directions(ps, o, true) | ◯ |
| arg_compare | arg_cmp | ◯ |
| arg_lower_bound / arg_upper_bound | ArgSorted::lower_bound / upper_bound | ◯ |
| exists_in_sector / exists_in_sector_circular | ArgSorted::exists / count（0 度をまたいでよい） | ◯ |

**ほかの 6 本**

| 旧 | 新 | 判定 |
|---|---|---|
| 259D_2.txt circle_relation（0 外部 / 1 外接 / 2 交差 / 3 内接 / 4 内部 / 5 一致） | circle_relation（値は共通接線の数 4 / 3 / 2 / 1 / 0、一致 5。**番号が違うので名前で比べる**） | ◯ |
| 275C.txt rot90(p, c, left) | rot90(p, c, k)（k = 1 左・−1 右）/ **Affine::rot90(c) / rot90_cw(c)** | ◯ |
| 275C.txt GridRange / TorusGrid / Neighbors | グラフ・グリッド第二世代 00_grid_base | 体系外（5-4） |
| 275C.txt Manhattan / Chebyshev / Circle の範囲 in() | in_diamond / chebyshev(c, p) ≤ D / contains(c, p) | ◯ |
| ManhattanTemplateTest.txt add_point / max_dist / max_pair / farthest_dist / farthest_point / diamond | ManhattanFarthest / DynamicManhattan（max_pair・farthest）/ in_diamond | ◯ |
| ManhattanMSTemplateTest.txt（＋ erase_point） | DynamicManhattan::erase | ◯ |
| RangeSplitTemplateTest.txt Range（1 次元）/ sign_split / cut_by_range | 区間の道具 | 体系外（5-4） |
| RangeSplitTemplateTest.txt Rect（閉区間。Rect(x1, x2, y1, y2) の順） | Rect（半開）/ rect_from_closed / **rect_from_ranges(x1, x2, y1, y2)**（同じ引数の順） | ◯ |
| RangeSplitTemplateTest.txt quadrant_split / sign_grid_split（番号 3i + j）/ split_by_x_axis / split_by_y_axis（原点だけ） | 同じ名前・同じ番号（中心を指定できる） | ◯ |
| RectangleDiagonalQueryTemplate.txt RectangleDiagonalQuery(L, R, D, U).plus(c) / minus(c) | diagonal_lattice(rect_from_closed(L, D, R, U), c, ±1) | ◯ |

### 7-3. 過去の解答・Dynalist の項目との照合

| 過去の解答 | 新しいテンプレートでの書き方 | 判定 |
|---|---|---|
| 266C.txt（凸四角形：自作の V と ccw） | ccw / is_convex（典型の照合で解いた） | ◯ |
| 374D.txt（順列と bit 全探索で線分をたどる距離） | dist（幾何の部分は距離だけ） | ◯ |
| 259D.txt / 259D_2.txt（円のつながり） | has_common_point ＋ Union-Find（典型の照合で解いた） | ◯ |
| 275C.txt（正方形の数） | count_squares / rot90（典型の照合で解いた） | ◯ |

| Dynalist の項目 | 新しい置き場所 |
|---|---|
| 09_幾何「AtCoder 幾何フルテンプレ（整数主導）」 | 00〜05・09 に分けた（7-2） |
| 09_幾何「円と円の関係」 | 05_circle circle_relation |
| 09_幾何「90度回転（任意中心・左右統合）」 | 00_geo_base rot90(p, c, k) / Affine |
| 09_幾何「範囲判定（動的サイズ対応）」 | グリッドの範囲はグラフ・グリッド 00_grid_base。幾何の範囲は in_diamond / contains / Rect |
| 09_幾何「マンハッタン距離 座標変換テンプレ（復元付き）」 | 07_manhattan manhattan_to_chebyshev / chebyshev_to_manhattan / ManhattanFarthest |
| 09_幾何「マンハッタン距離テンプレ（multiset対応・動的版）」 | 07_manhattan DynamicManhattan |
| 13・15「格子点・整数区間の分割をまとめたテンプレート」 | 08_rect quadrant_split など（1 次元の分割は体系外） |
| 13・15「矩形を斜めに切り取ったときの格子点」 | 08_rect diagonal_lattice |
| 07_データ構造「2D BIT」 | 長方形の中の点の数はこれを使う（幾何の索引から指す） |

旧ファイルの機能・過去の解答・Dynalist の項目で、新しいテンプレートにないものはない（体系外としたものを除く）。

### 7-4. 典型の問題（問題文を確かめて解き、サンプルを通した）

**Phase 2〜4 の照合の 33 問**（ABC033 D, ABC139 F, ABC144 D, ABC151 F, ABC157 F, ABC168 C, ABC181 C, ABC191 D, ABC197 D, ABC198 C, ABC207 D, ABC218 D, ABC224 C, ABC225 E, ABC226 D, ABC248 E, ABC250 F, ABC259 B, ABC259 D, ABC266 C, ABC275 C, ABC296 G, ABC304 C, ABC318 B, ABC343 E, ABC351 E, ABC362 B, ABC374 D, ARC042 B, 典型90 009 / 018 / 036 / 041）：問題文を開いて確かめ、テンプレートだけで解いてサンプルをすべて通した。最終照合で追加・修正した後のファイルで組み直して、**33 問ともサンプルを通ることを確かめ直した**（ABC343 E は答えが複数あるので専用の判定で確認）。

このとき見つかった使いにくさと直したこと：

| 見つかったこと | 直したこと |
|---|---|
| `rotated(Point, θ, PointD)` がコンパイルできない（ABC197 D。デモの例のとおりに書けない） | 中心の型を自由にした（2 次元・3 次元とも）。結果の −0 を 0 に |
| ABC157 F のデモの例 `Circle(p, T / c)` がコンパイルできない（整数の円に実数の半径） | `CircleD(to_d(p), T / c)` に直し、実数の判定の誤差で答えが 1e-9 程度ずれることを書いた |
| 入力で与えられた凸多角形を convex_contains にそのまま渡すと、h[0] が 180 度の頂点のとき誤る | 「convex_hull(P, true) を通す」をヘッダ・デモ・索引に書いた |
| 08_rect の rect_from_closed を連続な長方形（ABC318 B）に使うと 1 ずれる。入力 A B C D（x の範囲・y の範囲）を Rect(x1, y1, x2, y2) に渡し間違えやすい | 格子点を数えるときだけ rect_from_closed と書き直し、**rect_from_ranges(x1, x2, y1, y2)** を足した |
| 索引：ABC144 D（tan θ = h / w）の行がない、ABC139 F は専用の関数がないのに △ がない、ABC351 E は偶奇で分けることが書いていない、デモの ☆ 印の意味が書いていない | 索引の行を足し・直した。デモの ☆ 印を消した |

**最終照合で足した機能を使って新しく解いた 15 問**（問題文を開いて確かめ、テンプレートだけで解いた。サンプル・愚直解とのランダム比較・最大の入力での時間を確認）

| 問題 | 使った機能 | サンプル | ランダム比較 | 最大の入力（この環境） |
|---|---|---|---|---|
| ABC338 E | count_crossing_chords | 3/3 | 1000 件一致 | 0.26〜0.36 s |
| ABC405 F | chord_cross_counts | 2/2 | 1000 件一致 | 0.26 s |
| ABC424 F | ChordSet | 2/2 | 1000 件一致 | 0.14 s |
| 典型90 017 | count_crossing_chords | 5/5 | 1000 件一致 | 0.37 s |
| ABC402 D | chord_direction_key | 2/2 | 1000 件一致 | 0.32 s |
| ARC076 E | boundary_position（08_rect）＋ count_crossing_chords（05_circle。2 本を一緒に貼ってよいことも確認） | 4/4 | 1000 件一致 | 0.21 s |
| ARC051 A | contains(poly, c) / disk_contains(c, poly) / to_polygon | 5/5 | 1000 件＋全探索 59 万件一致 | － |
| ABC189 E | Affine（then の累積） | 2/2 | 400 件一致 | 0.16 s |
| ABC366 E | abs_sum_table | 3/3 | 300 件一致 | 1.0 s（大半は自分で書く sort） |
| ABC283 F | manhattan_nearest_each | 3/3 | 400 件一致 | 0.57 s |
| ABC250 F | PolygonPrefix | 3/3 | 420 件一致（頂点 10^5 でも全探索と一致） | 0.03 s |
| ABC016 D | intersects_strict | 2/2 | 300 件一致 | － |
| ABC197 D | rotated(Point, θ, PointD) | 2/2 | 300 件一致 | － |
| ARC029 B | rect_fits | 1/1 | （rect_fits 自体を角度の全探索とテストで比較） | － |
| ABC301 G | line_key（LineKey3）/ cross_point_frac | 2/2 | 有理数の愚直解と 150 件一致 | 0.02 s |

このときの指摘で直したこと：`ChordSet(n)` が位置 n を受け付けず 1〜N の番号で範囲外に書き込むことがあった → 位置 0〜n を受け付けるようにした。`boundary_position` に rect_from_closed で作った長方形を渡すと誤る → 「Rect(0, 0, R, C) をそのまま」と書いた。`manhattan_nearest_each` が番号だけで距離を自分で計算する必要があった（07_manhattan_mst には manhattan がない）→ `ManhattanNear{j, d}` を返すようにした。`reflect_x(p)` は「x = p で折り返す」で、x 軸での折り返しと取り違えやすい → デモに「x 軸で折り返すなら reflect_y(0)」と書いた。

### 7-5. 足したもの・直したもの（ファイル別）

| ファイル | 足したもの・直したもの |
|---|---|
| 00_geo_base | rotated / rotated_deg の中心の型を自由に・−0 を 0 に、AffineT（Affine / AffineD） |
| 02_line | segment_relation、半直線 RayT（on_ray / intersects / cross_point / cross_points / dist / dist2_cmp） |
| 03_polygon | 02_line を同梱。contains(poly, s) / contains(P, Q) / intersects(poly, s or l or r) / intersects(P, Q) / dist / boundary_dist、PolygonPrefix |
| 04_convex_hull | convex_contains(h, Q)、h[0] の注意。02_line を同梱（03_polygon が使うため） |
| 04_convex_polygon | convex_intersects、入力の凸多角形の注意 |
| 05_circle | in_circumcircle、円と半直線、disk_intersects / disk_contains / disk_dist / circumference_dist、弦（chords_cross / chord_direction_key / chord_cross_counts / count_crossing_chords / ChordSet）、デモの ABC157 F・☆ 印 |
| 05_circle_area | contains(poly, c) / disk_contains(c, poly) / disk_intersects(c, poly) / disk_dist(c, poly) |
| 07_manhattan | AbsSum / abs_sum_table |
| 07_manhattan_mst | manhattan_nearest_each（ManhattanNear） |
| 08_rect | rect_from_ranges / where_closed / intersection_dim（Rect・BoxK）/ boundary_position / to_polygon / rect_fits、rect_from_closed の説明 |
| 09_halfplane | 02_line を同梱（03_polygon が使うため） |
| 10_geo3d | line_relation / cross_point / cross_points / FracPoint3・cross_point_frac / direction_key / LineKey3・line_key / line_plane_relation / intersects・dist（直線・線分と平面）/ plane_key / plane_relation / dist（平面どうし）/ sphere_plane_relation / Circle3D・cross_circle / ball_intersects / ball_contains / ball_dist、rotated の型を自由に |
| 01_arg_sort・05_circle_cover・06_point_set・06_congruence | 変更なし（同梱している基盤・02_line の更新だけ） |

関数を足したファイルには、すべてランダムテスト（愚直な実装との比較。有理数・全探索・十分に長い線分での置き換え・角度の全探索など）と、使い方デモ（日本語の説明と類題）を足した。

### 7-6. 広い洗い出し（AtCoder の幾何 109 問）から決めたこと

Phase 4 の後で、AtCoder の幾何の問題 109 問（ABC 200 以降 65、それ以前 17、ARC 21、典型 90 6）を、すでに確かめた 34 問を除いて洗い出した。テンプレートでそのまま解ける ◯ が 28、数行足せば書ける △ が 40、足りない × が 32、幾何でない － が 9 だった。

| 足りなかったもの | 決めたこと |
|---|---|
| 円周の上の番号で与えられる弦（× の 8 問。最近の ABC の D〜G に 5 回） | 05_circle に弦の関数を足した（7-4 で 6 問を解いた） |
| 図形どうしの位置関係・包含・接触（長方形が辺で接する、円と長方形の包含、長方形が回して入るか、3 次元の直線の交点とキー、線分で多角形を切る） | 7-1 の表のとおり足した |
| 凸多角形の頂点の区間の面積・重心 | 03_polygon に PolygonPrefix |
| 1 次元の Σ w\|x − x_i\| の表 | 07_manhattan に AbsSum / abs_sum_table |
| マンハッタン距離で各点の最近点 | 07_manhattan_mst に manhattan_nearest_each |
| 整数のアフィン変換の合成 | 00_geo_base に Affine |
| 点の集合の範囲カウント（長方形の中の点の数） | 入れない（データ構造。07_データ構造 の 2D BIT を索引から指す） |
| 凸包の逐次・動的クエリ・Li Chao 木 | 入れない（データ構造。G / Ex 級） |
| 境界が有理数の領域の格子点（floor_sum） | 入れない（数学の道具） |
| 長方形の重なりの枚数ごとの面積（典型90 028） | 入れない（座標が小さいので 2 次元いもす法。グリッドの道具） |
| ユークリッド距離の最近点クエリ（ARC010 D）・障害物を避ける最短路・回転走査など | 入れない（出題が少なく、専用の構造が大きい） |
| 関数はあるのに索引から見つけにくいもの（halfplane の contains(h, p)、r.contains(r2)・bounding、dot128 系、45 度回転の使いどころ、時計回りの範囲、bisector_key・union_area・intersects_strict の使いどころ、convex_diameter・triangle_area2・arc_length など） | 索引に行を足した（7-7） |

### 7-7. 索引・Dynalist の項目の更新

- 索引（`dynalist_tool_index.txt` / `.opml`）：「▼ 図形の位置関係（早見）」（20 行）を「使い方」の次に新しく置いた。各大分類に、足した関数と見つけにくかった関数の行を足した（142 行 → 213 行）
- Dynalist の項目：`dynalist_items_phase2.txt`（00_geo_base）、`dynalist_items_phase3.txt`（02_line・03_polygon・04_convex_hull・05_circle・07_manhattan）、`dynalist_items_phase4.txt`（04_convex_polygon・05_circle_area・07_manhattan_mst・08_rect・09_halfplane・10_geo3d）の該当する項目を書き直した。Dynalist には、この 12 項目を入れ替える

### 7-8. 確認したこと

- 16 本すべて：g++（-Wall -Wextra -Wshadow -pedantic で警告 0）・clang（-pedantic-errors で警告 0）・g++ の C++20 / C++23（警告 0）・AddressSanitizer ＋ UndefinedBehaviorSanitizer ＋ _GLIBCXX_DEBUG・浮動小数点のオーバーフロー検出・Real を double にした版で、ランダムテストがすべて通った
- 貼り付けの確認：Main_base ＋ ACL ＋ グラフ・グリッド 46 本と一緒に貼ってもぶつからない（2 通りの順・gnu++2b / c++17・LOCAL あり / なし。足した関数も呼んで確かめた）、各ファイルを単独で Main_base に貼って動く（-Wshadow の警告 0。弦の関数の中の INF という名前が Main_base の INF と重なっていたので直した）、同梱しているブロック 48 個がすべて持ち主のファイルと同じ
- 参考：グラフ・グリッド第二世代の 2 か所（木の最小支配集合・最小費用流の slope）で、ローカル変数 INF が Main_base の INF を隠す -Wshadow の警告が出る（動作には影響しない。幾何の範囲外なので直していない）
- 典型 33 問を組み直して確認、新しい機能で 15 問（7-4）
- PC（`Desktop\claude\geometry_v2`、g++ 11.4）に保存して、16 本ともビルド（警告 0）・ランダムテストが通ることを確かめた
- 独立したレビューは、この回は使える回数の上限で別のエージェントを動かせなかったため、自分で辺の場合（長さ 0 の線分、始点が相手の上にある半直線、一直線の重なり、半径 0 の円・球、頂点 1〜2 個の多角形、番号が負・n 以上の PolygonPrefix など）を 1 つずつ見直し、実数の型でも整数と同じ判定になることをテストに足した。問題を解いた確認（7-4）では、2 つのエージェントがそれぞれ別の 7 問・6 問を解いて使いにくさを報告した

### 7-9. 数学第二世代との共通化（2026-10-10）

数学第二世代（`Desktop\claude\math_v2\design.md`）の Phase 1 の確認で、数学と幾何で同じ働きの関数を 1 つの名前にそろえることにした（数学の設計書 6-1 の確認 2：案 A）。

- `00_geo_base` の本体にあった `Int128` の typedef・`isqrt` / `isqrt_ceil`・`parse_fixed` / `read_fixed`・`int128_to_string` を、共通の基盤 `00_int_common.cpp`（ガード `GG_INT_COMMON`）に移した。16 本すべてが `GG_GEO_BASE` の前にこのブロックを同梱する。geometry_v2 にも math_v2 と同じ `00_int_common.cpp` を置いた（17 本。update_bundles がそのまま使える）
- 関数の名前・動きは同じ。変わったのは 2 つだけ：Int128 を `cout <<` / `cin >>` でそのまま出せる・読めるようになった（`uint128_to_string` / `parse_int128` なども足した）、`parse_fixed` が long long に収まらない値で止まるようになった（以前は黙って溢れていた）
- 見出しの「依存」に 00_int_common を足し、00_geo_base の説明・デモの「cout で直接出せない」を直した
- 確かめたこと：17 本とも g++（-O2 -Wall -Wextra -Wshadow -pedantic）・clang（-pedantic-errors）・C++20 / C++23 で警告なし・全テスト通過、ランダムテストの照合回数と使い方デモの出力が更新前と同じ（変えた説明の 2 行を除く）、`-D_GLIBCXX_DEBUG -fsanitize=address,undefined` で全体が通過、Real を double にした版も通過、Main_base（ACL 込み）にグラフ・グリッド 46 本・幾何 17 本・数学の本体を 2 通りの順で全部貼って g++・clang、gnu++2b / c++17、LOCAL あり / なしでコンパイル・実行、同梱部分 48 個が持ち主と一致
- **00_int_common を直すときは、math_v2 と geometry_v2 の両方の 00_int_common.cpp を同じにしてから**、それぞれのフォルダで update_bundles を実行する
