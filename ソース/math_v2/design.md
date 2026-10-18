# 数学テンプレート（幾何を除く） 第二世代 設計書

- 作成日：2026-10-10
- 状態：**Phase 1（体系設計）完了**（2026-10-10 ユーザー確認済み。6-1）。**Phase 2（基盤 00_int_common・00_math_base・00_modint）完了**（2026-10-10。2-P2 章）。**Phase 3（★ 10 本）完了**（2026-10-10。2-P3 章）。**Phase 4（☆ 13 本）完了**（2026-10-11。2-P4 章）。**最終照合完了**（2026-10-11。7 章：旧テンプレート・過去の解答・Dynalist・典型 48 問・AtCoder の数学 349 問の洗い出し）
- 方法：グラフ・グリッド第二世代（`Desktop\claude\graph_grid_v2\design.md`）・幾何第二世代（`Desktop\claude\geometry_v2\design.md`）と同じ方法・観点で作り直す
- 次回はこのファイルの「6. 未決事項と次にやること」と「7. 最終照合」から再開する

---

## 0. 目的・観点・これまでから引き継ぐこと

### 0-1. 目的

幾何を除いた数学（整数の基本・mod・桁と進数・素数と約数・合同式・組合せ・行列・数列の和・畳み込み・分数・ビット・ゲーム）のテンプレートを「第二世代」としてゼロから設計し直す。既存テンプレートの整理ではなく、新しい体系を作る。

成果物はテンプレートそのものと、**問題の特徴から使う道具に辿れる選択体系**（Dynalist の索引）。例えば次のように辿れる状態を目指す。

- 答えを 998244353 で割った余り → 00_modint（`using mint = ModInt998;` の 1 行）
- 「n 個から r 個選ぶ」が何度も出る → 04_binomial（`Binomial<mint> B(N); B.C(n, r)`）
- N 以下の数をすべて素因数分解する → 02_sieve（`Sieve S(N)`）。10^18 の数を 1 つなら 02_factorize（`factorize(n)`）
- Σ ⌊N/i⌋ を N = 10^12 で → 06_quotient（`quotient_blocks(N)`）
- x ≡ r (mod m) がいくつも与えられる → 03_mod_equation（`crt_all`）
- 漸化式の第 10^18 項 → 項が少なければ 05_matrix（`pow`）、多ければ 06_linear_recurrence
- a^b が 10^18 を超えるかを比べたい → 00_math_base（`pow_cap`）
- 負の数を割って切り捨てる → 00_math_base（`floor_div`）

### 0-2. 観点（グラフ・グリッド・幾何と同じ 3 つ）

| 観点 | 数学での具体化 |
|---|---|
| アクセスしやすさ | 「問題の特徴（答えを mod で・大きな数・素数・約数・合同式・数え上げ・漸化式・和・ビット…）→ ファイルと関数」で辿れる索引。ファイル名は分類番号＋名前。依存するファイルは同梱し、1 回のコピーで動く |
| 使いやすさ | modint・Int128・分数・多倍長を int / long long と同じ感覚で使える（`val()` なしで cout、負の数から作れる）。負の数・0・溢れで黙って間違えない（`floor_div`、`_cap` 版、溢れたら止める）。呼び方・戻り値の形をそろえる（解がないことがある関数は `ok` を持つ構造体、素因数分解は `{p, e}` の列） |
| 連携 | modint の型を基盤で統一し、組合せ・行列・畳み込み・漸化式がすべて同じ型を受け取る（`Binomial<mint>`、`Matrix<mint>`、`convolve(a, b)`）。グラフ・グリッド第二世代（46 本）・幾何第二世代（16 本）と同時に貼れる。グラフと組み合わせる形（行列木定理に辺のリスト、行列の累乗に隣接行列）をそろえ、組み合わせ方を索引に載せる |

### 0-3. 最重要の方針（グラフ・グリッドの 1-0 を最初から適用する）

第二世代のテンプレートは互いに関係し、旧テンプレートとの互換性もないため、**後から追加・修正しなくて済むこと**が最も重要。

- 追加でほしくなる可能性のある機能はすべて入れる。入れるか迷った機能は入れる
- 速さの改善が要らないくらい速くする
- 正しいこと（愚直解とのランダムテスト、サニタイザ、Main_base＋ACL に貼った状態でのコンパイル確認）
- 使いやすいこと

### 0-4. グラフ・グリッドと幾何での経験から、最初から取り入れること

| これまでに起きたこと | 数学での対応 |
|---|---|
| 1-0 の方針が Phase 2 の途中で決まった（グラフ） | 最初から適用する（0-3） |
| 最終照合で足りないものが見つかった（グラフは 7 本。幾何は 109 問の洗い出しで「足りない」が 32 問） | **典型との照合を Phase 1 で行い**（5-5）、ABC 300 以降の最近の問題も含めた。最終照合でも数学の問題を 100 問程度、広く洗い出す |
| 「数行で書けるが間違えやすい」ものをテンプレートにして役立った（ワーシャルフロイド、偏角の比較、小数の入力） | 最初から入れる：負の数の切り捨て・切り上げ・余り、溢れない累乗・lcm、誤差のない平方根、逆元を使わない等比数列の和、範囲の中の倍数の数、多くの数の lcm を mod で（mod を取りながら lcm を取ると誤る） |
| 使い方関数の例のとおりに書くとコンパイルできなかった（幾何の ABC197 D・ABC157 F。最終照合で見つかった） | 使い方関数の例はすべて実際にコンパイル・実行されるコードにする。**各 Phase の終わりに代表的な問題をテンプレートだけで解いて確かめる**（最終照合まで待たない） |
| std の lower_bound などがデバッグビルドで O(n) の検査をして終わらなくなった | 本体の二分探索・ヒープは自分で書く |
| MSVC ではビルドしない・`__int128` は使える（2026-10-08 確認） | `__int128` をそのまま使う |
| 名前の衝突を Main_base＋ACL＋グラフ・グリッドで確認した | **Main_base＋ACL＋グラフ・グリッド 46 本＋幾何 16 本**と一緒に貼って確認する。数学は ACL と名前が重なりやすい（modint・convolution・crt・floor_sum・pow_mod・inv_mod）。Main_base は `using namespace atcoder;` なので、同じ名前・同じ引数だと呼び出しがあいまいになりコンパイルエラーになる |
| 作ったものを見ていない別のエージェントのレビューで不具合が 4 件見つかった（幾何 Phase 4） | Phase 3・4 の終わりに毎回行う |
| 各 Phase で Dynalist の項目と索引を作った | 同じ（dynalist_items_phaseN.txt、dynalist_tool_index.txt / .opml） |

---

## 1. 決定事項

### 1-1. 進め方・既存テンプレートの扱い

- 既存テンプレートは「過去の実験結果」として扱う。既存の API を正解として引き継がない
- 参照したのは、新しい一覧と機能を決めた後（Phase 1 の終盤、5 章）だけ
- 原本は変更・削除しない。新しく作るものは `Desktop\claude\math_v2\` に出力する（最後に `C:\atcoder_template\典型90\ソース\math_v2\` へ置く想定。グラフ・グリッド・幾何と同じ）
- 既存の場所：`C:\atcoder_template\典型90\ソース`（5-0）
- Phase 1（体系設計）→ Phase 2（基盤）→ Phase 3（★）→ Phase 4（☆）→ 最終照合 の順に進め、区切りごとにユーザーの確認を取る
- **ACL は使わない**。modint・畳み込み・floor_sum・crt も std だけで自作する。速さは ACL と同じかそれ以上を目標にして、計測で比べる

### 1-2. 共通規約

グラフ・グリッド・幾何と同じもの：

| 項目 | 決定 |
|---|---|
| 優先順位 | 速度 → int / long long と変わらない使い勝手 → 正しさ |
| ファイル | 1テンプレート1ファイル。単体で実行できる .cpp |
| テスト | ランダムテスト関数と、全機能を使う使い方関数を同じファイルに入れ、main から呼ぶ（形は同じ：`xxx_test::random_test()`（Release は 8 倍の量）・`benchmark()`・`xxx_demo::demo()`・`setup_console()`） |
| 使い方の出力 | 関数名の表示に加え、何を返す・提供するかを日本語でも書く |
| 使うライブラリ | std のみ（ACL 不使用）。bits/stdc++.h に頼らず、使うヘッダを書く |
| 名前 | 型は PascalCase（`Sieve`, `Binomial`）、関数・メンバは snake_case（`floor_div`, `quotient_blocks`） |
| Main_base への依存 | テンプレートの中では Main_base の別名・マクロを使わない |
| 再帰 | 非再帰を基本にする。深さが log N で済むもの（分割統治の多項式の積、拡張ユークリッドなど）だけ再帰でよい |
| 「なし」 | 番号は -1。解がないことがある計算は `ok` を持つ構造体（下） |

数学で新しく決めるもの（★ は特に確認してほしいもの。6-1）：

| 項目 | 決定（案） |
|---|---|
| 整数の型 | 引数・戻り値は long long が基本。途中で溢れるところは 128 ビット（`Int128` = `__int128`）で計算する。答えが long long に収まらないことがあるものは、Int128 で返す版か、上限で止める版（名前に `_cap`）を用意する |
| 溢れ・範囲の保証 | 各関数は**引数の全範囲で正しい**か、正しい範囲をヘッダの一覧に書く（例：篩は N ≤ 5×10^7、ポラード・ローは n < 2^63）。**黙って溢れない**：厳密な値を返す関数（`pow_int`、分数の四則など）は溢れを検出したらメッセージを出して止める |
| 割り算・余り | 負の数があるところでは C++ の `/`・`%`（0 に向かって切り捨て）を使わず、`floor_div` / `ceil_div` / `mod_floor`（0 以上の余り）を使う。テンプレートの中も同じ |
| ★ modint | 定数の mod は `ModInt<MOD>`（別名 `ModInt998`・`ModInt107`）、入力で決まる mod は `DynModInt`（`DynModInt::set_mod(M)`）。mod は 1 以上 2^31 未満（素数でなくてよい、偶数でもよい）。逆元は拡張ユークリッドで求めるので、mod が合成数でも gcd = 1 なら求まる（なければメッセージを出して止める。止めずに調べるなら `try_inv`）。`val()` を書かずに cout できる。**`mint` という名前はテンプレートでは定義しない**（使う側で `using mint = ModInt998;` と 1 行書く）。数え上げ・行列・畳み込みのテンプレートは型を引数に取る（`Binomial<mint>`・`Matrix<mint>`）ので、ModInt と DynModInt のどちらでも動く |
| 64 ビットの mod | mod を引数で渡す関数（`mod_pow(a, e, m)`・`mod_inv(a, m)`・素数判定の内部）は 1 ≤ m < 2^63 で正しい（積は 128 ビット）。a は負でもよい |
| ★ 範囲の引数 | グラフ・幾何と同じく半開区間 [l, r)。問題で「L 以上 R 以下」と与えられることが多いもの（区間篩、倍数の数など）は閉区間版を `_closed` という名前で用意する。「N 以下」を表す関数は名前に `_up_to` を付ける（N を含む） |
| 解がないとき | 解がないことがある関数は `ok` を持つ構造体で返す（-1 などの番兵で返さない。-1 も答えになりうるため） |
| 素因数分解の形 | `vector<PrimePower>`（`PrimePower{p, e}`、p の昇順）。篩・単発の両方で同じ形 |
| 約数 | 昇順の `vector<long long>` |
| 表 | 「N まで」の表は長さ N + 1（添字 0〜N をそのまま使える） |
| 実数 | 実数の答えは long double。出し方（`fixed << setprecision(15)`）は使い方関数で示す |
| 名前の衝突 | Main_base（マクロ `rep` `erep` `rrep` `nall` `rall` `YES` `NO` `YN` `BUFF` `FIXED` `next_p` `debug`、名前 `ll` `ull` `sgn` `popcnt` `INF` `LINF` `PQ` `PQ_g` `chmin` `chmax` `sizei` `GraphIO`）・**ACL**（`using namespace atcoder;` で見える `modint` `static_modint` `dynamic_modint` `modint998244353` `modint1000000007` `convolution` `convolution_ll` `crt` `floor_sum` `pow_mod` `inv_mod` `internal` など）・std（`gcd` `lcm` `pow` `sqrt` `cbrt` `abs` `div` `exp` `log` `floor` `ceil` `round` `popcount` `bit_width` `bit_ceil` `bit_floor` `has_single_bit` `countl_zero` `countr_zero` `rotl` `rotr` `midpoint` `lerp` `beta` `iota` `accumulate` `reduce` `partial_sum` `transform`）・**グラフ・グリッド 46 本と幾何 16 本のグローバルな名前（508 個。`transpose` `doubling` `perm_compose` `congruence` `dot` `cross` `contains` など。2 章の案の名前と照らし合わせ、重なるのは幾何と同じにしたい `Int128` `isqrt` `isqrt_ceil` `int128_to_string` `parse_fixed` `read_fixed` だけ（6-1））**・よく使う変数名（N, M, K, Q, MOD, mod, P, A, B, C, D, X, Y, MAX, fact, inv, finv, primes, divs, dp, cnt, ans, res）・Windows のマクロ（`min` `max` `IN` `OUT`）と同じ名前をグローバルに置かない。配列によく使う名前（fact・primes など）は関数名にしない（`primes_up_to`・`divisors_of` のように前置詞などを入れる） |
| 出力 | テンプレートは cout しない（止めるときの cerr のメッセージを除く）。modint・Int128・分数・多倍長は `cout <<` で答えの形のまま出せる（幾何の Point の `(x, y)` のようなデバッグ用の形ではない） |

### 1-3. 共通API

```cpp
// ---- 使う側で 1 行（ModInt107・ModInt<任意の定数>・DynModInt も同じ使い方）----
using mint = ModInt998;
mint a = -3, b = 10;                  // 負の数・long long・Int128 から作れる
cout << a * b / 7 << '\n';            // val() なしで出せる
a.pow(K); a.inv(); mint::mod();
DynModInt::set_mod(M);                // 入力で mod が決まるとき（素数・奇数でなくてよい）

// ---- 整数の基本（00_math_base）----
floor_div(a, b); ceil_div(a, b); mod_floor(a, m);   // 負の数でも数学の切り捨て・切り上げ・余り
isqrt(n); kth_root_floor(n, k);                       // 誤差のない ⌊√n⌋・⌊n^(1/k)⌋
pow_cap(a, e, LIM);                                   // a^e（LIM を超えたら LIM。溢れない）
mod_pow(a, e, m); mod_inv(a, m);                      // m < 2^63・合成数の mod でも
Int128 x = ...; cout << x;                            // Int128 も cout で出せる

// ---- 前処理してから何度も聞くものは構造体 ----
Sieve S(10'000'000);   S.factorize(x);  S.divisors(x);  S.phi(x);
Binomial<mint> B(N);   B.C(n, r);  B.P(n, r);  B.H(n, r);
Matrix<mint> A(K, K);  auto Ae = A.pow(e);
XorBasis xb;           xb.insert(v);  xb.max_xor();

// ---- 1 回で終わる計算は関数 ----
auto fs = factorize(n);               // {{p, e}, ...}（n < 2^63。ポラード・ロー）
auto ds = divisors_of(n);             // 昇順
for (auto [q, l, r] : quotient_blocks(N)) { ... }   // ⌊N/i⌋ = q となる i ∈ [l, r)
auto c = crt_all(rs, ms);             // c.ok, c.r, c.m（mod は互いに素でなくてよい）
auto h = convolve(f, g);              // mint の列の積（mod は何でも）
```

- 1 回で終わる計算は関数（`factorize(n)`、`quotient_blocks(N)`）、前処理してから何度も聞くものは構造体（`Sieve S(N)`、`Binomial<mint> B(N)`）
- 複数の値を返す関数は、名前付きのメンバを持つ構造体で返す（構造化束縛でも受け取れる。pair・tuple で返さない）
- modint を受け取るテンプレートは ModInt・DynModInt のどちらでも動く。mod が素数でないと意味がない関数（Binomial の逆元など）は、ヘッダの一覧に書く

### 1-4. ファイルの形と共通化

- グラフ・グリッド・幾何と同じ：**.cpp 1 ファイル**、コピー範囲の目印（`ここからテンプレート本体（コピー範囲）` / `ここまでテンプレート本体`）、依存するテンプレートは**貼る範囲に同梱**して `#ifndef` ガードで重複を防ぐ
- ガード名は同じ付け方（`GG_` ＋ファイル名の大文字。例：`GG_MATH_BASE`、`GG_MODINT`）。`update_bundles.ps1 / .bat` を**そのまま**使える（math_v2 にコピーするだけ）
- ファイル名は「分類番号 2 桁＋名前」。**名前の部分はグラフ・グリッド・幾何のファイルと重ならない**ようにする（幾何に `06_congruence` があるので、合同式のファイルは `03_mod_equation` にした）
- 内部の名前空間は `math_internal`（`gg_internal`・`geo_internal` とは分ける）
- **すべてのファイルが 00_math_base を同梱する**（幾何の geo_base と同じ）
- 文字コードは UTF-8（BOM 付き）、main で Windows のコンソールの設定をする
- 幾何と共通の関数（`isqrt` / `isqrt_ceil`・`Int128` / `UInt128`・`int128_to_string`・`parse_fixed` / `read_fixed`）は、数学でも同じものが要る。幾何の基盤はこれらをガードなしで定義しているので、同じ名前のまま数学にも入れると、2 つを一緒に貼ったときに二重定義になる。扱いは 6-1 の確認 2 で決める
  - **案 A に決定（2026-10-10）**：この部分を小さな基盤 `00_int_common.cpp`（ガード `GG_INT_COMMON`）に切り出し、`00_math_base` と `00_geo_base` の両方が同梱する。geometry_v2 にも同じ `00_int_common.cpp` を置けば update_bundles でそのまま更新できる。幾何の 16 本は同梱部分が変わるだけ（関数の名前・動きは同じ）で、テストをやり直す。Int128 の `<<` / `>>` もここに入れる
  - （採らなかった案 B：数学の基盤に別の名前（例：`sqrt_floor`・`i128_to_string`）で入れる。幾何のファイルは変えないが、同じ働きの関数が 2 つの名前になる）

### 1-5. 正しさ・速さの確認のしかた（数学で追加するもの）

- **境界の値を必ず入れる**：0・±1・LLONG_MAX・LLONG_MIN（受け付ける関数）・2^31 ± 1・2^32・2^62・k² − 1 / k² / k² + 1（平方根）・k^3 の前後（立方根）・素数と合成数の境目・カーマイケル数と強擬素数（ミラー・ラビン）・mod = 1・2・合成数・偶数
- **愚直解**：小さい範囲の全探索、128 ビットや自作の小さな多倍長での計算、別のアルゴリズム（試し割り ⇔ ポラード・ロー、愚直な畳み込み ⇔ NTT、篩 ⇔ 素数の個数、行列の累乗 ⇔ 漸化式を回す、Berlekamp–Massey で元の漸化式に戻るか）
- **modint を受け取るテンプレート**は、ModInt と DynModInt（mod が素数・合成数）の両方でテストする
- **速さは ACL と比べる**：modint の積・畳み込み（長さ 2^20）・floor_sum・crt。ほかの目安：篩 N = 5×10^7、10^18 の数の素因数分解 × 10^4、Binomial の表 N = 10^7、100×100 の行列の 10^18 乗、quotient_blocks（N = 10^14）、prime_count（N = 10^11）
- **名前の衝突**：Main_base＋ACL に、グラフ・グリッド 46 本・幾何 16 本・数学の全ファイルの本体を 2 通りの順で全部貼り、gnu++2b / c++17・LOCAL あり / なしでコンパイル・実行する。単独のテンプレート（binary_search.cpp・static_seq_toolkit.cpp など）とも一緒に貼って確かめる
- **独立したレビュー**：Phase 3・4 の終わりに、作ったものを見ていない別のエージェントに本体を読ませ、全探索・厳密計算で確かめさせる
- **問題で確かめる**：Phase 3・4 の終わりに、5-5 の代表例の問題文を確かめてテンプレートだけで解き、サンプルを通す（幾何では最終照合まで待ったため、使いにくさの発見が遅れた）

---

## 2. 分類とテンプレート一覧（確定。2026-10-10）

記号：★＝優先して整備（Phase 2・3）、☆＝余裕が出てから（Phase 4。ただしグラフ・グリッド・幾何と同じく全部作る）、△＝専用テンプレートは作らず「型」として持つ

関数名は案。Phase 2〜4 の実装で確定し、この表に反映する（グラフ・グリッド・幾何と同じ）。

### 0. 基盤（Phase 2）

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `00_int_common` | **作成済み（Phase 2）**。数学・幾何で共通の整数の部分：`Int128` / `UInt128`（= __int128）、Int128・UInt128 の `<<` / `>>`、`int128_to_string` / `uint128_to_string` / `parse_int128`、`isqrt` / `isqrt_ceil`、`parse_fixed` / `read_fixed`（小数を 10^k 倍の整数に）。00_geo_base（幾何）と 00_math_base（数学）の両方が同梱する | なし | ABC169 C |
| ★ `00_math_base` | **作成済み（Phase 2）**。整数の基本（Int128・isqrt・parse_fixed は 00_int_common）。`Maybe<T>{ok, value}`。**割り算**：`floor_div` / `ceil_div` / `mod_floor`（0 以上の余り）/ `round_div`（四捨五入）、**追加：`ceil_multiple` / `floor_multiple` / `ceil_congruent` / `floor_congruent`（x 以上・以下で最も近い倍数・合同な数）**、`count_multiples(l, r, d)`（[l, r) の d の倍数の数。負でも）/ `count_multiples_closed` / `count_congruent(l, r, a, m)`（x ≡ a (mod m) となる x の数）。**累乗根**：`kth_root_floor(n, k)` / `kth_root_ceil` / `is_square` / `is_kth_power`。**累乗・溢れ**：`pow_int(a, e)`（溢れたら止める）/ `pow_cap(a, e, cap)` / `mul_cap` / `add_cap` / `mul_overflows(a, b)` / **`add_overflows`** / `ilog_floor(n, b)`（⌊log_b n⌋）/ **`ilog_ceil`** / `pow10_ll(k)`。**gcd・lcm**：`ext_gcd(a, b)`（`ExtGcd{g, x, y}`）/ `lcm_cap(a, b, cap)` / `gcd_of(v)` / `lcm_of_cap(v, cap)`。**64 ビットの mod**：`mod_mul(a, b, m)` / `mod_pow(a, e, m)`（e < 0 なら逆元の累乗）/ `mod_inv(a, m)`（合成数の mod でも）/ `try_mod_inv`、`geom_sum_mod(r, n, m)`（1 + r + … + r^(n−1)。逆元を使わないので m は何でも）/ `arith_sum(a, d, n)`（等差数列の和。Int128）/ `arith_sum_mod` / **`sum_range(l, r)` / `sum_range_closed(L, R)`**。**分数**：`frac_compare(a, b, c, d)`（a/b と c/d を厳密に）。**最終照合で追加・変更：`try_lcm(a, b)`（溢れたら ok = false）、`gcd_of` / `lcm_of_cap` を vector<int>・{…} でも、`sum_range_mod(l, r, m)` / `sum_range_closed_mod`** | int_common | ABC169 B、ABC227 C、典型90 020・022・038、ARC059 C |
| ★ `00_modint` | **作成済み（Phase 2）**。`ModInt<MOD>`（定数の mod。別名 `ModInt998` / `ModInt107`）、`DynModInt`（`set_mod(M)`。入力で決まる mod。`DynModIntT<ID>` で 2 つの mod を同時にも）。共通：int・long long・unsigned・Int128 から作る（負でも）/ `raw` / + − × ÷ 単項− ++ −−・複合代入（整数と混ざった式も）/ == != / `<`（map のキー用）/ `pow(e)`（e < 0 も）/ `inv()` / `try_inv` / `is_unit` / `val()` / `mod()` / `explicit operator bool` / `<<` / `>>`（長い数も読める）/ **`std::hash` の特殊化（`unordered_map<mint, V>` がそのまま使える。設計の ModIntHash から変更）** / `IS_PRIME` / `mod_is_prime()`。`rational_reconstruct(x)`（`RationalGuess{ok, num, den}`）（デバッグ用：mint の値を小さい分数 p/q に戻す。期待値などの答えを確かめる）。奇数の mod はモンゴメリ乗算（2^30 未満は遅延簡約）、偶数の mod は剰余（DynModInt は Barrett）。**最終照合で追加・変更：説明：`pow` の負の指数（N − 2 のような式）の注意** | math_base | 典型90 069、ABC280 E・ABC263 E（期待値を mod で。△） |

### 1. 桁・進数・回文数

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `01_digits` | **作成済み（Phase 3）**。**桁**：`digits_of(n, base)`（上の桁から。base ≤ 2^31 − 1）/ `digit_sum` / `num_digits` / `from_digits`。**進数**：`to_base(n, base)`（2〜36 進と −2〜−36 進。0・負の数も）/ `from_base(s, base)` / `from_base_cap(s, base, cap)` / `base_value_leq(s, base, limit)`（36 を超える進数も。溢れない）。**回文数**（k 番目は 0 から）：`is_palindrome_num(n, base)` / `kth_palindrome(k)`（文字列）/ `palindrome_rank` / `count_palindromes_up_to` / `_below` / `_of_length` / **`ceil_palindrome` / `floor_palindrome`**（設計の next / prev から改名）/ `make_palindrome(half, odd)`。**数をつなげる**：`concat_value`（Int128）/ `concat_mod` / `repeat_concat_mod` / `repunit_mod` / `pow10_mod`。**文字列の数**：`str_mod` / `str_num_compare` / **`substr_value`** / `DigitPrefix`。**数え上げ**：**`digit_length_blocks(l, r)`** / `total_digits_up_to` / `count_digit_up_to` / `digit_sum_total_up_to`（Int128）。**列挙**：`numbers_from_digits(ds, max_len)`。**最終照合で追加・変更：10^18 を超える回文数 `ceil_palindrome_str` / `floor_palindrome_str` / `count_palindromes_up_to_str` / `count_palindromes_below_str`（Int128。75 桁まで）/ `make_palindrome_str` / `is_palindrome_str`、`repeat_concat_value`、`str_mod(s, m, base)`、`carry_to_string(c, base)`（繰り上がり）** | math_base | ABC083 B、ABC105 C、ABC114 C、ABC186 C、ABC192 D、ABC336 C、ABC357 D、ABC363 D、ARC141 A、典型90 067・082 |

### 2. 素数・約数

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `02_sieve` | **作成済み（Phase 3）**。`Sieve S(N)`（線形篩。N ≤ 10^8）：is_prime / primes / spf / factorize / prime_factors / divisors / count_divisors / sum_divisors / phi / mu / omega / big_omega / largest_prime_factor / square_free_decompose / is_square_free / prime_count / count_primes(l, r) / kth_prime（0 番目が 2）/ lcm_mod。**N を超える数は 02_factorize に任せる**（同梱）。表：`phi_table` / `mu_table` / `divisor_count_table` / `divisor_sum_table` / `omega_table`。`primes_up_to(N)`（区間ごとのエラトステネス。奇数だけ。N < 2^31）。区間篩：`segment_primes(l, r)` / `segment_factorize(l, r, f)`（r ≤ 10^14、幅 ≤ 10^7。`_closed` 版も）。**最終照合で追加・変更：`Sieve(long long N)`、`S.primes()` と `primes_up_to` を vector<long long> に（`primes_up_to<int>` でメモリ半分）、`S.count_primes_closed`、`S.spf` の範囲の確認、表の関数は long long の N（0〜10^8 の外は止める）、`S.lcm_mod` を vector<int> でも、`divisor_count_table` / `divisor_sum_table` / `omega_table` を線形篩に（O(N)）** | math_base, **factorize** | ABC084 D、ABC177 E、ABC300 D、ABC383 D、ABC400 E、ABC412 E、典型90 030 |
| ★ `02_factorize` | **作成済み（Phase 3）**。1 つの大きな数（n < 2^63）：`is_prime(n)`（決定的ミラー・ラビン。モンゴメリ乗算）/ `factorize(n)`（ブレント版ポラード・ロー。`PrimePower{p, e}` の昇順）/ `factorize_trial` / `prime_factors` / `divisors_of` / `divisors_from_factors` / `for_each_divisor(fs, f)` / `count_divisors` / `sum_divisors`（Int128）/ `euler_phi` / `mobius` / `largest_prime_factor` / `square_free_decompose`（`SquareFreeParts{a, b}`）/ `is_square_free` / `radical` / **`is_prime_power`**。`factorial_exponent(n, p)` / **`factorial_max_power(n, k)`** / `min_factorial_multiple(k)`、`is_perfect_power` / `count_perfect_powers(N)`、`count_coprime_up_to(N, M)`、`lcm_mod_of(vs, m)`。**最終照合で追加・変更：`omega(n)` / `big_omega(n)`、`lcm_mod_of` を vector<int> でも** | math_base | ABC052 C、ABC169 D、ABC180 C、ABC190 D、ABC215 D、ABC254 D、ABC280 D、ABC284 D、ABC361 F、典型90 075・085 |
| ☆ `02_divisor_transform` | **作成済み（Phase 4）**。1〜N の上の `divisor_zeta(f, op)` / `divisor_mobius` / `multiple_zeta` / `multiple_mobius`（O(N log log N)。演算を渡せる）、`gcd_convolution` / `lcm_convolution`、`count_pairs_by_gcd(cnt)`、`phi_sum(N)`（Int128）/ `mertens(N)` / `count_coprime_pairs(N)`（杜教篩 O(N^{2/3})。N ≤ 10^13）、`sum_gcd_tuples(K, len, m)`（⌊K/d⌋ ごとにまとめて O(K^{2/3})）。**最終照合で追加・変更：`dirichlet_convolution(a, b)`、`DivisorsOf D(M)`（M ≤ 10^18 の約数の上の zeta / mobius / multiple_zeta / multiple_mobius）** | math_base, factorize, sieve | ABC162 E、ABC206 E、ABC428 G（ディリクレ積） |
| ☆ `02_prime_count` | **作成済み（Phase 4）**。`prime_count(N)`（Lucy。10^12 で約 1.3 秒）/ `prime_sum(N)`（Int128）/ `PrimeCountTable P(N)`（⌊N/i⌋ の形の v での π(v)）、`count_smooth(N, P)`（2 組に分けた半分全列挙。P ≤ 100 程度）。**最終照合で追加・変更：`multiplicative_sum<T>(N, coef, fpe)`（乗法的関数の和。Min_25 篩。N ≤ 10^12）** | math_base | ABC300 G、ABC172 D、ABC370 G（△） |

### 3. 合同式・mod の数論

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `03_mod_equation` | **作成済み（Phase 3）**。`linear_congruence(a, b, m)`（`ModSolution{ok, r, m, overflow}`）、`linear_diophantine(a, b, c)`（`Diophantine{ok, x0, y0, dx, dy}`、`d.x(k)` / `d.y(k)`。0 ≤ x0 < dx）/ `diophantine_count(a, b, c, xl, xr, yl, yr)` / **`diophantine_count_closed`** / `diophantine_min_x(a, b, c, x_min)`（`DiophantinePoint{ok, x, y}`）、`crt_merge` / `crt_all`（mod は互いに素でなくてよい。lcm が long long を超えたら overflow = true で、ok は「long long に入る解があるか」、r はその最小）、`garner(rs, ms, mod)`。使い方関数 `first_common_time`（周期的に繰り返す 2 つの区間が重なる最初の時刻） | math_base | ABC186 E、ABC193 E |
| ☆ `03_discrete_log` | **作成済み（Phase 4）**。`carmichael_lambda(n)`、`mult_order(a, m)` / `MultOrder ord(m)`（Maybe）、`primitive_root(p)`、`mod_sqrt(a, p)`（Maybe。小さい方）、`discrete_log(a, b, m)`（Maybe。m は合成数でもよい。gcd で減らしてから BSGS。m ≤ 10^13 程度）、`pow_tower_mod(as, m)` | math_base, factorize | ABC222 G、ABC228 E、ABC335 G |

### 4. 組合せ・順列

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `04_binomial` | **作成済み（Phase 3）**。`Binomial<mint> B(N)`（mod は素数。足りなければ 1.5 倍ずつ自動で伸ばす。DynModInt の mod が変わったら作り直す）：fact / inv_fact / inv / C / P / H / multinomial（vector か {…}）/ catalan / paths / paths_not_crossing（ずっと y ≤ x + k。鏡像法）。`binom_large_n<mint>(n, r)` / `binom_lucas(n, r, p)`（p ≤ 10^7 の素数）/ `binom_mod_any(n, r, m)`（素数べきごとの表と中国剰余定理。素数べきの部分 ≤ 10^7）/ `binom_exact` / `binom_cap` / `pascal_table<T>(N)`。**最終照合で追加・変更：`BinomialPrefix<mint>`（Σ_{i≤m} C(n, i) を n・m を ±1 して O(1)。Mo のアルゴリズム用）** | math_base, modint | ABC110 D、ABC132 D、ABC145 D、ABC156 D、ABC167 E、ABC202 D、ABC205 E、ABC425 E、ARC117 C、ABC463 G |
| ☆ `04_combinatorics` | **作成済み（Phase 4）**。`stirling2<mint>(n, k)` / `stirling2_row` / `stirling2_table<T>`、`stirling1_row` / `stirling1_table<T>`、`bell` / `bell_table<T>`、**`partition_row<mint>(N)`（級数の逆数で O(N log N)）** / `partition_table<T>` / `partition_count_table<T>`、`derangement_table<T>`、`surjections`、`balls_into_boxes(n, k, 玉, 箱, rule)`（写像 12 相）、`inclusion_exclusion<T>(n, f)` / **`inclusion_exclusion_count<mint>(n, g)`**、`count_divisible_by_any(l, r, ds)` / `_closed`。_table<T> の整数型は溢れたら止める。**最終照合で追加・変更：`binomial_transform` / `binomial_inverse` / `at_least_from_exactly` / `exactly_from_at_least`、`for_each_set_partition` / `for_each_integer_partition`、`q_binomial<mint>(n, k, q)`（q-Lucas）、`young_tableaux_count<mint>(shape)`（フック長公式）、inclusion_exclusion の向きの説明** | math_base, modint, binomial, convolution | ABC131 C、ABC242 F、ABC297 F、ABC390 D、ABC378 G（△） |
| ☆ `04_permutation` | **作成済み（Phase 4）**。`inversion_count(a)`、`perm_rank(p)`（Int128。n ≤ 33）/ `perm_rank_mod<mint>(p)` / `perm_unrank(n, k)`、`multiset_perm_rank(a)` / `multiset_perm_rank_mod<mint>(a)` / `multiset_perm_unrank(a, k)`、`kth_combination(n, r, k)` / `combination_rank(n, c)`。列は vector でも string でもよい。番号は 0 から。**最終照合で追加・変更：`perm_parity(a)`** | math_base, modint, binomial | ABC190 F、ABC202 D、ABC276 C、ABC470 F |

### 5. 行列・線形代数

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `05_matrix` | **作成済み（Phase 3）**。`Matrix<T>`（mint・long long（溢れたら止める）・double・半環）：identity / zero / (i, j) / [i][j] / + − × スカラー倍 / pow(e) / apply(v) / transpose / trace / `<<`。半環 `MinPlus<T>` / `MaxPlus<T>` / `XorAnd` / `OrAnd`。`MatrixFixed<T, K>`（同じ演算と det / inv / try_inv）。`companion_matrix(c)` / **`recurrence_nth_matrix(c, a, n)`**（漸化式の第 n 項。半環でも）。**最終照合で追加・変更：半環 `MinMax<T>` / `MaxMin<T>`、自作の半環の作り方の説明** | math_base, modint | ABC009 D、ABC293 E、EDPC R |
| ☆ `05_linear_algebra` | **作成済み（Phase 4）**。体（mint・double・Fraction。整数型はコンパイルエラー）の `mat_rank` / `mat_det` / `mat_inverse`（Maybe）/ `solve_linear(A, b)`（`LinearSolution{ok, x, kernel}`）、`det_mod_any(A, m)`、`count_spanning_trees<T>(n, us, vs)` / `count_arborescences<T>(n, us, vs, root)`、`char_poly(A)`（ヘッセンベルグ）/ **`det_affine(A, B)`（det(A + xB)。ABC323 G）** | math_base, modint, matrix | ABC253 Ex、ABC323 G |
| ☆ `05_xor_basis` | **作成済み（Phase 4）**。`XorBasis`（簡約した基底）：insert / contains / reduce / max_xor / min_xor / kth / **count_le** / rank / count / zero_by_nonempty / represent / merge / basis、`PrefixXorBasis`（push_back と [l, size()) のクエリ）、`Gf2System`（ビットで持つ。add / add_row / rank / solvable / free_count / solve / kernel） | math_base | ABC141 F、ABC223 H、ABC283 G |

### 6. 数列の和・漸化式

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `06_quotient` | **作成済み（Phase 3）**。`for_each_quotient(N, f)` / `quotient_blocks(N)`（`QuotientBlock{q, l, r}`）/ `quotient_blocks_ceil(N)` / `quotient_blocks2(N, M)`（i ≤ min(N, M)。`QuotientBlock2{q1, q2, l, r}`）。N < LLONG_MAX。`sum_floor_div(N)` / `sum_i_floor_div(N)` / `sum_mod_i(N)`（Int128）、`linear_floor_sum(n, m, a, b)`（a, b は負でも。Int128。n, m < 2^32 なら 64 ビットで速く）/ `count_lattice_under_line(a, b, c)`。**最終照合で追加・変更：`for_each_quotient` を途中で止める（f が false を返す）、`linear_mod_sum(n, m, a, b)`、`linear_floor_sum_ext<mint>`（Σ⌊⌋・Σ i⌊⌋・Σ⌊⌋²）、`universal_euclid`（万能ユークリッド）、引数の順が ACL と同じという説明** | math_base | ABC044 D、ABC172 D、ABC230 E、ABC372 G、ACL Practice C、ABC402 G・ABC429 G（universal_euclid で） |
| ★ `06_series` | **作成済み（Phase 3）**。`sum_powers<mint>(n, k)`（Σ_{i<n} i^k。0^0 = 1。O(k)。`(l, r, k)` 版と `sum_powers_closed`）、`lagrange_interpolate(ys, x)`（O(d)）/ `lagrange_interpolate_at(xs, ys, x)` / `poly_from_points(xs, ys)`（O(d²)）、`arith_geom_sum(r, n)`（Σ i·r^i。O(log n)。mod は何でも）、`poly_eval` / `poly_mul_naive` / `poly_divmod`（`PolyDivMod{q, r}`。整数なら割り切れないと止める） | math_base, modint | ABC208 F、ABC245 D、ARC033 D |
| ☆ `06_linear_recurrence` | **作成済み（Phase 4）**。`linear_recurrence_kth(a, c, K)`（Bostan–Mori。NTT が使える mod では Q(−x) を値の入れ替えで）/ `linear_recurrence_terms(a, c, n)` / `berlekamp_massey(s)` / `guess_kth(s, K)` | math_base, modint, convolution | ABC300 Ex |

### 7. 畳み込み・多項式

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ☆ `07_convolution` | **作成済み（Phase 4）**。`convolve(a, b)`（NTT は 4 進。NTT が使えない mod は 3 つの素数と Garner で自動に）、`convolve_ll` / `convolve_i128`（必要なら 5 つの素数）、`convolve_mod(a, b, m)`、`convolve_many(ps)`、`correlate(a, b)`、`ntt` / `intt`。**最終照合で追加・変更：`convolve_cyclic(a, b, n)`、`RelaxedConvolution<mint>`（オンラインの畳み込み）** | math_base, modint | ABC196 F、ABC291 G、ACL Practice F、ATC001 C、ABC357 G（オンラインの畳み込み） |
| ☆ `07_fps` | **作成済み（Phase 4）**。`Fps<mint>`（+ − × スカラー倍 / pre / eval / `<<`）、`fps_inv` / `fps_log` / `fps_exp` / `fps_pow`（k ≤ 10^18）/ `fps_sqrt`（Maybe）/ `fps_divmod`（`FpsDivMod{q, r}`）/ `fps_derivative` / `fps_integral` / `taylor_shift` / `fps_coeff_of_fraction(P, Q, K)` / `prod_linear`。**最終照合で追加・変更：`multipoint_eval` / `interpolate`、`power_sums(as, n, ws)`、`fps_sum_of_fractions`** | math_base, modint, convolution | ABC230 Ex |
| ☆ `07_subset_transform` | **作成済み（Phase 4）**。`subset_zeta` / `subset_mobius` / `superset_zeta` / `superset_mobius`（演算を渡せる）、`walsh_hadamard`、`or_convolution` / `and_convolution` / `xor_convolution`、`subset_convolution`（O(n^2 2^n)）。**最終照合で追加・変更：`subset_exp` / `subset_log`（割り算なし）** | math_base, modint | ABC212 H |

### 8. 分数・多倍長

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `08_fraction` | **作成済み（Phase 3）**。`FractionT<T>`（`Fraction` = long long、`Fraction128` = Int128）：常に既約・分母 > 0、+ − × ÷（途中は 128 ビット。結果が収まらなければ止める）、比較（long long は 128 ビットの掛け算、Int128 は掛け算をしない連分数の方法）/ `compare(x, y)`、floor / ceil / round / sign / abs / reciprocal / is_integer / to_real / `to_decimal(k)`（切り捨て）、`<<` / `>>` / `parse`、**`std::hash` の特殊化**（設計の FractionHash から変更）、`best_rational_approx(p, q, max_den)`。**最終照合で追加・変更：`rational_floor` / `rational_ceil`（分母が上限以下で以下の最大・以上の最小）、`simplest_between(x, y)`、`continued_fraction` / `from_continued_fraction`** | math_base | ABC308 C、ABC408 G |
| ☆ `08_bigint` | **作成済み（Phase 4）**。`BigInt`（10^9 進・符号つき）：整数・Int128・文字列から / + − × ÷ % ++ −−（× はカラツバ。長さが違うときは区切る。÷ は Knuth）/ 比較 / divmod / pow / isqrt（上半分から始めるニュートン法）/ gcd / abs / to_string / fits_ll / to_ll / to_i128 / to_ld / num_digits / digit_sum / mod_ll / `<<` / `>>`。**最終照合で追加・変更：`lcm(a, b)`（整数と混ぜても）** | math_base | ABC169 B |

### 9. ビット

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `09_bit` | **作成済み（Phase 3）**。`bit_count` / `top_bit` / `low_bit`（位置。0 なら −1）/ `bit_length` / `has_bit` / `bit_indices` / `for_each_bit` / `bit_reverse` / `mask_to_string`、`for_each_submask(mask, f, include_zero, include_self)`（小さい順）/ `for_each_supermask` / `for_each_mask_between(must, allowed, f)` と 1 の数 k を決める版 / `for_each_subset_k(n, k, f)`（ゴスパー）、`gray_code` / `gray_inverse`、`xor_up_to` / **`xor_range` / `xor_range_closed`** / `count_bit_up_to` / `popcount_sum_up_to`（Int128）、`xor_pair_sum` / `subarray_xor_sum`（Int128）。**最終照合で追加・変更：`xor_pair_sum` / `subarray_xor_sum` を vector<int> などでも** | math_base | ABC121 D、ABC147 D、ABC197 C、ABC238 D、ABC269 C、ABC356 D、ABC365 E |

### 10. ゲーム

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ☆ `10_game` | **作成済み（Phase 4）**。`mex(vs)`、`grundy_dag(n, moves)`（再帰しない。閉路は止める）、`nim_first_wins` / `misere_nim_first_wins`、`range_take_grundy(x, L, R)`、`subtraction_game_grundy(S, N)` / `SubtractionGame G(S)`（周期を見つけて巨大な x も）。**最終照合で追加・変更：`grundy_dag_multi(n, moves)`（1 手で複数のゲームに分かれる）、`game_outcome(n, moves)`（閉路のあるゲームの後退解析。`GameOutcome{result, steps}`）** | math_base | ABC206 F、ABC297 G、ABC368 F、ABC209 E |

### △（型）

| 型 | 置き場所 |
|---|---|
| 期待値・確率を mod で求める | 00_modint の使い方関数（ABC280 E・ABC263 E の形） |
| 包除原理（部分集合の符号付き和） | 04_combinatorics の `inclusion_exclusion` と使い方関数（ABC242 F・ABC297 F の形） |
| bit 全探索 | 09_bit の使い方関数（ABC197 C の形） |
| 答えの二分探索で値が大きくなる（溢れ） | 00_math_base の `pow_cap` / `mul_cap`（二分探索そのものは binary_search.cpp） |
| 周期的に繰り返す区間の重なり（時刻の合同式） | 03_mod_equation の使い方関数（ABC193 E の形） |

### 対象外（決定済み）

- 幾何（geometry_v2）
- 二分探索・三分探索（binary_search_template。三分探索も入っている）
- 累積和・いもす法・座標圧縮
- DP（ナップサック・桁 DP・bit DP の本体。部分集合の列挙だけは 09_bit）
- データ構造（BIT・セグメント木。転倒数・順列の番号に使う BIT はテンプレートの中に持つ）
- 文字列・区間・乱数・ハッシュ


## 2-P2. Phase 2（基盤）の実装メモ（2026-10-10）

### 作ったファイル（`Desktop\claude\math_v2\`）

| ファイル | 内容 | テスト（主な愚直解。Release の照合回数） |
|---|---|---|
| `00_int_common.cpp` | 数学・幾何で共通の整数の部分（Int128 / UInt128 の入出力、isqrt、parse_fixed） | 割り算を使わない 10 進変換、平方数の前後、long long の端の小数、ストリームの入出力（273 万回） |
| `00_math_base.cpp` | 数学基盤 | 128 ビットの愚直計算（引き算で余りをそろえる割り算、飽和する累乗、ループで数える倍数）、ext_gcd の解の大きさ、mod_pow は左から右の別の書き方、等比数列の和は小さい n をループ・大きい n を倍化の式で（1137 万回） |
| `00_modint.cpp` | ModInt / DynModInt / rational_reconstruct | 13 個の定数の mod（素数・合成数・偶数・2^30 の前後・1・2）と 78 個の入力の mod で全演算を 128 ビットの計算と比べる、2 万回の演算の連鎖（遅延簡約の値がずれないか）、ハッシュ・map・sort、2 つの mod を同時に、小さい分数の復元（253 万回） |
| `update_bundles.ps1` / `.bat` | graph_grid_v2 のものをそのままコピー | － |
| `dynalist_items_phase2.txt` | Dynalist に貼る項目（3 本） | － |

確認したこと：g++ 13（-O2 -Wall -Wextra -Wshadow -pedantic）と clang 18（-pedantic-errors）で C++17 / 20 / 23 とも警告なし・全テスト通過、`-D_GLIBCXX_DEBUG -fsanitize=address,undefined,float-cast-overflow` で全体（テスト・速度計測・デモ）が通過、Main_base.cpp（ACL 込み）に**グラフ・グリッド 46 本・幾何 16 本＋00_int_common・数学 3 本**の本体を 2 通りの順で全部貼り、gnu++2b / c++17・LOCAL あり / なしで g++・clang の両方でコンパイル・実行（名前の衝突なし。ACL の modint998244353・pow_mod・inv_mod・crt・floor_sum も今までどおり呼べる、`debug(mint)`・`debug(vector<mint>)` も可）、数学の各ファイルを 1 本ずつ Main_base に貼って -Wshadow の警告なし、ほかの単独テンプレート（binary_search.cpp・static_seq_toolkit.cpp・segtree・lazy_segtree・dynamic_unique_seq.cpp）のグローバルな名前とも重ならない、同梱部分が持ち主のファイルと全部一致（update_bundles -Check と同じ確認）。PC（`Desktop\claude\math_v2`・`geometry_v2`、g++ 11.4）に保存して、数学 3 本・幾何 17 本ともビルド（警告 0）・ランダムテストが通ることを確かめた。

保存したもの：更新前の幾何の 16 本と設計書は `Desktop\claude\geometry_v2_backup_20261010\` に残した。組み立て用の作業ファイル（本体・テスト・デモの元と build.py）は `math_v2\_work\` に置いた（update_bundles はサブフォルダを見ないので影響しない）。

### 幾何第二世代の更新（6-1 の確認 2：案 A）

- `00_geo_base` の本体から `Int128` の typedef・`isqrt` / `isqrt_ceil`・`parse_fixed` / `read_fixed`・`int128_to_string` を取り除き、16 本すべてに `GG_INT_COMMON` のブロックを `GG_GEO_BASE` の前に同梱した。geometry_v2 にも同じ `00_int_common.cpp` を置いた（17 本。update_bundles がそのまま使える）
- 関数の名前・動きは同じ（変わったのは、Int128 を `cout <<` で出せるようになったこと、`parse_fixed` が long long に収まらない値で止まるようになったこと）。各ファイルの見出しの「依存」に 00_int_common を足し、00_geo_base の説明・デモの「cout で直接出せない」を直した
- 確かめたこと：17 本とも g++（-pedantic）・clang（-pedantic-errors）・C++20 / C++23 で警告なし、ランダムテストの照合回数とデモの出力が更新前と完全に同じ（変えた 2 行の説明を除く）、サニタイザ ＋ _GLIBCXX_DEBUG、Real を double にした版、同梱部分 48 個が持ち主と一致

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）

| 処理 | 規模 | 時間 |
|---|---|---|
| isqrt / Int128 の出力 / Int128 の入力 / read_fixed(4) | 2×10^6 / 2×10^5 / 2×10^5 / 2×10^5 回 | 12 / 71 / 40 / 11 ms |
| floor_div + mod_floor / kth_root_floor(·, 3) / ext_gcd / pow_cap | 2×10^6 回 | 19 / 36 / 141 / 44 ms |
| mod_pow（m = 10^9+7）/ mod_pow（m ≈ 4.6×10^18）/ geom_sum_mod（指数は 63 ビットまで） | 5×10^5 回 | 86 / 173 / 211 ms |
| ModInt998：x = x·y + c の連鎖 / 内積 / 階乗と逆元の表 N = 10^7 / pow + inv | 2×10^7 / 1.6×10^7 / 1 / 10^6 回 | 100 / 27 / 109 / 211 ms |
| DynModInt（998244353）：同じ | 同じ | 84 / 24 / 121 / 208 ms |

ACL・旧テンプレートとの比較（同じ計測環境。ms）：

| 型 | 連鎖 1×10^8 | 内積 6.7×10^7 | O(n²) n = 4000 | 階乗の表 10^7 | pow 10^6 | inv 10^6 |
|---|---|---|---|---|---|---|
| ACL static_modint<998244353> | 527 | 106 | 19 | 110 | 105 | 118 |
| 旧 FastModInt<998244353> | 524 | 107 | 23 | 111 | 102 | 93 |
| **ModInt998** | 525 | 113 | 26 | 117 | 103 | 108 |
| ACL dynamic_modint（998244353） | 499 | 149 | 38 | 138 | 115 | 118 |
| 旧 DynModInt（998244353） | 441 | 95 | 19 | 108 | 121 | 98 |
| **DynModInt**（998244353） | 472 | 103 | 20 | 117 | 111 | 103 |
| ACL dynamic_modint（2^31 − 1） | 492 | 100 | 19 | 121 | 117 | 130 |
| **DynModInt**（2^31 − 1。Barrett） | 536 | 103 | 20 | 128 | 138 | 132 |

- 定数の mod はこの環境では ACL（コンパイラの定数の剰余）・旧 FastModInt（モンゴメリ）とほぼ同じ。旧テンプレートの記録（Skylake）ではモンゴメリが連鎖で 1.44 倍速かったので、モンゴメリのままにした
- 入力で決まる mod は ACL より速い（内積で約 1.4 倍、O(n²) の畳み込みで約 1.9 倍）。旧 DynModInt とほぼ同じ速さで、旧 DynModInt にはなかった偶数の mod・2^30 以上の mod にも対応した

### 実装で決めたこと

- **00_int_common** を作った（6-1 の確認 2）。Int128 / UInt128 は `__extension__ typedef`（-pedantic でも警告なし）。`operator<<` / `operator>>` を足した。Main_base の `debug(x)` は Int128 をそのままでは表示できない（debug の中の `<<` が後から書いたこの `<<` を見つけられない。`namespace std` に入れれば見つかるが規格違反なのでやめた）→ `debug(int128_to_string(x))` と説明に書いた
- **整数の関数の約束**：使える範囲の外の引数（d ≤ 0 の倍数、m ≤ 0 の mod、pow_cap に負の数など）はメッセージを出して止める。厳密な値の関数は溢れたら止める（pow_int・arith_sum・sum_range・ceil_multiple など）。LLONG_MIN / -1 は CPU が止める（説明に書いた）
- **ext_gcd** は |a|, |b| < 2^62 なら long long のまま（2^62 以上だけ Int128。Int128 の割り算は遅いため）。返す解は |x| ≤ max(1, |b|/g)、|y| ≤ max(1, |a|/g)
- **kth_root_floor** は double の pow で見当をつけ、溢れない掛け算で前後を直す（long double の pow は遅かった：2×10^6 回で 531 ms → 36 ms）
- **mod_pow / geom_sum_mod** は m ≤ 2^32 なら 64 ビットの掛け算、それより大きいときだけ 128 ビット
- **round_div** は「.5 は大きい方へ」（-2.5 → -2）。四捨五入の流儀は問題によるので説明に書いた
- **`Maybe<T>{ok, value}`** を基盤に置いた（try_mod_inv・try_inv が返す。ほかのテンプレートの「ないことがある値」もこれにそろえる）
- **ModInt** の内部：奇数の MOD はモンゴメリ乗算（MOD < 2^30 なら値を [0, 2·MOD) のまま持つ遅延簡約）、偶数の MOD は剰余。逆元は MOD が素数なら pow(MOD − 2)（拡張ユークリッドより速かった）、合成数なら拡張ユークリッド。0 の逆元・互いに素でない数の逆元はメッセージを出して止める（assert と違い、NDEBUG でも止まる）
- **DynModInt** の内部：mod が奇数で 2^30 未満ならモンゴメリ（遅延簡約）、それ以外は Barrett。状態は `inline static` のメンバ（関数の中の static と違い、使うたびの初期化の確認がない）。最初の mod は 998244353
- **ハッシュ**は `std::hash` の特殊化にした（`unordered_map<mint, V>` と書くだけで使える。設計の `ModIntHash` はやめた）
- **入力**は文字列で読んで 1 桁ずつ mod を取る（long long を超える長さの数も読める）
- 設計になかったが足したもの（迷ったら入れる）：`ceil_multiple` / `floor_multiple` / `ceil_congruent` / `floor_congruent`、`add_overflows`、`ilog_ceil`、`sum_range` / `sum_range_closed`、`parse_uint128` / `uint128_to_string`、`ModInt::IS_PRIME`、`DynModInt::mod_is_prime()`
- テストとデモの形はグラフ・グリッド・幾何と同じ（`xxx_test::random_test()`（Release は 8 倍の量）・`benchmark()`・`xxx_demo::demo()`・`setup_console()`）。テスト・デモの共通部品（CHECK・乱数・to_s・SHOW など）は組み立てのときに差し込む（作業用の仕組み。できた .cpp はほかの第二世代と同じ形）
- デモの例：ABC169 C（read_fixed）、ABC131 C（count_multiples_closed ＋ lcm_cap）、ABC227 C（kth_root_floor）、ABC169 B（mul_cap）、典型90 020（pow_cap）・038（lcm_cap）、ABC228 E（mod_pow）、ABC293 E・ABC357 D（geom_sum_mod）、ABC308 C（frac_compare）、典型90 069（pow）、ABC280 E（期待値と rational_reconstruct）

### 次回以降の作業のしかた

- 本体の正は各 .cpp の「ここからテンプレート本体（コピー範囲）」の中。基盤を直したら `update_bundles.bat` を実行する
- **00_int_common を直したときは、math_v2 と geometry_v2 の両方の 00_int_common.cpp を同じにしてから**、それぞれのフォルダで `update_bundles.bat` を実行する

## 2-P3. Phase 3（★ 10 本）の実装メモ（2026-10-10）

### 作ったファイル（`Desktop\claude\math_v2\`）

| ファイル | テスト（主な愚直解。Release の照合回数） |
|---|---|
| `01_digits.cpp` | 文字列での n 進変換・負の進数の値の計算、0〜2×10^6 の回文数を全部並べた表、128 ビットでの連結、1〜3×10^5 を全部書いた数え上げ、大きい N は 1 つ前との差（462 万回） |
| `02_factorize.cpp` | 試し割り・素数の表、sympy と同じ値の例、累乗数は全部数える、約数は全部割ってみる（336 万回） |
| `02_sieve.cpp` | 愚直な篩・試し割り、区間篩は 1 つずつ素数判定（10^14 の近くも）、ブロック（2^18）の境目（176 万回） |
| `03_mod_equation.cpp` | 解を全部調べる（小さい mod）、128 ビットで検算、答えを先に決めて作る CRT（lcm が long long を超える場合も）、互いに素な mod の garner（573 万回） |
| `04_binomial.cpp` | 128 ビットのパスカルの三角形（n ≤ 128）、DP の格子路、いろいろな素数の mod（2・3・13・2^31−1 も）、mod を変えた DynModInt、Lucas と合成数の mod を互いに（449 万回） |
| `05_matrix.cpp` | 愚直な掛け算、ちょうど K 歩の DP（最短・最長・行けるか）、漸化式を直接計算、行列式は置換で展開（66 万回） |
| `06_quotient.cpp` | 1 つずつ数える（N ≤ 2000）、大きい N は 1 つ前との差（約数の個数・約数の和）、floor_sum は 128 ビットで全部足す（147 万回） |
| `06_series.cpp` | 1 つずつ足すべき乗和、係数から作った値での補間、多項式の割り算は掛け戻して確かめる（68 万回） |
| `08_fraction.cpp` | 128 ビットの掛け算での検算、小数表示は愚直な割り算、近い分数は分母を全部試す（88 万回） |
| `09_bit.cpp` | 全部の mask を調べる列挙、0〜3×10^5 の累積、全部の組・区間の XOR の和（528 万回） |
| `dynalist_items_phase3.txt` | Dynalist に貼る項目（10 本） |
| `dynalist_tool_index.txt` / `.opml` | 数学道具選択の索引（Phase 2・3 の分。Phase 4 のファイルは「作成予定」と書いた。Phase 4 で完成） |

確認したこと（Phase 2 の 3 本も含めて全 13 本）：
- g++ 13 と clang 18（-O2 -Wall -Wextra -Wshadow -pedantic）で C++17 / 20 / 23 とも警告なし・全テスト通過
- `-D_GLIBCXX_DEBUG -fsanitize=address,undefined` と、テストの量を減らさない `-O2 -fsanitize=undefined` の両方で全体（テスト・速度計測・デモ）が通過（ここでテストの愚直解の溢れ 2 か所と、linear_floor_sum の速い経路の溢れ 1 か所を見つけて直した）
- Main_base.cpp（ACL 込み）にグラフ・グリッド 46 本・幾何 17 本・数学 13 本の本体を 2 通りの順で全部貼り、gnu++2b / c++17・LOCAL あり / なしで g++・clang の両方でコンパイル・実行（名前の衝突なし。Phase 3 の関数も呼んで確かめた）。数学の各ファイルを 1 本ずつ Main_base に貼っても警告なし
- 同梱部分 27 個が持ち主のファイルと全部一致
- **独立したレビュー**（別のエージェントが本体だけを読み、自分で書いた愚直解・Python の多倍長整数・sympy と比べた。2 つに分けて 5 本ずつ）：使える範囲の入力で間違った答えはなし。範囲の外の入力や long long の端の値で、メッセージを出して止まらずに未定義動作になるものなどを直した（下の表）
- **代表的な問題を Main_base ＋ テンプレートの本体だけで解いて、全部の例が合うことを確かめた**（24 問：ABC186 E、ABC193 E、ABC156 D、ABC205 E、ABC425 E、ARC117 C、EDPC R、ABC293 E、ABC009 D、ABC230 E、ACL Practice C、ABC172 D、ABC208 F、ABC245 D、ARC033 D、ABC308 C、ABC269 C、ABC356 D、ABC121 D、ABC365 E、ABC363 D、ABC084 D、ABC169 D、ABC110 D）。デモにもこのほか ABC083 B・ABC336 C・ABC186 C・ABC192 D・ABC357 D・ABC114 C・典型90 067・082・ABC132 D・ABC145 D・ABC167 E・ABC202 D・ABC044 D・ABC147 D・ABC197 C・ABC238 D の例を入れ、答えが問題の例と合うことを確かめた

### 独立したレビューで直したもの

| ファイル | 内容 | 直し方 |
|---|---|---|
| 02_sieve | `Sieve S(N)` で N > 10^8 のとき、止まる前に表を確保していた（N = 2^31−1 で溢れ） | 大きさを確かめてから確保 |
| 03_mod_equation | `garner` で 0 の mod があると、確かめる前に 0 で割っていた | 最初に全部の mod を確かめる |
| 01_digits | `DigitPrefix` / `concat_mod` / `repeat_concat_mod` で m = 0 のとき、確かめる前に % m | 最初に m ≥ 1 を確かめる |
| 01_digits | `numbers_from_digits` が 0〜9 以外の数字を受け付けていた／`substr_value`・`DigitPrefix::value` が範囲の外を読んでいた | 止める |
| 02_factorize・02_sieve | `lcm_mod_of` / `S.lcm_mod` が負の数を受け付けていた（LLONG_MIN で溢れ） | 負の数は止める |
| 06_quotient | `count_lattice_under_line(1, 1, LLONG_MAX)` で ⌊c/a⌋ + 1 が溢れた | x = 0 の列を分けて数える |
| 06_quotient | `for_each_quotient` などで N = LLONG_MAX のとき右端 N + 1 が溢れた／`quotient_blocks_ceil` の (N + l − 1)/l が溢れた／`sum_i_floor_div` で N > 6×10^18 のとき l + r − 1 が溢れた | N < LLONG_MAX と決めて止める／(N − 1)/l + 1／Int128 で |
| 05_matrix | `MatrixFixed<long long>` の引き算・スカラー倍が溢れを調べていなかった | Matrix と同じく調べる。`s * M`・`M *= s` も足した |
| 08_fraction | `Fraction128` の足し算で途中の値が −2^127 のとき符号を反転して溢れた | gcd の絶対値を符号なしで |
| 08_fraction | `to_decimal` が分母 > 2^127/10 で止まった／`parse` が "." や "+-1.5" を読んでいた／2^127 以上の UInt128 が負になった | 10 r を足し算で／形を確かめる／止める |
| 04_binomial | `binom_mod_any` で 10^7 を超える素因数がある m のとき、止まる前に √m まで試し割りした（p·p も溢れた） | 10^7 まで試したら止める |

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）

| 処理 | 時間 |
|---|---|
| is_prime（10^18 程度）10^5 回 / factorize（ランダム 10^18）10^4 回 / factorize（30 ビットの素数 2 つの積）1000 回 | 19 / 130 / 404 ms |
| Sieve(10^7) / S.factorize を 1〜10^7 全部 / primes_up_to(10^8) / segment_primes（10^14 から幅 10^7） | 75 / 751 / 164 / 151 ms |
| linear_congruence 10^6 回 / crt_all（10^6 個）＋ crt_merge 10^6 回 / garner（5 個）2×10^5 回 | 230 / 48 / 132 ms |
| Binomial(10^7) の表 / C 10^7 回（ランダムな位置）/ binom_large_n(10^18, 10^7) / binom_mod_any（m = 10^9）10^5 回 | 111 / 357 / 37 / 53 ms |
| Matrix<mint> 300×300 の掛け算 / 50×50 の 10^18 乗 / MatrixFixed<mint, 2> の累乗 10^5 回 | 37 / 14 / 182 ms |
| sum_floor_div(10^14) / sum_i_floor_div(10^14) / linear_floor_sum 10^6 回 | 27 / 169 / 206 ms |
| sum_powers(10^18, k = 2.5×10^6) / lagrange_interpolate（10^6 点） | 112 / 24 ms |
| Fraction 10^6 個の sort / 四則演算 10^6 回 / Fraction128（≈ 2^100）2×10^5 個の sort | 97 / 142 / 325 ms |
| 全 mask の for_each_submask（3^16 回）/ xor_pair_sum + subarray_xor_sum（n = 10^6） | 27 / 135 ms |

### 実装で決めたこと

- **02_sieve は 02_factorize を同梱する**（設計では math_base だけ）。`Sieve S(N)` に N を超える数が来たら factorize / is_prime に任せる（クエリの値が N を少し超えても止まらない）。factorize・divisors_of なども sieve を貼れば使える
- **素因数分解の形**は `PrimePower{p, e}` の p の昇順（`<<` で "p^e"）。平方因子の分解は `SquareFreeParts{a, b}`（n = a²·b）
- **回文数**：設計の `next_palindrome` / `prev_palindrome` は「n 以上・以下」（n 自身も含む）なので `ceil_palindrome` / `floor_palindrome` に改名した。`kth_palindrome` は 10^18 を超えるので文字列
- **digits_of** は桁を int で返すので base ≤ 2^31 − 1（digit_sum / num_digits は base がいくつでもよい）
- **crt** は lcm が long long を超えても止めない：`overflow = true` にして、ok は「long long に入る解があるか」、r はその最小の解（lcm が 2^63 を超えたあとは「今の r がそのまま解か」だけを見ればよいので厳密）
- **一次不定方程式**の解は x0 を 0 以上で最小にそろえた。y0 が long long に入らないときは linear_diophantine は止まるが、diophantine_count / diophantine_min_x は内部を 128 ビットで計算するので止まらない
- **Binomial<mint>** は mod が素数かを最初に確かめ、n ≥ mod なら止める（Lucas を案内）。DynModInt の mod が変わったら表を作り直す（1 回の比較なので速さは変わらない）
- **行列**は半環の型に `zero()` / `one()` を持たせ、普通の数（mint・long long・double）と同じ Matrix で扱う。long long は掛け算・足し算の溢れを調べて止める。漸化式の第 n 項 `recurrence_nth_matrix` を足した（ABC009 D がそのまま書ける）
- **linear_floor_sum** は ACL の floor_sum と同じ形。a, b は負でもよく、n, m < 2^32 なら 64 ビットで計算、それより大きいときは Int128 で途中の溢れを調べる
- **分数**は long long の比較を 128 ビットの掛け算、Int128 の比較を掛け算をしない連分数の方法にした（溢れも誤差もない）。gcd は 128 ビットでは二進 GCD（Int128 の割り算は遅い）。足し算・掛け算は先に約分して余分な gcd を省いた。ハッシュは `std::hash` の特殊化（設計の FractionHash はやめた。modint と同じ）
- **ビット**の列挙はすべて小さい順（ABC269 C の答えがそのまま出る）。`xor_range` / `xor_range_closed` を足した
- 設計になかったが足したもの：`substr_value`・`digit_length_blocks`（01）、`is_prime_power`・`factorial_max_power`（02）、`diophantine_count_closed`・`DiophantinePoint`（03）、`recurrence_nth_matrix`・`MatrixFixed` の `s * M`（05）、`sum_powers(l, r, k)`・`sum_powers_closed`（06）、`compare`・`is_integer`・`parse`（08）、`xor_range`・`xor_range_closed`・`popcount_sum_up_to`（09）

## 2-P4. Phase 4（☆ 13 本）の実装メモ（2026-10-11）

### 作ったファイル（`Desktop\claude\math_v2\`）

| ファイル | テスト（主な愚直解。Release の照合回数） |
|---|---|
| `07_convolution.cpp` | 愚直な掛け算（998244353・10^9+7・167772161、DynModInt で 2^31−1・10^9・2^30・7・2・1 などの合成数や偶数の mod）、128 ビットの愚直（convolve_i128・convolve_ll・convolve_mod）、convolve_many・correlate・ntt → intt（2.8 万回） |
| `07_fps.cpp` | 1 項ずつ決める漸化式（inv・log・exp）、f^k は掛け算の繰り返し、√ は 2 乗して戻す、割り算は掛け戻す、テイラーシフトは点での値、P/Q の係数は級数の展開とフィボナッチの倍化公式、小さい素数の mod（6,886 回。1 回が長さ数百の列） |
| `07_subset_transform.cpp` | 全部の組を調べる OR・AND・XOR・部分集合の畳み込み、ゼータ → メビウスで戻る（4 万回） |
| `02_divisor_transform.cpp` | 約数・倍数の和を全部、全部の組の gcd、φ・μ の表の和と 1 つ前との差、列を全部並べた gcd の和と φ の表での直接計算（2.7 万回） |
| `02_prime_count.cpp` | 篩で数えた値・よく知られた値（π(10^12) など）、試し割りでの smooth 数の数え上げ（4,800 回） |
| `03_discrete_log.cpp` | 1 つずつ掛ける位数と λ、平方根・離散対数は全部試す（p, m ≤ 600）、大きい m は答えから作る、累乗の塔は直接計算と周期（39 万回） |
| `04_combinatorics.cpp` | 順列・集合の分割を全部並べる（スターリング数・ベル数・完全順列）、ナップサックの分割数、写像を全部並べた 12 相（n, k ≤ 6）、包除は 1 つずつ数える・Int128 で 2 つ 3 つの包除（20 万回） |
| `04_permutation.cpp` | next_permutation で全部並べる（順列・重複のある並び）、組合せは全部列挙、n = 33 の端、O(N^2) の mod 版（58 万回） |
| `05_linear_algebra.cpp` | 置換の和の行列式（有理数・mint・double・合成数の mod）、小行列式での階数、全域木・有向木は辺の組を全部試す、特性多項式・det(A + xB) は点での行列式（8 万回） |
| `05_xor_basis.cpp` | 部分集合の XOR を全部（基底のすべての操作）、区間ごとに作り直した基底、mod 2 の方程式は全部の割り当て（39 万回） |
| `06_linear_recurrence.cpp` | 漸化式どおりの計算（いろいろな mod・合成数・NTT の分岐）、ずらした列で巨大な K、フィボナッチの倍化公式、BM は答えの漸化式で列を作って（8,336 回。1 回が長さ数百の列） |
| `08_bigint.cpp` | Int128 の計算と比較、大きい数は q·b + r = a と 3 つの素数の mod での掛け算の検算（15 万桁の 999…9 の 2 乗も）、よく知られた値（2^1000・100! の各桁の和）（27 万回） |
| `10_game.cpp` | 全部の手を試す Nim・逆形 Nim、DAG は小さい方から mex、石取りは定義どおりの表と周期（19 万回） |
| `dynalist_items_phase4.txt` | Dynalist に貼る項目（13 本） |
| `dynalist_tool_index.txt` / `.opml` | 数学道具選択の索引（完成版。「作成予定」をなくし、Phase 4 の関数名を入れた） |

確認したこと（全 26 本）：
- g++ 13 と clang 18（-O2 -Wall -Wextra -Wshadow -pedantic）で C++17 / 20 / 23 とも警告なし・全テスト通過
- `-D_GLIBCXX_DEBUG -fsanitize=address,undefined` と、テストの量を減らさない `-O2 -fsanitize=undefined` の両方で全体（テスト・速度計測・デモ）が通過（ここで bell_table<long long> の三角形の溢れと、convolve_many がデバッグ版で極端に遅いことを見つけて直した）
- Main_base.cpp（ACL 込み）にグラフ・グリッド 46 本・幾何 17 本・数学 26 本の本体を 2 通りの順で全部貼り、gnu++2b / c++17・LOCAL あり / なしで g++・clang の両方でコンパイル・実行（名前の衝突なし。BigInt の abs・pow・gcd・isqrt を足しても abs(-3)・pow(2.0, 10)・gcd(12LL, 18LL)・isqrt(99) は今までどおり）。数学の各ファイルを 1 本ずつ Main_base に貼っても警告なし
- 同梱部分 69 個が持ち主のファイルと全部一致
- **独立したレビュー**（別のエージェントが本体だけを読み、自分で書いた愚直解・Python の多倍長整数と比べた。6 本と 7 本に分けた）：下の表のものを直した
- **代表的な問題を Main_base ＋ テンプレートの本体だけで解いて、全部の例が合うことを確かめた**（27 問：ACL Practice F、ATC001 C、ABC196 F、ABC230 Ex、ABC212 H、ABC206 E、ABC162 E、ABC300 G、ABC222 G、ABC335 G、ABC228 E、ABC131 C、ABC242 F、ABC297 F、ABC190 F、ABC276 C、ABC202 D、ABC253 Ex、ABC323 G、ABC141 F、ABC223 H、ABC283 G、ABC300 Ex、ABC169 B、ABC206 F、ABC297 G、ABC368 F。Phase 3 の 24 問も合わせて 51 問）。デモにはこのほか ABC291 G の例も入れた

### 独立したレビュー・サニタイザで直したもの

| ファイル | 内容 | 直し方 |
|---|---|---|
| 05_linear_algebra | `Matrix<long long>` を渡すと割り算が切り捨てになり、行列式・階数・全域木の数を黙って間違えた | 整数型の行列はコンパイルエラーにした（mint・double を使う。合成数の mod は det_mod_any） |
| 04_combinatorics | long long の表（ベル数・スターリング数・完全順列・分割数）と inclusion_exclusion が黙って溢れた／ベルの三角形が B_N より先の B_{N+1} まで作って溢れた | 溢れたら止める（符号が交互の和は Int128 で足して最後に確かめる）。三角形は B_N に要る行まで |
| 04_combinatorics | mod ≤ 16 の DynModInt だと小さい n でも止まった（内部の階乗の表を最初に 16 まで作っていた） | 表は必要な分だけ作る |
| 04_combinatorics | 玉も箱も区別しない balls_into_boxes が上限の n = 2×10^6 で約 10 秒（五角数定理 O(n√n)） | `partition_row<mint>(N)`（級数の逆数で O(N log N)）を足して使う。上限 10^7 |
| 02_divisor_transform | count_pairs_by_gcd の m(m − 1)/2 が、答えは収まるのに途中で溢れた | 偶数の方を先に 2 で割り、溢れたら止める |
| 02_divisor_transform | sum_gcd_tuples が K = 10^8 で約 26 秒・0.9 GB／len = 0 で 0 を返さなかった | ⌊K/d⌋ が同じ d をまとめ、Σφ を杜教篩で（O(K^{2/3})。K ≤ 10^13）／0 を返す |
| 07_convolution | convolve_i128 / convolve_ll の短い列（愚直の経路）で、上限を確かめずに Int128 が溢れた | 上限の見積もりを愚直の前に |
| 07_convolution | convolve_many が priority_queue を使っていて、_GLIBCXX_DEBUG では 1 回ごとにヒープ全体を確かめるので 10^5 個で数十分かかった | 短い順に並べて隣どうしを掛ける形にした（Release の速さは同じ） |
| 07_fps | fps_sqrt が mod が奇素数かを確かめず、合成数の mod で無限ループ・誤答 | 確かめて止める |
| 07_fps | fps_exp・fps_pow が小さい素数の mod（項の数 < mod）でも、内部で 2 倍の長さの log を取って止まった | 最後の段は n 項まで |
| 07_fps | fps_coeff_of_fraction が次数 10^5・K = 10^18 で約 2.5 秒 | NTT が使える mod では Q(−x) を値の入れ替えで作る（約 1.1 秒） |
| 08_bigint | 長さが大きく違う掛け算（10^5 桁 × 400 桁）が、短い方を 0 で埋めてカラツバにして約 130 倍遅かった | 長い方を短い方の長さで区切って掛ける |
| 08_bigint | gcd(BigInt, int) がコンパイルできない（std::gcd が選ばれる）・後置の ++ −− がない／isqrt（10^5 桁）が 2.5 秒 | 整数との gcd と後置 ++ −− を足した／上半分の桁の平方根から始めるニュートン法にした（2×10^5 桁で約 1 秒） |
| そのほか | surjections(n, k) が k > n でも k! を作って止まった、balls_into_boxes が k が巨大（10^9）だと階乗の表で止まった、表の関数が N < 0 を受け付けた、count_smooth が大きい N・P でメモリを使い切った、phi_sum・mertens に N の上限がなかった | 0 を返す／n 個の積で計算／止める／途中の数が 4×10^7 個を超えたら止める／N ≤ 10^13 で止める |
| 説明 | 「long long は溢れたら止める」の表で、途中の値が溢れて止まる場合（ゼータ変換など）・xor_convolution は mod が奇数・fps_sqrt は mod が奇素数・count_smooth の P の目安・10^9+7 では Bostan–Mori が約 5 倍遅い、を書いていなかった | ヘッダーに書いた |

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）

| 処理 | 時間 |
|---|---|
| convolve<ModInt998>（5×10^5 + 5×10^5）/ convolve<ModInt107>（3 つの素数）/ convolve_i128（2×10^5）/ convolve_many（(1 + a x) を 10^5 個） | 64 / 236 / 109 / 110 ms |
| n = 5×10^5 の fps_inv / log / exp / pow（k = 10^18）/ sqrt / fps_coeff_of_fraction（次数 10^5、K = 10^18） | 111 / 168 / 437 / 620 / 307 / 1100 ms |
| xor_convolution + or_convolution（2^20）/ subset_convolution（2^18） | 140 / 368 ms |
| 約数・倍数のゼータ・メビウス 4 回（N = 10^7）/ phi_sum(10^11) / mertens(10^11) | 534 / 1620 / 1609 ms |
| prime_count(10^11) / prime_count(10^12) / prime_sum(10^11) / count_smooth(10^16, 97) | 290 / 1482 / 309 / 903 ms |
| MultOrder（m ≈ 10^13）2×10^5 回 / mod_sqrt（p ≈ 10^13）10^5 回 / discrete_log（m ≈ 10^12） | 271 / 274 / 178 ms |
| stirling2_row(5×10^5) / stirling1_row(2×10^5) / partition_row(10^6) / count_divisible_by_any（10^18、素数 20 個） | 106 / 226 / 528 / 58 ms |
| N = 10^6 の inversion_count / perm_rank_mod / perm_unrank / multiset_perm_rank_mod | 181 / 53 / 113 / 30 ms |
| n = 500 の mint の mat_det / mat_inverse / det_affine、double の solve_linear、det_mod_any（n = 200、m = 10^18） | 49 / 231 / 425 / 36 / 29 ms |
| XorBasis（insert 10^6 回 + クエリ 2×10^6 回）/ PrefixXorBasis（4×10^5 個 + クエリ 2×10^5）/ Gf2System 2000 × 2000 | 93 / 225 / 80 ms |
| linear_recurrence_kth（d = 5×10^4、K = 10^18）/ 同じく 10^9+7・d = 10^4 / linear_recurrence_terms（d = 5×10^4、n = 10^6）/ berlekamp_massey（4000 項） | 578 / 649 / 530 / 12 ms |
| BigInt 10^5 桁 × 10^5 桁 / 2×10^5 桁 ÷ 5×10^4 桁 / isqrt（2×10^5 桁）/ 10000! | 159 / 266 / 954 / 132 ms |
| grundy_dag（2×10^5 状態・2.2×10^6 辺）/ subtraction_game_grundy（3×10^6、|S| = 10） | 19 / 28 ms |

### 実装で決めたこと

- **依存**：07_subset_transform は modint を同梱する（数え上げで mint を使うことが多いので）。02_prime_count は篩を自分で持つので基盤だけを同梱する（設計では sieve）
- **NTT** は 4 進のバタフライ（ACL と同じ形の表）。表は mod ごとに作り、DynModInt の mod が変わったら作り直す。NTT が使えない mod（10^9+7・合成数・偶数）は 3 つの素数で計算してガーナーで戻す（自動）。convolve_i128 は結果の大きさの見積もりで 3 つか 5 つの素数を選ぶ。ACL の convolution とは名前が違うので一緒に貼れる
- **Bostan–Mori**（fps_coeff_of_fraction・linear_recurrence_kth）は、NTT の値の並びで j と j^1 が ±x の組になることを使い、Q(−x) を値の入れ替えで作る（1 段で NTT 4 回）
- **線形代数**の T は「割り算ができる型」。整数型は static_assert で止める。ABC323 G のために `det_affine(A, B)`（det(A + xB)）を足した（B が正則でなくてよい）。有向の全域木は `count_arborescences(n, us, vs, root)`（root から向かう向き）
- **数え上げ**：`inclusion_exclusion_count<mint>(n, g)`（対称な包除。ABC242 F）と `partition_row<mint>(N)` を足した。_table<T> は mint でも long long でも使え、整数型は溢れたら止める
- **順列**：列は vector でも string でもよい（値は比べられれば何でも。1〜n でも 0〜n−1 でも同じ番号）。厳密な番号は Int128（収まらなければ止める）、大きいときは _mod
- **XOR の基底**は簡約した形で持ち、先頭ビットの集合を 1 つの整数で持って insert・reduce・max_xor を速くした。x 以下の作れる値の個数 `count_le` を足した
- **BigInt** は 10^9 進。割り算は C++ と同じく 0 の方向に切り捨て。掛け算はカラツバで、途中の値が long long に収まるように 10^6 進（長さ 2^14 まで）か 10^3 進（2^23 まで）に直す。÷ は Knuth のアルゴリズム D
- **ゲーム**：grundy_dag は再帰しない（深い DAG でもスタックが溢れない）。`SubtractionGame` は直前 max(S) 個の Grundy 数の窓が同じになったところで周期を決める（周期の始まりはできるだけ前に）
- 設計になかったが足したもの：`count_le`（05_xor_basis）、`inclusion_exclusion_count`・`partition_row`（04_combinatorics）、`det_affine`（05_linear_algebra）、`perm_rank_mod`・`multiset_perm_rank_mod`（04_permutation）、`walsh_hadamard`（07_subset_transform）、`SubtractionGame`（10_game）

---

## 3. 依存関係

```
00_int_common（数学・幾何のすべてのファイルが同梱する）
00_math_base（数学のすべてのファイルが同梱する。00_int_common を同梱）
├─ 00_modint
├─ 01_digits
├─ 02_factorize
├─ 02_sieve              ← 02_factorize（Sieve S(N) が N を超える数を factorize に任せるため。Phase 3 で決定）
├─ 02_divisor_transform  ← 02_sieve（02_factorize も入る）
├─ 02_prime_count        （Phase 4 で決定：篩は自分で持つので基盤だけ）
├─ 03_mod_equation
├─ 03_discrete_log       ← 02_factorize
├─ 04_binomial           ← 00_modint
├─ 04_combinatorics      ← 00_modint, 04_binomial, 07_convolution
├─ 04_permutation        ← 00_modint, 04_binomial
├─ 05_matrix             ← 00_modint
├─ 05_linear_algebra     ← 00_modint, 05_matrix
├─ 05_xor_basis
├─ 06_quotient
├─ 06_series             ← 00_modint
├─ 06_linear_recurrence  ← 00_modint, 07_convolution
├─ 07_convolution        ← 00_modint
├─ 07_fps                ← 00_modint, 07_convolution
├─ 07_subset_transform  ← 00_modint（数え上げで mint を使うことが多いので。Phase 4 で決定）
├─ 08_fraction
├─ 08_bigint
├─ 09_bit
└─ 10_game
```

- ← の右に書いたファイル（と基盤）を同梱する。例：`04_binomial` には math_base と modint が入る（`Binomial<mint>` を 1 回のコピーで使える）
- `05_matrix` は modint を同梱する（`Matrix<mint>` がいちばんよく使う形なので）
- グラフ・グリッド・幾何には依存しない（同梱しない）。組み合わせ方は索引（4 章「グラフ・グリッド／幾何と組み合わせる」）と使い方関数のコメントに書く

### 作る順番

- Phase 2（基盤）：00_int_common、00_math_base、00_modint（＋ 幾何 16 本の同梱部分の更新と再テスト、update_bundles のコピー、Dynalist の項目）
- Phase 3（★）：01_digits、02_sieve、02_factorize、03_mod_equation、04_binomial、05_matrix、06_quotient、06_series、08_fraction、09_bit（10 本）
- Phase 4（☆）：07_convolution（ほかの ☆ が使うので最初）、07_fps、07_subset_transform、02_divisor_transform、02_prime_count、03_discrete_log、04_combinatorics、04_permutation、05_linear_algebra、05_xor_basis、06_linear_recurrence、08_bigint、10_game（13 本）
- 最終照合：旧ファイルを 1 機能ずつ、過去の解答、AtCoder の典型（問題文を確かめながら）と広い洗い出し、Dynalist の項目と突き合わせる

合計 26 ファイル（★ 13（00_int_common を含む）・☆ 13）。

---

## 4. 道具選択体系（案）

グラフ・グリッド・幾何と同じく、Dynalist の「00_典型発想インデックス」の中に 1 項目（▼ 数学道具選択（第二世代））として置く。辿り方は「問題の特徴（大分類）」→「聞かれていること（小分類）」→「ファイル名と関数」。各 Phase で関数が確定したら、`math_v2\dynalist_tool_index.txt` / `.opml` として完成させる。旧索引の「▼ 数学高速化」「▼ 重複を除いて数える」の各キーワード（mod inverse・CRT・φ関数・floor_sum・包除原理・Möbius）にも、ここの行を子として置く。

```
▼ 数学道具選択（第二世代）
	▼ 使い方
		1. 問題の特徴で大分類を選ぶ（mod／大きな数・割り算／桁・進数／素数・約数／gcd・lcm／合同式／累乗・位数／数え上げ／行列・漸化式／和／畳み込み／ビット／分数・多倍長／ゲーム）
		2. 聞かれていることで小分類を選ぶ
		3. 書いてあるファイルを貼る（依存するファイルは同梱済み。1回のコピーで動く）
		答えを mod で出すときは using mint = ModInt998;（10^9+7 なら ModInt107、入力で決まるなら DynModInt）を最初に書く
		範囲の引数は半開区間 [l, r)。「L 以上 R 以下」は _closed 版。k 番目は 0 から数える
	▼ 答えを mod で求める　00_modint
		998244353 / 10^9+7 で割った余り → using mint = ModInt998; ／ ModInt107
		mod が入力で与えられる（素数とは限らない） → DynModInt::set_mod(M)
		割り算・逆元（mod が合成数でも gcd = 1 なら） → a / b ／ a.inv()
		期待値・確率を mod で → 使い方デモ（ABC280 E）
		mod の答えから元の分数を推測する（デバッグ） → rational_reconstruct(x)
		64 ビットの mod の累乗・逆元 → 00_math_base　mod_pow(a, e, m) ／ mod_inv(a, m)
	▼ 大きな数・割り算・平方根　00_math_base
		負の数の切り捨て・切り上げ・余り → floor_div ／ ceil_div ／ mod_floor
		[l, r) の d の倍数の数 → count_multiples(l, r, d)（L 以上 R 以下なら _closed）
		平方根・k 乗根を誤差なく → isqrt ／ kth_root_floor(n, k)
		a^b が上限を超えるか・溢れずに累乗 → pow_cap(a, b, 上限)
		掛け算が溢れるか → mul_overflows(a, b) ／ mul_cap
		lcm が大きくなる → lcm_cap(a, b, 上限)
		Int128 の入出力 → cin >> x ／ cout << x
		小数の入力を誤差なく → parse_fixed(s, k) ／ read_fixed(k)
		分数の大小を誤差なく → frac_compare(a, b, c, d)（何度も使うなら 08_fraction）
		等差数列・等比数列の和（mod は何でも） → arith_sum ／ geom_sum_mod
	▼ 桁・進数・回文数　01_digits
		桁和・桁数・桁の配列 → digit_sum ／ num_digits ／ digits_of
		n 進数に直す・n 進数から直す（−2 進も） → to_base ／ from_base
		n 進数として読んだ値が上限以下か → base_value_leq(s, base, limit)
		回文数の k 番目（0 から）・N 以下の個数・N 以上で最小の回文数 → kth_palindrome ／ count_palindromes_up_to ／ ceil_palindrome ／ floor_palindrome
		数をつなげた値・同じ数を K 回つなげた値を mod で → concat_mod ／ repeat_concat_mod
		大きな数（文字列）の余り・大小 → str_mod ／ str_num_compare
		部分文字列の値を mod で何度も → DigitPrefix
		1〜N の桁数の和・数字 d の出現回数 → total_digits_up_to ／ count_digit_up_to
	▼ 素数・素因数分解・約数
		N 以下（N ≤ 5×10^7）の数を何度も分解・素数判定 → 02_sieve　Sieve S(N)
		1 つの数（10^18 まで）の素数判定・分解・約数 → 02_factorize　is_prime ／ factorize ／ divisors_of
		N 以下の素数の一覧 → 02_sieve　primes_up_to(N)
		[L, R) の素数・素因数分解（R ≤ 10^14、幅 ≤ 10^7） → 02_sieve　segment_primes ／ segment_factorize
		φ・μ・約数の個数・約数の和の表 → 02_sieve　phi_table ／ mu_table ／ divisor_count_table ／ divisor_sum_table
		異なる素因数の個数の表 → 02_sieve　omega_table
		n! が p で何回割れるか・k が n! を割り切る最小の n → 02_factorize　factorial_exponent ／ min_factorial_multiple
		平方数を掛けて平方数にする・平方因子を除く → square_free_decompose（Sieve のメンバにも）
		a^b（b ≥ 2）の形の数の個数 → 02_factorize　count_perfect_powers
		N 以下の素数の個数・素数の和（N ≤ 10^12） → 02_prime_count　prime_count ／ prime_sum ／ PrimeCountTable（⌊N/i⌋ ごとの値）
		素因数が P 以下の数の個数 → 02_prime_count　count_smooth
	▼ gcd・lcm・互いに素
		多くの数の lcm を mod で（mod を取りながら lcm を取ると誤る） → 02_factorize　lcm_mod_of ／ Sieve::lcm_mod
		N 以下で M と互いに素な数 → 02_factorize　count_coprime_up_to
		gcd が g になる組の数・列の gcd の和 → 02_divisor_transform　count_pairs_by_gcd ／ sum_gcd_tuples
		Σφ(i)・Σμ(i)（N ≤ 10^11）・互いに素な組の数 → 02_divisor_transform　phi_sum ／ mertens ／ count_coprime_pairs
		約数・倍数についての和（ゼータ・メビウス変換）・gcd ／ lcm の畳み込み → 02_divisor_transform　divisor_zeta ／ multiple_mobius ／ gcd_convolution
		いずれかで割り切れる数の個数（包除） → 04_combinatorics　count_divisible_by_any（L 以上 R 以下なら _closed）
	▼ 合同式・方程式　03_mod_equation
		ax ≡ b (mod m) → linear_congruence
		ax + by = c の整数解・範囲の中の解の数・x が最小の解 → linear_diophantine ／ diophantine_count ／ diophantine_min_x
		x ≡ r_i (mod m_i) をまとめる（mod が互いに素でなくても） → crt_merge ／ crt_all
		互いに素な mod での答えを別の mod で → garner
		周期的に繰り返す 2 つの区間が重なる最初の時刻 → 使い方デモ（ABC193 E）
	▼ 累乗・位数・離散対数　03_discrete_log
		a^k ≡ 1 となる最小の k（位数）・同じ mod で何度も → mult_order ／ MultOrder
		原始根・カーマイケル関数 λ(n) → primitive_root ／ carmichael_lambda
		a^x ≡ b となる x → discrete_log
		mod p の平方根 → mod_sqrt
		指数が巨大（a^(b^c)） → pow_tower_mod（素数 p なら指数を p − 1 で割った余りに。a が p の倍数のときに注意）
	▼ 数え上げ
		nCr・nPr・nHr を何度も → 04_binomial　Binomial<mint> B(N)
		n が大きく r が小さい nCr → 04_binomial　binom_large_n
		mod が合成数・mod なしの厳密な nCr → 04_binomial　pascal_table ／ binom_exact ／ binom_mod_any
		mod が小さい素数 → 04_binomial　binom_lucas
		格子路・境界を越えない格子路（カタラン） → 04_binomial　paths ／ paths_not_crossing
		スターリング数・ベル数・分割数・完全順列・全射 → 04_combinatorics　stirling2 ／ stirling1_row ／ bell ／ partition_row ／ derangement_table ／ surjections
		写像 12 相（玉と箱。区別する・しない、1 個以下・1 個以上） → 04_combinatorics　balls_into_boxes
		包除原理 → 04_combinatorics　inclusion_exclusion ／ inclusion_exclusion_count（使い方デモ：ABC242 F・ABC297 F）
		順列・重複のある並びの何番目か・k 番目 → 04_permutation　perm_rank ／ perm_unrank ／ multiset_perm_rank ／ multiset_perm_unrank
		k 個選ぶ組合せの k 番目・何番目か → 04_permutation　kth_combination ／ combination_rank
		転倒数 → 04_permutation　inversion_count
	▼ 行列・漸化式
		行列の累乗（項の少ない漸化式・K 歩の数え上げ） → 05_matrix　Matrix<mint>::pow
		ちょうど K 歩の最短路・K 歩で行けるか・XOR と AND の漸化式 → 05_matrix　MinPlus ／ OrAnd ／ XorAnd
		項の少ない線形漸化式の第 N 項 → 05_matrix　recurrence_nth_matrix(c, a, n)
		2×2・3×3 の行列を何度も・行列式・逆行列 → 05_matrix　MatrixFixed<T, K>
		連立一次方程式・行列式・逆行列・階数 → 05_linear_algebra　solve_linear ／ mat_det ／ mat_inverse ／ mat_rank
		合成数の mod の行列式 → 05_linear_algebra　det_mod_any
		特性多項式・det(A + xB) → 05_linear_algebra　char_poly ／ det_affine（使い方デモ：ABC323 G）
		全域木の数・有向木の数（行列木定理） → 05_linear_algebra　count_spanning_trees ／ count_arborescences
		XOR で作れる値（最大・k 番目・作れるか） → 05_xor_basis　XorBasis
		区間の要素の XOR で作れるか → 05_xor_basis　PrefixXorBasis
		mod 2 の連立方程式 → 05_xor_basis　Gf2System
		項の多い線形漸化式の第 K 項 → 06_linear_recurrence　linear_recurrence_kth
		最初の数項から漸化式を推測 → 06_linear_recurrence　berlekamp_massey ／ guess_kth
	▼ 和
		⌊N/i⌋ の和・⌊N/i⌋ が同じ i をまとめる → 06_quotient　quotient_blocks
		Σ i⌊N/i⌋・Σ (N mod i) → 06_quotient　sum_i_floor_div ／ sum_mod_i
		Σ ⌊(a i + b)/m⌋・直線の下の格子点 → 06_quotient　linear_floor_sum ／ count_lattice_under_line
		べき乗和 Σ i^k → 06_series　sum_powers<mint>
		Σ i·r^i → 06_series　arith_geom_sum
		多項式の値から別の点の値（ラグランジュ補間） → 06_series　lagrange_interpolate
		多項式の割り算（次数が小さい） → 06_series　poly_divmod
	▼ 畳み込み・多項式
		2 つの列の畳み込み（mod は何でも・整数のまま） → 07_convolution　convolve ／ convolve_ll ／ convolve_mod
		ずらして重ねた積の和・一致の数 → 07_convolution　correlate
		多くの多項式の積 → 07_convolution　convolve_many
		形式的冪級数（逆数・log・exp・累乗・平方根） → 07_fps　fps_inv ／ fps_log ／ fps_exp ／ fps_pow ／ fps_sqrt
		f(x + c)・多項式の割り算 → 07_fps　taylor_shift ／ fps_divmod
		P(x)/Q(x) の K 次の係数 → 07_fps　fps_coeff_of_fraction
		OR・AND・XOR の畳み込み・部分集合の和 → 07_subset_transform　or_convolution ／ and_convolution ／ xor_convolution ／ subset_zeta
		集合を 2 つに分ける数え上げ → 07_subset_transform　subset_convolution
	▼ ビット　09_bit
		部分集合の列挙・k 個選ぶ・必ず含む／含まないビット → for_each_submask ／ for_each_subset_k ／ for_each_mask_between
		立っているビットの位置・数 → bit_indices ／ for_each_bit ／ bit_count ／ top_bit
		グレイコード → gray_code
		0〜N の各ビットの 1 の個数・L ⊕ … ⊕ R → count_bit_up_to ／ popcount_sum_up_to ／ xor_range_closed
		全部の組・全部の区間の XOR の和 → xor_pair_sum ／ subarray_xor_sum
		bit 全探索 → 使い方デモ（ABC197 C）
	▼ 分数・多倍長
		分数を誤差なく（比で並べる・既約分数・小数で表示） → 08_fraction　Fraction ／ to_decimal
		分母が小さい近い分数 → 08_fraction　best_rational_approx
		long long・Int128 にも収まらない整数 → 08_bigint　BigInt
	▼ ゲーム　10_game
		Grundy 数（状態の DAG）・mex・Nim・逆形 Nim → grundy_dag ／ mex ／ nim_first_wins ／ misere_nim_first_wins
		L〜R 個取れる石取り・取れる数が決まった石取り（周期も） → range_take_grundy ／ subtraction_game_grundy ／ SubtractionGame
	▼ グラフ・グリッド／幾何と組み合わせる
		K 歩で行ける道の数 → graph_grid_v2 の隣接を行列にして 05_matrix　Matrix<mint>::pow
		ちょうど K 歩の最短路 → 05_matrix　Matrix<MinPlus<long long>>
		全域木の数 → graph_grid_v2 の辺（G.es）を us・vs にして 05_linear_algebra　count_spanning_trees へ
		最短路の本数を mod で → graph_grid_v2 02_shortest_path_dag　count_shortest_paths(G, r, mod)
		DAG の経路の数 → graph_grid_v2 04_topo_sort　dag_count_paths(G, s, mod)
		向き・傾きを gcd で正規化 → geometry_v2 00_geo_base　direction_key
```

---

## 5. 既存テンプレートとの照合（Phase 1-6）

### 5-0. 対象

Dynalist（dynalist-2026-10-7.opml）の「08_数学・数論」と、「14_テクニック・発想」「15_道具」に載っている数学関係の項目が指す旧ファイル（`C:\atcoder_template\典型90\ソース`）。

| 旧ファイル | Dynalist の項目 |
|---|---|
| PrimeFactorizeTemplateTest.txt / SPFTemplateTest.txt | 素因数分解系テンプレート（2 本が同じ名前） |
| 202D.txt / 425E.txt | パスカルの三角形（mod なし / mod あり）、modint（425E） |
| 368D.txt（368D_2.txt） | nCr・nHr・nPr・n が小さい場合の nCr・整数の累乗・累乗(mod)・逆元(mod)・逆元(非素数 mod) |
| 445E.txt | safe_gcd / safe_lcm / safe_lcm_mod |
| NumberTheoryTemplateTest.txt | 拡張ユークリッド・mod 逆元（合成数）・一次不定方程式・ax ≡ b (mod m)・CRT・φ 関数 |
| ARC141A.txt | 部分文字列の数値化・数値の結合 |
| MathSumTest.txt / 448E.txt | 数列の和・連結系テンプレ（非素数 mod 対応） |
| FastModIntTemplate.txt / DynModIntTemplate.txt / AnyModIntTemplate.txt / MyModintTemplateTest.txt | FastModInt・DynModInt・AnyModInt・FastModInt（定数 mod）・DynBarrettModInt |
| BigIntTemplateTest.txt | BigInt（多倍長整数） |
| FractionTemplateTest.txt | Fraction（有理数構造体） |
| CSR001J.txt | 転倒数を高速で取得 |
| PalindromeTemplateTest.txt / 363D.txt | Palindrome Library（回文数） |
| Base10TemplateTest.txt | n進数テンプレ（最小構成） |
| floor_sumTest.txt | ACL floor_sum 使用テンプレ |
| Mat2TemplateTest.txt | 2×2 行列テンプレート |
| LinearRecurrenceTemplate.txt | N項間線形漸化数列の第k項(mod) |
| 368F.txt | Grundy数(NIM)（14・15） |
| 238D.txt / BitSetGenerator.txt | 1 が立っている bit の index を列挙 / bitset ジェネレーター（15） |
| ARC59C.txt / 126C.txt | 繰り下がり、繰り上がり / 出力の桁数（14） |
| 114C.txt / 161D.txt | 状態遷移 DFS（N 以下の個数）/ 状態遷移 BFS（K 番目の値）（14） |

旧索引（00_典型発想インデックス）の数学関係：「▼ 数学高速化」（mod inverse・CRT・φ関数・floor_sum）、「▼ 重複を除いて数える」（包除原理・Möbius）、「▼ 状態遷移を考える」（bitDP・桁DP）。どのキーワードにも子項目がなく、テンプレートへのリンクになっていない。

### 5-1. 新設計にあって、既存にないもの

- 64 ビットの素数判定（ミラー・ラビン）・素因数分解（ポラード・ロー）、区間篩、累乗数の数、n! の中の素数の数、k が n! を割り切る最小の n
- 負の数の floor_div / ceil_div / mod_floor、溢れない累乗・掛け算（`_cap`）、k 乗根、Int128 の入出力、範囲の中の倍数の数
- 一般の N×N 行列と半環（最短路・XOR）、ガウスの消去法・行列式・行列木定理・特性多項式、XOR 基底、mod 2 の連立方程式
- ⌊N/i⌋ のブロックの列挙、ACL を使わない floor_sum、直線の下の格子点の数
- 二項係数の n が大きい版・リュカ・合成数の mod・格子路・鏡像法、スターリング数・ベル数・分割数・完全順列・全射・写像 12 相、包除原理、約数・倍数のゼータ／メビウス変換、gcd / lcm 畳み込み、φ・μ の和
- 順列・重複のある並び・組合せの番号と k 番目
- 位数・原始根・離散対数・mod の平方根・累乗の塔
- 任意の mod の畳み込み・FPS・Berlekamp–Massey、ラグランジュ補間・べき乗和、OR / AND / XOR 畳み込み・部分集合畳み込み
- ビットの補助（supermask・Gray コード・0〜N のビットの数・XOR の和）、Grundy の補助（mex・Nim・石取り）
- 分数の最良近似・小数の文字列から分数、`rational_reconstruct`

### 5-2. 既存にあって、新設計になかったもの → 追加した

| 追加先 | 追加した機能 | 参考にした既存 |
|---|---|---|
| 02_sieve | kth_prime（k 番目の素数）・count_primes(l, r)・largest_prime_factor・big_omega（重複も数える）・is_square_free・square_free_decompose（x = a²·b） | SPFTemplateTest.txt（kthPrime・range_prime_count・lpf・OmegaCount・isSquareFree・sqrt_decompose_fast） |
| 02_divisor_transform | count_coprime_pairs(N)（1 ≤ i, j ≤ N で互いに素な組の数） | SPFTemplateTest.txt（count_coprime_pairs） |
| 02_sieve / 02_factorize | 多くの数の lcm を mod で（lcm_mod / lcm_mod_of） | 445E.txt（safe_lcm_mod の「複数要素の lcm を取る場合は mod を取り続けると値がずれる」） |
| 00_modint | raw・explicit operator bool・is_unit / try_inv・ハッシュ・mod が偶数でも動く | FastModInt / AnyModInt |
| 01_digits | 回文数の一式（k 番目・番号・次・前・N 以下の個数・半分から作る・文字列版） | PalindromeTemplateTest.txt / 363D.txt |
| 01_digits | DigitPrefix（部分文字列の値）・concat | ARC141A.txt |
| 01_digits | base_value_leq（n 進数として読んだ値と上限の比較） | Base10TemplateTest.txt（leq_base） |
| 05_matrix | MatrixFixed（2×2 などを速く。det・inv も） | Mat2TemplateTest.txt |
| 08_fraction | to_decimal(k)・reciprocal | FractionTemplateTest.txt |
| 09_bit | for_each_mask_between（必ず 1・必ず 0 のビットを決めた列挙。個数 k の指定も）・for_each_subset_k | BitSetGenerator.txt（ForcedOne / ForcedZero / CombinationBits） |
| 04_binomial | pascal_table（合成数の mod・厳密な値） | 202D.txt / 425E.txt |
| 10_game | grundy_dag | 368F.txt（Grundy の add_edge / calc） |

### 5-3. 既存の不具合・使いにくさ（新設計で解消）

不具合の詳細は `Desktop\claude\atcoder_template_review_20260925.md` を参照。

| 既存 | 指摘 | 新設計での解消 |
|---|---|---|
| MathSumTest.txt（中-20） | 等差数列の和が __int128 に入る前に溢れる。a・d が負だと負の値を返す。等比数列も a < 0 で負 | Int128 で計算し、mod の結果は 0 以上にそろえる（負の値・大きい値をテストする） |
| Base10TemplateTest.txt（中-21） | 11 進以上に非対応、負の数で空文字列、digits(0) が空 | 0〜9・A〜Z、0 と負の数、負の進数 |
| NumberTheoryTemplateTest.txt（低） | safe_crt2 が lcm > 2^63 を検出しない | crt_all は `overflow` で知らせる |
| PrimeFactorizeTemplateTest.txt（3-1） | 試し割りだけで 10^18 に使えない | ミラー・ラビン・ポラード・ロー |
| Mat2TemplateTest.txt（3-1） | 2×2 だけ。inv のコメントが問題固有 | Matrix<T>（N×N・半環）と MatrixFixed |
| ARC59C.txt（3-1） | 切り上げを (sum + N − 1) / N で書いていて負の数で誤る。ヘッダにも floor_div / ceil_div がない | ceil_div / floor_div |
| floor_sumTest.txt（3-2） | ACL が必要。⌊N/i⌋ のブロックを列挙する関数がない | 自作の linear_floor_sum・quotient_blocks |
| LinearRecurrenceTemplate.txt（7-2） | ACL の convolution（NTT が使える mod だけ）が必要で、10^9+7 で使えない | 07_convolution の convolve（mod は何でも） |
| FractionTemplateTest.txt | 足し算・掛け算が long long で溢れても気づかない（比較は __int128 で正しい） | 途中を 128 ビットで計算し、収まらなければ止める。Int128 版も |
| 368D.txt / 445E.txt | mod_inv がフェルマー（mod が素数のとき）だけ。safe_lcm_mod は 2 数だけ | mod_inv は合成数の mod でも。lcm_mod_of は多くの数 |
| DynModIntTemplate.txt | mod は奇数かつ 2^30 以下 | DynModInt は 1 ≤ mod < 2^31（偶数も） |
| modint 4 種（Fast / Dyn / Any / MyModint） | 名前が重なる（FastModInt が 2 つ。見直し一覧 6-3）、どれを使うか迷う | ModInt・DynModInt の 2 種類に整理 |
| 368D.txt・425E.txt など | 問題の解答の中にあり、`using mint`・MOD・M などのグローバルな名前がぶつかる（6-3） | テンプレートにして名前を決める |
| Dynalist（6-8） | 素因数分解系テンプレートが同じ名前で 2 つ、368D の 8 項目が同じファイルを指す、「modint」が問題の解答（425E）を指す、索引の「包除原理」「Möbius」「φ関数」などにリンクがない | 新しい索引で整理する（4 章） |
| 見直し一覧 7-1・7-2 | 包除原理・Möbius のテンプレートがない。N×N 行列・ミラー・ラビン・floor_div・Berlekamp–Massey・任意 mod の畳み込みがない | 04_combinatorics・02_divisor_transform・05_matrix・02_factorize・00_math_base・06_linear_recurrence・07_convolution |

### 5-4. 体系の外に置くもの（今のまま残す）

- CSR001J.txt の BIT：データ構造。転倒数は 04_permutation の中に BIT を持つ
- RangeGCDQuerySegTree.txt：データ構造
- 114C.txt / 161D.txt（使える数字が決まった数を DFS / BFS で数える・列挙する）：探索の型として 14_テクニックに残す（多くは 01_digits の numbers_from_digits で書ける）
- 126C.txt（出力の桁数）：出力のしかた（setprecision）。使い方関数で示す
- DigitDPTemplateTest.txt（桁 DP）：DP
- 区間・二分探索・累積和：それぞれのテンプレート

### 5-5. AtCoder の典型との照合

グラフ・グリッドの最終照合で足りないものが見つかった反省から（幾何と同じく）、Phase 1 の時点で典型と突き合わせた。代表例は Phase 3・4 と最終照合で問題文を確かめ直す。

| 典型 | 代表例 | 担当 |
|---|---|---|
| 溢れの検出・巨大な累乗の比較 | ABC169 B、典型90 020 | 00_math_base pow_cap / mul_overflows |
| 小数の入力を誤差なく | ABC169 C | 00_math_base parse_fixed |
| lcm が大きくなる | 典型90 038 | 00_math_base lcm_cap |
| gcd | 典型90 022 | std::gcd、00_math_base gcd_of |
| 3 数の積が N 以下（立方根で範囲を決める） | ABC227 C | 00_math_base kth_root_floor |
| 切り上げ | ARC059 C | 00_math_base ceil_div |
| mod の累乗・逆元 | 典型90 069 | 00_modint |
| 期待値を mod で | ABC280 E、ABC263 E | 00_modint（△） |
| 桁和 | ABC083 B | 01_digits digit_sum |
| 進数変換（負の進数・8 進など） | ABC105 C、ABC186 C、ABC336 C、典型90 067 | 01_digits to_base / from_base |
| n 進数として読んだ値が上限以下 | ABC192 D | 01_digits base_value_leq |
| 1〜N の桁数の和 | 典型90 082 | 01_digits total_digits_up_to |
| 同じ数を K 回つなげた値・数をつなげる | ABC357 D、ARC141 A | 01_digits repeat_concat_mod / DigitPrefix |
| 回文数 | ABC363 D | 01_digits kth_palindrome |
| 素数の列挙 | ABC300 D、ABC383 D、ABC400 E | 02_sieve primes_up_to / Sieve |
| 篩で全部の数を分解・素因数の数の表 | ABC177 E、典型90 030 | 02_sieve Sieve / omega_table |
| 区間篩 | ABC412 E | 02_sieve segment_factorize |
| 素因数分解・約数 | ABC169 D、ABC180 C、ABC190 D、ABC215 D、典型90 075・085 | 02_factorize / 02_sieve |
| n! の中の素数の数 | ABC052 C、ABC280 D | 02_factorize factorial_exponent / min_factorial_multiple |
| 平方因子 | ABC254 D | 02_sieve square_free_decompose |
| 9×10^18 の数の分解（p²q の形） | ABC284 D | 02_factorize factorize |
| a^b の形の数の数 | ABC361 F | 02_factorize count_perfect_powers |
| 約数の個数の和 | ABC172 D | 06_quotient sum_i_floor_div（02_sieve divisor_count_table でも） |
| gcd の数え上げ（倍数のメビウス） | ABC162 E、ABC206 E | 02_divisor_transform |
| 素因数が P 以下の数 | ABC300 G | 02_prime_count count_smooth |
| 一次合同式 | ABC186 E | 03_mod_equation linear_congruence |
| 中国剰余定理 | ABC193 E | 03_mod_equation crt_merge（△ の型も） |
| 位数・離散対数 | ABC222 G、ABC335 G | 03_discrete_log mult_order / primitive_root |
| 指数が巨大な累乗（フェルマー） | ABC228 E | 03_discrete_log pow_tower_mod、00_math_base mod_pow |
| 二項係数 | ABC132 D、ABC145 D、ABC167 E、ABC262 E | 04_binomial Binomial |
| n が大きく r が小さい二項係数 | ABC156 D | 04_binomial binom_large_n |
| 重複組合せ | ABC110 D | 04_binomial H |
| 厳密な二項係数で k 番目の並び | ABC202 D | 04_binomial binom_exact、04_permutation multiset_perm_unrank |
| 鏡像法 | ABC205 E | 04_binomial paths_not_crossing |
| 合成数の mod のパスカルの三角形 | ABC425 E | 04_binomial pascal_table |
| リュカの定理（mod 3） | ARC117 C | 04_binomial binom_lucas |
| 包除原理 | ABC131 C、ABC242 F、ABC297 F | 04_combinatorics |
| 前の順列・順列の番号 | ABC276 C | 04_permutation |
| 転倒数 | ABC190 F | 04_permutation inversion_count |
| 行列の累乗 | EDPC R、ABC009 D（XOR と AND）、ABC293 E | 05_matrix |
| 行列木定理・行列式 | ABC253 Ex、ABC323 G | 05_linear_algebra |
| XOR 基底 | ABC141 F、ABC223 H、ABC283 G | 05_xor_basis |
| ⌊N/i⌋ のブロック | ABC044 D、ABC230 E | 06_quotient quotient_blocks |
| floor_sum・直線の下の格子点 | ABC372 G、ACL Practice C | 06_quotient linear_floor_sum / count_lattice_under_line |
| 等比数列の和（任意の mod） | ABC293 E | 00_math_base geom_sum_mod |
| ラグランジュ補間・べき乗和 | ABC208 F、ARC033 D | 06_series |
| 多項式の割り算 | ABC245 D | 06_series poly_divmod |
| 線形漸化式の第 K 項 | ABC300 Ex | 06_linear_recurrence |
| 畳み込み | ABC196 F、ABC291 G | 07_convolution |
| FPS | ABC230 Ex | 07_fps |
| XOR 畳み込み | ABC212 H | 07_subset_transform |
| 分数の比較 | ABC308 C | 08_fraction |
| ビット | ABC121 D、ABC147 D、ABC197 C、ABC238 D、ABC269 C、ABC356 D、ABC365 E | 09_bit |
| Grundy 数・石取り | ABC206 F、ABC297 G、ABC368 F | 10_game |

### 5-6. 追加しないもの（理由）

| 機能 | 理由 |
|---|---|
| 多点評価・多項式補間（O(N log² N)）・FPS の合成・逆関数 | ABC では G / Ex でもまれ。必要なら個別に書く |
| 乗法的関数の和の一般の方法（min_25 篩など） | 同上（素数の個数・和と φ・μ の和は入れる） |
| mod の k 乗根（k ≥ 3） | 同上（平方根は入れる） |
| 2^63 以上の数の素数判定・素因数分解 | 入力が long long の範囲 |
| 行列の高速な累乗（ケーリー・ハミルトン） | 06_linear_recurrence（特性多項式＋Bostan–Mori）で代わりになる |
| 確率・期待値の専用テンプレート | △（modint の使い方） |
| 数値計算（ニュートン法・数値積分） | 幾何・二分探索の範囲 |

---

## 6. 未決事項と次にやること

### 6-1. ユーザーへの確認（2026-10-10 回答済み）

| 確認したこと | 回答・決定 |
|---|---|
| 1. 一覧（2 章の 25 ファイル）と ★ / ☆ の分け方 | このままでよい。☆ も全部作る |
| 2. 幾何と共通の関数（isqrt / isqrt_ceil・Int128 / UInt128・int128_to_string・parse_fixed / read_fixed）の扱い | **案 A**：共通の基盤 `00_int_common`（ガード `GG_INT_COMMON`）に切り出し、数学・幾何の両方が同梱する。幾何の 16 本の同梱部分も更新して再テストする |
| 3. modint の名前と使い方 | `ModInt<MOD>`・`ModInt998`・`ModInt107`・`DynModInt`。`mint` はテンプレートで定義せず、使う側で `using mint = ModInt998;` と書く |
| 4. 範囲の引数 | 半開区間 [l, r)。閉区間は `_closed` の名前で用意する（グラフ・幾何と同じ） |
| 個別には聞かなかったもの | ACL を使わない、整数の型・溢れの方針（1-2）、名前の付け方。設計書の案のまま進める |

### 6-2. 次にやること

1. ~~6-1 の確認~~（済み）
2. ~~Phase 2：`00_int_common.cpp`・`00_math_base.cpp`・`00_modint.cpp`、幾何 16 本の同梱部分の更新・再テスト、update_bundles のコピー、`dynalist_items_phase2.txt`~~（済み。2-P2 章）。ユーザーが内容を確認し、Dynalist に `dynalist_items_phase2.txt` の項目を入れる
3. ~~Phase 3：★ の 10 ファイル（終わりに独立したレビューと、代表例の問題を解く確認）~~（済み。2-P3 章）。ユーザーが内容を確認し、Dynalist に `dynalist_items_phase3.txt` の項目と `dynalist_tool_index` を入れる
4. ~~Phase 4：☆ の 13 ファイル（同じ）~~（済み。2-P4 章）。ユーザーが内容を確認し、Dynalist に `dynalist_items_phase4.txt` の項目と `dynalist_tool_index`（完成版）を入れる
5. ~~最終照合：旧ファイル（5-0）を 1 機能ずつ、過去の解答、典型（5-5 の問題文を確かめながら）と数学の問題の広い洗い出し、Dynalist の項目~~（済み。7 章）。ユーザーが内容を確認し、Dynalist の項目（7-7 の 19 項目）と索引（`dynalist_tool_index.opml`）を入れ替える
6. Dynalist への反映はユーザーの環境で行う
7. 追加・修正は今までの成果物（既存のファイル・索引・Dynalist の項目・この設計書）に入れる（幾何の最終照合と同じ。新しいファイルは作っていない）

---

## 7. 最終照合（2026-10-11）

ユーザーの指示：「最終照合をお願いします。」（6-2 の 5：旧ファイルを 1 機能ずつ、過去の解答、典型（5-5 の問題文を確かめながら）と数学の問題の広い洗い出し、Dynalist の項目）

→ 足したもの・直したものは**すべて今までの 26 ファイル・索引・Dynalist の項目・この設計書に入れた**（幾何の最終照合と同じ。新しいファイルは作っていない。6-2 の 7）。

### 7-1. 旧テンプレートとの照合（1 機能ずつ）

記号：◯＝新しいテンプレートに同じかそれ以上の機能がある、追加＝この照合で足した、体系外＝グラフ・グリッドや自作の DP など別の体系にある。

| 旧ファイル | 旧の機能 | 新しいテンプレート | 判定 |
|---|---|---|---|
| PrimeFactorizeTemplateTest.txt | prime_factorize / enum_divisors / count_divisors / sum_divisors（試し割り）、primes_upto | 02_factorize の factorize（ポラード・ロー）/ factorize_trial / divisors_of / count_divisors / sum_divisors、02_sieve の primes_up_to | ◯ |
| 〃 | build_spf / prime_factorize_spf / enum_divisors_spf / count_divisors_spf / sum_divisors_spf | Sieve S(N)：factorize / divisors / count_divisors / sum_divisors | ◯ |
| SPFTemplateTest.txt | SPF（factorize / divisors / isPrime / primeCount / kthPrime / lpf） | Sieve（factorize / divisors / is_prime / prime_count / kth_prime / spf） | ◯ |
| 〃 | totient / mobius / divisorCount / divisorSum / omegaCount / OmegaCount / isSquareFree / sqrt_decompose_fast / is_square | S.phi / S.mu / S.count_divisors / S.sum_divisors / S.omega / S.big_omega / S.is_square_free / S.square_free_decompose / is_square（00_math_base）。1 つの数の ω・Ω も omega / big_omega（**追加**） | ◯ |
| 〃 | range_prime_count(l, r)（閉区間） | S.count_primes(l, r) / **S.count_primes_closed(L, R)（追加）** | ◯ |
| 〃 | count_coprime_pairs(n)（O(n)） | 02_divisor_transform の count_coprime_pairs（O(N^{2/3})） | ◯ |
| 202D.txt / 425E.txt | パスカルの三角形（mod なし / あり） | pascal_table<T>（long long・DynModInt） | ◯ |
| 368D.txt | power / mod_pow / mod_inv / comb / homo / perm / comb_small / extgcd / mod_inv_（合成数） | pow_int・mod_pow・mod_inv（00_math_base）、Binomial の C / H / P、binom_large_n、ext_gcd / try_mod_inv | ◯ |
| 445E.txt | safe_gcd / safe_lcm（溢れを調べる）/ safe_lcm_mod | gcd_of、**try_lcm（追加）**・lcm_cap / lcm_of_cap、lcm_mod_of / S.lcm_mod | ◯ |
| NumberTheoryTemplateTest.txt | 拡張ユークリッド・mod 逆元（合成数）・一次不定方程式・ax ≡ b・CRT・φ・phi_sieve | ext_gcd / try_mod_inv、linear_diophantine、linear_congruence、crt_merge / crt_all / garner、euler_phi / phi_table | ◯ |
| ARC141A.txt | parse（部分文字列の数値化）/ concat_number | substr_value / DigitPrefix、concat_value / concat_mod | ◯ |
| MathSumTest.txt / 448E.txt | sum_arithmetic / sum_geometric（合成数 mod）/ repunit_mod / digit_count / pow10_mod / concat / repeat_concat | arith_sum / arith_sum_mod / **sum_range_mod（追加）**、geom_sum_mod、repunit_mod / num_digits / pow10_mod / concat_mod / repeat_concat_mod、**repeat_concat_value（追加。mod なしの値）** | ◯ |
| FastModInt / DynModInt / AnyModInt / MyModint（DynBarrettModInt） | 定数 mod・入力の mod・合成数の mod（gcd_mod / is_unit / try_inv / order / solve_linear）、hash、Barrett | ModInt<MOD> / DynModInt（is_unit / try_inv / hash。偶数の mod は Barrett）、位数は mult_order（03_discrete_log）、合成数 mod の一次方程式は linear_congruence | ◯ |
| BigIntTemplateTest.txt | 四則演算・比較・gcd・lcm・pow・to_ll | BigInt（× はカラツバ、÷ は Knuth）、gcd・pow・to_ll、**lcm（追加）** | ◯ |
| FractionTemplateTest.txt | 四則演算・比較・sign / floor / ceil / to_ld / to_decimal / abs / reciprocal | Fraction（同じものと round / is_integer / hash / 入出力 / Fraction128） | ◯ |
| CSR001J.txt | BIT での転倒数 | inversion_count（値は何でも）、**perm_parity（追加）** | ◯ |
| PalindromeTemplateTest.txt / 363D.txt | kth_palindrome / palindrome_rank / prev・next_palindrome / is_palindrome / make_palindrome / count_palindrome_digits、文字列版（is / make / kth / next / prev / count_leq） | kth_palindrome（文字列）/ palindrome_rank / floor・ceil_palindrome / is_palindrome_num / make_palindrome / count_palindromes_of_length、**文字列版 ceil_palindrome_str / floor_palindrome_str / count_palindromes_up_to_str / make_palindrome_str / is_palindrome_str（追加）** | ◯ |
| Base10TemplateTest.txt | to_base / from_base / leq_base / digits / digit_sum | to_base / from_base（−2 進も）/ base_value_leq / digits_of / digit_sum、**str_mod の base（追加）** | ◯ |
| floor_sumTest.txt | ACL floor_sum の使い方（⌊ai/m⌋・⌊(ai+b)/m⌋・ずらし・まとめ方） | linear_floor_sum（同じ引数の順で範囲が広い）、**linear_mod_sum / linear_floor_sum_ext / universal_euclid（追加）**、quotient_blocks | ◯ |
| Mat2TemplateTest.txt | 2×2 行列（apply / pow / det / inv） | MatrixFixed<T, 2> | ◯ |
| LinearRecurrenceTemplate.txt | 線形漸化式の第 k 項（mod） | linear_recurrence_kth（Bostan–Mori）/ recurrence_nth_matrix | ◯ |
| 368F.txt | Grundy（辺を足してメモ化再帰） | grundy_dag（再帰しない）、**grundy_dag_multi / game_outcome（追加）** | ◯ |
| 238D.txt / BitSetGenerator.txt | 立っているビットの位置、CombinationBits / ForcedOne・ForcedZero（popcount K も） | bit_indices / for_each_bit、for_each_subset_k / for_each_mask_between(must, allowed[, k]) | ◯ |
| ARC59C.txt / 126C.txt | 繰り下がり・繰り上がり / 出力の桁数 | ceil_div・floor_div・mod_floor、**carry_to_string（追加）**、出力は FIXED（Main_base） | ◯ |
| 114C.txt / 161D.txt | 状態遷移の DFS / BFS（N 以下の数・K 番目） | numbers_from_digits（114C）、161D は BFS の型（典型の照合で解いた） | ◯（161D は体系外の型） |

旧ファイルの機能で、新しいテンプレートにないものはない。

### 7-2. 過去の解答・Dynalist の項目との照合

| 過去の解答（旧ファイルの問題） | 新しいテンプレートでの書き方 | 判定 |
|---|---|---|
| 202D（aab aba baa） | pascal_table / multiset_perm_unrank（Phase 4 の問題で解いた） | ◯ |
| 238D（AND and SUM） | bit の型（典型の照合で解いた） | ◯ |
| 363D（Palindromic Number） | kth_palindrome(N − 1) | ◯ |
| 368F（Dividing Game） | grundy_dag（約数へ）または S.big_omega の Nim | ◯ |
| 425E（Count Sequences 2） | pascal_table<DynModInt> | ◯ |
| 445E（Many LCMs） | Sieve で分解し、素数ごとに上位 2 つの指数（典型の照合で解いた） | ◯ |
| 448E（Simple Division） | repunit_mod / pow10_mod / mod_inv（典型の照合で解いた） | ◯ |
| ARC141 A（Periodic Number） | substr_value / repeat_concat の値（典型の照合で解いた） | ◯ |
| ARC059 C（Be Together） | ceil_div・floor_div（典型の照合で解いた） | ◯ |
| 114C（755）/ 161D（Lunlun Number） | numbers_from_digits / BFS の型（典型の照合で解いた） | ◯ |

Dynalist の項目（5-0 の表の 22 項目と旧索引の「数学高速化」「重複を除いて数える」「状態遷移を考える」）は、すべて新しい索引（`dynalist_tool_index`）の行から辿れる（mod inverse → 00_modint・00_math_base、CRT → 03_mod_equation、φ 関数 → 02_factorize・02_sieve、floor_sum → 06_quotient、包除原理 → 04_combinatorics、Möbius → 02_divisor_transform、bitDP・桁 DP は体系外の型）。

### 7-3. 典型の問題（問題文を確かめて解き、サンプルを通した）

- Phase 2〜4 で解いた 51 問に加え、5-5 の残りと過去の解答の問題 **48 問**（TYP020・022・030・038・067・069・075・082・085、ABC044 D・052 C・083 B・105 C・114 C・132 D・145 D・147 D・161 D・167 E・169 C・177 E・180 C・186 C・190 D・192 D・197 C・215 D・227 C・238 D・254 D・262 E・263 E・280 D・280 E・284 D・300 D・336 C・357 D・361 F・372 G・383 D・400 E・412 E・445 E・448 E、ARC059 C・141 A・175 A）を、問題文を開いて確かめ、Main_base ＋ テンプレートの本体だけで解いてサンプルをすべて通した（作ったものを見ていない別のエージェント 2 つに分けて解かせた）
- 最終照合で足した・直したあとのファイルで組み直して、**99 問ともサンプルを通ることを確かめ直した**
- 最終照合で足した機能で、さらに **17 問**を解いた：公式の例があるもの 3 問（ABC390 D：for_each_set_partition、ABC209 E：game_outcome、ABC172 D：multiplicative_sum）、自作の例（期待値は Python の愚直解で作った）14 問（多点評価・補間・集合べき級数の exp・連結な部分グラフの数（subset_log）・Σφ（multiplicative_sum、N = 10^10）・Σ⌊⌋²（linear_floor_sum_ext）・linear_mod_sum・オンラインの畳み込みでカタラン数・Σ gcd(x, M)（DivisorsOf、M ≤ 10^12）・q_binomial・Mo で Σ C(n, i)（BinomialPrefix）・大きな数の掛け算（convolve_ll ＋ carry_to_string）・次の回文数（100 桁）・置換の偶奇）

このとき見つかった使いにくさと直したこと：

| 見つかったこと | 直したこと |
|---|---|
| `primes_up_to` / `S.primes()` が vector<int> で、`for (auto p : ps) if (p * p > N)` が溢れる | vector<long long> にした（メモリを減らしたいときは `primes_up_to<int>`）。速さはほぼ同じ |
| `Sieve S(N)`・表の関数が int の N を受け取るので、long long の N を渡すと黙って切り詰められる | long long で受け取り、範囲の外（負・10^8 超え）は止める。`S.spf(x)` も範囲を確かめる |
| `gcd_of` / `lcm_of_cap` / `lcm_mod_of` / `S.lcm_mod` / `xor_pair_sum` / `subarray_xor_sum` が vector<long long> だけで、vector<int> を渡すとコンパイルエラー | テンプレートにして、vector<int>・{…} も受け取るようにした |
| 典型90 069 の `pow(N − 2)` は N = 1 で負の指数になり、エラーにならず逆元の累乗になる | 00_modint のヘッダ・デモに注意を書いた |
| `inclusion_exclusion` の f(S) が「条件 S を満たす」のか「破る」のか分かりにくい | ヘッダに「f(S) = 条件 S を全部満たすものの数 → どれも満たさないものの数」と書いた |
| `linear_floor_sum(n, m, a, b)` の引数の順が分かりにくい | 「ACL の floor_sum と同じ順」とヘッダ・デモ・索引に書いた |

### 7-4. 広い洗い出し（AtCoder の数学 349 問）から決めたこと

AtCoder の ABC200〜476 の問題のうち数学の道具を使うもの 349 問（C〜G 問題。典型の照合で確かめたものを除く）を、テンプレートでそのまま解ける ◯、数行足せば書ける △、足りない × に分けた：**◯ 85、△ 240、× 24**（解説を確かめたものを含む）。

| 足りなかったもの（× の問題） | 決めたこと |
|---|---|
| 集合の組分けを全部たどる（ABC390 D） | 04_combinatorics に for_each_set_partition（と for_each_integer_partition） |
| 区間の中で分母が最小の分数（ABC408 G） | 08_fraction に simplest_between（と rational_floor / rational_ceil・連分数） |
| Σ⌊⌋²・Σ i⌊⌋・モノイドでの床関数の和（ABC402 G・ABC429 G）、Σ (ai+b) mod m | 06_quotient に linear_floor_sum_ext・universal_euclid・linear_mod_sum |
| 乗法的関数の和（ABC370 G。Min_25） | 02_prime_count に multiplicative_sum（ABC370 G は素数を mod 3 で分けるので、そのままでは △） |
| ディリクレ積（ABC428 G） | 02_divisor_transform に dirichlet_convolution（と M の約数の上の変換 DivisorsOf） |
| オンラインの畳み込み（ABC357 G） | 07_convolution に RelaxedConvolution（と convolve_cyclic） |
| 二項係数の部分和を多数のクエリで（ABC463 G。Mo） | 04_binomial に BinomialPrefix |
| 置換の偶奇（ABC470 F）・巡回分解（ABC436 E） | 04_permutation に perm_parity（巡回分解はグラフ・グリッドの 06_perm_cycles を索引から指す） |
| フック長公式（ABC378 G） | 04_combinatorics に young_tableaux_count（問題は RSK も要るので △） |
| 分数の和 Σ a_k/(1 − r_k x)（ABC439 G の一部） | 07_fps に power_sums・fps_sum_of_fractions（と多点評価・補間） |
| power projection・合成逆（ABC345 G・ABC387 G・ABC439 G）、ガウス整数（ABC444 G）、F_{p^2} と chirp z（ABC381 G）、行列ゲームの LP（ABC448 G） | 入れない（出題が少なく、専用の道具が大きい） |

△ の中で、関数はあるのに索引から見つけにくかったもの・型として書いておくとよいもの（pairwise coprime、BEST 定理、min-max 包除、レピュニットと mult_order、鏡像法、(ai+b) mod m → 床関数の和、Π(1 + (k − 1)x) → stirling1_row、op を渡す divisor_zeta、XOR 基底の貪欲、期待値の線形性 など）は索引に行を足した（7-7）。ほかに、10^18 を超える回文数（文字列）・繰り上がり（carry_to_string）・ω / Ω・閉路のあるゲーム（game_outcome）・1 手で分かれるゲーム（grundy_dag_multi）・集合べき級数の exp / log・MinMax / MaxMin の半環・q 二項係数・二項変換も、旧ファイル・洗い出しで欲しくなったので足した。

### 7-5. 足したもの・直したもの（ファイル別）

| ファイル | 足したもの・直したもの |
|---|---|
| 00_math_base | try_lcm、gcd_of / lcm_of_cap を vector<int>・{…} でも、sum_range_mod / sum_range_closed_mod。列の要素を long long にするとき、収まらない値（unsigned の 2^63 以上・Int128）は止める（math_internal::to_ll_checked。列を受け取る関数はすべてこれを使う）、gcd_of の答えが 2^63 なら止める |
| 00_modint | pow の負の指数の注意（ヘッダ・デモ） |
| 01_digits | ceil_palindrome_str / floor_palindrome_str / count_palindromes_up_to_str / count_palindromes_below_str / make_palindrome_str / is_palindrome_str、repeat_concat_value、str_mod(s, m, base)、carry_to_string |
| 02_factorize | omega / big_omega、lcm_mod_of を vector<int> でも（先に全部の値を確かめてから 0 を見る） |
| 02_sieve | Sieve(long long N)（負なら止める）、S.primes() / primes_up_to を vector<long long> に（primes_up_to<int>）、S.count_primes_closed（R > N は止める）、S.spf の範囲の確認、表の関数を long long の N で（範囲の外は止める）、S.lcm_mod を vector<int> でも。**divisor_count_table / divisor_sum_table / omega_table を線形篩に**（O(N log N) で 10^8 に 12〜15 秒かかっていた → 1〜2 秒。10^7 なら 0.1〜0.2 秒） |
| 02_divisor_transform | dirichlet_convolution、DivisorsOf（zeta / mobius / multiple_zeta / multiple_mobius / divisors / index / factors） |
| 02_prime_count | multiplicative_sum（Min_25 篩。T は mint・long long・Int128。整数は溢れたら止める） |
| 04_binomial | BinomialPrefix（m ≥ n なら値は 2^n なので min(m, n) を渡すと速い、と説明）、multinomial の個数の和の溢れを調べる |
| 04_combinatorics | binomial_transform / binomial_inverse / at_least_from_exactly / exactly_from_at_least、for_each_set_partition / for_each_integer_partition（再帰しない）、q_binomial（q-Lucas の C は mod が小さいとき Lucas の定理で）、young_tableaux_count（マスの数の溢れを先に調べる）、inclusion_exclusion の説明 |
| 04_permutation | perm_parity |
| 05_matrix | MinMax / MaxMin、自作の半環の説明 |
| 06_quotient | for_each_quotient を途中で止める、linear_mod_sum、linear_floor_sum_ext（FloorSumExt）、universal_euclid（累乗の最後の使わない 2 乗をしない：溢れを調べるモノイドが止まらないように）、引数の順の説明 |
| 07_convolution | convolve_cyclic、RelaxedConvolution |
| 07_fps | multipoint_eval（f か点の一方が 64 個以下なら 1 点ずつ）、interpolate、power_sums、fps_sum_of_fractions（FpsFraction） |
| 07_subset_transform | subset_exp / subset_log |
| 08_fraction | rational_floor / rational_ceil、simplest_between、continued_fraction（項が 2^63 になる LLONG_MIN / −1 は止める）/ from_continued_fraction |
| 08_bigint | lcm（整数と混ぜても） |
| 09_bit | xor_pair_sum / subarray_xor_sum を vector<int> などでも |
| 10_game | grundy_dag_multi、game_outcome（GameOutcome・GAME_WIN / GAME_LOSE / GAME_DRAW） |
| 00_int_common・03_mod_equation・03_discrete_log・05_linear_algebra・05_xor_basis・06_series・06_linear_recurrence | 変更なし（同梱している基盤・ファイルの更新だけ） |

関数を足したファイルには、すべてランダムテスト（全探索・定義どおりの愚直な計算・128 ビットの厳密な値・別の方法との比較）と、使い方デモ（日本語の説明と類題）を足した。ヘッダ・本体の 1 行目の説明も足した機能に合わせて直した。

主な速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）：multiplicative_sum（Σφ を mint で）10^10 で 0.2 秒・10^11 で 1.1 秒・（Σd）10^12 で 3.5 秒、RelaxedConvolution 2×10^5 回で 0.37 秒、multipoint_eval（次数 10^5・10^5 点）0.5 秒・interpolate（10^5 点）0.76 秒、subset_exp（2^18）0.29 秒・subset_log（2^18）0.51 秒、primes_up_to(10^8) は vector<long long> で 0.29 秒（<int> で 0.21 秒）

### 7-6. 独立したレビューで直したもの

作ったものを見ていない別のエージェントに、足した関数の本体を読ませ、全探索・別の方法との比較・サニタイザで確かめさせた（25 本のテスト）。間違った答えを返すものはなかった。見つかった「黙って溢れる・止まらずに未定義動作になる・遅い」ものは次のとおり直した。

| 見つかったこと | 直したこと |
|---|---|
| `continued_fraction(LLONG_MIN, −1)` が黙って {LLONG_MIN} を返す（項 2^63） | 止める |
| `gcd_of({LLONG_MIN})` が負の値を返す（答えは 2^63） | 止める |
| 列を受け取る関数（gcd_of・lcm_of_cap・lcm_mod_of・xor_pair_sum など）に unsigned long long の 2^63 以上や Int128 を渡すと、long long に切り詰めて誤った答えになる | 要素を `math_internal::to_ll_checked` で変換し、収まらなければ止める |
| `S.count_primes_closed(0, LLONG_MAX)` で R + 1 が溢れる | R > N を先に止める |
| `q_binomial<DynModInt>`（mod = 7 など小さい素数）で q の位数 d が小さいと、C(⌊n/d⌋, ⌊k/d⌋) の n が mod 以上になって止まる | Lucas の定理で求める（mod 2・3・5・7・13 で q-パスカルの漸化式と比べるテストを足した） |
| `young_tableaux_count` の巨大な行で、マスの数の和が溢れてから大きさを確かめていた | 足す前に確かめる |
| `universal_euclid` の累乗が最後に使わない 2 乗をしていて、溢れを調べるモノイドだと正しい計算でも止まる | 最後の 2 乗をしない（テストを足した） |
| `Sieve(−5)` が黙って Sieve(1) になる | 負なら止める |
| 表の divisor_count_table / divisor_sum_table（O(N log N)）が 10^8 で 12〜15 秒 | 線形篩にした（1.2〜1.8 秒。omega_table も線形篩に） |
| `multipoint_eval` が 1 次式を 10^6 点で評価するときも積の木を使う（0.77 秒） | f か点の一方が 64 個以下なら 1 点ずつ |
| `for_each_integer_partition(10^6, f, 1)` が再帰が深すぎて落ちる | 再帰しない形にした |
| `lcm_mod_of({0, −5})` が負の値を確かめずに 0 を返す | 先に全部を確かめる |
| （自分の確認で）`DivisorsOf::index` が std::lower_bound を使っていて、_GLIBCXX_DEBUG では 1 回ごとに並びを全部確かめるため O(D) になる（約数 10 万個のテストがサニタイザのビルドで 17 分） | 自分で二分探索する（39 秒に） |
| （自分の確認で）`for_each_integer_partition` の `vector<Frame> st = {{…}}` に、サニタイザのビルドだけ GCC が -Warray-bounds の誤った警告を出す | push_back で入れる |

### 7-7. 索引・Dynalist の項目の更新

- 索引（`dynalist_tool_index.txt` / `.opml`）：足した関数の行と、7-4 の見つけにくかった使い方の行を足した（129 行 → 164 行）
- Dynalist の項目：`dynalist_items_phase2.txt`（00_math_base・00_modint）、`dynalist_items_phase3.txt`（01_digits・02_factorize・02_sieve・04_binomial・05_matrix・06_quotient・08_fraction・09_bit）、`dynalist_items_phase4.txt`（02_divisor_transform・02_prime_count・04_combinatorics・04_permutation・07_convolution・07_fps・07_subset_transform・08_bigint・10_game）の該当する項目を書き直した。Dynalist には、この 19 項目を入れ替える

### 7-8. 確認したこと

- 26 本すべて：g++ の C++17 / C++20 / C++23 と clang の C++17 / C++20 / C++23（-Wall -Wextra -Wshadow -pedantic で警告 0）、AddressSanitizer ＋ UndefinedBehaviorSanitizer ＋ _GLIBCXX_DEBUG、-O2 の UndefinedBehaviorSanitizer で、ランダムテスト・速度計測・デモがすべて通った（8 通り × 26 本）
- 貼り付けの確認：Main_base ＋ ACL ＋ グラフ・グリッド 46 本 ＋ 幾何 17 本と一緒に貼ってもぶつからない（2 通りの順・g++ / clang・gnu++2b / c++17・LOCAL あり / なし。最終照合で足した関数も呼んで確かめた。BigInt の lcm と std::lcm も区別される）、各ファイルを単独で Main_base に貼って動く（g++ / clang で警告 0。clang の 2 件は Main_base の使っていない INF・LINF）、同梱しているブロック 69 個がすべて持ち主のファイルと同じ、00_int_common は幾何と同じ
- 問題：Phase 2〜4 の 51 問、典型と過去の解答の 48 問、足した機能の 17 問の**合わせて 116 問がすべてサンプルを通る**（Main_base ＋ テンプレートの本体だけで）
- 索引・Dynalist の項目に書いた関数名が、すべてテンプレート（または graph_grid_v2・geometry_v2）の中にある
- `Desktop\claude\math_v2` に保存したファイルが作業環境のものと同じことを確かめ、ユーザーの PC の環境（g++ 11.4、C++17、-O2 -Wall -Wextra）でも 26 本がすべて警告 0 で通った
