# 文字列テンプレート 第二世代 設計書

- 作成日：2026-10-11
- 状態：**Phase 1（体系設計）完了**（2026-10-11 ユーザー確認済み。6-1）。**Phase 2（基盤 00_str_base）完了**（2026-10-12。2-P2 章）。**Phase 3（★ 11 本）完了**（2026-10-12。2-P3 章）。**Phase 4（☆ 7 本）完了**（2026-10-12。2-P4 章）。**最終照合 完了**（2026-10-13。7 章）
- 方法：グラフ・グリッド第二世代（`Desktop\claude\graph_grid_v2\design.md`）・幾何第二世代（`Desktop\claude\geometry_v2\design.md`）・数学第二世代（`Desktop\claude\math_v2\design.md`）と同じ方法・観点で作り直す
- 次回はこのファイルの「6. 未決事項と次にやること」から再開する

---

## 0. 目的・観点・これまでから引き継ぐこと

### 0-1. 目的

文字列（一致・出現、部分文字列の比較、接尾辞、回文、辞書順・回転、部分列、括弧列、接頭辞の集合、更新のある文字列）のテンプレートを「第二世代」としてゼロから設計し直す。既存テンプレートの整理ではなく、新しい体系を作る。

成果物はテンプレートそのものと、**問題の特徴から使う道具に辿れる選択体系**（Dynalist の索引）。例えば次のように辿れる状態を目指す。

- T が S のどこに何回現れるか → 02_pattern_match（`find_all(S, T)` / `count_occurrences(S, T)`）
- 2 つの部分文字列が等しいかを何度も聞かれる → 01_rolling_hash（`RollingHash H(S); H.get(l1, r1) == H.get(l2, r2)`）
- 区間 [l, r) が回文かを何度も聞かれる → 04_palindrome（`Palindromes P(S); P.is_palindrome(l, r)`）
- 異なる部分文字列の数・辞書順 k 番目の部分文字列 → 03_suffix_array（`SuffixArray SA(S); SA.count_distinct_substrings()`）
- 多くの文字列の共通接頭辞・LCP の和 → 08_trie（`Trie T; T.insert(s); T.sum_lcp_pairs()`）
- 長さ K の辞書順最小の部分列 → 06_subsequence（`lex_smallest_subsequence(S, K)`）
- 最長共通部分列（飛び飛び）と最長共通部分文字列（連続）→ 06_lcs_edit（`lcs_restore`）と 03_suffix_array（`longest_common_substring`）を名前で区別
- 正しい括弧列か・区間が正しい括弧列か → 07_bracket（`is_balanced(S)` / `BracketSegTree`）
- 前後への追加と反転を何度も → 09_reversible_deque（`ReversibleDeque<char>`）

### 0-2. 観点（グラフ・グリッド・幾何・数学と同じ 3 つ）

| 観点 | 文字列での具体化 |
|---|---|
| アクセスしやすさ | 「問題の特徴（一致・出現／部分文字列の比較／回文／辞書順・回転／部分列／括弧／接頭辞の集合／更新）→ ファイルと関数」で辿れる索引。ファイル名は分類番号＋名前。依存するファイルは同梱し、1 回のコピーで動く |
| 使いやすさ | `string` と `vector<整数>` のどちらもそのまま渡せる。範囲は [l, r)、「ない」は -1 にそろえる。文字の種類は `Alphabet` で指定し、範囲外の文字で黙って配列の外を読まない（旧 StringHelper・旧 Trie の不具合）。「部分文字列（連続）」と「部分列（飛び飛び）」を名前で区別する（旧 EDPCF の取り違え）。複数の値は名前付きのメンバを持つ構造体で返す |
| 連携 | mod の数え上げは型を引数に取る（math_v2 の `ModInt998` をそのまま使える：`count_distinct_subsequences<mint>(S)`）。KMP・Aho–Corasick のオートマトンを math_v2 の行列の累乗に、文字列の包含関係を graph_grid_v2 のマッチング・フローにつなぐ形を索引に載せる。グラフ・グリッド（46 本）・幾何（17 本）・数学（26 本）と同時に貼れる |

### 0-3. 最重要の方針（グラフ・グリッドの 1-0 を最初から適用する）

第二世代のテンプレートは互いに関係し、旧テンプレートとの互換性もないため、**後から追加・修正しなくて済むこと**が最も重要。

- 追加でほしくなる可能性のある機能はすべて入れる。入れるか迷った機能は入れる
- 速さの改善が要らないくらい速くする
- 正しいこと（愚直解とのランダムテスト、サニタイザ、Main_base＋ACL に貼った状態でのコンパイル確認）
- 使いやすいこと

### 0-4. これまでの経験から、最初から取り入れること

| これまでに起きたこと | 文字列での対応 |
|---|---|
| 最終照合で足りないものが見つかった（グラフは 7 本、幾何は 109 問の洗い出しで「足りない」が 32 問、数学は 349 問の洗い出しで機能を追加） | **典型との照合を Phase 1 で行う**（5-5。問題名は AtCoder の問題一覧で確かめた）。`Desktop\claude\AtCoder_文字列問題一覧.txt`（ABC200〜477 の 128 問）も使う。最終照合でも文字列の問題を広く洗い出す |
| 使い方関数の例のとおりに書くとコンパイルできなかった（幾何） | 使い方関数の例はすべて実際にコンパイル・実行されるコードにする。**各 Phase の終わりに代表的な問題をテンプレートだけで解いて確かめる** |
| std の lower_bound などがデバッグビルドで O(n) の検査をして終わらなくなった | 本体の二分探索・ヒープは自分で書く |
| MSVC ではビルドしない・`__int128` は使える | ハッシュの掛け算などに `unsigned __int128` をそのまま使う |
| 名前の衝突（ACL は `using namespace atcoder;` で見える） | **ACL の文字列の関数 `suffix_array` / `lcp_array` / `z_algorithm` と同じ名前にしない**（同じ引数だと呼び出しがあいまいになる）。グラフ・グリッド・幾何・数学のグローバルな名前（約 820 個。`CharSet` `trim` `contains` `transpose` `inversion_count` `is_palindrome_str` など）とも重ねない |
| 作ったものを見ていない別のエージェントのレビューで不具合が見つかった | Phase 3・4 の終わりに毎回行う |
| 旧テンプレートの不具合（ハッシュの基数を作り直して前の値と食い違う・範囲外の文字・閉区間と半開区間の混在・区切り文字が負の char になる） | 1-2 の規約で最初から防ぐ（ハッシュの基数は 1 回だけ決める、Alphabet で範囲を確かめる、[l, r) にそろえる、区切りは整数で持つ） |
| 各 Phase で Dynalist の項目と索引を作った | 同じ（dynalist_items_phaseN.txt、dynalist_tool_index.txt / .opml） |

---

## 1. 決定事項

### 1-1. 進め方・既存テンプレートの扱い

- 既存テンプレートは「過去の実験結果」として扱う。既存の API を正解として引き継がない
- 参照したのは、新しい一覧と機能を決めた後（Phase 1 の終盤、5 章）だけ
- 原本は変更・削除しない。新しく作るものは `Desktop\claude\string_v2\` に出力する（最後に `C:\atcoder_template\典型90\ソース\string_v2\` へ置く想定。ほかの第二世代と同じ）
- 既存の場所：`C:\atcoder_template\典型90\ソース`（5-0）
- Phase 1（体系設計）→ Phase 2（基盤）→ Phase 3（★）→ Phase 4（☆）→ 最終照合 の順に進め、区切りごとにユーザーの確認を取る
- **ACL は使わない**。接尾辞配列・Z も std だけで自作する。速さは ACL と同じかそれ以上を目標にして、計測で比べる

### 1-2. 共通規約

ほかの第二世代と同じもの：

| 項目 | 決定 |
|---|---|
| 優先順位 | 速度 → int / long long と変わらない使い勝手 → 正しさ |
| ファイル | 1テンプレート1ファイル。単体で実行できる .cpp |
| テスト | ランダムテスト関数と、全機能を使う使い方関数を同じファイルに入れ、main から呼ぶ（形は同じ：`xxx_test::random_test()`（Release は 8 倍の量）・`benchmark()`・`xxx_demo::demo()`・`setup_console()`） |
| 使い方の出力 | 関数名の表示に加え、何を返す・提供するかを日本語でも書く |
| 使うライブラリ | std のみ（ACL 不使用）。bits/stdc++.h に頼らず、使うヘッダを書く |
| 名前 | 型は PascalCase（`SuffixArray`, `Trie`）、関数・メンバは snake_case（`find_all`, `count_prefix`） |
| Main_base への依存 | テンプレートの中では Main_base の別名・マクロを使わない |
| 再帰 | 非再帰を基本にする（Trie・回文木・接尾辞オートマトンをたどる処理も）。深さが log N で済むもの（SA-IS の縮約など）だけ再帰でよい |
| 範囲 | 半開区間 [l, r)、0 始まり。k 番目は 0 から |
| 「なし」 | 位置・番号は -1。解がないことがある計算は `ok` を持つ構造体 |

文字列で新しく決めるもの（★ は特に確認してほしいもの。6-1）：

| 項目 | 決定（案） |
|---|---|
| 文字列と数の列 | 列を受け取る関数・構造体は **`string` と `vector<整数>`（int・long long）の両方**を受け取る（テンプレート）。数の列は値が大きくても負でもよい（接尾辞配列などは中で座標圧縮する）。char は中で unsigned char として扱う（128 以上の文字で負の添字にしない）。要素は整数型だけ（Phase 4 で `code_of` に static_assert を入れた。`"c"` のような文字列リテラル・ポインタ・小数はコンパイルエラー。unsigned long long は 2^63 未満の値で）。1 文字（char）も長さ 1 の列として渡せる（最終照合。`find_all(S, 'a')`） |
| ★ 文字の種類 | 子を配列で持つもの（Trie・Aho–Corasick・NextTable・CharPrefix・接尾辞オートマトン・回文木）は **`Alphabet` を引数に取る**（既定は英小文字）。`Alphabet::lower()` / `upper()` / `digits()` / `letters()`（A〜Z と a〜z の 52 種）/ `alnum()`（62 種）/ `ascii()` / `of("()[]")`（並べた文字だけ）/ `ints(K)`（整数 0〜K−1）。**範囲外の文字はメッセージを出して止める**（assert と違い NDEBUG でも止まる）。実行時に決まるので、K が入力で決まる数の列にも使える。表（(N+1)σ・ノード数 × σ）が 2×10^8 個（約 800MB）を超えるときは作らずに止める（Phase 4） |
| L 文字目〜R 文字目 | 問題の「L 文字目から R 文字目」（1 始まり・両端を含む）は `(L - 1, R)` で渡す。文字列の範囲は `_closed` 版を作らない（数学と違い、変換が 1 か所で間違えにくいため。使い方関数で示す） |
| 部分文字列の返し方 | 位置で返す：`StrSpan{l, r}`（`len()`）。文字列が欲しいときは `substr_of(S, sp)`（コピーを作らないため。最長の回文・k 番目の部分文字列など） |
| 個数 | 部分文字列・出現の数は long long（N²/2 で int を超える）。文字の数（count_chars・CharPrefix::count など）も long long で返す（最終照合で int から変えた。Main_base の chmax・min に ll の変数と並べて渡せる）。mod で数えるもの（異なる部分列の数など）は**型を引数に取る**：`count_distinct_subsequences<T>(S)`。T は math_v2 の `ModInt998` / `DynModInt`、long long、double など。整数型を渡して溢れたら止める（数学と同じ） |
| ★ ハッシュ | **mod 2^61 − 1 の 1 本**。基数は実行ごとの乱数で、**すべての文字列で共通**（別の文字列の部分文字列どうし・連結したものも比べられる）。累乗の表は基数を変えずに伸ばす（旧 RollingHashSegTree の不具合を防ぐ）。1 回の比較で衝突する確率は約 N / 2^61 |
| 比較の戻り値 | 辞書順の比較は -1 / 0 / 1（int）。「a が b の接頭辞」は a < b |
| 部分文字列と部分列 | 名前で区別する：連続するものは `substring`（`longest_common_substring`）、飛び飛びのものは `subsequence`（`lcs_restore` は「最長共通部分列」と説明に書く） |
| 名前の衝突 | Main_base（マクロ `rep` `erep` `rrep` `nall` `rall` `YES` `NO` `YN` `BUFF` `FIXED` `next_p` `debug`、名前 `ll` `ull` `sgn` `popcnt` `INF` `LINF` `PQ` `PQ_g` `chmin` `chmax` `sizei` `GraphIO`）・**ACL**（文字列の `suffix_array` `lcp_array` `z_algorithm`、ほかに `segtree` `fenwick_tree` `internal` など）・std（`hash` `rank` `next` `prev` `count` `find` `search` `mismatch` `reverse` `rotate` `replace` `size` `empty` `data` `tolower` `toupper` `isalpha` など。`using namespace std;` で見えるので、**グローバルな関数を `hash` `rank` `next` `prev` という名前にしない**（構造体のメンバならよい））・**グラフ・グリッド 46 本・幾何 17 本・数学 26 本のグローバルな名前（約 820 個）**・よく使う変数名（S, T, N, M, s, t, n, m, sa, lcp, z, pi, rk, cnt, pos, nxt, prv, dp, ans, res, rev, rle, runs, periods, lis, lcs, h, H, trie）・Windows のマクロ（`min` `max` `IN` `OUT`）と同じ名前をグローバルに置かない。配列によく使う名前は関数名にしない（`all_periods`・`maximal_repetitions` のようにする） |
| 出力 | テンプレートは cout しない（止めるときの cerr のメッセージを除く） |

### 1-3. 共通API

```cpp
// ---- 文字列も数の列もそのまま渡せる ----
string S, T;  vector<int> A, B;
auto z  = z_array(S);                       // z_array(A) も同じ
auto ps = find_all(S, T);                   // T の出現位置（重なりも）。find_all(A, B) も
long long c = count_occurrences(S, T);

// ---- 前処理してから何度も聞くものは構造体 ----
RollingHash H(S);   H.get(l, r) == H.get(l2, r2);  H.lcp(i, j);  H.compare(l1, r1, l2, r2);
SuffixArray SA(S);  SA.count_distinct_substrings();  SA.count(T);  SA.kth_substring(k);
Palindromes P(S);   P.is_palindrome(l, r);  auto sp = P.longest();  substr_of(S, sp);
Trie TR;            TR.insert(s);  TR.count_prefix(p);  TR.sum_lcp_pairs();
NextTable nx(S);    nx.next(i, 'c');  nx.is_subsequence(T);

// ---- 文字の種類（既定は英小文字。範囲外の文字が来たらメッセージを出して止まる）----
Trie TU(Alphabet::upper());                 // 英大文字
NextTable na(A, Alphabet::ints(K));         // 0〜K−1 の整数の列

// ---- mod で数える（math_v2 の ModInt を型に。long long でもよい）----
using mint = ModInt998;
mint d = count_distinct_subsequences<mint>(S);
```

- 1 回で終わる計算は関数（`find_all(S, T)`、`min_rotation(S)`）、前処理してから何度も聞くものは構造体（`RollingHash H(S)`、`SuffixArray SA(S)`）
- 複数の値を返す関数は、名前付きのメンバを持つ構造体で返す（構造化束縛でも受け取れる。pair・tuple で返さない）
- 構造体は作るときに必要なものだけを作る（旧 StringAlgo はコンストラクタで全部を作って重かった）。O(1) の LCP の疎テーブルなどは最初に使うときに作る

### 1-4. ファイルの形と共通化

- ほかの第二世代と同じ：**.cpp 1 ファイル**、コピー範囲の目印（`ここからテンプレート本体（コピー範囲）` / `ここまでテンプレート本体`）、依存するテンプレートは**貼る範囲に同梱**して `#ifndef` ガードで重複を防ぐ
- ガード名は同じ付け方（`GG_` ＋ファイル名の大文字。例：`GG_STR_BASE`、`GG_ROLLING_HASH`）。`update_bundles.ps1 / .bat` を**そのまま**使える（string_v2 にコピーするだけ）
- ファイル名は「分類番号 2 桁＋名前」。名前の部分はグラフ・グリッド・幾何・数学のファイルと重ならない（確認済み）
- 内部の名前空間は `str_internal`（`gg_internal`・`geo_internal`・`math_internal` とは分ける）
- **すべてのファイルが 00_str_base を同梱する**
- 数学の 00_int_common は同梱しない（文字列では Int128 の入出力を使わない。辞書の番号など 128 ビットが要るところは `unsigned __int128` / `__int128` を中で直接使う）
- 文字コードは UTF-8（BOM 付き）、main で Windows のコンソールの設定をする

### 1-5. 正しさ・速さの確認のしかた（文字列で追加するもの）

- **境界の値を必ず入れる**：空文字列（受け付ける関数）・長さ 1・全部同じ文字（"aaaa…"：Z・KMP・接尾辞配列・回文の最悪の形）・2 文字の繰り返し（"abab…"）・フィボナッチ文字列・Thue–Morse 列（周期・runs）・de Bruijn 列・文字の種類 σ = 1, 2, 3, 26・数の列に負の値と 10^18 級の値・Alphabet の範囲外の文字（止まることを確かめる）・128 以上の char
- **愚直解**：全部の部分文字列を `set<string>` に入れる（異なる部分文字列・k 番目）、`std::string::find` を全部の位置で（出現）、O(N²)〜O(N³) の回文判定・DP、全部の部分列の列挙（異なる部分列・最短の部分列でない文字列）、O(N²) の LIS、O(N M) の編集距離
- **ハッシュ**：全部の部分文字列の組の一致を `std::string` の比較と照らし合わせる（衝突しないことは確率なので、件数を多くする）
- **Alphabet を受け取るもの**は、英小文字・英大文字・並べた文字（括弧など）・整数（`ints(K)`）の全部でテストする
- **速さは ACL と旧テンプレートと比べる**：`suffix_array` / `lcp_array`（N = 10^6）、`z_algorithm`（N = 10^7）、旧 RollingHashTest2・ManacherTemplate・ALDS14C（Trie）・StringHelper・A24。ほかの目安：KMP・Manacher（N = 10^7）、ハッシュの取得 10^7 回、Aho–Corasick（パターンの長さの合計 10^6・本文 10^6）、Trie（10^6 文字）、接尾辞オートマトン・回文木（N = 10^6）、LCS のビット並列（5×10^4 × 5×10^4）、編集距離（10^4 × 10^4）、LIS（N = 10^6）、ワイルドカード一致（N = 5×10^5）、書き換えのある文字列（N = Q = 5×10^5）
- **名前の衝突**：Main_base＋ACL に、グラフ・グリッド 46 本・幾何 17 本・数学 26 本・文字列の全ファイルの本体を 2 通りの順で全部貼り、gnu++2b / c++17・LOCAL あり / なしでコンパイル・実行する。単独のテンプレート（binary_search.cpp・static_seq_toolkit.cpp（文字列の機能を持つ）・segtree・lazy_segtree・dynamic_unique_seq.cpp など）とも一緒に貼って確かめる
- **独立したレビュー**：Phase 3・4 の終わりに、作ったものを見ていない別のエージェントに本体を読ませ、愚直解で確かめさせる
- **問題で確かめる**：Phase 3・4 の終わりに、5-5 の代表例の問題文を確かめてテンプレートだけで解き、サンプルを通す

---

## 2. 分類とテンプレート一覧（確定。2026-10-11）

記号：★＝優先して整備（Phase 2・3）、☆＝余裕が出てから（Phase 4。ただしほかの第二世代と同じく全部作る）、△＝専用テンプレートは作らず「型」として持つ

関数名は案。Phase 2〜4 の実装で確定し、この表に反映する（ほかの第二世代と同じ）。代表的な問題の題名は AtCoder の問題一覧で確かめた（中身は Phase 3・4 で問題文を確かめ直す）。計算量の N・M は文字列（列）の長さ、σ は文字の種類の数。

### 0. 基盤（Phase 2）

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `00_str_base` | **作成済み（Phase 2）**。文字列の基本。**文字の種類**：`Alphabet`（`lower()`（既定）/ `upper()` / `digits()` / `letters()` / `alnum()` / `ascii()` / **`bytes()`**（0〜255）/ `of(chars)` / `ints(K)` / `ints(lo, hi)`）：`id(c)`（範囲外は止める）/ `contains(c)` / `value(i)` / `ch(i)` / `chars()` / `size()`、`to_ids(s, A)` / `from_ids(ids, A)`、**`compress_values(a)`**（`CompressedValues{ids, values}`。`.alphabet()` が `ints(K)`）。**文字の数・位置**：`count_chars(s, A)`、`CharPrefix cp(s, A)`（`count(l, r, c)` / `has` / `counts(l, r)` / `distinct(l, r)` / **`next(i, c)` / `prev(i, c)` / `kth(c, k)` / `total(c)`（どれも O(1)）**。メモリ (N+1)σ 個の int）、`CharPositions cps(s, A)`（`next` / `prev` / `count` / `has`（O(log N)）/ `kth` / `total` / `positions(c)` / **`occurrence_index(i)`**。メモリ O(N + σ)）。**ランレングス**：`run_length(s)`（`RunBlock<T>{value, len, start}`。len・start は long long）/ `run_length_decode(runs)`。**変換**：`to_lower` / `to_upper` / `swap_case` / `caesar_shift(s, k, A)` / `CharMap`（`replace(c, d)` / `set(c, d)` / `cm(c)` / `apply(s)` / `compose(g)` / `table()`）/ `reversed(s)` / `min_with_reverse(s)` / `sorted_chars(s)` / `repeat_string(s, k)`。**比較**：`common_prefix_length` / `common_suffix_length` / `hamming_distance` / `is_subsequence(t, s)` / `is_palindrome(s)`。**区切り**：`split(s, sep, skip_empty)` / `join(vs, sep)`（数の列も）。**共通の型**：`StrSpan{l, r}` と `substr_of(s, sp)`。`StrSpan`・`RunBlock` は `<<` で出せる（デバッグ用） **最終照合で足した**：`Alphabet::with`、`count_id` などの番号の版、1 文字の変換、`common_prefix_length(a, b, wild)`、`subsequence_positions`、`split(s, 文字列)`、`strip` / `strip_left` / `strip_right`、`split_runs_by`、`rearrange_no_adjacent_equal` / `lex_smallest_no_adjacent_equal`。数は long long で返す。1 文字も列として渡せる。 | なし | ABC219 C、ABC232 B、ABC259 C、ABC310 C、ABC329 C、ABC342 C、典型90 084、ABC380 D（△） |

### 1. ハッシュ

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `01_rolling_hash` | **作成済み（Phase 3）**。mod 2^61 − 1。基数は実行ごとの乱数・全部の文字列で共通。`RollingHash H(s)`（string・vector<整数>）：`get(l, r)`（`StrHash{value, len}`）/ `get_len(l, len)` / `prefix(r)` / `suffix(l)` / `all()` / `same(l1, r1, l2, r2)`（O(1)）/ `lcp(i, j)`・`lcp(i, H2, j)` / **`common_suffix(i, j)`・`common_suffix(i, H2, j)`** / `compare(l1, r1, l2, r2)`・`compare(l1, r1, H2, l2, r2)`（-1 / 0 / 1。O(log N)）/ `find_all(pat)` / `count(pat)`（pat は列でも StrHash でも）/ `code(i)` / `size()`。`StrHash`：`+` `+=`（連結）/ `== != <` / `std::hash` の特殊化 / `repeat(k)`（O(log k)）/ `<<`。`hash_of(s)` / `hash_of_char(c)` / `empty_hash()`。`RollingHash2D G(grid)`：`get(i1, j1, i2, j2)`（`GridHash{value, h, w}`）/ `find_all(pat)` / `all()` / `rows()` / `cols()` **最終照合で足した**：`center_palindrome_lengths_with_changes(s, k)`。 | str_base | ABC141 E、ABC284 F、典型90 047 |
| ☆ `01_set_hash` | **作成済み（Phase 4）**。集合・多重集合のハッシュ（値ごとに実行ごとの乱数を割り当てて mod 2^61 − 1 で足す。表を持たない）：`value_hash(x)`、`prefix_set_hashes(a)`（先頭 i 個の集合。重複は 1 回。O(N log N)）/ `prefix_multiset_hashes(a)`、`MultisetPrefix mp(a)`（`get(l, r)`（O(1)）/ `same(l1, r1, l2, r2)` / `size()`）、`MultisetHash`（`value` / `count`、`add(x, k)` / `remove(x, k)` / `merge(o)` / `==`）、`SetHash`（`value`、`add(x)` / `remove(x)`（なければ false）/ `count(x)` / `distinct()` / `==`） | str_base | ABC250 E、ABC367 F |

### 2. 一致・出現（パターン検索）

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `02_pattern_match` | **作成済み（Phase 3）**。Z と KMP（string・vector<整数>。text と pat の型を混ぜてもよい）：`z_array(s)`（z[0] は s の長さ）/ `prefix_function(s)` / `match_lengths(text, pat)` / `find_all(text, pat)`（重なりも含む開始位置。pat が空なら 0〜|text|）/ `find_first(text, pat, from)` / `find_last` / `count_occurrences` / `count_non_overlapping` / `replace_all(text, from, to)`。**周期・境界**：`all_borders(s)` / `min_period(s)` / `all_periods(s)` / `repetition_root(s)` / `prefix_min_periods(s)`。**回転**：`rotation_offsets(a, b)`（a, b がともに空なら {0}）/ `is_rotation(a, b)` / `find_all_cyclic(text, pat)`。**重なり**：`max_overlap(a, b)`。**オートマトン**：`KmpAutomaton KA(pat, A)`（`next(state, c)` / `next_id(state, id)` / `accepting(state)` / `size()` / `alphabet()`） **最終照合で足した**：`rewrite_until_stable` / `erase_until_stable`、`Occurrences`。 | str_base | ABC135 F、ABC150 F、ABC257 G、ABC284 F、ABC343 G、ABC430 E |
| ☆ `02_aho_corasick` | **作成済み（Phase 4）**。`AhoCorasick ac(A)`：`add(p)`（番号を返す。同じ文字列も可。空は止める）/ `build()` / `count_each(text)`（パターンごとの出現回数）/ `count_total(text)` / **`first_occurrence_each(text)`**（設計の first_end_each。最初の**開始**位置、なければ -1）/ `for_each_match(text, f)`（**f(開始位置, 番号)**。終わる位置の順）/ `contains_any(text)` / `next(state, c)` / `next_id(state, id)`（O(1)）/ `match_count(state)`（その状態で終わるパターンの数。接尾辞リンクの先も含む）/ `match_ids(state)`（長い順）/ `longest_match(state)` / `depth(state)` / `link(state)` / `bfs_order()` / `size()` / `pattern_count()` / `node_of(id)` / `pattern_length(id)`、`contained_patterns(id)`（パターン id の中に現れるほかのパターン。同じ文字列も入る。包含関係） | str_base | ABC268 Ex、ABC354 G、ABC362 G |
| ☆ `02_wildcard_match` | **作成済み（Phase 4）**。畳み込み（NTT を中に持つ。math_v2 には依存しない。text の長さは 2^23 まで）：`wildcard_match(text, pat)`（'?'。文字列だけ）/ `wildcard_match(text, pat, wild)`（両方の wild を任意の 1 文字として一致する開始位置。乱数・O(N log N)）/ `mismatch_counts(text, pat)` / `mismatch_counts(text, pat, wild)`（各ずらしでの不一致の数。正確・O(σ' N log N)、σ' は両方に現れる値の数）/ `match_with_k_mismatches(text, pat, k)` / `match_with_k_mismatches(text, pat, k, wild)` | str_base | ABC196 F |

### 3. 接尾辞

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `03_suffix_array` | **作成済み（Phase 3）**。`make_suffix_array(s)`（SA-IS。string・vector<整数>（座標圧縮））/ `make_lcp_array(s, sa)`（Kasai）。`SuffixArray SA(s)`：`sa` / `rank` / `lcp_adj`（隣どうし）/ `lcp(i, j)`（疎テーブルで O(1)。初回に作る）/ `compare(l1, r1, l2, r2)`（O(1)）/ `find_range(pat)` / `count(pat)` / `find_all(pat)`（昇順。pat が空なら置ける位置すべて）/ `count_occurrences(l, r)` / `occurrence_range(l, r)` / `count_distinct_substrings()` / `distinct_substrings_by_length()` / `kth_substring(k)`（`StrSpan`）/ `kth_substring_with_dup(k)` / `count_substrings_with_dup()` / `substring_rank(l, r)` / `longest_repeated()` / `longest_repeated_non_overlapping()` / `sum_lcp_with_all()` / `sum_lcp_all_pairs()` / `size()` / `length_at(p)` / `value(p)`。**複数の文字列**：`SuffixArray::of_strings(vs)`（区切りは整数で持つので本数に制限なし。`owner(p)`（区切りは -1）/ `offset(p)` / `position(t, off)`）、`longest_common_substring(a, b)`（`CommonSubstring{len, pos_a, pos_b}`）/ `longest_common_substrings_all(a, b)`（a での `StrSpan` の一覧）/ `longest_common_substring_k(vs)`（vs[0] での `StrSpan`） | str_base | ABC141 E、ABC213 F、ABC240 Ex、ABC272 F、ABC280 Ex、ABC362 G |
| ☆ `03_suffix_automaton` | **作成済み（Phase 4）**。`SuffixAutomaton sam(s, A)`（`extend(c)` で 1 文字ずつも。複数の文字列は `add_string`。extend は最後に足した文字列の続き）：`contains(t)` / `count_occurrences(t)`（空なら Σ(|s_i| + 1)）/ `first_occurrence(t)`（通し番号）/ `count_distinct_substrings()` / `total_length_distinct()`（足すたびに更新。O(1)）/ `count_substrings_with_dup()` / **`kth_substring(k, with_dup)`（string を返す）/ `kth_substring_ids(k, with_dup)`** / `longest_common_substring(t)`（t での `StrSpan`）/ 状態の情報（`size` / `len` / `link` / `next` / `next_id` / `endpos_count` / `first_end` / `order_by_len` / `last`）を DP 用に公開 | str_base | ABC362 G（03_suffix_array・02_aho_corasick の別解）、Library Checker Number of Substrings。ABC280 Ex（Q が多く答えが長い）は 03_suffix_array |

### 4. 回文

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `04_palindrome` | **作成済み（Phase 3）**。Manacher（string・vector<整数>）。`Palindromes P(s)`：`is_palindrome(l, r)`（O(1)）/ `radius_odd(i)` / `radius_even(i)` / **`length_at_center(c)` / `span_at_center(c)`**（中心 c = 0〜2N−2）/ `longest()`（`StrSpan`。同じ長さなら左）/ `count_all()` / `longest_prefix_palindrome()` / `longest_suffix_palindrome()` / `longest_starting_at()` / `longest_ending_at()` / `size()`。関数：**`center_palindrome_lengths(s)`**（設計の manacher_radii。2N − 1 個の中心ごとの最長の長さ）/ **`min_palindrome_with_prefix(s)` / `min_palindrome_with_suffix(s)`**（設計の min_append / min_prepend。足した後の回文そのものを返す）/ `min_changes_to_palindrome(s)` / `can_rearrange_to_palindrome(s)` / `longest_palindromic_subsequence(s)`（N ≤ 10^4）/ `longest_palindromic_subsequence_length(s)` / `palindrome_partition_min(s)`（`vector<StrSpan>`。O(N²)）/ `count_palindromic_subsequences<T>(s)`（O(N²)） **最終照合で足した**：`count_palindromes_leq<T>(s, A)`。 | str_base | ABC237 C、ABC320 B、ABC349 G、ABC398 F |
| ☆ `04_palindromic_tree` | **作成済み（Phase 4）**。回文木（eertree）。`PalindromicTree pt(s, A)` / `add(c)`：`count_distinct()` / `occurrences()`（それぞれの回文の出現回数。参照で返す）/ `node_of_prefix(i)`（s[0, i) の最長の回文接尾辞）/ `size()` / `length()` / `len(v)` / `link(v)` / `parent(v)`（両端を 1 文字ずつ除いた回文）/ `next(v, c)` / `first_end(v)` / `span(v)` / `palindromes()`（異なる回文の最初の出現）/ `count_total()` / `min_partition()` / `min_partition_spans()` / **`min_partition_prefix()`**（回文への最小の分割数。O(N log N)）/ `count_partitions<T>()` / **`count_partitions_prefix<T>()`**（回文への分割の数）。回文の包含関係は parent と link を辿れば全部わかる | str_base | ABC237 Ex、Library Checker Eertree |

### 5. 辞書順・回転

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `05_lex_order` | **作成済み（Phase 3）**。`concat_order(vs)`（a + b < b + a。同じなら元の順。空の列は最後。O(N log N × 最長の長さ)）/ `min_concat(vs)`、`lex_less_by(order)` / `lex_less_by(A)`（`LexLess`。`compare(a, b)` も）、`min_rotation(s)` / `max_rotation(s)`、**`bijective_to_string(n, A)` / `bijective_from_string(s, A)`**（文字の種類は Alphabet で。英大文字は `Alphabet::upper()`）、`dict_rank(s, max_len, A)` / `dict_unrank(k, max_len, A)`、`next_string(s, A)` **最終照合で足した**：`lex_rank_same_length<T>(s, A)`。 | str_base | ABC171 C、ABC219 C、ABC223 B、ABC225 F、ABC285 C |
| ☆ `05_lyndon` | **作成済み（Phase 4）**。`lyndon_factorization(s)`（Duval。各因子の `StrSpan`）/ `is_lyndon(s)` / `lyndon_array(s)`（各 i から始まる最長の Lyndon 語の長さ）/ `min_suffix_of_prefixes(s)`（各 i で s[0, i) の辞書順最小の接尾辞）/ `maximal_repetitions(s)`（runs：極大な繰り返し `Repetition{l, r, period}` の一覧。O(N log N)）/ `count_square_substrings(s)`（xx の形の部分文字列の数）。`Repetition` は `len()` / `==` / `<` / `<<` | str_base, **suffix_array** | Library Checker Lyndon Factorization・Enumerate Runs（ABC では G・Ex でまれ） |

### 6. 部分列

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `06_subsequence` | **作成済み（Phase 3）**。`NextTable nx(s, A)`（中は 00_str_base の CharPrefix。メモリ (N+1)σ 個の int）：`next(i, c)` / `prev(i, c)`（O(1)）/ `is_subsequence(t, from)` / `match_end(t, from)` / `match_length(t, from)` / `char_prefix()`。関数：`prefix_match_lengths(s, t)` / `suffix_match_lengths(s, t)`、`lex_smallest_subsequence(s, k)` / `lex_largest_subsequence(s, k)`、`shortest_non_subsequence(s, A)` / `shortest_non_subsequence_ids(s, A)`、`count_subsequence_occurrences<T>(s, t)`、`count_distinct_subsequences<T>(s, include_empty, min_step)`（O(N log N)）、`distinct_subsequences_by_length<T>(s)`（O(N²)） | str_base | 典型90 006・008、ARC081 E、ABC211 C、ABC214 F、ABC230 F、ABC324 E |
| ★ `06_lcs_edit` | **作成済み（Phase 3）**。**LCS**：`lcs_length(a, b)`（ビット並列 O(NM/64)）/ `lcs_restore(a, b)`（|a| × |b| が約 1.9×10^9 まで）/ `shortest_common_supersequence(a, b)`。**編集距離**：`edit_distance(a, b, ins, del, sub)`（費用は 0 以上。O(NM) 時間・O(M) メモリ）/ `edit_script(a, b, ins, del, sub)`（`EditOp{type, i, j}`。(N+1)(M+1) ≤ 2×10^8）/ `edit_distance_bounded(a, b, k)`（O((N + M) k)）/ `within_one_edit(a, b)` | str_base | EDPC F、ABC185 E、ABC324 C、ABC386 C・F |
| ★ `06_lis` | **作成済み（Phase 3）**。vector<整数>・string。strict = true で狭義、false で広義：`lis_length(a, strict)` / `lis_restore(a, strict)`（添字の列）/ `lis_ending_at(a, strict)` / `lis_starting_at(a, strict)` / `lis_members(a, strict)`（`LisMembers{length, in_some, in_all}`）/ `count_lis<T>(a, strict)` / `lds_length(a, strict)` / `min_increasing_cover(a, strict)` / `min_decreasing_cover(a, strict)` / `lis_2d(ps)`（x も y も < だけで比べる。double などでも） | str_base | ABC038 D、ABC134 E、ABC237 F（△）、ABC354 F、ABC360 G |

### 7. 括弧列

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `07_bracket` | **作成済み（Phase 3）**。括弧の種類は文字列 pairs で渡す（既定 "()"）。括弧以外の文字は無視する：`is_balanced(s, pairs)` / `match_brackets(s, pairs)` / `bracket_depths(s, pairs)` / `bracket_parents(s, pairs)`（囲む一番内側の組の開き括弧の位置）/ `erase_matched_pairs(s, pairs)`、`count_balanced_substrings(s, pairs)` / `longest_balanced_substring(s, pairs)`、`StaticBracket sb(s, pairs)`（`is_balanced(l, r)` O(1)・`partner(i)`。**何種類でも**）。1 種類の括弧だけ：`BracketBalance{need_open, need_close}`（`+` / `+=` / `balanced()` / `net()`）と `balance_of(s)` / `balance_of(c)`、`bracket_concat_order(vs)`（`BracketOrder{ok, order}`）/ `min_insertions_to_balance(s)` / `fix_brackets(s)`、`BracketSegTree bt(s)`（`set(i, c)` / `flip(i)` / `get(i)` / `is_balanced(l, r)` / `balance(l, r)` / `all_balanced()` / `to_string()`）、`for_each_balanced(n, f)` **最終照合で足した**：`count_balanced_completions<T>(s, wild, pairs)`。 | str_base | 典型90 002、ABC167 F、ABC223 F、ABC307 D、ABC350 F、ABC394 D、ABC428 C |

### 8. 文字列の集合（接頭辞）

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `08_trie` | **作成済み（Phase 3）**。`Trie TR(A)`（多重集合。子は 1 次元の配列）：`insert(s, cnt = 1)`（ノードを返す。cnt = 0 でノードだけ作る）/ `erase(s, cnt = 1)` / `erase_prefix(p)` / `clear()` / `count(s)` / `contains(s)` / `count_prefix(p)` / **`count_prefixes_of(s)`** / `size()` / `empty()` / `kth(k)` / `kth_ids(k)` / `rank(s)` / `words()` / `words_ids()` / `max_lcp(s, exclude_self)` / `sum_lcp_with(s)` / `sum_lcp_pairs()`。ノード：`find(s)` / **`path(s)`** / `child(v, c)` / `child_by_id(v, id)` / `parent(v)` / `depth(v)` / `label(v)` / `pass_count(v)` / `end_count(v)` / `string_of(v)` / `preorder()` / `node_count()` / `root()` / `alphabet()` | str_base | ABC268 G、ABC287 E、ABC353 E、ABC377 G、ABC403 E |

### 9. 操作・更新のある文字列

| ファイル | 責務・提供する機能（案） | 依存 | 代表的な問題 |
|---|---|---|---|
| ★ `09_reversible_deque` | **作成済み（Phase 3）**。`ReversibleDeque<T>`（T = char で文字列）：`push_front` / `push_back`（1 つ・列。自分の要素・自分自身も可）/ `pop_front` / `pop_back`（取り出して返す）/ `front` / `back` / `[]`（書き換えも）/ `reverse()`（O(1)）/ `reversed()` / `size` / `empty` / `clear` / `to_vector()` / `to_string()` / `substr(pos, len)`。`RotatableArray<T>`（大きさは固定）：`rotate_left(k)` / `rotate_right(k)` / `reverse()`（どれも O(1)）/ `[]` / `get(i)` / `set(i, v)` / `size()` / `to_vector()` / `to_string()` **最終照合で足した**：`reserve(n)`。 | str_base | ABC158 D、ABC199 C、ABC237 D |
| ☆ `09_mutable_string` | **作成済み（Phase 4）**。1 文字の書き換えがある文字列（ハッシュは 01_rolling_hash と同じ基数なので、ほかの文字列とも比べられる）。`MutableString ms(s, A)`：`set(i, c)` / `get(i)`（値。char は 0〜255）/ `size()` / `hash(l, r)` / **`reversed_hash(l, r)`**（`StrHash`）/ `same(l1, r1, l2, r2)` / `is_palindrome(l, r)`（O(log N)）/ `lcp(i, j)` / `compare(l1, r1, l2, r2)`（O(log² N)）/ `count(l, r, c)` / `counts(l, r)` / `is_sorted(l, r)`（Alphabet の順で広義単調増加か）/ `kth_char(l, r, k)`（区間を並べ替えた k 番目の値）/ `values()` / `to_string()` / `alphabet()` | str_base, rolling_hash | ABC285 F、ABC331 F |

### △（型）

| 型 | 置き場所 |
|---|---|
| 巨大な文字列の k 文字目（倍々に伸びる・置き換えを繰り返す） | 00_str_base の使い方関数（ABC380 D・ABC242 D の形） |
| パターンを含まない文字列の数え上げ（DP） | 02_pattern_match の `KmpAutomaton` の使い方関数（パターンが多ければ 02_aho_corasick、長さが 10^18 なら math_v2 05_matrix で累乗） |
| 書き換えのたびに「ABC」の数を数え直す | 書き換えた位置の前後だけ数え直す（02_pattern_match の使い方関数。ABC372 C の形） |
| 括弧列の数え上げ DP | 07_bracket の使い方関数（ABC312 D の形。DP の本体は DP 側） |
| 回文の条件から文字列を作る | 04_palindrome の使い方関数（ABC349 G の形。条件をまとめるのは graph_grid_v2 03_union_find） |
| 文字列を頂点にした最短路 | graph_grid_v2（00_state_id ＋ 01_bfs） |

### 対象外（決定済み）

- 数の文字列・回文数・大きな数の余り・部分文字列の値・数をつなげる（math_v2 01_digits）
- 挿入・削除・区間の反転・区間の代入がある文字列（DynamicString。データ構造）
- WaveletTree（文字列も可）・static_seq_toolkit の文字列の機能（データ構造。体系の外に今のまま残す）
- BinaryTrie（XOR の最小・最大。データ構造）
- HashMap / HashSet のラッパ（データ構造）
- 入力が文字列でも、本質が尺取り・累積和・一般の DP の問題（それぞれのテンプレート）
- 小さいグリッドの模様探し（graph_grid_v2 10_grid_transform の `find_subgrid`。大きいときは 01_rolling_hash の `RollingHash2D`）

## 2-P2. Phase 2（基盤）の実装メモ（2026-10-12）

### 作ったファイル（`Desktop\claude\string_v2\`）

| ファイル | 内容 | テスト（主な愚直解。Release の照合回数） |
|---|---|---|
| `00_str_base.cpp` | 文字列基盤（Alphabet・CharPrefix・CharPositions・ランレングス・変換・CharMap・比較・区切り・StrSpan） | 全 256 文字で Alphabet の番号・値・入っているかを map と比べる（of はランダムな並びで 2400 通り、ints は 10^18 級・負の範囲と境目の外側）、区間の文字の数・次 / 前の出現・k 番目・出現位置を全部の区間・全部の位置で愚直に数えて比べる（英小文字・括弧・128 以上の char・負の値の整数の列）、ランレングスは戻すと元に戻る・隣の run の値が違う、CharMap は replace を重ねたものと文字列を順に書き換えたものを比べる、部分列か は全部の部分集合で、split は手で区切ったものと（177 万回） |
| `update_bundles.ps1` / `.bat` | graph_grid_v2・math_v2 のものをそのままコピー | － |
| `dynalist_items_phase2.txt` | Dynalist に貼る項目（1 本） | － |
| `_work\` | 組み立て用の作業ファイル（本体・テスト・デモの元と build.py、問題の確認 problems\、名前の衝突の確認 tools\） | － |

確認したこと：g++ 13（-O2 -Wall -Wextra -Wshadow -pedantic）と clang 18（-pedantic-errors）で C++17 / 20 / 23 とも警告なし・全テスト通過・デモの出力が同じ、`-D_GLIBCXX_DEBUG -fsanitize=address,undefined,float-cast-overflow`（g++）と `-fsanitize=undefined,integer`（clang）で全体（テスト・速度計測・デモ）が通過、使えない引数で止まること（範囲外の文字・同じ文字が 2 回の of・範囲外の位置・l > r・長さの違う hamming_distance・負の k など 21 通り。どれもメッセージを出して abort）、Main_base.cpp（ACL 込み）に**グラフ・グリッド 46 本・幾何 17 本・数学 26 本・文字列 00_str_base** の本体を 2 通りの順で全部貼り、gnu++2b / c++17・LOCAL あり / なしで g++・clang の両方でコンパイル・実行（名前の衝突なし。ACL の `suffix_array` / `lcp_array` / `z_algorithm`、数学の `is_palindrome_str`、グラフの `trim` もそのまま呼べる。`debug(runs)`・`debug(StrSpan{…})` も可）、Main_base に 1 本だけ貼って -Wshadow の警告なし（g++・clang）、ほかの単独テンプレート（binary_search.cpp・static_seq_toolkit.cpp・dynamic_unique_seq.cpp・segtree・lazy_segtree）と一緒に貼ってコンパイルできる（static_seq_toolkit の `is_palindrome` / `run_length` はメンバなので重ならない）、代表的な問題 8 問（ABC219 C・ABC232 B・ABC259 C・ABC310 C・ABC329 C・ABC342 C・ABC380 D・典型90 084）を問題文と入出力例を確かめて Main_base ＋ 本体だけで解き、例 22 個がすべて合う。PC（`Desktop\claude\string_v2`、g++ 11.4）に保存して、ビルド（警告 0）・ランダムテストが通ることを確かめた。

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）

| 処理 | 規模 | 時間 |
|---|---|---|
| CharPrefix の構築 / count / next + prev | N = 10^6, σ = 26 / 10^7 回 / 10^6 回 | 63〜75 / 300〜380 / 40〜53 ms |
| CharPositions の構築 / next + prev | N = 10^6 / 10^6 回（N = 2×10^5 でも） | 9 / 280〜290 ms（約 200 ms） |
| run_length / to_ids / CharMap::apply | N = 10^7 | 120〜140 / 26 / 21 ms |
| compress_values | N = 10^6（long long） | 約 125 ms |

旧テンプレートとの比較（同じ計測環境。英大文字・小文字の 52 種。ms）：

| 処理 | N = 2×10^5 | N = 10^6 |
|---|---|---|
| 旧 StringHelper の構築 / next + prev × 10^6 | 56 / 25 | 280 / 42 |
| **CharPrefix** の構築 / next + prev × 10^6 | 24 / 20〜25 | 126 / 32〜40 |
| **CharPositions** の構築 / next + prev × 10^6 | 2 / 170 | 9 / 220 |

- CharPrefix は旧 StringHelper（nxt と prv の 2 つの表を vector<vector<int>> で持つ）と同じかそれより速く、メモリは約半分（(N+1)σ 個の int ＋ N）。旧版の「なし」は next が N、ここでは -1
- CharPositions の二分探索は std::lower_bound と同じ速さ（この環境の 1 回の探索は約 100 ns）。O(1) が要るときは CharPrefix（06_subsequence の NextTable も O(1) の予定）
- ランレングスは旧 RunLengthTemplateTest の run_length とほぼ同じ速さ（N = 10^7 で 36〜55 ms。開始位置も持つ）

### 実装で決めたこと

- **Alphabet は実行時の値**（6-1 の確認 3）。中身は 256 個の表（文字）か、lo と大きさ（整数）。`id` は表を 1 回引くだけ。`of` で並べた順が番号の順（= 辞書順）なので、独自の順の比較（ABC219 C）にも使える。設計になかった `bytes()`（0〜255。日本語の UTF-8 の文字列にも）を足した。ints は大きさ 10^9 まで（大きな値は `compress_values` で 0〜K-1 にしてから `ints(K)`）
- 範囲外の文字・位置は `str_internal::fail` で関数名と理由（どの文字か・どの範囲か）を出して abort する（数学と同じ。assert と違い NDEBUG でも止まる）
- **CharPrefix に next / prev / kth / total（O(1)）を足した**（設計では CharPositions だけが持つ予定だった。区間の数の表と位置の一覧があれば next(i, c) = 「c の (i より前の数) 番目」で O(1) になるため。迷ったら入れる）。CharPositions はメモリを抑えた版（O(N + σ)・O(log N)）として残した。(N+1)σ が 2×10^8 を超える CharPrefix は作らずに止める（Phase 4 で 6×10^8 から下げた。AtCoder の 1024MB に収めるため）
- 列は template で受け取り、`str_internal::seq` で文字列リテラル（"abc"）を std::string にする。要素は `code_of` で値にする（char は unsigned char として 0〜255。std::string の比較の順と同じ）。string と vector<int> を混ぜて比べてもよい（common_prefix_length(vector<int>, string) など）
- `RunBlock{value, len, start}` の len と start は long long（自分で作った巨大な run を `run_length_decode` に渡す場合も考えて）。decode は元の長さが 2×10^9 を超えると止める
- `caesar_shift` と `CharMap` は Alphabet にない文字をそのまま残す（空白・大文字などを含む文字列にそのまま使えるように。ほかの関数と違い止めない）
- `join` は string・char・整数・`<<` で出せるもの（数学の ModInt など）を受け取る（答えの空白区切りの出力に）
- `StrSpan` と `RunBlock` に `<<`（`[l, r)`・`{value, len, start}`）を付けた（Main_base の `debug` で見られるように。friend 関数なので名前は増えない）
- 二分探索は自分で書いた（分岐の少ない形。std::lower_bound のデバッグビルドでの検査を避ける）
- テストとデモの形はほかの第二世代と同じ。テスト用の文字列の作り方（一様・同じ文字の連続・周期・フィボナッチ・Thue–Morse・ほぼ同じ文字）を `_work\src\common_test.inc` に足した（Phase 3 以降も使う）
- デモの例：ABC219 C（Alphabet::of）、ABC259 C・典型90 084・ABC329 C（run_length）、ABC232 B（caesar_shift）、ABC342 C（CharMap。入力例 1）、ABC310 C（min_with_reverse）、ABC380 D（巨大な文字列の k 文字目。△ の型）

### 次回以降の作業のしかた

- 本体の正は各 .cpp の「ここからテンプレート本体（コピー範囲）」の中。基盤を直したら `update_bundles.bat` を実行する（Phase 3 以降、すべてのファイルが 00_str_base を同梱する）
- 組み立ては `_work\build.py`（src の head / body / rest から out に .cpp を作る）。名前の衝突は `_work\tools\check_all.py`（Main_base ＋ 全部の第二世代）と `check_singles.py`（単独テンプレート）、問題の確認は `_work\problems\run.py`

## 2-P3. Phase 3（★ 11 本）の実装メモ（2026-10-12）

### 作ったファイル（`Desktop\claude\string_v2\`）

| ファイル | テスト（主な愚直解。Release の照合回数） |
|---|---|
| `01_rolling_hash.cpp` | 全部の部分文字列の組で「ハッシュが等しい ⇔ 文字列が等しい」、lcp・common_suffix・compare を愚直に、連結・repeat を文字列を作って、2 次元は全部の長方形と模様探し（string・負と 10^18 級の数の列・混ぜたもの。619 万回） |
| `02_pattern_match.cpp` | std::string::find を全部の位置で、Z・π・周期・境界・回転・無限の繰り返し・重なりを O(N²) の定義どおりに、KmpAutomaton は全部の状態と文字で（11 万回） |
| `03_suffix_array.cpp` | 全部の接尾辞を並べたもの、全部の部分文字列の set（異なる数・k 番目・重複ありの k 番目・順位・出現回数）、LCP は全部の組、of_strings は区切りをまたがないこと、最長共通部分文字列は DP と（N < 400 で SA-IS の再帰も。25 万回） |
| `04_palindrome.cpp` | 全部の区間の回文判定、中心ごとの長さ、最長・個数・接頭辞・接尾辞、最短の回文（足す文字数を全部試す）、最長回文部分列・回文分割・回文の部分列の数は全部の部分列・O(N³) の DP と（4 万回） |
| `05_lex_order.cpp` | concat_order は全部の並べ方、回転は全部の回転、bijective・dict_rank / unrank・next_string は全部の文字列を並べたもの（独自の文字の順・ints も。41 万回） |
| `06_subsequence.cpp` | 全部の部分列の列挙（部分列か・異なる部分列の数・長さごと・部分列でない最短（BFS））、next / prev は全部の (i, c)、出現の数は O(NM) の DP と（1.2 万回） |
| `06_lcs_edit.cpp` | O(NM) の DP（LCS・重み付きの編集距離）、操作の列は当てはめて b になり費用が最小、k 以内は全部の k で（64 の倍数の前後の長さも。4.8 万回） |
| `06_lis.cpp` | 全部の部分集合の列挙（長さ・数・入りうる / 必ず入る位置）、O(N²) の DP（各位置で終わる / 始まる・数）、分割の最小数はバックトラック、2 次元は O(N²)（4.8 万回） |
| `07_bracket.cpp` | 「隣り合う組を消す」の愚直解（対応・正しいか・深さ・親）、全部の区間（StaticBracket・数・最長）、最小の挿入は区間 DP、concat は全部の並べ方、セグメント木はランダムな操作、全列挙は 2^n 通りから（5 万回） |
| `08_trie.cpp` | map で持った多重集合と、ランダムな insert / erase / erase_prefix のたびに全部の問い合わせを（英小文字・独自の順・負の整数・0/1。14 万回） |
| `09_reversible_deque.cpp` | std::deque（反転は std::reverse）・std::rotate と、ランダムな操作のたびに（1.3 万回） |
| `dynalist_items_phase3.txt` | Dynalist に貼る項目（11 本） |
| `dynalist_tool_index.txt` / `.opml` | 索引（4 章と同じ。☆ は Phase 4 で埋める） |

確認したこと：g++ 13 と clang 18 で C++17 / 20 / 23、-Wall -Wextra -Wshadow -pedantic（と -pedantic-errors）で警告なし。`-D_GLIBCXX_DEBUG -fsanitize=address,undefined`（g++）と `-fsanitize=undefined -fsanitize-trap`（clang）でテスト・速度計測・デモ全体が通過。使えない引数で止まること 38 通り（範囲外の位置・区間・Alphabet にない文字・不正な pairs・答えが long long に収まらない数え上げ・大きすぎる edit_script など。どれも関数名と理由を出して abort）。Main_base（ACL 込み）にグラフ・グリッド 46 本・幾何 17 本・数学 26 本・文字列 12 本の本体を 2 通りの順で全部貼り、gnu++2b / c++17・LOCAL あり / なし・g++ / clang で全部コンパイル・実行（ACL の suffix_array・lcp_array・z_algorithm と並べて呼べる。よく使う変数名 rank・hash・next・prev・count・runs・periods・lcs・sa・lcp・rev とも衝突しない）。1 本ずつ Main_base に貼って -Wshadow の警告なし。単独テンプレート（binary_search・static_seq_toolkit・dynamic_unique_seq・segtree・lazy_segtree）と一緒に貼ってコンパイルできる。グローバルな名前 143 個はほかの第二世代の約 900 個と重ならない。PC（`Desktop\claude\string_v2`、g++ 11.4）に保存して、12 本ともビルド（-Wall -Wextra -Wshadow で警告 0）・ランダムテストが通り、`_work\build.py` で組み立て直すと保存した .cpp と同じになることを確かめた。

### 代表的な問題（問題文と入出力例を確かめ、Main_base ＋ 本体だけで解いた。51 問・例 145 個がすべて合う）

ABC141 E（ハッシュと接尾辞配列の 2 通り）、ABC284 F、典型90 047、ABC257 G、ABC150 F、ABC135 F、ABC430 E、ABC343 G、ABC213 F、ABC272 F、ABC280 Ex、ABC237 C、ABC320 B、ABC398 F、ABC225 F、ABC219 C、ABC171 C、ABC285 C、ABC223 B、典型90 006・008、ABC211 C、ABC214 F、ABC230 F、ARC081 E、ABC324 E、EDPC F、ABC185 E、ABC324 C、ABC386 C・F、ABC038 D、ABC134 E、ABC354 F、ABC360 G、典型90 002、ABC167 F、ABC223 F、ABC307 D、ABC350 F、ABC394 D、ABC428 C、ABC268 G、ABC287 E、ABC353 E、ABC377 G、ABC403 E、ABC158 D、ABC199 C、ABC237 D。答えが複数ある問題（ABC284 F・EDPC F・ABC280 Ex）は答えが正しいかを確かめる関数で照合した。`_work\problems\phase3.py`。

- ABC240 Ex（接尾辞配列の上の DP。5-5 では 03_suffix_array の例）は解いていない（テンプレートの範囲ではなく難しい DP の問題。最終照合で扱いを決める）

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い。`_work\bench\cmp3.cpp`）

ACL・旧テンプレートとの比較（同じプログラムの中で 3 回の最小。ms）：

| 処理 | 規模 | 比べたもの | 新 |
|---|---|---|---|
| 接尾辞配列 | N = 10^6（σ = 26 / 2 / 1） | ACL suffix_array 84 / 79 / 19 | make_suffix_array 86 / 80 / 18 |
| 接尾辞配列（long long の列） | N = 10^6 | ACL 256 | 230 |
| LCP 配列 | N = 10^6（σ = 26） | ACL lcp_array 29 | make_lcp_array 30 |
| Z | N = 10^7（σ = 26 / 2 / 1） | ACL z_algorithm 89 / 118 / 95 | z_array 37 / 85 / 40 |
| ハッシュの構築 / get × 10^7 / lcp × 10^6 | N = 10^7 / 10^7 / 10^6 | 旧 RollingHashTest2 84 / 428 / 1085 | 76 / 375 / 273 |
| 回文の構築 / is_palindrome × 10^7 | N = 10^7 | 旧 ManacherTemplate 299 / 188 | 337 / 168 |
| Trie insert + 検索 / kth × 10^5 | 文字の合計 10^6 | 旧 ALDS14C 123 / 92 | 138 / 117 |
| LIS の復元 | N = 10^6 | 旧 A24 34 | 32 |

各ファイルの速度計測（Release。ms）：KMP の find_all（N = 10^7）118、match_lengths（10^7）43、SuffixArray の構築（10^6）177、lcp(i, j) × 10^6（表の構築込み）103、lcs_length（5×10^4 × 5×10^4）53、lcs_restore（3×10^4 × 3×10^4）103、edit_distance（10^4 × 10^4）243、edit_distance_bounded（5×10^5・k = 20）76、lis_members / count_lis / lis_2d（10^6）73 / 262 / 339、NextTable（10^6）+ match_length × 10^5 502、括弧の判定・対応・数（10^6）30、BracketSegTree（N = Q = 5×10^5）88、ReversibleDeque 10^7 回 140、RotatableArray 10^7 回 171。

- Trie は旧版より 1〜2 割遅い（Alphabet の確認・親と深さを持つぶん）。ほかは同じかそれより速い
- 速くするために変えたこと：Z・KMP は列を long long の列にコピーせずに読む（ACL の約 2 倍速い）。RollingHash は範囲の確認を軽くし（止めるときのメッセージ作りを別の関数に）、中では確認なしで引く（lcp が約 4 倍速い）。char の列では元の値を持たない（1 文字のハッシュ = 値 + 1 から戻す）。LIS の二分探索は「末尾より大きければすぐ足す」＋分岐ありの形（ばらばらな列で分岐なしより速かった）。SA-IS の L/S の型は vector<bool>（induce で飛び飛びに読むので小さい方が速い）

### 実装で決めたこと（設計からの変更を含む）

- 名前を変えたもの：`manacher_radii` → **`center_palindrome_lengths`**（半径ではなく長さを返すので）、`min_append_to_palindrome` / `min_prepend_to_palindrome` → **`min_palindrome_with_prefix` / `min_palindrome_with_suffix`**（足す文字ではなく回文そのものを返す。ABC398 F をそのまま出力できる）、`bijective_to_string(n, k, first)` → **`bijective_to_string(n, A)`**（文字の種類は Alphabet）
- 足したもの：`RollingHash::common_suffix`・別の文字列との compare、`Palindromes::length_at_center` / `span_at_center`、`NextTable::char_prefix()`、`shortest_non_subsequence_ids`、`LisMembers`、`StaticBracket::partner`、`BracketBalance::net()`、`Trie::count_prefixes_of`（ABC268 G）・`path(s)`（ノードに値を持たせる DP 用。ABC377 G・ABC403 E）・`child_by_id`・`label`・`string_of`・`clear`、`ReversibleDeque::reversed()`、`RotatableArray::get`
- NextTable は 00_str_base の CharPrefix を包むだけ（次 / 前の出現は CharPrefix が O(1) で持つようになったため）
- 00_str_base の中の補助（`str_internal`）を足した：`lower_index` / `upper_index`（自前の二分探索）、`checked_add`（整数型なら溢れたら止める）、`add_count`（途中で溢れても答えに使わなければ止めない）、`codes_of`、`lex_compare`、`range_fail`（範囲の確認を軽くする）。同じ補助をファイルごとに書くと一緒に貼ったときに重なるので、2 本以上で使うものは基盤に置く
- 数え上げの型 T：long long を渡したとき、答えが収まらなければ止める。途中の値だけが溢れる場合（count_lis のフェニック木の和・count_subsequence_occurrences の途中の dp）は止めない
- 括弧列の対応は「隣り合う同じ種類の開き・閉じを消すことを繰り返す」で一緒に消える組。種類の違う閉じ括弧が来たら、それより前の開き括弧は後ろと対応しない。区間 [l, r) が正しいか は「境目 l と r の状態の番号が等しいか」（開きを積むたびに新しい番号）で、何種類でも O(1)（設計の疎テーブル 2 つより簡単で速い）
- BracketBalance・concat・最小の挿入・BracketSegTree・全列挙は 1 種類の括弧だけ（何種類もの最小の挿入は O(N³) の区間 DP になるので作らない）
- SuffixArray で空の pat の出現は「置ける位置すべて」（1 つの文字列なら 0〜N の N + 1 個。of_strings では各文字列の 0〜長さ、合計 N 個）。02_pattern_match・01_rolling_hash と同じ考え方
- 編集距離の費用は 0 以上（負は止める）。足し算は上で止める（とても大きな費用でも溢れない）
- concat_order は O(N log N × 最長の長さ)。とても長い文字列が 1 本混ざると遅い（ABC225 F のような N・長さが小さい問題向け）

### 独立したレビュー（作ったものを見ていない 3 つのエージェントに、本体を読ませて愚直解で確かめさせた）

正しい入力で結果が間違うバグは 1 件（ReversibleDeque）。ほかは極端な入力・説明とのずれ。すべて直してテストを足した：

| ファイル | 指摘 | 直したこと |
|---|---|---|
| 09_reversible_deque | `d.push_back(d.front())` で、要素数が 16・32・… のときに置き場が広がり、解放した場所を読む（値が壊れる） | 足す値を先に写してから広げる。`d.push_back(d)`（自分自身）が止まらないのも直した |
| 06_lis | `lis_2d` が y を long long にしていた（double・2^63 以上の unsigned で誤る） | x・y とも `<` だけで比べ、y は順位にしてから LIS |
| 06_lis | `count_lis<long long>` が、答えは収まるのに途中の和の溢れで止まる | 溢れた印を持ち、答えに使われたときだけ止める（`add_count`） |
| 06_subsequence | `count_subsequence_occurrences<long long>` も同じ | 同じ |
| 04_palindrome | `count_palindromic_subsequences<long long>` が A + B − C の A + B で溢れる | A + (B − C) の順に（B ≥ C） |
| 06_lcs_edit | 大きな費用（10^18）で編集距離の足し算が溢れる。`edit_distance_bounded` に大きな k を渡すと k に比例するメモリ | 費用は 0 以上・足し算は上で止める。k は max(|a|, |b|) までに |
| 06_lcs_edit | `shortest_common_supersequence(string, vector<int>{1000})` が値を char に切り詰める | 収まらなければ止める。lcs_restore の上限の説明を「両方 4×10^4 くらい」に |
| 08_trie | `sum_lcp_pairs` の pass × (pass − 1) が数 3×10^9 で溢れる | 先に偶数の方を 2 で割り、溢れたら止める |
| 07_bracket | 引数なしで作った StaticBracket・BracketSegTree の判定で配列の外を読む。BracketSegTree に 256 以上の値を渡すと切り詰める | 空の文字列と同じにする。値が 0〜255 の外なら止める |
| 03_suffix_array | `offset(p)` を区切りの位置で呼ぶと配列の外を読む。owner・length_at・value・position に範囲の確認がない。空の pat の数え方がほかのファイルと違う | 確認を足し、区切りの offset は止める。空の pat は「置ける位置すべて」に |
| 01_rolling_hash | 引数なしで作った RollingHash・RollingHash2D の get で落ちる | 空の文字列と同じにする |
| 説明 | rotation_offsets("", "") が {0}・CommonSubstring の pos_a < pos_b・count_distinct_subsequences の計算量・concat_order の計算量 | 説明を直した |

### 次回以降の作業のしかた

- 2-P2 と同じ（本体の正は各 .cpp の中。組み立ては `_work\build.py`）。Phase 3 で足した確認の道具：`_work\tools\death.py`（止まることの確認）、`_work\tools\solv_check.cpp`（check_all.py で全部を貼ったときに呼ぶ solv）、`_work\bench\cmp3.cpp`（ACL・旧テンプレートとの速さの比較。旧テンプレートの部分を `C:\atcoder_template` から抜き出して作る。作り直しは `gen.py`）、`_work\problems\phase3.py`（代表的な問題）
- 基盤（00_str_base）を Phase 3 で直したので、PC に保存したあと `update_bundles.bat` を実行しなくても、今回の 12 本は同じ基盤を同梱している

## 2-P4. Phase 4（☆ 7 本）の実装メモ（2026-10-12）

### 作ったファイル（`Desktop\claude\string_v2\`）

| ファイル | テスト（主な愚直解。Release の照合回数） |
|---|---|
| `01_set_hash.cpp` | std::set・並べた列と比べて「ハッシュが等しい ⇔ 集合（多重集合）が等しい」を全部の区間の組で、SetHash・MultisetHash はランダムな add / remove / merge のたびに（string・負と 10^18 級の数の列・大きな個数と負の個数。3.2 万回） |
| `02_aho_corasick.cpp` | 全部の位置で全部のパターンを比べた出現（回数・最初・合計・for_each_match の順番）、t の各接頭辞を読んだ状態（match_count・match_ids・longest_match・depth）、bfs_order とリンクの深さ、包含関係は全部の組の find（同じ文字列・接頭辞・接尾辞のパターン、σ = 1〜3。1.6 万回） |
| `02_wildcard_match.cpp` | O(NM) の愚直解（wild なし・pat だけ・両方、数の列の wild、LLONG_MIN・LLONG_MAX の値、128 以上の char の wild、pat にしかない値、NTT の長さが変わる所。2.8 万回） |
| `03_suffix_automaton.cpp` | 全部の部分文字列の map（異なる数・長さの合計・重複ありの数・出現回数・最初の位置・k 番目（重複なし・あり）・最長共通部分文字列・状態ごとの回数と len / link）、1 文字ずつ足したものと同じ、複数の文字列、extend・add_string と問い合わせを混ぜる、数の列（2.3 万回） |
| `04_palindromic_tree.cpp` | 全部の部分文字列から集めた回文（種類・出現回数・最初の位置・link・parent・next・各接頭辞の最長の回文接尾辞）、分割は O(N²) の DP（最小・数・接頭辞ごと・mod）、問い合わせと add を混ぜる（2.1 万回） |
| `05_lyndon.cpp` | 定義どおりの Lyndon 語の判定（分解・各接頭辞の最小の接尾辞・Lyndon 配列）、runs は全部の区間の最小周期から、平方は全部の (i, L)（数の列・128 以上の char も。2.8 万回） |
| `09_mutable_string.cpp` | ランダムな set のたびに全部の問い合わせ（同じか・回文・lcp・compare・数・昇順か・k 番目）を愚直に（独自の文字の順・表を作る前の書き換えも）、ハッシュは RollingHash と同じ値か、大きな値の列はハッシュだけ（5600 回） |
| `dynalist_items_phase4.txt` | Dynalist に貼る項目（7 本） |
| `dynalist_tool_index.txt` / `.opml` | 索引（☆ の部分を埋めて完成。4 章と同じ） |

確認したこと：2-P3 と同じ（g++ 13 と clang 18 で C++17 / 20 / 23、-Wall -Wextra -Wshadow -pedantic（と -pedantic-errors）で警告なし。`-D_GLIBCXX_DEBUG -fsanitize=address,undefined`（g++）と `-fsanitize=undefined -fsanitize-trap`（clang）でテスト・速度計測・デモ全体が通過）。使えない引数で止まること 69 通り（Phase 4 で 31 通り足した：build の前後の使い方・空のパターン・範囲外の状態 / ノード / 番号 / k・Alphabet にない文字・NTT の長さの上限・表が 2×10^8 個を超える・答えが long long に収まらない分割の数など）。`static_assert` で止まること（数の列で wild を省いた wildcard_match、`ms.set(1, "c")` のような文字列リテラル、小数）。Main_base（ACL 込み）にグラフ・グリッド 46 本・幾何 17 本・数学 26 本・文字列 19 本の本体を 2 通りの順で全部貼り、gnu++2b / c++17・LOCAL あり / なし・g++ / clang で全部コンパイル・実行（よく使う変数名 link・len・occ・cnt・fail・go・wild・lyndon・word・period・node・state とも衝突しない）。1 本ずつ Main_base に貼って -Wshadow の警告なし。単独テンプレートと一緒に貼ってコンパイルできる。グローバルな名前は 170 個になり、ほかの第二世代の約 900 個と重ならない。PC（`Desktop\claude\string_v2`、g++ 11.4）に 19 本とも保存し直し、ビルド（-Wall -Wextra -Wshadow で警告 0）・ランダムテストが通り、`_work\build.py` で組み立て直すと保存した .cpp と同じになり、`phase4.py` の問題も PC で全部合うことを確かめた。

### 代表的な問題（問題文と入出力例を確かめ、Main_base ＋ 本体だけで解いた。14 通り・例 30 個がすべて合う）

ABC250 E、ABC367 F（01_set_hash）、ABC268 Ex、ABC354 G（ACL の mf_graph で最大流）、ABC362 G（02_aho_corasick）、ABC196 F（02_wildcard_match）、ABC362 G（03_suffix_automaton でも）、ABC237 Ex（04_palindromic_tree ＋ mf_graph で二部マッチング）、ABC285 F、ABC331 F（09_mutable_string）。Library Checker の Number of Substrings・Eertree・Lyndon Factorization・Enumerate Runs（入力例は公式の gen/example、出力は愚直解で作って照合）。`_work\problems\phase4.py`。ABC の 8 問は最大の大きさのランダムな入力でも 0.4 秒以内（この計測環境）。

- ABC280 Ex を接尾辞オートマトンの kth_substring で解くと、答えが長い問い合わせが多いと遅い（|答え| に比例）。03_suffix_array の kth_substring_with_dup（Phase 3 で解いた）を使う、とデモと索引に書いた

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い。Release。ms）

旧テンプレートに同じものがないので、各ファイルの速度計測だけ：

| 処理 | 規模 | 時間 |
|---|---|---|
| prefix_set_hashes / MultisetPrefix 作成 + get × 10^6 / SetHash 追加・削除 10^6 回 | N = 10^6 | 120 / 32〜60 / 180 |
| AhoCorasick build / count_each | パターンの長さの合計 10^6（σ = 4）・本文 10^6 | 67 / 46 |
| 同（σ = 26）build + count_each | 同 | 約 430（ノード 74 万） |
| mismatch_counts（σ = 2）/ wildcard_match（σ = 26）/ mismatch_counts（σ = 26） | 10^6 と 5×10^5 / 5×10^5 と 2×10^5 / 同 | 約 170〜200 / 約 170〜220 / 約 640〜700 |
| SuffixAutomaton 作成（σ = 26）/ 作成（σ = 2）/ count_occurrences × 10^5 | N = 5×10^5 / 10^6 | 100 / 130 / 20 |
| PalindromicTree 作成（σ = 2）/ min_partition / 全部 a の作成 + count_partitions | N = 10^6 | 20 / 22 / 115 |
| Lyndon 分解 + 各接頭辞の最小の接尾辞 / lyndon_array / maximal_repetitions | N = 10^6 | 25 / 110 / 700〜800 |
| MutableString 作成 / set + is_palindrome 10^6 回 / lcp × 10^5 | N = 10^6 | 37 / 190 / 190 |

メモリ（N = 10^6）：接尾辞オートマトン（σ = 26）約 170MB、回文木（全部 a・σ = 26）約 150MB、Aho–Corasick（ノード 70 万・σ = 26）約 115MB、maximal_repetitions 約 220MB。

### 実装で決めたこと（設計からの変更を含む）

- 名前・形を変えたもの：AhoCorasick の `first_end_each` → **`first_occurrence_each`**（終わりではなく開始位置。ほかのファイルの「位置」と同じ）、`for_each_match` の f は (開始位置, 番号)。SuffixAutomaton の `kth_substring` は string を返す（番号の列は `kth_substring_ids`。SuffixArray と違い位置を持たないため）。MutableString の `get` は値を返す。wildcard_match・mismatch_counts は wild を省く形と渡す形の 2 つ（数の列では wildcard_match の wild を必ず渡す）。`Repetition` の長さは `len()`（StrSpan と同じ）
- 足したもの：AhoCorasick の `next_id`・`link`・`bfs_order`・`pattern_count`・`node_of`・`pattern_length`、SuffixAutomaton の `kth_substring_ids`・`count_substrings_with_dup`・`next_id`・`first_end`・`order_by_len`・`last`、PalindromicTree の `next`・`span`・`size`・`length`・`min_partition_spans`・`min_partition_prefix`・`count_partitions_prefix<T>`、SetHash の `count`・`distinct`、MultisetPrefix の `same`・`size`、MutableString の `reversed_hash`・`values`・`to_string`・`alphabet`
- 05_lyndon は 03_suffix_array を同梱する（設計では 02_pattern_match）。Lyndon 配列は「次に小さい接尾辞」なので接尾辞配列の順位から O(N)、runs は 2 つの文字の順それぞれの Lyndon 根を前後の LCP（接尾辞配列 2 本の O(1) の LCP）で伸ばす。文字列は char のまま接尾辞配列に渡す（速い道）
- wildcard_match は乱数（Σ P T (P − T)^2 を 2 回）、mismatch_counts は正確（値ごとの畳み込みを周波数の側で足して逆変換は 1 回）。畳み込みの長さは text の長さ N 以上で足りる（読むのは M−1〜N−1 だけ）
- 接尾辞オートマトン：出現回数・k 番目は足した後の最初の呼び出しで数え直す。異なる部分文字列の数と長さの合計は足すたびに更新（複製では合計が変わらない）。add_string の後の extend はその文字列の続き。空文字列の数は Σ(|s_i| + 1)（SuffixArray::of_strings と同じ）。部分文字列の問い合わせの t は Alphabet の外の文字を含んでもよい（現れない文字として扱う）が、1 文字の遷移 next(v, c) は止める（数え上げの DP で番号を文字のつもりで渡す間違いを見つける）。Aho–Corasick の本文は Trie と同じく範囲外で止める
- 回文木の分割は直列リンク（diff・slink）で O(N log N)。数え上げの溢れは 06 系と同じく答えに使われたときだけ止める
- 子を配列で持つ構造（Aho–Corasick・接尾辞オートマトン・回文木）と表（CharPrefix・MutableString の文字の数）は、2×10^8 個（約 800MB）を超えるときはメッセージを出して止める（σ が大きい数の列で黙って MLE / bad_alloc にならないように）。列全体を渡すコンストラクタは先に領域を取る
- 集合のハッシュの unordered_map の鍵は、実行ごとの乱数を混ぜたハッシュにする（std::hash<long long> は値そのままで、狙った入力で遅くなる）。prefix_set_hashes は並べて最初の出現を探す（ハッシュ表を使わない）

### 独立したレビュー（作ったものを見ていない 3 つのエージェントに、本体を読ませて愚直解で確かめさせた。直した後にもう一度確かめさせた）

正しい入力で結果が間違うバグは 2 件（02_wildcard_match）。ほかは計算量・メモリの落とし穴と説明のずれ。すべて直してテストを足した：

| ファイル | 指摘 | 直したこと |
|---|---|---|
| 02_wildcard_match | 「wild なし」の目印に LLONG_MIN を使っていたので、値 LLONG_MIN どうしが一致にならない | 目印をやめ、wild を省く形と渡す形に分けた |
| 02_wildcard_match | 128 以上の char を wild に渡すと -1 になり、列の値（0〜255）と合わない | wild も code_of で値にする |
| 02_wildcard_match | 数の列でも既定の '?'（63）が wild になる | 数の列で wild を省くとコンパイルエラーに |
| 02_wildcard_match | pat にしかない値まで畳み込む・畳み込みの長さが倍 | 両方に現れる値だけ・長さ N 以上に（約 2 倍速い） |
| 01_set_hash | std::hash<long long> の unordered_set / map が狙った入力で 1000 倍遅くなる | 並べる形と、乱数を混ぜた鍵のハッシュに |
| 01_set_hash・09_mutable_string | same が長さの違いで範囲の確認をとばす。MultisetHash::remove の -k | 先に範囲を確かめる。remove は引き算で |
| 00_str_base（全部） | `ms.set(1, "c")` などポインタ・小数が黙って通る | code_of に static_assert |
| 02_aho_corasick | contained_patterns が呼ぶたびに O(パターン数)・リンクを何度もたどる（最悪 1.3 秒） | 印の番号を使い、一度見た状態の先は見ない（1〜10 ms） |
| 03_suffix_automaton | add_string の最後で last を根に戻していたので、last() と続きの extend が説明と違う。空文字列の数が of_strings と違う。next(v, 番号) が黙って -1 | 戻さない。Σ(|s_i| + 1) に。next は止める |
| 03_suffix_automaton | extend のたびに回数を聞くと O(N²)。デモの ABC280 Ex の解き方が本番の大きさで遅い | 説明に書き、異なる部分文字列の数は O(1) で更新に。ABC280 Ex は 03_suffix_array を案内 |
| 04_palindromic_tree | occurrences() が毎回写しを返す（ループの中で呼ぶと O(N²)）。各接頭辞の分割が取れない | 参照で返す。min_partition_prefix・count_partitions_prefix を足した |
| 共通 | σ が大きいとメッセージなしで bad_alloc / MLE。上限 6×10^8 は 1024MB を超える | 2×10^8 で止める。列を渡すコンストラクタは先に領域を取る |
| 05_lyndon | runs が文字列でも数の列の遅い道を通る。Lyndon 配列が要らない LCP も作る | char のまま渡す（約 1.4 倍速い）。順位だけ作る |
| 説明 | MutableString の set の計算量・lcp の計算量（設計書）・回文木の add の計算量・接尾辞オートマトンの各計算量・05_lyndon の依存（設計書） | 直した |

### 次回以降の作業のしかた

- 2-P3 と同じ。Phase 4 で足した確認の道具：`_work\problems\phase4.py`（代表的な問題）。`_work\tools\death.py` に Phase 4 の 31 通り、`_work\tools\solv_check.cpp` に Phase 4 の呼び出しを足した
- 基盤（00_str_base）を Phase 4 で直した（code_of の static_assert・CharPrefix の上限）ので、19 本とも新しい基盤を同梱し直した（PC には 19 本とも保存し直す）

---

## 3. 依存関係

```
00_str_base（文字列のすべてのファイルが同梱する）
├─ 01_rolling_hash
├─ 01_set_hash
├─ 02_pattern_match
├─ 02_aho_corasick
├─ 02_wildcard_match
├─ 03_suffix_array
├─ 03_suffix_automaton
├─ 04_palindrome
├─ 04_palindromic_tree
├─ 05_lex_order
├─ 05_lyndon            ← 03_suffix_array（Lyndon 配列と runs の LCP。Phase 4 で 02_pattern_match から変えた）
├─ 06_subsequence
├─ 06_lcs_edit
├─ 06_lis
├─ 07_bracket
├─ 08_trie
├─ 09_reversible_deque
└─ 09_mutable_string    ← 01_rolling_hash（同じ基数のハッシュ）
```

- ← の右に書いたファイル（と基盤）を同梱する
- グラフ・グリッド・幾何・数学には依存しない（同梱しない）。mod の数え上げは型を引数に取るので、math_v2 の 00_modint を一緒に貼ればそのまま使える。組み合わせ方は索引（4 章「ほかの第二世代と組み合わせる」）と使い方関数のコメントに書く

### 作る順番

- Phase 2（基盤）：00_str_base（＋ update_bundles のコピー、Dynalist の項目）
- Phase 3（★）：01_rolling_hash、02_pattern_match、03_suffix_array、04_palindrome、05_lex_order、06_subsequence、06_lcs_edit、06_lis、07_bracket、08_trie、09_reversible_deque（11 本）
- Phase 4（☆）：01_set_hash、02_aho_corasick、02_wildcard_match、03_suffix_automaton、04_palindromic_tree、05_lyndon、09_mutable_string（7 本）
- 最終照合：旧ファイルを 1 機能ずつ、過去の解答、AtCoder の典型（問題文を確かめながら）と広い洗い出し（`AtCoder_文字列問題一覧.txt` の 128 問を含む）、Dynalist の項目と突き合わせる

合計 19 ファイル（★ 12（00_str_base を含む）・☆ 7）。

---

## 4. 道具選択体系（Phase 3 で ★ の部分、Phase 4 で ☆ の部分を確定）

ほかの第二世代と同じく、Dynalist の「00_典型発想インデックス」の中に 1 項目（▼ 文字列道具選択（第二世代））として置き、テンプレートの項目は「文字列（第二世代）」に置く。辿り方は「問題の特徴（大分類）」→「聞かれていること（小分類）」→「ファイル名と関数」。`string_v2\dynalist_tool_index.txt` / `.opml` と同じ内容（Phase 4 で ☆ の関数名を入れて完成させた）。旧索引には文字列のキーワードがないので、旧「10_文字列」の項目は最終照合で新しい項目に置き換える。

```
▼ 文字列道具選択（第二世代）
	▼ 使い方
		1. 問題の特徴で大分類を選ぶ（文字・連続する文字／一致・出現／部分文字列の比較／回文／辞書順・回転／部分列／括弧／文字列の集合／操作・更新）
		2. 聞かれていることで小分類を選ぶ
		3. 書いてあるファイルを貼る（依存するファイルは同梱済み。1回のコピーで動く）
		範囲は半開区間 [l, r)（0 始まり）。「L 文字目から R 文字目」は (L - 1, R)。位置が「ない」ときは -1。k 番目は 0 から
		文字列も数の列（vector<int>・vector<long long>）もそのまま渡せる
		文字の種類は Alphabet で（既定は英小文字。Alphabet::upper()・Alphabet::of("()[]")・Alphabet::ints(K) など）
		記号・大文字が混ざる文字列は Alphabet::bytes()（番号 = 文字コード）か Alphabet::lower().with("@#")
		数（個数）は long long で返す（受け取る変数も ll に。chmax(int&, ll) はコンパイルエラー）。1 文字も長さ 1 の列として渡せる（find_all(S, 'a')）
		mod で数えるときは math_v2 の using mint = ModInt998;（10^9+7 なら ModInt107）を書いて count_xxx<mint>(…)
	▼ 文字・連続する文字　00_str_base
		文字の個数・区間 [l, r) の中の文字 c の個数 → count_chars ／ CharPrefix（count(l, r, c)。番号で回すなら count_id ／ total_id）
		区間の文字の数の差・偶奇が条件の部分文字列を数える → 先頭からの差（偶奇）を map で数える（ABC455 E・ABC418 D）
		位置 i 以降で最初の c・i より前で最後の c・k 番目の c → CharPrefix（next ／ prev。O(1)）／ CharPositions（next ／ prev ／ kth。省メモリ）
		連続する同じ文字をまとめる（ランレングス圧縮）・元に戻す → run_length ／ run_length_decode（ABC259 C・ABC329 C・典型90 084）
		「今ある c を全部 d に」を何度も → CharMap::replace ／ apply（ABC342 C）
		大文字・小文字・シーザー暗号 → to_lower ／ to_upper ／ swap_case ／ caesar_shift（1 文字の版も。ABC232 B）
		文字を番号に（'a' → 0）・独自の文字の順 → Alphabet::lower().id(c) ／ to_ids ／ Alphabet::of(順番)
		数の列を 0〜K−1 に → compress_values
		区切って分ける・つなげる（答えの出力） → split（1 文字・文字列で）／ join（vector・set など。区切りは文字列か 1 文字）
		前後の空白・特定の文字を除く → strip ／ strip_left ／ strip_right（ABC367 B）
		隣どうしの条件で区切る（隣が違う区間など） → split_runs_by（ABC456 C）
		隣が同じにならないように並べ替える・その辞書順最小 → rearrange_no_adjacent_equal ／ lex_smallest_no_adjacent_equal（ABC459 D）
		向きを区別しない（s と反転のうち小さい方） → min_with_reverse（ABC310 C）
		共通の接頭辞・接尾辞の長さ・違う位置の数 → common_prefix_length ／ common_suffix_length ／ hamming_distance
		巨大な文字列の k 文字目（倍々・置き換えの繰り返し） → 使い方デモ（ABC380 D・ABC242 D）
		S_n = S_{n-1} + S_{n-2} のような巨大な文字列の区間の文字の数 → 使い方デモ（長さと文字の数だけ持つ。ABC450 E）
	▼ 一致・出現（パターン検索）
		T が S のどこに現れるか・回数・最初・最後 → 02_pattern_match　find_all ／ count_occurrences ／ find_first ／ find_last
		重ならない出現の数・全部置き換える → count_non_overlapping ／ replace_all（1 回）
		置き換え・消すを現れなくなるまで繰り返す → rewrite_until_stable ／ erase_until_stable（ABC394 C・ABC328 D・ABC454 D）
		区間 [l, r) に T が現れるか・何回か（クエリ） → Occurrences::exists ／ count（ABC477 C）
		各位置から T と何文字一致するか（Z） → match_lengths ／ z_array（ABC257 G）
		数の列の一致（差や XOR の列にして） → find_all ／ rotation_offsets（ABC150 F）
		最小周期・s = t^k か・境界 → min_period ／ repetition_root ／ all_borders ／ all_periods ／ prefix_min_periods
		回転して一致するか・何回まわすか → rotation_offsets ／ is_rotation（ABC430 E）
		T を無限に繰り返した中での出現 → find_all_cyclic（ABC135 F）
		a の後ろと b の前が最大何文字重なるか → max_overlap（ABC343 G）
		T を含まない文字列を数える DP → KmpAutomaton
		多数のパターンの出現（それぞれの回数・最初の出現・どれかを含むか） → 02_aho_corasick　AhoCorasick::count_each ／ first_occurrence_each ／ contains_any（ABC362 G）
		多数のパターンを含まないようにする・DP → AhoCorasick::next ／ match_count（ABC268 Ex）
		パターンどうしの包含関係 → AhoCorasick::contained_patterns（ABC354 G）
		ワイルドカード（?）つきの一致 → 02_wildcard_match　wildcard_match
		? 込みで前から・後ろから何文字一致するか → 00_str_base　common_prefix_length(a, b, '?') ／ common_suffix_length（ABC287 D）
		ずらしごとの不一致の数・k 個以下の不一致での一致 → mismatch_counts ／ match_with_k_mismatches（ABC196 F）
		大きなグリッドの中で模様を探す → 01_rolling_hash　RollingHash2D::find_all（小さければ graph_grid_v2 find_subgrid）
	▼ 部分文字列の比較・同じ部分文字列
		2 つの部分文字列が等しいか（何度も） → 01_rolling_hash　RollingHash::get(l, r) ／ same（ABC284 F・典型90 047）
		長さを二分探索して同じ部分文字列を map で探す → RollingHash::get ＋ unordered_map<StrHash, …>（ABC141 E）
		部分文字列の LCP・辞書順の比較 → RollingHash::lcp ／ compare（O(log N)）、03_suffix_array　SuffixArray::lcp ／ compare（O(1)）
		つなげた・挿入した文字列のハッシュ → StrHash の + ／ repeat
		接尾辞配列そのもの（ACL と同じ） → make_suffix_array ／ make_lcp_array
		異なる部分文字列の数・辞書順 k 番目の部分文字列 → 03_suffix_array　count_distinct_substrings ／ kth_substring ／ kth_substring_with_dup（ABC280 Ex）
		2 回以上現れる最長の部分文字列（重なり可・不可） → longest_repeated ／ longest_repeated_non_overlapping（ABC141 E）
		2 つ・K 個の文字列の最長共通部分文字列（連続） → longest_common_substring ／ longest_common_substrings_all ／ longest_common_substring_k
		各接尾辞と全部の接尾辞の LCP の和 → sum_lcp_with_all（ABC213 F）
		パターンや部分文字列の出現回数を何度も → SuffixArray::count ／ count_occurrences(l, r)（ABC362 G）
		部分文字列を辞書順に処理する DP（辞書順に増える部分文字列の列など） → 03_suffix_array　rank・lcp_adj で（LCP が続く区間の左端, 長さ）を並べるキーに（ABC240 Ex）
		回転どうしの比較・複数の文字列をまとめて → make_suffix_array(S + S + '$' + T + T + '~')（ABC272 F）／ SuffixArray::of_strings（ABC280 Ex）
		部分文字列か・出現回数・1 文字足すたびに異なる部分文字列の数 → 03_suffix_automaton　SuffixAutomaton::contains ／ count_occurrences ／ count_distinct_substrings（ABC362 G）
		部分文字列を状態にした DP（出現位置の集合でまとめる） → SuffixAutomaton::len ／ link ／ next ／ endpos_count ／ order_by_len
		集合・多重集合が等しいか（並べ替えて同じか） → 01_set_hash　prefix_set_hashes ／ MultisetPrefix（ABC250 E・ABC367 F）
		追加・削除しながら集合・多重集合が等しいか → SetHash ／ MultisetHash
	▼ 回文
		区間 [l, r) が回文か（何度も） → 04_palindrome　Palindromes::is_palindrome
		最長の回文の部分文字列・回文の個数 → Palindromes::longest ／ count_all（ABC320 B）
		各位置を中心とする回文の長さ → radius_odd ／ radius_even ／ center_palindrome_lengths
		最長の回文の接頭辞・接尾辞・各位置から始まる最長 → longest_prefix_palindrome ／ longest_suffix_palindrome ／ longest_starting_at
		後ろ・前に足して回文にする最短 → min_palindrome_with_prefix ／ min_palindrome_with_suffix（ABC398 F）
		両端の文字を数えて間が回文か → Palindromes::is_palindrome（ABC237 C）
		回文にする最小の書き換え・並べ替えて回文にできるか → min_changes_to_palindrome ／ can_rearrange_to_palindrome
		最長回文部分列（飛び飛び）・回文の部分列の数 → longest_palindromic_subsequence ／ count_palindromic_subsequences<mint>
		k 文字まで書き換えてよい回文（中心ごとの最長・数） → 01_rolling_hash　center_palindrome_lengths_with_changes（ABC468 D）
		回転してから回文にする最小の書き換え → 回転ごとに min_changes_to_palindrome（ABC286 C）
		同じ長さで s 以下の回文の数（mod） → count_palindromes_leq<mint>(S, A)（ABC242 E）
		長さ 3 の回文の部分列（S_i = S_k）の数 → CharPrefix で真ん中ごとに 左の数 × 右の数（ABC375 D）
		回文に分割する最小の数・分け方の数 → palindrome_partition_min（N が大きいなら 04_palindromic_tree　PalindromicTree::min_partition ／ count_partitions<mint>）
		異なる回文の種類・それぞれの出現回数 → PalindromicTree::count_distinct ／ occurrences ／ count_total
		回文どうしの包含関係（含まない最大の選び方） → PalindromicTree::link ／ parent をたどる ＋ 二部マッチング（ABC237 Ex）
		書き換えがある回文判定 → 09_mutable_string　MutableString::is_palindrome（ABC331 F）
	▼ 辞書順・回転・名前の番号
		並べてつなげて辞書順最小 → 05_lex_order　concat_order ／ min_concat（ABC225 F）
		独自の文字の順で並べる → lex_less_by(order)（ABC219 C）
		a, …, z, aa, … の番号（表計算の列名） → bijective_to_string ／ bijective_from_string（ABC171 C・ABC285 C）
		同じ長さの文字列の中での辞書順の番号（mod） → lex_rank_same_length<mint>(S, A)
		同じ文字の個数が決まった並べ替えの辞書順 k 番目・何番目か → math_v2 04_permutation　multiset_perm_unrank ／ multiset_perm_rank（ABC202 D・ABC215 C）
		隣の入れ替えで T にする最小回数 → to_ids(S, Alphabet::of(T)) ＋ math_v2 04_permutation inversion_count（ABC264 D）
		長さ L 以下の全部の文字列の中での辞書順の番号・k 番目・次の文字列 → dict_rank ／ dict_unrank ／ next_string
		辞書順最小・最大の回転 → min_rotation ／ max_rotation（ABC223 B）
		Lyndon 分解・各接頭辞の最小の接尾辞・Lyndon 配列 → 05_lyndon　lyndon_factorization ／ min_suffix_of_prefixes ／ lyndon_array
		極大な繰り返し（runs）・xx の形の部分文字列の数 → maximal_repetitions ／ count_square_substrings
	▼ 部分列（飛び飛び）
		T が S の部分列か・先頭から何文字まで入るか → 06_subsequence　NextTable ／ prefix_match_lengths ／ suffix_match_lengths（ABC324 E）
		部分列として当てはめた位置 → 00_str_base　subsequence_positions(t, s)（ABC352 B。is_subsequence(t, s) と同じく拾う方 t が先）
		各開始位置から T を部分列として含む最短の終わり → NextTable::match_end（ABC452 D）
		次の文字 c の位置を何度も → NextTable::next ／ prev
		長さ K の辞書順最小・最大の部分列 → lex_smallest_subsequence ／ lex_largest_subsequence（典型90 006）
		部分列にならない最短の文字列 → shortest_non_subsequence（ARC081 E）
		T が部分列として何回現れるか → count_subsequence_occurrences<mint>（ABC211 C・典型90 008）
		異なる部分列の数（隣を選ばないなども） → count_distinct_subsequences<mint>（ABC214 F・ABC230 F）
		位置で区別して条件つきの部分列を数える（隣が同じでないなど） → 末尾の文字を状態にした DP（ABC456 D。count_distinct_subsequences は異なる文字列の数で別物）
		最長共通部分列（LCS）・復元・それを含む最短 → 06_lcs_edit　lcs_length ／ lcs_restore ／ shortest_common_supersequence（EDPC F）
		編集距離・操作の列 → edit_distance ／ edit_script（ABC185 E）
		1 回の操作で一致するか・K 回以内か → within_one_edit ／ edit_distance_bounded（ABC324 C・ABC386 C・ABC386 F）
		最長増加部分列（LIS）・復元・各位置で終わる / 始まる長さ → 06_lis　lis_length ／ lis_restore ／ lis_ending_at ／ lis_starting_at（ABC360 G）
		どれかの LIS・すべての LIS に入る位置 → lis_members（ABC354 F）
		LIS の数 → count_lis<mint>
		増加列（減少列）に分ける最小の数 → min_increasing_cover ／ min_decreasing_cover（ABC134 E）
		2 次元の LIS（箱の入れ子） → lis_2d（ABC038 D）
	▼ 括弧列　07_bracket
		正しい括弧列か（種類が複数でも） → is_balanced(s, "()[]<>")（ABC394 D）
		対応する括弧の位置・深さ・入れ子の木 → match_brackets ／ bracket_depths ／ bracket_parents（ABC350 F）
		対応する括弧を中身ごと消す → erase_matched_pairs（ABC307 D）
		1 文字ずつ足す・消すたびに正しいか → BracketBalance を積む（balance_of(c) を + する）（ABC428 C）
		括弧列を並べ替えてつなげる → bracket_concat_order（ABC167 F）
		最小の挿入で正しくする → min_insertions_to_balance ／ fix_brackets
		正しい括弧列の部分文字列の数・最長 → count_balanced_substrings ／ longest_balanced_substring
		区間が正しい括弧列か（書き換えなし・あり） → StaticBracket ／ BracketSegTree（ABC223 F）
		正しい括弧列をすべて辞書順に → for_each_balanced（典型90 002）
		? を開き・閉じに決めて正しい括弧列にする方法の数 → count_balanced_completions<mint>（ABC312 D）
		途中で閉じが開きより多くならないか（区間・書き換えあり。括弧でない 2 文字でも） → BracketSegTree(S, "AB").balance(l, r).need_open == 0（ABC473 F）
	▼ 文字列の集合（接頭辞）　08_trie
		p で始まる文字列の数・最長共通接頭辞 → Trie::count_prefix ／ max_lcp（ABC287 E）
		全部の組の LCP の和 → Trie::sum_lcp_pairs（ABC353 E）
		接頭辞である数・自分で始まる数（順位の期待値） → Trie::count_prefixes_of ／ count_prefix（ABC268 G）
		辞書順 k 番目・何番目か → Trie::kth ／ rank
		追加・削除しながら接頭辞を調べる・p で始まるものを全部消す → Trie::insert ／ erase ／ erase_prefix（ABC403 E）
		Trie の上で DP（各ノードの数・深さ） → find ／ path ／ child ／ pass_count ／ preorder（ABC377 G）
	▼ 操作・更新のある文字列
		前後に足す・反転を何度も → 09_reversible_deque　ReversibleDeque（ABC158 D・ABC237 D）
		回転（前半と後半の入れ替え）・反転を何度も → RotatableArray（ABC199 C）
		1 文字の書き換えと区間の比較・回文判定 → 09_mutable_string　MutableString::same ／ lcp ／ is_palindrome（ABC331 F）
		1 文字の書き換えと区間の文字の数・昇順か・k 番目の文字 → MutableString::count ／ counts ／ is_sorted ／ kth_char（ABC285 F）
		挿入・削除・区間の反転を何度も → 体系外（データ構造：DynamicString）
	▼ ほかの第二世代と組み合わせる
		mod で数える → math_v2 00_modint　count_xxx<mint>(…)
		パターンを含まない長さ N（10^18 まで）の文字列の数 → KmpAutomaton（複数なら AhoCorasick）の遷移を math_v2 05_matrix　Matrix<mint>::pow
		ある文字列が別の文字列に含まれる関係で最大の選び方（反鎖） → AhoCorasick::contained_patterns ／ PalindromicTree の link・parent ＋ graph_grid_v2 08_bipartite_matching ／ 08_max_flow（ABC354 G・ABC237 Ex）
		重なり（max_overlap）を重みにして bit DP → math_v2 09_bit（ABC343 G）
		周期的に一致する位置を辺にしてたどる → find_all_cyclic ＋ graph_grid_v2 06_functional_graph（ABC135 F）
		回文の条件から列を作る → graph_grid_v2 03_union_find で等しい組をまとめ、04_palindrome center_palindrome_lengths で確かめる（ABC349 G）
		文字列を頂点にした最短路 → graph_grid_v2 00_state_id ＋ 01_bfs
		数字の文字列・回文数・大きな数の余り → math_v2 01_digits
		答えが巨大な数（10^i の係数を繰り上げて出力） → math_v2 01_digits　carry_to_string（ABC379 E）
		縦書き（文字列の列を 90 度回す） → graph_grid_v2 10_grid_transform　rotate_cw（ABC366 B）
```

---

## 5. 既存テンプレートとの照合（Phase 1-6）

### 5-0. 対象

Dynalist（dynalist-2026-10-11.opml）の「10_文字列」と、「01_累積・前処理」「07_データ構造」「15_道具」に載っている文字列関係の項目が指す旧ファイル（`C:\atcoder_template\典型90\ソース`）。不具合の指摘は `Desktop\claude\atcoder_template_review_20260925.md`。

| 旧ファイル | Dynalist の項目 |
|---|---|
| StringTemplateTest.txt（[191]） | StringAlgo（string 系全部入り）：Z・SA / LCP・最長の繰り返し・異なる部分文字列の数・辞書順 K 番目（O(N) / O(log N)）・Rolling Hash・最長回文・2 つ・K 個の文字列の最長共通部分文字列（長さ・1 つ・全部） |
| RollingHashTest2.txt（[193]）/ RollingHashTest.txt（[194]） | RollingHash（concat・挿入・LCP・辞書順比較）、同じ名前の旧版（修正前） |
| RollingHashSegTree.txt（[118]） | ローリングハッシュ segtree（1 点更新＋部分文字列）。07_データ構造にも |
| PalindromeHash.txt（[195]） | 1 点更新回文判定 |
| ManacherTemplate.txt（[200]） | Manacher（各位置が中心のときの半径） |
| 468D.txt（[199]） | ExpandedString（文字の間に区切りを仮想的に入れる） |
| ALDS14C.txt（[198]） | Trie（文字列の multiset） |
| StringHelperTemplateTest.txt（[168]） | StringHelper（次・前の出現・区間にあるか・出現位置の列） |
| LazyBidirectionalStringTemplateTest.txt（[166]） | 前後から高速に文字列を構築できる構造体 |
| EDPCF.txt（[192]） | LCS 最長共通部分文字列（中身は最長共通部分列の DP） |
| A24.txt / A24_2.txt（[196][197]） | 最長増加部分列（LIS）（二分探索 / セグ木）。11_区間最適化にも |
| RunLengthTemplateTest.txt（[8]） | ランレングス圧縮（数列版 / string 版）（01・15） |
| BracketSegTree.txt（[117]） | 括弧列の区間判定（07_データ構造のセグ木） |
| WaveletTreeTemplateTest.txt（[139]）/ DynamicStringTemplate.txt（[62]） | WaveletTree（string 可能）/ 動的配列 文字列 全部入り（体系外。5-4） |

旧索引（00_典型発想インデックス）には文字列のキーワードがない。

### 5-1. 新設計にあって、既存にないもの

- KMP（prefix function）・出現の列挙と数・重ならない数・置き換え・周期・境界・回転の一致・無限に繰り返した中の出現・重なり・KMP オートマトン
- Aho–Corasick、ワイルドカード一致・ずらしごとの不一致の数（畳み込み）、接尾辞オートマトン、回文木、Lyndon 分解・runs
- 接尾辞配列の O(1) の LCP・部分文字列の比較・パターンの出現回数・重複ありの k 番目・重ならない最長の繰り返し・LCP の和・本数に制限のない複数の文字列
- 回文の個数・最長の回文接頭辞・接尾辞・足して回文にする最短・最長回文部分列・回文分割・回文の部分列の数
- 辞書順：つなげて最小・独自の文字の順・列名（a, …, z, aa）・辞書の番号・最小の回転
- 部分列：部分列にならない最短の文字列・異なる部分列の数・部分列としての出現回数・先頭から何文字入るか
- LCS のビット並列・編集距離（重み・制限つき・復元）・1 回の操作で一致するか
- LIS：各位置で終わる・始まる長さ・どれかの LIS／すべての LIS に入る位置・LIS の数・2 次元・増加列での被覆
- 括弧：複数の種類・対応・深さ・入れ子の木・つなげる順・最小の挿入・O(1) の区間判定
- Trie：rank・最大の LCP・LCP の和・接頭辞でまとめて消す・ノードの情報（DP 用）
- 集合・多重集合のハッシュ（見直し一覧 7-2 の Zobrist hash）、2 次元のハッシュ、ハッシュの k 回の連結
- O(1) の回転と反転（RotatableArray）、CharMap、文字の種類 Alphabet、CharPrefix

### 5-2. 既存にあって、新設計になかったもの → 追加した

| 追加先 | 追加した機能 | 参考にした既存 |
|---|---|---|
| 03_suffix_array | 2 つの文字列の最長共通部分文字列を全部（`longest_common_substrings_all`） | StringTemplateTest.txt（lcs_all） |
| 01_rolling_hash | 長さを指定して取得（`get_len`）・接尾辞（`suffix`）・空文字列（`empty_hash`） | RollingHashTest2.txt |
| 09_mutable_string | 書き換えがあるときの LCP・辞書順の比較 | RollingHashSegTree.txt（lcp / compare） |
| 09_reversible_deque | 文字列をまとめて前後に足す・部分文字列（`substr`） | LazyBidirectionalStringTemplateTest.txt（prepend(string) / substr） |
| 00_str_base | `CharPositions::positions(c)`（出現位置の列）・`CharPrefix::has`（区間にあるか） | StringHelperTemplateTest.txt（get_pos / exists） |
| 08_trie | `words()`（全部を辞書順に） | ALDS14C.txt（get_all） |
| 07_bracket | `BracketSegTree::get(i)` / `all_balanced()` | BracketSegTree.txt（get / all_ok） |
| 04_palindrome | 奇数長・偶数長の半径を別々に（`radius_odd` / `radius_even`） | ManacherTemplate.txt（odd / even）・468D.txt |

### 5-3. 既存の不具合・使いにくさ（新設計で解消）

| 既存 | 指摘 | 新設計での解消 |
|---|---|---|
| StringTemplateTest.txt（中-24） | K 個の文字列の LCS の区切り文字を `char('{' + i)` で作り、6 本目から負の char になって未定義動作 | 区切りを整数で持つ（`of_strings`・`longest_common_substring_k` は本数に制限なし） |
| StringTemplateTest.txt（4-3） | コンストラクタで Z・SA・LCP・kth・hash・Manacher を全部作って重い。BASE が固定。kth_substring_log が K の範囲を確かめない | 用途ごとの構造体に分け、必要なものだけ作る。基数は乱数。範囲外の k は止める |
| RollingHashSegTree.txt（中-14）/ PalindromeHash.txt（中-15） | 長い文字列で 2 本目を作ると基数と累乗の表を作り直し、1 本目と食い違う。S / op / e / SegTree のグローバルな名前が重なり 2 つを同時に使えない | 基数は 1 回だけ決め、表は基数を変えずに伸ばす。中の名前は `str_internal` に入れる |
| ALDS14C.txt（中-25） | コメントは英小文字だがコードは大文字前提で範囲外アクセス。入っていない文字列を消すとカウントが壊れる。kth が総数を超えると無限ループ。kth が 1 始まり | Alphabet で範囲を確かめる。erase は入っていなければ false で何も変えない。kth は 0 始まりで範囲外は止める |
| StringHelperTemplateTest.txt（4-3） | id() が英字以外を範囲外の添字にする。nxt / prv が (N+1)×52 の vector<vector<int>> で N = 2×10^5 だと約 85MB | Alphabet で範囲を確かめる。表は 1 次元の配列で、必要な文字の種類だけ（CharPositions は O(N)） |
| ManacherTemplate.txt（4-3） | is_palindrome(l, r) が閉区間。Manacher が 3 か所（[53][191][200]）に実装 | 半開区間にそろえ、04_palindrome の 1 つにまとめる |
| EDPCF.txt（4-3） | 名前は「最長共通部分文字列」だが中身は最長共通部分列。main に直書き | 部分列（`lcs_restore`）と部分文字列（`longest_common_substring`）を名前で区別し、関数にする |
| RunLengthTemplateTest.txt（4-3） | vector 版は個数が long long、string 版は int。開始位置と元に戻す関数がない | `RunBlock{value, len, start}`（len は long long）と `run_length_decode` |
| LazyBidirectionalStringTemplateTest.txt（3-2） | pop_front / pop_back がない | ReversibleDeque に入れる |
| RollingHashTest.txt（5-2） | 「修正前」の旧版が同じ名前で登録されている | 1 つにまとめる（Dynalist の項目を置き換える） |
| 468D.txt（4-3） | Main_base の sizei に依存して単体で動かない | 依存をなくし、Manacher の中心の番号の付け方（2N − 1 個）で同じことを表す |
| A24.txt / A24_2.txt | LIS の 2 つの版（二分探索・セグ木）があり、どちらも main の近くに置かれている。狭義・広義を切り替えられない | 1 つにまとめ、strict を引数に。各位置の長さ・数え上げも |

### 5-4. 体系の外に置くもの（今のまま残す）

- DynamicStringTemplate.txt（挿入・削除・区間の反転・ハッシュ）：データ構造。書き換えだけなら 09_mutable_string
- WaveletTreeTemplateTest.txt（文字列も可）・static_seq_toolkit.cpp（接尾辞配列・Z・回文など文字列の機能を持つ）：データ構造。名前の衝突だけ確かめる
- BinaryTrieTemplate.txt（XOR）・HashMap / HashSet：データ構造
- 部分文字列の数値化（ARC141A.txt）・Palindrome Library（回文数）：math_v2 01_digits に入っている

### 5-5. AtCoder の典型との照合

グラフ・グリッド・幾何の最終照合で足りないものが見つかった反省から（数学と同じく）、Phase 1 の時点で典型と突き合わせた。問題の題名は AtCoder の問題一覧で確かめた。代表例は Phase 3・4 と最終照合で問題文を確かめ直す。

| 典型 | 代表例 | 担当 |
|---|---|---|
| ランレングス圧縮 | ABC259 C、ABC329 C、典型90 084 | 00_str_base run_length |
| 文字の置き換えを重ねる | ABC342 C | 00_str_base CharMap |
| 反転を同一視する | ABC310 C | 00_str_base min_with_reverse |
| 巨大な文字列の k 文字目 | ABC380 D、ABC242 D | 00_str_base（△） |
| 部分文字列のハッシュ | ABC141 E、ABC284 F、典型90 047 | 01_rolling_hash |
| 集合・多重集合の一致 | ABC250 E、ABC367 F | 01_set_hash |
| Z（各位置の一致の長さ） | ABC257 G、ABC284 F | 02_pattern_match match_lengths / z_array |
| 数の列の一致（KMP） | ABC150 F | 02_pattern_match find_all |
| 無限に繰り返した中の出現 | ABC135 F | 02_pattern_match find_all_cyclic |
| 回転の一致 | ABC430 E、ABC223 B | 02_pattern_match rotation_offsets、05_lex_order min_rotation |
| 重なり・包含（bit DP） | ABC343 G | 02_pattern_match max_overlap |
| 多数のパターン | ABC268 Ex、ABC362 G | 02_aho_corasick |
| 文字列の包含関係で選ぶ（反鎖） | ABC354 G、ABC237 Ex | 02_aho_corasick contained_patterns、04_palindromic_tree |
| ずらしごとの不一致の数 | ABC196 F | 02_wildcard_match mismatch_counts |
| 接尾辞配列・LCP | ABC213 F、ABC272 F、ABC280 Ex、ABC362 G、ABC240 Ex | 03_suffix_array |
| 重ならない最長の繰り返し | ABC141 E | 03_suffix_array longest_repeated_non_overlapping |
| 回文の判定・作る | ABC237 C、ABC320 B、ABC398 F | 04_palindrome |
| 回文の条件から作る | ABC349 G | 04_palindrome（△ の型） |
| 異なる回文 | ABC237 Ex | 04_palindromic_tree |
| 書き換えのある回文判定 | ABC331 F | 09_mutable_string |
| つなげて辞書順最小 | ABC225 F | 05_lex_order concat_order |
| 独自の文字の順 | ABC219 C | 05_lex_order lex_less_by |
| 列名（a〜z, aa, …） | ABC171 C、ABC285 C | 05_lex_order bijective_to_string / bijective_from_string |
| 辞書順最小の部分列 | 典型90 006 | 06_subsequence lex_smallest_subsequence |
| 部分列として何回現れるか | 典型90 008、ABC211 C | 06_subsequence count_subsequence_occurrences |
| 異なる部分列の数 | ABC214 F、ABC230 F | 06_subsequence count_distinct_subsequences |
| 部分列にならない最短の文字列 | ARC081 E | 06_subsequence shortest_non_subsequence |
| 2 つをつなげて部分列になるか | ABC324 E | 06_subsequence prefix_match_lengths / suffix_match_lengths |
| LCS | EDPC F | 06_lcs_edit lcs_restore |
| 編集距離 | ABC185 E、ABC324 C、ABC386 C・F | 06_lcs_edit |
| LIS | ABC038 D、ABC134 E、ABC354 F、ABC360 G | 06_lis |
| 括弧列 | 典型90 002、ABC167 F、ABC223 F、ABC307 D、ABC350 F、ABC394 D、ABC428 C | 07_bracket |
| Trie | ABC268 G、ABC287 E、ABC353 E、ABC377 G、ABC403 E | 08_trie |
| 前後に足す・反転 | ABC158 D、ABC199 C、ABC237 D | 09_reversible_deque |
| 書き換えと区間の文字の数 | ABC285 F | 09_mutable_string |

### 5-6. 追加しないもの（理由）

| 機能 | 理由 |
|---|---|
| 接尾辞木 | 接尾辞配列＋LCP（と接尾辞オートマトン）で代わりになる |
| 挿入・削除がある文字列の検索（動的な接尾辞配列など） | データ構造（DynamicString）の範囲。ABC ではまれ |
| BWT・FM-index・LZ 系の圧縮 | AtCoder で使わない |
| 部分文字列の周期・Lyndon 分解を区間ごとに答えるクエリ | ABC では出ない（全体の runs・周期は入れる） |
| 二重ハッシュ | 2^61 − 1 の 1 本で衝突は十分少ない（6-1 の確認 4） |
| 正規表現 | 問題ごとの DP・オートマトン（KmpAutomaton・AhoCorasick）で書く |

---

## 6. 未決事項と次にやること

### 6-1. ユーザーへの確認（2026-10-11 回答済み）

| 確認したこと | 回答・決定 |
|---|---|
| 1. 一覧（2 章の 19 ファイル）と ★ / ☆ の分け方 | このままでよい。☆ も全部作る |
| 2. 範囲 | **提案どおり**：LIS（06_lis）・集合のハッシュ（01_set_hash）・前後に足して反転できる列（09_reversible_deque）を文字列の体系に入れる。BinaryTrie・挿入削除のある文字列・WaveletTree・HashMap は対象外（データ構造） |
| 3. 文字の種類の指定 | **実行時の `Alphabet`**：引数で渡す（既定は英小文字。`Alphabet::upper()`・`Alphabet::ints(K)` など）。範囲外の文字はメッセージを出して止める |
| 4. ハッシュ | **mod 2^61 − 1 の 1 本**（二重ハッシュはしない）。基数は実行ごとの乱数で、全部の文字列で共通 |
| 個別には聞かないもの | ACL を使わない、範囲は [l, r)（_closed 版は作らない）、string と vector<整数> の両方を受け取る、名前の付け方。設計書の案のまま進める |

### 6-2. 次にやること

1. ~~6-1 の確認~~（済み）
2. ~~Phase 2：`00_str_base.cpp`、update_bundles のコピー、`dynalist_items_phase2.txt`~~（済み。2-P2 章）。ユーザーが内容を確認し、Dynalist に `dynalist_items_phase2.txt` の項目を入れる
3. ~~Phase 3：★ の 11 ファイル（終わりに独立したレビューと、代表例の問題を解く確認）~~（済み。2-P3 章）。ユーザーが内容を確認し、Dynalist に `dynalist_items_phase3.txt` の項目と索引（`dynalist_tool_index.txt` / `.opml`）を入れる
4. ~~Phase 4：☆ の 7 ファイル（同じ）~~（済み。2-P4 章）。ユーザーが内容を確認し、Dynalist に `dynalist_items_phase4.txt` の項目を入れ、索引（`dynalist_tool_index.txt` / `.opml`。☆ の部分を埋めた版）を入れ直す
5. ~~最終照合：旧ファイル（5-0）を 1 機能ずつ、過去の解答、典型（5-5 の問題文を確かめながら）と文字列の問題の広い洗い出し、Dynalist の項目~~（済み。7 章）。ユーザーが Dynalist の項目（7-7）と索引を入れ替える
6. Dynalist への反映はユーザーの環境で行う
7. 追加・修正は今までの成果物（既存のファイル・索引・Dynalist の項目・この設計書）に入れる（幾何・数学の最終照合と同じ）

## 7. 最終照合（2026-10-12〜13）

ユーザーの指示：「最終照合をお願いします。」（6-2 の 5：旧ファイルを 1 機能ずつ、過去の解答、典型（5-5 の問題文を確かめながら）と文字列の問題の広い洗い出し、Dynalist の項目）

→ 足したもの・直したものは**すべて今までの 19 ファイル・索引・Dynalist の項目・この設計書に入れた**（新しいファイルは作っていない。数学・幾何と同じ）。

### 7-1. 旧テンプレートとの照合（1 機能ずつ）

記号：◯＝新しいテンプレートに同じかそれ以上の機能がある、追加＝この照合で足した、体系外＝データ構造など別の体系。

| 旧ファイル | 旧の機能 | 新しいテンプレート | 判定 |
|---|---|---|---|
| StringTemplateTest.txt（StringAlgo） | Z・SA / LCP / rank・最長の繰り返し・異なる部分文字列の数・辞書順 K 番目（O(N) / O(log N)）・get_hash・最長回文・2 つ / K 個の最長共通部分文字列（長さ・1 つ・全部） | z_array、SuffixArray（sa / rank / lcp_adj / longest_repeated / count_distinct_substrings / kth_substring）、RollingHash::get、Palindromes::longest、longest_common_substring / longest_common_substrings_all / longest_common_substring_k | ◯ |
| RollingHashTest2.txt / RollingHashTest.txt | get / get_len / get_from / empty_hash / concat / lcp(a, i, b, j) / less_than / print_hash | get / get_len / suffix / prefix / empty_hash / StrHash の + / lcp(i, H2, j) / compare / StrHash の << | ◯ |
| RollingHashSegTree.txt | set / get_hash / equal / lcp / compare（書き換えあり） | MutableString（set / hash / same / lcp / compare） | ◯ |
| PalindromeHash.txt | set / is_palindrome（書き換えあり） | MutableString::is_palindrome | ◯ |
| ManacherTemplate.txt | odd_radius / even_radius / odd_length / even_length / is_palindrome（閉区間） | radius_odd / radius_even / length_at_center / center_palindrome_lengths / is_palindrome（半開区間） | ◯ |
| 468D.txt（ExpandedString） | 文字の間に区切りを仮想的に入れて中心から広げる（ABC468 D の解答：食い違い 1 つまでの回文を中心ごとに数える。O(N^2)） | 中心の番号 2N − 1 個（center_palindrome_lengths）、**center_palindrome_lengths_with_changes(s, k)（追加。O(N k log N)）** | ◯ |
| ALDS14C.txt（Trie） | insert(s, ±1) / find / starts_with / kth（1 始まり）/ get_all | Trie（insert / erase / contains / count_prefix / kth（0 始まり）/ words） | ◯ |
| StringHelperTemplateTest.txt | id / next / prev / get_pos / exists、MultiStringHelper（StringHelper の vector） | CharPrefix（next / prev / has）・CharPositions（positions）、vector<CharPrefix> | ◯ |
| LazyBidirectionalStringTemplateTest.txt | prepend / append（1 文字・文字列）/ reverse / front / back / [] / substr / str / clear / size / empty / reserve | ReversibleDeque<char>（push_front / push_back / reverse / front / back / [] / substr / to_string / clear / size / empty）、**reserve（追加）** | ◯ |
| EDPCF.txt | LCS の DP と復元 | lcs_restore / lcs_length | ◯ |
| A24.txt / A24_2.txt | LIS（二分探索・セグ木で復元） | lis_length / lis_restore（strict 引数） | ◯ |
| RunLengthTemplateTest.txt | run_length（数列・string） | run_length / run_length_decode | ◯ |
| BracketSegTree.txt | set / range / ok / get / all_ok | BracketSegTree（set / balance / is_balanced / get / all_balanced） | ◯ |
| ARC141A.txt | 部分文字列の数値化 | math_v2 01_digits（substr_value） | ◯（数学） |
| WaveletTreeTemplateTest.txt / DynamicStringTemplate.txt | 静的区間クエリ / 挿入・削除のある文字列 | 体系外（5-4。名前の衝突だけ確かめた） | 体系外 |

旧ファイルの機能で、新しいテンプレートにないものはない。

### 7-2. 過去の解答・Dynalist の項目との照合

| 過去の解答（旧ファイルの問題） | 新しいテンプレートでの書き方 | 判定 |
|---|---|---|
| 468D（ABC468 D Pre-Palindrome） | center_palindrome_lengths_with_changes(S, 1) で中心ごとに数える（O(N log N)。旧は O(N^2)） | ◯（final.py） |
| ALDS14C（Trie の辞書） | Trie::insert / contains | ◯ |
| EDPCF（EDPC F LCS） | lcs_restore | ◯（phase3.py） |
| A24（LIS） | lis_length | ◯ |

Dynalist の「10_文字列」の 13 項目と、ほかの章の文字列関係の項目（01・15 のランレングス、07 の括弧列セグ木・ローリングハッシュ segtree・前後から構築・StringHelper・Wavelet・DynamicString・HashSet / HashMap、08・15 の部分文字列の数値化・回文数、11 の LIS）は、すべて新しい索引の行から辿れる（体系外のものは 5-4 のとおり）。

### 7-3. 典型の問題（5-5 の残り）と足した機能の問題

- 5-5 の代表例でまだ解いていなかった **ABC240 Ex**（接尾辞配列の rank・lcp_adj で部分文字列を辞書順に並べてフェニック木の DP。N = 2.5×10^4 で 0.4 秒）・**ABC349 G**（graph_grid_v2 03_union_find でまとめ、04_palindrome の center_palindrome_lengths で確かめる。N = 2×10^5 で 0.13 秒）・**ABC242 D**・**ABC362 G**（03_suffix_array の count。3 つ目の解き方）を、問題文を開いて確かめて解いた
- 最終照合で足した機能で 13 問を解いた（ABC242 E・ABC459 D・ABC287 D・ABC352 B・ABC456 C・ABC367 B・ABC468 D・ABC394 C・ABC328 D・ABC454 D・ABC477 C・ABC312 D・ABC473 F）
- `_work\problems\final.py`（17 通り・例 39 個がすべて合う）

### 7-4. 広い洗い出し（AtCoder_文字列問題一覧.txt の 128 問）

ABC200〜477 の B〜F の文字列の問題 128 問のうち、Phase 2〜4 で解いていた 22 問を除く **106 問**を、作ったものを見ていない 3 つのエージェントに分けて、問題文（必要なら解説）を開いて Main_base ＋ テンプレートの本体だけで解かせ、入出力例をすべて通した（`_work\problems\sweep1〜3.py`。例 292 個）。制約いっぱいの入力でも、重い問題はどれも 0.13 秒以内。

判定（テンプレートの関数がそのまま主な処理になる ◯、テンプレートを一部使うか std の数行で書ける △、足りない ×）：**◯ 40・△ 64・× 2**。

| 足りなかったもの・使いにくかったもの | 決めたこと |
|---|---|
| 隣が同じにならない並べ替え（ABC459 D。×） | 00_str_base に rearrange_no_adjacent_equal / lex_smallest_no_adjacent_equal |
| 置き換えを現れなくなるまで繰り返す（ABC454 D。×。ABC394 C・ABC328 D も。replace_all をループで回すと O(N^2)） | 02_pattern_match に rewrite_until_stable / erase_until_stable |
| 区間 [l, r) に T が現れるか（ABC477 C） | 02_pattern_match に Occurrences（count / exists / next / prev） |
| k 文字まで直してよい回文（ABC468 D。旧 468D も） | 01_rolling_hash に center_palindrome_lengths_with_changes |
| 同じ長さで s 以下の回文の数・同じ長さでの順位（ABC242 E。N = 10^6 で mod） | 04_palindrome に count_palindromes_leq<T>、05_lex_order に lex_rank_same_length<T> |
| ? を決めて正しい括弧列にする数（ABC312 D。索引が「型」を指していたのにデモがなかった） | 07_bracket に count_balanced_completions<T> |
| ? 込みの共通接頭辞・接尾辞（ABC287 D）・部分列の当てはめた位置（ABC352 B）・隣の条件で区切る（ABC456 C）・両端を除く（ABC367 B・ABC366 B） | 00_str_base に common_prefix_length(a, b, wild) / common_suffix_length(a, b, wild)、subsequence_positions、split_runs_by、strip / strip_left / strip_right |
| `from_ids(vector<ll>)`・`to_lower(char)`・`split(s, "#")`・`join(vs, '\n')`・`join(set)`・`find_all(S, 'a')`・`replace_all(s, ',', '.')` がコンパイルエラー | from_ids を整数の列なら何でも、1 文字の版の変換、split の文字列の区切り、join の 1 文字の区切りと任意のコンテナ、seq(char)（1 文字を長さ 1 の列に） |
| `count_chars` が vector<int> で `x * (x - 1) / 2` が溢れる。CharPrefix::count などが int で、Main_base の ll の変数と chmax / min で並べられない | 数は long long で返す（count_chars・CharPrefix・CharPositions・MutableString） |
| `rep(c, 26) cp.total(c)` が「値 0 は文字の種類に入っていない」で止まる | CharPrefix に番号の版（count_id / next_id / prev_id / kth_id / total_id） |
| 既定の英小文字のまま大文字・記号を渡して止まると直し方がわからない | Alphabet::id のメッセージに「upper()・bytes()・with(…) などを渡す」、Alphabet::with(extra) |
| 索引から見つけにくかった（ABC202 D・ABC215 C の multiset_perm_unrank、ABC264 D の inversion_count、ABC366 B の rotate_cw、ABC379 E の carry_to_string、ABC452 D の match_end、ABC473 F の need_open、ABC375 D・ABC455 E・ABC418 D・ABC456 D の型、ABC242 D・ABC450 E の巨大な文字列の型、is_subsequence の引数の順） | 索引に行を足した（7-7）。デモに ABC242 D・ABC450 E・ABC473 F の型を足した |
| 入れないもの：文字の集合としての Alphabet（同じ文字を許す）・run_length の長さの合計の上限つき和 | std::set・min で足りる |

### 7-5. 足したもの・直したもの（ファイル別）

| ファイル | 足したもの・直したもの |
|---|---|
| 00_str_base | seq(char)、Alphabet::with・id のメッセージ、from_ids（整数の列）、count_chars / CharPrefix / CharPositions の数を long long に、CharPrefix の番号の版、1 文字の to_lower / to_upper / swap_case / caesar_shift、common_prefix_length / common_suffix_length の wild、subsequence_positions、split（文字列の区切り）、strip / strip_left / strip_right、join（任意のコンテナ・1 文字の区切り・小数は fixed で 15 桁・string_view も）、split_runs_by、rearrange_no_adjacent_equal / lex_smallest_no_adjacent_equal、str_internal::checked_mul、デモに ABC242 D・ABC450 E の型 |
| 01_rolling_hash | center_palindrome_lengths_with_changes |
| 02_pattern_match | rewrite_until_stable / erase_until_stable、Occurrences |
| 04_palindrome | count_palindromes_leq<T> |
| 05_lex_order | lex_rank_same_length<T> |
| 07_bracket | count_balanced_completions<T>、デモに ABC473 F の型 |
| 09_mutable_string | count / counts を long long に |
| 09_reversible_deque | reserve |
| ほかの 11 本 | 変更なし（同梱している基盤の更新だけ） |

関数を足したファイルには、ランダムテスト（全探索・定義どおりの愚直解）と使い方デモ（日本語の説明と類題）を足し、ヘッダの 1 行目の説明も直した。止まることの確認を 13 通り足した（合わせて 82 通り）。

### 7-6. 独立したレビューで直したもの

作ったものを見ていない別のエージェントに、足した・変えた関数を読ませ、全探索・サニタイザで確かめさせた（lex_smallest_no_adjacent_equal は 4 万件・rewrite_until_stable は 31 万件で愚直解と一致）。間違った答えを返すものはなかった。

| 見つかったこと | 直したこと |
|---|---|
| `rewrite_until_stable("xyz", "a", "aa")` が「to が from を含む」で止まる（置き換えが起きないので終わる） | 最初の置き換えのときに確かめる |
| 伸び続ける置き換えを止める長さ 10^8 では、数の列だと 1024MB を超えてから止まる | 2×10^7 に |
| join が小数を 6 桁で出す（1234567.5 → 1.23457e+06） | 小数は fixed で 15 桁 |
| strip の既定が空白だけ（CRLF の \r が残る） | 空白・タブ・改行に |
| ヘッダの機能一覧に足した関数がない・01 と 02 の節の見出しが足した関数の上に来ていた | 直した |
| |to| ≥ |from| で置き換えの回数が多いと遅い・k 文字まで直してよい回文は k が大きいと遅い | 説明に書いた |
| 数を long long にしたので `int mx; chmax(mx, c[i])` がコンパイルエラーになる | 説明（ヘッダ・索引の使い方）に「受け取る変数も ll に」と書いた（ll にそろえるほうが Main_base の書き方に合う） |

### 7-7. 索引・Dynalist の項目の更新

- 索引（`dynalist_tool_index.txt` / `.opml`）：足した関数と、7-4 の見つけにくかった使い方の行を足した（122 行 → 146 行）
- Dynalist の項目：`dynalist_items_phase3.txt`（01_rolling_hash・02_pattern_match・04_palindrome・05_lex_order・06_subsequence・07_bracket・09_reversible_deque）、`dynalist_items_phase4.txt`（09_mutable_string）を書き直し、00_str_base の項目は `dynalist_items_final_check.txt` に書いた（`dynalist_items_phase2.txt` の 00_str_base の項目と入れ替える）。Dynalist には、この 9 項目を入れ替える
- 索引・Dynalist の項目に書いた関数名が、すべてテンプレート（または math_v2・graph_grid_v2・ACL）の中にあることを確かめた

### 7-8. 確認したこと

- 19 本すべて：g++ の C++17 / 20 / 23 と clang の C++17 / 20 / 23（-Wall -Wextra -Wshadow -pedantic と -pedantic-errors で警告 0）、AddressSanitizer ＋ UndefinedBehaviorSanitizer ＋ _GLIBCXX_DEBUG、clang の -O2 の UndefinedBehaviorSanitizer（trap）で、ランダムテスト・速度計測・デモがすべて通った
- 止まることの確認 82 通りがすべて関数名と理由を出して止まる。static_assert で止まるもの（数の列で wild を省いた wildcard_match、`set(1, "c")`、小数）
- 貼り付けの確認：Main_base ＋ ACL ＋ グラフ・グリッド 46 本 ＋ 幾何 17 本 ＋ 数学 26 本と一緒に貼ってもぶつからない（2 通りの順・g++ / clang・gnu++2b / c++17・LOCAL あり / なし。最終照合で足した関数も呼んで確かめた）。単独テンプレートと一緒に貼れる。1 本ずつ Main_base に貼って -Wshadow の警告 0。グローバルな名前は 184 個で、ほかの第二世代の約 900 個と重ならない
- 問題：Phase 2〜4 の 73 通り、最終照合の 17 通り、洗い出しの 106 問の**合わせて 196 通り（179 問）・入出力例 528 個がすべて合う**（Main_base ＋ テンプレートの本体だけで）。洗い出しの解答のうち 1 問（ABC381 E）は、数を long long にしたことで `int best; chmax(best, …)` がコンパイルエラーになったので、ll に直した
