# グラフ・グリッド テンプレート 第二世代 設計書

- 作成日：2026-10-05
- 状態：**Phase 4（☆）完了。最終照合（7 章）も完了し、照合で見つかった不足 7 ファイルを追加（計 46 ファイル）**（2026-10-07）。残りは 6 章（MSVC での確認など）
- 経緯：`Desktop\claude\graph_grid_template_review.md`（引き継ぎメモ）から開始
- 次回はこのファイルの「6. 未決事項と次にやること」から再開する

---

## 1. 決定事項

### 1-0. 最重要の方針（2026-10-06 ユーザー指示）

第二世代のテンプレートは互いに関係し、旧テンプレートとの互換性もないため、**後から追加・修正しなくて済むこと**が最も重要。Phase 2 に限らず、次を守る。

- 追加でほしくなる可能性のある機能はすべて入れる。入れるか迷った機能は入れる
- 速さの改善が要らないくらい速くする
- 正しいこと（愚直解とのランダムテスト、サニタイザ、Main_base＋ACL に貼った状態でのコンパイル確認）
- 使いやすいこと

### 1-1. 進め方・既存テンプレートの扱い（メモ3章・4章を承認）

- 既存テンプレートは「過去の実験結果」として扱う。参照するのは Phase 1 の終盤（1-6 照合）だけ
- 原本は変更・削除しない。新しく作るものは `Desktop\claude\graph_grid_v2\` に出力する
- 既存の場所：`C:\atcoder_template\典型90`
- Phase 1（体系設計）→ Phase 2（基盤）→ Phase 3（★）→ Phase 4（☆）の順に進め、区切りごとに確認を取る
- **ACL は使わない。すべて std だけで自作する**（SCC・2-SAT・フローも Phase 4 で自作）

### 1-2. 共通規約

| 項目 | 決定 |
|---|---|
| 優先順位 | 速度 → int / long long と変わらない使い勝手 → 正しさ |
| ファイル | 1テンプレート1ファイル。単体で実行できる .cpp |
| テスト | ランダムテスト関数と、全機能を使う使い方関数を同じファイルに入れ、main から呼ぶ |
| 使い方の出力 | 関数名の表示に加え、何を返す・提供するかを日本語でも書く |
| 使うライブラリ | std のみ（ACL 不使用） |
| 番号 | 頂点は 0-indexed。入力の 1-indexed は読み込み時に変換する |
| 区間 | 半開区間 [l, r)（例：部分木 = [in[v], out[v])） |
| 「なし」 | 頂点・親・辺番号は -1。到達不能は `DIST_INF`（= 1LL<<60、Main_base の LINF と同じ値）。経路がないときは空の vector |
| 重み・距離の型 | すべて `Cost`（= long long）。基盤の `using Cost = long long;` の1行で変更できる |
| 再帰 | DFS 系はすべて非再帰（Visual Studio のスタック 1MB 対策）。再帰DFSは△（型）として持つ |
| 名前 | 型は PascalCase（`Graph`, `RootedTree`）、関数・メンバは snake_case（`bfs`, `topo_sort`） |
| 名前の衝突 | Main_base の名前（`ll`, `rep`, `INF`, `LINF`, `chmin`, `GraphIO` など）、ACL の名前（`dsu`, `scc_graph`, `mf_graph` など）、よく使う変数名（N, M, H, W, G, K, Q, S, T, A, B など）と同じ名前をグローバルに置かない |
| Main_base への依存 | テンプレートの中では Main_base の別名・マクロ（`ll`, `rep`, `chmin` など）を使わない（単体で動かすため） |

### 1-3. 共通API

#### グラフ：辺を構造体にして1種類に統一

```cpp
using Cost = long long;
constexpr Cost DIST_INF = 1LL << 60;

struct Edge    { int to; int id; Cost w; };   // 隣接リストの要素。重みなしは w = 1
struct RawEdge { int from, to; Cost w; };     // 入力順の辺リストの要素

Graph G(N);
G.add_edge(u, v, w);        // 有向。辺番号を返す
G.add_undirected(u, v, w);  // 無向（両方向に同じ辺番号）
for (auto& e : G[v]) { e.to; e.w; e.id; }
G.es[id];                   // {from, to, w}
G.size();                   // 頂点数
G.reversed();               // 逆辺グラフ（有向グラフで t からの距離を求めるときなど）
G.sort_adj();               // 隣接を行き先の昇順に並べる（「小さい順に訪問」用）
```

入力は名前で分けた関数で読む（引数の順番の取り違えを防ぐ）。

```cpp
Graph G = read_undirected(N, M);     // 重みなし無向
Graph G = read_directed(N, M);       // 重みなし有向
Graph G = read_undirected_w(N, M);   // 重みあり無向
Graph G = read_directed_w(N, M);     // 重みあり有向
Graph G = read_tree(N);              // 辺 N-1 本の木（read_tree_w もあり）
Graph G = read_parents(N);           // 親の配列 p_2 .. p_N の形式の木
// どれも最後の引数 base（既定 1）で、入力の番号のずれを指定できる

// 特殊な与えられ方のときは、自分で作った隣接リストから Graph を作る（追加の入口）
vector<vector<int>> adj;                   // 重みなし
vector<vector<pair<int, long long>>> wadj; // 重みあり
Graph G = graph_from_adj(adj);             // 重みなしは w = 1
Graph G = graph_from_adj(wadj);
// 隣接リストの各要素を「u → v の有向辺1本」として取り込む。
// 無向グラフを両方向に入れた隣接リストなら、そのまま無向グラフとして使える
// （その場合、辺番号は向きごとに別になる）
```

基本は read 系の関数で読む。`graph_from_adj` は、行列や文字列で与えられるなど、read 系で読めない形のときの入口として追加する（ユーザーの要望）。

#### 探索の戻り値：距離と親をまとめて返す

```cpp
struct SearchResult {
    vector<Cost> dist;   // 距離（到達不能は DIST_INF）
    vector<int>  par;    // 直前の頂点（始点・未到達は -1）
    bool reached(int v) const;
    vector<int> path(int t) const;   // 始点 → t の頂点列（到達不能なら空）
};

auto r = dijkstra(G, s);   r.dist[t];  r.par[v];  r.path(t);
auto [dist, par] = bfs(G, s);         // 分けて受け取ることもできる
```

- BFS・0-1 BFS・ダイクストラ・グリッド探索はすべて `SearchResult` を返す。BFS とダイクストラをそのまま差し替えられる
- 構造化束縛（`auto [dist, par] = ...`）は、メンバ変数をちょうど2つにしておくことで使える。そのためメンバ変数はこの2つだけにする
- 辺番号の経路が要るとき（例：ABC252 E「Road Reduction」）は、基盤の `parent_edges(G, r)` で親の辺番号を求める（O(N + M)）

#### 探索の本体：「頂点数＋隣接を列挙する処理」を受け取る

BFS・0-1 BFS・ダイクストラ・連結成分は、本体（`*_core`）を1つだけ作る。通常のグラフ・グリッド・状態グラフは、本体を呼ぶ短い関数（ラッパー）にする。

```cpp
// 本体：adj(v, push) の中で、v の隣接を push(to) / push(to, w) で渡す
template<class F> SearchResult bfs_core(int n, const vector<int>& starts, F&& adj, int target = -1);   // 始点 1 つの int 版もある
template<class F> SearchResult dijkstra_core(int n, const vector<int>& starts, F&& adj, int target = -1); // {頂点, 初期コスト} の版もある

// ラッパー
auto r = bfs(G, s);                      // 多始点は bfs(G, {s1, s2})
auto r = grid_bfs(g, g.id(si, sj), "#X"); // 壁にする文字を指定（省略すると '#'）

// 状態グラフ（例：クーポン K 枚）は本体に直接ラムダを渡す
StateId S(N, K + 1);                     // (頂点, 状態) ↔ id
auto r = dijkstra_core(S.size(), S.id(s, 0), [&](int x, auto&& push) {   // 始点 1 つなら s、複数なら {s1, s2}
    auto [v, k] = S.get(x);
    for (auto& e : G[v]) {
        push(S.id(e.to, k), e.w);                // 普通に移動
        if (k < K) push(S.id(e.to, k + 1), 0);   // クーポンを使う
    }
});
```

#### 呼び方の基本形

- 1回で終わる計算は関数：`bfs(G, s)`、`topo_sort(G)`、`tree_diameter(G)`、`kruskal(G)`
- 前処理してから何度も聞くものは構造体：`RootedTree T(G, root)`、`LCA L(T)`、`FunctionalGraph F(nxt)`、`UnionFind uf(N)`
- 複数の値を返す関数は、名前付きのメンバを持つ構造体で返す（`pair` や `tuple` で返さない）

### 1-4. ファイルの形と共通化

- 保存形式：**.cpp 1ファイル**。貼る範囲を、既存の 9/23 版（union_find.cpp など）と同じ目印で囲む
  - 開始：`// ======================= ここからテンプレート本体（コピー範囲） =======================`
  - 終了：`// ======================= ここまでテンプレート本体 =======================`
- 依存するテンプレートは**貼る範囲に同梱**し、`#ifndef` ガードで重複を防ぐ。1回のコピーで動き、2つ貼っても重複定義にならない
- ガード名は `GG_` ＋ファイル名の大文字（例：`GG_GRAPH_BASE`、`GG_ROOTED_TREE`）
- ファイル名は「分類番号2桁＋名前」（例：`05_lca.cpp`）。分類と名前が対応し、並べると分類順になる
- 基盤を直したときは、同梱している全ファイルを更新する必要がある。一括更新の仕組みを Phase 2 で用意する
- アクセスツールは「目印の間だけをコピー」に対応する必要がある（ツールの管理はこの設計の対象外）

```cpp
// 05_lca.cpp ― LCA（HLD）、2点間距離、k個上の祖先、パス
// 依存：00_graph_base, 05_rooted_tree（同梱済み）
#include <bits/stdc++.h>
using namespace std;
// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_GRAPH_BASE
#define GG_GRAPH_BASE
  … グラフ基盤 …
#endif
#ifndef GG_ROOTED_TREE
#define GG_ROOTED_TREE
  … 根付き木化 …
#endif
#ifndef GG_LCA
#define GG_LCA
  … LCA 本体 …
#endif
// ======================= ここまでテンプレート本体 =======================
void random_test() { … }   // 愚直解との比較
void usage() { … }         // 全機能の使い方（日本語の説明つき）
int main() { random_test(); usage(); }
```

---

## 2. 分類とテンプレート一覧（確定）

記号：★＝優先して整備、☆＝余裕が出てから（Phase 4）、△＝専用テンプレートは作らず「型」として持つ

### 0. 基盤（Phase 2）

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ★ `00_graph_base` | **作成済み（Phase 2）**。`Cost`・`DIST_INF`・`EdgeDir`（`DIRECTED` / `UNDIRECTED`）・`Edge{to, id, w}`・`RawEdge{from, to, w}`・`Graph`（add_edge / add_undirected / add / add_vertex / [] / es / size / edge_count / reversed / sort_adj / in_degrees / out_degrees / adj_list / adj_list_w）、入力（read_undirected / read_directed / _w 版 / read_tree / read_tree_w / read_parents）、変換（graph_from_adj ×2 / graph_from_edges ×3 / graph_from_parents ×2）、`IntSpan`、`SearchResult`（dist / par / reached / path / nearest_source）、restore_path / parent_edges / path_edges、デバッグ表示（operator<<） | なし |
| ★ `00_grid_base` | **作成済み（Phase 2）**。`CharSet`（'#'・"#X"・string から暗黙に作れる、~ で補集合）、`Grid`（H / W / s、size / in / id / pos / row / col / (i, j) / at / []、count / find / find_all / mask、neighbor / for_each4 / for_each8（id 版・(i, j) 版、方向付きの f も可）/ for_each_offsets / トーラス版 3 つ / slide / wrap / wrap_id / padded / is_border / border_cells / path_to_dirs、static: DI / DJ / DIR_CHAR / dir_of（U・u・^ など）/ reverse_dir / turn_right / turn_left / offsets_manhattan / offsets_chebyshev / offsets_euclid_sq / offsets_knight）、read_grid(H, W) / read_grid() / read_matrix<T>、operator<< | なし |
| ★ `00_state_id` | **作成済み（Phase 2）**。`StateId`：2〜4 個の組 ↔ 連番（size / id / in / first / get / get3 / get4 / ids_of）。状態付き探索の土台 | なし |

### 1. 探索（Phase 3）

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ★ `01_bfs` | **作成済み（Phase 3）**。`bfs_core(n, s / starts, adj, target = -1)`（本体）、`bfs(G, s / starts, target = -1)`（多始点・t で打ち切り） | graph_base |
| ★ `01_zero_one_bfs` | **作成済み**。`zero_one_bfs_core(n, s / starts, adj, target)`（本体。push(u, w)、w は 0 / 1）、`zero_one_bfs(G, s / starts, target)` | graph_base |
| ★ `01_dfs` | **作成済み**。非再帰DFS。`dfs(G, s / starts)`、`dfs_all(G)`（全頂点）、`dfs_core(n, s / starts, adj)`。結果 `DfsResult`：pre / post / tour（戻りも含む列。ABC213 D）/ par / par_edge / depth / in / out / roots / reached | graph_base |
| ★ `01_components` | **作成済み**。`components_core(n, adj[, alive])`（本体。alive = false の頂点は対象外）、`components(G)`、`weak_components(G)`（有向の弱連結）、`component_edge_counts(G, c)`（ABC292 D）。結果 `Components`：count / id / group(k)（IntSpan、昇順）/ size_of / comp_size / same / sizes / max_size / groups | graph_base |
| ★ `01_bipartite` | **作成済み**。`bipartite(G)`、`xor_coloring(G)`（辺の重みが奇数なら両端の色が違う、偶数なら同じ。矛盾判定）。結果 `Bipartite`：ok / color / comp / comp_count / cnt[c]（成分ごとの色別の頂点数。ABC282 D）/ comp_ok / odd_cycle（矛盾の原因の閉路）。辺の向きは無視 | graph_base |
| ★ `01_cycle` | **作成済み**。`find_cycle_directed(G)`、`find_cycle_undirected(G)`（多重辺・自己ループ対応）、`shortest_cycle_directed(G, s)`（s を通る最短の有向閉路。ABC376 D）、`two_core(G)`（次数 1 以下を削り続けて残る頂点。なもりグラフの閉路）。結果 `Cycle`：vs / es / found() / size() | graph_base |
| △ 再帰DFS | バックトラック・全探索の型（例：ABC284 E）。N が小さいときだけ使う。型は dynalist_items_phase4.txt に載せた | － |

### 2. 最短路・状態グラフ

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ★ `02_dijkstra` | **作成済み**。`dijkstra_core(n, s / starts / {{s, 初期コスト}, ...}, adj, target)`（本体）、`dijkstra(G, 同じ)`、`dijkstra_dense(n, s / starts, cost(u, v))`（O(N²)、密なグラフ） | graph_base |
| ★ `02_shortest_path_dag` | **作成済み**。`count_shortest_paths(G, r, mod)`（BFS / 0-1 BFS / ダイクストラの結果から各頂点への本数。ABC211 D。重み 0 の閉路の先は -1）、`order_by_dist(r)`、`ShortestPathDag P(G, s, t, mod)`：best / reachable / count / on_vertex / on_edge / must_vertex / must_edge（すべての最短路が通るか）/ count_through_vertex / count_through_edge / path / ds / dt / cnt_s / cnt_t。重みは 0 以上 | dijkstra |
| ★ 状態付き探索 | 専用テンプレートは作らない。`00_state_id` ＋ `*_core` で組む。組み方は state_id の使い方関数に例を載せる | － |
| ☆ `02_bellman_ford` | **作成済み（Phase 4）**。`bellman_ford(G, s / starts)`・`bellman_ford_core(n, s, adj)`：dist（−∞ は -DIST_INF）/ par / par_edge / reached / neg_inf / path / negative_cycle。`find_negative_cycle(G)`（Cycle で返す）。牛ゲーの例つき | graph_base, cycle |
| ☆ `02_dijkstra_time` | **最終照合で追加**。`dijkstra_time(G, s / starts, cost(e, t))`（時刻 t に出て辺 e を通る時間が変わる最短路。発車時刻・混雑。ABC192 E・ABC204 E）、`dijkstra_time_core(n, s, adj(v, t, push))`、`next_time(t, k, r)`。結果 `TimeSearchResult`：SearchResult ＋ par_edge / path_edges | graph_base, dijkstra |
| ☆ `02_warshall_floyd` | **作成済み（△ からテンプレートにした）**。`WarshallFloyd W(G)` / `W(N)` + add_edge + run：W[i][j] / dist / reachable / neg_inf / negative_cycle / path / update_edge（辺を足して O(N^2) で更新。ABC375 F）/ matrix。−∞ の伝播・あふれ止めつき | graph_base |

### 3. Union-Find・最小全域木

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ★ `03_union_find` | **作成済み**。**9/23 版（union_find.cpp）を採用**。`UnionFind`：unite / find / same / size / group_count / groups / members / group_ids / roots / edge_count / unite(x, y, f)（併合時に呼ぶ関数）/ add / init。新体系に合わせて、ファイル名・`#ifndef` ガードと、**Graph から作る `union_find_from(G)`（G の全辺を unite）** を加える | graph_base |
| ★ `03_kruskal` | **作成済み**。`kruskal(G)` / `kruskal_max(G)`（最大全域木）/ `kruskal(G, uf)`・`kruskal_max(G, uf)`（uf で unite 済みの組から始める）。結果 `MstResult`：cost / edges（同じ重みは辺番号の小さい順）/ connected / components / tree(G)（選んだ辺だけのグラフ） | union_find |
| ☆ `03_weighted_union_find` | **作成済み**。9/23 版をそのまま採用 ＋ `weighted_union_find_from(G)`（辺 u → v（重み w）を「値[v] - 値[u] = w」として登録） | graph_base |
| ☆ `03_rollback_union_find` | **作成済み**。9/23 版 dynamic_connectivity_offline.cpp の RollbackUnionFind をそのまま取り出した（unite / undo / snapshot / rollback）。offline 版と両方を貼るときはこちらを貼らない | なし |
| ☆ `03_persistent_union_find` | **最終照合で追加**。`PersistentUnionFind`（部分永続）：unite（時刻が 1 進む）/ now / find(x, t) / same(x, y, t) / size(x, t) / group_count(t) / when(x, y)（初めてつながった時刻） | なし |
| ☆ `03_prim` | **作成済み**。`prim(G)` / `prim_max(G)`、`prim_dense(N, cost)` / `prim_dense_max(N, cost)`（全部の組に辺があるとき O(N^2)。ABC282 E）。結果 `PrimResult`：cost / par / par_edge / edges / connected / components | graph_base |

### 4. DAG

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ★ `04_topo_sort` | **作成済み**。`topo_sort(G)`（Kahn法）、`topo_sort_lex_min(G)`、`topo_sort_lex_max(G)`。結果 `TopoResult`：order / pos / is_dag / unique（ABC291 E）。DAG 上の DP：`dag_longest_path(G)`（DagLongest：dist / par / best / end / path）、`dag_count_paths(G, s, mod)` | graph_base |
| △ DAG上のDP | トポロジカル順に回す型。最長路・経路数は 04_topo_sort に関数として入れた。型は dynalist_items_phase4.txt に載せた | － |

### 5. 木

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ★ `05_rooted_tree`（Phase 2） | **作成済み**。`RootedTree T(G, root)` / `T(G, {根の一覧})` / `RootedTree::forest(G)`（非再帰。森・木でないグラフ（全域木）・有向グラフ（到達できる部分）にも対応）。メンバ：n / root / roots / par / par_edge / par_w / depth / wdepth / sz / in / out / order / root_of / edge_child。メソッド：size / reached / children（IntSpan）/ child_count / is_leaf / subtree（IntSpan）/ subtree_size / is_ancestor / path_to_root / leaves / bfs_order / postorder / euler_tour / heights / wheights / subtree_sum / subtree_fold / root_path_sum / root_path_fold。使い方関数に木DPの型・辺の彩色（ABC146 D）を載せた | graph_base |
| ★ `05_lca` | **作成済み**。**HLD で実装**（前処理 O(N)・メモリ O(N)・クエリ O(log N)）。`LCA L(T)` / `LCA L(G, root)`：lca(u, v) / lca({...}) / dist（辺の本数。別の木なら -1）/ wdist（重み）/ kth_ancestor / jump(u, v, k) / path / on_path(u, v, a) / is_ancestor / same_tree（森）/ subtree_range / path_segments(u, v, edge)（パスを HLD の区間に分ける。☆ HLD の土台）/ hld_pos / hld_order / head / heavy | rooted_tree |
| ★ `05_tree_diameter` | **作成済み**。`tree_diameter(G, start = 0)`：length / u / v / path / edges() / centers（1 つか 2 つ。重みありにも対応）/ ecc[x]（各頂点の離心数。ABC222 F 型）。非再帰の走査 3 回 | graph_base（rooted_tree は使わない） |
| △ 木DP | rooted_tree の order を逆順に回す型。型は dynalist_items_phase4.txt に載せた（全方位は ☆ 05_rerooting） | － |
| ☆ `05_hld` | **作成済み**。`TreeQuery Q(L, init, op, e, edge, commutative)`（1 点更新・path(u, v)・subtree(v)）、`TreeLazyQuery`（path_apply / subtree_apply も）。セグ木入り。非可換の演算でも path(u, v) は u → v の順。辺に値を載せる edge モード | lca |
| ☆ `05_rerooting` | **作成済み**。`rerooting(G, e, merge, put_edge, put_vertex)`：全頂点を根にした木DP の値を O(N)（非再帰。森も可） | graph_base |
| ☆ `05_centroid` | **作成済み**。`CentroidDecomposition C(G)`：par / level / order / roots / ancestor(v, k) / dist_to_ancestor(v, k) / each_component(c, f)。`tree_centroids(G)` | graph_base |
| ☆ `05_rerooting_edges` | **最終照合で追加**。`rerooting_edges(G, e, merge, put_edge, put_vertex)`：ans（rerooting と同じ）＋ 辺を切ったときの両側の dp（from_side / to_side） | graph_base |
| ☆ `05_virtual_tree` | **最終照合で追加**。`virtual_tree(L, vs)`：指定頂点とその LCA だけの圧縮した木（vs / par / dist / wdist / given / roots / index / graph）。O(K log N) | lca |
| ☆ `05_tree_classic` | **最終照合で追加**。木の最大（重み）独立集合・最小点被覆・最小支配集合・最大（重み）マッチング（復元つき）・独立集合の数え上げ（EDPC P） | rooted_tree |

### 6. 特殊グラフ

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ★ `06_functional_graph` | **作成済み**。ダブリングの表を作らない（サイクルは剰余、木の部分は深さごとの二分探索。前処理 O(N)・メモリ O(N)）。`read_next(N, base)` も入れた。`FunctionalGraph F(nxt)`：各頂点のサイクルまでの距離・入るサイクル・サイクルに入る頂点・サイクル内の位置・サイクルの一覧と長さ、`kth(v, K)`（K ≤ 10^18）、`reach(u, v)`（u→v に到達できるか、できるなら手数、できなければ -1）。**重みを渡したとき（`FunctionalGraph F(nxt, w)`）は、K 回移動したときの重みの合計も返す（kth_sum。例：ABC241 E）**。合計は long long に収まる前提 | graph_base（IntSpan） |
| ★ `06_perm_cycles` | **作成済み**。`PermCycles C(p)`：count / cycle(c) / cycle_len / len_of / cycle_id / cycle_pos / kth(i, K)（負も可）/ pow(K) / inverse / min_swaps / sort_swaps（手順。ABC350 C）/ parity / period(limit)（周期。超えたら -1）/ period_mod(mod)。`perm_compose(a, b)`、`read_perm(N, base)` | graph_base（IntSpan） |
| ☆ `06_doubling` | **最終照合で追加**。`Doubling D(nxt, val, op, e, max_k)` / `doubling(nxt)`：kth / fold（モノイド。非可換も可）/ walk / max_steps / first_true。行き先 -1（止まる）に対応。旧 Doubling の機能をすべて含む | なし |
| ☆ `06_euler_path` | **作成済み**。`euler_path_directed / euler_path_undirected / euler_circuit_directed / euler_circuit_undirected(G, s = -1, lex_min = false)`：ok / vs / es / closed()。辞書順最小にも対応 | graph_base |

### 7. 連結性（☆、Phase 4 で作成済み）

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ☆ `07_scc` | **作成済み**。`SCC S(G)` / `scc_core(n, adj)`：count / id（トポロジカル順）/ group / size_of / sizes / groups / same / dag(G, unique) / dag_edges(G, unique) | graph_base |
| ☆ `07_lowlink` | **作成済み**。`LowLink L(G)`：橋・関節点、二重辺連結成分（two_edge_id・bridge_tree）、二重頂点連結成分（block_of_edge・block_vertices・block_cut_tree）、頂点・辺を消したときの成分の数と大きさ（旧 BCC テンプレートの機能を引き継いだ）。多重辺・自己ループ対応 | graph_base |
| ☆ `07_two_sat` | **作成済み**。`TwoSat ts(n)`：add_clause / add_imply / add_nand / must / add_equal / add_not_equal / at_most_one / add_var / satisfiable / answer | graph_base, scc |

### 8. フロー（☆、Phase 4 で作成済み）

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ☆ `08_max_flow` | **作成済み**。`MaxFlow F(n)`：add_edge / flow(s, t, limit) / min_cut / cut_edges / get_edge / edges / change_edge（Dinic 法・非再帰。ACL と同じ速さ）。`LowerBoundFlow`（下限つき: max_flow / min_flow / feasible / flow_of） | なし |
| ☆ `08_bipartite_matching` | **作成済み**。`BipartiteMatching M(L, R)`：add_edge / solve（Hopcroft-Karp・非再帰）/ pairs / match_l / match_r / min_vertex_cover / max_independent_set / min_edge_cover | なし |
| ☆ `08_min_cost_flow` | **作成済み**。`MinCostFlow F(n)`：add_edge / flow(s, t, limit) / flow_min_cost / slope / get_edge / edges（主双対法。負のコストの辺も可） | なし |
| ☆ `08_hungarian` | **作成済み（追加）**。`hungarian(a)` / `hungarian_max(a)`：割り当て問題 O(n^3)。結果 cost / row_to_col / col_to_row | なし |

### 9. グリッド探索

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ★ `09_grid_bfs` | **作成済み**。`grid_bfs(g, s / starts, walls = '#')`、`grid_bfs8`、`grid_bfs_offsets(g, s, offs, walls)`、`enclosed_cells(g, walls)`。3 つ目に int（終点のつもり）を渡すとコンパイルエラー。使い方関数にワープ（ABC184 E）・向きを状態に持つ BFS・氷の上を滑る・トーラスの例 | grid_base, bfs |
| ★ `09_grid_01bfs` | **作成済み**。`grid_01bfs(g, s / starts, walls)`：壁に入るとコスト1。使い方関数に ABC176 D（ワープ魔法）・ABC213 E（パンチ）・曲がる回数の最小（ABC246 E 型）を本体に隣接処理を渡して解く例 | grid_base, zero_one_bfs |
| ★ `09_grid_components` | **作成済み**。`grid_components(g, walls = '#', eight = false)`（壁マスは -1）、`grid_components_same(g, eight)`（同じ文字どうし）。`grid_components(g, true)` はコンパイルエラー | grid_base, components |
| ☆ `09_grid_graph` | **最終照合で追加**。`grid_graph(g, walls, eight)` / `grid_graph_w(g, cost, walls, eight)` / `grid_graph_torus` / `grid_graph_offsets(g, offs, walls, torus)` / `grid_graph_directed(g, cost, walls, eight)`：グリッドを Graph にして、グラフの道具（最短路の本数・lowlink・閉路・最小全域木・フロー）を使う | graph_base, grid_base |
| ☆ → 使い方例（**反映済み**） | 向きを状態に持つBFS（state_id ＋ core）、外周から塗って囲まれた領域を判定（外周マスからの多始点BFS）、ワープ・トーラス（grid_base の wrap と for_each_offsets を使って隣接処理を書き換える）。専用ファイルは作らず、★テンプレートの使い方関数に例を載せる | |

### 10. グリッド変換

| ファイル | 責務・提供する機能 | 依存 |
|---|---|---|
| ★ `10_grid_transform` | **作成済み**。rotate_cw / rotate_ccw / rotate180 / rotate_k / all_rotations / transpose / flip_lr / flip_ud / shift_torus（ABC300 B）/ trim / match_rotation（ABC218 C）/ match_symmetry（裏返しも）/ find_subgrid（部分一致・ワイルドカード）/ rotate_point_cw・_ccw・rotate_point（ABC275 C）/ manhattan_to_chebyshev・chebyshev_to_manhattan（45度回転。☆ から前倒し）。`vector<string>` と `vector<vector<T>>` の両方に対応 | なし |

### 対象外（決定済み）

- 2次元集計（累積和・いもす法）
- 座標圧縮（汎用テンプレートとして別に持つ）

---

## 2-P2. Phase 2（基盤）の実装メモ（2026-10-06）

### 作ったファイル（`Desktop\claude\graph_grid_v2\`）

| ファイル | 内容 | テスト |
|---|---|---|
| `00_graph_base.cpp` | グラフ基盤 | ランダム 4,800 ケース（O2）。愚直なダイクストラで SearchResult 系も確認 |
| `00_grid_base.cpp` | グリッド基盤 | ランダム 12,000 ケース |
| `00_state_id.cpp` | StateId | ランダム 3,000 ケース（全組の番号の一致） |
| `05_rooted_tree.cpp` | 根付き木（グラフ基盤を同梱） | ランダム 18,400 ケース。木・森は再帰 DFS の愚直解と完全一致、一般のグラフ・有向グラフは全域木の性質を確認 |
| `update_bundles.ps1` / `update_bundles.bat` | 同梱部分の一括更新 | pwsh 7 で確認（更新・-Check・CRLF のファイル・BOM 保持・更新後のコンパイル） |
| `dynalist_items_phase2.txt` | Dynalist に貼る項目（説明・メモ・パス） | － |

各ファイルの確認：g++ 13 / clang 18 で -Wall -Wextra -Wshadow -pedantic 警告なし、C++17 / 20 / 23、`-D_GLIBCXX_DEBUG -fsanitize=address,undefined` で全テスト通過。
Main_base.cpp（ACL 込み）に全部の本体を 2 通りの順で貼り、gnu++2b / c++17・LOCAL あり / なしでコンパイルが通ること（名前の衝突なし、`debug(G)` などのデバッグ表示も可）を確認。
MSVC は手元にないため、MSVC で問題になりやすい書き方（bits/stdc++.h 依存、GCC 拡張、ラムダでの構造化束縛のキャプチャ（C++20）、空の vector の &v[0]）を避けた。

### 実装で決めたこと

- **Graph の内部表現は `vector<vector<Edge>>` のまま**（CSR にしない）。計測（N = M = 2×10^5）で構築 20ms 前後、探索は G[v] へのランダムアクセスが支配的で CSR にしても大差がないため、`G[v]` が普通の vector である使いやすさを優先した
- 向きの引数は `enum class EdgeDir`（`DIRECTED` / `UNDIRECTED`）。`base` と順番を取り違えるとコンパイルエラーになる
- 配列から作る変換（graph_from_edges / graph_from_parents）は `base` を必須の引数にした（読んだまま 1-indexed か、自分で作った 0-indexed かを毎回明示する）
- `SearchResult` のメンバ変数は dist と par の 2 つだけ（構造化束縛のため）。いちばん近い始点は `nearest_source()`、辺番号は `parent_edges` で後から求める
- `IntSpan`（連続した int の参照）を graph_base に置いた。`T.children(v)` などが返す。Phase 3・4 でグループの一覧などを返すときにも使う
- RootedTree の子の一覧は CSR（`IntSpan` で返す）。`euler_tour()`・`edge_child`・集計（subtree_sum など）も入れた
- RootedTree の速さ：N = 2×10^5 のランダムな番号の木で構築 40〜60ms（この計測環境）。親と深さだけを求める最小の走査が 20〜50ms で、差は書き込む配列の数による。メモリのランダムアクセスが支配的なので、これ以上の高速化は見込めない
- テストとデモの形は 9/23 版（union_find.cpp）にそろえた：`xxx_test::random_test()`（Release は 8 倍の量）、`benchmark()`、`xxx_demo::demo()`（API 一覧＋【番号】ごとの説明）、`setup_console()`
- Visual Studio を C++17 で使う場合の注意：自分のコードで `auto [v, k] = ...` の v, k をラムダでキャプチャするのは C++20 から（C++17 ではコンパイルエラーになることがある）

### 次回以降の作業のしかた

- 本体の正は各 .cpp の「ここからテンプレート本体（コピー範囲）」の中。基盤を直したら `update_bundles.bat` を実行する
- Phase 3 のファイルも、同梱する基盤は既存の .cpp から取り出して入れ、最後に `update_bundles.bat -Check` で一致を確認する

## 2-P3. Phase 3（★）の実装メモ（2026-10-06）

### 作ったファイル（`Desktop\claude\graph_grid_v2\`）

| ファイル | 同梱している本体 | テスト（Release のケース数） |
|---|---|---|
| `01_bfs.cpp` | graph_base, bfs | 12,000 |
| `01_zero_one_bfs.cpp` | graph_base, zero_one_bfs | 12,000 |
| `01_dfs.cpp` | graph_base, dfs | 12,000 |
| `01_components.cpp` | graph_base, components | 16,000 |
| `01_bipartite.cpp` | graph_base, bipartite | 20,000 |
| `01_cycle.cpp` | graph_base, cycle | 24,000 |
| `02_dijkstra.cpp` | graph_base, dijkstra | 12,000 |
| `02_shortest_path_dag.cpp` | graph_base, dijkstra, shortest_path_dag | 12,000 |
| `03_union_find.cpp` | graph_base, union_find | 9/23 版のテスト（8,832 ケース・97 万操作）＋ union_find_from 16,000 |
| `03_kruskal.cpp` | graph_base, union_find, kruskal | 9,600 |
| `04_topo_sort.cpp` | graph_base, topo_sort | 20,000（全順列と比較） |
| `05_lca.cpp` | graph_base, rooted_tree, lca | 12,000 |
| `05_tree_diameter.cpp` | graph_base, tree_diameter | 24,000 |
| `06_functional_graph.cpp` | graph_base, functional_graph | 24,000 |
| `06_perm_cycles.cpp` | graph_base, perm_cycles | 33,600（最小交換回数は全順列の BFS と比較） |
| `09_grid_bfs.cpp` | graph_base, grid_base, bfs, grid_bfs | 20,000（グリッドから作った Graph の bfs と比較） |
| `09_grid_01bfs.cpp` | graph_base, grid_base, zero_one_bfs, grid_01bfs | 20,000 |
| `09_grid_components.cpp` | graph_base, grid_base, components, grid_components | 24,000 |
| `10_grid_transform.cpp` | grid_transform のみ | 32,000 |

確認したこと（Phase 2 の 4 ファイルも含めた 23 ファイルすべて）：g++ 13（-O2 -Wall -Wextra -Wshadow -pedantic）で警告なし・全テスト通過、clang 18（-pedantic-errors -Wall -Wextra -Wshadow）で警告なし、C++20 / C++23 でコンパイル可、`-D_GLIBCXX_DEBUG -fsanitize=address,undefined` で全テスト通過、Main_base.cpp（ACL 込み）に 23 ファイルの本体を 2 通りの順で全部貼って gnu++2b / c++17・LOCAL あり / なしでコンパイル・実行（名前の衝突なし・`debug(...)` も可）、`update_bundles.ps1 -Check` で同梱部分が全ファイルで一致。

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）

| 処理 | 規模 | 時間 |
|---|---|---|
| bfs / 0-1 BFS / dijkstra | N = M = 2×10^5 | 12 / 41 / 78 ms |
| dfs_all / components / bipartite / find_cycle（有向＋無向） | N = M = 2×10^5 | 74 / 15 / 20 / 86 ms |
| ShortestPathDag（ダイクストラ 2 回＋3 つの mod で本数） | N = M = 2×10^5 | 240 ms |
| kruskal / topo_sort / topo_sort_lex_min / dag_longest_path | N = M = 2×10^5 | 32 / 9 / 23 / 14 ms |
| LCA 構築（RootedTree 込み）／ 2×10^5 クエリ | N = 2×10^5 | 54 / 68 ms |
| tree_diameter / FunctionalGraph 構築 / PermCycles 構築 | N = 2×10^5 | 41〜74 / 24 / 7 ms |
| grid_bfs / grid_bfs8 / grid_01bfs / grid_components | 2000 × 2000 | 110 / 220 / 210 / 110 ms |
| ABC176 D 型（5×5 ワープの 0-1 BFS） | 1000 × 1000 | 144 ms |

### 実装で決めたこと

- **探索の本体はすべて「始点 1 つ（int）」の版も持つ**：`bfs_core(n, s, adj)`。`{s}` と書いても動くが、clang の -Wall では「スカラーを波かっこで囲んでいる」警告が出るので、始点 1 つなら `s`、複数なら `{s1, s2}` と書く（デモもこの形）
- `dijkstra_core` / `dijkstra` は始点ごとの初期コスト `{{s1, c1}, {s2, c2}}` も受け取る
- **キューは生のポインタで回す**（bfs_core・components_core）。手書きの BFS と同じ速さになった（vector の push_back より約 2 割速い）
- **ダイクストラと辞書順のトポロジカルソートのヒープは自作の 2 分ヒープ**。std::priority_queue は `-D_GLIBCXX_DEBUG` のデバッグビルドで 1 回の操作ごとにヒープ全体を検査する（O(n)）ため、手元のデバッグ実行で N = 2×10^5 が終わらなくなる。自作のヒープは Release でも priority_queue と同じかそれより速い。`lower_bound` などの二分探索も同じ理由で本体では使っていない
- 利用者が受け取るフラグの配列は `vector<int>`（`vector<char>` はデバッグ表示で文字化けするため）。内部の作業用だけ `vector<char>`
- **LCA は HLD で実装**（ダブリングにしない）。前処理 O(N)・メモリ O(N) で、クエリも速い。`path_segments(u, v, edge)` でパスを HLD の区間に分けられるので、☆ の HLD（パス上の集計）はこれを土台にする
- **FunctionalGraph はダブリングの表を作らない**（N × 60 の表を作らず、サイクルは剰余、木の部分は深さごとの行きがけ順の二分探索）
- tree_diameter は rooted_tree を使わず、graph_base だけに依存する（非再帰の走査 3 回で、離心数 ecc も同時に求まる）
- functional_graph・perm_cycles は `IntSpan`（サイクルの頂点）のため graph_base を同梱する
- グリッドの探索は、3 つ目の引数に int を渡すと**コンパイルエラー**（`grid_bfs(g, s, t)` と終点のつもりで書く間違いを防ぐ。`= delete`）。`grid_components(g, true)` も同じ。`grid_components(g, "#X")` が bool 版に行かないよう `const char*` 版を置いた
- ShortestPathDag の「必ず通る」（must_vertex / must_edge）は、本数を 2 つの素数（998244353 と 10^9+7）で数えて「通る本数 = 全体の本数」で判定する（両方で偶然一致する確率は無視できる）
- kruskal は同じ重みなら辺番号の小さい順に選ぶ（結果が一意に決まる）
- 45 度回転（manhattan_to_chebyshev）は☆の予定だったが、数行なので 10_grid_transform に入れた
- Phase 2 の `00_state_id.cpp` は、説明の例（`dijkstra_core(S.size(), S.id(s, 0), ...)`）だけ直した。本体は変えていない

## 2-P4. Phase 4（☆）の実装メモ（2026-10-07）

### 作ったファイル（`Desktop\claude\graph_grid_v2\`）

| ファイル | 同梱している本体 | テスト（Release のケース数） |
|---|---|---|
| `02_bellman_ford.cpp` | graph_base, cycle, bellman_ford | 24,000（歩道の DP で −∞ まで確認） |
| `02_warshall_floyd.cpp` | graph_base, warshall_floyd | 24,000（update_edge は作り直したものと比較） |
| `03_weighted_union_find.cpp` | graph_base, weighted_union_find | 9/23 版のテスト（12,736 ケース・98 万操作）＋ weighted_union_find_from 16,000 |
| `03_rollback_union_find.cpp` | rollback_union_find のみ | 12,000 |
| `03_prim.cpp` | graph_base, prim | 24,000（クラスカル法と比較） |
| `05_hld.cpp` | graph_base, rooted_tree, lca, hld | 16,000（非可換の演算・遅延あり） |
| `05_rerooting.cpp` | graph_base, rerooting | 24,000（全部の根で木DP をやり直したものと比較） |
| `05_centroid.cpp` | graph_base, centroid | 24,000 |
| `06_euler_path.cpp` | graph_base, euler_path | 24,000（全部の一筆書きを列挙して比較。辞書順最小も） |
| `07_scc.cpp` | graph_base, scc | 24,000（到達可能性の行列と比較） |
| `07_lowlink.cpp` | graph_base, lowlink | 20,000（頂点・辺を 1 つずつ消して数え直したもの、単純閉路の列挙と比較） |
| `07_two_sat.cpp` | graph_base, scc, two_sat | 24,000（全部の割り当てと比較） |
| `08_max_flow.cpp` | max_flow のみ | 36,000（下限つきは全部の流し方の列挙と比較） |
| `08_bipartite_matching.cpp` | bipartite_matching のみ | 32,000 |
| `08_min_cost_flow.cpp` | min_cost_flow のみ | 26,400（全部の流し方の列挙・ベルマンフォードでの逐次最短路と比較） |
| `08_hungarian.cpp` | hungarian のみ | 24,000（全部の割り当てと比較） |

確認したこと（39 ファイルすべて）：g++ 13（-O2 -Wall -Wextra -Wshadow -pedantic）で警告なし・全テスト通過、clang 18（-pedantic-errors）で警告なし、C++20 / C++23 でコンパイル可、`-D_GLIBCXX_DEBUG -fsanitize=address,undefined` で全テスト通過、Main_base.cpp（ACL 込み）に 39 ファイルの本体を 2 通りの順で全部貼って gnu++2b / c++17・LOCAL あり / なしでコンパイル・実行（ACL の two_sat・mf_graph と同時に使えることも確認）、`update_bundles.ps1 -Check` で同梱部分が全ファイルで一致。Phase 2・3 の 23 ファイルは変更していない。

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）

| 処理 | 規模 | 時間 |
|---|---|---|
| bellman_ford / find_negative_cycle（最悪に近い形） | N = 2000, M = 4000 | 23 / 18 ms |
| WarshallFloyd / update_edge × N 回 | N = 400 | 49 / 58 ms |
| prim_dense（完全グラフ） | N = 5000（辺 1250 万本） | 131 ms |
| SCC / dag | N = M = 2×10^5 | 19 / 23 ms |
| LowLink（橋・関節点・ブロック・二重辺連結成分すべて） | N = 2×10^5, M = 3×10^5 | 96 ms |
| 2-SAT | 変数・条件 2×10^5 | 112 ms |
| 最大流（二部マッチングの形） / 二部マッチング（Hopcroft-Karp） | 左右 5×10^4、辺 2×10^5 | 121 / 29 ms |
| 最大流（グリッド 300 × 300、燃やす埋めるの形） | 9 万頂点 | 181 ms |
| 最小費用流（割り当て 300 × 300） / ハンガリアン法 1000 × 1000 | | 85 / 58 ms |
| HLD: 1 点更新 + パスの和 / パスに足す + パスの和（ランダムな木） | N = Q = 2×10^5 | 62 / 465 ms |
| 全方位木DP / 重心分解 / オイラー路 | N = 2×10^5 | 28 / 119 / 42 ms |

### 実装で決めたこと

- **ワーシャルフロイドは △（型）からテンプレートにした**。−∞ の伝播とあふれ止め、経路の復元、辺を足したときの O(N^2) の更新（辺を消していくクエリを逆順に解く）は、型として毎回書くには間違えやすいため
- 最大流は **AtCoder Library と同じ方式の Dinic 法を非再帰で**書いた（BFS は t に着いたら打ち切り、DFS は t から逆向きに 1 回で複数の道へ流す）。計測で ACL の maxflow と同じ速さ。下限つきの流れ（LowerBoundFlow）も同じファイルに入れた
- 最小費用流は主双対法。ダイクストラは t が確定したら打ち切り、確定した頂点だけポテンシャルを更新する（ACL と同じ方法。速さもほぼ同じ）。負のコストの辺があるときは最初だけベルマンフォード。`flow_min_cost`（コストが減る間だけ流す）を追加
- 二部マッチングは最大流とは別に Hopcroft-Karp で作った（最大流より 4 倍速い）。旧テンプレートの「座標圧縮 → グラフ構築 → 最大マッチング」の 3 段階をやめ、`add_edge` と `solve` だけにした。König の定理で最小点被覆・最大独立集合、最小辺被覆も返す
- **割り当て問題用にハンガリアン法（08_hungarian）を追加**（設計にはなかった。全部の組にコストがある割り当ては最小費用流より速く、書くのも楽）
- HLD は 05_lca の重軽分解（path_segments）の上にセグ木を載せた。非可換の演算に対応するため逆向きの畳み込みも持つ（commutative = true なら持たず 2〜3 割速い。ACL の遅延セグ木で書いた場合と同じ速さ）。辺に値を載せる edge モードつき
- 全方位木DP は「merge / put_edge / put_vertex」の 3 つの関数で書く汎用版。左右からの累積を使うので逆元が要らない。旧テンプレートの最大距離・距離の和は使い方デモで 3 つの関数を変えるだけで書ける
- 重心分解は、各頂点から重心の木の祖先までの距離を前もって持つ（「一番近い印までの距離」が O(log N)）。成分を列挙する each_component(c, f) は枝の番号も渡す（同じ枝どうしを引く数え上げに使う）
- lowlink は旧 BCC テンプレートの「頂点を消したときの成分の数・大きさ」も引き継いだ。二重辺連結成分・二重頂点連結成分（ブロック）・ブロックカット木も入れた（設計では「当面不要」だったが、迷ったら入れる方針で追加）
- オイラー路は辞書順最小（lex_min）にも対応（全部の一筆書きを列挙するテストで確認）
- 重み付き Union-Find・巻き戻せる Union-Find は 9/23 版をそのまま採用。巻き戻せる版は dynamic_connectivity_offline.cpp にも同じ構造体が入っているので、両方を貼るときはこちらを貼らない（ヘッダに書いた）
- 2-SAT は SCC（07_scc）を同梱して使う。at_most_one は補助の変数で O(k) 個の条件にする

## 2-F. 最終照合で追加したファイルの実装メモ（2026-10-07）

照合（7 章）で × になった機能を、既存の 39 ファイルを変えずに新しいファイル（新しいガード）として足した。

### 作ったファイル（`Desktop\claude\graph_grid_v2\`）

| ファイル | 同梱している本体 | テスト（Release のケース数） |
|---|---|---|
| `02_dijkstra_time.cpp` | graph_base, dijkstra, dijkstra_time | 12,000（「着く時刻」を何度も緩和する愚直と比較。発車時刻・ラッシュアワーの 2 種類） |
| `03_persistent_union_find.cpp` | persistent_union_find のみ | 24,000（時刻ごとに成分を全部持つ愚直と比較。when は全組） |
| `05_rerooting_edges.cpp` | graph_base, rerooting_edges | 18,400（辺を 1 本ずつ消して両側を再帰で木DP したものと比較。4 種類の DP） |
| `05_virtual_tree.cpp` | graph_base, rooted_tree, lca, virtual_tree | 16,000（全組の LCA を足した集合・親をたどった愚直と比較。森・届かない頂点も） |
| `05_tree_classic.cpp` | graph_base, rooted_tree, tree_classic | 16,000（頂点・辺の部分集合を全部調べたものと比較。重みあり・なし、選んだ集合の正しさも） |
| `06_doubling.cpp` | doubling のみ | 12,000（1 回ずつたどる愚直と比較。最小値・非可換（一次関数の合成）・行き先 -1） |
| `09_grid_graph.cpp` | graph_base, grid_base, grid_graph | 24,000（for_each4 / for_each8 / トーラス / 任意の移動の隣接と比較。BFS の距離も） |

確認したこと（46 ファイルすべて）：g++（-O2 -Wall -Wextra -Wshadow -pedantic）で警告なし・全テスト通過、clang（-pedantic-errors）で警告なし、C++20 / C++23 でコンパイル可、`-D_GLIBCXX_DEBUG -fsanitize=address,undefined` で全テスト通過、Main_base.cpp（ACL 込み）に 46 ファイルの本体を 2 通りの順で全部貼って gnu++2b / c++17・LOCAL あり / なしでコンパイル・実行、`update_bundles.ps1 -Check` で同梱部分が全ファイルで一致。**既存の 39 ファイルは 1 バイトも変えていない**（納品済みのものとバイト単位で一致を確認）。

### 速さ（この計測環境。AtCoder のジャッジはおよそ 2 倍速い）

| 処理 | 規模 | 時間 |
|---|---|---|
| dijkstra_time（発車時刻つき） / 参考: dijkstra | N = 2×10^5, M = 4×10^5 | 162 / 61 ms |
| PersistentUnionFind: unite / when + size + same | N = Q = 2×10^5 | 49 / 44 ms |
| rerooting_edges（距離の和）/ 星 | N = 2×10^5 | 80〜100 / 16 ms |
| virtual_tree（K の合計 2×10^5。K = 2 / 200 / 2×10^5） | N = 2×10^5 | 75 / 59 / 51 ms |
| tree_classic（最大独立集合 / 最小支配集合 / 最大重みマッチング / 独立集合の数） | N = 2×10^5 | 52 / 55 / 71 / 46 ms |
| Doubling（最小値・max_k = 10^18）: 作る / fold × 2×10^5 回、値なし: 作る / kth × 2×10^5 回 | N = 2×10^5 | 93 / 180〜230 / 35 / 145 ms |
| grid_graph（4 方向 / 8 方向）/ その上の BFS | 2000 × 2000 | 391 / 639 / 180 ms |

### 実装で決めたこと

- **02_dijkstra_time**：既存の `dijkstra_core` は adj に今の距離を渡さないので、時刻で重みが変わる問題（ABC192 E・ABC204 E）が書けなかった。adj(v, t, push) に今の時刻 t を渡す本体と、辺ごとに cost(e, t) を返すグラフ版を作った。ヒープは 02_dijkstra のものを共有（同梱）。「早く着いても損をしない」（FIFO）が前提で、待てる問題なら待ち時間も返せば成り立つ。ABC204 E の式（x + C + D / (x + 1) の最小）は問題に依存するので本体に入れず、使い方デモにコードを載せた
- **03_persistent_union_find**：時刻 = unite を呼んだ回数（併合しなかった回も 1 と数える。「k 本目の辺を足した直後」が時刻 k になる）。経路圧縮をしない union by size で、各頂点が親につながった時刻を持つ。when(x, y) は「つながった時刻の早いほうを親へ進める」と O(log n)
- **05_virtual_tree**：LCA の HLD 順（行きがけ順の一つ）で並べ、隣どうしの LCA を足してから、スタックで親を決める。L の根から届かない頂点は無視する（森にも対応）。元の頂点 → 添字は配列を作らず二分探索（クエリごとに O(N) の配列を作らないため）
- **05_rerooting_edges**：05_rerooting と同じ書き方で、辺を切った両側の dp（put_edge をかける前の値）を返す。木でないグラフでも落ちない（BFS の全域森で計算し、使わなかった辺は e のまま）。子への辺かどうかの判定は BFS のときに印を付けて O(1)（星でも O(N)）
- **06_doubling**：表は「段ごとに N 個」の並び。行き先 -1（止まる）に対応し、walk は止まったところまで進む。max_steps / first_true は単調な条件の二分探索（旧 Doubling の first_true の機能）。k ≥ 2^LOG は abort（黙って間違えないため）。値なし版は空の型で表を小さくした。旧 Doubling と名前が同じなので同時に貼らない（ヘッダに書いた）
- **05_tree_classic**：旧テンプレートの木DP「最大独立集合」「最大マッチング」が第二世代では「型で書く」だったので、よく出る最適化をまとめて関数にした。最小点被覆は最大重み独立集合の補集合（重み 0 の頂点で要らないものは外す）、最小支配集合は 3 状態の木DP（重みつき）、重みなしのマッチングは葉からの貪欲、重みつきは木DP。すべて選んだ頂点・辺を返す
- **09_grid_graph**：頂点 = マスの番号（壁のマスも辺のない頂点として残すので、番号がずれない）。辺は行優先で「右・下（8 方向なら右下・左下も）」の順に 1 本ずつ。トーラスで H・W が 2 以下のときは for_each4_torus と同じ数だけ辺ができる（自己ループは張らない）。`grid_graph(g, true)` は 2 つ目を壁の文字と取り違えるのでコンパイルエラーにした。速さはほぼ Graph の作成（頂点ごとの vector の確保）で決まる

## 3. 依存関係

```
00_graph_base ─┬─ 01_bfs ─────────────┬─ 09_grid_bfs ─────────┐
               ├─ 01_zero_one_bfs ────┼─ 09_grid_01bfs ───────┼─ 00_grid_base
               ├─ 01_components ──────┴─ 09_grid_components ──┘
               ├─ 01_dfs
               ├─ 01_bipartite
               ├─ 01_cycle ── 02_bellman_ford
               ├─ 02_dijkstra ─┬─ 02_shortest_path_dag
               │               └─ 02_dijkstra_time
               ├─ 02_warshall_floyd
               ├─ 03_union_find ── 03_kruskal
               ├─ 03_weighted_union_find
               ├─ 03_prim
               ├─ 04_topo_sort
               ├─ 05_rooted_tree ─┬─ 05_lca ─┬─ 05_hld
               │                  │          └─ 05_virtual_tree
               │                  └─ 05_tree_classic
               ├─ 05_tree_diameter / 05_rerooting / 05_rerooting_edges / 05_centroid
               ├─ 06_functional_graph / 06_perm_cycles / 06_euler_path
               ├─ 07_scc ── 07_two_sat
               ├─ 07_lowlink
               └─ 09_grid_graph（＋ 00_grid_base）

依存なし（単独で使う）：00_state_id, 10_grid_transform, 03_rollback_union_find, 03_persistent_union_find,
                       06_doubling, 08_max_flow, 08_bipartite_matching, 08_min_cost_flow, 08_hungarian
```

- 矢印の先（右側）のファイルは、左側を同梱する。例：`05_lca` には graph_base と rooted_tree が入る（Phase 3 で tree_diameter は graph_base だけに、functional_graph・perm_cycles は graph_base 同梱に変えた）
- `00_state_id` はどのファイルにも同梱しない（状態付き探索のときに一緒に貼る）
- 09 のグリッド探索は、`SearchResult` を使うため graph_base も同梱する

### 作る順番

- Phase 2（基盤）：00_graph_base、00_grid_base、00_state_id、05_rooted_tree
- Phase 3（★）：01_*、02_*、03_*、04_*、05_lca、05_tree_diameter、06_*、09_*、10_*
- Phase 4（☆）：7・8章と、各分類の☆（02_bellman_ford、02_warshall_floyd、03_weighted_union_find、03_rollback_union_find、03_prim、05_hld、05_rerooting、05_centroid、06_euler_path、07_*、08_*）
- 最終照合（7 章）で追加：02_dijkstra_time、03_persistent_union_find、05_rerooting_edges、05_virtual_tree、05_tree_classic、06_doubling、09_grid_graph

---

## 4. 道具選択体系（Phase 1-5、案）

アクセスツールは Dynalist の OPML を読み込む形なので、Dynalist に追加する内容として作る（ツールの管理はこの設計の対象外）。

- 本体：`graph_grid_v2\dynalist_tool_index.txt`（タブでの字下げ。Dynalist にそのまま貼り付けられる）と、同じ内容の `dynalist_tool_index.opml`（インポート用）
- 置き場所：Dynalist の「00_典型発想インデックス」の中に1項目として置く（承認済み）。既存のグラフ系の▼（最短距離・到達可能性、木を扱う、連結管理、K回操作・周期、グラフが特殊）は、新しいテンプレートがそろってから整理する
- 辿り方：「入力の形（大分類）」→「聞かれていること（小分類）」→「テンプレートのファイル名と関数」
- 大分類：グラフ（N 頂点 M 辺）／木／i → A_i・順列／グリッド／状態を持つ探索
- 記号：記号なし＝テンプレートあり、△＝型（Phase 4 で ☆ はすべて作成したので、☆ の印は索引から外した）
- 46 ファイル（★・☆・最終照合で追加した 7 本）はすべて、索引のどこかから辿れることを確認済み（Phase 3・4・最終照合で増えた関数も索引に追加した）

### Dynalist のテンプレート項目の形（Phase 2・3 で各テンプレートを作るときに追加する）

既存の項目と同じ「説明・メモ・ファイルパス」の形にする。

```
05_lca（LCA・2点間距離・k個上の祖先・パス）
    説明
        LCA L(T)：前処理 O(N log N)
        L.lca(u, v) / L.dist(u, v) / L.kth_ancestor(v, k) / L.path(u, v) / L.jump(u, v, k)
    メモ
        依存：00_graph_base, 05_rooted_tree（同梱済み）
    （ファイルパス）
```

## 5. 既存テンプレートとの照合（Phase 1-6、完了）

対象：Dynalist（dynalist-2026-9-26.opml）の 05_グラフ・06_フロー・13_グリッド・関連する道具の約40ファイル（`C:\atcoder_template\典型90\ソース`）と、過去の解答（209D・211D・226C・268D・270C・272D・309D・339D・350D・362D・383C・387D・400D・405D・70D）。不具合の詳細は `Desktop\claude\atcoder_template_review_20260925.md` を参照。

### 5-1. 新設計にあって、既存にないもの

- グリッドBFSの本体（過去の解答では毎回 dx / dy で書いている：383C・405D・339D）と、探索本体の共通化（グラフ・グリッド・状態グラフで同じ本体）
- 非再帰の根付き木化（親・深さ・サイズ・in / out をまとめて取る。レビュー 7-2 の TreeInfo に当たる）
- 連結成分（成分ごとの頂点リスト）、二部グラフ判定、非再帰DFS（行きがけ順・帰りがけ順・訪問列）
- 最短路DAG（本数・最短路に乗るかの判定）、多始点、parent_edges
- トポロジカルソートの queue 版と一意性の判定（既存は常に辞書順最小）
- LCA の jump、順列のサイクル、グリッド全体の回転・trim、StateId

### 5-2. 既存にあって、新設計になかったもの → 追加した（ユーザー承認済み）

| 追加先 | 追加した機能 | 参考にした既存 |
|---|---|---|
| 05_lca | on_path(u, v, a) | TreeTemplateTest.txt、70D.txt |
| 05_tree_diameter | 中心（重みありにも対応） | DiameterTemplate.txt（重みありは中心が空だった） |
| 06_functional_graph | K 回移動したときの重みの合計 | FunctionalGraphWeightedTemplateTest.txt |
| 00_grid_base | UDLR の文字列への変換 | 453D.txt |
| 00_grid_base | トーラスの座標の正規化 | 275C.txt（TorusGrid） |
| 00_grid_base | 任意の移動リストでの隣接、距離 D 以内の移動リストの生成 | 272D.txt（決まった距離へのジャンプ）、400D.txt（2マス先）、275C.txt（Manhattan / Chebyshev） |
| 10_grid_transform | 点の90度回転（任意の中心） | 275C.txt（rot90） |
| 00_graph_base | 隣接リストからの変換 graph_from_adj | GraphIO（303E.txt）の戻り値の形。ユーザーの要望 |

### 5-3. 既存で過剰・重複していたもの（新設計で解消）

- 木の統合版（LCA＋木DP＋全方位＋HLD を1つの struct に入れていた）→ rooted_tree・lca に分け、全方位・HLD は☆
- 木DP 9本（値ごとに別ファイル・同じ名前の struct TreeDP・再帰）→ rooted_tree 1本と「木DPの型」
- Functional graph の重みなし版と重みあり版の重複、Doubling の LOG=61 固定（N=2e5 で約50MB）→ 1本に統合し、必要な段数だけ持つ
- オイラーツアー系5ファイル → ★では rooted_tree の in / out に吸収。残り（出入り両方の列・HLD・部分木の遅延セグ木連携）は☆
- GraphIO の引数（indeg / outdeg のポインタ、directed が後ろ）→ 名前で分けた read 関数

### 5-4. 既存をそのまま採用するもの

- `union_find.cpp`（9/23 版）→ ★ 03_union_find。API はそのまま。ファイル名・ガード・`union_find_from(G)` だけ加える
- `weighted_union_find.cpp`（9/23 版）→ ☆ 重み付き Union-Find（Graph との連携は Phase 4 で合わせる）
- `dynamic_connectivity_offline.cpp` / `dynamic_connectivity_online.cpp`（9/23 版）→ 体系の外の上級テンプレートとしてそのまま残す（RollbackUnionFind は offline 版の中にある）

### 5-5. Phase 4（☆）で参考にした既存（反映済み）

| ☆ | 既存 | レビューでの指摘（要約） |
|---|---|---|
| SCC | SCCTest.txt（ACL の scc_graph を使用） | build 後に indeg が 0 になる。縮約DAG・サイズ・入次数・トポ順の API は参考になる |
| lowlink・BCC | BCCTemplateTest.txt | 成分サイズの計算が指数時間、多重辺を橋と誤判定 |
| 最大流・最小費用流 | FlowTemplateTest.txt | 最小費用流のオーバーフロー、辺番号・流量の取得・最小カットがない |
| 二部マッチング | BipartiteMatchingTemplate.txt | 3段階の呼び出しが必要、König の復元がない |
| HLD・オイラーツアー | HLDTemplate.txt、HLDWrapperTemplate.txt、EulerTour*.txt | 非可換の演算でパスの向きを考えていない、グローバル名の衝突 |
| 全方位木DP | MaxDistanceRerooting.txt、SumDistanceRerooting.txt | 固定の計算だけで、汎用（merge / add_root）がない |
| ベルマンフォード・ワーシャルフロイド | ShortestPathTemplateTest.txt | −∞ の伝播がない |

### 5-6. 体系の外に置くもの（今のまま残す。ユーザー承認済み）

- 汎用の Doubling（行き先 -1・first_true があり、グラフ以外にも使う）→ **最終照合で 06_doubling を追加した（旧 Doubling の機能をすべて含む。旧は使わなくてよい）**
- 状態遷移の列挙（161D：ルンルン数、114C：753 数。最短路ではなく列挙のひな形）
- マンハッタン距離テンプレ（ManhattanTemplateTest.txt など。45度回転の☆とは別に残す）
- Monoid DSU（420E.txt）→ 新しい UnionFind の unite(x, y, f) で代わりになる
- 木の辺彩色（146D.txt）は解答そのものなので、rooted_tree の使い方例に取り込む

## 6. 未決事項と次にやること

Phase 1（体系設計）・Phase 2（基盤）・Phase 3（★）・Phase 4（☆）・最終照合（7 章）はすべて完了（39 ファイル ＋ 照合で追加した 7 ファイル = 46 ファイル）。

1. MSVC（Visual Studio）と Windows PowerShell 5.1 での確認はユーザーの環境で行う（こちらに環境がないため未確認）
2. Dynalist に dynalist_items_phase2〜4.txt と dynalist_items_final_check.txt（最終照合で追加した 7 項目）の項目を入れ、索引（dynalist_tool_index.opml。最終照合で 13 行追加）を入れ替える
3. △（型）の 3 つ（再帰DFS・木DP・DAG上のDP）は dynalist_items_phase4.txt に型として載せた
4. 体系の外に置いたもの（5-6 章）と、テンプレートの管理方法はこの設計の対象外
5. 後から機能を足すときは、既存のファイルを直さずに新しいファイル（新しいガード）として足す（1-0 の方針）

## 7. 最終照合（2026-10-07）

対象：Dynalist（dynalist-2026-10-7.opml）のグラフ・グリッドに関係する全項目（05_グラフ・06_フロー・09_幾何・13_グリッド・14〜16 と、07_データ構造などに重複して載っているもの）が指す旧ファイル 43 本、過去の解答 18 本（209D・211D・211D_2・226C・268D_3・270C・272D・309D・339D・350D・362D・362D_2・383C・383C_2・387D・400D・405D・70D）、旧索引（00_典型発想インデックス）のグラフ関係の ▼ 7 つ。旧ファイルは 1 機能ずつ第二世代の API と突き合わせた。さらに AtCoder でよく出るグラフ・グリッドの典型とも突き合わせた。

判定：◯ ＝ 第二世代の API でそのまま（同等以上で）できる、△ ＝ 専用の API はないが数行で書ける（書き方を示す）／体系の外として残す決定済み、× ＝ できない・書くのが大変 → 新しいファイルとして追加した

### 7-1. 9/23 版（そのまま採用したもの）

コードの行単位で突き合わせた。旧のテンプレート本体の行はすべて第二世代に入っている（差は `#ifndef` ガードと `*_from(G)` の追加だけ）。

| 旧ファイル | 第二世代 | 判定 |
|---|---|---|
| union_find.cpp（133 行） | 03_union_find | ◯ 全行一致 |
| weighted_union_find.cpp（170 行） | 03_weighted_union_find | ◯ 全行一致 |
| dynamic_connectivity_offline.cpp の RollbackUnionFind | 03_rollback_union_find | ◯ 全行一致。OfflineDynamicConnectivity 本体は体系の外（5-4） |
| dynamic_connectivity_online.cpp | （体系の外） | △ そのまま残す（5-4）。第二世代と名前の衝突なし |

### 7-2. 旧テンプレート（1 機能ずつ）

| 旧ファイル | 旧の機能 | 第二世代の対応 | 判定 |
|---|---|---|---|
| 420E.txt（Monoid DSU） | 成分ごとの値をモノイドで併合・get(v)・root_value(v) | 03_union_find の `uf.unite(a, b, [&](int r, int c) { val[r] = op(val[r], val[c]); })`、`val[uf.find(v)]` | △（2 行） |
| WeightedUnionFindTemplateTest.txt（旧版） | find / weight / same / diff / size / unite / unite_check | 03_weighted_union_find の同名の関数 | ◯ |
| 〃 | consistent(x, y, w) が -1（未連結）/ 1 / 0 を返す | `!w.same(x, y)` で未連結、`w.consistent(x, y, w)` で矛盾の判定 | △（1 行） |
| SCCTest.txt | comp・sz・縮約DAG（重複辺なし）・入次数・トポロジカル順 | 07_scc：`S.id`（トポロジカル順の番号）/ `S.sizes()` / `S.dag(G)` / `S.dag(G).in_degrees()` / 番号順 = トポロジカル順 | ◯ |
| BCCTemplateTest.txt | 橋・関節点 | 07_lowlink：`L.bridges` / `L.articulations`（多重辺も正しく扱う） | ◯ |
| 〃 | 二重連結成分（辺の集合） | `L.block_of_edge[id]`（辺ごとのブロック）、`L.block_vertices(b)`。辺の一覧は `for (int id = 0; id < L.m; ++id) be[L.block_of_edge[id]].push_back(id);` | △（2 行） |
| 〃 | Block-Cut Tree（頂点 0..N-1・ブロック N..） | `L.block_cut_tree()`（同じ番号の付け方） | ◯ |
| 〃 | 頂点を消したときの各成分の大きさ・最大・成分数（全頂点分） | `L.pieces_after_removing_vertex(v)`（大きさの一覧。最大は max_element）、`L.components_after_removing_vertex(v)`（グラフ全体の成分数） | ◯ |
| FlowTemplateTest.txt | Dinic（add_edge / max_flow）、orig_cap で元の容量 | 08_max_flow：`F.add_edge` / `F.flow(s, t)` / `F.get_edge(id)`（容量と流量）。非再帰・ACL と同じ速さ | ◯ |
| 〃 | 最小費用流 min_cost_flow(s, t, f)（流せなければ -1） | 08_min_cost_flow：`auto [fl, c] = F.flow(s, t, f);` で `fl < f` なら流せない | △（1 行） |
| BipartiteMatchingTemplate.txt | 最大マッチング・組の復元・完全マッチング判定 | 08_bipartite_matching：`M.solve()` / `M.pairs()` / `M.solve() == L && L == R` | ◯ |
| 〃 | 頂点を long long の値で登録（内部で座標圧縮）・add_edge の引数の順は自由 | 番号は左右それぞれ 0 から。値は座標圧縮してから渡す（座標圧縮は汎用テンプレート。1-0 で対象外） | △ |
| ShortestPathTemplateTest.txt | dijkstra_with_prev / bfs01_with_prev / restore_path | 02_dijkstra `dijkstra(G, s)` / 01_zero_one_bfs / `r.path(t)` | ◯ |
| 〃 | bellman_ford（負閉路の有無） | 02_bellman_ford（−∞ の頂点・負閉路の取り出しまで） | ◯ |
| 〃 | warshall_floyd(N, 行列) | 02_warshall_floyd：`W(N)` に `add_edge` で入れて `run()`（行列からは 2 重ループで add_edge） | ◯ |
| 315E.txt | 有向・無向の閉路検出（頂点列） | 01_cycle：`find_cycle_directed / find_cycle_undirected`（辺番号も返す。多重辺・自己ループ対応） | ◯ |
| 〃 | トポロジカルソート（辞書順最小。DAG でなければ空） | 04_topo_sort：`topo_sort_lex_min(G)`（`is_dag` で判定） | ◯ |
| A67.txt | クラスカル（辺のリストから） | 03_kruskal：`kruskal(G)`（辺のリストは `graph_from_edges`） | ◯ |
| FunctionalGraph（重みなし・重みあり） | cycle 判定・entry・cycle の大きさ・位置・depth・k 回後・無限回後・k 回の重みの和・1 周の和 | 06_functional_graph：`on_cycle` / `entry`（無限回後に最初に着くサイクルの頂点）/ `cycle_len_of` / `cycle_pos` / `to_cycle` / `kth` / `kth_sum` / `cycle_weight` | ◯ |
| Doubling (FINAL COMPLETE VERSION) | f を k 回（行き先 -1 = 止まる）、first_true（条件を満たす最初の頂点。回数の上限つきも） | **なし → 06_doubling を追加**（行き先 -1・モノイドの集約・二分探索） | × |
| 303E.txt（GraphIO） | read_unweighted / read_weighted（有向・無向、入力の基数、入次数・出次数） | 00_graph_base：`read_undirected / read_directed / _w(N, M, base)`、`G.in_degrees()` / `G.out_degrees()`、`G.adj_list()` | ◯ |
| 146D.txt | 木の辺の彩色 | 05_rooted_tree の使い方デモ【11】 | ◯ |
| TreeTemplateTest.txt | LCA・距離・k 番目の祖先・on_path・パス・HLD の区間分割 | 05_lca：`lca / wdist / kth_ancestor / on_path / path / path_segments` | ◯ |
| 〃 | 木DP（根固定）・全方位木DP（★ を書き換える型） | 05_rooted_tree の型（order を逆順）・05_rerooting | ◯ |
| EulerTourBaseTemplate.txt | par / dep / tin / tout / rev（vertex_at）/ subtree_range / subtree_size / is_ancestor | 05_rooted_tree：`par / depth / in / out / order[pos] / subtree(v) / sz / is_ancestor` | ◯ |
| 〃 | LCA 用のオイラーツアー（2N-1）と深さの列 | `T.euler_tour()`、深さは `T.depth[x]` | ◯ |
| 〃 | first[v]（ツアーで最初に出る位置）、入りと出の列（~u） | `for (int i = 0; i < (int)et.size(); ++i) if (first[et[i]] < 0) first[et[i]] = i;`。入りと出の列を使う「根からの和 + 1 点更新」は `in / out` で区間加算・1 点取得に直せる | △ |
| EulerTourSubtreeTemplate.txt | 部分木の和・部分木への一括加算・1 点の取得・代入・加算 | 05_hld：`TreeLazyQuery` の `subtree / subtree_apply / get / set / path_apply(v, v, f)` | ◯ |
| 〃 | 位置の区間 [l, r) の集計・作用、prefix 探索 ge(x) | 位置は HLD 順（`L.hld_pos`）。二分探索まで要るときは ACL の lazy_segtree を `L.hld_pos` の順に作れば max_right が使える | △ |
| EulerTourLCATemplate.txt | Sparse Table で O(1) の LCA・dist・is_ancestor | 05_lca（HLD で O(log N)。前処理 O(N)・メモリ O(N)。速さは十分） | ◯ |
| HLDTemplate.txt / HLDWrapperTemplate.txt | lca・path_segments・query_path（和・任意の演算）・apply_path | 05_lca の `path_segments`、05_hld の `TreeQuery::path`（非可換の演算も向きどおり）/ `TreeLazyQuery::path_apply` | ◯ |
| DiameterTemplate.txt | 両端・長さ・パス・中心（重みなし・重みあり） | 05_tree_diameter（重みありでも中心を返す。離心数も） | ◯ |
| 木DP 9 本 | 部分木のサイズ / 深さ / 親 / 部分木の和 / 部分木の最大・最小 / 直径（down） | 05_rooted_tree：`sz / depth / par / subtree_sum / subtree_fold(val, max) / heights()`、05_tree_diameter | ◯ |
| 〃 | 部分木の中の距離の和 | `auto s = T.subtree_sum(T.wdepth);` → `s[v] - T.sz[v] * T.wdepth[v]` | △（2 行） |
| 〃 | 最大独立集合・最大マッチング | 木DPの型で書ける（4〜6 行）が、旧テンプレートでは専用の道具だった → **05_tree_classic を追加** | × → 追加 |
| MaxDistanceRerooting.txt / SumDistanceRerooting.txt | 各頂点から最も遠い距離・全頂点への距離の和 | 05_rerooting（デモ【1】【2】）、05_tree_diameter の `ecc` | ◯ |
| 275C.txt | rot90（任意の中心・左右） | 10_grid_transform：`rotate_point_cw / _ccw / rotate_point` | ◯ |
| 〃 | GridRange（in）・TorusGrid（norm）・Neighbors（4 / 8） | 00_grid_base：`g.in(i, j)` / `g.wrap(i, j)` / `for_each4 / for_each8`（数値のグリッドも `Grid g(H, W)` で形だけ持てる） | ◯ |
| 〃 | Manhattan / Chebyshev / Circle の範囲判定 | `Grid::offsets_manhattan(D) / offsets_chebyshev(D)`。円（di² + dj² ≤ R²）は offsets_chebyshev(R) を条件で絞る | △（1 行） |
| ManhattanTemplateTest.txt / ManhattanMSTemplateTest.txt | 最大マンハッタン距離・最遠点・ひし形の判定（点の追加・削除） | 座標の道具なので体系の外（5-6）。45 度回転は `manhattan_to_chebyshev` | △ |
| RangeSplitTemplateTest.txt / RectangleDiagonalQueryTemplate.txt | 区間・長方形の分割、斜めの直線との交点の格子点 | グラフ・グリッドの探索ではない座標計算の道具 → 体系の外（今のまま残す） | △ |
| 453D.txt | 状態付きグリッドの BFS（遷移を外から渡す）・UDLR の復元 | 00_state_id + 01_bfs の `bfs_core` + `r.path(t)` + `g.path_to_dirs(path)`。旧の restore は「最初に届いた状態」ではなく「番号最小の状態」を選んでいた（第二世代では状態ごとに dist を比べて選ぶ） | ◯ |
| 161D.txt / 114C.txt | 状態遷移の列挙（ルンルン数・753 数） | 最短路ではなく列挙のひな形 → 体系の外（5-6） | △ |

### 7-3. 過去の解答

| 解答 | していること | 第二世代での書き方 | 判定 |
|---|---|---|---|
| 209D | 木の 2 色塗り分けで偶奇 | `bipartite(G).color`（または `RootedTree` の `depth` の偶奇） | ◯ |
| 211D / 211D_2 | BFS で最短路の本数 | `count_shortest_paths(G, bfs(G, 0), MOD)` | ◯ |
| 226C | 隣接リストから BFS で到達する頂点 | `graph_from_adj`（入力は 1-indexed なので 1 引く）+ `bfs(G, N - 1).reached(v)` | ◯ |
| 268D_3 | 指定した頂点をすべて含む最小の部分木（葉を削る） | `RootedTree T(G, 指定頂点の 1 つ)`、`auto c = T.subtree_sum(mark);` で `c[v] > 0` の頂点を数える。多数のクエリなら **05_virtual_tree** | ◯ |
| 270C | 木の 2 頂点のパス | `LCA(G, 0).path(X, Y)`（または `bfs(G, X).path(Y)`） | ◯ |
| 272D | 決まった距離 √M のジャンプで BFS | `grid_bfs_offsets(g, s, Grid::offsets_euclid_sq(M))` | ◯ |
| 309D | 2 つの成分でそれぞれ最も遠い距離 | `bfs(G, 0).dist` の最大 | ◯ |
| 339D | 2 人が同時に動く BFS（N^4 状態） | `StateId S(N * N, N * N)` + `bfs_core` | ◯ |
| 350D | 成分ごとの頂点数と辺の本数 | `components(G)` + `component_edge_counts(G, c)`（または `union_find_from(G).edge_count`） | ◯ |
| 362D / 362D_2 | 頂点にも重みがあるダイクストラ・経路 | `dijkstra_core(N, {{0, A[0]}}, adj)`（`push(e.to, e.w + A[e.to])`）+ `r.path(t)` | ◯ |
| 383C / 383C_2 | グリッドの多始点 BFS | `grid_bfs(g, g.find_all('H'))` | ◯ |
| 387D | 縦・横を交互に動くグリッド BFS | `StateId S(H * W, 2)` + `bfs_core`（または始点の偶奇ごとに `grid_bfs_offsets` を 2 回） | ◯ |
| 400D | 前蹴り（壁を 2 マス壊す）の 0-1 BFS | `zero_one_bfs_core`（09_grid_01bfs のデモ【ABC213 E】と同じ形） | ◯ |
| 405D | 多始点 BFS で各マスから最寄りの出口への向き | `grid_bfs(g, g.find_all('E'))` の `par` から `g.path_to_dirs({v, par[v]})` → `"^>v<"[Grid::dir_of(c)]` | ◯ |
| 70D | LCA（距離・k 番目の祖先・パス上の判定・パス） | 05_lca | ◯ |

### 7-4. 旧索引・Dynalist の項目

- 旧索引の ▼ 最短距離（BFS・01BFS・Dijkstra・BellmanFord・WarshallFloyd）／▼ 木（Euler Tour・LCA・HLD・subtree query・Rerooting）／▼ 連結管理（UnionFind・重み付き・rollback）／▼ K回操作・周期（Doubling・Functional Graph・cycle detection）／▼ グラフが特殊（DAG DP・SCC・トポロジカルソート）は、すべて第二世代の索引（▼ グラフ・グリッド道具選択）に載っている道具に当たる。Doubling だけが第二世代になかった（→ 06_doubling）
- 第二世代の索引と Dynalist の第二世代の項目に書いた API 名（86 個・260 個）は、すべてテンプレートの中に実在することを確かめた
- Dynalist のグラフ関係の節で、上の旧ファイルの一覧から漏れている項目はなかった（Grundy 数・座標圧縮などはグラフではない）

### 7-5. AtCoder の典型との照合で見つかった不足

| 典型 | 第二世代（39 本）での状況 | 対応 |
|---|---|---|
| 時刻によってコストが変わる最短路（ABC192 E・ABC204 E。「今の時刻 t から次の発車まで待つ」） | `dijkstra_core` の adj に今の距離が渡らないので書けない | × → **02_dijkstra_time を追加** |
| グリッドをグラフにして、グラフの道具（最短路の本数・関節点・閉路・二部グラフ・最小全域木など）を使う | Graph に変える関数がない（毎回 for_each4 で辺を張る） | × → **09_grid_graph を追加** |
| 「いつ初めて u と v がつながるか」「時刻 t の時点で同じ成分か・成分の大きさ」 | 並列二分探索か、クエリをソートして Union-Find（オフライン）。オンラインでは書けない | × → **03_persistent_union_find を追加** |
| 指定した K 頂点だけの圧縮した木（虚树・auxiliary tree。ABC340 G・ABC359 G など） | LCA から作れるが 20 行ほどで間違えやすい | × → **05_virtual_tree を追加** |
| 辺を切った両側それぞれの木DP の値（全方位の「向き付きの辺」ごとの値） | rerooting は頂点ごとの答えしか返さない | × → **05_rerooting_edges を追加** |
| K 回移動したときの最小値・最大値（非可換・逆元なし）、「条件を満たす最初の回数」 | functional_graph は和だけ | × → **06_doubling を追加**（旧 Doubling の機能も含む） |
| 木の最大独立集合・最大マッチング・最小点被覆・最小支配集合・独立集合の数え上げ（EDPC P） | 木DPの型で書く | × → **05_tree_classic を追加** |
| 最短路の中で別の値を最大にする（2 番目のキー） | `ShortestPathDag` + `order_by_dist` の順に DP | △ |
| 最小全域木に辺 e を必ず入れたときの重み・e が何らかの MST に入るか | `kruskal` の木で `TreeQuery`（max）を作り、`cost + w - path(u, v)` | △ |
| 辺の本数が K 以下の最短路 | `StateId S(N, K + 1)` + `bfs_core` / `dijkstra_core`、またはベルマンフォードを K 回 | △ |
| 補グラフの BFS・K 番目の最短路・一般グラフの最大マッチング・木の同型判定・マンハッタン MST | ABC（D〜F）ではほとんど出ない。必要なら個別に書く | 追加しない |

### 7-6. 照合の結果として変えたもの

- 追加したファイル：02_dijkstra_time、03_persistent_union_find、05_rerooting_edges、05_virtual_tree、05_tree_classic、06_doubling、09_grid_graph（実装メモは 2-F 章）。既存の 39 ファイルは変えていない
- 索引（dynalist_tool_index.txt / .opml）に 13 行を足した（最短距離に dijkstra_time、到達・連結に persistent_union_find、木に rerooting_edges・virtual_tree（2 行）・tree_classic（2 行）、i → A_i に doubling（3 行）、グリッドに grid_graph（2 行）、状態を持つ探索に dijkstra_time_core）
- Dynalist の項目は `dynalist_items_final_check.txt`（7 項目）に書いた（dynalist_items_phase2〜4.txt は変えていない）
- 照合して △ にしたもの（数行で書ける・体系の外）は、書き方を 7-2・7-3 章の表に残した
