// =============================================================================
//  根付き木 RootedTree（親・深さ・部分木・行きがけ順・in/out・子の一覧・集計）   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 木（Phase 2 基盤。05_lca・05_tree_diameter・木DP・全方位木DP・HLD の土台）
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  特徴
//    - 非再帰（N = 10^6 の一直線の木でも Visual Studio のスタック 1MB で落ちない）
//    - 子は隣接リストの順に訪れる（G.sort_adj() しておけば番号の小さい子から）
//    - 森（複数の根・RootedTree::forest(G)）にも対応。木でないグラフを渡すと全域木（の 1 つ）になる
//    - 有向グラフを渡すと、根から辺の向きにたどれる頂点だけの木になる
//
//  高速化のポイント
//    - 1 回のスタック走査で親・深さ・重み付き深さ・行きがけ順を作り、部分木サイズは行きがけ順の逆順で 1 回
//    - 子の一覧は 1 本の配列（CSR）。vector<vector<int>> を作らないので速く、メモリも少ない
//
//  機能一覧                                                                       計算量
//    RootedTree T(G, root = 0)       根を決めて作る（root から届かない頂点は「未到達」）  O(N + M)
//    RootedTree T(G, {r1, r2, ...})  根を複数（森）                                     O(N + M)
//    RootedTree::forest(G)           全連結成分を、成分内で番号最小の頂点を根にして作る O(N + M)
//    T.par[v] / T.par_edge[v] / T.par_w[v]   親 / 親への辺の番号 / 親への辺の重み（根は -1, -1, 0）
//    T.depth[v] / T.wdepth[v]        根からの辺の本数 / 重みの和（未到達は -1 / 0）
//    T.sz[v]                         部分木の頂点数（T.subtree_size(v) も同じ）
//    T.order                         行きがけ順（逆順に回すと子 → 親の順。木DPはこれ）
//    T.in[v], T.out[v]               部分木 = order[in[v] .. out[v])（半開区間。BIT・セグ木に載せる）
//    T.root / T.roots / T.root_of[v] 根 / 根の一覧 / v の属する木の根
//    T.edge_child[id]                木の辺 id の子の側の頂点（木の辺でなければ -1）
//    T.children(v)                   子の一覧（for (int c : T.children(v))）。child_count(v) O(1)
//    T.subtree(v)                    部分木の頂点（行きがけ順。先頭が v）                O(1)
//    T.is_ancestor(u, v)             u が v の祖先か（u == v も true）                    O(1)
//    T.is_leaf(v) / T.reached(v) / T.size()                                               O(1)
//    T.path_to_root(v)               v, 親, ..., 根                                       O(深さ)
//    T.leaves() / T.bfs_order() / T.postorder()   葉 / 根に近い順 / 帰りがけ順             O(N)
//    T.euler_tour()                  入るときと子から戻るたびに頂点を書く列（2N - 1 個）    O(N)
//    T.heights() / T.wheights()      部分木の中で最も深い頂点までの辺の本数 / 重みの和      O(N)
//    T.subtree_sum(val)              部分木の値の和                                       O(N)
//    T.subtree_fold(val, op)         部分木の値を op（max・min など）でまとめる            O(N)
//    T.root_path_sum(val)            根から v までの値の和                                 O(N)
//    T.root_path_fold(val, op)       根から v までの値を op でまとめる                     O(N)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    rt_test::random_test()  … 愚直な実装（再帰 DFS など）と全機能を突き合わせるランダムテスト
//    rt_test::benchmark()    … 速度計測
//    rt_demo::demo()         … 全機能の使い方デモ（木DP の型・類題つき）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_GRAPH_BASE
#define GG_GRAPH_BASE
// ---- グラフ基盤（グラフ・グリッド第二世代）  Graph / 入力 / 変換 / SearchResult ----
using Cost = long long;                     // 重み・距離の型（この 1 行で変更できる）
constexpr Cost DIST_INF = Cost(1LL << 60);  // 到達不能（Main_base の LINF と同じ値）
enum class EdgeDir : unsigned char { Directed, Undirected };
constexpr EdgeDir DIRECTED = EdgeDir::Directed, UNDIRECTED = EdgeDir::Undirected;

struct Edge { int to; int id; Cost w; };   // 隣接リストの要素（重みなしの辺は w = 1）
struct RawEdge { int from, to; Cost w; };  // 辺リストの要素（G.es[id]。追加した順）

struct Graph {
    std::vector<std::vector<Edge>> adj;  // adj[v] = v から出る辺
    std::vector<RawEdge> es;             // es[id] = 辺 id

    Graph() {}
    explicit Graph(int n) : adj(n) {}
    int size() const { return (int)adj.size(); }       // 頂点数
    int edge_count() const { return (int)es.size(); }  // 辺の本数（無向辺も 1 本）
    int add_vertex() { adj.emplace_back(); return (int)adj.size() - 1; }  // 頂点を 1 つ増やし番号を返す
    int add_edge(int u, int v, Cost w = 1) {  // 有向辺 u → v。辺番号を返す
        const int id = (int)es.size();
        es.push_back({u, v, w});
        adj[u].push_back({v, id, w});
        return id;
    }
    int add_undirected(int u, int v, Cost w = 1) {  // 無向辺 u - v（両方向に同じ辺番号）
        const int id = (int)es.size();
        es.push_back({u, v, w});
        adj[u].push_back({v, id, w});
        adj[v].push_back({u, id, w});
        return id;
    }
    int add(int u, int v, Cost w, EdgeDir d) { return d == EdgeDir::Directed ? add_edge(u, v, w) : add_undirected(u, v, w); }
    std::vector<Edge>& operator[](int v) { return adj[v]; }
    const std::vector<Edge>& operator[](int v) const { return adj[v]; }

    Graph reversed() const {  // すべての辺の向きを逆にしたグラフ（辺番号はそのまま。無向辺は変わらない）
        const int n = size();
        Graph r(n);
        std::vector<int> deg(n, 0);
        for (int u = 0; u < n; ++u)
            for (const Edge& e : adj[u]) ++deg[e.to];
        for (int v = 0; v < n; ++v) r.adj[v].reserve(deg[v]);
        for (int u = 0; u < n; ++u)
            for (const Edge& e : adj[u]) r.adj[e.to].push_back({u, e.id, e.w});
        r.es.reserve(es.size());
        for (const RawEdge& e : es) r.es.push_back({e.to, e.from, e.w});
        return r;
    }
    void sort_adj() {  // 各頂点の隣接を (行き先, 重み, 辺番号) の昇順に並べる（小さい順に訪問したいとき）
        for (auto& a : adj)
            std::sort(a.begin(), a.end(), [](const Edge& x, const Edge& y) {
                return x.to != y.to ? x.to < y.to : x.w != y.w ? x.w < y.w : x.id < y.id;
            });
    }
    std::vector<int> in_degrees() const {  // 入次数（無向辺は両端に 1 ずつ、無向の自己ループは 2）
        std::vector<int> d(size(), 0);
        for (const auto& a : adj)
            for (const Edge& e : a) ++d[e.to];
        return d;
    }
    std::vector<int> out_degrees() const {  // 出次数 = G[v].size()
        std::vector<int> d(size());
        for (int v = 0; v < size(); ++v) d[v] = (int)adj[v].size();
        return d;
    }
    std::vector<std::vector<int>> adj_list() const {  // 行き先だけの隣接リスト
        std::vector<std::vector<int>> r(size());
        for (int v = 0; v < size(); ++v) {
            r[v].reserve(adj[v].size());
            for (const Edge& e : adj[v]) r[v].push_back(e.to);
        }
        return r;
    }
    std::vector<std::vector<std::pair<int, Cost>>> adj_list_w() const {  // {行き先, 重み} の隣接リスト
        std::vector<std::vector<std::pair<int, Cost>>> r(size());
        for (int v = 0; v < size(); ++v) {
            r[v].reserve(adj[v].size());
            for (const Edge& e : adj[v]) r[v].emplace_back(e.to, e.w);
        }
        return r;
    }
};

namespace gg_internal {
// 辺リストから隣接リストを作る（各頂点の容量を先に確保して再確保を避ける）
inline Graph build_graph(int n, std::vector<RawEdge>&& es, EdgeDir d) {
    Graph G(n);
    const bool und = (d == EdgeDir::Undirected);
    std::vector<int> deg(n, 0);
    for (const RawEdge& e : es) {
        ++deg[e.from];
        if (und) ++deg[e.to];
    }
    for (int v = 0; v < n; ++v) G.adj[v].reserve(deg[v]);
    for (int id = 0; id < (int)es.size(); ++id) {
        const RawEdge& e = es[id];
        G.adj[e.from].push_back({e.to, id, e.w});
        if (und) G.adj[e.to].push_back({e.from, id, e.w});
    }
    G.es = std::move(es);
    return G;
}
inline Graph read_edges(int n, int m, EdgeDir d, bool weighted, int base) {
    std::vector<RawEdge> es(m > 0 ? m : 0);
    for (RawEdge& e : es) {
        std::cin >> e.from >> e.to;
        e.from -= base, e.to -= base;
        if (weighted) std::cin >> e.w;
        else e.w = 1;
    }
    return build_graph(n, std::move(es), d);
}
}  // namespace gg_internal

// ---- 入力（base = 入力の番号の始まり。1-indexed なら 1、0-indexed なら 0） ----
inline Graph read_undirected(int n, int m, int base = 1) { return gg_internal::read_edges(n, m, UNDIRECTED, false, base); }
inline Graph read_directed(int n, int m, int base = 1) { return gg_internal::read_edges(n, m, DIRECTED, false, base); }
inline Graph read_undirected_w(int n, int m, int base = 1) { return gg_internal::read_edges(n, m, UNDIRECTED, true, base); }
inline Graph read_directed_w(int n, int m, int base = 1) { return gg_internal::read_edges(n, m, DIRECTED, true, base); }
inline Graph read_tree(int n, int base = 1) { return read_undirected(n, n - 1, base); }      // 辺 n-1 本
inline Graph read_tree_w(int n, int base = 1) { return read_undirected_w(n, n - 1, base); }  // 辺 n-1 本（重みあり）
inline Graph read_parents(int n, int base = 1) {  // p_2 .. p_n（頂点 i の親）の形の木。辺 id = i - 1（0-indexed の i の辺は i - 1）
    std::vector<RawEdge> es(n > 1 ? n - 1 : 0);
    for (int i = 1; i < n; ++i) {
        int p;
        std::cin >> p;
        es[i - 1] = {p - base, i, 1};
    }
    return gg_internal::build_graph(n, std::move(es), UNDIRECTED);
}

// ---- 変換（自分で作った配列から Graph を作る） ----
// 隣接リスト（0-indexed）から。各要素を「u → v の有向辺 1 本」として取り込む。
// 無向グラフを両方向に入れた隣接リストなら、そのまま無向グラフとして使える（辺番号は向きごとに別）
inline Graph graph_from_adj(const std::vector<std::vector<int>>& a) {
    const int n = (int)a.size();
    std::vector<RawEdge> es;
    size_t m = 0;
    for (const auto& r : a) m += r.size();
    es.reserve(m);
    for (int u = 0; u < n; ++u)
        for (int v : a[u]) es.push_back({u, v, 1});
    return gg_internal::build_graph(n, std::move(es), DIRECTED);
}
template <class W>
inline Graph graph_from_adj(const std::vector<std::vector<std::pair<int, W>>>& a) {  // {行き先, 重み}
    const int n = (int)a.size();
    std::vector<RawEdge> es;
    size_t m = 0;
    for (const auto& r : a) m += r.size();
    es.reserve(m);
    for (int u = 0; u < n; ++u)
        for (const auto& p : a[u]) es.push_back({u, p.first, (Cost)p.second});
    return gg_internal::build_graph(n, std::move(es), DIRECTED);
}
// 辺の配列から。us[i] - vs[i] が辺 i（base は必ず指定する。読んだままなら 1、0-indexed なら 0）
template <class T>
inline Graph graph_from_edges(int n, const std::vector<T>& us, const std::vector<T>& vs, EdgeDir d, int base) {
    std::vector<RawEdge> es(us.size());
    for (size_t i = 0; i < us.size(); ++i) es[i] = {(int)us[i] - base, (int)vs[i] - base, 1};
    return gg_internal::build_graph(n, std::move(es), d);
}
template <class T, class W>
inline Graph graph_from_edges(int n, const std::vector<T>& us, const std::vector<T>& vs, const std::vector<W>& ws, EdgeDir d, int base) {
    std::vector<RawEdge> es(us.size());
    for (size_t i = 0; i < us.size(); ++i) es[i] = {(int)us[i] - base, (int)vs[i] - base, (Cost)ws[i]};
    return gg_internal::build_graph(n, std::move(es), d);
}
template <class T>
inline Graph graph_from_edges(int n, const std::vector<std::pair<T, T>>& uv, EdgeDir d, int base) {
    std::vector<RawEdge> es(uv.size());
    for (size_t i = 0; i < uv.size(); ++i) es[i] = {(int)uv[i].first - base, (int)uv[i].second - base, 1};
    return gg_internal::build_graph(n, std::move(es), d);
}
// 親の配列から木（森）を作る。par[i] - base が i の親。負になる値（根の 0 や -1）は「親なし」。
// 辺は i の小さい順に番号が付く（親なしの頂点は飛ばす）
template <class T>
inline Graph graph_from_parents(const std::vector<T>& par, int base) {
    const int n = (int)par.size();
    std::vector<RawEdge> es;
    es.reserve(n);
    for (int i = 0; i < n; ++i) {
        const long long p = (long long)par[i] - base;
        if (p >= 0) es.push_back({(int)p, i, 1});
    }
    return gg_internal::build_graph(n, std::move(es), UNDIRECTED);
}
template <class T, class W>
inline Graph graph_from_parents(const std::vector<T>& par, const std::vector<W>& w, int base) {  // w[i] = i と親の辺の重み
    const int n = (int)par.size();
    std::vector<RawEdge> es;
    es.reserve(n);
    for (int i = 0; i < n; ++i) {
        const long long p = (long long)par[i] - base;
        if (p >= 0) es.push_back({(int)p, i, (Cost)w[i]});
    }
    return gg_internal::build_graph(n, std::move(es), UNDIRECTED);
}

// ---- 連続した int の並びの参照（コピーしない）。T.children(v) などが返す。for (int x : sp) / sp.size() / sp[i] ----
struct IntSpan {
    const int* b = nullptr;
    const int* e = nullptr;
    const int* begin() const { return b; }
    const int* end() const { return e; }
    int size() const { return (int)(e - b); }
    bool empty() const { return b == e; }
    int operator[](int i) const { return b[i]; }
    int front() const { return *b; }
    int back() const { return *(e - 1); }
    std::vector<int> to_vector() const { return std::vector<int>(b, e); }
};

// ---- 探索の結果（BFS・0-1 BFS・ダイクストラ・グリッド探索が共通で返す） ----
struct SearchResult {
    std::vector<Cost> dist;  // 距離（到達不能は DIST_INF）
    std::vector<int> par;    // 直前の頂点（始点・未到達は -1）

    bool reached(int v) const { return dist[v] != DIST_INF; }
    std::vector<int> path(int t) const {  // 始点 → t の頂点列（到達不能なら空）
        std::vector<int> p;
        if (!reached(t)) return p;
        for (int v = t; v != -1; v = par[v]) p.push_back(v);
        std::reverse(p.begin(), p.end());
        return p;
    }
    std::vector<int> nearest_source() const {  // 各頂点に最初に届いた始点（多始点で「いちばん近い始点」）。未到達は -1
        const int n = (int)par.size();
        std::vector<int> src(n, -2), st;
        for (int v = 0; v < n; ++v) {
            if (src[v] != -2) continue;
            if (!reached(v)) { src[v] = -1; continue; }
            int u = v;
            while (src[u] == -2 && par[u] != -1) st.push_back(u), u = par[u];
            const int s = (src[u] == -2 ? u : src[u]);
            src[u] = s;
            while (!st.empty()) src[st.back()] = s, st.pop_back();
        }
        return src;
    }
};

// 親の配列から s → t の頂点列（par[始点] = -1。t に到達していることが前提）
inline std::vector<int> restore_path(const std::vector<int>& par, int t) {
    std::vector<int> p;
    for (int v = t; v != -1; v = par[v]) p.push_back(v);
    std::reverse(p.begin(), p.end());
    return p;
}
// 各頂点に入るときに使った辺の番号（始点・未到達は -1）。r は G 上の BFS / 0-1 BFS / ダイクストラの結果。O(N + M)
inline std::vector<int> parent_edges(const Graph& G, const SearchResult& r) {
    const int n = G.size();
    std::vector<int> pe(n, -1);
    std::vector<Cost> best(n, 0);
    for (int u = 0; u < n; ++u) {
        if (!r.reached(u)) continue;
        for (const Edge& e : G[u])
            if (r.par[e.to] == u && (pe[e.to] == -1 || e.w < best[e.to])) pe[e.to] = e.id, best[e.to] = e.w;
    }
    return pe;
}
// 始点 → t の経路の辺番号の列（到達不能なら空）
inline std::vector<int> path_edges(const Graph& G, const SearchResult& r, int t) {
    std::vector<int> res;
    if (!r.reached(t)) return res;
    const std::vector<int> pe = parent_edges(G, r);
    for (int v = t; r.par[v] != -1; v = r.par[v]) res.push_back(pe[v]);
    std::reverse(res.begin(), res.end());
    return res;
}

// ---- デバッグ表示（debug(G) / debug(G[v]) / cout << e で使える） ----
inline std::ostream& operator<<(std::ostream& os, const Edge& e) { return os << "(to=" << e.to << " w=" << e.w << " id=" << e.id << ")"; }
inline std::ostream& operator<<(std::ostream& os, const RawEdge& e) { return os << "(" << e.from << "->" << e.to << " w=" << e.w << ")"; }
inline std::ostream& operator<<(std::ostream& os, const Graph& G) {
    for (int v = 0; v < G.size(); ++v) {
        os << "\n  " << v << ":";
        for (const Edge& e : G[v]) os << " " << e.to << "(w=" << e.w << ",id=" << e.id << ")";
    }
    return os;
}
#endif  // GG_GRAPH_BASE
#ifndef GG_ROOTED_TREE
#define GG_ROOTED_TREE
// ---- 根付き木（グラフ・グリッド第二世代）  非再帰。森・木でないグラフ（全域木になる）にも対応 ----
struct RootedTree {
    int n = 0;                    // 頂点数
    int root = -1;                // 根（複数あるときは roots[0]）
    std::vector<int> roots;       // 根の一覧（森のときは木ごとに 1 つ）
    std::vector<int> par;         // 親（根・未到達は -1）
    std::vector<int> par_edge;    // 親への辺の番号（根・未到達は -1）
    std::vector<Cost> par_w;      // 親への辺の重み（根・未到達は 0）
    std::vector<int> depth;       // 根からの辺の本数（未到達は -1）
    std::vector<Cost> wdepth;     // 根からの重みの和（未到達は 0）
    std::vector<int> sz;          // 部分木の頂点数（未到達は 0）
    std::vector<int> in, out;     // 行きがけ順の番号。v の部分木 = order[in[v] .. out[v])（未到達は -1）
    std::vector<int> order;       // 行きがけ順（DFS 順）。order[in[v]] = v。逆順に回すと子が親より先に来る
    std::vector<int> root_of;     // 属する木の根（未到達は -1）
    std::vector<int> edge_child;  // edge_child[辺 id] = 木の辺なら子の側の頂点、木の辺でなければ -1
    std::vector<int> child_start_, child_list_;  // 子の一覧（CSR）。children(v) で使う

    RootedTree() {}
    explicit RootedTree(const Graph& G, int r = 0) { build(G, std::vector<int>{r}, false); }
    RootedTree(const Graph& G, const std::vector<int>& rs) { build(G, rs, false); }  // 根を複数（森）
    static RootedTree forest(const Graph& G) {  // すべての連結成分を、成分内で番号最小の頂点を根にして作る
        RootedTree T;
        T.build(G, {}, true);
        return T;
    }

    int size() const { return n; }
    bool reached(int v) const { return depth[v] >= 0; }
    IntSpan children(int v) const { return {child_list_.data() + child_start_[v], child_list_.data() + child_start_[v + 1]}; }
    int child_count(int v) const { return child_start_[v + 1] - child_start_[v]; }
    bool is_leaf(int v) const { return depth[v] >= 0 && child_start_[v] == child_start_[v + 1]; }
    IntSpan subtree(int v) const {  // v の部分木の頂点（行きがけ順。先頭が v）
        if (in[v] < 0) return {};
        return {order.data() + in[v], order.data() + out[v]};
    }
    int subtree_size(int v) const { return sz[v]; }
    bool is_ancestor(int u, int v) const { return in[u] <= in[v] && in[v] < out[u]; }  // u が v の祖先か（u == v も true）
    std::vector<int> path_to_root(int v) const {  // v, 親, ..., 根
        std::vector<int> p;
        if (depth[v] < 0) return p;
        p.reserve(depth[v] + 1);
        for (; v != -1; v = par[v]) p.push_back(v);
        return p;
    }
    std::vector<int> leaves() const {  // 葉（子のない頂点。行きがけ順）
        std::vector<int> r;
        for (int v : order)
            if (child_start_[v] == child_start_[v + 1]) r.push_back(v);
        return r;
    }
    std::vector<int> bfs_order() const {  // 根に近い順（深さの昇順。同じ深さは行きがけ順の親の順）
        std::vector<int> q(roots);
        q.reserve(order.size());
        for (size_t h = 0; h < q.size(); ++h)
            for (int c : children(q[h])) q.push_back(c);
        return q;
    }
    std::vector<int> postorder() const {  // 帰りがけ順（子 → 親。兄弟は子の順）
        std::vector<int> res, st, it;
        res.reserve(order.size());
        for (int r : roots) {
            st.push_back(r), it.push_back(child_start_[r]);
            while (!st.empty()) {
                const int v = st.back();
                if (it.back() < child_start_[v + 1]) {
                    const int c = child_list_[it.back()++];
                    st.push_back(c), it.push_back(child_start_[c]);
                } else {
                    res.push_back(v), st.pop_back(), it.pop_back();
                }
            }
        }
        return res;
    }
    std::vector<int> euler_tour() const {  // 入るときと子から戻るたびに頂点を書く列（木ごとに長さ 2 * 頂点数 - 1。ABC213 D）
        std::vector<int> res, st, it;
        res.reserve(2 * order.size());
        for (int r : roots) {
            st.push_back(r), it.push_back(child_start_[r]), res.push_back(r);
            while (!st.empty()) {
                const int v = st.back();
                if (it.back() < child_start_[v + 1]) {
                    const int c = child_list_[it.back()++];
                    st.push_back(c), it.push_back(child_start_[c]), res.push_back(c);
                } else {
                    st.pop_back(), it.pop_back();
                    if (!st.empty()) res.push_back(st.back());
                }
            }
        }
        return res;
    }
    std::vector<int> heights() const {  // 高さ（部分木の中で v から最も深い頂点までの辺の本数。葉は 0、未到達は -1）
        std::vector<int> h(n, -1);
        for (int k = (int)order.size() - 1; k >= 0; --k) {
            const int v = order[k];
            if (h[v] < 0) h[v] = 0;
            if (par[v] != -1 && h[par[v]] < h[v] + 1) h[par[v]] = h[v] + 1;
        }
        return h;
    }
    std::vector<Cost> wheights() const {  // 重み付きの高さ（部分木の中で v から最も遠い頂点までの重みの和。未到達は 0）
        std::vector<Cost> h(n, 0);
        for (int k = (int)order.size() - 1; k >= 0; --k) {
            const int v = order[k];
            if (par[v] != -1 && h[par[v]] < h[v] + par_w[v]) h[par[v]] = h[v] + par_w[v];
        }
        return h;
    }
    // ---- 集計（木DPのよくある形）。未到達の頂点は val のまま ----
    template <class T>
    std::vector<T> subtree_sum(const std::vector<T>& val) const {  // 部分木の値の和
        std::vector<T> r(val);
        for (int k = (int)order.size() - 1; k > 0; --k) {
            const int v = order[k];
            if (par[v] != -1) r[par[v]] += r[v];
        }
        return r;
    }
    template <class T, class Op>
    std::vector<T> subtree_fold(const std::vector<T>& val, Op op) const {  // 部分木の値を op でまとめる（max・min・gcd など）
        std::vector<T> r(val);
        for (int k = (int)order.size() - 1; k > 0; --k) {
            const int v = order[k];
            if (par[v] != -1) r[par[v]] = op(r[par[v]], r[v]);
        }
        return r;
    }
    template <class T>
    std::vector<T> root_path_sum(const std::vector<T>& val) const {  // 根から v までの値の和（部分木への加算の累積にも使う。ABC138 D）
        std::vector<T> r(val);
        for (int v : order)
            if (par[v] != -1) r[v] = r[par[v]] + r[v];
        return r;
    }
    template <class T, class Op>
    std::vector<T> root_path_fold(const std::vector<T>& val, Op op) const {  // 根から v までの値を op でまとめる
        std::vector<T> r(val);
        for (int v : order)
            if (par[v] != -1) r[v] = op(r[par[v]], r[v]);
        return r;
    }

    void build(const Graph& G, const std::vector<int>& rs, bool all) {
        n = G.size();
        roots.clear(), order.clear();
        par.assign(n, -1), par_edge.assign(n, -1), par_w.assign(n, 0), depth.assign(n, -1), wdepth.assign(n, 0);
        sz.assign(n, 0), in.assign(n, -1), out.assign(n, -1), root_of.assign(n, -1), edge_child.assign(G.edge_count(), -1);
        order.reserve(n);
        std::vector<int> st;
        st.reserve(n);
        auto run = [&](int r) {
            if (r < 0 || r >= n || depth[r] != -1) return;
            roots.push_back(r), depth[r] = 0, root_of[r] = r, st.push_back(r);
            while (!st.empty()) {
                const int v = st.back();
                st.pop_back();
                in[v] = (int)order.size(), order.push_back(v);
                const std::vector<Edge>& a = G[v];
                for (int k = (int)a.size() - 1; k >= 0; --k) {  // 逆順に積むと、隣接の順に訪れる
                    const Edge& e = a[k];
                    const int u = e.to;
                    if (depth[u] != -1) continue;
                    depth[u] = depth[v] + 1, wdepth[u] = wdepth[v] + e.w;
                    par[u] = v, par_edge[u] = e.id, par_w[u] = e.w, root_of[u] = root_of[v];
                    st.push_back(u);
                }
            }
        };
        if (all) {
            for (int v = 0; v < n; ++v) run(v);
        } else {
            for (int r : rs) run(r);
        }
        for (int k = (int)order.size() - 1; k >= 0; --k) {
            const int v = order[k];
            sz[v] += 1;
            if (par[v] != -1) sz[par[v]] += sz[v];
        }
        child_start_.assign(n + 1, 0);
        for (int v : order) {
            out[v] = in[v] + sz[v];
            if (par[v] != -1) ++child_start_[par[v] + 1], edge_child[par_edge[v]] = v;
        }
        for (int v = 0; v < n; ++v) child_start_[v + 1] += child_start_[v];
        child_list_.assign(child_start_[n], 0);
        std::vector<int> pos(child_start_.begin(), child_start_.end() - 1);
        for (int v : order)
            if (par[v] != -1) child_list_[pos[par[v]]++] = v;
        root = roots.empty() ? -1 : roots[0];
    }
};
#endif  // GG_ROOTED_TREE
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace rt_test {

#define RT_CHECK(...)                                                                              \
    do {                                                                                           \
        if (!(__VA_ARGS__)) {                                                                      \
            std::cerr << "[rooted_tree] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                          \
        }                                                                                          \
    } while (0)

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
static const int TEST_SCALE = 1;
#else
static const int TEST_SCALE = 8;
#endif

static int rnd(std::mt19937_64& rng, int lo, int hi) { return lo + (int)(rng() % (unsigned long long)(hi - lo + 1)); }

// 愚直解（木・森専用）: 再帰 DFS で行きがけ順・帰りがけ順・オイラーツアーなどを作る
struct Brute {
    int n;
    std::vector<int> par, pe, dep, pre, post, tour, rootof;
    std::vector<Cost> pw, wdep;
    std::vector<std::vector<int>> ch;
    std::vector<int> roots;
    Brute(const Graph& G, const std::vector<int>& rs, bool all) : n(G.size()), par(n, -1), pe(n, -1), dep(n, -1), rootof(n, -1), pw(n, 0), wdep(n, 0), ch(n) {
        std::function<void(int)> dfs = [&](int v) {
            pre.push_back(v), tour.push_back(v);
            for (const Edge& e : G[v]) {
                if (dep[e.to] != -1) continue;
                dep[e.to] = dep[v] + 1, wdep[e.to] = wdep[v] + e.w, par[e.to] = v, pe[e.to] = e.id, pw[e.to] = e.w, rootof[e.to] = rootof[v];
                ch[v].push_back(e.to);
                dfs(e.to);
                tour.push_back(v);
            }
            post.push_back(v);
        };
        auto run = [&](int r) {
            if (r < 0 || r >= n || dep[r] != -1) return;
            roots.push_back(r), dep[r] = 0, rootof[r] = r;
            dfs(r);
        };
        if (all)
            for (int v = 0; v < n; ++v) run(v);
        else
            for (int r : rs) run(r);
    }
};

static Graph random_tree(std::mt19937_64& rng, int n, bool weighted, int shape) {
    // shape 0: ランダム、1: 一直線、2: 星、3: 2 分木風。辺の向き・順番もランダムにする
    std::vector<int> p(n, -1), perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    for (int i = 1; i < n; ++i) p[i] = shape == 0 ? rnd(rng, 0, i - 1) : shape == 1 ? i - 1 : shape == 2 ? 0 : (i - 1) / 2;
    std::vector<std::tuple<int, int, Cost>> es;
    for (int i = 1; i < n; ++i) {
        int a = perm[i], b = perm[p[i]];
        if (rng() & 1) std::swap(a, b);
        es.emplace_back(a, b, weighted ? (Cost)rnd(rng, 1, 20) : 1);
    }
    std::shuffle(es.begin(), es.end(), rng);
    Graph G(n);
    for (auto& [a, b, w] : es) G.add_undirected(a, b, w);
    return G;
}

static void check_common(const Graph& G, const RootedTree& T) {
    const int n = G.size();
    RT_CHECK(T.size() == n && (int)T.par.size() == n && (int)T.edge_child.size() == G.edge_count());
    // 子の一覧と親の整合、深さ、辺
    int reached = 0;
    for (int v = 0; v < n; ++v) {
        if (!T.reached(v)) {
            RT_CHECK(T.par[v] == -1 && T.depth[v] == -1 && T.in[v] == -1 && T.out[v] == -1 && T.sz[v] == 0 && T.root_of[v] == -1);
            RT_CHECK(T.children(v).empty() && T.subtree(v).empty() && !T.is_leaf(v) && T.path_to_root(v).empty());
            continue;
        }
        ++reached;
        const int p = T.par[v];
        if (p == -1) {
            RT_CHECK(T.depth[v] == 0 && T.wdepth[v] == 0 && T.par_edge[v] == -1 && T.root_of[v] == v && T.par_w[v] == 0);
        } else {
            const RawEdge& e = G.es[T.par_edge[v]];
            RT_CHECK((e.from == p && e.to == v) || (e.from == v && e.to == p));
            RT_CHECK(T.depth[v] == T.depth[p] + 1 && T.wdepth[v] == T.wdepth[p] + e.w && T.par_w[v] == e.w);
            RT_CHECK(T.root_of[v] == T.root_of[p] && T.edge_child[T.par_edge[v]] == v);
            bool found = false;
            for (int c : T.children(p)) found |= (c == v);
            RT_CHECK(found);
        }
        for (int c : T.children(v)) RT_CHECK(T.par[c] == v);
        RT_CHECK(T.child_count(v) == T.children(v).size() && T.is_leaf(v) == T.children(v).empty());
        // 子は行きがけ順に並ぶ
        const IntSpan cs = T.children(v);
        for (int k = 1; k < cs.size(); ++k) RT_CHECK(T.in[cs[k - 1]] < T.in[cs[k]]);
        if (!cs.empty()) RT_CHECK(cs.front() == cs[0] && cs.back() == cs[cs.size() - 1] && cs.to_vector().size() == (size_t)cs.size());
    }
    RT_CHECK((int)T.order.size() == reached);
    for (int id = 0; id < G.edge_count(); ++id) {
        const int c = T.edge_child[id];
        if (c != -1) RT_CHECK(T.par_edge[c] == id);
    }
    // 子孫の集合（愚直）と in / out / subtree / is_ancestor / sz
    std::vector<std::vector<char>> anc(n, std::vector<char>(n, 0));  // anc[u][v] = u が v の祖先
    for (int v = 0; v < n; ++v)
        if (T.reached(v))
            for (int u = v; u != -1; u = T.par[u]) anc[u][v] = 1;
    for (int u = 0; u < n; ++u) {
        if (!T.reached(u)) continue;
        RT_CHECK(T.order[T.in[u]] == u && T.out[u] - T.in[u] == T.sz[u] && T.subtree_size(u) == T.sz[u]);
        std::vector<int> ds;
        for (int v = 0; v < n; ++v) {
            RT_CHECK(T.is_ancestor(u, v) == (bool)anc[u][v]);
            if (anc[u][v]) ds.push_back(v);
        }
        std::vector<int> sub = T.subtree(u).to_vector();
        RT_CHECK(sub.size() == ds.size() && sub[0] == u);
        std::sort(sub.begin(), sub.end());
        RT_CHECK(sub == ds);
        const std::vector<int> pr = T.path_to_root(u);
        RT_CHECK((int)pr.size() == T.depth[u] + 1 && pr[0] == u && pr.back() == T.root_of[u]);
        for (size_t k = 1; k < pr.size(); ++k) RT_CHECK(T.par[pr[k - 1]] == pr[k]);
    }
    // 根・root
    RT_CHECK(T.root == (T.roots.empty() ? -1 : T.roots[0]));
    for (int r : T.roots) RT_CHECK(T.par[r] == -1 && T.reached(r));
    // bfs_order・postorder・leaves・euler_tour の形
    const std::vector<int> bo = T.bfs_order(), po = T.postorder(), lv = T.leaves(), et = T.euler_tour();
    RT_CHECK((int)bo.size() == reached && (int)po.size() == reached);
    for (size_t k = 1; k < bo.size(); ++k) RT_CHECK(T.depth[bo[k - 1]] <= T.depth[bo[k]]);
    std::vector<int> seen(n, 0);
    for (int v : po) {
        for (int c : T.children(v)) RT_CHECK(seen[c]);
        seen[v] = 1;
    }
    for (int v : lv) RT_CHECK(T.is_leaf(v));
    int nl = 0;
    for (int v = 0; v < n; ++v) nl += T.is_leaf(v);
    RT_CHECK((int)lv.size() == nl);
    RT_CHECK((int)et.size() == 2 * reached - (int)T.roots.size());
    // heights・wheights
    const std::vector<int> h = T.heights();
    const std::vector<Cost> wh = T.wheights();
    for (int v = 0; v < n; ++v) {
        if (!T.reached(v)) {
            RT_CHECK(h[v] == -1 && wh[v] == 0);
            continue;
        }
        int eh = 0;
        Cost ew = 0;
        for (int u : T.subtree(v)) eh = std::max(eh, T.depth[u] - T.depth[v]), ew = std::max(ew, T.wdepth[u] - T.wdepth[v]);
        RT_CHECK(h[v] == eh && wh[v] == ew);
    }
}

static void check_sums(const RootedTree& T, std::mt19937_64& rng) {
    const int n = T.size();
    std::vector<long long> val(n);
    for (auto& x : val) x = rnd(rng, -50, 50);
    const auto ss = T.subtree_sum(val);
    const auto sm = T.subtree_fold(val, [](long long a, long long b) { return std::max(a, b); });
    const auto rs = T.root_path_sum(val);
    const auto rm = T.root_path_fold(val, [](long long a, long long b) { return std::min(a, b); });
    for (int v = 0; v < n; ++v) {
        if (!T.reached(v)) {
            RT_CHECK(ss[v] == val[v] && sm[v] == val[v] && rs[v] == val[v] && rm[v] == val[v]);
            continue;
        }
        long long es = 0, em = val[v], er = 0, en = val[v];
        for (int u : T.subtree(v)) es += val[u], em = std::max(em, val[u]);
        for (int u : T.path_to_root(v)) er += val[u], en = std::min(en, val[u]);
        RT_CHECK(ss[v] == es && sm[v] == em && rs[v] == er && rm[v] == en);
    }
    const std::vector<int> vi(n, 1);
    const auto cnt = T.subtree_sum(vi);
    for (int v = 0; v < n; ++v)
        if (T.reached(v)) RT_CHECK(cnt[v] == T.sz[v]);
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    // 1. 木・森：再帰 DFS の愚直解と完全に一致する
    for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 25);
        Graph G = random_tree(rng, n, rng() & 1, rnd(rng, 0, 3));
        // ときどき辺を何本か抜いて森にする（辺を作り直す）
        if (it % 3 == 0 && G.edge_count() > 0) {
            Graph F(n);
            for (const RawEdge& e : G.es)
                if (rng() % 4) F.add_undirected(e.from, e.to, e.w);
            G = F;
        }
        if (rng() & 1) G.sort_adj();
        const int mode = rnd(rng, 0, 2);
        std::vector<int> rs;
        if (mode == 1) {
            const int k = rnd(rng, 1, 3);
            for (int i = 0; i < k; ++i) rs.push_back(rnd(rng, 0, n - 1));
        }
        const int r0 = rnd(rng, 0, n - 1);
        const RootedTree T = mode == 0 ? RootedTree(G, r0) : mode == 1 ? RootedTree(G, rs) : RootedTree::forest(G);
        const Brute B(G, mode == 0 ? std::vector<int>{r0} : rs, mode == 2);
        RT_CHECK(T.order == B.pre && T.par == B.par && T.par_edge == B.pe && T.depth == B.dep && T.wdepth == B.wdep && T.par_w == B.pw);
        RT_CHECK(T.roots == B.roots && T.root_of == B.rootof && T.postorder() == B.post && T.euler_tour() == B.tour);
        for (int v = 0; v < n; ++v) RT_CHECK(T.children(v).to_vector() == B.ch[v]);
        if (mode == 2) RT_CHECK((int)T.order.size() == n);
        check_common(G, T);
        check_sums(T, rng);
    }
    // 2. 木でないグラフ（無向・有向、多重辺・自己ループあり）：全域木としての性質を確かめる
    for (int it = 0; it < 800 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 20), m = rnd(rng, 0, 3 * n);
        const bool und = rng() & 1;
        Graph G(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            const Cost w = rnd(rng, 1, 9);
            und ? G.add_undirected(a, b, w) : G.add_edge(a, b, w);
        }
        const bool all = rng() & 1;
        const int r = rnd(rng, 0, n - 1);
        const RootedTree T = all ? RootedTree::forest(G) : RootedTree(G, r);
        check_common(G, T);
        check_sums(T, rng);
        // 到達できる頂点の集合が正しい
        std::vector<char> vis(n, 0);
        std::vector<int> q;
        if (all) {
            RT_CHECK((int)T.order.size() == n || !und);
        } else {
            vis[r] = 1, q.push_back(r);
            for (size_t h = 0; h < q.size(); ++h)
                for (const Edge& e : G[q[h]])
                    if (!vis[e.to]) vis[e.to] = 1, q.push_back(e.to);
            for (int v = 0; v < n; ++v) RT_CHECK(T.reached(v) == (bool)vis[v]);
        }
    }
    // 3. 端のケース
    {
        const Graph E;
        const RootedTree T0(E, 0);
        RT_CHECK(T0.size() == 0 && T0.root == -1 && T0.roots.empty() && T0.order.empty() && T0.euler_tour().empty());
        const RootedTree T1(Graph(1), 0);
        RT_CHECK(T1.root == 0 && T1.order == std::vector<int>{0} && T1.is_leaf(0) && T1.euler_tour() == std::vector<int>{0});
        const RootedTree T2(Graph(3), 5);  // 範囲外の根は無視
        RT_CHECK(T2.roots.empty() && T2.order.empty());
        const RootedTree T3(Graph(3), std::vector<int>{2, 2, 0});  // 同じ根を 2 回
        RT_CHECK(T3.roots == (std::vector<int>{2, 0}) && T3.order == (std::vector<int>{2, 0}));
    }
    std::cout << "[rooted_tree] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[rooted_tree] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        for (int shape : {0, 1}) {
            Graph G = random_tree(rng, n, true, shape);
            auto t0 = std::chrono::steady_clock::now();
            RootedTree T(G, 0);
            const long long tb = ms(t0);
            t0 = std::chrono::steady_clock::now();
            std::vector<long long> val(n, 1);
            const auto s = T.subtree_sum(val);
            const auto p = T.root_path_sum(val);
            const auto et = T.euler_tour();
            const auto po = T.postorder();
            long long chk = s[0] + p[n - 1] + (long long)et.size() + po[0];
            for (int k = n - 1; k >= 0; --k) {  // 木DP の型（子をなめる）
                const int v = T.order[k];
                for (int c : T.children(v)) chk += c & 1;
            }
            const long long tq = ms(t0);
            std::cout << "[rooted_tree] benchmark: N = " << n << (shape ? "（一直線）" : "（ランダム）") << "  build " << tb
                      << " ms, subtree_sum+root_path_sum+euler_tour+postorder+木DPの型 " << tq << " ms  (check " << chk << ")\n";
        }
    }
}

}  // namespace rt_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace rt_demo {

static std::string to_s(bool b) { return b ? "true" : "false"; }
template <class T, std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool>, int> = 0>
static std::string to_s(T v) { return std::to_string(v); }
template <class T>
static std::string to_s(const std::vector<T>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
static std::string to_s(const IntSpan& v) { return to_s(v.to_vector()); }
static std::string pad(const std::string& s, size_t w) {
    size_t cols = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = (unsigned char)s[i];
        if ((c & 0xC0) == 0x80) continue;
        cols += (c >= 0xE0 ? 2 : 1);
    }
    return cols >= w ? s : s + std::string(w - cols, ' ');
}
static void banner(const std::string& s) {
    std::cout << "\n=====================================================================\n"
              << "  " << s << "\n"
              << "=====================================================================\n";
}
static void api(const std::string& usage, const std::string& ret, const std::string& what, const std::string& cost) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << usage << "\n"
              << "    返り値 : " << ret << "\n"
              << "    内容   : " << what << "\n"
              << "    計算量 : " << cost << "\n";
}
static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
static void stmt(const std::string& s) {
    const bool block = !s.empty() && s.back() == '}' &&
                       (s.compare(0, 4, "for ") == 0 || s.compare(0, 3, "if ") == 0 || s.compare(0, 6, "while ") == 0);
    code(block ? s : s + ";");
}
template <class V>
static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 30) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 58) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) rt_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    rt_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void cheat_sheet() {
    struct Row { const char* use; const char* ret; const char* what; };
    const Row rows[] = {
        {"RootedTree T(G, root);", "", "根を決めて作る（O(N)）。T(G, {r1, r2}) は森、RootedTree::forest(G) は全成分"},
        {"T.par[v] / T.par_edge[v] / T.par_w[v]", "int / int / Cost", "親 / 親への辺の番号 / 重み（根は -1, -1, 0）"},
        {"T.depth[v] / T.wdepth[v]", "int / Cost", "根からの辺の本数 / 重みの和"},
        {"T.sz[v]", "int", "部分木の頂点数"},
        {"T.order", "vector<int>", "行きがけ順。逆順に回すと子が先（木DP）"},
        {"T.in[v], T.out[v]", "int", "部分木 = order[in[v] .. out[v])。BIT・セグ木で部分木クエリ"},
        {"T.children(v)", "IntSpan", "子の一覧。for (int c : T.children(v))"},
        {"T.subtree(v)", "IntSpan", "部分木の頂点（先頭が v）"},
        {"T.is_ancestor(u, v)", "bool", "u が v の祖先か（u == v も true）  O(1)"},
        {"T.edge_child[id]", "int", "木の辺 id の子の側（辺を切ると sz[c] と N - sz[c] に分かれる）"},
        {"T.root / T.roots / T.root_of[v]", "int / vector / int", "根 / 根の一覧 / 属する木の根"},
        {"T.path_to_root(v)", "vector<int>", "v から根まで"},
        {"T.leaves() / bfs_order() / postorder()", "vector<int>", "葉 / 根に近い順 / 帰りがけ順"},
        {"T.euler_tour()", "vector<int>", "入るときと子から戻るたびに書く列（2N - 1 個）"},
        {"T.heights() / T.wheights()", "vector", "部分木の中で最も深い頂点までの距離"},
        {"T.subtree_sum(val) / subtree_fold(val, op)", "vector<T>", "部分木の値の和 / op でまとめる"},
        {"T.root_path_sum(val) / root_path_fold(val, op)", "vector<T>", "根から v までの和 / op でまとめる"},
    };
    std::cout << "\n■ API 一覧（頂点は 0-indexed。根から届かない頂点は「未到達」で depth = -1）\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, 48) << pad(ret, 22) << r.what << "\n";
    }
}

static void demo() {
    banner("RootedTree 使い方デモ    見方:  コード  ->  返り値    // 意味");
    cheat_sheet();

    head("デモで使う木（辺の重みつき）");
    note("        0");
    note("      /   \\          辺: 0-1 (2), 0-2 (5), 1-3 (1), 1-4 (3), 2-5 (4)");
    note("     1     2");
    note("    / \\     \\");
    note("   3   4     5");
    {
        std::istringstream ss("6\n1 2 2\n1 3 5\n2 4 1\n2 5 3\n3 6 4\n");
        std::streambuf* old = std::cin.rdbuf(ss.rdbuf());
        code("// 入力: 6 / 1 2 2 / 1 3 5 / 2 4 1 / 2 5 3 / 3 6 4");
        RUN(int N; std::cin >> N);
        RUN(Graph G = read_tree_w(N));
        std::cin.rdbuf(old);

        api("【1】 RootedTree T(G, root);", "", "根を決めて、親・深さ・部分木などをまとめて作る。非再帰なので深い木でも落ちない", "O(N)");
        RUN(RootedTree T(G, 0));
        SHOW(T.par, "親（根は -1）");
        SHOW(T.par_edge, "親への辺の番号（入力の何本目か）");
        SHOW(T.par_w, "親への辺の重み");
        SHOW(T.depth, "根からの辺の本数");
        SHOW(T.wdepth, "根からの重みの和（= 根からの距離）");
        SHOW(T.sz, "部分木の頂点数");

        api("【2】 T.order / T.in[v] / T.out[v] / T.subtree(v)", "vector<int> / int / IntSpan",
            "行きがけ順。v の部分木は order の区間 [in[v], out[v])。部分木の和などを BIT・セグ木で求められる", "O(1)");
        SHOW(T.order, "行きがけ順（子は隣接の順）");
        SHOW(T.in, "order の何番目か");
        SHOW(T.out, "部分木の終わり（半開区間の右端）");
        SHOW(T.subtree(1), "頂点 1 の部分木 = order[1 .. 4)");
        SHOW(T.is_ancestor(1, 4), "1 は 4 の祖先");
        SHOW(T.is_ancestor(2, 4), "2 は 4 の祖先ではない");
        note("部分木に x を足す・部分木の和 → BIT の区間 [in[v], out[v]) に対する操作になる。");
        note("根から v までの値の和（点更新あり）→ BIT に +x を in[u]、-x を out[u] に足し、in[v] までの和を取る。");

        api("【3】 T.children(v) / T.child_count(v) / T.is_leaf(v) / T.leaves()", "IntSpan / int / bool / vector<int>",
            "子の一覧（コピーせずに参照する）。for (int c : T.children(v)) で回す", "O(1)（leaves は O(N)）");
        SHOW(T.children(1), "頂点 1 の子");
        SHOW(T.child_count(0), "根の子の数");
        SHOW(T.is_leaf(5), "5 は葉");
        SHOW(T.leaves(), "葉の一覧");

        api("【4】 T.edge_child[id]", "int", "木の辺 id の子の側の頂点。辺を切ると sz[c] 個と N - sz[c] 個に分かれる", "O(1)");
        SHOW(T.edge_child, "辺 0 (0-1) の子は 1、辺 2 (1-3) の子は 3 …");
        code("// 各辺を通る頂点の組の数 = sz[c] * (N - sz[c])");

        api("【5】 T.path_to_root(v) / T.bfs_order() / T.postorder() / T.euler_tour()", "vector<int>",
            "根までの頂点 / 根に近い順 / 帰りがけ順（子 → 親）/ オイラーツアー（入るときと戻るたびに書く。ABC213 D）", "O(深さ) / O(N)");
        SHOW(T.path_to_root(4), "4 → 1 → 0");
        SHOW(T.bfs_order(), "深さの順");
        SHOW(T.postorder(), "子を全部見てから親");
        SHOW(T.euler_tour(), "長さ 2N - 1");

        api("【6】 T.heights() / T.wheights()", "vector<int> / vector<Cost>", "部分木の中で、v から最も深い頂点までの辺の本数 / 重みの和", "O(N)");
        SHOW(T.heights(), "葉は 0");
        SHOW(T.wheights(), "重み付き");

        api("【7】 T.subtree_sum(val) / T.subtree_fold(val, op)", "vector<T>", "部分木の値の和 / 部分木の値を op（max・min・gcd など）でまとめる", "O(N)");
        RUN(std::vector<long long> val = {1, 10, 100, 1000, 10000, 100000});
        SHOW(T.subtree_sum(val), "部分木の和");
        RUN(auto mx = T.subtree_fold(val, [](long long a, long long b) { return std::max(a, b); }));
        SHOW(mx, "部分木の最大値");

        api("【8】 T.root_path_sum(val) / T.root_path_fold(val, op)", "vector<T>", "根から v までの値の和 / op でまとめる", "O(N)");
        SHOW(T.root_path_sum(val), "根から v までの和");
        note("ABC138 D（部分木に x を足す操作を Q 回）: add[p] += x を全部ためて、最後に root_path_sum(add)。");

        api("【9】 RootedTree T(G, {r1, r2, ...}) / RootedTree::forest(G)", "RootedTree", "森（連結でないグラフ）。roots に根の一覧、root_of[v] に v の属する木の根", "O(N + M)");
        RUN(Graph F(5));
        RUN(F.add_undirected(0, 1); F.add_undirected(3, 2));
        RUN(RootedTree TF = RootedTree::forest(F));
        SHOW(TF.roots, "成分ごとに番号最小の頂点が根");
        SHOW(TF.root_of, "属する木の根");
        RUN(RootedTree T1(F, 0));
        SHOW(T1.depth, "根 0 から届かない頂点は -1（未到達）");

        head("【10】 木DP の型（order を逆順に回すと、子が必ず先に計算済み）");
        code("// 例: 最大独立集合（隣り合わない頂点の最大個数）");
        code("std::vector<std::array<long long, 2>> dp(N);   // dp[v][0] = v を選ばない, dp[v][1] = 選ぶ");
        code("for (int k = N - 1; k >= 0; --k) {");
        code("    int v = T.order[k];");
        code("    dp[v] = {0, 1};");
        code("    for (int c : T.children(v)) {");
        code("        dp[v][0] += std::max(dp[c][0], dp[c][1]);");
        code("        dp[v][1] += dp[c][0];");
        code("    }");
        code("}");
        std::vector<std::array<long long, 2>> dp(N);
        for (int k = N - 1; k >= 0; --k) {
            const int v = T.order[k];
            dp[v] = {0, 1};
            for (int c : T.children(v)) {
                dp[v][0] += std::max(dp[c][0], dp[c][1]);
                dp[v][1] += dp[c][0];
            }
        }
        std::cout << "    答え = max(dp[0][0], dp[0][1]) = " << std::max(dp[0][0], dp[0][1]) << "（例: 0, 3, 4, 5）\n";
        note("根から下へ配る DP（親の値を使う）は order を前から回す。");

        head("【11】 例: 木の辺の彩色（ABC146 D）: 親への辺と違う色を、子に小さい順に配る");
        std::vector<int> color(G.edge_count(), 0);
        int K = 0;
        for (int v : T.order) {
            const int pc = (T.par_edge[v] == -1 ? 0 : color[T.par_edge[v]]);
            int c = 1;
            for (int ch : T.children(v)) {
                if (c == pc) ++c;
                color[T.par_edge[ch]] = c++;
            }
            K = std::max(K, T.child_count(v) + (T.par[v] != -1));
        }
        code("for (int v : T.order) {  int pc = (T.par_edge[v] == -1 ? 0 : color[T.par_edge[v]]);  int c = 1;");
        code("    for (int ch : T.children(v)) { if (c == pc) ++c;  color[T.par_edge[ch]] = c++; } }");
        std::cout << "    色の数 = " << K << ", 各辺の色 = " << to_s(color) << "\n";
    }

    head("典型的な使い方");
    note("木の入力           : Graph G = read_tree(N);  RootedTree T(G, 0);");
    note("親の配列の入力     : Graph G = read_parents(N);  RootedTree T(G, 0);");
    note("部分木クエリ       : BIT で [T.in[v], T.out[v]) を操作");
    note("2 点間の距離・LCA  : 05_lca の LCA L(T);");
    note("直径               : 05_tree_diameter の tree_diameter(G)");
    note("類題: ABC138 D（root_path_sum）、ABC146 D（辺の彩色）、ABC213 D（euler_tour）、ABC220 F（sz と全方位）");
}

}  // namespace rt_demo


#ifdef _WIN32
extern "C" __declspec(dllimport) int __stdcall SetConsoleOutputCP(unsigned int);
#endif
static void setup_console() {
#ifdef _WIN32
    const char* probe = "あ";  // UTF-8 なら 0xE3 0x81 0x82、CP932 なら 0x82 0xA0
    SetConsoleOutputCP((unsigned char)probe[0] == 0xE3 ? 65001u : 932u);
#endif
}

int main() {
    setup_console();
    rt_test::random_test();
    rt_test::benchmark();
    rt_demo::demo();
    return 0;
}
