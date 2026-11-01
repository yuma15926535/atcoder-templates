// =============================================================================
//  グラフ基盤（Graph / 入力 / 変換 / 探索の結果 SearchResult）   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 00 基盤（ほかのグラフ系テンプレートはすべてこれを同梱している）
//  依存       : なし
//  頂点番号   : 0-indexed（入力の 1-indexed は read 系の関数が 0-indexed に直す）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  共通の約束（第二世代のグラフ・グリッド系テンプレートすべてで同じ）
//    - 重み・距離の型は Cost（= long long）。到達不能は DIST_INF（= 1LL << 60 = Main_base の LINF）
//    - 「なし」は -1（頂点・親・辺番号）、経路がないときは空の vector、区間は半開区間 [l, r)
//    - DFS 系はすべて非再帰（Visual Studio のスタック 1MB でも落ちない）
//
//  高速化のポイント
//    - 辺は 16 バイトの構造体 Edge{to, id, w} 1 種類（重みなしも同じ型）
//    - read 系・変換系は辺を全部読んでから各頂点の容量を確保して詰める（vector の再確保なし）
//    - parent_edges は全辺を 1 回なめるだけ（O(N + M)。星形グラフでも遅くならない）
//
//  機能一覧                                                                    計算量
//    Graph G(n)                       頂点 0 ~ n-1、辺なしで作る                O(n)
//    G.add_edge(u, v, w = 1)  -> int  有向辺 u→v を追加し、辺番号を返す         ならし O(1)
//    G.add_undirected(u, v, w = 1)    無向辺 u-v（両方向に同じ辺番号）          ならし O(1)
//    G.add(u, v, w, DIRECTED / UNDIRECTED)   向きを引数で選ぶ版                 ならし O(1)
//    G.add_vertex()           -> int  頂点を 1 つ増やし番号を返す（超頂点用）   ならし O(1)
//    G[v]                             v から出る辺の一覧 vector<Edge>（e.to, e.w, e.id）
//    G.es[id]                         辺 id の {from, to, w}（追加した順）
//    G.size() / G.edge_count()        頂点数 / 辺の本数                         O(1)
//    G.reversed()                     全辺の向きを逆にしたグラフ（辺番号は同じ） O(N + M)
//    G.sort_adj()                     隣接を (行き先, 重み, 辺番号) の昇順に     O(M log M)
//    G.in_degrees() / G.out_degrees() 入次数 / 出次数の配列                     O(N + M)
//    G.adj_list() / G.adj_list_w()    vector<vector<int>> / {行き先, 重み} の形に変換  O(N + M)
//    read_undirected(n, m) / read_directed(n, m)          重みなしの辺 m 本を読む
//    read_undirected_w(n, m) / read_directed_w(n, m)      重みありの辺 m 本を読む（u v w）
//    read_tree(n) / read_tree_w(n)                        木（辺 n-1 本）を読む
//    read_parents(n)                                      p_2 .. p_n（親の番号）の形の木を読む
//        ※ どれも最後の引数 base（既定 1）で入力の番号の始まりを指定できる     O(N + M)
//    graph_from_adj(adj)              vector<vector<int>> / vector<vector<pair<int,W>>> から
//    graph_from_edges(n, us, vs, [ws,] DIRECTED / UNDIRECTED, base)   辺の配列から
//    graph_from_edges(n, uv_pairs, DIRECTED / UNDIRECTED, base)       pair の配列から
//    graph_from_parents(par, [w,] base)   親の配列から木（森）。par - base が負なら親なし
//    SearchResult r                   探索の結果。r.dist[v]（距離）, r.par[v]（直前の頂点）
//    r.reached(v) / r.path(t)         到達したか / 始点 → t の頂点列（到達不能なら空）
//    r.nearest_source()               各頂点に届いた始点（多始点で「いちばん近い始点」） O(N)
//    auto [dist, par] = bfs(G, s);    構造化束縛で分けて受け取れる
//    restore_path(par, t)             自分で作った親の配列から経路を復元        O(経路長)
//    parent_edges(G, r)               各頂点に入るときに使った辺の番号           O(N + M)
//    path_edges(G, r, t)              始点 → t の経路の辺番号の列                O(N + M)
//    cout << e / cout << G            デバッグ表示（Main_base の debug(G) でも使える）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    gb_test::random_test()  … 愚直な実装と全機能を突き合わせるランダムテスト
//    gb_test::benchmark()    … 速度計測
//    gb_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  テスト用の小道具（std::cin を文字列に差し替える）
// =============================================================================
struct CinFromString {
    std::istringstream ss;
    std::streambuf* old;
    explicit CinFromString(const std::string& s) : ss(s), old(std::cin.rdbuf(ss.rdbuf())) {}
    ~CinFromString() { std::cin.rdbuf(old); }
};

// =============================================================================
//  ランダムテスト
// =============================================================================
namespace gb_test {

#define GB_CHECK(...)                                                                           \
    do {                                                                                     \
        if (!(__VA_ARGS__)) {                                                                \
            std::cerr << "[graph_base] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                    \
        }                                                                                    \
    } while (0)

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
static const int TEST_SCALE = 1;
#else
static const int TEST_SCALE = 8;
#endif

static int rnd(std::mt19937_64& rng, int lo, int hi) { return lo + (int)(rng() % (unsigned long long)(hi - lo + 1)); }

struct E3 { int u, v; Cost w; };

static bool same_adj(const Graph& a, const Graph& b) {
    if (a.size() != b.size() || a.edge_count() != b.edge_count()) return false;
    for (int v = 0; v < a.size(); ++v) {
        if (a[v].size() != b[v].size()) return false;
        for (size_t k = 0; k < a[v].size(); ++k)
            if (a[v][k].to != b[v][k].to || a[v][k].id != b[v][k].id || a[v][k].w != b[v][k].w) return false;
    }
    for (int i = 0; i < a.edge_count(); ++i)
        if (a.es[i].from != b.es[i].from || a.es[i].to != b.es[i].to || a.es[i].w != b.es[i].w) return false;
    return true;
}
static std::vector<std::tuple<int, int, int, Cost>> edge_multiset(const Graph& g) {
    std::vector<std::tuple<int, int, int, Cost>> r;
    for (int v = 0; v < g.size(); ++v)
        for (const Edge& e : g[v]) r.emplace_back(v, e.to, e.id, e.w);
    std::sort(r.begin(), r.end());
    return r;
}

// テスト用の多始点ダイクストラ（O(N^2)。同じ距離なら番号の小さい頂点を先に確定）
static SearchResult naive_dijkstra(const Graph& G, const std::vector<int>& starts) {
    const int n = G.size();
    SearchResult r{std::vector<Cost>(n, DIST_INF), std::vector<int>(n, -1)};
    std::vector<char> done(n, 0);
    for (int s : starts) r.dist[s] = 0;
    for (int it = 0; it < n; ++it) {
        int u = -1;
        for (int v = 0; v < n; ++v)
            if (!done[v] && r.dist[v] != DIST_INF && (u == -1 || r.dist[v] < r.dist[u])) u = v;
        if (u == -1) break;
        done[u] = 1;
        for (const Edge& e : G[u])
            if (r.dist[u] + e.w < r.dist[e.to]) r.dist[e.to] = r.dist[u] + e.w, r.par[e.to] = u;
    }
    return r;
}

static void check_search_result(const Graph& G, const SearchResult& r, const std::vector<int>& starts, bool undirected) {
    const int n = G.size();
    std::vector<char> is_start(n, 0);
    for (int s : starts) is_start[s] = 1;
    const std::vector<int> src = r.nearest_source();
    const std::vector<int> pe = parent_edges(G, r);
    GB_CHECK((int)src.size() == n && (int)pe.size() == n);
    for (int t = 0; t < n; ++t) {
        const std::vector<int> p = r.path(t);
        if (!r.reached(t)) {
            GB_CHECK(p.empty() && src[t] == -1 && pe[t] == -1 && path_edges(G, r, t).empty());
            continue;
        }
        GB_CHECK(!p.empty() && p.back() == t && is_start[p[0]] && r.dist[p[0]] == 0);
        GB_CHECK(src[t] == p[0]);
        GB_CHECK(restore_path(r.par, t) == p);
        const std::vector<int> pes = path_edges(G, r, t);
        GB_CHECK(pes.size() + 1 == p.size());
        Cost sum = 0;
        for (size_t k = 0; k + 1 < p.size(); ++k) {
            const RawEdge& e = G.es[pes[k]];
            const bool ok = (e.from == p[k] && e.to == p[k + 1]) || (undirected && e.from == p[k + 1] && e.to == p[k]);
            GB_CHECK(ok);
            sum += e.w;
        }
        GB_CHECK(sum == r.dist[t]);
        if (r.par[t] == -1) GB_CHECK(pe[t] == -1);
        else GB_CHECK(r.dist[r.par[t]] + G.es[pe[t]].w == r.dist[t]);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    const int iters = 600 * TEST_SCALE;
    for (int it = 0; it < iters; ++it) {
        const int n = rnd(rng, 1, 12 + (it % 5) * 6);
        const int m = rnd(rng, 0, 3 * n);
        const bool und = rng() & 1;
        const bool wtd = rng() & 1;
        std::vector<E3> es(m);
        for (auto& e : es) e = {rnd(rng, 0, n - 1), rnd(rng, 0, n - 1), wtd ? (Cost)rnd(rng, 0, 20) : 1};

        // 1. add_edge / add_undirected / add と read 系・変換系が同じグラフを作る
        Graph A(n), B(n);
        for (int i = 0; i < m; ++i) {
            const int id = und ? A.add_undirected(es[i].u, es[i].v, es[i].w) : A.add_edge(es[i].u, es[i].v, es[i].w);
            GB_CHECK(id == i);
            GB_CHECK(B.add(es[i].u, es[i].v, es[i].w, und ? UNDIRECTED : DIRECTED) == i);
        }
        GB_CHECK(A.size() == n && A.edge_count() == m && same_adj(A, B));
        {
            const int base = rnd(rng, 0, 1);
            std::string in;
            for (auto& e : es) {
                in += std::to_string(e.u + base) + " " + std::to_string(e.v + base);
                if (wtd) in += " " + std::to_string(e.w);
                in += "\n";
            }
            CinFromString cf(in);
            Graph R = wtd ? (und ? read_undirected_w(n, m, base) : read_directed_w(n, m, base))
                          : (und ? read_undirected(n, m, base) : read_directed(n, m, base));
            GB_CHECK(same_adj(A, R));
        }
        {
            const int base = rnd(rng, 0, 1);
            std::vector<int> us(m), vs(m);
            std::vector<long long> ws(m);
            std::vector<std::pair<int, int>> uv(m);
            for (int i = 0; i < m; ++i) us[i] = es[i].u + base, vs[i] = es[i].v + base, ws[i] = es[i].w, uv[i] = {us[i], vs[i]};
            const EdgeDir d = und ? UNDIRECTED : DIRECTED;
            if (wtd) GB_CHECK(same_adj(A, graph_from_edges(n, us, vs, ws, d, base)));
            else {
                GB_CHECK(same_adj(A, graph_from_edges(n, us, vs, d, base)));
                GB_CHECK(same_adj(A, graph_from_edges(n, uv, d, base)));
            }
        }
        // 2. graph_from_adj（有向として取り込む）
        {
            std::vector<std::vector<int>> a(n);
            std::vector<std::vector<std::pair<int, long long>>> aw(n);
            for (int v = 0; v < n; ++v)
                for (const Edge& e : A[v]) a[v].push_back(e.to), aw[v].push_back({e.to, e.w});
            const Graph F = graph_from_adj(a), FW = graph_from_adj(aw);
            int id = 0;
            for (int v = 0; v < n; ++v) {
                GB_CHECK(F[v].size() == A[v].size() && FW[v].size() == A[v].size());
                for (size_t k = 0; k < A[v].size(); ++k, ++id) {
                    GB_CHECK(F[v][k].to == A[v][k].to && F[v][k].w == 1 && F[v][k].id == id);
                    GB_CHECK(FW[v][k].to == A[v][k].to && FW[v][k].w == A[v][k].w && FW[v][k].id == id);
                    GB_CHECK(F.es[id].from == v && F.es[id].to == A[v][k].to);
                }
            }
            GB_CHECK(F.edge_count() == id && F.adj_list() == a && FW.adj_list_w() == aw && A.adj_list() == a);
        }
        // 3. reversed / sort_adj / 次数
        {
            const Graph R = A.reversed();
            GB_CHECK(R.size() == n && R.edge_count() == m);
            std::vector<std::tuple<int, int, int, Cost>> ex;
            for (int v = 0; v < n; ++v)
                for (const Edge& e : A[v]) ex.emplace_back(e.to, v, e.id, e.w);
            std::sort(ex.begin(), ex.end());
            GB_CHECK(edge_multiset(R) == ex);
            for (int i = 0; i < m; ++i) GB_CHECK(R.es[i].from == A.es[i].to && R.es[i].to == A.es[i].from && R.es[i].w == A.es[i].w);
            if (und) GB_CHECK(edge_multiset(R) == edge_multiset(A));
            GB_CHECK(edge_multiset(R.reversed()) == edge_multiset(A));
            Graph S = A;
            S.sort_adj();
            GB_CHECK(edge_multiset(S) == edge_multiset(A));
            for (int v = 0; v < n; ++v)
                for (size_t k = 1; k < S[v].size(); ++k) {
                    const Edge &x = S[v][k - 1], &y = S[v][k];
                    GB_CHECK(std::make_tuple(x.to, x.w, x.id) <= std::make_tuple(y.to, y.w, y.id));  // 無向の自己ループは同じものが 2 つ
                }
            std::vector<int> indeg(n, 0), outdeg(n, 0);
            for (auto& e : es) {
                ++outdeg[e.u], ++indeg[e.v];
                if (und) ++outdeg[e.v], ++indeg[e.u];
            }
            GB_CHECK(A.in_degrees() == indeg && A.out_degrees() == outdeg);
        }
        // 4. add_vertex
        {
            Graph C = A;
            const int x = C.add_vertex();
            GB_CHECK(x == n && C.size() == n + 1 && C[x].empty());
            C.add_edge(x, 0, 0);
            GB_CHECK(C[x].size() == 1 && C.es.back().from == x);
        }
        // 5. SearchResult / parent_edges / path_edges / nearest_source / restore_path
        {
            std::vector<int> starts;
            const int k = rnd(rng, 1, std::min(n, 3));
            for (int i = 0; i < k; ++i) starts.push_back(rnd(rng, 0, n - 1));
            std::sort(starts.begin(), starts.end());
            starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
            const SearchResult r = naive_dijkstra(A, starts);
            check_search_result(A, r, starts, und);
            auto [dist, par] = naive_dijkstra(A, starts);  // 構造化束縛
            GB_CHECK(dist == r.dist && par == r.par);
        }
        // 6. 木の入力：read_parents / graph_from_parents
        {
            std::vector<int> p(n, -1);
            for (int i = 1; i < n; ++i) p[i] = rnd(rng, 0, i - 1);
            std::string in;
            for (int i = 1; i < n; ++i) in += std::to_string(p[i] + 1) + (i + 1 < n ? " " : "\n");
            CinFromString cf(in);
            const Graph T = read_parents(n, 1);
            GB_CHECK(T.size() == n && T.edge_count() == n - 1);
            for (int i = 1; i < n; ++i) GB_CHECK(T.es[i - 1].from == p[i] && T.es[i - 1].to == i && T.es[i - 1].w == 1);
            std::vector<long long> p1(n), w(n);
            for (int i = 0; i < n; ++i) p1[i] = (i == 0 ? (rng() & 1 ? 0 : -1) : p[i] + 1), w[i] = rnd(rng, 1, 9);
            const Graph T1 = graph_from_parents(p1, 1), T2 = graph_from_parents(p1, w, 1);
            GB_CHECK(same_adj(T, T1));
            GB_CHECK(T2.edge_count() == n - 1);
            for (int i = 1; i < n; ++i) GB_CHECK(T2.es[i - 1].from == p[i] && T2.es[i - 1].to == i && T2.es[i - 1].w == w[i]);
            std::string in2;
            for (int i = 1; i < n; ++i) in2 += std::to_string(i + 1) + " " + std::to_string(p[i] + 1) + " " + std::to_string(w[i]) + "\n";
            CinFromString cf2(in2);
            const Graph T3 = read_tree_w(n);
            GB_CHECK(T3.edge_count() == n - 1);
            for (int i = 1; i < n; ++i) GB_CHECK(T3.es[i - 1].from == i && T3.es[i - 1].to == p[i] && T3.es[i - 1].w == w[i]);
        }
    }
    // 7. 端のケース
    {
        Graph E;
        GB_CHECK(E.size() == 0 && E.edge_count() == 0 && E.reversed().size() == 0 && E.adj_list().empty());
        const Graph one = graph_from_edges(1, std::vector<int>{}, std::vector<int>{}, UNDIRECTED, 0);
        GB_CHECK(one.size() == 1 && one.edge_count() == 0);
        CinFromString cf("");
        const Graph t1 = read_tree(1);
        GB_CHECK(t1.size() == 1 && t1.edge_count() == 0);
        Graph L(2);
        L.add_undirected(1, 1, 5);  // 無向の自己ループは隣接に 2 回入る
        GB_CHECK(L[1].size() == 2 && L.in_degrees()[1] == 2 && L.reversed()[1].size() == 2);
    }
    std::cout << "[graph_base] random_test: OK (" << iters << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_base] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        const int m = n;
        std::string in;
        in.reserve((size_t)m * 24);
        std::vector<int> us(m), vs(m);
        std::vector<long long> ws(m);
        for (int i = 0; i < m; ++i) {
            us[i] = (int)(rng() % n), vs[i] = (int)(rng() % n), ws[i] = (long long)(rng() % 1000000000);
            in += std::to_string(us[i] + 1) + " " + std::to_string(vs[i] + 1) + " " + std::to_string(ws[i]) + "\n";
        }
        auto ms = [](std::chrono::steady_clock::time_point a) {
            return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
        };
        std::ios::sync_with_stdio(false);
        auto t0 = std::chrono::steady_clock::now();
        Graph G;
        {
            CinFromString cf(in);
            G = read_undirected_w(n, m);
        }
        const long long t_read = ms(t0);
        t0 = std::chrono::steady_clock::now();
        Graph H = graph_from_edges(n, us, vs, ws, UNDIRECTED, 0);
        const long long t_build = ms(t0);
        t0 = std::chrono::steady_clock::now();
        Graph R = G.reversed();
        const long long t_rev = ms(t0);
        t0 = std::chrono::steady_clock::now();
        Graph S = G;
        S.sort_adj();
        const long long t_sort = ms(t0);
        // 適当な探索木を作って parent_edges
        SearchResult r{std::vector<Cost>(n, DIST_INF), std::vector<int>(n, -1)};
        {
            std::vector<int> q{0};
            r.dist[0] = 0;
            for (size_t h = 0; h < q.size(); ++h)
                for (const Edge& e : G[q[h]])
                    if (r.dist[e.to] == DIST_INF) r.dist[e.to] = r.dist[q[h]] + 1, r.par[e.to] = q[h], q.push_back(e.to);
        }
        t0 = std::chrono::steady_clock::now();
        const std::vector<int> pe = parent_edges(G, r);
        const std::vector<int> src = r.nearest_source();
        const long long t_pe = ms(t0);
        std::cout << "[graph_base] benchmark: N = M = " << n << "  read_undirected_w(cin) " << t_read << " ms, graph_from_edges "
                  << t_build << " ms, reversed " << t_rev << " ms, sort_adj " << t_sort << " ms, parent_edges+nearest_source "
                  << t_pe << " ms  (check " << H.edge_count() + R.edge_count() + pe[n - 1] + src[n - 1] << ")\n";
    }
}

}  // namespace gb_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace gb_demo {

static std::string to_s(bool b) { return b ? "true" : "false"; }
static std::string to_s(int v) { return std::to_string(v); }
static std::string to_s(long long v) { return v == DIST_INF ? std::string("DIST_INF") : std::to_string(v); }
static std::string to_s(const std::vector<int>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + std::to_string(v[i]);
    return s + "}";
}
static std::string to_s(const std::vector<long long>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
static std::string to_s(const Edge& e) { return "{to=" + std::to_string(e.to) + ", id=" + std::to_string(e.id) + ", w=" + std::to_string(e.w) + "}"; }
static std::string to_s(const RawEdge& e) { return "{from=" + std::to_string(e.from) + ", to=" + std::to_string(e.to) + ", w=" + std::to_string(e.w) + "}"; }
static std::string to_s(const std::vector<Edge>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
static std::string to_s(const std::vector<std::vector<int>>& g) {
    std::string s = "{";
    for (size_t i = 0; i < g.size(); ++i) s += (i ? ", " : "") + to_s(g[i]);
    return s + "}";
}
static std::string pad(const std::string& s, size_t w) {
    size_t cols = 0;  // 全角を 2 文字分として数える
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
    std::string line = "  " + pad(expr, 34) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 52) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) gb_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    gb_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__
static void print_graph(const char* name, const Graph& G) {
    for (int v = 0; v < G.size(); ++v) std::cout << "    " << name << "[" << v << "] = " << to_s(G[v]) << "\n";
}

static void cheat_sheet() {
    struct Row { const char* use; const char* ret; const char* what; };
    const Row rows[] = {
        {"Graph G(n);", "", "頂点 0 ~ n-1・辺なし  O(n)"},
        {"G.add_edge(u, v, w = 1)", "int", "有向辺 u→v を追加。辺番号を返す"},
        {"G.add_undirected(u, v, w = 1)", "int", "無向辺 u-v を追加（両方向に同じ辺番号）"},
        {"G.add(u, v, w, UNDIRECTED)", "int", "向きを引数で選ぶ版（DIRECTED / UNDIRECTED）"},
        {"G.add_vertex()", "int", "頂点を 1 つ増やす（超頂点）"},
        {"G[v]", "vector<Edge>", "v から出る辺。e.to / e.w / e.id"},
        {"G.es[id]", "RawEdge", "辺 id の {from, to, w}（追加した順）"},
        {"G.size() / G.edge_count()", "int", "頂点数 / 辺の本数"},
        {"G.reversed()", "Graph", "全辺の向きを逆にしたグラフ（辺番号は同じ）"},
        {"G.sort_adj()", "void", "隣接を行き先の昇順に並べる"},
        {"G.in_degrees() / out_degrees()", "vector<int>", "入次数 / 出次数"},
        {"G.adj_list() / adj_list_w()", "vector<vector<..>>", "古い形の隣接リストに変換"},
        {"read_undirected(n, m) など", "Graph", "入力を読む（_w は重みあり。base 既定 1）"},
        {"read_tree(n) / read_parents(n)", "Graph", "木を読む（辺 n-1 本 / p_2..p_n）"},
        {"graph_from_adj(adj)", "Graph", "自分で作った隣接リストから"},
        {"graph_from_edges(n, us, vs, UNDIRECTED, base)", "Graph", "辺の配列から（ws 付き・pair 版あり）"},
        {"graph_from_parents(par, base)", "Graph", "親の配列から（w 付き版あり）"},
        {"r.dist[v] / r.par[v]", "Cost / int", "探索の結果：距離 / 直前の頂点"},
        {"r.reached(v) / r.path(t)", "bool / vector<int>", "到達したか / 始点 → t の頂点列"},
        {"r.nearest_source()", "vector<int>", "各頂点に届いた始点"},
        {"parent_edges(G, r) / path_edges(G, r, t)", "vector<int>", "使った辺の番号"},
        {"restore_path(par, t)", "vector<int>", "親の配列から経路"},
    };
    std::cout << "\n■ API 一覧（頂点は 0-indexed。到達不能は DIST_INF。詳しくは下の【1】~）\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, 47) << pad(ret, 23) << r.what << "\n";
    }
}

// デモ用の BFS（本物は 01_bfs の bfs）
static SearchResult demo_bfs(const Graph& G, const std::vector<int>& starts) {
    SearchResult r{std::vector<Cost>(G.size(), DIST_INF), std::vector<int>(G.size(), -1)};
    std::vector<int> q;
    for (int s : starts) r.dist[s] = 0, q.push_back(s);
    for (size_t h = 0; h < q.size(); ++h)
        for (const Edge& e : G[q[h]])
            if (r.dist[e.to] == DIST_INF) r.dist[e.to] = r.dist[q[h]] + 1, r.par[e.to] = q[h], q.push_back(e.to);
    return r;
}

static void demo() {
    banner("グラフ基盤 使い方デモ    見方:  コード  ->  返り値    // 意味");
    cheat_sheet();

    api("【1】 Graph G(n);  G.add_edge(u, v, w);  G.add_undirected(u, v, w);", "辺番号（int）",
        "辺を 1 本ずつ足してグラフを作る。重みを省略すると w = 1。辺番号は足した順に 0, 1, 2, ...", "ならし O(1)");
    RUN(Graph G(4));
    SHOW(G.add_undirected(0, 1, 5), "無向辺 0-1（重み 5）→ 辺番号 0");
    SHOW(G.add_undirected(1, 2), "無向辺 1-2（重み省略 = 1）→ 辺番号 1");
    SHOW(G.add_edge(2, 3, 7), "有向辺 2→3（重み 7）→ 辺番号 2");
    SHOW(G.add(3, 0, 2, DIRECTED), "向きを引数で選ぶ版");
    print_graph("G", G);
    note("無向辺は G[u] と G[v] の両方に、同じ辺番号 id で入る。");

    api("【2】 for (auto& e : G[v]) { e.to; e.w; e.id; }   G.es[id]", "Edge / RawEdge",
        "G[v] は v から出る辺の一覧。G.es[id] は辺 id の {from, to, w}（入力・追加した順）", "O(1)");
    SHOW(G[1], "頂点 1 から出る辺（0 へ重み 5、2 へ重み 1）");
    SHOW(G.es[0], "辺 0 は 0-1（向きは追加した通り）");
    SHOW(G.size(), "頂点数");
    SHOW(G.edge_count(), "辺の本数（無向辺も 1 本と数える）");
    note("辺に追加の情報（時刻・種類など）があるときは、別の配列 info[e.id] を辺番号で引く。");

    api("【3】 read_undirected(N, M) / read_directed / read_undirected_w / read_directed_w", "Graph",
        "「u v」または「u v w」が M 行の入力を読む。1-indexed を 0-indexed に直す（base 引数で変更可）", "O(N + M)");
    {
        CinFromString cf("1 2 10\n2 3 20\n1 3 5\n");
        code("// 入力:  1 2 10 / 2 3 20 / 1 3 5");
        RUN(Graph H = read_undirected_w(3, 3));
        print_graph("H", H);
    }
    {
        CinFromString cf("0 1\n1 2\n");
        code("// 入力が 0-indexed（0 1 / 1 2）なら最後の引数 base に 0");
        RUN(Graph H = read_directed(3, 2, 0));
        print_graph("H", H);
    }

    api("【4】 read_tree(N) / read_tree_w(N) / read_parents(N)", "Graph",
        "木を読む。read_tree は辺 N-1 本、read_parents は p_2 .. p_N（頂点 i の親）の形。どちらも無向辺で入る", "O(N)");
    {
        CinFromString cf("1 1 2 2\n");
        code("// 入力:  1 1 2 2   （頂点 2,3 の親が 1、頂点 4,5 の親が 2）");
        RUN(Graph T = read_parents(5));
        print_graph("T", T);
        note("辺番号は「0-indexed の頂点 i と親を結ぶ辺 = i - 1」。");
    }

    api("【5】 graph_from_adj(adj)", "Graph",
        "自分で作った隣接リスト vector<vector<int>> / vector<vector<pair<int, W>>> から作る。各要素を有向辺 1 本として取り込む",
        "O(N + M)");
    RUN(std::vector<std::vector<int>> adj = {{1, 2}, {0}, {0}});
    RUN(Graph A = graph_from_adj(adj));
    print_graph("A", A);
    RUN(std::vector<std::vector<std::pair<int, long long>>> wadj = {{{1, 3}}, {{0, 3}}});
    RUN(Graph AW = graph_from_adj(wadj));
    print_graph("AW", AW);
    note("両方向に入れた隣接リスト（無向グラフ）なら、そのまま無向グラフとして使える。辺番号は向きごとに別になる。");

    api("【6】 graph_from_edges(n, us, vs, [ws,] UNDIRECTED / DIRECTED, base)   graph_from_edges(n, pairs, 向き, base)",
        "Graph", "A_i, B_i が別々の行で与えられるときなど、辺の配列から作る。base は必ず指定（読んだままなら 1）", "O(N + M)");
    RUN(std::vector<int> As = {1, 2}, Bs = {2, 3});
    RUN(Graph E = graph_from_edges(3, As, Bs, UNDIRECTED, 1));
    print_graph("E", E);
    note("向きは DIRECTED / UNDIRECTED（型が別なので base と順番を間違えるとコンパイルエラーになる）。");

    api("【7】 graph_from_parents(par, [w,] base)", "Graph",
        "親の配列から木（森）を作る。par[i] - base が i の親。負になる値（根の 0 や -1）は「親なし」", "O(N)");
    RUN(std::vector<int> parent = {-1, 0, 0, 1});
    RUN(Graph P = graph_from_parents(parent, 0));
    print_graph("P", P);

    api("【8】 G.reversed() / G.sort_adj() / G.in_degrees() / G.out_degrees() / G.add_vertex()", "",
        "逆辺グラフ（有向グラフで「各頂点から t への距離」を t からの探索で求める）、隣接の整列、次数、頂点の追加",
        "O(N + M) / O(M log M) / O(N + M) / O(1)");
    RUN(Graph D(3));
    RUN(D.add_edge(0, 2); D.add_edge(0, 1); D.add_edge(1, 2));
    RUN(Graph Dr = D.reversed());
    print_graph("Dr", Dr);
    SHOW(D.in_degrees(), "入次数");
    SHOW(D.out_degrees(), "出次数");
    RUN(D.sort_adj());
    SHOW(D[0], "行き先の昇順になった（小さい番号から訪問したいとき）");
    SHOW(D.add_vertex(), "超頂点 3 を追加（例: 特別な頂点すべてと重み 0 でつなぐ）");
    SHOW(D.adj_list(), "vector<vector<int>> に変換（自分で書く探索に渡すとき）");

    api("【9】 SearchResult（BFS・0-1 BFS・ダイクストラ・グリッド探索が返す）", "dist / par",
        "r.dist[v] = 距離（到達不能は DIST_INF）、r.par[v] = 直前の頂点（始点・未到達は -1）", "");
    RUN(Graph Q(6));
    RUN(Q.add_undirected(0, 1); Q.add_undirected(1, 2); Q.add_undirected(3, 4); Q.add_undirected(2, 4));
    RUN(SearchResult r = demo_bfs(Q, {0}));
    note("（ここでは説明用の BFS を使っている。実際は 01_bfs の bfs(G, s) など）");
    SHOW(r.dist, "距離。頂点 5 は到達不能");
    SHOW(r.par, "直前の頂点");
    SHOW(r.reached(5), "到達したか");
    SHOW(r.path(3), "0 → 3 の頂点列");
    SHOW(r.path(5), "到達不能なら空");
    RUN(auto [dist, par] = demo_bfs(Q, {0}));
    SHOW(dist[4], "構造化束縛で分けて受け取ることもできる");

    api("【10】 r.nearest_source()", "vector<int>",
        "多始点の探索で、各頂点に最初に届いた始点（= いちばん近い始点）。未到達は -1", "O(N)");
    RUN(SearchResult r2 = demo_bfs(Q, {0, 3}));
    SHOW(r2.nearest_source(), "頂点 0,1,2 は始点 0、頂点 3,4 は始点 3 が近い");

    api("【11】 parent_edges(G, r) / path_edges(G, r, t)", "vector<int>",
        "各頂点に入るときに使った辺の番号 / 始点 → t の経路の辺番号の列。最短路木の辺を出力する問題で使う（例: ABC252 E）",
        "O(N + M)");
    SHOW(parent_edges(Q, r), "頂点 v に入った辺の番号（始点・未到達は -1）");
    SHOW(path_edges(Q, r, 3), "0 → 3 で通る辺の番号");

    api("【12】 restore_path(par, t)", "vector<int>", "自分で作った親の配列から経路を復元（par[始点] = -1）", "O(経路の長さ)");
    SHOW(restore_path(std::vector<int>{-1, 0, 1, 1}, 3), "0 → 1 → 3");

    head("典型的な使い方");
    note("無向グラフの最短路    : Graph G = read_undirected(N, M);  auto r = bfs(G, 0);  r.dist[N - 1]");
    note("各頂点から t へ（有向）: auto r = dijkstra(G.reversed(), t);");
    note("多始点                : bfs(G, {s1, s2, ...})  →  r.nearest_source() でどの始点が近いか");
    note("辺の削除・特殊な移動  : bfs_core(N, starts, [&](int v, auto&& push) { for (auto& e : G[v]) if (e.id != ban) push(e.to); });");
    note("類題: ABC252 E（最短路木の辺）、ABC211 D（最短路の本数）、ABC218 F（辺を 1 本ずつ消す）");
}

}  // namespace gb_demo


// =============================================================================
//  Windows のコンソールで日本語が文字化けしないようにする（Windows 以外では何もしない）
// =============================================================================
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
    gb_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    gb_test::benchmark();    // 速度計測
    gb_demo::demo();         // 使い方デモ
    return 0;
}
