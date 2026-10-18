// =============================================================================
//  BFS（幅優先探索）  bfs_core / bfs                              --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 01 探索
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed。結果は SearchResult（r.dist[v] / r.par[v] / r.path(t)。到達不能は DIST_INF）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  考え方
//    本体 bfs_core は「頂点数 + 隣接を列挙する処理（ラムダ）」を受け取る。グラフ・グリッド・状態グラフ・
//    辺を 1 本消したグラフ・数の操作など、隣接の書き方を変えるだけで同じ本体を使う。
//      bfs_core(n, s / starts, [&](int v, auto&& push) { ... v から 1 歩で行ける u ごとに push(u) ... });
//
//  高速化のポイント
//    - キューは配列 1 本（std::queue を使わない）。到達済みかどうかは dist だけで判定
//    - ラムダはテンプレートで受けるのでインライン展開される（関数ポインタ・std::function を使わない）
//
//  機能一覧                                                                    計算量
//    bfs(G, s) / bfs(G, {s1, s2, ...})   辺の本数での最短距離（多始点も可）     O(N + M)
//    bfs(G, s, target)                   target の距離が決まったら打ち切る       O(N + M)
//    bfs_core(n, s / starts, adj, target = -1)   隣接を自分で書く版                  O(N + 隣接の数)
//    （結果）r.dist / r.par / r.reached(v) / r.path(t) / r.nearest_source() / parent_edges(G, r)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    bfs_test::random_test()  … 愚直解と突き合わせるランダムテスト
//    bfs_test::benchmark()    … 速度計測
//    bfs_demo::demo()         … 使い方デモ（隣接を自分で書く例つき）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
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
#ifndef GG_BFS
#define GG_BFS
// ---- BFS（グラフ・グリッド第二世代）  本体 bfs_core は「頂点数 + 隣接を列挙する処理」を受け取る ----
// adj(v, push): v から 1 歩で行ける頂点 u ごとに push(u) を呼ぶ（到達済みかどうかは本体が調べる）。
// 通常のグラフ・グリッド・状態グラフ・辺を消したグラフなど、隣接の書き方を変えるだけで同じ本体を使える。
// target を指定すると、target の距離が決まった時点で打ち切る（ほかの頂点の距離は途中まで）。
template <class F>
SearchResult bfs_core(int n, const std::vector<int>& starts, F&& adj, int target = -1) {
    SearchResult r{std::vector<Cost>(n, DIST_INF), std::vector<int>(n, -1)};
    std::vector<int> q(n);  // 各頂点は高々 1 回しか入らないので n で足りる（生のポインタで回すほうが速い）
    int qt = 0;
    Cost* const dist = r.dist.data();
    int* const par = r.par.data();
    int* const qp = q.data();
    for (int s : starts)
        if (dist[s] != 0) dist[s] = 0, qp[qt++] = s;
    for (int h = 0; h < qt; ++h) {
        const int v = qp[h];
        if (v == target) break;
        const Cost nd = dist[v] + 1;
        adj(v, [&](int u) {
            if (dist[u] == DIST_INF) dist[u] = nd, par[u] = v, qp[qt++] = u;
        });
    }
    return r;
}
template <class F>
SearchResult bfs_core(int n, int s, F&& adj, int target = -1) { return bfs_core(n, std::vector<int>{s}, adj, target); }
// グラフの BFS（辺の重みは見ない。辺の本数での最短距離）
inline SearchResult bfs(const Graph& G, const std::vector<int>& starts, int target = -1) {
    return bfs_core(G.size(), starts, [&](int v, auto&& push) {
        for (const Edge& e : G[v]) push(e.to);
    }, target);
}
inline SearchResult bfs(const Graph& G, int s, int target = -1) { return bfs(G, std::vector<int>{s}, target); }
#endif  // GG_BFS
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "bfs";
// ---- テスト用の小道具 ----
#define CHECK(...)                                                                                \
    do {                                                                                          \
        if (!(__VA_ARGS__)) {                                                                     \
            std::cerr << "[" << TEST_TAG << "] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                         \
        }                                                                                         \
    } while (0)
// テスト量の倍率。デバッグビルド（最適化なし・範囲チェックあり）は遅いので自動で 1/8 にする
#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif
[[maybe_unused]] static int rnd(std::mt19937_64& rng, int lo, int hi) { return lo + (int)(rng() % (unsigned long long)(hi - lo + 1)); }
[[maybe_unused]] static long long ms_since(std::chrono::steady_clock::time_point t0) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
struct CinFromString {  // std::cin を文字列に差し替える
    std::istringstream ss;
    std::streambuf* old;
    explicit CinFromString(const std::string& s) : ss(s), old(std::cin.rdbuf(ss.rdbuf())) {}
    ~CinFromString() { std::cin.rdbuf(old); }
};

// =============================================================================
//  ランダムテスト
// =============================================================================
namespace bfs_test {

// 愚直解: 全部の辺を N 回緩和する（単位長）
static std::vector<Cost> naive(int n, const std::vector<std::vector<int>>& adj, const std::vector<int>& starts) {
    std::vector<Cost> d(n, DIST_INF);
    for (int s : starts) d[s] = 0;
    for (int it = 0; it < n; ++it)
        for (int v = 0; v < n; ++v)
            if (d[v] != DIST_INF)
                for (int u : adj[v]) d[u] = std::min(d[u], d[v] + 1);
    return d;
}
static void check_tree(const SearchResult& r, const std::vector<std::vector<int>>& adj, const std::vector<int>& starts) {
    const int n = (int)adj.size();
    std::vector<char> st(n, 0);
    for (int s : starts) st[s] = 1;
    for (int v = 0; v < n; ++v) {
        if (!r.reached(v)) {
            CHECK(r.par[v] == -1 && r.path(v).empty());
            continue;
        }
        if (r.par[v] == -1) {
            CHECK(st[v] && r.dist[v] == 0);
        } else {
            const int p = r.par[v];
            CHECK(r.dist[p] + 1 == r.dist[v]);
            CHECK(std::find(adj[p].begin(), adj[p].end(), v) != adj[p].end());
        }
        const std::vector<int> path = r.path(v);
        CHECK((Cost)path.size() == r.dist[v] + 1 && path.back() == v && st[path[0]]);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 20), m = rnd(rng, 0, 3 * n);
        const bool und = rng() & 1;
        Graph G(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            und ? G.add_undirected(a, b, rnd(rng, 0, 9)) : G.add_edge(a, b, rnd(rng, 0, 9));
        }
        const std::vector<std::vector<int>> adj = G.adj_list();
        std::vector<int> starts;
        const int k = rnd(rng, 0, 3);
        for (int i = 0; i < k; ++i) starts.push_back(rnd(rng, 0, n - 1));
        const std::vector<Cost> ex = naive(n, adj, starts);
        // 1. bfs(G, starts) / bfs(G, s)
        const SearchResult r = bfs(G, starts);
        CHECK(r.dist == ex);
        check_tree(r, adj, starts);
        if (k == 1) CHECK(bfs(G, starts[0]).dist == ex);
        // 2. target で打ち切っても target の距離と経路は正しい
        const int t = rnd(rng, 0, n - 1);
        const SearchResult rt = bfs(G, starts, t);
        CHECK(rt.dist[t] == ex[t]);
        if (rt.reached(t)) CHECK((Cost)rt.path(t).size() == ex[t] + 1);
        for (int v = 0; v < n; ++v) CHECK(rt.dist[v] == DIST_INF || rt.dist[v] == ex[v]);
        // 3. bfs_core に隣接を自分で書く: 辺 ban を使わない
        const int ban = m ? rnd(rng, 0, m - 1) : -1;
        const SearchResult rb = bfs_core(n, starts, [&](int v, auto&& push) {
            for (const Edge& e : G[v])
                if (e.id != ban) push(e.to);
        });
        std::vector<std::vector<int>> adj2(n);
        for (int v = 0; v < n; ++v)
            for (const Edge& e : G[v])
                if (e.id != ban) adj2[v].push_back(e.to);
        CHECK(rb.dist == naive(n, adj2, starts));
        check_tree(rb, adj2, starts);
        // 4. 状態グラフ（頂点, 歩数の偶奇）: (v, p) → (u, 1 - p)
        if (!starts.empty()) {
            const SearchResult rp = bfs_core(2 * n, starts[0] * 2, [&](int x, auto&& push) {
                for (int u : adj[x / 2]) push(u * 2 + (1 - x % 2));
            });
            std::vector<std::vector<int>> adj3(2 * n);
            for (int x = 0; x < 2 * n; ++x)
                for (int u : adj[x / 2]) adj3[x].push_back(u * 2 + (1 - x % 2));
            CHECK(rp.dist == naive(2 * n, adj3, {starts[0] * 2}));
        }
    }
    // 端のケース: 始点が重複、始点なし
    {
        Graph G(3);
        G.add_edge(0, 1);
        const SearchResult r = bfs(G, {0, 0, 1});
        CHECK(r.dist == (std::vector<Cost>{0, 0, DIST_INF}) && r.par == (std::vector<int>{-1, -1, -1}));
        const SearchResult e = bfs(G, std::vector<int>{});
        CHECK(e.dist == std::vector<Cost>(3, DIST_INF));
    }
    std::cout << "[bfs] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[bfs] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> us(n), vs(n);
        for (int i = 0; i < n; ++i) us[i] = (int)(rng() % n), vs[i] = (int)(rng() % n);
        const Graph G = graph_from_edges(n, us, vs, UNDIRECTED, 0);
        auto t0 = std::chrono::steady_clock::now();
        const SearchResult r = bfs(G, 0);
        const long long t1 = ms_since(t0);
        // 暗黙のグラフ（数 x → 2x mod n, x + 1 mod n）
        t0 = std::chrono::steady_clock::now();
        const SearchResult r2 = bfs_core(n, 1, [&](int x, auto&& push) {
            push((int)(2LL * x % n)), push((x + 1) % n);
        });
        const long long t2 = ms_since(t0);
        std::cout << "[bfs] benchmark: N = M = " << n << "  bfs(G, 0) " << t1 << " ms, bfs_core（数の操作のグラフ）" << t2 << " ms  (check "
                  << r.dist[n - 1] + r2.dist[n - 1] << ")\n";
    }
}

}  // namespace bfs_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace bfs_demo {
// ---- デモ用の小道具（見方: コード -> 返り値 // 意味） ----
template <class A, class B>
[[maybe_unused]] static std::string to_s(const std::pair<A, B>& p);
template <class T, size_t N>
[[maybe_unused]] static std::string to_s(const std::array<T, N>& a);
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v);
#ifdef GG_GRAPH_BASE
[[maybe_unused]] static std::string to_s(const Edge& e) { return "{to=" + std::to_string(e.to) + ", id=" + std::to_string(e.id) + ", w=" + std::to_string(e.w) + "}"; }
[[maybe_unused]] static std::string to_s(const RawEdge& e) { return "{" + std::to_string(e.from) + "-" + std::to_string(e.to) + ", w=" + std::to_string(e.w) + "}"; }
#endif
[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(char c) { return std::string("'") + c + "'"; }
[[maybe_unused]] static std::string to_s(const std::string& s) { return "\"" + s + "\""; }
[[maybe_unused]] static std::string to_s(const char* s) { return "\"" + std::string(s) + "\""; }
template <class T, std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool> && !std::is_same_v<T, char>, int> = 0>
[[maybe_unused]] static std::string to_s(T v) {
    if (std::is_same_v<T, long long> && (long long)v == (1LL << 60)) return "DIST_INF";
    return std::to_string(v);
}
template <class A, class B>
[[maybe_unused]] static std::string to_s(const std::pair<A, B>& p) { return "(" + to_s(p.first) + ", " + to_s(p.second) + ")"; }
template <class T, size_t N>
[[maybe_unused]] static std::string to_s(const std::array<T, N>& a) {
    std::string s = "(";
    for (size_t i = 0; i < N; ++i) s += (i ? ", " : "") + to_s(a[i]);
    return s + ")";
}
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
template <class T, class = decltype(std::declval<const T&>().to_vector())>
[[maybe_unused]] static std::string to_s(const T& v) { return to_s(v.to_vector()); }  // IntSpan など
[[maybe_unused]] static std::string pad(const std::string& s, size_t w) {  // 全角を 2 文字分として幅をそろえる
    size_t cols = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = (unsigned char)s[i];
        if ((c & 0xC0) == 0x80) continue;
        cols += (c >= 0xE0 ? 2 : 1);
    }
    return cols >= w ? s : s + std::string(w - cols, ' ');
}
[[maybe_unused]] static void banner(const std::string& s) {
    std::cout << "\n=====================================================================\n"
              << "  " << s << "\n"
              << "=====================================================================\n";
}
[[maybe_unused]] static void api(const std::string& usage, const std::string& ret, const std::string& what, const std::string& cost) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << usage << "\n"
              << "    返り値 : " << ret << "\n"
              << "    内容   : " << what << "\n"
              << "    計算量 : " << cost << "\n";
}
[[maybe_unused]] static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
[[maybe_unused]] static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
[[maybe_unused]] static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
[[maybe_unused]] static void stmt(const std::string& s) {
    const bool block = !s.empty() && s.back() == '}' &&
                       (s.compare(0, 4, "for ") == 0 || s.compare(0, 3, "if ") == 0 || s.compare(0, 6, "while ") == 0);
    code(block ? s : s + ";");
}
template <class V>
[[maybe_unused]] static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 36) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 56) + "  // " + comment;
    std::cout << line << "\n";
}
struct Row { const char* use; const char* ret; const char* what; };
template <size_t N>
[[maybe_unused]] static void cheat_sheet(const std::string& title, const Row (&rows)[N], size_t w1 = 46, size_t w2 = 22) {
    std::cout << "\n■ " << title << "\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, w1) << pad(ret, w2) << r.what << "\n";
    }
}
#define SHOW(expr, comment) show(#expr, (expr), comment)
#define RUN(...)             \
    stmt(#__VA_ARGS__);      \
    __VA_ARGS__

static void demo() {
    banner("BFS 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto r = bfs(G, s);", "SearchResult", "辺の本数での最短距離（重みは見ない）"},
        {"auto r = bfs(G, {s1, s2});", "SearchResult", "多始点（いちばん近い始点からの距離）"},
        {"auto r = bfs(G, s, t);", "SearchResult", "t の距離が決まったら打ち切る（t 以外は途中まで）"},
        {"bfs_core(n, s / starts, adj, t = -1)", "SearchResult", "隣接を自分で書く: adj(v, push) で push(u)"},
        {"r.dist[v] / r.par[v] / r.path(t)", "", "距離（届かなければ DIST_INF）/ 直前の頂点 / 経路"},
        {"auto [dist, par] = bfs(G, s);", "", "構造化束縛で受け取る"},
    };
    cheat_sheet("API 一覧", rows);

    head("デモで使うグラフ（無向）: 0-1, 1-2, 2-3, 0-4, 4-3, 5-6");
    RUN(Graph G(7));
    RUN(G.add_undirected(0, 1); G.add_undirected(1, 2); G.add_undirected(2, 3); G.add_undirected(0, 4); G.add_undirected(4, 3); G.add_undirected(5, 6));

    api("【1】 bfs(G, s)", "SearchResult", "s からの最短距離（辺の本数）。到達不能は DIST_INF。r.path(t) で経路", "O(N + M)");
    RUN(auto r = bfs(G, 0));
    SHOW(r.dist, "頂点 5, 6 には届かない");
    SHOW(r.par, "直前の頂点");
    SHOW(r.path(3), "0 → 3 の最短経路（0 → 4 → 3）");
    SHOW(r.reached(5), "届いたか");

    api("【2】 bfs(G, {s1, s2, ...})", "SearchResult", "多始点。各頂点の「いちばん近い始点までの距離」。r.nearest_source() でどの始点か", "O(N + M)");
    RUN(auto r2 = bfs(G, {1, 3}));
    SHOW(r2.dist, "始点 1 か 3 のどちらか近いほうからの距離");
    SHOW(r2.nearest_source(), "いちばん近い始点");

    api("【3】 bfs(G, s, t)", "SearchResult", "t の距離が決まった時点で打ち切る。t 以外の頂点の距離は途中まで（使わないこと）", "O(N + M)");
    SHOW(bfs(G, 0, 1).dist[1], "1 までの距離だけ欲しいとき");

    api("【4】 bfs_core(n, s / starts, adj)", "SearchResult",
        "隣接を自分で書く。adj(v, push) の中で、v から 1 歩で行ける u ごとに push(u) を呼ぶ（到達済みかどうかは本体が調べる）", "O(n + 隣接の数)");
    code("// 例 1: 辺 1 番（1-2）を使わない BFS（辺を 1 本ずつ消して調べる問題: ABC218 F）");
    RUN(int ban = 1);
    RUN(auto rb = bfs_core(G.size(), 0, [&](int v, auto&& push) { for (auto& e : G[v]) if (e.id != ban) push(e.to); }));
    SHOW(rb.dist, "辺 1（1-2）を使わずに");
    code("// 例 2: 状態つき（頂点, 歩数の偶奇）→ 番号 v * 2 + p。00_state_id の StateId を使うとわかりやすい");
    RUN(auto rp = bfs_core(G.size() * 2, 0 * 2 + 0, [&](int x, auto&& push) { for (auto& e : G[x / 2]) push(e.to * 2 + (1 - x % 2)); }));
    SHOW(rp.dist[3 * 2 + 1], "頂点 3 に奇数歩で着く最短（3 歩）");
    SHOW(rp.dist[3 * 2 + 0], "頂点 3 に偶数歩で着く最短（0-1-2-3 は 3 歩、0-4-3 は 2 歩）");
    code("// 例 3: 数の操作（x → x * 2 mod 10、x → x + 1 mod 10）で 1 から 7 にする最小回数");
    RUN(auto rn = bfs_core(10, 1, [&](int x, auto&& push) { push(x * 2 % 10); push((x + 1) % 10); }));
    SHOW(rn.dist[7], "1 → 2 → 3 → 6 → 7 など");
    SHOW(rn.path(7), "その手順");
    note("グリッドは 09_grid_bfs の grid_bfs（中身は bfs_core + g.for_each4）。");

    api("【5】 auto [dist, par] = bfs(G, s);", "", "SearchResult は dist と par の 2 つだけなので、構造化束縛で分けて受け取れる", "");
    RUN(auto [dist, par] = bfs(G, 5));
    SHOW(dist[6], "5 → 6");

    head("典型的な使い方");
    note("迷路・グラフの最短手数   : bfs(G, s).dist[t]（到達不能なら DIST_INF → -1 を出力など）");
    note("各頂点から最も近い○○     : bfs(G, 特別な頂点の一覧)  →  r.dist[v]、r.nearest_source()[v]");
    note("辺を消して調べる         : bfs_core で e.id を飛ばす");
    note("状態つき（鍵・偶奇など） : 00_state_id の StateId と bfs_core");
    note("最短路の本数             : 02_shortest_path_dag の count_shortest_paths(G, r, mod)");
}

}  // namespace bfs_demo


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
    bfs_test::random_test();
    bfs_test::benchmark();
    bfs_demo::demo();
    return 0;
}
