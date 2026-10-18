// =============================================================================
//  最短路 DAG（最短路の本数・最短路に乗る頂点と辺・必ず通る頂点と辺）   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 最短路・状態グラフ
//  依存       : 00_graph_base, 02_dijkstra（同梱済み）
//  頂点番号   : 0-indexed。重みは 0 以上（重み 0 の辺も可）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  「最短路に乗る辺」= dist[u] + w == dist[v] となる辺 u → v。これだけを集めると DAG になる（重み 0 の閉路がなければ）。
//  s → t の最短路に乗る ⇔ ds[u] + w + dt[v] == ds[t]（ds = s からの距離、dt = t までの距離）。
//  無向辺は向きごとに調べる（G の隣接 G[u] をたどるので、有向・無向を区別しなくてよい）。
//
//  機能一覧                                                                       計算量
//    count_shortest_paths(G, r, mod)  始点から各頂点への最短路の本数（mod。mod = 0 なら割らない）。
//                                     r は G の BFS / 0-1 BFS / ダイクストラの結果（重みと距離が合っていること）。
//                                     重み 0 の閉路の先（無限にある）は -1                    O(N + M)
//    order_by_dist(r)                 到達した頂点を距離の昇順に                              O(N log N)
//    ShortestPathDag P(G, s, t, mod = 998244353)   s → t についてまとめて求める              O((N + M) log M)
//      P.best / P.reachable() / P.count()          最短距離 / 届くか / 最短路の本数（mod）
//      P.on_vertex(v) / P.on_edge(id)              どれかの最短路に乗るか
//      P.must_vertex(v) / P.must_edge(id)          すべての最短路が通るか（消すと最短距離が変わる）
//      P.count_through_vertex(v) / P.count_through_edge(id)   そこを通る最短路の本数（mod）
//      P.path() / P.ds / P.dt / P.cnt_s / P.cnt_t   最短路を 1 本 / s から・t までの距離 / 本数
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤とダイクストラも入っている）
//    spd_test::random_test() / spd_test::benchmark() / spd_demo::demo()
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
#ifndef GG_DIJKSTRA
#define GG_DIJKSTRA
// ---- ダイクストラ（グラフ・グリッド第二世代）  重みが 0 以上の最短路。本体 dijkstra_core ----
// adj(v, push): v から行ける頂点 u ごとに push(u, w) を呼ぶ（w >= 0）。target は bfs_core と同じ。
// starts に {頂点, 初期コスト} の組を渡すと、始点ごとに初期コストを付けられる（「各頂点の値 + 距離」の最小など）
namespace gg_internal {
// (距離, 頂点) の最小ヒープ。std::priority_queue と同じ働きだが、-D_GLIBCXX_DEBUG のデバッグビルドでも
// 1 回の操作が O(log n) のまま（priority_queue はデバッグビルドだと毎回ヒープ全体を検査して O(n) になる）
struct DistHeap {
    std::vector<std::pair<Cost, int>> h;
    bool empty() const { return h.empty(); }
    void push(Cost d, int v) {
        size_t i = h.size();
        h.emplace_back();
        while (i > 0) {
            const size_t p = (i - 1) >> 1;
            if (h[p].first <= d) break;
            h[i] = h[p], i = p;
        }
        h[i] = {d, v};
    }
    std::pair<Cost, int> pop() {
        const std::pair<Cost, int> top = h[0], x = h.back();
        h.pop_back();
        const size_t n = h.size();
        if (n > 0) {
            size_t i = 0;
            for (size_t c = 1; c < n; c = 2 * i + 1) {
                if (c + 1 < n && h[c + 1].first < h[c].first) ++c;
                if (x.first <= h[c].first) break;
                h[i] = h[c], i = c;
            }
            h[i] = x;
        }
        return top;
    }
};
}  // namespace gg_internal
template <class F>
SearchResult dijkstra_core(int n, const std::vector<std::pair<int, Cost>>& starts, F&& adj, int target = -1) {
    SearchResult r{std::vector<Cost>(n, DIST_INF), std::vector<int>(n, -1)};
    Cost* const dist = r.dist.data();
    int* const par = r.par.data();
    gg_internal::DistHeap pq;
    for (const auto& s : starts)
        if (s.second < dist[s.first]) dist[s.first] = s.second, pq.push(s.second, s.first);
    while (!pq.empty()) {
        const std::pair<Cost, int> top = pq.pop();
        const Cost d = top.first;
        const int v = top.second;
        if (d != dist[v]) continue;
        if (v == target) break;
        adj(v, [&](int u, Cost w) {
            if (d + w < dist[u]) dist[u] = d + w, par[u] = v, pq.push(d + w, u);
        });
    }
    return r;
}
template <class F>
SearchResult dijkstra_core(int n, const std::vector<int>& starts, F&& adj, int target = -1) {
    std::vector<std::pair<int, Cost>> st;
    st.reserve(starts.size());
    for (int s : starts) st.emplace_back(s, 0);
    return dijkstra_core(n, st, adj, target);
}
template <class F>
SearchResult dijkstra_core(int n, int s, F&& adj, int target = -1) {  // 始点 1 つ（{s} と書いても曖昧にならないように）
    return dijkstra_core(n, std::vector<std::pair<int, Cost>>{{s, 0}}, adj, target);
}
// グラフのダイクストラ（辺の重み e.w >= 0）
inline SearchResult dijkstra(const Graph& G, const std::vector<std::pair<int, Cost>>& starts, int target = -1) {
    return dijkstra_core(G.size(), starts, [&](int v, auto&& push) {
        for (const Edge& e : G[v]) push(e.to, e.w);
    }, target);
}
inline SearchResult dijkstra(const Graph& G, const std::vector<int>& starts, int target = -1) {
    return dijkstra_core(G.size(), starts, [&](int v, auto&& push) {
        for (const Edge& e : G[v]) push(e.to, e.w);
    }, target);
}
inline SearchResult dijkstra(const Graph& G, int s, int target = -1) { return dijkstra(G, std::vector<int>{s}, target); }
// 密なグラフ用（辺が N^2 本近いとき）。cost(u, v) = u → v の重み（辺がなければ DIST_INF）。O(N^2)
template <class W>
SearchResult dijkstra_dense(int n, const std::vector<int>& starts, W&& cost) {
    SearchResult r{std::vector<Cost>(n, DIST_INF), std::vector<int>(n, -1)};
    std::vector<char> done(n, 0);
    for (int s : starts) r.dist[s] = 0;
    for (int it = 0; it < n; ++it) {
        int v = -1;
        for (int u = 0; u < n; ++u)
            if (!done[u] && r.dist[u] != DIST_INF && (v == -1 || r.dist[u] < r.dist[v])) v = u;
        if (v == -1) break;
        done[v] = 1;
        for (int u = 0; u < n; ++u) {
            if (done[u]) continue;
            const Cost w = cost(v, u);
            if (w != DIST_INF && r.dist[v] + w < r.dist[u]) r.dist[u] = r.dist[v] + w, r.par[u] = v;
        }
    }
    return r;
}
template <class W>
SearchResult dijkstra_dense(int n, int s, W&& cost) { return dijkstra_dense(n, std::vector<int>{s}, cost); }
#endif  // GG_DIJKSTRA
#ifndef GG_SHORTEST_PATH_DAG
#define GG_SHORTEST_PATH_DAG
// ---- 最短路 DAG（グラフ・グリッド第二世代）  最短路の本数・最短路に乗る頂点と辺・必ず通る頂点と辺 ----
// 「最短路に乗る辺」= dist[u] + w == dist[v] となる辺 u → v（無向辺は向きごとに調べる）。重みは 0 以上
// 到達した頂点を距離の昇順に（同じ距離なら番号順）
inline std::vector<int> order_by_dist(const SearchResult& r) {
    std::vector<int> v;
    for (int i = 0; i < (int)r.dist.size(); ++i)
        if (r.reached(i)) v.push_back(i);
    std::stable_sort(v.begin(), v.end(), [&](int a, int b) { return r.dist[a] < r.dist[b]; });
    return v;
}
namespace gg_internal {
// 最短路の本数を k 個（3 個まで）の mod でまとめて数える（Kahn 法 1 回）。数えられない頂点は -1
inline void count_paths_multi(const Graph& G, const SearchResult& r, const long long* mods, int k, std::vector<long long>* out) {
    const int n = G.size();
    std::vector<int> indeg(n, 0), q;
    for (int u = 0; u < n; ++u)
        if (r.reached(u))
            for (const Edge& e : G[u])
                if (e.to != u && r.dist[u] + e.w == r.dist[e.to]) ++indeg[e.to];
    std::vector<std::array<long long, 3>> c(n, {0, 0, 0});  // 3 つの mod の値を並べて持つ（メモリを 1 回で読む）
    q.reserve(n);
    for (int v = 0; v < n; ++v) {
        if (!r.reached(v)) continue;
        if (r.dist[v] == 0 && r.par[v] == -1) c[v] = {1, 1, 1};  // 始点
        if (indeg[v] == 0) q.push_back(v);
    }
    std::vector<char> done(n, 0);
    for (size_t h = 0; h < q.size(); ++h) {
        const int u = q[h];
        done[u] = 1;
        const std::array<long long, 3> cu = c[u];
        for (const Edge& e : G[u])
            if (e.to != u && r.dist[u] + e.w == r.dist[e.to]) {
                std::array<long long, 3>& cv = c[e.to];
                for (int j = 0; j < k; ++j) {
                    cv[j] += cu[j];
                    if (mods[j] && cv[j] >= mods[j]) cv[j] -= mods[j];  // どちらも mod 未満なので引けばよい
                }
                if (--indeg[e.to] == 0) q.push_back(e.to);
            }
    }
    for (int j = 0; j < k; ++j) {
        out[j].resize(n);
        for (int v = 0; v < n; ++v) out[j][v] = (r.reached(v) && !done[v]) ? -1 : c[v][j];
    }
}
}  // namespace gg_internal
// 始点（r.dist == 0 かつ r.par == -1 の頂点）から各頂点への最短路の本数（mod。mod = 0 なら割らない）。
// r は G の BFS / 0-1 BFS / ダイクストラの結果（多始点なら始点ごとに 1 本と数える）。未到達は 0。
// 重み 0 の閉路の先など数えられない（無限にある）頂点は -1
inline std::vector<long long> count_shortest_paths(const Graph& G, const SearchResult& r, long long mod) {
    std::vector<long long> c;
    gg_internal::count_paths_multi(G, r, &mod, 1, &c);
    return c;
}

// s → t の最短路の DAG をまとめて調べる（s から G で、t から G.reversed() でダイクストラ）
struct ShortestPathDag {
    int s = -1, t = -1;
    long long mod = 0;
    Cost best = DIST_INF;            // s → t の最短距離（届かなければ DIST_INF）
    SearchResult ds, dt;             // ds: s からの距離、dt: t までの距離（G.reversed() で t から）
    std::vector<long long> cnt_s;    // s → v の最短路の本数（mod）
    std::vector<long long> cnt_t;    // v → t の最短路の本数（mod）
    std::vector<int> on_v, on_e;     // 頂点・辺が s → t のどれかの最短路に乗るか（1 / 0）
    std::vector<long long> cnt_e;    // 辺を通る s → t の最短路の本数（mod）
    std::vector<int> must_v, must_e; // 頂点・辺が s → t のすべての最短路に乗るか（1 / 0。2 つの素数での本数で判定）

    ShortestPathDag() {}
    ShortestPathDag(const Graph& G, int s_, int t_, long long mod_ = 998244353) : s(s_), t(t_), mod(mod_) {
        const int n = G.size();
        const Graph R = G.reversed();
        ds = dijkstra(G, s), dt = dijkstra(R, t), best = ds.dist[t];
        const long long P1 = 998244353, P2 = 1000000007, mods[3] = {mod, P1, P2};
        std::vector<long long> cs[3], ct[3];  // [0] = mod, [1] = P1, [2] = P2（2 つの素数は「必ず通る」の判定用）
        gg_internal::count_paths_multi(G, ds, mods, 3, cs);
        gg_internal::count_paths_multi(R, dt, mods, 3, ct);
        cnt_s = std::move(cs[0]), cnt_t = std::move(ct[0]);
        const std::vector<long long>&a1 = cs[1], &b1 = ct[1], &a2 = cs[2], &b2 = ct[2];
        on_v.assign(n, 0), must_v.assign(n, 0), on_e.assign(G.edge_count(), 0), cnt_e.assign(G.edge_count(), 0), must_e.assign(G.edge_count(), 0);
        if (best == DIST_INF) return;
        auto mul = [](long long a, long long b, long long m) { return m ? (a % m) * (b % m) % m : a * b; };
        std::vector<long long> e1(G.edge_count(), 0), e2(G.edge_count(), 0);
        for (int v = 0; v < n; ++v) {
            if (!ds.reached(v) || !dt.reached(v) || ds.dist[v] + dt.dist[v] != best) continue;
            on_v[v] = 1;
            if (a1[v] >= 0 && b1[v] >= 0 && a1[t] >= 0)
                must_v[v] = (mul(a1[v], b1[v], P1) == a1[t] && mul(a2[v], b2[v], P2) == a2[t]);
            for (const Edge& e : G[v]) {
                const int u = e.to;
                if (u == v || !dt.reached(u) || ds.dist[v] + e.w + dt.dist[u] != best) continue;
                on_e[e.id] = 1;
                if (cnt_s[v] >= 0 && cnt_t[u] >= 0) {
                    cnt_e[e.id] += mul(cnt_s[v], cnt_t[u], mod);
                    if (mod) cnt_e[e.id] %= mod;
                }
                if (a1[v] >= 0 && b1[u] >= 0) e1[e.id] = (e1[e.id] + mul(a1[v], b1[u], P1)) % P1, e2[e.id] = (e2[e.id] + mul(a2[v], b2[u], P2)) % P2;
            }
        }
        if (a1[t] >= 0)
            for (int i = 0; i < G.edge_count(); ++i) must_e[i] = on_e[i] && e1[i] == a1[t] && e2[i] == a2[t];
    }
    bool reachable() const { return best != DIST_INF; }
    long long count() const { return reachable() ? cnt_s[t] : 0; }  // s → t の最短路の本数（mod。無限なら -1）
    bool on_vertex(int v) const { return on_v[v]; }
    bool on_edge(int id) const { return on_e[id]; }
    bool must_vertex(int v) const { return must_v[v]; }   // すべての最短路が v を通る（s と t も true）
    bool must_edge(int id) const { return must_e[id]; }   // すべての最短路が辺 id を通る（消すと最短距離が変わる）
    long long count_through_vertex(int v) const {        // v を通る最短路の本数（mod）
        if (!on_v[v] || cnt_s[v] < 0 || cnt_t[v] < 0) return on_v[v] ? -1 : 0;
        return mod ? cnt_s[v] % mod * (cnt_t[v] % mod) % mod : cnt_s[v] * cnt_t[v];
    }
    long long count_through_edge(int id) const { return cnt_e[id]; }  // 辺 id を通る最短路の本数（mod）
    std::vector<int> path() const { return ds.path(t); }               // 最短路を 1 本（頂点列）
};
#endif  // GG_SHORTEST_PATH_DAG
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "shortest_path_dag";
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

namespace spd_test {

// 愚直解: s から出る全部の経路（辺の列）を DFS で列挙し、t までの最短のものを集める（重みは正なので最短路は単純）
struct Brute {
    Cost best = DIST_INF;
    std::vector<std::vector<std::pair<int, int>>> paths;  // 最短路ごとの {頂点, 辺番号} の列
    Brute(const Graph& G, int s, int t) {
        std::vector<std::pair<int, int>> cur{{s, -1}};
        std::vector<char> used(G.size(), 0);
        used[s] = 1;
        std::function<void(int, Cost)> rec = [&](int v, Cost d) {
            if (v == t) {
                if (d < best) best = d, paths.clear();
                if (d == best) paths.push_back(cur);
                return;
            }
            for (const Edge& e : G[v]) {
                if (used[e.to]) continue;
                used[e.to] = 1, cur.push_back({e.to, e.id});
                rec(e.to, d + e.w);
                used[e.to] = 0, cur.pop_back();
            }
        };
        rec(s, 0);
    }
};

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    const long long MOD = 1000000007;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 7), m = rnd(rng, 0, 12);
        const bool und = rng() & 1;
        Graph G(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            const Cost w = rnd(rng, 1, 3);
            und ? G.add_undirected(a, b, w) : G.add_edge(a, b, w);
        }
        const int s = rnd(rng, 0, n - 1), t = rnd(rng, 0, n - 1);
        const Brute B(G, s, t);
        const ShortestPathDag P(G, s, t, MOD);
        CHECK(P.best == B.best && P.reachable() == (B.best != DIST_INF));
        CHECK(P.count() == (long long)B.paths.size());
        std::vector<int> on_v(n, 0), on_e(G.edge_count(), 0), must_v(n, 0), must_e(G.edge_count(), 0);
        std::vector<long long> thr_v(n, 0), thr_e(G.edge_count(), 0);
        if (!B.paths.empty()) {
            std::fill(must_v.begin(), must_v.end(), 1), std::fill(must_e.begin(), must_e.end(), 1);
            for (const auto& p : B.paths) {
                std::vector<char> pv(n, 0), pe(G.edge_count(), 0);
                for (auto [v, id] : p) {
                    pv[v] = 1, on_v[v] = 1, ++thr_v[v];
                    if (id >= 0) pe[id] = 1, on_e[id] = 1, ++thr_e[id];
                }
                for (int v = 0; v < n; ++v) must_v[v] &= pv[v];
                for (int i = 0; i < G.edge_count(); ++i) must_e[i] &= pe[i];
            }
        }
        for (int v = 0; v < n; ++v) {
            CHECK(P.on_vertex(v) == (bool)on_v[v] && P.must_vertex(v) == (bool)must_v[v]);
            CHECK(P.count_through_vertex(v) == thr_v[v]);
        }
        for (int i = 0; i < G.edge_count(); ++i) {
            CHECK(P.on_edge(i) == (bool)on_e[i] && P.must_edge(i) == (bool)must_e[i]);
            CHECK(P.count_through_edge(i) == thr_e[i]);
        }
        if (P.reachable()) {
            const std::vector<int> path = P.path();
            CHECK(path.front() == s && path.back() == t);
        }
        // count_shortest_paths（全頂点）と order_by_dist
        const SearchResult r = dijkstra(G, s);
        const std::vector<long long> cnt = count_shortest_paths(G, r, 0);
        for (int v = 0; v < n; ++v) CHECK(cnt[v] == (long long)Brute(G, s, v).paths.size());
        const std::vector<int> ord = order_by_dist(r);
        for (size_t k = 1; k < ord.size(); ++k)
            CHECK(r.dist[ord[k - 1]] < r.dist[ord[k]] || (r.dist[ord[k - 1]] == r.dist[ord[k]] && ord[k - 1] < ord[k]));
        int reached = 0;
        for (int v = 0; v < n; ++v) reached += r.reached(v);
        CHECK((int)ord.size() == reached);
    }
    // 重み 0 の辺: 閉路がなければ数えられる、閉路があれば -1
    {
        Graph G(4);
        G.add_edge(0, 1, 0), G.add_edge(0, 2, 0), G.add_edge(1, 3, 0), G.add_edge(2, 3, 0), G.add_edge(1, 2, 0);
        const auto c = count_shortest_paths(G, dijkstra(G, 0), 0);
        CHECK(c == (std::vector<long long>{1, 1, 2, 3}));
        G.add_edge(2, 1, 0);  // 1 ↔ 2 の重み 0 の閉路
        const auto c2 = count_shortest_paths(G, dijkstra(G, 0), 0);
        CHECK(c2[0] == 1 && c2[1] == -1 && c2[2] == -1 && c2[3] == -1);
    }
    // 多始点: 始点どうしが重み 0 でつながっていても、始点はそれぞれ 1 本と数える
    {
        Graph G(3);
        G.add_edge(0, 1, 0), G.add_edge(1, 2, 1), G.add_edge(0, 2, 1);
        const auto c = count_shortest_paths(G, dijkstra(G, std::vector<int>{0, 1}), 0);
        CHECK(c == (std::vector<long long>{1, 2, 3}));
    }
    std::cout << "[shortest_path_dag] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[shortest_path_dag] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    const int n = 200000;
    std::vector<int> us(n), vs(n);
    std::vector<long long> ws(n);
    for (int i = 0; i < n; ++i) us[i] = (int)(rng() % n), vs[i] = (int)(rng() % n), ws[i] = (long long)(rng() % 10 + 1);
    const Graph G = graph_from_edges(n, us, vs, ws, UNDIRECTED, 0);
    auto t0 = std::chrono::steady_clock::now();
    const ShortestPathDag P(G, 0, n - 1);
    std::cout << "[shortest_path_dag] benchmark: N = M = " << n << "  ShortestPathDag " << ms_since(t0) << " ms  (count " << P.count() << ")\n";
}

}  // namespace spd_test


namespace spd_demo {
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
    banner("最短路 DAG 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"count_shortest_paths(G, r, mod)", "vector<long long>", "始点から各頂点への最短路の本数（r は G の探索の結果）"},
        {"order_by_dist(r)", "vector<int>", "到達した頂点を距離の昇順に"},
        {"ShortestPathDag P(G, s, t, mod);", "", "s → t の最短路について、まとめて求める"},
        {"P.best / P.count()", "Cost / long long", "最短距離 / 最短路の本数"},
        {"P.on_vertex(v) / P.on_edge(id)", "bool", "どれかの最短路に乗るか"},
        {"P.must_vertex(v) / P.must_edge(id)", "bool", "すべての最短路が通るか"},
        {"P.count_through_vertex / _edge", "long long", "そこを通る最短路の本数"},
        {"P.path() / P.ds / P.dt", "", "最短路 1 本 / s からの距離 / t までの距離"},
    };
    cheat_sheet("API 一覧", rows, 38, 20);

    head("デモで使うグラフ（無向）: 0-1 (1), 0-2 (1), 1-3 (1), 2-3 (1), 3-4 (1), 0-4 (5)");
    note("    1");
    note("  /   \\");
    note(" 0     3 - 4      0 → 4 の最短は 3（0-1-3-4 と 0-2-3-4 の 2 本）");
    note("  \\   /");
    note("    2");
    RUN(Graph G(5));
    RUN(G.add_undirected(0, 1); G.add_undirected(0, 2); G.add_undirected(1, 3); G.add_undirected(2, 3); G.add_undirected(3, 4); G.add_undirected(0, 4, 5));

    api("【1】 count_shortest_paths(G, r, mod)", "vector<long long>",
        "始点から各頂点への最短路の本数（mod = 0 なら割らない）。r は G の BFS / ダイクストラの結果。ABC211 D", "O(N + M)");
    RUN(auto r = dijkstra(G, 0));
    SHOW(count_shortest_paths(G, r, 1000000007), "頂点 3, 4 へは 2 本");
    note("重みなしのグラフなら bfs(G, s) の結果を渡してもよい（G の重みが全部 1 のとき）。");

    api("【2】 ShortestPathDag P(G, s, t, mod = 998244353)", "", "s → t の最短路に乗る頂点・辺、必ず通る頂点・辺、通る本数をまとめて求める", "O((N + M) log M)");
    RUN(ShortestPathDag P(G, 0, 4));
    SHOW(P.best, "最短距離");
    SHOW(P.count(), "最短路の本数");
    SHOW(P.on_v, "最短路に乗る頂点（on_vertex(v)）");
    SHOW(P.on_e, "最短路に乗る辺（on_edge(id)）。辺 5（0-4）は乗らない");
    SHOW(P.must_v, "すべての最短路が通る頂点（0, 3, 4）");
    SHOW(P.must_e, "すべての最短路が通る辺 = 消すと最短距離が変わる辺（辺 4: 3-4）");
    SHOW(P.count_through_vertex(1), "頂点 1 を通る最短路の本数");
    SHOW(P.count_through_edge(4), "辺 4 を通る最短路の本数");
    SHOW(P.path(), "最短路を 1 本");
    note("「各辺を消したときの最短距離」: must_edge(id) でない辺は消しても変わらない（ABC375 G 型）。");

    api("【3】 order_by_dist(r)", "vector<int>", "到達した頂点を距離の昇順（同じなら番号順）に。距離の順に DP するとき", "O(N log N)");
    SHOW(order_by_dist(r), "距離 0, 1, 1, 2, 3 の順");
}

}  // namespace spd_demo

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
    spd_test::random_test();
    spd_test::benchmark();
    spd_demo::demo();
    return 0;
}
