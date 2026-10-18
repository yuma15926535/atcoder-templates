// =============================================================================
//  ベルマンフォード（負の重みがある最短路・−∞・負閉路）  bellman_ford / find_negative_cycle
//                                                                 --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 最短路・状態グラフ
//  依存       : 00_graph_base, 01_cycle（同梱済み。負閉路を Cycle で返すため）
//  頂点番号   : 0-indexed。重みは負でもよい
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  負閉路（一周すると重みの和が負になる閉路）があると、そこを何周もして距離をいくらでも小さくできる。
//  そうなる頂点の距離は −∞（-DIST_INF）にする。負閉路から行けない頂点の距離は正しく求まる。
//  無向グラフの負の辺は、行って戻るだけで負閉路になる（u → v → u）。
//
//  機能一覧                                                                       計算量
//    bellman_ford(G, s) / bellman_ford(G, {s1, s2, ...})   s からの最短距離（多始点も可）   O(NM)
//    bellman_ford_core(n, s / starts, adj)   隣接を自分で書く版（adj(v, push) で push(u, w)） O(n × 隣接の数)
//      結果 r: r.dist[v]（DIST_INF = 届かない、-DIST_INF = −∞）/ r.par[v] / r.par_edge[v]
//              r.reached(v) / r.neg_inf(v) / r.path(t) / r.negative_cycle（始点から行ける負閉路があるか）
//    find_negative_cycle(G)       グラフのどこかの負閉路を 1 つ（始点によらない）               O(NM)
//      結果 Cycle c: c.found() / c.vs（頂点の列）/ c.es（辺番号の列）
//
//  使い分け: 重みが 0 以上なら 02_dijkstra（ずっと速い）。全点対で N ≤ 400 程度なら 02_warshall_floyd
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤・閉路も入っている）
//    bf_test::random_test() / bf_test::benchmark() / bf_demo::demo()
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
#ifndef GG_CYCLE
#define GG_CYCLE
// ---- 閉路（グラフ・グリッド第二世代）  閉路を 1 つ見つける・最短の閉路・2-core ----
struct Cycle {
    std::vector<int> vs;  // 頂点 vs[0] → vs[1] → ... → vs[k-1] → vs[0]（同じ頂点は 2 回出てこない）
    std::vector<int> es;  // es[i] = vs[i] → vs[i+1] の辺番号（最後は vs[k-1] → vs[0]）
    bool found() const { return !vs.empty(); }
    int size() const { return (int)vs.size(); }
};
namespace gg_internal {
// 非再帰 DFS で閉路を探す。undirected なら「来た辺（同じ辺番号）」は戻らない
inline Cycle find_cycle_impl(const Graph& G, bool undirected) {
    const int n = G.size();
    std::vector<int> state(n, 0), par(n, -1), pe(n, -1), it(n, 0);  // state: 0 未訪問, 1 調べ中（スタック上）, 2 済み
    std::vector<int> st;
    for (int s = 0; s < n; ++s) {
        if (state[s]) continue;
        state[s] = 1, st.push_back(s);
        while (!st.empty()) {
            const int v = st.back();
            if (it[v] == (int)G[v].size()) {
                state[v] = 2, st.pop_back();
                continue;
            }
            const Edge& e = G[v][it[v]++];
            if (undirected && e.id == pe[v]) continue;
            const int u = e.to;
            if (state[u] == 0) {
                state[u] = 1, par[u] = v, pe[u] = e.id, st.push_back(u);
            } else if (state[u] == 1) {  // u はスタック上 → u → ... → v → u が閉路
                Cycle c;
                for (int x = v; x != u; x = par[x]) c.vs.push_back(x), c.es.push_back(pe[x]);
                c.vs.push_back(u);
                std::reverse(c.vs.begin(), c.vs.end());
                std::reverse(c.es.begin(), c.es.end());
                c.es.push_back(e.id);
                return c;
            }
        }
    }
    return Cycle{};
}
}  // namespace gg_internal
// 有向閉路を 1 つ（なければ空）。自己ループは長さ 1 の閉路
inline Cycle find_cycle_directed(const Graph& G) { return gg_internal::find_cycle_impl(G, false); }
// 無向グラフの閉路を 1 つ（なければ空 = 森）。自己ループは長さ 1、多重辺は長さ 2 の閉路
inline Cycle find_cycle_undirected(const Graph& G) { return gg_internal::find_cycle_impl(G, true); }
// s を通る最短（辺の本数）の有向閉路（なければ空）。ABC376 D
inline Cycle shortest_cycle_directed(const Graph& G, int s) {
    const int n = G.size();
    std::vector<int> dist(n, -1), par(n, -1), pe(n, -1), q{s};
    dist[s] = 0;
    for (size_t h = 0; h < q.size(); ++h) {
        const int v = q[h];
        for (const Edge& e : G[v]) {
            if (e.to == s) {  // v → s で閉じる（BFS 順なので最初に見つかったものが最短）
                Cycle c;
                for (int x = v; x != -1; x = par[x]) c.vs.push_back(x), c.es.push_back(pe[x]);
                std::reverse(c.vs.begin(), c.vs.end());
                std::reverse(c.es.begin(), c.es.end());
                c.es.erase(c.es.begin());  // 先頭は s の pe（-1）
                c.es.push_back(e.id);
                return c;
            }
            if (dist[e.to] == -1) dist[e.to] = dist[v] + 1, par[e.to] = v, pe[e.to] = e.id, q.push_back(e.to);
        }
    }
    return Cycle{};
}
// 無向グラフの 2-core: 次数 1 以下の頂点を削り続けて残る頂点に 1。閉路の上か、閉路どうしをつなぐ道の上の頂点。
// なもりグラフ（連結で辺 N 本）なら、ちょうど閉路の頂点になる。自己ループは次数 2、多重辺も別々に数える
inline std::vector<int> two_core(const Graph& G) {
    const int n = G.size();
    std::vector<int> deg(n, 0), q;
    for (const RawEdge& e : G.es) ++deg[e.from], ++deg[e.to];
    std::vector<int> in(n, 1);
    for (int v = 0; v < n; ++v)
        if (deg[v] <= 1) in[v] = 0, q.push_back(v);
    for (size_t h = 0; h < q.size(); ++h)
        for (const Edge& e : G[q[h]]) {
            const int u = e.to;
            if (in[u] && --deg[u] <= 1) in[u] = 0, q.push_back(u);
        }
    return in;
}
#endif  // GG_CYCLE
#ifndef GG_BELLMAN_FORD
#define GG_BELLMAN_FORD
// ---- ベルマンフォード（グラフ・グリッド第二世代）  負の重みがある最短路・負閉路の検出 ----
// 結果 r: r.dist[v] = 最短距離（届かなければ DIST_INF、負閉路を通っていくらでも小さくできるなら -DIST_INF = −∞）
// 無向グラフの負の辺は、行って戻るだけで負閉路になる（u → v → u）
struct BellmanFordResult {
    std::vector<Cost> dist;      // 最短距離。DIST_INF = 届かない、-DIST_INF = −∞
    std::vector<int> par;        // 直前の頂点（始点・未到達・−∞ は -1）
    std::vector<int> par_edge;   // 直前の辺の番号（同上は -1）
    bool negative_cycle = false; // 始点から行ける負閉路があるか（= どこかの頂点が −∞）
    bool reached(int v) const { return dist[v] != DIST_INF; }   // 届くか（−∞ も届く）
    bool neg_inf(int v) const { return dist[v] == -DIST_INF; }  // −∞ か
    std::vector<int> path(int t) const {  // 始点 → t の頂点列（届かない・−∞ なら空）
        if (!reached(t) || neg_inf(t)) return {};
        return restore_path(par, t);
    }
};
namespace gg_internal {
// adj(v, push): v から出る辺ごとに push(u, w, 辺番号) を呼ぶ
template <class F>
BellmanFordResult bellman_ford_impl(int n, const std::vector<int>& starts, F&& adj) {
    BellmanFordResult r;
    r.dist.assign(n, DIST_INF), r.par.assign(n, -1), r.par_edge.assign(n, -1);
    for (int s : starts) r.dist[s] = 0;
    std::vector<int> bad;  // n 回目でもまだ短くなる辺の行き先（負閉路の上か、その先）
    for (int round = 0; round < n; ++round) {
        bool updated = false;
        for (int v = 0; v < n; ++v) {
            const Cost dv = r.dist[v];
            if (dv == DIST_INF) continue;
            adj(v, [&](int u, Cost w, int id) {
                if (dv + w < r.dist[u]) {
                    updated = true;
                    if (round == n - 1) bad.push_back(u);
                    else r.dist[u] = dv + w, r.par[u] = v, r.par_edge[u] = id;
                }
            });
        }
        if (!updated) break;
    }
    if (bad.empty()) return r;
    // −∞ を、負閉路から行ける頂点すべてに広げる
    r.negative_cycle = true;
    std::vector<int> st;
    for (int u : bad)
        if (r.dist[u] != -DIST_INF) r.dist[u] = -DIST_INF, st.push_back(u);
    while (!st.empty()) {
        const int v = st.back();
        st.pop_back();
        adj(v, [&](int u, Cost, int) {
            if (r.dist[u] != -DIST_INF) r.dist[u] = -DIST_INF, st.push_back(u);
        });
    }
    for (int v = 0; v < n; ++v)
        if (r.dist[v] == -DIST_INF) r.par[v] = -1, r.par_edge[v] = -1;
    return r;
}
}  // namespace gg_internal
// 隣接を自分で書く版。adj(v, push) の中で push(u, w)（w は負でもよい）。O(n × 隣接の数)
template <class F>
BellmanFordResult bellman_ford_core(int n, const std::vector<int>& starts, F&& adj) {
    return gg_internal::bellman_ford_impl(n, starts, [&](int v, auto&& push) {
        adj(v, [&](int u, Cost w) { push(u, w, -1); });
    });
}
template <class F>
BellmanFordResult bellman_ford_core(int n, int s, F&& adj) { return bellman_ford_core(n, std::vector<int>{s}, adj); }
// グラフのベルマンフォード（辺の重み e.w は負でもよい）。O(NM)。負閉路がなければ N 回より早く終わることが多い
inline BellmanFordResult bellman_ford(const Graph& G, const std::vector<int>& starts) {
    return gg_internal::bellman_ford_impl(G.size(), starts, [&](int v, auto&& push) {
        for (const Edge& e : G[v]) push(e.to, e.w, e.id);
    });
}
inline BellmanFordResult bellman_ford(const Graph& G, int s) { return bellman_ford(G, std::vector<int>{s}); }
// グラフのどこかにある負閉路を 1 つ（始点によらない）。なければ c.found() == false。O(NM)
// c.vs = 頂点の列、c.es = 辺番号の列（vs[i] → vs[i+1]、最後は vs[k-1] → vs[0]）
inline Cycle find_negative_cycle(const Graph& G) {
    const int n = G.size();
    std::vector<Cost> dist(n, 0);  // すべての頂点から同時に始める（仮の始点からの重み 0 の辺）
    std::vector<int> par(n, -1), pe(n, -1);
    int last = -1;
    for (int round = 0; round < n; ++round) {
        last = -1;
        for (int v = 0; v < n; ++v)
            for (const Edge& e : G[v])
                if (dist[v] + e.w < dist[e.to]) dist[e.to] = dist[v] + e.w, par[e.to] = v, pe[e.to] = e.id, last = e.to;
        if (last == -1) return Cycle();
    }
    // n 回目でも更新された頂点から親を n 回たどると、必ず閉路の上に来る
    int x = last;
    for (int i = 0; i < n; ++i) x = par[x];
    Cycle c;
    for (int y = x;; y = par[y]) {
        c.vs.push_back(y), c.es.push_back(pe[y]);
        if (par[y] == x) break;
    }
    // ここまでで vs は x, par[x], ... の逆順。向きを直す（vs[i] → vs[i+1]）
    std::reverse(c.vs.begin(), c.vs.end());
    std::reverse(c.es.begin(), c.es.end());
    // es[i] は「vs[i] に入る辺」になっているので、1 つずらして「vs[i] から出る辺」にする
    std::rotate(c.es.begin(), c.es.begin() + 1, c.es.end());
    return c;
}
#endif  // GG_BELLMAN_FORD
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "bellman_ford";
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

namespace bf_test {

// 愚直: 辺の本数 k 以下の歩道の最小（k = 0 .. n-1）と、負の閉じた歩道（長さ n 以下）がある頂点
struct Brute {
    int n;
    std::vector<Cost> dist;      // 正解（−∞ は -DIST_INF）
    std::vector<int> on_neg;     // 負の閉じた歩道が通る頂点
};
static std::vector<int> neg_cycle_vertices(const Graph& G) {
    const int n = G.size();
    std::vector<int> res(n, 0);
    for (int c = 0; c < n; ++c) {
        std::vector<Cost> d(n, DIST_INF);
        d[c] = 0;
        Cost best = DIST_INF;
        for (int k = 0; k < n; ++k) {  // 長さ k + 1 の歩道
            std::vector<Cost> nd(n, DIST_INF);
            for (int v = 0; v < n; ++v)
                if (d[v] != DIST_INF)
                    for (const Edge& e : G[v]) nd[e.to] = std::min(nd[e.to], d[v] + e.w);
            d = nd;
            best = std::min(best, d[c]);
        }
        res[c] = (best < 0);
    }
    return res;
}
static std::vector<Cost> brute(const Graph& G, const std::vector<int>& starts) {
    const int n = G.size();
    std::vector<Cost> d(n, DIST_INF);
    for (int s : starts) d[s] = 0;
    for (int k = 0; k < n; ++k) {  // n - 1 回で十分（余分に回しても負閉路がなければ変わらない）
        std::vector<Cost> nd = d;
        for (int v = 0; v < n; ++v)
            if (d[v] != DIST_INF)
                for (const Edge& e : G[v]) nd[e.to] = std::min(nd[e.to], d[v] + e.w);
        if (k < n - 1) d = nd;
    }
    // 到達できる負閉路の頂点から行ける頂点は −∞
    const std::vector<int> neg = neg_cycle_vertices(G);
    std::vector<int> reach(n, 0), st;
    for (int s : starts)
        if (!reach[s]) reach[s] = 1, st.push_back(s);
    while (!st.empty()) {
        const int v = st.back();
        st.pop_back();
        for (const Edge& e : G[v])
            if (!reach[e.to]) reach[e.to] = 1, st.push_back(e.to);
    }
    std::vector<int> bad(n, 0);
    for (int c = 0; c < n; ++c)
        if (reach[c] && neg[c] && !bad[c]) bad[c] = 1, st.push_back(c);
    while (!st.empty()) {
        const int v = st.back();
        st.pop_back();
        for (const Edge& e : G[v])
            if (!bad[e.to]) bad[e.to] = 1, st.push_back(e.to);
    }
    for (int v = 0; v < n; ++v)
        if (bad[v]) d[v] = -DIST_INF;
    return d;
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 8), m = rnd(rng, 0, 14);
        const bool und = (it % 4 == 0);
        const int lo = (it % 3 == 0) ? -3 : (it % 3 == 1) ? -1 : 0;  // 負の辺の量を変える
        Graph G(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            const Cost w = rnd(rng, lo, 9);
            und ? G.add_undirected(a, b, w) : G.add_edge(a, b, w);
        }
        std::vector<int> starts;
        for (int k = rnd(rng, 1, 2); k > 0; --k) starts.push_back(rnd(rng, 0, n - 1));
        const std::vector<Cost> ex = brute(G, starts);
        const BellmanFordResult r = bellman_ford(G, starts);
        CHECK(r.dist == ex);
        bool any_neg = false;
        for (int v = 0; v < n; ++v) any_neg |= (ex[v] == -DIST_INF);
        CHECK(r.negative_cycle == any_neg);
        for (int v = 0; v < n; ++v) {
            CHECK(r.reached(v) == (ex[v] != DIST_INF) && r.neg_inf(v) == (ex[v] == -DIST_INF));
            if (!r.reached(v) || r.neg_inf(v)) {
                CHECK(r.path(v).empty() && r.par[v] == -1);
                continue;
            }
            // 経路: 始点から辺をたどり、重みの和が dist
            const std::vector<int> p = r.path(v);
            CHECK(!p.empty() && p.back() == v && std::find(starts.begin(), starts.end(), p[0]) != starts.end());
            Cost sum = 0;
            for (size_t k = 1; k < p.size(); ++k) {
                const int id = r.par_edge[p[k]];
                const RawEdge& e = G.es[id];
                CHECK(r.par[p[k]] == p[k - 1] && ((e.from == p[k - 1] && e.to == p[k]) || (und && e.to == p[k - 1] && e.from == p[k])));
                sum += e.w;
            }
            CHECK(sum == r.dist[v]);
        }
        // core 版も同じ
        if (starts.size() == 1) {
            const BellmanFordResult rc = bellman_ford_core(n, starts[0], [&](int v, auto&& push) {
                for (const Edge& e : G[v]) push(e.to, e.w);
            });
            CHECK(rc.dist == ex && rc.negative_cycle == r.negative_cycle);
            CHECK(bellman_ford(G, starts[0]).dist == ex);
        }
        // 負閉路を 1 つ
        const std::vector<int> neg = neg_cycle_vertices(G);
        const bool has = std::count(neg.begin(), neg.end(), 1) > 0;
        const Cycle c = find_negative_cycle(G);
        CHECK(c.found() == has);
        if (has) {
            const int k = c.size();
            CHECK((int)c.es.size() == k && k >= 1);
            Cost sum = 0;
            std::vector<int> seen(n, 0);
            for (int i = 0; i < k; ++i) {
                const int a = c.vs[i], b = c.vs[(i + 1) % k];
                const RawEdge& e = G.es[c.es[i]];
                CHECK((e.from == a && e.to == b) || (und && e.from == b && e.to == a));
                CHECK(!seen[a]);
                seen[a] = 1;
                sum += e.w;
            }
            CHECK(sum < 0);
        }
    }
    std::cout << "[bellman_ford] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[bellman_ford] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {2000, 5000}) {
        if (TEST_SCALE == 1 && n > 2000) break;
        // 最悪に近い形: 一直線の道を逆順に並べた辺（1 周で 1 つずつしか進まない）＋ ランダムな辺
        Graph G(n);
        for (int v = n - 2; v >= 0; --v) G.add_edge(v, v + 1, -1);
        for (int i = 0; i < n; ++i) G.add_edge((int)(rng() % n), (int)(rng() % n), (Cost)(rng() % 1000));
        auto t0 = std::chrono::steady_clock::now();
        const BellmanFordResult r = bellman_ford(G, 0);
        const long long t1 = ms_since(t0);
        G.add_edge(n - 1, 0, 0);  // 負閉路を作る
        t0 = std::chrono::steady_clock::now();
        const Cycle c = find_negative_cycle(G);
        std::cout << "[bellman_ford] benchmark: N = " << n << ", M = " << 2 * n << "（最悪に近い形）  bellman_ford " << t1 << " ms, find_negative_cycle " << ms_since(t0)
                  << " ms  (dist " << r.dist[n - 1] << ", cycle " << c.size() << ")\n";
    }
}

}  // namespace bf_test


namespace bf_demo {
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
struct DistVec {  // 距離の配列を「−∞」「DIST_INF」で表示するための入れ物（デモ用）
    std::vector<Cost> d;
};
[[maybe_unused]] static std::string to_s(const DistVec& x) {
    std::string s = "{";
    for (size_t i = 0; i < x.d.size(); ++i)
        s += (i ? ", " : "") + (x.d[i] == -DIST_INF ? std::string("-DIST_INF") : x.d[i] == DIST_INF ? std::string("DIST_INF") : std::to_string(x.d[i]));
    return s + "}";
}

static void demo() {
    banner("ベルマンフォード 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto r = bellman_ford(G, s);", "BellmanFordResult", "負の重みもある最短距離（多始点も可）"},
        {"r.dist[v] / r.path(t)", "", "距離（DIST_INF 届かない / -DIST_INF は −∞）/ 経路"},
        {"r.reached(v) / r.neg_inf(v)", "bool", "届くか / −∞ か"},
        {"r.negative_cycle", "bool", "始点から行ける負閉路があるか"},
        {"find_negative_cycle(G)", "Cycle", "どこかの負閉路を 1 つ（c.vs / c.es）"},
        {"bellman_ford_core(n, s, adj)", "BellmanFordResult", "隣接を自分で書く: push(u, w)"},
    };
    cheat_sheet("API 一覧", rows, 34, 20);

    head("デモで使うグラフ（有向）: 0→1 (4), 0→2 (5), 2→1 (-3), 1→3 (2), 3→4 (1), 4→3 (-2), 0→5 (1)");
    RUN(Graph G(6));
    RUN(G.add_edge(0, 1, 4); G.add_edge(0, 2, 5); G.add_edge(2, 1, -3); G.add_edge(1, 3, 2); G.add_edge(3, 4, 1); G.add_edge(4, 3, -2); G.add_edge(0, 5, 1));

    api("【1】 bellman_ford(G, s)", "BellmanFordResult", "負の辺があっても最短距離が求まる。負閉路から行ける頂点は −∞（-DIST_INF）", "O(NM)");
    RUN(auto r = bellman_ford(G, 0));
    show("r.dist", DistVec{r.dist}, "3 ↔ 4 が負閉路（1 + (-2) = -1）なので 3, 4 は −∞（-DIST_INF）");
    SHOW(r.dist[1], "0 → 2 → 1 で 5 + (-3) = 2");
    SHOW(r.path(1), "経路");
    SHOW(r.neg_inf(3), "−∞ か");
    SHOW(r.negative_cycle, "始点から行ける負閉路がある");
    note("出力が「−∞ なら inf と出す」なら if (r.neg_inf(t)) と書く（dist の値 -DIST_INF をそのまま使わない）。");

    api("【2】 find_negative_cycle(G)", "Cycle", "グラフのどこかにある負閉路を 1 つ（始点によらない）。なければ c.found() == false", "O(NM)");
    RUN(auto c = find_negative_cycle(G));
    SHOW(c.vs, "頂点");
    SHOW(c.es, "辺番号（vs[i] → vs[i+1]）");

    api("【3】 最大化は重みを負にして最小化（ABC061 D Score Attack）", "BellmanFordResult", "スコアの最大 = -(重みを負にした最短距離)。−∞（無限に増やせる）なら inf", "O(NM)");
    RUN(Graph S(3));
    RUN(S.add_edge(0, 1, -4); S.add_edge(1, 2, -3); S.add_edge(0, 2, -5));
    RUN(auto rs = bellman_ford(S, 0));
    SHOW(-rs.dist[2], "0 → 1 → 2 で 4 + 3 = 7");
    note("終点に行けない負閉路は関係ない（−∞ は「負閉路から行ける頂点」だけに付く）。");

    api("【4】 牛ゲー（差の制約の最大化）", "BellmanFordResult", "「x[v] - x[u] ≤ w」を辺 u → v（重み w）にすると、x[t] - x[s] の最大 = s から t の最短距離", "O(NM)");
    code("// 例: x1 - x0 ≤ 5, x2 - x1 ≤ 3, x2 - x0 ≤ 6, x0 - x2 ≤ -2（x2 ≥ x0 + 2）");
    RUN(Graph D(3));
    RUN(D.add_edge(0, 1, 5); D.add_edge(1, 2, 3); D.add_edge(0, 2, 6); D.add_edge(2, 0, -2));
    RUN(auto rd = bellman_ford(D, 0));
    SHOW(rd.dist[2], "x2 - x0 の最大");
    SHOW(rd.negative_cycle, "負閉路があれば条件を満たす x はない");
    note("「x[v] - x[u] ≥ w」は「x[u] - x[v] ≤ -w」に直して辺 v → u（重み -w）。等式は 2 本。");
    note("差が「ちょうど w」だけで不等式がないなら 03_weighted_union_find のほうが速い。");

    api("【5】 bellman_ford_core(n, s, adj)", "BellmanFordResult", "隣接を自分で書く。adj(v, push) の中で push(u, w)（w は負でもよい）", "O(n × 隣接の数)");
    RUN(auto rc = bellman_ford_core(G.size(), 0, [&](int v, auto&& push) { for (auto& e : G[v]) if (e.to != 4) push(e.to, e.w); }));
    SHOW(rc.dist, "頂点 4 に入らなければ負閉路がない");

    head("典型的な使い方");
    note("負の辺がある最短路     : bellman_ford(G, s).dist[t]（neg_inf(t) なら「いくらでも小さく」）");
    note("最大化（得点）         : 重みを負にして最短路。neg_inf なら「いくらでも大きく」");
    note("差の制約（牛ゲー）     : x[v] - x[u] ≤ w を辺 u → v（重み w）");
    note("負の辺が少しだけ・DAG  : 04_topo_sort の順に DP するほうが速い（DAG なら負の辺も問題ない）");
}

}  // namespace bf_demo

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
    bf_test::random_test();
    bf_test::benchmark();
    bf_demo::demo();
    return 0;
}
