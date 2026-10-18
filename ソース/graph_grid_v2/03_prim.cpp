// =============================================================================
//  プリム法（最小・最大全域木。密なグラフ用の O(N^2) 版つき）  prim / prim_dense   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 03 Union-Find・最小全域木
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed。辺の向きは無視する。非連結なら全域森（connected = false）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  使い分け: 辺が入力で与えられるなら 03_kruskal でよい。「全部の組の間に辺がある（重みは式で決まる）」
//            ときは prim_dense(N, cost) を使う。辺を作らないので N = 5000（辺 1250 万本）でも軽い
//
//  機能一覧                                                                       計算量
//    prim(G) / prim_max(G)                     最小 / 最大全域木（森）              O((N + M) log M)
//    prim_dense(N, cost) / prim_dense_max(N, cost)   cost(u, v) で重み（辺なしは DIST_INF）  O(N^2)
//    結果 PrimResult r:
//      r.cost        選んだ辺の重みの合計
//      r.par[v]      全域木での親（各成分の番号最小の頂点が根で -1）。木の辺は (v, par[v])
//      r.par_edge[v] / r.edges   親との辺の番号 / 選んだ辺の番号の一覧（prim(G) のとき）
//      r.connected / r.components   全域木になったか / 連結成分の数
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    prim_test::random_test() / prim_test::benchmark() / prim_demo::demo()
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
#ifndef GG_PRIM
#define GG_PRIM
// ---- プリム法（グラフ・グリッド第二世代）  最小（最大）全域木。とくに辺が N^2 本ある密なグラフ用 ----
// 非連結なら全域森（connected = false）。辺の向きは無視する
struct PrimResult {
    Cost cost = 0;                // 選んだ辺の重みの合計
    std::vector<int> par;         // 全域木での親（各成分の番号最小の頂点が根で -1）
    std::vector<int> par_edge;    // 親との辺の番号（prim(G) のとき。密な版は -1）
    std::vector<int> edges;       // 選んだ辺の番号（prim(G) のとき。選んだ順）
    bool connected = false;       // 全域木になったか
    int components = 0;           // 連結成分の数
};
namespace gg_internal {
struct PrimHeap {  // (重み, 頂点, 辺番号) の最小ヒープ（デバッグビルドでも O(log n)）
    struct Item { Cost w; int v, id; };
    std::vector<Item> h;
    bool empty() const { return h.empty(); }
    void push(Cost w, int v, int id) {
        size_t i = h.size();
        h.push_back({w, v, id});
        while (i > 0) {
            const size_t p = (i - 1) >> 1;
            if (h[p].w <= w) break;
            h[i] = h[p], i = p;
        }
        h[i] = {w, v, id};
    }
    Item pop() {
        const Item top = h[0], x = h.back();
        h.pop_back();
        const size_t n = h.size();
        if (n > 0) {
            size_t i = 0;
            for (size_t c = 1; c < n; c = 2 * i + 1) {
                if (c + 1 < n && h[c + 1].w < h[c].w) ++c;
                if (x.w <= h[c].w) break;
                h[i] = h[c], i = c;
            }
            h[i] = x;
        }
        return top;
    }
};
inline PrimResult prim_graph(const Graph& G, bool maximize) {
    const int n = G.size(), m = G.edge_count();
    PrimResult r;
    r.par.assign(n, -1), r.par_edge.assign(n, -1);
    // 辺の向きを無視した隣接（CSR。辺番号を両端に入れる）
    std::vector<int> start(n + 1, 0), ids(2 * (size_t)m);
    for (const RawEdge& e : G.es) ++start[e.from + 1], ++start[e.to + 1];
    for (int v = 0; v < n; ++v) start[v + 1] += start[v];
    {
        std::vector<int> pos(start.begin(), start.end() - 1);
        for (int id = 0; id < m; ++id) ids[pos[G.es[id].from]++] = id, ids[pos[G.es[id].to]++] = id;
    }
    std::vector<char> done(n, 0);
    PrimHeap pq;
    for (int s = 0; s < n; ++s) {
        if (done[s]) continue;
        ++r.components;
        pq.push(0, s, -1);
        while (!pq.empty()) {
            const PrimHeap::Item it = pq.pop();
            const int v = it.v;
            if (done[v]) continue;
            done[v] = 1;
            if (it.id >= 0) {
                const RawEdge& e = G.es[it.id];
                r.par[v] = (e.from == v) ? e.to : e.from, r.par_edge[v] = it.id, r.edges.push_back(it.id), r.cost += e.w;
            }
            for (int k = start[v]; k < start[v + 1]; ++k) {
                const RawEdge& e = G.es[ids[k]];
                const int u = (e.from == v) ? e.to : e.from;
                if (!done[u]) pq.push(maximize ? -e.w : e.w, u, ids[k]);
            }
        }
    }
    r.connected = (r.components <= 1);
    return r;
}
template <class F>
PrimResult prim_dense_impl(int n, F&& cost, bool maximize) {
    PrimResult r;
    r.par.assign(n, -1), r.par_edge.assign(n, -1);
    std::vector<Cost> best(n, DIST_INF);  // 木までの最良の重み（maximize なら符号を反転して持つ）
    std::vector<char> in(n, 0), has(n, 0);
    for (int it = 0; it < n; ++it) {
        int v = -1;
        for (int u = 0; u < n; ++u)
            if (!in[u] && has[u] && (v == -1 || best[u] < best[v])) v = u;
        if (v == -1) {  // 新しい成分を始める（番号最小の未使用の頂点）
            for (int u = 0; u < n; ++u)
                if (!in[u]) { v = u; break; }
            ++r.components;
        } else {
            r.cost += maximize ? -best[v] : best[v];
        }
        in[v] = 1;
        for (int u = 0; u < n; ++u) {
            if (in[u]) continue;
            const Cost w = cost(v, u);
            if (w == DIST_INF) continue;
            const Cost k = maximize ? -w : w;
            if (!has[u] || k < best[u]) best[u] = k, has[u] = 1, r.par[u] = v;
        }
    }
    r.connected = (r.components <= 1);
    return r;
}
}  // namespace gg_internal
// グラフの最小全域木（森）。O((N + M) log M)。辺が少ないなら 03_kruskal と同じ結果（どちらでもよい）
inline PrimResult prim(const Graph& G) { return gg_internal::prim_graph(G, false); }
inline PrimResult prim_max(const Graph& G) { return gg_internal::prim_graph(G, true); }
// 密なグラフ（全部の組に辺がある）用。cost(u, v) = u - v の重み（辺がなければ DIST_INF）。O(N^2)。
// 辺をグラフとして作らずに済む（N = 5000 なら辺 1250 万本でもメモリを使わない）。ABC282 E など
template <class F>
PrimResult prim_dense(int n, F&& cost) { return gg_internal::prim_dense_impl(n, cost, false); }
template <class F>
PrimResult prim_dense_max(int n, F&& cost) { return gg_internal::prim_dense_impl(n, cost, true); }
#endif  // GG_PRIM
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "prim";
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

namespace prim_test {

// 愚直: クラスカル法（重みの順に、ラベル書き換えでつなぐ）
static std::pair<Cost, int> kruskal_brute(int n, std::vector<std::tuple<Cost, int, int>> es, bool maximize) {
    std::sort(es.begin(), es.end());
    if (maximize) std::reverse(es.begin(), es.end());
    std::vector<int> lab(n);
    std::iota(lab.begin(), lab.end(), 0);
    Cost cost = 0;
    int comps = n;
    for (auto [w, a, b] : es) {
        if (lab[a] == lab[b]) continue;
        const int o = lab[b];
        for (int& x : lab)
            if (x == o) x = lab[a];
        cost += w, --comps;
    }
    return {cost, comps};
}
// par が全域森になっているか（根は各成分の番号最小の頂点）と、木の辺の重みの和
static bool forest_ok(int n, const std::vector<int>& par, const std::function<Cost(int, int)>& w, Cost cost, int comps) {
    Cost sum = 0;
    int roots = 0;
    for (int v = 0; v < n; ++v) {
        if (par[v] == -1) {
            ++roots;
            continue;
        }
        const Cost x = w(v, par[v]);
        if (x == DIST_INF) return false;
        sum += x;
    }
    if (sum != cost || roots != comps) return false;
    // 親をたどると必ず根に着き（閉路がない）、根はその成分で番号最小
    std::vector<int> root(n);
    for (int v = 0; v < n; ++v) {
        int x = v, k = 0;
        while (par[x] != -1 && k <= n) x = par[x], ++k;
        if (k > n) return false;
        root[v] = x;
    }
    for (int v = 0; v < n; ++v)
        if (root[v] > v) return false;
    return true;
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 10), m = rnd(rng, 0, 20);
        Graph G(n);
        std::vector<std::tuple<Cost, int, int>> es;
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            const Cost w = rnd(rng, -5, 9);
            (rng() & 1) ? G.add_undirected(a, b, w) : G.add_edge(a, b, w);
            es.emplace_back(w, a, b);
        }
        for (bool mx : {false, true}) {
            const PrimResult r = mx ? prim_max(G) : prim(G);
            const auto [ec, comps] = kruskal_brute(n, es, mx);
            CHECK(r.cost == ec && r.components == comps && r.connected == (comps <= 1));
            CHECK((int)r.edges.size() == n - comps);
            Cost sum = 0;
            for (int v = 0; v < n; ++v) {
                if (r.par[v] == -1) {
                    CHECK(r.par_edge[v] == -1);
                    continue;
                }
                const RawEdge& e = G.es[r.par_edge[v]];
                CHECK((e.from == v && e.to == r.par[v]) || (e.to == v && e.from == r.par[v]));
                sum += e.w;
            }
            CHECK(sum == r.cost);
            CHECK(forest_ok(n, r.par, [&](int a, int) { return G.es[r.par_edge[a]].w; }, r.cost, comps));
        }
        // 密な版: 隣接行列（一部は辺なし）
        std::vector<std::vector<Cost>> mat(n, std::vector<Cost>(n, DIST_INF));
        std::vector<std::tuple<Cost, int, int>> des;
        for (int a = 0; a < n; ++a)
            for (int b = a + 1; b < n; ++b)
                if (rng() % 4) {
                    const Cost w = rnd(rng, -5, 20);
                    mat[a][b] = mat[b][a] = w;
                    des.emplace_back(w, a, b);
                }
        for (bool mx : {false, true}) {
            auto cost = [&](int a, int b) { return mat[a][b]; };
            const PrimResult r = mx ? prim_dense_max(n, cost) : prim_dense(n, cost);
            const auto [ec, comps] = kruskal_brute(n, des, mx);
            CHECK(r.cost == ec && r.components == comps && r.connected == (comps <= 1) && r.edges.empty());
            CHECK(forest_ok(n, r.par, [&](int a, int b) { return mat[a][b]; }, r.cost, comps));
            for (int v = 0; v < n; ++v) CHECK(r.par_edge[v] == -1);
        }
    }
    {
        const PrimResult e = prim(Graph(0));
        CHECK(e.cost == 0 && e.components == 0 && e.connected);
    }
    std::cout << "[prim] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[prim] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    {
        const int n = 200000;
        std::vector<int> us(n), vs(n);
        std::vector<long long> ws(n);
        for (int i = 0; i < n; ++i) us[i] = (int)(rng() % n), vs[i] = (int)(rng() % n), ws[i] = (long long)(rng() % 1000000000);
        const Graph G = graph_from_edges(n, us, vs, ws, UNDIRECTED, 0);
        auto t0 = std::chrono::steady_clock::now();
        const PrimResult r = prim(G);
        std::cout << "[prim] benchmark: prim N = M = " << n << "  " << ms_since(t0) << " ms  (cost " << r.cost % 1000 << ")\n";
    }
    for (int n : {2000, 5000}) {
        if (TEST_SCALE == 1 && n > 2000) break;
        std::vector<long long> x(n), y(n);
        for (int i = 0; i < n; ++i) x[i] = (long long)(rng() % 1000000), y[i] = (long long)(rng() % 1000000);
        auto t0 = std::chrono::steady_clock::now();
        const PrimResult r = prim_dense(n, [&](int a, int b) { return std::abs(x[a] - x[b]) + std::abs(y[a] - y[b]); });
        std::cout << "[prim] benchmark: prim_dense N = " << n << "（完全グラフ、辺 " << (long long)n * (n - 1) / 2 << " 本）  " << ms_since(t0) << " ms  (cost "
                  << r.cost % 1000 << ")\n";
    }
}

}  // namespace prim_test


namespace prim_demo {
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
    banner("プリム法 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto r = prim(G); / prim_max(G);", "PrimResult", "最小 / 最大全域木（森）"},
        {"prim_dense(N, cost) / prim_dense_max(N, cost)", "PrimResult", "全部の組に辺（cost(u, v) で重み）。O(N^2)"},
        {"r.cost / r.par[v] / r.edges", "", "重みの合計 / 木での親 / 選んだ辺の番号"},
        {"r.connected / r.components", "", "全域木になったか / 連結成分の数"},
    };
    cheat_sheet("API 一覧", rows, 48, 14);

    api("【1】 prim_dense(N, cost)", "PrimResult", "全部の組 (u, v) の間に辺があり、重みが式で決まるとき。辺を作らずに O(N^2)", "O(N^2)");
    code("// 例: 点の間の重み = マンハッタン距離。4 点 (0,0) (2,0) (2,3) (5,3)");
    RUN(std::vector<long long> x = {0, 2, 2, 5}, y = {0, 0, 3, 3});
    RUN(auto r = prim_dense(4, [&](int a, int b) { return std::abs(x[a] - x[b]) + std::abs(y[a] - y[b]); }));
    SHOW(r.cost, "2 + 3 + 3");
    SHOW(r.par, "木の辺は (v, par[v])");
    note("辺がない組は DIST_INF を返す（「i と j がつながっていない」行列が与えられたときなど）。");

    api("【2】 prim_dense_max(N, cost)", "PrimResult", "最大全域木。ABC282 E: 重みが (A_i^A_j + A_j^A_i) mod M の完全グラフ", "O(N^2)");
    RUN(std::vector<long long> a = {3, 5, 6});
    RUN(auto rm = prim_dense_max(3, [&](int i, int j) { return a[i] ^ a[j]; }));
    SHOW(rm.cost, "XOR が大きい組: 3^5 = 6, 3^6 = 5, 5^6 = 3 → 6 + 5");

    api("【3】 prim(G) / prim_max(G)", "PrimResult", "辺が与えられたグラフ。結果は 03_kruskal と同じ（どちらでもよい）", "O((N + M) log M)");
    RUN(Graph G(4));
    RUN(G.add_undirected(0, 1, 3); G.add_undirected(1, 2, 1); G.add_undirected(0, 2, 2); G.add_undirected(2, 3, 5));
    RUN(auto rg = prim(G));
    SHOW(rg.cost, "1 + 2 + 5");
    SHOW(rg.edges, "選んだ辺の番号");
    SHOW(prim_max(G).cost, "最大全域木: 5 + 3 + 2");
    RUN(Graph H(4));
    RUN(H.add_undirected(0, 1, 1));
    SHOW(prim(H).components, "非連結なら全域森（成分 3 つ）");
    SHOW(prim(H).connected, "");

    head("典型的な使い方");
    note("点の間の距離が式で決まる完全グラフ : prim_dense(N, cost)");
    note("最大全域木                         : prim_dense_max / prim_max（03_kruskal なら kruskal_max）");
    note("辺が M 本で与えられる               : 03_kruskal の kruskal(G) でもよい（必ず使う辺があるなら kruskal(G, uf)）");
}

}  // namespace prim_demo

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
    prim_test::random_test();
    prim_test::benchmark();
    prim_demo::demo();
    return 0;
}
