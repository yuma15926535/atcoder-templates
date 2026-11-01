// =============================================================================
//  時刻で重みが変わるダイクストラ  dijkstra_time / dijkstra_time_core     --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 最短路・状態グラフ
//  依存       : 00_graph_base, 02_dijkstra（同梱済み。ヒープを共有する）
//  頂点番号   : 0-indexed（入力が 1-indexed なら read_directed(N, M) などで読むと 0-indexed になる）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  「辺を通るのにかかる時間が、その頂点を出る時刻で変わる」最短路（電車の発車時刻・時刻で混み具合が変わる道）。
//  02_dijkstra の dijkstra_core では、隣接を列挙する関数に今の時刻（距離）が渡らないので、こちらを使う。
//    - 早く着いても損をしない（早く着いたら待てばよい）とき正しい。待つと得なときは「待つ時間 + 通る時間」を返す
//    - 時刻を状態に持つ必要はない（頂点ごとに最も早く着く時刻だけ考えればよい）
//
//  機能一覧                                                                       計算量
//    dijkstra_time(G, s / starts / {{s, 時刻}, ...}, cost, target = -1)   cost(e, t) = 時刻 t に出て辺 e を通る時間   O((N + M) log N)
//      結果 r: r.dist[v]（着く時刻。届かなければ DIST_INF）/ r.par / r.path(t) / r.par_edge[v] / r.path_edges(t)
//    dijkstra_time_core(n, s / starts, adj, target = -1)   adj(v, t, push): 時刻 t に v にいるとき push(u, かかる時間)
//    next_time(t, k, r = 0)        時刻 t 以降で「k で割って r 余る」最初の時刻（k 分おきの発車）   O(1)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤・ダイクストラも入っている）
//    dt_test::random_test() / dt_test::benchmark() / dt_demo::demo()
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
#ifndef GG_DIJKSTRA_TIME
#define GG_DIJKSTRA_TIME
// ---- 時刻で重みが変わるダイクストラ（グラフ・グリッド第二世代）  02_dijkstra の本体に「その頂点に着いた時刻」を渡す版 ----
// adj(v, t, push): 頂点 v に時刻 t（= 始点からの最短距離）で着いたとき、行ける頂点 u ごとに push(u, w) を呼ぶ（w >= 0。u に着くのは時刻 t + w）。
// 「早く着いたほうが損をしない」（早く着いたら待てばよい = 遅く出て早く着く辺がない）なら正しい。待つと得な辺は、待つ時間も w に入れる。
// target・starts（{頂点, 始まりの時刻} の組も可）・返り値は dijkstra_core と同じ
template <class F>
SearchResult dijkstra_time_core(int n, const std::vector<std::pair<int, Cost>>& starts, F&& adj, int target = -1) {
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
        adj(v, d, [&](int u, Cost w) {
            if (w < DIST_INF && d + w < dist[u]) dist[u] = d + w, par[u] = v, pq.push(d + w, u);
        });
    }
    return r;
}
template <class F>
SearchResult dijkstra_time_core(int n, const std::vector<int>& starts, F&& adj, int target = -1) {
    std::vector<std::pair<int, Cost>> st;
    st.reserve(starts.size());
    for (int s : starts) st.emplace_back(s, 0);
    return dijkstra_time_core(n, st, adj, target);
}
template <class F>
SearchResult dijkstra_time_core(int n, int s, F&& adj, int target = -1) {
    return dijkstra_time_core(n, std::vector<std::pair<int, Cost>>{{s, 0}}, adj, target);
}
// グラフ版の結果: SearchResult（dist / par / path / reached）に、親から来た辺の番号 par_edge を足したもの
struct TimeSearchResult : SearchResult {
    std::vector<int> par_edge;  // par_edge[v] = 最短で v に着くときに最後に通った辺の番号（始点・未到達は -1）
    std::vector<int> path_edges(int t) const {  // 始点 → t で通る辺の番号の列（到達不能なら空）
        std::vector<int> es;
        if (dist[t] == DIST_INF) return es;
        for (int v = t; par_edge[v] != -1; v = par[v]) es.push_back(par_edge[v]);
        std::reverse(es.begin(), es.end());
        return es;
    }
};
// cost(e, t): 時刻 t に e の始点を出て、辺 e（G[v] の要素。e.to が行き先、e.id で辺の情報を引ける）を通って e.to に着くまでの時間
// （待ち時間を含む。0 以上。今は通れない・通らないなら DIST_INF）
template <class F>
TimeSearchResult dijkstra_time(const Graph& G, const std::vector<std::pair<int, Cost>>& starts, F&& cost, int target = -1) {
    const int n = G.size();
    TimeSearchResult r;
    r.dist.assign(n, DIST_INF), r.par.assign(n, -1), r.par_edge.assign(n, -1);
    gg_internal::DistHeap pq;
    for (const auto& s : starts)
        if (s.second < r.dist[s.first]) r.dist[s.first] = s.second, pq.push(s.second, s.first);
    while (!pq.empty()) {
        const std::pair<Cost, int> top = pq.pop();
        const Cost d = top.first;
        const int v = top.second;
        if (d != r.dist[v]) continue;
        if (v == target) break;
        for (const Edge& e : G[v]) {
            const Cost w = cost(e, d);
            if (w < DIST_INF && d + w < r.dist[e.to]) r.dist[e.to] = d + w, r.par[e.to] = v, r.par_edge[e.to] = e.id, pq.push(d + w, e.to);
        }
    }
    return r;
}
template <class F>
TimeSearchResult dijkstra_time(const Graph& G, const std::vector<int>& starts, F&& cost, int target = -1) {
    std::vector<std::pair<int, Cost>> st;
    st.reserve(starts.size());
    for (int s : starts) st.emplace_back(s, 0);
    return dijkstra_time(G, st, cost, target);
}
template <class F>
TimeSearchResult dijkstra_time(const Graph& G, int s, F&& cost, int target = -1) {
    return dijkstra_time(G, std::vector<std::pair<int, Cost>>{{s, 0}}, cost, target);
}
// 時刻 t 以降で「k で割ると r 余る」最初の時刻（k 分おきに出る電車の次の発車。ABC192 E）。k >= 1、t >= 0
inline Cost next_time(Cost t, Cost k, Cost r = 0) {
    r %= k;
    if (r < 0) r += k;
    const Cost q = t % k;
    return t + ((r - q) % k + k) % k;
}
#endif  // GG_DIJKSTRA_TIME
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "dijkstra_time";
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

namespace dt_test {

// 例 1: 辺 i は k[i] 分おきに発車して t[i] 分かかる（ABC192 E 型）
// 例 2: 時刻 x に出ると c[i] + d[i] / (x + 1) 分かかる。好きなだけ待ってよい（ABC204 E 型）
static Cost rush(Cost t, Cost c, Cost d) {  // 時刻 t 以降に出て、いちばん早く着く時刻
    Cost s = (Cost)std::sqrt((double)d);
    while (s * s > d) --s;
    while ((s + 1) * (s + 1) <= d) ++s;
    Cost best = t + c + d / (t + 1);
    for (Cost x1 = s - 1; x1 <= s + 2; ++x1) {
        const Cost x = std::max(t, x1 - 1);
        best = std::min(best, x + c + d / (x + 1));
    }
    return best;
}
static Cost rush_brute(Cost t, Cost c, Cost d) {
    Cost best = DIST_INF;
    for (Cost x = t; x <= t + d + 2; ++x) best = std::min(best, x + c + d / (x + 1));
    return best;
}

// 愚直: ベルマンフォードのように「着く時刻」を何度も緩和する（早く着いても損をしないので、これで最も早い時刻になる）
template <class F>
static std::vector<Cost> brute(const Graph& G, const std::vector<std::pair<int, Cost>>& st, F cost) {
    const int n = G.size();
    std::vector<Cost> d(n, DIST_INF);
    for (auto [s, t] : st) d[s] = std::min(d[s], t);
    for (int it = 0; it <= n; ++it)
        for (int v = 0; v < n; ++v) {
            if (d[v] == DIST_INF) continue;
            for (const Edge& e : G[v]) {
                const Cost w = cost(e, d[v]);
                if (w < DIST_INF) d[e.to] = std::min(d[e.to], d[v] + w);
            }
        }
    return d;
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int t = 0; t < 2000; ++t) CHECK(next_time(t, 7, 3) >= t && next_time(t, 7, 3) % 7 == 3 && next_time(t, 7, 3) - t < 7);
    for (int t = 0; t < 300; ++t)
        for (int c = 0; c < 3; ++c)
            for (int d = 0; d < 200; d += 7) CHECK(rush(t, c, d) == rush_brute(t, c, d));
    for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 9), m = rnd(rng, 0, 20);
        Graph G(n);
        std::vector<Cost> K(m), T(m), C(m), D(m);
        for (int i = 0; i < m; ++i) {
            const int u = rnd(rng, 0, n - 1), v = rnd(rng, 0, n - 1);
            if (it % 2) G.add_undirected(u, v, rnd(rng, 0, 9));
            else G.add_edge(u, v, rnd(rng, 0, 9));
            K[i] = rnd(rng, 1, 6), T[i] = rnd(rng, 0, 9), C[i] = rnd(rng, 0, 5), D[i] = rnd(rng, 0, 60);
        }
        std::vector<std::pair<int, Cost>> st;
        const int ns = rnd(rng, 1, 2);
        for (int k = 0; k < ns; ++k) st.push_back({rnd(rng, 0, n - 1), rnd(rng, 0, 5)});
        auto train = [&](const Edge& e, Cost t) { return next_time(t, K[e.id]) - t + T[e.id]; };
        auto rh = [&](const Edge& e, Cost t) { return rush(t, C[e.id], D[e.id]) - t; };
        auto closed = [&](const Edge& e, Cost t) { return (t % 4 == 3 && e.id % 2 == 0) ? Cost(DIST_INF) : Cost(e.w); };  // 時刻によって通れない（待てば通れる形ではないので FIFO でない。答えの一致は見ない）
        auto stat = [&](const Edge& e, Cost) { return e.w; };
        // (1) 発車時刻つき
        {
            const TimeSearchResult r = dijkstra_time(G, st, train);
            CHECK(r.dist == brute(G, st, train));
            for (int v = 0; v < n; ++v) {
                if (!r.reached(v)) continue;
                // 経路をたどって時刻を計算し直すと dist と一致する
                const std::vector<int> p = r.path(v), pe = r.path_edges(v);
                CHECK(p.size() == pe.size() + 1);
                Cost t0 = DIST_INF;
                for (auto [s, t] : st)
                    if (s == p[0]) t0 = std::min(t0, t);
                for (size_t k = 0; k < pe.size(); ++k) {
                    const RawEdge& e = G.es[pe[k]];
                    CHECK((e.from == p[k] && e.to == p[k + 1]) || (e.to == p[k] && e.from == p[k + 1]));
                    t0 += train(Edge{p[k + 1], pe[k], e.w}, t0);
                }
                CHECK(t0 == r.dist[v]);
            }
            // core 版も同じ
            const SearchResult rc = dijkstra_time_core(n, st, [&](int v, Cost t, auto&& push) {
                for (const Edge& e : G[v]) push(e.to, train(e, t));
            });
            CHECK(rc.dist == r.dist);
        }
        // (2) ラッシュアワー（待ってよい）
        {
            const TimeSearchResult r = dijkstra_time(G, st, rh);
            CHECK(r.dist == brute(G, st, rh));
        }
        // (3) 時刻によらない重みなら dijkstra と同じ
        {
            const TimeSearchResult r = dijkstra_time(G, st, stat);
            CHECK(r.dist == dijkstra(G, st).dist);
            const int tg = rnd(rng, 0, n - 1);
            const TimeSearchResult r2 = dijkstra_time(G, st, stat, tg);
            CHECK(r2.dist[tg] == r.dist[tg]);
        }
        // (4) DIST_INF を返す辺は通らない（落ちない・あふれない）
        {
            const TimeSearchResult r = dijkstra_time(G, st, closed);
            for (int v = 0; v < n; ++v) CHECK(r.dist[v] >= 0);
        }
        // (5) 始点 1 つの形
        {
            const int s = rnd(rng, 0, n - 1);
            CHECK(dijkstra_time(G, s, train).dist == brute(G, {{s, 0}}, train));
            CHECK(dijkstra_time(G, std::vector<int>{s}, train).dist == brute(G, {{s, 0}}, train));
            CHECK(dijkstra_time_core(n, s, [&](int v, Cost t, auto&& push) { for (const Edge& e : G[v]) push(e.to, train(e, t)); }).dist == brute(G, {{s, 0}}, train));
        }
    }
    std::cout << "[dijkstra_time] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[dijkstra_time] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    const int n = 200000, m = (TEST_SCALE == 1 ? 200000 : 400000);
    Graph G(n);
    std::vector<Cost> K(m), T(m);
    for (int i = 0; i < m; ++i) {
        G.add_undirected((int)(rng() % n), (int)(rng() % n), 1);
        K[i] = (Cost)(rng() % 100000) + 1, T[i] = (Cost)(rng() % 1000000000);
    }
    auto t0 = std::chrono::steady_clock::now();
    const TimeSearchResult r = dijkstra_time(G, 0, [&](const Edge& e, Cost t) { return next_time(t, K[e.id]) - t + T[e.id]; });
    const long long t1 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const SearchResult r2 = dijkstra(G, 0);
    std::cout << "[dijkstra_time] benchmark: N = " << n << ", M = " << m << "  発車時刻つき " << t1 << " ms（参考: 重み 1 の dijkstra " << ms_since(t0)
              << " ms）  (check " << r.dist[n - 1] % 1000 + r2.dist[n - 1] << ")\n";
}

}  // namespace dt_test


namespace dt_demo {
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
    banner("時刻で重みが変わるダイクストラ 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto r = dijkstra_time(G, s, cost);", "TimeSearchResult", "cost(e, t) = 時刻 t に出て辺 e を通るのにかかる時間"},
        {"dijkstra_time(G, {{s, 時刻}, ...}, cost)", "", "始まりの時刻を始点ごとに決める（多始点）"},
        {"r.dist[v] / r.reached(v)", "Cost / bool", "v に着くいちばん早い時刻（届かなければ DIST_INF）"},
        {"r.path(t) / r.path_edges(t)", "vector<int>", "通る頂点 / 通る辺の番号"},
        {"dijkstra_time_core(n, s, adj)", "SearchResult", "adj(v, t, push): 時刻 t に v にいるとき push(u, かかる時間)"},
        {"next_time(t, k, r = 0)", "Cost", "時刻 t 以降で k で割って r 余る最初の時刻"},
    };
    cheat_sheet("API 一覧", rows, 44, 21);

    head("デモで使うグラフ: 0 → 1 → 2、0 → 2（辺番号 0, 1, 2）");
    RUN(Graph G(3));
    RUN(G.add_edge(0, 1); G.add_edge(1, 2); G.add_edge(0, 2));

    api("【1】 電車: 辺 i は K[i] 分おきに発車し T[i] 分かかる（ABC192 E）", "TimeSearchResult", "待ち時間（次の発車まで）+ 乗る時間 を返す", "O((N + M) log N)");
    RUN(std::vector<Cost> K = {3, 5, 10}, T = {2, 2, 9});
    code("auto r = dijkstra_time(G, 0, [&](const Edge& e, Cost t) { return next_time(t, K[e.id]) - t + T[e.id]; });");
    auto r = dijkstra_time(G, 0, [&](const Edge& e, Cost t) { return next_time(t, K[e.id]) - t + T[e.id]; });
    SHOW(r.dist, "0 → 1: 時刻 2。1 → 2: 時刻 5 に出て 7。直通は 9");
    SHOW(r.path(2), "");
    SHOW(r.path_edges(2), "辺 0, 1 を使う");

    api("【2】 時刻 x に出ると C + floor(D / (x + 1)) かかる。好きなだけ待てる（ABC204 E）", "TimeSearchResult", "x ≥ t で x + C + D / (x + 1) を最小にする x を選ぶ（√D の近く）", "O((N + M) log N)");
    code("auto rush = [](Cost t, Cost c, Cost d) {   // 時刻 t 以降に出て、いちばん早く着く時刻");
    code("    Cost s = (Cost)std::sqrt((double)d);");
    code("    while (s * s > d) --s;");
    code("    while ((s + 1) * (s + 1) <= d) ++s;");
    code("    Cost best = t + c + d / (t + 1);");
    code("    for (Cost x1 = s - 1; x1 <= s + 2; ++x1) { Cost x = std::max(t, x1 - 1); best = std::min(best, x + c + d / (x + 1)); }");
    code("    return best;");
    code("};");
    RUN(std::vector<Cost> C = {1, 1, 1}, D = {100, 100, 1000});
    code("auto r2 = dijkstra_time(G, 0, [&](const Edge& e, Cost t) { return rush(t, C[e.id], D[e.id]) - t; });");
    auto r2 = dijkstra_time(G, 0, [&](const Edge& e, Cost t) { return dt_test::rush(t, C[e.id], D[e.id]) - t; });
    SHOW(r2.dist, "0 → 1: 時刻 9 に出て 20。1 → 2: すぐ出て 20 + 1 + 4 = 25（直通は 63）");

    api("【3】 本体 dijkstra_time_core（グリッド・状態グラフでも同じ）", "SearchResult", "adj(v, t, push) の t が今の時刻。push(u, w) で u に時刻 t + w に着く", "O((N + M) log N)");
    code("auto r3 = dijkstra_time_core(3, 0, [&](int v, Cost t, auto&& push) {");
    code("    for (const Edge& e : G[v]) push(e.to, next_time(t, K[e.id]) - t + T[e.id]);");
    code("});");
    auto r3 = dijkstra_time_core(3, 0, [&](int v, Cost t, auto&& push) {
        for (const Edge& e : G[v]) push(e.to, next_time(t, K[e.id]) - t + T[e.id]);
    });
    SHOW(r3.dist, "【1】と同じ");

    head("注意");
    note("「早く着いても損をしない」（遅く出たほうが早く着く辺がない）ことが前提。待てる問題なら、待つ時間も返せば必ず成り立つ。");
    note("通れない辺は DIST_INF を返す。時刻によらない重みなら 02_dijkstra の dijkstra でよい。");
    note("時刻 mod K で状態が変わる（その時刻にしか通れない・待てない）ときは 00_state_id で (頂点, 時刻 mod K) を状態にする。");
}

}  // namespace dt_demo

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
    dt_test::random_test();
    dt_test::benchmark();
    dt_demo::demo();
    return 0;
}
