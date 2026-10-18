// =============================================================================
//  トポロジカルソートと DAG 上の DP  topo_sort / topo_sort_lex_min / dag_longest_path / dag_count_paths
//                                                                 --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 04 DAG
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed。有向グラフ用（辺 u → v は「u が v より先」）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  機能一覧                                                                       計算量
//    topo_sort(G)              入次数 0 になった順（Kahn 法・キュー）                 O(N + M)
//    topo_sort_lex_min(G)      辞書順最小のトポロジカル順                             O((N + M) log N)
//    topo_sort_lex_max(G)      辞書順最大のトポロジカル順                             O((N + M) log N)
//      結果 TopoResult t: t.order（順序）, t.pos[v]（何番目か。並べられなければ -1）,
//                         t.is_dag（閉路がないか）, t.unique（順序が 1 通りか。ABC291 E）
//    dag_longest_path(G)       最長路（重みの和が最大。どこから始めてもよい）          O(N + M)
//      結果 DagLongest d: d.dist[v]（v で終わる最長路）, d.par, d.best, d.end, d.path(), d.path(v)
//    dag_count_paths(G, s, mod)  s から各頂点への路の本数（mod。mod = 0 なら割らない）   O(N + M)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    topo_test::random_test() / topo_test::benchmark() / topo_demo::demo()
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
#ifndef GG_TOPO_SORT
#define GG_TOPO_SORT
// ---- トポロジカルソート（Kahn 法）と DAG 上の DP（グラフ・グリッド第二世代）  有向グラフ用 ----
struct TopoResult {
    std::vector<int> order;  // トポロジカル順（辺 u → v なら u が先）。DAG でなければ並べられた頂点だけ
    std::vector<int> pos;    // pos[v] = order の何番目か（並べられなかった頂点 = 閉路の上か、閉路から来る頂点は -1）
    bool is_dag = true;      // 閉路がないか
    bool unique = true;      // 順序が 1 通りに決まるか（is_dag のときだけ意味がある。ABC291 E）
};
namespace gg_internal {
template <class Q>
TopoResult topo_impl(const Graph& G, Q& q) {  // q: push / top_pop / size を持つ入れ物
    const int n = G.size();
    TopoResult r;
    r.pos.assign(n, -1);
    r.order.reserve(n);
    std::vector<int> indeg = G.in_degrees();
    for (int v = 0; v < n; ++v)
        if (indeg[v] == 0) q.push(v);
    while (q.size()) {
        if (q.size() > 1) r.unique = false;
        const int v = q.top_pop();
        r.pos[v] = (int)r.order.size(), r.order.push_back(v);
        for (const Edge& e : G[v])
            if (--indeg[e.to] == 0) q.push(e.to);
    }
    r.is_dag = ((int)r.order.size() == n);
    if (!r.is_dag) r.unique = false;
    return r;
}
struct FifoQ {
    std::vector<int> a;
    size_t h = 0;
    void push(int v) { a.push_back(v); }
    size_t size() const { return a.size() - h; }
    int top_pop() { return a[h++]; }
};
// 2 分ヒープ（Cmp = std::less なら最小を、std::greater なら最大を取り出す）。
// std::priority_queue と違い、-D_GLIBCXX_DEBUG のデバッグビルドでも 1 回 O(log n)
template <class Cmp>
struct HeapQ {
    std::vector<int> h;
    Cmp before;  // before(a, b) = a を b より先に取り出す
    void push(int v) {
        size_t i = h.size();
        h.push_back(v);
        while (i > 0) {
            const size_t p = (i - 1) >> 1;
            if (!before(v, h[p])) break;
            h[i] = h[p], i = p;
        }
        h[i] = v;
    }
    size_t size() const { return h.size(); }
    int top_pop() {
        const int top = h[0], x = h.back();
        h.pop_back();
        const size_t n = h.size();
        if (n > 0) {
            size_t i = 0;
            for (size_t c = 1; c < n; c = 2 * i + 1) {
                if (c + 1 < n && before(h[c + 1], h[c])) ++c;
                if (!before(h[c], x)) break;
                h[i] = h[c], i = c;
            }
            h[i] = x;
        }
        return top;
    }
};
}  // namespace gg_internal
inline TopoResult topo_sort(const Graph& G) {  // 入次数 0 になった順（キュー）。O(N + M)
    gg_internal::FifoQ q;
    return gg_internal::topo_impl(G, q);
}
inline TopoResult topo_sort_lex_min(const Graph& G) {  // 辞書順最小（並べられる中で番号の小さい頂点から）。O((N + M) log N)
    gg_internal::HeapQ<std::less<int>> q;
    return gg_internal::topo_impl(G, q);
}
inline TopoResult topo_sort_lex_max(const Graph& G) {  // 辞書順最大
    gg_internal::HeapQ<std::greater<int>> q;
    return gg_internal::topo_impl(G, q);
}

// DAG の最長路（辺の重みの和が最大。始点・終点はどこでもよい）。dist[v] = v で終わる路の最大値（v 1 つだけなら 0）
struct DagLongest {
    bool is_dag = true;
    std::vector<Cost> dist;  // v で終わる最長路の長さ（重みの和。重みなしのグラフなら辺の本数）
    std::vector<int> par;    // 最長路での直前の頂点（-1 なら v が始まり）
    Cost best = 0;           // 全体の最長路の長さ
    int end = -1;            // 全体の最長路の終点
    std::vector<int> path(int v) const {  // v で終わる最長路の頂点列
        std::vector<int> p;
        for (; v != -1; v = par[v]) p.push_back(v);
        std::reverse(p.begin(), p.end());
        return p;
    }
    std::vector<int> path() const { return end < 0 ? std::vector<int>() : path(end); }  // 全体の最長路
};
inline DagLongest dag_longest_path(const Graph& G) {
    const int n = G.size();
    const TopoResult t = topo_sort(G);
    DagLongest d;
    d.is_dag = t.is_dag;
    d.dist.assign(n, 0), d.par.assign(n, -1);
    for (int u : t.order)
        for (const Edge& e : G[u])
            if (d.dist[u] + e.w > d.dist[e.to])
                d.dist[e.to] = d.dist[u] + e.w, d.par[e.to] = u;
    for (int v : t.order)
        if (d.end == -1 || d.dist[v] > d.best) d.best = d.dist[v], d.end = v;
    return d;
}
// DAG で s から各頂点への路の本数（mod。mod = 0 なら割らない）。DAG でない場合、閉路の先は正しくない
inline std::vector<long long> dag_count_paths(const Graph& G, int s, long long mod) {
    const TopoResult t = topo_sort(G);
    std::vector<long long> c(G.size(), 0);
    c[s] = 1;
    for (int u : t.order) {
        if (c[u] == 0) continue;
        for (const Edge& e : G[u]) {
            c[e.to] += c[u];
            if (mod) c[e.to] %= mod;
        }
    }
    return c;
}
#endif  // GG_TOPO_SORT
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "topo_sort";
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

namespace topo_test {

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 2500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 7), m = rnd(rng, 0, 10);
        const bool dag = rng() % 3 != 0;
        Graph G(n);
        std::vector<int> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        for (int i = 0; i < m; ++i) {
            int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            if (dag) {
                if (a == b) continue;
                if (a > b) std::swap(a, b);
                a = perm[a], b = perm[b];
            }
            G.add_edge(a, b, rnd(rng, -3, 9));
        }
        // 全順列を調べる
        std::vector<int> p(n);
        std::iota(p.begin(), p.end(), 0);
        std::vector<std::vector<int>> valid;
        do {
            std::vector<int> pos(n);
            for (int i = 0; i < n; ++i) pos[p[i]] = i;
            bool ok = true;
            for (const RawEdge& e : G.es) ok &= pos[e.from] < pos[e.to];
            if (ok) valid.push_back(p);
        } while (std::next_permutation(p.begin(), p.end()));
        const bool is_dag = !valid.empty();
        const TopoResult t = topo_sort(G), tmin = topo_sort_lex_min(G), tmax = topo_sort_lex_max(G);
        for (const TopoResult* r : {&t, &tmin, &tmax}) {
            CHECK(r->is_dag == is_dag && r->unique == (valid.size() == 1));
            for (int i = 0; i < (int)r->order.size(); ++i) CHECK(r->pos[r->order[i]] == i);
            int cnt = 0;
            for (int v = 0; v < n; ++v) cnt += (r->pos[v] >= 0);
            CHECK(cnt == (int)r->order.size());
            for (const RawEdge& e : G.es)  // 並べられた頂点の間の辺は順序を守る。並べられた頂点への辺は並べられた頂点から
                if (r->pos[e.to] >= 0) CHECK(r->pos[e.from] >= 0 && r->pos[e.from] < r->pos[e.to]);
            if (is_dag) CHECK((int)r->order.size() == n && std::find(valid.begin(), valid.end(), r->order) != valid.end());
        }
        if (is_dag) CHECK(tmin.order == valid.front() && tmax.order == valid.back());
        if (!is_dag) continue;
        // 最長路・路の本数: 全部の路を DFS で列挙（DAG なので有限）
        std::vector<Cost> best(n, 0);
        std::vector<long long> cnt(n, 0);
        const int s = rnd(rng, 0, n - 1);
        std::function<void(int, Cost, bool)> rec = [&](int v, Cost d, bool from_s) {
            best[v] = std::max(best[v], d);
            if (from_s) ++cnt[v];
            for (const Edge& e : G[v]) rec(e.to, d + e.w, from_s);
        };
        for (int v = 0; v < n; ++v) rec(v, 0, v == s);
        const DagLongest d = dag_longest_path(G);
        CHECK(d.is_dag && d.dist == best);
        Cost mx = 0;
        for (int v = 0; v < n; ++v) mx = std::max(mx, best[v]);
        CHECK(d.best == mx && d.dist[d.end] == mx);
        for (int v = 0; v < n; ++v) {  // 経路の重みの和が dist[v]
            const std::vector<int> path = d.path(v);
            CHECK(path.back() == v);
            Cost sum = 0;
            for (size_t k = 0; k + 1 < path.size(); ++k) {
                Cost w = -DIST_INF;
                for (const Edge& e : G[path[k]])
                    if (e.to == path[k + 1]) w = std::max(w, e.w);
                CHECK(w != -DIST_INF);
                sum += w;
            }
            CHECK(sum == d.dist[v]);
        }
        CHECK(dag_count_paths(G, s, 0) == cnt);
    }
    std::cout << "[topo_sort] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[topo_sort] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> us(n), vs(n);
        for (int i = 0; i < n; ++i) {
            int a = (int)(rng() % n), b = (int)(rng() % n);
            if (a > b) std::swap(a, b);
            if (a == b) b = (b + 1) % n, a = 0;
            us[i] = a, vs[i] = b;
        }
        const Graph G = graph_from_edges(n, us, vs, DIRECTED, 0);
        auto t0 = std::chrono::steady_clock::now();
        const TopoResult t = topo_sort(G);
        const long long t1 = ms_since(t0);
        t0 = std::chrono::steady_clock::now();
        const TopoResult t2 = topo_sort_lex_min(G);
        const long long t3 = ms_since(t0);
        t0 = std::chrono::steady_clock::now();
        const DagLongest d = dag_longest_path(G);
        std::cout << "[topo_sort] benchmark: N = M = " << n << "  topo_sort " << t1 << " ms, lex_min " << t3 << " ms, dag_longest_path " << ms_since(t0)
                  << " ms  (check " << t.is_dag + t2.is_dag + d.best << ")\n";
    }
}

}  // namespace topo_test


namespace topo_demo {
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
    banner("トポロジカルソート 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto t = topo_sort(G);", "TopoResult", "トポロジカル順（キュー）"},
        {"topo_sort_lex_min(G) / _lex_max(G)", "TopoResult", "辞書順最小 / 最大"},
        {"t.order / t.pos[v]", "vector<int>", "順序 / 何番目か（並べられなければ -1）"},
        {"t.is_dag / t.unique", "bool", "閉路がないか / 順序が 1 通りに決まるか"},
        {"auto d = dag_longest_path(G);", "DagLongest", "最長路（d.dist[v], d.best, d.path()）"},
        {"dag_count_paths(G, s, mod)", "vector<long long>", "s から各頂点への路の本数"},
    };
    cheat_sheet("API 一覧", rows, 38, 20);

    head("デモで使うグラフ（有向）: 3→1, 3→0, 1→2, 0→2, 2→4（「3 の後に 1 と 0」など）");
    RUN(Graph G(5));
    RUN(G.add_edge(3, 1); G.add_edge(3, 0); G.add_edge(1, 2); G.add_edge(0, 2); G.add_edge(2, 4));

    api("【1】 topo_sort(G)", "TopoResult", "辺 u → v なら u が先になる順。入次数 0 になった順に並べる", "O(N + M)");
    RUN(auto t = topo_sort(G));
    SHOW(t.order, "順序");
    SHOW(t.pos, "各頂点が何番目か");
    SHOW(t.is_dag, "閉路がない");
    SHOW(t.unique, "0 と 1 の順は決まらないので 1 通りではない");

    api("【2】 topo_sort_lex_min(G) / topo_sort_lex_max(G)", "TopoResult", "並べられる中で、番号の小さい（大きい）頂点を先にする（辞書順最小 / 最大）", "O((N + M) log N)");
    SHOW(topo_sort_lex_min(G).order, "辞書順最小");
    SHOW(topo_sort_lex_max(G).order, "辞書順最大");

    api("【3】 閉路があるとき", "TopoResult", "is_dag = false。order には並べられた頂点だけが入り、閉路の上と閉路から来る頂点は pos = -1", "O(N + M)");
    RUN(Graph C(4));
    RUN(C.add_edge(0, 1); C.add_edge(1, 2); C.add_edge(2, 1); C.add_edge(2, 3));
    RUN(auto tc = topo_sort(C));
    SHOW(tc.is_dag, "1 ↔ 2 が閉路");
    SHOW(tc.order, "並べられたのは 0 だけ");
    SHOW(tc.pos, "");

    api("【4】 dag_longest_path(G)", "DagLongest", "辺の重みの和が最大の路（どこから始めてもよい）。重みなしなら辺の本数", "O(N + M)");
    RUN(auto d = dag_longest_path(G));
    SHOW(d.dist, "v で終わる最長路の長さ");
    SHOW(d.best, "全体の最長路");
    SHOW(d.path(), "その路（3 → 1 → 2 → 4 など）");

    api("【5】 dag_count_paths(G, s, mod)", "vector<long long>", "s から各頂点への路の本数（mod = 0 なら割らない）", "O(N + M)");
    SHOW(dag_count_paths(G, 3, 0), "3 から 4 へは 2 本");

    head("DAG 上の DP の型（自分で書くとき）");
    code("auto t = topo_sort(G);");
    code("for (int v : t.order) for (auto& e : G[v]) dp[e.to] = max(dp[e.to], dp[v] + e.w);   // 前から配る");
    note("SCC で縮約した DAG（Phase 4）の上でも同じ形で使う。");
}

}  // namespace topo_demo

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
    topo_test::random_test();
    topo_test::benchmark();
    topo_demo::demo();
    return 0;
}
