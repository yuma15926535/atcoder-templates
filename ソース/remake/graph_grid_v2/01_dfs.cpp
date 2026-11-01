// =============================================================================
//  非再帰 DFS（深さ優先探索）  dfs / dfs_all / dfs_core           --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 01 探索
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  再帰 DFS と同じ順に訪れる（頂点に入ったときに訪問済みにし、隣接の順に進む）が、非再帰なので
//  N = 10^6 の一直線のグラフでも Visual Studio のスタック 1MB で落ちない。
//  すべての経路を調べる（戻るときに訪問済みを外す）バックトラックは、これではなく再帰で書く（デモ【5】の型）。
//
//  機能一覧                                                                       計算量
//    dfs(G, s) / dfs(G, {s1, s2, ...})   s から（starts の順に、未訪問なら）DFS        O(N + M)
//    dfs_all(G)                          番号の小さい頂点から全部の頂点を DFS            O(N + M)
//    dfs_core(n, s / starts, adj)            隣接を自分で書く版（adj(v, push) で push(u)）   O(N + 隣接の数)
//    結果 DfsResult d:
//      d.pre / d.post            行きがけ順 / 帰りがけ順
//      d.tour                    入るときと子から戻るたびに頂点を書く列（ABC213 D）
//      d.par / d.par_edge / d.depth   DFS 木の親 / 親から来た辺の番号 / 深さ（未到達は -1）
//      d.in / d.out              pre の何番目か。DFS 木の部分木 = pre[in[v] .. out[v])
//      d.roots / d.reached(v)    DFS を始めた頂点 / 訪れたか
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    dfs_test::random_test() / dfs_test::benchmark() / dfs_demo::demo()
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
#ifndef GG_DFS
#define GG_DFS
// ---- 非再帰 DFS（グラフ・グリッド第二世代）  再帰と同じ順に訪れる（入ったときに訪問済みにする） ----
struct DfsResult {
    std::vector<int> pre;       // 行きがけ順（訪れた順）
    std::vector<int> post;      // 帰りがけ順（その頂点から先を調べ終えた順）
    std::vector<int> tour;      // 入るときと、子から戻るたびに頂点を書く列（オイラーツアー。ABC213 D）
    std::vector<int> par;       // DFS 木の親（始点・未到達は -1）
    std::vector<int> par_edge;  // 親から来た辺の番号（始点・未到達は -1。dfs_core では -1）
    std::vector<int> depth;     // DFS 木の深さ（未到達は -1）
    std::vector<int> in, out;   // pre の何番目か。DFS 木の部分木 = pre[in[v] .. out[v])（未到達は -1）
    std::vector<int> roots;     // 実際に DFS を始めた頂点（訪問済みの始点は飛ばす）
    bool reached(int v) const { return depth[v] >= 0; }
};
// 本体: adj(v, push) で v の隣接 u ごとに push(u) を呼ぶ。starts の順に、未訪問なら DFS を始める
template <class F>
DfsResult dfs_core(int n, const std::vector<int>& starts, F&& adj) {
    DfsResult r;
    r.par.assign(n, -1), r.par_edge.assign(n, -1), r.depth.assign(n, -1), r.in.assign(n, -1), r.out.assign(n, -1);
    r.pre.reserve(n), r.post.reserve(n);
    struct Frame { int v, begin, cur, end; };  // buf[begin .. end) = v の隣接、cur = 次に見る位置
    std::vector<int> buf;
    std::vector<Frame> st;
    auto enter = [&](int v) {
        r.in[v] = (int)r.pre.size(), r.pre.push_back(v), r.tour.push_back(v);
        const int b = (int)buf.size();
        adj(v, [&](int u) { buf.push_back(u); });
        st.push_back({v, b, b, (int)buf.size()});
    };
    for (int s : starts) {
        if (r.depth[s] != -1) continue;
        r.roots.push_back(s);
        r.depth[s] = 0;
        enter(s);
        while (!st.empty()) {
            Frame& f = st.back();
            if (f.cur < f.end) {
                const int v = f.v, u = buf[f.cur++];
                if (r.depth[u] != -1) continue;
                r.depth[u] = r.depth[v] + 1, r.par[u] = v;
                enter(u);  // st が伸びるので f はこの後使わない
            } else {
                const int v = f.v;
                buf.resize(f.begin);
                r.out[v] = (int)r.pre.size(), r.post.push_back(v);
                st.pop_back();
                if (!st.empty()) r.tour.push_back(st.back().v);
            }
        }
    }
    return r;
}
template <class F>
DfsResult dfs_core(int n, int s, F&& adj) { return dfs_core(n, std::vector<int>{s}, adj); }
// グラフの DFS（隣接の順に訪れる。番号の小さい順にしたいときは先に G.sort_adj()）。par_edge も入る
inline DfsResult dfs(const Graph& G, const std::vector<int>& starts) {
    const int n = G.size();
    DfsResult r;
    r.par.assign(n, -1), r.par_edge.assign(n, -1), r.depth.assign(n, -1), r.in.assign(n, -1), r.out.assign(n, -1);
    r.pre.reserve(n), r.post.reserve(n);
    std::vector<int> st, it;
    for (int s : starts) {
        if (r.depth[s] != -1) continue;
        r.roots.push_back(s);
        r.depth[s] = 0, r.in[s] = (int)r.pre.size(), r.pre.push_back(s), r.tour.push_back(s);
        st.push_back(s), it.push_back(0);
        while (!st.empty()) {
            const int v = st.back();
            const std::vector<Edge>& a = G[v];
            if (it.back() < (int)a.size()) {
                const Edge& e = a[it.back()++];
                const int u = e.to;
                if (r.depth[u] != -1) continue;
                r.depth[u] = r.depth[v] + 1, r.par[u] = v, r.par_edge[u] = e.id, r.in[u] = (int)r.pre.size();
                r.pre.push_back(u), r.tour.push_back(u);
                st.push_back(u), it.push_back(0);
            } else {
                r.out[v] = (int)r.pre.size(), r.post.push_back(v);
                st.pop_back(), it.pop_back();
                if (!st.empty()) r.tour.push_back(st.back());
            }
        }
    }
    return r;
}
inline DfsResult dfs(const Graph& G, int s) { return dfs(G, std::vector<int>{s}); }
inline DfsResult dfs_all(const Graph& G) {  // 番号の小さい頂点から、未訪問なら DFS を始める（全頂点を訪れる）
    std::vector<int> all(G.size());
    for (int v = 0; v < G.size(); ++v) all[v] = v;
    return dfs(G, all);
}
#endif  // GG_DFS
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "dfs";
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

namespace dfs_test {

// 愚直解: 再帰 DFS
struct Brute {
    DfsResult r;
    Brute(const Graph& G, const std::vector<int>& starts, bool use_edge_id) {
        const int n = G.size();
        r.par.assign(n, -1), r.par_edge.assign(n, -1), r.depth.assign(n, -1), r.in.assign(n, -1), r.out.assign(n, -1);
        std::function<void(int)> rec = [&](int v) {
            r.in[v] = (int)r.pre.size(), r.pre.push_back(v), r.tour.push_back(v);
            for (const Edge& e : G[v]) {
                if (r.depth[e.to] != -1) continue;
                r.depth[e.to] = r.depth[v] + 1, r.par[e.to] = v;
                if (use_edge_id) r.par_edge[e.to] = e.id;
                rec(e.to);
                r.tour.push_back(v);
            }
            r.out[v] = (int)r.pre.size(), r.post.push_back(v);
        };
        for (int s : starts) {
            if (r.depth[s] != -1) continue;
            r.roots.push_back(s), r.depth[s] = 0;
            rec(s);
        }
    }
};
static bool same(const DfsResult& a, const DfsResult& b) {
    return a.pre == b.pre && a.post == b.post && a.tour == b.tour && a.par == b.par && a.par_edge == b.par_edge && a.depth == b.depth &&
           a.in == b.in && a.out == b.out && a.roots == b.roots;
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 25), m = rnd(rng, 0, 3 * n);
        const bool und = rng() & 1;
        Graph G(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            und ? G.add_undirected(a, b) : G.add_edge(a, b);
        }
        if (rng() & 1) G.sort_adj();
        std::vector<int> starts;
        for (int i = rnd(rng, 1, 4); i > 0; --i) starts.push_back(rnd(rng, 0, n - 1));
        const DfsResult d = dfs(G, starts);
        CHECK(same(d, Brute(G, starts, true).r));
        if (starts.size() == 1) CHECK(same(dfs(G, starts[0]), d));
        std::vector<int> all(n);
        std::iota(all.begin(), all.end(), 0);
        CHECK(same(dfs_all(G), Brute(G, all, true).r));
        // dfs_core は par_edge を持たない以外は同じ
        const DfsResult c = dfs_core(n, starts, [&](int v, auto&& push) {
            for (const Edge& e : G[v]) push(e.to);
        });
        CHECK(same(c, Brute(G, starts, false).r));
        for (int v = 0; v < n; ++v) CHECK(d.reached(v) == (d.depth[v] >= 0));
        CHECK((int)d.tour.size() == 2 * (int)d.pre.size() - (int)d.roots.size());
    }
    std::cout << "[dfs] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[dfs] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> us(n), vs(n);
        for (int i = 0; i < n; ++i) us[i] = (int)(rng() % n), vs[i] = (int)(rng() % n);
        const Graph G = graph_from_edges(n, us, vs, UNDIRECTED, 0);
        Graph P(n);  // 一直線
        for (int i = 0; i + 1 < n; ++i) P.add_undirected(i, i + 1);
        auto t0 = std::chrono::steady_clock::now();
        const DfsResult a = dfs_all(G);
        const long long t1 = ms_since(t0);
        t0 = std::chrono::steady_clock::now();
        const DfsResult b = dfs(P, 0);
        const long long t2 = ms_since(t0);
        std::cout << "[dfs] benchmark: N = M = " << n << "  dfs_all（ランダム）" << t1 << " ms, 一直線（深さ " << n - 1 << "）" << t2 << " ms  (check "
                  << a.post[0] + b.depth[n - 1] << ")\n";
    }
}

}  // namespace dfs_test


namespace dfs_demo {
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
    banner("非再帰 DFS 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto d = dfs(G, s);", "DfsResult", "s から DFS（隣接の順。番号順にしたいなら先に G.sort_adj()）"},
        {"dfs(G, {s1, s2}) / dfs_all(G)", "DfsResult", "複数の始点 / 全頂点"},
        {"dfs_core(n, s / starts, adj)", "DfsResult", "隣接を自分で書く: adj(v, push) で push(u)"},
        {"d.pre / d.post / d.tour", "vector<int>", "行きがけ順 / 帰りがけ順 / 戻りも書く列"},
        {"d.par / d.par_edge / d.depth", "vector<int>", "DFS 木の親 / 親から来た辺 / 深さ"},
        {"d.in / d.out", "vector<int>", "DFS 木の部分木 = pre[in[v] .. out[v])"},
        {"d.roots / d.reached(v)", "", "DFS を始めた頂点 / 訪れたか"},
    };
    cheat_sheet("API 一覧", rows, 40, 16);

    head("デモで使うグラフ（無向）: 0-1, 0-2, 1-3, 1-4, 2-4, 5-6");
    RUN(Graph G(7));
    RUN(G.add_undirected(0, 1); G.add_undirected(0, 2); G.add_undirected(1, 3); G.add_undirected(1, 4); G.add_undirected(2, 4); G.add_undirected(5, 6));

    api("【1】 dfs(G, s)", "DfsResult", "s から DFS。再帰と同じ順（入ったときに訪問済み・隣接の順）", "O(N + M)");
    RUN(auto d = dfs(G, 0));
    SHOW(d.pre, "行きがけ順: 0 → 1 → 3 → 4 → 2（2 は 4 から入る）");
    SHOW(d.post, "帰りがけ順");
    SHOW(d.tour, "戻りも書く（ABC213 D の出力そのもの。木なら長さ 2N - 1）");
    SHOW(d.par, "DFS 木の親（届かない 5, 6 は -1）");
    SHOW(d.par_edge, "親から来た辺の番号");
    SHOW(d.depth, "DFS 木の深さ");
    SHOW(d.in, "pre の何番目か");
    SHOW(d.out, "DFS 木の部分木の終わり");

    api("【2】 dfs(G, {s1, s2, ...}) / dfs_all(G)", "DfsResult", "starts の順に、まだ訪れていなければ DFS を始める / 全頂点（番号順）", "O(N + M)");
    RUN(auto da = dfs_all(G));
    SHOW(da.roots, "DFS を始めた頂点（成分ごとに 1 つ）");
    SHOW(da.pre, "全部の頂点の行きがけ順");

    api("【3】 dfs_core(n, s / starts, adj)", "DfsResult", "隣接を自分で書く（par_edge は -1 のまま）。adj(v, push) で push(u)", "O(n + 隣接の数)");
    code("// 例: 大きい番号の頂点から先に訪れる");
    RUN(auto dc = dfs_core(G.size(), 0, [&](int v, auto&& push) { for (int k = (int)G[v].size() - 1; k >= 0; --k) push(G[v][k].to); }));
    SHOW(dc.pre, "0 → 2 → 4 → 1 → 3");

    api("【4】 有向グラフ", "DfsResult", "辺の向きにたどる。帰りがけ順の逆はトポロジカル順（DAG のとき）", "O(N + M)");
    RUN(Graph D(4));
    RUN(D.add_edge(0, 1); D.add_edge(0, 2); D.add_edge(2, 1); D.add_edge(1, 3));
    RUN(auto dd = dfs_all(D));
    SHOW(dd.post, "帰りがけ順");
    note("トポロジカル順は 04_topo_sort の topo_sort(G) を使うのがわかりやすい。");

    head("【5】 すべての経路を調べる（バックトラック）は再帰で書く（N が小さいときだけ。ABC284 E）");
    code("std::vector<char> used(N);  long long cnt = 0;");
    code("auto rec = [&](auto&& self, int v) -> void {");
    code("    ++cnt;  used[v] = 1;                        // v で終わる単純な経路を 1 本数える");
    code("    for (auto& e : G[v]) if (!used[e.to]) self(self, e.to);");
    code("    used[v] = 0;                                // 戻るときに外す（ここが普通の DFS と違う）");
    code("};");
    code("rec(rec, 0);");
    {
        std::vector<char> used(G.size(), 0);
        long long cnt = 0;
        auto rec = [&](auto&& self, int v) -> void {
            ++cnt, used[v] = 1;
            for (auto& e : G[v])
                if (!used[e.to]) self(self, e.to);
            used[v] = 0;
        };
        rec(rec, 0);
        std::cout << "    頂点 0 から始まる単純な経路の数 = " << cnt << "\n";
    }
    note("Visual Studio で深い再帰をするときは、プロジェクトの設定でスタックサイズを増やす（リンカー → システム → スタックのサイズ）。");
}

}  // namespace dfs_demo

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
    dfs_test::random_test();
    dfs_test::benchmark();
    dfs_demo::demo();
    return 0;
}
