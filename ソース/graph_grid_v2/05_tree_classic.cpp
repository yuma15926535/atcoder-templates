// =============================================================================
//  木の上の有名な最適化  最大独立集合・最小点被覆・最小支配集合・最大マッチング・独立集合の数   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 木
//  依存       : 00_graph_base, 05_rooted_tree（同梱済み）
//  頂点番号   : 0-indexed。G は木（森でもよい。成分ごとに解いて合わせる）。read_tree(N) で読んだグラフをそのまま渡す
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  一般のグラフでは難しい（NP 困難な）問題も、木なら木DP・貪欲で O(N) で解ける。選んだ頂点・辺の復元つき。
//  頂点の重み w を省くとすべて 1（個数の最適化）。旧テンプレートの「最大独立集合」「最大マッチング」の木DP を置き換える
//
//  機能一覧                                                                       計算量
//    tree_max_independent_set(G, w = {})   隣り合わない頂点の集合で重みの和が最大（負の重みは選ばない）  O(N)
//    tree_min_vertex_cover(G, w = {})      全部の辺を覆う頂点の集合で重みの和が最小（w ≥ 0）             O(N)
//    tree_min_dominating_set(G, w = {})    全部の頂点が「自分か隣」に入る集合で重みの和が最小（w ≥ 0）   O(N)
//      結果 TreeVertexSet: value（重みの和・個数）/ vs（選んだ頂点。昇順）
//    tree_max_matching(G)                  端点を共有しない辺の集合で本数が最大                          O(N)
//    tree_max_weight_matching(G)           同じく辺の重み e.w の和が最大（負の辺は使わない）              O(N)
//      結果 TreeMatching: value / es（辺の番号。昇順）/ pairs（両端 (u, v)、u < v）
//    tree_count_independent_sets(G, mod)   独立集合（空も含む）の個数 mod（EDPC P）                       O(N)
//    （重みなしなら 最小点被覆の大きさ = 最大マッチングの本数 = N - 最大独立集合の大きさ）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤・根付き木も入っている）
//    tc_test::random_test() / tc_test::benchmark() / tc_demo::demo()
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
#ifndef GG_TREE_CLASSIC
#define GG_TREE_CLASSIC
// ---- 木の上の有名な最適化（グラフ・グリッド第二世代）  最大独立集合・最小点被覆・最小支配集合・最大マッチング・独立集合の数え上げ ----
// G は木（森でもよい。成分ごとに解いて合わせる）。すべて非再帰・O(N)。頂点の重み w を省くと、すべて 1（個数の最適化）
struct TreeVertexSet {
    Cost value = 0;       // 選んだ頂点の重みの和（重みなしなら個数）
    std::vector<int> vs;  // 選んだ頂点（昇順）
};
struct TreeMatching {
    Cost value = 0;                         // 選んだ辺の重みの和（重みなしなら本数）
    std::vector<int> es;                    // 選んだ辺の番号（昇順）
    std::vector<std::pair<int, int>> pairs; // 選んだ辺の両端 (u, v)（u < v。es の順）
};
namespace gg_internal {
inline std::vector<Cost> tree_weights(int n, const std::vector<Cost>& w) { return w.empty() ? std::vector<Cost>(n, 1) : w; }
inline TreeVertexSet tree_set_from(const std::vector<char>& in, const std::vector<Cost>& w) {
    TreeVertexSet r;
    for (int v = 0; v < (int)in.size(); ++v)
        if (in[v]) r.vs.push_back(v), r.value += w[v];
    return r;
}
}  // namespace gg_internal
// 最大（重み）独立集合: どの 2 つも隣り合わない頂点の集合で、重みの和が最大（負の重みの頂点は選ばない）
inline TreeVertexSet tree_max_independent_set(const Graph& G, const std::vector<Cost>& w_in = {}) {
    const int n = G.size();
    const std::vector<Cost> w = gg_internal::tree_weights(n, w_in);
    const RootedTree T = RootedTree::forest(G);
    std::vector<Cost> in(n, 0), out(n, 0);  // in[v] = v を選ぶときの部分木の最大、out[v] = 選ばないとき
    for (int k = n - 1; k >= 0; --k) {
        const int v = T.order[k];
        in[v] += std::max<Cost>(w[v], 0);
        if (T.par[v] >= 0) {
            const int p = T.par[v];
            in[p] += out[v], out[p] += std::max(in[v], out[v]);
        }
    }
    std::vector<char> take(n, 0);
    for (int v : T.order) {  // 根から: 親を選んだら選ばない、選ばなかったら良いほう
        const bool parent_taken = T.par[v] >= 0 && take[T.par[v]];
        take[v] = !parent_taken && w[v] > 0 && in[v] >= out[v];
    }
    return gg_internal::tree_set_from(take, w);
}
// 最小（重み）点被覆: すべての辺の端点のどちらかを含む頂点の集合で、重みの和が最小（重みは 0 以上）。
// 最大独立集合の補集合（重みなしなら 大きさ = N - 最大独立集合 = 最大マッチング）
inline TreeVertexSet tree_min_vertex_cover(const Graph& G, const std::vector<Cost>& w_in = {}) {
    const int n = G.size();
    const std::vector<Cost> w = gg_internal::tree_weights(n, w_in);
    const TreeVertexSet is = tree_max_independent_set(G, w);
    std::vector<char> in(n, 1);
    for (int v : is.vs) in[v] = 0;
    for (int v = 0; v < n; ++v)  // 重み 0 の頂点は独立集合に入れていないので、要らなければ外す
        if (in[v] && w[v] == 0) {
            bool need = false;
            for (const Edge& e : G[v])
                if (e.to != v && !in[e.to]) need = true;
            if (!need) in[v] = 0;
        }
    for (int v = 0; v < n; ++v)  // 自己ループは必ずその頂点で覆う
        for (const Edge& e : G[v])
            if (e.to == v) in[v] = 1;
    return gg_internal::tree_set_from(in, w);
}
// 最小（重み）支配集合: どの頂点も「自分か隣」が集合に入っている頂点の集合で、重みの和が最小（重みは 0 以上）
inline TreeVertexSet tree_min_dominating_set(const Graph& G, const std::vector<Cost>& w_in = {}) {
    const int n = G.size();
    const std::vector<Cost> w = gg_internal::tree_weights(n, w_in);
    const RootedTree T = RootedTree::forest(G);
    const Cost INF = DIST_INF;
    // a[v] = v を選ぶ、b[v] = v は選ばず子のどれかが選ばれている、c[v] = v は選ばず子も選ばれていない（親が選ぶ）
    std::vector<Cost> a(n), b(n), c(n);
    std::vector<int> best_child(n, -1);  // b[v] で「選ぶ子」にした子
    for (int k = n - 1; k >= 0; --k) {
        const int v = T.order[k];
        Cost sa = w[v], sb = 0, sc = 0, extra = INF;
        int bc = -1;
        for (int ch : T.children(v)) {
            const Cost m3 = std::min(a[ch], std::min(b[ch], c[ch]));
            const Cost m2 = std::min(a[ch], b[ch]);
            sa += m3, sb += m2, sc = (sc >= INF || b[ch] >= INF) ? INF : sc + b[ch];
            if (a[ch] - m2 < extra) extra = a[ch] - m2, bc = ch;
        }
        a[v] = sa, b[v] = (bc < 0 ? INF : sb + extra), c[v] = sc, best_child[v] = bc;
    }
    std::vector<int> st(n, -1);  // 0 = a, 1 = b, 2 = c
    std::vector<char> take(n, 0);
    for (int v : T.order) {
        const int p = T.par[v];
        if (p < 0) {
            st[v] = (a[v] <= b[v]) ? 0 : 1;
        } else if (st[p] == 0) {
            st[v] = (a[v] <= b[v] && a[v] <= c[v]) ? 0 : (b[v] <= c[v] ? 1 : 2);
        } else if (st[p] == 1 && best_child[p] == v) {
            st[v] = 0;
        } else if (st[p] == 1) {
            st[v] = (a[v] <= b[v]) ? 0 : 1;
        } else {  // st[p] == 2: 子はすべて b（選ばず、自分の子に支配される）
            st[v] = 1;
        }
        take[v] = (st[v] == 0);
    }
    return gg_internal::tree_set_from(take, w);
}
// 最大マッチング（重みなし。辺の数が最大）: 葉の側から貪欲に、親と組める頂点は組む
inline TreeMatching tree_max_matching(const Graph& G) {
    const int n = G.size();
    const RootedTree T = RootedTree::forest(G);
    std::vector<char> used(n, 0);
    TreeMatching r;
    for (int k = n - 1; k >= 0; --k) {
        const int v = T.order[k], p = T.par[v];
        if (p >= 0 && !used[v] && !used[p]) used[v] = used[p] = 1, r.es.push_back(T.par_edge[v]);
    }
    std::sort(r.es.begin(), r.es.end());
    for (int id : r.es) r.pairs.push_back({std::min(G.es[id].from, G.es[id].to), std::max(G.es[id].from, G.es[id].to)});
    r.value = (Cost)r.es.size();
    return r;
}
// 最大重みマッチング（辺の重み e.w の和が最大。負の重みの辺は使わない）
inline TreeMatching tree_max_weight_matching(const Graph& G) {
    const int n = G.size();
    const RootedTree T = RootedTree::forest(G);
    std::vector<Cost> f(n, 0), g(n, 0);  // f[v] = 部分木の最大（v は自由でも組んでもよい）、g[v] = v が子と組まないときの最大
    std::vector<int> mate(n, -1);         // f[v] で v と組む子（なければ -1）
    for (int k = n - 1; k >= 0; --k) {
        const int v = T.order[k];
        Cost sum = 0, gain = 0;
        int mc = -1;
        for (int ch : T.children(v)) sum += f[ch];
        for (int ch : T.children(v)) {
            const Cost d = T.par_w[ch] + g[ch] - f[ch];  // ch と組むと増える分
            if (d > gain) gain = d, mc = ch;
        }
        g[v] = sum, f[v] = sum + gain, mate[v] = mc;
    }
    std::vector<char> matched_up(n, 0);  // 親と組むことが決まっている
    TreeMatching r;
    for (int v : T.order) {
        if (matched_up[v]) continue;  // 親と組むので、子とは組まない
        if (mate[v] >= 0) matched_up[mate[v]] = 1, r.es.push_back(T.par_edge[mate[v]]);
    }
    std::sort(r.es.begin(), r.es.end());
    for (int id : r.es) {
        r.pairs.push_back({std::min(G.es[id].from, G.es[id].to), std::max(G.es[id].from, G.es[id].to)});
        r.value += G.es[id].w;
    }
    return r;
}
// 独立集合（空集合も含む）の個数 mod（EDPC P: 隣り合う 2 頂点をともに黒にしない塗り方の数）。森なら成分ごとの積
inline long long tree_count_independent_sets(const Graph& G, long long mod) {
    const int n = G.size();
    const RootedTree T = RootedTree::forest(G);
    std::vector<long long> in(n, 1), out(n, 1);
    for (int k = n - 1; k >= 0; --k) {
        const int v = T.order[k], p = T.par[v];
        if (p >= 0) in[p] = in[p] * out[v] % mod, out[p] = out[p] * ((in[v] + out[v]) % mod) % mod;
    }
    long long r = 1 % mod;
    for (int root : T.roots) r = r * ((in[root] + out[root]) % mod) % mod;
    return r;
}
#endif  // GG_TREE_CLASSIC
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "tree_classic";
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

namespace tc_test {

static Graph random_forest(std::mt19937_64& rng, int n, int shape) {
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    Graph G(n);
    for (int i = 1; i < n; ++i) {
        if (shape == 3 && rng() % 4 == 0) continue;  // 森
        const int p = shape == 1 ? i - 1 : shape == 2 ? 0 : rnd(rng, 0, i - 1);
        int a = perm[i], b = perm[p];
        if (rng() & 1) std::swap(a, b);
        G.add_undirected(a, b, rnd(rng, -3, 9));
    }
    return G;
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 11);
        const Graph G = random_forest(rng, n, rnd(rng, 0, 3));
        const int m = G.edge_count();
        std::vector<Cost> wa(n), wp(n);  // wa: 負もある、wp: 0 以上
        for (int v = 0; v < n; ++v) wa[v] = rnd(rng, -3, 9), wp[v] = rnd(rng, 0, 9);
        const bool unit = (it % 3 == 0);
        const std::vector<Cost> WA = unit ? std::vector<Cost>{} : wa, WP = unit ? std::vector<Cost>{} : wp;
        const std::vector<Cost> wa1 = unit ? std::vector<Cost>(n, 1) : wa, wp1 = unit ? std::vector<Cost>(n, 1) : wp;
        // 愚直: 頂点の部分集合を全部
        Cost best_is = 0, best_vc = DIST_INF, best_ds = DIST_INF;
        long long cnt_is = 0;
        for (int s = 0; s < (1 << n); ++s) {
            bool indep = true, cover = true;
            for (const RawEdge& e : G.es) {
                const bool a = (s >> e.from) & 1, b = (s >> e.to) & 1;
                if (a && b) indep = false;
                if (!a && !b) cover = false;
            }
            bool dom = true;
            for (int v = 0; v < n; ++v) {
                bool d = (s >> v) & 1;
                for (const Edge& e : G[v]) d = d || ((s >> e.to) & 1);
                if (!d) dom = false;
            }
            Cost sa = 0, sp = 0;
            for (int v = 0; v < n; ++v)
                if ((s >> v) & 1) sa += wa1[v], sp += wp1[v];
            if (indep) best_is = std::max(best_is, sa), ++cnt_is;
            if (cover) best_vc = std::min(best_vc, sp);
            if (dom) best_ds = std::min(best_ds, sp);
        }
        auto valid_set = [&](const TreeVertexSet& r, const std::vector<Cost>& w) {
            Cost s = 0;
            for (size_t i = 0; i < r.vs.size(); ++i) {
                CHECK(r.vs[i] >= 0 && r.vs[i] < n);
                if (i) CHECK(r.vs[i - 1] < r.vs[i]);
                s += w[r.vs[i]];
            }
            CHECK(s == r.value);
            std::vector<char> in(n, 0);
            for (int v : r.vs) in[v] = 1;
            return in;
        };
        {
            const TreeVertexSet r = tree_max_independent_set(G, WA);
            CHECK(r.value == best_is);
            const std::vector<char> in = valid_set(r, wa1);
            for (const RawEdge& e : G.es) CHECK(!(in[e.from] && in[e.to]));
        }
        {
            const TreeVertexSet r = tree_min_vertex_cover(G, WP);
            CHECK(r.value == best_vc);
            const std::vector<char> in = valid_set(r, wp1);
            for (const RawEdge& e : G.es) CHECK(in[e.from] || in[e.to]);
        }
        {
            const TreeVertexSet r = tree_min_dominating_set(G, WP);
            CHECK(r.value == best_ds);
            const std::vector<char> in = valid_set(r, wp1);
            for (int v = 0; v < n; ++v) {
                bool d = in[v];
                for (const Edge& e : G[v]) d = d || in[e.to];
                CHECK(d);
            }
        }
        CHECK(tree_count_independent_sets(G, 1000000007) == cnt_is % 1000000007);
        CHECK(tree_count_independent_sets(G, 7) == cnt_is % 7);
        // マッチング: 辺の部分集合を全部
        int best_m = 0;
        Cost best_wm = 0;
        for (int s = 0; s < (1 << m); ++s) {
            std::vector<int> deg(n, 0);
            bool ok = true;
            int c = 0;
            Cost w = 0;
            for (int id = 0; id < m; ++id)
                if ((s >> id) & 1) {
                    if (++deg[G.es[id].from] > 1 || ++deg[G.es[id].to] > 1) ok = false;
                    ++c, w += G.es[id].w;
                }
            if (ok) best_m = std::max(best_m, c), best_wm = std::max(best_wm, w);
        }
        auto valid_matching = [&](const TreeMatching& r, bool weighted) {
            std::vector<int> deg(n, 0);
            Cost w = 0;
            CHECK(r.es.size() == r.pairs.size());
            for (size_t i = 0; i < r.es.size(); ++i) {
                const RawEdge& e = G.es[r.es[i]];
                if (i) CHECK(r.es[i - 1] < r.es[i]);
                CHECK(++deg[e.from] == 1 && ++deg[e.to] == 1);
                CHECK(r.pairs[i] == std::make_pair(std::min(e.from, e.to), std::max(e.from, e.to)));
                w += weighted ? e.w : 1;
            }
            CHECK(w == r.value);
        };
        {
            const TreeMatching r = tree_max_matching(G);
            CHECK(r.value == best_m);
            valid_matching(r, false);
        }
        {
            const TreeMatching r = tree_max_weight_matching(G);
            CHECK(r.value == best_wm);
            valid_matching(r, true);
        }
        if (unit) {  // 重みなしの関係: 最小点被覆 = 最大マッチング = N - 最大独立集合
            CHECK(tree_min_vertex_cover(G).value == tree_max_matching(G).value);
            CHECK(tree_min_vertex_cover(G).value == n - tree_max_independent_set(G).value);
        }
    }
    std::cout << "[tree_classic] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[tree_classic] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    const int n = 200000;
    std::vector<int> us(n - 1), vs(n - 1);
    std::vector<long long> ws(n - 1);
    for (int i = 1; i < n; ++i) us[i - 1] = i, vs[i - 1] = (int)(rng() % i), ws[i - 1] = (long long)(rng() % 1000);
    const Graph G = graph_from_edges(n, us, vs, ws, UNDIRECTED, 0);
    std::vector<Cost> w(n);
    for (auto& x : w) x = (Cost)(rng() % 1000);
    auto t0 = std::chrono::steady_clock::now();
    const Cost a = tree_max_independent_set(G, w).value;
    const long long t1 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const Cost b = tree_min_dominating_set(G, w).value;
    const long long t2 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const Cost c = tree_max_weight_matching(G).value;
    const long long t3 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const long long d = tree_count_independent_sets(G, 998244353);
    std::cout << "[tree_classic] benchmark: N = " << n << "  最大独立集合 " << t1 << " ms, 最小支配集合 " << t2 << " ms, 最大重みマッチング " << t3 << " ms, 独立集合の数 "
              << ms_since(t0) << " ms  (check " << (a + b + c + d) % 1000 << ")\n";
}

}  // namespace tc_test


namespace tc_demo {
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
    banner("木の上の有名な最適化 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"tree_max_independent_set(G, w = {})", "TreeVertexSet", "隣り合わない頂点で重み最大 {value, vs}"},
        {"tree_min_vertex_cover(G, w = {})", "TreeVertexSet", "全部の辺を覆う頂点で重み最小"},
        {"tree_min_dominating_set(G, w = {})", "TreeVertexSet", "全部の頂点が自分か隣に入る、重み最小"},
        {"tree_max_matching(G)", "TreeMatching", "端点を共有しない辺の最大本数 {value, es, pairs}"},
        {"tree_max_weight_matching(G)", "TreeMatching", "辺の重みの和が最大"},
        {"tree_count_independent_sets(G, mod)", "long long", "独立集合の数（空集合も含む）"},
    };
    cheat_sheet("API 一覧", rows, 40, 18);

    head("デモで使う木: 0-1, 1-2, 1-3, 3-4, 3-5（重み 4, 1, 2, 5, 3）");
    RUN(Graph G(6));
    RUN(G.add_undirected(0, 1, 4); G.add_undirected(1, 2, 1); G.add_undirected(1, 3, 2); G.add_undirected(3, 4, 5); G.add_undirected(3, 5, 3));

    api("【1】 最大独立集合（重みなし = 個数）", "TreeVertexSet", "葉を選んでいくのが最適", "O(N)");
    RUN(auto is = tree_max_independent_set(G));
    SHOW(is.value, "");
    SHOW(is.vs, "");
    RUN(std::vector<Cost> w = {1, 10, 1, 1, 1, 1});
    SHOW(tree_max_independent_set(G, w).vs, "頂点の重みつき（1 が重い）");

    api("【2】 最小点被覆・最小支配集合", "TreeVertexSet", "重みは 0 以上", "O(N)");
    SHOW(tree_min_vertex_cover(G).vs, "大きさ = 最大マッチング");
    SHOW(tree_min_dominating_set(G).vs, "1 と 3 で全部が支配される");
    SHOW(tree_min_dominating_set(G, w).vs, "1 が重いと別の選び方");

    api("【3】 最大マッチング", "TreeMatching", "重みなし: 葉の側から貪欲。重みつき: 辺の重み e.w の和", "O(N)");
    RUN(auto mt = tree_max_matching(G));
    SHOW(mt.value, "");
    SHOW(mt.pairs, "");
    RUN(auto wm = tree_max_weight_matching(G));
    SHOW(wm.value, "4 + 5 = 9");
    SHOW(wm.es, "辺の番号");

    api("【4】 独立集合の数（EDPC P: 隣り合う 2 頂点をともに黒にしない塗り方）", "long long", "mod で数える", "O(N)");
    SHOW(tree_count_independent_sets(G, 1000000007), "");

    head("注意");
    note("G は木か森（閉路があると正しくない）。森は成分ごとに解いた結果をまとめて返す。");
    note("ほかの木DP は 05_rooted_tree の型（order を逆順に回す）で書く。全頂点を根にしたいなら 05_rerooting。");
}

}  // namespace tc_demo

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
    tc_test::random_test();
    tc_test::benchmark();
    tc_demo::demo();
    return 0;
}
