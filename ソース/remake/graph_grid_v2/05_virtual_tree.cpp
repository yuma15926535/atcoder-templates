// =============================================================================
//  虚树（auxiliary tree・圧縮した木）  virtual_tree / VirtualTree            --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 木
//  依存       : 00_graph_base, 05_rooted_tree, 05_lca（同梱済み）
//  頂点番号   : 0-indexed。G は木（森でもよい: LCA L(RootedTree::forest(G)) で作る）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  指定した K 頂点と、それらの LCA だけを残した木（頂点数 2K - 1 以下）を O(K log N) で作る。
//  祖先・子孫の関係と距離は元の木と同じ。クエリごとに頂点の集合が与えられ、その大きさの合計が小さい問題で、
//  圧縮した木の上で木DP をする（ABC340 G・ABC359 G など）。
//    - vs は元の木の行きがけ順（HLD 順）に並ぶので、i を大きい順に回すと子が親より先（木DP の順）
//    - 指定した頂点をすべて含む最小の部分木の辺の本数 = Σ dist[i]（ABC268 D 型の数え方）
//
//  機能一覧                                                                       計算量
//    VirtualTree T = virtual_tree(L, vs)    vs（重複可）とその LCA の木                   O(K log N)
//      T.size() / T.vs[i]（元の頂点）/ T.par[i]（親の添字、根は -1）/ T.roots
//      T.dist[i] / T.wdist[i]（親までの辺の本数 / 重みの和）/ T.given[i]（指定した頂点か）
//      T.index(v, L)（元の頂点 v の添字。なければ -1）                                      O(log K)
//      T.graph()（添字を頂点にした Graph。辺の重み = wdist）                                O(K)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤・根付き木・LCA も入っている）
//    vt_test::random_test() / vt_test::benchmark() / vt_demo::demo()
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
#ifndef GG_LCA
#define GG_LCA
// ---- LCA（グラフ・グリッド第二世代）  重軽分解（HLD）で求める。クエリ O(log N)・メモリ O(N)・非再帰 ----
// 森にも対応（別々の木の頂点どうしは lca = -1、dist = -1、wdist = DIST_INF）
struct LCA {
    int n = 0;
    std::vector<int> par, depth, sz, root_of;  // RootedTree からコピーしたもの
    std::vector<Cost> wdepth;
    std::vector<int> heavy;      // 部分木が最大の子（葉は -1）
    std::vector<int> head;       // v を含む heavy な道の先頭（根に近い側）
    std::vector<int> hld_pos;    // HLD 順での位置（heavy な道は連続、部分木も連続 = [hld_pos[v], hld_pos[v] + sz[v])）
    std::vector<int> hld_order;  // hld_order[hld_pos[v]] = v

    LCA() {}
    explicit LCA(const RootedTree& T) { build(T); }
    explicit LCA(const Graph& G, int root = 0) { build(RootedTree(G, root)); }
    void build(const RootedTree& T) {
        n = T.n, par = T.par, depth = T.depth, sz = T.sz, root_of = T.root_of, wdepth = T.wdepth;
        heavy.assign(n, -1), head.assign(n, -1), hld_pos.assign(n, -1), hld_order.clear();
        hld_order.reserve(T.order.size());
        for (int v : T.order)
            for (int c : T.children(v))
                if (heavy[v] == -1 || sz[c] > sz[heavy[v]]) heavy[v] = c;
        std::vector<int> st;
        for (int r : T.roots) {
            head[r] = r, st.push_back(r);
            while (!st.empty()) {
                const int v = st.back();
                st.pop_back();
                hld_pos[v] = (int)hld_order.size(), hld_order.push_back(v);
                const IntSpan cs = T.children(v);
                for (int k = cs.size() - 1; k >= 0; --k)  // 軽い子を先に積み、heavy を最後に積む（次に取り出す）
                    if (cs[k] != heavy[v]) head[cs[k]] = cs[k], st.push_back(cs[k]);
                if (heavy[v] != -1) head[heavy[v]] = head[v], st.push_back(heavy[v]);
            }
        }
    }
    bool same_tree(int u, int v) const { return root_of[u] != -1 && root_of[u] == root_of[v]; }
    int lca(int u, int v) const {  // 最小共通祖先（別の木・未到達なら -1）
        if (!same_tree(u, v)) return -1;
        while (head[u] != head[v]) {
            if (depth[head[u]] > depth[head[v]]) u = par[head[u]];
            else v = par[head[v]];
        }
        return depth[u] < depth[v] ? u : v;
    }
    int lca(const std::vector<int>& vs) const {  // 複数の頂点の LCA
        int a = vs.empty() ? -1 : vs[0];
        for (size_t i = 1; i < vs.size() && a != -1; ++i) a = lca(a, vs[i]);
        return a;
    }
    int dist(int u, int v) const {  // 辺の本数（別の木なら -1）
        const int a = lca(u, v);
        return a < 0 ? -1 : depth[u] + depth[v] - 2 * depth[a];
    }
    Cost wdist(int u, int v) const {  // 重みの和（別の木なら DIST_INF）
        const int a = lca(u, v);
        return a < 0 ? DIST_INF : wdepth[u] + wdepth[v] - 2 * wdepth[a];
    }
    int kth_ancestor(int v, int k) const {  // k 個上の祖先（k = 0 なら v。根より上なら -1）
        if (k < 0 || depth[v] < 0 || k > depth[v]) return -1;
        const int d = depth[v] - k;
        while (depth[head[v]] > d) v = par[head[v]];
        return hld_order[hld_pos[v] - (depth[v] - d)];
    }
    int jump(int u, int v, int k) const {  // u → v のパスで u から k 番目の頂点（k = 0 なら u。パスより長ければ -1）
        const int a = lca(u, v);
        if (a < 0 || k < 0) return -1;
        const int du = depth[u] - depth[a], dv = depth[v] - depth[a];
        if (k <= du) return kth_ancestor(u, k);
        if (k <= du + dv) return kth_ancestor(v, du + dv - k);
        return -1;
    }
    std::vector<int> path(int u, int v) const {  // u → v の頂点列（別の木なら空）
        const int a = lca(u, v);
        std::vector<int> p, q;
        if (a < 0) return p;
        for (int x = u; x != a; x = par[x]) p.push_back(x);
        p.push_back(a);
        for (int x = v; x != a; x = par[x]) q.push_back(x);
        p.insert(p.end(), q.rbegin(), q.rend());
        return p;
    }
    bool on_path(int u, int v, int a) const {  // a が u - v のパス上にあるか
        if (!same_tree(u, v) || !same_tree(u, a)) return false;
        return dist(u, a) + dist(a, v) == dist(u, v);
    }
    bool is_ancestor(int u, int v) const {  // u が v の祖先か（u == v も true）
        return same_tree(u, v) && hld_pos[u] <= hld_pos[v] && hld_pos[v] < hld_pos[u] + sz[u];
    }
    std::pair<int, int> subtree_range(int v) const { return {hld_pos[v], hld_pos[v] + sz[v]}; }  // HLD 順での部分木 [l, r)
    // u → v のパスを HLD 順の区間に分ける（パス上の値をセグ木で扱うとき）。
    // 返り値の各要素 (a, b) は「HLD 順の位置 a から b まで（両端を含む）を、この順にたどる」。a > b なら位置が減る向き。
    // 演算が可換なら [min(a, b), max(a, b) + 1) をそのまま使えばよい。
    // edge = true なら辺に載せた値用（各辺の値は子の側の頂点に置く）で、LCA の頂点を含めない
    std::vector<std::pair<int, int>> path_segments(int u, int v, bool edge = false) const {
        std::vector<std::pair<int, int>> up, down;
        if (!same_tree(u, v)) return up;
        while (head[u] != head[v]) {
            if (depth[head[u]] >= depth[head[v]]) up.push_back({hld_pos[u], hld_pos[head[u]]}), u = par[head[u]];
            else down.push_back({hld_pos[head[v]], hld_pos[v]}), v = par[head[v]];
        }
        const int e = edge ? 1 : 0;
        if (depth[u] >= depth[v]) {
            if (hld_pos[u] >= hld_pos[v] + e) up.push_back({hld_pos[u], hld_pos[v] + e});
        } else {
            if (hld_pos[v] >= hld_pos[u] + e) down.push_back({hld_pos[u] + e, hld_pos[v]});
        }
        up.insert(up.end(), down.rbegin(), down.rend());
        return up;
    }
};
#endif  // GG_LCA
#ifndef GG_VIRTUAL_TREE
#define GG_VIRTUAL_TREE
// ---- 虚树（auxiliary tree・圧縮した木）（グラフ・グリッド第二世代）  指定した K 頂点とその LCA だけの木を O(K log N) で作る ----
// 元の木の祖先・子孫の関係を保ったまま、指定した頂点どうしをつなぐのに要る頂点（指定した頂点と、それらの LCA）だけを残す。
// 頂点数は 2K - 1 以下。クエリごとに K 頂点が与えられ、K の合計が小さい問題（ABC340 G・ABC359 G など）で、木DP を K の大きさで回す。
// 05_lca の LCA を使う（根は LCA を作ったときの根。森なら木ごとに別の根になる）
struct VirtualTree {
    std::vector<int> vs;     // 頂点（元の番号）。元の木の行きがけ順（HLD 順）に並ぶ。par[i] < i
    std::vector<int> par;    // par[i] = 親の添字（根は -1）
    std::vector<int> dist;   // dist[i] = 親までの元の木での辺の本数（根は 0）
    std::vector<Cost> wdist; // wdist[i] = 親までの元の木での重みの和（根は 0）
    std::vector<int> given;  // given[i] = 1: vs[i] は指定した頂点 / 0: LCA として足した頂点
    std::vector<int> roots;  // 根の添字（森のときは木ごとに 1 つ）
    std::vector<int> pos_;   // pos_[i] = vs[i] の HLD 順の位置（index の二分探索用）
    int size() const { return (int)vs.size(); }
    // 元の頂点 v の添字（なければ -1）  O(log K)
    int index(int v, const LCA& L) const {
        const int p = L.hld_pos[v];
        int lo = 0, hi = size();
        while (lo < hi) {
            const int mid = (lo + hi) / 2;
            if (pos_[mid] < p) lo = mid + 1;
            else hi = mid;
        }
        return (lo < size() && pos_[lo] == p) ? lo : -1;
    }
    // 圧縮した木を Graph（頂点 = 添字、辺 i - 1 本目 = par[i] - i、重み = wdist[i]）にする。無向
    Graph graph() const {
        Graph T(size());
        for (int i = 0; i < size(); ++i)
            if (par[i] >= 0) T.add_undirected(par[i], i, wdist[i]);
        return T;
    }
};
// vs の頂点（重複してよい。L の根から届かない頂点は無視する）の虚树
inline VirtualTree virtual_tree(const LCA& L, const std::vector<int>& vs_in) {
    std::vector<std::pair<int, int>> a;  // (HLD 順の位置, 頂点)
    a.reserve(vs_in.size() * 2);
    for (int v : vs_in)
        if (L.root_of[v] != -1) a.push_back({L.hld_pos[v], v});
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
    const int k = (int)a.size();
    std::vector<int> given_v;
    given_v.reserve(k);
    for (const auto& p : a) given_v.push_back(p.second);
    for (int i = 0; i + 1 < k; ++i) {  // 行きがけ順で隣どうしの LCA を足せば、必要な LCA はすべてそろう
        const int c = L.lca(a[i].second, a[i + 1].second);
        if (c != -1) a.push_back({L.hld_pos[c], c});
    }
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
    VirtualTree T;
    const int m = (int)a.size();
    T.vs.resize(m), T.pos_.resize(m), T.par.assign(m, -1), T.dist.assign(m, 0), T.wdist.assign(m, 0), T.given.assign(m, 0);
    for (int i = 0; i < m; ++i) T.pos_[i] = a[i].first, T.vs[i] = a[i].second;
    {
        int j = 0;  // given_v も HLD 順に並んでいる
        for (int i = 0; i < m && j < k; ++i)
            if (T.vs[i] == given_v[j]) T.given[i] = 1, ++j;
    }
    std::vector<int> st;  // 今の頂点の祖先の添字（根に近い順）
    for (int i = 0; i < m; ++i) {
        const int v = T.vs[i];
        while (!st.empty() && !L.is_ancestor(T.vs[st.back()], v)) st.pop_back();
        if (st.empty()) {
            T.roots.push_back(i);
        } else {
            const int p = T.vs[st.back()];
            T.par[i] = st.back(), T.dist[i] = L.depth[v] - L.depth[p], T.wdist[i] = L.wdepth[v] - L.wdepth[p];
        }
        st.push_back(i);
    }
    return T;
}
#endif  // GG_VIRTUAL_TREE
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "virtual_tree";
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

namespace vt_test {

static Graph random_forest(std::mt19937_64& rng, int n, int shape) {
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    Graph G(n);
    for (int i = 1; i < n; ++i) {
        if (shape == 3 && rng() % 5 == 0) continue;  // 森にする
        const int p = shape == 1 ? i - 1 : shape == 2 ? 0 : rnd(rng, 0, i - 1);
        int a = perm[i], b = perm[p];
        if (rng() & 1) std::swap(a, b);
        G.add_undirected(a, b, rnd(rng, 0, 9));
    }
    return G;
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 30);
        const Graph G = random_forest(rng, n, rnd(rng, 0, 3));
        const RootedTree RT = (it & 1) ? RootedTree::forest(G) : RootedTree(G, rnd(rng, 0, n - 1));
        const LCA L(RT);
        std::vector<int> q;
        const int k = rnd(rng, 0, std::min(n, 8) + 2);
        for (int t = 0; t < k; ++t) q.push_back(rnd(rng, 0, n - 1));
        const VirtualTree T = virtual_tree(L, q);
        // 愚直: 頂点の集合 = 届く指定頂点と、同じ木の 2 つの LCA
        std::vector<char> in(n, 0), giv(n, 0);
        for (int v : q)
            if (L.root_of[v] != -1) in[v] = 1, giv[v] = 1;
        for (int a : q)
            for (int b : q) {
                const int c = L.lca(a, b);
                if (c != -1) in[c] = 1;
            }
        std::vector<int> exp;
        for (int v : L.hld_order)
            if (in[v]) exp.push_back(v);
        CHECK(T.vs == exp);
        CHECK(T.size() == (int)exp.size());
        CHECK((int)T.par.size() == T.size() && (int)T.dist.size() == T.size() && (int)T.wdist.size() == T.size() && (int)T.given.size() == T.size());
        int nroots = 0;
        for (int i = 0; i < T.size(); ++i) {
            const int v = T.vs[i];
            CHECK(T.given[i] == giv[v]);
            CHECK(T.index(v, L) == i);
            // 親 = 集合の中で最も近い真の祖先
            int p = RT.par[v];
            while (p != -1 && !in[p]) p = RT.par[p];
            if (p == -1) {
                CHECK(T.par[i] == -1 && T.dist[i] == 0 && T.wdist[i] == 0);
                ++nroots;
            } else {
                CHECK(T.par[i] >= 0 && T.par[i] < i && T.vs[T.par[i]] == p);
                CHECK(T.dist[i] == L.dist(v, p) && T.wdist[i] == L.wdist(v, p));
            }
        }
        CHECK((int)T.roots.size() == nroots);
        for (int r : T.roots) CHECK(T.par[r] == -1);
        for (int v = 0; v < n; ++v)
            if (!in[v]) CHECK(T.index(v, L) == -1);
        if (k > 0) CHECK(T.size() <= 2 * k);
        // graph(): 辺の本数 = 頂点数 - 根の数、重み = wdist
        const Graph TG = T.graph();
        CHECK(TG.size() == T.size() && TG.edge_count() == T.size() - nroots);
        for (const RawEdge& e : TG.es) CHECK(T.par[e.to] == e.from && T.wdist[e.to] == e.w);
        // 指定した頂点をすべて含む最小の部分木の辺の本数（1 つの木の中だけで指定したとき）
        if (it % 2 == 1 && k > 0 && nroots == 1) {
            std::vector<int> mark(n, 0);
            for (int v : q) mark[v] = 1;
            const RootedTree R2(G, T.vs[T.roots[0]]);
            const std::vector<int> c = R2.subtree_sum(mark);
            int cnt = 0;
            for (int v = 0; v < n; ++v) cnt += (R2.reached(v) && c[v] > 0);
            int sum = 0;
            for (int i = 0; i < T.size(); ++i) sum += T.dist[i];
            CHECK(sum + 1 == cnt);
        }
    }
    std::cout << "[virtual_tree] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[virtual_tree] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    const int n = 200000;
    Graph G(n);
    for (int i = 1; i < n; ++i) G.add_undirected(i, (int)(rng() % i));
    const LCA L(G, 0);
    for (int k : {2, 200, 200000}) {
        const int q = 200000 / k;
        auto t0 = std::chrono::steady_clock::now();
        long long chk = 0;
        for (int t = 0; t < q; ++t) {
            std::vector<int> vs(k);
            for (int& v : vs) v = (int)(rng() % n);
            const VirtualTree T = virtual_tree(L, vs);
            chk += T.size();
        }
        std::cout << "[virtual_tree] benchmark: N = " << n << "  K = " << k << " を " << q << " 回（K の合計 200000） " << ms_since(t0) << " ms  (check " << chk % 1000
                  << ")\n";
    }
}

}  // namespace vt_test


namespace vt_demo {
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
    banner("虚树（auxiliary tree）使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"VirtualTree T = virtual_tree(L, vs);", "VirtualTree", "vs とその LCA だけの木（L は 05_lca の LCA）"},
        {"T.size() / T.vs[i]", "int", "頂点数 / 添字 i の元の頂点（行きがけ順）"},
        {"T.par[i] / T.roots", "int", "親の添字（根は -1）/ 根の添字"},
        {"T.dist[i] / T.wdist[i]", "", "親までの辺の本数 / 重み"},
        {"T.given[i]", "int", "指定した頂点なら 1（LCA として足しただけなら 0）"},
        {"T.index(v, L)", "int", "元の頂点 v の添字（なければ -1）"},
        {"T.graph()", "Graph", "添字を頂点にしたグラフ（重み = wdist）"},
    };
    cheat_sheet("API 一覧", rows, 40, 16);

    head("デモで使う木（根 0）: 0-1, 0-2, 1-3, 1-4, 2-5, 4-6, 4-7");
    RUN(Graph G(8));
    RUN(for (auto [a, b] : std::vector<std::pair<int, int>>{{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}, {4, 6}, {4, 7}}) G.add_undirected(a, b));
    RUN(LCA L(G, 0));

    api("【1】 virtual_tree(L, {3, 6, 7, 5})", "VirtualTree", "3 と 6・7 の LCA は 1、6 と 7 の LCA は 4、全体の LCA は 0", "O(K log N)");
    RUN(VirtualTree T = virtual_tree(L, {3, 6, 7, 5}));
    SHOW(T.vs, "元の頂点（行きがけ順）");
    SHOW(T.par, "親の添字");
    SHOW(T.dist, "親までの辺の本数");
    SHOW(T.given, "指定した頂点か");
    SHOW(T.index(4, L), "頂点 4 の添字");

    api("【2】 指定した頂点をすべて結ぶ最小の部分木の辺の本数（ABC268 D 型）", "int", "Σ dist[i]（頂点数はそれ + 1）", "O(K)");
    code("int edges = 0; for (int i = 0; i < T.size(); ++i) edges += T.dist[i];");
    int edges = 0;
    for (int i = 0; i < T.size(); ++i) edges += T.dist[i];
    SHOW(edges, "0-1, 1-3, 1-4, 4-6, 4-7, 0-2, 2-5");

    api("【3】 圧縮した木の上の木DP: 添字を大きい順に回すと子が先", "vector<long long>", "例: 部分木に入る指定頂点の数 → 辺 i - par[i] を通るペアの数", "O(K)");
    code("std::vector<long long> cnt(T.size());");
    code("for (int i = T.size() - 1; i >= 0; --i) { cnt[i] += T.given[i]; if (T.par[i] >= 0) cnt[T.par[i]] += cnt[i]; }");
    code("long long sum = 0; for (int i = 0; i < T.size(); ++i) if (T.par[i] >= 0) sum += T.wdist[i] * cnt[i] * (K - cnt[i]);   // 指定頂点の全ペアの距離の和");
    std::vector<long long> cnt(T.size());
    for (int i = T.size() - 1; i >= 0; --i) {
        cnt[i] += T.given[i];
        if (T.par[i] >= 0) cnt[T.par[i]] += cnt[i];
    }
    long long sum = 0;
    for (int i = 0; i < T.size(); ++i)
        if (T.par[i] >= 0) sum += T.wdist[i] * cnt[i] * (4 - cnt[i]);
    SHOW(cnt, "");
    SHOW(sum, "d(3,6) + d(3,7) + d(3,5) + d(6,7) + d(6,5) + d(7,5) = 3 + 3 + 4 + 2 + 5 + 5");

    head("注意");
    note("L を作ったときの根から届かない頂点（森で別の木・RootedTree(G, r) で r と別の成分）は無視する。森なら LCA L(RootedTree::forest(G))。");
    note("クエリごとに作り直す（O(K log K + K log N)）。元の頂点 → 添字の対応が何度も要るなら index(v, L) か、配列に書いて使い終わったら戻す。");
}

}  // namespace vt_demo

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
    vt_test::random_test();
    vt_test::benchmark();
    vt_demo::demo();
    return 0;
}
