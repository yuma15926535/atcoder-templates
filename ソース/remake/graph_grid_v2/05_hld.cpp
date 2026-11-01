// =============================================================================
//  HLD（パス・部分木の集計と更新）  TreeQuery / TreeLazyQuery          --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 木
//  依存       : 00_graph_base, 05_rooted_tree, 05_lca（同梱済み。LCA の重軽分解をそのまま使う）
//  頂点番号   : 0-indexed。G は木（森でもよい: LCA L(RootedTree::forest(G)) で作る）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  木の頂点（または辺）に値を載せ、「u から v へのパス」「v の部分木」の集計と更新をする。セグ木も中に入っている。
//    - 演算 op は結合的であればよい（可換でなくてよい。path(u, v) は u から v への順に畳み込む）
//    - edge = true なら辺に値を載せる（辺の値は子の側の頂点に置く。パスの LCA の頂点は含めない）
//    - 部分木は HLD の順に畳み込む（可換な演算で使う）
//
//  機能一覧                                                                       計算量
//    TreeQuery Q(L, init, op, e, edge = false, commutative = false)   1 点の更新 + 集計（init[v] = 頂点 v の値）  O(N)
//      Q.set(v, x) / Q.get(v)                       頂点 v（edge なら v と親の辺）の値            O(log N)
//      Q.set_edge(u, v, x) / Q.get_edge(u, v)       辺 u - v の値（edge のとき）                  O(log N)
//      Q.path(u, v)                                 u → v のパスを順に畳み込む                    O(log^2 N)
//      Q.subtree(v)                                 v の部分木を畳み込む                          O(log N)
//    TreeLazyQuery Q(L, init, op, e, mapping, composition, id, edge = false, commutative = false)   区間への作用もできる版
//      Q.path_apply(u, v, f) / Q.subtree_apply(v, f) / Q.path / Q.subtree / Q.set / Q.get   O(log^2 N) / O(log N)
//      mapping(f, x) = 作用 f を値 x にかけたもの、composition(f, g) = g のあとに f、id = 何もしない作用（型が作用の型になる。0LL など）
//    commutative = true: 演算が可換（和・最大・最小など）なら付ける。逆向きの畳み込みを持たないので 2 〜 3 割速い
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤・根付き木・LCA も入っている）
//    hld_test::random_test() / hld_test::benchmark() / hld_demo::demo()
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
#ifndef GG_HLD
#define GG_HLD
// ---- HLD（重軽分解）でパス・部分木の集計と更新（グラフ・グリッド第二世代）  05_lca の path_segments の上にセグ木を載せる ----
// 頂点に値を載せる（edge = false）か、辺に値を載せる（edge = true。辺の値は子の側の頂点に置く）。
// 演算 op は結合的であればよい（可換でなくてよい。path(u, v) は u から v への順に畳み込む）。部分木は HLD の順に畳み込む
namespace gg_internal {
template <class T>
struct Identity { using type = T; };  // 型推論を止める（e に 0 と書いても S = long long のままにする）
// 前向き・後ろ向きの両方の畳み込みを持つセグ木（非可換の演算でパスを逆向きにたどるため）
template <class S, class Op>
struct BiSegTree {
    int n = 0, size = 1;
    bool both = true;     // 降順の畳み込みも持つか（可換な演算なら要らない）
    std::vector<S> f, b;  // f = 位置の昇順の畳み込み、b = 降順の畳み込み
    Op op;
    S e;
    BiSegTree(const std::vector<S>& v, Op op_, S e_, bool both_) : n((int)v.size()), both(both_), op(op_), e(e_) {
        while (size < n) size <<= 1;
        f.assign(2 * size, e);
        if (both) b.assign(2 * size, e);
        for (int i = 0; i < n; ++i) {
            f[size + i] = v[i];
            if (both) b[size + i] = v[i];
        }
        for (int k = size - 1; k >= 1; --k) pull(k);
    }
    void pull(int k) {
        f[k] = op(f[2 * k], f[2 * k + 1]);
        if (both) b[k] = op(b[2 * k + 1], b[2 * k]);
    }
    void set(int p, const S& x) {
        p += size, f[p] = x;
        if (both) b[p] = x;
        for (p >>= 1; p >= 1; p >>= 1) pull(p);
    }
    S get(int p) const { return f[p + size]; }
    // [l, r) を位置の昇順（rev なら降順）に畳み込む
    S fold(int l, int r, bool rev) const {
        S x = e, y = e;
        if (!rev) {
            for (l += size, r += size; l < r; l >>= 1, r >>= 1) {
                if (l & 1) x = op(x, f[l++]);
                if (r & 1) y = op(f[--r], y);
            }
            return op(x, y);
        }
        for (l += size, r += size; l < r; l >>= 1, r >>= 1) {
            if (l & 1) x = op(b[l++], x);
            if (r & 1) y = op(y, b[--r]);
        }
        return op(y, x);
    }
};
// 遅延伝播つき（区間に作用 f をまとめてかける）。前向き・後ろ向きの両方に同じ作用をかける
template <class S, class F, class Op, class Mapping, class Composition>
struct BiLazySegTree {
    int n = 0, size = 1, lg = 0;
    bool both = true;
    std::vector<S> f, b;
    std::vector<F> lz;
    Op op;
    S e;
    Mapping mapping;
    Composition composition;
    F id;
    BiLazySegTree(const std::vector<S>& v, Op op_, S e_, Mapping m_, Composition c_, F id_, bool both_)
        : n((int)v.size()), both(both_), op(op_), e(e_), mapping(m_), composition(c_), id(id_) {
        while (size < n) size <<= 1, ++lg;
        f.assign(2 * size, e), lz.assign(size, id);
        if (both) b.assign(2 * size, e);
        for (int i = 0; i < n; ++i) {
            f[size + i] = v[i];
            if (both) b[size + i] = v[i];
        }
        for (int k = size - 1; k >= 1; --k) pull(k);
    }
    void pull(int k) {
        f[k] = op(f[2 * k], f[2 * k + 1]);
        if (both) b[k] = op(b[2 * k + 1], b[2 * k]);
    }
    void all_apply(int k, const F& g) {
        f[k] = mapping(g, f[k]);
        if (both) b[k] = mapping(g, b[k]);
        if (k < size) lz[k] = composition(g, lz[k]);
    }
    void push(int k) {
        all_apply(2 * k, lz[k]), all_apply(2 * k + 1, lz[k]);
        lz[k] = id;
    }
    void set(int p, const S& x) {
        p += size;
        for (int i = lg; i >= 1; --i) push(p >> i);
        f[p] = x;
        if (both) b[p] = x;
        for (int i = 1; i <= lg; ++i) pull(p >> i);
    }
    S get(int p) {
        p += size;
        for (int i = lg; i >= 1; --i) push(p >> i);
        return f[p];
    }
    S fold(int l, int r, bool rev) {
        if (l == r) return e;
        l += size, r += size;
        for (int i = lg; i >= 1; --i) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        S x = e, y = e;
        if (!rev) {
            for (; l < r; l >>= 1, r >>= 1) {
                if (l & 1) x = op(x, f[l++]);
                if (r & 1) y = op(f[--r], y);
            }
            return op(x, y);
        }
        for (; l < r; l >>= 1, r >>= 1) {
            if (l & 1) x = op(b[l++], x);
            if (r & 1) y = op(y, b[--r]);
        }
        return op(y, x);
    }
    void apply(int l, int r, const F& g) {
        if (l == r) return;
        l += size, r += size;
        for (int i = lg; i >= 1; --i) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        {
            int l2 = l, r2 = r;
            for (; l < r; l >>= 1, r >>= 1) {
                if (l & 1) all_apply(l++, g);
                if (r & 1) all_apply(--r, g);
            }
            l = l2, r = r2;
        }
        for (int i = 1; i <= lg; ++i) {
            if (((l >> i) << i) != l) pull(l >> i);
            if (((r >> i) << i) != r) pull((r - 1) >> i);
        }
    }
};
}  // namespace gg_internal

// 1 点の更新 + パス・部分木の集計。O(log^2 N) / O(log N)
template <class S, class Op>
struct TreeQuery {
    LCA L;
    bool edge = false, commutative = false;
    gg_internal::BiSegTree<S, Op> seg;
    // init[v] = 頂点 v の値（edge なら「v と親を結ぶ辺の値」。根の分は使わない）。
    // commutative = true（演算が可換: 和・最大など）なら逆向きの畳み込みを持たず、少し速く省メモリになる
    TreeQuery(const LCA& L_, const std::vector<S>& init, Op op, typename gg_internal::Identity<S>::type e, bool edge_ = false, bool commutative_ = false)
        : L(L_), edge(edge_), commutative(commutative_), seg(arrange(L_, init, e), op, e, !commutative_) {}
    static std::vector<S> arrange(const LCA& L, const std::vector<S>& init, const S& e) {
        std::vector<S> v(L.n, e);
        for (int x = 0; x < L.n; ++x)
            if (L.hld_pos[x] >= 0) v[L.hld_pos[x]] = init[x];
        return v;
    }
    void set(int v, const S& x) { seg.set(L.hld_pos[v], x); }  // 頂点 v の値（edge なら v と親の辺の値）を x に
    S get(int v) const { return seg.get(L.hld_pos[v]); }
    // 辺 u - v の値を x に（edge のとき。u と v は隣どうし）
    void set_edge(int u, int v, const S& x) { set(L.par[u] == v ? u : v, x); }
    S get_edge(int u, int v) const { return get(L.par[u] == v ? u : v); }
    // u から v へのパスを、その順に op で畳み込む（edge なら通る辺の値）。別の木なら e
    S path(int u, int v) const {
        S r = seg.e;
        for (auto [a, b] : L.path_segments(u, v, edge)) r = seg.op(r, a <= b ? seg.fold(a, b + 1, false) : seg.fold(b, a + 1, !commutative));
        return r;
    }
    // v の部分木（edge なら部分木の中の辺）を HLD の順に畳み込む
    S subtree(int v) const {
        const auto [l, r] = L.subtree_range(v);
        return seg.fold(l + (edge ? 1 : 0), r, false);
    }
};
// 区間への作用（パス・部分木にまとめて足す・代入する など）+ 集計。O(log^2 N) / O(log N)
// mapping(f, x): 作用 f を値 x にかけた結果、composition(f, g): 「g のあとに f」をまとめた作用、id: 何もしない作用
template <class S, class F, class Op, class Mapping, class Composition>
struct TreeLazyQuery {
    LCA L;
    bool edge = false, commutative = false;
    gg_internal::BiLazySegTree<S, F, Op, Mapping, Composition> seg;
    // id（何もしない作用）の型が作用の型 F になる（long long なら 0LL と書く）
    TreeLazyQuery(const LCA& L_, const std::vector<S>& init, Op op, typename gg_internal::Identity<S>::type e, Mapping mapping, Composition composition,
                  F id, bool edge_ = false, bool commutative_ = false)
        : L(L_), edge(edge_), commutative(commutative_), seg(TreeQuery<S, Op>::arrange(L_, init, e), op, e, mapping, composition, id, !commutative_) {}
    void set(int v, const S& x) { seg.set(L.hld_pos[v], x); }
    S get(int v) { return seg.get(L.hld_pos[v]); }
    void set_edge(int u, int v, const S& x) { set(L.par[u] == v ? u : v, x); }
    S get_edge(int u, int v) { return get(L.par[u] == v ? u : v); }
    S path(int u, int v) {
        S r = seg.e;
        for (auto [a, b] : L.path_segments(u, v, edge)) r = seg.op(r, a <= b ? seg.fold(a, b + 1, false) : seg.fold(b, a + 1, !commutative));
        return r;
    }
    S subtree(int v) {
        const auto [l, r] = L.subtree_range(v);
        return seg.fold(l + (edge ? 1 : 0), r, false);
    }
    // u - v のパス上の頂点（edge なら辺）すべてに作用 f をかける
    void path_apply(int u, int v, const F& f) {
        for (auto [a, b] : L.path_segments(u, v, edge)) seg.apply(std::min(a, b), std::max(a, b) + 1, f);
    }
    // v の部分木（edge なら部分木の中の辺）すべてに作用 f をかける
    void subtree_apply(int v, const F& f) {
        const auto [l, r] = L.subtree_range(v);
        seg.apply(l + (edge ? 1 : 0), r, f);
    }
};
#endif  // GG_HLD
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "hld";
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

namespace hld_test {

const long long P = 1000000007;
using Aff = std::pair<long long, long long>;  // x → a x + b（mod P）
static Aff aff_op(Aff l, Aff r) { return {r.first * l.first % P, (r.first * l.second + r.second) % P}; }  // l のあとに r

struct Seq { long long first, last, sum; int cnt; };  // 列の先頭・末尾・和・個数（非可換）
static Seq seq_op(Seq a, Seq b) {
    if (a.cnt == 0) return b;
    if (b.cnt == 0) return a;
    return {a.first, b.last, a.sum + b.sum, a.cnt + b.cnt};
}
static bool operator==(const Seq& a, const Seq& b) { return a.first == b.first && a.last == b.last && a.sum == b.sum && a.cnt == b.cnt; }

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 25);
        Graph G(n);
        std::vector<int> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        for (int i = 1; i < n; ++i) {
            if (it % 6 == 0 && rng() % 4 == 0) continue;  // 森
            G.add_undirected(perm[i], perm[(it % 3 == 0) ? i - 1 : rnd(rng, 0, i - 1)]);
        }
        const RootedTree T = RootedTree::forest(G);
        const LCA L(T);
        const bool edge = (it % 2 == 1);
        // パス上の「値を置く頂点」の列（edge なら LCA を除く。u → v の順）
        auto path_nodes = [&](int u, int v) {
            std::vector<int> p = L.path(u, v);
            if (edge && !p.empty()) {
                const int a = L.lca(u, v);
                p.erase(std::find(p.begin(), p.end(), a));
            }
            return p;
        };
        auto sub_nodes = [&](int v) {
            std::vector<int> r;
            for (int x = 0; x < n; ++x)
                if (L.is_ancestor(v, x) && !(edge && x == v)) r.push_back(x);
            return r;
        };
        // (1) 非可換（一次関数の合成）、1 点更新
        {
            std::vector<Aff> val(n);
            for (auto& x : val) x = {rnd(rng, 0, 9), rnd(rng, 0, 9)};
            TreeQuery Q(L, val, aff_op, Aff{1, 0}, edge);
            for (int q = 0; q < 30; ++q) {
                const int t = rnd(rng, 0, 3), u = rnd(rng, 0, n - 1), v = rnd(rng, 0, n - 1);
                if (t == 0) {
                    const Aff x{rnd(rng, 0, 9), rnd(rng, 0, 9)};
                    if (edge && L.par[u] >= 0 && rng() % 2) Q.set_edge(L.par[u], u, x), val[u] = x;
                    else Q.set(u, x), val[u] = x;
                } else if (t == 1) {
                    Aff ex{1, 0};
                    if (L.same_tree(u, v))
                        for (int x : path_nodes(u, v)) ex = aff_op(ex, val[x]);
                    CHECK(Q.path(u, v) == ex);
                } else if (t == 2) {
                    CHECK(Q.get(u) == val[u]);
                    if (edge && L.par[u] >= 0) CHECK(Q.get_edge(u, L.par[u]) == val[u]);
                } else {
                    // 部分木は可換な演算で使うもの。ここでは HLD の順（hld_pos の昇順）に並べた愚直と比べる
                    std::vector<int> s = sub_nodes(u);
                    std::sort(s.begin(), s.end(), [&](int a, int b) { return L.hld_pos[a] < L.hld_pos[b]; });
                    Aff ex{1, 0};
                    for (int x : s) ex = aff_op(ex, val[x]);
                    CHECK(Q.subtree(u) == ex);
                }
            }
        }
        // (2) 遅延: 区間に足す。値は列の要約（先頭・末尾・和・個数）なので向きも確かめられる
        {
            std::vector<long long> val(n);
            for (auto& x : val) x = rnd(rng, -9, 9);
            std::vector<Seq> init(n);
            for (int x = 0; x < n; ++x) init[x] = {val[x], val[x], val[x], 1};
            auto mapping = [](long long f, Seq s) { return s.cnt == 0 ? s : Seq{s.first + f, s.last + f, s.sum + f * s.cnt, s.cnt}; };
            auto composition = [](long long f, long long g) { return f + g; };
            TreeLazyQuery Q(L, init, seq_op, Seq{0, 0, 0, 0}, mapping, composition, 0LL, edge);
            for (int q = 0; q < 30; ++q) {
                const int t = rnd(rng, 0, 5), u = rnd(rng, 0, n - 1), v = rnd(rng, 0, n - 1);
                const long long f = rnd(rng, -5, 5);
                if (t == 0) {
                    if (L.same_tree(u, v)) {
                        Q.path_apply(u, v, f);
                        for (int x : path_nodes(u, v)) val[x] += f;
                    }
                } else if (t == 1) {
                    Q.subtree_apply(u, f);
                    for (int x : sub_nodes(u)) val[x] += f;
                } else if (t == 2) {
                    Seq ex{0, 0, 0, 0};
                    if (L.same_tree(u, v))
                        for (int x : path_nodes(u, v)) ex = seq_op(ex, Seq{val[x], val[x], val[x], 1});
                    CHECK(Q.path(u, v) == ex);
                } else if (t == 3) {
                    long long sum = 0;
                    int cnt = 0;
                    for (int x : sub_nodes(u)) sum += val[x], ++cnt;
                    const Seq got = Q.subtree(u);
                    CHECK(got.sum == sum && got.cnt == cnt);
                } else if (t == 4) {
                    const long long x = rnd(rng, -9, 9);
                    Q.set(u, Seq{x, x, x, 1}), val[u] = x;
                } else {
                    const Seq g = Q.get(u);
                    CHECK(g.first == val[u] && g.cnt == 1);
                }
            }
        }
    }
    // (3) 可換な演算で commutative = true（和と最大）。1 点更新・区間に足す
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 25);
        Graph G(n);
        for (int i = 1; i < n; ++i) G.add_undirected(i, rnd(rng, 0, i - 1));
        const LCA L(G, rnd(rng, 0, n - 1));
        const bool edge = (it % 2 == 1);
        std::vector<long long> val(n);
        for (auto& x : val) x = rnd(rng, -9, 9);
        TreeQuery Q(L, val, [](long long a, long long b) { return std::max(a, b); }, -DIST_INF, edge, true);
        using S2 = std::pair<long long, long long>;
        std::vector<S2> init(n);
        for (int x = 0; x < n; ++x) init[x] = {val[x], 1};
        TreeLazyQuery Z(L, init, [](S2 a, S2 b) { return S2{a.first + b.first, a.second + b.second}; }, S2{0, 0},
                        [](long long f, S2 x) { return S2{x.first + f * x.second, x.second}; }, [](long long f, long long g) { return f + g; }, 0LL, edge, true);
        std::vector<long long> mx = val, sm = val;
        for (int q = 0; q < 30; ++q) {
            const int t = rnd(rng, 0, 3), u = rnd(rng, 0, n - 1), v = rnd(rng, 0, n - 1);
            std::vector<int> p = L.path(u, v);
            if (edge) p.erase(std::find(p.begin(), p.end(), L.lca(u, v)));
            if (t == 0) {
                const long long x = rnd(rng, -9, 9);
                Q.set(u, x), mx[u] = x;
            } else if (t == 1) {
                const long long f = rnd(rng, -3, 3);
                Z.path_apply(u, v, f);
                for (int x : p) sm[x] += f;
            } else {
                long long a = -DIST_INF, b = 0;
                for (int x : p) a = std::max(a, mx[x]), b += sm[x];
                CHECK(Q.path(u, v) == a && Z.path(u, v).first == b && Z.path(u, v).second == (long long)p.size());
            }
        }
    }
    std::cout << "[hld] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[hld] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    const int n = 200000;
    for (int line = 0; line < 2; ++line) {
        Graph G(n);
        for (int i = 1; i < n; ++i) G.add_undirected(i, line ? i - 1 : (int)(rng() % i));
        const LCA L(G, 0);
        std::vector<long long> init(n);
        for (auto& x : init) x = (long long)(rng() % 1000);
        auto t0 = std::chrono::steady_clock::now();
        TreeQuery Q(L, init, [](long long a, long long b) { return a + b; }, 0LL, false, true);
        long long chk = 0;
        for (int q = 0; q < n; ++q) {
            const int u = (int)(rng() % n), v = (int)(rng() % n);
            if (q & 1) Q.set(u, (long long)(rng() % 1000));
            else chk += Q.path(u, v);
        }
        const long long t1 = ms_since(t0);
        using S = std::pair<long long, long long>;  // {和, 個数}
        std::vector<S> init2(n);
        for (int x = 0; x < n; ++x) init2[x] = {init[x], 1};
        t0 = std::chrono::steady_clock::now();
        TreeLazyQuery Z(L, init2, [](S a, S b) { return S{a.first + b.first, a.second + b.second}; }, S{0, 0},
                        [](long long f, S s) { return S{s.first + f * s.second, s.second}; }, [](long long f, long long g) { return f + g; }, 0LL, false, true);
        for (int q = 0; q < n; ++q) {
            const int u = (int)(rng() % n), v = (int)(rng() % n);
            if (q & 1) Z.path_apply(u, v, (long long)(rng() % 10));
            else chk += Z.path(u, v).first;
        }
        std::cout << "[hld] benchmark: N = Q = " << n << (line ? "（一直線）" : "（ランダム）") << "  TreeQuery（1 点更新 + パスの和）" << t1 << " ms, TreeLazyQuery（パスに足す + パスの和）"
                  << ms_since(t0) << " ms  (check " << chk % 1000 << ")\n";
    }
}

}  // namespace hld_test


namespace hld_demo {
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
    banner("HLD（パス・部分木の集計と更新） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"LCA L(G, root); / LCA L(RootedTree::forest(G));", "", "まず LCA を作る（森なら forest）"},
        {"TreeQuery Q(L, init, op, e, edge = false, comm = false);", "", "1 点更新 + パス・部分木の集計"},
        {"Q.set(v, x) / Q.get(v) / Q.set_edge(u, v, x)", "", "頂点（辺）の値"},
        {"Q.path(u, v) / Q.subtree(v)", "S", "u → v の順に畳み込む / 部分木"},
        {"TreeLazyQuery Q(L, init, op, e, mapping, composition, id, edge, comm);", "", "パス・部分木への作用もできる版"},
        {"comm = true", "", "演算が可換（和・最大など）なら付けると速い"},
        {"Q.path_apply(u, v, f) / Q.subtree_apply(v, f)", "", "パス / 部分木のすべてに作用 f"},
    };
    cheat_sheet("API 一覧", rows, 64, 6);

    head("デモで使う木: 0-1, 1-2, 1-3, 3-4（根 0）、頂点の値 {5, 3, 7, 1, 4}");
    RUN(Graph G(5));
    RUN(G.add_undirected(0, 1); G.add_undirected(1, 2); G.add_undirected(1, 3); G.add_undirected(3, 4));
    RUN(LCA L(G, 0));

    api("【1】 TreeQuery（頂点の値の和・最大）", "long long", "パス上の頂点の値を畳み込む。1 点ずつ更新できる", "path O(log^2 N) / set O(log N)");
    RUN(TreeQuery Q(L, std::vector<long long>{5, 3, 7, 1, 4}, [](long long a, long long b) { return a + b; }, 0));
    SHOW(Q.path(2, 4), "2 → 1 → 3 → 4 で 7 + 3 + 1 + 4");
    SHOW(Q.subtree(1), "部分木 {1, 2, 3, 4}");
    RUN(Q.set(3, 10));
    SHOW(Q.path(2, 4), "頂点 3 を 10 にした");
    RUN(TreeQuery M(L, std::vector<long long>{5, 3, 7, 1, 4}, [](long long a, long long b) { return std::max(a, b); }, -DIST_INF));
    SHOW(M.path(4, 0), "パス上の最大");
    note("和・最大のように可換な演算なら、最後の引数 commutative を true にすると 2 〜 3 割速くなる: TreeQuery Q(L, init, op, 0, false, true);");

    api("【2】 辺に値を載せる（edge = true）", "long long", "辺の値は子の側の頂点に置く。重みつきの木の距離・パス上の辺の最大など（ABC294 G）", "同上");
    code("// 辺の重み: 0-1 (2), 1-2 (3), 1-3 (1), 3-4 (4) → 子の側の頂点に置く: init[1] = 2, init[2] = 3, init[3] = 1, init[4] = 4");
    RUN(TreeQuery E(L, std::vector<long long>{0, 2, 3, 1, 4}, [](long long a, long long b) { return a + b; }, 0, true));
    SHOW(E.path(2, 4), "辺 2-1, 1-3, 3-4 の和 = 3 + 1 + 4（LCA の頂点 1 は数えない）");
    RUN(E.set_edge(1, 3, 100));
    SHOW(E.path(2, 4), "辺 1-3 を 100 にした");

    api("【3】 非可換な演算（向きが大事なもの）", "S", "path(u, v) は u から v への順に畳み込むので、文字列の連結・一次関数の合成などにも使える", "同上");
    code("// 頂点に文字を置き、u → v のパスの文字列を作る");
    RUN(TreeQuery W(L, std::vector<std::string>{"a", "b", "c", "d", "e"}, [](std::string x, std::string y) { return x + y; }, std::string()));
    SHOW(W.path(2, 4), "c → b → d → e");
    SHOW(W.path(4, 2), "逆向き");

    api("【4】 TreeLazyQuery（パスにまとめて足す）", "S", "値を {和, 個数} にし、作用「全部に f を足す」を mapping で書く", "O(log^2 N)");
    code("using S = std::pair<long long, long long>;   // {和, 個数}");
    code("TreeLazyQuery Z(L, init, [](S a, S b) { return S{a.first + b.first, a.second + b.second}; }, S{0, 0},");
    code("    [](long long f, S s) { return S{s.first + f * s.second, s.second}; },   // mapping: 全部に f を足す");
    code("    [](long long f, long long g) { return f + g; }, 0LL, false, true);       // composition と id、辺ではない、可換");
    using S = std::pair<long long, long long>;
    std::vector<S> init = {{5, 1}, {3, 1}, {7, 1}, {1, 1}, {4, 1}};
    TreeLazyQuery Z(L, init, [](S a, S b) { return S{a.first + b.first, a.second + b.second}; }, S{0, 0},
                    [](long long f, S s) { return S{s.first + f * s.second, s.second}; }, [](long long f, long long g) { return f + g; }, 0LL, false, true);
    RUN(Z.path_apply(0, 4, 10));
    SHOW(Z.path(2, 4).first, "0 → 1 → 3 → 4 に 10 ずつ足した後の、2 → 4 の和（7 + 13 + 11 + 14）");
    RUN(Z.subtree_apply(3, -1));
    SHOW(Z.subtree(1).first, "部分木 {1, 2, 3, 4} の和");

    head("典型的な使い方");
    note("パス上の頂点の和・最大（1 点更新）: TreeQuery");
    note("辺の重みを変えながら 2 点の距離 : TreeQuery（edge = true）");
    note("パスにまとめて足す・代入する     : TreeLazyQuery");
    note("更新がなく、距離・LCA だけなら   : 05_lca の wdist / lca で十分");
}

}  // namespace hld_demo

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
    hld_test::random_test();
    hld_test::benchmark();
    hld_demo::demo();
    return 0;
}
