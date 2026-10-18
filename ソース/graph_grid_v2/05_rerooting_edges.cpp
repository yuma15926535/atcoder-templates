// =============================================================================
//  全方位木DP（辺を切った両側の値つき）  rerooting_edges                    --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 木
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed。G は木（森でもよい。成分ごとに計算する）。辺の番号は G に足した順（read_tree なら入力の順）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  05_rerooting と同じ書き方（merge / put_edge / put_vertex / e）で、頂点ごとの答えに加えて
//  「辺 i を切ったときの両側の木それぞれの dp（切った辺の端を根にする）」も O(N) で返す。
//    例: 辺を 1 本切ったときの両側の直径・大きさ・値の和、辺ごとに「両側の組み合わせ」を数える問題
//    dp[v] = put_vertex( merge( put_edge(dp[c1], 辺 v→c1), put_edge(dp[c2], 辺 v→c2), ... ), v )
//
//  機能一覧                                                                       計算量
//    auto R = rerooting_edges(G, e, merge, put_edge, put_vertex)                               O(N)
//      R.ans[v]        v を根にしたときの dp[v]（05_rerooting の返り値と同じ）
//      R.from_side[i]  辺 i（G.es[i] = {from, to}）を切ったときの from の側の dp（from が根）
//      R.to_side[i]    同じく to の側の dp（to が根）
//    木でないグラフを渡すと BFS で作った全域森で計算する（使わなかった辺の from_side / to_side は e）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    rre_test::random_test() / rre_test::benchmark() / rre_demo::demo()
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
#ifndef GG_REROOTING_EDGES
#define GG_REROOTING_EDGES
// ---- 全方位木DP（辺ごとの値つき）（グラフ・グリッド第二世代）  05_rerooting に「辺を切った両側の値」を足したもの（非再帰・O(N)） ----
// merge / put_edge / put_vertex / e の意味は 05_rerooting と同じ:
//   dp[v] = put_vertex( merge( put_edge(dp[c1], 辺 v→c1), put_edge(dp[c2], 辺 v→c2), ... ), v )
// 返り値:
//   ans[v]       = v を根にしたときの dp[v]（rerooting と同じ）
//   from_side[i] = 辺 i（G.es[i] = {from, to}）を切ったとき、from の側の木を from を根にして計算した dp
//   to_side[i]   = 同じく to の側の木を to を根にして計算した dp
//   （「辺を通して持ち上げた値」が要るときは put_edge(to_side[i], 辺) のように自分でかける）
// G は木（森でもよい）。木でないグラフなら BFS で作った全域森で計算し、使わなかった辺の from_side / to_side は e のまま
template <class T>
struct RerootingEdges {
    std::vector<T> ans, from_side, to_side;
};
template <class T, class Merge, class PutEdge, class PutVertex>
RerootingEdges<T> rerooting_edges(const Graph& G, const T& e, Merge merge, PutEdge put_edge, PutVertex put_vertex) {
    const int n = G.size(), m = G.edge_count();
    std::vector<int> order, par(n, -1), pidx(n, -1), off(n + 1, 0);  // pidx[v] = G[v] での親への辺の位置
    for (int v = 0; v < n; ++v) off[v + 1] = off[v] + (int)G[v].size();
    std::vector<char> child(off[n], 0);  // child[off[v] + i] = 1: G[v][i] は子への木の辺
    order.reserve(n);
    std::vector<char> seen(n, 0);
    for (int r = 0; r < n; ++r) {
        if (seen[r]) continue;
        seen[r] = 1;
        size_t h = order.size();
        order.push_back(r);
        for (; h < order.size(); ++h) {
            const int v = order[h];
            for (int i = 0; i < (int)G[v].size(); ++i) {
                const int u = G[v][i].to;
                if (seen[u]) continue;
                seen[u] = 1, par[u] = v, pidx[u] = -2 - G[v][i].id, child[off[v] + i] = 1, order.push_back(u);  // 親への辺の番号を仮に入れておく
            }
        }
    }
    for (int v = 0; v < n; ++v) {  // 親への辺の位置（G[v] の中で、親から来た辺と同じ番号の辺）
        if (par[v] < 0) continue;
        const int id = -2 - pidx[v];
        for (int i = 0; i < (int)G[v].size(); ++i)
            if (G[v][i].id == id) {
                pidx[v] = i;
                break;
            }
    }
    RerootingEdges<T> res;
    res.ans.assign(n, e), res.from_side.assign(m, e), res.to_side.assign(m, e);
    // 1 回目: 葉から根へ（部分木の値）
    std::vector<T> down(n, e);
    for (int k = n - 1; k >= 0; --k) {
        const int v = order[k];
        T acc = e;
        for (int i = 0; i < (int)G[v].size(); ++i)
            if (child[off[v] + i]) acc = merge(acc, put_edge(down[G[v][i].to], G[v][i]));
        down[v] = put_vertex(acc, v);
    }
    // 2 回目: 根から葉へ。左右からの累積で「自分以外」をまとめる
    std::vector<T> from_par(n, e), suf, c;
    std::vector<char> use;
    for (int k = 0; k < n; ++k) {
        const int v = order[k];
        const int d = (int)G[v].size();
        suf.assign(d + 1, e), c.assign(d, e), use.assign(d, 0);
        for (int i = 0; i < d; ++i) {
            const int u = G[v][i].to;
            if (i == pidx[v]) {
                c[i] = from_par[v], use[i] = 1;
            } else if (child[off[v] + i]) {  // 木の辺で u が子
                c[i] = put_edge(down[u], G[v][i]), use[i] = 1;
            }
        }
        for (int i = d - 1; i >= 0; --i) suf[i] = use[i] ? merge(c[i], suf[i + 1]) : suf[i + 1];
        res.ans[v] = put_vertex(suf[0], v);
        T pre = e;
        for (int i = 0; i < d; ++i) {
            if (!use[i]) continue;
            if (i != pidx[v]) {
                const int u = G[v][i].to;
                const T side = put_vertex(merge(pre, suf[i + 1]), v);  // 辺 v - u を切ったときの v の側（v が根）
                from_par[u] = put_edge(side, G[u][pidx[u]]);
                const int id = G[v][i].id;
                if (G.es[id].from == v) res.from_side[id] = side, res.to_side[id] = down[u];
                else res.to_side[id] = side, res.from_side[id] = down[u];
            }
            pre = merge(pre, c[i]);
        }
    }
    return res;
}
#endif  // GG_REROOTING_EDGES
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "rerooting_edges";
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

namespace rre_test {

// 愚直: 辺 skip を使わずに、v を根にして再帰で木DP（小さい木だけ）
template <class T, class M, class PE, class PV>
static T brute_dp(const Graph& G, int v, int p, int skip, const T& e, M merge, PE put_edge, PV put_vertex) {
    T acc = e;
    for (const Edge& ed : G[v])
        if (ed.to != p && ed.id != skip) acc = merge(acc, put_edge(brute_dp(G, ed.to, v, skip, e, merge, put_edge, put_vertex), ed));
    return put_vertex(acc, v);
}

template <class T, class M, class PE, class PV>
static void check_all(const Graph& G, const T& e, M merge, PE pe, PV pv) {
    const RerootingEdges<T> R = rerooting_edges(G, e, merge, pe, pv);
    const int n = G.size();
    CHECK((int)R.ans.size() == n && (int)R.from_side.size() == G.edge_count() && (int)R.to_side.size() == G.edge_count());
    for (int v = 0; v < n; ++v) CHECK(R.ans[v] == brute_dp(G, v, -1, -1, e, merge, pe, pv));
    for (int id = 0; id < G.edge_count(); ++id) {
        CHECK(R.from_side[id] == brute_dp(G, G.es[id].from, -1, id, e, merge, pe, pv));
        CHECK(R.to_side[id] == brute_dp(G, G.es[id].to, -1, id, e, merge, pe, pv));
    }
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    const long long MOD = 998244353;
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 12);
        Graph G(n);
        std::vector<int> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        for (int i = 1; i < n; ++i) {
            if (it % 5 == 0 && rng() % 4 == 0) continue;  // 森
            int a = perm[i], b = perm[rnd(rng, 0, i - 1)];
            if (rng() & 1) std::swap(a, b);
            G.add_undirected(a, b, rnd(rng, 0, 9));
        }
        std::vector<long long> val(n);
        for (auto& x : val) x = rnd(rng, 0, 9);
        // (1) いちばん遠い頂点までの距離
        check_all(G, 0LL, [](long long a, long long b) { return std::max(a, b); }, [](long long x, const Edge& e) { return x + e.w; }, [](long long x, int) { return x; });
        // (2) {距離の和, 頂点数}
        {
            using P = std::pair<long long, long long>;
            check_all(G, P{0, 0}, [](P a, P b) { return P{a.first + b.first, a.second + b.second}; }, [](P x, const Edge& e) { return P{x.first + x.second * e.w, x.second}; },
                      [](P x, int) { return P{x.first, x.second + 1}; });
        }
        // (3) 根を含む連結な部分集合の数（mod）と頂点の値
        check_all(G, 1LL, [&](long long a, long long b) { return a * b % MOD; }, [&](long long x, const Edge&) { return (x + 1) % MOD; },
                  [&](long long x, int v) { return x * (val[v] + 1) % MOD; });
        // (4) 直径（{根からの深さの上位 2 つ, 部分木の直径}）
        {
            struct D {
                long long a1, a2, dia;
                bool operator==(const D& o) const { return a1 == o.a1 && a2 == o.a2 && dia == o.dia; }
            };
            auto merge = [](D x, D y) {
                long long t[4] = {x.a1, x.a2, y.a1, y.a2};
                std::sort(t, t + 4, std::greater<long long>());
                return D{t[0], t[1], std::max(x.dia, y.dia)};
            };
            auto pe = [](D x, const Edge& e) { return D{x.a1 + e.w, 0, std::max(x.dia, x.a1 + x.a2)}; };
            auto pv = [](D x, int) { return x; };
            check_all(G, D{0, 0, 0}, merge, pe, pv);
        }
    }
    // 木でないグラフ（閉路・多重辺・自己ループ）でも落ちない。全域森の辺以外は e のまま
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 8), m = rnd(rng, 0, 12);
        Graph G(n);
        for (int i = 0; i < m; ++i) G.add_undirected(rnd(rng, 0, n - 1), rnd(rng, 0, n - 1));
        using P = std::pair<long long, long long>;
        const auto R = rerooting_edges(G, P{0, 0}, [](P a, P b) { return P{a.first + b.first, a.second + b.second}; },
                                       [](P x, const Edge& e) { return P{x.first + x.second * e.w, x.second}; }, [](P x, int) { return P{x.first, x.second + 1}; });
        long long tot = 0;
        for (int id = 0; id < m; ++id)
            if (R.from_side[id] != P{0, 0}) tot += R.from_side[id].second + R.to_side[id].second;  // 木の辺なら両側の頂点数の和 = 成分の頂点数
        CHECK(tot >= 0);
    }
    std::cout << "[rerooting_edges] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[rerooting_edges] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> us(n - 1), vs(n - 1);
        std::vector<long long> ws(n - 1);
        for (int i = 1; i < n; ++i) us[i - 1] = i, vs[i - 1] = (int)(rng() % i), ws[i - 1] = (long long)(rng() % 1000);
        const Graph G = graph_from_edges(n, us, vs, ws, UNDIRECTED, 0);
        using P = std::pair<long long, long long>;
        auto t0 = std::chrono::steady_clock::now();
        const auto R = rerooting_edges(G, P{0, 0}, [](P a, P b) { return P{a.first + b.first, a.second + b.second}; },
                                       [](P x, const Edge& e) { return P{x.first + x.second * e.w, x.second}; }, [](P x, int) { return P{x.first, x.second + 1}; });
        const long long t1 = ms_since(t0);
        Graph S(n);  // 星（中心 0）
        for (int v = 1; v < n; ++v) S.add_undirected(v, 0);
        t0 = std::chrono::steady_clock::now();
        const auto R2 = rerooting_edges(S, 0LL, [](long long a, long long b) { return std::max(a, b); }, [](long long x, const Edge& e) { return x + e.w; },
                                        [](long long x, int) { return x; });
        std::cout << "[rerooting_edges] benchmark: N = " << n << "  距離の和（辺ごとの両側つき） " << t1 << " ms, 星で最遠距離 " << ms_since(t0) << " ms  (check "
                  << (R.ans[0].first + R.to_side[0].first + R2.ans[1]) % 1000 << ")\n";
    }
}

}  // namespace rre_test


namespace rre_demo {
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
    banner("全方位木DP（辺ごとの値つき）使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto R = rerooting_edges(G, e, merge, put_edge, put_vertex);", "RerootingEdges<T>", ""},
        {"R.ans[v]", "T", "v を根にしたときの dp（05_rerooting と同じ）"},
        {"R.from_side[i] / R.to_side[i]", "T", "辺 i を切ったときの from 側 / to 側の dp（端が根）"},
    };
    cheat_sheet("API 一覧", rows, 62, 20);

    head("デモで使う木: 辺 0: 0-1, 辺 1: 1-2, 辺 2: 1-3, 辺 3: 3-4（重みなし）");
    RUN(Graph G(5));
    RUN(G.add_undirected(0, 1); G.add_undirected(1, 2); G.add_undirected(1, 3); G.add_undirected(3, 4));

    api("【1】 辺を切ったときの両側の頂点数", "RerootingEdges<long long>", "dp = 部分木の頂点数（merge = 和、put_vertex で + 1）", "O(N)");
    code("auto R = rerooting_edges(G, 0LL, [](long long a, long long b) { return a + b; },");
    code("    [](long long x, const Edge&) { return x; }, [](long long x, int) { return x + 1; });");
    auto R = rerooting_edges(G, 0LL, [](long long a, long long b) { return a + b; }, [](long long x, const Edge&) { return x; }, [](long long x, int) { return x + 1; });
    SHOW(R.from_side, "辺 i の from 側の頂点数");
    SHOW(R.to_side, "辺 i の to 側の頂点数");
    SHOW(R.ans, "どこを根にしても 5");

    api("【2】 辺を切ったときの両側の木の直径", "RerootingEdges<D>", "D = {根からの深さの上位 2 つ a1 ≥ a2, 直径}。put_edge で直径を確定させて持ち上げる", "O(N)");
    code("struct D { long long a1, a2, dia; };");
    code("auto merge = [](D x, D y) { long long t[4] = {x.a1, x.a2, y.a1, y.a2}; std::sort(t, t + 4, std::greater<long long>());");
    code("                            return D{t[0], t[1], std::max(x.dia, y.dia)}; };");
    code("auto put_edge = [](D x, const Edge& e) { return D{x.a1 + e.w, 0, std::max(x.dia, x.a1 + x.a2)}; };");
    code("auto put_vertex = [](D x, int) { return x; };");
    code("auto R2 = rerooting_edges(G, D{0, 0, 0}, merge, put_edge, put_vertex);");
    code("// 辺 i を切った両側の直径: max(from_side[i].dia, a1 + a2)（put_vertex の後はまだ直径を確定させていないので、ここで足す）");
    struct D { long long a1, a2, dia; };
    auto merge = [](D x, D y) {
        long long t[4] = {x.a1, x.a2, y.a1, y.a2};
        std::sort(t, t + 4, std::greater<long long>());
        return D{t[0], t[1], std::max(x.dia, y.dia)};
    };
    auto put_edge = [](D x, const Edge& e) { return D{x.a1 + e.w, 0, std::max(x.dia, x.a1 + x.a2)}; };
    auto put_vertex = [](D x, int) { return x; };
    auto R2 = rerooting_edges(G, D{0, 0, 0}, merge, put_edge, put_vertex);
    std::vector<std::pair<long long, long long>> dia;
    for (int i = 0; i < G.edge_count(); ++i) {
        const D a = R2.from_side[i], b = R2.to_side[i];
        dia.push_back({std::max(a.dia, a.a1 + a.a2), std::max(b.dia, b.a1 + b.a2)});
    }
    SHOW(dia, "辺 i を切ったときの (from 側の直径, to 側の直径)");

    head("書き方のコツ");
    note("dp の値は「辺 i の端を根にして、切った辺の向こうは無いもの」として計算される（put_edge はかかっていない）。");
    note("切った辺を通して持ち上げた値が要るなら put_edge(R.to_side[i], 辺) のように自分でかける。");
    note("頂点ごとの答えだけなら 05_rerooting（少し速い）。");
}

}  // namespace rre_demo

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
    rre_test::random_test();
    rre_test::benchmark();
    rre_demo::demo();
    return 0;
}
