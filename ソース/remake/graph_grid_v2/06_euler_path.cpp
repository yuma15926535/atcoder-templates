// =============================================================================
//  オイラー路（一筆書き）  euler_path_directed / euler_path_undirected / euler_circuit_*   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 06 特殊グラフ
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed。多重辺・自己ループ可。辺のない頂点（孤立点）は無視する
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  すべての辺をちょうど 1 回ずつ通る道（オイラー路）を求める。できる条件:
//    有向: 辺のある頂点がつながっていて、(出次数 - 入次数) が全部 0（閉路）か、+1 と -1 が 1 つずつで残りが 0
//    無向: 辺のある頂点がつながっていて、次数が奇数の頂点が 0 個（閉路）か 2 個（その 2 つが端）
//  有向版は G.es の向き（from → to）を、無向版は G.es の各辺を向きなしとして使う（add_edge で作っても無向として扱える）
//
//  機能一覧                                                                       計算量
//    euler_path_directed(G, s = -1, lex_min = false)     有向グラフの一筆書き           O(N + M)
//    euler_path_undirected(G, s = -1, lex_min = false)   無向グラフの一筆書き           O(N + M)
//    euler_circuit_directed / euler_circuit_undirected   始点に戻ってくるもの（オイラー閉路）
//      s を指定するとその頂点から始める。lex_min = true なら頂点の列が辞書順最小（O(M log M)）
//      結果 EulerPath r: r.ok（できたか）/ r.vs（頂点の列、辺の数 + 1 個）/ r.es（辺番号の列）/ r.closed()
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    ep_test::random_test() / ep_test::benchmark() / ep_demo::demo()
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
#ifndef GG_EULER_PATH
#define GG_EULER_PATH
// ---- オイラー路（一筆書き）（グラフ・グリッド第二世代）  Hierholzer 法（非再帰）。すべての辺をちょうど 1 回ずつ通る ----
// 有向版は G.es の向き（from → to）に、無向版は G.es の各辺を向きなしとして使う。孤立点は無視する
struct EulerPath {
    bool ok = false;        // 一筆書きできたか
    std::vector<int> vs;    // 通る頂点の列（辺の本数 + 1 個。辺がなければ始点 1 つ）
    std::vector<int> es;    // 通る辺の番号の列（es[i] は vs[i] → vs[i+1]）
    bool closed() const { return ok && !vs.empty() && vs.front() == vs.back(); }  // 始点に戻ってくるか（オイラー閉路）
};
namespace gg_internal {
inline EulerPath euler_impl(const Graph& G, int s, bool directed, bool lex_min, bool need_closed) {
    const int n = G.size(), m = G.edge_count();
    EulerPath r;
    if (n == 0) return r;
    // 隣接（CSR）: 有向なら from に、無向なら両端に辺番号を入れる
    std::vector<int> start(n + 1, 0), adj(directed ? m : 2 * (size_t)m);
    std::vector<int> out(n, 0), in(n, 0), deg(n, 0);
    for (const RawEdge& e : G.es) {
        ++out[e.from], ++in[e.to], deg[e.from] += 1, deg[e.to] += 1;
        ++start[e.from + 1];
        if (!directed) ++start[e.to + 1];
    }
    for (int v = 0; v < n; ++v) start[v + 1] += start[v];
    {
        std::vector<int> pos(start.begin(), start.end() - 1);
        for (int id = 0; id < m; ++id) {
            adj[pos[G.es[id].from]++] = id;
            if (!directed) adj[pos[G.es[id].to]++] = id;
        }
    }
    auto other = [&](int id, int v) { return G.es[id].from == v ? G.es[id].to : G.es[id].from; };
    if (lex_min)
        for (int v = 0; v < n; ++v)
            std::sort(adj.begin() + start[v], adj.begin() + start[v + 1], [&](int a, int b) {
                const int x = other(a, v), y = other(b, v);
                return x != y ? x < y : a < b;
            });
    // 始点を決める（次数の条件を確かめる）
    int first = -1, plus = -1, odd = 0;
    for (int v = 0; v < n; ++v) {
        if (deg[v] > 0 && first < 0) first = v;
        if (directed) {
            const int d = out[v] - in[v];
            if (d == 1) {
                if (plus >= 0) return r;
                plus = v;
            } else if (d == -1) {
                ++odd;
                if (odd > 1) return r;
            } else if (d != 0) {
                return r;
            }
        } else if (deg[v] % 2 == 1) {
            ++odd;
            if (plus < 0) plus = v;
        }
    }
    if (!directed && odd != 0 && odd != 2) return r;
    if (directed && (plus >= 0) != (odd == 1)) return r;
    const bool balanced = directed ? (plus < 0) : (odd == 0);
    if (need_closed && !balanced) return r;
    int st0;
    if (s >= 0) {
        if (!balanced) {
            if (directed && s != plus) return r;
            if (!directed && deg[s] % 2 == 0) return r;
        } else if (m > 0 && deg[s] == 0) {
            return r;
        }
        st0 = s;
    } else {
        st0 = !balanced ? plus : (first >= 0 ? first : 0);
    }
    // Hierholzer 法（非再帰）
    std::vector<int> it(start.begin(), start.end() - 1);
    std::vector<char> used(m, 0);
    std::vector<std::pair<int, int>> stk{{st0, -1}};  // (頂点, そこへ来た辺)
    while (!stk.empty()) {
        const int v = stk.back().first;
        while (it[v] < start[v + 1] && used[adj[it[v]]]) ++it[v];
        if (it[v] == start[v + 1]) {
            r.vs.push_back(v), r.es.push_back(stk.back().second);
            stk.pop_back();
            continue;
        }
        const int id = adj[it[v]++];
        used[id] = 1;
        stk.push_back({directed ? G.es[id].to : other(id, v), id});
    }
    if ((int)r.vs.size() != m + 1) {  // つながっていない辺がある
        r.vs.clear(), r.es.clear();
        return r;
    }
    std::reverse(r.vs.begin(), r.vs.end());
    std::reverse(r.es.begin(), r.es.end());
    r.es.erase(r.es.begin());  // 始点の「来た辺」(-1) を除く
    r.ok = true;
    return r;
}
}  // namespace gg_internal
// 有向グラフの一筆書き（s = -1 なら始点を自動で選ぶ）。lex_min = true なら頂点の列が辞書順最小のもの。O(N + M)（lex_min は O(M log M)）
inline EulerPath euler_path_directed(const Graph& G, int s = -1, bool lex_min = false) { return gg_internal::euler_impl(G, s, true, lex_min, false); }
// 無向グラフの一筆書き
inline EulerPath euler_path_undirected(const Graph& G, int s = -1, bool lex_min = false) { return gg_internal::euler_impl(G, s, false, lex_min, false); }
// 始点に戻ってくる一筆書き（オイラー閉路）
inline EulerPath euler_circuit_directed(const Graph& G, int s = -1, bool lex_min = false) { return gg_internal::euler_impl(G, s, true, lex_min, true); }
inline EulerPath euler_circuit_undirected(const Graph& G, int s = -1, bool lex_min = false) { return gg_internal::euler_impl(G, s, false, lex_min, true); }
#endif  // GG_EULER_PATH
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "euler_path";
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

namespace ep_test {

// 愚直: 全部の一筆書き（頂点の列）を列挙する
static std::vector<std::vector<int>> all_paths(const Graph& G, bool directed) {
    const int n = G.size(), m = G.edge_count();
    std::vector<std::vector<int>> res;
    std::vector<int> used(m, 0), path;
    std::function<void(int, int)> rec = [&](int v, int k) {
        if (k == m) {
            res.push_back(path);
            return;
        }
        for (int id = 0; id < m; ++id) {
            if (used[id]) continue;
            const RawEdge& e = G.es[id];
            int to = -1;
            if (e.from == v) to = e.to;
            else if (!directed && e.to == v) to = e.from;
            if (to < 0) continue;
            used[id] = 1, path.push_back(to);
            rec(to, k + 1);
            used[id] = 0, path.pop_back();
        }
    };
    for (int s = 0; s < n; ++s) {
        path.assign(1, s);
        rec(s, 0);
    }
    return res;
}

static void check_path(const Graph& G, const EulerPath& r, bool directed) {
    const int m = G.edge_count();
    CHECK((int)r.es.size() == m && (int)r.vs.size() == m + 1);
    std::vector<int> used(m, 0);
    for (int i = 0; i < m; ++i) {
        const RawEdge& e = G.es[r.es[i]];
        CHECK(!used[r.es[i]]);
        used[r.es[i]] = 1;
        const bool fw = (e.from == r.vs[i] && e.to == r.vs[i + 1]);
        const bool bw = (!directed && e.to == r.vs[i] && e.from == r.vs[i + 1]);
        CHECK(fw || bw);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 5), m = rnd(rng, 0, 7);
        Graph G(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            (rng() & 1) ? G.add_edge(a, b) : G.add_undirected(a, b);
        }
        for (bool directed : {true, false}) {
            const std::vector<std::vector<int>> all = all_paths(G, directed);
            for (int s = -1; s < n; ++s)
                for (bool closed : {false, true})
                    for (bool lex : {false, true}) {
                        std::vector<std::vector<int>> cand;
                        for (const auto& p : all)
                            if ((s < 0 || p[0] == s) && (!closed || p.front() == p.back())) cand.push_back(p);
                        if (m == 0 && (s >= 0 || n > 0)) cand.assign(1, std::vector<int>{s < 0 ? 0 : s});
                        const EulerPath r = directed ? (closed ? euler_circuit_directed(G, s, lex) : euler_path_directed(G, s, lex))
                                                     : (closed ? euler_circuit_undirected(G, s, lex) : euler_path_undirected(G, s, lex));
                        CHECK(r.ok == !cand.empty());
                        if (!r.ok) {
                            CHECK(r.vs.empty() && r.es.empty());
                            continue;
                        }
                        check_path(G, r, directed);
                        if (s >= 0) CHECK(r.vs[0] == s);
                        if (closed) CHECK(r.closed());
                        CHECK(std::find(cand.begin(), cand.end(), r.vs) != cand.end());
                        if (lex) CHECK(r.vs == *std::min_element(cand.begin(), cand.end()));
                    }
        }
    }
    std::cout << "[euler_path] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[euler_path] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        // ランダムな閉路を何本も重ねて、オイラー閉路があるグラフを作る
        Graph G(n);
        int v = 0;
        for (int k = 0; k < 2 * n; ++k) {
            const int u = (int)(rng() % n);
            G.add_edge(v, u), v = u;
        }
        G.add_edge(v, 0);
        auto t0 = std::chrono::steady_clock::now();
        const EulerPath r = euler_circuit_directed(G);
        const long long t1 = ms_since(t0);
        t0 = std::chrono::steady_clock::now();
        const EulerPath r2 = euler_path_undirected(G, -1, true);
        std::cout << "[euler_path] benchmark: N = " << n << ", M = " << G.edge_count() << "  有向 " << t1 << " ms, 無向 lex_min " << ms_since(t0) << " ms  (ok " << r.ok
                  << r2.ok << ")\n";
    }
}

}  // namespace ep_test


namespace ep_demo {
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
    banner("オイラー路（一筆書き） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto r = euler_path_undirected(G);", "EulerPath", "無向グラフの一筆書き"},
        {"auto r = euler_path_directed(G);", "EulerPath", "有向グラフの一筆書き"},
        {"euler_circuit_undirected(G) / _directed(G)", "EulerPath", "始点に戻ってくる一筆書き"},
        {"(G, s) / (G, -1, true)", "", "始点を s に / 辞書順最小"},
        {"r.ok / r.vs / r.es / r.closed()", "", "できたか / 頂点の列 / 辺番号の列 / 閉じているか"},
    };
    cheat_sheet("API 一覧", rows, 44, 14);

    head("デモで使うグラフ（無向）: 0-1, 1-2, 2-0, 2-3（次数が奇数なのは 2 と 3）");
    RUN(Graph G(4));
    RUN(G.add_undirected(0, 1); G.add_undirected(1, 2); G.add_undirected(2, 0); G.add_undirected(2, 3));

    api("【1】 euler_path_undirected(G)", "EulerPath", "すべての辺を 1 回ずつ通る道。次数が奇数の頂点が 2 つなら、そこが両端", "O(N + M)");
    RUN(auto r = euler_path_undirected(G));
    SHOW(r.ok, "一筆書きできる");
    SHOW(r.vs, "頂点の列");
    SHOW(r.es, "辺番号の列");
    SHOW(euler_path_undirected(G, 0).ok, "0 から始めることはできない（端は 2 か 3）");
    SHOW(euler_path_undirected(G, 3, true).vs, "3 から始める辞書順最小");
    SHOW(euler_circuit_undirected(G).ok, "奇数の頂点があるので閉路にはならない");

    api("【2】 euler_path_directed(G)", "EulerPath", "有向グラフ: 出次数 - 入次数 が全部 0 か、+1 と -1 が 1 つずつ", "O(N + M)");
    RUN(Graph D(3));
    RUN(D.add_edge(0, 1); D.add_edge(1, 2); D.add_edge(2, 0); D.add_edge(0, 2));
    SHOW(euler_path_directed(D).vs, "0 から出る辺が 1 本多いので 0 で始まり 2 で終わる");
    SHOW(euler_circuit_directed(D).ok, "");

    head("典型的な使い方");
    note("一筆書きできるか・その順番     : euler_path_undirected(G).ok / .vs");
    note("文字列のしりとりで全部つなぐ   : 文字を頂点、単語を辺（先頭の文字 → 末尾の文字）にして euler_path_directed");
    note("de Bruijn 列（全部の長さ k の文字列を含む最短の列）: 長さ k-1 の文字列を頂点にして euler_circuit_directed");
}

}  // namespace ep_demo

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
    ep_test::random_test();
    ep_test::benchmark();
    ep_demo::demo();
    return 0;
}
