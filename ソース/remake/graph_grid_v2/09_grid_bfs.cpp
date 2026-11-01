// =============================================================================
//  グリッドの BFS（迷路の最短路）  grid_bfs / grid_bfs8 / grid_bfs_offsets / enclosed_cells
//                                                                 --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 09 グリッド
//  依存       : 00_graph_base, 00_grid_base, 01_bfs（同梱済み）
//  マスの番号 : id = i * W + j（g.id(i, j) / g.pos(id)）。結果 r は r.dist[id]（届かなければ DIST_INF）/ r.par[id] / r.path(t)
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  壁の指定（3 つ目の引数。省略すると '#'）: '#' / "#X"（どれかの文字）/ ~CharSet(".SG")（これ以外すべて）
//    ※ 3 つ目は「終点」ではない。grid_bfs(g, s, t) と書くとコンパイルエラーにしてある（終点の距離は r.dist[t]）
//    ※ 始点は壁のマスでもよい（そこから壁でないマスへ出られる）
//
//  機能一覧                                                                       計算量
//    grid_bfs(g, s, walls = '#')          4 方向の最短手数（s は id。{s1, s2, ...} で多始点）   O(HW)
//    grid_bfs8(g, s, walls = '#')         8 方向（斜めも 1 手）                                 O(HW)
//    grid_bfs_offsets(g, s, offs, walls)  移動リスト offs で 1 手（桂馬・ジャンプ。途中のマスは見ない）  O(HW |offs|)
//    enclosed_cells(g, walls = '#')       外周から入れない（壁に囲まれた）マス e[id] = 1          O(HW)
//    （隣接を自分で書くとき）bfs_core(g.size(), starts, [&](int v, auto&& push) { ... g.for_each4(v, ...) ... })
//      → ワープ・向きつき・氷の上を滑る・トーラスの例は使い方デモにある
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤・グリッド基盤・BFS も入っている）
//    gbfs_test::random_test() / gbfs_test::benchmark() / gbfs_demo::demo()
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
#ifndef GG_GRID_BASE
#define GG_GRID_BASE
// ---- グリッド基盤（グラフ・グリッド第二世代）  Grid / CharSet / 隣接 / 移動 / 入力 ----
struct CharSet {  // 文字の集合（壁の文字などを O(1) で判定）。'#' / "#X" / std::string から暗黙に作れる
    std::array<bool, 256> f{};
    CharSet() {}
    CharSet(char c) { f[(unsigned char)c] = true; }
    CharSet(const char* s) { for (; *s; ++s) f[(unsigned char)*s] = true; }
    CharSet(const std::string& s) { for (char c : s) f[(unsigned char)c] = true; }
    bool operator()(char c) const { return f[(unsigned char)c]; }
    CharSet operator~() const {  // 補集合（例: ~CharSet(".") =「'.' 以外すべて」）
        CharSet r;
        for (int i = 0; i < 256; ++i) r.f[i] = !f[i];
        return r;
    }
};

namespace gg_internal {
template <class F>
inline void call_cell(F& f, int u, int k) {  // f(u) と f(u, 方向 or 番号) の両方を受け付ける
    if constexpr (std::is_invocable_v<F&, int, int>) f(u, k);
    else f(u);
}
template <class F>
inline void call_ij(F& f, int i, int j, int k) {  // f(i, j) と f(i, j, 方向 or 番号) の両方を受け付ける
    if constexpr (std::is_invocable_v<F&, int, int, int>) f(i, j, k);
    else f(i, j);
}
}  // namespace gg_internal

struct Grid {
    int H = 0, W = 0;
    std::vector<std::string> s;  // s[i][j] = マス (i, j) の文字
    // 方向: 0 上, 1 右, 2 下, 3 左（ここまで 4 方向）, 4 左上, 5 右上, 6 右下, 7 左下（8 方向）
    static constexpr int DI[8] = {-1, 0, 1, 0, -1, -1, 1, 1};
    static constexpr int DJ[8] = {0, 1, 0, -1, -1, 1, 1, -1};
    static constexpr char DIR_CHAR[5] = "URDL";  // DIR_CHAR[d] = 方向 d の文字

    Grid() {}
    Grid(int h, int w, char fill = '.') : H(h), W(w), s(h, std::string(w, fill)) {}
    explicit Grid(std::vector<std::string> rows) : H((int)rows.size()), W(rows.empty() ? 0 : (int)rows[0].size()), s(std::move(rows)) {}

    // ---- 座標と番号（番号 id = i * W + j） ----
    int size() const { return H * W; }
    bool in(int i, int j) const { return (unsigned)i < (unsigned)H && (unsigned)j < (unsigned)W; }
    bool in(std::pair<int, int> p) const { return in(p.first, p.second); }
    int id(int i, int j) const { return i * W + j; }
    int id(std::pair<int, int> p) const { return p.first * W + p.second; }
    std::pair<int, int> pos(int v) const { return {v / W, v % W}; }
    int row(int v) const { return v / W; }
    int col(int v) const { return v % W; }
    char& operator()(int i, int j) { return s[i][j]; }
    char operator()(int i, int j) const { return s[i][j]; }
    char& at(int v) { return s[v / W][v % W]; }
    char at(int v) const { return s[v / W][v % W]; }
    std::string& operator[](int i) { return s[i]; }
    const std::string& operator[](int i) const { return s[i]; }

    // ---- 文字の検索 ----
    int count(const CharSet& c) const {
        int r = 0;
        for (const std::string& t : s)
            for (char x : t) r += c(x);
        return r;
    }
    int find(const CharSet& c) const {  // 行優先で最初のマスの番号（なければ -1）
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j)
                if (c(s[i][j])) return i * W + j;
        return -1;
    }
    std::vector<int> find_all(const CharSet& c) const {  // 当てはまるマスの番号すべて（行優先）
        std::vector<int> r;
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j)
                if (c(s[i][j])) r.push_back(i * W + j);
        return r;
    }
    std::vector<char> mask(const CharSet& c) const {  // m[id] = 1（当てはまる） / 0
        std::vector<char> m(size());
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) m[i * W + j] = c(s[i][j]);
        return m;
    }

    // ---- 隣接の列挙（f は f(u) か f(u, 方向) / (i, j) 版は f(ni, nj) か f(ni, nj, 方向)） ----
    int neighbor(int v, int d) const {  // 方向 d の隣のマス（範囲外なら -1）
        const int i = v / W + DI[d], j = v % W + DJ[d];
        return in(i, j) ? i * W + j : -1;
    }
    template <class F>
    void for_each4(int v, F&& f) const {
        const int i = v / W, j = v - i * W;
        if (i > 0) gg_internal::call_cell(f, v - W, 0);
        if (j + 1 < W) gg_internal::call_cell(f, v + 1, 1);
        if (i + 1 < H) gg_internal::call_cell(f, v + W, 2);
        if (j > 0) gg_internal::call_cell(f, v - 1, 3);
    }
    template <class F>
    void for_each8(int v, F&& f) const {
        const int i = v / W, j = v - i * W;
        for (int d = 0; d < 8; ++d) {
            const int ni = i + DI[d], nj = j + DJ[d];
            if (in(ni, nj)) gg_internal::call_cell(f, ni * W + nj, d);
        }
    }
    template <class F>
    void for_each4(int i, int j, F&& f) const {
        for (int d = 0; d < 4; ++d) {
            const int ni = i + DI[d], nj = j + DJ[d];
            if (in(ni, nj)) gg_internal::call_ij(f, ni, nj, d);
        }
    }
    template <class F>
    void for_each8(int i, int j, F&& f) const {
        for (int d = 0; d < 8; ++d) {
            const int ni = i + DI[d], nj = j + DJ[d];
            if (in(ni, nj)) gg_internal::call_ij(f, ni, nj, d);
        }
    }
    // 任意の移動リスト（桂馬跳び・決まった距離へのジャンプなど）。f(u) か f(u, k)（k = offs の何番目か）
    template <class F>
    void for_each_offsets(int v, const std::vector<std::pair<int, int>>& offs, F&& f) const {
        const int i = v / W, j = v - i * W;
        for (int k = 0; k < (int)offs.size(); ++k) {
            const int ni = i + offs[k].first, nj = j + offs[k].second;
            if (in(ni, nj)) gg_internal::call_cell(f, ni * W + nj, k);
        }
    }
    template <class F>
    void for_each_offsets(int i, int j, const std::vector<std::pair<int, int>>& offs, F&& f) const {
        for (int k = 0; k < (int)offs.size(); ++k) {
            const int ni = i + offs[k].first, nj = j + offs[k].second;
            if (in(ni, nj)) gg_internal::call_ij(f, ni, nj, k);
        }
    }
    // トーラス（上下・左右の端がつながっている）版
    template <class F>
    void for_each4_torus(int v, F&& f) const {
        const int i = v / W, j = v - i * W;
        for (int d = 0; d < 4; ++d) gg_internal::call_cell(f, wrap_id(i + DI[d], j + DJ[d]), d);
    }
    template <class F>
    void for_each8_torus(int v, F&& f) const {
        const int i = v / W, j = v - i * W;
        for (int d = 0; d < 8; ++d) gg_internal::call_cell(f, wrap_id(i + DI[d], j + DJ[d]), d);
    }
    template <class F>
    void for_each_offsets_torus(int v, const std::vector<std::pair<int, int>>& offs, F&& f) const {
        const int i = v / W, j = v - i * W;
        for (int k = 0; k < (int)offs.size(); ++k) gg_internal::call_cell(f, wrap_id(i + offs[k].first, j + offs[k].second), k);
    }
    // 方向 d に、壁かグリッドの外の手前まで滑ったマス（すぐ前が壁なら v のまま）
    int slide(int v, int d, const CharSet& walls) const {
        int i = v / W, j = v % W;
        while (true) {
            const int ni = i + DI[d], nj = j + DJ[d];
            if (!in(ni, nj) || walls(s[ni][nj])) break;
            i = ni, j = nj;
        }
        return i * W + j;
    }
    // トーラスの座標の正規化（負の数・H 以上でもよい）
    std::pair<int, int> wrap(int i, int j) const { return {(i % H + H) % H, (j % W + W) % W}; }
    int wrap_id(int i, int j) const { return ((i % H + H) % H) * W + (j % W + W) % W; }

    // ---- 形を変えたコピー・外周 ----
    Grid padded(char c = '#') const {  // 外周を c で 1 マスずつ囲んだ (H+2)×(W+2)。元の (i, j) は (i+1, j+1) になる
        Grid g(H + 2, W + 2, c);
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) g.s[i + 1][j + 1] = s[i][j];
        return g;
    }
    bool is_border(int v) const {
        const int i = v / W, j = v % W;
        return i == 0 || j == 0 || i == H - 1 || j == W - 1;
    }
    std::vector<int> border_cells() const {  // 外周のマス（行優先・重複なし）。外から塗る多始点 BFS の始点に
        std::vector<int> r;
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j)
                if (i == 0 || j == 0 || i == H - 1 || j == W - 1) r.push_back(i * W + j);
        return r;
    }

    // ---- 方向 ----
    // 経路（マスの番号の列）を "URDL" の文字列にする。同じ行・列の何マス先への移動も向きの文字になる。
    // 斜めなどは '?'。torus = true なら端をまたぐ 1 マスの移動も正しい向きにする
    std::string path_to_dirs(const std::vector<int>& path, bool torus = false) const {
        std::string r;
        for (size_t k = 0; k + 1 < path.size(); ++k) {
            int di = path[k + 1] / W - path[k] / W, dj = path[k + 1] % W - path[k] % W;
            if (torus) {
                if (H > 2 && di == H - 1) di = -1;
                else if (H > 2 && di == -(H - 1)) di = 1;
                if (W > 2 && dj == W - 1) dj = -1;
                else if (W > 2 && dj == -(W - 1)) dj = 1;
            }
            const int si = (di > 0) - (di < 0), sj = (dj > 0) - (dj < 0);
            r += (si == -1 && sj == 0) ? 'U' : (si == 0 && sj == 1) ? 'R' : (si == 1 && sj == 0) ? 'D' : (si == 0 && sj == -1) ? 'L' : '?';
        }
        return r;
    }
    static int dir_of(char c) {  // 'U','R','D','L'（小文字・矢印 ^ > v < も可）→ 0 ~ 3。それ以外は -1
        switch (c) {
            case 'U': case 'u': case '^': return 0;
            case 'R': case 'r': case '>': return 1;
            case 'D': case 'd': case 'v': return 2;
            case 'L': case 'l': case '<': return 3;
            default: return -1;
        }
    }
    static int reverse_dir(int d) { return d < 4 ? (d + 2) & 3 : ((d - 2) & 3) + 4; }  // 逆向き（8 方向にも対応）
    static int turn_right(int d) { return (d + 1) & 3; }                                // 時計回りに 90 度
    static int turn_left(int d) { return (d + 3) & 3; }                                 // 反時計回りに 90 度

    // ---- 移動リストを作る（for_each_offsets に渡す） ----
    static std::vector<std::pair<int, int>> offsets_manhattan(int D, bool with_center = false) {  // |di| + |dj| <= D
        std::vector<std::pair<int, int>> r;
        for (int di = -D; di <= D; ++di) {
            const int rem = D - (di < 0 ? -di : di);
            for (int dj = -rem; dj <= rem; ++dj)
                if (with_center || di != 0 || dj != 0) r.emplace_back(di, dj);
        }
        return r;
    }
    static std::vector<std::pair<int, int>> offsets_chebyshev(int D, bool with_center = false) {  // max(|di|, |dj|) <= D
        std::vector<std::pair<int, int>> r;
        for (int di = -D; di <= D; ++di)
            for (int dj = -D; dj <= D; ++dj)
                if (with_center || di != 0 || dj != 0) r.emplace_back(di, dj);
        return r;
    }
    static std::vector<std::pair<int, int>> offsets_euclid_sq(long long d2) {  // di^2 + dj^2 == d2 ちょうど（ABC272 D 型）
        std::vector<std::pair<int, int>> r;
        if (d2 < 0) return r;
        auto isqrt = [](long long x) {
            long long y = (long long)std::sqrt((double)x);
            while (y * y > x) --y;
            while ((y + 1) * (y + 1) <= x) ++y;
            return y;
        };
        const long long R = isqrt(d2);
        for (long long di = -R; di <= R; ++di) {
            const long long rem = d2 - di * di, dj = isqrt(rem);
            if (dj * dj != rem) continue;
            if (dj == 0) r.emplace_back((int)di, 0);
            else r.emplace_back((int)di, (int)-dj), r.emplace_back((int)di, (int)dj);
        }
        return r;
    }
    static std::vector<std::pair<int, int>> offsets_knight() {  // 桂馬（チェスのナイト）の 8 方向
        return {{-2, -1}, {-2, 1}, {-1, -2}, {-1, 2}, {1, -2}, {1, 2}, {2, -1}, {2, 1}};
    }
};

// ---- 入力 ----
inline Grid read_grid(int H, int W) {  // H 行の文字列を読む
    Grid g(H, W);
    for (std::string& r : g.s) std::cin >> r;
    return g;
}
inline Grid read_grid() {  // 「H W」に続いて H 行を読む
    int H, W;
    std::cin >> H >> W;
    return read_grid(H, W);
}
template <class T>
inline std::vector<std::vector<T>> read_matrix(int H, int W) {  // 数値のグリッド（空白区切り）。形は Grid g(H, W) で扱える
    std::vector<std::vector<T>> a(H, std::vector<T>(W));
    for (auto& r : a)
        for (T& x : r) std::cin >> x;
    return a;
}
inline std::ostream& operator<<(std::ostream& os, const Grid& g) {  // 各行を改行区切りで出力（最後の行の後ろには改行なし）
    for (int i = 0; i < g.H; ++i) os << (i ? "\n" : "") << g.s[i];
    return os;
}
#endif  // GG_GRID_BASE
#ifndef GG_BFS
#define GG_BFS
// ---- BFS（グラフ・グリッド第二世代）  本体 bfs_core は「頂点数 + 隣接を列挙する処理」を受け取る ----
// adj(v, push): v から 1 歩で行ける頂点 u ごとに push(u) を呼ぶ（到達済みかどうかは本体が調べる）。
// 通常のグラフ・グリッド・状態グラフ・辺を消したグラフなど、隣接の書き方を変えるだけで同じ本体を使える。
// target を指定すると、target の距離が決まった時点で打ち切る（ほかの頂点の距離は途中まで）。
template <class F>
SearchResult bfs_core(int n, const std::vector<int>& starts, F&& adj, int target = -1) {
    SearchResult r{std::vector<Cost>(n, DIST_INF), std::vector<int>(n, -1)};
    std::vector<int> q(n);  // 各頂点は高々 1 回しか入らないので n で足りる（生のポインタで回すほうが速い）
    int qt = 0;
    Cost* const dist = r.dist.data();
    int* const par = r.par.data();
    int* const qp = q.data();
    for (int s : starts)
        if (dist[s] != 0) dist[s] = 0, qp[qt++] = s;
    for (int h = 0; h < qt; ++h) {
        const int v = qp[h];
        if (v == target) break;
        const Cost nd = dist[v] + 1;
        adj(v, [&](int u) {
            if (dist[u] == DIST_INF) dist[u] = nd, par[u] = v, qp[qt++] = u;
        });
    }
    return r;
}
template <class F>
SearchResult bfs_core(int n, int s, F&& adj, int target = -1) { return bfs_core(n, std::vector<int>{s}, adj, target); }
// グラフの BFS（辺の重みは見ない。辺の本数での最短距離）
inline SearchResult bfs(const Graph& G, const std::vector<int>& starts, int target = -1) {
    return bfs_core(G.size(), starts, [&](int v, auto&& push) {
        for (const Edge& e : G[v]) push(e.to);
    }, target);
}
inline SearchResult bfs(const Graph& G, int s, int target = -1) { return bfs(G, std::vector<int>{s}, target); }
#endif  // GG_BFS
#ifndef GG_GRID_BFS
#define GG_GRID_BFS
// ---- グリッドの BFS（グラフ・グリッド第二世代）  結果 r はマスの番号 id = i * W + j で引く ----
// walls = 入れないマスの文字（'#'、"#X"、~CharSet(".") など。省略すると '#'）。始点は壁でもよい
inline SearchResult grid_bfs(const Grid& g, const std::vector<int>& starts, const CharSet& walls = CharSet('#')) {
    const std::vector<char> wall = g.mask(walls);
    return bfs_core(g.size(), starts, [&](int v, auto&& push) {
        g.for_each4(v, [&](int u) {
            if (!wall[u]) push(u);
        });
    });
}
inline SearchResult grid_bfs(const Grid& g, int s, const CharSet& walls = CharSet('#')) { return grid_bfs(g, std::vector<int>{s}, walls); }
inline SearchResult grid_bfs(const Grid& g, int s, char wall) { return grid_bfs(g, std::vector<int>{s}, CharSet(wall)); }
inline SearchResult grid_bfs(const Grid& g, const std::vector<int>& starts, char wall) { return grid_bfs(g, starts, CharSet(wall)); }
SearchResult grid_bfs(const Grid&, int, int) = delete;  // 3 つ目は壁の文字（終点ではない）
SearchResult grid_bfs(const Grid&, const std::vector<int>&, int) = delete;
// 8 方向に動ける BFS
inline SearchResult grid_bfs8(const Grid& g, const std::vector<int>& starts, const CharSet& walls = CharSet('#')) {
    const std::vector<char> wall = g.mask(walls);
    return bfs_core(g.size(), starts, [&](int v, auto&& push) {
        g.for_each8(v, [&](int u) {
            if (!wall[u]) push(u);
        });
    });
}
inline SearchResult grid_bfs8(const Grid& g, int s, const CharSet& walls = CharSet('#')) { return grid_bfs8(g, std::vector<int>{s}, walls); }
inline SearchResult grid_bfs8(const Grid& g, int s, char wall) { return grid_bfs8(g, std::vector<int>{s}, CharSet(wall)); }
inline SearchResult grid_bfs8(const Grid& g, const std::vector<int>& starts, char wall) { return grid_bfs8(g, starts, CharSet(wall)); }
SearchResult grid_bfs8(const Grid&, int, int) = delete;
SearchResult grid_bfs8(const Grid&, const std::vector<int>&, int) = delete;
// 任意の移動リストで動ける BFS（1 回の移動が 1 手。行き先が壁でなければよく、途中のマスは見ない。桂馬・ジャンプなど）
inline SearchResult grid_bfs_offsets(const Grid& g, const std::vector<int>& starts, const std::vector<std::pair<int, int>>& offs, const CharSet& walls = CharSet('#')) {
    const std::vector<char> wall = g.mask(walls);
    return bfs_core(g.size(), starts, [&](int v, auto&& push) {
        g.for_each_offsets(v, offs, [&](int u) {
            if (!wall[u]) push(u);
        });
    });
}
inline SearchResult grid_bfs_offsets(const Grid& g, int s, const std::vector<std::pair<int, int>>& offs, const CharSet& walls = CharSet('#')) {
    return grid_bfs_offsets(g, std::vector<int>{s}, offs, walls);
}
// 囲まれたマス: 壁でないマスのうち、外周の壁でないマスから（4 方向に、壁を通らずに）たどり着けないもの。e[id] = 1
inline std::vector<int> enclosed_cells(const Grid& g, const CharSet& walls = CharSet('#')) {
    std::vector<int> starts;
    for (int v : g.border_cells())
        if (!walls(g.at(v))) starts.push_back(v);
    const SearchResult r = grid_bfs(g, starts, walls);
    std::vector<int> e(g.size(), 0);
    for (int v = 0; v < g.size(); ++v) e[v] = (!walls(g.at(v)) && !r.reached(v));
    return e;
}
#endif  // GG_GRID_BFS
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "grid_bfs";
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

namespace gbfs_test {

// 3 つ目の引数が終点（int）だとコンパイルエラーになること（= delete）を確かめる
template <class T, class = void>
struct CallableBfs : std::false_type {};
template <class T>
struct CallableBfs<T, std::void_t<decltype(grid_bfs(std::declval<const Grid&>(), 0, std::declval<T>()))>> : std::true_type {};
template <class T, class = void>
struct CallableBfsV : std::false_type {};
template <class T>
struct CallableBfsV<T, std::void_t<decltype(grid_bfs8(std::declval<const Grid&>(), std::vector<int>{}, std::declval<T>()))>> : std::true_type {};
static_assert(!CallableBfs<int>::value && !CallableBfsV<int>::value, "grid_bfs(g, s, t) は使えない");
static_assert(CallableBfs<char>::value && CallableBfs<const char*>::value && CallableBfs<std::string>::value && CallableBfs<CharSet>::value, "");
static_assert(CallableBfsV<char>::value && CallableBfsV<const char*>::value, "");

static Grid random_grid(std::mt19937_64& rng, int H, int W, const char* chars) {
    Grid g(H, W);
    const int k = (int)std::string(chars).size();
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) g.s[i][j] = chars[rng() % k];
    return g;
}
// 明示的なグラフ: v → u（u が壁でない、u は offs で行ける）
static Graph grid_graph(const Grid& g, const CharSet& walls, const std::vector<std::pair<int, int>>& offs) {
    Graph G(g.size());
    for (int i = 0; i < g.H; ++i)
        for (int j = 0; j < g.W; ++j)
            for (auto [di, dj] : offs) {
                const int ni = i + di, nj = j + dj;
                if (ni < 0 || ni >= g.H || nj < 0 || nj >= g.W || walls(g.s[ni][nj])) continue;
                G.add_edge(g.id(i, j), g.id(ni, nj));
            }
    return G;
}
static bool same_dist_valid_par(const Grid& g, const SearchResult& r, const SearchResult& ex, const std::vector<int>& starts,
                                const std::vector<std::pair<int, int>>& offs, const CharSet& walls) {
    if (r.dist != ex.dist) return false;
    for (int v = 0; v < g.size(); ++v) {
        if (!r.reached(v)) continue;
        if (r.dist[v] == 0) {
            if (std::find(starts.begin(), starts.end(), v) == starts.end() || r.par[v] != -1) return false;
            continue;
        }
        const int p = r.par[v];
        if (p < 0 || r.dist[p] + 1 != r.dist[v] || walls(g.at(v))) return false;
        bool adj = false;
        for (auto [di, dj] : offs) adj |= (g.row(p) + di == g.row(v) && g.col(p) + dj == g.col(v));
        if (!adj) return false;
    }
    return true;
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    const std::vector<std::pair<int, int>> o4 = {{-1, 0}, {0, 1}, {1, 0}, {0, -1}};
    const std::vector<std::pair<int, int>> o8 = {{-1, 0}, {0, 1}, {1, 0}, {0, -1}, {-1, -1}, {-1, 1}, {1, 1}, {1, -1}};
    for (int it = 0; it < 2500 * TEST_SCALE; ++it, ++cases) {
        const int H = rnd(rng, 1, 8), W = rnd(rng, 1, 8);
        const Grid g = random_grid(rng, H, W, (it % 3 == 0) ? "..#" : ".#X.S");
        const int ns = rnd(rng, 1, 3);
        std::vector<int> starts;
        for (int k = 0; k < ns; ++k) starts.push_back(rnd(rng, 0, g.size() - 1));  // 壁でもよい・重複してもよい
        const int s = starts[0];
        // 4 方向
        {
            const Graph G = grid_graph(g, CharSet('#'), o4);
            const SearchResult ex1 = bfs(G, s), exm = bfs(G, starts);
            const SearchResult r1 = grid_bfs(g, s), r2 = grid_bfs(g, s, '#'), r3 = grid_bfs(g, starts), r4 = grid_bfs(g, starts, '#');
            CHECK(same_dist_valid_par(g, r1, ex1, {s}, o4, CharSet('#')) && r2.dist == ex1.dist && r2.par == r1.par);
            CHECK(same_dist_valid_par(g, r3, exm, starts, o4, CharSet('#')) && r4.dist == exm.dist);
            // 経路は 4 方向に 1 マスずつ・壁を通らない
            for (int t = 0; t < g.size(); ++t) {
                if (!r1.reached(t)) continue;
                const std::vector<int> path = r1.path(t);
                CHECK(path.front() == s && path.back() == t && (Cost)path.size() == r1.dist[t] + 1);
                const std::string dirs = g.path_to_dirs(path);
                CHECK((Cost)dirs.size() == r1.dist[t]);
            }
        }
        {
            const CharSet walls("#X");
            const Graph G = grid_graph(g, walls, o4);
            CHECK(same_dist_valid_par(g, grid_bfs(g, starts, "#X"), bfs(G, starts), starts, o4, walls));
            const CharSet only_dot = ~CharSet(".");
            const Graph G2 = grid_graph(g, only_dot, o4);
            CHECK(same_dist_valid_par(g, grid_bfs(g, s, ~CharSet(".")), bfs(G2, s), {s}, o4, only_dot));
        }
        // 8 方向
        {
            const Graph G = grid_graph(g, CharSet('#'), o8);
            CHECK(same_dist_valid_par(g, grid_bfs8(g, s), bfs(G, s), {s}, o8, CharSet('#')));
            CHECK(same_dist_valid_par(g, grid_bfs8(g, starts, '#'), bfs(G, starts), starts, o8, CharSet('#')));
            const Graph GX = grid_graph(g, CharSet("#X"), o8);
            CHECK(same_dist_valid_par(g, grid_bfs8(g, s, "#X"), bfs(GX, s), {s}, o8, CharSet("#X")));
        }
        // 移動リスト
        {
            std::vector<std::pair<int, int>> offs;
            if (it % 2 == 0) offs = Grid::offsets_knight();
            else
                for (int k = rnd(rng, 1, 5); k > 0; --k) offs.emplace_back(rnd(rng, -2, 2), rnd(rng, -2, 2));
            const Graph G = grid_graph(g, CharSet('#'), offs);
            CHECK(same_dist_valid_par(g, grid_bfs_offsets(g, s, offs), bfs(G, s), {s}, offs, CharSet('#')));
            const Graph GX = grid_graph(g, CharSet("#X"), offs);
            CHECK(same_dist_valid_par(g, grid_bfs_offsets(g, starts, offs, "#X"), bfs(GX, starts), starts, offs, CharSet("#X")));
        }
        // 囲まれたマス: 外周の壁でないマスから塗る（愚直: 何度も広げる）
        for (const char* wc : {"#", "#X"}) {
            const CharSet walls(wc);
            std::vector<int> out(g.size(), 0);
            for (int v = 0; v < g.size(); ++v) out[v] = (!walls(g.at(v)) && (g.row(v) == 0 || g.row(v) == H - 1 || g.col(v) == 0 || g.col(v) == W - 1));
            for (bool ch = true; ch;) {
                ch = false;
                for (int v = 0; v < g.size(); ++v) {
                    if (out[v] || walls(g.at(v))) continue;
                    for (auto [di, dj] : o4) {
                        const int ni = g.row(v) + di, nj = g.col(v) + dj;
                        if (g.in(ni, nj) && out[g.id(ni, nj)]) {
                            out[v] = 1, ch = true;
                            break;
                        }
                    }
                }
            }
            std::vector<int> e(g.size());
            for (int v = 0; v < g.size(); ++v) e[v] = (!walls(g.at(v)) && !out[v]);
            CHECK(enclosed_cells(g, walls) == e);
            if (wc[1] == 0) CHECK(enclosed_cells(g) == e);
        }
    }
    std::cout << "[grid_bfs] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[grid_bfs] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {1000, 2000}) {
        if (TEST_SCALE == 1 && n > 1000) break;
        for (int pct : {0, 25}) {
            Grid g(n, n);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) g.s[i][j] = ((int)(rng() % 100) < pct) ? '#' : '.';
            g.s[n / 2][n / 2] = '.';
            const int s = g.id(n / 2, n / 2);
            auto t0 = std::chrono::steady_clock::now();
            const SearchResult r = grid_bfs(g, s);
            const long long t4 = ms_since(t0);
            t0 = std::chrono::steady_clock::now();
            const SearchResult r8 = grid_bfs8(g, s);
            const long long t8 = ms_since(t0);
            t0 = std::chrono::steady_clock::now();
            const std::vector<int> e = enclosed_cells(g);
            const long long te = ms_since(t0);
            int reached = 0;
            for (int v = 0; v < g.size(); ++v) reached += r.reached(v);
            std::cout << "[grid_bfs] benchmark: " << n << " × " << n << "（壁 " << pct << "%、到達 " << reached << " マス）  grid_bfs " << t4 << " ms, grid_bfs8 " << t8
                      << " ms, enclosed_cells " << te << " ms  (check " << (r8.dist[0] + e[n + 1]) % 1000 << ")\n";
        }
    }
}

}  // namespace gbfs_test


namespace gbfs_demo {
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

// 距離をグリッドの形で表示する（壁は #、届かないマスは -）
static void print_dist(const Grid& g, const SearchResult& r, const CharSet& walls = CharSet('#')) {
    for (int i = 0; i < g.H; ++i) {
        std::string line = "      ";
        for (int j = 0; j < g.W; ++j) {
            const int v = g.id(i, j);
            std::string c = walls(g.at(v)) && r.dist[v] != 0 ? "#" : r.reached(v) ? std::to_string(r.dist[v]) : "-";
            line += std::string(3 - c.size(), ' ') + c;
        }
        std::cout << line << "\n";
    }
}
static void print_grid(const Grid& g) {
    for (int i = 0; i < g.H; ++i) std::cout << "      " << g[i] << "\n";
}

static void demo() {
    banner("グリッドの BFS 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto r = grid_bfs(g, s);", "SearchResult", "4 方向の最短手数（壁は '#'）。s = g.id(i, j)"},
        {"grid_bfs(g, s, \"#X\") / (g, s, ~CharSet(\".\"))", "SearchResult", "壁の文字を指定（どれか / これ以外）"},
        {"grid_bfs(g, {s1, s2}) / grid_bfs(g, g.find_all('S'))", "SearchResult", "多始点"},
        {"grid_bfs8(g, s, walls)", "SearchResult", "8 方向"},
        {"grid_bfs_offsets(g, s, offs, walls)", "SearchResult", "移動リストで 1 手（桂馬など）"},
        {"enclosed_cells(g, walls)", "vector<int>", "壁に囲まれて外周から入れないマス"},
        {"r.dist[t] / r.path(t) / g.path_to_dirs(path)", "", "距離（DIST_INF なら届かない）/ 経路 / \"RRD..\""},
        {"bfs_core(g.size(), s, [&](int v, auto&& push) {..})", "SearchResult", "隣接を自分で書く（ワープ・向き・滑る・トーラス）"},
    };
    cheat_sheet("API 一覧", rows, 54, 16);

    head("デモで使う迷路（入力の形: H W と H 行）");
    CinFromString cf("4 6\nS..#..\n.#.#.#\n.#...G\n...##.\n");
    RUN(Grid g = read_grid());
    print_grid(g);
    RUN(int s = g.find('S'), t = g.find('G'));

    api("【1】 grid_bfs(g, s)", "SearchResult", "S から各マスへの最短手数（上下左右に 1 マスずつ、'#' には入れない）", "O(HW)");
    RUN(auto r = grid_bfs(g, s));
    SHOW(r.dist[t], "G までの最短手数（届かなければ DIST_INF）");
    note("各マスの距離（# は壁、- は届かない）:");
    print_dist(g, r);
    SHOW(r.path(t), "経路（マスの番号）");
    SHOW(g.path_to_dirs(r.path(t)), "経路を UDLR の文字列に");
    SHOW(g.pos(r.path(t)[3]), "3 手目のマスの (i, j)");

    api("【2】 grid_bfs(g, s, 壁の文字)", "SearchResult", "'#' 以外を壁にするとき。\"#X\" ならどちらも壁、~CharSet(\".SG\") なら . S G 以外すべてが壁", "O(HW)");
    SHOW(grid_bfs(g, s, "#.").reached(t), "'.' も壁にすると届かない");
    SHOW(grid_bfs(g, s, ~CharSet(".SG")).dist[t], "「. S G 以外は壁」（この迷路では '#' と同じ）");
    note("grid_bfs(g, s, t) と書くとコンパイルエラーになる（3 つ目は壁の文字。終点の距離は r.dist[t]）。");

    api("【3】 grid_bfs(g, {s1, s2, ...})", "SearchResult", "多始点。各マスの「いちばん近い始点までの距離」。r.nearest_source() でどの始点か", "O(HW)");
    RUN(auto rm = grid_bfs(g, {g.id(0, 0), g.id(3, 5)}));
    print_dist(g, rm);
    note("始点を文字で集めるなら grid_bfs(g, g.find_all('S'))。");

    api("【4】 grid_bfs8(g, s)", "SearchResult", "斜めも 1 手（8 方向）", "O(HW)");
    SHOW(grid_bfs8(g, s).dist[t], "斜めに抜けられるので短くなる");

    api("【5】 grid_bfs_offsets(g, s, offs)", "SearchResult", "移動リスト offs（{di, dj} の一覧）のどれかで 1 手。行き先が壁でなければよい（途中は見ない）", "O(HW |offs|)");
    SHOW(grid_bfs_offsets(g, s, Grid::offsets_knight()).dist[t], "桂馬（ナイト）の動きで");
    SHOW(grid_bfs_offsets(g, s, Grid::offsets_manhattan(2)).dist[t], "マンハッタン距離 2 以内へジャンプ（壁を飛び越えられる）");

    api("【6】 enclosed_cells(g, walls = '#')", "vector<int>", "壁でないマスのうち、外周の壁でないマスから入れない（壁に囲まれた）もの。e[id] = 1", "O(HW)");
    RUN(Grid h(std::vector<std::string>{"#####.", "#..#..", "#####.", "......"}));
    print_grid(h);
    RUN(auto e = enclosed_cells(h));
    SHOW(e[h.id(1, 1)], "(1, 1) は囲まれている");
    SHOW(e[h.id(1, 4)], "(1, 4) は外とつながっている");
    SHOW(std::count(e.begin(), e.end(), 1), "囲まれたマスの数");

    head("隣接を自分で書く（bfs_core に g.for_each4 などを渡す）");
    code("// 例 1: ワープ（同じ英小文字のマスへ 1 手で移れる。ABC184 E）");
    code("//   文字ごとのマスの一覧を作り、その文字に初めて着いたときに一覧のマスを全部 push して一覧を空にする");
    RUN(Grid w(std::vector<std::string>{"Sa#.G", "###.#", "...a#"}));
    print_grid(w);
    code("std::vector<std::vector<int>> cells(26);");
    code("for (int v = 0; v < w.size(); ++v) if ('a' <= w.at(v) && w.at(v) <= 'z') cells[w.at(v) - 'a'].push_back(v);");
    code("auto rw = bfs_core(w.size(), w.find('S'), [&](int v, auto&& push) {");
    code("    w.for_each4(v, [&](int u) { if (w.at(u) != '#') push(u); });");
    code("    const char c = w.at(v);");
    code("    if ('a' <= c && c <= 'z') { for (int u : cells[c - 'a']) push(u); cells[c - 'a'].clear(); }");
    code("});");
    std::vector<std::vector<int>> cells(26);
    for (int v = 0; v < w.size(); ++v)
        if ('a' <= w.at(v) && w.at(v) <= 'z') cells[w.at(v) - 'a'].push_back(v);
    auto rw = bfs_core(w.size(), w.find('S'), [&](int v, auto&& push) {
        w.for_each4(v, [&](int u) {
            if (w.at(u) != '#') push(u);
        });
        const char c = w.at(v);
        if ('a' <= c && c <= 'z') {
            for (int u : cells[c - 'a']) push(u);
            cells[c - 'a'].clear();
        }
    });
    SHOW(rw.dist[w.find('G')], "S → a → (ワープ) a → 上 → 上 → 右");

    code("// 例 2: 向きを状態に持つ（前に 1 マス進む・右を向く・左を向く がどれも 1 手）。状態 x = マス * 4 + 向き");
    code("//   00_state_id の StateId S(g.size(), 4) を使うと S.id(v, d) / S.get(x) で書ける");
    RUN(Grid m(std::vector<std::string>{"S..", ".#.", "..G"}));
    code("auto rd = bfs_core(m.size() * 4, m.find('S') * 4 + 1, [&](int x, auto&& push) {   // 右向き（方向 1）で出発");
    code("    const int v = x / 4, d = x % 4, u = m.neighbor(v, d);");
    code("    if (u != -1 && m.at(u) != '#') push(u * 4 + d);     // 前に進む");
    code("    push(v * 4 + Grid::turn_right(d));                  // 右を向く");
    code("    push(v * 4 + Grid::turn_left(d));                   // 左を向く");
    code("});");
    auto rd = bfs_core(m.size() * 4, m.find('S') * 4 + 1, [&](int x, auto&& push) {
        const int v = x / 4, d = x % 4, u = m.neighbor(v, d);
        if (u != -1 && m.at(u) != '#') push(u * 4 + d);
        push(v * 4 + Grid::turn_right(d));
        push(v * 4 + Grid::turn_left(d));
    });
    RUN(Cost best = DIST_INF);
    RUN(for (int d = 0; d < 4; ++d) best = std::min(best, rd.dist[m.find('G') * 4 + d]));
    SHOW(best, "右向きで出発: 右 2 → 右を向く → 下 2");

    code("// 例 3: 氷の上を滑る（壁にぶつかるまで止まれない）。g.slide(v, 方向, 壁) が止まるマス");
    RUN(Grid ice(std::vector<std::string>{"S...#", ".#...", "...#.", "#...G"}));
    print_grid(ice);
    code("auto ri = bfs_core(ice.size(), ice.find('S'), [&](int v, auto&& push) { for (int d = 0; d < 4; ++d) push(ice.slide(v, d, '#')); });");
    auto ri = bfs_core(ice.size(), ice.find('S'), [&](int v, auto&& push) {
        for (int d = 0; d < 4; ++d) push(ice.slide(v, d, '#'));
    });
    SHOW(ri.dist[ice.find('G')], "G に止まれる最小の回数");
    SHOW(ice.path_to_dirs(ri.path(ice.find('G'))), "滑った向き（同じ行・列の何マス先でも向きの文字になる）");

    code("// 例 4: トーラス（上下・左右がつながっている）。g.for_each4_torus");
    RUN(Grid tor(std::vector<std::string>{"S.#.", "###.", "G.#."}));
    print_grid(tor);
    code("auto rt = bfs_core(tor.size(), tor.find('S'), [&](int v, auto&& push) {");
    code("    tor.for_each4_torus(v, [&](int u) { if (tor.at(u) != '#') push(u); });");
    code("});");
    auto rt = bfs_core(tor.size(), tor.find('S'), [&](int v, auto&& push) {
        tor.for_each4_torus(v, [&](int u) {
            if (tor.at(u) != '#') push(u);
        });
    });
    SHOW(rt.dist[tor.find('G')], "S から上にはみ出すと下の行の G（ふつうのグリッドでは届かない）");
    SHOW(tor.path_to_dirs(rt.path(tor.find('G')), true), "端をまたぐ移動も向きにする（torus = true）");

    head("典型的な使い方");
    note("迷路の最短手数          : grid_bfs(g, g.find('S')).dist[g.find('G')]（DIST_INF なら -1 を出力）");
    note("各マスから最も近い○○   : grid_bfs(g, g.find_all('○'))  →  r.dist[id]");
    note("囲まれた領域・外から塗る : enclosed_cells(g)、または grid_bfs(g, 外周の壁でないマス)");
    note("壁を壊す最小回数         : 09_grid_01bfs の grid_01bfs");
    note("塗りつぶし・領域の数     : 09_grid_components の grid_components");
    note("数値のグリッド（高さ差など）: Grid g(H, W) で形だけ作り、bfs_core の中で A[i][j] を見て push する");
}

}  // namespace gbfs_demo

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
    gbfs_test::random_test();
    gbfs_test::benchmark();
    gbfs_demo::demo();
    return 0;
}
