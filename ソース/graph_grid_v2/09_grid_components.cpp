// =============================================================================
//  グリッドの連結成分（塗りつぶし）  grid_components / grid_components_same   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 09 グリッド
//  依存       : 00_graph_base, 00_grid_base, 01_components（同梱済み）
//  マスの番号 : id = i * W + j。結果 c は c.id[id]（成分番号。対象外のマスは -1）/ c.count / c.group(k) など
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  壁の指定は grid_bfs と同じ（'#' / "#X" / ~CharSet(".")）。「'#' の塊を数える」なら ~CharSet('#')（'#' 以外を壁にする）。
//  成分番号は行優先で最初に出てくるマスの成分から 0, 1, 2, ...（c.group(k) のマスは番号の昇順）。
//
//  機能一覧                                                                       計算量
//    grid_components(g, walls = '#', eight = false)   壁でないマスの連結成分（eight = true で 8 方向）   O(HW)
//    grid_components_same(g, eight = false)           同じ文字のマスどうしがつながる成分（全マスが対象）   O(HW)
//    （結果 Components c）c.count / c.id[id] / c.group(k) / c.size_of(k) / c.comp_size(id) / c.same(a, b)
//                          c.sizes() / c.max_size() / c.groups()
//    （隣接を自分で書くとき）components_core(g.size(), [&](int v, auto&& push) { ... }, alive)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤・グリッド基盤・連結成分も入っている）
//    gcc_test::random_test() / gcc_test::benchmark() / gcc_demo::demo()
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
#ifndef GG_COMPONENTS
#define GG_COMPONENTS
// ---- 連結成分（グラフ・グリッド第二世代）  本体 components_core は「頂点数 + 隣接を列挙する処理」を受け取る ----
struct Components {
    int count = 0;                  // 成分の数
    std::vector<int> id;            // 成分番号 0 ~ count-1（番号の小さい頂点を含む成分から順）。対象外の頂点は -1
    std::vector<int> start_, list_; // 成分ごとの頂点の一覧（CSR）。group(c) で使う
    IntSpan group(int c) const { return {list_.data() + start_[c], list_.data() + start_[c + 1]}; }  // 成分 c の頂点（昇順）
    int size_of(int c) const { return start_[c + 1] - start_[c]; }                                    // 成分 c の頂点数
    int comp_size(int v) const { return id[v] < 0 ? 0 : size_of(id[v]); }                             // v を含む成分の頂点数
    bool same(int u, int v) const { return id[u] >= 0 && id[u] == id[v]; }
    std::vector<int> sizes() const {
        std::vector<int> r(count);
        for (int c = 0; c < count; ++c) r[c] = size_of(c);
        return r;
    }
    int max_size() const {  // いちばん大きい成分の頂点数（成分がなければ 0）
        int m = 0;
        for (int c = 0; c < count; ++c) m = std::max(m, size_of(c));
        return m;
    }
    std::vector<std::vector<int>> groups() const {
        std::vector<std::vector<int>> r(count);
        for (int c = 0; c < count; ++c) r[c] = group(c).to_vector();
        return r;
    }
};
// adj(v, push): v の隣接 u ごとに push(u)。alive(v) が false の頂点は対象外（id = -1、たどらない）
template <class F, class P>
Components components_core(int n, F&& adj, P&& alive) {
    Components c;
    c.id.assign(n, -1);
    std::vector<int> q(n);  // 各頂点は高々 1 回しか入らない（生のポインタで回すほうが速い）
    int* const qp = q.data();
    int* const id = c.id.data();
    int qt = 0;
    for (int s = 0; s < n; ++s) {
        if (id[s] != -1 || !alive(s)) continue;
        const int k = c.count++;
        int h = qt;
        id[s] = k, qp[qt++] = s;
        for (; h < qt; ++h)
            adj(qp[h], [&](int u) {
                if (id[u] == -1 && alive(u)) id[u] = k, qp[qt++] = u;
            });
    }
    c.start_.assign(c.count + 1, 0);
    for (int v = 0; v < n; ++v)
        if (c.id[v] >= 0) ++c.start_[c.id[v] + 1];
    for (int k = 0; k < c.count; ++k) c.start_[k + 1] += c.start_[k];
    c.list_.assign(c.start_[c.count], 0);
    std::vector<int> pos(c.start_.begin(), c.start_.end() - 1);
    for (int v = 0; v < n; ++v)
        if (c.id[v] >= 0) c.list_[pos[c.id[v]]++] = v;
    return c;
}
template <class F>
Components components_core(int n, F&& adj) {
    return components_core(n, adj, [](int) { return true; });
}
// 無向グラフの連結成分（辺の向きのとおりにたどる。有向グラフの「向きを無視した成分」は weak_components）
inline Components components(const Graph& G) {
    return components_core(G.size(), [&](int v, auto&& push) {
        for (const Edge& e : G[v]) push(e.to);
    });
}
// 有向グラフの弱連結成分（辺の向きを無視してつながっている成分）
inline Components weak_components(const Graph& G) {
    const Graph R = G.reversed();
    return components_core(G.size(), [&](int v, auto&& push) {
        for (const Edge& e : G[v]) push(e.to);
        for (const Edge& e : R[v]) push(e.to);
    });
}
// 成分ごとの辺の本数（G.es で数える。多重辺・自己ループも 1 本。両端が対象外の辺は数えない）
inline std::vector<int> component_edge_counts(const Graph& G, const Components& c) {
    std::vector<int> r(c.count, 0);
    for (const RawEdge& e : G.es)
        if (c.id[e.from] >= 0) ++r[c.id[e.from]];
    return r;
}
#endif  // GG_COMPONENTS
#ifndef GG_GRID_COMPONENTS
#define GG_GRID_COMPONENTS
// ---- グリッドの連結成分（塗りつぶし）（グラフ・グリッド第二世代）  成分番号はマスの番号 id で引く ----
// 壁でないマスの連結成分（4 方向。eight = true なら 8 方向）。壁のマスは id = -1
inline Components grid_components(const Grid& g, const CharSet& walls = CharSet('#'), bool eight = false) {
    const std::vector<char> wall = g.mask(walls);
    auto alive = [&](int v) { return !wall[v]; };
    if (eight) return components_core(g.size(), [&](int v, auto&& push) { g.for_each8(v, push); }, alive);
    return components_core(g.size(), [&](int v, auto&& push) { g.for_each4(v, push); }, alive);
}
inline Components grid_components(const Grid& g, char wall, bool eight = false) { return grid_components(g, CharSet(wall), eight); }
inline Components grid_components(const Grid& g, const char* walls, bool eight = false) { return grid_components(g, CharSet(walls), eight); }  // "#X"（これがないと bool 版に行ってしまう）
Components grid_components(const Grid&, bool) = delete;  // 2 つ目は壁の文字（8 方向にするには grid_components(g, '#', true)）
// 同じ文字のマスどうしがつながる連結成分（色ごとの領域の数など）。全部のマスが対象
inline Components grid_components_same(const Grid& g, bool eight = false) {
    auto adj4 = [&](int v, auto&& push) {
        const char c = g.at(v);
        g.for_each4(v, [&](int u) {
            if (g.at(u) == c) push(u);
        });
    };
    auto adj8 = [&](int v, auto&& push) {
        const char c = g.at(v);
        g.for_each8(v, [&](int u) {
            if (g.at(u) == c) push(u);
        });
    };
    if (eight) return components_core(g.size(), adj8);
    return components_core(g.size(), adj4);
}
#endif  // GG_GRID_COMPONENTS
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "grid_components";
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

namespace gcc_test {

template <class T, class = void>
struct CallableBool : std::false_type {};
template <class T>
struct CallableBool<T, std::void_t<decltype(grid_components(std::declval<const Grid&>(), std::declval<T>()))>> : std::true_type {};
static_assert(!CallableBool<bool>::value, "grid_components(g, true) は使えない（8 方向は grid_components(g, '#', true)）");
static_assert(CallableBool<char>::value && CallableBool<const char*>::value && CallableBool<std::string>::value && CallableBool<CharSet>::value, "");

// 明示的なグラフで求めた成分を、対象のマスだけに番号を振り直したもの
static std::vector<int> expected_ids(const Grid& g, const std::function<bool(int)>& alive, const std::function<bool(int, int)>& link, bool eight, int& count) {
    Graph G(g.size());
    for (int v = 0; v < g.size(); ++v) {
        if (!alive(v)) continue;
        for (int d = 0; d < (eight ? 8 : 4); ++d) {
            const int u = g.neighbor(v, d);
            if (u != -1 && alive(u) && link(v, u) && v < u) G.add_undirected(v, u);
        }
    }
    const Components c = components(G);
    std::vector<int> re(c.count, -1), ids(g.size(), -1);
    count = 0;
    for (int v = 0; v < g.size(); ++v) {
        if (!alive(v)) continue;
        if (re[c.id[v]] == -1) re[c.id[v]] = count++;
        ids[v] = re[c.id[v]];
    }
    return ids;
}
static bool consistent(const Components& c) {
    const int n = (int)c.id.size();
    int total = 0;
    for (int k = 0; k < c.count; ++k) {
        const IntSpan gk = c.group(k);
        if (gk.size() != c.size_of(k) || gk.size() == 0) return false;
        for (int t = 0; t < gk.size(); ++t) {
            if (c.id[gk[t]] != k) return false;
            if (t > 0 && gk[t - 1] >= gk[t]) return false;
        }
        total += gk.size();
    }
    int alive = 0;
    for (int v = 0; v < n; ++v) alive += (c.id[v] >= 0);
    return total == alive;
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int H = rnd(rng, 1, 9), W = rnd(rng, 1, 9);
        Grid g(H, W);
        const char* chars = (it % 3 == 0) ? ".#" : (it % 3 == 1) ? "..#X" : "abc";
        const int k = (int)std::string(chars).size();
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) g.s[i][j] = chars[rng() % k];
        for (const char* wc : {"#", "#X", "a"}) {
            const CharSet walls(wc);
            auto alive = [&](int v) { return !walls(g.at(v)); };
            auto any = [](int, int) { return true; };
            for (bool eight : {false, true}) {
                int cnt = 0;
                const std::vector<int> ex = expected_ids(g, alive, any, eight, cnt);
                const Components c = grid_components(g, walls, eight);
                CHECK(c.count == cnt && c.id == ex && consistent(c));
                if (wc[1] == 0) {
                    const Components c2 = grid_components(g, wc[0], eight);
                    CHECK(c2.id == ex);
                }
            }
        }
        {
            const Components c = grid_components(g, "");  // 壁なし
            CHECK(c.count == 1 && c.size_of(0) == g.size() && consistent(c));
        }
        for (bool eight : {false, true}) {
            int cnt = 0;
            const std::vector<int> ex = expected_ids(g, [](int) { return true; }, [&](int a, int b) { return g.at(a) == g.at(b); }, eight, cnt);
            const Components c = eight ? grid_components_same(g, true) : grid_components_same(g);
            CHECK(c.count == cnt && c.id == ex && consistent(c));
        }
        // 既定の引数（'#'・4 方向）
        {
            int cnt = 0;
            const std::vector<int> ex = expected_ids(g, [&](int v) { return g.at(v) != '#'; }, [](int, int) { return true; }, false, cnt);
            CHECK(grid_components(g).id == ex);
        }
    }
    std::cout << "[grid_components] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[grid_components] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {1000, 2000}) {
        if (TEST_SCALE == 1 && n > 1000) break;
        Grid g(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) g.s[i][j] = "..#abc"[rng() % 6];
        auto t0 = std::chrono::steady_clock::now();
        const Components c = grid_components(g);
        const long long t1 = ms_since(t0);
        t0 = std::chrono::steady_clock::now();
        const Components c8 = grid_components(g, '#', true);
        const long long t2 = ms_since(t0);
        t0 = std::chrono::steady_clock::now();
        const Components cs = grid_components_same(g);
        std::cout << "[grid_components] benchmark: " << n << " × " << n << "  grid_components " << t1 << " ms, 8 方向 " << t2 << " ms, grid_components_same " << ms_since(t0)
                  << " ms  (count " << c.count << ", " << c8.count << ", " << cs.count << ")\n";
    }
}

}  // namespace gcc_test


namespace gcc_demo {
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

static void print_ids(const Grid& g, const Components& c) {
    for (int i = 0; i < g.H; ++i) {
        std::string line = "      ";
        for (int j = 0; j < g.W; ++j) {
            const int k = c.id[g.id(i, j)];
            const std::string s = k < 0 ? "#" : std::to_string(k);
            line += std::string(3 - s.size(), ' ') + s;
        }
        std::cout << line << "\n";
    }
}
static void print_grid(const Grid& g) {
    for (int i = 0; i < g.H; ++i) std::cout << "      " << g[i] << "\n";
}

static void demo() {
    banner("グリッドの連結成分（塗りつぶし） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto c = grid_components(g);", "Components", "'#' 以外のマスの成分（4 方向）"},
        {"grid_components(g, \"#X\") / (g, ~CharSet('#'))", "Components", "壁の文字を指定 / '#' の塊を数える"},
        {"grid_components(g, '#', true)", "Components", "8 方向"},
        {"grid_components_same(g, eight = false)", "Components", "同じ文字どうしがつながる成分"},
        {"c.count / c.id[id]", "int", "成分の数 / マスの成分番号（壁は -1）"},
        {"c.group(k) / c.size_of(k) / c.comp_size(id)", "IntSpan / int", "成分 k のマス / マスの数 / id の成分のマスの数"},
        {"c.same(a, b) / c.sizes() / c.max_size()", "", "同じ成分か / 全成分のマスの数 / 最大"},
    };
    cheat_sheet("API 一覧", rows, 48, 18);

    head("デモで使うグリッド");
    RUN(Grid g(std::vector<std::string>{"..#..", ".##.#", "#...#", "##.##"}));
    print_grid(g);

    api("【1】 grid_components(g, walls = '#')", "Components", "壁でないマスを上下左右でつないだ成分。番号は行優先で最初に出てくる成分から", "O(HW)");
    RUN(auto c = grid_components(g));
    SHOW(c.count, "島の数");
    note("各マスの成分番号（# は壁 = -1）:");
    print_ids(g, c);
    SHOW(c.sizes(), "各成分のマスの数");
    SHOW(c.max_size(), "いちばん大きい成分");
    SHOW(c.comp_size(g.id(0, 3)), "(0, 3) を含む成分のマスの数");
    SHOW(c.same(g.id(0, 0), g.id(3, 2)), "(0, 0) と (3, 2) は同じ成分か");
    SHOW(c.group(1), "成分 1 のマスの番号");
    note("for (int v : c.group(k)) { auto [i, j] = g.pos(v); ... } で成分ごとに処理できる。");

    api("【2】 grid_components(g, '#', true)", "Components", "8 方向（斜めもつながる）。ABC325 C（センサー）など", "O(HW)");
    RUN(auto c8 = grid_components(g, '#', true));
    SHOW(c8.count, "斜めでつながるので減る");
    print_ids(g, c8);
    note("grid_components(g, true) と書くとコンパイルエラーになる（2 つ目は壁の文字）。");

    api("【3】 grid_components(g, ~CharSet('#'))", "Components", "'#' 以外を壁にする = '#' の塊を数える", "O(HW)");
    RUN(auto cw = grid_components(g, ~CharSet('#')));
    SHOW(cw.count, "'#' の塊の数（4 方向）");
    SHOW(grid_components(g, ~CharSet('#'), true).count, "8 方向なら");

    api("【4】 grid_components_same(g, eight = false)", "Components", "同じ文字のマスどうしをつなぐ。全部のマスが対象（色ごとの領域）", "O(HW)");
    RUN(Grid p(std::vector<std::string>{"aab", "abb", "cca"}));
    print_grid(p);
    RUN(auto cs = grid_components_same(p));
    SHOW(cs.count, "領域の数");
    SHOW(cs.id, "");
    SHOW(p.at(cs.group(2)[0]), "成分 2 の文字");

    head("隣接を自分で書く（components_core に g.for_each4 などを渡す）");
    code("// 例: 数値のグリッドで、隣との差が 1 以下ならつながる。Grid shape(H, W) は形（番号の変換と隣接）のためだけに使う");
    RUN(std::vector<std::vector<int>> A = {{1, 2, 9}, {3, 2, 9}, {7, 7, 6}});
    RUN(Grid shape(3, 3));
    code("auto ca = components_core(shape.size(), [&](int v, auto&& push) {");
    code("    shape.for_each4(v, [&](int u) { if (std::abs(A[shape.row(u)][shape.col(u)] - A[shape.row(v)][shape.col(v)]) <= 1) push(u); });");
    code("});");
    auto ca = components_core(shape.size(), [&](int v, auto&& push) {
        shape.for_each4(v, [&](int u) {
            if (std::abs(A[shape.row(u)][shape.col(u)] - A[shape.row(v)][shape.col(v)]) <= 1) push(u);
        });
    });
    SHOW(ca.id, "{1, 2, 3, 2} / {9, 9} / {7, 7, 6}");
    note("マスを 1 つずつ足していく（後から道が増える）問題は 03_union_find の UnionFind で g.id(i, j) をつなぐ。");

    head("典型的な使い方");
    note("島・部屋の数             : grid_components(g).count");
    note("いちばん大きい領域       : grid_components(g).max_size()");
    note("2 マスが行き来できるか   : c.same(g.id(a), g.id(b))");
    note("塊ごとの集計（周の長さなど）: for (int v : c.group(k)) で各マスを見る");
}

}  // namespace gcc_demo

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
    gcc_test::random_test();
    gcc_test::benchmark();
    gcc_demo::demo();
    return 0;
}
