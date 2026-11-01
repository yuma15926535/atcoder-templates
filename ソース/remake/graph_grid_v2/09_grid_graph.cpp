// =============================================================================
//  グリッドをグラフにする  grid_graph / grid_graph_w / grid_graph_torus / grid_graph_offsets / grid_graph_directed   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 09 グリッド探索
//  依存       : 00_graph_base, 00_grid_base（同梱済み）
//  頂点番号   : 頂点 = マスの番号 g.id(i, j) = i * W + j（壁のマスも頂点として残る。辺がないだけ）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  グリッドを Graph にして、グラフの道具をそのまま使う:
//    最短路の本数（02_shortest_path_dag）・関節点 / 橋（07_lowlink）・閉路（01_cycle）・二部グラフ（01_bipartite）・
//    最小全域木（03_kruskal）・最大流 / マッチング（08_*）・連結成分（01_components）など。
//  最短距離・連結成分だけなら 09_grid_bfs / 09_grid_components のほうが速い（グラフを作らない）
//
//  機能一覧                                                                       計算量
//    grid_graph(g, walls = '#', eight = false)          壁でないマスの 4 方向（8 方向）の無向グラフ（重み 1）   O(HW)
//    grid_graph_w(g, cost, walls = '#', eight = false)  重みつきの無向グラフ（隣 u - v の重み cost(u, v)、u < v）  O(HW)
//    grid_graph_torus(g, walls = '#', eight = false)    トーラス版（端どうしもつながる）                        O(HW)
//    grid_graph_offsets(g, offs, walls = '#', torus = false)   任意の移動 (di, dj) の有向グラフ（重み 1）     O(HW |offs|)
//    grid_graph_directed(g, cost, walls = '#', eight = false)  隣への有向辺 u → v、重み cost(u, v)（DIST_INF なら張らない）  O(HW)
//    辺の番号: grid_graph / grid_graph_w は行優先で、各マスの「右・下（8 方向なら 右下・左下 も）」の順
//    grid_graph(g, true) はコンパイルエラー（2 つ目は壁の文字。8 方向は grid_graph(g, '#', true)）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤・グリッド基盤も入っている）
//    grg_test::random_test() / grg_test::benchmark() / grg_demo::demo()
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
#ifndef GG_GRID_GRAPH
#define GG_GRID_GRAPH
// ---- グリッドをグラフにする（グラフ・グリッド第二世代）  頂点 = マスの番号 g.id(i, j)（壁のマスは辺のない頂点として残る） ----
// グラフの道具（最短路の本数・lowlink・閉路・二部グラフ・最小全域木・フローなど）をグリッドにそのまま使うためのもの。
// 最短距離だけなら 09_grid_bfs のほうが速い（グラフを作らない）
namespace gg_internal {
template <class F>
Graph grid_graph_impl(const Grid& g, const CharSet& walls, bool eight, bool torus, F&& cost) {
    // 右・下（8 方向なら右下・左下も）への辺を 1 本ずつ張る。行優先の順
    static constexpr int FD[4][2] = {{0, 1}, {1, 0}, {1, 1}, {1, -1}};
    const int H = g.H, W = g.W, kinds = eight ? 4 : 2;
    std::vector<RawEdge> es;
    es.reserve((size_t)H * W * kinds);
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            if (walls(g.s[i][j])) continue;
            const int u = i * W + j;
            for (int k = 0; k < kinds; ++k) {
                int ni = i + FD[k][0], nj = j + FD[k][1];
                if (torus) {
                    ni = (ni % H + H) % H, nj = (nj % W + W) % W;
                } else if (!g.in(ni, nj)) {
                    continue;
                }
                const int v = ni * W + nj;
                if (v == u || walls(g.s[ni][nj])) continue;  // 自己ループは張らない（トーラスで H = 1 や W = 1 のとき）
                const Cost w = cost(u, v);
                if (w != DIST_INF) es.push_back({u, v, w});
            }
        }
    return build_graph(H * W, std::move(es), UNDIRECTED);
}
struct UnitCost {
    Cost operator()(int, int) const { return 1; }
};
}  // namespace gg_internal
// 壁でないマスどうしを 4 方向（eight なら 8 方向）で結んだ無向グラフ（重み 1）
inline Graph grid_graph(const Grid& g, const CharSet& walls = CharSet('#'), bool eight = false) {
    return gg_internal::grid_graph_impl(g, walls, eight, false, gg_internal::UnitCost());
}
Graph grid_graph(const Grid&, bool) = delete;  // 2 つ目は壁の文字（8 方向にするには grid_graph(g, '#', true)）
// トーラス版（上下・左右の端もつながる）。H や W が 2 以下だと同じ 2 マスの間に辺が 2 本できる（for_each4_torus と同じ数え方）
inline Graph grid_graph_torus(const Grid& g, const CharSet& walls = CharSet('#'), bool eight = false) {
    return gg_internal::grid_graph_impl(g, walls, eight, true, gg_internal::UnitCost());
}
Graph grid_graph_torus(const Grid&, bool) = delete;
// 重みつきの無向グラフ: 隣り合う壁でないマス u, v（u < v。トーラスでないとき）の辺の重みを cost(u, v) にする（DIST_INF なら張らない）。
// 「隣との高さの差」の最小全域木・最大値の最小化（03_kruskal / 03_union_find）など。辺の順は grid_graph と同じ
template <class F>
Graph grid_graph_w(const Grid& g, F&& cost, const CharSet& walls = CharSet('#'), bool eight = false) {
    return gg_internal::grid_graph_impl(g, walls, eight, false, cost);
}
// 任意の移動（offs の各 (di, dj)）で、壁でないマス u から壁でないマス u + (di, dj) への有向グラフ（重み 1）。torus なら端をまたいでよい
inline Graph grid_graph_offsets(const Grid& g, const std::vector<std::pair<int, int>>& offs, const CharSet& walls = CharSet('#'), bool torus = false) {
    const int H = g.H, W = g.W;
    std::vector<RawEdge> es;
    es.reserve((size_t)H * W * offs.size());
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            if (walls(g.s[i][j])) continue;
            for (const auto& o : offs) {
                int ni = i + o.first, nj = j + o.second;
                if (torus) {
                    ni = (ni % H + H) % H, nj = (nj % W + W) % W;
                } else if (!g.in(ni, nj)) {
                    continue;
                }
                if (walls(g.s[ni][nj])) continue;
                es.push_back({i * W + j, ni * W + nj, 1});
            }
        }
    return gg_internal::build_graph(H * W, std::move(es), DIRECTED);
}
// 重みつきの有向グラフ: 壁でないマス u から 4 方向（eight なら 8 方向）の壁でないマス v へ、重み cost(u, v) の辺
// （cost が DIST_INF なら辺を張らない）。「入るマスの数字がコスト」なら [&](int u, int v) { return A[v / W][v % W]; }
template <class F>
Graph grid_graph_directed(const Grid& g, F&& cost, const CharSet& walls = CharSet('#'), bool eight = false) {
    const int H = g.H, W = g.W;
    std::vector<RawEdge> es;
    es.reserve((size_t)H * W * (eight ? 8 : 4));
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            if (walls(g.s[i][j])) continue;
            const int u = i * W + j;
            for (int d = 0; d < (eight ? 8 : 4); ++d) {
                const int ni = i + Grid::DI[d], nj = j + Grid::DJ[d];
                if (!g.in(ni, nj) || walls(g.s[ni][nj])) continue;
                const int v = ni * W + nj;
                const Cost w = cost(u, v);
                if (w != DIST_INF) es.push_back({u, v, w});
            }
        }
    return gg_internal::build_graph(H * W, std::move(es), DIRECTED);
}
#endif  // GG_GRID_GRAPH
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "grid_graph";
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

namespace grg_test {

static Grid random_grid(std::mt19937_64& rng, int H, int W, int wall_pct) {
    Grid g(H, W);
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) g(i, j) = (rnd(rng, 0, 99) < wall_pct) ? '#' : (rnd(rng, 0, 9) == 0 ? 'X' : '.');
    return g;
}
// G[v] の行き先の一覧（ソート済み）
static std::vector<int> adj_sorted(const Graph& G, int v) {
    std::vector<int> a;
    for (const Edge& e : G[v]) a.push_back(e.to);
    std::sort(a.begin(), a.end());
    return a;
}
static std::vector<Cost> bfs_dist(const Graph& G, int s) {
    std::vector<Cost> d(G.size(), DIST_INF);
    std::vector<int> q{s};
    d[s] = 0;
    for (size_t h = 0; h < q.size(); ++h)
        for (const Edge& e : G[q[h]])
            if (d[e.to] == DIST_INF) d[e.to] = d[q[h]] + 1, q.push_back(e.to);
    return d;
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int H = rnd(rng, 1, 6), W = rnd(rng, 1, 6);
        const Grid g = random_grid(rng, H, W, rnd(rng, 0, 60));
        const CharSet walls = (it % 3 == 0) ? CharSet("#X") : CharSet('#');
        for (int eight = 0; eight < 2; ++eight) {
            // (1) 通常: 隣接が for_each4 / for_each8 と一致、辺の両端は壁でない、辺の本数 = 隣り合う組の数
            {
                const Graph G = grid_graph(g, walls, eight);
                CHECK(G.size() == H * W);
                int pairs = 0;
                for (int v = 0; v < H * W; ++v) {
                    std::vector<int> exp;
                    if (!walls(g.at(v))) {
                        auto f = [&](int u) { if (!walls(g.at(u))) exp.push_back(u); };
                        if (eight) g.for_each8(v, f);
                        else g.for_each4(v, f);
                    }
                    std::sort(exp.begin(), exp.end());
                    CHECK(adj_sorted(G, v) == exp);
                    pairs += (int)exp.size();
                }
                CHECK(G.edge_count() * 2 == pairs);
                for (const RawEdge& e : G.es) CHECK(!walls(g.at(e.from)) && !walls(g.at(e.to)) && e.w == 1 && e.from < e.to);
                for (int id = 1; id < G.edge_count(); ++id) CHECK(G.es[id - 1].from <= G.es[id].from);  // 行優先
                // BFS の距離が、グリッドの上で直接 BFS したものと一致
                const int s = rnd(rng, 0, H * W - 1);
                std::vector<Cost> d(H * W, DIST_INF);
                if (!walls(g.at(s))) {
                    std::vector<int> q{s};
                    d[s] = 0;
                    for (size_t h = 0; h < q.size(); ++h) {
                        auto f = [&](int u) { if (!walls(g.at(u)) && d[u] == DIST_INF) d[u] = d[q[h]] + 1, q.push_back(u); };
                        if (eight) g.for_each8(q[h], f);
                        else g.for_each4(q[h], f);
                    }
                } else {
                    d[s] = 0;
                }
                CHECK(bfs_dist(G, s) == d);
            }
            // (1') 重みつき無向: grid_graph と同じ辺（DIST_INF の辺を除く）で、重みが cost(u, v)
            {
                const Graph G0 = grid_graph(g, walls, eight);
                auto cost = [&](int u, int v) { return ((u * 7 + v * 3) % 5 == 0) ? DIST_INF : Cost((u * 31 + v * 17) % 10); };
                const Graph G = grid_graph_w(g, cost, walls, eight);
                std::vector<RawEdge> exp;
                for (const RawEdge& e : G0.es)
                    if (cost(e.from, e.to) != DIST_INF) exp.push_back({e.from, e.to, cost(e.from, e.to)});
                CHECK(G.edge_count() == (int)exp.size());
                for (int id = 0; id < G.edge_count(); ++id) CHECK(G.es[id].from == exp[id].from && G.es[id].to == exp[id].to && G.es[id].w == exp[id].w);
            }
            // (2) トーラス: 自分以外への隣接が for_each4_torus / for_each8_torus と一致（重複も数える）
            {
                const Graph G = grid_graph_torus(g, walls, eight);
                for (int v = 0; v < H * W; ++v) {
                    std::vector<int> exp;
                    if (!walls(g.at(v))) {
                        auto f = [&](int u) { if (u != v && !walls(g.at(u))) exp.push_back(u); };
                        if (eight) g.for_each8_torus(v, f);
                        else g.for_each4_torus(v, f);
                    }
                    std::sort(exp.begin(), exp.end());
                    if (eight == 0 || (H >= 3 && W >= 3)) CHECK(adj_sorted(G, v) == exp);  // 8 方向で H, W ≤ 2 は同じ組が何重にもなるので数え方が違う
                    else {
                        std::vector<int> a = adj_sorted(G, v);
                        a.erase(std::unique(a.begin(), a.end()), a.end());
                        exp.erase(std::unique(exp.begin(), exp.end()), exp.end());
                        CHECK(a == exp);  // 隣かどうかは一致
                    }
                }
                for (const RawEdge& e : G.es) CHECK(e.from != e.to);
            }
        }
        // (3) 任意の移動（有向）
        {
            std::vector<std::pair<int, int>> offs;
            const int k = rnd(rng, 0, 5);
            for (int t = 0; t < k; ++t) offs.push_back({rnd(rng, -3, 3), rnd(rng, -3, 3)});
            for (int torus = 0; torus < 2; ++torus) {
                const Graph G = grid_graph_offsets(g, offs, walls, torus);
                for (int v = 0; v < H * W; ++v) {
                    std::vector<int> exp;
                    if (!walls(g.at(v))) {
                        auto f = [&](int u) { if (!walls(g.at(u))) exp.push_back(u); };
                        if (torus) g.for_each_offsets_torus(v, offs, f);
                        else g.for_each_offsets(v, offs, f);
                    }
                    std::sort(exp.begin(), exp.end());
                    CHECK(adj_sorted(G, v) == exp);
                }
            }
        }
        // (4) 重みつき有向（入るマスの数字がコスト。'X' のマスには入れない）
        {
            std::vector<Cost> a(H * W);
            for (auto& x : a) x = rnd(rng, 0, 9);
            for (int eight = 0; eight < 2; ++eight) {
                const Graph G = grid_graph_directed(g, [&](int, int v) { return g.at(v) == 'X' ? DIST_INF : a[v]; }, '#', eight);
                for (int v = 0; v < H * W; ++v) {
                    std::vector<std::pair<int, Cost>> exp, got;
                    if (g.at(v) != '#') {
                        auto f = [&](int u) { if (g.at(u) == '.') exp.push_back({u, a[u]}); };
                        if (eight) g.for_each8(v, f);
                        else g.for_each4(v, f);
                    }
                    for (const Edge& e : G[v]) got.push_back({e.to, e.w});
                    std::sort(exp.begin(), exp.end()), std::sort(got.begin(), got.end());
                    CHECK(got == exp);
                }
            }
        }
    }
    std::cout << "[grid_graph] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[grid_graph] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    const int H = (TEST_SCALE == 1 ? 500 : 2000), W = H;
    const Grid g = random_grid(rng, H, W, 20);
    auto t0 = std::chrono::steady_clock::now();
    const Graph G = grid_graph(g);
    const long long t1 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const Graph G8 = grid_graph(g, '#', true);
    const long long t2 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<Cost> d = bfs_dist(G, 0);
    std::cout << "[grid_graph] benchmark: " << H << " × " << W << "  grid_graph " << t1 << " ms（辺 " << G.edge_count() << "）, 8 方向 " << t2 << " ms（辺 "
              << G8.edge_count() << "）, その上の BFS " << ms_since(t0) << " ms  (check " << d[H * W - 1] % 1000 << ")\n";
}

}  // namespace grg_test


namespace grg_demo {
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

static void print_grid(const Grid& g) {
    for (int i = 0; i < g.H; ++i) std::cout << "      " << g[i] << "\n";
}

static void demo() {
    banner("グリッドをグラフにする 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"Graph G = grid_graph(g);", "Graph", "壁 '#' でないマスの 4 方向の無向グラフ（頂点 = g.id(i, j)）"},
        {"grid_graph(g, \"#X\", true)", "Graph", "壁の文字を複数・8 方向"},
        {"grid_graph_w(g, cost, walls, eight)", "Graph", "隣 u - v の重みが cost(u, v) の無向グラフ（高さの差など）"},
        {"grid_graph_torus(g, walls, eight)", "Graph", "端どうしもつながる"},
        {"grid_graph_offsets(g, offs, walls, torus)", "Graph", "任意の移動 (di, dj) の有向グラフ（桂馬跳びなど）"},
        {"grid_graph_directed(g, cost, walls, eight)", "Graph", "u → v の重み cost(u, v)（入るマスのコストなど）"},
    };
    cheat_sheet("API 一覧", rows, 46, 12);

    head("デモで使うグリッド（3 × 4）");
    RUN(Grid g(std::vector<std::string>{"..#.", "....", ".#.."}));
    print_grid(g);

    api("【1】 grid_graph(g)", "Graph", "壁でない隣どうしを無向辺で結ぶ。壁のマスは辺のない頂点", "O(HW)");
    RUN(Graph G = grid_graph(g));
    SHOW(G.size(), "頂点数 = H × W");
    SHOW(G.edge_count(), "隣り合う（壁でない）組の数");
    SHOW(G[g.id(1, 1)], "マス (1, 1) の隣接（行き先・辺番号・重み 1）");
    SHOW(G[g.id(0, 2)].size(), "壁のマスには辺がない");

    api("【2】 grid_graph_directed(g, cost)", "Graph", "入るマスの数字がコストの有向グラフ（02_dijkstra に渡す）", "O(HW)");
    RUN(std::vector<std::vector<int>> A = {{1, 2, 0, 3}, {4, 1, 1, 1}, {9, 0, 5, 1}});
    RUN(Graph D = grid_graph_directed(g, [&](int, int v) { return (Cost)A[v / 4][v % 4]; }));
    SHOW(D[g.id(1, 1)], "(1, 1) から入れる隣と、そのマスの数字");

    api("【3】 grid_graph_w(g, cost)", "Graph", "隣どうしの重みが「数字の差の絶対値」の無向グラフ（最小全域木・最大値の最小化へ）", "O(HW)");
    RUN(Graph Wg = grid_graph_w(g, [&](int u, int v) { return (Cost)std::abs(A[u / 4][u % 4] - A[v / 4][v % 4]); }));
    SHOW(Wg[g.id(1, 1)], "(1, 1) の隣との差");

    api("【4】 grid_graph_offsets(g, Grid::offsets_knight())", "Graph", "桂馬跳びの有向グラフ", "O(HW × 8)");
    RUN(Graph K = grid_graph_offsets(g, Grid::offsets_knight()));
    SHOW(K[g.id(0, 0)], "(0, 0) から跳べるマス");

    head("ほかのテンプレートと組み合わせる（そのテンプレートを一緒に貼る）");
    code("auto r = bfs(G, s);  auto cnt = count_shortest_paths(G, r, MOD);   // 最短路の本数（02_shortest_path_dag）");
    code("LowLink L(G);  L.is_articulation[g.id(i, j)]                      // 消すと分かれるマス（07_lowlink）");
    code("Cycle c = find_cycle_undirected(G);                                // マスの閉路（01_cycle）");
    code("Bipartite b = bipartite(G);  // 市松模様の色。ドミノの最大数は 08_bipartite_matching へ");
    code("auto r = dijkstra(grid_graph_directed(g, cost), s);                // マスのコストの最短路（02_dijkstra）");
    code("MstResult t = kruskal(grid_graph_w(g, cost));                      // 隣どうしの重みの最小全域木（03_kruskal）");
    note("最短距離だけなら 09_grid_bfs（グラフを作らないぶん速い）。連結成分だけなら 09_grid_components。");
}

}  // namespace grg_demo

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
    grg_test::random_test();
    grg_test::benchmark();
    grg_demo::demo();
    return 0;
}
