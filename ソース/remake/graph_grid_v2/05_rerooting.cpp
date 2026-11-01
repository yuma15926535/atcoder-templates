// =============================================================================
//  全方位木DP  rerooting                                             --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 木
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed。G は木（森でもよい。成分ごとに計算する）。read_tree(N) で読んだグラフをそのまま渡す
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  「すべての頂点 v について、v を根にしたときの木DP の値」を O(N) で求める。木DP の形を 3 つの関数で書く:
//    dp[v] = put_vertex( merge( put_edge(dp[c1], 辺 v→c1), put_edge(dp[c2], 辺 v→c2), ... ), v )
//    merge(a, b)         子からの値をまとめる（結合的・可換）。e はその単位元（何もないとき）
//    put_edge(x, e)      子の側の値 x を、辺 e（e.to が子、e.w が重み）を通して持ち上げる
//    put_vertex(x, v)    まとめた値 x に頂点 v 自身の分を足す
//  （旧テンプレートの「最大距離」「距離の和」は、この 3 つを変えるだけで書ける。使い方デモを参照）
//
//  機能一覧                                                                       計算量
//    rerooting(G, e, merge, put_edge, put_vertex)   ans[v] = v を根にしたときの dp[v]    O(N)（merge などの呼び出し回数）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    rr_test::random_test() / rr_test::benchmark() / rr_demo::demo()
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
#ifndef GG_REROOTING
#define GG_REROOTING
// ---- 全方位木DP（グラフ・グリッド第二世代）  すべての頂点について「その頂点を根にしたときの木DP の値」を O(N) で求める（非再帰） ----
// 木DP の形（頂点 v を根とする部分木の値 dp[v]）:
//   dp[v] = put_vertex( merge( put_edge(dp[c1], 辺 v→c1), put_edge(dp[c2], 辺 v→c2), ... ), v )
//   merge(a, b): 子からの値をまとめる（結合的・可換。e は単位元）
//   put_edge(x, e): 辺 e（G[v] の要素。e.to が子の側）を通して、子の側の値 x を v から見た値にする（重みを足すなど）
//   put_vertex(x, v): まとめた値 x に頂点 v 自身の分を足す（葉なら x = e）
// 返り値 ans[v] = 頂点 v を根にしたときの dp[v]。森にも対応（成分ごとに計算する）
template <class T, class Merge, class PutEdge, class PutVertex>
std::vector<T> rerooting(const Graph& G, const T& e, Merge merge, PutEdge put_edge, PutVertex put_vertex) {
    const int n = G.size();
    std::vector<int> order, par(n, -1), pidx(n, -1);  // pidx[v] = G[v] での親への辺の位置
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
                seen[u] = 1, par[u] = v, order.push_back(u);
            }
        }
    }
    for (int v = 0; v < n; ++v) {
        if (par[v] < 0) continue;
        for (int i = 0; i < (int)G[v].size(); ++i)
            if (G[v][i].to == par[v]) {
                pidx[v] = i;
                break;
            }
    }
    // 1 回目: 葉から根へ（部分木の値）
    std::vector<T> down(n, e);
    for (int k = n - 1; k >= 0; --k) {
        const int v = order[k];
        T acc = e;
        for (int i = 0; i < (int)G[v].size(); ++i)
            if (i != pidx[v]) acc = merge(acc, put_edge(down[G[v][i].to], G[v][i]));
        down[v] = put_vertex(acc, v);
    }
    // 2 回目: 根から葉へ（親の側から来る値を配る）。左右からの累積で「自分以外」をまとめる
    std::vector<T> from_par(n, e), ans(n, e), suf, c;
    for (int k = 0; k < n; ++k) {
        const int v = order[k];
        const int d = (int)G[v].size();
        suf.assign(d + 1, e), c.assign(d, e);
        for (int i = 0; i < d; ++i) c[i] = (i == pidx[v]) ? from_par[v] : put_edge(down[G[v][i].to], G[v][i]);
        for (int i = d - 1; i >= 0; --i) suf[i] = merge(c[i], suf[i + 1]);
        ans[v] = put_vertex(suf[0], v);
        T pre = e;
        for (int i = 0; i < d; ++i) {
            if (i != pidx[v]) {
                const int u = G[v][i].to;
                from_par[u] = put_edge(put_vertex(merge(pre, suf[i + 1]), v), G[u][pidx[u]]);
            }
            pre = merge(pre, c[i]);
        }
    }
    return ans;
}
#endif  // GG_REROOTING
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "rerooting";
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

namespace rr_test {

// 愚直: 根を決めて再帰で木DP（小さい木だけ）
template <class T, class M, class PE, class PV>
static T brute_dp(const Graph& G, int v, int p, const T& e, M merge, PE put_edge, PV put_vertex) {
    T acc = e;
    for (const Edge& ed : G[v])
        if (ed.to != p) acc = merge(acc, put_edge(brute_dp(G, ed.to, v, e, merge, put_edge, put_vertex), ed));
    return put_vertex(acc, v);
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    const long long MOD = 998244353;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 12);
        Graph G(n);
        std::vector<int> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        for (int i = 1; i < n; ++i) {
            if (it % 5 == 0 && rng() % 4 == 0) continue;  // 森
            G.add_undirected(perm[i], perm[rnd(rng, 0, i - 1)], rnd(rng, 0, 9));
        }
        std::vector<long long> val(n);
        for (auto& x : val) x = rnd(rng, 0, 9);
        // (1) いちばん遠い頂点までの距離（重みつき）
        {
            auto merge = [](long long a, long long b) { return std::max(a, b); };
            auto pe = [](long long x, const Edge& e) { return x + e.w; };
            auto pv = [](long long x, int) { return x; };
            const std::vector<long long> ans = rerooting(G, 0LL, merge, pe, pv);
            for (int v = 0; v < n; ++v) CHECK(ans[v] == brute_dp(G, v, -1, 0LL, merge, pe, pv));
        }
        // (2) 全頂点への距離の和（{和, 頂点数}）
        {
            using P = std::pair<long long, long long>;
            auto merge = [](P a, P b) { return P{a.first + b.first, a.second + b.second}; };
            auto pe = [](P x, const Edge& e) { return P{x.first + x.second * e.w, x.second}; };
            auto pv = [](P x, int) { return P{x.first, x.second + 1}; };
            const std::vector<P> ans = rerooting(G, P{0, 0}, merge, pe, pv);
            for (int v = 0; v < n; ++v) CHECK(ans[v] == brute_dp(G, v, -1, P{0, 0}, merge, pe, pv));
        }
        // (3) v を含む連結な部分集合の数（EDPC V の形。mod）と頂点の値
        {
            auto merge = [&](long long a, long long b) { return a * b % MOD; };
            auto pe = [&](long long x, const Edge&) { return (x + 1) % MOD; };
            auto pv = [&](long long x, int v) { return x * (val[v] + 1) % MOD; };
            const std::vector<long long> ans = rerooting(G, 1LL, merge, pe, pv);
            for (int v = 0; v < n; ++v) CHECK(ans[v] == brute_dp(G, v, -1, 1LL, merge, pe, pv));
        }
        // (4) 頂点の値の最大値と、そこまでの距離の組（辺ごとに違う型の値を返す put_edge）
        {
            using P = std::pair<long long, long long>;  // {値の最大, 距離の最小（その値をとる頂点まで）}
            auto better = [](P a, P b) { return (a.first != b.first) ? (a.first > b.first ? a : b) : (a.second < b.second ? a : b); };
            auto pe = [](P x, const Edge& e) { return P{x.first, x.second + e.w}; };
            auto pv = [&](P x, int v) { return better(x, P{val[v], 0}); };
            const std::vector<P> ans = rerooting(G, P{-1, 0}, better, pe, pv);
            for (int v = 0; v < n; ++v) CHECK(ans[v] == brute_dp(G, v, -1, P{-1, 0}, better, pe, pv));
        }
    }
    std::cout << "[rerooting] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[rerooting] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> us(n - 1), vs(n - 1);
        std::vector<long long> ws(n - 1);
        for (int i = 1; i < n; ++i) us[i - 1] = i, vs[i - 1] = (int)(rng() % i), ws[i - 1] = (long long)(rng() % 1000);
        const Graph G = graph_from_edges(n, us, vs, ws, UNDIRECTED, 0);
        using P = std::pair<long long, long long>;
        auto t0 = std::chrono::steady_clock::now();
        const std::vector<P> ans = rerooting(G, P{0, 0}, [](P a, P b) { return P{a.first + b.first, a.second + b.second}; },
                                             [](P x, const Edge& e) { return P{x.first + x.second * e.w, x.second}; }, [](P x, int) { return P{x.first, x.second + 1}; });
        const long long t1 = ms_since(t0);
        Graph L(n);  // 一直線
        for (int v = 0; v + 1 < n; ++v) L.add_undirected(v, v + 1);
        t0 = std::chrono::steady_clock::now();
        const std::vector<long long> far = rerooting(L, 0LL, [](long long a, long long b) { return std::max(a, b); }, [](long long x, const Edge& e) { return x + e.w; },
                                                     [](long long x, int) { return x; });
        std::cout << "[rerooting] benchmark: N = " << n << "  距離の和 " << t1 << " ms, 一直線で最遠距離 " << ms_since(t0) << " ms  (check " << ans[0].first % 1000 + far[0]
                  << ")\n";
    }
}

}  // namespace rr_test


namespace rr_demo {
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
    banner("全方位木DP 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto ans = rerooting(G, e, merge, put_edge, put_vertex);", "vector<T>", "ans[v] = v を根にしたときの木DP の値"},
        {"merge(a, b) / e", "", "子からの値をまとめる / その単位元"},
        {"put_edge(x, e)", "", "子の側の値 x を辺 e（e.to が子）を通して持ち上げる"},
        {"put_vertex(x, v)", "", "まとめた値に頂点 v 自身の分を足す"},
    };
    cheat_sheet("API 一覧", rows, 58, 12);

    head("デモで使う木（重みつき）: 0-1 (2), 1-2 (3), 1-3 (1), 3-4 (4)");
    RUN(Graph G(5));
    RUN(G.add_undirected(0, 1, 2); G.add_undirected(1, 2, 3); G.add_undirected(1, 3, 1); G.add_undirected(3, 4, 4));

    api("【1】 各頂点からいちばん遠い頂点までの距離", "vector<long long>", "dp[v] = v から部分木の中でいちばん遠い距離。merge = max、辺で + w", "O(N)");
    code("auto far = rerooting(G, 0LL,");
    code("    [](long long a, long long b) { return std::max(a, b); },   // merge");
    code("    [](long long x, const Edge& e) { return x + e.w; },        // 辺を通ると + 重み");
    code("    [](long long x, int) { return x; });                       // 頂点の分はなし");
    auto far = rerooting(G, 0LL, [](long long a, long long b) { return std::max(a, b); }, [](long long x, const Edge& e) { return x + e.w; },
                         [](long long x, int) { return x; });
    SHOW(far, "05_tree_diameter の ecc と同じ");

    api("【2】 各頂点から全頂点への距離の和（ABC220 F）", "vector<pair<long long,long long>>", "値を {距離の和, 頂点数} の組にする。辺を通ると和が 頂点数 × 重み 増える", "O(N)");
    code("using P = std::pair<long long, long long>;   // {距離の和, 頂点数}");
    code("auto sum = rerooting(G, P{0, 0},");
    code("    [](P a, P b) { return P{a.first + b.first, a.second + b.second}; },");
    code("    [](P x, const Edge& e) { return P{x.first + x.second * e.w, x.second}; },");
    code("    [](P x, int) { return P{x.first, x.second + 1}; });   // 自分を 1 頂点として数える");
    using P = std::pair<long long, long long>;
    auto sum = rerooting(G, P{0, 0}, [](P a, P b) { return P{a.first + b.first, a.second + b.second}; },
                         [](P x, const Edge& e) { return P{x.first + x.second * e.w, x.second}; }, [](P x, int) { return P{x.first, x.second + 1}; });
    SHOW(sum, "{距離の和, 頂点数 = 5}");

    api("【3】 v を含む連結な部分集合（黒く塗る頂点）の数（EDPC V Subtree）", "vector<long long>", "dp[v] = Π (dp[子] + 1)。merge = 掛け算（単位元 1）", "O(N)");
    code("const long long M = 1000000007;");
    code("auto cnt = rerooting(G, 1LL, [&](long long a, long long b) { return a * b % M; },");
    code("    [&](long long x, const Edge&) { return (x + 1) % M; },   // 子の側を「全部白」にする 1 通りを足す");
    code("    [&](long long x, int) { return x; });");
    const long long M = 1000000007;
    auto cnt = rerooting(G, 1LL, [&](long long a, long long b) { return a * b % M; }, [&](long long x, const Edge&) { return (x + 1) % M; },
                         [&](long long x, int) { return x; });
    SHOW(cnt, "");
    note("merge は割り算のいらない形でよい（左右からの累積で「自分以外」をまとめるので、逆元がなくても動く）。");

    head("書き方のコツ");
    note("まず「根を固定した木DP」を、merge / put_edge / put_vertex の 3 つに分けて書く。");
    note("値が複数あるとき（和と個数など）は pair や struct にまとめる。");
    note("根を 1 つに決めてよいなら 05_rooted_tree の木DP の型（order を逆順に回す）で十分。");
}

}  // namespace rr_demo

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
    rr_test::random_test();
    rr_test::benchmark();
    rr_demo::demo();
    return 0;
}
