// =============================================================================
//  二部グラフ判定・2 色の塗り分け  bipartite / xor_coloring       --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 01 探索
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed。辺の向きは無視する（有向グラフを渡しても無向として扱う）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  bipartite(G)   : 隣り合う頂点を違う色（0 / 1）に塗れるか
//  xor_coloring(G): 辺 u-v の重み w について color[u] xor color[v] = w % 2 となるように塗れるか
//                   （「u と v は同じ / 違う」の条件が並ぶ問題。重み 0 = 同じ、1 = 違う）
//  どちらも、矛盾があれば矛盾する閉路（二部グラフなら奇数長の閉路）を 1 つ返す。
//
//  機能一覧                                                                       計算量
//    auto b = bipartite(G) / xor_coloring(G)                                       O(N + M)
//    b.ok                       全体が矛盾なく塗れるか
//    b.color[v]                 色 0 / 1（成分ごとに番号最小の頂点が 0）
//    b.comp[v] / b.comp_count() 成分番号 / 成分の数
//    b.cnt[c]                   成分 c の {色 0 の頂点数, 色 1 の頂点数}（ABC282 D など）
//    b.comp_ok[c]               成分 c が矛盾なく塗れるか
//    b.odd_cycle                矛盾する閉路の頂点列（ok のときは空）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    bip_test::random_test() / bip_test::benchmark() / bip_demo::demo()
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
#ifndef GG_BIPARTITE
#define GG_BIPARTITE
// ---- 二部グラフ判定・2 色の塗り分け（グラフ・グリッド第二世代）  辺の向きは無視する ----
struct Bipartite {
    bool ok = true;                        // 全体が二部グラフ（矛盾なく塗れる）か
    std::vector<int> color;                // 0 / 1（成分ごとに番号最小の頂点が 0）。矛盾があっても BFS で塗った色が入る
    std::vector<int> comp;                 // 成分番号（番号の小さい頂点を含む成分から 0, 1, ...）
    std::vector<std::array<int, 2>> cnt;   // cnt[c] = {成分 c の色 0 の頂点数, 色 1 の頂点数}
    std::vector<int> comp_ok;              // 成分ごとに矛盾なく塗れるか（1 / 0）
    std::vector<int> odd_cycle;            // ok でないとき、矛盾する閉路の頂点列を 1 つ（bipartite なら奇数長の閉路）
    int comp_count() const { return (int)cnt.size(); }
};
namespace gg_internal {
// 辺 u-v に「色の xor = p」の条件を付けて BFS で塗る（p は辺ごとに parity(e) で決める）
template <class Par>
Bipartite two_coloring(const Graph& G, Par parity) {
    const int n = G.size();
    // 向きを無視した隣接（CSR）: 辺 id ごとに両方向
    std::vector<int> st(n + 1, 0);
    for (const RawEdge& e : G.es) ++st[e.from + 1], ++st[e.to + 1];
    for (int v = 0; v < n; ++v) st[v + 1] += st[v];
    std::vector<int> to(st[n]), eid(st[n]), pos(st.begin(), st.end() - 1);
    for (int i = 0; i < G.edge_count(); ++i) {
        const RawEdge& e = G.es[i];
        to[pos[e.from]] = e.to, eid[pos[e.from]++] = i;
        to[pos[e.to]] = e.from, eid[pos[e.to]++] = i;
    }
    Bipartite b;
    b.color.assign(n, -1), b.comp.assign(n, -1);
    std::vector<int> par(n, -1), dep(n, 0), q;
    q.reserve(n);
    int bad_u = -1, bad_v = -1;
    for (int s = 0; s < n; ++s) {
        if (b.comp[s] != -1) continue;
        const int c = (int)b.cnt.size();
        b.cnt.push_back({0, 0}), b.comp_ok.push_back(1);
        b.comp[s] = c, b.color[s] = 0, q.clear(), q.push_back(s);
        for (size_t h = 0; h < q.size(); ++h) {
            const int v = q[h];
            ++b.cnt[c][b.color[v]];
            for (int k = st[v]; k < st[v + 1]; ++k) {
                const int u = to[k], p = parity(eid[k]) ? 1 : 0;
                if (b.comp[u] == -1) {
                    b.comp[u] = c, b.color[u] = b.color[v] ^ p, par[u] = v, dep[u] = dep[v] + 1, q.push_back(u);
                } else if ((b.color[u] ^ b.color[v]) != p) {
                    if (b.ok) bad_u = v, bad_v = u;
                    b.ok = false, b.comp_ok[c] = 0;
                }
            }
        }
    }
    if (!b.ok) {  // BFS 木の上の u → lca → v と、矛盾した辺 v-u で閉路になる
        int u = bad_u, v = bad_v;
        std::vector<int> a, d;
        while (dep[u] > dep[v]) a.push_back(u), u = par[u];
        while (dep[v] > dep[u]) d.push_back(v), v = par[v];
        while (u != v) a.push_back(u), d.push_back(v), u = par[u], v = par[v];
        a.push_back(u);
        b.odd_cycle = a;
        for (int i = (int)d.size() - 1; i >= 0; --i) b.odd_cycle.push_back(d[i]);
    }
    return b;
}
}  // namespace gg_internal
// 二部グラフ判定（隣り合う頂点を違う色に塗れるか）。ABC282 D など
inline Bipartite bipartite(const Graph& G) {
    return gg_internal::two_coloring(G, [](int) { return 1; });
}
// xor の条件での塗り分け: 辺 u-v の重み w について color[u] xor color[v] = w % 2（w が偶数なら同じ色、奇数なら違う色）
inline Bipartite xor_coloring(const Graph& G) {
    return gg_internal::two_coloring(G, [&](int id) { return (int)(G.es[id].w & 1); });
}
#endif  // GG_BIPARTITE
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "bipartite";
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

namespace bip_test {

// 愚直解: 成分ごとに全部の塗り方を試す（n <= 12）
static void check(const Graph& G, const Bipartite& b, bool use_w) {
    const int n = G.size();
    auto par = [&](int id) { return use_w ? (int)(G.es[id].w & 1) : 1; };
    // 成分（向きを無視）
    std::vector<int> comp(n, -1);
    int k = 0;
    for (int s = 0; s < n; ++s) {
        if (comp[s] != -1) continue;
        std::vector<int> q{s};
        comp[s] = k;
        for (size_t h = 0; h < q.size(); ++h)
            for (const RawEdge& e : G.es) {
                int o = -1;
                if (e.from == q[h]) o = e.to;
                else if (e.to == q[h]) o = e.from;
                if (o >= 0 && comp[o] == -1) comp[o] = k, q.push_back(o);
            }
        ++k;
    }
    CHECK(b.comp == comp && b.comp_count() == k && (int)b.comp_ok.size() == k);
    std::vector<char> ok(k, 0);
    for (int c = 0; c < k; ++c) {
        std::vector<int> vs;
        for (int v = 0; v < n; ++v)
            if (comp[v] == c) vs.push_back(v);
        for (int mask = 0; mask < (1 << vs.size()); ++mask) {
            std::vector<int> col(n, 0);
            for (size_t i = 0; i < vs.size(); ++i) col[vs[i]] = (mask >> i) & 1;
            bool good = true;
            for (int id = 0; id < G.edge_count() && good; ++id) {
                const RawEdge& e = G.es[id];
                if (comp[e.from] == c && (col[e.from] ^ col[e.to]) != par(id)) good = false;
            }
            if (good) {
                ok[c] = 1;
                break;
            }
        }
        CHECK((bool)b.comp_ok[c] == (bool)ok[c]);
    }
    bool all = true;
    for (int c = 0; c < k; ++c) all &= (bool)ok[c];
    CHECK(b.ok == all);
    // 色: ok の成分では条件を満たし、番号最小の頂点が 0。数も一致
    std::vector<std::array<int, 2>> cnt(k, {0, 0});
    std::vector<char> first(k, 1);
    for (int v = 0; v < n; ++v) {
        CHECK(b.color[v] == 0 || b.color[v] == 1);
        ++cnt[comp[v]][b.color[v]];
        if (first[comp[v]]) {
            CHECK(b.color[v] == 0);
            first[comp[v]] = 0;
        }
    }
    CHECK(b.cnt == cnt);
    for (int id = 0; id < G.edge_count(); ++id) {
        const RawEdge& e = G.es[id];
        if (ok[comp[e.from]]) CHECK((b.color[e.from] ^ b.color[e.to]) == par(id));
    }
    // 矛盾する閉路: 連続する頂点が辺で結ばれ、条件の xor の和が 1（同じ頂点は出てこない）
    if (b.ok) {
        CHECK(b.odd_cycle.empty());
    } else {
        const std::vector<int>& cy = b.odd_cycle;
        CHECK(!cy.empty());
        std::vector<int> sorted = cy;
        std::sort(sorted.begin(), sorted.end());
        CHECK(std::unique(sorted.begin(), sorted.end()) == sorted.end());
        // 閉路の各辺について、使える辺のうちパリティの合計が奇数になる選び方があること
        // （多重辺があるので、各区間の「使える辺のパリティの集合」で DP）
        int reach = 1;  // bit p = 合計のパリティ p が作れる
        const int L = (int)cy.size();
        for (int i = 0; i < L; ++i) {
            const int a = cy[i], c2 = cy[(i + 1) % L];
            int ps = 0;
            for (int id = 0; id < G.edge_count(); ++id) {
                const RawEdge& e = G.es[id];
                if ((e.from == a && e.to == c2) || (e.from == c2 && e.to == a)) ps |= 1 << par(id);
            }
            CHECK(ps != 0);
            int nr = 0;
            for (int p = 0; p < 2; ++p)
                if (reach >> p & 1)
                    for (int q = 0; q < 2; ++q)
                        if (ps >> q & 1) nr |= 1 << (p ^ q);
            reach = nr;
        }
        CHECK(reach >> 1 & 1);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 2500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 10), m = rnd(rng, 0, 2 * n);
        Graph G(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            const Cost w = rnd(rng, -3, 3);
            (rng() & 1) ? G.add_undirected(a, b, w) : G.add_edge(a, b, w);
        }
        check(G, bipartite(G), false);
        check(G, xor_coloring(G), true);
    }
    std::cout << "[bipartite] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[bipartite] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> us(n), vs(n);
        for (int i = 0; i < n; ++i) us[i] = (int)(rng() % (n / 2)) * 2, vs[i] = (int)(rng() % (n / 2)) * 2 + 1;  // 偶数-奇数の辺だけ
        const Graph G = graph_from_edges(n, us, vs, UNDIRECTED, 0);
        auto t0 = std::chrono::steady_clock::now();
        const Bipartite b = bipartite(G);
        std::cout << "[bipartite] benchmark: N = M = " << n << "  " << ms_since(t0) << " ms  (ok = " << b.ok << ")\n";
    }
}

}  // namespace bip_test


namespace bip_demo {
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
    banner("二部グラフ判定 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto b = bipartite(G);", "Bipartite", "隣り合う頂点を違う色に塗れるか"},
        {"auto b = xor_coloring(G);", "Bipartite", "color[u] xor color[v] = w % 2（同じ / 違う の条件）"},
        {"b.ok / b.comp_ok[c]", "bool", "全体 / 成分 c が矛盾なく塗れるか"},
        {"b.color[v]", "int", "色 0 / 1（成分ごとに番号最小の頂点が 0）"},
        {"b.comp[v] / b.comp_count()", "int", "成分番号 / 成分の数"},
        {"b.cnt[c]", "array<int, 2>", "成分 c の色 0 / 色 1 の頂点数"},
        {"b.odd_cycle", "vector<int>", "矛盾する閉路（奇数長の閉路）"},
    };
    cheat_sheet("API 一覧", rows, 34, 18);

    head("デモで使うグラフ（無向）: 0-1, 1-2, 3-4（成分 {0,1,2}, {3,4}, {5}）");
    RUN(Graph G(6));
    RUN(G.add_undirected(0, 1); G.add_undirected(1, 2); G.add_undirected(3, 4));

    api("【1】 bipartite(G)", "Bipartite", "隣り合う頂点を違う色（0 / 1）に塗れるか。色・成分・成分ごとの色別の数も返す", "O(N + M)");
    RUN(auto b = bipartite(G));
    SHOW(b.ok, "二部グラフ");
    SHOW(b.color, "色（成分ごとに番号最小の頂点が 0）");
    SHOW(b.comp, "成分番号");
    SHOW(b.cnt, "成分ごとの {色 0 の数, 色 1 の数}");
    code("// ABC282 D: 辺を 1 本足しても二部グラフのままの組の数 = N(N-1)/2 - M - Σ(同じ成分・同じ色の組)");
    code("long long ans = (long long)N * (N - 1) / 2 - M;");
    code("for (auto [a, c] : b.cnt) ans -= (long long)a * (a - 1) / 2 + (long long)c * (c - 1) / 2;");

    api("【2】 矛盾するとき: b.ok = false、b.odd_cycle", "", "奇数長の閉路があると二部グラフではない。その閉路を 1 つ返す", "O(N + M)");
    RUN(G.add_undirected(2, 0));
    RUN(auto b2 = bipartite(G));
    SHOW(b2.ok, "0-1-2-0 は長さ 3 の閉路");
    SHOW(b2.comp_ok, "成分ごと");
    SHOW(b2.odd_cycle, "奇数長の閉路");

    api("【3】 xor_coloring(G)", "Bipartite",
        "辺の重み w が偶数なら「同じ色」、奇数なら「違う色」の条件。条件がすべて満たせるか（重み付き Union-Find の xor 版の代わり）", "O(N + M)");
    code("// 例: 「A と B は同じ」「B と C は違う」「A と C は同じ」→ 矛盾");
    RUN(Graph X(3));
    RUN(X.add_undirected(0, 1, 0); X.add_undirected(1, 2, 1); X.add_undirected(0, 2, 0));
    RUN(auto bx = xor_coloring(X));
    SHOW(bx.ok, "矛盾");
    SHOW(bx.odd_cycle, "矛盾する条件の閉路");
    note("正直者・嘘つき系（ABC 「証言」）、スイッチの xor など。");
}

}  // namespace bip_demo

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
    bip_test::random_test();
    bip_test::benchmark();
    bip_demo::demo();
    return 0;
}
