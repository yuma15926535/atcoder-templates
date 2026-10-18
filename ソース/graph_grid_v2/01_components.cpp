// =============================================================================
//  連結成分  components / weak_components / components_core     --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 01 探索
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed。成分番号は「番号の小さい頂点を含む成分」から 0, 1, 2, ...
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  連結成分を 1 回だけ求めるならこれ（BFS で O(N + M)）。辺を足しながら何度も聞くなら 03_union_find。
//  グリッドの塗りつぶしは 09_grid_components（中身はこの components_core）。
//
//  高速化のポイント
//    - 成分ごとの頂点の一覧は 1 本の配列（CSR）。group(c) はコピーせずに参照（IntSpan）を返す
//
//  機能一覧                                                                       計算量
//    components(G)                    無向グラフの連結成分                           O(N + M)
//    weak_components(G)               有向グラフの弱連結成分（向きを無視）           O(N + M)
//    components_core(n, adj[, alive]) 隣接を自分で書く版。alive(v) = false の頂点は対象外（id = -1）
//    component_edge_counts(G, c)      成分ごとの辺の本数（ABC292 D）                 O(N + M)
//    結果 Components c:
//      c.count / c.id[v]              成分の数 / 成分番号
//      c.group(k) / c.size_of(k)      成分 k の頂点（昇順、IntSpan）/ 頂点数
//      c.comp_size(v) / c.same(u, v)  v を含む成分の頂点数 / 同じ成分か
//      c.sizes() / c.max_size() / c.groups()
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    cc_test::random_test() / cc_test::benchmark() / cc_demo::demo()
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
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "components";
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

namespace cc_test {

struct DSU {  // 愚直解用
    std::vector<int> p;
    explicit DSU(int n) : p(n) { std::iota(p.begin(), p.end(), 0); }
    int f(int x) { return p[x] == x ? x : p[x] = f(p[x]); }
    void u(int a, int b) { p[f(a)] = f(b); }
};
static void check(const Components& c, const std::vector<int>& label, const std::vector<char>& alive) {
    const int n = (int)label.size();
    // 番号の付け方: 番号の小さい頂点を含む成分から 0, 1, ...
    std::vector<int> ex(n, -1);
    int k = 0;
    std::vector<int> map(n, -1);
    for (int v = 0; v < n; ++v) {
        if (!alive[v]) continue;
        if (map[label[v]] == -1) map[label[v]] = k++;
        ex[v] = map[label[v]];
    }
    CHECK(c.id == ex && c.count == k);
    std::vector<std::vector<int>> g(k);
    for (int v = 0; v < n; ++v)
        if (ex[v] >= 0) g[ex[v]].push_back(v);
    CHECK(c.groups() == g);
    int mx = 0;
    for (int i = 0; i < k; ++i) {
        CHECK(c.group(i).to_vector() == g[i] && c.size_of(i) == (int)g[i].size());
        mx = std::max(mx, (int)g[i].size());
    }
    CHECK(c.max_size() == mx && (int)c.sizes().size() == k);
    for (int u = 0; u < n; ++u) {
        CHECK(c.comp_size(u) == (ex[u] < 0 ? 0 : (int)g[ex[u]].size()));
        for (int v = 0; v < n; ++v) CHECK(c.same(u, v) == (ex[u] >= 0 && ex[u] == ex[v]));
    }
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 20), m = rnd(rng, 0, 2 * n);
        const bool und = rng() & 1;
        Graph G(n);
        DSU d(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            und ? G.add_undirected(a, b) : G.add_edge(a, b);
            d.u(a, b);
        }
        std::vector<int> label(n);
        for (int v = 0; v < n; ++v) label[v] = d.f(v);
        const std::vector<char> all(n, 1);
        if (und) check(components(G), label, all);
        check(weak_components(G), label, all);
        // 辺の本数
        const Components c = weak_components(G);
        std::vector<int> ec(c.count, 0);
        for (const RawEdge& e : G.es) ++ec[c.id[e.from]];
        CHECK(component_edge_counts(G, c) == ec);
        // alive 付き
        std::vector<char> alive(n);
        for (auto& a : alive) a = (rng() % 4 != 0);
        DSU d2(n);
        for (const RawEdge& e : G.es)
            if (alive[e.from] && alive[e.to]) d2.u(e.from, e.to);
        std::vector<int> label2(n);
        for (int v = 0; v < n; ++v) label2[v] = d2.f(v);
        const Graph R = G.reversed();
        const Components ca = components_core(n, [&](int v, auto&& push) {
            for (const Edge& e : G[v]) push(e.to);
            for (const Edge& e : R[v]) push(e.to);
        }, [&](int v) { return (bool)alive[v]; });
        check(ca, label2, alive);
    }
    std::cout << "[components] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[components] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> us(n / 2), vs(n / 2);
        for (int i = 0; i < n / 2; ++i) us[i] = (int)(rng() % n), vs[i] = (int)(rng() % n);
        const Graph G = graph_from_edges(n, us, vs, UNDIRECTED, 0);
        auto t0 = std::chrono::steady_clock::now();
        const Components c = components(G);
        std::cout << "[components] benchmark: N = " << n << ", M = " << n / 2 << "  " << ms_since(t0) << " ms（成分 " << c.count << " 個）\n";
    }
}

}  // namespace cc_test


namespace cc_demo {
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
    banner("連結成分 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto c = components(G);", "Components", "無向グラフの連結成分"},
        {"auto c = weak_components(G);", "Components", "有向グラフで、向きを無視した成分"},
        {"components_core(n, adj[, alive])", "Components", "隣接を自分で書く。alive(v) が false は対象外"},
        {"c.count / c.id[v]", "int", "成分の数 / 成分番号（番号の小さい頂点の成分から 0, 1, ...）"},
        {"c.group(k) / c.size_of(k)", "IntSpan / int", "成分 k の頂点（昇順）/ 頂点数"},
        {"c.comp_size(v) / c.same(u, v)", "int / bool", "v の成分の頂点数 / 同じ成分か"},
        {"c.sizes() / c.max_size() / c.groups()", "", "全成分の頂点数 / 最大 / vector<vector<int>>"},
        {"component_edge_counts(G, c)", "vector<int>", "成分ごとの辺の本数"},
    };
    cheat_sheet("API 一覧", rows, 42, 16);

    head("デモで使うグラフ（無向）: 0-3, 3-5, 1-4, 2 は単独");
    RUN(Graph G(6));
    RUN(G.add_undirected(0, 3); G.add_undirected(3, 5); G.add_undirected(1, 4));

    api("【1】 components(G)", "Components", "連結成分。成分番号は番号の小さい頂点を含む成分から 0, 1, 2, ...", "O(N + M)");
    RUN(auto c = components(G));
    SHOW(c.count, "成分の数");
    SHOW(c.id, "各頂点の成分番号");
    SHOW(c.group(0), "成分 0 の頂点（昇順）");
    SHOW(c.size_of(1), "成分 1 の頂点数");
    SHOW(c.comp_size(5), "頂点 5 を含む成分の頂点数");
    SHOW(c.same(0, 5), "同じ成分か");
    SHOW(c.sizes(), "全成分の頂点数");
    SHOW(c.max_size(), "いちばん大きい成分");
    SHOW(c.groups(), "成分ごとの頂点の一覧");
    note("for (int v : c.group(k)) { ... } で成分ごとに処理できる。");
    note("「あと何本辺を足せば連結になるか」= c.count - 1。");

    api("【2】 component_edge_counts(G, c)", "vector<int>", "成分ごとの辺の本数（多重辺・自己ループも 1 本）。辺の本数 = 頂点数 なら閉路がちょうど 1 つ（ABC292 D）", "O(N + M)");
    SHOW(component_edge_counts(G, c), "成分 0 は 2 本、成分 1 は 1 本、成分 2 は 0 本");

    api("【3】 weak_components(G)", "Components", "有向グラフで、辺の向きを無視してつながっている成分（弱連結成分）", "O(N + M)");
    RUN(Graph D(4));
    RUN(D.add_edge(1, 0); D.add_edge(2, 3));
    SHOW(weak_components(D).id, "1→0 と 2→3");
    SHOW(components(D).id, "components は向きどおりにたどるので、有向グラフには使わない");

    api("【4】 components_core(n, adj, alive)", "Components", "隣接を自分で書く。alive(v) が false の頂点は対象外（id = -1、たどらない）", "O(n + 隣接の数)");
    code("// 例: 頂点 3 を取り除いたときの成分（頂点を 1 つずつ消して調べる問題）");
    RUN(auto c3 = components_core(G.size(), [&](int v, auto&& push) { for (auto& e : G[v]) push(e.to); }, [&](int v) { return v != 3; }));
    SHOW(c3.id, "3 は -1、0 と 5 は別の成分になる");
    SHOW(c3.count, "成分の数");
    note("グリッドの塗りつぶしは 09_grid_components の grid_components(g, '#')。");
}

}  // namespace cc_demo

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
    cc_test::random_test();
    cc_test::benchmark();
    cc_demo::demo();
    return 0;
}
