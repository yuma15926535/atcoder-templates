// =============================================================================
//  ワーシャルフロイド（全点対の最短距離・経路・負閉路・辺の追加による更新）  WarshallFloyd
//                                                                 --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 最短路・状態グラフ
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed。重みは負でもよい（負閉路があれば、そこを通れる組は −∞）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  高速化・安全のポイント
//    - 距離は 1 本の配列（N × N）に持つ（vector の vector より速い）。W[i][j] で読める
//    - DIST_INF どうしを足さない（あふれない）。負閉路で値が小さくなりすぎないよう -DIST_INF で止める
//
//  機能一覧                                                                       計算量
//    WarshallFloyd W(G)              G の辺を全部入れて計算まで行う                       O(N^3)
//    WarshallFloyd W(N); W.add_edge(u, v, w) / W.add_undirected(u, v, w); W.run();  （行列で与えられたときなど）
//    W[i][j] / W.dist(i, j)          最短距離（DIST_INF = 届かない、-DIST_INF = −∞）          O(1)
//    W.reachable(i, j) / W.neg_inf(i, j) / W.negative_cycle   届くか / −∞ か / 負閉路があるか
//    W.path(i, j)                    i → j の頂点列                                        O(経路の長さ)
//    W.update_edge(u, v, w) / W.update_undirected(u, v, w)   計算の後に辺を足して全点対を更新  O(N^2)
//    W.matrix()                      vector<vector<Cost>> にしたもの                          O(N^2)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    wf_test::random_test() / wf_test::benchmark() / wf_demo::demo()
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
#ifndef GG_WARSHALL_FLOYD
#define GG_WARSHALL_FLOYD
// ---- ワーシャルフロイド（グラフ・グリッド第二世代）  全点対の最短距離・経路・辺の追加による更新 ----
// W[i][j] = i → j の最短距離（届かなければ DIST_INF、負閉路を通っていくらでも小さくできるなら -DIST_INF = −∞）
struct WarshallFloyd {
    int n = 0;
    std::vector<Cost> d;    // d[i * n + j]（W[i][j] で読める）
    std::vector<int> nxt;   // nxt[i * n + j] = i → j の最短路で i の次の頂点（経路の復元用。なければ -1）
    bool negative_cycle = false;  // 負閉路があるか（run の後）

    WarshallFloyd() {}
    explicit WarshallFloyd(int n_) : n(n_), d((size_t)n_ * n_, DIST_INF), nxt((size_t)n_ * n_, -1) {
        for (int i = 0; i < n; ++i) d[(size_t)i * n + i] = 0, nxt[(size_t)i * n + i] = i;
    }
    explicit WarshallFloyd(const Graph& G) : WarshallFloyd(G.size()) {  // G の辺を全部入れて run() まで行う
        for (int v = 0; v < n; ++v)
            for (const Edge& e : G[v]) add_edge(v, e.to, e.w);
        run();
    }
    // 辺 u → v（重み w）を入れる（多重辺は小さいほう）。run() の前に使う
    void add_edge(int u, int v, Cost w) {
        Cost& x = d[(size_t)u * n + v];
        if (w < x) x = w, nxt[(size_t)u * n + v] = v;
    }
    void add_undirected(int u, int v, Cost w) { add_edge(u, v, w), add_edge(v, u, w); }
    // 全点対の最短距離を求める。O(N^3)
    void run() {
        const Cost LO = -DIST_INF;
        for (int k = 0; k < n; ++k) {
            const Cost* dk = &d[(size_t)k * n];
            for (int i = 0; i < n; ++i) {
                Cost* di = &d[(size_t)i * n];
                const Cost dik = di[k];
                if (dik == DIST_INF) continue;
                int* ni = &nxt[(size_t)i * n];
                const int nik = ni[k];
                for (int j = 0; j < n; ++j) {
                    if (dk[j] == DIST_INF) continue;
                    Cost c = dik + dk[j];
                    if (c < LO) c = LO;  // 負閉路で値が小さくなりすぎないように止める
                    if (c < di[j]) di[j] = c, ni[j] = nik;
                }
            }
        }
        // 負閉路の上の頂点 k（d[k][k] < 0）を通れる組は −∞
        negative_cycle = false;
        std::vector<int> neg;
        for (int k = 0; k < n; ++k)
            if (d[(size_t)k * n + k] < 0) neg.push_back(k);
        if (neg.empty()) return;
        negative_cycle = true;
        for (int k : neg)
            for (int i = 0; i < n; ++i) {
                if (d[(size_t)i * n + k] == DIST_INF) continue;
                for (int j = 0; j < n; ++j)
                    if (d[(size_t)k * n + j] != DIST_INF) d[(size_t)i * n + j] = -DIST_INF, nxt[(size_t)i * n + j] = -1;
            }
    }
    const Cost* operator[](int i) const { return &d[(size_t)i * n]; }  // W[i][j] で距離
    Cost dist(int i, int j) const { return d[(size_t)i * n + j]; }
    bool reachable(int i, int j) const { return d[(size_t)i * n + j] != DIST_INF; }
    bool neg_inf(int i, int j) const { return d[(size_t)i * n + j] == -DIST_INF; }
    // i → j の最短路の頂点列（届かない・−∞ なら空）
    std::vector<int> path(int i, int j) const {
        std::vector<int> p;
        if (!reachable(i, j) || neg_inf(i, j)) return p;
        for (int x = i; x != j; x = nxt[(size_t)x * n + j]) p.push_back(x);
        p.push_back(j);
        return p;
    }
    // run() の後に辺 u → v（重み w）を足して、全点対を更新する。O(N^2)。負閉路ができないこと
    // （辺を消していく問題を、逆順に「足していく」問題にして使う。ABC375 F など）
    void update_edge(int u, int v, Cost w) {
        const Cost* dv = &d[(size_t)v * n];
        for (int i = 0; i < n; ++i) {
            const Cost diu = d[(size_t)i * n + u];
            if (diu == DIST_INF) continue;
            Cost* di = &d[(size_t)i * n];
            int* ni = &nxt[(size_t)i * n];
            const int first = (i == u) ? v : ni[u];
            for (int j = 0; j < n; ++j) {
                if (dv[j] == DIST_INF) continue;
                const Cost c = diu + w + dv[j];
                if (c < di[j]) di[j] = c, ni[j] = first;
            }
        }
    }
    void update_undirected(int u, int v, Cost w) { update_edge(u, v, w), update_edge(v, u, w); }
    std::vector<std::vector<Cost>> matrix() const {  // vector<vector<Cost>> にしたもの
        std::vector<std::vector<Cost>> m(n, std::vector<Cost>(n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) m[i][j] = d[(size_t)i * n + j];
        return m;
    }
};
#endif  // GG_WARSHALL_FLOYD
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "warshall_floyd";
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

namespace wf_test {

// 愚直: 始点ごとに「辺の本数 n-1 以下の歩道の最小」と、負の閉じた歩道から行ける頂点（−∞）
static std::vector<std::vector<Cost>> brute(const Graph& G) {
    const int n = G.size();
    std::vector<int> neg(n, 0);
    for (int c = 0; c < n; ++c) {
        std::vector<Cost> d(n, DIST_INF);
        d[c] = 0;
        for (int k = 0; k < n; ++k) {
            std::vector<Cost> nd(n, DIST_INF);
            for (int v = 0; v < n; ++v)
                if (d[v] != DIST_INF)
                    for (const Edge& e : G[v]) nd[e.to] = std::min(nd[e.to], d[v] + e.w);
            d = nd;
            if (d[c] < 0) neg[c] = 1;
        }
    }
    std::vector<std::vector<int>> reach(n, std::vector<int>(n, 0));
    for (int s = 0; s < n; ++s) {
        std::vector<int> st{s};
        reach[s][s] = 1;
        while (!st.empty()) {
            const int v = st.back();
            st.pop_back();
            for (const Edge& e : G[v])
                if (!reach[s][e.to]) reach[s][e.to] = 1, st.push_back(e.to);
        }
    }
    std::vector<std::vector<Cost>> res(n, std::vector<Cost>(n, DIST_INF));
    for (int s = 0; s < n; ++s) {
        std::vector<Cost> d(n, DIST_INF);
        d[s] = 0;
        for (int k = 0; k + 1 < n; ++k) {
            std::vector<Cost> nd = d;
            for (int v = 0; v < n; ++v)
                if (d[v] != DIST_INF)
                    for (const Edge& e : G[v]) nd[e.to] = std::min(nd[e.to], d[v] + e.w);
            d = nd;
        }
        for (int t = 0; t < n; ++t) {
            res[s][t] = d[t];
            for (int c = 0; c < n; ++c)
                if (neg[c] && reach[s][c] && reach[c][t]) res[s][t] = -DIST_INF;
        }
    }
    return res;
}
static Cost min_edge(const Graph& G, int a, int b) {
    Cost w = DIST_INF;
    for (const Edge& e : G[a])
        if (e.to == b) w = std::min(w, e.w);
    return w;
}
static void check_paths(const WarshallFloyd& W, const Graph& G) {
    const int n = G.size();
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const std::vector<int> p = W.path(i, j);
            if (!W.reachable(i, j) || W.neg_inf(i, j)) {
                CHECK(p.empty());
                continue;
            }
            CHECK(!p.empty() && p.front() == i && p.back() == j && (int)p.size() <= n);
            Cost sum = 0;
            for (size_t k = 0; k + 1 < p.size(); ++k) {
                const Cost w = min_edge(G, p[k], p[k + 1]);
                CHECK(w != DIST_INF);
                sum += w;
            }
            CHECK(sum == W.dist(i, j));
        }
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 7), m = rnd(rng, 0, 14);
        const bool und = (it % 4 == 0);
        const int lo = (it % 3 == 0) ? -3 : (it % 3 == 1) ? -1 : 0;
        Graph G(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            const Cost w = rnd(rng, lo, 9);
            und ? G.add_undirected(a, b, w) : G.add_edge(a, b, w);
        }
        const auto ex = brute(G);
        const WarshallFloyd W(G);
        bool any = false;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                CHECK(W[i][j] == ex[i][j] && W.dist(i, j) == ex[i][j]);
                CHECK(W.reachable(i, j) == (ex[i][j] != DIST_INF) && W.neg_inf(i, j) == (ex[i][j] == -DIST_INF));
                any |= (ex[i][j] == -DIST_INF);
            }
        CHECK(W.negative_cycle == any && W.matrix() == ex);
        check_paths(W, G);
        // 手で辺を入れても同じ
        WarshallFloyd W2(n);
        for (const RawEdge& e : G.es) und ? W2.add_undirected(e.from, e.to, e.w) : W2.add_edge(e.from, e.to, e.w);
        W2.run();
        CHECK(W2.matrix() == ex);
        // update_edge: 重み 0 以上のグラフで、辺を 1 本ずつ足して作り直したものと比べる
        {
            Graph H(n);
            WarshallFloyd U(n);
            U.run();
            for (int k = 0; k < 6; ++k) {
                const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
                const Cost w = rnd(rng, 0, 9);
                if (und) H.add_undirected(a, b, w), U.update_undirected(a, b, w);
                else H.add_edge(a, b, w), U.update_edge(a, b, w);
                const WarshallFloyd F(H);
                CHECK(U.matrix() == F.matrix());
                check_paths(U, H);
            }
        }
    }
    std::cout << "[warshall_floyd] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[warshall_floyd] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {400, 800}) {
        if (TEST_SCALE == 1 && n > 400) break;
        Graph G(n);
        for (int i = 0; i < 4 * n; ++i) G.add_edge((int)(rng() % n), (int)(rng() % n), (Cost)(rng() % 1000000000));
        auto t0 = std::chrono::steady_clock::now();
        const WarshallFloyd W(G);
        const long long t1 = ms_since(t0);
        WarshallFloyd U = W;
        t0 = std::chrono::steady_clock::now();
        for (int k = 0; k < n; ++k) U.update_edge((int)(rng() % n), (int)(rng() % n), (Cost)(rng() % 1000));
        std::cout << "[warshall_floyd] benchmark: N = " << n << ", M = " << 4 * n << "  run " << t1 << " ms, update_edge " << n << " 回 " << ms_since(t0)
                  << " ms  (check " << W[0][n - 1] % 1000 + U[n - 1][0] % 1000 << ")\n";
    }
}

}  // namespace wf_test


namespace wf_demo {
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
struct Row2 {  // 行列の 1 行を「−∞」「DIST_INF」で表示するための入れ物（デモ用）
    std::vector<Cost> d;
};
[[maybe_unused]] static std::string to_s(const Row2& x) {
    std::string s = "{";
    for (size_t i = 0; i < x.d.size(); ++i)
        s += (i ? ", " : "") + (x.d[i] == -DIST_INF ? std::string("-DIST_INF") : x.d[i] == DIST_INF ? std::string("DIST_INF") : std::to_string(x.d[i]));
    return s + "}";
}

static void demo() {
    banner("ワーシャルフロイド 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"WarshallFloyd W(G);", "", "全点対の最短距離を計算する（O(N^3)）"},
        {"WarshallFloyd W(N); W.add_edge(u, v, w); W.run();", "", "辺を自分で入れる（行列で与えられたとき）"},
        {"W[i][j] / W.dist(i, j)", "Cost", "最短距離（DIST_INF 届かない / -DIST_INF は −∞）"},
        {"W.path(i, j)", "vector<int>", "経路の頂点列"},
        {"W.negative_cycle / W.neg_inf(i, j)", "bool", "負閉路があるか / i → j が −∞ か"},
        {"W.update_edge(u, v, w)", "", "計算の後に辺を足して更新（O(N^2)）"},
    };
    cheat_sheet("API 一覧", rows, 52, 14);

    head("デモで使うグラフ（無向）: 0-1 (4), 1-2 (1), 0-2 (7), 2-3 (2)");
    RUN(Graph G(4));
    RUN(G.add_undirected(0, 1, 4); G.add_undirected(1, 2, 1); G.add_undirected(0, 2, 7); G.add_undirected(2, 3, 2));

    api("【1】 WarshallFloyd W(G)", "WarshallFloyd", "全部の組 (i, j) の最短距離。W[i][j] で読む", "O(N^3)");
    RUN(WarshallFloyd W(G));
    SHOW(W[0][3], "0 → 1 → 2 → 3 で 4 + 1 + 2");
    SHOW(W.matrix(), "全部の組");
    SHOW(W.path(0, 3), "経路");
    note("N ≤ 400 程度なら全点対の距離はこれ。N が大きくて始点が少しならダイクストラを始点の数だけ。");

    api("【2】 WarshallFloyd W(N) と add_edge / run", "WarshallFloyd", "隣接行列や「i から j へコスト c」の形で与えられたとき。多重辺は小さいほう", "O(N^3)");
    RUN(WarshallFloyd M(3));
    RUN(M.add_edge(0, 1, 5); M.add_edge(1, 2, 3); M.add_edge(0, 2, 10); M.add_edge(0, 2, 9));
    RUN(M.run());
    SHOW(M[0][2], "0 → 1 → 2 で 8（直接は 9）");
    SHOW(M.reachable(2, 0), "有向なので 2 から 0 へは行けない");

    api("【3】 負の辺と負閉路", "bool / Cost", "負の辺があってもよい。負閉路を通れる組は −∞（-DIST_INF）。W.negative_cycle で有無", "O(N^3)");
    RUN(WarshallFloyd Q(3));
    RUN(Q.add_edge(0, 1, 1); Q.add_edge(1, 2, -3); Q.add_edge(2, 1, 1));
    RUN(Q.run());
    SHOW(Q.negative_cycle, "1 → 2 → 1 が -2 の負閉路");
    show("Q.matrix()[0]", Row2{Q.matrix()[0]}, "0 から: 1 と 2 は −∞");
    SHOW(Q.neg_inf(0, 2), "");

    api("【4】 W.update_edge(u, v, w)", "なし", "計算の後に辺を 1 本足して、全部の組を O(N^2) で更新する（負閉路ができないこと）", "O(N^2)");
    code("// 辺を消していくクエリ（ABC375 F）: 全部消した状態から始めて、クエリを逆順に見て辺を足す");
    RUN(W.update_undirected(0, 3, 1));
    SHOW(W[0][3], "新しい辺 0-3 (1)");
    SHOW(W[1][3], "1 → 0 → 3 で 4 + 1 = 5 より 1 → 2 → 3 の 3 が短い");
    SHOW(W.path(0, 2), "0 → 3 → 2 で 1 + 2 = 3");

    head("典型的な使い方");
    note("全点対の最短距離（N ≤ 400 程度）: WarshallFloyd W(G); W[i][j]");
    note("辺を消していくクエリ             : 逆順に update_edge で足していく");
    note("推移閉包（i から j に行けるか）  : 重み 1 で作って W.reachable(i, j)");
    note("経由地を全部回る（巡回）         : 距離行列 W[i][j] を求めてから bit DP");
}

}  // namespace wf_demo

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
    wf_test::random_test();
    wf_test::benchmark();
    wf_demo::demo();
    return 0;
}
