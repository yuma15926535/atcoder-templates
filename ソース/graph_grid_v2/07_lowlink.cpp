// =============================================================================
//  lowlink（橋・関節点・二重辺連結成分・二重頂点連結成分・消したときの成分）  LowLink
//                                                                 --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 07 連結性
//  依存       : 00_graph_base（同梱済み）
//  頂点番号   : 0-indexed。辺の向きは無視する（無向グラフ用）。多重辺・自己ループ・非連結に対応
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  橋     = 消すと連結成分が増える辺        関節点 = 消すと連結成分が増える頂点
//  二重辺連結成分 = 橋を全部消したときの連結成分（どの 1 辺を消してもつながったまま）
//  二重頂点連結成分（ブロック）= どの 1 頂点を消してもつながったままの、辺のまとまり
//  多重辺は辺番号で区別するので、2 本ある辺は橋にならない（DFS で「来た辺」を番号で判定）
//
//  機能一覧                                                                       計算量
//    LowLink L(G)                           前処理（非再帰の DFS 1 回）                O(N + M)
//    L.is_bridge[id] / L.bridges            橋か / 橋の辺番号の一覧
//    L.is_articulation[v] / L.articulations 関節点か / 関節点の一覧
//    L.components / L.comp[v] / L.comp_size[c]   連結成分の数 / 番号 / 頂点数
//    L.two_edge_count / L.two_edge_id[v]    二重辺連結成分の数 / 番号
//    L.bridge_tree()                        二重辺連結成分を橋で結んだ木（森）           O(N + M)
//    L.block_count / L.block_of_edge[id] / L.block_vertices(b)   ブロックの数 / 辺のブロック / ブロックの頂点
//    L.block_cut_tree()                     ブロックカット木（頂点 0 ~ N-1 = 元の頂点、N + b = ブロック b）
//    L.components_after_removing_vertex(v) / L.pieces_after_removing_vertex(v)   v を消したときの成分の数 / 各成分の頂点数
//    L.components_after_removing_edge(id) / L.pieces_after_removing_edge(id)    辺を消したとき          O(1)
//    L.ord / L.low / L.par / L.par_edge / L.sub   DFS の情報（自分で lowlink を使うとき）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    ll_test::random_test() / ll_test::benchmark() / ll_demo::demo()
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
#ifndef GG_LOWLINK
#define GG_LOWLINK
// ---- lowlink（グラフ・グリッド第二世代）  橋・関節点・二重辺連結成分・二重頂点連結成分（ブロック）・消したときの成分 ----
// 辺の向きは無視する（G.es の各辺を無向として扱う）。多重辺・自己ループ・非連結に対応。非再帰
struct LowLink {
    int n = 0, m = 0;
    int components = 0;                 // 連結成分の数
    std::vector<int> comp, comp_size;   // comp[v] = 連結成分の番号、comp_size[c] = 成分 c の頂点数
    std::vector<int> ord, low;          // DFS の訪問順・lowlink
    std::vector<int> par, par_edge;     // DFS 木の親と、親との辺の番号（根は -1）
    std::vector<int> sub;               // DFS 木での部分木の頂点数
    std::vector<int> is_bridge;         // is_bridge[辺番号] = 1（橋：消すと非連結になる辺）/ 0
    std::vector<int> is_articulation;   // is_articulation[v] = 1（関節点：消すと非連結になる頂点）/ 0
    std::vector<int> bridges;           // 橋の辺番号（昇順）
    std::vector<int> articulations;     // 関節点（昇順）
    int two_edge_count = 0;             // 二重辺連結成分の数（橋を消したときの連結成分）
    std::vector<int> two_edge_id;       // 頂点 v の二重辺連結成分の番号
    int block_count = 0;                // 二重頂点連結成分（ブロック）の数。辺を持たない孤立点はブロックに入らない
    std::vector<int> block_of_edge;     // 辺がどのブロックか（各辺はちょうど 1 つのブロック。自己ループは単独のブロック）
    std::vector<int> bstart_, blist_;   // ブロックの頂点（CSR）
    std::vector<int> start_, adj_;      // 無向の隣接（CSR。辺番号）
    std::vector<RawEdge> es_;

    LowLink() {}
    explicit LowLink(const Graph& G) { build(G); }
    void build(const Graph& G) {
        n = G.size(), m = G.edge_count(), es_ = G.es;
        start_.assign(n + 1, 0), adj_.assign(2 * (size_t)m, 0);
        for (const RawEdge& e : es_) ++start_[e.from + 1], ++start_[e.to + 1];
        for (int v = 0; v < n; ++v) start_[v + 1] += start_[v];
        {
            std::vector<int> pos(start_.begin(), start_.end() - 1);
            for (int id = 0; id < m; ++id) adj_[pos[es_[id].from]++] = id, adj_[pos[es_[id].to]++] = id;
        }
        ord.assign(n, -1), low.assign(n, 0), par.assign(n, -1), par_edge.assign(n, -1), sub.assign(n, 1);
        comp.assign(n, -1), comp_size.clear(), components = 0;
        is_bridge.assign(m, 0), is_articulation.assign(n, 0), bridges.clear(), articulations.clear();
        block_of_edge.assign(m, -1), block_count = 0;
        std::vector<int> it(n), child_cnt(n, 0), est, bverts, mark(n, -1);
        std::vector<int> bsz;  // ブロックごとの頂点数（あとで CSR に）
        std::vector<int> st;
        int k = 0;
        auto close_block = [&](int stop_edge) {  // stop_edge まで辺のスタックから取り出して 1 つのブロックにする
            const int b = block_count++;
            int cnt = 0;
            while (true) {
                const int id = est.back();
                est.pop_back();
                block_of_edge[id] = b;
                for (int x : {es_[id].from, es_[id].to})
                    if (mark[x] != b) mark[x] = b, bverts.push_back(x), ++cnt;
                if (id == stop_edge) break;
            }
            bsz.push_back(cnt);
        };
        for (int s = 0; s < n; ++s) {
            if (ord[s] != -1) continue;
            const int c = components++;
            comp_size.push_back(0);
            ord[s] = low[s] = k++, comp[s] = c, it[s] = start_[s], st.push_back(s);
            while (!st.empty()) {
                const int v = st.back();
                if (it[v] < start_[v + 1]) {
                    const int id = adj_[it[v]++];
                    const int u = es_[id].from == v ? es_[id].to : es_[id].from;
                    if (id == par_edge[v]) continue;
                    if (u == v) {  // 自己ループ: それだけで 1 つのブロック（隣接に 2 回出てくるので 1 回目だけ）
                        if (block_of_edge[id] == -1) est.push_back(id), close_block(id);
                        continue;
                    }
                    if (ord[u] == -1) {
                        par[u] = v, par_edge[u] = id, comp[u] = c, ++child_cnt[v];
                        ord[u] = low[u] = k++, it[u] = start_[u], est.push_back(id), st.push_back(u);
                    } else if (ord[u] < ord[v]) {  // 祖先への後退辺
                        low[v] = std::min(low[v], ord[u]), est.push_back(id);
                    }
                    continue;
                }
                st.pop_back();
                ++comp_size[c];
                const int p = par[v];
                if (p == -1) continue;
                low[p] = std::min(low[p], low[v]), sub[p] += sub[v];
                if (low[v] > ord[p]) is_bridge[par_edge[v]] = 1;
                if (low[v] >= ord[p]) {
                    if (par[p] != -1) is_articulation[p] = 1;
                    close_block(par_edge[v]);
                }
            }
            if (child_cnt[s] >= 2) is_articulation[s] = 1;
        }
        for (int id = 0; id < m; ++id)
            if (is_bridge[id]) bridges.push_back(id);
        for (int v = 0; v < n; ++v)
            if (is_articulation[v]) articulations.push_back(v);
        bstart_.assign(block_count + 1, 0);
        for (int b = 0; b < block_count; ++b) bstart_[b + 1] = bstart_[b] + bsz[b];
        blist_ = bverts;
        for (int b = 0; b < block_count; ++b) std::sort(blist_.begin() + bstart_[b], blist_.begin() + bstart_[b + 1]);
        // 二重辺連結成分: 橋を通らずに塗る
        two_edge_id.assign(n, -1), two_edge_count = 0;
        for (int s = 0; s < n; ++s) {
            if (two_edge_id[s] != -1) continue;
            const int c = two_edge_count++;
            two_edge_id[s] = c, st.push_back(s);
            while (!st.empty()) {
                const int v = st.back();
                st.pop_back();
                for (int q = start_[v]; q < start_[v + 1]; ++q) {
                    const int id = adj_[q];
                    if (is_bridge[id]) continue;
                    const int u = es_[id].from == v ? es_[id].to : es_[id].from;
                    if (two_edge_id[u] == -1) two_edge_id[u] = c, st.push_back(u);
                }
            }
        }
    }
    IntSpan block_vertices(int b) const { return {blist_.data() + bstart_[b], blist_.data() + bstart_[b + 1]}; }  // ブロック b の頂点（昇順）
    // 橋で結んだ木（森）: 頂点 = 二重辺連結成分、辺 k = bridges[k]（重みはそのまま）
    Graph bridge_tree() const {
        Graph T(two_edge_count);
        for (int id : bridges) T.add_undirected(two_edge_id[es_[id].from], two_edge_id[es_[id].to], es_[id].w);
        return T;
    }
    // ブロックカット木（森）: 頂点 0 ~ n-1 = 元の頂点、n + b = ブロック b。v がブロック b に入っていれば辺 v - (n + b)
    Graph block_cut_tree() const {
        Graph T(n + block_count);
        for (int b = 0; b < block_count; ++b)
            for (int v : block_vertices(b)) T.add_undirected(v, n + b);
        return T;
    }
    // 頂点 v（と v の辺）を消したときのグラフ全体の連結成分の数
    int components_after_removing_vertex(int v) const { return components - 1 + (int)pieces_after_removing_vertex(v).size(); }
    // 頂点 v を消したとき、v の連結成分だった部分がいくつに分かれ、それぞれ何頂点か（v が孤立点なら空）
    std::vector<int> pieces_after_removing_vertex(int v) const {
        std::vector<int> r;
        int rest = comp_size[comp[v]] - 1;
        for (int q = start_[v]; q < start_[v + 1]; ++q) {
            const int id = adj_[q];
            const int u = es_[id].from == v ? es_[id].to : es_[id].from;
            if (u != v && par[u] == v && par_edge[u] == id && low[u] >= ord[v]) r.push_back(sub[u]), rest -= sub[u];
        }
        if (rest > 0) r.push_back(rest);
        return r;
    }
    // 辺 id を消したときのグラフ全体の連結成分の数
    int components_after_removing_edge(int id) const { return components + is_bridge[id]; }
    // 辺 id を消したとき、その辺の連結成分だった部分の各成分の頂点数（橋なら 2 つ）
    std::vector<int> pieces_after_removing_edge(int id) const {
        const int a = es_[id].from, b = es_[id].to;
        const int total = comp_size[comp[a]];
        if (!is_bridge[id]) return {total};
        const int c = (par[a] == b && par_edge[a] == id) ? a : b;  // 子の側
        return {sub[c], total - sub[c]};
    }
};
#endif  // GG_LOWLINK
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "lowlink";
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

namespace ll_test {

// 愚直: 消す頂点・辺を指定して連結成分の頂点数を数える
static std::vector<int> pieces(const Graph& G, int skip_v, int skip_e) {
    const int n = G.size();
    std::vector<std::vector<int>> adj(n);
    for (int id = 0; id < G.edge_count(); ++id) {
        const RawEdge& e = G.es[id];
        if (id == skip_e || e.from == skip_v || e.to == skip_v) continue;
        adj[e.from].push_back(e.to), adj[e.to].push_back(e.from);
    }
    std::vector<int> seen(n, 0), res;
    for (int s = 0; s < n; ++s) {
        if (s == skip_v || seen[s]) continue;
        int cnt = 0;
        std::vector<int> st{s};
        seen[s] = 1;
        while (!st.empty()) {
            const int v = st.back();
            st.pop_back();
            ++cnt;
            for (int u : adj[v])
                if (!seen[u]) seen[u] = 1, st.push_back(u);
        }
        res.push_back(cnt);
    }
    return res;
}
static std::vector<int> labels(const Graph& G, const std::vector<int>& use_edge) {
    const int n = G.size();
    std::vector<int> lab(n);
    std::iota(lab.begin(), lab.end(), 0);
    for (bool ch = true; ch;) {
        ch = false;
        for (int id = 0; id < G.edge_count(); ++id) {
            if (!use_edge[id]) continue;
            const int a = G.es[id].from, b = G.es[id].to, mn = std::min(lab[a], lab[b]);
            if (lab[a] != mn || lab[b] != mn) lab[a] = lab[b] = mn, ch = true;
        }
    }
    return lab;
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 2500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 7), m = rnd(rng, 0, 9);
        Graph G(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            (rng() % 4 == 0) ? G.add_edge(a, b) : G.add_undirected(a, b);  // 有向で渡しても無向として扱う
        }
        const LowLink L(G);
        const std::vector<int> base = pieces(G, -1, -1);
        const int C = (int)base.size();
        CHECK(L.components == C);
        // 頂点を消す
        for (int v = 0; v < n; ++v) {
            std::vector<int> p = pieces(G, v, -1), q = L.pieces_after_removing_vertex(v);
            CHECK(L.components_after_removing_vertex(v) == (int)p.size());
            // q は v の成分だった部分だけ。p から「v を含まない成分（そのまま残る）」を除くと一致する
            std::vector<int> others;
            {
                std::vector<int> all = base;
                const int mine = L.comp_size[L.comp[v]];
                all.erase(std::find(all.begin(), all.end(), mine));
                others = all;
            }
            for (int x : others) p.erase(std::find(p.begin(), p.end(), x));
            std::sort(p.begin(), p.end()), std::sort(q.begin(), q.end());
            CHECK(p == q);
            CHECK(L.is_articulation[v] == (L.components_after_removing_vertex(v) > C));
        }
        // 辺を消す
        for (int id = 0; id < m; ++id) {
            const std::vector<int> p = pieces(G, -1, id);
            CHECK(L.components_after_removing_edge(id) == (int)p.size() && L.is_bridge[id] == ((int)p.size() > C));
            std::vector<int> q = L.pieces_after_removing_edge(id), pp;
            const int a = G.es[id].from;
            // a の成分（消す前）の頂点だけを数え直す
            std::vector<int> use(m, 1);
            use[id] = 0;
            const std::vector<int> lab = labels(G, use);
            const std::vector<int> lab0 = labels(G, std::vector<int>(m, 1));
            std::vector<int> cnt(n, 0);
            for (int v = 0; v < n; ++v)
                if (lab0[v] == lab0[a]) ++cnt[lab[v]];
            for (int x : cnt)
                if (x > 0) pp.push_back(x);
            std::sort(pp.begin(), pp.end()), std::sort(q.begin(), q.end());
            CHECK(pp == q);
        }
        // 一覧
        std::vector<int> br, ar;
        for (int id = 0; id < m; ++id)
            if (L.is_bridge[id]) br.push_back(id);
        for (int v = 0; v < n; ++v)
            if (L.is_articulation[v]) ar.push_back(v);
        CHECK(L.bridges == br && L.articulations == ar);
        // 二重辺連結成分 = 橋以外の辺でつながっている
        {
            std::vector<int> use(m);
            for (int id = 0; id < m; ++id) use[id] = !L.is_bridge[id];
            const std::vector<int> lab = labels(G, use);
            int cnt = 0;
            for (int v = 0; v < n; ++v) cnt += (lab[v] == v);
            CHECK(L.two_edge_count == cnt);
            for (int a = 0; a < n; ++a)
                for (int b = 0; b < n; ++b) CHECK((L.two_edge_id[a] == L.two_edge_id[b]) == (lab[a] == lab[b]));
            const Graph T = L.bridge_tree();
            CHECK(T.size() == cnt && T.edge_count() == (int)L.bridges.size());
            for (int k = 0; k < T.edge_count(); ++k) {
                const RawEdge& e = G.es[L.bridges[k]];
                const RawEdge& t = T.es[k];
                CHECK(t.from == L.two_edge_id[e.from] && t.to == L.two_edge_id[e.to]);
            }
        }
        // ブロック: 2 本の辺が同じブロック ⇔ 両方を含む単純閉路がある（閉路に入らない辺は 1 本で 1 ブロック）
        {
            std::vector<int> grp(m);
            std::iota(grp.begin(), grp.end(), 0);
            auto root = [&](int x) {
                while (grp[x] != x) x = grp[x];
                return x;
            };
            for (int mask = 1; mask < (1 << m); ++mask) {
                std::vector<int> deg(n, 0);
                for (int id = 0; id < m; ++id)
                    if (mask >> id & 1) ++deg[G.es[id].from], ++deg[G.es[id].to];
                bool ok = true;
                for (int v = 0; v < n; ++v) ok &= (deg[v] == 0 || deg[v] == 2);
                if (!ok) continue;
                std::vector<int> use(m, 0);
                int first = -1;
                for (int id = 0; id < m; ++id)
                    if (mask >> id & 1) use[id] = 1, first = (first < 0 ? id : first);
                const std::vector<int> lab = labels(G, use);
                bool conn = true;
                for (int v = 0; v < n; ++v)
                    if (deg[v]) conn &= (lab[v] == lab[G.es[first].from]);
                if (!conn) continue;
                for (int id = 0; id < m; ++id)
                    if (mask >> id & 1) grp[root(id)] = root(first);
            }
            for (int a = 0; a < m; ++a)
                for (int b = 0; b < m; ++b) CHECK((L.block_of_edge[a] == L.block_of_edge[b]) == (root(a) == root(b)));
            std::set<int> used;
            for (int id = 0; id < m; ++id) used.insert(L.block_of_edge[id]);
            CHECK((int)used.size() == L.block_count && (m == 0 || (*used.begin() == 0 && *used.rbegin() == L.block_count - 1)));
            for (int b = 0; b < L.block_count; ++b) {
                std::set<int> vs;
                for (int id = 0; id < m; ++id)
                    if (L.block_of_edge[id] == b) vs.insert(G.es[id].from), vs.insert(G.es[id].to);
                CHECK(L.block_vertices(b).to_vector() == std::vector<int>(vs.begin(), vs.end()));
            }
            const Graph BT = L.block_cut_tree();
            int sum = 0;
            for (int b = 0; b < L.block_count; ++b) sum += L.block_vertices(b).size();
            CHECK(BT.size() == n + L.block_count && BT.edge_count() == sum);
        }
    }
    std::cout << "[lowlink] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[lowlink] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> us, vs;
        for (int v = 1; v < n; ++v) us.push_back(v), vs.push_back((int)(rng() % v));  // 木
        for (int i = 0; i < n / 2; ++i) us.push_back((int)(rng() % n)), vs.push_back((int)(rng() % n));
        const Graph G = graph_from_edges(n, us, vs, UNDIRECTED, 0);
        auto t0 = std::chrono::steady_clock::now();
        const LowLink L(G);
        const long long t1 = ms_since(t0);
        Graph P(n);  // 一直線（深さ n の DFS）
        for (int v = 0; v + 1 < n; ++v) P.add_undirected(v, v + 1);
        t0 = std::chrono::steady_clock::now();
        const LowLink LP(P);
        std::cout << "[lowlink] benchmark: N = " << n << ", M = " << G.edge_count() << "  " << t1 << " ms, 一直線 " << ms_since(t0) << " ms  (bridges " << L.bridges.size()
                  << ", " << LP.bridges.size() << ", blocks " << L.block_count << ")\n";
    }
}

}  // namespace ll_test


namespace ll_demo {
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
    banner("lowlink（橋・関節点・二重連結成分） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"LowLink L(G);", "", "前処理（辺の向きは無視）"},
        {"L.bridges / L.is_bridge[id]", "", "橋（消すと非連結になる辺）"},
        {"L.articulations / L.is_articulation[v]", "", "関節点（消すと非連結になる頂点）"},
        {"L.two_edge_id[v] / L.bridge_tree()", "", "二重辺連結成分 / それを橋で結んだ木"},
        {"L.block_of_edge[id] / L.block_vertices(b)", "", "ブロック（二重頂点連結成分）"},
        {"L.block_cut_tree()", "Graph", "ブロックカット木"},
        {"L.components_after_removing_vertex(v)", "int", "v を消したときの連結成分の数"},
        {"L.pieces_after_removing_vertex(v)", "vector<int>", "v を消したときに分かれる各部分の頂点数"},
        {"L.components_after_removing_edge(id)", "int", "辺を消したときの連結成分の数"},
    };
    cheat_sheet("API 一覧", rows, 44, 14);

    head("デモで使うグラフ（無向）: 三角形 0-1-2、橋 2-3、三角形 3-4-5、橋 5-6、孤立点 7");
    RUN(Graph G(8));
    RUN(G.add_undirected(0, 1); G.add_undirected(1, 2); G.add_undirected(2, 0); G.add_undirected(2, 3));
    RUN(G.add_undirected(3, 4); G.add_undirected(4, 5); G.add_undirected(5, 3); G.add_undirected(5, 6));
    RUN(LowLink L(G));

    api("【1】 橋と関節点", "vector<int>", "橋 = 消すと連結成分が増える辺、関節点 = 消すと連結成分が増える頂点", "O(N + M)");
    SHOW(L.bridges, "辺 3（2-3）と辺 7（5-6）");
    SHOW(L.articulations, "2, 3, 5");
    SHOW(L.is_bridge[0], "三角形の辺は橋ではない");

    api("【2】 頂点・辺を消したとき", "int / vector<int>", "その頂点（辺）だけを消したときの連結成分の数と、分かれた部分の頂点数", "O(次数) / O(1)");
    SHOW(L.components, "今の連結成分の数（孤立点 7 を含む）");
    SHOW(L.components_after_removing_vertex(2), "{0, 1} と {3, 4, 5, 6} に分かれる");
    SHOW(L.pieces_after_removing_vertex(2), "分かれた部分の頂点数");
    SHOW(L.components_after_removing_vertex(0), "0 は関節点ではない");
    SHOW(L.pieces_after_removing_edge(3), "橋 2-3 を消すと 3 頂点と 4 頂点");
    note("「頂点 v を消したら残りの成分の大きさは」（ABC334 G など）は pieces_after_removing_vertex。");

    api("【3】 二重辺連結成分と橋の木", "vector<int> / Graph", "橋を全部消したときの成分。橋で結ぶと木（森）になる", "O(N + M)");
    SHOW(L.two_edge_id, "{0,1,2}, {3,4,5}, {6}, {7}");
    SHOW(L.bridge_tree().es, "成分を橋で結んだ木（辺 k は bridges[k]）");
    note("「2 点の間をどの 1 本を消しても行き来できるか」= 同じ二重辺連結成分か。");

    api("【4】 ブロック（二重頂点連結成分）とブロックカット木", "int / IntSpan / Graph",
        "どの 1 頂点を消してもつながったままの辺のまとまり。関節点は複数のブロックに入る", "O(N + M)");
    SHOW(L.block_count, "三角形 2 つと橋 2 本");
    SHOW(L.block_of_edge, "各辺のブロック");
    SHOW(L.block_vertices(L.block_of_edge[0]), "三角形 0-1-2");
    SHOW(L.block_cut_tree().size(), "元の 8 頂点 + ブロック 4 つ");
    note("ブロックカット木の上の 2 頂点のパスにある関節点 = その 2 頂点の間で必ず通る頂点。");

    head("典型的な使い方");
    note("消すと非連結になる辺・頂点       : L.bridges / L.articulations");
    note("辺を 1 本消しても行き来できる組  : L.two_edge_id[u] == L.two_edge_id[v]");
    note("頂点を消したときの成分の数・大きさ: components_after_removing_vertex / pieces_after_removing_vertex");
    note("必ず通る頂点（無向グラフ）       : L.block_cut_tree() の上で 05_lca の path（関節点だけ数える）");
}

}  // namespace ll_demo

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
    ll_test::random_test();
    ll_test::benchmark();
    ll_demo::demo();
    return 0;
}
