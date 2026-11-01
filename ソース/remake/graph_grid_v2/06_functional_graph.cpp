// =============================================================================
//  Functional graph（各頂点から辺が 1 本: v → nxt[v]）  FunctionalGraph   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 06 特殊グラフ
//  依存       : 00_graph_base（同梱済み。Cost・IntSpan のため）
//  頂点番号   : 0-indexed（入力が 1-indexed なら read_next(N) で読むと 0-indexed になる）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  どの頂点からも、何回か進むとサイクルに入る（サイクル + そこへ流れ込む木）。
//  「K 回移動した先（K ≤ 10^18）」「何回で初めて v に着くか」「K 回移動したときの値の合計」などを答える。
//
//  高速化のポイント
//    - ダブリングの表（N × 60）を作らない。サイクルは位置の剰余で O(1)、木の部分は
//      「深さごとに行きがけ順で並べた配列」の二分探索で O(log N)。前処理 O(N)・メモリ O(N)
//
//  機能一覧                                                                       計算量
//    FunctionalGraph F(nxt) / F(nxt, w)   作る（w[v] = v → nxt[v] の重み）            O(N)
//    read_next(N, base = 1)               A_1 .. A_N を読んで 0-indexed の nxt にする   O(N)
//    F.kth(v, K)                          v から K 回進んだ頂点                         O(log N)
//    F.kth_sum(v, K)                      v から K 回進むときの重みの和（ABC241 E）      O(log N)
//    F.reach(u, v)                        u から何回で初めて v に着くか（着かなければ -1）  O(1)
//    F.on_cycle(v) / F.cycle_id[v] / F.cycle_pos[v]   サイクル上か / 入るサイクル / サイクル内の位置
//    F.to_cycle[v] / F.entry[v] / F.to_cycle_w[v]     サイクルまでの回数 / 最初に着くサイクルの頂点 / その重み
//    F.cycle_count() / F.cycle(c) / F.cycle_len(c) / F.cycle_len_of(v) / F.cycle_weight(c)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    fg_test::random_test() / fg_test::benchmark() / fg_demo::demo()
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
#ifndef GG_FUNCTIONAL_GRAPH
#define GG_FUNCTIONAL_GRAPH
// ---- Functional graph（グラフ・グリッド第二世代）  各頂点から辺が 1 本（v → nxt[v]）。「K 回移動」「周期」 ----
// どの頂点からも、いくつか進むとサイクルに入る（サイクル + そこへ流れ込む木）
struct FunctionalGraph {
    int n = 0;
    std::vector<int> nxt;        // 行き先（0-indexed）
    std::vector<Cost> w;         // w[v] = v → nxt[v] の重み（重みなしで作ったら全部 1）
    std::vector<int> cycle_id;   // v が入っていくサイクルの番号
    std::vector<int> cycle_pos;  // サイクル上の頂点ならサイクル内の位置（サイクルの先頭から何番目か）。そうでなければ -1
    std::vector<int> to_cycle;   // サイクルに着くまでの回数（サイクル上は 0）
    std::vector<int> entry;      // 最初に着くサイクル上の頂点（サイクル上なら自分）
    std::vector<Cost> to_cycle_w;  // サイクルに着くまでの重みの和
    std::vector<int> cyc_start_, cyc_list_;  // サイクルの頂点（CSR。cycle(c) で使う。cyc_list_ は v → nxt[v] の順）
    std::vector<Cost> cyc_pre_;              // サイクルごとの重みの累積（cycle_pos の順。各サイクルは長さ + 1 個）
    std::vector<int> pre_start_;             // cyc_pre_ の各サイクルの開始位置
    std::vector<int> tin_, tout_;            // サイクルに流れ込む木（根 = サイクル上の頂点）の行きがけ順
    std::vector<int> lv_start_, lv_list_;    // 深さ（to_cycle）ごとの頂点を tin 順に並べたもの

    FunctionalGraph() {}
    explicit FunctionalGraph(const std::vector<int>& next) { build(next, std::vector<Cost>(next.size(), 1)); }
    template <class W>
    FunctionalGraph(const std::vector<int>& next, const std::vector<W>& weight) {
        build(next, std::vector<Cost>(weight.begin(), weight.end()));
    }

    int cycle_count() const { return (int)cyc_start_.size() - 1; }
    IntSpan cycle(int c) const { return {cyc_list_.data() + cyc_start_[c], cyc_list_.data() + cyc_start_[c + 1]}; }
    int cycle_len(int c) const { return cyc_start_[c + 1] - cyc_start_[c]; }
    int cycle_len_of(int v) const { return cycle_len(cycle_id[v]); }   // v が入っていくサイクルの長さ
    bool on_cycle(int v) const { return cycle_pos[v] >= 0; }
    Cost cycle_weight(int c) const { return cyc_pre_[pre_start_[c] + cycle_len(c)]; }  // サイクル 1 周の重みの和

    int kth(int v, long long K) const {  // v から K 回進んだ頂点（K >= 0。K は 10^18 程度まで）
        if (K <= to_cycle[v]) return ancestor(v, (int)K);
        const int c = cycle_id[v], L = cycle_len(c);
        const long long p = (cycle_pos[entry[v]] + (K - to_cycle[v])) % L;
        return cyc_list_[cyc_start_[c] + (int)p];
    }
    Cost kth_sum(int v, long long K) const {  // v から K 回進むときの重みの和（long long に収まる前提）
        if (K <= to_cycle[v]) return to_cycle_w[v] - to_cycle_w[ancestor(v, (int)K)];
        const int c = cycle_id[v], L = cycle_len(c), b = pre_start_[c];
        const long long r = K - to_cycle[v], q = r / L;
        const int p = cycle_pos[entry[v]], rem = (int)(r % L);
        Cost s = to_cycle_w[v] + (Cost)q * cyc_pre_[b + L];
        if (p + rem <= L) s += cyc_pre_[b + p + rem] - cyc_pre_[b + p];
        else s += (cyc_pre_[b + L] - cyc_pre_[b + p]) + cyc_pre_[b + p + rem - L];
        return s;
    }
    long long reach(int u, int v) const {  // u から何回進むと初めて v に着くか（着かなければ -1）
        if (on_cycle(v)) {
            if (cycle_id[u] != cycle_id[v]) return -1;
            const int L = cycle_len_of(v);
            return to_cycle[u] + ((cycle_pos[v] - cycle_pos[entry[u]]) % L + L) % L;
        }
        if (entry[u] != entry[v] || to_cycle[u] < to_cycle[v]) return -1;
        if (tin_[v] <= tin_[u] && tin_[u] < tout_[v]) return to_cycle[u] - to_cycle[v];
        return -1;
    }

    void build(const std::vector<int>& next, std::vector<Cost> weight) {
        n = (int)next.size(), nxt = next, w = std::move(weight);
        cycle_id.assign(n, -1), cycle_pos.assign(n, -1), to_cycle.assign(n, 0), entry.assign(n, -1), to_cycle_w.assign(n, 0);
        cyc_start_.assign(1, 0), cyc_list_.clear(), cyc_pre_.clear(), pre_start_.clear();
        // 1. サイクルを見つける（色: 0 未訪問, 1 今たどっている, 2 済み）
        std::vector<char> col(n, 0);
        std::vector<int> walk;
        for (int s = 0; s < n; ++s) {
            if (col[s]) continue;
            walk.clear();
            int v = s;
            while (col[v] == 0) col[v] = 1, walk.push_back(v), v = nxt[v];
            if (col[v] == 1) {  // v から始まるサイクル
                const int c = cycle_count();
                pre_start_.push_back((int)cyc_pre_.size());
                cyc_pre_.push_back(0);
                int x = v, k = 0;
                do {
                    cycle_id[x] = c, cycle_pos[x] = k++, entry[x] = x, cyc_list_.push_back(x);
                    cyc_pre_.push_back(cyc_pre_.back() + w[x]);
                    x = nxt[x];
                } while (x != v);
                cyc_start_.push_back((int)cyc_list_.size());
            }
            for (int x : walk) col[x] = 2;
        }
        // 2. 逆向きの辺（サイクル上の頂点への辺は、サイクルの外から来るものだけ）
        std::vector<int> rs(n + 1, 0);
        for (int v = 0; v < n; ++v)
            if (!on_cycle(v)) ++rs[nxt[v] + 1];
        for (int v = 0; v < n; ++v) rs[v + 1] += rs[v];
        std::vector<int> rl(rs[n]), pos(rs.begin(), rs.end() - 1);
        for (int v = 0; v < n; ++v)
            if (!on_cycle(v)) rl[pos[nxt[v]]++] = v;
        // 3. サイクル上の各頂点を根に、流れ込む木を DFS（非再帰）。行きがけ順と深さ
        tin_.assign(n, 0), tout_.assign(n, 0);
        std::vector<int> order, st;
        order.reserve(n);
        int maxd = 0;
        for (int r : cyc_list_) {
            st.push_back(r);
            while (!st.empty()) {
                const int v = st.back();
                if (v >= 0) {
                    tin_[v] = (int)order.size(), order.push_back(v);
                    st.back() = ~v;  // 帰りがけの印
                    for (int k = rs[v]; k < rs[v + 1]; ++k) {
                        const int u = rl[k];
                        to_cycle[u] = to_cycle[v] + 1, to_cycle_w[u] = to_cycle_w[v] + w[u], entry[u] = entry[v], cycle_id[u] = cycle_id[v];
                        maxd = std::max(maxd, to_cycle[u]);
                        st.push_back(u);
                    }
                } else {
                    tout_[~v] = (int)order.size(), st.pop_back();
                }
            }
        }
        // 4. 深さごとに tin の順で並べる（祖先を二分探索で求める）
        lv_start_.assign(maxd + 2, 0);
        for (int v = 0; v < n; ++v) ++lv_start_[to_cycle[v] + 1];
        for (int d = 0; d <= maxd; ++d) lv_start_[d + 1] += lv_start_[d];
        lv_list_.assign(n, 0);
        std::vector<int> lp(lv_start_.begin(), lv_start_.end() - 1);
        for (int v : order) lv_list_[lp[to_cycle[v]]++] = v;
    }

  private:
    int ancestor(int v, int k) const {  // 流れ込む木で v から k 回進んだ頂点（k <= to_cycle[v]）
        if (k == 0) return v;
        const int d = to_cycle[v] - k;
        const int* b = lv_list_.data() + lv_start_[d];
        const int* e = lv_list_.data() + lv_start_[d + 1];
        // 深さ d の頂点のうち、tin が tin[v] 以下で最大のもの
        int lo = 0, hi = (int)(e - b);
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            if (tin_[b[mid]] <= tin_[v]) lo = mid;
            else hi = mid;
        }
        return b[lo];
    }
};
// 「A_1 .. A_N（i の行き先）」を読む（base = 入力の番号の始まり）
inline std::vector<int> read_next(int n, int base = 1) {
    std::vector<int> a(n);
    for (int& x : a) std::cin >> x, x -= base;
    return a;
}
#endif  // GG_FUNCTIONAL_GRAPH
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "functional_graph";
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

namespace fg_test {

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 15);
        std::vector<int> nxt(n);
        std::vector<long long> w(n);
        for (int v = 0; v < n; ++v) nxt[v] = rnd(rng, 0, n - 1), w[v] = rnd(rng, -5, 9);
        if (it % 4 == 0) std::iota(nxt.begin(), nxt.end(), 0), std::shuffle(nxt.begin(), nxt.end(), rng);  // 順列
        const FunctionalGraph F(nxt, w), F1(nxt);
        // 構造: 愚直にたどる
        for (int v = 0; v < n; ++v) {
            std::vector<int> seen(n, -1), path;
            int x = v;
            while (seen[x] == -1) seen[x] = (int)path.size(), path.push_back(x), x = nxt[x];
            const int tail = seen[x], L = (int)path.size() - tail;  // x = 最初に 2 回目に来た頂点 = サイクルの入口
            CHECK(F.to_cycle[v] == tail && F.entry[v] == x && F.cycle_len_of(v) == L && F.on_cycle(v) == (tail == 0));
            Cost tw = 0;
            for (int k = 0; k < tail; ++k) tw += w[path[k]];
            CHECK(F.to_cycle_w[v] == tw);
            CHECK(F.cycle_id[v] == F.cycle_id[x]);
            // kth / kth_sum: 小さい K は愚直、大きい K は周期で縮めて比べる
            for (int t = 0; t < 8; ++t) {
                long long K = (t < 5) ? rnd(rng, 0, 3 * n) : (long long)(rng() % 1000000000000000000ULL);
                long long Kr = K;
                if (K > tail + L) Kr = tail + (K - tail) % L;
                int y = v;
                Cost s = 0;
                for (long long k = 0; k < Kr; ++k) s += w[y], y = nxt[y];
                CHECK(F.kth(v, K) == y && F1.kth(v, K) == y);
                // 和: 何周したかの分を足す
                Cost cyc = 0;
                for (int k = 0, z = x; k < L; ++k) cyc += w[z], z = nxt[z];
                const long long loops = (K > tail + L) ? (K - tail) / L - (Kr - tail) / L : 0;
                CHECK(F.kth_sum(v, K) == s + (Cost)loops * cyc);
                CHECK(F1.kth_sum(v, K) == K);
            }
            // reach: 2n 回までたどって最初に着く回数
            for (int u = 0; u < n; ++u) {
                long long ex = -1;
                int z = v;
                for (int k = 0; k <= 2 * n; ++k) {
                    if (z == u) {
                        ex = k;
                        break;
                    }
                    z = nxt[z];
                }
                CHECK(F.reach(v, u) == ex);
            }
        }
        // サイクルの一覧
        int on = 0;
        for (int c = 0; c < F.cycle_count(); ++c) {
            const IntSpan cy = F.cycle(c);
            CHECK(cy.size() == F.cycle_len(c));
            Cost cw = 0;
            for (int k = 0; k < cy.size(); ++k) {
                CHECK(nxt[cy[k]] == cy[(k + 1) % cy.size()] && F.cycle_pos[cy[k]] == k && F.cycle_id[cy[k]] == c);
                cw += w[cy[k]];
            }
            CHECK(F.cycle_weight(c) == cw);
            on += cy.size();
        }
        int on2 = 0;
        for (int v = 0; v < n; ++v) on2 += F.on_cycle(v);
        CHECK(on == on2);
    }
    {
        CinFromString cf("2 3 1\n");
        CHECK(read_next(3) == (std::vector<int>{1, 2, 0}));
    }
    std::cout << "[functional_graph] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[functional_graph] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> nxt(n);
        std::vector<long long> w(n);
        for (int v = 0; v < n; ++v) nxt[v] = (int)(rng() % n), w[v] = (long long)(rng() % 1000);
        auto t0 = std::chrono::steady_clock::now();
        const FunctionalGraph F(nxt, w);
        const long long tb = ms_since(t0);
        t0 = std::chrono::steady_clock::now();
        long long chk = 0;
        for (int q = 0; q < n; ++q) {
            const int v = (int)(rng() % n);
            const long long K = (long long)(rng() % 1000000000000000000ULL);
            chk += F.kth(v, K) + F.kth_sum(v, K % 1000000) + F.reach(v, (int)(rng() % n));
        }
        std::cout << "[functional_graph] benchmark: N = " << n << "  build " << tb << " ms, クエリ " << n << " 回（kth + kth_sum + reach）" << ms_since(t0)
                  << " ms  (check " << chk % 1000 << ")\n";
    }
}

}  // namespace fg_test


namespace fg_demo {
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
    banner("Functional graph 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"FunctionalGraph F(nxt) / F(nxt, w);", "", "v → nxt[v]（重み w[v]）"},
        {"read_next(N, base = 1)", "vector<int>", "A_1 .. A_N を読んで 0-indexed に"},
        {"F.kth(v, K) / F.kth_sum(v, K)", "int / Cost", "K 回進んだ頂点 / その間の重みの和"},
        {"F.reach(u, v)", "long long", "u から何回で初めて v に着くか（-1 なら着かない）"},
        {"F.on_cycle(v) / F.cycle_id[v] / F.cycle_pos[v]", "", "サイクル上か / 入るサイクル / 位置"},
        {"F.to_cycle[v] / F.entry[v]", "int", "サイクルまでの回数 / 最初に着くサイクルの頂点"},
        {"F.cycle(c) / F.cycle_len(c) / F.cycle_count()", "", "サイクルの頂点 / 長さ / 個数"},
    };
    cheat_sheet("API 一覧", rows, 48, 14);

    head("デモで使うグラフ: 0→1, 1→2, 2→3, 3→1, 4→2, 5→5（1→2→3→1 がサイクル、5 は自分に戻る）");
    RUN(std::vector<int> nxt = {1, 2, 3, 1, 2, 5});
    RUN(FunctionalGraph F(nxt));

    api("【1】 F.kth(v, K)", "int", "v から K 回進んだ頂点。K は 10^18 程度まで", "O(log N)");
    SHOW(F.kth(0, 1), "0 → 1");
    SHOW(F.kth(0, 4), "0 → 1 → 2 → 3 → 1");
    SHOW(F.kth(4, 1000000000000000000LL), "10^18 回でも一瞬");

    api("【2】 サイクルの情報", "", "各頂点が入っていくサイクル・サイクルまでの回数・最初に着くサイクルの頂点", "O(1)");
    SHOW(F.cycle_count(), "サイクルの数（1→2→3 と 5）");
    SHOW(F.cycle(0), "サイクル 0 の頂点（たどる順）");
    SHOW(F.to_cycle, "サイクルに着くまでの回数");
    SHOW(F.entry, "最初に着くサイクル上の頂点");
    SHOW(F.cycle_id, "入っていくサイクルの番号");
    SHOW(F.cycle_pos, "サイクル内の位置（サイクル上でなければ -1）");
    SHOW(F.cycle_len_of(0), "0 が入るサイクルの長さ");

    api("【3】 F.reach(u, v)", "long long", "u から何回進むと初めて v に着くか（着かなければ -1）", "O(1)");
    SHOW(F.reach(0, 3), "0 → 1 → 2 → 3");
    SHOW(F.reach(2, 1), "2 → 3 → 1（サイクルを回る）");
    SHOW(F.reach(1, 0), "サイクルから木には戻れない");
    SHOW(F.reach(4, 0), "別の枝には行けない");

    api("【4】 FunctionalGraph F(nxt, w) と F.kth_sum(v, K)", "Cost", "v から K 回進むときの重み w の和。サイクルは何周分かを掛け算で足す", "O(log N)");
    code("// ABC241 E: X に A[X mod N] を足す操作を K 回。頂点 i → (i + A[i]) mod N、重み A[i]");
    RUN(std::vector<long long> A = {2, 5, 7, 1, 3});
    RUN(int N = 5);
    RUN(std::vector<int> to(N));
    RUN(for (int i = 0; i < N; ++i) to[i] = (int)((i + A[i]) % N));
    RUN(FunctionalGraph FA(to, A));
    SHOW(FA.kth_sum(0, 3), "3 回で足される合計: 0 →(+2) 2 →(+7) 4 →(+3) 2 → 2 + 7 + 3");
    SHOW(FA.kth_sum(0, 1000000000000LL), "10^12 回でも一瞬");

    api("【5】 read_next(N, base = 1)", "vector<int>", "「A_1 A_2 ... A_N」（i の行き先、1-indexed）を読んで 0-indexed の nxt にする", "O(N)");
    {
        CinFromString cf("2 3 1\n");
        code("// 入力: 2 3 1");
        RUN(auto nx = read_next(3));
        SHOW(nx, "0-indexed になる");
    }
    note("順列（全頂点がサイクル上）なら 06_perm_cycles のほうが軽い（K 乗の順列・交換回数など）。");
}

}  // namespace fg_demo

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
    fg_test::random_test();
    fg_test::benchmark();
    fg_demo::demo();
    return 0;
}
