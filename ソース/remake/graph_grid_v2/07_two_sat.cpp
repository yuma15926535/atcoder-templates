// =============================================================================
//  2-SAT（「A または B」の条件を全部満たす真偽の割り当て）  TwoSat     --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 07 連結性
//  依存       : 00_graph_base, 07_scc（同梱済み）
//  変数       : 0-indexed の x_0 .. x_{n-1}（true / false）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  条件はすべて「(x_i == f) または (x_j == g)」の形に直して入れる（f, g は true / false）。
//  「ならば」「同時に成り立たない」「等しい」「違う」「高々 1 つ」は専用の関数がある。
//  中では「¬A → B」の含意グラフを作り、強連結成分分解（07_scc）で判定する。
//
//  機能一覧                                                                       計算量
//    TwoSat ts(n)                         変数 n 個で作る
//    ts.add_clause(i, f, j, g)            (x_i == f) または (x_j == g)                     O(1)
//    ts.add_imply(i, f, j, g)             (x_i == f) ならば (x_j == g)                     O(1)
//    ts.add_nand(i, f, j, g)              (x_i == f) と (x_j == g) は同時に成り立たない     O(1)
//    ts.must(i, f)                        x_i == f に決める                                O(1)
//    ts.add_equal(i, j) / ts.add_not_equal(i, j)   x_i == x_j / x_i != x_j                 O(1)
//    ts.at_most_one({{i1, f1}, {i2, f2}, ...})   成り立つのは高々 1 つ（補助の変数が増える）  O(k)
//    ts.add_var()                         変数を 1 つ増やす（番号を返す）                    O(1)
//    ts.satisfiable()                     全部満たせるか。満たせれば ts.answer[i] = 1 / 0    O(n + 条件の数)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤・SCC も入っている）
//    ts_test::random_test() / ts_test::benchmark() / ts_demo::demo()
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
#ifndef GG_SCC
#define GG_SCC
// ---- 強連結成分分解（グラフ・グリッド第二世代）  非再帰の Tarjan 法。成分番号はトポロジカル順 ----
// 互いに行き来できる頂点のまとまり（強連結成分）に分け、まとまりを 1 頂点に縮めた DAG を作る
struct SCC {
    int count = 0;                   // 成分の数
    std::vector<int> id;             // 成分番号 0 ~ count-1。トポロジカル順（辺 u → v なら id[u] <= id[v]）
    std::vector<int> start_, list_;  // 成分ごとの頂点（CSR）。group(c) で使う

    SCC() {}
    explicit SCC(const Graph& G) {
        build_core(G.size(), [&](int v, auto&& push) {
            for (const Edge& e : G[v]) push(e.to);
        });
    }
    // 隣接を自分で書く版。adj(v, push) の中で v から出る辺の行き先 u ごとに push(u)
    template <class F>
    void build_core(int n, F&& adj) {
        id.assign(n, -1), count = 0;
        std::vector<int> ord(n, -1), low(n, 0), st, buf;
        struct Frame { int v, begin, cur; };
        std::vector<Frame> cs;
        int k = 0;
        for (int s = 0; s < n; ++s) {
            if (ord[s] != -1) continue;
            auto enter = [&](int v) {
                ord[v] = low[v] = k++, st.push_back(v);
                const int b = (int)buf.size();
                adj(v, [&](int u) { buf.push_back(u); });
                cs.push_back({v, b, b});
            };
            enter(s);
            while (!cs.empty()) {
                Frame& f = cs.back();
                const int v = f.v;
                if (f.cur < (int)buf.size()) {
                    const int u = buf[f.cur++];
                    if (ord[u] == -1) enter(u);
                    else if (id[u] == -1) low[v] = std::min(low[v], ord[u]);  // u はスタック上
                    continue;
                }
                buf.resize(f.begin);
                cs.pop_back();
                if (low[v] == ord[v]) {  // v が成分の根
                    while (true) {
                        const int x = st.back();
                        st.pop_back();
                        id[x] = count;
                        if (x == v) break;
                    }
                    ++count;
                }
                if (!cs.empty()) low[cs.back().v] = std::min(low[cs.back().v], low[v]);
            }
        }
        // Tarjan 法は「行き先側の成分」から番号が付く（トポロジカル順の逆）ので反転する
        for (int& x : id) x = count - 1 - x;
        start_.assign(count + 1, 0);
        for (int v = 0; v < n; ++v) ++start_[id[v] + 1];
        for (int c = 0; c < count; ++c) start_[c + 1] += start_[c];
        list_.assign(n, 0);
        std::vector<int> pos(start_.begin(), start_.end() - 1);
        for (int v = 0; v < n; ++v) list_[pos[id[v]]++] = v;
    }
    IntSpan group(int c) const { return {list_.data() + start_[c], list_.data() + start_[c + 1]}; }  // 成分 c の頂点（昇順）
    int size_of(int c) const { return start_[c + 1] - start_[c]; }
    bool same(int u, int v) const { return id[u] == id[v]; }
    std::vector<int> sizes() const {
        std::vector<int> r(count);
        for (int c = 0; c < count; ++c) r[c] = size_of(c);
        return r;
    }
    std::vector<std::vector<int>> groups() const {
        std::vector<std::vector<int>> r(count);
        for (int c = 0; c < count; ++c) r[c] = group(c).to_vector();
        return r;
    }
    // 縮約した DAG で使う、元のグラフの辺の番号（成分をまたぐ辺だけ。unique なら同じ成分の組は最初の 1 本だけ）
    std::vector<int> dag_edges(const Graph& G, bool unique = true) const {
        std::vector<int> res;
        std::vector<int> last(count, -1);  // last[cv] = 直前にこの行き先を使った始点の成分（重複の判定）
        for (int c = 0; c < count; ++c)
            for (int v : group(c))
                for (const Edge& e : G[v]) {
                    const int d = id[e.to];
                    if (d == c) continue;
                    if (unique) {
                        if (last[d] == c) continue;
                        last[d] = c;
                    }
                    res.push_back(e.id);
                }
        return res;
    }
    // 縮約した DAG（頂点 = 成分、辺は dag_edges(G, unique) の順。重みはそのまま）。辺はすべて番号の小さい成分 → 大きい成分
    Graph dag(const Graph& G, bool unique = true) const {
        Graph D(count);
        for (int idx : dag_edges(G, unique)) {
            const RawEdge& e = G.es[idx];
            // 無向辺は同じ成分になるので、ここに来るのは有向辺（G[v] の向き = es の向き）
            D.add_edge(id[e.from], id[e.to], e.w);
        }
        return D;
    }
};
// 隣接を自分で書く版（グリッドや状態のグラフ）
template <class F>
SCC scc_core(int n, F&& adj) {
    SCC s;
    s.build_core(n, adj);
    return s;
}
#endif  // GG_SCC
#ifndef GG_TWO_SAT
#define GG_TWO_SAT
// ---- 2-SAT（グラフ・グリッド第二世代）  「(条件 A) または (条件 B)」の形の条件を全部満たす真偽の割り当て ----
// 変数 x_0 .. x_{n-1}（true / false）。条件は「x_i == f」の形（f は true / false）を 2 つ OR でつないだもの
struct TwoSat {
    int n = 0;                                // 変数の数（at_most_one で補助の変数が増える）
    std::vector<int> answer;                  // satisfiable() が true のときの割り当て（1 = true / 0 = false）
    std::vector<std::pair<int, int>> imp_;    // 含意の辺（リテラルの番号）

    TwoSat() {}
    explicit TwoSat(int n_) : n(n_) {}
    int add_var() { return n++; }  // 変数を 1 つ増やし、その番号を返す
    static int lit(int i, bool f) { return 2 * i + (f ? 0 : 1); }  // 「x_i == f」のリテラルの番号（否定は ^ 1）
    // (x_i == f) または (x_j == g)
    void add_clause(int i, bool f, int j, bool g) {
        const int a = lit(i, f), b = lit(j, g);
        imp_.emplace_back(a ^ 1, b), imp_.emplace_back(b ^ 1, a);  // ¬A → B、¬B → A
    }
    void add_imply(int i, bool f, int j, bool g) { add_clause(i, !f, j, g); }       // (x_i == f) ならば (x_j == g)
    void must(int i, bool f) { add_clause(i, f, i, f); }                           // x_i == f に決める
    void add_nand(int i, bool f, int j, bool g) { add_clause(i, !f, j, !g); }       // (x_i == f) と (x_j == g) は同時に成り立たない
    void add_equal(int i, int j) { add_clause(i, true, j, false), add_clause(i, false, j, true); }      // x_i == x_j
    void add_not_equal(int i, int j) { add_clause(i, true, j, true), add_clause(i, false, j, false); }  // x_i != x_j
    // 「x_i == f」の形の条件のうち、成り立つのは高々 1 つ（補助の変数を使って O(k) 個の条件にする）
    void at_most_one(const std::vector<std::pair<int, bool>>& lits) {
        const int k = (int)lits.size();
        if (k <= 1) return;
        int prev = -1;  // prev = 「ここまでのどれかが成り立つ」を表す補助の変数
        for (int t = 0; t < k; ++t) {
            const int i = lits[t].first;
            const bool f = lits[t].second;
            const int s = add_var();
            add_imply(i, f, s, true);  // lit_t → s_t
            if (prev >= 0) {
                add_imply(prev, true, s, true);   // s_{t-1} → s_t
                add_nand(prev, true, i, f);        // s_{t-1} → ¬lit_t
            }
            prev = s;
        }
    }
    // 全部の条件を満たす割り当てがあるか。あれば answer に入れる。O(n + 条件の数)
    bool satisfiable() {
        Graph G(2 * n);
        for (const auto& e : imp_) G.add_edge(e.first, e.second);
        const SCC S(G);
        answer.assign(n, 0);
        for (int i = 0; i < n; ++i) {
            const int t = S.id[lit(i, true)], f = S.id[lit(i, false)];
            if (t == f) {
                answer.clear();
                return false;
            }
            answer[i] = (t > f);  // トポロジカル順で後ろのほうを真にする
        }
        return true;
    }
};
#endif  // GG_TWO_SAT
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "two_sat";
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

namespace ts_test {

struct Cond {
    int type;  // 0 clause, 1 imply, 2 nand, 3 must, 4 equal, 5 not_equal, 6 at_most_one
    int i, j;
    bool f, g;
    std::vector<std::pair<int, bool>> lits;
};
static bool holds(const Cond& c, const std::vector<int>& x) {
    const bool a = (x[c.i] != 0) == c.f, b = (x[c.j] != 0) == c.g;
    switch (c.type) {
        case 0: return a || b;
        case 1: return !a || b;
        case 2: return !(a && b);
        case 3: return a;
        case 4: return x[c.i] == x[c.j];
        case 5: return x[c.i] != x[c.j];
        default: {
            int cnt = 0;
            for (auto [v, f] : c.lits) cnt += ((x[v] != 0) == f);
            return cnt <= 1;
        }
    }
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    int sat_cnt = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 9), q = rnd(rng, 0, 3 * n);
        TwoSat ts(n);
        std::vector<Cond> cs;
        for (int k = 0; k < q; ++k) {
            Cond c{rnd(rng, 0, 6), rnd(rng, 0, n - 1), rnd(rng, 0, n - 1), (rng() & 1) != 0, (rng() & 1) != 0, {}};
            switch (c.type) {
                case 0: ts.add_clause(c.i, c.f, c.j, c.g); break;
                case 1: ts.add_imply(c.i, c.f, c.j, c.g); break;
                case 2: ts.add_nand(c.i, c.f, c.j, c.g); break;
                case 3: ts.must(c.i, c.f); break;
                case 4: ts.add_equal(c.i, c.j); break;
                case 5: ts.add_not_equal(c.i, c.j); break;
                default:
                    for (int t = rnd(rng, 0, 4); t > 0; --t) c.lits.emplace_back(rnd(rng, 0, n - 1), (rng() & 1) != 0);
                    ts.at_most_one(c.lits);
            }
            cs.push_back(c);
        }
        bool ex = false;
        for (int mask = 0; mask < (1 << n) && !ex; ++mask) {
            std::vector<int> x(n);
            for (int i = 0; i < n; ++i) x[i] = mask >> i & 1;
            bool ok = true;
            for (const Cond& c : cs) ok &= holds(c, x);
            ex |= ok;
        }
        const bool got = ts.satisfiable();
        CHECK(got == ex);
        if (got) {
            ++sat_cnt;
            CHECK((int)ts.answer.size() == ts.n && ts.n >= n);
            for (const Cond& c : cs) CHECK(holds(c, ts.answer));
        } else {
            CHECK(ts.answer.empty());
        }
    }
    CHECK(sat_cnt > cases / 10 && sat_cnt < cases);  // 満たせる・満たせないの両方を十分に試している
    std::cout << "[two_sat] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[two_sat] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        // 答えを 1 つ決めて、それを満たす条件だけを入れる（満たせる最大級の入力）
        std::vector<int> ans(n);
        for (int& a : ans) a = (int)(rng() & 1);
        TwoSat ts(n);
        for (int k = 0; k < n; ++k) {
            const int i = (int)(rng() % n), j = (int)(rng() % n);
            const bool f = (rng() & 1) != 0;
            ts.add_clause(i, f, j, ans[j] != 0);  // 後ろが必ず成り立つ
        }
        auto t0 = std::chrono::steady_clock::now();
        const bool ok = ts.satisfiable();
        std::cout << "[two_sat] benchmark: 変数 " << n << ", 条件 " << n << "  " << ms_since(t0) << " ms  (satisfiable " << ok << ")\n";
    }
}

}  // namespace ts_test


namespace ts_demo {
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
    banner("2-SAT 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"TwoSat ts(n);", "", "変数 x_0 .. x_{n-1}"},
        {"ts.add_clause(i, f, j, g)", "", "(x_i == f) または (x_j == g)"},
        {"ts.add_imply(i, f, j, g)", "", "(x_i == f) ならば (x_j == g)"},
        {"ts.add_nand(i, f, j, g)", "", "(x_i == f) と (x_j == g) は同時に成り立たない"},
        {"ts.must(i, f)", "", "x_i == f に決める"},
        {"ts.add_equal(i, j) / add_not_equal(i, j)", "", "同じ / 違う"},
        {"ts.at_most_one({{i, f}, ...})", "", "成り立つのは高々 1 つ"},
        {"ts.satisfiable() / ts.answer", "bool / vector<int>", "満たせるか / 割り当て"},
    };
    cheat_sheet("API 一覧", rows, 42, 20);

    api("【1】 add_clause と satisfiable", "bool", "「(x_i == f) または (x_j == g)」を入れて、全部を満たす割り当てを探す", "O(n + 条件の数)");
    RUN(TwoSat ts(3));
    RUN(ts.add_clause(0, true, 1, true));    // x0 または x1
    RUN(ts.add_clause(0, false, 2, true));   // ¬x0 または x2
    RUN(ts.add_clause(1, false, 2, false));  // ¬x1 または ¬x2
    SHOW(ts.satisfiable(), "満たせる");
    SHOW(ts.answer, "x0 = 1, x1 = 0, x2 = 1 など（1 = true）");

    api("【2】 いろいろな条件", "なし", "ならば・同時に成り立たない・等しい・違う・決める。どれも add_clause に直して入れる", "O(1)");
    RUN(TwoSat t2(3));
    RUN(t2.add_imply(0, true, 1, true));  // x0 ならば x1
    RUN(t2.add_not_equal(1, 2));          // x1 != x2
    RUN(t2.must(2, true));                // x2 = true
    SHOW(t2.satisfiable(), "");
    SHOW(t2.answer, "x2 = 1 → x1 = 0 → x0 = 0");
    RUN(t2.must(0, true));
    SHOW(t2.satisfiable(), "x0 = true にすると x1 = true が必要で矛盾");

    api("【3】 at_most_one", "なし", "いくつかの条件のうち、成り立つのは高々 1 つ（補助の変数を足して O(k) 個の条件にする）", "O(k)");
    RUN(TwoSat t3(4));
    RUN(t3.at_most_one({{0, true}, {1, true}, {2, true}, {3, true}}));
    RUN(t3.add_clause(0, true, 1, true));  // x0 または x1
    SHOW(t3.satisfiable(), "");
    SHOW(t3.answer, "x0 か x1 のどちらか 1 つだけが 1（後ろは補助の変数）");

    head("典型的な使い方（「2 択」がたくさんあり、2 つずつの組の条件だけなら 2-SAT）");
    code("// 旗の置き場所（ACL Practice H）: 旗 i を置き場 A に置く = x_i。近すぎる組は同時に選べない");
    code("//   for 各組 (i, a), (j, b) で距離 < D:  ts.add_nand(i, a, j, b);");
    code("// カードの表裏（ABC210 F）: 素因数 p ごとに「p で割り切れる面」を集めて ts.at_most_one(その面の一覧)");
    code("// 2 色塗り分けで「この 2 つは違う色」: add_not_equal(i, j)（ただの二部グラフ判定なら 01_bipartite）");
    note("3 つ以上の OR（x または y または z）は 2-SAT では書けない。");
}

}  // namespace ts_demo

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
    ts_test::random_test();
    ts_test::benchmark();
    ts_demo::demo();
    return 0;
}
