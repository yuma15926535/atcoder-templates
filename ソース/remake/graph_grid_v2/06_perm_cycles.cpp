// =============================================================================
//  順列のサイクル分解  PermCycles                                    --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 06 特殊グラフ
//  依存       : 00_graph_base（同梱済み。IntSpan のため）
//  頂点番号   : 0-indexed。p は 0 ~ n-1 の並べ替え（i → p[i] の矢印と見る。入力が 1-indexed なら read_perm(N)）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  順列は「いくつかの輪（サイクル）」に分かれる。K 回あてた先・周期・入れ替え回数はサイクルごとに決まる。
//  （一般の「各頂点から 1 本」なら 06_functional_graph。順列ならこちらのほうが軽くて機能が多い）
//
//  機能一覧                                                                       計算量
//    PermCycles C(p)                作る                                               O(N)
//    read_perm(N, base = 1)         P_1 .. P_N を読んで 0-indexed にする                 O(N)
//    C.count() / C.cycle(c) / C.cycle_len(c) / C.len_of(i)   サイクルの数 / 要素（i → p[i] の順）/ 長さ
//    C.cycle_id[i] / C.cycle_pos[i]  i のサイクル番号（小さい要素を含むサイクルから 0, 1, ..）/ 位置
//    C.kth(i, K)                    i に p を K 回あてた値（K < 0 なら逆向き。|K| ≤ 10^18）  O(1)
//    C.pow(K)                       p を K 回合成した順列（r[i] = kth(i, K)）            O(N)
//    C.inverse()                    逆順列                                             O(N)
//    C.min_swaps()                  任意の 2 つの入れ替えで整列する最小回数（N - サイクル数）  O(1)
//    C.sort_swaps()                 その手順（A = p を A[i] = i にする位置の組。ABC350 C）    O(N)
//    C.parity()                     偶置換 0 / 奇置換 1                                  O(1)
//    C.period(limit = 4e18)         何回あてると元に戻るか（長さの最小公倍数。limit 超えは -1）  O(N log N)
//    C.period_mod(mod)              周期を mod で（大きすぎるとき）                        O(N)
//    perm_compose(a, b)             合成 r[i] = a[b[i]]（先に b、次に a）                  O(N)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    pc_test::random_test() / pc_test::benchmark() / pc_demo::demo()
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
#ifndef GG_PERM_CYCLES
#define GG_PERM_CYCLES
// ---- 順列のサイクル分解（グラフ・グリッド第二世代）  p は 0 ~ n-1 の並べ替え（i → p[i]） ----
struct PermCycles {
    int n = 0;
    std::vector<int> p;
    std::vector<int> cycle_id;    // i が属するサイクルの番号（番号の小さい要素を含むサイクルから 0, 1, ...）
    std::vector<int> cycle_pos;   // サイクル内の位置（cycle(c)[cycle_pos[i]] == i）
    std::vector<int> start_, list_;  // サイクルの要素（CSR）。各サイクルは i → p[i] → p[p[i]] → ... の順
    PermCycles() {}
    explicit PermCycles(const std::vector<int>& perm) : n((int)perm.size()), p(perm), cycle_id(n, -1), cycle_pos(n, -1) {
        start_.assign(1, 0), list_.reserve(n);
        for (int s = 0; s < n; ++s) {
            if (cycle_id[s] != -1) continue;
            const int c = (int)start_.size() - 1;
            int k = 0;
            for (int x = s; cycle_id[x] == -1; x = p[x]) cycle_id[x] = c, cycle_pos[x] = k++, list_.push_back(x);
            start_.push_back((int)list_.size());
        }
    }
    int count() const { return (int)start_.size() - 1; }  // サイクルの数（不動点も 1 つと数える）
    IntSpan cycle(int c) const { return {list_.data() + start_[c], list_.data() + start_[c + 1]}; }
    int cycle_len(int c) const { return start_[c + 1] - start_[c]; }
    int len_of(int i) const { return cycle_len(cycle_id[i]); }  // i が属するサイクルの長さ
    int kth(int i, long long K) const {  // i に p を K 回あてた値（p を K 回たどった先）。K < 0 なら逆向き
        const int c = cycle_id[i], L = cycle_len(c);
        return list_[start_[c] + (int)(((cycle_pos[i] + K % L) % L + L) % L)];
    }
    std::vector<int> pow(long long K) const {  // p を K 回合成した順列（r[i] = kth(i, K)）。O(n)
        std::vector<int> r(n);
        for (int i = 0; i < n; ++i) r[i] = kth(i, K);
        return r;
    }
    std::vector<int> inverse() const {
        std::vector<int> r(n);
        for (int i = 0; i < n; ++i) r[p[i]] = i;
        return r;
    }
    int min_swaps() const { return n - count(); }  // 任意の 2 つの入れ替えで恒等順列にする最小回数
    int parity() const { return min_swaps() & 1; }  // 偶置換なら 0、奇置換なら 1（転倒数の偶奇と同じ）
    // 配列 A = p（A[i] = p[i]）を、位置の入れ替え (i, j) を繰り返して A[i] = i にする手順（回数は min_swaps()。ABC350 C）
    std::vector<std::pair<int, int>> sort_swaps() const {
        std::vector<int> a(p), where(n);
        for (int i = 0; i < n; ++i) where[a[i]] = i;
        std::vector<std::pair<int, int>> ops;
        for (int i = 0; i < n; ++i)
            if (a[i] != i) {
                const int j = where[i];  // 値 i がある位置
                ops.emplace_back(i, j);
                where[a[i]] = j, where[i] = i;
                std::swap(a[i], a[j]);
            }
        return ops;
    }
    // 周期（何回あてると元に戻るか = サイクルの長さの最小公倍数）。limit を超えたら -1
    long long period(long long limit = (long long)4e18) const {
        long long l = 1;
        for (int c = 0; c < count(); ++c) {
            const long long L = cycle_len(c);
            long long a = l, b = L;
            while (b) a %= b, std::swap(a, b);
            const long long m = L / a;
            if (l > limit / m) return -1;
            l *= m;
        }
        return l;
    }
    // 周期を mod で（大きくなりすぎるとき）。各素因数の最大の指数から計算する
    long long period_mod(long long mod) const {
        std::vector<int> best(n + 1, 0);
        for (int c = 0; c < count(); ++c) {
            int L = cycle_len(c);
            for (int q = 2; q * q <= L; ++q) {
                int e = 0;
                while (L % q == 0) L /= q, ++e;
                best[q] = std::max(best[q], e);
            }
            if (L > 1) best[L] = std::max(best[L], 1);
        }
        long long r = 1 % mod;
        for (int q = 2; q <= n; ++q)
            for (int e = 0; e < best[q]; ++e) r = r * (q % mod) % mod;
        return r;
    }
};
// 合成 r[i] = a[b[i]]（先に b、次に a をあてる）
inline std::vector<int> perm_compose(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> r(b.size());
    for (size_t i = 0; i < b.size(); ++i) r[i] = a[b[i]];
    return r;
}
// P_1 .. P_N を読んで 0-indexed の順列にする（入力が 0-indexed なら base = 0）
inline std::vector<int> read_perm(int n, int base = 1) {
    std::vector<int> a(n);
    for (int& x : a) std::cin >> x, x -= base;
    return a;
}
#endif  // GG_PERM_CYCLES
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "perm_cycles";
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

namespace pc_test {

static long long gcd_ll(long long a, long long b) {
    while (b) a %= b, std::swap(a, b);
    return a;
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    // 最小交換回数の愚直: 恒等順列から転置で BFS（n ≤ 6）
    std::vector<std::vector<std::pair<std::vector<int>, int>>> dist_tab(7);
    for (int n = 0; n <= 6; ++n) {
        std::vector<int> id(n);
        std::iota(id.begin(), id.end(), 0);
        std::vector<std::vector<int>> q{id};
        std::vector<std::pair<std::vector<int>, int>> seen{{id, 0}};
        for (size_t h = 0; h < q.size(); ++h) {
            int d = -1;
            for (auto& pr : seen)
                if (pr.first == q[h]) d = pr.second;
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j) {
                    std::vector<int> r = q[h];
                    std::swap(r[i], r[j]);
                    bool found = false;
                    for (auto& pr : seen) found |= (pr.first == r);
                    if (!found) seen.push_back({r, d + 1}), q.push_back(r);
                }
        }
        dist_tab[n] = seen;
    }
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 0, 12);
        std::vector<int> p(n);
        std::iota(p.begin(), p.end(), 0);
        std::shuffle(p.begin(), p.end(), rng);
        if (it % 7 == 0) std::iota(p.begin(), p.end(), 0);  // 恒等順列
        const PermCycles C(p);
        CHECK(C.n == n && C.p == p);
        // サイクル
        int total = 0, prev_first = -1;
        for (int c = 0; c < C.count(); ++c) {
            const IntSpan cy = C.cycle(c);
            CHECK(cy.size() == C.cycle_len(c) && cy.size() >= 1);
            CHECK(cy[0] > prev_first);  // 小さい要素を含むサイクルから番号がつく
            prev_first = cy[0];
            for (int k = 0; k < cy.size(); ++k) {
                CHECK(p[cy[k]] == cy[(k + 1) % cy.size()] && C.cycle_id[cy[k]] == c && C.cycle_pos[cy[k]] == k);
                CHECK(cy[k] >= cy[0] && C.len_of(cy[k]) == cy.size());
            }
            total += cy.size();
        }
        CHECK(total == n);
        // inverse / kth / pow
        std::vector<int> inv(n);
        for (int i = 0; i < n; ++i) inv[p[i]] = i;
        CHECK(C.inverse() == inv);
        for (int t = 0; t < 6; ++t) {
            const long long K = (t < 4) ? rnd(rng, -3 * n - 3, 3 * n + 3) : (long long)(rng() % 2000000000000000000ULL) - 1000000000000000000LL;
            std::vector<int> r(n);
            for (int i = 0; i < n; ++i) {
                long long Kr = K;
                const int L = C.len_of(i);
                if (K > 3 * L || K < -3 * L) Kr = K % L;
                int x = i;
                for (long long k = 0; k < Kr; ++k) x = p[x];
                for (long long k = 0; k > Kr; --k) x = inv[x];
                r[i] = x;
                CHECK(C.kth(i, K) == x);
            }
            CHECK(C.pow(K) == r);
        }
        CHECK(C.pow(0) == C.pow((long long)C.period()) && C.pow(-1) == inv && C.pow(1) == p);
        // perm_compose
        CHECK(perm_compose(p, p) == C.pow(2) && perm_compose(p, inv) == C.pow(0) && perm_compose(inv, p) == C.pow(0));
        // min_swaps / sort_swaps / parity
        if (n <= 6) {
            int d = -1;
            for (auto& pr : dist_tab[n])
                if (pr.first == p) d = pr.second;
            CHECK(C.min_swaps() == d);
        }
        {
            std::vector<int> a = p;
            const auto ops = C.sort_swaps();
            CHECK((int)ops.size() == C.min_swaps());
            for (const auto& op : ops) {
                CHECK(0 <= op.first && op.first < n && 0 <= op.second && op.second < n && op.first != op.second);
                std::swap(a[op.first], a[op.second]);
            }
            for (int i = 0; i < n; ++i) CHECK(a[i] == i);
        }
        int invs = 0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) invs += (p[i] > p[j]);
        CHECK(C.parity() == (invs & 1));
        // period: 実際にあてて元に戻るまでの回数
        {
            std::vector<int> a = p, id(n);
            std::iota(id.begin(), id.end(), 0);
            long long per = 1;
            while (a != id) a = perm_compose(p, a), ++per;
            CHECK(C.period() == per);
            CHECK(C.period(per) == per && (per == 1 || C.period(per - 1) == -1));
            for (long long mod : {1LL, 2LL, 7LL, 998244353LL, (long long)rnd(rng, 1, 100)}) CHECK(C.period_mod(mod) == per % mod);
        }
    }
    // 大きい周期: 長さが相異なる素数（と 4, 9 など）のサイクルを並べる → 周期は積
    for (int it = 0; it < 200 * TEST_SCALE; ++it, ++cases) {
        std::vector<int> lens;
        for (int q = 2; q <= 300; ++q) {
            bool prime = true;
            for (int d = 2; d * d <= q; ++d) prime &= (q % d != 0);
            if (!prime || rng() % 3 == 0) continue;
            int L = q;
            if (q < 20 && rng() % 2 == 0) L *= q;  // 素数のべき
            lens.push_back(L);
            if (rng() % 4 == 0) lens.push_back(q);  // 同じ素因数の重複（最大の指数だけ効く）
        }
        std::shuffle(lens.begin(), lens.end(), rng);
        int n = 0;
        for (int L : lens) n += L;
        std::vector<int> pos(n), p(n);
        std::iota(pos.begin(), pos.end(), 0);
        std::shuffle(pos.begin(), pos.end(), rng);
        int at = 0;
        long long exact = 1;
        bool over = false;
        const long long LIM = (long long)4e18;
        for (int L : lens) {
            for (int k = 0; k < L; ++k) p[pos[at + k]] = pos[at + (k + 1) % L];
            at += L;
            const long long m = L / gcd_ll(exact, L);
            if (!over && exact > LIM / m) over = true;
            if (!over) exact *= m;
        }
        const PermCycles C(p);
        CHECK(C.count() == (int)lens.size());
        CHECK(C.period() == (over ? -1 : exact));
        // mod での期待値: 各素数の最大指数
        for (long long mod : {998244353LL, 1000000007LL, 2LL, 1LL}) {
            long long ex = 1 % mod;
            for (int q = 2; q <= 300; ++q) {
                bool prime = true;
                for (int d = 2; d * d <= q; ++d) prime &= (q % d != 0);
                if (!prime) continue;
                long long best = 1;
                for (int L : lens) {
                    long long pw = 1;
                    int x = L;
                    while (x % q == 0) x /= q, pw *= q;
                    best = std::max(best, pw);
                }
                ex = ex * (best % mod) % mod;
            }
            CHECK(C.period_mod(mod) == ex);
        }
    }
    {
        const PermCycles E(std::vector<int>{});
        CHECK(E.count() == 0 && E.min_swaps() == 0 && E.period() == 1 && E.period_mod(5) == 1 && E.pow(7).empty() && E.sort_swaps().empty());
        CinFromString cf("3 1 2\n");
        CHECK(read_perm(3) == (std::vector<int>{2, 0, 1}));
    }
    std::cout << "[perm_cycles] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[perm_cycles] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> p(n);
        std::iota(p.begin(), p.end(), 0);
        std::shuffle(p.begin(), p.end(), rng);
        auto t0 = std::chrono::steady_clock::now();
        const PermCycles C(p);
        const long long tb = ms_since(t0);
        t0 = std::chrono::steady_clock::now();
        const std::vector<int> r = C.pow(1000000000000000000LL);
        const auto ops = C.sort_swaps();
        const long long pm = C.period_mod(998244353);
        long long chk = 0;
        for (int q = 0; q < n; ++q) chk += C.kth((int)(rng() % n), (long long)(rng() % 1000000000000000000ULL));
        std::cout << "[perm_cycles] benchmark: N = " << n << "  build " << tb << " ms, pow + sort_swaps + period_mod + kth " << n << " 回 " << ms_since(t0)
                  << " ms  (check " << (r[0] + (long long)ops.size() + pm + chk) % 1000 << ")\n";
    }
}

}  // namespace pc_test


namespace pc_demo {
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
    banner("順列のサイクル分解 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"PermCycles C(p);", "", "p は 0 ~ N-1 の並べ替え（i → p[i]）"},
        {"read_perm(N, base = 1)", "vector<int>", "P_1 .. P_N を読んで 0-indexed に"},
        {"C.count() / C.cycle(c) / C.cycle_len(c)", "", "サイクルの数 / 要素 / 長さ"},
        {"C.cycle_id[i] / C.cycle_pos[i] / C.len_of(i)", "int", "i のサイクル番号 / 位置 / 長さ"},
        {"C.kth(i, K) / C.pow(K)", "int / vector<int>", "K 回あてた値 / K 乗の順列（K < 0 も可）"},
        {"C.inverse()", "vector<int>", "逆順列"},
        {"C.min_swaps() / C.sort_swaps()", "int / vector<pair>", "整列に必要な最小交換回数 / その手順"},
        {"C.parity()", "int", "偶置換 0 / 奇置換 1"},
        {"C.period(limit) / C.period_mod(mod)", "long long", "何回で元に戻るか（超えたら -1）/ mod で"},
        {"perm_compose(a, b)", "vector<int>", "r[i] = a[b[i]]"},
    };
    cheat_sheet("API 一覧", rows, 46, 22);

    head("デモで使う順列: p = {2, 0, 1, 4, 3, 5}（0 → 2 → 1 → 0、3 ↔ 4、5 はそのまま）");
    RUN(std::vector<int> p = {2, 0, 1, 4, 3, 5});
    RUN(PermCycles C(p));

    api("【1】 サイクル", "int / IntSpan / vector<int>", "i → p[i] の矢印でできる輪。番号は小さい要素を含むサイクルから。不動点も長さ 1 のサイクル", "O(N)");
    SHOW(C.count(), "{0, 2, 1}, {3, 4}, {5}");
    SHOW(C.cycle(0), "0 → 2 → 1 の順");
    SHOW(C.cycle(1), "");
    SHOW(C.cycle_id, "各要素のサイクル番号");
    SHOW(C.cycle_pos, "サイクル内の位置");
    SHOW(C.len_of(1), "1 が入っているサイクルの長さ");

    api("【2】 C.kth(i, K) / C.pow(K)", "int / vector<int>", "p を K 回あてた先。K は 10^18 でも負でもよい", "O(1) / O(N)");
    SHOW(C.kth(0, 1), "0 → 2");
    SHOW(C.kth(0, 2), "0 → 2 → 1");
    SHOW(C.kth(0, -1), "逆向きに 1 回（p[x] = 0 となる x）");
    SHOW(C.pow(2), "p を 2 回あてた順列");
    SHOW(C.pow(1000000000000000000LL), "10^18 回でも一瞬");
    SHOW(C.inverse(), "逆順列");
    code("// 「A_i ← A_{P_i} を K 回」の結果は A[C.kth(i, K)]。P が順列でないときは 06_functional_graph の F.kth（ABC367 E）");

    api("【3】 C.min_swaps() / C.sort_swaps()", "int / vector<pair<int,int>>", "好きな 2 つを入れ替えて昇順にする最小回数（N - サイクル数）と、その手順", "O(N)");
    SHOW(C.min_swaps(), "6 - 3");
    SHOW(C.sort_swaps(), "位置の組 (i, j) を順に入れ替えると A[i] = i になる");
    note("ABC350 C は 1-indexed で出力するので、各組に +1 して出す。");
    note("「隣どうしの入れ替え」の最小回数はサイクルではなく転倒数（BIT で数える）。偶奇だけなら C.parity()。");
    SHOW(C.parity(), "0 → 2 → 1 は偶、3 ↔ 4 は奇 → 全体で奇");

    api("【4】 C.period() / C.period_mod(mod)", "long long", "何回あてると元に戻るか（サイクルの長さの最小公倍数）", "O(N log N) / O(N)");
    SHOW(C.period(), "lcm(3, 2, 1)");
    SHOW(C.period(5), "5 を超えるので -1");
    SHOW(C.period_mod(4), "6 mod 4");
    note("N が大きいと周期は long long に収まらない（period は -1）。答えを mod で出すときは period_mod。");

    api("【5】 perm_compose(a, b)", "vector<int>", "r[i] = a[b[i]]（先に b、次に a をあてる）", "O(N)");
    SHOW(perm_compose(p, p), "C.pow(2) と同じ");
    SHOW(perm_compose(p, C.inverse()), "恒等順列");

    api("【6】 read_perm(N, base = 1)", "vector<int>", "P_1 .. P_N（1-indexed）を読んで 0-indexed にする", "O(N)");
    {
        CinFromString cf("3 1 2\n");
        code("// 入力: 3 1 2");
        RUN(auto q = read_perm(3));
        SHOW(q, "1 を引いた");
    }
}

}  // namespace pc_demo

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
    pc_test::random_test();
    pc_test::benchmark();
    pc_demo::demo();
    return 0;
}
