// =============================================================================
//  重み付き Union-Find（ポテンシャル付き Union-Find）            --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  各頂点が「値」を持ち、「値[y] - 値[x] = w」という 2 頂点の差の情報を登録・問い合わせできる。
//
//  ★ 向きの約束（ここだけは覚える）
//    unite(x, y, w) :  値[y] - 値[x] = w  を登録する（「y は x より w 大きい」）
//    diff(x, y)     :  値[y] - 値[x]      を返す
//    weight(x)      :  値[x] - 値[根]     （グループの根の値を 0 としたときの x の値）
//
//  頂点番号   : 0-indexed（0 ~ n-1）
//  重みの型 T : 既定は long long。T() が 0 で、a + b / a - b / a == b が使える型なら何でもよい
//               → int, long long, modint, 下の XorWeight（偶奇・2 色塗り分け・二部グラフ判定）、
//                 自作の 2 次元ベクトル型（相対位置）など。double は == が不安定なので非推奨
//  分類       : 03 Union-Find・最小全域木
//  依存       : 00_graph_base（同梱済み。weighted_union_find_from(G) のため）
//  由来       : 9/23 作成の weighted_union_find.cpp をそのまま採用し、第二世代の形（ファイル名・#ifndef ガード・
//               グラフ基盤の同梱）に合わせて weighted_union_find_from(G) を追加した
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る（<cassert> を使う）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  高速化のポイント
//    - union by size + 経路圧縮。union by size により木の高さは log2(n) 以下なので、
//      経路を固定長のローカル配列に積んで非再帰で圧縮する（vector 確保も再帰もなし）
//    - 親が根のとき（大半のケース）は即 return
//    - par_ 1 本で「親の番号 / 根ならサイズの負値」を兼用
//    - 同じグループの要素を循環リスト nxt_ でつなぎ、members(x) を O(そのグループのサイズ) に
//    - group_count() は O(1)、groups() は各グループのサイズ分を先に reserve
//
//  機能一覧                                                                     計算量
//    WeightedUnionFind<T> wuf(n)  頂点数 n で作る（<T> 省略で long long）       O(n)
//    wuf.find(x)                  x の根（グループの代表元）                   ならし O(α(n))
//    wuf.unite(x, y, w)           値[y] - 値[x] = w を登録して併合。新しく併合したら true
//                                 既に同じグループなら false（矛盾していたら consistent() が false に）
//    wuf.same(x, y)               同じグループなら true                        ならし O(α(n))
//    wuf.size(x)                  x のグループの要素数                         ならし O(α(n))
//    wuf.weight(x)                値[x] - 値[根]                               ならし O(α(n))
//    wuf.diff(x, y)               値[y] - 値[x]（x, y は同じグループであること） ならし O(α(n))
//    wuf.group_count()            グループの個数                               O(1)
//    wuf.groups()                 全グループ（最小要素の昇順・各グループ内も昇順）  O(n α(n))
//    wuf.members(x)               x と同じグループの要素一覧（順不同）         O(そのグループのサイズ)
//    wuf.consistent(x, y, w)      値[y] - 値[x] = w を足しても矛盾しないか（登録はしない）
//    wuf.unite_check(x, y, w)     矛盾しなければ登録して true、矛盾するなら何もせず false
//    ---- 追加機能 ----
//    wuf.consistent()             これまで unite に渡した情報が全て矛盾なしなら true   O(1)
//    wuf.group_ids()              各頂点のグループ番号 0 ~ k-1（groups() の添字と一致）  O(n α(n))
//    wuf.roots()                  各グループの根の一覧（昇順）                 O(n)
//    wuf.add()                    頂点を 1 個追加し、その番号を返す            ならし O(1)
//    wuf.init(n)                  頂点数 n で作り直す（全リセット）            O(n)
//    wuf.n()                      頂点数                                       O(1)
//    XorWeight                    重みを XOR で扱うための型（偶奇・二部グラフ判定用、任意で使う）
//    ---- 第二世代で追加 ----
//    weighted_union_find_from(G)  G の全辺 u → v（重み w）を「値[v] - 値[u] = w」として登録したもの    O(N + M α)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    wuf_test::random_test()  … 愚直解と全機能を突き合わせるランダムテスト（9/23 版のまま）
//    wuf_test::benchmark()    … 速度計測（9/23 版のまま）
//    wuf_demo::demo()         … 全機能の使い方デモ（9/23 版のまま）
//    wuff_test / wuff_demo    … weighted_union_find_from のテストとデモ（第二世代で追加）
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
#include <cassert>

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
#ifndef GG_WEIGHTED_UNION_FIND
#define GG_WEIGHTED_UNION_FIND
// ---- 重み付き Union-Find（9/23 版をそのまま採用。グラフ・グリッド第二世代では weighted_union_find_from(G) を追加） ----
// 必要なヘッダ: <vector>, <utility>, <cassert>
template <class T = long long>
struct WeightedUnionFind {
  public:
    WeightedUnionFind() : WeightedUnionFind(0) {}
    explicit WeightedUnionFind(int n) { init(n); }

    // 頂点数 n で作り直す（全頂点がバラバラ・矛盾フラグもリセット）  O(n)
    void init(int n) {
        n_ = n;
        cnt_ = n;
        ok_ = true;
        par_.assign(n, -1);
        w_.assign(n, T());
        nxt_.resize(n);
        for (int i = 0; i < n; ++i) nxt_[i] = i;
    }

    // 頂点を 1 個追加し、その番号（= 追加前の n()）を返す  ならし O(1)
    int add() {
        par_.push_back(-1);
        w_.push_back(T());
        nxt_.push_back(n_);
        ++cnt_;
        return n_++;
    }

    // 頂点数
    int n() const { return n_; }

    // x の根（代表元）。同じグループなら同じ値を返す  ならし O(α(n))
    // 呼んだ後は par_[x] が根（または x 自身が根）、w_[x] = 値[x] - 値[根] になっている
    int find(int x) {
        const int p = par_[x];
        if (p < 0) return x;          // x が根
        if (par_[p] < 0) return p;    // 親が根（圧縮済み）: 何もしなくてよい
        int path[64];                 // union by size なので木の高さ <= log2(n) < 64
        int len = 0;
        int r = x;
        while (par_[r] >= 0) {
            path[len++] = r;
            r = par_[r];
        }
        // path[len-1] は根の直接の子（重みは既に根からの差）。根に近い側から順に圧縮する
        for (int i = len - 2; i >= 0; --i) {
            const int v = path[i];
            w_[v] = w_[v] + w_[path[i + 1]];
            par_[v] = r;
        }
        return r;
    }

    // 「値[y] - 値[x] = w」を登録して併合する  ならし O(α(n))
    //   新しく併合したら true。既に同じグループなら false（何もしない）。
    //   既に同じグループで登録済みの情報と矛盾していた場合も false で、consistent() が false になる
    bool unite(int x, int y, T w) {
        const int rx = find(x), ry = find(y);
        if (rx == ry) {
            if (!(w_[y] - w_[x] == w)) ok_ = false;
            return false;
        }
        link(rx, ry, w + w_[x] - w_[y]);
        return true;
    }

    // 「値[y] - 値[x] = w」が矛盾しなければ登録して true、矛盾するなら何もせず false  ならし O(α(n))
    //   （既に同じグループで矛盾なしの場合も true。consistent() のフラグは変えない）
    bool unite_check(int x, int y, T w) {
        const int rx = find(x), ry = find(y);
        if (rx == ry) return w_[y] - w_[x] == w;
        link(rx, ry, w + w_[x] - w_[y]);
        return true;
    }

    // 「値[y] - 値[x] = w」を追加しても矛盾しないか（登録はしない）。別グループなら常に true
    bool consistent(int x, int y, T w) {
        const int rx = find(x), ry = find(y);
        return rx != ry || w_[y] - w_[x] == w;
    }

    // これまで unite に渡した情報が全て矛盾なしなら true  O(1)
    bool consistent() const { return ok_; }

    // x と y が同じグループなら true  ならし O(α(n))
    bool same(int x, int y) { return find(x) == find(y); }

    // x が属するグループの要素数  ならし O(α(n))
    int size(int x) { return -par_[find(x)]; }

    // 値[x] - 値[根]（根の値を 0 としたときの x の値）  ならし O(α(n))
    T weight(int x) {
        find(x);
        return w_[x];
    }

    // 値[y] - 値[x]。x と y は同じグループであること（別グループなら assert で停止）  ならし O(α(n))
    T diff(int x, int y) {
        const int rx = find(x), ry = find(y);
        assert(rx == ry);
        (void)rx;
        (void)ry;
        return w_[y] - w_[x];
    }

    // グループの個数  O(1)
    int group_count() const { return cnt_; }

    // x と同じグループの要素一覧（x から始まる順不同。昇順が必要なら sort）  O(サイズ)
    std::vector<int> members(int x) {
        std::vector<int> res;
        res.reserve(size(x));
        int v = x;
        do {
            res.push_back(v);
            v = nxt_[v];
        } while (v != x);
        return res;
    }

    // 全グループ。グループは最小要素の昇順、各グループ内も昇順  O(n α(n))
    std::vector<std::vector<int>> groups() {
        std::vector<int> id(n_, -1);
        std::vector<std::vector<int>> res;
        res.reserve(cnt_);
        for (int i = 0; i < n_; ++i) {
            const int r = find(i);
            if (id[r] < 0) {
                id[r] = (int)res.size();
                res.emplace_back();
                res.back().reserve(-par_[r]);
            }
            res[id[r]].push_back(i);
        }
        return res;
    }

    // 各頂点のグループ番号 0 ~ group_count()-1（groups() の添字と一致）  O(n α(n))
    std::vector<int> group_ids() {
        std::vector<int> id(n_, -1), res(n_);
        int k = 0;
        for (int i = 0; i < n_; ++i) {
            const int r = find(i);
            if (id[r] < 0) id[r] = k++;
            res[i] = id[r];
        }
        return res;
    }

    // 各グループの根の一覧（昇順）  O(n)
    std::vector<int> roots() const {
        std::vector<int> res;
        res.reserve(cnt_);
        for (int i = 0; i < n_; ++i)
            if (par_[i] < 0) res.push_back(i);
        return res;
    }

  private:
    int n_ = 0;              // 頂点数
    int cnt_ = 0;            // グループ数
    bool ok_ = true;         // これまでの unite に矛盾がないか
    std::vector<int> par_;   // 根なら -(サイズ)、それ以外は親の番号
    std::vector<T> w_;       // 値[x] - 値[親]（根では T()）
    std::vector<int> nxt_;   // 同じグループの要素をつなぐ循環リスト

    // rx, ry は異なる根で d = 値[ry] - 値[rx]。サイズの大きい方（同じなら rx）を根にする
    void link(int rx, int ry, T d) {
        if (par_[rx] > par_[ry]) {  // par_ はサイズの負値なので「>」で rx の方が小さい
            std::swap(rx, ry);
            d = T() - d;
        }
        par_[rx] += par_[ry];
        par_[ry] = rx;
        w_[ry] = d;
        std::swap(nxt_[rx], nxt_[ry]);  // 2 つの循環リストをつなげて 1 つにする
        --cnt_;
    }
};

// ---- 任意: XOR を「足し算」とみなす重み（偶奇・2 色塗り分け・二部グラフ判定用） ----
//   WeightedUnionFind<XorWeight> c(n);
//   c.unite(a, b, 1) :「a と b は違う色」   c.unite(a, b, 0) :「a と b は同じ色」
//   c.diff(a, b).v   : 0 なら同じ色、1 なら違う色（.v で中身の整数を取り出す）
struct XorWeight {
    unsigned long long v;
    XorWeight(unsigned long long v_ = 0) : v(v_) {}
    XorWeight operator+(const XorWeight& o) const { return XorWeight(v ^ o.v); }
    XorWeight operator-(const XorWeight& o) const { return XorWeight(v ^ o.v); }
    bool operator==(const XorWeight& o) const { return v == o.v; }
    bool operator!=(const XorWeight& o) const { return v != o.v; }
};
// グラフ G の全部の辺 u → v（重み w）を「値[v] - 値[u] = w」として登録した WeightedUnionFind<Cost>。
// 矛盾があれば consistent() が false。無向辺も「値[to] - 値[from] = w」（入力の向き）として扱う
inline WeightedUnionFind<Cost> weighted_union_find_from(const Graph& G) {
    WeightedUnionFind<Cost> w(G.size());
    for (const RawEdge& e : G.es) w.unite(e.from, e.to, e.w);
    return w;
}
#endif  // GG_WEIGHTED_UNION_FIND
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト: 愚直な実装（各頂点の値を直接持つ方式）と全機能を突き合わせる
// =============================================================================
namespace wuf_test {

// NDEBUG（VS の Release ビルド）でも無効にならないチェック
#define WUF_CHECK(cond)                                                                          \
    do {                                                                                         \
        if (!(cond)) {                                                                           \
            std::cerr << "[WeightedUnionFind] TEST FAILED  line " << __LINE__ << ": " #cond "\n"; \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)

// 愚直解: 各頂点の「値」val[v] とグループ番号 lab[v] をそのまま持つ。
//   併合時は小さい方のグループ全員の値をずらし、ラベルを書き換える。
template <class T>
struct Naive {
    std::vector<int> lab;
    std::vector<std::vector<int>> mem;  // ラベルごとの要素
    std::vector<T> val;
    int cnt = 0;
    bool ok = true;

    void init(int n) {
        lab.resize(n);
        mem.assign(n, std::vector<int>());
        for (int i = 0; i < n; ++i) lab[i] = i, mem[i].push_back(i);
        val.assign(n, T());
        cnt = n;
        ok = true;
    }
    int add() {
        const int v = (int)lab.size();
        lab.push_back(v);
        mem.push_back(std::vector<int>(1, v));
        val.push_back(T());
        ++cnt;
        return v;
    }
    int n() const { return (int)lab.size(); }
    bool same(int x, int y) const { return lab[x] == lab[y]; }
    int size(int x) const { return (int)mem[lab[x]].size(); }
    T diff(int x, int y) const { return val[y] - val[x]; }
    bool consistent(int x, int y, T w) const { return !same(x, y) || diff(x, y) == w; }
    void merge(int x, int y, T w) {  // 別グループの x, y に 値[y] - 値[x] = w を課す
        const int a = lab[x], b = lab[y];
        const T delta = val[x] + w - val[y];  // グループ b 全体に足せば条件を満たす
        if (mem[a].size() >= mem[b].size()) {
            for (int v : mem[b]) val[v] = val[v] + delta, lab[v] = a;
            mem[a].insert(mem[a].end(), mem[b].begin(), mem[b].end());
            mem[b].clear();
        } else {  // 逆にグループ a 全体を -delta ずらしても同じ
            const T nd = T() - delta;
            for (int v : mem[a]) val[v] = val[v] + nd, lab[v] = b;
            mem[b].insert(mem[b].end(), mem[a].begin(), mem[a].end());
            mem[a].clear();
        }
        --cnt;
    }
    bool unite(int x, int y, T w) {
        if (same(x, y)) {
            if (!(diff(x, y) == w)) ok = false;
            return false;
        }
        merge(x, y, w);
        return true;
    }
    bool unite_check(int x, int y, T w) {
        if (same(x, y)) return diff(x, y) == w;
        merge(x, y, w);
        return true;
    }
    std::vector<int> group_of(int x) const {  // 昇順
        std::vector<int> res = mem[lab[x]];
        std::sort(res.begin(), res.end());
        return res;
    }
    std::vector<std::vector<int>> groups() const {  // 最小要素の昇順
        std::vector<int> id(n(), -1);
        std::vector<std::vector<int>> res;
        for (int i = 0; i < n(); ++i) {
            if (id[lab[i]] < 0) {
                id[lab[i]] = (int)res.size();
                res.emplace_back();
            }
            res[id[lab[i]]].push_back(i);
        }
        return res;
    }
    std::vector<int> ids() const {
        std::vector<int> id(n(), -1), res(n());
        int k = 0;
        for (int i = 0; i < n(); ++i) {
            if (id[lab[i]] < 0) id[lab[i]] = k++;
            res[i] = id[lab[i]];
        }
        return res;
    }
};

// テスト量の倍率。デバッグビルド（最適化なし・範囲チェックあり）は非常に遅いので自動で 1/8 にする。
// もっと念入りに検査したいときはこの値を大きくする
#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
static const int TEST_SCALE = 1;
#else
static const int TEST_SCALE = 8;
#endif

static int rnd(std::mt19937_64& rng, int lo, int hi) {  // [lo, hi]
    return lo + (int)(rng() % (unsigned long long)(hi - lo + 1));
}

// 全機能の全頂点チェック。full = true なら全ペアの diff と全頂点の members も確認（O(n^2)）
template <class T>
static void check_all(WeightedUnionFind<T>& uf, const Naive<T>& nv, bool full, std::mt19937_64& rng) {
    const int n = nv.n();
    WUF_CHECK(uf.n() == n);
    WUF_CHECK(uf.group_count() == nv.cnt);
    WUF_CHECK(uf.consistent() == nv.ok);

    // find / size / weight
    std::vector<int> root_of_label(n, -1);
    for (int v = 0; v < n; ++v) {
        const int r = uf.find(v);
        WUF_CHECK(0 <= r && r < n);
        WUF_CHECK(uf.find(r) == r);
        WUF_CHECK(nv.lab[r] == nv.lab[v]);
        int& rl = root_of_label[nv.lab[v]];
        if (rl < 0) rl = r;
        WUF_CHECK(rl == r);
        WUF_CHECK(uf.size(v) == nv.size(v));
        WUF_CHECK(uf.weight(v) == nv.val[v] - nv.val[r]);
        WUF_CHECK(uf.weight(r) == T());
    }

    // groups / group_ids は並び順まで完全一致
    const std::vector<std::vector<int>> g = uf.groups();
    WUF_CHECK(g == nv.groups());
    const std::vector<int> ids = uf.group_ids();
    WUF_CHECK(ids == nv.ids());
    for (int i = 0; i < (int)g.size(); ++i)
        for (int v : g[i]) WUF_CHECK(ids[v] == i);

    // roots
    const std::vector<int> rs = uf.roots();
    WUF_CHECK((int)rs.size() == uf.group_count());
    WUF_CHECK(std::is_sorted(rs.begin(), rs.end()));
    std::vector<char> seen(n, 0);
    for (int r : rs) {
        WUF_CHECK(uf.find(r) == r);
        WUF_CHECK(!seen[nv.lab[r]]);
        seen[nv.lab[r]] = 1;
    }

    if (full) {
        for (int x = 0; x < n; ++x) {
            std::vector<int> m = uf.members(x);
            WUF_CHECK(!m.empty() && m[0] == x);
            std::sort(m.begin(), m.end());
            WUF_CHECK(m == nv.group_of(x));
            for (int y = 0; y < n; ++y) {
                WUF_CHECK(uf.same(x, y) == nv.same(x, y));
                if (nv.same(x, y)) {
                    WUF_CHECK(uf.diff(x, y) == nv.diff(x, y));
                    WUF_CHECK(uf.consistent(x, y, nv.diff(x, y)));
                }
            }
        }
    } else if (n > 0) {
        for (int k = 0; k < 2 * n; ++k) {
            const int x = rnd(rng, 0, n - 1);
            const std::vector<int>& m = nv.mem[nv.lab[x]];
            const int y = m[rnd(rng, 0, (int)m.size() - 1)];
            WUF_CHECK(uf.diff(x, y) == nv.diff(x, y));
        }
    }
}

// ランダム操作列を 1 本流す。gen(rng) で重みの乱数を作る
template <class T, class Gen>
static long long run_case(std::mt19937_64& rng, int n, int q, Gen gen) {
    WeightedUnionFind<T> uf(n);
    Naive<T> nv;
    nv.init(n);
    const bool small = (n <= 60);
    const int check_every = std::max(1, n / 2);
    // 同じグループなら半分の確率で「矛盾しない正しい差」、それ以外はランダム
    auto pick_w = [&](int x, int y) -> T {
        if (nv.same(x, y) && rnd(rng, 0, 1) == 0) return nv.diff(x, y);
        return gen(rng);
    };
    for (int t = 0; t < q; ++t) {
        const int cur = nv.n();
        const int type = rnd(rng, 0, 99);
        if (cur == 0 || type < 2) {  // add
            WUF_CHECK(uf.add() == nv.add());
            continue;
        }
        int x = rnd(rng, 0, cur - 1), y = rnd(rng, 0, cur - 1);
        if (rnd(rng, 0, 3) == 0) {  // 同じグループのペアも十分に作る
            const std::vector<int>& m = nv.mem[nv.lab[x]];
            y = m[rnd(rng, 0, (int)m.size() - 1)];
        }
        if (type < 30) {  // unite
            const int rx = uf.find(x), ry = uf.find(y);
            const int sx = uf.size(x), sy = uf.size(y);
            const T w = pick_w(x, y);
            const bool expect = nv.unite(x, y, w);
            WUF_CHECK(uf.unite(x, y, w) == expect);
            WUF_CHECK(uf.consistent() == nv.ok);
            if (expect) {
                WUF_CHECK(uf.find(x) == (sx >= sy ? rx : ry));  // union by size（同じなら x 側）
                WUF_CHECK(uf.size(y) == sx + sy);
                WUF_CHECK(uf.diff(x, y) == w);
            }
        } else if (type < 42) {  // unite_check
            const T w = pick_w(x, y);
            const bool before_ok = uf.consistent();
            const bool was_same = nv.same(x, y);
            const bool expect = nv.unite_check(x, y, w);
            WUF_CHECK(uf.unite_check(x, y, w) == expect);
            WUF_CHECK(uf.consistent() == before_ok);  // フラグは変わらない
            WUF_CHECK(uf.same(x, y) == (expect || was_same));
            if (expect) WUF_CHECK(uf.diff(x, y) == w);
        } else if (type < 52) {  // consistent(x, y, w)（登録しない）
            const T w = pick_w(x, y);
            WUF_CHECK(uf.consistent(x, y, w) == nv.consistent(x, y, w));
            WUF_CHECK(uf.same(x, y) == nv.same(x, y));
        } else if (type < 60) {  // weight
            WUF_CHECK(uf.weight(x) == nv.val[x] - nv.val[uf.find(x)]);
        } else if (type < 70) {  // diff
            if (nv.same(x, y)) WUF_CHECK(uf.diff(x, y) == nv.diff(x, y));
        } else if (type < 75) {  // same
            WUF_CHECK(uf.same(x, y) == nv.same(x, y));
        } else if (type < 79) {  // find
            WUF_CHECK((uf.find(x) == uf.find(y)) == nv.same(x, y));
            WUF_CHECK(nv.lab[uf.find(x)] == nv.lab[x]);
        } else if (type < 83) {  // size
            WUF_CHECK(uf.size(x) == nv.size(x));
        } else if (type < 88) {  // members
            std::vector<int> m = uf.members(x);
            WUF_CHECK(m[0] == x);
            std::sort(m.begin(), m.end());
            WUF_CHECK(m == nv.group_of(x));
        } else if (type < 90) {  // group_count
            WUF_CHECK(uf.group_count() == nv.cnt);
        } else if (type < 92) {  // groups
            WUF_CHECK(uf.groups() == nv.groups());
        } else if (type < 94) {  // group_ids
            WUF_CHECK(uf.group_ids() == nv.ids());
        } else if (type < 95) {  // init（作り直し）
            const int m = rnd(rng, 0, std::max(1, n));
            uf.init(m);
            nv.init(m);
        } else {  // n / roots / consistent()
            WUF_CHECK(uf.n() == nv.n());
            WUF_CHECK((int)uf.roots().size() == nv.cnt);
            WUF_CHECK(uf.consistent() == nv.ok);
        }
        if (t % check_every == 0) check_all(uf, nv, small, rng);
    }
    check_all(uf, nv, nv.n() <= 300, rng);
    return q;
}

// 最も木が深くなる併合順（同じサイズ同士を根で併合し続ける = 二項木）を作り、
// 一番深い頂点から問い合わせて、経路圧縮時の重みの足し込みを検証する
template <class T, class Gen>
static void deep_case(std::mt19937_64& rng, int k, Gen gen) {
    const int n = 1 << k;
    WeightedUnionFind<T> uf(n);
    Naive<T> nv;
    nv.init(n);
    for (int s = 1; s < n; s *= 2)
        for (int i = 0; i + s < n; i += 2 * s) {
            const T w = gen(rng);
            WUF_CHECK(nv.unite(i, i + s, w));
            WUF_CHECK(uf.unite(i, i + s, w));
        }
    // 頂点 n-1 の深さは k。最初の問い合わせで k 段の経路圧縮が起きる
    WUF_CHECK(uf.weight(n - 1) == nv.val[n - 1] - nv.val[0]);
    WUF_CHECK(uf.find(n - 1) == 0);
    // 深い頂点を後ろから順に（圧縮されていない経路が多く残る順に）
    for (int v = n - 1; v >= 0; v -= 3) WUF_CHECK(uf.diff(0, v) == nv.val[v] - nv.val[0]);
    check_all(uf, nv, n <= 64, rng);
}

// 形の決まったケース（n = 0、add だけで作る、自己ループ、一直線 など）
static void structured_cases(std::mt19937_64& rng) {
    {
        WeightedUnionFind<> uf;  // 型省略 = long long
        WUF_CHECK(uf.n() == 0 && uf.group_count() == 0 && uf.groups().empty() && uf.roots().empty());
        WUF_CHECK(uf.group_ids().empty() && uf.consistent());
        for (int i = 0; i < 5; ++i) WUF_CHECK(uf.add() == i);
        WUF_CHECK(uf.unite(0, 4, 10) && uf.unite(4, 2, -3) && !uf.unite(0, 2, 7) && uf.consistent());
        WUF_CHECK(uf.diff(0, 2) == 7 && uf.diff(2, 0) == -7 && uf.size(2) == 3);
        WUF_CHECK(!uf.unite(2, 2, 0) && uf.consistent());  // 自己ループ w = 0 は矛盾なし
        WUF_CHECK(!uf.unite(2, 2, 1) && !uf.consistent());  // 自己ループ w != 0 は矛盾
        uf.init(2);
        WUF_CHECK(uf.consistent() && uf.group_count() == 2);
    }
    for (int n : {1, 2, 3, 100, 1000}) {  // 一直線（常に大きい方に 1 頂点をつなぐ）
        WeightedUnionFind<long long> uf(n);
        Naive<long long> nv;
        nv.init(n);
        for (int i = 0; i + 1 < n; ++i) {
            const long long w = rnd(rng, -1000000000, 1000000000);
            WUF_CHECK(uf.unite(i + 1, i, w) == nv.unite(i + 1, i, w));
        }
        check_all(uf, nv, n <= 100, rng);
    }
}

static void random_test() {
    std::mt19937_64 rng(20260923);
    long long ops = 0;
    int cases = 0;
    auto gen_wide = [](std::mt19937_64& r) -> long long { return rnd(r, -1000000000, 1000000000); };
    auto gen_narrow = [](std::mt19937_64& r) -> long long { return rnd(r, -2, 2); };  // 偶然の一致も起こす
    auto gen_int = [](std::mt19937_64& r) -> int { return rnd(r, -3, 3); };
    auto gen_bit = [](std::mt19937_64& r) -> XorWeight { return XorWeight(r() & 1); };
    auto gen_x64 = [](std::mt19937_64& r) -> XorWeight { return XorWeight(r()); };

    // 小さい n を大量に（型・重みの範囲を変えて）
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = rnd(rng, 0, 12), q = rnd(rng, 0, 80);
        ops += run_case<long long>(rng, n, q, gen_wide);
        ops += run_case<long long>(rng, n, q, gen_narrow);
        ops += run_case<int>(rng, n, q, gen_int);
        ops += run_case<XorWeight>(rng, n, q, gen_bit);
        ops += run_case<XorWeight>(rng, n, q, gen_x64);
        cases += 5;
    }
    // 中くらいの n
    for (int it = 0; it < 30 * TEST_SCALE; ++it) {
        const int n = rnd(rng, 13, 150);
        ops += run_case<long long>(rng, n, 5 * n, gen_wide);
        ops += run_case<long long>(rng, n, 5 * n, gen_narrow);
        ops += run_case<XorWeight>(rng, n, 5 * n, gen_bit);
        cases += 3;
    }
    // 大きめの n
    for (int it = 0; it < TEST_SCALE; ++it) {
        const int n = rnd(rng, 1000, 5000);
        ops += run_case<long long>(rng, n, 4 * n, gen_wide);
        ops += run_case<long long>(rng, n, 4 * n, gen_narrow);
        cases += 2;
    }
    // 最も深い木での経路圧縮
    for (int k = 0; k <= 16; ++k) {
        deep_case<long long>(rng, k, gen_wide);
        deep_case<XorWeight>(rng, k, gen_x64);
    }
    structured_cases(rng);
    std::cout << "[WeightedUnionFind] random test: OK  (" << cases << " cases, " << ops
              << " ops + deepest-tree / structured cases)\n";
}

// 速度計測: n 頂点に対して「unite 1 回 : (same なら diff) 3 回」のランダム操作を q 回
static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[WeightedUnionFind] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int ns[] = {200000, 1000000};
    for (int n : ns) {
        if (TEST_SCALE == 1 && n > 200000) break;
        const int q = 4 * n;
        std::vector<int> a(q), b(q);
        std::vector<long long> w(q);
        for (int i = 0; i < q; ++i) a[i] = (int)(rng() % n), b[i] = (int)(rng() % n), w[i] = (long long)(rng() % 1000);
        const auto t0 = std::chrono::steady_clock::now();
        WeightedUnionFind<long long> uf(n);
        long long chk = 0;
        for (int i = 0; i < q; ++i) {
            if ((i & 3) == 0) chk += uf.unite(a[i], b[i], w[i]);
            else if (uf.same(a[i], b[i])) chk += uf.diff(a[i], b[i]);
        }
        const auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "[WeightedUnionFind] benchmark: n = " << n << ", " << q << " ops (unite:same+diff = 1:3)  "
                  << (long long)ms << " ms  (checksum " << chk << ")\n";
    }
}

}  // namespace wuf_test


// =============================================================================
//  使い方デモ: 全 API を 1 つずつ「返り値・内容・計算量・注意点」つきで実行して見せる。
//    出力の見方:   実際に書くコード  ->  返り値    // 何を返しているか（日本語）
//                  [状態] の行 … その時点のグループ分けと各頂点の weight（* が根）
// =============================================================================
namespace wuf_demo {

// 典型5 で使う自作の重みの型: 2 次元の相対位置 (x, y)
//   WeightedUnionFind の T に使うには「T() が 0」「+」「-」「==」があればよい
struct Vec2 {
    long long x = 0, y = 0;
    Vec2() {}
    Vec2(long long x_, long long y_) : x(x_), y(y_) {}
    Vec2 operator+(const Vec2& o) const { return Vec2(x + o.x, y + o.y); }
    Vec2 operator-(const Vec2& o) const { return Vec2(x - o.x, y - o.y); }
    bool operator==(const Vec2& o) const { return x == o.x && y == o.y; }
};

static std::string to_s(bool b) { return b ? "true" : "false"; }
static std::string to_s(int v) { return std::to_string(v); }
static std::string to_s(long long v) { return std::to_string(v); }
static std::string to_s(unsigned long long v) { return std::to_string(v); }
static std::string to_s(const std::string& s) { return s; }
static std::string to_s(const Vec2& p) { return "(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ")"; }
template <class V>
static std::string to_s(const std::vector<V>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
static std::string pad(const std::string& s, size_t w) { return s.size() >= w ? s : s + std::string(w - s.size(), ' '); }

static void banner(const std::string& s) {
    std::cout << "\n=====================================================================\n"
              << "  " << s << "\n"
              << "=====================================================================\n";
}
// API の見出し: 使い方 / 返り値 / 内容 / 計算量
static void api(const std::string& usage, const std::string& ret, const std::string& what, const std::string& cost) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << usage << "\n"
              << "    返り値 : " << ret << "\n"
              << "    内容   : " << what << "\n"
              << "    計算量 : " << cost << "\n";
}
static void head(const std::string& s) {  // 典型パターンなどの見出し
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
static void note(const std::string& s) { std::cout << "    " << s << "\n"; }  // 補足説明
static void code(const std::string& s) { std::cout << "  " << s << "\n"; }    // 実行したコード
// 文の末尾に ; を補う（for / if / while で } で終わる文には付けない）
static void stmt(const std::string& s) {
    const bool block = !s.empty() && s.back() == '}' &&
                       (s.compare(0, 4, "for ") == 0 || s.compare(0, 3, "if ") == 0 || s.compare(0, 6, "while ") == 0);
    code(block ? s : s + ";");
}
// 現在のグループ分けと weight を 1 行で表示（頂点:weight、* が根）
static void state(WeightedUnionFind<long long>& wuf) {
    std::string s;
    for (const auto& g : wuf.groups()) {
        const int r = wuf.find(g[0]);
        s += "{";
        for (size_t i = 0; i < g.size(); ++i)
            s += (i ? ", " : "") + std::to_string(g[i]) + (g[i] == r ? "*" : "") + ":" + std::to_string(wuf.weight(g[i]));
        s += "} ";
    }
    std::cout << "    [状態] " << s << " (グループ数 " << wuf.group_count() << ", consistent " << to_s(wuf.consistent())
              << ")\n";
}
template <class V>
static void show(const char* expr, const V& v, const std::string& comment) {  // コード -> 返り値 // 説明
    std::string line = "  " + pad(expr, 32) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 46) + "  // " + comment;
    std::cout << line << "\n";
}
// SHOW(式, "説明") : 式そのものと返り値と説明を表示
#define SHOW(expr, comment) wuf_demo::show(#expr, (expr), comment)
// RUN(文) : 文を表示してから実行（宣言にも使える）
#define RUN(...)                   \
    wuf_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

// 最初に表示する API 一覧（チートシート）
static void cheat_sheet() {
    struct Row {
        const char* use;
        const char* ret;
        const char* what;
    };
    const Row rows[] = {
        {"WeightedUnionFind<T> wuf(n);", "", "頂点 0 ~ n-1 で作る。<> で省略すると T = long long  O(n)"},
        {"wuf.n()", "int", "頂点数  O(1)"},
        {"wuf.unite(x, y, w)", "bool", "値[y] - 値[x] = w を登録して併合。新しく併合したら true  O(α(n))"},
        {"wuf.diff(x, y)", "T", "値[y] - 値[x]（x, y は同じグループであること）  O(α(n))"},
        {"wuf.weight(x)", "T", "値[x] - 値[根]（根を 0 としたときの x の値）  O(α(n))"},
        {"wuf.same(x, y)", "bool", "同じグループなら true  O(α(n))"},
        {"wuf.find(x)", "int", "x の根（代表元）  O(α(n))"},
        {"wuf.size(x)", "int", "x のグループの要素数  O(α(n))"},
        {"wuf.group_count()", "int", "グループの個数  O(1)"},
        {"wuf.groups()", "vector<vector<int>>", "全グループ（最小要素の昇順・中も昇順）  O(n α(n))"},
        {"wuf.members(x)", "vector<int>", "x と同じグループの要素（順不同）  O(サイズ)"},
        {"wuf.consistent(x, y, w)", "bool", "値[y] - 値[x] = w を足しても矛盾しないか（登録しない）  O(α(n))"},
        {"wuf.unite_check(x, y, w)", "bool", "矛盾しなければ登録して true / 矛盾なら何もせず false  O(α(n))"},
        {"wuf.consistent()", "bool", "これまでの unite に矛盾が無ければ true  O(1)"},
        {"wuf.group_ids()", "vector<int>", "各頂点のグループ番号 0 ~ k-1（groups() の添字）  O(n α(n))"},
        {"wuf.roots()", "vector<int>", "各グループの根（昇順）  O(n)"},
        {"wuf.add()", "int", "頂点を 1 個追加し、その番号を返す  ならし O(1)"},
        {"wuf.init(n)", "void", "頂点数 n で作り直す（矛盾フラグもリセット）  O(n)"},
    };
    std::cout << "\n■ API 一覧（x, y は頂点番号 0 ~ n-1。α(n) はほぼ定数。詳しくは下の【1】~【18】）\n";
    std::cout << "  ★ 向きの約束: unite(x, y, w) は「値[y] - 値[x] = w」（y は x より w 大きい）、diff(x, y) = 値[y] - 値[x]\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, 29) << pad(ret, 24) << r.what << "\n";
    }
}

static void demo() {
    banner("Weighted Union-Find 使い方デモ    見方:  コード  ->  返り値    // 意味");
    cheat_sheet();

    // ------------------------------------------------------------------ 【1】
    api("【1】 WeightedUnionFind<T> wuf(n);   /   WeightedUnionFind<T> wuf;（空で作って後から init(n) や add()）",
        "（コンストラクタ）", "頂点 0 ~ n-1 を作る。最初は全員が 1 人ずつ別々のグループ。T は値の差（重み）の型", "O(n)");
    note("T の例: long long（普通はこれ）/ int / modint / XorWeight（典型4）/ 自作の構造体（典型5）");
    note("WeightedUnionFind<> wuf(n); と <> だけ書くと T = long long");
    RUN(WeightedUnionFind<long long> wuf(6));
    state(wuf);
    RUN(WeightedUnionFind<> empty_wuf);
    SHOW(empty_wuf.n(), "空で作ると頂点数 0（T は long long）");
    note("1-indexed の入力は wuf.unite(a - 1, b - 1, w) のように 0-indexed に直す");

    // ------------------------------------------------------------------ 【2】
    api("【2】 wuf.n()", "int", "頂点数（add() で増える）", "O(1)");
    SHOW(wuf.n(), "頂点数");

    // ------------------------------------------------------------------ 【3】
    api("【3】 wuf.unite(x, y, w)", "bool … 新しく併合したら true / 既に同じグループなら false",
        "「値[y] - 値[x] = w」（y は x より w 大きい）を登録して、x と y のグループを 1 つにまとめる", "ならし O(α(n))");
    note("※ unite(x, y, w) と unite(y, x, -w) は同じ意味");
    note("※ 既に同じグループのときは何も変えずに false。登録済みの情報と矛盾していたら consistent() が false になる");
    SHOW(wuf.unite(0, 1, 5), "値[1] - 値[0] = 5（1 は 0 より 5 大きい）");
    state(wuf);
    SHOW(wuf.unite(2, 1, -3), "値[1] - 値[2] = -3（2 は 1 より 3 大きい）");
    SHOW(wuf.unite(3, 4, -2), "値[4] - 値[3] = -2（4 は 3 より 2 小さい）");
    state(wuf);
    SHOW(wuf.unite(0, 2, 8), "既に同じグループ → false（5 + 3 = 8 なので矛盾はない）");
    SHOW(wuf.consistent(), "矛盾していないので true のまま");
    state(wuf);
    note("[状態] の 頂点:数字 は weight（根 * の値を 0 としたときの値）");

    // ------------------------------------------------------------------ 【4】
    api("【4】 wuf.diff(x, y)", "T … 値[y] - 値[x]", "x から見た y の値の差。直接 unite していない 2 頂点でも、同じグループなら求まる",
        "ならし O(α(n))");
    note("※ x と y が別グループだと差は決まらないので assert で停止する。先に same で確認すること");
    SHOW(wuf.diff(0, 1), "値[1] - 値[0]");
    SHOW(wuf.diff(0, 2), "値[2] - 値[0] = 5 + 3（0 と 2 は直接 unite していない）");
    SHOW(wuf.diff(2, 0), "向きを逆にすると符号が逆");
    SHOW(wuf.diff(1, 1), "自分自身とは 0");
    SHOW(wuf.diff(3, 4), "値[4] - 値[3]");
    SHOW(wuf.same(0, 5), "別グループなので diff(0, 5) は呼んではいけない");
    note("安全な書き方:");
    RUN(if (wuf.same(0, 5)) std::cout << "      " << wuf.diff(0, 5) << "\n"; else std::cout << "      unknown\n");

    // ------------------------------------------------------------------ 【5】
    api("【5】 wuf.weight(x)", "T … 値[x] - 値[根]", "根の値を 0 としたときの x の値", "ならし O(α(n))");
    note("※ 根は unite で変わることがあり、そのとき weight も変わる。2 頂点を比べるときは diff を使う");
    note("※ 全頂点の weight(i) を並べると、登録した情報を全て満たす値の割り当ての 1 つになる（典型2）");
    SHOW(wuf.find(0), "{0,1,2} の根");
    SHOW(wuf.weight(0), "根の weight は常に 0");
    SHOW(wuf.weight(1), "値[1] - 値[0]");
    SHOW(wuf.weight(2), "値[2] - 値[0]");
    SHOW(wuf.weight(4), "{3,4} の根は 3 → 値[4] - 値[3]");
    SHOW(wuf.weight(5), "1 人だけのグループは自分が根なので 0");
    note("根が変わる例（別の WeightedUnionFind で確認）:");
    {
        RUN(WeightedUnionFind<long long> t(4));
        RUN(t.unite(2, 3, 10));
        SHOW(t.weight(3), "根は 2 → 値[3] - 値[2]");
        RUN(t.unite(0, 1, 1));
        RUN(t.unite(1, 2, 1));
        SHOW(t.find(3), "{0,1} と併合した結果、根が 0 に変わった");
        SHOW(t.weight(3), "値[3] - 値[0] に変わる（= 1 + 1 + 10）");
        SHOW(t.diff(2, 3), "diff は根に関係なく一定");
    }

    // ------------------------------------------------------------------ 【6】
    api("【6】 wuf.same(x, y)", "bool", "x と y が同じグループなら true（= 値の差が決まっているか）", "ならし O(α(n))");
    SHOW(wuf.same(0, 2), "0-1-2 はつながっている");
    SHOW(wuf.same(0, 3), "0 と 3 は別のグループ（値の差は不明）");
    SHOW(wuf.same(5, 5), "自分自身とは常に true");

    // ------------------------------------------------------------------ 【7】
    api("【7】 wuf.find(x)", "int … x の根（グループの代表元）の頂点番号",
        "同じグループなら同じ値が返る。weight(x) はこの根を基準にした値", "ならし O(α(n))");
    note("※ どの頂点が根になるかに意味はない。unite すると変わることがある（【5】の例）");
    SHOW(wuf.find(0), "{0,1,2} の根");
    SHOW(wuf.find(2), "同じグループなので find(0) と同じ値");
    SHOW(wuf.find(4), "{3,4} の根");
    SHOW(wuf.find(5) == 5, "find(x) == x なら x は根");

    // ------------------------------------------------------------------ 【8】
    api("【8】 wuf.size(x)", "int", "x が属するグループの要素数", "ならし O(α(n))");
    SHOW(wuf.size(0), "{0,1,2}");
    SHOW(wuf.size(4), "{3,4}");
    SHOW(wuf.size(5), "{5}");

    // ------------------------------------------------------------------ 【9】
    api("【9】 wuf.group_count()", "int", "グループの個数。unite が true になるたびに 1 減る", "O(1)");
    SHOW(wuf.group_count(), "{0,1,2} {3,4} {5} の 3 個");
    note("使用例: 全員の値の差が決まっているか → wuf.group_count() == 1");

    // ------------------------------------------------------------------ 【10】
    api("【10】 wuf.groups()", "vector<vector<int>> … グループごとの要素の一覧",
        "グループは最小要素の昇順、各グループの中も昇順に並ぶ", "O(n α(n))  ※ ループの中で何度も呼ばない");
    SHOW(wuf.groups(), "");
    note("使用例: グループごとに、その中の最大値と最小値の差を求める");
    code("for (const auto& g : wuf.groups()) {");
    code("    long long lo = 0, hi = 0;  // 根の weight は 0 なので 0 から始めてよい");
    code("    for (int v : g) { lo = std::min(lo, wuf.weight(v)); hi = std::max(hi, wuf.weight(v)); }");
    code("    std::cout << \"group of \" << g[0] << \" : \" << hi - lo << \"\\n\";");
    code("}");
    for (const auto& g : wuf.groups()) {
        long long lo = 0, hi = 0;
        for (int v : g) {
            lo = std::min(lo, wuf.weight(v));
            hi = std::max(hi, wuf.weight(v));
        }
        std::cout << "      group of " << g[0] << " : " << hi - lo << "\n";
    }

    // ------------------------------------------------------------------ 【11】
    api("【11】 wuf.members(x)", "vector<int> … x と同じグループの要素の一覧",
        "x のグループの要素だけを取り出す。※ 並び順は不定（x から始まる）。昇順が必要なら sort",
        "O(x のグループのサイズ)  ※ groups() と違って他のグループは見ないので速い");
    SHOW(wuf.members(0), "");
    SHOW(wuf.members(2), "どの頂点から呼ぶかで順番が変わる");
    RUN(std::vector<int> mem = wuf.members(1));
    RUN(std::sort(mem.begin(), mem.end()));
    SHOW(mem, "sort 後（昇順）");
    note("使用例: 頂点 0 から見た同じグループの各頂点の値の差");
    RUN(for (int v : wuf.members(0)) std::cout << "      diff(0, " << v << ") = " << wuf.diff(0, v) << "\n");

    // ------------------------------------------------------------------ 【12】
    api("【12】 wuf.consistent(x, y, w)", "bool … 矛盾しないなら true",
        "「値[y] - 値[x] = w」を足しても矛盾しないかを調べるだけ（登録はしない）", "ならし O(α(n))");
    note("同じグループなら diff(x, y) == w かどうか。別グループなら（まだ何も決まっていないので）常に true");
    SHOW(wuf.consistent(0, 2, 8), "値[2] - 値[0] = 8 なので矛盾しない");
    SHOW(wuf.consistent(0, 2, 7), "8 ≠ 7 なので矛盾する");
    SHOW(wuf.consistent(2, 0, -8), "向きを逆にしたもの");
    SHOW(wuf.consistent(0, 3, 100), "別グループ → true");
    SHOW(wuf.same(0, 3), "consistent は登録しないので別グループのまま");

    // ------------------------------------------------------------------ 【13】
    api("【13】 wuf.unite_check(x, y, w)", "bool … 登録できた（矛盾しない）なら true / 矛盾するなら false",
        "矛盾しなければ unite と同じく登録・併合する。矛盾するなら何もしない（consistent() のフラグも変えない）",
        "ならし O(α(n))");
    note("矛盾する情報は捨てて先へ進む問題向け（典型3）。unite との違いは返り値の意味:");
    note("  unite       … 新しく併合したら true（同じグループなら矛盾しなくても false）");
    note("  unite_check … 矛盾しなければ true（同じグループで矛盾しない場合も true）");
    SHOW(wuf.unite_check(2, 3, 10), "値[3] = 値[2] + 10 として併合");
    state(wuf);
    SHOW(wuf.unite_check(0, 4, 999), "値[4] - 値[0] = 8 + 10 - 2 = 16 ≠ 999 → 捨てる");
    SHOW(wuf.unite_check(0, 4, 16), "同じグループで矛盾しない → true");
    SHOW(wuf.unite(0, 4, 16), "同じ情報でも unite なら false（新しく併合していないため）");
    SHOW(wuf.consistent(), "unite_check で矛盾しても、フラグは true のまま");

    // ------------------------------------------------------------------ 【14】
    api("【14】 wuf.consistent()", "bool … これまで unite に渡した情報が全て矛盾なしなら true",
        "制約を全部 unite してから、最後に 1 回だけ確認する使い方ができる（典型1）", "O(1)");
    note("※ 矛盾した unite は false を返すだけで、値は何も変わらない（最初に登録された情報が残る）");
    SHOW(wuf.consistent(), "ここまでは矛盾なし");
    SHOW(wuf.unite(1, 4, 0), "値[4] - 値[1] = 11 ≠ 0 → 矛盾。false が返る");
    SHOW(wuf.consistent(), "矛盾があったことが記録された");
    SHOW(wuf.diff(1, 4), "値は元のまま");
    SHOW(wuf.unite(1, 4, 11), "その後に正しい情報を入れても");
    SHOW(wuf.consistent(), "一度 false になったら戻らない（init で戻る）");
    state(wuf);

    // ------------------------------------------------------------------ 【15】
    api("【15】 wuf.group_ids()", "vector<int> … ids[v] = 頂点 v のグループ番号（0 ~ group_count()-1）",
        "番号は groups() の添字と一致する。0 から詰めた番号なので配列の添字に使える", "O(n α(n))");
    SHOW(wuf.group_ids(), "頂点 0 ~ 5 のグループ番号");

    // ------------------------------------------------------------------ 【16】
    api("【16】 wuf.roots()", "vector<int> … 各グループの根（昇順）", "グループごとに 1 回ずつ処理したいときに使う", "O(n)");
    SHOW(wuf.roots(), "");
    RUN(for (int r : wuf.roots()) std::cout << "      root " << r << " : size " << wuf.size(r) << "\n");

    // ------------------------------------------------------------------ 【17】
    api("【17】 wuf.add()", "int … 追加した頂点の番号（= 追加前の n()）",
        "1 人だけのグループとして頂点を 1 個追加する（weight は 0）", "ならし O(1)");
    SHOW(wuf.add(), "頂点 6 を追加");
    SHOW(wuf.n(), "");
    SHOW(wuf.unite(5, 6, 3), "追加した頂点も普通に使える");
    SHOW(wuf.diff(5, 6), "");

    // ------------------------------------------------------------------ 【18】
    api("【18】 wuf.init(n)", "なし (void)", "頂点数 n で作り直す（全ての情報と consistent() のフラグをリセット）",
        "O(n)");
    RUN(wuf.init(3));
    SHOW(wuf.consistent(), "フラグも true に戻る");
    state(wuf);

    // ================================================================== 典型パターン
    head("典型1: 制約が全て同時に満たせるか（ABC 087 D People on a Line）");
    note("入力: N M のあと「L R D」が M 行（x[R] - x[L] = D、1-indexed）。全部 unite して最後に consistent()");
    note("提出用のコード:");
    code("int N, M;");
    code("std::cin >> N >> M;");
    code("WeightedUnionFind<long long> p(N);");
    code("for (int i = 0; i < M; ++i) {");
    code("    int l, r;");
    code("    long long d;");
    code("    std::cin >> l >> r >> d;");
    code("    p.unite(l - 1, r - 1, d);  // x[r] - x[l] = d（0-indexed に直す）");
    code("}");
    code("std::cout << (p.consistent() ? \"Yes\" : \"No\") << \"\\n\";");
    note("↓ 同じ処理を solve(入力の文字列) にして、2 通りの入力で実行した結果");
    auto solve = [](const std::string& input) -> std::string {
        std::istringstream in(input);
        int N, M;
        in >> N >> M;
        WeightedUnionFind<long long> p(N);
        for (int i = 0; i < M; ++i) {
            int l, r;
            long long d;
            in >> l >> r >> d;
            p.unite(l - 1, r - 1, d);
        }
        return p.consistent() ? "Yes" : "No";
    };
    SHOW(solve("3 3\n1 2 1\n2 3 1\n1 3 2\n"), "1 + 1 = 2 なので全て満たせる");
    SHOW(solve("3 3\n1 2 1\n2 3 1\n1 3 5\n"), "1 + 1 ≠ 5 なので満たせない");

    head("典型2: 制約を満たす値を 1 組作る（ABC 373 D Hidden Weights）");
    note("制約を全部 unite した後、各頂点 i について weight(i) を並べればよい");
    {
        RUN(WeightedUnionFind<long long> h(5));
        RUN(h.unite(0, 1, 3));
        RUN(h.unite(1, 2, -5));
        RUN(h.unite(3, 4, -2));
        RUN(std::vector<long long> ans(5));
        RUN(for (int i = 0; i < 5; ++i) ans[i] = h.weight(i));
        SHOW(ans, "値[1]-値[0]=3, 値[2]-値[1]=-5, 値[4]-値[3]=-2 を全て満たす");
        note("値を 0 以上にしたいとき: グループごとに最小値を引く（差は変わらない）");
        code("for (const auto& g : h.groups()) {");
        code("    long long mn = ans[g[0]];");
        code("    for (int v : g) mn = std::min(mn, ans[v]);");
        code("    for (int v : g) ans[v] -= mn;");
        code("}");
        for (const auto& g : h.groups()) {
            long long mn = ans[g[0]];
            for (int v : g) mn = std::min(mn, ans[v]);
            for (int v : g) ans[v] -= mn;
        }
        SHOW(ans, "各グループの最小値が 0 になった");
    }

    head("典型3: 矛盾する情報は捨てて、採用したものの番号を答える（ABC 328 F Good Set Query）");
    {
        RUN(WeightedUnionFind<long long> s(3));
        RUN(std::vector<std::array<long long, 3>> qs = {{0, 1, 2}, {1, 2, 3}, {0, 2, 4}, {0, 2, 5}});
        note("↑ {x, y, w} … 値[y] - 値[x] = w");
        RUN(std::vector<int> accepted);
        RUN(for (int i = 0; i < 4; ++i) if (s.unite_check((int)qs[i][0], (int)qs[i][1], qs[i][2])) accepted.push_back(i));
        SHOW(accepted, "番号 2（値[2]-値[0]=4）は 2 + 3 = 5 と矛盾するので不採用");
    }

    head("典型4: 重みを XorWeight にして、2 色の塗り分け・二部グラフ判定（ABC 327 D など）");
    note("unite(a, b, 1) =「a と b は違う色」、unite(a, b, 0) =「a と b は同じ色」");
    note("diff(a, b).v が 0 なら同じ色、1 なら違う色（.v で中身の整数を取り出す）");
    {
        RUN(WeightedUnionFind<XorWeight> col(4));
        SHOW(col.unite(0, 1, 1), "0 と 1 は違う色");
        SHOW(col.unite(1, 2, 1), "1 と 2 は違う色");
        SHOW(col.diff(0, 2).v, "0 → 0 と 2 は同じ色");
        SHOW(col.consistent(0, 2, 1), "0 と 2 を違う色にはできない");
        SHOW(col.unite(2, 3, 0), "2 と 3 は同じ色");
        SHOW(col.diff(0, 3).v, "0 と 3 も同じ色");
        SHOW(col.unite(3, 0, 1), "3 と 0 は違う色 → 矛盾");
        SHOW(col.consistent(), "2 色では塗り分けられない");
    }
    note("二部グラフ判定: 全ての辺 (a, b) を unite(a, b, 1) して consistent() が true なら二部グラフ");
    {
        RUN(std::vector<std::array<int, 2>> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
        RUN(WeightedUnionFind<XorWeight> bip(4));
        RUN(for (const auto& e : edges) bip.unite(e[0], e[1], 1));
        SHOW(bip.consistent(), "長さ 4 の閉路は二部グラフ");
        RUN(bip.unite(0, 2, 1));
        SHOW(bip.consistent(), "辺 0-2 を足すと長さ 3 の閉路ができて二部グラフでなくなる");
    }

    head("典型5: 自作の型で 2 次元の相対位置を扱う（ABC 320 D Relative Position）");
    note("T は「T() が 0」「+」「-」「==」があれば何でもよい。ここでは (x, y) の組 Vec2");
    note("（Vec2 の定義はこのファイルの wuf_demo::Vec2 を参照。modint も同じ条件を満たせば使える）");
    {
        RUN(WeightedUnionFind<Vec2> pos(4));
        SHOW(pos.unite(0, 1, Vec2(2, -1)), "人 1 は人 0 から見て (+2, -1)");
        SHOW(pos.unite(1, 2, Vec2(-3, 4)), "人 2 は人 1 から見て (-3, +4)");
        SHOW(pos.diff(0, 2), "人 0 から見た人 2 の位置");
        SHOW(pos.consistent(0, 2, Vec2(-1, 3)), "");
        SHOW(pos.same(0, 3), "人 3 の位置は決まらない");
        note("人 0 を原点としたとき各人の座標: same(0, i) なら diff(0, i)、そうでなければ undecidable");
    }

    head("類題メモ");
    note("制約の矛盾判定         : ABC 087 D People on a Line");
    note("値の構成               : ABC 373 D Hidden Weights");
    note("矛盾するものを捨てる   : ABC 328 F Good Set Query");
    note("2 次元の相対位置       : ABC 320 D Relative Position");
    note("2 色塗り分け (XOR)     : ABC 327 D Good Tuple Problem");
}

}


// =============================================================================
//  weighted_union_find_from のテストとデモ（第二世代で追加）
// =============================================================================
static const char* TEST_TAG = "weighted_union_find_from";
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

namespace wuff_test {
static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 12), m = rnd(rng, 0, 2 * n);
        // 半分は「値」を決めてから辺を作る（矛盾なし）、半分はでたらめ（矛盾することがある）
        const bool good = rng() & 1;
        std::vector<long long> val(n);
        for (auto& x : val) x = rnd(rng, -20, 20);
        Graph G(n);
        WeightedUnionFind<long long> ex(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            const long long w = good ? val[b] - val[a] : rnd(rng, -5, 5);
            (rng() & 1) ? G.add_undirected(a, b, w) : G.add_edge(a, b, w);
            ex.unite(a, b, w);
        }
        WeightedUnionFind<long long> w = weighted_union_find_from(G);
        CHECK(w.n() == n && w.group_count() == ex.group_count() && w.groups() == ex.groups() && w.consistent() == ex.consistent());
        if (good) CHECK(w.consistent());
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                if (w.same(a, b)) {
                    CHECK(w.diff(a, b) == ex.diff(a, b));
                    if (good) CHECK(w.diff(a, b) == val[b] - val[a]);
                }
    }
    std::cout << "[weighted_union_find_from] random_test: OK (" << cases << " cases)\n";
}
}  // namespace wuff_test

#undef SHOW
#undef RUN
namespace wuff_demo {
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
    banner("weighted_union_find_from（第二世代で追加）    見方:  コード  ->  返り値    // 意味");
    api("weighted_union_find_from(G)", "WeightedUnionFind<Cost>",
        "G の全辺 u → v（重み w）を「値[v] - 値[u] = w」として登録したもの。矛盾があれば consistent() が false", "O(N + M α(N))");
    code("// ABC087 D（条件「x[R] - x[L] = D」が M 個）: Graph G = read_directed_w(N, M); で読み、weighted_union_find_from(G).consistent()");
    RUN(Graph G(4));
    RUN(G.add_edge(0, 1, 3); G.add_edge(1, 2, -1); G.add_edge(0, 2, 2));
    RUN(auto w = weighted_union_find_from(G));
    SHOW(w.consistent(), "値[1] - 値[0] = 3、値[2] - 値[1] = -1、値[2] - 値[0] = 2 は矛盾しない");
    SHOW(w.diff(0, 2), "値[2] - 値[0]");
    SHOW(w.same(0, 3), "3 はどこともつながっていない");
    RUN(G.add_edge(2, 0, 5));
    SHOW(weighted_union_find_from(G).consistent(), "値[0] - 値[2] = 5 は矛盾");
    note("無向辺（add_undirected / read_undirected_w）も「値[to] - 値[from] = w」（入力の向き）として登録する。");
    note("差の条件が「≤」（不等式）なら 02_bellman_ford（牛ゲー）。");
}
}  // namespace wuff_demo

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
    wuf_test::random_test();   // ランダムテスト（失敗したらその行を表示して終了）
    wuf_test::benchmark();     // 速度計測
    wuff_test::random_test();  // weighted_union_find_from のテスト
    wuf_demo::demo();          // 使い方デモ（9/23 版）
    wuff_demo::demo();         // weighted_union_find_from のデモ
    return 0;
}
