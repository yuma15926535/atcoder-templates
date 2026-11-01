// =============================================================================
//  Union-Find（素集合データ構造）                                --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 03 Union-Find・最小全域木
//  依存       : 00_graph_base（同梱済み。union_find_from(G) のため）
//  由来       : 9/23 作成の union_find.cpp（union_find_templates）をそのまま採用し、
//               第二世代の形（ファイル名・#ifndef ガード・グラフ基盤の同梱）に合わせて union_find_from(G) を追加した
//  頂点番号   : 0-indexed（0 ~ n-1）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  高速化のポイント
//    - union by size + path halving（経路半分化）… ならし O(α(n))。再帰なし
//    - par_ 1 本で「親の番号 / 根ならサイズの負値」を兼用し、メモリアクセスを減らす
//    - 同じグループの要素を循環リスト nxt_ でつなぎ、members(x) を O(そのグループのサイズ) に
//    - group_count() は併合時に数えておくので O(1)
//    - groups() は各グループのサイズ分を先に reserve して再確保を避ける
//
//  機能一覧                                                         計算量
//    UnionFind uf(n)      頂点数 n で作る（最初は全員バラバラ）     O(n)
//    uf.find(x)           x の根（グループの代表元）                ならし O(α(n))
//    uf.unite(x, y)       併合。新しく併合したら true               ならし O(α(n))
//    uf.same(x, y)        同じグループなら true                     ならし O(α(n))
//    uf.size(x)           x のグループの要素数                      ならし O(α(n))
//    uf.group_count()     グループの個数                            O(1)
//    uf.groups()          全グループ（最小要素の昇順・各グループ内も昇順）  O(n α(n))
//    uf.members(x)        x と同じグループの要素一覧（順不同）      O(そのグループのサイズ)
//    ---- 追加機能 ----
//    uf.unite(x, y, f)    併合したときだけ f(新しい根, 吸収された根) を呼ぶ   ならし O(α(n))
//                         → グループごとのデータ（和・最小値など）を根で管理できる
//    uf.edge_count(x)     x のグループで呼ばれた unite の回数（= 辺の本数）   ならし O(α(n))
//                         → 「辺の本数 >= 頂点数」なら閉路あり、「== 頂点数 - 1」なら木
//    uf.group_ids()       各頂点のグループ番号 0 ~ k-1（groups() の添字と一致）  O(n α(n))
//    uf.roots()           各グループの根の一覧（昇順）              O(n)
//    uf.add()             頂点を 1 個追加し、その番号を返す         ならし O(1)
//    uf.init(n)           頂点数 n で作り直す（全リセット）         O(n)
//    uf.n()               頂点数                                    O(1)
//    ---- 第二世代で追加 ----
//    union_find_from(G)   グラフ G の全辺を unite した UnionFind（edge_count は成分ごとの辺の本数）  O(N + M α)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤も入っている）
//    uf_test::random_test()  … 愚直解と全機能を突き合わせるランダムテスト（9/23 版のまま）
//    uf_test::benchmark()    … 速度計測（9/23 版のまま）
//    uf_demo::demo()         … 全機能の使い方デモ（9/23 版のまま）
//    uff_test / uff_demo     … union_find_from のテストとデモ（第二世代で追加）
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
#ifndef GG_UNION_FIND
#define GG_UNION_FIND
// ---- Union-Find（9/23 版をそのまま採用。グラフ・グリッド第二世代では union_find_from(G) を追加） ----
// 必要なヘッダ: <vector>, <utility>
struct UnionFind {
  public:
    UnionFind() : UnionFind(0) {}
    explicit UnionFind(int n) { init(n); }

    // 頂点数 n で作り直す（全頂点がバラバラの状態に戻る）  O(n)
    void init(int n) {
        n_ = n;
        cnt_ = n;
        par_.assign(n, -1);
        edge_.assign(n, 0);
        nxt_.resize(n);
        for (int i = 0; i < n; ++i) nxt_[i] = i;
    }

    // 頂点を 1 個追加し、その番号（= 追加前の n()）を返す  ならし O(1)
    int add() {
        par_.push_back(-1);
        edge_.push_back(0);
        nxt_.push_back(n_);
        ++cnt_;
        return n_++;
    }

    // 頂点数
    int n() const { return n_; }

    // x の根（代表元）。同じグループなら同じ値を返す  ならし O(α(n))
    int find(int x) {
        // path halving: 親を祖父に付け替えながら登る（再帰なし・1 パス）
        while (par_[x] >= 0) {
            const int p = par_[x];
            const int g = par_[p];
            if (g < 0) return p;
            par_[x] = g;
            x = g;
        }
        return x;
    }

    // x と y のグループを併合する。新しく併合したら true、既に同じなら false  ならし O(α(n))
    bool unite(int x, int y) {
        x = find(x);
        y = find(y);
        if (x == y) {
            ++edge_[x];
            return false;
        }
        link(x, y);
        return true;
    }

    // 併合し、新しく併合したときだけ f(新しい根, 吸収された根) を呼ぶ  ならし O(α(n))
    //   例: auto merge = [&](int root, int child) { sum[root] += sum[child]; };
    //       uf.unite(a, b, merge);
    template <class F>
    bool unite(int x, int y, F&& f) {
        x = find(x);
        y = find(y);
        if (x == y) {
            ++edge_[x];
            return false;
        }
        link(x, y);  // x が新しい根、y が吸収された根になる
        f(x, y);
        return true;
    }

    // x と y が同じグループなら true  ならし O(α(n))
    bool same(int x, int y) { return find(x) == find(y); }

    // x が属するグループの要素数  ならし O(α(n))
    int size(int x) { return -par_[find(x)]; }

    // グループの個数  O(1)
    int group_count() const { return cnt_; }

    // x のグループで呼ばれた unite の回数（辺の本数。false だった unite・自己ループも数える）
    int edge_count(int x) { return edge_[find(x)]; }

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
    int n_ = 0;               // 頂点数
    int cnt_ = 0;             // グループ数
    std::vector<int> par_;    // 根なら -(サイズ)、それ以外は親の番号
    std::vector<int> nxt_;    // 同じグループの要素をつなぐ循環リスト
    std::vector<int> edge_;   // 根のみ有効: そのグループで呼ばれた unite の回数

    // x, y は異なる根。サイズの大きい方（同じなら x）を根にして併合する。
    // 呼び出し後 x = 新しい根、y = 吸収された根
    void link(int& x, int& y) {
        if (par_[x] > par_[y]) std::swap(x, y);  // par_ はサイズの負値なので「>」で小さい方
        par_[x] += par_[y];
        par_[y] = x;
        edge_[x] += edge_[y] + 1;
        std::swap(nxt_[x], nxt_[y]);  // 2 つの循環リストをつなげて 1 つにする
        --cnt_;
    }
};
// グラフの全辺を unite した UnionFind（edge_count は成分ごとの辺の本数になる。ABC292 D）
inline UnionFind union_find_from(const Graph& G) {
    UnionFind uf(G.size());
    for (const RawEdge& e : G.es) uf.unite(e.from, e.to);
    return uf;
}
#endif  // GG_UNION_FIND
// ======================= ここまでテンプレート本体 =======================



// =============================================================================
//  ランダムテスト: 愚直な実装（ラベル書き換え方式）と全機能を突き合わせる
// =============================================================================
namespace uf_test {

// NDEBUG（VS の Release ビルド）でも無効にならないチェック
#define UF_CHECK(cond)                                                                   \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            std::cerr << "[UnionFind] TEST FAILED  line " << __LINE__ << ": " #cond "\n"; \
            std::exit(1);                                                                \
        }                                                                                \
    } while (0)

// 愚直解: lab[v] = v のグループ番号、mem[ラベル] = そのグループの要素一覧をそのまま持つ。
//   併合時は小さい方のグループの要素を全部大きい方へ移し、ラベルを書き換える。
struct Naive {
    std::vector<int> lab;
    std::vector<std::vector<int>> mem;  // ラベルごとの要素
    std::vector<int> ecnt;              // ラベルごとの unite 回数
    std::vector<long long> sum;         // ラベルごとの値の合計（コールバック検証用）
    int cnt = 0;

    void init(int n, const std::vector<long long>& val) {
        lab.resize(n);
        mem.assign(n, std::vector<int>());
        for (int i = 0; i < n; ++i) lab[i] = i, mem[i].push_back(i);
        ecnt.assign(n, 0);
        sum = val;
        cnt = n;
    }
    int add(long long val) {
        const int v = (int)lab.size();
        lab.push_back(v);
        mem.push_back(std::vector<int>(1, v));
        ecnt.push_back(0);
        sum.push_back(val);
        ++cnt;
        return v;
    }
    int n() const { return (int)lab.size(); }
    bool same(int x, int y) const { return lab[x] == lab[y]; }
    int size(int x) const { return (int)mem[lab[x]].size(); }
    bool unite(int x, int y) {
        int a = lab[x], b = lab[y];
        if (a == b) {
            ++ecnt[a];
            return false;
        }
        if (mem[a].size() < mem[b].size()) std::swap(a, b);  // b を a へ移す
        for (int v : mem[b]) lab[v] = a;
        mem[a].insert(mem[a].end(), mem[b].begin(), mem[b].end());
        mem[b].clear();
        ecnt[a] += ecnt[b] + 1;
        sum[a] += sum[b];
        --cnt;
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

// 全機能の全頂点チェック。members_all = true なら全頂点の members も確認（O(n^2)）
static void check_all(UnionFind& uf, const Naive& nv, const std::vector<long long>* sum, bool members_all) {
    const int n = nv.n();
    UF_CHECK(uf.n() == n);
    UF_CHECK(uf.group_count() == nv.cnt);

    // find: 同じグループ ⇔ 同じ根。根は自分自身が根。
    std::vector<int> root_of_label(n, -1);
    for (int v = 0; v < n; ++v) {
        const int r = uf.find(v);
        UF_CHECK(0 <= r && r < n);
        UF_CHECK(uf.find(r) == r);
        UF_CHECK(nv.lab[r] == nv.lab[v]);
        int& rl = root_of_label[nv.lab[v]];
        if (rl < 0) rl = r;
        UF_CHECK(rl == r);
        UF_CHECK(uf.size(v) == nv.size(v));
        UF_CHECK(uf.edge_count(v) == nv.ecnt[nv.lab[v]]);
        if (sum) UF_CHECK((*sum)[r] == nv.sum[nv.lab[v]]);
    }

    // groups / group_ids は並び順まで完全一致
    const std::vector<std::vector<int>> g = uf.groups();
    UF_CHECK(g == nv.groups());
    UF_CHECK((int)g.size() == uf.group_count());
    const std::vector<int> ids = uf.group_ids();
    UF_CHECK(ids == nv.ids());
    for (int i = 0; i < (int)g.size(); ++i)
        for (int v : g[i]) UF_CHECK(ids[v] == i);

    // roots: 昇順・個数 = グループ数・全て根・全て別グループ
    const std::vector<int> rs = uf.roots();
    UF_CHECK((int)rs.size() == uf.group_count());
    UF_CHECK(std::is_sorted(rs.begin(), rs.end()));
    std::vector<char> seen(n, 0);
    for (int r : rs) {
        UF_CHECK(uf.find(r) == r);
        UF_CHECK(!seen[nv.lab[r]]);
        seen[nv.lab[r]] = 1;
    }

    // members: 集合として一致・先頭は x
    if (members_all) {
        for (int v = 0; v < n; ++v) {
            std::vector<int> m = uf.members(v);
            UF_CHECK(!m.empty() && m[0] == v);
            std::sort(m.begin(), m.end());
            UF_CHECK(m == nv.group_of(v));
        }
    }
}

// ランダム操作列を 1 本流す
//   use_cb = true なら全ての併合をコールバック版で行い、根で管理する合計値も検証する
static long long run_case(std::mt19937_64& rng, int n, int q, bool use_cb) {
    std::vector<long long> val(n);
    for (auto& x : val) x = rnd(rng, -1000000, 1000000);
    UnionFind uf(n);
    Naive nv;
    nv.init(n, val);
    std::vector<long long> sum = val;  // テスト対象側: sum[根] = グループの合計

    const bool small = (n <= 60);
    const int check_every = std::max(1, n / 2);
    for (int t = 0; t < q; ++t) {
        const int cur = nv.n();
        const int type = rnd(rng, 0, 99);
        if (cur == 0 || type < 2) {  // add
            const long long x = rnd(rng, -1000000, 1000000);
            const int a = uf.add(), b = nv.add(x);
            UF_CHECK(a == b);
            sum.push_back(x);
            continue;
        }
        // たまに同じグループから選んで「false の unite」も十分に起こす
        int x = rnd(rng, 0, cur - 1), y = rnd(rng, 0, cur - 1);
        if (rnd(rng, 0, 3) == 0) {
            const std::vector<int>& m = nv.mem[nv.lab[x]];
            y = m[rnd(rng, 0, (int)m.size() - 1)];
        }
        if (type < 45) {  // unite
            const int rx = uf.find(x), ry = uf.find(y);
            const int sx = uf.size(x), sy = uf.size(y);
            const bool expect = nv.unite(x, y);
            bool got;
            if (use_cb) {
                int calls = 0;
                if (t & 1) {  // ラムダを一時オブジェクトで渡す
                    got = uf.unite(x, y, [&](int r, int c) {
                        ++calls;
                        UF_CHECK(r != c && uf.find(r) == r && uf.find(c) == r);
                        sum[r] += sum[c];
                    });
                } else {  // ラムダを変数で渡す
                    auto f = [&](int r, int c) {
                        ++calls;
                        UF_CHECK(r != c && uf.find(r) == r && uf.find(c) == r);
                        sum[r] += sum[c];
                    };
                    got = uf.unite(x, y, f);
                }
                UF_CHECK(calls == (got ? 1 : 0));
            } else {
                got = uf.unite(x, y);
            }
            UF_CHECK(got == expect);
            if (got) {
                // union by size: 大きい方（同じなら x 側）の根が新しい根になる
                const int r = uf.find(x);
                UF_CHECK(r == (sx >= sy ? rx : ry));
                UF_CHECK(uf.size(y) == sx + sy);
            }
        } else if (type < 58) {  // same
            UF_CHECK(uf.same(x, y) == nv.same(x, y));
        } else if (type < 68) {  // find
            UF_CHECK((uf.find(x) == uf.find(y)) == nv.same(x, y));
            UF_CHECK(nv.lab[uf.find(x)] == nv.lab[x]);
        } else if (type < 76) {  // size
            UF_CHECK(uf.size(x) == nv.size(x));
        } else if (type < 82) {  // edge_count
            UF_CHECK(uf.edge_count(x) == nv.ecnt[nv.lab[x]]);
        } else if (type < 88) {  // members
            std::vector<int> m = uf.members(x);
            UF_CHECK(m[0] == x);
            std::sort(m.begin(), m.end());
            UF_CHECK(m == nv.group_of(x));
        } else if (type < 91) {  // group_count
            UF_CHECK(uf.group_count() == nv.cnt);
        } else if (type < 93) {  // groups
            UF_CHECK(uf.groups() == nv.groups());
        } else if (type < 95) {  // group_ids
            UF_CHECK(uf.group_ids() == nv.ids());
        } else if (type < 96) {  // init（作り直し）
            const int m = rnd(rng, 0, std::max(1, n));
            std::vector<long long> v2(m);
            for (auto& z : v2) z = rnd(rng, -1000000, 1000000);
            uf.init(m);
            nv.init(m, v2);
            sum = v2;
        } else {  // n / roots
            UF_CHECK(uf.n() == nv.n());
            UF_CHECK((int)uf.roots().size() == nv.cnt);
        }
        if (t % check_every == 0) check_all(uf, nv, use_cb ? &sum : nullptr, small);
    }
    check_all(uf, nv, use_cb ? &sum : nullptr, nv.n() <= 300);
    return q;
}

// 形の決まったケース: 一直線・スター・二項木型（最も深くなる併合順）・自己ループ
static void structured_cases() {
    for (int n : {1, 2, 3, 7, 8, 64, 1000, 1024}) {
        // 一直線に併合
        {
            UnionFind uf(n);
            Naive nv;
            nv.init(n, std::vector<long long>(n, 1));
            for (int i = 0; i + 1 < n; ++i) UF_CHECK(uf.unite(i, i + 1) == nv.unite(i, i + 1));
            check_all(uf, nv, nullptr, n <= 64);
            UF_CHECK(uf.group_count() == 1 && uf.size(n - 1) == n && uf.edge_count(0) == n - 1);
        }
        // 逆順に併合 + 全員を頂点 0 とつなぐスター
        {
            UnionFind uf(n);
            Naive nv;
            nv.init(n, std::vector<long long>(n, 1));
            for (int i = n - 1; i >= 1; --i) UF_CHECK(uf.unite(0, i) == nv.unite(0, i));
            check_all(uf, nv, nullptr, n <= 64);
        }
        // 二項木型: 同じサイズ同士を根で併合し続ける（木の高さが log2(n) になる）
        {
            UnionFind uf(n);
            Naive nv;
            nv.init(n, std::vector<long long>(n, 1));
            for (int s = 1; s < n; s *= 2)
                for (int i = 0; i + s < n; i += 2 * s) UF_CHECK(uf.unite(i, i + s) == nv.unite(i, i + s));
            check_all(uf, nv, nullptr, n <= 64);
        }
        // 自己ループ・多重辺は edge_count にだけ効く
        {
            UnionFind uf(n);
            Naive nv;
            nv.init(n, std::vector<long long>(n, 1));
            for (int i = 0; i < n; ++i) {
                UF_CHECK(uf.unite(i, i) == nv.unite(i, i));
                UF_CHECK(uf.unite(i, (i * 7 + 3) % n) == nv.unite(i, (i * 7 + 3) % n));
                UF_CHECK(uf.unite((i * 7 + 3) % n, i) == nv.unite((i * 7 + 3) % n, i));
            }
            check_all(uf, nv, nullptr, n <= 64);
        }
    }
    // n = 0 と、空から add だけで作るケース
    {
        UnionFind uf;
        UF_CHECK(uf.n() == 0 && uf.group_count() == 0 && uf.groups().empty() && uf.roots().empty());
        UF_CHECK(uf.group_ids().empty());
        UnionFind uf2(0);
        for (int i = 0; i < 5; ++i) UF_CHECK(uf2.add() == i);
        UF_CHECK(uf2.unite(0, 4) && uf2.unite(4, 2) && !uf2.unite(0, 2));
        UF_CHECK(uf2.group_count() == 3 && uf2.size(2) == 3 && uf2.edge_count(0) == 3);
    }
}

static void random_test() {
    std::mt19937_64 rng(20260923);
    long long ops = 0;
    int cases = 0;
    // 小さい n（衝突・同じグループ同士の unite が多い）を大量に
    for (int it = 0; it < 1000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 0, 12), q = rnd(rng, 0, 80);
        ops += run_case(rng, n, q, it & 1);
    }
    // 中くらいの n
    for (int it = 0; it < 100 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 13, 150);
        ops += run_case(rng, n, 5 * n, it & 1);
    }
    // 大きめの n
    for (int it = 0; it < 4 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1000, 4000);
        ops += run_case(rng, n, 4 * n, it & 1);
    }
    structured_cases();
    std::cout << "[UnionFind] random test: OK  (" << cases << " cases, " << ops << " ops + structured cases)\n";
}

// 速度計測: n 頂点に対して「unite 1 回 : same 3 回」のランダム操作を q 回
static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[UnionFind] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int ns[] = {200000, 1000000};
    for (int n : ns) {
        if (TEST_SCALE == 1 && n > 200000) break;
        const int q = 4 * n;
        std::vector<std::pair<int, int>> qs(q);
        for (auto& p : qs) p = std::make_pair((int)(rng() % n), (int)(rng() % n));
        const auto t0 = std::chrono::steady_clock::now();
        UnionFind uf(n);
        long long chk = 0;
        for (int i = 0; i < q; ++i) {
            if ((i & 3) == 0) chk += uf.unite(qs[i].first, qs[i].second);
            else chk += uf.same(qs[i].first, qs[i].second);
        }
        const auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "[UnionFind] benchmark: n = " << n << ", " << q << " ops (unite:same = 1:3)  " << (long long)ms
                  << " ms  (checksum " << chk << ")\n";
    }
}

}  // namespace uf_test


// =============================================================================
//  使い方デモ: 全 API を 1 つずつ「返り値・内容・計算量・注意点」つきで実行して見せる。
//    出力の見方:   実際に書くコード  ->  返り値    // 何を返しているか（日本語）
//                  [状態] の行 … その時点のグループ分け
// =============================================================================
namespace uf_demo {

static std::string to_s(bool b) { return b ? "true" : "false"; }
static std::string to_s(int v) { return std::to_string(v); }
static std::string to_s(long long v) { return std::to_string(v); }
static std::string to_s(const std::vector<int>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + std::to_string(v[i]);
    return s + "}";
}
static std::string to_s(const std::vector<std::vector<int>>& g) {
    std::string s = "{";
    for (size_t i = 0; i < g.size(); ++i) s += (i ? ", " : "") + to_s(g[i]);
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
// 現在のグループ分けを 1 行で表示
static void state(UnionFind& uf) {
    std::string s;
    for (const auto& g : uf.groups()) {
        s += "{";
        for (size_t i = 0; i < g.size(); ++i) s += (i ? "," : "") + std::to_string(g[i]);
        s += "} ";
    }
    std::cout << "    [状態] " << s << " (グループ数 " << uf.group_count() << ")\n";
}
template <class V>
static void show(const char* expr, const V& v, const std::string& comment) {  // コード -> 返り値 // 説明
    std::string line = "  " + pad(expr, 30) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 44) + "  // " + comment;
    std::cout << line << "\n";
}
// SHOW(式, "説明") : 式そのものと返り値と説明を表示
#define SHOW(expr, comment) uf_demo::show(#expr, (expr), comment)
// RUN(文) : 文を表示してから実行（宣言にも使える）
#define RUN(...)                  \
    uf_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

// 最初に表示する API 一覧（チートシート）
static void cheat_sheet() {
    struct Row {
        const char* use;
        const char* ret;
        const char* what;
    };
    const Row rows[] = {
        {"UnionFind uf(n);", "", "頂点 0 ~ n-1 で作る（最初は全員が別々のグループ）  O(n)"},
        {"uf.n()", "int", "頂点数  O(1)"},
        {"uf.unite(x, y)", "bool", "x と y のグループを併合。新しく併合したら true  O(α(n))"},
        {"uf.same(x, y)", "bool", "同じグループなら true  O(α(n))"},
        {"uf.find(x)", "int", "x の根（代表元）。同じグループなら同じ値  O(α(n))"},
        {"uf.size(x)", "int", "x のグループの要素数  O(α(n))"},
        {"uf.group_count()", "int", "グループの個数  O(1)"},
        {"uf.groups()", "vector<vector<int>>", "全グループ（最小要素の昇順・中も昇順）  O(n α(n))"},
        {"uf.members(x)", "vector<int>", "x と同じグループの要素（順不同）  O(サイズ)"},
        {"uf.edge_count(x)", "int", "x のグループで呼ばれた unite の回数 = 辺の本数  O(α(n))"},
        {"uf.group_ids()", "vector<int>", "各頂点のグループ番号 0 ~ k-1（groups() の添字）  O(n α(n))"},
        {"uf.roots()", "vector<int>", "各グループの根（昇順）  O(n)"},
        {"uf.unite(x, y, f)", "bool", "併合したときだけ f(新しい根, 吸収された根) を呼ぶ  O(α(n))"},
        {"uf.add()", "int", "頂点を 1 個追加し、その番号を返す  ならし O(1)"},
        {"uf.init(n)", "void", "頂点数 n で作り直す（全リセット）  O(n)"},
    };
    std::cout << "\n■ API 一覧（x, y は頂点番号 0 ~ n-1。α(n) はほぼ定数。詳しくは下の【1】~【15】）\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, 20) << pad(ret, 24) << r.what << "\n";
    }
}

static void demo() {
    banner("Union-Find 使い方デモ    見方:  コード  ->  返り値    // 意味");
    cheat_sheet();

    // ------------------------------------------------------------------ 【1】
    api("【1】 UnionFind uf(n);   /   UnionFind uf;（空で作って後から init(n) や add()）", "（コンストラクタ）",
        "頂点 0 ~ n-1 を作る。最初は全員が 1 人ずつ別々のグループ", "O(n)");
    RUN(UnionFind uf(8));
    state(uf);
    RUN(UnionFind empty_uf);
    SHOW(empty_uf.n(), "空で作ると頂点数 0");
    note("1-indexed の入力は uf.unite(a - 1, b - 1) のように 0-indexed に直すのが基本。");
    note("UnionFind uf(n + 1) にして 1 ~ n をそのまま使ってもよいが、使わない頂点 0 も");
    note("1 グループとして数えられるので、グループ数は group_count() - 1 になる点に注意");

    // ------------------------------------------------------------------ 【2】
    api("【2】 uf.n()", "int", "頂点数（add() で増える）", "O(1)");
    SHOW(uf.n(), "頂点数");

    // ------------------------------------------------------------------ 【3】
    api("【3】 uf.unite(x, y)", "bool … 新しく併合したら true / 既に同じグループなら false",
        "x のグループと y のグループを 1 つにまとめる（x と y の順番はどちらでもよい）", "ならし O(α(n))");
    SHOW(uf.unite(0, 1), "{0} と {1} を併合 → {0,1}");
    state(uf);
    SHOW(uf.unite(2, 1), "{2} と {0,1} を併合（順番は関係ない）");
    SHOW(uf.unite(3, 4), "{3} と {4} を併合");
    SHOW(uf.unite(5, 6), "{5} と {6} を併合");
    state(uf);
    SHOW(uf.unite(0, 2), "既に同じグループ → false（何も変わらない）");
    SHOW(uf.unite(7, 7), "自分自身との unite も false");
    state(uf);
    note("使い方: if (uf.unite(a, b)) { 新しくつながったときの処理 }   返り値は無視してもよい");

    // ------------------------------------------------------------------ 【4】
    api("【4】 uf.same(x, y)", "bool", "x と y が同じグループなら true（find(x) == find(y) と同じ）", "ならし O(α(n))");
    SHOW(uf.same(0, 2), "0-1-2 とつながっているので true");
    SHOW(uf.same(0, 3), "0 と 3 は別のグループ");
    SHOW(uf.same(4, 4), "自分自身とは常に true");

    // ------------------------------------------------------------------ 【5】
    api("【5】 uf.find(x)", "int … x の根（グループの代表元）の頂点番号",
        "同じグループなら同じ値が返る。値は「比較」か「配列の添字」に使う", "ならし O(α(n))");
    note("※ どの頂点が根になるかに意味はない（大きい方のグループの根が残る）");
    note("※ unite すると根が変わることがある。find の結果を保存したら unite の後に取り直すこと");
    SHOW(uf.find(0), "{0,1,2} の根");
    SHOW(uf.find(2), "同じグループなので find(0) と同じ値");
    SHOW(uf.find(4), "{3,4} の根");
    SHOW(uf.find(7), "1 人だけのグループは自分自身が根");
    SHOW(uf.find(3) == 3, "find(x) == x なら x は根");
    note("根が変わる例（別の UnionFind で確認）:");
    {
        RUN(UnionFind t(4));
        RUN(t.unite(2, 3));
        SHOW(t.find(3), "{2,3} の根は 2");
        RUN(t.unite(0, 1));
        RUN(t.unite(1, 3));
        SHOW(t.find(3), "{0,1} と併合した結果、根が 0 に変わった");
    }

    // ------------------------------------------------------------------ 【6】
    api("【6】 uf.size(x)", "int", "x が属するグループの要素数", "ならし O(α(n))");
    SHOW(uf.size(0), "{0,1,2}");
    SHOW(uf.size(2), "同じグループなら同じ値");
    SHOW(uf.size(4), "{3,4}");
    SHOW(uf.size(7), "{7}");
    note("使用例: 最大のグループの大きさ → for (int r : uf.roots()) best = max(best, uf.size(r));");

    // ------------------------------------------------------------------ 【7】
    api("【7】 uf.group_count()", "int", "グループの個数（= 連結成分の個数）。unite が true になるたびに 1 減る",
        "O(1)（数えておいた値を返すだけ）");
    SHOW(uf.group_count(), "{0,1,2} {3,4} {5,6} {7} の 4 個");
    SHOW(uf.group_count() - 1, "全体を連結にするのに足りない辺の本数");

    // ------------------------------------------------------------------ 【8】
    api("【8】 uf.groups()", "vector<vector<int>> … グループごとの要素の一覧",
        "グループは最小要素の昇順、各グループの中も昇順に並ぶ", "O(n α(n))  ※ 毎回全頂点を見るので、ループの中で何度も呼ばない");
    SHOW(uf.groups(), "");
    SHOW((int)uf.groups().size(), "= group_count()");
    note("使用例: グループごとにループする（g.size() が要素数、g[0] が最小の要素）");
    RUN(for (const auto& g : uf.groups()) std::cout << "      size " << g.size() << ", min " << g[0] << "\n");

    // ------------------------------------------------------------------ 【9】
    api("【9】 uf.members(x)", "vector<int> … x と同じグループの要素の一覧",
        "x のグループの要素だけを取り出す。※ 並び順は不定（x から始まる）。昇順が必要なら sort",
        "O(x のグループのサイズ)  ※ groups() と違って他のグループは見ないので速い");
    SHOW(uf.members(0), "");
    SHOW(uf.members(2), "同じグループでも、どの頂点から呼ぶかで順番が変わる");
    SHOW(uf.members(7), "1 人だけのグループ");
    RUN(std::vector<int> mem = uf.members(1));
    RUN(std::sort(mem.begin(), mem.end()));
    SHOW(mem, "sort 後（昇順）");
    note("使用例: for (int v : uf.members(x)) { x と同じグループの頂点 v ごとの処理 }");

    // ------------------------------------------------------------------ 【10】
    api("【10】 uf.edge_count(x)", "int … x のグループで呼ばれた unite の回数（= 辺の本数）",
        "false だった unite（既に同じグループ）・自己ループも 1 本と数える", "ならし O(α(n))");
    note("グラフの辺を全部 unite した後に:");
    note("  edge_count(x) == size(x) - 1  ならそのグループは木");
    note("  edge_count(x) >= size(x)      ならそのグループは閉路を含む");
    SHOW(uf.edge_count(0), "unite(0,1) unite(2,1) unite(0,2) の 3 本");
    SHOW(uf.size(0), "頂点 3 個・辺 3 本 → 閉路あり");
    SHOW(uf.edge_count(3), "unite(3,4) の 1 本・頂点 2 個 → 木");
    SHOW(uf.edge_count(7), "unite(7,7) の自己ループ 1 本 → 閉路あり扱い");
    SHOW(uf.edge_count(5), "");
    SHOW(uf.unite(6, 5), "多重辺（同じ辺をもう一度）も false");
    SHOW(uf.edge_count(5), "1 本増えた");

    // ------------------------------------------------------------------ 【11】
    api("【11】 uf.group_ids()", "vector<int> … ids[v] = 頂点 v のグループ番号（0 ~ group_count()-1）",
        "番号は groups() の添字と一致する（uf.groups()[ids[v]] に v が入っている）", "O(n α(n))");
    note("find と違って 0 から詰めた番号なので、そのまま配列の添字に使える（縮約グラフを作るときなど）");
    RUN(std::vector<int> ids = uf.group_ids());
    SHOW(ids, "頂点 0 ~ 7 のグループ番号");
    RUN(std::vector<int> cnt(uf.group_count(), 0));
    RUN(for (int v = 0; v < uf.n(); ++v) cnt[ids[v]]++);
    SHOW(cnt, "グループ番号ごとの頂点数");

    // ------------------------------------------------------------------ 【12】
    api("【12】 uf.roots()", "vector<int> … 各グループの根（昇順）",
        "グループごとに 1 回ずつ処理したいときに使う（要素数は group_count()）", "O(n)");
    SHOW(uf.roots(), "");
    note("使用例: 同じグループに属する 2 頂点の組の数 = Σ size * (size - 1) / 2");
    RUN(long long pairs = 0);
    RUN(for (int r : uf.roots()) pairs += (long long)uf.size(r) * (uf.size(r) - 1) / 2);
    SHOW(pairs, "{0,1,2} で 3 組 + {3,4} で 1 組 + {5,6} で 1 組");

    // ------------------------------------------------------------------ 【13】
    api("【13】 uf.unite(x, y, f)", "bool（unite(x, y) と同じ）",
        "新しく併合したときだけ f(新しい根, 吸収された根) を呼ぶ。グループごとのデータを根の位置で管理できる",
        "ならし O(α(n)) + f の計算量");
    note("f は [&](int root, int child) { ... } の形。root が残る根、child は根でなくなる頂点");
    note("データの読み出しは data[uf.find(x)]（根の位置の値だけが正しい）");
    {
        RUN(UnionFind g(5));
        RUN(std::vector<long long> sum = {10, 20, 30, 40, 50});
        RUN(std::vector<long long> mx = {10, 20, 30, 40, 50});
        RUN(auto merge = [&](int root, int child) { sum[root] += sum[child]; mx[root] = std::max(mx[root], mx[child]); });
        SHOW(g.unite(0, 1, merge), "併合 → merge(根, 吸収された根) が呼ばれる");
        SHOW(sum[g.find(1)], "{0,1} の合計 = 10 + 20");
        SHOW(g.unite(2, 3, merge), "");
        SHOW(g.unite(3, 4, merge), "");
        SHOW(sum[g.find(2)], "{2,3,4} の合計");
        SHOW(mx[g.find(2)], "{2,3,4} の最大値");
        SHOW(g.unite(1, 0, merge), "既に同じ → false。merge は呼ばれない");
        SHOW(g.unite(4, 0, merge), "{0,1} と {2,3,4} を併合");
        SHOW(sum[g.find(0)], "全体の合計");
        SHOW(mx[g.find(0)], "全体の最大値");
        note("ラムダを直接書いてもよい: uf.unite(a, b, [&](int r, int c) { sum[r] += sum[c]; });");
    }

    // ------------------------------------------------------------------ 【14】
    api("【14】 uf.add()", "int … 追加した頂点の番号（= 追加前の n()）",
        "1 人だけのグループとして頂点を 1 個追加する。途中で頂点が増える問題用", "ならし O(1)");
    SHOW(uf.add(), "頂点 8 を追加");
    SHOW(uf.n(), "");
    SHOW(uf.group_count(), "グループも 1 個増える");
    SHOW(uf.unite(8, 7), "追加した頂点も普通に使える");
    state(uf);

    // ------------------------------------------------------------------ 【15】
    api("【15】 uf.init(n)", "なし (void)",
        "頂点数 n で作り直す（併合も辺の本数も全てリセット）。テストケースが複数ある問題などで使う", "O(n)");
    RUN(uf.init(3));
    SHOW(uf.n(), "");
    state(uf);

    // ================================================================== 典型パターン
    head("典型1: 入力を読んで連結成分の個数を答える（頂点 1 ~ N、辺 M 本、1-indexed）");
    note("ここでは入力を istringstream で代用している。提出時は in を std::cin にするだけ");
    {
        RUN(std::istringstream in("5 3\n1 2\n2 3\n4 5\n"));
        RUN(int N, M);
        RUN(in >> N >> M);
        RUN(UnionFind g(N));
        RUN(for (int i = 0; i < M; ++i) { int a, b; in >> a >> b; g.unite(a - 1, b - 1); });
        SHOW(g.group_count(), "{1,2,3} {4,5} の 2 個");
        SHOW(g.same(0, 2), "頂点 1 と 3（0-indexed で 0 と 2）はつながっている");
    }

    head("典型2: 閉路検出 … 辺 (a, b) を足すとき unite が false なら、その辺で閉路ができる");
    {
        RUN(UnionFind c(4));
        SHOW(c.unite(0, 1), "");
        SHOW(c.unite(1, 2), "");
        SHOW(c.unite(2, 0), "false → 辺 2-0 で閉路ができた");
        SHOW(c.unite(2, 3), "");
    }

    head("典型3: クラスカル法（最小全域木）… 辺をコストの小さい順に見て unite が true なら採用");
    {
        RUN(std::vector<std::array<int, 3>> edges = {{4, 0, 1}, {1, 1, 2}, {3, 0, 2}, {2, 2, 3}, {7, 1, 3}});
        note("↑ {コスト, u, v}。コストを先頭に置くと sort だけでコスト順になる");
        RUN(std::sort(edges.begin(), edges.end()));
        RUN(UnionFind mst(4));
        RUN(long long total = 0);
        RUN(for (const auto& e : edges) if (mst.unite(e[1], e[2])) total += e[0]);
        SHOW(total, "採用: 1-2(1), 2-3(2), 0-2(3) → 最小全域木のコスト");
        SHOW(mst.group_count() == 1, "全体が連結か（false なら全域木は作れない）");
    }

    head("典型4: グリッドの連結成分 … マス (i, j) を頂点番号 i * W + j に対応させる");
    {
        RUN(std::vector<std::string> grid = {"##..#", "#..##", "..#..", "##..#"});
        for (const auto& row : grid) note("  " + row);
        RUN(const int H = 4, W = 5);
        RUN(UnionFind gu(H * W));
        code("for (int i = 0; i < H; ++i)");
        code("    for (int j = 0; j < W; ++j) {");
        code("        if (grid[i][j] != '#') continue;");
        code("        if (i + 1 < H && grid[i + 1][j] == '#') gu.unite(i * W + j, (i + 1) * W + j);  // 下とつなぐ");
        code("        if (j + 1 < W && grid[i][j + 1] == '#') gu.unite(i * W + j, i * W + j + 1);    // 右とつなぐ");
        code("    }");
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) {
                if (grid[i][j] != '#') continue;
                if (i + 1 < H && grid[i + 1][j] == '#') gu.unite(i * W + j, (i + 1) * W + j);
                if (j + 1 < W && grid[i][j + 1] == '#') gu.unite(i * W + j, i * W + j + 1);
            }
        note("'.' のマスも 1 人ずつのグループとして数えられるので、根が '#' のものだけ数える");
        RUN(int comp = 0, best = 0);
        RUN(for (int r : gu.roots()) if (grid[r / W][r % W] == '#') { ++comp; best = std::max(best, gu.size(r)); });
        SHOW(comp, "'#' の連結成分の個数");
        SHOW(best, "最大の成分のマス数");
        SHOW(gu.same(0 * W + 0, 1 * W + 0), "(0,0) と (1,0) は同じ成分");
    }

    head("類題メモ");
    note("edge_count で閉路判定 : ABC 292 D Unicyclic Components");
    note("add() で頂点を増やす  : ABC 279 F BOX");
}

}  // namespace uf_demo



// =============================================================================
//  union_find_from のテストとデモ（第二世代で追加）
// =============================================================================
static const char* TEST_TAG = "union_find_from";
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

namespace uff_test {
static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 15), m = rnd(rng, 0, 2 * n);
        Graph G(n);
        UnionFind ex(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            (rng() & 1) ? G.add_undirected(a, b) : G.add_edge(a, b);
            ex.unite(a, b);
        }
        UnionFind uf = union_find_from(G);
        CHECK(uf.n() == n && uf.group_count() == ex.group_count() && uf.groups() == ex.groups());
        for (int v = 0; v < n; ++v) CHECK(uf.edge_count(v) == ex.edge_count(v) && uf.size(v) == ex.size(v));
        int sum = 0;
        for (int r : uf.roots()) sum += uf.edge_count(r);
        CHECK(sum == m);
    }
    std::cout << "[union_find_from] random_test: OK (" << cases << " cases)\n";
}
}  // namespace uff_test

#undef SHOW
#undef RUN
namespace uff_demo {
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
    banner("union_find_from（第二世代で追加）    見方:  コード  ->  返り値    // 意味");
    api("union_find_from(G)", "UnionFind", "グラフ G の全辺（G.es）を unite した UnionFind。edge_count(x) は x の成分の辺の本数（多重辺・自己ループも数える）", "O(N + M α(N))");
    RUN(Graph G(5));
    RUN(G.add_undirected(0, 1); G.add_undirected(1, 2); G.add_undirected(2, 0); G.add_undirected(3, 4));
    RUN(UnionFind uf = union_find_from(G));
    SHOW(uf.group_count(), "成分の数");
    SHOW(uf.groups(), "成分");
    SHOW(uf.edge_count(0), "{0,1,2} の辺は 3 本 = 頂点数 → 閉路がちょうど 1 つ（ABC292 D）");
    SHOW(uf.edge_count(3), "{3,4} の辺は 1 本 = 頂点数 - 1 → 木");
    note("辺を 1 本ずつ足しながら答える問題は、uf.unite(a, b) を直接使う。");
    note("1 回だけ連結成分を求めるなら 01_components の components(G) でもよい（成分ごとの頂点の一覧が CSR で取れる）。");
}
}  // namespace uff_demo

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
    uf_test::random_test();   // ランダムテスト（失敗したらその行を表示して終了）
    uf_test::benchmark();     // 速度計測
    uff_test::random_test();  // union_find_from のテスト
    uf_demo::demo();          // 使い方デモ（9/23 版）
    uff_demo::demo();         // union_find_from のデモ
    return 0;
}
