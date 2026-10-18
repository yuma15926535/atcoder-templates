// =============================================================================
//  クラスカル法（最小全域木・最大全域木・全域森）  kruskal / kruskal_max   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 03 Union-Find・最小全域木
//  依存       : 00_graph_base, 03_union_find（同梱済み）
//  頂点番号   : 0-indexed。辺の向きは無視する（G.es の from - to を無向の辺として使う）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  機能一覧                                                                       計算量
//    kruskal(G) / kruskal_max(G)        最小 / 最大全域木（非連結なら全域森）         O(M log M)
//    kruskal(G, uf) / kruskal_max(G, uf)  uf で unite 済みの頂点はつながっているものとして始める
//    結果 MstResult m:
//      m.cost        選んだ辺の重みの合計
//      m.edges       選んだ辺の番号（重みの順。同じ重みは辺番号の小さい順）
//      m.connected   全域木になったか（= グラフが連結か）
//      m.components  選び終えたあとの連結成分の数
//      m.tree(G)     選んだ辺だけのグラフ（新しい辺番号 k は m.edges[k]）。MST 上で LCA などを使うとき
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（グラフ基盤と Union-Find も入っている）
//    mst_test::random_test() / mst_test::benchmark() / mst_demo::demo()
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
#ifndef GG_KRUSKAL
#define GG_KRUSKAL
// ---- クラスカル法（グラフ・グリッド第二世代）  最小（最大）全域木・全域森。辺の向きは無視する ----
struct MstResult {
    Cost cost = 0;             // 選んだ辺の重みの合計
    std::vector<int> edges;    // 選んだ辺の番号（選んだ順 = 重みの小さい順）
    bool connected = false;    // 全域木になったか（グラフ全体が連結か）
    int components = 0;        // 選び終えたあとの連結成分の数（1 なら全域木）
    // 選んだ辺だけのグラフ（無向）。新しい辺番号 k は edges[k] に対応する（MST 上で LCA などを使うとき）
    Graph tree(const Graph& G) const {
        Graph T(G.size());
        for (int id : edges) T.add_undirected(G.es[id].from, G.es[id].to, G.es[id].w);
        return T;
    }
};
namespace gg_internal {
inline MstResult kruskal_impl(const Graph& G, UnionFind uf, bool maximize) {
    std::vector<std::pair<Cost, int>> ord(G.edge_count());  // {重み, 辺番号}（同じ重みは辺番号の小さい順）
    for (int i = 0; i < G.edge_count(); ++i) ord[i] = {G.es[i].w, i};
    if (maximize) std::sort(ord.begin(), ord.end(), [](const std::pair<Cost, int>& a, const std::pair<Cost, int>& b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
    else std::sort(ord.begin(), ord.end());
    MstResult r;
    for (const auto& p : ord) {
        const RawEdge& e = G.es[p.second];
        if (uf.unite(e.from, e.to)) r.cost += e.w, r.edges.push_back(p.second);
    }
    r.components = uf.group_count();
    r.connected = (r.components <= 1);
    return r;
}
}  // namespace gg_internal
// 最小全域木（非連結なら最小全域森。connected = false）。同じ重みは辺番号の小さい順に選ぶ
inline MstResult kruskal(const Graph& G) { return gg_internal::kruskal_impl(G, UnionFind(G.size()), false); }
// 最大全域木
inline MstResult kruskal_max(const Graph& G) { return gg_internal::kruskal_impl(G, UnionFind(G.size()), true); }
// すでにつながっている頂点（uf で unite 済み）から始める版。必ず使う辺・最初からつながっている組がある問題で
inline MstResult kruskal(const Graph& G, const UnionFind& uf) { return gg_internal::kruskal_impl(G, uf, false); }
inline MstResult kruskal_max(const Graph& G, const UnionFind& uf) { return gg_internal::kruskal_impl(G, uf, true); }
#endif  // GG_KRUSKAL
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "kruskal";
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

namespace mst_test {

// 愚直解: 辺の部分集合を全部試す（森で、成分の数が最小のもののうち、重みの合計が最小 / 最大）
static void brute(const Graph& G, const std::vector<std::pair<int, int>>& pre, bool maximize, Cost& best, int& comps) {
    const int n = G.size(), m = G.edge_count();
    comps = n + 1;
    best = 0;
    for (int mask = 0; mask < (1 << m); ++mask) {
        std::vector<int> p(n);
        std::iota(p.begin(), p.end(), 0);
        std::function<int(int)> f = [&](int x) { return p[x] == x ? x : p[x] = f(p[x]); };
        int c = n;
        for (auto [a, b] : pre)
            if (f(a) != f(b)) p[f(a)] = f(b), --c;
        bool forest = true;
        Cost sum = 0;
        for (int i = 0; i < m && forest; ++i)
            if (mask >> i & 1) {
                const RawEdge& e = G.es[i];
                if (f(e.from) == f(e.to)) forest = false;
                else p[f(e.from)] = f(e.to), --c, sum += e.w;
            }
        if (!forest) continue;
        if (c < comps || (c == comps && (maximize ? sum > best : sum < best))) comps = c, best = sum;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 1200 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 7), m = rnd(rng, 0, 10);
        Graph G(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            const Cost w = rnd(rng, -5, 9);
            (rng() & 1) ? G.add_undirected(a, b, w) : G.add_edge(a, b, w);
        }
        std::vector<std::pair<int, int>> pre;
        UnionFind uf(n);
        if (it % 2)
            for (int k = rnd(rng, 0, 2); k > 0; --k) {
                const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
                pre.emplace_back(a, b), uf.unite(a, b);
            }
        for (int mx = 0; mx < 2; ++mx) {
            const MstResult r = pre.empty() ? (mx ? kruskal_max(G) : kruskal(G)) : (mx ? kruskal_max(G, uf) : kruskal(G, uf));
            Cost best;
            int comps;
            brute(G, pre, mx, best, comps);
            CHECK(r.cost == best && r.components == comps && r.connected == (comps <= 1));
            // 選んだ辺が森で、合計が cost
            UnionFind chk = uf;
            Cost sum = 0;
            for (int id : r.edges) {
                CHECK(chk.unite(G.es[id].from, G.es[id].to));
                sum += G.es[id].w;
            }
            CHECK(sum == r.cost && chk.group_count() == comps);
            for (size_t k = 1; k < r.edges.size(); ++k) {
                const RawEdge &a = G.es[r.edges[k - 1]], &b = G.es[r.edges[k]];
                CHECK(mx ? a.w >= b.w : a.w <= b.w);
            }
            const Graph T = r.tree(G);
            CHECK(T.edge_count() == (int)r.edges.size());
            for (int k = 0; k < T.edge_count(); ++k) CHECK(T.es[k].from == G.es[r.edges[k]].from && T.es[k].to == G.es[r.edges[k]].to && T.es[k].w == G.es[r.edges[k]].w);
        }
    }
    std::cout << "[kruskal] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[kruskal] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        std::vector<int> us(n), vs(n);
        std::vector<long long> ws(n);
        for (int i = 0; i < n; ++i) us[i] = (int)(rng() % n), vs[i] = (int)(rng() % n), ws[i] = (long long)(rng() % 1000000000);
        const Graph G = graph_from_edges(n, us, vs, ws, UNDIRECTED, 0);
        auto t0 = std::chrono::steady_clock::now();
        const MstResult r = kruskal(G);
        std::cout << "[kruskal] benchmark: N = M = " << n << "  " << ms_since(t0) << " ms  (edges " << r.edges.size() << ")\n";
    }
}

}  // namespace mst_test


namespace mst_demo {
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
    banner("クラスカル法 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto m = kruskal(G);", "MstResult", "最小全域木（非連結なら全域森）"},
        {"auto m = kruskal_max(G);", "MstResult", "最大全域木"},
        {"kruskal(G, uf) / kruskal_max(G, uf)", "MstResult", "uf でつないだ頂点はつながっているものとして始める"},
        {"m.cost / m.edges", "Cost / vector<int>", "重みの合計 / 選んだ辺の番号"},
        {"m.connected / m.components", "bool / int", "全域木になったか / 連結成分の数"},
        {"m.tree(G)", "Graph", "選んだ辺だけのグラフ"},
    };
    cheat_sheet("API 一覧", rows, 38, 22);

    head("デモで使うグラフ（無向）: 0-1 (4), 1-2 (2), 0-2 (3), 2-3 (5), 1-3 (7)");
    RUN(Graph G(4));
    RUN(G.add_undirected(0, 1, 4); G.add_undirected(1, 2, 2); G.add_undirected(0, 2, 3); G.add_undirected(2, 3, 5); G.add_undirected(1, 3, 7));

    api("【1】 kruskal(G)", "MstResult", "最小全域木。重みの小さい辺から、閉路にならないものを選ぶ", "O(M log M)");
    RUN(auto m = kruskal(G));
    SHOW(m.cost, "2 + 3 + 5");
    SHOW(m.edges, "選んだ辺の番号（重みの小さい順）");
    SHOW(m.connected, "連結（全域木になった）");

    api("【2】 kruskal_max(G)", "MstResult", "最大全域木（重みの大きい辺から）", "O(M log M)");
    SHOW(kruskal_max(G).cost, "7 + 5 + 4");

    api("【3】 kruskal(G, uf)", "MstResult", "uf で unite した頂点どうしは、最初からつながっているものとして始める（必ず使う辺・無料の道がある問題）", "O(M log M)");
    RUN(UnionFind uf(4));
    RUN(uf.unite(0, 3));
    SHOW(kruskal(G, uf).cost, "0 と 3 は無料でつながっている → 2 + 3");

    api("【4】 非連結のとき", "MstResult", "全域森（各成分の最小全域木）を返し、connected = false。components に成分の数", "O(M log M)");
    RUN(Graph F(5));
    RUN(F.add_undirected(0, 1, 1); F.add_undirected(2, 3, 1));
    RUN(auto mf = kruskal(F));
    SHOW(mf.connected, "連結でない");
    SHOW(mf.components, "成分の数（あと 2 本足せば連結）");

    api("【5】 m.tree(G)", "Graph", "選んだ辺だけのグラフ（無向）。新しい辺番号 k は m.edges[k]。MST 上の 2 点間の最大の辺などを LCA と組み合わせて", "O(N + M)");
    RUN(Graph T = m.tree(G));
    SHOW(T.edge_count(), "3 本");
    SHOW(T.es, "選んだ辺");

    head("典型的な使い方");
    note("全頂点を最小コストでつなぐ                  : kruskal(G).cost（連結でなければ -1 などを出力）");
    note("超頂点（どこからでも費用 c で行ける場所）   : G.add_vertex() で頂点を足し、各頂点と重み c でつないでから kruskal");
    note("すでにある道・必ず使う辺                    : UnionFind に先に unite してから kruskal(G, uf)");
}

}  // namespace mst_demo

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
    mst_test::random_test();
    mst_test::benchmark();
    mst_demo::demo();
    return 0;
}
