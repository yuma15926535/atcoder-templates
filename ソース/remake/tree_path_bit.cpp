// =============================================================================
//  【制限付き最速版】木のパスの和（オイラーツアー + BIT：1 点更新・パスの和 / パスに加算・1 点取得 を O(log N)）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : 演算は和・xor のように「引き算ができて順番によらない」ものだけ（min・max・行列の積などは不可）。
//                           TreePathSum は 1 点更新だけ（パスへの加算はない）、TreePathAdd はパスへの加算だけ（パスの和はない）
//    できること           : TreePathSum … 頂点（辺）の値の 1 点更新 + 2 点間のパスの和・根からの和・部分木の和
//                           TreePathAdd … パスの全頂点（全辺）への加算 + 1 頂点（1 辺）の値・部分木の和
//    一般版との比較       : graph_grid_v2/05_hld.cpp の TreeQuery / TreeLazyQuery・HLDWrapperTemplate.txt はパスが O(log² N)
//                           → この版はパスの質問・パスへの加算が O(log N)（LCA 1 回 + BIT を数回）
//    制限に合わないとき   : min・max・非可換な演算、パスへの加算とパスの和の両方が要るなら graph_grid_v2/05_hld.cpp
//
//  しくみ: 部分木が連続する順（HLD の行きがけ順）に並べ、頂点 v の値を「v に入る位置に +x、v の部分木を出る位置に −x」で BIT に置く。
//          位置 v までの累積和 = 根から v までの和。パス u–v の和 = R(u) + R(v) − 2 R(lca) + 値(lca)
//          パスへの加算は u・v に +x、lca と lca の親に −x を置き、頂点の値 = 部分木の中の置いた値の和
//    例: ABC294 G（Distance Queries on a Tree）… 辺の重みの変更と 2 点間の距離 → TreePathSum(…, edge = true)
//
//  機能一覧（頂点は 0〜N−1。g は vector<vector<int>> の隣接リスト、または (n, 辺の一覧)。木でなければ止める）     計算量
//    TreePathSum<T> ts(g, root = 0, init = {}, edge = false)  init[v] = 頂点 v の値（edge なら v と親の辺の値）  O(N)
//    TreePathSum<T> ts(n, edges, root = 0, init = {}, edge = false)
//      ts.add(v, x) / ts.set(v, x) / ts.get(v)       頂点 v（edge なら v と親の辺）の値に足す / 書く / 読む   O(log N) / O(log N) / O(1)
//      ts.add_edge(u, v, x) / set_edge / get_edge    辺 u–v の値（edge のとき。u, v は隣どうし）                O(log N)
//      ts.path(u, v)                                 u〜v のパスの和（edge なら辺の和）                        O(log N)
//      ts.root_path(v)                               根〜v の和                                                O(log N)
//      ts.subtree(v)                                 v の部分木の和（edge なら v より下の辺の和）              O(log N)
//    TreePathAdd<T> ta(g, root = 0, init = {}, edge = false) / ta(n, edges, …)
//      ta.path_add(u, v, x)                          u〜v のパスの全頂点（edge なら全辺）に x を足す            O(log N)
//      ta.get(v)                                     頂点 v（edge なら v と親の辺）の値                        O(log N)
//      ta.get_edge(u, v)                             辺 u–v の値（edge のとき）                                 O(log N)
//      ta.subtree_sum(v)                             v の部分木の値の和（edge なら v より下の辺の和）           O(log N)
//    共通: lca(u, v)（O(log N)）/ depth(v) / parent(v)（根は −1）/ size()
//
//  T は long long（ふつう）/ modint など（+, −, T(long long) との * ができる型。xor は T を自分で作る）
//  速さの目安: N = 5×10^5 のランダムな木で add と path を 5×10^5 回ずつ 約 0.3 秒（一般版の HLD + BIT（O(log² N)）は 約 0.44 秒、
//              完全二分木なら 約 0.39 秒 対 約 0.60 秒）。一直線に近い木では HLD の分解が 2 本ほどで済むので同じくらい
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    tpb_test::random_test()  … 親をたどってパスを作る愚直解との突き合わせ
//    tpb_test::benchmark()    … 速度計測（HLD + BIT との比較つき）
//    tpb_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace tpb_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 根付き木（非再帰）。pre は HLD の行きがけ順（部分木 v = [pre[v], pre[v] + sz[v])）。LCA は重軽分解
struct RootedTree {
    int n = 0, root = 0;
    std::vector<int> par, dep, sz, head, pre;
    void build(const std::vector<std::vector<int>>& g, int r) {
        n = (int)g.size(), root = r;
        if (r < 0 || r >= std::max(n, 1)) fail("tree_path_bit", "根が範囲外");
        par.assign(n, -1), dep.assign(n, 0), sz.assign(n, 1), head.assign(n, 0), pre.assign(n, -1);
        if (n == 0) return;
        std::vector<int> order;
        order.reserve(n);
        std::vector<char> seen(n, 0);
        order.push_back(r), seen[r] = 1;
        for (int k = 0; k < (int)order.size(); ++k) {
            const int v = order[k];
            bool skipped_parent = false;
            for (int u : g[v]) {
                if (u < 0 || u >= n) fail("tree_path_bit", "頂点の番号が範囲外");
                if (u == par[v] && !skipped_parent) {
                    skipped_parent = true;
                    continue;
                }
                if (seen[u]) fail("tree_path_bit", "木ではない（閉路・多重辺がある）");
                seen[u] = 1, par[u] = v, dep[u] = dep[v] + 1, order.push_back(u);
            }
        }
        if ((int)order.size() != n) fail("tree_path_bit", "木ではない（つながっていない）");
        std::vector<int> heavy(n, -1);
        for (int k = n - 1; k >= 1; --k) sz[par[order[k]]] += sz[order[k]];
        for (int k = n - 1; k >= 1; --k) {
            const int v = order[k], p = par[v];
            if (heavy[p] == -1 || sz[v] > sz[heavy[p]]) heavy[p] = v;
        }
        std::vector<int> st = {r};
        head[r] = r;
        int t = 0;
        while (!st.empty()) {
            const int v = st.back();
            st.pop_back();
            for (int x = v; x != -1; x = heavy[x]) {
                pre[x] = t++;
                for (int c : g[x])
                    if (c != par[x] && c != heavy[x]) head[c] = c, st.push_back(c);
                if (heavy[x] != -1) head[heavy[x]] = head[x];
            }
        }
    }
    int lca(int u, int v) const {
        while (head[u] != head[v]) {
            if (pre[head[u]] > pre[head[v]]) {
                u = par[head[u]];
            } else {
                v = par[head[v]];
            }
        }
        return pre[u] < pre[v] ? u : v;
    }
    int child_of_edge(int u, int v) const {  // 辺 u–v の子の側
        if (par[u] == v) return u;
        if (par[v] == u) return v;
        fail("tree_path_bit", "u と v は隣どうしではない");
    }
    void check(int v) const {
        if (v < 0 || v >= n) fail("tree_path_bit", "頂点が範囲外");
    }
};
inline std::vector<std::vector<int>> to_adj(int n, const std::vector<std::pair<int, int>>& edges) {
    std::vector<std::vector<int>> g(n);
    for (auto [u, v] : edges) {
        if (u < 0 || u >= n || v < 0 || v >= n) fail("tree_path_bit", "辺の端が範囲外");
        g[u].push_back(v), g[v].push_back(u);
    }
    return g;
}
template <class T>
struct Fenwick {  // 0 始まり。prefix(r) = [0, r) の和
    int n = 0;
    std::vector<T> t;
    void init(const std::vector<T>& a) {  // O(n)
        n = (int)a.size(), t.assign(n + 1, T(0));
        for (int i = 0; i < n; ++i) {
            t[i + 1] += a[i];
            const int j = (i + 1) + ((i + 1) & -(i + 1));
            if (j <= n) t[j] += t[i + 1];
        }
    }
    void add(int i, const T& x) {
        for (++i; i <= n; i += i & -i) t[i] += x;
    }
    T prefix(int r) const {
        T s(0);
        for (; r > 0; r -= r & -r) s += t[r];
        return s;
    }
};
}  // namespace tpb_internal

// ---- 1 点更新 + パスの和（根からの和・部分木の和も） ----
template <class T = long long>
struct TreePathSum {
    TreePathSum(const std::vector<std::vector<int>>& g, int root = 0, const std::vector<T>& init = {}, bool edge = false) : edge_(edge) { build(g, root, init); }
    TreePathSum(int n, const std::vector<std::pair<int, int>>& edges, int root = 0, const std::vector<T>& init = {}, bool edge = false) : edge_(edge) {
        build(tpb_internal::to_adj(n, edges), root, init);
    }
    void add(int v, const T& x) {  // 頂点 v（edge なら v と親の辺）の値に x
        tr_.check(v);
        if (edge_ && v == tr_.root) tpb_internal::fail("TreePathSum::add", "edge = true のとき根には辺がない");
        val_[v] += x;
        const int p = tr_.pre[v];
        a_.add(p, x), a_.add(p + tr_.sz[v], T(0) - x), b_.add(p, x);
    }
    void set(int v, const T& x) { add(v, x - get(v)); }
    const T& get(int v) const {
        tr_.check(v);
        return val_[v];
    }
    void add_edge(int u, int v, const T& x) { add(edge_child(u, v), x); }
    void set_edge(int u, int v, const T& x) { set(edge_child(u, v), x); }
    const T& get_edge(int u, int v) const { return get(edge_child(u, v)); }
    T root_path(int v) const {  // 根〜v の和
        tr_.check(v);
        return a_.prefix(tr_.pre[v] + 1);
    }
    T path(int u, int v) const {  // u〜v の和（edge なら辺の和）
        const int l = lca(u, v);
        const T rl = a_.prefix(tr_.pre[l] + 1);
        const T s = a_.prefix(tr_.pre[u] + 1) + a_.prefix(tr_.pre[v] + 1) - rl - rl;
        return edge_ ? s : s + val_[l];
    }
    T subtree(int v) const {  // v の部分木の和（edge なら v より下の辺の和）
        tr_.check(v);
        const int p = tr_.pre[v];
        const T s = b_.prefix(p + tr_.sz[v]) - b_.prefix(p);
        return edge_ ? s - val_[v] : s;
    }
    int lca(int u, int v) const {
        tr_.check(u), tr_.check(v);
        return tr_.lca(u, v);
    }
    int depth(int v) const { return tr_.check(v), tr_.dep[v]; }
    int parent(int v) const { return tr_.check(v), tr_.par[v]; }
    int size() const { return tr_.n; }

   private:
    bool edge_;
    tpb_internal::RootedTree tr_;
    std::vector<T> val_;
    tpb_internal::Fenwick<T> a_, b_;  // a_: 入る +x・出る −x（根からの和）、b_: 入る +x（部分木の和）
    void build(const std::vector<std::vector<int>>& g, int root, const std::vector<T>& init) {
        const int n = (int)g.size();
        if (!init.empty() && (int)init.size() != n) tpb_internal::fail("TreePathSum", "init の長さが頂点数と違う");
        tr_.build(g, root);
        val_ = init.empty() ? std::vector<T>(n, T(0)) : init;
        if (edge_ && n > 0) val_[root] = T(0);  // 根には辺がない
        std::vector<T> da(n + 1, T(0)), db(n, T(0));
        for (int v = 0; v < n; ++v) {
            const int p = tr_.pre[v];
            da[p] += val_[v], da[p + tr_.sz[v]] -= val_[v], db[p] += val_[v];
        }
        a_.init(da), b_.init(db);
    }
    int edge_child(int u, int v) const {
        if (!edge_) tpb_internal::fail("TreePathSum::*_edge", "edge = true で作ったときだけ使える");
        tr_.check(u), tr_.check(v);
        return tr_.child_of_edge(u, v);
    }
};

// ---- パスに加算 + 1 点取得（部分木の和も） ----
template <class T = long long>
struct TreePathAdd {
    TreePathAdd(const std::vector<std::vector<int>>& g, int root = 0, const std::vector<T>& init = {}, bool edge = false) : edge_(edge) { build(g, root, init); }
    TreePathAdd(int n, const std::vector<std::pair<int, int>>& edges, int root = 0, const std::vector<T>& init = {}, bool edge = false) : edge_(edge) {
        build(tpb_internal::to_adj(n, edges), root, init);
    }
    void path_add(int u, int v, const T& x) {  // u〜v の全頂点（edge なら全辺）に x
        const int l = lca(u, v);
        put(u, x), put(v, x);
        if (edge_) {
            put(l, T(0) - x - x);
        } else {
            put(l, T(0) - x);
            if (tr_.par[l] != -1) put(tr_.par[l], T(0) - x);
        }
    }
    T get(int v) const {  // 頂点 v（edge なら v と親の辺）の値
        tr_.check(v);
        const int p = tr_.pre[v];
        return init_[v] + (a_.prefix(p + tr_.sz[v]) - a_.prefix(p));
    }
    T get_edge(int u, int v) const {
        if (!edge_) tpb_internal::fail("TreePathAdd::get_edge", "edge = true で作ったときだけ使える");
        tr_.check(u), tr_.check(v);
        return get(tr_.child_of_edge(u, v));
    }
    T subtree_sum(int v) const {  // v の部分木の値の和（edge なら v より下の辺の和）
        tr_.check(v);
        const int p = tr_.pre[v], q = p + tr_.sz[v];
        const T da = a_.prefix(q) - a_.prefix(p), db = b_.prefix(q) - b_.prefix(p);
        const T s = (isum_[q] - isum_[p]) + db - T((long long)tr_.dep[v] - 1) * da;
        return edge_ ? s - get(v) : s;
    }
    int lca(int u, int v) const {
        tr_.check(u), tr_.check(v);
        return tr_.lca(u, v);
    }
    int depth(int v) const { return tr_.check(v), tr_.dep[v]; }
    int parent(int v) const { return tr_.check(v), tr_.par[v]; }
    int size() const { return tr_.n; }

   private:
    bool edge_;
    tpb_internal::RootedTree tr_;
    std::vector<T> init_, isum_;      // isum_ = 行きがけ順に並べた init の累積和
    tpb_internal::Fenwick<T> a_, b_;  // a_: 置いた値 d、b_: d × 深さ
    void put(int v, const T& x) {
        const int p = tr_.pre[v];
        a_.add(p, x), b_.add(p, T((long long)tr_.dep[v]) * x);
    }
    void build(const std::vector<std::vector<int>>& g, int root, const std::vector<T>& init) {
        const int n = (int)g.size();
        if (!init.empty() && (int)init.size() != n) tpb_internal::fail("TreePathAdd", "init の長さが頂点数と違う");
        tr_.build(g, root);
        init_ = init.empty() ? std::vector<T>(n, T(0)) : init;
        if (edge_ && n > 0) init_[root] = T(0);
        isum_.assign(n + 1, T(0));
        for (int v = 0; v < n; ++v) isum_[tr_.pre[v] + 1] = init_[v];
        for (int i = 0; i < n; ++i) isum_[i + 1] += isum_[i];
        a_.init(std::vector<T>(n, T(0))), b_.init(std::vector<T>(n, T(0)));
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace tpb_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define TPB_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[tree_path_bit] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                      \
        }                                                                                      \
    } while (0)

// [lo, hi] の一様乱数
[[maybe_unused]] static long long rnd(std::mt19937_64& rng, long long lo, long long hi) {
    const unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return (long long)((unsigned long long)lo + (w == 0 ? rng() : rng() % w));
}
[[maybe_unused]] static long long elapsed_ms(std::chrono::steady_clock::time_point a) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
}

using ll = long long;

// 愚直解：親をたどってパスの頂点を集める
struct BruteTree {
    int n;
    std::vector<int> par, dep;
    std::vector<std::vector<int>> g;
    BruteTree(const std::vector<std::vector<int>>& g_, int root) : n((int)g_.size()), par(n, -1), dep(n, 0), g(g_) {
        std::vector<int> st = {root};
        std::vector<char> seen(n, 0);
        seen[root] = 1;
        while (!st.empty()) {
            int v = st.back();
            st.pop_back();
            for (int u : g[v])
                if (!seen[u]) seen[u] = 1, par[u] = v, dep[u] = dep[v] + 1, st.push_back(u);
        }
    }
    std::vector<int> path(int u, int v) const {  // u〜v の頂点（lca を含む）
        std::vector<int> a, b;
        while (dep[u] > dep[v]) a.push_back(u), u = par[u];
        while (dep[v] > dep[u]) b.push_back(v), v = par[v];
        while (u != v) a.push_back(u), b.push_back(v), u = par[u], v = par[v];
        a.push_back(u);
        a.insert(a.end(), b.rbegin(), b.rend());
        return a;
    }
    bool in_sub(int w, int v) const {
        while (w != -1 && w != v) w = par[w];
        return w == v;
    }
};

static std::vector<std::vector<int>> random_tree(std::mt19937_64& rng, int n, int shape) {
    std::vector<std::vector<int>> g(n);
    std::vector<int> id(n);
    std::iota(id.begin(), id.end(), 0);
    std::shuffle(id.begin(), id.end(), rng);
    for (int i = 1; i < n; ++i) {
        const int p = shape == 0 ? (int)rnd(rng, 0, i - 1) : (shape == 1 ? i - 1 : (int)rnd(rng, std::max(0, i - 3), i - 1));
        g[id[i]].push_back(id[p]), g[id[p]].push_back(id[i]);
    }
    for (auto& a : g) std::shuffle(a.begin(), a.end(), rng);
    return g;
}

static void random_test() {
    std::mt19937_64 rng(20261101);
    long long cases = 0;
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 25);
        const auto g = random_tree(rng, n, (int)rnd(rng, 0, 2));
        const int root = (int)rnd(rng, 0, n - 1);
        const bool edge = rnd(rng, 0, 1);
        std::vector<ll> init(n);
        for (auto& v : init) v = rnd(rng, -9, 9);
        BruteTree bt(g, root);
        if (edge) init[root] = 0;
        // TreePathSum
        {
            TreePathSum<ll> ts(g, root, init, edge);
            std::vector<ll> val = init;
            for (int t = 0; t < 60; ++t) {
                const int k = (int)rnd(rng, 0, 4), u = (int)rnd(rng, 0, n - 1), v = (int)rnd(rng, 0, n - 1);
                const ll x = rnd(rng, -9, 9);
                if (k == 0 && !(edge && u == root)) {
                    ts.add(u, x), val[u] += x;
                } else if (k == 1 && !(edge && u == root)) {
                    ts.set(u, x), val[u] = x;
                } else if (k == 2 && edge && u != root) {
                    const int p = bt.par[u];
                    if (rnd(rng, 0, 1)) {
                        ts.add_edge(p, u, x), val[u] += x;
                    } else {
                        ts.set_edge(u, p, x), val[u] = x;
                    }
                    TPB_CHECK(ts.get_edge(u, p) == val[u]);
                }
                const auto ps = bt.path(u, v);
                const int l = ps[bt.dep[u] - std::min(bt.dep[u], bt.dep[v]) + 0 <= 0 ? 0 : 0];
                (void)l;
                ll s = 0;
                int lc = u;
                for (int w : ps)
                    if (bt.dep[w] < bt.dep[lc]) lc = w;
                for (int w : ps)
                    if (!edge || w != lc) s += val[w];
                TPB_CHECK(ts.path(u, v) == s && ts.lca(u, v) == lc);
                ll rs = 0;
                for (int w : bt.path(root, u))
                    if (!edge || w != root) rs += val[w];
                TPB_CHECK(ts.root_path(u) == rs);
                ll ss = 0;
                for (int w = 0; w < n; ++w)
                    if (bt.in_sub(w, u) && (!edge || w != u)) ss += val[w];
                TPB_CHECK(ts.subtree(u) == ss && ts.get(u) == val[u]);
                TPB_CHECK(ts.depth(u) == bt.dep[u] && ts.parent(u) == bt.par[u] && ts.size() == n);
            }
        }
        // TreePathAdd
        {
            TreePathAdd<ll> ta(g, root, init, edge);
            std::vector<ll> val = init;
            for (int t = 0; t < 60; ++t) {
                const int u = (int)rnd(rng, 0, n - 1), v = (int)rnd(rng, 0, n - 1);
                const ll x = rnd(rng, -9, 9);
                if (rnd(rng, 0, 1)) {
                    ta.path_add(u, v, x);
                    const auto ps = bt.path(u, v);
                    int lc = u;
                    for (int w : ps)
                        if (bt.dep[w] < bt.dep[lc]) lc = w;
                    for (int w : ps)
                        if (!edge || w != lc) val[w] += x;
                }
                for (int w = 0; w < n; ++w) TPB_CHECK(ta.get(w) == val[w]);
                if (edge && u != root) TPB_CHECK(ta.get_edge(bt.par[u], u) == val[u]);
                ll ss = 0;
                for (int w = 0; w < n; ++w)
                    if (bt.in_sub(w, u) && (!edge || w != u)) ss += val[w];
                TPB_CHECK(ta.subtree_sum(u) == ss);
                TPB_CHECK(ta.lca(u, v) == TreePathSum<ll>(g, root).lca(u, v) && ta.depth(u) == bt.dep[u] && ta.parent(u) == bt.par[u] && ta.size() == n);
            }
        }
        ++cases;
    }
    // 辺の一覧から作る形
    {
        TreePathSum<ll> ts(4, {{0, 1}, {1, 2}, {1, 3}}, 0, {}, true);
        ts.set_edge(1, 2, 5), ts.set_edge(3, 1, 7), ts.add_edge(0, 1, 1);
        TPB_CHECK(ts.path(2, 3) == 12 && ts.path(0, 2) == 6 && ts.subtree(1) == 12);
        TreePathAdd<ll> ta(4, {{0, 1}, {1, 2}, {1, 3}});
        ta.path_add(2, 3, 1);
        TPB_CHECK(ta.get(0) == 0 && ta.get(1) == 1 && ta.subtree_sum(1) == 3);
    }
    std::cout << "[tree_path_bit] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[tree_path_bit] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 50000 : 500000, Q = N;
    const auto g = random_tree(rng, N, 0);  // ランダムな木（親を自分より前からランダムに選ぶ）
    std::vector<std::array<int, 3>> qs(Q);
    for (auto& q : qs) q = {(int)rnd(rng, 0, N - 1), (int)rnd(rng, 0, N - 1), (int)rnd(rng, -1000, 1000)};
    auto t0 = std::chrono::steady_clock::now();
    ll chk = 0;
    {
        TreePathSum<ll> ts(g, 0);
        for (auto& q : qs) ts.add(q[0], q[2]), chk += ts.path(q[0], q[1]);
    }
    const ll ms1 = elapsed_ms(t0);
    // 一般版の形：HLD でパスを O(log N) 本の区間に分けて、それぞれを BIT で（O(log² N)）
    t0 = std::chrono::steady_clock::now();
    ll chk2 = 0;
    {
        tpb_internal::RootedTree tr;
        tr.build(g, 0);
        tpb_internal::Fenwick<ll> bit;
        bit.init(std::vector<ll>(N, 0));
        auto range = [&](int l, int r) { return bit.prefix(r) - bit.prefix(l); };
        for (auto& q : qs) {
            bit.add(tr.pre[q[0]], q[2]);
            int u = q[0], v = q[1];
            ll s = 0;
            while (tr.head[u] != tr.head[v]) {
                if (tr.pre[tr.head[u]] < tr.pre[tr.head[v]]) std::swap(u, v);
                s += range(tr.pre[tr.head[u]], tr.pre[u] + 1);
                u = tr.par[tr.head[u]];
            }
            if (tr.pre[u] > tr.pre[v]) std::swap(u, v);
            s += range(tr.pre[u], tr.pre[v] + 1);
            chk2 += s;
        }
    }
    const ll ms2 = elapsed_ms(t0);
    TPB_CHECK(chk == chk2);
    std::cout << "[tree_path_bit] benchmark: N = " << N << " の木で add と path を " << Q << " 回ずつ  TreePathSum " << ms1 << " ms / 一般版（HLD + BIT）" << ms2
              << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace tpb_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace tpb_demo {

using ll = long long;
[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(const std::string& s) { return "\"" + s + "\""; }
[[maybe_unused]] static std::string to_s(const char* s) { return to_s(std::string(s)); }
[[maybe_unused]] static std::string to_s(char c) { return std::string("'") + c + "'"; }
[[maybe_unused]] static std::string to_s(double v) {
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}
template <class T>
static std::string to_s(const T& v);
template <class A, class B>
static std::string to_s(const std::pair<A, B>& p);
template <class T>
static std::string to_s(const std::vector<T>& v);
template <class T>
static std::string to_s(const std::optional<T>& v);
template <class T>
static std::string to_s(const T& v) {  // << で出せるもの（整数など）
    std::ostringstream os;
    os << v;
    return os.str();
}
template <class A, class B>
static std::string to_s(const std::pair<A, B>& p) {
    return "(" + to_s(p.first) + ", " + to_s(p.second) + ")";
}
template <class T>
static std::string to_s(const std::vector<T>& v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "]";
}
template <class T>
static std::string to_s(const std::optional<T>& v) {
    return v ? to_s(*v) : std::string("nullopt");
}
[[maybe_unused]] static std::string pad(const std::string& s, size_t w) {
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
[[maybe_unused]] static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
[[maybe_unused]] static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
[[maybe_unused]] static void stmt(const std::string& s) { code(s + ";"); }
template <class V>
[[maybe_unused]] static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 38) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 66) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) tpb_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    tpb_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】木のパスの和（tree_path_bit） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 演算は和・xor など引き算ができて順番によらないものだけ。1 点更新 + パスの和 か、パスに加算 + 1 点取得 のどちらか\n"
              << "■ できること: パスの和・パスへの加算が O(log N)（一般版の HLD は O(log² N)）。根からの和・部分木の和も\n"
              << "■ API 一覧\n"
              << "  TreePathSum<T> ts(g / (n, edges), root, init, edge): add / set / get / add_edge / set_edge / get_edge / path / root_path / subtree / lca / depth / parent / size\n"
              << "  TreePathAdd<T> ta(g / (n, edges), root, init, edge): path_add / get / get_edge / subtree_sum / lca / depth / parent / size\n";

    // 木の形: 0 の子が 1, 2、1 の子が 3, 4、2 の子が 5
    RUN(std::vector<std::vector<int>> g = {{1, 2}, {0, 3, 4}, {0, 5}, {1}, {1}, {2}});
    api("【1】 TreePathSum<T> ts(g, root, init)  path(u, v) / root_path(v) / subtree(v)", "T", "パス u〜v の頂点の和 / 根〜v の和 / v の部分木の和", "O(log N)");
    RUN(TreePathSum<ll> ts(g, 0, {1, 2, 3, 4, 5, 6}));
    SHOW(ts.path(3, 5), "3-1-0-2-5 → 4 + 2 + 1 + 3 + 6");
    SHOW(ts.root_path(4), "0-1-4 → 1 + 2 + 5");
    SHOW(ts.subtree(1), "{1, 3, 4} → 2 + 4 + 5");
    SHOW(ts.lca(3, 4), "");

    api("【2】 add(v, x) / set(v, x) / get(v)", "void / void / const T&", "頂点 v の値に足す / 書く / 読む", "O(log N) / O(log N) / O(1)");
    RUN(ts.add(0, 100); ts.set(5, 0));
    SHOW(ts.path(3, 5), "根に +100、頂点 5 を 0 に");
    SHOW(ts.get(0), "");
    SHOW(ts.depth(5), "");
    SHOW(ts.parent(5), "");

    api("【3】 edge = true（辺に値を載せる）  set_edge / add_edge / get_edge / path", "void / void / const T& / T", "辺 u–v の値を書く / 足す / 読む、パスの辺の和（= 距離）", "O(log N)");
    RUN(TreePathSum<ll> dist(6, {{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}}, 0, {}, true));
    RUN(dist.set_edge(0, 1, 3); dist.set_edge(0, 2, 1); dist.set_edge(1, 3, 4); dist.set_edge(4, 1, 2); dist.set_edge(2, 5, 5));
    SHOW(dist.path(3, 5), "距離 4 + 3 + 1 + 5");
    RUN(dist.add_edge(1, 0, 10));
    SHOW(dist.path(3, 5), "辺 0–1 を +10（ABC294 G の形：辺の重みの変更と距離）");
    SHOW(dist.get_edge(0, 1), "");
    SHOW(dist.subtree(1), "頂点 1 より下の辺 {1–3, 1–4} の和");

    api("【4】 TreePathAdd<T> ta(g, root, init)  path_add(u, v, x) / get(v) / subtree_sum(v)", "void / T / T",
        "パス u〜v の全頂点に x / 頂点 v の値 / v の部分木の値の和", "O(log N)");
    RUN(TreePathAdd<ll> ta(g, 0));
    RUN(ta.path_add(3, 5, 1); ta.path_add(4, 4, 10); ta.path_add(1, 2, 100));
    SHOW(ta.get(0), "3〜5 で 1、1〜2 で 100");
    SHOW(ta.get(1), "");
    SHOW(ta.get(4), "");
    SHOW(ta.subtree_sum(1), "{1, 3, 4} → 101 + 1 + 10");
    RUN(TreePathAdd<ll> te(g, 0, {}, true));
    RUN(te.path_add(3, 5, 1));
    SHOW(te.get_edge(0, 1), "edge = true：パス 3〜5 の全部の辺に 1");
    SHOW(te.get(0), "根には辺がないので 0");
}

}  // namespace tpb_demo

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
    tpb_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    tpb_test::benchmark();    // 速度計測
    tpb_demo::demo();         // 使い方デモ
    return 0;
}
