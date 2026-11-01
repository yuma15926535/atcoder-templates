// =============================================================================
//  【制限付き最速版】木のいもす（パス・部分木・1 頂点に足す予約を O(1)（パスは + LCA）で、最後にまとめて各頂点の値にする）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : 足している途中で値を読めない（全部足してから build() で読む）。演算は加算（引き算ができて順番によらないもの）だけ。
//                           パスの和・最小などの質問はない
//    できること           : 「パス u〜v の全頂点（全辺）に x」「v の部分木の全頂点に x」「頂点 v に x」を何回でも予約して、
//                           build() で全頂点（全辺）の最終的な値を出す（build は何度でも呼べる）
//    一般版との比較       : graph_grid_v2/05_hld.cpp の TreeLazyQuery はパスへの加算 O(log² N)、
//                           EulerTourSubtreeTemplate.txt の部分木への加算は O(log N)
//                           → この版はパスへの加算 O(1) + LCA 1 回（O(log N)。実際はほぼ定数）、部分木・1 頂点は O(1)、build O(N)
//    制限に合わないとき   : 足すのと読むのがまざるなら tree_path_bit.cpp の TreePathAdd（O(log N)）か graph_grid_v2/05_hld.cpp
//
//  しくみ: パス u〜v に x は「u と v に +x、lca と lca の親に −x」を置き、最後に各頂点で部分木の中の合計を取る（葉から根へ）。
//          部分木 v に x は「v に +x」を置き、最後に根から v までの合計を取る（根から葉へ）
//    例: 「パス u〜v の頂点に 1 ずつ足す」を Q 回 → 各頂点が何回通られたか、「部分木に x を足す」を Q 回（ABC138 D）
//
//  機能一覧（頂点は 0〜N−1。g は vector<vector<int>> の隣接リスト、または (n, 辺の一覧)。木でなければ止める）     計算量
//    TreeImos<T> ti(g, root = 0, init = {}, edge = false)   init[v] = 頂点 v の初期値（edge なら v と親の辺の値）   O(N)
//    TreeImos<T> ti(n, edges, root = 0, init = {}, edge = false)
//      ti.path_add(u, v, x)              パス u〜v の全頂点（edge なら全辺）に x                                O(1) + LCA
//      ti.subtree_add(v, x)              v の部分木の全頂点（edge なら v より下の全辺）に x                        O(1)
//      ti.vertex_add(v, x)               頂点 v（edge なら v と親の辺）に x                                       O(1)
//      ti.build()                        全頂点の値を作って返す（予約は残る。続けて足して何度でも呼べる）          O(N)
//      ti.get(v) / ti.get_edge(u, v) / ti.to_vector()   build 後：頂点 v の値 / 辺 u–v の値（edge のとき）/ 全体  O(1)
//      ti.clear()                        予約を全部消す（初期値に戻る）                                            O(N)
//      ti.lca(u, v)（O(log N)）/ ti.depth(v) / ti.parent(v)（根は −1）/ ti.size()
//    edge = true のとき値は「v と親を結ぶ辺」の値を v に持つ（根の値は使わない・0）
//
//  T は long long（ふつう）/ modint など（+, − ができる型）
//  速さの目安: N = 10^6 のランダムな木で path_add 10^6 回 + build 約 0.26 秒（木を作るのに別に 約 0.15 秒）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    tim_test::random_test()  … パスの頂点を 1 つずつ足す愚直解との突き合わせ
//    tim_test::benchmark()    … 速度計測
//    tim_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
namespace tim_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 根付き木（非再帰）。order は根から近い順、pre は HLD の行きがけ順。LCA は重軽分解
struct RootedTree {
    int n = 0, root = 0;
    std::vector<int> par, dep, head, pre, order;
    void build(const std::vector<std::vector<int>>& g, int r) {
        n = (int)g.size(), root = r;
        if (r < 0 || r >= std::max(n, 1)) fail("tree_imos", "根が範囲外");
        par.assign(n, -1), dep.assign(n, 0), head.assign(n, 0), pre.assign(n, -1), order.clear();
        if (n == 0) return;
        order.reserve(n);
        std::vector<char> seen(n, 0);
        order.push_back(r), seen[r] = 1;
        for (int k = 0; k < (int)order.size(); ++k) {
            const int v = order[k];
            bool skipped_parent = false;
            for (int u : g[v]) {
                if (u < 0 || u >= n) fail("tree_imos", "頂点の番号が範囲外");
                if (u == par[v] && !skipped_parent) {
                    skipped_parent = true;
                    continue;
                }
                if (seen[u]) fail("tree_imos", "木ではない（閉路・多重辺がある）");
                seen[u] = 1, par[u] = v, dep[u] = dep[v] + 1, order.push_back(u);
            }
        }
        if ((int)order.size() != n) fail("tree_imos", "木ではない（つながっていない）");
        std::vector<int> sz(n, 1), heavy(n, -1);
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
    int child_of_edge(int u, int v) const {
        if (par[u] == v) return u;
        if (par[v] == u) return v;
        fail("tree_imos", "u と v は隣どうしではない");
    }
    void check(int v) const {
        if (v < 0 || v >= n) fail("tree_imos", "頂点が範囲外");
    }
};
inline std::vector<std::vector<int>> to_adj(int n, const std::vector<std::pair<int, int>>& edges) {
    std::vector<std::vector<int>> g(n);
    for (auto [u, v] : edges) {
        if (u < 0 || u >= n || v < 0 || v >= n) fail("tree_imos", "辺の端が範囲外");
        g[u].push_back(v), g[v].push_back(u);
    }
    return g;
}
}  // namespace tim_internal

template <class T = long long>
struct TreeImos {
    TreeImos(const std::vector<std::vector<int>>& g, int root = 0, const std::vector<T>& init = {}, bool edge = false) : edge_(edge) { build_tree(g, root, init); }
    TreeImos(int n, const std::vector<std::pair<int, int>>& edges, int root = 0, const std::vector<T>& init = {}, bool edge = false) : edge_(edge) {
        build_tree(tim_internal::to_adj(n, edges), root, init);
    }
    void path_add(int u, int v, const T& x) {  // パス u〜v の全頂点（edge なら全辺）に x
        const int l = lca(u, v);
        up_[u] += x, up_[v] += x;
        if (edge_) {
            up_[l] -= x, up_[l] -= x;
        } else {
            up_[l] -= x;
            if (tr_.par[l] != -1) up_[tr_.par[l]] -= x;
        }
    }
    void subtree_add(int v, const T& x) {  // v の部分木の全頂点（edge なら v より下の全辺）に x
        tr_.check(v);
        down_[v] += x;
        if (edge_) pt_[v] -= x;  // v と親の辺は部分木に入らない
    }
    void vertex_add(int v, const T& x) {  // 頂点 v（edge なら v と親の辺）に x
        tr_.check(v);
        if (edge_ && v == tr_.root) tim_internal::fail("TreeImos::vertex_add", "edge = true のとき根には辺がない");
        pt_[v] += x;
    }
    const std::vector<T>& build() {
        const int n = tr_.n;
        val_.assign(n, T(0));
        std::vector<T> s = up_;
        for (int k = n - 1; k >= 1; --k) {  // 葉から根へ：部分木の中の up_ の合計
            const int v = tr_.order[k];
            s[tr_.par[v]] += s[v];
        }
        std::vector<T> d = down_;
        for (int k = 1; k < n; ++k) {  // 根から葉へ：根から v までの down_ の合計
            const int v = tr_.order[k];
            d[v] += d[tr_.par[v]];
        }
        for (int v = 0; v < n; ++v) val_[v] = init_[v] + pt_[v] + s[v] + d[v];
        if (edge_ && n > 0) val_[tr_.root] = T(0);
        built_ = true;
        return val_;
    }
    const T& get(int v) const {
        check_built("TreeImos::get");
        tr_.check(v);
        return val_[v];
    }
    const T& get_edge(int u, int v) const {
        if (!edge_) tim_internal::fail("TreeImos::get_edge", "edge = true で作ったときだけ使える");
        check_built("TreeImos::get_edge");
        tr_.check(u), tr_.check(v);
        return val_[tr_.child_of_edge(u, v)];
    }
    const std::vector<T>& to_vector() const {
        check_built("TreeImos::to_vector");
        return val_;
    }
    void clear() {
        std::fill(up_.begin(), up_.end(), T(0)), std::fill(down_.begin(), down_.end(), T(0)), std::fill(pt_.begin(), pt_.end(), T(0));
        built_ = false;
    }
    int lca(int u, int v) const {
        tr_.check(u), tr_.check(v);
        return tr_.lca(u, v);
    }
    int depth(int v) const { return tr_.check(v), tr_.dep[v]; }
    int parent(int v) const { return tr_.check(v), tr_.par[v]; }
    int size() const { return tr_.n; }

   private:
    bool edge_, built_ = false;
    tim_internal::RootedTree tr_;
    std::vector<T> init_, up_, down_, pt_, val_;  // up_: パス用（部分木で合計）、down_: 部分木用（根から合計）、pt_: 1 頂点
    void build_tree(const std::vector<std::vector<int>>& g, int root, const std::vector<T>& init) {
        const int n = (int)g.size();
        if (!init.empty() && (int)init.size() != n) tim_internal::fail("TreeImos", "init の長さが頂点数と違う");
        tr_.build(g, root);
        init_ = init.empty() ? std::vector<T>(n, T(0)) : init;
        if (edge_ && n > 0) init_[root] = T(0);
        up_.assign(n, T(0)), down_.assign(n, T(0)), pt_.assign(n, T(0));
    }
    void check_built(const char* f) const {
        if (!built_) tim_internal::fail(f, "build() の前に呼んだ");
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace tim_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define TIM_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[tree_imos] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

struct BruteTree {
    int n;
    std::vector<int> par, dep;
    BruteTree(const std::vector<std::vector<int>>& g, int root) : n((int)g.size()), par(n, -1), dep(n, 0) {
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
    std::vector<int> path(int u, int v) const {
        std::vector<int> a;
        while (dep[u] > dep[v]) a.push_back(u), u = par[u];
        while (dep[v] > dep[u]) a.push_back(v), v = par[v];
        while (u != v) a.push_back(u), a.push_back(v), u = par[u], v = par[v];
        a.push_back(u);  // 最後が lca
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
        const int p = shape == 0 ? (int)rnd(rng, 0, i - 1) : (shape == 1 ? i - 1 : (i - 1) / 2);
        g[id[i]].push_back(id[p]), g[id[p]].push_back(id[i]);
    }
    for (auto& a : g) std::shuffle(a.begin(), a.end(), rng);
    return g;
}

static void random_test() {
    std::mt19937_64 rng(20261101);
    long long cases = 0;
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 25);
        const auto g = random_tree(rng, n, (int)rnd(rng, 0, 2));
        const int root = (int)rnd(rng, 0, n - 1);
        const bool edge = rnd(rng, 0, 1);
        std::vector<ll> init(n);
        for (auto& v : init) v = rnd(rng, -9, 9);
        TreeImos<ll> ti(g, root, init, edge);
        BruteTree bt(g, root);
        std::vector<ll> val = init;
        if (edge) val[root] = 0;
        for (int t = 0; t < 40; ++t) {
            const int k = (int)rnd(rng, 0, 2), u = (int)rnd(rng, 0, n - 1), v = (int)rnd(rng, 0, n - 1);
            const ll x = rnd(rng, -9, 9);
            if (k == 0) {
                ti.path_add(u, v, x);
                const auto ps = bt.path(u, v);
                for (size_t i = 0; i < ps.size(); ++i)
                    if (!edge || i + 1 != ps.size()) val[ps[i]] += x;
            } else if (k == 1) {
                ti.subtree_add(u, x);
                for (int w = 0; w < n; ++w)
                    if (bt.in_sub(w, u) && (!edge || w != u)) val[w] += x;
            } else if (!(edge && u == root)) {
                ti.vertex_add(u, x), val[u] += x;
            }
            if (rnd(rng, 0, 3) == 0) TIM_CHECK(ti.build() == val);
        }
        TIM_CHECK(ti.build() == val && ti.to_vector() == val);
        for (int v = 0; v < n; ++v) {
            TIM_CHECK(ti.get(v) == val[v] && ti.depth(v) == bt.dep[v] && ti.parent(v) == bt.par[v]);
            if (edge && v != root) TIM_CHECK(ti.get_edge(bt.par[v], v) == val[v] && ti.get_edge(v, bt.par[v]) == val[v]);
        }
        const int a = (int)rnd(rng, 0, n - 1), b = (int)rnd(rng, 0, n - 1);
        TIM_CHECK(ti.lca(a, b) == bt.path(a, b).back() && ti.size() == n);
        ti.clear();
        std::vector<ll> z = init;
        if (edge) z[root] = 0;
        TIM_CHECK(ti.build() == z);
        ++cases;
    }
    {
        TreeImos<ll> ti(4, {{0, 1}, {1, 2}, {1, 3}});
        ti.path_add(2, 3, 1), ti.subtree_add(1, 10);
        TIM_CHECK(ti.build() == std::vector<ll>({0, 11, 11, 11}));
    }
    std::cout << "[tree_imos] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[tree_imos] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 100000 : 1000000, Q = N;
    const auto g = random_tree(rng, N, 0);
    std::vector<std::array<int, 3>> qs(Q);
    for (auto& q : qs) q = {(int)rnd(rng, 0, N - 1), (int)rnd(rng, 0, N - 1), (int)rnd(rng, -1000, 1000)};
    auto t0 = std::chrono::steady_clock::now();
    TreeImos<ll> ti(g, 0);
    const ll ms0 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (auto& q : qs) ti.path_add(q[0], q[1], q[2]);
    const auto& v = ti.build();
    ll chk = 0;
    for (int i = 0; i < N; i += 7) chk += v[i];
    std::cout << "[tree_imos] benchmark: N = " << N << " のランダムな木で 作る " << ms0 << " ms、path_add " << Q << " 回 + build " << elapsed_ms(t0) << " ms  (check "
              << chk % 1000 << ")\n";
}

}  // namespace tim_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace tim_demo {

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
#define SHOW(expr, comment) tim_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    tim_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】木のいもす（tree_imos） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 足している途中で値を読めない（全部足してから build で読む）・加算だけ\n"
              << "■ できること: パス・部分木・1 頂点への加算を O(1)（パスは + LCA）で予約 → 全頂点（全辺）の値を O(N)\n"
              << "■ API 一覧\n"
              << "  TreeImos<T> ti(g / (n, edges), root, init, edge): path_add / subtree_add / vertex_add / build / get / get_edge / to_vector / clear / lca / depth / parent / size\n";

    // 木の形: 0 の子が 1, 2、1 の子が 3, 4、2 の子が 5
    RUN(std::vector<std::vector<int>> g = {{1, 2}, {0, 3, 4}, {0, 5}, {1}, {1}, {2}});
    api("【1】 TreeImos<T> ti(g, root)  path_add(u, v, x) → build()", "const vector<T>&", "パス u〜v の全頂点に x を足す予約を反映した、全頂点の値", "path_add O(1) + LCA、build O(N)");
    RUN(TreeImos<ll> ti(g, 0));
    RUN(ti.path_add(3, 5, 1); ti.path_add(4, 2, 10));
    SHOW(ti.build(), "3-1-0-2-5 に 1、4-1-0-2 に 10");
    SHOW(ti.get(1), "");

    api("【2】 subtree_add(v, x) / vertex_add(v, x)", "void", "v の部分木の全頂点に x / 頂点 v だけに x", "O(1)");
    RUN(ti.subtree_add(1, 100); ti.vertex_add(5, 1000));
    SHOW(ti.build(), "続けて足して build し直せる（{1, 3, 4} に 100、5 に 1000）");
    SHOW(ti.lca(3, 5), "");
    SHOW(ti.depth(5), "");
    SHOW(ti.parent(3), "");
    RUN(ti.clear());
    SHOW(ti.build(), "clear で予約を全部消す");

    api("【3】 edge = true（辺に値を載せる）  path_add / subtree_add / get_edge", "void / void / const T&", "パスの全辺に x / v より下の全辺に x / 辺 u–v の値",
        "O(1) + LCA / O(1) / O(1)");
    RUN(TreeImos<ll> te(6, {{0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}}, 0, {}, true));
    RUN(te.path_add(3, 4, 1); te.path_add(3, 5, 10); te.subtree_add(1, 100));
    RUN(te.build());
    SHOW(te.get_edge(1, 3), "3–4 と 3–5 の両方で通る + 部分木 1 の辺");
    SHOW(te.get_edge(0, 1), "3–5 だけ通る（辺 0–1 は部分木 1 に入らない）");
    SHOW(te.get_edge(5, 2), "");
    SHOW(te.to_vector(), "v と親の辺の値を v に持つ（根は 0）");
}

}  // namespace tim_demo

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
    tim_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    tim_test::benchmark();    // 速度計測
    tim_demo::demo();         // 使い方デモ
    return 0;
}
