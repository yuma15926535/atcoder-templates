// =============================================================================
//  二乗の木 DP（部分木のサイズの分だけ配列をマージする木 DP。全体で O(N^2)）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  dp[v][k]（k = 部分木の中で選んだ頂点の数・連結成分の数など、0〜部分木のサイズ）を、子の配列と 2 重ループで
//  マージしていく木 DP。どちらの配列も「部分木のサイズ + 1」の長さなら、全体の計算量は O(N^2) になる
//  （2 頂点の組が、LCA でちょうど 1 回だけ掛け合わされるため）。k を K 以下に切ると O(NK)
//    例: ABC287 F（Components）… 連結成分がちょうど k 個になる頂点の選び方 → tree_count_components<mint>(g)
//        「部分木から k 個選ぶ」「サイズ k の連結な部分の価値の最大」「親を選んだときだけ子を選べるナップサック」
//
//  機能一覧（g は木の隣接リスト vector<vector<int>>、頂点は 0〜N−1）                       計算量
//   一般形
//    square_tree_dp(g, root, init, merge, finish)   帰りがけ順に DP を作り、根の DP を返す            O(N^2)（merge の書き方による）
//      init(v) → DP                    頂点 v だけの DP（例: vector<T>{選ばない, 選ぶ}）
//      merge(DP& dv, DP&& dc, v, c)    子 c の DP を v の DP に混ぜる。dv.size() × dc.size() の 2 重ループで書く
//      finish(v, DP& dv)               v の子を全部混ぜた後に呼ばれる（親への辺の処理・答えの記録。省略可）
//    conv_add(a, b)                    c[i + j] += a[i] × b[j]（長さ |a| + |b| − 1）                       O(|a||b|)
//    conv_max(a, b, ng) / conv_min     c[i + j] = max / min (a[i] + b[j])（ng は「無い」を表す値。ng は足さない）
//   よく使う形（そのまま答えになる）
//    tree_count_components<T>(g)       ans[k] = 誘導部分グラフの連結成分がちょうど k 個になる頂点集合の数（k = 0〜N、ans[0] = 1 は空集合）
//    tree_count_independent<T>(g)      ans[k] = サイズ k の独立集合（どの 2 頂点も隣り合わない）の数（k = 0〜N）
//    tree_max_connected(g, val)        best[k] = サイズ k の連結な頂点集合の val の和の最大（k = 1〜N。best[0] = 0）
//    tree_knapsack_closed(g, root, w, v, W)   「子を選ぶなら親も選ぶ」頂点集合で、重さの和 ≤ W のときの価値の和の最大
//                                      （帰りがけ順の DP で O(N W)。W は 0 以上）
//
//  注意
//    再帰を使わないので N = 2×10^5 の深い木（道）でも落ちない。ただし O(N^2) なので N ≤ 5000 くらいが目安
//    merge では dv を新しい配列に作り直してよい（dv = std::move(新しい配列)）。dc は使い終わったら捨ててよい
//    T は mint（数え上げ）・long long など。long long の数え上げはすぐあふれるので mod を取る型にする
//  速さの目安: N = 5000 で tree_count_components 約 0.05 秒、tree_max_connected 約 0.02 秒（道・ランダムな木とも）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    treesq_test::random_test()  … 頂点集合を全部調べる愚直解との突き合わせ
//    treesq_test::benchmark()    … 速度計測
//    treesq_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace treesq_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 根からの行きがけ順と親（木でなければ止める）
inline void preorder(const std::vector<std::vector<int>>& g, int root, std::vector<int>& ord, std::vector<int>& par, const char* f) {
    const int n = (int)g.size();
    if (root < 0 || root >= n) fail(f, "root は 0〜N−1");
    ord.clear(), par.assign(n, -2);
    par[root] = -1;
    std::vector<int> st{root};
    long long edges = 0;
    for (const auto& r : g) edges += (long long)r.size();
    if (edges != 2LL * (n - 1)) fail(f, "g は木（辺が N − 1 本の無向グラフ。各辺を両向きに入れる）にする");
    while (!st.empty()) {
        const int v = st.back();
        st.pop_back();
        ord.push_back(v);
        for (int c : g[v]) {
            if (c < 0 || c >= n) fail(f, "頂点番号が範囲外");
            if (c == par[v]) continue;
            if (par[c] != -2) fail(f, "閉路がある（木ではない）");
            par[c] = v;
            st.push_back(c);
        }
    }
    if ((int)ord.size() != n) fail(f, "連結でない（木ではない）");
}
struct NoFinish {
    template <class DP>
    void operator()(int, DP&) const {}
};
}  // namespace treesq_internal

// ---- 一般形 ----
template <class Init, class Merge, class Finish = treesq_internal::NoFinish>
auto square_tree_dp(const std::vector<std::vector<int>>& g, int root, Init init, Merge merge, Finish finish = Finish()) {
    using DP = decltype(init(0));
    std::vector<int> ord, par;
    treesq_internal::preorder(g, root, ord, par, "square_tree_dp");
    const int n = (int)g.size();
    std::vector<DP> dp(n);
    for (int i = n - 1; i >= 0; --i) {  // 帰りがけ順（子が先）
        const int v = ord[i];
        DP cur = init(v);
        for (int c : g[v]) {
            if (c == par[v]) continue;
            merge(cur, std::move(dp[c]), v, c);
            DP().swap(dp[c]);  // 子の DP を捨てる（メモリを O(N) に保つ）
        }
        finish(v, cur);
        dp[v] = std::move(cur);
    }
    return std::move(dp[root]);
}
template <class T>
std::vector<T> conv_add(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.empty() || b.empty()) return {};
    std::vector<T> c(a.size() + b.size() - 1, T(0));
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j) c[i + j] += a[i] * b[j];
    return c;
}
template <class T>
std::vector<T> conv_max(const std::vector<T>& a, const std::vector<T>& b, T ng) {
    if (a.empty() || b.empty()) return {};
    std::vector<T> c(a.size() + b.size() - 1, ng);
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == ng) continue;
        for (size_t j = 0; j < b.size(); ++j) {
            if (b[j] == ng) continue;
            const T s = a[i] + b[j];
            if (c[i + j] == ng || c[i + j] < s) c[i + j] = s;
        }
    }
    return c;
}
template <class T>
std::vector<T> conv_min(const std::vector<T>& a, const std::vector<T>& b, T ng) {
    if (a.empty() || b.empty()) return {};
    std::vector<T> c(a.size() + b.size() - 1, ng);
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == ng) continue;
        for (size_t j = 0; j < b.size(); ++j) {
            if (b[j] == ng) continue;
            const T s = a[i] + b[j];
            if (c[i + j] == ng || s < c[i + j]) c[i + j] = s;
        }
    }
    return c;
}

// ---- よく使う形 ----
// 連結成分がちょうど k 個になる頂点集合の数（ABC287 F）
template <class T>
std::vector<T> tree_count_components(const std::vector<std::vector<int>>& g) {
    const int n = (int)g.size();
    if (n == 0) return {T(1)};
    using A = std::array<T, 2>;  // [0] = v を選ばない, [1] = v を選ぶ
    auto init = [](int) { return std::vector<A>{A{T(1), T(0)}, A{T(0), T(1)}}; };
    auto merge = [](std::vector<A>& dv, std::vector<A>&& dc, int, int) {
        std::vector<A> nd(dv.size() + dc.size() - 1, A{T(0), T(0)});
        for (size_t i = 0; i < dv.size(); ++i)
            for (size_t j = 0; j < dc.size(); ++j) {
                const T c0 = dc[j][0], c1 = dc[j][1];
                nd[i + j][0] += dv[i][0] * (c0 + c1);
                nd[i + j][1] += dv[i][1] * c0;
                if (i + j >= 1) nd[i + j - 1][1] += dv[i][1] * c1;  // v と子 c の成分がつながる
            }
        dv = std::move(nd);
    };
    const auto r = square_tree_dp(g, 0, init, merge);
    std::vector<T> ans(n + 1, T(0));
    for (size_t k = 0; k < r.size() && k <= (size_t)n; ++k) ans[k] = r[k][0] + r[k][1];
    return ans;
}
// サイズ k の独立集合の数
template <class T>
std::vector<T> tree_count_independent(const std::vector<std::vector<int>>& g) {
    const int n = (int)g.size();
    if (n == 0) return {T(1)};
    using A = std::array<T, 2>;
    auto init = [](int) { return std::vector<A>{A{T(1), T(0)}, A{T(0), T(1)}}; };
    auto merge = [](std::vector<A>& dv, std::vector<A>&& dc, int, int) {
        std::vector<A> nd(dv.size() + dc.size() - 1, A{T(0), T(0)});
        for (size_t i = 0; i < dv.size(); ++i)
            for (size_t j = 0; j < dc.size(); ++j) {
                nd[i + j][0] += dv[i][0] * (dc[j][0] + dc[j][1]);
                nd[i + j][1] += dv[i][1] * dc[j][0];
            }
        dv = std::move(nd);
    };
    const auto r = square_tree_dp(g, 0, init, merge);
    std::vector<T> ans(n + 1, T(0));
    for (size_t k = 0; k < r.size() && k <= (size_t)n; ++k) ans[k] = r[k][0] + r[k][1];
    return ans;
}
// サイズ k の連結な頂点集合の val の和の最大（best[0] = 0、作れない k は無い。N ≥ 1）
template <class T>
std::vector<T> tree_max_connected(const std::vector<std::vector<int>>& g, const std::vector<T>& val) {
    const int n = (int)g.size();
    if ((int)val.size() != n) treesq_internal::fail("tree_max_connected", "val の長さは N");
    if (n == 0) return {T(0)};
    // dp[v][k] = v を一番上に含む、v の部分木の中の連結な集合（サイズ k）の和の最大。k = 0 は「使わない」
    std::vector<T> best(n + 1);
    std::vector<char> has(n + 1, 0);
    best[0] = T(0), has[0] = 1;
    using P = std::pair<char, T>;  // (ある, 値)
    auto init = [&](int v) { return std::vector<P>{P{0, T(0)}, P{1, val[v]}}; };
    auto merge = [](std::vector<P>& dv, std::vector<P>&& dc, int, int) {
        std::vector<P> nd(dv.size() + dc.size() - 1, P{0, T(0)});
        for (size_t i = 1; i < dv.size(); ++i) {
            if (!dv[i].first) continue;
            if (!nd[i].first || nd[i].second < dv[i].second) nd[i] = dv[i];  // 子を使わない
            for (size_t j = 1; j < dc.size(); ++j) {
                if (!dc[j].first) continue;
                const T s = dv[i].second + dc[j].second;
                if (!nd[i + j].first || nd[i + j].second < s) nd[i + j] = {1, s};
            }
        }
        dv = std::move(nd);
    };
    auto finish = [&](int, std::vector<P>& dv) {
        for (size_t k = 1; k < dv.size(); ++k)
            if (dv[k].first && (!has[k] || best[k] < dv[k].second)) best[k] = dv[k].second, has[k] = 1;
    };
    square_tree_dp(g, 0, init, merge, finish);
    return best;
}
// 「子を選ぶなら親も選ぶ」頂点集合で、重さ w の和 ≤ W のときの価値 v の和の最大（空集合 0 も含む）
template <class T>
T tree_knapsack_closed(const std::vector<std::vector<int>>& g, int root, const std::vector<int>& w, const std::vector<T>& v, int W) {
    const int n = (int)g.size();
    if ((int)w.size() != n || (int)v.size() != n) treesq_internal::fail("tree_knapsack_closed", "w, v の長さは N");
    if (W < 0) treesq_internal::fail("tree_knapsack_closed", "W は 0 以上");
    if (n == 0) return T(0);
    std::vector<int> ord, par;
    treesq_internal::preorder(g, root, ord, par, "tree_knapsack_closed");
    std::vector<int> pos(n), sz(n, 1);
    for (int i = 0; i < n; ++i) pos[ord[i]] = i;
    for (int i = n - 1; i > 0; --i) sz[par[ord[i]]] += sz[ord[i]];
    // dp[i][c] = 行きがけ順で i 番目以降を決めるときの価値の最大（重さ c まで使える）。v を選ばないなら部分木ごと飛ばす
    const size_t Wn = (size_t)W + 1;
    std::vector<T> dp((size_t)(n + 1) * Wn, T(0));
    for (int i = n - 1; i >= 0; --i) {
        const int x = ord[i];
        T* cur = &dp[(size_t)i * Wn];
        const T* nxt = &dp[(size_t)(i + 1) * Wn];
        const T* skip = &dp[(size_t)(i + sz[x]) * Wn];
        for (int c = 0; c <= W; ++c) {
            T best = skip[c];
            if (w[x] <= c) {
                const T t = nxt[c - w[x]] + v[x];
                if (best < t) best = t;
            }
            cur[c] = best;
        }
    }
    return dp[W];
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace treesq_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define TS_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[dp_tree_square] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

struct Mint {  // テストと速度計測用の mod 998244353
    static constexpr unsigned MOD = 998244353;
    unsigned x = 0;
    Mint() {}
    Mint(long long v) : x((unsigned)(((v % (long long)MOD) + MOD) % MOD)) {}
    Mint& operator+=(const Mint& o) {
        x += o.x;
        if (x >= MOD) x -= MOD;
        return *this;
    }
    Mint& operator-=(const Mint& o) {
        x += MOD - o.x;
        if (x >= MOD) x -= MOD;
        return *this;
    }
    Mint operator+(const Mint& o) const { return Mint(*this) += o; }
    Mint operator-(const Mint& o) const { return Mint(*this) -= o; }
    Mint operator*(const Mint& o) const {
        Mint r;
        r.x = (unsigned)((unsigned long long)x * o.x % MOD);
        return r;
    }
    bool operator==(const Mint& o) const { return x == o.x; }
};

static std::vector<std::vector<int>> random_tree(std::mt19937_64& rng, int n, int shape) {
    std::vector<std::vector<int>> g(n);
    std::vector<int> label(n);
    for (int i = 0; i < n; ++i) label[i] = i;
    std::shuffle(label.begin(), label.end(), rng);
    for (int i = 1; i < n; ++i) {
        int p;
        if (shape == 0) p = (int)rnd(rng, 0, i - 1);              // ランダム
        else if (shape == 1) p = i - 1;                            // 道
        else if (shape == 2) p = 0;                                // 星
        else p = (int)rnd(rng, std::max(0, i - 3), i - 1);         // 細長い
        g[label[i]].push_back(label[p]), g[label[p]].push_back(label[i]);
    }
    return g;
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 11);
        const auto g = random_tree(rng, n, (int)rnd(rng, 0, 3));
        std::vector<ll> val(n);
        for (auto& x : val) x = rnd(rng, -20, 20);
        std::vector<ll> comp(n + 1, 0), indep(n + 1, 0);
        std::vector<ll> best(n + 1, std::numeric_limits<ll>::min());
        best[0] = 0;
        for (int m = 0; m < (1 << n); ++m) {
            // 連結成分の数（Union-Find の代わりに BFS）
            std::vector<int> seen(n, 0);
            int cc = 0, cnt = 0;
            bool ind = true;
            ll s = 0;
            for (int v = 0; v < n; ++v) {
                if (!(m >> v & 1)) continue;
                ++cnt, s += val[v];
                for (int u : g[v])
                    if (m >> u & 1) ind = false;
                if (seen[v]) continue;
                ++cc;
                std::vector<int> st{v};
                seen[v] = 1;
                while (!st.empty()) {
                    const int x = st.back();
                    st.pop_back();
                    for (int u : g[x])
                        if ((m >> u & 1) && !seen[u]) seen[u] = 1, st.push_back(u);
                }
            }
            comp[cc]++;
            if (ind) indep[cnt]++;
            if (cc == 1) best[cnt] = std::max(best[cnt], s);
        }
        const auto rc = tree_count_components<ll>(g);
        const auto ri = tree_count_independent<ll>(g);
        const auto rm = tree_max_connected(g, val);
        TS_CHECK((int)rc.size() == n + 1 && (int)ri.size() == n + 1 && (int)rm.size() == n + 1);
        for (int k = 0; k <= n; ++k) TS_CHECK(rc[k] == comp[k] && ri[k] == indep[k] && rm[k] == best[k]);
        // mod の型でも同じ
        const auto rcm = tree_count_components<Mint>(g);
        for (int k = 0; k <= n; ++k) TS_CHECK(rcm[k] == Mint(comp[k]));
        // ナップサック（根をいろいろ）
        const int root = (int)rnd(rng, 0, n - 1);
        std::vector<int> par(n, -1), st{root}, vis(n, 0);
        vis[root] = 1;
        while (!st.empty()) {
            const int x = st.back();
            st.pop_back();
            for (int u : g[x])
                if (!vis[u]) vis[u] = 1, par[u] = x, st.push_back(u);
        }
        std::vector<int> w(n);
        for (auto& x : w) x = (int)rnd(rng, 0, 6);
        const int W = (int)rnd(rng, 0, 20);
        ll kb = 0;
        for (int m = 0; m < (1 << n); ++m) {
            bool ok = true;
            int sw = 0;
            ll sv = 0;
            for (int v = 0; v < n && ok; ++v)
                if (m >> v & 1) {
                    if (par[v] >= 0 && !(m >> par[v] & 1)) ok = false;
                    sw += w[v], sv += val[v];
                }
            if (ok && sw <= W) kb = std::max(kb, sv);
        }
        TS_CHECK(tree_knapsack_closed(g, root, w, val, W) == kb);
        // 一般形：部分木のサイズ（init = {1}、merge = 足す）と、conv_add で「部分木から k 個選ぶ数」= C(サイズ, k)
        {
            std::vector<ll> szs(n, 0);
            auto r = square_tree_dp(
                g, root, [](int) { return std::vector<ll>{1, 1}; },
                [](std::vector<ll>& a, std::vector<ll>&& b, int, int) { a = conv_add(a, b); },
                [&](int v, std::vector<ll>& a) { szs[v] = (ll)a.size() - 1; });
            TS_CHECK((int)r.size() == n + 1);
            ll binom = 1;
            for (int k = 0; k <= n; ++k) {
                TS_CHECK(r[k] == binom);
                binom = binom * (n - k) / (k + 1);
            }
            ll tot = 0;
            for (ll x : szs) tot += x;
            ll depth_sum = 0;  // Σ サイズ = Σ (深さ + 1)
            std::vector<int> dep(n, 0), st2{root}, vis2(n, 0);
            vis2[root] = 1;
            while (!st2.empty()) {
                const int x = st2.back();
                st2.pop_back();
                depth_sum += dep[x] + 1;
                for (int u : g[x])
                    if (!vis2[u]) vis2[u] = 1, dep[u] = dep[x] + 1, st2.push_back(u);
            }
            TS_CHECK(tot == depth_sum);
        }
        // conv_max / conv_min
        {
            const ll NG = std::numeric_limits<ll>::min();
            std::vector<ll> a(rnd(rng, 1, 5)), b(rnd(rng, 1, 5));
            for (auto& x : a) x = rnd(rng, 0, 3) ? rnd(rng, -9, 9) : NG;
            for (auto& x : b) x = rnd(rng, 0, 3) ? rnd(rng, -9, 9) : NG;
            const auto cmx = conv_max(a, b, NG), cmn = conv_min(a, b, NG);
            for (size_t k = 0; k < a.size() + b.size() - 1; ++k) {
                ll mx = NG, mn = NG;
                for (size_t i = 0; i < a.size(); ++i) {
                    if (k < i || k - i >= b.size() || a[i] == NG || b[k - i] == NG) continue;
                    const ll s = a[i] + b[k - i];
                    if (mx == NG || s > mx) mx = s;
                    if (mn == NG || s < mn) mn = s;
                }
                TS_CHECK(cmx[k] == mx && cmn[k] == mn);
            }
        }
        cases += 4 * (n + 1);
    }
    std::cout << "[dp_tree_square] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[dp_tree_square] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 1000 : 5000;
    std::string out;
    unsigned long long chk = 0;
    for (int shape = 0; shape < 2; ++shape) {
        const auto g = random_tree(rng, N, shape);
        std::vector<ll> val(N);
        for (auto& x : val) x = rnd(rng, -1000000000, 1000000000);
        auto t0 = std::chrono::steady_clock::now();
        const auto r = tree_count_components<Mint>(g);
        const ll t1 = elapsed_ms(t0);
        t0 = std::chrono::steady_clock::now();
        const auto b = tree_max_connected(g, val);
        const ll t2 = elapsed_ms(t0);
        chk += r[N / 3].x + (unsigned long long)b[N / 2];
        out += std::string(shape == 0 ? " ランダムな木" : " 道") + ": tree_count_components " + std::to_string(t1) + " ms / tree_max_connected " +
               std::to_string(t2) + " ms";
    }
    {
        // 深さ 2×10^5 の道でも落ちない（N が大きいので DP は軽いものだけ）
        const int M = TEST_SCALE == 1 ? 20000 : 200000;
        const auto g = random_tree(rng, M, 1);
        const auto t0 = std::chrono::steady_clock::now();
        std::vector<int> w(M, 1);
        std::vector<ll> v(M, 1);
        chk += (unsigned long long)tree_knapsack_closed(g, 0, w, v, 10);
        out += " / 道 N = " + std::to_string(M) + " で tree_knapsack_closed(W = 10) " + std::to_string(elapsed_ms(t0)) + " ms";
    }
    std::cout << "[dp_tree_square] benchmark (N = " << N << "):" << out << "  (check " << chk % 1000 << ")\n";
}

}  // namespace treesq_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace treesq_demo {

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
#define SHOW(expr, comment) treesq_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    treesq_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("二乗の木 DP（dp_tree_square） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  square_tree_dp(g, root, init, merge, finish) / conv_add / conv_max / conv_min\n"
              << "  tree_count_components<T> / tree_count_independent<T> / tree_max_connected / tree_knapsack_closed\n";

    RUN(std::vector<std::vector<int>> path = {{1}, {0, 2}, {1, 3}, {2}});  // 0 - 1 - 2 - 3
    RUN(std::vector<std::vector<int>> star = {{1, 2, 3}, {0}, {0}, {0}});  // 0 が真ん中

    api("【1】 tree_count_components<T>(g)", "vector<T>（長さ N + 1）",
        "ans[k] = 選んだ頂点の誘導部分グラフの連結成分がちょうど k 個になる選び方の数（ABC287 F）。ans[0] = 1 は空集合", "O(N^2)");
    SHOW(tree_count_components<long long>(path), "4 頂点の道：k = 1 は区間 10 個、k = 2 は 5 個");
    SHOW(tree_count_components<long long>(star), "星：真ん中を選ぶと 1 つにまとまる");
    note("答えを mod で求めるときは T = mint（atcoder::modint998244353 など）");

    api("【2】 tree_count_independent<T>(g)", "vector<T>", "ans[k] = サイズ k の独立集合の数", "O(N^2)");
    SHOW(tree_count_independent<long long>(path), "道：{0, 2}, {0, 3}, {1, 3} の 3 個がサイズ 2");
    SHOW(tree_count_independent<long long>(star), "星：葉 3 つから選ぶ（と真ん中 1 つ）");

    api("【3】 tree_max_connected(g, val)", "vector<T>", "best[k] = サイズ k の連結な頂点集合の val の和の最大（k = 1〜N）", "O(N^2)");
    RUN(const std::vector<long long> val = {5, -2, 7, 3});
    SHOW(tree_max_connected(path, val), "道：k = 2 は {2, 3} の 10、k = 3 は {0, 1, 2} の 10");
    SHOW(tree_max_connected(star, val), "星：k ≥ 2 は真ん中 0 を含む");

    api("【4】 tree_knapsack_closed(g, root, w, v, W)", "T",
        "「子を選ぶなら親も選ぶ」（根に近い方から選ぶ）頂点集合で、重さ w の和 ≤ W のときの価値 v の和の最大", "O(N W)");
    RUN(const std::vector<int> w = {2, 1, 3, 1});
    SHOW(tree_knapsack_closed(path, 0, w, val, 3), "根 0：{0} の 5（1 を入れると −2）");
    SHOW(tree_knapsack_closed(path, 0, w, val, 6), "{0, 1, 2} = 5 − 2 + 7 = 10");
    SHOW(tree_knapsack_closed(star, 0, w, val, 6), "{0, 2, 3}：重さ 6・価値 15");

    api("【5】 square_tree_dp(g, root, init, merge, finish)：一般形", "DP（根の DP）",
        "init(v) で 1 頂点の DP、merge(dv, dc, v, c) で子を混ぜる（2 重ループ）、finish(v, dv) は子を混ぜ終わった後", "O(N^2)");
    code("// 例: 部分木から k 個選ぶ方法の数（= C(部分木のサイズ, k)）。conv_add で配列を掛ける");
    RUN(auto r = square_tree_dp(star, 0, [](int) { return std::vector<long long>{1, 1}; }, [](std::vector<long long>& a, std::vector<long long>&& b, int, int) { a = conv_add(a, b); }));
    SHOW(r, "C(4, k)");
    code("// finish で頂点ごとの答えを記録できる（ここでは部分木のサイズ）");
    RUN(std::vector<int> sz(4));
    RUN(square_tree_dp(path, 0, [](int) { return std::vector<long long>{1, 1}; }, [](std::vector<long long>& a, std::vector<long long>&& b, int, int) { a = conv_add(a, b); }, [&](int v, std::vector<long long>& a) { sz[v] = (int)a.size() - 1; }));
    SHOW(sz, "根 0 からの部分木のサイズ");

    api("【6】 conv_add(a, b) / conv_max(a, b, ng) / conv_min(a, b, ng)", "vector<T>", "c[i + j] = Σ a[i] b[j] / max (a[i] + b[j]) / min（ng は「無い」）",
        "O(|a||b|)");
    SHOW(conv_add(std::vector<long long>{1, 2}, std::vector<long long>{1, 3, 1}), "");
    SHOW(conv_max(std::vector<long long>{0, 5}, std::vector<long long>{0, -1, 4}, -1LL), "−1 は「無い」");

    api("【7】 K 個までで切ると O(NK)", "—", "merge で新しい配列の長さを min(|dv| + |dc| − 1, K + 1) にする（「k ≤ K だけ知りたい」とき）", "O(NK)");
}

}  // namespace treesq_demo

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
    treesq_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    treesq_test::benchmark();    // 速度計測
    treesq_demo::demo();         // 使い方デモ
    return 0;
}
