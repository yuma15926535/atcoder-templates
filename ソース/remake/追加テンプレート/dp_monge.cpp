// =============================================================================
//  Monge・決定単調性の DP（monotone minima・分割統治 DP・1 次元の DP・Alien DP（ラグランジュ緩和））   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  コスト c(i, j) が Monge（i < i'、j < j' で c(i, j) + c(i', j') ≤ c(i, j') + c(i', j)）なら、
//  「行ごとの最小を取る列」が行とともに右へずれる（決定単調性）ので、全部の組を見なくてよい。
//  例: c(i, j) = (区間 [i, j) の和)^2（値が 0 以上）、c(i, j) = f(j − i)（f が下に凸）
//  Alien DP は「ちょうど K 個選ぶ」の答え f(K) が K について下に凸のとき、1 個ごとにペナルティ λ を足して
//  個数の制約を外し、λ を二分探索する方法
//    例: ABC218 H（Red and Blue Lamps）… ちょうど R 個選ぶ → alien_dp
//
//  機能一覧（T は long long など。コストは「最小化」。最大化は −1 倍する）                       計算量
//    monotone_minima(H, W, f)           各行 i で f(i, j) を最小にする列 j（同じなら左）。f が Monge なら正しい   O((H + W) log H)
//    partition_monge(N, K, c)           [0, N) を空でない K 個の区間 [i, j) に分けるコスト c(i, j) の和の最小      O(K N log N)
//                                       （res[k] = ちょうど k 個に分けるときの最小。k = 0〜K）
//    monge_1d_dp(N, c)                  dp[0] = 0、dp[j] = min_{i<j} dp[i] + c(i, j)（j = 1〜N）を全部            O(N log N)
//    alien_dp(K, lo, hi, solve)         solve(λ) = {min (コスト + λ × 個数), その最小を取る中で一番少ない個数} を
//                                       返す関数。ちょうど K 個のときの最小コストを返す（λ は [lo, hi] の整数を二分探索）  O(log(hi − lo) × solve)
//    is_monge(H, W, f)                  f が Monge か（確かめる用）                                               O(HW)
//
//  注意
//    Monge でないコストに使うと、答えが正しくない（止まらない）。迷ったら小さい入力で愚直解と比べる
//    alien_dp：f(K) が下に凸、λ が整数で足りる（コストが整数）、solve(hi) の個数 ≤ K（hi を十分大きく）、
//              lo は −max_k (f(k − 1) − f(k)) より小さく（迷ったら lo = −(コストの絶対値の合計 + 1)、hi = −lo）。
//              solve は「同じ値なら個数の少ない方」を選ぶこと
//  速さの目安: monotone_minima（2×10^5 × 2×10^5）約 0.01 秒、partition_monge（N = 10^5、K = 10）約 0.03 秒、
//             monge_1d_dp（N = 2×10^5）約 0.02 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    monge_test::random_test()  … O(N^2)・O(KN^2) の愚直解との突き合わせ
//    monge_test::benchmark()    … 速度計測
//    monge_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
namespace monge_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace monge_internal

template <class F>
std::vector<int> monotone_minima(int H, int W, F f) {
    std::vector<int> res(H, 0);
    if (H <= 0 || W <= 0) return res;
    // 行 [r1, r2) の答えは列 [c1, c2] の中にある
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> st{{{0, H}, {0, W - 1}}};
    while (!st.empty()) {
        auto [rr, cc] = st.back();
        st.pop_back();
        const int r1 = rr.first, r2 = rr.second, c1 = cc.first, c2 = cc.second;
        if (r1 >= r2) continue;
        const int m = (r1 + r2) / 2;
        int best = c1;
        auto bv = f(m, c1);
        for (int j = c1 + 1; j <= c2; ++j) {
            auto v = f(m, j);
            if (v < bv) bv = v, best = j;
        }
        res[m] = best;
        st.push_back({{r1, m}, {c1, best}});
        st.push_back({{m + 1, r2}, {best, c2}});
    }
    return res;
}

template <class C>
auto partition_monge(int N, int K, C c) {
    using T = decltype(c(0, 1));
    if (N < 0 || K < 0) monge_internal::fail("partition_monge", "N, K は 0 以上");
    const T INF = std::numeric_limits<T>::max() / 4;
    std::vector<T> res(K + 1, INF);
    std::vector<T> dp(N + 1, INF), nd(N + 1);
    dp[0] = T(0);
    res[0] = N == 0 ? T(0) : INF;
    for (int k = 1; k <= K; ++k) {
        std::fill(nd.begin(), nd.end(), INF);
        // nd[j] = min_{i < j} dp[i] + c(i, j)（j = k〜N）。最小の i は j とともに増える → 分割統治
        std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> st{{{k, N + 1}, {k - 1, N - 1}}};
        while (!st.empty()) {
            auto [jr, ir] = st.back();
            st.pop_back();
            const int j1 = jr.first, j2 = jr.second, i1 = ir.first, i2 = ir.second;
            if (j1 >= j2) continue;
            const int m = (j1 + j2) / 2;
            int bi = -1;
            T bv = INF;
            for (int i = i1; i <= std::min(i2, m - 1); ++i) {
                if (!(dp[i] < INF)) continue;
                const T v = dp[i] + c(i, m);
                if (bi < 0 || v < bv) bv = v, bi = i;
            }
            nd[m] = bv;
            const int piv = bi < 0 ? i1 : bi;
            st.push_back({{j1, m}, {i1, piv}});
            st.push_back({{m + 1, j2}, {piv, i2}});
        }
        std::swap(dp, nd);
        res[k] = dp[N];
    }
    return res;
}

template <class C>
auto monge_1d_dp(int N, C c) {
    using T = decltype(c(0, 1));
    if (N < 0) monge_internal::fail("monge_1d_dp", "N は 0 以上");
    std::vector<T> dp(N + 1);
    dp[0] = T(0);
    if (N == 0) return dp;
    auto val = [&](int i, int j) { return dp[i] + c(i, j); };
    // 候補 i が一番よい j の範囲を [start, 次の start) として deque（配列 + 先頭位置）で持つ
    std::vector<int> cand, start;
    size_t head = 0;
    cand.push_back(0), start.push_back(1);
    for (int j = 1; j <= N; ++j) {
        while (head + 1 < cand.size() && start[head + 1] <= j) ++head;
        dp[j] = val(cand[head], j);
        if (j == N) break;
        // 候補 j を後ろに足す：後ろの候補より j が良くなる最初の位置を二分探索
        while (head < cand.size()) {
            const int s = std::max(start.back(), j + 1);
            if (val(j, s) <= val(cand.back(), s)) {
                if (cand.size() - head == 1) {  // 全部の範囲で j が良い
                    cand.back() = j, start.back() = j + 1;
                    break;
                }
                cand.pop_back(), start.pop_back();
                continue;
            }
            int lo = s, hi = N + 1;  // val(j, lo) > val(back, lo)、hi は j が良い最初（無ければ N + 1）
            while (hi - lo > 1) {
                const int mid = (lo + hi) / 2;
                if (val(j, mid) <= val(cand.back(), mid)) hi = mid;
                else lo = mid;
            }
            if (hi <= N) cand.push_back(j), start.push_back(hi);
            break;
        }
    }
    return dp;
}

// solve(λ) → {min (コスト + λ × 個数), そのときの最小の個数}
template <class Solve>
long long alien_dp(long long K, long long lo, long long hi, Solve solve) {
    if (lo > hi) monge_internal::fail("alien_dp", "lo ≤ hi にする");
    // 個数 ≤ K となる最小の λ を探す
    if (solve(hi).second > K) monge_internal::fail("alien_dp", "solve(hi) の個数が K より多い（hi を大きくする）");
    long long l = lo - 1, h = hi;  // solve(h) の個数 ≤ K
    while (h - l > 1) {
        const long long m = l + (h - l) / 2;
        if (solve(m).second <= K) h = m;
        else l = m;
    }
    const auto r = solve(h);
    return r.first - h * K;
}

template <class F>
bool is_monge(int H, int W, F f) {
    for (int i = 0; i + 1 < H; ++i)
        for (int j = 0; j + 1 < W; ++j)
            if (f(i, j) + f(i + 1, j + 1) > f(i, j + 1) + f(i + 1, j)) return false;
    return true;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace monge_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define MG_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[dp_monge] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
static const ll INF = std::numeric_limits<ll>::max() / 4;

// ちょうど k 個の「隣り合わない」要素を選ぶ和の最大（下に凸 = 最小化で −1 倍）の愚直 O(NK)
static std::vector<ll> brute_nonadj(const std::vector<ll>& a, int K) {
    const int n = (int)a.size();
    const ll NEG = std::numeric_limits<ll>::min() / 4;
    // dp[i][k][t]：i 番目まで、k 個、最後を選んだか
    std::vector<std::vector<std::array<ll, 2>>> dp(n + 1, std::vector<std::array<ll, 2>>(K + 1, {NEG, NEG}));
    dp[0][0][0] = 0;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k <= K; ++k)
            for (int t = 0; t < 2; ++t) {
                if (dp[i][k][t] == NEG) continue;
                dp[i + 1][k][0] = std::max(dp[i + 1][k][0], dp[i][k][t]);
                if (!t && k + 1 <= K) dp[i + 1][k + 1][1] = std::max(dp[i + 1][k + 1][1], dp[i][k][t] + a[i]);
            }
    std::vector<ll> r(K + 1);
    for (int k = 0; k <= K; ++k) r[k] = std::max(dp[n][k][0], dp[n][k][1]);
    return r;
}
// 隣り合わない要素を選ぶ：min Σ(−a_i) + λ × 個数、同じなら個数の少ない方
static std::pair<ll, ll> solve_nonadj(const std::vector<ll>& a, ll lam) {
    std::pair<ll, ll> no{0, 0}, yes{INF, 0};  // 最後を選ばない / 選んだ
    for (ll x : a) {
        const std::pair<ll, ll> ny{no.first - x + lam, no.second + 1};
        const std::pair<ll, ll> nn = std::min(no, yes);
        no = nn, yes = ny;
    }
    return std::min(no, yes);
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 14);
        std::vector<ll> a(n), pre(n + 1, 0);
        for (auto& x : a) x = rnd(rng, 0, 20);
        for (int i = 0; i < n; ++i) pre[i + 1] = pre[i] + a[i];
        const ll C = rnd(rng, 0, 30);
        auto c = [&](int i, int j) { return (pre[j] - pre[i]) * (pre[j] - pre[i]) + C; };
        MG_CHECK(is_monge(n + 1, n + 1, [&](int i, int j) { return i < j ? c(i, j) : (ll)1e12 * (i - j + 1); }));
        // monotone_minima：(x_i − y_j)^2 + g(i) + h(j)（x, y が増える）は Monge
        {
            const int H = (int)rnd(rng, 0, 12), W = (int)rnd(rng, 0, 12);
            std::vector<ll> x(H), y(W), g(H), h(W);
            for (auto& v : x) v = rnd(rng, -20, 20);
            for (auto& v : y) v = rnd(rng, -20, 20);
            std::sort(x.begin(), x.end()), std::sort(y.begin(), y.end());
            for (auto& v : g) v = rnd(rng, -50, 50);
            for (auto& v : h) v = rnd(rng, -50, 50);
            auto f = [&](int i, int j) { return (x[i] - y[j]) * (x[i] - y[j]) + g[i] + h[j]; };
            MG_CHECK(is_monge(H, W, f));
            const auto r = monotone_minima(H, W, f);
            for (int i = 0; i < H && W > 0; ++i) {
                int b = 0;
                for (int j = 1; j < W; ++j)
                    if (f(i, j) < f(i, b)) b = j;
                MG_CHECK(r[i] == b);
            }
            ++cases;
        }
        // partition_monge と O(K N^2)
        {
            const int K = (int)rnd(rng, 0, 6);
            const auto r = partition_monge(n, K, c);
            std::vector<ll> dp(n + 1, INF);
            dp[0] = 0;
            MG_CHECK(r[0] == (n == 0 ? 0 : INF));
            for (int k = 1; k <= K; ++k) {
                std::vector<ll> nd(n + 1, INF);
                for (int j = 1; j <= n; ++j)
                    for (int i = 0; i < j; ++i)
                        if (dp[i] < INF) nd[j] = std::min(nd[j], dp[i] + c(i, j));
                dp = nd;
                MG_CHECK(r[k] == dp[n]);
            }
            ++cases;
        }
        // monge_1d_dp と O(N^2)
        {
            const auto r = monge_1d_dp(n, c);
            std::vector<ll> dp(n + 1, INF);
            dp[0] = 0;
            for (int j = 1; j <= n; ++j)
                for (int i = 0; i < j; ++i) dp[j] = std::min(dp[j], dp[i] + c(i, j));
            MG_CHECK(r == dp);
            // f(j − i) が下に凸の別のコスト
            auto c2 = [&](int i, int j) { return (ll)(j - i) * (j - i) * 3 - 7 * (j - i); };
            const auto r2 = monge_1d_dp(n, c2);
            std::vector<ll> d2(n + 1, INF);
            d2[0] = 0;
            for (int j = 1; j <= n; ++j)
                for (int i = 0; i < j; ++i) d2[j] = std::min(d2[j], d2[i] + c2(i, j));
            MG_CHECK(r2 == d2);
            ++cases;
        }
        // alien_dp：ちょうど K 個の隣り合わない要素の和の最大
        {
            std::vector<ll> b(n);
            for (auto& x : b) x = rnd(rng, -10, 30);
            const int maxk = (n + 1) / 2;
            const auto br = brute_nonadj(b, maxk);
            for (int K = 0; K <= maxk; ++K) {
                const ll got = alien_dp(K, -10000, 10000, [&](ll lam) { return solve_nonadj(b, lam); });
                MG_CHECK(-got == br[K]);
                ++cases;
            }
        }
    }
    std::cout << "[dp_monge] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[dp_monge] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000, N2 = TEST_SCALE == 1 ? 10000 : 100000;
    std::vector<ll> a(N), pre(N + 1, 0);
    for (auto& x : a) x = rnd(rng, 0, 1000);
    for (int i = 0; i < N; ++i) pre[i + 1] = pre[i] + a[i];
    std::vector<ll> x(N), y(N);
    for (auto& v : x) v = rnd(rng, -1000000, 1000000);
    for (auto& v : y) v = rnd(rng, -1000000, 1000000);
    std::sort(x.begin(), x.end()), std::sort(y.begin(), y.end());
    ll chk = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int r : monotone_minima(N, N, [&](int i, int j) { return (x[i] - y[j]) * (x[i] - y[j]); })) chk += r;
    const ll t1 = elapsed_ms(t0);
    auto c = [&](int i, int j) { return (pre[j] - pre[i]) * (pre[j] - pre[i]); };
    t0 = std::chrono::steady_clock::now();
    chk += partition_monge(N2, 10, c)[10];
    const ll t2 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    chk += monge_1d_dp(N, [&](int i, int j) { return c(i, j) + 1000000; })[N];
    const ll t3 = elapsed_ms(t0);
    std::cout << "[dp_monge] benchmark: monotone_minima(" << N << " × " << N << ") " << t1 << " ms / partition_monge(N = " << N2 << ", K = 10) " << t2
              << " ms / monge_1d_dp(N = " << N << ") " << t3 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace monge_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace monge_demo {

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
#define SHOW(expr, comment) monge_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    monge_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("Monge・決定単調性の DP（dp_monge） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  monotone_minima(H, W, f) / partition_monge(N, K, c) / monge_1d_dp(N, c) / alien_dp(K, lo, hi, solve) / is_monge(H, W, f)\n";

    RUN(const std::vector<long long> a = {3, 1, 4, 1, 5, 9, 2, 6});
    RUN(std::vector<long long> pre(9, 0));
    RUN(for (int i = 0; i < 8; ++i) pre[i + 1] = pre[i] + a[i]);
    RUN(auto c = [&](int i, int j) { return (pre[j] - pre[i]) * (pre[j] - pre[i]); });  // 区間の和の 2 乗（Monge）

    api("【1】 partition_monge(N, K, c)", "vector<T>（長さ K + 1）", "[0, N) を空でない k 個の区間に分けて c(i, j) の和を最小に（k = 0〜K）",
        "O(K N log N)");
    SHOW(partition_monge(8, 4, c), "k = 0 は作れない（大きい値）。k = 1 は全体 31^2");

    api("【2】 monge_1d_dp(N, c)", "vector<T>（長さ N + 1）", "dp[0] = 0、dp[j] = min_{i<j} dp[i] + c(i, j)。区間の数は自由", "O(N log N)");
    SHOW(monge_1d_dp(8, [&](int i, int j) { return c(i, j) + 20; }), "1 区間ごとに 20 かかるときの最小");

    api("【3】 monotone_minima(H, W, f)", "vector<int>", "各行で f(i, j) を最小にする列（Monge なら O((H + W) log H)）", "O((H + W) log H)");
    RUN(const std::vector<long long> xs = {0, 3, 5, 10}, ys = {1, 2, 6, 9});
    SHOW(monotone_minima(4, 4, [&](int i, int j) { return (xs[i] - ys[j]) * (xs[i] - ys[j]); }), "各 xs[i] に一番近い ys の番号");

    api("【4】 alien_dp(K, lo, hi, solve)", "long long", "solve(λ) = {min (コスト + λ × 個数), 最小の個数}。ちょうど K 個のときの最小コスト",
        "O(log(hi − lo) × solve)");
    code("// 隣り合わない要素をちょうど K 個選んで和を最大に（コスト = −和 として最小化）");
    RUN(auto solve = [&](long long lam) { std::pair<long long, long long> no{0, 0}, yes{(long long)4e18, 0}; for (long long x : a) { std::pair<long long, long long> ny{no.first - x + lam, no.second + 1}; no = std::min(no, yes); yes = ny; } return std::min(no, yes); });
    SHOW(-alien_dp(2, -100, 100, solve), "K = 2：9 + 6 など");
    SHOW(-alien_dp(3, -100, 100, solve), "K = 3：4 + 9 + 6 など");
    note("λ の範囲は広めに（迷ったら ±(コストの絶対値の合計 + 1)）。狭いと答えがずれる");

    api("【5】 is_monge(H, W, f)", "bool", "Monge か確かめる（小さい入力で使う）", "O(HW)");
    SHOW(is_monge(9, 9, [&](int i, int j) { return i < j ? c(i, j) : 1000000000LL * (i - j + 1); }), "i < j だけ使うときは、それ以外を大きくして確かめる");
}

}  // namespace monge_demo

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
    monge_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    monge_test::benchmark();    // 速度計測
    monge_demo::demo();         // 使い方デモ
    return 0;
}
