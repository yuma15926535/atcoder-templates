// =============================================================================
//  CDQ 分割統治（列を半分に分け、「左から右への寄与」だけを掃いて求める。3 次元の支配の数え上げ・DP）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  1 つ目の軸で並べた列を半分に分けると、「左の要素 → 右の要素」の組は 1 つ目の軸の条件をもう満たしているので、
//  残りの 2 つの軸だけを掃く（2 つ目で並べ、3 つ目は BIT）。全体 O(N log^2 N)。
//  左を先に全部決めてから右へ寄与させる順番（左 → 横断 → 右）にすると、DP にも使える
//    例: 「a_j ≤ a_i かつ b_j ≤ b_i かつ c_j ≤ c_i となる j の数」（3 次元の偏序）、
//        dp[i] = 1 + max{dp[j] : j < i, a_j ≤ a_i, b_j ≤ b_i}（2 次元の LIS）、時刻つきの点の追加と長方形の質問
//
//  機能一覧                                                                                  計算量
//    dominance_sum_3d(points, queries)     points = {a, b, c, w}、queries = {a, b, c}。質問ごとに、a ≤ qa・b ≤ qb・c ≤ qc の
//                                          点の w の和                                                         O((N + Q) log^2 (N + Q))
//    count_dominated_3d(points)            points = {a, b, c}。点ごとに、自分以外で 3 つとも ≤ の点の数
//    count_dominated_3d_strict(points)     3 つとも < の点の数（座標は整数）
//    cdq_divide(n, cross)                  [0, n) を分割統治し、「左 → 横断 → 右」の順に cross(l, m, r) を呼ぶ
//                                          （[l, m) が全部決まってから [m, r) へ寄与させる。cross の中身は自分で書く）   O(n log n) 回の分割
//
//  注意: 座標は long long。重みの和は long long。cross の中で並べ替えるときは自分のコピーで（元の順番を崩さない）
//  速さの目安: count_dominated_3d（N = 2×10^5）約 0.4 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    cdq_test::random_test()  … 全部の組を調べる O(N^2) の愚直解との突き合わせ
//    cdq_test::benchmark()    … 速度計測
//    cdq_demo::demo()         … 全機能の使い方デモ（日本語の説明付き。2 次元の LIS を cdq_divide で）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace cdq_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace cdq_internal

template <class Cross>
void cdq_divide(int n, Cross cross) {
    if (n < 0) cdq_internal::fail("cdq_divide", "n は 0 以上");
    auto rec = [&](auto&& self, int l, int r) -> void {
        if (r - l <= 1) return;
        const int m = (l + r) / 2;
        self(self, l, m);
        cross(l, m, r);
        self(self, m, r);
    };
    rec(rec, 0, n);
}

inline std::vector<long long> dominance_sum_3d(const std::vector<std::array<long long, 4>>& points, const std::vector<std::array<long long, 3>>& queries) {
    struct Ev {
        long long a, b, c, w;
        int q;  // −1 = 点、それ以外 = 質問の番号
    };
    std::vector<Ev> E;
    E.reserve(points.size() + queries.size());
    for (const auto& p : points) E.push_back({p[0], p[1], p[2], p[3], -1});
    for (size_t i = 0; i < queries.size(); ++i) E.push_back({queries[i][0], queries[i][1], queries[i][2], 0, (int)i});
    // a の順（同じ a なら点が先：a_p ≤ a_q を数えるため）
    std::sort(E.begin(), E.end(), [](const Ev& x, const Ev& y) {
        if (x.a != y.a) return x.a < y.a;
        return (x.q < 0) > (y.q < 0);
    });
    std::vector<long long> cs;
    for (const auto& p : points) cs.push_back(p[2]);
    std::sort(cs.begin(), cs.end());
    cs.erase(std::unique(cs.begin(), cs.end()), cs.end());
    const int m = (int)cs.size();
    std::vector<int> cid(E.size());  // 点：c の番号、質問：c 以下の番号の数
    for (size_t i = 0; i < E.size(); ++i)
        cid[i] = E[i].q < 0 ? (int)(std::lower_bound(cs.begin(), cs.end(), E[i].c) - cs.begin()) : (int)(std::upper_bound(cs.begin(), cs.end(), E[i].c) - cs.begin());
    std::vector<long long> ans(queries.size(), 0), bit(m + 1, 0);
    std::vector<int> idx(E.size()), buf(E.size());
    for (size_t i = 0; i < E.size(); ++i) idx[i] = (int)i;
    auto bless = [&](int x, int y) {  // b の順（同じなら点が先）
        if (E[x].b != E[y].b) return E[x].b < E[y].b;
        return (E[x].q < 0) > (E[y].q < 0);
    };
    // マージソートの形：idx[l..r) を b の順に並べながら、左の点を右の質問に足す
    auto rec = [&](auto&& self, int l, int r) -> void {
        if (r - l <= 1) return;
        const int mid = (l + r) / 2;
        self(self, l, mid), self(self, mid, r);
        int i = l;
        for (int j = mid; j < r; ++j) {
            const int y = idx[j];
            if (E[y].q < 0) continue;
            while (i < mid && bless(idx[i], y)) {
                const int x = idx[i++];
                if (E[x].q < 0)
                    for (int k = cid[x] + 1; k <= m; k += k & -k) bit[k] += E[x].w;
            }
            long long s = 0;
            for (int k = cid[y]; k > 0; k -= k & -k) s += bit[k];
            ans[E[y].q] += s;
        }
        for (int t = l; t < i; ++t) {
            const int x = idx[t];
            if (E[x].q < 0)
                for (int k = cid[x] + 1; k <= m; k += k & -k) bit[k] -= E[x].w;
        }
        std::merge(idx.begin() + l, idx.begin() + mid, idx.begin() + mid, idx.begin() + r, buf.begin() + l, bless);
        std::copy(buf.begin() + l, buf.begin() + r, idx.begin() + l);
    };
    rec(rec, 0, (int)E.size());
    return ans;
}
inline std::vector<long long> count_dominated_3d(const std::vector<std::array<long long, 3>>& points) {
    std::vector<std::array<long long, 4>> p;
    for (const auto& x : points) p.push_back({x[0], x[1], x[2], 1});
    auto r = dominance_sum_3d(p, points);
    for (auto& v : r) --v;  // 自分を除く
    return r;
}
inline std::vector<long long> count_dominated_3d_strict(const std::vector<std::array<long long, 3>>& points) {
    std::vector<std::array<long long, 4>> p;
    std::vector<std::array<long long, 3>> q;
    for (const auto& x : points) p.push_back({x[0], x[1], x[2], 1}), q.push_back({x[0] - 1, x[1] - 1, x[2] - 1});
    return dominance_sum_3d(p, q);
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace cdq_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define CD_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_cdq] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 40), q = (int)rnd(rng, 0, 30);
        const ll V = rnd(rng, 1, 6);
        std::vector<std::array<ll, 4>> pts(n);
        std::vector<std::array<ll, 3>> qs(q), p3(n);
        for (int i = 0; i < n; ++i) pts[i] = {rnd(rng, 0, V), rnd(rng, 0, V), rnd(rng, 0, V), rnd(rng, -5, 5)}, p3[i] = {pts[i][0], pts[i][1], pts[i][2]};
        for (auto& x : qs) x = {rnd(rng, -1, V + 1), rnd(rng, -1, V + 1), rnd(rng, -1, V + 1)};
        const auto res = dominance_sum_3d(pts, qs);
        for (int i = 0; i < q; ++i) {
            ll s = 0;
            for (const auto& p : pts)
                if (p[0] <= qs[i][0] && p[1] <= qs[i][1] && p[2] <= qs[i][2]) s += p[3];
            CD_CHECK(res[i] == s);
            ++cases;
        }
        const auto c1 = count_dominated_3d(p3), c2 = count_dominated_3d_strict(p3);
        for (int i = 0; i < n; ++i) {
            ll a = 0, b = 0;
            for (int j = 0; j < n; ++j) {
                if (j != i && p3[j][0] <= p3[i][0] && p3[j][1] <= p3[i][1] && p3[j][2] <= p3[i][2]) ++a;
                if (p3[j][0] < p3[i][0] && p3[j][1] < p3[i][1] && p3[j][2] < p3[i][2]) ++b;
            }
            CD_CHECK(c1[i] == a && c2[i] == b);
            ++cases;
        }
        // cdq_divide で 2 次元の LIS：dp[i] = 1 + max{dp[j] : j < i, a_j ≤ a_i, b_j ≤ b_i}
        {
            std::vector<ll> A(n), B(n);
            for (int i = 0; i < n; ++i) A[i] = rnd(rng, 0, V), B[i] = rnd(rng, 0, V);
            std::vector<int> dp(n, 1);
            std::vector<int> bit(V + 2, 0);
            cdq_divide(n, [&](int l, int m, int r) {
                std::vector<int> L, R;
                for (int i = l; i < m; ++i) L.push_back(i);
                for (int i = m; i < r; ++i) R.push_back(i);
                std::sort(L.begin(), L.end(), [&](int x, int y) { return A[x] < A[y]; });
                std::sort(R.begin(), R.end(), [&](int x, int y) { return A[x] < A[y]; });
                size_t p = 0;
                for (int y : R) {
                    while (p < L.size() && A[L[p]] <= A[y]) {
                        for (ll k = B[L[p]] + 1; k <= V + 1; k += k & -k) bit[k] = std::max(bit[k], dp[L[p]]);
                        ++p;
                    }
                    int best = 0;
                    for (ll k = B[y] + 1; k > 0; k -= k & -k) best = std::max(best, bit[k]);
                    dp[y] = std::max(dp[y], best + 1);
                }
                for (size_t t = 0; t < p; ++t)
                    for (ll k = B[L[t]] + 1; k <= V + 1; k += k & -k) bit[k] = 0;
            });
            std::vector<int> exp(n, 1);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < i; ++j)
                    if (A[j] <= A[i] && B[j] <= B[i]) exp[i] = std::max(exp[i], exp[j] + 1);
            CD_CHECK(dp == exp);
            ++cases;
        }
    }
    std::cout << "[ds_cdq] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_cdq] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    std::vector<std::array<ll, 3>> p(N);
    for (auto& x : p) x = {rnd(rng, 0, 1000000000), rnd(rng, 0, 1000000000), rnd(rng, 0, 1000000000)};
    const auto t0 = std::chrono::steady_clock::now();
    ll chk = 0;
    for (ll v : count_dominated_3d(p)) chk += v;
    std::cout << "[ds_cdq] benchmark: count_dominated_3d（N = " << N << "） " << elapsed_ms(t0) << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace cdq_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace cdq_demo {

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
#define SHOW(expr, comment) cdq_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    cdq_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("CDQ 分割統治（ds_cdq） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  dominance_sum_3d(points{a, b, c, w}, queries{a, b, c}) / count_dominated_3d(points) / count_dominated_3d_strict(points)\n"
              << "  cdq_divide(n, cross)\n";

    api("【1】 count_dominated_3d(points) / count_dominated_3d_strict(points)", "vector<long long>",
        "点ごとに、自分以外で a・b・c が全部 ≤ の点の数 / 全部 < の点の数", "O(N log^2 N)");
    RUN(const std::vector<std::array<long long, 3>> p = {{1, 1, 1}, {2, 2, 2}, {2, 2, 2}, {3, 1, 5}, {0, 4, 0}});
    SHOW(count_dominated_3d(p), "");
    SHOW(count_dominated_3d_strict(p), "");

    api("【2】 dominance_sum_3d(points, queries)", "vector<long long>", "質問 (qa, qb, qc) ごとに、a ≤ qa・b ≤ qb・c ≤ qc の点の重みの和", "O((N + Q) log^2)");
    SHOW(dominance_sum_3d({{1, 1, 1, 10}, {2, 2, 2, 5}, {3, 1, 5, 1}}, {{2, 2, 2}, {3, 3, 9}, {0, 0, 0}}), "");

    api("【3】 cdq_divide(n, cross)", "void", "「左 → 横断 → 右」の順に cross(l, m, r)。左の答えが全部決まってから右に寄与させる（DP に使える）",
        "O(n log n) 回 × cross");
    {
        code("// 2 次元の LIS：dp[i] = 1 + max{dp[j] : j < i, a_j ≤ a_i, b_j ≤ b_i}。cross で左を a の順に BIT（b の max）へ入れる");
        const std::vector<int> A = {3, 1, 2, 2, 4, 1}, B = {1, 2, 2, 3, 3, 0};
        const int n = 6;
        std::vector<int> dp(n, 1), bit(6, 0);
        cdq_divide(n, [&](int l, int m, int r) {
            std::vector<int> L, R;
            for (int i = l; i < m; ++i) L.push_back(i);
            for (int i = m; i < r; ++i) R.push_back(i);
            std::sort(L.begin(), L.end(), [&](int x, int y) { return A[x] < A[y]; });
            std::sort(R.begin(), R.end(), [&](int x, int y) { return A[x] < A[y]; });
            size_t pi = 0;
            for (int y : R) {
                while (pi < L.size() && A[L[pi]] <= A[y]) {
                    for (int k = B[L[pi]] + 1; k < 6; k += k & -k) bit[k] = std::max(bit[k], dp[L[pi]]);
                    ++pi;
                }
                int best = 0;
                for (int k = B[y] + 1; k > 0; k -= k & -k) best = std::max(best, bit[k]);
                dp[y] = std::max(dp[y], best + 1);
            }
            for (size_t t = 0; t < pi; ++t)
                for (int k = B[L[t]] + 1; k < 6; k += k & -k) bit[k] = 0;
        });
        note("a = " + to_s(A) + "、b = " + to_s(B) + " → dp = " + to_s(dp));
    }
}

}  // namespace cdq_demo

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
    cdq_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    cdq_test::benchmark();    // 速度計測
    cdq_demo::demo();         // 使い方デモ
    return 0;
}
