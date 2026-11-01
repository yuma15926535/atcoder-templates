// =============================================================================
//  格子路の数え上げ（鏡像法：直線を越えない道・2 本の直線にはさまれる道・ballot 数・カタラン数）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  (0, 0) から (x, y) へ右（x + 1）か上（y + 1）に進む道で、途中ずっと y − x が [lo, hi] に入っているものの数を、
//  「境界に触れた道を鏡に映すと、映した点から出る道と 1 対 1 になる」ことを使って二項係数の足し引きで求める。
//  2 本のときは鏡に映すのを交互にくり返す（O((x + y) / (hi − lo)) 項）
//    例: 括弧列（カタラン数）、「A がずっと B 以上の票」（ballot）、「在庫が 0 未満にも上限超えにもならない並べ方」
//
//  機能一覧（答えは mod（素数。既定 998244353）。x, y, n ≤ max_n）                              計算量
//    LatticePath lp(max_n, mod)            階乗の表を作る                                              O(max_n)
//    lp.C(n, k)                            二項係数（範囲外は 0）                                       O(1)
//    lp.paths(x, y)                        (0, 0) → (x, y) の道の数 C(x + y, x)                         O(1)
//    lp.paths_below(x, y, k)               ずっと y − x ≤ k（k ≥ 0）の道の数（直線 y = x + k を越えない）  O(1)
//    lp.paths_between(x, y, lo, hi)        ずっと lo ≤ y − x ≤ hi（lo ≤ 0 ≤ hi）の道の数                O((x + y) / (hi − lo + 2) + 1)
//    lp.walks_in_range(n, D, lo, hi)       ±1 を n 回足して 0 から D へ、途中ずっと [lo, hi] にいる並びの数（上と同じもの）
//    lp.ballot(p, q)                       A が p 票・B が q 票（p > q）で、開票の途中ずっと A が B より多い順番の数（p ≤ q は 0）  O(1)
//    lp.catalan(n)                         カタラン数 C(2n, n) / (n + 1)                               O(1)
//
//  別の点から出る道は平行移動する：(a, b) → (c, d) で lo ≤ (y − x) ≤ hi なら、paths_between(c − a, d − b, lo − (b − a), hi − (b − a))
//  速さの目安: max_n = 10^6 の表 約 0.04 秒。paths_between は 1 回 O(n / 幅)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    lattice_test::random_test()  … 1 歩ずつ数える DP の愚直解との突き合わせ
//    lattice_test::benchmark()    … 速度計測
//    lattice_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
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
namespace lattice_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace lattice_internal

struct LatticePath {
    explicit LatticePath(int max_n = 0, long long mod = 998244353) : mod_(mod) {
        if (max_n < 0 || mod < 2) lattice_internal::fail("LatticePath", "max_n ≥ 0、mod は素数");
        if (max_n >= mod) lattice_internal::fail("LatticePath", "max_n は mod より小さくする（階乗が 0 になる）");
        const int n = 2 * max_n + 2;
        f_.assign(n + 1, 1), inv_f_.assign(n + 1, 1);
        for (int i = 1; i <= n; ++i) f_[i] = f_[i - 1] * (i % mod_) % mod_;
        inv_f_[n] = pw(f_[n], mod_ - 2);
        for (int i = n; i > 0; --i) inv_f_[i - 1] = inv_f_[i] * (i % mod_) % mod_;
    }
    long long C(long long n, long long k) const {
        if (k < 0 || n < 0 || k > n) return 0;
        if (n >= (long long)f_.size()) lattice_internal::fail("LatticePath::C", "n が max_n の 2 倍を超えた（表を大きく作る）");
        return f_[n] * inv_f_[k] % mod_ * inv_f_[n - k] % mod_;
    }
    long long paths(long long x, long long y) const { return (x < 0 || y < 0) ? 0 : C(x + y, x); }
    long long paths_below(long long x, long long y, long long k) const {  // ずっと y − x ≤ k
        if (k < 0) lattice_internal::fail("LatticePath::paths_below", "k ≥ 0（出発点が条件を満たす）");
        if (x < 0 || y < 0 || y - x > k) return 0;
        return (C(x + y, x) - C(x + y, x + k + 1) + mod_) % mod_;  // 直線 y = x + k + 1 で (0, 0) を映すと (−k−1, k+1)
    }
    long long walks_in_range(long long n, long long D, long long lo, long long hi) const {
        if (!(lo <= 0 && 0 <= hi)) lattice_internal::fail("LatticePath::walks_in_range", "lo ≤ 0 ≤ hi にする");
        if (n < 0 || D < lo || D > hi) return 0;
        const long long w = hi - lo + 2;  // 2 つの鏡（hi + 1 と lo − 1）の間隔 × 1
        auto N = [&](long long t) -> long long {  // 制約なしで 0 → t（n 歩）
            if (t < -n || t > n || ((n + t) & 1)) return 0;
            return C(n, (n + t) / 2);
        };
        long long res = 0;
        // Σ_k [N(D + 2kw) − N(2(hi + 1) − D + 2kw)]
        const long long kmax = n / (2 * w) + 2;
        for (long long k = -kmax; k <= kmax; ++k) {
            res += N(D + 2 * k * w);
            res -= N(2 * (hi + 1) - D + 2 * k * w);
            res %= mod_;
        }
        return (res % mod_ + mod_) % mod_;
    }
    long long paths_between(long long x, long long y, long long lo, long long hi) const {  // ずっと lo ≤ y − x ≤ hi
        if (x < 0 || y < 0) return 0;
        return walks_in_range(x + y, y - x, lo, hi);  // 上 = +1、右 = −1
    }
    long long ballot(long long p, long long q) const {  // (p − q) / (p + q) × C(p + q, p)
        if (p <= q || q < 0) return 0;
        return (C(p + q - 1, p - 1) - C(p + q - 1, p) + mod_) % mod_;  // 最初は A。そのあと「A − B ≥ 1 を保つ」道の数
    }
    long long catalan(long long n) const { return n < 0 ? 0 : (C(2 * n, n) - C(2 * n, n + 1) + mod_) % mod_; }
    long long mod() const { return mod_; }

   private:
    long long mod_;
    std::vector<long long> f_, inv_f_;
    long long pw(long long a, long long e) const {
        long long r = 1;
        a %= mod_;
        for (; e; e >>= 1, a = a * a % mod_)
            if (e & 1) r = r * a % mod_;
        return r;
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace lattice_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define LP_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[math_lattice_path] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    const ll MOD = 998244353;
    const LatticePath lp(200, MOD);
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const ll x = rnd(rng, 0, 25), y = rnd(rng, 0, 25);
        const ll lo = -rnd(rng, 0, 12), hi = rnd(rng, 0, 12);
        // DP：dp[i][j] = (0,0) から (i, j) まで、条件を満たす道の数
        std::vector<std::vector<ll>> dp(x + 1, std::vector<ll>(y + 1, 0)), all(x + 1, std::vector<ll>(y + 1, 0)), below(x + 1, std::vector<ll>(y + 1, 0));
        for (ll i = 0; i <= x; ++i)
            for (ll j = 0; j <= y; ++j) {
                const bool ok = lo <= j - i && j - i <= hi, okb = j - i <= hi;
                if (i == 0 && j == 0) {
                    dp[0][0] = ok, all[0][0] = 1, below[0][0] = okb;
                    continue;
                }
                ll a = 0, b = 0, c = 0;
                if (i > 0) a += dp[i - 1][j], b += all[i - 1][j], c += below[i - 1][j];
                if (j > 0) a += dp[i][j - 1], b += all[i][j - 1], c += below[i][j - 1];
                dp[i][j] = ok ? a % MOD : 0, all[i][j] = b % MOD, below[i][j] = okb ? c % MOD : 0;
            }
        LP_CHECK(lp.paths_between(x, y, lo, hi) == dp[x][y]);
        LP_CHECK(lp.paths(x, y) == all[x][y]);
        LP_CHECK(lp.paths_below(x, y, hi) == below[x][y]);
        cases += 3;
    }
    // ballot と カタラン数：小さい n で全部の並びを数える
    for (int p = 0; p <= 9; ++p)
        for (int q = 0; q <= 9; ++q) {
            ll cnt = 0;
            for (int m = 0; m < (1 << (p + q)); ++m) {
                if (__builtin_popcount(m) != p) continue;
                int d = 0;
                bool ok = true;
                for (int i = 0; i < p + q && ok; ++i) d += (m >> i & 1) ? 1 : -1, ok = d > 0;
                cnt += ok;
            }
            LP_CHECK(lp.ballot(p, q) == (p > q ? cnt : 0));  // p ≤ q は 0（p = q = 0 も 0 とする）
            ++cases;
        }
    const ll cat[] = {1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796};
    for (int n = 0; n <= 10; ++n) LP_CHECK(lp.catalan(n) == cat[n]);
    std::cout << "[math_lattice_path] random_test: OK (" << cases + 11 << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[math_lattice_path] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 100000 : 1000000;
    auto t0 = std::chrono::steady_clock::now();
    const LatticePath lp(N);
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    ll chk = 0;
    for (int i = 0; i < 1000; ++i) chk += lp.paths_between(N / 2, N / 2, -1000 - i, 1000);
    const ll t2 = elapsed_ms(t0);
    std::cout << "[math_lattice_path] benchmark: 表（max_n = " << N << "） " << t1 << " ms / paths_between（x = y = " << N / 2 << "、幅 2000）を 1000 回 "
              << t2 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace lattice_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace lattice_demo {

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
#define SHOW(expr, comment) lattice_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    lattice_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("格子路の数え上げ（math_lattice_path） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  LatticePath lp(max_n, mod): C / paths / paths_below / paths_between / walks_in_range / ballot / catalan\n";

    RUN(const LatticePath lp(100));

    api("【1】 paths(x, y) / paths_below(x, y, k)", "long long", "(0, 0) → (x, y) の道の数 / ずっと y − x ≤ k の道の数（直線を越えない）", "O(1)");
    SHOW(lp.paths(3, 3), "C(6, 3)");
    SHOW(lp.paths_below(3, 3, 0), "対角線より上に出ない = カタラン数 5");
    SHOW(lp.paths_below(4, 2, 1), "");

    api("【2】 paths_between(x, y, lo, hi)", "long long", "ずっと lo ≤ y − x ≤ hi の道の数（2 本の直線にはさまれる）", "O((x + y) / (hi − lo + 2) + 1)");
    SHOW(lp.paths_between(3, 3, -1, 1), "対角線から 1 より離れない");
    SHOW(lp.paths_between(5, 5, 0, 1), "y − x が 0 か 1 だけ：上・右を交互 → 1 通り");
    code("// 例: 在庫 0 から始め、+1（入荷）を a 回・−1（出荷）を b 回。在庫がずっと 0〜C のままの順番の数");
    SHOW(lp.walks_in_range(6, 0, 0, 2), "+1 と −1 を 3 回ずつ、在庫 0〜2");

    api("【3】 ballot(p, q) / catalan(n)", "long long", "A が p 票・B が q 票で、途中ずっと A が多い開票順 / カタラン数", "O(1)");
    SHOW(lp.ballot(3, 1), "(3 − 1) / (3 + 1) × C(4, 3) = 2");
    SHOW(lp.ballot(5, 3), "");
    SHOW(lp.catalan(10), "");
}

}  // namespace lattice_demo

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
    lattice_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    lattice_test::benchmark();    // 速度計測
    lattice_demo::demo();         // 使い方デモ
    return 0;
}
