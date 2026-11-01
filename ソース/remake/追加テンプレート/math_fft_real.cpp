// =============================================================================
//  実数の FFT（double の畳み込み：確率の分布の和・誤差を許す整数の畳み込み）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  c[k] = Σ_{i + j = k} a[i] b[j] を複素数の FFT で O(N log N)。mod の畳み込み（NTT）ではなく実数のまま計算するので、
//  確率・期待値の分布を足し合わせる（さいころの和など）ときに使う。整数でも結果が 10^13 くらいまでなら四捨五入で正しい
//    例: 「N 個のさいころの和の確率分布」、「2 つの集合から 1 つずつ選んだ和ごとの組の数（値が小さいとき）」
//
//  機能一覧                                                                                  計算量
//    convolve_real(a, b)                 vector<double> の畳み込み（長さ |a| + |b| − 1。どちらか空なら空）         O(N log N)
//    convolve_round(a, b)                vector<long long> の畳み込みを double で計算して四捨五入
//                                        （長さ 10^5 で結果の各値が 3×10^13 くらいまで正しいのを確認。超えるなら math_v2/07_convolution.cpp）
//    power_real(a, k)                    a を k 回畳み込んだもの（分布を k 回足す。max_len を超えた分は捨てる）        O(N log N log k)
//    fft_inplace(v, inverse)             複素数の FFT そのもの（v の長さは 2 の冪。inverse なら長さで割る）
//
//  注意: double の誤差は 最大の値 × 長さ × 10^−15 くらい。確率なら気にしなくてよい
//  速さの目安: 長さ 2^19 どうしの convolve_real 約 0.1 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    fft_test::random_test()  … O(N^2) の愚直な畳み込みとの突き合わせ
//    fft_test::benchmark()    … 速度計測
//    fft_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
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
namespace fft_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace fft_internal

inline void fft_inplace(std::vector<std::complex<double>>& v, bool inverse) {
    const int n = (int)v.size();
    if (n == 0) return;
    if (n & (n - 1)) fft_internal::fail("fft_inplace", "長さは 2 の冪");
    for (int i = 1, j = 0; i < n; ++i) {  // ビット反転の並べ替え
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(v[i], v[j]);
    }
    // 回転因子は 1 回だけ作る（長さ n/2 の表。誤差を小さくするため毎回 cos / sin で）
    static std::vector<std::complex<double>> rt;
    if ((int)rt.size() < n) {
        rt.assign(std::max(2, n), {1, 0});
        for (int k = 2; k < n; k <<= 1) {
            const double ang = std::acos(-1.0) / k;
            for (int i = k; i < 2 * k; ++i) rt[i] = (i & 1) ? std::polar(1.0, ang * (i - k)) : rt[i / 2];
        }
    }
    for (int len = 1; len < n; len <<= 1)
        for (int i = 0; i < n; i += 2 * len)
            for (int j = 0; j < len; ++j) {
                std::complex<double> w = rt[j + len];
                if (inverse) w = std::conj(w);
                const std::complex<double> z = v[i + j + len] * w;
                v[i + j + len] = v[i + j] - z;
                v[i + j] += z;
            }
    if (inverse)
        for (auto& x : v) x /= n;
}
inline std::vector<double> convolve_real(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.empty() || b.empty()) return {};
    const int need = (int)(a.size() + b.size() - 1);
    if (std::min(a.size(), b.size()) <= 32) {  // 短いなら愚直の方が速くて正確
        std::vector<double> c(need, 0.0);
        for (size_t i = 0; i < a.size(); ++i)
            for (size_t j = 0; j < b.size(); ++j) c[i + j] += a[i] * b[j];
        return c;
    }
    int n = 1;
    while (n < need) n <<= 1;
    // 1 回の FFT で 2 本：in = a + i b、out = in^2 の虚部 / 2 = a * b
    std::vector<std::complex<double>> in(n);
    for (size_t i = 0; i < a.size(); ++i) in[i].real(a[i]);
    for (size_t i = 0; i < b.size(); ++i) in[i].imag(b[i]);
    fft_inplace(in, false);
    for (auto& x : in) x *= x;
    fft_inplace(in, true);
    std::vector<double> c(need);
    for (int i = 0; i < need; ++i) c[i] = in[i].imag() / 2;
    return c;
}
inline std::vector<long long> convolve_round(const std::vector<long long>& a, const std::vector<long long>& b) {
    std::vector<double> x(a.begin(), a.end()), y(b.begin(), b.end());
    const auto c = convolve_real(x, y);
    std::vector<long long> r(c.size());
    for (size_t i = 0; i < c.size(); ++i) r[i] = (long long)std::llround(c[i]);
    return r;
}
inline std::vector<double> power_real(const std::vector<double>& a, long long k, size_t max_len = (size_t)1 << 22) {
    if (k < 0) fft_internal::fail("power_real", "k は 0 以上");
    std::vector<double> res{1.0}, base = a;
    if (base.size() > max_len) base.resize(max_len);
    while (k > 0) {
        if (k & 1) {
            res = convolve_real(res, base);
            if (res.size() > max_len) res.resize(max_len);
        }
        k >>= 1;
        if (k) {
            base = convolve_real(base, base);
            if (base.size() > max_len) base.resize(max_len);
        }
    }
    return res;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace fft_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define FT_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[math_fft_real] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    std::uniform_real_distribution<double> ud(-1.0, 1.0);
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 300), m = (int)rnd(rng, 0, 300);
        std::vector<double> a(n), b(m);
        for (auto& x : a) x = ud(rng);
        for (auto& x : b) x = ud(rng);
        const auto c = convolve_real(a, b);
        FT_CHECK((int)c.size() == (n && m ? n + m - 1 : 0));
        for (int k = 0; k < (int)c.size(); ++k) {
            double s = 0;
            for (int i = std::max(0, k - m + 1); i <= std::min(k, n - 1); ++i) s += a[i] * b[k - i];
            FT_CHECK(std::abs(c[k] - s) < 1e-9);
        }
        std::vector<ll> x(n), y(m);
        for (auto& v : x) v = rnd(rng, -1000000, 1000000);
        for (auto& v : y) v = rnd(rng, -1000000, 1000000);
        const auto z = convolve_round(x, y);
        for (int k = 0; k < (int)z.size(); ++k) {
            ll s = 0;
            for (int i = std::max(0, k - m + 1); i <= std::min(k, n - 1); ++i) s += x[i] * y[k - i];
            FT_CHECK(z[k] == s);
        }
        cases += 2;
    }
    // power_real：さいころ 3 個の和
    {
        std::vector<double> die(7, 1.0 / 6);
        die[0] = 0;
        const auto p = power_real(die, 3);
        FT_CHECK(p.size() == 19 && std::abs(p[10] - 27.0 / 216) < 1e-12 && std::abs(p[3] - 1.0 / 216) < 1e-12);
        double s = 0;
        for (double v : p) s += v;
        FT_CHECK(std::abs(s - 1) < 1e-12);
        cases += 2;
    }
    std::cout << "[math_fft_real] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[math_fft_real] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? (1 << 16) : (1 << 19);
    std::vector<double> a(N), b(N);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    for (auto& x : a) x = ud(rng);
    for (auto& x : b) x = ud(rng);
    const auto t0 = std::chrono::steady_clock::now();
    const auto c = convolve_real(a, b);
    std::cout << "[math_fft_real] benchmark: 長さ " << N << " どうしの convolve_real " << elapsed_ms(t0) << " ms  (check " << (long long)c[N / 2] % 1000 << ")\n";
}

}  // namespace fft_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace fft_demo {

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
#define SHOW(expr, comment) fft_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    fft_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("実数の FFT（math_fft_real） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  convolve_real(a, b) / convolve_round(a, b) / power_real(a, k, max_len) / fft_inplace(v, inverse)\n";

    api("【1】 convolve_real(a, b)", "vector<double>", "c[k] = Σ a[i] b[k − i]（実数のまま）", "O(N log N)");
    SHOW(convolve_real(std::vector<double>{0.5, 0.5}, std::vector<double>{0.5, 0.5}), "コインを 2 回投げた表の数の分布");

    api("【2】 power_real(a, k)", "vector<double>", "分布を k 回足した分布（さいころ k 個の和など）", "O(N log N log k)");
    RUN(std::vector<double> die = {0, 1.0 / 6, 1.0 / 6, 1.0 / 6, 1.0 / 6, 1.0 / 6, 1.0 / 6});
    RUN(const auto p = power_real(die, 2));
    SHOW(p[7], "2 個の和が 7 になる確率 = 6/36");
    SHOW(power_real(die, 100)[350], "100 個の和が 350 になる確率");

    api("【3】 convolve_round(a, b)", "vector<long long>", "整数の畳み込み（結果が 10^13 くらいまで）", "O(N log N)");
    SHOW(convolve_round(std::vector<long long>{1, 2, 3}, std::vector<long long>{4, 5}), "(1 + 2x + 3x²)(4 + 5x)");
}

}  // namespace fft_demo

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
    fft_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    fft_test::benchmark();    // 速度計測
    fft_demo::demo();         // 使い方デモ
    return 0;
}
