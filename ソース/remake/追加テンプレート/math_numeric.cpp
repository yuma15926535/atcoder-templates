// =============================================================================
//  数値計算（積分：適応シンプソン・ガウス求積・2 重積分 / 方程式：ニュートン法・二分法・解を全部）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  式で書けない面積・体積を ∫ f(x) dx で、解の公式が無い方程式 f(x) = 0 を反復で、double のまま求める。
//  積分は「区間を 2 つに分けても値がほとんど変わらなければ終わり」（適応シンプソン）、解は「接線で近づく」（ニュートン法）
//    例: 「曲線で囲まれた面積」「回転体の体積 ∫ π r(x)² dx」「円や放物線が重なった部分の面積」「x^x = 10 の解」
//
//  機能一覧（f は double を受け取って double を返す関数・ラムダ）                                   計算量（f を呼ぶ回数）
//    integrate(f, a, b, eps = 1e-10)          ∫_a^b f(x) dx（適応シンプソン。eps は誤差の目安：絶対誤差か相対誤差のどちらか） 滑らかなら数百回
//    integrate_n(f, a, b, n)                  区間を n 等分したシンプソン則（n は偶数に切り上げ。適応で止まらない形に）  n + 1 回
//    integrate_gauss(f, a, b, pieces = 1)     5 点のガウス求積を pieces 等分で（多項式なら 9 次まで正確）           5 × pieces 回
//    integrate_2d(f, xa, xb, ya, yb, eps)     ∫_xa^xb ∫_ya(x)^yb(x) f(x, y) dy dx（ya, yb は x の関数）             積分 × 積分
//    newton(f, df, x0)                        f(x) = 0 の解を x0 から（df は導関数）。止まらない・df = 0 なら nullopt       数回〜数十回
//    newton_numeric(f, x0)                    導関数を数値微分で（df を書かなくてよい）
//    bisect_root(f, lo, hi)                   f(lo) と f(hi) の符号が違うとき、間の解（100 回で幅 / 2^100）          100 回
//    find_roots(f, lo, hi, samples = 10000)   [lo, hi] の解を全部（samples 等分して符号が変わる所を二分法）       samples + 解の数 × 100 回
//    derivative(f, x)                         f'(x) を中心差分で（誤差 約 1e-10 × 大きさ）                         2 回
//
//  注意: 積分区間の中で f が飛ぶ（不連続・特異点）なら、その点で区間を分けて足す。周期的な関数で最初の 5 点が全部 0 だと
//        止まるので、最初に 2^3 等分してから始めている。find_roots は重解（x 軸に接するだけ）を見つけられない
//        最小値・最大値（凸関数）は binary_search.cpp の argmin_convex_real（黄金分割探索）を使う
//  速さの目安: integrate(sin(x) e^x, 0, 30) が f を 約 4000 回、球の体積の integrate_2d が f を 約 2.6 万回（どちらも 0.001 秒くらい）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    numeric_test::random_test()  … 原始関数が分かる多項式・解が分かる式との突き合わせ
//    numeric_test::benchmark()    … 速度計測
//    numeric_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
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
namespace numeric_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
template <class F>
double simpson_rec(F& f, double a, double b, double fa, double fm, double fb, double whole, double eps_abs, double eps_rel, int depth, int min_depth) {
    const double m = (a + b) / 2, lm = (a + m) / 2, rm = (m + b) / 2;
    const double flm = f(lm), frm = f(rm);
    const double left = (m - a) / 6 * (fa + 4 * flm + fm), right = (b - m) / 6 * (fm + 4 * frm + fb);
    const double delta = left + right - whole;
    // 絶対誤差 eps_abs（半分ずつに分ける）か、相対誤差 eps_rel のどちらかを満たせば終わり
    if (min_depth <= 0 && (depth <= 0 || std::abs(delta) <= 15 * std::max(eps_abs, eps_rel * std::abs(left + right)))) return left + right + delta / 15;
    return simpson_rec(f, a, m, fa, flm, fm, left, eps_abs / 2, eps_rel, depth - 1, min_depth - 1) +
           simpson_rec(f, m, b, fm, frm, fb, right, eps_abs / 2, eps_rel, depth - 1, min_depth - 1);
}
}  // namespace numeric_internal

template <class F>
double integrate(F f, double a, double b, double eps = 1e-10, int max_depth = 50) {
    if (a == b) return 0.0;
    const double fa = f(a), fb = f(b), fm = f((a + b) / 2);
    return numeric_internal::simpson_rec(f, a, b, fa, fm, fb, (b - a) / 6 * (fa + 4 * fm + fb), eps, eps, max_depth, 3);
}
template <class F>
double integrate_n(F f, double a, double b, int n) {
    if (n <= 0) numeric_internal::fail("integrate_n", "n は 1 以上");
    if (n & 1) ++n;
    const double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; ++i) s += f(a + h * i) * (i & 1 ? 4 : 2);
    return s * h / 3;
}
template <class F>
double integrate_gauss(F f, double a, double b, int pieces = 1) {
    if (pieces <= 0) numeric_internal::fail("integrate_gauss", "pieces は 1 以上");
    static const double X[5] = {0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640};
    static const double W[5] = {0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891};
    const double h = (b - a) / pieces;
    double s = 0;
    for (int p = 0; p < pieces; ++p) {
        const double l = a + h * p, c = l + h / 2;
        double t = 0;
        for (int k = 0; k < 5; ++k) t += W[k] * f(c + h / 2 * X[k]);
        s += t * h / 2;
    }
    return s;
}
template <class F, class YA, class YB>
double integrate_2d(F f, double xa, double xb, YA ya, YB yb, double eps = 1e-8) {
    auto inner = [&](double x) { return integrate([&](double y) { return f(x, y); }, ya(x), yb(x), eps, 40); };
    return integrate(inner, xa, xb, eps, 40);
}
template <class F, class DF>
std::optional<double> newton(F f, DF df, double x0, int max_iter = 100, double eps = 1e-12) {
    double x = x0;
    for (int it = 0; it < max_iter; ++it) {
        const double fx = f(x), d = df(x);
        if (fx == 0) return x;
        if (d == 0 || !std::isfinite(d) || !std::isfinite(fx)) return std::nullopt;
        const double nx = x - fx / d;
        if (!std::isfinite(nx)) return std::nullopt;
        if (std::abs(nx - x) <= eps * std::max(1.0, std::abs(x))) return nx;
        x = nx;
    }
    return std::nullopt;
}
template <class F>
double derivative(F f, double x) {
    const double h = 6e-6 * std::max(1.0, std::abs(x));  // 中心差分は h ≈ (double の誤差)^(1/3) がちょうどよい
    return (f(x + h) - f(x - h)) / (2 * h);
}
template <class F>
std::optional<double> newton_numeric(F f, double x0, int max_iter = 100, double eps = 1e-12) {
    return newton(f, [&](double x) { return derivative(f, x); }, x0, max_iter, eps);
}
template <class F>
double bisect_root(F f, double lo, double hi, int iters = 100) {
    double flo = f(lo), fhi = f(hi);
    if (flo == 0) return lo;
    if (fhi == 0) return hi;
    if ((flo > 0) == (fhi > 0)) numeric_internal::fail("bisect_root", "f(lo) と f(hi) の符号が同じ");
    for (int it = 0; it < iters; ++it) {
        const double mid = lo + (hi - lo) / 2;
        const double fm = f(mid);
        if (fm == 0) return mid;
        if ((fm > 0) == (flo > 0)) lo = mid, flo = fm;
        else hi = mid, fhi = fm;
    }
    return lo + (hi - lo) / 2;
}
template <class F>
std::vector<double> find_roots(F f, double lo, double hi, int samples = 10000) {
    if (samples <= 0) numeric_internal::fail("find_roots", "samples は 1 以上");
    std::vector<double> res;
    const double h = (hi - lo) / samples;
    double px = lo, pf = f(lo);
    if (pf == 0) res.push_back(lo);
    for (int i = 1; i <= samples; ++i) {
        const double x = i == samples ? hi : lo + h * i, fx = f(x);
        if (fx == 0) res.push_back(x);
        else if (pf != 0 && (pf > 0) != (fx > 0)) res.push_back(bisect_root(f, px, x));
        px = x, pf = fx;
    }
    return res;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace numeric_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define NU_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[math_numeric] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
using VD = std::vector<double>;

static double poly(const VD& c, double x) {
    double s = 0;
    for (int i = (int)c.size() - 1; i >= 0; --i) s = s * x + c[i];
    return s;
}
static double poly_int(const VD& c, double a, double b) {  // 原始関数で
    VD C(c.size() + 1, 0);
    for (size_t i = 0; i < c.size(); ++i) C[i + 1] = c[i] / (double)(i + 1);
    return poly(C, b) - poly(C, a);
}
static bool near(double x, double y, double tol) { return std::abs(x - y) <= tol * std::max(1.0, std::abs(y)); }

static void random_test() {
    std::mt19937_64 rng(20261031);
    std::uniform_real_distribution<double> ud(-3.0, 3.0);
    long long cases = 0;
    const double PI = std::acos(-1.0);
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int deg = (int)rnd(rng, 0, 9);
        VD c(deg + 1);
        for (auto& v : c) v = ud(rng);
        double a = ud(rng), b = ud(rng);
        if (it % 5 == 0) std::swap(a, b);  // a > b でもよい（符号が逆）
        auto f = [&](double x) { return poly(c, x); };
        const double ex = poly_int(c, a, b);
        NU_CHECK(near(integrate(f, a, b), ex, 1e-8));
        NU_CHECK(near(integrate_gauss(f, a, b), ex, 1e-9));
        NU_CHECK(near(integrate_n(f, a, b, 2000), ex, 1e-6));
        // 2 重積分：長方形の上の x^i y^j の和
        const int dx = (int)rnd(rng, 0, 3), dy = (int)rnd(rng, 0, 3);
        const double y0 = ud(rng), y1 = ud(rng), k = ud(rng);
        auto g = [&](double x, double y) { return k * std::pow(x, dx) * std::pow(y, dy); };
        const double ex2 = k * (std::pow(b, dx + 1) - std::pow(a, dx + 1)) / (dx + 1) * (std::pow(y1, dy + 1) - std::pow(y0, dy + 1)) / (dy + 1);
        if (it % 10 == 0) {
            NU_CHECK(near(integrate_2d(g, a, b, [&](double) { return y0; }, [&](double) { return y1; }), ex2, 1e-6));
            ++cases;
        }
        // 解：(x − r)(x² + p)（p > 0 で実数解は r だけ）
        const double r = ud(rng), p = std::abs(ud(rng)) + 0.1;
        auto h = [&](double x) { return (x - r) * (x * x + p); };
        auto dh = [&](double x) { return (x * x + p) + 2 * x * (x - r); };
        NU_CHECK(near(bisect_root(h, -10.0, 10.0), r, 1e-12));
        // 極値をまたぐと周期に入ることがある（3 次式のニュートン法）。解より外側（|x| が大きい側）から始めれば必ず近づく
        const double sg = r >= 0 ? 1.0 : -1.0;
        const auto nr = newton(h, dh, r + sg * std::abs(ud(rng)));
        NU_CHECK(nr && near(*nr, r, 1e-10));
        const auto nn = newton_numeric(h, r + sg * 0.5);
        NU_CHECK(nn && near(*nn, r, 1e-9));
        NU_CHECK(near(derivative(h, a), dh(a), 1e-8));
        // find_roots：0.2 以上離れた解を持つ多項式
        const int nroot = (int)rnd(rng, 0, 5);
        VD rs;
        while ((int)rs.size() < nroot) {
            const double t = ud(rng);
            bool ok = true;
            for (double u : rs) ok &= std::abs(u - t) > 0.2;
            if (ok) rs.push_back(t);
        }
        std::sort(rs.begin(), rs.end());
        auto q = [&](double x) {
            double s = 1;
            for (double u : rs) s *= x - u;
            return s * (x * x + 1);
        };
        const auto fr = find_roots(q, -3.5, 3.5, 1000);
        NU_CHECK(fr.size() == rs.size());
        for (size_t i = 0; i < fr.size() && i < rs.size(); ++i) NU_CHECK(std::abs(fr[i] - rs[i]) < 1e-9);
        cases += 9;
    }
    // 滑らかでない・特別な関数
    NU_CHECK(near(integrate([](double x) { return std::sin(x); }, 0.0, PI), 2.0, 1e-9));
    NU_CHECK(near(integrate([](double x) { return std::sqrt(std::max(0.0, 1 - x * x)); }, -1.0, 1.0), PI / 2, 1e-8));    // 端で傾きが無限
    NU_CHECK(near(integrate([](double x) { return std::sin(x) * std::sin(x); }, 0.0, 2 * PI * 8), 8 * PI, 1e-8));      // 周期的
    NU_CHECK(near(integrate([](double x) { return std::abs(x - 0.3); }, -1.0, 1.0), (1.3 * 1.3 + 0.7 * 0.7) / 2, 1e-9));  // 折れ曲がる
    NU_CHECK(near(integrate([](double x) { return 1 / (1 + x * x); }, 0.0, 1.0), PI / 4, 1e-10));
    NU_CHECK(near(integrate_gauss([](double x) { return std::exp(x); }, 0.0, 5.0, 100), std::exp(5.0) - 1, 1e-12));
    NU_CHECK(integrate([](double x) { return x; }, 2.0, 2.0) == 0.0);
    // 球の体積（上半分の高さ × 2 を円の上で）
    {
        const double R = 1.5;
        const double v = integrate_2d([&](double x, double y) { return 2 * std::sqrt(std::max(0.0, R * R - x * x - y * y)); }, -R, R,
                                      [&](double x) { return -std::sqrt(std::max(0.0, R * R - x * x)); }, [&](double x) { return std::sqrt(std::max(0.0, R * R - x * x)); });
        NU_CHECK(near(v, 4 * PI * R * R * R / 3, 1e-6));
    }
    // ニュートン法が止まらない・失敗する例
    NU_CHECK(!newton([](double x) { return x * x + 1; }, [](double x) { return 2 * x; }, 0.0).has_value());  // df = 0
    NU_CHECK(!newton([](double x) { return std::cbrt(x); }, [](double x) { return 1 / (3 * std::cbrt(x * x)); }, 1.0).has_value());  // 発散
    {
        const auto x = newton_numeric([](double t) { return std::pow(t, t) - 10; }, 2.0);
        NU_CHECK(x && std::abs(std::pow(*x, *x) - 10) < 1e-10);
    }
    cases += 13;
    std::cout << "[math_numeric] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[math_numeric] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    {
        long long calls = 0;
        const auto t0 = std::chrono::steady_clock::now();
        const double v = integrate(
            [&](double x) {
                ++calls;
                return std::sin(x) * std::exp(x);
            },
            0.0, 30.0);
        const double ex = std::exp(30.0) * (std::sin(30.0) - std::cos(30.0)) / 2 + 0.5;
        std::cout << "[math_numeric] benchmark: integrate(sin(x) e^x, 0, 30) " << elapsed_ms(t0) << " ms, f を " << calls << " 回  (相対誤差 " << std::abs(v - ex) / std::abs(ex)
                  << ")\n";
    }
    {
        const double R = 1.0;
        long long calls = 0;
        const auto t0 = std::chrono::steady_clock::now();
        const double v = integrate_2d(
            [&](double x, double y) {
                ++calls;
                return 2 * std::sqrt(std::max(0.0, R * R - x * x - y * y));
            },
            -R, R, [&](double x) { return -std::sqrt(std::max(0.0, R * R - x * x)); }, [&](double x) { return std::sqrt(std::max(0.0, R * R - x * x)); });
        std::cout << "[math_numeric] benchmark: 球の体積 integrate_2d " << elapsed_ms(t0) << " ms, f を " << calls << " 回  (誤差 " << std::abs(v - 4 * std::acos(-1.0) / 3)
                  << ")\n";
    }
}

}  // namespace numeric_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace numeric_demo {

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
#define SHOW(expr, comment) numeric_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    numeric_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("数値計算（math_numeric） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  integrate(f, a, b, eps) / integrate_n(f, a, b, n) / integrate_gauss(f, a, b, pieces) / integrate_2d(f, xa, xb, ya, yb, eps)\n"
              << "  newton(f, df, x0) / newton_numeric(f, x0) / bisect_root(f, lo, hi) / find_roots(f, lo, hi, samples) / derivative(f, x)\n";
    api("【1】 integrate(f, a, b, eps = 1e-10)", "double", "∫_a^b f(x) dx（適応シンプソン）", "滑らかなら f を数百回");
    SHOW(integrate([](double x) { return std::sqrt(1 - x * x); }, -1.0, 1.0), "半径 1 の半円の面積 = π / 2");
    SHOW(integrate([](double x) { return std::acos(-1.0) * x * x; }, 0.0, 3.0), "y = x を x 軸で回した円すい（底の半径 3、高さ 3）= 9π");
    SHOW(integrate([](double x) { return std::abs(std::sin(x)); }, 0.0, 10.0), "折れ曲がる関数も（区間を分けると速い）");

    api("【2】 integrate_n(f, a, b, n) / integrate_gauss(f, a, b, pieces)", "double", "n 等分のシンプソン則 / 5 点ガウス求積（多項式は 9 次まで正確）", "n + 1 回 / 5 × pieces 回");
    SHOW(integrate_n([](double x) { return std::exp(x); }, 0.0, 1.0, 100), "e − 1");
    SHOW(integrate_gauss([](double x) { return x * x * x * x; }, 0.0, 2.0), "32 / 5（1 区間で正確）");

    api("【3】 integrate_2d(f, xa, xb, ya(x), yb(x), eps)", "double", "y の範囲が x で変わる 2 重積分（体積など）", "積分 × 積分");
    code("// 半径 1 の球：高さ 2√(1 − x² − y²) を、円 x² + y² ≤ 1 の上で積分");
    RUN(auto top = [](double x, double y) { return 2 * std::sqrt(std::max(0.0, 1 - x * x - y * y)); });
    RUN(auto lo = [](double x) { return -std::sqrt(std::max(0.0, 1 - x * x)); });
    RUN(auto hi = [](double x) { return std::sqrt(std::max(0.0, 1 - x * x)); });
    SHOW(integrate_2d(top, -1.0, 1.0, lo, hi), "4π / 3 = 4.18879020479");

    api("【4】 newton(f, df, x0) / newton_numeric(f, x0)", "optional<double>", "接線で解に近づく（df を書かないなら numeric）", "数回〜数十回");
    SHOW(newton([](double x) { return x * x - 2; }, [](double x) { return 2 * x; }, 1.0), "√2");
    SHOW(newton_numeric([](double x) { return std::pow(x, x) - 10; }, 2.0), "x^x = 10");
    SHOW(newton([](double x) { return x * x + 1; }, [](double x) { return 2 * x; }, 0.0), "傾きが 0 で進めない → nullopt");

    api("【5】 bisect_root(f, lo, hi) / find_roots(f, lo, hi, samples)", "double / vector<double>", "符号が変わる所を二分法 / 区間の解を全部", "100 回 / samples 回 + 解の数 × 100 回");
    SHOW(bisect_root([](double x) { return std::cos(x) - x; }, 0.0, 1.0), "cos x = x");
    SHOW(find_roots([](double x) { return std::sin(x); }, -1.0, 10.0), "0, π, 2π, 3π");

    api("【6】 derivative(f, x)", "double", "f'(x) を中心差分で", "2 回");
    SHOW(derivative([](double x) { return std::exp(x); }, 1.0), "e");
}

}  // namespace numeric_demo

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
    numeric_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    numeric_test::benchmark();    // 速度計測
    numeric_demo::demo();         // 使い方デモ
    return 0;
}
