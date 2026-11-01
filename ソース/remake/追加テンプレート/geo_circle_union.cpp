// =============================================================================
//  円の和集合の面積・周の長さ・k 個以上の円に覆われた面積（円周を弧に分けてグリーンの定理）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  図形の面積は、周を反時計回りに 1 周する線積分 (1/2)∮(x dy − y dx)（グリーンの定理）。和集合の周は
//  「ほかのどの円にも覆われていない弧」の集まりなので、円ごとに覆われる角度の区間を集め、残った弧の線積分を足す。
//  弧を覆う円の数 c ごとに足すと「c + 1 個以上の円に覆われた部分」の面積も同時に出る
//    例: 「スプリンクラー N 個が水をまく範囲の面積」「ちょうど k 個の電波が届く範囲の面積」「円を重ねた図形の外周」
//
//  機能一覧（円は Disk {x, y, r}（long double）。N 個）                                         計算量
//    circle_union_area(cs)               和集合の面積（どれかの円に入る部分）                                O(N^2 log N)
//    circle_union_perimeter(cs)          和集合の周の長さ（穴の周も含む）                                      O(N^2 log N)
//    circle_cover_areas(cs)              res[k] = k 個以上の円に覆われた面積（k = 1〜N。res[0] は使わず 0）    O(N^2 log N)
//    circle_cover_exact(cs)              res[k] = ちょうど k 個の円に覆われた面積（k = 1〜N。res[0] は 0）     O(N^2 log N)
//
//  注意: 同じ円が何個あっても正しく数える（cover なら重ねた数だけ）。半径 0 以下の円は無視する
//        接する・含むの判定は 1e-12 × (座標と半径の大きさ) の誤差を許す。答えの誤差は 面積 × 1e-15 くらい
//        2 円の共通部分・円と多角形の共通部分は geometry_v2/05_circle_area.cpp
//  速さの目安: N = 2000（ランダム）の circle_union_area 約 0.4 秒（circle_cover_areas もほぼ同じ）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    cunion_test::random_test()  … x で輪切りにして長さを数値積分する愚直解との突き合わせ
//    cunion_test::benchmark()    … 速度計測
//    cunion_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
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
struct Disk {
    long double x, y, r;
};

namespace cunion_internal {
using LD = long double;
// 円 i の弧を、ほかの円に覆われている数 c ごとに f(i, c, a, b)（角度 a < b、反時計回り）で渡す
template <class F>
void for_each_arc(const std::vector<Disk>& all, F f) {
    std::vector<Disk> cs;
    for (const auto& c : all)
        if (c.r > 0) cs.push_back(c);
    const int n = (int)cs.size();
    const LD PI = std::acos((LD)-1), TAU = 2 * PI;
    std::vector<std::pair<LD, int>> ev;
    for (int i = 0; i < n; ++i) {
        const Disk& ci = cs[i];
        int base = 0;  // 円 i を丸ごと覆う円の数
        ev.clear();
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            const Disk& cj = cs[j];
            const LD dx = cj.x - ci.x, dy = cj.y - ci.y, d = std::hypot(dx, dy);
            const LD eps = 1e-12L * std::max({(LD)1, std::abs(ci.x), std::abs(ci.y), std::abs(cj.x), std::abs(cj.y), ci.r, cj.r});
            if (d <= eps && std::abs(ci.r - cj.r) <= eps) {  // 同じ円：番号の小さい方が上に重なっているとみなす
                if (j < i) ++base;
                continue;
            }
            if (d + ci.r <= cj.r + eps) {  // j が i を含む（内側で接するのも）
                ++base;
                continue;
            }
            if (d + cj.r <= ci.r + eps || d >= ci.r + cj.r - eps) continue;  // i が j を含む・離れている（接する）
            const LD phi = std::atan2(dy, dx);
            const LD ca = std::clamp((ci.r * ci.r + d * d - cj.r * cj.r) / (2 * ci.r * d), (LD)-1, (LD)1);
            const LD al = std::acos(ca);
            LD lo = phi - al, hi = phi + al;
            if (lo < 0) lo += TAU, hi += TAU;
            if (hi > TAU) {
                ev.push_back({lo, +1}), ev.push_back({TAU, -1});
                ev.push_back({0, +1}), ev.push_back({hi - TAU, -1});
            } else {
                ev.push_back({lo, +1}), ev.push_back({hi, -1});
            }
        }
        std::sort(ev.begin(), ev.end());
        LD prev = 0;
        int cur = base;
        for (const auto& [ang, dlt] : ev) {
            if (ang > prev) f(ci, cur, prev, ang);
            prev = std::max(prev, ang);
            cur += dlt;
        }
        if (TAU > prev) f(ci, cur, prev, TAU);
    }
}
// 円 c の弧 [a, b] に沿った (1/2)∫(x dy − y dx)
inline LD green(const Disk& c, LD a, LD b) {
    return (c.r * c.r * (b - a) + c.x * c.r * (std::sin(b) - std::sin(a)) - c.y * c.r * (std::cos(b) - std::cos(a))) / 2;
}
}  // namespace cunion_internal

inline long double circle_union_area(const std::vector<Disk>& cs) {
    long double s = 0;
    cunion_internal::for_each_arc(cs, [&](const Disk& c, int cov, long double a, long double b) {
        if (cov == 0) s += cunion_internal::green(c, a, b);
    });
    return s;
}
inline long double circle_union_perimeter(const std::vector<Disk>& cs) {
    long double s = 0;
    cunion_internal::for_each_arc(cs, [&](const Disk& c, int cov, long double a, long double b) {
        if (cov == 0) s += c.r * (b - a);
    });
    return s;
}
inline std::vector<long double> circle_cover_areas(const std::vector<Disk>& cs) {
    std::vector<long double> res(cs.size() + 2, 0);
    cunion_internal::for_each_arc(cs, [&](const Disk& c, int cov, long double a, long double b) { res[cov + 1] += cunion_internal::green(c, a, b); });
    res.resize(cs.size() + 1);
    return res;
}
inline std::vector<long double> circle_cover_exact(const std::vector<Disk>& cs) {
    auto at_least = circle_cover_areas(cs);
    std::vector<long double> res(at_least.size(), 0);
    for (size_t k = 1; k < at_least.size(); ++k) res[k] = at_least[k] - (k + 1 < at_least.size() ? at_least[k + 1] : 0);
    return res;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace cunion_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define CU_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[geo_circle_union] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
using LD = long double;

// 愚直：x で輪切りにして「k 個以上に覆われた y の長さ」を数値積分。
// 円の左右の端・2 円の交点の x で区間を分け、各区間は x = c + h sin θ と置いて端の √ を消してからシンプソン則
static std::vector<LD> brute(const std::vector<Disk>& cs) {
    const int n = (int)cs.size();
    std::vector<LD> res(n + 1, 0);
    std::vector<LD> xs;
    for (const auto& c : cs)
        if (c.r > 0) xs.push_back(c.x - c.r), xs.push_back(c.x + c.r);
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            const auto &a = cs[i], &b = cs[j];
            const LD d = std::hypot(b.x - a.x, b.y - a.y);
            if (d == 0 || d > a.r + b.r || d < std::abs(a.r - b.r)) continue;
            const LD t = (a.r * a.r - b.r * b.r + d * d) / (2 * d), h = std::sqrt(std::max((LD)0, a.r * a.r - t * t));
            const LD mx = a.x + (b.x - a.x) * t / d;
            xs.push_back(mx + h * (b.y - a.y) / d), xs.push_back(mx - h * (b.y - a.y) / d);
        }
    std::sort(xs.begin(), xs.end());
    std::vector<std::pair<LD, int>> ev;
    auto slice = [&](LD x, std::vector<LD>& acc, LD w) {
        ev.clear();
        for (const auto& c : cs) {
            const LD dx = x - c.x;
            if (c.r <= 0 || std::abs(dx) >= c.r) continue;
            const LD s = std::sqrt(c.r * c.r - dx * dx);
            ev.push_back({c.y - s, +1}), ev.push_back({c.y + s, -1});
        }
        std::sort(ev.begin(), ev.end());
        int cur = 0;
        for (size_t i = 0; i < ev.size(); ++i) {
            if (i > 0 && cur > 0) {
                const LD len = ev[i].first - ev[i - 1].first;
                for (int k = 1; k <= cur; ++k) acc[k] += w * len;
            }
            cur += ev[i].second;
        }
    };
    const LD PI = std::acos((LD)-1);
    const int M = 200;  // 偶数
    for (size_t p = 0; p + 1 < xs.size(); ++p) {
        const LD l = xs[p], r = xs[p + 1];
        if (r - l < 1e-12L) continue;
        const LD c = (l + r) / 2, h = (r - l) / 2, dt = PI / M;
        for (int i = 0; i <= M; ++i) {
            const LD th = -PI / 2 + dt * i;
            const LD w = (i == 0 || i == M ? 1 : (i & 1 ? 4 : 2)) * dt / 3 * h * std::cos(th);
            slice(c + h * std::sin(th), res, w);
        }
    }
    return res;
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    std::uniform_real_distribution<double> ud(-5.0, 5.0), ur(0.1, 4.0);
    long long cases = 0;
    const LD PI = std::acos((LD)-1);
    for (int it = 0; it < 250 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 6);
        std::vector<Disk> cs(n);
        const bool integral = it % 2 == 0;  // 整数の座標なら、同じ円・接する円・含む円がよく出る
        for (auto& c : cs) {
            if (integral) c = {(LD)rnd(rng, -3, 3), (LD)rnd(rng, -3, 3), (LD)rnd(rng, 1, 4)};
            else c = {(LD)ud(rng), (LD)ud(rng), (LD)ur(rng)};
        }
        if (n >= 2 && it % 7 == 0) cs[1] = cs[0];  // 同じ円
        if (n >= 1 && it % 11 == 0) cs[0].r = 0;   // 半径 0
        const auto br = brute(cs);
        const auto ca = circle_cover_areas(cs);
        const auto ex = circle_cover_exact(cs);
        const LD u = circle_union_area(cs);
        CU_CHECK((int)ca.size() == n + 1 && (int)ex.size() == n + 1);
        const LD tol = 1e-7L * std::max((LD)1, br.size() > 1 ? br[1] : 0);
        if (n >= 1) CU_CHECK(std::abs(u - br[1]) < tol);
        else CU_CHECK(u == 0);
        for (int k = 1; k <= n; ++k) {
            CU_CHECK(std::abs(ca[k] - br[k]) < tol);
            CU_CHECK(std::abs(ex[k] - (br[k] - (k < n ? br[k + 1] : 0))) < 2 * tol);
        }
        // 周の長さ：全部の半径を t だけ大きくしたときの面積の増え方（t → +0 の微分）= 周の長さ
        if (!integral) {
            const LD h = 1e-6L;
            auto grow = [&](LD t) {
                auto c2 = cs;
                for (auto& c : c2)
                    if (c.r > 0) c.r += t;
                return circle_union_area(c2);
            };
            const LD der = (-3 * grow(0) + 4 * grow(h) - grow(2 * h)) / (2 * h);
            CU_CHECK(std::abs(circle_union_perimeter(cs) - der) < 1e-4L * std::max((LD)1, der));
        }
        cases += 2 + 2 * n;
    }
    // 決まった形
    {
        CU_CHECK(std::abs(circle_union_area({{0, 0, 1}, {0, 0, 1}, {0, 0, 1}}) - PI) < 1e-15L);
        CU_CHECK(std::abs(circle_union_perimeter({{0, 0, 1}, {0, 0, 1}}) - 2 * PI) < 1e-15L);
        CU_CHECK(std::abs(circle_union_area({{0, 0, 3}, {1, 0, 1}}) - 9 * PI) < 1e-14L);                          // 内側で接する
        CU_CHECK(std::abs(circle_union_area({{0, 0, 1}, {2, 0, 1}}) - 2 * PI) < 1e-14L);                          // 外側で接する
        CU_CHECK(std::abs(circle_union_perimeter({{0, 0, 3}, {1, 0, 1}, {10, 0, 2}}) - 10 * PI) < 1e-13L);        // 含む・離れている
        const auto e = circle_cover_exact({{0, 0, 2}, {0, 0, 2}, {0, 0, 1}});
        CU_CHECK(std::abs(e[1]) < 1e-14L && std::abs(e[2] - 3 * PI) < 1e-13L && std::abs(e[3] - PI) < 1e-13L);  // 2 枚重ね + 小さい円
        // 2 つの単位円（中心の距離 1）：レンズ 2π/3 − √3/2
        const LD lens = 2 * PI / 3 - std::sqrt((LD)3) / 2;
        CU_CHECK(std::abs(circle_union_area({{0, 0, 1}, {1, 0, 1}}) - (2 * PI - lens)) < 1e-15L);
        CU_CHECK(std::abs(circle_cover_areas({{0, 0, 1}, {1, 0, 1}})[2] - lens) < 1e-15L);
        CU_CHECK(circle_union_area({}) == 0);
        cases += 9;
    }
    std::cout << "[geo_circle_union] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[geo_circle_union] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 300 : 2000;
    std::uniform_real_distribution<double> up(0.0, 1000.0), ur(1.0, 50.0);
    std::vector<Disk> cs(N);
    for (auto& c : cs) c = {(LD)up(rng), (LD)up(rng), (LD)ur(rng)};
    const auto t0 = std::chrono::steady_clock::now();
    const LD a = circle_union_area(cs);
    const auto t1 = elapsed_ms(t0);
    const auto t2s = std::chrono::steady_clock::now();
    const auto ca = circle_cover_areas(cs);
    std::cout << "[geo_circle_union] benchmark (N = " << N << "): circle_union_area " << t1 << " ms / circle_cover_areas " << elapsed_ms(t2s) << " ms  (面積 " << std::fixed
              << std::setprecision(3) << (double)a << "、2 個以上 " << (double)ca[2] << ")\n"
              << std::defaultfloat;
}

}  // namespace cunion_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace cunion_demo {

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
#define SHOW(expr, comment) cunion_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    cunion_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("円の和集合（geo_circle_union） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  circle_union_area(cs) / circle_union_perimeter(cs) / circle_cover_areas(cs) / circle_cover_exact(cs)   （cs は vector<Disk {x, y, r}>）\n";

    api("【1】 circle_union_area(cs) / circle_union_perimeter(cs)", "long double", "どれかの円に入る部分の面積 / その周の長さ", "O(N^2 log N)");
    RUN(const std::vector<Disk> cs = {{0, 0, 1}, {1, 0, 1}});
    SHOW((double)circle_union_area(cs), "中心の距離 1 の単位円 2 つ：2π − (2π/3 − √3/2) = 5.05481");
    SHOW((double)circle_union_perimeter(cs), "円周 2 つから 120° ずつ欠けた分：2 × 2π × (2/3) = 8.37758");
    SHOW((double)circle_union_area({{0, 0, 3}, {1, 0, 1}, {0, 0, 3}}), "含まれる円・同じ円は数えない：9π");

    api("【2】 circle_cover_areas(cs) / circle_cover_exact(cs)", "vector<long double>", "k 個以上 / ちょうど k 個の円に覆われた面積（添字 k = 1〜N）", "O(N^2 log N)");
    RUN(const std::vector<Disk> ts = {{0, 0, 2}, {2, 0, 2}, {1, 1.7, 2}});
    RUN(const auto at = circle_cover_areas(ts));
    SHOW((double)at[1], "3 つの円の和集合");
    SHOW((double)at[2], "2 つ以上が重なる部分");
    SHOW((double)at[3], "3 つとも重なる部分");
    RUN(const auto ex = circle_cover_exact(ts));
    SHOW((double)ex[2], "ちょうど 2 つ = at[2] − at[3]");
}

}  // namespace cunion_demo

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
    cunion_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    cunion_test::benchmark();    // 速度計測
    cunion_demo::demo();         // 使い方デモ
    return 0;
}
