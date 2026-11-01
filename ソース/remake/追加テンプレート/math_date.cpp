// =============================================================================
//  日付・曜日（うるう年・月の日数・日付 ↔ 通し日数・日数の差・曜日・何日後）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  グレゴリオ暦（4 で割り切れる年はうるう年。ただし 100 で割り切れて 400 で割り切れない年は平年）を、過去にも未来にも
//  そのまま延ばして使う。日付を「1970/01/01 からの日数」に直すと、差・何日後・曜日が足し算・割り算でできる
//    例: 「2 つの日付の間の日数」「N 日後の日付」「その日は何曜日か」「次の日付」「ある日付より前か」
//
//  機能一覧（年 y・月 m（1〜12）・日 d（1〜））                                                計算量
//    is_leap(y) / days_in_month(y, m) / days_in_year(y)    うるう年か / その月の日数 / その年の日数              O(1)
//    is_valid_date(y, m, d)              正しい日付か                                                      O(1)
//    date_to_days(y, m, d)               1970/01/01 を 0 とした通し日数（前は負）                          O(1)
//    days_to_date(z)                     通し日数 → Date {y, m, d}                                          O(1)
//    days_between(a, b)                  b − a（日数。Date どうし）                                         O(1)
//    weekday(y, m, d)                    曜日 0〜6（0 = 日曜, 1 = 月曜, …, 6 = 土曜）                        O(1)
//    add_days(date, k)                   k 日後（k < 0 なら前）の Date                                      O(1)
//    leap_years_up_to(y)                 1〜y 年のうるう年の数（y ≥ 0）                                     O(1)
//    Date {y, m, d}                      比較（< == など）、to_days()、next()、weekday()、to_string()（"2026/10/31"）、
//                                        parse_date("2026/10/31")（区切りは / か -）
//
//  注意: 年は long long（±10^15 年くらいまで）。0 年・負の年も天文学の数え方（0 年 = 紀元前 1 年）でそのまま計算する
//  速さの目安: date_to_days と days_to_date を 10^7 回ずつ 約 0.15 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    date_test::random_test()  … 1 日ずつ進める愚直解との突き合わせ（1 年〜2800 年の毎日）
//    date_test::benchmark()    … 速度計測
//    date_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace date_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
inline long long fdiv(long long a, long long b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }
}  // namespace date_internal

inline bool is_leap(long long y) { return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0; }
inline int days_in_month(long long y, int m) {
    static const int D[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (m < 1 || m > 12) date_internal::fail("days_in_month", "m は 1〜12");
    return m == 2 && is_leap(y) ? 29 : D[m];
}
inline int days_in_year(long long y) { return is_leap(y) ? 366 : 365; }
inline bool is_valid_date(long long y, int m, int d) { return 1 <= m && m <= 12 && 1 <= d && d <= days_in_month(y, m); }
// 1970/01/01 を 0 とする（Howard Hinnant の days_from_civil）
inline long long date_to_days(long long y, int m, int d) {
    if (!is_valid_date(y, m, d)) date_internal::fail("date_to_days", "正しくない日付");
    y -= m <= 2;
    const long long era = date_internal::fdiv(y, 400);
    const long long yoe = y - era * 400;                                     // [0, 399]
    const long long doy = (153LL * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;  // [0, 365]（3 月 1 日から）
    const long long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;             // [0, 146096]
    return era * 146097 + doe - 719468;
}
struct Date {
    long long y;
    int m, d;
    long long to_days() const { return date_to_days(y, m, d); }
    int weekday() const;
    Date next() const;
    std::string to_string() const {
        char b[64];
        std::snprintf(b, sizeof b, "%04lld/%02d/%02d", y, m, d);
        return b;
    }
    bool operator<(const Date& o) const { return std::tie(y, m, d) < std::tie(o.y, o.m, o.d); }
    bool operator>(const Date& o) const { return o < *this; }
    bool operator<=(const Date& o) const { return !(o < *this); }
    bool operator>=(const Date& o) const { return !(*this < o); }
    bool operator==(const Date& o) const { return y == o.y && m == o.m && d == o.d; }
    bool operator!=(const Date& o) const { return !(*this == o); }
};
inline Date days_to_date(long long z) {
    z += 719468;
    const long long era = date_internal::fdiv(z, 146097);
    const long long doe = z - era * 146097;
    const long long yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    const long long doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const long long mp = (5 * doy + 2) / 153;
    const int d = (int)(doy - (153 * mp + 2) / 5 + 1);
    const int m = (int)(mp < 10 ? mp + 3 : mp - 9);
    const long long y = yoe + era * 400 + (m <= 2);
    return {y, m, d};
}
inline long long days_between(const Date& a, const Date& b) { return b.to_days() - a.to_days(); }
inline int weekday(long long y, int m, int d) {
    const long long z = date_to_days(y, m, d);  // 1970/01/01 は木曜（4）
    return (int)(((z + 4) % 7 + 7) % 7);
}
inline Date add_days(const Date& a, long long k) { return days_to_date(a.to_days() + k); }
inline long long leap_years_up_to(long long y) {
    if (y < 0) date_internal::fail("leap_years_up_to", "y は 0 以上");
    return y / 4 - y / 100 + y / 400;
}
inline int Date::weekday() const { return ::weekday(y, m, d); }
inline Date Date::next() const { return add_days(*this, 1); }
inline Date parse_date(const std::string& s) {  // "2026/10/31" や "2026-10-31"
    long long y = 0;
    int m = 0, d = 0;
    char c1 = 0, c2 = 0;
    std::istringstream is(s);
    if (!(is >> y >> c1 >> m >> c2 >> d) || (c1 != '/' && c1 != '-') || c1 != c2 || !is_valid_date(y, m, d))
        date_internal::fail("parse_date", "\"YYYY/MM/DD\" の形の正しい日付にする");
    return {y, m, d};
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace date_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define DT_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[math_date] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    long long cases = 0;
    // 1 年 1 月 1 日から 1 日ずつ進める（月の日数を愚直に数える）
    Date cur{1, 1, 1};
    ll z = date_to_days(1, 1, 1);
    int wd = weekday(1, 1, 1);
    ll leaps = 0;
    for (ll y = 1; y <= 2800; ++y) {
        const bool lp = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
        DT_CHECK(is_leap(y) == lp && days_in_year(y) == (lp ? 366 : 365));
        leaps += lp;
        DT_CHECK(leap_years_up_to(y) == leaps);
        for (int m = 1; m <= 12; ++m) {
            const int dim = m == 2 ? (lp ? 29 : 28) : (m == 4 || m == 6 || m == 9 || m == 11 ? 30 : 31);
            DT_CHECK(days_in_month(y, m) == dim);
            DT_CHECK(!is_valid_date(y, m, 0) && !is_valid_date(y, m, dim + 1));
            for (int d = 1; d <= dim; ++d) {
                DT_CHECK(is_valid_date(y, m, d));
                DT_CHECK((cur == Date{y, m, d}));
                DT_CHECK(date_to_days(y, m, d) == z);
                DT_CHECK(days_to_date(z) == cur);
                DT_CHECK(weekday(y, m, d) == wd && cur.weekday() == wd);
                const Date nx = cur.next();
                DT_CHECK(cur < nx && nx.to_days() == z + 1);
                cur = nx, ++z, wd = (wd + 1) % 7;
                ++cases;
            }
        }
    }
    // 知っている曜日：1970/01/01 は木曜、2000/01/01 は土曜、2026/10/31 は土曜
    DT_CHECK(date_to_days(1970, 1, 1) == 0 && weekday(1970, 1, 1) == 4);
    DT_CHECK(weekday(2000, 1, 1) == 6 && weekday(2026, 10, 31) == 6);
    // 遠い年・負の年でも往復できる
    std::mt19937_64 rng(1);
    for (int i = 0; i < 100000; ++i) {
        const ll zz = rnd(rng, -400000000000000LL, 400000000000000LL);
        const Date a = days_to_date(zz);
        DT_CHECK(is_valid_date(a.y, a.m, a.d) && a.to_days() == zz);
        const ll k = rnd(rng, -1000000, 1000000);
        DT_CHECK(days_between(a, add_days(a, k)) == k);
        ++cases;
    }
    DT_CHECK((parse_date("2026/10/31") == Date{2026, 10, 31}) && (parse_date("1999-02-28") == Date{1999, 2, 28}));
    DT_CHECK(Date{2026, 1, 5}.to_string() == "2026/01/05");
    std::cout << "[math_date] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[math_date] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 1000000 : 10000000;
    const auto t0 = std::chrono::steady_clock::now();
    ll chk = 0;
    for (int i = 0; i < N; ++i) {
        const Date a = days_to_date(i * 37LL - 100000000);
        chk += a.d + date_to_days(a.y, a.m, a.d);
    }
    std::cout << "[math_date] benchmark: days_to_date + date_to_days を " << N << " 回 " << elapsed_ms(t0) << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace date_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace date_demo {

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
#define SHOW(expr, comment) date_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    date_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static std::string to_s(const Date& d) { return d.to_string(); }

static void demo() {
    banner("日付・曜日（math_date） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  is_leap / days_in_month / days_in_year / is_valid_date / date_to_days / days_to_date / days_between\n"
              << "  weekday / add_days / leap_years_up_to / parse_date      Date {y, m, d}: to_days / next / weekday / to_string / 比較\n";

    api("【1】 is_leap(y) / days_in_month(y, m) / is_valid_date(y, m, d)", "bool / int / bool", "うるう年か / 月の日数 / 正しい日付か", "O(1)");
    SHOW(is_leap(2024), "");
    SHOW(is_leap(1900), "100 で割り切れて 400 で割り切れない");
    SHOW(is_leap(2000), "400 で割り切れる");
    SHOW(days_in_month(2024, 2), "");
    SHOW(is_valid_date(2023, 2, 29), "");

    api("【2】 date_to_days(y, m, d) / days_to_date(z) / days_between(a, b)", "long long / Date / long long",
        "1970/01/01 = 0 の通し日数 / その逆 / b − a の日数", "O(1)");
    SHOW(date_to_days(2026, 10, 31), "");
    SHOW(to_s(days_to_date(20000)), "1970/01/01 の 20000 日後");
    SHOW((days_between(Date{2026, 1, 1}, Date{2026, 12, 25})), "");

    api("【3】 weekday(y, m, d) / add_days(date, k)", "int（0 = 日 … 6 = 土）/ Date", "曜日 / k 日後", "O(1)");
    SHOW(weekday(2026, 10, 31), "6 = 土曜");
    SHOW((to_s(add_days(Date{2024, 2, 28}, 1))), "うるう年");
    SHOW((to_s(add_days(Date{2026, 1, 1}, -1))), "前の日");

    api("【4】 Date の比較・next()・to_string() / parse_date(s)", "—", "日付の大小・次の日・\"YYYY/MM/DD\" の文字列 / 文字列から", "O(1)");
    RUN(const Date a = parse_date("2019/04/30"));
    SHOW((a <= Date{2019, 4, 30}), "「2019/04/30 以前か」の判定");
    SHOW(to_s(a.next()), "");
    SHOW(leap_years_up_to(2026), "1〜2026 年のうるう年の数");
}

}  // namespace date_demo

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
    date_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    date_test::benchmark();    // 速度計測
    date_demo::demo();         // 使い方デモ
    return 0;
}
