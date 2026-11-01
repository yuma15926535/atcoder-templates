// =============================================================================
//  整数の共通部分（Int128 / UInt128 とその入出力、整数の平方根、小数の入力）   --- 数学・幾何 第二世代で共通
// -----------------------------------------------------------------------------
//  分類       : 00 基盤（数学の 00_math_base と幾何の 00_geo_base が、これを同梱している）
//  依存       : なし
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  数学と幾何の両方で使う関数を 1 か所にまとめたもの。数学と幾何のテンプレートを一緒に貼っても
//  ガード GG_INT_COMMON で 1 回だけ定義される（どちらを先に貼ってもよい）。
//  このファイルは math_v2 と geometry_v2 の両方に同じものを置く（update_bundles がそれぞれのフォルダで使う）。
//
//  機能一覧                                                                    計算量
//    Int128 / UInt128                  128 ビット整数（= __int128 / unsigned __int128）
//    cout << x / cin >> x              Int128・UInt128 をそのまま出す・読む（範囲外の入力は止める）
//    int128_to_string(v) / uint128_to_string(u)   10 進の文字列に                          O(桁数)
//    parse_int128(s) / parse_uint128(s)           10 進の文字列から（読めない・範囲外なら止める）
//    isqrt(n) / isqrt_ceil(n)          ⌊√n⌋ / ⌈√n⌉（誤差なし。n ≤ 9.2×10^18。n ≤ 0 なら 0）     O(1)
//    parse_fixed(s, k)                 小数の文字列を 10^k 倍した整数に（"1.25", 2 → 125。誤差なし）
//    read_fixed(k)                     小数で与えられる値を 1 つ読んで 10^k 倍の整数に
//      parse_fixed は小数点以下に k 桁を超える 0 でない桁・数として読めない文字・long long に収まらない値で止める
//
//  注意
//    - Main_base の debug(x) は、Int128 をそのままでは表示できない（debug の中の << がこの << を見つけられない）。
//      debug(int128_to_string(x)) と書く
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    intc_test::random_test()  … 128 ビットの愚直計算と突き合わせるランダムテスト
//    intc_test::benchmark()    … 速度計測
//    intc_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_INT_COMMON
#define GG_INT_COMMON
// ---- 整数の共通部分（数学・幾何の第二世代で共通）  Int128 / UInt128 とその入出力、整数の平方根、小数の入力 ----
__extension__ typedef __int128 Int128;            // 128 ビット整数（途中の計算・大きな答え用。cout << / cin >> できる）
__extension__ typedef unsigned __int128 UInt128;  // 符号なし 128 ビット整数

// ---- 整数の平方根 ----
inline long long isqrt(long long n) {  // ⌊√n⌋（n ≤ 0 なら 0）
    if (n <= 0) return 0;
    long long r = (long long)std::sqrt((long double)n);
    while (r > 0 && (Int128)r * r > n) --r;
    while ((Int128)(r + 1) * (r + 1) <= n) ++r;
    return r;
}
inline long long isqrt_ceil(long long n) {  // ⌈√n⌉（n ≤ 0 なら 0）
    if (n <= 0) return 0;
    const long long r = isqrt(n);
    return (Int128)r * r == n ? r : r + 1;
}

// ---- Int128 / UInt128 の文字列化・読み込み・入出力 ----
inline std::string uint128_to_string(UInt128 u) {
    if (u == 0) return "0";
    std::string s;
    while (u) s += char('0' + (int)(u % 10)), u /= 10;
    std::reverse(s.begin(), s.end());
    return s;
}
inline std::string int128_to_string(Int128 v) {
    if (v == 0) return "0";
    const UInt128 u = v < 0 ? (UInt128)0 - (UInt128)v : (UInt128)v;
    return v < 0 ? "-" + uint128_to_string(u) : uint128_to_string(u);
}
// 10 進の文字列を Int128 に（先頭の + / - を受け付ける）。数として読めない・Int128 に収まらないときは止める
inline Int128 parse_int128(const std::string& s) {
    size_t i = 0;
    const bool neg = !s.empty() && s[0] == '-';
    if (!s.empty() && (s[0] == '+' || s[0] == '-')) ++i;
    if (i == s.size()) { std::cerr << "parse_int128: 数字がない: " << s << "\n"; std::abort(); }
    const UInt128 lim = neg ? (UInt128)1 << 127 : ((UInt128)1 << 127) - 1;
    UInt128 u = 0;
    for (; i < s.size(); ++i) {
        const char c = s[i];
        if (c < '0' || c > '9') { std::cerr << "parse_int128: 数として読めない: " << s << "\n"; std::abort(); }
        if (u > (lim - (UInt128)(c - '0')) / 10) { std::cerr << "parse_int128: Int128 に収まらない: " << s << "\n"; std::abort(); }
        u = u * 10 + (UInt128)(c - '0');
    }
    return neg ? (Int128)((UInt128)0 - u) : (Int128)u;
}
inline UInt128 parse_uint128(const std::string& s) {  // 先頭の + を受け付ける。読めない・収まらないときは止める
    size_t i = (!s.empty() && s[0] == '+') ? 1 : 0;
    if (i == s.size()) { std::cerr << "parse_uint128: 数字がない: " << s << "\n"; std::abort(); }
    UInt128 u = 0;
    for (; i < s.size(); ++i) {
        const char c = s[i];
        if (c < '0' || c > '9') { std::cerr << "parse_uint128: 数として読めない: " << s << "\n"; std::abort(); }
        if (u > (~(UInt128)0 - (UInt128)(c - '0')) / 10) { std::cerr << "parse_uint128: UInt128 に収まらない: " << s << "\n"; std::abort(); }
        u = u * 10 + (UInt128)(c - '0');
    }
    return u;
}
inline std::ostream& operator<<(std::ostream& os, Int128 v) { return os << int128_to_string(v); }
inline std::ostream& operator<<(std::ostream& os, UInt128 v) { return os << uint128_to_string(v); }
inline std::istream& operator>>(std::istream& is, Int128& v) {
    std::string s;
    if (is >> s) v = parse_int128(s);
    return is;
}
inline std::istream& operator>>(std::istream& is, UInt128& v) {
    std::string s;
    if (is >> s) v = parse_uint128(s);
    return is;
}

// ---- 小数の入力 ----
// 小数の文字列を 10^k 倍した整数にする（"-1.5", 2 → -150、"3", 2 → 300、".25", 2 → 25）。誤差なし
// 小数点以下に k 桁を超える 0 でない桁があったら止める（黙って切り捨てないため）。long long に収まらないときも止める
inline long long parse_fixed(const std::string& s, int k) {
    size_t i = 0;
    bool neg = false, any = false;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) neg = (s[i] == '-'), ++i;
    Int128 v = 0;
    const Int128 lim = (Int128)9223372036854775807LL + (neg ? 1 : 0);
    int frac = -1;  // 小数点以下の桁数（小数点がまだなら -1）
    for (; i < s.size(); ++i) {
        const char c = s[i];
        if (c == '.' && frac < 0) { frac = 0; continue; }
        if (c < '0' || c > '9') { std::cerr << "parse_fixed: 数として読めない: " << s << "\n"; std::abort(); }
        any = true;
        if (frac >= 0) {
            if (frac == k) {
                if (c != '0') { std::cerr << "parse_fixed: 小数点以下が " << k << " 桁を超えている: " << s << "\n"; std::abort(); }
                continue;
            }
            ++frac;
        }
        v = v * 10 + (c - '0');
        if (v > lim) { std::cerr << "parse_fixed: long long に収まらない: " << s << "\n"; std::abort(); }
    }
    if (!any) { std::cerr << "parse_fixed: 数字がない: " << s << "\n"; std::abort(); }
    for (int f = std::max(frac, 0); f < k; ++f) {
        v *= 10;
        if (v > lim) { std::cerr << "parse_fixed: long long に収まらない: " << s << "\n"; std::abort(); }
    }
    return (long long)(neg ? -v : v);
}
inline long long read_fixed(int k) {
    std::string s;
    std::cin >> s;
    return parse_fixed(s, k);
}
#endif  // GG_INT_COMMON
// ======================= ここまでテンプレート本体 =======================


// cin を文字列から読むようにする（テスト・デモ用）
struct CinFromString {
    std::istringstream ss;
    std::streambuf* old;
    explicit CinFromString(const std::string& s) : ss(s), old(std::cin.rdbuf(ss.rdbuf())) {}
    ~CinFromString() { std::cin.rdbuf(old); }
};

// =============================================================================
//  ランダムテスト
// =============================================================================
namespace intc_test {

#define INTC_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[int_common] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                       \
        }                                                                                       \
    } while (0)

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

[[maybe_unused]] static long long rnd(std::mt19937_64& rng, long long lo, long long hi) {  // [lo, hi] の一様乱数（hi - lo が 2^64 未満）
    const unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return w == 0 ? (long long)rng() : (long long)((unsigned long long)lo + rng() % w);
}
// 大きさがばらばらの long long（境目の値を多めに）
[[maybe_unused]] static long long rnd_ll(std::mt19937_64& rng) {
    static const long long special[] = {0, 1, -1, 2, -2, 3, -3, 7, -7, 10, -10, 2147483647LL, -2147483647LL, -2147483648LL, 2147483648LL,
                                        4294967295LL, 4294967296LL, -4294967296LL, 1000000000LL, -1000000000LL, 1000000007LL, 998244353LL,
                                        3037000499LL, 3037000500LL, -3037000500LL, 1000000000000000000LL, -1000000000000000000LL,
                                        4611686018427387904LL, -4611686018427387904LL, 9223372036854775807LL, -9223372036854775807LL - 1,
                                        9223372036854775806LL, -9223372036854775807LL};
    const int t = (int)(rng() % 8);
    if (t == 0) return special[rng() % (sizeof(special) / sizeof(special[0]))];
    const int bits = (int)(rng() % 64);
    long long v = (long long)(rng() >> (63 - bits) >> 1);
    if (t <= 1) v = (long long)(rng() % 21) - 10;
    return (rng() & 1) ? v : -v;
}

// 愚直：Int128 → 10 進（割り算を使わない別の方法：10 のべきを引いていく）
static std::string brute_to_string(Int128 v) {
    if (v == 0) return "0";
    std::string s;
    UInt128 u = v < 0 ? (UInt128)0 - (UInt128)v : (UInt128)v;
    UInt128 p = 1;
    int digits = 1;
    while (digits < 39 && u / p >= 10) p *= 10, ++digits;
    for (; digits > 0; --digits, p /= 10) {
        int d = 0;
        while (u >= p) u -= p, ++d;
        s += char('0' + d);
    }
    return (v < 0 ? "-" : "") + s;
}
static Int128 rnd128(std::mt19937_64& rng) {
    const int sh = (int)rnd(rng, 0, 64);
    UInt128 u = ((UInt128)rng() << 64) | rng();
    if (sh < 64) u >>= (int)rnd(rng, 0, 127);
    return (Int128)u;
}

static void test_isqrt(std::mt19937_64& rng, int& cases) {
    auto check = [&](long long n) {
        const long long r = isqrt(n), c = isqrt_ceil(n);
        if (n <= 0) {
            INTC_CHECK(r == 0 && c == 0);
            return;
        }
        INTC_CHECK((Int128)r * r <= n && (Int128)(r + 1) * (r + 1) > n);
        INTC_CHECK((Int128)c * c >= n && (Int128)(c - 1) * (c - 1) < n);
        ++cases;
    };
    for (long long n = -3; n <= 3000; ++n) check(n);
    for (int it = 0; it < 20000 * TEST_SCALE; ++it) {
        const int bits = (int)rnd(rng, 1, 63);
        const long long n = (long long)(rng() >> (64 - bits));
        check(n);
        const long long s = isqrt(n);
        for (long long d = -2; d <= 2; ++d) {
            const Int128 sq = (Int128)(s + d) * (s + d);
            if (sq >= 0 && sq <= (Int128)9223372036854775807LL) check((long long)sq), check((long long)sq - 1), check((long long)std::min<Int128>(sq + 1, (Int128)9223372036854775807LL));
        }
    }
    check(9223372036854775807LL);
    check(3037000499LL * 3037000499LL);
    check(3037000499LL * 3037000499LL - 1);
    check(-9223372036854775807LL - 1);
    INTC_CHECK(isqrt(9223372036854775807LL) == 3037000499LL && isqrt_ceil(9223372036854775807LL) == 3037000500LL);
}

static void test_int128_io(std::mt19937_64& rng, int& cases) {
    const Int128 mx = (Int128)(((UInt128)1 << 127) - 1), mn = -mx - 1;
    const UInt128 umx = ~(UInt128)0;
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        Int128 v = rnd128(rng);
        if (it % 50 == 0) v = it % 100 == 0 ? mx - (Int128)rnd(rng, 0, 5) : mn + (Int128)rnd(rng, 0, 5);
        const std::string s = int128_to_string(v);
        INTC_CHECK(s == brute_to_string(v));
        INTC_CHECK(parse_int128(s) == v);
        if (v >= 0) INTC_CHECK(parse_int128("+" + s) == v && parse_uint128(s) == (UInt128)v && uint128_to_string((UInt128)v) == s);
        const UInt128 u = (UInt128)v * 3 + 1;
        INTC_CHECK(parse_uint128(uint128_to_string(u)) == u);
        // ストリーム
        std::ostringstream os;
        os << v << ' ' << u << ' ' << (long long)(v >> 70) << ' ' << -7;
        std::istringstream is(os.str());
        Int128 a = 0;
        UInt128 b = 0;
        long long c = 0, d = 0;
        is >> a >> b >> c >> d;
        INTC_CHECK(a == v && b == u && c == (long long)(v >> 70) && d == -7);
    }
    INTC_CHECK(int128_to_string(mx) == "170141183460469231731687303715884105727");
    INTC_CHECK(int128_to_string(mn) == "-170141183460469231731687303715884105728");
    INTC_CHECK(parse_int128("-170141183460469231731687303715884105728") == mn && parse_int128("170141183460469231731687303715884105727") == mx);
    INTC_CHECK(uint128_to_string(umx) == "340282366920938463463374607431768211455" && parse_uint128("340282366920938463463374607431768211455") == umx);
    INTC_CHECK(int128_to_string(0) == "0" && int128_to_string(-7) == "-7" && parse_int128("-0") == 0 && parse_int128("007") == 7);
    {
        std::ostringstream os;
        os << std::setw(6) << (Int128)42 << '|' << (Int128)-5;  // setw も効く
        INTC_CHECK(os.str() == "    42|-5");
        // long long・int・char などは今までどおり（Int128 の << と取り違えない）
        std::ostringstream os2;
        os2 << 5LL << ' ' << 'a' << ' ' << 3u << ' ' << (short)-2 << ' ' << true;
        INTC_CHECK(os2.str() == "5 a 3 -2 1");
    }
    {
        CinFromString cf("123456789012345678901234567890 -5\n");
        Int128 a, b;
        std::cin >> a >> b;
        INTC_CHECK(a == parse_int128("123456789012345678901234567890") && b == -5);
    }
}

static void test_parse_fixed(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 6000 * TEST_SCALE; ++it, ++cases) {
        const int k = (int)rnd(rng, 0, 7);
        long long p10 = 1;
        for (int i = 0; i < k; ++i) p10 *= 10;
        const long long v = rnd(rng, -(it % 2 ? 1000LL : 1000000000000LL), (it % 2 ? 1000LL : 1000000000000LL));
        const long long a = std::llabs(v);
        std::string ip = std::to_string(a / p10), fp = std::to_string(a % p10);
        while ((int)fp.size() < k) fp = "0" + fp;
        if (k == 0) fp = "";
        const int style = (int)rnd(rng, 0, 6);
        std::string s;
        std::string f2 = fp;
        if (style == 1 || style == 2) {  // 末尾の 0 を消す
            while (!f2.empty() && f2.back() == '0') f2.pop_back();
        }
        if (style == 3) f2 += std::string((size_t)rnd(rng, 1, 3), '0');  // 余分な 0
        std::string i2 = ip;
        if (style == 4) i2 = std::string((size_t)rnd(rng, 1, 2), '0') + i2;  // 先頭の 0
        if (style == 5 && i2 == "0" && !f2.empty()) i2 = "";                  // ".25"
        s = i2;
        if (!f2.empty()) s += "." + f2;
        else if (style == 2) s += ".";  // "5."
        if (s.empty()) s = "0";
        if (v < 0) s = "-" + s;
        else if (style == 6) s = "+" + s;
        INTC_CHECK(parse_fixed(s, k) == v);
    }
    INTC_CHECK(parse_fixed("1.25", 2) == 125 && parse_fixed("-1.5", 2) == -150 && parse_fixed("3", 2) == 300);
    INTC_CHECK(parse_fixed(".25", 2) == 25 && parse_fixed("-0.0001", 4) == -1 && parse_fixed("2.5000", 1) == 25);
    INTC_CHECK(parse_fixed("-0", 3) == 0 && parse_fixed("100000.0000", 4) == 1000000000LL);
    // long long の端
    INTC_CHECK(parse_fixed("9223372036854775807", 0) == 9223372036854775807LL);
    INTC_CHECK(parse_fixed("-9223372036854775808", 0) == -9223372036854775807LL - 1);
    INTC_CHECK(parse_fixed("922337203685477580.7", 1) == 9223372036854775807LL);
    INTC_CHECK(parse_fixed("-92233720368547758.08", 2) == -9223372036854775807LL - 1);
    INTC_CHECK(parse_fixed("0000000000000000000000000001.5", 1) == 15);
    {
        CinFromString cf("1.5 -0.125 7\n");
        INTC_CHECK(read_fixed(2) == 150 && read_fixed(3) == -125 && read_fixed(0) == 7);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_isqrt(rng, cases);
    test_int128_io(rng, cases);
    test_parse_fixed(rng, cases);
    std::cout << "[int_common] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[int_common] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 2000000;
    std::vector<long long> xs(n);
    for (auto& x : xs) x = (long long)(rng() >> 1);
    auto t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (int i = 0; i < n; ++i) s1 += isqrt(xs[i]);
    const long long t_sqrt = ms(t0);
    const int m = 200000;
    std::vector<Int128> vs(m);
    for (auto& v : vs) v = (Int128)xs[rng() % n] * (long long)(rng() >> 2) * ((rng() & 1) ? 1 : -1);
    t0 = std::chrono::steady_clock::now();
    std::ostringstream os;
    for (int i = 0; i < m; ++i) os << vs[i] << '\n';
    const long long t_out = ms(t0);
    const std::string text = os.str();
    t0 = std::chrono::steady_clock::now();
    Int128 s2 = 0;
    {
        std::istringstream is(text);
        Int128 v = 0;
        for (int i = 0; i < m; ++i) is >> v, s2 ^= v;
    }
    const long long t_in = ms(t0);
    std::string fx;
    for (int i = 0; i < m; ++i) fx += std::to_string(xs[i] % 1000000) + "." + std::to_string(xs[i] % 10000) + " ";
    t0 = std::chrono::steady_clock::now();
    long long s3 = 0;
    {
        CinFromString cf(fx);
        for (int i = 0; i < m; ++i) s3 += read_fixed(4);
    }
    const long long t_fx = ms(t0);
    std::cout << "[int_common] benchmark: isqrt × " << n << " " << t_sqrt << " ms, Int128 の出力 × " << m << " " << t_out
              << " ms, 入力 × " << m << " " << t_in << " ms, read_fixed(4) × " << m << " " << t_fx << " ms  (check "
              << (s1 % 7) + (long long)(s2 % 7) + s3 % 7 << ")\n";
}

}  // namespace intc_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace intc_demo {

[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(long double v) {
    if (std::fabs(v) < 5e-13L) v = 0;
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}
[[maybe_unused]] static std::string to_s(double v) { return to_s((long double)v); }
[[maybe_unused]] static std::string to_s(const std::string& s) { return s; }
[[maybe_unused]] static std::string to_s(const char* s) { return s; }
template <class T>
[[maybe_unused]] static std::string to_s(const T& v) {  // << で出せるもの（整数・Int128・modint など）
    std::ostringstream os;
    os << v;
    return os.str();
}
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "]";
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
[[maybe_unused]] static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
[[maybe_unused]] static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
[[maybe_unused]] static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
[[maybe_unused]] static void stmt(const std::string& s) { code(s + ";"); }
template <class V>
[[maybe_unused]] static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 34) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 58) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) intc_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    intc_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("整数の共通部分（int_common） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Int128 / UInt128                       128 ビット整数（cout << x / cin >> x がそのまま使える）\n"
              << "  int128_to_string / uint128_to_string   10 進の文字列に\n"
              << "  parse_int128(s) / parse_uint128(s)     10 進の文字列から（読めない・範囲外なら止める）\n"
              << "  isqrt(n) / isqrt_ceil(n)               ⌊√n⌋ / ⌈√n⌉（誤差なし）\n"
              << "  parse_fixed(s, k) / read_fixed(k)      小数を 10^k 倍の整数に（誤差なし）\n"
              << "  （数学の 00_math_base・幾何の 00_geo_base はこれを同梱している。一緒に貼っても 1 回だけ定義される）\n";

    api("【1】 Int128 / UInt128、cout << x、cin >> x", "Int128 / UInt128",
        "long long（約 9.2×10^18）に収まらない途中の計算・答えに使う。cout でそのまま出せる", "O(1)（入出力は O(桁数)）");
    RUN(Int128 big = (Int128)1000000000000000000LL * 1000000000000000000LL);
    SHOW(big, "10^36（long long なら溢れる）");
    SHOW(-big / 7, "割り算・負の数も普通に");
    SHOW(int128_to_string(big + 1), "文字列にする（debug(...) で見るときなど）");
    SHOW(parse_int128("-123456789012345678901234567890"), "文字列から");
    {
        CinFromString cf("99999999999999999999999\n");
        RUN(Int128 x; std::cin >> x);
        SHOW(x, "cin >> x でそのまま読める（入力が 10^18 を超えるとき）");
    }
    note("範囲：Int128 は約 ±1.7×10^38、UInt128 は 0〜約 3.4×10^38。");
    note("Main_base の debug(x) は Int128 をそのまま表示できないので debug(int128_to_string(x)) と書く。");

    api("【2】 isqrt(n) / isqrt_ceil(n)", "long long", "⌊√n⌋ / ⌈√n⌉ を誤差なしで（sqrt は 10^18 付近で 1 ずれることがある）", "O(1)");
    SHOW(isqrt(99), "√99 = 9.94... の切り捨て");
    SHOW(isqrt_ceil(99), "切り上げ");
    SHOW(isqrt(999999999999999999LL), "10^18 − 1 → 999999999（sqrt だと 1000000000 になりうる）");
    SHOW(isqrt_ceil(1000000000000000000LL), "ちょうど平方数なら同じ値");
    SHOW(isqrt(-5), "0 以下は 0");
    note("平方数か調べる：isqrt(n) * isqrt(n) == n（00_math_base の is_square でも）");

    api("【3】 parse_fixed(s, k) / read_fixed(k)", "long long",
        "小数点以下 k 桁までの小数を 10^k 倍した整数にする（誤差なし）。小数の入力はこれで整数にしてから計算する", "O(桁数)");
    SHOW(parse_fixed("1.25", 2), "1.25 × 10^2");
    SHOW(parse_fixed("-0.5", 4), "負の数も");
    SHOW(parse_fixed("3", 4), "小数点がなくてもよい");
    SHOW(parse_fixed(".25", 2), "整数部がなくてもよい");
    {
        // ABC169 C: A × B（B は小数点以下 2 桁）の整数部分
        CinFromString cf("198 1.10\n");
        RUN(long long A; std::cin >> A; long long B100 = read_fixed(2));
        SHOW(A * B100 / 100, "ABC169 C: 198 × 1.10 = 217.8 → 217（double で掛けると誤差で 217 にならないことがある）");
    }
    note("小数点以下に k 桁を超える 0 でない桁・long long に収まらない値は止まる（黙って切り捨てない）。");
}

}  // namespace intc_demo

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
    intc_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    intc_test::benchmark();    // 速度計測
    intc_demo::demo();         // 使い方デモ
    return 0;
}
