// =============================================================================
//  べき乗和・ラグランジュ補間・次数の小さい多項式（Σ i^k、f(0)〜f(d) から f(x)、多項式の割り算）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 06 和・数列
//  依存       : 00_int_common, 00_math_base, 00_modint（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  使い方
//      using mint = ModInt107;
//      cout << sum_powers<mint>(N + 1, K) << '\n';           // 1^K + 2^K + … + N^K（K ≥ 1）
//      cout << lagrange_interpolate(ys, T) << '\n';           // ys = {f(0), …, f(d)} のとき f(T)
//
//  機能一覧                                                                    計算量
//    sum_powers<mint>(n, k)           Σ_{i=0}^{n−1} i^k（0^0 = 1。n ≤ 10^18。mod は k + 1 より大きい素数）    O(k)
//    sum_powers<mint>(l, r, k) / sum_powers_closed<mint>(L, R, k)    Σ_{i=l}^{r−1} / Σ_{i=L}^{R}
//    arith_geom_sum(r, n)             Σ_{i=0}^{n−1} i·r^i（r は mint。mod は何でもよい）                     O(log n)
//    lagrange_interpolate(ys, x)      f(0)〜f(d) が ys の d 次以下の多項式の f(x)（x は 10^18 でも負でも）     O(d)
//    lagrange_interpolate_at(xs, ys, x)   点 (xs[i], ys[i]) を通る多項式の f(x)（xs はばらばらの値）         O(d^2)
//    poly_from_points(xs, ys)         その多項式の係数（低い次数から）                                       O(d^2)
//   次数の小さい多項式（係数は低い次数から。T は mint・long long・double など）
//    poly_eval(c, x)                  c(x)                                                                   O(n)
//    poly_mul_naive(a, b)             a × b                                                                  O(|a||b|)
//    poly_divmod(a, b)                PolyDivMod{q, r}：a = q b + r（T が整数なら割り切れないと止める）         O(|a||b|)
//  （等差・等比数列の和は 00_math_base、次数の大きい多項式（畳み込み・FPS）は 07_convolution・07_fps）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    ser_test::random_test()  … 1 つずつ足す愚直解・係数からの計算との突き合わせ
//    ser_test::benchmark()    … 速度計測
//    ser_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
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
#ifndef GG_MATH_BASE
#define GG_MATH_BASE
// ---- 数学基盤（数学第二世代）  割り算・倍数・累乗根・累乗と溢れ・gcd / lcm・64 ビットの mod・数列の和・分数の比較 ----
template <class T>
struct Maybe {  // 「ないことがある値」：ok が false のとき value は使わない
    bool ok = false;
    T value{};
};

namespace math_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
inline Int128 floor_div128(Int128 a, Int128 b) {
    const Int128 q = a / b;
    return q - ((a % b != 0) && ((a < 0) != (b < 0)));
}
inline Int128 ceil_div128(Int128 a, Int128 b) {
    const Int128 q = a / b;
    return q + ((a % b != 0) && ((a < 0) == (b < 0)));
}
inline unsigned long long uabs(long long a) { return a < 0 ? 0ULL - (unsigned long long)a : (unsigned long long)a; }
inline unsigned long long ugcd(unsigned long long a, unsigned long long b) {
    while (b) {
        const unsigned long long t = a % b;
        a = b, b = t;
    }
    return a;
}
// 列の要素（int・unsigned・Int128 など）を long long に。収まらなければ止める（黙って切り詰めない）
template <class T>
long long to_ll_checked(const T& x, const char* func) {
    if constexpr (std::is_same<T, Int128>::value) {
        if (x > (Int128)9223372036854775807LL || x < -(Int128)9223372036854775807LL - 1) fail(func, "値が long long に収まらない");
    } else if constexpr (std::is_same<T, UInt128>::value) {
        if (x > (UInt128)9223372036854775807ULL) fail(func, "値が long long に収まらない");
    } else if constexpr (std::is_unsigned<T>::value && sizeof(T) >= sizeof(long long)) {
        if (x > (T)9223372036854775807ULL) fail(func, "値が long long に収まらない");
    }
    return (long long)x;
}
// y^k ≤ n か（y ≥ 0, k ≥ 1。n < 0 なら false。溢れない）
inline bool pow_le(long long y, int k, long long n) {
    if (n < 0) return false;
    unsigned long long r = 1;
    for (int i = 0; i < k; ++i)
        if (__builtin_mul_overflow(r, (unsigned long long)y, &r) || r > (unsigned long long)n) return false;
    return true;
}
}  // namespace math_internal

// ---- 割り算・余り（負の数でも数学の定義どおり。b ≠ 0。LLONG_MIN / -1 は答えが収まらない） ----
inline long long floor_div(long long a, long long b) {  // ⌊a / b⌋
    const long long q = a / b;
    return q - ((a % b != 0) && ((a ^ b) < 0));
}
inline long long ceil_div(long long a, long long b) {  // ⌈a / b⌉
    const long long q = a / b;
    return q + ((a % b != 0) && ((a ^ b) >= 0));
}
inline long long mod_floor(long long a, long long m) {  // a を m（> 0）で割った余りを [0, m) で（-7 mod 3 = 2）。m < 0 なら |m| で
    long long r = a % m;
    if (r < 0) r += (m < 0 ? -m : m);
    return r;
}
inline long long round_div(long long a, long long b) {  // a / b を四捨五入（ちょうど .5 は大きい方へ。-2.5 → -2）
    Int128 x = a, y = b;
    if (y < 0) x = -x, y = -y;
    return (long long)math_internal::floor_div128(2 * x + y, 2 * y);
}
inline long long ceil_multiple(long long x, long long d) {  // x 以上で最小の d（> 0）の倍数。long long に収まらなければ止める
    if (d <= 0) math_internal::fail("ceil_multiple", "d > 0 で使う");
    const Int128 r = (Int128)ceil_div(x, d) * d;
    if (r > (Int128)9223372036854775807LL) math_internal::fail("ceil_multiple", "long long に収まらない");
    return (long long)r;
}
inline long long floor_multiple(long long x, long long d) {  // x 以下で最大の d（> 0）の倍数
    if (d <= 0) math_internal::fail("floor_multiple", "d > 0 で使う");
    const Int128 r = (Int128)floor_div(x, d) * d;
    if (r < -(Int128)9223372036854775807LL - 1) math_internal::fail("floor_multiple", "long long に収まらない");
    return (long long)r;
}
inline long long ceil_congruent(long long x, long long a, long long m) {  // x 以上で最小の y ≡ a (mod m)（m > 0）
    if (m <= 0) math_internal::fail("ceil_congruent", "m > 0 で使う");
    const Int128 r = (Int128)x + (Int128)mod_floor((long long)(((Int128)a - x) % m), m);
    if (r > (Int128)9223372036854775807LL) math_internal::fail("ceil_congruent", "long long に収まらない");
    return (long long)r;
}
inline long long floor_congruent(long long x, long long a, long long m) {  // x 以下で最大の y ≡ a (mod m)（m > 0）
    if (m <= 0) math_internal::fail("floor_congruent", "m > 0 で使う");
    const Int128 r = (Int128)x - (Int128)mod_floor((long long)(((Int128)x - a) % m), m);
    if (r < -(Int128)9223372036854775807LL - 1) math_internal::fail("floor_congruent", "long long に収まらない");
    return (long long)r;
}
// [l, r) にある d（> 0）の倍数の数（l ≥ r なら 0）。閉区間 [L, R] は _closed
inline long long count_multiples(long long l, long long r, long long d) {
    if (d <= 0) math_internal::fail("count_multiples", "d > 0 で使う");
    return l >= r ? 0 : ceil_div(r, d) - ceil_div(l, d);
}
inline long long count_multiples_closed(long long L, long long R, long long d) {
    if (d <= 0) math_internal::fail("count_multiples_closed", "d > 0 で使う");
    return L > R ? 0 : floor_div(R, d) - ceil_div(L, d) + 1;
}
// [l, r) にある x ≡ a (mod m)（m > 0）の数。閉区間 [L, R] は _closed
inline long long count_congruent(long long l, long long r, long long a, long long m) {
    if (m <= 0) math_internal::fail("count_congruent", "m > 0 で使う");
    if (l >= r) return 0;
    using math_internal::ceil_div128;
    return (long long)(ceil_div128((Int128)r - a, m) - ceil_div128((Int128)l - a, m));
}
inline long long count_congruent_closed(long long L, long long R, long long a, long long m) {
    if (m <= 0) math_internal::fail("count_congruent_closed", "m > 0 で使う");
    if (L > R) return 0;
    using math_internal::ceil_div128;
    return (long long)(ceil_div128((Int128)R + 1 - a, m) - ceil_div128((Int128)L - a, m));
}

// ---- 累乗根（誤差なし。isqrt / isqrt_ceil は 00_int_common） ----
inline long long kth_root_floor(long long n, int k) {  // x^k ≤ n となる最大の x ≥ 0（n ≥ 0, k ≥ 1）
    if (n < 0 || k < 1) math_internal::fail("kth_root_floor", "n ≥ 0, k ≥ 1 で使う");
    if (k == 1 || n <= 1) return n;
    if (k == 2) return isqrt(n);
    if (k >= 63) return 1;  // 2^63 > n
    long long x = (long long)std::pow((double)n, 1.0 / k);  // 誤差は ±数（下で直す）
    while (x > 0 && !math_internal::pow_le(x, k, n)) --x;
    while (math_internal::pow_le(x + 1, k, n)) ++x;
    return x;
}
inline long long kth_root_ceil(long long n, int k) {  // x^k ≥ n となる最小の x ≥ 0（n ≥ 0, k ≥ 1）
    const long long x = kth_root_floor(n, k);
    return n > 0 && math_internal::pow_le(x, k, n - 1) ? x + 1 : x;  // x^k < n なら x + 1
}
inline bool is_kth_power(long long n, int k) {  // n = x^k となる整数 x ≥ 0 があるか（n ≥ 0）
    if (n < 0) return false;
    const long long x = kth_root_floor(n, k);
    return !math_internal::pow_le(x, k, n - 1);
}
inline bool is_square(long long n) {  // 平方数か（0 も平方数。負は false）
    if (n < 0) return false;
    const long long r = isqrt(n);
    return r * r == n;
}

// ---- 累乗・溢れ ----
inline bool mul_overflows(long long a, long long b) {  // a × b が long long に収まらないか
    long long r;
    return __builtin_mul_overflow(a, b, &r);
}
inline bool add_overflows(long long a, long long b) {  // a + b が long long に収まらないか
    long long r;
    return __builtin_add_overflow(a, b, &r);
}
inline long long mul_cap(long long a, long long b, long long cap) {  // min(a × b, cap)（a, b, cap ≥ 0。溢れない）
    if (a < 0 || b < 0 || cap < 0) math_internal::fail("mul_cap", "a, b, cap ≥ 0 で使う");
    const Int128 r = (Int128)a * b;
    return r > cap ? cap : (long long)r;
}
inline long long add_cap(long long a, long long b, long long cap) {  // min(a + b, cap)（a, b, cap ≥ 0。溢れない）
    if (a < 0 || b < 0 || cap < 0) math_internal::fail("add_cap", "a, b, cap ≥ 0 で使う");
    const Int128 r = (Int128)a + b;
    return r > cap ? cap : (long long)r;
}
inline long long pow_cap(long long a, long long e, long long cap) {  // min(a^e, cap)（a, e, cap ≥ 0。0^0 = 1。溢れない）
    if (a < 0 || e < 0 || cap < 0) math_internal::fail("pow_cap", "a, e, cap ≥ 0 で使う");
    long long r = 1, x = a;
    while (e > 0) {
        if (e & 1) r = mul_cap(r, x, cap);
        e >>= 1;
        if (e) x = mul_cap(x, x, cap);
        if (r == 0 || (r >= cap && x >= 1)) break;
    }
    return std::min(r, cap);
}
inline long long pow_int(long long a, long long e) {  // a^e を厳密に（e ≥ 0。a は負でもよい。long long に収まらなければ止める）
    if (e < 0) math_internal::fail("pow_int", "e ≥ 0 で使う");
    long long r = 1;
    while (e > 0) {
        if ((e & 1) && __builtin_mul_overflow(r, a, &r)) math_internal::fail("pow_int", "long long に収まらない");
        e >>= 1;
        if (e && __builtin_mul_overflow(a, a, &a)) math_internal::fail("pow_int", "long long に収まらない");
    }
    return r;
}
inline int ilog_floor(long long n, long long b) {  // b^k ≤ n となる最大の k（n ≥ 1, b ≥ 2）。10 進の桁数は ilog_floor(n, 10) + 1
    if (n < 1 || b < 2) math_internal::fail("ilog_floor", "n ≥ 1, b ≥ 2 で使う");
    int k = 0;
    for (long long p = 1; p <= n / b; p *= b) ++k;
    return k;
}
inline int ilog_ceil(long long n, long long b) {  // b^k ≥ n となる最小の k（n ≥ 1, b ≥ 2）
    const int k = ilog_floor(n, b);
    return pow_cap(b, k, n) == n ? k : k + 1;
}
inline long long pow10_ll(int k) {  // 10^k（0 ≤ k ≤ 18）
    static const long long t[19] = {1LL, 10LL, 100LL, 1000LL, 10000LL, 100000LL, 1000000LL, 10000000LL, 100000000LL, 1000000000LL,
                                    10000000000LL, 100000000000LL, 1000000000000LL, 10000000000000LL, 100000000000000LL,
                                    1000000000000000LL, 10000000000000000LL, 100000000000000000LL, 1000000000000000000LL};
    if (k < 0 || k > 18) math_internal::fail("pow10_ll", "0 ≤ k ≤ 18 で使う");
    return t[k];
}

// ---- gcd・lcm（gcd / lcm そのものは std::gcd / std::lcm） ----
struct ExtGcd {  // a x + b y = g（g = gcd(|a|, |b|) ≥ 0）
    long long g, x, y;
};
inline ExtGcd ext_gcd(long long a, long long b) {  // |x| ≤ max(1, |b| / g)、|y| ≤ max(1, |a| / g)。a = b = 0 なら {0, 0, 0}
    if (a == 0 && b == 0) return ExtGcd{0, 0, 0};
    auto run = [](auto r0, auto r1) {  // |a|, |b| < 2^62 なら long long のまま（速い）、それ以外は Int128
        decltype(r0) s0 = 1, s1 = 0, t0 = 0, t1 = 1;
        while (r1 != 0) {
            const auto q = r0 / r1;
            auto t = r0 - q * r1; r0 = r1, r1 = t;
            t = s0 - q * s1; s0 = s1, s1 = t;
            t = t0 - q * t1; t0 = t1, t1 = t;
        }
        if (r0 < 0) r0 = -r0, s0 = -s0, t0 = -t0;
        return ExtGcd{(long long)r0, (long long)s0, (long long)t0};
    };
    const long long L = 1LL << 62;
    if (-L < a && a < L && -L < b && b < L) return run(a, b);
    return run((Int128)a, (Int128)b);
}
inline long long lcm_cap(long long a, long long b, long long cap) {  // min(lcm(a, b), cap)（a, b, cap ≥ 0。lcm(0, x) = 0。溢れない）
    if (a < 0 || b < 0 || cap < 0) math_internal::fail("lcm_cap", "a, b, cap ≥ 0 で使う");
    if (a == 0 || b == 0) return 0;
    return mul_cap(a / (long long)math_internal::ugcd((unsigned long long)a, (unsigned long long)b), b, cap);
}
inline Maybe<long long> try_lcm(long long a, long long b) {  // lcm(|a|, |b|)。long long に収まらなければ ok = false（lcm(0, x) = 0）
    const unsigned long long x = math_internal::uabs(a), y = math_internal::uabs(b);
    if (x == 0 || y == 0) return {true, 0};
    const UInt128 l = (UInt128)(x / math_internal::ugcd(x, y)) * y;
    if (l > (UInt128)9223372036854775807ULL) return {false, 0};
    return {true, (long long)l};
}
template <class V>
long long gcd_of(const V& v) {  // 全部の gcd（絶対値で。空なら 0。vector<int> などどの整数の列でも）
    unsigned long long g = 0;
    for (const auto& x : v) g = math_internal::ugcd(g, math_internal::uabs(math_internal::to_ll_checked(x, "gcd_of")));
    if (g > 9223372036854775807ULL) math_internal::fail("gcd_of", "gcd = 2^63 は long long に収まらない");
    return (long long)g;
}
inline long long gcd_of(std::initializer_list<long long> v) { return gcd_of<std::initializer_list<long long>>(v); }  // gcd_of({a, b, c})
template <class V>
long long lcm_of_cap(const V& v, long long cap) {  // min(全部の lcm, cap)（絶対値で。空なら 1。0 があれば 0。cap ≥ 0）
    if (cap < 0) math_internal::fail("lcm_of_cap", "cap ≥ 0 で使う");
    long long r = 1;
    for (const auto& xx : v) {
        const unsigned long long u = math_internal::uabs(math_internal::to_ll_checked(xx, "lcm_of_cap"));
        if (u == 0) return 0;
        if (r < cap) r = u > (unsigned long long)cap ? cap : lcm_cap(r, (long long)u, cap);  // u > cap なら lcm も cap を超える
    }
    return std::min(r, cap);
}
inline long long lcm_of_cap(std::initializer_list<long long> v, long long cap) { return lcm_of_cap<std::initializer_list<long long>>(v, cap); }

// ---- 64 ビットの mod（m は 1 以上 2^63 未満。a は負でもよい。答えは [0, m)） ----
inline long long mod_mul(long long a, long long b, long long m) {
    if (m < 1) math_internal::fail("mod_mul", "m ≥ 1 で使う");
    Int128 r = (Int128)a * b % m;
    if (r < 0) r += m;
    return (long long)r;
}
inline Maybe<long long> try_mod_inv(long long a, long long m) {  // a x ≡ 1 (mod m) の x。gcd(a, m) ≠ 1 なら ok = false
    if (m < 1) math_internal::fail("try_mod_inv", "m ≥ 1 で使う");
    const ExtGcd e = ext_gcd(mod_floor(a, m), m);
    if (e.g != 1) return {false, 0};
    return {true, mod_floor(e.x, m)};
}
inline long long mod_inv(long long a, long long m) {  // a の逆元（m は合成数でもよい。gcd(a, m) ≠ 1 なら止める）
    const Maybe<long long> r = try_mod_inv(a, m);
    if (!r.ok) math_internal::fail("mod_inv", "逆元がない（gcd(a, m) ≠ 1）");
    return r.value;
}
inline long long mod_pow(long long a, long long e, long long m) {  // a^e mod m（e < 0 なら逆元の (-e) 乗。0^0 = 1）
    if (m < 1) math_internal::fail("mod_pow", "m ≥ 1 で使う");
    if (e < 0) a = mod_inv(a, m);
    unsigned long long ue = math_internal::uabs(e), x = (unsigned long long)mod_floor(a, m), r = 1 % (unsigned long long)m;
    const unsigned long long um = (unsigned long long)m;
    if (um <= (1ULL << 32)) {
        for (; ue; ue >>= 1, x = x * x % um)
            if (ue & 1) r = r * x % um;
    } else {
        for (; ue; ue >>= 1, x = (unsigned long long)((UInt128)x * x % um))
            if (ue & 1) r = (unsigned long long)((UInt128)r * x % um);
    }
    return (long long)r;
}
inline long long geom_sum_mod(long long r, long long n, long long m) {  // (1 + r + r^2 + … + r^(n-1)) mod m（n ≥ 0。逆元を使わないので m は何でも）
    if (m < 1 || n < 0) math_internal::fail("geom_sum_mod", "m ≥ 1, n ≥ 0 で使う");
    const unsigned long long um = (unsigned long long)m, ur = (unsigned long long)mod_floor(r, m);
    auto mul = [&](unsigned long long x, unsigned long long y) {
        return um <= (1ULL << 32) ? x * y % um : (unsigned long long)((UInt128)x * y % um);
    };
    unsigned long long s = 0, p = 1 % um;  // 長さ L の和 s と r^L
    for (int bit = 62; bit >= 0; --bit) {
        s = mul(s, (p + 1) % um), p = mul(p, p);  // L → 2L
        if ((n >> bit) & 1) s = (mul(s, ur) + 1) % um, p = mul(p, ur);  // L → L + 1
    }
    return (long long)s;
}

// ---- 数列の和 ----
inline Int128 arith_sum(long long a, long long d, long long n) {  // a + (a + d) + … + (a + (n-1)d)（n ≥ 0。Int128 に収まらなければ止める）
    if (n < 0) math_internal::fail("arith_sum", "n ≥ 0 で使う");
    const Int128 t = (Int128)n * (n - 1) / 2;
    Int128 x, y, z;
    if (__builtin_mul_overflow((Int128)n, (Int128)a, &x) || __builtin_mul_overflow(t, (Int128)d, &y) || __builtin_add_overflow(x, y, &z))
        math_internal::fail("arith_sum", "Int128 に収まらない（arith_sum_mod を使う）");
    return z;
}
inline long long arith_sum_mod(long long a, long long d, long long n, long long m) {  // arith_sum を mod m で（m は何でも）
    if (n < 0 || m < 1) math_internal::fail("arith_sum_mod", "n ≥ 0, m ≥ 1 で使う");
    const unsigned long long t = (unsigned long long)((UInt128)n * (UInt128)(n - (n > 0)) / 2 % (UInt128)m);
    return (long long)(((Int128)mod_floor(n, m) * mod_floor(a, m) + (Int128)mod_floor(d, m) * (Int128)t) % m);
}
inline Int128 sum_range(long long l, long long r) {  // l + (l+1) + … + (r-1)（l ≥ r なら 0）
    if (l >= r) return 0;
    const Int128 cnt = (Int128)r - l, s = (Int128)l + r - 1;
    Int128 res;
    if (__builtin_mul_overflow(cnt, s, &res)) math_internal::fail("sum_range", "Int128 に収まらない");
    return res / 2;
}
inline Int128 sum_range_closed(long long L, long long R) {  // L + (L+1) + … + R（L > R なら 0）
    if (L > R) return 0;
    const Int128 cnt = (Int128)R - L + 1, s = (Int128)L + R;
    Int128 res;
    if (__builtin_mul_overflow(cnt, s, &res)) math_internal::fail("sum_range_closed", "Int128 に収まらない");
    return res / 2;
}

inline long long sum_range_mod(long long l, long long r, long long m) {  // (l + (l+1) + … + (r-1)) mod m（[0, m)。l ≥ r なら 0）
    if (m < 1) math_internal::fail("sum_range_mod", "m ≥ 1 で使う");
    if (l >= r) return 0;
    Int128 cnt = (Int128)r - l, sm = (Int128)l + r - 1;  // 和 = cnt·sm/2（どちらかは偶数）
    if (cnt % 2 == 0) cnt /= 2;
    else sm /= 2;
    const Int128 a = ((cnt % m) + m) % m, b = ((sm % m) + m) % m;
    return (long long)(a * b % m);
}
inline long long sum_range_closed_mod(long long L, long long R, long long m) {  // (L + … + R) mod m（L > R なら 0）
    if (m < 1) math_internal::fail("sum_range_closed_mod", "m ≥ 1 で使う");
    if (L > R) return 0;
    if (R < 9223372036854775807LL) return sum_range_mod(L, R + 1, m);
    return (long long)(((Int128)sum_range_mod(L, R, m) + R % m) % m);
}

// ---- 分数の比較（誤差なし） ----
inline int frac_compare(long long a, long long b, long long c, long long d) {  // a/b と c/d の大小（-1 / 0 / 1。b, d ≠ 0）
    Int128 p = a, q = b, s = c, t = d;
    if (q < 0) p = -p, q = -q;
    if (t < 0) s = -s, t = -t;
    const Int128 l = p * t, r = s * q;
    return (l > r) - (l < r);
}
#endif  // GG_MATH_BASE
#ifndef GG_MODINT
#define GG_MODINT
// ---- modint（数学第二世代）  ModInt<MOD>（定数の mod）/ DynModInt（入力で決まる mod）。using mint = ModInt998; を使う側で書く ----
namespace math_internal {
constexpr unsigned mont_ninv(unsigned m) {  // -m^{-1} mod 2^32（m は奇数）
    unsigned x = m;
    for (int i = 0; i < 5; ++i) x *= 2u - m * x;
    return 0u - x;
}
constexpr bool is_prime_u32(unsigned n) {  // 決定的ミラー・ラビン（n < 2^32）
    if (n < 2) return false;
    if (n % 2 == 0) return n == 2;
    unsigned long long d = n - 1;
    int s = 0;
    while (d % 2 == 0) d /= 2, ++s;
    const unsigned bases[3] = {2, 7, 61};
    for (unsigned a : bases) {
        if (a % n == 0) continue;
        unsigned long long x = 1, b = a % n, e = d;
        for (; e; e >>= 1, b = b * b % n)
            if (e & 1) x = x * b % n;
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int i = 1; i < s && comp; ++i) {
            x = x * x % n;
            if (x == n - 1) comp = false;
        }
        if (comp) return false;
    }
    return true;
}
inline long long inv_by_gcd(long long a, long long m) {  // a の逆元 mod m（なければ -1）
    long long r0 = m, r1 = a, t0 = 0, t1 = 1;
    while (r1 != 0) {
        const long long q = r0 / r1;
        long long t = r0 - q * r1; r0 = r1, r1 = t;
        t = t0 - q * t1; t0 = t1, t1 = t;
    }
    if (r0 != 1) return -1;
    return t0 < 0 ? t0 + m : t0;
}
inline unsigned long long splitmix64(unsigned long long z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
template <class M>
M read_modint(std::istream& is) {  // 10 進の文字列を mod で読む（long long を超える長さでもよい）
    std::string s;
    M r = 0;
    if (!(is >> s)) return r;
    size_t i = (s[0] == '-' || s[0] == '+') ? 1 : 0;
    if (i == s.size()) fail("modint の入力", "数字がない");
    for (; i < s.size(); ++i) {
        if (s[i] < '0' || s[i] > '9') fail("modint の入力", "数として読めない");
        r = r * 10 + (s[i] - '0');
    }
    return s[0] == '-' ? -r : r;
}
}  // namespace math_internal

// 定数の mod（1 ≤ MOD < 2^31。素数でなくてよい、偶数でもよい）
// 奇数の MOD はモンゴメリ乗算（MOD < 2^30 なら遅延簡約でさらに速い）、偶数の MOD は普通の剰余
template <unsigned MOD>
struct ModInt {
    static_assert(1 <= MOD && MOD < (1u << 31), "ModInt: MOD は 1 以上 2^31 未満");

  private:
    using u32 = unsigned;
    using u64 = unsigned long long;
    static constexpr bool MONT = (MOD & 1) && MOD > 1;
    static constexpr bool LAZY = MONT && MOD < (1u << 30);  // 内部の値を [0, 2 MOD) のまま持つ
    static constexpr u32 VMAX = LAZY ? 2 * MOD : MOD;
    static constexpr u32 NINV = MONT ? math_internal::mont_ninv(MOD) : 0;
    static constexpr u32 R2 = MONT ? (u32)((~0ULL % MOD + 1) % MOD) : 0;  // 2^64 mod MOD
    static constexpr u32 ONE = MONT ? (u32)((1ULL << 32) % MOD) : 1 % MOD;  // 1 の内部表現
    u32 v_;  // MONT なら x·2^32 mod MOD（LAZY なら [0, 2 MOD)）、それ以外は x

    static u32 reduce(u64 t) { return (u32)((t + (u64)((u32)t * NINV) * MOD) >> 32); }  // t < MOD·2^32 → [0, 2 MOD)
    static u32 mul_in(u32 a, u32 b) {
        if constexpr (MONT) {
            const u32 r = reduce((u64)a * b);
            if constexpr (LAZY) return r;
            else return r >= MOD ? r - MOD : r;
        } else {
            return (u32)((u64)a * b % MOD);
        }
    }
    static u32 to_in(u32 x) {  // 0 ≤ x < MOD
        if constexpr (MONT) return mul_in(x, R2);
        else return x;
    }
    u32 norm_in() const { return v_ >= MOD ? v_ - MOD : v_; }

  public:
    static constexpr bool IS_PRIME = math_internal::is_prime_u32(MOD);
    static constexpr int mod() { return (int)MOD; }

    ModInt() : v_(0) {}
    template <class T, std::enable_if_t<std::is_integral_v<T> && std::is_signed_v<T>, int> = 0>
    ModInt(T x) {
        long long y = (long long)x % (long long)MOD;
        if (y < 0) y += MOD;
        v_ = to_in((u32)y);
    }
    template <class T, std::enable_if_t<std::is_integral_v<T> && std::is_unsigned_v<T>, int> = 0>
    ModInt(T x) : v_(to_in((u32)((u64)x % MOD))) {}
    ModInt(Int128 x) {
        Int128 y = x % (Int128)MOD;
        if (y < 0) y += MOD;
        v_ = to_in((u32)y);
    }
    ModInt(UInt128 x) : v_(to_in((u32)(x % MOD))) {}
    static ModInt raw(int x) {  // 0 ≤ x < MOD が分かっているとき（剰余を取らないので速い）
        ModInt r;
        r.v_ = to_in((u32)x);
        return r;
    }

    int val() const {
        if constexpr (MONT) {
            const u32 r = reduce(v_);
            return (int)(r >= MOD ? r - MOD : r);
        } else {
            return (int)v_;
        }
    }
    explicit operator bool() const { return norm_in() != 0; }

    ModInt& operator+=(const ModInt& o) {
        v_ += o.v_;
        if (v_ >= VMAX) v_ -= VMAX;
        return *this;
    }
    ModInt& operator-=(const ModInt& o) {
        if (v_ < o.v_) v_ += VMAX;
        v_ -= o.v_;
        return *this;
    }
    ModInt& operator*=(const ModInt& o) {
        v_ = mul_in(v_, o.v_);
        return *this;
    }
    ModInt& operator/=(const ModInt& o) { return *this *= o.inv(); }
    ModInt operator+() const { return *this; }
    ModInt operator-() const {
        ModInt r;
        r.v_ = v_ ? VMAX - v_ : 0;
        return r;
    }
    ModInt& operator++() {
        v_ += ONE;
        if (v_ >= VMAX) v_ -= VMAX;
        return *this;
    }
    ModInt& operator--() {
        if (v_ < ONE) v_ += VMAX;
        v_ -= ONE;
        return *this;
    }
    ModInt operator++(int) { ModInt t = *this; ++*this; return t; }
    ModInt operator--(int) { ModInt t = *this; --*this; return t; }
    friend ModInt operator+(ModInt a, const ModInt& b) { return a += b; }
    friend ModInt operator-(ModInt a, const ModInt& b) { return a -= b; }
    friend ModInt operator*(ModInt a, const ModInt& b) { return a *= b; }
    friend ModInt operator/(ModInt a, const ModInt& b) { return a /= b; }
    friend bool operator==(const ModInt& a, const ModInt& b) { return a.norm_in() == b.norm_in(); }
    friend bool operator!=(const ModInt& a, const ModInt& b) { return a.norm_in() != b.norm_in(); }
    friend bool operator<(const ModInt& a, const ModInt& b) { return a.val() < b.val(); }  // 値（0〜MOD-1）での比較。map のキー用
    friend bool operator>(const ModInt& a, const ModInt& b) { return a.val() > b.val(); }
    friend bool operator<=(const ModInt& a, const ModInt& b) { return a.val() <= b.val(); }
    friend bool operator>=(const ModInt& a, const ModInt& b) { return a.val() >= b.val(); }

    ModInt pow(long long e) const {  // e < 0 なら逆元の (-e) 乗
        ModInt x = e < 0 ? inv() : *this, r;
        r.v_ = ONE;
        for (unsigned long long u = math_internal::uabs(e); u; u >>= 1, x *= x)
            if (u & 1) r *= x;
        return r;
    }
    bool is_unit() const { return std::gcd((unsigned)val(), MOD) == 1; }  // 逆元があるか
    Maybe<ModInt> try_inv() const {
        if constexpr (IS_PRIME) {
            if (norm_in() == 0) return {false, ModInt()};
            return {true, pow(MOD - 2)};
        } else {
            const long long r = math_internal::inv_by_gcd(val(), MOD);
            if (r < 0) return {false, ModInt()};
            return {true, ModInt::raw((int)r)};
        }
    }
    ModInt inv() const {  // 逆元（MOD が合成数でも gcd = 1 なら。なければ止める）
        const Maybe<ModInt> r = try_inv();
        if (!r.ok) math_internal::fail("ModInt::inv", "逆元がない（0 か、mod と互いに素でない）");
        return r.value;
    }
    friend std::ostream& operator<<(std::ostream& os, const ModInt& a) { return os << a.val(); }
    friend std::istream& operator>>(std::istream& is, ModInt& a) {
        a = math_internal::read_modint<ModInt>(is);
        return is;
    }
};
using ModInt998 = ModInt<998244353>;
using ModInt107 = ModInt<1000000007>;

// 入力で決まる mod（DynModInt::set_mod(M)。1 ≤ M < 2^31。素数でなくてよい、偶数でもよい。最初は 998244353）
// 奇数で 2^30 未満の M はモンゴメリ乗算（遅延簡約）、それ以外は Barrett 簡約。ID を変えると別の mod を同時に持てる
template <int ID>
struct DynModIntT {
  private:
    using u32 = unsigned;
    using u64 = unsigned long long;
    struct State {  // 最初は 998244353
        u32 m = 998244353, vmax = 2 * 998244353u, ninv = math_internal::mont_ninv(998244353);
        u32 r2 = (u32)((~0ULL % 998244353 + 1) % 998244353), one = (u32)((1ULL << 32) % 998244353);
        u64 im = ~0ULL / 998244353 + 1;
        bool mont = true, prime = true;
    };
    inline static State s_{};
    static const State& st() { return s_; }
    u32 v_;

    static u32 mont_reduce(u64 t, u32 m, u32 ninv) { return (u32)((t + (u64)((u32)t * ninv) * m) >> 32); }
    static u32 barrett(u64 z, u32 m, u64 im) {  // z < m^2
        const u64 x = (u64)(((UInt128)z * im) >> 64);
        const u64 y = x * m;
        return (u32)(z - y + (z < y ? m : 0));
    }
    static u32 mul_in(u32 a, u32 b) {
        const State& s = st();
        if (s.mont) return mont_reduce((u64)a * b, s.m, s.ninv);
        return barrett((u64)a * b, s.m, s.im);
    }
    static u32 to_in(u32 x) {  // 0 ≤ x < m
        const State& s = st();
        return s.mont ? mont_reduce((u64)x * s.r2, s.m, s.ninv) : x;
    }
    u32 norm_in() const { return v_ >= st().m ? v_ - st().m : v_; }

  public:
    static void set_mod(long long m) {
        if (m < 1 || m >= (1LL << 31)) math_internal::fail("DynModInt::set_mod", "mod は 1 以上 2^31 未満");
        State& s = s_;
        s.m = (u32)m;
        s.mont = (m & 1) && m > 1 && m < (1LL << 30);
        s.vmax = s.mont ? 2 * s.m : s.m;
        s.im = ~0ULL / s.m + 1;
        s.prime = math_internal::is_prime_u32(s.m);
        if (s.mont) s.ninv = math_internal::mont_ninv(s.m), s.r2 = (u32)((~0ULL % s.m + 1) % s.m);
        s.one = s.mont ? (u32)((1ULL << 32) % s.m) : 1 % s.m;
    }
    static int mod() { return (int)st().m; }
    static bool mod_is_prime() { return st().prime; }

    DynModIntT() : v_(0) {}
    template <class T, std::enable_if_t<std::is_integral_v<T> && std::is_signed_v<T>, int> = 0>
    DynModIntT(T x) {
        const long long m = st().m;
        long long y = (long long)x % m;
        if (y < 0) y += m;
        v_ = to_in((u32)y);
    }
    template <class T, std::enable_if_t<std::is_integral_v<T> && std::is_unsigned_v<T>, int> = 0>
    DynModIntT(T x) : v_(to_in((u32)((u64)x % st().m))) {}
    DynModIntT(Int128 x) {
        Int128 y = x % (Int128)st().m;
        if (y < 0) y += st().m;
        v_ = to_in((u32)y);
    }
    DynModIntT(UInt128 x) : v_(to_in((u32)(x % st().m))) {}
    static DynModIntT raw(int x) {  // 0 ≤ x < mod が分かっているとき
        DynModIntT r;
        r.v_ = to_in((u32)x);
        return r;
    }

    int val() const {
        const State& s = st();
        if (!s.mont) return (int)v_;
        const u32 r = mont_reduce(v_, s.m, s.ninv);
        return (int)(r >= s.m ? r - s.m : r);
    }
    explicit operator bool() const { return norm_in() != 0; }

    DynModIntT& operator+=(const DynModIntT& o) {
        const u32 vmax = st().vmax;
        v_ += o.v_;
        if (v_ >= vmax) v_ -= vmax;
        return *this;
    }
    DynModIntT& operator-=(const DynModIntT& o) {
        if (v_ < o.v_) v_ += st().vmax;
        v_ -= o.v_;
        return *this;
    }
    DynModIntT& operator*=(const DynModIntT& o) {
        v_ = mul_in(v_, o.v_);
        return *this;
    }
    DynModIntT& operator/=(const DynModIntT& o) { return *this *= o.inv(); }
    DynModIntT operator+() const { return *this; }
    DynModIntT operator-() const {
        DynModIntT r;
        r.v_ = v_ ? st().vmax - v_ : 0;
        return r;
    }
    DynModIntT& operator++() {
        v_ += st().one;
        if (v_ >= st().vmax) v_ -= st().vmax;
        return *this;
    }
    DynModIntT& operator--() {
        if (v_ < st().one) v_ += st().vmax;
        v_ -= st().one;
        return *this;
    }
    DynModIntT operator++(int) { DynModIntT t = *this; ++*this; return t; }
    DynModIntT operator--(int) { DynModIntT t = *this; --*this; return t; }
    friend DynModIntT operator+(DynModIntT a, const DynModIntT& b) { return a += b; }
    friend DynModIntT operator-(DynModIntT a, const DynModIntT& b) { return a -= b; }
    friend DynModIntT operator*(DynModIntT a, const DynModIntT& b) { return a *= b; }
    friend DynModIntT operator/(DynModIntT a, const DynModIntT& b) { return a /= b; }
    friend bool operator==(const DynModIntT& a, const DynModIntT& b) { return a.norm_in() == b.norm_in(); }
    friend bool operator!=(const DynModIntT& a, const DynModIntT& b) { return a.norm_in() != b.norm_in(); }
    friend bool operator<(const DynModIntT& a, const DynModIntT& b) { return a.val() < b.val(); }
    friend bool operator>(const DynModIntT& a, const DynModIntT& b) { return a.val() > b.val(); }
    friend bool operator<=(const DynModIntT& a, const DynModIntT& b) { return a.val() <= b.val(); }
    friend bool operator>=(const DynModIntT& a, const DynModIntT& b) { return a.val() >= b.val(); }

    DynModIntT pow(long long e) const {
        DynModIntT x = e < 0 ? inv() : *this, r;
        r.v_ = st().one;
        for (unsigned long long u = math_internal::uabs(e); u; u >>= 1, x *= x)
            if (u & 1) r *= x;
        return r;
    }
    bool is_unit() const { return std::gcd((unsigned)val(), st().m) == 1; }
    Maybe<DynModIntT> try_inv() const {
        if (st().prime) {
            if (norm_in() == 0) return {false, DynModIntT()};
            return {true, pow(st().m - 2)};
        }
        const long long r = math_internal::inv_by_gcd(val(), st().m);
        if (r < 0) return {false, DynModIntT()};
        return {true, DynModIntT::raw((int)r)};
    }
    DynModIntT inv() const {
        const Maybe<DynModIntT> r = try_inv();
        if (!r.ok) math_internal::fail("DynModInt::inv", "逆元がない（0 か、mod と互いに素でない）");
        return r.value;
    }
    friend std::ostream& operator<<(std::ostream& os, const DynModIntT& a) { return os << a.val(); }
    friend std::istream& operator>>(std::istream& is, DynModIntT& a) {
        a = math_internal::read_modint<DynModIntT>(is);
        return is;
    }
};
using DynModInt = DynModIntT<0>;

// unordered_map<mint, V> / unordered_set<mint> がそのまま使えるようにする
namespace std {
template <unsigned MOD>
struct hash<ModInt<MOD>> {
    size_t operator()(const ModInt<MOD>& a) const { return (size_t)math_internal::splitmix64((unsigned long long)a.val()); }
};
template <int ID>
struct hash<DynModIntT<ID>> {
    size_t operator()(const DynModIntT<ID>& a) const { return (size_t)math_internal::splitmix64((unsigned long long)a.val()); }
};
}  // namespace std

// mint の値を小さい分数 num/den に戻す（デバッグ用：期待値などの答えが合っているか確かめる）
// |num|, den ≤ √(mod/2) の分数があれば見つける（なければ ok = false）。mod は素数でなくてもよい
struct RationalGuess {
    bool ok = false;
    long long num = 0, den = 1;
    friend std::ostream& operator<<(std::ostream& os, const RationalGuess& r) {
        if (!r.ok) return os << "(分数が見つからない)";
        return r.den == 1 ? os << r.num : os << r.num << "/" << r.den;
    }
};
template <class M>
RationalGuess rational_reconstruct(const M& x) {
    const long long m = M::mod(), bound = isqrt(m / 2);
    long long r0 = m, r1 = x.val(), t0 = 0, t1 = 1;
    while (r1 > bound) {
        const long long q = r0 / r1;
        long long t = r0 - q * r1; r0 = r1, r1 = t;
        t = t0 - q * t1; t0 = t1, t1 = t;
    }
    if (t1 < 0) r1 = -r1, t1 = -t1;
    if (t1 == 0 || t1 > bound || std::gcd(std::llabs(r1), t1) != 1) return RationalGuess{};
    return RationalGuess{true, r1, t1};
}
#endif  // GG_MODINT
#ifndef GG_SERIES
#define GG_SERIES
// ---- べき乗和・ラグランジュ補間・次数の小さい多項式（数学第二世代） ----
namespace math_internal {
template <class mint>
void check_prime_mod_over(long long need, const char* func) {  // mod が素数で need より大きいか（階乗の逆元を使うため）
    if (!is_prime_u32((unsigned)mint::mod())) fail(func, "mod が素数でない");
    if (need >= mint::mod()) fail(func, "点の数が mod 以上");
}
template <class mint>
std::vector<mint> inv_factorials(int n) {  // 1/0!〜1/n!
    std::vector<mint> f((size_t)n + 1, mint(1));
    for (int i = 1; i <= n; ++i) f[(size_t)i] = f[(size_t)i - 1] * mint(i);
    f[(size_t)n] = f[(size_t)n].inv();
    for (int i = n; i > 0; --i) f[(size_t)i - 1] = f[(size_t)i] * mint(i);
    return f;
}
}  // namespace math_internal

// ---- ラグランジュ補間 ----
// f(0), f(1), …, f(d) が ys のとき、d 次以下の多項式 f の f(x)（x は 10^18 でも負でもよい。mod は d より大きい素数）  O(d)
template <class mint>
mint lagrange_interpolate(const std::vector<mint>& ys, long long x) {
    const int n = (int)ys.size();  // 点の数 = d + 1
    if (n == 0) return mint(0);
    math_internal::check_prime_mod_over<mint>(n - 1, "lagrange_interpolate");
    if (0 <= x && x < n) return ys[(size_t)x];
    const mint X(x);
    std::vector<mint> pre((size_t)n + 1, mint(1)), suf((size_t)n + 1, mint(1));  // Π (x − j) の左右からの積
    for (int i = 0; i < n; ++i) pre[(size_t)i + 1] = pre[(size_t)i] * (X - mint(i));
    for (int i = n - 1; i >= 0; --i) suf[(size_t)i] = suf[(size_t)i + 1] * (X - mint(i));
    const std::vector<mint> fi = math_internal::inv_factorials<mint>(n - 1);
    mint r = 0;
    for (int i = 0; i < n; ++i) {
        mint t = ys[(size_t)i] * pre[(size_t)i] * suf[(size_t)i + 1] * fi[(size_t)i] * fi[(size_t)(n - 1 - i)];
        r += ((n - 1 - i) % 2 == 0) ? t : -t;
    }
    return r;
}
// 点 (xs[i], ys[i])（xs は互いに違う）を通る多項式の x での値（mod は素数）  O(d^2)
template <class mint>
mint lagrange_interpolate_at(const std::vector<mint>& xs, const std::vector<mint>& ys, mint x) {
    if (xs.size() != ys.size()) math_internal::fail("lagrange_interpolate_at", "xs と ys の長さが違う");
    const size_t n = xs.size();
    mint r = 0;
    for (size_t i = 0; i < n; ++i) {
        mint num = 1, den = 1;
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            if (xs[i] == xs[j]) math_internal::fail("lagrange_interpolate_at", "xs に同じ値がある");
            num *= x - xs[j], den *= xs[i] - xs[j];
        }
        r += ys[i] * num / den;
    }
    return r;
}
// 点 (xs[i], ys[i]) を通る多項式の係数（低い次数から。長さ n）  O(n^2)
template <class mint>
std::vector<mint> poly_from_points(const std::vector<mint>& xs, const std::vector<mint>& ys) {
    if (xs.size() != ys.size()) math_internal::fail("poly_from_points", "xs と ys の長さが違う");
    const size_t n = xs.size();
    std::vector<mint> all(n + 1, mint(0)), res(n, mint(0));  // all = Π (x − xs[j])
    all[0] = 1;
    for (size_t j = 0; j < n; ++j) {  // all ← all × (x − xs[j])（上の次数から更新）
        for (size_t k = j + 1; k > 0; --k) all[k] = all[k - 1] - xs[j] * all[k];
        all[0] = -xs[j] * all[0];
    }
    for (size_t i = 0; i < n; ++i) {
        mint den = 1;
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            if (xs[i] == xs[j]) math_internal::fail("poly_from_points", "xs に同じ値がある");
            den *= xs[i] - xs[j];
        }
        const mint c = ys[i] / den;
        // q = all / (x − xs[i])（組立除法）を c 倍して足す
        mint carry = 0;
        for (size_t k = n; k > 0; --k) {
            carry = all[k] + carry * xs[i];  // q の k−1 次の係数
            res[k - 1] += c * carry;
        }
    }
    return res;
}

// ---- べき乗和 ----
// Σ_{i=0}^{n−1} i^k（0^0 = 1。n ≤ 10^18、k ≤ 10^7 程度。mod は k + 1 より大きい素数）  O(k)
template <class mint>
mint sum_powers(long long n, int k) {
    if (n < 0 || k < 0) math_internal::fail("sum_powers", "n, k ≥ 0 で使う");
    // f(x) = Σ_{i<x} i^k は k + 1 次の多項式。f(0)〜f(k+1) を作って補間
    const int m = k + 2;
    if ((long long)m - 1 >= mint::mod()) math_internal::fail("sum_powers", "k + 1 が mod 以上");
    std::vector<mint> pw((size_t)m, mint(0));  // pw[i] = i^k（線形篩で：素数だけ累乗）
    std::vector<int> spf((size_t)m, 0), primes;
    if (m > 0) pw[0] = (k == 0 ? mint(1) : mint(0));
    if (m > 1) pw[1] = 1;
    for (int i = 2; i < m; ++i) {
        if (spf[(size_t)i] == 0) spf[(size_t)i] = i, primes.push_back(i), pw[(size_t)i] = mint(i).pow(k);
        for (int p : primes) {
            if (p > spf[(size_t)i] || (long long)p * i >= m) break;
            spf[(size_t)(p * i)] = p, pw[(size_t)(p * i)] = pw[(size_t)p] * pw[(size_t)i];
        }
    }
    std::vector<mint> ys((size_t)m, mint(0));
    for (int i = 1; i < m; ++i) ys[(size_t)i] = ys[(size_t)i - 1] + pw[(size_t)i - 1];
    return lagrange_interpolate(ys, n);
}
template <class mint>
mint sum_powers(long long l, long long r, int k) {  // Σ_{i=l}^{r−1} i^k（0 ≤ l。空なら 0）
    if (l < 0) math_internal::fail("sum_powers", "l ≥ 0 で使う");
    if (r <= l) return mint(0);
    return sum_powers<mint>(r, k) - sum_powers<mint>(l, k);
}
template <class mint>
mint sum_powers_closed(long long L, long long R, int k) {  // Σ_{i=L}^{R} i^k（0 ≤ L ≤ R < 2^63 − 1）
    if (L < 0) math_internal::fail("sum_powers_closed", "L ≥ 0 で使う");
    if (R < L) return mint(0);
    return sum_powers<mint>(L, R + 1, k);
}
// Σ_{i=0}^{n−1} i·r^i（mod は何でもよい）  O(log n)
template <class mint>
mint arith_geom_sum(mint r, long long n) {
    if (n < 0) math_internal::fail("arith_geom_sum", "n ≥ 0 で使う");
    // (len, S0 = Σ r^i, S1 = Σ i r^i, P = r^len) を上の桁から 2 倍・+1 で作る
    mint s0 = 0, s1 = 0, p = 1;
    long long len = 0;
    for (int b = 62; b >= 0; --b) {
        // 2 倍：[0, 2len) = [0, len) + len + [0, len)
        s1 = s1 + p * (s1 + mint(len) * s0), s0 = s0 * (mint(1) + p), p = p * p, len *= 2;
        if ((n >> b) & 1) s1 += mint(len) * p, s0 += p, p *= r, ++len;  // 1 項足す（i = len の項 len·r^len）
    }
    return s1;
}

// ---- 次数の小さい多項式（係数は低い次数から。T は mint・long long・double など） ----
template <class T>
T poly_eval(const std::vector<T>& c, const T& x) {  // c(x)（ホーナー法）  O(n)
    T r = T(0);
    for (size_t i = c.size(); i > 0; --i) r = r * x + c[i - 1];
    return r;
}
template <class T>
std::vector<T> poly_mul_naive(const std::vector<T>& a, const std::vector<T>& b) {  // a × b（長さ |a| + |b| − 1）  O(|a||b|)
    if (a.empty() || b.empty()) return {};
    std::vector<T> c(a.size() + b.size() - 1, T(0));
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j) c[i + j] = c[i + j] + a[i] * b[j];
    return c;
}
template <class T>
struct PolyDivMod {  // a = q × b + r（deg r < deg b）
    std::vector<T> q, r;
};
// a ÷ b（b の最高次の係数で割る。T が整数なら割り切れないと止める。商・余りの末尾の 0 は除く）  O(|a||b|)
template <class T>
PolyDivMod<T> poly_divmod(std::vector<T> a, std::vector<T> b) {
    auto trim = [](std::vector<T>& v) {
        while (!v.empty() && v.back() == T(0)) v.pop_back();
    };
    trim(a), trim(b);
    if (b.empty()) math_internal::fail("poly_divmod", "0 で割った");
    PolyDivMod<T> res;
    if (a.size() < b.size()) {
        res.r = a;
        return res;
    }
    const size_t n = a.size(), m = b.size();
    res.q.assign(n - m + 1, T(0));
    const T lc = b.back();
    for (size_t i = n - m + 1; i > 0; --i) {
        const size_t top = i - 1 + m - 1;  // 今の最高次
        T c;
        if constexpr (std::is_integral_v<T>) {
            if (a[top] % lc != 0) math_internal::fail("poly_divmod", "整数の多項式で割り切れない（最高次の係数で割れない）");
            c = a[top] / lc;
        } else {
            c = a[top] / lc;
        }
        res.q[i - 1] = c;
        for (size_t j = 0; j < m; ++j) a[i - 1 + j] = a[i - 1 + j] - c * b[j];
    }
    a.resize(m - 1);
    trim(a), trim(res.q);
    res.r = a;
    return res;
}
#endif  // GG_SERIES
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace ser_test {

#define SER_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[series] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

template <class mint>
static void test_sum_powers(std::mt19937_64& rng, int& cases, int kmax) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int k = (int)rnd(rng, 0, kmax);
        const long long n = rnd(rng, 0, 150);
        mint e = 0;
        for (long long i = 0; i < n; ++i) e += mint(i).pow(k);  // 0^0 = 1
        SER_CHECK(sum_powers<mint>(n, k) == e);
        const long long l = rnd(rng, 0, 150), r = rnd(rng, 0, 150);
        mint e2 = 0;
        for (long long i = l; i < r; ++i) e2 += mint(i).pow(k);
        SER_CHECK(sum_powers<mint>(l, r, k) == e2 && sum_powers_closed<mint>(l, r - 1, k) == e2);
    }
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {  // 大きい n：1 つ前との差、k = 1, 2 の公式
        const long long n = rnd(rng, 1, 1000000000000000000LL);
        const int k = (int)rnd(rng, 1, kmax);
        SER_CHECK(sum_powers<mint>(n + 1, k) - sum_powers<mint>(n, k) == mint(n).pow(k));
        const mint N(n);
        SER_CHECK(sum_powers<mint>(n, 1) * mint(2) == N * (N - mint(1)));
        SER_CHECK(sum_powers<mint>(n, 2) * mint(6) == (N - mint(1)) * N * (mint(2) * N - mint(1)));
    }
}

template <class mint>
static void test_arith_geom(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const mint r = it % 5 == 0 ? mint(rnd(rng, -2, 2)) : mint((long long)rng());
        const long long n = rnd(rng, 0, 60);
        mint e = 0, p = 1;
        for (long long i = 0; i < n; ++i) e += mint(i) * p, p *= r;
        SER_CHECK(arith_geom_sum(r, n) == e);
        const long long N = rnd(rng, 0, 1000000000000000000LL);
        SER_CHECK(arith_geom_sum(r, N + 1) - arith_geom_sum(r, N) == mint(N) * r.pow(N));
    }
}

template <class mint>
static void test_interpolate(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 5000 * TEST_SCALE; ++it, ++cases) {
        const int d = (int)rnd(rng, 0, 25);
        std::vector<mint> c((size_t)d + 1), ys((size_t)d + 1);
        for (auto& x : c) x = (long long)rng();
        for (int i = 0; i <= d; ++i) ys[(size_t)i] = poly_eval(c, mint(i));
        const long long x = it % 3 == 0 ? rnd(rng, -5, 40) : rnd_ll(rng);
        SER_CHECK(lagrange_interpolate(ys, x) == poly_eval(c, mint(x)));
        // ばらばらの点
        const int n = (int)rnd(rng, 1, 20);
        std::vector<mint> xs, vs;
        while ((int)xs.size() < n) {
            const mint v = it % 2 ? mint(rnd(rng, -30, 30)) : mint((long long)rng());
            if (std::find(xs.begin(), xs.end(), v) == xs.end()) xs.push_back(v), vs.push_back(mint((long long)rng()));
        }
        const std::vector<mint> co = poly_from_points(xs, vs);
        SER_CHECK((int)co.size() == n);
        for (int i = 0; i < n; ++i) SER_CHECK(poly_eval(co, xs[(size_t)i]) == vs[(size_t)i]);
        const mint t = (long long)rng();
        SER_CHECK(lagrange_interpolate_at(xs, vs, t) == poly_eval(co, t));
    }
    SER_CHECK(lagrange_interpolate(std::vector<mint>{}, 5) == mint(0) && lagrange_interpolate(std::vector<mint>{mint(7)}, 1000000000000000000LL) == mint(7));
}

static void test_poly(std::mt19937_64& rng, int& cases) {
    using mint = ModInt998;
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 12), m = (int)rnd(rng, 1, 8);
        std::vector<mint> a((size_t)n), b((size_t)m);
        for (auto& x : a) x = it % 3 ? mint((long long)rng()) : mint(rnd(rng, 0, 2));
        for (auto& x : b) x = it % 3 ? mint((long long)rng()) : mint(rnd(rng, 0, 2));
        if (b.back() == mint(0)) b.back() = 1;
        // 掛け算
        const std::vector<mint> c = poly_mul_naive(a, b);
        const mint x = (long long)rng();
        if (!a.empty()) SER_CHECK(poly_eval(c, x) == poly_eval(a, x) * poly_eval(b, x) && c.size() == a.size() + b.size() - 1);
        // 割り算：a = q b + r、deg r < deg b
        const PolyDivMod<mint> d = poly_divmod(a, b);
        std::vector<mint> back = poly_mul_naive(d.q, b);
        back.resize(std::max(back.size(), a.size()), mint(0));
        for (size_t i = 0; i < d.r.size(); ++i) back[i] += d.r[i];
        while (!back.empty() && back.back() == mint(0)) back.pop_back();
        std::vector<mint> at = a;
        while (!at.empty() && at.back() == mint(0)) at.pop_back();
        SER_CHECK(back == at && d.r.size() < b.size() && (d.r.empty() || d.r.back() != mint(0)));
        // 整数：b の最高次の係数が ±1 か、全部が割り切れるように作る
        std::vector<long long> qi((size_t)rnd(rng, 1, 6)), bi((size_t)m), ri((size_t)(m - 1));
        for (auto& v : qi) v = rnd(rng, -100, 100);
        for (auto& v : bi) v = rnd(rng, -100, 100);
        for (auto& v : ri) v = rnd(rng, -100, 100);
        bi.back() = it % 2 ? 1 : -1;
        if (qi.back() == 0) qi.back() = 3;
        std::vector<long long> ai = poly_mul_naive(qi, bi);
        for (size_t i = 0; i < ri.size(); ++i) ai[i] += ri[i];
        const PolyDivMod<long long> di = poly_divmod(ai, bi);
        while (!ri.empty() && ri.back() == 0) ri.pop_back();
        SER_CHECK(di.q == qi && di.r == ri);
    }
    // double
    const PolyDivMod<double> dd = poly_divmod<double>({-1, 0, 0, 1}, {-1, 1});  // x^3 − 1 = (x − 1)(x^2 + x + 1)
    SER_CHECK(dd.q == std::vector<double>({1, 1, 1}) && dd.r.empty());
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_sum_powers<ModInt998>(rng, cases, 12);
    test_sum_powers<ModInt107>(rng, cases, 12);
    test_sum_powers<ModInt<101>>(rng, cases, 12);
    test_arith_geom<ModInt998>(rng, cases);
    DynModInt::set_mod(1000000000);  // 合成数の mod でも使える
    test_arith_geom<DynModInt>(rng, cases);
    test_interpolate<ModInt998>(rng, cases);
    test_interpolate<ModInt107>(rng, cases);
    test_poly(rng, cases);
    std::cout << "[series] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[series] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    using mint = ModInt107;
    auto t0 = std::chrono::steady_clock::now();
    const mint s1 = sum_powers<mint>(1000000000000000000LL, 2500000);
    const long long t_sp = ms(t0);
    std::vector<mint> ys(1000001);
    for (auto& y : ys) y = (long long)rng();
    t0 = std::chrono::steady_clock::now();
    const mint s2 = lagrange_interpolate(ys, 1000000000000000000LL);
    const long long t_li = ms(t0);
    std::vector<mint> xs(2000), vs(2000);
    for (int i = 0; i < 2000; ++i) xs[(size_t)i] = mint(i * 7 + 3), vs[(size_t)i] = (long long)rng();
    t0 = std::chrono::steady_clock::now();
    const std::vector<mint> co = poly_from_points(xs, vs);
    const long long t_pp = ms(t0);
    std::vector<mint> a(4000), b(2000);
    for (auto& x : a) x = (long long)rng();
    for (auto& x : b) x = (long long)rng();
    t0 = std::chrono::steady_clock::now();
    const PolyDivMod<mint> d = poly_divmod(a, b);
    const long long t_div = ms(t0);
    std::cout << "[series] benchmark: sum_powers(10^18, k = 2.5×10^6) " << t_sp << " ms, lagrange_interpolate（10^6 + 1 点）" << t_li << " ms, poly_from_points（2000 点）" << t_pp
              << " ms, poly_divmod（4000 ÷ 2000）" << t_div << " ms  (check " << (s1 + s2 + co[5] + d.q[3]).val() % 7 << ")\n";
}

}  // namespace ser_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace ser_demo {

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
#define SHOW(expr, comment) ser_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    ser_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("べき乗和・ラグランジュ補間・小さい多項式（series） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  sum_powers<mint>(n, k) / sum_powers<mint>(l, r, k) / sum_powers_closed<mint>(L, R, k)   べき乗和\n"
              << "  arith_geom_sum(r, n)                       Σ i·r^i\n"
              << "  lagrange_interpolate(ys, x) / lagrange_interpolate_at(xs, ys, x) / poly_from_points(xs, ys)   補間\n"
              << "  poly_eval / poly_mul_naive / poly_divmod   次数の小さい多項式\n";
    note("using mint = ModInt107;");
    using mint = ModInt107;

    api("【1】 sum_powers<mint>(n, k) / sum_powers<mint>(l, r, k) / sum_powers_closed<mint>(L, R, k)", "mint",
        "Σ i^k（[0, n) / [l, r) / [L, R]）。k + 2 点の値を作ってラグランジュ補間", "O(k)");
    SHOW(sum_powers<mint>(11, 2), "0^2 + 1^2 + … + 10^2 = 385");
    SHOW(sum_powers_closed<mint>(1, 1000000000000000000LL, 3), "1^3 + … + (10^18)^3");
    SHOW(sum_powers<mint>(5, 8, 1), "5 + 6 + 7");
    SHOW(sum_powers<mint>(3, 0), "0^0 + 1^0 + 2^0（0^0 = 1）");

    api("【2】 arith_geom_sum(r, n)", "mint", "Σ_{i=0}^{n−1} i·r^i（等差 × 等比。mod は何でもよい）", "O(log n)");
    SHOW(arith_geom_sum(mint(2), 4), "0 + 2 + 8 + 24 = 34");
    SHOW(arith_geom_sum(mint(3), 1000000000000000000LL), "n = 10^18");

    api("【3】 lagrange_interpolate(ys, x)", "mint", "f(0)〜f(d) が ys の d 次以下の多項式の f(x)（x は 10^18 でもよい）", "O(d)");
    SHOW(lagrange_interpolate(std::vector<mint>{1, 3, 7}, 3), "ARC033 D 例 1：f(x) = x^2 + x + 1 → f(3) = 13");
    SHOW(lagrange_interpolate(std::vector<mint>{4, 16, 106, 484, 1624, 4384}, 1000000000), "ARC033 D 例 2");
    {
        // ABC208 F：f(n, m) は n の K + M 次の多項式。n = 0〜K+M の値を DP で作って補間
        auto solve = [](long long N, int M, int K) {
            const int D = K + M + 1;
            std::vector<mint> f((size_t)D);
            for (int n = 1; n < D; ++n) f[(size_t)n] = mint(n).pow(K);  // m = 0
            for (int m = 1; m <= M; ++m)
                for (int n = 1; n < D; ++n) f[(size_t)n] += f[(size_t)n - 1];  // f(n, m) = Σ_{j ≤ n} f(j, m − 1)
            return lagrange_interpolate(f, N);
        };
        note("ABC208 F（Cumulative Sum）: n = 0〜K+M の f(n, M) を累積和で作り、lagrange_interpolate(f, N)");
        std::cout << "    (3, 4, 2) → " << solve(3, 4, 2) << "、(0, 1, 2) → " << solve(0, 1, 2) << "、(10^18, 30, 123456) → " << solve(1000000000000000000LL, 30, 123456) << "\n";
    }

    api("【4】 lagrange_interpolate_at(xs, ys, x) / poly_from_points(xs, ys)", "mint / vector<mint>",
        "ばらばらの点 (xs[i], ys[i]) を通る多項式の f(x) / 係数（低い次数から）", "O(d^2)");
    RUN(std::vector<mint> xs = {1, 2, 4}, ys = {2, 5, 17});
    SHOW(poly_from_points(xs, ys), "f(x) = x^2 + 1 → 係数 [1, 0, 1]");
    SHOW(lagrange_interpolate_at(xs, ys, mint(10)), "f(10) = 101");

    api("【5】 poly_eval(c, x) / poly_mul_naive(a, b) / poly_divmod(a, b)", "T / vector<T> / PolyDivMod{q, r}",
        "係数は低い次数から。T は mint・long long・double。割り算は a = q b + r（整数なら割り切れないと止める）", "O(n) / O(nm)");
    SHOW(poly_eval<long long>({1, 0, 1}, 3), "x^2 + 1 で x = 3");
    SHOW(poly_mul_naive<long long>({2, 1}, {6, 4, 2}), "(x + 2)(2x^2 + 4x + 6)");
    RUN(PolyDivMod<long long> d = poly_divmod<long long>({5, 0, 0, 1}, {1, 1}));
    SHOW(d.q, "(x^3 + 5) ÷ (x + 1) の商 x^2 − x + 1");
    SHOW(d.r, "余り 4");
    note("ABC245 D（Polynomial division）: B = poly_divmod(C, A).q（係数は整数）");
    std::cout << "    例 1 → " << to_s(poly_divmod<long long>({12, 14, 8, 2}, {2, 1}).q) << "、例 2 → " << to_s(poly_divmod<long long>({10000, 0, -1}, {100, 1}).q) << "\n";
}

}  // namespace ser_demo

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
    ser_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ser_test::benchmark();    // 速度計測
    ser_demo::demo();         // 使い方デモ
    return 0;
}
