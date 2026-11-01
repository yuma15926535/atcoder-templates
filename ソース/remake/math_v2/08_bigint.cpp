// =============================================================================
//  多倍長整数（BigInt：long long・Int128 に収まらない整数の四則演算・比較・累乗・平方根・gcd・lcm・入出力）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 08 有理数・多倍長
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  使い方
//      BigInt a("123456789012345678901234567890"), b = 10;   // 文字列・整数から（cin >> a でもよい）
//      BigInt c = a * b + 1;                                   // int・long long と混ぜて書ける
//      cout << c / 7 << ' ' << c % 7 << '\n';
//  割り算は C++ の整数と同じ（0 の方向に切り捨て。% の符号は割られる数と同じ）
//  mod 998244353 などで足りるなら mint の方がずっと速い。答えが Int128 に収まるなら Int128 で
//
//  機能一覧                                                                    計算量（n, m は桁数）
//    BigInt(long long / Int128 / "文字列")   作る（文字列は先頭に + − を付けてよい）
//    + − × ÷ % 単項 −、+= など、++ −−、== < など、compare(a, b)
//                                     × は短ければ筆算、長ければカラツバ（10^5 桁同士で約 0.15 秒、10^6 桁同士は数秒）  × O(n^1.58)、÷ O(n m)
//    divmod(a, b)                     {商, 余り}
//    pow(a, e) / isqrt(a) / gcd(a, b) / abs(a)    累乗 / ⌊√a⌋（O(n^2)）/ 最大公約数（互除法 O(n^2)。2×10^4 桁で約 1 秒）/ 絶対値
//    lcm(a, b)                        最小公倍数（0 以上。どちらかが 0 なら 0。整数と混ぜてもよい）
//    a.to_string() / cout << a / cin >> a
//    a.fits_ll() / a.to_ll() / a.to_i128()   long long に収まるか / 変換（収まらなければ止める）
//    a.to_ld()                        おおよその値（long double）
//    a.num_digits() / a.digit_sum()   10 進の桁数 / 各桁の和
//    a.mod_ll(m)                      a mod m（[0, m)。巨大な数の mod を取るとき）              O(n)
//    a.sign() / a.is_zero() / a.is_even()
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    big_test::random_test()  … Int128 の計算・mod での検算との突き合わせ
//    big_test::benchmark()    … 速度計測
//    big_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
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
#ifndef GG_BIGINT
#define GG_BIGINT
// ---- 多倍長整数 BigInt（数学第二世代）  符号つき。10^9 進で持つ。× は桁が多いとカラツバ、÷ は筆算（Knuth） ----
// 割り算は C++ の整数と同じ（0 の方向に切り捨て。% の符号は割られる数と同じ）。0 で割ると止める
struct BigInt {
    static constexpr unsigned BASE = 1000000000U;
    static constexpr int BASE_DIGITS = 9;

    BigInt() {}
    template <class T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    BigInt(T v) {  // int・long long・unsigned など（Int128 も）
        if constexpr (std::is_signed<T>::value) {
            neg_ = v < 0;
            set_mag(neg_ ? (UInt128)0 - (UInt128)(Int128)v : (UInt128)(Int128)v);
        } else {
            set_mag((UInt128)v);
        }
    }
    BigInt(Int128 v) : neg_(v < 0) { set_mag(v < 0 ? (UInt128)0 - (UInt128)v : (UInt128)v); }
    BigInt(UInt128 v) { set_mag(v); }
    explicit BigInt(const std::string& s) { parse(s); }
    explicit BigInt(const char* s) { parse(std::string(s)); }

    // ---- 変換 ----
    std::string to_string() const {
        if (d_.empty()) return "0";
        std::string s = neg_ ? "-" : "";
        s += std::to_string(d_.back());
        char buf[16];
        for (size_t i = d_.size() - 1; i-- > 0;) {
            std::snprintf(buf, sizeof(buf), "%09u", d_[i]);
            s += buf;
        }
        return s;
    }
    bool fits_ll() const { return cmp_mag_u128((UInt128)9223372036854775807ULL + (neg_ ? 1 : 0)) <= 0; }
    long long to_ll() const {  // long long に収まらなければ止める
        if (!fits_ll()) math_internal::fail("BigInt::to_ll", "long long に収まらない");
        const UInt128 m = mag_u128();
        return neg_ ? (long long)((UInt128)0 - m) : (long long)m;
    }
    bool fits_i128() const { return cmp_mag_u128(((UInt128)1 << 127) - (neg_ ? 0 : 1)) <= 0; }
    Int128 to_i128() const {
        if (!fits_i128()) math_internal::fail("BigInt::to_i128", "Int128 に収まらない");
        const UInt128 m = mag_u128();
        return neg_ ? (Int128)((UInt128)0 - m) : (Int128)m;
    }
    long double to_ld() const {  // おおよその値
        long double r = 0;
        for (size_t i = d_.size(); i-- > 0;) r = r * BASE + d_[i];
        return neg_ ? -r : r;
    }

    // ---- 性質 ----
    bool is_zero() const { return d_.empty(); }
    int sign() const { return d_.empty() ? 0 : (neg_ ? -1 : 1); }
    bool is_even() const { return d_.empty() || d_[0] % 2 == 0; }
    int num_digits() const {  // 10 進の桁数（0 は 1 桁。符号は数えない）
        if (d_.empty()) return 1;
        int k = 0;
        for (unsigned t = d_.back(); t; t /= 10) ++k;
        return k + BASE_DIGITS * (int)(d_.size() - 1);
    }
    long long digit_sum() const {  // 10 進の各桁の和（符号は無視）
        long long s = 0;
        for (unsigned x : d_)
            for (; x; x /= 10) s += x % 10;
        return s;
    }
    long long mod_ll(long long m) const {  // この数を m（≥ 1）で割った余り（[0, m)。負の数も）  O(桁数)
        if (m < 1) math_internal::fail("BigInt::mod_ll", "m ≥ 1 で使う");
        UInt128 r = 0;
        for (size_t i = d_.size(); i-- > 0;) r = (r * BASE + d_[i]) % (UInt128)m;
        return neg_ && r ? m - (long long)r : (long long)r;
    }

    // ---- 比較 ----
    friend int compare(const BigInt& a, const BigInt& b) {  // −1 / 0 / 1
        if (a.neg_ != b.neg_) return a.neg_ ? -1 : 1;
        const int c = cmp_mag(a.d_, b.d_);
        return a.neg_ ? -c : c;
    }
    friend bool operator==(const BigInt& a, const BigInt& b) { return a.neg_ == b.neg_ && a.d_ == b.d_; }
    friend bool operator!=(const BigInt& a, const BigInt& b) { return !(a == b); }
    friend bool operator<(const BigInt& a, const BigInt& b) { return compare(a, b) < 0; }
    friend bool operator>(const BigInt& a, const BigInt& b) { return compare(a, b) > 0; }
    friend bool operator<=(const BigInt& a, const BigInt& b) { return compare(a, b) <= 0; }
    friend bool operator>=(const BigInt& a, const BigInt& b) { return compare(a, b) >= 0; }

    // ---- 四則演算 ----
    BigInt operator-() const {
        BigInt r = *this;
        if (!r.d_.empty()) r.neg_ = !r.neg_;
        return r;
    }
    BigInt operator+() const { return *this; }
    friend BigInt operator+(const BigInt& a, const BigInt& b) {
        if (a.neg_ == b.neg_) return make(add_mag(a.d_, b.d_), a.neg_);
        const int c = cmp_mag(a.d_, b.d_);
        if (c == 0) return BigInt();
        return c > 0 ? make(sub_mag(a.d_, b.d_), a.neg_) : make(sub_mag(b.d_, a.d_), b.neg_);
    }
    friend BigInt operator-(const BigInt& a, const BigInt& b) { return a + (-b); }
    friend BigInt operator*(const BigInt& a, const BigInt& b) { return make(mul_mag(a.d_, b.d_), a.neg_ != b.neg_); }
    friend BigInt operator/(const BigInt& a, const BigInt& b) { return divmod(a, b).first; }
    friend BigInt operator%(const BigInt& a, const BigInt& b) { return divmod(a, b).second; }
    BigInt& operator+=(const BigInt& b) { return *this = *this + b; }
    BigInt& operator-=(const BigInt& b) { return *this = *this - b; }
    BigInt& operator*=(const BigInt& b) { return *this = *this * b; }
    BigInt& operator/=(const BigInt& b) { return *this = *this / b; }
    BigInt& operator%=(const BigInt& b) { return *this = *this % b; }
    BigInt& operator++() { return *this += 1; }
    BigInt& operator--() { return *this -= 1; }
    BigInt operator++(int) {
        BigInt t = *this;
        *this += 1;
        return t;
    }
    BigInt operator--(int) {
        BigInt t = *this;
        *this -= 1;
        return t;
    }
    // 商と余り（0 の方向に切り捨て）  O(桁数 × 割る数の桁数)
    friend std::pair<BigInt, BigInt> divmod(const BigInt& a, const BigInt& b) {
        if (b.d_.empty()) math_internal::fail("BigInt::/", "0 で割った");
        std::vector<unsigned> q, r;
        div_mag(a.d_, b.d_, q, r);
        return {make(std::move(q), a.neg_ != b.neg_), make(std::move(r), a.neg_)};
    }

    // ⌊√n⌋（n ≥ 0）：上の半分の桁の平方根を初めの値にしたニュートン法（割り算 3 回ほど）  O(n^2)
    friend BigInt isqrt(const BigInt& n) {
        if (n.neg_) math_internal::fail("isqrt(BigInt)", "n ≥ 0 で使う");
        if (n.d_.size() <= 4) {  // 10^36 未満は UInt128 で
            const UInt128 v = n.mag_u128();
            UInt128 r = (UInt128)std::sqrt((long double)v);
            while (r > 0 && r * r > v) --r;
            while ((r + 1) * (r + 1) <= v) ++r;
            return BigInt(r);
        }
        const size_t s = n.d_.size() / 4;  // 下の 2s 桁（10^9 進）を落とした数の平方根 + 1 を BASE^s 倍すると √n 以上
        BigInt hi;
        hi.d_.assign(n.d_.begin() + (long long)(2 * s), n.d_.end());
        BigInt x = isqrt(hi) + 1;
        x.d_.insert(x.d_.begin(), s, 0U);
        while (true) {  // 上から近づく
            BigInt y = (x + n / x) / 2;
            if (y >= x) return x;
            x = std::move(y);
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const BigInt& a) { return os << a.to_string(); }
    friend std::istream& operator>>(std::istream& is, BigInt& a) {
        std::string s;
        if (is >> s) a = BigInt(s);
        return is;
    }

  private:
    std::vector<unsigned> d_;  // 10^9 進の桁（下の桁から。上の 0 は持たない。0 は空）
    bool neg_ = false;         // 0 のときは false

    void set_mag(UInt128 m) {
        d_.clear();
        for (; m; m /= BASE) d_.push_back((unsigned)(m % BASE));
        if (d_.empty()) neg_ = false;
    }
    UInt128 mag_u128() const {
        UInt128 r = 0;
        for (size_t i = d_.size(); i-- > 0;) r = r * BASE + d_[i];
        return r;
    }
    int cmp_mag_u128(UInt128 x) const {  // |this| と x の比較
        if (d_.size() > 5) return 1;  // 10^45 > 2^128
        BigInt t((UInt128)x);
        return cmp_mag(d_, t.d_);
    }
    void parse(const std::string& s) {
        size_t i = 0;
        bool ng = false;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) ng = s[i] == '-', ++i;
        if (i == s.size()) math_internal::fail("BigInt", "数字がない");
        for (size_t j = i; j < s.size(); ++j)
            if (s[j] < '0' || s[j] > '9') math_internal::fail("BigInt", "数として読めない文字がある");
        d_.clear();
        for (size_t end = s.size(); end > i;) {
            const size_t beg = end >= i + BASE_DIGITS ? end - BASE_DIGITS : i;
            unsigned x = 0;
            for (size_t j = beg; j < end; ++j) x = x * 10 + (unsigned)(s[j] - '0');
            d_.push_back(x), end = beg;
        }
        trim(d_);
        neg_ = ng && !d_.empty();
    }
    static void trim(std::vector<unsigned>& a) {
        while (!a.empty() && a.back() == 0) a.pop_back();
    }
    static BigInt make(std::vector<unsigned> m, bool ng) {
        trim(m);
        BigInt r;
        r.neg_ = ng && !m.empty();
        r.d_ = std::move(m);
        return r;
    }
    static int cmp_mag(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {
        if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
        for (size_t i = a.size(); i-- > 0;)
            if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
        return 0;
    }
    static std::vector<unsigned> add_mag(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {
        std::vector<unsigned> r(std::max(a.size(), b.size()) + 1, 0);
        unsigned c = 0;
        for (size_t i = 0; i + 1 < r.size(); ++i) {
            unsigned s = c + (i < a.size() ? a[i] : 0) + (i < b.size() ? b[i] : 0);
            c = s >= BASE;
            r[i] = c ? s - BASE : s;
        }
        r.back() = c;
        return r;
    }
    static std::vector<unsigned> sub_mag(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {  // |a| ≥ |b|
        std::vector<unsigned> r(a.size());
        long long br = 0;
        for (size_t i = 0; i < a.size(); ++i) {
            long long t = (long long)a[i] - br - (i < b.size() ? b[i] : 0);
            br = t < 0;
            r[i] = (unsigned)(br ? t + BASE : t);
        }
        return r;
    }
    // ---- 掛け算：短い方が 40 桁（10^9 進）以下なら筆算、それより長ければ 10^6 進にしてカラツバ ----
    static std::vector<unsigned> mul_mag(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {
        if (a.empty() || b.empty()) return {};
        if (std::min(a.size(), b.size()) <= 40) {
            std::vector<unsigned> r(a.size() + b.size(), 0);
            for (size_t i = 0; i < a.size(); ++i) {
                unsigned long long c = 0;
                for (size_t j = 0; j < b.size() || c; ++j) {
                    const unsigned long long cur = r[i + j] + c + (j < b.size() ? (unsigned long long)a[i] * b[j] : 0);
                    r[i + j] = (unsigned)(cur % BASE), c = cur / BASE;
                }
            }
            return r;
        }
        if (a.size() < b.size()) return mul_mag(b, a);
        if (a.size() >= 2 * b.size()) {  // 長さが大きく違う：長い方を短い方の長さで区切って掛け、ずらして足す
            std::vector<unsigned> r(a.size() + b.size() + 1, 0);
            for (size_t off = 0; off < a.size(); off += b.size()) {
                std::vector<unsigned> ch(a.begin() + (long long)off, a.begin() + (long long)std::min(off + b.size(), a.size()));
                trim(ch);
                std::vector<unsigned> pr = mul_mag(ch, b);
                trim(pr);  // カラツバの結果は後ろに 0 が付いていることがある
                unsigned c = 0;
                size_t i = 0;
                for (; i < pr.size() || c; ++i) {
                    const unsigned t = r[off + i] + c + (i < pr.size() ? pr[i] : 0);
                    c = t >= BASE;
                    r[off + i] = c ? t - BASE : t;
                }
            }
            return r;
        }
        // 10^k 進（k = 6、桁がとても多ければ 3）に直してカラツバ。途中の値が long long に収まるように k を選ぶ
        //   長さ n（2 の冪）のとき途中の値は n^2 (10^k)^2 / 32 以下：k = 6 なら n ≤ 2^14、k = 3 なら n ≤ 2^23
        int k = 6;
        size_t n = pow2_len(std::max(a.size(), b.size()), k);
        if (n > ((size_t)1 << 14)) k = 3, n = pow2_len(std::max(a.size(), b.size()), k);
        if (n > ((size_t)1 << 23)) math_internal::fail("BigInt::*", "桁が多すぎる（2×10^7 桁程度まで）");
        std::vector<long long> x = to_small(a, k), y = to_small(b, k);
        x.resize(n), y.resize(n);
        const std::vector<long long> z = karatsuba(x, y);
        const long long S = k == 6 ? 1000000 : 1000;
        std::vector<long long> ds;  // 10^k 進で繰り上げ
        long long c = 0;
        for (size_t i = 0; i < z.size() || c; ++i) {
            const long long cur = (i < z.size() ? z[i] : 0) + c;
            ds.push_back(cur % S), c = cur / S;
        }
        const size_t g = (size_t)(18 / k);  // 10^k 進 g 桁 = 10^18 未満 = 10^9 進 2 桁
        std::vector<unsigned> r;
        for (size_t i = 0; i < ds.size(); i += g) {
            unsigned long long v = 0;
            for (size_t t = g; t-- > 0;) v = v * (unsigned long long)S + (i + t < ds.size() ? (unsigned long long)ds[i + t] : 0);
            r.push_back((unsigned)(v % BASE)), r.push_back((unsigned)(v / BASE));
        }
        return r;
    }
    static size_t pow2_len(size_t limbs, int k) {  // 10^9 進 limbs 桁を 10^k 進にした長さ以上の 2 の冪
        const size_t len = (limbs + 1) / 2 * (size_t)(18 / k);
        size_t n = 1;
        while (n < len) n *= 2;
        return n;
    }
    static std::vector<long long> to_small(const std::vector<unsigned>& a, int k) {  // 10^9 進 2 桁（< 10^18）→ 10^k 進 18/k 桁
        const unsigned long long S = k == 6 ? 1000000 : 1000;
        std::vector<long long> r;
        for (size_t i = 0; i < a.size(); i += 2) {
            unsigned long long v = a[i] + (i + 1 < a.size() ? (unsigned long long)a[i + 1] * BASE : 0);
            for (int t = 0; t < 18 / k; ++t) r.push_back((long long)(v % S)), v /= S;
        }
        return r;
    }
    static std::vector<long long> karatsuba(const std::vector<long long>& a, const std::vector<long long>& b) {  // 長さは同じ 2 の冪。結果は長さ 2n
        const size_t n = a.size();
        std::vector<long long> r(2 * n, 0);
        if (n <= 32) {
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j) r[i + j] += a[i] * b[j];
            return r;
        }
        const size_t k = n / 2;
        const std::vector<long long> a1(a.begin(), a.begin() + (long long)k), a2(a.begin() + (long long)k, a.end());
        const std::vector<long long> b1(b.begin(), b.begin() + (long long)k), b2(b.begin() + (long long)k, b.end());
        const std::vector<long long> a1b1 = karatsuba(a1, b1), a2b2 = karatsuba(a2, b2);
        std::vector<long long> s1(k), s2(k);
        for (size_t i = 0; i < k; ++i) s1[i] = a1[i] + a2[i], s2[i] = b1[i] + b2[i];
        std::vector<long long> m = karatsuba(s1, s2);
        for (size_t i = 0; i < 2 * k; ++i) m[i] -= a1b1[i] + a2b2[i];
        for (size_t i = 0; i < 2 * k; ++i) r[i] += a1b1[i], r[i + k] += m[i], r[i + n] += a2b2[i];
        return r;
    }
    // ---- 割り算（Knuth のアルゴリズム D）----
    static void div_mag(const std::vector<unsigned>& a, const std::vector<unsigned>& b, std::vector<unsigned>& q, std::vector<unsigned>& r) {
        if (cmp_mag(a, b) < 0) {
            q.clear(), r = a;
            return;
        }
        if (b.size() == 1) {  // 1 桁で割る
            q.assign(a.size(), 0);
            unsigned long long rem = 0;
            for (size_t i = a.size(); i-- > 0;) {
                const unsigned long long cur = rem * BASE + a[i];
                q[i] = (unsigned)(cur / b[0]), rem = cur % b[0];
            }
            r.clear();
            if (rem) r.push_back((unsigned)rem);
            trim(q);
            return;
        }
        const size_t n = b.size(), m = a.size() - b.size();
        const unsigned long long norm = BASE / ((unsigned long long)b.back() + 1);  // 割る数の最上位を BASE/2 以上に
        std::vector<unsigned> u(a.size() + 1, 0), v(n, 0);
        {
            unsigned long long c = 0;
            for (size_t i = 0; i < a.size(); ++i) {
                const unsigned long long cur = a[i] * norm + c;
                u[i] = (unsigned)(cur % BASE), c = cur / BASE;
            }
            u[a.size()] = (unsigned)c, c = 0;
            for (size_t i = 0; i < n; ++i) {
                const unsigned long long cur = b[i] * norm + c;
                v[i] = (unsigned)(cur % BASE), c = cur / BASE;
            }
        }
        q.assign(m + 1, 0);
        for (size_t j = m + 1; j-- > 0;) {
            const unsigned long long num = (unsigned long long)u[j + n] * BASE + u[j + n - 1];
            unsigned long long qh = num / v[n - 1], rh = num % v[n - 1];
            while (qh >= BASE || qh * v[n - 2] > rh * BASE + u[j + n - 2]) {
                --qh, rh += v[n - 1];
                if (rh >= BASE) break;
            }
            // u[j..j+n] −= qh × v
            long long br = 0;
            unsigned long long c = 0;
            for (size_t i = 0; i < n; ++i) {
                const unsigned long long p = qh * v[i] + c;
                c = p / BASE;
                long long t = (long long)u[i + j] - (long long)(p % BASE) - br;
                br = t < 0;
                u[i + j] = (unsigned)(br ? t + BASE : t);
            }
            long long top = (long long)u[j + n] - (long long)c - br;
            if (top < 0) {  // 引きすぎ：1 回分戻す
                --qh;
                unsigned long long cc = 0;
                for (size_t i = 0; i < n; ++i) {
                    const unsigned long long s = (unsigned long long)u[i + j] + v[i] + cc;
                    u[i + j] = (unsigned)(s % BASE), cc = s / BASE;
                }
                top += (long long)cc;
            }
            u[j + n] = (unsigned)top;
            q[j] = (unsigned)qh;
        }
        trim(q);
        r.assign(n, 0);  // 余り = u の下 n 桁 / norm
        unsigned long long rem = 0;
        for (size_t i = n; i-- > 0;) {
            const unsigned long long cur = rem * BASE + u[i];
            r[i] = (unsigned)(cur / norm), rem = cur % norm;
        }
        trim(r);
    }
};
inline BigInt abs(const BigInt& a) { return a.sign() < 0 ? -a : a; }
inline BigInt pow(BigInt a, long long e) {  // a^e（e ≥ 0。0^0 = 1）
    if (e < 0) math_internal::fail("pow(BigInt)", "e ≥ 0 で使う");
    BigInt r = 1;
    for (; e > 0; e >>= 1) {
        if (e & 1) r *= a;
        if (e > 1) a *= a;
    }
    return r;
}
inline BigInt gcd(BigInt a, BigInt b) {  // 最大公約数（0 以上。gcd(0, 0) = 0）
    a = abs(a), b = abs(b);
    while (!b.is_zero()) {
        BigInt t = a % b;
        a = std::move(b), b = std::move(t);
    }
    return a;
}
template <class T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
inline BigInt gcd(const BigInt& a, T b) { return gcd(a, BigInt(b)); }
template <class T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
inline BigInt gcd(T a, const BigInt& b) { return gcd(BigInt(a), b); }
inline BigInt lcm(const BigInt& a, const BigInt& b) {  // 最小公倍数（0 以上。どちらかが 0 なら 0）
    if (a.is_zero() || b.is_zero()) return BigInt(0);
    return abs(a) / gcd(a, b) * abs(b);
}
template <class T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
inline BigInt lcm(const BigInt& a, T b) { return lcm(a, BigInt(b)); }
template <class T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
inline BigInt lcm(T a, const BigInt& b) { return lcm(BigInt(a), b); }
#endif  // GG_BIGINT
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace big_test {

#define BIG_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[bigint] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static Int128 rnd128(std::mt19937_64& rng) {  // 大きさがばらばらの Int128
    const int bits = (int)rnd(rng, 0, 126);
    UInt128 u = ((UInt128)rng() << 64) | rng();
    u = bits == 0 ? 0 : u >> (128 - bits);
    return (rng() & 1) ? -(Int128)u : (Int128)u;
}
static std::string rand_digits(std::mt19937_64& rng, int len) {  // len 桁（先頭は 0 でない。符号もランダム）
    std::string s = (rng() & 1) ? "-" : "";
    s += (char)('1' + rng() % 9);
    for (int i = 1; i < len; ++i) s += (char)('0' + (rng() % 5 == 0 ? (rng() % 2) * 9 : rng() % 10));  // 0 と 9 を多めに
    return s;
}
static const long long PRIMES[3] = {998244353, 1000000007, 4611686018427387847LL};

static void check_div(const BigInt& a, const BigInt& b) {
    const auto qr = divmod(a, b);
    const BigInt& q = qr.first;
    const BigInt& r = qr.second;
    BIG_CHECK(q * b + r == a && abs(r) < abs(b));
    BIG_CHECK(r.is_zero() || r.sign() == a.sign());
    BIG_CHECK(a / b == q && a % b == r);
}
static void check_mul(const BigInt& a, const BigInt& b) {
    const BigInt c = a * b;
    for (long long p : PRIMES) BIG_CHECK(c.mod_ll(p) == (long long)((Int128)a.mod_ll(p) * b.mod_ll(p) % p));
    BIG_CHECK(c.sign() == a.sign() * b.sign());
}

static void test_small(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const Int128 x = rnd128(rng), y = rnd128(rng);
        const BigInt a = x, b = y;
        BIG_CHECK(a.to_string() == int128_to_string(x) && BigInt(int128_to_string(x)) == a && a.to_i128() == x);
        BIG_CHECK((a < b) == (x < y) && (a == b) == (x == y) && (a >= b) == (x >= y) && compare(a, b) == (x < y ? -1 : x > y ? 1 : 0));
        // 足し算・引き算（Int128 で溢れない範囲だけ比べる：|x|, |y| < 2^126）
        BIG_CHECK((a + b).to_i128() == x + y && (a - b).to_i128() == x - y && (-a).to_i128() == -x);
        if (y != 0) BIG_CHECK((a / b).to_i128() == x / y && (a % b).to_i128() == x % y);
        const Int128 xs = x >> 64, ys = y >> 64;  // 掛け算は小さくして
        BIG_CHECK((BigInt(xs) * BigInt(ys)).to_i128() == xs * ys);
        const long long l = (long long)(x >> 64);
        BIG_CHECK(BigInt(l).fits_ll() && BigInt(l).to_ll() == l);
        BIG_CHECK(a.fits_ll() == (x >= LLONG_MIN && x <= LLONG_MAX));
        const long long m = rnd(rng, 1, 1000000000000000000LL);
        BIG_CHECK(a.mod_ll(m) == (long long)(((x % m) + m) % m));
        BIG_CHECK(a.num_digits() == (int)int128_to_string(x < 0 ? -x : x).size());
        if (x >= 0) {  // 平方根
            const BigInt s = isqrt(a);
            BIG_CHECK(s * s <= a && (s + 1) * (s + 1) > a);
        }
    }
    // 境目
    BIG_CHECK(BigInt(LLONG_MIN).to_ll() == LLONG_MIN && BigInt(LLONG_MAX).to_ll() == LLONG_MAX);
    BIG_CHECK(!(BigInt(LLONG_MAX) + 1).fits_ll() && !(BigInt(LLONG_MIN) - 1).fits_ll());
    const Int128 I128MIN = (Int128)((UInt128)1 << 127), I128MAX = (Int128)(((UInt128)1 << 127) - 1);
    BIG_CHECK(BigInt(I128MIN).to_i128() == I128MIN && BigInt(I128MAX).to_i128() == I128MAX && !(BigInt(I128MAX) + 1).fits_i128() && !(BigInt(I128MIN) - 1).fits_i128());
    BIG_CHECK(BigInt(I128MIN).to_string() == "-170141183460469231731687303715884105728");
    BIG_CHECK(BigInt((UInt128)-1).to_string() == "340282366920938463463374607431768211455" && BigInt(18446744073709551615ULL).to_string() == "18446744073709551615");
    BIG_CHECK(BigInt("-0") == BigInt(0) && BigInt("-0").sign() == 0 && BigInt("+000123") == 123 && BigInt("0000").to_string() == "0");
    BIG_CHECK((BigInt(-7) / 2) == -3 && (BigInt(-7) % 2) == -1 && (BigInt(7) / -2) == -3 && (BigInt(7) % -2) == 1);
    BIG_CHECK(BigInt(5) - 5 == 0 && (BigInt(5) - 5).sign() == 0 && -(BigInt(0)) == 0);
    BigInt t = 41;
    BIG_CHECK((++t) == 42 && (--t) == 41 && t.is_even() == false);
    cases += 10;
}

static void test_big(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int la = (int)rnd(rng, 1, it % 10 == 0 ? 2000 : 120), lb = (int)rnd(rng, 1, it % 3 == 0 ? la + 5 : 60);
        const BigInt a(rand_digits(rng, la)), b(rand_digits(rng, lb));
        BIG_CHECK(a.num_digits() == la && BigInt(a.to_string()) == a);
        check_mul(a, b);
        check_div(a, b);
        check_div(a * b + b / 3, b);
        BIG_CHECK((a * b) / b == a && (a * b) % b == 0 && (a + b) - b == a);
        // 9 で割った余りと各桁の和
        BIG_CHECK(a.digit_sum() % 9 == abs(a).mod_ll(9));
    }
    // 長い掛け算（カラツバ・10^3 進の方も）
    for (int len : {500, 3000, 20000, 100000, 120000}) {
        for (int t = 0; t < 2; ++t, ++cases) {
            const BigInt a(rand_digits(rng, len)), b(rand_digits(rng, len - t * 77));
            check_mul(a, b);
            if (len <= 20000) check_div(a * b + a, b);
        }
    }
    // 長さが大きく違う掛け算（短い方で区切る）
    for (int t = 0; t < 6; ++t, ++cases) {
        const BigInt a(rand_digits(rng, 50000 + t * 777)), b(rand_digits(rng, 300 + t * 401));
        check_mul(a, b);
        check_mul(b, a);
        BIG_CHECK(a * b == b * a && (a * b) / b == a);
    }
    // 大きい平方根・int との gcd・後置の ++ −−
    for (int len : {40, 41, 100, 1000, 5000, 20000}) {
        const BigInt n = abs(BigInt(rand_digits(rng, len)));
        const BigInt s = isqrt(n);
        BIG_CHECK(s * s <= n && (s + 1) * (s + 1) > n && isqrt(s * s) == s && isqrt(s * s - 1) == s - 1);
        ++cases;
    }
    {
        BigInt x = 5;
        BIG_CHECK(x++ == 5 && x == 6 && x-- == 6 && x == 5);
        BIG_CHECK(gcd(BigInt(12), 18) == 6 && gcd(30LL, BigInt(-45)) == 15 && gcd(BigInt(7), 0) == 7);
        cases += 2;
    }
    {  // 全部 9（繰り上がりがいちばん多い）
        const BigInt a(std::string(150000, '9'));
        const BigInt c = a * a;  // (10^n − 1)^2 = 10^{2n} − 2·10^n + 1
        BIG_CHECK(c == pow(BigInt(10), 300000) - 2 * pow(BigInt(10), 150000) + 1);
        ++cases;
    }
    // gcd・isqrt・pow
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const BigInt g(rand_digits(rng, (int)rnd(rng, 1, 40))), a(rand_digits(rng, (int)rnd(rng, 1, 60))), b(rand_digits(rng, (int)rnd(rng, 1, 60)));
        const BigInt h = gcd(a * g, b * g);
        BIG_CHECK(h.sign() > 0 && (a * g) % h == 0 && (b * g) % h == 0 && h % abs(g) == 0);
        BIG_CHECK(gcd((a * g) / h, (b * g) / h) == 1);
        const BigInt l = lcm(a * g, b * g);
        BIG_CHECK(l * h == abs(a * g * b * g) && (a.is_zero() || b.is_zero() || g.is_zero() || (l % (a * g) == 0 && l % (b * g) == 0)));
        const BigInt n = abs(a) * abs(b) + abs(g);
        const BigInt s = isqrt(n);
        BIG_CHECK(s * s <= n && (s + 1) * (s + 1) > n && isqrt(s * s) == s);
        const long long e = rnd(rng, 0, 30);
        BigInt p = 1;
        for (int i = 0; i < e; ++i) p *= g;
        BIG_CHECK(pow(g, e) == p);
    }
    BIG_CHECK(gcd(BigInt(0), BigInt(0)) == 0 && gcd(BigInt(0), BigInt(-12)) == 12 && pow(BigInt(0), 0) == 1);
    BIG_CHECK(lcm(BigInt(4), 6) == 12 && lcm(-4LL, BigInt(6)) == 12 && lcm(BigInt(0), 5) == 0 && lcm(BigInt(1000000007), BigInt(998244353)) == BigInt(1000000007LL * 998244353LL));
    // よく知られた値
    BIG_CHECK(pow(BigInt(2), 1000).digit_sum() == 1366);  // Project Euler 16
    BigInt f = 1;
    for (int i = 1; i <= 100; ++i) f *= i;
    BIG_CHECK(f.digit_sum() == 648 && f.num_digits() == 158);  // Project Euler 20
    BIG_CHECK(isqrt(pow(BigInt(10), 100)) == pow(BigInt(10), 50) && isqrt(pow(BigInt(10), 101)).to_string().substr(0, 10) == "3162277660");
    cases += 4;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_small(rng, cases);
    test_big(rng, cases);
    std::cout << "[bigint] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[bigint] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const BigInt a(rand_digits(rng, 100000)), b(rand_digits(rng, 100000)), c(rand_digits(rng, 50000));
    auto t0 = std::chrono::steady_clock::now();
    const BigInt p = a * b;
    const long long t_mul = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const BigInt q = p / c;
    const long long t_div = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const BigInt sq = isqrt(abs(p));
    const long long t_sqrt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    BigInt f = 1;
    for (int i = 1; i <= 10000; ++i) f *= i;
    const long long t_fact = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::string s = f.to_string();
    const BigInt g(s);
    const long long t_str = ms(t0);
    std::cout << "[bigint] benchmark: 10^5 桁 × 10^5 桁 " << t_mul << " ms, 2×10^5 桁 ÷ 5×10^4 桁 " << t_div << " ms, isqrt（2×10^5 桁）" << t_sqrt << " ms, 10000!（35660 桁）" << t_fact
              << " ms, 文字列との変換 " << t_str << " ms  (check " << (p.mod_ll(7) + q.mod_ll(7) + g.mod_ll(7) + sq.mod_ll(7)) % 7 << ")\n";
}

}  // namespace big_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace big_demo {

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
#define SHOW(expr, comment) big_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    big_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("多倍長整数（bigint） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  BigInt(long long / Int128 / \"文字列\")、+ − × ÷ %、比較、divmod / pow / isqrt / gcd / abs\n"
              << "  to_string / fits_ll / to_ll / to_i128 / to_ld / num_digits / digit_sum / mod_ll / sign / is_zero / is_even、cin >> / cout <<\n";

    api("【1】 BigInt a(\"…\") / BigInt b = 10 / + − × ÷ %", "BigInt", "符号つきの多倍長整数。int・long long と混ぜて書ける（÷ は 0 の方向に切り捨て）", "× O(n^1.58)、÷ O(n m)");
    RUN(BigInt a("123456789012345678901234567890"));
    SHOW(a * a, "");
    SHOW(a * 2 + 1, "");
    SHOW(a / 7, "");
    SHOW(a % 7, "");
    SHOW(BigInt(-7) / 2, "C++ と同じく 0 の方向へ");
    SHOW(BigInt(-7) % 2, "余りの符号は割られる数と同じ");
    SHOW(divmod(a, BigInt(1000000007)).second, "{商, 余り}");
    SHOW(a > BigInt(LLONG_MAX), "比較");

    api("【2】 pow(a, e) / isqrt(a) / gcd(a, b) / lcm(a, b) / abs(a)", "BigInt", "累乗 / ⌊√a⌋ / 最大公約数・最小公倍数（0 以上）/ 絶対値", "pow O(M(n) log e)");
    SHOW(pow(BigInt(2), 100), "2^100");
    SHOW(isqrt(BigInt("1000000000000000000000000000000")), "⌊√10^30⌋ = 10^15");
    SHOW(isqrt(pow(BigInt(10), 41)), "√10^41 = 316227766016837933199.88…");
    SHOW(gcd(pow(BigInt(6), 30), pow(BigInt(10), 20)), "2^20");
    {
        BigInt l = 1;
        for (int i = 1; i <= 100; ++i) l = lcm(l, i);
        std::cout << "  lcm(1, 2, …, 100)（BigInt で順に）-> " << l << "\n";
    }
    SHOW(abs(BigInt(-5)), "");

    api("【3】 to_string / fits_ll / to_ll / to_i128 / to_ld", "string / bool / long long / Int128 / long double", "変換（to_ll・to_i128 は収まらなければ止める）", "O(n)");
    SHOW(BigInt(LLONG_MAX).fits_ll(), "");
    SHOW((BigInt(LLONG_MAX) + 1).fits_ll(), "2^63 は long long に入らない");
    SHOW((BigInt(LLONG_MAX) + 1).to_i128(), "Int128 には入る");
    SHOW(pow(BigInt(3), 100).to_ld(), "おおよその値");

    api("【4】 num_digits / digit_sum / mod_ll(m)", "int / long long / long long", "10 進の桁数 / 各桁の和 / a mod m（[0, m)）", "O(n)");
    SHOW(pow(BigInt(2), 1000).num_digits(), "2^1000 は 302 桁");
    SHOW(pow(BigInt(2), 1000).digit_sum(), "各桁の和（Project Euler 16）");
    {
        BigInt f = 1;
        for (int i = 1; i <= 100; ++i) f *= i;
        note("100! = " + f.to_string().substr(0, 30) + "…（" + std::to_string(f.num_digits()) + " 桁、各桁の和 " + std::to_string(f.digit_sum()) + "）");
    }
    SHOW(BigInt(std::string(1000, '9')).mod_ll(998244353), "1000 桁の数 mod 998244353（文字列で来る巨大な N の mod）");
    note("比べるだけなら桁数 → 文字列の比較でもよい。答えが Int128（±1.7×10^38）に収まるなら Int128 の方が速い");
}

}  // namespace big_demo

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
    big_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    big_test::benchmark();    // 速度計測
    big_demo::demo();         // 使い方デモ
    return 0;
}
