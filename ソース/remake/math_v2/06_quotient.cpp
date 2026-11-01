// =============================================================================
//  商の列挙・床関数の和（⌊N/i⌋ が同じ i をまとめる、Σ⌊N/i⌋、Σ⌊(a i + b)/m⌋・Σ (a i + b) mod m・Σ⌊⌋²、直線の下の格子点、万能ユークリッド）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 06 和・数列
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  機能一覧                                                                    計算量
//    （N は 0 以上 LLONG_MAX 未満：区間の右端 r = N + 1 が long long に入るように）
//    for_each_quotient(N, f)          i = 1〜N を ⌊N/i⌋ = q が同じ区間 [l, r) に分けて f(q, l, r)（i の小さい順）   O(√N)
//    quotient_blocks(N)               上を vector<QuotientBlock{q, l, r}> で
//    quotient_blocks_ceil(N)          ⌈N/i⌉ = q が同じ区間に分ける
//    quotient_blocks2(N, M)           i = 1〜min(N, M) を ⌊N/i⌋・⌊M/i⌋ が両方同じ区間に（QuotientBlock2{q1, q2, l, r}）
//    sum_floor_div(N)                 Σ_{i=1}^{N} ⌊N/i⌋（1〜N の約数の個数の和。Int128）                   O(√N)
//    sum_i_floor_div(N)               Σ_{i=1}^{N} i·⌊N/i⌋（1〜N の約数の和の和。Int128）
//    sum_mod_i(N)                     Σ_{i=1}^{N} (N mod i)（Int128）
//    for_each_quotient で f が bool を返すと、false のところで止まる
//    linear_floor_sum(n, m, a, b)     Σ_{i=0}^{n−1} ⌊(a i + b) / m⌋（a, b は負でも。Int128。途中で溢れたら止める） O(log m)
//                                     （引数の順番は ACL の floor_sum と同じ：個数 n、割る数 m、傾き a、切片 b）
//    linear_mod_sum(n, m, a, b)       Σ_{i=0}^{n−1} ((a i + b) mod m)（a, b は負でも。Int128）                 O(log m)
//    linear_floor_sum_ext<mint>(n, m, a, b)   FloorSumExt{sum, sum_i, sum_sq}：Σ⌊⌋・Σ i⌊⌋・Σ⌊⌋^2 を mod で       O(log m)
//    universal_euclid(n, a, b, c, U, R)   万能ユークリッド：y = ⌊(a x + b)/c⌋ に沿って U（y が 1 増える）・R（x が 1 増える）を
//                                     掛けた積（x = 1〜n。a, b ≥ 0。M は結合的な * と単位元 M() を持つ型。可換でなくてよい）
//    count_lattice_under_line(a, b, c)    a x + b y ≤ c となる x, y ≥ 0 の組の数（a, b ≥ 1。Int128）      O(log)
//  （1〜N の約数の個数・約数の和の表は 02_sieve、数列の和の公式は 00_math_base の arith_sum / sum_range）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    quot_test::random_test()  … 1 つずつ数える愚直解との突き合わせ
//    quot_test::benchmark()    … 速度計測
//    quot_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
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
#ifndef GG_QUOTIENT
#define GG_QUOTIENT
// ---- 商の列挙・床関数の和（数学第二世代）  ⌊N/i⌋ が同じ i をまとめる、Σ⌊N/i⌋、Σ⌊(a i + b)/m⌋・Σ (a i + b) mod m・Σ⌊⌋²、万能ユークリッド ----
struct QuotientBlock {  // l ≤ i < r で ⌊N/i⌋（ceil 版は ⌈N/i⌉）が q
    long long q, l, r;
};
struct QuotientBlock2 {  // l ≤ i < r で ⌊N/i⌋ = q1、⌊M/i⌋ = q2
    long long q1, q2, l, r;
};

// i = 1, 2, …, N を ⌊N/i⌋ が同じ区間に分けて、i の小さい順に f(q, l, r)（区間の数は 2√N 以下）  O(√N)
//   f が bool を返すなら、false を返したところで止める
template <class F>
void for_each_quotient(long long N, F f) {  // （0 ≤ N < 2^63 − 1。区間の右端 r = N + 1 が long long に入るように）
    if (N < 0 || N == 9223372036854775807LL) math_internal::fail("for_each_quotient", "0 ≤ N < LLONG_MAX で使う");
    for (long long l = 1; l <= N;) {
        const long long q = N / l, r = N / q + 1;
        if constexpr (std::is_same<decltype(f(q, l, r)), bool>::value) {
            if (!f(q, l, r)) return;
        } else {
            f(q, l, r);
        }
        l = r;
    }
}
inline std::vector<QuotientBlock> quotient_blocks(long long N) {  // 上を vector で（q は大きい順）
    std::vector<QuotientBlock> v;
    for_each_quotient(N, [&](long long q, long long l, long long r) { v.push_back(QuotientBlock{q, l, r}); });
    return v;
}
inline std::vector<QuotientBlock> quotient_blocks_ceil(long long N) {  // i = 1〜N を ⌈N/i⌉ が同じ区間に分ける（N ≥ 1）
    if (N < 0 || N == 9223372036854775807LL) math_internal::fail("quotient_blocks_ceil", "0 ≤ N < LLONG_MAX で使う");
    std::vector<QuotientBlock> v;
    for (long long l = 1; l <= N;) {
        const long long q = (N - 1) / l + 1;              // ⌈N/l⌉（N ≥ 1。N + l − 1 は溢れることがある）
        const long long r = q == 1 ? N + 1 : (N - 1) / (q - 1) + 1;  // ⌈N/i⌉ = q となる最大の i は ⌊(N−1)/(q−1)⌋
        v.push_back(QuotientBlock{q, l, r});
        l = r;
    }
    return v;
}
inline std::vector<QuotientBlock2> quotient_blocks2(long long N, long long M) {  // i = 1〜min(N, M) を ⌊N/i⌋・⌊M/i⌋ が両方同じ区間に
    if (N < 0 || M < 0 || std::min(N, M) == 9223372036854775807LL) math_internal::fail("quotient_blocks2", "N, M ≥ 0、min(N, M) < LLONG_MAX で使う");
    std::vector<QuotientBlock2> v;
    const long long lim = std::min(N, M);
    for (long long l = 1; l <= lim;) {
        const long long q1 = N / l, q2 = M / l, r = std::min(N / q1, M / q2) + 1;
        v.push_back(QuotientBlock2{q1, q2, l, r});
        l = r;
    }
    return v;
}

// ---- ⌊N/i⌋ の和（N ≥ 0。答えは Int128） ----
inline Int128 sum_floor_div(long long N) {  // Σ_{i=1}^{N} ⌊N/i⌋（= 1〜N の約数の個数の和）  O(√N)
    if (N < 0) math_internal::fail("sum_floor_div", "N ≥ 0 で使う");
    const long long s = isqrt(N);
    Int128 r = 0;
    for (long long i = 1; i <= s; ++i) r += N / i;
    return 2 * r - (Int128)s * s;  // i ≤ √N と ⌊N/i⌋ ≤ √N の対称性
}
inline Int128 sum_i_floor_div(long long N) {  // Σ_{i=1}^{N} i·⌊N/i⌋（= 1〜N の約数の和の和）  O(√N)
    if (N < 0) math_internal::fail("sum_i_floor_div", "N ≥ 0 で使う");
    Int128 r = 0;
    for_each_quotient(N, [&](long long q, long long l, long long R) {
        const long long w = R - l;
        const Int128 t = (Int128)l + R - 1;  // l + … + (R − 1) = w·t / 2（どちらかは偶数。t は 2^63 を超えることがある）
        r += (Int128)q * (w % 2 == 0 ? (Int128)(w / 2) * t : (Int128)w * (t / 2));
    });
    return r;
}
inline Int128 sum_mod_i(long long N) {  // Σ_{i=1}^{N} (N mod i)  O(√N)
    return (Int128)N * N - sum_i_floor_div(N);
}

// ---- Σ_{i=0}^{n−1} ⌊(a i + b) / m⌋（a, b は負でもよい。n ≥ 0、m ≥ 1。答えは Int128。途中で溢れたら止める）  O(log m) ----
namespace math_internal {
inline Int128 qmul(Int128 x, Int128 y) {
    Int128 r;
    if (__builtin_mul_overflow(x, y, &r)) fail("linear_floor_sum", "Int128 に収まらない");
    return r;
}
inline Int128 qadd(Int128 x, Int128 y) {
    Int128 r;
    if (__builtin_add_overflow(x, y, &r)) fail("linear_floor_sum", "Int128 に収まらない");
    return r;
}
inline Int128 tri_count(Int128 n) { return (n % 2 == 0) ? qmul(n / 2, n - 1) : qmul(n, (n - 1) / 2); }  // n(n−1)/2
// a, b ≥ 0（ACL と同じユークリッドの互除法の形）。U は n, m < 2^32 なら unsigned long long（速い）、それ以外は Int128
template <class U>
inline Int128 floor_sum_nonneg(U n, U m, U a, U b) {
    Int128 ans = 0;
    while (true) {
        if constexpr (std::is_same_v<U, unsigned long long>) {  // n, m < 2^32 なら答えは 2^64 未満なので確認しない
            if (a >= m) ans += (Int128)(n * (n - 1) / 2) * (Int128)(a / m), a %= m;
            if (b >= m) ans += (Int128)n * (Int128)(b / m), b %= m;
        } else {
            if (a >= m) ans = qadd(ans, qmul(tri_count((Int128)n), (Int128)(a / m))), a %= m;
            if (b >= m) ans = qadd(ans, qmul((Int128)n, (Int128)(b / m))), b %= m;
        }
        const U y_max = a * n + b;  // a < m、b < m なので m(n+1) 未満（n は増えず、m は減る）
        if (y_max < m) break;
        n = y_max / m, b = y_max % m;
        std::swap(m, a);
    }
    return ans;
}
}  // namespace math_internal
inline Int128 linear_floor_sum(long long n, long long m, long long a, long long b) {
    if (n < 0 || m < 1) math_internal::fail("linear_floor_sum", "n ≥ 0、m ≥ 1 で使う");
    if (n == 0) return 0;
    if (n < (1LL << 32) && m < (1LL << 32)) {  // よくある大きさ：64 ビットで計算（答えは 2^127 未満なので溢れない）
        const long long ka = floor_div(a, m), kb = floor_div(b, m);
        const unsigned long long un = (unsigned long long)n;
        return (Int128)(un * (un - 1) / 2) * ka + (Int128)n * kb +
               math_internal::floor_sum_nonneg<unsigned long long>(un, (unsigned long long)m, (unsigned long long)mod_floor(a, m), (unsigned long long)mod_floor(b, m));
    }
    using math_internal::qadd;
    using math_internal::qmul;
    Int128 A = a, B = b, ans = 0;
    const Int128 N = n, M = m;
    if (A < 0 || A >= M) {  // a = a' + k m（0 ≤ a' < m）にして k·n(n−1)/2 を外に出す
        const Int128 k = math_internal::floor_div128(A, M);
        ans = qadd(ans, qmul(math_internal::tri_count(N), k)), A -= k * M;
    }
    if (B < 0 || B >= M) {
        const Int128 k = math_internal::floor_div128(B, M);
        ans = qadd(ans, qmul(N, k)), B -= k * M;
    }
    return qadd(ans, math_internal::floor_sum_nonneg<Int128>(N, M, A, B));
}
// a x + b y ≤ c となる x, y ≥ 0 の整数の組の数（a, b ≥ 1。c < 0 なら 0）  O(log)
inline Int128 count_lattice_under_line(long long a, long long b, long long c) {
    if (a < 1 || b < 1) math_internal::fail("count_lattice_under_line", "a, b ≥ 1 で使う");
    if (c < 0) return 0;
    // x = 0 は ⌊c/b⌋ + 1 個。x = 1〜q（q = ⌊c/a⌋）は Σ_{i=0}^{q−1} (⌊(c − a − a i)/b⌋ + 1)（q + 1 は溢れることがあるので分ける）
    const long long q = c / a;
    return linear_floor_sum(q, b, -a, c - a) + q + c / b + 1;
}
// Σ_{i=0}^{n−1} ((a i + b) mod m)（mod は 0〜m−1 の方。a, b は負でも。答えは n(m−1) 以下なので Int128）  O(log m)
inline Int128 linear_mod_sum(long long n, long long m, long long a, long long b) {
    if (n < 0 || m < 1) math_internal::fail("linear_mod_sum", "n ≥ 0、m ≥ 1 で使う");
    if (n == 0) return 0;
    const long long a2 = mod_floor(a, m), b2 = mod_floor(b, m);
    // Σ (a2 i + b2) − m Σ ⌊(a2 i + b2)/m⌋。途中は 2^128 を法として計算してよい（答えは 0 以上 2^126 未満）
    const UInt128 N = (UInt128)n, T = (N % 2 == 0) ? (N / 2) * (N - 1) : N * ((N - 1) / 2);
    const UInt128 F = (UInt128)linear_floor_sum(n, m, a2, b2);
    return (Int128)((UInt128)a2 * T + (UInt128)b2 * N - (UInt128)m * F);
}

// ---- 床関数の和の拡張：Σ⌊⌋・Σ i⌊⌋・Σ⌊⌋² を mint で（mod は奇数の素数）----
template <class T>
struct FloorSumExt {
    T sum, sum_i, sum_sq;  // Σ_{i=0}^{n−1} ⌊(a i + b)/m⌋、Σ i⌊…⌋、Σ ⌊…⌋^2
};
namespace math_internal {
template <class T>
T fs_s1(long long n) {  // 0 + 1 + … + n（n ≥ −1）
    long long x = n, y = n + 1;
    (x % 2 == 0 ? x : y) /= 2;
    return T(x) * T(y);
}
template <class T>
T fs_from128(Int128 v) {  // 0 ≤ v < 2^125 を T に（T は long long から作れればよい）
    return T((long long)(v >> 62)) * T(1LL << 62) + T((long long)(v & ((1LL << 62) - 1)));
}
template <class T>
T fs_s2(long long n) {  // 0^2 + 1^2 + … + n^2
    Int128 x = n, y = (Int128)n + 1, z = 2 * (Int128)n + 1;
    (x % 2 == 0 ? x : y) /= 2;
    (x % 3 == 0 ? x : (y % 3 == 0 ? y : z)) /= 3;
    return fs_from128<T>(x) * fs_from128<T>(y) * fs_from128<T>(z);
}
template <class T>
FloorSumExt<T> fs_ext(long long n, long long a, long long b, long long c) {  // x = 0〜n、a, b ≥ 0、c ≥ 1
    if (n < 0) return {T(0), T(0), T(0)};
    const T N(n), S1 = fs_s1<T>(n), S2 = fs_s2<T>(n);
    if (a == 0) {
        const T q(b / c);
        return {q * (N + T(1)), q * S1, q * q * (N + T(1))};
    }
    if (a >= c || b >= c) {
        const FloorSumExt<T> t = fs_ext<T>(n, a % c, b % c, c);
        const T qa(a / c), qb(b / c);
        FloorSumExt<T> r;
        r.sum = t.sum + qa * S1 + qb * (N + T(1));
        r.sum_i = t.sum_i + qa * S2 + qb * S1;
        r.sum_sq = qa * qa * S2 + qb * qb * (N + T(1)) + T(2) * qa * qb * S1 + t.sum_sq + T(2) * qb * t.sum + T(2) * qa * t.sum_i;
        return r;
    }
    const long long m = (long long)(((Int128)a * n + b) / c);  // a, b < c なので n + 1 未満
    const FloorSumExt<T> t = fs_ext<T>(m - 1, c, c - b - 1, a);
    const T M(m);
    FloorSumExt<T> r;
    r.sum = N * M - t.sum;
    r.sum_i = (M * N * (N + T(1)) - t.sum_sq - t.sum) / T(2);
    r.sum_sq = N * M * (M + T(1)) - T(2) * t.sum_i - T(2) * t.sum - r.sum;
    return r;
}
}  // namespace math_internal
// n ≥ 0、m ≥ 1、a, b は負でもよい。T は mint（mod は奇数の素数）  O(log m)
template <class T>
FloorSumExt<T> linear_floor_sum_ext(long long n, long long m, long long a, long long b) {
    if (n < 0 || m < 1) math_internal::fail("linear_floor_sum_ext", "n ≥ 0、m ≥ 1 で使う");
    if (n == 0) return {T(0), T(0), T(0)};
    const long long qa = floor_div(a, m), qb = floor_div(b, m);  // ⌊(a i + b)/m⌋ = qa i + qb + ⌊(a' i + b')/m⌋
    const FloorSumExt<T> t = math_internal::fs_ext<T>(n - 1, mod_floor(a, m), mod_floor(b, m), m);
    const T A(qa), B(qb), N(n), S1 = math_internal::fs_s1<T>(n - 1), S2 = math_internal::fs_s2<T>(n - 1);
    FloorSumExt<T> r;
    r.sum = t.sum + A * S1 + B * N;
    r.sum_i = t.sum_i + A * S2 + B * S1;
    r.sum_sq = A * A * S2 + B * B * N + T(2) * A * B * S1 + t.sum_sq + T(2) * B * t.sum + T(2) * A * t.sum_i;
    return r;
}

// ---- 万能ユークリッド（universal Euclid）----
// 直線 y = ⌊(a x + b)/c⌋ に沿って、x を 1 増やすたびに R、y が 1 増えるたびに U を掛けた積
//   U^{y(0)} · Π_{x=1}^{n} ( U^{y(x) − y(x−1)} · R )（a, b ≥ 0、c ≥ 1、n ≥ 0）を求める。
//   M は結合的な掛け算 operator* と単位元 M()（既定のコンストラクタ）を持つ型（行列・自作の struct など。可換でなくてよい）
//   Σ_x A^x B^{y(x)} のような「床関数を含む和」を、M に必要な量を持たせて求める  O(log max(a, c) × M の掛け算 × log)
namespace math_internal {
template <class M>
M ue_pow(M x, long long e) {  // 最後の 2 乗はしない（使わない値で溢れないように）
    M r;
    while (e > 0) {
        if (e & 1) r = r * x;
        e >>= 1;
        if (e > 0) x = x * x;
    }
    return r;
}
template <class M>
M ue_solve(long long p, long long q, long long r, long long l, const M& U, const M& R) {  // 0 ≤ r < q。y = ⌊(p x + r)/q⌋、x = 1〜l
    if (l == 0) return M();
    if (p >= q) return ue_solve<M>(p % q, q, r, l, U, ue_pow(U, p / q) * R);
    const long long m = (long long)(((Int128)l * p + r) / q);
    if (m == 0) return ue_pow(R, l);
    const long long cnt = l - (long long)(((Int128)q * m - r - 1) / p);
    return ue_pow(R, (q - r - 1) / p) * U * ue_solve<M>(q, p, (q - r - 1) % p, m - 1, R, U) * ue_pow(R, cnt);
}
}  // namespace math_internal
template <class M>
M universal_euclid(long long n, long long a, long long b, long long c, const M& U, const M& R) {
    if (n < 0 || a < 0 || b < 0 || c < 1) math_internal::fail("universal_euclid", "n, a, b ≥ 0、c ≥ 1 で使う");
    return math_internal::ue_pow(U, b / c) * math_internal::ue_solve<M>(a, c, b % c, n, U, R);
}
#endif  // GG_QUOTIENT
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace quot_test {

#define QUOT_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[quotient] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static const long long LMAX = std::numeric_limits<long long>::max();

static void check_blocks(long long N, int& cases) {  // 区間がすき間なく並び、各区間で値が一定で、区間が最大
    ++cases;
    const std::vector<QuotientBlock> f = quotient_blocks(N), c = quotient_blocks_ceil(N);
    long long x = 1;
    for (const auto& b : f) {
        QUOT_CHECK(b.l == x && b.l < b.r && N / b.l == b.q && N / (b.r - 1) == b.q && (b.r > N || N / b.r != b.q));
        x = b.r;
    }
    QUOT_CHECK(x == N + 1 || (N == 0 && f.empty()));
    QUOT_CHECK((long long)f.size() <= 2 * isqrt(N) + 1);
    x = 1;
    for (const auto& b : c) {
        auto cd = [&](long long i) { return (N + i - 1) / i; };
        QUOT_CHECK(b.l == x && b.l < b.r && cd(b.l) == b.q && cd(b.r - 1) == b.q && (b.r > N || cd(b.r) != b.q));
        x = b.r;
    }
    QUOT_CHECK(x == N + 1 || (N == 0 && c.empty()));
}

static void test_blocks(std::mt19937_64& rng, int& cases) {
    for (long long N = 0; N <= 2000; ++N) {
        check_blocks(N, cases);
        // 1 つずつ数える
        Int128 s1 = 0, s2 = 0, s3 = 0;
        for (long long i = 1; i <= N; ++i) s1 += N / i, s2 += (Int128)i * (N / i), s3 += N % i;
        QUOT_CHECK(sum_floor_div(N) == s1 && sum_i_floor_div(N) == s2 && sum_mod_i(N) == s3);
        long long cnt = 0;
        for_each_quotient(N, [&](long long q, long long l, long long r) {
            for (long long i = l; i < r; ++i) QUOT_CHECK(N / i == q);
            cnt += r - l;
        });
        QUOT_CHECK(cnt == N);
    }
    for (int it = 0; it < 100 * TEST_SCALE; ++it) check_blocks(rnd(rng, 0, 100000000LL), cases);
    check_blocks(1000000000000LL, cases), check_blocks(999999999999LL, cases);
    {  // N = 10^14（区間は 2×10^7 個）：for_each_quotient で区間がつながっていることと q の確認（一部）
        const long long N = 100000000000000LL;
        long long x = 1, k = 0;
        for_each_quotient(N, [&](long long q, long long l, long long r) {
            if (++k % 1000003 == 0 || l < 1000) QUOT_CHECK(l == x && N / l == q && N / (r - 1) == q && (r > N || N / r != q));
            x = r;
        });
        QUOT_CHECK(x == N + 1 && k <= 2 * isqrt(N));
    }
    // 2 つの商
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long N = rnd(rng, 0, it % 2 ? 300 : 1000000), M = rnd(rng, 0, it % 3 ? 300 : 1000000);
        long long x = 1;
        for (const auto& b : quotient_blocks2(N, M)) {
            QUOT_CHECK(b.l == x && b.l < b.r && N / b.l == b.q1 && M / b.l == b.q2 && N / (b.r - 1) == b.q1 && M / (b.r - 1) == b.q2);
            QUOT_CHECK(b.r > std::min(N, M) || N / b.r != b.q1 || M / b.r != b.q2);
            x = b.r;
        }
        QUOT_CHECK(x == std::min(N, M) + 1 || (std::min(N, M) == 0 && x == 1));
    }
    // 大きい N：1 つ前との差は約数の個数・約数の和
    for (int it = 0; it < 40 * TEST_SCALE; ++it, ++cases) {
        const long long N = rnd(rng, 2, it % 4 ? 1000000000LL : 100000000000LL);
        long long d = 0;
        Int128 sig = 0;
        for (long long i = 1; i * i <= N; ++i)
            if (N % i == 0) d += (i * i == N ? 1 : 2), sig += i + (i * i == N ? 0 : N / i);
        QUOT_CHECK(sum_floor_div(N) - sum_floor_div(N - 1) == d && sum_i_floor_div(N) - sum_i_floor_div(N - 1) == sig);
    }
}

static void test_floor_sum(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 40);
        const long long m = it % 3 == 0 ? rnd(rng, 1, 10) : rnd(rng, 1, it % 2 ? 1000000000LL : LMAX);
        const long long a = it % 4 == 0 ? rnd(rng, -20, 20) : rnd(rng, -1000000000000LL, 1000000000000LL);
        const long long b = it % 5 == 0 ? rnd(rng, -20, 20) : rnd_ll(rng);
        Int128 e = 0;
        for (long long i = 0; i < n; ++i) e += math_internal::floor_div128((Int128)a * i + b, m);
        QUOT_CHECK(linear_floor_sum(n, m, a, b) == e);
    }
    for (int it = 0; it < 50000 * TEST_SCALE; ++it, ++cases) {  // 大きい n：1 つ前との差
        const long long n = rnd(rng, 1, 1000000000000000000LL), m = rnd(rng, 1, 1000000000000000000LL);
        const long long a = rnd(rng, -1000000000LL, 1000000000LL), b = rnd(rng, -1000000000000000000LL, 1000000000000000000LL);
        QUOT_CHECK(linear_floor_sum(n, m, a, b) - linear_floor_sum(n - 1, m, a, b) == math_internal::floor_div128((Int128)a * (n - 1) + b, m));
    }
    QUOT_CHECK(linear_floor_sum(1000000000, 1000000000, 999999999, 999999999) == 499999999500000000LL);
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long a = rnd(rng, 1, 30), b = rnd(rng, 1, 30), c = rnd(rng, -5, 300);
        long long e = 0;
        for (long long x = 0; a * x <= c; ++x)
            for (long long y = 0; a * x + b * y <= c; ++y) ++e;
        QUOT_CHECK(count_lattice_under_line(a, b, c) == e);
    }
    QUOT_CHECK(count_lattice_under_line(1, 1, 1000000000000000000LL) == (Int128)1000000000000000001LL * 1000000000000000002LL / 2);
    QUOT_CHECK(count_lattice_under_line(1, 1, LMAX) == ((Int128)1 << 125) + ((Int128)1 << 62));  // c = LLONG_MAX でも（(c+1)(c+2)/2）
    for (int it = 0; it < 3000; ++it, ++cases) {  // c が LLONG_MAX 近く：x の数が少ないように a を大きくして 1 列ずつ足す。a と b を入れかえても同じ
        const long long c = rnd(rng, LMAX - 1000000, LMAX), a = rnd(rng, c / 1000, c), b = it % 2 ? rnd(rng, 1, 1000) : rnd(rng, 1, LMAX);
        Int128 e = 0;
        for (long long x = 0; x <= c / a; ++x) e += (c - a * x) / b + 1;
        QUOT_CHECK(count_lattice_under_line(a, b, c) == e && count_lattice_under_line(b, a, c) == e);
    }
}

struct ModP {  // テスト用の mod 型（998244353。割り算も）
    static constexpr long long MD = 998244353;
    long long v;
    ModP(long long x = 0) : v(x % MD) {
        if (v < 0) v += MD;
    }
    friend ModP operator+(ModP a, ModP b) { return ModP(a.v + b.v); }
    friend ModP operator-(ModP a, ModP b) { return ModP(a.v - b.v); }
    friend ModP operator*(ModP a, ModP b) { return ModP(a.v * b.v % MD); }
    friend ModP operator/(ModP a, ModP b) { return a * ModP(mod_pow(b.v, MD - 2, MD)); }
    friend bool operator==(ModP a, ModP b) { return a.v == b.v; }
};
static ModP m128(Int128 x) { return ModP((long long)(((x % ModP::MD) + ModP::MD) % ModP::MD)); }
struct UeNode {  // 万能ユークリッドで Σ y・Σ x y・Σ y^2 を求めるための値（x = R の数、y = U の数）
    Int128 x = 0, y = 0, s = 0, sx = 0, sy2 = 0;
    friend UeNode operator*(const UeNode& a, const UeNode& b) {  // a のあとに b（b の各点の x, y は a の分だけずれる）
        UeNode r;
        r.x = a.x + b.x, r.y = a.y + b.y;
        r.s = a.s + b.s + a.y * b.x;
        r.sx = a.sx + b.sx + a.x * b.s + a.y * (b.x * (b.x + 1) / 2) + a.x * a.y * b.x;
        r.sy2 = a.sy2 + b.sy2 + 2 * a.y * b.s + a.y * a.y * b.x;
        return r;
    }
};
struct StrMonoid {  // 可換でない積の確かめ：文字列をつなげる
    std::string s;
    friend StrMonoid operator*(const StrMonoid& a, const StrMonoid& b) { return StrMonoid{a.s + b.s}; }
};

static void test_extensions(std::mt19937_64& rng, int& cases) {
    // for_each_quotient を途中で止める
    for (int it = 0; it < 2000; ++it, ++cases) {
        const long long N = rnd(rng, 0, 100000), stop = rnd(rng, 0, 50);
        long long calls = 0, all = 0;
        for_each_quotient(N, [&](long long, long long, long long) { ++all; });
        for_each_quotient(N, [&](long long, long long, long long) { return ++calls < stop; });
        QUOT_CHECK(calls == std::min(all, std::max(stop, 1LL)));
    }
    // Σ (a i + b) mod m
    for (int it = 0; it < 50000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 60), m = it % 3 == 0 ? rnd(rng, 1, 10) : rnd(rng, 1, it % 2 ? 1000000000LL : LMAX);
        const long long a = it % 4 == 0 ? rnd(rng, -20, 20) : rnd_ll(rng), b = rnd_ll(rng);
        Int128 e = 0;
        for (long long i = 0; i < n; ++i) {
            const Int128 v = ((Int128)a * i + b) % m;
            e += v < 0 ? v + m : v;
        }
        QUOT_CHECK(linear_mod_sum(n, m, a, b) == e);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい n：m ごとの周期で
        const long long m = rnd(rng, 1, 300), n = rnd(rng, 0, LMAX), a = rnd_ll(rng), b = rnd_ll(rng);
        Int128 per = 0, rest = 0;
        for (long long i = 0; i < m; ++i) {
            const long long v = mod_floor((long long)(((Int128)mod_floor(a, m) * i + mod_floor(b, m)) % m), m);
            per += v;
            if (i < n % m) rest += v;
        }
        QUOT_CHECK(linear_mod_sum(n, m, a, b) == per * (n / m) + rest);
    }
    QUOT_CHECK(linear_mod_sum(LMAX, LMAX, LMAX - 1, LMAX - 1) > 0 && linear_mod_sum(0, 5, 1, 1) == 0);
    // Σ⌊⌋・Σ i⌊⌋・Σ⌊⌋^2（mod p）
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 50), m = it % 3 == 0 ? rnd(rng, 1, 10) : rnd(rng, 1, it % 2 ? 1000000000LL : LMAX);
        const long long a = it % 4 == 0 ? rnd(rng, -20, 20) : rnd_ll(rng), b = it % 5 == 0 ? rnd(rng, -20, 20) : rnd_ll(rng);
        Int128 s = 0, si = 0;
        ModP sq = 0;
        for (long long i = 0; i < n; ++i) {
            const Int128 f = math_internal::floor_div128((Int128)a * i + b, m);
            s += f, si += f * i, sq = sq + m128(f) * m128(f);
        }
        const FloorSumExt<ModP> r = linear_floor_sum_ext<ModP>(n, m, a, b);
        QUOT_CHECK(r.sum == m128(s) && r.sum_i == m128(si) && r.sum_sq == sq);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい n：1 つ前との差
        const long long n = rnd(rng, 1, LMAX), m = rnd(rng, 1, LMAX), a = rnd_ll(rng), b = rnd_ll(rng);
        const FloorSumExt<ModP> r = linear_floor_sum_ext<ModP>(n, m, a, b), q = linear_floor_sum_ext<ModP>(n - 1, m, a, b);
        const ModP f = m128(math_internal::floor_div128((Int128)a * (n - 1) + b, m));
        QUOT_CHECK(r.sum - q.sum == f && r.sum_i - q.sum_i == f * ModP(n - 1) && r.sum_sq - q.sum_sq == f * f);
        if (m >= 1000000000LL && a >= -1000000000LL && a <= 1000000000LL) QUOT_CHECK(r.sum == m128(linear_floor_sum(n, m, a, b)));  // Int128 に収まる範囲で
    }
    // 万能ユークリッド：文字列（順番）と Σ の値
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 30), a = rnd(rng, 0, 20), b = rnd(rng, 0, 40), c = rnd(rng, 1, 20);
        std::string e(std::string((size_t)(b / c), 'U'));
        for (long long x = 1; x <= n; ++x) e += std::string((size_t)((a * x + b) / c - (a * (x - 1) + b) / c), 'U') + "R";
        QUOT_CHECK(universal_euclid(n, a, b, c, StrMonoid{"U"}, StrMonoid{"R"}).s == e);
        UeNode U, R;
        U.y = 1, R.x = 1;
        const UeNode got = universal_euclid(n, a, b, c, U, R);
        Int128 s = 0, sx = 0, sy2 = 0;
        for (long long x = 1; x <= n; ++x) {
            const Int128 y = (a * x + b) / c;
            s += y, sx += x * y, sy2 += y * y;
        }
        QUOT_CHECK(got.x == n && got.y == (a * n + b) / c && got.s == s && got.sx == sx && got.sy2 == sy2);
    }
    {  // 溢れを調べる掛け算のモノイド：使わない最後の 2 乗をしないこと（2^40 は収まるが 2^64 は収まらない）
        struct ProdM {
            long long x = 1;
            ProdM operator*(const ProdM& o) const {
                ProdM r;
                if (__builtin_mul_overflow(x, o.x, &r.x)) math_internal::fail("ProdM", "溢れた");
                return r;
            }
        };
        ProdM U, R;
        U.x = 3, R.x = 2;
        QUOT_CHECK(universal_euclid(40, 0, 0, 1, U, R).x == (1LL << 40) && universal_euclid(10, 1, 0, 1, U, R).x == 60466176LL);  // 6^10 = 3^10 2^10
        cases += 2;
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい値：Σ y は linear_floor_sum と同じ
        const long long n = rnd(rng, 0, 1000000000LL), a = rnd(rng, 0, 1000000000LL), b = rnd(rng, 0, 1000000000LL), c = rnd(rng, 1, 1000000000LL);
        UeNode U, R;
        U.y = 1, R.x = 1;
        QUOT_CHECK(universal_euclid(n, a, b, c, U, R).s == linear_floor_sum(n + 1, c, a, b) - b / c);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_blocks(rng, cases);
    test_floor_sum(rng, cases);
    test_extensions(rng, cases);
    std::cout << "[quotient] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[quotient] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    const Int128 s1 = sum_floor_div(100000000000000LL);
    const long long t_sum = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Int128 s2 = sum_i_floor_div(100000000000000LL);
    const long long t_isum = ms(t0);
    std::vector<std::array<long long, 4>> qs(1000000);
    for (auto& q : qs) q = {rnd(rng, 1, 1000000000), rnd(rng, 1, 1000000000), rnd(rng, 0, 1000000000), rnd(rng, 0, 1000000000)};
    t0 = std::chrono::steady_clock::now();
    Int128 s3 = 0;
    for (const auto& q : qs) s3 += linear_floor_sum(q[0], q[1], q[2], q[3]);
    const long long t_fs = ms(t0);
    std::cout << "[quotient] benchmark: sum_floor_div(10^14) " << t_sum << " ms, sum_i_floor_div(10^14) " << t_isum << " ms, linear_floor_sum × 10^6 " << t_fs
              << " ms  (check " << (long long)((s1 + s2 + s3) % 7) << ")\n";
}

}  // namespace quot_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace quot_demo {

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
#define SHOW(expr, comment) quot_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    quot_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("商の列挙・床関数の和（quotient） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  for_each_quotient(N, f) / quotient_blocks(N) / quotient_blocks_ceil(N) / quotient_blocks2(N, M)   ⌊N/i⌋ が同じ i をまとめる\n"
              << "  sum_floor_div(N) / sum_i_floor_div(N) / sum_mod_i(N)        Σ⌊N/i⌋ / Σ i⌊N/i⌋ / Σ (N mod i)\n"
              << "  linear_floor_sum(n, m, a, b) / count_lattice_under_line(a, b, c)   Σ⌊(a i + b)/m⌋ / 直線の下の格子点\n";

    api("【1】 for_each_quotient(N, f) / quotient_blocks(N)", "なし（f(q, l, r) を呼ぶ）/ vector<QuotientBlock{q, l, r}>",
        "i = 1〜N を ⌊N/i⌋ = q が同じ区間 [l, r) に分ける（区間は 2√N 個以下）", "O(√N)");
    RUN(std::vector<std::string> v; for_each_quotient(10, [&](long long q, long long l, long long r) { v.push_back("q=" + std::to_string(q) + ":[" + std::to_string(l) + "," + std::to_string(r) + ")"); }));
    SHOW(v, "10/1 = 10、10/2 = 5、10/3 = 3、10/4 = 10/5 = 2、10/6〜10/10 = 1");
    SHOW(quotient_blocks(1000000000000LL).size(), "N = 10^12 でも 2×10^6 個くらい");
    SHOW(quotient_blocks_ceil(10).size(), "⌈10/i⌉：10, 5, 4, 3, 2（i = 5〜9）, 1（i = 10）");
    SHOW(quotient_blocks2(10, 7).size(), "⌊10/i⌋ と ⌊7/i⌋ が両方同じ区間（i ≤ 7）");
    {
        // ABC172 D：Σ_{K ≤ N} K·d(K) = Σ_i i · T(⌊N/i⌋)（T(m) = m(m+1)/2）
        auto solve = [](long long N) {
            Int128 s = 0;
            for_each_quotient(N, [&](long long q, long long l, long long r) { s += sum_range(l, r) * (Int128)(q * (q + 1) / 2); });
            return s;
        };
        note("ABC172 D（Sum of Divisors）: 区間 [l, r) の i の和（sum_range）× T(q)。O(√N)");
        std::cout << "    N = 4 → " << solve(4) << "、N = 100 → " << solve(100) << "、N = 10^7 → " << solve(10000000) << "\n";
    }
    {
        // ABC044 D：b ≤ √n は直接、b > √n は 2 桁で ⌊n/b⌋ = q が一定の区間ごとに b = (n − s)/q + 1
        auto solve = [](long long n, long long s) {
            if (s == n) return n + 1;  // b ≤ n なら 2 桁以上で桁和 < n なので、b = n + 1（1 桁）
            long long best = -1;
            for (long long b = 2; b * b <= n && best < 0; ++b) {
                long long t = 0;
                for (long long x = n; x; x /= b) t += x % b;
                if (t == s) best = b;
            }
            if (best >= 0) return best;
            for_each_quotient(n, [&](long long q, long long l, long long r) {
                if (best >= 0 || n < s || (n - s) % q != 0) return;
                const long long b = (n - s) / q + 1;  // 桁和 = q + (n − q b)
                if (l <= b && b < r && b > isqrt(n) && b >= 2) best = b;  // b·b は溢れることがあるので isqrt
            });
            return best;
        };
        note("ABC044 D（桁和）: b ≤ √n は直接。b > √n は ⌊n/b⌋ = q が同じ区間で n − q(b − 1) = s を解く");
        std::cout << "    (87654, 30) → " << solve(87654, 30) << "、(87654, 138) → " << solve(87654, 138) << "、(87654, 45678) → " << solve(87654, 45678)
                  << "、(31415926535, 1) → " << solve(31415926535LL, 1) << "、(1, 31415926535) → " << solve(1, 31415926535LL) << "\n";
    }

    api("【2】 sum_floor_div(N) / sum_i_floor_div(N) / sum_mod_i(N)", "Int128",
        "Σ⌊N/i⌋（約数の個数の和）/ Σ i⌊N/i⌋（約数の和の和）/ Σ (N mod i)", "O(√N)");
    SHOW(sum_floor_div(3), "ABC230 E: 3 + 1 + 1");
    SHOW(sum_floor_div(10000000000LL), "ABC230 E: N = 10^10");
    SHOW(sum_i_floor_div(4), "σ(1) + σ(2) + σ(3) + σ(4) = 1 + 3 + 4 + 7");
    SHOW(sum_mod_i(10), "0 + 0 + 1 + 2 + 0 + 4 + 3 + 2 + 1 + 0");

    api("【3】 linear_floor_sum(n, m, a, b)", "Int128", "Σ_{i=0}^{n−1} ⌊(a i + b) / m⌋（a, b は負でもよい。ACL の floor_sum と同じで範囲が広い）", "O(log m)");
    SHOW(linear_floor_sum(4, 10, 6, 3), "ACL Practice C 例 1");
    SHOW(linear_floor_sum(6, 5, 4, 3), "例 2");
    SHOW(linear_floor_sum(1, 1, 0, 0), "例 3");
    SHOW(linear_floor_sum(31415, 92653, 58979, 32384), "例 4");
    SHOW(linear_floor_sum(1000000000, 1000000000, 999999999, 999999999), "例 5");
    SHOW(linear_floor_sum(5, 3, -2, 1), "負の係数も：⌊1/3⌋ + ⌊−1/3⌋ + ⌊−3/3⌋ + ⌊−5/3⌋ + ⌊−7/3⌋ = 0 − 1 − 1 − 2 − 3");
    note("引数の順番は ACL の floor_sum(n, m, a, b) と同じ（個数 n、割る数 m、傾き a、切片 b）");
    SHOW(linear_mod_sum(5, 7, 3, 1), "Σ (3i + 1) mod 7 = 1 + 4 + 0 + 3 + 6");
    note("(a i + b) mod m の和は「Σ (a i + b) − m Σ ⌊(a i + b)/m⌋」なので linear_mod_sum で O(log m)");

    api("【4】 count_lattice_under_line(a, b, c)", "Int128", "a x + b y ≤ c となる x, y ≥ 0 の組の数（三角形の中の格子点）", "O(log)");
    SHOW(count_lattice_under_line(2, 3, 6), "(0,0) (0,1) (0,2) (1,0) (1,1) (2,0) (3,0)");
    SHOW(count_lattice_under_line(1, 1, 1000000000000LL), "(10^12 + 1)(10^12 + 2) / 2");

    api("【5】 linear_floor_sum_ext<mint>(n, m, a, b)", "FloorSumExt<mint>{sum, sum_i, sum_sq}",
        "Σ⌊(a i + b)/m⌋・Σ i⌊…⌋・Σ⌊…⌋^2（i = 0〜n−1）を mod で（mod は奇数の素数）", "O(log m)");
    {
        const auto r = linear_floor_sum_ext<quot_test::ModP>(4, 3, 2, 1);  // ⌊1/3⌋, ⌊3/3⌋, ⌊5/3⌋, ⌊7/3⌋ = 0, 1, 1, 2
        std::cout << "  linear_floor_sum_ext<mint>(4, 3, 2, 1) -> {" << r.sum.v << ", " << r.sum_i.v << ", " << r.sum_sq.v
                  << "}   // 値 0, 1, 1, 2：和 4、Σ i·値 = 0 + 1 + 2 + 6 = 9、Σ 値^2 = 6\n";
    }

    api("【6】 universal_euclid(n, a, b, c, U, R)", "M（自作の型）",
        "直線 y = ⌊(a x + b)/c⌋ に沿って、y が 1 増えるたびに U、x が 1 増えるたびに R を掛けた積（x = 1〜n。M は結合的な * と単位元 M()）",
        "O(log × M の掛け算 × log)");
    note("Σ_x A^x B^{⌊(ax+b)/c⌋}（行列）のような床関数を含む和に。M に「U の数・R の数・欲しい和」を持たせて * を書く");
    {
        const auto s = universal_euclid(4, 2, 1, 3, quot_test::StrMonoid{"U"}, quot_test::StrMonoid{"R"}).s;
        std::cout << "  universal_euclid(4, 2, 1, 3, \"U\", \"R\") -> " << s << "   // y(x) = ⌊(2x+1)/3⌋ = 1, 1, 2, 3（x = 1〜4）\n";
    }
}

}  // namespace quot_demo

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
    quot_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    quot_test::benchmark();    // 速度計測
    quot_demo::demo();         // 使い方デモ
    return 0;
}
