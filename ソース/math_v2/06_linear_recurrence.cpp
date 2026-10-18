// =============================================================================
//  線形漸化式（項の多い漸化式の第 K 項・最初の n 項、最初の数項から漸化式を求める Berlekamp–Massey）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 06 数列・和
//  依存       : 00_int_common, 00_math_base, 00_modint, 07_convolution（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  漸化式 a_i = c_1 a_{i−1} + c_2 a_{i−2} + … + c_d a_{i−d}
//    a = {a_0, …, a_{d−1}}（最初の d 項）、c = {c_1, …, c_d}（c[0] が 1 つ前の項の係数）
//  使い方
//      using mint = ModInt998;
//      mint f = linear_recurrence_kth<mint>({0, 1}, {1, 1}, K);   // フィボナッチ数 F_K（K ≤ 10^18）
//      mint g = guess_kth(first_terms, K);                        // DP で出した最初の数十項から第 K 項を推測
//
//  機能一覧                                                                    計算量
//    linear_recurrence_kth(a, c, K)   第 K 項（Bostan–Mori。mod は何でも。d = 5×10^4 で約 0.5 秒。
//                                     10^9+7 など NTT が使えない mod は 3 つの素数の畳み込みなので約 5 倍遅い）  O(d log d log K)
//    linear_recurrence_terms(a, c, n) 最初の n 項                                                 O(n log n)
//    berlekamp_massey(s)              列 s をみたす最短の漸化式の c（mod は素数。次数 d には 2d 項以上いる）  O(N^2)
//    guess_kth(s, K)                  berlekamp_massey + linear_recurrence_kth
//  （項が少ない（d ≤ 10 程度）なら 05_matrix の行列累乗でもよい。P(x)/Q(x) の第 K 係数は 07_fps の fps_coeff_of_fraction）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    linrec_test::random_test()  … 漸化式どおりに計算した値との突き合わせ
//    linrec_test::benchmark()    … 速度計測
//    linrec_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#include <queue>
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
#ifndef GG_CONVOLUTION
#define GG_CONVOLUTION
// ---- 畳み込み（数学第二世代）  mint の列（NTT・任意の mod）、整数のまま（long long・Int128）、多くの多項式の積、相関、巡回、オンライン ----
namespace math_internal {
template <class mint>
struct NttTable {  // mod ごとの 1 の冪根の表（ModInt は 1 回、DynModInt は mod が変わったら作り直す）
    int mod = 0, rank2 = 0;
    std::vector<mint> root, iroot, rate2, irate2, rate3, irate3;
};
inline unsigned primitive_root_u32(unsigned p) {  // 素数 p の原始根（いちばん小さいもの）
    if (p == 2) return 1;
    std::vector<unsigned> qs;
    unsigned x = p - 1;
    for (unsigned q = 2; (unsigned long long)q * q <= x; ++q)
        if (x % q == 0) {
            qs.push_back(q);
            while (x % q == 0) x /= q;
        }
    if (x > 1) qs.push_back(x);
    for (unsigned g = 2;; ++g) {
        bool ok = true;
        for (unsigned q : qs) {
            unsigned long long r = 1, b = g, e = (p - 1) / q;
            for (; e; e >>= 1, b = b * b % p)
                if (e & 1) r = r * b % p;
            if (r == 1) { ok = false; break; }
        }
        if (ok) return g;
    }
}
template <class mint>
const NttTable<mint>& ntt_table() {
    static NttTable<mint> t;
    if (t.mod == mint::mod()) return t;
    const unsigned p = (unsigned)mint::mod();
    if (!is_prime_u32(p)) fail("ntt", "mod が素数でない（convolve なら自動で 3 つの素数を使う）");
    t = NttTable<mint>();
    t.mod = (int)p;
    t.rank2 = __builtin_ctz(p - 1);
    const int r = t.rank2;
    t.root.assign(r + 1, mint(1)), t.iroot.assign(r + 1, mint(1));
    t.rate2.assign(std::max(0, r - 1), mint(1)), t.irate2.assign(std::max(0, r - 1), mint(1));
    t.rate3.assign(std::max(0, r - 2), mint(1)), t.irate3.assign(std::max(0, r - 2), mint(1));
    const mint g = mint((long long)primitive_root_u32(p));
    t.root[r] = g.pow((long long)((p - 1) >> r));
    t.iroot[r] = t.root[r].inv();
    for (int i = r - 1; i >= 0; --i) t.root[i] = t.root[i + 1] * t.root[i + 1], t.iroot[i] = t.iroot[i + 1] * t.iroot[i + 1];
    mint prod = 1, iprod = 1;
    for (int i = 0; i + 2 <= r; ++i) {
        t.rate2[i] = t.root[i + 2] * prod, t.irate2[i] = t.iroot[i + 2] * iprod;
        prod *= t.iroot[i + 2], iprod *= t.root[i + 2];
    }
    prod = 1, iprod = 1;
    for (int i = 0; i + 3 <= r; ++i) {
        t.rate3[i] = t.root[i + 3] * prod, t.irate3[i] = t.iroot[i + 3] * iprod;
        prod *= t.iroot[i + 3], iprod *= t.root[i + 3];
    }
    return t;
}
template <class mint>
bool ntt_usable(size_t len) {  // mod が素数で、長さ len（2 の冪）の NTT ができるか
    const unsigned p = (unsigned)mint::mod();
    if (p < 3 || !is_prime_u32(p)) return false;
    return (size_t)1 << __builtin_ctz(p - 1) >= len;
}
inline int ceil_log2(size_t n) {
    int h = 0;
    while (((size_t)1 << h) < n) ++h;
    return h;
}
}  // namespace math_internal

// 長さ 2^h の NTT（結果はビットを逆にした順。intt と組で、点ごとの積をとるのに使う。mod は素数で 2^h | mod − 1）
template <class mint>
void ntt(std::vector<mint>& a) {
    const int n = (int)a.size(), h = math_internal::ceil_log2(a.size());
    if (n == 0) return;
    if ((1 << h) != n) math_internal::fail("ntt", "長さは 2 の冪");
    const auto& t = math_internal::ntt_table<mint>();
    if (h > t.rank2) math_internal::fail("ntt", "この mod では長さが大きすぎる");
    int len = 0;
    while (len < h) {
        if (h - len == 1) {
            const int p = 1 << (h - len - 1);
            mint rot = 1;
            for (int s = 0; s < (1 << len); ++s) {
                const int off = s << (h - len);
                for (int i = 0; i < p; ++i) {
                    const mint l = a[i + off], r = a[i + off + p] * rot;
                    a[i + off] = l + r, a[i + off + p] = l - r;
                }
                if (s + 1 != (1 << len)) rot *= t.rate2[__builtin_ctz(~(unsigned)s)];
            }
            ++len;
        } else {
            const int p = 1 << (h - len - 2);
            mint rot = 1;
            const mint imag = t.root[2];
            for (int s = 0; s < (1 << len); ++s) {
                const mint rot2 = rot * rot, rot3 = rot2 * rot;
                const int off = s << (h - len);
                for (int i = 0; i < p; ++i) {
                    const mint a0 = a[i + off], a1 = a[i + off + p] * rot, a2 = a[i + off + 2 * p] * rot2, a3 = a[i + off + 3 * p] * rot3;
                    const mint x = (a1 - a3) * imag;
                    a[i + off] = a0 + a2 + a1 + a3;
                    a[i + off + p] = a0 + a2 - a1 - a3;
                    a[i + off + 2 * p] = a0 - a2 + x;
                    a[i + off + 3 * p] = a0 - a2 - x;
                }
                if (s + 1 != (1 << len)) rot *= t.rate3[__builtin_ctz(~(unsigned)s)];
            }
            len += 2;
        }
    }
}
template <class mint>
void intt(std::vector<mint>& a, bool divide = true) {  // ntt の逆（divide なら 1/n 倍まで）
    const int n = (int)a.size(), h = math_internal::ceil_log2(a.size());
    if (n == 0) return;
    if ((1 << h) != n) math_internal::fail("intt", "長さは 2 の冪");
    const auto& t = math_internal::ntt_table<mint>();
    if (h > t.rank2) math_internal::fail("intt", "この mod では長さが大きすぎる");
    int len = h;
    while (len) {
        if (len == 1) {
            const int p = 1 << (h - len);
            mint irot = 1;
            for (int s = 0; s < (1 << (len - 1)); ++s) {
                const int off = s << (h - len + 1);
                for (int i = 0; i < p; ++i) {
                    const mint l = a[i + off], r = a[i + off + p];
                    a[i + off] = l + r, a[i + off + p] = (l - r) * irot;
                }
                if (s + 1 != (1 << (len - 1))) irot *= t.irate2[__builtin_ctz(~(unsigned)s)];
            }
            --len;
        } else {
            const int p = 1 << (h - len);
            mint irot = 1;
            const mint iimag = t.iroot[2];
            for (int s = 0; s < (1 << (len - 2)); ++s) {
                const mint irot2 = irot * irot, irot3 = irot2 * irot;
                const int off = s << (h - len + 2);
                for (int i = 0; i < p; ++i) {
                    const mint a0 = a[i + off], a1 = a[i + off + p], a2 = a[i + off + 2 * p], a3 = a[i + off + 3 * p];
                    const mint x = (a2 - a3) * iimag;
                    a[i + off] = a0 + a1 + a2 + a3;
                    a[i + off + p] = (a0 - a1 + x) * irot;
                    a[i + off + 2 * p] = (a0 + a1 - a2 - a3) * irot2;
                    a[i + off + 3 * p] = (a0 - a1 - x) * irot3;
                }
                if (s + 1 != (1 << (len - 2))) irot *= t.irate3[__builtin_ctz(~(unsigned)s)];
            }
            len -= 2;
        }
    }
    if (divide) {
        const mint inv_n = mint(n).inv();
        for (mint& x : a) x *= inv_n;
    }
}

namespace math_internal {
template <class mint>
std::vector<mint> conv_naive(const std::vector<mint>& a, const std::vector<mint>& b) {
    std::vector<mint> c(a.size() + b.size() - 1);
    if (a.size() < b.size()) {
        for (size_t j = 0; j < b.size(); ++j)
            for (size_t i = 0; i < a.size(); ++i) c[i + j] += a[i] * b[j];
    } else {
        for (size_t i = 0; i < a.size(); ++i)
            for (size_t j = 0; j < b.size(); ++j) c[i + j] += a[i] * b[j];
    }
    return c;
}
template <class mint>
std::vector<mint> conv_ntt(std::vector<mint> a, std::vector<mint> b) {  // mod が NTT に使える素数のとき
    const size_t z = a.size() + b.size() - 1, n = (size_t)1 << ceil_log2(z);
    const bool same = a == b;
    const mint inv_n = mint((long long)n).inv();  // intt の 1/n は点ごとの積に混ぜる
    a.resize(n), ntt(a);
    if (same) {
        for (size_t i = 0; i < n; ++i) a[i] *= a[i] * inv_n;
    } else {
        b.resize(n), ntt(b);
        for (size_t i = 0; i < n; ++i) a[i] *= b[i] * inv_n;
    }
    intt(a, false);
    a.resize(z);
    return a;
}
// 3〜5 個の NTT 素数（2^23 以上の 2 の冪で割れる）での畳み込みの余り
constexpr unsigned NTT_P0 = 754974721, NTT_P1 = 167772161, NTT_P2 = 469762049, NTT_P3 = 998244353, NTT_P4 = 2013265921;
template <unsigned P>
std::vector<unsigned> conv_residue(const std::vector<long long>& a, const std::vector<long long>& b) {
    using M = ModInt<P>;
    std::vector<M> x(a.size()), y(b.size());
    for (size_t i = 0; i < a.size(); ++i) x[i] = M(a[i]);
    for (size_t i = 0; i < b.size(); ++i) y[i] = M(b[i]);
    const std::vector<M> z = (std::min(a.size(), b.size()) <= 60) ? conv_naive(x, y) : conv_ntt(std::move(x), std::move(y));
    std::vector<unsigned> r(z.size());
    for (size_t i = 0; i < z.size(); ++i) r[i] = (unsigned)z[i].val();
    return r;
}
// 素数 ps[0..k) での余り rs から、−M/2 < x < M/2（M = Π ps）の x を 2^128 で割った余りとして戻す（つり合った混合基数のガーナー）
struct GarnerBalanced {
    int k;
    long long p[5], inv[5], pm[5][5];  // inv[j] = (p_0 … p_{j−1})^{−1} mod p_j、pm[i][j] = p_0 … p_{i−1} mod p_j
    UInt128 base[5];                   // p_0 … p_{j−1}（2^128 で割った余り）
    GarnerBalanced(const unsigned* ps, int k_) : k(k_) {
        for (int j = 0; j < k; ++j) p[j] = ps[j];
        for (int j = 0; j < k; ++j) {
            long long pr = 1;
            for (int i = 0; i <= j; ++i) pm[i][j] = pr, pr = pr * p[i] % p[j];
            inv[j] = mod_inv(pm[j][j], p[j]);
        }
        UInt128 b = 1;
        for (int j = 0; j < k; ++j) base[j] = b, b *= (UInt128)p[j];
    }
    Int128 operator()(const unsigned* rs) const {
        long long t[5];
        UInt128 x = 0;
        for (int j = 0; j < k; ++j) {
            long long cur = 0;  // すでに決めた部分 Σ_{i<j} t_i (p_0 … p_{i−1}) を p_j で割った余り
            for (int i = 0; i < j; ++i) cur = (cur + (t[i] % p[j] + p[j]) % p[j] * pm[i][j]) % p[j];
            long long v = ((long long)rs[j] - cur) % p[j];
            if (v < 0) v += p[j];
            v = v * inv[j] % p[j];
            if (v > p[j] / 2) v -= p[j];  // (−p/2, p/2] に
            t[j] = v;
            x += (UInt128)(Int128)v * base[j];
        }
        return (Int128)x;
    }
};
}  // namespace math_internal

// mint の列の畳み込み c[k] = Σ_{i+j=k} a[i] b[j]。mod が 998244353 などの NTT 素数なら NTT、それ以外（10^9+7・合成数）は
// 3 つの素数で計算してガーナーで戻す（自動）。短い方が 60 以下なら愚直  O((n + m) log(n + m))
template <class mint>
std::vector<mint> convolve(const std::vector<mint>& a, const std::vector<mint>& b) {
    if (a.empty() || b.empty()) return {};
    if (std::min(a.size(), b.size()) <= 60) return math_internal::conv_naive(a, b);
    const size_t z = a.size() + b.size() - 1;
    if (math_internal::ntt_usable<mint>((size_t)1 << math_internal::ceil_log2(z))) return math_internal::conv_ntt(a, b);
    if (z > (1u << 24)) math_internal::fail("convolve", "長さが 2^24 を超える");
    using namespace math_internal;
    std::vector<long long> x(a.size()), y(b.size());
    for (size_t i = 0; i < a.size(); ++i) x[i] = a[i].val();
    for (size_t i = 0; i < b.size(); ++i) y[i] = b[i].val();
    const std::vector<unsigned> r0 = conv_residue<NTT_P0>(x, y), r1 = conv_residue<NTT_P1>(x, y), r2 = conv_residue<NTT_P2>(x, y);
    // 値は 0 以上 min(n, m)·(mod − 1)^2 < P0 P1 P2 なので、ふつうのガーナーで mod の値に
    const long long m = mint::mod();
    const long long i01 = mod_inv(NTT_P0, NTT_P1), i012 = mod_inv((long long)NTT_P0 * NTT_P1 % NTT_P2, NTT_P2);
    const long long p01_m = (long long)NTT_P0 * NTT_P1 % m;
    std::vector<mint> c(z);
    for (size_t i = 0; i < z; ++i) {
        const long long t0 = r0[i];
        const long long t1 = (long long)(r1[i] + NTT_P1 - t0 % NTT_P1) % NTT_P1 * i01 % NTT_P1;
        const long long cur2 = (t0 + (long long)NTT_P0 % NTT_P2 * t1) % NTT_P2;
        const long long t2 = ((long long)r2[i] - cur2 % NTT_P2 + NTT_P2) % NTT_P2 * i012 % NTT_P2;
        c[i] = mint((t0 + (Int128)NTT_P0 * t1 + (Int128)p01_m * t2) % m);
    }
    return c;
}
namespace math_internal {
// 整数のまま畳み込む（結果の絶対値が 2^126 未満と分かっているときだけ。足りるだけの素数を使う）
inline std::vector<Int128> convolve_exact(const std::vector<long long>& a, const std::vector<long long>& b, const char* func) {
    if (a.empty() || b.empty()) return {};
    const size_t z = a.size() + b.size() - 1;
    Int128 ma = 0, mb = 0;
    for (long long v : a) ma = std::max<Int128>(ma, v < 0 ? -(Int128)v : (Int128)v);
    for (long long v : b) mb = std::max<Int128>(mb, v < 0 ? -(Int128)v : (Int128)v);
    // 係数の大きさの上限 ≤ min(n, m)·max|a|·max|b|（long double で見積もる）
    const long double bound = (long double)std::min(a.size(), b.size()) * (long double)ma * (long double)mb;
    if (bound >= 8.0e37L) fail(func, "結果が Int128 に収まらないかもしれない");
    if (std::min(a.size(), b.size()) <= 60) {  // 愚直（128 ビットで。上の見積もりで溢れない）
        std::vector<Int128> c(z, 0);
        for (size_t i = 0; i < a.size(); ++i)
            for (size_t j = 0; j < b.size(); ++j) c[i + j] += (Int128)a[i] * b[j];
        return c;
    }
    const int k = bound < 4.0e24L ? 3 : 5;  // 3 つの素数の積 ≈ 5.9×10^25、5 つ ≈ 1.2×10^44（2^126 ≈ 8.5×10^37 まで）
    if (z > (1u << (k == 3 ? 24 : 23))) fail(func, "長さが大きすぎる");
    const unsigned ps[5] = {NTT_P0, NTT_P1, NTT_P2, NTT_P3, NTT_P4};
    std::vector<unsigned> rs[5];
    rs[0] = conv_residue<NTT_P0>(a, b), rs[1] = conv_residue<NTT_P1>(a, b), rs[2] = conv_residue<NTT_P2>(a, b);
    if (k == 5) rs[3] = conv_residue<NTT_P3>(a, b), rs[4] = conv_residue<NTT_P4>(a, b);
    std::vector<Int128> c(z);
    if (k == 3) {  // 3 つの素数：ふつうのガーナーで [0, M) にしてから M/2 を超えたら M を引く
        const long long P0 = NTT_P0, P1 = NTT_P1, P2 = NTT_P2;
        const long long i01 = mod_inv(P0, P1), i012 = mod_inv(P0 * P1 % P2, P2), p0_2 = P0 % P2;
        const Int128 M = (Int128)P0 * P1 * P2, P01 = (Int128)P0 * P1;
        for (size_t i = 0; i < z; ++i) {
            const long long t0 = rs[0][i];
            const long long t1 = (rs[1][i] - t0 % P1 + P1) % P1 * i01 % P1;
            const long long t2 = ((long long)rs[2][i] - (t0 + p0_2 * t1) % P2 + P2) % P2 * i012 % P2;
            Int128 x = t0 + (Int128)P0 * t1 + P01 * t2;
            if (x > M / 2) x -= M;
            c[i] = x;
        }
        return c;
    }
    const GarnerBalanced g(ps, k);
    for (size_t i = 0; i < z; ++i) {
        unsigned r[5];
        for (int j = 0; j < k; ++j) r[j] = rs[j][i];
        c[i] = g(r);
    }
    return c;
}
}  // namespace math_internal
// 整数の列の畳み込み（厳密。負の数も。結果が long long に収まらなければ止める）
inline std::vector<long long> convolve_ll(const std::vector<long long>& a, const std::vector<long long>& b) {
    const std::vector<Int128> c = math_internal::convolve_exact(a, b, "convolve_ll");
    std::vector<long long> r(c.size());
    for (size_t i = 0; i < c.size(); ++i) {
        if (c[i] > (Int128)9223372036854775807LL || c[i] < -(Int128)9223372036854775807LL - 1) math_internal::fail("convolve_ll", "結果が long long に収まらない（convolve_i128 を使う）");
        r[i] = (long long)c[i];
    }
    return r;
}
// 整数の列の畳み込みを Int128 で（結果の絶対値は 2^126 程度まで）
inline std::vector<Int128> convolve_i128(const std::vector<long long>& a, const std::vector<long long>& b) { return math_internal::convolve_exact(a, b, "convolve_i128"); }
// long long の列を mod m（1 ≤ m < 2^31）で畳み込む
inline std::vector<long long> convolve_mod(const std::vector<long long>& a, const std::vector<long long>& b, long long m) {
    if (m < 1 || m >= (1LL << 31)) math_internal::fail("convolve_mod", "1 ≤ m < 2^31 で使う");
    using M = DynModIntT<-1>;  // 利用者の DynModInt とは別の mod を持つ
    const long long keep = M::mod();
    M::set_mod(m);
    std::vector<M> x(a.size()), y(b.size());
    for (size_t i = 0; i < a.size(); ++i) x[i] = M(a[i]);
    for (size_t i = 0; i < b.size(); ++i) y[i] = M(b[i]);
    const std::vector<M> z = convolve(x, y);
    std::vector<long long> r(z.size());
    for (size_t i = 0; i < z.size(); ++i) r[i] = z[i].val();
    M::set_mod(keep);
    return r;
}
// 多くの多項式の積（短いものから 2 つずつ掛ける）  O(N log^2 N)
template <class mint>
std::vector<mint> convolve_many(std::vector<std::vector<mint>> ps) {
    if (ps.empty()) return {mint(1)};
    for (const auto& p : ps)
        if (p.empty()) return {};
    while (ps.size() > 1) {  // 短い順に並べ、隣どうしを掛けて半分にする（短いものから掛けていく）
        std::sort(ps.begin(), ps.end(), [](const std::vector<mint>& x, const std::vector<mint>& y) { return x.size() < y.size(); });
        std::vector<std::vector<mint>> nx;
        nx.reserve(ps.size() / 2 + 1);
        for (size_t i = 0; i + 1 < ps.size(); i += 2) nx.push_back(convolve(ps[i], ps[i + 1]));
        if (ps.size() % 2) nx.push_back(std::move(ps.back()));
        ps.swap(nx);
    }
    return ps[0];
}
// 相関 c[k] = Σ_{i<|b|} a[k+i]·b[i]（k = 0〜|a|−|b|。a の中で b をずらして重ねた積の和。|a| ≥ |b| ≥ 1）
template <class mint>
std::vector<mint> correlate(const std::vector<mint>& a, const std::vector<mint>& b) {
    if (b.empty() || a.size() < b.size()) math_internal::fail("correlate", "|a| ≥ |b| ≥ 1 で使う");
    const std::vector<mint> rb(b.rbegin(), b.rend());
    const std::vector<mint> c = convolve(a, rb);
    return std::vector<mint>(c.begin() + (long)(b.size() - 1), c.begin() + (long)a.size());
}
// 巡回畳み込み c[k] = Σ_{(i+j) mod n = k} a[i] b[j]（長さ n。a, b の添字は mod n で見る）
template <class mint>
std::vector<mint> convolve_cyclic(const std::vector<mint>& a, const std::vector<mint>& b, size_t n) {
    if (n == 0) math_internal::fail("convolve_cyclic", "n ≥ 1 で使う");
    std::vector<mint> c(n, mint(0));
    const std::vector<mint> d = convolve(a, b);
    for (size_t k = 0; k < d.size(); ++k) c[k % n] += d[k];
    return c;
}
// オンラインの畳み込み：a_n, b_n を 1 つずつ渡すと、その時点で c_n = Σ_{i+j=n} a_i b_j を返す（後の値を知らずに）
//   f_{n+1} = Σ_{i=0}^{n} f_i g_{n−i} のように、答えを次の入力に使う DP に。全体で O(N log^2 N)
template <class mint>
struct RelaxedConvolution {
    mint push(const mint& x, const mint& y) {  // a_n = x、b_n = y を加えて c_n を返す（n = これまでに渡した数）
        const size_t n = a_.size();
        a_.push_back(x), b_.push_back(y);
        if (c_.size() < 2 * n + 2) c_.resize(2 * n + 2, mint(0));
        c_[n] += n == 0 ? a_[0] * b_[0] : a_[n] * b_[0] + a_[0] * b_[n];
        // 添字が 1 以上の部分を大きさ s の正方形に分けて、使う値がそろった時点で掛ける（[s, 2s) × [m s, (m+1) s) と、m ≥ 2 ならその逆）
        for (size_t s = 1; 2 * s <= n + 1 && (n + 1) % s == 0; s *= 2) {
            const size_t m = (n + 1) / s - 1;
            add_block(a_, s, b_, m * s, s, (m + 1) * s);
            if (m >= 2) add_block(b_, s, a_, m * s, s, (m + 1) * s);
        }
        return c_[n];
    }
    size_t size() const { return a_.size(); }  // これまでに渡した数
    mint get(size_t n) const {                // 確定した c_n（n < size()）
        if (n >= a_.size()) math_internal::fail("RelaxedConvolution::get", "n < size() で使う");
        return c_[n];
    }

  private:
    std::vector<mint> a_, b_, c_;
    void add_block(const std::vector<mint>& p, size_t pl, const std::vector<mint>& q, size_t ql, size_t s, size_t out) {  // p[pl, pl+s) × q[ql, ql+s) を c[out + …] に
        if (s <= 64) {
            for (size_t i = 0; i < s; ++i)
                for (size_t j = 0; j < s; ++j) c_[out + i + j] += p[pl + i] * q[ql + j];
            return;
        }
        const std::vector<mint> u(p.begin() + (long)pl, p.begin() + (long)(pl + s)), v(q.begin() + (long)ql, q.begin() + (long)(ql + s));
        const std::vector<mint> w = convolve(u, v);
        for (size_t k = 0; k < w.size(); ++k) c_[out + k] += w[k];
    }
};
#endif  // GG_CONVOLUTION
#ifndef GG_LINEAR_RECURRENCE
#define GG_LINEAR_RECURRENCE
// ---- 線形漸化式（数学第二世代）  a_i = c_1 a_{i−1} + … + c_d a_{i−d} の第 K 項（Bostan–Mori）・最初の数項から漸化式を求める（Berlekamp–Massey） ----
// c は {c_1, c_2, …, c_d}（c[0] が 1 つ前の項の係数）。a は最初の d 項 a_0〜a_{d−1}
namespace math_internal {
template <class mint>
mint linrec_bostan_mori(std::vector<mint> P, std::vector<mint> Q, unsigned long long K) {  // [x^K] P/Q（deg P < deg Q = d、Q[0] = 1）
    const size_t d = Q.size() - 1;
    P.resize(d);
    while (K > 0) {
        const size_t n = (size_t)1 << ceil_log2(2 * d + 1);
        if (d >= 16 && ntt_usable<mint>(n)) {  // NTT の値の並びでは j と j^1 が ±x の組（Q(−x) の値は隣と入れ替えるだけ）
            P.resize(n), Q.resize(n);
            ntt(P), ntt(Q);
            for (size_t j = 0; j < n; j += 2) {
                const mint q0 = Q[j], q1 = Q[j + 1];
                P[j] *= q1, P[j + 1] *= q0;
                Q[j] = Q[j + 1] = q0 * q1;
            }
            intt(P), intt(Q);
        } else {
            std::vector<mint> Qm(Q);
            for (size_t i = 1; i < Qm.size(); i += 2) Qm[i] = -Qm[i];  // Q(−x)
            P = d == 0 ? std::vector<mint>() : convolve(P, Qm), Q = convolve(Q, Qm);
        }
        // P(x) Q(−x) の K の偶奇と同じ次数の項、Q(x) Q(−x) の偶数次の項
        std::vector<mint> np, nq;
        np.reserve(d), nq.reserve(d + 1);
        for (size_t i = (size_t)(K & 1); i < P.size() && np.size() < d; i += 2) np.push_back(P[i]);
        for (size_t i = 0; i < Q.size() && nq.size() < d + 1; i += 2) nq.push_back(Q[i]);
        np.resize(d), nq.resize(d + 1);
        P.swap(np), Q.swap(nq), K >>= 1;
    }
    return d == 0 ? mint(0) : P[0];
}
template <class mint>
std::vector<mint> linrec_series_inv(const std::vector<mint>& Q, size_t n) {  // 1/Q の最初の n 項（Q[0] = 1。掛け算だけのニュートン法なので mod は何でも）
    std::vector<mint> g{mint(1)};
    for (size_t k = 1; k < n;) {
        k *= 2;
        std::vector<mint> q(Q.begin(), Q.begin() + (long long)std::min(Q.size(), k));
        std::vector<mint> t = convolve(q, g);  // g ← g (2 − Q g)
        t.resize(k);
        for (auto& x : t) x = -x;
        t[0] += mint(2);
        g = convolve(g, t);
        g.resize(k);
    }
    g.resize(n);
    return g;
}
template <class mint>
void linrec_check(const std::vector<mint>& a, const std::vector<mint>& c, const char* func) {
    if (a.size() < c.size()) fail(func, "a は c と同じ数（d 個）以上の項が必要");
}
}  // namespace math_internal

// 第 K 項 a_K（a は最初の d 項以上、c = {c_1, …, c_d}。K ≤ 10^18。mod は何でも）  O(d log d log K)
template <class mint>
mint linear_recurrence_kth(const std::vector<mint>& a, const std::vector<mint>& c, long long K) {
    math_internal::linrec_check(a, c, "linear_recurrence_kth");
    if (K < 0) math_internal::fail("linear_recurrence_kth", "K ≥ 0 で使う");
    if (K < (long long)a.size()) return a[(size_t)K];
    const size_t d = c.size();
    if (d == 0) return mint(0);
    std::vector<mint> Q(d + 1), A(a.begin(), a.begin() + (long long)d);  // Q = 1 − c_1 x − … − c_d x^d、P = A Q mod x^d
    Q[0] = 1;
    for (size_t j = 0; j < d; ++j) Q[j + 1] = -c[j];
    std::vector<mint> P = convolve(A, Q);
    P.resize(d);
    return math_internal::linrec_bostan_mori(P, Q, (unsigned long long)K);
}
// 最初の n 項 a_0〜a_{n−1}  O(n log n)（d が小さければ O(n d)）
template <class mint>
std::vector<mint> linear_recurrence_terms(const std::vector<mint>& a, const std::vector<mint>& c, long long n) {
    math_internal::linrec_check(a, c, "linear_recurrence_terms");
    if (n < 0) math_internal::fail("linear_recurrence_terms", "n ≥ 0 で使う");
    const size_t d = c.size(), N = (size_t)n;
    std::vector<mint> res(a.begin(), a.begin() + (long long)std::min(a.size(), N));
    if (res.size() == N) return res;
    if (d == 0) {
        res.resize(N);
        return res;
    }
    if ((double)d * (double)N <= 3e7 || d <= 64) {  // 漸化式どおりに
        res.resize(N);
        for (size_t i = a.size(); i < N; ++i) {
            mint s = 0;
            for (size_t j = 0; j < d; ++j) s += c[j] * res[i - 1 - j];
            res[i] = s;
        }
        return res;
    }
    std::vector<mint> Q(d + 1), A(a.begin(), a.begin() + (long long)d);
    Q[0] = 1;
    for (size_t j = 0; j < d; ++j) Q[j + 1] = -c[j];
    std::vector<mint> P = convolve(A, Q);
    P.resize(d);
    std::vector<mint> r = convolve(P, math_internal::linrec_series_inv(Q, N));
    r.resize(N);
    return r;
}
// 列 s をみたす最短の線形漸化式の c = {c_1, …, c_d}（mod は素数。2d 項以上あれば正しい漸化式が出る）  O(N^2)
template <class mint>
std::vector<mint> berlekamp_massey(const std::vector<mint>& s) {
    std::vector<mint> C{mint(1)}, B{mint(1)};
    mint b = 1;
    size_t L = 0, m = 1;
    for (size_t n = 0; n < s.size(); ++n) {
        mint dlt = s[n];
        for (size_t i = 1; i <= L; ++i) dlt += C[i] * s[n - i];
        if (dlt == mint(0)) {
            ++m;
            continue;
        }
        const mint coef = dlt / b;
        if (2 * L <= n) {
            const std::vector<mint> T = C;
            if (C.size() < B.size() + m) C.resize(B.size() + m, mint(0));
            for (size_t i = 0; i < B.size(); ++i) C[i + m] -= coef * B[i];
            L = n + 1 - L, B = T, b = dlt, m = 1;
        } else {
            if (C.size() < B.size() + m) C.resize(B.size() + m, mint(0));
            for (size_t i = 0; i < B.size(); ++i) C[i + m] -= coef * B[i];
            ++m;
        }
    }
    C.resize(L + 1, mint(0));
    std::vector<mint> c(L);
    for (size_t j = 1; j <= L; ++j) c[j - 1] = -C[j];
    return c;
}
// 最初の数項 s から漸化式を推測して第 K 項（漸化式の次数 d に対して 2d 項以上いる。K ≤ 10^18）
template <class mint>
mint guess_kth(const std::vector<mint>& s, long long K) {
    if (K < 0) math_internal::fail("guess_kth", "K ≥ 0 で使う");
    if (K < (long long)s.size()) return s[(size_t)K];
    return linear_recurrence_kth(s, berlekamp_massey(s), K);
}
#endif  // GG_LINEAR_RECURRENCE
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace linrec_test {

#define LINREC_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[linear_recurrence] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
static std::vector<mint> naive_terms(const std::vector<mint>& a, const std::vector<mint>& c, size_t n) {
    std::vector<mint> r(a.begin(), a.begin() + (long long)std::min(a.size(), n));
    r.resize(n);
    for (size_t i = a.size(); i < n; ++i) {
        mint s = 0;
        for (size_t j = 0; j < c.size(); ++j) s += c[j] * r[i - 1 - j];
        r[i] = s;
    }
    return r;
}
template <class mint>
static std::vector<mint> rand_vec(std::mt19937_64& rng, size_t n) {
    std::vector<mint> v(n);
    for (auto& x : v) x = rng() % 4 == 0 ? mint(rnd(rng, -2, 2)) : mint((long long)(rng() >> 2));
    return v;
}

template <class mint>
static void test_kth(std::mt19937_64& rng, int& cases, int iters) {
    for (int it = 0; it < iters; ++it) {
        const size_t d = (size_t)(it % 3 == 0 ? rnd(rng, 0, 6) : rnd(rng, 1, it % 2 ? 40 : 200));
        const std::vector<mint> c = rand_vec<mint>(rng, d), a = naive_terms(rand_vec<mint>(rng, d), c, d + (size_t)rnd(rng, 0, 2));  // 多めに渡してもよい
        const size_t N = d * 3 + 50;
        const std::vector<mint> t = naive_terms(a, c, N);
        for (int q = 0; q < 6; ++q, ++cases) {
            const long long K = rnd(rng, 0, (long long)N - 1);
            LINREC_CHECK(linear_recurrence_kth(a, c, K) == t[(size_t)K]);
        }
        ++cases;
        LINREC_CHECK(linear_recurrence_terms(a, c, (long long)N) == t);
        // K が巨大：t 項ずらした列で同じ値になるか
        const long long K = rnd(rng, 0, 1000000000000000000LL), sh = rnd(rng, 0, (long long)N - (long long)d);
        const std::vector<mint> a2(t.begin() + sh, t.begin() + sh + (long long)d);
        if (K >= sh) LINREC_CHECK(linear_recurrence_kth(a2, c, K - sh) == linear_recurrence_kth(a, c, K));
    }
}

template <class mint>
static void test_bm(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 80 * TEST_SCALE; ++it, ++cases) {
        const size_t d = (size_t)rnd(rng, 0, it % 2 ? 10 : 80);
        const std::vector<mint> c = rand_vec<mint>(rng, d), a = rand_vec<mint>(rng, d);
        const std::vector<mint> s = naive_terms(a, c, 2 * d + (size_t)rnd(rng, 0, 5));
        const std::vector<mint> g = berlekamp_massey(s);
        LINREC_CHECK(g.size() <= d);
        const std::vector<mint> t = naive_terms(a, c, 2 * d + 100);
        LINREC_CHECK(naive_terms(std::vector<mint>(s.begin(), s.begin() + (long long)g.size()), g, t.size()) == t);
        const long long K = rnd(rng, 0, 1000000000000000000LL);
        LINREC_CHECK(guess_kth(s, K) == linear_recurrence_kth(a, c, K));
    }
    // n^3 は (1 − x)^4 の漸化式
    std::vector<mint> cube;
    for (long long n = 0; n < 12; ++n) cube.push_back(mint(n * n * n));
    LINREC_CHECK(berlekamp_massey(cube) == std::vector<mint>({mint(4), mint(-6), mint(4), mint(-1)}));
    LINREC_CHECK(guess_kth(cube, 1000000) == mint(1000000).pow(3));
    LINREC_CHECK(berlekamp_massey(std::vector<mint>(10, mint(0))).empty() && berlekamp_massey(std::vector<mint>()).empty());
    cases += 3;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_kth<ModInt998>(rng, cases, 60 * TEST_SCALE);
    test_kth<ModInt107>(rng, cases, 25 * TEST_SCALE);  // NTT が使えない素数
    DynModInt::set_mod(1000000000);
    test_kth<DynModInt>(rng, cases, 15 * TEST_SCALE);  // 合成数の mod でも第 K 項・最初の n 項は使える
    DynModInt::set_mod(998244353);
    test_kth<DynModInt>(rng, cases, 8 * TEST_SCALE);
    test_bm<ModInt998>(rng, cases);
    test_bm<ModInt107>(rng, cases);
    {  // フィボナッチ数：倍数公式と
        using mint = ModInt998;
        auto fib = [](unsigned long long n) {  // (F_n, F_{n+1})
            mint x = 0, y = 1;
            for (int b = 63; b >= 0; --b) {
                const mint x2 = x * (y * mint(2) - x), y2 = x * x + y * y;
                x = x2, y = y2;
                if (n >> b & 1) {
                    const mint t = x + y;
                    x = y, y = t;
                }
            }
            return x;
        };
        for (int it = 0; it < 1000; ++it, ++cases) {
            const long long K = rnd(rng, 0, 1000000000000000000LL);
            LINREC_CHECK(linear_recurrence_kth<mint>({0, 1}, {1, 1}, K) == fib((unsigned long long)K));
        }
    }
    {  // 長い列（O(n log n) の方）：漸化式どおりの値と
        using mint = ModInt998;
        for (int it = 0; it < 2; ++it, ++cases) {
            const size_t d = 100, N = 310000;
            const std::vector<mint> c = rand_vec<mint>(rng, d), a = rand_vec<mint>(rng, d);
            LINREC_CHECK(linear_recurrence_terms(a, c, (long long)N) == naive_terms(a, c, N));
        }
        const std::vector<ModInt107> c7 = rand_vec<ModInt107>(rng, 120), a7 = rand_vec<ModInt107>(rng, 120);
        LINREC_CHECK(linear_recurrence_terms(a7, c7, 260000) == naive_terms(a7, c7, 260000));
    }
    std::cout << "[linear_recurrence] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[linear_recurrence] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    using mint = ModInt998;
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::vector<mint> c = rand_vec<mint>(rng, 50000), a = rand_vec<mint>(rng, 50000);
    auto t0 = std::chrono::steady_clock::now();
    const mint x = linear_recurrence_kth(a, c, 1000000000000000000LL);
    const long long t_k = ms(t0);
    const std::vector<ModInt107> c7 = rand_vec<ModInt107>(rng, 10000), a7 = rand_vec<ModInt107>(rng, 10000);
    t0 = std::chrono::steady_clock::now();
    const ModInt107 y = linear_recurrence_kth(a7, c7, 1000000000000000000LL);
    const long long t_k7 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<mint> t = linear_recurrence_terms(a, c, 1000000);
    const long long t_t = ms(t0);
    const std::vector<mint> s = linear_recurrence_terms(std::vector<mint>(a.begin(), a.begin() + 2000), std::vector<mint>(c.begin(), c.begin() + 2000), 4000);
    t0 = std::chrono::steady_clock::now();
    const std::vector<mint> g = berlekamp_massey(s);
    const long long t_bm = ms(t0);
    std::cout << "[linear_recurrence] benchmark: kth（d = 5×10^4, K = 10^18）" << t_k << " ms, kth（10^9+7, d = 10^4）" << t_k7 << " ms, terms（d = 5×10^4, n = 10^6）"
              << t_t << " ms, berlekamp_massey（4000 項）" << t_bm << " ms  (check " << (x.val() + y.val() + t.back().val() + g.size()) % 7 << ")\n";
}

}  // namespace linrec_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace linrec_demo {

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
#define SHOW(expr, comment) linrec_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    linrec_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    using mint = ModInt998;
    banner("線形漸化式（linear_recurrence） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  linear_recurrence_kth(a, c, K) / linear_recurrence_terms(a, c, n) / berlekamp_massey(s) / guess_kth(s, K)\n"
              << "  a_i = c_1 a_{i−1} + … + c_d a_{i−d}、a = {a_0, …, a_{d−1}}、c = {c_1, …, c_d}（この中では using mint = ModInt998;）\n";

    api("【1】 linear_recurrence_kth(a, c, K)", "mint", "第 K 項（K ≤ 10^18。Bostan–Mori。mod は何でも）", "O(d log d log K)");
    SHOW(linear_recurrence_kth<mint>({0, 1}, {1, 1}, 10), "フィボナッチ数 F_10");
    SHOW(linear_recurrence_kth<mint>({0, 1}, {1, 1}, 1000000000000000000LL), "F_{10^18} mod 998244353");
    SHOW(linear_recurrence_kth<mint>({0, 0, 1}, {1, 1, 1}, 30), "トリボナッチ数");
    SHOW(linear_recurrence_kth<ModInt107>({1, 3}, {4, -1}, 1000000000000000000LL), "a_i = 4a_{i−1} − a_{i−2}（mod 10^9+7）");

    api("【2】 linear_recurrence_terms(a, c, n)", "vector<mint>", "最初の n 項", "O(n log n)（d が小さいときは O(n d)）");
    SHOW(linear_recurrence_terms<mint>({0, 1}, {1, 1}, 12), "F_0〜F_11");
    SHOW(linear_recurrence_terms<mint>({1, 1, 1}, {0, 1, 1}, 12), "パドヴァン数列 a_i = a_{i−2} + a_{i−3}");

    api("【3】 berlekamp_massey(s) / guess_kth(s, K)", "vector<mint> / mint",
        "列をみたす最短の漸化式の c（次数 d には 2d 項以上いる）/ それで第 K 項", "O(N^2) + O(d log d log K)");
    SHOW(berlekamp_massey(std::vector<mint>({0, 1, 1, 2, 3, 5, 8, 13})), "フィボナッチ：c = {1, 1}");
    SHOW(berlekamp_massey(std::vector<mint>({0, 1, 8, 27, 64, 125, 216, 343, 512, 729})), "n^3：c = {4, −6, 4, −1}（(1 − x)^4）");
    {
        // 3 × n の盤面をドミノで敷き詰める数：列ごとの DP で最初の 40 項を出し、guess_kth で n = 10^18
        std::vector<mint> f;
        for (int n = 0; n < 40; ++n) {
            std::vector<mint> dp(8, mint(0));  // 列の境目で、次の列にはみ出しているマス（3 ビット）
            dp[0] = 1;
            for (int col = 0; col < n; ++col) {
                std::vector<mint> nd(8, mint(0));
                for (int s = 0; s < 8; ++s) {
                    if (dp[(size_t)s] == mint(0)) continue;
                    std::function<void(int, int, int)> fill = [&](int r, int cur, int nxt) {  // r 行目から置く。cur は埋まった、nxt ははみ出し
                        if (r == 3) { nd[(size_t)nxt] += dp[(size_t)s]; return; }
                        if (cur >> r & 1) { fill(r + 1, cur, nxt); return; }
                        fill(r + 1, cur | 1 << r, nxt | 1 << r);  // 横に置く
                        if (r + 1 < 3 && !(cur >> (r + 1) & 1)) fill(r + 2, cur | 3 << r, nxt);  // 縦に置く
                    };
                    fill(0, s, 0);
                }
                dp.swap(nd);
            }
            f.push_back(dp[0]);
        }
        note("3 × n をドミノで敷き詰める数：DP で最初の 40 項 → berlekamp_massey → guess_kth で n = 10^18");
        std::cout << "    最初の項 " << to_s(std::vector<mint>(f.begin(), f.begin() + 9)) << " …\n";
        std::cout << "    漸化式 c = " << to_s(berlekamp_massey(f)) << "（a_n = 4a_{n−2} − a_{n−4}）\n";
        std::cout << "    n = 10^18 → " << guess_kth(f, 1000000000000000000LL) << "\n";
    }
    {
        // ABC300 Ex（Fibonacci: Revisited）：Σ_{m ⊆ N} a_m。P/Q の Bostan–Mori で、N のビットが 1 なら偶奇両方の項を足す
        auto solve = [](int K, long long N) {
            std::vector<mint> Q((size_t)K + 1), A((size_t)K, mint(1));  // Q = 1 − x − … − x^K、P = A Q mod x^K
            Q[0] = 1;
            for (int j = 1; j <= K; ++j) Q[(size_t)j] = -1;
            std::vector<mint> P = convolve(A, Q);
            P.resize((size_t)K);
            for (; N > 0; N >>= 1) {
                std::vector<mint> Qm(Q);
                for (size_t i = 1; i < Qm.size(); i += 2) Qm[i] = -Qm[i];
                const std::vector<mint> U = convolve(P, Qm), V = convolve(Q, Qm);
                std::vector<mint> nP((size_t)K, mint(0)), nQ((size_t)K + 1, mint(0));
                for (size_t i = 0; i < U.size(); ++i)
                    if (i % 2 == 0 || (N & 1)) nP[i / 2] += U[i];  // N のビットが 0 なら偶数次だけ、1 なら両方
                for (size_t i = 0; i < V.size(); i += 2) nQ[i / 2] = V[i];
                P.swap(nP), Q.swap(nQ);
            }
            return P[0] / Q[0];
        };
        note("ABC300 Ex（Fibonacci: Revisited）: Bostan–Mori で、N のビットが 1 の段は偶数次と奇数次の項を両方足す（convolve で）");
        std::cout << "    (2, 6) → " << solve(2, 6) << "、(2, 8) → " << solve(2, 8) << "、(1, 123456789) → " << solve(1, 123456789) << "、(300, 20230429) → "
                  << solve(300, 20230429) << "\n    (42923, 999999999558876113) → " << solve(42923, 999999999558876113LL) << "\n";
    }
}

}  // namespace linrec_demo

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
    linrec_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    linrec_test::benchmark();    // 速度計測
    linrec_demo::demo();         // 使い方デモ
    return 0;
}
