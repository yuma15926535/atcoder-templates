// =============================================================================
//  合同式・一次不定方程式・中国剰余定理（a x ≡ b (mod m)・a x + b y = c・x ≡ r_i (mod m_i)）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 03 合同式・方程式
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  答えの形
//    ModSolution{ok, r, m}            解は x ≡ r (mod m)（0 ≤ r < m）。ok が false なら解なし
//    Diophantine{ok, x0, y0, dx, dy}  a x + b y = c の解すべて：x = x0 + k·dx、y = y0 − k·dy（d.x(k) / d.y(k) でも）
//  （1 つの合同式の逆元・累乗・個数は 00_math_base：mod_inv / mod_pow / ceil_congruent / count_congruent。このファイルに入っている）
//
//  機能一覧                                                                    計算量
//    linear_congruence(a, b, m)       a x ≡ b (mod m) の解を x ≡ r (mod m / gcd(a, m)) で（a, b は負でも）      O(log m)
//    linear_diophantine(a, b, c)      a x + b y = c の解すべて（0 ≤ x0 < dx = |b| / g）                       O(log)
//    diophantine_count(a, b, c, xl, xr, yl, yr)        xl ≤ x < xr かつ yl ≤ y < yr の解の数
//    diophantine_count_closed(a, b, c, XL, XR, YL, YR)  XL ≤ x ≤ XR かつ YL ≤ y ≤ YR の解の数（x, y ≥ 0 なら範囲は 0〜LLONG_MAX）
//                                     （2 つとも、個数が long long に収まらなければ止める）
//    diophantine_min_x(a, b, c, x_min)                 x ≥ x_min で x が最小の解 DiophantinePoint{ok, x, y}
//    crt_merge(r1, m1, r2, m2)        x ≡ r1 (mod m1) かつ x ≡ r2 (mod m2)（互いに素でなくてよい）          O(log)
//    crt_all(rs, ms)                  x ≡ rs[i] (mod ms[i]) をすべて。答えの r は 0 以上で最小の解            O(n log)
//                                     lcm が long long を超えたら overflow = true（ok は「long long に入る解があるか」、r はその最小）
//    garner(rs, ms, mod)              互いに素な ms での答え x（0 ≤ x < Π ms。long long を超えてよい）を mod で  O(n^2)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    meq_test::random_test()  … 全部調べる愚直解・128 ビットでの検算との突き合わせ
//    meq_test::benchmark()    … 速度計測
//    meq_demo::demo()         … 全機能の使い方デモ（日本語の説明付き。周期的に繰り返す 2 つの区間が重なる最初の時刻も）
// =============================================================================
#include <algorithm>
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
#ifndef GG_MOD_EQUATION
#define GG_MOD_EQUATION
// ---- 合同式・一次不定方程式・中国剰余定理（数学第二世代） ----
struct ModSolution {  // 解は x ≡ r (mod m)（0 ≤ r < m）。ok が false なら解なし
    bool ok = false;
    long long r = 0, m = 1;
    bool overflow = false;  // crt で lcm が long long を超えた。そのとき ok は「long long に入る解があるか」、r はその最小。m は 0（使わない）
};
struct Diophantine {  // a x + b y = c の整数解すべて：x = x0 + k·dx、y = y0 − k·dy（k は整数）。ok が false なら解なし
    bool ok = false;
    long long x0 = 0, y0 = 0, dx = 0, dy = 0;  // dx ≥ 0。b ≠ 0 なら 0 ≤ x0 < dx（x0 は 0 以上で最小の x）
    long long x(long long k) const { return check(x0 + (Int128)k * dx); }  // k 番目の解の x（収まらなければ止める）
    long long y(long long k) const { return check(y0 - (Int128)k * dy); }

  private:
    static long long check(Int128 v) {
        if (v > (Int128)9223372036854775807LL || v < -(Int128)9223372036854775807LL - 1) math_internal::fail("Diophantine", "long long に収まらない");
        return (long long)v;
    }
};
struct DiophantinePoint {  // a x + b y = c の 1 つの解（ok が false なら条件を満たす解なし）
    bool ok = false;
    long long x = 0, y = 0;
};

namespace math_internal {
inline bool fits_ll(Int128 v) { return -(Int128)9223372036854775807LL - 1 <= v && v <= (Int128)9223372036854775807LL; }
// x ≡ r1 (mod m1) と x ≡ r2 (mod m2)（0 ≤ r1 < m1 ≤ 2^63 − 1、0 ≤ r2 < m2 ≤ 2^63 − 1）を合わせる。r < m = lcm < 2^126
struct Crt128 {
    bool ok;
    Int128 r, m;
};
inline Crt128 crt_merge128(long long r1, long long m1, long long r2, long long m2) {
    if (m1 < m2) std::swap(r1, r2), std::swap(m1, m2);
    if (m1 % m2 == 0) return Crt128{r1 % m2 == r2, r1, m1};  // m2 | m1（よくある形を速く）
    const ExtGcd e = ext_gcd(m1, m2);  // m1·e.x + m2·e.y = g なので (m1 / g)·e.x ≡ 1 (mod m2 / g)
    const long long g = e.g, diff = r2 - r1;  // |diff| < 2^63
    if (diff % g != 0) return Crt128{false, 0, 1};
    const long long mm = m2 / g;  // r1 + m1·t ≡ r2 (mod m2) ⇔ (m1 / g)·t ≡ diff / g (mod mm)
    const long long t = (long long)((Int128)mod_floor(diff / g, mm) * mod_floor(e.x, mm) % mm);
    return Crt128{true, (Int128)r1 + (Int128)m1 * t, (Int128)m1 * mm};
}
// 以下 crt_all の本体：(r, m) に (r2, m2) を足していく。lcm が long long を超えたら「r がそのまま解か」だけを見る
struct CrtState {
    bool ok = true, overflow = false;
    long long r = 0, m = 1;
    void add(long long r2, long long m2, const char* func) {
        if (m2 < 1) fail(func, "mod は 1 以上で使う");
        if (!ok) return;
        r2 = mod_floor(r2, m2);
        if (overflow) {  // 解は r + t·(lcm)（lcm > 2^63）。r が合わなければ次の解はもう long long に入らない
            if (r % m2 != r2) ok = false;
            return;
        }
        const Crt128 c = crt_merge128(r, m, r2, m2);
        if (!c.ok) { ok = false; return; }
        if (c.m > (Int128)9223372036854775807LL) {
            overflow = true;
            if (c.r > (Int128)9223372036854775807LL) { ok = false; return; }
            r = (long long)c.r, m = 0;
            return;
        }
        r = (long long)c.r, m = (long long)c.m;
    }
    ModSolution result() const {
        if (!ok) return ModSolution{false, 0, overflow ? 0 : 1, overflow};
        return ModSolution{true, r, overflow ? 0 : m, overflow};
    }
};
}  // namespace math_internal

// ---- 一次合同式 a x ≡ b (mod m) ----
inline ModSolution linear_congruence(long long a, long long b, long long m) {  // 解を x ≡ r (mod m / g) の形で（g = gcd(a, m)。a, b は負でも）
    if (m < 1) math_internal::fail("linear_congruence", "m ≥ 1 で使う");
    a = mod_floor(a, m), b = mod_floor(b, m);
    const ExtGcd e = ext_gcd(a, m);  // a·e.x + m·e.y = g なので (a / g)·e.x ≡ 1 (mod m / g)。a = 0 なら g = m
    const long long g = e.g;
    if (b % g != 0) return ModSolution{false, 0, 1, false};
    const long long m2 = m / g;
    return ModSolution{true, (long long)((Int128)(b / g) * mod_floor(e.x, m2) % m2), m2, false};
}

// ---- 一次不定方程式 a x + b y = c ----
namespace math_internal {
struct Dio128 {  // Diophantine と同じ形（128 ビット。y0 が long long に入らないときにも使える）
    bool ok = false;
    Int128 x0 = 0, y0 = 0, dx = 0, dy = 0;
};
inline Dio128 dio_solve(long long a, long long b, long long c, const char* func) {
    const long long LMIN = -9223372036854775807LL - 1;
    if (a == LMIN || b == LMIN) fail(func, "a, b に LLONG_MIN は使えない");
    Dio128 d;
    if (a == 0 && b == 0) {
        d.ok = c == 0;
        return d;
    }
    if (b == 0) {
        if ((Int128)c % a != 0) return d;
        d.ok = true, d.x0 = (Int128)c / a, d.y0 = 0, d.dx = 0, d.dy = -1;
        return d;
    }
    const ExtGcd e = ext_gcd(a, b);  // a·e.x + b·e.y = g
    const long long g = e.g;
    if (c % g != 0) return d;
    d.ok = true;
    d.dx = (b < 0 ? -(Int128)b : (Int128)b) / g;
    d.dy = (Int128)(a / g) * (b < 0 ? -1 : 1);
    d.x0 = (Int128)e.x * (c / g) % d.dx;  // |e.x| ≤ |b| / g なので収まる
    if (d.x0 < 0) d.x0 += d.dx;
    d.y0 = ((Int128)c - (Int128)a * d.x0) / b;  // 割り切れる
    return d;
}
}  // namespace math_internal
// 解を x = x0 + k·dx、y = y0 − k·dy で（g = gcd(a, b)、dx = |b| / g、dy = a / g × (b の符号)）。0 ≤ x0 < dx
//   b = 0 なら x は 1 つに決まり y は何でもよい（dx = 0、dy = −1、y0 = 0）。a = 0 なら x が何でもよい（dx = 1、dy = 0）
//   a = b = 0 なら c = 0 のとき x, y は何でもよい（ok = true、dx = dy = 0）。a, b に LLONG_MIN は使えない
inline Diophantine linear_diophantine(long long a, long long b, long long c) {
    const math_internal::Dio128 d = math_internal::dio_solve(a, b, c, "linear_diophantine");
    if (!math_internal::fits_ll(d.x0) || !math_internal::fits_ll(d.y0)) math_internal::fail("linear_diophantine", "x0 / y0 が long long に収まらない");
    return Diophantine{d.ok, (long long)d.x0, (long long)d.y0, (long long)d.dx, (long long)d.dy};
}
namespace math_internal {
// a x + b y = c、xl ≤ x < xr、yl ≤ y < yr の解の数（2^100 を超えたら 2^100）
inline Int128 dio_count(long long a, long long b, long long c, Int128 xl, Int128 xr, Int128 yl, Int128 yr, const char* func) {
    const Dio128 d = dio_solve(a, b, c, func);
    if (xl >= xr || yl >= yr || !d.ok) return 0;
    const Int128 INF = (Int128)1 << 100;
    if (d.dx == 0 && d.dy == 0) {  // a = b = 0, c = 0：どの組も解
        const Int128 w = xr - xl, h = yr - yl;
        return w > INF / h ? INF : w * h;
    }
    Int128 lo = -INF, hi = INF;  // k の範囲 [lo, hi)
    if (d.dx == 0) {
        if (!(xl <= d.x0 && d.x0 < xr)) return 0;
    } else {  // xl ≤ x0 + k dx < xr
        lo = std::max(lo, ceil_div128(xl - d.x0, d.dx));
        hi = std::min(hi, ceil_div128(xr - d.x0, d.dx));
    }
    if (d.dy == 0) {
        if (!(yl <= d.y0 && d.y0 < yr)) return 0;
    } else if (d.dy > 0) {  // yl ≤ y0 − k dy < yr ⇔ (y0 − yr) / dy < k ≤ (y0 − yl) / dy
        lo = std::max(lo, floor_div128(d.y0 - yr, d.dy) + 1);
        hi = std::min(hi, floor_div128(d.y0 - yl, d.dy) + 1);
    } else {  // y = y0 + k e（e = −dy > 0）
        lo = std::max(lo, ceil_div128(yl - d.y0, -d.dy));
        hi = std::min(hi, ceil_div128(yr - d.y0, -d.dy));
    }
    return hi > lo ? hi - lo : 0;
}
}  // namespace math_internal
inline long long diophantine_count(long long a, long long b, long long c, long long xl, long long xr, long long yl, long long yr) {  // xl ≤ x < xr、yl ≤ y < yr の解の数
    const Int128 n = math_internal::dio_count(a, b, c, xl, xr, yl, yr, "diophantine_count");
    if (n > (Int128)9223372036854775807LL) math_internal::fail("diophantine_count", "個数が long long に収まらない");
    return (long long)n;
}
inline long long diophantine_count_closed(long long a, long long b, long long c, long long XL, long long XR, long long YL, long long YR) {  // XL ≤ x ≤ XR、YL ≤ y ≤ YR
    const Int128 n = math_internal::dio_count(a, b, c, XL, (Int128)XR + 1, YL, (Int128)YR + 1, "diophantine_count_closed");
    if (n > (Int128)9223372036854775807LL) math_internal::fail("diophantine_count_closed", "個数が long long に収まらない");
    return (long long)n;
}
// x ≥ x_min で x が最小の解（b = 0 なら y = 0 を返す。x が long long に入らなければ ok = false。y が入らなければ止める）
inline DiophantinePoint diophantine_min_x(long long a, long long b, long long c, long long x_min) {
    const math_internal::Dio128 d = math_internal::dio_solve(a, b, c, "diophantine_min_x");
    if (!d.ok) return DiophantinePoint{};
    Int128 x = x_min, y = 0;  // a = b = c = 0 ならこのまま
    if (d.dx == 0 && d.dy != 0) {  // b = 0：x は 1 つ
        if (d.x0 < x_min) return DiophantinePoint{};
        x = d.x0;
    } else if (d.dx != 0) {
        const Int128 k = math_internal::ceil_div128((Int128)x_min - d.x0, d.dx);
        x = d.x0 + k * d.dx, y = d.y0 - k * d.dy;
    }
    if (x > (Int128)9223372036854775807LL) return DiophantinePoint{};  // long long に入る解がない
    if (!math_internal::fits_ll(y)) math_internal::fail("diophantine_min_x", "y が long long に収まらない");
    return DiophantinePoint{true, (long long)x, (long long)y};
}

// ---- 中国剰余定理（mod は互いに素でなくてよい。答えは 0 ≤ r < lcm で最小） ----
inline ModSolution crt_merge(long long r1, long long m1, long long r2, long long m2) {  // x ≡ r1 (mod m1) かつ x ≡ r2 (mod m2)
    math_internal::CrtState s;
    s.add(r1, m1, "crt_merge");
    s.add(r2, m2, "crt_merge");
    return s.result();
}
inline ModSolution crt_all(const std::vector<long long>& rs, const std::vector<long long>& ms) {  // x ≡ rs[i] (mod ms[i]) をすべて（空なら x ≡ 0 (mod 1)）
    if (rs.size() != ms.size()) math_internal::fail("crt_all", "rs と ms の長さが違う");
    math_internal::CrtState s;
    for (size_t i = 0; i < rs.size(); ++i) s.add(rs[i], ms[i], "crt_all");
    return s.result();
}
// x ≡ rs[i] (mod ms[i])（ms は互いに素）を満たす 0 ≤ x < Π ms の x を mod で（Π ms が long long を超えてよい）  O(n^2)
inline long long garner(const std::vector<long long>& rs, const std::vector<long long>& ms, long long mod) {
    if (rs.size() != ms.size()) math_internal::fail("garner", "rs と ms の長さが違う");
    if (mod < 1) math_internal::fail("garner", "mod ≥ 1 で使う");
    for (long long m : ms)
        if (m < 1) math_internal::fail("garner", "mod は 1 以上で使う");
    const size_t n = ms.size();
    std::vector<long long> coef(n + 1, 1), cons(n + 1, 0), mm(ms);
    mm.push_back(mod);
    for (size_t i = 0; i < n; ++i) {
        const Maybe<long long> inv = try_mod_inv(coef[i], ms[i]);
        if (!inv.ok) math_internal::fail("garner", "ms が互いに素でない");
        const long long t = mod_mul(mod_floor(rs[i], ms[i]) - cons[i], inv.value, ms[i]);
        for (size_t j = i + 1; j <= n; ++j) {
            cons[j] = (long long)(((Int128)cons[j] + (Int128)t * coef[j]) % mm[j]);
            coef[j] = (long long)((Int128)coef[j] * ms[i] % mm[j]);
        }
    }
    return cons[n];
}
#endif  // GG_MOD_EQUATION
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace meq_test {

#define MEQ_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[mod_equation] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
static const long long LMIN = std::numeric_limits<long long>::min();
static unsigned long long bf_gcd(unsigned long long a, unsigned long long b) {
    while (b) {
        const unsigned long long t = a % b;
        a = b, b = t;
    }
    return a;
}
static unsigned long long uab(long long a) { return a < 0 ? 0ULL - (unsigned long long)a : (unsigned long long)a; }
static Int128 mod128(Int128 a, Int128 m) {
    a %= m;
    return a < 0 ? a + m : a;
}
static long long no_lmin(long long v) { return v == LMIN ? LMIN + 1 : v; }

static void test_linear_congruence(std::mt19937_64& rng, int& cases) {
    for (long long m = 1; m <= 40; ++m)
        for (int it = 0; it < 300; ++it, ++cases) {
            const long long a = it % 7 == 0 ? rnd_ll(rng) : rnd(rng, -100, 100), b = it % 5 == 0 ? rnd_ll(rng) : rnd(rng, -100, 100);
            const ModSolution s = linear_congruence(a, b, m);
            std::vector<long long> sols, expect;
            for (long long x = 0; x < m; ++x)
                if (mod128((Int128)a * x - b, m) == 0) sols.push_back(x);
            MEQ_CHECK(s.ok == !sols.empty() && !s.overflow);
            if (!s.ok) continue;
            MEQ_CHECK(0 <= s.r && s.r < s.m && m % s.m == 0);
            for (long long x = 0; x < m; ++x)
                if (x % s.m == s.r) expect.push_back(x);
            MEQ_CHECK(sols == expect);
        }
    for (int it = 0; it < 300000 * TEST_SCALE; ++it, ++cases) {
        const long long m = it % 3 == 0 ? rnd(rng, 1, LMAX) : (it % 3 == 1 ? rnd(rng, 1, 1000000000) : std::max<long long>(1, (long long)uab(no_lmin(rnd_ll(rng)))));
        long long a = rnd_ll(rng), b = rnd_ll(rng);
        if (it % 4 == 0) a = a / 1000 * 6, b = b / 1000 * 6;  // 公約数を持たせる
        const ModSolution s = linear_congruence(a, b, m);
        const unsigned long long g = bf_gcd((unsigned long long)mod128(a, m), (unsigned long long)m);
        MEQ_CHECK(s.ok == (mod128(b, m) % (Int128)g == 0));
        if (s.ok) MEQ_CHECK(s.m == (long long)((unsigned long long)m / g) && 0 <= s.r && s.r < s.m && mod128((Int128)a * s.r - b, m) == 0);
    }
}

static void test_diophantine(std::mt19937_64& rng, int& cases) {
    // 小さい値：窓の中の解を全部調べる
    for (int it = 0; it < 25000 * TEST_SCALE; ++it, ++cases) {
        const long long a = rnd(rng, -12, 12), b = rnd(rng, -12, 12), c = rnd(rng, -60, 60);
        const Diophantine d = linear_diophantine(a, b, c);
        if (a == 0 && b == 0) {
            MEQ_CHECK(d.ok == (c == 0) && d.dx == 0 && d.dy == 0);
        } else if (b == 0) {
            MEQ_CHECK(d.ok == (c % a == 0));
            if (d.ok) MEQ_CHECK(d.x0 == c / a && d.dx == 0 && d.dy == -1 && d.y(5) == 5 && d.x(7) == c / a);
        } else {
            std::vector<long long> xs;
            for (long long x = -150; x <= 150; ++x)
                if ((c - a * x) % b == 0) xs.push_back(x);
            MEQ_CHECK(d.ok == !xs.empty());
            if (d.ok) {
                MEQ_CHECK(0 <= d.x0 && d.x0 < d.dx && d.dx == std::llabs(b) / (long long)bf_gcd(uab(a), uab(b)));
                std::vector<long long> got;
                for (long long k = -170; k <= 170; ++k) {  // x0 < 12、dx ≥ 1 なので窓 [−150, 150] を覆う
                    const long long x = d.x(k);
                    if (-150 <= x && x <= 150) {
                        got.push_back(x);
                        MEQ_CHECK(d.y(k) == (c - a * x) / b);
                    }
                }
                MEQ_CHECK(got == xs);
            }
        }
        // 範囲の中の解の数・x が最小の解
        const long long xl = rnd(rng, -25, 25), xr = rnd(rng, -25, 25), yl = rnd(rng, -25, 25), yr = rnd(rng, -25, 25);
        long long cnt = 0, cnt_closed = 0;
        for (long long x = xl; x <= xr; ++x)
            for (long long y = yl; y <= yr; ++y)
                if (a * x + b * y == c) cnt_closed++, cnt += (x < xr && y < yr);
        MEQ_CHECK(diophantine_count(a, b, c, xl, xr, yl, yr) == cnt && diophantine_count_closed(a, b, c, xl, xr, yl, yr) == cnt_closed);
        const DiophantinePoint p = diophantine_min_x(a, b, c, xl);
        long long bx = LMAX;
        for (long long x = xl; x < xl + 200 && bx == LMAX; ++x)
            if (b == 0 ? a * x == c : (c - a * x) % b == 0) bx = x;
        MEQ_CHECK(p.ok == (bx != LMAX));
        if (p.ok) MEQ_CHECK(p.x == bx && a * p.x + b * p.y == c && (b != 0 || p.y == 0));
    }
    // 大きい値：128 ビットで検算
    for (int it = 0; it < 150000 * TEST_SCALE; ++it, ++cases) {
        long long a = no_lmin(rnd_ll(rng)), b = no_lmin(rnd_ll(rng)), c = rnd_ll(rng);
        if (it % 3 == 0) {  // 解を持つように作る
            const long long x = rnd(rng, -(1LL << 26), 1LL << 26), y = rnd(rng, -(1LL << 26), 1LL << 26);  // |a x| < 2^61
            a = rnd(rng, -(1LL << 31), 1LL << 31) * (it % 2 ? 1 : 12), b = rnd(rng, -(1LL << 31), 1LL << 31) * (it % 2 ? 1 : 18);
            c = a * x + b * y;
        }
        const unsigned long long g = bf_gcd(uab(a), uab(b));
        const bool solvable = g == 0 ? c == 0 : (Int128)c % (Int128)g == 0;
        // y0 が long long に入るときだけ linear_diophantine を呼ぶ（入らなければ止まる仕様）
        const bool fits = b == 0 ? !(a == -1 && c == LMIN) : (Int128)uab(c) / uab(b) + (Int128)uab(a) / g < (Int128)LMAX;
        if (fits) {
            const Diophantine d = linear_diophantine(a, b, c);
            MEQ_CHECK(d.ok == solvable);
            if (d.ok && b != 0) {
                MEQ_CHECK((Int128)a * d.x0 + (Int128)b * d.y0 == c && 0 <= d.x0 && d.x0 < d.dx && (unsigned long long)d.dx == uab(b) / g);
                MEQ_CHECK((Int128)d.dy == (Int128)(a / (long long)g) * (b < 0 ? -1 : 1));
                const long long k = rnd(rng, -3, 3);
                const Int128 X = (Int128)d.x0 + (Int128)k * d.dx, Y = (Int128)d.y0 - (Int128)k * d.dy;
                const Int128 L = (Int128)LMAX;
                if (-L <= X && X <= L && -L <= Y && Y <= L) MEQ_CHECK((Int128)a * X + (Int128)b * Y == c);  // 128 ビットに収まるときだけ
            }
        }
        // 範囲の中の解の数：x の範囲が狭いときに 1 つずつ調べる
        const long long xl = rnd_ll(rng), w = rnd(rng, 0, 100);
        const long long xr = xl > LMAX - w ? LMAX : xl + w;
        long long yl = rnd_ll(rng), yr = rnd_ll(rng);
        if (yl > yr) std::swap(yl, yr);
        if (it % 4 == 0) yl = LMIN, yr = LMAX;
        long long cnt = 0;
        long long first = LMAX;
        bool first_found = false;
        if (b == 0) continue;  // b = 0 は小さい値のテストで
        Int128 first_y = 0;
        for (long long x = xl; x < xr; ++x) {
            const Int128 rest = (Int128)c - (Int128)a * x;
            if (rest % b != 0) continue;
            const Int128 y = rest / b;
            if (!first_found) first = x, first_y = y, first_found = true;
            if (yl <= y && y < yr) cnt++;
        }
        MEQ_CHECK(diophantine_count(a, b, c, xl, xr, yl, yr) == cnt);
        if (first_found && uab(b) / g <= 100 && first_y <= (Int128)LMAX && first_y >= (Int128)LMIN) {  // 窓の中に最小の解がある（y が収まるとき）
            const DiophantinePoint p = diophantine_min_x(a, b, c, xl);
            MEQ_CHECK(p.ok && p.x == first && (Int128)p.y == first_y);
        }
    }
    // 端
    MEQ_CHECK(diophantine_count_closed(3, 5, 100, 0, LMAX, 0, LMAX) == 7);  // (0,20) (5,17) … (30,2)
    MEQ_CHECK(diophantine_count_closed(1, 1, LMAX - 1, 0, LMAX, 0, LMAX) == LMAX);
    MEQ_CHECK(diophantine_count(0, 0, 0, 0, 3, 0, 4) == 12 && diophantine_count(0, 0, 1, 0, 3, 0, 4) == 0);
    MEQ_CHECK(!diophantine_min_x(2, 4, 5, 0).ok && diophantine_min_x(1, 0, 7, 8).ok == false && diophantine_min_x(1, 0, 7, 7).x == 7);
    {
        const Diophantine d = linear_diophantine(LMAX, LMAX - 1, 1);  // x = 1, y = −1
        MEQ_CHECK(d.ok && d.x0 == 1 && d.y0 == -1 && d.dx == LMAX - 1 && d.dy == LMAX);
    }
}

static void test_crt(std::mt19937_64& rng, int& cases) {
    // 小さい値：lcm まで全部調べる
    for (int it = 0; it < 40000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 4);
        std::vector<long long> rs(n), ms(n);
        long long L = 1;
        for (int i = 0; i < n; ++i) {
            ms[i] = rnd(rng, 1, 30), rs[i] = it % 2 ? rnd(rng, -100, 100) : rnd_ll(rng);
            L = L / std::gcd(L, ms[i]) * ms[i];
        }
        if (L > 3000) continue;  // 全部調べるのは lcm が小さいときだけ
        long long best = -1;
        for (long long x = 0; x < L && best < 0; ++x) {
            bool ok = true;
            for (int i = 0; i < n && ok; ++i) ok = mod128(x - (Int128)rs[i], ms[i]) == 0;
            if (ok) best = x;
        }
        const ModSolution s = crt_all(rs, ms);
        MEQ_CHECK(!s.overflow && s.ok == (best >= 0));
        if (s.ok) MEQ_CHECK(s.r == best && s.m == L);
        if (n == 2) {
            const ModSolution t = crt_merge(rs[0], ms[0], rs[1], ms[1]);
            MEQ_CHECK(t.ok == s.ok && t.r == s.r && t.m == s.m && !t.overflow);
        }
    }
    // 大きい値：答え x を先に決めて作る
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 6);
        std::vector<long long> rs(n), ms(n);
        const Int128 x = it % 3 == 0 ? (Int128)rnd(rng, 0, LMAX) : (Int128)rnd(rng, 0, 1LL << 40);
        Int128 L = 1;
        bool big = false;
        for (int i = 0; i < n; ++i) {
            const int t = (int)rnd(rng, 0, 3);
            ms[i] = t == 0 ? rnd(rng, 1, LMAX) : (t == 1 ? rnd(rng, 1, 1000000000) : (t == 2 ? rnd(rng, 1, 64) : (1LL << rnd(rng, 0, 62))));
            rs[i] = (long long)mod128(x, ms[i]);  // 余りに mod の倍数を足してもよい
            if (it % 4 == 0 && rs[i] <= LMAX - ms[i]) rs[i] += ms[i];
            if (it % 4 == 1) rs[i] -= ms[i];
            if (!big) {
                const unsigned long long g = bf_gcd((unsigned long long)(L % ms[i]), (unsigned long long)ms[i]);
                L = L / g * ms[i];
                if (L > (Int128)LMAX) big = true;
            }
        }
        const ModSolution s = crt_all(rs, ms);
        MEQ_CHECK(s.overflow == big);
        if (!big) {
            MEQ_CHECK(s.ok && s.m == (long long)L && (Int128)s.r == mod128(x, L));
        } else {  // lcm > 2^63 > x なので x が最小の解
            MEQ_CHECK(s.ok && (Int128)s.r == x && s.m == 0);
        }
        // 1 つの余りをずらす → 解がないか、あっても lcm が収まらないときは最小の解が変わる
        if (n >= 2) {
            std::vector<long long> rs2 = rs;
            const int i = (int)rnd(rng, 0, n - 1);
            if (ms[i] == 1) continue;
            rs2[i] = (long long)mod128((Int128)rs2[i] + 1, ms[i]);
            const ModSolution u = crt_all(rs2, ms);
            // 1 つずつ確かめる：u.ok なら u.r はすべての式を満たす
            if (u.ok) {
                for (int j = 0; j < n; ++j) MEQ_CHECK(mod128((Int128)u.r - rs2[j], ms[j]) == 0);
                if (!u.overflow) MEQ_CHECK(u.r < u.m && u.m == s.m);
            }
        }
    }
    // lcm が収まらず、最小の解も long long を超える
    {
        const long long p = 3037000507LL, q = 3037000537LL, r2 = 4611686014132420609LL;  // p, q は素数（p·q > 2^63）
        const Int128 x = (Int128)LMAX + 12345;
        const ModSolution s = crt_merge((long long)(x % p), p, (long long)(x % q), q);
        MEQ_CHECK(s.overflow && !s.ok);
        const ModSolution t = crt_merge(5, p, 5, q);
        MEQ_CHECK(t.overflow && t.ok && t.r == 5);
        const ModSolution u = crt_all({5, 5, 6}, {p, q, 7});  // 5 ≢ 6 (mod 7) なので long long に入る解はない
        MEQ_CHECK(u.overflow && !u.ok);
        MEQ_CHECK(crt_all({0}, {r2}).ok && crt_all({}, {}).ok && crt_all({}, {}).m == 1);
    }
}

static void test_garner(std::mt19937_64& rng, int& cases) {
    const std::vector<long long> primes = {998244353, 1000000007, 1000000009, 167772161, 469762049, 754974721, 2147483647, 1000000000000000003LL, 999999999999999989LL};
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        std::vector<long long> ms;
        Int128 P = 1;
        // 互いに素な mod：素数から選ぶか、小さい数を互いに素になるように選ぶ
        if (it % 2 == 0) {
            std::vector<long long> pool = primes;
            std::shuffle(pool.begin(), pool.end(), rng);
            for (long long p : pool) {
                if (P > ((Int128)1 << 125) / p) break;
                if (ms.size() >= 4 || rng() % 4 == 0) break;
                ms.push_back(p), P *= p;
            }
        } else {
            for (int k = 0; k < 8; ++k) {
                const long long m = rnd(rng, 1, 1000);
                bool ok = true;
                for (long long q : ms) ok &= std::gcd(q, m) == 1;
                if (ok && P <= ((Int128)1 << 110) / m) ms.push_back(m), P *= m;
            }
        }
        const UInt128 ux = ((UInt128)rng() << 64) | (UInt128)rng();
        const Int128 x = (Int128)(ux % (UInt128)P);
        std::vector<long long> rs;
        for (long long m : ms) rs.push_back((long long)(x % m) - (it % 3 == 0 ? m : 0));
        const long long mod = it % 3 == 0 ? 998244353 : rnd(rng, 1, LMAX);
        MEQ_CHECK(garner(rs, ms, mod) == (long long)(x % mod));
    }
    MEQ_CHECK(garner({}, {}, 7) == 0 && garner({3}, {1}, 10) == 0);
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_linear_congruence(rng, cases);
    test_diophantine(rng, cases);
    test_crt(rng, cases);
    test_garner(rng, cases);
    std::cout << "[mod_equation] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[mod_equation] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 1000000;
    std::vector<long long> xs(n), ys(n);
    for (int i = 0; i < n; ++i) xs[i] = rnd(rng, 1, 1000000000000000000LL), ys[i] = rnd(rng, 1, 1000000000000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (int i = 0; i < n; ++i) {
        const ModSolution s = linear_congruence(xs[i], ys[i], ys[(i + 1) % n]);
        s1 += s.ok ? s.r % 7 : 1;
    }
    const long long t_lc = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for (int i = 0; i < n; ++i) {
        const Diophantine d = linear_diophantine(xs[i] % 1000000000, ys[i] % 1000000000, xs[(i + 1) % n]);
        s2 += d.ok ? d.x0 % 7 : 1;
    }
    const long long t_dio = ms(t0);
    // crt_all：10^6 個の式（mod は 720720 の約数。答えはすべての式を満たす）
    std::vector<long long> divs;
    for (long long d = 1; d <= 720720; ++d)
        if (720720 % d == 0) divs.push_back(d);
    std::vector<long long> rs(n), mods(n);
    const long long x = 123456;
    for (int i = 0; i < n; ++i) mods[i] = divs[rng() % divs.size()], rs[i] = x % mods[i];
    t0 = std::chrono::steady_clock::now();
    const ModSolution c = crt_all(rs, mods);
    long long s3 = c.r;
    // 互いに素な 2 つ（素数の積は 10^18 程度）× 10^6
    for (int i = 0; i < n; ++i) s3 += crt_merge(xs[i], 1000000007, ys[i], 999999937).r % 7;
    const long long t_crt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s4 = 0;
    const std::vector<long long> gm = {998244353, 1000000007, 1000000009, 167772161, 469762049};
    for (int i = 0; i < n / 5; ++i) s4 += garner({xs[i] % gm[0], xs[i] % gm[1], ys[i] % gm[2], ys[i] % gm[3], xs[i] % gm[4]}, gm, 998244353);
    const long long t_gar = ms(t0);
    std::cout << "[mod_equation] benchmark: linear_congruence × 10^6 " << t_lc << " ms, linear_diophantine × 10^6 " << t_dio << " ms, crt_all（10^6 個）+ crt_merge × 10^6 "
              << t_crt << " ms, garner（5 個）× 2×10^5 " << t_gar << " ms  (check " << (s1 + s2 + s3 + s4) % 7 << ")\n";
}

}  // namespace meq_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace meq_demo {

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
#define SHOW(expr, comment) meq_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    meq_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

// 使い方関数：周期的に繰り返す 2 つの区間が重なる最初の時刻（ABC193 E の形）
//   t mod T1 ∈ [l1, r1) かつ t mod T2 ∈ [l2, r2) となる最小の t ≥ 0（なければ −1）。(r1 − l1)(r2 − l2) 回の crt_merge
static long long first_common_time(long long T1, long long l1, long long r1, long long T2, long long l2, long long r2) {
    long long best = -1;
    for (long long a = l1; a < r1; ++a)
        for (long long b = l2; b < r2; ++b) {
            const ModSolution s = crt_merge(a, T1, b, T2);
            if (s.ok && (best < 0 || s.r < best)) best = s.r;
        }
    return best;
}

static void demo() {
    banner("合同式・一次不定方程式・中国剰余定理（mod_equation） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  linear_congruence(a, b, m)                    a x ≡ b (mod m) → ModSolution{ok, r, m}（x ≡ r (mod m)）\n"
              << "  linear_diophantine(a, b, c)                   a x + b y = c → Diophantine{ok, x0, y0, dx, dy}\n"
              << "  diophantine_count / _closed / diophantine_min_x   範囲の中の解の数 / x ≥ x_min で最小の解\n"
              << "  crt_merge(r1, m1, r2, m2) / crt_all(rs, ms)   中国剰余定理（互いに素でなくてよい）\n"
              << "  garner(rs, ms, mod)                           互いに素な mod での答えを別の mod で\n";

    api("【1】 linear_congruence(a, b, m)", "ModSolution{ok, r, m}",
        "a x ≡ b (mod m) の解。x ≡ r (mod m')（m' = m / gcd(a, m)、0 ≤ r < m'）。解がなければ ok = false", "O(log m)");
    RUN(ModSolution s = linear_congruence(3, 6, 10));
    SHOW(s.ok, "解がある");
    SHOW(s.r, "x ≡ 2 (mod 10)");
    SHOW(linear_congruence(4, 6, 10).r, "4x ≡ 6 (mod 10) → x ≡ 4 (mod 5)");
    SHOW(linear_congruence(4, 6, 10).m, "解の mod は 10 / gcd(4, 10) = 5");
    SHOW(linear_congruence(2, 1, 4).ok, "2x ≡ 1 (mod 4) は解なし");
    {
        // ABC186 E（Throne）：S + K x ≡ 0 (mod N) となる最小の x ≥ 0
        note("ABC186 E（Throne）: K x ≡ −S (mod N) を解く。r が答え（解がなければ −1）");
        const long long cases[4][3] = {{10, 4, 3}, {1000, 11, 2}, {998244353, 897581057, 595591169}, {10000, 6, 14}};
        for (auto& t : cases) {
            const ModSolution r = linear_congruence(t[2], -t[1], t[0]);
            std::cout << "    N = " << t[0] << ", S = " << t[1] << ", K = " << t[2] << " → " << (r.ok ? r.r : -1) << "\n";
        }
    }

    api("【2】 linear_diophantine(a, b, c)", "Diophantine{ok, x0, y0, dx, dy}",
        "a x + b y = c の整数解すべて：x = x0 + k·dx、y = y0 − k·dy（d.x(k) / d.y(k)）。0 ≤ x0 < dx（x0 は 0 以上で最小の x）", "O(log)");
    RUN(Diophantine d = linear_diophantine(3, 5, 100));
    SHOW(d.x0, "x = 0 + 5k");
    SHOW(d.y0, "y = 20 − 3k");
    SHOW(d.x(2), "k = 2 なら x = 10");
    SHOW(d.y(2), "y = 14（3·10 + 5·14 = 100）");
    SHOW(linear_diophantine(4, 6, 7).ok, "4x + 6y = 7 は解なし（gcd = 2 が 7 を割らない）");

    api("【3】 diophantine_count(a, b, c, xl, xr, yl, yr) / diophantine_count_closed(…) / diophantine_min_x(a, b, c, x_min)", "long long / DiophantinePoint{ok, x, y}",
        "x, y の範囲の中の解の数（[xl, xr) / 閉区間）。x ≥ x_min で x が最小の解", "O(log)");
    SHOW(diophantine_count_closed(3, 5, 100, 0, LLONG_MAX, 0, LLONG_MAX), "x, y ≥ 0 の解：(0,20) (5,17) … (30,2)");
    SHOW(diophantine_count(1, 1, 10, 0, 5, 0, 100), "x ∈ [0, 5) で x + y = 10");
    RUN(DiophantinePoint p = diophantine_min_x(3, 5, 100, 1));
    SHOW(p.x, "x ≥ 1 で最小の x");
    SHOW(p.y, "そのときの y");
    note("例：a 円玉と b 円玉でちょうど c 円にする組み合わせの数 = diophantine_count_closed(a, b, c, 0, LLONG_MAX, 0, LLONG_MAX)");

    api("【4】 crt_merge(r1, m1, r2, m2) / crt_all(rs, ms)", "ModSolution{ok, r, m, overflow}",
        "x ≡ r_i (mod m_i) をすべて満たす x ≡ r (mod lcm)。r は 0 以上で最小の解。mod は互いに素でなくてよい", "O(n log)");
    SHOW(crt_merge(2, 3, 3, 5).r, "x ≡ 2 (mod 3)、x ≡ 3 (mod 5) → 8");
    SHOW(crt_merge(2, 3, 3, 5).m, "lcm = 15");
    SHOW(crt_all({1, 3, 5}, {4, 6, 8}).ok, "互いに素でなくても、公約数の分で矛盾しなければ解がある");
    SHOW(crt_all({1, 3, 5}, {4, 6, 8}).r, "21 ≡ 1 (mod 4)、≡ 3 (mod 6)、≡ 5 (mod 8)");
    SHOW(crt_all({0, 1}, {4, 6}).ok, "0 と 1 は mod 2 で違うので解なし");
    RUN(ModSolution big = crt_merge(5, 3037000507LL, 5, 3037000537LL));
    SHOW(big.overflow, "lcm が long long を超えた");
    SHOW(big.ok, "それでも long long に入る最小の解はわかる");
    SHOW(big.r, "5");
    note("ABC193 E（Train）：使い方関数 first_common_time（下の例）。t mod (2X+2Y) ∈ [X, X+Y) かつ t mod (P+Q) ∈ [P, P+Q)");
    {
        const long long cases[3][4] = {{5, 2, 7, 6}, {1, 1, 3, 1}, {999999999, 1, 1000000000, 1}};
        for (auto& c : cases) {
            const long long X = c[0], Y = c[1], P = c[2], Q = c[3];
            const long long t = first_common_time(2 * X + 2 * Y, X, X + Y, P + Q, P, P + Q);
            std::cout << "    X = " << X << ", Y = " << Y << ", P = " << P << ", Q = " << Q << " → " << (t < 0 ? std::string("infinity") : std::to_string(t)) << "\n";
        }
    }

    api("【5】 garner(rs, ms, mod)", "long long",
        "互いに素な ms での余りから、0 ≤ x < Π ms の x（Π ms は long long を超えてよい）を mod で", "O(n^2)");
    SHOW(garner({2, 3, 2}, {3, 5, 7}, 1000000007), "x = 23（2 mod 3、3 mod 5、2 mod 7）");
    RUN(const long long X = 123456789012345678LL);
    SHOW(garner({X % 998244353, X % 1000000007, X % 1000000009}, {998244353, 1000000007, 1000000009}, 1000000000000LL), "X mod 10^12 = 789012345678");
    note("いくつかの素数 mod で計算した答えから本当の値（の mod）を戻すとき（任意 mod の畳み込みなど）に使う");
}

}  // namespace meq_demo

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
    meq_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    meq_test::benchmark();    // 速度計測
    meq_demo::demo();         // 使い方デモ
    return 0;
}
