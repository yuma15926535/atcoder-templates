// =============================================================================
//  素数の個数・和（N が大きい π(N)・素数の和、⌊N/i⌋ ごとの π の表、素因数がすべて P 以下の数の個数、乗法的関数の和（Min_25））   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 素数・約数
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  機能一覧                                                                    計算量
//    prime_count(N)                   N 以下の素数の個数（N = 10^11 で約 0.3 秒、10^12 で約 1.3 秒）       O(N^{3/4} / log N)
//    prime_sum(N)                     N 以下の素数の和（Int128）
//    PrimeCountTable P(N)             ⌊N/i⌋ の形のすべての v での π(v)：P.count(v) / P(v)（v ≤ √N でもよい）
//    count_smooth(N, P)               N 以下で素因数がすべて P 以下の数の個数（1 も数える。P ≤ 100 程度。(10^16, 97) で約 0.8 秒。
//                                     素数を 2 組に分けて作った数を合わせる。作る数が 4×10^7 個を超えると止める）
//    multiplicative_sum<T>(N, coef, fpe)   Σ_{i≤N} f(i)（f は乗法的。Min_25 篩）。素数で f(p) = Σ coef[k] p^k（3 次まで）、
//                                     素数の累乗で f(p^e) = fpe(p, e)。T は mint・long long・Int128（整数は溢れたら止める）
//                                     例 d：{2}, e+1 ／ φ：{−1, 1}, p^e − p^(e−1) ／ μ：{−1}, (e == 1 ? −1 : 0)
//                                     N ≤ 10^12（mint で 10^10 なら約 0.2 秒、10^11 で約 1 秒）                  O(N^{3/4} / log N) 程度
//  （N ≤ 10^8 なら 02_sieve の Sieve / primes_up_to の方が速く、区間の素数は segment_primes）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    pcnt_test::random_test()  … 篩で数えた値・よく知られた値との突き合わせ
//    pcnt_test::benchmark()    … 速度計測
//    pcnt_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
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
#ifndef GG_PRIME_COUNT
#define GG_PRIME_COUNT
// ---- 素数の個数・和（数学第二世代）  N が大きい（10^11〜10^12）π(N)・素数の和、⌊N/i⌋ ごとの表、P-smooth な数の個数、乗法的関数の和（Min_25） ----
// ⌊N/i⌋ の形の値 v ごとの π(v)（Lucy の方法）。P.count(v) は v ≤ √N か v = ⌊N/i⌋ のとき  O(N^{3/4} / log N)
struct PrimeCountTable {
    explicit PrimeCountTable(long long N) : n_(N) {
        if (N < 0) math_internal::fail("PrimeCountTable", "N ≥ 0 で使う");
        r_ = isqrt(N);
        lo_.assign((size_t)r_ + 1, 0), hi_.assign((size_t)r_ + 1, 0);
        for (long long v = 1; v <= r_; ++v) lo_[(size_t)v] = v - 1;    // lo[v] = π(v) の候補（はじめは 2〜v の個数）
        for (long long i = 1; i <= r_; ++i) hi_[(size_t)i] = N / i - 1;  // hi[i] = π(N/i) の候補
        for (long long p = 2; p <= r_; ++p) {
            if (lo_[(size_t)p] == lo_[(size_t)p - 1]) continue;  // p は合成数
            const long long pc = lo_[(size_t)p - 1], p2 = p * p, lim = std::min(r_, N / p2);
            for (long long i = 1; i <= lim; ++i) {  // N/i ≥ p^2 のもの：p を最小の素因数に持つ数を除く
                const long long d = i * p;
                hi_[(size_t)i] -= (d <= r_ ? hi_[(size_t)d] : lo_[(size_t)(N / d)]) - pc;
            }
            for (long long v = r_; v >= p2; --v) lo_[(size_t)v] -= lo_[(size_t)(v / p)] - pc;
        }
    }
    long long count(long long v) const {  // π(v)（v ≤ √N か v = ⌊N/i⌋）
        if (v < 0) return 0;
        if (v <= r_) return lo_[(size_t)v];
        if (v > n_ || n_ / (n_ / v) != v) math_internal::fail("PrimeCountTable::count", "v は √N 以下か ⌊N/i⌋ の形");
        return hi_[(size_t)(n_ / v)];
    }
    long long operator()(long long v) const { return count(v); }
    long long n() const { return n_; }

  private:
    long long n_, r_;
    std::vector<long long> lo_, hi_;
};
inline long long prime_count(long long N) { return N < 2 ? 0 : PrimeCountTable(N).count(N); }  // N 以下の素数の個数

// N 以下の素数の和（Int128。N ≤ 10^12 程度）  O(N^{3/4} / log N)
inline Int128 prime_sum(long long N) {
    if (N < 2) return 0;
    const long long r = isqrt(N);
    std::vector<Int128> lo((size_t)r + 1, 0), hi((size_t)r + 1, 0);
    auto tri = [](long long v) { return (Int128)v * (v + 1) / 2 - 1; };  // 2 + 3 + … + v
    for (long long v = 1; v <= r; ++v) lo[(size_t)v] = tri(v);
    for (long long i = 1; i <= r; ++i) hi[(size_t)i] = tri(N / i);
    for (long long p = 2; p <= r; ++p) {
        if (lo[(size_t)p] == lo[(size_t)p - 1]) continue;
        const Int128 pc = lo[(size_t)p - 1];
        const long long p2 = p * p, lim = std::min(r, N / p2);
        for (long long i = 1; i <= lim; ++i) {
            const long long d = i * p;
            hi[(size_t)i] -= (Int128)p * ((d <= r ? hi[(size_t)d] : lo[(size_t)(N / d)]) - pc);
        }
        for (long long v = r; v >= p2; --v) lo[(size_t)v] -= (Int128)p * (lo[(size_t)(v / p)] - pc);
    }
    return hi[1];
}

// N 以下の正の整数で、素因数がすべて P 以下のものの個数（1 も数える。P ≤ 100 程度・N ≤ 10^18。素数を 2 組に分けて数え合わせる）
inline long long count_smooth(long long N, long long P) {
    if (N < 1) return 0;
    if (P > 1000) math_internal::fail("count_smooth", "P ≤ 1000 で使う（N が大きいときは P ≤ 100 程度）");
    std::vector<long long> ps;
    for (long long p = 2; p <= P; ++p) {
        bool pr = true;
        for (long long q = 2; q * q <= p; ++q) pr &= p % q != 0;
        if (pr) ps.push_back(p);
    }
    std::vector<long long> g1, g2;  // 素数を交互に 2 組に
    for (size_t i = 0; i < ps.size(); ++i) (i % 2 ? g2 : g1).push_back(ps[i]);
    size_t total = 0;
    auto enumerate = [&](const std::vector<long long>& g) {  // g の素数だけでできた N 以下の数（1 も）
        std::vector<long long> v = {1};
        for (long long p : g) {
            const size_t m = v.size();
            for (size_t i = 0; i < m; ++i)
                for (long long x = v[i]; x <= N / p;) {
                    x *= p, v.push_back(x);
                    if (++total > 40000000) math_internal::fail("count_smooth", "途中で作る数が多すぎる（N か P が大きすぎる。(10^18, 97) くらいまで）");
                }
        }
        std::sort(v.begin(), v.end());
        return v;
    };
    const std::vector<long long> A = enumerate(g1), B = enumerate(g2);
    long long cnt = 0;
    long long j = (long long)B.size() - 1;
    for (long long a : A) {  // a·b ≤ N となる b の数（a が増えると b の上限は減る）
        while (j >= 0 && B[(size_t)j] > N / a) --j;
        if (j < 0) break;
        cnt += j + 1;
    }
    return cnt;
}

// ---- 乗法的関数の和（Min_25 篩）----
namespace math_internal {
template <class T>
constexpr bool pc_is_int = std::is_integral<T>::value || std::is_same<T, Int128>::value;
template <class T>
inline T pc_add(const T& a, const T& b) {
    if constexpr (pc_is_int<T>) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail("multiplicative_sum", "値が整数の型に収まらない（mint か Int128 で）");
        return r;
    } else {
        return a + b;
    }
}
template <class T>
inline T pc_sub(const T& a, const T& b) {
    if constexpr (pc_is_int<T>) {
        T r;
        if (__builtin_sub_overflow(a, b, &r)) fail("multiplicative_sum", "値が整数の型に収まらない（mint か Int128 で）");
        return r;
    } else {
        return a - b;
    }
}
template <class T>
inline T pc_mul(const T& a, const T& b) {
    if constexpr (pc_is_int<T>) {
        T r;
        if (__builtin_mul_overflow(a, b, &r)) fail("multiplicative_sum", "値が整数の型に収まらない（mint か Int128 で）");
        return r;
    } else {
        return a * b;
    }
}
template <class T>
T pc_power_sum(long long v, int k) {  // 1^k + 2^k + … + v^k（k ≤ 3。割り算は先に整数で済ませる）
    if (k == 0) return T(v);
    long long a = v, b = v + 1;
    (a % 2 == 0 ? a : b) /= 2;
    if (k == 1) return pc_mul(T(a), T(b));
    if (k == 3) {
        const T s = pc_mul(T(a), T(b));
        return pc_mul(s, s);
    }
    long long x = v, y = v + 1, z = 2 * v + 1;  // k = 2：v(v+1)(2v+1)/6
    (x % 2 == 0 ? x : y) /= 2;
    (x % 3 == 0 ? x : (y % 3 == 0 ? y : z)) /= 3;
    return pc_mul(pc_mul(T(x), T(y)), T(z));
}
}  // namespace math_internal
// Σ_{i=1}^{N} f(i)：f は乗法的（f(1) = 1）で、素数 p では f(p) = coef[0] + coef[1] p + coef[2] p^2 + coef[3] p^3（3 次まで）、
// 素数の累乗では f(p^e) = fpe(p, e)（e ≥ 1。fpe(p, 1) は f(p) と同じ値にする）。T は mint・long long・Int128（整数は溢れたら止める）
//   例：約数の個数 d … coef {2}、fpe = e + 1 ／ φ … coef {−1, 1}、fpe = p^e − p^(e−1) ／ μ … coef {−1}、fpe = (e == 1 ? −1 : 0)
//   N ≤ 10^12（mint で 10^10 なら約 0.2 秒、10^11 で約 1 秒、10^12 で約 4 秒）  O(N^{3/4} / log N) 程度
template <class T, class Fpe>
T multiplicative_sum(long long N, const std::vector<T>& coef, Fpe fpe) {
    using namespace math_internal;
    if (N < 1) return T(0);
    if (coef.empty() || coef.size() > 4) fail("multiplicative_sum", "coef は 1〜4 個（f(p) は p の 3 次まで）");
    if (N > 1000000000000LL) fail("multiplicative_sum", "N ≤ 10^12 で使う");
    const long long r = isqrt(N);
    std::vector<long long> ps;  // √N 以下の素数
    {
        std::vector<char> comp((size_t)r + 1, 0);
        for (long long p = 2; p <= r; ++p)
            if (!comp[(size_t)p]) {
                ps.push_back(p);
                for (long long q = p * p; q <= r; q += p) comp[(size_t)q] = 1;
            }
    }
    // ⌊N/i⌋ の形の値 v での Σ_{素数 q ≤ v} f(q)：v ≤ r は Glo[v]、v > r は Ghi[N/v]（Lucy の方法を p^k ごとに）
    std::vector<T> Glo((size_t)r + 1, T(0)), Ghi((size_t)r + 1, T(0));
    for (int k = 0; k < (int)coef.size(); ++k) {
        if (coef[(size_t)k] == T(0)) continue;
        std::vector<T> lo((size_t)r + 1), hi((size_t)r + 1);
        lo[0] = T(0);
        for (long long v = 1; v <= r; ++v) lo[(size_t)v] = pc_sub(pc_power_sum<T>(v, k), T(1));
        for (long long i = 1; i <= r; ++i) hi[(size_t)i] = pc_sub(pc_power_sum<T>(N / i, k), T(1));
        for (const long long p : ps) {
            T pw = T(1);
            for (int t = 0; t < k; ++t) pw = pc_mul(pw, T(p));
            const T pc = lo[(size_t)p - 1];  // p より小さい素数の q^k の和
            const long long p2 = p * p, lim = std::min(r, N / p2);
            if (k == 0) {  // p^0 = 1 なので掛け算なし
                for (long long i = 1; i <= lim; ++i) {
                    const long long d = i * p;
                    hi[(size_t)i] = pc_sub(hi[(size_t)i], pc_sub(d <= r ? hi[(size_t)d] : lo[(size_t)(N / d)], pc));
                }
                for (long long v = r; v >= p2; --v) lo[(size_t)v] = pc_sub(lo[(size_t)v], pc_sub(lo[(size_t)(v / p)], pc));
            } else {
                for (long long i = 1; i <= lim; ++i) {
                    const long long d = i * p;
                    hi[(size_t)i] = pc_sub(hi[(size_t)i], pc_mul(pw, pc_sub(d <= r ? hi[(size_t)d] : lo[(size_t)(N / d)], pc)));
                }
                for (long long v = r; v >= p2; --v) lo[(size_t)v] = pc_sub(lo[(size_t)v], pc_mul(pw, pc_sub(lo[(size_t)(v / p)], pc)));
            }
        }
        for (long long v = 1; v <= r; ++v) Glo[(size_t)v] = pc_add(Glo[(size_t)v], pc_mul(coef[(size_t)k], lo[(size_t)v]));
        for (long long i = 1; i <= r; ++i) Ghi[(size_t)i] = pc_add(Ghi[(size_t)i], pc_mul(coef[(size_t)k], hi[(size_t)i]));
    }
    const int P = (int)ps.size();
    // S(v, j) = Σ_{2 ≤ x ≤ v、x の最小の素因数 ≥ ps[j]} f(x)
    auto S = [&](auto&& self, long long v, int j) -> T {
        const long long last = j == 0 ? 1 : ps[(size_t)j - 1];
        if (v <= last) return T(0);
        T res = pc_sub(v <= r ? Glo[(size_t)v] : Ghi[(size_t)(N / v)], Glo[(size_t)last]);  // ps[j] 以上の素数
        for (int k = j; k < P && ps[(size_t)k] * ps[(size_t)k] <= v; ++k) {
            const long long p = ps[(size_t)k];
            long long pe = p;
            T fe = T(fpe(p, 1));
            for (int e = 1; pe <= v / p; ++e, pe *= p) {
                const T fe1 = T(fpe(p, e + 1));
                res = pc_add(res, pc_add(pc_mul(fe, self(self, v / pe, k + 1)), fe1));
                fe = fe1;
            }
        }
        return res;
    };
    return pc_add(S(S, N, 0), T(1));
}
#endif  // GG_PRIME_COUNT
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace pcnt_test {

#define PCNT_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[prime_count] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

struct Mod7 {  // テスト用の小さな mod 型（mint の代わり。10^9+7）
    static constexpr unsigned long long MD = 1000000007;
    long long v;
    Mod7(long long x = 0) : v(x % (long long)MD) {
        if (v < 0) v += MD;
    }
    friend Mod7 operator+(Mod7 a, Mod7 b) {
        a.v += b.v;
        if (a.v >= (long long)MD) a.v -= MD;
        return a;
    }
    friend Mod7 operator-(Mod7 a, Mod7 b) {
        a.v -= b.v;
        if (a.v < 0) a.v += MD;
        return a;
    }
    friend Mod7 operator*(Mod7 a, Mod7 b) {
        a.v = (long long)((unsigned long long)a.v * (unsigned long long)b.v % MD);
        return a;
    }
    friend bool operator==(Mod7 a, Mod7 b) { return a.v == b.v; }
};
static long long ipow(long long p, int e) {
    long long r = 1;
    while (e-- > 0) r *= p;
    return r;
}
static void test_multiplicative(std::mt19937_64& rng, int& cases) {
    // 小さい N：最小の素因数の表から各関数を作って累積和と比べる
    const int M = 200000;
    std::vector<int> spf((size_t)M + 1, 0);
    for (int i = 2; i <= M; ++i)
        if (!spf[(size_t)i])
            for (int j = i; j <= M; j += i)
                if (!spf[(size_t)j]) spf[(size_t)j] = i;
    // 0: d, 1: σ1, 2: φ, 3: μ, 4: μ², 5: σ2, 6: i^3
    const int F = 7;
    std::vector<std::vector<Int128>> pre(F, std::vector<Int128>((size_t)M + 1, 0));
    for (int x = 1; x <= M; ++x) {
        Int128 v[F] = {1, 1, 1, 1, 1, 1, 1};
        for (int y = x; y > 1;) {
            const int p = spf[(size_t)y];
            int e = 0;
            while (y % p == 0) y /= p, ++e;
            const Int128 pe = ipow(p, e);
            Int128 s1 = 0, s2 = 0;
            for (int t = 0; t <= e; ++t) s1 += ipow(p, t), s2 += (Int128)ipow(p, t) * ipow(p, t);
            v[0] *= e + 1, v[1] *= s1, v[2] *= pe - pe / p, v[3] *= e == 1 ? -1 : 0, v[4] *= e == 1 ? 1 : 0, v[5] *= s2, v[6] *= pe * pe * pe;
        }
        for (int f = 0; f < F; ++f) pre[(size_t)f][(size_t)x] = pre[(size_t)f][(size_t)x - 1] + v[f];
    }
    auto sig = [](long long p, int e, int k) {  // 1 + p^k + … + p^{ke}
        Int128 s = 0, t = 1;
        for (int i = 0; i <= e; ++i) s += t, t *= ipow(p, k);
        return s;
    };
    auto calc = [&](int f, long long N) -> Int128 {
        switch (f) {
            case 0: return multiplicative_sum<Int128>(N, {2}, [](long long, int e) { return (Int128)e + 1; });
            case 1: return multiplicative_sum<Int128>(N, {1, 1}, [&](long long p, int e) { return sig(p, e, 1); });
            case 2: return multiplicative_sum<Int128>(N, {-1, 1}, [](long long p, int e) { return (Int128)ipow(p, e) - ipow(p, e - 1); });
            case 3: return multiplicative_sum<Int128>(N, {-1}, [](long long, int e) { return (Int128)(e == 1 ? -1 : 0); });
            case 4: return multiplicative_sum<Int128>(N, {1}, [](long long, int e) { return (Int128)(e == 1); });
            case 5: return multiplicative_sum<Int128>(N, {1, 0, 1}, [&](long long p, int e) { return sig(p, e, 2); });
            default: return multiplicative_sum<Int128>(N, {0, 0, 0, 1}, [](long long p, int e) { const Int128 q = ipow(p, e); return q * q * q; });
        }
    };
    for (int it = 0; it < 200 * TEST_SCALE; ++it, ++cases) {
        const long long N = it < 40 ? it : rnd(rng, 0, it % 3 ? 3000 : M);
        const int f = it % F;
        PCNT_CHECK(calc(f, N) == pre[(size_t)f][(size_t)N]);
        if (f == 0) PCNT_CHECK(multiplicative_sum<long long>(N, {2}, [](long long, int e) { return e + 1LL; }) == (long long)pre[0][(size_t)N]);
        if (f == 2) {  // mod 型でも
            const Mod7 m = multiplicative_sum<Mod7>(N, {Mod7(-1), Mod7(1)}, [](long long p, int e) { return Mod7(ipow(p, e) - ipow(p, e - 1)); });
            PCNT_CHECK(m.v == (long long)(pre[2][(size_t)N] % 1000000007));
        }
    }
    // 大きい N：別の方法で
    const long long NB = TEST_SCALE > 1 ? 10000000000LL : 100000000LL;
    {
        Int128 dsum = 0, s1 = 0;  // Σ d(i) = Σ ⌊N/i⌋、Σ σ1(i) = Σ i ⌊N/i⌋
        for (long long i = 1; i <= NB;) {
            const long long q = NB / i, j = NB / q;
            dsum += (Int128)q * (j - i + 1), s1 += (Int128)q * ((Int128)(i + j) * (j - i + 1) / 2);
            i = j + 1;
        }
        PCNT_CHECK(calc(0, NB) == dsum && calc(1, NB) == s1);
        const Int128 n = NB / 100;  // Σ i^3 = (n(n+1)/2)^2
        PCNT_CHECK(calc(6, NB / 100) == (n * (n + 1) / 2) * (n * (n + 1) / 2));
        cases += 3;
    }
    {  // Σφ・Σμ をメモ化再帰で（Φ(N) = N(N+1)/2 − Σ_{d≥2} Φ(N/d)、M(N) = 1 − Σ_{d≥2} M(N/d)）
        const long long N = TEST_SCALE > 1 ? 1000000000LL : 10000000LL;
        std::map<long long, std::pair<Int128, Int128>> memo;
        auto rec = [&](auto&& self, long long v) -> std::pair<Int128, Int128> {
            if (v <= M) return {pre[2][(size_t)v], pre[3][(size_t)v]};
            auto itr = memo.find(v);
            if (itr != memo.end()) return itr->second;
            Int128 a = (Int128)v * (v + 1) / 2, b = 1;
            for (long long d = 2; d <= v;) {
                const long long q = v / d, d2 = v / q;
                const auto s = self(self, q);
                a -= (Int128)(d2 - d + 1) * s.first, b -= (Int128)(d2 - d + 1) * s.second;
                d = d2 + 1;
            }
            return memo[v] = {a, b};
        };
        const auto ex = rec(rec, N);
        PCNT_CHECK(calc(2, N) == ex.first && calc(3, N) == ex.second);
        cases += 2;
    }
    PCNT_CHECK(calc(4, 100) == 61 && calc(3, 1) == 1 && calc(0, 0) == 0);  // 100 以下の平方因子なしは 61 個
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    const int M = 3000000;
    std::vector<char> comp((size_t)M + 1, 0);
    std::vector<long long> pc((size_t)M + 1, 0);
    std::vector<Int128> ps((size_t)M + 1, 0);
    for (int i = 2; i <= M; ++i) {
        if (!comp[(size_t)i])
            for (long long j = (long long)i * i; j <= M; j += i) comp[(size_t)j] = 1;
        pc[(size_t)i] = pc[(size_t)i - 1] + !comp[(size_t)i];
        ps[(size_t)i] = ps[(size_t)i - 1] + (comp[(size_t)i] ? 0 : i);
    }
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const long long N = it < 50 ? it : rnd(rng, 0, it % 3 ? 100000 : M);
        PCNT_CHECK(prime_count(N) == pc[(size_t)N] && prime_sum(N) == ps[(size_t)N]);
        if (it % 10 == 0) {  // 表のすべての値
            const PrimeCountTable T(N);
            for (long long i = 1; i <= N; i = N / (N / i) + 1) PCNT_CHECK(T.count(N / i) == pc[(size_t)(N / i)]);
            for (long long v = 0; v * v <= N; ++v) PCNT_CHECK(T(v) == pc[(size_t)v]);
        }
    }
    // よく知られた値
    PCNT_CHECK(prime_count(10000000000LL) == 455052511 && prime_count(100000000000LL) == 4118054813LL);
    PCNT_CHECK(prime_sum(2000000) == 142913828922LL);  // 200 万未満の素数の和
    // P-smooth：小さい N は全部調べる
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const long long N = rnd(rng, 0, 20000), P = rnd(rng, 1, 60);
        long long e = 0;
        for (long long x = 1; x <= N; ++x) {
            long long y = x;
            for (long long p = 2; p <= P && y > 1; ++p)
                while (y % p == 0) y /= p;
            e += y == 1;
        }
        PCNT_CHECK(count_smooth(N, P) == e);
    }
    PCNT_CHECK(count_smooth(36, 3) == 14 && count_smooth(1000000000000000000LL, 2) == 60 && count_smooth(LLONG_MAX, 1) == 1);
    test_multiplicative(rng, cases);
    std::cout << "[prime_count] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[prime_count] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    const long long c11 = prime_count(100000000000LL);
    const long long t11 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const long long c12 = prime_count(1000000000000LL);
    const long long t12 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Int128 s11 = prime_sum(100000000000LL);
    const long long ts = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const long long sm = count_smooth(10000000000000000LL, 97);
    const long long tsm = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Mod7 m10 = multiplicative_sum<Mod7>(10000000000LL, {Mod7(-1), Mod7(1)}, [](long long p, int e) { return Mod7(ipow(p, e) - ipow(p, e - 1)); });
    const long long tm10 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Mod7 m11 = multiplicative_sum<Mod7>(100000000000LL, {Mod7(-1), Mod7(1)}, [](long long p, int e) { return Mod7(ipow(p, e) - ipow(p, e - 1)); });
    const long long tm11 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Mod7 m12 = multiplicative_sum<Mod7>(1000000000000LL, {Mod7(2)}, [](long long, int e) { return Mod7(e + 1); });
    const long long tm12 = ms(t0);
    std::cout << "[prime_count] benchmark: prime_count(10^11) " << t11 << " ms, prime_count(10^12) " << t12 << " ms, prime_sum(10^11) " << ts
              << " ms, count_smooth(10^16, 97) " << tsm << " ms\n"
              << "[prime_count] benchmark: multiplicative_sum（Σφ mod p）10^10 " << tm10 << " ms, 10^11 " << tm11 << " ms, （Σd mod p）10^12 " << tm12
              << " ms  (check " << (long long)((c11 + c12 + s11 + sm + m10.v + m11.v + m12.v) % 7) << ")\n";
}

}  // namespace pcnt_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace pcnt_demo {

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
#define SHOW(expr, comment) pcnt_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    pcnt_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("素数の個数・和（prime_count） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  prime_count(N) / prime_sum(N)        N 以下の素数の個数 / 和（N ≤ 10^12 程度）\n"
              << "  PrimeCountTable P(N)：P.count(v)       ⌊N/i⌋ の形の v での π(v)\n"
              << "  count_smooth(N, P)                   素因数がすべて P 以下の数の個数\n"
              << "  multiplicative_sum<T>(N, coef, fpe)  乗法的関数の和 Σ_{i≤N} f(i)（Min_25 篩。N ≤ 10^12）\n";

    api("【1】 prime_count(N) / prime_sum(N)", "long long / Int128", "N 以下の素数の個数 / 和（篩では足りない大きな N で。Lucy の方法）", "O(N^{3/4} / log N)");
    SHOW(prime_count(100), "2, 3, 5, …, 97");
    SHOW(prime_count(100000000000LL), "π(10^11)");
    SHOW(prime_sum(100), "100 以下の素数の和");
    SHOW(prime_sum(10000000000LL), "10^10 以下の素数の和（Int128）");

    api("【2】 PrimeCountTable P(N); P.count(v)", "long long", "⌊N/i⌋ の形のすべての v での π(v)（√N 以下の v も）。乗法的関数の和を求める DP の部品", "前処理 O(N^{3/4} / log N)・O(1)");
    RUN(PrimeCountTable P(100));
    SHOW(P.count(100 / 3), "π(33)");
    SHOW(P(7), "π(7)");
    note("例：N 以下で「素数 p < q の積 p·q」になる数 = Σ_{p ≤ √N} (π(N/p) − π(p))（N/p は ⌊N/i⌋ の形なので表にある）");
    {
        auto semiprimes = [](long long N) {
            const PrimeCountTable T(N);
            long long c = 0;
            for (long long p = 2; p * p <= N; ++p)
                if (T(p) != T(p - 1)) c += T(N / p) - T(p);
            return c;
        };
        std::cout << "    N = 100 → " << semiprimes(100) << "（6, 10, 14, 15, …, 95）、N = 10^12 → " << semiprimes(1000000000000LL) << "\n";
    }

    api("【3】 count_smooth(N, P)", "long long", "N 以下で素因数がすべて P 以下の数の個数（1 も数える。素数を 2 組に分けて数え合わせる）", "O(組の大きさ log)");
    SHOW(count_smooth(36, 3), "ABC300 G 例 1：2^a 3^b ≤ 36 は 14 個");
    SHOW(count_smooth(10000000000000000LL, 97), "ABC300 G 例 2");

    api("【4】 multiplicative_sum<T>(N, coef, fpe)", "T（mint・long long・Int128）",
        "Σ_{i=1}^{N} f(i)（f は乗法的）。素数で f(p) = coef[0] + coef[1] p + …（3 次まで）、素数の累乗で f(p^e) = fpe(p, e)。mod なら T = mint",
        "O(N^{3/4} / log N) 程度（10^11 で約 1 秒）");
    SHOW(multiplicative_sum<long long>(10000000000LL, {2}, [](long long, int e) { return e + 1LL; }), "Σ d(i)（約数の個数）、N = 10^10");
    SHOW(multiplicative_sum<long long>(10000000000LL, {-1}, [](long long, int e) { return e == 1 ? -1LL : 0LL; }), "Σ μ(i)（メルテンス関数）");
    SHOW(multiplicative_sum<long long>(10000000000LL, {1}, [](long long, int e) { return e == 1 ? 1LL : 0LL; }), "平方因子を持たない数の個数（μ(i)^2 の和）");
    SHOW(multiplicative_sum<Int128>(10000000000LL, {-1, 1}, [](long long p, int e) { return (Int128)pow_int(p, e) - pow_int(p, e - 1); }), "Σ φ(i)：f(p) = p − 1 → coef {−1, 1}");
    SHOW(multiplicative_sum<Int128>(1000000000LL, {1, 1}, [](long long p, int e) { return (Int128)(pow_int(p, e + 1) - 1) / (p - 1); }), "Σ σ(i)（約数の和）、N = 10^9");
    note("fpe(p, 1) は f(p)（coef の多項式）と同じ値にする。整数の T は溢れたら止まるので、mod の答えなら 00_modint の mint を T に");
}

}  // namespace pcnt_demo

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
    pcnt_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pcnt_test::benchmark();    // 速度計測
    pcnt_demo::demo();         // 使い方デモ
    return 0;
}
