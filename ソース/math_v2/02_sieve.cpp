// =============================================================================
//  篩（Sieve：N 以下の数の素因数分解・素数判定を何度も、φ・μ・約数の個数などの表、素数の一覧、区間篩）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 素数・約数
//  依存       : 00_int_common, 00_math_base, 02_factorize（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  Sieve S(N) は N を超える数も受け付ける（02_factorize の is_prime / factorize に任せる。このファイルに入っている）
//
//  機能一覧                                                                    計算量
//    Sieve S(N)                       線形篩（最小の素因数の表）。N ≤ 5×10^7 程度（メモリ 4N バイト。上限 10^8）   O(N)
//      S.is_prime(x) / S.spf(x)       素数か / 最小の素因数（x ≤ N なら O(1)）
//      S.primes()                     N 以下の素数（昇順の vector<long long>。p * p も溢れない）
//      S.factorize(x)                 素因数分解 vector<PrimePower>（x ≤ N なら O(log x)）
//      S.prime_factors(x) / S.divisors(x)          異なる素因数 / 約数（昇順）
//      S.count_divisors(x) / S.sum_divisors(x)     約数の個数 / 約数の和（Int128）
//      S.phi(x) / S.mu(x) / S.omega(x) / S.big_omega(x)   φ / μ / 異なる素因数の数 / 素因数の数（重複も）
//      S.largest_prime_factor(x) / S.square_free_decompose(x) / S.is_square_free(x)
//      S.prime_count(x)               x 以下の素数の数（x ≤ N）                                     O(log N)
//      S.count_primes(l, r) / S.count_primes_closed(L, R)   [l, r) / [L, R] の素数の数（r ≤ N + 1）
//      S.kth_prime(k)                 小さい方から k 番目の素数（0 番目が 2）
//      S.lcm_mod(vs, m)               多くの数の lcm を mod m で（vector<int> も可）
//    （以下の表は 0 ≤ N ≤ 10^8。添字がそのまま x）
//    phi_table(N) / mu_table(N)       0〜N の φ / μ の表                                            O(N)
//    divisor_count_table(N) / divisor_sum_table(N)    0〜N の約数の個数 / 約数の和の表（線形篩）          O(N)
//    omega_table(N)                   0〜N の異なる素因数の数の表                                    O(N)
//    （表はどれも 10^7 で 0.1〜0.2 秒、10^8 で 1〜2 秒。10^8 のときメモリは phi_table・divisor_sum_table が約 900MB、ほかは約 500MB）
//    primes_up_to(N)                  N 以下の素数 vector<long long>（N < 2^31。区間ごとのエラトステネス。10^8 で約 0.3 秒）
//    primes_up_to<int>(N)             同じものを vector<int> で（メモリ半分。N = 10^9 なら 5×10^7 個で 200MB）
//    segment_primes(l, r) / segment_primes_closed(L, R)     [l, r) / [L, R] の素数（r ≤ 10^14、幅 ≤ 10^7 程度）
//    segment_factorize(l, r, f) / segment_factorize_closed(L, R, f)
//                                     [l, r) の各数 x を小さい順に f(x, 素因数分解) へ渡す（幅 ≤ 10^7 程度）
//    （1 つの大きな数の素数判定・素因数分解は 02_factorize：is_prime / factorize / divisors_of など。このファイルに入っている）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    sieve_test::random_test()  … 愚直な篩・試し割りとの突き合わせ
//    sieve_test::benchmark()    … 速度計測
//    sieve_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
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
#ifndef GG_FACTORIZE
#define GG_FACTORIZE
// ---- 素因数分解（数学第二世代）  1 つの大きな数（n < 2^63）の素数判定・素因数分解・約数・φ・μ など ----
struct PrimePower {  // 素因数分解の 1 項 p^e（素因数分解は p の昇順の vector<PrimePower>）
    long long p;
    int e;
    friend bool operator==(const PrimePower& a, const PrimePower& b) { return a.p == b.p && a.e == b.e; }
    friend bool operator!=(const PrimePower& a, const PrimePower& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const PrimePower& f) { return f.e == 1 ? os << f.p : os << f.p << "^" << f.e; }
};
struct SquareFreeParts {  // n = a^2 · b（b は平方因子を持たない）
    long long a, b;
};

namespace math_internal {
struct Mont64 {  // 64 ビットのモンゴメリ乗算（n は奇数。値は [0, n) で持つ）
    unsigned long long n, ninv, r2;
    explicit Mont64(unsigned long long n_) : n(n_), ninv(n_), r2(0) {
        for (int i = 0; i < 6; ++i) ninv *= 2 - n * ninv;  // n^{-1} mod 2^64
        r2 = (unsigned long long)(((UInt128)0 - (UInt128)n) % n);  // 2^128 mod n
    }
    unsigned long long reduce(UInt128 t) const {  // t < n·2^64 → t·2^{-64} mod n
        const unsigned long long q = (unsigned long long)t * ninv;
        const unsigned long long hi = (unsigned long long)(t >> 64), m = (unsigned long long)(((UInt128)q * n) >> 64);
        return hi >= m ? hi - m : hi - m + n;
    }
    unsigned long long mul(unsigned long long a, unsigned long long b) const { return reduce((UInt128)a * b); }
    unsigned long long to(unsigned long long a) const { return mul(a % n, r2); }
    unsigned long long from(unsigned long long a) const { return reduce(a); }
    unsigned long long pow(unsigned long long a, unsigned long long e) const {  // a はモンゴメリ表現
        unsigned long long r = to(1);
        for (; e; e >>= 1, a = mul(a, a))
            if (e & 1) r = mul(r, a);
        return r;
    }
};
inline bool miller_rabin(unsigned long long n) {  // n は 5 以上の奇数
    unsigned long long d = n - 1;
    int s = 0;
    while (!(d & 1)) d >>= 1, ++s;
    if (n < (1ULL << 32)) {
        for (unsigned long long a : {2ULL, 7ULL, 61ULL}) {
            if (a % n == 0) continue;
            unsigned long long x = 1, b = a % n;
            for (unsigned long long e = d; e; e >>= 1, b = b * b % n)
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
    const Mont64 mt(n);
    const unsigned long long one = mt.to(1), mone = mt.to(n - 1);
    for (unsigned long long a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (a % n == 0) continue;
        unsigned long long x = mt.pow(mt.to(a), d);
        if (x == one || x == mone) continue;
        bool comp = true;
        for (int i = 1; i < s && comp; ++i) {
            x = mt.mul(x, x);
            if (x == mone) comp = false;
        }
        if (comp) return false;
    }
    return true;
}
inline unsigned long long find_factor(unsigned long long n) {  // 合成数 n の 1 でも n でもない約数（ポラード・ロー、Brent 版）
    if (n % 2 == 0) return 2;
    const Mont64 mt(n);
    for (unsigned long long c0 = 1;; ++c0) {
        const unsigned long long c = mt.to(c0);
        auto f = [&](unsigned long long x) {
            const unsigned long long y = mt.mul(x, x) + c;
            return y >= n ? y - n : y;
        };
        unsigned long long x = 0, y = mt.to(2), ys = 0, q = mt.to(1), g = 1;
        const unsigned long long M = 128;
        for (unsigned long long r = 1; g == 1; r <<= 1) {
            x = y;
            for (unsigned long long i = 0; i < r; ++i) y = f(y);
            for (unsigned long long k = 0; k < r && g == 1; k += M) {
                ys = y;
                for (unsigned long long i = 0; i < M && i < r - k; ++i) y = f(y), q = mt.mul(q, x > y ? x - y : y - x);
                g = ugcd(q, n);
            }
        }
        if (g == n) {
            do {
                ys = f(ys);
                g = ugcd(x > ys ? x - ys : ys - x, n);
            } while (g == 1);
        }
        if (g != n) return g;
    }
}
inline std::vector<PrimePower> group_primes(std::vector<long long>& ps) {
    std::sort(ps.begin(), ps.end());
    std::vector<PrimePower> fs;
    for (long long p : ps) {
        if (!fs.empty() && fs.back().p == p) ++fs.back().e;
        else fs.push_back({p, 1});
    }
    return fs;
}
inline void check_positive(long long n, const char* func) {
    if (n < 1) fail(func, "n ≥ 1 で使う");
}
}  // namespace math_internal

inline bool is_prime(long long n) {  // 素数か（n < 2^63。決定的ミラー・ラビン）
    if (n < 2) return false;
    for (long long p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        if (n % p == 0) return n == p;
    }
    if (n < 37 * 37) return true;
    return math_internal::miller_rabin((unsigned long long)n);
}
inline std::vector<PrimePower> factorize(long long n) {  // 素因数分解（n ≥ 1。p の昇順。n = 1 なら空。ポラード・ロー）
    math_internal::check_positive(n, "factorize");
    std::vector<long long> ps;
    for (long long p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47})
        while (n % p == 0) ps.push_back(p), n /= p;
    std::vector<unsigned long long> st;
    if (n > 1) st.push_back((unsigned long long)n);
    while (!st.empty()) {
        const unsigned long long m = st.back();
        st.pop_back();
        if (m < 53 * 53 || is_prime((long long)m)) {
            ps.push_back((long long)m);
            continue;
        }
        const unsigned long long d = math_internal::find_factor(m);
        st.push_back(d), st.push_back(m / d);
    }
    return math_internal::group_primes(ps);
}
inline std::vector<PrimePower> factorize_trial(long long n) {  // 試し割りの素因数分解（n ≥ 1。O(√n)。確かめ用）
    math_internal::check_positive(n, "factorize_trial");
    std::vector<PrimePower> fs;
    for (long long p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
        if (n % p) continue;
        int e = 0;
        while (n % p == 0) n /= p, ++e;
        fs.push_back({p, e});
    }
    if (n > 1) fs.push_back({n, 1});
    return fs;
}
inline std::vector<long long> prime_factors(long long n) {  // 異なる素因数（昇順）
    std::vector<long long> r;
    for (const PrimePower& f : factorize(n)) r.push_back(f.p);
    return r;
}
template <class F>
void for_each_divisor(const std::vector<PrimePower>& fs, F f) {  // 素因数分解から全部の約数を f(d) に渡す（順番は決まっていない）
    std::vector<long long> ds = {1};
    for (const PrimePower& pe : fs) {
        const size_t sz = ds.size();
        long long pk = 1;
        for (int k = 1; k <= pe.e; ++k) {
            pk *= pe.p;
            for (size_t i = 0; i < sz; ++i) ds.push_back(ds[i] * pk);
        }
    }
    for (long long d : ds) f(d);
}
inline std::vector<long long> divisors_from_factors(const std::vector<PrimePower>& fs) {  // 素因数分解から約数（昇順）
    std::vector<long long> ds;
    for_each_divisor(fs, [&](long long d) { ds.push_back(d); });
    std::sort(ds.begin(), ds.end());
    return ds;
}
inline std::vector<long long> divisors_of(long long n) { return divisors_from_factors(factorize(n)); }  // 約数（昇順。n ≥ 1）
inline long long count_divisors(long long n) {  // 約数の個数
    long long c = 1;
    for (const PrimePower& f : factorize(n)) c *= f.e + 1;
    return c;
}
inline Int128 sum_divisors(long long n) {  // 約数の和（Int128）
    Int128 s = 1;
    for (const PrimePower& f : factorize(n)) {
        Int128 t = 1, pk = 1;
        for (int k = 0; k < f.e; ++k) pk *= f.p, t += pk;
        s *= t;
    }
    return s;
}
inline long long euler_phi(long long n) {  // オイラーの φ(n)（1 以上 n 以下で n と互いに素な数の個数）
    long long r = n;
    for (const PrimePower& f : factorize(n)) r = r / f.p * (f.p - 1);
    return r;
}
inline int mobius(long long n) {  // メビウスの μ(n)（平方因子があれば 0、なければ (-1)^(素因数の数)）
    const std::vector<PrimePower> fs = factorize(n);
    for (const PrimePower& f : fs)
        if (f.e >= 2) return 0;
    return fs.size() % 2 ? -1 : 1;
}
inline int omega(long long n) { return (int)factorize(n).size(); }  // 異なる素因数の個数 ω(n)（12 = 2^2·3 → 2）
inline int big_omega(long long n) {  // 重複を含めた素因数の個数 Ω(n)（12 = 2^2·3 → 3）
    int c = 0;
    for (const PrimePower& f : factorize(n)) c += f.e;
    return c;
}
inline long long largest_prime_factor(long long n) {  // 最大の素因数（n = 1 なら 1）
    const std::vector<PrimePower> fs = factorize(n);
    return fs.empty() ? 1 : fs.back().p;
}
inline SquareFreeParts square_free_decompose(long long n) {  // n = a^2 · b（b は平方因子を持たない）
    long long a = 1, b = 1;
    for (const PrimePower& f : factorize(n)) {
        for (int k = 0; k < f.e / 2; ++k) a *= f.p;
        if (f.e % 2) b *= f.p;
    }
    return {a, b};
}
inline bool is_square_free(long long n) { return square_free_decompose(n).a == 1; }  // 平方因子を持たないか
inline long long radical(long long n) {  // 異なる素因数の積
    long long r = 1;
    for (const PrimePower& f : factorize(n)) r *= f.p;
    return r;
}
inline bool is_prime_power(long long n) {  // n = p^k（p は素数、k ≥ 1）か
    return n >= 2 && factorize(n).size() == 1;
}
inline long long factorial_exponent(long long n, long long p) {  // n! が素数 p で何回割れるか（ルジャンドル。n ≥ 0）
    if (p < 2) math_internal::fail("factorial_exponent", "p ≥ 2 で使う");
    long long c = 0;
    for (long long m = n / p; m > 0; m /= p) c += m;
    return c;
}
inline long long factorial_max_power(long long n, long long k) {  // k^e が n! を割り切る最大の e（k ≥ 2 は合成数でもよい）
    if (k < 2) math_internal::fail("factorial_max_power", "k ≥ 2 で使う");
    long long best = std::numeric_limits<long long>::max();
    for (const PrimePower& f : factorize(k)) best = std::min(best, factorial_exponent(n, f.p) / f.e);
    return best;
}
inline long long min_factorial_multiple(long long k) {  // k が n! を割り切る最小の n ≥ 1（k ≥ 1）
    long long ans = 1;
    for (const PrimePower& f : factorize(k)) {
        long long n = 0, cnt = 0;
        while (cnt < f.e) {
            n += f.p;
            for (long long m = n; m % f.p == 0; m /= f.p) ++cnt;
        }
        ans = std::max(ans, n);
    }
    return ans;
}
inline bool is_perfect_power(long long n) {  // n = a^b（a ≥ 1, b ≥ 2）と書けるか（1 は true。0 以下は false）
    if (n < 1) return false;
    if (n == 1) return true;
    for (int b = 2; b < 63; ++b) {
        const long long a = kth_root_floor(n, b);
        if (a < 2) break;
        if (is_kth_power(n, b)) return true;
    }
    return false;
}
inline long long count_perfect_powers(long long N) {  // 1 以上 N 以下の、a^b（a ≥ 1, b ≥ 2）と書ける数の個数（1 も数える）
    if (N < 1) return 0;
    std::vector<long long> v;  // 平方数でない 3 乗以上の累乗
    for (long long a = 2; a <= N / a / a; ++a) {
        if (is_square(a)) continue;  // a が平方数なら a^b は (√a)^(2b) で、下の数え方に含まれる
        Int128 x = (Int128)a * a * a;
        for (int b = 3; x <= N; ++b, x *= a)
            if (b % 2 == 1) v.push_back((long long)x);  // b が偶数なら平方数
    }
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return isqrt(N) + (long long)v.size();
}
inline long long count_coprime_up_to(long long N, long long M) {  // 1 以上 N 以下で M（≥ 1）と互いに素な数の個数（包除原理）
    if (N < 1) return 0;
    const std::vector<long long> ps = prime_factors(M);
    const int k = (int)ps.size();
    long long r = 0;
    for (int s = 0; s < (1 << k); ++s) {
        Int128 prod = 1;
        int bits = 0;
        for (int i = 0; i < k; ++i)
            if (s >> i & 1) prod *= ps[i], ++bits;
        if (prod > N) continue;
        r += (bits % 2 ? -1 : 1) * (N / (long long)prod);
    }
    return r;
}
template <class V>
long long lcm_mod_of(const V& vs, long long m) {  // 全部の数（≥ 1）の lcm を mod m で（m は何でも。0 があれば 0。vector<int> なども可）
    if (m < 1) math_internal::fail("lcm_mod_of", "m ≥ 1 で使う");
    std::vector<long long> xs;
    for (const auto& xv : vs) {  // 先に全部を確かめる
        const long long x = math_internal::to_ll_checked(xv, "lcm_mod_of");
        if (x < 0) math_internal::fail("lcm_mod_of", "数は 0 以上で使う");
        xs.push_back(x);
    }
    std::vector<PrimePower> best;
    for (const long long x : xs) {
        if (x == 0) return 0;
        for (const PrimePower& f : factorize(x)) best.push_back(f);
    }
    std::sort(best.begin(), best.end(), [](const PrimePower& a, const PrimePower& b) { return a.p != b.p ? a.p < b.p : a.e > b.e; });
    long long r = 1 % m;
    for (size_t i = 0; i < best.size(); ++i)
        if (i == 0 || best[i].p != best[i - 1].p) r = mod_mul(r, mod_pow(best[i].p, best[i].e, m), m);
    return r;
}
inline long long lcm_mod_of(std::initializer_list<long long> vs, long long m) { return lcm_mod_of<std::initializer_list<long long>>(vs, m); }
#endif  // GG_FACTORIZE
#ifndef GG_SIEVE
#define GG_SIEVE
// ---- 篩（数学第二世代）  Sieve S(N)：N 以下の数を何度も素因数分解・素数判定、φ・μ などの表、区間篩 ----
// Sieve は N を超える数も受け付ける（02_factorize の is_prime / factorize に任せる）
namespace math_internal {
inline int checked_table_size(long long N, const char* func) {  // 0〜N の表を作る前に大きさを確かめる
    if (N < 0) fail(func, "N ≥ 0 で使う");
    if (N > 100000000) fail(func, "N が大きすぎる（10^8 まで）");
    return (int)N;
}
}  // namespace math_internal
struct Sieve {
  private:
    int n_;
    std::vector<int> spf_;
    std::vector<long long> primes_;
    static int checked_size(long long N) {  // 表を作る前に大きさを確かめる（N = 0 は 1 として作る）
        if (N < 0) math_internal::fail("Sieve", "N ≥ 0 で使う");
        if (N > 100000000) math_internal::fail("Sieve", "N が大きすぎる（10^8 まで。素数の一覧だけなら primes_up_to）");
        return (int)std::max(N, 1LL);
    }
    std::vector<PrimePower> small_factorize(long long x) const {
        std::vector<PrimePower> fs;
        while (x > 1) {
            const int p = spf_[x];
            int e = 0;
            while (x % p == 0) x /= p, ++e;
            fs.push_back({p, e});
        }
        return fs;
    }

  public:
    explicit Sieve(long long N) : n_(checked_size(N)), spf_((size_t)n_ + 1, 0) {  // 線形篩 O(N)。N ≤ 5×10^7 程度（メモリ 4N バイト）
        spf_[1] = 1;
        for (int i = 2; i <= n_; ++i) {
            if (spf_[i] == 0) spf_[i] = i, primes_.push_back(i);
            for (const long long p : primes_) {
                if (p > spf_[i] || p * i > n_) break;
                spf_[p * i] = (int)p;
            }
        }
    }
    int size() const { return n_; }
    const std::vector<long long>& primes() const { return primes_; }  // N 以下の素数（昇順。long long なので p * p も溢れない）
    int spf(long long x) const {                                        // 最小の素因数（2 ≤ x ≤ N。spf(1) = 1）
        if (x < 1 || x > n_) math_internal::fail("Sieve::spf", "1 ≤ x ≤ N で使う");
        return spf_[x];
    }
    bool is_prime(long long x) const { return x <= n_ ? (x >= 2 && spf_[x] == x) : ::is_prime(x); }
    std::vector<PrimePower> factorize(long long x) const {  // x ≥ 1（N を超えたら 02_factorize の factorize）
        math_internal::check_positive(x, "Sieve::factorize");
        return x <= n_ ? small_factorize(x) : ::factorize(x);
    }
    std::vector<long long> prime_factors(long long x) const {
        std::vector<long long> r;
        for (const PrimePower& f : factorize(x)) r.push_back(f.p);
        return r;
    }
    std::vector<long long> divisors(long long x) const { return divisors_from_factors(factorize(x)); }  // 昇順
    long long count_divisors(long long x) const {
        long long c = 1;
        for (const PrimePower& f : factorize(x)) c *= f.e + 1;
        return c;
    }
    Int128 sum_divisors(long long x) const {
        Int128 s = 1;
        for (const PrimePower& f : factorize(x)) {
            Int128 t = 1, pk = 1;
            for (int k = 0; k < f.e; ++k) pk *= f.p, t += pk;
            s *= t;
        }
        return s;
    }
    long long phi(long long x) const {
        long long r = x;
        for (const PrimePower& f : factorize(x)) r = r / f.p * (f.p - 1);
        return r;
    }
    int mu(long long x) const {
        const std::vector<PrimePower> fs = factorize(x);
        for (const PrimePower& f : fs)
            if (f.e >= 2) return 0;
        return fs.size() % 2 ? -1 : 1;
    }
    int omega(long long x) const { return (int)factorize(x).size(); }  // 異なる素因数の数
    int big_omega(long long x) const {                                  // 素因数の数（重複も数える）
        int c = 0;
        for (const PrimePower& f : factorize(x)) c += f.e;
        return c;
    }
    long long largest_prime_factor(long long x) const {  // x = 1 なら 1
        const std::vector<PrimePower> fs = factorize(x);
        return fs.empty() ? 1 : fs.back().p;
    }
    SquareFreeParts square_free_decompose(long long x) const {  // x = a^2 · b
        long long a = 1, b = 1;
        for (const PrimePower& f : factorize(x)) {
            for (int k = 0; k < f.e / 2; ++k) a *= f.p;
            if (f.e % 2) b *= f.p;
        }
        return {a, b};
    }
    bool is_square_free(long long x) const { return square_free_decompose(x).a == 1; }
    long long prime_count(long long x) const {  // x 以下の素数の数（x ≤ N）
        if (x > n_) math_internal::fail("Sieve::prime_count", "x ≤ N で使う（大きい x は 02_prime_count）");
        int lo = 0, hi = (int)primes_.size();  // primes_[lo..hi) で x を超える最初の位置を二分探索
        while (lo < hi) {
            const int mid = (lo + hi) / 2;
            if (primes_[mid] <= x) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
    long long count_primes(long long l, long long r) const {  // [l, r) の素数の数（r ≤ N + 1）
        if (l >= r) return 0;
        return prime_count(r - 1) - (l >= 1 ? prime_count(l - 1) : 0);
    }
    long long count_primes_closed(long long L, long long R) const {  // [L, R] の素数の数（R ≤ N）
        if (L > R) return 0;
        if (R > n_) math_internal::fail("Sieve::count_primes_closed", "R ≤ N で使う（大きい R は 02_prime_count）");
        return count_primes(L, R + 1);
    }
    long long kth_prime(long long k) const {  // 小さい方から k 番目の素数（0 番目が 2。N 以下の素数の中で）
        if (k < 0 || k >= (long long)primes_.size()) math_internal::fail("Sieve::kth_prime", "k が範囲の外（N を大きくする）");
        return primes_[k];
    }
    template <class V>
    long long lcm_mod(const V& vs, long long m) const {  // 全部の数（≥ 1）の lcm を mod m で（vector<int> なども可）
        if (m < 1) math_internal::fail("Sieve::lcm_mod", "m ≥ 1 で使う");
        std::vector<long long> xs;
        for (const auto& xv : vs) {  // 先に全部を確かめる
            const long long x = math_internal::to_ll_checked(xv, "Sieve::lcm_mod");
            if (x < 0) math_internal::fail("Sieve::lcm_mod", "数は 0 以上で使う");
            xs.push_back(x);
        }
        std::vector<PrimePower> all;
        for (const long long x : xs) {
            if (x == 0) return 0;
            for (const PrimePower& f : factorize(x)) all.push_back(f);
        }
        std::sort(all.begin(), all.end(), [](const PrimePower& a, const PrimePower& b) { return a.p != b.p ? a.p < b.p : a.e > b.e; });
        long long r = 1 % m;
        for (size_t i = 0; i < all.size(); ++i)
            if (i == 0 || all[i].p != all[i - 1].p) r = mod_mul(r, mod_pow(all[i].p, all[i].e, m), m);
        return r;
    }
    long long lcm_mod(std::initializer_list<long long> vs, long long m) const { return lcm_mod<std::initializer_list<long long>>(vs, m); }
};

// ---- 0〜N の表（O(N)。0 ≤ N ≤ 10^8）----
inline std::vector<long long> phi_table(long long N_) {  // phi[x] = φ(x)（phi[0] = 0）
    const int N = math_internal::checked_table_size(N_, "phi_table");
    std::vector<long long> phi(N + 1);
    std::vector<int> ps;
    for (int i = 0; i <= N; ++i) phi[i] = i;
    std::vector<char> comp(N + 1, 0);
    for (int i = 2; i <= N; ++i) {
        if (!comp[i]) ps.push_back(i), phi[i] = i - 1;
        for (int p : ps) {
            if ((long long)p * i > N) break;
            comp[p * i] = 1;
            if (i % p == 0) {
                phi[p * i] = phi[i] * p;
                break;
            }
            phi[p * i] = phi[i] * (p - 1);
        }
    }
    return phi;
}
inline std::vector<int> mu_table(long long N_) {  // mu[x] = μ(x)（mu[0] = 0）
    const int N = math_internal::checked_table_size(N_, "mu_table");
    std::vector<int> mu(N + 1, 0), ps;
    std::vector<char> comp(N + 1, 0);
    if (N >= 1) mu[1] = 1;
    for (int i = 2; i <= N; ++i) {
        if (!comp[i]) ps.push_back(i), mu[i] = -1;
        for (int p : ps) {
            if ((long long)p * i > N) break;
            comp[p * i] = 1;
            if (i % p == 0) {
                mu[p * i] = 0;
                break;
            }
            mu[p * i] = -mu[i];
        }
    }
    return mu;
}
// 以下も線形篩（x = p·i（p は x の最小の素因数）を i から決める。e[x] は最小の素因数の指数）
inline std::vector<int> divisor_count_table(long long N_) {  // d[x] = x の約数の個数（d[0] = 0）  O(N)
    const int N = math_internal::checked_table_size(N_, "divisor_count_table");
    std::vector<int> d(N + 1, 0), ps;
    std::vector<unsigned char> e(N + 1, 0);
    if (N >= 1) d[1] = 1;
    for (int i = 2; i <= N; ++i) {
        if (d[i] == 0) ps.push_back(i), d[i] = 2, e[i] = 1;  // まだ決まっていない = 素数
        for (int p : ps) {
            if ((long long)p * i > N) break;
            const int j = p * i;
            if (i % p == 0) {
                e[j] = (unsigned char)(e[i] + 1), d[j] = d[i] / (e[i] + 1) * (e[i] + 2);
                break;
            }
            e[j] = 1, d[j] = d[i] * 2;
        }
    }
    return d;
}
inline std::vector<long long> divisor_sum_table(long long N_) {  // s[x] = x の約数の和（s[0] = 0）  O(N)
    const int N = math_internal::checked_table_size(N_, "divisor_sum_table");
    std::vector<long long> s(N + 1, 0);
    std::vector<unsigned char> e(N + 1, 0);
    std::vector<int> ps;
    if (N >= 1) s[1] = 1;
    for (int i = 2; i <= N; ++i) {
        if (s[i] == 0) ps.push_back(i), s[i] = i + 1, e[i] = 1;
        for (int p : ps) {
            if ((long long)p * i > N) break;
            const int j = p * i;
            if (i % p == 0) {  // i = p^e · rest → j = p^{e+1} · rest
                long long pk = 1, term = 1, pw = 1;
                for (int t = 0; t < e[i]; ++t) pk *= p;
                for (int t = 0; t <= e[i]; ++t) pw *= p, term += pw;
                e[j] = (unsigned char)(e[i] + 1), s[j] = s[i / pk] * term;
                break;
            }
            e[j] = 1, s[j] = s[i] * (p + 1);
        }
    }
    return s;
}
inline std::vector<int> omega_table(long long N_) {  // w[x] = x の異なる素因数の数（w[0] = w[1] = 0）  O(N)
    const int N = math_internal::checked_table_size(N_, "omega_table");
    std::vector<int> w(N + 1, 0), ps;
    std::vector<char> comp(N + 1, 0);
    for (int i = 2; i <= N; ++i) {
        if (!comp[i]) ps.push_back(i), w[i] = 1;
        for (int p : ps) {
            if ((long long)p * i > N) break;
            const int j = p * i;
            comp[j] = 1;
            if (i % p == 0) {
                w[j] = w[i];
                break;
            }
            w[j] = w[i] + 1;
        }
    }
    return w;
}
// N 以下の素数（エラトステネス。奇数だけ・区間ごと。N < 2^31。10^9 なら 5×10^7 個・約 3 秒）
// 返り値は vector<long long>（p * p も溢れない）。メモリを半分にするなら primes_up_to<int>(N)（N = 10^9 で 400MB → 200MB）
template <class T = long long>
std::vector<T> primes_up_to(long long N) {
    if (N > 2147483647LL) math_internal::fail("primes_up_to", "N < 2^31 で使う");
    std::vector<T> ps;
    if (N < 2) return ps;
    ps.push_back(2);
    const int n = (int)N, sq = (int)isqrt(N);
    std::vector<char> small(sq + 1, 1);
    std::vector<int> base;  // √N 以下の奇素数
    for (int i = 3; i <= sq; i += 2)
        if (small[i]) {
            base.push_back(i);
            for (long long j = (long long)i * i; j <= sq; j += 2 * i) small[j] = 0;
        }
    const int B = 1 << 16;  // 奇数 B 個ずつ
    std::vector<char> blk(B);
    std::vector<long long> nxt(base.size());  // 次に消す奇数の倍数
    for (size_t i = 0; i < base.size(); ++i) nxt[i] = (long long)base[i] * base[i];
    for (long long lo = 3; lo <= n; lo += 2LL * B) {  // [lo, lo + 2B) の奇数
        const long long hi = std::min<long long>(lo + 2LL * B, (long long)n + 1);
        std::fill(blk.begin(), blk.end(), 1);
        for (size_t i = 0; i < base.size(); ++i) {
            long long j = nxt[i];
            const long long step = 2LL * base[i];
            for (; j < hi; j += step) blk[(j - lo) >> 1] = 0;
            nxt[i] = j;
        }
        for (long long x = lo; x < hi; x += 2)
            if (blk[(x - lo) >> 1]) ps.push_back((T)x);
    }
    return ps;
}

// ---- 区間篩（[l, r) の数。r ≤ 10^14 程度、r − l ≤ 10^7 程度。√r 以下の素数で割る）----
inline std::vector<long long> segment_primes(long long l, long long r) {  // [l, r) の素数（昇順）
    std::vector<long long> res;
    l = std::max(l, 2LL);
    if (l >= r) return res;
    const long long sq = isqrt(r - 1);
    std::vector<char> small(sq + 1, 1), seg((size_t)(r - l), 1);
    for (long long p = 2; p <= sq; ++p) {
        if (!small[p]) continue;
        for (long long j = p * p; j <= sq; j += p) small[j] = 0;
        for (long long j = std::max(p * p, ceil_multiple(l, p)); j < r; j += p) seg[(size_t)(j - l)] = 0;
    }
    for (long long x = l; x < r; ++x)
        if (seg[(size_t)(x - l)]) res.push_back(x);
    return res;
}
inline std::vector<long long> segment_primes_closed(long long L, long long R) { return segment_primes(L, R + 1); }
// [l, r) の各数 x（≥ 1）の素因数分解を、x の小さい順に f(x, fs) へ渡す（fs は p の昇順の vector<PrimePower>）
template <class F>
void segment_factorize(long long l, long long r, F f) {
    l = std::max(l, 1LL);
    if (l >= r) return;
    const long long sq = isqrt(r - 1);
    std::vector<char> small(sq + 1, 1);
    std::vector<long long> base;
    for (long long p = 2; p <= sq; ++p) {
        if (!small[p]) continue;
        base.push_back(p);
        for (long long j = p * p; j <= sq; j += p) small[j] = 0;
    }
    const long long B = 1 << 18;
    std::vector<long long> rem, nxt(base.size());
    for (size_t k = 0; k < base.size(); ++k) nxt[k] = ceil_multiple(l, base[k]);  // 次に割る倍数
    std::vector<std::vector<PrimePower>> fs((size_t)std::min(B, r - l));
    for (long long lo = l; lo < r; lo += B) {
        const long long hi = std::min(lo + B, r);
        rem.resize((size_t)(hi - lo));
        for (long long x = lo; x < hi; ++x) rem[(size_t)(x - lo)] = x, fs[(size_t)(x - lo)].clear();
        for (size_t k = 0; k < base.size(); ++k) {
            const long long p = base[k];
            long long j = nxt[k];
            for (; j < hi; j += p) {
                long long& v = rem[(size_t)(j - lo)];
                int e = 0;
                while (v % p == 0) v /= p, ++e;
                fs[(size_t)(j - lo)].push_back({p, e});
            }
            nxt[k] = j;
        }
        for (long long x = lo; x < hi; ++x) {
            const size_t i = (size_t)(x - lo);
            if (rem[i] > 1) fs[i].push_back({rem[i], 1});
            f(x, (const std::vector<PrimePower>&)fs[i]);
        }
    }
}
template <class F>
void segment_factorize_closed(long long L, long long R, F f) { segment_factorize(L, R + 1, f); }
#endif  // GG_SIEVE
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace sieve_test {

#define SIEVE_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[sieve] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static std::vector<char> brute_sieve(int N) {
    std::vector<char> ip(N + 1, 1);
    ip[0] = 0;
    if (N >= 1) ip[1] = 0;
    for (int i = 2; (long long)i * i <= N; ++i)
        if (ip[i])
            for (int j = i * i; j <= N; j += i) ip[j] = 0;
    return ip;
}

static void check_members(const Sieve& S, long long x) {
    const std::vector<PrimePower> fs = factorize_trial(x);
    SIEVE_CHECK(S.factorize(x) == fs);
    long long cd = 1, phi = x, a = 1, b = 1, lp = 1;
    Int128 sd = 1;
    int mu = 1, big = 0;
    for (const PrimePower& f : fs) {
        cd *= f.e + 1, phi = phi / f.p * (f.p - 1), big += f.e, lp = f.p;
        Int128 t = 1, pk = 1;
        for (int k = 0; k < f.e; ++k) pk *= f.p, t += pk;
        sd *= t;
        mu = f.e >= 2 ? 0 : -mu;
        for (int k = 0; k < f.e / 2; ++k) a *= f.p;
        if (f.e % 2) b *= f.p;
    }
    SIEVE_CHECK(S.count_divisors(x) == cd && S.sum_divisors(x) == sd && S.phi(x) == phi && S.mu(x) == mu);
    SIEVE_CHECK(S.omega(x) == (int)fs.size() && S.big_omega(x) == big && S.largest_prime_factor(x) == lp);
    const SquareFreeParts sp = S.square_free_decompose(x);
    SIEVE_CHECK(sp.a == a && sp.b == b && S.is_square_free(x) == (a == 1) && S.prime_factors(x).size() == fs.size());
    SIEVE_CHECK(S.is_prime(x) == (fs.size() == 1 && fs[0].e == 1));
}

static void test_sieve(std::mt19937_64& rng, int& cases) {
    for (int N : {1, 2, 3, 4, 10, 97, 100, 1000, 65537, 200000}) {
        const Sieve S(N);
        const std::vector<char> ip = brute_sieve(N);
        std::vector<long long> ps;
        for (int x = 0; x <= N; ++x)
            if (ip[x]) ps.push_back(x);
        SIEVE_CHECK(S.primes() == ps && S.size() == N);
        long long cnt = 0;
        for (int x = 0; x <= N; ++x, ++cases) {
            cnt += ip[x];
            SIEVE_CHECK(S.is_prime(x) == (bool)ip[x] && S.prime_count(x) == cnt);
            if (x >= 2) SIEVE_CHECK(S.spf(x) == factorize_trial(x)[0].p);
            if (x >= 1 && (N <= 1000 || x % 7 == 0)) check_members(S, x);
        }
        for (int k = 0; k < (int)ps.size(); ++k) SIEVE_CHECK(S.kth_prime(k) == ps[k]);
        for (int it = 0; it < 200; ++it, ++cases) {
            const long long l = rnd(rng, -3, N + 1), r = rnd(rng, -3, N + 1);
            long long c = 0;
            for (long long x = std::max(0LL, l); x < r; ++x) c += ip[x];
            SIEVE_CHECK(S.count_primes(l, r) == c && S.count_primes_closed(l, r - 1) == c);
        }
        // N を超える数は 02_factorize に任せる
        for (int it = 0; it < 300; ++it, ++cases) {
            const long long x = rnd(rng, N + 1, it % 2 ? 10000000000LL : 1000000000000000000LL);
            SIEVE_CHECK(S.factorize(x) == factorize(x) && S.is_prime(x) == is_prime(x));
            if (x <= 10000000000LL) check_members(S, x);
        }
    }
    // lcm_mod
    const Sieve S(100000);
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        std::vector<long long> vs((size_t)rnd(rng, 0, 6));
        for (auto& v : vs) v = rnd(rng, 1, it % 2 ? 100000 : 1000000000000LL);
        const long long m = rnd(rng, 1, 1000000007);
        SIEVE_CHECK(S.lcm_mod(vs, m) == lcm_mod_of(vs, m));
    }
    SIEVE_CHECK(S.lcm_mod(std::vector<int>{4, 6, 10}, 1000) == 60 && S.lcm_mod({8, 12}, 1000) == 24 && S.count_primes_closed(2, 2) == 1 && S.count_primes_closed(5, 4) == 0);
    {  // long long の N でも作れる・primes() は long long
        const long long NL = 1000;
        const Sieve S2(NL);
        long long sq = 0;
        for (auto p : S2.primes()) sq += p * p;
        SIEVE_CHECK(S2.size() == 1000 && sq > 0 && std::is_same<std::decay_t<decltype(S2.primes()[0])>, long long>::value);
    }
}

static void test_tables(std::mt19937_64& rng, int& cases) {
    for (int N : {0, 1, 2, 3, 10, 1000, 100000}) {
        const Sieve S(std::max(N, 1));
        const std::vector<long long> phi = phi_table(N), ds = divisor_sum_table(N);
        const std::vector<int> mu = mu_table(N), dc = divisor_count_table(N), om = omega_table(N);
        SIEVE_CHECK((int)phi.size() == N + 1 && (int)mu.size() == N + 1 && (int)ds.size() == N + 1 && (int)dc.size() == N + 1 && (int)om.size() == N + 1);
        SIEVE_CHECK(phi[0] == 0 && mu[0] == 0 && dc[0] == 0 && ds[0] == 0 && om[0] == 0);
        for (int x = 1; x <= N; ++x, ++cases)
            SIEVE_CHECK(phi[x] == S.phi(x) && mu[x] == S.mu(x) && dc[x] == S.count_divisors(x) && ds[x] == (long long)S.sum_divisors(x) && om[x] == S.omega(x));
    }
    // primes_up_to（区間の境目の前後も）
    const int M = 3000000;
    const std::vector<char> ip = brute_sieve(M);
    std::vector<int> all;
    for (int x = 0; x <= M; ++x)
        if (ip[x]) all.push_back(x);
    std::vector<long long> Ns = {0, 1, 2, 3, 4, 5, 9, 25, 131073, 131074, 131075, 131076, 262147, 262149, M};
    for (int it = 0; it < 20; ++it) Ns.push_back(rnd(rng, 0, M));
    for (long long N : Ns) {
        const std::vector<int> ps = primes_up_to<int>(N);
        const std::vector<long long> pl = primes_up_to(N);
        size_t k = 0;
        while (k < all.size() && all[k] <= N) ++k;
        SIEVE_CHECK(ps == std::vector<int>(all.begin(), all.begin() + (long)k) && pl == std::vector<long long>(all.begin(), all.begin() + (long)k));
        ++cases;
    }
    SIEVE_CHECK(primes_up_to(20000000).size() == 1270607 && primes_up_to<int>(20000000).size() == 1270607);
    SIEVE_CHECK(phi_table(3LL).size() == 4 && omega_table(0).size() == 1);
    {  // 大きめの N で、約数を全部たどる方法と比べる
        const int N = 2000000;
        std::vector<int> dc(N + 1, 0), om(N + 1, 0);
        std::vector<long long> ds(N + 1, 0);
        for (int i = 1; i <= N; ++i)
            for (int j = i; j <= N; j += i) ++dc[j], ds[j] += i;
        for (int p = 2; p <= N; ++p)
            if (om[p] == 0 && dc[p] == 2)
                for (int j = p; j <= N; j += p) ++om[j];
        SIEVE_CHECK(divisor_count_table(N) == dc && divisor_sum_table(N) == ds && omega_table(N) == om);
        cases += 3;
    }
}

static void test_segment(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 60 * TEST_SCALE; ++it) {
        long long l;  // √r までの篩を毎回作るので、大きい l（10^14 付近）ははじめの 24 回だけ
        switch (it < 24 ? it % 8 : it % 6) {
            case 0: case 1: case 2: l = rnd(rng, -5, 100); break;
            case 3: case 4: case 5: l = rnd(rng, 1, 1000000000LL); break;
            case 6: l = rnd(rng, 1, 100000000000000LL); break;
            default: l = 100000000000000LL - rnd(rng, 0, 3000); break;
        }
        const long long r = l + rnd(rng, 0, 2000);
        std::vector<long long> expect;
        for (long long x = std::max(l, 0LL); x < r; ++x)
            if (is_prime(x)) expect.push_back(x);
        SIEVE_CHECK(segment_primes(l, r) == expect && segment_primes_closed(l, r - 1) == expect);
        long long next = std::max(l, 1LL);
        segment_factorize(l, r, [&](long long x, const std::vector<PrimePower>& fs) {
            SIEVE_CHECK(x == next && fs == factorize(x));
            ++next, ++cases;
        });
        SIEVE_CHECK(next == std::max(r, std::max(l, 1LL)));
    }
    // ブロック（2^18 個）をまたぐ
    long long next = 1;
    segment_factorize_closed(1, 600000, [&](long long x, const std::vector<PrimePower>& fs) {
        SIEVE_CHECK(x == next && fs == factorize(x));
        ++next, ++cases;
    });
    SIEVE_CHECK(next == 600001);
    next = 99999999700000LL;  // 10^14 の手前（√r ≈ 10^7）
    segment_factorize(99999999700000LL, 100000000000000LL, [&](long long x, const std::vector<PrimePower>& fs) {
        SIEVE_CHECK(x == next++);
        if (x % 97 == 0) SIEVE_CHECK(fs == factorize(x));
        Int128 p = 1;
        for (const PrimePower& f : fs) {
            for (int k = 0; k < f.e; ++k) p *= f.p;
        }
        SIEVE_CHECK(p == x);
        ++cases;
    });
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_sieve(rng, cases);
    test_tables(rng, cases);
    test_segment(rng, cases);
    std::cout << "[sieve] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[sieve] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    const Sieve S(10000000);
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (int x = 1; x <= 10000000; ++x) s1 += (long long)S.factorize(x).size();
    const long long t_fac = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<long long> phi = phi_table(10000000);
    const long long t_phi = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<int> om = omega_table(10000000);
    const long long t_om = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const size_t np = primes_up_to(100000000).size();
    const long long t_pr = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const size_t ns = segment_primes(100000000000000LL, 100000000000000LL + 10000000).size();
    const long long t_seg = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    segment_factorize(100000000000000LL, 100000000000000LL + 1000000, [&](long long, const std::vector<PrimePower>& fs) { s2 += (long long)fs.size(); });
    const long long t_segf = ms(t0);
    std::cout << "[sieve] benchmark: Sieve(10^7) " << t_build << " ms, S.factorize を 1〜10^7 全部 " << t_fac << " ms, phi_table(10^7) " << t_phi
              << " ms, omega_table(10^7) " << t_om << " ms\n"
              << "[sieve] benchmark: primes_up_to(10^8)（" << np << " 個）" << t_pr << " ms, segment_primes（10^14 から幅 10^7。" << ns << " 個）" << t_seg
              << " ms, segment_factorize（10^14 から幅 10^6）" << t_segf << " ms  (check " << (s1 + s2 + phi[12345] + om[30]) % 7 << ")\n";
}

}  // namespace sieve_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace sieve_demo {

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
#define SHOW(expr, comment) sieve_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    sieve_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("篩（sieve） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Sieve S(N)                          N 以下の最小の素因数の表（N ≤ 5×10^7 程度）\n"
              << "    is_prime / spf / primes / factorize / prime_factors / divisors / count_divisors / sum_divisors\n"
              << "    phi / mu / omega / big_omega / largest_prime_factor / square_free_decompose / is_square_free\n"
              << "    prime_count(x) / count_primes(l, r) / count_primes_closed(L, R) / kth_prime(k) / lcm_mod(vs, m)   （N を超える数は 02_factorize に任せる）\n"
              << "  phi_table / mu_table / divisor_count_table / divisor_sum_table / omega_table   0〜N の表（N ≤ 10^8）\n"
              << "  primes_up_to(N) / primes_up_to<int>(N)   N 以下の素数 vector<long long> / メモリ半分の vector<int>（N < 2^31）\n"
              << "  segment_primes(l, r) / segment_factorize(l, r, f)（_closed も）   区間篩（r ≤ 10^14）\n";

    api("【1】 Sieve S(N);  S.is_prime(x) / S.factorize(x) / S.divisors(x)", "bool / vector<PrimePower> / vector<long long>",
        "最小の素因数の表を作り、N 以下の数の素数判定・素因数分解を O(log x) で。N を超える数を渡してもよい（ポラード・ロー）", "構築 O(N)、1 回 O(log x)");
    RUN(Sieve S(1000000));
    SHOW(S.is_prime(999983), "100 万以下で最大の素数");
    SHOW(S.factorize(720720), "2^4 · 3^2 · 5 · 7 · 11 · 13");
    SHOW(S.divisors(60), "昇順");
    SHOW(S.phi(36), "φ(36)");
    SHOW(S.prime_count(1000000), "100 万以下の素数の数");
    SHOW(S.count_primes_closed(100, 200), "[100, 200] の素数の数（count_primes(l, r) は [l, r)）");
    SHOW(S.kth_prime(0), "0 番目の素数");
    SHOW(S.factorize(1000000016000000063LL), "N を超えたら 02_factorize（ポラード・ロー）");
    {
        // ABC177 E: すべての組が互いに素（pairwise）/ 全体の gcd が 1（setwise）/ どちらでもない
        const std::vector<int> A = {6, 10, 15};
        std::vector<int> seen(1000001, 0);
        bool pairwise = true;
        for (int a : A)
            for (long long p : S.prime_factors(a)) pairwise &= ++seen[p] == 1;
        long long g = 0;
        for (int a : A) g = std::gcd(g, (long long)a);
        note("ABC177 E（Coprime）: 各素因数が 2 つ以上の数に出てこなければ pairwise。そうでなく全体の gcd が 1 なら setwise");
        std::cout << "    {6, 10, 15} → " << (pairwise ? "pairwise coprime" : g == 1 ? "setwise coprime" : "not coprime") << "\n";
    }
    {
        // ABC254 D: i·j が平方数になる 1 ≤ i, j ≤ N の組の数 = Σ_b (平方因子を除いた部分が b の数)^2
        const int N = 4;
        std::map<long long, long long> cnt;
        for (int i = 1; i <= N; ++i) cnt[S.square_free_decompose(i).b]++;
        long long ans = 0;
        for (auto& kv : cnt) ans += kv.second * kv.second;
        note("ABC254 D（Together Square）: 平方因子を除いた部分 b が同じ組なら積が平方数");
        std::cout << "    N = 4 → " << ans << "\n";
    }

    api("【2】 phi_table(N) / mu_table(N) / divisor_count_table(N) / divisor_sum_table(N) / omega_table(N)", "vector",
        "0〜N の値の表（添字がそのまま x）。たくさんの x の φ・μ などが要るとき", "O(N)（線形篩）");
    RUN(auto om = omega_table(15));
    SHOW(om, "異なる素因数の数");
    {
        // 典型90 030: N 以下で異なる素因数を K 個以上持つ数の個数
        const int N = 15, K = 2;
        int c = 0;
        for (int x = 1; x <= N; ++x) c += om[x] >= K;
        note("典型90 030（K Factors）: omega_table(N) で数える");
        std::cout << "    N = 15, K = 2 → " << c << "（6, 10, 12, 14, 15）\n";
    }
    SHOW(phi_table(10), "φ(0)〜φ(10)");
    SHOW(mu_table(10), "μ(0)〜μ(10)");

    api("【3】 primes_up_to(N) / primes_up_to<int>(N)", "vector<long long> / vector<int>",
        "N 以下の素数（昇順）。Sieve より速く省メモリ（10^8 で約 0.3 秒・素数 576 万個）。long long なので p * p も溢れない。<int> はメモリ半分", "O(N log log N)");
    RUN(auto ps = primes_up_to(30));
    SHOW(ps, "30 以下の素数");
    {
        // ABC383 D: N 以下で約数がちょうど 9 個の数 = p^8 または p^2 q^2（p < q）
        const long long N = 200;
        const std::vector<long long> qs = primes_up_to(isqrt(N));
        long long c = 0;
        for (size_t i = 0; i < qs.size(); ++i) {
            if (pow_cap(qs[i], 8, N + 1) <= N) ++c;
            for (size_t j = i + 1; j < qs.size() && qs[i] * qs[j] <= isqrt(N); ++j) ++c;
        }
        note("ABC383 D（9 Divisors）: p^8 か p²q²。√N 以下の素数で数える");
        std::cout << "    N = 200 → " << c << "（36, 100, 196）\n";
    }

    api("【4】 segment_primes(l, r) / segment_factorize(l, r, f)（_closed も）", "vector<long long> / なし",
        "[l, r) の素数 / 各数の素因数分解（r ≤ 10^14、幅 ≤ 10^7 程度。√r 以下の素数だけで割る）", "O((r − l) log log r + √r)");
    SHOW(segment_primes_closed(1000000000000LL, 1000000000100LL), "10^12 〜 10^12 + 100 の素数");
    {
        // ABC412 E: A_n = lcm(1, …, n)。A_L, …, A_R の値の種類 = 1 + (L, R] の素数の累乗の数
        const long long L = 4, R = 12;
        long long c = 1;
        segment_factorize_closed(L + 1, R, [&](long long, const std::vector<PrimePower>& fs) { c += fs.size() == 1; });
        note("ABC412 E（LCM Sequence）: lcm(1..n) は n が素数の累乗のときだけ変わる。(L, R] の素数の累乗を segment_factorize で数える");
        std::cout << "    L = 4, R = 12 → " << c << "（5, 7, 8, 9, 11 で変わる）\n";
    }
    note("x の小さい順に f(x, fs) が呼ばれる。fs は p の昇順の vector<PrimePower>（1 なら空）。");
}

}  // namespace sieve_demo

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
    sieve_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    sieve_test::benchmark();    // 速度計測
    sieve_demo::demo();         // 使い方デモ
    return 0;
}
