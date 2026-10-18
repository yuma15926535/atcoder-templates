// =============================================================================
//  素因数分解（1 つの大きな数の素数判定・素因数分解・約数・φ・μ、n! の素因数、累乗数の数など）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 素数・約数
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  使い分け：N 以下（N ≤ 5×10^7 程度）の数を何度も分解するなら 02_sieve の Sieve（このファイルを同梱している）。
//            1 つ・少数の大きな数（10^18 まで）ならこのファイル
//
//  機能一覧（n < 2^63）                                                       計算量
//    PrimePower{p, e}                      素因数分解の 1 項 p^e（cout << で "p^e"、e = 1 なら "p"）
//    is_prime(n)                           素数か（決定的ミラー・ラビン。誤判定なし）               O(log n)
//    factorize(n)                          素因数分解 vector<PrimePower>（p の昇順。n ≥ 1。1 なら空）  期待 O(n^{1/4})
//    factorize_trial(n)                    試し割りの素因数分解（確かめ用）                          O(√n)
//    prime_factors(n)                      異なる素因数（昇順）
//    divisors_of(n)                        約数（昇順）
//    divisors_from_factors(fs) / for_each_divisor(fs, f)   素因数分解から約数（昇順）/ 約数を f(d) に渡す
//    count_divisors(n) / sum_divisors(n)   約数の個数 / 約数の和（Int128）
//    euler_phi(n) / mobius(n)              オイラーの φ / メビウスの μ
//    omega(n) / big_omega(n)               異なる素因数の個数 ω / 重複を含めた素因数の個数 Ω（12 → 2 / 3）
//    largest_prime_factor(n) / radical(n)  最大の素因数 / 異なる素因数の積
//    square_free_decompose(n)              SquareFreeParts{a, b}：n = a^2 · b（b は平方因子なし）
//    is_square_free(n) / is_prime_power(n) 平方因子を持たないか / 素数の累乗か
//    factorial_exponent(n, p)              n! が素数 p で何回割れるか（ルジャンドル）                O(log n)
//    factorial_max_power(n, k)             k^e が n! を割り切る最大の e（k は合成数でもよい）
//    min_factorial_multiple(k)             k が n! を割り切る最小の n
//    is_perfect_power(n)                   a^b（a ≥ 1, b ≥ 2）と書けるか（1 は true）
//    count_perfect_powers(N)               1 以上 N 以下の a^b（b ≥ 2）の形の数の個数（1 も数える）   O(N^{1/3} log)
//    count_coprime_up_to(N, M)             1 以上 N 以下で M と互いに素な数の個数（包除）           O(2^ω(M))
//    lcm_mod_of(vs, m)                     多くの数の lcm を mod m で（mod を取りながら lcm を取ると誤るため。vector<int> も可）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    fac_test::random_test()  … 試し割り・篩・作り方の分かっている数との突き合わせ
//    fac_test::benchmark()    … 速度計測
//    fac_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace fac_test {

#define FAC_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[factorize] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 愚直な篩（N 以下の素数判定）
static std::vector<char> brute_sieve(int N) {
    std::vector<char> ip(N + 1, 1);
    ip[0] = 0;
    if (N >= 1) ip[1] = 0;
    for (int i = 2; (long long)i * i <= N; ++i)
        if (ip[i])
            for (int j = i * i; j <= N; j += i) ip[j] = 0;
    return ip;
}
static long long prod_of(const std::vector<PrimePower>& fs) {
    Int128 r = 1;
    for (const PrimePower& f : fs)
        for (int k = 0; k < f.e; ++k) r *= f.p;
    return (long long)r;
}
static bool well_formed(const std::vector<PrimePower>& fs) {  // p が昇順・重複なし・すべて素数・e ≥ 1
    for (size_t i = 0; i < fs.size(); ++i) {
        if (fs[i].e < 1 || !is_prime(fs[i].p)) return false;
        if (i && fs[i - 1].p >= fs[i].p) return false;
    }
    return true;
}
static long long rnd_prime(std::mt19937_64& rng, long long lo, long long hi) {
    while (true) {
        const long long x = rnd(rng, lo, hi);
        if (is_prime(x)) return x;
    }
}

static void test_is_prime(std::mt19937_64& rng, int& cases) {
    const int N = 2000000;
    const std::vector<char> ip = brute_sieve(N);
    for (int x = -5; x <= N; ++x, ++cases) FAC_CHECK(is_prime(x) == (x >= 0 && ip[x]));
    // 試し割りで確かめられる大きさ
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const long long x = rnd(rng, 1, 10000000000LL);
        bool pr = x >= 2;
        for (long long d = 2; d * d <= x && pr; ++d) pr = x % d != 0;
        FAC_CHECK(is_prime(x) == pr);
    }
    // 強擬素数・カーマイケル数（合成数）と、知られている大きな素数
    for (long long c : {561LL, 41041LL, 825265LL, 321197185LL, 5394826801LL, 232250619601LL, 9746347772161LL, 3215031751LL, 2152302898747LL,
                        3474749660383LL, 341550071728321LL, 3825123056546413051LL, 318665857834031151LL, 7999252175582851LL,
                        55245642489451LL, 1122004669633LL, 25326001LL, 4759123141LL, 1000000016000000063LL})
        FAC_CHECK(!is_prime(c));
    for (long long p : {2305843009213693951LL, 1000000000000000003LL, 1000000000000000009LL, 998244353LL, 1000000007LL,
                        4611686018427387847LL, 9223372036854775783LL, 4294967311LL, 4294967291LL, 2147483647LL})
        FAC_CHECK(is_prime(p));
    FAC_CHECK(!is_prime(9223372036854775807LL) && !is_prime(-7) && !is_prime(0) && !is_prime(1));
}

static void test_factorize(std::mt19937_64& rng, int& cases) {
    // 小さい数は全部、試し割りと
    for (long long n = 1; n <= 200000; ++n, ++cases) FAC_CHECK(factorize(n) == factorize_trial(n));
    for (int it = 0; it < 1000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 1, it % 4 ? 1000000 : 10000000000LL);
        FAC_CHECK(factorize(n) == factorize_trial(n));
    }
    // 大きい数は「作り方の分かっている数」と
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        std::vector<long long> ps;
        Int128 n = 1;
        const int kind = it % 5;
        if (kind == 0) {  // 2 つの大きい素数の積（ポラード・ローの苦手な形）
            const long long p = rnd_prime(rng, 1000000000LL, 3037000000LL), q = rnd_prime(rng, 2, 9223372036854775807LL / p);
            ps = {p, q};
        } else if (kind == 1) {  // 素数の累乗
            const long long p = rnd_prime(rng, 2, 3000000);
            for (Int128 x = p; x * p <= (Int128)9223372036854775807LL && (ps.size() < 2 || rng() % 3); x *= p) ps.push_back(p);
            if (ps.empty()) ps.push_back(p);
        } else {  // 大きさのばらばらな素数をいくつか
            while (true) {
                const int bits = (int)rnd(rng, 2, 40);
                const long long p = rnd_prime(rng, 2, (1LL << bits));
                if (n * p > (Int128)9223372036854775807LL) break;
                n *= p, ps.push_back(p);
                if (rng() % 4 == 0) break;
            }
            if (ps.empty()) ps.push_back(rnd_prime(rng, 2, 1000000000000000000LL));
        }
        n = 1;
        for (long long p : ps) n *= p;
        std::sort(ps.begin(), ps.end());
        std::vector<PrimePower> expect;
        for (long long p : ps) {
            if (!expect.empty() && expect.back().p == p) ++expect.back().e;
            else expect.push_back({p, 1});
        }
        FAC_CHECK(factorize((long long)n) == expect);
    }
    // ランダムな 64 ビットの数：積が戻ること・素数であること・昇順
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 1, 9223372036854775807LL);
        const std::vector<PrimePower> fs = factorize(n);
        FAC_CHECK(prod_of(fs) == n && well_formed(fs));
    }
    FAC_CHECK(factorize(1).empty() && factorize(9223372036854775807LL) == (std::vector<PrimePower>{{7, 2}, {73, 1}, {127, 1}, {337, 1}, {92737, 1}, {649657, 1}}));
    FAC_CHECK(factorize(1000000016000000063LL) == (std::vector<PrimePower>{{1000000007, 1}, {1000000009, 1}}));
}

static void test_divisor_functions(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 + 300 * TEST_SCALE; ++it, ++cases) {
        const long long n = it < 2000 ? it + 1 : rnd(rng, 1, it % 2 ? 1000000 : 10000000000LL);
        std::vector<long long> ds;
        for (long long d = 1; d * d <= n; ++d)
            if (n % d == 0) {
                ds.push_back(d);
                if (d != n / d) ds.push_back(n / d);
            }
        std::sort(ds.begin(), ds.end());
        FAC_CHECK(divisors_of(n) == ds && count_divisors(n) == (long long)ds.size());
        Int128 s = 0;
        for (long long d : ds) s += d;
        FAC_CHECK(sum_divisors(n) == s);
        // φ・μ・最大の素因数・平方因子・根基を試し割りの分解から
        const std::vector<PrimePower> fs = factorize_trial(n);
        long long phi = n, rad = 1, a = 1, b = 1;
        int mu = 1;
        for (const PrimePower& f : fs) {
            phi = phi / f.p * (f.p - 1), rad *= f.p;
            mu = f.e >= 2 ? 0 : -mu;
            for (int k = 0; k < f.e / 2; ++k) a *= f.p;
            if (f.e % 2) b *= f.p;
        }
        FAC_CHECK(euler_phi(n) == phi && mobius(n) == mu && radical(n) == rad && largest_prime_factor(n) == (fs.empty() ? 1 : fs.back().p));
        const SquareFreeParts sp = square_free_decompose(n);
        FAC_CHECK(sp.a == a && sp.b == b && is_square_free(n) == (a == 1) && is_prime_power(n) == (fs.size() == 1));
        FAC_CHECK(prime_factors(n).size() == fs.size() && omega(n) == (int)fs.size());
        {
            int bo = 0;
            for (const PrimePower& f : fs) bo += f.e;
            FAC_CHECK(big_omega(n) == bo);
        }
        // for_each_divisor は同じ集合
        std::vector<long long> d2;
        for_each_divisor(fs, [&](long long d) { d2.push_back(d); });
        std::sort(d2.begin(), d2.end());
        FAC_CHECK(d2 == ds && divisors_from_factors(fs) == ds);
        if (n <= 100000) {  // φ を数える
            long long c = 0;
            for (long long x = 1; x <= n; ++x) c += std::gcd(x, n) == 1;
            FAC_CHECK(c == phi);
        }
    }
}

static void test_factorial_and_powers(std::mt19937_64& rng, int& cases) {
    // ルジャンドル
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 3000), p = rnd_prime(rng, 2, 60);
        long long c = 0;
        for (long long i = 1; i <= n; ++i)
            for (long long m = i; m % p == 0; m /= p) ++c;
        FAC_CHECK(factorial_exponent(n, p) == c);
        const long long k = rnd(rng, 2, 200);
        long long e = 0;
        while (true) {  // k^(e+1) | n! か（素因数ごとに）
            bool ok = true;
            for (const PrimePower& f : factorize_trial(k)) ok &= factorial_exponent(n, f.p) >= (long long)f.e * (e + 1);
            if (!ok) break;
            ++e;
        }
        FAC_CHECK(factorial_max_power(n, k) == e);
    }
    for (long long k = 1; k <= 3000; ++k, ++cases) {  // k | n! となる最小の n
        long long n = 1;
        Int128 f = 1;
        while ((f % k) != 0) ++n, f = f * n % ((Int128)k * 1000000007);
        FAC_CHECK(min_factorial_multiple(k) == n);
    }
    FAC_CHECK(min_factorial_multiple(1000000007LL * 2) == 1000000007LL && min_factorial_multiple(1LL << 40) == 44);
    FAC_CHECK(factorial_exponent(1000000000000000000LL, 2) == 1000000000000000000LL - __builtin_popcountll(1000000000000000000ULL));
    // 累乗数
    {
        const long long LIM = 1000000;
        std::vector<char> pp(LIM + 1, 0);
        pp[1] = 1;
        for (long long a = 2; a * a <= LIM; ++a)
            for (long long x = a * a; x <= LIM; x *= a) pp[x] = 1;
        long long c = 0;
        for (long long x = 1; x <= LIM; ++x) {
            c += pp[x];
            FAC_CHECK(is_perfect_power(x) == (bool)pp[x]);
            if (x % 997 == 0 || x < 3000) FAC_CHECK(count_perfect_powers(x) == c);
            ++cases;
        }
    }
    // 大きい N はメビウスの式 1 − Σ_{k ≥ 2} μ(k)(⌊N^{1/k}⌋ − 1) と
    for (int it = 0; it < 4 * TEST_SCALE; ++it, ++cases) {
        const long long N = it == 0 ? 1000000000000000000LL : rnd(rng, 1, 9223372036854775807LL >> (rng() % 40));
        long long expect = 1;
        for (int k = 2; k < 63; ++k) expect -= mobius(k) * (kth_root_floor(N, k) - 1);
        FAC_CHECK(count_perfect_powers(N) == expect);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい累乗数
        const long long a = rnd(rng, 2, 3000000), b = rnd(rng, 2, 40);
        const long long x = pow_cap(a, b, 9223372036854775807LL);
        if (x == 9223372036854775807LL) continue;
        FAC_CHECK(is_perfect_power(x));
        if (b % 2 == 1 && !is_square(a)) FAC_CHECK(!is_square(x));
    }
}

static void test_misc(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long N = rnd(rng, 0, 800), M = it % 2 ? rnd(rng, 1, 3000) : rnd(rng, 1, 1000000000000000000LL);
        long long c = 0;
        for (long long x = 1; x <= N; ++x) c += std::gcd(x, M) == 1;
        FAC_CHECK(count_coprime_up_to(N, M) == c);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 6);
        std::vector<long long> vs(n);
        Int128 l = 1;
        bool big = false;
        for (auto& v : vs) {
            v = it % 2 ? rnd(rng, 1, 60) : rnd(rng, 1, 1000000000000LL);
            if (big) continue;
            const Int128 t = l / std::gcd((long long)(l % v), v);
            if (t > ((Int128)1 << 100) / v) big = true;
            else l = t * v;
        }
        const long long m = rnd(rng, 1, 1000000007);
        if (!big) FAC_CHECK(lcm_mod_of(vs, m) == (long long)(l % m));
    }
    FAC_CHECK(lcm_mod_of({0, 5}, 7) == 0 && lcm_mod_of({}, 7) == 1 && lcm_mod_of({4, 6}, 1) == 0 && lcm_mod_of(std::vector<int>{4, 6, 10}, 1000) == 60);
    FAC_CHECK(omega(1) == 0 && big_omega(1) == 0 && omega(720) == 3 && big_omega(720) == 7 && big_omega(1LL << 62) == 62);
    std::ostringstream os;
    for (const PrimePower& f : factorize(360)) os << f << ' ';
    FAC_CHECK(os.str() == "2^3 3^2 5 ");
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_is_prime(rng, cases);
    test_factorize(rng, cases);
    test_divisor_functions(rng, cases);
    test_factorial_and_powers(rng, cases);
    test_misc(rng, cases);
    std::cout << "[factorize] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[factorize] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 100000;
    std::vector<long long> xs(n);
    for (auto& x : xs) x = rnd(rng, 1, 1000000000000000000LL) | 1;
    auto t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (long long x : xs) s1 += is_prime(x);
    const long long t_pr = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for (int i = 0; i < 10000; ++i) s2 += (long long)factorize(xs[i]).size();
    const long long t_fac = ms(t0);
    std::vector<long long> semi(1000);
    for (auto& x : semi) x = rnd_prime(rng, 1000000000LL, 3000000000LL) * rnd_prime(rng, 1000000000LL, 3000000000LL);
    t0 = std::chrono::steady_clock::now();
    long long s3 = 0;
    for (long long x : semi) s3 += factorize(x)[0].p;
    const long long t_semi = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const long long s4 = count_perfect_powers(1000000000000000000LL);
    const long long t_pp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const long long s5 = (long long)divisors_of(963761198400LL).size();
    const long long t_div = ms(t0);
    std::cout << "[factorize] benchmark: is_prime（10^18 程度）× " << n << " " << t_pr << " ms, factorize（ランダム 10^18）× 10^4 " << t_fac
              << " ms, factorize（約 30 ビットの素数 2 つの積）× 1000 " << t_semi << " ms\n"
              << "[factorize] benchmark: count_perfect_powers(10^18) " << t_pp << " ms, divisors_of(963761198400)（約数 6720 個）" << t_div
              << " ms  (check " << (s1 + s2 + s3 + s4 + s5) % 7 << ")\n";
}

}  // namespace fac_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace fac_demo {

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
#define SHOW(expr, comment) fac_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    fac_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("素因数分解（factorize） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧（n < 2^63）\n"
              << "  is_prime(n) / factorize(n) / prime_factors(n)        素数判定 / 素因数分解（{p, e} の昇順）/ 異なる素因数\n"
              << "  divisors_of(n) / divisors_from_factors(fs) / for_each_divisor(fs, f)   約数\n"
              << "  count_divisors / sum_divisors / euler_phi / mobius / omega / big_omega / largest_prime_factor / radical\n"
              << "  square_free_decompose(n) / is_square_free(n) / is_prime_power(n)          n = a²·b、平方因子なし、素数の累乗\n"
              << "  factorial_exponent(n, p) / factorial_max_power(n, k) / min_factorial_multiple(k)   n! の素因数\n"
              << "  is_perfect_power(n) / count_perfect_powers(N)        a^b の形の数\n"
              << "  count_coprime_up_to(N, M) / lcm_mod_of(vs, m)        M と互いに素な数の個数 / 多くの数の lcm を mod で\n"
              << "  （N 以下の数を何度も分解するなら 02_sieve の Sieve。このファイルの関数も使える）\n";

    api("【1】 is_prime(n) / factorize(n)", "bool / vector<PrimePower>",
        "10^18 の数でも速い（ミラー・ラビンとポラード・ロー）。factorize は p の昇順の {p, e}（PrimePower）。cout << で p^e と出る", "O(log n) / 期待 O(n^{1/4})");
    SHOW(is_prime(1000000007), "素数");
    SHOW(is_prime(1000000016000000063LL), "10^9+7 と 10^9+9 の積");
    SHOW(factorize(360), "2^3 · 3^2 · 5");
    SHOW(factorize(1000000016000000063LL), "大きい素数 2 つの積も速い");
    RUN(for (auto [p, e] : factorize(720)) std::cout << "  " << p << "^" << e; std::cout << "\n");
    {
        // ABC169 D: N = p1^e1 ... に、z = p^k（使っていないもの）で割る操作を何回できるか。各 e について 1 + 2 + … + k ≤ e となる最大の k
        const long long N = 24;
        long long ans = 0;
        for (auto [p, e] : factorize(N))
            for (long long k = 1, rest = e; rest >= k; rest -= k, ++k) ++ans;
        note("ABC169 D（Div Game）: 各素因数の指数 e について 1 + 2 + … ≤ e となる回数を足す");
        std::cout << "    N = 24 → " << ans << "\n";
    }
    {
        // ABC284 D: N = p^2 q（N ≤ 9×10^18）の p と q
        const long long N = 2023;
        auto fs = factorize(N);
        const long long p = fs[0].e == 2 ? fs[0].p : fs[1].p, q = fs[0].e == 2 ? fs[1].p : fs[0].p;
        note("ABC284 D（N = p²q）: factorize(N) で分解して指数 2 の方が p");
        std::cout << "    N = 2023 → p = " << p << ", q = " << q << "\n";
    }

    api("【2】 divisors_of(n) / count_divisors(n) / sum_divisors(n)", "vector<long long> / long long / Int128",
        "約数（昇順）・約数の個数・約数の和。素因数分解が分かっているなら divisors_from_factors(fs)", "O(約数の個数 · log)");
    SHOW(divisors_of(36), "昇順");
    SHOW(count_divisors(963761198400LL), "10^12 以下で約数が最も多い数の 1 つ");
    SHOW(sum_divisors(28), "完全数（28 + 28）");
    {
        // ABC190 D: 和が N となる連続する整数の列の数 = 2 × (N の奇数の約数の個数)
        long long N = 12, odd = N;
        while (odd % 2 == 0) odd /= 2;
        note("ABC190 D（和が N の等差数列の数）: 2 × count_divisors(N から 2 を除いた数)");
        std::cout << "    N = 12 → " << 2 * count_divisors(odd) << "\n";
    }

    api("【3】 euler_phi(n) / mobius(n) / omega(n) / big_omega(n) / largest_prime_factor(n) / radical(n) / is_prime_power(n)", "long long / int / bool",
        "φ(n)（n 以下で n と互いに素な数の個数）・μ(n)・最大の素因数・異なる素因数の積・素数の累乗か", "factorize と同じ");
    SHOW(euler_phi(36), "36 と互いに素な 1〜36 の数");
    SHOW(mobius(30), "2·3·5 → (-1)^3");
    SHOW(omega(720), "720 = 2^4·3^2·5：異なる素因数は 3 個");
    SHOW(big_omega(720), "重複を含めると 4 + 2 + 1 = 7 個");
    SHOW(mobius(12), "平方因子があれば 0");
    SHOW(largest_prime_factor(1001), "7·11·13");
    SHOW(is_prime_power(243), "3^5");

    api("【4】 square_free_decompose(n) / is_square_free(n)", "SquareFreeParts{a, b}", "n = a²·b（b は平方因子を持たない）。「掛けて平方数になる組」は b が同じ組", "factorize と同じ");
    RUN(auto [a, b] = square_free_decompose(72));
    SHOW(a, "72 = 6² · 2");
    SHOW(b, "平方因子を除いた部分");
    note("ABC254 D（i·j が平方数になる組）: i·j が平方数 ⇔ i と j の b が同じ（N 以下を全部なら 02_sieve の Sieve::square_free_decompose）");

    api("【5】 factorial_exponent(n, p) / factorial_max_power(n, k) / min_factorial_multiple(k)", "long long",
        "n! が素数 p で何回割れるか（ルジャンドル）/ k^e | n! の最大の e / k | n! となる最小の n", "O(log n) / factorize(k)");
    SHOW(factorial_exponent(100, 5), "100! の末尾の 0 の数（5 の方が少ない）");
    SHOW(factorial_max_power(10, 12), "12^e | 10! の最大の e（2^8 · 3^4 → 4）");
    SHOW(min_factorial_multiple(30), "ABC280 D: 30 | 5!");
    SHOW(min_factorial_multiple(123456789011LL), "大きい素数なら k そのもの");
    note("ABC052 C（N! の約数の個数）: 各素数 p ≤ N について factorial_exponent(N, p) + 1 を掛ける");

    api("【6】 is_perfect_power(n) / count_perfect_powers(N)", "bool / long long", "a^b（a ≥ 1, b ≥ 2）の形か / 1〜N にいくつあるか（1 も数える）", "O(N^{1/3} log)");
    SHOW(is_perfect_power(1000000), "10^6 = 1000^2 = 100^3");
    SHOW(count_perfect_powers(99), "ABC361 F: 1, 4, 8, 9, 16, 25, 27, 32, 36, 49, 64, 81");
    SHOW(count_perfect_powers(1000000000000000000LL), "N = 10^18 でも速い");

    api("【7】 count_coprime_up_to(N, M) / lcm_mod_of(vs, m)", "long long", "1〜N で M と互いに素な数の個数（包除）/ 多くの数の lcm を mod で", "O(2^(M の素因数の数)) / 分解");
    SHOW(count_coprime_up_to(100, 30), "2・3・5 のどれでも割り切れない数");
    SHOW(lcm_mod_of({1000000007, 1000000009, 6}, 998244353), "lcm は約 6×10^18。mod を取りながら lcm を取ると誤る");
    note("N 以下の数をたくさん分解するときは 02_sieve の Sieve S(N)（このファイルの関数も一緒に入っている）。");
}

}  // namespace fac_demo

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
    fac_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    fac_test::benchmark();    // 速度計測
    fac_demo::demo();         // 使い方デモ
    return 0;
}
