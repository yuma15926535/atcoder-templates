// =============================================================================
//  順列・組合せの番号（辞書順で何番目か・k 番目の並び・重複のある並び・転倒数・置換の偶奇）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 04 組合せ
//  依存       : 00_int_common, 00_math_base, 00_modint, 04_binomial（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  番号（何番目か・k 番目）はすべて 0 から。列は vector でも string でもよい
//  厳密な番号は Int128（cout << できる。収まらなければ止める）。大きいときは _mod（using mint = ModInt998; など）
//
//  機能一覧                                                                    計算量
//    inversion_count(a)               転倒数：i < j かつ a_i > a_j の組の数（値は比べられれば何でも）     O(N log N)
//    perm_parity(a)                   置換の偶奇（偶 0 / 奇 1。値がすべて異なる列。サイクルの数から）        O(N)〜O(N log N)
//    perm_rank(p)                     値がすべて異なる列 p が、その並べ替え全体の辞書順で何番目か（n ≤ 33）  O(N log N)
//    perm_rank_mod<mint>(p)           同じく mod で（n はいくつでも）
//    perm_unrank(n, k)                0〜n−1 の順列で辞書順 k 番目のもの（0 ≤ k < n!）                    O(N log N)
//    multiset_perm_rank(a)            同じ値があってよい列 a が、その並べ替え（同じものは区別しない）の何番目か  O(N log N)
//    multiset_perm_rank_mod<mint>(a)  同じく mod で（mod は n より大きい素数）
//    multiset_perm_unrank(a, k)       a の並べ替えで辞書順 k 番目のもの（"aabb" の 3 番目 → "abba"）        O(N × 値の種類)
//    kth_combination(n, r, k)         {0..n−1} から r 個選ぶ組（小さい順の列）で辞書順 k 番目                O(n × min(r, 130))
//    combination_rank(n, c)           組 c が辞書順で何番目か
//  （次・前の順列は std::next_permutation / std::prev_permutation。並べ方の数・nCr は 04_binomial）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    perm_test::random_test()  … next_permutation で全部並べた愚直解との突き合わせ
//    perm_test::benchmark()    … 速度計測
//    perm_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <initializer_list>
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
#ifndef GG_BINOMIAL
#define GG_BINOMIAL
// ---- 二項係数（数学第二世代）  Binomial<mint>（階乗の表）・n が大きい nCr・Lucas・合成数の mod・厳密な値・パスカルの三角形・Σ_{i≤m} C(n, i)（Mo 用） ----
// Binomial<mint> B(N)：0〜N の階乗・階乗の逆元の表（mint の mod は素数。N を超える n が来たら自動で表を伸ばす）
template <class mint>
struct Binomial {
    explicit Binomial(long long n = 0) { ensure(std::max(n, 1LL)); }
    mint fact(long long n) const {  // n!
        if (n < 0) math_internal::fail("Binomial::fact", "n ≥ 0 で使う");
        ensure(n);
        return f_[(size_t)n];
    }
    mint inv_fact(long long n) const {  // 1 / n!
        if (n < 0) math_internal::fail("Binomial::inv_fact", "n ≥ 0 で使う");
        ensure(n);
        return fi_[(size_t)n];
    }
    mint inv(long long n) const {  // 1 / n（n ≥ 1）
        if (n < 1) math_internal::fail("Binomial::inv", "n ≥ 1 で使う");
        ensure(n);
        return fi_[(size_t)n] * f_[(size_t)n - 1];
    }
    mint C(long long n, long long r) const {  // nCr（r < 0・r > n・n < 0 なら 0）
        if (r < 0 || r > n) return mint(0);
        ensure(n);
        return f_[(size_t)n] * fi_[(size_t)r] * fi_[(size_t)(n - r)];
    }
    mint P(long long n, long long r) const {  // nPr = n! / (n − r)!（範囲外は 0）
        if (r < 0 || r > n) return mint(0);
        ensure(n);
        return f_[(size_t)n] * fi_[(size_t)(n - r)];
    }
    mint H(long long n, long long r) const {  // nHr = (n + r − 1)Cr（n 種類から重複を許して r 個。H(n, 0) = 1）
        if (n < 0 || r < 0) return mint(0);
        if (r == 0) return mint(1);
        return C(n + r - 1, r);
    }
    template <class V>
    mint multinomial(const V& ks) const {  // (Σ k)! / Π k!（負の k があれば 0）
        long long s = 0;
        for (auto kv : ks) {
            const long long k = math_internal::to_ll_checked(kv, "Binomial::multinomial");
            if (k < 0) return mint(0);
            if (__builtin_add_overflow(s, k, &s)) math_internal::fail("Binomial::multinomial", "個数の和が long long に収まらない");
        }
        ensure(s);
        mint r = f_[(size_t)s];
        for (auto k : ks) r *= fi_[(size_t)k];
        return r;
    }
    mint multinomial(std::initializer_list<long long> ks) const { return multinomial<std::initializer_list<long long>>(ks); }
    mint catalan(long long n) const {  // カタラン数 C(2n, n) / (n + 1)
        if (n < 0) return mint(0);
        return C(2 * n, n) - C(2 * n, n + 1);
    }
    mint paths(long long x, long long y) const { return x < 0 || y < 0 ? mint(0) : C(x + y, x); }  // (0, 0) → (x, y) の右・上だけの道
    // (0, 0) → (x, y) の右・上だけの道で、ずっと y ≤ x + k のもの（k ≥ 0。鏡像法：C(x + y, x) − C(x + y, y − k − 1)）
    mint paths_not_crossing(long long x, long long y, long long k) const {
        if (x < 0 || y < 0 || k < 0 || y > x + k) return mint(0);
        return C(x + y, x) - C(x + y, y - k - 1);
    }

  private:
    mutable std::vector<mint> f_, fi_;
    mutable int built_mod_ = 0;
    void ensure(long long n) const {
        if (n >= (long long)f_.size() || built_mod_ != mint::mod()) grow(n);
    }
    void grow(long long n) const {
        const long long M = mint::mod();
        if (built_mod_ != M) {  // はじめて、または DynModInt の mod が変わった
            if (!math_internal::is_prime_u32((unsigned)M)) math_internal::fail("Binomial", "mod が素数でない（合成数の mod は pascal_table か binom_mod_any）");
            f_.clear(), fi_.clear(), built_mod_ = (int)M;
        }
        if (n >= M) math_internal::fail("Binomial", "n ≥ mod では n! が 0 になる（binom_lucas を使う）");
        if (n > 100000000) math_internal::fail("Binomial", "n が大きすぎる（表は 10^8 まで。r が小さいなら binom_large_n）");
        const size_t old = f_.size();
        if ((long long)old > n) return;
        size_t sz = std::max<size_t>((size_t)n + 1, std::min<size_t>(old + old / 2, (size_t)M));  // 1.5 倍ずつ伸ばす
        sz = std::min<size_t>(sz, (size_t)100000001);
        f_.resize(sz), fi_.resize(sz);
        if (old == 0) f_[0] = 1;
        for (size_t i = std::max<size_t>(old, 1); i < sz; ++i) f_[i] = f_[i - 1] * mint((long long)i);
        fi_[sz - 1] = f_[sz - 1].inv();
        for (size_t i = sz - 1; i > old; --i) fi_[i - 1] = fi_[i] * mint((long long)i);
    }
};

// S(n, m) = Σ_{i=0}^{m} C(n, i)（n 個から m 個以下を選ぶ数）を、n・m を 1 ずつ動かしながら O(1) で保つ（Mo のアルゴリズムで多くのクエリに）
//   n ≥ 0。m < 0 なら 0、m ≥ n なら 2^n。mod は 3 以上の素数（n を減らすとき 2 で割る）
template <class mint>
struct BinomialPrefix {
    explicit BinomialPrefix(long long max_n = 0) : B_(std::max(max_n, 1LL)) {
        if (mint::mod() == 2) math_internal::fail("BinomialPrefix", "mod = 2 では 2 で割れない");
    }
    long long n() const { return n_; }
    long long m() const { return m_; }
    mint value() const { return s_; }  // Σ_{i=0}^{m} C(n, i)
    void inc_n() { s_ = s_ + s_ - B_.C(n_, m_), ++n_; }  // S(n+1, m) = 2 S(n, m) − C(n, m)
    void dec_n() {                                       // S(n−1, m) = (S(n, m) + C(n−1, m)) / 2
        if (n_ == 0) math_internal::fail("BinomialPrefix::dec_n", "n ≥ 0 のまま使う");
        --n_, s_ = (s_ + B_.C(n_, m_)) * inv2();
    }
    void inc_m() { ++m_, s_ += B_.C(n_, m_); }
    void dec_m() { s_ -= B_.C(n_, m_), --m_; }
    mint move_to(long long n, long long m) {  // (n, m) まで 1 ずつ動かして値を返す（|Δn| + |Δm| 回。m ≥ n なら値は 2^n なので、m は min(m, n) にして渡すと速い）
        if (n < 0) math_internal::fail("BinomialPrefix::move_to", "n ≥ 0 で使う");
        while (n_ < n) inc_n();
        while (m_ > m) dec_m();
        while (n_ > n) dec_n();
        while (m_ < m) inc_m();
        return s_;
    }

  private:
    Binomial<mint> B_;
    long long n_ = 0, m_ = 0;
    mint s_ = mint(1);  // S(0, 0) = 1
    mint inv2() const {
        static long long cached_mod = -1;
        static mint v;
        if (cached_mod != (long long)mint::mod()) cached_mod = mint::mod(), v = mint(2).inv();
        return v;
    }
};

// n が大きく r が小さい nCr（n ≤ 10^18、r ≤ 10^7 程度、r < mod。mod は素数）  O(r)
template <class mint>
mint binom_large_n(long long n, long long r) {
    if (r < 0 || r > n) return mint(0);
    if (r >= mint::mod()) math_internal::fail("binom_large_n", "r ≥ mod（binom_lucas を使う）");
    mint num = 1, den = 1;
    for (long long i = 0; i < r; ++i) num *= mint(n - i), den *= mint(i + 1);
    return num / den;
}

// nCr mod p（p は素数で 10^7 以下。n, r ≤ 10^18。Lucas の定理）  前処理 O(p)（同じ p なら 1 回）・O(log_p n)
inline long long binom_lucas(long long n, long long r, long long p) {
    if (p < 2 || p > 10000000 || !math_internal::is_prime_u32((unsigned)p)) math_internal::fail("binom_lucas", "p は 10^7 以下の素数");
    if (r < 0 || r > n) return 0;
    static long long cached = 0;
    static std::vector<int> f, fi;
    if (cached != p) {
        f.assign((size_t)p, 1), fi.assign((size_t)p, 1);
        for (long long i = 1; i < p; ++i) f[(size_t)i] = (int)(f[(size_t)i - 1] * i % p);
        fi[(size_t)p - 1] = (int)mod_pow(f[(size_t)p - 1], p - 2, p);
        for (long long i = p - 1; i > 0; --i) fi[(size_t)i - 1] = (int)(fi[(size_t)i] * i % p);
        cached = p;
    }
    long long res = 1;
    while (r > 0) {
        const long long ni = n % p, ri = r % p;
        if (ri > ni) return 0;
        res = res * f[(size_t)ni] % p * fi[(size_t)ri] % p * fi[(size_t)(ni - ri)] % p;
        n /= p, r /= p;
    }
    return res;
}

namespace math_internal {
struct BinomPrimePower {  // nCr mod p^e のための表（p で割れる数を除いた積）
    long long p, pe;
    int e;
    std::vector<int> prod;  // prod[i] = Π_{j ≤ i, p ∤ j} j mod p^e
    long long fact_wo_p(long long n) const {  // n! から p をすべて除いたもの mod p^e
        long long r = 1;
        while (n > 0) {
            r = r * mod_pow(prod[(size_t)pe], n / pe, pe) % pe * prod[(size_t)(n % pe)] % pe;
            n /= p;
        }
        return r;
    }
    long long legendre(long long n) const {  // n! の中の p の個数
        long long c = 0;
        while (n > 0) n /= p, c += n;
        return c;
    }
    long long binom(long long n, long long r) const {
        const long long k = legendre(n) - legendre(r) - legendre(n - r);
        if (k >= e) return 0;
        long long v = fact_wo_p(n) * mod_inv(fact_wo_p(r) * fact_wo_p(n - r) % pe, pe) % pe;
        for (long long i = 0; i < k; ++i) v = v * p % pe;
        return v;
    }
};
}  // namespace math_internal
// nCr mod m（m は合成数でもよい。m の素数べきの部分はそれぞれ 10^7 以下。n, r ≤ 10^18）  前処理 O(Σ p^e)（同じ m なら 1 回）
inline long long binom_mod_any(long long n, long long r, long long m) {
    if (m < 1) math_internal::fail("binom_mod_any", "m ≥ 1 で使う");
    if (r < 0 || r > n) return 0;
    if (m == 1) return 0;
    static long long cached = 0;
    static std::vector<math_internal::BinomPrimePower> parts;
    if (cached != m) {
        parts.clear();
        long long x = m;
        for (long long p = 2; p <= 10000000 && p * p <= x; ++p) {  // 10^7 より大きい素因数は下で止める
            if (x % p) continue;
            math_internal::BinomPrimePower b{p, 1, 0, {}};
            while (x % p == 0) x /= p, b.pe *= p, b.e++;
            parts.push_back(b);
        }
        if (x > 1) parts.push_back(math_internal::BinomPrimePower{x, x, 1, {}});
        for (auto& b : parts) {
            if (b.pe > 10000000) math_internal::fail("binom_mod_any", "m の素数べきの部分が 10^7 を超える");
            b.prod.assign((size_t)b.pe + 1, 1);
            for (long long i = 1; i <= b.pe; ++i) b.prod[(size_t)i] = (int)(i % b.p == 0 ? b.prod[(size_t)i - 1] : b.prod[(size_t)i - 1] * i % b.pe);
        }
        cached = m;
    }
    long long res = 0, M = 1;  // 中国剰余定理で合わせる（素数べきは互いに素）
    for (const auto& b : parts) {
        const long long a = b.binom(n, r);
        const long long t = mod_mul(a - res % b.pe, mod_inv(M % b.pe, b.pe), b.pe);
        res += M * t, M *= b.pe;
    }
    return res;
}

// 厳密な nCr（long long に収まらなければ止める）/ min(nCr, cap)（溢れない）  O(min(r, n − r))（値が大きくなる前に終わる）
namespace math_internal {
inline Int128 binom_sat(long long n, long long r, Int128 cap) {  // min(nCr, cap)（0 ≤ r ≤ n、cap ≤ 2^63）
    r = std::min(r, n - r);
    Int128 v = 1;
    for (long long i = 1; i <= r; ++i) {  // v = C(n − r + i, i)。i が増えると増えるので、cap を超えたら終わり
        v = v * (n - r + i) / i;
        if (v > cap) return cap;
    }
    return std::min(v, cap);  // r = 0 で cap = 0 のとき
}
}  // namespace math_internal
inline long long binom_exact(long long n, long long r) {
    if (r < 0 || r > n) return 0;
    const Int128 L = (Int128)9223372036854775807LL;
    const Int128 v = math_internal::binom_sat(n, r, L + 1);
    if (v > L) math_internal::fail("binom_exact", "nCr が long long に収まらない（binom_cap か mod で）");
    return (long long)v;
}
inline long long binom_cap(long long n, long long r, long long cap) {
    if (cap < 0) math_internal::fail("binom_cap", "cap ≥ 0 で使う");
    if (r < 0 || r > n) return 0;
    return (long long)math_internal::binom_sat(n, r, cap);
}

// パスカルの三角形 t[n][r] = nCr（0 ≤ r ≤ n ≤ N）。T は何でも（mint・DynModInt・double・long long は溢れたら止める）  O(N^2)
template <class T = long long>
std::vector<std::vector<T>> pascal_table(int N) {
    if (N < 0) math_internal::fail("pascal_table", "N ≥ 0 で使う");
    std::vector<std::vector<T>> t((size_t)N + 1);
    for (int n = 0; n <= N; ++n) {
        t[(size_t)n].assign((size_t)n + 1, T(1));
        for (int r = 1; r < n; ++r) {
            if constexpr (std::is_integral_v<T>) {
                if (__builtin_add_overflow(t[(size_t)n - 1][(size_t)r - 1], t[(size_t)n - 1][(size_t)r], &t[(size_t)n][(size_t)r]))
                    math_internal::fail("pascal_table", "値が型に収まらない（long long なら N ≤ 66）");
            } else {
                t[(size_t)n][(size_t)r] = t[(size_t)n - 1][(size_t)r - 1] + t[(size_t)n - 1][(size_t)r];
            }
        }
    }
    return t;
}
#endif  // GG_BINOMIAL
#ifndef GG_PERMUTATION
#define GG_PERMUTATION
// ---- 順列・組合せの番号（数学第二世代）  辞書順で何番目か・k 番目の並び（番号は 0 から）・転倒数・置換の偶奇 ----
// 列は vector でも string でもよい。厳密な番号は Int128（収まらなければ止める。そのときは _mod）
namespace math_internal {
struct PermFenwick {  // 0/1 の出現数を数える BIT
    explicit PermFenwick(int n) : n_(n), t_((size_t)n + 1, 0) {}
    void add(int i, int v) {
        for (++i; i <= n_; i += i & -i) t_[(size_t)i] += v;
    }
    int prefix(int i) const {  // [0, i) の和
        int s = 0;
        for (; i > 0; i -= i & -i) s += t_[(size_t)i];
        return s;
    }
    int kth(int k) const {  // 和が k + 1 になる最小の位置（k 番目の 1。0 から）
        int pos = 0, lg = 1;
        while (lg * 2 <= n_) lg *= 2;
        for (; lg > 0; lg /= 2)
            if (pos + lg <= n_ && t_[(size_t)(pos + lg)] <= k) pos += lg, k -= t_[(size_t)pos];
        return pos;
    }

  private:
    int n_;
    std::vector<int> t_;
};
inline Int128 gcd128(Int128 a, Int128 b) {  // a, b ≥ 0
    while (b) {
        const Int128 t = a % b;
        a = b, b = t;
    }
    return a;
}
const Int128 PERM_INF = (Int128)(((UInt128)1 << 127) - 1);  // Int128 の最大値（「とても大きい」の印にも使う）
// M × a / b（割り切れると分かっているとき）。PERM_INF を超えたら PERM_INF
inline Int128 mul_div_exact(Int128 M, Int128 a, Int128 b) {
    const Int128 g = gcd128(M, b);
    M /= g, b /= g;
    a /= b;  // gcd(M/g, b/g) = 1 なので b/g は a を割り切る
    if (a != 0 && M > PERM_INF / a) return PERM_INF;
    return M * a;
}
inline Int128 binom128(long long n, long long r) {  // min(nCr, PERM_INF)（PERM_INF は nCr になれないので「超えた」の印）
    if (r < 0 || r > n) return 0;
    r = std::min(r, n - r);
    Int128 v = 1;
    for (long long i = 1; i <= r; ++i) {  // v = C(n − r + i, i)：増えていくので PERM_INF になったら終わり
        v = mul_div_exact(v, n - r + i, i);
        if (v == PERM_INF) return v;
    }
    return v;
}
template <class C>
std::vector<int> perm_compress(const C& a, int& kinds) {  // 値を 0, 1, … に（同じ値は同じ番号）
    using T = typename std::decay<decltype(a[0])>::type;
    const size_t n = a.size();
    std::vector<int> r(n);
    kinds = 0;
    if (n == 0) return r;
    if constexpr (std::is_integral<T>::value) {  // 値の幅が狭い（1〜n など）なら数えるだけ  O(N)
        const auto mm = std::minmax_element(a.begin(), a.end());
        const T lo = *mm.first, hi = *mm.second;
        if ((UInt128)((Int128)hi - (Int128)lo) < (UInt128)2 * n + 64) {
            std::vector<int> id((size_t)((Int128)hi - (Int128)lo) + 1, 0);
            for (size_t i = 0; i < n; ++i) id[(size_t)((Int128)a[i] - (Int128)lo)] = 1;
            for (int& x : id) x = x ? kinds++ : 0;
            for (size_t i = 0; i < n; ++i) r[i] = id[(size_t)((Int128)a[i] - (Int128)lo)];
            return r;
        }
    }
    std::vector<std::pair<T, int>> v(n);
    for (size_t i = 0; i < n; ++i) v[i] = {a[i], (int)i};
    std::sort(v.begin(), v.end());
    for (size_t i = 0; i < n; ++i) {
        if (i && v[i - 1].first < v[i].first) ++kinds;
        r[(size_t)v[i].second] = kinds;
    }
    ++kinds;
    return r;
}
}  // namespace math_internal

// 転倒数：i < j かつ a_i > a_j の組の数（値は比べられれば何でも。同じ値は数えない）  O(N log N)
template <class C>
long long inversion_count(const C& a) {
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(a, k);
    math_internal::PermFenwick bit(k);
    long long s = 0;
    for (size_t i = 0; i < c.size(); ++i) {
        s += (long long)i - bit.prefix(c[i] + 1);
        bit.add(c[i], 1);
    }
    return s;
}
// 置換の偶奇：値がすべて異なる列 a をソートする置換が偶置換なら 0、奇置換なら 1（= 転倒数の偶奇 = (N − サイクル数) mod 2）  O(N log N)（値が 0〜N 程度なら O(N)）
template <class C>
int perm_parity(const C& a) {
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(a, k);
    if (k != (int)c.size()) math_internal::fail("perm_parity", "値がすべて異なる列で使う（同じ値があると偶奇が決まらない）");
    std::vector<char> seen(c.size(), 0);
    int cycles = 0;
    for (size_t i = 0; i < c.size(); ++i)
        if (!seen[i]) {
            ++cycles;
            for (size_t j = i; !seen[j]; j = (size_t)c[j]) seen[j] = 1;
        }
    return ((int)c.size() - cycles) % 2;
}

// ---- 順列（値がすべて異なる列）----
// p が「同じ値の並べ替え全体」の辞書順で何番目か（0 から）。値は何でもよい（1〜n でも 0〜n−1 でも）。n ≤ 33  O(N log N)
template <class C>
Int128 perm_rank(const C& p) {
    const int n = (int)p.size();
    if (n > 33) math_internal::fail("perm_rank", "n ≤ 33 で使う（番号が Int128 に収まらない。大きいときは perm_rank_mod）");
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(p, k);
    if (k != n) math_internal::fail("perm_rank", "値がすべて異なる列で使う（同じ値があるなら multiset_perm_rank）");
    math_internal::PermFenwick bit(n);
    Int128 r = 0;
    for (int i = 0; i < n; ++i) {
        const int less = c[(size_t)i] - bit.prefix(c[(size_t)i]);  // まだ使っていない、より小さい値の数
        r = r * (n - i) + less;                                     // ホーナー法：Σ less_i (n−1−i)!
        bit.add(c[(size_t)i], 1);
    }
    return r;
}
template <class mint, class C>
mint perm_rank_mod(const C& p) {  // perm_rank を mod で（n はいくつでも）  O(N log N)
    const int n = (int)p.size();
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(p, k);
    if (k != n) math_internal::fail("perm_rank_mod", "値がすべて異なる列で使う（同じ値があるなら multiset_perm_rank_mod）");
    math_internal::PermFenwick bit(n);
    mint r = 0;
    for (int i = 0; i < n; ++i) {
        r = r * mint(n - i) + mint(c[(size_t)i] - bit.prefix(c[(size_t)i]));
        bit.add(c[(size_t)i], 1);
    }
    return r;
}
// 0〜n−1 の順列で辞書順 k 番目（0 から）のもの（0 ≤ k < n!）  O(N log N)
inline std::vector<int> perm_unrank(int n, Int128 k) {
    if (n < 0 || k < 0) math_internal::fail("perm_unrank", "n ≥ 0、k ≥ 0 で使う");
    std::vector<Int128> f((size_t)n + 1, 1);  // f[m] = m!（34 以上は PERM_INF）
    for (int m = 1; m <= n; ++m) f[(size_t)m] = m >= 34 ? math_internal::PERM_INF : f[(size_t)m - 1] * m;
    if (n <= 33 && k >= f[(size_t)n]) math_internal::fail("perm_unrank", "k < n! で使う");
    math_internal::PermFenwick bit(n);
    for (int i = 0; i < n; ++i) bit.add(i, 1);
    std::vector<int> p((size_t)n);
    for (int i = 0; i < n; ++i) {
        const Int128 fm = f[(size_t)(n - 1 - i)];
        const int d = fm == math_internal::PERM_INF ? 0 : (int)(k / fm);
        if (fm != math_internal::PERM_INF) k %= fm;
        p[(size_t)i] = bit.kth(d), bit.add(p[(size_t)i], -1);
    }
    return p;
}

// ---- 重複のある並び（同じ値があってよい。"aabbc" の並べ替えなど）----
namespace math_internal {
inline Int128 multiset_count(const std::vector<int>& cnt, const char* func) {  // 並べ方の数 (Σ c)! / Π c!（Int128 に収まらなければ止める）
    Int128 M = 1;
    long long s = 0;
    for (int c : cnt)
        for (int j = 1; j <= c; ++j) {
            ++s;
            M = mul_div_exact(M, s, j);  // 1 個ずつ足す：M × s / (その値の個数)
            if (M == PERM_INF) fail(func, "並べ方の数が Int128 に収まらない（_mod を使う）");
        }
    return M;
}
}  // namespace math_internal
// a が「a の並べ替え（同じものは区別しない）」全体の辞書順で何番目か（0 から）  O(N log N)
template <class C>
Int128 multiset_perm_rank(const C& a) {
    const int n = (int)a.size();
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(a, k);
    std::vector<int> cnt((size_t)k, 0);
    for (int x : c) ++cnt[(size_t)x];
    Int128 M = math_internal::multiset_count(cnt, "multiset_perm_rank");  // 残りの並べ方の数
    math_internal::PermFenwick bit(k);
    for (int v = 0; v < k; ++v) bit.add(v, cnt[(size_t)v]);
    Int128 r = 0;
    for (int i = 0; i < n; ++i) {
        const long long s = n - i;
        const int u = c[(size_t)i], less = bit.prefix(u);  // 先頭に u より小さい値を置く並べ方：M × less / s
        r += math_internal::mul_div_exact(M, less, s);
        M = math_internal::mul_div_exact(M, cnt[(size_t)u], s);
        --cnt[(size_t)u], bit.add(u, -1);
    }
    return r;
}
template <class mint, class C>
mint multiset_perm_rank_mod(const C& a) {  // multiset_perm_rank を mod で（mod は素数で n より大きい）  O(N log N)
    const int n = (int)a.size();
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(a, k);
    std::vector<int> cnt((size_t)k, 0);
    for (int x : c) ++cnt[(size_t)x];
    Binomial<mint> B(n);
    mint M = B.fact(n);
    for (int x : cnt) M *= B.inv_fact(x);
    math_internal::PermFenwick bit(k);
    for (int v = 0; v < k; ++v) bit.add(v, cnt[(size_t)v]);
    mint r = 0;
    for (int i = 0; i < n; ++i) {
        const int u = c[(size_t)i];
        const mint is = B.inv(n - i);
        r += M * mint(bit.prefix(u)) * is;
        M *= mint(cnt[(size_t)u]) * is;
        --cnt[(size_t)u], bit.add(u, -1);
    }
    return r;
}
// a の並べ替え（同じものは区別しない）で辞書順 k 番目（0 から）のもの（a の順番は何でもよい）  O(N × 値の種類)
template <class C>
C multiset_perm_unrank(C a, Int128 k) {
    std::sort(a.begin(), a.end());
    const int n = (int)a.size();
    int kinds = 0;
    const std::vector<int> c = math_internal::perm_compress(a, kinds);
    std::vector<int> cnt((size_t)kinds, 0), first((size_t)kinds, 0);
    for (int i = n - 1; i >= 0; --i) ++cnt[(size_t)c[(size_t)i]], first[(size_t)c[(size_t)i]] = i;
    Int128 M = math_internal::multiset_count(cnt, "multiset_perm_unrank");
    if (k < 0 || k >= M) math_internal::fail("multiset_perm_unrank", "0 ≤ k < 並べ方の数 で使う");
    C res = a;
    for (int i = 0; i < n; ++i) {
        const long long s = n - i;
        for (int v = 0; v < kinds; ++v) {
            if (cnt[(size_t)v] == 0) continue;
            const Int128 blk = math_internal::mul_div_exact(M, cnt[(size_t)v], s);  // 先頭が v の並べ方
            if (k < blk) {
                res[(size_t)i] = a[(size_t)first[(size_t)v]], M = blk, --cnt[(size_t)v];
                break;
            }
            k -= blk;
        }
    }
    return res;
}

// ---- 組合せ（{0, …, n−1} から r 個。小さい順に並べた列の辞書順）----
// k 番目（0 から）の r 個の組（小さい順）。0 ≤ k < C(n, r)  O(n × min(r, 130))
inline std::vector<int> kth_combination(int n, int r, Int128 k) {
    if (n < 0 || r < 0 || r > n) math_internal::fail("kth_combination", "0 ≤ r ≤ n で使う");
    const Int128 tot = math_internal::binom128(n, r);
    if (k < 0 || (tot != math_internal::PERM_INF && k >= tot)) math_internal::fail("kth_combination", "0 ≤ k < C(n, r) で使う");
    std::vector<int> res;
    for (int x = 0; x < n && (int)res.size() < r; ++x) {
        const int need = r - (int)res.size();
        const Int128 c = math_internal::binom128(n - 1 - x, need - 1);  // x を選ぶ組の数（PERM_INF は「k より大きい」）
        if (c == math_internal::PERM_INF || k < c) res.push_back(x);
        else k -= c;
    }
    return res;
}
// c（{0, …, n−1} の異なる r 個。順番は何でもよい）が辞書順で何番目か（0 から。C(n, r) が Int128 に収まること）  O(n + r log r)
template <class C>
Int128 combination_rank(int n, const C& c) {
    std::vector<long long> v(c.begin(), c.end());
    std::sort(v.begin(), v.end());
    const int r = (int)v.size();
    for (int i = 0; i < r; ++i)
        if (v[(size_t)i] < 0 || v[(size_t)i] >= n || (i && v[(size_t)i] == v[(size_t)i - 1])) math_internal::fail("combination_rank", "0〜n−1 の異なる値で使う");
    if (math_internal::binom128(n, r) == math_internal::PERM_INF) math_internal::fail("combination_rank", "C(n, r) が Int128 に収まらない");
    Int128 res = 0;
    long long prev = -1;
    for (int i = 0; i < r; ++i) {  // prev+1 〜 v_i − 1 を i 番目に置く組：Σ_x C(n−1−x, r−i−1) = C(n−prev−1, r−i) − C(n−v_i, r−i)
        res += math_internal::binom128(n - prev - 1, r - i) - math_internal::binom128(n - v[(size_t)i], r - i);
        prev = v[(size_t)i];
    }
    return res;
}
#endif  // GG_PERMUTATION
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace perm_test {

#define PERM_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[permutation] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

template <class C>
static long long brute_inv(const C& a) {
    long long s = 0;
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = i + 1; j < a.size(); ++j) s += a[j] < a[i];
    return s;
}
static Int128 rnd128(std::mt19937_64& rng, Int128 hi) {  // [0, hi)
    const UInt128 x = ((UInt128)rng() << 64) | rng();
    return (Int128)(x % (UInt128)hi);
}

static void test_perm(std::mt19937_64& rng, int& cases) {
    using mint = ModInt998;
    for (int n = 0; n <= 7; ++n) {
        std::vector<int> p((size_t)n);
        std::iota(p.begin(), p.end(), 0);
        long long idx = 0;
        const long long off = rnd(rng, -5, 5), mul = rnd(rng, 1, 4);
        do {
            ++cases;
            PERM_CHECK(perm_rank(p) == idx && perm_unrank(n, idx) == p && perm_rank_mod<mint>(p) == mint(idx));
            std::vector<long long> q;  // 値をずらしても同じ番号
            for (int x : p) q.push_back(x * mul + off);
            PERM_CHECK(perm_rank(q) == idx);
            PERM_CHECK(inversion_count(p) == brute_inv(p) && perm_parity(p) == brute_inv(p) % 2 && perm_parity(q) == brute_inv(p) % 2);
            ++idx;
        } while (std::next_permutation(p.begin(), p.end()));
    }
    // n = 33：番号が Int128 いっぱいまで
    Int128 f33 = 1;
    for (int i = 1; i <= 33; ++i) f33 *= i;
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 33);
        Int128 fn = 1;
        for (int i = 1; i <= n; ++i) fn *= i;
        const Int128 k = it == 0 ? fn - 1 : rnd128(rng, fn);
        const std::vector<int> p = perm_unrank(n, k);
        PERM_CHECK(perm_rank(p) == k && perm_rank_mod<mint>(p) == mint((long long)(k % 998244353)));
    }
    // n > 33：前の方は動かない（k < 34!）
    for (int it = 0; it < 30 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 34, 300);
        const Int128 k = rnd128(rng, math_internal::PERM_INF);
        const std::vector<int> p = perm_unrank(n, k);
        for (int i = 0; i < n - 34; ++i) PERM_CHECK(p[(size_t)i] == i);
        PERM_CHECK(p[(size_t)n - 34] == n - 34 + (int)(k / f33));
        const std::vector<int> tail(p.end() - 33, p.end());
        PERM_CHECK(perm_rank(tail) == k % f33);
    }
    // 大きい n の mod 版：O(n^2) の愚直と
    for (int it = 0; it < 20 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 2000);
        std::vector<int> p((size_t)n);
        std::iota(p.begin(), p.end(), 0);
        std::shuffle(p.begin(), p.end(), rng);
        mint r = 0;
        ModInt107 r7 = 0;
        for (int i = 0; i < n; ++i) {
            int less = 0;
            for (int j = i + 1; j < n; ++j) less += p[(size_t)j] < p[(size_t)i];
            r = r * mint(n - i) + mint(less);
            r7 = r7 * ModInt107(n - i) + ModInt107(less);
        }
        PERM_CHECK(perm_rank_mod<mint>(p) == r && perm_rank_mod<ModInt107>(p) == r7);
    }
    // 転倒数：同じ値・文字列・負の値
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 60);
        std::vector<long long> a((size_t)n);
        for (auto& x : a) x = it % 2 ? rnd(rng, -3, 3) : rnd_ll(rng);
        PERM_CHECK(inversion_count(a) == brute_inv(a));
        std::string s;
        for (int i = 0; i < n; ++i) s += (char)('a' + rnd(rng, 0, 4));
        PERM_CHECK(inversion_count(s) == brute_inv(s));
    }
}

static void test_multiset(std::mt19937_64& rng, int& cases) {
    using mint = ModInt998;
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 8);
        std::vector<int> a((size_t)n);
        for (auto& x : a) x = (int)rnd(rng, 0, 3) * 7 - 5;
        std::sort(a.begin(), a.end());
        std::vector<int> orig = a;
        std::shuffle(orig.begin(), orig.end(), rng);
        long long idx = 0;
        do {
            ++cases;
            PERM_CHECK(multiset_perm_rank(a) == idx && multiset_perm_unrank(orig, idx) == a && multiset_perm_rank_mod<mint>(a) == mint(idx));
            ++idx;
        } while (std::next_permutation(a.begin(), a.end()));
    }
    {  // 文字列
        std::string s = "aabbc";
        long long idx = 0;
        do {
            ++cases;
            PERM_CHECK(multiset_perm_rank(s) == idx && multiset_perm_unrank(std::string("cbaba"), idx) == s);
            ++idx;
        } while (std::next_permutation(s.begin(), s.end()));
        PERM_CHECK(idx == 30);
    }
    // 大きい：60 文字（並べ方 ≈ 5.8×10^26）で行って戻る
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        std::string s;
        const int n = (int)rnd(rng, 1, 60), kinds = (int)rnd(rng, 1, 3);
        for (int i = 0; i < n; ++i) s += (char)('a' + rnd(rng, 0, kinds - 1));
        const Int128 r = multiset_perm_rank(s);
        PERM_CHECK(multiset_perm_unrank(s, r) == s && multiset_perm_rank_mod<mint>(s) == mint((long long)(r % 998244353)));
    }
    PERM_CHECK(multiset_perm_unrank(std::string(30, 'a') + std::string(30, 'b'), (Int128)118264581564861423LL) == std::string(30, 'b') + std::string(30, 'a'));
    // mod 版は n が大きくても（Horner の O(N^2) と）
    for (int it = 0; it < 10 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 1500);
        std::vector<int> a((size_t)n);
        for (auto& x : a) x = (int)rnd(rng, 0, 5);
        // 残りの並べ方 M を mint で持ち、M × (小さい値の個数) / 残りの長さ を足す（逆元は mint で）
        std::vector<int> cnt(6, 0);
        for (int x : a) ++cnt[(size_t)x];
        mint M = 1;
        for (int i = 1; i <= n; ++i) M *= mint(i);
        for (int c : cnt)
            for (int i = 1; i <= c; ++i) M /= mint(i);
        mint r = 0;
        for (int i = 0; i < n; ++i) {
            int less = 0;
            for (int v = 0; v < a[(size_t)i]; ++v) less += cnt[(size_t)v];
            r += M * mint(less) / mint(n - i);
            M = M * mint(cnt[(size_t)a[(size_t)i]]) / mint(n - i);
            --cnt[(size_t)a[(size_t)i]];
        }
        PERM_CHECK(multiset_perm_rank_mod<mint>(a) == r);
    }
}

static void test_comb(std::mt19937_64& rng, int& cases) {
    for (int n = 0; n <= 12; ++n)
        for (int r = 0; r <= n; ++r) {
            std::vector<std::vector<int>> all;
            std::vector<int> cur;
            std::function<void(int)> rec = [&](int x) {
                if ((int)cur.size() == r) { all.push_back(cur); return; }
                for (int y = x; y < n; ++y) cur.push_back(y), rec(y + 1), cur.pop_back();
            };
            rec(0);
            for (size_t i = 0; i < all.size(); ++i, ++cases) {
                PERM_CHECK(kth_combination(n, r, (Int128)i) == all[i] && combination_rank(n, all[i]) == (Int128)i);
                std::vector<int> sh = all[i];
                std::shuffle(sh.begin(), sh.end(), rng);
                PERM_CHECK(combination_rank(n, sh) == (Int128)i);
            }
        }
    // C(120, 60) ≈ 10^35：行って戻る
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 120), r = (int)rnd(rng, 0, n);
        const Int128 tot = math_internal::binom128(n, r);
        const Int128 k = it % 5 == 0 ? tot - 1 : rnd128(rng, tot);
        const std::vector<int> c = kth_combination(n, r, k);
        PERM_CHECK((int)c.size() == r && std::is_sorted(c.begin(), c.end()) && combination_rank(n, c) == k);
    }
    // C(n, r) が Int128 を超えても kth_combination は使える（前の方から決まる）
    for (int it = 0; it < 50 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 200, 3000), r = (int)rnd(rng, 100, n - 50);
        const Int128 k = rnd128(rng, math_internal::PERM_INF);
        const std::vector<int> c = kth_combination(n, r, k);
        PERM_CHECK((int)c.size() == r && std::adjacent_find(c.begin(), c.end(), std::greater_equal<int>()) == c.end() && c.back() < n);
        // 後ろの方を除いた部分（右端の 60 個の位置から選ぶ部分）は C(60 + …) で確かめにくいので、先頭が 0, 1, 2, … と始まることだけ
        PERM_CHECK(c[0] == 0);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_perm(rng, cases);
    test_multiset(rng, cases);
    test_comb(rng, cases);
    std::cout << "[permutation] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[permutation] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    using mint = ModInt998;
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 1000000;
    std::vector<int> p((size_t)N);
    std::iota(p.begin(), p.end(), 0);
    std::shuffle(p.begin(), p.end(), rng);
    std::vector<long long> a((size_t)N);
    for (auto& x : a) x = rnd(rng, 0, 1000000000);
    auto t0 = std::chrono::steady_clock::now();
    const long long inv = inversion_count(a);
    const long long t_inv = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const mint r = perm_rank_mod<mint>(p);
    const long long t_rank = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<int> q = perm_unrank(N, math_internal::PERM_INF - 12345);
    const long long t_un = ms(t0);
    std::vector<int> b((size_t)N);
    for (auto& x : b) x = (int)rnd(rng, 0, 9);
    t0 = std::chrono::steady_clock::now();
    const mint mr = multiset_perm_rank_mod<mint>(b);
    const long long t_ms = ms(t0);
    std::cout << "[permutation] benchmark: N = 10^6 で inversion_count " << t_inv << " ms, perm_rank_mod " << t_rank << " ms, perm_unrank " << t_un
              << " ms, multiset_perm_rank_mod " << t_ms << " ms  (check " << (inv + r.val() + q[777777] + mr.val()) % 7 << ")\n";
}

}  // namespace perm_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace perm_demo {

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
#define SHOW(expr, comment) perm_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    perm_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    using mint = ModInt998;
    banner("順列・組合せの番号（permutation） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  inversion_count / perm_parity / perm_rank / perm_rank_mod / perm_unrank\n"
              << "  multiset_perm_rank / multiset_perm_rank_mod / multiset_perm_unrank / kth_combination / combination_rank\n"
              << "  （番号はすべて 0 から。この中では using mint = ModInt998; としている）\n";

    api("【1】 inversion_count(a) / perm_parity(a)", "long long / int", "転倒数：i < j かつ a_i > a_j の組の数（値は何でも。同じ値の組は数えない。隣どうしの交換でソートする回数）", "O(N log N)");
    SHOW(inversion_count(std::vector<int>{3, 1, 2}), "(3,1)・(3,2)");
    SHOW(inversion_count(std::string("banana")), "文字列も");
    SHOW(perm_parity(std::vector<int>{3, 1, 2}), "置換の偶奇（0 = 偶置換）。転倒数 2 の偶奇。値がすべて異なる列で O(N)");
    SHOW(perm_parity(std::string("dcba")), "転倒数 6 → 偶");
    note("「2 つを入れ替える」操作で並べ替えられるかは偶奇が決め手（15 パズル・3 つの巡回は偶置換、など）");
    {
        // ABC190 F（Shift and Inversions）：先頭の a_0 を末尾に回すと、転倒数は −a_0 + (N − 1 − a_0) 変わる
        auto solve = [](std::vector<long long> a) {
            const long long n = (long long)a.size();
            long long c = inversion_count(a);
            std::vector<long long> res;
            for (long long k = 0; k < n; ++k) res.push_back(c), c += (n - 1 - a[(size_t)k]) - a[(size_t)k];
            return res;
        };
        note("ABC190 F（Shift and Inversions）: 最初だけ inversion_count、あとは先頭を回すごとに差分で");
        std::cout << "    例 1 → " << to_s(solve({0, 1, 2, 3})) << "\n    例 2 → " << to_s(solve({0, 3, 1, 5, 4, 2, 9, 6, 8, 7})) << "\n";
    }

    api("【2】 perm_rank(p) / perm_rank_mod<mint>(p) / perm_unrank(n, k)", "Int128 / mint / vector<int>",
        "値がすべて異なる列が辞書順で何番目か（値は何でも。n ≤ 33。大きいなら mod）/ 0〜n−1 の順列の k 番目", "O(N log N)");
    SHOW(perm_rank(std::vector<int>{0, 1, 2}), "いちばん最初");
    SHOW(perm_rank(std::vector<int>{3, 1, 2}), "1〜3 の順列：123, 132, 213, 231, 312 → 4 番目（0 から）");
    SHOW(perm_unrank(4, 23), "0〜3 の順列の最後（4! − 1）");
    SHOW(perm_rank_mod<mint>(perm_unrank(100, 1000000007)), "n が大きいときは mod で");
    {
        // ABC276 C（Previous Permutation）：perm_unrank(N, perm_rank(P) − 1) に 1 を足す
        auto solve = [](const std::vector<int>& P) {
            std::vector<int> q = perm_unrank((int)P.size(), perm_rank(P) - 1);
            for (int& x : q) ++x;
            return q;
        };
        note("ABC276 C（Previous Permutation）: perm_unrank(N, perm_rank(P) − 1)（N ≤ 33 のとき。std::prev_permutation でもよい）");
        std::cout << "    例 1 → " << to_s(solve({3, 1, 2})) << "\n    例 2 → " << to_s(solve({9, 8, 6, 5, 10, 3, 1, 2, 4, 7})) << "\n";
    }

    api("【3】 multiset_perm_rank(a) / multiset_perm_rank_mod<mint>(a) / multiset_perm_unrank(a, k)", "Int128 / mint / a と同じ型",
        "同じ値があってよい列の並べ替え（同じものは区別しない）で何番目か / k 番目の並び（a の順番は何でもよい）", "O(N log N) / O(N × 値の種類)");
    SHOW(multiset_perm_rank(std::string("abba")), "aabb, abab, abba → 2");
    SHOW(multiset_perm_unrank(std::string("aabb"), 3), "4 番目（0 から 3）");
    SHOW(multiset_perm_rank(std::vector<int>{2, 2, 1, 1, 1}), "最後：5! / (2! 3!) − 1 = 9");
    {
        // ABC202 D（aab aba baa）：a が A 個、b が B 個の文字列で辞書順 K 番目（1 から）
        auto solve = [](int A, int B, long long K) { return multiset_perm_unrank(std::string(A, 'a') + std::string(B, 'b'), K - 1); };
        note("ABC202 D（aab aba baa）: multiset_perm_unrank(\"a…ab…b\", K − 1)");
        std::cout << "    (2, 2, 4) → " << solve(2, 2, 4) << "\n    (30, 30, 118264581564861424) → " << solve(30, 30, 118264581564861424LL) << "\n";
    }

    api("【4】 kth_combination(n, r, k) / combination_rank(n, c)", "vector<int> / Int128",
        "{0..n−1} から r 個選ぶ組（小さい順の列）の辞書順 k 番目 / 組 c が何番目か", "O(n × min(r, 130)) / O(n + r log r)");
    SHOW(kth_combination(5, 3, 0), "最初");
    SHOW(kth_combination(5, 3, 9), "最後（C(5, 3) − 1 = 9）");
    SHOW(combination_rank(5, std::vector<int>{0, 3, 4}), "012, 013, 014, 023, 024, 034 → 5");
}

}  // namespace perm_demo

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
    perm_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    perm_test::benchmark();    // 速度計測
    perm_demo::demo();         // 使い方デモ
    return 0;
}
