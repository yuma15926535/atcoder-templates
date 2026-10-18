// =============================================================================
//  ビット（1 の数・位置、部分集合の列挙、グレイコード、0〜n の XOR・ビットの数、XOR の和）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 09 ビット
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  x は 64 ビットとして扱う（負の数は 2 の補数のビット）。ビットの位置 i は 0〜63（0 が最下位）
//  集合 {0, …, n−1} の部分集合を n ビットの数で表す（i ビット目が 1 ⇔ i を含む）。列挙はすべて小さい順
//
//  機能一覧                                                                    計算量
//   1 つの数
//    bit_count(x) / top_bit(x) / low_bit(x)     1 の数 / 最上位・最下位の 1 の位置（0 なら −1）       O(1)
//    bit_length(x) / has_bit(x, i)              2 進数の桁数（0 なら 0）/ i ビット目が 1 か
//    bit_indices(x) / for_each_bit(x, f)        1 のビットの位置（小さい順）
//    bit_reverse(x, n) / mask_to_string(mask, n)    下 n ビットの左右反転 / "0101" の形（上の桁が左）
//   部分集合の列挙（小さい順に f(s)）
//    for_each_submask(mask, f, include_zero = true, include_self = true)   s ⊆ mask（全 mask で合計 3^n）  O(2^|mask|)
//    for_each_supermask(mask, n, f)             mask ⊆ s ⊆ {0, …, n−1}
//    for_each_mask_between(must, allowed, f)    must ⊆ s ⊆ allowed（必ず 1 のビットと、使ってよいビット）
//    for_each_mask_between(must, allowed, k, f) そのうち 1 の数がちょうど k のもの                     O(個数 × k)
//    for_each_subset_k(n, k, f)                 {0, …, n−1} の k 個の部分集合（ゴスパーの方法）       O(C(n, k))
//   グレイコード
//    gray_code(i) / gray_inverse(g)             i 番目（隣どうしは 1 ビットだけ違う）/ 何番目か
//   0〜n の XOR・ビットの数（n ≥ 0）
//    xor_up_to(n) / xor_range(l, r) / xor_range_closed(L, R)    0 ⊕ … ⊕ n / l ⊕ … ⊕ (r − 1) / L ⊕ … ⊕ R   O(1)
//    count_bit_up_to(n, k)                      0〜n のうち k ビット目が 1 の数                        O(1)
//    popcount_sum_up_to(n)                      Σ_{x=0}^{n} bit_count(x)（Int128）                     O(63)
//   XOR の和（a_i ≥ 0。Int128。a は vector<long long>・vector<int> など）
//    xor_pair_sum(a)                            Σ_{i<j} (a_i ⊕ a_j)                                    O(63 n)
//    subarray_xor_sum(a)                        すべての区間（長さ 1 以上）の XOR の和                 O(63 n)
//  （ビット DP・高速ゼータ変換は 07_subset_transform、XOR 基底は 05_xor_basis）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    bit_test::random_test()  … 全部調べる愚直解との突き合わせ
//    bit_test::benchmark()    … 速度計測
//    bit_demo::demo()         … 全機能の使い方デモ（日本語の説明付き。bit 全探索の形も）
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
#ifndef GG_BIT
#define GG_BIT
// ---- ビット（数学第二世代）  ビットの数・位置、部分集合の列挙、グレイコード、0〜n の XOR・ビットの数、XOR の和 ----
//   x は 64 ビットとして扱う（負の数は 2 の補数のビット）。位置 i は 0〜63（0 が最下位）
inline int bit_count(long long x) { return __builtin_popcountll((unsigned long long)x); }                // 1 の数
inline int top_bit(long long x) { return x == 0 ? -1 : 63 - __builtin_clzll((unsigned long long)x); }   // 最上位の 1 の位置（0 なら −1）
inline int low_bit(long long x) { return x == 0 ? -1 : __builtin_ctzll((unsigned long long)x); }       // 最下位の 1 の位置（0 なら −1。値がほしいなら x & -x）
inline int bit_length(long long x) { return top_bit(x) + 1; }                                          // 2 進数の桁数（0 なら 0）
inline bool has_bit(long long x, int i) {                                                              // i ビット目が 1 か
    if (i < 0 || i > 63) math_internal::fail("has_bit", "i は 0〜63");
    return ((unsigned long long)x >> i) & 1;
}
inline std::vector<int> bit_indices(long long x) {  // 1 のビットの位置（小さい順）
    std::vector<int> r;
    for (unsigned long long u = (unsigned long long)x; u; u &= u - 1) r.push_back(__builtin_ctzll(u));
    return r;
}
template <class F>
void for_each_bit(long long x, F f) {  // 1 のビットの位置 i を小さい順に f(i)
    for (unsigned long long u = (unsigned long long)x; u; u &= u - 1) f(__builtin_ctzll(u));
}
inline long long bit_reverse(long long x, int n) {  // 下 n ビットを左右反転（n は 0〜63。上のビットは 0 になる）
    if (n < 0 || n > 63) math_internal::fail("bit_reverse", "n は 0〜63");
    unsigned long long r = 0;
    for (int i = 0; i < n; ++i) r |= (((unsigned long long)x >> i) & 1) << (n - 1 - i);
    return (long long)r;
}
inline std::string mask_to_string(long long mask, int n) {  // 下 n ビットを "0101" の形に（上の桁が左。n は 0〜64）
    if (n < 0 || n > 64) math_internal::fail("mask_to_string", "n は 0〜64");
    std::string s((size_t)n, '0');
    for (int i = 0; i < n; ++i)
        if (((unsigned long long)mask >> i) & 1) s[(size_t)(n - 1 - i)] = '1';
    return s;
}

// ---- 部分集合の列挙（小さい順） ----
namespace math_internal {
inline void check_width(int n, const char* func) {
    if (n < 0 || n > 63) fail(func, "n は 0〜63");
}
}  // namespace math_internal
template <class F>
void for_each_submask(long long mask, F f, bool include_zero = true, bool include_self = true) {  // s ⊆ mask を小さい順に f(s)（3^n 通りの部分集合 DP にも）
    if (mask < 0) math_internal::fail("for_each_submask", "mask ≥ 0 で使う");
    const unsigned long long m = (unsigned long long)mask;
    unsigned long long s = 0;
    do {
        if ((s != 0 || include_zero) && (s != m || include_self)) f((long long)s);
        s = (s - m) & m;  // 次に大きい部分集合
    } while (s != 0);
}
template <class F>
void for_each_supermask(long long mask, int n, F f) {  // mask ⊆ s ⊆ {0, …, n−1} を小さい順に f(s)
    math_internal::check_width(n, "for_each_supermask");
    if (mask < 0 || (unsigned long long)mask >> n) math_internal::fail("for_each_supermask", "mask は 0 以上 2^n 未満");
    const unsigned long long full = (1ULL << n) - 1;
    for (unsigned long long s = (unsigned long long)mask;; s = (s + 1) | (unsigned long long)mask) {
        f((long long)s);
        if (s == full) break;
    }
}
template <class F>
void for_each_mask_between(long long must, long long allowed, F f) {  // must ⊆ s ⊆ allowed を小さい順に f(s)（必ず 1 のビットと使ってよいビット。must ⊄ allowed なら何もしない）
    if (must < 0 || allowed < 0) math_internal::fail("for_each_mask_between", "must, allowed ≥ 0 で使う");
    if (must & ~allowed) return;
    for_each_submask(allowed & ~must, [&](long long s) { f(must | s); });
}
template <class F>
void for_each_mask_between(long long must, long long allowed, int k, F f) {  // さらに 1 の数がちょうど k のものだけ（小さい順）
    if (must < 0 || allowed < 0) math_internal::fail("for_each_mask_between", "must, allowed ≥ 0 で使う");
    if (must & ~allowed) return;
    const std::vector<int> pos = bit_indices(allowed & ~must);  // 自由なビット
    const int c = (int)pos.size(), t = k - bit_count(must);
    if (t < 0 || t > c) return;
    if (t == 0) { f(must); return; }
    // 自由なビットの番号の上でゴスパーの方法（allowed ≥ 0 なので c ≤ 63）
    const unsigned long long end = 1ULL << c;
    for (unsigned long long x = (1ULL << t) - 1; x < end;) {
        unsigned long long s = (unsigned long long)must;
        for (unsigned long long u = x; u; u &= u - 1) s |= 1ULL << pos[(size_t)__builtin_ctzll(u)];
        f((long long)s);
        const unsigned long long lo = x & (0 - x), y = x + lo;
        x = (((x ^ y) >> 2) / lo) | y;
    }
}
template <class F>
void for_each_subset_k(int n, int k, F f) {  // {0, …, n−1} の k 個の部分集合（1 の数が k の n ビットの数）を小さい順に f(s)
    math_internal::check_width(n, "for_each_subset_k");
    if (k < 0 || k > n) return;
    if (k == 0) { f(0LL); return; }
    const unsigned long long end = 1ULL << n;
    for (unsigned long long x = (1ULL << k) - 1; x < end;) {
        f((long long)x);
        const unsigned long long lo = x & (0 - x), y = x + lo;
        x = (((x ^ y) >> 2) / lo) | y;
    }
}

// ---- グレイコード（隣どうしが 1 ビットだけ違う並び） ----
inline long long gray_code(long long i) { return i ^ (long long)((unsigned long long)i >> 1); }  // i 番目のグレイコード
inline long long gray_inverse(long long g) {                                                       // gray_code の逆（何番目か）
    unsigned long long u = (unsigned long long)g;
    for (int s = 1; s < 64; s <<= 1) u ^= u >> s;
    return (long long)u;
}

// ---- 0〜n の XOR・ビットの数（n ≥ 0） ----
inline long long xor_up_to(long long n) {  // 0 ⊕ 1 ⊕ … ⊕ n（n < 0 なら 0）
    if (n < 0) return 0;
    switch (n & 3) {
        case 0: return n;
        case 1: return 1;
        case 2: return n + 1;
        default: return 0;
    }
}
inline long long xor_range(long long l, long long r) {  // l ⊕ … ⊕ (r − 1)（0 ≤ l。空なら 0）
    if (l < 0) math_internal::fail("xor_range", "l ≥ 0 で使う");
    return r <= l ? 0 : xor_up_to(r - 1) ^ xor_up_to(l - 1);
}
inline long long xor_range_closed(long long L, long long R) {  // L ⊕ … ⊕ R（0 ≤ L。ABC121 D）
    if (L < 0) math_internal::fail("xor_range_closed", "L ≥ 0 で使う");
    return R < L ? 0 : xor_up_to(R) ^ xor_up_to(L - 1);
}
inline long long count_bit_up_to(long long n, int k) {  // 0〜n のうち k ビット目が 1 の数（k は 0〜63）
    if (k < 0 || k > 63) math_internal::fail("count_bit_up_to", "k は 0〜63");
    if (n < 0 || k == 63) return 0;
    const unsigned long long m = (unsigned long long)n + 1, p = 1ULL << k;  // 0〜n は m 個。周期 2p のうち後半 p 個が 1
    const unsigned long long full = (m >> (k + 1)) << k, rest = m & (2 * p - 1);
    return (long long)(full + (rest > p ? rest - p : 0));
}
inline Int128 popcount_sum_up_to(long long n) {  // Σ_{x=0}^{n} bit_count(x)（n < 0 なら 0）
    if (n < 0) return 0;
    const unsigned long long m = (unsigned long long)n + 1;  // 0〜m−1 を、m の上の桁から 1 のビットごとに区切って数える
    UInt128 s = 0;
    unsigned long long ones = 0;  // ここより上で 1 のビットの数
    for (unsigned long long u = m; u; ++ones) {
        const int k = 63 - __builtin_clzll(u);
        u ^= 1ULL << k;
        s += ((UInt128)(unsigned long long)k << k >> 1) + ((UInt128)ones << k);  // 下 k ビットが全部の形（k·2^(k−1) 個）＋ 上の 1 の分
    }
    return (Int128)s;
}

// ---- XOR の和（a_i ≥ 0） ----
// （a は vector<long long>・vector<int> など。{1, 2, 3} とも書ける）
template <class V>
Int128 xor_pair_sum(const V& a) {  // Σ_{i<j} (a_i ⊕ a_j)（ABC147 D）  O(63 n)
    long long cnt[63] = {};
    for (const auto& x : a) {
        const long long v = math_internal::to_ll_checked(x, "xor_pair_sum");
        if (v < 0) math_internal::fail("xor_pair_sum", "a_i ≥ 0 で使う");
        for (unsigned long long u = (unsigned long long)v; u; u &= u - 1) cnt[__builtin_ctzll(u)]++;
    }
    const long long n = (long long)a.size();
    Int128 s = 0;
    for (int k = 0; k < 63; ++k) s += (Int128)cnt[k] * (n - cnt[k]) << k;
    return s;
}
inline Int128 xor_pair_sum(std::initializer_list<long long> a) { return xor_pair_sum<std::initializer_list<long long>>(a); }
template <class V>
Int128 subarray_xor_sum(const V& a) {  // すべての区間 [l, r)（長さ 1 以上）の XOR の和（ABC365 E）  O(63 n)
    for (const auto& x : a)
        if (math_internal::to_ll_checked(x, "subarray_xor_sum") < 0) math_internal::fail("subarray_xor_sum", "a_i ≥ 0 で使う");
    Int128 s = 0;
    for (int k = 0; k < 63; ++k) {
        long long c[2] = {1, 0};  // 累積 XOR の k ビット目が 0 / 1 の個数（はじめの 0 を含む）
        int cur = 0;
        for (const auto& x : a) cur ^= (int)(((long long)x >> k) & 1), c[cur]++;
        s += (Int128)c[0] * c[1] << k;
    }
    return s;
}
inline Int128 subarray_xor_sum(std::initializer_list<long long> a) { return subarray_xor_sum<std::initializer_list<long long>>(a); }
#endif  // GG_BIT
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace bit_test {

#define BIT_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[bit] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static void test_single(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300000 * TEST_SCALE; ++it, ++cases) {
        const long long x = rnd_ll(rng);
        const unsigned long long u = (unsigned long long)x;
        int cnt = 0, top = -1, low = -1;
        std::vector<int> idx;
        for (int i = 0; i < 64; ++i)
            if ((u >> i) & 1) {
                cnt++, top = i, idx.push_back(i);
                if (low < 0) low = i;
            }
        BIT_CHECK(bit_count(x) == cnt && top_bit(x) == top && low_bit(x) == low && bit_length(x) == top + 1 && bit_indices(x) == idx);
        std::vector<int> got;
        for_each_bit(x, [&](int i) { got.push_back(i); });
        BIT_CHECK(got == idx);
        const int i = (int)rnd(rng, 0, 63);
        BIT_CHECK(has_bit(x, i) == (((u >> i) & 1) == 1));
        const int n = (int)rnd(rng, 0, 63);
        unsigned long long rev = 0;
        std::string s;
        for (int j = n - 1; j >= 0; --j) s += ((u >> j) & 1) ? '1' : '0';
        for (int j = 0; j < n; ++j)
            if ((u >> j) & 1) rev |= 1ULL << (n - 1 - j);
        BIT_CHECK(bit_reverse(x, n) == (long long)rev && mask_to_string(x, n) == s);
        // グレイコード
        const long long g = gray_code(x);
        BIT_CHECK(gray_inverse(g) == x && bit_count(gray_code(x) ^ gray_code((long long)(u + 1))) == 1);
    }
    BIT_CHECK(mask_to_string(-1, 64) == std::string(64, '1') && mask_to_string(5, 0) == "" && bit_count(-1) == 64 && top_bit(-1) == 63);
}

static void test_enumerate(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 11);
        const long long full = (1LL << n) - 1;
        const long long mask = (long long)(rng() & (unsigned long long)full), mask2 = (long long)(rng() & (unsigned long long)full);
        // 部分集合（4 通りのフラグ）
        for (int fl = 0; fl < 4; ++fl) {
            const bool z = fl & 1, self = fl & 2;
            std::vector<long long> expect, got;
            for (long long s = 0; s <= full; ++s)
                if ((s & ~mask) == 0 && (s != 0 || z) && (s != mask || self)) expect.push_back(s);
            for_each_submask(mask, [&](long long s) { got.push_back(s); }, z, self);
            BIT_CHECK(got == expect);
        }
        {
            std::vector<long long> expect, got;
            for (long long s = 0; s <= full; ++s)
                if ((s & mask) == mask) expect.push_back(s);
            for_each_supermask(mask, n, [&](long long s) { got.push_back(s); });
            BIT_CHECK(got == expect);
        }
        // must ⊆ s ⊆ allowed（must = mask & mask2 なら必ず含まれる。そうでないときは空）
        const long long must = it % 3 == 0 ? (mask & mask2) : (long long)(rng() & (unsigned long long)full), allowed = mask;
        const int k = (int)rnd(rng, -1, n + 1);
        std::vector<long long> e1, e2, g1, g2;
        for (long long s = 0; s <= full; ++s)
            if ((s & must) == must && (s & ~allowed) == 0) {
                e1.push_back(s);
                if (bit_count(s) == k) e2.push_back(s);
            }
        for_each_mask_between(must, allowed, [&](long long s) { g1.push_back(s); });
        for_each_mask_between(must, allowed, k, [&](long long s) { g2.push_back(s); });
        BIT_CHECK(g1 == e1 && g2 == e2);
        // k 個の部分集合
        std::vector<long long> e3, g3;
        for (long long s = 0; s <= full; ++s)
            if (bit_count(s) == k) e3.push_back(s);
        for_each_subset_k(n, k, [&](long long s) { g3.push_back(s); });
        BIT_CHECK(g3 == e3);
    }
    // 大きいビット：63 ビット目まで
    long long c = 0, last = 0;
    for_each_subset_k(63, 62, [&](long long s) { c++, last = s; });
    BIT_CHECK(c == 63 && last == LMAX - 1);  // 1〜62 ビット目
    c = 0;
    for_each_subset_k(63, 1, [&](long long s) { c++, last = s; });
    BIT_CHECK(c == 63 && last == (1LL << 62));
    c = 0;
    for_each_mask_between(1LL << 62, LMAX, 2, [&](long long s) { c++, last = s; });
    BIT_CHECK(c == 62 && last == (3LL << 61));
    c = 0;
    for_each_supermask(LMAX - 1, 63, [&](long long s) { c++, last = s; });
    BIT_CHECK(c == 2 && last == LMAX);
    c = 0;
    for_each_submask(1LL << 62 | 1, [&](long long s) { c++, last = s; });
    BIT_CHECK(c == 4 && last == (1LL << 62 | 1));
}

static void test_counting(std::mt19937_64& rng, int& cases) {
    long long x = 0;
    Int128 pc = 0;
    std::vector<long long> cnt(63, 0);
    for (long long n = 0; n <= 300000; ++n, ++cases) {
        x ^= n, pc += bit_count(n);
        for (int k = 0; k < 20; ++k) cnt[k] += (n >> k) & 1;
        BIT_CHECK(xor_up_to(n) == x && popcount_sum_up_to(n) == pc);
        const int k = (int)(n % 23);
        BIT_CHECK(count_bit_up_to(n, k) == cnt[k]);
    }
    BIT_CHECK(xor_up_to(-1) == 0 && popcount_sum_up_to(-5) == 0 && count_bit_up_to(-1, 3) == 0 && count_bit_up_to(LMAX, 63) == 0);
    for (int it = 0; it < 300000 * TEST_SCALE; ++it, ++cases) {  // 大きい n は 1 つ前との差
        const long long n = rnd(rng, 1, LMAX);
        const int k = (int)rnd(rng, 0, 63);
        BIT_CHECK((xor_up_to(n) ^ xor_up_to(n - 1)) == n);
        BIT_CHECK(count_bit_up_to(n, k) - count_bit_up_to(n - 1, k) == (long long)((n >> k) & 1));
        BIT_CHECK(popcount_sum_up_to(n) - popcount_sum_up_to(n - 1) == bit_count(n));
        // 区間の XOR
        const long long l = rnd(rng, 0, 300), r = rnd(rng, 0, 300);
        long long e = 0;
        for (long long v = l; v < r; ++v) e ^= v;
        BIT_CHECK(xor_range(l, r) == e && xor_range_closed(l, r - 1) == e);
        const long long L = rnd(rng, 0, LMAX), R = L + std::min<long long>(LMAX - L, rnd(rng, 0, 50));
        long long e2 = 0;
        for (long long v = L;; ++v) {
            e2 ^= v;
            if (v == R) break;
        }
        BIT_CHECK(xor_range_closed(L, R) == e2);
    }
    for (int k = 0; k < 63; ++k) BIT_CHECK(count_bit_up_to(LMAX, k) == (1LL << 62));
    BIT_CHECK(popcount_sum_up_to(LMAX) == (Int128)63 << 62);
}

static void test_xor_sums(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 30);
        std::vector<long long> a(n);
        for (auto& v : a) v = it % 2 ? rnd(rng, 0, 15) : rnd(rng, 0, LMAX);
        Int128 e1 = 0, e2 = 0;
        for (int i = 0; i < n; ++i) {
            long long cur = 0;
            for (int j = i; j < n; ++j) {
                cur ^= a[j], e2 += cur;
                if (j > i) e1 += a[i] ^ a[j];
            }
        }
        BIT_CHECK(xor_pair_sum(a) == e1 && subarray_xor_sum(a) == e2);
        bool small = true;
        for (long long v : a) small &= v <= INT_MAX;
        if (small) {  // vector<int> でも
            const std::vector<int> ai(a.begin(), a.end());
            BIT_CHECK(xor_pair_sum(ai) == e1 && subarray_xor_sum(ai) == e2);
        }
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_single(rng, cases);
    test_enumerate(rng, cases);
    test_counting(rng, cases);
    test_xor_sums(rng, cases);
    std::cout << "[bit] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[bit] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    const int N = 16;
    for (long long m = 0; m < (1LL << N); ++m) for_each_submask(m, [&](long long s) { s1 += s; });  // 3^16 ≈ 4.3×10^7
    const long long t_sub = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for_each_subset_k(40, 6, [&](long long s) { s2 ^= s; });  // C(40, 6) ≈ 3.8×10^6
    for_each_mask_between(1, (1LL << 40) - 1, 6, [&](long long s) { s2 += s & 7; });
    const long long t_k = ms(t0);
    const int n = 1000000;
    std::vector<long long> a(n);
    for (auto& v : a) v = rnd(rng, 0, (1LL << 60) - 1);
    t0 = std::chrono::steady_clock::now();
    const Int128 s3 = xor_pair_sum(a) + subarray_xor_sum(a);
    const long long t_sum = ms(t0);
    t0 = std::chrono::steady_clock::now();
    Int128 s4 = 0;
    for (int i = 0; i < n; ++i) s4 += popcount_sum_up_to(a[i]) + count_bit_up_to(a[i], i % 60) + xor_up_to(a[i]);
    const long long t_cnt = ms(t0);
    std::cout << "[bit] benchmark: 全 mask の for_each_submask（n = 16、3^16 回）" << t_sub << " ms, for_each_subset_k(40, 6) + for_each_mask_between(k = 6) "
              << t_k << " ms\n"
              << "[bit] benchmark: xor_pair_sum + subarray_xor_sum（n = 10^6）" << t_sum << " ms, popcount_sum_up_to + count_bit_up_to + xor_up_to × 10^6 " << t_cnt
              << " ms  (check " << (long long)((s3 + s4 + s1 + s2) % 7) << ")\n";
}

}  // namespace bit_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace bit_demo {

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
#define SHOW(expr, comment) bit_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    bit_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("ビット（bit） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  bit_count / top_bit / low_bit / bit_length / has_bit / bit_indices / for_each_bit   1 つの数のビット\n"
              << "  bit_reverse / mask_to_string                                  反転 / \"0101\" の形\n"
              << "  for_each_submask / for_each_supermask / for_each_mask_between / for_each_subset_k   部分集合の列挙（小さい順）\n"
              << "  gray_code / gray_inverse                                      グレイコード\n"
              << "  xor_up_to / xor_range / xor_range_closed / count_bit_up_to / popcount_sum_up_to   0〜n の XOR・ビットの数\n"
              << "  xor_pair_sum / subarray_xor_sum                               XOR の和\n";

    api("【1】 bit_count(x) / top_bit(x) / low_bit(x) / bit_length(x) / has_bit(x, i) / bit_indices(x)", "int / bool / vector<int>",
        "1 の数 / 最上位・最下位の 1 の位置（0 なら −1）/ 桁数 / i ビット目 / 1 の位置の一覧", "O(1)");
    SHOW(bit_count(13), "13 = 1101");
    SHOW(top_bit(13), "最上位は 3 ビット目");
    SHOW(low_bit(12), "12 = 1100 の最下位の 1 は 2 ビット目（値なら 12 & -12 = 4）");
    SHOW(bit_length(0), "0 は 0 桁");
    SHOW(has_bit(13, 1), "1 ビット目は 0");
    SHOW(bit_indices(13), "1 のビットの位置");
    SHOW(mask_to_string(13, 6), "上の桁が左");
    SHOW(bit_reverse(1, 4), "0001 → 1000");

    api("【2】 for_each_submask(mask, f, include_zero, include_self) / for_each_supermask(mask, n, f)", "なし（f(s) を呼ぶ）",
        "mask の部分集合 / mask を含む n ビットの集合を小さい順に。0 や mask 自身を除くこともできる", "O(個数)");
    RUN(std::vector<long long> v; for_each_submask(11, [&](long long s) { v.push_back(s); }));
    SHOW(v, "ABC269 C（Submask）: 11 = 1011 の部分集合を小さい順に（答えそのもの）");
    RUN(v.clear(); for_each_submask(6, [&](long long s) { v.push_back(s); }, false, false));
    SHOW(v, "0 と自分自身を除く（空でない真部分集合）");
    RUN(v.clear(); for_each_supermask(5, 4, [&](long long s) { v.push_back(s); }));
    SHOW(v, "0101 を含む 4 ビットの集合");
    note("部分集合 DP：for (m = 0; m < 2^n; ++m) for_each_submask(m, …) は全部で 3^n 回（n = 16 で約 4×10^7）");

    api("【3】 for_each_mask_between(must, allowed, f) / for_each_mask_between(must, allowed, k, f) / for_each_subset_k(n, k, f)", "なし（f(s) を呼ぶ）",
        "必ず 1 のビット must と使ってよいビット allowed を決めた列挙（1 の数 k も指定できる）/ n 個から k 個選ぶ", "O(個数) / O(個数 × k)");
    RUN(v.clear(); for_each_mask_between(0b0001, 0b1011, [&](long long s) { v.push_back(s); }));
    SHOW(v, "0 ビット目は必ず 1、2 ビット目は必ず 0");
    RUN(v.clear(); for_each_mask_between(0b0001, 0b1111, 2, [&](long long s) { v.push_back(s); }));
    SHOW(v, "さらに 1 の数がちょうど 2");
    RUN(v.clear(); for_each_subset_k(4, 2, [&](long long s) { v.push_back(s); }));
    SHOW(v, "4 個から 2 個選ぶ（C(4, 2) = 6 通り）");

    head("【4】 bit 全探索（使い方関数の形）  ABC197 C（ORXOR）：区切り方を 2^(N−1) 通り全部試す");
    {
        auto solve = [](const std::vector<long long>& a) {
            const int n = (int)a.size();
            long long best = LLONG_MAX;
            for (long long cut = 0; cut < (1LL << (n - 1)); ++cut) {  // i ビット目が 1 なら a[i] と a[i+1] のあいだで区切る
                long long x = 0, o = 0;
                for (int i = 0; i < n; ++i) {
                    o |= a[i];
                    if (i == n - 1 || has_bit(cut, i)) x ^= o, o = 0;
                }
                best = std::min(best, x);
            }
            return best;
        };
        code("for (long long cut = 0; cut < (1LL << (n - 1)); ++cut)  // has_bit(cut, i) なら i と i+1 のあいだで区切る");
        std::cout << "    A = {1, 5, 7} → " << solve({1, 5, 7}) << "、A = {10, 10, 10} → " << solve({10, 10, 10}) << "、A = {1, 3, 3, 1} → " << solve({1, 3, 3, 1}) << "\n";
    }

    api("【5】 gray_code(i) / gray_inverse(g)", "long long", "i 番目のグレイコード（隣どうしは 1 ビットだけ違う）/ その逆", "O(1)");
    RUN(v.clear(); for (int i = 0; i < 8; ++i) v.push_back(gray_code(i)));
    SHOW(v, "0, 1, 3, 2, 6, 7, 5, 4");
    SHOW(gray_inverse(5), "5 は 6 番目");
    note("集合を 1 つずつ足し引きしながら全部の部分集合を回る（差分更新の bit 全探索）ときに使う");

    api("【6】 xor_up_to(n) / xor_range(l, r) / xor_range_closed(L, R) / count_bit_up_to(n, k) / popcount_sum_up_to(n)", "long long / Int128",
        "0 ⊕ … ⊕ n / 区間の XOR / 0〜n で k ビット目が 1 の数 / 0〜n の 1 の数の合計", "O(1) / O(63)");
    SHOW(xor_range_closed(2, 4), "ABC121 D: 2 ⊕ 3 ⊕ 4 = 5");
    SHOW(xor_range_closed(123, 456), "ABC121 D: 435");
    SHOW(xor_range(2, 5), "[2, 5) も同じ");
    SHOW(count_bit_up_to(10, 1), "0〜10 で 1 ビット目が 1：2, 3, 6, 7, 10");
    SHOW(popcount_sum_up_to(7), "0〜7 の 1 の数の合計 = 12");
    {
        // ABC356 D：Σ_{k=0}^{N} popcount(k & M) mod 998244353 = Σ_{M の i ビット目が 1} count_bit_up_to(N, i)
        auto solve = [](long long N, long long M) {
            long long s = 0;
            for_each_bit(M, [&](int i) { s = (s + count_bit_up_to(N, i)) % 998244353; });
            return s;
        };
        note("ABC356 D（Masked Popcount）: M の 1 のビット i について count_bit_up_to(N, i) を足す");
        std::cout << "    N = 4, M = 3 → " << solve(4, 3) << "、N = M = 2^60 − 1 → " << solve((1LL << 60) - 1, (1LL << 60) - 1) << "\n";
    }

    api("【7】 xor_pair_sum(a) / subarray_xor_sum(a)", "Int128",
        "Σ_{i<j} (a_i ⊕ a_j) / すべての区間の XOR の和（a_i ≥ 0。ビットごとに数える）", "O(63 n)");
    SHOW(xor_pair_sum({1, 2, 3}), "ABC147 D: 3 + 2 + 1 = 6（答えは % 1000000007）");
    SHOW((long long)(xor_pair_sum({3, 14, 159, 2653, 58979, 323846, 2643383, 27950288, 419716939, 9375105820LL}) % 1000000007), "ABC147 D の例 3");
    SHOW(subarray_xor_sum({1, 3, 2}), "長さ 1 の区間も含む：1 + 3 + 2 + 2 + 1 + 0 = 9");
    {
        const std::vector<long long> a = {2, 5, 6, 5, 2, 1, 7};
        Int128 s = subarray_xor_sum(a);
        for (long long x : a) s -= x;
        note("ABC365 E（Xor Sigma Problem）: 長さ 2 以上なので subarray_xor_sum(a) − Σ a_i");
        std::cout << "    A = {1, 3, 2} → " << subarray_xor_sum({1, 3, 2}) - 6 << "、A = {2, 5, 6, 5, 2, 1, 7} → " << s << "\n";
    }
    head("【8】 ビットの式の例  ABC238 D（AND and SUM）：x & y = a、x + y = s となる x, y ≥ 0 があるか");
    {
        auto ok = [](long long a, long long s) { return s - 2 * a >= 0 && ((s - 2 * a) & a) == 0; };  // x + y = 2(x & y) + (x ⊕ y)
        code("s − 2a ≥ 0 かつ ((s − 2a) & a) == 0   // x + y = 2(x & y) + (x ⊕ y)");
        std::cout << "    (a, s) = (1, 8) → " << (ok(1, 8) ? "Yes" : "No") << "、(4, 2) → " << (ok(4, 2) ? "Yes" : "No") << "\n";
    }
}

}  // namespace bit_demo

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
    bit_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    bit_test::benchmark();    // 速度計測
    bit_demo::demo();         // 使い方デモ
    return 0;
}
