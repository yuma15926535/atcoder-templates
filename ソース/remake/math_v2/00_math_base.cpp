// =============================================================================
//  数学基盤（割り算・倍数・累乗根・累乗と溢れ・gcd / lcm・64 ビットの mod・数列の和・分数の比較）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 00 基盤（ほかの数学テンプレートはすべてこれを同梱している）
//  依存       : 00_int_common（同梱済み。Int128・isqrt・parse_fixed などは幾何第二世代と共通）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  共通の約束（第二世代の数学テンプレートすべてで同じ）
//    - 整数は long long が基本。途中で溢れるところは Int128 で計算する。答えが long long に収まらないことが
//      あるものは Int128 で返すか、上限で止める版（名前に _cap）を用意する
//    - 黙って溢れない：厳密な値を返す関数（pow_int・arith_sum など）は、溢れるときメッセージを出して止める
//    - 使える範囲の外の引数（倍数の d ≤ 0、mod の m ≤ 0 など）も、メッセージを出して止める（黙って変な値を返さない）
//    - 負の数の割り算・余りは floor_div / ceil_div / mod_floor を使う（C++ の / と % は 0 に向かって切り捨てる）
//    - 範囲の引数は半開区間 [l, r)。「L 以上 R 以下」で使うことが多いものは _closed 版もある
//    - 解がないことがある関数は ok を持つ構造体で返す（Maybe<T>{ok, value} など。-1 などの番兵で返さない）
//    - modint は 00_modint（using mint = ModInt998; を使う側で書く）
//
//  機能一覧                                                                    計算量
//    Maybe<T>{ok, value}                ないことがある値（ok が false なら value は使わない）
//   割り算・倍数（b ≠ 0。LLONG_MIN / -1 は答えが収まらない）
//    floor_div(a, b) / ceil_div(a, b)   ⌊a/b⌋ / ⌈a/b⌉（負の数でも）                        O(1)
//    mod_floor(a, m)                    a を m で割った余り（0 以上 m 未満。-7 mod 3 = 2）
//    round_div(a, b)                    a/b の四捨五入（.5 は大きい方へ）
//    ceil_multiple(x, d) / floor_multiple(x, d)          x 以上で最小 / x 以下で最大の d の倍数
//    ceil_congruent(x, a, m) / floor_congruent(x, a, m)  x 以上で最小 / x 以下で最大の y ≡ a (mod m)
//    count_multiples(l, r, d) / count_multiples_closed(L, R, d)        [l, r) / [L, R] にある d の倍数の数
//    count_congruent(l, r, a, m) / count_congruent_closed(L, R, a, m)  [l, r) / [L, R] にある x ≡ a (mod m) の数
//   累乗根（誤差なし）
//    isqrt(n) / isqrt_ceil(n)           ⌊√n⌋ / ⌈√n⌉（00_int_common）                       O(1)
//    kth_root_floor(n, k) / kth_root_ceil(n, k)   ⌊n^(1/k)⌋ / ⌈n^(1/k)⌉（n ≥ 0, k ≥ 1）       O(k)
//    is_square(n) / is_kth_power(n, k)  平方数か / k 乗数か
//   累乗・溢れ
//    pow_int(a, e)                      a^e を厳密に（溢れたら止める）                       O(log e)
//    pow_cap(a, e, cap)                 min(a^e, cap)（a ≥ 0。溢れない。「a^e が cap を超えるか」に）  O(log e)
//    mul_cap(a, b, cap) / add_cap(a, b, cap)      min(a·b, cap) / min(a+b, cap)（0 以上の値）      O(1)
//    mul_overflows(a, b) / add_overflows(a, b)    a·b / a+b が long long に収まらないか           O(1)
//    ilog_floor(n, b) / ilog_ceil(n, b) b^k ≤ n となる最大の k / b^k ≥ n となる最小の k          O(log n)
//    pow10_ll(k)                        10^k（0 ≤ k ≤ 18）
//   gcd・lcm（gcd / lcm そのものは std::gcd / std::lcm）
//    ext_gcd(a, b)                      ExtGcd{g, x, y}：a x + b y = g = gcd(|a|, |b|)          O(log)
//    lcm_cap(a, b, cap)                 min(lcm(a, b), cap)（溢れない）
//    try_lcm(a, b)                      lcm（Maybe：long long に収まらなければ ok = false）
//    gcd_of(v) / lcm_of_cap(v, cap)     全部の gcd / min(全部の lcm, cap)（vector<int> でも、gcd_of({a, b, c}) でも）  O(N log)
//   64 ビットの mod（1 ≤ m < 2^63。a は負でもよい。答えは [0, m)）
//    mod_mul(a, b, m)                   a·b mod m（128 ビットで計算）                           O(1)
//    mod_pow(a, e, m)                   a^e mod m（e < 0 なら逆元の累乗）                       O(log e)
//    mod_inv(a, m) / try_mod_inv(a, m)  逆元（m は合成数でもよい。ないと止まる / Maybe で返す）     O(log m)
//    geom_sum_mod(r, n, m)              1 + r + … + r^(n-1) mod m（逆元を使わないので m は何でも）   O(log n)
//   数列の和
//    arith_sum(a, d, n)                 a + (a+d) + … + (a+(n-1)d)（Int128。収まらなければ止める）   O(1)
//    arith_sum_mod(a, d, n, m)          同じ和を mod m で（m は何でも）
//    sum_range(l, r) / sum_range_closed(L, R)     l + … + (r-1) / L + … + R（Int128）
//    sum_range_mod(l, r, m) / sum_range_closed_mod(L, R, m)   同じ和を mod m で（[0, m)）
//   分数
//    frac_compare(a, b, c, d)           a/b と c/d の大小（-1 / 0 / 1。b, d ≠ 0。誤差なし）        O(1)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mathb_test::random_test()  … 128 ビット・ループの愚直計算と全機能を突き合わせるランダムテスト
//    mathb_test::benchmark()    … 速度計測
//    mathb_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace mathb_test {

#define MATHB_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[math_base] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static const long long LMAX = std::numeric_limits<long long>::max(), LMIN = std::numeric_limits<long long>::min();
static int sgn128(Int128 v) { return (v > 0) - (v < 0); }
// 愚直な ⌊a/b⌋（引き算で余りを 0 以上にそろえる）
static Int128 bf_floor(Int128 a, Int128 b) {
    if (b < 0) a = -a, b = -b;
    Int128 q = a / b, r = a - q * b;
    while (r < 0) r += b, --q;
    while (r >= b) r -= b, ++q;
    return q;
}
static Int128 bf_ceil(Int128 a, Int128 b) { return -bf_floor(-a, b); }
static unsigned long long bf_gcd(unsigned long long a, unsigned long long b) { return b ? bf_gcd(b, a % b) : a; }
static Int128 bf_gcd128(Int128 a, Int128 b) { return b ? bf_gcd128(b, a % b) : a; }
// y^k を ±2^100 で止める（|y| < 2^64）
static Int128 bf_pow_sat(Int128 y, long long k) {
    const Int128 SAT = (Int128)1 << 100;
    Int128 r = 1;
    const Int128 ay = y < 0 ? -y : y;
    for (long long i = 0; i < k; ++i) {
        if (r == 0) return 0;
        const Int128 ar = r < 0 ? -r : r;
        if (ay > 1 && ar > SAT / ay) return ((r < 0) != (y < 0)) ? -SAT : SAT;
        r *= y;
        if (ay <= 1 && (r == 1 || r == 0)) return (y == -1 && (k - i - 1) % 2) ? -r : r;
        if (ay <= 1 && r == -1) return ((k - i - 1) % 2) ? 1 : -1;
    }
    return r;
}
static long long bf_mod_pow(long long a, unsigned long long e, long long m) {  // 別の書き方（右から左ではなく左から右）
    Int128 r = 1 % m, x = ((Int128)a % m + m) % m;
    for (int bit = 63; bit >= 0; --bit) {
        r = r * r % m;
        if ((e >> bit) & 1) r = r * x % m;
    }
    return (long long)r;
}

static void test_division(std::mt19937_64& rng, int& cases) {
    auto check = [&](long long a, long long b) {
        if (b == 0 || (a == LMIN && b == -1)) return;
        ++cases;
        MATHB_CHECK(floor_div(a, b) == (long long)bf_floor(a, b));
        MATHB_CHECK(ceil_div(a, b) == (long long)bf_ceil(a, b));
        if (b != LMIN) {
            const long long r = mod_floor(a, b), ab = b < 0 ? -b : b;
            MATHB_CHECK(0 <= r && r < ab && ((Int128)a - r) % ab == 0);
        }
        // 四捨五入：-1/2 ≤ a/b - q < 1/2（ちょうど 1/2 のときは大きい方）
        const long long q = round_div(a, b);
        Int128 x = a, y = b;
        if (y < 0) x = -x, y = -y;
        const Int128 diff2 = 2 * (x - (Int128)q * y);  // 2(a/b - q)·y
        MATHB_CHECK(-y <= diff2 && diff2 < y);
    };
    for (long long a = -40; a <= 40; ++a)
        for (long long b = -12; b <= 12; ++b) check(a, b);
    for (int it = 0; it < 200000 * TEST_SCALE; ++it) check(rnd_ll(rng), rnd_ll(rng));
    MATHB_CHECK(floor_div(-7, 2) == -4 && ceil_div(-7, 2) == -3 && mod_floor(-7, 3) == 2 && round_div(-5, 2) == -2 && round_div(5, 2) == 3);
    MATHB_CHECK(floor_div(LMIN, 1) == LMIN && ceil_div(LMAX, -1) == -LMAX && mod_floor(LMIN, 7) == (long long)(((Int128)LMIN % 7 + 7) % 7));
}

static void test_multiples(std::mt19937_64& rng, int& cases) {
    // 小さい範囲は数える
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long l = rnd(rng, -60, 60), r = rnd(rng, -60, 60), d = rnd(rng, 1, 13), a = rnd(rng, -30, 30);
        long long c1 = 0, c2 = 0, c3 = 0, c4 = 0;
        for (long long x = l; x < r; ++x) c1 += (x % d == 0), c3 += (((x - a) % d) == 0);
        for (long long x = l; x <= r; ++x) c2 += (x % d == 0), c4 += (((x - a) % d) == 0);
        MATHB_CHECK(count_multiples(l, r, d) == c1 && count_multiples_closed(l, r, d) == c2);
        MATHB_CHECK(count_congruent(l, r, a, d) == c3 && count_congruent_closed(l, r, a, d) == c4);
    }
    // 大きい範囲は別の式（⌊(r-1)/d⌋ - ⌊(l-1)/d⌋）と性質で
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        long long l = rnd_ll(rng), r = rnd_ll(rng);
        if (l > r) std::swap(l, r);
        while ((Int128)r - l > ((Int128)1 << 62)) l /= 2, r /= 2;  // 個数が long long に収まる範囲
        long long d = rnd_ll(rng);
        if (d <= 0) d = d == LMIN ? 3 : -d + (d == 0);
        MATHB_CHECK((Int128)count_multiples(l, r, d) == bf_floor((Int128)r - 1, d) - bf_floor((Int128)l - 1, d));
        MATHB_CHECK((Int128)count_multiples_closed(l, r, d) == bf_floor(r, d) - bf_floor((Int128)l - 1, d));
        const long long a = rnd_ll(rng);
        MATHB_CHECK((Int128)count_congruent(l, r, a, d) == bf_floor((Int128)r - 1 - a, d) - bf_floor((Int128)l - 1 - a, d));
        MATHB_CHECK((Int128)count_congruent_closed(l, r, a, d) == bf_floor((Int128)r - a, d) - bf_floor((Int128)l - 1 - a, d));
        // 倍数・合同な値の切り上げ・切り捨て
        const long long x = rnd_ll(rng) / 2;
        d = d % (1LL << 61) + 1;  // 答えが long long に収まる範囲
        const long long cm = ceil_multiple(x, d), fm = floor_multiple(x, d);
        MATHB_CHECK(cm % d == 0 && cm >= x && (Int128)cm - d < x && fm % d == 0 && fm <= x && (Int128)fm + d > x);
        const long long cc = ceil_congruent(x, a, d), fc = floor_congruent(x, a, d);
        MATHB_CHECK(((Int128)cc - a) % d == 0 && cc >= x && (Int128)cc - d < x && ((Int128)fc - a) % d == 0 && fc <= x && (Int128)fc + d > x);
    }
}

static void test_roots(std::mt19937_64& rng, int& cases) {
    auto check = [&](long long n, int k) {
        if (n < 0) return;
        ++cases;
        const long long x = kth_root_floor(n, k), c = kth_root_ceil(n, k);
        MATHB_CHECK(x >= 0 && bf_pow_sat(x, k) <= n && bf_pow_sat((Int128)x + 1, k) > n);
        MATHB_CHECK(c >= 0 && bf_pow_sat(c, k) >= n && (c == 0 || bf_pow_sat((Int128)c - 1, k) < n));
        MATHB_CHECK(is_kth_power(n, k) == (bf_pow_sat(x, k) == n));
        if (k == 2) MATHB_CHECK(is_square(n) == (bf_pow_sat(x, 2) == n) && x == isqrt(n) && c == isqrt_ceil(n));
    };
    for (long long n = 0; n <= 3000; ++n)
        for (int k = 1; k <= 13; ++k) check(n, k);
    for (int it = 0; it < 30000 * TEST_SCALE; ++it) {
        const int k = (int)rnd(rng, 1, 66);
        long long n = rnd_ll(rng);
        if (n < 0) n = n == LMIN ? LMAX : -n;
        check(n, k);
        // 累乗数の前後
        const long long base = kth_root_floor(n, k);
        for (Int128 y = std::max(0LL, base - 1); y <= (Int128)base + 1; ++y) {
            const Int128 p = bf_pow_sat(y, k);
            if (p > LMAX) continue;
            check((long long)p, k);
            if (p > 0) check((long long)p - 1, k);
            if (p < LMAX) check((long long)p + 1, k);
        }
    }
    check(LMAX, 2), check(LMAX, 3), check(LMAX, 62), check(LMAX, 63), check(LMAX, 64);
    MATHB_CHECK(kth_root_floor(1000000000000000000LL, 3) == 1000000 && kth_root_floor(999999999999999999LL, 3) == 999999);
    MATHB_CHECK(kth_root_floor(LMAX, 62) == 2 && kth_root_floor(LMAX, 63) == 1 && kth_root_floor(4611686018427387904LL, 62) == 2);
    MATHB_CHECK(!is_square(-4) && is_square(0) && is_square(1) && !is_square(2) && is_kth_power(0, 5) && !is_kth_power(-8, 3));
}

static void test_pow(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 200000 * TEST_SCALE; ++it, ++cases) {
        long long a = rnd_ll(rng), b = rnd_ll(rng);
        MATHB_CHECK(mul_overflows(a, b) == ((Int128)a * b > LMAX || (Int128)a * b < LMIN));
        MATHB_CHECK(add_overflows(a, b) == ((Int128)a + b > LMAX || (Int128)a + b < LMIN));
        a = a < 0 ? (a == LMIN ? LMAX : -a) : a, b = b < 0 ? (b == LMIN ? LMAX : -b) : b;
        long long cap = rnd_ll(rng);
        cap = cap < 0 ? (cap == LMIN ? 0 : -cap) : cap;
        MATHB_CHECK((Int128)mul_cap(a, b, cap) == std::min<Int128>((Int128)a * b, cap));
        MATHB_CHECK((Int128)add_cap(a, b, cap) == std::min<Int128>((Int128)a + b, cap));
        // 累乗
        const long long base = it % 3 == 0 ? rnd(rng, -12, 12) : (it % 3 == 1 ? rnd(rng, -3100000000LL, 3100000000LL) : rnd_ll(rng));
        const long long e = it % 4 == 0 ? rnd(rng, 0, 70) : (it % 4 == 1 ? rnd(rng, 0, 3) : rnd(rng, 0, LMAX));
        const Int128 p = bf_pow_sat(base, std::min<long long>(e, 200));
        const long long ab = base < 0 ? (base == LMIN ? LMAX : -base) : base;
        const Int128 pa = bf_pow_sat(ab, std::min<long long>(e, 200));
        MATHB_CHECK((Int128)pow_cap(ab, e, cap) == std::min<Int128>(pa, cap));
        if (p >= LMIN && p <= LMAX && e <= 200) MATHB_CHECK(pow_int(base, e) == (long long)p);
    }
    MATHB_CHECK(pow_int(-2, 63) == LMIN && pow_int(2, 62) == 4611686018427387904LL && pow_int(0, 0) == 1 && pow_int(-1, LMAX) == -1);
    MATHB_CHECK(pow_cap(0, 0, 5) == 1 && pow_cap(0, 0, 0) == 0 && pow_cap(2, 1000000000000000000LL, LMAX) == LMAX && pow_cap(1, LMAX, 7) == 1);
    MATHB_CHECK(pow_cap(10, 18, LMAX) == 1000000000000000000LL && pow_cap(10, 19, LMAX) == LMAX);
    for (long long n = 1; n <= 5000; ++n)
        for (long long b = 2; b <= 12; ++b) {
            int f = 0, c = 0;
            for (long long p = 1; p * b <= n; p *= b) ++f;
            for (long long p = 1; p < n; p *= b) ++c;
            MATHB_CHECK(ilog_floor(n, b) == f && ilog_ceil(n, b) == c);
            ++cases;
        }
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        long long n = rnd_ll(rng);
        n = n <= 0 ? (n <= -LMAX ? LMAX : -n + 1) : n;
        const long long b = it % 2 ? rnd(rng, 2, 40) : (rnd_ll(rng) & LMAX) | 2;
        const int f = ilog_floor(n, b), c = ilog_ceil(n, b);
        MATHB_CHECK(bf_pow_sat(b, f) <= n && bf_pow_sat(b, f + 1) > n && bf_pow_sat(b, c) >= n && (c == 0 || bf_pow_sat(b, c - 1) < n));
    }
    for (int k = 0; k <= 18; ++k) MATHB_CHECK((Int128)pow10_ll(k) == bf_pow_sat(10, k));
    MATHB_CHECK(ilog_floor(LMAX, 2) == 62 && ilog_floor(999, 10) == 2 && ilog_floor(1000, 10) == 3 && ilog_ceil(1001, 10) == 4);
}

static void test_gcd(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 200000 * TEST_SCALE; ++it, ++cases) {
        long long a = rnd_ll(rng), b = rnd_ll(rng);
        if (it % 5 == 0) b = a / rnd(rng, 3, 100) * rnd(rng, -3, 3);  // 公約数の多い組
        const unsigned long long ua = a < 0 ? 0ULL - (unsigned long long)a : (unsigned long long)a, ub = b < 0 ? 0ULL - (unsigned long long)b : (unsigned long long)b;
        const unsigned long long g = bf_gcd(ua, ub);
        if (g > (unsigned long long)LMAX) continue;  // gcd = 2^63（LLONG_MIN と 0 など）は収まらない
        const ExtGcd e = ext_gcd(a, b);
        MATHB_CHECK((unsigned long long)e.g == g && (Int128)a * e.x + (Int128)b * e.y == e.g);
        if (g) {
            MATHB_CHECK((Int128)std::llabs(e.x) <= std::max<Int128>(1, (Int128)ub / g) && (Int128)std::llabs(e.y) <= std::max<Int128>(1, (Int128)ua / g));
        } else {
            MATHB_CHECK(e.x == 0 && e.y == 0);
        }
        // lcm
        const long long pa = (long long)(ua & (unsigned long long)LMAX), pb = (long long)(ub & (unsigned long long)LMAX);
        long long cap = rnd_ll(rng);
        cap = cap < 0 ? (cap == LMIN ? 0 : -cap) : cap;
        const unsigned long long gg = bf_gcd(pa, pb);
        const Int128 l = (pa == 0 || pb == 0) ? 0 : (Int128)(pa / gg) * pb;
        MATHB_CHECK((Int128)lcm_cap(pa, pb, cap) == std::min<Int128>(l, cap));
    }
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 8);
        std::vector<long long> v(n);
        for (auto& x : v) x = it % 2 ? rnd(rng, -60, 60) : rnd_ll(rng) / (long long)rnd(rng, 1, 1000000);
        unsigned long long g = 0;
        Int128 l = 1;
        bool zero = false;
        const Int128 SAT = (Int128)1 << 100;
        for (long long x : v) {
            const unsigned long long u = x < 0 ? 0ULL - (unsigned long long)x : (unsigned long long)x;
            g = bf_gcd(g, u);
            if (u == 0) zero = true;
            else if (l < SAT) {
                const Int128 t = l / bf_gcd128(l, (Int128)u);
                l = t > SAT / (Int128)u ? SAT : std::min<Int128>(SAT, t * (Int128)u);
            }
        }
        long long cap = rnd(rng, 0, it % 3 ? 1000000 : LMAX);
        MATHB_CHECK((unsigned long long)gcd_of(v) == g);
        MATHB_CHECK((Int128)lcm_of_cap(v, cap) == (zero ? 0 : std::min<Int128>(l, cap)));
        if (it % 2) {  // vector<int> でも
            std::vector<int> vi;
            for (long long x : v) vi.push_back((int)(x % 1000));
            unsigned long long gi = 0;
            for (int x : vi) gi = bf_gcd(gi, (unsigned long long)(x < 0 ? -x : x));
            MATHB_CHECK((unsigned long long)gcd_of(vi) == gi);
        }
        // try_lcm：2 つの数の lcm を 128 ビットで
        const long long a = it % 3 ? rnd(rng, -100000, 100000) : rnd_ll(rng), b = it % 3 ? rnd(rng, -100000, 100000) : rnd_ll(rng);
        const UInt128 ua = a < 0 ? (UInt128)0 - (UInt128)(Int128)a : (UInt128)a, ub = b < 0 ? (UInt128)0 - (UInt128)(Int128)b : (UInt128)b;
        const Maybe<long long> tl = try_lcm(a, b);
        if (ua == 0 || ub == 0) MATHB_CHECK(tl.ok && tl.value == 0);
        else {
            const UInt128 L = ua / bf_gcd128((Int128)ua, (Int128)ub) * ub;
            MATHB_CHECK(tl.ok == (L <= (UInt128)LMAX) && (!tl.ok || (UInt128)tl.value == L));
        }
    }
    MATHB_CHECK(gcd_of({12, -18, 30}) == 6 && lcm_of_cap({4, 6}, 100) == 12 && gcd_of(std::vector<int>{}) == 0 && lcm_of_cap(std::vector<unsigned>{4u, 6u}, 5) == 5);
    MATHB_CHECK(!try_lcm(LMAX, LMAX - 1).ok && try_lcm(LMIN + 1, 1).value == LMAX && try_lcm(-4, 6).value == 12);
    cases += 2;
}

static void test_mod(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        long long m;
        switch (it % 5) {
            case 0: m = rnd(rng, 1, 30); break;
            case 1: m = rnd(rng, 1, 4294967296LL); break;
            case 2: m = rnd(rng, 4294967297LL, LMAX); break;
            case 3: m = it % 2 ? 998244353 : 1000000007; break;
            default: m = LMAX - rnd(rng, 0, 3); break;
        }
        const long long a = rnd_ll(rng), b = rnd_ll(rng);
        MATHB_CHECK(mod_mul(a, b, m) == (long long)((((Int128)a * b) % m + m) % m));
        const unsigned long long e = it % 3 ? (unsigned long long)rnd(rng, 0, 20) : rng() >> 1;
        const long long p = mod_pow(a, (long long)e, m);
        MATHB_CHECK(p == bf_mod_pow(a, e, m));
        const unsigned long long g = bf_gcd((unsigned long long)mod_floor(a, m), (unsigned long long)m);
        const Maybe<long long> inv = try_mod_inv(a, m);
        MATHB_CHECK(inv.ok == (g == 1));
        if (inv.ok) {
            MATHB_CHECK(0 <= inv.value && inv.value < m && mod_mul(a, inv.value, m) == 1 % m && mod_inv(a, m) == inv.value);
            MATHB_CHECK(mod_mul(mod_pow(a, -(long long)e, m), p, m) == 1 % m);
        }
    }
    MATHB_CHECK(mod_pow(0, 0, 7) == 1 && mod_pow(5, 0, 1) == 0 && mod_pow(-2, 3, 5) == 2 && mod_pow(3, -1, 7) == 5);
    MATHB_CHECK(mod_pow(2, LMAX, LMAX) == bf_mod_pow(2, (unsigned long long)LMAX, LMAX) && mod_inv(LMAX - 1, LMAX) == LMAX - 1);
    // 等比数列の和
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long m = it % 3 == 0 ? rnd(rng, 1, 50) : (it % 3 == 1 ? rnd(rng, 1, 4294967296LL) : rnd(rng, 1, LMAX));
        const long long r = rnd_ll(rng), n = rnd(rng, 0, 300);
        Int128 s = 0, pw = 1 % m;
        const Int128 rr = ((Int128)r % m + m) % m;
        for (long long i = 0; i < n; ++i) s = (s + pw) % m, pw = pw * rr % m;
        MATHB_CHECK(geom_sum_mod(r, n, m) == (long long)s);
        // 大きい n：S(2n) = S(n)(1 + r^n)、S(n+1) = 1 + r S(n)
        const long long N = rnd(rng, 0, LMAX / 2 - 1);
        const long long sn = geom_sum_mod(r, N, m), rn = mod_pow(r, N, m);
        MATHB_CHECK(geom_sum_mod(r, 2 * N, m) == mod_mul(sn, (long long)(((Int128)rn + 1) % m), m));
        MATHB_CHECK(geom_sum_mod(r, N + 1, m) == (long long)(((Int128)mod_mul(r, sn, m) + 1) % m));
    }
    MATHB_CHECK(geom_sum_mod(2, 10, 1000000) == 1023 && geom_sum_mod(1, 1000000000000000000LL, 998244353) == 1000000000000000000LL % 998244353);
}

static void test_sums(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long a = it % 2 ? rnd(rng, -1000, 1000) : rnd_ll(rng), d = it % 2 ? rnd(rng, -1000, 1000) : rnd_ll(rng), n = rnd(rng, 0, 60);
        Int128 s = 0;
        for (long long i = 0; i < n; ++i) s += (Int128)a + (Int128)i * d;
        MATHB_CHECK(arith_sum(a, d, n) == s);
        const long long m = it % 3 ? rnd(rng, 1, 1000000007) : rnd(rng, 1, LMAX);
        MATHB_CHECK(arith_sum_mod(a, d, n, m) == (long long)((s % m + m) % m));
        const long long l = it % 2 ? rnd(rng, -100, 100) : rnd_ll(rng) / 4, r = l + rnd(rng, -3, 80);
        Int128 t = 0, tc = 0;
        for (long long x = l; x < r; ++x) t += x;
        for (long long x = l; x <= r; ++x) tc += x;
        MATHB_CHECK(sum_range(l, r) == t && sum_range_closed(l, r) == tc);
        MATHB_CHECK(sum_range_mod(l, r, m) == (long long)((t % m + m) % m) && sum_range_closed_mod(l, r, m) == (long long)((tc % m + m) % m));
    }
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {  // 端の値：sum_range（Int128）と
        long long l = rnd_ll(rng), r = rnd_ll(rng);
        if (l > r) std::swap(l, r);
        const long long m = it % 2 ? rnd(rng, 1, 1000) : rnd(rng, 1, LMAX);
        if ((Int128)r - l < ((Int128)1 << 62)) {
            const Int128 t = sum_range(l, r);
            MATHB_CHECK(sum_range_mod(l, r, m) == (long long)((t % m + m) % m));
        }
        MATHB_CHECK(sum_range_closed_mod(l, r, m) == (long long)(((Int128)sum_range_mod(l, r, m) + mod_floor(r, m)) % m));
    }
    MATHB_CHECK(sum_range_closed_mod(LMIN, LMAX, 1000000007) == (long long)((-((Int128)1 << 63) % 1000000007 + 1000000007) % 1000000007));
    // 大きい n は式どうしで
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 4000000000000000000LL), a = rnd(rng, -1000000000000000000LL, 1000000000000000000LL), d = rnd(rng, -10, 10);
        const Int128 expect = (Int128)n * a + (Int128)d * ((Int128)n * (n - 1) / 2);
        MATHB_CHECK(arith_sum(a, d, n) == expect);
        const long long m = rnd(rng, 1, LMAX);
        MATHB_CHECK(arith_sum_mod(a, d, n, m) == (long long)((expect % m + m) % m));
        const long long dm = rnd_ll(rng), nm = rnd(rng, 0, LMAX);
        // n(n-1)/2 を 2 つに分けて計算した値と比べる
        const Int128 h = nm % 2 == 0 ? (Int128)(nm / 2) % m * ((nm - (nm > 0)) % m) : (Int128)(nm % m) * (((nm - 1) / 2) % m);
        MATHB_CHECK(arith_sum_mod(0, dm, nm, m) == (long long)(((Int128)mod_floor(dm, m) * (h % m)) % m));
    }
    MATHB_CHECK(sum_range(LMIN, LMAX) == -((Int128)1 << 64) + 1 && sum_range_closed(LMIN, LMAX) == -((Int128)1 << 63));
    MATHB_CHECK(sum_range_closed(1, 100) == 5050 && sum_range(5, 5) == 0 && sum_range_closed(3, 2) == 0 && arith_sum(1, 2, 3) == 9);
}

static void test_frac(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 200000 * TEST_SCALE; ++it, ++cases) {
        long long a = rnd_ll(rng), b = rnd_ll(rng), c = rnd_ll(rng), d = rnd_ll(rng);
        if (b == 0) b = 1;
        if (d == 0) d = -1;
        if (it % 4 == 0) {  // 等しい分数を作る
            a /= 4, b /= 4;
            if (b == 0) b = 1;
            const long long k = rnd(rng, 1, 3) * ((rng() & 1) ? 1 : -1);
            c = a * k, d = b * k;
        }
        // 愚直：符号を b·d の符号で合わせて a·d と c·b を比べる
        const int s = sgn128((Int128)b) * sgn128((Int128)d);
        const Int128 l = (Int128)a * d * s, r = (Int128)c * b * s;
        MATHB_CHECK(frac_compare(a, b, c, d) == (l > r) - (l < r));
        MATHB_CHECK(frac_compare(c, d, a, b) == -frac_compare(a, b, c, d));
    }
    MATHB_CHECK(frac_compare(1, 3, 2, 6) == 0 && frac_compare(-1, 2, 1, -2) == 0 && frac_compare(1, -3, 0, 5) == -1 && frac_compare(LMAX, LMAX - 1, LMAX - 1, LMAX - 2) == -1);
}

static void test_misc(int& cases) {
    const Maybe<long long> m1 = try_mod_inv(4, 8), m2 = try_mod_inv(3, 8);
    MATHB_CHECK(!m1.ok && m2.ok && m2.value == 3);
    auto [ok, v] = try_mod_inv(2, 9);
    MATHB_CHECK(ok && v == 5);
    const Maybe<int> none;
    MATHB_CHECK(!none.ok && none.value == 0);
    ++cases;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_division(rng, cases);
    test_multiples(rng, cases);
    test_roots(rng, cases);
    test_pow(rng, cases);
    test_gcd(rng, cases);
    test_mod(rng, cases);
    test_sums(rng, cases);
    test_frac(rng, cases);
    test_misc(cases);
    std::cout << "[math_base] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[math_base] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 2000000;
    std::vector<long long> xs(n), ys(n);
    for (int i = 0; i < n; ++i) xs[i] = rnd_ll(rng), ys[i] = rnd_ll(rng) | 1, ys[i] = ys[i] == -1 ? 3 : ys[i];  // LLONG_MIN / -1 を避ける
    auto t0 = std::chrono::steady_clock::now();
    unsigned long long s1 = 0;
    for (int i = 0; i < n; ++i) s1 += (unsigned long long)(floor_div(xs[i], ys[i]) + mod_floor(xs[i], ys[i]));
    const long long t_div = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s2 = 0;
    for (int i = 0; i < n / 4; ++i) s2 += (unsigned long long)(mod_pow(xs[i], ys[i] & LMAX, 1000000007));
    const long long t_pow32 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s3 = 0;
    for (int i = 0; i < n / 4; ++i) s3 += (unsigned long long)(mod_pow(xs[i], ys[i] & LMAX, 4611686018427387847LL));
    const long long t_pow64 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s4 = 0;
    for (int i = 0; i < n; ++i) s4 += (unsigned long long)(kth_root_floor(xs[i] & LMAX, 3));
    const long long t_cbrt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s5 = 0;
    for (int i = 0; i < n; ++i) s5 += (unsigned long long)(ext_gcd(xs[i], ys[i]).x);
    const long long t_gcd = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s6 = 0;
    for (int i = 0; i < n / 4; ++i) s6 += (unsigned long long)(geom_sum_mod(xs[i], ys[i] & LMAX, 998244353));
    const long long t_geom = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s7 = 0;
    for (int i = 0; i < n; ++i) s7 += (unsigned long long)(pow_cap(xs[i] & 1023, ys[i] & 63, 1000000000000000000LL));
    const long long t_cap = ms(t0);
    std::cout << "[math_base] benchmark: floor_div + mod_floor × " << n << " " << t_div << " ms, kth_root_floor(·, 3) × " << n << " " << t_cbrt
              << " ms, ext_gcd × " << n << " " << t_gcd << " ms, pow_cap × " << n << " " << t_cap << " ms\n"
              << "[math_base] benchmark: mod_pow（m = 10^9+7）× " << n / 4 << " " << t_pow32 << " ms, mod_pow（m ≈ 4.6×10^18）× " << n / 4 << " "
              << t_pow64 << " ms, geom_sum_mod × " << n / 4 << " " << t_geom << " ms  (check "
              << (s1 + s2 + s3 + s4 + s5 + s6 + s7) % 7 << ")\n";
}

}  // namespace mathb_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace mathb_demo {

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
#define SHOW(expr, comment) mathb_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    mathb_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("数学基盤（math_base） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  floor_div / ceil_div / mod_floor / round_div      負の数でも正しい切り捨て・切り上げ・余り・四捨五入\n"
              << "  ceil_multiple / floor_multiple / ceil_congruent / floor_congruent   x 以上・以下で最も近い倍数・合同な数\n"
              << "  count_multiples(_closed) / count_congruent(_closed)                  範囲の中の倍数・合同な数の個数\n"
              << "  isqrt / isqrt_ceil / kth_root_floor / kth_root_ceil / is_square / is_kth_power   累乗根（誤差なし）\n"
              << "  pow_int / pow_cap / mul_cap / add_cap / mul_overflows / add_overflows  累乗と溢れ\n"
              << "  ilog_floor / ilog_ceil / pow10_ll                 対数の切り捨て・切り上げ / 10 の累乗\n"
              << "  ext_gcd / lcm_cap / try_lcm / gcd_of / lcm_of_cap  拡張ユークリッド / 溢れない lcm / 全部の gcd・lcm\n"
              << "  mod_mul / mod_pow / mod_inv / try_mod_inv         64 ビットの mod（合成数の mod でも）\n"
              << "  geom_sum_mod / arith_sum / arith_sum_mod / sum_range(_closed)(_mod)  等比・等差数列の和\n"
              << "  frac_compare / Maybe<T>                           分数の大小（誤差なし）/ ないことがある値\n"
              << "  （Int128・isqrt・parse_fixed は 00_int_common。cout << Int128 もそのまま使える）\n";

    api("【1】 floor_div(a, b) / ceil_div(a, b) / mod_floor(a, m) / round_div(a, b)", "long long",
        "⌊a/b⌋・⌈a/b⌉・0 以上の余り・四捨五入。C++ の / と % は 0 に向かって切り捨てるので、負の数ではこちらを使う", "O(1)");
    SHOW(-7 / 2, "C++ の割り算（0 に向かって切り捨て）");
    SHOW(floor_div(-7, 2), "⌊-3.5⌋ = -4");
    SHOW(ceil_div(-7, 2), "⌈-3.5⌉ = -3");
    SHOW((-7 + 3 - 1) / 3, "よくある切り上げ (a + b - 1) / b は負の数で誤る（⌈-7/3⌉ = -2）");
    SHOW(ceil_div(-7, 3), "正しい切り上げ");
    SHOW(-7 % 3, "C++ の余り（負になる）");
    SHOW(mod_floor(-7, 3), "0 以上の余り");
    SHOW(round_div(7, 2), "3.5 → 4");
    SHOW(round_div(-7, 2), "-3.5 → -3（.5 は大きい方へ）");

    api("【2】 count_multiples(l, r, d) / count_multiples_closed(L, R, d) / count_congruent(l, r, a, m)", "long long",
        "[l, r)（_closed なら [L, R]）にある d の倍数・x ≡ a (mod m) の個数。負の範囲でもよい", "O(1)");
    SHOW(count_multiples_closed(1, 100, 7), "1〜100 の 7 の倍数");
    SHOW(count_multiples(-10, 10, 3), "[-10, 10) の 3 の倍数（-9, -6, -3, 0, 3, 6, 9）");
    SHOW(count_congruent_closed(1, 100, 2, 7), "1〜100 で 7 で割って 2 余る数");
    {
        // ABC131 C: A 以上 B 以下で、C でも D でも割り切れない数
        const long long A = 4, B = 9, C = 2, D = 3;
        const long long cd = lcm_cap(C, D, B + 1);  // lcm が B を超えたら、どちらでも割り切れる数はない
        const long long ans = (B - A + 1) - count_multiples_closed(A, B, C) - count_multiples_closed(A, B, D) + count_multiples_closed(A, B, cd);
        note("ABC131 C（C でも D でも割り切れない数）: 全体 − C の倍数 − D の倍数 + lcm(C, D) の倍数");
        std::cout << "    A = 4, B = 9, C = 2, D = 3 → " << ans << "（5, 7）\n";
    }
    SHOW(ceil_multiple(-7, 3), "-7 以上で最小の 3 の倍数");
    SHOW(floor_multiple(-7, 3), "-7 以下で最大の 3 の倍数");
    SHOW(ceil_congruent(10, 2, 7), "10 以上で 7 で割って 2 余る最小の数");
    SHOW(floor_congruent(10, 2, 7), "10 以下で最大");

    api("【3】 isqrt(n) / kth_root_floor(n, k) / kth_root_ceil(n, k) / is_square(n) / is_kth_power(n, k)", "long long / bool",
        "誤差のない ⌊√n⌋・⌊n^(1/k)⌋・⌈n^(1/k)⌉。std::pow / std::cbrt は 10^18 付近で 1 ずれることがある", "O(k)");
    SHOW(kth_root_floor(1000000000000000000LL, 3), "10^18 の立方根");
    SHOW(kth_root_floor(999999999999999999LL, 3), "10^18 − 1 → 999999（cbrt だと 1000000 になりうる）");
    SHOW(kth_root_ceil(10, 3), "x^3 ≥ 10 となる最小の x");
    SHOW(is_square(144), "12^2");
    SHOW(is_kth_power(1024, 5), "4^5");
    {
        // ABC227 C: A ≤ B ≤ C かつ ABC ≤ N となる組の数
        const long long N = 100;
        long long cnt = 0;
        for (long long a = 1; a <= kth_root_floor(N, 3); ++a)
            for (long long b = a; b <= isqrt(N / a); ++b) cnt += N / (a * b) - b + 1;
        note("ABC227 C（ABC ≤ N の数え上げ）: A ≤ ∛N、B ≤ √(N/A) の範囲で C の個数 ⌊N/(AB)⌋ − B + 1 を足す");
        std::cout << "    N = 100 → " << cnt << "\n";
    }

    api("【4】 pow_cap(a, e, cap) / mul_cap / add_cap / pow_int(a, e) / mul_overflows / add_overflows", "long long / bool",
        "溢れずに「cap を超えるか」を調べる（_cap は min(値, cap)）。pow_int は厳密な値（溢れたら止める）", "O(log e) / O(1)");
    SHOW(pow_cap(10, 30, 1000000000000000001LL), "10^30 は 10^18 + 1 で止まる → 10^18 を超える");
    SHOW(pow_int(-3, 5), "厳密な値");
    SHOW(mul_overflows(3000000000LL, 4000000000LL), "1.2×10^19 は long long に収まらない");
    {
        // ABC169 B: A1 × … × AN。10^18 を超えたら -1（0 があれば 0）
        const std::vector<long long> A = {1000000000LL, 1000000000LL, 10};
        const long long LIM = 1000000000000000000LL;
        long long prod = 1;
        bool zero = false;
        for (long long a : A) zero |= (a == 0);
        for (long long a : A) prod = mul_cap(prod, a, LIM + 1);
        note("ABC169 B（積が 10^18 を超えたら -1）: mul_cap(積, a, 10^18 + 1) で溢れずに掛ける。0 があれば先に 0");
        std::cout << "    {10^9, 10^9, 10} → " << (zero ? 0 : (prod > LIM ? -1 : prod)) << "\n";
    }
    {
        // 典型90 020: log2(a) < b·log2(c) ⇔ a < c^b
        const long long a = 4, b = 3, c = 2;
        note("典型90 020（log2 a < b log2 c）: a < c^b を pow_cap(c, b, a + 1) > a で（実数の log を使わない）");
        std::cout << "    a = 4, b = 3, c = 2 → " << (pow_cap(c, b, a + 1) > a ? "Yes" : "No") << "\n";
    }

    api("【5】 ilog_floor(n, b) / ilog_ceil(n, b) / pow10_ll(k)", "int / long long", "b^k ≤ n となる最大の k / b^k ≥ n となる最小の k / 10^k", "O(log n)");
    SHOW(ilog_floor(999, 10) + 1, "10 進の桁数");
    SHOW(ilog_ceil(1000000, 2), "2^k ≥ 10^6 となる最小の k（ダブリングの段数など）");
    SHOW(pow10_ll(18), "10^18");

    api("【6】 ext_gcd(a, b) / lcm_cap(a, b, cap) / try_lcm(a, b) / gcd_of(v) / lcm_of_cap(v, cap)", "ExtGcd{g, x, y} / long long / Maybe",
        "a x + b y = gcd の解（負の数でも）/ 溢れない lcm / 全部の gcd・lcm（gcd・lcm 1 つずつなら std::gcd / std::lcm）", "O(log)");
    RUN(ExtGcd e = ext_gcd(240, 46));
    SHOW(e.g, "gcd(240, 46)");
    SHOW(240 * e.x + 46 * e.y, "240 x + 46 y = 2 を満たす");
    SHOW(lcm_cap(1000000000000LL, 999999999999LL, 1000000000000000001LL), "lcm ≈ 10^24 は 10^18 + 1 で止まる（典型90 038 の Large 判定）");
    SHOW(lcm_cap(1000000000LL, 999999999LL, 1000000000000000001LL), "超えなければ lcm そのもの");
    SHOW(gcd_of({12, -18, 30}), "全部の gcd（絶対値で。vector<int> などでも）");
    SHOW(lcm_of_cap({4, 6, 10}, 1000), "全部の lcm");
    SHOW(try_lcm(1000000000000LL, 999999999999LL).ok, "lcm が long long に収まらない（Maybe）");
    SHOW(try_lcm(-4, 6).value, "収まれば lcm（絶対値で）");
    note("多くの数の lcm を mod で求めるときは、mod を取りながら lcm を取ってはいけない（02_factorize の lcm_mod_of）。");

    api("【7】 mod_pow(a, e, m) / mod_inv(a, m) / try_mod_inv(a, m) / mod_mul(a, b, m)", "long long / Maybe<long long>",
        "64 ビットの mod（m < 2^63）。逆元は合成数の mod でも gcd(a, m) = 1 なら求まる。mod が 2^31 未満の定数なら 00_modint の方が速い", "O(log e)");
    SHOW(mod_pow(2, 100, 1000000007), "2^100 mod (10^9+7)");
    SHOW(mod_pow(3, -1, 10), "3 の逆元 mod 10（3 × 7 = 21 ≡ 1）");
    SHOW(mod_inv(3, 10), "同じ");
    SHOW(try_mod_inv(4, 10).ok, "gcd(4, 10) = 2 なので逆元はない");
    SHOW(mod_mul(4000000000000000000LL, 3, 9000000000000000000LL), "1.2×10^19 mod 9×10^18（128 ビットで計算）");
    {
        // ABC228 E: M^(K^N) mod p（p = 998244353）。指数は p - 1 で割った余りに。ただし M が p の倍数なら 0
        const long long N = 2, K = 2, M = 2, p = 998244353;
        const long long ans = M % p == 0 ? 0 : mod_pow(M, mod_pow(K, N, p - 1), p);
        note("ABC228 E（M^(K^N) mod p）: フェルマーで指数を p − 1 で割った余りにする。M ≡ 0 のときは答えが 0（指数が 0 になっても 1 にしない）");
        std::cout << "    N = 2, K = 2, M = 2 → 2^4 = " << ans << "\n";
    }

    api("【8】 geom_sum_mod(r, n, m)", "long long", "1 + r + r^2 + … + r^(n-1) を mod m で。逆元を使わないので m が合成数・r - 1 が割れないときも正しい", "O(log n)");
    SHOW(geom_sum_mod(3, 4, 7), "ABC293 E: 1 + 3 + 9 + 27 = 40 ≡ 5 (mod 7)");
    {
        // ABC357 D: N を N 回つなげた数 mod 998244353 = N × (1 + 10^d + 10^2d + …)
        const long long N = 5, p = 998244353;
        const int d = ilog_floor(N, 10) + 1;
        note("ABC357 D（N を N 回つなげた数）: N × geom_sum_mod(10^d, N, p)（d は N の桁数。10^d は mod_pow(10, d, p)）");
        std::cout << "    N = 5 → " << mod_mul(N % p, geom_sum_mod(mod_pow(10, d, p), N, p), p) << "（55555）\n";
    }

    api("【9】 arith_sum(a, d, n) / arith_sum_mod(a, d, n, m) / sum_range(l, r) / sum_range_closed(L, R) / sum_range_mod(l, r, m)", "Int128 / long long",
        "等差数列の和（Int128 に収まらなければ止める。mod で欲しいなら _mod）。sum_range は l 以上 r 未満の整数の和", "O(1)");
    SHOW(arith_sum(1, 2, 5), "1 + 3 + 5 + 7 + 9");
    SHOW(sum_range_closed(1, 1000000000000LL), "1 + … + 10^12（Int128。cout でそのまま出せる）");
    SHOW(arith_sum_mod(1, 1, 1000000000000000000LL, 998244353), "1 + … + 10^18 を mod で");
    SHOW(sum_range_mod(100, 1000000000000000000LL, 998244353), "100 + … + (10^18 − 1) を mod で（引数は区間の端）");

    api("【10】 frac_compare(a, b, c, d)", "int（-1 / 0 / 1）", "a/b と c/d を誤差なしで比べる（分母は負でもよい。128 ビットで掛け算）", "O(1)");
    SHOW(frac_compare(1, 3, 333333333, 1000000000), "1/3 > 0.333333333");
    SHOW(frac_compare(2, -4, -1, 2), "-1/2 = -1/2");
    {
        // ABC308 C: 成功率 A/(A+B) の大きい順（同じなら番号の小さい順）
        std::vector<long long> A = {1, 3, 2}, B = {3, 1, 2};
        std::vector<int> idx = {0, 1, 2};
        std::stable_sort(idx.begin(), idx.end(), [&](int i, int j) { return frac_compare(A[i], A[i] + B[i], A[j], A[j] + B[j]) > 0; });
        note("ABC308 C（成功率で並べる）: stable_sort の比較に frac_compare(A_i, A_i + B_i, A_j, A_j + B_j) > 0");
        std::cout << "    → " << idx[0] + 1 << " " << idx[1] + 1 << " " << idx[2] + 1 << "\n";
    }
    note("何度も分数を足したり掛けたりするなら 08_fraction の Fraction。");

    api("【11】 Maybe<T>{ok, value}", "構造体", "ないことがある値。try_mod_inv などが返す。ok を確かめてから value を使う（構造化束縛でも）", "O(1)");
    RUN(auto [ok, inv] = try_mod_inv(7, 12));
    SHOW(ok, "7 と 12 は互いに素");
    SHOW(inv, "7 × 7 = 49 ≡ 1 (mod 12)");
}

}  // namespace mathb_demo

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
    mathb_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mathb_test::benchmark();    // 速度計測
    mathb_demo::demo();         // 使い方デモ
    return 0;
}
