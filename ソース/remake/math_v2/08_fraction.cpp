// =============================================================================
//  分数（Fraction：いつも既約・分母 > 0。誤差なしの四則演算・比較・切り捨て・小数表示・近い分数・あいだの最も簡単な分数・連分数）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 08 分数・大きな数
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  使い方
//      Fraction a(1, 3), b = 2;             // 1/3 と 2（long long。分子・分母が 10^18 を超えるなら Fraction128）
//      cout << a + b << '\n';               // 7/3
//      sort(v.begin(), v.end());            // 比率で誤差なく並べる（< が使える）
//  途中の計算は 128 ビット。結果が型に収まらなければメッセージを出して止める
//
//  機能一覧                                                                    計算量
//    Fraction x(n, d) / Fraction x = n    n/d（約分する。d ≠ 0）/ 整数。x.num・x.den で分子・分母    O(log)
//    + − × ÷ 単項− += −= ×= ÷=            整数と混ざった式も（x * 2 + 1 など）                       O(log)
//    == != < > <= >= / compare(x, y)      誤差なしの比較 / −1・0・1                                  O(1)（Fraction128 は O(log)）
//    x.floor() / x.ceil() / x.round()     ⌊x⌋ / ⌈x⌉ / 四捨五入（ちょうど .5 は 0 から遠い方）
//    x.sign() / x.abs() / x.reciprocal() / x.is_integer()
//    x.to_real()                          long double に
//    x.to_decimal(k)                      小数第 k 位まで（それより下は切り捨て）の文字列                 O(k)
//    cout << x / cin >> x                 "p/q"（整数なら "p"）で出す / "p/q"・"p"・"-1.25" を読む
//    unordered_map<Fraction, V>           std::hash を用意してある
//    best_rational_approx(p, q, max_den)  p/q にいちばん近い、分母が max_den 以下の分数（同じ近さなら分母、値の小さい方） O(log)
//    rational_floor(p, q, max_den) / rational_ceil(p, q, max_den)   分母が max_den 以下で p/q 以下の最大 / 以上の最小     O(log)
//    simplest_between(x, y)           x < r < y となる分数 r で分母が最小のもの（|分子| も最小。開区間）               O(log)
//    continued_fraction(p, q) / from_continued_fraction(cf)   連分数展開 [a0; a1, …]（vector<long long>）/ その値         O(log)
//    Fraction128                          分子・分母が Int128 のもの（使い方は同じ）
//  （分数 a/b と c/d の大小を 1 回だけ比べるなら 00_math_base の frac_compare でもよい）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    frac_test::random_test()  … 128 ビットの掛け算での検算・全部調べる愚直解との突き合わせ
//    frac_test::benchmark()    … 速度計測
//    frac_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
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
#ifndef GG_FRACTION
#define GG_FRACTION
// ---- 分数（数学第二世代）  FractionT<T>：いつも既約・分母 > 0。Fraction = long long、Fraction128 = Int128 ----
namespace math_internal {
inline Int128 frac_abs(Int128 a) { return a < 0 ? -a : a; }
inline int frac_ctz(UInt128 x) {  // x ≠ 0
    const unsigned long long lo = (unsigned long long)x;
    return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((unsigned long long)(x >> 64));
}
inline UInt128 frac_uabs(Int128 a) { return a < 0 ? (UInt128)0 - (UInt128)a : (UInt128)a; }  // −2^127 でも溢れない
inline Int128 frac_gcd(Int128 sa, Int128 sb) {  // gcd(|a|, |b|)（128 ビットは二進 GCD、64 ビットに入ったら速い方で。結果 < 2^127 のときに使う）
    UInt128 a = frac_uabs(sa), b = frac_uabs(sb);
    if (a == 0) return (Int128)b;
    if (b == 0) return (Int128)a;
    const int za = frac_ctz(a), zb = frac_ctz(b), k = za < zb ? za : zb;
    a >>= za, b >>= zb;
    while (true) {  // a, b は奇数
        if ((a >> 64) == 0 && (b >> 64) == 0) return (Int128)((UInt128)ugcd((unsigned long long)a, (unsigned long long)b) << k);
        if (a > b) std::swap(a, b);
        b -= a;
        if (b == 0) return (Int128)(a << k);
        b >>= frac_ctz(b);
    }
}
inline Int128 frac_mul(Int128 a, Int128 b) {
    Int128 r;
    if (__builtin_mul_overflow(a, b, &r)) fail("Fraction", "途中の値が 128 ビットに収まらない");
    return r;
}
inline Int128 frac_add(Int128 a, Int128 b) {
    Int128 r;
    if (__builtin_add_overflow(a, b, &r)) fail("Fraction", "途中の値が 128 ビットに収まらない");
    return r;
}
// a/b と c/d（b, d > 0）の大小（-1 / 0 / 1）。掛け算をしない（連分数のように整数部分を比べていく）
inline int frac_cmp_euclid(Int128 a, Int128 b, Int128 c, Int128 d) {
    int sign = 1;
    while (true) {
        Int128 q1 = a / b, r1 = a % b, q2 = c / d, r2 = c % d;  // 掛け算をしないので溢れない
        if (r1 < 0) q1 -= 1, r1 += b;                          // 床関数の商と 0 ≤ r1 < b に
        if (r2 < 0) q2 -= 1, r2 += d;
        if (q1 != q2) return q1 < q2 ? -sign : sign;
        if (r1 == 0 || r2 == 0) return r1 == r2 ? 0 : (r1 == 0 ? -sign : sign);
        a = b, b = r1, c = d, d = r2;  // r1/b と r2/d の比較は逆数 b/r1 と d/r2 の比較（向きが逆）
        sign = -sign;
    }
}
}  // namespace math_internal

template <class T>
struct FractionT {
    static_assert(std::is_same_v<T, long long> || std::is_same_v<T, Int128>, "FractionT は long long か Int128");
    T num, den;  // いつも既約で den > 0（直接書きかえない）

    FractionT() : num(0), den(1) {}
    template <class I, std::enable_if_t<std::is_integral_v<I> || std::is_same_v<I, Int128>, int> = 0>
    FractionT(I n) : num(0), den(1) {  // 整数から
        *this = make(to_i128(n), 1);
    }
    template <class I, class J, std::enable_if_t<(std::is_integral_v<I> || std::is_same_v<I, Int128>) && (std::is_integral_v<J> || std::is_same_v<J, Int128>), int> = 0>
    FractionT(I n, J d) : num(0), den(1) {  // n / d（d ≠ 0。約分する）
        if (d == 0) math_internal::fail("Fraction", "分母が 0");
        *this = make(to_i128(n), to_i128(d));
    }

    friend FractionT operator+(const FractionT& x, const FractionT& y) {
        using namespace math_internal;
        const Int128 g = frac_gcd(x.den, y.den);
        const Int128 t = frac_add(frac_mul(x.num, y.den / g), frac_mul(y.num, x.den / g));
        if (g == 1) return reduced(t, frac_mul(x.den, y.den));  // 分母が互いに素なら約分は要らない
        const Int128 g2 = frac_gcd(t, g);                      // gcd(t, 分母) = gcd(t, g)
        return reduced(t / g2, frac_mul(x.den / g, y.den / g2));
    }
    friend FractionT operator-(const FractionT& x, const FractionT& y) { return x + (-y); }
    friend FractionT operator*(const FractionT& x, const FractionT& y) {
        using namespace math_internal;
        const Int128 g1 = frac_gcd(x.num, y.den), g2 = frac_gcd(y.num, x.den);
        if (g1 == 0 || g2 == 0) return FractionT();  // どちらかが 0
        return reduced(frac_mul(x.num / g1, y.num / g2), frac_mul(x.den / g2, y.den / g1));  // 先に約分したので既約
    }
    friend FractionT operator/(const FractionT& x, const FractionT& y) { return x * y.reciprocal(); }
    FractionT operator-() const {
        FractionT r;
        r.num = -num, r.den = den;  // num は T の最小値にならないようにしてある
        return r;
    }
    FractionT operator+() const { return *this; }
    FractionT& operator+=(const FractionT& y) { return *this = *this + y; }
    FractionT& operator-=(const FractionT& y) { return *this = *this - y; }
    FractionT& operator*=(const FractionT& y) { return *this = *this * y; }
    FractionT& operator/=(const FractionT& y) { return *this = *this / y; }

    friend int compare(const FractionT& x, const FractionT& y) {  // −1 / 0 / 1
        if constexpr (std::is_same_v<T, long long>) {
            const Int128 l = (Int128)x.num * y.den, r = (Int128)y.num * x.den;
            return l < r ? -1 : (l > r ? 1 : 0);
        } else {
            const Int128 L = (Int128)1 << 63;  // 4 つとも 2^63 未満なら掛け算で
            if (-L < x.num && x.num < L && x.den < L && -L < y.num && y.num < L && y.den < L) {
                const Int128 l = x.num * y.den, r = y.num * x.den;
                return l < r ? -1 : (l > r ? 1 : 0);
            }
            return math_internal::frac_cmp_euclid(x.num, x.den, y.num, y.den);
        }
    }
    friend bool operator==(const FractionT& x, const FractionT& y) { return x.num == y.num && x.den == y.den; }
    friend bool operator!=(const FractionT& x, const FractionT& y) { return !(x == y); }
    friend bool operator<(const FractionT& x, const FractionT& y) { return compare(x, y) < 0; }
    friend bool operator>(const FractionT& x, const FractionT& y) { return compare(x, y) > 0; }
    friend bool operator<=(const FractionT& x, const FractionT& y) { return compare(x, y) <= 0; }
    friend bool operator>=(const FractionT& x, const FractionT& y) { return compare(x, y) >= 0; }

    T floor() const { return (T)math_internal::floor_div128(num, den); }  // ⌊x⌋
    T ceil() const { return (T)math_internal::ceil_div128(num, den); }    // ⌈x⌉
    T round() const {                                                    // 四捨五入（ちょうど .5 は 0 から遠い方）
        const T q = num / den, r = num % den;                             // 0 に向かって切り捨て
        const T ar = r < 0 ? -r : r;
        if (ar >= den - ar) return r < 0 ? q - 1 : q + 1;
        return q;
    }
    int sign() const { return num > 0 ? 1 : (num < 0 ? -1 : 0); }
    FractionT abs() const { return num < 0 ? -*this : *this; }
    FractionT reciprocal() const {  // 1/x（0 なら止める）
        if (num == 0) math_internal::fail("Fraction::reciprocal", "0 の逆数");
        FractionT r;
        r.num = num < 0 ? -den : den, r.den = num < 0 ? -num : num;
        return r;
    }
    bool is_integer() const { return den == 1; }
    long double to_real() const { return (long double)num / (long double)den; }
    std::string to_decimal(int k) const {  // 小数第 k 位まで（それより下は切り捨て。k = 0 なら整数部分だけ）
        if (k < 0) math_internal::fail("Fraction::to_decimal", "k ≥ 0 で使う");
        const Int128 a = math_internal::frac_abs(num);
        const Int128 ip = a / den;
        Int128 r = a % den;
        std::string digits;
        const UInt128 D = (UInt128)den;
        for (int i = 0; i < k; ++i) {  // 10 r を den で割る（r を 10 回足しながら引く。10 r が 128 ビットを超えても大丈夫）
            UInt128 acc = 0, R = (UInt128)r;
            int dg = 0;
            for (int t = 0; t < 10; ++t) {
                acc += R;  // acc < den、R < den なので 2^128 未満
                if (acc >= D) acc -= D, ++dg;
            }
            digits += (char)('0' + dg), r = (Int128)acc;
        }
        const bool zero = ip == 0 && digits.find_first_not_of('0') == std::string::npos;
        return (num < 0 && !zero ? "-" : "") + int128_to_string(ip) + (k ? "." + digits : "");
    }
    friend std::ostream& operator<<(std::ostream& os, const FractionT& x) {  // "p/q"（整数なら "p"）
        os << int128_to_string(x.num);
        if (x.den != 1) os << '/' << int128_to_string(x.den);
        return os;
    }
    friend std::istream& operator>>(std::istream& is, FractionT& x) {  // "p/q"・"p"・"-1.25" を読む
        std::string s;
        if (!(is >> s)) return is;
        x = parse(s);
        return is;
    }
    static FractionT parse(const std::string& s) {  // "p/q"・"p"・"-1.25"・".5" を読む
        using namespace math_internal;
        const size_t slash = s.find('/');
        if (slash != std::string::npos) return FractionT(parse_int128(s.substr(0, slash)), parse_int128(s.substr(slash + 1)));
        const size_t dot = s.find('.');
        if (dot == std::string::npos) return FractionT(parse_int128(s));
        const size_t st = (!s.empty() && (s[0] == '-' || s[0] == '+')) ? 1 : 0;
        Int128 n = 0, d = 1;
        int cnt = 0;  // 数字の数
        for (size_t i = st; i < s.size(); ++i) {
            if (i == dot) continue;
            if (s[i] < '0' || s[i] > '9') fail("Fraction::parse", "数の形が正しくない");
            n = frac_add(frac_mul(n, 10), s[i] - '0'), ++cnt;
            if (i > dot) d = frac_mul(d, 10);
        }
        if (cnt == 0) fail("Fraction::parse", "数字がない");
        return FractionT(st == 1 && s[0] == '-' ? -n : n, d);
    }

  private:
    template <class I>
    static Int128 to_i128(I v) {  // 整数を Int128 に（2^127 以上の UInt128 は止める）
        if constexpr (sizeof(I) == 16 && !std::is_same_v<I, Int128>) {
            if (v > (I)((UInt128)1 << 127) - 1) math_internal::fail("Fraction", "値が型に収まらない");
        }
        return (Int128)v;
    }
    static FractionT reduced(Int128 n, Int128 d) {  // 既約で d > 0 と分かっているもの（T に入るかだけ確かめる）
        if constexpr (std::is_same_v<T, long long>) {
            const Int128 L = (Int128)9223372036854775807LL;
            if (n > L || n < -L || d > L) math_internal::fail("Fraction", "値が long long に収まらない（Fraction128 を使う）");
        } else {
            const Int128 MIN128 = (Int128)((UInt128)1 << 127);
            if (n == MIN128) math_internal::fail("Fraction", "値が型に収まらない");
        }
        FractionT r;
        r.num = (T)n, r.den = (T)d;
        return r;
    }
    static FractionT make(Int128 n, Int128 d) {  // 約分して分母を正に。T に入らなければ止める
        using namespace math_internal;
        if (d == 0) fail("Fraction", "分母が 0");
        const Int128 MIN128 = (Int128)((UInt128)1 << 127);
        if (n == MIN128 || d == MIN128) fail("Fraction", "値が型に収まらない");
        if (d < 0) n = -n, d = -d;
        const Int128 g = frac_gcd(n, d);
        if (g > 1) n /= g, d /= g;
        if (n == 0) d = 1;
        if constexpr (std::is_same_v<T, long long>) {
            const Int128 L = (Int128)9223372036854775807LL;
            if (n > L || n < -L || d > L) fail("Fraction", "値が long long に収まらない（Fraction128 を使う）");
        }
        FractionT r;
        r.num = (T)n, r.den = (T)d;
        return r;
    }
};
using Fraction = FractionT<long long>;
using Fraction128 = FractionT<Int128>;
namespace std {
template <class T>
struct hash<FractionT<T>> {  // unordered_map<Fraction, V> でそのまま使える
    size_t operator()(const FractionT<T>& x) const {
        const unsigned long long a = (unsigned long long)(UInt128)x.num, b = (unsigned long long)((UInt128)x.num >> 64);
        const unsigned long long c = (unsigned long long)(UInt128)x.den, d = (unsigned long long)((UInt128)x.den >> 64);
        unsigned long long h = a * 0x9E3779B97F4A7C15ULL ^ (b + 0x632BE59BD9B4E019ULL);
        h = (h ^ (h >> 29)) * 0xBF58476D1CE4E5B9ULL ^ c;
        h = (h ^ (h >> 31)) * 0x94D049BB133111EBULL ^ (d + 0x7F4A7C15ULL);
        return (size_t)(h ^ (h >> 32));
    }
};
}  // namespace std

// p/q（q ≠ 0）にいちばん近い、分母が max_den 以下の分数（同じ近さなら分母の小さい方、分母も同じなら値の小さい方）  O(log)
inline Fraction best_rational_approx(long long p, long long q, long long max_den) {
    if (q == 0 || max_den < 1) math_internal::fail("best_rational_approx", "q ≠ 0、max_den ≥ 1 で使う");
    Int128 P = p, Q = q;
    if (Q < 0) P = -P, Q = -Q;
    {
        const Int128 g = math_internal::frac_gcd(P, Q);
        P /= g, Q /= g;
        if (Q <= max_den) return Fraction(P, Q);
    }
    // 連分数の近似分数 h/k を分母が max_den を超える手前まで
    Int128 h0 = 0, h1 = 1, k0 = 1, k1 = 0, a = P, b = Q;
    while (b != 0) {
        const Int128 t = math_internal::floor_div128(a, b);
        const Int128 k2 = t * k1 + k0;
        if (k2 > max_den) break;
        const Int128 h2 = t * h1 + h0;
        h0 = h1, h1 = h2, k0 = k1, k1 = k2;
        const Int128 r = a - t * b;
        a = b, b = r;
    }
    // 候補：最後の近似分数 h1/k1 と、中間近似分数 (h0 + m h1)/(k0 + m k1)（m は分母が収まる最大）
    const Int128 m = (max_den - k0) / k1;
    const Int128 hs = h0 + m * h1, ks = k0 + m * k1;
    const Int128 d1 = math_internal::frac_abs(P * k1 - h1 * Q), d2 = math_internal::frac_abs(P * ks - hs * Q);  // 距離 × Q × 分母
    const Int128 l = d1 * ks, r = d2 * k1;  // d1/k1 と d2/ks を比べる
    if (l < r || (l == r && k1 <= ks)) return Fraction(h1, k1);
    return Fraction(hs, ks);
}
namespace math_internal {
inline void frac_bounded(long long p, long long q, long long max_den, const char* func, Int128& lp, Int128& lq, Int128& rp, Int128& rq) {
    if (q == 0 || max_den < 1) fail(func, "q ≠ 0、max_den ≥ 1 で使う");
    Int128 P = p, Q = q;
    if (Q < 0) P = -P, Q = -Q;
    const Int128 g = frac_gcd(P, Q);
    P /= g, Q /= g;
    if (Q <= max_den) {  // x そのもの
        lp = rp = P, lq = rq = Q;
        return;
    }
    // L < x < R を隣り合う分数（Stern–Brocot 木）のまま、分母 max_den 以下で x に近づける
    lp = floor_div128(P, Q), lq = 1, rp = lp + 1, rq = 1;
    while (true) {
        Int128 k1 = (P * lq - lp * Q) / (rp * Q - P * rq);  // L + k R ≤ x となる最大の k
        k1 = std::min(k1, (max_den - lq) / rq);
        lp += k1 * rp, lq += k1 * rq;
        Int128 k2 = (rp * Q - P * rq) / (P * lq - lp * Q);  // R + k L ≥ x となる最大の k
        k2 = std::min(k2, (max_den - rq) / lq);
        rp += k2 * lp, rq += k2 * lq;
        if (k1 == 0 && k2 == 0) return;
    }
}
}  // namespace math_internal
// p/q 以下で最大の、分母が max_den 以下の分数 / p/q 以上で最小のもの（p/q の分母が max_den 以下なら p/q）  O(log)
inline Fraction rational_floor(long long p, long long q, long long max_den) {
    Int128 lp, lq, rp, rq;
    math_internal::frac_bounded(p, q, max_den, "rational_floor", lp, lq, rp, rq);
    return Fraction(lp, lq);
}
inline Fraction rational_ceil(long long p, long long q, long long max_den) {
    Int128 lp, lq, rp, rq;
    math_internal::frac_bounded(p, q, max_den, "rational_ceil", lp, lq, rp, rq);
    return Fraction(rp, rq);
}
// x < p/q < y となる分数のうち、分母が最小のもの（分母が同じなら |分子| が最小。開区間。Stern–Brocot 木でいちばん上）  O(log)
inline Fraction simplest_between(const Fraction& x, const Fraction& y) {
    if (!(x < y)) math_internal::fail("simplest_between", "x < y で使う");
    if (x.sign() < 0 && y.sign() > 0) return Fraction(0);
    if (y.sign() <= 0) return -simplest_between(-y, -x);
    // 0 ≤ a/b < c/d の開区間で（分母最小 = 分子も最小）
    auto rec = [](auto&& self, Int128 a, Int128 b, Int128 c, Int128 d) -> std::pair<Int128, Int128> {
        const Int128 k = a / b;
        if ((k + 1) * d < c) return {k + 1, 1};  // 整数 k + 1 が入る
        a -= k * b, c -= k * d;                  // 0 ≤ a/b < c/d ≤ 1 に
        if (a == 0) {                            // (0, c/d) で分母最小：1/q（q = ⌊d/c⌋ + 1）
            const Int128 q = d / c + 1;
            return {k * q + 1, q};
        }
        const std::pair<Int128, Int128> r = self(self, d, c, b, a);  // 逆数をとる：d/c < P/Q < b/a
        return {k * r.first + r.second, r.first};
    };
    const std::pair<Int128, Int128> r = rec(rec, x.num, x.den, y.num, y.den);
    return Fraction(r.first, r.second);
}
// 連分数展開 p/q = a0 + 1/(a1 + 1/(a2 + …))（q ≠ 0。a0 は負もある。最後の項は 2 以上（p/q が整数なら [a0] だけ））  O(log)
inline std::vector<long long> continued_fraction(long long p, long long q) {
    if (q == 0) math_internal::fail("continued_fraction", "q ≠ 0 で使う");
    Int128 a = p, b = q;
    if (b < 0) a = -a, b = -b;
    std::vector<long long> cf;
    while (b != 0) {
        const Int128 t = math_internal::floor_div128(a, b);
        if (t > (Int128)LLONG_MAX) math_internal::fail("continued_fraction", "項が long long に収まらない（p/q = 2^63）");
        cf.push_back((long long)t);
        const Int128 r = a - t * b;
        a = b, b = r;
    }
    return cf;
}
// 連分数 [a0; a1, a2, …] の値（a1 以降は 1 以上。溢れたら止める）  O(長さ)
template <class V>
Fraction from_continued_fraction(const V& cfv) {
    std::vector<long long> cf;
    for (const auto& x : cfv) cf.push_back(math_internal::to_ll_checked(x, "from_continued_fraction"));
    if (cf.empty()) math_internal::fail("from_continued_fraction", "空の連分数");
    Fraction x(cf.back());
    for (size_t i = cf.size() - 1; i > 0; --i) {
        if (cf[i] < 1) math_internal::fail("from_continued_fraction", "a1 以降は 1 以上で使う");
        x = Fraction(cf[i - 1]) + x.reciprocal();
    }
    return x;
}
inline Fraction from_continued_fraction(std::initializer_list<long long> cf) { return from_continued_fraction<std::initializer_list<long long>>(cf); }
#endif  // GG_FRACTION
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace frac_test {

#define FRAC_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[fraction] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
static Int128 g128(Int128 a, Int128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        const Int128 t = a % b;
        a = b, b = t;
    }
    return a;
}
static long long nz(std::mt19937_64& rng, long long lo, long long hi) {  // 0 でない乱数
    long long v;
    do v = rnd(rng, lo, hi);
    while (v == 0);
    return v;
}
template <class F>
static bool same_value(const F& x, Int128 n, Int128 d) {  // x == n/d（d ≠ 0）か（既約・分母 > 0 も確かめる）
    if (x.den <= 0 || g128(x.num, x.den) != 1) return false;
    return (Int128)x.num * d == n * (Int128)x.den;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 50000 * TEST_SCALE; ++it, ++cases) {
        const long long R = it % 3 == 0 ? 20 : (it % 3 == 1 ? 1000000000 : LMAX);
        const long long a = rnd(rng, -R, R), b = nz(rng, -R, R);
        const Fraction x(a, b);
        FRAC_CHECK(same_value(x, a, b));
        // 四則演算（値が long long に収まる大きさで）
        const long long S = 1000000000;
        const long long p = rnd(rng, -S, S), q = nz(rng, -S, S), r = rnd(rng, -S, S), s = nz(rng, -S, S);
        const Fraction u(p, q), v(r, s);
        FRAC_CHECK(same_value(u + v, (Int128)p * s + (Int128)r * q, (Int128)q * s));
        FRAC_CHECK(same_value(u - v, (Int128)p * s - (Int128)r * q, (Int128)q * s));
        FRAC_CHECK(same_value(u * v, (Int128)p * r, (Int128)q * s));
        if (r != 0) FRAC_CHECK(same_value(u / v, (Int128)p * s, (Int128)q * r));
        FRAC_CHECK(same_value(-u, -p, q) && same_value(u * 3 + 1, (Int128)3 * p + q, q) && same_value(2 - u, (Int128)2 * q - p, q));
        // 比較（大きい値でも）
        const long long c = rnd(rng, -R, R), d = nz(rng, -R, R);
        const Fraction y(c, d);
        const Int128 lhs = (Int128)a * d * (b < 0 ? -1 : 1) * (d < 0 ? -1 : 1), rhs = (Int128)c * b * (b < 0 ? -1 : 1) * (d < 0 ? -1 : 1);
        const int e = lhs < rhs ? -1 : (lhs > rhs ? 1 : 0);
        FRAC_CHECK(compare(x, y) == e && (x < y) == (e < 0) && (x == y) == (e == 0) && (x >= y) == (e >= 0));
        // 切り捨て・切り上げ・四捨五入
        FRAC_CHECK((Int128)x.floor() == math_internal::floor_div128(a, b) && (Int128)x.ceil() == math_internal::ceil_div128(a, b));
        const Int128 fl = math_internal::floor_div128(a, b), two_rem = 2 * ((Int128)x.num - fl * x.den);  // 2 × (x − ⌊x⌋) × den
        Int128 rd = two_rem < x.den ? fl : (two_rem > x.den ? fl + 1 : (x.num >= 0 ? fl + 1 : fl));
        FRAC_CHECK((Int128)x.round() == rd);
        FRAC_CHECK(x.sign() == (a == 0 ? 0 : ((a < 0) != (b < 0) ? -1 : 1)) && x.abs() >= Fraction(0) && (x.is_integer() == (a % b == 0)));
        if (a != 0) FRAC_CHECK(x * x.reciprocal() == Fraction(1));
        FRAC_CHECK(std::fabs((double)x.to_real() - (double)a / (double)b) <= 1e-9 * (1 + std::fabs((double)a / (double)b)));
    }
}

static void test_text(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long a = rnd(rng, -LMAX, LMAX), b = nz(rng, 1, it % 2 ? 1000 : LMAX);
        const Fraction x(a, b);
        // 出して読み戻す
        std::ostringstream os;
        os << x;
        std::istringstream is(os.str());
        Fraction y;
        is >> y;
        FRAC_CHECK(x == y && (x.den != 1 || os.str().find('/') == std::string::npos));
        // 小数 k 桁（切り捨て）：愚直な割り算
        const int k = (int)rnd(rng, 0, 25);
        const Int128 A = a < 0 ? -(Int128)a : a;
        Int128 rem = A % b;
        std::string digits;
        for (int i = 0; i < k; ++i) rem *= 10, digits += (char)('0' + (int)(rem / b)), rem %= b;
        const std::string body = int128_to_string(A / b) + (k ? "." + digits : "");
        const bool zero = A / b == 0 && digits.find_first_not_of('0') == std::string::npos;
        const std::string e = (a < 0 && !zero ? "-" : "") + body;
        FRAC_CHECK(Fraction(a, b).to_decimal(k) == e);
        // 小数を読む：to_decimal の文字列はその値（切り捨て）を表す
        Int128 p10 = 1;
        for (int i = 0; i < k; ++i) p10 *= 10;
        if (k > 0 && k <= 18 && (A / b + 1) * p10 < (Int128)LMAX) {  // 読んだ値が long long に収まるとき
            const Fraction t = Fraction::parse(e);
            FRAC_CHECK((Int128)t.num * p10 == (parse_int128(int128_to_string(A / b) + digits) * (a < 0 && !zero ? -1 : 1)) * (Int128)t.den);
        }
    }
    FRAC_CHECK(Fraction::parse("-1.25") == Fraction(-5, 4) && Fraction::parse("3/-6") == Fraction(-1, 2) && Fraction::parse(".5") == Fraction(1, 2) && Fraction::parse("7") == Fraction(7));
    FRAC_CHECK(Fraction(-1, 3).to_decimal(5) == "-0.33333" && Fraction(-1, 3000000).to_decimal(5) == "0.00000" && Fraction(7, 2).to_decimal(0) == "3");
    // ハッシュ
    std::unordered_map<Fraction, int> mp;
    for (int i = 1; i <= 100; ++i) mp[Fraction(i, 2 * i)]++, mp[Fraction(-i, 3 * i)]++;
    FRAC_CHECK(mp.size() == 2 && mp[Fraction(1, 2)] == 100);
}

static void test_128(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 15000 * TEST_SCALE; ++it, ++cases) {
        const long long R = 1LL << 60;
        const long long a = rnd(rng, -R, R), b = nz(rng, -R, R), c = rnd(rng, -R, R), d = nz(rng, -R, R);
        const Fraction128 x(a, b), y(c, d);
        FRAC_CHECK(same_value(x, a, b) && x + y == Fraction128((Int128)a * d + (Int128)c * b, (Int128)b * d) && x * y == Fraction128((Int128)a * c, (Int128)b * d));
        FRAC_CHECK(x - y == Fraction128((Int128)a * d - (Int128)c * b, (Int128)b * d) && (c == 0 || x / y == Fraction128((Int128)a * d, (Int128)b * c)));
        const Int128 lhs = (Int128)a * d * (b < 0 ? -1 : 1) * (d < 0 ? -1 : 1), rhs = (Int128)c * b * (b < 0 ? -1 : 1) * (d < 0 ? -1 : 1);
        const int e = lhs < rhs ? -1 : (lhs > rhs ? 1 : 0);
        FRAC_CHECK(compare(x, y) == e && compare(y, x) == -e);
        // 大きい分子・分母：比較を long double と（十分離れているとき）
        const Int128 P60 = (Int128)1 << 60;
        const Fraction128 big1 = Fraction128((Int128)rnd(rng, -LMAX, LMAX) * P60 + rnd_ll(rng), (Int128)rnd(rng, 1, LMAX) * P60 + 1);
        const Fraction128 big2 = Fraction128((Int128)rnd(rng, -LMAX, LMAX) * P60 + rnd_ll(rng), (Int128)rnd(rng, 1, LMAX) * P60 + 3);
        const long double r1 = big1.to_real(), r2 = big2.to_real();
        if (std::fabs((double)(r1 - r2)) > 1e-6 * (std::fabs((double)r1) + std::fabs((double)r2) + 1)) FRAC_CHECK((big1 < big2) == (r1 < r2));
        FRAC_CHECK(big1.to_decimal(30).size() > 30 && Fraction128::parse(big1.to_decimal(0)) == Fraction128(big1.num / big1.den));
        FRAC_CHECK(compare(big1, big1) == 0 && compare(big1 + Fraction128(1, 3), big1) == 1 && compare(big1, big1 + Fraction128(1, 3)) == -1);
    }
}

static void test_edge128(int& cases) {  // 途中の値が −2^127 になる足し算、分母が 2^125 を超える小数表示
    const Int128 P126 = (Int128)1 << 126;
    const Fraction128 x(-P126 - 1, 2), y(-P126 + 1, 2), z(P126 - 1, 2);
    FRAC_CHECK(x + y == Fraction128(-P126) && x - z == Fraction128(-P126));
    const Int128 d = (Int128)3 << 125;  // (d − 1)/d = 0.999…
    FRAC_CHECK(Fraction128(d - 1, d).to_decimal(3) == "0.999" && Fraction128(-(d - 1), d).to_decimal(5) == "-0.99999" && Fraction128(1, d).to_decimal(3) == "0.000");
    cases += 3;
}

static void test_approx(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 15000 * TEST_SCALE; ++it, ++cases) {
        const long long q = nz(rng, -1000000, 1000000), p = it % 2 ? rnd(rng, -1000000, 1000000) : rnd(rng, -50, 50);
        const long long M = rnd(rng, 1, 60);
        // 愚直：分母 k = 1〜M、分子は p k / q の前後
        Int128 bh = 0, bk = 0;
        for (long long k = 1; k <= M; ++k) {
            const Int128 f = math_internal::floor_div128((Int128)p * k, q);
            for (Int128 h = f; h <= f + 1; ++h) {
                // |p/q − h/k| = |p k − h q| / (|q| k)。前の候補 bh/bk と比べる
                const Int128 dn = (Int128)p * k - h * q, dist = dn < 0 ? -dn : dn;
                if (bk == 0) {
                    bh = h, bk = k;
                    continue;
                }
                const Int128 bdn = (Int128)p * bk - bh * q, bdist = bdn < 0 ? -bdn : bdn;
                if (dist * bk < bdist * k) bh = h, bk = k;
            }
        }
        const Fraction got = best_rational_approx(p, q, M);
        FRAC_CHECK(got == Fraction(bh, bk));
    }
    FRAC_CHECK(best_rational_approx(314159265358979LL, 100000000000000LL, 1000) == Fraction(355, 113) && best_rational_approx(3, 4, 2) == Fraction(1));
    // 分母が max_den 以下で p/q 以下の最大・以上の最小
    for (int it = 0; it < 15000 * TEST_SCALE; ++it, ++cases) {
        const long long q = nz(rng, -100000, 100000), p = it % 2 ? rnd(rng, -1000000, 1000000) : rnd(rng, -50, 50);
        const long long M = rnd(rng, 1, 60);
        const Fraction x(p, q);
        Fraction lo, hi;
        bool has = false;
        for (long long k = 1; k <= M; ++k) {
            const Fraction f(math_internal::floor_div128((Int128)p * k, q), k), c(math_internal::ceil_div128((Int128)p * k, q), k);
            if (!has || lo < f) lo = f;
            if (!has || c < hi) hi = c;
            has = true;
        }
        FRAC_CHECK(rational_floor(p, q, M) == lo && rational_ceil(p, q, M) == hi && lo <= x && x <= hi);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい値：隣り合う（右 − 左 = 1/(分母の積)）こと、分母の和が max_den を超えること
        const long long q = nz(rng, -LLONG_MAX, LLONG_MAX), p = rnd(rng, -LLONG_MAX, LLONG_MAX), M = rnd(rng, 1, it % 2 ? 1000000 : LLONG_MAX);  // Fraction の分子は −2^63 を使わない
        const Fraction l = rational_floor(p, q, M), r = rational_ceil(p, q, M);
        const Fraction128 x(p, q);
        FRAC_CHECK(l.den <= M && r.den <= M && Fraction128(l.num, l.den) <= x && x <= Fraction128(r.num, r.den));
        if (!(l == r)) FRAC_CHECK((Int128)r.num * l.den - (Int128)l.num * r.den == 1 && (Int128)l.den + r.den > M);
    }
    // いちばん簡単な分数（開区間）
    for (int it = 0; it < 15000 * TEST_SCALE; ++it, ++cases) {
        Fraction x(rnd(rng, -30, 30), rnd(rng, 1, 30)), y(rnd(rng, -30, 30), rnd(rng, 1, 30));
        if (x == y) continue;
        if (y < x) std::swap(x, y);
        Fraction best;
        bool found = false;
        for (long long k = 1; k <= 61 && !found; ++k) {  // 分母の小さい順、分子は絶対値の小さい順
            for (long long a = 0; a <= 2000 && !found; ++a)
                for (long long s : {1, -1}) {
                    const Fraction f(s * a, k);
                    if (x < f && f < y) {
                        best = f, found = true;
                        break;
                    }
                }
        }
        FRAC_CHECK(found && simplest_between(x, y) == best);
    }
    FRAC_CHECK(simplest_between(Fraction(1, 3), Fraction(1, 2)) == Fraction(2, 5) && simplest_between(Fraction(-7, 2), Fraction(-3)) == Fraction(-10, 3));
    FRAC_CHECK(simplest_between(Fraction(999999999999999999LL, 1000000000000000000LL), Fraction(1)).den == 1000000000000000001LL);
    // 連分数
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long q = nz(rng, -LLONG_MAX, LLONG_MAX), p = rnd(rng, -LLONG_MAX, LLONG_MAX);
        const std::vector<long long> cf = continued_fraction(p, q);
        for (size_t i = 1; i < cf.size(); ++i) FRAC_CHECK(cf[i] >= 1 && (i + 1 < cf.size() || cf[i] >= 2));
        FRAC_CHECK(from_continued_fraction(cf) == Fraction(p, q));
    }
    FRAC_CHECK(continued_fraction(355, 113) == std::vector<long long>({3, 7, 16}) && continued_fraction(-7, 3) == std::vector<long long>({-3, 1, 2}) && from_continued_fraction({1, 2, 2, 2}) == Fraction(17, 12));
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_basic(rng, cases);
    test_text(rng, cases);
    test_128(rng, cases);
    test_edge128(cases);
    test_approx(rng, cases);
    std::cout << "[fraction] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[fraction] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 1000000;
    std::vector<Fraction> v(n);
    for (auto& x : v) x = Fraction(rnd(rng, 0, 1000000000), rnd(rng, 1, 1000000000));
    auto t0 = std::chrono::steady_clock::now();
    std::sort(v.begin(), v.end());
    const long long t_sort = ms(t0);
    t0 = std::chrono::steady_clock::now();
    Fraction s = 0;
    for (int i = 0; i < n; ++i) s = Fraction(v[(size_t)i].num % 1000 + 1, v[(size_t)i].den % 1000 + 1) * Fraction(i % 7 + 1, i % 5 + 1) + Fraction(1, i % 9 + 1);
    const long long t_ops = ms(t0);
    std::vector<Fraction128> w(n / 5);
    for (auto& x : w) x = Fraction128(((Int128)rnd(rng, 0, LMAX) << 40) + 1, ((Int128)rnd(rng, 1, LMAX) << 40) + 1);
    t0 = std::chrono::steady_clock::now();
    std::sort(w.begin(), w.end());
    const long long t_sort128 = ms(t0);
    std::cout << "[fraction] benchmark: Fraction 10^6 個の sort " << t_sort << " ms, 四則演算 × 10^6 " << t_ops << " ms, Fraction128（分子・分母 ≈ 2^100）2×10^5 個の sort "
              << t_sort128 << " ms  (check " << (long long)(s.num % 7 + w[5].num % 7) << ")\n";
}

}  // namespace frac_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace frac_demo {

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
#define SHOW(expr, comment) frac_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    frac_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("分数（fraction） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Fraction x(n, d) / x.num / x.den / + − × ÷ / == < など / compare(x, y)\n"
              << "  x.floor() / x.ceil() / x.round() / x.sign() / x.abs() / x.reciprocal() / x.is_integer() / x.to_real() / x.to_decimal(k)\n"
              << "  cout << x / cin >> x / Fraction::parse(s) / unordered_map<Fraction, V>\n"
              << "  best_rational_approx(p, q, max_den)      分母が小さい近い分数\n"
              << "  Fraction128                              分子・分母が Int128\n";

    api("【1】 Fraction x(n, d) / + − × ÷ / 比較", "Fraction",
        "いつも既約で分母 > 0。途中は 128 ビットで計算し、結果が long long に収まらなければ止める", "O(log)");
    RUN(Fraction a(1, 3), b(-2, 4));
    SHOW(b, "約分され、符号は分子に");
    SHOW(a + b, "1/3 − 1/2");
    SHOW(a * b, "掛け算");
    SHOW(a / b, "割り算");
    SHOW(a * 3, "整数と混ぜても");
    SHOW(a < b, "誤差なしの比較");
    SHOW(compare(Fraction(2, 6), a), "等しいので 0");
    SHOW(a.num, "分子");
    SHOW(a.den, "分母");
    {
        // ABC308 C：A/(A+B) の大きい順、同じなら番号の小さい順（stable_sort + 誤差なしの比較）
        auto solve = [](const std::vector<std::pair<long long, long long>>& ab) {
            std::vector<int> id(ab.size());
            std::iota(id.begin(), id.end(), 0);
            std::stable_sort(id.begin(), id.end(), [&](int i, int j) {
                return Fraction(ab[(size_t)i].first, ab[(size_t)i].first + ab[(size_t)i].second) > Fraction(ab[(size_t)j].first, ab[(size_t)j].first + ab[(size_t)j].second);
            });
            std::string s;
            for (int i : id) s += (s.empty() ? "" : " ") + std::to_string(i + 1);
            return s;
        };
        note("ABC308 C（Standings）: 成功率 A/(A+B) を Fraction で比べて stable_sort");
        std::cout << "    例 1 → " << solve({{1, 3}, {3, 1}, {2, 2}}) << "、例 2 → " << solve({{1, 3}, {2, 6}})
                  << "、例 3 → " << solve({{999999999, 1000000000}, {333333333, 999999999}, {1000000000, 999999997}, {999999998, 1000000000}}) << "\n";
    }

    api("【2】 x.floor() / x.ceil() / x.round() / x.sign() / x.abs() / x.reciprocal() / x.to_real() / x.to_decimal(k)", "long long / Fraction / long double / string",
        "切り捨て・切り上げ・四捨五入（.5 は 0 から遠い方）/ 符号 / 絶対値 / 逆数 / 実数 / 小数第 k 位まで（切り捨て）", "O(1) / O(k)");
    RUN(Fraction c(-7, 2));
    SHOW(c.floor(), "⌊−3.5⌋");
    SHOW(c.ceil(), "⌈−3.5⌉");
    SHOW(c.round(), "−3.5 の四捨五入（0 から遠い方）");
    SHOW(c.abs(), "絶対値");
    SHOW(c.reciprocal(), "逆数");
    SHOW(Fraction(1, 7).to_decimal(20), "1/7 を小数第 20 位まで（誤差なし）");
    note("出力の桁を決めて実数で出すなら cout << fixed << setprecision(15) << (double)x.to_real()");

    api("【3】 cout << x / cin >> x / Fraction::parse(s) / unordered_map<Fraction, V>", "", "\"p/q\"（整数なら \"p\"）で出す / \"p/q\"・\"p\"・\"-1.25\" を読む / ハッシュ", "O(桁数)");
    SHOW(Fraction::parse("-1.25"), "小数を読む");
    SHOW(Fraction::parse("6/-4"), "分数を読む");
    RUN(std::unordered_map<Fraction, int> cnt; cnt[Fraction(1, 2)]++; cnt[Fraction(2, 4)]++);
    SHOW(cnt[Fraction(1, 2)], "同じ値は同じキー（傾きを数えるときなど）");

    api("【4】 best_rational_approx(p, q, max_den)", "Fraction", "p/q にいちばん近い、分母が max_den 以下の分数（同じ近さなら分母の小さい方）", "O(log)");
    SHOW(best_rational_approx(314159265358979LL, 100000000000000LL, 1000), "π ≈ 355/113");
    SHOW(best_rational_approx(314159265358979LL, 100000000000000LL, 100), "分母 100 以下なら 311/99");
    SHOW(best_rational_approx(1, 3, 2), "1/3 に近い分母 2 以下：1/2 と 0 は 1/6 と 1/3 離れている");

    api("【5】 Fraction128", "Fraction128", "分子・分母が Int128（10^38 まで）。比較は掛け算をしない方法で誤差も溢れもない", "O(log)");
    RUN(Fraction128 h(1000000000000000000LL, 3));
    SHOW(h * h, "10^36 / 9（long long には入らない）");
    SHOW(h * h > Fraction128((Int128)111111111111111111LL * 1000000000000000000LL), "比較");

    api("【6】 rational_floor(p, q, N) / rational_ceil(p, q, N) / simplest_between(x, y) / continued_fraction(p, q) / from_continued_fraction(cf)",
        "Fraction / vector<long long>", "分母 N 以下で p/q 以下の最大・以上の最小 / x < r < y で分母最小の r / 連分数展開とその逆", "O(log)");
    SHOW(rational_floor(314159265358979LL, 100000000000000LL, 100), "π 以下で分母 100 以下の最大");
    SHOW(rational_ceil(314159265358979LL, 100000000000000LL, 100), "π 以上で分母 100 以下の最小");
    SHOW(simplest_between(Fraction(1, 3), Fraction(1, 2)), "1/3 と 1/2 のあいだで分母最小");
    SHOW(simplest_between(Fraction(-7, 2), Fraction(-3)), "負でも");
    SHOW(continued_fraction(355, 113), "355/113 = 3 + 1/(7 + 1/16)");
    SHOW(from_continued_fraction({1, 2, 2, 2}), "√2 の近似 17/12");
    note("区間 (x, y) に入る分数の分母の最小値は simplest_between(x, y).den（それより小さい分母の分数は入らない）");
}

}  // namespace frac_demo

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
    frac_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    frac_test::benchmark();    // 速度計測
    frac_demo::demo();         // 使い方デモ
    return 0;
}
