// =============================================================================
//  XOR の基底（XorBasis：XOR で作れる値の最大・k 番目・作れるか・どの要素で作るか、区間の要素の基底、mod 2 の連立一次方程式）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 行列・線形代数
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  「作れる値」= いくつかの要素（0 個でもよい）の XOR。値は unsigned long long（64 ビット）
//
//  機能一覧                                                                    計算量
//    XorBasis B;                      空の基底
//      B.insert(x)                    x を加える（作れる値が増えたら true）                         O(64)
//      B.contains(x)                  x を作れるか
//      B.max_xor(x = 0) / B.min_xor(x = 0)   x XOR (作れる値) の最大 / 最小
//      B.kth(k)                       作れる値を小さい順に並べた k 番目（0 から。0 番目は 0。k < 2^rank）
//      B.count_le(x)                  x 以下の作れる値の個数（Int128）
//      B.rank() / B.count()           基底の大きさ / 作れる値の個数 2^rank（Int128）
//      B.zero_by_nonempty()           1 個以上の要素の XOR で 0 を作れるか（insert で増えなかった要素がある）
//      B.represent(x)                 x を作る要素の番号（insert した順に 0, 1, …）。Maybe
//      B.merge(C) / B.basis()         C の分も加える / 基底（簡約済み。大きい順）
//    PrefixXorBasis P;                a_0, a_1, … を P.push_back(a) で順に足し、今の末尾までの [l, size()) について
//      P.contains(l, x) / P.max_xor(l, x = 0) / P.rank(l) / P.basis(l)    （クエリは右端の順に並べて答える） O(64)
//    Gf2System S(m);                  mod 2 の連立一次方程式（変数 x_0〜x_{m−1}）
//      S.add(vars, rhs) / S.add_row(coef, rhs)   vars の XOR = rhs / Σ coef_j x_j = rhs を加える
//      S.rank() / S.solvable() / S.free_count()  階数 / 解があるか / 自由な変数の数（解の数は 2^free_count）
//      S.solve() / S.kernel()         解の 1 つ（Maybe）/ 右辺 0 の解の基底                         O(式 × rank × m / 64)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    xorb_test::random_test()  … 部分集合を全部試す愚直解との突き合わせ
//    xorb_test::benchmark()    … 速度計測
//    xorb_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
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
#ifndef GG_XOR_BASIS
#define GG_XOR_BASIS
// ---- XOR の基底（数学第二世代）  XorBasis（64 ビットの値の XOR で作れる値）・PrefixXorBasis（区間の要素）・Gf2System（mod 2 の連立方程式） ----
// 「作れる値」= いくつかの要素（0 個でもよい）の XOR。0 個なら 0
struct XorBasis {
    using u64 = unsigned long long;
    XorBasis() {
        for (int i = 0; i < 64; ++i) b_[i] = 0, comp_[i] = 0;
    }
    // x を加える。新しく作れる値が増えたら true（既に作れるなら false）  O(64)
    bool insert(u64 x) {
        const int id = inserted_++;
        u64 c = 0;  // 引いた基底がどの slot の XOR か
        // 簡約した基底は他の基底の先頭ビットを持たないので、x の「基底の先頭ビット」の位置の基底を全部 XOR すればよい
        for (u64 t = x & have_; t; t &= t - 1) {
            const int i = __builtin_ctzll(t);
            x ^= b_[i], c ^= comp_[i];
        }
        if (x == 0) return false;
        const int h = 63 - __builtin_clzll(x), s = rank_++;
        slot_id_[s] = id;
        c ^= 1ULL << s;
        for (u64 t = have_ >> h >> 1; t; t &= t - 1) {  // 上の基底から x の先頭ビットを消す（簡約を保つ）
            const int i = h + 1 + __builtin_ctzll(t);
            if (b_[i] >> h & 1) b_[i] ^= x, comp_[i] ^= c;
        }
        b_[h] = x, comp_[h] = c, have_ |= 1ULL << h;
        return true;
    }
    bool contains(u64 x) const { return reduce(x) == 0; }  // x が作れるか
    u64 reduce(u64 x) const {                              // min_{s は作れる値} (x XOR s)
        for (u64 t = x & have_; t; t &= t - 1) x ^= b_[__builtin_ctzll(t)];
        return x;
    }
    u64 min_xor(u64 x = 0) const { return reduce(x); }  // x XOR (作れる値) の最小（x = 0 なら 0）
    u64 max_xor(u64 x = 0) const {                      // x XOR (作れる値) の最大（先頭ビットが 0 のところの基底を全部 XOR）
        for (u64 t = ~x & have_; t; t &= t - 1) x ^= b_[__builtin_ctzll(t)];
        return x;
    }
    // 作れる値を小さい順に並べた k 番目（0 から。0 番目は 0。k < 2^rank）
    u64 kth(u64 k) const {
        if (rank_ < 64 && (k >> rank_) != 0) math_internal::fail("XorBasis::kth", "k < 2^rank で使う");
        u64 x = 0;
        int j = 0;
        for (u64 t = have_; t; t &= t - 1, ++j)
            if (k >> j & 1) x ^= b_[__builtin_ctzll(t)];
        return x;
    }
    Int128 count_le(u64 x) const {  // x 以下の作れる値の個数（kth が増えていくことを使い、k を上のビットから決める）  O(64)
        u64 k = 0, cur = 0;
        int j = rank_;
        for (int i = 63; i >= 0; --i) {
            if (!b_[i]) continue;
            --j;
            if ((cur ^ b_[i]) <= x) cur ^= b_[i], k |= 1ULL << j;  // k のビット j を立てて、下のビットは 0 の値
        }
        return (Int128)k + 1;
    }
    int rank() const { return rank_; }
    Int128 count() const { return (Int128)1 << rank_; }  // 作れる値の個数 2^rank
    int inserted() const { return inserted_; }          // insert を呼んだ回数
    // 0 個でない要素の XOR で 0 を作れるか（insert したのに増えなかった要素があるか）
    bool zero_by_nonempty() const { return inserted_ > rank_; }
    void merge(const XorBasis& o) {  // o の作れる値も作れるようにする（represent の番号は o の基底の分だけ増える）
        for (int i = 0; i < 64; ++i)
            if (o.b_[i]) insert(o.b_[i]);
    }
    std::vector<u64> basis() const {  // 基底（簡約済み。先頭ビットの大きい順）
        std::vector<u64> r;
        for (int i = 63; i >= 0; --i)
            if (b_[i]) r.push_back(b_[i]);
        return r;
    }
    // x を作る要素の番号（insert を呼んだ順に 0, 1, …）。作れなければ ok = false
    Maybe<std::vector<int>> represent(u64 x) const {
        u64 c = 0;
        for (int i = 63; i >= 0; --i)
            if (x >> i & 1) {
                if (!b_[i]) return {false, {}};
                x ^= b_[i], c ^= comp_[i];
            }
        std::vector<int> r;
        for (int s = 0; s < rank_; ++s)
            if (c >> s & 1) r.push_back(slot_id_[s]);
        std::sort(r.begin(), r.end());
        return {true, r};
    }

  private:
    u64 b_[64];        // b_[i]：先頭ビットが i の基底（なければ 0）
    u64 comp_[64];     // b_[i] がどの slot の XOR か
    int slot_id_[64];  // slot s の要素が何番目に insert されたか
    u64 have_ = 0;     // 基底がある先頭ビットの集合
    int rank_ = 0, inserted_ = 0;
};

// ---- 区間の要素の XOR：a_0, a_1, … を順に push_back して、今の末尾までの [l, size()) について答える ----
// クエリは右端の小さい順に並べ替えて答える（各ビットの基底に「いちばん右の番号」を残す）  push_back・各クエリ O(64)
struct PrefixXorBasis {
    using u64 = unsigned long long;
    PrefixXorBasis() {
        for (int i = 0; i < 64; ++i) b_[i] = 0, pos_[i] = -1;
    }
    void push_back(u64 x) {  // a_{size()} = x を足す
        int p = n_++;
        for (int i = 63; i >= 0 && x; --i) {
            if (!(x >> i & 1)) continue;
            if (!b_[i]) {
                b_[i] = x, pos_[i] = p;
                return;
            }
            if (pos_[i] < p) std::swap(b_[i], x), std::swap(pos_[i], p);  // 新しい（右の）方を基底に残す
            x ^= b_[i];
        }
    }
    int size() const { return n_; }
    bool contains(int l, u64 x) const {  // [l, size()) の要素（0 個でもよい）の XOR で x を作れるか
        check(l, "PrefixXorBasis::contains");
        for (int i = 63; i >= 0; --i)
            if ((x >> i & 1) && b_[i] && pos_[i] >= l) x ^= b_[i];
        return x == 0;
    }
    u64 max_xor(int l, u64 x = 0) const {  // x XOR ([l, size()) で作れる値) の最大
        check(l, "PrefixXorBasis::max_xor");
        for (int i = 63; i >= 0; --i)
            if (b_[i] && pos_[i] >= l && !(x >> i & 1)) x ^= b_[i];
        return x;
    }
    int rank(int l) const {  // [l, size()) の要素の基底の大きさ
        check(l, "PrefixXorBasis::rank");
        int r = 0;
        for (int i = 0; i < 64; ++i) r += b_[i] && pos_[i] >= l;
        return r;
    }
    XorBasis basis(int l) const {  // [l, size()) の要素の基底（XorBasis として。kth なども使える）
        check(l, "PrefixXorBasis::basis");
        XorBasis r;
        for (int i = 63; i >= 0; --i)
            if (b_[i] && pos_[i] >= l) r.insert(b_[i]);
        return r;
    }

  private:
    u64 b_[64];
    int pos_[64];
    int n_ = 0;
    void check(int l, const char* f) const {
        if (l < 0 || l > n_) math_internal::fail(f, "0 ≤ l ≤ size() で使う");
    }
};

// ---- mod 2 の連立一次方程式（変数 x_0〜x_{m−1} ∈ {0, 1}。ビットで持つ）----
// 式を add したあと rank / solve / kernel。消去は O(式の数 × rank × m / 64)
struct Gf2System {
    using u64 = unsigned long long;
    explicit Gf2System(int m) : m_(m), w_((m + 64) / 64) {
        if (m < 0) math_internal::fail("Gf2System", "m ≥ 0 で使う");
    }
    void add(const std::vector<int>& vars, int rhs) {  // vars の変数の XOR = rhs（同じ変数が 2 回あれば消える）
        std::vector<u64> row((size_t)w_, 0);
        for (int v : vars) {
            if (v < 0 || v >= m_) math_internal::fail("Gf2System::add", "変数の番号は 0〜m−1");
            row[(size_t)(v >> 6)] ^= 1ULL << (v & 63);
        }
        if (rhs & 1) row[(size_t)(m_ >> 6)] ^= 1ULL << (m_ & 63);
        rows_.push_back(std::move(row)), dirty_ = true;
    }
    void add_row(const std::vector<int>& coef, int rhs) {  // Σ coef_j x_j = rhs（coef は長さ m の 0/1）
        if ((int)coef.size() != m_) math_internal::fail("Gf2System::add_row", "係数の長さは m");
        std::vector<int> vars;
        for (int j = 0; j < m_; ++j)
            if (coef[(size_t)j] & 1) vars.push_back(j);
        add(vars, rhs);
    }
    int rank() { return eliminate(), (int)piv_.size(); }
    bool solvable() { return eliminate(), ok_; }
    int free_count() { return eliminate(), m_ - (int)piv_.size(); }  // 解があれば解の数は 2^free_count
    Maybe<std::vector<int>> solve() {                                // 解の 1 つ（自由な変数は 0）
        eliminate();
        if (!ok_) return {false, {}};
        std::vector<int> x((size_t)m_, 0);
        for (size_t i = 0; i < piv_.size(); ++i) x[(size_t)piv_[i]] = bit(rows_[i], m_);
        return {true, x};
    }
    std::vector<std::vector<int>> kernel() {  // 右辺を 0 にした方程式の解の基底（解は solve() + kernel の XOR で全部）
        eliminate();
        std::vector<bool> is_piv((size_t)m_, false);
        for (int c : piv_) is_piv[(size_t)c] = true;
        std::vector<std::vector<int>> res;
        for (int f = 0; f < m_; ++f) {
            if (is_piv[(size_t)f]) continue;
            std::vector<int> v((size_t)m_, 0);
            v[(size_t)f] = 1;
            for (size_t i = 0; i < piv_.size(); ++i) v[(size_t)piv_[i]] = bit(rows_[i], f);
            res.push_back(v);
        }
        return res;
    }
    int equations() const { return (int)rows_.size(); }

  private:
    int m_, w_;
    std::vector<std::vector<u64>> rows_;
    std::vector<int> piv_;
    bool dirty_ = true, ok_ = true;
    static int bit(const std::vector<u64>& r, int j) { return (int)(r[(size_t)(j >> 6)] >> (j & 63) & 1); }
    void eliminate() {  // 行の簡約階段形に
        if (!dirty_) return;
        dirty_ = false, piv_.clear();
        const int n = (int)rows_.size();
        int r = 0;
        for (int c = 0; c < m_ && r < n; ++c) {
            const int wc = c >> 6;
            const u64 bc = 1ULL << (c & 63);
            int p = r;
            while (p < n && !(rows_[(size_t)p][(size_t)wc] & bc)) ++p;
            if (p == n) continue;
            std::swap(rows_[(size_t)p], rows_[(size_t)r]);
            const u64* pr = rows_[(size_t)r].data();
            for (int i = 0; i < n; ++i) {
                if (i == r || !(rows_[(size_t)i][(size_t)wc] & bc)) continue;
                u64* ri = rows_[(size_t)i].data();
                for (int k = wc; k < w_; ++k) ri[k] ^= pr[k];
            }
            piv_.push_back(c), ++r;
        }
        ok_ = true;
        for (int i = r; i < n; ++i)
            if (bit(rows_[(size_t)i], m_)) ok_ = false;
    }
};
#endif  // GG_XOR_BASIS
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace xorb_test {

#define XORB_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[xor_basis] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

using u64 = unsigned long long;

static std::vector<u64> span_of(const std::vector<u64>& a) {  // 部分集合の XOR を全部（小さい順・重複なし）
    std::set<u64> s{0};
    for (u64 x : a) {
        std::set<u64> t = s;
        for (u64 y : s) t.insert(x ^ y);
        s.swap(t);
    }
    return std::vector<u64>(s.begin(), s.end());
}

static void test_basis(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 800 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 12), bits = (int)rnd(rng, 1, 7);
        std::vector<u64> a((size_t)n);
        for (auto& x : a) x = rng() & ((1ULL << bits) - 1);
        XorBasis B;
        std::vector<u64> pre;
        for (u64 x : a) {
            const size_t before = span_of(pre).size();
            pre.push_back(x);
            XORB_CHECK(B.insert(x) == (span_of(pre).size() > before));
        }
        const std::vector<u64> sp = span_of(a);
        ++cases;
        XORB_CHECK(B.count() == (Int128)sp.size() && (1ULL << B.rank()) == sp.size() && B.inserted() == n);
        XORB_CHECK(B.zero_by_nonempty() == (n > B.rank()));
        for (size_t k = 0; k < sp.size(); ++k) XORB_CHECK(B.kth(k) == sp[k]);
        for (u64 x = 0; x < (1ULL << bits) + 3; ++x, ++cases) {
            u64 mn = ~0ULL, mx = 0;
            for (u64 s : sp) mn = std::min(mn, x ^ s), mx = std::max(mx, x ^ s);
            XORB_CHECK(B.contains(x) == std::binary_search(sp.begin(), sp.end(), x));
            XORB_CHECK(B.min_xor(x) == mn && B.reduce(x) == mn && B.max_xor(x) == mx);
            XORB_CHECK(B.count_le(x) == (Int128)(std::upper_bound(sp.begin(), sp.end(), x) - sp.begin()));
            const Maybe<std::vector<int>> r = B.represent(x);
            XORB_CHECK(r.ok == B.contains(x));
            if (r.ok) {
                u64 y = 0;
                for (size_t i = 0; i < r.value.size(); ++i) {
                    XORB_CHECK(r.value[i] >= 0 && r.value[i] < n && (i == 0 || r.value[i - 1] < r.value[i]));
                    y ^= a[(size_t)r.value[i]];
                }
                XORB_CHECK(y == x);
            }
        }
        // merge：2 つに分けて作っても同じ
        XorBasis C, D;
        for (int i = 0; i < n; ++i) (i % 2 ? C : D).insert(a[(size_t)i]);
        C.merge(D);
        XORB_CHECK(C.basis() == B.basis());
    }
    // 64 ビットいっぱいの値
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        XorBasis B;
        const int n = (int)rnd(rng, 0, 80);
        std::vector<u64> a((size_t)n);
        for (auto& x : a) x = it % 3 ? rng() : (rng() & rng() & rng());
        for (u64 x : a) B.insert(x);
        const std::vector<u64> bs = B.basis();
        XORB_CHECK((int)bs.size() == B.rank());
        for (size_t i = 0; i < bs.size(); ++i) {  // 簡約：先頭ビットは他の基底に出てこない
            const int h = 63 - __builtin_clzll(bs[i]);
            for (size_t j = 0; j < bs.size(); ++j) XORB_CHECK(i == j || !(bs[j] >> h & 1));
        }
        const u64 x = rng();
        const Int128 c = B.count_le(x);
        XORB_CHECK(c >= 1 && B.kth((u64)(c - 1)) <= x && (c == B.count() || B.kth((u64)c) > x));
        const u64 k = B.rank() == 64 ? rng() : rng() % (1ULL << B.rank());
        XORB_CHECK(B.count_le(B.kth(k)) == (Int128)k + 1);
        u64 s = 0;  // 作れる値 → represent
        for (u64 v : a)
            if (rng() & 1) s ^= v;
        const Maybe<std::vector<int>> r = B.represent(s);
        XORB_CHECK(r.ok);
        u64 y = 0;
        for (int i : r.value) y ^= a[(size_t)i];
        XORB_CHECK(y == s && B.contains(s) && B.max_xor(s) >= s);
    }
}

static void test_prefix(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 40 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 40), bits = (int)rnd(rng, 1, 64);
        std::vector<u64> a((size_t)n);
        for (auto& x : a) x = bits == 64 ? rng() : rng() & ((1ULL << bits) - 1);
        PrefixXorBasis P;
        for (int r = 1; r <= n; ++r) {
            P.push_back(a[(size_t)r - 1]);
            XORB_CHECK(P.size() == r);
            for (int l = 0; l <= r; ++l, ++cases) {
                XorBasis B;
                for (int i = l; i < r; ++i) B.insert(a[(size_t)i]);
                const u64 x = rng() % 3 ? (rng() & (bits == 64 ? ~0ULL : (1ULL << bits) - 1)) : B.kth(B.rank() == 64 ? rng() : rng() % (1ULL << B.rank()));
                XORB_CHECK(P.contains(l, x) == B.contains(x) && P.max_xor(l, x) == B.max_xor(x) && P.rank(l) == B.rank());
                XORB_CHECK(P.basis(l).basis() == B.basis());
            }
        }
    }
}

static void test_gf2(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int m = (int)rnd(rng, 0, 10), n = (int)rnd(rng, 0, 12);
        Gf2System S(m);
        std::vector<std::vector<int>> eq((size_t)n);
        std::vector<int> rhs((size_t)n);
        for (int i = 0; i < n; ++i) {
            eq[(size_t)i].assign((size_t)m, 0);
            for (auto& c : eq[(size_t)i]) c = (int)(rng() % 3 == 0);
            rhs[(size_t)i] = (int)(rng() & 1);
            if (i % 2) S.add_row(eq[(size_t)i], rhs[(size_t)i]);
            else {
                std::vector<int> vars;
                for (int j = 0; j < m; ++j)
                    if (eq[(size_t)i][(size_t)j]) vars.push_back(j), vars.push_back(j), vars.push_back(j);  // 3 回 = 1 回
                S.add(vars, rhs[(size_t)i]);
            }
        }
        auto sat = [&](const std::vector<int>& x, bool homo) {
            for (int i = 0; i < n; ++i) {
                int s = 0;
                for (int j = 0; j < m; ++j) s ^= eq[(size_t)i][(size_t)j] & x[(size_t)j];
                if (s != (homo ? 0 : rhs[(size_t)i])) return false;
            }
            return true;
        };
        long long cnt = 0;
        for (int mask = 0; mask < (1 << m); ++mask) {
            std::vector<int> x((size_t)m);
            for (int j = 0; j < m; ++j) x[(size_t)j] = mask >> j & 1;
            cnt += sat(x, false);
        }
        XORB_CHECK(S.solvable() == (cnt > 0) && S.equations() == n);
        if (cnt > 0) XORB_CHECK(cnt == (1LL << S.free_count()) && S.rank() == m - S.free_count());
        const Maybe<std::vector<int>> x = S.solve();
        XORB_CHECK(x.ok == (cnt > 0) && (!x.ok || sat(x.value, false)));
        const auto K = S.kernel();
        XORB_CHECK((int)K.size() == m - S.rank());
        XorBasis kb;  // 基底が一次独立か（m ≤ 10 なのでビットにして）
        for (const auto& v : K) {
            XORB_CHECK(sat(v, true));
            u64 b = 0;
            for (int j = 0; j < m; ++j) b |= (u64)v[(size_t)j] << j;
            XORB_CHECK(kb.insert(b));
        }
    }
    // 大きい：答えを決めてから式を作る
    for (int it = 0; it < 20 * TEST_SCALE; ++it, ++cases) {
        const int m = (int)rnd(rng, 1, 400), n = (int)rnd(rng, 1, 450);
        std::vector<int> ans((size_t)m);
        for (auto& v : ans) v = (int)(rng() & 1);
        Gf2System S(m);
        std::vector<std::vector<int>> eq;
        std::vector<int> rhs;
        for (int i = 0; i < n; ++i) {
            std::vector<int> c((size_t)m);
            int s = 0;
            for (int j = 0; j < m; ++j) c[(size_t)j] = (int)(rng() & 1), s ^= c[(size_t)j] & ans[(size_t)j];
            S.add_row(c, s), eq.push_back(c), rhs.push_back(s);
        }
        const Maybe<std::vector<int>> x = S.solve();
        XORB_CHECK(x.ok);
        for (int i = 0; i < n; ++i) {
            int s = 0;
            for (int j = 0; j < m; ++j) s ^= eq[(size_t)i][(size_t)j] & x.value[(size_t)j];
            XORB_CHECK(s == rhs[(size_t)i]);
        }
        S.add({0}, 1 - x.value[0]);  // 式を足すと解き直す
        XORB_CHECK(S.equations() == n + 1);
        if (S.solvable()) XORB_CHECK(S.solve().value[0] == 1 - x.value[0]);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_basis(rng, cases);
    test_prefix(rng, cases);
    test_gf2(rng, cases);
    std::cout << "[xor_basis] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[xor_basis] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    XorBasis B;
    for (int i = 0; i < 1000000; ++i) B.insert(rng() >> 4);
    u64 s1 = 0;
    for (int i = 0; i < 1000000; ++i) s1 ^= B.max_xor(rng()) ^ B.reduce(rng());
    const long long t_b = ms(t0);
    t0 = std::chrono::steady_clock::now();
    PrefixXorBasis P;
    u64 s2 = 0;
    for (int i = 0; i < 400000; ++i) {
        P.push_back(rng() >> 4);
        if (i % 2) s2 += P.contains((int)rnd(rng, 0, i), rng() >> 4);
    }
    const long long t_p = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const int M = 2000;
    Gf2System S(M);
    for (int i = 0; i < M; ++i) {
        std::vector<int> c((size_t)M);
        for (auto& v : c) v = (int)(rng() & 1);
        S.add_row(c, (int)(rng() & 1));
    }
    const int rk = S.rank();
    const long long t_g = ms(t0);
    std::cout << "[xor_basis] benchmark: XorBasis insert × 10^6 + クエリ × 2×10^6 " << t_b << " ms, PrefixXorBasis（4×10^5 個 + クエリ 2×10^5）" << t_p
              << " ms, Gf2System 2000 × 2000 " << t_g << " ms  (check " << (s1 + s2 + (u64)rk) % 7 << ")\n";
}

}  // namespace xorb_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace xorb_demo {

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
#define SHOW(expr, comment) xorb_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    xorb_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("XOR の基底（xor_basis） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  XorBasis：insert / contains / max_xor / min_xor / kth / count_le / rank / count / zero_by_nonempty / represent / merge / basis\n"
              << "  PrefixXorBasis：push_back / contains(l, x) / max_xor(l, x) / rank(l) / basis(l)\n"
              << "  Gf2System：add / add_row / rank / solvable / free_count / solve / kernel\n";

    api("【1】 XorBasis B; B.insert(x) / contains / max_xor / min_xor / rank / count", "bool / bool / u64 / u64 / int / Int128",
        "XOR で作れる値（0 個の XOR = 0 も含む）。insert は作れる値が増えたら true", "O(64)");
    RUN(XorBasis B);
    SHOW(B.insert(2), "");
    SHOW(B.insert(21), "");
    SHOW(B.insert(17), "");
    SHOW(B.insert(21), "21 は作れる（増えない）");
    SHOW(B.contains(6), "2 ^ 21 ^ 17 = 6");
    SHOW(B.contains(1), "");
    SHOW(B.max_xor(), "作れる値の最大");
    SHOW(B.max_xor(8), "8 XOR (作れる値) の最大");
    SHOW(B.min_xor(23), "23 XOR (作れる値) の最小");
    SHOW(B.rank(), "");
    SHOW(B.count(), "2^3 = 8 個");
    SHOW(B.zero_by_nonempty(), "21 ^ 21 = 0（1 個以上の XOR で 0 を作れる）");
    {
        // ABC141 F（Xor Sum 3）：全体の XOR T で立っているビットは必ず 1 回だけ足される。残りのビットは 2 倍 → T を消した値の max_xor
        auto solve = [](const std::vector<unsigned long long>& A) {
            unsigned long long T = 0;
            for (auto a : A) T ^= a;
            XorBasis Bs;
            for (auto a : A) Bs.insert(a & ~T);
            return T + 2 * Bs.max_xor();
        };
        note("ABC141 F（Xor Sum 3）: 全体の XOR を T として、各 A_i から T のビットを消した基底の max_xor × 2 + T");
        std::cout << "    例 1 → " << solve({3, 6, 5}) << "、例 2 → " << solve({23, 36, 66, 65}) << "、例 3 → "
                  << solve({1008288677408720767ULL, 539403903321871999ULL, 1044301017184589821ULL, 215886900497862655ULL, 504277496111605629ULL,
                            972104334925272829ULL, 792625803473366909ULL, 972333547668684797ULL, 467386965442856573ULL, 755861732751878143ULL,
                            1151846447448561405ULL, 467257771752201853ULL, 683930041385277311ULL, 432010719984459389ULL, 319104378117934975ULL,
                            611451291444233983ULL, 647509226592964607ULL, 251832107792119421ULL, 827811265410084479ULL, 864032478037725181ULL})
                  << "\n";
    }

    api("【2】 B.kth(k) / B.count_le(x)", "u64 / Int128", "作れる値を小さい順に並べた k 番目（0 から。0 番目は 0）/ x 以下の作れる値の個数", "O(64)");
    SHOW(B.kth(0), "");
    SHOW(B.kth(7), "最大");
    SHOW(B.count_le(18), "0, 2, 4, 6, 17 の 5 個");
    {
        // ABC283 G（Partial Xor Enumeration）：s_L〜s_R（1 から）= kth(L−1)〜kth(R−1)
        auto solve = [](long long L, long long R, const std::vector<unsigned long long>& A) {
            XorBasis Bs;
            for (auto a : A) Bs.insert(a);
            std::vector<unsigned long long> res;
            for (long long k = L - 1; k < R; ++k) res.push_back(Bs.kth((unsigned long long)k));
            return res;
        };
        note("ABC283 G（Partial Xor Enumeration）: kth(L − 1)〜kth(R − 1)");
        std::cout << "    例 1 → " << to_s(solve(1, 8, {2, 21, 17, 21})) << "\n    例 2 → " << to_s(solve(3, 7, {2, 21, 17, 21})) << "\n    例 3 → "
                  << to_s(solve(1, 1, {0, 0, 0, 0, 0})) << "、例 4 → " << to_s(solve(21, 21, {88, 44, 82, 110, 121, 80})) << "\n";
    }

    api("【3】 B.represent(x) / B.basis() / B.merge(C)", "Maybe<vector<int>> / vector<u64> / なし",
        "x を作る要素の番号（insert を呼んだ順に 0 から）/ 基底（簡約済み。大きい順）/ C で作れる値も加える", "O(64)");
    SHOW(B.represent(6).value, "2 ^ 21 ^ 17 → insert の 0, 1, 2 番目");
    SHOW(B.represent(1).ok, "作れない");
    SHOW(B.basis(), "簡約した基底");

    api("【4】 PrefixXorBasis P; P.push_back(a) / P.contains(l, x) / P.max_xor(l, x) / P.rank(l)", "bool / u64 / int",
        "a_0, a_1, … を順に足し、今の末尾 r = size() までの区間 [l, r) の要素で（クエリは r の順に並べて答える）", "O(64)");
    {
        // ABC223 H（Xor Query）：クエリを R の順に並べ、A_R まで足したところで contains(L − 1, X)
        auto solve = [](const std::vector<unsigned long long>& A, const std::vector<std::array<long long, 3>>& Q) {
            std::vector<int> ord(Q.size());
            std::iota(ord.begin(), ord.end(), 0);
            std::sort(ord.begin(), ord.end(), [&](int a, int b) { return Q[(size_t)a][1] < Q[(size_t)b][1]; });
            PrefixXorBasis P;
            std::vector<std::string> ans(Q.size());
            for (int q : ord) {
                while (P.size() < Q[(size_t)q][1]) P.push_back(A[(size_t)P.size()]);
                ans[(size_t)q] = P.contains((int)Q[(size_t)q][0] - 1, (unsigned long long)Q[(size_t)q][2]) ? "Yes" : "No";
            }
            return ans;
        };
        note("ABC223 H（Xor Query）: クエリを R の順に並べて、A_1〜A_R を push_back したところで contains(L − 1, X)");
        std::cout << "    例 1 → " << to_s(solve({3, 1, 4, 1, 5}, {{1, 3, 7}, {2, 5, 7}})) << "\n    例 2 → "
                  << to_s(solve({8, 45, 56, 9, 38, 28, 33, 5, 15, 19}, {{10, 10, 53}, {3, 8, 60}, {1, 10, 29}, {5, 7, 62}, {3, 7, 51}, {8, 8, 52}, {1, 4, 60}, {6, 8, 32}, {4, 8, 58}, {5, 9, 2}}))
                  << "\n";
    }

    api("【5】 Gf2System S(m); S.add(vars, rhs) / S.solve() / S.kernel() / S.free_count()", "なし / Maybe<vector<int>> / vector<vector<int>> / int",
        "mod 2 の連立一次方程式（vars の変数の XOR = rhs）。解の 1 つと、右辺 0 の解の基底", "O(式 × rank × m / 64)");
    RUN(Gf2System S(3));
    RUN(S.add({0, 1}, 1));
    RUN(S.add({1, 2}, 0));
    SHOW(S.solve().value, "x0 ^ x1 = 1、x1 ^ x2 = 0 の解の 1 つ（x2 = 0）");
    SHOW(S.kernel(), "(1, 1, 1) を XOR しても解");
    SHOW(S.free_count(), "解は 2^1 = 2 個");
    RUN(S.add({0, 2}, 0));
    SHOW(S.solvable(), "x0 ^ x2 = 0 を足すと x0 ^ x1 ^ x1 ^ x2 = 1 と矛盾");
}

}  // namespace xorb_demo

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
    xorb_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    xorb_test::benchmark();    // 速度計測
    xorb_demo::demo();         // 使い方デモ
    return 0;
}
