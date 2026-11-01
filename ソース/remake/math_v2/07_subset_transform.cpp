// =============================================================================
//  部分集合の変換（高速ゼータ・メビウス変換、OR・AND・XOR の畳み込み、部分集合の畳み込み・exp・log）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 07 畳み込み・多項式
//  依存       : 00_int_common, 00_math_base, 00_modint（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  配列の長さは 2^n（i 番目が集合 i：i ビット目が 1 ⇔ 要素 i を含む）。T は mint・long long（溢れたら止める）・double など
//  （modint は数え上げで使うことが多いので同梱している）
//  long long は途中の値（ゼータ変換の和など）が溢れても止まる。答えだけ収まる場合もあるので、大きくなりそうなら mint で
//  xor_convolution・walsh_hadamard(f, true) は 2^n で割るので、mint の mod は奇数（DynModInt の偶数の mod では使えない）
//
//  機能一覧                                                                    計算量
//    subset_zeta(f, op) / subset_mobius(f, inv_op)       f[S] ← Σ_{T ⊆ S} f[T] / その逆（op を渡すと max・gcd などでも） O(n 2^n)
//    superset_zeta(f, op) / superset_mobius(f, inv_op)   f[S] ← Σ_{T ⊇ S} f[T] / その逆
//    walsh_hadamard(f, inverse)       アダマール変換（XOR の畳み込み用。inverse なら 2^n で割る）
//    or_convolution(a, b)             c[k] = Σ_{i | j = k} a[i] b[j]（長さは 2 の冪にそろえる）
//    and_convolution(a, b)            c[k] = Σ_{i & j = k} a[i] b[j]
//    xor_convolution(a, b)            c[k] = Σ_{i ^ j = k} a[i] b[j]
//    subset_convolution(a, b)         c[S] = Σ_{T ⊆ S} a[T] b[S∖T]（集合を 2 つに分ける）                O(n^2 2^n)
//    subset_exp(f) / subset_log(f)    g[S] = Σ_{S の組分け} Π f[組]（f[∅] = 0）/ その逆（f[∅] = 1）。連結な部分の数え上げ  O(n^2 2^n)
//                                     （割り算を使わないので long long・合成数の mod でも。2^18 で exp 約 0.3 秒、log 約 0.5 秒）
//  （部分集合の列挙は 09_bit、約数・倍数のゼータ変換は 02_divisor_transform）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    subt_test::random_test()  … 全部の組を調べる愚直解との突き合わせ
//    subt_test::benchmark()    … 速度計測
//    subt_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_SUBSET_TRANSFORM
#define GG_SUBSET_TRANSFORM
// ---- 部分集合の変換（数学第二世代）  高速ゼータ・メビウス変換、OR・AND・XOR の畳み込み、部分集合の畳み込み・exp・log ----
// 配列の長さは 2^n（i 番目が集合 i）。T は mint・long long（溢れたら止める）・double など
namespace math_internal {
template <class T>
inline T st_add(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail("subset_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a + b;
    }
}
template <class T>
inline T st_sub(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_sub_overflow(a, b, &r)) fail("subset_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a - b;
    }
}
template <class T>
inline T st_mul(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_mul_overflow(a, b, &r)) fail("subset_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a * b;
    }
}
inline int st_log2(size_t n, const char* func) {  // n = 2^k の k（2 の冪でなければ止める）
    if (n == 0 || (n & (n - 1))) fail(func, "長さは 2 の冪（2^n）");
    return __builtin_ctzll((unsigned long long)n);
}
struct StAdd {
    template <class T>
    T operator()(const T& a, const T& b) const { return st_add(a, b); }
};
struct StSub {
    template <class T>
    T operator()(const T& a, const T& b) const { return st_sub(a, b); }
};
}  // namespace math_internal

// f[S] ← Σ_{T ⊆ S} f[T]（op を渡すと Σ の代わりに op でまとめる：max・min・gcd・OR など）  O(n 2^n)
template <class T, class Op = math_internal::StAdd>
void subset_zeta(std::vector<T>& f, Op op = Op{}) {
    const size_t n = f.size();
    math_internal::st_log2(n, "subset_zeta");
    for (size_t w = 1; w < n; w <<= 1)
        for (size_t s = 0; s < n; ++s)
            if (s & w) f[s] = op(f[s], f[s ^ w]);
}
// subset_zeta の逆（f[S] ← Σ_{T ⊆ S} (−1)^{|S∖T|} f[T]。inv_op は op の逆）
template <class T, class Op = math_internal::StSub>
void subset_mobius(std::vector<T>& f, Op inv_op = Op{}) {
    const size_t n = f.size();
    math_internal::st_log2(n, "subset_mobius");
    for (size_t w = 1; w < n; w <<= 1)
        for (size_t s = 0; s < n; ++s)
            if (s & w) f[s] = inv_op(f[s], f[s ^ w]);
}
// f[S] ← Σ_{T ⊇ S} f[T]
template <class T, class Op = math_internal::StAdd>
void superset_zeta(std::vector<T>& f, Op op = Op{}) {
    const size_t n = f.size();
    math_internal::st_log2(n, "superset_zeta");
    for (size_t w = 1; w < n; w <<= 1)
        for (size_t s = 0; s < n; ++s)
            if (!(s & w)) f[s] = op(f[s], f[s | w]);
}
template <class T, class Op = math_internal::StSub>
void superset_mobius(std::vector<T>& f, Op inv_op = Op{}) {
    const size_t n = f.size();
    math_internal::st_log2(n, "superset_mobius");
    for (size_t w = 1; w < n; w <<= 1)
        for (size_t s = 0; s < n; ++s)
            if (!(s & w)) f[s] = inv_op(f[s], f[s | w]);
}
// アダマール変換（XOR の畳み込み用。inverse なら 2^n で割る：T が整数なら割り切れることを確かめる）
template <class T>
void walsh_hadamard(std::vector<T>& f, bool inverse = false) {
    const size_t n = f.size();
    math_internal::st_log2(n, "walsh_hadamard");
    for (size_t w = 1; w < n; w <<= 1)
        for (size_t s = 0; s < n; ++s)
            if (!(s & w)) {
                const T a = f[s], b = f[s | w];
                f[s] = math_internal::st_add(a, b), f[s | w] = math_internal::st_sub(a, b);
            }
    if (inverse) {
        if constexpr (std::is_integral_v<T>) {
            for (T& x : f) {
                if (x % (T)n != 0) math_internal::fail("walsh_hadamard", "2^n で割り切れない");
                x /= (T)n;
            }
        } else {
            const T inv = T(1) / T((long long)n);
            for (T& x : f) x = x * inv;
        }
    }
}
namespace math_internal {
template <class T>
void st_pad(std::vector<T>& a, std::vector<T>& b) {  // 同じ 2 の冪の長さにそろえる
    size_t n = 1;
    while (n < std::max(a.size(), b.size())) n <<= 1;
    a.resize(n, T(0)), b.resize(n, T(0));
}
}  // namespace math_internal
// c[k] = Σ_{i | j = k} a[i] b[j]（長さは 2 の冪にそろえる）  O(n 2^n)
template <class T>
std::vector<T> or_convolution(std::vector<T> a, std::vector<T> b) {
    if (a.empty() || b.empty()) return {};
    math_internal::st_pad(a, b);
    subset_zeta(a), subset_zeta(b);
    for (size_t i = 0; i < a.size(); ++i) a[i] = math_internal::st_mul(a[i], b[i]);
    subset_mobius(a);
    return a;
}
// c[k] = Σ_{i & j = k} a[i] b[j]
template <class T>
std::vector<T> and_convolution(std::vector<T> a, std::vector<T> b) {
    if (a.empty() || b.empty()) return {};
    math_internal::st_pad(a, b);
    superset_zeta(a), superset_zeta(b);
    for (size_t i = 0; i < a.size(); ++i) a[i] = math_internal::st_mul(a[i], b[i]);
    superset_mobius(a);
    return a;
}
// c[k] = Σ_{i ^ j = k} a[i] b[j]
template <class T>
std::vector<T> xor_convolution(std::vector<T> a, std::vector<T> b) {
    if (a.empty() || b.empty()) return {};
    math_internal::st_pad(a, b);
    walsh_hadamard(a), walsh_hadamard(b);
    for (size_t i = 0; i < a.size(); ++i) a[i] = math_internal::st_mul(a[i], b[i]);
    walsh_hadamard(a, true);
    return a;
}
// 部分集合の畳み込み c[S] = Σ_{T ⊆ S} a[T] b[S∖T]（ランク付きのゼータ変換）  O(n^2 2^n)
template <class T>
std::vector<T> subset_convolution(std::vector<T> a, std::vector<T> b) {
    if (a.empty() || b.empty()) return {};
    math_internal::st_pad(a, b);
    const size_t N = a.size();
    const int n = math_internal::st_log2(N, "subset_convolution");
    std::vector<std::vector<T>> ra((size_t)n + 1, std::vector<T>(N, T(0))), rb = ra, rc = ra;
    for (size_t s = 0; s < N; ++s) {
        const int k = __builtin_popcountll((unsigned long long)s);
        ra[(size_t)k][s] = a[s], rb[(size_t)k][s] = b[s];
    }
    for (int k = 0; k <= n; ++k) subset_zeta(ra[(size_t)k]), subset_zeta(rb[(size_t)k]);
    for (int i = 0; i <= n; ++i)
        for (int j = 0; i + j <= n; ++j)
            for (size_t s = 0; s < N; ++s) rc[(size_t)(i + j)][s] = math_internal::st_add(rc[(size_t)(i + j)][s], math_internal::st_mul(ra[(size_t)i][s], rb[(size_t)j][s]));
    for (int k = 0; k <= n; ++k) subset_mobius(rc[(size_t)k]);
    std::vector<T> c(N, T(0));
    for (size_t s = 0; s < N; ++s) c[s] = rc[(size_t)__builtin_popcountll((unsigned long long)s)][s];
    return c;
}
namespace math_internal {
template <class T>
std::vector<T> st_subset_inv(const std::vector<T>& a) {  // subset_convolution の逆元（a[∅] = 1。割り算はいらない）  O(n^2 2^n)
    const size_t N = a.size();
    const int n = st_log2(N, "subset_log");
    std::vector<std::vector<T>> r((size_t)n + 1, std::vector<T>(N, T(0)));
    for (size_t s = 0; s < N; ++s) r[(size_t)__builtin_popcountll((unsigned long long)s)][s] = a[s];
    for (int k = 0; k <= n; ++k) subset_zeta(r[(size_t)k]);
    std::vector<T> q((size_t)n + 1);
    for (size_t s = 0; s < N; ++s) {  // 各 s で多項式 Σ_k r[k][s] z^k（定数項 1）の逆数
        q[0] = T(1);
        for (int m = 1; m <= n; ++m) {
            T acc = T(0);
            for (int k = 1; k <= m; ++k) acc = st_add(acc, st_mul(r[(size_t)k][s], q[(size_t)(m - k)]));
            q[(size_t)m] = st_sub(T(0), acc);
        }
        for (int m = 0; m <= n; ++m) r[(size_t)m][s] = q[(size_t)m];
    }
    for (int k = 0; k <= n; ++k) subset_mobius(r[(size_t)k]);
    std::vector<T> c(N);
    for (size_t s = 0; s < N; ++s) c[s] = r[(size_t)__builtin_popcountll((unsigned long long)s)][s];
    return c;
}
}  // namespace math_internal
// 集合のべき級数の exp：g[S] = Σ_{S の組分け {B_1, …, B_k}} f[B_1] … f[B_k]（f[∅] = 0、g[∅] = 1。長さ 2^n）
//   「つながった部分の数え方」f から「全体の数え方」g へ。割り算を使わないので long long・合成数の mod でも  O(n^2 2^n)
template <class T>
std::vector<T> subset_exp(const std::vector<T>& f) {
    const int n = math_internal::st_log2(f.size(), "subset_exp");
    if (!(f[0] == T(0))) math_internal::fail("subset_exp", "f[∅] = 0 で使う");
    std::vector<T> g = {T(1)};
    for (int k = 0; k < n; ++k) {  // 要素 k を含む集合：k を含む組 {k} ∪ T と残りの組分け
        const std::vector<T> h(f.begin() + (1LL << k), f.begin() + (2LL << k));
        const std::vector<T> hi = subset_convolution(h, g);
        g.insert(g.end(), hi.begin(), hi.end());
    }
    return g;
}
// subset_exp の逆：f[∅] = 1 の f から g（f = subset_exp(g)、g[∅] = 0）。「全体」から「つながったもの」を数える  O(n^2 2^n)
template <class T>
std::vector<T> subset_log(const std::vector<T>& f) {
    const int n = math_internal::st_log2(f.size(), "subset_log");
    if (!(f[0] == T(1))) math_internal::fail("subset_log", "f[∅] = 1 で使う");
    std::vector<T> g = {T(0)};
    for (int k = 0; k < n; ++k) {  // f の上半分 = g の上半分 * f の下半分 なので、f の下半分の逆元を掛ける
        const std::vector<T> lo(f.begin(), f.begin() + (1LL << k)), hi(f.begin() + (1LL << k), f.begin() + (2LL << k));
        const std::vector<T> gh = subset_convolution(hi, math_internal::st_subset_inv(lo));
        g.insert(g.end(), gh.begin(), gh.end());
    }
    return g;
}
#endif  // GG_SUBSET_TRANSFORM
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace subt_test {

#define SUBT_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[subset_transform] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

template <class T, class Gen>
static void test_type(std::mt19937_64& rng, int& cases, int iters, Gen gen) {
    for (int it = 0; it < iters; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 7);
        const size_t N = (size_t)1 << n;
        std::vector<T> a(N), b(N);
        for (auto& x : a) x = gen(rng);
        for (auto& x : b) x = gen(rng);
        // ゼータ・メビウス
        std::vector<T> z = a, zs = a, e(N, T(0)), es(N, T(0));
        for (size_t s = 0; s < N; ++s)
            for (size_t t = 0; t < N; ++t) {
                if ((t & s) == t) e[s] = e[s] + a[t];
                if ((t & s) == s) es[s] = es[s] + a[t];
            }
        subset_zeta(z), superset_zeta(zs);
        SUBT_CHECK(z == e && zs == es);
        subset_mobius(z), superset_mobius(zs);
        SUBT_CHECK(z == a && zs == a);
        // 畳み込み
        std::vector<T> co(N, T(0)), ca(N, T(0)), cx(N, T(0)), cs(N, T(0));
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j) {
                co[i | j] = co[i | j] + a[i] * b[j], ca[i & j] = ca[i & j] + a[i] * b[j], cx[i ^ j] = cx[i ^ j] + a[i] * b[j];
                if ((i & j) == 0) cs[i | j] = cs[i | j] + a[i] * b[j];
            }
        SUBT_CHECK(or_convolution(a, b) == co && and_convolution(a, b) == ca && xor_convolution(a, b) == cx && subset_convolution(a, b) == cs);
    }
}

template <class T>
static std::vector<T> brute_exp(const std::vector<T>& f) {  // g[S] = Σ_{T ∋ S の最小の要素、T ⊆ S} f[T] g[S∖T]
    const size_t N = f.size();
    std::vector<T> g(N, T(0));
    g[0] = T(1);
    for (size_t s = 1; s < N; ++s) {
        const size_t low = s & (~s + 1), rest = s ^ low;
        for (size_t t = rest;; t = (t - 1) & rest) {  // 最小の要素を含む組 low | t
            g[s] = g[s] + f[low | t] * g[rest ^ t];
            if (t == 0) break;
        }
    }
    return g;
}
static void test_exp_log(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const size_t N = (size_t)1 << rnd(rng, 0, 9);
        std::vector<ModInt998> f(N);
        for (auto& x : f) x = (long long)rng();
        f[0] = 0;
        const std::vector<ModInt998> g = subset_exp(f);
        SUBT_CHECK(g == brute_exp(f) && subset_log(g) == f);
        std::vector<long long> fl(N);  // long long でも（小さい値）
        for (auto& x : fl) x = rnd(rng, -3, 3);
        fl[0] = 0;
        const std::vector<long long> gl = subset_exp(fl);
        SUBT_CHECK(gl == brute_exp(fl) && subset_log(gl) == fl);
    }
    // 全部 1 の exp はベル数、2^{C(|S|,2)}（グラフの数）の log は連結なグラフの数
    const int n = 6;
    std::vector<long long> one((size_t)1 << n, 1), graphs((size_t)1 << n);
    one[0] = 0;
    for (size_t s = 0; s < graphs.size(); ++s) {
        const int k = __builtin_popcountll(s);
        graphs[s] = 1LL << (k * (k - 1) / 2);
    }
    const std::vector<long long> bell = subset_exp(one), conn = subset_log(graphs);
    SUBT_CHECK(bell.back() == 203 && conn[1] == 1 && conn[3] == 1 && conn[7] == 4 && conn[15] == 38 && conn[31] == 728 && conn.back() == 26704);
    DynModInt::set_mod(1000000000);  // 合成数の mod でも（割り算を使わない）
    std::vector<DynModInt> d((size_t)1 << 7);
    for (auto& x : d) x = (long long)rng();
    d[0] = 0;
    SUBT_CHECK(subset_log(subset_exp(d)) == d);
    cases += 2;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_type<ModInt998>(rng, cases, 2000 * TEST_SCALE, [](std::mt19937_64& r) { return ModInt998((long long)r()); });
    test_type<long long>(rng, cases, 2000 * TEST_SCALE, [](std::mt19937_64& r) { return rnd(r, -1000, 1000); });
    DynModInt::set_mod(1000000000);  // 合成数・偶数の mod でも（XOR の逆変換は 2^n の逆元が要るので除く）
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {
        const size_t N = (size_t)1 << rnd(rng, 0, 6);
        std::vector<DynModInt> a(N), b(N), co(N), cs(N);
        for (auto& x : a) x = (long long)rng();
        for (auto& x : b) x = (long long)rng();
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j) {
                co[i | j] += a[i] * b[j];
                if ((i & j) == 0) cs[i | j] += a[i] * b[j];
            }
        SUBT_CHECK(or_convolution(a, b) == co && subset_convolution(a, b) == cs);
    }
    // 演算を渡す：部分集合の最大値
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {
        const size_t N = (size_t)1 << rnd(rng, 0, 8);
        std::vector<long long> a(N), e(N, LLONG_MIN);
        for (auto& x : a) x = rnd(rng, -100, 100);
        for (size_t s = 0; s < N; ++s)
            for (size_t t = 0; t < N; ++t)
                if ((t & s) == t) e[s] = std::max(e[s], a[t]);
        subset_zeta(a, [](long long x, long long y) { return std::max(x, y); });
        SUBT_CHECK(a == e);
    }
    // 長さがちがうものは 2 の冪にそろえる
    SUBT_CHECK(xor_convolution<long long>({1, 2, 3}, {4, 5}) == std::vector<long long>({4 + 10, 5 + 8, 12 + 0, 15 + 0}));
    test_exp_log(rng, cases);
    std::cout << "[subset_transform] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[subset_transform] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const size_t N = (size_t)1 << 20;
    std::vector<ModInt998> a(N), b(N);
    for (auto& x : a) x = (long long)rng();
    for (auto& x : b) x = (long long)rng();
    auto t0 = std::chrono::steady_clock::now();
    const auto c1 = xor_convolution(a, b);
    const auto c2 = or_convolution(a, b);
    const long long t1 = ms(t0);
    std::vector<ModInt998> sa(a.begin(), a.begin() + (1 << 18)), sb(b.begin(), b.begin() + (1 << 18));
    t0 = std::chrono::steady_clock::now();
    const auto c3 = subset_convolution(sa, sb);
    const long long t2 = ms(t0);
    std::vector<ModInt998> e(sa);
    e[0] = 0;
    t0 = std::chrono::steady_clock::now();
    const auto ex = subset_exp(e);
    const long long t3 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto lg = subset_log(ex);
    const long long t4 = ms(t0);
    std::cout << "[subset_transform] benchmark: xor_convolution + or_convolution（2^20）" << t1 << " ms, subset_convolution（2^18）" << t2 << " ms, subset_exp（2^18）"
              << t3 << " ms, subset_log（2^18）" << t4 << " ms  (check " << (c1[5] + c2[7] + c3[9]).val() % 7 << (lg == e ? "" : " log(exp) ≠ f") << ")\n";
}

}  // namespace subt_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace subt_demo {

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
#define SHOW(expr, comment) subt_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    subt_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("部分集合の変換（subset_transform） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  subset_zeta / subset_mobius / superset_zeta / superset_mobius（f を書きかえる。op を渡せる）\n"
              << "  walsh_hadamard / or_convolution / and_convolution / xor_convolution / subset_convolution\n";

    api("【1】 subset_zeta(f) / subset_mobius(f) / superset_zeta(f) / superset_mobius(f)", "なし（f を書きかえる）",
        "f[S] ← Σ_{T ⊆ S} f[T] / その逆 / ⊇ の版。op を渡すと max・min・gcd などでまとめる", "O(n 2^n)");
    RUN(std::vector<long long> f = {1, 2, 3, 4}; subset_zeta(f));
    SHOW(f, "f[3] = f[0] + f[1] + f[2] + f[3]");
    RUN(subset_mobius(f));
    SHOW(f, "元に戻る");
    RUN(std::vector<long long> g = {5, 1, 7, 2}; subset_zeta(g, [](long long x, long long y) { return std::max(x, y); }));
    SHOW(g, "部分集合の最大値");
    note("bit DP：「S の部分集合すべての和」を全 S で求めるときに。superset は「S を含む集合すべての和」");

    api("【2】 or_convolution(a, b) / and_convolution(a, b) / xor_convolution(a, b)", "vector<T>",
        "c[k] = Σ_{i ∘ j = k} a[i] b[j]（∘ は OR・AND・XOR）", "O(n 2^n)");
    SHOW(or_convolution<long long>({1, 1}, {1, 1}), "{0, 1} から 2 つ選んだ OR：0 が 1 通り、1 が 3 通り");
    SHOW(xor_convolution<long long>({0, 1, 1, 0}, {0, 1, 1, 0}), "{1, 2} から 2 つ選んだ XOR");
    {
        // ABC212 H：A の数を 1〜N 個並べて XOR ≠ 0 になる数。XOR = 0 の数をアダマール変換の上で等比数列の和に
        auto solve = [](long long N, const std::vector<int>& A) {
            using mint = ModInt998;
            std::vector<mint> h(1 << 16);
            for (int a : A) h[(size_t)a] = 1;
            walsh_hadamard(h);
            for (mint& x : h) x = (x == mint(1)) ? mint(N) : x * (x.pow(N) - 1) / (x - 1);  // x + x^2 + … + x^N
            walsh_hadamard(h, true);
            const mint K = (long long)A.size();
            const mint all = (K == mint(1)) ? mint(N) : K * (K.pow(N) - 1) / (K - 1);
            return all - h[0];
        };
        note("ABC212 H（Nim Counting）: アダマール変換して点ごとに x + x^2 + … + x^N、戻して XOR = 0 の数を引く");
        std::cout << "    例 1 → " << solve(2, {1, 2}) << "、例 2 → " << solve(100, {3, 5, 7}) << "\n";
    }

    api("【3】 subset_convolution(a, b)", "vector<T>", "c[S] = Σ_{T ⊆ S} a[T] b[S∖T]（集合 S を 2 つに分けるすべての方法）", "O(n^2 2^n)");
    SHOW(subset_convolution<long long>({1, 1, 1, 1}, {1, 1, 1, 1}), "{0, 1} の部分集合 S を 2 つに分ける方法の数 = 2^|S|");
    note("集合を k 個のグループに分ける数え上げ（bit DP の O(3^n) を O(n^2 2^n) に）で使う");

    api("【4】 subset_exp(f) / subset_log(f)", "vector<T>",
        "g[S] = Σ_{S の組分け} Π f[組]（f[∅] = 0）/ その逆（f[∅] = 1）。割り算なしなので long long・合成数の mod でも", "O(n^2 2^n)");
    SHOW(subset_exp<long long>({0, 1, 1, 1, 1, 1, 1, 1}), "どの組も 1 通り：組分けの数 = ベル数（|S| = 3 で 5）");
    {
        // 頂点の集合ごとの「辺の数 e(S)」から、連結な部分グラフ（全頂点を使う辺の選び方）の数：全体 2^{e(S)} の log
        const int n = 4;
        const std::vector<std::pair<int, int>> es = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}};
        std::vector<ModInt998> all((size_t)1 << n);
        for (size_t s = 0; s < all.size(); ++s) {
            int e = 0;
            for (auto [u, v] : es) e += (s >> u & 1) && (s >> v & 1);
            all[s] = ModInt998(2).pow(e);  // S の中の辺を自由に選ぶ
        }
        const std::vector<ModInt998> conn = subset_log(all);
        note("連結の数え上げ：「S の中の辺を自由に選ぶ」2^{e(S)} を subset_log すると「S が連結になる選び方」");
        std::cout << "    4 頂点（辺 01 12 23 30 02）で全頂点が連結になる辺の選び方 → " << conn.back() << "\n";
    }
}

}  // namespace subt_demo

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
    subt_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    subt_test::benchmark();    // 速度計測
    subt_demo::demo();         // 使い方デモ
    return 0;
}
