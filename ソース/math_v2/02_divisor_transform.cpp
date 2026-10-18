// =============================================================================
//  約数・倍数の変換（約数・倍数のゼータ・メビウス変換、gcd・lcm の畳み込み、gcd ごとの組の数、Σφ・Σμ、ディリクレ積、大きな M の約数の上の変換）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 素数・約数
//  依存       : 00_int_common, 00_math_base, 02_factorize, 02_sieve（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  配列は 1〜N を使う（長さ N + 1。f[0] は使わない）。T は long long（溢れたら止める）・double・mint（00_modint を貼れば）など
//  （long long は途中の値（ゼータ変換の和など）が溢れても止まる。答えだけ収まる場合もあるので、大きくなりそうなら mint で）
//
//  機能一覧                                                                    計算量
//    divisor_zeta(f, op) / divisor_mobius(f, inv_op)     g[n] = Σ_{d | n} f[d] / その逆（op を渡すと max・gcd なども） O(N log log N)
//    multiple_zeta(f, op) / multiple_mobius(f, inv_op)   g[d] = Σ_{d | n} f[n]（倍数の和）/ その逆
//    gcd_convolution(a, b) / lcm_convolution(a, b)       c[k] = Σ_{gcd(i,j)=k} a[i] b[j] / lcm の版（k ≤ N）
//    count_pairs_by_gcd(cnt)          値 v の個数 cnt[v] から、gcd がちょうど g の組 {i < j} の数（全部の g）
//    dirichlet_convolution(a, b)      ディリクレ積 c[n] = Σ_{d|n} a[d] b[n/d]（n = 1〜N）                          O(N log N)
//    DivisorsOf D(M)                  M（≤ 10^18）の約数の上の変換。f は約数の番号で引く（長さ D.size()）          構築 分解 + O(D ω)
//      D.divisors() / D.size() / D[i] / D.index(d) / D.factors()     約数（昇順）/ 個数 / i 番目 / d の番号（約数でなければ -1）/ 素因数分解
//      D.zeta(f, op) / D.mobius(f, inv_op)                g[d] = Σ_{e|d} f[e] / その逆                                   O(D ω(M))
//      D.multiple_zeta(f, op) / D.multiple_mobius(f, inv_op)   g[d] = Σ_{d|e|M} f[e] / その逆
//    phi_sum(N) / mertens(N)          Σ_{i≤N} φ(i)（Int128）/ Σ_{i≤N} μ(i)（N = 10^10 で約 0.3 秒、10^11 で約 1.5 秒。N ≤ 10^13） O(N^{2/3})
//    count_coprime_pairs(N)           1 ≤ a, b ≤ N で gcd = 1 の組の数（順序あり。Int128）
//    sum_gcd_tuples(K, len, m)        値が 1〜K で長さ len の列すべての gcd の和 mod m（ABC162 E。K ≤ 10^13） O(K^{2/3} + √K log len)
//  （φ・μ・約数の個数の表は 02_sieve：このファイルに入っている）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    divt_test::random_test()  … 全部調べる愚直解・1 つ前との差（φ(N)・μ(N)）との突き合わせ
//    divt_test::benchmark()    … 速度計測
//    divt_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_DIVISOR_TRANSFORM
#define GG_DIVISOR_TRANSFORM
// ---- 約数・倍数の変換（数学第二世代）  約数・倍数のゼータ・メビウス変換、gcd・lcm の畳み込み、gcd ごとの組の数、Σφ・Σμ、ディリクレ積、M の約数の上の変換 ----
// 配列は 1〜N を使う（f[0] は使わない・変えない）。T は mint・long long（溢れたら止める）・double など
namespace math_internal {
template <class T>
inline T dt_add(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail("divisor_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a + b;
    }
}
template <class T>
inline T dt_sub(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_sub_overflow(a, b, &r)) fail("divisor_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a - b;
    }
}
template <class T>
inline T dt_mul(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_mul_overflow(a, b, &r)) fail("divisor_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a * b;
    }
}
struct DtAdd {
    template <class T>
    T operator()(const T& a, const T& b) const { return dt_add(a, b); }
};
struct DtSub {
    template <class T>
    T operator()(const T& a, const T& b) const { return dt_sub(a, b); }
};
inline const std::vector<int>& dt_primes(int N) {  // N 以下の素数（前より大きい N が来たら作り直す）
    static std::vector<int> ps;
    static int built = -1;
    if (built < N) ps = primes_up_to<int>(std::max(N, 2)), built = std::max(N, 2);
    return ps;
}
}  // namespace math_internal

// g[n] = Σ_{d | n} f[d]（n = 1〜N。op を渡すと Σ の代わり）  O(N log log N)
template <class T, class Op = math_internal::DtAdd>
void divisor_zeta(std::vector<T>& f, Op op = Op{}) {
    const int N = (int)f.size() - 1;
    for (int p : math_internal::dt_primes(N)) {
        if (p > N) break;
        for (int i = 1; (long long)i * p <= N; ++i) f[(size_t)i * p] = op(f[(size_t)i * p], f[(size_t)i]);
    }
}
// divisor_zeta の逆（f[n] = Σ_{d | n} μ(n/d) g[d]）
template <class T, class Op = math_internal::DtSub>
void divisor_mobius(std::vector<T>& f, Op inv_op = Op{}) {
    const int N = (int)f.size() - 1;
    for (int p : math_internal::dt_primes(N)) {
        if (p > N) break;
        for (int i = N / p; i >= 1; --i) f[(size_t)i * p] = inv_op(f[(size_t)i * p], f[(size_t)i]);
    }
}
// g[d] = Σ_{d | n, n ≤ N} f[n]（倍数の和）
template <class T, class Op = math_internal::DtAdd>
void multiple_zeta(std::vector<T>& f, Op op = Op{}) {
    const int N = (int)f.size() - 1;
    for (int p : math_internal::dt_primes(N)) {
        if (p > N) break;
        for (int i = N / p; i >= 1; --i) f[(size_t)i] = op(f[(size_t)i], f[(size_t)i * p]);
    }
}
// multiple_zeta の逆
template <class T, class Op = math_internal::DtSub>
void multiple_mobius(std::vector<T>& f, Op inv_op = Op{}) {
    const int N = (int)f.size() - 1;
    for (int p : math_internal::dt_primes(N)) {
        if (p > N) break;
        for (int i = 1; (long long)i * p <= N; ++i) f[(size_t)i] = inv_op(f[(size_t)i], f[(size_t)i * p]);
    }
}
// c[k] = Σ_{gcd(i, j) = k} a[i] b[j]（i, j, k = 1〜N。長さは長い方にそろえる）  O(N log log N)
template <class T>
std::vector<T> gcd_convolution(std::vector<T> a, std::vector<T> b) {
    const size_t n = std::max(a.size(), b.size());
    a.resize(n, T(0)), b.resize(n, T(0));
    multiple_zeta(a), multiple_zeta(b);
    for (size_t i = 1; i < n; ++i) a[i] = math_internal::dt_mul(a[i], b[i]);
    multiple_mobius(a);
    if (n) a[0] = T(0);
    return a;
}
// c[k] = Σ_{lcm(i, j) = k} a[i] b[j]（k ≤ N のものだけ）
template <class T>
std::vector<T> lcm_convolution(std::vector<T> a, std::vector<T> b) {
    const size_t n = std::max(a.size(), b.size());
    a.resize(n, T(0)), b.resize(n, T(0));
    divisor_zeta(a), divisor_zeta(b);
    for (size_t i = 1; i < n; ++i) a[i] = math_internal::dt_mul(a[i], b[i]);
    divisor_mobius(a);
    if (n) a[0] = T(0);
    return a;
}
// cnt[v] = 値 v（1〜N）の個数のとき、gcd がちょうど g になる組 {i < j} の数（g = 1〜N）  O(N log log N)
inline std::vector<long long> count_pairs_by_gcd(const std::vector<long long>& cnt) {
    std::vector<long long> m(cnt);
    if (!m.empty()) m[0] = 0;
    multiple_zeta(m);  // m[d] = d の倍数の個数
    for (size_t d = 1; d < m.size(); ++d) {  // C(m, 2)（偶数の方を先に 2 で割る。溢れたら止める）
        const long long x = m[d], a = x % 2 == 0 ? x / 2 : x, b = x % 2 == 0 ? x - 1 : (x - 1) / 2;
        if (__builtin_mul_overflow(a, b, &m[d])) math_internal::fail("count_pairs_by_gcd", "組の数が long long に収まらない");
    }
    multiple_mobius(m);
    return m;
}
// ディリクレ積 c[n] = Σ_{d | n} a[d] b[n/d]（n = 1〜N。長さは長い方にそろえる）  O(N log N)
template <class T>
std::vector<T> dirichlet_convolution(const std::vector<T>& a, const std::vector<T>& b) {
    const size_t n = std::max(a.size(), b.size());
    std::vector<T> c(n, T(0));
    for (size_t i = 1; i < a.size(); ++i)
        for (size_t j = 1; j < b.size() && i * j < n; ++j) c[i * j] = math_internal::dt_add(c[i * j], math_internal::dt_mul(a[i], b[j]));
    return c;
}
// ---- M の約数の上での変換（M ≤ 10^18。約数の数 D ≤ 103680）----
// DivisorsOf D(M)：D.divisors() は M の約数（昇順）。配列 f は約数の番号（D.index(d)）で引く（長さ D.size()）
struct DivisorsOf {
  private:
    long long M_;
    std::vector<PrimePower> fs_;
    std::vector<long long> ds_;
    std::vector<std::vector<std::pair<int, int>>> steps_;  // 素数ごとに (i, j)：ds_[j] = ds_[i] · p（i の昇順）
    template <class T>
    void check_size(const std::vector<T>& f, const char* func) const {
        if (f.size() != ds_.size()) math_internal::fail(func, "配列の長さを約数の数 size() にする");
    }

  public:
    explicit DivisorsOf(long long M) : M_(M) {
        math_internal::check_positive(M, "DivisorsOf");
        fs_ = ::factorize(M), ds_ = divisors_from_factors(fs_);
        for (const PrimePower& f : fs_) {
            std::vector<std::pair<int, int>> st;
            size_t j = 0;
            for (size_t i = 0; i < ds_.size(); ++i) {
                if ((M / ds_[i]) % f.p != 0) continue;  // ds_[i] · p が M を割り切るものだけ
                const long long t = ds_[i] * f.p;
                while (ds_[j] < t) ++j;
                st.push_back({(int)i, (int)j});
            }
            steps_.push_back(std::move(st));
        }
    }
    long long value() const { return M_; }
    int size() const { return (int)ds_.size(); }                          // 約数の数
    const std::vector<long long>& divisors() const { return ds_; }       // 約数（昇順）
    const std::vector<PrimePower>& factors() const { return fs_; }       // M の素因数分解
    long long operator[](int i) const { return ds_[(size_t)i]; }          // i 番目の約数
    int index(long long d) const {  // 約数 d の番号（d が M の約数でなければ -1）。二分探索（std::lower_bound はデバッグビルドで O(D) になるので自分で）
        size_t lo = 0, hi = ds_.size();
        while (lo < hi) {
            const size_t mid = (lo + hi) / 2;
            if (ds_[mid] < d) lo = mid + 1;
            else hi = mid;
        }
        return lo < ds_.size() && ds_[lo] == d ? (int)lo : -1;
    }
    // g[d] = Σ_{e | d} f[e]（d | M）  O(D · ω(M))
    template <class T, class Op = math_internal::DtAdd>
    void zeta(std::vector<T>& f, Op op = Op{}) const {
        check_size(f, "DivisorsOf::zeta");
        for (const auto& st : steps_)
            for (const auto& ij : st) f[(size_t)ij.second] = op(f[(size_t)ij.second], f[(size_t)ij.first]);
    }
    template <class T, class Op = math_internal::DtSub>
    void mobius(std::vector<T>& f, Op inv_op = Op{}) const {  // zeta の逆
        check_size(f, "DivisorsOf::mobius");
        for (const auto& st : steps_)
            for (auto it = st.rbegin(); it != st.rend(); ++it) f[(size_t)it->second] = inv_op(f[(size_t)it->second], f[(size_t)it->first]);
    }
    // g[d] = Σ_{d | e | M} f[e]（倍数の和）
    template <class T, class Op = math_internal::DtAdd>
    void multiple_zeta(std::vector<T>& f, Op op = Op{}) const {
        check_size(f, "DivisorsOf::multiple_zeta");
        for (const auto& st : steps_)
            for (auto it = st.rbegin(); it != st.rend(); ++it) f[(size_t)it->first] = op(f[(size_t)it->first], f[(size_t)it->second]);
    }
    template <class T, class Op = math_internal::DtSub>
    void multiple_mobius(std::vector<T>& f, Op inv_op = Op{}) const {  // multiple_zeta の逆
        check_size(f, "DivisorsOf::multiple_mobius");
        for (const auto& st : steps_)
            for (const auto& ij : st) f[(size_t)ij.first] = inv_op(f[(size_t)ij.first], f[(size_t)ij.second]);
    }
};
// Σ_{i=1}^{N} φ(i)・Σ μ(i)（N ≤ 10^11 程度。杜教篩：小さいところは表、大きいところは ⌊N/k⌋ ごとに）  O(N^{2/3})
namespace math_internal {
struct DtPrefix {  // ⌊N/i⌋ の形の v での和 S(v)：v ≤ T は表 small、それより大きい v = ⌊N/k⌋ は big[k]
    long long N = 0, T = 0;
    std::vector<long long> small;
    std::vector<Int128> big;
    Int128 operator()(long long v) const { return v <= T ? (Int128)small[(size_t)v] : big[(size_t)(N / v)]; }
};
template <class Base>
void dt_du_sieve(DtPrefix& P, Base base) {  // S(x) = base(x) − Σ_{d=2}^{x} S(⌊x/d⌋)
    const long long N = P.N, T = P.T;
    if (N <= T) return;
    const long long K = N / T;
    P.big.assign((size_t)K + 1, 0);
    for (long long k = K; k >= 1; --k) {
        const long long x = N / k;
        Int128 r = base(x);
        for (long long d = 2; d <= x;) {
            const long long q = x / d, d2 = x / q;
            const Int128 v = q <= T ? (Int128)P.small[(size_t)q] : P.big[(size_t)(k * d)];
            r -= (Int128)(d2 - d + 1) * v;
            d = d2 + 1;
        }
        P.big[(size_t)k] = r;
    }
}
inline long long dt_threshold(long long N) {  // 表を作る大きさ（N^{2/3}。10^7 まで）
    long long t = (long long)std::cbrt((long double)N);
    t = t * t;
    return std::max<long long>(1, std::min<long long>({N, t, 10000000LL}));
}
inline DtPrefix dt_phi_prefix(long long N, const char* func) {  // Σφ を ⌊N/i⌋ の形の値すべてで
    if (N > 10000000000000LL) fail(func, "N ≤ 10^13 で使う（表の大きさ。10^11 で約 1.5 秒）");
    DtPrefix P;
    P.N = std::max(N, 0LL), P.T = dt_threshold(std::max(N, 1LL));
    P.small = phi_table((int)P.T);
    for (size_t i = 1; i < P.small.size(); ++i) P.small[i] += P.small[i - 1];
    dt_du_sieve(P, [](long long x) { return (Int128)x * (x + 1) / 2; });
    return P;
}
}  // namespace math_internal
inline Int128 phi_sum(long long N) {
    if (N <= 0) return 0;
    return math_internal::dt_phi_prefix(N, "phi_sum")(N);
}
// Σ_{i=1}^{N} μ(i)（メルテンス関数。N ≤ 10^11 程度）  O(N^{2/3})
inline long long mertens(long long N) {
    if (N <= 0) return 0;
    if (N > 10000000000000LL) math_internal::fail("mertens", "N ≤ 10^13 で使う（表の大きさ。10^11 で約 1.5 秒）");
    math_internal::DtPrefix P;
    P.N = N, P.T = math_internal::dt_threshold(N);
    const std::vector<int> mu = mu_table((int)P.T);
    P.small.assign(mu.size(), 0);
    for (size_t i = 1; i < P.small.size(); ++i) P.small[i] = P.small[i - 1] + mu[i];
    math_internal::dt_du_sieve(P, [](long long) { return (Int128)1; });
    return (long long)P(N);
}
// 1 ≤ a, b ≤ N で gcd(a, b) = 1 となる組 (a, b) の数（順序あり）= 2 Σφ − 1
inline Int128 count_coprime_pairs(long long N) { return N <= 0 ? 0 : 2 * phi_sum(N) - 1; }
// 長さ len、値が 1〜K の列すべての gcd の和 mod m（ABC162 E）= Σ_{d ≤ K} φ(d) ⌊K/d⌋^len（⌊K/d⌋ が同じ d をまとめる。len = 0 なら 0）
//   K ≤ 10^13 程度  O(K^{2/3} + √K log len)
inline long long sum_gcd_tuples(long long K, long long len, long long m) {
    if (K < 0 || len < 0 || m < 1) math_internal::fail("sum_gcd_tuples", "K, len ≥ 0、m ≥ 1 で使う");
    if (K == 0 || len == 0) return 0;  // 列がない / 空の列の gcd は 0
    const math_internal::DtPrefix P = math_internal::dt_phi_prefix(K, "sum_gcd_tuples");
    long long s = 0;
    Int128 prev = 0;  // Φ(d − 1)
    for (long long d = 1; d <= K;) {
        const long long q = K / d, d2 = K / q;
        const Int128 cur = P(d2);
        const long long blk = (long long)((cur - prev) % m);
        s = (long long)(((Int128)s + (Int128)blk * mod_pow(q, len, m)) % m);
        prev = cur, d = d2 + 1;
    }
    return s;
}
#endif  // GG_DIVISOR_TRANSFORM
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace divt_test {

#define DIVT_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[divisor_transform] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static long long bgcd(long long a, long long b) { return b ? bgcd(b, a % b) : a; }

static void test_transforms(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int N = (int)rnd(rng, 0, it % 2 ? 40 : 300);
        std::vector<long long> f((size_t)N + 1), g((size_t)N + 1);
        for (auto& x : f) x = rnd(rng, -1000, 1000);
        for (auto& x : g) x = rnd(rng, -1000, 1000);
        std::vector<long long> dz((size_t)N + 1, 0), mz((size_t)N + 1, 0);
        dz[0] = f[0], mz[0] = f[0];
        for (int n = 1; n <= N; ++n)
            for (int d = 1; d <= N; ++d) {
                if (n % d == 0) dz[(size_t)n] += f[(size_t)d];
                if (d % n == 0) mz[(size_t)n] += f[(size_t)d];
            }
        std::vector<long long> a = f, b = f;
        divisor_zeta(a), multiple_zeta(b);
        DIVT_CHECK(a == dz && b == mz);
        divisor_mobius(a), multiple_mobius(b);
        DIVT_CHECK(a == f && b == f);
        // gcd・lcm の畳み込み
        if (N <= 80) {
            std::vector<long long> cg((size_t)N + 1, 0), cl((size_t)N + 1, 0);
            for (int i = 1; i <= N; ++i)
                for (int j = 1; j <= N; ++j) {
                    const long long gg = bgcd(i, j), l = (long long)i / gg * j;
                    cg[(size_t)gg] += f[(size_t)i] * g[(size_t)j];
                    if (l <= N) cl[(size_t)l] += f[(size_t)i] * g[(size_t)j];
                }
            DIVT_CHECK(gcd_convolution(f, g) == cg && lcm_convolution(f, g) == cl);
            // gcd ごとの組の数
            std::vector<long long> cnt((size_t)N + 1, 0), e((size_t)N + 1, 0);
            std::vector<long long> vals;
            for (int t = 0; t < (int)rnd(rng, 0, 30) && N >= 1; ++t) {
                const long long v = rnd(rng, 1, N);
                cnt[(size_t)v]++, vals.push_back(v);
            }
            for (size_t i = 0; i < vals.size(); ++i)
                for (size_t j = i + 1; j < vals.size(); ++j) e[(size_t)bgcd(vals[i], vals[j])]++;
            DIVT_CHECK(count_pairs_by_gcd(cnt) == e);
        }
    }
    // 演算を渡す：約数の最大値
    std::vector<long long> f = {0, 5, 1, 7, 2, 3, 9};
    divisor_zeta(f, [](long long x, long long y) { return std::max(x, y); });
    DIVT_CHECK(f == std::vector<long long>({0, 5, 5, 7, 5, 5, 9}));
    // ディリクレ積
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int na = (int)rnd(rng, 0, 120), nb = it % 3 ? na : (int)rnd(rng, 0, 120);
        std::vector<long long> a((size_t)na), b((size_t)nb);
        for (auto& x : a) x = rnd(rng, -1000, 1000);
        for (auto& x : b) x = rnd(rng, -1000, 1000);
        const size_t n = std::max(a.size(), b.size());
        std::vector<long long> c(n, 0);
        for (size_t k = 1; k < n; ++k)
            for (size_t d = 1; d <= k; ++d)
                if (k % d == 0 && d < a.size() && k / d < b.size()) c[k] += a[d] * b[k / d];
        DIVT_CHECK(dirichlet_convolution(a, b) == c);
    }
    {  // μ * 1 = ε、φ * 1 = id
        const int N = 1000;
        const std::vector<int> mu = mu_table(N);
        const std::vector<long long> phi = phi_table(N);
        const std::vector<long long> one((size_t)N + 1, 1), mul(mu.begin(), mu.end());
        const std::vector<long long> e = dirichlet_convolution(mul, one), id = dirichlet_convolution(phi, one);
        for (int k = 1; k <= N; ++k) DIVT_CHECK(e[(size_t)k] == (k == 1) && id[(size_t)k] == k);
    }
    // M の約数の上での変換
    std::vector<long long> Ms = {1, 2, 12, 1024, 720720, 735134400LL, 963761198400LL, 1000000007LL, 999999999999999989LL, 897612484786617600LL};
    for (int it = 0; it < 300 * TEST_SCALE; ++it) Ms.push_back(rnd(rng, 1, it % 3 == 0 ? 1000000000000000000LL : (it % 3 == 1 ? 1000000 : 1000)));
    for (const long long M : Ms) {
        const DivisorsOf D(M);
        const std::vector<long long> ds = divisors_of(M);
        const int K = D.size();
        DIVT_CHECK(D.divisors() == ds && K == (int)ds.size() && D.value() == M && D.factors() == factorize(M));
        for (int i = 0; i < K; ++i) DIVT_CHECK(D.index(ds[(size_t)i]) == i && D[i] == ds[(size_t)i]);
        DIVT_CHECK(D.index(M + 1) == -1 && (M < 3 || M % 3 == 0 || D.index(3) == -1));
        std::vector<long long> fv((size_t)K);
        for (auto& x : fv) x = rnd(rng, -1000000, 1000000);
        std::vector<long long> a = fv, b = fv;
        D.zeta(a), D.multiple_zeta(b);
        if (K <= 1500) {
            for (int i = 0; i < K; ++i, ++cases) {
                long long sa = 0, sb = 0;
                for (int j = 0; j < K; ++j) {
                    if (ds[(size_t)i] % ds[(size_t)j] == 0) sa += fv[(size_t)j];
                    if (ds[(size_t)j] % ds[(size_t)i] == 0) sb += fv[(size_t)j];
                }
                DIVT_CHECK(a[(size_t)i] == sa && b[(size_t)i] == sb);
            }
        } else {  // 大きいときは全部 1 の変換 = 約数の個数 / 倍数の個数
            std::vector<long long> o((size_t)K, 1), o2((size_t)K, 1);
            D.zeta(o), D.multiple_zeta(o2);
            for (int i = 0; i < K; i += 37, ++cases) DIVT_CHECK(o[(size_t)i] == count_divisors(ds[(size_t)i]) && o2[(size_t)i] == count_divisors(M / ds[(size_t)i]));
        }
        D.mobius(a), D.multiple_mobius(b);
        DIVT_CHECK(a == fv && b == fv);
    }
    {  // 演算を渡す：倍数の最小値
        const DivisorsOf D(12);  // 1 2 3 4 6 12
        std::vector<long long> g = {5, 3, 9, 8, 1, 7};
        D.multiple_zeta(g, [](long long x, long long y) { return std::min(x, y); });
        DIVT_CHECK(g == std::vector<long long>({1, 1, 1, 7, 1, 7}));
    }
}

static void test_sums(std::mt19937_64& rng, int& cases) {
    const int M = 2000000;
    const std::vector<long long> phi = phi_table(M);
    const std::vector<int> mu = mu_table(M);
    std::vector<Int128> ps((size_t)M + 1, 0);
    std::vector<long long> ms((size_t)M + 1, 0);
    for (int i = 1; i <= M; ++i) ps[(size_t)i] = ps[(size_t)i - 1] + phi[(size_t)i], ms[(size_t)i] = ms[(size_t)i - 1] + mu[(size_t)i];
    for (int it = 0; it < 60 * TEST_SCALE; ++it, ++cases) {
        const long long N = it < 20 ? it : rnd(rng, 1, M);
        DIVT_CHECK(phi_sum(N) == ps[(size_t)N] && mertens(N) == ms[(size_t)N] && count_coprime_pairs(N) == (N ? 2 * ps[(size_t)N] - 1 : 0));
    }
    for (int it = 0; it < 4; ++it, ++cases) {  // 大きい N：1 つ前との差は φ(N)・μ(N)
        const long long N = rnd(rng, 3000000000LL, it % 2 ? 10000000000LL : 3000000000LL + 100);
        DIVT_CHECK(phi_sum(N) - phi_sum(N - 1) == euler_phi(N) && mertens(N) - mertens(N - 1) == mobius(N));
    }
    // gcd の和：小さい場合を全部数える
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const long long K = rnd(rng, 1, 6), len = rnd(rng, 1, 5), m = rnd(rng, 1, 1000000007);
        long long total = 0;
        std::vector<long long> v((size_t)len, 1);
        while (true) {
            long long gg = 0;
            for (long long x : v) gg = bgcd(gg, x);
            total += gg;
            size_t i = 0;
            while (i < v.size() && v[i] == K) v[i] = 1, ++i;
            if (i == v.size()) break;
            v[i]++;
        }
        DIVT_CHECK(sum_gcd_tuples(K, len, m) == total % m);
    }
    // 大きい K：φ の表で Σ φ(d) ⌊K/d⌋^len を直接
    {
        const std::vector<long long> ph = phi_table(3000000);
        for (int it = 0; it < 6; ++it, ++cases) {
            const long long K = it == 0 ? 3000000 : rnd(rng, 1, 300000), len = it == 0 ? 3 : rnd(rng, 1, 1000000000000000000LL), m = rnd(rng, 1, 1000000000000LL);
            long long e = 0;
            for (long long d = 1; d <= K; ++d) e = (long long)(((Int128)e + (Int128)(ph[(size_t)d] % m) * mod_pow(K / d, len, m)) % m);
            DIVT_CHECK(sum_gcd_tuples(K, len, m) == e);
        }
    }
    DIVT_CHECK(sum_gcd_tuples(5, 0, 7) == 0 && sum_gcd_tuples(0, 3, 7) == 0);
    DIVT_CHECK(count_pairs_by_gcd({0, 3000000000LL})[1] == 4499999998500000000LL);
    cases += 2;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_transforms(rng, cases);
    test_sums(rng, cases);
    std::cout << "[divisor_transform] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[divisor_transform] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 10000000;
    std::vector<long long> f((size_t)N + 1);
    for (auto& x : f) x = rnd(rng, 0, 1000);
    auto t0 = std::chrono::steady_clock::now();
    divisor_zeta(f), divisor_mobius(f), multiple_zeta(f), multiple_mobius(f);
    const long long t_tr = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Int128 s1 = phi_sum(100000000000LL);
    const long long t_phi = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const long long s2 = mertens(100000000000LL);
    const long long t_mu = ms(t0);
    std::cout << "[divisor_transform] benchmark: 約数・倍数のゼータ・メビウス 4 回（N = 10^7）" << t_tr << " ms, phi_sum(10^11) " << t_phi << " ms, mertens(10^11) " << t_mu
              << " ms  (check " << (long long)((s1 + s2 + f[12345]) % 7 + 7) % 7 << ")\n";
}

}  // namespace divt_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace divt_demo {

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
#define SHOW(expr, comment) divt_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    divt_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("約数・倍数の変換（divisor_transform） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  divisor_zeta / divisor_mobius / multiple_zeta / multiple_mobius（f[1..N] を書きかえる。op を渡せる）\n"
              << "  gcd_convolution / lcm_convolution / count_pairs_by_gcd / dirichlet_convolution\n"
              << "  phi_sum / mertens / count_coprime_pairs / sum_gcd_tuples\n"
              << "  DivisorsOf D(M)：M の約数の上で zeta / mobius / multiple_zeta / multiple_mobius（M ≤ 10^18）\n";

    api("【1】 divisor_zeta(f) / divisor_mobius(f) / multiple_zeta(f) / multiple_mobius(f)", "なし（f を書きかえる）",
        "約数の和 g[n] = Σ_{d|n} f[d] / 倍数の和 g[d] = Σ_{d|n} f[n] / それぞれの逆（f[0] は使わない）", "O(N log log N)");
    RUN(std::vector<long long> f = {0, 1, 1, 1, 1, 1, 1}; divisor_zeta(f));
    SHOW(f, "1〜6 の約数の個数");
    RUN(std::vector<long long> g = {0, 1, 1, 1, 1, 1, 1}; multiple_zeta(g));
    SHOW(g, "6 以下の d の倍数の個数 = ⌊6/d⌋");
    note("「gcd がちょうど g」の数は、まず「g で割り切れる」で数えて multiple_mobius（包除）で戻す");

    api("【2】 gcd_convolution(a, b) / lcm_convolution(a, b) / count_pairs_by_gcd(cnt)", "vector",
        "c[k] = Σ_{gcd(i,j)=k} a[i] b[j] / lcm の版 / gcd がちょうど g の組 {i < j} の数", "O(N log log N)");
    SHOW(gcd_convolution<long long>({0, 1, 1, 1, 1}, {0, 1, 1, 1, 1}), "1〜4 から 2 つ（順序あり）選んだ gcd ごとの数");
    SHOW(count_pairs_by_gcd({0, 0, 2, 1, 1}), "{2, 2, 3, 4} の 2 つ組の gcd：gcd 1 が 3 組、gcd 2 が 3 組");
    {
        // ABC206 E：[L, R] の組 (x, y) で gcd = g ≥ 2 かつ x ≠ g かつ y ≠ g の数
        auto solve = [](long long L, long long R) {
            std::vector<long long> F((size_t)R + 1, 0);
            for (long long d = 1; d <= R; ++d) {
                const long long c = R / d - (L - 1) / d;  // [L, R] の d の倍数
                F[(size_t)d] = c * c;
            }
            multiple_mobius(F);  // gcd がちょうど d の組の数
            long long ans = 0;
            for (long long d = 2; d <= R; ++d) {
                ans += F[(size_t)d];
                if (d >= L) ans -= 2 * (R / d - (L - 1) / d) - 1;  // x = d または y = d
            }
            return ans;
        };
        note("ABC206 E（Divide Both）: 「g の倍数どうしの組」c_g^2 を multiple_mobius して gcd ちょうど g の数に");
        std::cout << "    (3, 7) → " << solve(3, 7) << "、(4, 10) → " << solve(4, 10) << "、(1, 10^6) → " << solve(1, 1000000) << "\n";
    }

    api("【3】 phi_sum(N) / mertens(N) / count_coprime_pairs(N) / sum_gcd_tuples(K, len, m)", "Int128 / long long",
        "Σφ（Int128）/ Σμ（N ≤ 10^11 程度）/ 互いに素な組の数 / 列の gcd の和 mod m", "O(N^{2/3}) / O(K log len)");
    SHOW(phi_sum(10), "1 + 1 + 2 + 2 + 4 + 2 + 6 + 4 + 6 + 4");
    SHOW(phi_sum(10000000000LL), "N = 10^10");
    SHOW(mertens(1000000000LL), "Σμ(i)（i ≤ 10^9）");
    SHOW(count_coprime_pairs(3), "(1,1) (1,2) (2,1) (1,3) (3,1) (2,3) (3,2)");
    note("ABC162 E（Sum of gcd of Tuples (Hard)）: sum_gcd_tuples(K, N, 10^9+7)");
    std::cout << "    (N, K) = (3, 2) → " << sum_gcd_tuples(2, 3, 1000000007) << "、(3, 200) → " << sum_gcd_tuples(200, 3, 1000000007) << "、(10^5, 10^5) → "
              << sum_gcd_tuples(100000, 100000, 1000000007) << "\n";

    api("【4】 dirichlet_convolution(a, b)", "vector", "ディリクレ積 c[n] = Σ_{d|n} a[d] b[n/d]（n = 1〜N。乗法的関数の組み合わせ）", "O(N log N)");
    SHOW(dirichlet_convolution<long long>({0, 1, 1, 1, 1, 1, 1}, {0, 1, 1, 1, 1, 1, 1}), "1 * 1 = 約数の個数");
    SHOW(dirichlet_convolution<long long>({0, 1, 2, 3, 4, 5, 6}, {0, 1, -1, -1, 0, -1, 1}), "id * μ = φ");

    api("【5】 DivisorsOf D(M);  D.zeta(f) / D.mobius(f) / D.multiple_zeta(f) / D.multiple_mobius(f)", "なし（f を書きかえる）",
        "M（≤ 10^18）の約数だけの上の約数の和・倍数の和とその逆。f は約数の番号で引く（長さ D.size()、番号は D.index(d)、値は D[i]）",
        "O(D · ω(M))（D は約数の数）");
    RUN(const DivisorsOf D(12));
    SHOW(D.divisors(), "約数（昇順）");
    SHOW(D.index(6), "6 は 4 番目（0 から）");
    {
        // 1〜M の x で gcd(x, M) = d になる数：まず「d の倍数」M/d 個を数えて multiple_mobius
        const long long M = 1000000000000LL;
        const DivisorsOf E(M);
        std::vector<long long> c((size_t)E.size());
        for (int i = 0; i < E.size(); ++i) c[(size_t)i] = M / E[i];
        E.multiple_mobius(c);  // c[i] = gcd(x, M) = E[i] になる x の数（= φ(M / E[i])）
        Int128 s = 0;
        for (int i = 0; i < E.size(); ++i) s += (Int128)E[i] * c[(size_t)i];
        note("Σ_{x=1}^{M} gcd(x, M)：gcd ごとの個数を「倍数の個数 M/d」の multiple_mobius で");
        std::cout << "    M = 10^12（約数 " << E.size() << " 個）→ " << s << "\n";
    }
    {
        // M の約数 a, b で lcm(a, b) = M の組 (a, b) の数：「d の約数どうし」(約数の数)^2 を mobius
        const long long M = 720720;
        const DivisorsOf E(M);
        std::vector<long long> cnt((size_t)E.size(), 1);
        E.zeta(cnt);  // cnt[i] = E[i] の約数の個数
        for (auto& x : cnt) x *= x;
        E.mobius(cnt);  // cnt[i] = lcm がちょうど E[i] の組の数
        note("lcm(a, b) = M となる約数の組 (a, b)：「どちらも d の約数」の数を mobius で「lcm がちょうど d」に");
        std::cout << "    M = 720720 → " << cnt.back() << "\n";
    }
}

}  // namespace divt_demo

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
    divt_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    divt_test::benchmark();    // 速度計測
    divt_demo::demo();         // 使い方デモ
    return 0;
}
