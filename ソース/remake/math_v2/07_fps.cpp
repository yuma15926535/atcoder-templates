// =============================================================================
//  形式的冪級数（Fps<mint>：逆数・log・exp・累乗・平方根・多項式の割り算・テイラーシフト・P/Q の K 次の係数・多点評価・補間・べき乗和）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 07 畳み込み・多項式
//  依存       : 00_int_common, 00_math_base, 00_modint, 07_convolution（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  使い方
//      using mint = ModInt998;               // 998244353 がいちばん速い（ほかの mod でも 3 つの素数の畳み込みで動く）
//      Fps<mint> f = {1, 1};                 // 係数は低い次数から（1 + x）
//      Fps<mint> g = fps_inv(f, 10);         // 1/f の最初の 10 項（x^0〜x^9）
//
//  機能一覧（n を渡すものは最初の n 項を返す）                                    計算量
//    Fps<mint>                        vector<mint> と同じ。+ − ×（畳み込み）・スカラー倍・単項 −、f.pre(n)（最初の n 項）、f.eval(x)、
//                                     cout << f（係数を空白区切りで）
//    fps_inv(f, n)                    1/f（f[0] ≠ 0）                                                      O(n log n)
//    fps_log(f, n) / fps_exp(f, n)    log f（f[0] = 1）/ exp f（f[0] = 0）
//    fps_pow(f, k, n)                 f^k（k ≤ 10^18。f[0] = 0 でもよい）
//    fps_sqrt(f, n)                   √f（Maybe：ok が false なら平方根なし）
//    fps_divmod(f, g)                 多項式の割り算 FpsDivMod{q, r}（f = q g + r）
//    fps_derivative(f) / fps_integral(f)    微分 / 積分（定数項 0）                                      O(n)
//    taylor_shift(f, c)               f(x + c) の係数                                                       O(n log n)
//    fps_coeff_of_fraction(P, Q, K)   [x^K] P/Q（K ≤ 10^18。Bostan–Mori。線形漸化式の第 K 項）            O(d log d log K)
//    prod_linear(as)                  Π (1 + a_i x)                                                         O(n log^2 n)
//    multipoint_eval(f, xs)           f(x_0), …, f(x_{m−1})（多点評価。10^5 点で約 0.4 秒）                  O(n log n + m log^2 m)
//    interpolate(xs, ys)              点 (x_i, y_i)（x は異なる）を通る次数 n − 1 以下の多項式                 O(n log^2 n)
//    power_sums(as, n, ws)            べき乗和 p_k = Σ w_i a_i^k（k = 0〜n−1。ws を省くと 1）                 O(m log^2 m + n log n)
//    fps_sum_of_fractions(nums, dens) Σ nums[i] / dens[i] を 1 つの分数 FpsFraction{num, den} に（通分だけ）  O(N log^2 N)
//  log・exp・pow・integral・taylor_shift は mod が項の数より大きい素数のとき（998244353 など）、sqrt は mod が奇素数のとき
//  （漸化式の推測と第 K 項は 06_linear_recurrence、次数の小さい多項式の補間は 06_series）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    fps_test::random_test()  … 1 項ずつ決める愚直な漸化式との突き合わせ
//    fps_test::benchmark()    … 速度計測
//    fps_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#include <queue>
#include <random>
#include <set>
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
#ifndef GG_CONVOLUTION
#define GG_CONVOLUTION
// ---- 畳み込み（数学第二世代）  mint の列（NTT・任意の mod）、整数のまま（long long・Int128）、多くの多項式の積、相関、巡回、オンライン ----
namespace math_internal {
template <class mint>
struct NttTable {  // mod ごとの 1 の冪根の表（ModInt は 1 回、DynModInt は mod が変わったら作り直す）
    int mod = 0, rank2 = 0;
    std::vector<mint> root, iroot, rate2, irate2, rate3, irate3;
};
inline unsigned primitive_root_u32(unsigned p) {  // 素数 p の原始根（いちばん小さいもの）
    if (p == 2) return 1;
    std::vector<unsigned> qs;
    unsigned x = p - 1;
    for (unsigned q = 2; (unsigned long long)q * q <= x; ++q)
        if (x % q == 0) {
            qs.push_back(q);
            while (x % q == 0) x /= q;
        }
    if (x > 1) qs.push_back(x);
    for (unsigned g = 2;; ++g) {
        bool ok = true;
        for (unsigned q : qs) {
            unsigned long long r = 1, b = g, e = (p - 1) / q;
            for (; e; e >>= 1, b = b * b % p)
                if (e & 1) r = r * b % p;
            if (r == 1) { ok = false; break; }
        }
        if (ok) return g;
    }
}
template <class mint>
const NttTable<mint>& ntt_table() {
    static NttTable<mint> t;
    if (t.mod == mint::mod()) return t;
    const unsigned p = (unsigned)mint::mod();
    if (!is_prime_u32(p)) fail("ntt", "mod が素数でない（convolve なら自動で 3 つの素数を使う）");
    t = NttTable<mint>();
    t.mod = (int)p;
    t.rank2 = __builtin_ctz(p - 1);
    const int r = t.rank2;
    t.root.assign(r + 1, mint(1)), t.iroot.assign(r + 1, mint(1));
    t.rate2.assign(std::max(0, r - 1), mint(1)), t.irate2.assign(std::max(0, r - 1), mint(1));
    t.rate3.assign(std::max(0, r - 2), mint(1)), t.irate3.assign(std::max(0, r - 2), mint(1));
    const mint g = mint((long long)primitive_root_u32(p));
    t.root[r] = g.pow((long long)((p - 1) >> r));
    t.iroot[r] = t.root[r].inv();
    for (int i = r - 1; i >= 0; --i) t.root[i] = t.root[i + 1] * t.root[i + 1], t.iroot[i] = t.iroot[i + 1] * t.iroot[i + 1];
    mint prod = 1, iprod = 1;
    for (int i = 0; i + 2 <= r; ++i) {
        t.rate2[i] = t.root[i + 2] * prod, t.irate2[i] = t.iroot[i + 2] * iprod;
        prod *= t.iroot[i + 2], iprod *= t.root[i + 2];
    }
    prod = 1, iprod = 1;
    for (int i = 0; i + 3 <= r; ++i) {
        t.rate3[i] = t.root[i + 3] * prod, t.irate3[i] = t.iroot[i + 3] * iprod;
        prod *= t.iroot[i + 3], iprod *= t.root[i + 3];
    }
    return t;
}
template <class mint>
bool ntt_usable(size_t len) {  // mod が素数で、長さ len（2 の冪）の NTT ができるか
    const unsigned p = (unsigned)mint::mod();
    if (p < 3 || !is_prime_u32(p)) return false;
    return (size_t)1 << __builtin_ctz(p - 1) >= len;
}
inline int ceil_log2(size_t n) {
    int h = 0;
    while (((size_t)1 << h) < n) ++h;
    return h;
}
}  // namespace math_internal

// 長さ 2^h の NTT（結果はビットを逆にした順。intt と組で、点ごとの積をとるのに使う。mod は素数で 2^h | mod − 1）
template <class mint>
void ntt(std::vector<mint>& a) {
    const int n = (int)a.size(), h = math_internal::ceil_log2(a.size());
    if (n == 0) return;
    if ((1 << h) != n) math_internal::fail("ntt", "長さは 2 の冪");
    const auto& t = math_internal::ntt_table<mint>();
    if (h > t.rank2) math_internal::fail("ntt", "この mod では長さが大きすぎる");
    int len = 0;
    while (len < h) {
        if (h - len == 1) {
            const int p = 1 << (h - len - 1);
            mint rot = 1;
            for (int s = 0; s < (1 << len); ++s) {
                const int off = s << (h - len);
                for (int i = 0; i < p; ++i) {
                    const mint l = a[i + off], r = a[i + off + p] * rot;
                    a[i + off] = l + r, a[i + off + p] = l - r;
                }
                if (s + 1 != (1 << len)) rot *= t.rate2[__builtin_ctz(~(unsigned)s)];
            }
            ++len;
        } else {
            const int p = 1 << (h - len - 2);
            mint rot = 1;
            const mint imag = t.root[2];
            for (int s = 0; s < (1 << len); ++s) {
                const mint rot2 = rot * rot, rot3 = rot2 * rot;
                const int off = s << (h - len);
                for (int i = 0; i < p; ++i) {
                    const mint a0 = a[i + off], a1 = a[i + off + p] * rot, a2 = a[i + off + 2 * p] * rot2, a3 = a[i + off + 3 * p] * rot3;
                    const mint x = (a1 - a3) * imag;
                    a[i + off] = a0 + a2 + a1 + a3;
                    a[i + off + p] = a0 + a2 - a1 - a3;
                    a[i + off + 2 * p] = a0 - a2 + x;
                    a[i + off + 3 * p] = a0 - a2 - x;
                }
                if (s + 1 != (1 << len)) rot *= t.rate3[__builtin_ctz(~(unsigned)s)];
            }
            len += 2;
        }
    }
}
template <class mint>
void intt(std::vector<mint>& a, bool divide = true) {  // ntt の逆（divide なら 1/n 倍まで）
    const int n = (int)a.size(), h = math_internal::ceil_log2(a.size());
    if (n == 0) return;
    if ((1 << h) != n) math_internal::fail("intt", "長さは 2 の冪");
    const auto& t = math_internal::ntt_table<mint>();
    if (h > t.rank2) math_internal::fail("intt", "この mod では長さが大きすぎる");
    int len = h;
    while (len) {
        if (len == 1) {
            const int p = 1 << (h - len);
            mint irot = 1;
            for (int s = 0; s < (1 << (len - 1)); ++s) {
                const int off = s << (h - len + 1);
                for (int i = 0; i < p; ++i) {
                    const mint l = a[i + off], r = a[i + off + p];
                    a[i + off] = l + r, a[i + off + p] = (l - r) * irot;
                }
                if (s + 1 != (1 << (len - 1))) irot *= t.irate2[__builtin_ctz(~(unsigned)s)];
            }
            --len;
        } else {
            const int p = 1 << (h - len);
            mint irot = 1;
            const mint iimag = t.iroot[2];
            for (int s = 0; s < (1 << (len - 2)); ++s) {
                const mint irot2 = irot * irot, irot3 = irot2 * irot;
                const int off = s << (h - len + 2);
                for (int i = 0; i < p; ++i) {
                    const mint a0 = a[i + off], a1 = a[i + off + p], a2 = a[i + off + 2 * p], a3 = a[i + off + 3 * p];
                    const mint x = (a2 - a3) * iimag;
                    a[i + off] = a0 + a1 + a2 + a3;
                    a[i + off + p] = (a0 - a1 + x) * irot;
                    a[i + off + 2 * p] = (a0 + a1 - a2 - a3) * irot2;
                    a[i + off + 3 * p] = (a0 - a1 - x) * irot3;
                }
                if (s + 1 != (1 << (len - 2))) irot *= t.irate3[__builtin_ctz(~(unsigned)s)];
            }
            len -= 2;
        }
    }
    if (divide) {
        const mint inv_n = mint(n).inv();
        for (mint& x : a) x *= inv_n;
    }
}

namespace math_internal {
template <class mint>
std::vector<mint> conv_naive(const std::vector<mint>& a, const std::vector<mint>& b) {
    std::vector<mint> c(a.size() + b.size() - 1);
    if (a.size() < b.size()) {
        for (size_t j = 0; j < b.size(); ++j)
            for (size_t i = 0; i < a.size(); ++i) c[i + j] += a[i] * b[j];
    } else {
        for (size_t i = 0; i < a.size(); ++i)
            for (size_t j = 0; j < b.size(); ++j) c[i + j] += a[i] * b[j];
    }
    return c;
}
template <class mint>
std::vector<mint> conv_ntt(std::vector<mint> a, std::vector<mint> b) {  // mod が NTT に使える素数のとき
    const size_t z = a.size() + b.size() - 1, n = (size_t)1 << ceil_log2(z);
    const bool same = a == b;
    const mint inv_n = mint((long long)n).inv();  // intt の 1/n は点ごとの積に混ぜる
    a.resize(n), ntt(a);
    if (same) {
        for (size_t i = 0; i < n; ++i) a[i] *= a[i] * inv_n;
    } else {
        b.resize(n), ntt(b);
        for (size_t i = 0; i < n; ++i) a[i] *= b[i] * inv_n;
    }
    intt(a, false);
    a.resize(z);
    return a;
}
// 3〜5 個の NTT 素数（2^23 以上の 2 の冪で割れる）での畳み込みの余り
constexpr unsigned NTT_P0 = 754974721, NTT_P1 = 167772161, NTT_P2 = 469762049, NTT_P3 = 998244353, NTT_P4 = 2013265921;
template <unsigned P>
std::vector<unsigned> conv_residue(const std::vector<long long>& a, const std::vector<long long>& b) {
    using M = ModInt<P>;
    std::vector<M> x(a.size()), y(b.size());
    for (size_t i = 0; i < a.size(); ++i) x[i] = M(a[i]);
    for (size_t i = 0; i < b.size(); ++i) y[i] = M(b[i]);
    const std::vector<M> z = (std::min(a.size(), b.size()) <= 60) ? conv_naive(x, y) : conv_ntt(std::move(x), std::move(y));
    std::vector<unsigned> r(z.size());
    for (size_t i = 0; i < z.size(); ++i) r[i] = (unsigned)z[i].val();
    return r;
}
// 素数 ps[0..k) での余り rs から、−M/2 < x < M/2（M = Π ps）の x を 2^128 で割った余りとして戻す（つり合った混合基数のガーナー）
struct GarnerBalanced {
    int k;
    long long p[5], inv[5], pm[5][5];  // inv[j] = (p_0 … p_{j−1})^{−1} mod p_j、pm[i][j] = p_0 … p_{i−1} mod p_j
    UInt128 base[5];                   // p_0 … p_{j−1}（2^128 で割った余り）
    GarnerBalanced(const unsigned* ps, int k_) : k(k_) {
        for (int j = 0; j < k; ++j) p[j] = ps[j];
        for (int j = 0; j < k; ++j) {
            long long pr = 1;
            for (int i = 0; i <= j; ++i) pm[i][j] = pr, pr = pr * p[i] % p[j];
            inv[j] = mod_inv(pm[j][j], p[j]);
        }
        UInt128 b = 1;
        for (int j = 0; j < k; ++j) base[j] = b, b *= (UInt128)p[j];
    }
    Int128 operator()(const unsigned* rs) const {
        long long t[5];
        UInt128 x = 0;
        for (int j = 0; j < k; ++j) {
            long long cur = 0;  // すでに決めた部分 Σ_{i<j} t_i (p_0 … p_{i−1}) を p_j で割った余り
            for (int i = 0; i < j; ++i) cur = (cur + (t[i] % p[j] + p[j]) % p[j] * pm[i][j]) % p[j];
            long long v = ((long long)rs[j] - cur) % p[j];
            if (v < 0) v += p[j];
            v = v * inv[j] % p[j];
            if (v > p[j] / 2) v -= p[j];  // (−p/2, p/2] に
            t[j] = v;
            x += (UInt128)(Int128)v * base[j];
        }
        return (Int128)x;
    }
};
}  // namespace math_internal

// mint の列の畳み込み c[k] = Σ_{i+j=k} a[i] b[j]。mod が 998244353 などの NTT 素数なら NTT、それ以外（10^9+7・合成数）は
// 3 つの素数で計算してガーナーで戻す（自動）。短い方が 60 以下なら愚直  O((n + m) log(n + m))
template <class mint>
std::vector<mint> convolve(const std::vector<mint>& a, const std::vector<mint>& b) {
    if (a.empty() || b.empty()) return {};
    if (std::min(a.size(), b.size()) <= 60) return math_internal::conv_naive(a, b);
    const size_t z = a.size() + b.size() - 1;
    if (math_internal::ntt_usable<mint>((size_t)1 << math_internal::ceil_log2(z))) return math_internal::conv_ntt(a, b);
    if (z > (1u << 24)) math_internal::fail("convolve", "長さが 2^24 を超える");
    using namespace math_internal;
    std::vector<long long> x(a.size()), y(b.size());
    for (size_t i = 0; i < a.size(); ++i) x[i] = a[i].val();
    for (size_t i = 0; i < b.size(); ++i) y[i] = b[i].val();
    const std::vector<unsigned> r0 = conv_residue<NTT_P0>(x, y), r1 = conv_residue<NTT_P1>(x, y), r2 = conv_residue<NTT_P2>(x, y);
    // 値は 0 以上 min(n, m)·(mod − 1)^2 < P0 P1 P2 なので、ふつうのガーナーで mod の値に
    const long long m = mint::mod();
    const long long i01 = mod_inv(NTT_P0, NTT_P1), i012 = mod_inv((long long)NTT_P0 * NTT_P1 % NTT_P2, NTT_P2);
    const long long p01_m = (long long)NTT_P0 * NTT_P1 % m;
    std::vector<mint> c(z);
    for (size_t i = 0; i < z; ++i) {
        const long long t0 = r0[i];
        const long long t1 = (long long)(r1[i] + NTT_P1 - t0 % NTT_P1) % NTT_P1 * i01 % NTT_P1;
        const long long cur2 = (t0 + (long long)NTT_P0 % NTT_P2 * t1) % NTT_P2;
        const long long t2 = ((long long)r2[i] - cur2 % NTT_P2 + NTT_P2) % NTT_P2 * i012 % NTT_P2;
        c[i] = mint((t0 + (Int128)NTT_P0 * t1 + (Int128)p01_m * t2) % m);
    }
    return c;
}
namespace math_internal {
// 整数のまま畳み込む（結果の絶対値が 2^126 未満と分かっているときだけ。足りるだけの素数を使う）
inline std::vector<Int128> convolve_exact(const std::vector<long long>& a, const std::vector<long long>& b, const char* func) {
    if (a.empty() || b.empty()) return {};
    const size_t z = a.size() + b.size() - 1;
    Int128 ma = 0, mb = 0;
    for (long long v : a) ma = std::max<Int128>(ma, v < 0 ? -(Int128)v : (Int128)v);
    for (long long v : b) mb = std::max<Int128>(mb, v < 0 ? -(Int128)v : (Int128)v);
    // 係数の大きさの上限 ≤ min(n, m)·max|a|·max|b|（long double で見積もる）
    const long double bound = (long double)std::min(a.size(), b.size()) * (long double)ma * (long double)mb;
    if (bound >= 8.0e37L) fail(func, "結果が Int128 に収まらないかもしれない");
    if (std::min(a.size(), b.size()) <= 60) {  // 愚直（128 ビットで。上の見積もりで溢れない）
        std::vector<Int128> c(z, 0);
        for (size_t i = 0; i < a.size(); ++i)
            for (size_t j = 0; j < b.size(); ++j) c[i + j] += (Int128)a[i] * b[j];
        return c;
    }
    const int k = bound < 4.0e24L ? 3 : 5;  // 3 つの素数の積 ≈ 5.9×10^25、5 つ ≈ 1.2×10^44（2^126 ≈ 8.5×10^37 まで）
    if (z > (1u << (k == 3 ? 24 : 23))) fail(func, "長さが大きすぎる");
    const unsigned ps[5] = {NTT_P0, NTT_P1, NTT_P2, NTT_P3, NTT_P4};
    std::vector<unsigned> rs[5];
    rs[0] = conv_residue<NTT_P0>(a, b), rs[1] = conv_residue<NTT_P1>(a, b), rs[2] = conv_residue<NTT_P2>(a, b);
    if (k == 5) rs[3] = conv_residue<NTT_P3>(a, b), rs[4] = conv_residue<NTT_P4>(a, b);
    std::vector<Int128> c(z);
    if (k == 3) {  // 3 つの素数：ふつうのガーナーで [0, M) にしてから M/2 を超えたら M を引く
        const long long P0 = NTT_P0, P1 = NTT_P1, P2 = NTT_P2;
        const long long i01 = mod_inv(P0, P1), i012 = mod_inv(P0 * P1 % P2, P2), p0_2 = P0 % P2;
        const Int128 M = (Int128)P0 * P1 * P2, P01 = (Int128)P0 * P1;
        for (size_t i = 0; i < z; ++i) {
            const long long t0 = rs[0][i];
            const long long t1 = (rs[1][i] - t0 % P1 + P1) % P1 * i01 % P1;
            const long long t2 = ((long long)rs[2][i] - (t0 + p0_2 * t1) % P2 + P2) % P2 * i012 % P2;
            Int128 x = t0 + (Int128)P0 * t1 + P01 * t2;
            if (x > M / 2) x -= M;
            c[i] = x;
        }
        return c;
    }
    const GarnerBalanced g(ps, k);
    for (size_t i = 0; i < z; ++i) {
        unsigned r[5];
        for (int j = 0; j < k; ++j) r[j] = rs[j][i];
        c[i] = g(r);
    }
    return c;
}
}  // namespace math_internal
// 整数の列の畳み込み（厳密。負の数も。結果が long long に収まらなければ止める）
inline std::vector<long long> convolve_ll(const std::vector<long long>& a, const std::vector<long long>& b) {
    const std::vector<Int128> c = math_internal::convolve_exact(a, b, "convolve_ll");
    std::vector<long long> r(c.size());
    for (size_t i = 0; i < c.size(); ++i) {
        if (c[i] > (Int128)9223372036854775807LL || c[i] < -(Int128)9223372036854775807LL - 1) math_internal::fail("convolve_ll", "結果が long long に収まらない（convolve_i128 を使う）");
        r[i] = (long long)c[i];
    }
    return r;
}
// 整数の列の畳み込みを Int128 で（結果の絶対値は 2^126 程度まで）
inline std::vector<Int128> convolve_i128(const std::vector<long long>& a, const std::vector<long long>& b) { return math_internal::convolve_exact(a, b, "convolve_i128"); }
// long long の列を mod m（1 ≤ m < 2^31）で畳み込む
inline std::vector<long long> convolve_mod(const std::vector<long long>& a, const std::vector<long long>& b, long long m) {
    if (m < 1 || m >= (1LL << 31)) math_internal::fail("convolve_mod", "1 ≤ m < 2^31 で使う");
    using M = DynModIntT<-1>;  // 利用者の DynModInt とは別の mod を持つ
    const long long keep = M::mod();
    M::set_mod(m);
    std::vector<M> x(a.size()), y(b.size());
    for (size_t i = 0; i < a.size(); ++i) x[i] = M(a[i]);
    for (size_t i = 0; i < b.size(); ++i) y[i] = M(b[i]);
    const std::vector<M> z = convolve(x, y);
    std::vector<long long> r(z.size());
    for (size_t i = 0; i < z.size(); ++i) r[i] = z[i].val();
    M::set_mod(keep);
    return r;
}
// 多くの多項式の積（短いものから 2 つずつ掛ける）  O(N log^2 N)
template <class mint>
std::vector<mint> convolve_many(std::vector<std::vector<mint>> ps) {
    if (ps.empty()) return {mint(1)};
    for (const auto& p : ps)
        if (p.empty()) return {};
    while (ps.size() > 1) {  // 短い順に並べ、隣どうしを掛けて半分にする（短いものから掛けていく）
        std::sort(ps.begin(), ps.end(), [](const std::vector<mint>& x, const std::vector<mint>& y) { return x.size() < y.size(); });
        std::vector<std::vector<mint>> nx;
        nx.reserve(ps.size() / 2 + 1);
        for (size_t i = 0; i + 1 < ps.size(); i += 2) nx.push_back(convolve(ps[i], ps[i + 1]));
        if (ps.size() % 2) nx.push_back(std::move(ps.back()));
        ps.swap(nx);
    }
    return ps[0];
}
// 相関 c[k] = Σ_{i<|b|} a[k+i]·b[i]（k = 0〜|a|−|b|。a の中で b をずらして重ねた積の和。|a| ≥ |b| ≥ 1）
template <class mint>
std::vector<mint> correlate(const std::vector<mint>& a, const std::vector<mint>& b) {
    if (b.empty() || a.size() < b.size()) math_internal::fail("correlate", "|a| ≥ |b| ≥ 1 で使う");
    const std::vector<mint> rb(b.rbegin(), b.rend());
    const std::vector<mint> c = convolve(a, rb);
    return std::vector<mint>(c.begin() + (long)(b.size() - 1), c.begin() + (long)a.size());
}
// 巡回畳み込み c[k] = Σ_{(i+j) mod n = k} a[i] b[j]（長さ n。a, b の添字は mod n で見る）
template <class mint>
std::vector<mint> convolve_cyclic(const std::vector<mint>& a, const std::vector<mint>& b, size_t n) {
    if (n == 0) math_internal::fail("convolve_cyclic", "n ≥ 1 で使う");
    std::vector<mint> c(n, mint(0));
    const std::vector<mint> d = convolve(a, b);
    for (size_t k = 0; k < d.size(); ++k) c[k % n] += d[k];
    return c;
}
// オンラインの畳み込み：a_n, b_n を 1 つずつ渡すと、その時点で c_n = Σ_{i+j=n} a_i b_j を返す（後の値を知らずに）
//   f_{n+1} = Σ_{i=0}^{n} f_i g_{n−i} のように、答えを次の入力に使う DP に。全体で O(N log^2 N)
template <class mint>
struct RelaxedConvolution {
    mint push(const mint& x, const mint& y) {  // a_n = x、b_n = y を加えて c_n を返す（n = これまでに渡した数）
        const size_t n = a_.size();
        a_.push_back(x), b_.push_back(y);
        if (c_.size() < 2 * n + 2) c_.resize(2 * n + 2, mint(0));
        c_[n] += n == 0 ? a_[0] * b_[0] : a_[n] * b_[0] + a_[0] * b_[n];
        // 添字が 1 以上の部分を大きさ s の正方形に分けて、使う値がそろった時点で掛ける（[s, 2s) × [m s, (m+1) s) と、m ≥ 2 ならその逆）
        for (size_t s = 1; 2 * s <= n + 1 && (n + 1) % s == 0; s *= 2) {
            const size_t m = (n + 1) / s - 1;
            add_block(a_, s, b_, m * s, s, (m + 1) * s);
            if (m >= 2) add_block(b_, s, a_, m * s, s, (m + 1) * s);
        }
        return c_[n];
    }
    size_t size() const { return a_.size(); }  // これまでに渡した数
    mint get(size_t n) const {                // 確定した c_n（n < size()）
        if (n >= a_.size()) math_internal::fail("RelaxedConvolution::get", "n < size() で使う");
        return c_[n];
    }

  private:
    std::vector<mint> a_, b_, c_;
    void add_block(const std::vector<mint>& p, size_t pl, const std::vector<mint>& q, size_t ql, size_t s, size_t out) {  // p[pl, pl+s) × q[ql, ql+s) を c[out + …] に
        if (s <= 64) {
            for (size_t i = 0; i < s; ++i)
                for (size_t j = 0; j < s; ++j) c_[out + i + j] += p[pl + i] * q[ql + j];
            return;
        }
        const std::vector<mint> u(p.begin() + (long)pl, p.begin() + (long)(pl + s)), v(q.begin() + (long)ql, q.begin() + (long)(ql + s));
        const std::vector<mint> w = convolve(u, v);
        for (size_t k = 0; k < w.size(); ++k) c_[out + k] += w[k];
    }
};
#endif  // GG_CONVOLUTION
#ifndef GG_FPS
#define GG_FPS
// ---- 形式的冪級数（数学第二世代）  Fps<mint>：逆数・log・exp・累乗・平方根・割り算・テイラーシフト・P/Q の K 次の係数・多点評価・補間・べき乗和 ----
// Fps<mint> は vector<mint>（低い次数から）に + − × を足したもの。n を渡す関数は「最初の n 項（x^0〜x^{n−1}）」を返す
template <class mint>
struct Fps : std::vector<mint> {
    using std::vector<mint>::vector;
    Fps() = default;
    Fps(const std::vector<mint>& v) : std::vector<mint>(v) {}
    Fps(std::vector<mint>&& v) : std::vector<mint>(std::move(v)) {}
    Fps pre(size_t n) const {  // 最初の n 項（足りなければ 0 を足す）
        Fps r(this->begin(), this->begin() + (long)std::min(n, this->size()));
        r.resize(n);
        return r;
    }
    mint eval(const mint& x) const {  // f(x)
        mint r = 0;
        for (size_t i = this->size(); i > 0; --i) r = r * x + (*this)[i - 1];
        return r;
    }
    Fps& operator+=(const Fps& g) {
        if (g.size() > this->size()) this->resize(g.size());
        for (size_t i = 0; i < g.size(); ++i) (*this)[i] += g[i];
        return *this;
    }
    Fps& operator-=(const Fps& g) {
        if (g.size() > this->size()) this->resize(g.size());
        for (size_t i = 0; i < g.size(); ++i) (*this)[i] -= g[i];
        return *this;
    }
    Fps& operator*=(const Fps& g) { return *this = Fps(convolve<mint>(*this, g)); }
    Fps& operator*=(const mint& s) {
        for (mint& x : *this) x *= s;
        return *this;
    }
    friend Fps operator+(Fps f, const Fps& g) { return f += g; }
    friend Fps operator-(Fps f, const Fps& g) { return f -= g; }
    friend Fps operator*(const Fps& f, const Fps& g) { return Fps(convolve<mint>(f, g)); }
    friend Fps operator*(Fps f, const mint& s) { return f *= s; }
    friend Fps operator*(const mint& s, Fps f) { return f *= s; }
    Fps operator-() const {
        Fps r(*this);
        for (mint& x : r) x = -x;
        return r;
    }
    friend std::ostream& operator<<(std::ostream& os, const Fps& f) {  // 係数を空白区切りで（そのまま答えの出力に使える）
        for (size_t i = 0; i < f.size(); ++i) os << (i ? " " : "") << f[i];
        return os;
    }
};

namespace math_internal {
template <class mint>
void fps_need_prime(long long n, const char* func) {  // 1〜n の逆数を使う関数（mod は n より大きい素数）
    if (!is_prime_u32((unsigned)mint::mod()) || n >= mint::mod()) fail(func, "mod が項の数より大きい素数のときに使う");
}
template <class mint>
std::vector<mint> fps_inverses(int n) {  // 1/0（使わない）, 1/1, …, 1/n
    std::vector<mint> iv((size_t)n + 1, mint(1));
    const long long p = mint::mod();
    for (int i = 2; i <= n; ++i) iv[(size_t)i] = -iv[(size_t)(p % i)] * mint(p / i);
    return iv;
}
}  // namespace math_internal

template <class mint>
Fps<mint> fps_derivative(const std::vector<mint>& f) {  // f'
    Fps<mint> r(f.size() > 1 ? f.size() - 1 : 0);
    for (size_t i = 1; i < f.size(); ++i) r[i - 1] = f[i] * mint((long long)i);
    return r;
}
template <class mint>
Fps<mint> fps_integral(const std::vector<mint>& f) {  // ∫f（定数項 0）
    math_internal::fps_need_prime<mint>((long long)f.size(), "fps_integral");
    const std::vector<mint> iv = math_internal::fps_inverses<mint>((int)f.size());
    Fps<mint> r(f.size() + 1);
    for (size_t i = 0; i < f.size(); ++i) r[i + 1] = f[i] * iv[i + 1];
    return r;
}
// 1/f の最初の n 項（f[0] ≠ 0。f[0] が mod と互いに素なら合成数の mod でも）  O(n log n)
template <class mint>
Fps<mint> fps_inv(const std::vector<mint>& f, int n) {
    if (n < 0) math_internal::fail("fps_inv", "n ≥ 0 で使う");
    if (n == 0) return {};
    if (f.empty() || f[0] == mint(0)) math_internal::fail("fps_inv", "f[0] が 0");
    Fps<mint> g{f[0].inv()};
    const bool fast = math_internal::ntt_usable<mint>((size_t)1 << math_internal::ceil_log2((size_t)2 * n));
    for (int m = 1; m < n; m *= 2) {
        if (fast) {  // NTT 5 回：g ← g − g (f g − 1)
            std::vector<mint> F((size_t)2 * m), G((size_t)2 * m);
            for (int i = 0; i < 2 * m && i < (int)f.size(); ++i) F[i] = f[i];
            for (int i = 0; i < m; ++i) G[i] = g[i];
            ntt(F), ntt(G);
            for (int i = 0; i < 2 * m; ++i) F[i] *= G[i];
            intt(F);
            for (int i = 0; i < m; ++i) F[i] = 0;  // f g の最初の m 項は 1, 0, …（使わない）
            ntt(F);
            for (int i = 0; i < 2 * m; ++i) F[i] *= G[i];
            intt(F);
            g.resize((size_t)2 * m);
            for (int i = m; i < 2 * m; ++i) g[i] = -F[i];
        } else {  // どの mod でも：g ← g (2 − f g)
            Fps<mint> h = Fps<mint>(f).pre((size_t)2 * m) * g;
            h.resize((size_t)2 * m);
            for (mint& x : h) x = -x;
            h[0] += 2;
            g = (g * h).pre((size_t)2 * m);
        }
    }
    return g.pre((size_t)n);
}
// log f の最初の n 項（f[0] = 1。mod は n より大きい素数）  O(n log n)
template <class mint>
Fps<mint> fps_log(const std::vector<mint>& f, int n) {
    if (n < 0) math_internal::fail("fps_log", "n ≥ 0 で使う");
    if (n == 0) return {};
    if (f.empty() || f[0] != mint(1)) math_internal::fail("fps_log", "f[0] = 1 で使う");
    Fps<mint> d = fps_derivative(Fps<mint>(f).pre((size_t)n)) * fps_inv(f, n);
    d.resize((size_t)n - 1);
    return fps_integral(d).pre((size_t)n);
}
// exp f の最初の n 項（f[0] = 0。mod は n より大きい素数）  O(n log n)
template <class mint>
Fps<mint> fps_exp(const std::vector<mint>& f, int n) {
    if (n < 0) math_internal::fail("fps_exp", "n ≥ 0 で使う");
    if (n == 0) return {};
    if (!f.empty() && f[0] != mint(0)) math_internal::fail("fps_exp", "f[0] = 0 で使う");
    math_internal::fps_need_prime<mint>(n, "fps_exp");
    Fps<mint> g{mint(1)};
    for (int m = 1; m < n; m *= 2) {  // g ← g (1 − log g + f)（最後の段は n 項まで：mod が小さい素数でも 1/n までしか使わない）
        const int len = std::min(2 * m, n);
        Fps<mint> h = fps_log(g, len);
        for (int i = 0; i < len; ++i) h[i] = (i < (int)f.size() ? f[i] : mint(0)) - h[i];
        h[0] += 1;
        g = (g * h).pre((size_t)len);
    }
    return g.pre((size_t)n);
}
// f^k の最初の n 項（k ≥ 0。k は 10^18 でもよい。f[0] = 0 でもよい。mod は n より大きい素数）  O(n log n)
template <class mint>
Fps<mint> fps_pow(const std::vector<mint>& f, long long k, int n) {
    if (n < 0 || k < 0) math_internal::fail("fps_pow", "n, k ≥ 0 で使う");
    Fps<mint> r((size_t)n);
    if (n == 0) return r;
    if (k == 0) {
        r[0] = 1;
        return r;
    }
    size_t d = 0;
    while (d < f.size() && f[d] == mint(0)) ++d;
    if (d == f.size() || (Int128)d * k >= n) return r;  // 全部 0
    const size_t shift = (size_t)d * (size_t)k, m = (size_t)n - shift;
    const mint c = f[d], ic = c.inv();
    Fps<mint> g(f.begin() + (long)d, f.begin() + (long)std::min(f.size(), d + m));
    for (mint& x : g) x *= ic;  // g = f / (c x^d)、g[0] = 1
    Fps<mint> lg = fps_log(g, (int)m);
    for (mint& x : lg) x *= mint(k);
    const Fps<mint> e = fps_exp(lg, (int)m);
    const mint ck = c.pow(k);
    for (size_t i = 0; i < m; ++i) r[shift + i] = e[i] * ck;
    return r;
}
namespace math_internal {
template <class mint>
Maybe<mint> fps_mod_sqrt(const mint& a) {  // mod は奇素数。トネリ・シャンクス
    const long long p = mint::mod();
    if (a == mint(0)) return {true, mint(0)};
    if (p == 2) return {true, a};
    if (a.pow((p - 1) / 2) != mint(1)) return {false, mint(0)};
    long long q = p - 1;
    int s = 0;
    while (q % 2 == 0) q /= 2, ++s;
    mint z = 2;
    while (z.pow((p - 1) / 2) == mint(1)) z += 1;
    mint c = z.pow(q), x = a.pow((q + 1) / 2), t = a.pow(q);
    int m = s;
    while (t != mint(1)) {
        int i = 0;
        mint tt = t;
        while (tt != mint(1)) tt *= tt, ++i;
        mint b = c;
        for (int j = 0; j < m - i - 1; ++j) b *= b;
        x *= b, c = b * b, t *= c, m = i;
    }
    return {true, x};
}
}  // namespace math_internal
// √f の最初の n 項（なければ ok = false。f[0] = 0 でもよい。mod は奇素数）  O(n log n)
template <class mint>
Maybe<Fps<mint>> fps_sqrt(const std::vector<mint>& f, int n) {
    if (n < 0) math_internal::fail("fps_sqrt", "n ≥ 0 で使う");
    if (mint::mod() == 2 || !math_internal::is_prime_u32((unsigned)mint::mod())) math_internal::fail("fps_sqrt", "mod が奇素数のときに使う");
    Fps<mint> r((size_t)n);
    if (n == 0) return {true, r};
    size_t d = 0;
    while (d < f.size() && f[d] == mint(0)) ++d;
    if (d == f.size()) return {true, r};  // f = 0
    if (d % 2 == 1) return {false, Fps<mint>()};
    const Maybe<mint> s0 = math_internal::fps_mod_sqrt(f[d]);
    if (!s0.ok) return {false, Fps<mint>()};
    if (d / 2 >= (size_t)n) return {true, r};  // 最初の n 項は 0
    const size_t m = (size_t)n - d / 2;
    const Fps<mint> g(f.begin() + (long)d, f.begin() + (long)std::min(f.size(), d + m));
    Fps<mint> h{s0.value};
    const mint inv2 = mint(2).inv();
    for (size_t k = 1; k < m; k *= 2) {  // h ← (h + g/h) / 2
        Fps<mint> t = (g.pre(2 * k) * fps_inv(h, (int)(2 * k))).pre(2 * k);
        h.resize(2 * k);
        for (size_t i = 0; i < 2 * k; ++i) h[i] = (h[i] + t[i]) * inv2;
    }
    for (size_t i = 0; i < m; ++i) r[d / 2 + i] = h[i];
    return {true, r};
}
template <class mint>
struct FpsDivMod {  // f = q g + r（deg r < deg g）
    Fps<mint> q, r;
};
// 多項式の割り算（g の最高次の係数で割る。末尾の 0 は除く）  O(n log n)
template <class mint>
FpsDivMod<mint> fps_divmod(std::vector<mint> f, std::vector<mint> g) {
    while (!f.empty() && f.back() == mint(0)) f.pop_back();
    while (!g.empty() && g.back() == mint(0)) g.pop_back();
    if (g.empty()) math_internal::fail("fps_divmod", "0 で割った");
    FpsDivMod<mint> res;
    if (f.size() < g.size()) {
        res.r = f;
        return res;
    }
    const size_t k = f.size() - g.size() + 1;
    Fps<mint> rf(f.rbegin(), f.rend()), rg(g.rbegin(), g.rend());
    Fps<mint> q = (rf.pre(k) * fps_inv(rg, (int)k)).pre(k);
    std::reverse(q.begin(), q.end());
    Fps<mint> r = Fps<mint>(f) - Fps<mint>(g) * q;
    r.resize(g.size() - 1);
    while (!r.empty() && r.back() == mint(0)) r.pop_back();
    while (!q.empty() && q.back() == mint(0)) q.pop_back();
    res.q = q, res.r = r;
    return res;
}
// f(x + c) の係数（mod は f の長さより大きい素数）  O(n log n)
template <class mint>
Fps<mint> taylor_shift(const std::vector<mint>& f, const mint& c) {
    const int n = (int)f.size();
    if (n == 0) return {};
    math_internal::fps_need_prime<mint>(n, "taylor_shift");
    std::vector<mint> fact((size_t)n, mint(1)), ifact((size_t)n, mint(1));
    for (int i = 1; i < n; ++i) fact[i] = fact[i - 1] * mint(i);
    ifact[n - 1] = fact[n - 1].inv();
    for (int i = n - 1; i > 0; --i) ifact[i - 1] = ifact[i] * mint(i);
    std::vector<mint> a((size_t)n), b((size_t)n);
    mint pw = 1;
    for (int i = 0; i < n; ++i) a[i] = f[i] * fact[i], b[(size_t)(n - 1 - i)] = pw * ifact[i], pw *= c;
    const std::vector<mint> conv = convolve(a, b);
    Fps<mint> r((size_t)n);
    for (int j = 0; j < n; ++j) r[j] = conv[(size_t)(n - 1 + j)] * ifact[j];
    return r;
}
// [x^K] P(x)/Q(x)（Q[0] ≠ 0。K ≤ 10^18。Bostan–Mori）  O(d log d log K)（d = 次数）
template <class mint>
mint fps_coeff_of_fraction(std::vector<mint> P, std::vector<mint> Q, long long K) {
    if (K < 0) return mint(0);
    if (Q.empty() || Q[0] == mint(0)) math_internal::fail("fps_coeff_of_fraction", "Q[0] が 0");
    while (K > 0) {
        if (P.empty()) return mint(0);
        const size_t up = P.size() + Q.size() - 1, vq = 2 * Q.size() - 1;
        const size_t n = (size_t)1 << math_internal::ceil_log2(std::max(up, vq));
        if (Q.size() >= 16 && math_internal::ntt_usable<mint>(n)) {  // NTT の値の並びでは j と j^1 が ±x の組（Q(−x) の値は隣と入れ替えるだけ）
            P.resize(n), Q.resize(n);
            ntt(P), ntt(Q);
            for (size_t j = 0; j < n; j += 2) {
                const mint q0 = Q[j], q1 = Q[j + 1];
                P[j] *= q1, P[j + 1] *= q0;
                Q[j] = Q[j + 1] = q0 * q1;
            }
            intt(P), intt(Q);
            P.resize(up), Q.resize(vq);
        } else {
            std::vector<mint> Qm(Q);
            for (size_t i = 1; i < Qm.size(); i += 2) Qm[i] = -Qm[i];  // Q(−x)
            P = convolve(P, Qm), Q = convolve(Q, Qm);
        }
        std::vector<mint> np, nq;
        for (size_t i = (size_t)(K & 1); i < P.size(); i += 2) np.push_back(P[i]);
        for (size_t i = 0; i < Q.size(); i += 2) nq.push_back(Q[i]);
        P.swap(np), Q.swap(nq), K >>= 1;
    }
    return P.empty() ? mint(0) : P[0] / Q[0];
}
// Π (1 + a_i x)（分割統治）  O(n log^2 n)
template <class mint>
Fps<mint> prod_linear(const std::vector<mint>& as) {
    std::vector<std::vector<mint>> ps;
    ps.reserve(as.size());
    for (const mint& a : as) ps.push_back({mint(1), a});
    return Fps<mint>(convolve_many(std::move(ps)));
}

// ---- 多点評価・補間・分数の和 ----
namespace math_internal {
template <class mint>
struct FpsProductTree {  // 区間 [l, r) の点の Π (x − x_i) を持つ木（葉の近くは 32 点ずつまとめる）
    static constexpr int LEAF = 32;
    const std::vector<mint>& xs;
    std::vector<Fps<mint>> t;
    explicit FpsProductTree(const std::vector<mint>& x) : xs(x), t(4 * std::max<size_t>(x.size(), 1)) { build(1, 0, (int)x.size()); }
    void build(int k, int l, int r) {
        if (r - l <= LEAF) {
            Fps<mint> p = {mint(1)};
            for (int i = l; i < r; ++i) {  // p ← p (x − x_i)
                p.push_back(mint(0));
                for (size_t j = p.size() - 1; j > 0; --j) p[j] = p[j - 1] - p[j] * xs[(size_t)i];
                p[0] = -p[0] * xs[(size_t)i];
            }
            t[(size_t)k] = p;
            return;
        }
        const int m = (l + r) / 2;
        build(2 * k, l, m), build(2 * k + 1, m, r);
        t[(size_t)k] = t[(size_t)(2 * k)] * t[(size_t)(2 * k + 1)];
    }
    void eval(int k, int l, int r, Fps<mint> f, std::vector<mint>& out) const {  // f の次数は r − l 未満にしてから
        if (r - l <= LEAF) {
            for (int i = l; i < r; ++i) out[(size_t)i] = f.eval(xs[(size_t)i]);
            return;
        }
        const int m = (l + r) / 2;
        eval(2 * k, l, m, fps_divmod(f, t[(size_t)(2 * k)]).r, out);
        eval(2 * k + 1, m, r, fps_divmod(f, t[(size_t)(2 * k + 1)]).r, out);
    }
    Fps<mint> combine(int k, int l, int r, const std::vector<mint>& c) const {  // Σ c_i Π_{j ≠ i} (x − x_j)
        if (r - l <= LEAF) {
            Fps<mint> res((size_t)(r - l), mint(0));
            for (int i = l; i < r; ++i) {  // t[k] / (x − x_i) を組み立て除法で
                const Fps<mint>& p = t[(size_t)k];
                mint carry = 0;
                for (size_t j = p.size() - 1; j > 0; --j) {
                    carry = p[j] + carry * xs[(size_t)i];
                    res[j - 1] += carry * c[(size_t)i];
                }
            }
            return res;
        }
        const int m = (l + r) / 2;
        return combine(2 * k, l, m, c) * t[(size_t)(2 * k + 1)] + combine(2 * k + 1, m, r, c) * t[(size_t)(2 * k)];
    }
};
}  // namespace math_internal
// 多点評価：f(x_0), f(x_1), …, f(x_{m−1})  O(n log n + m log^2 m)（n = |f|, m = |xs|）
template <class mint>
std::vector<mint> multipoint_eval(const std::vector<mint>& f, const std::vector<mint>& xs) {
    std::vector<mint> out(xs.size());
    if (xs.empty()) return out;
    if ((long long)f.size() * (long long)xs.size() <= 1000000 || std::min(f.size(), xs.size()) <= 64) {  // 小さいとき・片方が短いときは 1 点ずつ
        const Fps<mint> F(f);
        for (size_t i = 0; i < xs.size(); ++i) out[i] = F.eval(xs[i]);
        return out;
    }
    const math_internal::FpsProductTree<mint> T(xs);
    T.eval(1, 0, (int)xs.size(), fps_divmod(f, T.t[1]).r, out);
    return out;
}
// 補間：x_i（すべて異なる）で y_i をとる次数 n − 1 以下の多項式（係数は低い次数から。長さ n）  O(n log^2 n)
template <class mint>
Fps<mint> interpolate(const std::vector<mint>& xs, const std::vector<mint>& ys) {
    if (xs.size() != ys.size()) math_internal::fail("interpolate", "xs と ys の長さが違う");
    const int n = (int)xs.size();
    if (n == 0) return {};
    const math_internal::FpsProductTree<mint> T(xs);
    const std::vector<mint> d = multipoint_eval<mint>(fps_derivative(T.t[1]), xs);  // M'(x_i) = Π_{j ≠ i} (x_i − x_j)
    std::vector<mint> c((size_t)n);
    for (int i = 0; i < n; ++i) {
        if (d[(size_t)i] == mint(0)) math_internal::fail("interpolate", "x が重複している");
        c[(size_t)i] = ys[(size_t)i] / d[(size_t)i];
    }
    Fps<mint> r = T.combine(1, 0, n, c);
    r.resize((size_t)n);
    return r;
}
template <class mint>
struct FpsFraction {  // num / den
    Fps<mint> num, den;
};
// 分数の和 Σ nums[i] / dens[i] を 1 つの分数に（分割統治。通分するだけで約分はしない）  O(N log^2 N)（N = 次数の合計）
template <class mint>
FpsFraction<mint> fps_sum_of_fractions(const std::vector<Fps<mint>>& nums, const std::vector<Fps<mint>>& dens) {
    if (nums.size() != dens.size()) math_internal::fail("fps_sum_of_fractions", "nums と dens の長さが違う");
    if (nums.empty()) return {Fps<mint>{}, Fps<mint>{mint(1)}};
    auto rec = [&](auto&& self, size_t l, size_t r) -> FpsFraction<mint> {
        if (r - l == 1) return {nums[l], dens[l]};
        const size_t m = (l + r) / 2;
        const FpsFraction<mint> a = self(self, l, m), b = self(self, m, r);
        return {a.num * b.den + b.num * a.den, a.den * b.den};
    };
    return rec(rec, 0, nums.size());
}
// べき乗和 p_k = Σ_i w_i a_i^k（k = 0〜n−1。ws を省くと w_i = 1）= Σ w_i / (1 − a_i x) の係数  O(m log^2 m + n log n)
template <class mint>
Fps<mint> power_sums(const std::vector<mint>& as, int n, const std::vector<mint>& ws = {}) {
    if (!ws.empty() && ws.size() != as.size()) math_internal::fail("power_sums", "ws と as の長さが違う");
    if (n <= 0) return {};
    if (as.empty()) return Fps<mint>((size_t)n, mint(0));
    std::vector<Fps<mint>> nums(as.size()), dens(as.size());
    for (size_t i = 0; i < as.size(); ++i) nums[i] = {ws.empty() ? mint(1) : ws[i]}, dens[i] = {mint(1), -as[i]};
    const FpsFraction<mint> fr = fps_sum_of_fractions(nums, dens);
    return (fr.num.pre((size_t)n) * fps_inv(fr.den, n)).pre((size_t)n);
}
#endif  // GG_FPS
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace fps_test {

#define FPS_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[fps] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

template <class mint>
static std::vector<mint> rand_vec(std::mt19937_64& rng, size_t n) {
    std::vector<mint> v(n);
    for (auto& x : v) x = rng() % 5 == 0 ? mint(rnd(rng, 0, 3)) : mint((long long)rng());
    return v;
}
template <class mint>
static std::vector<mint> mul_trunc(const std::vector<mint>& a, const std::vector<mint>& b, size_t n) {
    std::vector<mint> c(n);
    for (size_t i = 0; i < a.size() && i < n; ++i)
        for (size_t j = 0; j < b.size() && i + j < n; ++j) c[i + j] += a[i] * b[j];
    return c;
}
template <class mint>
static std::vector<mint> brute_inv(const std::vector<mint>& f, size_t n) {  // g_k = −(Σ_{i≥1} f_i g_{k−i}) / f_0
    std::vector<mint> g(n);
    const mint i0 = f[0].inv();
    for (size_t k = 0; k < n; ++k) {
        mint s = k == 0 ? mint(1) : mint(0);
        for (size_t i = 1; i <= k && i < f.size(); ++i) s -= f[i] * g[k - i];
        g[k] = s * i0;
    }
    return g;
}
template <class mint>
static std::vector<mint> brute_exp(const std::vector<mint>& f, size_t n) {  // g_k = (1/k) Σ_{i=1}^{k} i f_i g_{k−i}
    std::vector<mint> g(n);
    if (n) g[0] = 1;
    for (size_t k = 1; k < n; ++k) {
        mint s = 0;
        for (size_t i = 1; i <= k && i < f.size(); ++i) s += mint((long long)i) * f[i] * g[k - i];
        g[k] = s / mint((long long)k);
    }
    return g;
}
template <class mint>
static std::vector<mint> brute_log(const std::vector<mint>& f, size_t n) {  // f L' = f'
    std::vector<mint> d(n), L(n);
    for (size_t k = 0; k + 1 < n; ++k) {
        mint s = (k + 1 < f.size()) ? f[k + 1] * mint((long long)(k + 1)) : mint(0);
        for (size_t i = 1; i <= k && i < f.size(); ++i) s -= f[i] * d[k - i];
        d[k] = s;
        L[k + 1] = d[k] / mint((long long)(k + 1));
    }
    return L;
}

template <class mint>
static void test_basic(std::mt19937_64& rng, int& cases, int iters, bool prime_mod) {
    for (int it = 0; it < iters; ++it, ++cases) {
        const size_t n = (size_t)(it % 3 == 0 ? rnd(rng, 0, 8) : rnd(rng, 1, it % 2 ? 80 : 400));
        std::vector<mint> f = rand_vec<mint>(rng, (size_t)rnd(rng, 1, (long long)n + 3));
        if (f[0] == mint(0)) f[0] = 1;
        if (!prime_mod)
            while (!f[0].is_unit()) f[0] += 1;
        FPS_CHECK(fps_inv(f, (int)n) == Fps<mint>(brute_inv(f, n)));
        if (!prime_mod) continue;
        // log・exp
        std::vector<mint> g = f;
        g[0] = 1;
        FPS_CHECK(fps_log(g, (int)n) == Fps<mint>(brute_log(g, n)));
        g[0] = 0;
        FPS_CHECK(fps_exp(g, (int)n) == Fps<mint>(brute_exp(g, n)));
        // 累乗（先頭に 0 があってもよい）
        std::vector<mint> h(f);
        if (it % 3 == 0) h.insert(h.begin(), (size_t)rnd(rng, 1, 4), mint(0));
        const long long k = rnd(rng, 0, 5);
        std::vector<mint> e(n);
        if (n) e[0] = 1;
        for (long long t = 0; t < k; ++t) e = mul_trunc(e, h, n);
        FPS_CHECK(fps_pow(h, k, (int)n) == Fps<mint>(e));
        const long long k1 = rnd(rng, 0, 1000000000000000000LL), k2 = rnd(rng, 0, 1000);
        if (it % 4 == 0) FPS_CHECK(Fps<mint>(mul_trunc<mint>(fps_pow(f, k1, (int)n), fps_pow(f, k2, (int)n), n)) == fps_pow(f, k1 + k2, (int)n));
        // 平方根：g^2 の平方根の 2 乗がもとに戻る（先頭の 0 は偶数個）
        std::vector<mint> sq = mul_trunc(h, h, n + 6);
        const Maybe<Fps<mint>> s = fps_sqrt(sq, (int)n);
        FPS_CHECK(s.ok && Fps<mint>(mul_trunc<mint>(s.value, s.value, n)) == Fps<mint>(sq).pre(n));
        // 微分・積分
        FPS_CHECK(fps_derivative(fps_integral(f)) == Fps<mint>(f));
    }
}

template <class mint>
static void test_poly(std::mt19937_64& rng, int& cases, int iters) {
    for (int it = 0; it < iters; ++it, ++cases) {
        // 割り算
        std::vector<mint> f = rand_vec<mint>(rng, (size_t)rnd(rng, 0, 300)), g = rand_vec<mint>(rng, (size_t)rnd(rng, 1, 150));
        if (g.back() == mint(0)) g.back() = 1;
        const FpsDivMod<mint> d = fps_divmod(f, g);
        Fps<mint> back = Fps<mint>(g) * d.q + d.r;
        while (!back.empty() && back.back() == mint(0)) back.pop_back();
        while (!f.empty() && f.back() == mint(0)) f.pop_back();
        FPS_CHECK(back == Fps<mint>(f) && d.r.size() < g.size() && (d.r.empty() || d.r.back() != mint(0)));
        // テイラーシフト：いくつかの点で f(x + c) と比べる
        const Fps<mint> p(rand_vec<mint>(rng, (size_t)rnd(rng, 0, 200)));
        const mint c = (long long)rng();
        const Fps<mint> sh = taylor_shift(p, c);
        for (int t = 0; t < 3; ++t) {
            const mint x = (long long)rng();
            FPS_CHECK(sh.eval(x) == p.eval(x + c));
        }
        // P/Q の K 次の係数：級数を展開したものと
        std::vector<mint> P = rand_vec<mint>(rng, (size_t)rnd(rng, 0, 20)), Q = rand_vec<mint>(rng, (size_t)rnd(rng, 1, 20));
        if (Q[0] == mint(0)) Q[0] = 1;
        const size_t K = (size_t)rnd(rng, 0, 300);
        const std::vector<mint> ser = mul_trunc(P, brute_inv(Q, K + 1), K + 1);
        FPS_CHECK(fps_coeff_of_fraction(P, Q, (long long)K) == ser[K]);
        if (it % 4 == 0) {  // 次数が大きい（NTT で Q(−x) を作る方）・P の方が長い
            std::vector<mint> P2 = rand_vec<mint>(rng, (size_t)rnd(rng, 0, 300)), Q2 = rand_vec<mint>(rng, (size_t)rnd(rng, 16, 200));
            Q2[0] = 1;
            const size_t K2 = (size_t)rnd(rng, 0, 1500);
            const std::vector<mint> ser2 = mul_trunc(P2, brute_inv(Q2, K2 + 1), K2 + 1);
            FPS_CHECK(fps_coeff_of_fraction(P2, Q2, (long long)K2) == ser2[K2]);
        }
        // Π (1 + a_i x)
        const std::vector<mint> as = rand_vec<mint>(rng, (size_t)rnd(rng, 0, 100));
        std::vector<mint> e = {mint(1)};
        for (const mint& a : as) e = mul_trunc<mint>(e, {mint(1), a}, e.size() + 1);
        FPS_CHECK(prod_linear(as) == Fps<mint>(e));
    }
    // フィボナッチ：x / (1 − x − x^2) の K 次の係数を倍化の公式と
    for (int it = 0; it < 200; ++it, ++cases) {
        const long long K = rnd(rng, 0, 1000000000000000000LL);
        mint a = 0, b = 1;  // (F_k, F_{k+1}) を上の桁から
        for (int bit = 62; bit >= 0; --bit) {
            const mint c2 = a * (b * 2 - a), d2 = a * a + b * b;
            a = c2, b = d2;
            if ((K >> bit) & 1) {
                const mint t = a + b;
                a = b, b = t;
            }
        }
        FPS_CHECK(fps_coeff_of_fraction<mint>({0, 1}, {1, -1, -1}, K) == a);
    }
}

static void test_edges(int& cases) {
    using mint = ModInt998;
    FPS_CHECK(fps_inv<mint>({mint(2)}, 0).empty() && fps_pow<mint>({0, 0, 3}, 2, 4) == Fps<mint>{0, 0, 0, 0} && fps_pow<mint>({0, 0, 3}, 1, 4) == Fps<mint>{0, 0, 3, 0});
    FPS_CHECK(fps_pow<mint>({}, 0, 3) == Fps<mint>{1, 0, 0} && fps_pow<mint>({}, 5, 3) == Fps<mint>{0, 0, 0});
    FPS_CHECK(!fps_sqrt<mint>({0, 1}, 5).ok && !fps_sqrt<mint>({3}, 2).ok && fps_sqrt<mint>({0, 0, 4}, 3).value.pre(2) * fps_sqrt<mint>({0, 0, 4}, 3).value.pre(2) == Fps<mint>{0, 0, 4});  // 3 は平方剰余でない
    FPS_CHECK(fps_coeff_of_fraction<mint>({}, {1}, 5) == mint(0) && fps_coeff_of_fraction<mint>({7}, {1}, 0) == mint(7));
    {  // 小さい素数の mod：項の数 < mod なら exp・pow も使える
        DynModInt::set_mod(7);
        using D = DynModInt;
        const Fps<D> e = fps_exp(Fps<D>{0, 1}, 5);
        FPS_CHECK(e == Fps<D>{1, 1, D(2).inv(), D(6).inv(), D(24).inv()});
        FPS_CHECK(fps_pow(Fps<D>{1, 1}, 3, 5) == Fps<D>{1, 3, 3, 1, 0} && fps_exp(Fps<D>{0, 3}, 6).size() == 6);
        DynModInt::set_mod(998244353);
    }
    cases += 6;
}

template <class mint>
static void test_points(std::mt19937_64& rng, int& cases, int iters) {
    for (int it = 0; it < iters; ++it, ++cases) {
        // 多点評価（小さいときは 1 点ずつ、大きいときは木）
        const int n = (int)rnd(rng, 0, it % 4 ? 60 : 3000), m = (int)rnd(rng, 0, it % 4 ? 60 : 3000);
        Fps<mint> f((size_t)n);
        for (auto& x : f) x = (long long)rng();
        std::vector<mint> xs((size_t)m);
        for (auto& x : xs) x = it % 3 ? mint((long long)rng()) : mint(rnd(rng, 0, 5));  // 同じ点が何度も出てもよい
        const std::vector<mint> ev = multipoint_eval<mint>(f, xs);
        for (int i = 0; i < m; i += std::max(1, m / 50)) FPS_CHECK(ev[(size_t)i] == f.eval(xs[(size_t)i]));
        // 補間（異なる点）
        const int k = (int)rnd(rng, 0, it % 4 ? 70 : 2000);
        std::vector<mint> px, py;
        std::set<long long> used;
        while ((int)px.size() < k) {
            const long long v = it % 2 ? rnd(rng, 0, 1000000) : (long long)(rng() % mint::mod());
            if (used.insert(v).second) px.push_back(mint(v)), py.push_back(mint((long long)rng()));
        }
        const Fps<mint> P = interpolate(px, py);
        FPS_CHECK((int)P.size() == k);
        for (int i = 0; i < k; i += std::max(1, k / 50)) FPS_CHECK(P.eval(px[(size_t)i]) == py[(size_t)i]);
        // べき乗和
        const int na = (int)rnd(rng, 0, it % 4 ? 30 : 1000), L = (int)rnd(rng, 0, 40);
        std::vector<mint> as((size_t)na), ws((size_t)na);
        for (auto& x : as) x = (long long)rng();
        for (auto& x : ws) x = (long long)rng();
        const Fps<mint> p1 = power_sums(as, L), pw = power_sums(as, L, ws);
        FPS_CHECK((int)p1.size() == L && (int)pw.size() == L);
        for (int j = 0; j < L; j += std::max(1, L / 10)) {
            mint e1 = 0, ew = 0;
            for (int i = 0; i < na; ++i) e1 += as[(size_t)i].pow(j), ew += ws[(size_t)i] * as[(size_t)i].pow(j);
            FPS_CHECK(p1[(size_t)j] == e1 && pw[(size_t)j] == ew);
        }
    }
    // 分数の和：Σ P_i/Q_i を 1 つにまとめた N/D が、それぞれの和と同じ級数になる
    for (int it = 0; it < 50 * TEST_SCALE; ++it, ++cases) {
        const int c = (int)rnd(rng, 1, 8), L = 20;
        std::vector<Fps<mint>> nums((size_t)c), dens((size_t)c);
        Fps<mint> e((size_t)L, mint(0));
        for (int i = 0; i < c; ++i) {
            nums[(size_t)i].resize((size_t)rnd(rng, 0, 4)), dens[(size_t)i].resize((size_t)rnd(rng, 1, 4));
            for (auto& x : nums[(size_t)i]) x = (long long)rng();
            for (auto& x : dens[(size_t)i]) x = (long long)rng();
            dens[(size_t)i][0] = rnd(rng, 1, 100);
            e += (nums[(size_t)i] * fps_inv(dens[(size_t)i], L)).pre(L);
        }
        const FpsFraction<mint> fr = fps_sum_of_fractions(nums, dens);
        FPS_CHECK((fr.num * fps_inv(fr.den, L)).pre(L) == e);
    }
    {
        const std::vector<mint> xs = {1, 2, 3};
        FPS_CHECK(interpolate<mint>(xs, {2, 4, 6}) == Fps<mint>({0, 2, 0}) && multipoint_eval<mint>({}, xs) == std::vector<mint>(3, mint(0)));
        FPS_CHECK(interpolate<mint>({}, {}).empty() && power_sums<mint>({}, 3) == Fps<mint>(3, mint(0)));
        cases += 2;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_points<ModInt998>(rng, cases, 60 * TEST_SCALE);
    test_points<ModInt107>(rng, cases, 20 * TEST_SCALE);
    test_basic<ModInt998>(rng, cases, 300 * TEST_SCALE, true);
    test_basic<ModInt107>(rng, cases, 150 * TEST_SCALE, true);  // NTT が使えない素数（3 つの素数の畳み込み）
    DynModInt::set_mod(1000000000);
    test_basic<DynModInt>(rng, cases, 100 * TEST_SCALE, false);  // 合成数の mod でも fps_inv は使える
    test_poly<ModInt998>(rng, cases, 200 * TEST_SCALE);
    test_poly<ModInt107>(rng, cases, 60 * TEST_SCALE);
    test_edges(cases);
    std::cout << "[fps] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[fps] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    using mint = ModInt998;
    const int n = 500000;
    std::vector<mint> f(n);
    for (auto& x : f) x = (long long)rng();
    f[0] = 1;
    auto t0 = std::chrono::steady_clock::now();
    const auto gi = fps_inv(f, n);
    const long long t_inv = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto gl = fps_log(f, n);
    const long long t_log = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto ge = fps_exp(gl, n);
    const long long t_exp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto gp = fps_pow(f, 1000000000000000000LL, n);
    const long long t_pow = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto gs = fps_sqrt(f, n);
    const long long t_sqrt = ms(t0);
    std::vector<mint> Q(100001);
    for (auto& x : Q) x = (long long)rng();
    Q[0] = 1;
    t0 = std::chrono::steady_clock::now();
    const mint cf = fps_coeff_of_fraction<mint>({1}, Q, 1000000000000000000LL);
    const long long t_bm = ms(t0);
    std::vector<mint> xs(100000), ys(100000);
    for (size_t i = 0; i < xs.size(); ++i) xs[i] = mint((long long)i * 7 + 3), ys[i] = (long long)rng();
    const std::vector<mint> f1(f.begin(), f.begin() + 100000);
    t0 = std::chrono::steady_clock::now();
    const auto ev = multipoint_eval(f1, xs);
    const long long t_mp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto ip = interpolate(xs, ys);
    const long long t_ip = ms(t0);
    std::cout << "[fps] benchmark: n = 5×10^5 で inv " << t_inv << " ms, log " << t_log << " ms, exp " << t_exp << " ms, pow(k = 10^18) " << t_pow << " ms, sqrt " << t_sqrt
              << " ms, coeff_of_fraction（次数 10^5、K = 10^18）" << t_bm << " ms\n"
              << "[fps] benchmark: multipoint_eval（次数 10^5、10^5 点）" << t_mp << " ms, interpolate（10^5 点）" << t_ip << " ms  (check "
              << (gi[5] + ge[7] + gp[9] + cf + ev[77] + ip[55] + (gs.ok ? gs.value[3] : mint(0))).val() % 7 << (ge == Fps<mint>(f) ? "" : " exp(log f) ≠ f") << ")\n";
}

}  // namespace fps_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace fps_demo {

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
#define SHOW(expr, comment) fps_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    fps_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("形式的冪級数（fps） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Fps<mint>：+ − × / f.pre(n) / f.eval(x)\n"
              << "  fps_inv / fps_log / fps_exp / fps_pow / fps_sqrt（最初の n 項）\n"
              << "  fps_divmod / fps_derivative / fps_integral / taylor_shift\n"
              << "  fps_coeff_of_fraction(P, Q, K) / prod_linear(as)\n";
    note("using mint = ModInt998;");
    using mint = ModInt998;

    api("【1】 Fps<mint> f = {…} / f * g / fps_inv(f, n)", "Fps<mint>", "係数は低い次数から。× は畳み込み。1/f の最初の n 項（f[0] ≠ 0）", "O(n log n)");
    RUN(Fps<mint> f = {1, -1, -1});
    SHOW(fps_inv(f, 10), "1/(1 − x − x^2) = フィボナッチ数 1, 1, 2, 3, 5, …");
    SHOW(f * f, "(1 − x − x^2)^2");
    {
        // 分割数 p(n)：1/Π(1 − x^k) の係数
        const int N = 101;
        Fps<mint> e(N);
        e[0] = 1;
        for (long long j = 1;; ++j) {  // Π(1 − x^k) = 1 + Σ_{j≥1} (−1)^j (x^{j(3j−1)/2} + x^{j(3j+1)/2})（五角数定理）
            const long long a = j * (3 * j - 1) / 2, b = j * (3 * j + 1) / 2;
            if (a >= N) break;
            e[a] += (j % 2 ? mint(-1) : mint(1));
            if (b < N) e[b] += (j % 2 ? mint(-1) : mint(1));
        }
        const Fps<mint> p = fps_inv(e, N);
        note("分割数 p(n) = [x^n] 1/Π(1 − x^k)（Π は五角数定理で作る）");
        std::cout << "    p(10) = " << p[10] << "、p(100) = " << p[100] << "\n";
    }

    api("【2】 fps_log(f, n) / fps_exp(f, n) / fps_pow(f, k, n) / fps_sqrt(f, n)", "Fps<mint> / Maybe<Fps<mint>>",
        "log f（f[0] = 1）/ exp f（f[0] = 0）/ f^k（k ≤ 10^18）/ √f（ないときは ok = false）", "O(n log n)");
    {
        // ベル数：exp(e^x − 1) の係数 × n!
        const int N = 11;
        Fps<mint> ex(N);
        mint fact = 1;
        for (int i = 1; i < N; ++i) fact *= mint(i), ex[i] = fact.inv();  // e^x − 1
        const Fps<mint> b = fps_exp(ex, N);
        mint f10 = 1;
        for (int i = 1; i <= 10; ++i) f10 *= mint(i);
        note("ベル数 B_n = n! [x^n] exp(e^x − 1)");
        std::cout << "    B_10 = " << b[10] * f10 << "（115975）\n";
    }
    SHOW(fps_pow<mint>({1, 1}, 5, 7), "(1 + x)^5");
    SHOW(fps_pow<mint>({0, 1, 1}, 1000000000000000000LL, 3), "(x + x^2)^{10^18} は x^{10^18} から始まるので最初の 3 項は 0");
    {
        // カタラン数：C(x) = (1 − √(1 − 4x)) / (2x)
        const Maybe<Fps<mint>> s = fps_sqrt<mint>({1, -4}, 12);
        Fps<mint> c(11);
        for (int i = 0; i < 11; ++i) c[i] = -s.value[i + 1] / 2;
        note("カタラン数 = (1 − √(1 − 4x)) / (2x) の係数（fps_sqrt）");
        std::cout << "    " << to_s(std::vector<mint>(c.begin(), c.end())) << "\n";
    }
    SHOW(fps_sqrt<mint>({0, 1}, 3).ok, "x の平方根は冪級数にならない");

    api("【3】 fps_divmod(f, g) / fps_derivative(f) / fps_integral(f) / taylor_shift(f, c)", "FpsDivMod{q, r} / Fps<mint>",
        "多項式の割り算 f = q g + r / 微分 / 積分 / f(x + c)", "O(n log n)");
    RUN(FpsDivMod<mint> d = fps_divmod<mint>({5, 0, 0, 1}, {1, 1}));
    SHOW(d.q, "(x^3 + 5) ÷ (x + 1) の商 x^2 − x + 1（−1 は mod で 998244352）");
    SHOW(d.r, "余り 4");
    SHOW(fps_derivative<mint>({1, 2, 3}), "1 + 2x + 3x^2 の微分");
    SHOW(taylor_shift<mint>({0, 0, 1}, mint(1)), "x^2 → (x + 1)^2 = 1 + 2x + x^2");

    api("【4】 fps_coeff_of_fraction(P, Q, K) / prod_linear(as)", "mint / Fps<mint>",
        "[x^K] P/Q（線形漸化式の第 K 項。K ≤ 10^18）/ Π (1 + a_i x)", "O(d log d log K) / O(n log^2 n)");
    SHOW(fps_coeff_of_fraction<mint>({0, 1}, {1, -1, -1}, 10), "フィボナッチ F_10 = 55（x/(1 − x − x^2)）");
    SHOW(fps_coeff_of_fraction<mint>({0, 1}, {1, -1, -1}, 1000000000000000000LL), "F_{10^18} mod 998244353");
    SHOW(prod_linear<mint>({1, 2, 3}), "(1 + x)(1 + 2x)(1 + 3x)：基本対称式");
    note("a_i を 1 つ以上選んだ積の和など、k 個選ぶ積の和は prod_linear の k 次の係数");

    api("【5】 multipoint_eval(f, xs) / interpolate(xs, ys)", "vector<mint> / Fps<mint>",
        "f(x_0), …, f(x_{m−1}) をまとめて / 点 (x_i, y_i)（x は異なる）を通る次数 n − 1 以下の多項式", "O(n log^2 n)");
    SHOW(multipoint_eval<mint>({1, 0, 1}, {0, 1, 2, 3}), "1 + x^2 の x = 0, 1, 2, 3 での値");
    SHOW(interpolate<mint>({0, 1, 2}, {1, 2, 5}), "(0,1) (1,2) (2,5) を通る：1 + x^2");
    note("x が 0, 1, …, n−1 の等間隔で 1 点の値だけ欲しいなら 06_series の lagrange_interpolate（O(n)）が速い");

    api("【6】 power_sums(as, n, ws) / fps_sum_of_fractions(nums, dens)", "Fps<mint> / FpsFraction<mint>{num, den}",
        "p_k = Σ w_i a_i^k（k = 0〜n−1。ws を省くと 1）/ Σ P_i/Q_i を 1 つの分数に（分割統治で通分）", "O(m log^2 m + n log n)");
    SHOW(power_sums<mint>({1, 2, 3}, 4), "Σ a^k：3, 6, 14, 36");
    SHOW(power_sums<mint>({2, 3}, 3, {10, 1}), "10·2^k + 3^k：11, 23, 49");
    note("Σ_i a_i^k をすべての k = 0〜n−1 について一度に（1/(1 − a_i x) の和の係数）");
}

}  // namespace fps_demo

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
    fps_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    fps_test::benchmark();    // 速度計測
    fps_demo::demo();         // 使い方デモ
    return 0;
}
