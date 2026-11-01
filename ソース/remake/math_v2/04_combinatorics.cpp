// =============================================================================
//  数え上げ（スターリング数・ベル数・分割数・完全順列・全射・写像 12 相・包除原理・二項変換・組分けの列挙・q 二項係数・ヤング盤）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 04 組合せ
//  依存       : 00_int_common, 00_math_base, 00_modint, 04_binomial, 07_convolution（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  使い方
//      using mint = ModInt998;                               // mod は素数（10^9+7 でもよい）
//      mint s = stirling2<mint>(10, 3);                      // S(10, 3)
//      mint w = balls_into_boxes<mint>(5, 3, true, false, 2); // 区別する玉 5 個を区別しない箱 3 個へ、どの箱も 1 個以上
//
//  機能一覧（T は mint でも long long でもよい表（long long は溢れたら止める）。mint だけのものは <mint>）   計算量
//    stirling2<mint>(n, k)            第 2 種スターリング数 S(n, k)：n 個を k 個の空でない組に分ける数（n ≤ 10^18）  O(k log n)
//    stirling2_row<mint>(n)           S(n, 0), …, S(n, n)                                                 O(n log n)
//    stirling2_table<T>(N)            t[n][k]（0 ≤ k ≤ n ≤ N）                                            O(N^2)
//    stirling1_row<mint>(n)           第 1 種（符号なし）c(n, k)：サイクルが k 個の順列の数。k = 0〜n         O(n log^2 n)
//    stirling1_table<T>(N)            t[n][k]（符号つきは (−1)^{n−k} 倍）                                 O(N^2)
//    bell<mint>(n) / bell_table<T>(N) ベル数 B_n（n 個を空でない組にいくつかに分ける数）                     O(n log n) / O(N^2)
//    partition_row<mint>(N)           分割数 p(0), …, p(N)（n を正の整数の和に。順番は区別しない）            O(N log N)
//    partition_table<T>(N)            同じく五角数定理で（long long でも使える）                             O(N √N)
//    partition_count_table<T>(N, K)   t[k][n]：n を k 個以下に分ける数（k ≤ K, n ≤ N）                     O(N K)
//    derangement_table<T>(N)          完全順列（どの i も i 番目に来ない順列）の数 D_0, …, D_N               O(N)
//    surjections<mint>(n, k)          n 個から k 個への全射の数 k! S(n, k)                                  O(k log n)
//    balls_into_boxes<mint>(n, k, 玉を区別, 箱を区別, rule)                                                写像 12 相を全部
//                                     rule：0 = 何個でもよい、1 = 各箱 1 個以下、2 = 各箱 1 個以上
//                                     （玉も箱も区別しないときは分割数：k ≥ n なら n ≤ 10^7、k < n なら n·k ≤ 3×10^8）
//    inclusion_exclusion<T>(n, f)     Σ_{S ⊆ {0..n−1}} (−1)^{|S|} f(S)（n ≤ 30。S はビット）                       O(2^n)
//                                     f(S) =「条件 S を全部満たすものの数」を渡すと「どの条件も満たさないものの数」になる
//                                     （「少なくとも 1 つ満たす」数は 全体 − これ）
//    inclusion_exclusion_count<mint>(n, g)  Σ_i (−1)^i C(n, i) g(n − i)（g(m)：m 個だけ使ってよいときの数 → n 個全部使う数）  O(n)
//    count_divisible_by_any(l, r, ds) [l, r) で ds のどれかで割り切れる数の個数（閉区間は _closed）        O(2^|ds|) 以下
//    binomial_transform<mint>(g) / binomial_inverse<mint>(f)    f(n) = Σ_k C(n, k) g(k) / その逆（(−1)^{n−k} を付ける）  O(n log n)
//    at_least_from_exactly<mint>(E) / exactly_from_at_least<mint>(A)
//                                     A[k] = Σ_{j≥k} C(j, k) E[j]（k 個を選んで満たす、の和）/ その逆（ちょうど k 個 E を A から）
//    for_each_set_partition(n, f)     {0..n−1} の組分けを全部 f(g, k)（g[i] = 組の番号、k = 組の数。n ≤ 12 程度）  O(B_n · n)
//    for_each_integer_partition(n, f, max_part)   n の分割を全部 f(parts)（大きい順。各部分 ≤ max_part）     O(p(n) · n)
//    q_binomial<mint>(n, k, q)        q 二項係数（F_q^n の k 次元部分空間の数。q の位数が小さくても q-Lucas で。n ≤ 10^7）  O(n)
//    young_tableaux_count<mint>(shape)  形 shape（行の長さ。大きい順）の標準ヤング盤の数（フック長公式）       O(マスの数)
//  （nCr・カタラン数・格子路は 04_binomial、順列の番号・転倒数は 04_permutation）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    comb_test::random_test()  … 全部並べて数える愚直解・漸化式との突き合わせ
//    comb_test::benchmark()    … 速度計測
//    comb_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_COMBINATORICS
#define GG_COMBINATORICS
// ---- 数え上げ（数学第二世代）  スターリング数・ベル数・分割数・完全順列・全射・写像 12 相・包除原理・二項変換・組分けの列挙・q 二項係数・ヤング盤 ----
// mint の mod は素数（Binomial を使う）。_table は O(N^2) の DP、_row は畳み込みで O(n log n)
namespace math_internal {
template <class mint>
Binomial<mint>& comb_binom() {  // ファイルの中で使い回す階乗の表（必要なら自動で伸びる）
    static Binomial<mint> B(0);
    return B;
}
// T が整数型なら溢れを調べる足し算・掛け算（long long の表で黙って溢れないように）
template <class T>
T cb_add(const T& a, const T& b, const char* f) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(f, "値が整数型に収まらない（mint で使う）");
        return r;
    } else {
        return a + b;
    }
}
template <class T>
T cb_mul(const T& a, const T& b, const char* f) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_mul_overflow(a, b, &r)) fail(f, "値が整数型に収まらない（mint で使う）");
        return r;
    } else {
        return a * b;
    }
}
template <class T>
using CbAcc = typename std::conditional<std::is_integral<T>::value, Int128, T>::type;  // 符号が交互の和は Int128 で足してから戻す
template <class T>
T cb_from_acc(const CbAcc<T>& v, const char* f) {
    if constexpr (std::is_integral<T>::value) {
        if (v < (Int128)std::numeric_limits<T>::min() || v > (Int128)std::numeric_limits<T>::max()) fail(f, "値が整数型に収まらない（mint で使う）");
        return (T)v;
    } else {
        return v;
    }
}
inline void cb_check_n(long long N, const char* f) {
    if (N < 0) fail(f, "N ≥ 0 で使う");
}
template <class mint>
std::vector<mint> comb_series_inv(const std::vector<mint>& Q, size_t n) {  // 1/Q の最初の n 項（Q[0] = 1。掛け算だけのニュートン法）
    std::vector<mint> g{mint(1)};
    for (size_t k = 1; k < n;) {
        k *= 2;
        std::vector<mint> q(Q.begin(), Q.begin() + (long long)std::min(Q.size(), k));
        std::vector<mint> t = convolve(q, g);  // g ← g (2 − Q g)
        t.resize(k);
        for (auto& x : t) x = -x;
        t[0] += mint(2);
        g = convolve(g, t);
        g.resize(k);
    }
    g.resize(n);
    return g;
}
template <class mint>
mint comb_falling(long long a, long long n) {  // a (a−1) … (a−n+1)
    mint r = 1;
    for (long long i = 0; i < n; ++i) r *= mint(a - i);
    return r;
}
}  // namespace math_internal

// ---- 第 2 種スターリング数 S(n, k)：n 個の区別できるものを k 個の空でない（区別しない）組に分ける数 ----
template <class mint>
mint stirling2(long long n, long long k) {  // 1 つだけ：Σ_{i=0}^{k} (−1)^{k−i} C(k, i) i^n / k!  O(k log n)
    if (n < 0 || k < 0 || k > n) return mint(0);
    auto& B = math_internal::comb_binom<mint>();
    mint s = 0;
    for (long long i = 0; i <= k; ++i) {
        const mint t = B.C(k, i) * mint(i).pow(n);
        s += ((k - i) % 2 ? -t : t);
    }
    return s * B.inv_fact(k);
}
template <class mint>
std::vector<mint> stirling2_row(int n) {  // S(n, 0), …, S(n, n)（畳み込み 1 回）  O(n log n)
    math_internal::cb_check_n(n, "stirling2_row");
    auto& B = math_internal::comb_binom<mint>();
    std::vector<mint> a((size_t)n + 1), b((size_t)n + 1);
    for (int i = 0; i <= n; ++i) a[(size_t)i] = mint(i).pow(n) * B.inv_fact(i), b[(size_t)i] = (i % 2 ? -B.inv_fact(i) : B.inv_fact(i));
    std::vector<mint> c = convolve(a, b);
    c.resize((size_t)n + 1);
    return c;
}
template <class T>
std::vector<std::vector<T>> stirling2_table(int N) {  // t[n][k]（0 ≤ k ≤ n ≤ N。S(n, k) = k S(n−1, k) + S(n−1, k−1)。整数型は溢れたら止める）  O(N^2)
    math_internal::cb_check_n(N, "stirling2_table");
    std::vector<std::vector<T>> t((size_t)N + 1);
    for (int n = 0; n <= N; ++n) {
        t[(size_t)n].assign((size_t)n + 1, T(0));
        t[(size_t)n][0] = (n == 0 ? T(1) : T(0));
        for (int k = 1; k <= n; ++k)
            t[(size_t)n][(size_t)k] = math_internal::cb_add(k < n ? math_internal::cb_mul(t[(size_t)n - 1][(size_t)k], T(k), "stirling2_table") : T(0), t[(size_t)n - 1][(size_t)k - 1], "stirling2_table");
    }
    return t;
}
// ---- 第 1 種スターリング数（符号なし）c(n, k)：n 個の順列でサイクルがちょうど k 個のものの数（符号つきは (−1)^{n−k} 倍） ----
template <class mint>
std::vector<mint> stirling1_row(int n) {  // c(n, 0), …, c(n, n) = x(x+1)…(x+n−1) の係数（分割統治）  O(n log^2 n)
    math_internal::cb_check_n(n, "stirling1_row");
    std::vector<std::vector<mint>> ps;
    for (int i = 0; i < n; ++i) ps.push_back({mint(i), mint(1)});
    return convolve_many(std::move(ps));
}
template <class T>
std::vector<std::vector<T>> stirling1_table(int N) {  // t[n][k] = c(n, k)（c(n, k) = (n−1) c(n−1, k) + c(n−1, k−1)）  O(N^2)
    math_internal::cb_check_n(N, "stirling1_table");
    std::vector<std::vector<T>> t((size_t)N + 1);
    for (int n = 0; n <= N; ++n) {
        t[(size_t)n].assign((size_t)n + 1, T(0));
        t[(size_t)n][0] = (n == 0 ? T(1) : T(0));
        for (int k = 1; k <= n; ++k)
            t[(size_t)n][(size_t)k] = math_internal::cb_add(k < n ? math_internal::cb_mul(t[(size_t)n - 1][(size_t)k], T(n - 1), "stirling1_table") : T(0), t[(size_t)n - 1][(size_t)k - 1], "stirling1_table");
    }
    return t;
}
// ---- ベル数 B_n（n 個を空でない組にいくつかに分ける数）----
template <class mint>
mint bell(int n) {  // B_n = Σ_k S(n, k)  O(n log n)
    mint s = 0;
    for (const mint& x : stirling2_row<mint>(n)) s += x;
    return s;
}
template <class T>
std::vector<T> bell_table(int N) {  // B_0, …, B_N（ベルの三角形。n − 1 行目の最後が B_n）  O(N^2)
    math_internal::cb_check_n(N, "bell_table");
    std::vector<T> res((size_t)N + 1), row{T(1)};
    res[0] = T(1);
    for (int n = 1; n <= N; ++n) {
        res[(size_t)n] = row.back();
        if (n == N) break;
        std::vector<T> nxt{row.back()};
        for (const T& x : row) nxt.push_back(math_internal::cb_add(nxt.back(), x, "bell_table"));
        row.swap(nxt);
    }
    return res;
}
// ---- 分割数 ----
template <class T>
std::vector<T> partition_table(int N) {  // p(0), …, p(N)：n を正の整数の和に分ける数（順番は区別しない。五角数定理）  O(N √N)
    math_internal::cb_check_n(N, "partition_table");
    std::vector<T> p((size_t)N + 1, T(0));
    p[0] = T(1);
    for (int n = 1; n <= N; ++n) {
        math_internal::CbAcc<T> s = 0;
        for (long long j = 1;; ++j) {
            const long long a = j * (3 * j - 1) / 2, b = j * (3 * j + 1) / 2;
            if (a > n) break;
            const bool plus = j % 2 == 1;
            s = plus ? s + p[(size_t)(n - a)] : s - p[(size_t)(n - a)];
            if (b <= n) s = plus ? s + p[(size_t)(n - b)] : s - p[(size_t)(n - b)];
        }
        p[(size_t)n] = math_internal::cb_from_acc<T>(s, "partition_table");
    }
    return p;
}
template <class mint>
std::vector<mint> partition_row(int N) {  // p(0), …, p(N) を 1/Π(1 − x^k)（五角数の級数の逆数）で  O(N log N)
    math_internal::cb_check_n(N, "partition_row");
    std::vector<mint> Q((size_t)N + 1, mint(0));  // Π(1 − x^k) = Σ_j (−1)^j x^{j(3j−1)/2}（j は ±）
    Q[0] = 1;
    for (long long j = 1;; ++j) {
        const long long a = j * (3 * j - 1) / 2, b = j * (3 * j + 1) / 2;
        if (a > N) break;
        const mint sg = j % 2 ? mint(-1) : mint(1);
        Q[(size_t)a] += sg;
        if (b <= N) Q[(size_t)b] += sg;
    }
    return math_internal::comb_series_inv(Q, (size_t)N + 1);
}
template <class T>
std::vector<std::vector<T>> partition_count_table(int N, int K) {  // t[k][n]：n を k 個以下に分ける数（= k 以下の数の和にする数）  O(N K)
    math_internal::cb_check_n(N, "partition_count_table"), math_internal::cb_check_n(K, "partition_count_table");
    std::vector<std::vector<T>> t((size_t)K + 1, std::vector<T>((size_t)N + 1, T(0)));
    t[0][0] = T(1);
    for (int k = 1; k <= K; ++k)
        for (int n = 0; n <= N; ++n) t[(size_t)k][(size_t)n] = math_internal::cb_add(t[(size_t)k - 1][(size_t)n], n >= k ? t[(size_t)k][(size_t)(n - k)] : T(0), "partition_count_table");
    return t;
}
// ---- 完全順列・全射 ----
template <class T>
std::vector<T> derangement_table(int N) {  // D_0, …, D_N：どの i も i 番目に来ない順列の数（D_n = (n−1)(D_{n−1} + D_{n−2})）
    math_internal::cb_check_n(N, "derangement_table");
    std::vector<T> d((size_t)std::max(N + 1, 2), T(0));
    d[0] = T(1), d[1] = T(0);
    for (int n = 2; n <= N; ++n)
        d[(size_t)n] = math_internal::cb_mul(T(n - 1), math_internal::cb_add(d[(size_t)n - 1], d[(size_t)n - 2], "derangement_table"), "derangement_table");
    d.resize((size_t)N + 1);
    return d;
}
template <class mint>
mint surjections(long long n, long long k) {  // n 個から k 個への全射の数 k! S(n, k)
    if (n < 0 || k < 0 || k > n) return mint(0);
    return math_internal::comb_binom<mint>().fact(k) * stirling2<mint>(n, k);
}

// ---- 写像 12 相：n 個の玉を k 個の箱に入れる方法の数 ----
// balls_distinct / boxes_distinct：区別するか。rule：0 = 何個でもよい、1 = 各箱に 1 個以下、2 = 各箱に 1 個以上
template <class mint>
mint balls_into_boxes(long long n, long long k, bool balls_distinct, bool boxes_distinct, int rule) {
    if (n < 0 || k < 0 || rule < 0 || rule > 2) math_internal::fail("balls_into_boxes", "n, k ≥ 0、rule は 0・1・2");
    auto& B = math_internal::comb_binom<mint>();
    if (rule == 1 && n > k) return mint(0);
    if (rule == 2 && n < k) return mint(0);
    const bool big = n + k > 20000000;  // k が巨大なら階乗の表を作らず、n 個の積で
    if (balls_distinct && boxes_distinct) {
        if (rule == 0) return mint(k).pow(n);
        if (rule == 1) return big ? math_internal::comb_falling<mint>(k, n) : B.P(k, n);
        return surjections<mint>(n, k);
    }
    if (!balls_distinct && boxes_distinct) {
        if (rule == 0) {
            if (n == 0) return mint(1);
            return big ? math_internal::comb_falling<mint>(n + k - 1, n) * B.inv_fact(n) : B.H(k, n);
        }
        if (rule == 1) return big ? math_internal::comb_falling<mint>(k, n) * B.inv_fact(n) : B.C(k, n);
        return n == 0 ? mint(k == 0 ? 1 : 0) : B.C(n - 1, k - 1);
    }
    if (balls_distinct && !boxes_distinct) {
        if (rule == 1) return mint(1);
        if (rule == 2) return stirling2<mint>(n, k);
        // Σ_{i ≤ k} S(n, i) = Σ_j j^n / j! · Σ_{t ≤ k−j} (−1)^t / t!  O(k log n)
        const long long kk = std::min(n, k);
        std::vector<mint> e((size_t)kk + 1);  // e[t] = Σ_{s ≤ t} (−1)^s / s!
        for (long long t = 0; t <= kk; ++t) e[(size_t)t] = (t ? e[(size_t)t - 1] : mint(0)) + (t % 2 ? -B.inv_fact(t) : B.inv_fact(t));
        mint s = 0;
        for (long long j = 0; j <= kk; ++j) s += mint(j).pow(n) * B.inv_fact(j) * e[(size_t)(kk - j)];
        return s;
    }
    // 玉も箱も区別しない：分割数（m を k 個以下に分ける数）
    if (rule == 1) return mint(1);
    const long long m = rule == 2 ? n - k : n;  // 1 個以上なら、先に 1 個ずつ入れて残りを k 個以下に
    if (k == 0) return mint(m == 0 ? 1 : 0);
    if (k >= m) {  // 制限なし：分割数 p(m)
        if (m > 10000000) math_internal::fail("balls_into_boxes", "玉も箱も区別しない場合は n ≤ 10^7 程度");
        return partition_row<mint>((int)m)[(size_t)m];
    }
    if ((Int128)(m + 1) * k > 300000000) math_internal::fail("balls_into_boxes", "玉も箱も区別せず k < n の場合は n·k ≤ 3×10^8 程度");
    std::vector<mint> dp((size_t)m + 1, mint(0));  // 1〜k の数を何個でも使って和を作る
    dp[0] = 1;
    for (long long j = 1; j <= k; ++j)
        for (long long v = j; v <= m; ++v) dp[(size_t)v] += dp[(size_t)(v - j)];
    return dp[(size_t)m];
}

// ---- 包除原理 ----
// Σ_{S ⊆ {0, …, n−1}} (−1)^{|S|} f(S)（n ≤ 30。f は集合 S（ビット）を受け取る。「どの条件も満たさない数」を「条件 S を全部満たす数」から）  O(2^n)
template <class T, class F>
T inclusion_exclusion(int n, F f) {
    if (n < 0 || n > 30) math_internal::fail("inclusion_exclusion", "0 ≤ n ≤ 30 で使う");
    math_internal::CbAcc<T> s = 0;  // 整数型なら Int128 で足して、最後に収まるか確かめる
    for (long long S = 0; S < (1LL << n); ++S) {
        const T v = f(S);
        s = (__builtin_popcountll((unsigned long long)S) % 2) ? s - v : s + v;
    }
    return math_internal::cb_from_acc<T>(s, "inclusion_exclusion");
}
// 対称な包除原理：g(m) = 「n 個の候補のうち、ある m 個だけを使ってよい」ときの数 → n 個を全部使う数 = Σ_{i=0}^{n} (−1)^i C(n, i) g(n − i)  O(n)
template <class mint, class G>
mint inclusion_exclusion_count(int n, G g) {
    math_internal::cb_check_n(n, "inclusion_exclusion_count");
    auto& B = math_internal::comb_binom<mint>();
    mint s = 0;
    for (int i = 0; i <= n; ++i) {
        const mint t = B.C(n, i) * mint(g(n - i));
        s += (i % 2 ? -t : t);
    }
    return s;
}
// ---- 二項変換・「ちょうど k 個」と「k 個を選んで満たす」の変換（mint。畳み込み 1 回）----
namespace math_internal {
template <class mint>
std::vector<mint> cb_egf_mul(const std::vector<mint>& g, bool alt) {  // f(n) = Σ_k C(n, k) g(k) (−1)^{(n−k)·alt}
    const size_t n = g.size();
    if (n == 0) return {};
    auto& B = comb_binom<mint>();
    std::vector<mint> a(n), b(n);
    for (size_t k = 0; k < n; ++k) a[k] = g[k] * B.inv_fact((long long)k), b[k] = (alt && k % 2) ? -B.inv_fact((long long)k) : B.inv_fact((long long)k);
    std::vector<mint> c = convolve(a, b);
    c.resize(n);
    for (size_t k = 0; k < n; ++k) c[k] *= B.fact((long long)k);
    return c;
}
template <class mint>
std::vector<mint> cb_upper(const std::vector<mint>& E, bool alt) {  // A[k] = Σ_{j≥k} C(j, k) E[j] (−1)^{(j−k)·alt}
    const size_t n = E.size();
    if (n == 0) return {};
    auto& B = comb_binom<mint>();
    std::vector<mint> R(n), b(n);
    for (size_t i = 0; i < n; ++i) R[i] = B.fact((long long)(n - 1 - i)) * E[n - 1 - i], b[i] = (alt && i % 2) ? -B.inv_fact((long long)i) : B.inv_fact((long long)i);
    const std::vector<mint> c = convolve(R, b);
    std::vector<mint> A(n);
    for (size_t k = 0; k < n; ++k) A[k] = c[n - 1 - k] * B.inv_fact((long long)k);
    return A;
}
}  // namespace math_internal
// 二項変換 f(n) = Σ_{k=0}^{n} C(n, k) g(k)（n = 0〜|g|−1）  O(n log n)
template <class mint>
std::vector<mint> binomial_transform(const std::vector<mint>& g) { return math_internal::cb_egf_mul(g, false); }
// その逆 g(n) = Σ_{k=0}^{n} (−1)^{n−k} C(n, k) f(k)（「n 個の候補のうち使ってよいもの」f から「n 個全部使う」g へ。包除原理）
template <class mint>
std::vector<mint> binomial_inverse(const std::vector<mint>& f) { return math_internal::cb_egf_mul(f, true); }
// E[j] = ちょうど j 個の条件を満たすものの数 → A[k] = Σ_{j≥k} C(j, k) E[j]（「k 個の条件を選び、それを全部満たす」の和）
template <class mint>
std::vector<mint> at_least_from_exactly(const std::vector<mint>& E) { return math_internal::cb_upper(E, false); }
// その逆：A[k]（k 個の条件を選んで全部満たす、の和。数えやすい）→ E[k] = Σ_{j≥k} (−1)^{j−k} C(j, k) A[j]（ちょうど k 個）
template <class mint>
std::vector<mint> exactly_from_at_least(const std::vector<mint>& A) { return math_internal::cb_upper(A, true); }

// [l, r) の整数で、ds のどれかで割り切れるものの数（ds ≥ 1。包除原理。lcm が範囲を超えた組はまとめて打ち切る）  O(2^|ds|) 以下
inline long long count_divisible_by_any(long long l, long long r, std::vector<long long> ds) {
    for (long long d : ds)
        if (d < 1) math_internal::fail("count_divisible_by_any", "ds ≥ 1 で使う");
    if (r <= l || ds.empty()) return 0;
    // 他の数の倍数になっている数はいらない（その倍数は全部もう数えている）
    std::sort(ds.begin(), ds.end());
    ds.erase(std::unique(ds.begin(), ds.end()), ds.end());
    {
        std::vector<long long> keep;
        for (long long d : ds) {
            bool red = false;
            for (long long e : keep)
                if (d % e == 0) { red = true; break; }
            if (!red) keep.push_back(d);
        }
        ds.swap(keep);
    }
    using U = unsigned long long;
    const U ub = std::max(math_internal::uabs(l), math_internal::uabs(r));  // lcm がこれより大きいと、範囲にある倍数は 0 だけ
    const int n = (int)ds.size();
    const long long zero = (l <= 0 && 0 < r) ? 1 : 0;
    Int128 total = 0;
    auto dfs = [&](auto&& self, int i, U lc, int cnt) -> void {
        if (i == n) {
            if (cnt > 0) {
                const Int128 c = math_internal::ceil_div128(r, (Int128)lc) - math_internal::ceil_div128(l, (Int128)lc);
                total += (cnt % 2 ? c : -c);
            }
            return;
        }
        self(self, i + 1, lc, cnt);
        const U d = (U)ds[(size_t)i], g = math_internal::ugcd(lc, d);
        const UInt128 nl = (UInt128)(lc / g) * d;
        if (nl > ub) {  // この先の組はみな「0 だけ」なので、残りを足し引きすると打ち消し合う（最後の 1 個のときだけ残る）
            if (i + 1 == n) total += ((cnt + 1) % 2 ? zero : -zero);
            return;
        }
        self(self, i + 1, (U)nl, cnt + 1);
    };
    dfs(dfs, 0, 1, 0);
    if (total > (Int128)LLONG_MAX) math_internal::fail("count_divisible_by_any", "答えが long long に収まらない");
    return (long long)total;
}
inline long long count_divisible_by_any_closed(long long L, long long R, const std::vector<long long>& ds) {  // [L, R]
    if (R < L) return 0;
    if (R < LLONG_MAX) return count_divisible_by_any(L, R + 1, ds);
    long long c = count_divisible_by_any(L, R, ds);
    for (long long d : ds)
        if (d >= 1 && R % d == 0) return c + 1;
    return c;
}

// ---- 分け方を全部たどる ----
// {0, …, n−1} を空でない組に分ける方法をすべて f(g, k) に渡す（g[i] = i の組の番号（初めて出る順に 0, 1, …）、k = 組の数）
//   ベル数 B_n 通り（n = 10 で 11 万、12 で 420 万、13 で 2760 万）  O(B_n · n)
template <class F>
void for_each_set_partition(int n, F f) {
    if (n < 0 || n > 15) math_internal::fail("for_each_set_partition", "0 ≤ n ≤ 15 で使う（B_15 は 13 億）");
    std::vector<int> g((size_t)n, 0);
    auto rec = [&](auto&& self, int i, int k) -> void {
        if (i == n) {
            f((const std::vector<int>&)g, k);
            return;
        }
        for (int b = 0; b <= k; ++b) g[(size_t)i] = b, self(self, i + 1, b == k ? k + 1 : k);
    };
    rec(rec, 0, 0);
}
// n を正の整数の和に分ける方法（順番は区別しない）をすべて f(parts) に渡す（parts は大きい順。max_part を渡すと各部分は max_part 以下）
//   分割数 p(n) 通り（n = 50 で 20 万、70 で 400 万）  O(p(n) · n)
template <class F>
void for_each_integer_partition(int n, F f, int max_part = -1) {
    if (n < 0) math_internal::fail("for_each_integer_partition", "n ≥ 0 で使う");
    if (max_part < 0 || max_part > n) max_part = n;
    std::vector<int> parts;
    if (n == 0) {
        f((const std::vector<int>&)parts);
        return;
    }
    struct Frame {
        int rem, x;  // 残り rem を、次は x 以下の部分から試す
    };
    std::vector<Frame> st;  // 再帰しない（部分の数が多くても深くならない）
    st.push_back(Frame{n, std::min(n, max_part)});
    while (!st.empty()) {
        Frame& t = st.back();
        if (t.x == 0) {  // この段は全部試した：ひとつ前の段で選んだ部分を戻す
            st.pop_back();
            if (!parts.empty()) parts.pop_back();
            continue;
        }
        const int x = t.x--, rem = t.rem - x;
        parts.push_back(x);
        if (rem == 0) {
            f((const std::vector<int>&)parts);
            parts.pop_back();
            continue;
        }
        st.push_back({rem, std::min(rem, x)});
    }
}

// ---- q 二項係数・ヤング盤 ----
// q 二項係数 [n, k]_q = Π_{i=0}^{k−1} (1 − q^{n−i}) / (1 − q^{i+1})（q が素数べきなら F_q^n の k 次元部分空間の数）
//   q は何でも（q^d ≡ 1 となる d があれば q-Lucas の定理で）。n ≤ 10^7  O(n)
template <class mint>
mint q_binomial(long long n, long long k, long long q) {
    if (n < 0 || k < 0 || k > n) return mint(0);
    if (n > 10000000) math_internal::fail("q_binomial", "n ≤ 10^7 で使う");
    k = std::min(k, n - k);
    const mint Q(q);
    if (Q == mint(0)) return mint(1);  // [n, k]_0 = 1
    long long d = 0;                   // q の位数（q^d = 1 となる最小の d ≤ n。なければ 0）
    {
        mint t = 1;
        for (long long j = 1; j <= n; ++j)
            if ((t *= Q) == mint(1)) {
                d = j;
                break;
            }
    }
    auto direct = [&](long long a, long long b) {  // [a, b]_q（b ≤ a で分母が 0 にならないとき）
        mint num = 1, den = 1, x = Q.pow(a - b + 1), y = Q;
        for (long long i = 0; i < b; ++i) num *= mint(1) - x, den *= mint(1) - y, x *= Q, y *= Q;
        return num / den;
    };
    if (d == 0) return direct(n, k);
    // q-Lucas：[n, k]_q ≡ C(⌊n/d⌋, ⌊k/d⌋) · [n mod d, k mod d]_q。C は mod が小さいときは Lucas の定理で
    auto& B = math_internal::comb_binom<mint>();
    const long long p = (long long)mint::mod();
    mint c = 1;
    for (long long a = n / d, b = k / d; (a > 0 || b > 0) && c != mint(0); a /= p, b /= p) c *= B.C(a % p, b % p);
    return c == mint(0) ? c : c * direct(n % d, k % d);
}
// 形 shape（各行の長さ。大きい順）の標準ヤング盤の数 = n! / Π（各マスのフックの長さ）（フック長公式。n はマスの数）  O(n)
template <class mint, class V>
mint young_tableaux_count(const V& shape) {
    std::vector<long long> row;
    for (const auto& x : shape) row.push_back(math_internal::to_ll_checked(x, "young_tableaux_count"));
    while (!row.empty() && row.back() == 0) row.pop_back();
    long long n = 0;
    for (size_t i = 0; i < row.size(); ++i) {
        if (row[i] < 0 || (i > 0 && row[i] > row[i - 1])) math_internal::fail("young_tableaux_count", "行の長さは 0 以上で大きい順に並べる");
        if (row[i] > 100000000 - n) math_internal::fail("young_tableaux_count", "マスが多すぎる（10^8 まで）");
        n += row[i];
    }
    std::vector<long long> col(row.empty() ? 0 : (size_t)row[0], 0);  // col[j] = j 列目の長さ
    for (long long r : row)
        for (long long j = 0; j < r; ++j) ++col[(size_t)j];
    auto& B = math_internal::comb_binom<mint>();
    mint prod = 1;
    for (size_t i = 0; i < row.size(); ++i)
        for (long long j = 0; j < row[i]; ++j) prod *= mint((row[i] - j - 1) + (col[(size_t)j] - (long long)i - 1) + 1);
    return B.fact(n) / prod;
}
template <class mint>
mint young_tableaux_count(std::initializer_list<long long> shape) { return young_tableaux_count<mint, std::initializer_list<long long>>(shape); }
#endif  // GG_COMBINATORICS
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace comb_test {

#define COMB_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[combinatorics] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 全部の写像 [n] → [k] を並べて、玉・箱の区別に合わせて同じものをまとめて数える
static long long brute_twelve(int n, int k, bool bd, bool xd, int rule) {
    std::set<std::vector<int>> seen;
    std::vector<int> f((size_t)n, 0);
    long long total = 1;
    for (int i = 0; i < n; ++i) total *= k;
    if (k == 0) total = (n == 0 ? 1 : 0);
    for (long long code = 0; code < total; ++code) {
        long long c = code;
        for (int i = 0; i < n; ++i) f[(size_t)i] = (int)(c % k), c /= k;
        std::vector<int> cnt((size_t)k, 0);
        for (int x : f) ++cnt[(size_t)x];
        bool ok = true;
        for (int x : cnt)
            if ((rule == 1 && x > 1) || (rule == 2 && x == 0)) ok = false;
        if (!ok) continue;
        std::vector<int> key;
        if (bd && xd) key = f;
        else if (!bd && xd) key = cnt;
        else if (bd && !xd) {  // 箱を出てきた順に番号を振り直す
            std::vector<int> lab((size_t)k, -1);
            int nx = 0;
            for (int x : f) {
                if (lab[(size_t)x] < 0) lab[(size_t)x] = nx++;
                key.push_back(lab[(size_t)x]);
            }
        } else {
            key = cnt;
            std::sort(key.begin(), key.end());
        }
        seen.insert(key);
    }
    return (long long)seen.size();
}
static std::vector<std::vector<int>> all_perms(int n) {
    std::vector<int> p((size_t)n);
    std::iota(p.begin(), p.end(), 0);
    std::vector<std::vector<int>> res;
    do res.push_back(p);
    while (std::next_permutation(p.begin(), p.end()));
    return res;
}

template <class mint>
static void test_stirling(std::mt19937_64& rng, int& cases) {
    const int N = 70;
    const auto t2 = stirling2_table<mint>(N);
    const auto t1 = stirling1_table<mint>(N);
    for (int n = 0; n <= N; ++n) {
        const auto r2 = stirling2_row<mint>(n);
        const auto r1 = stirling1_row<mint>(n);
        COMB_CHECK(r2.size() == (size_t)n + 1 && r1.size() == (size_t)n + 1);
        for (int k = 0; k <= n; ++k, ++cases) {
            COMB_CHECK(r2[(size_t)k] == t2[(size_t)n][(size_t)k] && stirling2<mint>(n, k) == t2[(size_t)n][(size_t)k]);
            COMB_CHECK(r1[(size_t)k] == t1[(size_t)n][(size_t)k]);
        }
        COMB_CHECK(stirling2<mint>(n, n + 1) == mint(0) && stirling2<mint>(n, -1) == mint(0));
        mint b = 0;
        for (int k = 0; k <= n; ++k) b += t2[(size_t)n][(size_t)k];
        COMB_CHECK(bell<mint>(n) == b && bell_table<mint>(N)[(size_t)n] == b);
        // Σ_k c(n, k) = n!、x^n = Σ_k S(n, k) x(x−1)…(x−k+1)
        mint s1 = 0, f = 1;
        for (int k = 0; k <= n; ++k) s1 += t1[(size_t)n][(size_t)k];
        for (int i = 1; i <= n; ++i) f *= mint(i);
        COMB_CHECK(s1 == f);
        const long long x = rnd(rng, 0, 1000);
        mint lhs = mint(x).pow(n), rhs = 0, ff = 1;
        for (int k = 0; k <= n; ++k) rhs += t2[(size_t)n][(size_t)k] * ff, ff *= mint(x - k);
        COMB_CHECK(lhs == rhs);
    }
    // 大きい n：1 つだけ版と行の版
    for (int it = 0; it < 4 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1000, 3000);
        const auto r2 = stirling2_row<mint>(n);
        for (int t = 0; t < 20; ++t) {
            const int k = (int)rnd(rng, 0, n);
            COMB_CHECK(stirling2<mint>(n, k) == r2[(size_t)k]);
        }
        const auto r1 = stirling1_row<mint>(n);  // c(n, k) = (n−1) c(n−1, k) + c(n−1, k−1) を 1 段だけ確かめる
        const auto q1 = stirling1_row<mint>(n - 1);
        for (int k = 1; k < n; ++k) COMB_CHECK(r1[(size_t)k] == mint(n - 1) * q1[(size_t)k] + q1[(size_t)k - 1]);
    }
}

static void test_small_exact(int& cases) {
    // 順列を全部並べて：サイクルの数・完全順列
    const auto t1 = stirling1_table<long long>(8);
    const auto der = derangement_table<long long>(8);
    for (int n = 0; n <= 8; ++n) {
        std::vector<long long> cyc((size_t)n + 1, 0);
        long long d = 0;
        for (const auto& p : all_perms(n)) {
            std::vector<bool> vis((size_t)n, false);
            int c = 0;
            bool fixed = false;
            for (int i = 0; i < n; ++i) {
                if (p[(size_t)i] == i) fixed = true;
                if (vis[(size_t)i]) continue;
                ++c;
                for (int j = i; !vis[(size_t)j]; j = p[(size_t)j]) vis[(size_t)j] = true;
            }
            ++cyc[(size_t)c];
            if (!fixed) ++d;
        }
        for (int k = 0; k <= n; ++k, ++cases) COMB_CHECK(t1[(size_t)n][(size_t)k] == cyc[(size_t)k]);
        COMB_CHECK(der[(size_t)n] == d);
    }
    // 集合の分割を全部並べて（制限成長列）：S(n, k)・ベル数
    const auto t2 = stirling2_table<long long>(10);
    const auto bt = bell_table<long long>(10);
    for (int n = 0; n <= 10; ++n) {
        std::vector<long long> cnt((size_t)n + 1, 0);
        std::vector<int> a((size_t)n, 0);
        std::function<void(int, int)> rec = [&](int i, int mx) {
            if (i == n) { ++cnt[(size_t)mx]; return; }
            for (int v = 0; v <= mx; ++v) a[(size_t)i] = v, rec(i + 1, std::max(mx, v + 1));
        };
        rec(0, 0);
        long long b = 0;
        for (int k = 0; k <= n; ++k, ++cases) {
            COMB_CHECK(t2[(size_t)n][(size_t)k] == cnt[(size_t)k]);
            b += cnt[(size_t)k];
        }
        COMB_CHECK(bt[(size_t)n] == b);
    }
    COMB_CHECK(bell_table<long long>(25)[25] == 4638590332229999353LL);
    // 分割数：ナップサックの DP と
    const int N = 400;
    const auto p = partition_table<long long>(N);
    const auto pc = partition_count_table<long long>(N, 30);
    for (int k = 0; k <= 30; ++k) {  // k 以下の数を何個でも使うナップサック
        std::vector<long long> q((size_t)N + 1, 0);
        q[0] = 1;
        for (int j = 1; j <= k; ++j)
            for (int v = j; v <= N; ++v) q[(size_t)v] += q[(size_t)(v - j)];
        for (int v = 0; v <= N; ++v, ++cases) COMB_CHECK(pc[(size_t)k][(size_t)v] == q[(size_t)v]);
    }
    {
        std::vector<long long> q((size_t)N + 1, 0);
        q[0] = 1;
        for (int j = 1; j <= N; ++j)
            for (int v = j; v <= N; ++v) q[(size_t)v] += q[(size_t)(v - j)];
        for (int v = 0; v <= N; ++v, ++cases) COMB_CHECK(p[(size_t)v] == q[(size_t)v]);
    }
    COMB_CHECK(partition_table<long long>(100)[100] == 190569292LL);
    COMB_CHECK(derangement_table<long long>(0).size() == 1 && derangement_table<long long>(0)[0] == 1);
}

template <class mint>
static void test_twelve(std::mt19937_64& rng, int& cases) {
    for (int n = 0; n <= 6; ++n)
        for (int k = 0; k <= 6; ++k)
            for (int m = 0; m < 4; ++m)
                for (int rule = 0; rule < 3; ++rule, ++cases) {
                    const bool bd = m & 1, xd = m & 2;
                    COMB_CHECK(balls_into_boxes<mint>(n, k, bd, xd, rule) == mint(brute_twelve(n, k, bd, xd, rule)));
                }
    // 大きい n：Σ_{i ≤ k} S(n, i) と、全射の数（包除原理）と
    for (int it = 0; it < 200 * TEST_SCALE; ++it, ++cases) {
        const long long n = it % 2 ? rnd(rng, 0, 1000000000000000000LL) : rnd(rng, 0, 60), k = rnd(rng, 0, 40);
        mint s = 0;
        for (long long i = 0; i <= k; ++i) s += stirling2<mint>(n, i);
        COMB_CHECK(balls_into_boxes<mint>(n, k, true, false, 0) == s);
        const mint sur = inclusion_exclusion_count<mint>((int)k, [&](int m) { return mint(m).pow(n); });
        COMB_CHECK(surjections<mint>(n, k) == sur && balls_into_boxes<mint>(n, k, true, true, 2) == sur);
    }
    // 玉も箱も区別しない：分割数の表と
    const auto pc = partition_count_table<mint>(300, 300);
    for (int it = 0; it < 300; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 300), k = (int)rnd(rng, 0, 300);
        COMB_CHECK(balls_into_boxes<mint>(n, k, false, false, 0) == pc[(size_t)k][(size_t)n]);
        COMB_CHECK(balls_into_boxes<mint>(n, k, false, false, 2) == (n >= k ? pc[(size_t)k][(size_t)(n - k)] : mint(0)));
    }
}

static void test_inclusion(std::mt19937_64& rng, int& cases) {
    // [1, X] で、どの d_i でも割り切れない数（inclusion_exclusion）
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 8);
        const long long X = rnd(rng, 0, 3000);
        std::vector<long long> d((size_t)n);
        for (auto& x : d) x = rnd(rng, 1, 40);
        const long long got = inclusion_exclusion<long long>(n, [&](long long S) {
            long long l = 1;
            for (int i = 0; i < n; ++i)
                if (S >> i & 1) l = lcm_cap(l, d[(size_t)i], 1000000);
            return X / l;
        });
        long long e = 0;
        for (long long v = 1; v <= X; ++v) {
            bool any = false;
            for (long long x : d) any = any || v % x == 0;
            e += !any;
        }
        COMB_CHECK(got == e);
        // count_divisible_by_any（負の範囲も）
        const long long l = rnd(rng, -2000, 2000), r = rnd(rng, -2000, 2000);
        long long e2 = 0;
        for (long long v = l; v < r; ++v) {
            bool any = false;
            for (long long x : d) any = any || v % x == 0;
            e2 += any;
        }
        COMB_CHECK(count_divisible_by_any(l, r, d) == e2);
        COMB_CHECK(count_divisible_by_any_closed(l, r, d) == e2 + (r >= l && [&] { for (long long x : d) if (r % x == 0) return true; return false; }() ? 1 : 0));
    }
    // 大きい範囲・大きい d：Int128 で 2 つ・3 つの包除を直接
    using math_internal::ceil_div128;
    auto cm = [&](long long l, long long r, Int128 d) { return ceil_div128(r, d) - ceil_div128(l, d); };
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        long long l = rnd_ll(rng), r = rnd_ll(rng);
        if (l > r) std::swap(l, r);
        const long long a = it % 3 ? rnd(rng, 1, 1000000000000LL) : std::max(1LL, std::llabs(rnd_ll(rng) / 2)), b = rnd(rng, 2, 1000000), c = rnd(rng, 2, 1000);
        auto lc = [](Int128 x, Int128 y) {
            Int128 g = x, h = y;
            while (h) { const Int128 t = g % h; g = h, h = t; }
            return x / g * y;
        };
        const Int128 e = cm(l, r, a) + cm(l, r, b) + cm(l, r, c) - cm(l, r, lc(a, b)) - cm(l, r, lc(a, c)) - cm(l, r, lc(b, c)) + cm(l, r, lc(lc(a, b), c));
        if (e > (Int128)LLONG_MAX) continue;
        COMB_CHECK((Int128)count_divisible_by_any(l, r, {a, b, c}) == e);
    }
    COMB_CHECK((Int128)count_divisible_by_any(LLONG_MIN, LLONG_MAX, {4, 6}) == cm(LLONG_MIN, LLONG_MAX, 4) + cm(LLONG_MIN, LLONG_MAX, 6) - cm(LLONG_MIN, LLONG_MAX, 12));
    COMB_CHECK(count_divisible_by_any_closed(LLONG_MAX, LLONG_MAX, {7, 73}) == 1 && count_divisible_by_any_closed(LLONG_MAX, LLONG_MAX, {2}) == 0);
    COMB_CHECK(count_divisible_by_any(0, 1, {1000000007, 998244353}) == 1 && count_divisible_by_any(1, 1, {1}) == 0 && count_divisible_by_any(5, 9, {}) == 0);
    // 素数を 40 個：lcm がすぐ範囲を超えるので速い
    {
        std::vector<long long> ps;
        for (long long v = 2; ps.size() < 40; ++v) {
            bool pr = true;
            for (long long q = 2; q * q <= v; ++q) pr = pr && v % q != 0;
            if (pr) ps.push_back(v);
        }
        const long long X = 100000;
        long long e = 0;
        for (long long v = 1; v <= X; ++v) {
            bool any = false;
            for (long long p : ps) any = any || v % p == 0;
            e += any;
        }
        COMB_CHECK(count_divisible_by_any_closed(1, X, ps) == e);
        COMB_CHECK(count_divisible_by_any_closed(1, X, {2, 4, 6, 3, 9, 2}) == X / 2 + X / 3 - X / 6);
        cases += 2;
    }
}

template <class mint>
static void test_transforms(std::mt19937_64& rng, int& cases) {
    const Binomial<mint> B(500);
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, it % 4 ? 30 : 300);
        std::vector<mint> g((size_t)n);
        for (auto& x : g) x = mint(rnd(rng, 0, 1000000000));
        std::vector<mint> f((size_t)n, 0), gi((size_t)n, 0), A((size_t)n, 0), E((size_t)n, 0);
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b) {
                if (b <= a) f[(size_t)a] += B.C(a, b) * g[(size_t)b], gi[(size_t)a] += ((a - b) % 2 ? -B.C(a, b) : B.C(a, b)) * g[(size_t)b];
                if (b >= a) A[(size_t)a] += B.C(b, a) * g[(size_t)b], E[(size_t)a] += ((b - a) % 2 ? -B.C(b, a) : B.C(b, a)) * g[(size_t)b];
            }
        COMB_CHECK(binomial_transform(g) == f && binomial_inverse(g) == gi && at_least_from_exactly(g) == A && exactly_from_at_least(g) == E);
        COMB_CHECK(binomial_inverse(binomial_transform(g)) == g && exactly_from_at_least(at_least_from_exactly(g)) == g);
    }
    // q 二項係数：q-パスカルの漸化式 [n, k] = [n−1, k−1] + q^k [n−1, k] と比べる（位数の小さい q も）
    std::vector<long long> qs = {0, 1, 2, 3, -1, (long long)mint::mod() - 1, (long long)mint::mod() + 1, 1000000};
    if (mint::mod() == 998244353)
        for (long long d : {2, 4, 7, 8, 14, 17}) qs.push_back((long long)mint(3).pow((998244353 - 1) / d).val());  // 位数 d
    if (mint::mod() == 1000000007)
        for (long long d : {2, 3, 6, 41}) qs.push_back((long long)mint(5).pow((1000000007 - 1) / d).val());
    for (long long q : qs) {
        const int N = 60;
        std::vector<std::vector<mint>> t(N + 1, std::vector<mint>(N + 1, 0));
        const mint Q(q);
        for (int n = 0; n <= N; ++n) {
            t[(size_t)n][0] = 1;
            for (int k = 1; k <= n; ++k) t[(size_t)n][(size_t)k] = t[(size_t)n - 1][(size_t)k - 1] + Q.pow(k) * t[(size_t)n - 1][(size_t)k];
        }
        for (int n = 0; n <= N; ++n)
            for (int k = -1; k <= n + 1; ++k, ++cases) COMB_CHECK(q_binomial<mint>(n, k, q) == (k < 0 || k > n ? mint(0) : t[(size_t)n][(size_t)k]));
    }
}
static void test_q_binomial_small_mod(int& cases) {  // mod が小さい素数（q の位数 d で ⌊n/d⌋ ≥ mod になる）でも
    for (int p : {2, 3, 5, 7, 13}) {
        DynModInt::set_mod(p);
        const int N = 45;
        for (long long q = -2; q <= p + 1; ++q) {
            const DynModInt Q(q);
            std::vector<std::vector<DynModInt>> t(N + 1, std::vector<DynModInt>(N + 1, DynModInt(0)));
            for (int n = 0; n <= N; ++n) {
                t[(size_t)n][0] = 1;
                for (int k = 1; k <= n; ++k) t[(size_t)n][(size_t)k] = t[(size_t)n - 1][(size_t)k - 1] + Q.pow(k) * t[(size_t)n - 1][(size_t)k];
            }
            for (int n = 0; n <= N; ++n)
                for (int k = 0; k <= n; ++k, ++cases) COMB_CHECK(q_binomial<DynModInt>(n, k, q) == t[(size_t)n][(size_t)k]);
        }
    }
    DynModInt::set_mod(998244353);
}

static long long brute_syt(std::vector<int> sh, std::map<std::vector<int>, long long>& memo) {  // 角を 1 つずつ取る
    while (!sh.empty() && sh.back() == 0) sh.pop_back();
    if (sh.empty()) return 1;
    auto it = memo.find(sh);
    if (it != memo.end()) return it->second;
    long long r = 0;
    for (size_t i = 0; i < sh.size(); ++i)
        if (i + 1 == sh.size() || sh[i + 1] < sh[i]) {
            --sh[i];
            r += brute_syt(sh, memo);
            ++sh[i];
        }
    return memo[sh] = r;
}

static void test_enumerate(std::mt19937_64& rng, int& cases) {
    const std::vector<long long> bell_ll = bell_table<long long>(12), part_ll = partition_table<long long>(40);
    for (int n = 0; n <= 10; ++n) {
        long long cnt = 0;
        std::set<std::vector<int>> seen;
        for_each_set_partition(n, [&](const std::vector<int>& g, int k) {
            ++cnt;
            int mx = -1;
            for (int x : g) {  // 初めて出る順の番号
                COMB_CHECK(x <= mx + 1);
                mx = std::max(mx, x);
            }
            COMB_CHECK(mx + 1 == k && (int)g.size() == n);
            if (n <= 7) seen.insert(g);
        });
        COMB_CHECK(cnt == bell_ll[(size_t)n] && (n > 7 || (long long)seen.size() == cnt));
        ++cases;
    }
    const auto pct = partition_count_table<long long>(40, 40);  // pct[k][n]：n を k 個以下に分ける数 = 最大の部分が k 以下
    for (int n = 0; n <= 40; ++n) {
        long long cnt = 0;
        for_each_integer_partition(n, [&](const std::vector<int>& p) {
            ++cnt;
            int s = 0;
            for (size_t i = 0; i < p.size(); ++i) {
                s += p[i];
                COMB_CHECK(p[i] >= 1 && (i == 0 || p[i] <= p[i - 1]));
            }
            COMB_CHECK(s == n);
        });
        COMB_CHECK(cnt == part_ll[(size_t)n]);
        const int mp = (int)rnd(rng, 0, n + 1);
        long long c2 = 0;
        std::vector<int> prev;
        for_each_integer_partition(
            n, [&](const std::vector<int>& p) {
                ++c2;
                COMB_CHECK(p.empty() || p[0] <= mp);
                COMB_CHECK(prev.empty() || p < prev);  // 辞書順で大きい方から
                prev = p;
            },
            mp);
        COMB_CHECK(c2 == (mp == 0 ? (n == 0) : pct[(size_t)std::min(mp, n)][(size_t)n]));
        cases += 2;
    }
    // フック長公式
    std::map<std::vector<int>, long long> memo;
    for (int it = 0; it < 500; ++it, ++cases) {
        std::vector<int> sh;
        int rem = (int)rnd(rng, 0, 14), mx = rem;
        while (rem > 0) {
            const int x = (int)rnd(rng, 1, std::min(rem, mx));
            sh.push_back(x), rem -= x, mx = x;
        }
        COMB_CHECK(young_tableaux_count<ModInt998>(sh) == ModInt998(brute_syt(sh, memo)));
    }
    COMB_CHECK(young_tableaux_count<ModInt998>({3, 2}) == ModInt998(5) && young_tableaux_count<ModInt998>(std::vector<int>{}) == ModInt998(1));
    COMB_CHECK(young_tableaux_count<ModInt998>({5, 5}) == ModInt998(42));  // 2 行 n 列はカタラン数
    {  // 部分の数が多くても再帰しない
        long long calls = 0;
        size_t len = 0;
        for_each_integer_partition(1000000, [&](const std::vector<int>& p) { ++calls, len = p.size(); }, 1);
        COMB_CHECK(calls == 1 && len == 1000000);
        long long c0 = 0;
        for_each_integer_partition(5, [&](const std::vector<int>&) { ++c0; }, 0);
        for_each_integer_partition(0, [&](const std::vector<int>& p) { c0 += 10 + (long long)p.size(); });
        COMB_CHECK(c0 == 10);
        cases += 2;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_transforms<ModInt998>(rng, cases);
    test_transforms<ModInt107>(rng, cases);
    test_enumerate(rng, cases);
    test_q_binomial_small_mod(cases);
    test_stirling<ModInt998>(rng, cases);
    test_stirling<ModInt107>(rng, cases);  // NTT が使えない素数
    test_small_exact(cases);
    test_twelve<ModInt998>(rng, cases);
    test_twelve<ModInt107>(rng, cases);
    test_inclusion(rng, cases);
    {  // partition_row（級数の逆数）と五角数定理の表
        const auto a = partition_row<ModInt998>(5000), b = partition_table<ModInt998>(5000);
        const auto c = partition_row<ModInt107>(3000), d = partition_table<ModInt107>(3000);
        COMB_CHECK(a == b && c == d && partition_row<ModInt998>(0).size() == 1);
        cases += 3;
    }
    {  // 小さい素数の mod（DynModInt）：mod 7 でも n < 7 なら使える
        DynModInt::set_mod(7);
        COMB_CHECK(stirling2<DynModInt>(3, 2) == DynModInt(3) && surjections<DynModInt>(3, 2) == DynModInt(6) && bell<DynModInt>(4) == DynModInt(15 % 7));
        COMB_CHECK(balls_into_boxes<DynModInt>(3, 2, true, false, 0) == DynModInt(4) && surjections<DynModInt>(3, 2000000000) == DynModInt(0));
        DynModInt::set_mod(998244353);
        // k が巨大：階乗の表を作らずに
        using mint = ModInt998;
        COMB_CHECK(balls_into_boxes<mint>(2, 1000000000, true, true, 1) == mint(1000000000) * mint(999999999));
        COMB_CHECK(balls_into_boxes<mint>(3, 1000000000, false, true, 1) == mint(1000000000) * mint(999999999) * mint(999999998) / mint(6));
        COMB_CHECK(balls_into_boxes<mint>(2, 1000000000, false, true, 0) == mint(1000000001) * mint(1000000000) / mint(2));
        COMB_CHECK(balls_into_boxes<mint>(100000, 100000, false, false, 0) == partition_row<mint>(100000)[100000]);
        cases += 4;
    }
    std::cout << "[combinatorics] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[combinatorics] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    using mint = ModInt998;
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    const auto r2 = stirling2_row<mint>(500000);
    const long long t_s2 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto r1 = stirling1_row<mint>(200000);
    const long long t_s1 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto p = partition_row<mint>(1000000);
    const long long t_p = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto t2 = stirling2_table<mint>(3000);
    const long long t_tab = ms(t0);
    t0 = std::chrono::steady_clock::now();
    std::vector<long long> ds;
    for (long long v = 2; ds.size() < 20; ++v) {
        bool pr = true;
        for (long long q = 2; q * q <= v; ++q) pr = pr && v % q != 0;
        if (pr) ds.push_back(v);
    }
    const long long c = count_divisible_by_any(0, 1000000000000000000LL, ds);
    const long long t_ie = ms(t0);
    std::cout << "[combinatorics] benchmark: stirling2_row(5×10^5) " << t_s2 << " ms, stirling1_row(2×10^5) " << t_s1 << " ms, partition_row(10^6) " << t_p
              << " ms, stirling2_table(3000) " << t_tab << " ms, count_divisible_by_any（10^18, 素数 20 個）" << t_ie << " ms  (check "
              << (r2[1234] + r1[777] + p.back() + t2[3000][1500] + mint(c)).val() % 7 << ")\n";
}

}  // namespace comb_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace comb_demo {

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
#define SHOW(expr, comment) comb_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    comb_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    using mint = ModInt998;
    banner("数え上げ（combinatorics） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  stirling2 / stirling2_row / stirling2_table / stirling1_row / stirling1_table / bell / bell_table\n"
              << "  partition_table / partition_count_table / derangement_table / surjections / balls_into_boxes\n"
              << "  inclusion_exclusion / inclusion_exclusion_count / count_divisible_by_any(_closed)\n"
              << "  （この中では using mint = ModInt998; としている）\n";

    api("【1】 stirling2<mint>(n, k) / stirling2_row<mint>(n) / stirling2_table<T>(N)", "mint / vector<mint> / vector<vector<T>>",
        "第 2 種スターリング数 S(n, k)：n 個の区別できるものを k 個の空でない（区別しない）組に分ける数", "O(k log n) / O(n log n) / O(N^2)");
    SHOW(stirling2<mint>(5, 2), "{1..5} を 2 組に：2^4 − 1 = 15");
    SHOW(stirling2_row<mint>(5), "S(5, 0..5)");
    SHOW(stirling2<mint>(1000000000000000000LL, 3), "n が 10^18 でも k が小さければ");
    SHOW(stirling2_table<long long>(4)[4], "long long の表も（mod なし）");

    api("【2】 stirling1_row<mint>(n) / stirling1_table<T>(N)", "vector<mint> / vector<vector<T>>",
        "第 1 種スターリング数（符号なし）c(n, k)：n 個の順列でサイクルがちょうど k 個のものの数。x(x+1)…(x+n−1) の係数", "O(n log^2 n) / O(N^2)");
    SHOW(stirling1_row<mint>(4), "c(4, 0..4)。合計 24 = 4!");
    note("符号つき s(n, k) は (−1)^{n−k} c(n, k)");

    api("【3】 bell<mint>(n) / bell_table<T>(N)", "mint / vector<T>", "ベル数 B_n：n 個を空でない組にいくつかに分ける数（= Σ_k S(n, k)）", "O(n log n) / O(N^2)");
    SHOW(bell_table<long long>(8), "B_0..B_8");
    SHOW(bell<mint>(100000), "B_100000 mod 998244353");

    api("【4】 partition_row<mint>(N) / partition_table<T>(N) / partition_count_table<T>(N, K)", "vector<mint> / vector<T> / vector<vector<T>>",
        "分割数 p(n)：n を正の整数の和に分ける数（順番は区別しない）/ 同じく long long でも使える版 / t[k][n]：k 個以下に分ける数", "O(N log N) / O(N √N) / O(N K)");
    SHOW(partition_row<mint>(100000)[100000], "p(10^5) mod 998244353");
    SHOW(partition_table<long long>(10), "p(0..10)：5 = 5 = 4+1 = 3+2 = 3+1+1 = 2+2+1 = 2+1+1+1 = 1+1+1+1+1 で 7");
    SHOW(partition_count_table<long long>(10, 3)[3][10], "10 を 3 個以下に分ける");

    api("【5】 derangement_table<T>(N) / surjections<mint>(n, k)", "vector<T> / mint",
        "完全順列（どの i も i 番目に来ない順列）の数 / n 個から k 個への全射（どの行き先にも 1 個以上）の数", "O(N) / O(k log n)");
    SHOW(derangement_table<long long>(6), "D_0..D_6（プレゼント交換で全員が他人のものをもらう）");
    SHOW(surjections<mint>(5, 3), "3! S(5, 3) = 6 × 25");

    api("【6】 balls_into_boxes<mint>(n, k, 玉を区別するか, 箱を区別するか, rule)", "mint",
        "写像 12 相：n 個の玉を k 個の箱に入れる方法の数。rule：0 = 何個でもよい、1 = 各箱 1 個以下、2 = 各箱 1 個以上", "O(k log n) など（玉も箱も区別しないときは分割数）");
    {
        const char* names[4] = {"玉 区別しない・箱 区別しない", "玉 区別する  ・箱 区別しない", "玉 区別しない・箱 区別する  ", "玉 区別する  ・箱 区別する  "};
        note("n = 5、k = 3 の 12 通り：                  何個でも   1 個以下   1 個以上");
        for (int m = 0; m < 4; ++m) {
            std::cout << "    " << names[m];
            for (int rule = 0; rule < 3; ++rule) std::cout << std::setw(11) << balls_into_boxes<mint>(5, 3, m & 1, m & 2, rule).val();
            std::cout << "\n";
        }
        note("（3^5 = 243、S(5,1)+S(5,2)+S(5,3) = 41、p(5 を 3 個以下) = 5、7H5 = C(7, 5) = 21、全射 150、C(4, 2) = 6 …）");
    }
    SHOW(balls_into_boxes<mint>(2, 3, true, true, 1), "「1 個以下」は n ≤ k のときだけ 0 でない：3P2 = 6");

    api("【7】 inclusion_exclusion<T>(n, f)", "T", "Σ_{S ⊆ {0..n−1}} (−1)^{|S|} f(S)：f(S) =「条件 S を全部破る」数 → どの条件も破らない数", "O(2^n × f)");
    {
        // ABC297 F（Minimum Bounding Box 2）：h × w の長方形がちょうど外接長方形になる選び方を「4 辺のどれかが空」の包除で
        auto solve = [](long long H, long long W, long long K) {
            Binomial<mint> B(H * W);
            mint s = 0;
            for (long long h = 1; h <= H; ++h)
                for (long long w = 1; w <= W; ++w) {
                    const mint ways = inclusion_exclusion<mint>(4, [&](long long S) {  // ビット 0,1 = 上・下の行、2,3 = 左・右の列が空
                        const long long hh = std::max(0LL, h - (S & 1) - (S >> 1 & 1)), ww = std::max(0LL, w - (S >> 2 & 1) - (S >> 3 & 1));
                        return B.C(hh * ww, K);
                    });
                    s += ways * mint(h * w) * mint((H - h + 1) * (W - w + 1));
                }
            return s / B.C(H * W, K);
        };
        note("ABC297 F（Minimum Bounding Box 2）: 外接長方形がちょうど h × w になる選び方を「上下左右の 4 辺のどれかが空」の包除で");
        std::cout << "    (2, 2, 2) → " << solve(2, 2, 2) << "、(10, 10, 1) → " << solve(10, 10, 1) << "、(314, 159, 2653) → " << solve(314, 159, 2653) << "\n";
    }

    api("【8】 inclusion_exclusion_count<mint>(n, g)", "mint", "Σ_i (−1)^i C(n, i) g(n − i)：g(m) =「n 個のうち決まった m 個だけを使ってよい」数 → n 個を全部使う数", "O(n × g)");
    SHOW(inclusion_exclusion_count<mint>(3, [](int m) { return mint(m).pow(5); }), "5 個を 3 個へ、3 個とも使う（全射 150）");
    {
        // ABC242 F（Black and White Rooks）：r 行 c 列を全部使って k 個を置く数 f(r, c, k) を行・列の 2 重の包除で
        auto solve = [](int N, int M, int Bk, int Wk) {
            Binomial<mint> B(N * M + N + M);
            auto f = [&](int k) {
                std::vector<std::vector<mint>> t((size_t)N + 1, std::vector<mint>((size_t)M + 1));
                for (int r = 0; r <= N; ++r)
                    for (int c = 0; c <= M; ++c)
                        t[(size_t)r][(size_t)c] = inclusion_exclusion_count<mint>(r, [&](int rr) { return inclusion_exclusion_count<mint>(c, [&](int cc) { return B.C(rr * cc, k); }); });
                return t;
            };
            const auto fb = f(Bk), fw = f(Wk);
            mint s = 0;
            for (int r1 = 0; r1 <= N; ++r1)
                for (int c1 = 0; c1 <= M; ++c1) {
                    if (fb[(size_t)r1][(size_t)c1] == mint(0)) continue;
                    for (int r2 = 0; r1 + r2 <= N; ++r2)
                        for (int c2 = 0; c1 + c2 <= M; ++c2)
                            s += B.C(N, r1) * B.C(N - r1, r2) * B.C(M, c1) * B.C(M - c1, c2) * fb[(size_t)r1][(size_t)c1] * fw[(size_t)r2][(size_t)c2];
                }
            return s;
        };
        note("ABC242 F（Black and White Rooks）: 「r 行 c 列を全部使う」置き方を inclusion_exclusion_count の 2 重で");
        std::cout << "    (2, 2, 1, 1) → " << solve(2, 2, 1, 1) << "、(1, 2, 1, 1) → " << solve(1, 2, 1, 1) << "、(40, 40, 30, 30) → " << solve(40, 40, 30, 30) << "\n";
    }

    api("【9】 count_divisible_by_any(l, r, ds) / _closed(L, R, ds)", "long long", "[l, r)（[L, R]）の整数で ds のどれかで割り切れるものの数（包除原理。lcm が大きすぎる組は打ち切る）", "O(2^|ds|) 以下");
    SHOW(count_divisible_by_any_closed(1, 100, {2, 3, 5}), "1〜100 で 2・3・5 のどれかの倍数");
    SHOW(count_divisible_by_any(-10, 10, {4, 6}), "負の範囲も（0 も数える）");
    {
        // ABC131 C（Anti-Division）：A〜B で C でも D でも割り切れない数
        auto solve = [](long long A, long long B, long long C, long long D) { return (B - A + 1) - count_divisible_by_any_closed(A, B, {C, D}); };
        note("ABC131 C（Anti-Division）: (B − A + 1) − count_divisible_by_any_closed(A, B, {C, D})");
        std::cout << "    (4, 9, 2, 3) → " << solve(4, 9, 2, 3) << "、(10, 40, 6, 8) → " << solve(10, 40, 6, 8) << "、例 3 → "
                  << solve(314159265358979323LL, 846264338327950288LL, 419716939, 937510582) << "\n";
    }

    api("【10】 binomial_transform / binomial_inverse / at_least_from_exactly / exactly_from_at_least（<mint> の vector）", "vector<mint>",
        "f(n) = Σ_k C(n, k) g(k) とその逆 / A[k] = Σ_{j≥k} C(j, k) E[j] とその逆（「k 個を選んで満たす」の和 A から「ちょうど k 個」E へ）", "O(n log n)");
    note("「ちょうど k 個の条件を満たす数」は直接数えにくい → 「k 個の条件を選び、それを全部満たす」組の数 A[k] を数えて exactly_from_at_least");
    {
        // 例：1〜n の順列で、ちょうど k 個の i が p_i = i（不動点）になる数。A[k] = C(n, k) (n − k)!（k 個を選んで固定、残りは自由）
        const int n = 5;
        std::vector<mint> A(n + 1);
        const Binomial<mint> B(n);
        for (int k = 0; k <= n; ++k) A[(size_t)k] = B.C(n, k) * B.fact(n - k);
        std::cout << "    n = 5 の順列の不動点の数ごとの個数 → ";
        for (const mint& x : exactly_from_at_least(A)) std::cout << x << " ";
        std::cout << "（0 個は完全順列 44）\n";
    }
    SHOW(binomial_inverse<ModInt998>({1, 2, 4, 8, 16}), "2^n = Σ C(n, k)·1 なので 1 に戻る");

    api("【11】 for_each_set_partition(n, f) / for_each_integer_partition(n, f, max_part)", "なし（f を呼ぶ）",
        "{0..n−1} の組分けを全部 f(g, k)（g[i] = i の組の番号、k = 組の数）/ n の分割を全部 f(parts)（大きい順）", "O(B_n · n) / O(p(n) · n)");
    {
        long long c = 0;
        for_each_set_partition(4, [&](const std::vector<int>&, int k) { c += k == 2; });
        std::cout << "  for_each_set_partition(4, …)：4 個を 2 組に分ける方法 → " << c << "（S(4, 2) = 7）\n";
        // ABC390 D：袋の石をいくつかの組にまとめ、組ごとの和の XOR として出てくる値の種類（N ≤ 12）
        const std::vector<long long> a = {2, 5, 7};
        std::set<long long> vals;
        std::vector<long long> sum(a.size());
        for_each_set_partition((int)a.size(), [&](const std::vector<int>& g, int k) {
            std::fill(sum.begin(), sum.begin() + k, 0);
            for (size_t i = 0; i < a.size(); ++i) sum[(size_t)g[i]] += a[i];
            long long x = 0;
            for (int j = 0; j < k; ++j) x ^= sum[(size_t)j];
            vals.insert(x);
        });
        note("ABC390 D（Stone XOR）: N ≤ 12 の組分けを全部（B_12 = 420 万通り）たどって XOR の種類を数える");
        std::cout << "    {2, 5, 7} → " << vals.size() << "\n";
        std::cout << "  for_each_integer_partition(5, …) → ";
        for_each_integer_partition(5, [&](const std::vector<int>& p) {
            std::cout << "{";
            for (size_t i = 0; i < p.size(); ++i) std::cout << p[i] << (i + 1 < p.size() ? "," : "");
            std::cout << "} ";
        });
        std::cout << " // 大きい順、p(5) = 7 通り\n";
    }

    api("【12】 q_binomial<mint>(n, k, q) / young_tableaux_count<mint>(shape)", "mint",
        "q 二項係数 [n, k]_q（F_q^n の k 次元部分空間の数。q^d ≡ 1 でも q-Lucas で正しく）/ 標準ヤング盤の数（フック長公式）", "O(n)");
    SHOW(q_binomial<ModInt998>(4, 2, 2), "F_2^4 の 2 次元部分空間の数");
    SHOW(q_binomial<ModInt998>(10, 3, 1), "q = 1 なら普通の C(10, 3)");
    SHOW(young_tableaux_count<ModInt998>({3, 2, 1}), "形 (3, 2, 1) に 1〜6 を行・列とも増えるように入れる数");
    SHOW(young_tableaux_count<ModInt998>({10, 10}), "2 行 10 列 = カタラン数 C_10");
}

}  // namespace comb_demo

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
    comb_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    comb_test::benchmark();    // 速度計測
    comb_demo::demo();         // 使い方デモ
    return 0;
}
