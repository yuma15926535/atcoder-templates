// =============================================================================
//  桁・進数・回文数（桁和・n 進数・回文数（文字列の大きなものも）・数をつなげる・大きな数の文字列・繰り上がり・桁の数え上げ）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 01 桁・進数
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  機能一覧                                                                    計算量
//   桁（base は 2 以上。負の数は絶対値で）
//    digits_of(n, base = 10)            各桁 vector<int>（上の桁から。0 なら {0}。base ≤ 2^31 − 1）  O(桁数)
//    digit_sum(n, base) / num_digits(n, base)   桁和 / 桁数（0 は 1 桁）
//    from_digits(ds, base)              桁の並びから数に（溢れたら止める）
//   n 進数（2〜36 進と −2〜−36 進。10 以上の桁は A〜Z）
//    to_base(n, base)                   文字列に（0・負の数も。負の進数なら符号なしで表せる）
//    from_base(s, base)                 文字列から（小文字も可。溢れたら止める）
//    from_base_cap(s, base, cap)        base は 2 以上いくつでも：min(値, cap)（溢れない）
//    base_value_leq(s, base, limit)     s を base 進数で読んだ値 ≤ limit か（溢れない。ABC192 D の二分探索に）
//   回文数（10 進。0 以上。k 番目は 0 から数える：0 番目が 0）
//    is_palindrome_num(n, base = 10)    回文か
//    kth_palindrome(k)                  k 番目の回文数（10^18 を超えるので文字列）                O(桁数)
//    palindrome_rank(n)                 回文数 n が何番目か（kth_palindrome の逆）
//    count_palindromes_up_to(n) / count_palindromes_below(n)     [0, n] / [0, n) の回文数の個数
//    count_palindromes_of_length(len)   len 桁の回文数の個数
//    ceil_palindrome(n) / floor_palindrome(n)    n 以上で最小 / n 以下で最大の回文数
//    make_palindrome(half, odd)         左半分から作る（123 → odd: 12321 / even: 123321）
//   大きな回文数（10 進の文字列。long long を超える数に）
//    ceil_palindrome_str(s) / floor_palindrome_str(s)    s 以上で最小 / s 以下で最大の回文数（文字列）   O(|s|)
//    count_palindromes_up_to_str(s) / count_palindromes_below_str(s)   [0, s] / [0, s) の回文数の個数（Int128。75 桁まで）
//    make_palindrome_str(half, odd) / is_palindrome_str(s)   左半分から作る / 文字列が回文か
//   数をつなげる（0 以上）
//    concat_value(a, b)                 a の後ろに b（Int128）
//    repeat_concat_value(x, k)          x を k 回つなげた数（Int128。収まらなければ止める）
//    concat_mod(a, b, m) / repeat_concat_mod(x, k, m)    つなげた数 mod m / x を k 回つなげた数 mod m   O(log k)
//    repunit_mod(k, m) / pow10_mod(k, m)                 11…1（k 個）mod m / 10^k mod m
//   大きな数の文字列（10 進）
//    str_mod(s, m, base = 10) / str_num_compare(a, b)    値 mod m（base 進数 2〜36 でも）/ 大小（-1 / 0 / 1。符号・先頭の 0 も扱う）  O(桁数)
//    carry_to_string(c, base = 10)      繰り上がり：c[i] が base^i の係数（下の桁から。畳み込みの結果など）の数を文字列に  O(|c| + 桁数)
//    substr_value(s, l, r)              s[l, r) の値（溢れたら止める）
//    DigitPrefix P(s, m)                P.value(l, r)：s[l, r) の値 mod m を O(1)                 前処理 O(|s|)
//   桁の数え上げ
//    digit_length_blocks(l, r)          [l, r) を桁数ごとの区間 DigitLengthBlock{len, lo, hi} に分ける
//    total_digits_up_to(N)              1〜N の桁数の合計（Int128）
//    count_digit_up_to(N, d)            1〜N を書いたときの数字 d の個数（Int128。先頭の 0 は数えない）  O(桁数)
//    digit_sum_total_up_to(N)           1〜N の桁和の合計（Int128）
//    numbers_from_digits(ds, max_len)   数字 ds だけでできた 1〜max_len 桁の数（小さい順）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    dig_test::random_test()  … 文字列・128 ビット・全部数える愚直解との突き合わせ
//    dig_test::benchmark()    … 速度計測
//    dig_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_DIGITS
#define GG_DIGITS
// ---- 桁・進数・回文数（数学第二世代）  桁和・n 進数・回文数（文字列の大きなものも）・数をつなげる・大きな数の文字列・繰り上がり・桁の数え上げ ----
namespace math_internal {
inline int digit_value(char c) {  // '0'〜'9'・'A'〜'Z'・'a'〜'z' の値（それ以外は -1）
    if ('0' <= c && c <= '9') return c - '0';
    if ('A' <= c && c <= 'Z') return c - 'A' + 10;
    if ('a' <= c && c <= 'z') return c - 'a' + 10;
    return -1;
}
inline void check_base(long long base, const char* func) {
    if (base < 2) fail(func, "base ≥ 2 で使う");
}
inline long long checked_mod(long long m, const char* func) {  // m ≥ 1 を確かめて返す（% m の前に）
    if (m < 1) fail(func, "m ≥ 1 で使う");
    return m;
}
inline Int128 pow10_128(int k) {
    Int128 r = 1;
    while (k-- > 0) r *= 10;
    return r;
}
inline std::string normalize_nonneg_str(const std::string& s, const char* func) {  // 先頭の + と余分な 0 を取る（0 以上の 10 進でなければ止める）
    size_t i = (!s.empty() && s[0] == '+') ? 1 : 0;
    if (i == s.size()) fail(func, "数字がない");
    for (size_t j = i; j < s.size(); ++j)
        if (s[j] < '0' || s[j] > '9') fail(func, "0 以上の 10 進の数として読めない");
    while (i + 1 < s.size() && s[i] == '0') ++i;
    return s.substr(i);
}
inline std::string str_inc(std::string s) {  // 0 以上の 10 進の文字列 + 1
    int i = (int)s.size() - 1;
    while (i >= 0 && s[i] == '9') s[i--] = '0';
    if (i < 0) s.insert(s.begin(), '1');
    else ++s[i];
    return s;
}
inline std::string str_dec(std::string s) {  // 1 以上の 10 進の文字列 − 1（先頭の 0 を取る）
    int i = (int)s.size() - 1;
    while (s[i] == '0') s[i--] = '9';
    --s[i];
    if (s.size() > 1 && s[0] == '0') s.erase(s.begin());
    return s;
}
inline Int128 str_to_128(const std::string& s) {  // 38 桁までの 0 以上の 10 進
    Int128 v = 0;
    for (char c : s) v = v * 10 + (c - '0');
    return v;
}
}  // namespace math_internal

// ---- 桁 ----
inline std::vector<int> digits_of(long long n, long long base = 10) {  // |n| の各桁（上の桁から。0 なら {0}。base ≤ 2^31 − 1）
    math_internal::check_base(base, "digits_of");
    if (base > 2147483647LL) math_internal::fail("digits_of", "桁を int で返すので base ≤ 2^31 − 1 で使う");
    unsigned long long u = math_internal::uabs(n);
    std::vector<int> ds;
    do ds.push_back((int)(u % (unsigned long long)base)), u /= (unsigned long long)base;
    while (u);
    std::reverse(ds.begin(), ds.end());
    return ds;
}
inline long long digit_sum(long long n, long long base = 10) {  // |n| の桁和
    math_internal::check_base(base, "digit_sum");
    unsigned long long u = math_internal::uabs(n), s = 0;
    for (; u; u /= (unsigned long long)base) s += u % (unsigned long long)base;
    return (long long)s;
}
inline int num_digits(long long n, long long base = 10) {  // |n| の桁数（0 は 1 桁）
    math_internal::check_base(base, "num_digits");
    unsigned long long u = math_internal::uabs(n);
    int d = 1;
    while (u >= (unsigned long long)base) u /= (unsigned long long)base, ++d;
    return d;
}
inline long long from_digits(const std::vector<int>& ds, long long base = 10) {  // 上の桁からの並びを数に（溢れたら止める）
    math_internal::check_base(base, "from_digits");
    long long v = 0;
    for (int d : ds) {
        if (d < 0 || d >= base) math_internal::fail("from_digits", "桁が 0 以上 base 未満でない");
        if (__builtin_mul_overflow(v, base, &v) || __builtin_add_overflow(v, (long long)d, &v)) math_internal::fail("from_digits", "long long に収まらない");
    }
    return v;
}

// ---- n 進数（2〜36 進と −2〜−36 進。10 以上の桁は A〜Z） ----
inline std::string to_base(long long n, int base) {
    if (!((2 <= base && base <= 36) || (-36 <= base && base <= -2))) math_internal::fail("to_base", "base は 2〜36 か −36〜−2");
    static const char* D = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    if (n == 0) return "0";
    std::string s;
    if (base > 0) {
        unsigned long long u = math_internal::uabs(n);
        for (; u; u /= (unsigned long long)base) s += D[u % (unsigned long long)base];
        if (n < 0) s += '-';
    } else {  // 負の進数：余りを 0〜|base|-1 にそろえる
        Int128 x = n;
        while (x != 0) {
            Int128 r = x % base;
            x /= base;
            if (r < 0) r -= base, x += 1;
            s += D[(int)r];
        }
    }
    std::reverse(s.begin(), s.end());
    return s;
}
inline long long from_base(const std::string& s, int base) {  // to_base の逆（小文字も可。正の進数なら先頭の - も可。溢れたら止める）
    if (!((2 <= base && base <= 36) || (-36 <= base && base <= -2))) math_internal::fail("from_base", "base は 2〜36 か −36〜−2");
    size_t i = 0;
    const bool neg = base > 0 && !s.empty() && s[0] == '-';
    if (neg || (!s.empty() && s[0] == '+')) ++i;
    if (i == s.size()) math_internal::fail("from_base", "数字がない");
    Int128 v = 0;
    const Int128 lim = (Int128)9223372036854775807LL + 1;
    for (; i < s.size(); ++i) {
        const int d = math_internal::digit_value(s[i]);
        if (d < 0 || d >= (base < 0 ? -base : base)) math_internal::fail("from_base", "その進数の桁でない文字がある");
        v = v * base + d;
        if (v > lim || v < -lim) math_internal::fail("from_base", "long long に収まらない");
    }
    if (neg) v = -v;
    if (v > (Int128)9223372036854775807LL || v < -lim) math_internal::fail("from_base", "long long に収まらない");
    return (long long)v;
}
// s（0〜9・A〜Z）を base 進数（2 以上。36 を超えてもよい）として読んだ値と cap の小さい方（溢れない。0 以上の数として読む）
inline long long from_base_cap(const std::string& s, long long base, long long cap) {
    math_internal::check_base(base, "from_base_cap");
    if (cap < 0) math_internal::fail("from_base_cap", "cap ≥ 0 で使う");
    if (s.empty()) math_internal::fail("from_base_cap", "数字がない");
    Int128 v = 0;
    for (char c : s) {
        const int d = math_internal::digit_value(c);
        if (d < 0 || d >= base) math_internal::fail("from_base_cap", "その進数の桁でない文字がある");
        if (v < cap) v = std::min<Int128>(v * base + d, cap);  // cap に届いたらそのまま（桁の確認だけ続ける）
    }
    return (long long)v;
}
inline bool base_value_leq(const std::string& s, long long base, long long limit) {  // s を base 進数で読んだ値 ≤ limit か（溢れない。ABC192 D）
    math_internal::check_base(base, "base_value_leq");
    if (s.empty()) math_internal::fail("base_value_leq", "数字がない");
    Int128 v = 0;
    bool over = limit < 0;
    for (char c : s) {
        const int d = math_internal::digit_value(c);
        if (d < 0 || d >= base) math_internal::fail("base_value_leq", "その進数の桁でない文字がある");
        if (!over) v = v * base + d, over = v > limit;
    }
    return !over;
}

// ---- 回文数（0 以上の整数。k 番目は 0 から数える：0 番目が 0） ----
inline bool is_palindrome_num(long long n, long long base = 10) {  // n ≥ 0 が base 進数で回文か
    if (n < 0) return false;
    const std::vector<int> d = digits_of(n, base);
    for (size_t i = 0, j = d.size() - 1; i < j; ++i, --j)
        if (d[i] != d[j]) return false;
    return true;
}
inline long long count_palindromes_of_length(int len) {  // len 桁の回文数の個数（1 桁は 0〜9 の 10 個。len ≤ 38 まで）
    if (len < 1) return 0;
    if (len == 1) return 10;
    return 9 * pow10_ll((len + 1) / 2 - 1);
}
inline std::string kth_palindrome(long long k) {  // 小さい方から k 番目の回文数（0 番目は 0）。10^18 を超えるので文字列
    if (k < 0) math_internal::fail("kth_palindrome", "k ≥ 0 で使う");
    if (k < 10) return std::string(1, char('0' + k));
    k -= 10;
    for (int len = 2;; ++len) {
        const long long c = count_palindromes_of_length(len);
        if (k < c) {
            const std::string half = std::to_string(pow10_ll((len + 1) / 2 - 1) + k);
            std::string rev(half.rbegin() + (len % 2), half.rend());
            return half + rev;
        }
        k -= c;
    }
}
inline long long make_palindrome(long long half, bool odd) {  // half（≥ 0）を左半分にした回文数（odd なら真ん中の桁を共有：123 → 12321 / 123321）
    if (half < 0) math_internal::fail("make_palindrome", "half ≥ 0 で使う");
    Int128 v = half;
    for (long long t = odd ? half / 10 : half; t > 0; t /= 10) {
        v = v * 10 + t % 10;
        if (v > (Int128)9223372036854775807LL) math_internal::fail("make_palindrome", "long long に収まらない");
    }
    return (long long)v;
}
inline long long count_palindromes_up_to(long long n) {  // 0 以上 n 以下の回文数の個数（n < 0 なら 0）
    if (n < 0) return 0;
    const int len = num_digits(n);
    long long c = 0;
    for (int l = 1; l < len; ++l) c += count_palindromes_of_length(l);
    if (len == 1) return n + 1;
    const long long hdiv = pow10_ll(len / 2), half = n / hdiv, lo = pow10_ll((len + 1) / 2 - 1);
    c += half - lo;  // 左半分が half より小さいもの
    if (make_palindrome(half, len % 2) <= n) ++c;
    return c;
}
inline long long count_palindromes_below(long long n) { return n <= 0 ? 0 : count_palindromes_up_to(n - 1); }  // [0, n) の回文数の個数
inline long long palindrome_rank(long long n) {  // 回文数 n が小さい方から何番目か（0 から数える。kth_palindrome の逆）
    if (!is_palindrome_num(n)) math_internal::fail("palindrome_rank", "n が回文数でない");
    return count_palindromes_below(n);
}
inline long long ceil_palindrome(long long n) {  // n 以上で最小の回文数（long long に収まらなければ止める）
    if (n <= 0) return 0;
    const int len = num_digits(n);
    if (len == 1) return n;
    const long long half = n / pow10_ll(len / 2);
    const long long p = make_palindrome(half, len % 2);
    if (p >= n) return p;
    if (half + 1 == pow10_ll(num_digits(half))) {  // 99…9 → 次の桁数の最小 10…01
        if (len >= 19) math_internal::fail("ceil_palindrome", "long long に収まらない");
        return pow10_ll(len) + 1;
    }
    return make_palindrome(half + 1, len % 2);
}
inline long long floor_palindrome(long long n) {  // n 以下で最大の回文数（n ≥ 0）
    if (n < 0) math_internal::fail("floor_palindrome", "n ≥ 0 で使う");
    const int len = num_digits(n);
    if (len == 1) return n;
    const long long half = n / pow10_ll(len / 2);
    const long long p = make_palindrome(half, len % 2);
    if (p <= n) return p;
    if (half == pow10_ll(num_digits(half) - 1)) return pow10_ll(len - 1) - 1;  // 10…0 → 前の桁数の最大 99…9
    return make_palindrome(half - 1, len % 2);
}

// ---- 大きな回文数（10 進の文字列。0 以上。先頭の 0 は無視。long long を超える数に） ----
inline bool is_palindrome_str(const std::string& s) {  // 文字列 s が回文か（数でなくてもよい）
    for (size_t i = 0, j = s.size(); i + 1 < j; ++i, --j)
        if (s[i] != s[j - 1]) return false;
    return true;
}
inline std::string make_palindrome_str(const std::string& half, bool odd) {  // half を左半分にした回文（odd なら真ん中の文字を共有："123" → "12321" / "123321"）
    if (half.empty()) math_internal::fail("make_palindrome_str", "half が空");
    std::string r = half;
    r.append(half.rbegin() + (odd ? 1 : 0), half.rend());
    return r;
}
inline std::string ceil_palindrome_str(const std::string& n) {  // n 以上で最小の回文数
    const std::string s = math_internal::normalize_nonneg_str(n, "ceil_palindrome_str");
    const size_t len = s.size();
    if (len == 1) return s;
    const std::string half = s.substr(0, (len + 1) / 2);
    const std::string p = make_palindrome_str(half, len % 2);
    if (p >= s) return p;  // 同じ長さなので文字列の比較 = 数の比較
    if (half.find_first_not_of('9') == std::string::npos) return "1" + std::string(len - 1, '0') + "1";  // 99…9 → 10…01
    return make_palindrome_str(math_internal::str_inc(half), len % 2);
}
inline std::string floor_palindrome_str(const std::string& n) {  // n 以下で最大の回文数
    const std::string s = math_internal::normalize_nonneg_str(n, "floor_palindrome_str");
    const size_t len = s.size();
    if (len == 1) return s;
    const std::string half = s.substr(0, (len + 1) / 2);
    const std::string p = make_palindrome_str(half, len % 2);
    if (p <= s) return p;
    if (half == "1" + std::string(half.size() - 1, '0')) return std::string(len - 1, '9');  // 10…0 → 99…9
    return make_palindrome_str(math_internal::str_dec(half), len % 2);
}
inline Int128 count_palindromes_up_to_str(const std::string& n) {  // [0, n] の回文数の個数（n は 75 桁まで）
    const std::string s = math_internal::normalize_nonneg_str(n, "count_palindromes_up_to_str");
    const int len = (int)s.size();
    if (len > 75) math_internal::fail("count_palindromes_up_to_str", "Int128 に収まるよう 75 桁まで");
    if (len == 1) return s[0] - '0' + 1;
    Int128 c = 10;
    for (int l = 2; l < len; ++l) c += 9 * math_internal::pow10_128((l + 1) / 2 - 1);
    const std::string half = s.substr(0, (len + 1) / 2);
    c += math_internal::str_to_128(half) - math_internal::pow10_128((len + 1) / 2 - 1);  // 左半分が half より小さいもの
    if (make_palindrome_str(half, len % 2) <= s) ++c;
    return c;
}
inline Int128 count_palindromes_below_str(const std::string& n) {  // [0, n) の回文数の個数（n は 75 桁まで）
    const std::string s = math_internal::normalize_nonneg_str(n, "count_palindromes_below_str");
    return count_palindromes_up_to_str(s) - (is_palindrome_str(s) ? 1 : 0);
}

// ---- 数をつなげる（a, b ≥ 0） ----
inline Int128 concat_value(long long a, long long b) {  // a の後ろに b をつなげた数（"12" + "345" → 12345）
    if (a < 0 || b < 0) math_internal::fail("concat_value", "a, b ≥ 0 で使う");
    return (Int128)a * math_internal::pow10_128(num_digits(b)) + b;
}
inline Int128 repeat_concat_value(long long x, long long k) {  // x（≥ 0）を k 回つなげた数（Int128 に収まらなければ止める。k = 0 なら 0）
    if (x < 0 || k < 0) math_internal::fail("repeat_concat_value", "x, k ≥ 0 で使う");
    if (x == 0) return 0;
    const Int128 p = math_internal::pow10_128(num_digits(x)), LIM = (Int128)(~(UInt128)0 >> 1);
    Int128 v = 0;
    for (long long i = 0; i < k; ++i) {
        if (v > (LIM - x) / p) math_internal::fail("repeat_concat_value", "Int128 に収まらない");
        v = v * p + x;
    }
    return v;
}
inline long long pow10_mod(long long k, long long m) { return mod_pow(10, k, m); }  // 10^k mod m
inline long long concat_mod(long long a, long long b, long long m) {  // (a の後ろに b をつなげた数) mod m
    if (a < 0 || b < 0) math_internal::fail("concat_mod", "a, b ≥ 0 で使う");
    math_internal::checked_mod(m, "concat_mod");
    return (long long)(((Int128)mod_mul(a, pow10_mod(num_digits(b), m), m) + b % m) % m);
}
inline long long repeat_concat_mod(long long x, long long k, long long m) {  // x（≥ 0）を k 回つなげた数 mod m（ABC357 D）
    if (x < 0 || k < 0) math_internal::fail("repeat_concat_mod", "x, k ≥ 0 で使う");
    math_internal::checked_mod(m, "repeat_concat_mod");
    return mod_mul(x % m, geom_sum_mod(pow10_mod(num_digits(x), m), k, m), m);
}
inline long long repunit_mod(long long k, long long m) { return geom_sum_mod(10, k, m); }  // 1 を k 個並べた数 mod m

// ---- 大きな数の文字列（10 進。先頭の - と 0 があってもよい） ----
inline long long str_mod(const std::string& s, long long m, int base = 10) {  // 値 mod m（[0, m)）。base 進数（2〜36。10 以上の桁は A〜Z / a〜z）としても読める
    if (m < 1) math_internal::fail("str_mod", "m ≥ 1 で使う");
    if (base < 2 || base > 36) math_internal::fail("str_mod", "base は 2〜36");
    size_t i = (!s.empty() && (s[0] == '-' || s[0] == '+')) ? 1 : 0;
    if (i == s.size()) math_internal::fail("str_mod", "数字がない");
    long long r = 0;
    for (; i < s.size(); ++i) {
        const int d = math_internal::digit_value(s[i]);
        if (d < 0 || d >= base) math_internal::fail("str_mod", "数として読めない");
        r = (long long)(((Int128)r * base + d) % m);
    }
    return s[0] == '-' ? (m - r) % m : r;
}
// 繰り上がり：c[i] が base^i の係数（下の桁から。畳み込みの結果など。負もよいが全体は 0 以上）の数を base 進の文字列（上の桁から）に
template <class V>
std::string carry_to_string(const V& c, int base = 10) {
    if (base < 2 || base > 36) math_internal::fail("carry_to_string", "base は 2〜36");
    static const char* D = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    std::string s;
    Int128 carry = 0;
    auto push = [&](Int128 v) {  // v を足して 1 桁出す（余りは 0〜base−1、商は切り捨て）
        Int128 r = v % base, q = v / base;
        if (r < 0) r += base, --q;
        s += D[(int)r];
        carry = q;
    };
    for (const auto& x : c) push(carry + (Int128)x);
    while (carry > 0) push(carry);
    if (carry < 0) math_internal::fail("carry_to_string", "値が負");
    while (s.size() > 1 && s.back() == '0') s.pop_back();
    if (s.empty()) s = "0";
    std::reverse(s.begin(), s.end());
    return s;
}
inline std::string carry_to_string(std::initializer_list<long long> c, int base = 10) { return carry_to_string<std::initializer_list<long long>>(c, base); }
inline int str_num_compare(const std::string& a, const std::string& b) {  // 10 進の整数の文字列の大小（-1 / 0 / 1）
    auto parse = [](const std::string& s, bool& neg) {
        size_t i = (!s.empty() && (s[0] == '-' || s[0] == '+')) ? 1 : 0;
        neg = !s.empty() && s[0] == '-';
        while (i + 1 < s.size() && s[i] == '0') ++i;
        std::string d = s.substr(i);
        if (d.empty()) math_internal::fail("str_num_compare", "数字がない");
        for (char c : d)
            if (c < '0' || c > '9') math_internal::fail("str_num_compare", "数として読めない");
        if (d == "0") neg = false;
        return d;
    };
    bool na, nb;
    const std::string x = parse(a, na), y = parse(b, nb);
    if (na != nb) return na ? -1 : 1;
    int c = x.size() != y.size() ? (x.size() < y.size() ? -1 : 1) : (x < y ? -1 : (x > y ? 1 : 0));
    return na ? -c : c;
}
inline long long substr_value(const std::string& s, int l, int r) {  // 数字の文字列 s の [l, r) の値（l ≥ r なら 0。溢れたら止める）
    if (l < r && (l < 0 || r > (int)s.size())) math_internal::fail("substr_value", "範囲が文字列の外");
    long long v = 0;
    for (int i = l; i < r; ++i) {
        if (s[i] < '0' || s[i] > '9') math_internal::fail("substr_value", "数字でない文字がある");
        if (__builtin_mul_overflow(v, 10LL, &v) || __builtin_add_overflow(v, (long long)(s[i] - '0'), &v)) math_internal::fail("substr_value", "long long に収まらない");
    }
    return v;
}
struct DigitPrefix {  // 数字の文字列 s の部分文字列 [l, r) の値を mod m で O(1)
    long long m;
    std::vector<long long> pre, pw;  // pre[i] = s[0, i) の値 mod m、pw[i] = 10^i mod m
    DigitPrefix(const std::string& s, long long m_) : m(math_internal::checked_mod(m_, "DigitPrefix")), pre(s.size() + 1, 0), pw(s.size() + 1, 1 % m) {
        for (size_t i = 0; i < s.size(); ++i) {
            if (s[i] < '0' || s[i] > '9') math_internal::fail("DigitPrefix", "数字でない文字がある");
            pre[i + 1] = (long long)(((Int128)pre[i] * 10 + (s[i] - '0')) % m);
            pw[i + 1] = (long long)((Int128)pw[i] * 10 % m);
        }
    }
    long long value(int l, int r) const {  // s[l, r) の値 mod m（l ≥ r なら 0。0 ≤ l、r ≤ |s|）
        if (l >= r) return 0;
        if (l < 0 || r >= (int)pre.size()) math_internal::fail("DigitPrefix::value", "範囲が文字列の外");
        return (long long)((((Int128)pre[r] - (Int128)pre[l] * pw[r - l]) % m + m) % m);
    }
};

// ---- 桁の数え上げ（1〜N） ----
struct DigitLengthBlock {  // 桁数が len の数の範囲 [lo, hi)
    int len;
    long long lo, hi;
};
inline std::vector<DigitLengthBlock> digit_length_blocks(long long l, long long r) {  // [l, r)（0 ≤ l）を桁数ごとに分ける
    std::vector<DigitLengthBlock> res;
    if (l < 0) math_internal::fail("digit_length_blocks", "l ≥ 0 で使う");
    for (long long x = l; x < r;) {
        const int len = num_digits(x);
        const long long end = len >= 19 ? r : std::min(r, pow10_ll(len));
        res.push_back({len, x, end});
        x = end;
    }
    return res;
}
inline Int128 total_digits_up_to(long long N) {  // 1〜N の桁数の合計（"123…N" の長さ）
    Int128 s = 0;
    if (N < 1) return 0;
    for (const DigitLengthBlock& b : digit_length_blocks(1, N + (N < std::numeric_limits<long long>::max()))) s += (Int128)b.len * (b.hi - b.lo);
    if (N == std::numeric_limits<long long>::max()) s += 19;
    return s;
}
inline Int128 count_digit_up_to(long long N, int d) {  // 1〜N を書いたとき数字 d（0〜9）が出る回数（先頭の 0 は数えない）
    if (d < 0 || d > 9) math_internal::fail("count_digit_up_to", "d は 0〜9");
    if (N < 1) return 0;
    Int128 c = 0;
    for (Int128 p = 1; p <= N; p *= 10) {  // 各位について
        const Int128 high = N / (p * 10), cur = N / p % 10, low = N % p;
        if (d == 0) {
            if (high == 0) continue;  // 先頭の位
            c += (high - 1) * p + (cur > 0 ? p : low + 1);
        } else {
            c += high * p + (cur > d ? p : (cur == d ? low + 1 : 0));
        }
    }
    return c;
}
inline Int128 digit_sum_total_up_to(long long N) {  // 1〜N の桁和の合計
    Int128 s = 0;
    for (int d = 1; d <= 9; ++d) s += count_digit_up_to(N, d) * d;
    return s;
}

// ---- 使える数字が決まった数の列挙 ----
// 数字の集合 ds（0〜9）だけを使う 1〜max_len 桁の数（先頭は 0 以外。0 が ds にあれば 0 も）を小さい順に。max_len ≤ 18
inline std::vector<long long> numbers_from_digits(std::vector<int> ds, int max_len) {
    if (max_len < 1 || max_len > 18) math_internal::fail("numbers_from_digits", "1 ≤ max_len ≤ 18 で使う");
    for (int d : ds)
        if (d < 0 || d > 9) math_internal::fail("numbers_from_digits", "数字は 0〜9");
    std::sort(ds.begin(), ds.end());
    ds.erase(std::unique(ds.begin(), ds.end()), ds.end());
    std::vector<long long> res, cur;
    if (!ds.empty() && ds[0] == 0) res.push_back(0);
    for (int d : ds)
        if (d > 0) cur.push_back(d);
    for (int len = 1; len <= max_len && !cur.empty(); ++len) {
        res.insert(res.end(), cur.begin(), cur.end());
        if (len == max_len) break;
        if (cur.size() * ds.size() > 50000000) math_internal::fail("numbers_from_digits", "数が多すぎる（5×10^7 個まで）");
        std::vector<long long> nxt;
        nxt.reserve(cur.size() * ds.size());
        for (long long x : cur)
            for (int d : ds) nxt.push_back(x * 10 + d);
        cur.swap(nxt);
    }
    return res;
}
#endif  // GG_DIGITS
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace dig_test {

#define DIG_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[digits] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
static std::string brute_base(unsigned long long u, unsigned long long b) {  // 正の進数（桁は 0-9A-Z。36 まで）
    static const char* D = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    if (u == 0) return "0";
    std::string s;
    while (u) s = D[u % b] + s, u /= b;
    return s;
}
static bool brute_pal(const std::string& s) { return std::string(s.rbegin(), s.rend()) == s; }
static Int128 parse128(const std::string& s) {  // 10 進（先頭の - / + 可）
    size_t i = (s[0] == '-' || s[0] == '+') ? 1 : 0;
    Int128 v = 0;
    for (; i < s.size(); ++i) v = v * 10 + (s[i] - '0');
    return s[0] == '-' ? -v : v;
}

static void test_digits_and_base(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd_ll(rng);
        const long long base = it % 3 == 0 ? 10 : (it % 3 == 1 ? rnd(rng, 2, 36) : rnd(rng, 2, LMAX));
        const unsigned long long u = n < 0 ? 0ULL - (unsigned long long)n : (unsigned long long)n;
        std::vector<int> ds;
        unsigned long long t = u, s = 0;
        do ds.insert(ds.begin(), (int)(t % (unsigned long long)base)), s += t % (unsigned long long)base, t /= (unsigned long long)base;
        while (t);
        DIG_CHECK(digit_sum(n, base) == (long long)s && num_digits(n, base) == (int)ds.size());
        if (base <= 2147483647LL) {
            DIG_CHECK(digits_of(n, base) == ds);
            if (u <= (unsigned long long)LMAX) DIG_CHECK(from_digits(ds, base) == (long long)u);
        }
        if (base <= 36) {
            const std::string b = brute_base(u, (unsigned long long)base);
            DIG_CHECK(to_base(n, (int)base) == (n < 0 ? "-" : "") + b);
            DIG_CHECK(from_base(to_base(n, (int)base), (int)base) == n);
            std::string low = b;
            for (char& c : low) c = (char)std::tolower(c);
            if (n >= 0) DIG_CHECK(from_base(low, (int)base) == n);
        }
        // 負の進数：値・桁の範囲・先頭が 0 でないこと
        const int nb = -(int)rnd(rng, 2, 36);
        const std::string ns = to_base(n, nb);
        Int128 v = 0;
        for (char c : ns) {
            const int d = c <= '9' ? c - '0' : c - 'A' + 10;
            DIG_CHECK(0 <= d && d < -nb);
            v = v * nb + d;
        }
        DIG_CHECK(v == n && (ns == "0" || ns[0] != '0') && from_base(ns, nb) == n);
    }
    DIG_CHECK(to_base(-9, -2) == "1011" && to_base(0, -2) == "0" && to_base(255, 16) == "FF" && to_base(-255, 16) == "-FF");
    DIG_CHECK(to_base(std::numeric_limits<long long>::min(), 2) == "-1" + std::string(63, '0') && from_base("-1" + std::string(63, '0'), 2) == std::numeric_limits<long long>::min());
    DIG_CHECK(digits_of(0) == std::vector<int>{0} && num_digits(0) == 1 && digit_sum(-123) == 6);
    // from_base_cap / base_value_leq
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        const int len = (int)rnd(rng, 1, 25);
        const long long base = it % 2 ? rnd(rng, 2, 36) : rnd(rng, 10, it % 4 == 0 ? LMAX : 100000);
        std::string s;
        for (int i = 0; i < len; ++i) s += "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[rng() % std::min<long long>(base, it % 3 ? 10 : 36)];
        const Int128 SAT = (Int128)1 << 100;
        Int128 v = 0;
        for (char c : s) {
            const int d = c <= '9' ? c - '0' : c - 'A' + 10;
            v = v > SAT / base ? SAT : std::min<Int128>(v * base + d, SAT);  // 愚直解も 2^100 で止める
        }
        long long cap = rnd_ll(rng);
        cap = cap < 0 ? (cap == std::numeric_limits<long long>::min() ? 0 : -cap) : cap;
        DIG_CHECK((Int128)from_base_cap(s, base, cap) == std::min<Int128>(v, cap));
        long long lim = it % 2 ? rnd_ll(rng) : (v <= LMAX ? (long long)v + rnd(rng, -1, 1) : LMAX);
        DIG_CHECK(base_value_leq(s, base, lim) == (v <= lim));
    }
}

static void test_palindrome(std::mt19937_64& rng, int& cases) {
    const long long M = 2000000;
    std::vector<long long> pals;
    for (long long x = 0; x <= M; ++x)
        if (brute_pal(std::to_string(x))) pals.push_back(x);
    long long cnt = 0;  // [0, x] の回文数の数
    size_t k = 0;
    for (long long x = 0; x <= M; ++x, ++cases) {
        const bool p = k < pals.size() && pals[k] == x;
        cnt += p;
        DIG_CHECK(is_palindrome_num(x) == p && count_palindromes_up_to(x) == cnt && count_palindromes_below(x) == cnt - p);
        if (x < 1900000) DIG_CHECK(ceil_palindrome(x) == pals[k]);  // k はまだ x 以上で最初の回文数の番号
        if (p) {
            DIG_CHECK(palindrome_rank(x) == (long long)k && kth_palindrome((long long)k) == std::to_string(x));
            ++k;
        }
        DIG_CHECK(floor_palindrome(x) == pals[k - 1]);  // k − 1 は x 以下で最後の回文数の番号
        if (x < 3000) DIG_CHECK(is_palindrome_num(x, 2) == brute_pal(brute_base((unsigned long long)x, 2)));
    }
    // 大きい数：kth と rank、ceil / floor のあいだに回文数がないこと
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long kk = rnd(rng, 0, it % 2 ? 1000000000LL : 1900000000LL);
        const std::string s = kth_palindrome(kk);
        DIG_CHECK(brute_pal(s) && (s.size() == 1 || s[0] != '0'));
        DIG_CHECK(palindrome_rank(std::stoll(s)) == kk);
        const long long n = rnd(rng, 0, it % 3 ? 1000000000000000LL : 9000000000000000000LL);
        const long long c = ceil_palindrome(n), f = floor_palindrome(n);
        DIG_CHECK(c >= n && f <= n && is_palindrome_num(c) && is_palindrome_num(f));
        DIG_CHECK(count_palindromes_below(c) == count_palindromes_below(n) && count_palindromes_up_to(f) == count_palindromes_up_to(n));
        const long long half = rnd(rng, 1, 999999999);
        std::string hs = std::to_string(half);
        DIG_CHECK(std::to_string(make_palindrome(half, true)) == hs + std::string(hs.rbegin() + 1, hs.rend()) && std::to_string(make_palindrome(half, false)) == hs + std::string(hs.rbegin(), hs.rend()));
    }
    // 長い回文数（文字列）：桁数ごとの個数と合うこと
    long long base = 0;
    for (int len = 1; len <= 36; ++len) {
        const long long c = count_palindromes_of_length(len);
        DIG_CHECK(kth_palindrome(base).size() == (size_t)len && kth_palindrome(base + c - 1) == std::string((size_t)len, '9'));
        if (__builtin_add_overflow(base, c, &base)) break;
        ++cases;
    }
    DIG_CHECK(kth_palindrome(45) == "363" && ceil_palindrome(LMAX - 1000000000000LL) == 9223371037301733229LL && floor_palindrome(LMAX) == 9223372036302733229LL);
    // 文字列の回文数：小さい数は全部、long long の範囲は long long 版と、長い数は性質で
    {
        size_t kk = 0;
        Int128 c = 0;
        for (long long x = 0; x <= 300000; ++x, ++cases) {
            const std::string s = std::to_string(x);
            while (pals[kk] < x) ++kk;
            const bool p = pals[kk] == x;
            c += p;
            DIG_CHECK(ceil_palindrome_str(s) == std::to_string(pals[kk]) && floor_palindrome_str(s) == std::to_string(pals[kk - !p]));
            DIG_CHECK(count_palindromes_up_to_str(s) == c && count_palindromes_below_str(s) == c - p && is_palindrome_str(s) == p);
            if (x % 1000 == 7) DIG_CHECK(ceil_palindrome_str("00" + s) == ceil_palindrome_str(s) && count_palindromes_up_to_str("+0" + s) == c);
        }
    }
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, it % 2 ? 1000000000000LL : 9000000000000000000LL);
        const std::string s = std::to_string(n);
        DIG_CHECK(ceil_palindrome_str(s) == std::to_string(ceil_palindrome(n)) && floor_palindrome_str(s) == std::to_string(floor_palindrome(n)));
        DIG_CHECK(count_palindromes_up_to_str(s) == count_palindromes_up_to(n) && count_palindromes_below_str(s) == count_palindromes_below(n));
        // 長い数（0・1・9 が多いと境目を通る）
        const int len = (int)rnd(rng, 1, 75);
        std::string t(1, char('1' + rng() % 9));
        const char* pool = it % 3 == 0 ? "019" : (it % 3 == 1 ? "0123456789" : "09");
        const int psz = (int)std::char_traits<char>::length(pool);
        while ((int)t.size() < len) t += pool[rng() % psz];
        if (it % 5 == 0) t = make_palindrome_str(t.substr(0, (len + 1) / 2), len % 2) + "";  // 回文そのもの
        if (it % 7 == 0) t = std::string((size_t)len, '9');
        if (it % 11 == 0) t = "1" + std::string((size_t)len - 1, '0');
        const std::string ce = ceil_palindrome_str(t), fl = floor_palindrome_str(t);
        DIG_CHECK(brute_pal(ce) && brute_pal(fl) && str_num_compare(fl, t) <= 0 && str_num_compare(t, ce) <= 0);
        DIG_CHECK((ce.size() == 1 || ce[0] != '0') && (fl.size() == 1 || fl[0] != '0'));
        if (t.size() <= 74) {
            const Int128 up = count_palindromes_up_to_str(t), below = count_palindromes_below_str(t);
            DIG_CHECK(count_palindromes_below_str(ce) == below && count_palindromes_up_to_str(fl) == up);  // あいだに回文数がない
            DIG_CHECK(up - below == (brute_pal(t) ? 1 : 0) && count_palindromes_up_to_str(ce) == below + 1);
        }
    }
    DIG_CHECK(ceil_palindrome_str("100000000000000000000000") == "100000000000000000000001" && floor_palindrome_str("100000000000000000000000") == "99999999999999999999999");
    DIG_CHECK(ceil_palindrome_str("123456789012345678901234567890") == "123456789012346643210987654321" && ceil_palindrome_str("199") == "202");
    DIG_CHECK(count_palindromes_up_to_str(kth_palindrome(1000000000000000000LL)) == (Int128)1000000000000000001LL && make_palindrome_str("123", true) == "12321");
    DIG_CHECK(count_palindromes_up_to_str(std::string(75, '9')) > 0 && floor_palindrome_str("0") == "0" && ceil_palindrome_str("000") == "0");
}

static void test_concat_and_strings(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long a = rnd(rng, 0, it % 2 ? 1000 : LMAX), b = rnd(rng, 0, it % 3 ? 1000 : LMAX), m = rnd(rng, 1, it % 2 ? 1000000007 : LMAX);
        const std::string s = std::to_string(a) + std::to_string(b);
        DIG_CHECK(concat_value(a, b) == parse128(s));
        DIG_CHECK(concat_mod(a, b, m) == (long long)(parse128(s) % m) && str_mod(s, m) == (long long)(parse128(s) % m));
        const long long k = rnd(rng, 0, 12), x = rnd(rng, 0, it % 2 ? 99 : 999999);
        long long r = 0, rep = 0;
        const long long p = pow10_mod(num_digits(x), m);
        for (long long i = 0; i < k; ++i) r = (long long)(((Int128)r * p + x) % m), rep = (long long)(((Int128)rep * 10 + 1) % m);
        DIG_CHECK(repeat_concat_mod(x, k, m) == r && repunit_mod(k, m) == rep);
        {
            std::string rs;
            for (long long i = 0; i < k; ++i) rs += std::to_string(x);
            if (rs.size() <= 38) DIG_CHECK(repeat_concat_value(x, k) == (rs.empty() ? 0 : parse128(rs)));
        }
        // base 進数の文字列の mod（36 進 20 桁まで）
        {
            const int bs = (int)rnd(rng, 2, 36), len = (int)rnd(rng, 1, 20);
            std::string t;
            Int128 val = 0;
            for (int i = 0; i < len; ++i) {
                const int d = (int)(rng() % bs);
                t += (d < 10 ? char('0' + d) : char((rng() % 2 ? 'A' : 'a') + d - 10));
                val = val * bs + d;
            }
            const bool neg = rng() % 2;
            DIG_CHECK(str_mod((neg ? "-" : "") + t, m, bs) == (long long)(((neg ? -val : val) % m + m) % m));
            const long long y = rnd_ll(rng);
            DIG_CHECK(str_mod(to_base(y, bs), m, bs) == (long long)(((Int128)y % m + m) % m));
        }
        // 繰り上がり
        {
            const int bs = it % 2 ? 10 : (int)rnd(rng, 2, 36), len = (int)rnd(rng, 0, 8);
            std::vector<long long> co(len);
            Int128 val = 0, pw = 1;
            for (int i = 0; i < len; ++i) co[i] = rnd(rng, it % 3 == 0 ? -1000 : 0, it % 4 == 0 ? 1000000000000LL : 50), val += co[i] * pw, pw *= bs;
            if (val >= 0) {
                std::string expect;
                for (Int128 v = val; v > 0; v /= bs) expect += "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[(int)(v % bs)];
                if (expect.empty()) expect = "0";
                std::reverse(expect.begin(), expect.end());
                DIG_CHECK(carry_to_string(co, bs) == expect);
            }
            // 桁の畳み込み = 掛け算
            const long long A = rnd(rng, 0, LMAX), B = rnd(rng, 0, LMAX);
            std::vector<int> da = digits_of(A), db = digits_of(B);
            std::reverse(da.begin(), da.end()), std::reverse(db.begin(), db.end());
            std::vector<long long> cv(da.size() + db.size() - 1, 0);
            for (size_t i = 0; i < da.size(); ++i)
                for (size_t j = 0; j < db.size(); ++j) cv[i + j] += da[i] * db[j];
            std::ostringstream os;
            os << (Int128)A * B;
            DIG_CHECK(carry_to_string(cv) == os.str());
        }
        // 文字列の比較・mod（36 桁まで。符号・先頭の 0）
        auto rnds = [&]() {
            std::string t = std::to_string(rnd(rng, 0, 9));
            const int len = (int)rnd(rng, 0, 35);
            for (int i = 0; i < len; ++i) t += char('0' + rng() % 10);
            if (rng() % 4 == 0) t = "00" + t;
            if (rng() % 2) t = "-" + t;
            return t;
        };
        const std::string u = rnds(), v = rnds();
        const Int128 iu = parse128(u), iv = parse128(v);
        DIG_CHECK(str_num_compare(u, v) == (iu > iv) - (iu < iv) && str_num_compare(u, u) == 0);
        DIG_CHECK(str_mod(u, m) == (long long)((iu % m + m) % m));
        // 部分文字列
        std::string w;
        for (int i = 0; i < 40; ++i) w += char('0' + rng() % 10);
        const DigitPrefix P(w, m);
        const int l = (int)rnd(rng, 0, 40), rr = (int)rnd(rng, l, std::min(40, l + 18));
        DIG_CHECK(substr_value(w, l, rr) == (rr > l ? std::stoll(w.substr(l, rr - l)) : 0) && P.value(l, rr) == substr_value(w, l, rr) % m);
        const int l2 = (int)rnd(rng, 0, 40), r2 = (int)rnd(rng, l2, 40);
        DIG_CHECK(P.value(l2, r2) == str_mod(r2 > l2 ? w.substr(l2, r2 - l2) : "0", m));
    }
    DIG_CHECK(str_num_compare("-0", "0") == 0 && str_num_compare("007", "7") == 0 && str_num_compare("-10", "-9") == -1);
    DIG_CHECK(repeat_concat_value(9, 38) == parse128(std::string(38, '9')) && repeat_concat_value(0, LMAX) == 0 && repeat_concat_value(123, 0) == 0);
    DIG_CHECK(carry_to_string({12, 3}) == "42" && carry_to_string(std::vector<int>{}) == "0" && carry_to_string({0, 0, 0}) == "0" && carry_to_string({-1, 1}) == "9");
}

static void test_counting(std::mt19937_64& rng, int& cases) {
    Int128 td = 0, dsum_total = 0;
    Int128 cnt[10] = {0};
    for (long long x = 1; x <= 300000; ++x, ++cases) {
        const std::string s = std::to_string(x);
        td += (Int128)s.size();
        for (char c : s) cnt[c - '0']++, dsum_total += c - '0';
        DIG_CHECK(total_digits_up_to(x) == td && digit_sum_total_up_to(x) == dsum_total);
        if (x % 13 == 0 || x < 2000)
            for (int d = 0; d <= 9; ++d) DIG_CHECK(count_digit_up_to(x, d) == cnt[d]);
    }
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {  // 大きい N は 1 つ前との差
        const long long N = rnd(rng, 2, it % 2 ? 1000000000000LL : LMAX);
        const std::string s = std::to_string(N);
        DIG_CHECK(total_digits_up_to(N) - total_digits_up_to(N - 1) == (Int128)s.size());
        long long dsum = 0;
        for (char c : s) dsum += c - '0';
        DIG_CHECK(digit_sum_total_up_to(N) - digit_sum_total_up_to(N - 1) == dsum);
        const int d = (int)rnd(rng, 0, 9);
        DIG_CHECK(count_digit_up_to(N, d) - count_digit_up_to(N - 1, d) == (Int128)std::count(s.begin(), s.end(), char('0' + d)));
        // 桁数ごとの区間
        const long long l = rnd(rng, 0, N), r = rnd(rng, l, LMAX);
        long long x = l;
        for (const DigitLengthBlock& b : digit_length_blocks(l, r)) {
            DIG_CHECK(b.lo == x && b.lo < b.hi && num_digits(b.lo) == b.len && num_digits(b.hi - 1) == b.len);
            x = b.hi;
        }
        DIG_CHECK(x == std::max(l, r) || (l == r && x == l));
    }
    DIG_CHECK(total_digits_up_to(LMAX) - total_digits_up_to(LMAX - 1) == 19);
    // 使える数字が決まった数
    for (int it = 0; it < 200; ++it, ++cases) {
        std::vector<int> ds;
        for (int d = 0; d <= 9; ++d)
            if (rng() % 3 == 0) ds.push_back(d);
        const int L = (int)rnd(rng, 1, 5);
        std::vector<long long> expect;
        for (long long y = 0; y < pow10_ll(L); ++y) {
            bool ok = true;
            for (char c : std::to_string(y)) ok &= std::find(ds.begin(), ds.end(), c - '0') != ds.end();
            if (ok) expect.push_back(y);
        }
        DIG_CHECK(numbers_from_digits(ds, L) == expect);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_digits_and_base(rng, cases);
    test_palindrome(rng, cases);
    test_concat_and_strings(rng, cases);
    test_counting(rng, cases);
    std::cout << "[digits] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[digits] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 1000000;
    std::vector<long long> xs(n);
    for (auto& x : xs) x = rnd(rng, 0, LMAX);
    auto t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (long long x : xs) s1 += digit_sum(x) + num_digits(x);
    const long long t_ds = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for (int i = 0; i < n / 4; ++i) s2 += (long long)to_base(xs[i], -2).size() + from_base(to_base(xs[i], 7), 7) % 3;
    const long long t_base = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s3 = 0;
    for (int i = 0; i < n / 4; ++i) s3 += ceil_palindrome(xs[i] / 2) % 7 + count_palindromes_up_to(xs[i]) % 7 + (long long)kth_palindrome(xs[i] / 10).size();
    const long long t_pal = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s4 = 0;
    for (int i = 0; i < n / 4; ++i) s4 += (long long)(count_digit_up_to(xs[i], i % 10) % 7) + repeat_concat_mod(xs[i] % 1000000, xs[i], 998244353);
    const long long t_cnt = ms(t0);
    std::string big(1000000, '7');
    t0 = std::chrono::steady_clock::now();
    const DigitPrefix P(big, 998244353);
    long long s5 = str_mod(big, 1000000007);
    for (int i = 0; i < n; ++i) s5 += P.value((int)(xs[i] % 500000), (int)(500000 + xs[i] % 500000));
    const long long t_str = ms(t0);
    std::cout << "[digits] benchmark: digit_sum + num_digits × 10^6 " << t_ds << " ms, to_base(−2) + to_base / from_base(7) × 2.5×10^5 " << t_base
              << " ms, 回文数 3 種 × 2.5×10^5 " << t_pal << " ms\n"
              << "[digits] benchmark: count_digit_up_to + repeat_concat_mod × 2.5×10^5 " << t_cnt << " ms, 10^6 桁の str_mod + DigitPrefix の前処理と 10^6 回の value " << t_str
              << " ms  (check " << (s1 + s2 + s3 + s4 + s5) % 7 << ")\n";
}

}  // namespace dig_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace dig_demo {

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
#define SHOW(expr, comment) dig_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    dig_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("桁・進数・回文数（digits） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  digits_of / digit_sum / num_digits / from_digits              桁（base を指定できる）\n"
              << "  to_base / from_base（±2〜36 進）/ from_base_cap / base_value_leq（36 を超える進数も）\n"
              << "  is_palindrome_num / kth_palindrome / palindrome_rank / count_palindromes_up_to / _below / _of_length\n"
              << "  ceil_palindrome / floor_palindrome / make_palindrome           回文数\n"
              << "  ceil_palindrome_str / floor_palindrome_str / count_palindromes_up_to_str / _below_str / make_palindrome_str / is_palindrome_str   大きな回文数\n"
              << "  concat_value / repeat_concat_value / concat_mod / repeat_concat_mod / repunit_mod / pow10_mod   数をつなげる\n"
              << "  str_mod / str_num_compare / substr_value / DigitPrefix / carry_to_string   大きな数の文字列\n"
              << "  digit_length_blocks / total_digits_up_to / count_digit_up_to / digit_sum_total_up_to   桁の数え上げ\n"
              << "  numbers_from_digits                                            使える数字が決まった数の列挙\n";

    api("【1】 digits_of(n, base) / digit_sum(n, base) / num_digits(n, base) / from_digits(ds, base)", "vector<int> / long long / int",
        "各桁（上の桁から）/ 桁和 / 桁数。base は 2 以上（digits_of は 2^31 − 1 まで）。負の数は絶対値で", "O(桁数)");
    SHOW(digits_of(2026), "上の桁から");
    SHOW(digits_of(10, 2), "2 進数の桁");
    SHOW(digit_sum(999), "桁和");
    SHOW(num_digits(1000000), "桁数");
    SHOW(from_digits({1, 0, 1}, 2), "桁の並びから（2 進の 101 = 5）");
    {
        // ABC083 B: 1〜N で桁和が A 以上 B 以下の数の和
        const long long N = 20, A = 2, B = 5;
        long long s = 0;
        for (long long x = 1; x <= N; ++x)
            if (A <= digit_sum(x) && digit_sum(x) <= B) s += x;
        note("ABC083 B（Some Sums）: digit_sum で桁和を調べる");
        std::cout << "    N = 20, A = 2, B = 5 → " << s << "\n";
    }

    api("【2】 to_base(n, base) / from_base(s, base)", "string / long long", "2〜36 進と −2〜−36 進の文字列との変換（10 以上の桁は A〜Z。0・負の数も）", "O(桁数)");
    SHOW(to_base(-9, -2), "ABC105 C: −2 進数");
    SHOW(to_base(255, 16), "16 進数");
    SHOW(from_base("ff", 16), "小文字も読める");
    SHOW(to_base(-10, 3), "正の進数の負の数は - が付く");
    {
        // 典型90 067: 8 進数の N を 9 進数にして 8 を 5 に置き換える、を K 回
        std::string s = "21";
        const int K = 1;
        for (int k = 0; k < K; ++k) {
            s = to_base(from_base(s, 8), 9);
            for (char& c : s)
                if (c == '8') c = '5';
        }
        note("典型90 067（Base 8 to 9）: from_base(s, 8) → to_base(·, 9) → 8 を 5 に");
        std::cout << "    N = 21, K = 1 → " << s << "\n";
    }
    {
        // ABC336 C: 各桁が偶数の数の N 番目（0 が 1 番目）= N − 1 の 5 進数の各桁を 2 倍
        const long long N = 8;
        std::string s = to_base(N - 1, 5);
        for (char& c : s) c = char('0' + (c - '0') * 2);
        note("ABC336 C（Even Digits）: N − 1 を 5 進数にして各桁を 2 倍");
        std::cout << "    N = 8 → " << s << "\n";
    }
    {
        // ABC186 C: N 以下で 10 進でも 8 進でも 7 を含まない数
        const long long N = 20;
        long long c = 0;
        for (long long x = 1; x <= N; ++x) {
            const std::string a = to_base(x, 10), b = to_base(x, 8);
            c += a.find('7') == std::string::npos && b.find('7') == std::string::npos;
        }
        note("ABC186 C（Unlucky 7）: to_base(x, 8) で 8 進数の文字列にして 7 を探す");
        std::cout << "    N = 20 → " << c << "\n";
    }

    api("【3】 from_base_cap(s, base, cap) / base_value_leq(s, base, limit)", "long long / bool",
        "base は 2 以上いくつでも（10^18 でも）。値が cap・limit を超えるかを溢れずに調べる", "O(|s|)");
    SHOW(from_base_cap("99", 1000000000000LL, 1000000000000000000LL), "99 × 10^12 + 9 は上限 10^18 で止まる");
    SHOW(base_value_leq("22", 4, 10), "22 (4 進) = 10 ≤ 10");
    {
        // ABC192 D: X を n 進数（n > X の最大の桁）として読んだ値が M 以下になる n の個数（二分探索）
        const std::string X = "22";
        const long long M = 10;
        long long dmax = 0;
        for (char c : X) dmax = std::max<long long>(dmax, c - '0');
        long long ans;
        if (X.size() == 1) {
            ans = (X[0] - '0') <= M ? 1 : 0;
        } else {
            long long lo = dmax, hi = M + 1;  // lo は条件を満たす最大の進数（dmax は「満たさない」の番兵）
            while (hi - lo > 1) {
                const long long mid = lo + (hi - lo) / 2;
                (base_value_leq(X, mid, M) ? lo : hi) = mid;
            }
            ans = lo - dmax;
        }
        note("ABC192 D（Base n）: base_value_leq(X, n, M) で二分探索（X が 1 桁なら別扱い）");
        std::cout << "    X = 22, M = 10 → " << ans << "（3 進と 4 進）\n";
    }

    api("【4】 kth_palindrome(k) / palindrome_rank(n) / count_palindromes_up_to(n) / ceil_palindrome(n) / floor_palindrome(n)", "string / long long",
        "回文数（10 進）。k 番目は 0 から数える（0 番目が 0）。kth_palindrome は 10^18 を超えるので文字列", "O(桁数)");
    SHOW(kth_palindrome(45), "ABC363 D: 46 番目（N = 46 なら kth_palindrome(N − 1)）");
    SHOW(kth_palindrome(1000000000000000000LL), "とても大きい k");
    SHOW(palindrome_rank(363), "363 は 45 番目");
    SHOW(count_palindromes_up_to(1000), "0〜1000 の回文数");
    SHOW(ceil_palindrome(1234), "1234 以上で最小");
    SHOW(floor_palindrome(1000), "1000 以下で最大");
    SHOW(make_palindrome(123, true), "左半分から（真ん中を共有）");
    SHOW(is_palindrome_num(9, 2), "2 進数で 1001");
    note("long long を超える数は文字列版（_str）。count_palindromes_up_to_str は Int128 を返す（75 桁まで）");
    SHOW(ceil_palindrome_str("123456789012345678901234567890"), "30 桁");
    SHOW(floor_palindrome_str("1000000000000000000000"), "10^21 以下で最大");
    SHOW(count_palindromes_up_to_str("100000000000000000000000000000"), "0〜10^29 の回文数");
    SHOW(make_palindrome_str("abc", false), "文字列でも作れる");
    SHOW(is_palindrome_str("racecar"), "文字列が回文か");

    api("【5】 concat_value(a, b) / concat_mod(a, b, m) / repeat_concat_mod(x, k, m) / repunit_mod(k, m)", "Int128 / long long",
        "数をつなげた値 / その mod / x を k 回つなげた数の mod（k ≤ 10^18）/ 11…1 の mod", "O(1) / O(log k)");
    SHOW(concat_value(123, 45), "12345");
    SHOW(repeat_concat_value(2026, 3), "202620262026");
    SHOW(concat_mod(1000000000LL, 1000000000LL, 998244353), "10^9 の後ろに 10^9（19 桁）を mod で");
    SHOW(repeat_concat_mod(5, 5, 998244353), "ABC357 D: 5 を 5 回 → 55555");
    SHOW(repeat_concat_mod(10, 3, 1000000007), "101010");
    SHOW(repunit_mod(1000000000000000000LL, 998244353), "1 を 10^18 個並べた数 mod p");

    api("【6】 str_mod(s, m, base) / str_num_compare(a, b) / substr_value(s, l, r) / DigitPrefix P(s, m) / carry_to_string(c, base)", "long long / int / string",
        "long long に入らない数（文字列）の mod・大小 / 部分文字列の値（P.value(l, r) は mod で O(1)）/ 係数の繰り上がり", "O(|s|) / O(1)");
    SHOW(str_mod("123456789012345678901234567890", 9), "桁和と同じ（9 で割った余り）");
    SHOW(str_num_compare("100000000000000000000", "-5"), "符号も見る");
    SHOW(str_num_compare("0012", "12"), "先頭の 0 は無視");
    SHOW(substr_value("20261010", 4, 8), "部分文字列 [4, 8)");
    RUN(DigitPrefix P("31415926535", 1000000007));
    SHOW(P.value(2, 7), "\"41592\"");
    SHOW(str_mod("ZZ", 1000, 36), "36 進数の ZZ = 1295 の mod");
    note("carry_to_string：桁ごとの係数（下の桁から）を繰り上げて文字列に。大きな数の掛け算を畳み込みでするときの最後の一歩");
    SHOW(carry_to_string({18, 27, 28, 13, 4}), "123 × 456：桁 {3, 2, 1} と {6, 5, 4} の畳み込み → 56088");
    SHOW(carry_to_string({1, 1, 0, 1}, 2), "2 進数の 1011");

    api("【7】 digit_length_blocks(l, r) / total_digits_up_to(N) / count_digit_up_to(N, d) / digit_sum_total_up_to(N)", "vector / Int128",
        "[l, r) を桁数ごとの区間に分ける / 1〜N の桁数の合計・数字 d の個数・桁和の合計", "O(桁数)");
    SHOW(total_digits_up_to(100), "\"123…100\" の長さ（9 + 180 + 3）");
    SHOW(count_digit_up_to(100, 0), "1〜100 の 0 の個数（10, 20, …, 100）");
    SHOW(digit_sum_total_up_to(10), "1 + 2 + … + 9 + 1");
    {
        // 典型90 082: Σ_{x=L}^{R} x × (x の桁数) mod 10^9+7
        const long long L = 98, R = 100, MOD = 1000000007;
        long long ans = 0;
        for (const DigitLengthBlock& b : digit_length_blocks(L, R + 1))
            ans = (ans + mod_mul(b.len, arith_sum_mod(b.lo, 1, b.hi - b.lo, MOD), MOD)) % MOD;
        note("典型90 082（Counting Numbers）: 桁数ごとの区間に分けて、桁数 × 区間の和（arith_sum_mod）");
        std::cout << "    L = 98, R = 100 → " << ans << "（98·2 + 99·2 + 100·3）\n";
    }

    api("【8】 numbers_from_digits(ds, max_len)", "vector<long long>", "数字 ds だけでできた 1〜max_len 桁の数（小さい順。0 が ds にあれば 0 も）", "O(個数)");
    SHOW(numbers_from_digits({1, 3}, 2), "1, 3, 11, 13, 31, 33");
    {
        // ABC114 C: N 以下で、数字が 3・5・7 だけで、3・5・7 をすべて含む数（七五三数）
        const long long N = 575;
        long long c = 0;
        for (long long x : numbers_from_digits({3, 5, 7}, num_digits(N))) {
            const std::string s = std::to_string(x);
            if (x <= N && s.find('3') != std::string::npos && s.find('5') != std::string::npos && s.find('7') != std::string::npos) ++c;
        }
        note("ABC114 C（755）: numbers_from_digits({3, 5, 7}, 桁数) から条件を満たすものを数える");
        std::cout << "    N = 575 → " << c << "\n";
    }
}

}  // namespace dig_demo

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
    dig_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    dig_test::benchmark();    // 速度計測
    dig_demo::demo();         // 使い方デモ
    return 0;
}
