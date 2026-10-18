// =============================================================================
//  ゲーム（mex・Grundy 数：状態の DAG・石取り（取れる数が L〜R・集合 S、周期も）、分かれる手、Nim・逆形 Nim、閉路のあるゲームの後退解析）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 10 ゲーム
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  「動けなくなった方が負け」の 2 人ゲーム。いくつかのゲームを並べ、手番ではどれか 1 つを 1 手進めるとき
//  先手必勝 ⇔ 各ゲームの Grundy 数の XOR ≠ 0（Grundy 数 = 1 手で行ける状態の Grundy 数の mex）
//
//  機能一覧                                                                    計算量
//    mex(vs)                          vs に含まれない最小の 0 以上の整数                            O(|vs|)
//    grundy_dag(n, moves)             状態 0〜n−1、moves[v] = v から 1 手で行ける状態。各状態の Grundy 数  O(n + 辺)
//    nim_first_wins(piles)            Nim（最後に取った方が勝ち）で先手必勝か（XOR ≠ 0）
//    misere_nim_first_wins(piles)     逆形 Nim（最後に取った方が負け）で先手必勝か
//    range_take_grundy(x, L, R)       1 回に L〜R 個取れる石取りの x 個の山の Grundy 数               O(1)
//    subtraction_game_grundy(S, N)    1 回に取れる数が集合 S の石取りの Grundy 数 g(0)〜g(N)          O(N |S|)
//    SubtractionGame G(S)             同じく x が巨大なとき：周期を見つけて G.grundy(x)。G.preperiod() / G.period()
//    grundy_dag_multi(n, moves)       1 手で複数のゲームに分かれるとき（moves[v] = 手の一覧、手 = 行き先の状態の列。値は XOR）  O(n + 行き先の総数)
//    game_outcome(n, moves)           閉路があってもよいゲームの勝ち負け（後退解析）：GameOutcome{result, steps}         O(n + 辺)
//                                     result[v] = GAME_WIN / GAME_LOSE / GAME_DRAW（v で手番の人から）、steps[v] = 決着までの手数
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    game_test::random_test()  … 全部の手を試すメモ化探索との突き合わせ
//    game_test::benchmark()    … 速度計測
//    game_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_GAME
#define GG_GAME
// ---- ゲーム（数学第二世代）  mex・Grundy 数（DAG・石取り・分かれる手）・Nim・閉路のあるゲームの後退解析 ----
// 「動けなくなった方が負け」の 2 人ゲームで、いくつかのゲームを並べたもの（どれか 1 つで 1 手）の先手必勝 ⇔ 各 Grundy 数の XOR ≠ 0
template <class T>
int mex(const std::vector<T>& vs) {  // vs に含まれない最小の 0 以上の整数（負の値は無視）  O(|vs|)
    std::vector<char> seen(vs.size() + 1, 0);
    for (const T& v : vs) {
        if constexpr (std::is_signed<T>::value)
            if (v < 0) continue;
        if ((unsigned long long)v <= vs.size()) seen[(size_t)v] = 1;
    }
    int m = 0;
    while (seen[(size_t)m]) ++m;
    return m;
}
// 状態 0〜n−1、moves[v] = v から 1 手で行ける状態。各状態の Grundy 数（閉路があれば止める）  O(n + 辺の数)
inline std::vector<int> grundy_dag(int n, const std::vector<std::vector<int>>& moves) {
    if (n < 0 || (int)moves.size() != n) math_internal::fail("grundy_dag", "moves の大きさは n");
    std::vector<int> g((size_t)n, -1), st;
    std::vector<char> state((size_t)n, 0);  // 0 = まだ、1 = 調べている途中、2 = 決まった
    std::vector<size_t> it((size_t)n, 0);
    std::vector<char> seen;
    for (int s = 0; s < n; ++s) {
        if (state[(size_t)s]) continue;
        st.push_back(s), state[(size_t)s] = 1;
        while (!st.empty()) {  // 行き先を全部決めてから自分を決める（再帰を使わない深さ優先）
            const int v = st.back();
            const std::vector<int>& mv = moves[(size_t)v];
            if (it[(size_t)v] < mv.size()) {
                const int w = mv[it[(size_t)v]++];
                if (w < 0 || w >= n) math_internal::fail("grundy_dag", "行き先の番号は 0〜n−1");
                if (state[(size_t)w] == 1) math_internal::fail("grundy_dag", "閉路がある（同じ状態に戻れる）");
                if (state[(size_t)w] == 0) st.push_back(w), state[(size_t)w] = 1;
                continue;
            }
            seen.assign(mv.size() + 1, 0);
            for (int w : mv)
                if ((size_t)g[(size_t)w] <= mv.size()) seen[(size_t)g[(size_t)w]] = 1;
            int m = 0;
            while (seen[(size_t)m]) ++m;
            g[(size_t)v] = m, state[(size_t)v] = 2, st.pop_back();
        }
    }
    return g;
}
// Nim：好きな山から 1 個以上取る。最後に取った方が勝ち ⇔ XOR ≠ 0 なら先手必勝
template <class T>
bool nim_first_wins(const std::vector<T>& piles) {
    unsigned long long x = 0;
    for (const T& a : piles) x ^= (unsigned long long)a;
    return x != 0;
}
// 逆形 Nim（最後に取った方が負け）：全部の山が 1 以下なら「1 の山が偶数個」、そうでなければ XOR ≠ 0 なら先手必勝
template <class T>
bool misere_nim_first_wins(const std::vector<T>& piles) {
    unsigned long long x = 0;
    bool big = false;
    long long ones = 0;
    for (const T& a : piles) x ^= (unsigned long long)a, big = big || a > 1, ones += a == 1;
    return big ? x != 0 : ones % 2 == 0;
}
// 1 回に L〜R 個取れる石取り（L ≥ 1）の、x 個の山の Grundy 数 ⌊(x mod (L + R)) / L⌋  O(1)
inline long long range_take_grundy(long long x, long long L, long long R) {
    if (x < 0 || L < 1 || R < L) math_internal::fail("range_take_grundy", "x ≥ 0、1 ≤ L ≤ R で使う");
    return (long long)(((UInt128)x % ((UInt128)L + (UInt128)R)) / (UInt128)L);
}
// 1 回に取れる数が集合 S（S の要素は 1 以上）の石取りの Grundy 数 g(0), …, g(N)  O(N |S|)
inline std::vector<int> subtraction_game_grundy(const std::vector<long long>& S, int N) {
    for (long long s : S)
        if (s < 1) math_internal::fail("subtraction_game_grundy", "取れる数は 1 以上");
    if (N < 0) return {};
    std::vector<int> g((size_t)N + 1, 0);
    std::vector<int> mark(S.size() + 2, -1);
    for (int x = 0; x <= N; ++x) {
        for (long long s : S)
            if (s <= x && (size_t)g[(size_t)(x - s)] <= S.size()) mark[(size_t)g[(size_t)(x - s)]] = x;
        int m = 0;
        while (mark[(size_t)m] == x) ++m;
        g[(size_t)x] = m;
    }
    return g;
}
// 同じ石取りで x が巨大（10^18）なとき：Grundy 数の列はいつか周期的になるので、周期を見つけて使う
//   SubtractionGame G(S); G.grundy(x) / G.preperiod()（周期が始まる位置）/ G.period()
//   周期は limit までに見つからなければ止める。O(limit × |S|)
struct SubtractionGame {
    explicit SubtractionGame(const std::vector<long long>& S, int limit = 2000000) {
        if (S.empty()) {
            pre_ = 0, per_ = 1, g_ = {0};
            return;
        }
        long long smax = 0;
        for (long long s : S) smax = std::max(smax, s);
        if (smax > limit) math_internal::fail("SubtractionGame", "取れる数が limit より大きい");
        const int w = (int)smax;  // 直前 w 個の Grundy 数が同じなら、その先は全部同じ
        const unsigned long long B = 0x9E3779B97F4A7C15ULL;
        unsigned long long Bw = 1;
        for (int i = 0; i < w; ++i) Bw *= B;
        std::unordered_map<unsigned long long, int> first;  // 窓のハッシュ → 位置
        std::vector<int> mark(S.size() + 2, -1);
        unsigned long long h = 0;
        for (int x = 0; x <= limit; ++x) {
            if (x >= w) {  // 窓 g[x−w, x) を調べる
                auto itr = first.find(h);
                if (itr != first.end() && std::equal(g_.begin() + (itr->second - w), g_.begin() + itr->second, g_.begin() + (x - w))) {
                    per_ = x - itr->second, pre_ = itr->second - w;
                    while (pre_ > 0 && g_[(size_t)pre_ - 1] == g_[(size_t)(pre_ - 1 + per_)]) --pre_;  // 周期をもっと前から
                    return;
                }
                first.emplace(h, x);
            }
            for (long long s : S)
                if (s <= x && (size_t)g_[(size_t)(x - s)] <= S.size()) mark[(size_t)g_[(size_t)(x - s)]] = x;
            int m = 0;
            while (mark[(size_t)m] == x) ++m;
            g_.push_back(m);
            h = h * B + (unsigned long long)(m + 1);
            if (x >= w) h -= Bw * (unsigned long long)(g_[(size_t)(x - w)] + 1);
        }
        math_internal::fail("SubtractionGame", "limit までに周期が見つからない");
    }
    int grundy(long long x) const {
        if (x < 0) math_internal::fail("SubtractionGame::grundy", "x ≥ 0 で使う");
        if (x < (long long)g_.size()) return g_[(size_t)x];
        return g_[(size_t)(pre_ + (x - pre_) % per_)];
    }
    long long preperiod() const { return pre_; }
    long long period() const { return per_; }

  private:
    std::vector<int> g_;
    long long pre_ = 0, per_ = 1;
};

// ---- 1 手で「いくつかのゲームに分かれる」ときの Grundy 数 ----
// moves[v] = v の手の一覧。1 つの手は「行き先の状態の列」で、その手の値は行き先の Grundy 数の XOR
//   （山を 2 つに分ける：{a, b}、山を取り除くだけ：{}（値 0））。閉路があれば止める  O(n + 行き先の総数)
inline std::vector<int> grundy_dag_multi(int n, const std::vector<std::vector<std::vector<int>>>& moves) {
    if (n < 0 || (int)moves.size() != n) math_internal::fail("grundy_dag_multi", "moves の大きさは n");
    std::vector<int> g((size_t)n, -1), st;
    std::vector<char> state((size_t)n, 0);  // 0 = まだ、1 = 調べている途中、2 = 決まった
    std::vector<size_t> it_move((size_t)n, 0), it_pos((size_t)n, 0);
    std::vector<char> seen;
    for (int s = 0; s < n; ++s) {
        if (state[(size_t)s]) continue;
        st.push_back(s), state[(size_t)s] = 1;
        while (!st.empty()) {
            const int v = st.back();
            const auto& mv = moves[(size_t)v];
            bool pushed = false;
            while (it_move[(size_t)v] < mv.size()) {  // 次に調べる行き先
                const auto& to = mv[it_move[(size_t)v]];
                if (it_pos[(size_t)v] == to.size()) {
                    ++it_move[(size_t)v], it_pos[(size_t)v] = 0;
                    continue;
                }
                const int w = to[it_pos[(size_t)v]++];
                if (w < 0 || w >= n) math_internal::fail("grundy_dag_multi", "行き先の番号は 0〜n−1");
                if (state[(size_t)w] == 1) math_internal::fail("grundy_dag_multi", "閉路がある（同じ状態に戻れる）");
                if (state[(size_t)w] == 0) {
                    st.push_back(w), state[(size_t)w] = 1, pushed = true;
                    break;
                }
            }
            if (pushed) continue;
            seen.assign(mv.size() + 1, 0);
            for (const auto& to : mv) {
                int x = 0;
                for (int w : to) x ^= g[(size_t)w];
                if ((size_t)x <= mv.size()) seen[(size_t)x] = 1;
            }
            int m = 0;
            while (seen[(size_t)m]) ++m;
            g[(size_t)v] = m, state[(size_t)v] = 2, st.pop_back();
        }
    }
    return g;
}

// ---- 閉路があってもよいゲームの勝ち負け（後退解析）----
// 状態 0〜n−1、moves[v] = v から 1 手で行ける状態。動けない状態は手番の人の負け。
// result[v] = GAME_WIN（手番の人が勝つ）/ GAME_LOSE / GAME_DRAW（どちらも勝てず続く）、
// steps[v] = 決着までの手数（勝つ方は最短、負ける方は最長で。引き分けは −1）  O(n + 辺の数)
enum GameResultKind : int { GAME_LOSE = 0, GAME_WIN = 1, GAME_DRAW = 2 };
struct GameOutcome {
    std::vector<int> result, steps;
};
inline GameOutcome game_outcome(int n, const std::vector<std::vector<int>>& moves) {
    if (n < 0 || (int)moves.size() != n) math_internal::fail("game_outcome", "moves の大きさは n");
    std::vector<int> cnt((size_t)n), head((size_t)n + 1, 0), radj;  // 逆向きの辺（CSR）
    for (int v = 0; v < n; ++v) {
        cnt[(size_t)v] = (int)moves[(size_t)v].size();
        for (int w : moves[(size_t)v]) {
            if (w < 0 || w >= n) math_internal::fail("game_outcome", "行き先の番号は 0〜n−1");
            ++head[(size_t)w + 1];
        }
    }
    for (int v = 0; v < n; ++v) head[(size_t)v + 1] += head[(size_t)v];
    radj.resize((size_t)head[(size_t)n]);
    {
        std::vector<int> pos(head.begin(), head.end() - 1);
        for (int v = 0; v < n; ++v)
            for (int w : moves[(size_t)v]) radj[(size_t)pos[(size_t)w]++] = v;
    }
    GameOutcome r;
    r.result.assign((size_t)n, GAME_DRAW), r.steps.assign((size_t)n, -1);
    std::vector<int> q;
    q.reserve((size_t)n);
    for (int v = 0; v < n; ++v)
        if (cnt[(size_t)v] == 0) r.result[(size_t)v] = GAME_LOSE, r.steps[(size_t)v] = 0, q.push_back(v);
    for (size_t h = 0; h < q.size(); ++h) {  // 手数の小さい順に決まる
        const int v = q[h];
        for (int i = head[(size_t)v]; i < head[(size_t)v + 1]; ++i) {
            const int u = radj[(size_t)i];
            if (r.result[(size_t)u] != GAME_DRAW) continue;
            if (r.result[(size_t)v] == GAME_LOSE) {  // 負けの状態に行ければ勝ち
                r.result[(size_t)u] = GAME_WIN, r.steps[(size_t)u] = r.steps[(size_t)v] + 1, q.push_back(u);
            } else if (--cnt[(size_t)u] == 0) {  // どこへ行っても相手の勝ち
                r.result[(size_t)u] = GAME_LOSE, r.steps[(size_t)u] = r.steps[(size_t)v] + 1, q.push_back(u);
            }
        }
    }
    return r;
}
#endif  // GG_GAME
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace game_test {

#define GAME_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[game] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 山の組を状態にして、全部の手を試す（勝ち負けだけ）
static bool brute_win(std::vector<int> piles, const std::vector<int>& S, bool misere, std::map<std::vector<int>, bool>& memo) {
    std::sort(piles.begin(), piles.end());
    auto it = memo.find(piles);
    if (it != memo.end()) return it->second;
    bool any_move = false, win = false;
    for (size_t i = 0; i < piles.size() && !win; ++i)
        for (int s : S) {
            if (s > piles[i]) continue;
            any_move = true;
            std::vector<int> q = piles;
            q[i] -= s;
            if (!brute_win(q, S, misere, memo)) { win = true; break; }
        }
    if (!any_move) win = misere;  // 動けない：普通は負け、逆形は（最後に取ったのが相手なので）勝ち
    return memo[piles] = win;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 10);
        std::vector<long long> v((size_t)n);
        for (auto& x : v) x = rnd(rng, -2, 12);
        int e = 0;
        while (std::find(v.begin(), v.end(), (long long)e) != v.end()) ++e;
        GAME_CHECK(mex(v) == e);
        std::vector<unsigned> u;
        for (long long x : v)
            if (x >= 0) u.push_back((unsigned)x);
        GAME_CHECK(mex(u) == e);
    }
    // Nim・逆形 Nim：全部の手を試す
    std::map<std::vector<int>, bool> m1, m2;
    const std::vector<int> all{1, 2, 3, 4, 5, 6, 7};
    for (int it = 0; it < 3000; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 4);
        std::vector<int> p((size_t)n);
        for (auto& x : p) x = (int)rnd(rng, 0, 7);
        GAME_CHECK(nim_first_wins(p) == brute_win(p, all, false, m1));
        GAME_CHECK(misere_nim_first_wins(p) == brute_win(p, all, true, m2));
    }
}

static void test_dag(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 40);
        std::vector<int> perm((size_t)n);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        std::vector<std::vector<int>> mv((size_t)n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j)
                if (rng() % 4 == 0) mv[(size_t)perm[(size_t)i]].push_back(perm[(size_t)j]);  // 番号の大きい方から小さい方へ（並べ替えた番号で）
        for (auto& v : mv)
            if (!v.empty() && rng() % 3 == 0) v.push_back(v[0]);  // 同じ行き先が 2 回
        const std::vector<int> g = grundy_dag(n, mv);
        std::vector<int> e((size_t)n, -1);  // 小さい方から順に
        for (int i = 0; i < n; ++i) {
            const int v = perm[(size_t)i];
            std::vector<int> ws;
            for (int w : mv[(size_t)v]) ws.push_back(e[(size_t)w]);
            e[(size_t)v] = mex(ws);
        }
        GAME_CHECK(g == e);
    }
    // 長い 1 本道（再帰しないので深くても大丈夫）
    const int N = 1000000;
    std::vector<std::vector<int>> mv((size_t)N);
    for (int i = 1; i < N; ++i) mv[(size_t)i].push_back(i - 1);
    const std::vector<int> g = grundy_dag(N, mv);
    GAME_CHECK(g[0] == 0 && g[1] == 1 && g[N - 1] == (N - 1) % 2);
    ++cases;
}

static void test_multi_and_outcome(std::mt19937_64& rng, int& cases) {
    // 1 手で複数のゲームに分かれる
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 30);
        std::vector<std::vector<std::vector<int>>> mv((size_t)n);
        std::vector<std::vector<int>> single((size_t)n);
        bool all_single = it % 3 == 0;
        for (int v = 1; v < n; ++v) {  // 小さい番号へだけ
            const int k = (int)rnd(rng, 0, 4);
            for (int t = 0; t < k; ++t) {
                std::vector<int> to;
                const int len = all_single ? 1 : (int)rnd(rng, 0, 3);
                for (int j = 0; j < len; ++j) to.push_back((int)rnd(rng, 0, v - 1));
                mv[(size_t)v].push_back(to);
                if (all_single) single[(size_t)v].push_back(to[0]);
            }
        }
        std::vector<int> e((size_t)n);
        for (int v = 0; v < n; ++v) {
            std::vector<int> vals;
            for (const auto& to : mv[(size_t)v]) {
                int x = 0;
                for (int w : to) x ^= e[(size_t)w];
                vals.push_back(x);
            }
            e[(size_t)v] = mex(vals);
        }
        GAME_CHECK(grundy_dag_multi(n, mv) == e && (!all_single || grundy_dag(n, single) == e));
    }
    // 後退解析：手数ごとに決める愚直解と（閉路あり）
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 25);
        std::vector<std::vector<int>> mv((size_t)n);
        for (int v = 0; v < n; ++v)
            for (int w = 0; w < n; ++w)
                if (rng() % (it % 2 ? 4 : 8) == 0) mv[(size_t)v].push_back(w);
        std::vector<int> res((size_t)n, GAME_DRAW), st((size_t)n, -1);
        for (int v = 0; v < n; ++v)
            if (mv[(size_t)v].empty()) res[(size_t)v] = GAME_LOSE, st[(size_t)v] = 0;
        for (int k = 1; k <= 2 * n + 2; ++k) {  // 手数 k で決まるもの
            std::vector<int> nr = res, ns = st;
            for (int u = 0; u < n; ++u) {
                if (res[(size_t)u] != GAME_DRAW) continue;
                bool win = false, all_win = true;
                for (int w : mv[(size_t)u]) {
                    win |= res[(size_t)w] == GAME_LOSE && st[(size_t)w] == k - 1;
                    all_win &= res[(size_t)w] == GAME_WIN;
                }
                if (win) nr[(size_t)u] = GAME_WIN, ns[(size_t)u] = k;
                else if (all_win) nr[(size_t)u] = GAME_LOSE, ns[(size_t)u] = k;
            }
            res = nr, st = ns;
        }
        const GameOutcome o = game_outcome(n, mv);
        GAME_CHECK(o.result == res && o.steps == st);
    }
}

static void test_take(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 400 * TEST_SCALE; ++it, ++cases) {
        const long long L = rnd(rng, 1, 20), R = rnd(rng, L, 40);
        std::vector<long long> S;
        for (long long s = L; s <= R; ++s) S.push_back(s);
        const std::vector<int> g = subtraction_game_grundy(S, 300);
        for (int x = 0; x <= 300; ++x) GAME_CHECK(range_take_grundy(x, L, R) == g[(size_t)x]);
        const long long X = rnd(rng, 0, 1000000000000000000LL);
        const SubtractionGame G(S);
        GAME_CHECK(G.grundy(X) == range_take_grundy(X, L, R) && G.period() == L + R);
    }
    GAME_CHECK(range_take_grundy(LLONG_MAX, LLONG_MAX, LLONG_MAX) == 1 && range_take_grundy(5, 1, 1) == 1);
    // 集合 S：メモ化探索の Grundy 数と、周期
    for (int it = 0; it < 400 * TEST_SCALE; ++it, ++cases) {
        const int k = (int)rnd(rng, 1, 5);
        std::vector<long long> S;
        for (int i = 0; i < k; ++i) S.push_back(rnd(rng, 1, it % 2 ? 10 : 40));
        const int N = 3000;
        const std::vector<int> g = subtraction_game_grundy(S, N);
        for (int x = 0; x <= N; ++x) {  // 定義どおり
            std::vector<int> ws;
            for (long long s : S)
                if (s <= x) ws.push_back(g[(size_t)(x - s)]);
            GAME_CHECK(g[(size_t)x] == mex(ws));
        }
        const SubtractionGame G(S);
        for (int x = 0; x <= N; ++x) GAME_CHECK(G.grundy(x) == g[(size_t)x]);
        const long long P = G.preperiod(), T = G.period();
        GAME_CHECK(T >= 1 && P >= 0);
        if (P + 2 * T <= N) {
            for (long long x = P; x + T <= N; ++x) GAME_CHECK(g[(size_t)x] == g[(size_t)(x + T)]);
            GAME_CHECK(P == 0 || g[(size_t)P - 1] != g[(size_t)(P - 1 + T)]);  // 周期はできるだけ前から
        }
        const long long X = rnd(rng, P + T, 1000000000000000000LL);
        GAME_CHECK(G.grundy(X) == G.grundy(X - T));
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_basic(rng, cases);
    test_dag(rng, cases);
    test_multi_and_outcome(rng, cases);
    test_take(rng, cases);
    std::cout << "[game] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[game] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 200000;  // 約数へ動けるゲーム（ABC368 F の形）
    std::vector<std::vector<int>> mv((size_t)N + 1);
    for (int d = 1; d <= N; ++d)
        for (int m = 2 * d; m <= N; m += d) mv[(size_t)m].push_back(d);
    auto t0 = std::chrono::steady_clock::now();
    const std::vector<int> g = grundy_dag(N + 1, mv);
    const long long t_dag = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<int> h = subtraction_game_grundy({1, 3, 7, 19, 31, 47, 53, 61, 97, 101}, 3000000);
    const long long t_sub = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const SubtractionGame G({2, 5, 11, 23, 37, 89, 123});
    const long long t_per = ms(t0);
    std::cout << "[game] benchmark: grundy_dag（2×10^5 状態, 約 2.2×10^6 辺）" << t_dag << " ms, subtraction_game_grundy（3×10^6, |S| = 10）" << t_sub
              << " ms, SubtractionGame の周期 " << t_per << " ms（周期 " << G.period() << "）  (check " << (g[N] + h.back() + G.grundy(1000000000000LL)) % 7 << ")\n";
}

}  // namespace game_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace game_demo {

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
#define SHOW(expr, comment) game_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    game_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("ゲーム（game） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  mex / grundy_dag / nim_first_wins / misere_nim_first_wins / range_take_grundy / subtraction_game_grundy / SubtractionGame\n"
              << "  先手必勝 ⇔ 並べたゲームの Grundy 数の XOR ≠ 0\n";

    api("【1】 mex(vs)", "int", "vs に含まれない最小の 0 以上の整数（Grundy 数 = 行き先の Grundy 数の mex）", "O(|vs|)");
    SHOW(mex(std::vector<int>{0, 1, 3}), "");
    SHOW(mex(std::vector<int>{1, 2}), "");
    SHOW(mex(std::vector<int>{}), "行き先がない（負けの状態）→ 0");
    {
        // ABC206 F（Interval Game 2）：g[l][r] = mex{ g[l][L_i] ^ g[R_i][r] : [L_i, R_i) ⊆ [l, r) }
        auto solve = [](const std::vector<std::pair<int, int>>& iv) {
            std::vector<std::vector<int>> g(101, std::vector<int>(101, 0));
            for (int len = 1; len <= 99; ++len)
                for (int l = 1; l + len <= 100; ++l) {
                    const int r = l + len;
                    std::vector<int> ws;
                    for (const auto& p : iv)
                        if (l <= p.first && p.second <= r) ws.push_back(g[(size_t)l][(size_t)p.first] ^ g[(size_t)p.second][(size_t)r]);
                    g[(size_t)l][(size_t)r] = mex(ws);
                }
            return std::string(g[1][100] ? "Alice" : "Bob");
        };
        note("ABC206 F（Interval Game 2）: 区間 [l, r) の Grundy 数 = mex（選んだ区間の左右の Grundy 数の XOR）");
        std::cout << "    例 → " << solve({{53, 98}, {8, 43}, {12, 53}}) << ", "
                  << solve({{4, 7}, {5, 7}, {3, 7}, {4, 5}, {5, 8}, {6, 9}, {4, 8}, {5, 10}, {1, 9}, {5, 10}}) << ", " << solve({{58, 98}, {11, 29}}) << ", "
                  << solve({{79, 83}, {44, 83}, {38, 74}, {49, 88}, {18, 45}, {64, 99}}) << ", " << solve({{5, 9}}) << "\n";
    }

    api("【2】 grundy_dag(n, moves)", "vector<int>", "状態 0〜n−1、moves[v] = v から 1 手で行ける状態。各状態の Grundy 数（閉路があれば止める。再帰しない）", "O(n + 辺の数)");
    RUN(std::vector<std::vector<int>> mv = {{}, {0}, {0, 1}, {1, 2}});
    SHOW(grundy_dag(4, mv), "0 → 行き先なし、1 → {0}、2 → {0, 1}、3 → {1, 2}");
    {
        // ABC368 F（Dividing Game）：A_i を真の約数にする。状態 = 値、行き先 = 真の約数
        auto solve = [](const std::vector<int>& A) {
            const int M = *std::max_element(A.begin(), A.end());
            std::vector<std::vector<int>> moves((size_t)M + 1);
            for (int d = 1; d <= M; ++d)
                for (int m = 2 * d; m <= M; m += d) moves[(size_t)m].push_back(d);
            const std::vector<int> g = grundy_dag(M + 1, moves);
            int x = 0;
            for (int a : A) x ^= g[(size_t)a];
            return std::string(x ? "Anna" : "Bruno");
        };
        note("ABC368 F（Dividing Game）: 状態 = 値、行き先 = 真の約数 で grundy_dag、XOR");
        std::cout << "    例 1 → " << solve({2, 3, 4}) << "、例 2 → " << solve({2, 3, 4, 6}) << "\n";
    }

    api("【3】 nim_first_wins(piles) / misere_nim_first_wins(piles)", "bool", "Nim（最後に取った方が勝ち）/ 逆形 Nim（最後に取った方が負け）で先手必勝か", "O(N)");
    SHOW(nim_first_wins(std::vector<int>{1, 2, 3}), "1 ^ 2 ^ 3 = 0 → 後手必勝");
    SHOW(nim_first_wins(std::vector<int>{3, 4, 5}), "");
    SHOW(misere_nim_first_wins(std::vector<int>{1, 1}), "全部 1：1 の山が偶数個なら先手必勝");
    SHOW(misere_nim_first_wins(std::vector<int>{1, 2, 3}), "2 以上がある：普通の Nim と同じ");

    api("【4】 range_take_grundy(x, L, R)", "long long", "1 回に L〜R 個取れる石取りの x 個の山の Grundy 数 ⌊(x mod (L + R)) / L⌋", "O(1)");
    SHOW(range_take_grundy(10, 1, 2), "10 mod 3 = 1");
    SHOW(range_take_grundy(7, 3, 5), "7 mod 8 = 7 → ⌊7/3⌋");
    {
        // ABC297 G（Constrained Nim 2）
        auto solve = [](long long L, long long R, const std::vector<long long>& A) {
            long long x = 0;
            for (long long a : A) x ^= range_take_grundy(a, L, R);
            return std::string(x ? "First" : "Second");
        };
        note("ABC297 G（Constrained Nim 2）: 各山の range_take_grundy の XOR");
        std::cout << "    例 1 → " << solve(1, 2, {2, 3, 3}) << "、例 2 → " << solve(1, 1, {3, 1, 4, 1, 5}) << "、例 3 → " << solve(3, 14, {10, 20, 30, 40, 50, 60, 70}) << "\n";
    }

    api("【5】 subtraction_game_grundy(S, N) / SubtractionGame G(S)", "vector<int> / 構造体",
        "1 回に取れる数が集合 S の石取りの Grundy 数 g(0)〜g(N) / x が巨大なとき：周期を見つけて G.grundy(x)", "O(N |S|) / O(周期が見つかるまで × |S|)");
    SHOW(subtraction_game_grundy({1, 3, 4}, 14), "取れる数 {1, 3, 4}");
    RUN(SubtractionGame G({1, 3, 4}));
    SHOW(G.period(), "周期 7（0 1 0 1 2 3 2 のくり返し）");
    SHOW(G.preperiod(), "最初から周期的");
    SHOW(G.grundy(1000000000000000000LL), "10^18 mod 7 = 1 → g(1)");
    SHOW(subtraction_game_grundy({2, 4, 7}, 16), "取れる数 {2, 4, 7}");
    RUN(SubtractionGame H({2, 4, 7}));
    SHOW(H.preperiod(), "g(8) から周期的になる（途中から周期的な例）");
    SHOW(H.period(), "1 0 2 のくり返し");

    api("【6】 grundy_dag_multi(n, moves)", "vector<int>",
        "1 手で「いくつかのゲームに分かれる」とき：moves[v] は手の一覧、1 つの手は行き先の状態の列（その手の値 = 行き先の Grundy 数の XOR）", "O(n + 行き先の総数)");
    {
        // グランディのゲーム：山を大きさの違う 2 つの山に分ける
        const int N = 16;
        std::vector<std::vector<std::vector<int>>> moves2(N);
        for (int x = 0; x < N; ++x)
            for (int a = 1; 2 * a < x; ++a) moves2[(size_t)x].push_back({a, x - a});
        const std::vector<int> g = grundy_dag_multi(N, moves2);
        std::cout << "  グランディのゲーム（山を大きさの違う 2 つに分ける）g(0)〜g(15) -> ";
        for (int x : g) std::cout << x << " ";
        std::cout << "\n";
    }

    api("【7】 game_outcome(n, moves)", "GameOutcome{result, steps}",
        "閉路があってもよいゲームの後退解析。result[v] = GAME_WIN / GAME_LOSE / GAME_DRAW（v で手番の人から見て）、steps[v] = 決着までの手数（勝つ方は最短、負ける方は最長）",
        "O(n + 辺の数)");
    {
        // ABC209 E（Shiritori）：単語 = 状態、末尾 3 文字 = 次の単語の先頭 3 文字なら行ける。高橋くんが単語 i を言うと青木くんの手番
        const std::vector<std::string> w = {"abcd", "bcda", "ada"};
        const int n = (int)w.size();
        std::vector<std::vector<int>> moves2(n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (w[(size_t)i].substr(w[(size_t)i].size() - 3) == w[(size_t)j].substr(0, 3)) moves2[(size_t)i].push_back(j);
        const GameOutcome o = game_outcome(n, moves2);
        note("ABC209 E（Shiritori）: 単語ごとの状態で後退解析（実際は「先頭 3 文字」の状態にまとめて辺を減らす）");
        std::cout << "    ";
        for (int i = 0; i < n; ++i) std::cout << (o.result[(size_t)i] == GAME_WIN ? "Aoki" : o.result[(size_t)i] == GAME_LOSE ? "Takahashi" : "Draw") << " ";
        std::cout << "（steps: " << o.steps[0] << " " << o.steps[1] << " " << o.steps[2] << "）\n";
    }
}

}  // namespace game_demo

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
    game_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    game_test::benchmark();    // 速度計測
    game_demo::demo();         // 使い方デモ
    return 0;
}
