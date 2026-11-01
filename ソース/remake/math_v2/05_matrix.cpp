// =============================================================================
//  行列（Matrix<T>：mint・long long・double・半環（MinPlus・MinMax など）、行列の累乗、漸化式の行列、MatrixFixed<T, K>）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 行列・線形代数
//  依存       : 00_int_common, 00_math_base, 00_modint（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  使い方
//      using mint = ModInt107;
//      Matrix<mint> A = {{1, 1}, {1, 0}};     // 2×2（Matrix<mint> A(n, m) なら全部 0）
//      cout << A.pow(K)(0, 1) << '\n';         // フィボナッチ数 F_K
//  T が long long のときは溢れたら止める（数え上げは mint で）
//
//  機能一覧                                                                    計算量
//    Matrix<T> A(n, m) / A = {{…}, {…}} / Matrix<T>::identity(n)   0 の行列 / 値を並べて / 単位行列
//      A(i, j) / A[i][j] / A.rows() / A.cols()                    成分 / 大きさ
//      A + B / A - B / A * B / A * s（スカラー）/ += -= *=         O(n m k)
//      A.pow(e)                       A^e（e ≤ 10^18）                                          O(n^3 log e)
//      A.apply(v)                     A v（行列 × 縦ベクトル。vector<T>）                         O(n m)
//      A.transpose() / A.trace()      転置 / 対角成分の和
//      cout << A                      1 行ずつ空白区切りで
//    半環（行列の + と × を別のものにする。Matrix<MinPlus<long long>> のように使う）
//      MinPlus<T>                     + が min、× が +：ちょうど K 歩の最短路（道がないのは INF）
//      MaxPlus<T>                     + が max、× が +：ちょうど K 歩の最長路
//      MinMax<T> / MaxMin<T>          + が min、× が max：通る辺の最大を最小に（ボトルネック）/ その逆（最大容量路）
//      XorAnd                         + が XOR、× が AND（unsigned long long。ABC009 D）
//      OrAnd                          + が OR、× が AND：ちょうど K 歩で行けるか
//      自作の半環：operator+ と operator* と static zero()（+ の単位元・× で 0 になる）/ static one()（× の単位元）を
//                 持つ struct を作れば Matrix<自作> で使える（== と << もあると便利）
//    companion_matrix(c)              a_i = c_0 a_{i−1} + … + c_{k−1} a_{i−k} の行列
//    recurrence_nth_matrix(c, a, n)   その漸化式の a_n（a = a_0〜a_{k−1}。半環でもよい。long long は        O(k^3 log n)
//                                     途中の行列の値が溢れると、a_n が収まっても止まる → mint で）
//    MatrixFixed<T, K>                K×K の小さい行列（2×2・3×3 などを速く）。M(i, j) / M[i][j]、+ − × スカラー倍、
//                                     pow / apply（std::array）/ transpose / trace と
//      M.det() / M.inv() / M.try_inv()    行列式 / 逆行列（T は mint か double）                       O(K^3)
//  （N×N の行列式・逆行列・連立一次方程式・ランクは 05_linear_algebra）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mat_test::random_test()  … 愚直な掛け算・DP・余因子展開との突き合わせ
//    mat_test::benchmark()    … 速度計測
//    mat_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
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
#ifndef GG_MATRIX
#define GG_MATRIX
// ---- 行列（数学第二世代）  Matrix<T>（mint・long long・double・半環）・MatrixFixed<T, K>（小さい行列を速く）・漸化式の行列 ----
namespace math_internal {
template <class T, class = void>
struct has_semiring_one : std::false_type {};
template <class T>
struct has_semiring_one<T, std::void_t<decltype(T::one())>> : std::true_type {};
template <class T>
T sr_zero() {  // 足し算の単位元（普通の数は 0、半環は T::zero()）
    if constexpr (has_semiring_one<T>::value) return T::zero();
    else return T(0);
}
template <class T>
T sr_one() {  // 掛け算の単位元（普通の数は 1、半環は T::one()）
    if constexpr (has_semiring_one<T>::value) return T::one();
    else return T(1);
}
template <class T>
inline void mat_fma(T& acc, const T& a, const T& b) {  // acc += a * b（整数の型は溢れたら止める）
    if constexpr (std::is_integral_v<T>) {
        T p;
        if (__builtin_mul_overflow(a, b, &p) || __builtin_add_overflow(acc, p, &acc)) fail("Matrix", "値が整数の型に収まらない（数え上げは mint で）");
    } else {
        acc = acc + a * b;
    }
}
template <class T>
inline T mat_sub(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T s;
        if (__builtin_sub_overflow(a, b, &s)) fail("Matrix", "値が整数の型に収まらない（数え上げは mint で）");
        return s;
    } else {
        return a - b;
    }
}
template <class T>
inline T mat_mul(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T s;
        if (__builtin_mul_overflow(a, b, &s)) fail("Matrix", "値が整数の型に収まらない（数え上げは mint で）");
        return s;
    } else {
        return a * b;
    }
}
template <class T>
inline T mat_add(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T s;
        if (__builtin_add_overflow(a, b, &s)) fail("Matrix", "値が整数の型に収まらない（数え上げは mint で）");
        return s;
    } else {
        return a + b;
    }
}
}  // namespace math_internal

// ---- 半環（行列の足し算・掛け算を別のものに変える） ----
template <class T>
struct MinPlus {  // 足し算 = min、掛け算 = +（ちょうど K 歩の最短路）。道がないのは INF
    static constexpr T INF = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();
    T v;
    MinPlus() : v(INF) {}
    MinPlus(T x) : v(x) {}
    static MinPlus zero() { return MinPlus(); }
    static MinPlus one() { return MinPlus(T(0)); }
    bool is_inf() const { return v == INF; }
    friend MinPlus operator+(const MinPlus& a, const MinPlus& b) { return MinPlus(a.v < b.v ? a.v : b.v); }
    friend MinPlus operator*(const MinPlus& a, const MinPlus& b) {
        if (a.v == INF || b.v == INF) return MinPlus();
        if constexpr (std::is_integral_v<T>) {
            T s;
            if (__builtin_add_overflow(a.v, b.v, &s) || s == INF) math_internal::fail("MinPlus", "長さが型に収まらない");
            return MinPlus(s);
        } else {
            return MinPlus(a.v + b.v);
        }
    }
    friend bool operator==(const MinPlus& a, const MinPlus& b) { return a.v == b.v; }
    friend bool operator!=(const MinPlus& a, const MinPlus& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const MinPlus& a) { return a.v == INF ? os << "INF" : os << a.v; }
};
template <class T>
struct MaxPlus {  // 足し算 = max、掛け算 = +（ちょうど K 歩の最長路）。道がないのは −INF
    static constexpr T NEG_INF = std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest();
    T v;
    MaxPlus() : v(NEG_INF) {}
    MaxPlus(T x) : v(x) {}
    static MaxPlus zero() { return MaxPlus(); }
    static MaxPlus one() { return MaxPlus(T(0)); }
    bool is_inf() const { return v == NEG_INF; }
    friend MaxPlus operator+(const MaxPlus& a, const MaxPlus& b) { return MaxPlus(a.v > b.v ? a.v : b.v); }
    friend MaxPlus operator*(const MaxPlus& a, const MaxPlus& b) {
        if (a.v == NEG_INF || b.v == NEG_INF) return MaxPlus();
        if constexpr (std::is_integral_v<T>) {
            T s;
            if (__builtin_add_overflow(a.v, b.v, &s) || s == NEG_INF) math_internal::fail("MaxPlus", "長さが型に収まらない");
            return MaxPlus(s);
        } else {
            return MaxPlus(a.v + b.v);
        }
    }
    friend bool operator==(const MaxPlus& a, const MaxPlus& b) { return a.v == b.v; }
    friend bool operator!=(const MaxPlus& a, const MaxPlus& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const MaxPlus& a) { return a.v == NEG_INF ? os << "-INF" : os << a.v; }
};
template <class T>
struct MinMax {  // 足し算 = min、掛け算 = max（ちょうど K 歩で、通る辺の重みの最大を最小に。ボトルネック）。道がないのは INF
    static constexpr T INF = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();
    static constexpr T NEG_INF = std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest();
    T v;
    MinMax() : v(INF) {}
    MinMax(T x) : v(x) {}
    static MinMax zero() { return MinMax(); }         // 道がない
    static MinMax one() { return MinMax(NEG_INF); }   // 0 歩（辺を通らない）
    bool is_inf() const { return v == INF; }
    friend MinMax operator+(const MinMax& a, const MinMax& b) { return MinMax(a.v < b.v ? a.v : b.v); }
    friend MinMax operator*(const MinMax& a, const MinMax& b) { return MinMax(a.v > b.v ? a.v : b.v); }
    friend bool operator==(const MinMax& a, const MinMax& b) { return a.v == b.v; }
    friend bool operator!=(const MinMax& a, const MinMax& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const MinMax& a) { return a.v == INF ? os << "INF" : (a.v == NEG_INF ? os << "-INF" : os << a.v); }
};
template <class T>
struct MaxMin {  // 足し算 = max、掛け算 = min（ちょうど K 歩で、通る辺の重みの最小を最大に。最大容量路）。道がないのは −INF
    static constexpr T INF = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();
    static constexpr T NEG_INF = std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest();
    T v;
    MaxMin() : v(NEG_INF) {}
    MaxMin(T x) : v(x) {}
    static MaxMin zero() { return MaxMin(); }      // 道がない
    static MaxMin one() { return MaxMin(INF); }    // 0 歩
    bool is_inf() const { return v == NEG_INF; }
    friend MaxMin operator+(const MaxMin& a, const MaxMin& b) { return MaxMin(a.v > b.v ? a.v : b.v); }
    friend MaxMin operator*(const MaxMin& a, const MaxMin& b) { return MaxMin(a.v < b.v ? a.v : b.v); }
    friend bool operator==(const MaxMin& a, const MaxMin& b) { return a.v == b.v; }
    friend bool operator!=(const MaxMin& a, const MaxMin& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const MaxMin& a) { return a.v == NEG_INF ? os << "-INF" : (a.v == INF ? os << "INF" : os << a.v); }
};
struct XorAnd {  // 足し算 = XOR、掛け算 = AND（ビットごとの漸化式。ABC009 D）
    unsigned long long v;
    XorAnd(unsigned long long x = 0) : v(x) {}
    static XorAnd zero() { return XorAnd(0); }
    static XorAnd one() { return XorAnd(~0ULL); }
    friend XorAnd operator+(const XorAnd& a, const XorAnd& b) { return XorAnd(a.v ^ b.v); }
    friend XorAnd operator*(const XorAnd& a, const XorAnd& b) { return XorAnd(a.v & b.v); }
    friend bool operator==(const XorAnd& a, const XorAnd& b) { return a.v == b.v; }
    friend bool operator!=(const XorAnd& a, const XorAnd& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const XorAnd& a) { return os << a.v; }
};
struct OrAnd {  // 足し算 = OR、掛け算 = AND（ちょうど K 歩で行けるか）
    bool v;
    OrAnd(bool x = false) : v(x) {}
    static OrAnd zero() { return OrAnd(false); }
    static OrAnd one() { return OrAnd(true); }
    friend OrAnd operator+(const OrAnd& a, const OrAnd& b) { return OrAnd(a.v || b.v); }
    friend OrAnd operator*(const OrAnd& a, const OrAnd& b) { return OrAnd(a.v && b.v); }
    friend bool operator==(const OrAnd& a, const OrAnd& b) { return a.v == b.v; }
    friend bool operator!=(const OrAnd& a, const OrAnd& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const OrAnd& a) { return os << (a.v ? 1 : 0); }
};

// ---- 行列 Matrix<T>（n 行 m 列。T は mint・long long・double・上の半環） ----
template <class T>
struct Matrix {
    Matrix() : n_(0), m_(0) {}
    Matrix(int n, int m) : n_(n), m_(m) {  // すべて 0（半環なら zero）
        if (n < 0 || m < 0) math_internal::fail("Matrix", "大きさは 0 以上");
        a_.assign((size_t)n * (size_t)m, math_internal::sr_zero<T>());
    }
    explicit Matrix(int n) : Matrix(n, n) {}
    Matrix(const std::vector<std::vector<T>>& v) : Matrix((int)v.size(), v.empty() ? 0 : (int)v[0].size()) {
        for (int i = 0; i < n_; ++i) {
            if ((int)v[(size_t)i].size() != m_) math_internal::fail("Matrix", "行の長さがそろっていない");
            for (int j = 0; j < m_; ++j) (*this)(i, j) = v[(size_t)i][(size_t)j];
        }
    }
    Matrix(std::initializer_list<std::initializer_list<T>> v) : Matrix(std::vector<std::vector<T>>(v.begin(), v.end())) {}
    static Matrix identity(int n) {  // 単位行列
        Matrix r(n, n);
        for (int i = 0; i < n; ++i) r(i, i) = math_internal::sr_one<T>();
        return r;
    }
    static Matrix zero(int n, int m) { return Matrix(n, m); }
    int rows() const { return n_; }
    int cols() const { return m_; }
    T& operator()(int i, int j) { return a_[(size_t)i * (size_t)m_ + (size_t)j]; }
    const T& operator()(int i, int j) const { return a_[(size_t)i * (size_t)m_ + (size_t)j]; }
    T* operator[](int i) { return a_.data() + (size_t)i * (size_t)m_; }  // A[i][j] とも書ける
    const T* operator[](int i) const { return a_.data() + (size_t)i * (size_t)m_; }

    Matrix& operator+=(const Matrix& b) {
        same_size(b, "Matrix::+");
        for (size_t i = 0; i < a_.size(); ++i) a_[i] = math_internal::mat_add(a_[i], b.a_[i]);
        return *this;
    }
    Matrix& operator-=(const Matrix& b) {  // 半環では使えない
        same_size(b, "Matrix::-");
        for (size_t i = 0; i < a_.size(); ++i) a_[i] = math_internal::mat_sub(a_[i], b.a_[i]);
        return *this;
    }
    Matrix& operator*=(const Matrix& b) { return *this = *this * b; }
    Matrix& operator*=(const T& s) {  // スカラー倍
        for (T& x : a_) x = math_internal::mat_mul(x, s);
        return *this;
    }
    friend Matrix operator+(Matrix a, const Matrix& b) { return a += b; }
    friend Matrix operator-(Matrix a, const Matrix& b) { return a -= b; }
    friend Matrix operator*(Matrix a, const T& s) { return a *= s; }
    friend Matrix operator*(const T& s, Matrix a) { return a *= s; }
    friend Matrix operator*(const Matrix& a, const Matrix& b) {  // O(n m k)
        if (a.m_ != b.n_) math_internal::fail("Matrix::*", "左の列数と右の行数が違う");
        const T z = math_internal::sr_zero<T>();
        Matrix c(a.n_, b.m_);
        for (int i = 0; i < a.n_; ++i) {
            T* ci = c[i];
            for (int k = 0; k < a.m_; ++k) {
                const T aik = a(i, k);
                if (aik == z) continue;
                const T* bk = b[k];
                for (int j = 0; j < b.m_; ++j) math_internal::mat_fma(ci[j], aik, bk[j]);
            }
        }
        return c;
    }
    Matrix pow(long long e) const {  // A^e（e ≥ 0。e ≤ 10^18 でよい）  O(n^3 log e)
        if (n_ != m_) math_internal::fail("Matrix::pow", "正方行列で使う");
        if (e < 0) math_internal::fail("Matrix::pow", "e ≥ 0 で使う");
        Matrix r = identity(n_), b = *this;
        for (; e > 0; e >>= 1) {
            if (e & 1) r = r * b;
            if (e > 1) b = b * b;
        }
        return r;
    }
    std::vector<T> apply(const std::vector<T>& v) const {  // A v（行列 × 縦ベクトル）  O(n m)
        if ((int)v.size() != m_) math_internal::fail("Matrix::apply", "ベクトルの長さが列数と違う");
        std::vector<T> r((size_t)n_, math_internal::sr_zero<T>());
        for (int i = 0; i < n_; ++i)
            for (int j = 0; j < m_; ++j) math_internal::mat_fma(r[(size_t)i], (*this)(i, j), v[(size_t)j]);
        return r;
    }
    Matrix transpose() const {
        Matrix r(m_, n_);
        for (int i = 0; i < n_; ++i)
            for (int j = 0; j < m_; ++j) r(j, i) = (*this)(i, j);
        return r;
    }
    T trace() const {  // 対角成分の和
        if (n_ != m_) math_internal::fail("Matrix::trace", "正方行列で使う");
        T s = math_internal::sr_zero<T>();
        for (int i = 0; i < n_; ++i) s = math_internal::mat_add(s, (*this)(i, i));
        return s;
    }
    friend bool operator==(const Matrix& a, const Matrix& b) { return a.n_ == b.n_ && a.m_ == b.m_ && a.a_ == b.a_; }
    friend bool operator!=(const Matrix& a, const Matrix& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const Matrix& a) {  // 1 行ずつ（行の中は空白区切り）
        for (int i = 0; i < a.n_; ++i) {
            for (int j = 0; j < a.m_; ++j) os << (j ? " " : "") << a(i, j);
            if (i + 1 < a.n_) os << '\n';
        }
        return os;
    }

  private:
    int n_, m_;
    std::vector<T> a_;
    void same_size(const Matrix& b, const char* f) const {
        if (n_ != b.n_ || m_ != b.m_) math_internal::fail(f, "大きさが違う");
    }
};

// ---- 線形漸化式の行列 ----
// a_i = c_0 a_{i−1} + c_1 a_{i−2} + … + c_{k−1} a_{i−k} の行列（1 行目が c、その下が 1 つずらす形）
//   (a_{i}, a_{i−1}, …, a_{i−k+1})^T = M (a_{i−1}, …, a_{i−k})^T
template <class T>
Matrix<T> companion_matrix(const std::vector<T>& c) {
    const int k = (int)c.size();
    if (k == 0) math_internal::fail("companion_matrix", "c が空");
    Matrix<T> m(k, k);
    for (int j = 0; j < k; ++j) m(0, j) = c[(size_t)j];
    for (int i = 1; i < k; ++i) m(i, i - 1) = math_internal::sr_one<T>();
    return m;
}
// 上の漸化式で a_0〜a_{k−1} が a のとき a_n（n は 0 から数える。n ≤ 10^18。半環でもよい）  O(k^3 log n)
template <class T>
T recurrence_nth_matrix(const std::vector<T>& c, const std::vector<T>& a, long long n) {
    const int k = (int)c.size();
    if ((int)a.size() != k) math_internal::fail("recurrence_nth_matrix", "c と a の長さが違う");
    if (n < 0) math_internal::fail("recurrence_nth_matrix", "n ≥ 0 で使う");
    if (n < k) return a[(size_t)n];
    std::vector<T> s(a.rbegin(), a.rend());  // (a_{k−1}, …, a_0)
    return companion_matrix(c).pow(n - k + 1).apply(s)[0];
}

// ---- 小さい行列を速く MatrixFixed<T, K>（K×K。大きさがコンパイル時に決まるので展開される） ----
template <class T, int K>
struct MatrixFixed {
    std::array<std::array<T, K>, K> a;
    MatrixFixed() {
        for (auto& row : a) row.fill(math_internal::sr_zero<T>());
    }
    MatrixFixed(std::initializer_list<std::initializer_list<T>> v) : MatrixFixed() {
        if ((int)v.size() != K) math_internal::fail("MatrixFixed", "行の数が K でない");
        int i = 0;
        for (const auto& row : v) {
            if ((int)row.size() != K) math_internal::fail("MatrixFixed", "列の数が K でない");
            int j = 0;
            for (const T& x : row) a[(size_t)i][(size_t)j++] = x;
            ++i;
        }
    }
    static MatrixFixed identity() {
        MatrixFixed r;
        for (int i = 0; i < K; ++i) r.a[(size_t)i][(size_t)i] = math_internal::sr_one<T>();
        return r;
    }
    T& operator()(int i, int j) { return a[(size_t)i][(size_t)j]; }
    const T& operator()(int i, int j) const { return a[(size_t)i][(size_t)j]; }
    std::array<T, K>& operator[](int i) { return a[(size_t)i]; }
    const std::array<T, K>& operator[](int i) const { return a[(size_t)i]; }
    friend MatrixFixed operator+(MatrixFixed x, const MatrixFixed& y) {
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) x.a[i][j] = math_internal::mat_add(x.a[i][j], y.a[i][j]);
        return x;
    }
    friend MatrixFixed operator-(MatrixFixed x, const MatrixFixed& y) {  // 半環では使えない
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) x.a[i][j] = math_internal::mat_sub(x.a[i][j], y.a[i][j]);
        return x;
    }
    friend MatrixFixed operator*(const MatrixFixed& x, const MatrixFixed& y) {
        MatrixFixed r;
        for (int i = 0; i < K; ++i)
            for (int k = 0; k < K; ++k)
                for (int j = 0; j < K; ++j) math_internal::mat_fma(r.a[i][j], x.a[i][k], y.a[k][j]);
        return r;
    }
    friend MatrixFixed operator*(MatrixFixed x, const T& s) {  // スカラー倍
        for (auto& row : x.a)
            for (T& v : row) v = math_internal::mat_mul(v, s);
        return x;
    }
    friend MatrixFixed operator*(const T& s, const MatrixFixed& x) { return x * s; }
    MatrixFixed& operator*=(const T& s) { return *this = *this * s; }
    MatrixFixed& operator+=(const MatrixFixed& y) { return *this = *this + y; }
    MatrixFixed& operator-=(const MatrixFixed& y) { return *this = *this - y; }
    MatrixFixed& operator*=(const MatrixFixed& y) { return *this = *this * y; }
    MatrixFixed pow(long long e) const {
        if (e < 0) math_internal::fail("MatrixFixed::pow", "e ≥ 0 で使う");
        MatrixFixed r = identity(), b = *this;
        for (; e > 0; e >>= 1) {
            if (e & 1) r = r * b;
            if (e > 1) b = b * b;
        }
        return r;
    }
    std::array<T, K> apply(const std::array<T, K>& v) const {
        std::array<T, K> r;
        r.fill(math_internal::sr_zero<T>());
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) math_internal::mat_fma(r[(size_t)i], a[i][j], v[(size_t)j]);
        return r;
    }
    MatrixFixed transpose() const {
        MatrixFixed r;
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) r.a[j][i] = a[i][j];
        return r;
    }
    T trace() const {
        T s = math_internal::sr_zero<T>();
        for (int i = 0; i < K; ++i) s = math_internal::mat_add(s, a[i][i]);
        return s;
    }
    T det() const {  // 行列式（T は mint（mod は素数）か double）  O(K^3)
        static_assert(!std::is_integral_v<T>, "det は mint か double で使う");
        MatrixFixed b = *this;
        T d = T(1);
        for (int c = 0; c < K; ++c) {
            const int p = pivot(b, c);
            if (p < 0) return T(0);
            if (p != c) std::swap(b.a[p], b.a[c]), d = -d;
            d = d * b.a[c][c];
            const T iv = T(1) / b.a[c][c];
            for (int r = c + 1; r < K; ++r) {
                const T f = b.a[r][c] * iv;
                for (int j = c; j < K; ++j) b.a[r][j] = b.a[r][j] - f * b.a[c][j];
            }
        }
        return d;
    }
    Maybe<MatrixFixed> try_inv() const {  // 逆行列（なければ ok = false。T は mint（mod は素数）か double）
        static_assert(!std::is_integral_v<T>, "inv は mint か double で使う");
        MatrixFixed b = *this, r = identity();
        for (int c = 0; c < K; ++c) {
            const int p = pivot(b, c);
            if (p < 0) return Maybe<MatrixFixed>{false, MatrixFixed()};
            std::swap(b.a[p], b.a[c]), std::swap(r.a[p], r.a[c]);
            const T iv = T(1) / b.a[c][c];
            for (int j = 0; j < K; ++j) b.a[c][j] = b.a[c][j] * iv, r.a[c][j] = r.a[c][j] * iv;
            for (int i = 0; i < K; ++i) {
                if (i == c) continue;
                const T f = b.a[i][c];
                for (int j = 0; j < K; ++j) b.a[i][j] = b.a[i][j] - f * b.a[c][j], r.a[i][j] = r.a[i][j] - f * r.a[c][j];
            }
        }
        return Maybe<MatrixFixed>{true, r};
    }
    MatrixFixed inv() const {  // 逆行列（なければ止める）
        const Maybe<MatrixFixed> r = try_inv();
        if (!r.ok) math_internal::fail("MatrixFixed::inv", "逆行列がない（行列式が 0）");
        return r.value;
    }
    friend bool operator==(const MatrixFixed& x, const MatrixFixed& y) { return x.a == y.a; }
    friend bool operator!=(const MatrixFixed& x, const MatrixFixed& y) { return !(x == y); }
    friend std::ostream& operator<<(std::ostream& os, const MatrixFixed& x) {
        for (int i = 0; i < K; ++i) {
            for (int j = 0; j < K; ++j) os << (j ? " " : "") << x.a[i][j];
            if (i + 1 < K) os << '\n';
        }
        return os;
    }

  private:
    static int pivot(const MatrixFixed& b, int c) {  // c 列目で c 行目以降のピボット（double は絶対値が最大のもの。なければ −1）
        if constexpr (std::is_floating_point_v<T>) {
            int p = c;
            for (int r = c + 1; r < K; ++r)
                if (std::abs(b.a[r][c]) > std::abs(b.a[p][c])) p = r;
            return std::abs(b.a[p][c]) < 1e-12 ? -1 : p;
        } else {
            for (int r = c; r < K; ++r)
                if (b.a[r][c] != T(0)) return r;
            return -1;
        }
    }
};
#endif  // GG_MATRIX
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace mat_test {

#define MAT_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[matrix] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

using mint = ModInt998;
template <class T>
static Matrix<T> naive_mul(const Matrix<T>& a, const Matrix<T>& b) {  // 愚直（i, j, k の順）
    Matrix<T> c(a.rows(), b.cols());
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < b.cols(); ++j) {
            T s = math_internal::sr_zero<T>();
            for (int k = 0; k < a.cols(); ++k) s = s + a(i, k) * b(k, j);
            c(i, j) = s;
        }
    return c;
}
static Matrix<mint> rand_mat(std::mt19937_64& rng, int n, int m) {
    Matrix<mint> a(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j) a(i, j) = rng() % 4 == 0 ? mint(0) : mint((long long)rng());
    return a;
}

static void test_mint(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 6), m = (int)rnd(rng, 0, 6), k = (int)rnd(rng, 0, 6);
        const Matrix<mint> A = rand_mat(rng, n, m), B = rand_mat(rng, m, k), C = rand_mat(rng, n, m);
        MAT_CHECK(A * B == naive_mul(A, B));
        MAT_CHECK((A * B).transpose() == B.transpose() * A.transpose());
        const mint s = (long long)rng();
        Matrix<mint> D = A + C, E = A - C, F = A * s, G = s * A;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) MAT_CHECK(D(i, j) == A(i, j) + C(i, j) && E[i][j] == A[i][j] - C[i][j] && F(i, j) == A(i, j) * s && G(i, j) == F(i, j));
        // apply は列ベクトルとの掛け算
        std::vector<mint> v((size_t)m);
        Matrix<mint> V(m, 1);
        for (int j = 0; j < m; ++j) v[(size_t)j] = (long long)rng(), V(j, 0) = v[(size_t)j];
        const std::vector<mint> av = A.apply(v);
        const Matrix<mint> AV = A * V;
        for (int i = 0; i < n; ++i) MAT_CHECK(av[(size_t)i] == AV(i, 0));
        // 累乗
        const Matrix<mint> S = rand_mat(rng, n, n);
        Matrix<mint> P = Matrix<mint>::identity(n);
        const int e = (int)rnd(rng, 0, 12);
        for (int t = 0; t < e; ++t) P = P * S;
        MAT_CHECK(S.pow(e) == P);
        const long long e1 = rnd(rng, 0, 1000000000000000000LL), e2 = rnd(rng, 0, 1000000000000000000LL);
        if (it % 8 == 0) MAT_CHECK(S.pow(e1) * S.pow(e2) == S.pow(e1 + e2));
        mint tr = 0;
        for (int i = 0; i < n; ++i) tr += S(i, i);
        MAT_CHECK(S.trace() == tr);
    }
    MAT_CHECK(Matrix<mint>::identity(3) * Matrix<mint>({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}) == Matrix<mint>({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}));
}

static void test_integer_and_real(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 5), m = (int)rnd(rng, 1, 5), k = (int)rnd(rng, 1, 5);
        Matrix<long long> A(n, m), B(m, k);
        Matrix<double> Ad(n, m), Bd(m, k);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) A(i, j) = rnd(rng, -1000000, 1000000), Ad(i, j) = (double)A(i, j);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < k; ++j) B(i, j) = rnd(rng, -1000000, 1000000), Bd(i, j) = (double)B(i, j);
        const Matrix<long long> C = A * B;
        const Matrix<double> Cd = Ad * Bd;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < k; ++j) {
                Int128 s = 0;
                for (int t = 0; t < m; ++t) s += (Int128)A(i, t) * B(t, j);
                MAT_CHECK((Int128)C(i, j) == s && std::abs(Cd(i, j) - (double)s) <= 1e-6 * (1 + std::abs((double)s)));
            }
    }
    // フィボナッチ：long long で F_90 まで厳密に
    const Matrix<long long> F = Matrix<long long>({{1, 1}, {1, 0}}).pow(90);
    MAT_CHECK(F(0, 1) == 2880067194370816120LL);
}

template <class S>
static void test_semiring_paths(std::mt19937_64& rng, int& cases, bool minimize) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 6), K = (int)rnd(rng, 0, 12);
        Matrix<S> A(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (rng() % 3) A(i, j) = S(rnd(rng, -20, 100));
        // DP：ちょうど t 歩
        const long long NONE = std::numeric_limits<long long>::min();
        std::vector<std::vector<long long>> d((size_t)n, std::vector<long long>((size_t)n, NONE));
        for (int i = 0; i < n; ++i) d[(size_t)i][(size_t)i] = 0;
        for (int t = 0; t < K; ++t) {
            std::vector<std::vector<long long>> nd((size_t)n, std::vector<long long>((size_t)n, NONE));
            for (int i = 0; i < n; ++i)
                for (int u = 0; u < n; ++u) {
                    if (d[(size_t)i][(size_t)u] == NONE) continue;
                    for (int v = 0; v < n; ++v) {
                        if (A(u, v).is_inf()) continue;
                        const long long c = d[(size_t)i][(size_t)u] + A(u, v).v;
                        long long& x = nd[(size_t)i][(size_t)v];
                        if (x == NONE || (minimize ? c < x : c > x)) x = c;
                    }
                }
            d = nd;
        }
        const Matrix<S> P = A.pow(K);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) MAT_CHECK(d[(size_t)i][(size_t)j] == NONE ? P(i, j).is_inf() : (!P(i, j).is_inf() && P(i, j).v == d[(size_t)i][(size_t)j]));
    }
}

template <class S>
static void test_bottleneck(std::mt19937_64& rng, int& cases, bool min_of_max) {  // MinMax / MaxMin
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 6), K = (int)rnd(rng, 0, 12);
        Matrix<S> A(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (rng() % 3) A(i, j) = S(rnd(rng, -50, 50));
        const long long NONE = LLONG_MIN, START = min_of_max ? LLONG_MIN + 1 : LLONG_MAX;  // 0 歩は「辺なし」
        std::vector<std::vector<long long>> d((size_t)n, std::vector<long long>((size_t)n, NONE));
        for (int i = 0; i < n; ++i) d[(size_t)i][(size_t)i] = START;
        for (int t = 0; t < K; ++t) {
            std::vector<std::vector<long long>> nd((size_t)n, std::vector<long long>((size_t)n, NONE));
            for (int i = 0; i < n; ++i)
                for (int u = 0; u < n; ++u) {
                    if (d[(size_t)i][(size_t)u] == NONE) continue;
                    for (int v = 0; v < n; ++v) {
                        if (A(u, v).is_inf()) continue;
                        const long long w = A(u, v).v, c = min_of_max ? std::max(d[(size_t)i][(size_t)u], w) : std::min(d[(size_t)i][(size_t)u], w);
                        long long& x = nd[(size_t)i][(size_t)v];
                        if (x == NONE || (min_of_max ? c < x : c > x)) x = c;
                    }
                }
            d = nd;
        }
        const Matrix<S> P = A.pow(K);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (d[(size_t)i][(size_t)j] == NONE) MAT_CHECK(P(i, j).is_inf());
                else if (d[(size_t)i][(size_t)j] == START) MAT_CHECK(P(i, j) == S::one());
                else MAT_CHECK(!P(i, j).is_inf() && P(i, j).v == d[(size_t)i][(size_t)j]);
            }
    }
}

static void test_bool_and_xor(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {  // OrAnd：ちょうど K 歩で行けるか
        const int n = (int)rnd(rng, 1, 7), K = (int)rnd(rng, 0, 15);
        Matrix<OrAnd> A(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) A(i, j) = rng() % 4 == 0;
        std::vector<std::vector<bool>> r((size_t)n, std::vector<bool>((size_t)n, false));
        for (int i = 0; i < n; ++i) r[(size_t)i][(size_t)i] = true;
        for (int t = 0; t < K; ++t) {
            std::vector<std::vector<bool>> nr((size_t)n, std::vector<bool>((size_t)n, false));
            for (int i = 0; i < n; ++i)
                for (int u = 0; u < n; ++u)
                    if (r[(size_t)i][(size_t)u])
                        for (int v = 0; v < n; ++v)
                            if (A(u, v).v) nr[(size_t)i][(size_t)v] = true;
            r = nr;
        }
        const Matrix<OrAnd> P = A.pow(K);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) MAT_CHECK(P(i, j).v == r[(size_t)i][(size_t)j]);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // XorAnd と mint の漸化式：直接計算と比べる
        const int k = (int)rnd(rng, 1, 6);
        std::vector<XorAnd> c((size_t)k), a((size_t)k);
        std::vector<mint> cm((size_t)k), am((size_t)k);
        for (int i = 0; i < k; ++i) c[(size_t)i] = rng(), a[(size_t)i] = rng(), cm[(size_t)i] = (long long)rng(), am[(size_t)i] = (long long)rng();
        std::vector<XorAnd> seq(a);
        std::vector<mint> sm(am);
        const int N = 150;
        while ((int)seq.size() <= N) {
            XorAnd x = 0;
            mint y = 0;
            const size_t i = seq.size();
            for (int j = 0; j < k; ++j) x = x + c[(size_t)j] * seq[i - 1 - (size_t)j], y += cm[(size_t)j] * sm[i - 1 - (size_t)j];
            seq.push_back(x), sm.push_back(y);
        }
        const long long n = rnd(rng, 0, N);
        MAT_CHECK(recurrence_nth_matrix(c, a, n) == seq[(size_t)n] && recurrence_nth_matrix(cm, am, n) == sm[(size_t)n]);
    }
}

template <int K>
static void test_fixed(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 5000 * TEST_SCALE; ++it, ++cases) {
        MatrixFixed<mint, K> A, B;
        Matrix<mint> a(K, K), b(K, K);
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) {
                const long long x = it % 3 == 0 ? rnd(rng, 0, 2) : (long long)rng(), y = (long long)rng();
                A(i, j) = x, a(i, j) = x, B[i][j] = y, b(i, j) = y;
            }
        const MatrixFixed<mint, K> C = A * B, D = A + B, E = A - B, P = A.pow(it % 50);
        const Matrix<mint> c = a * b, p = a.pow(it % 50);
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) MAT_CHECK(C(i, j) == c(i, j) && D(i, j) == a(i, j) + b(i, j) && E(i, j) == a(i, j) - b(i, j) && P(i, j) == p(i, j));
        // 行列式：置換で展開
        std::vector<int> perm(K);
        std::iota(perm.begin(), perm.end(), 0);
        mint det = 0;
        do {
            int inv = 0;
            for (int i = 0; i < K; ++i)
                for (int j = i + 1; j < K; ++j) inv += perm[(size_t)i] > perm[(size_t)j];
            mint t = inv % 2 ? mint(-1) : mint(1);
            for (int i = 0; i < K; ++i) t *= A(i, perm[(size_t)i]);
            det += t;
        } while (std::next_permutation(perm.begin(), perm.end()));
        MAT_CHECK(A.det() == det && (A * B).det() == A.det() * B.det());
        const Maybe<MatrixFixed<mint, K>> iv = A.try_inv();
        MAT_CHECK(iv.ok == (det != mint(0)));
        if (iv.ok) MAT_CHECK(A * iv.value == (MatrixFixed<mint, K>::identity()) && iv.value * A == (MatrixFixed<mint, K>::identity()));
        std::array<mint, K> v;
        for (auto& x : v) x = (long long)rng();
        const std::array<mint, K> av = A.apply(v);
        for (int i = 0; i < K; ++i) {
            mint s = 0;
            for (int j = 0; j < K; ++j) s += A(i, j) * v[(size_t)j];
            MAT_CHECK(av[(size_t)i] == s);
        }
        MAT_CHECK(A.transpose().det() == det && A.trace() == a.trace());
        const mint sc = (long long)rng();
        MatrixFixed<mint, K> S1 = A * sc, S2 = sc * A, S3 = A;
        S3 *= sc;
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) MAT_CHECK(S1(i, j) == A(i, j) * sc && S2 == S1 && S3 == S1);
    }
    // double
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        MatrixFixed<double, K> A;
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) A(i, j) = (double)rnd(rng, -100, 100) / 7;
        const double d = A.det();
        const Maybe<MatrixFixed<double, K>> iv = A.try_inv();
        if (std::abs(d) > 1e-6) {
            MAT_CHECK(iv.ok);
            const MatrixFixed<double, K> I = A * iv.value;
            for (int i = 0; i < K; ++i)
                for (int j = 0; j < K; ++j) MAT_CHECK(std::abs(I(i, j) - (i == j ? 1.0 : 0.0)) < 1e-8);
            MAT_CHECK(std::abs(iv.value.det() * d - 1) < 1e-8);
        }
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_mint(rng, cases);
    test_integer_and_real(rng, cases);
    test_semiring_paths<MinPlus<long long>>(rng, cases, true);
    test_semiring_paths<MaxPlus<long long>>(rng, cases, false);
    test_bottleneck<MinMax<long long>>(rng, cases, true);
    test_bottleneck<MaxMin<long long>>(rng, cases, false);
    test_bool_and_xor(rng, cases);
    test_fixed<1>(rng, cases);
    test_fixed<2>(rng, cases);
    test_fixed<3>(rng, cases);
    test_fixed<4>(rng, cases);
    std::cout << "[matrix] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[matrix] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    using M107 = ModInt107;
    Matrix<M107> A(300, 300), B(300, 300);
    for (int i = 0; i < 300; ++i)
        for (int j = 0; j < 300; ++j) A(i, j) = (long long)rng(), B(i, j) = (long long)rng();
    auto t0 = std::chrono::steady_clock::now();
    const Matrix<M107> C = A * B;
    const long long t_mul = ms(t0);
    Matrix<M107> G(50, 50);
    for (int i = 0; i < 50; ++i)
        for (int j = 0; j < 50; ++j) G(i, j) = (long long)(rng() % 2);
    t0 = std::chrono::steady_clock::now();
    const Matrix<M107> GP = G.pow(1000000000000000000LL);
    const long long t_pow = ms(t0);
    t0 = std::chrono::steady_clock::now();
    MatrixFixed<M107, 2> F = {{1, 1}, {1, 0}};
    M107 s = 0;
    for (int i = 0; i < 100000; ++i) s += F.pow((long long)(rng() >> 4))(0, 1);
    const long long t_fixed = ms(t0);
    Matrix<MinPlus<long long>> D(100, 100);
    for (int i = 0; i < 100; ++i)
        for (int j = 0; j < 100; ++j)
            if (rng() % 2) D(i, j) = rnd(rng, 1, 1000000);
    t0 = std::chrono::steady_clock::now();
    const Matrix<MinPlus<long long>> DP = D.pow(1000);
    const long long t_mp = ms(t0);
    std::cout << "[matrix] benchmark: Matrix<mint> 300×300 の掛け算 " << t_mul << " ms, 50×50 の 10^18 乗 " << t_pow << " ms, MatrixFixed<mint, 2> の累乗 × 10^5 " << t_fixed
              << " ms, Matrix<MinPlus> 100×100 の 1000 乗 " << t_mp << " ms  (check " << (C(1, 2) + GP(3, 4) + s).val() % 7 + DP(0, 0).v % 7 << ")\n";
}

}  // namespace mat_test


// デモで行列を 1 行に表示する（行を " / " でつなぐ。テンプレートの cout << A は 1 行ずつ改行）
template <class T>
static std::string to_s(const Matrix<T>& a) {
    std::ostringstream os;
    os << "[";
    for (int i = 0; i < a.rows(); ++i) {
        os << (i ? " / " : "");
        for (int j = 0; j < a.cols(); ++j) os << (j ? " " : "") << a(i, j);
    }
    return os.str() + "]";
}
template <class T, int K>
static std::string to_s(const MatrixFixed<T, K>& a) {
    std::ostringstream os;
    os << "[";
    for (int i = 0; i < K; ++i) {
        os << (i ? " / " : "");
        for (int j = 0; j < K; ++j) os << (j ? " " : "") << a(i, j);
    }
    return os.str() + "]";
}

// =============================================================================
//  使い方デモ
// =============================================================================
namespace mat_demo {

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
#define SHOW(expr, comment) mat_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    mat_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("行列（matrix） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Matrix<T>：A(i, j) / A[i][j] / + - * / A.pow(e) / A.apply(v) / A.transpose() / A.trace() / identity(n)\n"
              << "  半環：MinPlus<T> / MaxPlus<T> / MinMax<T> / MaxMin<T> / XorAnd / OrAnd（Matrix<MinPlus<long long>> のように）\n"
              << "  companion_matrix(c) / recurrence_nth_matrix(c, a, n)    線形漸化式\n"
              << "  MatrixFixed<T, K>：K×K の小さい行列（det / inv も）\n";

    api("【1】 Matrix<T> A = {{…}, {…}} / A * B / A.pow(e) / A.apply(v)", "Matrix<T> / vector<T>",
        "行列の計算。T は mint・long long（溢れたら止める）・double など", "掛け算 O(n^3)・累乗 O(n^3 log e)");
    note("using mint = ModInt107;");
    using mint = ModInt107;
    RUN(Matrix<mint> A = {{1, 1}, {1, 0}});
    SHOW(A.pow(10), "フィボナッチ：右上が F_10 = 55");
    SHOW(A.pow(1000000000000000000LL)(0, 1), "F_{10^18} mod 10^9+7");
    SHOW(A * A, "掛け算");
    SHOW((A.apply({mint(1), mint(0)})), "A × (1, 0)^T");
    SHOW(Matrix<long long>({{1, 2, 3}, {4, 5, 6}}).transpose(), "転置（2×3 → 3×2）");
    SHOW(Matrix<long long>({{1, 2}, {3, 4}}).trace(), "対角成分の和");
    {
        // EDPC R：隣接行列の K 乗の全成分の和 = 長さ K の道の数
        auto walks = [](const std::vector<std::vector<int>>& adj, long long K) {
            const int n = (int)adj.size();
            Matrix<mint> G(n, n);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) G(i, j) = adj[(size_t)i][(size_t)j];
            const Matrix<mint> P = G.pow(K);
            mint s = 0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) s += P(i, j);
            return s;
        };
        note("EDPC R（Walk）: 隣接行列の K 乗の全成分の和（グラフは graph_grid_v2 で読んで隣接行列にしてもよい）");
        std::cout << "    例 1 → " << walks({{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {1, 0, 0, 0}}, 2) << "、例 2 → " << walks({{0, 1, 0}, {1, 0, 1}, {0, 0, 0}}, 3)
                  << "、例 3 → " << walks({{0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0}}, 2)
                  << "、例 4 → " << walks({{0}}, 1) << "\n";
    }
    {
        // ABC293 E：[[A, 1], [0, 1]]^X の右上 = Σ_{i<X} A^i（mod M は合成数のこともある → DynModInt）
        auto solve = [](long long a, long long X, long long M) {
            DynModInt::set_mod(M);
            Matrix<DynModInt> T = {{DynModInt(a), DynModInt(1)}, {DynModInt(0), DynModInt(1)}};
            return T.pow(X)(0, 1);
        };
        note("ABC293 E（Geometric Sequence Sum）: [[A, 1], [0, 1]]^X の右上（Matrix<DynModInt>。00_math_base の geom_sum_mod でも）");
        std::cout << "    (3, 4, 7) → " << solve(3, 4, 7) << "、(8, 10, 9) → " << solve(8, 10, 9) << "、(10^9, 10^12, 998244353) → " << solve(1000000000, 1000000000000LL, 998244353) << "\n";
    }

    api("【2】 companion_matrix(c) / recurrence_nth_matrix(c, a, n)", "Matrix<T> / T",
        "a_i = c_0 a_{i−1} + … + c_{k−1} a_{i−k}（a_0〜a_{k−1} が a）の a_n。n は 0 から数える", "O(k^3 log n)");
    SHOW(recurrence_nth_matrix<mint>({1, 1}, {0, 1}, 90), "フィボナッチ F_90 mod 10^9+7");
    SHOW(recurrence_nth_matrix<mint>({1, 1, 1}, {0, 0, 1}, 1000000000000000000LL), "トリボナッチの 10^18 項目");
    SHOW(companion_matrix<long long>({2, 3}), "a_i = 2 a_{i−1} + 3 a_{i−2} の行列");
    note("k が大きい（数百以上）ときは 06_linear_recurrence（きたまさ法・Bostan–Mori）のほうが速い");

    api("【3】 半環：Matrix<MinPlus<long long>> / Matrix<OrAnd> / Matrix<XorAnd>", "Matrix<半環>",
        "+ と × を min と +（最短路）、OR と AND（行けるか）、XOR と AND（ビットの漸化式）にした行列の累乗", "O(n^3 log K)");
    RUN(Matrix<MinPlus<long long>> D(3, 3); D(0, 1) = 5; D(1, 2) = 1; D(2, 0) = 2; D(1, 0) = 10);
    SHOW(D.pow(4)(0, 1), "0 → 1 にちょうど 4 歩で行く最短：0→1→2→0→1 = 5+1+2+5");
    SHOW(D.pow(2)(0, 2), "ちょうど 2 歩：0→1→2 = 6");
    SHOW(D.pow(2)(0, 0), "ちょうど 2 歩で 0 → 0：0→1→0 = 15");
    SHOW(D.pow(1)(0, 2), "1 歩では行けない");
    RUN(Matrix<MinMax<int>> W(3, 3); W(0, 1) = 7; W(1, 2) = 3; W(0, 2) = 9; W(2, 2) = 1);
    SHOW(W.pow(2)(0, 2), "MinMax：ちょうど 2 歩で 0 → 2、通る辺の最大を最小に（0→1→2 なら max(7, 3) = 7）");
    RUN(Matrix<MaxMin<int>> Cap(3, 3); Cap(0, 1) = 7; Cap(1, 2) = 3; Cap(0, 2) = 9; Cap(2, 2) = 1);
    SHOW(Cap.pow(2)(0, 2), "MaxMin：通る辺の最小を最大に（0→1→2 なら 3、0→2→2 なら 1）");
    RUN(Matrix<OrAnd> R = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}});
    SHOW(R.pow(5)(0, 2), "ちょうど 5 歩で 0 → 2 に行けるか（3 歩で一周）");
    {
        // ABC009 D：A_{N+K} = XOR_j (C_j AND A_{N+K−j}) の M 番目（1 から数える）
        auto solve = [](const std::vector<unsigned long long>& as, const std::vector<unsigned long long>& cs, long long M) {
            std::vector<XorAnd> a(as.begin(), as.end()), c(cs.begin(), cs.end());
            return recurrence_nth_matrix(c, a, M - 1).v;
        };
        note("ABC009 D（漸化式）: XorAnd で recurrence_nth_matrix(C, A, M − 1)");
        std::cout << "    例 1 → " << solve({10, 20, 30}, {7, 19, 13}, 5) << "、例 2 → "
                  << solve({2345678901ULL, 1001001001ULL, 3333333333ULL, 3141592653ULL, 1234567890ULL}, {2147483648ULL, 2147483647ULL, 4294967295ULL, 4294967294ULL, 3434343434ULL}, 100) << "\n";
    }

    api("【4】 MatrixFixed<T, K>：M.pow(e) / M.det() / M.inv() / M.try_inv()", "MatrixFixed / T / Maybe",
        "K×K の小さい行列（大きさがコンパイル時に決まるので速い）。行列式・逆行列は T が mint か double", "O(K^3)");
    RUN(MatrixFixed<mint, 2> F = {{1, 1}, {1, 0}});
    SHOW(F.pow(10)(0, 1), "F_10");
    SHOW(F.det(), "1·0 − 1·1 = −1（mod 10^9+7）");
    SHOW(F.inv(), "逆行列 [[0, 1], [1, −1]]");
    RUN(MatrixFixed<double, 3> R3 = {{2, 0, 0}, {0, 3, 0}, {1, 0, 4}});
    SHOW(R3.det(), "2·3·4");
    SHOW(R3.try_inv().ok, "逆行列がある");
    SHOW((MatrixFixed<double, 2>{{1, 2}, {2, 4}}.try_inv().ok), "行列式 0 なので逆行列なし");
}

}  // namespace mat_demo

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
    mat_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mat_test::benchmark();    // 速度計測
    mat_demo::demo();         // 使い方デモ
    return 0;
}
