// =============================================================================
//  modint（ModInt<MOD>：定数の mod / DynModInt：入力で決まる mod）   --- 数学第二世代
// -----------------------------------------------------------------------------
//  分類       : 00 基盤（組合せ・行列・畳み込みなど mod で数えるテンプレートはこれを同梱している）
//  依存       : 00_int_common, 00_math_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。ACL は使わない
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  使い方（mint という名前はテンプレートでは定義しない。使う側で 1 行書く）
//      using mint = ModInt998;          // 998244353。10^9+7 なら ModInt107、ほかの定数なら ModInt<M>
//      using mint = DynModInt;          // mod が入力で決まるとき。mint::set_mod(M); を最初に呼ぶ
//      mint a = -3, b = 10;  cout << a * b / 7 << '\n';   // val() を書かずに出せる
//
//  機能一覧（ModInt と DynModInt で同じ）                                      計算量
//    mint x = 整数                   int・long long・unsigned・Int128 から（負の数でもよい）     O(1)
//    mint::raw(x)                    0 ≤ x < mod が分かっているときの速い作り方
//    + - * / 単項- ++ -- += -= *= /= 整数と混ざった式も（a * 2 + 1 など）                      O(1)（/ は逆元）
//    == != < > <= >=                 値（0〜mod-1）での比較（< は map のキー・sort 用）
//    x.val() / mint::mod()           値（int）/ mod（int）
//    x.pow(e)                        x^e（e は long long。負なら逆元の (-e) 乗）                O(log |e|)
//                                    （指数が N − 2 などの式で負になりうるときは注意：エラーにならず逆元の累乗になる）
//    x.inv() / x.try_inv()           逆元（mod が合成数でも gcd = 1 なら。ないと止まる / Maybe で返す）  O(log mod)
//    x.is_unit()                     逆元があるか（gcd(x, mod) = 1）
//    if (x)                          0 でないか
//    cout << x / cin >> x            そのまま出す / 読む（long long を超える長さの数も読める）
//    unordered_map<mint, V>          std::hash を用意してある
//    ModInt<M>::IS_PRIME             M が素数か（コンパイル時の定数）
//    DynModInt::set_mod(M)           mod を設定（1 ≤ M < 2^31。素数でなくてよい・偶数でもよい）
//    DynModInt::mod_is_prime()       今の mod が素数か
//    DynModIntT<ID>                  ID を変えると別の mod を同時に持てる（DynModInt = DynModIntT<0>）
//    rational_reconstruct(x)         RationalGuess{ok, num, den}：mint の値を小さい分数に戻す（デバッグ用）
//
//  注意
//    - mod は 1 以上 2^31 未満（int の範囲）。それより大きい mod は 00_math_base の mod_mul / mod_pow を使う
//    - 内部はモンゴメリ乗算（奇数の mod）。偶数の mod は普通の剰余（DynModInt は Barrett）で、どちらも正しく動く
//    - DynModInt の set_mod を途中で変えると、それまでに作った値は意味がなくなる（作り直す）
//    - 0 の逆元・mod と互いに素でない数の逆元は、メッセージを出して止まる（try_inv なら止まらない）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mint_test::random_test()  … 128 ビットの愚直計算と全機能を突き合わせるランダムテスト（いろいろな mod で）
//    mint_test::benchmark()    … 速度計測
//    mint_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
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
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace mint_test {

#define MINT_CHECK(...)                                                                          \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[modint] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static long long refmod(Int128 x, long long m) {
    x %= m;
    if (x < 0) x += m;
    return (long long)x;
}
static long long ref_pow(long long a, unsigned long long e, long long m) {
    Int128 r = 1 % m, x = refmod(a, m);
    for (; e; e >>= 1, x = x * x % m)
        if (e & 1) r = r * x % m;
    return (long long)r;
}
static bool brute_is_prime(long long n) {
    if (n < 2) return false;
    for (long long d = 2; d * d <= n; ++d)
        if (n % d == 0) return false;
    return true;
}
static std::string rnd_digits(std::mt19937_64& rng, int len) {
    std::string s;
    for (int i = 0; i < len; ++i) s += char('0' + rng() % 10);
    return s;
}

// M の全機能を、128 ビットの愚直計算と突き合わせる（M::mod() は呼ぶ前に決まっていること）
template <class M>
static void test_ops(std::mt19937_64& rng, int& cases, int iters) {
    const long long m = M::mod();
    for (int it = 0; it < iters; ++it, ++cases) {
        const long long x = rnd_ll(rng), y = it % 7 == 0 ? x / 4 * (long long)rnd(rng, -2, 2) : rnd_ll(rng);  // 7 回に 1 回は x と関係のある値
        const M a = x, b = y;
        const long long rx = refmod(x, m), ry = refmod(y, m);
        MINT_CHECK(a.val() == rx && b.val() == ry && 0 <= a.val() && a.val() < m);
        MINT_CHECK((a + b).val() == refmod((Int128)x + y, m) && (a - b).val() == refmod((Int128)x - y, m));
        MINT_CHECK((a * b).val() == refmod((Int128)rx * ry, m) && (-a).val() == refmod(-(Int128)x, m) && (+a) == a);
        MINT_CHECK((a * 3 + 1).val() == refmod((Int128)rx * 3 + 1, m) && (2 - a).val() == refmod(2 - (Int128)x, m) && (a == b) == (rx == ry));
        MINT_CHECK((a != b) == (rx != ry) && (a < b) == (rx < ry) && (a <= b) == (rx <= ry) && (a > b) == (rx > ry) && (a >= b) == (rx >= ry));
        MINT_CHECK((bool)a == (rx != 0) && (a == x) && (b == y));
        M c = a;
        c += b, c -= a, c *= a, c -= b * a;
        MINT_CHECK(c == 0);
        c = a;
        MINT_CHECK((c++).val() == rx && c.val() == refmod((Int128)x + 1, m) && (--c).val() == rx && (++c).val() == refmod((Int128)x + 1, m) && (c--) == a + 1 && c == a);
        // 別の型から
        const unsigned long long ux = (unsigned long long)x;
        const Int128 big = (Int128)x * (long long)(rng() >> 1) + (Int128)y;
        const UInt128 ubig = ((UInt128)rng() << 64) | rng();
        MINT_CHECK(M(ux).val() == (long long)(ux % (unsigned long long)m) && M(big).val() == refmod(big, m) && M(ubig).val() == (long long)(ubig % (UInt128)m));
        MINT_CHECK(M((int)x).val() == refmod((int)x, m) && M((unsigned)x).val() == (long long)((unsigned)x % (unsigned long long)m) && M(true).val() == 1 % m && M((char)7).val() == 7 % m);
        MINT_CHECK(M::raw((int)rx) == a);
        // 累乗・逆元
        const unsigned long long e = it % 3 ? (unsigned long long)rnd(rng, 0, 30) : rng() >> 1;
        MINT_CHECK(a.pow((long long)e).val() == ref_pow(rx, e, m));
        const bool unit = std::gcd(rx, m) == 1;
        MINT_CHECK(a.is_unit() == unit);
        const Maybe<M> inv = a.try_inv();
        MINT_CHECK(inv.ok == unit);
        if (unit) {
            MINT_CHECK((a * inv.value).val() == 1 % m && a.inv() == inv.value && (b / a * a) == b && (b / a).val() == refmod((Int128)ry * inv.value.val(), m));
            MINT_CHECK(a.pow(-(long long)e) * a.pow((long long)e) == M(1) && a.pow(-1) == inv.value);
        }
        // 入出力
        std::ostringstream os;
        os << a << ' ' << -b;
        MINT_CHECK(os.str() == std::to_string(rx) + " " + std::to_string(refmod(-(Int128)y, m)));
        const std::string d = rnd_digits(rng, (int)rnd(rng, 1, 36));
        Int128 dv = 0;
        for (char ch : d) dv = dv * 10 + (ch - '0');
        std::istringstream is(d + " -" + d + " +" + d);
        M r1, r2, r3;
        is >> r1 >> r2 >> r3;
        MINT_CHECK(r1.val() == refmod(dv, m) && r2.val() == refmod(-dv, m) && r3 == r1);
    }
    // 長い計算の連鎖（遅延簡約の値がずれないか）
    {
        M acc = 1, acc2 = 0;
        long long ra = 1 % m, rb = 0;
        for (int i = 0; i < 20000; ++i, ++cases) {
            const long long v = rnd_ll(rng);
            const long long rv = refmod(v, m);
            switch (rng() % 6) {
                case 0: acc *= v, ra = refmod((Int128)ra * rv, m); break;
                case 1: acc += v, ra = refmod((Int128)ra + rv, m); break;
                case 2: acc -= v, ra = refmod((Int128)ra - rv, m); break;
                case 3: acc = -acc, ra = refmod(-(Int128)ra, m); break;
                case 4: acc2 += acc * acc, rb = refmod((Int128)rb + (Int128)ra * ra, m); break;
                default: acc2 = acc2 * acc2 - acc, rb = refmod((Int128)rb * rb - ra, m); break;
            }
            MINT_CHECK(acc.val() == ra && acc2.val() == rb);
        }
    }
    // ハッシュ・map・sort
    {
        std::unordered_set<M> hs;
        std::map<M, int> mp;
        std::vector<long long> vals;
        for (int i = 0; i < 2000; ++i) {
            const long long v = rnd(rng, -100000, 100000);
            hs.insert(M(v)), mp[M(v)]++, vals.push_back(refmod(v, m));
        }
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
        MINT_CHECK(hs.size() == vals.size() && mp.size() == vals.size() && mp.begin()->first.val() == vals[0]);
        ++cases;
    }
}

template <unsigned MOD>
static void test_static(std::mt19937_64& rng, int& cases) {
    static_assert(ModInt<MOD>::mod() == (int)MOD, "mod");
    MINT_CHECK(ModInt<MOD>::IS_PRIME == brute_is_prime(MOD));
    test_ops<ModInt<MOD>>(rng, cases, 3000 * TEST_SCALE);
}

static void test_dynamic(std::mt19937_64& rng, int& cases) {
    std::vector<long long> mods = {998244353, 1000000007, 1, 2, 3, 4, 6, 7, 2147483647, 2147483646, 1073741823, 1073741824, 1073741825,
                                   1073741789, 1000000000, 999999999, 65536, 167772161};
    for (int i = 0; i < 60; ++i) mods.push_back(rnd(rng, 1, i % 3 == 0 ? 100 : 2147483647));
    for (long long m : mods) {
        DynModInt::set_mod(m);
        MINT_CHECK(DynModInt::mod() == m && DynModInt::mod_is_prime() == brute_is_prime(m));
        test_ops<DynModInt>(rng, cases, 500 * TEST_SCALE);
    }
    // 2 つの mod を同時に
    DynModInt::set_mod(1000000007);
    DynModIntT<1>::set_mod(998244353);
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const long long x = rnd_ll(rng), y = rnd_ll(rng);
        const DynModInt a = x, b = y;
        const DynModIntT<1> c = x, d = y;
        MINT_CHECK((a * b).val() == refmod((Int128)refmod(x, 1000000007) * refmod(y, 1000000007), 1000000007));
        MINT_CHECK((c * d).val() == refmod((Int128)refmod(x, 998244353) * refmod(y, 998244353), 998244353));
        MINT_CHECK((a * b).val() == (ModInt107(x) * ModInt107(y)).val() && (c + d).val() == (ModInt998(x) + ModInt998(y)).val());
    }
    DynModInt::set_mod(998244353);
}

template <class M>
static void test_reconstruct(std::mt19937_64& rng, int& cases) {
    const long long m = M::mod(), B = isqrt(m / 2);
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long lim = it % 2 ? std::min<long long>(B, 1000) : B;
        long long p = rnd(rng, -lim, lim), q = rnd(rng, 1, std::max<long long>(1, lim));
        const long long g = std::gcd(std::llabs(p), q);
        p /= g, q /= g;
        if (std::gcd(q, m) != 1) continue;
        const RationalGuess r = rational_reconstruct(M(p) / M(q));
        MINT_CHECK(r.ok && r.num == p && r.den == q);
        // 小さい分数でないかもしれない値：見つかったら条件を満たすこと
        const M z = rnd_ll(rng);
        const RationalGuess s = rational_reconstruct(z);
        if (s.ok) MINT_CHECK(s.den >= 1 && s.den <= B && std::llabs(s.num) <= B && M(s.num) == z * M(s.den));
    }
    std::ostringstream os;
    os << rational_reconstruct(M(1) / M(3)) << ' ' << rational_reconstruct(M(5));
    MINT_CHECK(os.str() == "1/3 5");
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_static<998244353>(rng, cases);
    test_static<1000000007>(rng, cases);
    test_static<167772161>(rng, cases);
    test_static<2147483647>(rng, cases);   // 素数・2^30 以上（モンゴメリ・遅延簡約なし）
    test_static<1073741823>(rng, cases);   // 合成数・2^30 未満（遅延簡約）
    test_static<1073741825>(rng, cases);   // 合成数・2^30 以上
    test_static<999999999>(rng, cases);
    test_static<1000000000>(rng, cases);   // 偶数
    test_static<65536>(rng, cases);
    test_static<7>(rng, cases);
    test_static<3>(rng, cases);
    test_static<2>(rng, cases);
    test_static<1>(rng, cases);
    test_dynamic(rng, cases);
    test_reconstruct<ModInt998>(rng, cases);
    test_reconstruct<ModInt107>(rng, cases);
    DynModInt::set_mod(1000000000);
    test_reconstruct<DynModInt>(rng, cases);
    DynModInt::set_mod(998244353);
    std::cout << "[modint] random_test: OK (" << cases << " cases)\n";
}

template <class M>
static void bench_one(const char* name) {
    std::mt19937_64 rng(1);
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    unsigned long long chk = 0;
    auto t0 = std::chrono::steady_clock::now();
    {
        M x = 12345, y = 678901, c = 3;
        for (int i = 0; i < 20000000; ++i) x = x * y + c;
        chk += (unsigned long long)x.val();
    }
    const long long t_chain = ms(t0);
    const int n = 1 << 20;
    std::vector<M> a(n), b(n);
    for (int i = 0; i < n; ++i) a[i] = (long long)(rng() >> 2), b[i] = (long long)(rng() >> 2);
    t0 = std::chrono::steady_clock::now();
    {
        M s = 0;
        for (int r = 0; r < 16; ++r)
            for (int i = 0; i < n; ++i) s += a[i] * b[i];
        chk += (unsigned long long)s.val();
    }
    const long long t_dot = ms(t0);
    t0 = std::chrono::steady_clock::now();
    {
        const int N = 10000000;
        std::vector<M> f(N + 1), g(N + 1);
        f[0] = 1;
        for (int i = 1; i <= N; ++i) f[i] = f[i - 1] * i;
        g[N] = f[N].inv();
        for (int i = N; i > 0; --i) g[i - 1] = g[i] * i;
        chk += (unsigned long long)g[7].val();
    }
    const long long t_fact = ms(t0);
    t0 = std::chrono::steady_clock::now();
    {
        M s = 0;
        for (int i = 1; i < 1000000; ++i) s += a[i].pow(1000000000 + i) + a[i].inv();
        chk += (unsigned long long)s.val();
    }
    const long long t_pow = ms(t0);
    std::cout << "[modint] benchmark: " << name << "  x = x*y+c × 2×10^7 " << t_chain << " ms, 内積 1.6×10^7 回 " << t_dot
              << " ms, 階乗と逆元の表 N = 10^7 " << t_fact << " ms, pow + inv × 10^6 " << t_pow << " ms  (check " << chk % 7 << ")\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[modint] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    bench_one<ModInt998>("ModInt998    ");
    bench_one<ModInt107>("ModInt107    ");
    DynModInt::set_mod(998244353);
    bench_one<DynModInt>("DynModInt(998244353)");
    DynModInt::set_mod(2147483647);
    bench_one<DynModInt>("DynModInt(2^31-1)");
    DynModInt::set_mod(998244353);
}

}  // namespace mint_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace mint_demo {

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
#define SHOW(expr, comment) mint_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    mint_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("modint 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧（ModInt と DynModInt で同じ）\n"
              << "  using mint = ModInt998;  / ModInt107 / ModInt<M> / DynModInt（mint::set_mod(M)）   使う側で 1 行書く\n"
              << "  mint x = 整数 / mint::raw(x)            作る（負の数・long long・Int128 から）\n"
              << "  + - * / ++ -- 比較                      int と同じように（整数と混ざってもよい）\n"
              << "  x.val() / mint::mod()                   値 / mod\n"
              << "  x.pow(e) / x.inv() / x.try_inv() / x.is_unit()   累乗（負の指数も）/ 逆元 / 止まらない逆元 / 逆元があるか\n"
              << "  cout << x / cin >> x / unordered_map<mint, V>    入出力 / ハッシュ\n"
              << "  rational_reconstruct(x)                 mint の値を小さい分数に戻す（デバッグ用）\n";

    api("【1】 using mint = ModInt998;  mint a = x;  + - * /  cout << a", "mint",
        "998244353 で割った余りで計算する整数。int・long long と同じように書ける。val() なしで cout できる", "O(1)（/ は O(log mod)）");
    using mint = ModInt998;
    RUN(mint a = -3; mint b = 10);
    SHOW(a, "負の数は mod を足した値（998244350）");
    SHOW(a * b + 5, "整数と混ざってもよい");
    SHOW(b / 4, "割り算は逆元を掛ける（10 / 4 = 5/2）");
    SHOW(b / 4 * 4, "掛け戻すと 10");
    SHOW(mint(1000000000000000000LL), "long long から（10^18 mod 998244353）");
    SHOW(mint::mod(), "mod");
    RUN(mint c = 5; c *= c; c -= 1; ++c);
    SHOW(c, "複合代入・++ も");
    SHOW(c.val() + 1LL, "値を整数として使うとき val()");

    api("【2】 x.pow(e) / x.inv()", "mint", "累乗（e は long long。負なら逆元の累乗）/ 逆元", "O(log e) / O(log mod)");
    SHOW(mint(2).pow(10), "1024");
    SHOW(mint(2).pow(-1), "2 の逆元（499122177）");
    SHOW(mint(3).inv() * 3, "逆元を掛けると 1");
    {
        // 典型90 069: N 個のブロックを K 色で、距離 2 以内の色がすべて違う塗り方。K × (K-1) × (K-2)^(N-2)
        using mint107 = ModInt107;
        const long long N = 1000000000000000000LL, K = 3;
        const mint107 ans = N == 1 ? mint107(K) : mint107(K) * (K - 1) * mint107(K - 2).pow(N - 2);
        note("典型90 069（N ≤ 10^18 の累乗）: K × (K − 1) × (K − 2)^(N − 2) を mod 10^9+7 で（ModInt107）");
        std::cout << "    N = 10^18, K = 3 → " << ans << "\n";
        note("注意：pow の指数が負だと逆元の累乗になる。N − 2 のような式が負になりうるとき（N = 1）は、上のように場合分けする");
    }

    api("【3】 DynModInt::set_mod(M)", "なし", "mod が入力で決まるとき。素数でなくてよい・偶数でもよい（逆元は mod と互いに素なときだけ）", "O(1)");
    {
        using dmint = DynModInt;
        RUN(dmint::set_mod(100); dmint x = 7; dmint y = 15);
        SHOW(x * y, "105 mod 100");
        SHOW(x.inv(), "7 × 43 = 301 ≡ 1 (mod 100)");
        SHOW(y.is_unit(), "15 と 100 は互いに素でない");
        SHOW(y.try_inv().ok, "逆元はない（inv() を呼ぶと止まる）");
        SHOW(dmint::mod_is_prime(), "100 は素数でない");
        dmint::set_mod(998244353);
        note("ABC425 E のように mod が素数でないときは、割り算を使わない形（パスカルの三角形など。04_binomial）で数える。");
        note("2 つの mod を同時に使うなら DynModIntT<1>::set_mod(M2) で別の型を作る。");
    }

    api("【4】 期待値・確率を mod で", "mint", "確率 p/100 は mint(p) / 100。分数のまま計算した値と同じものが mod で求まる", "-");
    {
        // ABC280 E: 体力 N。確率 P% で 2、それ以外で 1 のダメージ。体力が 0 以下になるまでの攻撃回数の期待値
        const int N = 3, P = 10;
        const mint p = mint(P) / 100;
        std::vector<mint> dp(N + 1);
        for (int i = 1; i <= N; ++i) dp[i] = 1 + p * (i >= 2 ? dp[i - 2] : mint(0)) + (1 - p) * dp[i - 1];
        note("ABC280 E: dp[i] = 1 + p·dp[i−2] + (1 − p)·dp[i−1]（p = P/100 を mint で）");
        std::cout << "    N = 3, P = 10 → " << dp[N] << "\n";
        SHOW(rational_reconstruct(dp[N]), "分数に戻すと 281/100 = 2.81（手計算と合うか確かめられる）");
    }

    api("【5】 rational_reconstruct(x)", "RationalGuess{ok, num, den}", "mint の値を |分子|, 分母 ≤ √(mod/2) の分数に戻す（デバッグ用。見つからなければ ok = false）", "O(log mod)");
    SHOW(rational_reconstruct(mint(1) / 3), "1/3");
    SHOW(rational_reconstruct(mint(-7) / 12), "-7/12");
    SHOW(rational_reconstruct(mint(123456789)), "小さい分数がない値");

    api("【6】 比較・ハッシュ・入力", "-", "比較は 0〜mod-1 の値で（sort・map に使える）。unordered_map<mint, V> もそのまま。cin は長い数も読める", "-");
    {
        RUN(std::unordered_map<mint, int> cnt; cnt[mint(5)]++; cnt[mint(998244358)]++);
        SHOW(cnt[mint(5)], "998244358 ≡ 5 なので同じキー");
        std::istringstream is("123456789012345678901234567890");
        mint big;
        is >> big;
        SHOW(big, "30 桁の数を読んで mod を取る（long long に入らない入力）");
        SHOW(mint(3) < mint(-1), "値での比較（-1 は 998244352）");
        SHOW(ModInt998::IS_PRIME, "mod が素数か（コンパイル時）");
        SHOW(mint::raw(5), "0 ≤ x < mod と分かっているとき（剰余を取らないので少し速い）");
    }
    note("組合せ（nCr）は 04_binomial の Binomial<mint>、行列は 05_matrix の Matrix<mint>、畳み込みは 07_convolution。");
    note("mod が 2^31 以上（10^18 など）のときは modint ではなく 00_math_base の mod_mul / mod_pow を使う。");
}

}  // namespace mint_demo

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
    mint_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mint_test::benchmark();    // 速度計測
    mint_demo::demo();         // 使い方デモ
    return 0;
}
