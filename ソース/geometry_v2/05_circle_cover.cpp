// =============================================================================
//  点を円で囲む（最小包含円・半径 r の円で最大何点を囲めるか）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 円
//  依存       : 00_int_common, 00_geo_base, 02_line, 05_circle（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - 答えは実数（CircleD / PointD）。「円の内側」は円周の上を含み、半径に合わせた相対誤差 GEO_EPS で判定する
//      （座標が 10^9 でも、近い点どうしの小さい円を正しく求める。Real が double なら座標の大きさ × 10^-14 程度の誤差は残る）
//    - 点は Point（整数）でも PointD（実数）でもよい。整数なら外接円の中心は 128 ビットで求める
//
//  機能一覧                                                                    計算量
//    min_enclosing_circle(ps)          全部の点を囲む最小の円（CircleD。点がなければ r = -1）       期待 O(N)
//    max_points_in_circle(ps, r)       半径 r の円で囲める点の最大の数と、その円の中心              O(N^2 log N)
//                                      CircleCover{count, center}（点がなければ count = 0）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    cv_test::random_test()  … 全部の候補の円を調べる愚直な実装と突き合わせるランダムテスト
//    cv_test::benchmark()    … 速度計測
//    cv_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
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
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
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
#ifndef GG_GEO_BASE
#define GG_GEO_BASE
// ---- 幾何基盤（幾何第二世代）  点 Point（整数）/ PointD（実数）、直線・線分・円、向き、回転、角度、入力 ----
using Coord = long long;    // 整数の座標（|x|, |y| ≤ 10^9 まで正しい）
using Real = long double;   // 実数（double に変えてもよい）
constexpr Real GEO_EPS = 1e-9;  // 実数の比較の相対誤差（大きさ 1 未満の値は絶対誤差）
constexpr Real GEO_PI = 3.141592653589793238462643383279502884L;

namespace geo_internal {
__extension__ typedef unsigned __int128 UInt128;
template <class T, class... A>
constexpr bool coord_args_ok = (std::is_arithmetic_v<A> && ...) && (std::is_floating_point_v<T> || (std::is_integral_v<A> && ...));
}  // namespace geo_internal

// 実数の符号（|x| ≤ GEO_EPS · max(1, |scale|) なら 0）。scale には比べている量の大きさを渡す
inline int real_sign(Real x, Real scale = 1) {
    const Real e = GEO_EPS * std::max<Real>(1, std::fabs(scale));
    return x > e ? 1 : (x < -e ? -1 : 0);
}
inline int real_cmp(Real a, Real b) { return real_sign(a - b, std::max(std::fabs(a), std::fabs(b))); }  // -1 / 0 / 1
inline bool real_eq(Real a, Real b) { return real_cmp(a, b) == 0; }

template <class T>
struct PointT {
    T x = 0, y = 0;
    constexpr PointT() {}
    template <class A, class B, std::enable_if_t<geo_internal::coord_args_ok<T, A, B>, int> = 0>
    constexpr PointT(A x_, B y_) : x(static_cast<T>(x_)), y(static_cast<T>(y_)) {}
    template <class U, std::enable_if_t<!std::is_same_v<T, U> && std::is_floating_point_v<T>, int> = 0>
    constexpr explicit PointT(const PointT<U>& p) : x(static_cast<T>(p.x)), y(static_cast<T>(p.y)) {}

    constexpr PointT operator+() const { return *this; }
    constexpr PointT operator-() const { return PointT(-x, -y); }
    constexpr PointT& operator+=(const PointT& p) { x += p.x, y += p.y; return *this; }
    constexpr PointT& operator-=(const PointT& p) { x -= p.x, y -= p.y; return *this; }
    constexpr PointT& operator*=(T k) { x *= k, y *= k; return *this; }
    constexpr PointT& operator/=(T k) { x /= k, y /= k; return *this; }  // 整数の点は切り捨て（割り切れるときに使う）
    template <class K, std::enable_if_t<std::is_floating_point_v<K> && std::is_integral_v<T>, int> = 0>
    PointT& operator*=(K) = delete;  // 整数の点に実数は掛けない（to_d(p) を使う）
    template <class K, std::enable_if_t<std::is_floating_point_v<K> && std::is_integral_v<T>, int> = 0>
    PointT& operator/=(K) = delete;

    friend constexpr PointT operator+(PointT a, const PointT& b) { return a += b; }
    friend constexpr PointT operator-(PointT a, const PointT& b) { return a -= b; }
    friend constexpr PointT operator*(PointT a, T k) { return a *= k; }
    friend constexpr PointT operator*(T k, PointT a) { return a *= k; }
    friend constexpr PointT operator/(PointT a, T k) { return a /= k; }
    template <class K>  // 整数の点に実数は掛けない・割らない（コンパイルエラーにする。to_d(p) を使う）
    friend std::enable_if_t<std::is_floating_point_v<K> && std::is_integral_v<T>, PointT> operator*(const PointT&, K) = delete;
    template <class K>
    friend std::enable_if_t<std::is_floating_point_v<K> && std::is_integral_v<T>, PointT> operator*(K, const PointT&) = delete;
    template <class K>
    friend std::enable_if_t<std::is_floating_point_v<K> && std::is_integral_v<T>, PointT> operator/(const PointT&, K) = delete;

    friend constexpr bool operator==(const PointT& a, const PointT& b) { return a.x == b.x && a.y == b.y; }
    friend constexpr bool operator!=(const PointT& a, const PointT& b) { return !(a == b); }
    friend constexpr bool operator<(const PointT& a, const PointT& b) { return a.x != b.x ? a.x < b.x : a.y < b.y; }
    friend constexpr bool operator>(const PointT& a, const PointT& b) { return b < a; }
    friend constexpr bool operator<=(const PointT& a, const PointT& b) { return !(b < a); }
    friend constexpr bool operator>=(const PointT& a, const PointT& b) { return !(a < b); }
};
using Point = PointT<Coord>;
using PointD = PointT<Real>;

template <class T>
std::ostream& operator<<(std::ostream& os, const PointT<T>& p) { return os << "(" << p.x << ", " << p.y << ")"; }
template <class T>
std::istream& operator>>(std::istream& is, PointT<T>& p) { return is >> p.x >> p.y; }

// ---- 内積・外積・長さ ----
template <class T>
constexpr T dot(const PointT<T>& a, const PointT<T>& b) { return a.x * b.x + a.y * b.y; }
template <class T>
constexpr T cross(const PointT<T>& a, const PointT<T>& b) { return a.x * b.y - a.y * b.x; }
template <class T>
constexpr T dot(const PointT<T>& o, const PointT<T>& a, const PointT<T>& b) { return dot(a - o, b - o); }
template <class T>
constexpr T cross(const PointT<T>& o, const PointT<T>& a, const PointT<T>& b) { return cross(a - o, b - o); }
template <class T>
constexpr T norm2(const PointT<T>& a) { return a.x * a.x + a.y * a.y; }
inline Int128 dot128(const Point& a, const Point& b) { return (Int128)a.x * b.x + (Int128)a.y * b.y; }
inline Int128 cross128(const Point& a, const Point& b) { return (Int128)a.x * b.y - (Int128)a.y * b.x; }
inline Int128 norm2_128(const Point& a) { return (Int128)a.x * a.x + (Int128)a.y * a.y; }
template <class T>
Real len(const PointT<T>& a) { return std::sqrt((Real)a.x * (Real)a.x + (Real)a.y * (Real)a.y); }
template <class T>
constexpr T dist2(const PointT<T>& a, const PointT<T>& b) { return norm2(a - b); }
template <class T>
Real dist(const PointT<T>& a, const PointT<T>& b) { return len(a - b); }

// 実数の点どうしの比較（座標の大きさに合わせた相対誤差）
inline bool real_eq(const PointD& a, const PointD& b) {
    const Real s = std::max(std::max(std::fabs(a.x), std::fabs(a.y)), std::max(std::fabs(b.x), std::fabs(b.y)));
    return real_sign(a.x - b.x, s) == 0 && real_sign(a.y - b.y, s) == 0;
}

// ---- 向き ----
enum CcwDetail : int { CLOCKWISE = -1, ON_SEGMENT = 0, COUNTER_CLOCKWISE = 1, ONLINE_FRONT = -2, ONLINE_BACK = 2 };
enum Where : int { OUTSIDE = 0, ON_BORDER = 1, INSIDE = 2 };  // 内外判定の結果

// +1: a → b → c で左に曲がる（反時計回り）/ -1: 右 / 0: 一直線。整数なら厳密
template <class T>
int ccw(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c) {
    const PointT<T> u = b - a, v = c - a;
    if constexpr (std::is_integral_v<T>) {
        const T s = cross(u, v);
        return (s > 0) - (s < 0);
    } else {
        return real_sign(cross(u, v), len(u) * len(v));
    }
}
// 一直線のときの位置まで返す（a ≠ b が前提。a == b なら c == a で ON_SEGMENT、それ以外は ONLINE_FRONT）
template <class T>
int ccw_detail(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c) {
    const int s = ccw(a, b, c);
    if (s != 0) return s;
    const PointT<T> u = b - a, v = c - a;
    if constexpr (std::is_integral_v<T>) {
        if (dot(u, v) < 0) return ONLINE_BACK;
        if (norm2(u) < norm2(v)) return ONLINE_FRONT;
    } else {
        const Real nu = norm2(u), nv = norm2(v);
        if (real_sign(dot(u, v), std::sqrt(nu * nv)) < 0) return ONLINE_BACK;
        if (real_sign(nv - nu, std::max(nu, nv)) > 0) return ONLINE_FRONT;
    }
    return ON_SEGMENT;
}
template <class T>
bool collinear(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c) { return ccw(a, b, c) == 0; }
namespace geo_internal {
template <class T>
int sgn_cross(const PointT<T>& u, const PointT<T>& v) {  // 外積の符号（実数なら |u||v| に合わせた相対誤差）
    if constexpr (std::is_integral_v<T>) {
        const T c = cross(u, v);
        return (c > 0) - (c < 0);
    } else {
        return real_sign(cross(u, v), len(u) * len(v));
    }
}
template <class T>
int sgn_dot(const PointT<T>& u, const PointT<T>& v) {  // 内積の符号（同上）
    if constexpr (std::is_integral_v<T>) {
        const T c = dot(u, v);
        return (c > 0) - (c < 0);
    } else {
        return real_sign(dot(u, v), len(u) * len(v));
    }
}
}  // namespace geo_internal

// ---- 回転 ----
template <class T>
constexpr PointT<T> rot90(const PointT<T>& p) { return PointT<T>(-p.y, p.x); }  // 反時計回りに 90 度
template <class T>
constexpr PointT<T> rot90_cw(const PointT<T>& p) { return PointT<T>(p.y, -p.x); }  // 時計回りに 90 度
template <class T>
constexpr PointT<T> rot180(const PointT<T>& p) { return -p; }
template <class T>
constexpr PointT<T> rot90(const PointT<T>& p, const PointT<T>& c, long long k) {  // c のまわりに反時計回りに 90 度を k 回
    PointT<T> d = p - c;
    const long long r = (k % 4 + 4) % 4;
    if (r == 1) d = rot90(d);
    if (r == 2) d = -d;
    if (r == 3) d = rot90_cw(d);
    return d + c;
}
template <class T>
constexpr PointD to_d(const PointT<T>& p) { return PointD(p.x, p.y); }
template <class T>
PointD rotated(const PointT<T>& p, Real theta) {  // 原点のまわりに反時計回りに theta（ラジアン）
    const Real c = std::cos(theta), s = std::sin(theta);
    const PointD q = to_d(p);
    return PointD(q.x * c - q.y * s + (Real)0, q.x * s + q.y * c + (Real)0);  // + 0 で -0 を 0 に
}
// center のまわりに。p と center の型は違ってよい（Point と PointD など）
template <class T, class U>
PointD rotated(const PointT<T>& p, Real theta, const PointT<U>& center) {
    const PointD r = rotated(to_d(p) - to_d(center), theta) + to_d(center);
    return PointD(r.x + (Real)0, r.y + (Real)0);
}
template <class T>
PointD rotated_deg(const PointT<T>& p, Real deg) {  // 度で。90 度の倍数なら誤差なし
    const Real q = deg / 90;
    if (q == std::floor(q) && std::fabs(q) < 1e18L) {
        PointT<T> r = p;
        const long long k = ((long long)q % 4 + 4) % 4;
        for (long long i = 0; i < k; ++i) r = rot90(r);
        return PointD((Real)r.x + (Real)0, (Real)r.y + (Real)0);
    }
    return rotated(p, deg * GEO_PI / 180);
}
template <class T, class U>
PointD rotated_deg(const PointT<T>& p, Real deg, const PointT<U>& center) {
    const PointD r = rotated_deg(to_d(p) - to_d(center), deg) + to_d(center);
    return PointD(r.x + (Real)0, r.y + (Real)0);
}
template <class T>
constexpr PointT<T> complex_mul(const PointT<T>& a, const PointT<T>& b) { return PointT<T>(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
template <class T>
constexpr PointT<T> conjugate(const PointT<T>& a) { return PointT<T>(a.x, -a.y); }

// ---- アフィン変換 f(p) = (a x + b y + c, d x + e y + f)。回転・折り返し・平行移動を合成して 1 つにまとめる ----
// 整数の Affine なら 90 度の回転・折り返し・平行移動・整数倍の合成は整数のまま（誤差なし）。任意の角は AffineD::rotate
template <class T>
struct AffineT {
    T a = 1, b = 0, c = 0, d = 0, e = 1, f = 0;  // 既定は恒等変換（何もしない）
    constexpr AffineT() {}
    constexpr AffineT(T a_, T b_, T c_, T d_, T e_, T f_) : a(a_), b(b_), c(c_), d(d_), e(e_), f(f_) {}
    constexpr PointT<T> operator()(const PointT<T>& p) const { return PointT<T>(a * p.x + b * p.y + c, d * p.x + e * p.y + f); }
    constexpr AffineT then(const AffineT& g) const {  // この変換をしてから g をする変換（g ∘ this）
        return AffineT(g.a * a + g.b * d, g.a * b + g.b * e, g.a * c + g.b * f + g.c, g.d * a + g.e * d, g.d * b + g.e * e,
                       g.d * c + g.e * f + g.f);
    }
    static constexpr AffineT translate(const PointT<T>& v) { return AffineT(1, 0, v.x, 0, 1, v.y); }  // v だけ平行移動
    static constexpr AffineT rot90(const PointT<T>& o = PointT<T>()) { return AffineT(0, -1, o.x + o.y, 1, 0, o.y - o.x); }  // o のまわりに反時計回り 90 度
    static constexpr AffineT rot90_cw(const PointT<T>& o = PointT<T>()) { return AffineT(0, 1, o.x - o.y, -1, 0, o.y + o.x); }  // 時計回り 90 度
    static constexpr AffineT rot180(const PointT<T>& o = PointT<T>()) { return AffineT(-1, 0, 2 * o.x, 0, -1, 2 * o.y); }
    static constexpr AffineT reflect_x(T p) { return AffineT(-1, 0, 2 * p, 0, 1, 0); }  // 直線 x = p について折り返す
    static constexpr AffineT reflect_y(T p) { return AffineT(1, 0, 0, 0, -1, 2 * p); }  // 直線 y = p について折り返す
    static constexpr AffineT scale(T k, const PointT<T>& o = PointT<T>()) { return AffineT(k, 0, o.x - k * o.x, 0, k, o.y - k * o.y); }  // o を中心に k 倍
    static AffineT rotate(Real theta, const PointT<T>& o = PointT<T>()) {  // o のまわりに反時計回りに theta（AffineD だけ）
        static_assert(std::is_floating_point_v<T>, "rotate は AffineD で使う（整数の Affine は rot90 / rot90_cw / rot180）");
        const T cs = std::cos(theta), sn = std::sin(theta);
        return AffineT(cs, -sn, o.x - cs * o.x + sn * o.y, sn, cs, o.y - sn * o.x - cs * o.y);
    }
    friend constexpr bool operator==(const AffineT& p, const AffineT& q) {
        return p.a == q.a && p.b == q.b && p.c == q.c && p.d == q.d && p.e == q.e && p.f == q.f;
    }
    friend constexpr bool operator!=(const AffineT& p, const AffineT& q) { return !(p == q); }
};
using Affine = AffineT<Coord>;
using AffineD = AffineT<Real>;

// ---- 実数の点・角度 ----
template <class T>
PointD midpoint(const PointT<T>& a, const PointT<T>& b) { return PointD(((Real)a.x + (Real)b.x) / 2, ((Real)a.y + (Real)b.y) / 2); }
template <class T>
PointD unit(const PointT<T>& v) {  // 長さ 1 に（0 ベクトルは 0 のまま）
    const Real l = len(v);
    return l == 0 ? PointD() : PointD((Real)v.x / l, (Real)v.y / l);
}
inline PointD from_polar(Real r, Real theta) { return PointD(r * std::cos(theta), r * std::sin(theta)); }
inline Real deg_to_rad(Real d) { return d * GEO_PI / 180; }
inline Real rad_to_deg(Real r) { return r * 180 / GEO_PI; }
template <class T>
Real angle_of(const PointT<T>& v) {  // 偏角 [0, 2π)。0 ベクトルは 0
    Real t = std::atan2((Real)v.y, (Real)v.x);
    if (t < 0) t += 2 * GEO_PI;
    if (t >= 2 * GEO_PI || t == 0) t = 0;  // -0 や、2π に丸められた値を 0 に
    return t;
}
template <class T>
Real angle_of_deg(const PointT<T>& v) { return rad_to_deg(angle_of(v)); }
template <class T>
Real angle_between(const PointT<T>& a, const PointT<T>& b) {  // なす角 [0, π]（どちらかが 0 ベクトルなら 0）
    if (a == PointT<T>() || b == PointT<T>()) return 0;  // atan2(0, -0) = π にならないように
    const Real c = (Real)a.x * (Real)b.y - (Real)a.y * (Real)b.x, d = (Real)a.x * (Real)b.x + (Real)a.y * (Real)b.y;
    return std::atan2(std::fabs(c), d);
}
template <class T>
Real angle_between_deg(const PointT<T>& a, const PointT<T>& b) { return rad_to_deg(angle_between(a, b)); }

// ---- 向きの正規化（同じ向き・同じ傾きをまとめるキー） ----
inline Point direction_key(const Point& v) {  // gcd で割る。向きを区別する（(2,4) と (-2,-4) は別）。0 ベクトルは (0, 0)
    const Coord g = std::gcd(v.x < 0 ? -v.x : v.x, v.y < 0 ? -v.y : v.y);
    return g == 0 ? Point() : Point(v.x / g, v.y / g);
}
inline Point line_direction_key(const Point& v) {  // 向きを区別しない（x > 0、または x == 0 かつ y > 0 にそろえる）
    Point d = direction_key(v);
    if (d.x < 0 || (d.x == 0 && d.y < 0)) d = -d;
    return d;
}

// ---- 入力 ----
inline Point read_point() { Point p; std::cin >> p.x >> p.y; return p; }
inline PointD read_point_d() { PointD p; std::cin >> p.x >> p.y; return p; }
inline Point read_point_fixed(int k) { const long long x = read_fixed(k); return Point(x, read_fixed(k)); }
inline std::vector<Point> read_points(int n) {
    std::vector<Point> ps(n);
    for (auto& p : ps) std::cin >> p.x >> p.y;
    return ps;
}
inline std::vector<PointD> read_points_d(int n) {
    std::vector<PointD> ps(n);
    for (auto& p : ps) std::cin >> p.x >> p.y;
    return ps;
}
inline std::vector<Point> read_points_fixed(int n, int k) {
    std::vector<Point> ps(n);
    for (auto& p : ps) p = read_point_fixed(k);
    return ps;
}

struct PointHash {  // unordered_set<Point, PointHash> / unordered_map<Point, V, PointHash>
    static unsigned long long mix(unsigned long long z) {
        z += 0x9e3779b97f4a7c15ULL;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    size_t operator()(const Point& p) const { return (size_t)mix(mix((unsigned long long)p.x) ^ (unsigned long long)p.y); }
};

// ---- 直線・線分・円の型（アルゴリズムは 02_line / 05_circle など） ----
template <class T>
struct LineT {  // 直線 a → b（a ≠ b。向きを持つ：左側 = cross(b - a, p - a) > 0）
    PointT<T> a, b;
    constexpr LineT() {}
    constexpr LineT(const PointT<T>& a_, const PointT<T>& b_) : a(a_), b(b_) {}
    template <class A, class B, class C, class D, std::enable_if_t<geo_internal::coord_args_ok<T, A, B, C, D>, int> = 0>
    constexpr LineT(A x1, B y1, C x2, D y2) : a(x1, y1), b(x2, y2) {}
    constexpr PointT<T> dir() const { return b - a; }
    friend constexpr bool operator==(const LineT& l, const LineT& m) { return l.a == m.a && l.b == m.b; }
    friend constexpr bool operator!=(const LineT& l, const LineT& m) { return !(l == m); }
};
template <class T>
struct SegmentT {  // 線分 a - b（a == b の点も可）
    PointT<T> a, b;
    constexpr SegmentT() {}
    constexpr SegmentT(const PointT<T>& a_, const PointT<T>& b_) : a(a_), b(b_) {}
    template <class A, class B, class C, class D, std::enable_if_t<geo_internal::coord_args_ok<T, A, B, C, D>, int> = 0>
    constexpr SegmentT(A x1, B y1, C x2, D y2) : a(x1, y1), b(x2, y2) {}
    constexpr PointT<T> dir() const { return b - a; }
    friend constexpr bool operator==(const SegmentT& s, const SegmentT& t) { return s.a == t.a && s.b == t.b; }
    friend constexpr bool operator!=(const SegmentT& s, const SegmentT& t) { return !(s == t); }
};
template <class T>
struct CircleT {  // 円（中心 c・半径 r）
    PointT<T> c;
    T r = 0;
    constexpr CircleT() {}
    template <class R, std::enable_if_t<geo_internal::coord_args_ok<T, R>, int> = 0>
    constexpr CircleT(const PointT<T>& c_, R r_) : c(c_), r(static_cast<T>(r_)) {}
    template <class A, class B, class R, std::enable_if_t<geo_internal::coord_args_ok<T, A, B, R>, int> = 0>
    constexpr CircleT(A x, B y, R r_) : c(x, y), r(static_cast<T>(r_)) {}
    friend constexpr bool operator==(const CircleT& p, const CircleT& q) { return p.c == q.c && p.r == q.r; }
    friend constexpr bool operator!=(const CircleT& p, const CircleT& q) { return !(p == q); }
};
using Line = LineT<Coord>;
using LineD = LineT<Real>;
using Segment = SegmentT<Coord>;
using SegmentD = SegmentT<Real>;
using Circle = CircleT<Coord>;
using CircleD = CircleT<Real>;

template <class T>
constexpr LineD to_d(const LineT<T>& l) { return LineD(to_d(l.a), to_d(l.b)); }
template <class T>
constexpr SegmentD to_d(const SegmentT<T>& s) { return SegmentD(to_d(s.a), to_d(s.b)); }
template <class T>
constexpr CircleD to_d(const CircleT<T>& c) { return CircleD(to_d(c.c), c.r); }
template <class T>
std::ostream& operator<<(std::ostream& os, const LineT<T>& l) { return os << "Line[" << l.a << " -> " << l.b << "]"; }
template <class T>
std::ostream& operator<<(std::ostream& os, const SegmentT<T>& s) { return os << "Segment[" << s.a << " - " << s.b << "]"; }
template <class T>
std::ostream& operator<<(std::ostream& os, const CircleT<T>& c) { return os << "Circle[c=" << c.c << ", r=" << c.r << "]"; }
inline Segment read_segment() { const Point a = read_point(); return Segment(a, read_point()); }
inline Line read_line() { const Point a = read_point(); return Line(a, read_point()); }
inline Circle read_circle() { const Point c = read_point(); Coord r; std::cin >> r; return Circle(c, r); }
#endif  // GG_GEO_BASE
#ifndef GG_LINE
#define GG_LINE
// ---- 直線・線分（幾何第二世代）  判定は整数で厳密、交点・距離・射影は実数（PointD / Real） ----
enum LineRelation : int { LINE_CROSS = 0, LINE_PARALLEL = 1, LINE_SAME = 2 };
// 線分どうしの位置関係：共有点なし / 1 点で接する（端点・T 字・一直線で端どうし）/ 真に交差 / 一直線で長さのある部分が重なる
enum SegmentRelation : int { SEGMENT_SEPARATE = 0, SEGMENT_TOUCH = 1, SEGMENT_CROSS = 2, SEGMENT_OVERLAP = 3 };
template <class T>
struct RayT {  // 半直線：a から b の向きへどこまでも（a ≠ b）。a は含む
    PointT<T> a, b;
    constexpr RayT() {}
    constexpr RayT(const PointT<T>& a_, const PointT<T>& b_) : a(a_), b(b_) {}
    template <class A, class B, class C, class D, std::enable_if_t<geo_internal::coord_args_ok<T, A, B, C, D>, int> = 0>
    constexpr RayT(A x1, B y1, C x2, D y2) : a(x1, y1), b(x2, y2) {}
    constexpr PointT<T> dir() const { return b - a; }
    friend constexpr bool operator==(const RayT& r, const RayT& q) { return r.a == q.a && r.b == q.b; }
    friend constexpr bool operator!=(const RayT& r, const RayT& q) { return !(r == q); }
};
using Ray = RayT<Coord>;
using RayD = RayT<Real>;
template <class T>
constexpr RayD to_d(const RayT<T>& r) { return RayD(to_d(r.a), to_d(r.b)); }
template <class T>
std::ostream& operator<<(std::ostream& os, const RayT<T>& r) { return os << "Ray[" << r.a << " -> " << r.b << " ...]"; }
namespace geo_internal {
template <class X>
struct linear_t : std::false_type {};
template <class T>
struct linear_t<LineT<T>> : std::true_type {};
template <class T>
struct linear_t<SegmentT<T>> : std::true_type {};
template <class T>
struct linear_t<RayT<T>> : std::true_type {};
template <class X>
constexpr bool is_linear_v = linear_t<X>::value;
// base + dir * (num / den)。整数なら分子 base * den + dir * num を 128 ビットで正確に求めてから割る
template <class T>
PointD affine(const PointT<T>& base, const PointT<T>& dir, Int128 num, Int128 den) {
    if constexpr (std::is_integral_v<T>) {
        const Int128 x = (Int128)base.x * den + (Int128)dir.x * num, y = (Int128)base.y * den + (Int128)dir.y * num;
        return PointD((Real)x / (Real)den, (Real)y / (Real)den);
    } else {
        return base + dir * ((Real)num / (Real)den);
    }
}
inline Int128 abs128(Int128 v) { return v < 0 ? -v : v; }
inline Int128 gcd128(Int128 a, Int128 b) {
    a = abs128(a), b = abs128(b);
    while (b != 0) {
        const Int128 t = a % b;
        a = b, b = t;
    }
    return a;
}
}  // namespace geo_internal

// ---- 厳密な判定（整数の型なら誤差なし） ----
template <class L1, class L2, std::enable_if_t<geo_internal::is_linear_v<L1> && geo_internal::is_linear_v<L2>, int> = 0>
bool is_parallel(const L1& l, const L2& m) { return geo_internal::sgn_cross(l.dir(), m.dir()) == 0; }
template <class L1, class L2, std::enable_if_t<geo_internal::is_linear_v<L1> && geo_internal::is_linear_v<L2>, int> = 0>
bool is_orthogonal(const L1& l, const L2& m) { return geo_internal::sgn_dot(l.dir(), m.dir()) == 0; }
template <class T>
bool on_line(const LineT<T>& l, const PointT<T>& p) { return ccw(l.a, l.b, p) == 0; }
template <class T>
bool on_segment(const SegmentT<T>& s, const PointT<T>& p) { return ccw_detail(s.a, s.b, p) == ON_SEGMENT; }
template <class T>
int side(const LineT<T>& l, const PointT<T>& p) { return ccw(l.a, l.b, p); }  // 左 +1 / 右 -1 / 直線の上 0
// 線分どうしが共有点を持つか（端点で触れる・一直線で重なるも含む）
template <class T>
bool intersects(const SegmentT<T>& s, const SegmentT<T>& t) {
    if (s.a == s.b) return on_segment(t, s.a);
    if (t.a == t.b) return on_segment(s, t.a);
    return ccw_detail(s.a, s.b, t.a) * ccw_detail(s.a, s.b, t.b) <= 0 && ccw_detail(t.a, t.b, s.a) * ccw_detail(t.a, t.b, s.b) <= 0;
}
// 真に交差（1 点で交わり、その点がどちらの線分の端点でもない）
template <class T>
bool intersects_strict(const SegmentT<T>& s, const SegmentT<T>& t) {
    return ccw(s.a, s.b, t.a) * ccw(s.a, s.b, t.b) < 0 && ccw(t.a, t.b, s.a) * ccw(t.a, t.b, s.b) < 0;
}
template <class T>
bool intersects(const LineT<T>& l, const SegmentT<T>& s) { return ccw(l.a, l.b, s.a) * ccw(l.a, l.b, s.b) <= 0; }
template <class T>
bool intersects(const SegmentT<T>& s, const LineT<T>& l) { return intersects(l, s); }
template <class T>
int line_relation(const LineT<T>& l, const LineT<T>& m) {
    if (!is_parallel(l, m)) return LINE_CROSS;
    return on_line(l, m.a) ? LINE_SAME : LINE_PARALLEL;
}
template <class T>
bool intersects(const LineT<T>& l, const LineT<T>& m) { return line_relation(l, m) != LINE_PARALLEL; }

// ---- 射影・距離（実数） ----
template <class T>
PointD projection(const LineT<T>& l, const PointT<T>& p) {  // 垂線の足
    const PointT<T> d = l.dir(), w = p - l.a;
    if constexpr (std::is_integral_v<T>) return geo_internal::affine(l.a, d, dot128(d, w), norm2_128(d));
    else return l.a + d * (dot(d, w) / norm2(d));
}
template <class T>
PointD reflection(const LineT<T>& l, const PointT<T>& p) {  // 直線について線対称な点
    if constexpr (std::is_integral_v<T>) {
        const Point d = l.dir(), w = p - l.a;
        const Int128 den = norm2_128(d), num = dot128(d, w);
        const Int128 x = 2 * ((Int128)l.a.x * den + (Int128)d.x * num) - (Int128)p.x * den;
        const Int128 y = 2 * ((Int128)l.a.y * den + (Int128)d.y * num) - (Int128)p.y * den;
        return PointD((Real)x / (Real)den, (Real)y / (Real)den);
    } else {
        return projection(l, p) * (Real)2 - p;
    }
}
template <class T>
Real dist(const LineT<T>& l, const PointT<T>& p) {
    const PointT<T> d = l.dir();
    if constexpr (std::is_integral_v<T>) {
        const T c = cross(d, p - l.a);
        return (Real)(c < 0 ? -c : c) / len(d);
    } else {
        return std::fabs(cross(d, p - l.a)) / len(d);
    }
}
template <class T>
Real dist(const SegmentT<T>& s, const PointT<T>& p) {
    const PointT<T> d = s.dir();
    if (geo_internal::sgn_dot(d, p - s.a) <= 0) return dist(s.a, p);
    if (geo_internal::sgn_dot(-d, p - s.b) <= 0) return dist(s.b, p);
    return dist(LineT<T>(s.a, s.b), p);
}
template <class T>
Real dist(const PointT<T>& p, const LineT<T>& l) { return dist(l, p); }
template <class T>
Real dist(const PointT<T>& p, const SegmentT<T>& s) { return dist(s, p); }
template <class T>
Real dist(const SegmentT<T>& s, const SegmentT<T>& t) {
    if (intersects(s, t)) return 0;
    return std::min(std::min(dist(s, t.a), dist(s, t.b)), std::min(dist(t, s.a), dist(t, s.b)));
}
template <class T>
Real dist(const LineT<T>& l, const SegmentT<T>& s) { return intersects(l, s) ? 0 : std::min(dist(l, s.a), dist(l, s.b)); }
template <class T>
Real dist(const SegmentT<T>& s, const LineT<T>& l) { return dist(l, s); }
template <class T>
Real dist(const LineT<T>& l, const LineT<T>& m) { return is_parallel(l, m) ? dist(l, m.a) : 0; }
// 距離の 2 乗と整数 D2 を誤差なしで比べる（-1: 距離² < D2 / 0: 等しい / 1: 大きい）
inline int dist2_cmp(const Line& l, const Point& p, Coord D2) {
    const Point d = l.dir();
    const Int128 c = cross128(d, p - l.a), lhs = c * c, rhs = (Int128)D2 * norm2_128(d);
    return (lhs > rhs) - (lhs < rhs);
}
inline int dist2_cmp(const Segment& s, const Point& p, Coord D2) {
    const Point d = s.dir();
    const bool near_a = dot128(d, p - s.a) <= 0, near_b = !near_a && dot128(-d, p - s.b) <= 0;
    if (!near_a && !near_b) return dist2_cmp(Line(s.a, s.b), p, D2);
    const Int128 v = norm2_128(p - (near_a ? s.a : s.b));
    return (v > D2) - (v < D2);
}

// ---- 交点 ----
template <class T>
std::vector<PointD> cross_point(const LineT<T>& l, const LineT<T>& m) {  // 直線の交点（平行・同じ直線なら空）
    const PointT<T> d1 = l.dir(), d2 = m.dir();
    if (geo_internal::sgn_cross(d1, d2) == 0) return {};
    if constexpr (std::is_integral_v<T>) return {geo_internal::affine(l.a, d1, (Int128)cross(m.a - l.a, d2), (Int128)cross(d1, d2))};
    else return {l.a + d1 * (cross(m.a - l.a, d2) / cross(d1, d2))};
}
template <class T>
std::vector<PointD> cross_point(const LineT<T>& l, const SegmentT<T>& s) {  // 直線と線分（線分が直線の上なら両端）
    if (!intersects(l, s)) return {};
    if (s.a == s.b) return {to_d(s.a)};
    if (geo_internal::sgn_cross(l.dir(), s.dir()) == 0) return {to_d(s.a), to_d(s.b)};
    return cross_point(l, LineT<T>(s.a, s.b));
}
template <class T>
std::vector<PointD> cross_point(const SegmentT<T>& s, const LineT<T>& l) { return cross_point(l, s); }
// 線分の交点。重なるときは重なりの両端（s の向きの順）、1 点で重なるなら 1 個
template <class T>
std::vector<PointD> cross_points(const SegmentT<T>& s, const SegmentT<T>& t) {
    if (!intersects(s, t)) return {};
    if (s.a == s.b) return {to_d(s.a)};
    if (t.a == t.b) return {to_d(t.a)};
    const PointT<T> d = s.dir();
    if (geo_internal::sgn_cross(d, t.dir()) != 0) return cross_point(LineT<T>(s.a, s.b), LineT<T>(t.a, t.b));
    auto key = [&](const PointT<T>& p) { return dot(d, p - s.a); };  // s の向きに沿った位置
    PointT<T> lo = t.a, hi = t.b;
    if (key(lo) > key(hi)) std::swap(lo, hi);
    const PointT<T> L = key(lo) > key(s.a) ? lo : s.a, R = key(hi) < key(s.b) ? hi : s.b;
    bool same;
    if constexpr (std::is_integral_v<T>) same = (L == R);
    else same = real_eq(L, R);
    if (same) return {to_d(L)};
    return {to_d(L), to_d(R)};
}

// 線分どうしの位置関係（SEGMENT_SEPARATE / SEGMENT_TOUCH / SEGMENT_CROSS / SEGMENT_OVERLAP。整数なら厳密）
template <class T>
int segment_relation(const SegmentT<T>& s, const SegmentT<T>& t) {
    if (!intersects(s, t)) return SEGMENT_SEPARATE;
    if (intersects_strict(s, t)) return SEGMENT_CROSS;
    return cross_points(s, t).size() == 2 ? SEGMENT_OVERLAP : SEGMENT_TOUCH;
}

// 有理数の点 (x / d, y / d)（d > 0、既約）。同じ交点かどうかを誤差なしで判定するのに使う
struct FracPoint {
    Int128 x = 0, y = 0, d = 1;
    PointD point_d() const { return PointD((Real)x / (Real)d, (Real)y / (Real)d); }
    bool is_integer() const { return d == 1; }
    friend bool operator==(const FracPoint& p, const FracPoint& q) { return p.x == q.x && p.y == q.y && p.d == q.d; }
    friend bool operator!=(const FracPoint& p, const FracPoint& q) { return !(p == q); }
    friend bool operator<(const FracPoint& p, const FracPoint& q) {  // set / map 用の順（大きさの順ではない）
        if (p.x != q.x) return p.x < q.x;
        if (p.y != q.y) return p.y < q.y;
        return p.d < q.d;
    }
};
struct FracPointHash {
    size_t operator()(const FracPoint& p) const {
        const PointHash h;
        const size_t a = h(Point((Coord)(p.x >> 64) ^ (Coord)p.x, (Coord)(p.y >> 64) ^ (Coord)p.y));
        return a ^ (size_t)PointHash::mix((unsigned long long)((Coord)(p.d >> 64) ^ (Coord)p.d));
    }
};
inline std::vector<FracPoint> cross_point_frac(const Line& l, const Line& m) {  // 直線の交点を既約分数で（平行・同じ直線なら空）
    const Point d1 = l.dir(), d2 = m.dir();
    Int128 den = cross128(d1, d2), num = cross128(m.a - l.a, d2);
    if (den == 0) return {};
    if (den < 0) den = -den, num = -num;
    Int128 x = (Int128)l.a.x * den + (Int128)d1.x * num, y = (Int128)l.a.y * den + (Int128)d1.y * num;
    const Int128 g = geo_internal::gcd128(geo_internal::gcd128(x, y), den);
    return {FracPoint{x / g, y / g, den / g}};
}

// ---- 直線を作る ----
template <class T>
LineT<T> perpendicular(const LineT<T>& l, const PointT<T>& p) { return LineT<T>(p, p + rot90(l.dir())); }  // p を通る垂線
template <class T>
LineT<T> parallel_through(const LineT<T>& l, const PointT<T>& p) { return LineT<T>(p, p + l.dir()); }  // p を通る平行線
template <class T>
LineD bisector(const PointT<T>& a, const PointT<T>& b) {  // 垂直二等分線（a ≠ b）
    const PointD m = midpoint(a, b);
    return LineD(m, m + to_d(rot90(b - a)));
}

// 直線 a x + b y = c を、(a, b, c) を gcd で割って a > 0（a == 0 なら b > 0）にそろえたキー。同じ直線なら同じキー
struct LineKey {
    Coord a = 0, b = 0, c = 0;
    LineKey() {}
    LineKey(Coord a_, Coord b_, Coord c_) {  // (a, b) ≠ (0, 0)
        const Coord g = std::gcd(std::gcd(a_ < 0 ? -a_ : a_, b_ < 0 ? -b_ : b_), c_ < 0 ? -c_ : c_);
        a = a_ / g, b = b_ / g, c = c_ / g;
        if (a < 0 || (a == 0 && b < 0)) a = -a, b = -b, c = -c;
    }
    bool contains(const Point& p) const { return (Int128)a * p.x + (Int128)b * p.y == (Int128)c; }
    Point dir_key() const { return line_direction_key(Point(b, -a)); }  // 傾きのキー（平行なら同じ。line_direction_key と同じ形）
    LineD line() const {                                                // 実数の直線に
        const Real n = (Real)a * a + (Real)b * b;
        const PointD p((Real)a * c / n, (Real)b * c / n);
        return LineD(p, p + PointD(-b, a));
    }
    friend bool operator==(const LineKey& p, const LineKey& q) { return p.a == q.a && p.b == q.b && p.c == q.c; }
    friend bool operator!=(const LineKey& p, const LineKey& q) { return !(p == q); }
    friend bool operator<(const LineKey& p, const LineKey& q) {
        if (p.a != q.a) return p.a < q.a;
        if (p.b != q.b) return p.b < q.b;
        return p.c < q.c;
    }
};
struct LineKeyHash {
    size_t operator()(const LineKey& k) const { return PointHash()(Point(k.a, k.b)) ^ (size_t)PointHash::mix((unsigned long long)k.c + 0x51ed27ULL); }
};
inline LineKey line_key(const Point& p, const Point& q) {  // p と q（p ≠ q）を通る直線
    const Point d = q - p;
    return LineKey(d.y, -d.x, (Coord)((Int128)d.y * p.x - (Int128)d.x * p.y));
}
inline LineKey line_key(const Line& l) { return line_key(l.a, l.b); }
inline LineKey bisector_key(const Point& p, const Point& q) {  // p と q（p ≠ q）の垂直二等分線
    return LineKey(2 * (q.x - p.x), 2 * (q.y - p.y), (Coord)(norm2_128(q) - norm2_128(p)));
}

// ---- 半直線 Ray r(a, b)（a から b の向きへ。a を含む）。判定は整数なら厳密 ----
template <class T>
bool on_ray(const RayT<T>& r, const PointT<T>& p) {  // 半直線の上か（始点 a を含む）
    const int c = ccw_detail(r.a, r.b, p);
    return c == ON_SEGMENT || c == ONLINE_FRONT;
}
template <class T>
bool intersects(const RayT<T>& r, const LineT<T>& l) {  // 半直線と直線が共有点を持つか
    const int s = ccw(l.a, l.b, r.a);  // 始点が直線のどちら側か
    return s == 0 || s * geo_internal::sgn_cross(l.dir(), r.dir()) < 0;
}
template <class T>
bool intersects(const LineT<T>& l, const RayT<T>& r) { return intersects(r, l); }
template <class T>
bool intersects(const RayT<T>& r, const SegmentT<T>& s) {  // 半直線と線分が共有点を持つか（端点で触れる・重なるも含む）
    if (s.a == s.b) return on_ray(r, s.a);
    const int c1 = ccw(r.a, r.b, s.a), c2 = ccw(r.a, r.b, s.b);
    if (c1 * c2 > 0) return false;
    if (c1 == 0 && c2 == 0) return on_ray(r, s.a) || on_ray(r, s.b);  // 一直線
    // 線分は半直線を含む直線と 1 点で交わる。その点が始点より前か（t ≥ 0）
    return ccw(s.a, s.b, r.a) * geo_internal::sgn_cross(r.dir(), s.dir()) >= 0;
}
template <class T>
bool intersects(const SegmentT<T>& s, const RayT<T>& r) { return intersects(r, s); }
template <class T>
bool intersects(const RayT<T>& r, const RayT<T>& q) {  // 半直線どうしが共有点を持つか
    if (geo_internal::sgn_cross(r.dir(), q.dir()) == 0) {
        if (ccw(r.a, r.b, q.a) != 0) return false;  // 平行で別の直線
        return on_ray(r, q.a) || on_ray(q, r.a);
    }
    return intersects(r, LineT<T>(q.a, q.b)) && intersects(q, LineT<T>(r.a, r.b));
}
// 交点がちょうど 1 つのときだけ返す（共有点なし・一直線で重なるときは空。重なりは intersects で分かる）
template <class T>
std::vector<PointD> cross_point(const RayT<T>& r, const LineT<T>& l) {
    if (geo_internal::sgn_cross(r.dir(), l.dir()) == 0 || !intersects(r, l)) return {};
    return cross_point(LineT<T>(r.a, r.b), l);
}
template <class T>
std::vector<PointD> cross_point(const LineT<T>& l, const RayT<T>& r) { return cross_point(r, l); }
template <class T>
std::vector<PointD> cross_point(const RayT<T>& r, const RayT<T>& q) {
    if (!intersects(r, q)) return {};
    if (geo_internal::sgn_cross(r.dir(), q.dir()) != 0) return cross_point(LineT<T>(r.a, r.b), LineT<T>(q.a, q.b));
    if (r.a == q.a && geo_internal::sgn_dot(r.dir(), q.dir()) < 0) return {to_d(r.a)};  // 逆向きで始点だけ共有
    return {};
}
// 半直線と線分の交点（一直線で重なるときは重なりの両端を半直線の向きの順に。1 点なら 1 個）
template <class T>
std::vector<PointD> cross_points(const RayT<T>& r, const SegmentT<T>& s) {
    if (!intersects(r, s)) return {};
    if (s.a == s.b) return {to_d(s.a)};
    const PointT<T> d = r.dir();
    if (geo_internal::sgn_cross(d, s.dir()) != 0) return cross_point(LineT<T>(r.a, r.b), LineT<T>(s.a, s.b));
    auto key = [&](const PointT<T>& p) { return dot(d, p - r.a); };  // 半直線の向きに沿った位置
    PointT<T> lo = s.a, hi = s.b;
    if (key(lo) > key(hi)) std::swap(lo, hi);
    const PointT<T> L = on_ray(r, lo) ? lo : r.a;
    bool same;
    if constexpr (std::is_integral_v<T>) same = (L == hi);
    else same = real_eq(L, hi);
    if (same) return {to_d(L)};
    return {to_d(L), to_d(hi)};
}
template <class T>
std::vector<PointD> cross_points(const SegmentT<T>& s, const RayT<T>& r) { return cross_points(r, s); }
template <class T>
Real dist(const RayT<T>& r, const PointT<T>& p) {
    if (geo_internal::sgn_dot(r.dir(), p - r.a) <= 0) return dist(r.a, p);
    return dist(LineT<T>(r.a, r.b), p);
}
template <class T>
Real dist(const PointT<T>& p, const RayT<T>& r) { return dist(r, p); }
template <class T>
Real dist(const RayT<T>& r, const SegmentT<T>& s) {
    if (intersects(r, s)) return 0;
    return std::min({dist(r, s.a), dist(r, s.b), dist(s, r.a)});
}
template <class T>
Real dist(const SegmentT<T>& s, const RayT<T>& r) { return dist(r, s); }
template <class T>
Real dist(const RayT<T>& r, const RayT<T>& q) { return intersects(r, q) ? 0 : std::min(dist(r, q.a), dist(q, r.a)); }
template <class T>
Real dist(const RayT<T>& r, const LineT<T>& l) { return intersects(r, l) ? 0 : dist(l, r.a); }
template <class T>
Real dist(const LineT<T>& l, const RayT<T>& r) { return dist(r, l); }
inline int dist2_cmp(const Ray& r, const Point& p, Coord D2) {  // 距離² と整数 D2 の比較（-1 / 0 / 1。誤差なし）
    if (dot128(r.dir(), p - r.a) <= 0) {
        const Int128 v = norm2_128(p - r.a);
        return (v > D2) - (v < D2);
    }
    return dist2_cmp(Line(r.a, r.b), p, D2);
}

// ---- 格子点 ----
inline Coord lattice_points_on(const Segment& s) {  // 線分の上（両端を含む）の格子点の数
    const Point d = s.dir();
    return std::gcd(d.x < 0 ? -d.x : d.x, d.y < 0 ? -d.y : d.y) + 1;
}
inline Coord cells_crossed(const Segment& s) {  // 線分が内部を通る 1×1 のマスの数（端点は格子点）
    const Point d = s.dir();
    const Coord dx = d.x < 0 ? -d.x : d.x, dy = d.y < 0 ? -d.y : d.y;
    return dx + dy - std::gcd(dx, dy);
}
#endif  // GG_LINE
#ifndef GG_CIRCLE
#define GG_CIRCLE
// ---- 円（幾何第二世代）  位置関係・交点の数は整数で厳密、交点・接点・中心は実数 ----
enum CircleRelation : int { CIRCLE_CONTAIN = 0, CIRCLE_INTERNAL = 1, CIRCLE_INTERSECT = 2, CIRCLE_EXTERNAL = 3, CIRCLE_SEPARATE = 4, CIRCLE_SAME = 5 };

template <class T>
Where contains(const CircleT<T>& c, const PointT<T>& p) {  // 円の内外（円周の上は ON_BORDER）
    if constexpr (std::is_integral_v<T>) {
        const Int128 d2 = norm2_128(p - c.c), r2 = (Int128)c.r * c.r;
        return d2 < r2 ? INSIDE : (d2 == r2 ? ON_BORDER : OUTSIDE);
    } else {
        const int s = real_cmp(dist(c.c, p), c.r);
        return s < 0 ? INSIDE : (s == 0 ? ON_BORDER : OUTSIDE);
    }
}
// 2 円の位置関係（値は共通接線の数。同じ円は CIRCLE_SAME = 5）
template <class T>
int circle_relation(const CircleT<T>& a, const CircleT<T>& b) {
    if constexpr (std::is_integral_v<T>) {
        const Int128 d2 = norm2_128(b.c - a.c), s = (Int128)a.r + b.r, df = a.r > b.r ? (Int128)a.r - b.r : (Int128)b.r - a.r;
        if (d2 == 0 && a.r == b.r) return CIRCLE_SAME;
        if (d2 > s * s) return CIRCLE_SEPARATE;
        if (d2 == s * s) return CIRCLE_EXTERNAL;
        if (d2 > df * df) return CIRCLE_INTERSECT;
        if (d2 == df * df) return CIRCLE_INTERNAL;
        return CIRCLE_CONTAIN;
    } else {
        const Real d = dist(a.c, b.c), s = a.r + b.r, df = std::fabs(a.r - b.r), sc = std::max(d, s);
        if (real_sign(d, sc) == 0 && real_sign(a.r - b.r, sc) == 0) return CIRCLE_SAME;
        const int c1 = real_sign(d - s, sc);
        if (c1 > 0) return CIRCLE_SEPARATE;
        if (c1 == 0) return CIRCLE_EXTERNAL;
        const int c2 = real_sign(d - df, sc);
        if (c2 > 0) return CIRCLE_INTERSECT;
        if (c2 == 0) return CIRCLE_INTERNAL;
        return CIRCLE_CONTAIN;
    }
}
template <class T>
bool has_common_point(const CircleT<T>& a, const CircleT<T>& b) {  // 円周どうしが共有点を持つか
    const int r = circle_relation(a, b);
    return r == CIRCLE_EXTERNAL || r == CIRCLE_INTERSECT || r == CIRCLE_INTERNAL || r == CIRCLE_SAME;
}
template <class T>
int count_cross(const CircleT<T>& c, const LineT<T>& l) {  // 円周と直線の交点の数（0 / 1（接する）/ 2）
    const PointT<T> d = l.dir();
    if constexpr (std::is_integral_v<T>) {
        const Int128 cr = cross128(d, c.c - l.a), lhs = cr * cr, rhs = (Int128)c.r * c.r * norm2_128(d);
        return lhs > rhs ? 0 : (lhs == rhs ? 1 : 2);
    } else {
        const int s = real_cmp(dist(l, c.c), c.r);
        return s > 0 ? 0 : (s == 0 ? 1 : 2);
    }
}
namespace geo_internal {
// 線分 s.a + t (s.b - s.a)（0 ≤ t ≤ 1）と円周の交点が、二次方程式 A t² + 2B t + C = 0 の解 t1 ≤ t2 のどれか
// （bit0: t1、bit1: t2。接する（重解）なら bit0 だけ）
template <class T>
int circle_segment_mask(const CircleT<T>& c, const SegmentT<T>& s) {
    const PointT<T> d = s.b - s.a, w = s.a - c.c;
    if constexpr (std::is_integral_v<T>) {
        const Int128 A = norm2_128(d), B = dot128(d, w), C = norm2_128(w) - (Int128)c.r * c.r;
        if (A == 0) return C == 0 ? 1 : 0;
        const Int128 D = B * B - A * C, f0 = C, f1 = A + 2 * B + C;
        if (D < 0 || (f0 < 0 && f1 < 0)) return 0;
        if (f0 < 0) return 2;
        if (f1 < 0) return 1;
        if (-B < 0) return f0 == 0 ? 2 : 0;
        if (-B > A) return f1 == 0 ? 1 : 0;
        return D == 0 ? 1 : 3;
    } else {
        const Real A = norm2(d), B = dot(d, w), nw = norm2(w), r2 = c.r * c.r, C = nw - r2;
        const Real sc = std::max({A, std::fabs(B), nw, r2});
        if (real_sign(A, sc) == 0) return real_sign(C, sc) == 0 ? 1 : 0;
        // C = |w|² - r² は桁落ちするので、誤差の目安は |w|² と r² の大きい方。判別式 B² - AC もそれに合わせる
        const Real cs = std::max(nw, r2);
        const int sD = real_sign(B * B - A * C, std::max(B * B, A * cs)), s0 = real_sign(C, cs), s1 = real_sign(A + 2 * B + C, sc);
        if (sD < 0 || (s0 < 0 && s1 < 0)) return 0;
        if (s0 < 0) return 2;
        if (s1 < 0) return 1;
        if (real_sign(-B, sc) < 0) return s0 == 0 ? 2 : 0;
        if (real_sign(-B - A, sc) > 0) return s1 == 0 ? 1 : 0;
        return sD == 0 ? 1 : 3;
    }
}
}  // namespace geo_internal
template <class T>
int count_cross(const CircleT<T>& c, const SegmentT<T>& s) {  // 円周と線分の交点の数（0〜2。端点も含む）
    const int m = geo_internal::circle_segment_mask(c, s);
    return (m & 1) + (m >> 1 & 1);
}

// ---- 交点・接点（実数） ----
template <class T>
std::vector<PointD> cross_points(const CircleT<T>& c, const LineT<T>& l) {  // 円周と直線の交点（直線の向き a → b の順）
    const int k = count_cross(c, l);
    if (k == 0) return {};
    const PointD h = projection(l, c.c);
    if (k == 1) return {h};
    const PointT<T> d = l.dir();
    Real half2;  // 半弦の長さの 2 乗 = r² - (中心と直線の距離)²
    if constexpr (std::is_integral_v<T>) {
        const Int128 cr = cross128(d, c.c - l.a), n2 = norm2_128(d);
        half2 = (Real)((Int128)c.r * c.r * n2 - cr * cr) / (Real)n2;
    } else {
        const Real e = dist(l, c.c);
        half2 = c.r * c.r - e * e;
    }
    const PointD u = unit(d) * std::sqrt(std::max<Real>(0, half2));
    return {h - u, h + u};
}
template <class T>
std::vector<PointD> cross_points(const CircleT<T>& c, const SegmentT<T>& s) {  // 円周と線分の交点（a → b の順）
    const int m = geo_internal::circle_segment_mask(c, s);
    if (m == 0) return {};
    if (s.a == s.b) return {to_d(s.a)};
    const std::vector<PointD> v = cross_points(c, LineT<T>(s.a, s.b));
    if (v.size() <= 1) return v;
    std::vector<PointD> r;
    if (m & 1) r.push_back(v[0]);
    if (m & 2) r.push_back(v[1]);
    return r;
}
// 2 円の交点（c1 → c2 の向きの左側の点が先。接するなら 1 個、同じ円・交わらないなら空）
template <class T>
std::vector<PointD> cross_points(const CircleT<T>& a, const CircleT<T>& b) {
    const int rel = circle_relation(a, b);
    if (rel == CIRCLE_SAME || rel == CIRCLE_SEPARATE || rel == CIRCLE_CONTAIN) return {};
    const PointT<T> v = b.c - a.c;
    PointD base, off;
    if constexpr (std::is_integral_v<T>) {
        const Int128 D2 = norm2_128(v), K = (Int128)a.r * a.r - (Int128)b.r * b.r + D2;  // K = 2 · d · (a.c から交点の弦までの距離)
        base = geo_internal::affine(a.c, v, K, 2 * D2);
        if (rel == CIRCLE_EXTERNAL || rel == CIRCLE_INTERNAL) return {base};
        const Int128 H = 4 * D2 * a.r * a.r - K * K;  // (2 · d · 半弦)²
        off = to_d(rot90(v)) * (std::sqrt((Real)(H < 0 ? 0 : H)) / (2 * (Real)D2));
    } else {
        const Real D2 = norm2(v), K = a.r * a.r - b.r * b.r + D2;
        base = a.c + v * (K / (2 * D2));
        if (rel == CIRCLE_EXTERNAL || rel == CIRCLE_INTERNAL) return {base};
        const Real H = 4 * D2 * a.r * a.r - K * K;
        off = rot90(v) * (std::sqrt(std::max<Real>(0, H)) / (2 * D2));
    }
    return {base + off, base - off};
}
// p から円に引いた接線の接点（p が円周の上なら {p}、内側なら空）
template <class T>
std::vector<PointD> tangent_points(const CircleT<T>& c, const PointT<T>& p) {
    const Where w = contains(c, p);
    if (w == INSIDE) return {};
    if (w == ON_BORDER) return {to_d(p)};
    const PointT<T> v = p - c.c;
    PointD base, off;
    if constexpr (std::is_integral_v<T>) {
        const Int128 D2 = norm2_128(v), r2 = (Int128)c.r * c.r;
        base = geo_internal::affine(c.c, v, r2, D2);
        off = to_d(rot90(v)) * ((Real)c.r * std::sqrt((Real)(D2 - r2)) / (Real)D2);
    } else {
        const Real D2 = norm2(v), r2 = c.r * c.r;
        base = c.c + v * (r2 / D2);
        off = rot90(v) * (c.r * std::sqrt(std::max<Real>(0, D2 - r2)) / D2);
    }
    return {base + off, base - off};
}
// 共通接線。a は c1 との接点、b は c2 との接点（接点が同じなら a + 接線の向き）。外側の接線 → 内側の接線の順
// 本数は circle_relation の値（4 / 3 / 2 / 1 / 0）。同じ円・中心が同じなら空。半径は正
template <class T>
std::vector<LineD> common_tangents(const CircleT<T>& a, const CircleT<T>& b) {
    const int rel = circle_relation(a, b);
    std::vector<LineD> res;
    if (rel == CIRCLE_SAME || rel == CIRCLE_CONTAIN) return res;
    const PointD c1 = to_d(a.c), c2 = to_d(b.c), v = to_d(b.c - a.c), w = rot90(v);
    const Real r1 = (Real)a.r, r2 = (Real)b.r;
    Real D2, e2, i2;  // |v|², |v|² - (r1 - r2)², |v|² - (r1 + r2)²
    if constexpr (std::is_integral_v<T>) {
        const Int128 d2 = norm2_128(b.c - a.c), dm = (Int128)a.r - b.r, dp = (Int128)a.r + b.r;
        D2 = (Real)d2, e2 = (Real)(d2 - dm * dm), i2 = (Real)(d2 - dp * dp);
    } else {
        D2 = norm2(v), e2 = D2 - (r1 - r2) * (r1 - r2), i2 = D2 - (r1 + r2) * (r1 + r2);
    }
    auto add = [&](Real cosn, Real sin2, int sgn, bool inner) {  // 法線 n = (v · cosn ± w · √sin2) / |v|²
        const PointD n = (v * cosn + w * (sgn * std::sqrt(std::max<Real>(0, sin2)))) / D2;
        const PointD q1 = c1 + n * r1, q2 = inner ? c2 - n * r2 : c2 + n * r2;
        res.push_back(real_eq(q1, q2) ? LineD(q1, q1 + rot90(n)) : LineD(q1, q2));
    };
    add(r1 - r2, e2, 1, false);  // 外側の接線
    if (rel != CIRCLE_INTERNAL) add(r1 - r2, e2, -1, false);
    if (rel == CIRCLE_EXTERNAL || rel == CIRCLE_SEPARATE) {  // 内側の接線
        add(r1 + r2, i2, 1, true);
        if (rel == CIRCLE_SEPARATE) add(r1 + r2, i2, -1, true);
    }
    return res;
}

// ---- 三角形の円・2 点を通る円 ----
template <class T>
CircleD circumcircle(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c) {  // 外接円（一直線なら r = -1）
    if (ccw(a, b, c) == 0) return CircleD(to_d(a), -1);
    const PointT<T> u = b - a, v = c - a;
    PointD o;
    if constexpr (std::is_integral_v<T>) {
        const Int128 D = 2 * cross128(u, v), nu = norm2_128(u), nv = norm2_128(v);
        const Int128 X = (Int128)v.y * nu - (Int128)u.y * nv, Y = (Int128)u.x * nv - (Int128)v.x * nu;
        o = PointD((Real)((Int128)a.x * D + X) / (Real)D, (Real)((Int128)a.y * D + Y) / (Real)D);
    } else {
        const Real D = 2 * cross(u, v), nu = norm2(u), nv = norm2(v);
        o = a + PointD(v.y * nu - u.y * nv, u.x * nv - v.x * nu) / D;
    }
    return CircleD(o, dist(o, to_d(a)));
}
namespace geo_internal {
inline int sign_sum_mul3(const Int128 (&a)[3], const long long (&b)[3]) {  // Σ a_i b_i の符号（|a_i| < 2^120。積は 192 ビットまで誤差なし）
    const UInt128 mask = ((UInt128)1 << 64) - 1;
    Int128 H = 0, L = 0;  // 和 = H · 2^64 + L
    for (int i = 0; i < 3; ++i) {
        const Int128 P = (Int128)((UInt128)a[i] & mask) * b[i];
        H += (a[i] >> 64) * b[i] + (P >> 64), L += (Int128)((UInt128)P & mask);
    }
    H += L >> 64, L &= (Int128)mask;
    return H > 0 ? 1 : (H < 0 ? -1 : (L > 0 ? 1 : 0));
}
}  // namespace geo_internal
// 点 p が 3 点 a, b, c を通る円（外接円）の内側か（OUTSIDE / ON_BORDER / INSIDE）。整数なら厳密。a, b, c が一直線なら OUTSIDE
// 4 点が同じ円の上にあるか = in_circumcircle(a, b, c, d) == ON_BORDER
template <class T>
Where in_circumcircle(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c, const PointT<T>& p) {
    const int o = ccw(a, b, c);
    if (o == 0) return OUTSIDE;
    const PointT<T> u = a - p, v = b - p, w = c - p;
    int s;
    if constexpr (std::is_integral_v<T>) {
        const Int128 nu = norm2_128(u), nv = norm2_128(v), nw = norm2_128(w);
        // 行列式 | u.x u.y |u|² ; v.x v.y |v|² ; w.x w.y |w|² | を 1 行目で展開（|u|² ≤ 8×10^18 は long long に収まる）
        const Int128 A[3] = {(Int128)v.y * nw - (Int128)w.y * nv, (Int128)w.x * nv - (Int128)v.x * nw, cross128(v, w)};
        const long long B[3] = {u.x, u.y, (long long)nu};
        s = geo_internal::sign_sum_mul3(A, B);
    } else {
        const Real nu = norm2(u), nv = norm2(v), nw = norm2(w), m = std::max({nu, nv, nw});
        s = real_sign(u.x * (v.y * nw - w.y * nv) - u.y * (v.x * nw - w.x * nv) + nu * cross(v, w), m * m);
    }
    s *= o;
    return s > 0 ? INSIDE : (s == 0 ? ON_BORDER : OUTSIDE);
}
template <class T>
CircleD incircle(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c) {  // 内接円（一直線なら r = -1）
    if (ccw(a, b, c) == 0) return CircleD(to_d(a), -1);
    const Real la = dist(b, c), lb = dist(c, a), lc = dist(a, b), s = la + lb + lc;
    const PointD o = (to_d(a) * la + to_d(b) * lb + to_d(c) * lc) / s;
    const T cr = cross(b - a, c - a);
    return CircleD(o, std::fabs((Real)cr) / s);
}
template <class T>
CircleD circle_with_diameter(const PointT<T>& a, const PointT<T>& b) { return CircleD(midpoint(a, b), dist(a, b) / 2); }
// a と b（a ≠ b）を通る半径 r の円の中心（0〜2 個）。整数の点と整数の r なら個数は厳密
template <class T, class R>
std::vector<PointD> circles_through(const PointT<T>& a, const PointT<T>& b, R r) {
    const PointT<T> v = b - a;
    if (a == b || r < 0) return {};
    int k;
    Real h2;  // (2h)² = 4r² - |v|²（h は中点から中心までの距離）
    if constexpr (std::is_integral_v<T> && std::is_integral_v<R>) {
        const Int128 d2 = norm2_128(v), f = 4 * (Int128)r * r;
        k = d2 > f ? 0 : (d2 == f ? 1 : 2);
        h2 = (Real)(f - d2);
    } else {
        const Real d = len(v);
        const int s = real_cmp(d, 2 * (Real)r);
        k = s > 0 ? 0 : (s == 0 ? 1 : 2);
        h2 = 4 * (Real)r * (Real)r - d * d;
    }
    if (k == 0) return {};
    const PointD m = midpoint(a, b);
    if (k == 1) return {m};
    const PointD off = unit(rot90(v)) * (std::sqrt(std::max<Real>(0, h2)) / 2);
    return {m + off, m - off};
}
// ---- 円と半直線（円周との交点） ----
template <class T>
int count_cross(const CircleT<T>& c, const RayT<T>& r) {  // 円周と半直線の交点の数（0〜2。始点も含む。整数なら厳密）
    // r.a + t d（t ≥ 0）を代入した A t² + 2B t + C = 0 の解のうち t ≥ 0 のものの数
    const PointT<T> d = r.dir(), w = r.a - c.c;
    int sD, sB, sC;
    if constexpr (std::is_integral_v<T>) {
        const Int128 A = norm2_128(d), B = dot128(d, w), C = norm2_128(w) - (Int128)c.r * c.r, D = B * B - A * C;
        sD = (D > 0) - (D < 0), sB = (B > 0) - (B < 0), sC = (C > 0) - (C < 0);
    } else {
        const Real A = norm2(d), B = dot(d, w), nw = norm2(w), r2 = c.r * c.r, C = nw - r2, cs = std::max(nw, r2);
        sD = real_sign(B * B - A * C, std::max(B * B, A * cs)), sB = real_sign(B, std::sqrt(A * nw)), sC = real_sign(C, cs);
    }
    if (sD < 0) return 0;
    if (sD == 0) return sB <= 0 ? 1 : 0;
    return (sB <= 0 || sC <= 0) + (sB <= 0 && sC >= 0);
}
template <class T>
std::vector<PointD> cross_points(const CircleT<T>& c, const RayT<T>& r) {  // 円周と半直線の交点（始点に近い順）
    const int k = count_cross(c, r);
    if (k == 0) return {};
    std::vector<PointD> v = cross_points(c, LineT<T>(r.a, r.b));
    if (v.size() <= 1 || k == 2) return v;
    return {v[1]};
}

// ---- 円板（円の内側。円周を含む）としての位置関係・距離。整数なら厳密 ----
namespace geo_internal {
template <class T, class X>
int disk_cmp(const CircleT<T>& c, const X& x) {  // (中心と x の距離) と r の比較（-1 / 0 / 1）
    if constexpr (std::is_integral_v<T>) return dist2_cmp(x, c.c, c.r * c.r);
    else return real_cmp(dist(x, c.c), c.r);
}
}  // namespace geo_internal
template <class T>
bool disk_intersects(const CircleT<T>& c, const SegmentT<T>& s) { return geo_internal::disk_cmp(c, s) <= 0; }  // 円板と線分が共有点を持つか
template <class T>
bool disk_intersects(const CircleT<T>& c, const LineT<T>& l) { return geo_internal::disk_cmp(c, l) <= 0; }
template <class T>
bool disk_intersects(const CircleT<T>& c, const RayT<T>& r) { return geo_internal::disk_cmp(c, r) <= 0; }
template <class T>
bool disk_intersects(const CircleT<T>& a, const CircleT<T>& b) {  // 2 つの円板が共有点を持つか（中心の距離 ≤ r1 + r2）
    if constexpr (std::is_integral_v<T>) return norm2_128(b.c - a.c) <= ((Int128)a.r + b.r) * ((Int128)a.r + b.r);
    else return real_cmp(dist(a.c, b.c), a.r + b.r) <= 0;
}
template <class T>
bool disk_contains(const CircleT<T>& c, const SegmentT<T>& s) {  // 線分が円板に入っているか（両端が入っていればよい）
    return contains(c, s.a) != OUTSIDE && contains(c, s.b) != OUTSIDE;
}
template <class T>
bool disk_contains(const CircleT<T>& a, const CircleT<T>& b) {  // 円板 b が円板 a に入っているか（内接・同じ円も含む）
    if constexpr (std::is_integral_v<T>) {
        if (a.r < b.r) return false;
        const Int128 df = (Int128)a.r - b.r;
        return norm2_128(b.c - a.c) <= df * df;
    } else {
        return real_cmp(dist(a.c, b.c) + b.r, a.r) <= 0;
    }
}
// 円板までの距離（円板の中・共有点があれば 0）
template <class T>
Real disk_dist(const CircleT<T>& c, const PointT<T>& p) { return contains(c, p) != OUTSIDE ? 0 : dist(c.c, p) - (Real)c.r; }
template <class T>
Real disk_dist(const CircleT<T>& c, const SegmentT<T>& s) { return disk_intersects(c, s) ? 0 : dist(s, c.c) - (Real)c.r; }
template <class T>
Real disk_dist(const CircleT<T>& c, const LineT<T>& l) { return disk_intersects(c, l) ? 0 : dist(l, c.c) - (Real)c.r; }
template <class T>
Real disk_dist(const CircleT<T>& c, const RayT<T>& r) { return disk_intersects(c, r) ? 0 : dist(r, c.c) - (Real)c.r; }
template <class T>
Real disk_dist(const CircleT<T>& a, const CircleT<T>& b) { return disk_intersects(a, b) ? 0 : dist(a.c, b.c) - (Real)a.r - (Real)b.r; }
// 円周までの距離（中の点でも円周まで）。2 円は円周どうしの距離（共有点があれば 0。片方が中にあれば内側のすき間）
template <class T>
Real circumference_dist(const CircleT<T>& c, const PointT<T>& p) { return std::fabs(dist(c.c, p) - (Real)c.r); }
template <class T>
Real circumference_dist(const CircleT<T>& a, const CircleT<T>& b) {
    const int rel = circle_relation(a, b);
    if (rel == CIRCLE_SEPARATE) return dist(a.c, b.c) - (Real)a.r - (Real)b.r;
    if (rel == CIRCLE_CONTAIN) return std::fabs((Real)a.r - (Real)b.r) - dist(a.c, b.c);
    return 0;
}

// ---- 円周の上の点を結ぶ弦（点は円周に沿った位置で表す：番号・周に沿った座標など、一周の順に並んだ整数） ----
// 弦 ab と cd が円の内部の 1 点で交わるか（端点を共有する・長さ 0 なら false）
inline bool chords_cross(long long a, long long b, long long c, long long d) {
    if (a > b) std::swap(a, b);
    if (c > d) std::swap(c, d);
    return (a < c && c < b && b < d) || (c < a && a < d && d < b);
}
// 円周を n 等分した点 a, b を通る直線の傾きのキー（2 本が平行 ⇔ キーが同じ）。点の番号は 0 〜 n-1 でも 1 〜 n でもよい
inline long long chord_direction_key(long long a, long long b, long long n) { return ((a + b) % n + n) % n; }
// 各 query の弦が、chords のうち何本と交わるか（chords_cross と同じ意味）。オフライン O((N + Q) log N)
inline std::vector<long long> chord_cross_counts(const std::vector<std::pair<long long, long long>>& chords,
                                                 const std::vector<std::pair<long long, long long>>& queries) {
    std::vector<std::pair<long long, long long>> ch;
    for (auto [a, b] : chords) {
        if (a > b) std::swap(a, b);
        if (a < b) ch.emplace_back(a, b);
    }
    std::sort(ch.begin(), ch.end());
    std::vector<long long> bs;
    for (const auto& e : ch) bs.push_back(e.second);
    std::sort(bs.begin(), bs.end());
    const int m = (int)bs.size();
    std::vector<int> bit(m + 1, 0);
    auto add = [&](long long b) {
        for (int i = (int)(std::lower_bound(bs.begin(), bs.end(), b) - bs.begin()) + 1; i <= m; i += i & -i) ++bit[i];
    };
    auto count_less = [&](long long b) {  // 入れた弦のうち他端 < b の数
        long long r = 0;
        for (int i = (int)(std::lower_bound(bs.begin(), bs.end(), b) - bs.begin()); i > 0; i -= i & -i) r += bit[i];
        return r;
    };
    // 交わる ⇔ (c < a < d かつ b > d) または (a < c かつ c < b < d)。「a < X の弦のうち…」の差で数える
    struct Ev {
        long long X, lo, hi;  // a < X の弦のうち lo ≤ b < hi のもの
        int q, sgn;
    };
    std::vector<Ev> ev;
    const long long no_upper = std::numeric_limits<long long>::max();  // 上限なし
    for (int q = 0; q < (int)queries.size(); ++q) {
        auto [c, d] = queries[q];
        if (c > d) std::swap(c, d);
        if (c == d) continue;
        ev.push_back({d, d + 1, no_upper, q, +1}), ev.push_back({c + 1, d + 1, no_upper, q, -1});
        ev.push_back({c, c + 1, d, q, +1});
    }
    std::sort(ev.begin(), ev.end(), [](const Ev& x, const Ev& y) { return x.X < y.X; });
    std::vector<long long> res(queries.size(), 0);
    size_t k = 0;
    for (const Ev& e : ev) {
        while (k < ch.size() && ch[k].first < e.X) add(ch[k++].second);
        const long long hi = e.hi == no_upper ? (long long)k : count_less(e.hi);
        res[e.q] += e.sgn * (hi - count_less(e.lo));
    }
    return res;
}
// 交わる弦の組の数（端点を共有する組は数えない）。O(N log N)
inline long long count_crossing_chords(const std::vector<std::pair<long long, long long>>& chords) {
    long long s = 0;
    for (long long v : chord_cross_counts(chords, chords)) s += v;
    return s / 2;
}
// 弦を足しながら「弦 ab が今ある弦のどれかと交わるか」を調べる（位置は 0 〜 n。1 〜 n の番号もそのまま使える）。O(log n)
struct ChordSet {
    int sz = 1;
    std::vector<long long> mx, mn;  // 位置 a に端がある弦の、もう一方の端の最大（a が小さい側）/ 最小（a が大きい側）
    explicit ChordSet(int n) {
        while (sz < n + 1) sz <<= 1;
        mx.assign(2 * sz, std::numeric_limits<long long>::min()), mn.assign(2 * sz, std::numeric_limits<long long>::max());
    }
    void add(long long a, long long b) {  // 弦 ab を足す
        if (a > b) std::swap(a, b);
        if (a == b) return;
        for (long long i = a + sz; i > 0; i >>= 1) mx[i] = std::max(mx[i], b);
        for (long long i = b + sz; i > 0; i >>= 1) mn[i] = std::min(mn[i], a);
    }
    bool crosses(long long c, long long d) const {  // 弦 cd が今ある弦のどれかと円の内部で交わるか
        if (c > d) std::swap(c, d);
        long long hi = std::numeric_limits<long long>::min(), lo = std::numeric_limits<long long>::max();
        for (long long l = c + 1 + sz, r = d + sz; l < r; l >>= 1, r >>= 1) {  // 位置 (c, d) の範囲
            if (l & 1) hi = std::max(hi, mx[l]), lo = std::min(lo, mn[l]), ++l;
            if (r & 1) --r, hi = std::max(hi, mx[r]), lo = std::min(lo, mn[r]);
        }
        return hi > d || lo < c;
    }
};

inline Real arc_length(Real r, Real theta) { return r * theta; }
inline Real sector_area(Real r, Real theta) { return r * r * theta / 2; }

// 円の中の格子点 (X, Y) の数。中心・半径を scale 倍した整数 cx, cy, r で渡す（(X·scale - cx)² + (Y·scale - cy)² ≤ r²）
// 小数点以下 k 桁の入力なら parse_fixed(s, k) / read_fixed(k) と scale = 10^k。include_border = false なら円周の上を除く。O(r / scale)
inline long long lattice_points_in_circle(Coord cx, Coord cy, Coord r, Coord scale = 1, bool include_border = true) {
    if (r < 0) return 0;
    auto fdiv = [](Coord a, Coord b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); };
    auto cdiv = [&](Coord a, Coord b) { return -fdiv(-a, b); };
    const Int128 r2 = (Int128)r * r;
    long long cnt = 0;
    for (Coord X = cdiv(cx - r, scale), xe = fdiv(cx + r, scale); X <= xe; ++X) {
        const Coord dx = X * scale - cx;
        const Int128 rem = r2 - (Int128)dx * dx;
        if (rem < 0 || (!include_border && rem == 0)) continue;
        const Coord s = isqrt((long long)(include_border ? rem : rem - 1));  // |Y·scale - cy| ≤ s
        const Coord lo = cdiv(cy - s, scale), hi = fdiv(cy + s, scale);
        if (hi >= lo) cnt += hi - lo + 1;
    }
    return cnt;
}
#endif  // GG_CIRCLE
#ifndef GG_CIRCLE_COVER
#define GG_CIRCLE_COVER
// ---- 点を円で囲む（幾何第二世代）  最小包含円・半径 r の円で囲める最大の点の数（答えは実数） ----
struct CircleCover {  // 囲める点の数と、そのときの円の中心
    int count = 0;
    PointD center;
};

namespace geo_internal {
// 円の内側か円周の上。誤差は「半径 × GEO_EPS」と「中心の座標を Real で表すときの丸め（座標の大きさ × 64 ε）」の大きい方
//（座標の大きさ × GEO_EPS にすると、座標 10^9 で 1 だけ離れた点を同じ点とみなしてしまう）
inline bool cover_in(const CircleD& c, const PointD& p) {
    const Real tol = std::max(GEO_EPS * std::max<Real>(1, c.r), 64 * std::numeric_limits<Real>::epsilon() * std::max(std::fabs(c.c.x), std::fabs(c.c.y)));
    return dist(c.c, p) <= c.r + tol;
}
template <class T>
CircleD cover_circle3(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c) {  // 3 点を通る円（一直線なら最も遠い 2 点を直径に）
    const CircleD cc = circumcircle(a, b, c);
    if (cc.r >= 0) return cc;
    const Real ab = dist2(to_d(a), to_d(b)), bc = dist2(to_d(b), to_d(c)), ca = dist2(to_d(c), to_d(a));
    if (ab >= bc && ab >= ca) return circle_with_diameter(a, b);
    return bc >= ca ? circle_with_diameter(b, c) : circle_with_diameter(c, a);
}
}  // namespace geo_internal

// 全部の点を囲む最小の円（Welzl の方法。点の順番を乱数で混ぜて期待 O(N)）。点がなければ r = -1
template <class T>
CircleD min_enclosing_circle(const std::vector<PointT<T>>& ps) {
    const int n = (int)ps.size();
    if (n == 0) return CircleD(PointD(), -1);
    std::vector<int> id(n);
    for (int i = 0; i < n; ++i) id[i] = i;
    std::mt19937 rng(20261008u);
    std::shuffle(id.begin(), id.end(), rng);
    auto P = [&](int k) { return ps[id[k]]; };
    CircleD c(to_d(P(0)), 0);
    for (int i = 1; i < n; ++i) {
        if (geo_internal::cover_in(c, to_d(P(i)))) continue;
        c = CircleD(to_d(P(i)), 0);
        for (int j = 0; j < i; ++j) {
            if (geo_internal::cover_in(c, to_d(P(j)))) continue;
            c = circle_with_diameter(P(i), P(j));
            for (int k = 0; k < j; ++k)
                if (!geo_internal::cover_in(c, to_d(P(k)))) c = geo_internal::cover_circle3(P(i), P(j), P(k));
        }
    }
    return c;
}

// 半径 r の円で囲める点の最大の数（円周の上も含む）と、その円の中心。どれか 1 点を円周の上に置いて偏角で走査する
template <class T, class R>
CircleCover max_points_in_circle(const std::vector<PointT<T>>& ps, R r_in) {
    const int n = (int)ps.size();
    const Real r = (Real)r_in;
    if (n == 0 || r < 0) return CircleCover();
    CircleCover best{1, to_d(ps[0])};
    std::vector<std::pair<Real, int>> ev;
    for (int i = 0; i < n; ++i) {
        ev.clear();
        int cur = 1;
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            if (ps[j] == ps[i]) { ++cur; continue; }
            const PointD v = to_d(ps[j]) - to_d(ps[i]);
            const Real d = len(v);
            if (d > 2 * r && real_sign(d - 2 * r, d) > 0) continue;
            // 中心 c = ps[i] + r (cos t, sin t) のうち ps[j] を囲むのは t ∈ [θ - φ, θ + φ]
            const Real th = angle_of(v), ph = std::acos(std::min<Real>(1, d / (2 * r))) + GEO_EPS;
            Real s = th - ph, e = th + ph;
            if (s < 0) s += 2 * GEO_PI, e += 2 * GEO_PI;
            ev.push_back({s, -1});  // 同じ角度なら入る方を先に（-1 < +1 の順に並べる）
            if (e >= 2 * GEO_PI) ++cur, ev.push_back({e - 2 * GEO_PI, +1});
            else ev.push_back({e, +1});
        }
        std::sort(ev.begin(), ev.end());
        auto center = [&](Real t) { return to_d(ps[i]) + PointD(r * std::cos(t), r * std::sin(t)); };
        if (cur > best.count) best = CircleCover{cur, center(0)};
        for (const auto& [t, k] : ev) {
            cur -= k;  // 入る（-1）なら +1、出る（+1）なら -1
            if (cur > best.count) best = CircleCover{cur, center(t)};
        }
    }
    return best;
}
#endif  // GG_CIRCLE_COVER
// ======================= ここまでテンプレート本体 =======================


// cin を文字列から読むようにする（テスト・デモ用）
struct CinFromString {
    std::istringstream ss;
    std::streambuf* old;
    explicit CinFromString(const std::string& s) : ss(s), old(std::cin.rdbuf(ss.rdbuf())) {}
    ~CinFromString() { std::cin.rdbuf(old); }
};

// =============================================================================
//  ランダムテスト
// =============================================================================
namespace cv_test {

#define CV_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[circle_cover] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK CV_CHECK

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
static const int TEST_SCALE = 1;
#else
static const int TEST_SCALE = 8;
#endif
// 実数の許容誤差（Real の精度に合わせる。long double なら 1e-15、double なら 1e-12 程度）
[[maybe_unused]] static const Real EPSR = std::numeric_limits<Real>::epsilon();
[[maybe_unused]] static const Real TOL = std::max<Real>(EPSR * 4096, (Real)1e-15L);

[[maybe_unused]] static long long rnd(std::mt19937_64& rng, long long lo, long long hi) { return lo + (long long)(rng() % (unsigned long long)(hi - lo + 1)); }
[[maybe_unused]] static Real rndr(std::mt19937_64& rng, Real lo, Real hi) { return lo + (hi - lo) * (Real)(rng() >> 11) / (Real)(1ULL << 53); }
[[maybe_unused]] static Point rndp(std::mt19937_64& rng, long long c) { return Point(rnd(rng, -c, c), rnd(rng, -c, c)); }
[[maybe_unused]] static std::vector<Point> rndps(std::mt19937_64& rng, int n, long long c) {
    std::vector<Point> ps(n);
    for (auto& p : ps) p = rndp(rng, c);
    return ps;
}
[[maybe_unused]] static int sgn128(Int128 v) { return (v > 0) - (v < 0); }
[[maybe_unused]] static bool near(Real a, Real b, Real tol) { return std::fabs(a - b) <= tol * std::max<Real>(1, std::max(std::fabs(a), std::fabs(b))); }
[[maybe_unused]] static bool nearp(const PointD& a, const PointD& b, Real tol) {
    const Real s = std::max<Real>(1, std::max(std::max(std::fabs(a.x), std::fabs(a.y)), std::max(std::fabs(b.x), std::fabs(b.y))));
    return std::fabs(a.x - b.x) <= tol * s && std::fabs(a.y - b.y) <= tol * s;
}
[[maybe_unused]] static long long ms_since(std::chrono::steady_clock::time_point a) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
}

static bool b_in(const CircleD& c, const PointD& p, Real tol) { return dist(c.c, p) <= c.r + tol * std::max<Real>(1, std::max({c.r, std::fabs(c.c.x), std::fabs(c.c.y)})); }

static void test_mec(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 6000 * TEST_SCALE; ++it, ++cases) {
        const bool big = it % 3 == 0;
        const long long C = big ? 1000000000LL : (it % 3 == 1 ? 2 : 20);
        const int n = (int)rnd(rng, 1, 9);
        std::vector<Point> ps = rndps(rng, n, C);
        if (it % 7 == 0)  // 一直線・同じ点
            for (auto& p : ps) p = Point(1, 2) + Point(3, -1) * rnd(rng, -2, 2);
        const CircleD c = min_enclosing_circle(ps);
        for (const Point& p : ps) CHECK(b_in(c, to_d(p), 1e-12L));
        // 愚直：1 点・2 点を直径・3 点の外接円の中で、全部を囲む最小の円
        Real best = std::numeric_limits<Real>::infinity();
        auto consider = [&](const CircleD& cc) {
            if (cc.r < 0) return;
            for (const Point& p : ps)
                if (!b_in(cc, to_d(p), 1e-12L)) return;
            best = std::min(best, cc.r);
        };
        for (int i = 0; i < n; ++i) {
            consider(CircleD(to_d(ps[i]), 0));
            for (int j = i + 1; j < n; ++j) {
                consider(circle_with_diameter(ps[i], ps[j]));
                for (int k = j + 1; k < n; ++k) consider(circumcircle(ps[i], ps[j], ps[k]));
            }
        }
        CHECK(std::fabs(c.r - best) <= 1e-12L * std::max<Real>(1, (Real)C));
        // 実数の点でも同じ
        std::vector<PointD> pd;
        for (const Point& p : ps) pd.push_back(to_d(p));
        CHECK(std::fabs(min_enclosing_circle(pd).r - best) <= 1e-9L * std::max<Real>(1, (Real)C));
    }
    CHECK(min_enclosing_circle(std::vector<Point>{}).r == -1);
    // 原点から遠い（座標 10^9 近く）小さい点の集まり：半径は小さいまま正しく
    CHECK(near(min_enclosing_circle(std::vector<Point>{{1000000000, 1000000000}, {1000000000, 999999999}}).r, 0.5L, TOL * 1e9L));
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const Point o(rnd(rng, 999999000, 1000000000) * (rng() % 2 ? 1 : -1), rnd(rng, 999999000, 1000000000) * (rng() % 2 ? 1 : -1));
        const int n = (int)rnd(rng, 1, 8);
        std::vector<Point> ps;
        for (int k = 0; k < n; ++k) {
            const Point q = o + rndp(rng, 6);
            if (std::llabs(q.x) <= 1000000000 && std::llabs(q.y) <= 1000000000) ps.push_back(q);
        }
        if (ps.empty()) continue;
        std::vector<Point> qs = ps;  // 原点の近くへ平行移動したもの（こちらは誤差が小さい）
        for (Point& q : qs) q -= o;
        const CircleD c = min_enclosing_circle(ps), cq = min_enclosing_circle(qs);
        const Real tol = std::max<Real>(1e-9L, EPSR * 1e9L * 256);
        CHECK(std::fabs(c.r - cq.r) <= tol && nearp(c.c - to_d(o), cq.c, tol));
        for (const Point& p : ps) CHECK(dist(c.c, to_d(p)) <= c.r + tol);
    }
}

static void test_max_points(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 3 : 8;
        const int n = (int)rnd(rng, 0, 12);
        std::vector<Point> ps = rndps(rng, n, C);
        // 半径：整数・半整数・√整数（ちょうど円周の上に乗る点が出るように）
        Real r;
        const int t = it % 3;
        if (t == 0) r = (Real)rnd(rng, 0, 4);
        else if (t == 1) r = (Real)rnd(rng, 1, 9) / 2;
        else r = std::sqrt((Real)rnd(rng, 1, 30));
        const CircleCover got = max_points_in_circle(ps, r);
        // 愚直：中心の候補は各点と、2 点を通る半径 r の円の中心
        int best = 0;
        auto count = [&](const PointD& c, Real tol) {
            int k = 0;
            for (const Point& p : ps) k += dist(c, to_d(p)) <= r + tol;
            return k;
        };
        for (int i = 0; i < n; ++i) {
            best = std::max(best, count(to_d(ps[i]), 1e-9L));
            for (int j = i + 1; j < n; ++j)
                for (const PointD& c : circles_through(to_d(ps[i]), to_d(ps[j]), r)) best = std::max(best, count(c, 1e-9L));
        }
        CHECK(got.count == best);
        if (n > 0) CHECK(count(got.center, 1e-7L) >= got.count);
    }
    // 大きい座標・整数の半径
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 30);
        const std::vector<Point> ps = rndps(rng, n, 1000000000LL);
        const long long r = rnd(rng, 0, 1500000000LL);
        const CircleCover got = max_points_in_circle(ps, r);
        int best = 1;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                for (const PointD& c : circles_through(ps[i], ps[j], r)) {
                    int k = 0;
                    for (const Point& p : ps) k += dist(c, to_d(p)) <= (Real)r * (1 + 1e-12L);
                    best = std::max(best, k);
                }
        CHECK(got.count == best);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_mec(rng, cases);
    test_max_points(rng, cases);
    std::cout << "[circle_cover] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[circle_cover] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 200000;
    const std::vector<Point> ps = rndps(rng, n, 1000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    const CircleD c = min_enclosing_circle(ps);
    const long long t1 = ms_since(t0);
    const int m = 2000;
    const std::vector<Point> qs = rndps(rng, m, 1000);
    t0 = std::chrono::steady_clock::now();
    const CircleCover cv = max_points_in_circle(qs, 100);
    const long long t2 = ms_since(t0);
    std::cout << "[circle_cover] benchmark: min_enclosing_circle N = " << n << " " << t1 << " ms, max_points_in_circle N = " << m << " " << t2
              << " ms  (check " << ((long long)c.r + cv.count) % 7 << ")\n";
}

}  // namespace cv_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace cv_demo {

[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
template <class T, std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool>, int> = 0>
[[maybe_unused]] static std::string to_s(T v) { return std::to_string(v); }
[[maybe_unused]] static std::string to_s(Int128 v) { return int128_to_string(v); }
[[maybe_unused]] static std::string to_s(Real v) {
    if (std::fabs(v) < 5e-7L) v = 0;
    std::ostringstream os;
    os << std::fixed << std::setprecision(6) << v;
    return os.str();
}
template <class T, std::enable_if_t<std::is_floating_point_v<T> && !std::is_same_v<T, Real>, int> = 0>
[[maybe_unused]] static std::string to_s(T v) { return to_s((Real)v); }
[[maybe_unused]] static std::string to_s(const Point& p) { return "(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ")"; }
[[maybe_unused]] static std::string to_s(const PointD& p) { return "(" + to_s(p.x) + ", " + to_s(p.y) + ")"; }
[[maybe_unused]] static std::string to_s(const std::string& s) { return s; }
[[maybe_unused]] static std::string to_s(const char* s) { return s; }
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
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
[[maybe_unused]] static void show(const char* expr, const std::string& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 34) + " -> " + v;
    if (!comment.empty()) line = pad(line, 58) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) cv_demo::show(#expr, cv_demo::to_s(expr), comment)
#define RUN(...)                   \
    cv_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("点を円で囲む 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  min_enclosing_circle(ps)          全部の点を囲む最小の円（CircleD）\n"
              << "  max_points_in_circle(ps, r)       半径 r の円で囲める最大の点の数と中心（CircleCover{count, center}）\n";

    api("【1】 min_enclosing_circle(ps)", "CircleD", "全部の点を囲む最小の円（円周の上も可）。点がなければ r = -1", "期待 O(N)");
    RUN(std::vector<Point> ps = {{0, 0}, {4, 0}, {2, 1}, {1, 3}});
    RUN(CircleD c = min_enclosing_circle(ps));
    SHOW(c.c, "中心");
    SHOW(c.r, "半径");
    SHOW(min_enclosing_circle(std::vector<Point>{{0, 0}, {6, 0}, {3, 1}}).r, "鈍角三角形 → 最も遠い 2 点を直径に（半径 3）");
    note("ABC151 F（Enclose All）: 半径をそのまま出力する（cout << fixed << setprecision(15) << c.r）");

    api("【2】 max_points_in_circle(ps, r)", "CircleCover{count, center}", "半径 r の円で囲める点の最大の数（円周の上も数える）と、その中心", "O(N^2 log N)");
    RUN(std::vector<Point> qs = {{0, 0}, {2, 0}, {1, 1}, {5, 5}, {6, 5}, {9, 9}});
    RUN(CircleCover cv = max_points_in_circle(qs, 1));
    SHOW(cv.count, "(0,0), (2,0), (1,1) は中心 (1, 0) の半径 1 の円で囲める");
    SHOW(cv.center, "中心の一つ（円周の上にどれかの点がある）");
    SHOW(max_points_in_circle(qs, 0.5).count, "半径 0.5 なら (5,5), (6,5) の 2 点");
    note("考え方：最適な円は動かしていくと、どれかの点が円周の上に来る → 各点を円周の上に置いて、中心の偏角の区間を走査する");
}

}  // namespace cv_demo


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
    cv_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    cv_test::benchmark();    // 速度計測
    cv_demo::demo();         // 使い方デモ
    return 0;
}

