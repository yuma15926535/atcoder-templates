// =============================================================================
//  合同（回転・平行移動（・裏返し）で重なるか。点の集合・多角形）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 06 点の集合
//  依存       : 00_int_common, 00_geo_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - 判定は整数で厳密（128 ビット）。回転の角は 90 度の倍数でなくてよい（(3, 4) / 5 のような回転も見つける）
//    - 重ね方：B[to[i]] = R(angle) · F(A[i]) + shift。F は flipped なら x 軸について折り返し (x, y) → (x, -y)、そうでなければそのまま
//    - 点の集合は同じ点を含んでもよい（同じ点どうしの対応は好きに選ぶ）。座標は 10^9 まで、点の数は 2×10^5 程度まで
//
//  機能一覧                                                                    計算量
//    congruence(A, B, allow_flip = false)      点の集合 A を回転・平行移動（allow_flip なら裏返しも）で B に重ねられるか  O(N log N)
//                                              Congruence{ok, flipped, to, angle, shift}
//    congruent(A, B, allow_flip = false)       同じ判定の bool だけ
//    polygon_congruence(P, Q, allow_flip = false)  多角形（頂点の順番つき。どちら回りでもよい・始点はどこでもよい）  O(N)
//                                              to[i] は P[i] が重なる Q の頂点の番号。隣り合う頂点は同じ点でないこと
//    polygon_congruent(P, Q, allow_flip = false)   同じ判定の bool だけ
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    cg_test::random_test()  … 2 点の対応から変換を作って全部確かめる愚直な実装と突き合わせるランダムテスト
//    cg_test::benchmark()    … 速度計測
//    cg_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_CONGRUENCE
#define GG_CONGRUENCE
// ---- 合同（幾何第二世代）  回転・平行移動（・裏返し）で重なるか。整数で厳密 ----
struct Congruence {
    bool ok = false;       // 重ねられるか
    bool flipped = false;  // 裏返したか（x 軸について折り返してから回す）
    std::vector<int> to;   // A[i] が B[to[i]] に重なる
    Real angle = 0;        // 回転の角（反時計回り、[0, 2π)）
    PointD shift;          // 平行移動（B[to[i]] = R(angle) · F(A[i]) + shift）
};

namespace geo_internal {
struct CgVec {  // 128 ビットの点
    Int128 x, y;
};
inline Int128 cg_cross(const CgVec& a, const CgVec& b) { return a.x * b.y - a.y * b.x; }
inline Int128 cg_dot(const CgVec& a, const CgVec& b) { return a.x * b.x + a.y * b.y; }
struct CgTok {  // 回転で変わらない量：長さの 2 乗と、次との外積・内積
    Int128 n, c, d;
    bool operator==(const CgTok& o) const { return n == o.n && c == o.c && d == o.d; }
};
// pat が text（巡回）のどこから始まるか（KMP）。見つからなければ -1
inline int cg_find_cyclic(const std::vector<CgTok>& pat, const std::vector<CgTok>& text) {
    const int m = (int)pat.size();
    if ((int)text.size() != m) return -1;
    if (m == 0) return 0;
    std::vector<int> f(m, 0);
    for (int i = 1, k = 0; i < m; ++i) {
        while (k > 0 && !(pat[i] == pat[k])) k = f[k - 1];
        if (pat[i] == pat[k]) ++k;
        f[i] = k;
    }
    for (int i = 0, k = 0; i < 2 * m - 1; ++i) {
        const CgTok& t = text[i % m];
        while (k > 0 && !(t == pat[k])) k = f[k - 1];
        if (t == pat[k]) ++k;
        if (k == m) return i - m + 1;
    }
    return -1;
}
// 重心を原点に移して n 倍した点を偏角（同じ向きは近い順）の順に並べる。原点の点は別に
inline void cg_prepare(const std::vector<Point>& ps, bool flip, std::vector<int>& order, std::vector<CgVec>& vs, std::vector<int>& zeros) {
    const long long n = (long long)ps.size();
    Int128 sx = 0, sy = 0;
    for (const Point& p : ps) sx += p.x, sy += flip ? -p.y : p.y;
    std::vector<CgVec> all(n);
    for (int i = 0; i < n; ++i) all[i] = CgVec{(Int128)ps[i].x * n - sx, (Int128)(flip ? -ps[i].y : ps[i].y) * n - sy};
    order.clear(), zeros.clear();
    for (int i = 0; i < n; ++i) (all[i].x == 0 && all[i].y == 0 ? zeros : order).push_back(i);
    auto half = [](const CgVec& v) { return (v.y < 0 || (v.y == 0 && v.x < 0)) ? 1 : 0; };
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        const CgVec &u = all[a], &v = all[b];
        const int hu = half(u), hv = half(v);
        if (hu != hv) return hu < hv;
        const Int128 c = cg_cross(u, v);
        if (c != 0) return c > 0;
        return cg_dot(u, u) < cg_dot(v, v);
    });
    vs.clear();
    for (int i : order) vs.push_back(all[i]);
}
inline std::vector<CgTok> cg_tokens(const std::vector<CgVec>& v) {
    const int m = (int)v.size();
    std::vector<CgTok> t(m);
    for (int k = 0; k < m; ++k) {
        const CgVec &a = v[k], &b = v[(k + 1) % m];
        t[k] = CgTok{cg_dot(a, a), cg_cross(a, b), cg_dot(a, b)};
    }
    return t;
}
inline void cg_transform(const std::vector<Point>& A, const std::vector<Point>& B, Congruence& r) {
    // 角は「重心から最も遠い点 A[i0] → B[to[i0]]」の向きの差から、重心を n 倍した整数のベクトルで誤差なく近く求める
    const int n = (int)A.size();
    if (n == 0) return;
    auto f = [&](const Point& p) { return r.flipped ? Point(p.x, -p.y) : p; };
    Int128 sx = 0, sy = 0, tx = 0, ty = 0;
    for (const Point& p : A) sx += f(p).x, sy += f(p).y;
    for (const Point& p : B) tx += p.x, ty += p.y;
    int i0 = 0;
    Int128 best = -1;
    for (int i = 0; i < n; ++i) {
        const Int128 ux = (Int128)f(A[i]).x * n - sx, uy = (Int128)f(A[i]).y * n - sy, m = ux * ux + uy * uy;
        if (m > best) best = m, i0 = i;
    }
    const Int128 ux = (Int128)f(A[i0]).x * n - sx, uy = (Int128)f(A[i0]).y * n - sy;
    const Int128 vx = (Int128)B[r.to[i0]].x * n - tx, vy = (Int128)B[r.to[i0]].y * n - ty;
    Real t = 0;
    if (best > 0) {
        t = std::atan2((Real)(ux * vy - uy * vx), (Real)(ux * vx + uy * vy));
        if (t < 0) t += 2 * GEO_PI;
        if (t >= 2 * GEO_PI) t = 0;
    }
    const PointD ga((Real)sx / n, (Real)sy / n), gb((Real)tx / n, (Real)ty / n);
    r.angle = t;
    r.shift = gb - rotated(ga, t);
}
}  // namespace geo_internal

// 点の集合 A を回転・平行移動（allow_flip なら裏返しも）で B に重ねられるか。重ね方も返す
inline Congruence congruence(const std::vector<Point>& A, const std::vector<Point>& B, bool allow_flip = false) {
    Congruence r;
    if (A.size() != B.size()) return r;
    const int n = (int)A.size();
    std::vector<int> ob, zb, oa, za;
    std::vector<geo_internal::CgVec> vb, va;
    geo_internal::cg_prepare(B, false, ob, vb, zb);
    const std::vector<geo_internal::CgTok> tb = geo_internal::cg_tokens(vb);
    for (int f = 0; f < (allow_flip ? 2 : 1); ++f) {
        geo_internal::cg_prepare(A, f == 1, oa, va, za);
        if (oa.size() != ob.size()) return r;
        const int sh = geo_internal::cg_find_cyclic(geo_internal::cg_tokens(va), tb);
        if (sh < 0) continue;
        const int m = (int)oa.size();
        r.ok = true, r.flipped = f == 1, r.to.assign(n, -1);
        for (int k = 0; k < m; ++k) r.to[oa[k]] = ob[(k + sh) % m];
        for (size_t k = 0; k < za.size(); ++k) r.to[za[k]] = zb[k];
        geo_internal::cg_transform(A, B, r);
        return r;
    }
    return r;
}
inline bool congruent(const std::vector<Point>& A, const std::vector<Point>& B, bool allow_flip = false) { return congruence(A, B, allow_flip).ok; }

// 多角形 P（頂点の順番つき）を回転・平行移動（allow_flip なら裏返しも）で Q に重ねられるか。
// どちら回りに並んでいてもよく、始点もどこでもよい。to[i] は P[i] が重なる Q の頂点の番号
inline Congruence polygon_congruence(const std::vector<Point>& P, const std::vector<Point>& Q, bool allow_flip = false) {
    Congruence r;
    if (P.size() != Q.size()) return r;
    const int n = (int)P.size();
    auto edges = [&](const std::vector<Point>& v, bool flip, bool rev) {  // 辺のベクトルの列（裏返し・逆回り）
        std::vector<geo_internal::CgVec> e(n);
        for (int k = 0; k < n; ++k) {
            const int i = rev ? n - 1 - k : k, j = rev ? (2 * n - 2 - k) % n : (k + 1) % n;
            const Point d = v[j] - v[i];
            e[k] = geo_internal::CgVec{d.x, flip ? -d.y : d.y};
        }
        return geo_internal::cg_tokens(e);
    };
    if (n == 0) {
        r.ok = true;
        return r;
    }
    for (int f = 0; f < (allow_flip ? 2 : 1); ++f) {
        const std::vector<geo_internal::CgTok> tp = edges(P, f == 1, false);
        for (int rev = 0; rev < 2; ++rev) {
            const int sh = geo_internal::cg_find_cyclic(tp, edges(Q, false, rev == 1));
            if (sh < 0) continue;
            r.ok = true, r.flipped = f == 1, r.to.assign(n, -1);
            for (int k = 0; k < n; ++k) r.to[k] = rev ? (n - 1 - (k + sh) % n) : (k + sh) % n;
            geo_internal::cg_transform(P, Q, r);
            return r;
        }
    }
    return r;
}
inline bool polygon_congruent(const std::vector<Point>& P, const std::vector<Point>& Q, bool allow_flip = false) { return polygon_congruence(P, Q, allow_flip).ok; }
#endif  // GG_CONGRUENCE
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
namespace cg_test {

#define CG_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[congruence] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK CG_CHECK

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

// 愚直：A の 2 点の行き先を全部試し、変換（複素数の比。有理数で厳密）で全部が B に重なるか
static bool brute_map(const std::vector<Point>& A, const std::vector<Point>& B, bool flip) {
    const int n = (int)A.size();
    if (n != (int)B.size()) return false;
    if (n == 0) return true;
    std::vector<Point> a = A;
    if (flip)
        for (Point& p : a) p.y = -p.y;
    int j = -1;
    for (int k = 1; k < n; ++k)
        if (a[k] != a[0]) { j = k; break; }
    std::multiset<Point> sb(B.begin(), B.end());
    if (j < 0) {  // 全部が同じ点
        for (const Point& q : B)
            if (q != B[0]) return false;
        return true;
    }
    const Point u = a[j] - a[0];
    const Int128 uu = norm2_128(u);
    for (int k = 0; k < n; ++k)
        for (int l = 0; l < n; ++l) {
            const Point w = B[l] - B[k];
            if (norm2_128(w) != uu) continue;
            // p → B[k] + (p - a0) · w · conj(u) / |u|²
            std::multiset<Point> s;
            bool ok = true;
            for (const Point& p : a) {
                const Point d = p - a[0];
                const Int128 re = ((Int128)d.x * w.x - (Int128)d.y * w.y) * u.x + ((Int128)d.x * w.y + (Int128)d.y * w.x) * u.y;
                const Int128 im = ((Int128)d.x * w.y + (Int128)d.y * w.x) * u.x - ((Int128)d.x * w.x - (Int128)d.y * w.y) * u.y;
                if (re % uu != 0 || im % uu != 0) { ok = false; break; }
                s.insert(Point((Coord)(re / uu) + B[k].x, (Coord)(im / uu) + B[k].y));
            }
            if (ok && s == sb) return true;
        }
    return false;
}
// 結果の重ね方が正しいか：to が置換で、距離が全部同じ・向き（外積の符号）が同じ（裏返しなら逆）・変換で重なる
static bool check_map(const std::vector<Point>& A, const std::vector<Point>& B, const Congruence& r) {
    const int n = (int)A.size();
    std::vector<int> s = r.to;
    std::sort(s.begin(), s.end());
    for (int i = 0; i < n; ++i)
        if (s[i] != i) return false;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (dist2(A[i], A[j]) != dist2(B[r.to[i]], B[r.to[j]])) return false;
            for (int k = 0; k < std::min(n, 6); ++k) {
                const Int128 c1 = cross128(A[j] - A[i], A[k] - A[i]), c2 = cross128(B[r.to[j]] - B[r.to[i]], B[r.to[k]] - B[r.to[i]]);
                if (c1 != (r.flipped ? -c2 : c2)) return false;
            }
        }
    Coord C = 1;
    for (const Point& p : B) C = std::max({C, std::llabs(p.x), std::llabs(p.y)});
    for (int i = 0; i < n; ++i) {
        const Point f = r.flipped ? Point(A[i].x, -A[i].y) : A[i];
        if (!nearp(rotated(f, r.angle) + r.shift, to_d(B[r.to[i]]), 1e-12L * (Real)n)) return false;
    }
    return true;
}
static std::vector<Point> rnd_move(std::mt19937_64& rng, std::vector<Point> A, long long C, bool& flipped) {
    flipped = rng() % 2;
    if (flipped)
        for (Point& p : A) p.y = -p.y;
    const long long k = rnd(rng, 0, 3);
    const Point t = rndp(rng, C);
    for (Point& p : A) p = rot90(p, Point(), k) + t;
    std::shuffle(A.begin(), A.end(), rng);
    return A;
}

static void test_points(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 7);
        const long long C = it % 3 == 0 ? 1 : 3;
        std::vector<Point> A = rndps(rng, n, C), B;
        const int t = it % 4;
        bool fl = false;
        if (t == 0) B = rndps(rng, n, C);                      // ばらばら
        else if (t == 1) B = rnd_move(rng, A, 5, fl);          // 90 度の倍数の回転・裏返し
        else if (t == 2) {                                     // (3, 4) / 5 の回転（A を 5 倍しておく）
            for (Point& p : A) p *= 5;
            B = A;
            for (Point& p : B) p = complex_mul(p, Point(3, 4)) / 5 + Point(1, -2);
            if (rng() % 2)
                for (Point& p : B) p.y = -p.y, fl = true;
            std::shuffle(B.begin(), B.end(), rng);
        } else {  // 1 点だけ少しずらす
            B = rnd_move(rng, A, 5, fl);
            if (n > 0) B[rng() % n] += Point(rnd(rng, -1, 1), rnd(rng, -1, 1));
        }
        for (int f = 0; f < 2; ++f) {
            const bool e = brute_map(A, B, false) || (f == 1 && brute_map(A, B, true));
            const Congruence r = congruence(A, B, f == 1);
            CHECK(r.ok == e && congruent(A, B, f == 1) == e);
            if (r.ok) CHECK(check_map(A, B, r) && (f == 1 || !r.flipped));
        }
    }
    // 大きい座標・多い点
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 60);
        std::vector<Point> A = rndps(rng, n, 100000000LL);
        bool fl;
        std::vector<Point> B = rnd_move(rng, A, 800000000LL, fl);
        Congruence r = congruence(A, B, true);
        CHECK(r.ok && check_map(A, B, r));
        CHECK(congruent(A, B, false) == (!fl || brute_map(A, B, false)));
        B[0] += Point(1, 0);
        CHECK(!congruent(A, B, true) || brute_map(A, B, false) || brute_map(A, B, true));
    }
}

static void test_polygons(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 7);
        const long long C = it % 3 == 0 ? 1 : 3;
        std::vector<Point> P = rndps(rng, n, C);
        for (int k = 0; k < n; ++k)  // 隣り合う頂点は違う点に
            if (n > 1 && P[k] == P[(k + 1) % n]) P[k] += Point(C + 1, 0);
        std::vector<Point> Q = P;
        const int t = it % 3;
        if (t == 0) Q = rndps(rng, n, C);
        else {
            if (rng() % 2)
                for (Point& p : Q) p.y = -p.y;
            const long long k = rnd(rng, 0, 3);
            for (Point& p : Q) p = rot90(p, Point(), k) + Point(2, -1);
            std::rotate(Q.begin(), Q.begin() + rng() % n, Q.end());
            if (rng() % 2) std::reverse(Q.begin(), Q.end());
            if (t == 2) Q[rng() % n] += Point(0, 1);
        }
        bool okq = true;
        for (int k = 0; k < n; ++k) okq &= n == 1 || Q[k] != Q[(k + 1) % n];
        if (!okq) continue;
        for (int f = 0; f < 2; ++f) {
            // 愚直：始点と向き（と裏返し）を全部試し、頂点の対応で brute_map と同じ変換を確かめる
            bool e = false;
            for (int fl = 0; fl < (f == 1 ? 2 : 1) && !e; ++fl)
                for (int s = 0; s < n && !e; ++s)
                    for (int rev = 0; rev < 2 && !e; ++rev) {
                        std::vector<Point> q2(n);
                        for (int k = 0; k < n; ++k) q2[k] = Q[rev ? ((s - k) % n + n) % n : (s + k) % n];
                        // P[k] → q2[k] の対応で重なるか：距離・外積がすべて同じ
                        bool ok = true;
                        for (int i = 0; i < n && ok; ++i)
                            for (int j = 0; j < n && ok; ++j) {
                                ok &= dist2(P[i], P[j]) == dist2(q2[i], q2[j]);
                                for (int k = 0; k < n && ok; ++k) {
                                    const Int128 c1 = cross128(P[j] - P[i], P[k] - P[i]), c2 = cross128(q2[j] - q2[i], q2[k] - q2[i]);
                                    ok &= c1 == (fl ? -c2 : c2);
                                }
                            }
                        e = ok;
                    }
            const Congruence r = polygon_congruence(P, Q, f == 1);
            CHECK(r.ok == e && polygon_congruent(P, Q, f == 1) == e);
            if (r.ok) {
                if (!check_map(P, Q, r)) {
                    std::cerr << "P:"; for (auto& p : P) std::cerr << " " << p; std::cerr << "\nQ:"; for (auto& q : Q) std::cerr << " " << q;
                    std::cerr << "\nto:"; for (int x : r.to) std::cerr << " " << x; std::cerr << " flipped=" << r.flipped << " f=" << f << "\n";
                }
                CHECK(check_map(P, Q, r));
                // 隣り合う頂点は隣り合う頂点に
                for (int k = 0; k < n; ++k) {
                    const int d = ((r.to[(k + 1) % n] - r.to[k]) % n + n) % n;
                    CHECK(n <= 2 || d == 1 || d == n - 1);
                }
            }
        }
    }
}

// 重ね方の精度：重心のすぐ近くの点しかない・座標が 10^9 近く・点が多い
static void test_precision(std::mt19937_64& rng, int& cases) {
    auto max_err = [](const std::vector<Point>& A, const std::vector<Point>& B, const Congruence& r) {
        Real e = 0;
        for (size_t i = 0; i < A.size(); ++i) {
            const Point f = r.flipped ? Point(A[i].x, -A[i].y) : A[i];
            const PointD q = rotated(f, r.angle) + r.shift;
            e = std::max({e, std::fabs(q.x - (Real)B[r.to[i]].x), std::fabs(q.y - (Real)B[r.to[i]].y)});
        }
        return e;
    };
    {  // 3 点・座標 10^9 近く（90 度回して平行移動）
        const std::vector<Point> A = {{443501150, 443501150}, {943168399, -53790020}, {-56166098, 940792322}};
        std::vector<Point> B;
        for (const Point& p : A) B.push_back(rot90(p) + Point(123456789, -387654321));
        const Congruence r = congruence(A, B);
        CHECK(r.ok && std::fabs(r.angle - GEO_PI / 2) <= TOL && max_err(A, B, r) <= TOL * 1e9L * 16);
        ++cases;
    }
    for (int it = 0; it < 40 * TEST_SCALE; ++it, ++cases) {  // (c, c) について対称な組 + 重心から 1 だけずれた点
        const Coord c = 476543217 + rnd(rng, -1000, 1000);
        const int m = it % 2 ? 10 : 1000;
        std::vector<Point> A = {Point(c, c), Point(c + 1, c + 2)};
        for (int k = 0; k < m; ++k) {
            const Point d(rnd(rng, -500000000, 500000000), rnd(rng, -500000000, 500000000));
            A.push_back(Point(c, c) + d), A.push_back(Point(c, c) - d);
        }
        const long long k = rnd(rng, 0, 3);
        std::vector<Point> B;
        for (const Point& p : A) B.push_back(rot90(p, Point(), k) + Point(-7, 11));
        std::shuffle(B.begin(), B.end(), rng);
        const Congruence r = congruence(A, B);
        CHECK(r.ok && max_err(A, B, r) <= TOL * 1e9L * 16);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_precision(rng, cases);
    test_points(rng, cases);
    test_polygons(rng, cases);
    std::cout << "[congruence] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[congruence] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 200000;
    const std::vector<Point> A = rndps(rng, n, 1000000000LL / 5);
    std::vector<Point> B = A;
    for (Point& p : B) p = complex_mul(p, Point(3, 4)) / 1 + Point(7, 7);  // 5 倍の大きさ（合同ではない）
    std::vector<Point> C = A;
    for (Point& p : C) p = rot90(p, Point(1, 1), 1);
    std::shuffle(C.begin(), C.end(), rng);
    auto t0 = std::chrono::steady_clock::now();
    const bool r1 = congruent(A, B, true), r2 = congruent(A, C);
    const long long t1 = ms_since(t0);
    std::vector<Point> P = A, Q = A;
    std::rotate(Q.begin(), Q.begin() + 12345, Q.end());
    t0 = std::chrono::steady_clock::now();
    const bool r3 = polygon_congruent(P, Q);
    const long long t2 = ms_since(t0);
    std::cout << "[congruence] benchmark: N = " << n << "  congruence × 2 " << t1 << " ms, polygon_congruence " << t2 << " ms  (check " << r1 << r2 << r3 << ")\n";
}

}  // namespace cg_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace cg_demo {

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
#define SHOW(expr, comment) cg_demo::show(#expr, cg_demo::to_s(expr), comment)
#define RUN(...)                   \
    cg_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("合同 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  congruence(A, B, allow_flip = false) / congruent(...)                 点の集合（順番なし）\n"
              << "  polygon_congruence(P, Q, allow_flip = false) / polygon_congruent(...)  多角形（頂点の順番つき）\n"
              << "  結果 Congruence{ok, flipped, to, angle, shift}：B[to[i]] = R(angle)·F(A[i]) + shift\n";

    api("【1】 congruence(A, B, allow_flip) / congruent(A, B, allow_flip)", "Congruence / bool",
        "点の集合 A を回転・平行移動（allow_flip なら裏返しも）して B に重ねられるか。整数で厳密", "O(N log N)");
    RUN(std::vector<Point> A = {{0, 0}, {2, 0}, {0, 1}});
    RUN(std::vector<Point> B = {{5, 5}, {5, 7}, {4, 5}});
    RUN(Congruence r = congruence(A, B));
    SHOW(r.ok, "90 度回して (5, 5) へ");
    SHOW(r.to, "A[i] → B[to[i]]");
    SHOW(rad_to_deg(r.angle), "回転の角（度）");
    SHOW(r.shift, "平行移動");
    SHOW(congruent(A, std::vector<Point>{{0, 0}, {2, 0}, {0, -1}}), "裏返しは不可（allow_flip = false）");
    SHOW(congruent(A, std::vector<Point>{{0, 0}, {2, 0}, {0, -1}}, true), "allow_flip = true なら重なる");
    SHOW(congruent(std::vector<Point>{{0, 0}, {5, 0}}, std::vector<Point>{{1, 1}, {4, 5}}), "(3, 4) の向きへの回転も見つける");
    note("ABC207 D（Congruence Points）: congruent(S, T) で Yes / No（回転の角は 90 度の倍数とは限らない）");
    note("考え方：重心を原点に移し、偏角の順に並べて（長さの 2 乗, 次との外積, 内積）の列を作る → 巡回して一致するか（KMP）");

    api("【2】 polygon_congruence(P, Q, allow_flip) / polygon_congruent(P, Q, allow_flip)", "Congruence / bool",
        "多角形として重なるか（頂点の順番を保つ。始点・回り方は自由）。to[i] は P[i] が重なる Q の頂点", "O(N)");
    RUN(std::vector<Point> P = {{0, 0}, {3, 0}, {3, 1}, {0, 2}});
    RUN(std::vector<Point> Q = {{10, 3}, {10, 0}, {8, 0}, {9, 3}});  // P を 90 度回して (10, 0) だけ動かし、始点と回り方を変えた
    SHOW(polygon_congruent(P, Q), "");
    SHOW(polygon_congruence(P, Q).to, "P[i] → Q[to[i]]（(0,0) → Q[1] = (10,0) など）");
    SHOW(polygon_congruent(P, std::vector<Point>{{0, 0}, {3, 0}, {3, 2}, {0, 1}}), "裏返し（台形の向きが逆）は不可");
    SHOW(polygon_congruent(P, std::vector<Point>{{0, 0}, {3, 0}, {3, 2}, {0, 1}}, true), "allow_flip = true なら重なる");
    note("点の集合として同じでも、頂点の順番（辺のつながり）が違えば多角形としては別");
}

}  // namespace cg_demo


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
    cg_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    cg_test::benchmark();    // 速度計測
    cg_demo::demo();         // 使い方デモ
    return 0;
}

