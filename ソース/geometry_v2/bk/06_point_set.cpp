// =============================================================================
//  点の集合（最近点対・距離 D 以内の組・一直線・三角形・直線の数・角の分類・最大の角・長方形 / 正方形 / 平行四辺形 / 台形の数）
//                                                                                    --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 06 点の集合
//  依存       : 00_geo_base, 02_line, 01_arg_sort（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - 点は vector<Point>。番号は入力の順（0-indexed）。すべて整数で厳密（|x|, |y| ≤ 10^9）
//    - 三角形・直線の数え上げは「番号の組」で数える（同じ点が 2 つあれば別の点として数える。三角形なら面積 0 なので数えない）
//    - 長方形・正方形・平行四辺形・台形の数は、同じ点を 1 つにまとめてから数える
//
//  機能一覧                                                                    計算量
//    closest_pair(ps)                  最も近い 2 点 ClosestPair{i, j, d2}（i < j。2 点未満なら -1）   O(N log N)
//    pairs_within(ps, D) / pairs_within_d2(ps, D2)   距離 ≤ D（距離² ≤ D2）の組 (i, j)（i < j）の一覧   O(N log N + 答えの数)
//    count_collinear_triples(ps)       一直線に並ぶ 3 点の組の数                                        O(N^2 log N)
//    count_triangles(ps)               面積が正の三角形の数（= C(N, 3) - 一直線の組）                   O(N^2 log N)
//    max_collinear(ps)                 1 本の直線の上に並ぶ最大の点の数 MaxCollinear{count, i, j}       O(N^2 log N)
//    lines_through(ps)                 2 点以上を通る直線ごとの点の数 vector<LineCount{key, count}>      O(N^2 log N)
//    count_lines_with_at_least(ps, K)  K 点以上を通る直線の数（無限にあるなら -1）                     O(N^2 log N)
//    count_triangles_by_angle(ps)      鋭角・直角・鈍角三角形の数 TriangleCounts{acute, right, obtuse}  O(N^2 log N)
//    max_angle(ps)                     最大の角 ∠p_i p_j p_k（j が頂点）MaxAngle{angle（ラジアン）, i, j, k}  O(N^2 log N)
//    count_axis_rectangles(ps)         4 頂点が全部ある、辺が軸に平行な長方形の数                        O(Σ 列の点の数²)
//    count_squares(ps)                 4 頂点が全部ある正方形の数（傾いてもよい）                        O(N^2)
//    count_parallelograms(ps)          4 頂点が全部ある平行四辺形の数（つぶれたものは除く）              O(N^2 log N)
//    count_trapezoids(ps)              4 頂点が全部ある台形の数（1 組以上の対辺が平行。平行四辺形も含む） O(N^2 log N)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    ps_test::random_test()  … 全探索（全部の組・3 点・4 点の組）と突き合わせるランダムテスト
//    ps_test::benchmark()    … 速度計測
//    ps_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_GEO_BASE
#define GG_GEO_BASE
// ---- 幾何基盤（幾何第二世代）  点 Point（整数）/ PointD（実数）、直線・線分・円、向き、回転、角度、入力 ----
using Coord = long long;    // 整数の座標（|x|, |y| ≤ 10^9 まで正しい）
using Real = long double;   // 実数（double に変えてもよい）
__extension__ typedef __int128 Int128;  // 途中の計算用（出力は int128_to_string）
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

// ---- 入力 ----
// 小数の文字列を 10^k 倍した整数にする（"-1.5", 2 → -150、"3", 2 → 300、".25", 2 → 25）。誤差なし
// 小数点以下に k 桁を超える 0 でない桁があったら止める（黙って切り捨てないため）。結果は long long に収まる前提
inline long long parse_fixed(const std::string& s, int k) {
    size_t i = 0;
    bool neg = false, any = false;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) neg = (s[i] == '-'), ++i;
    long long v = 0;
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
    }
    if (!any) { std::cerr << "parse_fixed: 数字がない: " << s << "\n"; std::abort(); }
    for (int f = std::max(frac, 0); f < k; ++f) v *= 10;
    return neg ? -v : v;
}
inline long long read_fixed(int k) {
    std::string s;
    std::cin >> s;
    return parse_fixed(s, k);
}
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

inline std::string int128_to_string(Int128 v) {
    if (v == 0) return "0";
    using U = geo_internal::UInt128;
    U u = v < 0 ? (U)0 - (U)v : (U)v;
    std::string s;
    while (u) s += char('0' + (int)(u % 10)), u /= 10;
    if (v < 0) s += '-';
    std::reverse(s.begin(), s.end());
    return s;
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
#ifndef GG_ARG_SORT
#define GG_ARG_SORT
// ---- 偏角（幾何第二世代）  偏角の比較・並べ替え・向きの範囲の数え上げ・角の大きさの比較（すべて整数で厳密） ----
// 偏角は x 軸の正の向きから反時計回りに [0, 2π)。0 ベクトルはいちばん小さい
inline int arg_half(const Point& v) {  // 0: [0, π) / 1: [π, 2π) / -1: 0 ベクトル
    if (v.x == 0 && v.y == 0) return -1;
    return (v.y > 0 || (v.y == 0 && v.x > 0)) ? 0 : 1;
}
inline int arg_cmp(const Point& a, const Point& b) {  // -1: a の偏角が小さい / 0: 同じ向き / 1: a が大きい
    const int ha = arg_half(a), hb = arg_half(b);
    if (ha != hb) return ha < hb ? -1 : 1;
    if (ha < 0) return 0;
    const Coord c = cross(a, b);
    return c > 0 ? -1 : (c < 0 ? 1 : 0);
}
inline bool arg_less(const Point& a, const Point& b) { return arg_cmp(a, b) < 0; }
inline bool same_direction(const Point& a, const Point& b) { return arg_cmp(a, b) == 0; }

// o から見た偏角の順に並べた番号。o と同じ点が先頭、同じ向きは o に近い順（far_first = true なら遠い順）、さらに番号の順
inline std::vector<int> arg_sort_idx(const std::vector<Point>& ps, const Point& o = Point(), bool far_first = false) {
    const int n = (int)ps.size();
    std::vector<int> idx(n);
    std::vector<Point> v(n);
    std::vector<signed char> h(n);
    std::vector<Coord> nr(n);
    for (int i = 0; i < n; ++i) idx[i] = i, v[i] = ps[i] - o, h[i] = (signed char)arg_half(v[i]), nr[i] = norm2(v[i]);
    std::sort(idx.begin(), idx.end(), [&](int i, int j) {
        if (h[i] != h[j]) return h[i] < h[j];
        if (h[i] >= 0) {
            const Coord c = cross(v[i], v[j]);
            if (c != 0) return c > 0;
        }
        if (nr[i] != nr[j]) return far_first ? nr[i] > nr[j] : nr[i] < nr[j];
        return i < j;
    });
    return idx;
}
inline std::vector<Point> arg_sort(const std::vector<Point>& ps, const Point& o = Point(), bool far_first = false) {
    std::vector<Point> r;
    r.reserve(ps.size());
    for (int i : arg_sort_idx(ps, o, far_first)) r.push_back(ps[i]);
    return r;
}
// 向き start（0 ベクトル不可）から反時計回りに一周する順（偏角が start 以上の点が先、start より小さい点が後）。o と同じ点は先頭
inline std::vector<int> arg_sort_idx_from(const std::vector<Point>& ps, const Point& o, const Point& start, bool far_first = false) {
    std::vector<int> idx = arg_sort_idx(ps, o, far_first);
    size_t z = 0;
    while (z < idx.size() && ps[idx[z]] == o) ++z;
    size_t k = z;
    while (k < idx.size() && arg_cmp(ps[idx[k]] - o, start) < 0) ++k;
    std::rotate(idx.begin() + (long)z, idx.begin() + (long)k, idx.end());
    return idx;
}
inline std::vector<Point> arg_sort_from(const std::vector<Point>& ps, const Point& o, const Point& start, bool far_first = false) {
    std::vector<Point> r;
    r.reserve(ps.size());
    for (int i : arg_sort_idx_from(ps, o, start, far_first)) r.push_back(ps[i]);
    return r;
}
// 向きごとに 1 点だけ残す（keep_far = false なら o に近い点、true なら遠い点）。偏角の順。o と同じ点は除く
inline std::vector<int> unique_directions_idx(const std::vector<Point>& ps, const Point& o = Point(), bool keep_far = false) {
    std::vector<int> r;
    Point last;
    for (int i : arg_sort_idx(ps, o, keep_far)) {
        const Point v = ps[i] - o;
        if (v == Point()) continue;
        if (!r.empty() && arg_cmp(last, v) == 0) continue;
        r.push_back(i), last = v;
    }
    return r;
}
inline std::vector<Point> unique_directions(const std::vector<Point>& ps, const Point& o = Point(), bool keep_far = false) {
    std::vector<Point> r;
    for (int i : unique_directions_idx(ps, o, keep_far)) r.push_back(ps[i]);
    return r;
}

// 点を o から見た偏角の順に並べ、向きの範囲に入る点を二分探索で数える
struct ArgSorted {
    Point o;
    std::vector<int> order;  // 方向のある点の番号（偏角の順。同じ向きは o に近い順）
    std::vector<Point> vec;  // vec[k] = ps[order[k]] - o
    std::vector<int> zeros;  // o と同じ点の番号
    ArgSorted() {}
    explicit ArgSorted(const std::vector<Point>& ps, const Point& o_ = Point()) : o(o_) {
        for (int i : arg_sort_idx(ps, o)) {
            const Point v = ps[i] - o;
            if (v == Point()) zeros.push_back(i);
            else order.push_back(i), vec.push_back(v);
        }
    }
    int size() const { return (int)order.size(); }
    int lower_bound(const Point& dir) const {  // 偏角が dir 以上の最初の位置（0 ~ size()）
        int lo = 0, hi = size();
        while (lo < hi) {
            const int m = (lo + hi) / 2;
            if (arg_cmp(vec[m], dir) < 0) lo = m + 1;
            else hi = m;
        }
        return lo;
    }
    int upper_bound(const Point& dir) const {  // 偏角が dir より大きい最初の位置
        int lo = 0, hi = size();
        while (lo < hi) {
            const int m = (lo + hi) / 2;
            if (arg_cmp(vec[m], dir) <= 0) lo = m + 1;
            else hi = m;
        }
        return lo;
    }
    // 向き l から反時計回りに向き r までの範囲（2π をまたいでよい）にある点の数。端を含むかは include_l / include_r
    // l と r が同じ向きなら、両端を含むときはその向きの点の数、それ以外は 0（一周ぶんは size()）。l, r は 0 ベクトル不可
    int count(const Point& l, const Point& r, bool include_l = true, bool include_r = false) const {
        const int c = arg_cmp(l, r);
        if (c == 0) return include_l && include_r ? upper_bound(l) - lower_bound(l) : 0;
        const int a = include_l ? lower_bound(l) : upper_bound(l);
        const int b = include_r ? upper_bound(r) : lower_bound(r);
        if (c < 0) return std::max(0, b - a);
        return (size() - a) + b;
    }
    bool exists(const Point& l, const Point& r, bool include_l = true, bool include_r = false) const { return count(l, r, include_l, include_r) > 0; }
    std::vector<std::vector<int>> groups() const {  // 同じ向きの点の番号のまとまり（偏角の順。中は o に近い順）
        std::vector<std::vector<int>> g;
        for (int k = 0; k < size(); ++k) {
            if (k == 0 || arg_cmp(vec[k - 1], vec[k]) != 0) g.emplace_back();
            g.back().push_back(order[k]);
        }
        return g;
    }
};

// ∠(a, b)（ベクトル a と b のなす角 [0, π]）と ∠(c, d) を厳密に比べる。-1: ∠(a, b) が小さい / 0: 等しい / 1: 大きい
// 0 ベクトルを含む角は 0 とする。∠AOB と ∠COD を比べるなら angle_cmp(A - O, B - O, C - O, D - O)
inline int angle_cmp(const Point& a, const Point& b, const Point& c, const Point& d) {
    auto key = [](const Point& u, const Point& v, Int128& x, Int128& y) {  // なす角 = (内積, |外積|) の偏角
        if (u == Point() || v == Point()) { x = 1, y = 0; return; }
        x = dot128(u, v), y = cross128(u, v);
        if (y < 0) y = -y;
    };
    Int128 x1, y1, x2, y2;
    key(a, b, x1, y1), key(c, d, x2, y2);
    const int h1 = (y1 == 0 && x1 < 0), h2 = (y2 == 0 && x2 < 0);  // 1: ちょうど π
    if (h1 != h2) return h1 < h2 ? -1 : 1;
    if (h1) return 0;
    const Int128 cr = x1 * y2 - y1 * x2;
    return cr > 0 ? -1 : (cr < 0 ? 1 : 0);
}
inline bool angle_less(const Point& a, const Point& b, const Point& c, const Point& d) { return angle_cmp(a, b, c, d) < 0; }

// 各点 p について、o から見た偏角が (θp, θp + π) にある点の数（= cross(p - o, q - o) > 0 の q の数）。o と同じ点は 0
inline std::vector<int> count_in_open_half(const std::vector<Point>& ps, const Point& o = Point()) {
    const ArgSorted S(ps, o);
    std::vector<int> res(ps.size(), 0);
    for (int k = 0; k < S.size(); ++k) res[S.order[k]] = S.count(S.vec[k], -S.vec[k], false, false);
    return res;
}
#endif  // GG_ARG_SORT
#ifndef GG_POINT_SET
#define GG_POINT_SET
// ---- 点の集合（幾何第二世代）  最近点対・距離 D 以内の組・一直線・三角形・直線・角・長方形などの数え上げ（整数で厳密） ----
struct ClosestPair {  // 最も近い 2 点の番号（i < j）と距離の 2 乗。2 点未満なら -1
    int i = -1, j = -1;
    Coord d2 = -1;
};
struct MaxCollinear {  // 1 本の直線の上に並ぶ最大の点の数と、その直線の上の異なる 2 点（全部が同じ点なら j = -1）
    int count = 0, i = -1, j = -1;
};
struct LineCount {  // 直線（LineKey）と、その上にある点の数
    LineKey key;
    int count = 0;
};
struct TriangleCounts {  // 鋭角・直角・鈍角三角形の数
    long long acute = 0, right = 0, obtuse = 0;
};
struct MaxAngle {  // 最大の角 ∠p_i p_j p_k（j が頂点、ラジアン [0, π]）。3 点がとれなければ angle = -1
    Real angle = -1;
    int i = -1, j = -1, k = -1;
};

namespace geo_internal {
// 点 i から、番号が i より大きい点への向きを「直線の向き」（v と -v を同じ）でまとめる。zero = i と同じ点の数
// f(まとまりの点の数, まとまりの最初の点の番号) をまとまりごとに呼ぶ
template <class F>
void line_groups(const std::vector<Point>& ps, int i, std::vector<std::pair<Point, int>>& buf, int& zero, F&& f) {
    const int n = (int)ps.size();
    buf.clear(), zero = 0;
    for (int j = i + 1; j < n; ++j) {
        Point v = ps[j] - ps[i];
        if (v == Point()) { ++zero; continue; }
        if (v.y < 0 || (v.y == 0 && v.x < 0)) v = -v;  // 偏角を [0, π) にそろえる
        buf.push_back({v, j});
    }
    std::sort(buf.begin(), buf.end(), [](const std::pair<Point, int>& a, const std::pair<Point, int>& b) { return cross(a.first, b.first) > 0; });
    for (size_t k = 0; k < buf.size();) {
        size_t l = k;
        while (l < buf.size() && cross(buf[k].first, buf[l].first) == 0) ++l;
        f((int)(l - k), buf[k].second);
        k = l;
    }
}
inline std::vector<Point> unique_points(std::vector<Point> s) {
    std::sort(s.begin(), s.end());
    s.erase(std::unique(s.begin(), s.end()), s.end());
    return s;
}
inline Point upper_dir(Point v) {  // v と -v のうち偏角が [0, π) の方
    if (v.y < 0 || (v.y == 0 && v.x < 0)) v = -v;
    return v;
}
inline long long choose2(long long c) { return c * (c - 1) / 2; }
}  // namespace geo_internal

// 最も近い 2 点（分割統治。再帰の深さは log N）
inline ClosestPair closest_pair(const std::vector<Point>& ps) {
    const int n = (int)ps.size();
    if (n < 2) return ClosestPair();
    std::vector<int> id(n), buf(n), strip(n);
    for (int i = 0; i < n; ++i) id[i] = i;
    std::sort(id.begin(), id.end(), [&](int a, int b) { return ps[a] < ps[b]; });
    ClosestPair best{std::min(id[0], id[1]), std::max(id[0], id[1]), dist2(ps[id[0]], ps[id[1]])};
    auto upd = [&](int a, int b) {
        const Coord d = dist2(ps[a], ps[b]);
        if (d < best.d2) best = ClosestPair{std::min(a, b), std::max(a, b), d};
    };
    auto by_y = [&](int a, int b) { return ps[a].y < ps[b].y; };
    auto rec = [&](auto&& self, int l, int r) -> void {  // id[l, r) は x の順 → 終わると y の順
        if (r - l <= 3) {
            for (int a = l; a < r; ++a)
                for (int b = a + 1; b < r; ++b) upd(id[a], id[b]);
            std::sort(id.begin() + l, id.begin() + r, by_y);
            return;
        }
        const int m = (l + r) / 2;
        const Coord midx = ps[id[m]].x;
        self(self, l, m), self(self, m, r);
        std::merge(id.begin() + l, id.begin() + m, id.begin() + m, id.begin() + r, buf.begin() + l, by_y);
        std::copy(buf.begin() + l, buf.begin() + r, id.begin() + l);
        int sz = 0;
        for (int k = l; k < r; ++k) {
            const Point& p = ps[id[k]];
            const Coord dx = p.x - midx;
            if (dx * dx >= best.d2) continue;
            for (int t = sz - 1; t >= 0; --t) {
                const Coord dy = p.y - ps[strip[t]].y;
                if (dy * dy >= best.d2) break;
                upd(id[k], strip[t]);
            }
            strip[sz++] = id[k];
        }
    };
    rec(rec, 0, n);
    return best;
}
// 距離の 2 乗が D2 以下の組 (i, j)（i < j）をすべて（番号の順）。マス分け（1 辺 ⌈√D2⌉ のマスと隣のマスだけ調べる）
inline std::vector<std::pair<int, int>> pairs_within_d2(const std::vector<Point>& ps, Coord D2) {
    const int n = (int)ps.size();
    std::vector<std::pair<int, int>> res;
    if (D2 < 0) return res;
    const Coord s = std::max<Coord>(1, isqrt_ceil(D2));
    auto fdiv = [](Coord a, Coord b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); };
    std::vector<std::pair<Point, int>> cell(n);
    for (int i = 0; i < n; ++i) cell[i] = {Point(fdiv(ps[i].x, s), fdiv(ps[i].y, s)), i};
    std::sort(cell.begin(), cell.end());
    std::unordered_map<Point, std::pair<int, int>, PointHash> range;  // マス → cell の中の範囲 [l, r)
    range.reserve((size_t)n * 2 + 1);
    for (int l = 0, r; l < n; l = r) {
        for (r = l; r < n && cell[r].first == cell[l].first; ++r) {}
        range[cell[l].first] = {l, r};
    }
    auto add = [&](int a, int b) {
        if (dist2(ps[a], ps[b]) <= D2) res.push_back({std::min(a, b), std::max(a, b)});
    };
    const Point nb[4] = {Point(1, -1), Point(1, 0), Point(1, 1), Point(0, 1)};  // 隣のマスは片側の 4 つだけ見る（同じ組を 2 回見ない）
    for (const auto& [c, lr] : range) {
        for (int a = lr.first; a < lr.second; ++a)
            for (int b = a + 1; b < lr.second; ++b) add(cell[a].second, cell[b].second);
        for (const Point& d : nb) {
            const auto it = range.find(c + d);
            if (it == range.end()) continue;
            for (int a = lr.first; a < lr.second; ++a)
                for (int b = it->second.first; b < it->second.second; ++b) add(cell[a].second, cell[b].second);
        }
    }
    std::sort(res.begin(), res.end());
    return res;
}
inline std::vector<std::pair<int, int>> pairs_within(const std::vector<Point>& ps, Coord D) {  // 距離 ≤ D（D は 3×10^9 まで）
    return D < 0 ? std::vector<std::pair<int, int>>() : pairs_within_d2(ps, D * D);
}

// ---- 一直線・三角形・直線 ----
inline long long count_collinear_triples(const std::vector<Point>& ps) {  // 一直線に並ぶ 3 点の組（番号の組）の数
    const int n = (int)ps.size();
    long long tot = 0;
    std::vector<std::pair<Point, int>> buf;
    for (int i = 0; i < n; ++i) {
        int z;
        long long acc = 0;
        geo_internal::line_groups(ps, i, buf, z, [&](int g, int) { acc += geo_internal::choose2(g); });
        const long long m = n - 1 - i;
        tot += acc + geo_internal::choose2(z) + (long long)z * (m - z);
    }
    return tot;
}
inline long long count_triangles(const std::vector<Point>& ps) {  // 面積が正の三角形の数
    const long long n = (long long)ps.size();
    return n * (n - 1) * (n - 2) / 6 - count_collinear_triples(ps);
}
inline MaxCollinear max_collinear(const std::vector<Point>& ps) {  // 1 本の直線の上に並ぶ最大の点の数
    const int n = (int)ps.size();
    MaxCollinear best;
    std::vector<std::pair<Point, int>> buf;
    for (int i = 0; i < n; ++i) {
        int z;
        bool any = false;
        geo_internal::line_groups(ps, i, buf, z, [&](int g, int j) {
            any = true;
            if (g + z + 1 > best.count) best = MaxCollinear{g + z + 1, i, j};
        });
        if (!any && z + 1 > best.count) best = MaxCollinear{z + 1, i, -1};
    }
    return best;
}
// 2 つ以上の異なる点を通る直線ごとに、その上の点の数（同じ点は重ねて数える）。LineKey の順
inline std::vector<LineCount> lines_through(const std::vector<Point>& ps) {
    const int n = (int)ps.size();
    std::vector<std::pair<LineKey, int>> all;
    std::vector<std::pair<Point, int>> buf;
    for (int i = 0; i < n; ++i) {
        int z;
        geo_internal::line_groups(ps, i, buf, z, [&](int g, int j) { all.push_back({line_key(ps[i], ps[j]), g + z + 1}); });
    }
    std::sort(all.begin(), all.end(), [](const std::pair<LineKey, int>& a, const std::pair<LineKey, int>& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second > b.second;
    });
    std::vector<LineCount> res;
    for (size_t k = 0; k < all.size(); ++k)
        if (k == 0 || all[k].first != all[k - 1].first) res.push_back(LineCount{all[k].first, all[k].second});
    return res;
}
// K 点以上を通る直線の数。K ≤ 1、または同じ点が K 個以上あるとき（その点を通る直線が全部当てはまる）は無限なので -1
inline long long count_lines_with_at_least(const std::vector<Point>& ps, int K) {
    if (K <= 1) return -1;
    const std::vector<Point> s = [&] { std::vector<Point> t = ps; std::sort(t.begin(), t.end()); return t; }();
    for (size_t k = 0, l; k < s.size(); k = l) {
        for (l = k; l < s.size() && s[l] == s[k]; ++l) {}
        if ((int)(l - k) >= K) return -1;
    }
    const int n = (int)ps.size();
    std::vector<LineKey> keys;  // K 点以上の直線（同じ直線は何度も出てくる）
    std::vector<std::pair<Point, int>> buf;
    for (int i = 0; i < n; ++i) {
        int z;
        geo_internal::line_groups(ps, i, buf, z, [&](int g, int j) {
            if (g + z + 1 >= K) keys.push_back(line_key(ps[i], ps[j]));
        });
    }
    std::sort(keys.begin(), keys.end());
    return (long long)(std::unique(keys.begin(), keys.end()) - keys.begin());
}

// ---- 角度 ----
namespace geo_internal {
// ps[i] から見た、ps[i] と違う点への向き（偏角の順）と番号
inline void sorted_dirs(const std::vector<Point>& ps, int i, std::vector<std::pair<Point, int>>& d) {
    d.clear();
    for (int j = 0; j < (int)ps.size(); ++j)
        if (ps[j] != ps[i]) d.push_back({ps[j] - ps[i], j});
    // 上半分 [0, π) と下半分 [π, 2π) に分けてから、それぞれを外積で並べる
    const auto mid = std::partition(d.begin(), d.end(), [](const std::pair<Point, int>& a) { return arg_half(a.first) == 0; });
    auto by_cross = [](const std::pair<Point, int>& a, const std::pair<Point, int>& b) { return cross(a.first, b.first) > 0; };
    std::sort(d.begin(), mid, by_cross);
    std::sort(mid, d.end(), by_cross);
}
// w の、v から反時計回りに測った角の区分：0: [0, π/2) / 1: π/2 / 2: (π/2, π) / 3: [π, 2π)
inline int quarter(const Point& v, const Point& w) {
    const Coord c = cross(v, w), t = dot(v, w);
    if (c < 0 || (c == 0 && t < 0)) return 3;
    if (t > 0) return 0;
    return t == 0 ? 1 : 2;
}
}  // namespace geo_internal
// 鋭角・直角・鈍角三角形の数（面積が正の三角形だけ。番号の組で数える）
inline TriangleCounts count_triangles_by_angle(const std::vector<Point>& ps) {
    const int n = (int)ps.size();
    TriangleCounts r;
    std::vector<std::pair<Point, int>> d;
    for (int i = 0; i < n; ++i) {  // 頂点 i の角が直角・鈍角になる組を数える（1 つの三角形で直角・鈍角の頂点は 1 つだけ）
        geo_internal::sorted_dirs(ps, i, d);
        const int m = (int)d.size();
        int a = 0, b = 0, c = 0, g = 0;  // k から反時計回りに並べたとき、区分が 1 以上 / 2 以上 / 3 の最初の位置。g は k と同じ向きの最初
        for (int k = 0; k < m; ++k) {
            const Point& v = d[k].first;
            if (k > 0 && !same_direction(d[k - 1].first, v)) g = k;
            const int e = g + m;  // 見る範囲は [k, e)（k より前の同じ向きの点を除く）
            a = std::max(a, k), b = std::max(b, a), c = std::max(c, b);
            while (a < e && geo_internal::quarter(v, d[a % m].first) < 1) ++a;
            b = std::max(b, a);
            while (b < e && geo_internal::quarter(v, d[b % m].first) < 2) ++b;
            c = std::max(c, b);
            while (c < e && geo_internal::quarter(v, d[c % m].first) < 3) ++c;
            r.right += b - a, r.obtuse += c - b;
        }
    }
    r.acute = count_triangles(ps) - r.right - r.obtuse;
    return r;
}
// 最大の角 ∠p_i p_j p_k（j が頂点。i, k は j と違う点）。典型90 009
inline MaxAngle max_angle(const std::vector<Point>& ps) {
    const int n = (int)ps.size();
    MaxAngle best;
    Point ba, bb;  // いちばん大きい角の 2 本のベクトル
    std::vector<std::pair<Point, int>> d;
    for (int j = 0; j < n; ++j) {
        geo_internal::sorted_dirs(ps, j, d);
        const int m = (int)d.size();
        if (m < 2) continue;
        int c = 0, g = 0;  // k から反時計回りに並べたとき、π 以上の最初の位置（反対の向きにいちばん近い点の候補は c と c - 1）
        for (int k = 0; k < m; ++k) {
            if (k > 0 && !same_direction(d[k - 1].first, d[k].first)) g = k;
            const int e = g + m;  // 見る範囲は [k, e)
            c = std::max(c, k + 1);
            while (c < e && geo_internal::quarter(d[k].first, d[c % m].first) < 3) ++c;
            for (int t : {c, c - 1}) {
                if (t <= k || t >= e) continue;
                const Point& w = d[t % m].first;
                if (best.i < 0 || angle_cmp(d[k].first, w, ba, bb) > 0) best = MaxAngle{0, d[k].second, j, d[t % m].second}, ba = d[k].first, bb = w;
            }
        }
    }
    if (best.i >= 0) best.angle = angle_between(ba, bb);
    return best;
}

// ---- 図形の数（同じ点は 1 つにまとめる） ----
inline long long count_axis_rectangles(const std::vector<Point>& ps) {  // 辺が軸に平行な長方形（面積 > 0）
    const std::vector<Point> s = geo_internal::unique_points(ps);
    std::vector<std::pair<Coord, Coord>> keys;  // 同じ x の列の 2 点の (y1, y2)
    for (size_t a = 0, b; a < s.size(); a = b) {
        for (b = a; b < s.size() && s[b].x == s[a].x; ++b) {}
        for (size_t p = a; p < b; ++p)
            for (size_t q = p + 1; q < b; ++q) keys.push_back({s[p].y, s[q].y});
    }
    std::sort(keys.begin(), keys.end());
    long long res = 0;
    for (size_t a = 0, b; a < keys.size(); a = b) {
        for (b = a; b < keys.size() && keys[b] == keys[a]; ++b) {}
        res += geo_internal::choose2((long long)(b - a));
    }
    return res;
}
inline long long count_squares(const std::vector<Point>& ps) {  // 正方形（傾いてもよい、面積 > 0）
    const std::vector<Point> s = geo_internal::unique_points(ps);
    const std::unordered_set<Point, PointHash> st(s.begin(), s.end());
    long long c = 0;
    for (const Point& a : s)
        for (const Point& b : s) {
            if (a == b) continue;
            const Point r = rot90(b - a);  // a → b の左側に正方形
            if (st.count(b + r) && st.count(a + r)) ++c;
        }
    return c / 4;
}
inline long long count_parallelograms(const std::vector<Point>& ps) {  // 平行四辺形（4 点が一直線のつぶれたものは除く）
    const std::vector<Point> s = geo_internal::unique_points(ps);
    const int m = (int)s.size();
    std::vector<std::pair<Point, Point>> es;  // (2 点の和 = 中点の 2 倍, 向き)
    es.reserve((size_t)m * (m - 1) / 2);
    for (int i = 0; i < m; ++i)
        for (int j = i + 1; j < m; ++j) es.push_back({s[i] + s[j], geo_internal::upper_dir(s[j] - s[i])});
    std::sort(es.begin(), es.end(), [](const std::pair<Point, Point>& a, const std::pair<Point, Point>& b) {
        if (a.first != b.first) return a.first < b.first;
        return cross(a.second, b.second) > 0;
    });
    long long res = 0;
    for (size_t a = 0, b; a < es.size(); a = b) {
        for (b = a; b < es.size() && es[b].first == es[a].first; ++b) {}
        res += geo_internal::choose2((long long)(b - a));  // 対角線の中点が同じ 2 組
        for (size_t c = a, d; c < b; c = d) {           // 同じ直線の上の 2 組はつぶれている
            for (d = c; d < b && cross(es[c].second, es[d].second) == 0; ++d) {}
            res -= geo_internal::choose2((long long)(d - c));
        }
    }
    return res;
}
// 台形（1 組以上の対辺が平行な四角形。平行四辺形も含む）。ちょうど 1 組だけなら count_trapezoids - count_parallelograms
inline long long count_trapezoids(const std::vector<Point>& ps) {
    const std::vector<Point> s = geo_internal::unique_points(ps);
    const int m = (int)s.size();
    std::vector<std::pair<Point, int>> es;  // (向き, 線分の端点の番号)
    es.reserve((size_t)m * (m - 1) / 2);
    for (int i = 0; i < m; ++i)
        for (int j = i + 1; j < m; ++j) es.push_back({geo_internal::upper_dir(s[j] - s[i]), i});
    std::sort(es.begin(), es.end(), [](const std::pair<Point, int>& a, const std::pair<Point, int>& b) { return cross(a.first, b.first) > 0; });
    long long par = 0;  // 平行で別の直線の上にある線分の組の数
    std::vector<Coord> key;
    for (size_t a = 0, b; a < es.size(); a = b) {
        for (b = a; b < es.size() && cross(es[a].first, es[b].first) == 0; ++b) {}
        par += geo_internal::choose2((long long)(b - a));
        key.clear();
        for (size_t c = a; c < b; ++c) key.push_back(cross(es[a].first, s[es[c].second]));  // 平行な直線を区別する値
        std::sort(key.begin(), key.end());
        for (size_t c = 0, d; c < key.size(); c = d) {
            for (d = c; d < key.size() && key[d] == key[c]; ++d) {}
            par -= geo_internal::choose2((long long)(d - c));
        }
    }
    return par - count_parallelograms(s);  // 平行四辺形は 2 回数えている
}
#endif  // GG_POINT_SET
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
namespace ps_test {

#define PS_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[point_set] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK PS_CHECK

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

static std::vector<Point> rnd_set(std::mt19937_64& rng, int n, long long C, bool dup) {
    std::vector<Point> ps = rndps(rng, n, C);
    if (dup && n > 1)
        for (int k = 0; k < n / 3; ++k) ps[rng() % n] = ps[rng() % n];
    return ps;
}
static bool col3(const Point& a, const Point& b, const Point& c) { return cross128(b - a, c - a) == 0; }

static void test_closest(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, it % 10 == 0 ? 300 : 30);
        const long long C = it % 3 == 0 ? 1000000000LL : (it % 3 == 1 ? 5 : 1000);
        const std::vector<Point> ps = rnd_set(rng, n, C, it % 4 == 0);
        const ClosestPair r = closest_pair(ps);
        if (n < 2) {
            CHECK(r.i == -1 && r.d2 == -1);
        } else {
            Coord best = std::numeric_limits<Coord>::max();
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j) best = std::min(best, dist2(ps[i], ps[j]));
            CHECK(r.d2 == best && r.i < r.j && dist2(ps[r.i], ps[r.j]) == best);
        }
        // 距離 D 以内
        const Coord D = it % 5 == 0 ? 0 : rnd(rng, 0, std::min<Coord>(C * 2, 3000000000LL));
        std::vector<std::pair<int, int>> e;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if ((Int128)dist2(ps[i], ps[j]) <= (Int128)D * D) e.push_back({i, j});
        CHECK(pairs_within(ps, D) == e);
        const Coord D2 = rnd(rng, 0, 2 * C * C);
        e.clear();
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (dist2(ps[i], ps[j]) <= D2) e.push_back({i, j});
        CHECK(pairs_within_d2(ps, D2) == e);
    }
}

static void test_lines(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 14);
        const bool big = it % 4 == 0;
        std::vector<Point> ps = rnd_set(rng, n, big ? 1000000000LL : 2, it % 3 == 0);
        if (big && n > 2) {  // 大きい座標で一直線の点を作る
            const Point d = rndp(rng, 1000), o = rndp(rng, 500000000LL);
            for (int k = 0; k < n / 2; ++k) ps[rng() % n] = o + d * rnd(rng, -100000, 100000);
        }
        long long col = 0;
        for (int a = 0; a < n; ++a)
            for (int b = a + 1; b < n; ++b)
                for (int c = b + 1; c < n; ++c) col += col3(ps[a], ps[b], ps[c]);
        CHECK(count_collinear_triples(ps) == col && count_triangles(ps) == (long long)n * (n - 1) * (n - 2) / 6 - col);
        // 直線ごとの点の数
        std::map<LineKey, int> lines;
        int mx = n > 0 ? 1 : 0;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b) {
                if (ps[a] == ps[b]) continue;
                int c = 0;
                for (int k = 0; k < n; ++k) c += col3(ps[a], ps[b], ps[k]);
                lines[line_key(ps[a], ps[b])] = c;
                mx = std::max(mx, c);
            }
        int mult = 0;
        for (int a = 0; a < n; ++a) {
            int c = 0;
            for (int b = 0; b < n; ++b) c += ps[a] == ps[b];
            mult = std::max(mult, c);
        }
        if (lines.empty()) mx = std::max(mx, mult);
        const MaxCollinear m = max_collinear(ps);
        CHECK(m.count == mx);
        if (m.j >= 0) {
            int c = 0;
            for (int k = 0; k < n; ++k) c += col3(ps[m.i], ps[m.j], ps[k]);
            CHECK(ps[m.i] != ps[m.j] && c == mx);
        }
        const std::vector<LineCount> lt = lines_through(ps);
        CHECK(lt.size() == lines.size());
        size_t k = 0;
        for (const auto& [key, c] : lines) {
            CHECK(lt[k].key == key && lt[k].count == c);
            ++k;
        }
        for (int K = 0; K <= 5; ++K) {
            long long e = 0;
            for (const auto& kv : lines) e += kv.second >= K;
            if (K <= 1 || mult >= K) e = -1;
            CHECK(count_lines_with_at_least(ps, K) == e);
        }
    }
}

static void test_angles(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 12);
        const bool big = it % 4 == 0;
        std::vector<Point> ps = rnd_set(rng, n, big ? 1000000000LL : 3, it % 3 == 0);
        if (big && n >= 3) {  // 直角三角形を混ぜる
            const Point a = rndp(rng, 400000000LL), v = rndp(rng, 400000000LL);
            ps[0] = a, ps[1] = a + v, ps[2] = a + rot90(v);
        }
        TriangleCounts e;
        for (int a = 0; a < n; ++a)
            for (int b = a + 1; b < n; ++b)
                for (int c = b + 1; c < n; ++c) {
                    if (col3(ps[a], ps[b], ps[c])) continue;
                    Int128 s[3] = {norm2_128(ps[b] - ps[c]), norm2_128(ps[c] - ps[a]), norm2_128(ps[a] - ps[b])};
                    std::sort(s, s + 3);
                    if (s[2] == s[0] + s[1]) ++e.right;
                    else if (s[2] > s[0] + s[1]) ++e.obtuse;
                    else ++e.acute;
                }
        const TriangleCounts r = count_triangles_by_angle(ps);
        CHECK(r.acute == e.acute && r.right == e.right && r.obtuse == e.obtuse);
        // 最大の角
        Real best = -1;
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < n; ++k)
                    if (ps[i] != ps[j] && ps[k] != ps[j] && i != k) best = std::max(best, angle_between(ps[i] - ps[j], ps[k] - ps[j]));
        const MaxAngle ma = max_angle(ps);
        CHECK(near(ma.angle, best, TOL * 16));
        if (best >= 0) CHECK(ps[ma.i] != ps[ma.j] && ps[ma.k] != ps[ma.j] && ma.i != ma.k && near(angle_between(ps[ma.i] - ps[ma.j], ps[ma.k] - ps[ma.j]), best, TOL * 16));
        else CHECK(ma.i == -1);
    }
}

static void test_shapes(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 11);
        const long long C = it % 2 ? 2 : 3;
        std::vector<Point> ps = rnd_set(rng, n, C, it % 3 == 0);
        const std::vector<Point> s = geo_internal::unique_points(ps);
        const int m = (int)s.size();
        long long rect = 0, sq = 0, par = 0, trap = 0;
        for (int a = 0; a < m; ++a)
            for (int b = a + 1; b < m; ++b)
                for (int c = b + 1; c < m; ++c)
                    for (int d = c + 1; d < m; ++d) {
                        const Point q[4] = {s[a], s[b], s[c], s[d]};
                        // 軸に平行な長方形：x が 2 種類・y が 2 種類で、4 つの組み合わせが全部ある
                        std::set<Coord> xs, ys;
                        for (const Point& p : q) xs.insert(p.x), ys.insert(p.y);
                        if (xs.size() == 2 && ys.size() == 2) ++rect;
                        // 正方形：6 本の距離のうち小さい 4 本が等しく、大きい 2 本がその 2 倍
                        std::vector<Coord> dd;
                        for (int i = 0; i < 4; ++i)
                            for (int j = i + 1; j < 4; ++j) dd.push_back(dist2(q[i], q[j]));
                        std::sort(dd.begin(), dd.end());
                        if (dd[0] > 0 && dd[0] == dd[3] && dd[4] == dd[5] && dd[4] == 2 * dd[0]) ++sq;
                        // 3 通りの 2 組への分け方
                        const int pr[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
                        const bool allcol = col3(q[0], q[1], q[2]) && col3(q[0], q[1], q[3]);
                        bool isp = false, ist = false;
                        for (const auto& P : pr) {
                            const Point A = q[P[0]], B = q[P[1]], Cc = q[P[2]], D = q[P[3]];
                            if (!allcol && A + B == Cc + D) isp = true;                          // 対角線の中点が同じ
                            if (cross128(B - A, D - Cc) == 0 && !col3(A, B, Cc)) ist = true;   // 平行で別の直線
                        }
                        par += isp, trap += ist;
                    }
        CHECK(count_axis_rectangles(ps) == rect && count_squares(ps) == sq);
        CHECK(count_parallelograms(ps) == par && count_trapezoids(ps) == trap);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_closest(rng, cases);
    test_lines(rng, cases);
    test_angles(rng, cases);
    test_shapes(rng, cases);
    std::cout << "[point_set] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[point_set] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 200000;
    const std::vector<Point> ps = rndps(rng, n, 1000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    const ClosestPair cp = closest_pair(ps);
    const long long t1 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const auto pw = pairs_within(ps, 5000000);
    const long long t2 = ms_since(t0);
    const int m = 2000;
    const std::vector<Point> qs = rndps(rng, m, 30);  // 一直線・同じ向きが多い
    t0 = std::chrono::steady_clock::now();
    const long long col = count_collinear_triples(qs);
    const long long t3 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const long long lk = count_lines_with_at_least(rndps(rng, m, 1000000000LL), 2);
    const long long t4 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const TriangleCounts tc = count_triangles_by_angle(qs);
    const long long t5 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const MaxAngle ma = max_angle(qs);
    const long long t6 = ms_since(t0);
    const std::vector<Point> gs = rndps(rng, m, 40);
    t0 = std::chrono::steady_clock::now();
    const long long sh = count_axis_rectangles(gs) + count_squares(gs) + count_parallelograms(gs) + count_trapezoids(gs);
    const long long t7 = ms_since(t0);
    std::cout << "[point_set] benchmark: N = " << n << "  closest_pair " << t1 << " ms, pairs_within " << t2 << " ms\n"
              << "[point_set] benchmark: N = " << m << "  count_collinear_triples " << t3 << " ms, count_lines_with_at_least " << t4
              << " ms, count_triangles_by_angle " << t5 << " ms, max_angle " << t6 << " ms, 長方形+正方形+平行四辺形+台形 " << t7 << " ms  (check "
              << cp.d2 % 7 + (long long)pw.size() % 7 + col % 7 + lk % 7 + tc.obtuse % 7 + ma.i + sh % 7 << ")\n";
}

}  // namespace ps_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace ps_demo {

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
#define SHOW(expr, comment) ps_demo::show(#expr, ps_demo::to_s(expr), comment)
#define RUN(...)                   \
    ps_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("点の集合 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  closest_pair(ps) / pairs_within(ps, D) / pairs_within_d2(ps, D2)          最も近い 2 点 / 距離 D 以内の組\n"
              << "  count_collinear_triples(ps) / count_triangles(ps) / max_collinear(ps)    一直線の 3 点 / 三角形の数 / 一直線の最大の点の数\n"
              << "  lines_through(ps) / count_lines_with_at_least(ps, K)                       直線ごとの点の数 / K 点以上を通る直線の数\n"
              << "  count_triangles_by_angle(ps) / max_angle(ps)                             鋭角・直角・鈍角の数 / 最大の角\n"
              << "  count_axis_rectangles / count_squares / count_parallelograms / count_trapezoids   図形の数\n";

    api("【1】 closest_pair(ps) / pairs_within(ps, D)", "ClosestPair{i, j, d2} / vector<pair<int, int>>",
        "最も近い 2 点（番号 i < j と距離の 2 乗）/ 距離 D 以内の組を全部（pairs_within_d2 なら距離² ≤ D2）", "O(N log N) / O(N log N + 答えの数)");
    RUN(std::vector<Point> ps = {{0, 0}, {5, 5}, {1, 2}, {9, 9}, {6, 4}});
    RUN(ClosestPair c = closest_pair(ps));
    std::cout << "  c.i, c.j, c.d2                     -> " << c.i << ", " << c.j << ", " << c.d2 << "         // (5,5) と (6,4)：距離² = 2\n";
    SHOW(pairs_within(ps, 3).size(), "距離 3 以内の組の数");
    note("ABC304 C（距離 D 以内の人にうつる）: for (auto [i, j] : pairs_within(ps, D)) uf.unite(i, j);  → 0 とつながっているか");

    api("【2】 count_collinear_triples(ps) / count_triangles(ps) / max_collinear(ps)", "long long / MaxCollinear{count, i, j}",
        "一直線に並ぶ 3 点の組の数 / 面積が正の三角形の数 / 1 本の直線に並ぶ最大の点の数", "O(N^2 log N)");
    RUN(std::vector<Point> qs = {{0, 0}, {1, 1}, {2, 2}, {3, 0}, {0, 3}, {1, 2}});
    SHOW(count_collinear_triples(qs), "y = x の 3 点、x + y = 3 の 3 点 (3,0) (1,2) (0,3)");
    SHOW(count_triangles(qs), "C(6, 3) - 2");
    SHOW(max_collinear(qs).count, "1 本の直線の上に 3 点まで");
    note("ABC181 C（3 点が一直線か）: count_collinear_triples(ps) > 0。ABC224 C（三角形の数）: count_triangles(ps)");

    api("【3】 lines_through(ps) / count_lines_with_at_least(ps, K)", "vector<LineCount{key, count}> / long long",
        "2 点以上を通る直線ごとの点の数 / K 点以上を通る直線の数（K ≤ 1 など無限にあるときは -1）", "O(N^2 log N)");
    SHOW(lines_through(qs).size(), "異なる直線の数");
    SHOW(count_lines_with_at_least(qs, 3), "3 点以上を通る直線");
    SHOW(count_lines_with_at_least(qs, 1), "無限 → -1");
    note("ABC248 E（K 点以上を通る直線）: auto r = count_lines_with_at_least(ps, K); r == -1 なら Infinity");

    api("【4】 count_triangles_by_angle(ps) / max_angle(ps)", "TriangleCounts{acute, right, obtuse} / MaxAngle{angle, i, j, k}",
        "鋭角・直角・鈍角三角形の数（ABC033 D）/ 最大の角 ∠p_i p_j p_k（j が頂点、ラジアン。典型90 009）", "O(N^2 log N)");
    RUN(TriangleCounts t = count_triangles_by_angle(std::vector<Point>{{0, 0}, {4, 0}, {0, 3}, {4, 3}, {1, 1}}));
    std::cout << "  t.acute, t.right, t.obtuse         -> " << t.acute << ", " << t.right << ", " << t.obtuse << "\n";
    RUN(MaxAngle a = max_angle(std::vector<Point>{{0, 0}, {10, 0}, {5, 1}}));
    SHOW(rad_to_deg(a.angle), "頂点 (5, 1) の角（度）");
    std::cout << "  a.i, a.j, a.k                      -> " << a.i << ", " << a.j << ", " << a.k << "          // j が頂点\n";

    api("【5】 count_axis_rectangles / count_squares / count_parallelograms / count_trapezoids", "long long",
        "4 頂点が全部ある長方形（軸に平行）/ 正方形（傾いてもよい）/ 平行四辺形 / 台形（平行四辺形を含む）の数", "O(N^2) 〜 O(N^2 log N)");
    RUN(std::vector<Point> g;  for (Coord x = 0; x < 3; ++x) for (Coord y = 0; y < 3; ++y) g.push_back(Point(x, y)));
    SHOW(count_axis_rectangles(g), "3×3 の格子点：軸に平行な長方形");
    SHOW(count_squares(g), "正方形（傾いた 1 個を含む）");
    SHOW(count_parallelograms(g), "平行四辺形");
    SHOW(count_trapezoids(g), "台形（平行四辺形を含む）");
    note("ABC218 D（長方形の数）: count_axis_rectangles(ps)。ABC275 C（9×9 の盤の正方形）: '#' のマスを点にして count_squares");
    note("ちょうど 1 組だけ平行な台形 = count_trapezoids(ps) - count_parallelograms(ps)");
}

}  // namespace ps_demo


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
    ps_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ps_test::benchmark();    // 速度計測
    ps_demo::demo();         // 使い方デモ
    return 0;
}

