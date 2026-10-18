// =============================================================================
//  直線・線分（交差判定・交点・距離・射影・直線の正規化・格子点）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 02 直線・線分
//  依存       : 00_geo_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - Line l(a, b) は直線 a → b（a ≠ b。左側 = ccw(a, b, p) > 0）、Segment s(a, b) は線分（a == b の点も可）
//    - 整数の型（Line / Segment / Point）なら判定は厳密。実数の型（LineD ...）なら相対誤差 GEO_EPS
//    - 交点・射影・距離は実数で返す（交点は vector<PointD>。なしは空）。整数の入力なら分子を 128 ビットで正確に求めてから割る
//    - 直線の関数は a と向き b - a しか使わないので、perpendicular などで作った直線の b は座標の範囲を少し超えてもよい
//
//  機能一覧                                                                    計算量
//    is_parallel(l1, l2) / is_orthogonal(l1, l2)      平行 / 垂直（直線・線分どちらでも）      O(1)
//    on_line(l, p) / on_segment(s, p) / side(l, p)    直線の上 / 線分の上（端点を含む）/ 左 +1・右 -1・上 0
//    intersects(s1, s2)                               線分どうしが共有点を持つか（端点で触れる・重なるも含む）
//    intersects_strict(s1, s2)                        真に交差（1 点で交わり、その点がどちらの端点でもない）
//    intersects(l, s) / intersects(s, l) / intersects(l1, l2)   直線と線分 / 直線どうし
//    line_relation(l1, l2)                            LINE_CROSS / LINE_PARALLEL / LINE_SAME（同じ直線）
//    segment_relation(s1, s2)                         SEGMENT_SEPARATE（共有点なし）/ SEGMENT_TOUCH（1 点で接する：端点・T 字）
//                                                     / SEGMENT_CROSS（真に交差）/ SEGMENT_OVERLAP（一直線で長さのある部分が重なる）
//    projection(l, p) / reflection(l, p)              垂線の足 / 線対称な点（PointD）
//    dist(l, p) / dist(s, p) / dist(s1, s2) / dist(l, s) / dist(l1, l2)   距離（Real。引数の順は逆でもよい）
//    dist2_cmp(l or s, p, D2)                         距離の 2 乗と整数 D2 の比較（-1 / 0 / 1。誤差なし）
//    cross_point(l1, l2)                              直線の交点（平行・同じ直線なら空）
//    cross_point(l, s) / cross_point(s, l)            直線と線分の交点（線分が直線の上なら両端）
//    cross_points(s1, s2)                             線分の交点（重なるときは重なりの両端。1 点に重なれば 1 個）
//    cross_point_frac(l1, l2)                         交点を既約分数 FracPoint{x, y, d}（= (x/d, y/d)）で。同じ交点の判定に
//    perpendicular(l, p) / parallel_through(l, p)     p を通る垂線 / p を通る平行線（整数のまま）
//    bisector(a, b)                                   垂直二等分線（LineD）
//    LineKey line_key(p, q) / line_key(l)             直線 a x + b y = c を正規化したキー（同じ直線なら同じ。==・<・LineKeyHash）
//      k.contains(p) / k.dir_key() / k.line()         点が直線の上か / 傾きのキー / 実数の直線に
//    bisector_key(p, q)                               垂直二等分線の LineKey
//    lattice_points_on(s) / cells_crossed(s)          線分の上の格子点の数（gcd + 1）/ 線分が内部を通るマスの数
//    Ray r(a, b) / RayD                               半直線（a から b の向きへ。a を含む）。is_parallel / is_orthogonal も使える
//      on_ray(r, p) / intersects(r, 直線・線分・半直線)   半直線の上か / 共有点を持つか（引数の順は逆でもよい。整数なら厳密）
//      cross_point(r, l) / cross_point(r, r2)         交点がちょうど 1 つのときその点（なし・一直線で重なるなら空）
//      cross_points(r, s)                             半直線と線分の交点（重なるなら重なりの両端）
//      dist(r, 点・直線・線分・半直線) / dist2_cmp(r, p, D2)   距離（Real）/ 距離² と整数の比較
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    ln_test::random_test()  … 愚直な実装（端点の判定・有理数・サンプリング）と突き合わせるランダムテスト
//    ln_test::benchmark()    … 速度計測
//    ln_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
namespace ln_test {

#define LN_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[line] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK LN_CHECK

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

static bool b_on_seg(const Point& a, const Point& b, const Point& p) {
    return cross128(b - a, p - a) == 0 && std::min(a.x, b.x) <= p.x && p.x <= std::max(a.x, b.x) && std::min(a.y, b.y) <= p.y &&
           p.y <= std::max(a.y, b.y);
}
static bool b_strict(const Segment& s, const Segment& t) {
    return sgn128(cross128(s.b - s.a, t.a - s.a)) * sgn128(cross128(s.b - s.a, t.b - s.a)) < 0 &&
           sgn128(cross128(t.b - t.a, s.a - t.a)) * sgn128(cross128(t.b - t.a, s.b - t.a)) < 0;
}
static bool b_inter(const Segment& s, const Segment& t) {
    return b_on_seg(s.a, s.b, t.a) || b_on_seg(s.a, s.b, t.b) || b_on_seg(t.a, t.b, s.a) || b_on_seg(t.a, t.b, s.b) || b_strict(s, t);
}
// 点と線分の距離の 2 乗（有理数 num / den）
static void b_dist2(const Segment& s, const Point& p, Int128& num, Int128& den) {
    const Point d = s.b - s.a;
    const Int128 t = dot128(d, p - s.a), n = norm2_128(d);
    if (n == 0 || t <= 0) { num = norm2_128(p - s.a), den = 1; return; }
    if (t >= n) { num = norm2_128(p - s.b), den = 1; return; }
    const Int128 c = cross128(d, p - s.a);
    num = c * c, den = n;
}
static Segment rnd_seg(std::mt19937_64& rng, long long C) {
    Segment s(rndp(rng, C), rndp(rng, C));
    if (rng() % 8 == 0) s.b = s.a;
    return s;
}

static void test_predicates(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 2 : 5;
        const Segment s = rnd_seg(rng, C), t = rnd_seg(rng, C);
        const Point p = rndp(rng, C);
        CHECK(on_segment(s, p) == b_on_seg(s.a, s.b, p));
        CHECK(intersects(s, t) == b_inter(s, t) && intersects(t, s) == b_inter(s, t));
        CHECK(intersects_strict(s, t) == b_strict(s, t));
        CHECK(is_parallel(s, t) == (cross128(s.b - s.a, t.b - t.a) == 0) && is_orthogonal(s, t) == (dot128(s.b - s.a, t.b - t.a) == 0));
        if (s.a == s.b || t.a == t.b) continue;
        const Line l(s.a, s.b), m(t.a, t.b);
        CHECK(on_line(l, p) == (cross128(l.b - l.a, p - l.a) == 0) && side(l, p) == sgn128(cross128(l.b - l.a, p - l.a)));
        const bool lt = sgn128(cross128(l.b - l.a, t.a - l.a)) * sgn128(cross128(l.b - l.a, t.b - l.a)) <= 0;
        CHECK(intersects(l, t) == lt && intersects(t, l) == lt);
        const int rel = cross128(l.b - l.a, m.b - m.a) != 0 ? LINE_CROSS : (cross128(l.b - l.a, m.a - l.a) == 0 ? LINE_SAME : LINE_PARALLEL);
        CHECK(line_relation(l, m) == rel && intersects(l, m) == (rel != LINE_PARALLEL));
        CHECK(is_parallel(l, t) == is_parallel(s, m));
    }
    // 大きい座標：端点で触れる・一直線で重なるケースを a + k d で作る
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const Point d = rndp(rng, 50);
        if (d == Point()) continue;
        const Point a = rndp(rng, 900000000LL);
        const long long K = 1000000;
        Segment s(a + d * rnd(rng, -K, K), a + d * rnd(rng, -K, K)), t(a + d * rnd(rng, -K, K), a + d * rnd(rng, -K, K));
        if (it % 3 == 0) t.b = rndp(rng, 1000000000LL);
        if (it % 5 == 0) t.a = s.b;
        CHECK(intersects(s, t) == b_inter(s, t) && intersects_strict(s, t) == b_strict(s, t));
        const Segment u(rndp(rng, 1000000000LL), rndp(rng, 1000000000LL)), v(rndp(rng, 1000000000LL), rndp(rng, 1000000000LL));
        CHECK(intersects(u, v) == b_inter(u, v) && intersects_strict(u, v) == b_strict(u, v));
        CHECK(on_segment(s, t.a) == b_on_seg(s.a, s.b, t.a));
    }
}

static void test_metric(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 6 : 1000000000LL;
        const Segment s = rnd_seg(rng, C), t = rnd_seg(rng, C);
        const Point p = rndp(rng, C);
        // 点と線分の距離
        Int128 num, den;
        b_dist2(s, p, num, den);
        CHECK(near(dist(s, p), std::sqrt((Real)num / (Real)den), TOL) && dist(p, s) == dist(s, p));
        for (int q = 0; q < 3; ++q) {
            Coord D2 = (Coord)std::llround((Real)num / (Real)den) + rnd(rng, -1, 1);
            if (D2 < 0) D2 = 0;
            CHECK(dist2_cmp(s, p, D2) == sgn128(num - (Int128)D2 * den));
        }
        // 線分どうし
        Real e = 1e30L;
        for (const auto& [sg, q] : {std::make_pair(s, t.a), std::make_pair(s, t.b), std::make_pair(t, s.a), std::make_pair(t, s.b)}) {
            Int128 a2, b2;
            b_dist2(sg, q, a2, b2);
            e = std::min(e, std::sqrt((Real)a2 / (Real)b2));
        }
        if (b_inter(s, t)) e = 0;
        CHECK(near(dist(s, t), e, TOL));
        if (s.a == s.b) continue;
        const Line l(s.a, s.b);
        // 直線と点
        const Int128 c = cross128(l.b - l.a, p - l.a);
        CHECK(near(dist(l, p), std::fabs((Real)c) / len(l.b - l.a), TOL) && dist(p, l) == dist(l, p));
        for (int q = 0; q < 3; ++q) {
            const Int128 n2 = norm2_128(l.b - l.a);
            Coord D2 = (Coord)std::llround((Real)(c * c) / (Real)n2) + rnd(rng, -1, 1);
            if (D2 < 0) D2 = 0;
            CHECK(dist2_cmp(l, p, D2) == sgn128(c * c - (Int128)D2 * n2));
        }
        // 射影・反射
        const PointD h = projection(l, p), r = reflection(l, p);
        const Real sc = std::max<Real>(1, std::max(len(l.b - l.a), dist(l.a, p)));
        CHECK(std::fabs(cross(to_d(l.b - l.a), h - to_d(l.a))) <= TOL * sc * len(l.b - l.a) * 4);
        CHECK(std::fabs(dot(to_d(l.b - l.a), to_d(p) - h)) <= TOL * sc * len(l.b - l.a) * 4);
        CHECK(nearp(midpoint(to_d(p), r), h, TOL * 4) && near(dist(h, to_d(p)), dist(l, p), TOL * 4));
        // 直線と線分・直線どうし
        const Real lt = intersects(l, t) ? 0 : std::min(dist(l, t.a), dist(l, t.b));
        CHECK(dist(l, t) == lt && dist(t, l) == lt);
        if (t.a != t.b) {
            const Line m(t.a, t.b);
            CHECK(dist(l, m) == (is_parallel(l, m) ? dist(l, m.a) : 0));
        }
    }
}

static void test_cross(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 3 : (it % 3 == 1 ? 8 : 1000000000LL);
        Segment s = rnd_seg(rng, C), t = rnd_seg(rng, C);
        if (it % 4 == 0 && C > 10) {  // 一直線で重なる
            const Point d = rndp(rng, 30), a = rndp(rng, 900000000LL);
            const long long K = 3000000;
            s = Segment(a + d * rnd(rng, -K, K), a + d * rnd(rng, -K, K)), t = Segment(a + d * rnd(rng, -K, K), a + d * rnd(rng, -K, K));
        }
        const std::vector<PointD> got = cross_points(s, t);
        if (!b_inter(s, t)) {
            CHECK(got.empty());
            continue;
        }
        const Point ds = s.b - s.a, dt = t.b - t.a;
        if (s.a != s.b && t.a != t.b && cross128(ds, dt) != 0) {
            CHECK(got.size() == 1);
            // 有理数の交点
            const Int128 den = cross128(ds, dt), num = cross128(t.a - s.a, dt);
            const PointD e((Real)((Int128)s.a.x * den + (Int128)ds.x * num) / (Real)den, (Real)((Int128)s.a.y * den + (Int128)ds.y * num) / (Real)den);
            CHECK(nearp(got[0], e, TOL));
            const std::vector<FracPoint> f = cross_point_frac(Line(s.a, s.b), Line(t.a, t.b));
            CHECK(f.size() == 1 && f[0].d > 0 && nearp(f[0].point_d(), e, TOL));
            CHECK(geo_internal::gcd128(geo_internal::gcd128(f[0].x, f[0].y), f[0].d) == 1);
            CHECK(f[0] == cross_point_frac(Line(t.b, t.a), Line(s.b, s.a))[0]);
            CHECK(f[0] == cross_point_frac(Line(s.a, s.a + ds * 3), Line(t.b, t.b - dt * 2))[0]);
            const std::vector<PointD> lp = cross_point(Line(s.a, s.b), Line(t.a, t.b));
            CHECK(lp.size() == 1 && nearp(lp[0], e, TOL));
            continue;
        }
        // 重なり：もう一方の上にある端点のうち、両端
        std::vector<Point> on;
        for (const Point& q : {s.a, s.b})
            if (b_on_seg(t.a, t.b, q)) on.push_back(q);
        for (const Point& q : {t.a, t.b})
            if (b_on_seg(s.a, s.b, q)) on.push_back(q);
        Point dir = ds != Point() ? ds : dt;
        std::sort(on.begin(), on.end(), [&](const Point& u, const Point& v) { return dot128(dir, u - s.a) < dot128(dir, v - s.a); });
        on.erase(std::unique(on.begin(), on.end()), on.end());
        std::vector<Point> e = {on.front()};
        if (on.back() != on.front()) e.push_back(on.back());
        CHECK(got.size() == e.size());
        for (size_t k = 0; k < e.size(); ++k) CHECK(got[k] == to_d(e[k]));
    }
    // 直線どうし・直線と線分
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 5 : 1000000000LL;
        const Point a = rndp(rng, C), b = rndp(rng, C);
        if (a == b) continue;
        const Line l(a, b);
        const Segment t = rnd_seg(rng, C);
        const std::vector<PointD> g = cross_point(l, t);
        if (!intersects(l, t)) CHECK(g.empty());
        else if (t.a == t.b) CHECK(g.size() == 1 && g[0] == to_d(t.a));
        else if (cross128(l.b - l.a, t.b - t.a) == 0) CHECK(g.size() == 2 && g[0] == to_d(t.a) && g[1] == to_d(t.b));
        else {
            CHECK(g.size() == 1);
            const Real sc = std::max<Real>(1, std::max(std::fabs(g[0].x), std::fabs(g[0].y)));
            CHECK(dist(to_d(l), g[0]) <= TOL * 64 * sc && dist(to_d(t), g[0]) <= TOL * 64 * sc);
        }
        CHECK(cross_point(t, l).size() == g.size());
    }
}

static void test_keys(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 4 : 999000000LL;
        Point p = rndp(rng, C), q = rndp(rng, C), r = rndp(rng, C), s = rndp(rng, C);
        if (it % 3 == 0) {  // 同じ直線の上
            const Point d = rndp(rng, 3);
            r = p + d * rnd(rng, -2, 2), s = p + d * rnd(rng, -2, 2), q = p + d * rnd(rng, -2, 2);
            if (C > 10) r = p + d * rnd(rng, -100000, 100000), s = p - d * rnd(rng, -100000, 100000);
        }
        if (p == q || r == s) continue;
        const LineKey k1 = line_key(p, q), k2 = line_key(r, s);
        const bool same = cross128(q - p, r - p) == 0 && cross128(q - p, s - p) == 0;
        CHECK((k1 == k2) == same && (k1 < k2 || k2 < k1) == !same);
        CHECK(k1 == line_key(q, p) && k1 == line_key(Line(p, q)));
        CHECK(k1.a > 0 || (k1.a == 0 && k1.b > 0));
        CHECK(std::gcd(std::gcd(std::llabs(k1.a), std::llabs(k1.b)), std::llabs(k1.c)) == 1);
        CHECK(k1.contains(p) && k1.contains(q) && k1.contains(r) == (cross128(q - p, r - p) == 0));
        CHECK((k1.dir_key() == k2.dir_key()) == (cross128(q - p, s - r) == 0) && k1.dir_key() == line_direction_key(q - p));
        CHECK(LineKeyHash()(k1) == LineKeyHash()(line_key(q, p)));
        const LineD ld = k1.line();
        const Real sc = std::max<Real>(1, len(p));
        CHECK(dist(ld, to_d(p)) <= TOL * 64 * sc && dist(ld, to_d(q)) <= TOL * 64 * std::max<Real>(1, len(q)));
        // 垂直二等分線
        const LineKey b = bisector_key(p, q);
        CHECK(b.contains(r) == (norm2_128(r - p) == norm2_128(r - q)));
        const Point mid2 = p + q;  // 中点の 2 倍
        if ((mid2.x % 2 == 0) && (mid2.y % 2 == 0)) CHECK(b.contains(Point(mid2.x / 2, mid2.y / 2)));
        const LineD bl = bisector(p, q);
        CHECK(near(dist(bl.a, to_d(p)), dist(bl.a, to_d(q)), TOL * 4) && near(dist(bl.b, to_d(p)), dist(bl.b, to_d(q)), TOL * 4));
        // 垂線・平行線
        const Line pl = perpendicular(Line(p, q), r), pa = parallel_through(Line(p, q), r);
        CHECK(pl.a == r && is_orthogonal(pl, Line(p, q)) && pa.a == r && is_parallel(pa, Line(p, q)));
    }
    // 小さい座標で垂直二等分線の上の格子点を全部調べる
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {
        const Point p = rndp(rng, 4), q = rndp(rng, 4);
        if (p == q) continue;
        const LineKey b = bisector_key(p, q);
        for (Coord x = -8; x <= 8; ++x)
            for (Coord y = -8; y <= 8; ++y) CHECK(b.contains(Point(x, y)) == (dist2(Point(x, y), p) == dist2(Point(x, y), q)));
    }
}

static void test_lattice(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const Segment s(rndp(rng, 6), rndp(rng, 6));
        // 線分の上の格子点
        Coord c = 0;
        for (Coord x = -6; x <= 6; ++x)
            for (Coord y = -6; y <= 6; ++y) c += b_on_seg(s.a, s.b, Point(x, y));
        CHECK(lattice_points_on(s) == c);
        // 内部を通るマス：t = (2k + 1) / (2M) の点がどのマスの内部にあるか（格子線の上は数えない）
        const long long M = 4000;
        std::set<Point> cells;
        const Point d = s.b - s.a;
        auto fl = [](long long a, long long b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); };
        for (long long k = 0; k < M; ++k) {
            const long long nx = s.a.x * 2 * M + d.x * (2 * k + 1), ny = s.a.y * 2 * M + d.y * (2 * k + 1);
            if (nx % (2 * M) == 0 || ny % (2 * M) == 0) continue;
            cells.insert(Point(fl(nx, 2 * M), fl(ny, 2 * M)));
        }
        CHECK(cells_crossed(s) == (Coord)cells.size());
    }
    CHECK(lattice_points_on(Segment(0, 0, 1000000000, 1000000000)) == 1000000001 && cells_crossed(Segment(0, 0, 4, 6)) == 8);
}

static void test_real(std::mt19937_64& rng, int& cases) {
    // 実数の型：整数の図形を回転・拡大・平行移動しても、判定が変わらない
    for (int it = 0; it < 10000 * TEST_SCALE; ++it, ++cases) {
        const Real th = rndr(rng, 0, 2 * GEO_PI), sc = rndr(rng, 0.01L, 100);
        const PointD sh(rndr(rng, -1000, 1000), rndr(rng, -1000, 1000));
        auto tf = [&](const Point& p) { return rotated(p, th) * sc + sh; };
        const Segment s = rnd_seg(rng, 5), t = rnd_seg(rng, 5);
        const SegmentD sd(tf(s.a), tf(s.b)), td(tf(t.a), tf(t.b));
        CHECK(intersects(sd, td) == intersects(s, t) && intersects_strict(sd, td) == intersects_strict(s, t));
        const Point p = rndp(rng, 5);
        CHECK(on_segment(sd, tf(p)) == on_segment(s, p));
        CHECK(near(dist(sd, tf(p)), dist(s, p) * sc, TOL * 1e3L));
        const std::vector<PointD> a = cross_points(s, t), b = cross_points(sd, td);
        CHECK(a.size() == b.size());
        for (size_t k = 0; k < a.size(); ++k) CHECK(nearp(rotated(a[k], th) * sc + sh, b[k], TOL * 1e3L));
        if (s.a == s.b || t.a == t.b) continue;
        const LineD l(sd.a, sd.b), m(td.a, td.b);
        CHECK(line_relation(l, m) == line_relation(Line(s.a, s.b), Line(t.a, t.b)));
        CHECK(is_parallel(l, m) == is_parallel(s, t) && is_orthogonal(l, m) == is_orthogonal(s, t));
        CHECK(nearp(projection(l, tf(p)), rotated(projection(Line(s.a, s.b), p), th) * sc + sh, TOL * 1e3L));
        CHECK(nearp(reflection(l, tf(p)), rotated(reflection(Line(s.a, s.b), p), th) * sc + sh, TOL * 1e3L));
        // 半直線・線分の位置関係
        const Ray r(s.a, s.b), q(t.a, t.b);
        const RayD rd(sd.a, sd.b), qd(td.a, td.b);
        CHECK(segment_relation(sd, td) == segment_relation(s, t) && on_ray(rd, tf(p)) == on_ray(r, p));
        CHECK(intersects(rd, qd) == intersects(r, q) && intersects(rd, td) == intersects(r, t) && intersects(rd, m) == intersects(r, Line(t.a, t.b)));
        CHECK(near(dist(rd, tf(p)), dist(r, p) * sc, TOL * 1e3L) && near(dist(rd, td), dist(r, t) * sc, TOL * 1e3L));
        const std::vector<PointD> u = cross_points(r, t), v = cross_points(rd, td);
        CHECK(u.size() == v.size());
        for (size_t k = 0; k < u.size(); ++k) CHECK(nearp(rotated(u[k], th) * sc + sh, v[k], TOL * 1e3L));
    }
}

// 半直線を十分に長い線分で置き換えた愚直な判定（座標が小さいので交点は始点から 10^6 倍の範囲に入る）
static Segment b_ray_seg(const Ray& r) { return Segment(r.a, r.a + r.dir() * 1000000); }
static Segment b_line_seg(const Line& l) { return Segment(l.a - l.dir() * 1000000, l.a + l.dir() * 1000000); }
static int b_seg_rel(const Segment& s, const Segment& t) {
    if (!b_inter(s, t)) return SEGMENT_SEPARATE;
    if (b_strict(s, t)) return SEGMENT_CROSS;
    if (s.a == s.b || t.a == t.b || cross128(s.b - s.a, t.b - t.a) != 0 || cross128(s.b - s.a, t.a - s.a) != 0) return SEGMENT_TOUCH;
    const Point d = s.b - s.a;
    Int128 a0 = dot128(d, s.a), a1 = dot128(d, s.b), b0 = dot128(d, t.a), b1 = dot128(d, t.b);
    if (a0 > a1) std::swap(a0, a1);
    if (b0 > b1) std::swap(b0, b1);
    return std::max(a0, b0) < std::min(a1, b1) ? SEGMENT_OVERLAP : SEGMENT_TOUCH;
}
static Ray rnd_ray(std::mt19937_64& rng, long long C) {
    while (true) {
        const Point a = rndp(rng, C), b = rndp(rng, C);
        if (a != b) return Ray(a, b);
    }
}

static void test_ray(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 2 : 5;
        const Ray r = rnd_ray(rng, C), q = rnd_ray(rng, C);
        const Segment s = rnd_seg(rng, C), t = rnd_seg(rng, C);
        const Point p = rndp(rng, C);
        const Line l(q.a, q.b);
        const Segment rs = b_ray_seg(r), qs = b_ray_seg(q), ls = b_line_seg(l);
        // 線分どうしの位置関係
        const int rel = segment_relation(s, t);
        CHECK(rel == b_seg_rel(s, t) && segment_relation(t, s) == rel);
        // 判定
        CHECK(on_ray(r, p) == b_on_seg(rs.a, rs.b, p));
        CHECK(intersects(r, s) == b_inter(rs, s) && intersects(s, r) == intersects(r, s));
        CHECK(intersects(r, l) == b_inter(rs, ls) && intersects(l, r) == intersects(r, l));
        CHECK(intersects(r, q) == b_inter(rs, qs) && intersects(q, r) == intersects(r, q));
        CHECK(is_parallel(r, q) == (cross128(r.dir(), q.dir()) == 0) && is_orthogonal(r, s) == (dot128(r.dir(), s.dir()) == 0));
        // 交点
        const std::vector<PointD> g1 = cross_points(r, s), e1 = cross_points(rs, s);
        CHECK(g1.size() == e1.size() && cross_points(s, r).size() == g1.size());
        for (size_t k = 0; k < g1.size(); ++k) CHECK(nearp(g1[k], e1[k], TOL));
        const std::vector<PointD> g2 = cross_point(r, l);
        if (b_inter(rs, ls) && cross128(r.dir(), l.dir()) != 0) {
            CHECK(g2.size() == 1 && nearp(g2[0], cross_point(Line(r.a, r.b), l)[0], TOL));
        } else {
            CHECK(g2.empty());
        }
        CHECK(cross_point(l, r).size() == g2.size());
        const std::vector<PointD> g3 = cross_point(r, q);
        if (!b_inter(rs, qs)) CHECK(g3.empty());
        else if (cross128(r.dir(), q.dir()) != 0) CHECK(g3.size() == 1 && nearp(g3[0], cross_point(Line(r.a, r.b), Line(q.a, q.b))[0], TOL));
        else if (r.a == q.a && dot128(r.dir(), q.dir()) < 0) CHECK(g3.size() == 1 && g3[0] == to_d(r.a));
        else CHECK(g3.empty());
        // 距離（長い線分での距離と同じ）
        CHECK(near(dist(r, p), dist(rs, p), TOL) && dist(p, r) == dist(r, p));
        CHECK(near(dist(r, s), dist(rs, s), TOL) && dist(s, r) == dist(r, s));
        CHECK(near(dist(r, q), dist(rs, qs), TOL));
        CHECK(near(dist(r, l), dist(rs, ls), TOL) && dist(l, r) == dist(r, l));
        for (int k = 0; k < 3; ++k) {
            Int128 num, den;
            b_dist2(rs, p, num, den);
            Coord D2 = (Coord)((Real)num / (Real)den) + rnd(rng, -1, 1);
            if (D2 < 0) D2 = 0;
            CHECK(dist2_cmp(r, p, D2) == sgn128(num - (Int128)D2 * den));
        }
    }
    // 大きい座標：一直線の上で始点・端点が触れる場合
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const Point d = rndp(rng, 30);
        if (d == Point()) continue;
        const Point a = rndp(rng, 900000000LL);
        const long long K = 1000000;
        const Ray r(a + d * rnd(rng, -K, K), a + d * rnd(rng, -K, K));
        if (r.a == r.b) continue;
        Segment s(a + d * rnd(rng, -K, K), a + d * rnd(rng, -K, K));
        if (it % 3 == 0) s.b = rndp(rng, 1000000000LL);
        if (it % 5 == 0) s.a = r.a;
        // 一直線なら区間の重なりで判定できる
        bool e;
        if (cross128(r.dir(), s.b - s.a) == 0 && cross128(r.dir(), s.a - r.a) == 0) {
            e = dot128(r.dir(), s.a - r.a) >= 0 || dot128(r.dir(), s.b - r.a) >= 0;
        } else {
            e = intersects(r, s);  // 一直線でない場合は下の別の式と比べる
            const Point sd = s.b - s.a;
            const Int128 den = cross128(r.dir(), sd);
            if (s.a == s.b) CHECK(e == (cross128(r.dir(), s.a - r.a) == 0 && dot128(r.dir(), s.a - r.a) >= 0));
            else if (den == 0) CHECK(!e);
            else {
                const Int128 tn = cross128(s.a - r.a, sd), un = cross128(s.a - r.a, r.dir());  // t = tn / den（半直線）、u = un / den（線分）
                const bool ok = (den > 0 ? tn >= 0 && 0 <= un && un <= den : tn <= 0 && den <= un && un <= 0);
                CHECK(e == ok);
            }
        }
        CHECK(intersects(r, s) == e);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_predicates(rng, cases);
    test_ray(rng, cases);
    test_metric(rng, cases);
    test_cross(rng, cases);
    test_keys(rng, cases);
    test_lattice(rng, cases);
    test_real(rng, cases);
    std::cout << "[line] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[line] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 2000000;
    std::vector<Segment> ss(n);
    for (auto& s : ss) s = Segment(rndp(rng, 1000000000LL), rndp(rng, 1000000000LL));
    auto t0 = std::chrono::steady_clock::now();
    long long c1 = 0;
    for (int i = 0; i + 1 < n; ++i) c1 += intersects(ss[i], ss[i + 1]);
    const long long t_int = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    Real c2 = 0;
    for (int i = 0; i + 1 < n; ++i) {
        const auto v = cross_point(Line(ss[i].a, ss[i].b), Line(ss[i + 1].a, ss[i + 1].b));
        if (!v.empty()) c2 += v[0].x;
    }
    const long long t_cp = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    Real c3 = 0;
    for (int i = 0; i + 1 < n; ++i) c3 += dist(ss[i], ss[i + 1].a);
    const long long t_d = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    long long c4 = 0;
    for (int i = 0; i < n; ++i) c4 += line_key(ss[i].a, ss[i].b).c & 7;
    const long long t_k = ms_since(t0);
    std::cout << "[line] benchmark: " << n << " 回  intersects " << t_int << " ms, cross_point " << t_cp << " ms, dist(線分, 点) " << t_d
              << " ms, line_key " << t_k << " ms  (check " << c1 + (long long)c2 % 7 + (long long)c3 % 7 + c4 << ")\n";
}

}  // namespace ln_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace ln_demo {

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
#define SHOW(expr, comment) ln_demo::show(#expr, ln_demo::to_s(expr), comment)
#define RUN(...)                   \
    ln_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("直線・線分 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Line l(a, b)（直線 a→b）/ Segment s(a, b)（線分）/ LineD / SegmentD（実数）         型は 00_geo_base\n"
              << "  is_parallel / is_orthogonal / on_line / on_segment / side                       判定（整数なら厳密）\n"
              << "  intersects(s1, s2) / intersects_strict(s1, s2) / intersects(l, s) / line_relation(l1, l2)\n"
              << "  projection(l, p) / reflection(l, p) / dist(図形, 図形) / dist2_cmp(l or s, p, D2)\n"
              << "  cross_point(l1, l2) / cross_point(l, s) / cross_points(s1, s2) / cross_point_frac(l1, l2)\n"
              << "  perpendicular(l, p) / parallel_through(l, p) / bisector(a, b)\n"
              << "  line_key(p, q) / bisector_key(p, q)（LineKey: contains / dir_key / line）/ lattice_points_on / cells_crossed\n"
              << "  segment_relation(s1, s2)（SEGMENT_SEPARATE / TOUCH / CROSS / OVERLAP）\n"
              << "  Ray r(a, b)（半直線）: on_ray / intersects / cross_point / cross_points / dist / dist2_cmp\n";

    api("【1】 on_line(l, p) / on_segment(s, p) / side(l, p) / is_parallel / is_orthogonal", "bool / int",
        "直線の上か / 線分の上か（端点を含む）/ どちら側か（左 +1・右 -1・上 0）/ 平行 / 垂直", "O(1)");
    RUN(Line l(Point(0, 0), Point(4, 2)); Segment s(Point(0, 0), Point(4, 2)));
    SHOW(on_line(l, Point(8, 4)), "直線の延長の上");
    SHOW(on_segment(s, Point(8, 4)), "線分の外");
    SHOW(on_segment(s, Point(2, 1)), "線分の上（中点）");
    SHOW(side(l, Point(0, 5)), "左側");
    SHOW(is_parallel(l, Line(Point(1, 0), Point(3, 1))), "傾き 1/2 どうし");
    SHOW(is_orthogonal(l, Segment(0, 0, -1, 2)), "内積 0");

    api("【2】 intersects(s1, s2) / intersects_strict(s1, s2) / line_relation(l1, l2)", "bool / int",
        "線分どうしが共有点を持つか（端点で触れる・重なるも含む）/ 真に交差するか / 直線の関係", "O(1)");
    SHOW(intersects(Segment(0, 0, 2, 2), Segment(0, 2, 2, 0)), "×の形");
    SHOW(intersects(Segment(0, 0, 2, 2), Segment(2, 2, 3, 0)), "端点で触れる");
    SHOW(intersects_strict(Segment(0, 0, 2, 2), Segment(2, 2, 3, 0)), "端点で触れるだけなら false");
    SHOW(intersects(Segment(0, 0, 3, 0), Segment(2, 0, 5, 0)), "一直線で重なる");
    SHOW(line_relation(Line(0, 0, 1, 1), Line(2, 2, 5, 5)) == LINE_SAME, "同じ直線");
    SHOW(line_relation(Line(0, 0, 1, 1), Line(1, 0, 2, 1)) == LINE_PARALLEL, "平行（交わらない）");

    api("【3】 cross_point(l1, l2) / cross_point(l, s) / cross_points(s1, s2)", "vector<PointD>",
        "交点（なしは空）。線分が重なるときは重なりの両端。整数の入力なら分子を 128 ビットで正確に求めてから割る", "O(1)");
    SHOW(cross_point(Line(0, 0, 4, 4), Line(0, 4, 4, 0)), "(2, 2)");
    SHOW(cross_point(Line(0, 0, 1, 1), Line(1, 0, 2, 1)), "平行なので空");
    SHOW(cross_points(Segment(0, 0, 4, 0), Segment(1, 0, 7, 0)), "重なりの両端 (1,0)〜(4,0)");
    SHOW(cross_points(Segment(0, 0, 2, 0), Segment(2, 0, 2, 5)), "端点で触れる");
    SHOW(cross_point(Line(0, 1, 1, 1), Segment(3, 0, 3, 7)), "直線 y = 1 と線分 x = 3");

    api("【4】 cross_point_frac(l1, l2)", "vector<FracPoint>（x / d, y / d。d > 0、既約）",
        "交点を誤差なしの分数で。「異なる交点の数」などを set<FracPoint> / unordered_set<FracPoint, FracPointHash> で数える", "O(log)");
    {
        const FracPoint f = cross_point_frac(Line(0, 0, 3, 1), Line(0, 1, 1, 0))[0];
        std::cout << "  cross_point_frac(Line(0,0,3,1), Line(0,1,1,0)) -> x/d = " << int128_to_string(f.x) << "/" << int128_to_string(f.d)
                  << ", y/d = " << int128_to_string(f.y) << "/" << int128_to_string(f.d) << "   // (3/4, 1/4)\n";
        SHOW(f.point_d(), "実数に");
    }

    api("【5】 projection(l, p) / reflection(l, p) / dist(...) / dist2_cmp(l or s, p, D2)", "PointD / Real / int",
        "垂線の足 / 線対称な点 / 距離（点・直線・線分のどの組でも、順は逆でもよい）/ 距離² と整数の比較（誤差なし）", "O(1)");
    SHOW(projection(Line(0, 0, 4, 0), Point(1, 3)), "垂線の足");
    SHOW(reflection(Line(0, 0, 1, 1), Point(3, 0)), "y = x について対称");
    SHOW(dist(Line(0, 0, 4, 0), Point(1, 3)), "直線と点");
    SHOW(dist(Segment(0, 0, 4, 0), Point(7, 4)), "線分と点（端点 (4, 0) が近い）");
    SHOW(dist(Segment(0, 0, 4, 0), Segment(1, 2, 3, 5)), "線分どうし");
    SHOW(dist2_cmp(Segment(0, 0, 4, 0), Point(7, 4), 25), "距離² = 25 と等しい → 0");
    SHOW(dist2_cmp(Line(0, 0, 3, 1), Point(0, 2), 3), "距離² = 3.6 > 3 → 1");
    note("ARC042 B（多角形の内側の点から辺までの最短距離）: for (i) ans = min(ans, dist(Segment(P[i], P[(i+1)%N]), p));");

    api("【6】 perpendicular(l, p) / parallel_through(l, p) / bisector(a, b)", "Line / LineD",
        "p を通る垂線・平行線（整数のまま）/ 垂直二等分線（実数）", "O(1)");
    SHOW(perpendicular(Line(0, 0, 4, 2), Point(1, 1)).b, "向き (4,2) を 90 度回した (-2, 4) だけ進んだ点");
    SHOW(bisector(Point(0, 0), Point(4, 0)).a, "中点 (2, 0) を通る");

    api("【7】 LineKey k = line_key(p, q) / bisector_key(p, q)", "LineKey{a, b, c}（a x + b y = c）",
        "同じ直線なら同じキー（gcd と符号でそろえる）。set / map / LineKeyHash で直線ごとに数える", "O(log)");
    RUN(LineKey k = line_key(Point(1, 1), Point(3, 5)));
    std::cout << "  k.a, k.b, k.c                      -> " << k.a << ", " << k.b << ", " << k.c << "        // 2x - y = 1\n";
    SHOW(k == line_key(Point(5, 9), Point(2, 3)), "同じ直線の別の 2 点");
    SHOW(k.contains(Point(0, -1)), "直線の上の点か");
    SHOW(k.dir_key(), "傾きのキー（平行な直線で同じ）");
    SHOW(bisector_key(Point(0, 0), Point(2, 2)).contains(Point(2, 0)), "(0,0) と (2,2) から等距離");
    note("ABC248 E（K 点以上を通る直線の数）は 06_point_set の count_lines_with_at_least が line_key を使う");

    api("【8】 lattice_points_on(s) / cells_crossed(s)", "Coord", "線分の上の格子点の数（gcd(|dx|, |dy|) + 1）/ 線分が内部を通るマスの数（|dx| + |dy| - gcd）", "O(log)");
    SHOW(lattice_points_on(Segment(0, 0, 6, 4)), "gcd(6, 4) + 1");
    SHOW(cells_crossed(Segment(0, 0, 6, 4)), "6 + 4 - 2");

    api("【9】 segment_relation(s1, s2)", "int（SEGMENT_SEPARATE 0 / SEGMENT_TOUCH 1 / SEGMENT_CROSS 2 / SEGMENT_OVERLAP 3）",
        "線分どうしの位置関係を 1 回で。TOUCH は共有点が 1 点だけで交差していない（端点どうし・T 字・一直線で端が触れる）", "O(1)");
    SHOW(segment_relation(Segment(0, 0, 2, 2), Segment(0, 2, 2, 0)) == SEGMENT_CROSS, "×の形");
    SHOW(segment_relation(Segment(0, 0, 2, 0), Segment(1, 0, 1, 3)) == SEGMENT_TOUCH, "T 字（端点が辺の途中に乗る）");
    SHOW(segment_relation(Segment(0, 0, 2, 0), Segment(2, 0, 5, 0)) == SEGMENT_TOUCH, "一直線で端どうしが触れる");
    SHOW(segment_relation(Segment(0, 0, 3, 0), Segment(2, 0, 5, 0)) == SEGMENT_OVERLAP, "一直線で重なる（長さ 1）");
    SHOW(segment_relation(Segment(0, 0, 1, 0), Segment(2, 0, 5, 0)) == SEGMENT_SEPARATE, "離れている");

    api("【10】 Ray r(a, b)：on_ray / intersects / cross_point / cross_points / dist / dist2_cmp", "bool / vector<PointD> / Real",
        "半直線（a から b の向きへ、a を含む）。直線・線分・半直線との判定は整数なら厳密（引数の順は逆でもよい）", "O(1)");
    RUN(Ray r(Point(0, 0), Point(1, 1)));
    SHOW(on_ray(r, Point(5, 5)), "向きの先");
    SHOW(on_ray(r, Point(-1, -1)), "始点の後ろは含まない");
    SHOW(intersects(r, Segment(3, 0, 3, 5)), "x = 3 の線分と交わる");
    SHOW(intersects(r, Segment(-3, 0, -3, -5)), "後ろの線分とは交わらない");
    SHOW(intersects(r, Line(0, 2, 1, 2)), "直線 y = 2 と交わる");
    SHOW(intersects(r, Ray(Point(4, 0), Point(4, 1))), "半直線どうし（(4, 4) で交わる）");
    SHOW(cross_point(r, Line(0, 2, 1, 2)), "交点 (2, 2)");
    SHOW(cross_points(r, Segment(-2, -2, 3, 3)), "一直線で重なる部分の両端 (0,0)〜(3,3)");
    SHOW(dist(r, Point(-3, -4)), "始点 (0, 0) が最も近い");
    SHOW(dist(r, Point(0, 2)), "途中の点 (1, 1) が最も近い（√2）");
    note("ABC016 D 型（線分が多角形の辺と何回交わるか）は intersects_strict。点から伸ばした半直線で何回横切るかは intersects(r, 辺) を数える。");
}

}  // namespace ln_demo


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
    ln_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ln_test::benchmark();    // 速度計測
    ln_demo::demo();         // 使い方デモ
    return 0;
}

