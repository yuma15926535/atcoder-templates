// =============================================================================
//  軸に平行な長方形・直方体（共通部分・位置関係・和集合の面積・k 次元の箱・象限などへの分割・傾き ±1 の直線の上の格子点）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 08 軸に平行な長方形・直方体
//  依存       : 00_geo_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - Rect(x1, y1, x2, y2) は半開 [x1, x2) × [y1, y2)（グラフ・グリッドの区間と同じ）。面積 = (x2 - x1)(y2 - y1)
//      格子点（マス）で考えると x1 ≤ x < x2 の整数 → 個数も area()。閉区間 L ≤ x ≤ R の「格子点・マス」を数えるときだけ rect_from_closed で直す
//      連続な図形（A ≤ x ≤ B, C ≤ y ≤ D の長方形の面積・重なり）はそのまま Rect(A, C, B, D)（面積 (B - A)(D - C)）。入力が A B C D の順なら rect_from_ranges
//    - x1 ≥ x2 などは空の長方形（面積 0）。座標は 10^9 まで（面積は long long、k 次元の体積は Int128）
//
//  機能一覧                                                                    計算量
//    Rect r(x1, y1, x2, y2)            r.width() / r.height() / r.area() / r.empty() / r.contains(p) / r.contains(r2)
//    rect_from_points(a, b)            2 頂点 a, b（向きは自由）を対角とする長方形                   O(1)
//    rect_from_closed(L, D, R, U)      L ≤ x ≤ R, D ≤ y ≤ U の格子点 → Rect(L, D, R + 1, U + 1)          O(1)
//    rect_from_ranges(x1, x2, y1, y2)  x の範囲・y の範囲の順で作る（= Rect(x1, y1, x2, y2)。入力が x1 x2 y1 y2 の順のとき）  O(1)
//    lattice_points_closed(r)          閉じた長方形 [x1, x2] × [y1, y2] の周を含む格子点の数          O(1)
//   閉じた長方形 [x1, x2] × [y1, y2] としての位置関係（隣り合う・角で接するの区別）
//    where_closed(r, p)                点の内外 OUTSIDE / ON_BORDER / INSIDE                          O(1)
//    intersection_dim(a, b)            共通部分の次元：-1 共有点なし / 0 角の 1 点 / 1 辺を共有（隣り合う）/ 2 面積のある重なり  O(1)
//                                      BoxK どうしも同じ（3 次元なら 2 は面を共有、3 は体積のある重なり）       O(K)
//    boundary_position(r, p)           周に沿った位置（(x1, y1) から反時計回り。周の上でなければ -1）。周の上の点を結ぶ弦の問題に  O(1)
//                                      ※ 0 ≤ x ≤ R, 0 ≤ y ≤ C の長方形の周なら Rect(0, 0, R, C) をそのまま渡す（rect_from_closed は使わない）
//    to_polygon(r)                     4 頂点の多角形（反時計回り）。03_polygon・05_circle_area の関数に渡すときに   O(1)
//    rect_fits(p, q, a, b)             p × q の長方形を、回して（斜めも可）a × b の長方形に入れられるか（辺に触れてもよい。厳密）  O(1)
//    intersect(a, b) / overlaps(a, b) / bounding(a, b)   共通部分 / 面積が正の重なりがあるか / 両方を含む最小の長方形
//    union_area(rects)                 長方形の和集合の面積（走査＋セグメント木）                     O(N log N)
//    BoxK<K> b(lo, hi)                 k 次元の箱 [lo_i, hi_i)。b.volume()（Int128。辺が 2×10^9 までなら K ≤ 4 で収まる）/ b.empty() / b.contains(p)
//    intersect(b1, b2) / box_intersection_volume(b1, b2)   箱の共通部分 / その体積（Box3 = BoxK<3> など）  O(K)
//    quadrant_split(r, c = (0, 0))     象限で 4 つに分ける：[0] x≥cx,y≥cy [1] x<cx,y≥cy [2] x<cx,y<cy [3] x≥cx,y<cy  O(1)
//    sign_grid_split(r, c = (0, 0))    x < = > cx と y < = > cy で 9 つに（番号 3·(x の側) + (y の側)。側は < = > の順に 0, 1, 2）
//    split_by_x_axis(r, cy = 0) / split_by_y_axis(r, cx = 0)   y < = > cy / x < = > cx で 3 つに          O(1)
//    diagonal_lattice(r, c, slope)     直線 y = slope·x + c（slope = 1 か -1）の上の r の中の格子点
//                                      DiagonalLattice{exist, first, last, count}（first は x が最小）         O(1)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    rc_test::random_test()  … マスを全部数える愚直な実装と突き合わせるランダムテスト
//    rc_test::benchmark()    … 速度計測
//    rc_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_RECT
#define GG_RECT
// ---- 軸に平行な長方形・直方体（幾何第二世代）  半開 [x1, x2) × [y1, y2) ----
struct Rect {
    Coord x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    constexpr Rect() {}
    constexpr Rect(Coord x1_, Coord y1_, Coord x2_, Coord y2_) : x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}
    constexpr Coord width() const { return x2 > x1 ? x2 - x1 : 0; }
    constexpr Coord height() const { return y2 > y1 ? y2 - y1 : 0; }
    constexpr Coord area() const { return width() * height(); }  // 格子点（マス）の数でもある
    constexpr bool empty() const { return x1 >= x2 || y1 >= y2; }
    constexpr bool contains(const Point& p) const { return x1 <= p.x && p.x < x2 && y1 <= p.y && p.y < y2; }
    constexpr bool contains(const Rect& r) const { return r.empty() || (x1 <= r.x1 && r.x2 <= x2 && y1 <= r.y1 && r.y2 <= y2); }
    friend constexpr bool operator==(const Rect& a, const Rect& b) { return a.x1 == b.x1 && a.y1 == b.y1 && a.x2 == b.x2 && a.y2 == b.y2; }
    friend constexpr bool operator!=(const Rect& a, const Rect& b) { return !(a == b); }
};
inline std::ostream& operator<<(std::ostream& os, const Rect& r) { return os << "Rect[" << r.x1 << ", " << r.x2 << ") x [" << r.y1 << ", " << r.y2 << ")"; }
inline Rect rect_from_points(const Point& a, const Point& b) { return Rect(std::min(a.x, b.x), std::min(a.y, b.y), std::max(a.x, b.x), std::max(a.y, b.y)); }
inline Rect rect_from_closed(Coord L, Coord D, Coord R, Coord U) { return Rect(L, D, R + 1, U + 1); }  // L ≤ x ≤ R, D ≤ y ≤ U
inline Coord lattice_points_closed(const Rect& r) {  // [x1, x2] × [y1, y2]（周を含む）の格子点の数
    return (r.x1 <= r.x2 && r.y1 <= r.y2) ? (r.x2 - r.x1 + 1) * (r.y2 - r.y1 + 1) : 0;
}
inline Rect intersect(const Rect& a, const Rect& b) {  // 共通部分（なければ空の長方形）
    return Rect(std::max(a.x1, b.x1), std::max(a.y1, b.y1), std::min(a.x2, b.x2), std::min(a.y2, b.y2));
}
inline bool overlaps(const Rect& a, const Rect& b) { return !intersect(a, b).empty(); }  // 面積が正の重なりがあるか
inline Rect bounding(const Rect& a, const Rect& b) {  // 両方を含む最小の長方形（空のものは無視）
    if (a.empty()) return b;
    if (b.empty()) return a;
    return Rect(std::min(a.x1, b.x1), std::min(a.y1, b.y1), std::max(a.x2, b.x2), std::max(a.y2, b.y2));
}

inline Rect rect_from_ranges(Coord x1, Coord x2, Coord y1, Coord y2) { return Rect(x1, y1, x2, y2); }  // x の範囲 → y の範囲の順で
// ---- 閉じた長方形 [x1, x2] × [y1, y2] としての位置関係（x1 ≤ x2, y1 ≤ y2 とする） ----
inline Where where_closed(const Rect& r, const Point& p) {  // 点の内外（周の上は ON_BORDER）
    if (p.x < r.x1 || r.x2 < p.x || p.y < r.y1 || r.y2 < p.y) return OUTSIDE;
    return (p.x == r.x1 || p.x == r.x2 || p.y == r.y1 || p.y == r.y2) ? ON_BORDER : INSIDE;
}
// 共通部分の次元：-1 共有点なし / 0 角の 1 点だけ / 1 辺の一部を共有（隣り合う）/ 2 面積のある重なり
inline int intersection_dim(const Rect& a, const Rect& b) {
    const Coord w = std::min(a.x2, b.x2) - std::max(a.x1, b.x1), h = std::min(a.y2, b.y2) - std::max(a.y1, b.y1);
    if (w < 0 || h < 0) return -1;
    return (w > 0) + (h > 0);
}
// 周に沿った位置：(x1, y1) から反時計回り（下の辺 → 右の辺 → 上の辺 → 左の辺）。周の上でなければ -1（幅・高さは正とする）
// 閉じた長方形 [x1, x2] × [y1, y2] の周。0 ≤ x ≤ R, 0 ≤ y ≤ C なら Rect(0, 0, R, C)（rect_from_closed で + 1 しない）
inline Coord boundary_position(const Rect& r, const Point& p) {
    if (where_closed(r, p) != ON_BORDER) return -1;
    const Coord W = r.x2 - r.x1, H = r.y2 - r.y1;
    if (p.y == r.y1) return p.x - r.x1;
    if (p.x == r.x2) return W + (p.y - r.y1);
    if (p.y == r.y2) return W + H + (r.x2 - p.x);
    return 2 * W + H + (r.y2 - p.y);
}
inline std::vector<Point> to_polygon(const Rect& r) { return {Point(r.x1, r.y1), Point(r.x2, r.y1), Point(r.x2, r.y2), Point(r.x1, r.y2)}; }
// 辺の長さ p, q の長方形を、回して（斜めにしてもよい）辺の長さ a, b の長方形の中に置けるか（辺に触れてもよい）。整数で厳密
// 長い辺どうし・短い辺どうしで入らないときは、斜めに置ける条件 ((a+b)/(p+q))² + ((a-b)/(p-q))² ≥ 2 で判定する
inline bool rect_fits(Coord p, Coord q, Coord a, Coord b) {
    if (p < q) std::swap(p, q);
    if (a < b) std::swap(a, b);
    if (p <= a && q <= b) return true;
    if (p <= a || q > b) return false;  // ここから p > a ≥ b ≥ q（p > q）
    const Int128 s = (Int128)(a + b) * (a + b), d = (Int128)(a - b) * (a - b), P = (Int128)(p + q) * (p + q), Q = (Int128)(p - q) * (p - q);
    return s * Q + d * P >= 2 * P * Q;
}

// 長方形の和集合の面積（x で走査し、y を座標圧縮したセグメント木で覆われた長さを持つ）
inline Coord union_area(const std::vector<Rect>& rects) {
    std::vector<Coord> ys;
    struct Ev {
        Coord x;
        int a, b, d;  // y の区間 [ys[a], ys[b]) に d を足す
    };
    for (const Rect& r : rects)
        if (!r.empty()) ys.push_back(r.y1), ys.push_back(r.y2);
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    const int m = (int)ys.size() - 1;
    if (m <= 0) return 0;
    auto idx = [&](Coord y) {  // ys の中の位置（二分探索）
        int lo = 0, hi = (int)ys.size() - 1;
        while (lo < hi) {
            const int md = (lo + hi) / 2;
            (ys[md] < y ? lo : hi) = md + (ys[md] < y);
        }
        return lo;
    };
    std::vector<Ev> ev;
    for (const Rect& r : rects)
        if (!r.empty()) {
            const int a = idx(r.y1), b = idx(r.y2);
            ev.push_back({r.x1, a, b, 1}), ev.push_back({r.x2, a, b, -1});
        }
    std::sort(ev.begin(), ev.end(), [](const Ev& p, const Ev& q) { return p.x < q.x; });
    // 節点 k：覆う y の長さ full、区間に丸ごと足された回数 cnt、覆われている長さ len（非再帰のセグメント木）
    int sz = 1;
    while (sz < m) sz <<= 1;
    struct Node {
        Coord len, full;
        int cnt;
    };
    std::vector<Node> t(2 * sz, Node{0, 0, 0});
    for (int i = 0; i < m; ++i) t[sz + i].full = ys[i + 1] - ys[i];
    for (int k = sz - 1; k >= 1; --k) t[k].full = t[2 * k].full + t[2 * k + 1].full;
    auto pull = [&](int k) { t[k].len = t[k].cnt > 0 ? t[k].full : (k >= sz ? 0 : t[2 * k].len + t[2 * k + 1].len); };
    auto update = [&](int a, int b, int d) {  // y の区間 [a, b) に d を足す
        int l = a + sz, r = b + sz;
        int k1 = l >> 1, k2 = (r - 1) >> 1;
        while (l < r) {
            if (l & 1) t[l].cnt += d, pull(l), ++l;
            if (r & 1) --r, t[r].cnt += d, pull(r);
            l >>= 1, r >>= 1;
        }
        for (; k1 >= 1; k1 >>= 1, k2 >>= 1) {  // 両端の祖先を下から（同じ深さなので並べて進める）
            pull(k1);
            if (k2 != k1) pull(k2);
        }
    };
    Coord res = 0;
    for (size_t i = 0; i < ev.size(); ++i) {
        if (i > 0) res += t[1].len * (ev[i].x - ev[i - 1].x);
        update(ev[i].a, ev[i].b, ev[i].d);
    }
    return res;
}

// ---- k 次元の箱 [lo_i, hi_i) ----
template <int K>
struct BoxK {
    std::array<Coord, K> lo{}, hi{};
    constexpr BoxK() {}
    constexpr BoxK(const std::array<Coord, K>& lo_, const std::array<Coord, K>& hi_) : lo(lo_), hi(hi_) {}
    bool empty() const {
        for (int i = 0; i < K; ++i)
            if (lo[i] >= hi[i]) return true;
        return false;
    }
    Int128 volume() const {  // 体積（空なら 0）。出力は int128_to_string か (long long)。辺が 2×10^9 までなら K ≤ 4 で Int128 に収まる
        if (empty()) return 0;
        Int128 v = 1;
        for (int i = 0; i < K; ++i) v *= hi[i] - lo[i];
        return v;
    }
    bool contains(const std::array<Coord, K>& p) const {
        for (int i = 0; i < K; ++i)
            if (p[i] < lo[i] || hi[i] <= p[i]) return false;
        return true;
    }
};
using Box3 = BoxK<3>;
template <int K>
BoxK<K> intersect(const BoxK<K>& a, const BoxK<K>& b) {
    BoxK<K> r;
    for (int i = 0; i < K; ++i) r.lo[i] = std::max(a.lo[i], b.lo[i]), r.hi[i] = std::min(a.hi[i], b.hi[i]);
    return r;
}
template <int K>
Int128 box_intersection_volume(const BoxK<K>& a, const BoxK<K>& b) { return intersect(a, b).volume(); }
// 閉じた箱 [lo_i, hi_i] として共通部分の次元（-1 共有点なし / 0 頂点 / 1 辺 / 2 面 / ... / K 体積のある重なり）
template <int K>
int intersection_dim(const BoxK<K>& a, const BoxK<K>& b) {
    int d = 0;
    for (int i = 0; i < K; ++i) {
        const Coord w = std::min(a.hi[i], b.hi[i]) - std::max(a.lo[i], b.lo[i]);
        if (w < 0) return -1;
        d += w > 0;
    }
    return d;
}

// ---- 格子点の長方形の分割（中心 c のまわり。半開なので「x = cx」は [cx, cx + 1)） ----
namespace geo_internal {
inline std::array<Rect, 3> rect_split3(const Rect& r, Coord c, bool along_x) {  // along_x：x で分ける
    std::array<Rect, 3> s{r, r, r};
    if (along_x) {
        s[0].x2 = std::min(r.x2, c), s[1].x1 = std::max(r.x1, c), s[1].x2 = std::min(r.x2, c + 1), s[2].x1 = std::max(r.x1, c + 1);
    } else {
        s[0].y2 = std::min(r.y2, c), s[1].y1 = std::max(r.y1, c), s[1].y2 = std::min(r.y2, c + 1), s[2].y1 = std::max(r.y1, c + 1);
    }
    return s;
}
}  // namespace geo_internal
// [0] x ≥ cx, y ≥ cy（第 1 象限）/ [1] x < cx, y ≥ cy / [2] x < cx, y < cy / [3] x ≥ cx, y < cy
inline std::array<Rect, 4> quadrant_split(const Rect& r, const Point& c = Point()) {
    const Coord lx = std::min(r.x2, c.x), hx = std::max(r.x1, c.x), ly = std::min(r.y2, c.y), hy = std::max(r.y1, c.y);
    return {Rect(hx, hy, r.x2, r.y2), Rect(r.x1, hy, lx, r.y2), Rect(r.x1, r.y1, lx, ly), Rect(hx, r.y1, r.x2, ly)};
}
// 番号 3·i + j：i = x の側（0: x < cx、1: x = cx、2: x > cx）、j = y の側（同様）
inline std::array<Rect, 9> sign_grid_split(const Rect& r, const Point& c = Point()) {
    const std::array<Rect, 3> xs = geo_internal::rect_split3(r, c.x, true);
    std::array<Rect, 9> res;
    for (int i = 0; i < 3; ++i) {
        const std::array<Rect, 3> ys = geo_internal::rect_split3(xs[i], c.y, false);
        for (int j = 0; j < 3; ++j) res[3 * i + j] = ys[j];
    }
    return res;
}
inline std::array<Rect, 3> split_by_x_axis(const Rect& r, Coord cy = 0) { return geo_internal::rect_split3(r, cy, false); }  // y < = > cy
inline std::array<Rect, 3> split_by_y_axis(const Rect& r, Coord cx = 0) { return geo_internal::rect_split3(r, cx, true); }   // x < = > cx

// ---- 傾き ±1 の直線の上の格子点 ----
struct DiagonalLattice {
    bool exist = false;
    Point first, last;  // x が最小・最大の格子点
    Coord count = 0;
};
// r の中（半開）の格子点のうち、直線 y = slope·x + c（slope = 1 か -1）の上にあるもの
inline DiagonalLattice diagonal_lattice(const Rect& r, Coord c, int slope) {
    // x の範囲：x1 ≤ x < x2 かつ y1 ≤ slope·x + c < y2
    Coord lo = r.x1, hi = r.x2 - 1;
    if (slope == 1) lo = std::max(lo, r.y1 - c), hi = std::min(hi, r.y2 - 1 - c);
    else lo = std::max(lo, c - (r.y2 - 1)), hi = std::min(hi, c - r.y1);
    DiagonalLattice d;
    if (r.empty() || lo > hi) return d;
    d.exist = true, d.count = hi - lo + 1;
    d.first = Point(lo, slope * lo + c), d.last = Point(hi, slope * hi + c);
    return d;
}
#endif  // GG_RECT
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
namespace rc_test {

#define RC_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[rect] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK RC_CHECK

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

static Rect rnd_rect(std::mt19937_64& rng, long long C) {
    Coord a = rnd(rng, -C, C), b = rnd(rng, -C, C), c = rnd(rng, -C, C), d = rnd(rng, -C, C);
    if (rng() % 8 != 0) {  // たいていは正しい向き
        if (a > b) std::swap(a, b);
        if (c > d) std::swap(c, d);
    }
    return Rect(a, c, b, d);
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const Rect a = rnd_rect(rng, 4), b = rnd_rect(rng, 4);
        // マスを全部数える
        Coord ca = 0, ci = 0, cu = 0;
        bool sub = true;
        for (Coord x = -5; x <= 5; ++x)
            for (Coord y = -5; y <= 5; ++y) {
                const Point p(x, y);
                const bool ia = a.contains(p), ib = b.contains(p);
                ca += ia, ci += ia && ib, cu += ia || ib;
                if (ib && !ia) sub = false;
                CHECK(intersect(a, b).contains(p) == (ia && ib));
                if (!a.empty() && !b.empty()) CHECK(!(ia || ib) || bounding(a, b).contains(p));
            }
        CHECK(a.area() == ca && intersect(a, b).area() == ci && overlaps(a, b) == (ci > 0));
        CHECK(union_area({a, b}) == cu && a.contains(b) == sub);
        CHECK(rect_from_points(Point(a.x2, a.y1), Point(a.x1, a.y2)).area() == std::llabs((a.x2 - a.x1) * (a.y2 - a.y1)));
        // 閉区間
        const Coord L = rnd(rng, -3, 3), R = rnd(rng, -3, 3), D = rnd(rng, -3, 3), U = rnd(rng, -3, 3);
        Coord cc = 0;
        for (Coord x = L; x <= R; ++x)
            for (Coord y = D; y <= U; ++y) ++cc;
        CHECK(rect_from_closed(L, D, R, U).area() == cc && lattice_points_closed(Rect(L, D, R, U)) == cc);
    }
}

static void test_union(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 12);
        std::vector<Rect> rs(n);
        for (auto& r : rs) r = rnd_rect(rng, 6);
        Coord c = 0;
        for (Coord x = -7; x <= 7; ++x)
            for (Coord y = -7; y <= 7; ++y) {
                bool in = false;
                for (const Rect& r : rs) in |= r.contains(Point(x, y));
                c += in;
            }
        CHECK(union_area(rs) == c);
    }
    // 大きい座標：座標圧縮したマスごとに数える
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 25);
        std::vector<Rect> rs(n);
        for (auto& r : rs) r = rnd_rect(rng, 1000000000LL);
        std::vector<Coord> xs, ys;
        for (const Rect& r : rs) xs.push_back(r.x1), xs.push_back(r.x2), ys.push_back(r.y1), ys.push_back(r.y2);
        std::sort(xs.begin(), xs.end()), std::sort(ys.begin(), ys.end());
        Coord c = 0;
        for (size_t i = 0; i + 1 < xs.size(); ++i)
            for (size_t j = 0; j + 1 < ys.size(); ++j) {
                if (xs[i] == xs[i + 1] || ys[j] == ys[j + 1]) continue;
                bool in = false;
                for (const Rect& r : rs) in |= r.contains(Point(xs[i], ys[j]));
                if (in) c += (xs[i + 1] - xs[i]) * (ys[j + 1] - ys[j]);
            }
        CHECK(union_area(rs) == c);
    }
}

static void test_box_split(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        Box3 a, b;
        for (int k = 0; k < 3; ++k) {
            a.lo[k] = rnd(rng, -3, 3), a.hi[k] = rnd(rng, -3, 3), b.lo[k] = rnd(rng, -3, 3), b.hi[k] = rnd(rng, -3, 3);
            if (rng() % 6) {
                if (a.lo[k] > a.hi[k]) std::swap(a.lo[k], a.hi[k]);
                if (b.lo[k] > b.hi[k]) std::swap(b.lo[k], b.hi[k]);
            }
        }
        long long va = 0, vi = 0;
        for (Coord x = -3; x <= 3; ++x)
            for (Coord y = -3; y <= 3; ++y)
                for (Coord z = -3; z <= 3; ++z) {
                    const std::array<Coord, 3> p = {x, y, z};
                    va += a.contains(p), vi += a.contains(p) && b.contains(p);
                }
        CHECK(a.volume() == va && box_intersection_volume(a, b) == vi && a.empty() == (va == 0));
        // 分割：全部の格子点がちょうど 1 つの部分に、正しい側に入る
        const Rect r = rnd_rect(rng, 4);
        const Point c = rndp(rng, 3);
        const auto q = quadrant_split(r, c);
        const auto g = sign_grid_split(r, c);
        const auto sx = split_by_x_axis(r, c.y), sy = split_by_y_axis(r, c.x);
        for (Coord x = -5; x <= 5; ++x)
            for (Coord y = -5; y <= 5; ++y) {
                const Point p(x, y);
                const bool in = r.contains(p);
                int k4 = 0, k9 = 0, k3x = 0, k3y = 0;
                for (int t = 0; t < 4; ++t)
                    if (q[t].contains(p)) {
                        ++k4;
                        CHECK(t == (x >= c.x ? (y >= c.y ? 0 : 3) : (y >= c.y ? 1 : 2)));
                    }
                for (int t = 0; t < 9; ++t)
                    if (g[t].contains(p)) {
                        ++k9;
                        const int i = x < c.x ? 0 : (x == c.x ? 1 : 2), j = y < c.y ? 0 : (y == c.y ? 1 : 2);
                        CHECK(t == 3 * i + j);
                    }
                for (int t = 0; t < 3; ++t) {
                    if (sx[t].contains(p)) {
                        ++k3x;
                        CHECK(t == (y < c.y ? 0 : (y == c.y ? 1 : 2)));
                    }
                    if (sy[t].contains(p)) {
                        ++k3y;
                        CHECK(t == (x < c.x ? 0 : (x == c.x ? 1 : 2)));
                    }
                }
                CHECK(k4 == in && k9 == in && k3x == in && k3y == in);
            }
        // 傾き ±1 の直線の上の格子点
        for (int s : {1, -1}) {
            const Coord cc = rnd(rng, -8, 8);
            std::vector<Point> v;
            for (Coord x = -5; x <= 5; ++x) {
                const Point p(x, s * x + cc);
                if (r.contains(p)) v.push_back(p);
            }
            const DiagonalLattice d = diagonal_lattice(r, cc, s);
            CHECK(d.exist == !v.empty() && d.count == (Coord)v.size());
            if (!v.empty()) CHECK(d.first == v.front() && d.last == v.back());
        }
    }
}

static void test_closed(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        auto rr = [&](long long C) {
            Coord a = rnd(rng, -C, C), b = rnd(rng, -C, C), c = rnd(rng, -C, C), d = rnd(rng, -C, C);
            if (a > b) std::swap(a, b);
            if (c > d) std::swap(c, d);
            return Rect(a, c, b, d);
        };
        const Rect a = rr(4), b = rr(4);
        // 閉じた長方形の共通部分を 0.5 刻みの点で調べる（座標を 2 倍）
        std::set<Coord> xs, ys;
        for (Coord X = -10; X <= 10; ++X)
            for (Coord Y = -10; Y <= 10; ++Y)
                if (2 * a.x1 <= X && X <= 2 * a.x2 && 2 * a.y1 <= Y && Y <= 2 * a.y2 && 2 * b.x1 <= X && X <= 2 * b.x2 && 2 * b.y1 <= Y && Y <= 2 * b.y2)
                    xs.insert(X), ys.insert(Y);
        const int e = xs.empty() ? -1 : (xs.size() > 1) + (ys.size() > 1);
        CHECK(intersection_dim(a, b) == e && intersection_dim(b, a) == e);
        const Point p = rndp(rng, 5);
        const bool in = a.x1 <= p.x && p.x <= a.x2 && a.y1 <= p.y && p.y <= a.y2;
        const bool strict = a.x1 < p.x && p.x < a.x2 && a.y1 < p.y && p.y < a.y2;
        CHECK(where_closed(a, p) == (strict ? INSIDE : (in ? ON_BORDER : OUTSIDE)));
        CHECK(rect_from_ranges(a.x1, a.x2, a.y1, a.y2) == a);
        const std::vector<Point> poly = to_polygon(a);
        CHECK(poly.size() == 4 && cross128(poly[1] - poly[0], poly[2] - poly[0]) + cross128(poly[2] - poly[0], poly[3] - poly[0]) == 2 * (Int128)a.area());
        // 周に沿った位置：(x1, y1) から反時計回りに 1 ずつ歩いた順
        if (a.x1 < a.x2 && a.y1 < a.y2) {
            std::vector<Point> walk;
            Point q(a.x1, a.y1);
            const Point dirs[4] = {Point(1, 0), Point(0, 1), Point(-1, 0), Point(0, -1)};
            for (int k = 0; k < 4; ++k) {
                const Coord len = k % 2 == 0 ? a.x2 - a.x1 : a.y2 - a.y1;
                for (Coord t = 0; t < len; ++t) walk.push_back(q), q += dirs[k];
            }
            for (Coord X = -5; X <= 5; ++X)
                for (Coord Y = -5; Y <= 5; ++Y) {
                    const auto f = std::find(walk.begin(), walk.end(), Point(X, Y));
                    CHECK(boundary_position(a, Point(X, Y)) == (f == walk.end() ? -1 : (Coord)(f - walk.begin())));
                }
        }
        // 3 次元の箱
        Box3 A, B;
        for (int i = 0; i < 3; ++i) {
            A.lo[i] = rnd(rng, -2, 2), A.hi[i] = A.lo[i] + rnd(rng, 0, 3);
            B.lo[i] = rnd(rng, -2, 2), B.hi[i] = B.lo[i] + rnd(rng, 0, 3);
        }
        int dim = 0;
        bool none = false;
        for (int i = 0; i < 3; ++i) {
            std::set<Coord> vs;
            for (Coord X = -10; X <= 12; ++X)
                if (2 * A.lo[i] <= X && X <= 2 * A.hi[i] && 2 * B.lo[i] <= X && X <= 2 * B.hi[i]) vs.insert(X);
            if (vs.empty()) none = true;
            dim += vs.size() > 1;
        }
        CHECK(intersection_dim(A, B) == (none ? -1 : dim));
    }
}

static void test_fits(std::mt19937_64& rng, int& cases) {
    // 角 θ を細かく動かして、回した長方形の外接長方形 (p cos θ + q sin θ, p sin θ + q cos θ) が a × b に入るかを調べる
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const Coord p = rnd(rng, 1, 30), q = rnd(rng, 1, 30), a = rnd(rng, 1, 30), b = rnd(rng, 1, 30);
        bool e = false;
        const int M = 4000;
        for (int k = 0; k <= M && !e; ++k) {
            const Real t = GEO_PI / 2 * k / M, c = std::cos(t), sn = std::sin(t);
            e = (p * c + q * sn <= a + 1e-9L && p * sn + q * c <= b + 1e-9L);
        }
        const bool f = rect_fits(p, q, a, b);
        if (f != e) {  // 境界ぎりぎりは細かくして確かめる
            const int M2 = 400000;
            e = false;
            for (int k = 0; k <= M2 && !e; ++k) {
                const Real t = GEO_PI / 2 * k / M2, c = std::cos(t), sn = std::sin(t);
                e = (p * c + q * sn <= a + 1e-9L && p * sn + q * c <= b + 1e-9L);
            }
        }
        CHECK(f == e && rect_fits(q, p, b, a) == f && rect_fits(p, q, b, a) == f);
    }
    // 大きい値（細い長方形を正方形の対角線の向きに）：対角線 8×10^8·√2 ≈ 1.13×10^9 には入り、7×10^8·√2 ≈ 0.99×10^9 には入らない
    CHECK(rect_fits(1000000000, 1, 1000000000, 1) && rect_fits(1000000000, 1, 800000000, 800000000) && !rect_fits(1000000000, 1, 700000000, 700000000));
    CHECK(rect_fits(10, 1, 8, 8) && !rect_fits(10, 1, 7, 7) && rect_fits(1, 1, 1, 1) && !rect_fits(2, 1, 1, 1));
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_basic(rng, cases);
    test_fits(rng, cases);
    test_closed(rng, cases);
    test_union(rng, cases);
    test_box_split(rng, cases);
    std::cout << "[rect] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[rect] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 200000;
    std::vector<Rect> rs(n);
    for (auto& r : rs) r = rnd_rect(rng, 1000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    const Coord u = union_area(rs);
    const long long t1 = ms_since(t0);
    std::cout << "[rect] benchmark: union_area N = " << n << " " << t1 << " ms  (check " << u % 7 << ")\n";
}

}  // namespace rc_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace rc_demo {

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
#define SHOW(expr, comment) rc_demo::show(#expr, rc_demo::to_s(expr), comment)
#define RUN(...)                   \
    rc_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static std::string to_s(const Rect& r) {
    std::ostringstream os;
    os << r;
    return os.str();
}
static void demo() {
    banner("軸に平行な長方形・直方体 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧（Rect は半開 [x1, x2) × [y1, y2)）\n"
              << "  Rect r(x1, y1, x2, y2);  r.area() / r.empty() / r.contains(p or r2) / r.width() / r.height()\n"
              << "  rect_from_points(a, b) / rect_from_closed(L, D, R, U) / lattice_points_closed(r)\n"
              << "  intersect(a, b) / overlaps(a, b) / bounding(a, b) / union_area(rects)\n"
              << "  BoxK<K>（Box3）：volume / empty / contains、intersect / box_intersection_volume\n"
              << "  quadrant_split / sign_grid_split / split_by_x_axis / split_by_y_axis / diagonal_lattice\n"
              << "  rect_from_ranges(x1, x2, y1, y2)（入力が x1 x2 y1 y2 の順）/ to_polygon(r) / rect_fits(p, q, a, b)（回して入るか）\n"
              << "  閉じた長方形として: where_closed(r, p) / intersection_dim(a, b)（BoxK も）/ boundary_position(r, p)\n";

    api("【1】 Rect(x1, y1, x2, y2) / intersect / overlaps / union_area", "Rect / bool / Coord",
        "半開の長方形。面積はマスの数と同じ。和集合の面積は走査＋セグメント木", "O(1) / O(N log N)");
    RUN(Rect a(0, 0, 4, 3), b(2, 1, 6, 5));
    SHOW(a.area(), "4 × 3");
    SHOW(to_s(intersect(a, b)), "共通部分");
    SHOW(intersect(a, b).area(), "");
    SHOW(overlaps(a, Rect(4, 0, 5, 1)), "辺で接するだけ → 重なりなし");
    SHOW(union_area({a, b}), "12 + 16 - 4");
    SHOW(a.contains(Point(4, 0)), "x = 4 は含まない（半開）");
    note("ABC318 B（シートの重なり）: union_area でそのまま。N が小さければ 2 次元のマスを塗ってもよい");
    SHOW(rect_from_closed(1, 1, 3, 2).area(), "1 ≤ x ≤ 3, 1 ≤ y ≤ 2 の格子点は 6 個");
    SHOW(lattice_points_closed(Rect(0, 0, 2, 2)), "[0, 2] × [0, 2] の周を含む格子点は 9 個");

    api("【2】 BoxK<K> / Box3 / box_intersection_volume(b1, b2)", "Int128",
        "k 次元の箱 [lo_i, hi_i) の体積・共通部分（体積は Int128。小さいなら (long long) で出力）", "O(K)");
    RUN(Box3 c1({0, 0, 0}, {7, 7, 7}), c2({3, 0, 0}, {10, 7, 7}));
    SHOW(c1.volume(), "7³");
    SHOW(box_intersection_volume(c1, c2), "4 × 7 × 7");
    SHOW(intersect(c1, intersect(c2, Box3({0, 5, 0}, {7, 12, 7}))).volume(), "3 つの共通部分");
    note("ABC343 E（7 × 7 × 7 の立方体 3 つ）: 位置を全探索して、ちょうど 1 / 2 / 3 個に含まれる体積を包除原理で");

    api("【3】 quadrant_split / sign_grid_split / split_by_x_axis / split_by_y_axis", "array<Rect, 4 / 9 / 3>",
        "格子点の長方形を、中心 c（既定は原点）の象限・符号・軸の上下で分ける（半開なので「x = cx」は [cx, cx + 1)）", "O(1)");
    RUN(Rect r(-3, -2, 5, 6));
    RUN(auto q = quadrant_split(r));
    for (int t = 0; t < 4; ++t) std::cout << "  q[" << t << "]                               -> " << to_s(q[t]) << "  （格子点 " << q[t].area() << "）\n";
    SHOW(sign_grid_split(r)[4].area(), "(0, 0) の 1 点");
    SHOW(to_s(split_by_x_axis(r)[2]), "y > 0");

    api("【4】 diagonal_lattice(r, c, slope)", "DiagonalLattice{exist, first, last, count}",
        "直線 y = x + c（slope = 1）/ y = -x + c（slope = -1）の上の、r の中の格子点（first は x が最小）", "O(1)");
    RUN(DiagonalLattice d = diagonal_lattice(rect_from_closed(-2, -1, 4, 3), 1, 1));
    std::cout << "  d.exist, d.count, d.first, d.last  -> " << d.exist << ", " << d.count << ", " << d.first << ", " << d.last << "   // -2 ≤ x ≤ 4, -1 ≤ y ≤ 3 で y = x + 1\n";
    SHOW(diagonal_lattice(rect_from_closed(-2, -1, 4, 3), 10, -1).exist, "y = -x + 10 は通らない");
    note("45 度回したマンハッタン距離の問題（ひし形の辺の上の格子点など）で使う");

    api("【5】 where_closed(r, p) / intersection_dim(a, b) / boundary_position(r, p) / to_polygon(r) / rect_from_ranges / rect_fits", "Where / int / Coord / vector<Point> / bool",
        "閉じた長方形 [x1, x2] × [y1, y2] としての位置関係。隣り合う（辺を共有）・角で接するを区別する", "O(1)");
    RUN(Rect sq(0, 0, 2, 2));
    SHOW(where_closed(sq, Point(2, 1)) == ON_BORDER, "右の辺の上");
    SHOW(intersection_dim(sq, Rect(2, 0, 4, 2)), "辺 x = 2 を共有 → 1（隣り合う）");
    SHOW(intersection_dim(sq, Rect(2, 2, 3, 3)), "角 (2, 2) だけ → 0");
    SHOW(intersection_dim(sq, Rect(1, 1, 3, 3)), "面積のある重なり → 2");
    SHOW(intersection_dim(sq, Rect(3, 0, 4, 2)), "離れている → -1");
    SHOW(intersection_dim(Box3({0, 0, 0}, {1, 1, 1}), Box3({1, 0, 0}, {2, 1, 1})), "立方体が面を共有 → 2（ABC312 E）");
    SHOW(boundary_position(Rect(0, 0, 4, 3), Point(4, 1)), "下の辺 4 + 右の辺を 1 上がる");
    SHOW(boundary_position(Rect(0, 0, 4, 3), Point(1, 1)), "周の上でない");
    SHOW(to_polygon(Rect(0, 0, 4, 3)), "4 頂点（反時計回り）");
    SHOW(rect_from_ranges(1, 5, 2, 3), "入力が x1 x2 y1 y2 の順（ABC318 B など）");
    SHOW(rect_fits(10, 1, 8, 8), "rect_fits: 10 × 1 は 8 × 8 に斜め（45 度）なら入る");
    SHOW(rect_fits(10, 1, 7, 7), "7 × 7 には入らない");
    note("ARC029 B（紙を封筒に入れられるか）: rect_fits(A, B, C, D)");
    note("ABC312 E（直方体どうしが面を共有する数）: マスに分けて数えるのが速いが、少数なら intersection_dim == 2 で判定できる");
    note("ARC076 E 型: Rect board(0, 0, R, C);（rect_from_closed は使わない）両端が周の上にある線分を boundary_position の弦にして");
    note("    05_circle の count_crossing_chords(弦) == 0 か（08_rect と 05_circle を両方貼ってよい）");
}

}  // namespace rc_demo


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
    rc_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    rc_test::benchmark();    // 速度計測
    rc_demo::demo();         // 使い方デモ
    return 0;
}

