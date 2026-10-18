// =============================================================================
//  凸包（凸包・一直線の点を残す凸包・上下の凸包・最も遠い 2 点・凸多角形の内外判定 O(log N)）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 04 凸包・凸多角形
//  依存       : 00_geo_base, 02_line, 03_polygon（同梱済み。凸包の面積・ピックの定理・多角形の位置関係もこのファイルだけで使える）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - 凸包は反時計回りで、(x, y) が最小の点から始まる。同じ点は 1 つ（番号の小さい方）
//    - すべて整数で厳密（|x|, |y| ≤ 10^9）
//
//  機能一覧                                                                    計算量
//    convex_hull(ps, keep_collinear = false)      凸包の頂点（keep_collinear = true なら辺の上の点も）   O(N log N)
//    convex_hull_idx(ps, keep_collinear = false)  同じものを番号で
//        ※ 全部が一直線なら：false → 両端の 2 点、true → 端から端へ全部の点。1 点なら 1 点
//    lower_hull(ps) / upper_hull(ps)              下側 / 上側の凸包（x の小さい順。一直線の点なし）       O(N log N)
//    farthest_pair(ps)                            最も遠い 2 点 FarthestPair{i, j, d2}（番号は ps の中、i < j）  O(N log N)
//    convex_diameter(h)                           凸多角形 h の最も遠い 2 頂点（番号は h の中。回転キャリパー）  O(N)
//    convex_contains(h, p)                        凸多角形 h の内外（OUTSIDE / ON_BORDER / INSIDE）       O(log N)
//    convex_contains(h, Q)                        多角形 Q（凸でなくてよい）が凸多角形 h（周を含む）に入っているか  O(M log N)
//        ※ h は反時計回りで h[0] が角（180 度でない頂点）であること。convex_hull の結果ならそうなっている
//          入力で与えられた凸多角形をそのまま使うなら convex_hull(P, true) を通す（向き・始点・一直線の頂点をそろえる）
//    lattice_points_in_convex(h)                  凸多角形 h の内部と周の上の格子点の数（頂点 0〜2 個も可）  O(N log)
//        ※ h は convex_hull の結果（反時計回り。一直線の頂点があってもよい。頂点 0〜2 個の退化した形も可）
//    （同梱の 03_polygon）area2(h) / area(h) / interior_lattice_points(h) / boundary_lattice_points(h) / contains / is_convex ...
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    ch_test::random_test()  … 愚直な実装（三角形に含まれるか・支える直線・全組の距離）と突き合わせるランダムテスト
//    ch_test::benchmark()    … 速度計測
//    ch_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_POLYGON
#define GG_POLYGON
// ---- 多角形・三角形（幾何第二世代）  面積・内外判定・凸・単純・重心・格子点・三角形の種類 ----
enum TriangleKind : int { TRIANGLE_DEGENERATE = 0, TRIANGLE_ACUTE = 1, TRIANGLE_RIGHT = 2, TRIANGLE_OBTUSE = 3 };
namespace geo_internal {
template <class T>
int poly_arg_cmp(const PointT<T>& a, const PointT<T>& b) {  // 偏角の比較（凸判定で一周を数えるため）
    auto half = [](const PointT<T>& v) { return (v.y > 0 || (v.y == 0 && v.x > 0)) ? 0 : 1; };
    const int ha = half(a), hb = half(b);
    if (ha != hb) return ha < hb ? -1 : 1;
    return -sgn_cross(a, b);
}
template <class T>
bool poly_seg_meet(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c, const PointT<T>& d) {  // 線分 ab と cd が共有点を持つか
    if (a == b) return ccw_detail(c, d, a) == ON_SEGMENT;
    if (c == d) return ccw_detail(a, b, c) == ON_SEGMENT;
    return ccw_detail(a, b, c) * ccw_detail(a, b, d) <= 0 && ccw_detail(c, d, a) * ccw_detail(c, d, b) <= 0;
}
}  // namespace geo_internal

// 符号付き面積の 2 倍（反時計回りなら正）。整数の点なら厳密（途中は 128 ビット）
template <class T>
T area2(const std::vector<PointT<T>>& poly) {
    const int n = (int)poly.size();
    if (n < 3) return 0;
    if constexpr (std::is_integral_v<T>) {
        Int128 s = 0;
        for (int i = 1; i + 1 < n; ++i) s += cross128(poly[i] - poly[0], poly[i + 1] - poly[0]);
        return (T)s;
    } else {
        T s = 0;
        for (int i = 1; i + 1 < n; ++i) s += cross(poly[i] - poly[0], poly[i + 1] - poly[0]);
        return s;
    }
}
template <class T>
Real area(const std::vector<PointT<T>>& poly) {
    const T s = area2(poly);
    return (Real)(s < 0 ? -s : s) / 2;
}
template <class T>
bool is_ccw(const std::vector<PointT<T>>& poly) {
    const T s = area2(poly);
    if constexpr (std::is_integral_v<T>) return s > 0;
    else {
        Real sc = 0;
        for (const auto& p : poly) sc = std::max(sc, norm2(p - poly[0]));
        return real_sign(s, sc) > 0;
    }
}
template <class T>
std::vector<PointT<T>> make_ccw(std::vector<PointT<T>> poly) {  // 時計回りなら逆順にする
    if (area2(poly) < 0) std::reverse(poly.begin(), poly.end());
    return poly;
}
// 点 p が多角形の内側か（OUTSIDE / ON_BORDER / INSIDE）。凸でなくてよい（単純多角形）
template <class T>
Where contains(const std::vector<PointT<T>>& poly, const PointT<T>& p) {
    const int n = (int)poly.size();
    bool in = false;
    for (int i = 0; i < n; ++i) {
        const PointT<T>& a = poly[i];
        const PointT<T>& b = poly[i + 1 == n ? 0 : i + 1];
        if (ccw_detail(a, b, p) == ON_SEGMENT) return ON_BORDER;
        PointT<T> u = a - p, v = b - p;
        if (u.y > v.y) std::swap(u, v);
        if (u.y <= 0 && 0 < v.y && geo_internal::sgn_cross(u, v) > 0) in = !in;
    }
    return in ? INSIDE : OUTSIDE;
}
// 凸多角形か。strict = false なら 180 度の頂点（一直線に並ぶ 3 点）と、続けて同じ点を許す
template <class T>
bool is_convex(const std::vector<PointT<T>>& poly_in, bool strict = false) {
    std::vector<PointT<T>> poly;
    for (const auto& p : poly_in)
        if (poly.empty() || poly.back() != p) poly.push_back(p);
    while (poly.size() > 1 && poly.back() == poly.front()) poly.pop_back();
    if (strict && poly.size() != poly_in.size()) return false;
    const int n = (int)poly.size();
    if (n < 3) return false;
    int dir = 0;
    for (int i = 0; i < n; ++i) {
        const PointT<T>&a = poly[i], &b = poly[(i + 1) % n], &c = poly[(i + 2) % n];
        const int s = ccw(a, b, c);
        if (s == 0) {
            if (strict || ccw_detail(a, b, c) != ONLINE_FRONT) return false;  // 折り返し
            continue;
        }
        if (dir == 0) dir = s;
        else if (s != dir) return false;
    }
    if (dir == 0) return false;
    int wraps = 0;  // 辺の向きがちょうど一周するか（星形五角形などを除く）
    for (int i = 0; i < n; ++i) {
        const int c = geo_internal::poly_arg_cmp(poly[(i + 1) % n] - poly[i], poly[(i + 2) % n] - poly[(i + 1) % n]);
        if (dir > 0 ? c > 0 : c < 0) ++wraps;
    }
    return wraps == 1;
}
// 単純多角形か（長さ 0 の辺がなく、隣り合う辺は共有する頂点だけで接し、隣り合わない辺は交わらない）
template <class T>
bool is_simple(const std::vector<PointT<T>>& poly) {
    const int n = (int)poly.size();
    if (n < 3) return false;
    for (int i = 0; i < n; ++i)
        if (poly[i] == poly[(i + 1) % n]) return false;
    for (int i = 0; i < n; ++i) {
        const PointT<T>&a = poly[i], &b = poly[(i + 1) % n];
        for (int j = i + 1; j < n; ++j) {
            const PointT<T>&c = poly[j], &d = poly[(j + 1) % n];
            if (j == i + 1) {  // 辺 ab と bd（共有点 b）。d が直線 ab の上で a の側にあると重なる
                if (ccw(a, b, d) == 0 && ccw_detail(a, b, d) != ONLINE_FRONT) return false;
                continue;
            }
            if (i == 0 && j == n - 1) {  // 辺 ca と ab（共有点 a）
                if (ccw(c, a, b) == 0 && ccw_detail(c, a, b) != ONLINE_FRONT) return false;
                continue;
            }
            if (geo_internal::poly_seg_meet(a, b, c, d)) return false;
        }
    }
    return true;
}
template <class T>
Real perimeter(const std::vector<PointT<T>>& poly) {
    const int n = (int)poly.size();
    Real s = 0;
    for (int i = 0; i < n && n >= 2; ++i) s += dist(poly[i], poly[(i + 1) % n]);
    return s;
}
// 重心（面積の重心）。面積が 0 なら頂点の平均
template <class T>
PointD centroid(const std::vector<PointT<T>>& poly) {
    const int n = (int)poly.size();
    if (n == 0) return PointD();
    const T a2 = area2(poly);
    bool zero;
    if constexpr (std::is_integral_v<T>) zero = (a2 == 0);
    else zero = (std::fabs(a2) == 0);
    if (zero) {
        Real x = 0, y = 0;
        for (const auto& p : poly) x += (Real)p.x, y += (Real)p.y;
        return PointD(x / n, y / n);
    }
    if constexpr (std::is_integral_v<T>) {
        Int128 sx = 0, sy = 0;
        for (int i = 1; i + 1 < n; ++i) {
            const Point u = poly[i] - poly[0], v = poly[i + 1] - poly[0];
            const Int128 c = cross128(u, v);
            sx += c * (u.x + v.x), sy += c * (u.y + v.y);
        }
        return PointD((Real)poly[0].x + (Real)sx / ((Real)a2 * 3), (Real)poly[0].y + (Real)sy / ((Real)a2 * 3));
    } else {
        Real sx = 0, sy = 0;
        for (int i = 1; i + 1 < n; ++i) {
            const PointD u = poly[i] - poly[0], v = poly[i + 1] - poly[0];
            const Real c = cross(u, v);
            sx += c * (u.x + v.x), sy += c * (u.y + v.y);
        }
        return PointD(poly[0].x + sx / (a2 * 3), poly[0].y + sy / (a2 * 3));
    }
}
// 周の上の格子点の数（頂点は格子点）
inline Coord boundary_lattice_points(const std::vector<Point>& poly) {
    const int n = (int)poly.size();
    if (n == 1) return 1;
    Coord b = 0;
    for (int i = 0; i < n; ++i) {
        const Point d = poly[(i + 1) % n] - poly[i];
        b += std::gcd(d.x < 0 ? -d.x : d.x, d.y < 0 ? -d.y : d.y);
    }
    return b;
}
// 内部（周を除く）の格子点の数。ピックの定理 S = I + B / 2 - 1 より I = (2S - B + 2) / 2（単純多角形、面積 > 0）
inline Coord interior_lattice_points(const std::vector<Point>& poly) {
    Coord a2 = area2(poly);
    if (a2 < 0) a2 = -a2;
    return (a2 - boundary_lattice_points(poly) + 2) / 2;
}
template <class T>
T triangle_area2(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c) { return cross(b - a, c - a); }
// 三角形の種類（整数なら厳密）。3 点が一直線（同じ点を含む）なら TRIANGLE_DEGENERATE
template <class T>
int triangle_kind(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c) {
    if (ccw(a, b, c) == 0) return TRIANGLE_DEGENERATE;
    const int s0 = geo_internal::sgn_dot(b - a, c - a), s1 = geo_internal::sgn_dot(a - b, c - b), s2 = geo_internal::sgn_dot(a - c, b - c);
    if (s0 == 0 || s1 == 0 || s2 == 0) return TRIANGLE_RIGHT;
    if (s0 < 0 || s1 < 0 || s2 < 0) return TRIANGLE_OBTUSE;
    return TRIANGLE_ACUTE;
}
template <class T>
Where in_triangle(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c, const PointT<T>& p) {
    return contains(std::vector<PointT<T>>{a, b, c}, p);
}

// ---- 多角形と線分・直線・半直線・多角形の位置関係（多角形は周を含む閉じた領域。単純多角形。整数なら厳密） ----
namespace geo_internal {
// 2 倍した点 (X, Y) の内外（辺の中点のような半整数の点を 128 ビットで厳密に）
inline Where contains_doubled(const std::vector<Point>& poly, Int128 X, Int128 Y) {
    const int n = (int)poly.size();
    bool in = false;
    for (int i = 0; i < n; ++i) {
        const Point& a = poly[i];
        const Point& b = poly[i + 1 == n ? 0 : i + 1];
        Int128 ux = 2 * (Int128)a.x - X, uy = 2 * (Int128)a.y - Y, vx = 2 * (Int128)b.x - X, vy = 2 * (Int128)b.y - Y;
        if (ux * vy - uy * vx == 0 && ux * vx + uy * vy <= 0) return ON_BORDER;
        if (uy > vy) std::swap(ux, vx), std::swap(uy, vy);
        if (uy <= 0 && 0 < vy && ux * vy - uy * vx > 0) in = !in;
    }
    return in ? INSIDE : OUTSIDE;
}
template <class T>
SegmentT<T> poly_edge(const std::vector<PointT<T>>& poly, int i) { return SegmentT<T>(poly[i], poly[i + 1 == (int)poly.size() ? 0 : i + 1]); }
}  // namespace geo_internal
// 線分 s が多角形（周を含む）に入っているか。凸でなくてよい。O(N log N)
template <class T>
bool contains(const std::vector<PointT<T>>& poly, const SegmentT<T>& s) {
    if (contains(poly, s.a) == OUTSIDE || contains(poly, s.b) == OUTSIDE) return false;
    if (s.a == s.b) return true;
    const int n = (int)poly.size();
    std::vector<PointT<T>> pts = {s.a, s.b};  // s の上の「区切りの点」：両端と、s の上にある頂点
    for (int i = 0; i < n; ++i) {
        if (intersects_strict(s, geo_internal::poly_edge(poly, i))) return false;  // 辺を横切ると外に出る
        if (on_segment(s, poly[i])) pts.push_back(poly[i]);
    }
    const PointT<T> d = s.dir();
    std::sort(pts.begin(), pts.end(), [&](const PointT<T>& p, const PointT<T>& q) { return dot(d, p - s.a) < dot(d, q - s.a); });
    for (size_t k = 0; k + 1 < pts.size(); ++k) {  // 区切りの間は内外が変わらないので、中点で調べる
        if (pts[k] == pts[k + 1]) continue;
        if constexpr (std::is_integral_v<T>) {
            if (geo_internal::contains_doubled(poly, (Int128)pts[k].x + pts[k + 1].x, (Int128)pts[k].y + pts[k + 1].y) == OUTSIDE) return false;
        } else {
            if (contains(poly, midpoint(pts[k], pts[k + 1])) == OUTSIDE) return false;
        }
    }
    return true;
}
// 線分・直線・半直線と多角形（周を含む）が共有点を持つか。O(N)
template <class T>
bool intersects(const std::vector<PointT<T>>& poly, const SegmentT<T>& s) {
    if (poly.empty()) return false;
    if (contains(poly, s.a) != OUTSIDE) return true;
    for (int i = 0; i < (int)poly.size(); ++i)
        if (intersects(s, geo_internal::poly_edge(poly, i))) return true;
    return false;
}
template <class T>
bool intersects(const std::vector<PointT<T>>& poly, const LineT<T>& l) {
    bool pos = false, neg = false;
    for (const auto& p : poly) {
        const int c = ccw(l.a, l.b, p);
        if (c == 0) return true;
        (c > 0 ? pos : neg) = true;
    }
    return pos && neg;
}
template <class T>
bool intersects(const std::vector<PointT<T>>& poly, const RayT<T>& r) {
    if (poly.empty()) return false;
    if (contains(poly, r.a) != OUTSIDE) return true;
    for (int i = 0; i < (int)poly.size(); ++i)
        if (intersects(r, geo_internal::poly_edge(poly, i))) return true;
    return false;
}
// 多角形 Q が多角形 P（周を含む）に入っているか（Q の周が全部 P に入っていればよい）。O(N M log N)
template <class T>
bool contains(const std::vector<PointT<T>>& P, const std::vector<PointT<T>>& Q) {
    if (Q.empty()) return true;
    for (int i = 0; i < (int)Q.size(); ++i)
        if (!contains(P, geo_internal::poly_edge(Q, i))) return false;
    return true;
}
// 2 つの多角形（周を含む）が共有点を持つか（辺が交わる・接する・片方がもう片方に入っている）。O(N M)
template <class T>
bool intersects(const std::vector<PointT<T>>& P, const std::vector<PointT<T>>& Q) {
    if (P.empty() || Q.empty()) return false;
    if (contains(P, Q[0]) != OUTSIDE || contains(Q, P[0]) != OUTSIDE) return true;
    for (int i = 0; i < (int)P.size(); ++i)
        for (int j = 0; j < (int)Q.size(); ++j)
            if (intersects(geo_internal::poly_edge(P, i), geo_internal::poly_edge(Q, j))) return true;
    return false;
}
// 距離（多角形は中身を含む。中にある・共有点があれば 0）。boundary_dist は周までの距離（中の点でも周までの距離）
template <class T>
Real boundary_dist(const std::vector<PointT<T>>& poly, const PointT<T>& p) {
    Real r = std::numeric_limits<Real>::infinity();
    for (int i = 0; i < (int)poly.size(); ++i) r = std::min(r, dist(geo_internal::poly_edge(poly, i), p));
    return r;
}
template <class T>
Real dist(const std::vector<PointT<T>>& poly, const PointT<T>& p) { return contains(poly, p) != OUTSIDE ? 0 : boundary_dist(poly, p); }
template <class T>
Real dist(const std::vector<PointT<T>>& poly, const SegmentT<T>& s) {
    if (intersects(poly, s)) return 0;
    Real r = std::numeric_limits<Real>::infinity();
    for (int i = 0; i < (int)poly.size(); ++i) r = std::min(r, dist(geo_internal::poly_edge(poly, i), s));
    return r;
}
template <class T>
Real dist(const std::vector<PointT<T>>& P, const std::vector<PointT<T>>& Q) {
    if (intersects(P, Q)) return 0;
    Real r = std::numeric_limits<Real>::infinity();
    for (int i = 0; i < (int)P.size(); ++i)
        for (int j = 0; j < (int)Q.size(); ++j) r = std::min(r, dist(geo_internal::poly_edge(P, i), geo_internal::poly_edge(Q, j)));
    return r;
}

// ---- 連続する頂点で囲む部分の面積・重心を O(1) で（凸多角形を対角線で切る問題など） ----
// 部分 (i, j) = 頂点 v[i], v[i+1], ..., v[j] を順に結んで v[i] に戻る多角形（番号は n で割った余り。i から j へ番号が増える向き）
// 頂点の数は ((j - i) mod n) + 1。i == j なら 1 点（面積 0）、j = i - 1 なら多角形全体
struct PolygonPrefix {
    int n = 0;
    std::vector<Point> v;
    std::vector<Int128> cr, mx, my, sx, sy;  // 2 周分の累積和：cross(v[t], v[t+1]) / それに (v[t] + v[t+1]) を掛けたもの / 頂点の座標
    PolygonPrefix() {}
    explicit PolygonPrefix(const std::vector<Point>& poly) : n((int)poly.size()), v(poly) {
        cr.assign(2 * n + 1, 0), mx = my = sx = sy = cr;
        for (int t = 0; t < 2 * n; ++t) {
            const Point &p = v[t % n], &q = v[(t + 1) % n];
            const Int128 c = cross128(p, q);
            cr[t + 1] = cr[t] + c, mx[t + 1] = mx[t] + c * (p.x + q.x), my[t + 1] = my[t] + c * (p.y + q.y);
            sx[t + 1] = sx[t] + p.x, sy[t + 1] = sy[t] + p.y;
        }
    }
    // 符号付き面積の 2 倍（多角形が反時計回りなら正）。Int128 で返す（4 倍して比べるなどしても溢れない）
    Int128 area2(int i, int j) const {
        int a, b;
        range(i, j, a, b);
        return cr[b] - cr[a] + cross128(v[b % n], v[a]);
    }
    Real area(int i, int j) const {  // 面積（|area2| / 2）
        const Int128 s = area2(i, j);
        return (Real)(s < 0 ? -s : s) / 2;
    }
    PointD centroid(int i, int j) const {  // 重心（面積の重心。面積 0 なら頂点の平均）
        int a, b;
        range(i, j, a, b);
        const Point &p = v[b % n], &q = v[a];
        const Int128 c = cross128(p, q), A = cr[b] - cr[a] + c;
        if (A == 0) {
            const Real k = (Real)(b - a + 1);
            return PointD((Real)(sx[b + 1] - sx[a]) / k, (Real)(sy[b + 1] - sy[a]) / k);
        }
        const Int128 X = mx[b] - mx[a] + c * (p.x + q.x), Y = my[b] - my[a] + c * (p.y + q.y);
        return PointD((Real)X / ((Real)A * 3), (Real)Y / ((Real)A * 3));
    }

  private:
    void range(int i, int j, int& a, int& b) const {  // 累積和の位置 a ≤ b（a < n）
        a = ((i % n) + n) % n, b = a + (((j - i) % n) + n) % n;
    }
};
#endif  // GG_POLYGON
#ifndef GG_CONVEX_HULL
#define GG_CONVEX_HULL
// ---- 凸包（幾何第二世代）  凸包・上下の凸包・最も遠い 2 点・凸多角形の内外判定（すべて整数で厳密） ----
struct FarthestPair {  // 2 点の番号（i < j。点が 1 つだけなら i == j）と距離の 2 乗。点がなければ -1
    int i = -1, j = -1;
    Coord d2 = -1;
};

// 凸包の頂点の番号（反時計回り、(x, y) が最小の点から）。同じ点は番号の小さい 1 つだけ
// keep_collinear = true なら辺の上の点も含める。全部が一直線なら false → 両端の 2 点、true → 端から端へ全部
inline std::vector<int> convex_hull_idx(const std::vector<Point>& ps, bool keep_collinear = false) {
    const int n = (int)ps.size();
    std::vector<int> id(n);
    for (int i = 0; i < n; ++i) id[i] = i;
    std::sort(id.begin(), id.end(), [&](int a, int b) { return ps[a] != ps[b] ? ps[a] < ps[b] : a < b; });
    std::vector<int> u;  // 同じ点を除いた番号（(x, y) の順）
    for (int i : id)
        if (u.empty() || ps[u.back()] != ps[i]) u.push_back(i);
    const int m = (int)u.size();
    if (m <= 1) return u;
    bool col = true;
    for (int k = 2; k < m && col; ++k) col = cross(ps[u[1]] - ps[u[0]], ps[u[k]] - ps[u[0]]) == 0;
    if (col) return keep_collinear ? u : std::vector<int>{u[0], u[m - 1]};
    auto bad = [&](int o, int a, int b) {  // o → a → b が右に曲がる（keep_collinear = false なら一直線も）
        const Coord c = cross(ps[a] - ps[o], ps[b] - ps[a]);
        return keep_collinear ? c < 0 : c <= 0;
    };
    std::vector<int> h(2 * m);
    int k = 0;
    for (int i = 0; i < m; ++i) {
        while (k >= 2 && bad(h[k - 2], h[k - 1], u[i])) --k;
        h[k++] = u[i];
    }
    for (int i = m - 2, t = k + 1; i >= 0; --i) {
        while (k >= t && bad(h[k - 2], h[k - 1], u[i])) --k;
        h[k++] = u[i];
    }
    h.resize(k - 1);
    return h;
}
inline std::vector<Point> convex_hull(const std::vector<Point>& ps, bool keep_collinear = false) {
    std::vector<Point> r;
    for (int i : convex_hull_idx(ps, keep_collinear)) r.push_back(ps[i]);
    return r;
}
namespace geo_internal {
inline std::vector<Point> half_hull(const std::vector<Point>& ps, bool lower) {
    std::vector<Point> s = ps;
    std::sort(s.begin(), s.end());
    std::vector<Point> u;  // 各 x で、下側なら y が最小、上側なら最大の点だけ
    for (const Point& p : s) {
        if (!u.empty() && u.back().x == p.x) {
            if (!lower) u.back() = p;
            continue;
        }
        u.push_back(p);
    }
    std::vector<Point> h;
    for (const Point& p : u) {
        while (h.size() >= 2) {
            const Coord c = cross(h[h.size() - 1] - h[h.size() - 2], p - h[h.size() - 1]);
            if (lower ? c <= 0 : c >= 0) h.pop_back();
            else break;
        }
        h.push_back(p);
    }
    return h;
}
}  // namespace geo_internal
inline std::vector<Point> lower_hull(const std::vector<Point>& ps) { return geo_internal::half_hull(ps, true); }   // 下側（x の小さい順）
inline std::vector<Point> upper_hull(const std::vector<Point>& ps) { return geo_internal::half_hull(ps, false); }  // 上側（x の小さい順）

// 凸多角形 h（反時計回り・一直線の頂点なし）の最も遠い 2 頂点（番号は h の中）。回転キャリパー
inline FarthestPair convex_diameter(const std::vector<Point>& h) {
    const int n = (int)h.size();
    if (n == 0) return FarthestPair();
    if (n == 1) return FarthestPair{0, 0, 0};
    FarthestPair best{0, 1, dist2(h[0], h[1])};
    if (n == 2) return best;
    auto upd = [&](int a, int b) {
        const Coord d = dist2(h[a], h[b]);
        if (d > best.d2) best = FarthestPair{std::min(a, b), std::max(a, b), d};
    };
    int j = 1;
    for (int i = 0; i < n; ++i) {
        const int ni = (i + 1) % n;
        const Point e = h[ni] - h[i];
        while (cross(e, h[(j + 1) % n] - h[j]) > 0) j = (j + 1) % n;
        upd(i, j), upd(ni, j);
    }
    return best;
}
// 点の集合で最も遠い 2 点（番号は ps の中、i < j）
inline FarthestPair farthest_pair(const std::vector<Point>& ps) {
    const int n = (int)ps.size();
    if (n == 0) return FarthestPair();
    const std::vector<int> idx = convex_hull_idx(ps);
    if (idx.size() == 1) return n == 1 ? FarthestPair{0, 0, 0} : FarthestPair{0, 1, 0};
    std::vector<Point> h;
    for (int i : idx) h.push_back(ps[i]);
    const FarthestPair r = convex_diameter(h);
    const int a = idx[r.i], b = idx[r.j];
    return FarthestPair{std::min(a, b), std::max(a, b), r.d2};
}
// 凸多角形 h（反時計回り。convex_hull の結果）の内外。O(log N)
inline Where convex_contains(const std::vector<Point>& h, const Point& p) {
    const int n = (int)h.size();
    if (n == 0) return OUTSIDE;
    if (n == 1) return h[0] == p ? ON_BORDER : OUTSIDE;
    if (n == 2) return ccw_detail(h[0], h[1], p) == ON_SEGMENT ? ON_BORDER : OUTSIDE;
    const Point a = h[0];
    const Coord s1 = cross(h[1] - a, p - a), sn = cross(h[n - 1] - a, p - a);
    if (s1 < 0 || sn > 0) return OUTSIDE;
    if (s1 == 0) {  // a → h[1] の直線の上：その向きの最も遠い頂点まで（一直線の頂点が続いていてもよい）
        int lo = 1, hi = n;
        while (hi - lo > 1) {
            const int m = (lo + hi) / 2;
            if (cross(h[1] - a, h[m] - a) == 0) lo = m;
            else hi = m;
        }
        return ccw_detail(a, h[lo], p) == ON_SEGMENT ? ON_BORDER : OUTSIDE;
    }
    if (sn == 0) {  // a → h[n-1] の直線の上
        int lo = 0, hi = n - 1;
        while (hi - lo > 1) {
            const int m = (lo + hi) / 2;
            if (cross(h[n - 1] - a, h[m] - a) == 0) hi = m;
            else lo = m;
        }
        return ccw_detail(a, h[hi], p) == ON_SEGMENT ? ON_BORDER : OUTSIDE;
    }
    int lo = 1, hi = n - 1;  // cross(h[lo] - a, p - a) > 0、cross(h[hi] - a, p - a) < 0 の間で二分探索
    while (hi - lo > 1) {
        const int m = (lo + hi) / 2;
        if (cross(h[m] - a, p - a) >= 0) lo = m;
        else hi = m;
    }
    const Coord c = cross(h[lo + 1] - h[lo], p - h[lo]);
    return c > 0 ? INSIDE : (c == 0 ? ON_BORDER : OUTSIDE);
}
// 多角形 Q（凸でなくてよい）が凸多角形 h（周を含む）に入っているか（Q の頂点が全部入っていればよい）。O(M log N)
inline bool convex_contains(const std::vector<Point>& h, const std::vector<Point>& Q) {
    for (const Point& q : Q)
        if (convex_contains(h, q) == OUTSIDE) return false;
    return true;
}
// 凸多角形 h（convex_hull の結果）の内部と周の上の格子点の数（頂点は格子点）。頂点 0〜2 個の退化した形も可
inline Coord lattice_points_in_convex(const std::vector<Point>& h) {
    const int n = (int)h.size();
    if (n <= 1) return n;
    if (n == 2) {
        const Point d = h[1] - h[0];
        return std::gcd(d.x < 0 ? -d.x : d.x, d.y < 0 ? -d.y : d.y) + 1;
    }
    return interior_lattice_points(h) + boundary_lattice_points(h);
}
#endif  // GG_CONVEX_HULL
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
namespace ch_test {

#define CH_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[convex_hull] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK CH_CHECK

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

static bool b_in_tri(const Point& a, const Point& b, const Point& c, const Point& p) {  // 閉じた三角形（一直線なら線分）に含まれるか
    const Int128 s1 = cross128(b - a, p - a), s2 = cross128(c - b, p - b), s3 = cross128(a - c, p - c);
    const bool neg = s1 < 0 || s2 < 0 || s3 < 0, pos = s1 > 0 || s2 > 0 || s3 > 0;
    if (neg && pos) return false;
    if (cross128(b - a, c - a) != 0) return true;
    // 一直線：3 点の作る線分の上か
    Point lo = a, hi = a;
    for (const Point& q : {b, c}) lo = std::min(lo, q), hi = std::max(hi, q);
    return cross128(hi - lo, p - lo) == 0 && lo <= p && p <= hi;
}
// 期待する凸包の頂点の集合（番号は同じ点の中で最小）
static std::vector<int> b_hull_set(const std::vector<Point>& ps, bool keep) {
    const int n = (int)ps.size();
    std::vector<int> d;  // 同じ点の代表
    for (int i = 0; i < n; ++i) {
        bool first = true;
        for (int j = 0; j < i; ++j) first &= ps[j] != ps[i];
        if (first) d.push_back(i);
    }
    std::vector<int> r;
    for (int i : d) {
        bool in;
        if (!keep) {  // ほかの 3 点（2 点でも可）の作る三角形に含まれなければ頂点
            in = false;
            for (int a : d)
                for (int b : d)
                    for (int c : d)
                        if (a != i && b != i && c != i && !in && b_in_tri(ps[a], ps[b], ps[c], ps[i])) in = true;
        } else {  // 支える直線（i を通り、全部の点が左側か直線の上）があれば周の上
            in = true;
            if (d.size() == 1) in = false;
            for (int q : d) {
                if (q == i || !in) continue;
                bool ok = true;
                for (int s : d) ok &= cross128(ps[q] - ps[i], ps[s] - ps[i]) >= 0;
                if (ok) in = false;
            }
        }
        if (!in) r.push_back(i);
    }
    std::sort(r.begin(), r.end());
    return r;
}

static void test_hull(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 6000 * TEST_SCALE; ++it, ++cases) {
        const bool big = it % 4 == 0;
        const int n = (int)rnd(rng, 0, big ? 25 : 11);
        std::vector<Point> ps = rndps(rng, n, big ? 1000000000LL : (it % 4 == 1 ? 1 : 3));
        if (it % 5 == 0 && n > 0) {  // 一直線の点を混ぜる
            const Point d = rndp(rng, 2), o = rndp(rng, big ? 900000000LL : 1);
            for (int k = 0; k < n / 2; ++k) ps[rng() % n] = o + d * rnd(rng, -3, 3);
        }
        for (int keep = 0; keep < 2; ++keep) {
            const std::vector<int> h = convex_hull_idx(ps, keep);
            std::vector<int> s = h;
            std::sort(s.begin(), s.end());
            CHECK(std::unique(s.begin(), s.end()) == s.end());
            // 全部が一直線（または 1 点）か
            std::vector<Point> dp = ps;
            std::sort(dp.begin(), dp.end());
            dp.erase(std::unique(dp.begin(), dp.end()), dp.end());
            bool col = true;
            for (size_t k = 2; k < dp.size(); ++k) col &= cross128(dp[1] - dp[0], dp[k] - dp[0]) == 0;
            if (dp.size() <= 1 || col) {
                std::vector<int> e;  // 端から端へ（(x, y) の順）
                for (const Point& p : dp)
                    for (int i = 0; i < n; ++i)
                        if (ps[i] == p) { e.push_back(i); break; }
                if (!keep && e.size() > 2) e = {e.front(), e.back()};
                CHECK(h == e);
                continue;
            }
            CHECK(s == b_hull_set(ps, keep));
            // 反時計回り・(x, y) 最小の点から・凸
            for (int i : h) CHECK(!(ps[i] < ps[h[0]]));
            const int m = (int)h.size();
            for (int k = 0; k < m; ++k) {
                const Int128 c = cross128(ps[h[(k + 1) % m]] - ps[h[k]], ps[h[(k + 2) % m]] - ps[h[(k + 1) % m]]);
                CHECK(keep ? c >= 0 : c > 0);
            }
            std::vector<Point> hp;
            for (int i : h) hp.push_back(ps[i]);
            CHECK(hp == convex_hull(ps, keep) && is_convex(hp, !keep) && area2(hp) > 0);
            // 内外判定（凸包の中の点・外の点・周の上の点）
            for (int q = 0; q < 20; ++q) {
                Point p = rndp(rng, big ? 1000000000LL : 4);
                if (q % 4 == 0) p = ps[rng() % n];
                CHECK(convex_contains(hp, p) == contains(hp, p));
            }
        }
    }
    // 頂点 0〜2 個の退化した凸包の内外
    CHECK(convex_contains({}, Point(0, 0)) == OUTSIDE && convex_contains({Point(1, 2)}, Point(1, 2)) == ON_BORDER);
    CHECK(convex_contains({Point(0, 0), Point(4, 2)}, Point(2, 1)) == ON_BORDER && convex_contains({Point(0, 0), Point(4, 2)}, Point(6, 3)) == OUTSIDE);
}

static void test_half_far(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 6000 * TEST_SCALE; ++it, ++cases) {
        const bool big = it % 3 == 0;
        const int n = (int)rnd(rng, 1, big ? 40 : 10);
        const std::vector<Point> ps = rndps(rng, n, big ? 1000000000LL : 3);
        // 下側・上側：凸包を (x 最小で y 最小) から (x 最大で y 最小) まで反時計回りにたどったもの
        std::vector<Point> dp = ps;
        std::sort(dp.begin(), dp.end());
        dp.erase(std::unique(dp.begin(), dp.end()), dp.end());
        const Coord xmin = dp.front().x, xmax = dp.back().x;
        const Point ll = dp.front(), ur = dp.back();  // 左下・右上
        Point lr = ur, ul = ll;                        // 右下・左上
        for (const Point& p : dp) {
            if (p.x == xmax && p.y < lr.y) lr = p;
            if (p.x == xmin && p.y > ul.y) ul = p;
        }
        std::vector<Point> lo_e, up_e;
        if (xmin == xmax) lo_e = {ll}, up_e = {ur};
        else {
            const std::vector<Point> h = convex_hull(ps);
            const int m = (int)h.size();
            int k = 0;
            while (h[k] != ll) ++k;
            for (;; k = (k + 1) % m) { lo_e.push_back(h[k]); if (h[k] == lr) break; }
            while (h[k] != ur) k = (k + 1) % m;
            for (;; k = (k + 1) % m) { up_e.push_back(h[k]); if (h[k] == ul) break; }
            std::reverse(up_e.begin(), up_e.end());
        }
        CHECK(lower_hull(ps) == lo_e && upper_hull(ps) == up_e);
        // 最も遠い 2 点
        Coord best = -1;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) best = std::max(best, dist2(ps[i], ps[j]));
        const FarthestPair f = farthest_pair(ps);
        if (n == 1) CHECK(f.i == 0 && f.j == 0 && f.d2 == 0);
        else CHECK(f.d2 == best && f.i < f.j && dist2(ps[f.i], ps[f.j]) == best);
    }
    CHECK(farthest_pair({}).i == -1 && farthest_pair({}).d2 == -1);
    // 凸包の中の格子点（小さい座標で全部数える）
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const std::vector<Point> ps = rndps(rng, (int)rnd(rng, 0, 6), it % 2 ? 2 : 5);
        const std::vector<Point> h = convex_hull(ps);
        Coord c = 0;
        for (Coord x = -5; x <= 5; ++x)
            for (Coord y = -5; y <= 5; ++y) c += convex_contains(h, Point(x, y)) != OUTSIDE;
        CHECK(lattice_points_in_convex(h) == c && lattice_points_in_convex(convex_hull(ps, true)) == c);
    }
}

// 頂点数 n 前後の凸多角形（辺ベクトルの成分は C 以下）
static std::vector<Point> make_convex_polygon(std::mt19937_64& rng, int n, long long C) {
    std::vector<Point> vs;
    for (int k = 0; k < n / 2; ++k) {
        const Point v = rndp(rng, C);
        if (v != Point()) vs.push_back(v), vs.push_back(-v);
    }
    auto h = [](const Point& v) { return (v.y > 0 || (v.y == 0 && v.x > 0)) ? 0 : 1; };
    std::sort(vs.begin(), vs.end(), [&](const Point& a, const Point& b) {
        if (h(a) != h(b)) return h(a) < h(b);
        return cross(a, b) > 0;
    });
    std::vector<Point> P;
    Point cur;
    for (size_t k = 0; k < vs.size(); ++k) {
        if (k > 0 && h(vs[k - 1]) == h(vs[k]) && cross(vs[k - 1], vs[k]) == 0) { cur += vs[k]; continue; }
        P.push_back(cur), cur += vs[k];
    }
    return P;
}

static void test_large_polygon(std::mt19937_64& rng, int& cases) {
    // 頂点の多い凸多角形：凸包が全頂点、内外判定は O(N) の contains と一致、最遠点対は全組と一致
    for (int it = 0; it < 200 * TEST_SCALE; ++it, ++cases) {
        const std::vector<Point> P = make_convex_polygon(rng, (int)rnd(rng, 3, 120), it % 2 ? 5 : 100000);
        if (P.size() < 3) continue;
        std::vector<Point> ps = P;
        std::shuffle(ps.begin(), ps.end(), rng);
        const std::vector<Point> h = convex_hull(ps);
        CHECK(h.size() == P.size() && is_convex(h, true));
        for (int q = 0; q < 50; ++q) {
            Point p = P[rng() % P.size()] + rndp(rng, it % 2 ? 2 : 50000);
            CHECK(convex_contains(h, p) == contains(h, p));
        }
        Coord best = 0;
        for (size_t i = 0; i < ps.size(); ++i)
            for (size_t j = i + 1; j < ps.size(); ++j) best = std::max(best, dist2(ps[i], ps[j]));
        CHECK(farthest_pair(ps).d2 == best && convex_diameter(h).d2 == best);
    }
}

static void test_contains_polygon(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 5000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 4 : 1000000000LL;
        const std::vector<Point> h = convex_hull(rndps(rng, (int)rnd(rng, 1, 12), C), it % 3 == 0);
        std::vector<Point> Q = rndps(rng, (int)rnd(rng, 1, 6), C);
        if (it % 4 == 0) Q = h;
        if (it % 4 == 1 && h.size() >= 3) Q = {h[0], h[1], h[2]};
        CHECK(convex_contains(h, Q) == contains(h, Q));  // 03_polygon の一般の多角形の版（辺が全部入るか）と同じ
    }
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_hull(rng, cases);
    test_contains_polygon(rng, cases);
    test_half_far(rng, cases);
    test_large_polygon(rng, cases);
    std::cout << "[convex_hull] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[convex_hull] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 200000;
    const std::vector<Point> ps = rndps(rng, n, 1000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    const std::vector<Point> h1 = convex_hull(ps);
    const long long t1 = ms_since(t0);
    std::vector<Point> qs = make_convex_polygon(rng, n, 5000);  // 全部が凸包の頂点
    std::shuffle(qs.begin(), qs.end(), rng);
    t0 = std::chrono::steady_clock::now();
    const std::vector<Point> h2 = convex_hull(qs);
    const long long t2 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const FarthestPair f = farthest_pair(qs);
    const long long t3 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    int in = 0;
    for (int k = 0; k < n; ++k) in += convex_contains(h2, qs[k % qs.size()] + rndp(rng, 3)) == INSIDE;
    const long long t4 = ms_since(t0);
    std::cout << "[convex_hull] benchmark: N = " << n << "  convex_hull（ランダム）" << t1 << " ms, convex_hull（全部が頂点 " << h2.size()
              << " 個）" << t2 << " ms, farthest_pair " << t3 << " ms, convex_contains × N " << t4 << " ms  (check "
              << h1.size() + f.d2 % 7 + in % 7 << ")\n";
}

}  // namespace ch_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace ch_demo {

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
#define SHOW(expr, comment) ch_demo::show(#expr, ch_demo::to_s(expr), comment)
#define RUN(...)                   \
    ch_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("凸包 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  convex_hull(ps, keep_collinear) / convex_hull_idx(ps, keep_collinear)   凸包（反時計回り、(x, y) 最小の点から）\n"
              << "  lower_hull(ps) / upper_hull(ps)                                           下側 / 上側の凸包（x の小さい順）\n"
              << "  farthest_pair(ps) / convex_diameter(h)                                    最も遠い 2 点（FarthestPair{i, j, d2}）\n"
              << "  convex_contains(h, p) / convex_contains(h, Q)                             凸多角形の内外 O(log N) / 多角形 Q が入るか\n"
              << "  （同梱の 03_polygon）area2 / area / interior_lattice_points / boundary_lattice_points / contains / is_convex\n";

    api("【1】 convex_hull(ps, keep_collinear = false) / convex_hull_idx(ps, keep_collinear)", "vector<Point> / vector<int>",
        "凸包の頂点（反時計回り、(x, y) が最小の点から）。keep_collinear = true なら辺の上の点も。同じ点は 1 つ", "O(N log N)");
    RUN(std::vector<Point> ps = {{0, 0}, {2, 0}, {4, 0}, {4, 4}, {0, 4}, {1, 1}, {2, 3}, {0, 2}, {4, 4}});
    SHOW(convex_hull(ps), "四角形の 4 頂点");
    SHOW(convex_hull_idx(ps), "番号で（同じ点 (4,4) は番号の小さい 3）");
    SHOW(convex_hull(ps, true), "辺の上の (2,0) と (0,2) も");
    SHOW(convex_hull(std::vector<Point>{{3, 3}, {1, 1}, {2, 2}}), "全部が一直線：両端の 2 点");
    SHOW(convex_hull(std::vector<Point>{{3, 3}, {1, 1}, {2, 2}}, true), "keep_collinear なら端から端へ全部");

    api("【2】 凸包と 03_polygon（同梱）", "Coord / Real",
        "凸包の面積・周の長さ・格子点の数は、凸包を多角形として 03_polygon の関数に渡す", "O(N)");
    RUN(auto h = convex_hull(ps));
    SHOW(area2(h), "面積の 2 倍（4 × 4 × 2）");
    SHOW(interior_lattice_points(h) + boundary_lattice_points(h), "凸包の中（周を含む）の格子点の数 = 25");
    SHOW(lattice_points_in_convex(h), "同じもの（凸包が点・線分の退化した形でも正しく数える）");
    note("典型90 041（新しく打てる杭の数）: lattice_points_in_convex(convex_hull(ps)) - N");

    api("【3】 farthest_pair(ps) / convex_diameter(h)", "FarthestPair{i, j, d2}",
        "最も遠い 2 点の番号（ps の中、i < j）と距離の 2 乗。凸包をとってから回転キャリパー", "O(N log N) / O(N)");
    RUN(FarthestPair f = farthest_pair(ps));
    std::cout << "  f.i, f.j, f.d2                     -> " << f.i << ", " << f.j << ", " << f.d2 << "        // (0,0) と (4,4)：距離² = 32\n";
    note("距離そのものは std::sqrt((Real)f.d2)。点が 1 つなら i == j、なければ -1");

    api("【4】 convex_contains(h, p) / convex_contains(h, Q)", "Where（OUTSIDE / ON_BORDER / INSIDE）/ bool",
        "凸多角形 h（convex_hull の結果）の内外を O(log N) で。クエリが多いとき（ABC296 G）", "O(log N)");
    SHOW(convex_contains(h, Point(1, 3)) == INSIDE, "内側");
    SHOW(convex_contains(h, Point(4, 1)) == ON_BORDER, "辺の上");
    SHOW(convex_contains(h, Point(5, 2)) == OUTSIDE, "外");
    SHOW(convex_contains(h, std::vector<Point>{Point(1, 1), Point(3, 1), Point(2, 3)}), "三角形が全部入る（O(M log N)）");
    note("凸でない多角形なら 03_polygon の contains(poly, p)（O(N)）");
    note("h[0] は角（180 度でない頂点）であること。入力の凸多角形をそのまま使うなら h = convex_hull(P, true) にしてから");

    api("【5】 lower_hull(ps) / upper_hull(ps)", "vector<Point>",
        "下側・上側の凸包（x の小さい順。各 x で一番下 / 一番上の点だけ、一直線の点なし）", "O(N log N)");
    SHOW(lower_hull(ps), "下側");
    SHOW(upper_hull(ps), "上側");
    note("点 (x_i, y_i) の a·x_i + y_i の最小は、lower_hull の上で a に応じた頂点（二分探索・尺取り）。クエリが多いなら 04_convex_polygon の extreme_point（O(log N)）");
}

}  // namespace ch_demo


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
    ch_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ch_test::benchmark();    // 速度計測
    ch_demo::demo();         // 使い方デモ
    return 0;
}

