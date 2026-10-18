// =============================================================================
//  マンハッタン距離の最小全域木（候補の辺 4N 本・最小全域木・チェビシェフ距離版）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 07 マンハッタン距離・チェビシェフ距離
//  依存       : 00_geo_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - 点の番号は入力順（0-indexed）。辺は ManhattanEdges{us, vs, ws}（us[k] - vs[k]、重み ws[k]）
//      graph_grid_v2 の graph_from_edges(N, e.us, e.vs, e.ws, UNDIRECTED, 0) にそのまま渡せる
//    - 座標は 10^9 まで（同じ点があってもよい）
//
//  機能一覧                                                                    計算量
//    manhattan_mst_candidates(ps)      最小全域木に使われうる辺（4N 本以下）。全部の組の代わりにグラフに使う   O(N log N)
//    manhattan_mst(ps)                 マンハッタン距離の最小全域木 ManhattanMST{cost, edges}（edges は重みの小さい順）  O(N log N)
//    chebyshev_mst(ps)                 チェビシェフ距離 max(|dx|, |dy|) の最小全域木（45 度回して同じ方法で）   O(N log N)
//    manhattan_nearest_each(ps)        各点からマンハッタン距離で最も近いほかの点 ManhattanNear{j, d}（番号・距離。点が 1 つなら -1）  O(N log N)
//                                      ※ 候補の辺（右半分の 4 方向）だけでは各点の最近点は求まらない。これは 8 方向を見る
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mm_test::random_test()  … 全部の組を使うプリム法と突き合わせるランダムテスト
//    mm_test::benchmark()    … 速度計測
//    mm_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_MANHATTAN_MST
#define GG_MANHATTAN_MST
// ---- マンハッタン距離の最小全域木（幾何第二世代）  候補の辺 4N 本から ----
struct ManhattanEdges {  // 辺 k：us[k] - vs[k]、重み ws[k]（graph_from_edges(N, us, vs, ws, UNDIRECTED, 0) に渡せる）
    std::vector<int> us, vs;
    std::vector<Coord> ws;
    int size() const { return (int)us.size(); }
    void add(int u, int v, Coord w) { us.push_back(u), vs.push_back(v), ws.push_back(w); }
};
struct ManhattanMST {
    Coord cost = 0;        // 重みの合計
    ManhattanEdges edges;  // 選んだ辺（N - 1 本。重みの小さい順）
};

// 最小全域木に使われうる辺（各点から、45 度ずつに分けた右半分の 4 方向それぞれで最も近い点への辺。4N 本以下）
// 方向 [45°, 90°]（dx ≥ 0, dy ≥ dx）で最も近い点 = x + y が最小の点。x の大きい順に見て、y - x が自分以上の点の中の最小を BIT で
// ほかの 3 方向は座標を入れ替え・反転して同じことをする（並べ替えは 4 回で済むよう共有する）
inline ManhattanEdges manhattan_mst_candidates(const std::vector<Point>& ps) {
    const int n = (int)ps.size();
    ManhattanEdges E;
    E.us.reserve(4 * n), E.vs.reserve(4 * n), E.ws.reserve(4 * n);
    // 並べ方：ox = (x, y) の小さい順、oy = (y, x) の小さい順
    std::vector<std::pair<std::pair<Coord, Coord>, int>> buf(n);
    std::vector<int> ox(n), oy(n);
    for (int i = 0; i < n; ++i) buf[i] = {{ps[i].x, ps[i].y}, i};
    std::sort(buf.begin(), buf.end());
    for (int i = 0; i < n; ++i) ox[i] = buf[i].second;
    for (int i = 0; i < n; ++i) buf[i] = {{ps[i].y, ps[i].x}, i};
    std::sort(buf.begin(), buf.end());
    for (int i = 0; i < n; ++i) oy[i] = buf[i].second;
    // 順位：rd = y - x の小さい順、rs = x + y の小さい順（同じ値は同じ順位）。md, ms は値の種類の数
    std::vector<std::pair<Coord, int>> kv(n);
    std::vector<int> rd(n), rs(n);
    auto ranks = [&](std::vector<int>& rk) {
        std::sort(kv.begin(), kv.end());
        int r = 0;
        for (int k = 0; k < n; ++k) {
            if (k > 0 && kv[k].first != kv[k - 1].first) ++r;
            rk[kv[k].second] = r;
        }
        return n == 0 ? 0 : r + 1;
    };
    for (int i = 0; i < n; ++i) kv[i] = {ps[i].y - ps[i].x, i};
    const int md = ranks(rd);
    for (int i = 0; i < n; ++i) kv[i] = {ps[i].x + ps[i].y, i};
    const int ms = ranks(rs);
    auto rev_groups = [&](std::vector<int>& o, bool byx) {  // 同じ x（y）の中だけ逆順にする
        for (int k = 0; k < n;) {
            int l = k;
            while (l < n && (byx ? ps[o[l]].x == ps[o[k]].x : ps[o[l]].y == ps[o[k]].y)) ++l;
            std::reverse(o.begin() + k, o.begin() + l);
            k = l;
        }
    };
    std::vector<int> bi(n + 1), ord(n), key(n);
    std::vector<Coord> bv(n + 1), val(n);
    const Coord INF_ = std::numeric_limits<Coord>::max();
    for (int dir = 0; dir < 4; ++dir) {
        // dir 0：(x, y)、1：(y, x)、2：(-y, x)、3：(x, -y) と見て、方向 [45°, 90°] を探す
        if (dir == 0) ord.assign(ox.rbegin(), ox.rend());               // x 大 → y 大
        if (dir == 1) ord.assign(oy.rbegin(), oy.rend());               // y 大 → x 大
        if (dir == 2) ord = oy, rev_groups(ord, false);                 // y 小 → x 大
        if (dir == 3) ord.assign(ox.rbegin(), ox.rend()), rev_groups(ord, true);  // x 大 → y 小
        for (int i = 0; i < n; ++i) {
            const Coord x = ps[i].x, y = ps[i].y;
            if (dir == 0) key[i] = md - rd[i], val[i] = x + y;  // 変換後の y - x が大きいほど小さい番号（1 から）
            if (dir == 1) key[i] = rd[i] + 1, val[i] = x + y;
            if (dir == 2) key[i] = ms - rs[i], val[i] = x - y;
            if (dir == 3) key[i] = rs[i] + 1, val[i] = x - y;
        }
        std::fill(bv.begin(), bv.end(), INF_), std::fill(bi.begin(), bi.end(), -1);
        for (const int i : ord) {
            Coord best = INF_;
            int bj = -1;
            for (int p = key[i]; p > 0; p -= p & -p)
                if (bv[p] < best) best = bv[p], bj = bi[p];
            if (bj >= 0) E.add(i, bj, best - val[i]);
            for (int p = key[i]; p <= n; p += p & -p)
                if (val[i] < bv[p]) bv[p] = val[i], bi[p] = i;
        }
    }
    return E;
}

namespace geo_internal {
inline ManhattanMST mst_from_candidates(int n, const ManhattanEdges& E) {
    std::vector<std::pair<Coord, int>> ord(E.size());
    for (int k = 0; k < E.size(); ++k) ord[k] = {E.ws[k], k};
    std::sort(ord.begin(), ord.end());
    std::vector<int> par(n);
    for (int i = 0; i < n; ++i) par[i] = i;
    auto find = [&](int x) {
        while (par[x] != x) x = par[x] = par[par[x]];
        return x;
    };
    ManhattanMST r;
    for (const auto& wk : ord) {
        const int k = wk.second, a = find(E.us[k]), b = find(E.vs[k]);
        if (a == b) continue;
        par[a] = b;
        r.cost += E.ws[k], r.edges.add(E.us[k], E.vs[k], E.ws[k]);
        if (r.edges.size() == n - 1) break;
    }
    return r;
}
}  // namespace geo_internal

// マンハッタン距離の最小全域木（edges は重みの小さい順、N - 1 本）
inline ManhattanMST manhattan_mst(const std::vector<Point>& ps) {
    return geo_internal::mst_from_candidates((int)ps.size(), manhattan_mst_candidates(ps));
}
// チェビシェフ距離の最小全域木（45 度回す：チェビシェフ距離 = 回した点のマンハッタン距離 / 2）
inline ManhattanMST chebyshev_mst(const std::vector<Point>& ps) {
    std::vector<Point> q;
    q.reserve(ps.size());
    for (const Point& p : ps) q.push_back(Point(p.x + p.y, p.x - p.y));
    ManhattanMST r = manhattan_mst(q);
    r.cost /= 2;
    for (Coord& w : r.edges.ws) w /= 2;
    return r;
}
struct ManhattanNear {  // 最も近いほかの点の番号 j と距離 d（点が 1 つなら j = -1, d = -1）
    int j = -1;
    Coord d = -1;
};
// 各点からマンハッタン距離で最も近いほかの点（同じ距離ならどれか）。O(N log N)
// 候補の辺は右半分の 4 方向しか見ないので、点を 180 度回した集合でも同じことをして 8 方向にする
inline std::vector<ManhattanNear> manhattan_nearest_each(const std::vector<Point>& ps) {
    const int n = (int)ps.size();
    std::vector<ManhattanNear> best(n);
    std::vector<Point> rev(n);
    for (int i = 0; i < n; ++i) rev[i] = -ps[i];
    for (int t = 0; t < 2; ++t) {
        const ManhattanEdges E = manhattan_mst_candidates(t == 0 ? ps : rev);
        for (int k = 0; k < E.size(); ++k) {
            ManhattanNear& b = best[E.us[k]];
            if (b.j < 0 || E.ws[k] < b.d) b.j = E.vs[k], b.d = E.ws[k];
        }
    }
    std::vector<int> id(n);  // 同じ点があれば距離 0 で、それが最も近い
    for (int i = 0; i < n; ++i) id[i] = i;
    std::sort(id.begin(), id.end(), [&](int a, int b) { return ps[a] < ps[b]; });
    for (int k = 0; k + 1 < n; ++k)
        if (ps[id[k]] == ps[id[k + 1]]) best[id[k]] = {id[k + 1], 0}, best[id[k + 1]] = {id[k], 0};
    return best;
}
#endif  // GG_MANHATTAN_MST
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
namespace mm_test {

#define MM_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[manhattan_mst] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK MM_CHECK

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

static Coord mdist(const Point& a, const Point& b) { return std::llabs(a.x - b.x) + std::llabs(a.y - b.y); }
static Coord cdist(const Point& a, const Point& b) { return std::max(std::llabs(a.x - b.x), std::llabs(a.y - b.y)); }
static Coord prim(const std::vector<Point>& ps, bool cheb) {
    const int n = (int)ps.size();
    if (n == 0) return 0;
    auto w = [&](int i, int j) {
        const Coord dx = std::llabs(ps[i].x - ps[j].x), dy = std::llabs(ps[i].y - ps[j].y);
        return cheb ? std::max(dx, dy) : dx + dy;
    };
    std::vector<Coord> best(n, std::numeric_limits<Coord>::max());
    std::vector<char> used(n, 0);
    best[0] = 0;
    Coord total = 0;
    for (int it = 0; it < n; ++it) {
        int v = -1;
        for (int i = 0; i < n; ++i)
            if (!used[i] && (v < 0 || best[i] < best[v])) v = i;
        used[v] = 1, total += best[v];
        for (int i = 0; i < n; ++i)
            if (!used[i]) best[i] = std::min(best[i], w(v, i));
    }
    return total;
}
static bool is_spanning_tree(int n, const ManhattanEdges& e) {
    if (e.size() != std::max(0, n - 1)) return false;
    std::vector<int> par(n);
    for (int i = 0; i < n; ++i) par[i] = i;
    std::function<int(int)> f = [&](int x) { return par[x] == x ? x : par[x] = f(par[x]); };
    for (int k = 0; k < e.size(); ++k) {
        const int a = f(e.us[k]), b = f(e.vs[k]);
        if (a == b) return false;
        par[a] = b;
    }
    return true;
}

static void test_mst(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 5000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, it % 5 == 0 ? 60 : 12);
        const long long C = it % 3 == 0 ? 1000000000LL : (it % 3 == 1 ? 2 : 30);
        std::vector<Point> ps = rndps(rng, n, C);
        if (it % 7 == 0 && n > 2) ps[1] = ps[0], ps[2] = ps[0];  // 同じ点
        const ManhattanEdges cand = manhattan_mst_candidates(ps);
        CHECK(cand.size() <= 4 * n);
        for (int k = 0; k < cand.size(); ++k) CHECK(cand.us[k] != cand.vs[k] && cand.ws[k] == mdist(ps[cand.us[k]], ps[cand.vs[k]]));
        for (int cheb = 0; cheb < 2; ++cheb) {
            const ManhattanMST m = cheb ? chebyshev_mst(ps) : manhattan_mst(ps);
            CHECK(m.cost == prim(ps, cheb) && is_spanning_tree(n, m.edges));
            Coord s = 0;
            for (int k = 0; k < m.edges.size(); ++k) {
                const Point a = ps[m.edges.us[k]], b = ps[m.edges.vs[k]];
                CHECK(m.edges.ws[k] == (cheb ? cdist(a, b) : mdist(a, b)));
                CHECK(k == 0 || m.edges.ws[k - 1] <= m.edges.ws[k]);
                s += m.edges.ws[k];
            }
            CHECK(s == m.cost);
        }
    }
}

static void test_nearest(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 5000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, it % 5 == 0 ? 60 : 12);
        const long long C = it % 3 == 0 ? 1000000000LL : (it % 3 == 1 ? 2 : 30);
        std::vector<Point> ps = rndps(rng, n, C);
        if (it % 7 == 0 && n > 2) ps[n - 1] = ps[0];  // 同じ点
        const std::vector<ManhattanNear> nb = manhattan_nearest_each(ps);
        CHECK((int)nb.size() == n);
        for (int i = 0; i < n; ++i) {
            Coord e = std::numeric_limits<Coord>::max();
            for (int j = 0; j < n; ++j)
                if (j != i) e = std::min(e, mdist(ps[i], ps[j]));
            if (n == 1) CHECK(nb[i].j == -1 && nb[i].d == -1);
            else CHECK(nb[i].j >= 0 && nb[i].j != i && mdist(ps[i], ps[nb[i].j]) == e && nb[i].d == e);
        }
    }
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_mst(rng, cases);
    test_nearest(rng, cases);
    std::cout << "[manhattan_mst] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[manhattan_mst] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 200000;
    const std::vector<Point> ps = rndps(rng, n, 1000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    const ManhattanMST m = manhattan_mst(ps);
    const long long t1 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const ManhattanMST c = chebyshev_mst(ps);
    const long long t2 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<ManhattanNear> nb = manhattan_nearest_each(ps);
    const long long t3 = ms_since(t0);
    std::cout << "[manhattan_mst] benchmark: N = " << n << "  manhattan_mst " << t1 << " ms, chebyshev_mst " << t2 << " ms, manhattan_nearest_each " << t3
              << " ms  (check " << (m.cost + c.cost + nb[0].d) % 7 << ")\n";
}

}  // namespace mm_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace mm_demo {

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
#define SHOW(expr, comment) mm_demo::show(#expr, mm_demo::to_s(expr), comment)
#define RUN(...)                   \
    mm_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("マンハッタン距離の最小全域木 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  manhattan_mst_candidates(ps)      候補の辺（4N 本以下）ManhattanEdges{us, vs, ws}\n"
              << "  manhattan_mst(ps)                 最小全域木 ManhattanMST{cost, edges}\n"
              << "  chebyshev_mst(ps)                 チェビシェフ距離の最小全域木\n"
              << "  manhattan_nearest_each(ps)        各点から最も近いほかの点 ManhattanNear{j, d}\n";

    RUN(std::vector<Point> ps = {{0, 0}, {3, 1}, {1, 4}, {5, 5}, {-2, 3}});

    api("【1】 manhattan_mst(ps)", "ManhattanMST{cost, edges}", "全部の点をマンハッタン距離の辺でつなぐ最小のコストと、選んだ辺（重みの小さい順）", "O(N log N)");
    RUN(ManhattanMST m = manhattan_mst(ps));
    SHOW(m.cost, "重みの合計");
    for (int k = 0; k < m.edges.size(); ++k)
        std::cout << "  m.edges の " << k << " 本目: " << m.edges.us[k] << " - " << m.edges.vs[k] << "（重み " << m.edges.ws[k] << "）\n";
    note("全部の組（N(N-1)/2 本）を作ると N = 2×10^5 では間に合わない → 候補 4N 本だけで最小全域木になる");

    api("【2】 manhattan_mst_candidates(ps)", "ManhattanEdges{us, vs, ws}", "最小全域木に使われうる辺（45 度ずつの 8 方向それぞれで最も近い点への辺。半分の 4 方向だけ見れば足りる）。4N 本以下", "O(N log N)");
    RUN(ManhattanEdges e = manhattan_mst_candidates(ps));
    SHOW(e.size(), "候補の数");
    note("graph_grid_v2 と一緒に：Graph G = graph_from_edges(N, e.us, e.vs, e.ws, UNDIRECTED, 0); auto r = kruskal(G);");
    note("はじめからつながっている点があるとき：つながっている点を unite した UnionFind uf で kruskal(G, uf) とすればよい（使えない辺がある問題には使えない）");

    api("【3】 chebyshev_mst(ps)", "ManhattanMST", "チェビシェフ距離 max(|dx|, |dy|) の最小全域木（(x + y, x - y) に回してマンハッタン距離で）", "O(N log N)");
    SHOW(chebyshev_mst(ps).cost, "");

    api("【4】 manhattan_nearest_each(ps)", "vector<ManhattanNear{j, d}>", "各点からマンハッタン距離で最も近いほかの点の番号 j と距離 d（同じ距離ならどれか。点が 1 つなら -1）", "O(N log N)");
    RUN(std::vector<ManhattanNear> nb = manhattan_nearest_each(ps));
    SHOW(nb[0].j, "0 から最も近いのは 1");
    SHOW(nb[0].d, "距離 4");
    note("ABC283 F（各 i で min |P_i - P_j| + |i - j|）: 点 (i, P_i) にして manhattan_nearest_each(ps)[i].d");
    note("候補の辺 manhattan_mst_candidates は右半分の 4 方向だけなので、各点の最近点を求めるのには使えない（これを使う）");
}

}  // namespace mm_demo


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
    mm_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mm_test::benchmark();    // 速度計測
    mm_demo::demo();         // 使い方デモ
    return 0;
}

