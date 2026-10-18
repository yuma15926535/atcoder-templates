// =============================================================================
//  幾何基盤（点 Point / PointD・直線・線分・円の型、向き、回転、角度、入力）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 00 基盤（ほかの幾何テンプレートはすべてこれを同梱している）
//  依存       : なし
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  共通の約束（第二世代の幾何テンプレートすべてで同じ）
//    - 点は 2 種類だけ：整数の Point（Coord = long long）と実数の PointD（Real = long double）
//    - 判定（向き・交差・内外・円の位置関係・偏角）と、面積の 2 倍・距離の 2 乗は整数で厳密に求める（誤差なし）
//    - 整数の座標は |x|, |y| ≤ 10^9（半径も 10^9）まで正しい。はみ出す途中の計算は Int128（= __int128）
//    - 答えが整数にならないもの（交点・回転・角度など）だけ PointD / Real で返す。実数の比較は相対誤差 GEO_EPS
//      実数の点の誤差は「座標の大きさ × 1e-19 程度」（long double。Real を double にすると 1e-16 程度）
//    - 反時計回りが正（cross(a, b) > 0 は b が a の左）。偏角は x 軸の正の向きから反時計回りに [0, 2π)
//    - 返す多角形は反時計回り。交点の「なし」は空の vector、番号の「なし」は -1
//    - 内外判定は OUTSIDE / ON_BORDER / INSIDE（Windows のマクロ IN・OUT とぶつからない名前）
//    - 座標の向き：幾何の (x, y) は「横・縦（上が正）」。グリッドの (i, j)（行・列。下が正）とは
//      回転の向きが逆になる。グリッドのマスを幾何で扱うなら Point{j, -i} にすると見た目と向きがそろう
//
//  機能一覧                                                                    計算量
//    Coord / Real / Int128 / GEO_EPS / GEO_PI   座標の型 / 実数の型 / 128 ビット整数 / 相対誤差 1e-9 / π
//    Point p(x, y) / PointD q(x, y)    整数の点 / 実数の点（PointT<T> の別名）。p.x, p.y
//      + - 単項- *k k* /k += -= *= /= == != < > <= >=（x → y の辞書順）
//      整数の点に実数を掛ける・割るとコンパイルエラー（黙って切り捨てないため）。実数にするには to_d(p)
//    dot(a, b) / cross(a, b) / norm2(a)           内積 / 外積 / 長さの 2 乗（整数のまま）     O(1)
//    dot(o, a, b) / cross(o, a, b)                o を始点にした内積・外積（cross(o, a, b) = cross(a - o, b - o)）
//    dot128 / cross128 / norm2_128                Int128 で計算する版（大きい座標の積の比較に）
//    len(a) / dist2(a, b) / dist(a, b)            長さ（Real）/ 距離の 2 乗（整数）/ 距離（Real）
//    ccw(a, b, c)                                 +1 左に曲がる / -1 右 / 0 一直線（整数なら厳密）
//    ccw_detail(a, b, c)                          一直線のときの位置まで（AOJ CGL_1_C と同じ 5 分類）
//        COUNTER_CLOCKWISE(1) / CLOCKWISE(-1) / ONLINE_BACK(2: c-a-b) / ONLINE_FRONT(-2: a-b-c) / ON_SEGMENT(0: a-c-b)
//    collinear(a, b, c)                           3 点が一直線か
//    rot90(p) / rot90_cw(p) / rot180(p)           反時計回り 90 度 / 時計回り 90 度 / 180 度（整数のまま）
//    rot90(p, c, k)                               c のまわりに反時計回りに 90 度を k 回（負なら時計回り）
//    rotated(p, θ [, c]) / rotated_deg(p, 度 [, c])   任意の角だけ反時計回りに回す（PointD。90 度の倍数は厳密。c は PointD でもよい）
//    complex_mul(a, b) / conjugate(a)             複素数としての積 / 共役（回転・拡大を整数のまま）
//    Affine f / AffineD                           アフィン変換 (a x + b y + c, d x + e y + f)。f(p) で適用、f.then(g) で「f のあと g」
//      Affine::rot90(o) / rot90_cw(o) / rot180(o) / reflect_x(p)（直線 x = p で折り返し）/ reflect_y(p) / translate(v) / scale(k, o)
//      AffineD::rotate(θ, o)                      操作の列を前から合成しておけば「i 回目の操作のあとの位置」が O(1)（ABC189 E）
//    to_d(p / l / s / c)                          実数の型に直す
//    midpoint(a, b) / unit(v) / from_polar(r, θ)  中点 / 長さ 1 のベクトル / 極座標から（PointD）
//    angle_of(v) / angle_of_deg(v)                偏角 [0, 2π) / [0, 360)（0 ベクトルは 0）
//    angle_between(a, b) / angle_between_deg      2 つのベクトルのなす角 [0, π] / [0, 180]
//    deg_to_rad(d) / rad_to_deg(r)                度 ↔ ラジアン
//    real_sign(x [, scale]) / real_cmp(a, b) / real_eq(a, b)   実数の符号・比較（相対誤差。点どうしも可）
//    direction_key(v)                             gcd で割った向き（(2,4) → (1,2)、(-2,-4) → (-1,-2)）
//    line_direction_key(v)                        向きを区別しない（v と -v を同じにする。傾きのキー）
//    isqrt(n) / isqrt_ceil(n)                     ⌊√n⌋ / ⌈√n⌉（誤差なし。n ≤ 9.2×10^18）     O(1)
//    parse_fixed(s, k)                            小数の文字列を 10^k 倍した整数に（"1.25", 2 → 125。誤差なし）
//    read_fixed(k) / read_point_fixed(k) / read_points_fixed(n, k)   小数で与えられる値・点を読む
//    read_point() / read_points(n) / read_point_d() / read_points_d(n)   点を読む（x y）
//    read_segment() / read_line() / read_circle()                         x1 y1 x2 y2 / x1 y1 x2 y2 / x y r
//    int128_to_string(v)                          Int128 を文字列に（cout で直接出せないため）
//    PointHash                                    unordered_set<Point, PointHash> などに使うハッシュ
//    Line l(a, b) / Segment s(a, b) / Circle c(中心, r)   直線 a→b（向きあり）/ 線分 / 円（LineD / SegmentD / CircleD も）
//    Where（OUTSIDE / ON_BORDER / INSIDE）        内外判定の結果（多角形・円・凸多角形のテンプレートで使う）
//    cout << p                                    "(x, y)"（デバッグ表示。Main_base の debug(p) でも使える）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    geob_test::random_test()  … 128 ビットの愚直計算・別の式と全機能を突き合わせるランダムテスト
//    geob_test::benchmark()    … 速度計測
//    geob_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
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
namespace geob_test {

#define GEOB_CHECK(...)                                                                         \
    do {                                                                                        \
        if (!(__VA_ARGS__)) {                                                                   \
            std::cerr << "[geo_base] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                       \
        }                                                                                       \
    } while (0)

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
static const int TEST_SCALE = 1;
#else
static const int TEST_SCALE = 8;
#endif

// 実数の許容誤差（Real の精度に合わせる。long double なら 1e-15、double なら 1e-12 程度）
static const Real EPSR = std::numeric_limits<Real>::epsilon();
static const Real TOL = std::max<Real>(EPSR * 4096, (Real)1e-15L);

static long long rnd(std::mt19937_64& rng, long long lo, long long hi) { return lo + (long long)(rng() % (unsigned long long)(hi - lo + 1)); }
static Real rndr(std::mt19937_64& rng, Real lo, Real hi) { return lo + (hi - lo) * (Real)(rng() >> 11) / (Real)(1ULL << 53); }
static Point rndp(std::mt19937_64& rng, long long c) { return Point(rnd(rng, -c, c), rnd(rng, -c, c)); }
static int sgn128(Int128 v) { return (v > 0) - (v < 0); }
static bool near(Real a, Real b, Real tol) { return std::fabs(a - b) <= tol * std::max<Real>(1, std::max(std::fabs(a), std::fabs(b))); }
static bool nearp(const PointD& a, const PointD& b, Real tol) {
    const Real s = std::max<Real>(1, std::max(std::max(std::fabs(a.x), std::fabs(a.y)), std::max(std::fabs(b.x), std::fabs(b.y))));
    return std::fabs(a.x - b.x) <= tol * s && std::fabs(a.y - b.y) <= tol * s;
}
// 愚直な 5 分類（128 ビットで計算）
static int brute_ccw_detail(const Point& a, const Point& b, const Point& c) {
    const Point u = b - a, v = c - a;
    const int s = sgn128(cross128(u, v));
    if (s != 0) return s;
    if (dot128(u, v) < 0) return ONLINE_BACK;
    if (norm2_128(u) < norm2_128(v)) return ONLINE_FRONT;
    return ON_SEGMENT;
}

static void test_point_ops(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 1000000000LL : 10;
        const Point a = rndp(rng, C), b = rndp(rng, C);
        const long long k = rnd(rng, -3, 3);
        GEOB_CHECK((a + b).x == a.x + b.x && (a + b).y == a.y + b.y);
        GEOB_CHECK((a - b).x == a.x - b.x && (a - b).y == a.y - b.y);
        GEOB_CHECK((-a).x == -a.x && (-a).y == -a.y && +a == a);
        GEOB_CHECK(a * k == Point(a.x * k, a.y * k) && k * a == a * k);
        if (k != 0) GEOB_CHECK((a * k) / k == a);
        Point c = a;
        c += b, c -= b, c *= 2, c /= 2;
        GEOB_CHECK(c == a);
        GEOB_CHECK((a == b) == (a.x == b.x && a.y == b.y) && (a != b) == !(a == b));
        GEOB_CHECK((a < b) == (std::make_pair(a.x, a.y) < std::make_pair(b.x, b.y)));
        GEOB_CHECK((a > b) == (b < a) && (a <= b) == !(b < a) && (a >= b) == !(a < b));
        // 内積・外積・距離（|座標| ≤ 10^9 なら差の積でも long long に収まる）
        const Point o = rndp(rng, C);
        GEOB_CHECK((Int128)dot(a - o, b - o) == dot128(a - o, b - o) && dot(o, a, b) == dot(a - o, b - o));
        GEOB_CHECK((Int128)cross(a - o, b - o) == cross128(a - o, b - o) && cross(o, a, b) == cross(a - o, b - o));
        GEOB_CHECK((Int128)norm2(a - b) == norm2_128(a - b) && dist2(a, b) == norm2(a - b));
        GEOB_CHECK(near(dist(a, b), std::sqrt((Real)norm2_128(a - b)), TOL));
        GEOB_CHECK(near(len(a), std::hypot((Real)a.x, (Real)a.y), TOL));
        // 実数の点
        const PointD da = to_d(a), db = to_d(b);
        GEOB_CHECK(da.x == (Real)a.x && da.y == (Real)a.y && PointD(a) == da);
        const Real lab = std::max<Real>(1, len(a) * len(b));  // 内積・外積の誤差は |a||b| に比例する
        GEOB_CHECK(std::fabs(dot(da, db) - (Real)dot128(a, b)) <= TOL * lab && std::fabs(cross(da, db) - (Real)cross128(a, b)) <= TOL * lab);
        const PointD m = midpoint(a, b);
        GEOB_CHECK(m.x * 2 == (Real)a.x + (Real)b.x && m.y * 2 == (Real)a.y + (Real)b.y);
        const PointD half = da / (Real)2;
        GEOB_CHECK(half.x * 2 == da.x && (da * (Real)0.5) == half && ((Real)0.5 * da) == half);
        // 複素数としての積
        const Point p = rndp(rng, 30000), q = rndp(rng, 30000);
        const std::complex<long long> cp(p.x, p.y), cq(q.x, q.y), cr = cp * cq;
        GEOB_CHECK(complex_mul(p, q) == Point(cr.real(), cr.imag()) && conjugate(p) == Point(p.x, -p.y));
        GEOB_CHECK(norm2(complex_mul(p, q)) == norm2(p) * norm2(q));
    }
    // 型と変換
    const Point p(3, -4);
    const PointD q = to_d(p);
    GEOB_CHECK(q.x == 3 && q.y == -4 && len(p) == 5 && norm2(p) == 25);
    constexpr Point cp(1, 2);
    static_assert(cp.x == 1 && rot90(cp) == Point(-2, 1), "constexpr");
    auto [sx, sy] = p;
    GEOB_CHECK(sx == 3 && sy == -4);
    static_assert(!std::is_constructible_v<Point, double, double>, "整数の点を実数から作れない");
    static_assert(std::is_constructible_v<PointD, long long, int>, "実数の点は整数から作れる");
    static_assert(!std::is_convertible_v<Point, PointD>, "Point → PointD は明示（to_d）");
}

static void test_ccw(std::mt19937_64& rng, int& cases) {
    // 小さい座標：全部の退化ケースを含む
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 2 : 6;
        const Point a = rndp(rng, C), b = rndp(rng, C), c = rndp(rng, C);
        GEOB_CHECK(ccw(a, b, c) == sgn128(cross128(b - a, c - a)));
        if (a != b) GEOB_CHECK(ccw_detail(a, b, c) == brute_ccw_detail(a, b, c));
        GEOB_CHECK(collinear(a, b, c) == (cross128(b - a, c - a) == 0));
        GEOB_CHECK(geo_internal::sgn_cross(b - a, c - a) == sgn128(cross128(b - a, c - a)) && geo_internal::sgn_dot(b - a, c - a) == sgn128(dot128(b - a, c - a)));
        GEOB_CHECK(geo_internal::sgn_cross(to_d(b - a), to_d(c - a)) == sgn128(cross128(b - a, c - a)) && geo_internal::sgn_dot(to_d(b - a), to_d(c - a)) == sgn128(dot128(b - a, c - a)));
        GEOB_CHECK(ccw(a, b, c) == -ccw(a, c, b) && ccw(a, b, c) == ccw(b, c, a));
    }
    GEOB_CHECK(ccw_detail(Point(1, 1), Point(1, 1), Point(1, 1)) == ON_SEGMENT);
    GEOB_CHECK(ccw_detail(Point(1, 1), Point(1, 1), Point(2, 1)) == ONLINE_FRONT);
    // 大きい座標：一直線の 3 点を a + k d で作る
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const Point d = rndp(rng, it % 2 ? 3 : 1000);
        if (d == Point()) continue;
        const long long lim = 1000000000LL / std::max<long long>(1, std::max(std::llabs(d.x), std::llabs(d.y)));
        const Point a = Point(rnd(rng, -1000000000LL + lim * std::llabs(d.x), 1000000000LL - lim * std::llabs(d.x)),
                              rnd(rng, -1000000000LL + lim * std::llabs(d.y), 1000000000LL - lim * std::llabs(d.y)));
        long long k1 = rnd(rng, -lim, lim), k2 = rnd(rng, -lim, lim);
        if (it % 5 == 0) k2 = k1;
        if (it % 7 == 0) k2 = 0;
        if (k1 == 0) k1 = 1;
        const Point b = a + d * k1, c = a + d * k2;
        GEOB_CHECK(ccw(a, b, c) == 0);
        const int expect = (k1 > 0 ? k2 < 0 : k2 > 0) ? ONLINE_BACK : (std::llabs(k2) > std::llabs(k1) ? ONLINE_FRONT : ON_SEGMENT);
        GEOB_CHECK(ccw_detail(a, b, c) == expect);
        // 1 だけずらすと一直線でなくなる（外積の符号は 128 ビットと一致）
        const Point e = c + Point(rnd(rng, -1, 1), rnd(rng, -1, 1));
        if (std::llabs(e.x) <= 1000000000LL && std::llabs(e.y) <= 1000000000LL) {
            GEOB_CHECK(ccw(a, b, e) == sgn128(cross128(b - a, e - a)));
            GEOB_CHECK(ccw_detail(a, b, e) == brute_ccw_detail(a, b, e));
        }
        // 端の値どうし
        const Point f = rndp(rng, 1000000000LL), g = rndp(rng, 1000000000LL), h = rndp(rng, 1000000000LL);
        GEOB_CHECK(ccw(f, g, h) == sgn128(cross128(g - f, h - f)));
        if (f != g) GEOB_CHECK(ccw_detail(f, g, h) == brute_ccw_detail(f, g, h));
    }
    {
        const Point a(-1000000000LL, -1000000000LL), b(1000000000LL, 1000000000LL), c(1000000000LL, -1000000000LL), d(-1000000000LL, 1000000000LL);
        GEOB_CHECK(ccw(a, b, c) == -1 && ccw(a, b, d) == 1 && ccw(c, d, a) == ccw(c, d, b) * -1 + 0 * 0 - 0);
        GEOB_CHECK(ccw(a, b, Point(0, 0)) == 0 && ccw_detail(a, b, Point(0, 0)) == ON_SEGMENT);
        GEOB_CHECK(ccw_detail(Point(0, 0), b, a) == ONLINE_BACK && ccw_detail(a, Point(0, 0), b) == ONLINE_FRONT);
    }
    // 実数の点：回転・平行移動しても向きが変わらない
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 1000000 : 1000000000LL;
        const Real th = rndr(rng, 0, 2 * GEO_PI), sc = rndr(rng, 0.001L, 1000);
        const PointD sh(rndr(rng, -1e6L, 1e6L), rndr(rng, -1e6L, 1e6L));
        auto tf = [&](const Point& p) { return rotated(p, th) * sc + sh; };
        Point a = rndp(rng, C), b = rndp(rng, C), c = rndp(rng, C);
        if (it % 3 == 0) {  // 一直線
            const Point d = rndp(rng, 1000);
            const long long k1 = rnd(rng, -100, 100), k2 = rnd(rng, -100, 100);
            a = rndp(rng, C / 2), b = a + d * k1, c = a + d * k2;
        }
        if (a == b) continue;
        {  // Real の精度で区別できる場合だけ（座標の大きさ × 精度 が、ベクトルの長さ × GEO_EPS より十分小さい）
            const Real big = std::max({std::fabs(tf(a).x), std::fabs(tf(a).y), std::fabs(tf(b).x), std::fabs(tf(b).y), std::fabs(tf(c).x), std::fabs(tf(c).y)});
            const Real small = std::min(len(tf(b) - tf(a)), c == a ? len(tf(b) - tf(a)) : len(tf(c) - tf(a)));
            if (EPSR * 64 * std::max<Real>(big, 1) > GEO_EPS * 0.01L * small) continue;
        }
        const Int128 cr = cross128(b - a, c - a);
        const Real rel = (Real)cr / (len(b - a) * std::max<Real>(len(c - a), 1));
        if (cr == 0 || std::fabs(rel) > 1e-6L) {
            GEOB_CHECK(ccw(tf(a), tf(b), tf(c)) == sgn128(cr));
            GEOB_CHECK(ccw_detail(tf(a), tf(b), tf(c)) == brute_ccw_detail(a, b, c) || (c == a || c == b));
        }
    }
}

static void test_rotation(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 1000000000LL : 20;
        const Point p = rndp(rng, C / 2), c = rndp(rng, C / 2);
        GEOB_CHECK(rot90(p) == Point(-p.y, p.x) && rot90_cw(p) == Point(p.y, -p.x) && rot180(p) == -p);
        GEOB_CHECK(rot90(rot90_cw(p)) == p && rot90(rot90(p)) == rot180(p));
        GEOB_CHECK(dot(p, rot90(p)) == 0 && cross(p, rot90(p)) == norm2(p));
        const long long k = rnd(rng, -9, 9);
        Point e = p;
        for (long long t = 0; t < ((k % 4) + 4) % 4; ++t) e = rot90(e - c) + c;
        GEOB_CHECK(rot90(p, c, k) == e);
        GEOB_CHECK(rot90(rot90(p, c, k), c, -k) == p);
        // 任意の角：複素数の積と比べる
        const Real th = rndr(rng, -10, 10);
        const std::complex<Real> z = std::complex<Real>((Real)p.x, (Real)p.y) * std::polar((Real)1, th);
        GEOB_CHECK(nearp(rotated(p, th), PointD(z.real(), z.imag()), TOL));
        const std::complex<Real> w = (std::complex<Real>((Real)(p.x - c.x), (Real)(p.y - c.y))) * std::polar((Real)1, th) +
                                     std::complex<Real>((Real)c.x, (Real)c.y);
        GEOB_CHECK(nearp(rotated(p, th, c), PointD(w.real(), w.imag()), TOL));
        GEOB_CHECK(near(len(rotated(p, th)), len(p), TOL));
        // 度：90 度の倍数は厳密、それ以外は rotated と同じ
        const long long kk = rnd(rng, -8, 8);
        GEOB_CHECK(rotated_deg(p, (Real)(90 * kk)) == to_d(rot90(p, Point(), kk)));
        GEOB_CHECK(rotated_deg(p, (Real)(90 * kk), c) == to_d(rot90(p, c, kk)));
        const Real deg = rndr(rng, -720, 720);
        GEOB_CHECK(nearp(rotated_deg(p, deg), rotated(p, deg_to_rad(deg)), TOL));
        const PointD pd = to_d(p);
        GEOB_CHECK(rotated_deg(pd, 270) == PointD(pd.y, -pd.x));
        // 中心の型が違ってもよい（Point を PointD のまわりに）
        const PointD cd = to_d(c) + PointD(0.5L, -0.25L);
        GEOB_CHECK(nearp(rotated(p, th, cd), rotated(pd, th, cd), TOL) && nearp(rotated(pd, th, c), rotated(p, th, c), TOL));
        const std::complex<Real> w2 = (std::complex<Real>((Real)p.x - cd.x, (Real)p.y - cd.y)) * std::polar((Real)1, th) +
                                      std::complex<Real>(cd.x, cd.y);
        GEOB_CHECK(nearp(rotated(p, th, cd), PointD(w2.real(), w2.imag()), TOL));
        GEOB_CHECK(rotated_deg(p, (Real)(90 * kk), cd) == rot90(pd, cd, kk) + PointD(0, 0));
        GEOB_CHECK(nearp(rotated_deg(p, deg, cd), rotated(p, deg_to_rad(deg), cd), TOL));
    }
    // -0 を出さない（cout で "-0" と表示されないように）
    for (const PointD& q : {rotated(Point(0, 0), 1), rotated(Point(1, 0), GEO_PI / 2), rotated_deg(PointD(1, 0), 90), rotated_deg(PointD(0, -0.0L), 180),
                            rotated_deg(Point(1, 0), 90, PointD(0, 0)), rotated(PointD(-0.0L, 0), 0)}) {
        GEOB_CHECK(!std::signbit(q.x) || q.x != 0);
        GEOB_CHECK(!std::signbit(q.y) || q.y != 0);
    }
}

static void test_affine(std::mt19937_64& rng, int& cases) {
    // 整数のアフィン変換：合成と、1 つずつ適用した結果が一致する
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int m = (int)rnd(rng, 0, 8);
        const long long C = it % 2 ? 1000 : 3;
        std::vector<Affine> ops;
        std::vector<int> kind;
        std::vector<Point> arg;
        for (int k = 0; k < m; ++k) {
            const int t = (int)rnd(rng, 0, 6);
            const Point o = rndp(rng, C);
            kind.push_back(t), arg.push_back(o);
            if (t == 0) ops.push_back(Affine::rot90(o));
            if (t == 1) ops.push_back(Affine::rot90_cw(o));
            if (t == 2) ops.push_back(Affine::rot180(o));
            if (t == 3) ops.push_back(Affine::reflect_x(o.x));
            if (t == 4) ops.push_back(Affine::reflect_y(o.y));
            if (t == 5) ops.push_back(Affine::translate(o));
            if (t == 6) ops.push_back(Affine::scale(rnd(rng, -3, 3), o));
        }
        Affine f;
        for (const Affine& g : ops) f = f.then(g);
        const Point p = rndp(rng, C);
        Point q = p;
        for (int k = 0; k < m; ++k) {
            const Point o = arg[k];
            const int t = kind[k];
            if (t == 0) GEOB_CHECK(ops[k](q) == rot90(q, o, 1));
            if (t == 1) GEOB_CHECK(ops[k](q) == rot90(q, o, -1));
            if (t == 2) GEOB_CHECK(ops[k](q) == rot90(q, o, 2));
            if (t == 3) GEOB_CHECK(ops[k](q) == Point(2 * o.x - q.x, q.y));
            if (t == 4) GEOB_CHECK(ops[k](q) == Point(q.x, 2 * o.y - q.y));
            if (t == 5) GEOB_CHECK(ops[k](q) == q + o);
            if (t == 6) GEOB_CHECK(ops[k](q) - o == (q - o) * ops[k].a);
            q = ops[k](q);
        }
        GEOB_CHECK(f(p) == q);
        GEOB_CHECK(Affine()(p) == p && Affine().then(f) == f && f.then(Affine()) == f);
        // 結合法則
        if (m >= 3) GEOB_CHECK(ops[0].then(ops[1]).then(ops[2]) == ops[0].then(ops[1].then(ops[2])));
    }
    // 実数：任意の角の回転は rotated と同じ
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const Real th = rndr(rng, -10, 10), th2 = rndr(rng, -10, 10);
        const PointD o(rndr(rng, -100, 100), rndr(rng, -100, 100)), p(rndr(rng, -100, 100), rndr(rng, -100, 100));
        const AffineD f = AffineD::rotate(th, o), g = AffineD::rotate(th2), h = f.then(g).then(AffineD::scale(2, o));
        GEOB_CHECK(nearp(f(p), rotated(p, th, o), TOL * 100));
        GEOB_CHECK(nearp(h(p), (rotated(rotated(p, th, o), th2) - o) * (Real)2 + o, TOL * 1000));
    }
}

static void test_angles(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 3 : (it % 3 == 1 ? 1000 : 1000000000LL);
        const Point a = rndp(rng, C), b = rndp(rng, C);
        const Real t = angle_of(a);
        GEOB_CHECK(0 <= t && t < 2 * GEO_PI);
        if (a != Point()) {
            Real e = std::arg(std::complex<Real>((Real)a.x, (Real)a.y));
            if (e < 0) e += 2 * GEO_PI;
            GEOB_CHECK(std::fabs(t - e) < TOL || std::fabs(t - e) > 2 * GEO_PI - TOL);
            GEOB_CHECK(nearp(from_polar(len(a), t), to_d(a), TOL));
            GEOB_CHECK(near(len(unit(a)), 1, TOL) && ccw(PointD(), to_d(a), unit(a)) == 0 && dot(to_d(a), unit(a)) > 0);
        } else {
            GEOB_CHECK(t == 0 && unit(a) == PointD());
        }
        GEOB_CHECK(near(angle_of_deg(a), rad_to_deg(t), TOL));
        const Real ab = angle_between(a, b);
        GEOB_CHECK(0 <= ab && ab <= GEO_PI && angle_between(b, a) == ab);
        if (a == Point() || b == Point()) GEOB_CHECK(ab == 0 && angle_between(to_d(a), to_d(b)) == 0);
        if (a != Point() && b != Point()) {
            const Real cs = (Real)dot128(a, b) / (len(a) * len(b));
            if (std::fabs(cs) < 0.9L) GEOB_CHECK(near(ab, std::acos(cs), 1e-12L));
            if (cross128(a, b) == 0) GEOB_CHECK(dot128(a, b) > 0 ? std::fabs(ab) <= TOL : near(ab, GEO_PI, TOL));
            if (dot128(a, b) == 0) GEOB_CHECK(near(ab, GEO_PI / 2, TOL));
        }
        GEOB_CHECK(near(angle_between_deg(a, b), rad_to_deg(ab), TOL));
        const Real d = rndr(rng, -1000, 1000);
        GEOB_CHECK(near(rad_to_deg(deg_to_rad(d)), d, TOL));
    }
    GEOB_CHECK(angle_of(Point(1, 0)) == 0 && near(angle_of(Point(0, 1)), GEO_PI / 2, TOL) && near(angle_of(Point(-1, 0)), GEO_PI, TOL));
    GEOB_CHECK(near(angle_of(Point(0, -1)), 3 * GEO_PI / 2, TOL) && angle_of(PointD(1, -0.0L)) == 0);
    GEOB_CHECK(angle_of(PointD(1, -1e-30L)) < 2 * GEO_PI);
}

static void test_keys(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 6 : 1000000000LL;
        Point u = rndp(rng, C), v = rndp(rng, C);
        if (it % 4 == 0) {
            const Point d = rndp(rng, 5);
            u = d * rnd(rng, -100, 100), v = d * rnd(rng, -100, 100);
        }
        const Point ku = direction_key(u), kv = direction_key(v);
        if (u == Point()) {
            GEOB_CHECK(ku == Point() && line_direction_key(u) == Point());
            continue;
        }
        GEOB_CHECK(std::gcd(std::llabs(ku.x), std::llabs(ku.y)) == 1 && cross128(ku, u) == 0 && dot128(ku, u) > 0);
        if (v == Point()) continue;
        GEOB_CHECK((ku == kv) == (cross128(u, v) == 0 && dot128(u, v) > 0));
        GEOB_CHECK((line_direction_key(u) == line_direction_key(v)) == (cross128(u, v) == 0));
        const Point lk = line_direction_key(u);
        GEOB_CHECK(lk.x > 0 || (lk.x == 0 && lk.y > 0));
        GEOB_CHECK(direction_key(u * 7) == ku && direction_key(-u) == -ku);
    }
}

static void test_isqrt(std::mt19937_64& rng, int& cases) {
    auto check = [&](long long n) {
        const long long r = isqrt(n), c = isqrt_ceil(n);
        if (n <= 0) {
            GEOB_CHECK(r == 0 && c == 0);
            return;
        }
        GEOB_CHECK((Int128)r * r <= n && (Int128)(r + 1) * (r + 1) > n);
        GEOB_CHECK((Int128)c * c >= n && (Int128)(c - 1) * (c - 1) < n);
        ++cases;
    };
    for (long long n = -3; n <= 2000; ++n) check(n);
    for (int it = 0; it < 20000 * TEST_SCALE; ++it) {
        const int bits = (int)rnd(rng, 1, 63);
        const long long n = (long long)(rng() >> (64 - bits));
        check(n);
        const long long s = isqrt(n);
        for (long long d = -2; d <= 2; ++d) {
            const Int128 sq = (Int128)(s + d) * (s + d);
            if (sq >= 0 && sq <= (Int128)9223372036854775807LL) check((long long)sq), check((long long)sq - 1), check((long long)std::min<Int128>(sq + 1, (Int128)9223372036854775807LL));
        }
    }
    check(9223372036854775807LL);
    check(3037000499LL * 3037000499LL);
    check(3037000499LL * 3037000499LL - 1);
    GEOB_CHECK(isqrt(9223372036854775807LL) == 3037000499LL && isqrt_ceil(9223372036854775807LL) == 3037000500LL);
}

static void test_parse(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 6000 * TEST_SCALE; ++it, ++cases) {
        const int k = (int)rnd(rng, 0, 7);
        long long p10 = 1;
        for (int i = 0; i < k; ++i) p10 *= 10;
        const long long v = rnd(rng, -(it % 2 ? 1000LL : 1000000000000LL) * (k ? 1 : 1), (it % 2 ? 1000LL : 1000000000000LL));
        const long long a = std::llabs(v);
        std::string ip = std::to_string(a / p10), fp = std::to_string(a % p10);
        while ((int)fp.size() < k) fp = "0" + fp;
        if (k == 0) fp = "";
        const int style = (int)rnd(rng, 0, 6);
        std::string s;
        std::string f2 = fp;
        if (style == 1 || style == 2) {  // 末尾の 0 を消す
            while (!f2.empty() && f2.back() == '0') f2.pop_back();
        }
        if (style == 3) f2 += std::string((size_t)rnd(rng, 1, 3), '0');  // 余分な 0
        std::string i2 = ip;
        if (style == 4) i2 = std::string((size_t)rnd(rng, 1, 2), '0') + i2;  // 先頭の 0
        if (style == 5 && i2 == "0" && !f2.empty()) i2 = "";                  // ".25"
        s = i2;
        if (!f2.empty()) s += "." + f2;
        else if (style == 2) s += ".";  // "5."
        if (s.empty()) s = "0";
        if (v < 0) s = "-" + s;
        else if (style == 6) s = "+" + s;
        GEOB_CHECK(parse_fixed(s, k) == v);
    }
    GEOB_CHECK(parse_fixed("1.25", 2) == 125 && parse_fixed("-1.5", 2) == -150 && parse_fixed("3", 2) == 300);
    GEOB_CHECK(parse_fixed(".25", 2) == 25 && parse_fixed("-0.0001", 4) == -1 && parse_fixed("2.5000", 1) == 25);
    GEOB_CHECK(parse_fixed("-0", 3) == 0 && parse_fixed("100000.0000", 4) == 1000000000LL);
}

static void test_io_misc(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 20), k = (int)rnd(rng, 0, 4);
        std::vector<Point> ps(n), fs(n);
        std::string in;
        for (auto& p : ps) p = rndp(rng, 1000000000LL), in += std::to_string(p.x) + " " + std::to_string(p.y) + "\n";
        std::string in2;
        for (auto& p : fs) {
            p = rndp(rng, 1000000);
            auto fmt = [&](long long v) {
                long long p10 = 1;
                for (int i = 0; i < k; ++i) p10 *= 10;
                std::string f = std::to_string(std::llabs(v) % p10);
                while ((int)f.size() < k) f = "0" + f;
                return std::string(v < 0 ? "-" : "") + std::to_string(std::llabs(v) / p10) + (k ? "." + f : "");
            };
            in2 += fmt(p.x) + " " + fmt(p.y) + "\n";
        }
        {
            CinFromString cf(in);
            GEOB_CHECK(read_points(n) == ps);
        }
        {
            CinFromString cf(in);
            const std::vector<PointD> d = read_points_d(n);
            for (int i = 0; i < n; ++i) GEOB_CHECK(d[i] == to_d(ps[i]));
        }
        {
            CinFromString cf(in2);
            GEOB_CHECK(read_points_fixed(n, k) == fs);
        }
    }
    {
        CinFromString cf("1 2 3 4\n5 6 7 8\n-1 -2 10\n1.5 -2.25\n0.125\n3 4\n2.5 3.5\n");
        GEOB_CHECK(read_segment() == Segment(1, 2, 3, 4) && read_line() == Line(Point(5, 6), Point(7, 8)));
        GEOB_CHECK(read_circle() == Circle(-1, -2, 10));
        GEOB_CHECK(read_point_fixed(2) == Point(150, -225) && read_fixed(3) == 125 && read_point() == Point(3, 4));
        GEOB_CHECK(read_point_d() == PointD(2.5L, 3.5L));
    }
    // Int128 の文字列化
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        using U = geo_internal::UInt128;
        const Int128 v = (Int128)(((U)rng() << (int)rnd(rng, 0, 64)) ^ (U)rng());
        std::string s = int128_to_string(v);
        Int128 back = 0;
        for (size_t i = (s[0] == '-'); i < s.size(); ++i) back = back * 10 + (s[i] - '0');
        if (s[0] == '-') back = -back;
        GEOB_CHECK(back == v);
    }
    const Int128 mx = (Int128)(((geo_internal::UInt128)1 << 127) - 1);
    GEOB_CHECK(int128_to_string(mx) == "170141183460469231731687303715884105727");
    GEOB_CHECK(int128_to_string(-mx - 1) == "-170141183460469231731687303715884105728");
    GEOB_CHECK(int128_to_string(0) == "0" && int128_to_string(-7) == "-7");
    // ハッシュ
    {
        std::unordered_set<Point, PointHash> hs;
        std::set<Point> ss;
        for (int i = 0; i < 20000; ++i) {
            const Point p = rndp(rng, i % 2 ? 50 : 1000000000LL);
            hs.insert(p), ss.insert(p);
        }
        GEOB_CHECK(hs.size() == ss.size());
        for (const Point& p : ss) GEOB_CHECK(hs.count(p) == 1);
        ++cases;
    }
    // 実数の比較
    GEOB_CHECK(real_sign(1e-10L) == 0 && real_sign(1e-8L) == 1 && real_sign(-1e-8L) == -1);
    GEOB_CHECK(real_sign(1e-3L, 1e9L) == 0 && real_sign(2.0L, 1e9L) == 1);
    GEOB_CHECK(real_cmp(1e9L, 1e9L + 0.5L) == 0 && real_cmp(1, 1.001L) == -1 && real_eq(0.1L + 0.2L, 0.3L));
    GEOB_CHECK(real_eq(PointD(1e9L, 0.5L), PointD(1e9L + 0.1L, 0.5L)) && !real_eq(PointD(1, 2), PointD(1, 2.001L)));
    // 直線・線分・円の型
    {
        const Line l(Point(1, 2), Point(3, 5));
        const Segment s(0, 0, 4, 4);
        const Circle c(Point(1, 1), 3), c2(1, 1, 3);
        GEOB_CHECK(l.dir() == Point(2, 3) && s.dir() == Point(4, 4) && c == c2 && c.r == 3);
        GEOB_CHECK(to_d(l).a == PointD(1, 2) && to_d(s).b == PointD(4, 4) && to_d(c).r == 3 && to_d(c).c == PointD(1, 1));
        const LineD ld(0.5L, 1, 2, 3.5L);
        const CircleD cd(PointD(0.5L, 0.5L), 1.5L);
        GEOB_CHECK(ld.b == PointD(2, 3.5L) && cd.r == 1.5L);
        static_assert(!std::is_constructible_v<Circle, Point, double>, "整数の円に実数の半径は渡せない");
        static_assert(!std::is_constructible_v<Segment, double, int, int, int>, "整数の線分に実数は渡せない");
        std::ostringstream os;
        os << Point(1, -2) << " " << l << " " << s << " " << c;
        GEOB_CHECK(os.str() == "(1, -2) Line[(1, 2) -> (3, 5)] Segment[(0, 0) - (4, 4)] Circle[c=(1, 1), r=3]");
        GEOB_CHECK(OUTSIDE == 0 && ON_BORDER == 1 && INSIDE == 2);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_point_ops(rng, cases);
    test_ccw(rng, cases);
    test_rotation(rng, cases);
    test_affine(rng, cases);
    test_angles(rng, cases);
    test_keys(rng, cases);
    test_isqrt(rng, cases);
    test_parse(rng, cases);
    test_io_misc(rng, cases);
    std::cout << "[geo_base] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[geo_base] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 2000000;
    std::vector<Point> ps(n);
    for (auto& p : ps) p = rndp(rng, 1000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (int i = 0; i + 2 < n; ++i) s1 += ccw(ps[i], ps[i + 1], ps[i + 2]);
    const long long t_ccw = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for (int i = 0; i + 2 < n; ++i) s2 += ccw_detail(ps[i], ps[i + 1], ps[i + 2]);
    const long long t_ccwd = ms(t0);
    t0 = std::chrono::steady_clock::now();
    Real s3 = 0;
    for (int i = 0; i < n; ++i) s3 += rotated(ps[i], (Real)i).x;
    const long long t_rot = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s4 = 0;
    for (int i = 0; i < n; ++i) s4 += isqrt(ps[i].x * ps[i].x);
    const long long t_sqrt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s5 = 0;
    for (int i = 0; i + 1 < n; ++i) s5 += direction_key(ps[i + 1] - ps[i]).x;
    const long long t_key = ms(t0);
    const int m = 200000;
    std::string in;
    for (int i = 0; i < m; ++i) in += std::to_string(ps[i].x) + " " + std::to_string(ps[i].y) + "\n";
    std::string inf;
    for (int i = 0; i < m; ++i) inf += std::to_string(ps[i].x % 1000) + "." + std::to_string(std::llabs(ps[i].y) % 10000) + " 1.5\n";
    std::ios::sync_with_stdio(false);
    t0 = std::chrono::steady_clock::now();
    std::vector<Point> rp;
    {
        CinFromString cf(in);
        rp = read_points(m);
    }
    const long long t_read = ms(t0);
    t0 = std::chrono::steady_clock::now();
    std::vector<Point> rf;
    {
        CinFromString cf(inf);
        rf = read_points_fixed(m, 4);
    }
    const long long t_readf = ms(t0);
    std::cout << "[geo_base] benchmark: " << n << " 回  ccw " << t_ccw << " ms, ccw_detail " << t_ccwd << " ms, rotated " << t_rot
              << " ms, isqrt " << t_sqrt << " ms, direction_key " << t_key << " ms\n"
              << "[geo_base] benchmark: N = " << m << "  read_points " << t_read << " ms, read_points_fixed(N, 4) " << t_readf
              << " ms  (check " << s1 + s2 + (long long)s3 % 7 + s4 % 7 + s5 % 7 + rp[m - 1].x % 7 + rf[m - 1].y % 7 << ")\n";
}

}  // namespace geob_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace geob_demo {

static std::string to_s(bool b) { return b ? "true" : "false"; }
template <class T, std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool>, int> = 0>
static std::string to_s(T v) { return std::to_string(v); }
static std::string to_s(Real v) {
    if (std::fabs(v) < 5e-7L) v = 0;
    std::ostringstream os;
    os << std::fixed << std::setprecision(6) << v;
    return os.str();
}
static std::string to_s(const Point& p) { return "(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ")"; }
static std::string to_s(const PointD& p) { return "(" + to_s(p.x) + ", " + to_s(p.y) + ")"; }
static std::string to_s(const std::string& s) { return s; }
static std::string pad(const std::string& s, size_t w) {
    size_t cols = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = (unsigned char)s[i];
        if ((c & 0xC0) == 0x80) continue;
        cols += (c >= 0xE0 ? 2 : 1);
    }
    return cols >= w ? s : s + std::string(w - cols, ' ');
}
static void banner(const std::string& s) {
    std::cout << "\n=====================================================================\n"
              << "  " << s << "\n"
              << "=====================================================================\n";
}
static void api(const std::string& usage, const std::string& ret, const std::string& what, const std::string& cost) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << usage << "\n"
              << "    返り値 : " << ret << "\n"
              << "    内容   : " << what << "\n"
              << "    計算量 : " << cost << "\n";
}
static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
static void stmt(const std::string& s) { code(s + ";"); }
template <class V>
static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 34) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 58) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) geob_demo::show(#expr, (expr), comment)
#define RUN(...)                    \
    geob_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("幾何基盤 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧（型）\n"
              << "  Point p(x, y)  / PointD q(x, y)          整数の点（long long）/ 実数の点（long double）\n"
              << "  Line l(a, b) / Segment s(a, b) / Circle c(中心, r)   直線 a→b / 線分 / 円（LineD / SegmentD / CircleD も）\n"
              << "  to_d(p) / to_d(l) / to_d(s) / to_d(c)    実数の型に直す（Point → PointD は自動では変わらない）\n"
              << "  Where: OUTSIDE / ON_BORDER / INSIDE       内外判定の結果（03_polygon・05_circle などで使う）\n"
              << "■ API 一覧（関数）\n"
              << "  dot / cross / norm2 / dist2              内積 / 外積 / 長さの 2 乗 / 距離の 2 乗（整数のまま）\n"
              << "  len / dist                               長さ / 距離（Real）\n"
              << "  ccw(a, b, c) / ccw_detail / collinear    向き（+1 左 / -1 右 / 0 一直線）/ 5 分類 / 一直線か\n"
              << "  rot90 / rot90_cw / rot180 / rot90(p, c, k)     90 度回転（整数のまま）\n"
              << "  rotated(p, θ[, c]) / rotated_deg(p, 度[, c])   任意の角の回転（PointD。c は PointD でもよい）\n"
              << "  Affine / AffineD（f(p)・f.then(g)・rot90 / rot90_cw / rot180 / reflect_x / reflect_y / translate / scale / rotate）\n"
              << "  angle_of / angle_between（_deg も）     偏角 [0, 2π) / なす角 [0, π]\n"
              << "  direction_key / line_direction_key       向き / 傾きのキー（gcd で割る）\n"
              << "  isqrt / isqrt_ceil                       整数の平方根（切り捨て / 切り上げ）\n"
              << "  parse_fixed(s, k) / read_fixed(k)        小数を 10^k 倍の整数に（誤差なし）\n"
              << "  read_point / read_points(n) / read_points_fixed(n, k) / read_segment / read_circle   入力\n"
              << "  real_sign / real_cmp / real_eq           実数の比較（相対誤差 GEO_EPS = 1e-9）\n"
              << "  int128_to_string / PointHash             Int128 の出力 / unordered_set 用のハッシュ\n";

    api("【1】 Point p(x, y);  p.x, p.y;  + - * / == <", "Point（整数）",
        "整数の点。足し算・引き算・整数倍・比較（< は x → y の辞書順。sort や set に使える）", "O(1)");
    RUN(Point a(3, 1); Point b(1, 2));
    SHOW(a + b, "足し算");
    SHOW(a - b, "引き算（b から a へのベクトル）");
    SHOW(a * 3, "整数倍（a * 1.5 はコンパイルエラー。実数にするなら to_d(a) * 1.5）");
    SHOW(b < a, "x が小さい方が小さい");
    SHOW(to_d(a) * (Real)0.5, "実数の点にしてから半分");
    note("構造化束縛もできる: auto [x, y] = p;");

    api("【2】 dot(a, b) / cross(a, b) / norm2(a) / dist2(a, b) / dist(a, b)", "Coord（整数）/ Real",
        "内積・外積・長さの 2 乗・距離の 2 乗は整数のまま（誤差なし）。距離そのものは Real", "O(1)");
    SHOW(dot(a, b), "3*1 + 1*2");
    SHOW(cross(a, b), "3*2 - 1*1。正なら b は a の左（反時計回り側）");
    SHOW(cross(Point(0, 0), a, b), "cross(o, a, b) = cross(a - o, b - o)。三角形 o a b の面積の 2 倍");
    SHOW(dist2(a, b), "距離の 2 乗（比較はこれで。誤差なし）");
    SHOW(dist(a, b), "√5");
    note("|座標| ≤ 10^9 なら差の積（cross(b - a, c - a) など）も long long に収まる。");
    note("それ以上の大きさの積を比べたいときは dot128 / cross128 / norm2_128（Int128）を使う。");

    api("【3】 ccw(a, b, c) / ccw_detail(a, b, c) / collinear(a, b, c)", "int",
        "a → b → c が左に曲がる +1 / 右 -1 / 一直線 0。ccw_detail は一直線のときの並び方まで返す", "O(1)");
    SHOW(ccw(Point(0, 0), Point(2, 0), Point(1, 1)), "左（反時計回り）");
    SHOW(ccw(Point(0, 0), Point(2, 0), Point(1, -1)), "右（時計回り）");
    SHOW(ccw(Point(0, 0), Point(2, 0), Point(5, 0)), "一直線");
    SHOW(ccw_detail(Point(0, 0), Point(2, 0), Point(5, 0)) == ONLINE_FRONT, "a - b - c の順（c は b の先）");
    SHOW(ccw_detail(Point(0, 0), Point(2, 0), Point(-1, 0)) == ONLINE_BACK, "c - a - b の順（c は a の後ろ）");
    SHOW(ccw_detail(Point(0, 0), Point(2, 0), Point(1, 0)) == ON_SEGMENT, "c は線分 ab の上（端点も含む）");
    note("ccw_detail の値: COUNTER_CLOCKWISE(1) / CLOCKWISE(-1) / ONLINE_BACK(2) / ONLINE_FRONT(-2) / ON_SEGMENT(0)");
    note("ccw_detail(...) > 0 と書くと ONLINE_BACK も入るので、名前で比べる。向きだけなら ccw を使う。");
    note("類題: ABC266 C（凸四角形：4 つの角で ccw がすべて +1）、ABC181 C（3 点が一直線か）");

    api("【4】 rot90(p) / rot90_cw(p) / rot90(p, c, k) / rotated(p, θ) / rotated_deg(p, 度, c)", "Point / PointD",
        "90 度の回転は整数のまま。任意の角は PointD（90 度の倍数を度で渡すと誤差なし）", "O(1)");
    SHOW(rot90(Point(2, 1)), "反時計回りに 90 度");
    SHOW(rot90_cw(Point(2, 1)), "時計回りに 90 度");
    SHOW(rot90(Point(3, 1), Point(1, 1), 2), "(1, 1) のまわりに 180 度");
    SHOW(rotated_deg(Point(2, 2), 180), "ABC259 B: (2, 2) を 180 度 → (-2, -2)（誤差なし）");
    SHOW(rotated_deg(Point(5, 0), 120), "ABC259 B: (5, 0) を 120 度");
    SHOW(rotated(Point(1, 0), GEO_PI / 3, Point(0, 0)), "60 度（ラジアン）");
    note("ABC259 B の出力: auto q = rotated_deg(Point(a, b), d); cout << fixed << setprecision(15) << q.x << ' ' << q.y;");
    SHOW(rotated(Point(1, 1), GEO_PI / 2, PointD(0.5L, 0.5L)), "中心は PointD でもよい（(0.5, 0.5) のまわりに 90 度）");
    note("ABC197 D（正 N 角形の反対の頂点）: PointD c = midpoint(p0, pN2); PointD p1 = rotated(p0, 2 * GEO_PI / N, c);");
    note("    0 に近い値が -0.000000 と表示されることがあるが、誤差の範囲なので正解になる（90 度の倍数なら rotated_deg が厳密）");

    api("【4-2】 Affine f / AffineD（f(p)・f.then(g)・Affine::rot90(o) / rot90_cw(o) / rot180(o) / reflect_x(p) / reflect_y(p) / translate(v) / scale(k, o)）",
        "Affine（a x + b y + c, d x + e y + f）", "回転・折り返し・平行移動を 1 つの変換にまとめる。整数のまま合成できる（誤差なし）", "O(1)");
    RUN(Affine f = Affine::rot90_cw().then(Affine::reflect_x(3)));
    SHOW(f(Point(1, 2)), "時計回り 90 度 → (2, -1)、x = 3 で折り返し → (4, -1)");
    SHOW(Affine::rot90(Point(1, 1))(Point(3, 1)), "(1, 1) のまわりに反時計回り 90 度");
    SHOW(AffineD::rotate(GEO_PI / 2, PointD(1, 1))(PointD(3, 1)), "実数は任意の角（AffineD::rotate）");
    note("ABC189 E（操作の列とクエリ「A 回目の操作のあとの点 B」）: pre[0] = Affine(); pre[i + 1] = pre[i].then(操作 i);");
    note("    答えは pre[A](P[B])。操作 1 = Affine::rot90_cw()、2 = Affine::rot90()、3 = reflect_x(p)、4 = reflect_y(p)");
    note("reflect_x(p) は「直線 x = p で折り返す」（x → 2p - x）。x 軸で折り返す（y → -y）なら reflect_y(0)");

    api("【5】 angle_of(v) / angle_between(a, b)（_deg も）/ from_polar(r, θ) / unit(v)", "Real / PointD",
        "偏角 [0, 2π) / 2 つのベクトルのなす角 [0, π] / 長さ r・角 θ の点 / 長さ 1 のベクトル", "O(1)");
    SHOW(angle_of_deg(Point(-1, 1)), "135 度");
    SHOW(angle_of_deg(Point(0, -2)), "270 度（負にならない）");
    SHOW(angle_between_deg(Point(1, 0), Point(-1, 1)), "135 度");
    SHOW(unit(Point(3, 4)), "長さ 1");
    {
        // ABC168 C: 長さ A の時針と長さ B の分針の先端の距離（H 時 M 分）
        const Real A = 3, B = 4;
        const int H = 9, M = 0;
        const PointD h = from_polar(A, deg_to_rad(90 - 30 * H - 0.5L * M)), m = from_polar(B, deg_to_rad(90 - 6 * M));
        note("ABC168 C（時計の針）: h = from_polar(A, deg_to_rad(90 - 30H - 0.5M)), m = from_polar(B, deg_to_rad(90 - 6M)); dist(h, m)");
        std::cout << "    A = 3, B = 4, 9 時 0 分 → 距離 " << to_s(dist(h, m)) << "（直角なので 5）\n";
    }

    api("【6】 direction_key(v) / line_direction_key(v)", "Point",
        "向き・傾きが同じベクトルを同じキーにする（gcd で割る）。set や map で数えるときに", "O(log)");
    SHOW(direction_key(Point(4, -6)), "gcd 2 で割る");
    SHOW(direction_key(Point(-4, 6)), "向きが逆なら別のキー");
    SHOW(line_direction_key(Point(-4, 6)), "向きを区別しない（傾きのキー。x > 0 にそろえる）");
    note("ABC226 D: 全部の組 (i, j) の direction_key(ps[j] - ps[i]) を set に入れた個数が答え");
    note("同じ直線（傾きだけでなく位置も）のキーは 02_line の line_key(a, b)");

    api("【7】 isqrt(n) / isqrt_ceil(n)", "long long", "⌊√n⌋ / ⌈√n⌉ を誤差なしで。距離の比較・切り上げを整数で", "O(1)");
    SHOW(isqrt(99), "√99 = 9.94... の切り捨て");
    SHOW(isqrt_ceil(99), "切り上げ");
    SHOW(isqrt(1000000000000000000LL), "10^18 でも誤差なし");
    {
        // ABC198 C: 1 歩 R で (0, 0) → (X, Y) へ。最小の歩数
        const long long R = 5, X = 15, Y = 0;
        const long long d2 = X * X + Y * Y;
        const long long steps = d2 == 0 ? 0 : (d2 < R * R ? 2 : (isqrt_ceil(d2) + R - 1) / R);
        note("ABC198 C（1 歩 R で (X, Y) へ）: d2 = X² + Y²。d2 < R² なら 2、それ以外は ⌈⌈√d2⌉ / R⌉");
        std::cout << "    R = 5, (X, Y) = (15, 0) → " << steps << " 歩\n";
    }

    api("【8】 parse_fixed(s, k) / read_fixed(k) / read_points_fixed(n, k)", "long long / Point",
        "小数点以下 k 桁までの小数を 10^k 倍した整数にする（誤差なし）。小数の入力はこれで整数にしてから解く", "O(桁数)");
    SHOW(parse_fixed("1.25", 4), "1.25 × 10^4");
    SHOW(parse_fixed("-0.5", 4), "負の数も");
    SHOW(parse_fixed("3", 4), "小数点がなくてもよい");
    note("ABC191 D（円の中の格子点）: X, Y, R が小数点以下 4 桁 → X = read_fixed(4) などで 10^4 倍の整数にする。");
    note("    格子点 (x, y) は 10^4 の倍数の点になる（数えるのは 05_circle の lattice_points_in_circle）");
    note("小数点以下に k 桁を超える 0 でない桁があると止まる（黙って切り捨てない）。");

    api("【9】 Line / Segment / Circle、to_d、read_segment / read_line / read_circle", "型",
        "直線は a → b の向きを持つ（左側 = cross(b - a, p - a) > 0）。アルゴリズムは 02_line・05_circle など", "O(1)");
    RUN(Line l(Point(0, 0), Point(4, 2)); Segment s(1, 1, 3, 5); Circle c(Point(0, 0), 5));
    SHOW(l.dir(), "向きのベクトル b - a");
    SHOW(to_d(c).r, "実数の円に");
    std::cout << "  cout << l << s << c  ->  " << l << " " << s << " " << c << "\n";
    note("入力が x1 y1 x2 y2 なら read_segment()、x y r なら read_circle()。");

    api("【10】 Int128 / int128_to_string(v) / PointHash", "Int128 / string / ハッシュ",
        "Int128 は __int128（cout で直接出せないので文字列にする）。点を unordered_set に入れるハッシュ", "O(1)");
    SHOW(int128_to_string(cross128(Point(1000000000, -1000000000), Point(1000000000, 1000000000)) * 1000000000), "2×10^27");
    RUN(std::unordered_set<Point, PointHash> seen; seen.insert(Point(1, 2)); seen.insert(Point(1, 2)));
    SHOW(seen.size(), "同じ点は 1 つ");

    api("【11】 real_sign(x, scale) / real_cmp(a, b) / real_eq(a, b)", "int / bool",
        "実数の比較（相対誤差 GEO_EPS = 1e-9。大きさ 1 未満は絶対誤差）。点どうしの real_eq もある", "O(1)");
    SHOW(real_eq(0.1L + 0.2L, 0.3L), "誤差を無視して等しい");
    SHOW(real_cmp(1e9L, 1e9L + 0.5L), "10^9 に対して 0.5 は誤差の範囲（相対 1e-9）");
    SHOW(real_sign(-1e-12L), "0 とみなす");
    note("実数の点の向き ccw(PointD...) も、ベクトルの長さに合わせた相対誤差で判定する。");

    head("【12】 グリッドの (i, j) と幾何の (x, y)・実数の出力");
    note("幾何の (x, y) は「横・縦（上が正）」、グリッドの (i, j) は「行・列（下が正）」。");
    note("グリッドのマスを幾何で扱うときは Point{j, -i} にすると、見た目の左回り = 反時計回りになる。");
    note("マスのまま 90 度回すならグラフ・グリッドの 10_grid_transform（rotate_point_cw など）。");
    note("実数の答えは cout << fixed << setprecision(15) << x;（Main_base なら cout << FIXED(15) << x;）");
}

}  // namespace geob_demo


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
    geob_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    geob_test::benchmark();    // 速度計測
    geob_demo::demo();         // 使い方デモ
    return 0;
}
