// =============================================================================
//  3 次元（点・ベクトル・四面体の体積・角度・回転・直線・線分・平面・球）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 10 3 次元
//  依存       : 00_int_common, 00_geo_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - 点は整数の Point3（Coord = long long）と実数の Point3D（Real）。2 次元の Point / PointD と同じ使い方（+ - * / == < など）
//    - 整数の型なら判定（向き・同じ平面・一直線・交わるか・内外・球の位置関係・交点の数）は誤差なし。|x|, |y|, |z| ≤ 10^9 まで正しい
//      ただし dot / norm2 / dist2 は long long で返すので、座標の差が 1.7×10^9 を超えるとはみ出す → dot128 / norm2_128 / dist2_128
//    - 交点・射影・距離・角度・体積は実数で返す（交点は vector<Point3D>。なしは空）。実数の型は相対誤差 GEO_EPS で判定する
//    - 直線 Line3 l(a, b) は a → b（a ≠ b）、線分 Segment3 s(a, b)（a == b の点も可）
//    - 平面 Plane3 pl(a, b, c, d) は a x + b y + c z = d（法線 n = (a, b, c) ≠ 0。整数なら d は Int128）
//      plane_from_points(p, q, r)：3 点を通る平面（p → q → r が反時計回りに見える側が表）。表側（side = +1）は法線の向き
//    - 回転はすべて右ねじの向き（軸の先から見て反時計回り）
//
//  機能一覧                                                                    計算量
//    Point3 p(x, y, z) / Point3D       演算子、cout << p は (x, y, z)。to_d / midpoint / unit
//    dot / cross（ベクトル）/ norm2 / dist2 / len / dist、dot128 / norm2_128 / dist2_128（Int128）     O(1)
//    volume6(a, b, c, d)               符号付きの四面体の体積の 6 倍（整数なら Int128 で厳密）/ tetra_volume（Real）
//    orient3d(a, b, c, d)              d が平面 abc の表側 +1 / 裏側 -1 / 同じ平面 0。coplanar(a, b, c, d) / collinear(a, b, c)
//    triangle_area(a, b, c)            三角形の面積（Real）
//    in_tetrahedron(a, b, c, d, p)     OUTSIDE / ON_BORDER / INSIDE（四面体はつぶれていないこと）
//    angle_between(u, v) / _deg        ベクトルのなす角 [0, π]。直線どうし・平面どうし・直線と平面は [0, π/2]
//    elevation_angle(v) / _deg         仰角（xy 平面から上へ測った角。[-π/2, π/2]。下向きなら負 = 俯角）
//    rot90_x / rot90_y / rot90_z(p)    座標軸のまわりに 90 度（整数のまま）/ rotation24(p, k)：立方体の 24 通りの向きの k 番目
//    rotated(p, θ, 軸 [, 中心]) / rotated_deg(p, 度, 軸 [, 中心])   任意の軸のまわりの回転（Point3D。座標軸の 90 度の倍数は誤差なし。型が違ってもよい）
//    on_line / on_segment / is_parallel / is_orthogonal / intersects（直線・線分のどの組でも）
//    projection(l or pl, p) / reflection(l or pl, p)    垂線の足 / 対称な点（Point3D）
//    dist(図形, 図形)                  点と点・直線・線分・平面、直線・線分どうし（Real。引数の順は逆でもよい）
//    side(pl, p) / on_plane(pl, p)     平面の表 +1 / 裏 -1 / 上 0
//    cross_point(l or s, pl)           直線・線分と平面の交点（直線が平面に平行・平面の上なら空。線分が平面の上なら両端）
//    cross_line(pl1, pl2)              2 平面の交線（vector<Line3D>。平行なら空。a は原点に最も近い点、向きは pl1.n × pl2.n）
//    Sphere sp(c, r)                   contains / sphere_relation（SPHERE_*）/ cross_points(sp, l) / intersection_volume(sp1, sp2)
//    sphere_cap_volume(r, h)           高さ h の球冠の体積
//   位置関係（整数なら厳密）
//    line_relation(l, m)               LINE3_CROSS（交わる）/ LINE3_PARALLEL / LINE3_SAME（同じ直線）/ LINE3_SKEW（ねじれの位置）
//    cross_point(l, m) / cross_point(l, s) / cross_points(s, t)   直線・線分どうしの交点（重なるなら両端）
//    cross_point_frac(l, m)            交点を既約分数 FracPoint3{x, y, z, d} で（異なる交点を set / FracPoint3Hash で数える）
//    direction_key(v) / line_direction_key(v)   gcd で割った向き / 向きを区別しない（傾きのキー）
//    line_key(p, q) / line_key(l)      3 次元の直線のキー LineKey3{d, m}（同じ直線なら同じ。==・<・LineKey3Hash）
//    line_plane_relation(l, pl)        LINE_PLANE_CROSS / LINE_PLANE_PARALLEL / LINE_PLANE_ON（直線が平面の上）
//    intersects(l or s, pl) / dist(l or s, pl)   直線・線分と平面が共有点を持つか / 距離
//    plane_relation(p, q) / plane_key(pl) / dist(p, q)   PLANE_CROSS / PLANE_PARALLEL / PLANE_SAME / 同じ平面なら同じキー / 平行な平面の間隔
//    sphere_plane_relation(sp, pl)     球面と平面：0 離れる / 1 接する / 2 交わる（整数なら 256 ビットで厳密）
//    cross_circle(sp, pl) / cross_circle(sp1, sp2)   交わりの円 Circle3D{c, n（単位法線）, r}（接するなら r = 0、なければ空）
//    ball_intersects(sp, 線分・直線・平面・球) / ball_contains(sp, 線分・球) / ball_dist(sp, 点・線分・直線・平面・球)   球体（中身を含む）として
//    read_point3() / read_points3(N) / read_point3_d() / read_points3_d(N) / Point3Hash
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    g3_test::random_test()  … 行列式・連立方程式（有理数）・三分探索・数値積分の愚直解と突き合わせるランダムテスト
//    g3_test::benchmark()    … 速度計測
//    g3_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_GEO3D
#define GG_GEO3D
// ---- 3 次元（幾何第二世代）  Point3（整数）/ Point3D（実数）、直線・線分・平面・球。判定は整数で厳密 ----
template <class T>
struct Point3T {
    T x = 0, y = 0, z = 0;
    constexpr Point3T() {}
    template <class A, class B, class C, std::enable_if_t<geo_internal::coord_args_ok<T, A, B, C>, int> = 0>
    constexpr Point3T(A x_, B y_, C z_) : x(static_cast<T>(x_)), y(static_cast<T>(y_)), z(static_cast<T>(z_)) {}
    template <class U, std::enable_if_t<!std::is_same_v<T, U> && std::is_floating_point_v<T>, int> = 0>
    constexpr explicit Point3T(const Point3T<U>& p) : x(static_cast<T>(p.x)), y(static_cast<T>(p.y)), z(static_cast<T>(p.z)) {}

    constexpr Point3T operator+() const { return *this; }
    constexpr Point3T operator-() const { return Point3T(-x, -y, -z); }
    constexpr Point3T& operator+=(const Point3T& p) { x += p.x, y += p.y, z += p.z; return *this; }
    constexpr Point3T& operator-=(const Point3T& p) { x -= p.x, y -= p.y, z -= p.z; return *this; }
    constexpr Point3T& operator*=(T k) { x *= k, y *= k, z *= k; return *this; }
    constexpr Point3T& operator/=(T k) { x /= k, y /= k, z /= k; return *this; }  // 整数の点は切り捨て（割り切れるときに使う）
    template <class K, std::enable_if_t<std::is_floating_point_v<K> && std::is_integral_v<T>, int> = 0>
    Point3T& operator*=(K) = delete;  // 整数の点に実数は掛けない（to_d(p) を使う）
    template <class K, std::enable_if_t<std::is_floating_point_v<K> && std::is_integral_v<T>, int> = 0>
    Point3T& operator/=(K) = delete;

    friend constexpr Point3T operator+(Point3T a, const Point3T& b) { return a += b; }
    friend constexpr Point3T operator-(Point3T a, const Point3T& b) { return a -= b; }
    friend constexpr Point3T operator*(Point3T a, T k) { return a *= k; }
    friend constexpr Point3T operator*(T k, Point3T a) { return a *= k; }
    friend constexpr Point3T operator/(Point3T a, T k) { return a /= k; }
    template <class K>  // 整数の点に実数は掛けない・割らない（コンパイルエラーにする。to_d(p) を使う）
    friend std::enable_if_t<std::is_floating_point_v<K> && std::is_integral_v<T>, Point3T> operator*(const Point3T&, K) = delete;
    template <class K>
    friend std::enable_if_t<std::is_floating_point_v<K> && std::is_integral_v<T>, Point3T> operator*(K, const Point3T&) = delete;
    template <class K>
    friend std::enable_if_t<std::is_floating_point_v<K> && std::is_integral_v<T>, Point3T> operator/(const Point3T&, K) = delete;

    friend constexpr bool operator==(const Point3T& a, const Point3T& b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
    friend constexpr bool operator!=(const Point3T& a, const Point3T& b) { return !(a == b); }
    friend constexpr bool operator<(const Point3T& a, const Point3T& b) { return a.x != b.x ? a.x < b.x : (a.y != b.y ? a.y < b.y : a.z < b.z); }
    friend constexpr bool operator>(const Point3T& a, const Point3T& b) { return b < a; }
    friend constexpr bool operator<=(const Point3T& a, const Point3T& b) { return !(b < a); }
    friend constexpr bool operator>=(const Point3T& a, const Point3T& b) { return !(a < b); }
};
using Point3 = Point3T<Coord>;
using Point3D = Point3T<Real>;
template <class T>
std::ostream& operator<<(std::ostream& os, const Point3T<T>& p) { return os << "(" << p.x << ", " << p.y << ", " << p.z << ")"; }
template <class T>
std::istream& operator>>(std::istream& is, Point3T<T>& p) { return is >> p.x >> p.y >> p.z; }

// ---- 内積・外積・長さ ----
template <class T>
constexpr T dot(const Point3T<T>& a, const Point3T<T>& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
template <class T>
constexpr Point3T<T> cross(const Point3T<T>& a, const Point3T<T>& b) {  // 外積（ベクトル）。整数なら成分は差が 2×10^9 まで long long に収まる
    return Point3T<T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
template <class T>
constexpr T norm2(const Point3T<T>& a) { return dot(a, a); }
template <class T>
constexpr T dist2(const Point3T<T>& a, const Point3T<T>& b) { return norm2(a - b); }
inline Int128 dot128(const Point3& a, const Point3& b) { return (Int128)a.x * b.x + (Int128)a.y * b.y + (Int128)a.z * b.z; }
inline Int128 norm2_128(const Point3& a) { return dot128(a, a); }
inline Int128 dist2_128(const Point3& a, const Point3& b) { return norm2_128(a - b); }
template <class T>
Real len(const Point3T<T>& a) { return std::sqrt((Real)a.x * (Real)a.x + (Real)a.y * (Real)a.y + (Real)a.z * (Real)a.z); }
template <class T>
Real dist(const Point3T<T>& a, const Point3T<T>& b) { return len(a - b); }
template <class T>
constexpr Point3D to_d(const Point3T<T>& p) { return Point3D(p.x, p.y, p.z); }
template <class T>
Point3D midpoint(const Point3T<T>& a, const Point3T<T>& b) { return Point3D(((Real)a.x + (Real)b.x) / 2, ((Real)a.y + (Real)b.y) / 2, ((Real)a.z + (Real)b.z) / 2); }
template <class T>
Point3D unit(const Point3T<T>& v) {  // 長さ 1 に（0 ベクトルは 0 のまま）
    const Real l = len(v);
    return l == 0 ? Point3D() : to_d(v) / l;
}

namespace geo_internal {
template <class T>
using wide3_t = std::conditional_t<std::is_integral_v<T>, Int128, Real>;  // 途中の計算の型（整数なら 128 ビット）
template <class T>
wide3_t<T> dotw(const Point3T<T>& a, const Point3T<T>& b) {
    using W = wide3_t<T>;
    return (W)a.x * b.x + (W)a.y * b.y + (W)a.z * b.z;
}
template <class T>
bool is_zero3(const Point3T<T>& v, Real scale) {  // 0 ベクトルか（実数は scale に合わせた相対誤差）
    if constexpr (std::is_integral_v<T>) return v.x == 0 && v.y == 0 && v.z == 0;
    else return real_sign(len(v), scale) == 0;
}
inline Real fold_half_pi(Real a) { return std::min(a, GEO_PI - a); }  // [0, π] → [0, π/2]（向きを区別しない角）
}  // namespace geo_internal

// ---- 体積・向き ----
// 符号付きの四面体の体積の 6 倍 = (b - a, c - a, d - a) の行列式。d が平面 abc の表側（a → b → c が反時計回りに見える側）なら正
template <class T>
geo_internal::wide3_t<T> volume6(const Point3T<T>& a, const Point3T<T>& b, const Point3T<T>& c, const Point3T<T>& d) {
    return geo_internal::dotw(cross(b - a, c - a), d - a);
}
template <class T>
Real tetra_volume(const Point3T<T>& a, const Point3T<T>& b, const Point3T<T>& c, const Point3T<T>& d) { return std::fabs((Real)volume6(a, b, c, d)) / 6; }
template <class T>
int orient3d(const Point3T<T>& a, const Point3T<T>& b, const Point3T<T>& c, const Point3T<T>& d) {  // +1 表側 / -1 裏側 / 0 同じ平面
    const auto v = volume6(a, b, c, d);
    if constexpr (std::is_integral_v<T>) return (v > 0) - (v < 0);
    else return real_sign(v, len(b - a) * len(c - a) * len(d - a));
}
template <class T>
bool coplanar(const Point3T<T>& a, const Point3T<T>& b, const Point3T<T>& c, const Point3T<T>& d) { return orient3d(a, b, c, d) == 0; }
template <class T>
bool collinear(const Point3T<T>& a, const Point3T<T>& b, const Point3T<T>& c) {  // 3 点が一直線（同じ点を含むときも）
    const Point3T<T> u = b - a, v = c - a;
    return geo_internal::is_zero3(cross(u, v), len(u) * len(v));
}
template <class T>
Real triangle_area(const Point3T<T>& a, const Point3T<T>& b, const Point3T<T>& c) { return len(cross(to_d(b - a), to_d(c - a))) / 2; }
// 四面体 abcd（体積が 0 でないこと。0 なら OUTSIDE）の内外。面・辺・頂点の上は ON_BORDER
template <class T>
Where in_tetrahedron(const Point3T<T>& a, const Point3T<T>& b, const Point3T<T>& c, const Point3T<T>& d, const Point3T<T>& p) {
    const int o = orient3d(a, b, c, d);
    if (o == 0) return OUTSIDE;
    const int s[4] = {orient3d(p, b, c, d) * o, orient3d(a, p, c, d) * o, orient3d(a, b, p, d) * o, orient3d(a, b, c, p) * o};
    bool zero = false;
    for (int k = 0; k < 4; ++k) {
        if (s[k] < 0) return OUTSIDE;
        zero |= s[k] == 0;
    }
    return zero ? ON_BORDER : INSIDE;
}

// ---- 角度 ----
template <class T>
Real angle_between(const Point3T<T>& u, const Point3T<T>& v) {  // なす角 [0, π]（どちらかが 0 ベクトルなら 0）
    if (u == Point3T<T>() || v == Point3T<T>()) return 0;  // atan2(0, -0) = π にならないように
    const Point3D a = to_d(u), b = to_d(v);
    return std::atan2(len(cross(a, b)), dot(a, b));
}
template <class T>
Real angle_between_deg(const Point3T<T>& u, const Point3T<T>& v) { return rad_to_deg(angle_between(u, v)); }
template <class T>
Real elevation_angle(const Point3T<T>& v) { return std::atan2((Real)v.z, std::hypot((Real)v.x, (Real)v.y)); }  // 仰角 [-π/2, π/2]
template <class T>
Real elevation_angle_deg(const Point3T<T>& v) { return rad_to_deg(elevation_angle(v)); }

// ---- 回転（右ねじの向き：軸の先から見て反時計回り） ----
template <class T>
constexpr Point3T<T> rot90_x(const Point3T<T>& p) { return Point3T<T>(p.x, -p.z, p.y); }  // x 軸のまわりに 90 度（y → z）
template <class T>
constexpr Point3T<T> rot90_y(const Point3T<T>& p) { return Point3T<T>(p.z, p.y, -p.x); }  // y 軸のまわりに 90 度（z → x）
template <class T>
constexpr Point3T<T> rot90_z(const Point3T<T>& p) { return Point3T<T>(-p.y, p.x, p.z); }  // z 軸のまわりに 90 度（x → y）
// 立方体の 24 通りの回転の k 番目（0 ≤ k < 24。k = 0 はそのまま）。x 軸のまわりに k % 4 回まわしてから、x 軸の正の向きを 6 方向の k / 4 番目へ
template <class T>
constexpr Point3T<T> rotation24(Point3T<T> p, int k) {
    for (int i = 0; i < k % 4; ++i) p = rot90_x(p);
    const int f = k / 4;
    if (f == 1) p = rot90_z(p);                              // x → y
    else if (f == 2) p = rot90_z(rot90_z(p));                // x → -x
    else if (f == 3) p = rot90_z(rot90_z(rot90_z(p)));       // x → -y
    else if (f == 4) p = rot90_y(rot90_y(rot90_y(p)));       // x → z
    else if (f == 5) p = rot90_y(p);                         // x → -z
    return p;
}
// 原点を通る軸 axis（0 でない）のまわりに theta（ラジアン）回す（ロドリゲスの回転公式）
// p・軸・中心の型は違ってよい（Point3 と Point3D など）。結果の -0 は 0 にする
namespace geo_internal {
inline Point3D fix_zero3(const Point3D& p) { return Point3D(p.x + (Real)0, p.y + (Real)0, p.z + (Real)0); }
}  // namespace geo_internal
template <class T, class U>
Point3D rotated(const Point3T<T>& p, Real theta, const Point3T<U>& axis) {
    const Point3D k = unit(axis), v = to_d(p);
    const Real c = std::cos(theta), s = std::sin(theta);
    return geo_internal::fix_zero3(v * c + cross(k, v) * s + k * (dot(k, v) * (1 - c)));
}
template <class T, class U, class V>
Point3D rotated(const Point3T<T>& p, Real theta, const Point3T<U>& axis, const Point3T<V>& center) {
    return geo_internal::fix_zero3(rotated(to_d(p) - to_d(center), theta, axis) + to_d(center));
}
template <class T, class U>
Point3D rotated_deg(const Point3T<T>& p, Real deg, const Point3T<U>& axis) {  // 度で。座標軸のまわりの 90 度の倍数なら誤差なし
    const Real q = deg / 90;
    const int nz = (axis.x != 0) + (axis.y != 0) + (axis.z != 0);
    if (nz == 1 && q == std::floor(q) && std::fabs(q) < 1e18L) {
        long long k = ((long long)q % 4 + 4) % 4;
        if (axis.x + axis.y + axis.z < 0) k = (4 - k) % 4;  // 軸が負の向き → 逆回り
        Point3T<T> r = p;
        for (long long i = 0; i < k; ++i) r = axis.x != 0 ? rot90_x(r) : (axis.y != 0 ? rot90_y(r) : rot90_z(r));
        return geo_internal::fix_zero3(to_d(r));
    }
    return rotated(p, deg_to_rad(deg), axis);
}
template <class T, class U, class V>
Point3D rotated_deg(const Point3T<T>& p, Real deg, const Point3T<U>& axis, const Point3T<V>& center) {
    return geo_internal::fix_zero3(rotated_deg(to_d(p) - to_d(center), deg, axis) + to_d(center));
}

// ---- 直線・線分・平面・球の型 ----
template <class T>
struct Line3T {  // 直線 a → b（a ≠ b）
    Point3T<T> a, b;
    constexpr Line3T() {}
    constexpr Line3T(const Point3T<T>& a_, const Point3T<T>& b_) : a(a_), b(b_) {}
    constexpr Point3T<T> dir() const { return b - a; }
    friend constexpr bool operator==(const Line3T& l, const Line3T& m) { return l.a == m.a && l.b == m.b; }
    friend constexpr bool operator!=(const Line3T& l, const Line3T& m) { return !(l == m); }
};
template <class T>
struct Segment3T {  // 線分 a - b（a == b の点も可）
    Point3T<T> a, b;
    constexpr Segment3T() {}
    constexpr Segment3T(const Point3T<T>& a_, const Point3T<T>& b_) : a(a_), b(b_) {}
    constexpr Point3T<T> dir() const { return b - a; }
    friend constexpr bool operator==(const Segment3T& s, const Segment3T& t) { return s.a == t.a && s.b == t.b; }
    friend constexpr bool operator!=(const Segment3T& s, const Segment3T& t) { return !(s == t); }
};
template <class T>
struct Plane3T {  // 平面 n · x = d（法線 n = (a, b, c) ≠ 0。整数なら d は Int128）。表側（side = +1）は法線の向き
    Point3T<T> n;
    geo_internal::wide3_t<T> d = 0;
    constexpr Plane3T() {}
    constexpr Plane3T(const Point3T<T>& n_, geo_internal::wide3_t<T> d_) : n(n_), d(d_) {}
    template <class A, class B, class C, class D,
              std::enable_if_t<geo_internal::coord_args_ok<T, A, B, C> && (std::is_arithmetic_v<D> || std::is_same_v<D, Int128>) &&
                                   (std::is_floating_point_v<T> || !std::is_floating_point_v<D>),
                               int> = 0>
    constexpr Plane3T(A a, B b, C c, D d_) : n(a, b, c), d(static_cast<geo_internal::wide3_t<T>>(d_)) {}  // a x + b y + c z = d
    friend constexpr bool operator==(const Plane3T& p, const Plane3T& q) { return p.n == q.n && p.d == q.d; }
    friend constexpr bool operator!=(const Plane3T& p, const Plane3T& q) { return !(p == q); }
    friend constexpr bool operator<(const Plane3T& p, const Plane3T& q) { return p.n != q.n ? p.n < q.n : p.d < q.d; }  // set / map 用
};
template <class T>
struct SphereT {  // 球（中心 c・半径 r）
    Point3T<T> c;
    T r = 0;
    constexpr SphereT() {}
    template <class R, std::enable_if_t<geo_internal::coord_args_ok<T, R>, int> = 0>
    constexpr SphereT(const Point3T<T>& c_, R r_) : c(c_), r(static_cast<T>(r_)) {}
    template <class A, class B, class C, class R, std::enable_if_t<geo_internal::coord_args_ok<T, A, B, C, R>, int> = 0>
    constexpr SphereT(A x, B y, C z, R r_) : c(x, y, z), r(static_cast<T>(r_)) {}
    friend constexpr bool operator==(const SphereT& p, const SphereT& q) { return p.c == q.c && p.r == q.r; }
    friend constexpr bool operator!=(const SphereT& p, const SphereT& q) { return !(p == q); }
};
using Line3 = Line3T<Coord>;
using Line3D = Line3T<Real>;
using Segment3 = Segment3T<Coord>;
using Segment3D = Segment3T<Real>;
using Plane3 = Plane3T<Coord>;
using Plane3D = Plane3T<Real>;
using Sphere = SphereT<Coord>;
using SphereD = SphereT<Real>;
template <class T>
constexpr Line3D to_d(const Line3T<T>& l) { return Line3D(to_d(l.a), to_d(l.b)); }
template <class T>
constexpr Segment3D to_d(const Segment3T<T>& s) { return Segment3D(to_d(s.a), to_d(s.b)); }
template <class T>
Plane3D to_d(const Plane3T<T>& p) { return Plane3D(to_d(p.n), (Real)p.d); }
template <class T>
constexpr SphereD to_d(const SphereT<T>& s) { return SphereD(to_d(s.c), s.r); }
template <class T>
std::ostream& operator<<(std::ostream& os, const Line3T<T>& l) { return os << "Line3[" << l.a << " -> " << l.b << "]"; }
template <class T>
std::ostream& operator<<(std::ostream& os, const Segment3T<T>& s) { return os << "Segment3[" << s.a << " - " << s.b << "]"; }
template <class T>
std::ostream& operator<<(std::ostream& os, const Plane3T<T>& p) {
    os << "Plane3[" << p.n.x << " x + " << p.n.y << " y + " << p.n.z << " z = ";
    if constexpr (std::is_integral_v<T>) os << int128_to_string(p.d);
    else os << p.d;
    return os << "]";
}
template <class T>
std::ostream& operator<<(std::ostream& os, const SphereT<T>& s) { return os << "Sphere[c=" << s.c << ", r=" << s.r << "]"; }
// 3 点 p, q, r（一直線でない）を通る平面。p → q → r が反時計回りに見える側が表（法線 = (q - p) × (r - p)）
template <class T>
Plane3T<T> plane_from_points(const Point3T<T>& p, const Point3T<T>& q, const Point3T<T>& r) {
    const Point3T<T> n = cross(q - p, r - p);
    return Plane3T<T>(n, geo_internal::dotw(n, p));
}

namespace geo_internal {
template <class X>
struct linear3_t : std::false_type {};
template <class T>
struct linear3_t<Line3T<T>> : std::true_type {};
template <class T>
struct linear3_t<Segment3T<T>> : std::true_type {};
template <class X>
constexpr bool is_linear3_v = linear3_t<X>::value;
template <class T>
Real scale3(const Point3T<T>& a, const Point3T<T>& b, const Point3T<T>& c, const Point3T<T>& d) {  // 実数の判定の誤差の目安（座標の大きさ）
    return std::max(std::max(len(a), len(b)), std::max(len(c), len(d)));
}
// 直線 a1 + s d1 と a2 + t d2（n = d1 × d2 ≠ 0）で最も近い点の s, t。s = s_num / den、t = t_num / den（den = |n|² > 0）
template <class T>
void closest_params(const Point3T<T>& a1, const Point3T<T>& d1, const Point3T<T>& a2, const Point3T<T>& d2, const Point3T<T>& n, wide3_t<T>& s_num,
                    wide3_t<T>& t_num, wide3_t<T>& den) {
    const Point3T<T> w = a2 - a1;
    s_num = dotw(cross(w, d2), n), t_num = dotw(cross(w, d1), n), den = dotw(n, n);
}
}  // namespace geo_internal

// ---- 直線・線分 ----
template <class L1, class L2, std::enable_if_t<geo_internal::is_linear3_v<L1> && geo_internal::is_linear3_v<L2>, int> = 0>
bool is_parallel(const L1& l, const L2& m) { return geo_internal::is_zero3(cross(l.dir(), m.dir()), len(l.dir()) * len(m.dir())); }
template <class L1, class L2, std::enable_if_t<geo_internal::is_linear3_v<L1> && geo_internal::is_linear3_v<L2>, int> = 0>
bool is_orthogonal(const L1& l, const L2& m) {
    const auto v = geo_internal::dotw(l.dir(), m.dir());
    if constexpr (std::is_integral_v<std::decay_t<decltype(l.a.x)>>) return v == 0;
    else return real_sign(v, len(l.dir()) * len(m.dir())) == 0;
}
template <class T>
bool on_line(const Line3T<T>& l, const Point3T<T>& p) { return collinear(l.a, l.b, p); }
template <class T>
Point3D projection(const Line3T<T>& l, const Point3T<T>& p) {  // 垂線の足
    const Point3T<T> d = l.dir();
    return to_d(l.a) + to_d(d) * ((Real)geo_internal::dotw(p - l.a, d) / (Real)geo_internal::dotw(d, d));
}
template <class T>
Point3D reflection(const Line3T<T>& l, const Point3T<T>& p) { return projection(l, p) * 2 - to_d(p); }  // 直線について対称な点
template <class T>
Real dist(const Line3T<T>& l, const Point3T<T>& p) {
    const Point3D d = to_d(l.dir());
    return len(cross(d, to_d(p - l.a))) / len(d);
}
template <class T>
Real dist(const Segment3T<T>& s, const Point3T<T>& p) {
    const Point3T<T> d = s.dir();
    if (geo_internal::dotw(p - s.a, d) <= 0) return dist(p, s.a);
    if (geo_internal::dotw(p - s.b, d) >= 0) return dist(p, s.b);
    return dist(Line3T<T>(s.a, s.b), p);
}
template <class T>
Real dist(const Point3T<T>& p, const Line3T<T>& l) { return dist(l, p); }
template <class T>
Real dist(const Point3T<T>& p, const Segment3T<T>& s) { return dist(s, p); }
template <class T>
bool on_segment(const Segment3T<T>& s, const Point3T<T>& p) {  // 線分の上（端点を含む）
    if constexpr (std::is_integral_v<T>) return collinear(s.a, s.b, p) && dot128(s.a - p, s.b - p) <= 0;
    else return real_sign(dist(s, p), geo_internal::scale3(s.a, s.b, p, p)) == 0;
}
template <class T>
Real dist(const Line3T<T>& l, const Line3T<T>& m) {  // 直線どうし（ねじれの位置でも）
    const Point3T<T> u = l.dir(), v = m.dir(), n = cross(u, v);
    if (geo_internal::is_zero3(n, len(u) * len(v))) return dist(l, m.a);
    return std::fabs((Real)geo_internal::dotw(n, m.a - l.a)) / len(n);
}
template <class T>
Real dist(const Segment3T<T>& s, const Segment3T<T>& t) {  // 線分どうし
    const Point3T<T> d1 = s.dir(), d2 = t.dir(), n = cross(d1, d2);
    if (!geo_internal::is_zero3(n, len(d1) * len(d2))) {  // 最も近い点が両方の線分の内側なら、直線どうしの距離
        geo_internal::wide3_t<T> ps, pt, den;
        geo_internal::closest_params(s.a, d1, t.a, d2, n, ps, pt, den);
        if (0 <= ps && ps <= den && 0 <= pt && pt <= den) return std::fabs((Real)geo_internal::dotw(n, t.a - s.a)) / len(n);
    }
    return std::min(std::min(dist(s, t.a), dist(s, t.b)), std::min(dist(t, s.a), dist(t, s.b)));
}
template <class T>
Real dist(const Line3T<T>& l, const Segment3T<T>& s) {  // 直線と線分
    const Point3T<T> d1 = l.dir(), d2 = s.dir(), n = cross(d1, d2);
    if (!geo_internal::is_zero3(n, len(d1) * len(d2))) {
        geo_internal::wide3_t<T> ps, pt, den;
        geo_internal::closest_params(l.a, d1, s.a, d2, n, ps, pt, den);
        if (0 <= pt && pt <= den) return std::fabs((Real)geo_internal::dotw(n, s.a - l.a)) / len(n);
    }
    return std::min(dist(l, s.a), dist(l, s.b));
}
template <class T>
Real dist(const Segment3T<T>& s, const Line3T<T>& l) { return dist(l, s); }
// 共有点を持つか（整数なら厳密。線分は端点で触れる・重なるも含む）
template <class T>
bool intersects(const Segment3T<T>& s, const Segment3T<T>& t) {
    if constexpr (std::is_integral_v<T>) {
        const Point3T<T> d1 = s.dir(), d2 = t.dir(), n = cross(d1, d2);
        if (d1 == Point3T<T>()) return on_segment(t, s.a);
        if (d2 == Point3T<T>()) return on_segment(s, t.a);
        if (n == Point3T<T>()) {  // 平行：同じ直線の上で重なるか
            if (!collinear(s.a, s.b, t.a)) return false;
            return on_segment(s, t.a) || on_segment(s, t.b) || on_segment(t, s.a) || on_segment(t, s.b);
        }
        if (dot128(n, t.a - s.a) != 0) return false;  // 同じ平面にない（ねじれ）
        Int128 ps, pt, den;
        geo_internal::closest_params(s.a, d1, t.a, d2, n, ps, pt, den);
        return 0 <= ps && ps <= den && 0 <= pt && pt <= den;
    } else {
        return real_sign(dist(s, t), geo_internal::scale3(s.a, s.b, t.a, t.b)) == 0;
    }
}
template <class T>
bool intersects(const Line3T<T>& l, const Line3T<T>& m) {
    if constexpr (std::is_integral_v<T>) {
        const Point3T<T> n = cross(l.dir(), m.dir());
        if (n == Point3T<T>()) return collinear(l.a, l.b, m.a);  // 平行：同じ直線か
        return dot128(n, m.a - l.a) == 0;                          // 同じ平面にあれば交わる
    } else {
        return real_sign(dist(l, m), geo_internal::scale3(l.a, l.b, m.a, m.b)) == 0;
    }
}
template <class T>
bool intersects(const Line3T<T>& l, const Segment3T<T>& s) {
    if constexpr (std::is_integral_v<T>) {
        const Point3T<T> d1 = l.dir(), d2 = s.dir(), n = cross(d1, d2);
        if (n == Point3T<T>()) return on_line(l, s.a);  // 平行（線分が点のときも）
        if (dot128(n, s.a - l.a) != 0) return false;
        Int128 ps, pt, den;
        geo_internal::closest_params(l.a, d1, s.a, d2, n, ps, pt, den);
        return 0 <= pt && pt <= den;
    } else {
        return real_sign(dist(l, s), geo_internal::scale3(l.a, l.b, s.a, s.b)) == 0;
    }
}
template <class T>
bool intersects(const Segment3T<T>& s, const Line3T<T>& l) { return intersects(l, s); }
// 2 直線の位置関係：交わる / 平行 / 同じ直線 / ねじれの位置（整数なら厳密）
enum Line3Relation : int { LINE3_CROSS = 0, LINE3_PARALLEL = 1, LINE3_SAME = 2, LINE3_SKEW = 3 };
template <class T>
int line_relation(const Line3T<T>& l, const Line3T<T>& m) {
    if (is_parallel(l, m)) return on_line(l, m.a) ? LINE3_SAME : LINE3_PARALLEL;
    return intersects(l, m) ? LINE3_CROSS : LINE3_SKEW;
}
// 交点（実数）。直線どうしは 1 点で交わるときだけ（平行・同じ直線・ねじれなら空）
template <class T>
std::vector<Point3D> cross_point(const Line3T<T>& l, const Line3T<T>& m) {
    if (line_relation(l, m) != LINE3_CROSS) return {};
    const Point3T<T> d1 = l.dir(), d2 = m.dir();
    geo_internal::wide3_t<T> ps, pt, den;
    geo_internal::closest_params(l.a, d1, m.a, d2, cross(d1, d2), ps, pt, den);
    return {to_d(l.a) + to_d(d1) * ((Real)ps / (Real)den)};
}
// 直線と線分（線分が直線の上なら両端）
template <class T>
std::vector<Point3D> cross_point(const Line3T<T>& l, const Segment3T<T>& s) {
    if (!intersects(l, s)) return {};
    if (s.a == s.b) return {to_d(s.a)};
    if (is_parallel(l, s)) return {to_d(s.a), to_d(s.b)};
    return cross_point(l, Line3T<T>(s.a, s.b));
}
template <class T>
std::vector<Point3D> cross_point(const Segment3T<T>& s, const Line3T<T>& l) { return cross_point(l, s); }
// 線分どうしの交点（一直線で重なるときは重なりの両端を s の向きの順に。1 点なら 1 個）
template <class T>
std::vector<Point3D> cross_points(const Segment3T<T>& s, const Segment3T<T>& t) {
    if (!intersects(s, t)) return {};
    if (s.a == s.b) return {to_d(s.a)};
    if (t.a == t.b) return {to_d(t.a)};
    if (!is_parallel(s, t)) return cross_point(Line3T<T>(s.a, s.b), Line3T<T>(t.a, t.b));
    const Point3T<T> d = s.dir();
    auto key = [&](const Point3T<T>& p) { return geo_internal::dotw(d, p - s.a); };  // s の向きに沿った位置
    Point3T<T> lo = t.a, hi = t.b;
    if (key(lo) > key(hi)) std::swap(lo, hi);
    const Point3T<T> L = key(lo) > key(s.a) ? lo : s.a, R = key(hi) < key(s.b) ? hi : s.b;
    bool same;
    if constexpr (std::is_integral_v<T>) same = (L == R);
    else same = real_sign(len(R - L), geo_internal::scale3(s.a, s.b, t.a, t.b)) == 0;
    if (same) return {to_d(L)};
    return {to_d(L), to_d(R)};
}

// 有理数の点 (x / d, y / d, z / d)（d > 0、既約）。同じ交点かどうかを誤差なしで判定する
struct FracPoint3 {
    Int128 x = 0, y = 0, z = 0, d = 1;
    Point3D point_d() const { return Point3D((Real)x / (Real)d, (Real)y / (Real)d, (Real)z / (Real)d); }
    bool is_integer() const { return d == 1; }
    friend bool operator==(const FracPoint3& p, const FracPoint3& q) { return p.x == q.x && p.y == q.y && p.z == q.z && p.d == q.d; }
    friend bool operator!=(const FracPoint3& p, const FracPoint3& q) { return !(p == q); }
    friend bool operator<(const FracPoint3& p, const FracPoint3& q) {  // set / map 用の順
        if (p.x != q.x) return p.x < q.x;
        if (p.y != q.y) return p.y < q.y;
        if (p.z != q.z) return p.z < q.z;
        return p.d < q.d;
    }
};
namespace geo_internal {
inline Int128 abs128_3(Int128 v) { return v < 0 ? -v : v; }
inline Int128 gcd128_3(Int128 a, Int128 b) {
    a = abs128_3(a), b = abs128_3(b);
    while (b != 0) {
        const Int128 t = a % b;
        a = b, b = t;
    }
    return a;
}
}  // namespace geo_internal
// 2 直線の交点を既約分数で（1 点で交わるときだけ。それ以外は空）
inline std::vector<FracPoint3> cross_point_frac(const Line3& l, const Line3& m) {
    if (line_relation(l, m) != LINE3_CROSS) return {};
    const Point3 d1 = l.dir(), d2 = m.dir(), n = cross(d1, d2), w = m.a - l.a;
    // n の 0 でない成分の座標を消した 2 次元で解く：s = cross2(w, d2) / cross2(d1, d2)（cross2(d1, d2) は n のその成分）
    Int128 den, num;
    if (n.z != 0) den = n.z, num = (Int128)w.x * d2.y - (Int128)w.y * d2.x;
    else if (n.x != 0) den = n.x, num = (Int128)w.y * d2.z - (Int128)w.z * d2.y;
    else den = n.y, num = (Int128)w.z * d2.x - (Int128)w.x * d2.z;
    if (den < 0) den = -den, num = -num;
    const Int128 x = (Int128)l.a.x * den + (Int128)d1.x * num, y = (Int128)l.a.y * den + (Int128)d1.y * num, z = (Int128)l.a.z * den + (Int128)d1.z * num;
    const Int128 g = geo_internal::gcd128_3(geo_internal::gcd128_3(geo_internal::gcd128_3(x, y), z), den);
    return {FracPoint3{x / g, y / g, z / g, den / g}};
}
// 向きの正規化（gcd で割る）。direction_key は向きを区別する、line_direction_key は v と -v を同じに（最初の 0 でない成分を正）
inline Point3 direction_key(const Point3& v) {
    const Coord g = std::gcd(std::gcd(v.x < 0 ? -v.x : v.x, v.y < 0 ? -v.y : v.y), v.z < 0 ? -v.z : v.z);
    return g == 0 ? Point3() : Point3(v.x / g, v.y / g, v.z / g);
}
inline Point3 line_direction_key(const Point3& v) {
    Point3 d = direction_key(v);
    if (d.x < 0 || (d.x == 0 && (d.y < 0 || (d.y == 0 && d.z < 0)))) d = -d;
    return d;
}
// 3 次元の直線のキー（同じ直線なら同じ）。d = line_direction_key(向き)、m = p × d（直線の上のどの点 p でも同じ値）
struct LineKey3 {
    Point3 d, m;
    friend bool operator==(const LineKey3& a, const LineKey3& b) { return a.d == b.d && a.m == b.m; }
    friend bool operator!=(const LineKey3& a, const LineKey3& b) { return !(a == b); }
    friend bool operator<(const LineKey3& a, const LineKey3& b) { return a.d != b.d ? a.d < b.d : a.m < b.m; }
};
inline LineKey3 line_key(const Point3& p, const Point3& q) {  // p と q（p ≠ q）を通る直線
    const Point3 d = line_direction_key(q - p);
    return LineKey3{d, cross(p, d)};
}
inline LineKey3 line_key(const Line3& l) { return line_key(l.a, l.b); }

// ---- 平面 ----
template <class T>
int side(const Plane3T<T>& pl, const Point3T<T>& p) {  // 表側 +1 / 裏側 -1 / 平面の上 0
    const auto v = geo_internal::dotw(pl.n, p) - pl.d;
    if constexpr (std::is_integral_v<T>) return (v > 0) - (v < 0);
    else return real_sign(v, len(pl.n) * len(p) + std::fabs(pl.d));
}
template <class T>
bool on_plane(const Plane3T<T>& pl, const Point3T<T>& p) { return side(pl, p) == 0; }
template <class T>
Real dist(const Plane3T<T>& pl, const Point3T<T>& p) { return std::fabs((Real)(geo_internal::dotw(pl.n, p) - pl.d)) / len(pl.n); }
template <class T>
Real dist(const Point3T<T>& p, const Plane3T<T>& pl) { return dist(pl, p); }
template <class T>
Point3D projection(const Plane3T<T>& pl, const Point3T<T>& p) {  // 垂線の足
    return to_d(p) - to_d(pl.n) * ((Real)(geo_internal::dotw(pl.n, p) - pl.d) / (Real)geo_internal::dotw(pl.n, pl.n));
}
template <class T>
Point3D reflection(const Plane3T<T>& pl, const Point3T<T>& p) { return projection(pl, p) * 2 - to_d(p); }  // 平面について対称な点
template <class T>
std::vector<Point3D> cross_point(const Line3T<T>& l, const Plane3T<T>& pl) {  // 直線と平面の交点（平行・平面の上なら空）
    const Point3T<T> d = l.dir();
    const auto den = geo_internal::dotw(pl.n, d);
    if constexpr (std::is_integral_v<T>) {
        if (den == 0) return {};
    } else {
        if (real_sign(den, len(pl.n) * len(d)) == 0) return {};
    }
    const auto num = pl.d - geo_internal::dotw(pl.n, l.a);
    return {to_d(l.a) + to_d(d) * ((Real)num / (Real)den)};
}
template <class T>
std::vector<Point3D> cross_point(const Plane3T<T>& pl, const Line3T<T>& l) { return cross_point(l, pl); }
template <class T>
std::vector<Point3D> cross_point(const Segment3T<T>& s, const Plane3T<T>& pl) {  // 線分と平面の交点（線分が平面の上なら両端）
    const int sa = side(pl, s.a), sb = side(pl, s.b);
    if (sa * sb > 0) return {};
    if (sa == 0 && sb == 0) return s.a == s.b ? std::vector<Point3D>{to_d(s.a)} : std::vector<Point3D>{to_d(s.a), to_d(s.b)};
    if (sa == 0) return {to_d(s.a)};
    if (sb == 0) return {to_d(s.b)};
    const auto va = geo_internal::dotw(pl.n, s.a) - pl.d, vb = geo_internal::dotw(pl.n, s.b) - pl.d;  // 符号が逆
    return {to_d(s.a) + to_d(s.dir()) * ((Real)va / (Real)(va - vb))};
}
template <class T>
std::vector<Point3D> cross_point(const Plane3T<T>& pl, const Segment3T<T>& s) { return cross_point(s, pl); }
// 直線と平面の位置関係：1 点で交わる / 平行（共有点なし）/ 直線が平面の上
enum LinePlaneRelation : int { LINE_PLANE_CROSS = 0, LINE_PLANE_PARALLEL = 1, LINE_PLANE_ON = 2 };
template <class T>
int line_plane_relation(const Line3T<T>& l, const Plane3T<T>& pl) {
    const auto den = geo_internal::dotw(pl.n, l.dir());
    bool par;
    if constexpr (std::is_integral_v<T>) par = (den == 0);
    else par = real_sign(den, len(pl.n) * len(l.dir())) == 0;
    if (!par) return LINE_PLANE_CROSS;
    return on_plane(pl, l.a) ? LINE_PLANE_ON : LINE_PLANE_PARALLEL;
}
template <class T>
bool intersects(const Line3T<T>& l, const Plane3T<T>& pl) { return line_plane_relation(l, pl) != LINE_PLANE_PARALLEL; }
template <class T>
bool intersects(const Segment3T<T>& s, const Plane3T<T>& pl) { return side(pl, s.a) * side(pl, s.b) <= 0; }
template <class T>
Real dist(const Line3T<T>& l, const Plane3T<T>& pl) { return line_plane_relation(l, pl) == LINE_PLANE_PARALLEL ? dist(pl, l.a) : 0; }
template <class T>
Real dist(const Segment3T<T>& s, const Plane3T<T>& pl) { return intersects(s, pl) ? 0 : std::min(dist(pl, s.a), dist(pl, s.b)); }
// 平面のキー：(a, b, c, d) を gcd で割り、(a, b, c) の最初の 0 でない成分を正にしたもの（同じ平面なら同じ。整数だけ）
inline Plane3 plane_key(const Plane3& pl) {
    using geo_internal::gcd128_3;
    Int128 g = gcd128_3(gcd128_3(gcd128_3(pl.n.x, pl.n.y), pl.n.z), pl.d);
    const Point3& n = pl.n;
    if (n.x < 0 || (n.x == 0 && (n.y < 0 || (n.y == 0 && n.z < 0)))) g = -g;
    return Plane3(Point3((Coord)(n.x / g), (Coord)(n.y / g), (Coord)(n.z / g)), pl.d / g);
}
// 2 平面の位置関係：交わる（交線がある）/ 平行 / 同じ平面（整数なら厳密）
enum PlaneRelation : int { PLANE_CROSS = 0, PLANE_PARALLEL = 1, PLANE_SAME = 2 };
template <class T>
int plane_relation(const Plane3T<T>& p, const Plane3T<T>& q) {
    if constexpr (std::is_integral_v<T>) {
        const Plane3 a = plane_key(p), b = plane_key(q);
        if (line_direction_key(a.n) != line_direction_key(b.n)) return PLANE_CROSS;
        return a == b ? PLANE_SAME : PLANE_PARALLEL;
    } else {
        if (real_sign(len(cross(p.n, q.n)), len(p.n) * len(q.n)) != 0) return PLANE_CROSS;
        const Point3D x0 = p.n * (p.d / norm2(p.n));  // p の上の点
        return on_plane(q, x0) ? PLANE_SAME : PLANE_PARALLEL;
    }
}
template <class T>
Real dist(const Plane3T<T>& p, const Plane3T<T>& q) {  // 平行なら間の距離、交わるなら 0
    if (plane_relation(p, q) != PLANE_PARALLEL) return 0;
    const Point3D np = to_d(p.n), nq = to_d(q.n);
    const Real a = (Real)p.d / len(np), b = (Real)q.d / len(nq);  // 原点から各平面までの符号付き距離（法線の向きで）
    return std::fabs(dot(np, nq) > 0 ? a - b : a + b);
}
namespace geo_internal {
// a·b - c·e（a, c は Int128、b, e は long long。積は 192 ビットまで）を誤差なく求めてから Real に直す
inline Real mul_sub_real(Int128 a, long long b, Int128 c, long long e) {
    const UInt128 mask = ((UInt128)1 << 64) - 1;
    auto split = [&](Int128 x, long long y, Int128& H, Int128& L) {  // x·y = H·2^64 + L（0 ≤ L < 2^64）
        const Int128 xl = (Int128)((UInt128)x & mask), P = xl * y;  // 0 ≤ xl < 2^64 なので |P| < 2^127
        H = (x >> 64) * y + (P >> 64), L = (Int128)((UInt128)P & mask);
    };
    Int128 H1, L1, H2, L2;
    split(a, b, H1, L1), split(c, e, H2, L2);
    Int128 H = H1 - H2, L = L1 - L2;  // -2^64 < L < 2^64
    // |L| < 2^63 にそろえる（H·2^64 と L が打ち消し合って、Real に直すときの丸めが目立たないように）
    const Int128 half = (Int128)1 << 63;
    if (L >= half) L -= (Int128)1 << 64, ++H;
    else if (L < -half) L += (Int128)1 << 64, --H;
    return std::ldexp((Real)H, 64) + (Real)L;
}
}  // namespace geo_internal
// 2 平面の交線（平行・同じ平面なら空）。a は交線の上で原点に最も近い点、向き b - a は p.n × q.n
// a = (w × u) / |u|²（u = p.n × q.n、w = p.d q.n - q.d p.n）。整数なら u と w を誤差なく求めてから割る（ほぼ平行な平面でも正しく）
template <class T>
std::vector<Line3D> cross_line(const Plane3T<T>& p, const Plane3T<T>& q) {
    Point3D u, w;
    if constexpr (std::is_integral_v<T>) {
        const Int128 cx = (Int128)p.n.y * q.n.z - (Int128)p.n.z * q.n.y, cy = (Int128)p.n.z * q.n.x - (Int128)p.n.x * q.n.z,
                     cz = (Int128)p.n.x * q.n.y - (Int128)p.n.y * q.n.x;
        if (cx == 0 && cy == 0 && cz == 0) return {};
        u = Point3D((Real)cx, (Real)cy, (Real)cz);
        w = Point3D(geo_internal::mul_sub_real(p.d, q.n.x, q.d, p.n.x), geo_internal::mul_sub_real(p.d, q.n.y, q.d, p.n.y),
                    geo_internal::mul_sub_real(p.d, q.n.z, q.d, p.n.z));
    } else {
        u = cross(p.n, q.n);
        if (real_sign(len(u), len(p.n) * len(q.n)) == 0) return {};
        w = q.n * p.d - p.n * q.d;
    }
    const Point3D x0 = cross(w, u) / norm2(u);
    return {Line3D(x0, x0 + u)};
}
template <class T>
Real angle_between(const Line3T<T>& l, const Line3T<T>& m) { return geo_internal::fold_half_pi(angle_between(l.dir(), m.dir())); }  // [0, π/2]
template <class T>
Real angle_between(const Plane3T<T>& p, const Plane3T<T>& q) { return geo_internal::fold_half_pi(angle_between(p.n, q.n)); }  // 2 面角 [0, π/2]
template <class T>
Real angle_between(const Line3T<T>& l, const Plane3T<T>& p) { return GEO_PI / 2 - geo_internal::fold_half_pi(angle_between(l.dir(), p.n)); }
template <class T>
Real angle_between(const Plane3T<T>& p, const Line3T<T>& l) { return angle_between(l, p); }

// ---- 球 ----
enum SphereRelation : int { SPHERE_CONTAIN = 0, SPHERE_INTERNAL = 1, SPHERE_INTERSECT = 2, SPHERE_EXTERNAL = 3, SPHERE_SEPARATE = 4, SPHERE_SAME = 5 };
template <class T>
Where contains(const SphereT<T>& s, const Point3T<T>& p) {  // 球の内外（球面の上は ON_BORDER）
    if constexpr (std::is_integral_v<T>) {
        const Int128 d2 = dist2_128(p, s.c), r2 = (Int128)s.r * s.r;
        return d2 < r2 ? INSIDE : (d2 == r2 ? ON_BORDER : OUTSIDE);
    } else {
        const int c = real_cmp(dist(s.c, p), s.r);
        return c < 0 ? INSIDE : (c == 0 ? ON_BORDER : OUTSIDE);
    }
}
// 2 つの球の位置関係（05_circle の CIRCLE_* と同じ値）：離れている 4 / 外接 3 / 交わる 2 / 内接 1 / 内側 0 / 同じ球 5
template <class T>
int sphere_relation(const SphereT<T>& a, const SphereT<T>& b) {
    if constexpr (std::is_integral_v<T>) {
        const Int128 d2 = dist2_128(a.c, b.c), s = (Int128)a.r + b.r, df = a.r > b.r ? (Int128)a.r - b.r : (Int128)b.r - a.r;
        if (d2 == 0 && a.r == b.r) return SPHERE_SAME;
        if (d2 > s * s) return SPHERE_SEPARATE;
        if (d2 == s * s) return SPHERE_EXTERNAL;
        if (d2 > df * df) return SPHERE_INTERSECT;
        if (d2 == df * df) return SPHERE_INTERNAL;
        return SPHERE_CONTAIN;
    } else {
        const Real d = dist(a.c, b.c), s = a.r + b.r, df = std::fabs(a.r - b.r), sc = std::max(d, s);
        if (real_sign(d, sc) == 0 && real_sign(a.r - b.r, sc) == 0) return SPHERE_SAME;
        const int c1 = real_sign(d - s, sc);
        if (c1 > 0) return SPHERE_SEPARATE;
        if (c1 == 0) return SPHERE_EXTERNAL;
        const int c2 = real_sign(d - df, sc);
        if (c2 > 0) return SPHERE_INTERSECT;
        if (c2 == 0) return SPHERE_INTERNAL;
        return SPHERE_CONTAIN;
    }
}
// 球面と直線の交点（l の向きの順。接するなら 1 個、交わらなければ空）。整数なら個数は厳密
template <class T>
std::vector<Point3D> cross_points(const SphereT<T>& s, const Line3T<T>& l) {
    const Point3T<T> d = l.dir(), w = l.a - s.c;
    Real A, B, C, sq = 0;
    int sD;
    if constexpr (std::is_integral_v<T>) {
        const Int128 A1 = norm2_128(d), B1 = dot128(d, w), C1 = norm2_128(w) - (Int128)s.r * s.r, D1 = B1 * B1 - A1 * C1;
        A = (Real)A1, B = (Real)B1, C = (Real)C1, sD = (D1 > 0) - (D1 < 0);
        if (D1 > 0) sq = std::sqrt((Real)D1);
    } else {
        const Real nw = norm2(w), r2 = s.r * s.r;
        A = norm2(d), B = dot(d, w), C = nw - r2;
        const Real disc = B * B - A * C;
        sD = real_sign(disc, std::max(B * B, A * std::max(nw, r2)));
        if (sD > 0) sq = std::sqrt(disc);
    }
    const Point3D a = to_d(l.a), dd = to_d(d);
    if (sD < 0) return {};
    if (sD == 0) return {a + dd * (-B / A)};
    const Real q = -(B + (B >= 0 ? sq : -sq));  // 桁落ちを避ける：t1 = q / A、t2 = C / q
    Real t1 = q / A, t2 = C / q;
    if (t1 > t2) std::swap(t1, t2);
    return {a + dd * t1, a + dd * t2};
}
inline Real sphere_cap_volume(Real r, Real h) { return GEO_PI * h * h * (3 * r - h) / 3; }  // 半径 r の球の、高さ h（0 ≤ h ≤ 2r）の部分の体積
// 2 つの球の共通部分の体積（2 つの球冠の和。接する直前でも桁落ちしにくい式で）
template <class T>
Real intersection_volume(const SphereT<T>& a, const SphereT<T>& b) {
    const int rel = sphere_relation(a, b);
    if (rel == SPHERE_SEPARATE || rel == SPHERE_EXTERNAL) return 0;
    const Real r1 = a.r, r2 = b.r;
    if (rel != SPHERE_INTERSECT) {
        const Real m = std::min(r1, r2);
        return 4 * GEO_PI * m * m * m / 3;
    }
    // d = 中心の距離。P = r1 + r2 - d、Q1 = d + r2 - r1、Q2 = d + r1 - r2（引き算で桁落ちする方は 2 乗の差 / 和で）
    Real d, P, Q1, Q2;
    if constexpr (std::is_integral_v<T>) {
        const Int128 d2 = dist2_128(a.c, b.c), s = (Int128)a.r + b.r, df = (Int128)a.r - b.r;
        d = std::sqrt((Real)d2);
        P = (Real)(s * s - d2) / ((Real)s + d);
        const Real q = (Real)(d2 - df * df) / (d + std::fabs((Real)df));
        if (df >= 0) Q1 = q, Q2 = d + (Real)df;
        else Q1 = d - (Real)df, Q2 = q;
    } else {
        const Real d2 = dist2(a.c, b.c), s = r1 + r2, df = r1 - r2;
        d = std::sqrt(d2);
        P = (s * s - d2) / (s + d);
        const Real q = (d2 - df * df) / (d + std::fabs(df));
        if (df >= 0) Q1 = q, Q2 = d + df;
        else Q1 = d - df, Q2 = q;
    }
    const Real h1 = P * Q1 / (2 * d), h2 = P * Q2 / (2 * d);  // 球 1・球 2 の、相手の中にある球冠の高さ
    return sphere_cap_volume(r1, h1) + sphere_cap_volume(r2, h2);
}

// ---- 球と平面・球どうしの交わり（円）、中身の詰まった球（球体）としての位置関係 ----
struct Circle3D {  // 3 次元の円：中心 c、円のある平面の単位法線 n、半径 r
    Point3D c, n;
    Real r = 0;
};
namespace geo_internal {
inline void mul_u128(UInt128 a, UInt128 b, UInt128& hi, UInt128& lo) {  // a · b = hi · 2^128 + lo
    const UInt128 M = ((UInt128)1 << 64) - 1;
    const UInt128 a0 = a & M, a1 = a >> 64, b0 = b & M, b1 = b >> 64;
    const UInt128 p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    const UInt128 mid = (p00 >> 64) + (p01 & M) + (p10 & M);
    lo = (p00 & M) | (mid << 64);
    hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
}
// (球の中心と平面の距離) と r の比較（-1 / 0 / 1）。整数なら (n·c - d)² と r²|n|² を 256 ビットで比べる
template <class T>
int sphere_plane_cmp(const SphereT<T>& sp, const Plane3T<T>& pl) {
    if constexpr (std::is_integral_v<T>) {
        const Int128 v = dotw(pl.n, sp.c) - pl.d;
        const UInt128 av = v < 0 ? (UInt128)0 - (UInt128)v : (UInt128)v;
        const UInt128 N2 = (UInt128)((Int128)pl.n.x * pl.n.x) + (UInt128)((Int128)pl.n.y * pl.n.y) + (UInt128)((Int128)pl.n.z * pl.n.z);
        UInt128 h1, l1, h2, l2;
        mul_u128(av, av, h1, l1), mul_u128((UInt128)((Int128)sp.r * sp.r), N2, h2, l2);
        if (h1 != h2) return h1 < h2 ? -1 : 1;
        return l1 < l2 ? -1 : (l1 > l2 ? 1 : 0);
    } else {
        return real_cmp(dist(pl, sp.c), sp.r);
    }
}
template <class T>
int ball_cmp(const SphereT<T>& sp, const Segment3T<T>& s) {  // (球の中心と線分の距離) と r の比較
    if constexpr (std::is_integral_v<T>) {
        const Point3 d = s.dir(), w = sp.c - s.a;
        const Int128 t = dot128(d, w), L = norm2_128(d), r2 = (Int128)sp.r * sp.r;
        Int128 lhs, rhs;
        if (L == 0 || t <= 0) lhs = norm2_128(w), rhs = r2;
        else if (t >= L) lhs = dist2_128(sp.c, s.b), rhs = r2;
        else {
            const Int128 cx = (Int128)d.y * w.z - (Int128)d.z * w.y, cy = (Int128)d.z * w.x - (Int128)d.x * w.z, cz = (Int128)d.x * w.y - (Int128)d.y * w.x;
            lhs = cx * cx + cy * cy + cz * cz, rhs = r2 * L;  // |d × w|² ≤ |d|²|w|² < 1.5×10^38
        }
        return (lhs > rhs) - (lhs < rhs);
    } else {
        return real_cmp(dist(s, sp.c), sp.r);
    }
}
template <class T>
int ball_cmp(const SphereT<T>& sp, const Line3T<T>& l) {
    if constexpr (std::is_integral_v<T>) {
        const Point3 d = l.dir(), w = sp.c - l.a;
        const Int128 cx = (Int128)d.y * w.z - (Int128)d.z * w.y, cy = (Int128)d.z * w.x - (Int128)d.x * w.z, cz = (Int128)d.x * w.y - (Int128)d.y * w.x;
        const Int128 lhs = cx * cx + cy * cy + cz * cz, rhs = (Int128)sp.r * sp.r * norm2_128(d);
        return (lhs > rhs) - (lhs < rhs);
    } else {
        return real_cmp(dist(l, sp.c), sp.r);
    }
}
}  // namespace geo_internal
// 球面と平面：0 交わらない / 1 接する / 2 円で交わる（整数なら厳密）
template <class T>
int sphere_plane_relation(const SphereT<T>& sp, const Plane3T<T>& pl) { return 1 - geo_internal::sphere_plane_cmp(sp, pl); }
// 球面と平面の交わりの円（接するなら r = 0、交わらなければ空）
template <class T>
std::vector<Circle3D> cross_circle(const SphereT<T>& sp, const Plane3T<T>& pl) {
    const int k = sphere_plane_relation(sp, pl);
    if (k == 0) return {};
    const Point3D c = projection(pl, sp.c);
    Real r = 0;
    if (k == 2) {
        const Real e = dist(pl, sp.c);
        r = std::sqrt(std::max<Real>(0, ((Real)sp.r - e) * ((Real)sp.r + e)));
    }
    return {Circle3D{c, unit(pl.n), r}};
}
// 2 つの球面の交わりの円（接するなら r = 0。同じ球・交わらないなら空）。法線は a.c → b.c の向き
template <class T>
std::vector<Circle3D> cross_circle(const SphereT<T>& a, const SphereT<T>& b) {
    const int rel = sphere_relation(a, b);
    if (rel == SPHERE_SAME || rel == SPHERE_SEPARATE || rel == SPHERE_CONTAIN) return {};
    const Point3D u = to_d(b.c - a.c);
    Real t, r = 0;  // 中心 = a.c + u · t
    if constexpr (std::is_integral_v<T>) {
        const Int128 D2 = dist2_128(a.c, b.c), K = (Int128)a.r * a.r - (Int128)b.r * b.r + D2;  // K = 2 · d · (a.c から円の中心まで)
        t = (Real)K / (2 * (Real)D2);
        if (rel == SPHERE_INTERSECT) {
            const Int128 H = 4 * D2 * a.r * a.r - K * K;  // (2 · d · r)²
            r = std::sqrt((Real)(H < 0 ? 0 : H)) / (2 * std::sqrt((Real)D2));
        }
    } else {
        const Real D2 = norm2(u), K = a.r * a.r - b.r * b.r + D2;
        t = K / (2 * D2);
        if (rel == SPHERE_INTERSECT) r = std::sqrt(std::max<Real>(0, 4 * D2 * a.r * a.r - K * K)) / (2 * std::sqrt(D2));
    }
    return {Circle3D{to_d(a.c) + u * t, unit(u), r}};
}
// 球体（中身の詰まった球。球面を含む）としての位置関係・距離。整数なら判定は厳密
template <class T>
bool ball_intersects(const SphereT<T>& sp, const Segment3T<T>& s) { return geo_internal::ball_cmp(sp, s) <= 0; }
template <class T>
bool ball_intersects(const SphereT<T>& sp, const Line3T<T>& l) { return geo_internal::ball_cmp(sp, l) <= 0; }
template <class T>
bool ball_intersects(const SphereT<T>& sp, const Plane3T<T>& pl) { return geo_internal::sphere_plane_cmp(sp, pl) <= 0; }
template <class T>
bool ball_intersects(const SphereT<T>& a, const SphereT<T>& b) {  // 中心の距離 ≤ r1 + r2
    if constexpr (std::is_integral_v<T>) return dist2_128(a.c, b.c) <= ((Int128)a.r + b.r) * ((Int128)a.r + b.r);
    else return real_cmp(dist(a.c, b.c), a.r + b.r) <= 0;
}
template <class T>
bool ball_contains(const SphereT<T>& sp, const Segment3T<T>& s) { return contains(sp, s.a) != OUTSIDE && contains(sp, s.b) != OUTSIDE; }
template <class T>
bool ball_contains(const SphereT<T>& a, const SphereT<T>& b) {  // 球体 b が球体 a に入っているか（内接・同じ球も含む）
    if constexpr (std::is_integral_v<T>) {
        if (a.r < b.r) return false;
        const Int128 df = (Int128)a.r - b.r;
        return dist2_128(a.c, b.c) <= df * df;
    } else {
        return real_cmp(dist(a.c, b.c) + b.r, a.r) <= 0;
    }
}
template <class T>
Real ball_dist(const SphereT<T>& sp, const Point3T<T>& p) { return contains(sp, p) != OUTSIDE ? 0 : dist(sp.c, p) - (Real)sp.r; }
template <class T>
Real ball_dist(const SphereT<T>& sp, const Segment3T<T>& s) { return ball_intersects(sp, s) ? 0 : dist(s, sp.c) - (Real)sp.r; }
template <class T>
Real ball_dist(const SphereT<T>& sp, const Line3T<T>& l) { return ball_intersects(sp, l) ? 0 : dist(l, sp.c) - (Real)sp.r; }
template <class T>
Real ball_dist(const SphereT<T>& sp, const Plane3T<T>& pl) { return ball_intersects(sp, pl) ? 0 : dist(pl, sp.c) - (Real)sp.r; }
template <class T>
Real ball_dist(const SphereT<T>& a, const SphereT<T>& b) { return ball_intersects(a, b) ? 0 : dist(a.c, b.c) - (Real)a.r - (Real)b.r; }

// ---- 入力・ハッシュ ----
inline Point3 read_point3() { Point3 p; std::cin >> p.x >> p.y >> p.z; return p; }
inline Point3D read_point3_d() { Point3D p; std::cin >> p.x >> p.y >> p.z; return p; }
inline std::vector<Point3> read_points3(int n) {
    std::vector<Point3> ps(n);
    for (auto& p : ps) std::cin >> p.x >> p.y >> p.z;
    return ps;
}
inline std::vector<Point3D> read_points3_d(int n) {
    std::vector<Point3D> ps(n);
    for (auto& p : ps) std::cin >> p.x >> p.y >> p.z;
    return ps;
}
struct Point3Hash {  // unordered_set<Point3, Point3Hash> / unordered_map<Point3, V, Point3Hash>
    size_t operator()(const Point3& p) const {
        return (size_t)PointHash::mix(PointHash::mix(PointHash::mix((unsigned long long)p.x) ^ (unsigned long long)p.y) ^ (unsigned long long)p.z);
    }
};
struct LineKey3Hash {  // unordered_set<LineKey3, LineKey3Hash>
    size_t operator()(const LineKey3& k) const { return Point3Hash()(k.d) ^ (size_t)PointHash::mix((unsigned long long)Point3Hash()(k.m) + 0x9e37ULL); }
};
struct FracPoint3Hash {  // unordered_set<FracPoint3, FracPoint3Hash>
    size_t operator()(const FracPoint3& p) const {
        auto f = [](Int128 v) { return (unsigned long long)(Coord)(v >> 64) ^ (unsigned long long)(Coord)v; };
        return (size_t)PointHash::mix(PointHash::mix(PointHash::mix(f(p.x)) ^ f(p.y)) ^ f(p.z) ^ PointHash::mix(f(p.d)));
    }
};
#endif  // GG_GEO3D
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
namespace g3_test {

#define G3_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[geo3d] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK G3_CHECK

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

[[maybe_unused]] static Point3 rndp3(std::mt19937_64& rng, long long c) { return Point3(rnd(rng, -c, c), rnd(rng, -c, c), rnd(rng, -c, c)); }
[[maybe_unused]] static bool nearp3(const Point3D& a, const Point3D& b, Real tol) {
    const Real s = std::max<Real>(1, std::max({std::fabs(a.x), std::fabs(a.y), std::fabs(a.z), std::fabs(b.x), std::fabs(b.y), std::fabs(b.z)}));
    return std::fabs(a.x - b.x) <= tol * s && std::fabs(a.y - b.y) <= tol * s && std::fabs(a.z - b.z) <= tol * s;
}
// ---- 愚直な計算（128 ビット） ----
struct W3 {
    Int128 x, y, z;
};
static W3 w3(const Point3& p) { return W3{p.x, p.y, p.z}; }
static W3 wsub(const W3& a, const W3& b) { return W3{a.x - b.x, a.y - b.y, a.z - b.z}; }
static Int128 wdot(const W3& a, const W3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static W3 wcross(const W3& a, const W3& b) { return W3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
static bool wzero(const W3& a) { return a.x == 0 && a.y == 0 && a.z == 0; }
// 3×3 の行列式（余因子展開）
static Int128 det3(const W3& r0, const W3& r1, const W3& r2) {
    return r0.x * (r1.y * r2.z - r1.z * r2.y) - r0.y * (r1.x * r2.z - r1.z * r2.x) + r0.z * (r1.x * r2.y - r1.y * r2.x);
}
static Int128 b_vol6(const Point3& a, const Point3& b, const Point3& c, const Point3& d) { return det3(wsub(w3(b), w3(a)), wsub(w3(c), w3(a)), wsub(w3(d), w3(a))); }
static bool b_on_seg(const Point3& a, const Point3& b, const Point3& p) {  // p が線分 ab の上（a == b も）
    const W3 u = wsub(w3(b), w3(a)), v = wsub(w3(p), w3(a));
    if (!wzero(wcross(u, v))) return false;
    return wdot(wsub(w3(a), w3(p)), wsub(w3(b), w3(p))) <= 0;
}
// a1 + s d1 = a2 + t d2 を解く（2 つの成分の連立方程式をクラメルの公式で解き、残りの成分で確かめる）
// 戻り値：0 解なし、1 ただ 1 つ（s = sn / den, t = tn / den、den > 0）、2 平行（d1 × d2 = 0）
static int b_solve(const Point3& a1, const Point3& d1, const Point3& a2, const Point3& d2, Int128& sn, Int128& tn, Int128& den) {
    const Int128 D1[3] = {d1.x, d1.y, d1.z}, D2[3] = {d2.x, d2.y, d2.z}, Wv[3] = {(Int128)a2.x - a1.x, (Int128)a2.y - a1.y, (Int128)a2.z - a1.z};
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j) {
            // s D1[i] - t D2[i] = W[i]、s D1[j] - t D2[j] = W[j]
            const Int128 det = D1[i] * (-D2[j]) - (-D2[i]) * D1[j];
            if (det == 0) continue;
            Int128 s = Wv[i] * (-D2[j]) - (-D2[i]) * Wv[j], t = D1[i] * Wv[j] - Wv[i] * D1[j], dd = det;
            if (dd < 0) dd = -dd, s = -s, t = -t;
            const int k = 3 - i - j;
            if (s * D1[k] - t * D2[k] != Wv[k] * dd) return 0;
            sn = s, tn = t, den = dd;
            return 1;
        }
    return 2;
}
// 線分・直線の共有点（愚直）。kind：0 線分どうし、1 直線どうし、2 直線（1 つめ）と線分（2 つめ）
static bool b_meet(const Point3& a1, const Point3& b1, const Point3& a2, const Point3& b2, int kind) {
    const Point3 d1 = b1 - a1, d2 = b2 - a2;
    if (kind == 0 && d1 == Point3()) return b_on_seg(a2, b2, a1);
    if (kind != 1 && d2 == Point3()) {
        if (kind == 0) return b_on_seg(a1, b1, a2);
        return wzero(wcross(w3(d1), wsub(w3(a2), w3(a1))));
    }
    Int128 sn, tn, den;
    const int r = b_solve(a1, d1, a2, d2, sn, tn, den);
    if (r == 0) return false;
    if (r == 1) {
        const bool s_ok = kind != 0 || (0 <= sn && sn <= den), t_ok = kind == 1 || (0 <= tn && tn <= den);
        return s_ok && t_ok;
    }
    // 平行：同じ直線の上か
    if (!wzero(wcross(w3(d1), wsub(w3(a2), w3(a1))))) return false;
    if (kind != 0) return true;
    return b_on_seg(a1, b1, a2) || b_on_seg(a1, b1, b2) || b_on_seg(a2, b2, a1) || b_on_seg(a2, b2, b1);
}
// 直線どうしの距離の 2 乗（正規方程式を有理数で解く）。lines：直線どうし
static Real b_line_dist(const Point3& a1, const Point3& d1, const Point3& a2, const Point3& d2) {
    const Int128 A = norm2_128(d1), B = dot128(d1, d2), C = norm2_128(d2), D = dot128(d1, a1 - a2), E = dot128(d2, a1 - a2);
    const Int128 den = A * C - B * B;
    if (den == 0) {  // 平行：a2 から直線 1 への距離
        const Real t = -(Real)D / (Real)A;
        const Point3D p = to_d(a1) + to_d(d1) * t;
        return dist(p, to_d(a2));
    }
    // A s - B t = -D、B s - C t = -E（係数の行列式は B² - AC = -den）
    const Int128 sn = (-D) * (-C) - (-B) * (-E), tn = A * (-E) - (-D) * B;
    const Real s = -(Real)sn / (Real)den, t = -(Real)tn / (Real)den;
    return dist(to_d(a1) + to_d(d1) * s, to_d(a2) + to_d(d2) * t);
}
// 点と線分の距離（t を [0, 1] に切り詰める）
static Real b_seg_point(const Point3D& a, const Point3D& b, const Point3D& p) {
    const Point3D d = b - a;
    const Real dd = dot(d, d);
    Real t = dd == 0 ? 0 : dot(p - a, d) / dd;
    t = std::max<Real>(0, std::min<Real>(1, t));
    return dist(a + d * t, p);
}
// 凸関数の最小（三分探索）
template <class F>
static Real ternary_min(const F& f, Real lo, Real hi) {
    for (int it = 0; it < 200; ++it) {
        const Real m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
        if (f(m1) < f(m2)) hi = m2;
        else lo = m1;
    }
    return std::min({f(lo), f((lo + hi) / 2), f(hi)});
}
// 適応シンプソン
template <class F>
static Real simpson_rec(const F& f, Real a, Real b, Real fa, Real fm, Real fb, Real whole, Real eps, int depth) {
    const Real m = (a + b) / 2, lm = (a + m) / 2, rm = (m + b) / 2;
    const Real flm = f(lm), frm = f(rm);
    const Real left = (m - a) / 6 * (fa + 4 * flm + fm), right = (b - m) / 6 * (fm + 4 * frm + fb);
    if (depth <= 0 || std::fabs(left + right - whole) <= 15 * eps) return left + right + (left + right - whole) / 15;
    return simpson_rec(f, a, m, fa, flm, fm, left, eps / 2, depth - 1) + simpson_rec(f, m, b, fm, frm, fb, right, eps / 2, depth - 1);
}
template <class F>
static Real integrate(const F& f, Real a, Real b, Real eps) {
    Real s = 0;
    const int K = 32;
    for (int k = 0; k < K; ++k) {
        const Real x0 = a + (b - a) * k / K, x1 = a + (b - a) * (k + 1) / K, xm = (x0 + x1) / 2;
        const Real f0 = f(x0), fm = f(xm), f1 = f(x1);
        s += simpson_rec(f, x0, x1, f0, fm, f1, (x1 - x0) / 6 * (f0 + 4 * fm + f1), eps / K, 40);
    }
    return s;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 1000000000LL : (it % 3 == 1 ? 3 : 1000000);
        const Point3 a = rndp3(rng, C), b = rndp3(rng, C);
        // 演算子
        CHECK(a + b - b == a && -(-a) == a && a * 3 / 3 == a && 2 * a == a + a);
        Point3 t = a;
        t += b, t -= a;
        CHECK(t == b);
        CHECK((a < b) == (std::make_tuple(a.x, a.y, a.z) < std::make_tuple(b.x, b.y, b.z)) && (a <= b) == !(b < a) && (a == b) == !(a != b));
        // 128 ビットの版は愚直と同じ。差が 1.7×10^9 以下なら long long の版も同じ
        CHECK(dot128(a, b) == wdot(w3(a), w3(b)) && norm2_128(a - b) == wdot(wsub(w3(a), w3(b)), wsub(w3(a), w3(b))) && dist2_128(a, b) == norm2_128(a - b));
        const Point3 c = cross(a - b, b);  // 差が 2×10^9 まででも外積の成分は long long に収まる
        const W3 cw = wcross(wsub(w3(a), w3(b)), w3(b));
        CHECK(c.x == cw.x && c.y == cw.y && c.z == cw.z);
        if (C <= 1000000) CHECK(dot(a, b) == (Coord)dot128(a, b) && norm2(a) == (Coord)norm2_128(a) && dist2(a, b) == (Coord)dist2_128(a, b));
        CHECK(wdot(w3(c), wsub(w3(a), w3(b))) == 0 && wdot(w3(c), w3(b)) == 0);  // 外積は両方に垂直
        CHECK(near(len(a - b), std::sqrt((Real)dist2_128(a, b)), TOL) && near(dist(a, b), len(b - a), TOL));
        // 実数
        const Point3D ad = to_d(a), bd = to_d(b);
        CHECK(nearp3(midpoint(a, b), (ad + bd) / 2, TOL) && nearp3(Point3D(a), ad, 0));
        if (a != Point3()) CHECK(near(len(unit(a)), 1, TOL) && nearp3(unit(a) * len(a), ad, TOL));
        CHECK(unit(Point3()) == Point3D());
        CHECK(Point3Hash()(a) == Point3Hash()(Point3(a.x, a.y, a.z)));
    }
}

static void test_volume(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 1000000000LL : (it % 3 == 1 ? 1 : 4);
        Point3 a = rndp3(rng, C), b = rndp3(rng, C), c = rndp3(rng, C), d = rndp3(rng, C);
        if (it % 5 == 0 && C > 4) {  // 同じ平面・一直線の点を大きい座標で作る（a + s u + t v）
            const Point3 u = rndp3(rng, 1000), v = rndp3(rng, 1000);
            a = rndp3(rng, 800000000LL);
            b = a + u * rnd(rng, -100000, 100000), c = a + v * rnd(rng, -100000, 100000);
            d = a + u * rnd(rng, -100000, 100000) + v * rnd(rng, -100000, 100000);
            if (it % 10 == 0) d += Point3(rnd(rng, -1, 1), rnd(rng, -1, 1), rnd(rng, -1, 1));
            if (it % 15 == 0) c = a + u * rnd(rng, -100000, 100000);  // a, b, c が一直線
        }
        const Int128 v = b_vol6(a, b, c, d);
        CHECK(volume6(a, b, c, d) == v && volume6(a, c, b, d) == -v && volume6(b, c, d, a) == -v);
        CHECK(orient3d(a, b, c, d) == sgn128(v) && coplanar(a, b, c, d) == (v == 0));
        CHECK(near(tetra_volume(a, b, c, d), std::fabs((Real)v) / 6, TOL));
        const bool col = wzero(wcross(wsub(w3(b), w3(a)), wsub(w3(c), w3(a))));
        CHECK(collinear(a, b, c) == col && collinear(c, a, b) == col);
        CHECK(near(triangle_area(a, b, c), std::sqrt((Real)wdot(wcross(wsub(w3(b), w3(a)), wsub(w3(c), w3(a))), wcross(wsub(w3(b), w3(a)), wsub(w3(c), w3(a))))) / 2, TOL));
        // 表側：plane_from_points の side と orient3d は同じ
        if (!col) {
            const Plane3 pl = plane_from_points(a, b, c);
            CHECK(side(pl, d) == orient3d(a, b, c, d) && side(pl, a) == 0 && side(pl, b) == 0 && side(pl, c) == 0);
        }
        // 実数の点：小さい座標なら同じ。回転・拡大・平行移動しても同じ
        if (C <= 4) {
            const Point3D ad = to_d(a), bd = to_d(b), cd = to_d(c), dd = to_d(d);
            CHECK(orient3d(ad, bd, cd, dd) == sgn128(v) && collinear(ad, bd, cd) == col);
            const Point3D ax = Point3D(rndr(rng, -1, 1), rndr(rng, -1, 1), rndr(rng, 0.1L, 1));
            const Real th = rndr(rng, 0, 7), sc = std::pow((Real)10, rndr(rng, 0, 6));
            const Point3D sh(rndr(rng, -1e6L, 1e6L), rndr(rng, -1e6L, 1e6L), rndr(rng, -1e6L, 1e6L));
            auto f = [&](const Point3D& p) { return rotated(p, th, ax) * sc + sh; };
            CHECK(orient3d(f(ad), f(bd), f(cd), f(dd)) == sgn128(v) && collinear(f(ad), f(bd), f(cd)) == col);
        }
    }
    // 四面体の内外
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 2 : 1000000000LL / 4;
        Point3 a = rndp3(rng, C), b = rndp3(rng, C), c = rndp3(rng, C), d = rndp3(rng, C);
        Point3 p = rndp3(rng, C + 1);
        const int mode = it % 5;
        if (mode == 1) a *= 4, b *= 4, c *= 4, d *= 4, p = (a + b + c + d) / 4;  // 重心
        if (mode == 2) a *= 2, b *= 2, c *= 2, d *= 2, p = (b + c) / 2;          // 辺の中点
        if (mode == 3) a *= 3, b *= 3, c *= 3, d *= 3, p = (b + c + d) / 3;      // 面の重心
        if (mode == 4) p = a;                                                    // 頂点
        const Int128 V = b_vol6(a, b, c, d);
        const Where got = in_tetrahedron(a, b, c, d, p);
        if (V == 0) {
            CHECK(got == OUTSIDE);
        } else {
            // 重心座標の符号（愚直：行列式）
            const Int128 s0 = b_vol6(p, b, c, d), s1 = b_vol6(a, p, c, d), s2 = b_vol6(a, b, p, d), s3 = b_vol6(a, b, c, p);
            CHECK(s0 + s1 + s2 + s3 == V);
            const int o = sgn128(V);
            const int t[4] = {sgn128(s0) * o, sgn128(s1) * o, sgn128(s2) * o, sgn128(s3) * o};
            const bool out = t[0] < 0 || t[1] < 0 || t[2] < 0 || t[3] < 0, zero = t[0] == 0 || t[1] == 0 || t[2] == 0 || t[3] == 0;
            CHECK(got == (out ? OUTSIDE : (zero ? ON_BORDER : INSIDE)));
            CHECK(in_tetrahedron(b, a, c, d, p) == got && in_tetrahedron(a, b, c, d, a) == ON_BORDER);
            if (C <= 2) CHECK(in_tetrahedron(to_d(a), to_d(b), to_d(c), to_d(d), to_d(p)) == got);
        }
    }
}

static void test_angle_rotation(std::mt19937_64& rng, int& cases) {
    const Point3 ex(1, 0, 0), ey(0, 1, 0), ez(0, 0, 1);
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 1000000000LL : 5;
        const Point3 u = rndp3(rng, C), v = rndp3(rng, C);
        // なす角：カハンの式 2 atan2(| |v|u - |u|v |, | |v|u + |u|v |) と比べる
        const Point3D ud = to_d(u), vd = to_d(v);
        const Real lu = len(u), lv = len(v);
        if (lu > 0 && lv > 0) {
            const Real k = 2 * std::atan2(len(ud * lv - vd * lu), len(ud * lv + vd * lu));
            CHECK(std::fabs(angle_between(u, v) - k) <= 1e-12L);
            CHECK(near(angle_between_deg(u, v), rad_to_deg(k), 1e-12L));
            CHECK(std::fabs(elevation_angle(u) - (GEO_PI / 2 - angle_between(u, ez))) <= 1e-12L);
        }
        CHECK(angle_between(u, Point3()) == 0 && angle_between(Point3(), u) == 0 && angle_between(to_d(u), Point3D()) == 0 && elevation_angle(Point3()) == 0);
        CHECK(elevation_angle(Point3(-1, -1, 0)) == 0 && elevation_angle(Point3D(-1, -1, -0.0L)) <= 0);
        // 90 度の回転：4 回でもとに戻る・内積と外積を保つ（向きを変えない回転）
        CHECK(rot90_x(rot90_x(rot90_x(rot90_x(u)))) == u && rot90_y(rot90_y(rot90_y(rot90_y(u)))) == u && rot90_z(rot90_z(rot90_z(rot90_z(u)))) == u);
        if (C <= 5) {
            for (int k = 0; k < 24; ++k) {
                const Point3 ru = rotation24(u, k), rv = rotation24(v, k);
                CHECK(dot(ru, rv) == dot(u, v) && cross(ru, rv) == rotation24(cross(u, v), k));
            }
            CHECK(rot90_x(ey) == ez && rot90_y(ez) == ex && rot90_z(ex) == ey);
        }
        // 任意の軸のまわりの回転：長さ・軸方向の成分を保つ、逆に回すと戻る、z 軸なら 2 次元の回転と同じ、右ねじ
        Point3 ax = rndp3(rng, 5);
        if (ax == Point3()) ax = ez;
        const Real th = rndr(rng, -7, 7);
        const Point3D r = rotated(u, th, ax);
        CHECK(near(len(r), len(u), TOL));
        CHECK(std::fabs(dot(r, unit(ax)) - dot(ud, unit(ax))) <= TOL * std::max<Real>(1, lu));
        CHECK(nearp3(rotated(r, -th, to_d(ax)), ud, TOL));
        const PointD q2 = rotated(Point(u.x, u.y), th);
        CHECK(nearp3(rotated(u, th, ez), Point3D(q2.x, q2.y, (Real)u.z), TOL));
        CHECK(nearp3(rotated(ex, GEO_PI / 2, ez), to_d(ey), TOL) && nearp3(rotated(ey, GEO_PI / 2, ex), to_d(ez), TOL));
        const Point3 cen = rndp3(rng, C);
        CHECK(nearp3(rotated(u, th, ax, cen), rotated(u - cen, th, ax) + to_d(cen), TOL));
        // 度で：座標軸のまわりの 90 度の倍数は誤差なし
        const long long q = rnd(rng, -9, 9);
        for (int axis = 0; axis < 3; ++axis)
            for (int sg = -1; sg <= 1; sg += 2) {
                const Point3 A = axis == 0 ? ex * sg : (axis == 1 ? ey * sg : ez * sg);
                Point3 e = u, e2 = u - cen;
                const long long k = ((q * sg) % 4 + 4) % 4;
                for (long long i = 0; i < k; ++i) {
                    e = axis == 0 ? rot90_x(e) : (axis == 1 ? rot90_y(e) : rot90_z(e));
                    e2 = axis == 0 ? rot90_x(e2) : (axis == 1 ? rot90_y(e2) : rot90_z(e2));
                }
                CHECK(rotated_deg(u, (Real)(90 * q), A) == to_d(e));
                CHECK(rotated_deg(u, (Real)(90 * q), A * 7, cen) == to_d(e2 + cen));
            }
        const Real deg = rndr(rng, -400, 400);
        CHECK(nearp3(rotated_deg(u, deg, ax), rotated(u, deg_to_rad(deg), ax), TOL));
        CHECK(nearp3(rotated_deg(u, deg, ax, cen), rotated(u, deg_to_rad(deg), ax, cen), TOL));
    }
    // rotation24：24 通りが全部違う
    {
        std::vector<std::pair<Point3, Point3>> seen;
        for (int k = 0; k < 24; ++k) seen.push_back({rotation24(Point3(1, 0, 0), k), rotation24(Point3(0, 1, 0), k)});
        std::sort(seen.begin(), seen.end());
        CHECK(std::unique(seen.begin(), seen.end()) == seen.end() && rotation24(Point3(1, 2, 3), 0) == Point3(1, 2, 3));
        ++cases;
    }
}

static void test_lines(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 40000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 4 == 0 ? 1 : (it % 4 == 1 ? 2 : 5);
        Point3 a1 = rndp3(rng, C), b1 = rndp3(rng, C), a2 = rndp3(rng, C), b2 = rndp3(rng, C);
        if (it % 7 == 0) b2 = a2 + (b1 - a1) * rnd(rng, -2, 2);  // 平行
        if (it % 11 == 0) a2 = a1 + (b1 - a1) * rnd(rng, -2, 2);  // 同じ直線の上から
        const Segment3 s(a1, b1), t(a2, b2);
        // 交わるか（線分どうし・直線どうし・直線と線分）
        CHECK(intersects(s, t) == b_meet(a1, b1, a2, b2, 0) && intersects(t, s) == intersects(s, t));
        if (a1 != b1 && a2 != b2) {
            const Line3 l(a1, b1), m(a2, b2);
            CHECK(intersects(l, m) == b_meet(a1, b1, a2, b2, 1));
            CHECK(intersects(l, t) == b_meet(a1, b1, a2, b2, 2) && intersects(t, l) == intersects(l, t));
            CHECK(is_parallel(l, m) == wzero(wcross(w3(b1 - a1), w3(b2 - a2))) && is_parallel(s, m) == is_parallel(l, t));
            CHECK(is_orthogonal(l, m) == (dot128(b1 - a1, b2 - a2) == 0) && is_orthogonal(s, t) == is_orthogonal(l, m));
            // 実数の型でも同じ（小さい座標なので）
            CHECK(intersects(to_d(l), to_d(m)) == intersects(l, m) && intersects(to_d(l), to_d(t)) == intersects(l, t));
            CHECK(is_parallel(to_d(l), to_d(m)) == is_parallel(l, m) && is_orthogonal(to_d(l), to_d(t)) == is_orthogonal(l, m));
            // 距離：直線どうしは正規方程式、直線と線分は三分探索
            CHECK(std::fabs(dist(l, m) - b_line_dist(a1, b1 - a1, a2, b2 - a2)) <= 1e-12L);
            const Real e_ls = ternary_min([&](Real u) { return dist(to_d(l), to_d(a2) + to_d(b2 - a2) * u); }, 0, 1);
            CHECK(std::fabs(dist(l, t) - e_ls) <= 1e-9L && dist(t, l) == dist(l, t));
            CHECK(std::fabs(dist(to_d(l), to_d(m)) - dist(l, m)) <= 1e-12L && std::fabs(dist(to_d(l), to_d(t)) - dist(l, t)) <= 1e-12L);
            // 点：直線の上か・射影・対称な点・距離
            const Point3 p = rndp3(rng, C + 1);
            CHECK(on_line(l, p) == wzero(wcross(w3(b1 - a1), w3(p - a1))) && on_line(to_d(l), to_d(p)) == on_line(l, p));
            const Point3D pr = projection(l, p);
            CHECK(std::fabs(dot(to_d(p) - pr, to_d(b1 - a1))) <= 1e-12L && dist(to_d(l), pr) <= 1e-12L);
            CHECK(nearp3(midpoint(to_d(p), reflection(l, p)), pr, TOL));
            CHECK(std::fabs(dist(l, p) - dist(pr, to_d(p))) <= 1e-12L && dist(p, l) == dist(l, p));
        }
        // 点と線分
        const Point3 p = rndp3(rng, C + 1);
        CHECK(on_segment(s, p) == b_on_seg(a1, b1, p) && on_segment(to_d(s), to_d(p)) == on_segment(s, p));
        CHECK(std::fabs(dist(s, p) - b_seg_point(to_d(a1), to_d(b1), to_d(p))) <= 1e-12L && dist(p, s) == dist(s, p));
        // 線分どうしの距離：三分探索（s の上の点から t までの距離は凸）
        const Real e_ss = ternary_min([&](Real u) { return b_seg_point(to_d(a2), to_d(b2), to_d(a1) + to_d(b1 - a1) * u); }, 0, 1);
        CHECK(std::fabs(dist(s, t) - e_ss) <= 1e-9L && dist(t, s) == dist(s, t));
        CHECK((dist(s, t) == 0) == intersects(s, t));
        CHECK(std::fabs(dist(to_d(s), to_d(t)) - dist(s, t)) <= 1e-12L && intersects(to_d(s), to_d(t)) == intersects(s, t));
    }
    // 大きい座標：X で交わる線分を作る → 交わる。1 ずらす → 愚直と同じ
    for (int it = 0; it < 10000 * TEST_SCALE; ++it, ++cases) {
        const Point3 X = rndp3(rng, 900000000LL), u = rndp3(rng, 1000), v = rndp3(rng, 1000);
        Point3 a1 = X - u * rnd(rng, 0, 50000), b1 = X + u * rnd(rng, 0, 50000), a2 = X - v * rnd(rng, 0, 50000), b2 = X + v * rnd(rng, 0, 50000);
        CHECK(intersects(Segment3(a1, b1), Segment3(a2, b2)) && dist(Segment3(a1, b1), Segment3(a2, b2)) == 0);
        if (u != Point3() && v != Point3() && !wzero(wcross(w3(u), w3(v)))) {
            CHECK(intersects(Line3(a1, a1 + u), Line3(a2, a2 + v)));
            // 平面の上の 2 直線の交点 = X：平面 (a1, u, v) で切ってもよい
            const Plane3 pl = plane_from_points(X, X + u, X + v);
            CHECK(on_plane(pl, a1) && on_plane(pl, b2));
        }
        a2 += Point3(rnd(rng, -1, 1), rnd(rng, -1, 1), rnd(rng, -1, 1));
        CHECK(intersects(Segment3(a1, b1), Segment3(a2, b2)) == b_meet(a1, b1, a2, b2, 0));
        if (a1 != b1 && a2 != b2) {
            CHECK(intersects(Line3(a1, b1), Line3(a2, b2)) == b_meet(a1, b1, a2, b2, 1));
            CHECK(intersects(Line3(a1, b1), Segment3(a2, b2)) == b_meet(a1, b1, a2, b2, 2));
            // 距離（相対誤差）：直線どうしは正規方程式の値と
            const Real e = b_line_dist(a1, b1 - a1, a2, b2 - a2), g = dist(Line3(a1, b1), Line3(a2, b2));
            CHECK(std::fabs(g - e) <= 1e-9L * e + TOL * 1e9L);
        }
    }
    // 端の座標（10^9）でも 128 ビットに収まる
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long C = 1000000000LL;
        const Point3 a1 = rndp3(rng, C), b1 = rndp3(rng, C), a2 = rndp3(rng, C), b2 = rndp3(rng, C);
        const Segment3 s(a1, b1), t(a2, b2);
        CHECK(intersects(s, t) == b_meet(a1, b1, a2, b2, 0));
        const Real e = ternary_min([&](Real u) { return b_seg_point(to_d(a2), to_d(b2), to_d(a1) + to_d(b1 - a1) * u); }, 0, 1);
        CHECK(std::fabs(dist(s, t) - e) <= 1e-9L * e + TOL * 1e10L);
    }
}

static void test_planes(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 1000000000LL : (it % 3 == 1 ? 2 : 6);
        const Point3 p = rndp3(rng, C), q = rndp3(rng, C), r = rndp3(rng, C);
        if (collinear(p, q, r)) continue;
        const Plane3 pl = plane_from_points(p, q, r);
        const Point3 x = rndp3(rng, C);
        // 表裏・距離（体積 / 底面積）・射影・対称な点
        CHECK(side(pl, x) == sgn128(b_vol6(p, q, r, x)) && on_plane(pl, x) == (b_vol6(p, q, r, x) == 0));
        const Real e = std::fabs((Real)b_vol6(p, q, r, x)) / len(to_d(cross(q - p, r - p)));
        CHECK(near(dist(pl, x), e, 1e-15L) && dist(x, pl) == dist(pl, x));
        const Point3D pr = projection(pl, x), rf = reflection(pl, x);
        const Real sc0 = std::max<Real>(1, (Real)C);
        CHECK(dist(to_d(pl), pr) <= TOL * 16 * sc0 && nearp3((to_d(x) + rf) / 2, pr, TOL));
        if (dist(pl, x) > 1e-6L * sc0) CHECK(len(cross(unit(to_d(x) - pr), unit(to_d(pl.n)))) <= 1e-9L);
        // 方程式で作った平面
        const Point3 n = rndp3(rng, 3);
        if (n != Point3()) {
            const Coord d = rnd(rng, -10, 10);
            const Plane3 pe(n.x, n.y, n.z, d);
            CHECK(side(pe, x) == sgn128(dot128(n, x) - d));
            if (C <= 6) CHECK(side(to_d(pe), to_d(x)) == side(pe, x));
            CHECK(near(dist(pe, x), std::fabs((Real)(dot128(n, x) - d)) / len(n), 1e-15L));
        }
        // 直線・線分と平面の交点
        Point3 a = rndp3(rng, C), b = rndp3(rng, C);
        if (it % 5 == 0) b = a + (q - p) * rnd(rng, -1, 1);  // 平面に平行
        if (a == b) continue;
        const Line3 l(a, b);
        const Int128 den = dot128(pl.n, b - a), num = (Int128)0 + pl.d - dot128(pl.n, a);
        const std::vector<Point3D> cp = cross_point(l, pl);
        CHECK(cp.size() == (den != 0 ? 1u : 0u));
        if (den != 0) {
            const Point3D X = cp[0];
            const Real sc = std::max<Real>(std::max<Real>(1, (Real)C), len(X));
            CHECK(dist(to_d(pl), X) <= TOL * 64 * sc && dist(to_d(l), X) <= TOL * 64 * sc);
            CHECK(nearp3(X, to_d(a) + to_d(b - a) * ((Real)num / (Real)den), TOL));
            CHECK(cross_point(pl, l).size() == 1);
        }
        const Segment3 sg(a, b);
        const int sa = side(pl, a), sb = side(pl, b);
        const std::vector<Point3D> cs = cross_point(sg, pl);
        const size_t exp = sa * sb > 0 ? 0 : (sa == 0 && sb == 0 ? 2 : 1);
        CHECK(cs.size() == exp && cross_point(pl, sg).size() == exp);
        if (exp == 1 && sa != 0 && sb != 0) CHECK(nearp3(cs[0], cp[0], TOL * 16));
        // 実数の平面でも同じ（小さい座標）
        if (C <= 6) {
            CHECK(side(to_d(pl), to_d(x)) == side(pl, x) && cross_point(to_d(l), to_d(pl)).size() == cp.size() && cross_point(to_d(sg), to_d(pl)).size() == exp);
        }
    }
    // 2 平面の交線：同じ 2 点 P, Q を通る 2 平面（ほぼ平行なものも）の交線は P, Q を通る
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 1000000000LL : (it % 3 == 1 ? 4 : 1000);
        const Point3 P = rndp3(rng, C), Q = rndp3(rng, C), R1 = rndp3(rng, C);
        Point3 R2 = rndp3(rng, C);
        if (it % 2 == 0) {  // R2 を R1 のすぐ近くに（ほぼ平行な平面）
            R2 = R1 + Point3(rnd(rng, -1, 1), rnd(rng, -1, 1), rnd(rng, -1, 1));
            if (std::llabs(R2.x) > C || std::llabs(R2.y) > C || std::llabs(R2.z) > C) R2 = R1;
        }
        if (collinear(P, Q, R1) || collinear(P, Q, R2)) continue;
        const Plane3 A = plane_from_points(P, Q, R1), B = plane_from_points(P, Q, R2);
        const std::vector<Line3D> L = cross_line(A, B);
        const bool par = wzero(wcross(w3(A.n), w3(B.n)));
        CHECK(L.size() == (par ? 0u : 1u));
        if (L.empty()) continue;
        const Real sc = std::max({(Real)1, len(to_d(P)), len(to_d(Q))});
        CHECK(dist(L[0], to_d(P)) <= TOL * 64 * sc && dist(L[0], to_d(Q)) <= TOL * 64 * sc);
        CHECK(len(cross(unit(L[0].dir()), unit(to_d(Q - P)))) <= TOL * 64);
        CHECK(dot(L[0].dir(), to_d(cross(Point3D(A.n), Point3D(B.n)))) > 0);
        CHECK(std::fabs(dot(L[0].a, unit(L[0].dir()))) <= TOL * 64 * std::max<Real>(1, len(L[0].a)));  // a は原点に最も近い点
        CHECK(len(L[0].a) <= len(to_d(P)) * (1 + TOL * 64) + TOL);
        if (C <= 4) {  // 実数の平面でも同じ
            const std::vector<Line3D> Ld = cross_line(to_d(A), to_d(B));
            CHECK(Ld.size() == 1 && nearp3(Ld[0].a, L[0].a, TOL * 64) && len(cross(unit(Ld[0].dir()), unit(L[0].dir()))) <= TOL * 64);
        }
    }
    // 2 平面の角度・直線と平面の角度
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 200000000LL : 4;
        const Point3 p = rndp3(rng, C), q = rndp3(rng, C), r = rndp3(rng, C);
        Point3 u = rndp3(rng, C), v = rndp3(rng, C), w = rndp3(rng, C);
        if (it % 4 == 0) u = p + (q - p) * 2, v = q + (r - p) * 3 - (q - p), w = p + (r - p) - (q - p) * 2;  // 同じ平面
        if (it % 4 == 1) u = p + Point3(0, 0, 1), v = q + Point3(0, 0, 1), w = r + Point3(0, 0, 1);  // 平行移動
        if (collinear(p, q, r) || collinear(u, v, w)) continue;
        const Plane3 A = plane_from_points(p, q, r), B = plane_from_points(u, v, w);
        CHECK(cross_line(A, B).empty() == wzero(wcross(w3(A.n), w3(B.n))));
        const Real ab = angle_between(A, B);
        CHECK(0 <= ab && ab <= GEO_PI / 2 + TOL && std::fabs(ab - std::min(angle_between(A.n, B.n), GEO_PI - angle_between(A.n, B.n))) <= TOL);
        const Line3 l(p, u == p ? p + Point3(1, 0, 0) : u);
        const Real lp = angle_between(l, B);
        CHECK(std::fabs(lp - (GEO_PI / 2 - std::min(angle_between(l.dir(), B.n), GEO_PI - angle_between(l.dir(), B.n)))) <= TOL && angle_between(B, l) == lp);
        CHECK(std::fabs(angle_between(l, Line3(v, w)) - std::min(angle_between(l.dir(), w - v), GEO_PI - angle_between(l.dir(), w - v))) <= TOL);
    }
    CHECK(std::fabs(angle_between(Plane3(0, 0, 1, 0), Plane3(1, 0, 0, 5)) - GEO_PI / 2) <= TOL && angle_between(Line3(Point3(), Point3(0, 0, 1)), Plane3(0, 0, 1, 3)) == GEO_PI / 2);
}

// 192 ビットの a·b - c·e：打ち消し合っても正しく Real に直る
static void test_mul_sub(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 200000 * TEST_SCALE; ++it, ++cases) {
        const int sh = (int)rnd(rng, 0, 95);
        const Int128 a = (Int128)rnd(rng, -4000000000000000000LL, 4000000000000000000LL) * ((Int128)1 << (sh / 3)) + rnd(rng, -5, 5);
        const long long b = rnd(rng, -4000000000000000000LL, 4000000000000000000LL) >> rnd(rng, 0, 61);
        const Int128 d = rnd(rng, -1000000, 1000000);
        // a·b - a·(b - k) = a·k、(a + d)·b - a·b = d·b（積は 2^127 を超えうるが、差は小さい）
        const long long k = rnd(rng, -3, 3);
        if ((b > 0 && b - k < -4000000000000000000LL) || (b < 0 && b - k > 4000000000000000000LL)) continue;
        CHECK(near(geo_internal::mul_sub_real(a, b, a, b - k), (Real)(a * k), TOL));
        CHECK(near(geo_internal::mul_sub_real(a + d, b, a, b), (Real)(d * b), TOL));
        CHECK(geo_internal::mul_sub_real(a, b, a, b) == 0);
        if (std::llabs(b) < (1LL << 30) && a < ((Int128)1 << 95) && a > -((Int128)1 << 95)) CHECK(near(geo_internal::mul_sub_real(a, b, 0, 0), (Real)(a * b), TOL));
        CHECK(near(geo_internal::mul_sub_real(a, b, d, k), -geo_internal::mul_sub_real(d, k, a, b), 0));
    }
}

static void test_spheres(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const bool big = it % 3 == 0;
        const long long C = big ? 1000000000LL : 3;
        const Sphere a(rndp3(rng, C), rnd(rng, 0, big ? 1000000000LL : 4)), b(rndp3(rng, C), rnd(rng, 0, big ? 1000000000LL : 4));
        // 内外
        const Point3 x = rndp3(rng, C + 1);
        const Int128 dx = dist2_128(x, a.c), r2 = (Int128)a.r * a.r;
        CHECK(contains(a, x) == (dx < r2 ? INSIDE : (dx == r2 ? ON_BORDER : OUTSIDE)));
        if (!big) CHECK(contains(to_d(a), to_d(x)) == contains(a, x));
        // 位置関係（実数で求めた距離と比べる。小さい座標なら実数の版とも同じ）
        const int rel = sphere_relation(a, b);
        const Int128 d2 = dist2_128(a.c, b.c), s2 = ((Int128)a.r + b.r) * ((Int128)a.r + b.r), f2 = ((Int128)a.r - b.r) * ((Int128)a.r - b.r);
        const int e = (d2 == 0 && a.r == b.r) ? SPHERE_SAME
                                               : (d2 > s2 ? SPHERE_SEPARATE : (d2 == s2 ? SPHERE_EXTERNAL : (d2 > f2 ? SPHERE_INTERSECT : (d2 == f2 ? SPHERE_INTERNAL : SPHERE_CONTAIN))));
        CHECK(rel == e && sphere_relation(b, a) == rel);
        if (!big) CHECK(sphere_relation(to_d(a), to_d(b)) == rel);
        // 共通部分の体積：中心を結ぶ軸に沿って断面（同心の円 2 つ）の面積を数値積分
        const Real got = intersection_volume(a, b);
        CHECK(near(got, intersection_volume(b, a), 1e-12L));
        if (!big) {
            const Real d = std::sqrt((Real)d2), R1 = a.r, R2 = b.r;
            auto f = [&](Real t) {  // 中心 a から軸に沿って t の位置
                const Real q1 = R1 * R1 - t * t, q2 = R2 * R2 - (t - d) * (t - d);
                return q1 > 0 && q2 > 0 ? GEO_PI * std::min(q1, q2) : (Real)0;
            };
            const Real lo = std::max(-R1, d - R2), hi = std::min(R1, d + R2);
            Real ev = 0;
            if (hi > lo) {  // 断面が入れ替わる位置（2 つの球面が交わる平面）で区切る
                const Real m = d > 0 ? (d * d + R1 * R1 - R2 * R2) / (2 * d) : lo;
                if (lo < m && m < hi) ev = integrate(f, lo, m, 1e-13L) + integrate(f, m, hi, 1e-13L);
                else ev = integrate(f, lo, hi, 1e-13L);
            }
            CHECK(std::fabs(got - ev) <= 1e-9L * std::max<Real>(1, ev));
            CHECK(near(intersection_volume(to_d(a), to_d(b)), got, 1e-12L));
        } else {
            CHECK(got >= 0 && got <= 4 * GEO_PI / 3 * std::pow((Real)std::min(a.r, b.r), 3) * (1 + TOL));
        }
        // 直線との交点
        Point3 p = rndp3(rng, C + 2), q = rndp3(rng, C + 2);
        if (it % 4 == 0 && !big) p = a.c + Point3(a.r, rnd(rng, -2, 2), rnd(rng, -2, 2)), q = p + Point3(0, rnd(rng, -2, 2), rnd(rng, -2, 2));  // 接する直線など
        if (p == q) continue;
        const Line3 l(p, q);
        const std::vector<Point3D> cps = cross_points(a, l);
        const Point3 dd = q - p, ww = p - a.c;
        const Int128 A = norm2_128(dd), B = dot128(dd, ww), Cc = norm2_128(ww) - (Int128)a.r * a.r, D = B * B - A * Cc;
        CHECK((int)cps.size() == (D < 0 ? 0 : (D == 0 ? 1 : 2)));
        const Real sc = std::max<Real>(1, (Real)std::max<long long>(C, a.r));
        for (const Point3D& X : cps) CHECK(std::fabs(dist(to_d(a.c), X) - (Real)a.r) <= TOL * 1000 * sc && dist(to_d(l), X) <= TOL * 1000 * sc);
        if (cps.size() == 2) CHECK(dot(cps[1] - cps[0], to_d(dd)) > 0);
        if (!big) CHECK(cross_points(to_d(a), to_d(l)).size() == cps.size());
    }
    // 接する直前（細い共通部分）：正で、小さい
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const Coord R = rnd(rng, 1, 1000000000LL), r = rnd(rng, 1, 1000000000LL);
        const Sphere a(Point3(0, 0, 0), R), b(Point3(R + r - 1, 0, 0), r);
        const Real v = intersection_volume(a, b);
        // 高さ h ≈ 1 程度の球冠 2 つ：π h² r の程度
        CHECK(v > 0 && v < 2 * GEO_PI * (Real)std::min(R, r) * 1.01L);
        if (r < R) {
            const Sphere c(Point3(R - r + 1, 0, 0), r);  // 内接より 1 だけ外
            const Real full = 4 * GEO_PI / 3 * (Real)r * r * r, vi = intersection_volume(a, c);
            CHECK(vi <= full * (1 + TOL) && vi > full * 0.5L);
        }
    }
    CHECK(near(sphere_cap_volume(2, 4), 4 * GEO_PI / 3 * 8, TOL) && near(sphere_cap_volume(2, 2), 2 * GEO_PI / 3 * 8, TOL));
}

static void test_input() {
    {
        CinFromString in("1 2 3\n4 5 6 -7 8 9\n1.5 2 -3\n0.5 0 0 1 1 1\n");
        const Point3 a = read_point3();
        const std::vector<Point3> v = read_points3(2);
        const Point3D d = read_point3_d();
        const std::vector<Point3D> w = read_points3_d(2);
        CHECK(a == Point3(1, 2, 3) && v[0] == Point3(4, 5, 6) && v[1] == Point3(-7, 8, 9) && d == Point3D(1.5L, 2, -3) && w[1] == Point3D(1, 1, 1));
    }
    std::ostringstream os;
    os << Point3(1, -2, 3) << " " << Plane3(1, 2, 3, 4) << " " << Sphere(0, 0, 0, 5);
    CHECK(os.str() == "(1, -2, 3) Plane3[1 x + 2 y + 3 z = 4] Sphere[c=(0, 0, 0), r=5]");
}

// ---- 位置関係（直線どうし・交点・キー・平面どうし・直線と平面・球と平面・球体） ----
static Plane3 rnd_plane(std::mt19937_64& rng, long long C, Point3 pts[3]) {
    while (true) {
        for (int k = 0; k < 3; ++k) pts[k] = rndp3(rng, C);
        if (!wzero(wcross(wsub(w3(pts[1]), w3(pts[0])), wsub(w3(pts[2]), w3(pts[0]))))) return plane_from_points(pts[0], pts[1], pts[2]);
    }
}
static void test_relations3(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 2 : 5;
        Point3 a = rndp3(rng, C), b = rndp3(rng, C), c = rndp3(rng, C), d = rndp3(rng, C);
        if (it % 4 == 0) c = a + (b - a) * rnd(rng, -2, 2), d = c + (b - a) * rnd(rng, -2, 2) + (it % 8 == 0 ? rndp3(rng, 1) : Point3());  // 平行・同じ直線
        if (a == b || c == d) continue;
        const Line3 l(a, b), m(c, d);
        // 直線どうしの位置関係（愚直：連立方程式を有理数で）
        Int128 sn, tn, den;
        const int r = b_solve(a, b - a, c, d - c, sn, tn, den);
        const bool same = r == 2 && wzero(wcross(w3(b - a), wsub(w3(c), w3(a))));
        const int e = r == 1 ? LINE3_CROSS : (r == 2 ? (same ? LINE3_SAME : LINE3_PARALLEL) : LINE3_SKEW);
        CHECK(line_relation(l, m) == e && line_relation(m, l) == e);
        const std::vector<Point3D> cp = cross_point(l, m);
        const std::vector<FracPoint3> fp = cross_point_frac(l, m);
        if (e == LINE3_CROSS) {
            const Point3D x = to_d(a) + to_d(b - a) * ((Real)sn / (Real)den);
            CHECK(cp.size() == 1 && nearp3(cp[0], x, TOL));
            CHECK(fp.size() == 1 && fp[0].d > 0 && nearp3(fp[0].point_d(), x, TOL));
            // 既約で、両方の直線の上に誤差なくある
            const FracPoint3 f = fp[0];
            CHECK(geo_internal::gcd128_3(geo_internal::gcd128_3(geo_internal::gcd128_3(f.x, f.y), f.z), f.d) == 1);
            const W3 P{f.x, f.y, f.z};
            for (const Line3& g : {l, m}) {
                const W3 A{(Int128)g.a.x * f.d, (Int128)g.a.y * f.d, (Int128)g.a.z * f.d};
                CHECK(wzero(wcross(w3(g.dir()), wsub(P, A))));
            }
            // 別の表し方の直線でも同じ交点
            CHECK(cross_point_frac(Line3(b, a), Line3(d + (d - c) * 2, c))[0] == f && cross_point_frac(m, l)[0] == f);
            CHECK(FracPoint3Hash()(f) == FracPoint3Hash()(cross_point_frac(m, l)[0]));
        } else {
            CHECK(cp.empty() && fp.empty());
        }
        // 直線のキー
        const LineKey3 k1 = line_key(a, b), k2 = line_key(c, d);
        CHECK((k1 == k2) == same && (k1 < k2 || k2 < k1) == !same && k1 == line_key(b, a) && k1 == line_key(Line3(a + (b - a) * 3, a - (b - a) * 2)));
        CHECK(LineKey3Hash()(k1) == LineKey3Hash()(line_key(b, a)));
        CHECK(direction_key(b - a) == direction_key((b - a) * 3) && line_direction_key(b - a) == line_direction_key(a - b));
        CHECK((direction_key(b - a) == direction_key(d - c)) == (r == 2 && wdot(w3(b - a), w3(d - c)) > 0));
        // 線分どうしの交点
        const Segment3 s1(a, b), s2(c, d);
        const std::vector<Point3D> sp = cross_points(s1, s2);
        CHECK(!sp.empty() == b_meet(a, b, c, d, 0));
        for (const Point3D& q : sp) CHECK(b_seg_point(to_d(a), to_d(b), q) <= TOL * 16 && b_seg_point(to_d(c), to_d(d), q) <= TOL * 16);
        if (sp.size() == 2) CHECK(e == LINE3_SAME && dot(to_d(b - a), sp[1] - sp[0]) > 0);
        if (e == LINE3_SAME && b_meet(a, b, c, d, 0)) {  // 重なりの両端：もう一方の上にある端点のうち、s の向きで最小・最大
            Int128 lo = 0, hi = 0;
            bool first = true;
            for (const Point3& q : {a, b, c, d})
                if (b_on_seg(a, b, q) && b_on_seg(c, d, q)) {
                    const Int128 t = wdot(w3(b - a), wsub(w3(q), w3(a)));
                    if (first) lo = hi = t, first = false;
                    lo = std::min(lo, t), hi = std::max(hi, t);
                }
            CHECK(sp.size() == (lo == hi ? 1u : 2u));
        }
        const std::vector<Point3D> lp = cross_point(l, s2);
        CHECK(!lp.empty() == b_meet(a, b, c, d, 2) && cross_point(s2, l).size() == lp.size());
        // 平面
        Point3 P1[3], P2[3];
        const Plane3 p = rnd_plane(rng, C, P1);
        Plane3 q = rnd_plane(rng, C, P2);
        if (it % 5 == 0) q = Plane3(p.n * rnd(rng, -3, 3), p.d * rnd(rng, -3, 3)), P2[0] = Point3(1 << 20, 0, 0);  // 平行・同じ平面
        if (q.n == Point3()) continue;
        const bool par = wzero(wcross(w3(p.n), w3(q.n)));
        // 同じ平面：p の上の 3 点が q の上にある
        bool sm = par;
        for (int k = 0; k < 3; ++k) sm = sm && wdot(w3(q.n), w3(P1[k])) == q.d;
        CHECK(plane_relation(p, q) == (!par ? PLANE_CROSS : (sm ? PLANE_SAME : PLANE_PARALLEL)) && plane_relation(q, p) == plane_relation(p, q));
        CHECK((plane_key(p) == plane_key(q)) == sm && plane_key(Plane3(p.n * -2, p.d * -2)) == plane_key(p));
        if (par && !sm) CHECK(near(dist(p, q), dist(q, P1[0]), TOL * 16) && near(dist(q, p), dist(p, q), TOL * 16));
        else CHECK(dist(p, q) == 0);
        // 直線と平面
        const Int128 dn = wdot(w3(p.n), w3(b - a));
        const int lpr = dn != 0 ? LINE_PLANE_CROSS : (wdot(w3(p.n), w3(a)) == p.d ? LINE_PLANE_ON : LINE_PLANE_PARALLEL);
        CHECK(line_plane_relation(l, p) == lpr && intersects(l, p) == (lpr != LINE_PLANE_PARALLEL) && cross_point(l, p).size() == (lpr == LINE_PLANE_CROSS ? 1u : 0u));
        CHECK(near(dist(l, p), lpr == LINE_PLANE_PARALLEL ? dist(p, a) : 0, TOL));
        const Int128 va = wdot(w3(p.n), w3(a)) - p.d, vb = wdot(w3(p.n), w3(b)) - p.d;
        const bool si = (va <= 0 && vb >= 0) || (va >= 0 && vb <= 0);
        CHECK(intersects(s1, p) == si && near(dist(s1, p), si ? 0 : std::min(dist(p, a), dist(p, b)), TOL));
        // 球と平面・球体
        const Sphere S(rndp3(rng, C), rnd(rng, 0, C)), S2(rndp3(rng, C), rnd(rng, 0, C));
        const Int128 v = wdot(w3(p.n), w3(S.c)) - p.d, lhs = v * v, rhs = (Int128)S.r * S.r * norm2_128(p.n);
        const int spr = lhs > rhs ? 0 : (lhs == rhs ? 1 : 2);
        CHECK(sphere_plane_relation(S, p) == spr && ball_intersects(S, p) == (spr > 0));
        const std::vector<Circle3D> cc = cross_circle(S, p);
        CHECK(cc.size() == (spr > 0 ? 1u : 0u));
        if (!cc.empty()) {
            CHECK(on_plane(to_d(p), cc[0].c) && near(len(cc[0].n), 1, TOL) && near(std::fabs(dot(cc[0].n, unit(p.n))), 1, TOL));
            if (spr == 1) CHECK(cc[0].r == 0);
            // 円の上の点は球面と平面の上
            const Point3D u = unit(cross(cc[0].n, std::fabs(cc[0].n.x) < 0.5L ? Point3D(1, 0, 0) : Point3D(0, 1, 0)));
            const Point3D x = cc[0].c + u * cc[0].r;
            CHECK(near(dist(x, to_d(S.c)), (Real)S.r, TOL * 64) && dist(to_d(p), x) <= TOL * 64 * (1 + len(x)));
        }
        const std::vector<Circle3D> c2 = cross_circle(S, S2);
        const int rel = sphere_relation(S, S2);
        CHECK(c2.size() == ((rel == SPHERE_INTERSECT || rel == SPHERE_EXTERNAL || rel == SPHERE_INTERNAL) ? 1u : 0u));
        if (!c2.empty()) {
            const Point3D u = unit(cross(c2[0].n, std::fabs(c2[0].n.x) < 0.5L ? Point3D(1, 0, 0) : Point3D(0, 1, 0)));
            const Point3D x = c2[0].c + u * c2[0].r;
            CHECK(near(dist(x, to_d(S.c)), (Real)S.r, TOL * 64) && near(dist(x, to_d(S2.c)), (Real)S2.r, TOL * 64));
            if (rel != SPHERE_INTERSECT) CHECK(c2[0].r == 0);
        }
        // 球体と線分・直線（中心との距離の 2 乗を有理数で）
        auto seg_d2 = [&](const Point3& A, const Point3& B, bool line, Int128& num, Int128& dd) {
            const W3 dv = wsub(w3(B), w3(A)), w = wsub(w3(S.c), w3(A));
            const Int128 t = wdot(dv, w), L = wdot(dv, dv);
            if (!line && (L == 0 || t <= 0)) { num = wdot(w, w), dd = 1; return; }
            if (!line && t >= L) { const W3 w2 = wsub(w3(S.c), w3(B)); num = wdot(w2, w2), dd = 1; return; }
            const W3 cr = wcross(dv, w);
            num = wdot(cr, cr), dd = L;
        };
        Int128 nm, dd;
        const Int128 r2 = (Int128)S.r * S.r;
        seg_d2(a, b, false, nm, dd);
        CHECK(ball_intersects(S, s1) == (nm <= r2 * dd) && near(ball_dist(S, s1), std::max<Real>(0, std::sqrt((Real)nm / (Real)dd) - (Real)S.r), TOL));
        seg_d2(a, b, true, nm, dd);
        CHECK(ball_intersects(S, l) == (nm <= r2 * dd) && near(ball_dist(S, l), std::max<Real>(0, std::sqrt((Real)nm / (Real)dd) - (Real)S.r), TOL));
        CHECK(ball_contains(S, s1) == (dist2_128(a, S.c) <= r2 && dist2_128(b, S.c) <= r2));
        const Int128 cd2 = dist2_128(S.c, S2.c);
        CHECK(ball_intersects(S, S2) == (cd2 <= (Int128)(S.r + S2.r) * (S.r + S2.r)));
        CHECK(ball_contains(S, S2) == (S.r >= S2.r && cd2 <= (Int128)(S.r - S2.r) * (S.r - S2.r)));
        CHECK(near(ball_dist(S, a), std::max<Real>(0, dist(to_d(S.c), to_d(a)) - (Real)S.r), TOL));
        CHECK(near(ball_dist(S, S2), std::max<Real>(0, dist(to_d(S.c), to_d(S2.c)) - (Real)(S.r + S2.r)), TOL));
        CHECK(near(ball_dist(S, p), std::max<Real>(0, dist(p, S.c) - (Real)S.r), TOL));
        // 実数の型でも同じ判定
        CHECK(line_relation(to_d(l), to_d(m)) == e && line_plane_relation(to_d(l), to_d(p)) == lpr && plane_relation(to_d(p), to_d(q)) == plane_relation(p, q));
        CHECK(sphere_plane_relation(to_d(S), to_d(p)) == spr && ball_intersects(to_d(S), to_d(s1)) == ball_intersects(S, s1));
    }
    // 大きい座標：接する平面（法線が大きい）・交わる直線
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const Point3 c = rndp3(rng, 400000000LL);
        const long long k = rnd(rng, 1, 100000000), M = rnd(rng, 1, 1800000000000000000LL);
        const Sphere S(c, 5 * k);
        const Point3 n(3 * M, 4 * M, 0);
        const Int128 d0 = wdot(w3(n), w3(c)) + (Int128)25 * k * M;  // 点 c + (3k, 4k, 0) を通る接平面
        CHECK(sphere_plane_relation(S, Plane3(n, d0)) == 1 && sphere_plane_relation(S, Plane3(n, d0 - 1)) == 2 && sphere_plane_relation(S, Plane3(n, d0 + 1)) == 0);
        CHECK(sphere_plane_relation(S, Plane3(-n, -d0)) == 1 && ball_intersects(S, Plane3(n, d0)) && !ball_intersects(S, Plane3(n, d0 + 1)));
        // 交わる 2 直線：交点 X を通る
        const Point3 X = rndp3(rng, 900000000LL), u = rndp3(rng, 1000), v = rndp3(rng, 1000);
        if (u == Point3() || v == Point3() || cross(u, v) == Point3()) continue;
        const Line3 l(X - u * rnd(rng, 1, 100000), X + u * rnd(rng, 1, 100000)), m(X + v * rnd(rng, 1, 100000), X - v * rnd(rng, 1, 100000));
        CHECK(line_relation(l, m) == LINE3_CROSS && line_relation(l, Line3(m.a + Point3(0, 0, 1), m.b + Point3(0, 0, 1))) != LINE3_SAME);
        const std::vector<FracPoint3> f = cross_point_frac(l, m);
        CHECK(f.size() == 1 && f[0].is_integer() && f[0].x == X.x && f[0].y == X.y && f[0].z == X.z);
        CHECK(line_key(l.a, l.b) == line_key(X, X + u) && plane_relation(plane_from_points(l.a, l.b, m.a), plane_from_points(m.b, l.a, X)) == PLANE_SAME);
    }
    // 回転：型の違う中心・軸、-0 を出さない
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const Point3 p = rndp3(rng, 100), ax = rndp3(rng, 5);
        if (ax == Point3()) continue;
        const Point3D cd(rndr(rng, -10, 10), rndr(rng, -10, 10), rndr(rng, -10, 10));
        const Real th = rndr(rng, -7, 7);
        CHECK(nearp3(rotated(p, th, ax, cd), rotated(to_d(p), th, to_d(ax), cd), TOL * 16));
        CHECK(nearp3(rotated_deg(p, 90, ax, cd), rotated(to_d(p), GEO_PI / 2, ax, cd), TOL * 64));
    }
    for (const Point3D& q : {rotated_deg(Point3(1, 0, 0), 90, Point3(0, 0, 1)), rotated(Point3(0, 0, 0), 1, Point3(1, 1, 1)), rotated_deg(Point3D(0, -0.0L, 0), 180, Point3(1, 0, 0))})
        CHECK(!(std::signbit(q.x) && q.x == 0) && !(std::signbit(q.y) && q.y == 0) && !(std::signbit(q.z) && q.z == 0));
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_basic(rng, cases);
    test_relations3(rng, cases);
    test_volume(rng, cases);
    test_angle_rotation(rng, cases);
    test_lines(rng, cases);
    test_planes(rng, cases);
    test_mul_sub(rng, cases);
    test_spheres(rng, cases);
    test_input();
    std::cout << "[geo3d] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[geo3d] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 2000000;
    std::vector<Point3> ps(n + 3);
    for (auto& p : ps) p = rndp3(rng, 1000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    long long s = 0;
    for (int i = 0; i < n; ++i) s += orient3d(ps[i], ps[i + 1], ps[i + 2], ps[i + 3]);
    const long long t1 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    Real acc = 0;
    for (int i = 0; i + 3 < n; i += 2) acc += dist(Segment3(ps[i], ps[i + 1]), Segment3(ps[i + 2], ps[i + 3]));
    const long long t2 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i + 3 < n; i += 2) s += intersects(Segment3(ps[i], ps[i + 1]), Segment3(ps[i + 2], ps[i + 3]));
    const long long t3 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i + 3 < n; i += 2) {
        acc += intersection_volume(Sphere(ps[i], std::llabs(ps[i + 1].x)), Sphere(ps[i + 2], std::llabs(ps[i + 3].x)));
        s += (long long)cross_points(Sphere(ps[i], std::llabs(ps[i + 1].x)), Line3(ps[i + 2], ps[i + 3])).size();
    }
    const long long t4 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i + 1 < n; i += 2) acc += rotated(ps[i], (Real)i, ps[i + 1]).x;
    const long long t5 = ms_since(t0);
    std::cout << "[geo3d] benchmark: orient3d × " << n << " " << t1 << " ms, dist(線分, 線分) × " << n / 2 << " " << t2 << " ms, intersects(線分, 線分) × " << n / 2
              << " " << t3 << " ms, 球の共通部分 + 球と直線 × " << n / 2 << " " << t4 << " ms, rotated × " << n / 2 << " " << t5 << " ms  (check "
              << (s + (long long)std::fmod(acc, (Real)7)) % 7 << ")\n";
}

}  // namespace g3_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace g3_demo {

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
#define SHOW(expr, comment) g3_demo::show(#expr, g3_demo::to_s(expr), comment)
#define RUN(...)                   \
    g3_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static std::string to_s(const Point3& p) { return "(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ", " + std::to_string(p.z) + ")"; }
static std::string to_s(const Point3D& p) { return "(" + to_s(p.x) + ", " + to_s(p.y) + ", " + to_s(p.z) + ")"; }
static std::string to_s(const std::vector<Point3D>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
static std::string to_s(const Line3D& l) { return "Line3D[" + to_s(l.a) + " -> " + to_s(l.b) + "]"; }
static std::string where_s(Where w) { return w == INSIDE ? "INSIDE" : (w == ON_BORDER ? "ON_BORDER" : "OUTSIDE"); }
static void demo() {
    banner("3 次元 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧（Point3 は整数、Point3D は実数。判定は整数なら誤差なし）\n"
              << "  Point3 p(x, y, z);  dot / cross / norm2 / dist2 / len / dist / dot128 / norm2_128 / dist2_128 / to_d / midpoint / unit\n"
              << "  volume6 / tetra_volume / orient3d / coplanar / collinear / triangle_area / in_tetrahedron\n"
              << "  angle_between(_deg) / elevation_angle(_deg) / rot90_x / rot90_y / rot90_z / rotation24 / rotated / rotated_deg\n"
              << "  Line3 / Segment3：on_line / on_segment / is_parallel / is_orthogonal / intersects / projection / reflection / dist\n"
              << "  Plane3(a, b, c, d) / plane_from_points：side / on_plane / dist / projection / reflection / cross_point / cross_line / angle_between\n"
              << "  Sphere(c, r)：contains / sphere_relation / cross_points / intersection_volume / sphere_cap_volume\n"
              << "  read_point3 / read_points3 / read_point3_d / read_points3_d / Point3Hash\n"
              << "  位置関係: line_relation（LINE3_*）/ cross_point(l, m) / cross_points(s, t) / cross_point_frac / line_key（LineKey3）/ direction_key\n"
              << "            line_plane_relation（LINE_PLANE_*）/ plane_relation（PLANE_*）/ plane_key / dist(平面, 平面・直線・線分)\n"
              << "            sphere_plane_relation / cross_circle(球, 平面 or 球)（Circle3D）/ ball_intersects / ball_contains / ball_dist\n";

    api("【1】 Point3 / dot / cross / norm2 / dist2 / dist", "Coord / Point3 / Real", "2 次元の Point と同じ使い方。cross はベクトル（両方に垂直）", "O(1)");
    RUN(Point3 a(1, 0, 0), b(0, 2, 0), c(1, 2, 3));
    SHOW(dot(a, c), "1·1 + 0·2 + 0·3");
    SHOW(cross(a, b), "x 軸 × y 軸の向き = z 軸の向き");
    SHOW(dist2(a, c), "距離の 2 乗（整数のまま）");
    SHOW(dist(a, c), "距離");
    SHOW(c * 2 - a, "演算子");
    note("座標が 10^9 で差が大きいとき dot / norm2 / dist2 は long long をはみ出す → dot128 / norm2_128 / dist2_128（Int128）");

    api("【2】 volume6 / tetra_volume / orient3d / coplanar / collinear / triangle_area / in_tetrahedron", "Int128 / Real / int / bool / Where",
        "体積の 6 倍（符号付き・厳密）/ 体積 / d が平面 abc の表 +1・裏 -1・同じ平面 0 / 4 点が同じ平面 / 3 点が一直線 / 面積 / 四面体の内外", "O(1)");
    RUN(Point3 o(0, 0, 0), ex(1, 0, 0), ey(0, 1, 0), ez(0, 0, 1));
    SHOW(volume6(o, ex, ey, ez), "6 × 1/6");
    SHOW(tetra_volume(o, ex, ey, ez), "1/6");
    SHOW(orient3d(o, ex, ey, ez), "o → ex → ey が反時計回りに見える側（z > 0）に ez");
    SHOW(orient3d(o, ex, ey, -ez), "裏側");
    SHOW(coplanar(o, ex, ey, ex + ey * 5), "z = 0 の平面の上");
    SHOW(collinear(o, Point3(1, 2, 3), Point3(-2, -4, -6)), "");
    SHOW(triangle_area(ex, ey, ez), "正三角形（1 辺 √2）= √3/2");
    SHOW(where_s(in_tetrahedron(o, ex * 4, ey * 4, ez * 4, Point3(1, 1, 1))), "x + y + z < 4");
    SHOW(where_s(in_tetrahedron(o, ex * 4, ey * 4, ez * 4, Point3(2, 2, 0))), "面 z = 0 の上");

    api("【3】 angle_between(u, v) / elevation_angle(v)（_deg で度）", "Real", "ベクトルのなす角 [0, π] / 仰角（xy 平面から上へ測った角 [-π/2, π/2]）", "O(1)");
    SHOW(angle_between_deg(Point3(1, 0, 0), Point3(1, 1, 0)), "45 度");
    SHOW(angle_between_deg(Point3(1, 1, 1), Point3(-1, -1, -1)), "反対向き");
    SHOW(elevation_angle_deg(Point3(1, 0, 1)), "45 度上向き");
    SHOW(elevation_angle_deg(Point3(0, 3, -4)), "下向きは負（俯角は符号を取る）");
    note("典型90 018（Statue of Chokudai）: 観覧車の位置 P（Point3D）と像の位置 S から、見下ろす角 = elevation_angle_deg(P - S)");
    {
        // 例：直径 L、周期 T の観覧車（中心 (0, 0, L/2)、y-z 平面を回る）。時刻 e の位置
        const Real T = 4, L = 2, X = 1, Y = 1, e = 1;
        const Point3D P(0, -L / 2 * std::sin(2 * GEO_PI * e / T), L / 2 - L / 2 * std::cos(2 * GEO_PI * e / T)), S(X, Y, 0);
        code("Point3D P(0, -L/2·sin(2πe/T), L/2 - L/2·cos(2πe/T)), S(X, Y, 0);   // T = 4, L = 2, X = Y = 1, e = 1");
        SHOW(P, "観覧車の位置");
        SHOW(elevation_angle_deg(P - S), "像から見上げる角 = 観覧車から見下ろす角（度）");
    }

    api("【4】 rot90_x / rot90_y / rot90_z / rotation24 / rotated(p, θ, 軸 [, 中心]) / rotated_deg", "Point3 / Point3D",
        "右ねじの向き（軸の先から見て反時計回り）に回す。rotation24(p, k) は立方体の 24 通りの向き（k = 0〜23）", "O(1)");
    SHOW(rot90_z(Point3(1, 2, 3)), "z 軸のまわりに 90 度（x → y）");
    SHOW(rot90_x(Point3(1, 2, 3)), "x 軸のまわり（y → z）");
    SHOW(rotation24(Point3(1, 2, 3), 7), "24 通りのうち 7 番目");
    SHOW(rotated_deg(Point3(1, 0, 0), 90, Point3(0, 0, 1)), "座標軸の 90 度の倍数は誤差なし");
    SHOW(rotated_deg(Point3(1, 0, 0), 120, Point3(1, 1, 1)), "(1,1,1) のまわりに 120 度 → x → y");
    SHOW(rotated(Point3(2, 0, 0), GEO_PI, Point3(0, 1, 0), Point3(1, 0, 0)), "中心 (1, 0, 0)・y 軸の向きのまわりに 180 度");
    note("立体のブロックを回して重ねる問題：全部の点に rotation24(p, k) を k = 0〜23 で（裏返しも許すなら x を -x にしたものも）");

    api("【5】 Line3 / Segment3：intersects / on_segment / dist / projection", "bool / Real / Point3D",
        "交わるか（整数なら厳密。ねじれの位置も正しく）/ 線分の上か / 距離（点・直線・線分のどの組でも）/ 垂線の足", "O(1)");
    RUN(Segment3 s1(Point3(0, 0, 0), Point3(2, 2, 0)), s2(Point3(0, 2, 0), Point3(2, 0, 0)), s3(Point3(0, 2, 1), Point3(2, 0, 1)));
    SHOW(intersects(s1, s2), "(1, 1, 0) で交わる");
    SHOW(intersects(s1, s3), "高さが違う（ねじれの位置）");
    SHOW(dist(s1, s3), "ねじれの位置の最短距離");
    SHOW(dist(Line3(Point3(0, 0, 0), Point3(1, 0, 0)), Point3(3, 4, 0)), "x 軸から (3, 4, 0) まで");
    SHOW(dist(s1, Point3(5, 5, 0)), "線分の端 (2, 2, 0) が最も近い");
    SHOW(projection(Line3(Point3(0, 0, 0), Point3(1, 1, 1)), Point3(3, 0, 0)), "垂線の足");
    SHOW(on_segment(s1, Point3(1, 1, 0)), "");

    api("【6】 Plane3(a, b, c, d) / plane_from_points(p, q, r)：side / dist / projection / cross_point / cross_line / angle_between", "int / Real / Point3D / vector",
        "a x + b y + c z = d の平面。表 +1・裏 -1（法線の向きが表）/ 距離 / 垂線の足 / 直線・線分との交点 / 2 平面の交線 / 角度", "O(1)");
    RUN(Plane3 pl = plane_from_points(Point3(1, 0, 0), Point3(0, 1, 0), Point3(0, 0, 1)));
    std::cout << "  pl                                 -> " << pl << "   // x + y + z = 1\n";
    SHOW(side(pl, Point3(1, 1, 1)), "表（原点と反対側）");
    SHOW(dist(pl, Point3(0, 0, 0)), "1/√3");
    SHOW(projection(pl, Point3(0, 0, 0)), "(1/3, 1/3, 1/3)");
    SHOW(cross_point(Line3(Point3(0, 0, 0), Point3(1, 0, 0)), pl), "x 軸との交点");
    SHOW(cross_point(Segment3(Point3(0, 0, 0), Point3(0, 0, -5)), pl), "交わらない → 空");
    SHOW(to_s(cross_line(Plane3(0, 0, 1, 0), Plane3(1, 0, 0, 2))[0]), "z = 0 と x = 2 の交線（a は原点に最も近い点、向きは n1 × n2）");
    SHOW(rad_to_deg(angle_between(Plane3(0, 0, 1, 0), pl)), "2 面角（度）");
    note("平面の式が a x + b y + c z + d = 0 の形なら Plane3(a, b, c, -d)");

    api("【7】 Sphere(c, r)：contains / sphere_relation / cross_points / intersection_volume", "Where / int / vector<Point3D> / Real",
        "球の内外 / 2 球の位置関係（SPHERE_SEPARATE(4) / EXTERNAL(3) / INTERSECT(2) / INTERNAL(1) / CONTAIN(0) / SAME(5)）/ 直線との交点 / 共通部分の体積", "O(1)");
    RUN(Sphere A(Point3(0, 0, 0), 2), B(Point3(3, 0, 0), 2));
    SHOW(where_s(contains(A, Point3(0, 2, 0))), "球面の上");
    SHOW(sphere_relation(A, B), "SPHERE_INTERSECT = 2");
    SHOW(sphere_relation(A, Sphere(Point3(4, 0, 0), 2)), "外接 SPHERE_EXTERNAL = 3");
    SHOW(cross_points(A, Line3(Point3(-5, 0, 0), Point3(5, 0, 0))), "x 軸と 2 点で（直線の向きの順）");
    SHOW(intersection_volume(A, B), "共通部分（高さ 0.5 の球冠 2 つ）");
    SHOW(sphere_cap_volume(2, 0.5), "その片方");

    api("【8】 read_point3() / read_points3(N) / read_point3_d() / Point3Hash", "Point3 / vector<Point3> / Point3D", "入力（x y z の順）。unordered_set<Point3, Point3Hash>", "O(N)");
    {
        CinFromString in("2\n1 2 3\n-4 5 6\n");
        int N;
        std::cin >> N;
        const std::vector<Point3> ps = read_points3(N);
        code("int N; cin >> N; vector<Point3> ps = read_points3(N);   // 入力 \"2 / 1 2 3 / -4 5 6\"");
        SHOW(ps[1], "");
        std::unordered_set<Point3, Point3Hash> st(ps.begin(), ps.end());
        code("unordered_set<Point3, Point3Hash> st(ps.begin(), ps.end());");
        SHOW(st.count(Point3(1, 2, 3)), "");
    }

    api("【9】 line_relation(l, m) / cross_point(l, m) / cross_points(s, t) / cross_point_frac(l, m) / line_key(p, q)",
        "int（LINE3_CROSS 0 / LINE3_PARALLEL 1 / LINE3_SAME 2 / LINE3_SKEW 3）/ vector<Point3D> / vector<FracPoint3> / LineKey3",
        "2 直線の位置関係（ねじれの位置も）/ 交点 / 交点を既約分数で / 同じ直線なら同じキー（整数なら厳密）", "O(1)");
    RUN(Line3 l(Point3(0, 0, 0), Point3(1, 1, 0)); Line3 m(Point3(2, 0, 0), Point3(2, 1, 0)));
    SHOW(line_relation(l, m) == LINE3_CROSS, "同じ平面 z = 0 で交わる");
    SHOW(to_s(cross_point(l, m)), "交点 (2, 2, 0)");
    SHOW(line_relation(l, Line3(Point3(0, 0, 1), Point3(0, 1, 1))) == LINE3_SKEW, "z = 1 の y 軸に平行な直線とはねじれ");
    SHOW(line_relation(l, Line3(Point3(3, 3, 0), Point3(5, 5, 0))) == LINE3_SAME, "同じ直線");
    {
        const FracPoint3 f = cross_point_frac(Line3(Point3(0, 0, 0), Point3(3, 1, 0)), Line3(Point3(0, 1, 0), Point3(1, 0, 0)))[0];
        std::cout << "  cross_point_frac(...)              -> (" << int128_to_string(f.x) << ", " << int128_to_string(f.y) << ", " << int128_to_string(f.z) << ") / "
                  << int128_to_string(f.d) << "   // (3/4, 1/4, 0)。set<FracPoint3> で異なる交点を数える\n";
    }
    SHOW(line_key(Point3(1, 1, 1), Point3(3, 5, 7)) == line_key(Point3(5, 9, 13), Point3(2, 3, 4)), "同じ直線の別の 2 点 → 同じキー");
    SHOW(to_s(cross_points(Segment3(Point3(0, 0, 0), Point3(4, 0, 0)), Segment3(Point3(2, 0, 0), Point3(9, 0, 0)))), "一直線で重なる部分の両端");
    note("ABC301 G 型（3 次元の直線の交点・同じ直線をまとめる）: set<LineKey3> と cross_point_frac");

    api("【10】 line_plane_relation(l, pl) / plane_relation(p, q) / plane_key(pl) / dist(平面, 平面) / intersects(線分, 平面)",
        "int（LINE_PLANE_CROSS / PARALLEL / ON、PLANE_CROSS / PARALLEL / SAME）/ Plane3 / Real / bool",
        "直線と平面・平面どうしの位置関係（整数なら厳密）/ 同じ平面なら同じキー（gcd で割って向きをそろえる）/ 平行な平面の間隔", "O(1)");
    RUN(Plane3 p(0, 0, 1, 2));
    note("p = 平面 z = 2");
    SHOW(line_plane_relation(Line3(Point3(0, 0, 0), Point3(1, 0, 0)), p) == LINE_PLANE_PARALLEL, "x 軸は平行（交わらない）");
    SHOW(line_plane_relation(Line3(Point3(0, 0, 2), Point3(1, 5, 2)), p) == LINE_PLANE_ON, "平面の上の直線");
    SHOW(plane_relation(p, Plane3(0, 0, -3, -6)) == PLANE_SAME, "-3z = -6 は同じ平面");
    SHOW(plane_relation(p, Plane3(0, 0, 2, 9)) == PLANE_PARALLEL, "2z = 9 は平行");
    SHOW(dist(p, Plane3(0, 0, 2, 9)), "z = 2 と z = 4.5 の間隔");
    {
        std::ostringstream os;
        os << plane_key(Plane3(0, 0, -3, -6));
        std::cout << "  plane_key(Plane3(0, 0, -3, -6))    -> " << os.str() << "   // z = 2 にそろう\n";
    }
    SHOW(intersects(Segment3(Point3(0, 0, 0), Point3(0, 0, 5)), p), "線分が平面を貫く");

    api("【11】 sphere_plane_relation(sp, pl) / cross_circle(sp, pl or sp2) / ball_intersects / ball_contains / ball_dist",
        "int（0 離れる / 1 接する / 2 交わる）/ vector<Circle3D{c, n, r}> / bool / Real",
        "球面と平面の位置関係（整数なら 256 ビットで厳密）/ 交わりの円 / 球体（中身の詰まった球）としての共有点・包含・距離", "O(1)");
    RUN(Sphere S(Point3(0, 0, 0), 5));
    SHOW(sphere_plane_relation(S, Plane3(0, 0, 1, 5)), "z = 5 で接する");
    SHOW(sphere_plane_relation(S, Plane3(0, 0, 1, 3)), "z = 3 で交わる");
    {
        const Circle3D cc = cross_circle(S, Plane3(0, 0, 1, 3))[0];
        std::cout << "  cross_circle(S, z = 3)             -> c = " << to_s(cc.c) << ", n = " << to_s(cc.n) << ", r = " << to_s(cc.r) << "   // 半径 4 の円\n";
        const Circle3D e = cross_circle(S, Sphere(Point3(6, 0, 0), 5))[0];
        std::cout << "  cross_circle(S, 球 (6,0,0) r=5)    -> c = " << to_s(e.c) << ", r = " << to_s(e.r) << "   // x = 3 の平面の半径 4 の円\n";
    }
    SHOW(ball_intersects(S, Segment3(Point3(1, 1, 1), Point3(2, 2, 2))), "球体の中の線分（球面とは交わらない）");
    SHOW(ball_contains(S, Segment3(Point3(3, 4, 0), Point3(0, 0, -5))), "両端が球面の上でも入る");
    SHOW(ball_dist(S, Line3(Point3(8, 0, -1), Point3(8, 0, 1))), "直線 x = 8 まで 3");
    SHOW(ball_contains(S, Sphere(Point3(1, 0, 0), 4)), "内接する球体は入る");
}

}  // namespace g3_demo


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
    g3_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    g3_test::benchmark();    // 速度計測
    g3_demo::demo();         // 使い方デモ
    return 0;
}

