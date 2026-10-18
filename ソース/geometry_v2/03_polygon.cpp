// =============================================================================
//  多角形・三角形（面積・内外判定・位置関係・凸か・単純か・重心・格子点（ピックの定理）・三角形の種類）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 03 多角形・三角形
//  依存       : 00_int_common, 00_geo_base, 02_line（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - 多角形は vector<Point>（頂点を順に。始点を末尾に繰り返さない）。向きは時計回りでも反時計回りでもよい
//    - 整数の点なら判定・面積の 2 倍は厳密（|x|, |y| ≤ 10^9 の単純多角形なら面積の 2 倍は long long に収まる）
//    - contains / is_convex などは凸でない多角形にも使える（辺が交わらない「単純多角形」が前提。is_simple で確かめられる）
//
//  機能一覧                                                                    計算量
//    area2(poly)                      符号付き面積の 2 倍（反時計回りなら正。整数のまま）       O(N)
//    area(poly)                       面積（Real。|area2| / 2）                                O(N)
//    is_ccw(poly) / make_ccw(poly)    反時計回りか / 反時計回りにした多角形                    O(N)
//    contains(poly, p)                OUTSIDE / ON_BORDER / INSIDE（凸でなくてよい）           O(N)
//    is_convex(poly, strict = false)  凸か（strict = true なら 180 度の頂点・同じ点も不可）    O(N)
//    is_simple(poly)                  辺どうしが（隣の辺との共有点以外で）交わらないか         O(N^2)
//    perimeter(poly) / centroid(poly) 周の長さ / 重心（面積の重心。面積 0 なら頂点の平均）     O(N)
//    boundary_lattice_points(poly)    周の上の格子点の数（Σ gcd）                              O(N log)
//    interior_lattice_points(poly)    内部の格子点の数（ピックの定理。単純多角形）             O(N log)
//    triangle_area2(a, b, c)          三角形の符号付き面積の 2 倍（= cross(a, b, c)）          O(1)
//    triangle_kind(a, b, c)           TRIANGLE_ACUTE / TRIANGLE_RIGHT / TRIANGLE_OBTUSE / TRIANGLE_DEGENERATE（厳密）
//    in_triangle(a, b, c, p)          三角形の内外（OUTSIDE / ON_BORDER / INSIDE）             O(1)
//   位置関係（多角形は周を含む閉じた領域として。凸でなくてよい。整数なら厳密）
//    contains(poly, s)                線分 s が多角形に入っているか（周の上を通ってもよい）     O(N log N)
//    contains(P, Q)                   多角形 Q が多角形 P に入っているか                       O(N M log N)
//    intersects(poly, s / l / r)      線分 / 直線 / 半直線と多角形が共有点を持つか（接するのも含む）  O(N)
//    intersects(P, Q)                 2 つの多角形が共有点を持つか（辺が交わる・接する・片方が中に入る）  O(N M)
//    dist(poly, p / s) / dist(P, Q)   距離（中にある・共有点があれば 0）/ boundary_dist(poly, p)：周までの距離  O(N) / O(N M)
//    PolygonPrefix pp(poly)           連続する頂点 v[i..j] で囲む部分の面積・重心を O(1) で（凸多角形を対角線で切る問題に）
//      pp.area2(i, j)（Int128）/ pp.area(i, j) / pp.centroid(i, j)   番号は n で割った余り。j = i - 1 なら全体    O(1)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    pg_test::random_test()  … 愚直な実装（角度の和・全探索・三角形分割）と突き合わせるランダムテスト
//    pg_test::benchmark()    … 速度計測
//    pg_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
namespace pg_test {

#define PG_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[polygon] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK PG_CHECK

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
static bool b_meet(const Point& a, const Point& b, const Point& c, const Point& d) {
    if (b_on_seg(a, b, c) || b_on_seg(a, b, d) || b_on_seg(c, d, a) || b_on_seg(c, d, b)) return true;
    return sgn128(cross128(b - a, c - a)) * sgn128(cross128(b - a, d - a)) < 0 && sgn128(cross128(d - c, a - c)) * sgn128(cross128(d - c, b - c)) < 0;
}
static bool b_simple(const std::vector<Point>& P) {
    const int n = (int)P.size();
    if (n < 3) return false;
    for (int i = 0; i < n; ++i) {
        const Point a = P[i], b = P[(i + 1) % n], c = P[(i + 2) % n];
        if (a == b) return false;
        if (n >= 3 && cross128(a - b, c - b) == 0 && dot128(a - b, c - b) > 0) return false;  // 隣の辺と重なる
    }
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            const bool adj = (j == i + 1) || (i == 0 && j == n - 1);
            if (adj) continue;
            if (b_meet(P[i], P[(i + 1) % n], P[j], P[(j + 1) % n])) return false;
        }
    return true;
}
static Int128 b_area2(const std::vector<Point>& P) {
    Int128 s = 0;
    const int n = (int)P.size();
    for (int i = 0; i < n; ++i) s += cross128(P[i], P[(i + 1) % n]);
    return s;
}
// 角度の和で巻き数を求める（周の上は先に調べる）
static Where b_contains(const std::vector<Point>& P, const Point& p) {
    const int n = (int)P.size();
    for (int i = 0; i < n; ++i)
        if (b_on_seg(P[i], P[(i + 1) % n], p)) return ON_BORDER;
    Real s = 0;
    for (int i = 0; i < n; ++i) {
        const Point u = P[i] - p, v = P[(i + 1) % n] - p;
        s += std::atan2((Real)cross128(u, v), (Real)dot128(u, v));
    }
    return std::fabs(s) > GEO_PI ? INSIDE : OUTSIDE;
}
static std::vector<Point> rnd_star(std::mt19937_64& rng, int n, long long C) {
    // 原点のまわりに偏角の順で並べた多角形（同じ向きの点は 1 つにする）。単純とは限らないので使う側で確かめる
    std::vector<Point> ps;
    for (int k = 0; k < n; ++k) {
        Point p = rndp(rng, C);
        if (p != Point()) ps.push_back(p);
    }
    std::sort(ps.begin(), ps.end(), [](const Point& a, const Point& b) {
        auto h = [](const Point& v) { return (v.y > 0 || (v.y == 0 && v.x > 0)) ? 0 : 1; };
        if (h(a) != h(b)) return h(a) < h(b);
        const Int128 c = cross128(a, b);
        if (c != 0) return c > 0;
        return norm2_128(a) < norm2_128(b);
    });
    std::vector<Point> r;
    for (const Point& p : ps)
        if (r.empty() || !(cross128(r.back(), p) == 0 && dot128(r.back(), p) > 0)) r.push_back(p);
        else if (rng() % 2) r.back() = p;
    if (rng() % 2) std::reverse(r.begin(), r.end());
    if (!r.empty()) std::rotate(r.begin(), r.begin() + (long)(rng() % r.size()), r.end());
    return r;
}

static void test_area_contains(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const bool big = it % 3 == 0;
        const long long C = big ? 1000000000LL : (it % 3 == 1 ? 3 : 8);
        std::vector<Point> P = (it % 4 == 0) ? rndps(rng, (int)rnd(rng, 0, 7), C) : rnd_star(rng, (int)rnd(rng, 3, 12), C);
        if (it % 5 == 0 && !P.empty()) P.insert(P.begin() + (long)(rng() % P.size()), P[rng() % P.size()]);  // 同じ点・とげ
        const bool simple = b_simple(P);
        if (big && !simple) continue;  // 大きい座標は単純多角形だけ（面積の 2 倍が long long に収まる）
        CHECK((Int128)area2(P) == ((int)P.size() < 3 ? 0 : b_area2(P)));
        CHECK(near(area(P), std::fabs((Real)b_area2(P)) / 2, TOL) && is_ccw(P) == (b_area2(P) > 0));
        const std::vector<Point> Q = make_ccw(P);
        CHECK(area2(Q) >= 0 && (Q == P || (area2(P) < 0 && std::equal(Q.begin(), Q.end(), P.rbegin()))));
        CHECK(is_simple(P) == simple);
        if (!simple) continue;
        for (int q = 0; q < 30; ++q) {
            Point p = rndp(rng, C + 1);
            if (q % 5 == 0) p = P[rng() % P.size()];
            if (q % 5 == 1 && !big) {  // 辺の上
                const size_t i = rng() % P.size();
                const Point d = P[(i + 1) % P.size()] - P[i];
                const Coord g = std::gcd(std::llabs(d.x), std::llabs(d.y));
                if (g > 0) p = P[i] + Point(d.x / g, d.y / g) * rnd(rng, 0, g);
            }
            CHECK(contains(P, p) == b_contains(P, p));
        }
        // 周の長さ・重心
        Real per = 0;
        for (size_t i = 0; i < P.size(); ++i) per += std::sqrt((Real)norm2_128(P[(i + 1) % P.size()] - P[i]));
        CHECK(near(perimeter(P), per, TOL * 4));
        Real cx = 0, cy = 0, aw = 0;
        for (size_t i = 1; i + 1 < P.size(); ++i) {
            const Real w = (Real)cross128(P[i] - P[0], P[i + 1] - P[0]);
            cx += w * ((Real)P[0].x + (Real)P[i].x + (Real)P[i + 1].x) / 3, cy += w * ((Real)P[0].y + (Real)P[i].y + (Real)P[i + 1].y) / 3, aw += w;
        }
        const PointD g = centroid(P);
        if (aw != 0) CHECK(nearp(g, PointD(cx / aw, cy / aw), TOL * 1e3L));
        // 内側の格子点（小さい座標で全部数える）
        if (!big) {
            Coord in = 0, on = 0;
            for (Coord x = -C; x <= C; ++x)
                for (Coord y = -C; y <= C; ++y) {
                    const Where w = b_contains(P, Point(x, y));
                    in += (w == INSIDE), on += (w == ON_BORDER);
                }
            CHECK(boundary_lattice_points(P) == on && interior_lattice_points(P) == in);
        }
    }
    // 特別な形
    const std::vector<Point> sq = {{0, 0}, {4, 0}, {4, 4}, {0, 4}};
    CHECK(area2(sq) == 32 && contains(sq, Point(2, 2)) == INSIDE && contains(sq, Point(4, 1)) == ON_BORDER && contains(sq, Point(5, 1)) == OUTSIDE);
    CHECK(interior_lattice_points(sq) == 9 && boundary_lattice_points(sq) == 16);
    CHECK(contains(std::vector<Point>{}, Point(0, 0)) == OUTSIDE && contains(std::vector<Point>{{1, 1}}, Point(1, 1)) == ON_BORDER);
}

static void test_convex(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 1 : (it % 3 == 1 ? 2 : 1000000000LL);
        std::vector<Point> P = (it % 2) ? rndps(rng, (int)rnd(rng, 0, 6), C) : rnd_star(rng, (int)rnd(rng, 3, 8), C);
        if (it % 7 == 0 && !P.empty()) {  // 続けて同じ点
            const size_t i = rng() % P.size();
            P.insert(P.begin() + (long)i, P[i]);
        }
        for (int strict = 0; strict < 2; ++strict) {
            // 愚直：続けて同じ点をまとめてから（strict なら同じ点があれば凸でない）、単純・面積 > 0・どの辺についても全頂点が同じ側
            std::vector<Point> Q;
            for (const Point& p : P)
                if (Q.empty() || Q.back() != p) Q.push_back(p);
            while (Q.size() > 1 && Q.back() == Q.front()) Q.pop_back();
            bool e = Q.size() >= 3 && !(strict && Q.size() != P.size()) && b_simple(Q) && b_area2(Q) != 0;
            if (e) {
                const int sg = b_area2(Q) > 0 ? 1 : -1;
                for (size_t i = 0; i < Q.size() && e; ++i)
                    for (size_t j = 0; j < Q.size(); ++j) {
                        const int s = sgn128(cross128(Q[(i + 1) % Q.size()] - Q[i], Q[j] - Q[i]));
                        if (s == -sg) e = false;
                    }
                if (strict)
                    for (size_t i = 0; i < Q.size(); ++i)
                        if (cross128(Q[(i + 1) % Q.size()] - Q[i], Q[(i + 2) % Q.size()] - Q[i]) == 0) e = false;
            }
            CHECK(is_convex(P, strict) == e);
        }
    }
    const std::vector<Point> star = {{0, 3}, {2, -3}, {-3, 1}, {3, 1}, {-2, -3}};  // 星形五角形（どの角も同じ向きに曲がる）
    CHECK(!is_convex(star) && !is_simple(star));
}

static void test_triangle(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 3 : 1000000000LL;
        const Point a = rndp(rng, C), b = rndp(rng, C);
        Point c = rndp(rng, C);
        if (it % 4 == 0) c = a + rot90(b - a) * rnd(rng, -1, 1);  // 直角
        // 愚直：一番長い辺の 2 乗と、ほかの 2 辺の 2 乗の和を比べる
        int e;
        if (cross128(b - a, c - a) == 0) e = TRIANGLE_DEGENERATE;
        else {
            Int128 s[3] = {norm2_128(b - c), norm2_128(c - a), norm2_128(a - b)};
            std::sort(s, s + 3);
            e = s[2] == s[0] + s[1] ? TRIANGLE_RIGHT : (s[2] > s[0] + s[1] ? TRIANGLE_OBTUSE : TRIANGLE_ACUTE);
        }
        CHECK(triangle_kind(a, b, c) == e && triangle_area2(a, b, c) == (Coord)cross128(b - a, c - a));
        const Point p = rndp(rng, C);
        CHECK(in_triangle(a, b, c, p) == b_contains({a, b, c}, p) || cross128(b - a, c - a) == 0);
        if (cross128(b - a, c - a) == 0) CHECK(in_triangle(a, b, c, p) == ((b_on_seg(a, b, p) || b_on_seg(b, c, p) || b_on_seg(c, a, p)) ? ON_BORDER : OUTSIDE));
        if (C < 10) CHECK(triangle_kind(to_d(a), to_d(b), to_d(c)) == e);
    }
}

static void test_real(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const std::vector<Point> P = rnd_star(rng, (int)rnd(rng, 3, 10), 6);
        if (!b_simple(P)) continue;
        const Real th = rndr(rng, 0, 2 * GEO_PI), sc = rndr(rng, 0.01L, 100);
        const PointD sh(rndr(rng, -1000, 1000), rndr(rng, -1000, 1000));
        auto tf = [&](const Point& p) { return rotated(p, th) * sc + sh; };
        std::vector<PointD> D;
        for (const Point& p : P) D.push_back(tf(p));
        CHECK(near(area(D), area(P) * sc * sc, TOL * 1e4L) && is_ccw(D) == is_ccw(P));
        CHECK(is_convex(D) == is_convex(P) && is_convex(D, true) == is_convex(P, true) && is_simple(D) == is_simple(P));
        CHECK(nearp(centroid(D), rotated(centroid(P), th) * sc + sh, TOL * 1e4L));
        for (int q = 0; q < 20; ++q) {
            const Point p = rndp(rng, 7);
            CHECK(contains(D, tf(p)) == contains(P, p));
        }
        // 位置関係（実数の型でも同じ判定）
        for (int q = 0; q < 6; ++q) {
            Segment sg(rndp(rng, 7), rndp(rng, 7));
            if (q < 2) sg = Segment(P[rng() % P.size()], P[rng() % P.size()]);
            const SegmentD sd(tf(sg.a), tf(sg.b));
            CHECK(contains(D, sd) == contains(P, sg) && intersects(D, sd) == intersects(P, sg));
            if (sg.a != sg.b) CHECK(intersects(D, LineD(sd.a, sd.b)) == intersects(P, Line(sg.a, sg.b)) && intersects(D, RayD(sd.a, sd.b)) == intersects(P, Ray(sg.a, sg.b)));
            CHECK(near(dist(D, sd), dist(P, sg) * sc, TOL * 1e4L));
        }
    }
}

// ---- 位置関係の愚直な判定：線分の上の「区切りの点」（両端・頂点・辺との交点）とその中点を、有理数のまま調べる ----
struct BFrac {
    Int128 n, d;  // n / d（d > 0）
};
static Where b_where_at(const std::vector<Point>& P, const Segment& s, const BFrac& t) {  // 点 s.a + (s.b - s.a)·t の内外
    std::vector<Point> S;
    for (const Point& p : P) S.push_back(Point((Coord)(p.x * t.d), (Coord)(p.y * t.d)));
    const Point d = s.b - s.a;
    return b_contains(S, Point((Coord)(s.a.x * t.d + d.x * t.n), (Coord)(s.a.y * t.d + d.y * t.n)));
}
// {線分の全部が多角形に入る, 線分のどこかが多角形に入る}
static std::pair<bool, bool> b_seg_in(const std::vector<Point>& P, const Segment& s) {
    if (s.a == s.b) {
        const bool w = b_contains(P, s.a) != OUTSIDE;
        return {w, w};
    }
    const Point d = s.b - s.a;
    const int n = (int)P.size();
    std::vector<BFrac> ts = {{0, 1}, {1, 1}};
    for (int i = 0; i < n; ++i) {
        const Point u = P[i], w = P[(i + 1) % n];
        if (b_on_seg(s.a, s.b, u)) ts.push_back({dot128(d, u - s.a), norm2_128(d)});
        Int128 den = cross128(d, w - u), tn = cross128(u - s.a, w - u), un = cross128(u - s.a, d);
        if (den == 0) continue;
        if (den < 0) den = -den, tn = -tn, un = -un;
        if (0 <= tn && tn <= den && 0 <= un && un <= den) ts.push_back({tn, den});
    }
    std::sort(ts.begin(), ts.end(), [](const BFrac& a, const BFrac& b) { return a.n * b.d < b.n * a.d; });
    std::vector<BFrac> all;
    for (size_t k = 0; k < ts.size(); ++k) {
        all.push_back(ts[k]);
        if (k + 1 < ts.size() && ts[k].n * ts[k + 1].d != ts[k + 1].n * ts[k].d)
            all.push_back({ts[k].n * ts[k + 1].d + ts[k + 1].n * ts[k].d, 2 * ts[k].d * ts[k + 1].d});
    }
    bool every = true, some = false;
    for (const BFrac& t : all) {
        const Int128 g = geo_internal::gcd128(t.n, t.d);
        const bool in = b_where_at(P, s, {t.n / g, t.d / g}) != OUTSIDE;
        every = every && in, some = some || in;
    }
    return {every, some};
}
static Real b_seg_dist(const Segment& s, const Segment& t) {
    if (b_meet(s.a, s.b, t.a, t.b)) return 0;
    auto pd = [](const Segment& g, const Point& p) {
        const Point d = g.b - g.a;
        const Int128 tt = dot128(d, p - g.a), nn = norm2_128(d);
        if (nn == 0 || tt <= 0) return std::sqrt((Real)norm2_128(p - g.a));
        if (tt >= nn) return std::sqrt((Real)norm2_128(p - g.b));
        return std::fabs((Real)cross128(d, p - g.a)) / std::sqrt((Real)nn);
    };
    return std::min({pd(s, t.a), pd(s, t.b), pd(t, s.a), pd(t, s.b)});
}
static Segment b_edge(const std::vector<Point>& P, int i) { return Segment(P[i], P[(i + 1) % P.size()]); }
static std::vector<Point> rnd_simple(std::mt19937_64& rng, int lo, int hi, long long C) {
    while (true) {
        const std::vector<Point> P = rnd_star(rng, (int)rnd(rng, lo, hi), C);
        if (b_simple(P)) return P;
    }
}

static void test_relations(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 3 : 5;
        const std::vector<Point> P = rnd_simple(rng, 3, 9, C);
        const int n = (int)P.size();
        // 線分
        for (int q = 0; q < 6; ++q) {
            Segment s(rndp(rng, C + 1), rndp(rng, C + 1));
            if (q == 1) s = Segment(P[rng() % n], P[rng() % n]);  // 頂点どうし（対角線・辺）
            if (q == 2) s.a = P[rng() % n];
            if (q == 3) s.b = s.a;
            const auto [all, some] = b_seg_in(P, s);
            CHECK(contains(P, s) == all && intersects(P, s) == some);
            const Point p = rndp(rng, C + 1);
            const Real dp = b_contains(P, p) != OUTSIDE ? 0 : [&] {
                Real r = 1e30L;
                for (int i = 0; i < n; ++i) r = std::min(r, b_seg_dist(b_edge(P, i), Segment(p, p)));
                return r;
            }();
            CHECK(near(dist(P, p), dp, TOL));
            Real bd = 1e30L, ds = 1e30L;
            for (int i = 0; i < n; ++i) bd = std::min(bd, b_seg_dist(b_edge(P, i), Segment(p, p))), ds = std::min(ds, b_seg_dist(b_edge(P, i), s));
            CHECK(near(boundary_dist(P, p), bd, TOL));
            CHECK(near(dist(P, s), some ? 0 : ds, TOL));
            // 直線・半直線（座標が小さいので、十分に長い線分の辺との交差で判定できる）
            if (s.a != s.b) {
                const Line l(s.a, s.b);
                const Ray r(s.a, s.b);
                const Segment ls(s.a - s.dir() * 1000000, s.a + s.dir() * 1000000), rs(s.a, s.a + s.dir() * 1000000);
                bool el = false, er = b_contains(P, r.a) != OUTSIDE;
                for (int i = 0; i < n; ++i) el = el || b_meet(ls.a, ls.b, P[i], P[(i + 1) % n]), er = er || b_meet(rs.a, rs.b, P[i], P[(i + 1) % n]);
                CHECK(intersects(P, l) == el && intersects(P, r) == er);
            }
        }
        // 多角形どうし
        std::vector<Point> Q = rnd_simple(rng, 3, 6, C);
        if (it % 3 == 0) {  // 小さく縮めて中に入りやすくする
            for (Point& p : Q) p = Point(p.x / 2, p.y / 2);
            if (!b_simple(Q)) continue;
        }
        bool all = true, some = b_contains(Q, P[0]) != OUTSIDE;
        Real dq = 1e30L;
        for (int j = 0; j < (int)Q.size(); ++j) {
            const auto [a, b] = b_seg_in(P, b_edge(Q, j));
            all = all && a, some = some || b;
            for (int i = 0; i < n; ++i) dq = std::min(dq, b_seg_dist(b_edge(P, i), b_edge(Q, j)));
        }
        CHECK(contains(P, Q) == all && intersects(P, Q) == some && intersects(Q, P) == some);
        CHECK(near(dist(P, Q), some ? 0 : dq, TOL) && near(dist(Q, P), dist(P, Q), TOL));
        if (all) CHECK(area(Q) <= area(P));
        CHECK(contains(P, P));
    }
    // 大きい座標：小さい図形を K 倍して平行移動しても結果が変わらない（128 ビットの中点の判定など）
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long C = 5, K = rnd(rng, 1, 150000000);
        const Point sh = rndp(rng, 200000000);
        auto tf = [&](const Point& p) { return p * K + sh; };
        const std::vector<Point> P = rnd_simple(rng, 3, 9, C);
        std::vector<Point> PB;
        for (const Point& p : P) PB.push_back(tf(p));
        for (int q = 0; q < 6; ++q) {
            Segment s(rndp(rng, C + 1), rndp(rng, C + 1));
            if (q < 2) s = Segment(P[rng() % P.size()], P[rng() % P.size()]);
            const Segment sb(tf(s.a), tf(s.b));
            CHECK(contains(PB, sb) == contains(P, s) && intersects(PB, sb) == intersects(P, s));
            if (s.a != s.b) CHECK(intersects(PB, Line(sb.a, sb.b)) == intersects(P, Line(s.a, s.b)) && intersects(PB, Ray(sb.a, sb.b)) == intersects(P, Ray(s.a, s.b)));
            CHECK(near(dist(PB, sb), dist(P, s) * K, TOL * 16));
        }
    }
}

static void test_prefix(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 1000000000LL : 6;
        std::vector<Point> P = rndps(rng, (int)rnd(rng, 1, 9), C);
        if (it % 3 == 0) {
            const std::vector<Point> S = rnd_star(rng, (int)rnd(rng, 3, 9), C);
            if (!S.empty()) P = S;
        }
        const int n = (int)P.size();
        const PolygonPrefix pp(P);
        for (int q = 0; q < 10; ++q) {
            const int i = (int)rnd(rng, -20, 20), j = (int)rnd(rng, -20, 20);
            std::vector<Point> sub;
            const int a = ((i % n) + n) % n, cnt = (((j - i) % n) + n) % n + 1;
            for (int k = 0; k < cnt; ++k) sub.push_back(P[(a + k) % n]);
            const Int128 e = (int)sub.size() < 3 ? 0 : b_area2(sub);
            CHECK(pp.area2(i, j) == e && near(pp.area(i, j), std::fabs((Real)e) / 2, TOL));
            const PointD g = pp.centroid(i, j), h = centroid(sub);
            const Real sc = std::max<Real>(1, (Real)C);
            CHECK(std::fabs(g.x - h.x) <= TOL * sc * 64 && std::fabs(g.y - h.y) <= TOL * sc * 64);
        }
        CHECK(pp.area2(0, n - 1) == (n < 3 ? 0 : b_area2(P)) && pp.area2(3, 2) == pp.area2(0, n - 1));
    }
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_area_contains(rng, cases);
    test_relations(rng, cases);
    test_prefix(rng, cases);
    test_convex(rng, cases);
    test_triangle(rng, cases);
    test_real(rng, cases);
    std::cout << "[polygon] random_test: OK (" << cases << " cases)\n";
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

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[polygon] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    // 凸多角形（N = 2×10^5）：向きの違う辺ベクトルとその逆を偏角の順に並べて足していく
    const int n = 200000;
    std::vector<Point> P = make_convex_polygon(rng, n, 5000);
    auto t0 = std::chrono::steady_clock::now();
    Coord a2 = 0;
    for (int k = 0; k < 10; ++k) a2 ^= area2(P);
    const long long t_area = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    int in = 0;
    for (int k = 0; k < 100; ++k) in += contains(P, rndp(rng, 1000000000LL)) == INSIDE;
    const long long t_cont = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const bool cv = is_convex(P);
    const long long t_conv = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const PointD g = centroid(P);
    const Coord I = interior_lattice_points(P);
    const long long t_cent = ms_since(t0);
    std::vector<Point> Q(P.begin(), P.begin() + 2000);
    t0 = std::chrono::steady_clock::now();
    const bool sp = is_simple(Q);
    const long long t_simple = ms_since(t0);
    std::cout << "[polygon] benchmark: N = " << P.size() << "  area2 × 10 " << t_area << " ms, contains × 100 " << t_cont << " ms, is_convex "
              << t_conv << " ms, centroid + interior_lattice_points " << t_cent << " ms, is_simple(N = 2000) " << t_simple << " ms  (check "
              << (a2 & 7) + in + cv + sp + (long long)g.x % 7 + I % 7 << ")\n";
}

}  // namespace pg_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace pg_demo {

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
#define SHOW(expr, comment) pg_demo::show(#expr, pg_demo::to_s(expr), comment)
#define RUN(...)                   \
    pg_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("多角形・三角形 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  area2(poly) / area(poly) / is_ccw(poly) / make_ccw(poly)        面積の 2 倍（符号付き・整数）/ 面積 / 向き\n"
              << "  contains(poly, p)                                              OUTSIDE / ON_BORDER / INSIDE（凸でなくてよい）\n"
              << "  is_convex(poly, strict) / is_simple(poly)                      凸か / 辺が交わらないか\n"
              << "  perimeter(poly) / centroid(poly)                               周の長さ / 重心\n"
              << "  boundary_lattice_points(poly) / interior_lattice_points(poly)  周の上 / 内部の格子点の数（ピックの定理）\n"
              << "  triangle_area2(a, b, c) / triangle_kind(a, b, c) / in_triangle(a, b, c, p)\n"
              << "  contains(poly, s) / contains(P, Q) / intersects(poly, s or l or r) / intersects(P, Q)   位置関係（周を含む）\n"
              << "  dist(poly, p or s) / dist(P, Q) / boundary_dist(poly, p)          距離（中なら 0）/ 周までの距離\n"
              << "  PolygonPrefix pp(poly): pp.area2(i, j) / pp.area(i, j) / pp.centroid(i, j)   頂点 i..j の部分を O(1) で\n";

    api("【1】 area2(poly) / area(poly) / is_ccw(poly) / make_ccw(poly)", "Coord / Real / bool / vector<Point>",
        "符号付き面積の 2 倍（反時計回りなら正。整数のまま誤差なし）/ 面積 / 反時計回りか / 反時計回りにそろえる", "O(N)");
    RUN(std::vector<Point> poly = {{0, 0}, {4, 0}, {4, 3}, {2, 1}, {0, 3}});
    SHOW(area2(poly), "凹んだ五角形。面積の 2 倍");
    SHOW(area(poly), "面積");
    SHOW(is_ccw(poly), "反時計回り");
    SHOW(make_ccw(std::vector<Point>{{0, 0}, {0, 2}, {2, 0}}), "時計回りの三角形を反時計回りに");
    note("面積の 2 倍で比べる・出力する（area2 が奇数なら答えは .5）と誤差がない");

    api("【2】 contains(poly, p)", "Where（OUTSIDE / ON_BORDER / INSIDE）", "点が多角形の内側か。凸でなくてよい（凸でクエリが多いなら 04_convex_hull の convex_contains）", "O(N)");
    SHOW(contains(poly, Point(1, 1)) == INSIDE, "内側");
    SHOW(contains(poly, Point(2, 2)) == OUTSIDE, "凹んだところは外");
    SHOW(contains(poly, Point(3, 2)) == ON_BORDER, "辺 (4,3)-(2,1) の上");

    api("【3】 is_convex(poly, strict = false) / is_simple(poly)", "bool",
        "凸か（strict なら 180 度の頂点を許さない。向きはどちらでもよい）/ 辺が交わらない単純多角形か", "O(N) / O(N^2)");
    SHOW(is_convex(poly), "凹んでいる");
    SHOW(is_convex(std::vector<Point>{{0, 0}, {2, 0}, {4, 0}, {4, 4}}), "180 度の頂点 (2,0) があっても凸");
    SHOW(is_convex(std::vector<Point>{{0, 0}, {2, 0}, {4, 0}, {4, 4}}, true), "strict なら false");
    SHOW(is_simple(std::vector<Point>{{0, 0}, {2, 2}, {2, 0}, {0, 2}}), "砂時計の形（辺が交わる）");
    note("ABC266 C（凸四角形か）: is_convex({A, B, C, D}, true)。四角形が反時計回りで与えられるなら ccw がすべて +1 でもよい");

    api("【4】 perimeter(poly) / centroid(poly)", "Real / PointD", "周の長さ / 重心（面積の重心。面積 0 なら頂点の平均）", "O(N)");
    SHOW(perimeter(std::vector<Point>{{0, 0}, {3, 0}, {3, 4}}), "3 + 4 + 5");
    SHOW(centroid(std::vector<Point>{{0, 0}, {6, 0}, {0, 6}}), "三角形の重心 (2, 2)");

    api("【5】 boundary_lattice_points(poly) / interior_lattice_points(poly)", "Coord",
        "周の上の格子点の数（Σ gcd(|dx|, |dy|)）/ 内部の格子点の数（ピックの定理 I = (2S - B + 2) / 2）", "O(N log)");
    SHOW(boundary_lattice_points(poly), "周の上");
    SHOW(interior_lattice_points(poly), "内部");
    note("典型90 041（凸包の中の格子点）: auto h = convex_hull(ps);（04_convex_hull）→ interior + boundary - N");

    api("【6】 triangle_kind(a, b, c) / triangle_area2(a, b, c) / in_triangle(a, b, c, p)", "int / Coord / Where",
        "三角形の種類（誤差なし）/ 符号付き面積の 2 倍 / 三角形の内外", "O(1)");
    SHOW(triangle_kind(Point(0, 0), Point(4, 0), Point(0, 3)) == TRIANGLE_RIGHT, "直角三角形（ABC362 B）");
    SHOW(triangle_kind(Point(0, 0), Point(4, 0), Point(-1, 3)) == TRIANGLE_OBTUSE, "鈍角");
    SHOW(triangle_kind(Point(0, 0), Point(4, 0), Point(2, 3)) == TRIANGLE_ACUTE, "鋭角");
    SHOW(triangle_kind(Point(0, 0), Point(1, 1), Point(2, 2)) == TRIANGLE_DEGENERATE, "一直線");
    SHOW(in_triangle(Point(0, 0), Point(4, 0), Point(0, 4), Point(1, 1)) == INSIDE, "内側");
    note("値の名前は TRIANGLE_ACUTE / TRIANGLE_RIGHT / TRIANGLE_OBTUSE / TRIANGLE_DEGENERATE（RIGHT などとぶつからない）");

    api("【7】 contains(poly, s) / contains(P, Q) / intersects(poly, s / l / r) / intersects(P, Q)", "bool",
        "多角形（周を含む閉じた領域）と線分・直線・半直線・多角形の位置関係。凸でなくてよい（単純多角形）。整数なら厳密", "O(N log N) / O(N M log N) / O(N) / O(N M)");
    RUN(std::vector<Point> L = {Point(0, 0), Point(4, 0), Point(4, 2), Point(2, 2), Point(2, 4), Point(0, 4)});
    note("L 字形（凹んだ角が (2, 2)）");
    SHOW(contains(L, Segment(1, 1, 3, 1)), "中を通る");
    SHOW(contains(L, Segment(1, 3, 3, 1)), "(2, 2) の角をかすめる（周の上を通るのは OK）");
    SHOW(contains(L, Segment(1, 3, 3, 2)), "(1,3)→(3,2)：凹んだ部分（外）を通る");
    SHOW(contains(L, Segment(2, 4, 2, 2)), "周の上だけを通る線分も中");
    SHOW(contains(L, Segment(3, 3, 3, 3)), "外の点（長さ 0 の線分）");
    SHOW(intersects(L, Segment(3, 3, 5, 5)), "外の線分で共有点なし");
    SHOW(intersects(L, Line(3, 3, 4, 3)), "直線 y = 3 は L の左の部分を通る");
    SHOW(intersects(L, Ray(Point(3, 3), Point(4, 3))), "右向きの半直線は当たらない");
    SHOW(contains(L, std::vector<Point>{Point(0, 0), Point(2, 0), Point(2, 2), Point(0, 2)}), "左下の正方形は L に入る");
    SHOW(intersects(L, std::vector<Point>{Point(2, 2), Point(3, 3), Point(2, 3)}), "(2, 2) で接する三角形");
    note("凸多角形どうしの共通部分そのもの・面積は 04_convex_polygon の convex_intersection、円と多角形は 05_circle_area");

    api("【8】 dist(poly, p) / dist(poly, s) / dist(P, Q) / boundary_dist(poly, p)", "Real",
        "距離（多角形は中身を含む。中にある・共有点があれば 0）/ 周までの距離（中の点でも周まで）", "O(N) / O(N M)");
    SHOW(dist(L, Point(3, 3)), "凹んだ部分の点から L まで（1）");
    SHOW(boundary_dist(L, Point(1, 1)), "中の点から周まで（ARC042 B）");
    SHOW(dist(L, Segment(3, 3, 5, 5)), "線分と L");

    api("【9】 PolygonPrefix pp(poly) / pp.area2(i, j) / pp.area(i, j) / pp.centroid(i, j)", "Int128 / Real / PointD",
        "連続する頂点 v[i], v[i+1], ..., v[j]（番号は n で割った余り）で囲む部分の符号付き面積の 2 倍・面積・重心", "前処理 O(N)、O(1)");
    RUN(PolygonPrefix pp(std::vector<Point>{Point(0, 0), Point(4, 0), Point(4, 4), Point(0, 4)}));
    SHOW(int128_to_string(pp.area2(0, 2)), "三角形 (0,0) (4,0) (4,4) の面積の 2 倍");
    SHOW(int128_to_string(pp.area2(3, 1)), "番号は一周してよい：(0,4) (0,0) (4,0)");
    SHOW(pp.centroid(0, 3), "全体の重心");
    note("ABC250 F（凸多角形を対角線で切って面積 S/4 に近づける）: 尺取りで j を進めながら pp.area2(i, j) を O(1) で求める");
}

}  // namespace pg_demo


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
    pg_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pg_test::benchmark();    // 速度計測
    pg_demo::demo();         // 使い方デモ
    return 0;
}

