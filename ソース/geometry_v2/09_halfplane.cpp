// =============================================================================
//  半平面交差（a x + b y ≤ c をすべて満たす領域・有界か）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 09 半平面交差
//  依存       : 00_int_common, 00_geo_base, 02_line, 03_polygon（同梱済み。このファイルだけ貼れば動く。結果に area などをそのまま使える）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - 半平面 HalfPlane h(a, b, c) は a x + b y ≤ c（境界の直線の上を含む）。HalfPlane(l) は直線 l の左側
//    - 整数の半平面（HalfPlane / Line）なら、どの交点を残すかの判定は 128 ビットで誤差なし
//      （|a|, |b| ≤ 2×10^9、|c| ≤ 4×10^18 まで。座標 10^9 までの点で作った Line はこの範囲に入る）
//    - 答えの頂点だけ実数（PointD）。交点を整数に丸めない
//    - 実数の半平面（HalfPlaneD / LineD）は相対誤差 GEO_EPS で判定する。答えの頂点は各半平面を誤差の範囲で満たす
//      （面積が 0 に近い答え（線分・点）は空になることがある。面積 0 を区別したいときは整数の版を使う）
//
//  機能一覧                                                                    計算量
//    HalfPlane h(a, b, c) / HalfPlane(l) / HalfPlaneD   a x + b y ≤ c / 直線 l の左側。h.dir()（境界の向き。左側が領域）
//    contains(h, p)                    OUTSIDE / ON_BORDER / INSIDE                                 O(1)
//    halfplane_intersection(hs, box)   共通部分の凸多角形（vector<PointD>、反時計回り、(x, y) 最小の頂点（Real に直した値で比べる）から）  O(N log N)
//                                      空なら空。面積 0 の線分なら 2 点、1 点なら 1 点。同じ点・一直線の頂点は含まない
//                                      有界でない（境界の向きが 180 度以上あく）ときだけ |x|, |y| ≤ box の正方形で切る（box の既定は 2×10^18）
//                                      有界なら box は使わない（ほぼ平行な境界で頂点が 10^18 より遠くなっても正しい）
//    halfplane_intersection(lines, box)  各直線の左側の共通部分（vector<Line> / vector<LineD> のまま渡せる）
//    halfplane_unbounded(hs)           共通部分が空でないとき、有界でない（どこまでも続く）か      O(N log N)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    hp_test::random_test()  … 全部の交点から有理数で作った凸包と突き合わせるランダムテスト
//    hp_test::benchmark()    … 速度計測
//    hp_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_HALFPLANE
#define GG_HALFPLANE
// ---- 半平面交差（幾何第二世代）  a x + b y ≤ c の共通部分。整数の入力なら判定は誤差なし、頂点だけ実数 ----
template <class T>
struct HalfPlaneT {  // a x + b y ≤ c（境界の直線の上を含む）
    T a = 0, b = 0, c = 0;
    constexpr HalfPlaneT() {}
    template <class A, class B, class C, std::enable_if_t<geo_internal::coord_args_ok<T, A, B, C>, int> = 0>
    constexpr HalfPlaneT(A a_, B b_, C c_) : a(static_cast<T>(a_)), b(static_cast<T>(b_)), c(static_cast<T>(c_)) {}
    constexpr explicit HalfPlaneT(const LineT<T>& l)  // 直線 l（a → b）の左側
        : a(l.b.y - l.a.y), b(l.a.x - l.b.x), c((l.b.y - l.a.y) * l.a.x + (l.a.x - l.b.x) * l.a.y) {}
    constexpr PointT<T> dir() const { return PointT<T>(-b, a); }  // 境界の直線の向き（左側が領域）
};
using HalfPlane = HalfPlaneT<Coord>;
using HalfPlaneD = HalfPlaneT<Real>;
template <class T>
std::ostream& operator<<(std::ostream& os, const HalfPlaneT<T>& h) { return os << "HalfPlane[" << h.a << " x + " << h.b << " y <= " << h.c << "]"; }

template <class T>
Where contains(const HalfPlaneT<T>& h, const PointT<T>& p) {  // OUTSIDE / ON_BORDER（境界の直線の上）/ INSIDE
    int s;
    if constexpr (std::is_integral_v<T>) {
        const Int128 v = (Int128)h.a * p.x + (Int128)h.b * p.y - h.c;
        s = (v > 0) - (v < 0);
    } else {
        s = real_sign(h.a * p.x + h.b * p.y - h.c, std::fabs(h.a * p.x) + std::fabs(h.b * p.y) + std::fabs(h.c));
    }
    return s < 0 ? INSIDE : (s == 0 ? ON_BORDER : OUTSIDE);
}

namespace geo_internal {
template <class T>
using hp_wide = std::conditional_t<std::is_integral_v<T>, Int128, Real>;
template <class T>
struct HpCross {  // 境界の直線 2 本の交点 (X / D, Y / D)
    hp_wide<T> X, Y, D;
};
template <class T>
HpCross<T> hp_cross(const HalfPlaneT<T>& i, const HalfPlaneT<T>& j) {
    using W = hp_wide<T>;
    return HpCross<T>{(W)i.c * j.b - (W)j.c * i.b, (W)i.a * j.c - (W)j.a * i.c, (W)i.a * j.b - (W)j.a * i.b};
}
// 交点 q が半平面 k の外なら 1、境界の上なら 0、内側なら -1
template <class T>
int hp_side(const HalfPlaneT<T>& k, const HpCross<T>& q) {
    using W = hp_wide<T>;
    const W u = (W)k.a * q.X, v = (W)k.b * q.Y, w = (W)k.c * q.D;
    int s;
    if constexpr (std::is_integral_v<T>) {
        const W t = u + v - w;
        s = (t > 0) - (t < 0);
    } else {
        s = real_sign(u + v - w, std::fabs(u) + std::fabs(v) + std::fabs(w));
    }
    return q.D > 0 ? s : -s;
}
// 半平面 h を足すとき、交点 q を作る境界を取り除くか。整数は q が h の外（境界の上は残す。面積 0 の答えも正しく出る）
// 実数は q が h の外か境界の上（誤差の範囲）でも取り除く（ほぼ平行な境界で頂点が逆戻りしないように）
template <class T>
bool hp_pop(const HalfPlaneT<T>& h, const HpCross<T>& q) {
    if constexpr (std::is_integral_v<T>) return hp_side(h, q) > 0;
    else return hp_side(h, q) >= 0;
}
template <class T>
int hp_dir_cross(const HalfPlaneT<T>& i, const HalfPlaneT<T>& j) {  // 境界の向きの外積の符号
    if constexpr (std::is_integral_v<T>) {
        const Int128 c = (Int128)i.a * j.b - (Int128)i.b * j.a;  // cross((-bi, ai), (-bj, aj))
        return (c > 0) - (c < 0);
    } else {
        return sgn_cross(i.dir(), j.dir());
    }
}
template <class T>
bool hp_same_dir(const HalfPlaneT<T>& i, const HalfPlaneT<T>& j) {  // 境界の向きが同じ（平行で同じ向き）
    if constexpr (std::is_integral_v<T>) return hp_dir_cross(i, j) == 0 && (Int128)i.a * j.a + (Int128)i.b * j.b > 0;
    else return hp_dir_cross(i, j) == 0 && sgn_dot(i.dir(), j.dir()) > 0;
}
template <class T>
bool hp_stricter(const HalfPlaneT<T>& j, const HalfPlaneT<T>& i) {  // 同じ向きの 2 つで、j の方が狭い（真に）
    if constexpr (std::is_integral_v<T>) {
        const bool useA = i.a != 0;
        const Int128 ai = useA ? (i.a < 0 ? -(Int128)i.a : (Int128)i.a) : (i.b < 0 ? -(Int128)i.b : (Int128)i.b);
        const Int128 aj = useA ? (j.a < 0 ? -(Int128)j.a : (Int128)j.a) : (j.b < 0 ? -(Int128)j.b : (Int128)j.b);
        return (Int128)j.c * ai < (Int128)i.c * aj;
    } else {
        return j.c / std::sqrt(j.a * j.a + j.b * j.b) < i.c / std::sqrt(i.a * i.a + i.b * i.b);
    }
}
// 境界の向きの偏角の順に並べる（整数なら外積で厳密に、実数なら atan2 で）
template <class T>
void hp_sort(std::vector<HalfPlaneT<T>>& v) {
    if constexpr (std::is_integral_v<T>) {
        auto half = [](const HalfPlaneT<T>& h) { const Point d = h.dir(); return (d.y < 0 || (d.y == 0 && d.x < 0)) ? 1 : 0; };
        std::sort(v.begin(), v.end(), [&](const HalfPlaneT<T>& i, const HalfPlaneT<T>& j) {
            const int hi = half(i), hj = half(j);
            return hi != hj ? hi < hj : hp_dir_cross(i, j) > 0;
        });
    } else {
        std::vector<std::pair<Real, int>> key(v.size());
        for (size_t k = 0; k < v.size(); ++k) key[k] = {std::atan2(v[k].a, -v[k].b), (int)k};
        std::sort(key.begin(), key.end());
        std::vector<HalfPlaneT<T>> s(v.size());
        for (size_t k = 0; k < v.size(); ++k) s[k] = v[key[k].second];
        v.swap(s);
    }
}
// v（偏角の順）の境界の向きが 180 度以上あくか（あくなら、共通部分は空でなければ有界でない）
template <class T>
bool hp_open(const std::vector<HalfPlaneT<T>>& v) {
    std::vector<HalfPlaneT<T>> u;  // 向きの違うものだけ
    for (const auto& h : v)
        if (u.empty() || !hp_same_dir(u.back(), h)) u.push_back(h);
    while (u.size() >= 2 && hp_same_dir(u.back(), u.front())) u.pop_back();
    const int n = (int)u.size();
    if (n <= 2) return true;
    for (int i = 0; i < n; ++i) {
        const auto& p = u[i];
        const auto& q = u[(i + 1) % n];
        const int c = hp_dir_cross(p, q);
        if (c < 0 || (c == 0 && !hp_same_dir(p, q))) return true;  // 次の向きまで 180 度以上
    }
    return false;
}
template <class T>
std::vector<PointD> hp_run(const std::vector<HalfPlaneT<T>>& in, T box) {
    std::vector<HalfPlaneT<T>> v;
    v.reserve(in.size() + 4);
    for (const auto& h : in) {
        if (h.a == 0 && h.b == 0) {  // 0 ≤ c
            if (h.c < 0) return {};
            continue;
        }
        v.push_back(h);
    }
    hp_sort(v);
    if (hp_open(v)) {  // 有界でないかもしれないときだけ箱で切る
        v.push_back(HalfPlaneT<T>(1, 0, box)), v.push_back(HalfPlaneT<T>(0, 1, box));
        v.push_back(HalfPlaneT<T>(-1, 0, box)), v.push_back(HalfPlaneT<T>(0, -1, box));
        hp_sort(v);
    }
    // 同じ向きは最も狭いものだけ残す
    std::vector<HalfPlaneT<T>> u;
    for (const auto& h : v) {
        if (!u.empty() && hp_same_dir(u.back(), h)) {
            if (hp_stricter(h, u.back())) u.back() = h;
            continue;
        }
        u.push_back(h);
    }
    if (u.size() >= 2 && hp_same_dir(u.back(), u.front())) {  // 実数で偏角が 0 と 2π に分かれた同じ向き
        if (hp_stricter(u.back(), u.front())) u.front() = u.back();
        u.pop_back();
    }
    const int n = (int)u.size();
    std::vector<HalfPlaneT<T>> dq(n);
    int fr = 0, bk = 0;
    for (const auto& h : u) {
        while (bk - fr >= 2 && hp_pop(h, hp_cross(dq[bk - 2], dq[bk - 1]))) --bk;
        while (bk - fr >= 2 && hp_pop(h, hp_cross(dq[fr], dq[fr + 1]))) ++fr;
        if (bk - fr >= 1 && hp_dir_cross(dq[bk - 1], h) <= 0) return {};  // 180 度以上曲がる → 空
        dq[bk++] = h;
    }
    while (bk - fr >= 3 && hp_pop(dq[fr], hp_cross(dq[bk - 2], dq[bk - 1]))) --bk;
    while (bk - fr >= 3 && hp_pop(dq[bk - 1], hp_cross(dq[fr], dq[fr + 1]))) ++fr;
    const int m = bk - fr;
    if (m < 3 || hp_dir_cross(dq[bk - 1], dq[fr]) <= 0) return {};
    // 頂点 q_i = (i 番目と i + 1 番目の境界の交点)。q_i が i + 2 番目の境界の上なら q_{i+1} と同じ点
    std::vector<PointD> res;
    for (int i = 0; i < m; ++i) {
        const HpCross<T> q = hp_cross(dq[fr + i], dq[fr + (i + 1) % m]);
        if constexpr (std::is_integral_v<T>)  // 整数：q が次の次の境界の上なら次の頂点と同じ点（誤差なし）
            if (hp_side(dq[fr + (i + 2) % m], q) == 0 && !(i == m - 1 && res.empty())) continue;
        res.push_back(PointD((Real)q.X / (Real)q.D, (Real)q.Y / (Real)q.D));
    }
    if constexpr (!std::is_integral_v<T>) {  // 実数：誤差の範囲で同じ点を除く
        std::vector<PointD> r2;
        for (const PointD& p : res)
            if (r2.empty() || !real_eq(r2.back(), p)) r2.push_back(p);
        while (r2.size() >= 2 && real_eq(r2.back(), r2.front())) r2.pop_back();
        res.swap(r2);
    }
    std::rotate(res.begin(), std::min_element(res.begin(), res.end()), res.end());
    return res;
}
}  // namespace geo_internal

// 半平面の共通部分（反時計回りの凸多角形。(x, y) 最小の頂点から）。空なら空、線分なら 2 点、1 点なら 1 点
// 有界でないときは |x|, |y| ≤ box の正方形で切る
template <class T>
std::vector<PointD> halfplane_intersection(const std::vector<HalfPlaneT<T>>& hs, T box = static_cast<T>(2000000000000000000LL)) {
    return geo_internal::hp_run(hs, box);
}
template <class T>
std::vector<PointD> halfplane_intersection(const std::vector<LineT<T>>& lines, T box = static_cast<T>(2000000000000000000LL)) {  // 各直線の左側
    std::vector<HalfPlaneT<T>> hs;
    hs.reserve(lines.size());
    for (const auto& l : lines) hs.push_back(HalfPlaneT<T>(l));
    return geo_internal::hp_run(hs, box);
}
// 共通部分が空でないとき、有界でないか（境界の向きが半円（180 度以内）に収まるとき）
template <class T>
bool halfplane_unbounded(const std::vector<HalfPlaneT<T>>& hs) {
    std::vector<HalfPlaneT<T>> v;  // 境界の向き（偏角の順。整数なら厳密）
    for (const auto& h : hs)
        if (h.a != 0 || h.b != 0) v.push_back(h);
    geo_internal::hp_sort(v);
    return geo_internal::hp_open(v);
}
template <class T>
bool halfplane_unbounded(const std::vector<LineT<T>>& lines) {
    std::vector<HalfPlaneT<T>> hs;
    for (const auto& l : lines) hs.push_back(HalfPlaneT<T>(l));
    return halfplane_unbounded(hs);
}
#endif  // GG_HALFPLANE
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
namespace hp_test {

#define HP_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[halfplane] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK HP_CHECK

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

// ---- 愚直解：全部の交点のうち全部の半平面を満たすものの凸包（有理数で厳密） ----
struct BQ {  // (x / d, y / d)、d > 0
    Int128 x, y, d;
};
static bool bq_less(const BQ& p, const BQ& q) {
    const Int128 a = p.x * q.d, b = q.x * p.d;
    if (a != b) return a < b;
    return p.y * q.d < q.y * p.d;
}
static bool bq_eq(const BQ& p, const BQ& q) { return p.x * q.d == q.x * p.d && p.y * q.d == q.y * p.d; }
static int bq_ccw(const BQ& a, const BQ& b, const BQ& c) {
    const Int128 ux = b.x * a.d - a.x * b.d, uy = b.y * a.d - a.y * b.d, vx = c.x * a.d - a.x * c.d, vy = c.y * a.d - a.y * c.d;
    const Int128 s = ux * vy - uy * vx;
    return (s > 0) - (s < 0);
}
static std::vector<PointD> brute(std::vector<HalfPlane> hs, Coord box) {
    for (const auto& h : hs)
        if (h.a == 0 && h.b == 0 && h.c < 0) return {};
    std::vector<HalfPlane> v;
    for (const auto& h : hs)
        if (h.a != 0 || h.b != 0) v.push_back(h);
    v.push_back(HalfPlane(1, 0, box)), v.push_back(HalfPlane(-1, 0, box)), v.push_back(HalfPlane(0, 1, box)), v.push_back(HalfPlane(0, -1, box));
    std::vector<BQ> cand;
    for (size_t i = 0; i < v.size(); ++i)
        for (size_t j = i + 1; j < v.size(); ++j) {
            Int128 D = (Int128)v[i].a * v[j].b - (Int128)v[j].a * v[i].b;
            if (D == 0) continue;
            Int128 X = (Int128)v[i].c * v[j].b - (Int128)v[j].c * v[i].b, Y = (Int128)v[i].a * v[j].c - (Int128)v[j].a * v[i].c;
            if (D < 0) D = -D, X = -X, Y = -Y;
            bool ok = true;
            for (const auto& h : v) ok &= (Int128)h.a * X + (Int128)h.b * Y <= (Int128)h.c * D;
            if (ok) cand.push_back(BQ{X, Y, D});
        }
    std::sort(cand.begin(), cand.end(), bq_less);
    std::vector<BQ> u;
    for (const BQ& q : cand)
        if (u.empty() || !bq_eq(u.back(), q)) u.push_back(q);
    std::vector<BQ> h;
    if (u.size() <= 2) h = u;
    else {
        std::vector<BQ> st(2 * u.size());
        int k = 0;
        for (size_t i = 0; i < u.size(); ++i) {
            while (k >= 2 && bq_ccw(st[k - 2], st[k - 1], u[i]) <= 0) --k;
            st[k++] = u[i];
        }
        for (int i = (int)u.size() - 2, t = k + 1; i >= 0; --i) {
            while (k >= t && bq_ccw(st[k - 2], st[k - 1], u[i]) <= 0) --k;
            st[k++] = u[i];
        }
        st.resize(k - 1);
        h = st;
        if (h.size() == 1) h = {u.front(), u.back()};  // 全部が一直線
    }
    std::vector<PointD> r;
    for (const BQ& q : h) r.push_back(PointD((Real)q.x / (Real)q.d, (Real)q.y / (Real)q.d));
    return r;
}
static bool same_poly(const std::vector<PointD>& a, const std::vector<PointD>& b, Real tol) {  // 始点のずれは許す（同じ x の頂点は丸めで順番が入れ替わりうる）
    if (a.size() != b.size()) return false;
    const size_t n = a.size();
    for (size_t s = 0; s < std::max<size_t>(n, 1); ++s) {
        bool ok = true;
        for (size_t i = 0; i < n && ok; ++i) ok = nearp(a[i], b[(i + s) % n], tol);
        if (ok) return true;
    }
    return false;
}
static HalfPlane rnd_hp(std::mt19937_64& rng, long long C) {
    if (rng() % 3 == 0) return HalfPlane(rnd(rng, -3, 3), rnd(rng, -3, 3), rnd(rng, -2 * C, 2 * C));
    Point a = rndp(rng, C), b = rndp(rng, C);
    while (b == a) b = rndp(rng, C);
    return HalfPlane(Line(a, b));
}

static void test_small(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 2 : 4;
        const int n = (int)rnd(rng, 0, it % 4 == 0 ? 10 : 6);
        std::vector<HalfPlane> hs;
        for (int k = 0; k < n; ++k) {
            const int t = (int)(rng() % 6);
            if (t == 0 && !hs.empty()) hs.push_back(hs[rng() % hs.size()]);  // 同じもの
            else if (t == 1 && !hs.empty()) {                                // 逆向き（境界が同じ・平行）
                const HalfPlane h = hs[rng() % hs.size()];
                hs.push_back(HalfPlane(-h.a, -h.b, -h.c + rnd(rng, -1, 1) * (rng() % 2)));
            } else if (t == 2 && !hs.empty()) {  // 同じ向きで定数だけ違う・係数を何倍かしたもの
                const HalfPlane h = hs[rng() % hs.size()];
                const Coord s = rnd(rng, 1, 3);
                hs.push_back(HalfPlane(h.a * s, h.b * s, h.c * s + rnd(rng, -2, 2)));
            } else hs.push_back(rnd_hp(rng, C));
        }
        if (it % 50 == 0) hs.push_back(HalfPlane(0, 0, rnd(rng, -1, 1)));
        const Coord box = it % 2 ? 6 : 50;
        // 有界なら箱は使わない（箱を大きくしても変わらないもの）、有界でなければ箱で切ったもの
        const std::vector<PointD> b1 = brute(hs, 5000), b2 = brute(hs, 9000);
        const bool unb = !same_poly(b1, b2, 1e-12L);
        const std::vector<PointD> exp = unb ? brute(hs, box) : b1, got = halfplane_intersection(hs, box);
        if (!same_poly(exp, got, 1e-12L)) {
            std::cerr << "n=" << hs.size() << " box=" << box << "\n";
            for (auto& h : hs) std::cerr << "  " << h << "\n";
            std::cerr << "exp:";
            for (auto& p : exp) std::cerr << " " << p;
            std::cerr << "\ngot:";
            for (auto& p : got) std::cerr << " " << p;
            std::cerr << "\n";
        }
        CHECK(same_poly(exp, got, 1e-12L));
        // 実数の版（誤差つき）でも同じ（小さい座標なので）
        std::vector<HalfPlaneD> hd;
        for (const auto& h : hs) hd.push_back(HalfPlaneD((Real)h.a, (Real)h.b, (Real)h.c));
        const std::vector<PointD> gd = halfplane_intersection(hd, (Real)box);
        if (exp.size() >= 3) CHECK(same_poly(exp, gd, 1e-9L));
        else CHECK(gd.empty() || same_poly(exp, gd, 1e-9L));  // 面積 0 は実数では空になることがある
        // 内外：領域の頂点はすべての半平面の内側か境界の上
        for (const PointD& p : got)
            for (const auto& h : hd) CHECK(contains(h, p) != OUTSIDE);
        // 有界か：箱の大きさを変えても同じなら有界
        if (!b1.empty()) CHECK(halfplane_unbounded(hs) == unb);
    }
    // 直線の左側として
    for (int it = 0; it < 5000 * TEST_SCALE; ++it, ++cases) {
        std::vector<Line> ls;
        std::vector<HalfPlane> hs;
        for (int k = 0, n = (int)rnd(rng, 0, 6); k < n; ++k) {
            Point a = rndp(rng, 3), b = rndp(rng, 3);
            while (b == a) b = rndp(rng, 3);
            ls.push_back(Line(a, b)), hs.push_back(HalfPlane(Line(a, b)));
        }
        const std::vector<PointD> b1 = brute(hs, 5000);
        const std::vector<PointD> e = same_poly(b1, brute(hs, 9000), 1e-12L) ? b1 : brute(hs, 7);
        CHECK(same_poly(halfplane_intersection(ls, 7LL), e, 1e-12L));
        CHECK(halfplane_unbounded(ls) == halfplane_unbounded(hs));
    }
}

// 大きい座標：凸多角形の辺と、外側を通る余分な直線 → 答えは元の凸多角形
static void test_big(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 1000000000LL : 300;
        std::vector<Point> ps = rndps(rng, (int)rnd(rng, 3, 12), C);
        // 凸包（素朴に：Andrew）
        std::sort(ps.begin(), ps.end());
        ps.erase(std::unique(ps.begin(), ps.end()), ps.end());
        if (ps.size() < 3) continue;
        std::vector<Point> h(2 * ps.size());
        int k = 0;
        for (size_t i = 0; i < ps.size(); ++i) {
            while (k >= 2 && cross128(h[k - 1] - h[k - 2], ps[i] - h[k - 1]) <= 0) --k;
            h[k++] = ps[i];
        }
        for (int i = (int)ps.size() - 2, t = k + 1; i >= 0; --i) {
            while (k >= t && cross128(h[k - 1] - h[k - 2], ps[i] - h[k - 1]) <= 0) --k;
            h[k++] = ps[i];
        }
        h.resize(k - 1);
        if (h.size() < 3) continue;
        std::vector<Line> ls;
        const int m = (int)h.size();
        for (int i = 0; i < m; ++i) ls.push_back(Line(h[i], h[(i + 1) % m]));
        for (int q = 0; q < 6; ++q) {  // 頂点を通る余分な直線（多角形の外側へ向けて）
            const int i = (int)(rng() % m);
            const Point d = h[(i + 1) % m] - h[(i + m - 1) % m];  // 頂点 i の前後を結ぶ向き（この向きの直線は多角形を左に見る）
            ls.push_back(Line(h[i], h[i] + d));
        }
        std::shuffle(ls.begin(), ls.end(), rng);
        const std::vector<PointD> got = halfplane_intersection(ls);
        CHECK((int)got.size() == m);
        for (int i = 0; i < m; ++i) CHECK(nearp(got[i], to_d(h[i]), 1e-15L));
        CHECK(!halfplane_unbounded(ls));
        // 辺 0 の逆向きを足す → 辺 0 の線分だけ。外へ 1 ずらした逆向き → 空
        std::vector<HalfPlane> hs;
        for (const Line& l : ls) hs.push_back(HalfPlane(l));
        const HalfPlane e0(Line(h[0], h[1]));
        hs.push_back(HalfPlane(-e0.a, -e0.b, -e0.c));
        const std::vector<PointD> seg = halfplane_intersection(hs);
        CHECK(seg.size() == 2 && nearp(seg[0], to_d(std::min(h[0], h[1])), 1e-15L) && nearp(seg[1], to_d(std::max(h[0], h[1])), 1e-15L));
        hs.back().c -= 1;  // a x + b y ≥ c + 1
        CHECK(halfplane_intersection(hs).empty());
        hs.back().c += 2;  // a x + b y ≥ c - 1（少し内側まで）→ 細い多角形
        CHECK(halfplane_intersection(hs).size() >= 3);
        // 1 点で接する：辺 0 と辺 m-1 の逆向き → 頂点 0 だけ
        hs.back() = HalfPlane(Line(h[1], h[0]));
        hs.push_back(HalfPlane(Line(h[0], h[m - 1])));
        const std::vector<PointD> pt = halfplane_intersection(hs);
        CHECK(pt.size() == 1 && nearp(pt[0], to_d(h[0]), 1e-15L));
    }
}

// 隣り合う境界の向きの差が 10^-19 ラジアン程度（実数の偏角では並べる順を間違える）でも正しく
// 有界で、頂点が 10^18 より遠い（箱 2×10^18 の外）三角形も正しく
static void test_tiny_angle(int& cases) {
    {
        const std::vector<Line> ls = {Line(Point(-1000000000, 0), Point(1000000000, 1)), Line(Point(999999999, 2), Point(-1000000000, 1)),
                                      Line(Point(1000000000, 0), Point(1000000000, 1))};
        std::vector<HalfPlane> hs;
        for (const Line& l : ls) hs.push_back(HalfPlane(l));
        // 有界な三角形なので、頂点は 3 本の直線どうしの交点 3 つ（brute の有理数は 128 ビットを超えるので使わない）
        std::vector<PointD> exp;
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 3; ++j) {
                const auto q = geo_internal::hp_cross(hs[i], hs[j]);
                exp.push_back(PointD((Real)q.X / (Real)q.D, (Real)q.Y / (Real)q.D));
            }
        const std::vector<PointD> got = halfplane_intersection(ls);
        CHECK(!halfplane_unbounded(ls) && got.size() == 3);
        for (const PointD& e : exp) CHECK(nearp(e, got[0], TOL) || nearp(e, got[1], TOL) || nearp(e, got[2], TOL));
        CHECK(near(got[0].x, -3999999999000000000.0L, TOL) && near(area(got), area(exp), 1e-9L) && area(got) > 1e18L);
        ++cases;
    }
    const std::vector<Point> q = {{-999999999, -999999999}, {0, -1}, {1000000000, 999999998}, {-999999999, 999999998}};
    CHECK(cross128(q[1] - q[0], q[2] - q[1]) == 1);
    std::vector<Line> ls;
    for (int i = 0; i < 4; ++i) ls.push_back(Line(q[i], q[(i + 1) % 4]));
    for (int rep = 0; rep < 4; ++rep, ++cases) {
        std::rotate(ls.begin(), ls.begin() + 1, ls.end());
        CHECK(!halfplane_unbounded(ls));
        const std::vector<PointD> got = halfplane_intersection(ls);
        CHECK(got.size() == 4);
        for (int i = 0; i < 4; ++i) CHECK(nearp(got[i], to_d(q[i]), TOL));
    }
    std::vector<Line> ls3 = {ls[0], ls[1], ls[2]};  // 1 本除くと有界でない
    CHECK(halfplane_unbounded(ls3));
}

// たくさんのほぼ平行な境界（円に接する直線）：実数の版の頂点が各半平面を誤差の範囲で満たし、面積が整数の版と合う
static void test_dense(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 6; ++it, ++cases) {
        const int n = it < 3 ? 3000 : 300;
        const Real R = it % 2 ? 1e9L : 1e5L;
        std::vector<Line> ls;
        for (int i = 0; i < n; ++i) {
            const Real t = rndr(rng, 0, 2 * GEO_PI);
            const Point p((Coord)(R * std::cos(t)), (Coord)(R * std::sin(t)));
            ls.push_back(Line(p, p + rot90(Point((Coord)(R / 1000 * std::cos(t)) + 1, (Coord)(R / 1000 * std::sin(t))))));
        }
        std::vector<LineD> ld;
        for (const auto& l : ls) ld.push_back(to_d(l));
        const std::vector<PointD> r = halfplane_intersection(ls), rd = halfplane_intersection(ld);
        CHECK(r.size() >= 3 && rd.size() >= 3 && near(area(r), area(rd), 1e-12L));
        for (const PointD& p : rd)
            for (const LineD& l : ld) CHECK(cross(l.dir(), p - l.a) / len(l.dir()) >= -R * 1e-8L);  // 距離で R の 1e-8 倍まで
        for (const PointD& p : r)
            for (const LineD& l : ld) CHECK(cross(l.dir(), p - l.a) / len(l.dir()) >= -R * 1e-15L);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_small(rng, cases);
    test_big(rng, cases);
    test_dense(rng, cases);
    test_tiny_angle(cases);
    CHECK(halfplane_intersection(std::vector<HalfPlane>{}, 3LL).size() == 4 && halfplane_unbounded(std::vector<HalfPlane>{}));
    CHECK(halfplane_intersection(std::vector<HalfPlane>{HalfPlane(0, 0, -1)}).empty());
    std::cout << "[halfplane] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[halfplane] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 200000;
    std::vector<Line> ls;  // 半径 10^9 の円に接する直線（向きはばらばら）＋原点の近くを通る直線
    for (int i = 0; i < n; ++i) {
        const Real t = rndr(rng, 0, 2 * GEO_PI);
        const Point p((Coord)(1e9L * std::cos(t)), (Coord)(1e9L * std::sin(t)));
        ls.push_back(Line(p, p + rot90(Point((Coord)(1e6L * std::cos(t)), (Coord)(1e6L * std::sin(t))))));
    }
    auto t0 = std::chrono::steady_clock::now();
    const std::vector<PointD> r = halfplane_intersection(ls);
    const long long t1 = ms_since(t0);
    std::vector<LineD> ld;
    for (const auto& l : ls) ld.push_back(to_d(l));
    t0 = std::chrono::steady_clock::now();
    const std::vector<PointD> rd = halfplane_intersection(ld);
    const long long t2 = ms_since(t0);
    std::cout << "[halfplane] benchmark: N = " << n << "  整数 " << t1 << " ms（頂点 " << r.size() << "）, 実数 " << t2 << " ms（頂点 " << rd.size() << "）\n";
}

}  // namespace hp_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace hp_demo {

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
#define SHOW(expr, comment) hp_demo::show(#expr, hp_demo::to_s(expr), comment)
#define RUN(...)                   \
    hp_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("半平面交差 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  HalfPlane h(a, b, c)  /  HalfPlane(l)  /  HalfPlaneD               a x + b y ≤ c / 直線 l の左側（境界を含む）\n"
              << "  contains(h, p)                                                      OUTSIDE / ON_BORDER / INSIDE\n"
              << "  halfplane_intersection(hs or lines, box = 2e18)                     共通部分の凸多角形（vector<PointD>）\n"
              << "  halfplane_unbounded(hs or lines)                                    （空でないとき）有界でないか\n";

    api("【1】 HalfPlane h(a, b, c)  /  HalfPlane(Line(p, q))", "HalfPlane", "a x + b y ≤ c（境界を含む）。直線から作ると左側", "O(1)");
    RUN(HalfPlane h(1, 1, 4));
    SHOW(contains(h, Point(1, 2)) == INSIDE, "1 + 2 ≤ 4");
    SHOW(contains(h, Point(2, 2)) == ON_BORDER, "境界の上（2 + 2 = 4）");
    SHOW(contains(HalfPlane(Line(Point(0, 0), Point(1, 0))), Point(5, -1)) == OUTSIDE, "x 軸の左側（上側）ではない");
    note("x ≥ 0 は HalfPlane(-1, 0, 0)、y ≥ 1 は HalfPlane(0, -1, -1) のように ≤ にそろえて書く");

    api("【2】 halfplane_intersection(hs)", "vector<PointD>",
        "共通部分（反時計回り、(x, y) 最小の頂点から）。空なら空、線分なら 2 点、1 点なら 1 点", "O(N log N)");
    RUN(std::vector<HalfPlane> hs = {HalfPlane(1, 1, 4), HalfPlane(-1, 0, 0), HalfPlane(0, -1, 0), HalfPlane(1, -1, 2)});
    RUN(std::vector<PointD> poly = halfplane_intersection(hs));
    SHOW(poly, "x + y ≤ 4, x ≥ 0, y ≥ 0, x - y ≤ 2");
    SHOW(area(poly), "03_polygon の area がそのまま使える");
    SHOW(halfplane_intersection(std::vector<HalfPlane>{HalfPlane(1, 0, 1), HalfPlane(-1, 0, -2)}).size(), "x ≤ 1 かつ x ≥ 2 → 空");
    SHOW(halfplane_intersection(std::vector<HalfPlane>{HalfPlane(1, 0, 1), HalfPlane(-1, 0, -1), HalfPlane(0, 1, 3), HalfPlane(0, -1, 0)}),
         "x = 1 の線分（0 ≤ y ≤ 3）");
    note("条件をすべて満たす点があるか → !halfplane_intersection(hs).empty()（境界の上だけでもよい。整数なら誤差なし）");
    note("a x + b y の最大（線形計画）→ 有界なら、どれかの頂点で最大になる");

    api("【3】 halfplane_intersection(lines)  /  halfplane_unbounded(hs)", "vector<PointD> / bool",
        "各直線の左側の共通部分 / 空でないとき、どこまでも続くか（大きな箱で切られた形になる）", "O(N log N)");
    RUN(std::vector<Line> ls = {Line(Point(0, 0), Point(4, 0)), Line(Point(4, 0), Point(0, 4))});
    SHOW(halfplane_unbounded(ls), "2 本では有界にならない");
    RUN(ls.push_back(Line(Point(0, 4), Point(0, 0))));
    SHOW(halfplane_intersection(ls), "3 本で三角形");
    SHOW(halfplane_unbounded(ls), "");
    note("多角形の核（すべての辺を内側から見渡せる点の集まり）= 反時計回りの多角形の各辺の左側の共通部分");

    api("【4】 実数：HalfPlaneD / LineD", "vector<PointD>", "相対誤差 GEO_EPS で判定する。境界を r だけ内側へずらす（c - r·|(a, b)|）などに使う", "O(N log N)");
    note("凸多角形に入る最大の円（半径を二分探索）:");
    RUN(std::vector<Point> tri = {{0, 0}, {6, 0}, {0, 8}});
    Real lo = 0, hi = 10;
    for (int it = 0; it < 100; ++it) {
        const Real r = (lo + hi) / 2;
        std::vector<HalfPlaneD> hd;
        for (int i = 0; i < 3; ++i) {
            HalfPlaneD e(to_d(Line(tri[i], tri[(i + 1) % 3])));
            e.c -= r * std::sqrt(e.a * e.a + e.b * e.b);  // 境界を r だけ内側へ
            hd.push_back(e);
        }
        (halfplane_intersection(hd).empty() ? hi : lo) = r;
    }
    code("Real lo = 0, hi = 10;  // 100 回：各辺を r だけ内側へずらした共通部分が空でなければ r は入る");
    SHOW(lo, "3-4-5 の三角形（2 倍）の内接円の半径 = 2");
}

}  // namespace hp_demo


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
    hp_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    hp_test::benchmark();    // 速度計測
    hp_demo::demo();         // 使い方デモ
    return 0;
}

