// =============================================================================
//  円の共通部分の面積（2 円・円と多角形・円と直線の左側・弓形）と、円と多角形の位置関係   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 05 円
//  依存       : 00_int_common, 00_geo_base, 02_line, 03_polygon, 05_circle（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - 面積は実数（Real）。整数の円・多角形なら、場合分け（交わるか・接するか）と判別式は 128 ビットで厳密
//    - 接する直前のような細い部分でも精度が落ちないように、弓形は角が小さいとき級数で計算する
//
//  機能一覧                                                                    計算量
//    circular_segment_area(r, θ)       弓形の面積（中心角 θ ∈ [0, 2π]。弦と弧で囲まれた部分）      O(1)
//    intersection_area(c1, c2)         2 円の共通部分の面積                                        O(1)
//    intersection_area(c, poly)        円と多角形（凸でなくてよい・向きはどちらでもよい）の共通部分の面積  O(N)
//    circle_cut_area(c, l)             円の、直線 l の左側の部分の面積                            O(1)
//   円と多角形の位置関係（多角形は周を含む閉じた領域・単純多角形。円は円板。整数なら厳密）
//    contains(poly, c)                 円が多角形に入っているか（bool。円周が周に接してもよい。点の contains(poly, p) は Where）  O(N)
//    disk_contains(c, poly)            多角形が円板に入っているか（全部の頂点が円板に入ればよい）   O(N)
//    disk_intersects(c, poly)          円板と多角形が共有点を持つか                                O(N)
//    disk_dist(c, poly)                円板と多角形の距離（共有点があれば 0）                      O(N)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    ca_test::random_test()  … 数値積分（適応シンプソン）と突き合わせるランダムテスト
//    ca_test::benchmark()    … 速度計測
//    ca_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_CIRCLE
#define GG_CIRCLE
// ---- 円（幾何第二世代）  位置関係・交点の数は整数で厳密、交点・接点・中心は実数 ----
enum CircleRelation : int { CIRCLE_CONTAIN = 0, CIRCLE_INTERNAL = 1, CIRCLE_INTERSECT = 2, CIRCLE_EXTERNAL = 3, CIRCLE_SEPARATE = 4, CIRCLE_SAME = 5 };

template <class T>
Where contains(const CircleT<T>& c, const PointT<T>& p) {  // 円の内外（円周の上は ON_BORDER）
    if constexpr (std::is_integral_v<T>) {
        const Int128 d2 = norm2_128(p - c.c), r2 = (Int128)c.r * c.r;
        return d2 < r2 ? INSIDE : (d2 == r2 ? ON_BORDER : OUTSIDE);
    } else {
        const int s = real_cmp(dist(c.c, p), c.r);
        return s < 0 ? INSIDE : (s == 0 ? ON_BORDER : OUTSIDE);
    }
}
// 2 円の位置関係（値は共通接線の数。同じ円は CIRCLE_SAME = 5）
template <class T>
int circle_relation(const CircleT<T>& a, const CircleT<T>& b) {
    if constexpr (std::is_integral_v<T>) {
        const Int128 d2 = norm2_128(b.c - a.c), s = (Int128)a.r + b.r, df = a.r > b.r ? (Int128)a.r - b.r : (Int128)b.r - a.r;
        if (d2 == 0 && a.r == b.r) return CIRCLE_SAME;
        if (d2 > s * s) return CIRCLE_SEPARATE;
        if (d2 == s * s) return CIRCLE_EXTERNAL;
        if (d2 > df * df) return CIRCLE_INTERSECT;
        if (d2 == df * df) return CIRCLE_INTERNAL;
        return CIRCLE_CONTAIN;
    } else {
        const Real d = dist(a.c, b.c), s = a.r + b.r, df = std::fabs(a.r - b.r), sc = std::max(d, s);
        if (real_sign(d, sc) == 0 && real_sign(a.r - b.r, sc) == 0) return CIRCLE_SAME;
        const int c1 = real_sign(d - s, sc);
        if (c1 > 0) return CIRCLE_SEPARATE;
        if (c1 == 0) return CIRCLE_EXTERNAL;
        const int c2 = real_sign(d - df, sc);
        if (c2 > 0) return CIRCLE_INTERSECT;
        if (c2 == 0) return CIRCLE_INTERNAL;
        return CIRCLE_CONTAIN;
    }
}
template <class T>
bool has_common_point(const CircleT<T>& a, const CircleT<T>& b) {  // 円周どうしが共有点を持つか
    const int r = circle_relation(a, b);
    return r == CIRCLE_EXTERNAL || r == CIRCLE_INTERSECT || r == CIRCLE_INTERNAL || r == CIRCLE_SAME;
}
template <class T>
int count_cross(const CircleT<T>& c, const LineT<T>& l) {  // 円周と直線の交点の数（0 / 1（接する）/ 2）
    const PointT<T> d = l.dir();
    if constexpr (std::is_integral_v<T>) {
        const Int128 cr = cross128(d, c.c - l.a), lhs = cr * cr, rhs = (Int128)c.r * c.r * norm2_128(d);
        return lhs > rhs ? 0 : (lhs == rhs ? 1 : 2);
    } else {
        const int s = real_cmp(dist(l, c.c), c.r);
        return s > 0 ? 0 : (s == 0 ? 1 : 2);
    }
}
namespace geo_internal {
// 線分 s.a + t (s.b - s.a)（0 ≤ t ≤ 1）と円周の交点が、二次方程式 A t² + 2B t + C = 0 の解 t1 ≤ t2 のどれか
// （bit0: t1、bit1: t2。接する（重解）なら bit0 だけ）
template <class T>
int circle_segment_mask(const CircleT<T>& c, const SegmentT<T>& s) {
    const PointT<T> d = s.b - s.a, w = s.a - c.c;
    if constexpr (std::is_integral_v<T>) {
        const Int128 A = norm2_128(d), B = dot128(d, w), C = norm2_128(w) - (Int128)c.r * c.r;
        if (A == 0) return C == 0 ? 1 : 0;
        const Int128 D = B * B - A * C, f0 = C, f1 = A + 2 * B + C;
        if (D < 0 || (f0 < 0 && f1 < 0)) return 0;
        if (f0 < 0) return 2;
        if (f1 < 0) return 1;
        if (-B < 0) return f0 == 0 ? 2 : 0;
        if (-B > A) return f1 == 0 ? 1 : 0;
        return D == 0 ? 1 : 3;
    } else {
        const Real A = norm2(d), B = dot(d, w), nw = norm2(w), r2 = c.r * c.r, C = nw - r2;
        const Real sc = std::max({A, std::fabs(B), nw, r2});
        if (real_sign(A, sc) == 0) return real_sign(C, sc) == 0 ? 1 : 0;
        // C = |w|² - r² は桁落ちするので、誤差の目安は |w|² と r² の大きい方。判別式 B² - AC もそれに合わせる
        const Real cs = std::max(nw, r2);
        const int sD = real_sign(B * B - A * C, std::max(B * B, A * cs)), s0 = real_sign(C, cs), s1 = real_sign(A + 2 * B + C, sc);
        if (sD < 0 || (s0 < 0 && s1 < 0)) return 0;
        if (s0 < 0) return 2;
        if (s1 < 0) return 1;
        if (real_sign(-B, sc) < 0) return s0 == 0 ? 2 : 0;
        if (real_sign(-B - A, sc) > 0) return s1 == 0 ? 1 : 0;
        return sD == 0 ? 1 : 3;
    }
}
}  // namespace geo_internal
template <class T>
int count_cross(const CircleT<T>& c, const SegmentT<T>& s) {  // 円周と線分の交点の数（0〜2。端点も含む）
    const int m = geo_internal::circle_segment_mask(c, s);
    return (m & 1) + (m >> 1 & 1);
}

// ---- 交点・接点（実数） ----
template <class T>
std::vector<PointD> cross_points(const CircleT<T>& c, const LineT<T>& l) {  // 円周と直線の交点（直線の向き a → b の順）
    const int k = count_cross(c, l);
    if (k == 0) return {};
    const PointD h = projection(l, c.c);
    if (k == 1) return {h};
    const PointT<T> d = l.dir();
    Real half2;  // 半弦の長さの 2 乗 = r² - (中心と直線の距離)²
    if constexpr (std::is_integral_v<T>) {
        const Int128 cr = cross128(d, c.c - l.a), n2 = norm2_128(d);
        half2 = (Real)((Int128)c.r * c.r * n2 - cr * cr) / (Real)n2;
    } else {
        const Real e = dist(l, c.c);
        half2 = c.r * c.r - e * e;
    }
    const PointD u = unit(d) * std::sqrt(std::max<Real>(0, half2));
    return {h - u, h + u};
}
template <class T>
std::vector<PointD> cross_points(const CircleT<T>& c, const SegmentT<T>& s) {  // 円周と線分の交点（a → b の順）
    const int m = geo_internal::circle_segment_mask(c, s);
    if (m == 0) return {};
    if (s.a == s.b) return {to_d(s.a)};
    const std::vector<PointD> v = cross_points(c, LineT<T>(s.a, s.b));
    if (v.size() <= 1) return v;
    std::vector<PointD> r;
    if (m & 1) r.push_back(v[0]);
    if (m & 2) r.push_back(v[1]);
    return r;
}
// 2 円の交点（c1 → c2 の向きの左側の点が先。接するなら 1 個、同じ円・交わらないなら空）
template <class T>
std::vector<PointD> cross_points(const CircleT<T>& a, const CircleT<T>& b) {
    const int rel = circle_relation(a, b);
    if (rel == CIRCLE_SAME || rel == CIRCLE_SEPARATE || rel == CIRCLE_CONTAIN) return {};
    const PointT<T> v = b.c - a.c;
    PointD base, off;
    if constexpr (std::is_integral_v<T>) {
        const Int128 D2 = norm2_128(v), K = (Int128)a.r * a.r - (Int128)b.r * b.r + D2;  // K = 2 · d · (a.c から交点の弦までの距離)
        base = geo_internal::affine(a.c, v, K, 2 * D2);
        if (rel == CIRCLE_EXTERNAL || rel == CIRCLE_INTERNAL) return {base};
        const Int128 H = 4 * D2 * a.r * a.r - K * K;  // (2 · d · 半弦)²
        off = to_d(rot90(v)) * (std::sqrt((Real)(H < 0 ? 0 : H)) / (2 * (Real)D2));
    } else {
        const Real D2 = norm2(v), K = a.r * a.r - b.r * b.r + D2;
        base = a.c + v * (K / (2 * D2));
        if (rel == CIRCLE_EXTERNAL || rel == CIRCLE_INTERNAL) return {base};
        const Real H = 4 * D2 * a.r * a.r - K * K;
        off = rot90(v) * (std::sqrt(std::max<Real>(0, H)) / (2 * D2));
    }
    return {base + off, base - off};
}
// p から円に引いた接線の接点（p が円周の上なら {p}、内側なら空）
template <class T>
std::vector<PointD> tangent_points(const CircleT<T>& c, const PointT<T>& p) {
    const Where w = contains(c, p);
    if (w == INSIDE) return {};
    if (w == ON_BORDER) return {to_d(p)};
    const PointT<T> v = p - c.c;
    PointD base, off;
    if constexpr (std::is_integral_v<T>) {
        const Int128 D2 = norm2_128(v), r2 = (Int128)c.r * c.r;
        base = geo_internal::affine(c.c, v, r2, D2);
        off = to_d(rot90(v)) * ((Real)c.r * std::sqrt((Real)(D2 - r2)) / (Real)D2);
    } else {
        const Real D2 = norm2(v), r2 = c.r * c.r;
        base = c.c + v * (r2 / D2);
        off = rot90(v) * (c.r * std::sqrt(std::max<Real>(0, D2 - r2)) / D2);
    }
    return {base + off, base - off};
}
// 共通接線。a は c1 との接点、b は c2 との接点（接点が同じなら a + 接線の向き）。外側の接線 → 内側の接線の順
// 本数は circle_relation の値（4 / 3 / 2 / 1 / 0）。同じ円・中心が同じなら空。半径は正
template <class T>
std::vector<LineD> common_tangents(const CircleT<T>& a, const CircleT<T>& b) {
    const int rel = circle_relation(a, b);
    std::vector<LineD> res;
    if (rel == CIRCLE_SAME || rel == CIRCLE_CONTAIN) return res;
    const PointD c1 = to_d(a.c), c2 = to_d(b.c), v = to_d(b.c - a.c), w = rot90(v);
    const Real r1 = (Real)a.r, r2 = (Real)b.r;
    Real D2, e2, i2;  // |v|², |v|² - (r1 - r2)², |v|² - (r1 + r2)²
    if constexpr (std::is_integral_v<T>) {
        const Int128 d2 = norm2_128(b.c - a.c), dm = (Int128)a.r - b.r, dp = (Int128)a.r + b.r;
        D2 = (Real)d2, e2 = (Real)(d2 - dm * dm), i2 = (Real)(d2 - dp * dp);
    } else {
        D2 = norm2(v), e2 = D2 - (r1 - r2) * (r1 - r2), i2 = D2 - (r1 + r2) * (r1 + r2);
    }
    auto add = [&](Real cosn, Real sin2, int sgn, bool inner) {  // 法線 n = (v · cosn ± w · √sin2) / |v|²
        const PointD n = (v * cosn + w * (sgn * std::sqrt(std::max<Real>(0, sin2)))) / D2;
        const PointD q1 = c1 + n * r1, q2 = inner ? c2 - n * r2 : c2 + n * r2;
        res.push_back(real_eq(q1, q2) ? LineD(q1, q1 + rot90(n)) : LineD(q1, q2));
    };
    add(r1 - r2, e2, 1, false);  // 外側の接線
    if (rel != CIRCLE_INTERNAL) add(r1 - r2, e2, -1, false);
    if (rel == CIRCLE_EXTERNAL || rel == CIRCLE_SEPARATE) {  // 内側の接線
        add(r1 + r2, i2, 1, true);
        if (rel == CIRCLE_SEPARATE) add(r1 + r2, i2, -1, true);
    }
    return res;
}

// ---- 三角形の円・2 点を通る円 ----
template <class T>
CircleD circumcircle(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c) {  // 外接円（一直線なら r = -1）
    if (ccw(a, b, c) == 0) return CircleD(to_d(a), -1);
    const PointT<T> u = b - a, v = c - a;
    PointD o;
    if constexpr (std::is_integral_v<T>) {
        const Int128 D = 2 * cross128(u, v), nu = norm2_128(u), nv = norm2_128(v);
        const Int128 X = (Int128)v.y * nu - (Int128)u.y * nv, Y = (Int128)u.x * nv - (Int128)v.x * nu;
        o = PointD((Real)((Int128)a.x * D + X) / (Real)D, (Real)((Int128)a.y * D + Y) / (Real)D);
    } else {
        const Real D = 2 * cross(u, v), nu = norm2(u), nv = norm2(v);
        o = a + PointD(v.y * nu - u.y * nv, u.x * nv - v.x * nu) / D;
    }
    return CircleD(o, dist(o, to_d(a)));
}
namespace geo_internal {
inline int sign_sum_mul3(const Int128 (&a)[3], const long long (&b)[3]) {  // Σ a_i b_i の符号（|a_i| < 2^120。積は 192 ビットまで誤差なし）
    const UInt128 mask = ((UInt128)1 << 64) - 1;
    Int128 H = 0, L = 0;  // 和 = H · 2^64 + L
    for (int i = 0; i < 3; ++i) {
        const Int128 P = (Int128)((UInt128)a[i] & mask) * b[i];
        H += (a[i] >> 64) * b[i] + (P >> 64), L += (Int128)((UInt128)P & mask);
    }
    H += L >> 64, L &= (Int128)mask;
    return H > 0 ? 1 : (H < 0 ? -1 : (L > 0 ? 1 : 0));
}
}  // namespace geo_internal
// 点 p が 3 点 a, b, c を通る円（外接円）の内側か（OUTSIDE / ON_BORDER / INSIDE）。整数なら厳密。a, b, c が一直線なら OUTSIDE
// 4 点が同じ円の上にあるか = in_circumcircle(a, b, c, d) == ON_BORDER
template <class T>
Where in_circumcircle(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c, const PointT<T>& p) {
    const int o = ccw(a, b, c);
    if (o == 0) return OUTSIDE;
    const PointT<T> u = a - p, v = b - p, w = c - p;
    int s;
    if constexpr (std::is_integral_v<T>) {
        const Int128 nu = norm2_128(u), nv = norm2_128(v), nw = norm2_128(w);
        // 行列式 | u.x u.y |u|² ; v.x v.y |v|² ; w.x w.y |w|² | を 1 行目で展開（|u|² ≤ 8×10^18 は long long に収まる）
        const Int128 A[3] = {(Int128)v.y * nw - (Int128)w.y * nv, (Int128)w.x * nv - (Int128)v.x * nw, cross128(v, w)};
        const long long B[3] = {u.x, u.y, (long long)nu};
        s = geo_internal::sign_sum_mul3(A, B);
    } else {
        const Real nu = norm2(u), nv = norm2(v), nw = norm2(w), m = std::max({nu, nv, nw});
        s = real_sign(u.x * (v.y * nw - w.y * nv) - u.y * (v.x * nw - w.x * nv) + nu * cross(v, w), m * m);
    }
    s *= o;
    return s > 0 ? INSIDE : (s == 0 ? ON_BORDER : OUTSIDE);
}
template <class T>
CircleD incircle(const PointT<T>& a, const PointT<T>& b, const PointT<T>& c) {  // 内接円（一直線なら r = -1）
    if (ccw(a, b, c) == 0) return CircleD(to_d(a), -1);
    const Real la = dist(b, c), lb = dist(c, a), lc = dist(a, b), s = la + lb + lc;
    const PointD o = (to_d(a) * la + to_d(b) * lb + to_d(c) * lc) / s;
    const T cr = cross(b - a, c - a);
    return CircleD(o, std::fabs((Real)cr) / s);
}
template <class T>
CircleD circle_with_diameter(const PointT<T>& a, const PointT<T>& b) { return CircleD(midpoint(a, b), dist(a, b) / 2); }
// a と b（a ≠ b）を通る半径 r の円の中心（0〜2 個）。整数の点と整数の r なら個数は厳密
template <class T, class R>
std::vector<PointD> circles_through(const PointT<T>& a, const PointT<T>& b, R r) {
    const PointT<T> v = b - a;
    if (a == b || r < 0) return {};
    int k;
    Real h2;  // (2h)² = 4r² - |v|²（h は中点から中心までの距離）
    if constexpr (std::is_integral_v<T> && std::is_integral_v<R>) {
        const Int128 d2 = norm2_128(v), f = 4 * (Int128)r * r;
        k = d2 > f ? 0 : (d2 == f ? 1 : 2);
        h2 = (Real)(f - d2);
    } else {
        const Real d = len(v);
        const int s = real_cmp(d, 2 * (Real)r);
        k = s > 0 ? 0 : (s == 0 ? 1 : 2);
        h2 = 4 * (Real)r * (Real)r - d * d;
    }
    if (k == 0) return {};
    const PointD m = midpoint(a, b);
    if (k == 1) return {m};
    const PointD off = unit(rot90(v)) * (std::sqrt(std::max<Real>(0, h2)) / 2);
    return {m + off, m - off};
}
// ---- 円と半直線（円周との交点） ----
template <class T>
int count_cross(const CircleT<T>& c, const RayT<T>& r) {  // 円周と半直線の交点の数（0〜2。始点も含む。整数なら厳密）
    // r.a + t d（t ≥ 0）を代入した A t² + 2B t + C = 0 の解のうち t ≥ 0 のものの数
    const PointT<T> d = r.dir(), w = r.a - c.c;
    int sD, sB, sC;
    if constexpr (std::is_integral_v<T>) {
        const Int128 A = norm2_128(d), B = dot128(d, w), C = norm2_128(w) - (Int128)c.r * c.r, D = B * B - A * C;
        sD = (D > 0) - (D < 0), sB = (B > 0) - (B < 0), sC = (C > 0) - (C < 0);
    } else {
        const Real A = norm2(d), B = dot(d, w), nw = norm2(w), r2 = c.r * c.r, C = nw - r2, cs = std::max(nw, r2);
        sD = real_sign(B * B - A * C, std::max(B * B, A * cs)), sB = real_sign(B, std::sqrt(A * nw)), sC = real_sign(C, cs);
    }
    if (sD < 0) return 0;
    if (sD == 0) return sB <= 0 ? 1 : 0;
    return (sB <= 0 || sC <= 0) + (sB <= 0 && sC >= 0);
}
template <class T>
std::vector<PointD> cross_points(const CircleT<T>& c, const RayT<T>& r) {  // 円周と半直線の交点（始点に近い順）
    const int k = count_cross(c, r);
    if (k == 0) return {};
    std::vector<PointD> v = cross_points(c, LineT<T>(r.a, r.b));
    if (v.size() <= 1 || k == 2) return v;
    return {v[1]};
}

// ---- 円板（円の内側。円周を含む）としての位置関係・距離。整数なら厳密 ----
namespace geo_internal {
template <class T, class X>
int disk_cmp(const CircleT<T>& c, const X& x) {  // (中心と x の距離) と r の比較（-1 / 0 / 1）
    if constexpr (std::is_integral_v<T>) return dist2_cmp(x, c.c, c.r * c.r);
    else return real_cmp(dist(x, c.c), c.r);
}
}  // namespace geo_internal
template <class T>
bool disk_intersects(const CircleT<T>& c, const SegmentT<T>& s) { return geo_internal::disk_cmp(c, s) <= 0; }  // 円板と線分が共有点を持つか
template <class T>
bool disk_intersects(const CircleT<T>& c, const LineT<T>& l) { return geo_internal::disk_cmp(c, l) <= 0; }
template <class T>
bool disk_intersects(const CircleT<T>& c, const RayT<T>& r) { return geo_internal::disk_cmp(c, r) <= 0; }
template <class T>
bool disk_intersects(const CircleT<T>& a, const CircleT<T>& b) {  // 2 つの円板が共有点を持つか（中心の距離 ≤ r1 + r2）
    if constexpr (std::is_integral_v<T>) return norm2_128(b.c - a.c) <= ((Int128)a.r + b.r) * ((Int128)a.r + b.r);
    else return real_cmp(dist(a.c, b.c), a.r + b.r) <= 0;
}
template <class T>
bool disk_contains(const CircleT<T>& c, const SegmentT<T>& s) {  // 線分が円板に入っているか（両端が入っていればよい）
    return contains(c, s.a) != OUTSIDE && contains(c, s.b) != OUTSIDE;
}
template <class T>
bool disk_contains(const CircleT<T>& a, const CircleT<T>& b) {  // 円板 b が円板 a に入っているか（内接・同じ円も含む）
    if constexpr (std::is_integral_v<T>) {
        if (a.r < b.r) return false;
        const Int128 df = (Int128)a.r - b.r;
        return norm2_128(b.c - a.c) <= df * df;
    } else {
        return real_cmp(dist(a.c, b.c) + b.r, a.r) <= 0;
    }
}
// 円板までの距離（円板の中・共有点があれば 0）
template <class T>
Real disk_dist(const CircleT<T>& c, const PointT<T>& p) { return contains(c, p) != OUTSIDE ? 0 : dist(c.c, p) - (Real)c.r; }
template <class T>
Real disk_dist(const CircleT<T>& c, const SegmentT<T>& s) { return disk_intersects(c, s) ? 0 : dist(s, c.c) - (Real)c.r; }
template <class T>
Real disk_dist(const CircleT<T>& c, const LineT<T>& l) { return disk_intersects(c, l) ? 0 : dist(l, c.c) - (Real)c.r; }
template <class T>
Real disk_dist(const CircleT<T>& c, const RayT<T>& r) { return disk_intersects(c, r) ? 0 : dist(r, c.c) - (Real)c.r; }
template <class T>
Real disk_dist(const CircleT<T>& a, const CircleT<T>& b) { return disk_intersects(a, b) ? 0 : dist(a.c, b.c) - (Real)a.r - (Real)b.r; }
// 円周までの距離（中の点でも円周まで）。2 円は円周どうしの距離（共有点があれば 0。片方が中にあれば内側のすき間）
template <class T>
Real circumference_dist(const CircleT<T>& c, const PointT<T>& p) { return std::fabs(dist(c.c, p) - (Real)c.r); }
template <class T>
Real circumference_dist(const CircleT<T>& a, const CircleT<T>& b) {
    const int rel = circle_relation(a, b);
    if (rel == CIRCLE_SEPARATE) return dist(a.c, b.c) - (Real)a.r - (Real)b.r;
    if (rel == CIRCLE_CONTAIN) return std::fabs((Real)a.r - (Real)b.r) - dist(a.c, b.c);
    return 0;
}

// ---- 円周の上の点を結ぶ弦（点は円周に沿った位置で表す：番号・周に沿った座標など、一周の順に並んだ整数） ----
// 弦 ab と cd が円の内部の 1 点で交わるか（端点を共有する・長さ 0 なら false）
inline bool chords_cross(long long a, long long b, long long c, long long d) {
    if (a > b) std::swap(a, b);
    if (c > d) std::swap(c, d);
    return (a < c && c < b && b < d) || (c < a && a < d && d < b);
}
// 円周を n 等分した点 a, b を通る直線の傾きのキー（2 本が平行 ⇔ キーが同じ）。点の番号は 0 〜 n-1 でも 1 〜 n でもよい
inline long long chord_direction_key(long long a, long long b, long long n) { return ((a + b) % n + n) % n; }
// 各 query の弦が、chords のうち何本と交わるか（chords_cross と同じ意味）。オフライン O((N + Q) log N)
inline std::vector<long long> chord_cross_counts(const std::vector<std::pair<long long, long long>>& chords,
                                                 const std::vector<std::pair<long long, long long>>& queries) {
    std::vector<std::pair<long long, long long>> ch;
    for (auto [a, b] : chords) {
        if (a > b) std::swap(a, b);
        if (a < b) ch.emplace_back(a, b);
    }
    std::sort(ch.begin(), ch.end());
    std::vector<long long> bs;
    for (const auto& e : ch) bs.push_back(e.second);
    std::sort(bs.begin(), bs.end());
    const int m = (int)bs.size();
    std::vector<int> bit(m + 1, 0);
    auto add = [&](long long b) {
        for (int i = (int)(std::lower_bound(bs.begin(), bs.end(), b) - bs.begin()) + 1; i <= m; i += i & -i) ++bit[i];
    };
    auto count_less = [&](long long b) {  // 入れた弦のうち他端 < b の数
        long long r = 0;
        for (int i = (int)(std::lower_bound(bs.begin(), bs.end(), b) - bs.begin()); i > 0; i -= i & -i) r += bit[i];
        return r;
    };
    // 交わる ⇔ (c < a < d かつ b > d) または (a < c かつ c < b < d)。「a < X の弦のうち…」の差で数える
    struct Ev {
        long long X, lo, hi;  // a < X の弦のうち lo ≤ b < hi のもの
        int q, sgn;
    };
    std::vector<Ev> ev;
    const long long no_upper = std::numeric_limits<long long>::max();  // 上限なし
    for (int q = 0; q < (int)queries.size(); ++q) {
        auto [c, d] = queries[q];
        if (c > d) std::swap(c, d);
        if (c == d) continue;
        ev.push_back({d, d + 1, no_upper, q, +1}), ev.push_back({c + 1, d + 1, no_upper, q, -1});
        ev.push_back({c, c + 1, d, q, +1});
    }
    std::sort(ev.begin(), ev.end(), [](const Ev& x, const Ev& y) { return x.X < y.X; });
    std::vector<long long> res(queries.size(), 0);
    size_t k = 0;
    for (const Ev& e : ev) {
        while (k < ch.size() && ch[k].first < e.X) add(ch[k++].second);
        const long long hi = e.hi == no_upper ? (long long)k : count_less(e.hi);
        res[e.q] += e.sgn * (hi - count_less(e.lo));
    }
    return res;
}
// 交わる弦の組の数（端点を共有する組は数えない）。O(N log N)
inline long long count_crossing_chords(const std::vector<std::pair<long long, long long>>& chords) {
    long long s = 0;
    for (long long v : chord_cross_counts(chords, chords)) s += v;
    return s / 2;
}
// 弦を足しながら「弦 ab が今ある弦のどれかと交わるか」を調べる（位置は 0 〜 n。1 〜 n の番号もそのまま使える）。O(log n)
struct ChordSet {
    int sz = 1;
    std::vector<long long> mx, mn;  // 位置 a に端がある弦の、もう一方の端の最大（a が小さい側）/ 最小（a が大きい側）
    explicit ChordSet(int n) {
        while (sz < n + 1) sz <<= 1;
        mx.assign(2 * sz, std::numeric_limits<long long>::min()), mn.assign(2 * sz, std::numeric_limits<long long>::max());
    }
    void add(long long a, long long b) {  // 弦 ab を足す
        if (a > b) std::swap(a, b);
        if (a == b) return;
        for (long long i = a + sz; i > 0; i >>= 1) mx[i] = std::max(mx[i], b);
        for (long long i = b + sz; i > 0; i >>= 1) mn[i] = std::min(mn[i], a);
    }
    bool crosses(long long c, long long d) const {  // 弦 cd が今ある弦のどれかと円の内部で交わるか
        if (c > d) std::swap(c, d);
        long long hi = std::numeric_limits<long long>::min(), lo = std::numeric_limits<long long>::max();
        for (long long l = c + 1 + sz, r = d + sz; l < r; l >>= 1, r >>= 1) {  // 位置 (c, d) の範囲
            if (l & 1) hi = std::max(hi, mx[l]), lo = std::min(lo, mn[l]), ++l;
            if (r & 1) --r, hi = std::max(hi, mx[r]), lo = std::min(lo, mn[r]);
        }
        return hi > d || lo < c;
    }
};

inline Real arc_length(Real r, Real theta) { return r * theta; }
inline Real sector_area(Real r, Real theta) { return r * r * theta / 2; }

// 円の中の格子点 (X, Y) の数。中心・半径を scale 倍した整数 cx, cy, r で渡す（(X·scale - cx)² + (Y·scale - cy)² ≤ r²）
// 小数点以下 k 桁の入力なら parse_fixed(s, k) / read_fixed(k) と scale = 10^k。include_border = false なら円周の上を除く。O(r / scale)
inline long long lattice_points_in_circle(Coord cx, Coord cy, Coord r, Coord scale = 1, bool include_border = true) {
    if (r < 0) return 0;
    auto fdiv = [](Coord a, Coord b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); };
    auto cdiv = [&](Coord a, Coord b) { return -fdiv(-a, b); };
    const Int128 r2 = (Int128)r * r;
    long long cnt = 0;
    for (Coord X = cdiv(cx - r, scale), xe = fdiv(cx + r, scale); X <= xe; ++X) {
        const Coord dx = X * scale - cx;
        const Int128 rem = r2 - (Int128)dx * dx;
        if (rem < 0 || (!include_border && rem == 0)) continue;
        const Coord s = isqrt((long long)(include_border ? rem : rem - 1));  // |Y·scale - cy| ≤ s
        const Coord lo = cdiv(cy - s, scale), hi = fdiv(cy + s, scale);
        if (hi >= lo) cnt += hi - lo + 1;
    }
    return cnt;
}
#endif  // GG_CIRCLE
#ifndef GG_CIRCLE_AREA
#define GG_CIRCLE_AREA
// ---- 円の共通部分の面積（幾何第二世代）  2 円・円と多角形・円と直線の左側・弓形 ----
// 弓形の面積 r²/2 (θ - sin θ)。θ が小さいときは桁落ちしないよう級数で
inline Real circular_segment_area(Real r, Real theta) {
    Real f;
    if (theta < 0.01L) {
        const Real t2 = theta * theta;
        f = theta * t2 / 6 * (1 - t2 / 20 * (1 - t2 / 42 * (1 - t2 / 72)));  // θ³/6 - θ⁵/120 + θ⁷/5040 - θ⁹/362880
    } else {
        f = theta - std::sin(theta);
    }
    return r * r / 2 * f;
}

// 2 円の共通部分の面積
template <class T>
Real intersection_area(const CircleT<T>& a, const CircleT<T>& b) {
    const int rel = circle_relation(a, b);
    if (rel == CIRCLE_SEPARATE || rel == CIRCLE_EXTERNAL) return 0;
    const Real r1 = a.r, r2 = b.r;
    if (rel != CIRCLE_INTERSECT) return GEO_PI * std::min(r1, r2) * std::min(r1, r2);  // 片方がもう片方の内側
    // 中心どうしの距離 d、三角形 (d, r1, r2) の面積の 16 倍 K = 4 d² r1² - (d² + r1² - r2²)²
    Real K, c1, c2;  // c1 = d² + r1² - r2²、c2 = d² + r2² - r1²
    if constexpr (std::is_integral_v<T>) {
        const Int128 d2 = norm2_128(b.c - a.c), s1 = (Int128)a.r * a.r, s2 = (Int128)b.r * b.r;
        const Int128 e1 = d2 + s1 - s2;
        K = (Real)(4 * d2 * s1 - e1 * e1), c1 = (Real)e1, c2 = (Real)(d2 + s2 - s1);
    } else {
        const Real d2 = norm2(b.c - a.c);
        c1 = d2 + r1 * r1 - r2 * r2, c2 = d2 + r2 * r2 - r1 * r1;
        K = 4 * d2 * r1 * r1 - c1 * c1;
    }
    const Real sk = std::sqrt(std::max<Real>(0, K));
    const Real a1 = std::atan2(sk, c1), a2 = std::atan2(sk, c2);  // 交点を見込む半分の角
    return circular_segment_area(r1, 2 * a1) + circular_segment_area(r2, 2 * a2);
}

namespace geo_internal {
// 中心を原点に移した円（半径 r）と三角形 (0, a, b) の共通部分の符号付き面積
template <class T>
Real circle_tri_area(const PointT<T>& a, const PointT<T>& b, T r) {
    const PointT<T> d = b - a;
    // |a + t d|² = r² → A t² + 2 B t + C = 0
    Real ts[2];
    int k = 0;
    Real A, B, C, disc;
    if constexpr (std::is_integral_v<T>) {
        const Int128 A1 = norm2_128(d), B1 = dot128(a, d), C1 = norm2_128(a) - (Int128)r * r, D1 = B1 * B1 - A1 * C1;
        A = (Real)A1, B = (Real)B1, C = (Real)C1, disc = D1 > 0 && A1 > 0 ? (Real)D1 : -1;
    } else {
        A = norm2(d), B = dot(a, d), C = norm2(a) - r * r;
        disc = A > 0 ? B * B - A * C : -1;
    }
    // 線分を円周との交点 t1 < t2 で切り、(t1, t2) の部分は円の内側（三角形）、それ以外は外側（扇形）
    // 交点が 2 つないとき（接する・離れている）は全部が外側
    Real t1 = 2, t2 = 2;
    if (disc > 0) {
        const Real sd = std::sqrt(disc);
        const Real q = -(B + (B >= 0 ? sd : -sd));  // 桁落ちを避ける：t1 = q / A、t2 = C / q（t1 t2 = C / A）
        t1 = q / A, t2 = C / q;
        if (t1 > t2) std::swap(t1, t2);
        if (0 < t1 && t1 < 1) ts[k++] = t1;
        if (0 < t2 && t2 < 1) ts[k++] = t2;
    }
    const PointD ad = to_d(a), dd = to_d(d);
    const Real rr = (Real)r * (Real)r;
    Real res = 0, s = 0;
    for (int i = 0; i <= k; ++i) {
        const Real e = i < k ? ts[i] : 1;
        if (e <= s) continue;
        const PointD p = ad + dd * s, q = ad + dd * e;
        const Real mid = (s + e) / 2;
        if (t1 < mid && mid < t2) res += cross(p, q) / 2;                  // 円の内側：三角形
        else res += rr * std::atan2(cross(p, q), dot(p, q)) / 2;            // 外側：扇形
        s = e;
    }
    return res;
}
}  // namespace geo_internal

// 円と多角形（単純多角形。凸でなくてよい・向きはどちらでもよい）の共通部分の面積
template <class T>
Real intersection_area(const CircleT<T>& c, const std::vector<PointT<T>>& poly) {
    const int n = (int)poly.size();
    if (n < 3 || c.r <= 0) return 0;
    Real s = 0;
    for (int i = 0; i < n; ++i) s += geo_internal::circle_tri_area(poly[i] - c.c, poly[(i + 1) % n] - c.c, c.r);
    return std::fabs(s);
}

// 円の、直線 l の左側の部分の面積
template <class T>
Real circle_cut_area(const CircleT<T>& c, const LineT<T>& l) {
    const Real full = GEO_PI * (Real)c.r * (Real)c.r;
    if (c.r <= 0) return 0;
    // 中心から直線までの符号付き距離 s = cross / |dir|（左が正）。|s| < r なら弦がある
    Real cr, h;  // cr = cross(dir, c - a)、h = √(r²|dir|² - cr²)（弦の半分 × |dir|）
    int side;    // |s| と r の比較
    if constexpr (std::is_integral_v<T>) {
        const Int128 x = cross128(l.dir(), c.c - l.a), w = (Int128)c.r * c.r * norm2_128(l.dir()) - x * x;
        cr = (Real)x, h = w > 0 ? std::sqrt((Real)w) : 0, side = w > 0 ? -1 : (w == 0 ? 0 : 1);
    } else {
        const Real x = cross(l.dir(), c.c - l.a), dl = len(l.dir()), w = c.r * c.r * dl * dl - x * x;
        cr = x, h = w > 0 ? std::sqrt(w) : 0, side = real_sign(std::fabs(x) / dl - c.r, c.r) >= 0 ? 1 : -1;
    }
    if (side >= 0) return cr > 0 ? full : 0;  // 直線が円の外（接する）
    const Real seg = circular_segment_area(c.r, 2 * std::atan2(h, std::fabs(cr)));  // 中心から遠い側の弓形
    return cr >= 0 ? full - seg : seg;
}

// ---- 円と多角形の位置関係（多角形は周を含む閉じた領域。円は円板。整数なら厳密） ----
template <class T>
bool contains(const std::vector<PointT<T>>& poly, const CircleT<T>& c) {  // 円が多角形に入っているか
    if (contains(poly, c.c) == OUTSIDE) return false;
    for (int i = 0; i < (int)poly.size(); ++i)
        if (geo_internal::disk_cmp(c, geo_internal::poly_edge(poly, i)) < 0) return false;  // 中心と辺の距離 < r なら辺の外へはみ出す
    return true;
}
template <class T>
bool disk_contains(const CircleT<T>& c, const std::vector<PointT<T>>& poly) {  // 多角形が円板に入っているか
    for (const auto& p : poly)
        if (contains(c, p) == OUTSIDE) return false;
    return true;
}
template <class T>
bool disk_intersects(const CircleT<T>& c, const std::vector<PointT<T>>& poly) {  // 円板と多角形が共有点を持つか
    if (poly.empty()) return false;
    if (contains(poly, c.c) != OUTSIDE) return true;
    for (int i = 0; i < (int)poly.size(); ++i)
        if (geo_internal::disk_cmp(c, geo_internal::poly_edge(poly, i)) <= 0) return true;
    return false;
}
template <class T>
Real disk_dist(const CircleT<T>& c, const std::vector<PointT<T>>& poly) { return disk_intersects(c, poly) ? 0 : dist(poly, c.c) - (Real)c.r; }
#endif  // GG_CIRCLE_AREA
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
namespace ca_test {

#define CA_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[circle_area] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK CA_CHECK

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

// ---- 数値積分（適応シンプソン） ----
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
    // 折れ目でつまずかないように 64 個に分けてから
    Real s = 0;
    const int K = 64;
    for (int k = 0; k < K; ++k) {
        const Real x0 = a + (b - a) * k / K, x1 = a + (b - a) * (k + 1) / K, xm = (x0 + x1) / 2;
        const Real f0 = f(x0), fm = f(xm), f1 = f(x1);
        s += simpson_rec(f, x0, x1, f0, fm, f1, (x1 - x0) / 6 * (f0 + 4 * fm + f1), eps / K, 40);
    }
    return s;
}
// 円の縦の区間 [cy - h, cy + h]（x で）
static bool circle_span(const CircleD& c, Real x, Real& lo, Real& hi) {
    const Real dx = x - c.c.x, w = c.r * c.r - dx * dx;
    if (w <= 0) return false;
    lo = c.c.y - std::sqrt(w), hi = c.c.y + std::sqrt(w);
    return true;
}
// 多角形の x での縦の区間（偶奇の規則）と円の区間の共通部分の長さ
static Real poly_circle_len(const std::vector<PointD>& poly, const CircleD& c, Real x) {
    Real clo, chi;
    if (!circle_span(c, x, clo, chi)) return 0;
    std::vector<Real> ys;
    const int n = (int)poly.size();
    for (int i = 0; i < n; ++i) {
        const PointD p = poly[i], q = poly[(i + 1) % n];
        if ((p.x <= x && x < q.x) || (q.x <= x && x < p.x)) ys.push_back(p.y + (q.y - p.y) * (x - p.x) / (q.x - p.x));
    }
    std::sort(ys.begin(), ys.end());
    Real s = 0;
    for (size_t k = 0; k + 1 < ys.size(); k += 2) s += std::max<Real>(0, std::min(ys[k + 1], chi) - std::max(ys[k], clo));
    return s;
}
// 星形の多角形（中心のまわりに角度の順に点を並べる。単純でないこともあるので is_simple で確かめて使う）
static std::vector<Point> rnd_star(std::mt19937_64& rng, int n, long long C) {
    std::vector<std::pair<Real, Point>> v;
    const Point o = rndp(rng, C / 3);
    for (int k = 0; k < n; ++k) {
        Point p = rndp(rng, C);
        if (p == o) continue;
        v.push_back({angle_of(p - o), p});
    }
    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<Point> r;
    for (size_t k = 0; k < v.size(); ++k)
        if (k == 0 || std::fabs(v[k].first - v[k - 1].first) > 1e-12L) r.push_back(v[k].second);  // 同じ向きの点は 1 つ
    if (rng() % 2) std::reverse(r.begin(), r.end());
    return r;
}

static void test_segment(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const Real r = rndr(rng, 0.1L, 10), th = it % 3 ? rndr(rng, 0, 2 * GEO_PI) : std::pow((Real)10, -rndr(rng, 1, 8));
        const Real a = circular_segment_area(r, th);
        if (th >= 0.01L) CHECK(near(a, r * r / 2 * (th - std::sin(th)), 1e-15L));
        else CHECK(near(a, r * r * th * th * th / 12, th * th));  // 先頭の項
        CHECK(near(circular_segment_area(r, th) + circular_segment_area(r, 2 * GEO_PI - th), GEO_PI * r * r, 1e-15L));
    }
}

static void test_two_circles(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const Circle a(rndp(rng, 4), rnd(rng, 0, 4)), b(rndp(rng, 4), rnd(rng, 0, 4));
        const Real got = intersection_area(a, b);
        CHECK(near(got, intersection_area(b, a), 1e-15L));
        // 数値積分
        const CircleD ad = to_d(a), bd = to_d(b);
        auto f = [&](Real x) {
            Real l1, h1, l2, h2;
            if (!circle_span(ad, x, l1, h1) || !circle_span(bd, x, l2, h2)) return (Real)0;
            return std::max<Real>(0, std::min(h1, h2) - std::max(l1, l2));
        };
        const Real lo = std::max(ad.c.x - ad.r, bd.c.x - bd.r), hi = std::min(ad.c.x + ad.r, bd.c.x + bd.r);
        const Real exp = hi > lo ? integrate(f, lo, hi, 1e-11L) : 0;
        CHECK(std::fabs(got - exp) <= 1e-8L);
        // 実数の円でも同じ・大きく拡大しても比が同じ（10^8 倍 → 面積 10^16 倍）
        CHECK(near(intersection_area(ad, bd), got, 1e-12L));
        const Coord s = 100000000;
        const Circle A(a.c * s, a.r * s), B(b.c * s, b.r * s);
        CHECK(near(intersection_area(A, B) / ((Real)s * s), got, 1e-12L));
    }
    // 接する直前（細い共通部分）：拡大しても比が保たれる
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const Coord R = rnd(rng, 1, 1000000000LL), r = rnd(rng, 1, 1000000000LL);
        const Circle a(Point(0, 0), R), b(Point(R + r - 1, 0), r);  // 外接より 1 だけ近い
        const Real got = intersection_area(a, b);
        // 弓形 2 つ：h = 半弦、角は 2 asin(h / r) を級数で
        CHECK(got > 0 && got < 2 * (Real)std::min(R, r) * 2);
        const Circle c(Point(R - r + 1, 0), r);  // 内接より 1 だけ外
        if (r < R) CHECK(intersection_area(a, c) < GEO_PI * (Real)r * r && intersection_area(a, c) > GEO_PI * (Real)r * r * 0.9L);
    }
}

static void test_polygon(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const std::vector<Point> poly = rnd_star(rng, (int)rnd(rng, 3, 9), 6);
        if (poly.size() < 3 || !is_simple(poly)) continue;
        const Circle c(rndp(rng, 5), rnd(rng, 0, 6));
        const Real got = intersection_area(c, poly);
        std::vector<PointD> pd;
        for (const Point& p : poly) pd.push_back(to_d(p));
        const CircleD cd = to_d(c);
        // 被積分関数が滑らかでない x（頂点・辺と円周の交点）で区切って積分する
        std::vector<Real> xs = {cd.c.x - cd.r, cd.c.x + cd.r};
        for (size_t i = 0; i < pd.size(); ++i) {
            xs.push_back(pd[i].x);
            for (const PointD& q : cross_points(cd, SegmentD(pd[i], pd[(i + 1) % pd.size()]))) xs.push_back(q.x);
        }
        std::sort(xs.begin(), xs.end());
        Real exp = 0;
        for (size_t k = 0; k + 1 < xs.size(); ++k) {
            const Real lo = std::max(xs[k], cd.c.x - cd.r), hi = std::min(xs[k + 1], cd.c.x + cd.r);
            if (c.r > 0 && hi > lo) exp += integrate([&](Real x) { return poly_circle_len(pd, cd, x); }, lo, hi, 1e-12L);
        }
        if (std::fabs(got - exp) > 1e-7L) {
            std::cerr << std::setprecision(17) << "c=" << c << " got=" << got << " exp=" << exp << " poly:";
            for (auto& p : poly) std::cerr << " " << p;
            std::cerr << "\n";
        }
        CHECK(std::fabs(got - exp) <= 1e-7L);
        CHECK(near(intersection_area(cd, pd), got, 1e-12L));
        // 円が大きければ多角形の面積、多角形が大きければ円の面積
        CHECK(near(intersection_area(Circle(c.c, 100), poly), area(poly), 1e-15L));
        const Coord s = 100000000;
        std::vector<Point> big = poly;
        for (Point& p : big) p *= s;
        CHECK(near(intersection_area(Circle(c.c * s, c.r * s), big) / ((Real)s * s), got, 1e-12L));
    }
    CHECK(near(intersection_area(Circle(0, 0, 1), std::vector<Point>{{-5, -5}, {5, -5}, {5, 5}, {-5, 5}}), GEO_PI, 1e-15L));
}

static void test_cut(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const Circle c(rndp(rng, 5), rnd(rng, 0, 5));
        Point a = rndp(rng, 6), b = rndp(rng, 6);
        while (b == a) b = rndp(rng, 6);
        const Line l(a, b);
        const Real got = circle_cut_area(c, l);
        // 左側 = 大きな正方形を l で切った多角形との共通部分
        const Real exp = intersection_area(to_d(c), [&] {
            std::vector<PointD> sq = {{-100, -100}, {100, -100}, {100, 100}, {-100, 100}}, r;
            const LineD ld = to_d(l);
            for (int i = 0; i < 4; ++i) {  // 凸多角形を直線で切る（素朴に）
                const PointD p = sq[i], q = sq[(i + 1) % 4];
                const Real sp = cross(ld.dir(), p - ld.a), sq2 = cross(ld.dir(), q - ld.a);
                if (sp >= 0) r.push_back(p);
                if ((sp > 0 && sq2 < 0) || (sp < 0 && sq2 > 0)) r.push_back(p + (q - p) * (sp / (sp - sq2)));
            }
            return r;
        }());
        CHECK(std::fabs(got - exp) <= 1e-9L);
        CHECK(near(got + circle_cut_area(c, Line(b, a)), GEO_PI * (Real)c.r * c.r, 1e-15L));
        CHECK(near(circle_cut_area(to_d(c), to_d(l)), got, 1e-12L));
    }
}

static void test_relation(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 6000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 4 : 8;
        std::vector<Point> P = rnd_star(rng, (int)rnd(rng, 3, 9), C);
        if (!is_simple(P)) continue;
        const Circle c(rndp(rng, C), rnd(rng, 0, C));
        const Int128 r2 = (Int128)c.r * c.r;
        const int n = (int)P.size();
        // 中心と各辺の距離の 2 乗（有理数）と r² の比較
        int mn = 1;  // 最小の辺までの距離² と r² の比較
        for (int i = 0; i < n; ++i) {
            const Point a = P[i], b = P[(i + 1) % n], d = b - a;
            const Int128 t = dot128(d, c.c - a), L = norm2_128(d);
            Int128 num, den;
            if (t <= 0) num = norm2_128(c.c - a), den = 1;
            else if (t >= L) num = norm2_128(c.c - b), den = 1;
            else num = cross128(d, c.c - a) * cross128(d, c.c - a), den = L;
            const Int128 v = num - r2 * den;
            mn = std::min(mn, (v > 0) - (v < 0));
        }
        const Where w = contains(P, c.c);
        CHECK(contains(P, c) == (w != OUTSIDE && mn >= 0));
        bool all = true;
        for (const Point& p : P) all = all && norm2_128(p - c.c) <= r2;
        CHECK(disk_contains(c, P) == all);
        CHECK(disk_intersects(c, P) == (w != OUTSIDE || mn <= 0));
        // 円周の上の点を細かく取って確かめる（はっきりした場合だけ）
        const std::vector<PointD> PD = [&] { std::vector<PointD> v; for (const Point& p : P) v.push_back(to_d(p)); return v; }();
        bool out = false;
        for (int k = 0; k < 360; ++k) {
            const PointD q = to_d(c.c) + from_polar((Real)c.r, 2 * GEO_PI * k / 360);
            if (contains(PD, q) == OUTSIDE && boundary_dist(PD, q) > 1e-6L) out = true;
        }
        if (out) CHECK(!contains(P, c));
        if (contains(P, c)) CHECK(!out);
        // 距離：共有点がなければ（多角形までの距離）- r
        const Real e = disk_intersects(c, P) ? 0 : dist(PD, to_d(c.c)) - (Real)c.r;
        CHECK(near(disk_dist(c, P), e, TOL) && disk_dist(c, P) >= 0);
        // 実数の型でも同じ（回転・拡大・平行移動）
        const Real th = rndr(rng, 0, 2 * GEO_PI), sc = rndr(rng, 0.1L, 10);
        const PointD sh(rndr(rng, -100, 100), rndr(rng, -100, 100));
        std::vector<PointD> Q;
        for (const Point& p : P) Q.push_back(rotated(p, th) * sc + sh);
        const CircleD cd(rotated(c.c, th) * sc + sh, (Real)c.r * sc);
        CHECK(contains(Q, cd) == contains(P, c) && disk_contains(cd, Q) == disk_contains(c, P) && disk_intersects(cd, Q) == disk_intersects(c, P));
    }
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_relation(rng, cases);
    test_segment(rng, cases);
    test_two_circles(rng, cases);
    test_polygon(rng, cases);
    test_cut(rng, cases);
    std::cout << "[circle_area] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[circle_area] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int q = 1000000;
    auto t0 = std::chrono::steady_clock::now();
    Real s = 0;
    for (int i = 0; i < q; ++i) s += intersection_area(Circle(rndp(rng, 1000000000LL), rnd(rng, 0, 1000000000LL)), Circle(rndp(rng, 1000000000LL), rnd(rng, 0, 1000000000LL)));
    const long long t1 = ms_since(t0);
    const int n = 200000;
    std::vector<Point> poly;
    for (int i = 0; i < n; ++i) {
        const Real t = 2 * GEO_PI * i / n;
        poly.push_back(Point((Coord)(1e9L * std::cos(t)), (Coord)(1e9L * std::sin(t))));
    }
    t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < 10; ++k) s += intersection_area(Circle(rndp(rng, 500000000LL), rnd(rng, 0, 1000000000LL)), poly);
    const long long t2 = ms_since(t0);
    std::cout << "[circle_area] benchmark: intersection_area(円, 円) × " << q << " " << t1 << " ms, intersection_area(円, 多角形 N = " << n << ") × 10 " << t2
              << " ms  (check " << (long long)std::fmod(s, (Real)7) << ")\n";
}

}  // namespace ca_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace ca_demo {

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
#define SHOW(expr, comment) ca_demo::show(#expr, ca_demo::to_s(expr), comment)
#define RUN(...)                   \
    ca_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("円の共通部分の面積 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  circular_segment_area(r, θ)     弓形の面積\n"
              << "  intersection_area(c1, c2)       2 円の共通部分の面積\n"
              << "  intersection_area(c, poly)      円と多角形の共通部分の面積\n"
              << "  circle_cut_area(c, l)           円の、直線 l の左側の部分の面積\n"
              << "  contains(poly, c) / disk_contains(c, poly) / disk_intersects(c, poly) / disk_dist(c, poly)   円と多角形の位置関係\n";

    api("【1】 circular_segment_area(r, θ)", "Real", "中心角 θ（ラジアン）の弓形（弦と弧で囲まれた部分）の面積 r²/2 (θ - sin θ)", "O(1)");
    SHOW(circular_segment_area(1, GEO_PI), "半円 = π/2");
    SHOW(circular_segment_area(2, GEO_PI / 2), "4/2 (π/2 - 1)");

    api("【2】 intersection_area(c1, c2)", "Real", "2 円の共通部分（レンズ形）の面積。整数の円なら場合分けは厳密", "O(1)");
    SHOW(intersection_area(Circle(0, 0, 2), Circle(2, 0, 2)), "半径 2 の円を 2 だけずらす");
    SHOW(intersection_area(Circle(0, 0, 5), Circle(1, 1, 2)), "内側 → 小さい円の面積 4π");
    SHOW(intersection_area(Circle(0, 0, 1), Circle(2, 0, 1)), "外接 → 0");
    note("2 円の和集合の面積 = π r1² + π r2² - intersection_area(c1, c2)");

    api("【3】 intersection_area(c, poly)", "Real", "円と多角形（凸でなくてよい・向きはどちらでもよい）の共通部分の面積", "O(N)");
    RUN(std::vector<Point> sq = {{0, 0}, {2, 0}, {2, 2}, {0, 2}});
    SHOW(intersection_area(Circle(0, 0, 2), sq), "正方形の中の 1/4 円 = π");
    SHOW(intersection_area(Circle(1, 1, 1), sq), "内接円 = π");
    SHOW(intersection_area(Circle(1, 1, 5), sq), "正方形全体 = 4");

    api("【4】 circle_cut_area(c, l)", "Real", "円の、直線 l の左側（l の向きに対して左）の部分の面積", "O(1)");
    SHOW(circle_cut_area(Circle(0, 0, 1), Line(Point(0, 0), Point(1, 0))), "x 軸の上側 = π/2");
    SHOW(circle_cut_area(Circle(0, 0, 2), Line(Point(1, 0), Point(1, 1))), "x = 1 の左側（大きい方）");
    SHOW(circle_cut_area(Circle(0, 0, 2), Line(Point(1, 1), Point(1, 0))), "x = 1 の右側（弓形）");

    api("【5】 contains(poly, c) / disk_contains(c, poly) / disk_intersects(c, poly) / disk_dist(c, poly)", "bool / Real",
        "円（円板）と多角形（周を含む）の位置関係。円が多角形に入る / 多角形が円板に入る / 共有点を持つ / 距離（整数なら判定は厳密）", "O(N)");
    RUN(std::vector<Point> box = {Point(0, 0), Point(10, 0), Point(10, 10), Point(0, 10)});
    SHOW(contains(box, Circle(5, 5, 5)), "内接円は入る（周に接してよい）");
    SHOW(contains(box, Circle(5, 5, 6)), "はみ出す");
    SHOW(disk_contains(Circle(5, 5, 8), box), "4 頂点が半径 8 に入る（5√2 ≈ 7.07）");
    SHOW(disk_contains(Circle(5, 5, 7), box), "頂点がはみ出す");
    SHOW(disk_intersects(Circle(15, 5, 5), box), "辺 x = 10 に接する");
    SHOW(disk_dist(Circle(20, 5, 5), box), "10 - 5 = 5");
    note("ARC051 A（円と長方形：赤だけ・青だけの部分があるか）: 長方形を 4 頂点の多角形にして");
    note("    赤だけがある ⇔ !contains(rect, c)、青だけがある ⇔ !disk_contains(c, rect)（08_rect の to_polygon(r) でも作れる）");
}

}  // namespace ca_demo


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
    ca_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ca_test::benchmark();    // 速度計測
    ca_demo::demo();         // 使い方デモ
    return 0;
}

