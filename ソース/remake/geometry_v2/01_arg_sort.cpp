// =============================================================================
//  偏角（偏角の比較・並べ替え・向きの範囲の数え上げ・角の大きさの比較）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 01 角度・偏角
//  依存       : 00_int_common, 00_geo_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  約束
//    - 偏角は x 軸の正の向きから反時計回りに [0, 2π)。0 ベクトル（原点そのもの）はいちばん小さい
//    - すべて整数で厳密（atan2 を使わない）。ベクトルの成分は 2×10^9 まで（= 座標 10^9 どうしの差）
//    - 実数の点を偏角で並べたいときは 00_geo_base の angle_of(v) をキーにする
//
//  機能一覧                                                                    計算量
//    arg_half(v)                     0: [0, π) / 1: [π, 2π) / -1: 0 ベクトル            O(1)
//    arg_cmp(a, b) / arg_less(a, b)  偏角の比較（-1 / 0 / 1）/ a の偏角が小さいか         O(1)
//    same_direction(a, b)            同じ向きか（どちらも 0 ベクトルでも true）          O(1)
//    arg_sort(ps, o = {0,0}, far_first = false)      o から見た偏角の順に並べた点       O(N log N)
//    arg_sort_idx(ps, o, far_first)                  同じ順の番号。o と同じ点が先頭、同じ向きは近い順（far_first なら遠い順）
//    arg_sort_from(ps, o, start) / arg_sort_idx_from 向き start から反時計回りに一周する順  O(N log N)
//    unique_directions(ps, o, keep_far = false) / unique_directions_idx   向きごとに 1 点（近い点 / 遠い点）。o と同じ点は除く
//    ArgSorted S(ps, o)              偏角の順に並べて、向きの範囲を二分探索で数える      O(N log N)
//      S.size() / S.order[k] / S.vec[k] / S.zeros   方向のある点の数 / k 番目の番号 / そのベクトル / o と同じ点の番号
//      S.lower_bound(dir) / S.upper_bound(dir)       偏角が dir 以上 / dir より大きい最初の位置   O(log N)
//      S.count(l, r, include_l = true, include_r = false)   l から反時計回りに r までの点の数（2π をまたいでよい）  O(log N)
//      S.exists(l, r, ...) / S.groups()              あるか / 同じ向きの点の番号のまとまり
//    angle_cmp(a, b, c, d) / angle_less(a, b, c, d)  ∠(a, b) と ∠(c, d)（なす角 [0, π]）を厳密に比べる   O(1)
//    count_in_open_half(ps, o)       各点 p について、偏角が (θp, θp + π) にある点の数（o から見て p の左側）  O(N log N)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    as_test::random_test()  … 愚直な実装（atan2・全探索）と突き合わせるランダムテスト
//    as_test::benchmark()    … 速度計測
//    as_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
namespace as_test {

#define AS_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[arg_sort] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK AS_CHECK

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

// 愚直：小さい座標なら atan2 で偏角を区別できる（同じ向きかどうかは外積と内積で厳密に）
static Real ang(const Point& v) {
    if (v == Point()) return -1;
    Real t = std::atan2((Real)v.y, (Real)v.x);
    if (t < 0) t += 2 * GEO_PI;
    if (t >= 2 * GEO_PI) t = 0;
    return t;
}
static bool brute_same(const Point& a, const Point& b) {
    if (a == Point() || b == Point()) return a == b;
    return cross128(a, b) == 0 && dot128(a, b) > 0;
}
static int brute_cmp(const Point& a, const Point& b) {
    if (brute_same(a, b)) return 0;
    return ang(a) < ang(b) ? -1 : 1;
}
// 向き l から反時計回りに測った角（l と同じ向きなら 0）
static Real rel_ang(const Point& l, const Point& v) {
    if (brute_same(l, v)) return 0;
    Real d = ang(v) - ang(l);
    if (d < 0) d += 2 * GEO_PI;
    return d;
}
static int brute_count(const std::vector<Point>& ps, const Point& o, const Point& l, const Point& r, bool il, bool ir) {
    int c = 0;
    const bool lr = brute_same(l, r);
    const Real dr = rel_ang(l, r);
    for (const Point& p : ps) {
        const Point v = p - o;
        if (v == Point()) continue;
        const Real d = rel_ang(l, v);
        const bool atl = brute_same(l, v), atr = brute_same(r, v);
        bool in;
        if (lr) in = il && ir && atl;
        else if (atl) in = il;
        else if (atr) in = ir;
        else in = d < dr;
        if (in) ++c;
    }
    return c;
}
static int brute_angle_cmp(const Point& a, const Point& b, const Point& c, const Point& d) {
    auto an = [](const Point& u, const Point& v) -> Real {
        if (u == Point() || v == Point()) return 0;
        return std::atan2(std::fabs((Real)cross128(u, v)), (Real)dot128(u, v));
    };
    const Real x = an(a, b), y = an(c, d);
    if (std::fabs(x - y) < 1e-12L) return 0;
    return x < y ? -1 : 1;
}

static void test_cmp(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 3 == 0 ? 2 : 6;
        Point a = rndp(rng, C), b = rndp(rng, C);
        if (it % 5 == 0) b = a * rnd(rng, 1, 3);
        CHECK(arg_half(a) == (a == Point() ? -1 : (ang(a) < GEO_PI ? 0 : 1)));
        CHECK(arg_cmp(a, b) == brute_cmp(a, b));
        CHECK(arg_less(a, b) == (brute_cmp(a, b) < 0) && same_direction(a, b) == brute_same(a, b));
        CHECK(arg_cmp(b, a) == -arg_cmp(a, b));
    }
    // 大きい成分：同じ向きは整数倍で作る。違う向きは 128 ビットの外積で判定
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const Point a = rndp(rng, 2000000000LL), b = rndp(rng, 2000000000LL);
        if (a == Point() || b == Point()) continue;
        const Real x = ang(a), y = ang(b);
        if (std::fabs(x - y) > 1e-12L) CHECK(arg_cmp(a, b) == (x < y ? -1 : 1));
        const Point d = rndp(rng, 1000);
        if (d == Point()) continue;
        const long long k1 = rnd(rng, 1, 2000000), k2 = rnd(rng, 1, 2000000);
        CHECK(arg_cmp(d * k1, d * k2) == 0 && arg_cmp(d * k1, -d * k2) != 0);
        CHECK(arg_cmp(d * k1, rot90(d) * k2) == (arg_half(d) == arg_half(rot90(d)) ? -1 : (arg_half(d) == 0 ? -1 : 1)));
    }
}

static void test_sort(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 14);
        const long long C = it % 2 ? 3 : 8;
        const std::vector<Point> ps = rndps(rng, n, C);
        const Point o = rndp(rng, C);
        for (int far = 0; far < 2; ++far) {
            std::vector<int> e(n);
            std::iota(e.begin(), e.end(), 0);
            std::sort(e.begin(), e.end(), [&](int i, int j) {
                const Point u = ps[i] - o, v = ps[j] - o;
                const int c = brute_cmp(u, v);
                if (c != 0) return c < 0;
                if (norm2(u) != norm2(v)) return far ? norm2(u) > norm2(v) : norm2(u) < norm2(v);
                return i < j;
            });
            CHECK(arg_sort_idx(ps, o, far) == e);
            const std::vector<Point> sp = arg_sort(ps, o, far);
            for (int k = 0; k < n; ++k) CHECK(sp[k] == ps[e[k]]);
            // 始める向き
            Point st = rndp(rng, C);
            if (st == Point()) st = Point(1, 0);
            std::vector<int> f;
            for (int i : e)
                if (ps[i] == o) f.push_back(i);
            for (int i : e)
                if (ps[i] != o && brute_cmp(ps[i] - o, st) >= 0) f.push_back(i);
            for (int i : e)
                if (ps[i] != o && brute_cmp(ps[i] - o, st) < 0) f.push_back(i);
            CHECK(arg_sort_idx_from(ps, o, st, far) == f);
            const std::vector<Point> fp = arg_sort_from(ps, o, st, far);
            for (int k = 0; k < n; ++k) CHECK(fp[k] == ps[f[k]]);
            // 向きごとに 1 点
            std::vector<int> u;
            for (int i : e) {
                if (ps[i] == o) continue;
                if (!u.empty() && brute_same(ps[u.back()] - o, ps[i] - o)) continue;
                u.push_back(i);
            }
            CHECK(unique_directions_idx(ps, o, far) == u);
            const std::vector<Point> up = unique_directions(ps, o, far);
            CHECK(up.size() == u.size());
            for (size_t k = 0; k < u.size(); ++k) CHECK(up[k] == ps[u[k]]);
        }
    }
    // 大きい座標：並びが偏角の順になっていること
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 300);
        std::vector<Point> ps = rndps(rng, n, 1000000000LL);
        const Point o = rndp(rng, 1000000000LL);
        for (int k = 0; k < n / 3; ++k) ps[rnd(rng, 0, n - 1)] = o + (ps[rnd(rng, 0, n - 1)] - o);  // 重複
        const std::vector<int> idx = arg_sort_idx(ps, o);
        std::vector<int> seen(n, 0);
        for (int i : idx) seen[i]++;
        for (int i = 0; i < n; ++i) CHECK(seen[i] == 1);
        for (int k = 0; k + 1 < n; ++k) {
            const int c = arg_cmp(ps[idx[k]] - o, ps[idx[k + 1]] - o);
            CHECK(c <= 0);
            if (c == 0) CHECK(norm2_128(ps[idx[k]] - o) <= norm2_128(ps[idx[k + 1]] - o));
        }
    }
}

static void test_arg_sorted(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 15);
        const long long C = it % 2 ? 2 : 5;
        const std::vector<Point> ps = rndps(rng, n, C);
        const Point o = rndp(rng, 1);
        const ArgSorted S(ps, o);
        int z = 0;
        for (const Point& p : ps) z += (p == o);
        CHECK((int)S.zeros.size() == z && S.size() == n - z);
        for (int k = 0; k < S.size(); ++k) CHECK(S.vec[k] == ps[S.order[k]] - o);
        for (int q = 0; q < 20; ++q) {
            Point l = rndp(rng, C), r = rndp(rng, C);
            if (l == Point()) l = Point(1, 0);
            if (r == Point()) r = Point(0, -1);
            if (q % 4 == 0) r = l * rnd(rng, 1, 2);
            const bool il = rng() & 1, ir = rng() & 1;
            CHECK(S.count(l, r, il, ir) == brute_count(ps, o, l, r, il, ir));
            CHECK(S.exists(l, r, il, ir) == (brute_count(ps, o, l, r, il, ir) > 0));
            int lb = 0, ub = 0;
            for (int k = 0; k < S.size(); ++k) lb += brute_cmp(S.vec[k], l) < 0, ub += brute_cmp(S.vec[k], l) <= 0;
            CHECK(S.lower_bound(l) == lb && S.upper_bound(l) == ub);
        }
        const auto g = S.groups();
        int tot = 0;
        for (size_t a = 0; a < g.size(); ++a) {
            tot += (int)g[a].size();
            for (size_t b = 0; b < g[a].size(); ++b) CHECK(brute_same(ps[g[a][b]] - o, ps[g[a][0]] - o));
            if (a + 1 < g.size()) CHECK(brute_cmp(ps[g[a][0]] - o, ps[g[a + 1][0]] - o) < 0);
        }
        CHECK(tot == S.size());
        // 半平面の点の数
        const std::vector<int> h = count_in_open_half(ps, o);
        for (int i = 0; i < n; ++i) {
            int c = 0;
            if (ps[i] != o)
                for (int j = 0; j < n; ++j) c += cross128(ps[i] - o, ps[j] - o) > 0;
            CHECK(h[i] == c);
        }
    }
    // 大きい座標の半平面
    for (int it = 0; it < 200 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 60);
        std::vector<Point> ps = rndps(rng, n, 1000000000LL);
        const Point o = rndp(rng, 1000000000LL);
        if (it % 2) {  // o を通る直線の上に点を並べる
            const Point d = rndp(rng, 5);
            for (int k = 0; k < n / 2 && d != Point(); ++k) ps[k] = o + d * rnd(rng, -100000, 100000);
        }
        const std::vector<int> h = count_in_open_half(ps, o);
        for (int i = 0; i < n; ++i) {
            int c = 0;
            if (ps[i] != o)
                for (int j = 0; j < n; ++j) c += cross128(ps[i] - o, ps[j] - o) > 0;
            CHECK(h[i] == c);
        }
    }
}

static void test_angle_cmp(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 3 : 6;
        const Point a = rndp(rng, C), b = rndp(rng, C), c = rndp(rng, C), d = rndp(rng, C);
        CHECK(angle_cmp(a, b, c, d) == brute_angle_cmp(a, b, c, d));
        CHECK(angle_less(a, b, c, d) == (brute_angle_cmp(a, b, c, d) < 0));
    }
    // 大きい座標：等しい角を回転・拡大・入れ替えで作る
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const Point a = rndp(rng, 1000000000LL), b = rndp(rng, 1000000000LL);
        if (a == Point() || b == Point()) continue;
        const long long k = rnd(rng, 1, 2);
        CHECK(angle_cmp(a, b, rot90(a) * k, rot90(b)) == 0 && angle_cmp(a, b, b, a) == 0 && angle_cmp(a, b, -a, -b) == 0);
        CHECK(angle_cmp(a, b, a, -b) == brute_angle_cmp(a, b, a, -b) || std::fabs(angle_between(a, b) - GEO_PI / 2) < 1e-12L);
        const Point c = rndp(rng, 2000000000LL), d = rndp(rng, 2000000000LL);
        const Real x = angle_between(a, b), y = angle_between(c, d);
        if (c != Point() && d != Point() && std::fabs(x - y) > 1e-12L) CHECK(angle_cmp(a, b, c, d) == (x < y ? -1 : 1));
        CHECK(angle_cmp(a, a * 2, Point(1, 0), Point(5, 0)) == 0 && angle_cmp(a, -a, Point(0, 1), Point(0, -3)) == 0);
        CHECK(angle_cmp(a, -a, a, b) >= 0);
    }
    CHECK(angle_cmp(Point(0, 0), Point(1, 0), Point(1, 0), Point(1, 0)) == 0);
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_cmp(rng, cases);
    test_sort(rng, cases);
    test_arg_sorted(rng, cases);
    test_angle_cmp(rng, cases);
    std::cout << "[arg_sort] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[arg_sort] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 200000;
    const std::vector<Point> ps = rndps(rng, n, 1000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    const std::vector<int> idx = arg_sort_idx(ps);
    const long long t_sort = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const ArgSorted S(ps);
    long long s = 0;
    for (int i = 0; i < n; ++i) s += S.count(ps[i], ps[(i + 1) % n]);
    const long long t_cnt = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<int> h = count_in_open_half(ps);
    const long long t_half = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for (int i = 0; i + 3 < n; ++i) s2 += angle_cmp(ps[i], ps[i + 1], ps[i + 2], ps[i + 3]);
    const long long t_ang = ms_since(t0);
    std::cout << "[arg_sort] benchmark: N = " << n << "  arg_sort_idx " << t_sort << " ms, ArgSorted + count × N " << t_cnt
              << " ms, count_in_open_half " << t_half << " ms, angle_cmp × N " << t_ang << " ms  (check " << idx[0] + s % 7 + h[0] % 7 + s2 % 7 << ")\n";
}

}  // namespace as_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace as_demo {

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
#define SHOW(expr, comment) as_demo::show(#expr, as_demo::to_s(expr), comment)
#define RUN(...)                   \
    as_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("偏角 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  arg_cmp(a, b) / arg_less(a, b) / same_direction(a, b) / arg_half(v)   偏角の比較（整数で厳密）\n"
              << "  arg_sort(ps, o, far_first) / arg_sort_idx(...)        o から見た偏角の順（同じ向きは近い順）\n"
              << "  arg_sort_from(ps, o, start) / arg_sort_idx_from(...)  向き start から一周する順\n"
              << "  unique_directions(ps, o, keep_far) / _idx             向きごとに 1 点\n"
              << "  ArgSorted S(ps, o);  S.count(l, r, 端を含むか...) / S.lower_bound(dir) / S.groups()   向きの範囲の点の数\n"
              << "  angle_cmp(a, b, c, d) / angle_less                    ∠(a, b) と ∠(c, d) を厳密に比べる\n"
              << "  count_in_open_half(ps, o)                             各点の左側（(θ, θ + π)）にある点の数\n";

    api("【1】 arg_cmp(a, b) / arg_less(a, b) / same_direction(a, b)", "int / bool",
        "x 軸の正の向きから反時計回りに測った偏角 [0, 2π) を、atan2 を使わず整数で比べる", "O(1)");
    SHOW(arg_cmp(Point(1, 1), Point(-1, 1)), "45 度 < 135 度");
    SHOW(arg_cmp(Point(0, -1), Point(1, 0)), "270 度 > 0 度");
    SHOW(same_direction(Point(2, 4), Point(1, 2)), "同じ向き（長さは関係ない）");
    SHOW(same_direction(Point(2, 4), Point(-1, -2)), "逆向きは別");
    note("sort(v.begin(), v.end(), arg_less) で偏角の順に並ぶ（同じ向きの順は決まらない。決めたいなら arg_sort）");

    api("【2】 arg_sort(ps, o, far_first) / arg_sort_idx(ps, o, far_first)", "vector<Point> / vector<int>",
        "o から見た偏角の順に並べる。o と同じ点が先頭、同じ向きは o に近い順（far_first = true なら遠い順）", "O(N log N)");
    RUN(std::vector<Point> ps = {{2, 0}, {0, 1}, {-1, -1}, {1, 0}, {0, 0}, {-3, 2}});
    SHOW(arg_sort(ps), "原点 (0, 0) が先頭、(1, 0) と (2, 0) は近い順");
    SHOW(arg_sort_idx(ps), "番号で");
    SHOW(arg_sort(ps, Point(0, 0), true), "同じ向きは遠い順");

    api("【3】 arg_sort_from(ps, o, start) / unique_directions(ps, o, keep_far)", "vector<Point>",
        "向き start から反時計回りに一周する順 / 向きごとに 1 点（o と同じ点は除く）", "O(N log N)");
    SHOW(arg_sort_from(ps, Point(0, 0), Point(-1, 0)), "180 度の向きから始める");
    SHOW(unique_directions(ps), "(1, 0) と (2, 0) は近い方だけ");
    SHOW(unique_directions(ps, Point(0, 0), true), "遠い方を残す");

    api("【4】 ArgSorted S(ps, o);  S.count(l, r, include_l = true, include_r = false)", "int",
        "向き l から反時計回りに向き r までにある点の数（2π をまたいでもよい）。扇形の中の点の数", "前処理 O(N log N)、1 回 O(log N)");
    RUN(ArgSorted S(ps));
    SHOW(S.size(), "方向のある点の数（原点と同じ点は S.zeros へ）");
    SHOW(S.count(Point(1, 0), Point(0, 1)), "[0 度, 90 度)：(1,0) と (2,0)");
    SHOW(S.count(Point(1, 0), Point(0, 1), true, true), "[0 度, 90 度]：(0,1) も入る");
    SHOW(S.count(Point(0, -1), Point(1, 1)), "270 度 → 45 度（0 度をまたぐ）");
    SHOW(S.count(Point(1, 0), Point(1, 0), true, true), "l と r が同じ向き：その向きの点だけ（一周なら S.size()）");
    SHOW(S.groups().size(), "向きの種類の数");

    api("【5】 angle_cmp(a, b, c, d) / angle_less(a, b, c, d)", "int / bool",
        "∠(a, b)（a と b のなす角 [0, π]）と ∠(c, d) を誤差なしで比べる。∠AOB なら (A - O, B - O)", "O(1)");
    SHOW(angle_cmp(Point(1, 0), Point(1, 1), Point(0, 1), Point(-1, 0)), "45 度 < 90 度");
    SHOW(angle_cmp(Point(3, 4), Point(-4, 3), Point(1, 0), Point(0, -7)), "どちらも 90 度");
    SHOW(angle_less(Point(1, 0), Point(-1, 1), Point(1, 0), Point(-1, -1)), "135 度と 135 度（下向きでもなす角は同じ）");

    api("【6】 count_in_open_half(ps, o)", "vector<int>",
        "各点 p について、o から見て p の左側（偏角が (θp, θp + π)）にある点の数", "O(N log N)");
    {
        const std::vector<Point> qs = {{1, 0}, {0, 1}, {-1, -1}, {2, -1}, {-1, 2}};
        SHOW(count_in_open_half(qs), "(1,0) の左側は (0,1) と (-1,2)");
        note("原点を内側に含む三角形の数 = C(N, 3) - Σ C(h_i, 2)（原点を通る直線の上に 2 点が並ばないとき）");
        long long n = (long long)qs.size(), t = n * (n - 1) * (n - 2) / 6;
        for (int h : count_in_open_half(qs)) t -= (long long)h * (h - 1) / 2;
        std::cout << "    この例では C(5, 3) = 10 個のうち " << t << " 個\n";
    }

    head("【7】 例: ABC139 F（ベクトルをいくつか選んで和の長さを最大に）");
    note("偏角の順に並べると、選ぶべきベクトルは連続した区間（一周してよい）になる。全区間を試す");
    {
        const std::vector<Point> vs = {{0, 10}, {5, -5}, {-5, -5}, {3, 3}};
        const std::vector<Point> s = arg_sort(vs);
        const int n = (int)s.size();
        Coord best = 0;
        for (int l = 0; l < n; ++l) {
            Point sum;
            for (int len = 1; len <= n; ++len) sum += s[(l + len - 1) % n], best = std::max(best, norm2(sum));
        }
        code("auto s = arg_sort(vs);  for (l) { Point sum; for (len = 1..N) { sum += s[(l + len - 1) % N]; best = max(best, norm2(sum)); } }");
        std::cout << "    (0,10) (5,-5) (-5,-5) (3,3) → 最大の長さ = " << to_s(std::sqrt((Real)best)) << "（(0,10) + (3,3) = (3,13)）\n";
    }

    head("【8】 例: 角度の区間を重ならないように選ぶ（ABC225 E の形）");
    note("区間を [l_i, r_i]（向きのベクトル、反時計回りに l → r）で表し、r の偏角の順に並べて貪欲に選ぶ");
    code("sort(id, [&](int i, int j) { return arg_less(R[i], R[j]); });");
    code("Point last = 最初は選んだものなし;  for (int i : id) if (なし || arg_cmp(last, L[i]) <= 0) { ++ans; last = R[i]; }");
    note("端が接してよいなら arg_cmp(last, L[i]) <= 0、接してはいけないなら < 0。区間が 0 度をまたがない形に直してから使う");
    note("類題: 典型90 009（最大の角）は 06_point_set の max_angle、ABC033 D（三角形の分類）は count_triangles_by_angle");
}

}  // namespace as_demo


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
    as_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    as_test::benchmark();    // 速度計測
    as_demo::demo();         // 使い方デモ
    return 0;
}

