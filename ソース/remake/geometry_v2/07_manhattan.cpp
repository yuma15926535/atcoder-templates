// =============================================================================
//  マンハッタン距離・チェビシェフ距離（45 度回転・最も遠い点（静的・追加削除あり）・距離の和・中央値）   --- 幾何第二世代
// -----------------------------------------------------------------------------
//  分類       : 07 マンハッタン距離・チェビシェフ距離
//  依存       : 00_int_common, 00_geo_base（同梱済み。このファイルだけ貼れば動く）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。__int128 を使う（GCC / Clang）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  考え方
//    u = x + y, v = x - y と 45 度回すと、マンハッタン距離 |dx| + |dy| = max(|du|, |dv|)（チェビシェフ距離）になる
//    チェビシェフ距離 max(|dx|, |dy|) = (|du| + |dv|) / 2
//
//  機能一覧                                                                    計算量
//    manhattan(a, b) / chebyshev(a, b)               マンハッタン距離 / チェビシェフ距離           O(1)
//    manhattan_to_chebyshev(p) / chebyshev_to_manhattan(q)   (x + y, x - y) / 元に戻す（q.x と q.y の偶奇が同じとき）
//        ※ グラフ・グリッドの 10_grid_transform の (x, y) 版と同じ名前（Point 版の重ね定義。両方貼ってよい）
//    in_diamond(c, D, p)                             マンハッタン距離 ≤ D か（ひし形の内側）         O(1)
//    ManhattanFarthest F(ps)                         点の集合（変わらない）                          O(N)
//      F.max_pair()                                  最も遠い 2 点 ManhattanPair{i, j, d}（i < j）   O(1)
//      F.farthest(q)                                 q から最も遠い点 ManhattanFar{i, d}             O(1)
//    DynamicManhattan M                              点の追加・削除がある場合
//      M.add(p) → 番号 / M.erase(番号)（消えていたら false）/ M.size() / M.alive(番号) / M.point(番号)   O(log N)
//      M.max_pair() / M.farthest(q)                  同上（今ある点だけ）                            O(1)
//    sum_manhattan_pairs(ps) / sum_chebyshev_pairs(ps)   全部の組の距離の和（Int128）                 O(N log N)
//    sum_manhattan_from_each(ps)                     各点から全点へのマンハッタン距離の和            O(N log N)
//    manhattan_median(ps)                            距離の和が最小になる格子点（x, y それぞれの中央値）  O(N)
//    AbsSum S(xs [, ws])                             1 次元の f(x) = Σ w_i |x - x_i|。S(x)（Int128）は O(log N)、S.argmin() は重み付きの中央値
//    abs_sum_table(xs, lo, hi)                       x = lo..hi の f(x) = Σ |x - x_i| の表（vector<long long>）  O(N log N + (hi - lo))
//                                                    マンハッタン距離の和 ≤ D の格子点を数える（ABC366 E）：x と y の表を作って尺取り
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mh_test::random_test()  … 全部の組を調べる愚直な実装と突き合わせるランダムテスト
//    mh_test::benchmark()    … 速度計測
//    mh_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#ifndef GG_MANHATTAN
#define GG_MANHATTAN
// ---- マンハッタン距離・チェビシェフ距離（幾何第二世代）  u = x + y, v = x - y で 45 度回す ----
struct ManhattanPair {  // 2 点の番号（i < j。点が 1 つだけなら i == j）と距離。点がなければ -1
    int i = -1, j = -1;
    Coord d = -1;
};
struct ManhattanFar {  // いちばん遠い点の番号と距離。点がなければ -1
    int i = -1;
    Coord d = -1;
};
inline Coord manhattan(const Point& a, const Point& b) {
    const Point d = a - b;
    return (d.x < 0 ? -d.x : d.x) + (d.y < 0 ? -d.y : d.y);
}
inline Coord chebyshev(const Point& a, const Point& b) {
    const Point d = a - b;
    return std::max(d.x < 0 ? -d.x : d.x, d.y < 0 ? -d.y : d.y);
}
inline Point manhattan_to_chebyshev(const Point& p) { return Point(p.x + p.y, p.x - p.y); }
inline Point chebyshev_to_manhattan(const Point& q) { return Point((q.x + q.y) / 2, (q.x - q.y) / 2); }  // q.x と q.y の偶奇が同じとき
inline bool in_diamond(const Point& c, Coord D, const Point& p) { return manhattan(c, p) <= D; }

namespace geo_internal {
// u・v の最小・最大（番号つき）から、最も遠い 2 点・q から最も遠い点を求める
inline ManhattanPair mh_max_pair(int n, Coord umin, int iumin, Coord umax, int iumax, Coord vmin, int ivmin, Coord vmax, int ivmax, int other) {
    if (n == 0) return ManhattanPair();
    const Coord du = umax - umin, dv = vmax - vmin;
    int a = du >= dv ? iumin : ivmin, b = du >= dv ? iumax : ivmax;
    if (a == b) {  // 全部が同じ点
        if (n == 1) return ManhattanPair{a, a, 0};
        b = other;
    }
    return ManhattanPair{std::min(a, b), std::max(a, b), std::max(du, dv)};
}
inline ManhattanFar mh_farthest(int n, const Point& q, Coord umin, int iumin, Coord umax, int iumax, Coord vmin, int ivmin, Coord vmax, int ivmax) {
    if (n == 0) return ManhattanFar();
    const Coord qu = q.x + q.y, qv = q.x - q.y;
    ManhattanFar r{iumin, qu - umin};
    if (umax - qu > r.d) r = ManhattanFar{iumax, umax - qu};
    if (qv - vmin > r.d) r = ManhattanFar{ivmin, qv - vmin};
    if (vmax - qv > r.d) r = ManhattanFar{ivmax, vmax - qv};
    return r;
}
}  // namespace geo_internal

// 変わらない点の集合で、最も遠い 2 点・ある点から最も遠い点（典型90 036）
struct ManhattanFarthest {
    int n = 0;
    Coord umin = 0, umax = 0, vmin = 0, vmax = 0;
    int iumin = -1, iumax = -1, ivmin = -1, ivmax = -1, second = -1;  // second：0 番と違う番号（全部が同じ点のとき用）
    ManhattanFarthest() {}
    explicit ManhattanFarthest(const std::vector<Point>& ps) : n((int)ps.size()) {
        for (int i = 0; i < n; ++i) {
            const Coord u = ps[i].x + ps[i].y, v = ps[i].x - ps[i].y;
            if (i == 0 || u < umin) umin = u, iumin = i;
            if (i == 0 || u > umax) umax = u, iumax = i;
            if (i == 0 || v < vmin) vmin = v, ivmin = i;
            if (i == 0 || v > vmax) vmax = v, ivmax = i;
        }
        second = n >= 2 ? 1 : -1;
    }
    ManhattanPair max_pair() const {  // 最も遠い 2 点（i < j）と距離
        return geo_internal::mh_max_pair(n, umin, iumin, umax, iumax, vmin, ivmin, vmax, ivmax, iumin == 0 ? second : 0);
    }
    ManhattanFar farthest(const Point& q) const { return geo_internal::mh_farthest(n, q, umin, iumin, umax, iumax, vmin, ivmin, vmax, ivmax); }
};

// 点の追加・削除があるとき（u, v を set で持つ）
struct DynamicManhattan {
    std::vector<Point> pts;
    std::vector<char> alive_;
    std::set<std::pair<Coord, int>> su, sv;
    int add(const Point& p) {  // 点を足して番号（0, 1, 2, ... の順）を返す
        const int id = (int)pts.size();
        pts.push_back(p), alive_.push_back(1);
        su.insert({p.x + p.y, id}), sv.insert({p.x - p.y, id});
        return id;
    }
    bool erase(int id) {  // 消す（もう消えている・ない番号なら何もせず false）
        if (id < 0 || id >= (int)pts.size() || !alive_[id]) return false;
        alive_[id] = 0;
        su.erase({pts[id].x + pts[id].y, id}), sv.erase({pts[id].x - pts[id].y, id});
        return true;
    }
    int size() const { return (int)su.size(); }
    bool alive(int id) const { return 0 <= id && id < (int)pts.size() && alive_[id]; }
    Point point(int id) const { return pts[id]; }
    ManhattanPair max_pair() const {
        if (su.empty()) return ManhattanPair();
        // set は (値, 番号) の順なので、全部が同じ点でも最小と最大の番号は違う（1 点だけのときは同じ）
        return geo_internal::mh_max_pair(size(), su.begin()->first, su.begin()->second, su.rbegin()->first, su.rbegin()->second, sv.begin()->first,
                                         sv.begin()->second, sv.rbegin()->first, sv.rbegin()->second, -1);
    }
    ManhattanFar farthest(const Point& q) const {
        if (su.empty()) return ManhattanFar();
        return geo_internal::mh_farthest(size(), q, su.begin()->first, su.begin()->second, su.rbegin()->first, su.rbegin()->second, sv.begin()->first,
                                         sv.begin()->second, sv.rbegin()->first, sv.rbegin()->second);
    }
};

namespace geo_internal {
inline Int128 sum_abs_pairs(std::vector<Coord> a) {  // Σ_{i<j} |a_i - a_j|
    std::sort(a.begin(), a.end());
    Int128 s = 0;
    const long long n = (long long)a.size();
    for (long long k = 0; k < n; ++k) s += (Int128)a[k] * (2 * k - n + 1);
    return s;
}
}  // namespace geo_internal
inline Int128 sum_manhattan_pairs(const std::vector<Point>& ps) {  // 全部の組のマンハッタン距離の和
    std::vector<Coord> xs, ys;
    for (const Point& p : ps) xs.push_back(p.x), ys.push_back(p.y);
    return geo_internal::sum_abs_pairs(xs) + geo_internal::sum_abs_pairs(ys);
}
inline Int128 sum_chebyshev_pairs(const std::vector<Point>& ps) {  // 全部の組のチェビシェフ距離の和（= (Σ|du| + Σ|dv|) / 2）
    std::vector<Coord> us, vs;
    for (const Point& p : ps) us.push_back(p.x + p.y), vs.push_back(p.x - p.y);
    return (geo_internal::sum_abs_pairs(us) + geo_internal::sum_abs_pairs(vs)) / 2;
}
// 各点から全点へのマンハッタン距離の和（N · 4×10^9 程度まで long long に収まる）
inline std::vector<long long> sum_manhattan_from_each(const std::vector<Point>& ps) {
    const int n = (int)ps.size();
    std::vector<long long> res(n, 0);
    for (int axis = 0; axis < 2; ++axis) {
        std::vector<int> id(n);
        for (int i = 0; i < n; ++i) id[i] = i;
        auto key = [&](int i) { return axis == 0 ? ps[i].x : ps[i].y; };
        std::sort(id.begin(), id.end(), [&](int a, int b) { return key(a) < key(b); });
        std::vector<long long> pre(n + 1, 0);
        for (int k = 0; k < n; ++k) pre[k + 1] = pre[k] + key(id[k]);
        for (int k = 0; k < n; ++k) {
            const long long x = key(id[k]);
            res[id[k]] += x * k - pre[k] + (pre[n] - pre[k + 1]) - x * (n - 1 - k);
        }
    }
    return res;
}
// マンハッタン距離の和が最小になる格子点（x, y それぞれの下側の中央値）。点がなければ (0, 0)
inline Point manhattan_median(const std::vector<Point>& ps) {
    if (ps.empty()) return Point();
    std::vector<Coord> xs, ys;
    for (const Point& p : ps) xs.push_back(p.x), ys.push_back(p.y);
    const size_t m = (ps.size() - 1) / 2;
    std::nth_element(xs.begin(), xs.begin() + (long)m, xs.end());
    std::nth_element(ys.begin(), ys.begin() + (long)m, ys.end());
    return Point(xs[m], ys[m]);
}

// ---- 1 次元の距離の和 f(x) = Σ w_i |x - x_i|（マンハッタン距離の和は x と y に分けて足す） ----
struct AbsSum {
    std::vector<Coord> xs;        // x_i を小さい順に
    std::vector<Int128> pw, pwx;  // 重みの累積和・重み × 座標の累積和
    AbsSum() : pw(1, 0), pwx(1, 0) {}
    explicit AbsSum(const std::vector<Coord>& x, const std::vector<Coord>& w = {}) {  // w を省くと重みは全部 1（重みは 0 以上）
        const int n = (int)x.size();
        std::vector<int> id(n);
        for (int i = 0; i < n; ++i) id[i] = i;
        std::sort(id.begin(), id.end(), [&](int a, int b) { return x[a] < x[b]; });
        xs.resize(n), pw.assign(n + 1, 0), pwx.assign(n + 1, 0);
        for (int k = 0; k < n; ++k) {
            const Coord wk = w.empty() ? 1 : w[id[k]];
            xs[k] = x[id[k]], pw[k + 1] = pw[k] + wk, pwx[k + 1] = pwx[k] + (Int128)wk * xs[k];
        }
    }
    Int128 operator()(Coord x) const {  // f(x)。O(log N)
        const int n = (int)xs.size(), k = (int)(std::upper_bound(xs.begin(), xs.end(), x) - xs.begin());
        return (Int128)x * pw[k] - pwx[k] + (pwx[n] - pwx[k]) - (Int128)x * (pw[n] - pw[k]);
    }
    Coord argmin() const {  // f を最小にする x のうち最小のもの（重み付きの中央値）。点がなければ 0
        const int n = (int)xs.size();
        for (int k = 0; k < n; ++k)
            if (2 * pw[k + 1] >= pw[n]) return xs[k];
        return 0;
    }
};
// x = lo, lo + 1, ..., hi での f(x) = Σ |x - x_i|（重みなし）の表。O(N log N + (hi - lo))。値は long long に収まる範囲で
inline std::vector<long long> abs_sum_table(std::vector<Coord> xs, Coord lo, Coord hi) {
    if (hi < lo) return {};
    std::sort(xs.begin(), xs.end());
    const long long n = (long long)xs.size();
    std::vector<long long> res((size_t)(hi - lo + 1));
    long long f = 0;
    for (const Coord v : xs) f += v > lo ? v - lo : lo - v;
    size_t k = std::upper_bound(xs.begin(), xs.end(), lo) - xs.begin();  // x_i ≤ x の数
    for (Coord x = lo;; ++x) {
        res[(size_t)(x - lo)] = f;
        if (x == hi) break;
        f += 2 * (long long)k - n;  // x → x + 1 で、x_i ≤ x の点は 1 遠く、ほかは 1 近くなる
        while (k < xs.size() && xs[k] <= x + 1) ++k;
    }
    return res;
}
#endif  // GG_MANHATTAN
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
namespace mh_test {

#define MH_CHECK(...)                                                                                 \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[manhattan] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n";    \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)
#define CHECK MH_CHECK

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

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long C = it % 2 ? 5 : 1000000000LL;
        const Point a = rndp(rng, C), b = rndp(rng, C), c = rndp(rng, C);
        CHECK(manhattan(a, b) == std::llabs(a.x - b.x) + std::llabs(a.y - b.y) && chebyshev(a, b) == std::max(std::llabs(a.x - b.x), std::llabs(a.y - b.y)));
        const Point ua = manhattan_to_chebyshev(a), ub = manhattan_to_chebyshev(b);
        CHECK(chebyshev(ua, ub) == manhattan(a, b) && chebyshev_to_manhattan(ua) == a);
        CHECK(manhattan(ua, ub) == 2 * chebyshev(a, b));
        const Coord D = rnd(rng, 0, 2 * C);
        CHECK(in_diamond(c, D, a) == (manhattan(c, a) <= D));
    }
}

static void test_farthest(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 30);
        const long long C = it % 3 == 0 ? 1000000000LL : (it % 3 == 1 ? 2 : 50);
        std::vector<Point> ps = rndps(rng, n, C);
        if (it % 7 == 0)
            for (auto& p : ps) p = Point(3, 3);  // 全部が同じ点
        const ManhattanFarthest F(ps);
        const ManhattanPair mp = F.max_pair();
        Coord best = -1;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) best = std::max(best, manhattan(ps[i], ps[j]));
        if (n == 0) CHECK(mp.i == -1 && mp.d == -1);
        else if (n == 1) CHECK(mp.i == 0 && mp.j == 0 && mp.d == 0);
        else CHECK(mp.d == best && mp.i < mp.j && manhattan(ps[mp.i], ps[mp.j]) == best);
        for (int q = 0; q < 10; ++q) {
            const Point p = rndp(rng, C);
            const ManhattanFar f = F.farthest(p);
            Coord e = -1;
            for (int i = 0; i < n; ++i) e = std::max(e, manhattan(p, ps[i]));
            CHECK(f.d == e && (n == 0 ? f.i == -1 : manhattan(p, ps[f.i]) == e));
        }
    }
    // 追加・削除
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        DynamicManhattan M;
        std::vector<Point> all;
        std::vector<int> alive;
        const long long C = it % 2 ? 3 : 1000000000LL;
        for (int op = 0; op < 200; ++op) {
            const int t = (int)rnd(rng, 0, 3);
            if (t <= 1) {
                const Point p = it % 5 == 0 ? Point(1, 1) : rndp(rng, C);
                CHECK(M.add(p) == (int)all.size());
                all.push_back(p), alive.push_back(1);
            } else if (t == 2) {
                const int id = (int)rnd(rng, -1, (int)all.size());
                const bool e = 0 <= id && id < (int)all.size() && alive[id];
                CHECK(M.erase(id) == e);
                if (e) alive[id] = 0;
            }
            int cnt = 0;
            Coord best = -1;
            for (size_t i = 0; i < all.size(); ++i) {
                if (!alive[i]) continue;
                ++cnt;
                CHECK(M.alive((int)i) && M.point((int)i) == all[i]);
                for (size_t j = i + 1; j < all.size(); ++j)
                    if (alive[j]) best = std::max(best, manhattan(all[i], all[j]));
            }
            CHECK(M.size() == cnt);
            const ManhattanPair mp = M.max_pair();
            if (cnt == 0) CHECK(mp.i == -1);
            else if (cnt == 1) CHECK(mp.i == mp.j && alive[mp.i] && mp.d == 0);
            else CHECK(mp.d == best && mp.i < mp.j && alive[mp.i] && alive[mp.j] && manhattan(all[mp.i], all[mp.j]) == best);
            const Point q = rndp(rng, C);
            Coord e = -1;
            for (size_t i = 0; i < all.size(); ++i)
                if (alive[i]) e = std::max(e, manhattan(q, all[i]));
            const ManhattanFar f = M.farthest(q);
            CHECK(f.d == e && (cnt == 0 ? f.i == -1 : (alive[f.i] && manhattan(q, all[f.i]) == e)));
        }
    }
}

static void test_sums(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 40);
        const long long C = it % 2 ? 1000000000LL : 6;
        const std::vector<Point> ps = rndps(rng, n, C);
        Int128 sm = 0, sc = 0;
        std::vector<long long> each(n, 0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (i < j) sm += manhattan(ps[i], ps[j]), sc += chebyshev(ps[i], ps[j]);
                each[i] += manhattan(ps[i], ps[j]);
            }
        CHECK(sum_manhattan_pairs(ps) == sm && sum_chebyshev_pairs(ps) == sc && sum_manhattan_from_each(ps) == each);
        // 中央値：距離の和が最小（小さい座標なら範囲を全部調べる）
        if (C < 10 && n > 0) {
            const Point m = manhattan_median(ps);
            auto f = [&](const Point& p) {
                long long s = 0;
                for (const Point& q : ps) s += manhattan(p, q);
                return s;
            };
            long long best = std::numeric_limits<long long>::max();
            for (Coord x = -C; x <= C; ++x)
                for (Coord y = -C; y <= C; ++y) best = std::min(best, f(Point(x, y)));
            CHECK(f(m) == best);
        }
    }
    CHECK(manhattan_median({}) == Point() && sum_manhattan_pairs({}) == 0);
}

static void test_abs_sum(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 20);
        const bool big = it % 3 == 0;
        std::vector<Coord> xs(n), ws(n);
        for (auto& x : xs) x = big ? rnd(rng, -1000000000LL, 1000000000LL) : rnd(rng, -8, 8);
        for (auto& w : ws) w = big ? rnd(rng, 0, 1000000000LL) : rnd(rng, 0, 3);
        const AbsSum S(xs), W(xs, ws);
        for (int q = 0; q < 10; ++q) {
            const Coord x = big ? rnd(rng, -1000000000LL, 1000000000LL) : rnd(rng, -10, 10);
            Int128 e = 0, ew = 0;
            for (int i = 0; i < n; ++i) e += xs[i] > x ? xs[i] - x : x - xs[i], ew += (Int128)ws[i] * (xs[i] > x ? xs[i] - x : x - xs[i]);
            CHECK(S(x) == e && W(x) == ew);
        }
        // 最小になる x（最小のもの）
        if (!big && n > 0) {
            Int128 best = -1;
            Coord bx = 0;
            for (Coord x = -10; x <= 10; ++x)
                if (best < 0 || W(x) < best) best = W(x), bx = x;
            const Coord m = W.argmin();
            CHECK(W(m) == best);
            bool any_w = false;
            for (const Coord w : ws) any_w = any_w || w > 0;
            if (any_w) CHECK(m == bx);  // 最小にする x のうち最小のもの
            Int128 bs = -1;
            for (Coord x = -10; x <= 10; ++x)
                if (bs < 0 || S(x) < bs) bs = S(x);
            CHECK(S(S.argmin()) == bs);
        }
        // 表
        const Coord lo = big ? rnd(rng, -1000000000LL, 1000000000LL - 50) : rnd(rng, -12, 12), hi = lo + rnd(rng, -1, 30);
        const std::vector<long long> tb = abs_sum_table(xs, lo, hi);
        CHECK((Coord)tb.size() == std::max<Coord>(0, hi - lo + 1));
        for (Coord x = lo; x <= hi; ++x) CHECK(tb[(size_t)(x - lo)] == (long long)S(x));
    }
    CHECK(AbsSum()(5) == 0 && AbsSum().argmin() == 0 && abs_sum_table({}, 0, 2) == std::vector<long long>(3, 0));
}

static void random_test() {
    std::mt19937_64 rng(20261008);
    int cases = 0;
    test_basic(rng, cases);
    test_abs_sum(rng, cases);
    test_farthest(rng, cases);
    test_sums(rng, cases);
    std::cout << "[manhattan] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[manhattan] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = 200000;
    const std::vector<Point> ps = rndps(rng, n, 1000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    const ManhattanFarthest F(ps);
    long long s1 = F.max_pair().d;
    for (int i = 0; i < n; ++i) s1 += F.farthest(ps[i]).d & 7;
    const long long t1 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    DynamicManhattan M;
    long long s2 = 0;
    for (int i = 0; i < n; ++i) {
        M.add(ps[i]);
        if (i % 3 == 2) M.erase(i - 1);
        s2 += M.farthest(ps[i]).d & 7;
    }
    const long long t2 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const Int128 s3 = sum_manhattan_pairs(ps) + sum_chebyshev_pairs(ps);
    const std::vector<long long> each = sum_manhattan_from_each(ps);
    const long long t3 = ms_since(t0);
    std::cout << "[manhattan] benchmark: N = " << n << "  ManhattanFarthest + farthest × N " << t1 << " ms, DynamicManhattan（追加 + 削除 + farthest）× N " << t2
              << " ms, 距離の和 3 種 " << t3 << " ms  (check " << s1 % 7 + s2 % 7 + (long long)(s3 % 7) + each[0] % 7 << ")\n";
}

}  // namespace mh_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace mh_demo {

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
#define SHOW(expr, comment) mh_demo::show(#expr, mh_demo::to_s(expr), comment)
#define RUN(...)                   \
    mh_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("マンハッタン距離 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  manhattan(a, b) / chebyshev(a, b) / in_diamond(c, D, p)                       距離 / ひし形の内側か\n"
              << "  manhattan_to_chebyshev(p) / chebyshev_to_manhattan(q)                         45 度回転（u = x + y, v = x - y）\n"
              << "  ManhattanFarthest F(ps);  F.max_pair() / F.farthest(q)                        最も遠い 2 点 / q から最も遠い点\n"
              << "  DynamicManhattan M;  M.add(p) / M.erase(id) / M.max_pair() / M.farthest(q)  点の追加・削除があるとき\n"
              << "  sum_manhattan_pairs(ps) / sum_chebyshev_pairs(ps) / sum_manhattan_from_each(ps) / manhattan_median(ps)\n"
              << "  AbsSum S(xs [, ws]): S(x) / S.argmin()  ・  abs_sum_table(xs, lo, hi)           1 次元の距離の和 Σ w|x - x_i|\n";

    api("【1】 manhattan(a, b) / chebyshev(a, b) / manhattan_to_chebyshev(p)", "Coord / Point",
        "マンハッタン距離 |dx| + |dy| / チェビシェフ距離 max(|dx|, |dy|) / 45 度回転 (x + y, x - y)", "O(1)");
    SHOW(manhattan(Point(1, 2), Point(4, -2)), "3 + 4");
    SHOW(chebyshev(Point(1, 2), Point(4, -2)), "max(3, 4)");
    SHOW(manhattan_to_chebyshev(Point(1, 2)), "(x + y, x - y)");
    SHOW(chebyshev(manhattan_to_chebyshev(Point(1, 2)), manhattan_to_chebyshev(Point(4, -2))), "回した後のチェビシェフ距離 = 元のマンハッタン距離");
    SHOW(chebyshev_to_manhattan(Point(3, -1)), "元に戻す");
    SHOW(in_diamond(Point(0, 0), 3, Point(2, 1)), "|2| + |1| ≤ 3");
    note("マンハッタン距離が D 以内の点を数える → 45 度回して「長方形の中の点」にする（2 次元の数え上げは BIT など）");

    api("【2】 ManhattanFarthest F(ps);  F.max_pair() / F.farthest(q)", "ManhattanPair{i, j, d} / ManhattanFar{i, d}",
        "最も遠い 2 点（i < j）/ q から最も遠い点（番号と距離）。u, v の最小・最大だけ覚える", "前処理 O(N)、1 回 O(1)");
    RUN(std::vector<Point> ps = {{1, 1}, {3, 4}, {-2, 0}, {5, -1}});
    RUN(ManhattanFarthest F(ps));
    RUN(ManhattanPair mp = F.max_pair());
    std::cout << "  mp.i, mp.j, mp.d                   -> " << mp.i << ", " << mp.j << ", " << mp.d << "          // (3,4) と (-2,0)：5 + 4 = 9\n";
    RUN(ManhattanFar f = F.farthest(Point(0, 0)));
    std::cout << "  f.i, f.d                           -> " << f.i << ", " << f.d << "             // (3,4)：距離 7\n";
    note("典型90 036: 各クエリ q で F.farthest(ps[q]).d");

    api("【3】 DynamicManhattan M;  M.add(p) / M.erase(id) / M.max_pair() / M.farthest(q)", "int / bool / ManhattanPair / ManhattanFar",
        "点の追加・削除がある場合（番号は add した順。消した点をもう一度消しても壊れない）", "追加・削除 O(log N)、質問 O(1)");
    RUN(DynamicManhattan M; M.add(Point(0, 0)); M.add(Point(5, 5)); M.add(Point(1, 1)));
    SHOW(M.max_pair().d, "(0,0) と (5,5)");
    RUN(M.erase(1));
    SHOW(M.max_pair().d, "(5,5) を消すと (0,0) と (1,1)");
    SHOW(M.erase(1), "もう消えている → false");

    api("【4】 sum_manhattan_pairs(ps) / sum_chebyshev_pairs(ps) / sum_manhattan_from_each(ps)", "Int128 / Int128 / vector<long long>",
        "全部の組の距離の和（大きくなるので Int128。出力は int128_to_string）/ 各点から全点への距離の和", "O(N log N)");
    SHOW(sum_manhattan_pairs(ps), "6 組の和");
    SHOW(sum_chebyshev_pairs(ps), "チェビシェフ距離の和");
    SHOW(sum_manhattan_from_each(ps), "各点から");
    note("ABC351 E（斜めに跳ぶうさぎ）: x + y の偶奇で 2 組に分け、それぞれ sum_chebyshev_pairs（偶奇が違う組は 0）");

    api("【5】 manhattan_median(ps)", "Point", "マンハッタン距離の和が最小になる格子点（x, y それぞれの中央値）", "O(N)");
    SHOW(manhattan_median(ps), "x の中央値 1、y の中央値 0");

    api("【6】 AbsSum S(xs [, ws]);  S(x) / S.argmin()  ・  abs_sum_table(xs, lo, hi)", "Int128 / Coord / vector<long long>",
        "1 次元の距離の和 f(x) = Σ w_i |x - x_i|（重みなしなら w_i = 1）/ 最小にする x / x = lo..hi の f(x) の表", "O(log N) / O(N) / O(N log N + (hi - lo))");
    RUN(AbsSum S(std::vector<Coord>{1, 4, 6}));
    SHOW(int128_to_string(S(3)), "|3-1| + |3-4| + |3-6| = 6");
    SHOW(S.argmin(), "中央値 4");
    SHOW(abs_sum_table({1, 4, 6}, 0, 6), "x = 0..6 の f(x)");
    SHOW(int128_to_string(AbsSum({0, 10}, {3, 1})(4)), "重み付き：3·4 + 1·6");
    note("ABC366 E（Σ|x - x_i| + Σ|y - y_i| ≤ D の格子点の数）: fx = abs_sum_table(xs, -2·10^6, 2·10^6)、fy も同じ。");
    note("    fy を並べて、fx[i] ごとに fy ≤ D - fx[i] の個数を二分探索（または尺取り）で足す");
}

}  // namespace mh_demo


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
    mh_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mh_test::benchmark();    // 速度計測
    mh_demo::demo();         // 使い方デモ
    return 0;
}

