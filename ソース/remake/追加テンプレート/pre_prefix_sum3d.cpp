// =============================================================================
//  3 次元累積和・3 次元いもす（直方体の和を O(1)・直方体にまとめて足す）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  X × Y × Z の 3 次元配列で、直方体 [x1, x2) × [y1, y2) × [z1, z2) の和を包除（8 項）で O(1) に答える。
//  逆に「直方体に足す」を何回もしてから各マスの値にするのが 3 次元いもす
//    例: ABC366 D（Cuboid Sum Query）… N ≤ 100 の立方体で直方体の和を Q 回 → PrefixSum3D の sum_closed
//
//  機能一覧（0-indexed・半開区間。閉区間は _closed 版。範囲のはみ出しは捨てる）                 計算量
//    PrefixSum3D<T> ps(a)               a[x][y][z]（vector<vector<vector<…>>>）から作る            O(XYZ)
//    ps.sum(x1, y1, z1, x2, y2, z2)     [x1, x2) × [y1, y2) × [z1, z2) の和                       O(1)
//    ps.sum_closed(x1, y1, z1, x2, y2, z2)   [x1, x2] × [y1, y2] × [z1, z2] の和
//    Imos3D<T> im(X, Y, Z)              全部 0                                                     O(XYZ)
//    im.add(x1, y1, z1, x2, y2, z2, v)  直方体 [x1, x2) × … に v を足す（予約）                     O(1)
//    im.add_closed(…, v)                閉区間版
//    im.build() / im.get(x, y, z)       各マスの値にする / マスの値                                 O(XYZ) / O(1)
//    im.to_prefix_sum()                 build の後：その値の PrefixSum3D を作る（足してから直方体の和を聞く）  O(XYZ)
//
//  注意: メモリは (X + 1)(Y + 1)(Z + 1) 個の T。X = Y = Z = 400 の long long で 約 0.5 GB なので、大きいときは int を使う
//  速さの目安: X = Y = Z = 100 で PrefixSum3D の構築 + sum を 2×10^5 回 約 0.03 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    ps3d_test::random_test()  … 全部足す愚直解との突き合わせ
//    ps3d_test::benchmark()    … 速度計測
//    ps3d_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace ps3d_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace ps3d_internal

template <class T = long long>
struct PrefixSum3D {
    PrefixSum3D() {}
    template <class A>
    explicit PrefixSum3D(const std::vector<std::vector<std::vector<A>>>& a) {
        X_ = (int)a.size(), Y_ = X_ ? (int)a[0].size() : 0, Z_ = (X_ && Y_) ? (int)a[0][0].size() : 0;
        s_.assign((size_t)(X_ + 1) * (Y_ + 1) * (Z_ + 1), T(0));
        for (int x = 0; x < X_; ++x) {
            if ((int)a[x].size() != Y_) ps3d_internal::fail("PrefixSum3D", "大きさがそろっていない");
            for (int y = 0; y < Y_; ++y) {
                if ((int)a[x][y].size() != Z_) ps3d_internal::fail("PrefixSum3D", "大きさがそろっていない");
                for (int z = 0; z < Z_; ++z) at(x + 1, y + 1, z + 1) = T(a[x][y][z]);
            }
        }
        build_prefix();
    }
    T sum(long long x1, long long y1, long long z1, long long x2, long long y2, long long z2) const {
        x1 = std::max(x1, 0LL), y1 = std::max(y1, 0LL), z1 = std::max(z1, 0LL);
        x2 = std::min(x2, (long long)X_), y2 = std::min(y2, (long long)Y_), z2 = std::min(z2, (long long)Z_);
        if (x1 >= x2 || y1 >= y2 || z1 >= z2) return T(0);
        return cat(x2, y2, z2) - cat(x1, y2, z2) - cat(x2, y1, z2) - cat(x2, y2, z1) + cat(x1, y1, z2) + cat(x1, y2, z1) + cat(x2, y1, z1) -
               cat(x1, y1, z1);
    }
    T sum_closed(long long x1, long long y1, long long z1, long long x2, long long y2, long long z2) const {
        return sum(x1, y1, z1, x2 + 1, y2 + 1, z2 + 1);
    }

   private:
    template <class>
    friend struct Imos3D;
    int X_ = 0, Y_ = 0, Z_ = 0;
    std::vector<T> s_;  // s_[x][y][z] = [0, x) × [0, y) × [0, z) の和
    size_t id(long long x, long long y, long long z) const { return ((size_t)x * (Y_ + 1) + (size_t)y) * (Z_ + 1) + (size_t)z; }
    T& at(long long x, long long y, long long z) { return s_[id(x, y, z)]; }
    T cat(long long x, long long y, long long z) const { return s_[id(x, y, z)]; }
    void build_prefix() {  // 3 方向の累積和
        for (int x = 0; x <= X_; ++x)
            for (int y = 0; y <= Y_; ++y)
                for (int z = 1; z <= Z_; ++z) at(x, y, z) += at(x, y, z - 1);
        for (int x = 0; x <= X_; ++x)
            for (int y = 1; y <= Y_; ++y)
                for (int z = 0; z <= Z_; ++z) at(x, y, z) += at(x, y - 1, z);
        for (int x = 1; x <= X_; ++x)
            for (int y = 0; y <= Y_; ++y)
                for (int z = 0; z <= Z_; ++z) at(x, y, z) += at(x - 1, y, z);
    }
};

template <class T = long long>
struct Imos3D {
    Imos3D(int X = 0, int Y = 0, int Z = 0) : X_(X), Y_(Y), Z_(Z), d_((size_t)(X + 1) * (Y + 1) * (Z + 1), T(0)) {
        if (X < 0 || Y < 0 || Z < 0) ps3d_internal::fail("Imos3D", "大きさは 0 以上");
    }
    void add(long long x1, long long y1, long long z1, long long x2, long long y2, long long z2, T v) {
        if (built_) ps3d_internal::fail("Imos3D::add", "build() の後には足せない");
        x1 = std::max(x1, 0LL), y1 = std::max(y1, 0LL), z1 = std::max(z1, 0LL);
        x2 = std::min(x2, (long long)X_), y2 = std::min(y2, (long long)Y_), z2 = std::min(z2, (long long)Z_);
        if (x1 >= x2 || y1 >= y2 || z1 >= z2) return;
        at(x1, y1, z1) += v, at(x2, y1, z1) -= v, at(x1, y2, z1) -= v, at(x1, y1, z2) -= v;
        at(x2, y2, z1) += v, at(x2, y1, z2) += v, at(x1, y2, z2) += v, at(x2, y2, z2) -= v;
    }
    void add_closed(long long x1, long long y1, long long z1, long long x2, long long y2, long long z2, T v) {
        add(x1, y1, z1, x2 + 1, y2 + 1, z2 + 1, v);
    }
    void build() {
        if (built_) ps3d_internal::fail("Imos3D::build", "build() は 1 回だけ");
        for (int x = 0; x <= X_; ++x)
            for (int y = 0; y <= Y_; ++y)
                for (int z = 1; z <= Z_; ++z) at(x, y, z) += at(x, y, z - 1);
        for (int x = 0; x <= X_; ++x)
            for (int y = 1; y <= Y_; ++y)
                for (int z = 0; z <= Z_; ++z) at(x, y, z) += at(x, y - 1, z);
        for (int x = 1; x <= X_; ++x)
            for (int y = 0; y <= Y_; ++y)
                for (int z = 0; z <= Z_; ++z) at(x, y, z) += at(x - 1, y, z);
        built_ = true;
    }
    T get(int x, int y, int z) const {
        if (!built_) ps3d_internal::fail("Imos3D::get", "先に build() を呼ぶ");
        if (x < 0 || x >= X_ || y < 0 || y >= Y_ || z < 0 || z >= Z_) ps3d_internal::fail("Imos3D::get", "範囲外");
        return d_[id(x, y, z)];
    }
    PrefixSum3D<T> to_prefix_sum() const {
        if (!built_) ps3d_internal::fail("Imos3D::to_prefix_sum", "先に build() を呼ぶ");
        PrefixSum3D<T> ps;
        ps.X_ = X_, ps.Y_ = Y_, ps.Z_ = Z_;
        ps.s_.assign(d_.size(), T(0));
        for (int x = 0; x < X_; ++x)
            for (int y = 0; y < Y_; ++y)
                for (int z = 0; z < Z_; ++z) ps.at(x + 1, y + 1, z + 1) = d_[id(x, y, z)];
        ps.build_prefix();
        return ps;
    }

   private:
    int X_, Y_, Z_;
    std::vector<T> d_;
    bool built_ = false;
    size_t id(long long x, long long y, long long z) const { return ((size_t)x * (Y_ + 1) + (size_t)y) * (Z_ + 1) + (size_t)z; }
    T& at(long long x, long long y, long long z) { return d_[id(x, y, z)]; }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace ps3d_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define PS_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[pre_prefix_sum3d] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                      \
        }                                                                                      \
    } while (0)

// [lo, hi] の一様乱数
[[maybe_unused]] static long long rnd(std::mt19937_64& rng, long long lo, long long hi) {
    const unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return (long long)((unsigned long long)lo + (w == 0 ? rng() : rng() % w));
}
[[maybe_unused]] static long long elapsed_ms(std::chrono::steady_clock::time_point a) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
}

using ll = long long;
using V3 = std::vector<std::vector<std::vector<ll>>>;

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 600 * TEST_SCALE; ++it) {
        const int X = (int)rnd(rng, 0, 5), Y = (int)rnd(rng, 0, 5), Z = (int)rnd(rng, 0, 5);
        V3 a(X, std::vector<std::vector<ll>>(Y, std::vector<ll>(Z)));
        for (auto& p : a)
            for (auto& q : p)
                for (auto& v : q) v = rnd(rng, -9, 9);
        auto box = [&](const V3& g, ll x1, ll y1, ll z1, ll x2, ll y2, ll z2) {
            ll s = 0;
            for (ll x = std::max(0LL, x1); x < std::min<ll>(X, x2); ++x)
                for (ll y = std::max(0LL, y1); y < std::min<ll>(Y, y2); ++y)
                    for (ll z = std::max(0LL, z1); z < std::min<ll>(Z, z2); ++z) s += g[x][y][z];
            return s;
        };
        const PrefixSum3D<ll> ps(a);
        for (int q = 0; q < 30; ++q) {
            const ll x1 = rnd(rng, -1, X + 1), x2 = rnd(rng, -1, X + 1), y1 = rnd(rng, -1, Y + 1), y2 = rnd(rng, -1, Y + 1), z1 = rnd(rng, -1, Z + 1),
                     z2 = rnd(rng, -1, Z + 1);
            PS_CHECK(ps.sum(x1, y1, z1, x2, y2, z2) == box(a, x1, y1, z1, x2, y2, z2));
            PS_CHECK(ps.sum_closed(x1, y1, z1, x2, y2, z2) == box(a, x1, y1, z1, x2 + 1, y2 + 1, z2 + 1));
            cases += 2;
        }
        Imos3D<ll> im(X, Y, Z);
        V3 g(X, std::vector<std::vector<ll>>(Y, std::vector<ll>(Z, 0)));
        for (int t = 0; t < 15; ++t) {
            ll x1 = rnd(rng, -1, X + 1), x2 = rnd(rng, -1, X + 1), y1 = rnd(rng, -1, Y + 1), y2 = rnd(rng, -1, Y + 1), z1 = rnd(rng, -1, Z + 1),
               z2 = rnd(rng, -1, Z + 1);
            const ll v = rnd(rng, -5, 5);
            const bool closed = rnd(rng, 0, 1);
            if (closed) im.add_closed(x1, y1, z1, x2, y2, z2, v), ++x2, ++y2, ++z2;
            else im.add(x1, y1, z1, x2, y2, z2, v);
            for (ll x = std::max(0LL, x1); x < std::min<ll>(X, x2); ++x)
                for (ll y = std::max(0LL, y1); y < std::min<ll>(Y, y2); ++y)
                    for (ll z = std::max(0LL, z1); z < std::min<ll>(Z, z2); ++z) g[x][y][z] += v;
        }
        im.build();
        for (int x = 0; x < X; ++x)
            for (int y = 0; y < Y; ++y)
                for (int z = 0; z < Z; ++z) PS_CHECK(im.get(x, y, z) == g[x][y][z]);
        const auto ps2 = im.to_prefix_sum();
        PS_CHECK(ps2.sum(0, 0, 0, X, Y, Z) == box(g, 0, 0, 0, X, Y, Z));
        cases += (long long)X * Y * Z + 1;
    }
    std::cout << "[pre_prefix_sum3d] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[pre_prefix_sum3d] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 40 : 100, Q = 200000;
    V3 a(N, std::vector<std::vector<ll>>(N, std::vector<ll>(N)));
    for (auto& p : a)
        for (auto& q : p)
            for (auto& v : q) v = rnd(rng, 0, 999);
    const auto t0 = std::chrono::steady_clock::now();
    const PrefixSum3D<ll> ps(a);
    ll chk = 0;
    for (int i = 0; i < Q; ++i) {
        ll l[3], r[3];
        for (int k = 0; k < 3; ++k) {
            l[k] = rnd(rng, 0, N - 1), r[k] = rnd(rng, 0, N - 1);
            if (l[k] > r[k]) std::swap(l[k], r[k]);
        }
        chk += ps.sum_closed(l[0], l[1], l[2], r[0], r[1], r[2]);
    }
    std::cout << "[pre_prefix_sum3d] benchmark (N = " << N << "、構築 + sum_closed " << Q << " 回): " << elapsed_ms(t0) << " ms  (check " << chk % 1000
              << ")\n";
}

}  // namespace ps3d_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace ps3d_demo {

[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(const std::string& s) { return "\"" + s + "\""; }
[[maybe_unused]] static std::string to_s(const char* s) { return to_s(std::string(s)); }
[[maybe_unused]] static std::string to_s(char c) { return std::string("'") + c + "'"; }
[[maybe_unused]] static std::string to_s(double v) {
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}
template <class T>
static std::string to_s(const T& v);
template <class A, class B>
static std::string to_s(const std::pair<A, B>& p);
template <class T>
static std::string to_s(const std::vector<T>& v);
template <class T>
static std::string to_s(const std::optional<T>& v);
template <class T>
static std::string to_s(const T& v) {  // << で出せるもの（整数など）
    std::ostringstream os;
    os << v;
    return os.str();
}
template <class A, class B>
static std::string to_s(const std::pair<A, B>& p) {
    return "(" + to_s(p.first) + ", " + to_s(p.second) + ")";
}
template <class T>
static std::string to_s(const std::vector<T>& v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "]";
}
template <class T>
static std::string to_s(const std::optional<T>& v) {
    return v ? to_s(*v) : std::string("nullopt");
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
[[maybe_unused]] static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
[[maybe_unused]] static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
[[maybe_unused]] static void stmt(const std::string& s) { code(s + ";"); }
template <class V>
[[maybe_unused]] static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 38) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 66) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) ps3d_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    ps3d_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("3 次元累積和・3 次元いもす（pre_prefix_sum3d） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  PrefixSum3D<T> ps(a): sum / sum_closed\n"
              << "  Imos3D<T> im(X, Y, Z): add / add_closed → build → get / to_prefix_sum\n";

    api("【1】 PrefixSum3D<> ps(a) / sum(x1, y1, z1, x2, y2, z2) / sum_closed(…)", "— / T / T",
        "直方体 [x1, x2) × [y1, y2) × [z1, z2) の和（閉区間は _closed）", "O(XYZ) / O(1)");
    RUN(std::vector<std::vector<std::vector<long long>>> a = {{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}});  // 2 × 2 × 2
    RUN(const PrefixSum3D<> ps(a));
    SHOW(ps.sum(0, 0, 0, 2, 2, 2), "全部 = 36");
    SHOW(ps.sum(1, 0, 0, 2, 2, 2), "x = 1 の面 = 5 + 6 + 7 + 8");
    SHOW(ps.sum_closed(0, 1, 1, 1, 1, 1), "(0,1,1) と (1,1,1) = 4 + 8");
    note("ABC366 D（1-indexed の閉区間 [Lx, Rx] …）は ps.sum_closed(Lx − 1, Ly − 1, Lz − 1, Rx − 1, Ry − 1, Rz − 1)");

    api("【2】 Imos3D<> im(X, Y, Z) / add(…, v) / build() / get(x, y, z)", "— / void / void / T", "直方体にまとめて足してから各マスの値にする",
        "O(XYZ) / O(1) / O(XYZ) / O(1)");
    RUN(Imos3D<> im(3, 3, 3));
    RUN(im.add(0, 0, 0, 2, 2, 2, 1));
    RUN(im.add_closed(1, 1, 1, 2, 2, 2, 10));
    RUN(im.build());
    SHOW(im.get(0, 0, 0), "1 つ目だけ");
    SHOW(im.get(1, 1, 1), "両方");
    SHOW(im.get(2, 2, 2), "2 つ目だけ");
    SHOW(im.to_prefix_sum().sum(0, 0, 0, 3, 3, 3), "全体の和 = 8 × 1 + 8 × 10");
}

}  // namespace ps3d_demo

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
    ps3d_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ps3d_test::benchmark();    // 速度計測
    ps3d_demo::demo();         // 使い方デモ
    return 0;
}
