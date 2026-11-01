// =============================================================================
//  ひし形・三角形・斜めのいもす（マンハッタン距離 r 以内・直角二等辺三角形・斜めの線分と帯・長方形にまとめて足す）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  H × W のマス（x が 0〜H−1、y が 0〜W−1）に「(x, y) からマンハッタン距離 r 以内（ひし形）」「斜めの線分」に
//  値を足す操作を何回もしてから、build() で各マスの値にする。
//  45 度回転 (u, w) = (x + y, x − y + W − 1) すると、ひし形も斜めの線分も長方形になるので 2 次元いもすで足せる。
//  三角形は、横の差分を「縦に流す分」と「斜めに流す分」に分けて持つ（JOI の「釘」の方法）
//    例: 「爆弾が距離 r 以内を壊す」「見張りがマンハッタン距離 r 以内を見る」を全部重ねて、各マスの回数を数える
//
//  機能一覧（座標は long long。マスの外にはみ出した部分は捨てる）                              計算量
//    ImosDiamond<T> im(H, W)            全部 0（T の既定は long long）                         O((H + W)^2)
//    add_diamond(x, y, r, v)            |i − x| + |j − y| ≤ r のマス (i, j) に v（r < 0 なら何もしない）  O(1)
//    add_square(x, y, r, v)             max(|i − x|, |j − y|) ≤ r（チェビシェフ距離）のマスに v   O(1)
//    add_rect(x1, y1, x2, y2, v)        [x1, x2) × [y1, y2) に v                                O(1)
//    add_diag(x, y, len, v)             (x + k, y + k)（k = 0〜len−1。右下向きの斜め）に v      O(1)
//    add_anti(x, y, len, v)             (x + k, y − k)（左下向きの斜め）に v                     O(1)
//    add_triangle(x, y, k, dx, dy, v)   直角の頂点 (x, y)、向き (dx, dy)（±1）、足の長さ k の直角二等辺三角形：
//                                       (x + dx·a, y + dy·b)（a, b ≥ 0、a + b ≤ k）に v（k + 1 マスの段）     O(1)
//    add_anti_band(s1, s2, v)           s1 ≤ i + j < s2 のマス全部に v（左下向きの帯）           O(1)
//    add_diag_band(d1, d2, v)           d1 ≤ i − j < d2 のマス全部に v（右下向きの帯）           O(1)
//    build()                            まとめて各マスの値にする                                O((H + W)^2)
//    get(x, y) / to_vector()            build の後：マスの値 / 全マス                           O(1) / O(HW)
//
//  注意
//    回転した表は (H + W − 1)^2 マス。H = W = 2000 の long long で 約 130 MB（int なら半分）
//    add_triangle を使った向きごとに (H + 1) × W の表を 2 枚足す（H = W = 2000 の long long で 1 つの向きにつき 約 64 MB。使わなければ作らない。build で捨てる）
//    build の後に add すると止める。get を build の前に呼んでも止める
//  速さの目安: H = W = 2000 で add_diamond を 10^6 回 + build 約 0.25 秒、add_triangle（4 つの向き）10^6 回 + build 約 0.65 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    imosdia_test::random_test()  … 1 マスずつ足す愚直解との突き合わせ
//    imosdia_test::benchmark()    … 速度計測
//    imosdia_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
namespace imosdia_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace imosdia_internal

template <class T = long long>
struct ImosDiamond {
    ImosDiamond(int H = 0, int W = 0) : H_(H), W_(W), S_(H + W - 1 > 0 ? H + W - 1 : 0) {
        if (H < 0 || W < 0) imosdia_internal::fail("ImosDiamond", "H, W は 0 以上");
        rect_.assign((size_t)(H + 1) * (W + 1), T(0));
        rot_.assign((size_t)(S_ + 1) * (S_ + 1), T(0));
    }
    void add_rect(long long x1, long long y1, long long x2, long long y2, T v) {
        need_not_built("ImosDiamond::add_rect");
        rect_add(rect_, H_, W_, x1, y1, x2, y2, v);
    }
    void add_square(long long x, long long y, long long r, T v) {
        if (r < 0) return;
        add_rect(x - r, y - r, x + r + 1, y + r + 1, v);
    }
    void add_diamond(long long x, long long y, long long r, T v) {
        need_not_built("ImosDiamond::add_diamond");
        if (r < 0) return;
        const long long u = x + y, w = x - y + (W_ - 1);
        rect_add(rot_, S_, S_, u - r, w - r, u + r + 1, w + r + 1, v);
    }
    void add_diag(long long x, long long y, long long len, T v) {  // (x + k, y + k)：w が同じで u が 2 ずつ増える
        need_not_built("ImosDiamond::add_diag");
        if (len <= 0) return;
        const long long u = x + y, w = x - y + (W_ - 1);
        rect_add(rot_, S_, S_, u, w, u + 2 * len - 1, w + 1, v);
    }
    void add_anti(long long x, long long y, long long len, T v) {  // (x + k, y − k)：u が同じで w が 2 ずつ増える
        need_not_built("ImosDiamond::add_anti");
        if (len <= 0) return;
        const long long u = x + y, w = x - y + (W_ - 1);
        rect_add(rot_, S_, S_, u, w, u + 1, w + 2 * len - 1, v);
    }
    // 直角二等辺三角形：直角の頂点 (x, y) から向き (dx, dy)（それぞれ +1 か −1）へ、(x + dx·a, y + dy·b)（a, b ≥ 0、a + b ≤ k）に v
    void add_triangle(long long x, long long y, long long k, int dx, int dy, T v) {
        need_not_built("ImosDiamond::add_triangle");
        if ((dx != 1 && dx != -1) || (dy != 1 && dy != -1)) imosdia_internal::fail("ImosDiamond::add_triangle", "dx, dy は +1 か −1");
        if (k < 0 || H_ == 0 || W_ == 0) return;
        Tri& t = tri_[(dx < 0) * 2 + (dy < 0)];  // 向きごとに、上下・左右を裏返して (+1, +1) の向きにした表
        if (t.V.empty()) t.V.assign((size_t)(H_ + 1) * W_, T(0)), t.A.assign((size_t)(H_ + 1) * W_, T(0));
        tri_add(t, dx > 0 ? x : (H_ - 1) - x, dy > 0 ? y : (W_ - 1) - y, k, v);
    }
    // 斜めの帯：x + y が [s1, s2) のマス全部 / x − y が [d1, d2) のマス全部に v
    void add_anti_band(long long s1, long long s2, T v) {
        need_not_built("ImosDiamond::add_anti_band");
        rect_add(rot_, S_, S_, s1, 0, s2, S_, v);
    }
    void add_diag_band(long long d1, long long d2, T v) {
        need_not_built("ImosDiamond::add_diag_band");
        rect_add(rot_, S_, S_, 0, d1 + (W_ - 1), S_, d2 + (W_ - 1), v);
    }
    void build() {
        need_not_built("ImosDiamond::build");
        prefix(rect_, H_, W_);
        prefix(rot_, S_, S_);
        for (int o = 0; o < 4; ++o) {
            Tri& t = tri_[o];
            if (t.V.empty()) continue;
            const size_t w = W_;
            for (int i = 1; i < H_; ++i)
                for (int j = 0; j < W_; ++j) {
                    t.V[i * w + j] += t.V[(i - 1) * w + j];                     // 縦に流す
                    if (j + 1 < W_) t.A[i * w + j] += t.A[(i - 1) * w + j + 1];  // 左下へ斜めに流す
                }
            for (int i = 0; i < H_; ++i) {
                T run = T(0);
                const int x = (o & 2) ? H_ - 1 - i : i;
                for (int j = 0; j < W_; ++j) {
                    run += t.V[i * w + j] + t.A[i * w + j];  // 横に流す
                    const int y = (o & 1) ? W_ - 1 - j : j;
                    rect_[(size_t)x * (W_ + 1) + y] += run;
                }
            }
            std::vector<T>().swap(t.V), std::vector<T>().swap(t.A);
        }
        built_ = true;
    }
    T get(int x, int y) const {
        if (!built_) imosdia_internal::fail("ImosDiamond::get", "先に build() を呼ぶ");
        if (x < 0 || x >= H_ || y < 0 || y >= W_) imosdia_internal::fail("ImosDiamond::get", "(x, y) が範囲外");
        const int u = x + y, w = x - y + (W_ - 1);
        return rect_[(size_t)x * (W_ + 1) + y] + rot_[(size_t)u * (S_ + 1) + w];
    }
    std::vector<std::vector<T>> to_vector() const {
        std::vector<std::vector<T>> g(H_, std::vector<T>(W_));
        for (int x = 0; x < H_; ++x)
            for (int y = 0; y < W_; ++y) g[x][y] = get(x, y);
        return g;
    }

   private:
    int H_, W_, S_;
    std::vector<T> rect_, rot_;  // 普通の差分 (H+1)×(W+1)、回転した差分 (S+1)×(S+1)
    struct Tri {
        std::vector<T> V, A;  // 横の差分を、縦に流す分 V と斜めに流す分 A に分けて持つ（(H+1)×W）
    };
    Tri tri_[4];
    bool built_ = false;
    // 向き (+1, +1) の三角形 i ≥ x, j ≥ y, i + j ≤ s（s = x + y + k）。行 i の横の差分は「列 y に +v、列 s − i + 1 に −v」。
    // +v は列 y を縦に、−v は (i, s − i + 1) を斜めに並ぶので、縦・斜めの差分にする。列が負の分は列 0 に、W 以上は捨てる
    void tri_add(Tri& t, long long x, long long y, long long k, T v) {
        const long long s = x + y + k, i1 = std::max(x, 0LL), i2 = std::min(x + k, (long long)H_ - 1);
        if (i1 > i2) return;
        auto at = [&](std::vector<T>& d, long long i, long long j) -> T& { return d[(size_t)i * W_ + (size_t)j]; };
        if (y < W_) {
            const long long cy = std::max(y, 0LL);
            at(t.V, i1, cy) += v, at(t.V, i2 + 1, cy) -= v;
        }
        const long long a1 = std::max(i1, s + 2 - W_), a2 = std::min(i2, s + 1);  // −v の列が [0, W) の行
        if (a1 <= a2) {
            at(t.A, a1, s - a1 + 1) -= v;
            if (s - a2 >= 0) at(t.A, a2 + 1, s - a2) += v;
        }
        const long long b1 = std::max(i1, s + 2);  // −v の列が負の行 → 列 0
        if (b1 <= i2) at(t.V, b1, 0) -= v, at(t.V, i2 + 1, 0) += v;
    }
    void need_not_built(const char* f) const {
        if (built_) imosdia_internal::fail(f, "build() の後には足せない");
    }
    static void rect_add(std::vector<T>& d, int h, int w, long long x1, long long y1, long long x2, long long y2, T v) {
        x1 = std::max(x1, 0LL), y1 = std::max(y1, 0LL), x2 = std::min(x2, (long long)h), y2 = std::min(y2, (long long)w);
        if (x1 >= x2 || y1 >= y2) return;
        const size_t s = (size_t)w + 1;
        d[x1 * s + y1] += v, d[x1 * s + y2] -= v, d[x2 * s + y1] -= v, d[x2 * s + y2] += v;
    }
    static void prefix(std::vector<T>& d, int h, int w) {
        const size_t s = (size_t)w + 1;
        for (int x = 0; x <= h; ++x)
            for (int y = 1; y <= w; ++y) d[x * s + y] += d[x * s + y - 1];
        for (int x = 1; x <= h; ++x)
            for (int y = 0; y <= w; ++y) d[x * s + y] += d[(x - 1) * s + y];
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace imosdia_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define ID_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[pre_imos_diamond] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 800 * TEST_SCALE; ++it) {
        const int H = (int)rnd(rng, 0, 9), W = (int)rnd(rng, 0, 9);
        ImosDiamond<ll> im(H, W);
        std::vector<std::vector<ll>> g(H, std::vector<ll>(W, 0));
        auto put = [&](ll i, ll j, ll v) {
            if (0 <= i && i < H && 0 <= j && j < W) g[i][j] += v;
        };
        const int ops = (int)rnd(rng, 0, 25);
        for (int t = 0; t < ops; ++t) {
            const ll x = rnd(rng, -4, H + 4), y = rnd(rng, -4, W + 4), r = rnd(rng, -1, 8), v = rnd(rng, -5, 5);
            const int kind = (int)rnd(rng, 0, 7);
            if (kind == 0) {
                im.add_diamond(x, y, r, v);
                for (ll i = -15; i < 25; ++i)
                    for (ll j = -15; j < 25; ++j)
                        if (std::abs(i - x) + std::abs(j - y) <= r) put(i, j, v);
            } else if (kind == 1) {
                im.add_square(x, y, r, v);
                for (ll i = x - r; i <= x + r; ++i)
                    for (ll j = y - r; j <= y + r; ++j) put(i, j, v);
            } else if (kind == 2) {
                const ll x2 = rnd(rng, -4, H + 4), y2 = rnd(rng, -4, W + 4);
                im.add_rect(x, y, x2, y2, v);
                for (ll i = x; i < x2; ++i)
                    for (ll j = y; j < y2; ++j) put(i, j, v);
            } else if (kind == 3) {
                im.add_diag(x, y, r, v);
                for (ll k = 0; k < r; ++k) put(x + k, y + k, v);
            } else if (kind == 4) {
                im.add_anti(x, y, r, v);
                for (ll k = 0; k < r; ++k) put(x + k, y - k, v);
            } else if (kind == 5) {
                const int dx = rnd(rng, 0, 1) ? 1 : -1, dy = rnd(rng, 0, 1) ? 1 : -1;
                const ll k = rnd(rng, -1, 12);
                im.add_triangle(x, y, k, dx, dy, v);
                for (ll a = 0; a <= k; ++a)
                    for (ll b = 0; a + b <= k; ++b) put(x + dx * a, y + dy * b, v);
            } else if (kind == 6) {
                const ll s1 = rnd(rng, -3, H + W + 2), s2 = rnd(rng, -3, H + W + 2);
                im.add_anti_band(s1, s2, v);
                for (ll i = 0; i < H; ++i)
                    for (ll j = 0; j < W; ++j)
                        if (s1 <= i + j && i + j < s2) put(i, j, v);
            } else {
                const ll d1 = rnd(rng, -W - 2, H + 2), d2 = rnd(rng, -W - 2, H + 2);
                im.add_diag_band(d1, d2, v);
                for (ll i = 0; i < H; ++i)
                    for (ll j = 0; j < W; ++j)
                        if (d1 <= i - j && i - j < d2) put(i, j, v);
            }
        }
        im.build();
        ID_CHECK(im.to_vector() == g);
        cases += (long long)H * W + 1;
    }
    std::cout << "[pre_imos_diamond] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[pre_imos_diamond] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int H = TEST_SCALE == 1 ? 400 : 2000, Q = TEST_SCALE == 1 ? 100000 : 1000000;
    const auto t0 = std::chrono::steady_clock::now();
    ImosDiamond<ll> im(H, H);
    for (int i = 0; i < Q; ++i) im.add_diamond(rnd(rng, 0, H - 1), rnd(rng, 0, H - 1), rnd(rng, 0, H / 4), 1);
    im.build();
    ll chk = 0;
    for (int i = 0; i < H; i += 7) chk += im.get(i, H - 1 - i);
    std::cout << "[pre_imos_diamond] benchmark (H = W = " << H << "、add_diamond " << Q << " 回 + build): " << elapsed_ms(t0) << " ms  (check "
              << chk % 1000 << ")\n";
    const auto t1 = std::chrono::steady_clock::now();
    ImosDiamond<ll> tr(H, H);
    for (int i = 0; i < Q; ++i) tr.add_triangle(rnd(rng, 0, H - 1), rnd(rng, 0, H - 1), rnd(rng, 0, H / 4), rnd(rng, 0, 1) ? 1 : -1, rnd(rng, 0, 1) ? 1 : -1, 1);
    tr.build();
    chk = 0;
    for (int i = 0; i < H; i += 7) chk += tr.get(i, H - 1 - i);
    std::cout << "[pre_imos_diamond] benchmark (H = W = " << H << "、add_triangle " << Q << " 回 + build): " << elapsed_ms(t1) << " ms  (check "
              << chk % 1000 << ")\n";
}

}  // namespace imosdia_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace imosdia_demo {

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
#define SHOW(expr, comment) imosdia_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    imosdia_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void show_grid(const std::vector<std::vector<long long>>& g) {
    for (const auto& r : g) {
        std::string s = "      ";
        for (long long v : r) s += std::to_string(v) + " ";
        std::cout << s << "\n";
    }
}

static void demo() {
    banner("ひし形・斜めのいもす（pre_imos_diamond） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  ImosDiamond<T> im(H, W)   add_diamond / add_square / add_rect / add_diag / add_anti / add_triangle / add_anti_band / add_diag_band\n"
              << "                            → build → get / to_vector\n";

    api("【1】 add_diamond(x, y, r, v) / build() / get(x, y)", "void / void / T", "マンハッタン距離 r 以内のマスに v を足す（はみ出しは捨てる）",
        "O(1) / O((H + W)^2) / O(1)");
    RUN(ImosDiamond<> im(5, 7));
    RUN(im.add_diamond(2, 3, 2, 1));
    RUN(im.add_diamond(4, 0, 1, 10));  // 下の端（はみ出す）
    RUN(im.build());
    SHOW(im.get(2, 3), "中心");
    SHOW(im.get(0, 3), "距離 2");
    SHOW(im.get(0, 2), "距離 3 なので入らない");
    note("im.to_vector() =");
    show_grid(im.to_vector());

    api("【2】 add_square(x, y, r, v) / add_rect(x1, y1, x2, y2, v)", "void", "チェビシェフ距離 r 以内（正方形）/ [x1, x2) × [y1, y2)", "O(1)");
    RUN(ImosDiamond<> sq(4, 4));
    RUN(sq.add_square(1, 1, 1, 1));
    RUN(sq.add_rect(2, 0, 4, 4, 5));
    RUN(sq.build());
    note("sq.to_vector() =");
    show_grid(sq.to_vector());

    api("【3】 add_diag(x, y, len, v) / add_anti(x, y, len, v)", "void", "(x + k, y + k) / (x + k, y − k)（k = 0〜len−1）に v。斜めの線分", "O(1)");
    RUN(ImosDiamond<> dg(4, 5));
    RUN(dg.add_diag(0, 0, 4, 1));
    RUN(dg.add_anti(0, 4, 10, 2));  // はみ出した部分は捨てる
    RUN(dg.build());
    note("dg.to_vector() =");
    show_grid(dg.to_vector());

    api("【4】 add_triangle(x, y, k, dx, dy, v)", "void", "直角の頂点 (x, y) から向き (dx, dy) に、足の長さ k の直角二等辺三角形に v", "O(1)");
    RUN(ImosDiamond<> tg(5, 6));
    RUN(tg.add_triangle(0, 0, 3, +1, +1, 1));  // 左上が直角
    RUN(tg.add_triangle(4, 5, 2, -1, -1, 5));  // 右下が直角（上・左へ）
    RUN(tg.build());
    note("tg.to_vector() =");
    show_grid(tg.to_vector());

    api("【5】 add_anti_band(s1, s2, v) / add_diag_band(d1, d2, v)", "void", "i + j が [s1, s2) / i − j が [d1, d2) のマス全部に v（斜めの帯）", "O(1)");
    RUN(ImosDiamond<> bd(4, 5));
    RUN(bd.add_anti_band(2, 4, 1));
    RUN(bd.add_diag_band(0, 1, 10));
    RUN(bd.build());
    note("bd.to_vector() =");
    show_grid(bd.to_vector());
}

}  // namespace imosdia_demo

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
    imosdia_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    imosdia_test::benchmark();    // 速度計測
    imosdia_demo::demo();         // 使い方デモ
    return 0;
}
