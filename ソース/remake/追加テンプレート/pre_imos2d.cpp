// =============================================================================
//  2 次元いもす（2 次元の差分配列：長方形にまとめて足してから、各マスの値にする）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  H × W のマス（x が 1 つ目で 0〜H−1、y が 2 つ目で 0〜W−1）に、長方形ごとに値を足す操作を何回もしてから、
//  build() で各マスの値にまとめる。1 回の add は O(1)、build は O(HW)
//    例: 典型90 028（Cluttered Paper）… 紙ごとに add(lx, ly, rx, ry, 1) → build → 値が k のマスを数える
//  長方形がマスの外にはみ出していたら、はみ出した部分は捨てる（座標は long long で渡してよい）
//
//  機能一覧（長方形は半開区間 [x1, x2) × [y1, y2)。閉区間は _closed 版）              計算量
//    Imos2D<T> im(H, W)                    全部 0 の H × W（T の既定は long long）            O(HW)
//    add(x1, y1, x2, y2, v)                [x1, x2) × [y1, y2) の全マスに v を足す（予約）     O(1)
//    add_closed(x1, y1, x2, y2, v)         [x1, x2] × [y1, y2] の全マスに v を足す（予約）     O(1)
//    add_point(x, y, v)                    1 マス (x, y) に v を足す（予約）                   O(1)
//    build()                               予約した add をまとめて、各マスの値にする           O(HW)
//    --- 以下は build() の後 ---
//    get(x, y)                             マス (x, y) の値                                   O(1)
//    to_vector()                           全マスの値（vector<vector<T>>）                     O(HW)
//    count_if(f)                           f(値) が true になるマスの数                        O(HW)
//    max_value() / min_value()             全マスの最大値 / 最小値                             O(HW)
//    build_sum()                           長方形の和を答える表を作る（メモリが 2 倍になる）   O(HW)
//    sum(x1, y1, x2, y2) / sum_closed(…)   build_sum の後：長方形の値の和                     O(1)
//    --- いつでも ---
//    height() / width() / is_built()       H / W / build 済みか
//    clear()                               全部 0・build 前に戻す                              O(HW)
//
//  注意
//    build() の後に add すると止まる（もう一度足すなら clear() してからやり直す）
//    get・sum などを build() の前に呼ぶと止まる。範囲外の (x, y) も止まる
//    T は long long（ふつう）/ int（マスが多くてメモリを減らしたいとき）/ double・mint など足し算できる型
//  速さの目安: H = W = 2000、add を 10^6 回、build、build_sum で 合計 約 0.07 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    imos2d_test::random_test()  … 1 マスずつ足す愚直解との突き合わせ
//    imos2d_test::benchmark()    … 速度計測
//    imos2d_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
namespace imos2d_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace imos2d_internal

template <class T = long long>
struct Imos2D {
    Imos2D(int H = 0, int W = 0) : H_(H), W_(W), d_((size_t)(H + 1) * (W + 1), T(0)) {
        if (H < 0 || W < 0) imos2d_internal::fail("Imos2D", "H, W は 0 以上");
    }
    int height() const { return H_; }
    int width() const { return W_; }
    bool is_built() const { return built_; }

    // [x1, x2) × [y1, y2) に v を足す（はみ出した部分は捨てる）
    void add(long long x1, long long y1, long long x2, long long y2, T v) {
        if (built_) imos2d_internal::fail("Imos2D::add", "build() の後には足せない（clear() してからやり直す）");
        x1 = std::max(x1, 0LL), y1 = std::max(y1, 0LL);
        x2 = std::min(x2, (long long)H_), y2 = std::min(y2, (long long)W_);
        if (x1 >= x2 || y1 >= y2) return;
        at(x1, y1) += v;
        at(x1, y2) -= v;
        at(x2, y1) -= v;
        at(x2, y2) += v;
    }
    // [x1, x2] × [y1, y2] に v を足す
    void add_closed(long long x1, long long y1, long long x2, long long y2, T v) {
        if (x2 < x1 || y2 < y1) return;
        add(x1, y1, x2 + 1, y2 + 1, v);
    }
    void add_point(long long x, long long y, T v) { add(x, y, x + 1, y + 1, v); }

    // 各マスの値にする（2 次元の累積和）
    void build() {
        if (built_) imos2d_internal::fail("Imos2D::build", "build() は 1 回だけ（足し直すなら clear()）");
        const size_t w = (size_t)W_ + 1;
        for (int x = 0; x <= H_; ++x) {
            T* r = &d_[(size_t)x * w];
            for (int y = 1; y <= W_; ++y) r[y] += r[y - 1];
        }
        for (int x = 1; x <= H_; ++x) {
            T* r = &d_[(size_t)x * w];
            const T* p = &d_[(size_t)(x - 1) * w];
            for (int y = 0; y <= W_; ++y) r[y] += p[y];
        }
        built_ = true;
    }
    T get(int x, int y) const {
        need_built("Imos2D::get");
        if (x < 0 || x >= H_ || y < 0 || y >= W_) imos2d_internal::fail("Imos2D::get", "(x, y) が範囲外");
        return d_[(size_t)x * (W_ + 1) + y];
    }
    std::vector<std::vector<T>> to_vector() const {
        need_built("Imos2D::to_vector");
        std::vector<std::vector<T>> g(H_, std::vector<T>(W_));
        for (int x = 0; x < H_; ++x)
            for (int y = 0; y < W_; ++y) g[x][y] = d_[(size_t)x * (W_ + 1) + y];
        return g;
    }
    template <class F>
    long long count_if(F f) const {
        need_built("Imos2D::count_if");
        long long c = 0;
        for (int x = 0; x < H_; ++x)
            for (int y = 0; y < W_; ++y) c += f(d_[(size_t)x * (W_ + 1) + y]) ? 1 : 0;
        return c;
    }
    T max_value() const { return extreme(true); }
    T min_value() const { return extreme(false); }

    // 長方形の和を答える表を作る（build の後）
    void build_sum() {
        need_built("Imos2D::build_sum");
        const size_t w = (size_t)W_ + 1;
        s_.assign((size_t)(H_ + 1) * w, T(0));
        for (int x = 0; x < H_; ++x) {
            T row = T(0);
            for (int y = 0; y < W_; ++y) {
                row += d_[(size_t)x * w + y];
                s_[(size_t)(x + 1) * w + y + 1] = s_[(size_t)x * w + y + 1] + row;
            }
        }
        has_sum_ = true;
    }
    // [x1, x2) × [y1, y2) の値の和（はみ出しは捨てる）
    T sum(long long x1, long long y1, long long x2, long long y2) const {
        if (!has_sum_) imos2d_internal::fail("Imos2D::sum", "先に build() と build_sum() を呼ぶ");
        x1 = std::max(x1, 0LL), y1 = std::max(y1, 0LL);
        x2 = std::min(x2, (long long)H_), y2 = std::min(y2, (long long)W_);
        if (x1 >= x2 || y1 >= y2) return T(0);
        const size_t w = (size_t)W_ + 1;
        return s_[x2 * w + y2] - s_[x1 * w + y2] - s_[x2 * w + y1] + s_[x1 * w + y1];
    }
    T sum_closed(long long x1, long long y1, long long x2, long long y2) const {
        if (x2 < x1 || y2 < y1) return T(0);
        return sum(x1, y1, x2 + 1, y2 + 1);
    }
    void clear() {
        std::fill(d_.begin(), d_.end(), T(0));
        s_.clear(), s_.shrink_to_fit();
        built_ = has_sum_ = false;
    }

   private:
    int H_, W_;
    std::vector<T> d_, s_;
    bool built_ = false, has_sum_ = false;
    T& at(long long x, long long y) { return d_[(size_t)x * (W_ + 1) + (size_t)y]; }
    void need_built(const char* f) const {
        if (!built_) imos2d_internal::fail(f, "先に build() を呼ぶ");
    }
    T extreme(bool mx) const {
        need_built(mx ? "Imos2D::max_value" : "Imos2D::min_value");
        if (H_ == 0 || W_ == 0) imos2d_internal::fail(mx ? "Imos2D::max_value" : "Imos2D::min_value", "マスが 0 個");
        T r = d_[0];
        for (int x = 0; x < H_; ++x)
            for (int y = 0; y < W_; ++y) {
                const T& v = d_[(size_t)x * (W_ + 1) + y];
                if (mx ? r < v : v < r) r = v;
            }
        return r;
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace imos2d_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define IM_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[pre_imos2d] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int H = (int)rnd(rng, 0, 12), W = (int)rnd(rng, 0, 12);
        Imos2D<long long> im(H, W);
        std::vector<std::vector<long long>> g(H, std::vector<long long>(W, 0));
        const int ops = (int)rnd(rng, 0, 30);
        for (int t = 0; t < ops; ++t) {
            // はみ出し・空の長方形も出す
            long long x1 = rnd(rng, -3, H + 3), x2 = rnd(rng, -3, H + 3), y1 = rnd(rng, -3, W + 3), y2 = rnd(rng, -3, W + 3);
            const long long v = rnd(rng, -5, 5);
            const int kind = (int)rnd(rng, 0, 2);
            if (kind == 0) {
                if (rnd(rng, 0, 4) && x1 > x2) std::swap(x1, x2);
                if (rnd(rng, 0, 4) && y1 > y2) std::swap(y1, y2);
                im.add(x1, y1, x2, y2, v);
                for (long long x = std::max(0LL, x1); x < std::min<long long>(H, x2); ++x)
                    for (long long y = std::max(0LL, y1); y < std::min<long long>(W, y2); ++y) g[x][y] += v;
            } else if (kind == 1) {
                if (rnd(rng, 0, 4) && x1 > x2) std::swap(x1, x2);
                if (rnd(rng, 0, 4) && y1 > y2) std::swap(y1, y2);
                im.add_closed(x1, y1, x2, y2, v);
                for (long long x = std::max(0LL, x1); x <= std::min<long long>(H - 1, x2); ++x)
                    for (long long y = std::max(0LL, y1); y <= std::min<long long>(W - 1, y2); ++y) g[x][y] += v;
            } else {
                im.add_point(x1, y1, v);
                if (0 <= x1 && x1 < H && 0 <= y1 && y1 < W) g[x1][y1] += v;
            }
        }
        IM_CHECK(!im.is_built());
        im.build();
        IM_CHECK(im.is_built() && im.height() == H && im.width() == W);
        IM_CHECK(im.to_vector() == g);
        long long pos = 0, mx = 0, mn = 0;
        for (int x = 0; x < H; ++x)
            for (int y = 0; y < W; ++y) {
                IM_CHECK(im.get(x, y) == g[x][y]);
                pos += g[x][y] > 0;
                if (x == 0 && y == 0) mx = mn = g[x][y];
                mx = std::max(mx, g[x][y]), mn = std::min(mn, g[x][y]);
                ++cases;
            }
        IM_CHECK(im.count_if([](long long v) { return v > 0; }) == pos);
        if (H > 0 && W > 0) IM_CHECK(im.max_value() == mx && im.min_value() == mn);
        im.build_sum();
        for (int q = 0; q < 30; ++q) {
            const long long x1 = rnd(rng, -2, H + 2), x2 = rnd(rng, -2, H + 2), y1 = rnd(rng, -2, W + 2), y2 = rnd(rng, -2, W + 2);
            long long s = 0, sc = 0;
            for (long long x = std::max(0LL, x1); x < std::min<long long>(H, x2); ++x)
                for (long long y = std::max(0LL, y1); y < std::min<long long>(W, y2); ++y) s += g[x][y];
            for (long long x = std::max(0LL, x1); x <= std::min<long long>(H - 1, x2); ++x)
                for (long long y = std::max(0LL, y1); y <= std::min<long long>(W - 1, y2); ++y) sc += g[x][y];
            IM_CHECK(im.sum(x1, y1, x2, y2) == s);
            IM_CHECK(im.sum_closed(x1, y1, x2, y2) == sc);
            ++cases;
        }
        // clear してやり直せる
        im.clear();
        IM_CHECK(!im.is_built());
        im.add(0, 0, H, W, 7);
        im.build();
        if (H > 0 && W > 0) IM_CHECK(im.get(H - 1, W - 1) == 7 && im.get(0, 0) == 7);
    }
    // 大きい座標（long long）でもはみ出しを捨てられる
    {
        Imos2D<int> im(3, 4);
        im.add(-1000000000000LL, 1, 1000000000000LL, 3, 2);
        im.add_closed(1, -5, 1, 1000000000000000000LL, 1);
        im.build();
        IM_CHECK(im.get(0, 1) == 2 && im.get(2, 2) == 2 && im.get(1, 0) == 1 && im.get(1, 1) == 3 && im.get(0, 3) == 0);
        cases += 5;
    }
    std::cout << "[pre_imos2d] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[pre_imos2d] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int H = TEST_SCALE == 1 ? 500 : 2000, W = H, Q = TEST_SCALE == 1 ? 100000 : 1000000;
    std::vector<int> a(Q), b(Q), c(Q), d(Q);
    for (int i = 0; i < Q; ++i) {
        a[i] = (int)rnd(rng, 0, H), b[i] = (int)rnd(rng, 0, H), c[i] = (int)rnd(rng, 0, W), d[i] = (int)rnd(rng, 0, W);
        if (a[i] > b[i]) std::swap(a[i], b[i]);
        if (c[i] > d[i]) std::swap(c[i], d[i]);
    }
    const auto t0 = std::chrono::steady_clock::now();
    Imos2D<long long> im(H, W);
    for (int i = 0; i < Q; ++i) im.add(a[i], c[i], b[i], d[i], 1);
    im.build();
    im.build_sum();
    long long chk = im.sum(0, 0, H, W) + im.count_if([](long long v) { return v % 2 == 0; });
    std::cout << "[pre_imos2d] benchmark (H = W = " << H << "、add " << Q << " 回 + build + build_sum): " << elapsed_ms(t0)
              << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace imos2d_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace imos2d_demo {

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
#define SHOW(expr, comment) imos2d_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    imos2d_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("2 次元いもす（pre_imos2d） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Imos2D<T> im(H, W)   add / add_closed / add_point → build → get / to_vector / count_if / max_value / min_value\n"
              << "  build_sum → sum / sum_closed      height / width / is_built / clear\n";

    api("【1】 Imos2D<> im(H, W) / add(x1, y1, x2, y2, v) / build()", "— / void / void",
        "[x1, x2) × [y1, y2) に v を足す予約を何回もして、最後に build() で各マスの値にする", "O(HW) / O(1) / O(HW)");
    RUN(Imos2D<> im(4, 5));
    RUN(im.add(0, 0, 2, 3, 1));    // x = 0〜1, y = 0〜2 に 1
    RUN(im.add(1, 1, 4, 5, 10));   // x = 1〜3, y = 1〜4 に 10
    RUN(im.build());
    SHOW(im.get(0, 0), "1 つ目の長方形だけ");
    SHOW(im.get(1, 2), "両方が重なる → 11");
    SHOW(im.get(3, 4), "2 つ目だけ");
    SHOW(im.get(3, 0), "どちらも無い");
    SHOW(im.to_vector(), "全マス（x ごとの行）");

    api("【2】 add_closed(x1, y1, x2, y2, v) / add_point(x, y, v)", "void",
        "閉区間 [x1, x2] × [y1, y2] / 1 マス。はみ出した部分は捨てる（座標は long long でよい）", "O(1)");
    RUN(Imos2D<> g(3, 3));
    RUN(g.add_closed(0, 0, 1, 1, 5));
    RUN(g.add_point(2, 2, 7));
    RUN(g.add(-100, 2, 100, 100, 1));   // y = 2 の列全部（はみ出しは捨てる）
    RUN(g.build());
    SHOW(g.to_vector(), "");

    api("【3】 count_if(f) / max_value() / min_value()", "long long / T / T", "値の条件に合うマスの数 / 最大値 / 最小値", "O(HW)");
    SHOW(im.count_if([](long long v) { return v >= 10; }), "im で 10 以上のマス");
    SHOW(im.max_value(), "最大値");
    SHOW(im.min_value(), "最小値");

    api("【4】 build_sum() / sum(x1, y1, x2, y2) / sum_closed(…)", "void / T / T",
        "build の後に長方形の和の表を作る / [x1, x2) × [y1, y2) の和（閉区間は _closed）", "O(HW) / O(1)");
    RUN(im.build_sum());
    SHOW(im.sum(0, 0, 4, 5), "全体の和");
    SHOW(im.sum(1, 1, 2, 3), "マス (1, 1)・(1, 2) の和 = 11 + 11");
    SHOW(im.sum_closed(3, 4, 3, 4), "1 マス (3, 4)");

    api("【5】 clear()", "void", "全部 0・build 前に戻す（足し直すとき）", "O(HW)");
    RUN(im.clear());
    RUN(im.add(0, 0, 1, 1, 3));
    RUN(im.build());
    SHOW(im.get(0, 0), "");

    api("【6】 使い方の例: 典型90 028（Cluttered Paper）の形", "—",
        "紙 i は [lx, rx) × [ly, ry)。ちょうど k 枚重なっている面積を k = 1〜N で答える（下は自作の小さい例）", "O(N + HW)");
    {
        code("Imos2D<int> im(1000, 1000);");
        code("for (紙 (lx, ly, rx, ry)) im.add(lx, ly, rx, ry, 1);");
        code("im.build();  for (k = 1..N) 答え[k] = im.count_if([&](int v) { return v == k; });");
        const std::vector<std::vector<int>> papers = {{1, 1, 3, 4}, {2, 2, 4, 3}, {0, 0, 5, 1}};
        Imos2D<int> p(6, 6);
        for (const auto& q : papers) p.add(q[0], q[1], q[2], q[3], 1);
        p.build();
        std::vector<long long> ans;
        for (int k = 1; k <= (int)papers.size(); ++k) ans.push_back(p.count_if([&](int v) { return v == k; }));
        note("紙 = {[1,3)×[1,4), [2,4)×[2,3), [0,5)×[0,1)} → k = 1, 2, 3 の面積 = " + to_s(ans));
        note("マスが多い（1000 × 1000 など）ときは Imos2D<int> にするとメモリが半分になる");
    }
}

}  // namespace imos2d_demo

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
    imos2d_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    imos2d_test::benchmark();    // 速度計測
    imos2d_demo::demo();         // 使い方デモ
    return 0;
}
