// =============================================================================
//  Slope trick（下に凸な区分線形関数を 2 つのヒープで持ち、|x − a| を足す・累積 min・平行移動する）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  f(x) = min_f + Σ_{l ∈ L} max(0, l − x) + Σ_{r ∈ R} max(0, x − r) の形で、傾きが変わる点を
//  左側 L（最大ヒープ）と右側 R（最小ヒープ）に入れておく。dp[i][x] の x についての関数がこの形のまま更新できるとき、
//  O(N log N) で DP ができる
//    例: 「数列を広義単調増加にする Σ|a_i − b_i| の最小」、ABC217 H（Snuketoon）など「各時刻に位置の関数を足し、1 秒で ±1 動ける」
//
//  機能一覧（値・座標は long long）                                                           計算量
//    SlopeTrick st                       f(x) = 0
//    st.add_all(a)                       f(x) += a                                                 O(1)
//    st.add_x_minus_a(a)                 f(x) += max(0, x − a)（a より右で増える）                    O(log N)
//    st.add_a_minus_x(a)                 f(x) += max(0, a − x)（a より左で増える）                    O(log N)
//    st.add_abs(a)                       f(x) += |x − a|                                           O(log N)
//    st.clear_right()                    f(x) ← min_{y ≤ x} f(y)（累積 min。右側を平らに）            O(1)（ならし）
//    st.clear_left()                     f(x) ← min_{y ≥ x} f(y)（左側を平らに）
//    st.shift(a)                         f(x) ← f(x − a)（右へ a ずらす）                            O(1)
//    st.slide_min(a, b)                  f(x) ← min_{x − b ≤ y ≤ x − a} f(y)（a ≤ b。1 歩で a〜b 動ける）  O(1)
//    st.get_min() / st.argmin_range()    最小値 / 最小になる x の範囲 {lo, hi}（端が無いときは ∓INF）    O(1)
//    st.eval(x)                          f(x)（ヒープを全部見る）                                     O(N)
//
//  注意: 値はあふれないように（min_f は long long）。INF = 4×10^18 を「端が無い」の意味で使う
//  速さの目安: N = 2×10^5 で add_abs と clear_right を N 回ずつ 約 0.02 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    slope_test::random_test()  … 整数の点ごとに関数の値を持つ愚直解との突き合わせ
//    slope_test::benchmark()    … 速度計測
//    slope_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
struct SlopeTrick {
    static constexpr long long INF = 4000000000000000000LL;
    void add_all(long long a) { min_f_ += a; }
    void add_x_minus_a(long long a) {  // max(0, x − a)
        const long long l = top_L();
        if (l > a) {  // 最小の範囲より左に a がある → 最小値が増える
            min_f_ += l - a;
            push_L(a);
            push_R(pop_L());
        } else {
            push_R(a);
        }
    }
    void add_a_minus_x(long long a) {  // max(0, a − x)
        const long long r = top_R();
        if (r < a) {
            min_f_ += a - r;
            push_R(a);
            push_L(pop_R());
        } else {
            push_L(a);
        }
    }
    void add_abs(long long a) {
        add_x_minus_a(a);
        add_a_minus_x(a);
    }
    void clear_right() {
        while (!R_.empty()) R_.pop();
    }
    void clear_left() {
        while (!L_.empty()) L_.pop();
    }
    void shift(long long a) { add_L_ += a, add_R_ += a; }
    void slide_min(long long a, long long b) {  // f(x) ← min_{x−b ≤ y ≤ x−a} f(y)
        if (a > b) std::swap(a, b);
        add_L_ += a, add_R_ += b;
    }
    long long get_min() const { return min_f_; }
    std::pair<long long, long long> argmin_range() const { return {top_L(), top_R()}; }
    long long eval(long long x) const {
        long long s = min_f_;
        auto l = L_;
        while (!l.empty()) s += std::max(0LL, l.top() + add_L_ - x), l.pop();
        auto r = R_;
        while (!r.empty()) s += std::max(0LL, x - (r.top() + add_R_)), r.pop();
        return s;
    }
    size_t size() const { return L_.size() + R_.size(); }

   private:
    long long min_f_ = 0, add_L_ = 0, add_R_ = 0;
    std::priority_queue<long long> L_;                                              // 中身 + add_L_ が本当の値
    std::priority_queue<long long, std::vector<long long>, std::greater<long long>> R_;  // 中身 + add_R_
    long long top_L() const { return L_.empty() ? -INF : L_.top() + add_L_; }
    long long top_R() const { return R_.empty() ? INF : R_.top() + add_R_; }
    void push_L(long long a) { L_.push(a - add_L_); }
    void push_R(long long a) { R_.push(a - add_R_); }
    long long pop_L() {
        const long long v = top_L();
        if (!L_.empty()) L_.pop();
        return v;
    }
    long long pop_R() {
        const long long v = top_R();
        if (!R_.empty()) R_.pop();
        return v;
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace slope_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define SL_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_slope_trick] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    const int D = 1000, W = 300;  // 関数は [−D, D] の整数で持ち、[−W, W] で比べる（折れ目は ±200 くらいに収まる）
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        SlopeTrick st;
        std::vector<ll> f(2 * D + 1, 0);
        auto at = [&](int x) -> ll& { return f[x + D]; };
        const int ops = (int)rnd(rng, 1, 20);
        for (int t = 0; t < ops; ++t) {
            const int k = (int)rnd(rng, 0, 8);
            const ll a = rnd(rng, -10, 10);
            std::vector<ll> g = f;
            if (k == 0) {
                st.add_all(a);
                for (auto& v : f) v += a;
            } else if (k == 1) {
                st.add_x_minus_a(a);
                for (int x = -D; x <= D; ++x) at(x) += std::max(0LL, x - a);
            } else if (k == 2) {
                st.add_a_minus_x(a);
                for (int x = -D; x <= D; ++x) at(x) += std::max(0LL, a - x);
            } else if (k == 3) {
                st.add_abs(a);
                for (int x = -D; x <= D; ++x) at(x) += std::abs(x - a);
            } else if (k == 4) {
                st.clear_right();
                for (int x = -D + 1; x <= D; ++x) at(x) = std::min(at(x), at(x - 1));
            } else if (k == 5) {
                st.clear_left();
                for (int x = D - 1; x >= -D; --x) at(x) = std::min(at(x), at(x + 1));
            } else if (k == 6) {
                st.shift(a);
                for (int x = -D; x <= D; ++x) {
                    const int y = (int)(x - a);
                    at(x) = (y >= -D && y <= D) ? g[y + D] : (ll)4e18;
                }
            } else {
                const ll b = a + rnd(rng, 0, 5);
                st.slide_min(a, b);
                for (int x = -D; x <= D; ++x) {
                    ll m = (ll)4e18;
                    for (ll y = x - b; y <= x - a; ++y)
                        if (y >= -D && y <= D) m = std::min(m, g[y + D]);
                    at(x) = m;
                }
            }
            // 比べる
            ll mn = (ll)4e18;
            for (int x = -W; x <= W; ++x) mn = std::min(mn, at(x));
            SL_CHECK(st.get_min() == mn);
            const auto [lo, hi] = st.argmin_range();
            for (int x = -W; x <= W; x += 7) SL_CHECK(st.eval(x) == at(x));
            for (int x = -W; x <= W; ++x) SL_CHECK((at(x) == mn) == (lo <= x && x <= hi));
            cases += 3;
        }
    }
    std::cout << "[ds_slope_trick] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_slope_trick] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    const auto t0 = std::chrono::steady_clock::now();
    SlopeTrick st;
    for (int i = 0; i < N; ++i) st.add_abs(rnd(rng, -1000000000, 1000000000)), st.clear_right();
    std::cout << "[ds_slope_trick] benchmark: add_abs と clear_right を " << N << " 回ずつ " << elapsed_ms(t0) << " ms  (check " << st.get_min() % 1000 << ")\n";
}

}  // namespace slope_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace slope_demo {

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
#define SHOW(expr, comment) slope_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    slope_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("Slope trick（ds_slope_trick） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  SlopeTrick st: add_all / add_x_minus_a / add_a_minus_x / add_abs / clear_right / clear_left / shift / slide_min\n"
              << "                 get_min / argmin_range / eval / size\n";

    api("【1】 add_abs(a) / get_min() / argmin_range()", "void / long long / pair", "f(x) += |x − a| / 最小値 / 最小になる x の範囲", "O(log N) / O(1) / O(1)");
    RUN(SlopeTrick st);
    RUN(st.add_abs(1); st.add_abs(5); st.add_abs(3));
    SHOW(st.get_min(), "x = 3 で 2 + 0 + 2");
    SHOW(st.argmin_range(), "");
    RUN(st.add_abs(10));
    SHOW(st.argmin_range(), "4 点の中央値の範囲");
    SHOW(st.eval(0), "f(0) = 1 + 5 + 3 + 10");

    api("【2】 clear_right()：数列を広義単調増加にする Σ|a_i − b_i| の最小", "—",
        "dp_i(x) = min_{y ≤ x} dp_{i−1}(y) + |x − a_i|。各 i で clear_right() してから add_abs(a_i)", "O(N log N)");
    {
        const std::vector<long long> a = {3, 1, 4, 1, 5, 9, 2, 6};
        SlopeTrick t;
        for (long long x : a) t.clear_right(), t.add_abs(x);
        note("a = " + to_s(a) + " → 最小 " + to_s(t.get_min()));
    }

    api("【3】 slide_min(a, b) / shift(a)", "void", "1 歩で a〜b 動ける（f(x) ← min_{x−b ≤ y ≤ x−a} f(y)）/ 右へ a ずらす", "O(1)");
    RUN(SlopeTrick s);
    RUN(s.add_abs(0));        // 最初は x = 0 にいたい
    RUN(s.slide_min(-1, 1));  // 1 秒で ±1 動ける
    SHOW(s.argmin_range(), "どこでも 0 のまま行ける範囲が [−1, 1] に広がる");
    RUN(s.add_abs(5));        // 次の時刻は 5 にいたい
    SHOW(s.get_min(), "1 まで動いて |1 − 5| = 4");
}

}  // namespace slope_demo

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
    slope_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    slope_test::benchmark();    // 速度計測
    slope_demo::demo();         // 使い方デモ
    return 0;
}
