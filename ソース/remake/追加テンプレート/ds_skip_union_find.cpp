// =============================================================================
//  区間を飛ばす Union-Find（「i 以上でまだ処理していない最小の位置」を O(α)。各位置を 1 回だけ処理する）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  処理した位置 i を「i + 1 へ飛ぶ」矢印にして、根をたどる（経路圧縮）と、次のまだ処理していない位置が見つかる。
//  区間 [l, r) のまだの位置を全部処理するのが、全体で O((N + Q) α(N))
//    例: 区間を塗る（一度塗ったマスはもう塗らない／後ろの質問から処理すると「最後に塗った色」）、
//        「まだ使っていない一番小さい番号以上の席」、BFS で「まだ訪れていない頂点」を区間でまとめて取る
//
//  機能一覧（位置は 0〜n−1。区間は半開区間 [l, r)）                                              計算量
//    SkipUnionFind sk(n)                 全部まだ                                                     O(n)
//    sk.next(i)                          i 以上でまだの最小の位置（無ければ n）                         ならし O(α(n))
//    sk.prev(i)                          i 以下でまだの最大の位置（無ければ −1）
//    sk.erase(i)                         位置 i を処理済みにする（もう済みなら何もしない）
//    sk.contains(i)                      位置 i がまだか                                              O(1)
//    sk.take_range(l, r, f)              [l, r) のまだの位置を小さい順に f(i) してから処理済みにする     全体 ならし O((n + Q) α(n))
//    sk.count()                          まだの位置の数                                               O(1)
//
//  注意: next(n) は n、prev(−1) は −1 を返す（範囲の端として使ってよい）
//  速さの目安: n = 10^6 で take_range と next を 10^6 回ずつ（ばらばらの区間）約 0.1 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    skipuf_test::random_test()  … std::set を使う愚直解との突き合わせ
//    skipuf_test::benchmark()    … 速度計測
//    skipuf_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace skipuf_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace skipuf_internal

struct SkipUnionFind {
    explicit SkipUnionFind(int n = 0) : n_(n), cnt_(n), nx_(n + 1), pv_(n + 1), alive_(n, 1) {
        if (n < 0) skipuf_internal::fail("SkipUnionFind", "n は 0 以上");
        std::iota(nx_.begin(), nx_.end(), 0);  // nx_[i] = i（まだ）。済みなら i + 1 へ
        std::iota(pv_.begin(), pv_.end(), 0);  // pv_[i + 1] で位置 i を表す（pv_[0] = −1 の番兵）。済みなら i へ
    }
    int next(int i) {
        if (i < 0) i = 0;
        if (i >= n_) return n_;
        return find(nx_, i);
    }
    int prev(int i) {
        if (i < 0) return -1;
        if (i >= n_) i = n_ - 1;
        return find(pv_, i + 1) - 1;
    }
    bool contains(int i) const { return 0 <= i && i < n_ && alive_[i]; }
    void erase(int i) {
        if (i < 0 || i >= n_) skipuf_internal::fail("SkipUnionFind::erase", "位置が範囲外");
        if (!alive_[i]) return;
        alive_[i] = 0, --cnt_;
        nx_[i] = i + 1;
        pv_[i + 1] = i;
    }
    template <class F>
    void take_range(int l, int r, F f) {
        if (l < 0 || r > n_ || l > r) skipuf_internal::fail("SkipUnionFind::take_range", "0 ≤ l ≤ r ≤ n の半開区間 [l, r)");
        for (int i = next(l); i < r; i = next(i + 1)) {
            f(i);
            erase(i);
        }
    }
    int count() const { return cnt_; }
    int size() const { return n_; }

   private:
    int n_, cnt_;
    std::vector<int> nx_, pv_;
    std::vector<char> alive_;
    static int find(std::vector<int>& p, int x) {
        int r = x;
        while (p[r] != r) r = p[r];
        while (p[x] != r) {  // 経路圧縮
            const int y = p[x];
            p[x] = r, x = y;
        }
        return r;
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace skipuf_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define SK_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_skip_union_find] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 40);
        SkipUnionFind sk(n);
        std::set<int> alive;
        for (int i = 0; i < n; ++i) alive.insert(i);
        for (int t = 0; t < 80; ++t) {
            const int k = (int)rnd(rng, 0, 4);
            const int i = (int)rnd(rng, -2, n + 2);
            if (k == 0 && 0 <= i && i < n) {
                sk.erase(i), alive.erase(i);
            } else if (k == 1) {
                int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
                if (l > r) std::swap(l, r);
                std::vector<int> got, exp;
                sk.take_range(l, r, [&](int x) { got.push_back(x); });
                for (auto it2 = alive.lower_bound(l); it2 != alive.end() && *it2 < r;) exp.push_back(*it2), it2 = alive.erase(it2);
                SK_CHECK(got == exp);
            } else {
                auto a = alive.lower_bound(std::max(i, 0));
                SK_CHECK(sk.next(i) == (i >= n ? n : (a == alive.end() ? n : *a)));
                auto b = alive.upper_bound(std::min(i, n - 1));
                SK_CHECK(sk.prev(i) == (i < 0 ? -1 : (b == alive.begin() ? -1 : *std::prev(b))));
                SK_CHECK(sk.contains(i) == (alive.count(i) > 0));
            }
            SK_CHECK(sk.count() == (int)alive.size());
            ++cases;
        }
    }
    std::cout << "[ds_skip_union_find] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_skip_union_find] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 100000 : 1000000;
    const auto t0 = std::chrono::steady_clock::now();
    SkipUnionFind sk(N);
    ll chk = 0;
    for (int i = 0; i < N; ++i) {
        int l = (int)rnd(rng, 0, N), r = (int)rnd(rng, 0, N);
        if (l > r) std::swap(l, r);
        sk.take_range(l, r, [&](int x) { chk += x; });
        chk += sk.next((int)rnd(rng, 0, N - 1));
    }
    std::cout << "[ds_skip_union_find] benchmark: n = " << N << " で take_range と next を " << N << " 回ずつ " << elapsed_ms(t0) << " ms  (check " << chk % 1000
              << ")\n";
}

}  // namespace skipuf_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace skipuf_demo {

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
#define SHOW(expr, comment) skipuf_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    skipuf_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("区間を飛ばす Union-Find（ds_skip_union_find） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  SkipUnionFind sk(n): next / prev / erase / contains / take_range(l, r, f) / count / size\n";

    api("【1】 next(i) / prev(i) / erase(i)", "int / int / void", "i 以上でまだの最小 / i 以下でまだの最大 / 処理済みにする", "ならし O(α(n))");
    RUN(SkipUnionFind sk(8));
    RUN(sk.erase(3); sk.erase(4); sk.erase(6));
    SHOW(sk.next(3), "3, 4 は済み → 5");
    SHOW(sk.next(6), "");
    SHOW(sk.prev(4), "");
    SHOW(sk.count(), "");

    api("【2】 take_range(l, r, f)", "void", "[l, r) のまだの位置を小さい順に f してから処理済みに（各位置は全体で 1 回だけ）", "全体 ならし O((n + Q) α(n))");
    {
        code("// 区間を塗る：後ろの質問から処理して、まだのマスだけ塗る → 各マスは「最後に塗った色」になる");
        const int n = 10;
        const std::vector<std::array<int, 3>> paint = {{0, 6, 1}, {3, 9, 2}, {5, 7, 3}};  // [l, r) を色 c に
        SkipUnionFind s(n);
        std::vector<int> col(n, 0);
        for (int k = (int)paint.size() - 1; k >= 0; --k) s.take_range(paint[k][0], paint[k][1], [&](int i) { col[i] = paint[k][2]; });
        note("[0,6) を 1、[3,9) を 2、[5,7) を 3 の順に塗る → " + to_s(col));
    }
}

}  // namespace skipuf_demo

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
    skipuf_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    skipuf_test::benchmark();    // 速度計測
    skipuf_demo::demo();         // 使い方デモ
    return 0;
}
