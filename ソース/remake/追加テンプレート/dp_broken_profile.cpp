// =============================================================================
//  枠 DP（broken profile：グリッドを 1 マスずつ進み、境界の W マス分の状態を bit で持つ DP）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  幅 W が小さい（W ≤ 20 くらい）グリッドで、左上から 1 マスずつ進む。マス (i, j) を見るときの状態 mask の
//  bit k は「k < j なら (i, k)、k ≥ j なら (i − 1, k) の何か（塗ったか・埋まっているか）」を表す（ちょうど境界の W マス）。
//  1 マスごとに 2^W 個の状態を更新するので O(H W 2^W)。行ごとの bit DP（O(H 4^W)）より速い
//    例: ドミノの敷き詰め方の数、隣り合わないようにマスを選ぶ方法の数・重みの最大
//
//  機能一覧（grid は vector<string>。'#' は使えないマス、それ以外は使えるマス）                 計算量
//    broken_profile_dp(H, W, init, zero, trans, combine)   一般形。init[mask] から始め、マスごとに
//        trans(i, j, mask, value, emit) を呼ぶ。emit(新しい mask, 新しい値) で次の状態へ。同じ mask に来た値は
//        combine(今の値, 来た値) でまとめる（数え上げなら足す、最大なら max）。zero は「来ていない」値。最後の表を返す  O(H W 2^W)
//    count_domino_tilings<T>(grid)        使えるマスを 1 × 2 のドミノで隙間なく埋める方法の数                 O(H W 2^W)
//    count_independent_cells<T>(grid)     使えるマスから、上下左右に隣り合わないようにマスを選ぶ方法の数（空も 1 通り）
//    max_independent_cells(w, grid)       同じ条件で選んだマスの重み w[i][j] の和の最大
//
//  注意: W が大きく H が小さいときは転置して W を小さくする（状態は 2^min(H, W)）。W ≤ 25 くらいまで
//  速さの目安: count_domino_tilings（16 × 16）約 0.015 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    profile_test::random_test()  … 全部の置き方・選び方を調べる愚直解との突き合わせ
//    profile_test::benchmark()    … 速度計測
//    profile_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace profile_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
inline void check_grid(const std::vector<std::string>& g, const char* f) {
    const size_t W = g.empty() ? 0 : g[0].size();
    for (const auto& r : g)
        if (r.size() != W) fail(f, "行の長さがそろっていない");
    if (W > 25) fail(f, "W が大きすぎる（転置して W を小さくする）");
}
}  // namespace profile_internal

template <class T, class Trans, class Combine>
std::vector<T> broken_profile_dp(int H, int W, const std::vector<T>& init, T zero, Trans trans, Combine combine) {
    if (W < 0 || W > 25 || H < 0) profile_internal::fail("broken_profile_dp", "0 ≤ W ≤ 25、H ≥ 0 にする");
    const size_t S = (size_t)1 << W;
    if (init.size() != S) profile_internal::fail("broken_profile_dp", "init の長さは 2^W");
    std::vector<T> cur = init, nxt(S);
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            std::fill(nxt.begin(), nxt.end(), zero);
            for (size_t m = 0; m < S; ++m) {
                if (cur[m] == zero) continue;
                trans(i, j, (unsigned)m, cur[m], [&](unsigned nm, const T& v) {
                    if (nm >= S) profile_internal::fail("broken_profile_dp", "emit した mask が 2^W 以上");
                    nxt[nm] = nxt[nm] == zero ? v : combine(nxt[nm], v);
                });
            }
            std::swap(cur, nxt);
        }
    return cur;
}

// ドミノ：bit j = 「(i, j) はもう埋まっている」（上からの縦置き・左からの横置き）
template <class T>
T count_domino_tilings(const std::vector<std::string>& g) {
    profile_internal::check_grid(g, "count_domino_tilings");
    const int H = (int)g.size(), W = H ? (int)g[0].size() : 0;
    std::vector<T> init((size_t)1 << W, T(0));
    init[0] = T(1);
    auto free_cell = [&](int i, int j) { return 0 <= i && i < H && 0 <= j && j < W && g[i][j] != '#'; };
    auto res = broken_profile_dp<T>(
        H, W, init, T(0),
        [&](int i, int j, unsigned m, const T& v, auto emit) {
            const unsigned b = 1u << j;
            if (m & b) {  // もう埋まっている → 下の (i + 1, j) はまだ
                emit(m & ~b, v);
                return;
            }
            if (!free_cell(i, j)) {
                emit(m, v);
                return;
            }
            if (free_cell(i + 1, j)) emit(m | b, v);                                     // 縦置き：(i + 1, j) を埋める
            if (free_cell(i, j + 1) && !(m & (b << 1))) emit((m & ~b) | (b << 1), v);    // 横置き：(i, j + 1) を埋める
        },
        [](const T& a, const T& c) { return a + c; });
    return res[0];
}
// 隣り合わない選び方：bit j = 「境界のそのマスを選んだ」
template <class T>
T count_independent_cells(const std::vector<std::string>& g) {
    profile_internal::check_grid(g, "count_independent_cells");
    const int H = (int)g.size(), W = H ? (int)g[0].size() : 0;
    std::vector<T> init((size_t)1 << W, T(0));
    init[0] = T(1);
    auto res = broken_profile_dp<T>(
        H, W, init, T(0),
        [&](int i, int j, unsigned m, const T& v, auto emit) {
            const unsigned b = 1u << j;
            emit(m & ~b, v);  // 選ばない
            const bool up = m & b, left = j > 0 && (m & (b >> 1));
            if (g[i][j] != '#' && !up && !left) emit(m | b, v);
        },
        [](const T& a, const T& c) { return a + c; });
    T s = T(0);
    for (const T& x : res) s += x;
    return s;
}
inline long long max_independent_cells(const std::vector<std::vector<long long>>& w, const std::vector<std::string>& g) {
    profile_internal::check_grid(g, "max_independent_cells");
    const int H = (int)g.size(), W = H ? (int)g[0].size() : 0;
    if ((int)w.size() != H) profile_internal::fail("max_independent_cells", "w と grid の大きさが違う");
    const long long NG = std::numeric_limits<long long>::min();
    std::vector<long long> init((size_t)1 << W, NG);
    init[0] = 0;
    auto res = broken_profile_dp<long long>(
        H, W, init, NG,
        [&](int i, int j, unsigned m, const long long& v, auto emit) {
            const unsigned b = 1u << j;
            emit(m & ~b, v);
            const bool up = m & b, left = j > 0 && (m & (b >> 1));
            if (g[i][j] != '#' && !up && !left) emit(m | b, v + w[i][j]);
        },
        [](long long a, long long c) { return std::max(a, c); });
    return *std::max_element(res.begin(), res.end());
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace profile_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define BP_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[dp_broken_profile] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static ll brute_domino(std::vector<std::string> g) {
    const int H = (int)g.size(), W = H ? (int)g[0].size() : 0;
    auto rec = [&](auto&& self) -> ll {
        int p = -1;
        for (int k = 0; k < H * W && p < 0; ++k)
            if (g[k / W][k % W] == '.') p = k;
        if (p < 0) return 1;
        const int i = p / W, j = p % W;
        ll c = 0;
        g[i][j] = 'x';
        if (j + 1 < W && g[i][j + 1] == '.') g[i][j + 1] = 'x', c += self(self), g[i][j + 1] = '.';
        if (i + 1 < H && g[i + 1][j] == '.') g[i + 1][j] = 'x', c += self(self), g[i + 1][j] = '.';
        g[i][j] = '.';
        return c;
    };
    return rec(rec);
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int H = (int)rnd(rng, 0, 5), W = (int)rnd(rng, 0, 4);
        std::vector<std::string> g(H, std::string(W, '.'));
        const int wall = (int)rnd(rng, 0, 3);
        for (auto& r : g)
            for (auto& c : r)
                if ((int)rnd(rng, 0, 9) < wall) c = '#';
        std::vector<std::vector<ll>> w(H, std::vector<ll>(W));
        for (auto& r : w)
            for (auto& x : r) x = rnd(rng, -5, 20);
        BP_CHECK(count_domino_tilings<ll>(g) == brute_domino(g));
        ll cnt = 0, best = 0;
        const int C = H * W;
        for (int m = 0; m < (1 << C); ++m) {
            bool ok = true;
            ll s = 0;
            for (int k = 0; k < C && ok; ++k) {
                if (!(m >> k & 1)) continue;
                const int i = k / W, j = k % W;
                if (g[i][j] == '#') ok = false;
                if (j + 1 < W && (m >> (k + 1) & 1)) ok = false;
                if (i + 1 < H && (m >> (k + W) & 1)) ok = false;
                s += w[i][j];
            }
            if (ok) ++cnt, best = std::max(best, s);
        }
        BP_CHECK(count_independent_cells<ll>(g) == cnt);
        BP_CHECK(max_independent_cells(w, g) == best);
        cases += 3;
    }
    // よく知られた値：8 × 8 のドミノの敷き詰めは 12988816 通り
    BP_CHECK(count_domino_tilings<ll>(std::vector<std::string>(8, std::string(8, '.'))) == 12988816);
    std::cout << "[dp_broken_profile] random_test: OK (" << cases + 1 << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[dp_broken_profile] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 10 : 16;
    const auto t0 = std::chrono::steady_clock::now();
    const unsigned long long r = count_domino_tilings<unsigned long long>(std::vector<std::string>(N, std::string(N, '.')));
    std::cout << "[dp_broken_profile] benchmark: count_domino_tilings(" << N << " × " << N << ") " << elapsed_ms(t0) << " ms  (check " << r % 1000
              << ")\n";
}

}  // namespace profile_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace profile_demo {

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
#define SHOW(expr, comment) profile_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    profile_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("枠 DP（dp_broken_profile） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  broken_profile_dp(H, W, init, zero, trans, combine)\n"
              << "  count_domino_tilings<T>(grid) / count_independent_cells<T>(grid) / max_independent_cells(w, grid)\n";

    api("【1】 count_domino_tilings<T>(grid)", "T", "'#' 以外のマスを 1 × 2 のドミノで隙間なく埋める方法の数（mod なら T = mint）", "O(H W 2^W)");
    SHOW(count_domino_tilings<long long>(std::vector<std::string>{"..", ".."}), "2 × 2 は縦 2 本か横 2 本");
    SHOW(count_domino_tilings<long long>(std::vector<std::string>(8, "........")), "8 × 8（チェス盤）");
    SHOW(count_domino_tilings<long long>(std::vector<std::string>{"#...", "....", "...#"}), "角が 2 つ欠けた 3 × 4");

    api("【2】 count_independent_cells<T>(grid) / max_independent_cells(w, grid)", "T / long long",
        "上下左右に隣り合わないようにマスを選ぶ方法の数（空も 1 通り）/ 選んだ重みの和の最大", "O(H W 2^W)");
    SHOW(count_independent_cells<long long>(std::vector<std::string>{"...", "..."}), "2 × 3");
    RUN(const std::vector<std::vector<long long>> w = {{5, 1, 5}, {1, 9, 1}});
    SHOW(max_independent_cells(w, std::vector<std::string>{"...", "..."}), "斜めは隣ではないので 5 + 5 + 9");
    SHOW(max_independent_cells(w, std::vector<std::string>{"...", ".#."}), "9 のマスが使えない → 5 + 5");

    api("【3】 broken_profile_dp(H, W, init, zero, trans, combine)：一般形", "vector<T>（最後の 2^W 個の値）",
        "trans(i, j, mask, value, emit) で次の状態を emit(新しい mask, 値)。bit k は k < j なら (i, k)、k ≥ j なら (i − 1, k)",
        "O(H W 2^W)");
    code("// 例: 2 × 3 の各マスを白か黒に塗る。縦に並ぶ 2 マスがどちらも黒にならない塗り方（bit = 黒）");
    RUN(std::vector<long long> init(8, 0));
    RUN(init[0] = 1);
    RUN(auto res = broken_profile_dp<long long>(2, 3, init, 0, [](int, int j, unsigned m, long long v, auto emit) { emit(m & ~(1u << j), v); if (!(m >> j & 1)) emit(m | (1u << j), v); }, [](long long a, long long b) { return a + b; }));
    SHOW(res, "最後の行の塗り方ごとの数");
    {
        long long s = 0;
        for (long long x : res) s += x;
        note("→ 合計 " + to_s(s) + " 通り（1 列ごとに 3 通り → 3^3）");
    }
}

}  // namespace profile_demo

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
    profile_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    profile_test::benchmark();    // 速度計測
    profile_demo::demo();         // 使い方デモ
    return 0;
}
