// =============================================================================
//  最小平均閉路（辺の重みの平均が最小の閉路。Karp の方法 O(NM)、閉路も復元）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  d[k][v] = ちょうど k 本の辺で v に着く歩道の重みの最小（どこから出てもよい）とすると、
//  最小平均 = min_v max_{0≤k<N} (d[N][v] − d[k][v]) / (N − k)（Karp）。分数のまま比べるので誤差がない
//    例: 「1 周あたりの平均コストが最小のループ」、「(利益の和) / (時間の和) が最大の周回」は二分探索 + 負閉路でもよいが、
//        平均（÷ 辺の数）ならこれ。最大平均は重みを −1 倍する
//
//  機能一覧（有向グラフ。辺は {from, to, w}。頂点は 0〜N−1）                                    計算量
//    min_mean_cycle(n, edges)            optional<MeanCycle>。閉路が無ければ nullopt                     O(NM)
//                                          MeanCycle {num, den, cycle}：平均 = num / den（den は閉路の辺の数、約分しない）、
//                                          cycle = 閉路の頂点の列（cycle[0] → cycle[1] → … → cycle[0]）
//    max_mean_cycle(n, edges)            平均が最大の閉路（重みを −1 倍して求め、戻したもの）
//
//  注意: 重みは long long（N × |w| が long long に収まること）。自己ループも長さ 1 の閉路として数える
//  速さの目安: N = 2000、M = 10^4 で 約 0.15 秒（表は (N + 1) × N 個）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mmc_test::random_test()  … 単純な閉路を全部並べる愚直解との突き合わせ
//    mmc_test::benchmark()    … 速度計測
//    mmc_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
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
namespace mmc_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
inline bool frac_less(long long a, long long b, long long c, long long d) {  // a/b < c/d（b, d > 0）
    return (__int128)a * d < (__int128)c * b;
}
}  // namespace mmc_internal

struct MeanCycle {
    long long num = 0, den = 1;
    std::vector<int> cycle;
};

inline std::optional<MeanCycle> min_mean_cycle(int n, const std::vector<std::array<long long, 3>>& edges) {
    using mmc_internal::frac_less;
    if (n <= 0) return std::nullopt;
    for (const auto& e : edges)
        if (e[0] < 0 || e[0] >= n || e[1] < 0 || e[1] >= n) mmc_internal::fail("min_mean_cycle", "頂点番号が範囲外");
    const long long INF = std::numeric_limits<long long>::max() / 4;
    // d[k][v]、どの辺から来たか pe[k][v]
    std::vector<std::vector<long long>> d(n + 1, std::vector<long long>(n, INF));
    std::vector<std::vector<int>> pe(n + 1, std::vector<int>(n, -1));
    for (int v = 0; v < n; ++v) d[0][v] = 0;
    for (int k = 0; k < n; ++k)
        for (int i = 0; i < (int)edges.size(); ++i) {
            const int a = (int)edges[i][0], b = (int)edges[i][1];
            if (d[k][a] >= INF) continue;
            const long long nd = d[k][a] + edges[i][2];
            if (nd < d[k + 1][b]) d[k + 1][b] = nd, pe[k + 1][b] = i;
        }
    // Karp の式
    bool found = false;
    long long bn = 0, bd = 1;
    int bv = -1;
    for (int v = 0; v < n; ++v) {
        if (d[n][v] >= INF) continue;
        long long wn = 0, wd = 1;  // max_k (d[n][v] − d[k][v]) / (n − k)
        bool any = false;
        for (int k = 0; k < n; ++k) {
            if (d[k][v] >= INF) continue;
            const long long a = d[n][v] - d[k][v], b = n - k;
            if (!any || frac_less(wn, wd, a, b)) wn = a, wd = b, any = true;
        }
        if (any && (!found || frac_less(wn, wd, bn, bd))) bn = wn, bd = wd, bv = v, found = true;
    }
    if (!found) return std::nullopt;
    // d[n][bv] の歩道を戻り、出てきた閉路のうち平均が最小のもの（= 最小平均）を取る
    std::vector<int> walk{bv};  // 逆順の頂点
    std::vector<int> wedge;
    for (int k = n, v = bv; k > 0; --k) {
        const int i = pe[k][v];
        wedge.push_back(i);
        v = (int)edges[i][0];
        walk.push_back(v);
    }
    std::reverse(walk.begin(), walk.end());
    std::reverse(wedge.begin(), wedge.end());  // wedge[j] は walk[j] → walk[j+1]
    MeanCycle best;
    bool have = false;
    std::vector<int> pos(n, -1), stk_v, stk_e;
    std::vector<long long> pre{0};
    for (int j = 0; j <= n; ++j) {
        const int v = walk[j];
        if (pos[v] >= 0) {  // 閉路：stk_v[pos[v]] … v
            const int s = pos[v];
            const long long w = pre.back() - pre[s];
            const long long len = (long long)stk_v.size() - s;
            if (!have || frac_less(w, len, best.num, best.den)) {
                best.num = w, best.den = len, best.cycle.assign(stk_v.begin() + s, stk_v.end());
                have = true;
            }
            for (size_t t = s; t < stk_v.size(); ++t) pos[stk_v[t]] = -1;
            stk_v.resize(s), pre.resize(s + 1);
        }
        pos[v] = (int)stk_v.size();
        stk_v.push_back(v);
        if (j < n) pre.push_back(pre.back() + edges[wedge[j]][2]);
    }
    if (!have || frac_less(bn, bd, best.num, best.den) || frac_less(best.num, best.den, bn, bd))
        mmc_internal::fail("min_mean_cycle", "閉路の復元に失敗（内部エラー）");
    return best;
}
inline std::optional<MeanCycle> max_mean_cycle(int n, std::vector<std::array<long long, 3>> edges) {
    for (auto& e : edges) e[2] = -e[2];
    auto r = min_mean_cycle(n, edges);
    if (r) r->num = -r->num;
    return r;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace mmc_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define MM_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_min_mean_cycle] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 6), m = (int)rnd(rng, 0, 12);
        std::vector<std::array<ll, 3>> es(m);
        for (auto& e : es) e = {rnd(rng, 0, n - 1), rnd(rng, 0, n - 1), rnd(rng, -10, 10)};
        // 愚直：単純な閉路を全部（多重辺は重みの小さい方 / 大きい方）
        bool found = false;
        ll bn = 0, bd = 1, xn = 0, xd = 1;
        std::vector<char> used(n, 0);
        auto dfs = [&](auto&& self, int start, int v, ll w, ll len) -> void {
            for (const auto& e : es) {
                if (e[0] != v) continue;
                const int u = (int)e[1];
                if (u == start) {
                    const ll W = w + e[2], L = len + 1;
                    if (!found || mmc_internal::frac_less(W, L, bn, bd)) bn = W, bd = L;
                    if (!found || mmc_internal::frac_less(xn, xd, W, L)) xn = W, xd = L;
                    found = true;
                } else if (u > start && !used[u]) {
                    used[u] = 1;
                    self(self, start, u, w + e[2], len + 1);
                    used[u] = 0;
                }
            }
        };
        for (int s = 0; s < n; ++s) used[s] = 1, dfs(dfs, s, s, 0, 0), used[s] = 0;
        const auto r = min_mean_cycle(n, es), R = max_mean_cycle(n, es);
        MM_CHECK(r.has_value() == found && R.has_value() == found);
        if (found) {
            MM_CHECK((__int128)r->num * bd == (__int128)bn * r->den);
            MM_CHECK((__int128)R->num * xd == (__int128)xn * R->den);
            // 復元した閉路が本当にあり、その重みの和と長さ
            const auto& c = r->cycle;
            MM_CHECK((ll)c.size() == r->den);
            ll w = 0;
            for (size_t i = 0; i < c.size(); ++i) {
                const int a = c[i], b = c[(i + 1) % c.size()];
                ll mn = std::numeric_limits<ll>::max();
                for (const auto& e : es)
                    if (e[0] == a && e[1] == b) mn = std::min(mn, e[2]);
                MM_CHECK(mn != std::numeric_limits<ll>::max());
                w += mn;
            }
            MM_CHECK(w <= r->num);  // 同じ頂点の間に安い辺があればもっと小さくなるが、平均が最小なので等しいはず
            MM_CHECK(w == r->num);
        }
        ++cases;
    }
    std::cout << "[graph_min_mean_cycle] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_min_mean_cycle] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 500 : 2000, M = 5 * N;
    std::vector<std::array<ll, 3>> es(M);
    for (auto& e : es) e = {rnd(rng, 0, N - 1), rnd(rng, 0, N - 1), rnd(rng, -1000000, 1000000)};
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = min_mean_cycle(N, es);
    std::cout << "[graph_min_mean_cycle] benchmark (N = " << N << ", M = " << M << "): " << elapsed_ms(t0) << " ms  (check " << (r ? r->num % 1000 : -1) << ")\n";
}

}  // namespace mmc_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace mmc_demo {

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
#define SHOW(expr, comment) mmc_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    mmc_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("最小平均閉路（graph_min_mean_cycle） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  min_mean_cycle(n, edges{from, to, w}) / max_mean_cycle(n, edges) → optional<MeanCycle {num, den, cycle}>\n";

    api("【1】 min_mean_cycle(n, edges)", "optional<MeanCycle>", "辺の重みの平均 num / den が最小の閉路（cycle は頂点の列）", "O(NM)");
    RUN(const std::vector<std::array<long long, 3>> es = {{0, 1, 2}, {1, 0, 4}, {1, 2, 1}, {2, 3, 1}, {3, 1, 2}, {3, 3, 5}});
    RUN(const auto r = min_mean_cycle(4, es));
    SHOW(r->num, "");
    SHOW(r->den, "1 → 2 → 3 → 1 は (1 + 1 + 2) / 3");
    SHOW(r->cycle, "");

    api("【2】 max_mean_cycle(n, edges)", "optional<MeanCycle>", "平均が最大の閉路", "O(NM)");
    RUN(const auto R = max_mean_cycle(4, es));
    SHOW(R->num, "自己ループ 3 → 3 の 5 / 1");
    SHOW(R->den, "");
    SHOW(min_mean_cycle(3, {{0, 1, 1}, {1, 2, 1}}).has_value(), "閉路が無い");
}

}  // namespace mmc_demo

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
    mmc_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mmc_test::benchmark();    // 速度計測
    mmc_demo::demo();         // 使い方デモ
    return 0;
}
