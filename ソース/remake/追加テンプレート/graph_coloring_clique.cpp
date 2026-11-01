// =============================================================================
//  彩色数・最大独立集合・最大クリーク（頂点数が小さい一般グラフ）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  彩色数：独立集合の数 I(S) を全部の部分集合で数え、「k 個の独立集合で全体を覆える ⇔ Σ_S (−1)^{N−|S|} I(S)^k > 0」
//          （包除原理）で k を 1 から試す。O(2^N N)。mod を 2 つ使って数えるので、まず間違えない
//  最大独立集合：次数 0・1 の頂点は必ず入れてよい、それ以外は「次数の大きい頂点を入れる / 入れない」で分ける分枝。N ≤ 50 くらい
//    例: ABC187 F（Close Group）… 補グラフの彩色数（= クリークに分ける最小の数）→ min_clique_cover(g)
//        「仲の悪い人どうしを同じ部屋にしない最小の部屋数」= 彩色数、「互いに知り合いの最大の人数」= 最大クリーク
//
//  機能一覧（g は無向グラフの隣接リスト vector<vector<int>>。頂点は 0〜N−1）                          計算量
//    chromatic_number(g)                 隣り合う頂点が違う色になる塗り方の最小の色の数（N ≤ 25 くらい）       O(2^N N)
//    min_clique_cover(g)                 頂点をクリーク（全部つながっている組）に分ける最小の数 = 補グラフの彩色数
//    max_independent_set(g)              どの 2 つも隣り合わない最大の頂点集合（頂点の番号の一覧。N ≤ 64）       O(1.4^N くらい)
//    max_clique(g)                       どの 2 つも隣り合う最大の頂点集合（補グラフの最大独立集合）
//
//  注意: 多重辺・自己ループは無視する（自己ループがあっても色は塗れるものとする）。N = 0 の彩色数は 0
//  速さの目安: chromatic_number（N = 20）約 0.1 秒、max_independent_set（N = 50、辺の確率 0.3）約 0.001 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    color_test::random_test()  … 塗り方・部分集合を全部調べる愚直解との突き合わせ
//    color_test::benchmark()    … 速度計測
//    color_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdint>
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
namespace color_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
inline std::vector<uint64_t> to_masks(const std::vector<std::vector<int>>& g, bool complement, const char* f) {
    const int n = (int)g.size();
    if (n > 64) fail(f, "頂点が多すぎる（64 まで）");
    std::vector<uint64_t> adj(n, 0);
    for (int v = 0; v < n; ++v)
        for (int u : g[v]) {
            if (u < 0 || u >= n) fail(f, "頂点番号が範囲外");
            if (u != v) adj[v] |= 1ULL << u, adj[u] |= 1ULL << v;
        }
    if (complement) {
        const uint64_t all = n == 64 ? ~0ULL : ((1ULL << n) - 1);
        for (int v = 0; v < n; ++v) adj[v] = ~adj[v] & all & ~(1ULL << v);
    }
    return adj;
}
inline int chromatic(const std::vector<uint64_t>& adj) {
    const int n = (int)adj.size();
    if (n == 0) return 0;
    if (n > 27) fail("chromatic_number", "頂点が多すぎる（2^N の表を作る）");
    const size_t full = (size_t)1 << n;
    // I[S] = S に含まれる独立集合の数（空集合も含む）：I[S] = I[S − v] + I[S − N[v] − v]。素数 2 つの mod で数える
    const uint64_t mods[2] = {1000000007ULL, 998244353ULL};
    std::vector<std::vector<uint64_t>> Im(2, std::vector<uint64_t>(full));
    for (int t = 0; t < 2; ++t) {
        auto& J = Im[t];
        J[0] = 1;
        for (size_t S = 1; S < full; ++S) {
            const int v = __builtin_ctzll(S);
            J[S] = (J[S & (S - 1)] + J[(S & (S - 1)) & ~(size_t)adj[v]]) % mods[t];
        }
    }
    std::vector<std::vector<uint64_t>> pw = Im;  // I(S)^k
    for (int k = 1; k <= n; ++k) {
        bool ok = false;
        for (int t = 0; t < 2 && !ok; ++t) {
            uint64_t s = 0;
            for (size_t S = 0; S < full; ++S) {
                const bool neg = ((n - __builtin_popcountll(S)) & 1);
                s = (s + (neg ? mods[t] - pw[t][S] : pw[t][S])) % mods[t];
            }
            ok = s != 0;
        }
        if (ok) return k;
        for (int t = 0; t < 2; ++t)
            for (size_t S = 0; S < full; ++S) pw[t][S] = pw[t][S] * Im[t][S] % mods[t];
    }
    return n;
}
// 最大独立集合（分枝）。cand の中から選ぶ
inline void mis_rec(const std::vector<uint64_t>& adj, uint64_t cand, uint64_t cur, int cnt, int& best, uint64_t& best_set) {
    if (cand == 0) {
        if (cnt > best) best = cnt, best_set = cur;
        return;
    }
    if (cnt + __builtin_popcountll(cand) <= best) return;
    // 次数 0・1 の頂点は入れてよい。次数が最大の頂点を探す
    int pick = -1, maxdeg = -1;
    for (uint64_t c = cand; c; c &= c - 1) {
        const int v = __builtin_ctzll(c);
        const int d = __builtin_popcountll(adj[v] & cand);
        if (d <= 1) {
            mis_rec(adj, cand & ~(adj[v] | (1ULL << v)), cur | (1ULL << v), cnt + 1, best, best_set);
            return;
        }
        if (d > maxdeg) maxdeg = d, pick = v;
    }
    // pick を入れる / 入れない
    mis_rec(adj, cand & ~(adj[pick] | (1ULL << pick)), cur | (1ULL << pick), cnt + 1, best, best_set);
    mis_rec(adj, cand & ~(1ULL << pick), cur, cnt, best, best_set);
}
inline std::vector<int> mis(const std::vector<uint64_t>& adj) {
    const int n = (int)adj.size();
    int best = -1;
    uint64_t bs = 0;
    const uint64_t all = n == 64 ? ~0ULL : ((1ULL << n) - 1);
    mis_rec(adj, all, 0, 0, best, bs);
    std::vector<int> res;
    for (int v = 0; v < n; ++v)
        if (bs >> v & 1) res.push_back(v);
    return res;
}
}  // namespace color_internal

inline int chromatic_number(const std::vector<std::vector<int>>& g) {
    return color_internal::chromatic(color_internal::to_masks(g, false, "chromatic_number"));
}
inline int min_clique_cover(const std::vector<std::vector<int>>& g) {
    return color_internal::chromatic(color_internal::to_masks(g, true, "min_clique_cover"));
}
inline std::vector<int> max_independent_set(const std::vector<std::vector<int>>& g) {
    return color_internal::mis(color_internal::to_masks(g, false, "max_independent_set"));
}
inline std::vector<int> max_clique(const std::vector<std::vector<int>>& g) {
    return color_internal::mis(color_internal::to_masks(g, true, "max_clique"));
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace color_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define CC_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_coloring_clique] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 600 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 9);
        const int p = (int)rnd(rng, 0, 10);
        std::vector<std::vector<int>> g(n);
        std::vector<std::vector<char>> a(n, std::vector<char>(n, 0));
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if ((int)rnd(rng, 0, 9) < p) g[i].push_back(j), a[i][j] = a[j][i] = 1;
        if (n > 0 && rnd(rng, 0, 3) == 0) g[0].push_back(0);  // 自己ループ（無視される）
        // 彩色数：k 色で塗れるかを全部試す
        auto colorable = [&](int k, bool comp) {
            std::vector<int> c(n, 0);
            auto rec = [&](auto&& self, int v) -> bool {
                if (v == n) return true;
                for (int x = 0; x < k; ++x) {
                    bool ok = true;
                    for (int u = 0; u < v && ok; ++u)
                        if ((comp ? (!a[u][v]) : a[u][v]) && c[u] == x) ok = false;
                    if (ok) {
                        c[v] = x;
                        if (self(self, v + 1)) return true;
                    }
                }
                return false;
            };
            return rec(rec, 0);
        };
        int chi = 0, cov = 0;
        while (!colorable(chi, false)) ++chi;
        while (!colorable(cov, true)) ++cov;
        CC_CHECK(chromatic_number(g) == chi);
        CC_CHECK(min_clique_cover(g) == cov);
        // 最大独立集合・最大クリーク
        int mi = 0, mc = 0;
        for (int m = 0; m < (1 << n); ++m) {
            bool ind = true, cl = true;
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j)
                    if ((m >> i & 1) && (m >> j & 1)) ind = ind && !a[i][j], cl = cl && a[i][j];
            if (ind) mi = std::max(mi, __builtin_popcount(m));
            if (cl) mc = std::max(mc, __builtin_popcount(m));
        }
        const auto I = max_independent_set(g), C = max_clique(g);
        CC_CHECK((int)I.size() == mi && (int)C.size() == mc);
        for (size_t x = 0; x < I.size(); ++x)
            for (size_t y = x + 1; y < I.size(); ++y) CC_CHECK(!a[I[x]][I[y]]);
        for (size_t x = 0; x < C.size(); ++x)
            for (size_t y = x + 1; y < C.size(); ++y) CC_CHECK(a[C[x]][C[y]]);
        cases += 4;
    }
    std::cout << "[graph_coloring_clique] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_coloring_clique] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n1 = TEST_SCALE == 1 ? 14 : 20, n2 = TEST_SCALE == 1 ? 30 : 50;
    std::vector<std::vector<int>> g(n1), h(n2);
    for (int i = 0; i < n1; ++i)
        for (int j = i + 1; j < n1; ++j)
            if (rnd(rng, 0, 1)) g[i].push_back(j);
    for (int i = 0; i < n2; ++i)
        for (int j = i + 1; j < n2; ++j)
            if (rnd(rng, 0, 9) < 3) h[i].push_back(j);
    auto t0 = std::chrono::steady_clock::now();
    ll chk = chromatic_number(g);
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    chk += (ll)max_independent_set(h).size();
    const ll t2 = elapsed_ms(t0);
    std::cout << "[graph_coloring_clique] benchmark: chromatic_number(N = " << n1 << ") " << t1 << " ms / max_independent_set(N = " << n2 << ", 辺の確率 0.3) "
              << t2 << " ms  (check " << chk << ")\n";
}

}  // namespace color_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace color_demo {

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
#define SHOW(expr, comment) color_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    color_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("彩色数・最大独立集合・最大クリーク（graph_coloring_clique） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  chromatic_number(g) / min_clique_cover(g) / max_independent_set(g) / max_clique(g)\n";

    RUN(const std::vector<std::vector<int>> c5 = {{1}, {2}, {3}, {4}, {0}});  // 5 頂点の閉路
    RUN(const std::vector<std::vector<int>> k4 = {{1, 2, 3}, {2, 3}, {3}, {}});  // 4 頂点の完全グラフ

    api("【1】 chromatic_number(g)", "int", "隣り合う頂点が違う色になる最小の色の数", "O(2^N N)");
    SHOW(chromatic_number(c5), "奇数の閉路は 3 色");
    SHOW(chromatic_number(k4), "完全グラフ K4 は 4 色");

    api("【2】 min_clique_cover(g)", "int", "頂点をクリークに分ける最小の数（ABC187 F の形）", "O(2^N N)");
    SHOW(min_clique_cover(c5), "辺 2 本 + 頂点 1 つ");
    SHOW(min_clique_cover(k4), "全体が 1 つのクリーク");

    api("【3】 max_independent_set(g) / max_clique(g)", "vector<int>", "どの 2 つも隣り合わない / 隣り合う最大の頂点の組", "O(1.4^N くらい)");
    SHOW(max_independent_set(c5), "");
    SHOW(max_clique(k4), "");
    SHOW(max_clique(c5), "三角形が無いので 2");
}

}  // namespace color_demo

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
    color_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    color_test::benchmark();    // 速度計測
    color_demo::demo();         // 使い方デモ
    return 0;
}
