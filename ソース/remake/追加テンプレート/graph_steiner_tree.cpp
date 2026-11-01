// =============================================================================
//  最小シュタイナー木（指定した K 個の頂点をすべてつなぐ最小コストの木。部分集合の DP + ダイクストラ）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  dp[S][v] = 「ターミナルの集合 S と頂点 v をつなぐ木」の最小コスト。
//    ① dp[S][v] = min_{T ⊂ S} dp[T][v] + dp[S∖T][v]（v で 2 つの木をくっつける）
//    ② dp[S][*] を初期値にしてダイクストラ（木を辺に沿って伸ばす）
//  の順に S を大きくしていく。O(3^K N + 2^K (N + M) log N)
//    例: ABC364 G（Last Major City）… K − 1 個の都市を固定し、各都市 v について dp[全部][v] を答える → steiner_table
//
//  機能一覧（g は重み付き無向グラフ vector<vector<pair<int, T>>>（両向きに入れる）。重みは 0 以上）       計算量
//    steiner_tree(g, terminals)          terminals をすべてつなぐ木の最小コスト（つなげないなら steiner_inf<T>()）   O(3^K N + 2^K M log N)
//    steiner_table(g, terminals)         dp[S][v]（S は terminals の部分集合の bit。S = 0 は 0 ではなく「v だけ」= 0）
//    steiner_inf<T>()                    「つなげない」を表す大きい値（T の最大 / 4）
//
//  注意: K は 15 くらいまで（3^K）。グリッドは頂点を i × W + j にしてグラフに直す。頂点に重みがあるときは
//        「頂点に入るときに重みを足す」有向の辺にする（そのとき出発点の重みは別に足す）
//  速さの目安: N = 10^4、M = 4×10^4、K = 8 で 約 1.2 秒（ダイクストラ 2^K 回が重い）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    steiner_test::random_test()  … 辺の部分集合を全部調べる愚直解との突き合わせ
//    steiner_test::benchmark()    … 速度計測
//    steiner_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace steiner_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace steiner_internal

template <class T>
constexpr T steiner_inf() {
    return std::numeric_limits<T>::max() / 4;
}

template <class T>
std::vector<std::vector<T>> steiner_table(const std::vector<std::vector<std::pair<int, T>>>& g, const std::vector<int>& terminals) {
    const int n = (int)g.size(), K = (int)terminals.size();
    if (K > 20) steiner_internal::fail("steiner_table", "ターミナルが多すぎる（3^K の DP）");
    const T INF = steiner_inf<T>();
    const int full = 1 << K;
    std::vector<std::vector<T>> dp(full, std::vector<T>(n, INF));
    for (int v = 0; v < n; ++v) dp[0][v] = T(0);
    for (int i = 0; i < K; ++i) {
        if (terminals[i] < 0 || terminals[i] >= n) steiner_internal::fail("steiner_table", "ターミナルの番号が範囲外");
        dp[1 << i][terminals[i]] = T(0);
    }
    using P = std::pair<T, int>;
    for (int S = 1; S < full; ++S) {
        auto& cur = dp[S];
        // ① 2 つに分けてくっつける（T は S の真部分集合で、一番小さいビットを含むものだけ見れば十分）
        const int low = S & -S;
        for (int sub = (S - 1) & S; sub > 0; sub = (sub - 1) & S) {
            if (!(sub & low)) continue;
            const auto& a = dp[sub];
            const auto& b = dp[S ^ sub];
            for (int v = 0; v < n; ++v)
                if (a[v] < INF && b[v] < INF && a[v] + b[v] < cur[v]) cur[v] = a[v] + b[v];
        }
        // ② ダイクストラで伸ばす
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        for (int v = 0; v < n; ++v)
            if (cur[v] < INF) pq.push({cur[v], v});
        while (!pq.empty()) {
            auto [d, v] = pq.top();
            pq.pop();
            if (d != cur[v]) continue;
            for (auto [u, w] : g[v]) {
                if (w < T(0)) steiner_internal::fail("steiner_table", "重みは 0 以上");
                if (d + w < cur[u]) cur[u] = d + w, pq.push({cur[u], u});
            }
        }
    }
    return dp;
}
template <class T>
T steiner_tree(const std::vector<std::vector<std::pair<int, T>>>& g, const std::vector<int>& terminals) {
    if (terminals.empty()) return T(0);
    const auto dp = steiner_table(g, terminals);
    return *std::min_element(dp.back().begin(), dp.back().end());
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace steiner_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define SN_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_steiner_tree] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
using G = std::vector<std::vector<std::pair<int, ll>>>;

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    const ll INF = steiner_inf<ll>();
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 7);
        std::vector<std::array<ll, 3>> es;
        const int m = (int)rnd(rng, 0, 9);
        for (int i = 0; i < m; ++i) es.push_back({rnd(rng, 0, n - 1), rnd(rng, 0, n - 1), rnd(rng, 0, 9)});
        G g(n);
        for (auto [a, b, w] : es) g[a].push_back({(int)b, w}), g[b].push_back({(int)a, w});
        std::vector<int> term;
        for (int v = 0; v < n; ++v)
            if (rnd(rng, 0, 2) == 0) term.push_back(v);
        if (rnd(rng, 0, 4) == 0 && n > 0) term.push_back(term.empty() ? 0 : term[0]);  // 同じ頂点が 2 回
        const int K = (int)term.size();
        const auto dp = steiner_table(g, term);
        // 愚直：辺の部分集合ごとに、各頂点と同じ成分にある端末の集合 mask を求め、bv[mask][v] にコストの最小を入れる。
        // 最後に「S ⊆ mask なら S もつながる」ので上位集合の最小を下ろす
        std::vector<std::vector<ll>> bv((size_t)1 << K, std::vector<ll>(n, INF));
        std::vector<int> par(n), tm(n);
        auto root = [&](int x) {
            while (par[x] != x) x = par[x] = par[par[x]];
            return x;
        };
        for (int sm = 0; sm < (1 << m); ++sm) {
            std::iota(par.begin(), par.end(), 0);
            ll cost = 0;
            for (int i = 0; i < m; ++i)
                if (sm >> i & 1) par[root((int)es[i][0])] = root((int)es[i][1]), cost += es[i][2];
            std::fill(tm.begin(), tm.end(), 0);
            for (int i = 0; i < K; ++i) tm[root(term[i])] |= 1 << i;
            for (int v = 0; v < n; ++v) {
                ll& x = bv[tm[root(v)]][v];
                x = std::min(x, cost);
            }
        }
        for (int i = 0; i < K; ++i)
            for (int S = 0; S < (1 << K); ++S)
                if (!(S >> i & 1))
                    for (int v = 0; v < n; ++v) bv[S][v] = std::min(bv[S][v], bv[S | 1 << i][v]);
        for (int S = 0; S < (1 << K); ++S)
            for (int v = 0; v < n; ++v) SN_CHECK(dp[S][v] == bv[S][v]);
        ll mn = INF;
        for (int v = 0; v < n; ++v) mn = std::min(mn, bv[(1 << K) - 1][v]);
        SN_CHECK(steiner_tree(g, term) == (K == 0 ? 0 : mn));
        cases += (1 << K) * n + 1;
    }
    std::cout << "[graph_steiner_tree] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_steiner_tree] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 200 : 10000, M = 4 * N, K = 8;  // デバッグビルドは priority_queue の検査が重いので小さく
    G g(N);
    for (int i = 1; i < N; ++i) {
        const int p = (int)rnd(rng, 0, i - 1);
        const ll w = rnd(rng, 1, 1000000000);
        g[i].push_back({p, w}), g[p].push_back({i, w});
    }
    for (int i = N - 1; i < M; ++i) {
        const int a = (int)rnd(rng, 0, N - 1), b = (int)rnd(rng, 0, N - 1);
        const ll w = rnd(rng, 1, 1000000000);
        g[a].push_back({b, w}), g[b].push_back({a, w});
    }
    std::vector<int> term(K);
    for (auto& t : term) t = (int)rnd(rng, 0, N - 1);
    const auto t0 = std::chrono::steady_clock::now();
    const ll r = steiner_tree(g, term);
    std::cout << "[graph_steiner_tree] benchmark (N = " << N << ", M = " << M << ", K = " << K << "): " << elapsed_ms(t0) << " ms  (check " << r % 1000 << ")\n";
}

}  // namespace steiner_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace steiner_demo {

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
#define SHOW(expr, comment) steiner_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    steiner_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("最小シュタイナー木（graph_steiner_tree） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  steiner_tree(g, terminals) / steiner_table(g, terminals) / steiner_inf<T>()\n";

    RUN(std::vector<std::vector<std::pair<int, long long>>> g(6));
    RUN(auto add = [&](int a, int b, long long w) { g[a].push_back({b, w}); g[b].push_back({a, w}); });
    RUN(add(0, 4, 3); add(1, 4, 3); add(2, 4, 3); add(0, 1, 5); add(1, 2, 5); add(2, 3, 1); add(4, 5, 1));

    api("【1】 steiner_tree(g, terminals)", "T", "terminals をすべてつなぐ木の最小コスト", "O(3^K N + 2^K M log N)");
    SHOW(steiner_tree(g, std::vector<int>{0, 1, 2}), "真ん中の 4 を経由する 3 + 3 + 3（直接つなぐと 5 + 5）");
    SHOW(steiner_tree(g, std::vector<int>{0, 3}), "0 → 4 → 2 → 3 = 3 + 3 + 1");
    SHOW(steiner_tree(g, std::vector<int>{1}), "1 つだけなら 0");

    api("【2】 steiner_table(g, terminals)", "vector<vector<T>>", "dp[S][v] = ターミナルの集合 S（bit）と頂点 v をつなぐ最小コスト", "同上");
    RUN(const auto dp = steiner_table(g, std::vector<int>{0, 1}));
    SHOW(dp[3], "0 と 1 と、各頂点 v をつなぐ最小（ABC364 G の形）");
}

}  // namespace steiner_demo

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
    steiner_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    steiner_test::benchmark();    // 速度計測
    steiner_demo::demo();         // 使い方デモ
    return 0;
}
