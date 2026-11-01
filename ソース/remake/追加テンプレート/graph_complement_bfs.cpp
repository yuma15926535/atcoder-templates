// =============================================================================
//  補グラフの BFS・連結成分（「辺が無い組」が与えられ、補グラフの辺が多すぎるとき O(N + M)）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  補グラフ（元のグラフで辺の無い 2 頂点の間に辺がある）は辺が N^2 本近くあるので、そのまま作れない。
//  まだ訪れていない頂点の一覧を持ち、頂点 v から進むとき「一覧のうち、元のグラフで v と隣り合っていないもの」を全部取り出す。
//  残った頂点は元の辺 1 本ずつに対応するので、全体 O(N + M)
//    例: 「M 組の仲の悪い人以外は 1 歩で行ける」ときの最短距離・グループの数、補グラフが連結か
//
//  機能一覧（g は元のグラフの隣接リスト vector<vector<int>>（無向。両向きに入れる））                 計算量
//    complement_bfs(g, s)                補グラフでの s からの最短距離（辺の数。届かないなら −1）          O(N + M)
//    complement_components(g)            補グラフの連結成分の番号（0, 1, …）。成分の数は最大の番号 + 1     O(N + M)
//    complement_component_count(g)       補グラフの連結成分の数
//
//  注意: 自己ループ・多重辺があってもよい（補グラフに自己ループは無いものとする）
//  速さの目安: N = 2×10^5、M = 2×10^5 で complement_bfs と complement_component_count 合わせて 約 0.01 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    compbfs_test::random_test()  … 補グラフを作って BFS する愚直解との突き合わせ
//    compbfs_test::benchmark()    … 速度計測
//    compbfs_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
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
namespace compbfs_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 補グラフで s から広げる。dist[v] を書き、訪れた頂点を返す。rest = まだの頂点の一覧（書き換える）
inline void spread(const std::vector<std::vector<int>>& g, int s, std::vector<int>& rest, std::vector<int>& dist, std::vector<int>& mark, int& stamp,
                   int label, bool use_label) {
    std::vector<int> q{s};
    dist[s] = use_label ? label : 0;
    for (size_t k = 0; k < q.size(); ++k) {
        const int v = q[k];
        ++stamp;
        for (int u : g[v]) mark[u] = stamp;  // 元のグラフで v と隣り合う
        std::vector<int> keep;
        for (int u : rest) {
            if (mark[u] == stamp) keep.push_back(u);
            else dist[u] = use_label ? label : dist[v] + 1, q.push_back(u);
        }
        rest.swap(keep);
    }
}
}  // namespace compbfs_internal

inline std::vector<int> complement_bfs(const std::vector<std::vector<int>>& g, int s) {
    const int n = (int)g.size();
    if (s < 0 || s >= n) compbfs_internal::fail("complement_bfs", "s が範囲外");
    std::vector<int> dist(n, -1), mark(n, 0), rest;
    for (int v = 0; v < n; ++v)
        if (v != s) rest.push_back(v);
    int stamp = 0;
    compbfs_internal::spread(g, s, rest, dist, mark, stamp, 0, false);
    return dist;
}
inline std::vector<int> complement_components(const std::vector<std::vector<int>>& g) {
    const int n = (int)g.size();
    std::vector<int> comp(n, -1), mark(n, 0), rest;
    for (int v = 0; v < n; ++v) rest.push_back(v);
    int stamp = 0, c = 0;
    while (!rest.empty()) {
        const int s = rest.back();
        rest.pop_back();
        compbfs_internal::spread(g, s, rest, comp, mark, stamp, c++, true);
    }
    return comp;
}
inline int complement_component_count(const std::vector<std::vector<int>>& g) {
    const auto c = complement_components(g);
    return c.empty() ? 0 : *std::max_element(c.begin(), c.end()) + 1;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace compbfs_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define CB_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_complement_bfs] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 1000 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 12), p = (int)rnd(rng, 0, 10);
        std::vector<std::vector<int>> g(n);
        std::vector<std::vector<char>> a(n, std::vector<char>(n, 0));
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if ((int)rnd(rng, 0, 9) < p) {
                    g[i].push_back(j), g[j].push_back(i), a[i][j] = a[j][i] = 1;
                    if (rnd(rng, 0, 5) == 0) g[i].push_back(j), g[j].push_back(i);  // 多重辺
                }
        if (rnd(rng, 0, 3) == 0) g[0].push_back(0);  // 自己ループ
        const int s = (int)rnd(rng, 0, n - 1);
        // 愚直：補グラフを作って BFS
        std::vector<int> d(n, -1);
        std::queue<int> q;
        d[s] = 0, q.push(s);
        while (!q.empty()) {
            const int v = q.front();
            q.pop();
            for (int u = 0; u < n; ++u)
                if (u != v && !a[v][u] && d[u] < 0) d[u] = d[v] + 1, q.push(u);
        }
        CB_CHECK(complement_bfs(g, s) == d);
        // 成分：同じ成分 ⇔ 補グラフで届く
        const auto c = complement_components(g);
        for (int x = 0; x < n; ++x) {
            std::vector<int> dd(n, -1);
            std::queue<int> qq;
            dd[x] = 0, qq.push(x);
            while (!qq.empty()) {
                const int v = qq.front();
                qq.pop();
                for (int u = 0; u < n; ++u)
                    if (u != v && !a[v][u] && dd[u] < 0) dd[u] = 0, qq.push(u);
            }
            for (int y = 0; y < n; ++y) CB_CHECK((c[x] == c[y]) == (dd[y] >= 0));
        }
        std::vector<int> cs = c;
        std::sort(cs.begin(), cs.end());
        cs.erase(std::unique(cs.begin(), cs.end()), cs.end());
        CB_CHECK(complement_component_count(g) == (int)cs.size() && cs.back() == (int)cs.size() - 1);
        cases += 3;
    }
    std::cout << "[graph_complement_bfs] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_complement_bfs] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000, M = N;
    std::vector<std::vector<int>> g(N);
    for (int i = 0; i < M; ++i) {
        const int a = (int)rnd(rng, 0, N - 1), b = (int)rnd(rng, 0, N - 1);
        g[a].push_back(b), g[b].push_back(a);
    }
    const auto t0 = std::chrono::steady_clock::now();
    const auto d = complement_bfs(g, 0);
    const int c = complement_component_count(g);
    ll chk = c;
    for (int x : d) chk += x;
    std::cout << "[graph_complement_bfs] benchmark (N = " << N << ", M = " << M << "): complement_bfs + complement_component_count " << elapsed_ms(t0)
              << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace compbfs_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace compbfs_demo {

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
#define SHOW(expr, comment) compbfs_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    compbfs_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("補グラフの BFS・連結成分（graph_complement_bfs） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  complement_bfs(g, s) / complement_components(g) / complement_component_count(g)\n";

    api("【1】 complement_bfs(g, s)", "vector<int>", "元のグラフ g で隣り合っていない頂点どうしが 1 歩で行けるときの最短距離", "O(N + M)");
    RUN(const std::vector<std::vector<int>> g = {{1, 2, 3}, {0}, {0}, {0}});  // 0 は 1, 2, 3 と仲が悪い
    SHOW(complement_bfs(g, 0), "0 からはどこへも行けない");
    SHOW(complement_bfs(g, 1), "1 → 2, 3 は 1 歩");

    api("【2】 complement_components(g) / complement_component_count(g)", "vector<int> / int", "補グラフの連結成分の番号 / 数", "O(N + M)");
    SHOW(complement_components(g), "");
    SHOW(complement_component_count(g), "{0} と {1, 2, 3}");
}

}  // namespace compbfs_demo

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
    compbfs_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    compbfs_test::benchmark();    // 速度計測
    compbfs_demo::demo();         // 使い方デモ
    return 0;
}
