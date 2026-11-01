// =============================================================================
//  辺彩色（二部グラフは最大次数 Δ 色・一般の単純グラフは Δ + 1 色で、同じ頂点の辺が同じ色にならないように塗る）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  二部グラフ：辺 u–v を足すとき、u で空いている色 a・v で空いている色 b が違えば、v から a, b を交互にたどる道の色を
//              入れ替えて a を空ける（二部なので u には戻らない）。ちょうど Δ 色（König の定理）。多重辺があってもよい
//  一般のグラフ：Misra–Gries の方法（扇を回して色をずらす）で Δ + 1 色（Vizing の定理）。単純グラフのみ
//    例: 「各人が 1 日に 1 試合だけの総当たり戦を最短何日で組めるか」、「2 つのグループの組を、1 日に 1 人 1 組ずつで何日か」
//
//  機能一覧（辺の色は 0 から。答えは edges の順に色を並べたもの）                                計算量
//    edge_coloring_bipartite(L, R, edges)   左 L 個・右 R 個の二部グラフ（辺は {左の番号, 右の番号}）を Δ 色で   O(M (L + R))
//    edge_coloring(n, edges)                単純な無向グラフ（辺は {u, v}）を Δ + 1 色以下で                  O(N M)
//
//  注意: 頂点の数 × 色の数 の表を作る（N = 2×10^3、Δ = 10^3 で 約 8 MB）。edge_coloring は多重辺・自己ループがあると止める
//        木（二部グラフ）なら edge_coloring_bipartite で Δ 色（146D.txt の木の辺彩色と同じ答えの数）
//  速さの目安: edge_coloring_bipartite（左右 1000 頂点、M = 10^5）約 0.03 秒、edge_coloring（N = 1000、M = 10^4）約 0.002 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    edgecolor_test::random_test()  … 色が正しく塗れているか・色の数を調べる
//    edgecolor_test::benchmark()    … 速度計測
//    edgecolor_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace edgecolor_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace edgecolor_internal

inline std::vector<int> edge_coloring_bipartite(int L, int R, const std::vector<std::pair<int, int>>& edges) {
    const int n = L + R, m = (int)edges.size();
    std::vector<int> deg(n, 0);
    for (auto [a, b] : edges) {
        if (a < 0 || a >= L || b < 0 || b >= R) edgecolor_internal::fail("edge_coloring_bipartite", "頂点番号が範囲外");
        ++deg[a], ++deg[L + b];
    }
    const int D = n ? *std::max_element(deg.begin(), deg.end()) : 0;
    std::vector<int> col(m, -1);
    if (D == 0) return col;
    std::vector<int> at((size_t)n * D, -1);  // at[v * D + c] = v で色 c の辺の番号
    auto other = [&](int e, int v) { return edges[e].first == v ? L + edges[e].second : edges[e].first; };
    auto slot = [&](int v, int c) -> int& { return at[(size_t)v * D + c]; };
    auto free_color = [&](int v) {
        for (int c = 0; c < D; ++c)
            if (slot(v, c) < 0) return c;
        edgecolor_internal::fail("edge_coloring_bipartite", "空いている色が無い（内部エラー）");
    };
    std::vector<int> path;
    for (int e = 0; e < m; ++e) {
        const int u = edges[e].first, v = L + edges[e].second;
        const int a = free_color(u), b = free_color(v);
        if (slot(v, a) >= 0) {
            // v から a, b を交互にたどる道（色 a の辺から）の色を入れ替える
            path.clear();
            int cur = v, c = a;
            while (slot(cur, c) >= 0) {
                const int f = slot(cur, c);
                path.push_back(f);
                cur = other(f, cur);
                c = c == a ? b : a;
            }
            for (int f : path) {  // いったん外す
                const int x = edges[f].first, y = L + edges[f].second;
                slot(x, col[f]) = -1, slot(y, col[f]) = -1;
            }
            for (int f : path) {
                col[f] = col[f] == a ? b : a;
                const int x = edges[f].first, y = L + edges[f].second;
                slot(x, col[f]) = f, slot(y, col[f]) = f;
            }
        }
        col[e] = a;
        slot(u, a) = e, slot(v, a) = e;
    }
    return col;
}

// Misra–Gries（KACTL の EdgeColoring と同じ流れ）
inline std::vector<int> edge_coloring(int N, const std::vector<std::pair<int, int>>& eds) {
    {
        std::set<std::pair<int, int>> seen;
        for (auto [a, b] : eds) {
            if (a < 0 || a >= N || b < 0 || b >= N) edgecolor_internal::fail("edge_coloring", "頂点番号が範囲外");
            if (a == b) edgecolor_internal::fail("edge_coloring", "自己ループがある（単純グラフにする）");
            if (!seen.insert({std::min(a, b), std::max(a, b)}).second) edgecolor_internal::fail("edge_coloring", "多重辺がある（単純グラフにする）");
        }
    }
    std::vector<int> cc(N + 1, 0), ret(eds.size(), 0), fan(N + 1), freec(N, 0), loc;
    for (auto [a, b] : eds) ++cc[a], ++cc[b];
    if (eds.empty()) return ret;
    const int ncols = *std::max_element(cc.begin(), cc.end()) + 1;
    cc.assign(std::max(N + 1, ncols + 1), 0);
    std::vector<std::vector<int>> adj(N, std::vector<int>(ncols, -1));  // adj[v][c] = 色 c の辺の相手
    for (auto [u0, v0] : eds) {
        int u = u0, v = v0;
        fan[0] = v;
        loc.assign(ncols, 0);
        int at = u, end = u, d, c = freec[u], ind = 0, i = 0;
        while (d = freec[v], !loc[d] && (v = adj[u][d]) != -1) loc[d] = ++ind, cc[ind] = d, fan[ind] = v;
        cc[loc[d]] = c;
        for (int cd = d; at != -1; cd ^= c ^ d, at = adj[at][cd]) std::swap(adj[at][cd], adj[end = at][cd ^ c ^ d]);
        while (adj[fan[i]][d] != -1) {
            const int left = fan[i], right = fan[++i], e = cc[i];
            adj[u][e] = left, adj[left][e] = u, adj[right][e] = -1, freec[right] = e;
        }
        adj[u][d] = fan[i], adj[fan[i]][d] = u;
        for (int y : {fan[0], u, end}) {
            int& z = freec[y] = 0;
            while (adj[y][z] != -1) ++z;
        }
    }
    for (size_t k = 0; k < eds.size(); ++k) {
        const int u = eds[k].first, v = eds[k].second;
        while (adj[u][ret[k]] != v) ++ret[k];
    }
    return ret;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace edgecolor_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define EC_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_edge_coloring] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
        // 二部グラフ（多重辺あり）
        {
            const int L = (int)rnd(rng, 0, 8), R = (int)rnd(rng, 0, 8), m = (L && R) ? (int)rnd(rng, 0, 30) : 0;
            std::vector<std::pair<int, int>> es(m);
            for (auto& e : es) e = {(int)rnd(rng, 0, L - 1), (int)rnd(rng, 0, R - 1)};
            const auto col = edge_coloring_bipartite(L, R, es);
            std::vector<int> deg(L + R, 0);
            for (auto [a, b] : es) ++deg[a], ++deg[L + b];
            const int D = deg.empty() ? 0 : *std::max_element(deg.begin(), deg.end());
            std::set<std::pair<int, int>> used;
            for (int e = 0; e < m; ++e) {
                EC_CHECK(0 <= col[e] && col[e] < D);
                EC_CHECK(used.insert({es[e].first, col[e]}).second && used.insert({L + es[e].second, col[e]}).second);
            }
            ++cases;
        }
        // 一般の単純グラフ
        {
            const int n = (int)rnd(rng, 0, 10), p = (int)rnd(rng, 0, 10);
            std::vector<std::pair<int, int>> es;
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j)
                    if ((int)rnd(rng, 0, 9) < p) es.push_back(rnd(rng, 0, 1) ? std::make_pair(i, j) : std::make_pair(j, i));
            std::shuffle(es.begin(), es.end(), rng);
            const auto col = edge_coloring(n, es);
            std::vector<int> deg(n, 0);
            for (auto [a, b] : es) ++deg[a], ++deg[b];
            const int D = n ? *std::max_element(deg.begin(), deg.end()) : 0;
            std::set<std::pair<int, int>> used;
            for (size_t e = 0; e < es.size(); ++e) {
                EC_CHECK(0 <= col[e] && col[e] <= D);
                EC_CHECK(used.insert({es[e].first, col[e]}).second && used.insert({es[e].second, col[e]}).second);
            }
            ++cases;
        }
    }
    std::cout << "[graph_edge_coloring] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_edge_coloring] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int L = TEST_SCALE == 1 ? 300 : 1000, M = TEST_SCALE == 1 ? 10000 : 100000;
    std::vector<std::pair<int, int>> es(M);
    for (auto& e : es) e = {(int)rnd(rng, 0, L - 1), (int)rnd(rng, 0, L - 1)};
    auto t0 = std::chrono::steady_clock::now();
    const auto c1 = edge_coloring_bipartite(L, L, es);
    const ll t1 = elapsed_ms(t0);
    const int n = TEST_SCALE == 1 ? 300 : 1000, m = 10 * n;
    std::set<std::pair<int, int>> s;
    while ((int)s.size() < m) {
        int a = (int)rnd(rng, 0, n - 1), b = (int)rnd(rng, 0, n - 1);
        if (a == b) continue;
        s.insert({std::min(a, b), std::max(a, b)});
    }
    std::vector<std::pair<int, int>> ge(s.begin(), s.end());
    t0 = std::chrono::steady_clock::now();
    const auto c2 = edge_coloring(n, ge);
    const ll t2 = elapsed_ms(t0);
    std::cout << "[graph_edge_coloring] benchmark: edge_coloring_bipartite(左右 " << L << ", M = " << M << ") " << t1 << " ms（" << *std::max_element(c1.begin(), c1.end()) + 1
              << " 色） / edge_coloring(N = " << n << ", M = " << m << ") " << t2 << " ms（" << *std::max_element(c2.begin(), c2.end()) + 1 << " 色）\n";
}

}  // namespace edgecolor_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace edgecolor_demo {

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
#define SHOW(expr, comment) edgecolor_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    edgecolor_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("辺彩色（graph_edge_coloring） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  edge_coloring_bipartite(L, R, edges{左, 右}) / edge_coloring(n, edges{u, v})   → 辺ごとの色（0 から）\n";

    api("【1】 edge_coloring_bipartite(L, R, edges)", "vector<int>", "二部グラフを最大次数 Δ 色で（多重辺もよい）", "O(M (L + R))");
    RUN(const std::vector<std::pair<int, int>> be = {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {2, 1}});
    SHOW(edge_coloring_bipartite(3, 2, be), "右の 1 は次数 3 → 3 色");

    api("【2】 edge_coloring(n, edges)", "vector<int>", "単純な一般グラフを Δ + 1 色以下で", "O(N M)");
    code("// 総当たり戦（完全グラフ K4）：1 日に 1 人 1 試合 → 色 = 何日目");
    RUN(const std::vector<std::pair<int, int>> k4 = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}});
    SHOW(edge_coloring(4, k4), "Δ = 3。K4 は 3 色（3 日）で組める");
    RUN(const std::vector<std::pair<int, int>> c3 = {{0, 1}, {1, 2}, {2, 0}});
    SHOW(edge_coloring(3, c3), "三角形は Δ + 1 = 3 色が必要");
}

}  // namespace edgecolor_demo

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
    edgecolor_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    edgecolor_test::benchmark();    // 速度計測
    edgecolor_demo::demo();         // 使い方デモ
    return 0;
}
