// =============================================================================
//  支配木（Dominator tree：始点から各頂点へ行くとき必ず通る頂点。Lengauer–Tarjan）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  有向グラフで、始点 r から v へのどの道も頂点 d を通るとき「d は v を支配する」。v を支配する頂点のうち v に一番近いもの
//  idom(v) を親にすると木になる（支配木）。v の支配木での祖先 = v へ行くのに必ず通る頂点。
//  DFS の順番と「半支配点」を使う Lengauer–Tarjan の方法。O(M log N)。再帰を使わない
//    例: 「どの頂点を壊すとゴールに行けなくなるか」= ゴールの支配木での祖先、「必ず通る中継点の数」
//
//  機能一覧（g は有向グラフの隣接リスト vector<vector<int>>。頂点は 0〜N−1）                          計算量
//    dominator_tree(g, root)             idom[v]（根は root 自身、root から届かない頂点は −1）               O(M log N)
//    dominators_of(idom, root, v)        v を支配する頂点の一覧（v から根へ。v と根も含む。届かないなら空）   O(深さ)
//
//  注意: 無向グラフの「この頂点を消すと分かれる」は関節点（07_lowlink.cpp）。最短路で必ず通る頂点は 02_shortest_path_dag.cpp
//  速さの目安: N = 2×10^5、M = 5×10^5 で 約 0.2 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    domtree_test::random_test()  … 頂点を 1 つずつ消して届くか調べる愚直解との突き合わせ
//    domtree_test::benchmark()    … 速度計測
//    domtree_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace domtree_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace domtree_internal

inline std::vector<int> dominator_tree(const std::vector<std::vector<int>>& g, int root) {
    const int n = (int)g.size();
    if (root < 0 || root >= n) domtree_internal::fail("dominator_tree", "root が範囲外");
    std::vector<int> dfn(n, -1), ord, par(n, -1);
    {  // 反復の DFS（行きがけ順）
        std::vector<std::pair<int, int>> st{{root, 0}};
        dfn[root] = 0, ord.push_back(root);
        while (!st.empty()) {
            auto& [v, k] = st.back();
            if (k < (int)g[v].size()) {
                const int u = g[v][k++];
                if (u < 0 || u >= n) domtree_internal::fail("dominator_tree", "頂点番号が範囲外");
                if (dfn[u] < 0) dfn[u] = (int)ord.size(), ord.push_back(u), par[u] = v, st.push_back({u, 0});
            } else {
                st.pop_back();
            }
        }
    }
    const int m = (int)ord.size();
    std::vector<std::vector<int>> rg(n), bucket(n);
    for (int v = 0; v < n; ++v)
        if (dfn[v] >= 0)
            for (int u : g[v])
                if (dfn[u] >= 0) rg[u].push_back(v);
    std::vector<int> semi(n), idom(n, -1), anc(n, -1), label(n);
    for (int v = 0; v < n; ++v) semi[v] = dfn[v], label[v] = v;
    std::vector<int> stk;
    auto eval = [&](int v) {  // 森の根の手前までで semi が最小の頂点（経路圧縮。再帰なし）
        if (anc[v] < 0) return v;
        stk.clear();
        for (int x = v; anc[anc[x]] >= 0; x = anc[x]) stk.push_back(x);
        for (size_t k = stk.size(); k-- > 0;) {
            const int x = stk[k], a = anc[x];
            if (semi[label[a]] < semi[label[x]]) label[x] = label[a];
            anc[x] = anc[a];
        }
        return label[v];
    };
    for (int i = m - 1; i >= 1; --i) {
        const int w = ord[i];
        for (int v : rg[w]) {
            const int u = eval(v);
            if (semi[u] < semi[w]) semi[w] = semi[u];
        }
        bucket[ord[semi[w]]].push_back(w);
        anc[w] = par[w];
        for (int v : bucket[par[w]]) {
            const int u = eval(v);
            idom[v] = semi[u] < semi[v] ? u : par[w];
        }
        bucket[par[w]].clear();
    }
    for (int i = 1; i < m; ++i) {
        const int w = ord[i];
        if (idom[w] != ord[semi[w]]) idom[w] = idom[idom[w]];
    }
    idom[root] = root;
    return idom;
}
inline std::vector<int> dominators_of(const std::vector<int>& idom, int root, int v) {
    std::vector<int> res;
    if (v < 0 || v >= (int)idom.size() || idom[v] < 0) return res;
    for (int x = v;; x = idom[x]) {
        res.push_back(x);
        if (x == root) break;
    }
    return res;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace domtree_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define DT_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_dominator_tree] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
        const int n = (int)rnd(rng, 1, 10), m = (int)rnd(rng, 0, 25);
        std::vector<std::vector<int>> g(n);
        for (int i = 0; i < m; ++i) g[rnd(rng, 0, n - 1)].push_back((int)rnd(rng, 0, n - 1));
        const int root = (int)rnd(rng, 0, n - 1);
        const auto idom = dominator_tree(g, root);
        // 愚直：頂点 d を消して届かなくなる頂点 = d が支配する頂点
        auto reach = [&](int banned) {
            std::vector<char> vis(n, 0);
            if (banned == root) return vis;
            std::vector<int> st{root};
            vis[root] = 1;
            while (!st.empty()) {
                const int v = st.back();
                st.pop_back();
                for (int u : g[v])
                    if (u != banned && !vis[u]) vis[u] = 1, st.push_back(u);
            }
            return vis;
        };
        const auto all = reach(-1);
        std::vector<std::vector<char>> dom(n, std::vector<char>(n, 0));  // dom[v][d] = d が v を支配
        for (int d = 0; d < n; ++d) {
            const auto r = reach(d);
            for (int v = 0; v < n; ++v)
                if (all[v] && (v == d || !r[v])) dom[v][d] = 1;
        }
        for (int v = 0; v < n; ++v) {
            if (!all[v]) {
                DT_CHECK(idom[v] == -1 && dominators_of(idom, root, v).empty());
                continue;
            }
            if (v == root) {
                DT_CHECK(idom[v] == root);
                continue;
            }
            // idom(v) = 自分以外の支配点で、他の支配点すべてに支配されるもの
            int exp = -1;
            for (int d = 0; d < n; ++d) {
                if (d == v || !dom[v][d]) continue;
                bool ok = true;
                for (int e = 0; e < n && ok; ++e)
                    if (e != v && dom[v][e] && !dom[d][e]) ok = false;
                if (ok) exp = d;
            }
            DT_CHECK(idom[v] == exp);
            auto ds = dominators_of(idom, root, v);
            std::sort(ds.begin(), ds.end());
            std::vector<int> ed;
            for (int d = 0; d < n; ++d)
                if (dom[v][d]) ed.push_back(d);
            DT_CHECK(ds == ed);
            ++cases;
        }
    }
    std::cout << "[graph_dominator_tree] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_dominator_tree] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000, M = 5 * N / 2;
    std::vector<std::vector<int>> g(N);
    for (int i = 1; i < N; ++i) g[i - 1].push_back(i);  // 長い道（深い DFS）
    for (int i = N - 1; i < M; ++i) g[rnd(rng, 0, N - 1)].push_back((int)rnd(rng, 0, N - 1));
    const auto t0 = std::chrono::steady_clock::now();
    const auto idom = dominator_tree(g, 0);
    ll chk = 0;
    for (int x : idom) chk += x;
    std::cout << "[graph_dominator_tree] benchmark (N = " << N << ", M = " << M << "): " << elapsed_ms(t0) << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace domtree_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace domtree_demo {

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
#define SHOW(expr, comment) domtree_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    domtree_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("支配木（graph_dominator_tree） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  dominator_tree(g, root) / dominators_of(idom, root, v)\n";

    api("【1】 dominator_tree(g, root)", "vector<int>", "idom[v] = v へ行くとき必ず通る頂点のうち、v に一番近いもの（根は自分、届かないなら −1）",
        "O(M log N)");
    RUN(const std::vector<std::vector<int>> g = {{1, 2}, {3}, {3}, {4}, {}, {4}});  // 0→1→3, 0→2→3, 3→4、5 は届かない
    RUN(const auto idom = dominator_tree(g, 0));
    SHOW(idom, "3 へは 1 か 2 を通るので必ず通るのは 0 だけ。4 へは 3 を必ず通る");

    api("【2】 dominators_of(idom, root, v)", "vector<int>", "v へ行くとき必ず通る頂点（v から根へ）", "O(深さ)");
    SHOW(dominators_of(idom, 0, 4), "4 は 3 と 0 を必ず通る");
    SHOW(dominators_of(idom, 0, 5), "届かない");
}

}  // namespace domtree_demo

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
    domtree_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    domtree_test::benchmark();    // 速度計測
    domtree_demo::demo();         // 使い方デモ
    return 0;
}
