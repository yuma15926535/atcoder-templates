// =============================================================================
//  一般グラフの最大マッチング（Edmonds の花：二部グラフでなくてよい）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  マッチング（どの頂点も高々 1 本の辺にしか使わない辺の組）の最大の大きさ。二部グラフなら最大流で求まるが、
//  奇数の閉路があると増加路の探索がうまくいかないので、奇数の閉路（花）を 1 つの頂点に縮めながら探す（Edmonds）。O(N^3)
//    例: 「仲のよい 2 人組を最大何組作れるか」（男女の区別が無い）、グリッドの隣り合うマスの組でも、斜めありなどで二部でないとき
//
//  機能一覧（頂点は 0〜N−1。無向グラフ）                                                        計算量
//    GeneralMatching gm(n)               頂点 n 個、辺なし
//    gm.add_edge(u, v)                   辺 u–v を足す（多重辺・自己ループは無視される）                   O(1)
//    gm.max_matching()                   最大マッチングの大きさ（何度呼んでもよい。辺を足したら続きから）  O(N^3)
//    gm.mate(v)                          max_matching の後：v の相手（いなければ −1）                      O(1)
//    gm.pairs()                          組の一覧 {u, v}（u < v）                                         O(N)
//
//  注意: N ≤ 500 くらい（O(N^3)）。重み付き（最大重みマッチング）はできない（二部グラフなら 08_hungarian.cpp・最小費用流）
//  速さの目安: N = 500 のランダムなグラフで 約 0.001 秒（最悪は O(N^3)）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    gmatch_test::random_test()  … bit DP の愚直解との突き合わせ
//    gmatch_test::benchmark()    … 速度計測
//    gmatch_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
namespace gmatch_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace gmatch_internal

struct GeneralMatching {
    explicit GeneralMatching(int n = 0) : n_(n), g_(n), match_(n, -1) {}
    void add_edge(int u, int v) {
        if (u < 0 || u >= n_ || v < 0 || v >= n_) gmatch_internal::fail("GeneralMatching::add_edge", "頂点番号が範囲外");
        if (u == v) return;
        g_[u].push_back(v), g_[v].push_back(u);
    }
    int max_matching() {
        for (int v = 0; v < n_; ++v)
            if (match_[v] == -1) {
                const int t = find_path(v);
                if (t == -1) continue;
                // 増加路に沿って付け替える
                for (int x = t; x != -1;) {
                    const int pv = p_[x], ppv = match_[pv];
                    match_[x] = pv, match_[pv] = x;
                    x = ppv;
                }
            }
        int c = 0;
        for (int v = 0; v < n_; ++v) c += match_[v] > v;
        return c;
    }
    int mate(int v) const {
        if (v < 0 || v >= n_) gmatch_internal::fail("GeneralMatching::mate", "頂点番号が範囲外");
        return match_[v];
    }
    std::vector<std::pair<int, int>> pairs() const {
        std::vector<std::pair<int, int>> r;
        for (int v = 0; v < n_; ++v)
            if (match_[v] > v) r.push_back({v, match_[v]});
        return r;
    }

   private:
    int n_;
    std::vector<std::vector<int>> g_;
    std::vector<int> match_, p_, base_, q_;
    std::vector<char> used_, blossom_;
    int lca(int a, int b) {
        std::vector<char> seen(n_, 0);
        while (true) {
            a = base_[a];
            seen[a] = 1;
            if (match_[a] == -1) break;
            a = p_[match_[a]];
        }
        while (true) {
            b = base_[b];
            if (seen[b]) return b;
            b = p_[match_[b]];
        }
    }
    void mark_path(int v, int b, int child) {
        while (base_[v] != b) {
            blossom_[base_[v]] = blossom_[base_[match_[v]]] = 1;
            p_[v] = child;
            child = match_[v];
            v = p_[match_[v]];
        }
    }
    int find_path(int root) {
        used_.assign(n_, 0), p_.assign(n_, -1), base_.resize(n_);
        std::iota(base_.begin(), base_.end(), 0);
        used_[root] = 1;
        q_.assign(1, root);
        for (size_t qh = 0; qh < q_.size(); ++qh) {
            const int v = q_[qh];
            for (int to : g_[v]) {
                if (base_[v] == base_[to] || match_[v] == to) continue;
                if (to == root || (match_[to] != -1 && p_[match_[to]] != -1)) {  // 奇数の閉路（花）を縮める
                    const int cb = lca(v, to);
                    blossom_.assign(n_, 0);
                    mark_path(v, cb, to);
                    mark_path(to, cb, v);
                    for (int i = 0; i < n_; ++i)
                        if (blossom_[base_[i]]) {
                            base_[i] = cb;
                            if (!used_[i]) used_[i] = 1, q_.push_back(i);
                        }
                } else if (p_[to] == -1) {
                    p_[to] = v;
                    if (match_[to] == -1) return to;  // 増加路が見つかった
                    used_[match_[to]] = 1;
                    q_.push_back(match_[to]);
                }
            }
        }
        return -1;
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace gmatch_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define GM_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_general_matching] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
        const int n = (int)rnd(rng, 0, 12), p = (int)rnd(rng, 0, 10);
        std::vector<std::vector<char>> a(n, std::vector<char>(n, 0));
        GeneralMatching gm(n);
        std::vector<std::pair<int, int>> later;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if ((int)rnd(rng, 0, 19) < p) {
                    a[i][j] = a[j][i] = 1;
                    if (rnd(rng, 0, 3) == 0) later.push_back({i, j});
                    else gm.add_edge(i, j);
                }
        // 途中で一度求め、辺を足して続きから
        gm.max_matching();
        for (auto [u, v] : later) gm.add_edge(v, u);
        const int got = gm.max_matching();
        // 愚直：dp[mask] = mask の頂点だけでの最大マッチング（一番小さい頂点を使わないか、誰かと組む）
        std::vector<int> dp(1 << n, 0);
        for (int m = 1; m < (1 << n); ++m) {
            const int v = __builtin_ctz(m);
            const int rest = m & (m - 1);
            dp[m] = dp[rest];
            for (int u = v + 1; u < n; ++u)
                if ((rest >> u & 1) && a[v][u]) dp[m] = std::max(dp[m], 1 + dp[rest & ~(1 << u)]);
        }
        GM_CHECK(got == dp[(1 << n) - 1]);
        const auto ps = gm.pairs();
        GM_CHECK((int)ps.size() == got);
        std::vector<int> deg(n, 0);
        for (auto [u, v] : ps) {
            GM_CHECK(u < v && a[u][v] && gm.mate(u) == v && gm.mate(v) == u);
            ++deg[u], ++deg[v];
        }
        for (int v = 0; v < n; ++v) GM_CHECK(deg[v] <= 1 && (deg[v] == 1) == (gm.mate(v) != -1));
        ++cases;
    }
    // 奇数の閉路（三角形 2 つを辺でつないだもの）
    {
        GeneralMatching gm(6);
        for (auto [u, v] : std::vector<std::pair<int, int>>{{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {2, 3}}) gm.add_edge(u, v);
        GM_CHECK(gm.max_matching() == 3);
        ++cases;
    }
    std::cout << "[graph_general_matching] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_general_matching] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 200 : 500, M = N;
    GeneralMatching gm(N);
    for (int i = 0; i < M; ++i) gm.add_edge((int)rnd(rng, 0, N - 1), (int)rnd(rng, 0, N - 1));
    const auto t0 = std::chrono::steady_clock::now();
    const int r = gm.max_matching();
    std::cout << "[graph_general_matching] benchmark (N = " << N << ", M = " << M << "): " << elapsed_ms(t0) << " ms  (マッチング " << r << ")\n";
}

}  // namespace gmatch_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace gmatch_demo {

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
#define SHOW(expr, comment) gmatch_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    gmatch_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("一般グラフの最大マッチング（graph_general_matching） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  GeneralMatching gm(n): add_edge(u, v) / max_matching() / mate(v) / pairs()\n";

    api("【1】 add_edge(u, v) / max_matching()", "void / int", "二部グラフでなくてよい最大マッチング", "O(1) / O(N^3)");
    RUN(GeneralMatching gm(5));
    RUN(gm.add_edge(0, 1); gm.add_edge(1, 2); gm.add_edge(2, 0); gm.add_edge(2, 3); gm.add_edge(3, 4));  // 三角形 0-1-2 と 2-3-4
    SHOW(gm.max_matching(), "5 頂点なので最大 2 組");

    api("【2】 mate(v) / pairs()", "int / vector<pair<int, int>>", "v の相手（−1 は相手なし）/ 組の一覧", "O(1) / O(N)");
    SHOW(gm.pairs(), "");
    SHOW(gm.mate(2), "");
}

}  // namespace gmatch_demo

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
    gmatch_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    gmatch_test::benchmark();    // 速度計測
    gmatch_demo::demo();         // 使い方デモ
    return 0;
}
