// =============================================================================
//  Gomory–Hu 木（無向グラフの全点対の最小カットを、N − 1 回の最大流で）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  無向グラフで「u と v を分けるのに切る辺の容量の和の最小（= u–v の最大流）」を全部の組について答える。
//  Gusfield の方法：頂点 i（1〜N−1）と今の親 p[i] の間で最大流を流し、i 側に残った頂点の親を i にする。
//  できた木の上で、u–v のパスの辺の重みの最小が u–v の最小カット。最大流（Dinic）N − 1 回
//    例: 「全部の組の最小カットの和」「最小カットが k 以上の組の数」「どの 2 頂点もつながりの強さが k 以上か」
//
//  機能一覧（頂点は 0〜N−1。容量は 0 以上の long long）                                         計算量
//    GomoryHu gh(n)                      頂点 n 個、辺なし
//    gh.add_edge(u, v, cap)              無向の辺 u–v（容量 cap）を足す                                  O(1)
//    gh.build()                          木を作る（add_edge の後に 1 回）                                 O(N × 最大流)
//    gh.parent(v) / gh.weight(v)         木での親（根 0 は −1）/ 辺 v–parent(v) の重み                   O(1)
//    gh.min_cut(u, v)                    u と v の最小カット（u = v なら 0 ではなく「無限」= −1 を返す）   O(N)
//    gh.all_pairs()                      全部の組の最小カットの表（対角は −1）                             O(N^2)
//
//  注意: 連結でない組の最小カットは 0。N ≤ 500 くらい（最大流 N 回）。有向グラフには使えない
//  速さの目安: N = 200、M = 2000 で build + all_pairs 約 0.02 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    gomoryhu_test::random_test()  … 頂点の分け方を全部試す愚直解との突き合わせ
//    gomoryhu_test::benchmark()    … 速度計測
//    gomoryhu_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace gomoryhu_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 無向の辺の Dinic（毎回容量を元に戻して使う）
struct Dinic {
    struct E {
        int to, rev;
        long long cap, orig;
    };
    std::vector<std::vector<E>> g;
    std::vector<int> level, it;
    explicit Dinic(int n) : g(n), level(n), it(n) {}
    void add_undirected(int u, int v, long long c) {
        g[u].push_back({v, (int)g[v].size(), c, c});
        g[v].push_back({u, (int)g[u].size() - 1, c, c});
    }
    void reset() {
        for (auto& es : g)
            for (auto& e : es) e.cap = e.orig;
    }
    bool bfs(int s, int t) {
        std::fill(level.begin(), level.end(), -1);
        std::queue<int> q;
        level[s] = 0, q.push(s);
        while (!q.empty()) {
            const int v = q.front();
            q.pop();
            for (const E& e : g[v])
                if (e.cap > 0 && level[e.to] < 0) level[e.to] = level[v] + 1, q.push(e.to);
        }
        return level[t] >= 0;
    }
    long long dfs(int v, int t, long long f) {  // 深さは BFS の段の数（≤ N）
        if (v == t) return f;
        for (int& i = it[v]; i < (int)g[v].size(); ++i) {
            E& e = g[v][i];
            if (e.cap <= 0 || level[e.to] != level[v] + 1) continue;
            const long long d = dfs(e.to, t, std::min(f, e.cap));
            if (d > 0) {
                e.cap -= d, g[e.to][e.rev].cap += d;
                return d;
            }
        }
        return 0;
    }
    long long max_flow(int s, int t) {
        long long flow = 0;
        while (bfs(s, t)) {
            std::fill(it.begin(), it.end(), 0);
            long long f;
            while ((f = dfs(s, t, std::numeric_limits<long long>::max())) > 0) flow += f;
        }
        return flow;
    }
};
}  // namespace gomoryhu_internal

struct GomoryHu {
    explicit GomoryHu(int n = 0) : n_(n), flow_(n) {}
    void add_edge(int u, int v, long long cap) {
        if (u < 0 || u >= n_ || v < 0 || v >= n_) gomoryhu_internal::fail("GomoryHu::add_edge", "頂点番号が範囲外");
        if (cap < 0) gomoryhu_internal::fail("GomoryHu::add_edge", "容量は 0 以上");
        if (u == v) return;
        flow_.add_undirected(u, v, cap);
    }
    void build() {
        par_.assign(n_, 0), w_.assign(n_, 0);
        if (n_ > 0) par_[0] = -1;
        for (int i = 1; i < n_; ++i) {
            flow_.reset();
            w_[i] = flow_.max_flow(i, par_[i]);
            // 残余グラフで i から届く側（level ≥ 0）にいる、同じ親の頂点の親を i に
            for (int j = i + 1; j < n_; ++j)
                if (par_[j] == par_[i] && flow_.level[j] >= 0) par_[j] = i;
        }
        built_ = true;
    }
    int parent(int v) const {
        need("GomoryHu::parent");
        return par_[v];
    }
    long long weight(int v) const {
        need("GomoryHu::weight");
        return v == 0 ? -1 : w_[v];
    }
    long long min_cut(int u, int v) const {
        need("GomoryHu::min_cut");
        if (u == v) return -1;
        // 木の上で u–v のパスの最小（深さをそろえて上がる）
        auto dep = [&](int x) {
            int d = 0;
            for (int y = x; par_[y] >= 0; y = par_[y]) ++d;
            return d;
        };
        int du = dep(u), dv = dep(v);
        long long best = std::numeric_limits<long long>::max();
        while (du > dv) best = std::min(best, w_[u]), u = par_[u], --du;
        while (dv > du) best = std::min(best, w_[v]), v = par_[v], --dv;
        while (u != v) best = std::min({best, w_[u], w_[v]}), u = par_[u], v = par_[v];
        return best;
    }
    std::vector<std::vector<long long>> all_pairs() const {
        need("GomoryHu::all_pairs");
        std::vector<std::vector<long long>> res(n_, std::vector<long long>(n_, -1));
        // 木の上で各頂点から広げる
        std::vector<std::vector<std::pair<int, long long>>> t(n_);
        for (int v = 1; v < n_; ++v) t[v].push_back({par_[v], w_[v]}), t[par_[v]].push_back({v, w_[v]});
        for (int s = 0; s < n_; ++s) {
            std::vector<int> st{s};
            std::vector<char> vis(n_, 0);
            vis[s] = 1;
            std::vector<long long> mn(n_, std::numeric_limits<long long>::max());
            while (!st.empty()) {
                const int x = st.back();
                st.pop_back();
                for (auto [y, w] : t[x])
                    if (!vis[y]) vis[y] = 1, mn[y] = std::min(mn[x], w), res[s][y] = mn[y], st.push_back(y);
            }
        }
        return res;
    }

   private:
    int n_;
    gomoryhu_internal::Dinic flow_;
    std::vector<int> par_;
    std::vector<long long> w_;
    bool built_ = false;
    void need(const char* f) const {
        if (!built_) gomoryhu_internal::fail(f, "先に build() を呼ぶ");
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace gomoryhu_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define GH_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_gomory_hu] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
        const int n = (int)rnd(rng, 1, 8), m = (int)rnd(rng, 0, 14);
        std::vector<std::array<ll, 3>> es(m);
        GomoryHu gh(n);
        for (auto& e : es) {
            e = {rnd(rng, 0, n - 1), rnd(rng, 0, n - 1), rnd(rng, 0, 9)};
            gh.add_edge((int)e[0], (int)e[1], e[2]);
        }
        gh.build();
        const auto ap = gh.all_pairs();
        for (int u = 0; u < n; ++u)
            for (int v = 0; v < n; ++v) {
                if (u == v) {
                    GH_CHECK(gh.min_cut(u, v) == -1 && ap[u][v] == -1);
                    continue;
                }
                // 愚直：u を含み v を含まない頂点集合を全部
                ll best = -1;
                for (int S = 0; S < (1 << n); ++S) {
                    if (!(S >> u & 1) || (S >> v & 1)) continue;
                    ll c = 0;
                    for (auto& e : es)
                        if (((S >> e[0]) & 1) != ((S >> e[1]) & 1)) c += e[2];
                    if (best < 0 || c < best) best = c;
                }
                GH_CHECK(gh.min_cut(u, v) == best && ap[u][v] == best);
                ++cases;
            }
    }
    std::cout << "[graph_gomory_hu] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_gomory_hu] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 60 : 200, M = 10 * N;
    GomoryHu gh(N);
    for (int i = 0; i < M; ++i) gh.add_edge((int)rnd(rng, 0, N - 1), (int)rnd(rng, 0, N - 1), rnd(rng, 1, 1000000));
    const auto t0 = std::chrono::steady_clock::now();
    gh.build();
    ll chk = 0;
    for (const auto& r : gh.all_pairs())
        for (ll x : r) chk += x;
    std::cout << "[graph_gomory_hu] benchmark (N = " << N << ", M = " << M << "): build + all_pairs " << elapsed_ms(t0) << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace gomoryhu_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace gomoryhu_demo {

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
#define SHOW(expr, comment) gomoryhu_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    gomoryhu_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("Gomory–Hu 木（graph_gomory_hu） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  GomoryHu gh(n): add_edge(u, v, cap) → build() → min_cut(u, v) / all_pairs() / parent(v) / weight(v)\n";

    api("【1】 add_edge(u, v, cap) / build() / min_cut(u, v)", "void / void / long long", "全部の組の最小カット（= 最大流）を N − 1 回の最大流で", "O(N × 最大流)");
    RUN(GomoryHu gh(4));
    RUN(gh.add_edge(0, 1, 3); gh.add_edge(1, 2, 1); gh.add_edge(2, 3, 4); gh.add_edge(0, 2, 2));
    RUN(gh.build());
    SHOW(gh.min_cut(0, 3), "{0, 1} と {2, 3} に分ける：1–2 と 0–2 で 1 + 2");
    SHOW(gh.min_cut(0, 1), "1 は 0–1 と 1–2 で 3 + 1");
    SHOW(gh.min_cut(1, 2), "");

    api("【2】 all_pairs() / parent(v) / weight(v)", "表 / int / long long", "全部の組の表 / 木の親と辺の重み", "O(N^2) / O(1)");
    SHOW(gh.all_pairs(), "対角は −1");
    SHOW(gh.parent(3), "");
    SHOW(gh.weight(3), "");
}

}  // namespace gomoryhu_demo

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
    gomoryhu_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    gomoryhu_test::benchmark();    // 速度計測
    gomoryhu_demo::demo();         // 使い方デモ
    return 0;
}
