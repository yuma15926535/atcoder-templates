// =============================================================================
//  【制限付き最速版】重みが小さい整数のダイクストラ（Dial 法：バケットを輪にして使う。O(N + M + 最短距離の最大値)）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : 辺の重みは 0〜C の整数で、C が小さいこと（目安 C ≤ 100 程度。最短距離の最大値 D の分だけバケットを見るので、
//                           C が大きいと遅くなる）。負の重み・実数の重みは不可
//    できること           : 単一始点・多始点の最短距離と経路（隣接リスト / 隣接を関数で渡す / グリッド）
//    一般版との比較       : graph_grid_v2/02_dijkstra.cpp は O(M log N)（ヒープ）
//                           → この版は O(N + M + D)（D = 最短距離の最大値 ≤ C (N − 1)）。ヒープの log が消える
//                           （重みが 0 か 1 だけなら 01_zero_one_bfs.cpp、全部 1 なら 01_bfs.cpp でも同じオーダー）
//    制限に合わないとき   : 重みが大きい・実数なら graph_grid_v2/02_dijkstra.cpp
//
//  しくみ: 距離 d の頂点をバケット d mod (C + 1) に入れる。距離の小さい順にバケットを見ていけば、ヒープなしで小さい順に取り出せる
//          （今の距離 d から見て、まだのバケットには d〜d + C の距離しか入らないので、C + 1 個を輪にして使い回せる）
//    例: マスに入るコストが 1〜9 のグリッドの最短路（マスが 10^6〜10^7 個）、重みが 1〜10 のグラフ
//
//  機能一覧（頂点は 0〜N−1。到達できない頂点の距離は DIAL_INF = 2^60）                                     計算量
//    dial_dijkstra(g, s) / dial_dijkstra(g, {s1, s2, …})   g = vector<vector<pair<int, W>>>（W は整数型）。C は自動     O(N + M + D)
//    dial_dijkstra_core(n, starts, C, for_each_edge)      隣接を関数で渡す版：for_each_edge(u, push) の中で push(v, w)   O(N + M + D)
//      返り値 DialResult r : r.dist[v] / r.par[v]（最短路の 1 つ前。始点・届かないなら −1）/ r.reached(v) / r.path(t)（始点〜t の頂点列）
//    dial_grid(cost, si, sj) / dial_grid(cost, {{i1, j1}, …})   cost[i][j] = マス (i, j) に入るコスト（0〜C の整数、負なら壁）。上下左右
//      返り値 vector<vector<long long>>（各マスの最短距離。届かなければ DIAL_INF。始点のマスのコストは数えない）   O(HW + D)
//
//  速さの目安: 3000 × 3000 のグリッド（コスト 1〜9）で dial_grid 約 0.36 秒（一般版の priority_queue のダイクストラは 約 1.1 秒）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    dial_test::random_test()  … priority_queue のダイクストラとの突き合わせ
//    dial_test::benchmark()    … 速度計測（priority_queue のダイクストラとの比較つき）
//    dial_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace dial_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 本体：dist（と par）を埋める
template <class F>
void run(int n, const std::vector<int>& starts, int C, F&& for_each_edge, std::vector<long long>& dist, std::vector<int>* par) {
    if (C < 0) fail("dial_dijkstra", "重みの上限 C は 0 以上");
    const long long INF = 1LL << 60;
    dist.assign(n, INF);
    if (par) par->assign(n, -1);
    std::vector<std::vector<int>> bk(C + 1);
    long long pending = 0;
    for (int s : starts) {
        if (s < 0 || s >= n) fail("dial_dijkstra", "始点が範囲外");
        if (dist[s] != 0) dist[s] = 0, bk[0].push_back(s), ++pending;
    }
    long long d = 0;
    int u = -1, cur = 0;  // cur = d mod (C + 1)
    auto push = [&](int v, long long w) {
        if (w < 0 || w > C) fail("dial_dijkstra", "辺の重みが 0〜C の範囲外");
        const long long nd = d + w;
        if (nd < dist[v]) {
            dist[v] = nd;
            if (par) (*par)[v] = u;
            int idx = cur + (int)w;
            if (idx > C) idx -= C + 1;
            bk[idx].push_back(v), ++pending;
        }
    };
    for (; pending > 0; ++d, cur = (cur == C ? 0 : cur + 1)) {
        std::vector<int>& b = bk[cur];
        for (size_t i = 0; i < b.size(); ++i) {  // 重み 0 の辺で同じバケットに増えることがあるので、毎回 size を見る
            u = b[i];
            if (dist[u] != d) continue;  // もっと短い距離で取り出し済み
            for_each_edge(u, push);
        }
        pending -= (long long)b.size();
        b.clear();
    }
}
}  // namespace dial_internal

constexpr long long DIAL_INF = 1LL << 60;

struct DialResult {
    std::vector<long long> dist;  // 届かなければ DIAL_INF
    std::vector<int> par;         // 最短路での 1 つ前（始点・届かないなら −1）
    bool reached(int v) const { return dist[v] < DIAL_INF; }
    std::vector<int> path(int t) const {  // 始点〜t の頂点列（届かなければ空）
        if (!reached(t)) return {};
        std::vector<int> p;
        for (int v = t; v != -1; v = par[v]) p.push_back(v);
        std::reverse(p.begin(), p.end());
        return p;
    }
};

template <class F>
DialResult dial_dijkstra_core(int n, const std::vector<int>& starts, int C, F&& for_each_edge) {
    DialResult r;
    dial_internal::run(n, starts, C, for_each_edge, r.dist, &r.par);
    return r;
}
template <class W>
DialResult dial_dijkstra(const std::vector<std::vector<std::pair<int, W>>>& g, const std::vector<int>& starts) {
    long long C = 0;
    for (const auto& es : g)
        for (const auto& e : es) {
            if (e.second < 0) dial_internal::fail("dial_dijkstra", "負の重みがある");
            C = std::max<long long>(C, (long long)e.second);
        }
    if (C > 100000000) dial_internal::fail("dial_dijkstra", "重みが大きすぎる（graph_grid_v2/02_dijkstra.cpp を使う）");
    return dial_dijkstra_core((int)g.size(), starts, (int)C, [&](int u, auto& push) {
        for (const auto& e : g[u]) push(e.first, (long long)e.second);
    });
}
template <class W>
DialResult dial_dijkstra(const std::vector<std::vector<std::pair<int, W>>>& g, int s) {
    return dial_dijkstra(g, std::vector<int>{s});
}

inline std::vector<std::vector<long long>> dial_grid(const std::vector<std::vector<int>>& cost, const std::vector<std::pair<int, int>>& starts) {
    const int H = (int)cost.size(), W = H ? (int)cost[0].size() : 0;
    int C = 0;
    for (const auto& row : cost) {
        if ((int)row.size() != W) dial_internal::fail("dial_grid", "行の長さがそろっていない");
        for (int c : row) C = std::max(C, c);
    }
    std::vector<int> flat((size_t)H * W), st;
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) flat[(size_t)i * W + j] = cost[i][j];
    for (auto [i, j] : starts) {
        if (i < 0 || i >= H || j < 0 || j >= W) dial_internal::fail("dial_grid", "始点が範囲外");
        st.push_back(i * W + j);
    }
    std::vector<long long> dist;
    dial_internal::run(H * W, st, C, [&](int u, auto& push) {
        const int i = u / W, j = u - i * W;
        if (i > 0 && flat[u - W] >= 0) push(u - W, flat[u - W]);
        if (i + 1 < H && flat[u + W] >= 0) push(u + W, flat[u + W]);
        if (j > 0 && flat[u - 1] >= 0) push(u - 1, flat[u - 1]);
        if (j + 1 < W && flat[u + 1] >= 0) push(u + 1, flat[u + 1]);
    }, dist, nullptr);
    std::vector<std::vector<long long>> res(H, std::vector<long long>(W));
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) res[i][j] = dist[(size_t)i * W + j];
    return res;
}
inline std::vector<std::vector<long long>> dial_grid(const std::vector<std::vector<int>>& cost, int si, int sj) {
    return dial_grid(cost, std::vector<std::pair<int, int>>{{si, sj}});
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace dial_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define DIAL_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_dial] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 一般版：priority_queue のダイクストラ
static std::vector<ll> pq_dijkstra(const std::vector<std::vector<std::pair<int, ll>>>& g, const std::vector<int>& starts) {
    std::vector<ll> dist(g.size(), DIAL_INF);
    std::priority_queue<std::pair<ll, int>, std::vector<std::pair<ll, int>>, std::greater<>> pq;
    for (int s : starts)
        if (dist[s] != 0) dist[s] = 0, pq.push({0, s});
    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (d != dist[u]) continue;
        for (auto [v, w] : g[u])
            if (d + w < dist[v]) dist[v] = d + w, pq.push({dist[v], v});
    }
    return dist;
}

static void random_test() {
    std::mt19937_64 rng(20261101);
    long long cases = 0;
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 30), m = (int)rnd(rng, 0, 80), C = (int)rnd(rng, 0, 6);
        std::vector<std::vector<std::pair<int, ll>>> g(n);
        for (int e = 0; e < m; ++e) g[rnd(rng, 0, n - 1)].push_back({(int)rnd(rng, 0, n - 1), rnd(rng, 0, C)});
        std::vector<int> st;
        const int k = (int)rnd(rng, 1, 3);
        for (int i = 0; i < k; ++i) st.push_back((int)rnd(rng, 0, n - 1));
        const auto exp = pq_dijkstra(g, st);
        const DialResult r = (k == 1 && rnd(rng, 0, 1)) ? dial_dijkstra(g, st[0]) : dial_dijkstra(g, st);
        DIAL_CHECK(r.dist == exp);
        for (int v = 0; v < n; ++v) {
            DIAL_CHECK(r.reached(v) == (exp[v] < DIAL_INF));
            const auto p = r.path(v);
            if (!r.reached(v)) {
                DIAL_CHECK(p.empty() && r.par[v] == -1);
                continue;
            }
            DIAL_CHECK(!p.empty() && p.back() == v && std::count(st.begin(), st.end(), p.front()) > 0);
            ll len = 0;
            for (size_t i = 0; i + 1 < p.size(); ++i) {  // 経路の辺が実在して長さが合う
                ll best = DIAL_INF;
                for (auto [to, w] : g[p[i]])
                    if (to == p[i + 1]) best = std::min(best, w);
                DIAL_CHECK(best < DIAL_INF);
                len += best;
            }
            DIAL_CHECK(len == exp[v]);
        }
        // グリッド
        const int H = (int)rnd(rng, 1, 7), W = (int)rnd(rng, 1, 7);
        std::vector<std::vector<int>> cost(H, std::vector<int>(W));
        for (auto& row : cost)
            for (auto& c : row) c = (int)rnd(rng, -1, 5);
        const int si = (int)rnd(rng, 0, H - 1), sj = (int)rnd(rng, 0, W - 1);
        const auto gd = dial_grid(cost, si, sj);
        std::vector<std::vector<std::pair<int, ll>>> gg(H * W);
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) {
                const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
                for (int d = 0; d < 4; ++d) {
                    const int a = i + di[d], b = j + dj[d];
                    if (0 <= a && a < H && 0 <= b && b < W && cost[a][b] >= 0) gg[i * W + j].push_back({a * W + b, cost[a][b]});
                }
            }
        const auto e2 = pq_dijkstra(gg, {si * W + sj});
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) DIAL_CHECK(gd[i][j] == e2[i * W + j]);
        ++cases;
    }
    std::cout << "[graph_dial] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_dial] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int H = TEST_SCALE == 1 ? 500 : 3000, W = H;
    std::vector<std::vector<int>> cost(H, std::vector<int>(W));
    for (auto& row : cost)
        for (auto& c : row) c = (int)rnd(rng, 1, 9);
    auto t0 = std::chrono::steady_clock::now();
    const auto d1 = dial_grid(cost, 0, 0);
    const ll ms1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    std::vector<ll> dist((size_t)H * W, DIAL_INF);
    {
        std::priority_queue<std::pair<ll, int>, std::vector<std::pair<ll, int>>, std::greater<>> pq;
        dist[0] = 0, pq.push({0, 0});
        while (!pq.empty()) {
            const ll d = pq.top().first;
            const int u = pq.top().second;
            pq.pop();
            if (d != dist[u]) continue;
            const int i = u / W, j = u % W;
            auto relax = [&](int v, int a, int b) {
                if (d + cost[a][b] < dist[v]) dist[v] = d + cost[a][b], pq.push({dist[v], v});
            };
            if (i > 0) relax(u - W, i - 1, j);
            if (i + 1 < H) relax(u + W, i + 1, j);
            if (j > 0) relax(u - 1, i, j - 1);
            if (j + 1 < W) relax(u + 1, i, j + 1);
        }
    }
    const ll ms2 = elapsed_ms(t0);
    DIAL_CHECK(d1[H - 1][W - 1] == dist[(size_t)H * W - 1]);
    std::cout << "[graph_dial] benchmark: " << H << " × " << W << " のグリッド（コスト 1〜9）  dial_grid " << ms1 << " ms / 一般版（priority_queue）" << ms2
              << " ms  (check " << d1[H - 1][W - 1] << ")\n";
}

}  // namespace dial_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace dial_demo {

using ll = long long;
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
#define SHOW(expr, comment) dial_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    dial_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】重みが小さい整数のダイクストラ（graph_dial） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 重みは 0〜C の小さい整数だけ（負・実数・大きい重みは不可）\n"
              << "■ できること: 最短距離・経路を O(N + M + D)（ヒープの log なし）。隣接リスト / 関数で隣接を渡す / グリッド\n"
              << "■ API 一覧\n"
              << "  dial_dijkstra(g, s / starts) / dial_dijkstra_core(n, starts, C, for_each_edge) → DialResult（dist / par / reached / path）\n"
              << "  dial_grid(cost, si, sj / starts) → vector<vector<long long>>\n";

    api("【1】 dial_dijkstra(g, s)", "DialResult", "g = vector<vector<pair<int, W>>>（重みは 0〜C の整数。C は自動で決まる）の最短距離・経路", "O(N + M + D)");
    RUN(std::vector<std::vector<std::pair<int, ll>>> g = {{{1, 4}, {2, 1}}, {{3, 1}}, {{1, 2}, {3, 5}}, {}, {{0, 1}}});
    RUN(DialResult r = dial_dijkstra(g, 0));
    SHOW(r.dist[3], "0 → 2 → 1 → 3 で 1 + 2 + 1");
    SHOW(r.path(3), "経路");
    SHOW(r.par[1], "1 の前は 2");
    SHOW(r.reached(4), "4 へは行けない");
    SHOW(r.dist[4] == DIAL_INF, "届かないときは DIAL_INF");

    api("【2】 dial_dijkstra(g, {s1, s2, …})", "DialResult", "多始点（一番近い始点からの距離）", "O(N + M + D)");
    SHOW(dial_dijkstra(g, std::vector<int>{0, 4}).dist, "");

    api("【3】 dial_dijkstra_core(n, starts, C, for_each_edge)", "DialResult", "for_each_edge(u, push) の中で push(v, w) と隣接を渡す（辺を作らなくてよい）", "O(N + M + D)");
    {
        code("// 例：頂点 v から v + 1 へ重み 2、v から 2v へ重み 3（頂点は 0〜20）。1 から 20 への最短");
        auto r2 = dial_dijkstra_core(21, {1}, 3, [&](int v, auto& push) {
            if (v + 1 <= 20) push(v + 1, 2);
            if (2 * v <= 20) push(2 * v, 3);
        });
        note("dist[20] = " + to_s(r2.dist[20]) + "、経路 " + to_s(r2.path(20)));
    }

    api("【4】 dial_grid(cost, si, sj)", "vector<vector<long long>>", "cost[i][j] = マスに入るコスト（負なら壁）。上下左右に動く最短距離（始点のコストは数えない）", "O(HW + D)");
    RUN(std::vector<std::vector<int>> cost = {{1, 9, 1}, {1, -1, 1}, {1, 1, 1}});
    SHOW(dial_grid(cost, 0, 0), "真ん中は壁（DIAL_INF）。右上へは下を回る方が近い");
    SHOW(dial_grid(cost, {{0, 0}, {0, 2}})[2][1], "多始点（左上と右上）");
}

}  // namespace dial_demo

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
    dial_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    dial_test::benchmark();    // 速度計測
    dial_demo::demo();         // 使い方デモ
    return 0;
}
