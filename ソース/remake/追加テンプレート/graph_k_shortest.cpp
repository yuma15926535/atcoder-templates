// =============================================================================
//  k 番目の最短路（同じ頂点を通ってよい歩道の k 番目まで・単純な道の k 番目まで（Yen））   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  歩道（同じ頂点・辺を何度通ってもよい）なら、ダイクストラで各頂点を「K 回まで」確定させると、
//  t が確定した i 回目の距離が i 番目に短い歩道の長さになる。O(K M log(K M))
//  単純な道（同じ頂点を通らない）は Yen の方法：1 番目の道から少しずつずらした候補を作り、短い順に取る。O(K N × ダイクストラ)
//    例: 「2 番目に短い経路の長さ」（同じ長さの別の道でもよい）、「上位 K 本の経路の長さの和」
//
//  機能一覧（g は重み付き有向グラフ vector<vector<pair<int, T>>>。無向なら両向きに入れる。重みは 0 以上）   計算量
//    k_shortest_walks(g, s, t, K)        s → t の歩道の長さを短い順に K 個まで（足りなければ短くなる）        O(K M log(K M))
//    k_shortest_walks_all(g, s, K)       全部の頂点について、s からの歩道の長さを短い順に K 個まで
//    k_shortest_paths(g, s, t, K)        s → t の単純な道を短い順に K 本まで：{長さ, 頂点の列}（同じ長さは順不同）   O(K N (M log N))
//
//  注意: 歩道は同じ長さでも別の道なら別に数える（辺の並びが違えば別）。重みが 0 の閉路があると同じ長さが何個も出る
//  速さの目安: k_shortest_walks_all（N = 10^5、M = 2×10^5、K = 10）約 0.4 秒、k_shortest_paths（N = 200、M = 2000、K = 20）約 0.005 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    kshort_test::random_test()  … 歩道・単純な道を全部並べる愚直解との突き合わせ
//    kshort_test::benchmark()    … 速度計測
//    kshort_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace kshort_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace kshort_internal

template <class T>
std::vector<std::vector<T>> k_shortest_walks_all(const std::vector<std::vector<std::pair<int, T>>>& g, int s, int K) {
    const int n = (int)g.size();
    if (s < 0 || s >= n) kshort_internal::fail("k_shortest_walks_all", "s が範囲外");
    std::vector<std::vector<T>> res(n);
    if (K <= 0) return res;
    using P = std::pair<T, int>;
    std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
    pq.push({T(0), s});
    while (!pq.empty()) {
        auto [d, v] = pq.top();
        pq.pop();
        if ((int)res[v].size() >= K) continue;
        res[v].push_back(d);
        for (auto [u, w] : g[v]) {
            if (w < T(0)) kshort_internal::fail("k_shortest_walks_all", "重みは 0 以上");
            if ((int)res[u].size() < K) pq.push({d + w, u});
        }
    }
    return res;
}
template <class T>
std::vector<T> k_shortest_walks(const std::vector<std::vector<std::pair<int, T>>>& g, int s, int t, int K) {
    if (t < 0 || t >= (int)g.size()) kshort_internal::fail("k_shortest_walks", "t が範囲外");
    return k_shortest_walks_all(g, s, K)[t];
}

namespace kshort_internal {
// banned_v[v] の頂点と banned_e の辺（from, 何番目）を使わないダイクストラ。s → t の {長さ, 頂点の列}（無ければ空）
template <class T>
std::pair<T, std::vector<int>> dijkstra_path(const std::vector<std::vector<std::pair<int, T>>>& g, int s, int t, const std::vector<char>& banned_v,
                                             const std::set<std::pair<int, int>>& banned_e) {
    const int n = (int)g.size();
    const T INF = std::numeric_limits<T>::max() / 4;
    std::vector<T> dist(n, INF);
    std::vector<int> prv(n, -1);
    using P = std::pair<T, int>;
    std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
    dist[s] = T(0), pq.push({T(0), s});
    while (!pq.empty()) {
        auto [d, v] = pq.top();
        pq.pop();
        if (d != dist[v]) continue;
        if (v == t) break;
        for (int k = 0; k < (int)g[v].size(); ++k) {
            const auto [u, w] = g[v][k];
            if (banned_v[u] || banned_e.count({v, k})) continue;
            if (d + w < dist[u]) dist[u] = d + w, prv[u] = v, pq.push({dist[u], u});
        }
    }
    if (dist[t] >= INF) return {INF, {}};
    std::vector<int> path;
    for (int x = t; x != -1; x = prv[x]) path.push_back(x);
    std::reverse(path.begin(), path.end());
    return {dist[t], path};
}
}  // namespace kshort_internal

template <class T>
std::vector<std::pair<T, std::vector<int>>> k_shortest_paths(const std::vector<std::vector<std::pair<int, T>>>& g, int s, int t, int K) {
    const int n = (int)g.size();
    if (s < 0 || s >= n || t < 0 || t >= n) kshort_internal::fail("k_shortest_paths", "s, t が範囲外");
    std::vector<std::pair<T, std::vector<int>>> A;
    if (K <= 0) return A;
    std::vector<char> nob(n, 0);
    auto first = kshort_internal::dijkstra_path(g, s, t, nob, {});
    if (first.second.empty()) return A;
    A.push_back(first);
    std::set<std::pair<T, std::vector<int>>> B;  // 候補（重複なし）
    // 道の長さを計算（隣り合う頂点の間の最小の辺の重み）
    auto edge_w = [&](int a, int b) {
        T best = std::numeric_limits<T>::max();
        for (auto [u, w] : g[a])
            if (u == b && w < best) best = w;
        return best;
    };
    while ((int)A.size() < K) {
        const auto& last = A.back().second;
        T root_cost = T(0);
        for (size_t i = 0; i + 1 < last.size(); ++i) {
            const int spur = last[i];
            // 同じ根を持つ確定済みの道の、次の辺を禁止する
            std::set<std::pair<int, int>> be;
            for (const auto& [c, p] : A)
                if (p.size() > i && std::equal(last.begin(), last.begin() + i + 1, p.begin()))
                    for (int k = 0; k < (int)g[spur].size(); ++k)
                        if (g[spur][k].first == p[i + 1]) be.insert({spur, k});
            std::vector<char> bv(n, 0);
            for (size_t j = 0; j < i; ++j) bv[last[j]] = 1;  // 根の頂点（spur 以外）は使わない
            auto sp = kshort_internal::dijkstra_path(g, spur, t, bv, be);
            if (!sp.second.empty()) {
                std::vector<int> path(last.begin(), last.begin() + i);
                path.insert(path.end(), sp.second.begin(), sp.second.end());
                B.insert({root_cost + sp.first, path});
            }
            root_cost += edge_w(last[i], last[i + 1]);
        }
        // 確定済みと同じ道は飛ばす
        while (!B.empty()) {
            auto it = B.begin();
            bool dup = false;
            for (const auto& a : A)
                if (a.second == it->second) dup = true;
            if (!dup) break;
            B.erase(it);
        }
        if (B.empty()) break;
        A.push_back(*B.begin());
        B.erase(B.begin());
    }
    return A;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace kshort_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define KS_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_k_shortest] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 6), m = (int)rnd(rng, 0, 10);
        G g(n);
        for (int i = 0; i < m; ++i) g[rnd(rng, 0, n - 1)].push_back({(int)rnd(rng, 0, n - 1), rnd(rng, 1, 5)});
        const int s = (int)rnd(rng, 0, n - 1), t = (int)rnd(rng, 0, n - 1), K = (int)rnd(rng, 1, 8);
        // 歩道：長さ d ちょうどの歩道の数 cnt[v][d] を DP（K 個あれば十分なので K で打ち切る）
        const ll L = 12;
        std::vector<std::vector<ll>> cnt(n, std::vector<ll>(L + 1, 0)), walks(n);
        cnt[s][0] = 1;
        for (ll d = 0; d <= L; ++d)
            for (int v = 0; v < n; ++v) {
                if (!cnt[v][d]) continue;
                for (auto [u, w] : g[v])
                    if (d + w <= L) cnt[u][d + w] = std::min<ll>(K, cnt[u][d + w] + cnt[v][d]);
            }
        for (int v = 0; v < n; ++v)
            for (ll d = 0; d <= L; ++d)
                for (ll c = 0; c < cnt[v][d]; ++c) walks[v].push_back(d);
        const auto all = k_shortest_walks_all(g, s, K);
        for (int v = 0; v < n; ++v) {
            auto w = walks[v];
            std::sort(w.begin(), w.end());
            std::vector<ll> exp;
            for (ll x : w)
                if ((int)exp.size() < K) exp.push_back(x);
            // 長さ L を超えるものは愚直が知らないので、L 以下の分だけ比べる
            std::vector<ll> got;
            for (ll x : all[v])
                if (x <= L) got.push_back(x);
            if ((int)exp.size() == K || (int)got.size() < K) KS_CHECK(got == exp);
            else KS_CHECK(std::equal(got.begin(), got.end(), exp.begin()));
            ++cases;
        }
        KS_CHECK(k_shortest_walks(g, s, t, K) == all[t]);
        // 単純な道：全部並べる
        std::vector<std::pair<ll, std::vector<int>>> paths;
        std::vector<char> used(n, 0);
        std::vector<int> cur{s};
        used[s] = 1;
        auto dfs2 = [&](auto&& self, int v, ll d) -> void {
            if (v == t) {
                paths.push_back({d, cur});
                return;
            }
            for (auto [u, w] : g[v])
                if (!used[u]) used[u] = 1, cur.push_back(u), self(self, u, d + w), cur.pop_back(), used[u] = 0;
        };
        dfs2(dfs2, s, 0);
        // 同じ頂点の列は（多重辺で）1 つにまとめ、短い方
        std::sort(paths.begin(), paths.end(), [](const auto& a, const auto& b) { return a.second != b.second ? a.second < b.second : a.first < b.first; });
        paths.erase(std::unique(paths.begin(), paths.end(), [](const auto& a, const auto& b) { return a.second == b.second; }), paths.end());
        std::vector<ll> costs;
        for (auto& p : paths) costs.push_back(p.first);
        std::sort(costs.begin(), costs.end());
        const auto got = k_shortest_paths(g, s, t, K);
        KS_CHECK((int)got.size() == std::min(K, (int)costs.size()));
        for (size_t i = 0; i < got.size(); ++i) {
            KS_CHECK(got[i].first == costs[i]);
            // 本当に s → t の単純な道で、その長さ
            const auto& p = got[i].second;
            KS_CHECK(p.front() == s && p.back() == t);
            std::vector<int> sp = p;
            std::sort(sp.begin(), sp.end());
            KS_CHECK(std::unique(sp.begin(), sp.end()) == sp.end());
            bool found = false;
            for (auto& q : paths)
                if (q.second == p) found = q.first == got[i].first;
            KS_CHECK(found);
        }
        for (size_t i = 0; i < got.size(); ++i)
            for (size_t j = i + 1; j < got.size(); ++j) KS_CHECK(got[i].second != got[j].second);
        ++cases;
    }
    std::cout << "[graph_k_shortest] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_k_shortest] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 300 : 100000, M = 2 * N;  // デバッグビルドは priority_queue の検査（毎回 O(サイズ)）が重いので小さく
    G g(N);
    for (int i = 0; i < M; ++i) g[rnd(rng, 0, N - 1)].push_back({(int)rnd(rng, 0, N - 1), rnd(rng, 1, 1000000)});
    auto t0 = std::chrono::steady_clock::now();
    const auto r = k_shortest_walks_all(g, 0, 10);
    ll chk = 0;
    for (auto& v : r) chk += (ll)v.size();
    const ll t1 = elapsed_ms(t0);
    const int n2 = TEST_SCALE == 1 ? 40 : 200, m2 = TEST_SCALE == 1 ? 200 : 2000;
    G h(n2);
    for (int i = 0; i < m2; ++i) h[rnd(rng, 0, n2 - 1)].push_back({(int)rnd(rng, 0, n2 - 1), rnd(rng, 1, 1000)});
    t0 = std::chrono::steady_clock::now();
    chk += (ll)k_shortest_paths(h, 0, n2 - 1, 20).size();
    const ll t2 = elapsed_ms(t0);
    std::cout << "[graph_k_shortest] benchmark: k_shortest_walks_all(N = " << N << ", M = " << M << ", K = 10) " << t1 << " ms / k_shortest_paths(N = " << n2
              << ", M = " << m2 << ", K = 20) " << t2 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace kshort_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace kshort_demo {

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
#define SHOW(expr, comment) kshort_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    kshort_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("k 番目の最短路（graph_k_shortest） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  k_shortest_walks(g, s, t, K) / k_shortest_walks_all(g, s, K) / k_shortest_paths(g, s, t, K)\n";

    RUN(std::vector<std::vector<std::pair<int, long long>>> g(4));
    RUN(auto add = [&](int a, int b, long long w) { g[a].push_back({b, w}); g[b].push_back({a, w}); });  // 無向
    RUN(add(0, 1, 1); add(1, 3, 1); add(0, 2, 2); add(2, 3, 2); add(1, 2, 1));

    api("【1】 k_shortest_walks(g, s, t, K)", "vector<T>", "s → t の歩道（同じ頂点を何度通ってもよい）の長さを短い順に K 個", "O(K M log(K M))");
    SHOW(k_shortest_walks(g, 0, 3, 5), "0-1-3 が 2、そのあと 0-1-0-1-3 なども");

    api("【2】 k_shortest_paths(g, s, t, K)", "vector<pair<T, vector<int>>>", "s → t の単純な道（同じ頂点を通らない）を短い順に K 本", "O(K N (M log N))");
    RUN(const auto ps = k_shortest_paths(g, 0, 3, 4));
    for (const auto& p : ps) note("長さ " + to_s(p.first) + "：" + to_s(p.second));

    api("【3】 k_shortest_walks_all(g, s, K)", "vector<vector<T>>", "全部の頂点について、短い順に K 個", "O(K M log(K M))");
    SHOW(k_shortest_walks_all(g, 0, 2), "各頂点の 1 番目と 2 番目");
}

}  // namespace kshort_demo

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
    kshort_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    kshort_test::benchmark();    // 速度計測
    kshort_demo::demo();         // 使い方デモ
    return 0;
}
