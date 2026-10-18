// =============================================================================
//  最小費用流  MinCostFlow                                          --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 08 フロー
//  依存       : なし（単独で貼る）
//  頂点番号   : 0-indexed。容量・コストは long long。負のコストの辺も可（ただし負閉路は不可）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  主双対法: ポテンシャル（各頂点の最短距離）でコストを 0 以上に直し、ダイクストラで最短路を見つけて流す。
//  負のコストの辺があるときは、最初だけベルマンフォードでポテンシャルを作る（「最大」重みの問題も符号を変えて解ける）。
//  注意: 流量 × コストが long long に収まること。容量無限の道があるなら limit を必ず指定する
//
//  機能一覧                                                                       計算量
//    MinCostFlow F(n) / F.add_vertex()              頂点 n 個で作る / 頂点を足す
//    F.add_edge(from, to, cap, cost)  -> int       辺を足す（辺番号を返す）
//    F.flow(s, t, limit = CAP_INF)    -> {流量, コスト}   最大 limit まで、コスト最小で流す   O(F (N + M) log M)
//    F.flow_min_cost(s, t)            -> {流量, コスト}   流量は自由で、コストの合計を最小にする
//    F.slope(s, t, limit)             -> 折れ線の頂点 (流量, コスト) の列（流量ごとの最小コスト）
//    F.get_edge(id) / F.edges()       -> {from, to, cap, flow, cost}（流し方の復元）
//    （流量をちょうど K にしたいとき: auto [f, c] = F.flow(s, t, K); if (f < K) 不可能）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mcf_test::random_test() / mcf_test::benchmark() / mcf_demo::demo()
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_MIN_COST_FLOW
#define GG_MIN_COST_FLOW
// ---- 最小費用流（グラフ・グリッド第二世代）  主双対法（ポテンシャル + ダイクストラ）。負のコストの辺も可（負閉路は不可） ----
// 頂点は 0-indexed。容量・コストは long long。流量 × コストが long long に収まること
struct MinCostFlow {
    static constexpr long long CAP_INF = (long long)4e18;
    struct EdgeInfo { int from, to; long long cap, flow, cost; };
    struct E { int to, rev; long long cap, cost; };
    std::vector<std::vector<E>> g;
    std::vector<std::pair<int, int>> pos_;

    MinCostFlow() {}
    explicit MinCostFlow(int n) : g(n) {}
    int size() const { return (int)g.size(); }
    int add_vertex() { g.emplace_back(); return (int)g.size() - 1; }
    // 辺 from → to（容量 cap、1 流すごとにコスト cost）。辺番号を返す
    int add_edge(int from, int to, long long cap, long long cost) {
        const int id = (int)pos_.size();
        pos_.emplace_back(from, (int)g[from].size());
        const int rf = (int)g[from].size(), rt = (int)g[to].size() + (from == to ? 1 : 0);
        g[from].push_back({to, rt, cap, cost});
        g[to].push_back({from, rf, 0, -cost});
        return id;
    }
    // s から t へ最大 limit 流す（コストが最小になるように）。{流した量, コストの合計} を返す
    std::pair<long long, long long> flow(int s, int t, long long limit = CAP_INF) { return slope(s, t, limit).back(); }
    // 流す量を自由に決めて、コストの合計を最小にする（1 単位あたりのコストが負の間だけ流す）。{流量, コスト}
    std::pair<long long, long long> flow_min_cost(int s, int t, long long limit = CAP_INF) { return slope_impl(s, t, limit, true).back(); }
    // 流量とコストの関係（折れ線の頂点 (流量, コスト) の列。最初は (0, 0)）。この間は直線で、傾きは増えていく
    std::vector<std::pair<long long, long long>> slope(int s, int t, long long limit = CAP_INF) { return slope_impl(s, t, limit, false); }
    EdgeInfo get_edge(int id) const {
        const E& e = g[pos_[id].first][pos_[id].second];
        const E& r = g[e.to][e.rev];
        return {pos_[id].first, e.to, e.cap + r.cap, r.cap, e.cost};
    }
    std::vector<EdgeInfo> edges() const {
        std::vector<EdgeInfo> r;
        for (int id = 0; id < (int)pos_.size(); ++id) r.push_back(get_edge(id));
        return r;
    }
    int edge_count() const { return (int)pos_.size(); }

    std::vector<std::pair<long long, long long>> slope_impl(int s, int t, long long limit, bool only_negative) {
        const int n = size();
        const long long INF = (long long)4e18;
        std::vector<std::pair<long long, long long>> res{{0, 0}};
        if (s == t) return res;
        // ポテンシャル h: 残余グラフに負のコストの辺があればベルマンフォードで初期化（なければ 0）
        std::vector<long long> h(n, 0), d(n);
        bool neg = false;
        for (int v = 0; v < n && !neg; ++v)
            for (const E& e : g[v])
                if (e.cap > 0 && e.cost < 0) { neg = true; break; }
        if (neg) {
            std::fill(d.begin(), d.end(), INF);
            d[s] = 0;
            for (int round = 0; round < n; ++round) {
                bool upd = false;
                for (int v = 0; v < n; ++v) {
                    if (d[v] == INF) continue;
                    for (const E& e : g[v])
                        if (e.cap > 0 && d[v] + e.cost < d[e.to]) d[e.to] = d[v] + e.cost, upd = true;
                }
                if (!upd) break;
            }
            for (int v = 0; v < n; ++v) h[v] = (d[v] == INF) ? 0 : d[v];
        }
        std::vector<int> pv(n), pe(n);
        std::vector<char> vis(n);
        struct Item { long long d; int v; };
        std::vector<Item> heap;
        auto push = [&](long long dist, int v) {
            size_t i = heap.size();
            heap.push_back({dist, v});
            while (i > 0) {
                const size_t p = (i - 1) >> 1;
                if (heap[p].d <= dist) break;
                heap[i] = heap[p], i = p;
            }
            heap[i] = {dist, v};
        };
        auto pop = [&]() {
            const Item top = heap[0], x = heap.back();
            heap.pop_back();
            const size_t m = heap.size();
            if (m > 0) {
                size_t i = 0;
                for (size_t c = 1; c < m; c = 2 * i + 1) {
                    if (c + 1 < m && heap[c + 1].d < heap[c].d) ++c;
                    if (x.d <= heap[c].d) break;
                    heap[i] = heap[c], i = c;
                }
                heap[i] = x;
            }
            return top;
        };
        long long fl = 0, cost = 0, prev_unit = 0;
        bool first = true;
        while (fl < limit) {
            // 削減コスト cost + h[u] - h[v]（0 以上）でダイクストラ
            std::fill(d.begin(), d.end(), INF), std::fill(vis.begin(), vis.end(), 0);
            heap.clear();
            d[s] = 0, push(0, s);
            while (!heap.empty()) {
                const Item it = pop();
                const int v = it.v;
                if (vis[v]) continue;
                vis[v] = 1;
                if (v == t) break;  // t の距離が決まったら打ち切る（ポテンシャルの更新で正しさを保つ）
                for (int i = 0; i < (int)g[v].size(); ++i) {
                    const E& e = g[v][i];
                    if (e.cap == 0 || vis[e.to]) continue;
                    const long long nd = d[v] + e.cost + h[v] - h[e.to];
                    if (nd < d[e.to]) d[e.to] = nd, pv[e.to] = v, pe[e.to] = i, push(nd, e.to);
                }
            }
            if (!vis[t]) break;
            // 確定した頂点だけ h[v] += d[v] - d[t]（全体に定数を足しても削減コストは変わらないので、これで 0 以上が保たれる）
            for (int v = 0; v < n; ++v)
                if (vis[v]) h[v] -= d[t] - d[v];
            const long long unit = h[t] - h[s];  // 1 流すごとのコスト
            if (only_negative && unit >= 0) break;
            long long f = limit - fl;
            for (int v = t; v != s; v = pv[v]) f = std::min(f, g[pv[v]][pe[v]].cap);
            for (int v = t; v != s; v = pv[v]) {
                E& e = g[pv[v]][pe[v]];
                e.cap -= f, g[v][e.rev].cap += f;
            }
            fl += f, cost += f * unit;
            if (!first && unit == prev_unit) res.pop_back();
            res.push_back({fl, cost});
            prev_unit = unit, first = false;
        }
        return res;
    }
};
#endif  // GG_MIN_COST_FLOW
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "min_cost_flow";
// ---- テスト用の小道具 ----
#define CHECK(...)                                                                                \
    do {                                                                                          \
        if (!(__VA_ARGS__)) {                                                                     \
            std::cerr << "[" << TEST_TAG << "] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                         \
        }                                                                                         \
    } while (0)
// テスト量の倍率。デバッグビルド（最適化なし・範囲チェックあり）は遅いので自動で 1/8 にする
#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif
[[maybe_unused]] static int rnd(std::mt19937_64& rng, int lo, int hi) { return lo + (int)(rng() % (unsigned long long)(hi - lo + 1)); }
[[maybe_unused]] static long long ms_since(std::chrono::steady_clock::time_point t0) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
struct CinFromString {  // std::cin を文字列に差し替える
    std::istringstream ss;
    std::streambuf* old;
    explicit CinFromString(const std::string& s) : ss(s), old(std::cin.rdbuf(ss.rdbuf())) {}
    ~CinFromString() { std::cin.rdbuf(old); }
};

namespace mcf_test {

struct Ed { int a, b; long long cap, cost; };
// 愚直: 全部の流し方を数え上げ、流量ごとの最小コスト（流せない流量は大きい値）
static std::vector<long long> brute(int n, const std::vector<Ed>& es, int s, int t) {
    const long long BIG = 1LL << 60;
    std::vector<long long> best(64, BIG);
    const int m = (int)es.size();
    std::vector<long long> f(m, 0);
    std::function<void(int)> rec = [&](int i) {
        if (i == m) {
            std::vector<long long> net(n, 0);
            long long c = 0;
            for (int k = 0; k < m; ++k) net[es[k].a] -= f[k], net[es[k].b] += f[k], c += f[k] * es[k].cost;
            for (int v = 0; v < n; ++v)
                if (v != s && v != t && net[v] != 0) return;
            if (net[t] < 0) return;
            best[net[t]] = std::min(best[net[t]], c);
            return;
        }
        for (long long x = 0; x <= es[i].cap; ++x) f[i] = x, rec(i + 1);
    };
    rec(0);
    return best;
}
static long long cost_at(const std::vector<std::pair<long long, long long>>& sl, long long k) {
    for (size_t i = 1; i < sl.size(); ++i)
        if (k <= sl[i].first) {
            const long long f0 = sl[i - 1].first, c0 = sl[i - 1].second, f1 = sl[i].first, c1 = sl[i].second;
            return c0 + (c1 - c0) / (f1 - f0) * (k - f0);
        }
    return sl[0].second;  // k == 0
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 2, 4), m = rnd(rng, 0, 6);
        std::vector<long long> p(n);  // 負閉路ができないよう、コスト = p[b] - p[a] + 0 以上
        for (auto& x : p) x = (it % 2) ? rnd(rng, -4, 4) : 0;
        std::vector<Ed> es;
        MinCostFlow F(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            const long long cap = rnd(rng, 0, 2), cost = p[b] - p[a] + rnd(rng, 0, 5);
            if (a == b && cost < 0) continue;
            es.push_back({a, b, cap, cost});
            CHECK(F.add_edge(a, b, cap, cost) == (int)es.size() - 1);
        }
        int s = rnd(rng, 0, n - 1), t = rnd(rng, 0, n - 1);
        if (s == t) t = (s + 1) % n;
        const std::vector<long long> best = brute(n, es, s, t);
        long long maxf = 0;
        while (maxf + 1 < 64 && best[maxf + 1] < (1LL << 60)) ++maxf;
        const int mode = it % 3;
        if (mode == 0) {
            MinCostFlow G = F;
            const auto sl = G.slope(s, t);
            CHECK(sl.front() == std::make_pair(0LL, 0LL) && sl.back().first == maxf);
            for (long long k = 0; k <= maxf; ++k) CHECK(cost_at(sl, k) == best[k]);
            for (size_t i = 2; i < sl.size(); ++i)  // 傾きは真に増える（同じ傾きはまとめる）
                CHECK((sl[i].second - sl[i - 1].second) * (sl[i - 1].first - sl[i - 2].first) > (sl[i - 1].second - sl[i - 2].second) * (sl[i].first - sl[i - 1].first));
        } else if (mode == 1) {
            const long long lim = rnd(rng, 0, 4);
            const auto r = F.flow(s, t, lim);
            const long long k = std::min(lim, maxf);
            CHECK(r.first == k && r.second == best[k]);
            // 流し方の復元
            std::vector<long long> net(n, 0);
            long long c = 0;
            for (const auto& e : F.edges()) {
                CHECK(0 <= e.flow && e.flow <= e.cap);
                net[e.from] -= e.flow, net[e.to] += e.flow, c += e.flow * e.cost;
            }
            for (int v = 0; v < n; ++v)
                if (v != s && v != t) CHECK(net[v] == 0);
            CHECK(net[t] == k && c == r.second);
            // 続けて残りを流す（2 回目の呼び出しでは残余グラフに負の辺があるのでベルマンフォードから始まる）
            const auto r2 = F.flow(s, t);
            CHECK(r.first + r2.first == maxf && r.second + r2.second == best[maxf]);
        } else {
            const auto r = F.flow_min_cost(s, t);
            long long mc = 0, mk = 0;
            for (long long k = 0; k <= maxf; ++k)
                if (best[k] < mc) mc = best[k], mk = k;
            CHECK(r.second == mc && r.first == mk);
        }
    }
    // 少し大きいグラフ: 愚直な「ベルマンフォードで最短路を見つけて流す」と比べる
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 2, 25), m = rnd(rng, 0, 80);
        std::vector<long long> p(n);
        for (auto& x : p) x = rnd(rng, -50, 50);
        struct BE { int to, rev; long long cap, cost; };
        std::vector<std::vector<BE>> bg(n);
        MinCostFlow F(n);
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            if (a == b) continue;
            const long long cap = rnd(rng, 1, 10), cost = p[b] - p[a] + rnd(rng, 0, 30);
            F.add_edge(a, b, cap, cost);
            bg[a].push_back({b, (int)bg[b].size(), cap, cost});
            bg[b].push_back({a, (int)bg[a].size() - 1, 0, -cost});
        }
        const int s = 0, t = n - 1;
        const long long lim = (it % 2) ? rnd(rng, 0, 30) : (1LL << 60);
        long long bf = 0, bc = 0;
        while (bf < lim) {
            std::vector<long long> d(n, 1LL << 60);
            std::vector<int> pv(n, -1), pe(n, -1);
            d[s] = 0;
            for (int round = 0; round < n; ++round)
                for (int v = 0; v < n; ++v)
                    if (d[v] < (1LL << 60))
                        for (int k = 0; k < (int)bg[v].size(); ++k)
                            if (bg[v][k].cap > 0 && d[v] + bg[v][k].cost < d[bg[v][k].to]) d[bg[v][k].to] = d[v] + bg[v][k].cost, pv[bg[v][k].to] = v, pe[bg[v][k].to] = k;
            if (d[t] >= (1LL << 60)) break;
            long long f = lim - bf;
            for (int v = t; v != s; v = pv[v]) f = std::min(f, bg[pv[v]][pe[v]].cap);
            for (int v = t; v != s; v = pv[v]) {
                BE& e = bg[pv[v]][pe[v]];
                e.cap -= f, bg[v][e.rev].cap += f;
            }
            bf += f, bc += f * d[t];
        }
        const auto r = (it % 2) ? F.flow(s, t, lim) : F.flow(s, t);
        CHECK(r.first == bf && r.second == bc);
    }
    std::cout << "[min_cost_flow] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[min_cost_flow] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {100, 300}) {
        if (TEST_SCALE == 1 && n > 100) break;
        // 割り当て問題（左 n・右 n、全部の組に辺、コスト 0 ~ 10^9）
        MinCostFlow F(2 * n + 2);
        const int s = 2 * n, t = s + 1;
        for (int i = 0; i < n; ++i) F.add_edge(s, i, 1, 0), F.add_edge(n + i, t, 1, 0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) F.add_edge(i, n + j, 1, (long long)(rng() % 1000000000));
        auto t0 = std::chrono::steady_clock::now();
        const auto r = F.flow(s, t);
        std::cout << "[min_cost_flow] benchmark: 割り当て " << n << " × " << n << "（辺 " << n * n << "）  " << ms_since(t0) << " ms  (flow " << r.first << ", cost "
                  << r.second % 1000 << ")\n";
    }
    {
        // 負のコストつき（最大重みマッチングを符号を変えて）
        const int n = 200;
        MinCostFlow F(2 * n + 2);
        const int s = 2 * n, t = s + 1;
        for (int i = 0; i < n; ++i) F.add_edge(s, i, 1, 0), F.add_edge(n + i, t, 1, 0);
        for (int k = 0; k < 20 * n; ++k) F.add_edge((int)(rng() % n), n + (int)(rng() % n), 1, -(long long)(rng() % 1000000000));
        auto t0 = std::chrono::steady_clock::now();
        const auto r = F.flow_min_cost(s, t);
        std::cout << "[min_cost_flow] benchmark: 負のコスト（最大重みマッチング " << n << " × " << n << "、辺 " << 20 * n << "）  " << ms_since(t0) << " ms  (flow " << r.first
                  << ")\n";
    }
}

}  // namespace mcf_test


namespace mcf_demo {
[[maybe_unused]] static std::string to_s(const MinCostFlow::EdgeInfo& e) {  // demo_common の vector 版から見えるよう先に置く
    return "{" + std::to_string(e.from) + "→" + std::to_string(e.to) + ", cap=" + std::to_string(e.cap) + ", flow=" + std::to_string(e.flow) +
           ", cost=" + std::to_string(e.cost) + "}";
}
// ---- デモ用の小道具（見方: コード -> 返り値 // 意味） ----
template <class A, class B>
[[maybe_unused]] static std::string to_s(const std::pair<A, B>& p);
template <class T, size_t N>
[[maybe_unused]] static std::string to_s(const std::array<T, N>& a);
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v);
#ifdef GG_GRAPH_BASE
[[maybe_unused]] static std::string to_s(const Edge& e) { return "{to=" + std::to_string(e.to) + ", id=" + std::to_string(e.id) + ", w=" + std::to_string(e.w) + "}"; }
[[maybe_unused]] static std::string to_s(const RawEdge& e) { return "{" + std::to_string(e.from) + "-" + std::to_string(e.to) + ", w=" + std::to_string(e.w) + "}"; }
#endif
[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(char c) { return std::string("'") + c + "'"; }
[[maybe_unused]] static std::string to_s(const std::string& s) { return "\"" + s + "\""; }
[[maybe_unused]] static std::string to_s(const char* s) { return "\"" + std::string(s) + "\""; }
template <class T, std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool> && !std::is_same_v<T, char>, int> = 0>
[[maybe_unused]] static std::string to_s(T v) {
    if (std::is_same_v<T, long long> && (long long)v == (1LL << 60)) return "DIST_INF";
    return std::to_string(v);
}
template <class A, class B>
[[maybe_unused]] static std::string to_s(const std::pair<A, B>& p) { return "(" + to_s(p.first) + ", " + to_s(p.second) + ")"; }
template <class T, size_t N>
[[maybe_unused]] static std::string to_s(const std::array<T, N>& a) {
    std::string s = "(";
    for (size_t i = 0; i < N; ++i) s += (i ? ", " : "") + to_s(a[i]);
    return s + ")";
}
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
template <class T, class = decltype(std::declval<const T&>().to_vector())>
[[maybe_unused]] static std::string to_s(const T& v) { return to_s(v.to_vector()); }  // IntSpan など
[[maybe_unused]] static std::string pad(const std::string& s, size_t w) {  // 全角を 2 文字分として幅をそろえる
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
[[maybe_unused]] static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
[[maybe_unused]] static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
[[maybe_unused]] static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
[[maybe_unused]] static void stmt(const std::string& s) {
    const bool block = !s.empty() && s.back() == '}' &&
                       (s.compare(0, 4, "for ") == 0 || s.compare(0, 3, "if ") == 0 || s.compare(0, 6, "while ") == 0);
    code(block ? s : s + ";");
}
template <class V>
[[maybe_unused]] static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 36) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 56) + "  // " + comment;
    std::cout << line << "\n";
}
struct Row { const char* use; const char* ret; const char* what; };
template <size_t N>
[[maybe_unused]] static void cheat_sheet(const std::string& title, const Row (&rows)[N], size_t w1 = 46, size_t w2 = 22) {
    std::cout << "\n■ " << title << "\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, w1) << pad(ret, w2) << r.what << "\n";
    }
}
#define SHOW(expr, comment) show(#expr, (expr), comment)
#define RUN(...)             \
    stmt(#__VA_ARGS__);      \
    __VA_ARGS__

static void demo() {
    banner("最小費用流 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"MinCostFlow F(n);", "", "頂点 n 個"},
        {"F.add_edge(from, to, cap, cost)", "int", "辺（1 流すごとに cost。負も可）"},
        {"auto [f, c] = F.flow(s, t, limit);", "pair", "最大 limit まで、コスト最小で流す"},
        {"F.flow_min_cost(s, t)", "pair", "流量は自由でコストを最小に"},
        {"F.slope(s, t)", "vector<pair>", "流量ごとの最小コスト（折れ線）"},
        {"F.get_edge(id) / F.edges()", "EdgeInfo", "{from, to, cap, flow, cost}"},
    };
    cheat_sheet("API 一覧", rows, 36, 16);

    head("デモで使うグラフ: 0→1 (容量 2, コスト 1), 0→2 (1, 2), 1→2 (1, 1), 1→3 (1, 3), 2→3 (2, 1)  s = 0, t = 3");
    RUN(MinCostFlow F(4));
    RUN(F.add_edge(0, 1, 2, 1); F.add_edge(0, 2, 1, 2); F.add_edge(1, 2, 1, 1); F.add_edge(1, 3, 1, 3); F.add_edge(2, 3, 2, 1));

    api("【1】 F.slope(s, t)", "vector<pair<long long,long long>>", "流量 k ごとの最小コスト。折れ線の頂点の列（あいだは直線）", "O(F (N + M) log M)");
    RUN(MinCostFlow G = F);
    SHOW(G.slope(0, 3), "1 流すと 3、2 流すと 6、3 流すと 10");

    api("【2】 F.flow(s, t, limit)", "pair<long long,long long>", "最大 limit まで、コストが最小になるように流す。{流量, コスト}", "O(F (N + M) log M)");
    RUN(auto [f, c] = F.flow(0, 3, 2));
    SHOW(f, "流した量");
    SHOW(c, "コストの合計");
    SHOW(F.edges(), "各辺の流量");
    note("流量をちょうど K にしたいとき: flow(s, t, K) の流量が K 未満なら不可能。");

    api("【3】 最大重み（得点）の問題は、コストの符号を変える", "pair<long long,long long>", "負のコストも入れてよい（中で最初だけベルマンフォード）。flow_min_cost なら得点が増える間だけ流す", "O(F (N + M) log M) + 最初に O(NM)");
    code("// 例: 左 2 人・右 2 つの仕事、得点 a[i][j]（3 5 / 4 1）。全員を割り当てて得点の合計を最大に");
    RUN(std::vector<std::vector<long long>> a = {{3, 5}, {4, 1}});
    RUN(MinCostFlow M(6));
    RUN(for (int i = 0; i < 2; ++i) { M.add_edge(4, i, 1, 0); M.add_edge(2 + i, 5, 1, 0); for (int j = 0; j < 2; ++j) M.add_edge(i, 2 + j, 1, -a[i][j]); });
    SHOW(-M.flow(4, 5).second, "5 + 4");
    note("負のコストを避けたいなら「大きい定数 B - 得点」をコストにして、最後に B × 流量 を足し戻してもよい。");
    note("割り当て問題（全部の組に辺がある正方行列）なら 08_hungarian のほうが速い。");

    head("典型的な使い方");
    note("運ぶ量を決めてコスト最小       : F.flow(s, t, K)");
    note("最大重みマッチング（一部だけ）  : コストを -重み にして flow_min_cost(s, t)");
    note("K 個選んで得点最大（区間など）  : 頂点を一直線に並べ、容量 K の流れで表す（区間スケジューリングの形）");
}

}  // namespace mcf_demo

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
    mcf_test::random_test();
    mcf_test::benchmark();
    mcf_demo::demo();
    return 0;
}
