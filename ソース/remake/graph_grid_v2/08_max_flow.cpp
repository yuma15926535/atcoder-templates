// =============================================================================
//  最大流（最小カット・下限つきの流れ）  MaxFlow / LowerBoundFlow     --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 08 フロー
//  依存       : なし（単独で貼る）
//  頂点番号   : 0-indexed。容量は long long（「無限」は MaxFlow::CAP_INF = 4×10^18）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  高速化のポイント
//    - Dinic 法を非再帰で書いた（Visual Studio のスタック 1MB でも深いグラフで落ちない）
//    - BFS は t に着いたら打ち切り、DFS は t から逆向きにたどって 1 回でいくつもの道にまとめて流す
//      （AtCoder Library の maxflow と同じ方式。計測でも同じくらいの速さ）
//
//  機能一覧                                                                       計算量
//    MaxFlow F(n) / F.add_vertex()           頂点 n 個で作る / 頂点を足す
//    F.add_edge(from, to, cap)  -> int       辺を足す（辺番号を返す）
//    F.flow(s, t, limit = CAP_INF)           s → t に最大で limit 流す（流した量）     O(N^2 M)（実際はずっと速い）
//    F.min_cut(s) / F.cut_edges(s)           最小カットの s 側（1 / 0）/ カットの辺の番号     O(N + M)
//    F.get_edge(id) / F.edges()              {from, to, cap, flow}（流量の復元）
//    F.change_edge(id, cap, flow)            辺の容量・流量を書き換える
//    LowerBoundFlow L(n); L.add_edge(from, to, lo, hi)   流量を lo 以上 hi 以下にする辺
//      L.max_flow(s, t) / L.min_flow(s, t) / L.feasible()  最大 / 最小 / 循環が作れるか（解なしは -1、L.ok = false）
//      L.flow_of(id)                         各辺の流量
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mf_test::random_test() / mf_test::benchmark() / mf_demo::demo()
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
#ifndef GG_MAX_FLOW
#define GG_MAX_FLOW
// ---- 最大流（グラフ・グリッド第二世代）  Dinic 法（非再帰）・最小カット・下限つきの流れ ----
// 頂点は 0-indexed。容量は long long。add_edge が返す辺番号で、あとから流量を読める
struct MaxFlow {
    static constexpr long long CAP_INF = (long long)4e18;  // 「容量無限」の辺に使う値・flow の上限の既定値
    struct EdgeInfo { int from, to; long long cap, flow; };
    struct E { int to, rev; long long cap; };
    std::vector<std::vector<E>> g;
    std::vector<std::pair<int, int>> pos_;  // 辺番号 → (始点, g[始点] での位置)

    MaxFlow() {}
    explicit MaxFlow(int n) : g(n) {}
    int size() const { return (int)g.size(); }
    int add_vertex() { g.emplace_back(); return (int)g.size() - 1; }
    // 辺 from → to（容量 cap）。辺番号を返す
    int add_edge(int from, int to, long long cap) {
        const int id = (int)pos_.size();
        pos_.emplace_back(from, (int)g[from].size());
        const int rf = (int)g[from].size(), rt = (int)g[to].size() + (from == to ? 1 : 0);
        g[from].push_back({to, rt, cap});
        g[to].push_back({from, rf, 0});
        return id;
    }
    // s から t へ最大で limit 流す（流した量を返す）。続けて呼ぶと、今の流れに足していく。O(N^2 M)（実際はずっと速い）
    long long flow(int s, int t, long long limit = CAP_INF) {
        if (s == t) return 0;
        const int n = size();
        std::vector<int> level(n), it(n), q(n);
        struct Frame { int v; long long up, res; };
        std::vector<Frame> st;
        long long total = 0;
        while (total < limit) {
            // 残余グラフで s からの距離（レベル）。t に着いたら打ち切る
            std::fill(level.begin(), level.end(), -1);
            int qh = 0, qt = 0;
            level[s] = 0, q[qt++] = s;
            while (qh < qt && level[t] < 0) {
                const int v = q[qh++];
                for (const E& e : g[v])
                    if (e.cap > 0 && level[e.to] < 0) {
                        level[e.to] = level[v] + 1, q[qt++] = e.to;
                        if (e.to == t) break;
                    }
            }
            if (level[t] < 0) break;
            std::fill(it.begin(), it.end(), 0);
            // t から s へ、レベルが 1 ずつ減る辺を逆にたどり、1 回の探索でできるだけ多く流す（非再帰）
            while (total < limit) {
                long long ret = 0;
                bool back = false;
                st.clear();
                st.push_back({t, limit - total, 0});
                while (!st.empty()) {
                    const int k = (int)st.size() - 1;
                    const int v = st[k].v;
                    if (v == s) {
                        ret = st[k].up, st.pop_back(), back = true;
                        continue;
                    }
                    if (back) {  // 子（g[v][it[v]] の行き先）から ret 流せた
                        back = false;
                        if (ret > 0) {
                            E& e = g[v][it[v]];
                            e.cap += ret, g[e.to][e.rev].cap -= ret, st[k].res += ret;
                            if (st[k].res == st[k].up) {
                                ret = st[k].res, st.pop_back(), back = true;
                                continue;
                            }
                        }
                        ++it[v];
                    }
                    bool pushed = false;
                    for (; it[v] < (int)g[v].size(); ++it[v]) {
                        const E& e = g[v][it[v]];
                        const long long c = g[e.to][e.rev].cap;  // e.to → v の残り容量
                        if (level[v] <= level[e.to] || c == 0) continue;
                        st.push_back({e.to, std::min(st[k].up - st[k].res, c), 0});
                        pushed = true;
                        break;
                    }
                    if (pushed) continue;
                    level[v] = n;  // これ以上流せない
                    ret = st[k].res, st.pop_back(), back = true;
                }
                if (ret == 0) break;
                total += ret;
            }
        }
        return total;
    }
    // flow の後: s から残余グラフで行ける頂点（= 最小カットの s 側）なら 1
    std::vector<int> min_cut(int s) const {
        std::vector<int> vis(size(), 0), st{s};
        vis[s] = 1;
        while (!st.empty()) {
            const int v = st.back();
            st.pop_back();
            for (const E& e : g[v])
                if (e.cap > 0 && !vis[e.to]) vis[e.to] = 1, st.push_back(e.to);
        }
        return vis;
    }
    // flow の後: 最小カットの辺（s 側 → t 側の辺の番号。容量の和 = 最大流）
    std::vector<int> cut_edges(int s) const {
        const std::vector<int> side = min_cut(s);
        std::vector<int> r;
        for (int id = 0; id < (int)pos_.size(); ++id) {
            const EdgeInfo e = get_edge(id);
            if (side[e.from] && !side[e.to]) r.push_back(id);
        }
        return r;
    }
    EdgeInfo get_edge(int id) const {  // 辺 id の {始点, 終点, 容量, 流量}
        const E& e = g[pos_[id].first][pos_[id].second];
        const E& r = g[e.to][e.rev];
        return {pos_[id].first, e.to, e.cap + r.cap, r.cap};
    }
    std::vector<EdgeInfo> edges() const {
        std::vector<EdgeInfo> r;
        for (int id = 0; id < (int)pos_.size(); ++id) r.push_back(get_edge(id));
        return r;
    }
    int edge_count() const { return (int)pos_.size(); }
    // 辺 id の容量と流量を書き換える（0 ≤ new_flow ≤ new_cap）
    void change_edge(int id, long long new_cap, long long new_flow) {
        E& e = g[pos_[id].first][pos_[id].second];
        E& r = g[e.to][e.rev];
        e.cap = new_cap - new_flow, r.cap = new_flow;
    }
};
// 下限つきの流れ: 各辺の流量を lo 以上 hi 以下にする。max_flow / min_flow / feasible は何回呼んでもよい（毎回作り直す）
// 条件を満たせないときは -1 を返し ok = false（正味の流量が負になりうるグラフでは ok で判定する）
struct LowerBoundFlow {
    int n = 0;
    bool ok = false;               // 直前の計算で条件を満たせたか
    struct Ed { int from, to; long long lo, hi; };
    std::vector<Ed> es;
    std::vector<long long> flow_;  // 解いた後の各辺の流量
    LowerBoundFlow() {}
    explicit LowerBoundFlow(int n_) : n(n_) {}
    int add_edge(int from, int to, long long lo, long long hi) {  // lo ≤ 流量 ≤ hi。辺番号を返す
        es.push_back({from, to, lo, hi});
        return (int)es.size() - 1;
    }
    // 循環（すべての頂点で入る量 = 出る量）として条件を満たせるか
    bool feasible() { return solve(-1, -1, 0) >= 0; }
    // s → t の流量の最大（条件を満たせなければ -1）
    long long max_flow(int s, int t) { return solve(s, t, 1); }
    // s → t の流量の最小（条件を満たせなければ -1）
    long long min_flow(int s, int t) { return solve(s, t, 2); }
    long long flow_of(int id) const { return flow_[id]; }  // 解いた後の辺 id の流量

    long long solve(int s, int t, int mode) {
        ok = false;
        const int S = n, T = n + 1;
        MaxFlow F(n + 2);
        std::vector<long long> ex(n, 0);
        std::vector<int> ids(es.size());
        for (size_t i = 0; i < es.size(); ++i) {
            const Ed& e = es[i];
            if (e.lo > e.hi) return -1;
            ids[i] = F.add_edge(e.from, e.to, e.hi - e.lo);
            ex[e.to] += e.lo, ex[e.from] -= e.lo;
        }
        int ts = -1, st = -1;  // s と t の間を無限の容量で両向きにつなぐ（正味の流量が負の場合も含めて循環にする）
        if (mode != 0) ts = F.add_edge(t, s, MaxFlow::CAP_INF), st = F.add_edge(s, t, MaxFlow::CAP_INF);
        long long need = 0;
        for (int v = 0; v < n; ++v) {
            if (ex[v] > 0) F.add_edge(S, v, ex[v]), need += ex[v];
            if (ex[v] < 0) F.add_edge(v, T, -ex[v]);
        }
        if (F.flow(S, T) != need) return -1;
        long long ans = 0;
        if (mode != 0) {
            ans = F.get_edge(ts).flow - F.get_edge(st).flow;
            F.change_edge(ts, 0, 0), F.change_edge(st, 0, 0);
            if (mode == 1) ans += F.flow(s, t);
            else ans -= F.flow(t, s);
        }
        flow_.assign(es.size(), 0);
        for (size_t i = 0; i < es.size(); ++i) flow_[i] = es[i].lo + F.get_edge(ids[i]).flow;
        ok = true;
        return ans;
    }
};
#endif  // GG_MAX_FLOW
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "max_flow";
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

namespace mf_test {

// 愚直: 隣接行列で Edmonds-Karp
static long long brute_flow(int n, const std::vector<std::array<long long, 3>>& es, int s, int t) {
    std::vector<std::vector<long long>> c(n, std::vector<long long>(n, 0));
    for (auto& e : es) c[e[0]][e[1]] += e[2];
    long long total = 0;
    while (true) {
        std::vector<int> pre(n, -1);
        pre[s] = s;
        std::vector<int> q{s};
        for (size_t h = 0; h < q.size(); ++h)
            for (int u = 0; u < n; ++u)
                if (pre[u] < 0 && c[q[h]][u] > 0) pre[u] = q[h], q.push_back(u);
        if (pre[t] < 0) break;
        long long f = (long long)4e18;
        for (int v = t; v != s; v = pre[v]) f = std::min(f, c[pre[v]][v]);
        for (int v = t; v != s; v = pre[v]) c[pre[v]][v] -= f, c[v][pre[v]] += f;
        total += f;
    }
    return total;
}
static void check_flow(const MaxFlow& F, int n, int s, int t, long long value) {
    std::vector<long long> net(n, 0);
    for (const auto& e : F.edges()) {
        CHECK(0 <= e.flow && e.flow <= e.cap);
        net[e.from] -= e.flow, net[e.to] += e.flow;
    }
    for (int v = 0; v < n; ++v)
        if (v != s && v != t) CHECK(net[v] == 0);
    CHECK(net[t] == value && net[s] == -value);
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 2, 8), m = rnd(rng, 0, 16);
        MaxFlow F(n);
        std::vector<std::array<long long, 3>> es;
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            const long long c = (it % 5 == 0) ? (long long)(rng() % 1000000000000ULL) : rnd(rng, 0, 6);
            CHECK(F.add_edge(a, b, c) == i);
            es.push_back({a, b, c});
        }
        int s = rnd(rng, 0, n - 1), t = rnd(rng, 0, n - 1);
        if (s == t) t = (s + 1) % n;
        const long long ex = brute_flow(n, es, s, t);
        if (it % 3 == 0) {  // limit つきで 2 回に分けて流す
            const long long lim = rnd(rng, 0, 5);
            const long long a = F.flow(s, t, lim);
            CHECK(a == std::min(ex, lim));
            const long long b = F.flow(s, t);
            CHECK(a + b == ex);
        } else {
            CHECK(F.flow(s, t) == ex);
        }
        check_flow(F, n, s, t, ex);
        // 最小カット
        const std::vector<int> side = F.min_cut(s);
        CHECK(side[s] == 1 && side[t] == 0);
        long long cut = 0;
        for (int id : F.cut_edges(s)) {
            const auto e = F.get_edge(id);
            CHECK(side[e.from] && !side[e.to] && e.flow == e.cap);
            cut += e.cap;
        }
        CHECK(cut == ex);
        for (int id = 0; id < m; ++id) {
            const auto e = F.get_edge(id);
            CHECK(e.from == es[id][0] && e.to == es[id][1] && e.cap == es[id][2]);
        }
        // change_edge: 全部の流れを消してから流し直しても同じ
        for (int id = 0; id < m; ++id) F.change_edge(id, es[id][2], 0);
        CHECK(F.flow(s, t) == ex);
        CHECK(F.flow(s, s) == 0);
    }
    // 下限つき: 全部の流し方を数え上げて比べる
    for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 2, 4), m = rnd(rng, 1, 5);
        LowerBoundFlow L(n);
        std::vector<std::array<int, 4>> es;
        for (int i = 0; i < m; ++i) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1), lo = rnd(rng, 0, 2), hi = lo + rnd(rng, 0, 2);
            L.add_edge(a, b, lo, hi);
            es.push_back({a, b, lo, hi});
        }
        int s = rnd(rng, 0, n - 1), t = rnd(rng, 0, n - 1);
        if (s == t) t = (s + 1) % n;
        bool feas = false, feas_c = false;
        long long mx = -(1LL << 60), mn = 1LL << 60;
        std::vector<int> f(m);
        std::function<void(int)> rec = [&](int i) {
            if (i == m) {
                std::vector<int> net(n, 0);
                for (int k = 0; k < m; ++k) net[es[k][0]] -= f[k], net[es[k][1]] += f[k];
                bool ok = true, okc = true;
                for (int v = 0; v < n; ++v) {
                    if (net[v] != 0) okc = false;
                    if (v != s && v != t && net[v] != 0) ok = false;
                }
                if (okc) feas_c = true;
                if (ok) feas = true, mx = std::max(mx, (long long)net[t]), mn = std::min(mn, (long long)net[t]);
                return;
            }
            for (int x = es[i][2]; x <= es[i][3]; ++x) f[i] = x, rec(i + 1);
        };
        rec(0);
        const auto check_sol = [&](bool circ) {
            std::vector<long long> net(n, 0);
            for (int k = 0; k < m; ++k) {
                CHECK(es[k][2] <= L.flow_of(k) && L.flow_of(k) <= es[k][3]);
                net[es[k][0]] -= L.flow_of(k), net[es[k][1]] += L.flow_of(k);
            }
            for (int v = 0; v < n; ++v)
                if (circ || (v != s && v != t)) CHECK(net[v] == 0);
            return net[t];
        };
        CHECK(L.feasible() == feas_c && L.ok == feas_c);
        if (feas_c) check_sol(true);
        const long long a = L.max_flow(s, t);
        CHECK(L.ok == feas);
        if (feas) CHECK(a == mx && check_sol(false) == mx);
        else CHECK(a == -1);
        const long long b = L.min_flow(s, t);
        CHECK(L.ok == feas);
        if (feas) CHECK(b == mn && check_sol(false) == mn);
    }
    std::cout << "[max_flow] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[max_flow] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    {  // 二部マッチング（左右 5 万、辺 20 万）
        const int L = 50000, R = 50000, M = 200000;
        MaxFlow F(L + R + 2);
        const int s = L + R, t = s + 1;
        for (int i = 0; i < L; ++i) F.add_edge(s, i, 1);
        for (int j = 0; j < R; ++j) F.add_edge(L + j, t, 1);
        for (int k = 0; k < M; ++k) F.add_edge((int)(rng() % L), L + (int)(rng() % R), 1);
        auto t0 = std::chrono::steady_clock::now();
        const long long f = F.flow(s, t);
        std::cout << "[max_flow] benchmark: 二部マッチング 左右 " << L << "・辺 " << M << "  " << ms_since(t0) << " ms  (flow " << f << ")\n";
    }
    {  // ランダムなグラフ（N = 5000, M = 50000, 容量ランダム）
        const int n = 5000, m = 50000;
        MaxFlow F(n);
        for (int k = 0; k < m; ++k) F.add_edge((int)(rng() % n), (int)(rng() % n), (long long)(rng() % 1000000));
        auto t0 = std::chrono::steady_clock::now();
        const long long f = F.flow(0, n - 1);
        std::cout << "[max_flow] benchmark: ランダム N = " << n << ", M = " << m << "  " << ms_since(t0) << " ms  (flow " << f << ")\n";
    }
    {  // グリッド（燃やす埋めるの形。300 × 300 マス、隣どうしの容量ランダム）
        const int H = 300, W = 300, n = H * W;
        MaxFlow F(n + 2);
        const int s = n, t = n + 1;
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) {
                const int v = i * W + j;
                F.add_edge(s, v, (long long)(rng() % 100)), F.add_edge(v, t, (long long)(rng() % 100));
                if (i + 1 < H) F.add_edge(v, v + W, (long long)(rng() % 50)), F.add_edge(v + W, v, (long long)(rng() % 50));
                if (j + 1 < W) F.add_edge(v, v + 1, (long long)(rng() % 50)), F.add_edge(v + 1, v, (long long)(rng() % 50));
            }
        auto t0 = std::chrono::steady_clock::now();
        const long long f = F.flow(s, t);
        std::cout << "[max_flow] benchmark: グリッド " << H << " × " << W << "（燃やす埋めるの形）  " << ms_since(t0) << " ms  (flow " << f << ")\n";
    }
}

}  // namespace mf_test


namespace mf_demo {
[[maybe_unused]] static std::string to_s(const MaxFlow::EdgeInfo& e) {  // demo_common の vector 版から見えるよう先に置く
    return "{" + std::to_string(e.from) + "→" + std::to_string(e.to) + ", cap=" + std::to_string(e.cap) + ", flow=" + std::to_string(e.flow) + "}";
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
    banner("最大流 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"MaxFlow F(n);", "", "頂点 n 個"},
        {"int id = F.add_edge(from, to, cap);", "int", "辺を足す（辺番号）"},
        {"F.flow(s, t) / F.flow(s, t, limit)", "long long", "最大流 / 最大で limit まで"},
        {"F.min_cut(s) / F.cut_edges(s)", "vector<int>", "最小カットの s 側 / カットの辺"},
        {"F.get_edge(id) / F.edges()", "EdgeInfo", "{from, to, cap, flow}"},
        {"F.change_edge(id, cap, flow)", "", "辺を書き換える"},
        {"LowerBoundFlow L(n); L.add_edge(u, v, lo, hi);", "", "下限つき（lo ≤ 流量 ≤ hi）"},
        {"L.max_flow(s, t) / L.min_flow(s, t) / L.feasible()", "long long / bool", "解なしは -1（L.ok = false）"},
    };
    cheat_sheet("API 一覧", rows, 50, 18);

    head("デモで使うグラフ: 0→1 (3), 0→2 (2), 1→2 (1), 1→3 (2), 2→3 (3)  s = 0, t = 3");
    RUN(MaxFlow F(4));
    RUN(F.add_edge(0, 1, 3); F.add_edge(0, 2, 2); F.add_edge(1, 2, 1); F.add_edge(1, 3, 2); F.add_edge(2, 3, 3));

    api("【1】 F.flow(s, t)", "long long", "s から t へ流せる量の最大", "O(N^2 M)（実際はずっと速い）");
    SHOW(F.flow(0, 3), "");
    SHOW(F.edges(), "各辺の流量（流し方の復元）");

    api("【2】 F.min_cut(s) / F.cut_edges(s)", "vector<int>", "最小カット（容量の和が最小になる「s 側と t 側の分け方」）。容量の和 = 最大流", "O(N + M)");
    SHOW(F.min_cut(0), "s 側なら 1");
    SHOW(F.cut_edges(0), "s 側 → t 側の辺の番号");

    api("【3】 二部マッチング（最大流で解く形）", "long long", "左 L 人・右 R 個。s → 左（容量 1）、左 → 右（容量 1）、右 → t（容量 1）", "O(M √N)");
    RUN(int L = 3, R = 3);
    RUN(MaxFlow B(L + R + 2));
    RUN(int s = L + R, t = s + 1);
    RUN(for (int i = 0; i < L; ++i) B.add_edge(s, i, 1));
    RUN(for (int j = 0; j < R; ++j) B.add_edge(L + j, t, 1));
    RUN(std::vector<int> ids = {B.add_edge(0, L + 0, 1), B.add_edge(0, L + 1, 1), B.add_edge(1, L + 0, 1), B.add_edge(2, L + 0, 1)});
    SHOW(B.flow(s, t), "最大マッチング");
    note("組を出力するなら、左 → 右の辺のうち get_edge(id).flow == 1 のもの。08_bipartite_matching のほうが速くて便利。");

    api("【4】 燃やす埋める（最小カットで 2 択の最適化）", "long long", "「v を A にすると a、B にすると b 失う」「u が A で v が B なら c 失う」の最小", "最大流と同じ");
    code("// 頂点 v: s 側 = A、t 側 = B。s → v（容量 b）: v を B にすると切れる。v → t（容量 a）: v を A にすると切れる");
    code("// 「u が A かつ v が B なら c」: u → v（容量 c）。「得る」は先に全部足しておき、失う量として入れる");
    RUN(MaxFlow P(2 + 2));
    RUN(int S = 2, T = 3);
    RUN(P.add_edge(S, 0, 5); P.add_edge(0, T, 1); P.add_edge(S, 1, 2); P.add_edge(1, T, 4); P.add_edge(0, 1, 3));
    SHOW(P.flow(S, T), "失う量の最小");
    SHOW(P.min_cut(S), "0, 1 がそれぞれ A（s 側 = 1）か B か");

    api("【5】 LowerBoundFlow（流量に下限がある）", "long long", "各辺の流量を lo 以上 hi 以下にする。解なしは -1（L.ok = false）", "最大流 2 回分");
    RUN(LowerBoundFlow LB(3));
    RUN(LB.add_edge(0, 1, 2, 5); LB.add_edge(1, 2, 0, 3); LB.add_edge(0, 2, 1, 4));
    SHOW(LB.max_flow(0, 2), "1 → 2 が 3 までなので 3 + 4");
    SHOW(LB.min_flow(0, 2), "0 → 1 に 2 以上、0 → 2 に 1 以上流す必要がある");
    SHOW(LB.flow_of(0), "最後に解いた流し方での辺 0 の流量");
    RUN(LowerBoundFlow LC(2));
    RUN(LC.add_edge(0, 1, 3, 3));
    SHOW(LC.feasible(), "0 → 1 に 3 流すと戻る道がないので循環にならない");

    head("典型的な使い方");
    note("最大マッチング             : 08_bipartite_matching（速い。最小点被覆・最大独立集合も）");
    note("辺・頂点を何本消せば分断   : 最小カット（頂点は v_in → v_out の容量 1 の辺に分ける）");
    note("2 択の最適化（燃やす埋める）: 最小カット");
    note("コストも最小にしたい       : 08_min_cost_flow");
}

}  // namespace mf_demo

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
    mf_test::random_test();
    mf_test::benchmark();
    mf_demo::demo();
    return 0;
}
