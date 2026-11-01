// =============================================================================
//  bit DP（巡回セールスマン・ハミルトン路・割り当ての数え上げ・集合をいくつかの組に分ける）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  集合 {0, …, n−1} の部分集合を n ビットの数 mask で表し（i ビット目が 1 ⇔ i を含む）、mask を状態にする DP
//    例: ABC180 E（Traveling Salesman among Aerial Cities）… tsp_cycle(dist, 0)
//        EDPC O（Matching）… count_assignments<mint>(ok)
//        ABC187 F（Close Group）… min_partition_count(n, clique_table(adj))
//
//  機能一覧（dist は n × n の距離。辺が無いところは bitdp_inf<T>() 以上の値）          計算量
//   巡回セールスマン・ハミルトン路
//    bitdp_inf<T>()                        「無い・届かない」を表す大きい値（T の最大値 / 4）
//    tsp_table(dist, start)                dp[mask][v]：start から出て mask の頂点をちょうど 1 回ずつ通り、
//                                          今 v にいる最小コスト（届かないなら bitdp_inf。start = −1 ならどこから出てもよい）  O(2^n n^2)
//    tsp_cycle(dist, start = 0)            全頂点を 1 回ずつ通って start に戻る最小コスト（無理なら bitdp_inf）  O(2^n n^2)
//    tsp_path(dist, start = −1, goal = −1) 全頂点を 1 回ずつ通る道の最小コスト（−1 はどこでもよい）
//    tsp_cycle_restore(dist, start = 0)    {コスト, 回る順番（start, …, start の n + 1 個）}。無理なら {bitdp_inf, {}}
//    tsp_path_restore(dist, start, goal)   {コスト, 通る順番（n 個）}
//   割り当ての数え上げ
//    count_assignments<T>(ok)              n 人を m 個の仕事（n ≤ m）に 1 人 1 つずつ、どの 2 人も同じ仕事にならないように
//                                          割り当てる方法の数（ok[i][j] が 0 でないとき i は j をできる）。T は mint など  O(2^m m)
//   集合を組に分ける（どの要素もちょうど 1 つの組に入る。組の順番は区別しない）
//    subset_sum_table(a)                   s[mask] = mask に入る a_i の和                                         O(2^n)
//    clique_table(adj) / independent_table(adj)   ok[mask] = mask がクリーク / 独立集合か（adj は n × n の 0 / 1）  O(2^n)
//    min_partition_cost(n, cost)           組 g のコスト cost[g]（g ≠ 0）の和の最小。無理（どの分け方も bitdp_inf を含む）なら bitdp_inf  O(3^n)
//    min_partition_restore(n, cost)        {最小コスト, 組の一覧（mask）}
//    min_partition_cost_k(n, cost, K)      ちょうど K 組に分けるときの最小（無理なら bitdp_inf）               O(K 3^n)
//    min_partition_count(n, ok)            ok[g] が true の組だけで分けるときの最小の組数（無理なら −1）       O(3^n)
//    count_partitions<T>(n, ok)            ok[g] が true の組だけでの分け方の数                                  O(3^n)
//
//  注意
//    n の目安：tsp は n ≤ 20（dp 表が 2^n × n 個。n = 20 の long long で 約 170 MB）、count_assignments は m ≤ 22、
//    分割（3^n）は n ≤ 18（3^18 / 2 ≈ 2 × 10^8 回）
//    コストは T（long long など）。足し算は bitdp_inf を超えたら打ち切るので、bitdp_inf 未満の値どうしならあふれない
//    最大化したいときはコストを −1 倍して最小化し、答えを −1 倍する（tsp・分割とも）
//    n = 1 の tsp_cycle は 0（動かない）
//  速さの目安: tsp_cycle（n = 17）約 0.07 秒、min_partition_count（n = 18）約 0.2 秒、count_assignments（m = 21）約 0.03 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    bitdp_test::random_test()  … 順列・分け方を全部調べる愚直解との突き合わせ
//    bitdp_test::benchmark()    … 速度計測
//    bitdp_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace bitdp_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
template <class T>
void check_square(const std::vector<std::vector<T>>& d, const char* f) {
    const size_t n = d.size();
    for (const auto& r : d)
        if (r.size() != n) fail(f, "dist は n × n にする");
    if (n > 25) fail(f, "頂点が多すぎる（2^n × n の表を作る）");
}
}  // namespace bitdp_internal

template <class T>
constexpr T bitdp_inf() {
    return std::numeric_limits<T>::max() / 4;
}

// ---- 巡回セールスマン・ハミルトン路 ----
template <class T>
std::vector<std::vector<T>> tsp_table(const std::vector<std::vector<T>>& dist, int start) {
    bitdp_internal::check_square(dist, "tsp_table");
    const int n = (int)dist.size();
    if (start < -1 || start >= n) bitdp_internal::fail("tsp_table", "start は −1 か 0〜n−1");
    const T INF = bitdp_inf<T>();
    std::vector<std::vector<T>> dp((size_t)1 << n, std::vector<T>(n, INF));
    if (n == 0) return dp;
    if (start >= 0) dp[(size_t)1 << start][start] = T(0);
    else
        for (int v = 0; v < n; ++v) dp[(size_t)1 << v][v] = T(0);
    const size_t full = (size_t)1 << n;
    for (size_t mask = 1; mask < full; ++mask) {
        if (start >= 0 && !(mask >> start & 1)) continue;
        const auto& cur = dp[mask];
        for (int v = 0; v < n; ++v) {
            if (!(mask >> v & 1) || !(cur[v] < INF)) continue;
            const T base = cur[v];
            const auto& dv = dist[v];
            for (int u = 0; u < n; ++u) {
                if (mask >> u & 1) continue;
                if (!(dv[u] < INF)) continue;
                T nv = base + dv[u];
                if (!(nv < INF)) continue;
                T& to = dp[mask | ((size_t)1 << u)][u];
                if (nv < to) to = nv;
            }
        }
    }
    return dp;
}
template <class T>
T tsp_cycle(const std::vector<std::vector<T>>& dist, int start = 0) {
    const int n = (int)dist.size();
    if (n == 0) return T(0);
    if (start < 0 || start >= n) bitdp_internal::fail("tsp_cycle", "start は 0〜n−1");
    if (n == 1) return T(0);
    const auto dp = tsp_table(dist, start);
    const T INF = bitdp_inf<T>();
    T best = INF;
    const auto& last = dp[((size_t)1 << n) - 1];
    for (int v = 0; v < n; ++v) {
        if (v == start || !(last[v] < INF) || !(dist[v][start] < INF)) continue;
        const T c = last[v] + dist[v][start];
        if (c < best) best = c;
    }
    return best;
}
template <class T>
T tsp_path(const std::vector<std::vector<T>>& dist, int start = -1, int goal = -1) {
    const int n = (int)dist.size();
    if (n == 0) return T(0);
    if (goal < -1 || goal >= n) bitdp_internal::fail("tsp_path", "goal は −1 か 0〜n−1");
    const auto dp = tsp_table(dist, start);
    const auto& last = dp[((size_t)1 << n) - 1];
    if (goal >= 0) return last[goal];
    return *std::min_element(last.begin(), last.end());
}
namespace bitdp_internal {
// dp 表から (mask, v) までの順番を戻す
template <class T>
std::vector<int> tsp_back(const std::vector<std::vector<T>>& dp, const std::vector<std::vector<T>>& dist, size_t mask, int v) {
    const int n = (int)dist.size();
    std::vector<int> order{v};
    while ((mask & (mask - 1)) != 0) {  // 2 つ以上
        const size_t pm = mask ^ ((size_t)1 << v);
        int prv = -1;
        for (int u = 0; u < n && prv < 0; ++u)
            if ((pm >> u & 1) && dp[pm][u] < bitdp_inf<T>() && dist[u][v] < bitdp_inf<T>() && dp[pm][u] + dist[u][v] == dp[mask][v]) prv = u;
        if (prv < 0) fail("tsp_restore", "復元できない（内部エラー）");
        order.push_back(prv);
        mask = pm, v = prv;
    }
    std::reverse(order.begin(), order.end());
    return order;
}
}  // namespace bitdp_internal
template <class T>
std::pair<T, std::vector<int>> tsp_cycle_restore(const std::vector<std::vector<T>>& dist, int start = 0) {
    const int n = (int)dist.size();
    if (n == 0) return {T(0), {}};
    if (start < 0 || start >= n) bitdp_internal::fail("tsp_cycle_restore", "start は 0〜n−1");
    if (n == 1) return {T(0), {start, start}};
    const auto dp = tsp_table(dist, start);
    const T INF = bitdp_inf<T>();
    const size_t full = ((size_t)1 << n) - 1;
    T best = INF;
    int bv = -1;
    for (int v = 0; v < n; ++v) {
        if (v == start || !(dp[full][v] < INF) || !(dist[v][start] < INF)) continue;
        const T c = dp[full][v] + dist[v][start];
        if (c < best) best = c, bv = v;
    }
    if (bv < 0) return {INF, {}};
    auto order = bitdp_internal::tsp_back(dp, dist, full, bv);
    order.push_back(start);
    return {best, order};
}
template <class T>
std::pair<T, std::vector<int>> tsp_path_restore(const std::vector<std::vector<T>>& dist, int start = -1, int goal = -1) {
    const int n = (int)dist.size();
    if (n == 0) return {T(0), {}};
    if (goal < -1 || goal >= n) bitdp_internal::fail("tsp_path_restore", "goal は −1 か 0〜n−1");
    const auto dp = tsp_table(dist, start);
    const T INF = bitdp_inf<T>();
    const size_t full = ((size_t)1 << n) - 1;
    int bv = goal;
    if (goal < 0) bv = (int)(std::min_element(dp[full].begin(), dp[full].end()) - dp[full].begin());
    if (!(dp[full][bv] < INF)) return {INF, {}};
    return {dp[full][bv], bitdp_internal::tsp_back(dp, dist, full, bv)};
}

// ---- 割り当ての数え上げ ----
template <class T, class Mat>
T count_assignments(const Mat& ok) {
    const int n = (int)ok.size();
    if (n == 0) return T(1);
    const int m = (int)ok[0].size();
    for (const auto& r : ok)
        if ((int)r.size() != m) bitdp_internal::fail("count_assignments", "ok の行の長さがそろっていない");
    if (n > m) return T(0);
    if (m > 26) bitdp_internal::fail("count_assignments", "仕事が多すぎる（2^m の表を作る）");
    std::vector<unsigned long long> okm(n, 0);  // 人 i ができる仕事の集合
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            if (ok[i][j]) okm[i] |= 1ULL << j;
    std::vector<T> dp((size_t)1 << m, T(0));
    dp[0] = T(1);
    T total = T(0);
    for (size_t mask = 0; mask < ((size_t)1 << m); ++mask) {
        const int i = __builtin_popcountll((unsigned long long)mask);  // 次に割り当てる人
        if (i == n) {
            total += dp[mask];
            continue;
        }
        if (i > n) continue;
        const T cur = dp[mask];
        for (unsigned long long av = okm[i] & ~(unsigned long long)mask; av; av &= av - 1)
            dp[mask | ((size_t)1 << __builtin_ctzll(av))] += cur;
    }
    return total;
}

// ---- 集合を組に分ける ----
template <class T>
std::vector<T> subset_sum_table(const std::vector<T>& a) {
    const int n = (int)a.size();
    if (n > 30) bitdp_internal::fail("subset_sum_table", "要素が多すぎる");
    std::vector<T> s((size_t)1 << n, T(0));
    for (int i = 0; i < n; ++i)
        for (size_t m = 0; m < ((size_t)1 << i); ++m) s[m | ((size_t)1 << i)] = s[m] + a[i];
    return s;
}
namespace bitdp_internal {
template <class Mat>
std::vector<char> clique_like(const Mat& adj, bool want_edge, const char* f) {
    const int n = (int)adj.size();
    for (const auto& r : adj)
        if ((int)r.size() != n) fail(f, "adj は n × n にする");
    if (n > 30) fail(f, "頂点が多すぎる");
    std::vector<unsigned> nb(n, 0);  // i とつながっている（独立集合なら、つながっていない）頂点
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && ((adj[i][j] != 0) == want_edge)) nb[i] |= 1u << j;
    std::vector<char> ok((size_t)1 << n, 0);
    ok[0] = 1;
    for (size_t m = 1; m < ((size_t)1 << n); ++m) {
        const int v = __builtin_ctzll((unsigned long long)m);
        const size_t rest = m & (m - 1);
        ok[m] = ok[rest] && ((rest & ~(size_t)nb[v]) == 0);
    }
    return ok;
}
}  // namespace bitdp_internal
template <class Mat>
std::vector<char> clique_table(const Mat& adj) {
    return bitdp_internal::clique_like(adj, true, "clique_table");
}
template <class Mat>
std::vector<char> independent_table(const Mat& adj) {
    return bitdp_internal::clique_like(adj, false, "independent_table");
}
namespace bitdp_internal {
template <class V>
int check_table(int n, const V& t, const char* f) {
    if (n < 0 || n > 26) fail(f, "n は 0〜26");
    if (t.size() != ((size_t)1 << n)) fail(f, "表の長さは 2^n にする");
    return n;
}
}  // namespace bitdp_internal
template <class T>
T min_partition_cost(int n, const std::vector<T>& cost) {
    bitdp_internal::check_table(n, cost, "min_partition_cost");
    const T INF = bitdp_inf<T>();
    const size_t full = ((size_t)1 << n) - 1;
    std::vector<T> dp(full + 1, INF);
    dp[0] = T(0);
    for (size_t S = 1; S <= full; ++S) {
        const size_t low = S & (~S + 1), rest = S ^ low;
        T best = INF;
        for (size_t sub = rest;; sub = (sub - 1) & rest) {  // 一番小さい要素 low を含む組 g = sub | low
            const size_t g = sub | low;
            if (cost[g] < INF && dp[S ^ g] < INF) {
                const T c = dp[S ^ g] + cost[g];
                if (c < best) best = c;
            }
            if (sub == 0) break;
        }
        dp[S] = best < INF ? best : INF;
    }
    return dp[full];
}
template <class T>
std::pair<T, std::vector<unsigned long long>> min_partition_restore(int n, const std::vector<T>& cost) {
    bitdp_internal::check_table(n, cost, "min_partition_restore");
    const T INF = bitdp_inf<T>();
    const size_t full = ((size_t)1 << n) - 1;
    std::vector<T> dp(full + 1, INF);
    std::vector<size_t> pick(full + 1, 0);
    dp[0] = T(0);
    for (size_t S = 1; S <= full; ++S) {
        const size_t low = S & (~S + 1), rest = S ^ low;
        for (size_t sub = rest;; sub = (sub - 1) & rest) {
            const size_t g = sub | low;
            if (cost[g] < INF && dp[S ^ g] < INF) {
                const T c = dp[S ^ g] + cost[g];
                if (c < dp[S]) dp[S] = c, pick[S] = g;
            }
            if (sub == 0) break;
        }
    }
    if (!(dp[full] < INF)) return {INF, {}};
    std::vector<unsigned long long> groups;
    for (size_t S = full; S != 0; S ^= pick[S]) groups.push_back(pick[S]);
    return {dp[full], groups};
}
template <class T>
T min_partition_cost_k(int n, const std::vector<T>& cost, int K) {
    bitdp_internal::check_table(n, cost, "min_partition_cost_k");
    const T INF = bitdp_inf<T>();
    const size_t full = ((size_t)1 << n) - 1;
    if (K < 0) return INF;
    if (n == 0) return K == 0 ? T(0) : INF;
    if (K == 0 || K > n) return INF;
    std::vector<T> dp(full + 1, INF), nd(full + 1);
    dp[0] = T(0);
    for (int k = 1; k <= K; ++k) {
        std::fill(nd.begin(), nd.end(), INF);
        for (size_t S = 1; S <= full; ++S) {
            const size_t low = S & (~S + 1), rest = S ^ low;
            T best = INF;
            for (size_t sub = rest;; sub = (sub - 1) & rest) {
                const size_t g = sub | low;
                if (cost[g] < INF && dp[S ^ g] < INF) {
                    const T c = dp[S ^ g] + cost[g];
                    if (c < best) best = c;
                }
                if (sub == 0) break;
            }
            nd[S] = best;
        }
        std::swap(dp, nd);
    }
    return dp[full];
}
template <class V>
int min_partition_count(int n, const V& ok) {
    bitdp_internal::check_table(n, ok, "min_partition_count");
    const size_t full = ((size_t)1 << n) - 1;
    const int BIG = 1 << 29;
    std::vector<int> dp(full + 1, BIG);
    dp[0] = 0;
    for (size_t S = 1; S <= full; ++S) {
        if (ok[S]) {
            dp[S] = 1;
            continue;
        }
        const size_t low = S & (~S + 1), rest = S ^ low;
        int best = BIG;
        for (size_t sub = rest;; sub = (sub - 1) & rest) {
            const size_t g = sub | low;
            if (ok[g] && dp[S ^ g] + 1 < best) best = dp[S ^ g] + 1;
            if (sub == 0) break;
        }
        dp[S] = best;
    }
    return dp[full] >= BIG ? -1 : dp[full];
}
template <class T, class V>
T count_partitions(int n, const V& ok) {
    bitdp_internal::check_table(n, ok, "count_partitions");
    const size_t full = ((size_t)1 << n) - 1;
    std::vector<T> dp(full + 1, T(0));
    dp[0] = T(1);
    for (size_t S = 1; S <= full; ++S) {
        const size_t low = S & (~S + 1), rest = S ^ low;
        T c = T(0);
        for (size_t sub = rest;; sub = (sub - 1) & rest) {
            const size_t g = sub | low;
            if (ok[g]) c += dp[S ^ g];
            if (sub == 0) break;
        }
        dp[S] = c;
    }
    return dp[full];
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace bitdp_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define BD_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[dp_bitdp] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
static const ll INF = bitdp_inf<ll>();

static void test_tsp(std::mt19937_64& rng, long long& cases) {
    const int n = (int)rnd(rng, 1, 7);
    std::vector<std::vector<ll>> d(n, std::vector<ll>(n));
    const int none = (int)rnd(rng, 0, 3);  // 辺が無い割合
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) d[i][j] = (int)rnd(rng, 0, 9) < none ? (rnd(rng, 0, 1) ? INF : INF + 5) : rnd(rng, 0, 20);
    std::vector<int> p(n);
    auto cost_of = [&](const std::vector<int>& ord) {  // 順番どおりに通るコスト（無理なら INF）
        ll c = 0;
        for (size_t i = 0; i + 1 < ord.size(); ++i) {
            if (d[ord[i]][ord[i + 1]] >= INF) return INF;
            c += d[ord[i]][ord[i + 1]];
        }
        return c;
    };
    for (int s = 0; s < n; ++s) {
        // 巡回
        ll best = n == 1 ? 0 : INF;
        if (n > 1) {
            std::vector<int> rest;
            for (int v = 0; v < n; ++v)
                if (v != s) rest.push_back(v);
            std::sort(rest.begin(), rest.end());
            do {
                std::vector<int> ord{s};
                ord.insert(ord.end(), rest.begin(), rest.end());
                ord.push_back(s);
                best = std::min(best, cost_of(ord));
            } while (std::next_permutation(rest.begin(), rest.end()));
        }
        BD_CHECK(tsp_cycle(d, s) == best);
        const auto rc = tsp_cycle_restore(d, s);
        BD_CHECK(rc.first == best);
        if (best < INF) {
            BD_CHECK((int)rc.second.size() == n + 1 && rc.second.front() == s && rc.second.back() == s);
            std::vector<int> seen(rc.second.begin(), rc.second.end() - 1);
            std::sort(seen.begin(), seen.end());
            for (int v = 0; v < n; ++v) BD_CHECK(seen[v] == v);
            BD_CHECK(n == 1 || cost_of(rc.second) == best);  // n = 1 は動かないので 0
        } else {
            BD_CHECK(rc.second.empty());
        }
        ++cases;
    }
    // 道（start・goal は決める・決めない）
    for (int s = -1; s < n; ++s)
        for (int g = -1; g < n; ++g) {
            ll best = INF;
            for (int v = 0; v < n; ++v) p[v] = v;
            do {
                if (s >= 0 && p[0] != s) continue;
                if (g >= 0 && p[n - 1] != g) continue;
                best = std::min(best, cost_of(p));
            } while (std::next_permutation(p.begin(), p.end()));
            BD_CHECK(tsp_path(d, s, g) == best);
            const auto rp = tsp_path_restore(d, s, g);
            BD_CHECK(rp.first == best);
            if (best < INF) {
                BD_CHECK((int)rp.second.size() == n && cost_of(rp.second) == best);
                if (s >= 0) BD_CHECK(rp.second.front() == s);
                if (g >= 0) BD_CHECK(rp.second.back() == g);
                std::vector<int> seen = rp.second;
                std::sort(seen.begin(), seen.end());
                for (int v = 0; v < n; ++v) BD_CHECK(seen[v] == v);
            }
            ++cases;
        }
}

static void test_assign(std::mt19937_64& rng, long long& cases) {
    const int n = (int)rnd(rng, 0, 6), m = (int)rnd(rng, 0, 7);
    std::vector<std::vector<int>> ok(n, std::vector<int>(m));
    for (auto& r : ok)
        for (auto& x : r) x = (int)rnd(rng, 0, 2) > 0;
    // 愚直：仕事を全部の順列で並べ、先頭 n 個を人に（同じ割り当てが (m − n)! 回出る）
    ll c = 0;
    if (n <= m) {
        std::vector<int> p(m);
        for (int j = 0; j < m; ++j) p[j] = j;
        ll dup = 1;
        for (int k = 1; k <= m - n; ++k) dup *= k;
        do {
            bool g = true;
            for (int i = 0; i < n && g; ++i) g = ok[i][p[i]];
            c += g;
        } while (std::next_permutation(p.begin(), p.end()));
        c /= dup;
    }
    if (n == 0) c = 1;
    BD_CHECK(count_assignments<ll>(ok) == c);
    ++cases;
}

// 全部の分け方を列挙（各組の mask の一覧）
static void all_partitions(int n, int i, std::vector<unsigned>& cur, std::vector<std::vector<unsigned>>& out) {
    if (i == n) {
        out.push_back(cur);
        return;
    }
    for (size_t k = 0; k < cur.size(); ++k) {
        cur[k] |= 1u << i;
        all_partitions(n, i + 1, cur, out);
        cur[k] ^= 1u << i;
    }
    cur.push_back(1u << i);
    all_partitions(n, i + 1, cur, out);
    cur.pop_back();
}

static void test_partition(std::mt19937_64& rng, long long& cases) {
    const int n = (int)rnd(rng, 0, 7);
    std::vector<ll> cost((size_t)1 << n);
    for (auto& c : cost) c = rnd(rng, 0, 4) == 0 ? INF : rnd(rng, -10, 30);
    std::vector<std::vector<int>> adj(n, std::vector<int>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) adj[i][j] = adj[j][i] = (int)rnd(rng, 0, 1);
    const auto cl = clique_table(adj), in = independent_table(adj);
    for (size_t m = 0; m < cost.size(); ++m) {
        bool c = true, d = true;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if ((m >> i & 1) && (m >> j & 1)) c = c && adj[i][j], d = d && !adj[i][j];
        BD_CHECK(cl[m] == c && in[m] == d);
    }
    std::vector<char> ok((size_t)1 << n);
    for (auto& x : ok) x = (char)(rnd(rng, 0, 2) > 0);
    std::vector<unsigned> cur;
    std::vector<std::vector<unsigned>> parts;
    all_partitions(n, 0, cur, parts);
    ll best = INF, cnt = 0;
    int bestk = -1;
    std::vector<ll> bestK(n + 2, INF);
    for (const auto& pt : parts) {
        ll c = 0;
        bool fine = true, okall = true;
        for (unsigned g : pt) {
            if (cost[g] >= INF) fine = false;
            else c += cost[g];
            okall = okall && ok[g];
        }
        if (fine) best = std::min(best, c), bestK[pt.size()] = std::min(bestK[pt.size()], c);
        if (okall) {
            ++cnt;
            if (bestk < 0 || (int)pt.size() < bestk) bestk = (int)pt.size();
        }
    }
    if (n == 0) bestk = 0;
    BD_CHECK(min_partition_cost(n, cost) == best);
    const auto rs = min_partition_restore(n, cost);
    BD_CHECK(rs.first == best);
    if (best < INF) {
        unsigned all = 0;
        ll c = 0;
        for (auto g : rs.second) {
            BD_CHECK((all & g) == 0 && g != 0);
            all |= (unsigned)g, c += cost[g];
        }
        BD_CHECK(all == (1u << n) - 1 && c == best);
    }
    for (int k = 0; k <= n + 1; ++k) BD_CHECK(min_partition_cost_k(n, cost, k) == (n == 0 ? (k == 0 ? 0 : INF) : (k == 0 ? INF : bestK[k])));
    BD_CHECK(min_partition_count(n, ok) == bestk);
    BD_CHECK(count_partitions<ll>(n, ok) == cnt);
    std::vector<ll> a(n);
    for (auto& x : a) x = rnd(rng, -5, 5);
    const auto st = subset_sum_table(a);
    for (size_t m = 0; m < st.size(); ++m) {
        ll s = 0;
        for (int i = 0; i < n; ++i)
            if (m >> i & 1) s += a[i];
        BD_CHECK(st[m] == s);
    }
    cases += 6;
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 150 * TEST_SCALE; ++it) test_tsp(rng, cases);
    for (int it = 0; it < 300 * TEST_SCALE; ++it) test_assign(rng, cases);
    for (int it = 0; it < 150 * TEST_SCALE; ++it) test_partition(rng, cases);
    std::cout << "[dp_bitdp] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[dp_bitdp] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = TEST_SCALE == 1 ? 12 : 17, pn = TEST_SCALE == 1 ? 12 : 18, am = TEST_SCALE == 1 ? 14 : 21;
    std::vector<std::vector<ll>> d(n, std::vector<ll>(n));
    for (auto& r : d)
        for (auto& x : r) x = rnd(rng, 1, 1000000);
    auto t0 = std::chrono::steady_clock::now();
    ll chk = tsp_cycle(d, 0);
    const ll t1 = elapsed_ms(t0);
    std::vector<std::vector<int>> adj(pn, std::vector<int>(pn, 0));
    for (int i = 0; i < pn; ++i)
        for (int j = i + 1; j < pn; ++j) adj[i][j] = adj[j][i] = rnd(rng, 0, 9) < 7;
    t0 = std::chrono::steady_clock::now();
    chk += min_partition_count(pn, clique_table(adj));
    const ll t2 = elapsed_ms(t0);
    std::vector<std::vector<int>> ok(am, std::vector<int>(am));
    for (auto& r : ok)
        for (auto& x : r) x = rnd(rng, 0, 1);
    t0 = std::chrono::steady_clock::now();
    chk += (ll)(count_assignments<unsigned long long>(ok) % 1000);
    const ll t3 = elapsed_ms(t0);
    std::cout << "[dp_bitdp] benchmark: tsp_cycle(n = " << n << ") " << t1 << " ms / min_partition_count(n = " << pn << ") " << t2
              << " ms / count_assignments(m = " << am << ") " << t3 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace bitdp_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace bitdp_demo {

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
#define SHOW(expr, comment) bitdp_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    bitdp_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("bit DP（dp_bitdp） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  bitdp_inf<T> / tsp_table / tsp_cycle / tsp_path / tsp_cycle_restore / tsp_path_restore\n"
              << "  count_assignments<T>\n"
              << "  subset_sum_table / clique_table / independent_table\n"
              << "  min_partition_cost / min_partition_restore / min_partition_cost_k / min_partition_count / count_partitions<T>\n";

    api("【1】 tsp_cycle(dist, start) / tsp_cycle_restore(dist, start)", "T / pair<T, vector<int>>",
        "全頂点を 1 回ずつ通って start に戻る最小コスト / その順番（start, …, start）", "O(2^n n^2)");
    RUN(const long long X = bitdp_inf<long long>());  // 辺が無い
    RUN(const std::vector<std::vector<long long>> d = {{0, 10, 15, 20}, {5, 0, 9, 10}, {6, 13, 0, 12}, {8, 8, 9, 0}});
    SHOW(tsp_cycle(d, 0), "0 → 1 → 3 → 2 → 0 = 10 + 10 + 9 + 6");
    SHOW(tsp_cycle_restore(d, 0), "");
    RUN(const std::vector<std::vector<long long>> d2 = {{0, 1, X}, {X, 0, 1}, {X, X, 0}});
    SHOW(tsp_cycle(d2, 0) == X, "戻れない → bitdp_inf");

    api("【2】 tsp_path(dist, start, goal) / tsp_path_restore(dist, start, goal)", "T / pair<T, vector<int>>",
        "全頂点を 1 回ずつ通る道の最小（start・goal は −1 ならどこでもよい）/ その順番", "O(2^n n^2)");
    SHOW(tsp_path(d), "どこから出てどこで終わってもよい");
    SHOW(tsp_path_restore(d, 2, -1), "2 から出る");
    SHOW(tsp_path(d2, 0, 2), "0 → 1 → 2");

    api("【3】 tsp_table(dist, start)", "vector<vector<T>>", "dp[mask][v]：start から mask を 1 回ずつ通って v にいる最小（一部だけ回る問題に）",
        "O(2^n n^2)");
    RUN(const auto tb = tsp_table(d, 0));
    SHOW(tb[0b1011][3], "{0, 1, 3} を通って 3 にいる：0 → 1 → 3");
    note("例: 「時間 T 以内に回れる都市の数の最大」→ dp[mask][v] + dist[v][goal] ≤ T となる mask の 1 の数の最大");

    api("【4】 count_assignments<T>(ok)", "T", "n 人を m 個の仕事に 1 人 1 つずつ、仕事が重ならないように割り当てる方法の数", "O(2^m m)");
    RUN(const std::vector<std::vector<int>> ok = {{0, 1, 1}, {1, 0, 1}, {1, 1, 1}});
    SHOW(count_assignments<long long>(ok), "3 人・3 つの仕事（EDPC O の形。答えは mod なら T = mint）");
    SHOW(count_assignments<long long>(std::vector<std::vector<int>>{{1, 1, 1}, {1, 1, 1}}), "2 人を 3 つに：3 × 2");

    api("【5】 clique_table(adj) / independent_table(adj) / subset_sum_table(a)", "vector<char> / vector<char> / vector<T>",
        "mask がクリークか / 独立集合か / mask に入る a の和（組のコストを作るのに使う）", "O(2^n)");
    RUN(const std::vector<std::vector<int>> adj = {{0, 1, 1, 0}, {1, 0, 1, 0}, {1, 1, 0, 1}, {0, 0, 1, 0}});
    RUN(const auto cl = clique_table(adj));
    SHOW((int)cl[0b0111], "{0, 1, 2} は三角形");
    SHOW((int)cl[0b1011], "{0, 1, 3}：0 と 3 がつながっていない");
    SHOW(subset_sum_table(std::vector<long long>{3, 5, 9}), "mask = 0〜7 の順");

    api("【6】 min_partition_count(n, ok) / count_partitions<T>(n, ok)", "int / T",
        "ok な組だけで全体を分けるときの最小の組数（無理なら −1）/ 分け方の数", "O(3^n)");
    SHOW(min_partition_count(4, cl), "クリークに分ける最小の数（ABC187 F の形）：{0, 1}, {2, 3} など");
    SHOW(count_partitions<long long>(4, cl), "クリークへの分け方の数");

    api("【7】 min_partition_cost(n, cost) / min_partition_restore(n, cost) / min_partition_cost_k(n, cost, K)", "T / pair<T, vector<mask>> / T",
        "組 g のコスト cost[g] の和を最小に分ける / その組 / ちょうど K 組で", "O(3^n) / O(3^n) / O(K 3^n)");
    {
        code("// 4 人を組に分ける。組のコスト = (人数 − 2)²（2 人組が一番よい）");
        std::vector<long long> cost(16);
        for (int m = 1; m < 16; ++m) cost[m] = (long long)(__builtin_popcount(m) - 2) * (__builtin_popcount(m) - 2);
        SHOW(min_partition_cost(4, cost), "2 人 + 2 人");
        SHOW(min_partition_restore(4, cost), "組の mask の一覧");
        SHOW(min_partition_cost_k(4, cost, 1), "1 組（4 人）なら (4 − 2)²");
        SHOW(min_partition_cost_k(4, cost, 4), "4 組（1 人ずつ）なら 1 × 4");
    }
}

}  // namespace bitdp_demo

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
    bitdp_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    bitdp_test::benchmark();    // 速度計測
    bitdp_demo::demo();         // 使い方デモ
    return 0;
}
