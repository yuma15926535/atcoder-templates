// =============================================================================
//  集める（Gather）：数直線・円周上の点を 1 か所 / K 個連続に集める最小コスト（中央値・重み付き中央値）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  距離の和 Σ|x_i − p| を最小にする p は中央値（重みがあれば重み付き中央値）。K 個を隣り合わせ（t, t+1, …, t+K−1）に
//  集めるときは、小さい順に並べて y_i = x_i − i とすると「y を 1 か所に集める」問題になる
//    例: 0 / 1 の列で、隣どうしの入れ替えで 1 を K 個つなげる最小回数 → gather_consecutive(1 の位置, K).cost
//        家を 1 か所に集める・会議の場所を決める（重み = 人数）・円いテーブルで隣に集まる
//
//  機能一覧（座標は long long の整数。cost は移動距離の和）                                       計算量
//   1 か所に集める
//    gather_point(xs)                     {cost, pos}：Σ|x − pos| が最小の pos（中央値。2 つあれば小さい方）    O(N)
//    gather_point_weighted<T>(xs, ws)     {cost, pos}：Σ w|x − pos| の最小（重み付き中央値。w ≥ 0）              O(N log N)
//    gather_point_square<T>(xs)           {cost, pos}：Σ(x − pos)² の最小（整数の pos。平均に近い整数）          O(N)
//    gather_point_2d(xs, ys)              {cost, x, y}：マンハッタン距離の和の最小（x と y を別々に）             O(N)
//    gather_point_circular(xs, L)         {cost, pos}：周の長さ L の円周上（座標 0〜L−1）で、円周上の距離の和の最小  O(N log N)
//    GatherCost<T> gc(xs) / gc(xs, ws)    Σ w|x − p| をいろいろな p について何度も答える                        O(N log N)
//      gc.cost(p) / gc.best()             p に集めるコスト / {cost, pos}                                         O(log N) / O(1)
//   K 個を選んで集める
//    gather_k_to_point(xs, K)             {cost, pos, first}：K 個を選んで 1 か所に集める最小（小さい順で first 番目から K 個）  O(N log N)
//    gather_consecutive(xs, K)            {cost, first, start}：K 個を選んで start, start+1, …, start+K−1 に 1 個ずつ集める最小
//                                         （xs は相異なる整数。K を省略すると全部）                              O(N log N)
//    gather_consecutive_circular(xs, L, K)  円周（長さ L）で K 個を連続する K マスに集める最小（start は 0〜L−1）  O(N log N)
//
//  注意
//    コストは long long（T を指定できるものは T）。|x| ≤ 10^9・N ≤ 2×10^5 なら距離の和は long long に収まる。
//    重み付き・2 乗は __int128 などを T に指定する（例: gather_point_weighted<__int128>(xs, ws)）
//    xs が空のときは cost = 0、pos = 0
//    gather_consecutive の first は「xs を小さい順に並べたとき」の番号（円周版は first から輪のように K 個）
//  速さの目安: N = 2×10^5 で gather_consecutive 約 0.01 秒、gather_point_circular 約 0.02 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    gather_test::random_test()  … 集める場所・選び方を全部調べる愚直解との突き合わせ
//    gather_test::benchmark()    … 速度計測
//    gather_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace gather_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
inline long long floor_div_ll(long long a, long long b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }
}  // namespace gather_internal

template <class T = long long>
struct GatherPoint {
    T cost;
    long long pos;
};
struct GatherPoint2D {
    long long cost, x, y;
};
struct GatherK {
    long long cost, pos;
    int first;  // 小さい順に並べたときの、選んだ K 個の先頭の番号
};
struct GatherRun {
    long long cost;
    int first;        // 小さい順に並べたときの、選んだ K 個の先頭の番号
    long long start;  // 集めた先の左端（start, start + 1, …, start + K − 1）
};

// ---- 1 か所に集める ----
inline GatherPoint<long long> gather_point(std::vector<long long> xs) {
    if (xs.empty()) return {0, 0};
    const size_t m = (xs.size() - 1) / 2;  // 下側の中央値
    std::nth_element(xs.begin(), xs.begin() + m, xs.end());
    const long long p = xs[m];
    long long c = 0;
    for (long long x : xs) c += x < p ? p - x : x - p;
    return {c, p};
}
template <class T = long long>
GatherPoint<T> gather_point_weighted(const std::vector<long long>& xs, const std::vector<long long>& ws) {
    if (xs.size() != ws.size()) gather_internal::fail("gather_point_weighted", "xs と ws の長さが違う");
    for (long long w : ws)
        if (w < 0) gather_internal::fail("gather_point_weighted", "重みは 0 以上");
    const size_t n = xs.size();
    if (n == 0) return {T(0), 0};
    std::vector<size_t> id(n);
    std::iota(id.begin(), id.end(), 0);
    std::sort(id.begin(), id.end(), [&](size_t a, size_t b) { return xs[a] < xs[b]; });
    T total = T(0);
    for (long long w : ws) total += T(w);
    // 左側の重みの和 × 2 ≥ 全体 となる最初の点（重みが全部 0 なら一番左）
    T acc = T(0);
    long long p = xs[id[0]];
    for (size_t k = 0; k < n; ++k) {
        acc += T(ws[id[k]]);
        if (!(acc * T(2) < total)) {
            p = xs[id[k]];
            break;
        }
    }
    T c = T(0);
    for (size_t i = 0; i < n; ++i) c += T(ws[i]) * T(xs[i] < p ? p - xs[i] : xs[i] - p);
    return {c, p};
}
template <class T = long long>
GatherPoint<T> gather_point_square(const std::vector<long long>& xs) {
    const long long n = (long long)xs.size();
    if (n == 0) return {T(0), 0};
    T s = T(0), s2 = T(0);
    long long sum_ll = 0;  // 平均の候補を出すため（|x| ≤ 4.6×10^13・N ≤ 2×10^5 程度まで）
    for (long long x : xs) s += T(x), s2 += T(x) * T(x), sum_ll += x;
    const long long f = gather_internal::floor_div_ll(sum_ll, n);
    auto cost = [&](long long p) { return s2 - T(2) * T(p) * s + T(n) * T(p) * T(p); };
    const T c0 = cost(f), c1 = cost(f + 1);
    if (c1 < c0) return {c1, f + 1};
    return {c0, f};
}
inline GatherPoint2D gather_point_2d(const std::vector<long long>& xs, const std::vector<long long>& ys) {
    if (xs.size() != ys.size()) gather_internal::fail("gather_point_2d", "xs と ys の長さが違う");
    const auto a = gather_point(xs), b = gather_point(ys);
    return {a.cost + b.cost, a.pos, b.pos};
}

template <class T = long long>
struct GatherCost {
    GatherCost() {}
    explicit GatherCost(const std::vector<long long>& xs) : GatherCost(xs, std::vector<long long>(xs.size(), 1)) {}
    GatherCost(const std::vector<long long>& xs, const std::vector<long long>& ws) {
        if (xs.size() != ws.size()) gather_internal::fail("GatherCost", "xs と ws の長さが違う");
        const size_t n = xs.size();
        std::vector<size_t> id(n);
        std::iota(id.begin(), id.end(), 0);
        std::sort(id.begin(), id.end(), [&](size_t a, size_t b) { return xs[a] < xs[b]; });
        x_.resize(n), pw_.assign(n + 1, T(0)), pwx_.assign(n + 1, T(0));
        for (size_t k = 0; k < n; ++k) {
            if (ws[id[k]] < 0) gather_internal::fail("GatherCost", "重みは 0 以上");
            x_[k] = xs[id[k]];
            pw_[k + 1] = pw_[k] + T(ws[id[k]]);
            pwx_[k + 1] = pwx_[k] + T(ws[id[k]]) * T(xs[id[k]]);
        }
        // 一番よい場所（重み付き中央値）
        best_pos_ = n ? x_[0] : 0;
        for (size_t k = 0; k < n; ++k)
            if (!(pw_[k + 1] * T(2) < pw_[n])) {
                best_pos_ = x_[k];
                break;
            }
    }
    T cost(long long p) const {
        const size_t n = x_.size();
        const size_t k = (size_t)(std::lower_bound(x_.begin(), x_.end(), p) - x_.begin());  // x < p の数
        return T(p) * pw_[k] - pwx_[k] + (pwx_[n] - pwx_[k]) - T(p) * (pw_[n] - pw_[k]);
    }
    GatherPoint<T> best() const { return {cost(best_pos_), best_pos_}; }

   private:
    std::vector<long long> x_;
    std::vector<T> pw_, pwx_;  // 重みの累積和・重み × 座標の累積和
    long long best_pos_ = 0;
};

// 周の長さ L の円周上（座標 0〜L−1）。円周上の距離 min(|a − b|, L − |a − b|) の和が最小の点（候補は xs のどれか）
inline GatherPoint<long long> gather_point_circular(const std::vector<long long>& xs, long long L) {
    if (L <= 0) gather_internal::fail("gather_point_circular", "L は 1 以上");
    for (long long x : xs)
        if (x < 0 || x >= L) gather_internal::fail("gather_point_circular", "座標は 0〜L−1");
    const size_t n = xs.size();
    if (n == 0) return {0, 0};
    std::vector<long long> s = xs;
    std::sort(s.begin(), s.end());
    std::vector<long long> v(3 * n);  // x − L, x, x + L
    for (size_t i = 0; i < n; ++i) v[i] = s[i] - L, v[n + i] = s[i], v[2 * n + i] = s[i] + L;
    std::vector<long long> pre(3 * n + 1, 0);
    for (size_t i = 0; i < 3 * n; ++i) pre[i + 1] = pre[i] + v[i];
    GatherPoint<long long> best{std::numeric_limits<long long>::max(), 0};
    size_t lo = 0, mid = 0;  // p が増えると lo・mid も増える（尺取り）
    for (size_t i = 0; i < n; ++i) {
        const long long p = s[i];
        // p − L/2 < v ≤ p + L/2 にある n 個（2v − 2p ∈ (−L, L]）。lo = その先頭、mid = v ≥ p の先頭
        while (2 * (v[lo] - p) <= -L) ++lo;
        while (v[mid] < p) ++mid;
        const size_t hi = lo + n;
        const long long c = p * (long long)(mid - lo) - (pre[mid] - pre[lo]) + (pre[hi] - pre[mid]) - p * (long long)(hi - mid);
        if (c < best.cost) best = {c, p};
    }
    return best;
}

// ---- K 個を選んで集める ----
inline GatherK gather_k_to_point(const std::vector<long long>& xs, int K) {
    const int n = (int)xs.size();
    if (K < 1 || K > n) gather_internal::fail("gather_k_to_point", "1 ≤ K ≤ N にする");
    std::vector<long long> s = xs;
    std::sort(s.begin(), s.end());
    std::vector<long long> pre(n + 1, 0);
    for (int i = 0; i < n; ++i) pre[i + 1] = pre[i] + s[i];
    GatherK best{std::numeric_limits<long long>::max(), 0, 0};
    for (int i = 0; i + K <= n; ++i) {
        const int m = i + (K - 1) / 2;
        const long long p = s[m];
        const long long c = p * (m - i) - (pre[m] - pre[i]) + (pre[i + K] - pre[m]) - p * (i + K - m);
        if (c < best.cost) best = {c, p, i};
    }
    return best;
}
namespace gather_internal {
// y（小さい順）の窓 [i, i + K) を中央値に集めるコストの最小。返り値 {cost, i, 中央値}
inline std::pair<long long, std::pair<int, long long>> best_window(const std::vector<long long>& y, int starts, int K) {
    std::vector<long long> pre(y.size() + 1, 0);
    for (size_t i = 0; i < y.size(); ++i) pre[i + 1] = pre[i] + y[i];
    long long bc = std::numeric_limits<long long>::max();
    int bi = 0;
    long long bm = 0;
    for (int i = 0; i < starts; ++i) {
        const int m = i + (K - 1) / 2;
        const long long p = y[m];
        const long long c = p * (m - i) - (pre[m] - pre[i]) + (pre[i + K] - pre[m]) - p * (i + K - m);
        if (c < bc) bc = c, bi = i, bm = p;
    }
    return {bc, {bi, bm}};
}
}  // namespace gather_internal
inline GatherRun gather_consecutive(const std::vector<long long>& xs, int K) {
    const int n = (int)xs.size();
    if (K < 0 || K > n) gather_internal::fail("gather_consecutive", "0 ≤ K ≤ N にする");
    if (K == 0) return {0, 0, 0};
    std::vector<long long> s = xs;
    std::sort(s.begin(), s.end());
    for (int i = 0; i + 1 < n; ++i)
        if (s[i] == s[i + 1]) gather_internal::fail("gather_consecutive", "xs は相異なる整数にする");
    std::vector<long long> y(n);
    for (int i = 0; i < n; ++i) y[i] = s[i] - i;  // 小さい順（相異なるので減らない）
    const auto r = gather_internal::best_window(y, n - K + 1, K);
    const int i = r.second.first;
    return {r.first, i, r.second.second + i};
}
inline GatherRun gather_consecutive(const std::vector<long long>& xs) { return gather_consecutive(xs, (int)xs.size()); }
inline GatherRun gather_consecutive_circular(const std::vector<long long>& xs, long long L, int K) {
    const int n = (int)xs.size();
    if (L <= 0) gather_internal::fail("gather_consecutive_circular", "L は 1 以上");
    if (K < 0 || K > n || K > L) gather_internal::fail("gather_consecutive_circular", "0 ≤ K ≤ min(N, L) にする");
    for (long long x : xs)
        if (x < 0 || x >= L) gather_internal::fail("gather_consecutive_circular", "座標は 0〜L−1");
    if (K == 0) return {0, 0, 0};
    std::vector<long long> s = xs;
    std::sort(s.begin(), s.end());
    for (int i = 0; i + 1 < n; ++i)
        if (s[i] == s[i + 1]) gather_internal::fail("gather_consecutive_circular", "xs は相異なる整数にする");
    // 2 周分に伸ばし（x, x + L）、窓の先頭を 0〜N−1 で動かす
    std::vector<long long> y(2 * n);
    for (int i = 0; i < 2 * n; ++i) y[i] = (i < n ? s[i] : s[i - n] + L) - i;
    const auto r = gather_internal::best_window(y, n, K);
    const int i = r.second.first;
    long long st = (r.second.second + i) % L;
    if (st < 0) st += L;
    return {r.first, i, st};
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace gather_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define GA_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[gather_cost] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static ll cdist(ll a, ll b, ll L) {
    const ll d = a > b ? a - b : b - a;
    return std::min(d, L - d);
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 8);
        const ll V = rnd(rng, 1, 12);
        std::vector<ll> xs(n), ws(n), ys(n);
        for (auto& x : xs) x = rnd(rng, -V, V);
        for (auto& w : ws) w = rnd(rng, 0, 3) ? rnd(rng, 0, 5) : 0;
        for (auto& y : ys) y = rnd(rng, -V, V);
        // 1 か所
        auto lin = [&](ll p) {
            ll c = 0;
            for (ll x : xs) c += std::abs(x - p);
            return c;
        };
        auto wlin = [&](ll p) {
            ll c = 0;
            for (int i = 0; i < n; ++i) c += ws[i] * std::abs(xs[i] - p);
            return c;
        };
        auto sq = [&](ll p) {
            ll c = 0;
            for (ll x : xs) c += (x - p) * (x - p);
            return c;
        };
        ll b1 = 0, b2 = 0, b3 = 0;
        if (n > 0) {
            b1 = b2 = b3 = std::numeric_limits<ll>::max();
            for (ll p = -V - 2; p <= V + 2; ++p) b1 = std::min(b1, lin(p)), b2 = std::min(b2, wlin(p)), b3 = std::min(b3, sq(p));
        }
        const auto g1 = gather_point(xs);
        GA_CHECK(g1.cost == b1 && lin(g1.pos) == b1);
        const auto g2 = gather_point_weighted(xs, ws);
        GA_CHECK(g2.cost == b2 && (n == 0 || wlin(g2.pos) == b2));
        const auto g3 = gather_point_square(xs);
        GA_CHECK(g3.cost == b3 && (n == 0 || sq(g3.pos) == b3));
        const auto g4 = gather_point_weighted<__int128>(xs, ws);
        GA_CHECK((ll)g4.cost == b2);
        const auto g2d = gather_point_2d(xs, ys);
        ll by = 0;
        if (n > 0) {
            by = std::numeric_limits<ll>::max();
            for (ll p = -V - 2; p <= V + 2; ++p) {
                ll c = 0;
                for (ll y : ys) c += std::abs(y - p);
                by = std::min(by, c);
            }
        }
        GA_CHECK(g2d.cost == b1 + by);
        const GatherCost<ll> gc(xs, ws), gu(xs);
        for (ll p = -V - 3; p <= V + 3; ++p) GA_CHECK(gc.cost(p) == wlin(p) && gu.cost(p) == lin(p));
        GA_CHECK(gc.best().cost == b2 && gu.best().cost == b1);
        cases += 8;
        // K 個を 1 か所に
        for (int K = 1; K <= n; ++K) {
            ll best = std::numeric_limits<ll>::max();
            for (int m = 0; m < (1 << n); ++m) {
                if (__builtin_popcount(m) != K) continue;
                for (ll p = -V; p <= V; ++p) {
                    ll c = 0;
                    for (int i = 0; i < n; ++i)
                        if (m >> i & 1) c += std::abs(xs[i] - p);
                    best = std::min(best, c);
                }
            }
            const auto r = gather_k_to_point(xs, K);
            GA_CHECK(r.cost == best);
            std::vector<ll> s = xs;
            std::sort(s.begin(), s.end());
            ll c = 0;
            for (int i = r.first; i < r.first + K; ++i) c += std::abs(s[i] - r.pos);
            GA_CHECK(c == best);
            ++cases;
        }
        // 円周の 1 か所
        {
            const ll L = rnd(rng, 1, 12);
            std::vector<ll> cs(n);
            for (auto& x : cs) x = rnd(rng, 0, L - 1);
            ll best = 0;
            if (n > 0) {
                best = std::numeric_limits<ll>::max();
                for (ll p = 0; p < L; ++p) {
                    ll c = 0;
                    for (ll x : cs) c += cdist(x, p, L);
                    best = std::min(best, c);
                }
            }
            const auto r = gather_point_circular(cs, L);
            GA_CHECK(r.cost == best);
            if (n > 0) {
                ll c = 0;
                for (ll x : cs) c += cdist(x, r.pos, L);
                GA_CHECK(c == best && 0 <= r.pos && r.pos < L);
            }
            ++cases;
        }
    }
    // K 個を連続に（相異なる点。直線・円周）
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const ll L = rnd(rng, 1, 9);
        const int n = (int)rnd(rng, 0, std::min<ll>(L, 6));
        std::vector<ll> all(L);
        std::iota(all.begin(), all.end(), 0);
        std::shuffle(all.begin(), all.end(), rng);
        std::vector<ll> xs(all.begin(), all.begin() + n);
        const ll shift = rnd(rng, -5, 5);
        std::vector<ll> lx = xs;
        for (auto& x : lx) x += shift;
        for (int K = 0; K <= n; ++K) {
            // 直線：選ぶ K 個・行き先の左端 t・対応のさせ方を全部
            ll best = K == 0 ? 0 : std::numeric_limits<ll>::max(), bestc = best;
            if (K > 0) {
                std::vector<int> idx(n);
                std::iota(idx.begin(), idx.end(), 0);
                for (int m = 0; m < (1 << n); ++m) {
                    if (__builtin_popcount(m) != K) continue;
                    std::vector<int> sel;
                    for (int i = 0; i < n; ++i)
                        if (m >> i & 1) sel.push_back(i);
                    std::sort(sel.begin(), sel.end());
                    do {
                        for (ll t = -16; t <= 16; ++t) {
                            ll c = 0;
                            for (int j = 0; j < K; ++j) c += std::abs(lx[sel[j]] - (t + j));
                            best = std::min(best, c);
                        }
                        for (ll t = 0; t < L; ++t) {
                            ll c = 0;
                            for (int j = 0; j < K; ++j) c += cdist(xs[sel[j]], (t + j) % L, L);
                            bestc = std::min(bestc, c);
                        }
                    } while (std::next_permutation(sel.begin(), sel.end()));
                }
            }
            const auto r = gather_consecutive(lx, K);
            GA_CHECK(r.cost == best);
            if (K > 0) {
                std::vector<ll> s = lx;
                std::sort(s.begin(), s.end());
                ll c = 0;
                for (int j = 0; j < K; ++j) c += std::abs(s[r.first + j] - (r.start + j));
                GA_CHECK(c == best);
            }
            const auto rc = gather_consecutive_circular(xs, L, K);
            GA_CHECK(rc.cost == bestc);
            if (K > 0) {
                std::vector<ll> s = xs;
                std::sort(s.begin(), s.end());
                ll c = 0;
                for (int j = 0; j < K; ++j) c += cdist(s[(rc.first + j) % n], (rc.start + j) % L, L);
                GA_CHECK(c == bestc && 0 <= rc.start && rc.start < L);
            }
            cases += 2;
        }
        if (n > 0) GA_CHECK(gather_consecutive(lx).cost == gather_consecutive(lx, n).cost);
    }
    std::cout << "[gather_cost] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[gather_cost] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    const ll L = 1000000000;
    std::vector<ll> xs(N);
    {
        std::vector<ll> t;
        while ((int)t.size() < N) {
            for (int i = 0; i < N; ++i) t.push_back(rnd(rng, 0, L - 1));
            std::sort(t.begin(), t.end());
            t.erase(std::unique(t.begin(), t.end()), t.end());
        }
        t.resize(N);
        std::shuffle(t.begin(), t.end(), rng);
        xs = t;
    }
    ll chk = 0;
    auto t0 = std::chrono::steady_clock::now();
    chk += gather_consecutive(xs, N / 2).cost;
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    chk += gather_point_circular(xs, L).cost;
    const ll t2 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    chk += gather_consecutive_circular(xs, L, N / 2).cost;
    const ll t3 = elapsed_ms(t0);
    std::cout << "[gather_cost] benchmark (N = " << N << "): gather_consecutive " << t1 << " ms / gather_point_circular " << t2
              << " ms / gather_consecutive_circular " << t3 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace gather_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace gather_demo {

static std::string to_s(const GatherPoint<long long>& g);
static std::string to_s(const GatherPoint2D& g);
static std::string to_s(const GatherK& g);
static std::string to_s(const GatherRun& g);

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
#define SHOW(expr, comment) gather_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    gather_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static std::string to_s(const GatherPoint<long long>& g) { return "{cost = " + to_s(g.cost) + ", pos = " + to_s(g.pos) + "}"; }
static std::string to_s(const GatherPoint2D& g) {
    return "{cost = " + to_s(g.cost) + ", x = " + to_s(g.x) + ", y = " + to_s(g.y) + "}";
}
static std::string to_s(const GatherK& g) {
    return "{cost = " + to_s(g.cost) + ", pos = " + to_s(g.pos) + ", first = " + to_s(g.first) + "}";
}
static std::string to_s(const GatherRun& g) {
    return "{cost = " + to_s(g.cost) + ", first = " + to_s(g.first) + ", start = " + to_s(g.start) + "}";
}

static void demo() {
    banner("集める（gather_cost） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  gather_point / gather_point_weighted<T> / gather_point_square<T> / gather_point_2d / gather_point_circular\n"
              << "  GatherCost<T>: cost(p) / best()\n"
              << "  gather_k_to_point / gather_consecutive / gather_consecutive_circular\n";

    api("【1】 gather_point(xs)", "GatherPoint {cost, pos}", "Σ|x − pos| が最小の場所（中央値）とそのコスト", "O(N)");
    RUN(const std::vector<long long> xs = {1, 4, 9, 10, 30});
    SHOW(gather_point(xs), "中央値 9：8 + 5 + 0 + 1 + 21");

    api("【2】 gather_point_weighted<T>(xs, ws)", "GatherPoint<T>", "Σ w|x − pos| の最小（重み付き中央値）。大きくなるなら T = __int128", "O(N log N)");
    RUN(const std::vector<long long> ws = {1, 1, 1, 1, 10});
    SHOW(gather_point_weighted(xs, ws), "30 に 10 人いるので 30 に集まる");

    api("【3】 gather_point_square<T>(xs) / gather_point_2d(xs, ys)", "GatherPoint<T> / GatherPoint2D {cost, x, y}",
        "Σ(x − pos)² の最小（平均に近い整数）/ マンハッタン距離の和の最小（x と y を別々に中央値）", "O(N)");
    SHOW(gather_point_square(xs), "平均 10.8 → 11");
    SHOW(gather_point_2d(std::vector<long long>{0, 2, 10}, std::vector<long long>{5, 1, 3}), "(2, 3) に集める");

    api("【4】 gather_point_circular(xs, L)", "GatherPoint", "周の長さ L の円周上（0〜L−1）で、円周上の距離の和の最小", "O(N log N)");
    SHOW(gather_point_circular(std::vector<long long>{1, 2, 9}, 10LL), "9 → 1 は 2 歩（0 をまたぐ）");

    api("【5】 GatherCost<T> gc(xs) / gc(xs, ws) → gc.cost(p) / gc.best()", "— / T / GatherPoint<T>", "いろいろな p で Σ w|x − p| を何度も答える",
        "O(N log N) / O(log N)");
    RUN(const GatherCost<> gc(xs));
    SHOW(gc.cost(0), "");
    SHOW(gc.cost(9), "");
    SHOW(gc.best(), "");

    api("【6】 gather_k_to_point(xs, K)", "GatherK {cost, pos, first}", "K 個を選んで 1 か所に集める最小（小さい順で first 番目から K 個を pos へ）",
        "O(N log N)");
    SHOW(gather_k_to_point(xs, 3), "{4, 9, 10} を 9 へ");

    api("【7】 gather_consecutive(xs, K)", "GatherRun {cost, first, start}",
        "K 個を選んで start, start+1, …, start+K−1 に 1 個ずつ集める最小（xs は相異なる）。K を省略すると全部", "O(N log N)");
    SHOW(gather_consecutive(xs, 3), "{4, 9, 10} → {8, 9, 10}：4 + 0 + 0");
    SHOW(gather_consecutive(xs), "全部を 5 マス連続に");
    {
        code("// 0 / 1 の列で、隣どうしの入れ替えで 1 を全部つなげる最小回数 = 1 の位置を連続に集めるコスト");
        const std::string s = "1001011";
        std::vector<long long> pos;
        for (int i = 0; i < (int)s.size(); ++i)
            if (s[i] == '1') pos.push_back(i);
        note("s = \"1001011\" → 1 の位置 " + to_s(pos) + " → gather_consecutive(pos).cost = " + to_s(gather_consecutive(pos).cost));
    }

    api("【8】 gather_consecutive_circular(xs, L, K)", "GatherRun", "円周（長さ L）で K 個を連続する K マスに集める最小（start は 0〜L−1）",
        "O(N log N)");
    SHOW(gather_consecutive_circular(std::vector<long long>{0, 1, 8, 9}, 10LL, 4), "8, 9, 0, 1 はもう連続（0 をまたぐ）");
    SHOW(gather_consecutive_circular(std::vector<long long>{0, 3, 6}, 10LL, 3), "2, 3, 4 に集める：0 → 2、6 → 4");
}

}  // namespace gather_demo

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
    gather_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    gather_test::benchmark();    // 速度計測
    gather_demo::demo();         // 使い方デモ
    return 0;
}
