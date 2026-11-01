// =============================================================================
//  反転（Flip）：端から決める貪欲（長さ K の区間の反転・+1 mod m・長方形の反転）とライツアウト   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  「長さ K の区間をまとめて反転する」操作は、順番を入れ替えても結果が同じで、同じ区間を 2 回やるのは無駄。
//  だから左端のマスを直せるのは「そのマスから始まる区間」だけになり、左から順に決まる（貪欲）。
//  反転した回数はいもす（差分）で持つので全体 O(N)。2 次元でも左上から同じように決まる
//    例: 「プログラミングコンテストチャレンジブック」の Face The Right Way（POJ 3276）… flip_k_best
//        Fliptile（POJ 3279）… lights_out（1 行目の押し方を全部試すと、残りの行は決まる）
//
//  機能一覧（a は 0 / 1 の vector<int>。target = 全部をそろえたい値（0 か 1））                    計算量
//   1 次元
//    flip_k_min(a, K, target = 0)          長さ K の区間の反転で全部 target にする最小回数（無理なら −1）        O(N)
//    flip_k_starts(a, K, target = 0)       そのとき反転する区間の左端の一覧（小さい順）。無理なら nullopt      O(N)
//    flip_k_best(a, target = 0)            {K, 回数}：回数が最小になる K（同じなら小さい K）。どの K でも無理なら {−1, −1}  O(N^2)
//    add_k_mod_min(a, K, m)                長さ K の区間に +1（mod m）する操作で全部 0（mod m）にする最小回数（無理なら −1）
//                                          a は 0 以上の整数でなくてもよい（mod m で考える）                    O(N)
//   2 次元（g は H × W の 0 / 1）
//    flip_rect_min(g, h, w, target = 0)    h × w の長方形の反転で全部 target にする最小回数（無理なら −1）      O(HW)
//    flip_rect_plan(g, h, w, target = 0)   press[i][j] = 1 なら左上 (i, j) の長方形を反転。無理なら nullopt     O(HW)
//    lights_out(g, target = 0)             マスを押すと自分と上下左右が反転。{最小の押す回数, 押し方 press}。
//                                          最小の中で press を行ごとに並べた列が辞書順最小のもの。無理なら {−1, {}}  O(2^W H)（W ≤ 25 くらい）
//
//  注意
//    K < 1 や K > N は止める。h, w も同じ（1 ≤ h ≤ H、1 ≤ w ≤ W）
//    lights_out は 1 行目を全部（2^W 通り）試す。W が大きく H が小さいときは転置してから渡す（W ≤ 63）
//    add_k_mod_min の回数は long long（最大 N × (m − 1)）
//  速さの目安: flip_k_min（N = 10^6）約 0.006 秒、flip_k_best（N = 5000）約 0.04 秒、lights_out（20 × 20）約 0.1 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    flip_test::random_test()  … 押し方を全部試す・幅優先探索の愚直解との突き合わせ
//    flip_test::benchmark()    … 速度計測
//    flip_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace flip_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
inline void check01(const std::vector<int>& a, int target, const char* f) {
    for (int x : a)
        if (x != 0 && x != 1) fail(f, "a は 0 か 1 だけ");
    if (target != 0 && target != 1) fail(f, "target は 0 か 1");
}
// 長さ K の反転で全部 target。starts に左端を入れる（nullptr なら入れない）。無理なら −1
inline int flip_core(const std::vector<int>& a, int K, int target, std::vector<int>* starts, const char* f) {
    const int n = (int)a.size();
    if (K < 1 || K > n) fail(f, "1 ≤ K ≤ N にする");
    std::vector<char> ends(n + 1, 0);  // ends[i] = i で反転の効果が終わる区間の数（の偶奇）
    int par = 0, cnt = 0;
    for (int i = 0; i < n; ++i) {
        par ^= ends[i];
        if ((a[i] ^ par) != target) {
            if (i + K > n) return -1;
            par ^= 1, ends[i + K] ^= 1, ++cnt;
            if (starts) starts->push_back(i);
        }
    }
    return cnt;
}
}  // namespace flip_internal

// ---- 1 次元 ----
inline int flip_k_min(const std::vector<int>& a, int K, int target = 0) {
    flip_internal::check01(a, target, "flip_k_min");
    return flip_internal::flip_core(a, K, target, nullptr, "flip_k_min");
}
inline std::optional<std::vector<int>> flip_k_starts(const std::vector<int>& a, int K, int target = 0) {
    flip_internal::check01(a, target, "flip_k_starts");
    std::vector<int> st;
    if (flip_internal::flip_core(a, K, target, &st, "flip_k_starts") < 0) return std::nullopt;
    return st;
}
inline std::pair<int, int> flip_k_best(const std::vector<int>& a, int target = 0) {
    flip_internal::check01(a, target, "flip_k_best");
    std::pair<int, int> best{-1, -1};
    for (int K = 1; K <= (int)a.size(); ++K) {
        const int c = flip_internal::flip_core(a, K, target, nullptr, "flip_k_best");
        if (c >= 0 && (best.second < 0 || c < best.second)) best = {K, c};
    }
    return best;
}
inline long long add_k_mod_min(const std::vector<long long>& a, int K, long long m) {
    const int n = (int)a.size();
    if (K < 1 || K > n) flip_internal::fail("add_k_mod_min", "1 ≤ K ≤ N にする");
    if (m < 1) flip_internal::fail("add_k_mod_min", "m は 1 以上");
    std::vector<long long> ends(n + 1, 0);  // ends[i] = i で効果が終わる足し算の合計（mod m）
    long long add = 0, cnt = 0;
    for (int i = 0; i < n; ++i) {
        add = ((add - ends[i]) % m + m) % m;
        const long long cur = ((a[i] % m + m) % m + add) % m;
        if (cur != 0) {
            if (i + K > n) return -1;
            const long long t = m - cur;  // この位置から始まる区間に t 回
            cnt += t, add = (add + t) % m, ends[i + K] = (ends[i + K] + t) % m;
        }
    }
    return cnt;
}

// ---- 2 次元 ----
inline std::optional<std::vector<std::vector<int>>> flip_rect_plan(const std::vector<std::vector<int>>& g, int h, int w, int target = 0) {
    const int H = (int)g.size(), W = H ? (int)g[0].size() : 0;
    for (const auto& r : g) {
        if ((int)r.size() != W) flip_internal::fail("flip_rect_plan", "行の長さがそろっていない");
        flip_internal::check01(r, target, "flip_rect_plan");
    }
    if (H == 0 || W == 0) return std::vector<std::vector<int>>(H, std::vector<int>(W, 0));
    if (h < 1 || h > H || w < 1 || w > W) flip_internal::fail("flip_rect_plan", "1 ≤ h ≤ H、1 ≤ w ≤ W にする");
    std::vector<std::vector<char>> d(H + 1, std::vector<char>(W + 1, 0)), P(H + 1, std::vector<char>(W + 1, 0));
    std::vector<std::vector<int>> press(H, std::vector<int>(W, 0));
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) {
            // P[i+1][j+1] = 左上 (0, 0)〜(i, j) の d の XOR = (i, j) を覆う反転の回数の偶奇
            char p = d[i][j] ^ P[i][j + 1] ^ P[i + 1][j] ^ P[i][j];
            if ((g[i][j] ^ p) != target) {
                if (i + h > H || j + w > W) return std::nullopt;
                press[i][j] = 1;
                d[i][j] ^= 1, d[i][j + w] ^= 1, d[i + h][j] ^= 1, d[i + h][j + w] ^= 1;
                p ^= 1;
            }
            P[i + 1][j + 1] = p;
        }
    return press;
}
inline long long flip_rect_min(const std::vector<std::vector<int>>& g, int h, int w, int target = 0) {
    const auto pl = flip_rect_plan(g, h, w, target);
    if (!pl) return -1;
    long long c = 0;
    for (const auto& r : *pl)
        for (int x : r) c += x;
    return c;
}
inline std::pair<int, std::vector<std::vector<int>>> lights_out(const std::vector<std::vector<int>>& g, int target = 0) {
    const int H = (int)g.size(), W = H ? (int)g[0].size() : 0;
    for (const auto& r : g) {
        if ((int)r.size() != W) flip_internal::fail("lights_out", "行の長さがそろっていない");
        flip_internal::check01(r, target, "lights_out");
    }
    if (H == 0 || W == 0) return {0, std::vector<std::vector<int>>(H, std::vector<int>(W, 0))};
    if (W > 30) flip_internal::fail("lights_out", "W が大きすぎる（2^W 通り試す。転置して W を小さくする）");
    const uint64_t full = (W == 64) ? ~0ULL : ((1ULL << W) - 1);
    // 行をビットで（列 j を W − 1 − j ビット目に：列 0 が一番上の桁 → 1 行目を小さい順に試すと辞書順）
    std::vector<uint64_t> row(H, 0);
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j)
            if (g[i][j] != target) row[i] |= 1ULL << (W - 1 - j);
    int best = -1;
    uint64_t best_first = 0;
    std::vector<uint64_t> cur(H), pr(H);
    for (uint64_t first = 0; first <= full; ++first) {
        cur = row;
        int cnt = 0;
        for (int i = 0; i < H; ++i) {
            const uint64_t p = i == 0 ? first : cur[i - 1];  // 上の行で点いているところの下を押す
            pr[i] = p;
            cnt += __builtin_popcountll(p);
            if (best >= 0 && cnt >= best) break;
            cur[i] ^= (p ^ (p << 1) ^ (p >> 1)) & full;
            if (i > 0) cur[i - 1] ^= p;
            if (i + 1 < H) cur[i + 1] ^= p;
        }
        if ((best < 0 || cnt < best) && cur[H - 1] == 0) {
            bool ok = true;
            for (int i = 0; i < H && ok; ++i) ok = cur[i] == 0;
            if (ok) best = cnt, best_first = first;
        }
        if (first == full) break;
    }
    if (best < 0) return {-1, {}};
    // 押し方を作り直す
    std::vector<std::vector<int>> press(H, std::vector<int>(W, 0));
    cur = row;
    for (int i = 0; i < H; ++i) {
        const uint64_t p = i == 0 ? best_first : cur[i - 1];
        for (int j = 0; j < W; ++j) press[i][j] = (int)(p >> (W - 1 - j) & 1);
        cur[i] ^= (p ^ (p << 1) ^ (p >> 1)) & full;
        if (i > 0) cur[i - 1] ^= p;
        if (i + 1 < H) cur[i + 1] ^= p;
    }
    return {best, press};
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace flip_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define FL_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[flip_greedy] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static void test_1d(std::mt19937_64& rng, long long& cases) {
    const int n = (int)rnd(rng, 1, 10);
    std::vector<int> a(n);
    for (auto& x : a) x = (int)rnd(rng, 0, 1);
    for (int target = 0; target <= 1; ++target) {
        std::pair<int, int> best{-1, -1};
        for (int K = 1; K <= n; ++K) {
            // 幅優先探索（状態 = ビット列）
            const int full = (1 << n) - 1;
            int s0 = 0, goal = target ? full : 0;
            for (int i = 0; i < n; ++i) s0 |= a[i] << i;
            std::vector<int> dist(1 << n, -1);
            std::queue<int> q;
            dist[s0] = 0, q.push(s0);
            while (!q.empty()) {
                const int s = q.front();
                q.pop();
                for (int l = 0; l + K <= n; ++l) {
                    const int t = s ^ (((1 << K) - 1) << l);
                    if (dist[t] < 0) dist[t] = dist[s] + 1, q.push(t);
                }
            }
            const int r = flip_k_min(a, K, target);
            FL_CHECK(r == dist[goal]);
            const auto st = flip_k_starts(a, K, target);
            FL_CHECK(st.has_value() == (r >= 0));
            if (st) {
                FL_CHECK((int)st->size() == r);
                std::vector<int> b = a;
                for (int l : *st)
                    for (int i = l; i < l + K; ++i) b[i] ^= 1;
                for (int x : b) FL_CHECK(x == target);
                FL_CHECK(std::is_sorted(st->begin(), st->end()));
            }
            if (r >= 0 && (best.second < 0 || r < best.second)) best = {K, r};
            ++cases;
        }
        FL_CHECK(flip_k_best(a, target) == best);
    }
    // +1 mod m
    const ll m = rnd(rng, 1, 4);
    const int n2 = (int)rnd(rng, 1, 6);
    std::vector<ll> b(n2);
    for (auto& x : b) x = rnd(rng, -10, 10);
    for (int K = 1; K <= n2; ++K) {
        // 状態 = 各位置の値 mod m（m^n2 ≤ 4^6）
        int S = 1;
        for (int i = 0; i < n2; ++i) S *= (int)m;
        auto enc = [&](const std::vector<int>& v) {
            int c = 0;
            for (int i = n2 - 1; i >= 0; --i) c = c * (int)m + v[i];
            return c;
        };
        std::vector<int> v0(n2);
        for (int i = 0; i < n2; ++i) v0[i] = (int)(((b[i] % m) + m) % m);
        std::vector<int> dist(S, -1);
        std::queue<int> q;
        dist[enc(v0)] = 0, q.push(enc(v0));
        while (!q.empty()) {
            const int s = q.front();
            q.pop();
            std::vector<int> v(n2);
            int t = s;
            for (int i = 0; i < n2; ++i) v[i] = t % (int)m, t /= (int)m;
            for (int l = 0; l + K <= n2; ++l) {
                std::vector<int> u = v;
                for (int i = l; i < l + K; ++i) u[i] = (u[i] + 1) % (int)m;
                const int e = enc(u);
                if (dist[e] < 0) dist[e] = dist[s] + 1, q.push(e);
            }
        }
        FL_CHECK(add_k_mod_min(b, K, m) == dist[0]);
        ++cases;
    }
}

static void test_2d(std::mt19937_64& rng, long long& cases) {
    const int H = (int)rnd(rng, 1, 4), W = (int)rnd(rng, 1, 4);
    std::vector<std::vector<int>> g(H, std::vector<int>(W));
    for (auto& r : g)
        for (auto& x : r) x = (int)rnd(rng, 0, 1);
    for (int target = 0; target <= 1; ++target) {
        // 長方形：使う左上の集合を全部
        for (int h = 1; h <= H; ++h)
            for (int w = 1; w <= W; ++w) {
                const int P = (H - h + 1) * (W - w + 1);
                int best = -1;
                for (int m = 0; m < (1 << P); ++m) {
                    std::vector<std::vector<int>> b = g;
                    for (int k = 0; k < P; ++k)
                        if (m >> k & 1) {
                            const int i0 = k / (W - w + 1), j0 = k % (W - w + 1);
                            for (int i = i0; i < i0 + h; ++i)
                                for (int j = j0; j < j0 + w; ++j) b[i][j] ^= 1;
                        }
                    bool ok = true;
                    for (auto& r : b)
                        for (int x : r) ok = ok && x == target;
                    if (ok && (best < 0 || __builtin_popcount(m) < best)) best = __builtin_popcount(m);
                }
                FL_CHECK(flip_rect_min(g, h, w, target) == best);
                const auto pl = flip_rect_plan(g, h, w, target);
                FL_CHECK(pl.has_value() == (best >= 0));
                if (pl) {
                    std::vector<std::vector<int>> b = g;
                    for (int i0 = 0; i0 < H; ++i0)
                        for (int j0 = 0; j0 < W; ++j0)
                            if ((*pl)[i0][j0]) {
                                FL_CHECK(i0 + h <= H && j0 + w <= W);
                                for (int i = i0; i < i0 + h; ++i)
                                    for (int j = j0; j < j0 + w; ++j) b[i][j] ^= 1;
                            }
                    for (auto& r : b)
                        for (int x : r) FL_CHECK(x == target);
                }
                ++cases;
            }
        // ライツアウト：押すマスの集合を全部（辞書順最小も確かめる）
        const int C = H * W;
        int best = -1;
        std::vector<int> bestp;
        for (int m = 0; m < (1 << C); ++m) {
            std::vector<std::vector<int>> b = g;
            for (int k = 0; k < C; ++k)
                if (m >> k & 1) {
                    const int i = k / W, j = k % W;
                    b[i][j] ^= 1;
                    if (i > 0) b[i - 1][j] ^= 1;
                    if (i + 1 < H) b[i + 1][j] ^= 1;
                    if (j > 0) b[i][j - 1] ^= 1;
                    if (j + 1 < W) b[i][j + 1] ^= 1;
                }
            bool ok = true;
            for (auto& r : b)
                for (int x : r) ok = ok && x == target;
            if (!ok) continue;
            std::vector<int> pv(C);
            for (int k = 0; k < C; ++k) pv[k] = m >> k & 1;
            const int c = __builtin_popcount(m);
            if (best < 0 || c < best || (c == best && pv < bestp)) best = c, bestp = pv;
        }
        const auto lo = lights_out(g, target);
        FL_CHECK(lo.first == best);
        if (best >= 0) {
            std::vector<int> pv;
            for (auto& r : lo.second)
                for (int x : r) pv.push_back(x);
            FL_CHECK(pv == bestp);
        }
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 150 * TEST_SCALE; ++it) test_1d(rng, cases);
    for (int it = 0; it < 150 * TEST_SCALE; ++it) test_2d(rng, cases);
    std::cout << "[flip_greedy] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[flip_greedy] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 100000 : 1000000, M = TEST_SCALE == 1 ? 800 : 5000, G = TEST_SCALE == 1 ? 12 : 20;
    std::vector<int> a(N), b(M);
    for (auto& x : a) x = (int)rnd(rng, 0, 1);
    for (auto& x : b) x = (int)rnd(rng, 0, 1);
    std::vector<std::vector<int>> g(G, std::vector<int>(G));
    for (auto& r : g)
        for (auto& x : r) x = (int)rnd(rng, 0, 1);
    ll chk = 0;
    auto t0 = std::chrono::steady_clock::now();
    chk += flip_k_min(a, 3);
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    chk += flip_k_best(b).second;
    const ll t2 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    chk += lights_out(g).first;
    const ll t3 = elapsed_ms(t0);
    std::cout << "[flip_greedy] benchmark: flip_k_min(N = " << N << ") " << t1 << " ms / flip_k_best(N = " << M << ") " << t2 << " ms / lights_out("
              << G << " × " << G << ") " << t3 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace flip_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace flip_demo {

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
#define SHOW(expr, comment) flip_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    flip_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("反転（flip_greedy） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  flip_k_min / flip_k_starts / flip_k_best / add_k_mod_min\n"
              << "  flip_rect_min / flip_rect_plan / lights_out\n";

    api("【1】 flip_k_min(a, K, target) / flip_k_starts(a, K, target)", "int / optional<vector<int>>",
        "長さ K の区間をまとめて反転して全部 target にする最小回数（無理なら −1）/ 反転する区間の左端", "O(N)");
    RUN(const std::vector<int> a = {1, 1, 0, 1, 0, 1, 1});  // 1 = 後ろ向きの牛
    SHOW(flip_k_min(a, 3), "全部 0 にする");
    SHOW(flip_k_starts(a, 3), "[0, 3) [2, 5) [4, 7) を反転");
    SHOW(flip_k_min(a, 2), "K = 2 では無理");
    SHOW(flip_k_min(a, 3, 1), "全部 1 にするのは無理");

    api("【2】 flip_k_best(a, target)", "pair<int, int>", "{K, 回数}：回数が最小になる K（同じなら小さい K）。Face The Right Way の形", "O(N^2)");
    SHOW(flip_k_best(a), "");

    api("【3】 add_k_mod_min(a, K, m)", "long long", "長さ K の区間に +1（mod m）する操作で全部 0（mod m）にする最小回数（無理なら −1）", "O(N)");
    SHOW(add_k_mod_min(std::vector<long long>{2, 0, 0, 2}, 2, 3LL), "左から：+1 を [0,2) に 1 回、[1,3) に 2 回、[2,4) に 1 回");
    SHOW(add_k_mod_min(std::vector<long long>{1, 0, 0}, 2, 3LL), "最後が合わない → −1");

    api("【4】 flip_rect_min(g, h, w, target) / flip_rect_plan(g, h, w, target)", "long long / optional<vector<vector<int>>>",
        "h × w の長方形の反転で全部 target にする最小回数 / 反転する長方形の左上（press[i][j] = 1）", "O(HW)");
    RUN(const std::vector<std::vector<int>> g = {{1, 1, 0}, {1, 0, 1}, {0, 1, 1}});
    SHOW(flip_rect_min(g, 2, 2), "左上 (0,0) と (1,1) の 2 × 2");
    SHOW(flip_rect_plan(g, 2, 2), "");
    SHOW(flip_rect_min(g, 3, 2), "無理");

    api("【5】 lights_out(g, target)", "pair<int, vector<vector<int>>>",
        "マスを押すと自分と上下左右が反転。全部 target にする最小の押す回数と押し方（最小の中で辞書順最小）。Fliptile の形", "O(2^W H)");
    RUN(const std::vector<std::vector<int>> L = {{1, 0, 0, 1}, {0, 1, 1, 0}, {0, 1, 1, 0}, {1, 0, 0, 1}});
    SHOW(lights_out(L), "4 × 4");
    SHOW(lights_out(std::vector<std::vector<int>>{{1, 1, 0}, {1, 0, 0}, {0, 0, 0}}), "左上を 1 回押す");
}

}  // namespace flip_demo

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
    flip_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    flip_test::benchmark();    // 速度計測
    flip_demo::demo();         // 使い方デモ
    return 0;
}
