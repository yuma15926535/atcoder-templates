// =============================================================================
//  ダブリング（K 回たどる・値の集約・二分探索）  Doubling / doubling        --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 06 特殊グラフ
//  依存       : なし（単独で貼れる）
//  頂点番号   : 0-indexed。nxt[v] = v の次の頂点（-1 なら先がない = そこで止まる）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  v → nxt[v] を K 回（K ≤ 10^18）たどる。1 回の移動ごとの値 val[v] を、モノイド op で移動の順に集約できる
//  （和・最小・最大・一次関数の合成など。可換でなくてよい）。さらに「条件を満たす最初の回数」を二分探索で求める。
//    - 06_functional_graph との使い分け: K 回後の頂点・重みの和だけなら functional_graph（表を作らないので速い・省メモリ）。
//      最小・最大・非可換の集約、先がない頂点（-1）、「何回目に初めて〜になるか」はこちら
//    - 表は LOG × N（max_k = 10^18 なら LOG = 60。N = 2×10^5 で値が long long なら約 140 MB）。K の上限が小さいなら max_k を小さくする
//    - 旧テンプレートの Doubling（first_true つき）の機能はすべて入っている（旧と同時に貼らない。名前が同じ）
//
//  機能一覧                                                                       計算量
//    Doubling D(nxt, val, op, e, max_k = 10^18)   値つき（op(先, 後) で集約、e は単位元）         O(N LOG)
//    auto D = doubling(nxt, max_k = 10^18)        値なし（行き先だけ）                             O(N LOG)
//    D.kth(v, k)           -> int      ちょうど k 回進んだ頂点（途中で先がなければ -1）           O(LOG)
//    D.fold(v, k)          -> S        k 回分の値の集約（途中で止まるならそこまで）               O(LOG)
//    D.walk(v, k)          -> Walk     最大 k 回進む: {to（着いた頂点）, steps（進んだ回数）, value}  O(LOG)
//    D.max_steps(v, limit, ok)  -> Walk   ok(頂点, 値) が true のまま進める最大の回数（ok は true → false と変わる）  O(LOG)
//    D.first_true(v, limit, f)  -> Walk   f(頂点, 値) が初めて true になる回数（f は false → true。なければ steps = -1）  O(LOG)
//    k・limit は 0 以上 2^LOG 未満（max_k 以下ならよい。超えると abort）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    db_test::random_test() / db_test::benchmark() / db_demo::demo()
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
#ifndef GG_DOUBLING
#define GG_DOUBLING
// ---- ダブリング（グラフ・グリッド第二世代）  v → nxt[v] を K 回たどる。移動の値の集約（モノイド）と二分探索つき ----
// nxt[v] = -1 は「先がない（そこで止まる）」。val[v] = v → nxt[v] の 1 回の移動の値。
// op(a, b) は結合的（可換でなくてよい。移動した順に op(先の分, 後の分) で畳み込む）。e は op の単位元。
// 表は LOG × N（LOG = max_k のビット数。max_k = 10^18 なら 60）。K 回後の頂点・重みの和だけなら 06_functional_graph のほうが速く省メモリ
// 必要なヘッダ: <vector>, <cstdlib>
template <class S, class Op>
struct Doubling {
    struct Walk {
        int to = -1;          // 着いた頂点
        long long steps = 0;  // 実際に進んだ回数（先がなくなって止まると、頼んだ回数より少なくなる）
        S value;              // 進んだ分の値の集約
    };
    int n = 0, lg = 0;
    long long max_k = 0;
    std::vector<int> to_;  // to_[k * n + v] = v から 2^k 回進んだ頂点（進めなければ -1）
    std::vector<S> val_;   // val_[k * n + v] = その 2^k 回分の値
    Op op;
    S e;
    // max_k = 一度に進む回数の上限（これ以下の k だけ使える）
    Doubling(const std::vector<int>& nxt, const std::vector<S>& val, Op op_, typename std::common_type<S>::type e_, long long max_k_ = 1000000000000000000LL)
        : n((int)nxt.size()), max_k(max_k_), op(op_), e(e_) {
        while (lg < 62 && (1LL << lg) <= max_k) ++lg;
        if (lg == 0) lg = 1;
        to_.assign((size_t)lg * n, -1), val_.assign((size_t)lg * n, e);
        for (int v = 0; v < n; ++v) to_[v] = nxt[v], val_[v] = (nxt[v] < 0 ? e : val[v]);
        for (int k = 1; k < lg; ++k) {
            const int* pt = &to_[(size_t)(k - 1) * n];
            const S* pv = &val_[(size_t)(k - 1) * n];
            int* ct = &to_[(size_t)k * n];
            S* cv = &val_[(size_t)k * n];
            for (int v = 0; v < n; ++v) {
                const int mid = pt[v];
                if (mid < 0 || pt[mid] < 0) continue;  // 2^k 回は進めない
                ct[v] = pt[mid], cv[v] = op(pv[v], pv[mid]);
            }
        }
    }
    // v から最大 k 回進む（先がなくなったらそこで止まる）  O(LOG)
    Walk walk(int v, long long k) const {
        check(k);
        Walk w{v, 0, e};
        for (int b = lg - 1; b >= 0; --b) {
            if (!((k >> b) & 1)) continue;
            const size_t i = (size_t)b * n + w.to;
            if (to_[i] >= 0) {
                w.value = op(w.value, val_[i]), w.to = to_[i], w.steps += 1LL << b;
            } else {  // 2^b 回は進めない: 細かく刻んで、進めるところまで進む
                for (int c = b - 1; c >= 0; --c) {
                    const size_t j = (size_t)c * n + w.to;
                    if (to_[j] >= 0) w.value = op(w.value, val_[j]), w.to = to_[j], w.steps += 1LL << c;
                }
                return w;
            }
        }
        return w;
    }
    int kth(int v, long long k) const {  // ちょうど k 回進んだ頂点（途中で先がなくなれば -1）
        const Walk w = walk(v, k);
        return w.steps == k ? w.to : -1;
    }
    S fold(int v, long long k) const { return walk(v, k).value; }  // k 回分の値（途中で止まるならそこまで）
    // ok(頂点, 値) が、進んだ回数 m について「m が小さいうちは true、ある所から先は false」のとき、
    // ok が true のままでいられる最大の m（limit 以下・先がなくなるまで）。m = 0（ok(v, e)）は true とみなす  O(LOG)
    template <class F>
    Walk max_steps(int v, long long limit, F ok) const {
        check(limit);
        Walk w{v, 0, e};
        for (int b = lg - 1; b >= 0; --b) {
            if (w.steps + (1LL << b) > limit) continue;
            const size_t i = (size_t)b * n + w.to;
            if (to_[i] < 0) continue;
            const S nv = op(w.value, val_[i]);
            if (ok(to_[i], nv)) w.value = nv, w.to = to_[i], w.steps += 1LL << b;
        }
        return w;
    }
    // f(頂点, 値) が「false, ..., false, true, ..., true」と変わるとき、初めて true になる回数 m（0 ≤ m ≤ limit）。
    // limit 回以内（先がなくなるまで）に true にならなければ steps = -1  O(LOG)
    template <class F>
    Walk first_true(int v, long long limit, F f) const {
        if (f(v, e)) return Walk{v, 0, e};
        Walk w = max_steps(v, limit, [&](int u, const S& x) { return !f(u, x); });
        if (w.steps < limit && to_[w.to] >= 0) {
            const S nv = op(w.value, val_[w.to]);
            if (f(to_[w.to], nv)) return Walk{to_[w.to], w.steps + 1, nv};
        }
        return Walk{-1, -1, e};
    }

  private:
    void check(long long k) const {
        if (k < 0 || (k >> lg) != 0) std::abort();  // k は 0 以上 2^LOG 未満（max_k を大きくして作り直す）
    }
};
namespace gg_internal {
struct NoValue {};
struct NoValueOp {
    NoValue operator()(NoValue, NoValue) const { return {}; }
};
}  // namespace gg_internal
// 値を使わない（行き先だけの）ダブリング: auto D = doubling(nxt); D.kth(v, K) / D.walk(v, K).steps / D.first_true(v, K, f)
inline Doubling<gg_internal::NoValue, gg_internal::NoValueOp> doubling(const std::vector<int>& nxt, long long max_k = 1000000000000000000LL) {
    return Doubling<gg_internal::NoValue, gg_internal::NoValueOp>(nxt, std::vector<gg_internal::NoValue>(nxt.size()), gg_internal::NoValueOp(), gg_internal::NoValue(), max_k);
}
#endif  // GG_DOUBLING
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "doubling";
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

namespace db_test {

static const long long MOD = 998244353;
struct Aff {  // x → a x + b（mod）。op(f, g) = 「f のあと g」
    long long a, b;
    bool operator==(const Aff& o) const { return a == o.a && b == o.b; }
};

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 10);
        std::vector<int> nxt(n);
        for (int v = 0; v < n; ++v) nxt[v] = (it % 3 == 0 && rng() % 4 == 0) ? -1 : rnd(rng, 0, n - 1);
        std::vector<long long> w(n);
        std::vector<Aff> f(n);
        for (int v = 0; v < n; ++v) w[v] = rnd(rng, 0, 9), f[v] = {rnd(rng, 1, 5), rnd(rng, 0, 5)};
        const long long max_k = (it % 2) ? 1000000000000000000LL : rnd(rng, 0, 40);
        auto mn = [](long long a, long long b) { return std::min(a, b); };
        auto comp = [](Aff p, Aff q) { return Aff{q.a * p.a % MOD, (q.a * p.b + q.b) % MOD}; };
        const Doubling Dsum(nxt, w, [](long long a, long long b) { return a + b; }, 0, max_k);
        const Doubling Dmin(nxt, w, mn, (long long)4e18, max_k);
        const Doubling Daff(nxt, f, comp, Aff{1, 0}, max_k);
        const auto Dp = doubling(nxt, max_k);
        for (int q = 0; q < 30; ++q) {
            const int v = rnd(rng, 0, n - 1);
            const long long k = (rng() % 3 == 0) ? (long long)(rng() % (unsigned long long)(max_k + 1)) : std::min<long long>(max_k, rnd(rng, 0, 30));
            // 愚直（k が大きいときはサイクルで同じことの繰り返しになるので、30 + n 回で打ち切って比べられる範囲だけ）
            if (k <= 60) {
                int x = v;
                long long steps = 0, s = 0, m = (long long)4e18;
                Aff af{1, 0};
                for (; steps < k && nxt[x] >= 0; ++steps) s += w[x], m = std::min(m, w[x]), af = comp(af, f[x]), x = nxt[x];
                const auto ws = Dsum.walk(v, k);
                CHECK(ws.to == x && ws.steps == steps && ws.value == s);
                CHECK(Dmin.fold(v, k) == m);
                CHECK(Daff.fold(v, k) == af);
                CHECK(Dsum.kth(v, k) == (steps == k ? x : -1));
                CHECK(Dp.kth(v, k) == (steps == k ? x : -1) && Dp.walk(v, k).steps == steps);
                // max_steps: 和が X 以下でいられる最大の回数
                const long long X = rnd(rng, 0, 40);
                const auto ms = Dsum.max_steps(v, k, [&](int, long long val) { return val <= X; });
                {
                    int y = v;
                    long long st = 0, sm = 0;
                    while (st < k && nxt[y] >= 0 && sm + w[y] <= X) sm += w[y], y = nxt[y], ++st;
                    CHECK(ms.steps == st && ms.to == y && ms.value == sm);
                }
                // first_true: 頂点 t に初めて着く回数（頂点の条件は単調ではないので、「t に着いたことがある」を値で持つ）
                const int tg = rnd(rng, 0, n - 1);
                std::vector<long long> hit(n);
                for (int u = 0; u < n; ++u) hit[u] = (nxt[u] == tg);
                const Doubling Dor(nxt, hit, [](long long a, long long b) { return a | b; }, 0, max_k);
                const auto ft = Dor.first_true(v, k, [&](int u, long long val) { return u == tg || val != 0; });
                {
                    int y = v;
                    long long st = 0, ans = -1;
                    while (true) {
                        if (y == tg) { ans = st; break; }
                        if (st >= k || nxt[y] < 0) break;
                        y = nxt[y], ++st;
                    }
                    CHECK(ft.steps == ans);
                    if (ans >= 0) CHECK(ft.to == tg);
                }
                // 値なし版の first_true（いつまでも false なら見つからない）
                const auto fp = Dp.first_true(v, k, [&](int, gg_internal::NoValue) { return false; });
                CHECK(fp.steps == -1 && fp.to == -1);
            } else {
                // 大きい k: 結果は「k mod 周期」で決まる。walk の steps と kth の整合だけ見る
                const auto ws = Dsum.walk(v, k);
                CHECK(ws.steps <= k);
                CHECK(Dsum.kth(v, k) == (ws.steps == k ? ws.to : -1));
                CHECK(Dp.walk(v, k).to == ws.to);
                // k = a + b と分けて進んでも同じ
                const long long a = (long long)(rng() % (unsigned long long)(k + 1));
                const auto w1 = Dsum.walk(v, a);
                if (w1.steps == a) {
                    const auto w2 = Dsum.walk(w1.to, k - a);
                    CHECK(w2.to == ws.to && w1.steps + w2.steps == ws.steps && w1.value + w2.value == ws.value);
                }
            }
        }
    }
    std::cout << "[doubling] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[doubling] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    const int n = 200000, q = 200000;
    std::vector<int> nxt(n);
    std::vector<long long> w(n);
    for (int v = 0; v < n; ++v) nxt[v] = (int)(rng() % n), w[v] = (long long)(rng() % 1000000000);
    auto t0 = std::chrono::steady_clock::now();
    const Doubling D(nxt, w, [](long long a, long long b) { return std::min(a, b); }, (long long)4e18);
    const long long t1 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    long long chk = 0;
    for (int i = 0; i < q; ++i) chk += D.fold((int)(rng() % n), (long long)(rng() % 1000000000000000000ULL)) % 7;
    const long long t2 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    const auto P = doubling(nxt);
    const long long t3 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < q; ++i) chk += P.kth((int)(rng() % n), (long long)(rng() % 1000000000000000000ULL)) % 7;
    std::cout << "[doubling] benchmark: N = " << n << "  作る（最小値・10^18）" << t1 << " ms, fold " << q << " 回 " << t2 << " ms, 値なしを作る " << t3 << " ms, kth " << q
              << " 回 " << ms_since(t0) << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace db_test


namespace db_demo {
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
    banner("ダブリング 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"Doubling D(nxt, val, op, e, max_k);", "", "nxt[v] = 次（-1 で止まる）、val[v] = その移動の値"},
        {"auto D = doubling(nxt, max_k);", "", "値なし（行き先だけ）"},
        {"D.kth(v, k)", "int", "ちょうど k 回後の頂点（進めなければ -1）"},
        {"D.fold(v, k)", "S", "k 回分の値を移動の順に op でまとめたもの"},
        {"D.walk(v, k)", "Walk", "{to, steps, value}: 最大 k 回（止まるまで）進む"},
        {"D.max_steps(v, limit, ok)", "Walk", "ok(頂点, 値) が true のまま進める最大の回数"},
        {"D.first_true(v, limit, f)", "Walk", "f(頂点, 値) が初めて true になる回数（なければ steps = -1）"},
    };
    cheat_sheet("API 一覧", rows, 40, 10);

    head("デモで使う移動: 0 → 1 → 2 → 3 → 1（1, 2, 3 が周期 3 のサイクル）、4 → -1（止まる）");
    RUN(std::vector<int> nxt = {1, 2, 3, 1, -1});
    RUN(std::vector<long long> w = {5, 2, 7, 1, 9});

    api("【1】 K 回移動したときの最小値（逆元がないので functional_graph の和ではできない）", "long long", "op = min、e = とても大きい値", "作る O(N LOG)、1 回 O(LOG)");
    RUN(Doubling D(nxt, w, [](long long a, long long b) { return std::min(a, b); }, (long long)4e18));
    SHOW(D.fold(0, 2), "0 → 1 → 2 の辺の値 5, 2 の最小");
    SHOW(D.kth(0, 1000000000000000000LL), "10^18 回後");
    SHOW(D.fold(2, 1), "");

    api("【2】 先がない頂点（-1）", "int / Walk", "kth は -1、walk は止まったところまで", "O(LOG)");
    SHOW(D.kth(4, 1), "4 の先はない");
    SHOW(D.walk(4, 5).steps, "1 回も進めない");
    SHOW(D.kth(4, 0), "0 回なら自分");

    api("【3】 何回目に初めて和が X を超えるか（first_true）", "Walk", "f(頂点, 値の和) が false → true と変わる。値は 0 以上なので和は増えていく", "O(LOG)");
    RUN(Doubling S(nxt, w, [](long long a, long long b) { return a + b; }, 0));
    code("auto r = S.first_true(0, 100, [](int, long long s) { return s > 20; });   // 5 + 2 + 7 + 1 + 2 + 7 = 24 > 20");
    auto r = S.first_true(0, 100, [](int, long long s) { return s > 20; });
    SHOW(r.steps, "6 回目");
    SHOW(r.to, "そのとき着いている頂点");
    SHOW(r.value, "そのときの和");
    code("auto m = S.max_steps(0, 100, [](int, long long s) { return s <= 20; });   // 和が 20 以下でいられる最大の回数");
    auto m = S.max_steps(0, 100, [](int, long long s) { return s <= 20; });
    SHOW(m.steps, "5 回（和 17）");

    api("【4】 非可換の値: 一次関数の合成（移動の順に f を重ねる）", "S", "op(f, g) = f のあと g（値の順が大事なときも正しい）", "O(LOG)");
    code("struct F { long long a, b; };   // x → a x + b");
    code("Doubling A(nxt, fs, [](F p, F q) { return F{q.a * p.a, q.a * p.b + q.b}; }, F{1, 0});");
    note("A.fold(v, K) が「v から K 回の移動の関数を順に合成したもの」になる（mod を取るなら op の中で取る）。");

    head("使い分け");
    note("K 回後の頂点・重みの和（逆元のある演算）だけなら 06_functional_graph（表を作らないので速い・省メモリ）。");
    note("木の上で親へ K 回（k 番目の祖先）は 05_lca の kth_ancestor。パス上の最小・最大は 05_hld。");
    note("max_k を問題の K の上限にすると表が小さくなる（K ≤ N なら max_k = N）。");
}

}  // namespace db_demo

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
    db_test::random_test();
    db_test::benchmark();
    db_demo::demo();
    return 0;
}
