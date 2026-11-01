// =============================================================================
//  ヒューリスティック（AHC）の型：速い乱数・タイマー・焼きなまし・ビームサーチ   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  AtCoder Heuristic Contest（AHC）でよく書く部品。焼きなましは「1 手変えたときのスコアの差 delta を計算 →
//  受け入れるなら確定、受け入れないなら元に戻す」を時間いっぱい繰り返す型にしてある（差分計算は自分で書く）
//
//  機能一覧                                                                                  計算量
//    Xorshift rng(seed)                 速い乱数（xorshift64）
//      rng() / rng.next_int(lo, hi)     64 ビットの乱数 / [lo, hi] の整数                         O(1)
//      rng.next_double() / rng.shuffle(v)   [0, 1) の実数 / vector を混ぜる
//    HeuristicTimer timer               作った時から測る。timer.sec() / timer.ms() / timer.reset()
//    annealing(time_limit, T0, T1, propose, accept, reject, maximize = true, seed)
//        time_limit 秒まで：delta = propose()（1 手を準備してスコアの差を返す）→ 受け入れるなら accept()、
//        受け入れないなら reject()。温度は T0 から T1 へ指数的に下げる（T = T0 (T1/T0)^進み具合）。
//        maximize = true なら delta ≥ 0 は必ず受け入れ、delta < 0 は確率 exp(delta / T) で受け入れる
//        返り値 AnnealStats {試した回数, 受け入れた回数}
//    beam_search(init, depth, width, expand, score)
//        各段で expand(状態, 子の vector) で子を作り、score（大きいほど良い）の上位 width 個を残す。depth 段後の一番良い状態
//
//  注意
//    一番良い状態は自分で覚えておく（accept の中で「今のスコア > 最高」なら保存する）
//    T0 は「受け入れてよい悪化の大きさ」くらい、T1 はその 1/100〜1/1000 くらいから調整する
//    時間は 256 回ごとに測る（propose が重いときは time_limit を少し短めに）
//  速さの目安: annealing の枠だけで 1 秒に 約 10^8 回（propose・accept が軽いとき）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    heuristic_test::random_test()  … 乱数の範囲、全探索で答えが分かる小さい問題で焼きなまし・ビームサーチを確認
//    heuristic_test::benchmark()    … 速度計測
//    heuristic_demo::demo()         … 全機能の使い方デモ（日本語の説明付き。巡回セールスマンの 2-opt）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
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
struct Xorshift {
    uint64_t x;
    explicit Xorshift(uint64_t seed = 88172645463325252ULL) : x(seed ? seed : 88172645463325252ULL) {}
    uint64_t operator()() {
        x ^= x << 7;
        x ^= x >> 9;
        return x;
    }
    long long next_int(long long lo, long long hi) {  // [lo, hi]
        const uint64_t w = (uint64_t)(hi - lo) + 1;
        return lo + (long long)(w == 0 ? (*this)() : (*this)() % w);
    }
    double next_double() { return (double)((*this)() >> 11) * (1.0 / 9007199254740992.0); }  // [0, 1)
    template <class T>
    void shuffle(std::vector<T>& v) {
        for (size_t i = v.size(); i > 1; --i) std::swap(v[i - 1], v[(size_t)next_int(0, (long long)i - 1)]);
    }
};

struct HeuristicTimer {
    std::chrono::steady_clock::time_point st = std::chrono::steady_clock::now();
    void reset() { st = std::chrono::steady_clock::now(); }
    double sec() const { return std::chrono::duration<double>(std::chrono::steady_clock::now() - st).count(); }
    double ms() const { return sec() * 1000.0; }
};

struct AnnealStats {
    long long iterations = 0, accepted = 0;
};
template <class Propose, class Accept, class Reject>
AnnealStats annealing(double time_limit, double T0, double T1, Propose propose, Accept accept, Reject reject, bool maximize = true,
                      uint64_t seed = 12345) {
    Xorshift rng(seed);
    HeuristicTimer timer;
    AnnealStats st;
    double T = T0, progress = 0;
    const double ratio = (T1 > 0 && T0 > 0) ? T1 / T0 : 1.0;
    while (true) {
        if ((st.iterations & 255) == 0) {
            progress = timer.sec() / time_limit;
            if (progress >= 1.0) break;
            T = T0 * std::pow(ratio, progress);
        }
        ++st.iterations;
        double delta = propose();
        if (!maximize) delta = -delta;
        if (delta >= 0 || (T > 0 && rng.next_double() < std::exp(delta / T))) {
            accept();
            ++st.accepted;
        } else {
            reject();
        }
    }
    return st;
}

template <class State, class Expand, class Score>
State beam_search(const State& init, int depth, int width, Expand expand, Score score) {
    std::vector<State> cur{init}, nxt;
    for (int d = 0; d < depth; ++d) {
        nxt.clear();
        for (const State& s : cur) expand(s, nxt);
        if (nxt.empty()) break;
        if ((int)nxt.size() > width) {
            std::vector<std::pair<double, int>> key(nxt.size());
            for (size_t i = 0; i < nxt.size(); ++i) key[i] = {-(double)score(nxt[i]), (int)i};
            std::nth_element(key.begin(), key.begin() + width, key.end());
            std::vector<State> keep;
            keep.reserve(width);
            for (int i = 0; i < width; ++i) keep.push_back(std::move(nxt[key[i].second]));
            nxt.swap(keep);
        }
        cur.swap(nxt);
    }
    size_t best = 0;
    for (size_t i = 1; i < cur.size(); ++i)
        if (score(cur[best]) < score(cur[i])) best = i;
    return cur[best];
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace heuristic_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define HE_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[tool_heuristic] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static double tour_len(const std::vector<int>& p, const std::vector<double>& x, const std::vector<double>& y) {
    double s = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        const int a = p[i], b = p[(i + 1) % p.size()];
        s += std::hypot(x[a] - x[b], y[a] - y[b]);
    }
    return s;
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    // 乱数の範囲
    {
        Xorshift xr(1);
        for (int i = 0; i < 100000; ++i) {
            const ll lo = rnd(rng, -10, 10), hi = lo + rnd(rng, 0, 20);
            const ll v = xr.next_int(lo, hi);
            HE_CHECK(lo <= v && v <= hi);
            const double d = xr.next_double();
            HE_CHECK(0 <= d && d < 1);
            ++cases;
        }
        std::vector<int> v(50);
        std::iota(v.begin(), v.end(), 0);
        xr.shuffle(v);
        std::vector<int> s = v;
        std::sort(s.begin(), s.end());
        for (int i = 0; i < 50; ++i) HE_CHECK(s[i] == i);
    }
    // 焼きなまし：8 都市の巡回セールスマン（全探索の最短と同じになるはず）。2-opt の差分
    for (int it = 0; it < 5 * TEST_SCALE; ++it) {
        const int n = 8;
        std::vector<double> x(n), y(n);
        for (int i = 0; i < n; ++i) x[i] = (double)rnd(rng, 0, 1000), y[i] = (double)rnd(rng, 0, 1000);
        std::vector<int> p(n);
        std::iota(p.begin(), p.end(), 0);
        double best = 1e18;
        {
            std::vector<int> q(p.begin() + 1, p.end());
            do {
                std::vector<int> t{0};
                t.insert(t.end(), q.begin(), q.end());
                best = std::min(best, tour_len(t, x, y));
            } while (std::next_permutation(q.begin(), q.end()));
        }
        Xorshift xr(7 + it);
        double cur = tour_len(p, x, y), found = cur;
        int a = 0, b = 0;
        auto dist = [&](int u, int v) { return std::hypot(x[u] - x[v], y[u] - y[v]); };
        const auto st = annealing(
            0.05, 300, 1,
            [&]() {  // 区間 [a, b] を反転したときの長さの変化
                a = (int)xr.next_int(1, n - 1), b = (int)xr.next_int(1, n - 1);
                if (a > b) std::swap(a, b);
                const int pa = p[a - 1], qa = p[a], qb = p[b], nb = p[(b + 1) % n];
                return dist(pa, qa) + dist(qb, nb) - dist(pa, qb) - dist(qa, nb);  // 短くなる分（大きいほど良い）
            },
            [&]() {
                std::reverse(p.begin() + a, p.begin() + b + 1);
                cur = tour_len(p, x, y);
                found = std::min(found, cur);
            },
            [&]() {});
        HE_CHECK(st.iterations > 1000 && st.accepted > 0);
        HE_CHECK(std::abs(found - best) < 1e-6);
        ++cases;
    }
    // 時間を守る
    {
        HeuristicTimer t;
        annealing(0.02, 1, 0.01, [] { return -1.0; }, [] {}, [] {});
        HE_CHECK(t.sec() < 0.2);
        ++cases;
    }
    // ビームサーチ：長さ 6 の 0〜3 の列、スコア = Σ c[i][a_i] − 隣が同じなら 5。幅を全部にすれば全探索と同じ
    for (int it = 0; it < 50 * TEST_SCALE; ++it) {
        const int D = 6, K = 4;
        std::vector<std::vector<int>> c(D, std::vector<int>(K));
        for (auto& r : c)
            for (auto& v : r) v = (int)rnd(rng, 0, 9);
        struct St {
            std::vector<int> a;
            int sc = 0;
        };
        auto expand = [&](const St& s, std::vector<St>& out) {
            const int i = (int)s.a.size();
            for (int v = 0; v < K; ++v) {
                St t = s;
                t.a.push_back(v);
                t.sc += c[i][v] - (i > 0 && s.a.back() == v ? 5 : 0);
                out.push_back(t);
            }
        };
        const St r = beam_search(St{}, D, 4096, expand, [](const St& s) { return s.sc; });
        int best = -1000;
        for (int m = 0; m < 4096; ++m) {
            int sc = 0, prev = -1, t = m;
            for (int i = 0; i < D; ++i) {
                const int v = t % K;
                t /= K;
                sc += c[i][v] - (v == prev ? 5 : 0);
                prev = v;
            }
            best = std::max(best, sc);
        }
        HE_CHECK(r.sc == best && (int)r.a.size() == D);
        const St r2 = beam_search(St{}, D, 3, expand, [](const St& s) { return s.sc; });
        HE_CHECK(r2.sc <= best);
        ++cases;
    }
    std::cout << "[tool_heuristic] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[tool_heuristic] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    long long sum = 0;
    Xorshift xr(3);
    const auto st = annealing(0.2, 10, 0.1, [&]() { return (double)((long long)(xr() & 15) - 8); }, [&]() { ++sum; }, [&]() {});
    std::cout << "[tool_heuristic] benchmark: annealing の枠 0.2 秒で " << st.iterations << " 回（受け入れ " << st.accepted << " 回）  (check "
              << sum % 1000 << ")\n";
}

}  // namespace heuristic_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace heuristic_demo {

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
#define SHOW(expr, comment) heuristic_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    heuristic_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("ヒューリスティック（tool_heuristic） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Xorshift rng(seed): rng() / next_int(lo, hi) / next_double() / shuffle(v)\n"
              << "  HeuristicTimer: sec() / ms() / reset()\n"
              << "  annealing(time_limit, T0, T1, propose, accept, reject, maximize, seed) → AnnealStats {iterations, accepted}\n"
              << "  beam_search(init, depth, width, expand, score)\n";

    api("【1】 Xorshift rng(seed)", "—", "速い乱数。next_int(lo, hi) は [lo, hi]、next_double() は [0, 1)", "O(1)");
    RUN(Xorshift rng(2026));
    SHOW(rng.next_int(1, 6), "さいころ");
    SHOW(rng.next_int(1, 6), "");
    RUN(std::vector<int> v = {1, 2, 3, 4, 5});
    RUN(rng.shuffle(v));
    SHOW(v, "");

    api("【2】 annealing(time_limit, T0, T1, propose, accept, reject)", "AnnealStats",
        "propose で 1 手を準備してスコアの差を返す → 受け入れたら accept（確定）、だめなら reject（戻す）", "時間いっぱい");
    {
        code("// 20 都市の巡回セールスマン：2-opt（区間の反転）で短くする。差 = 短くなる長さ");
        const int n = 20;
        Xorshift r(1);
        std::vector<double> x(n), y(n);
        for (int i = 0; i < n; ++i) x[i] = (double)r.next_int(0, 100), y[i] = (double)r.next_int(0, 100);
        std::vector<int> p(n);
        std::iota(p.begin(), p.end(), 0);
        auto d = [&](int a, int b) { return std::hypot(x[a] - x[b], y[a] - y[b]); };
        auto len = [&]() {
            double s = 0;
            for (int i = 0; i < n; ++i) s += d(p[i], p[(i + 1) % n]);
            return s;
        };
        const double before = len();
        int a = 0, b = 0;
        const auto st = annealing(
            0.1, 20, 0.1,
            [&]() {
                a = (int)r.next_int(1, n - 1), b = (int)r.next_int(1, n - 1);
                if (a > b) std::swap(a, b);
                return d(p[a - 1], p[a]) + d(p[b], p[(b + 1) % n]) - d(p[a - 1], p[b]) - d(p[a], p[(b + 1) % n]);
            },
            [&]() { std::reverse(p.begin() + a, p.begin() + b + 1); }, [&]() {});
        note("長さ " + to_s(before) + " → " + to_s(len()) + "（" + to_s(st.iterations) + " 回試して " + to_s(st.accepted) + " 回受け入れ）");
    }

    api("【3】 beam_search(init, depth, width, expand, score)", "State", "各段で子を作り、score の上位 width 個だけ残して depth 段進む", "O(depth × width × 子の数)");
    {
        code("// 長さ 5 の 0〜9 の列。スコア = 各桁の値の和 − 隣と同じなら 8（自作の例）");
        struct St {
            std::vector<int> a;
            int sc = 0;
        };
        auto expand = [](const St& s, std::vector<St>& out) {
            for (int dg = 0; dg < 10; ++dg) {
                St t = s;
                t.sc += dg - (!s.a.empty() && s.a.back() == dg ? 8 : 0);
                t.a.push_back(dg);
                out.push_back(t);
            }
        };
        const St r = beam_search(St{}, 5, 10, expand, [](const St& s) { return s.sc; });
        note("→ " + to_s(r.a) + "（スコア " + to_s(r.sc) + "）");
    }
}

}  // namespace heuristic_demo

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
    heuristic_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    heuristic_test::benchmark();    // 速度計測
    heuristic_demo::demo();         // 使い方デモ
    return 0;
}
