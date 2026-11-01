// =============================================================================
//  周期検出（状態が 1 手ずつ決まるシミュレーションで、K 手後の状態・K 手分の和を答える）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  状態 s から次の状態 next(s) が 1 つに決まるとき、いつか前に出た状態に戻り、そこからは同じ並びを繰り返す。
//    s0 → s1 → … → s_μ → … → s_(μ+λ−1) → s_μ → …（μ = 周期に入るまでの手数、λ = 周期の長さ）
//  状態は map に入れられる型なら何でもよい（整数・pair・vector・string・盤面など）
//    例: ABC179 E（Sequence Sum）… A_(n+1) = A_n² mod M の先頭 N 項の和（N ≤ 10^10）
//        「盤面を 10^18 回動かした後」「操作を K 回した後の位置」
//
//  機能一覧（k は手数 0 以上。10^18 でもよい）                                       計算量（L = μ + λ）
//    CycleFinder<S> cf(s0, next)          s0 から next を繰り返して周期を見つける             O(L log L)
//    CycleFinder<S, Map> cf(s0, next)     Map = unordered_map<S, int> なども使える（速い）    O(L)
//    preperiod() / period()               μ（周期に入るまでの手数）/ λ（周期の長さ）           O(1)
//    distinct_count()                     出てくる状態の種類 L = μ + λ
//    at(k)                                k 手後の状態（at(0) = s0）                           O(1)
//    index(k)                             k 手後の状態が states() の何番目か（0〜L−1）         O(1)
//    states()                             s0, s1, …, s_(L−1)（どれも違う状態）
//    first_time(s)                        初めて状態 s になる手数（ならなければ −1）           O(log L)
//    visit_count(s, k)                    0〜k−1 手後のうち状態が s になる回数                 O(log L)
//    make_sum<T>(value)                   value(状態) の和を答える表 CycleSum<T> を作る        O(L)
//      .sum(k)                            value(s0) + … + value(s_(k−1))（先頭 k 個）          O(1)
//      .sum(l, r) / .sum_closed(l, r)     l〜r−1 手後 / l〜r 手後 の value の和
//   メモリを使わない版（状態が大きくて全部は持てないとき）
//    brent_cycle(s0, next)                {μ, λ} を Brent の方法で（状態は == で比べる）       O(L) 時間・O(1) メモリ
//    state_after(s0, next, k, mu, lam)    k 手後の状態を、周期を使って min(k, μ + (k−μ) mod λ) 手だけ動かして求める
//
//  注意
//    next は同じ状態に対していつも同じ状態を返すこと（乱数や外の変数で変わると壊れる）
//    L が limit（既定 10^7）を超えたら止める。CycleFinder は L 個の状態を 2 回（states と map）持つので、状態が大きいときは注意
//    整数の状態 0〜N−1 で、いろいろな出発点について何度も聞くなら graph_grid_v2/06_functional_graph.cpp の方が速い
//    sum の T は long long・mint など。T(k) × T で掛け算する（k ≤ 10^18。long long であふれないかは自分で確認）
//  速さの目安: L = 3.3×10^5（long long の状態）で std::map 約 0.3 秒、unordered_map 約 0.05 秒、brent_cycle 約 0.003 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    cycle_test::random_test()  … 1 手ずつ動かす愚直解との突き合わせ
//    cycle_test::benchmark()    … 速度計測
//    cycle_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace cycle_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace cycle_internal

// value の和を答える表（CycleFinder::make_sum で作る）
template <class T>
struct CycleSum {
    long long mu = 0, lam = 1;
    std::vector<T> pre;  // pre[i] = value(s0) + … + value(s_(i−1))  (i = 0〜L)
    // 先頭 k 個（0〜k−1 手後）の和
    T sum(long long k) const {
        if (k < 0) cycle_internal::fail("CycleSum::sum", "k は 0 以上");
        const long long L = mu + lam;
        if (k <= L) return pre[(size_t)k];
        const long long full = (k - mu) / lam, rem = (k - mu) % lam;
        const T cyc = pre[(size_t)L] - pre[(size_t)mu];
        return pre[(size_t)(mu + rem)] + T(full) * cyc;
    }
    T sum(long long l, long long r) const {  // l〜r−1 手後
        if (l < 0 || r < l) cycle_internal::fail("CycleSum::sum", "0 ≤ l ≤ r にする");
        return sum(r) - sum(l);
    }
    T sum_closed(long long l, long long r) const {  // l〜r 手後
        if (l < 0 || r + 1 < l) cycle_internal::fail("CycleSum::sum_closed", "0 ≤ l ≤ r + 1 にする");
        return sum(r + 1) - sum(l);
    }
};

template <class S, class Map = std::map<S, int>>
struct CycleFinder {
    template <class Next>
    CycleFinder(const S& s0, Next next, long long limit = 10000000) {
        S cur = s0;
        while (true) {
            auto it = id_.find(cur);
            if (it != id_.end()) {
                mu_ = it->second;
                lam_ = (long long)seq_.size() - mu_;
                break;
            }
            if ((long long)seq_.size() >= limit) cycle_internal::fail("CycleFinder", "limit 手までに同じ状態に戻らなかった");
            id_.emplace(cur, (int)seq_.size());
            seq_.push_back(cur);
            cur = next(static_cast<const S&>(seq_.back()));
        }
    }
    long long preperiod() const { return mu_; }
    long long period() const { return lam_; }
    long long distinct_count() const { return (long long)seq_.size(); }
    const std::vector<S>& states() const { return seq_; }
    size_t index(long long k) const {
        if (k < 0) cycle_internal::fail("CycleFinder::index", "k は 0 以上");
        if (k < mu_) return (size_t)k;
        return (size_t)(mu_ + (k - mu_) % lam_);
    }
    const S& at(long long k) const { return seq_[index(k)]; }
    long long first_time(const S& s) const {
        auto it = id_.find(s);
        return it == id_.end() ? -1 : (long long)it->second;
    }
    // 0〜k−1 手後のうち状態が s になる回数
    long long visit_count(const S& s, long long k) const {
        if (k < 0) cycle_internal::fail("CycleFinder::visit_count", "k は 0 以上");
        auto it = id_.find(s);
        if (it == id_.end()) return 0;
        const long long t = it->second;
        if (t >= k) return 0;
        if (t < mu_) return 1;
        return (k - 1 - t) / lam_ + 1;
    }
    template <class T, class F>
    CycleSum<T> make_sum(F value) const {
        CycleSum<T> cs;
        cs.mu = mu_, cs.lam = lam_;
        cs.pre.assign(seq_.size() + 1, T(0));
        for (size_t i = 0; i < seq_.size(); ++i) cs.pre[i + 1] = cs.pre[i] + T(value(seq_[i]));
        return cs;
    }

   private:
    std::vector<S> seq_;
    Map id_;
    long long mu_ = 0, lam_ = 1;
};

// Brent の方法：{μ, λ}（メモリ O(1)。状態を == で比べる）
template <class S, class Next>
std::pair<long long, long long> brent_cycle(const S& s0, Next next, long long limit = 4000000000000000000LL) {
    long long power = 1, lam = 1;
    S tortoise = s0, hare = next(s0);
    while (!(tortoise == hare)) {
        if (power == lam) tortoise = hare, power *= 2, lam = 0;
        hare = next(hare);
        if (++lam > limit) cycle_internal::fail("brent_cycle", "limit 手までに同じ状態に戻らなかった");
    }
    long long mu = 0;
    tortoise = s0, hare = s0;
    for (long long i = 0; i < lam; ++i) hare = next(hare);
    while (!(tortoise == hare)) tortoise = next(tortoise), hare = next(hare), ++mu;
    return {mu, lam};
}
// k 手後の状態（mu, lam は brent_cycle などで求めたもの）
template <class S, class Next>
S state_after(const S& s0, Next next, long long k, long long mu, long long lam) {
    if (k < 0 || lam <= 0) cycle_internal::fail("state_after", "k ≥ 0, lam ≥ 1 にする");
    const long long steps = k < mu ? k : mu + (k - mu) % lam;
    S s = s0;
    for (long long i = 0; i < steps; ++i) s = next(s);
    return s;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace cycle_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define CY_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[search_cycle_detect] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 関数グラフ f（0〜n−1）で、いろいろな出発点を試す
template <class Map>
static void test_int(std::mt19937_64& rng, long long& cases) {
    const int n = (int)rnd(rng, 1, 40);
    std::vector<int> f(n);
    for (auto& x : f) x = (int)rnd(rng, 0, n - 1);
    if (rnd(rng, 0, 3) == 0)  // 長い道 + 閉路になりやすく
        for (int i = 0; i < n; ++i) f[i] = i + 1 < n ? i + 1 : (int)rnd(rng, 0, n - 1);
    std::vector<ll> val(n);
    for (auto& v : val) v = rnd(rng, -100, 100);
    const int s0 = (int)rnd(rng, 0, n - 1);
    auto nx = [&](const int& x) { return f[x]; };
    CycleFinder<int, Map> cf(s0, nx);
    const auto cs = cf.template make_sum<ll>([&](int x) { return val[x]; });
    const auto br = brent_cycle(s0, nx);
    CY_CHECK(br.first == cf.preperiod() && br.second == cf.period());
    CY_CHECK(cf.distinct_count() == cf.preperiod() + cf.period());
    // 愚直に 300 手
    std::vector<int> sim{s0};
    for (int k = 0; k < 300; ++k) sim.push_back(f[sim.back()]);
    std::vector<ll> pre{0};
    for (int k = 0; k < 300; ++k) pre.push_back(pre.back() + val[sim[k]]);
    // μ, λ の定義どおりか
    {
        std::map<int, int> first;
        int mu = -1, lam = -1;
        for (int k = 0; k <= 300; ++k) {
            if (first.count(sim[k])) {
                mu = first[sim[k]], lam = k - mu;
                break;
            }
            first[sim[k]] = k;
        }
        CY_CHECK(mu == cf.preperiod() && lam == cf.period());
        for (int x = 0; x < n; ++x) CY_CHECK(cf.first_time(x) == (first.count(x) ? first[x] : -1));
    }
    for (int k = 0; k < 300; ++k) {
        CY_CHECK(cf.at(k) == sim[k]);
        CY_CHECK(cs.sum(k) == pre[k]);
        CY_CHECK(state_after(s0, nx, k, br.first, br.second) == sim[k]);
        ++cases;
    }
    for (int q = 0; q < 30; ++q) {
        const ll l = rnd(rng, 0, 299), r = rnd(rng, l, 300);
        CY_CHECK(cs.sum(l, r) == pre[r] - pre[l]);
        if (r < 300) CY_CHECK(cs.sum_closed(l, r) == pre[r + 1] - pre[l]);
        const int x = (int)rnd(rng, 0, n - 1);
        const ll k = rnd(rng, 0, 300);
        CY_CHECK(cf.visit_count(x, k) == (ll)std::count(sim.begin(), sim.begin() + k, x));
        cases += 3;
    }
    // 大きい k：周期で割った位置と同じ
    for (int q = 0; q < 10; ++q) {
        const ll k = rnd(rng, 0, rnd(rng, 0, 1) ? 1000000000000LL : 1000000000000000000LL), mu = cf.preperiod(), lam = cf.period();
        const ll k2 = k < mu ? k : mu + (k - mu) % lam;
        CY_CHECK(cf.at(k) == cf.at(k2));
        // 和：k = k2 + t λ なら sum(k) = sum(k2) + t × (1 周の和)
        if (k >= mu && k <= 1000000000000LL) {  // 和は |値| ≤ 100・L ≤ 40 なので 10^12 手ならあふれない
            const ll t = (k - k2) / lam;
            const ll cyc = cs.sum(mu + lam) - cs.sum(mu);
            CY_CHECK(cs.sum(k) == cs.sum(k2) + t * cyc);
        }
        ++cases;
    }
}

// 状態が vector（置換を何回もかける：盤面の形）
static void test_vector(std::mt19937_64& rng, long long& cases) {
    const int n = (int)rnd(rng, 1, 7);
    std::vector<int> p(n);
    for (int i = 0; i < n; ++i) p[i] = i;
    std::shuffle(p.begin(), p.end(), rng);
    const int m = (int)rnd(rng, 1, 3);
    auto nx = [&](const std::vector<int>& s) {  // 並べ替えてから先頭を +1 mod m（前周期ができることもある）
        std::vector<int> t(n);
        for (int i = 0; i < n; ++i) t[i] = s[p[i]];
        t[0] = std::min(m - 1, t[0] + 1);
        return t;
    };
    std::vector<int> s0(n);
    for (auto& x : s0) x = (int)rnd(rng, 0, m - 1);
    CycleFinder<std::vector<int>> cf(s0, nx);
    const auto br = brent_cycle(s0, nx);
    CY_CHECK(br.first == cf.preperiod() && br.second == cf.period());
    std::vector<int> cur = s0;
    for (int k = 0; k < 200; ++k) {
        CY_CHECK(cf.at(k) == cur);
        cur = nx(cur);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 300 * TEST_SCALE; ++it) test_int<std::map<int, int>>(rng, cases);
    for (int it = 0; it < 300 * TEST_SCALE; ++it) test_int<std::unordered_map<int, int>>(rng, cases);
    for (int it = 0; it < 300 * TEST_SCALE; ++it) test_vector(rng, cases);
    // pair の状態・mint 風の型でも和が取れる
    {
        CycleFinder<std::pair<int, int>> cf({0, 1}, [](const std::pair<int, int>& s) { return std::make_pair(s.second, (s.first + s.second) % 10); });
        CY_CHECK(cf.preperiod() == 0 && cf.period() == 60);  // フィボナッチ mod 10 の周期（ピサノ周期）
        const auto cs = cf.make_sum<ll>([](const std::pair<int, int>& s) { return s.first; });
        ll brute = 0;
        std::pair<int, int> s{0, 1};
        for (int k = 0; k < 1000; ++k) brute += s.first, s = {s.second, (s.first + s.second) % 10};
        CY_CHECK(cs.sum(1000) == brute);
        cases += 2;
    }
    std::cout << "[search_cycle_detect] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[search_cycle_detect] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const ll M = TEST_SCALE == 1 ? 100003 : 1000003;  // x → 3x + 1 mod M（L は 33 万くらい）
    auto nx = [M](const ll& x) { return (x * 3 + 1) % M; };
    auto t0 = std::chrono::steady_clock::now();
    CycleFinder<ll> cf(1, nx);
    const auto cs = cf.make_sum<ll>([](ll x) { return x; });
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    CycleFinder<ll, std::unordered_map<ll, int>> cf2(1, nx);
    const ll t2 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto br = brent_cycle(1LL, nx);
    const ll t3 = elapsed_ms(t0);
    std::cout << "[search_cycle_detect] benchmark (L = " << cf.distinct_count() << "): map " << t1 << " ms / unordered_map " << t2
              << " ms / brent_cycle " << t3 << " ms  (check " << (cs.sum(1000000000000LL) + cf2.period() + br.second) % 1000 << ")\n";
}

}  // namespace cycle_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace cycle_demo {

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
#define SHOW(expr, comment) cycle_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    cycle_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("周期検出（search_cycle_detect） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  CycleFinder<S, Map = map<S, int>> cf(s0, next)\n"
              << "    preperiod / period / distinct_count / at / index / states / first_time / visit_count / make_sum<T>(value)\n"
              << "  CycleSum<T>: sum(k) / sum(l, r) / sum_closed(l, r)\n"
              << "  brent_cycle(s0, next) / state_after(s0, next, k, mu, lam)\n";

    api("【1】 CycleFinder<S> cf(s0, next)", "—", "s0 から next を繰り返し、前に出た状態に戻るまで動かす", "O(L log L)（L = μ + λ）");
    code("// x → x² mod 1001 を 2 から（ABC179 E の数列）");
    RUN(CycleFinder<long long> cf(2, [](long long x) { return x * x % 1001; }));
    SHOW(cf.states(), "s0, s1, …（どれも違う状態）");
    SHOW(cf.preperiod(), "μ：周期に入るまでの手数");
    SHOW(cf.period(), "λ：周期の長さ");
    SHOW(cf.distinct_count(), "出てくる状態の種類 μ + λ");

    api("【2】 at(k) / index(k)", "const S& / size_t", "k 手後の状態 / それが states() の何番目か（k は 10^18 でもよい）", "O(1)");
    SHOW(cf.at(0), "0 手後 = s0");
    SHOW(cf.at(5), "5 手後");
    SHOW(cf.at(1000000000000000000LL), "10^18 手後");
    SHOW(cf.index(1000000000000000000LL), "states() の何番目か");

    api("【3】 first_time(s) / visit_count(s, k)", "long long", "初めて s になる手数（ならなければ −1）/ 0〜k−1 手後に s になる回数", "O(log L)");
    SHOW(cf.first_time(256), "");
    SHOW(cf.first_time(3), "3 にはならない");
    SHOW(cf.visit_count(cf.at(cf.preperiod()), 100), "周期の中の状態は λ 手ごとに来る");

    api("【4】 make_sum<T>(value) → .sum(k) / .sum(l, r) / .sum_closed(l, r)", "CycleSum<T> / T",
        "value(状態) の和。sum(k) は先頭 k 個（0〜k−1 手後）、sum(l, r) は l〜r−1 手後", "O(L) / O(1)");
    RUN(const auto cs = cf.make_sum<long long>([](long long x) { return x; }));
    SHOW(cs.sum(6), "2 + 4 + 16 + 256 + 471 + 620（ABC179 E で N = 6, X = 2, M = 1001）");
    SHOW(cs.sum(10000000000LL), "N = 10^10 でもすぐ");
    SHOW(cs.sum(1, 3), "1, 2 手後 = 4 + 16");
    SHOW(cs.sum_closed(1, 3), "1〜3 手後 = 4 + 16 + 256");

    api("【5】 状態が盤面（vector<string>）でも使える", "—", "map に入れられる型なら何でも。速くしたいなら Map に unordered_map とハッシュを渡す", "—");
    {
        code("// 1 行の盤面を右へ 1 つ回し、先頭 2 マスが両方 'o' なら先頭を消す（自作の例）");
        RUN(auto step = [](const std::string& s) { std::string t = s.back() + s.substr(0, s.size() - 1); if (t[0] == 'o' && t[1] == 'o') t[0] = '.'; return t; });
        RUN(CycleFinder<std::string> board("oo..o", step));
        SHOW(board.states(), "");
        SHOW(board.preperiod(), "");
        SHOW(board.period(), "");
        SHOW(board.at(1000000000000000000LL), "10^18 手後の盤面");
    }

    api("【6】 brent_cycle(s0, next) / state_after(s0, next, k, mu, lam)", "pair<long long, long long> / S",
        "状態を持たずに {μ, λ} を求める（状態が大きくて全部持てないとき）/ それを使って k 手後を動かして求める", "O(L) 時間・O(1) メモリ");
    RUN(auto f = [](long long x) { return x * x % 1001; });
    RUN(const auto mt = brent_cycle(2LL, f));
    SHOW(mt, "{μ, λ}（CycleFinder と同じ）");
    SHOW(state_after(2LL, f, 1000000000000000000LL, mt.first, mt.second), "10^18 手後（μ + λ 手以内だけ動かす）");
}

}  // namespace cycle_demo

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
    cycle_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    cycle_test::benchmark();    // 速度計測
    cycle_demo::demo();         // 使い方デモ
    return 0;
}
