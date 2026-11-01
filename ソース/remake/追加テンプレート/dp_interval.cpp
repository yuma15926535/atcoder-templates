// =============================================================================
//  区間 DP（dp[l][r] を短い区間から埋める型・区間を 2 つに分ける型・Knuth 最適化）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  区間 [l, r)（半開区間。0 ≤ l ≤ r ≤ n）ごとの値 dp[l][r] を、長さの短い区間から順に求める
//    例: EDPC N（Slimes）… 隣どうしを合体させるコストの最小 → split（O(n^3)）/ split_knuth（O(n^2)）
//        行列の積の順番、回文にするための最小の挿入、区間を取り除く操作の最小回数 など
//
//  機能一覧（dp は (n + 1) × (n + 1) 個。区間は半開区間 [l, r)、閉区間は _closed 版）    計算量
//    IntervalDP<T> dp(n)                 長さ n の列の区間 DP（T は long long など）             O(n^2)
//    run(base, trans)                    一般形：dp[i][i] = base(i)（空の区間）、dp[l][r] = trans(l, r, dp)。
//                                        trans の中では dp.get(a, b) で [l, r) に含まれる、より短い区間の値を読める
//                                        （l の大きい順・r の小さい順に求める）                   O(n^2 × trans)
//    split(leaf, cost, is_min = true)    2 つに分ける形：dp[i][i+1] = leaf(i)、
//                                        dp[l][r] = min（is_min = false なら max）_{l<m<r} dp[l][m] + dp[m][r] + cost(l, m, r)  O(n^3)
//    split_knuth(leaf, w)                split の cost が w(l, r)（m によらない）で、w が四角不等式
//                                        w(a, c) + w(b, d) ≤ w(a, d) + w(b, c)（a ≤ b ≤ c ≤ d）と単調性 w(b, c) ≤ w(a, d) を
//                                        満たすときの最小化（例: w = 区間の和（0 以上）、区間の和の 2 乗）       O(n^2)
//    get(l, r) / get_closed(l, r)        dp[l][r]（[l, r)）/ [l, r]                               O(1)
//    all()                               dp[0][n]
//    split_point(l, r)                   split / split_knuth で [l, r) を分けた位置 m              O(1)
//    merge_order()                       [0, n) を分ける木の、合体の順番 (l, m, r) の一覧（短い方から）  O(n)
//    size()                              n
//
//  注意
//    メモリは (n + 1)^2 × T（split 系はさらに int の (n + 1)^2）。n = 5000 の long long で 約 200 MB（split 系は 約 300 MB）
//    run の trans で [l, r) に含まれない区間・同じ長さの区間を get すると止める（まだ求めていない値を読む誤りを防ぐ）
//    円環の区間 DP は、列を 2 周分（長さ 2n）つないで作り、長さ n の区間 get(i, i + n) の最小を取る
//    split で cost が m によらないなら cost(l, m, r) の m は無視してよい（例: [](int l, int, int r) { … }）
//  速さの目安: split（n = 400）約 0.01 秒、split_knuth（n = 5000）約 0.3 秒、run（n = 5000、trans が O(1)）約 0.13 秒（多くは表を作る時間）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    intervaldp_test::random_test()  … メモ化再帰の愚直解との突き合わせ
//    intervaldp_test::benchmark()    … 速度計測
//    intervaldp_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace intervaldp_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace intervaldp_internal

template <class T>
struct IntervalDP {
    explicit IntervalDP(int n = 0) : n_(n), dp_((size_t)(n + 1) * (n + 1), T{}) {
        if (n < 0) intervaldp_internal::fail("IntervalDP", "n は 0 以上");
    }
    int size() const { return n_; }

    // 一般形：dp[i][i] = base(i)、dp[l][r] = trans(l, r, *this)（長さの短い順）
    template <class Base, class Trans>
    void run(Base base, Trans trans) {
        for (int i = 0; i <= n_; ++i) at(i, i) = base(i);
        running_ = true;
        for (int l = n_ - 1; l >= 0; --l)  // [l, r) に含まれる短い区間は先に求まっている順番（メモリを順に読むので速い）
            for (int r = l + 1; r <= n_; ++r) {
                cur_l_ = l, cur_r_ = r;
                at(l, r) = trans(l, r, static_cast<const IntervalDP&>(*this));
            }
        running_ = false;
        arg_.clear();
    }
    // 2 つに分ける形：dp[i][i+1] = leaf(i)、dp[l][r] = best_m dp[l][m] + dp[m][r] + cost(l, m, r)
    template <class Leaf, class Cost>
    void split(Leaf leaf, Cost cost, bool is_min = true) {
        arg_.assign((size_t)(n_ + 1) * (n_ + 1), -1);
        for (int i = 0; i < n_; ++i) at(i, i + 1) = leaf(i);
        for (int l = n_ - 2; l >= 0; --l)
            for (int r = l + 2; r <= n_; ++r) {
                T best{};
                int bm = -1;
                for (int m = l + 1; m < r; ++m) {
                    const T c = at(l, m) + at(m, r) + cost(l, m, r);
                    if (bm < 0 || (is_min ? c < best : best < c)) best = c, bm = m;
                }
                at(l, r) = best, arg_[idx(l, r)] = bm;
            }
    }
    // Knuth 最適化（最小化。w は (l, r) だけで決まり、四角不等式と単調性を満たす）
    template <class Leaf, class W>
    void split_knuth(Leaf leaf, W w) {
        arg_.assign((size_t)(n_ + 1) * (n_ + 1), -1);
        for (int i = 0; i < n_; ++i) at(i, i + 1) = leaf(i), arg_[idx(i, i + 1)] = i;
        for (int l = n_ - 2; l >= 0; --l)
            for (int r = l + 2; r <= n_; ++r) {
                const int lo = std::max(l + 1, arg_[idx(l, r - 1)]), hi = std::min(r - 1, arg_[idx(l + 1, r)]);
                T best{};
                int bm = -1;
                for (int m = lo; m <= hi; ++m) {
                    const T c = at(l, m) + at(m, r);
                    if (bm < 0 || c < best) best = c, bm = m;
                }
                if (bm < 0) intervaldp_internal::fail("IntervalDP::split_knuth", "分ける位置が見つからない（内部エラー）");
                at(l, r) = best + w(l, r), arg_[idx(l, r)] = bm;
            }
    }

    T get(int l, int r) const {
        if (l < 0 || r > n_ || l > r) intervaldp_internal::fail("IntervalDP::get", "0 ≤ l ≤ r ≤ n にする");
        if (running_ && !(cur_l_ <= l && r <= cur_r_ && r - l < cur_r_ - cur_l_))
            intervaldp_internal::fail("IntervalDP::get", "run の中では今の区間 [l, r) に含まれる、より短い区間しか読めない");
        return dp_[idx(l, r)];
    }
    T get_closed(int l, int r) const { return get(l, r + 1); }
    T all() const { return get(0, n_); }
    int split_point(int l, int r) const {
        if (arg_.empty()) intervaldp_internal::fail("IntervalDP::split_point", "split / split_knuth の後に使う");
        if (l < 0 || r > n_ || r - l < 2) intervaldp_internal::fail("IntervalDP::split_point", "長さ 2 以上の区間 [l, r) にする");
        return arg_[idx(l, r)];
    }
    // [0, n) を分ける木で、合体 (l, m, r) を子が先になる順に
    std::vector<std::tuple<int, int, int>> merge_order() const {
        std::vector<std::tuple<int, int, int>> out;
        if (n_ < 2) return out;
        std::vector<std::pair<int, int>> st{{0, n_}};
        while (!st.empty()) {
            auto [l, r] = st.back();
            st.pop_back();
            const int m = split_point(l, r);
            out.push_back({l, m, r});
            if (m - l >= 2) st.push_back({l, m});
            if (r - m >= 2) st.push_back({m, r});
        }
        std::reverse(out.begin(), out.end());
        return out;
    }

   private:
    int n_;
    std::vector<T> dp_;
    std::vector<int> arg_;
    int cur_l_ = 0, cur_r_ = 0;
    bool running_ = false;
    size_t idx(int l, int r) const { return (size_t)l * (n_ + 1) + r; }
    T& at(int l, int r) { return dp_[idx(l, r)]; }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace intervaldp_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define ID_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[dp_interval] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 12);
        std::vector<ll> a(n), d(n + 1);
        for (auto& x : a) x = rnd(rng, 0, 30);
        for (auto& x : d) x = rnd(rng, 1, 20);
        std::vector<ll> pre(n + 1, 0);
        for (int i = 0; i < n; ++i) pre[i + 1] = pre[i] + a[i];
        auto sum = [&](int l, int r) { return pre[r] - pre[l]; };
        // 愚直（メモ化なしの再帰。n ≤ 12 なので間に合う程度に小さく）
        std::vector<ll> memo((n + 1) * (n + 1), -1);
        auto brute = [&](auto&& self, int l, int r, int kind) -> ll {
            if (r - l == 1) return kind == 2 ? a[l] : 0;
            ll& mm = memo[l * (n + 1) + r];
            if (mm >= 0) return mm;
            ll best = -1;
            for (int m = l + 1; m < r; ++m) {
                ll c = self(self, l, m, kind) + self(self, m, r, kind);
                if (kind == 0) c += sum(l, r);                 // Slimes
                else if (kind == 1) c += d[l] * d[m] * d[r];   // 行列の積
                else c += sum(l, r) * sum(l, r) % 7 * (m - l); // 最大化（適当な式）
                if (best < 0 || (kind == 2 ? c > best : c < best)) best = c;
            }
            return mm = best;
        };
        for (int kind = 0; kind < 3; ++kind) {
            std::fill(memo.begin(), memo.end(), -1);
            IntervalDP<ll> dp(n);
            if (kind == 0) dp.split([](int) { return 0LL; }, [&](int l, int, int r) { return sum(l, r); });
            else if (kind == 1) dp.split([](int) { return 0LL; }, [&](int l, int m, int r) { return d[l] * d[m] * d[r]; });
            else dp.split([&](int i) { return a[i]; }, [&](int l, int m, int r) { return sum(l, r) * sum(l, r) % 7 * (m - l); }, false);
            ID_CHECK(dp.size() == n);
            for (int l = 0; l < n; ++l)
                for (int r = l + 1; r <= n; ++r) {
                    ID_CHECK(dp.get(l, r) == brute(brute, l, r, kind));
                    ID_CHECK(dp.get_closed(l, r - 1) == dp.get(l, r));
                    if (r - l >= 2) {  // 分けた位置で値が合う
                        const int m = dp.split_point(l, r);
                        ID_CHECK(l < m && m < r);
                        ll c = dp.get(l, m) + dp.get(m, r);
                        if (kind == 0) c += sum(l, r);
                        else if (kind == 1) c += d[l] * d[m] * d[r];
                        else c += sum(l, r) * sum(l, r) % 7 * (m - l);
                        ID_CHECK(c == dp.get(l, r));
                    }
                    ++cases;
                }
            if (n >= 1) ID_CHECK(dp.all() == dp.get(0, n));
            // merge_order：n − 1 回の合体で、子が先
            const auto ord = dp.merge_order();
            ID_CHECK((int)ord.size() == std::max(0, n - 1));
            std::vector<std::vector<char>> done(n + 1, std::vector<char>(n + 1, 0));
            for (int i = 0; i < n; ++i) done[i][i + 1] = 1;
            for (auto [l, m, r] : ord) {
                ID_CHECK(done[l][m] && done[m][r] && !done[l][r]);
                done[l][r] = 1;
            }
            if (n >= 1) ID_CHECK(done[0][n]);
        }
        // Knuth：w = 区間の和、区間の和の 2 乗
        for (int kind = 0; kind < 2; ++kind) {
            auto w = [&](int l, int r) { return kind == 0 ? sum(l, r) : sum(l, r) * sum(l, r); };
            IntervalDP<ll> p(n), q(n);
            p.split([](int) { return 0LL; }, [&](int l, int, int r) { return w(l, r); });
            q.split_knuth([](int) { return 0LL; }, w);
            for (int l = 0; l < n; ++l)
                for (int r = l + 1; r <= n; ++r) {
                    ID_CHECK(p.get(l, r) == q.get(l, r));
                    if (r - l >= 2) {
                        const int m = q.split_point(l, r);
                        ID_CHECK(q.get(l, m) + q.get(m, r) + w(l, r) == q.get(l, r));
                    }
                    ++cases;
                }
        }
        // run：回文にする最小の挿入回数（文字は 0〜2）と、愚直の再帰
        {
            std::vector<int> s(n);
            for (auto& c : s) c = (int)rnd(rng, 0, 2);
            IntervalDP<int> p(n);
            p.run([](int) { return 0; }, [&](int l, int r, const IntervalDP<int>& dp) {
                if (r - l == 1) return 0;
                if (s[l] == s[r - 1]) return dp.get(l + 1, r - 1);
                return 1 + std::min(dp.get(l + 1, r), dp.get(l, r - 1));
            });
            auto br = [&](auto&& self, int l, int r) -> int {
                if (r - l <= 1) return 0;
                if (s[l] == s[r - 1]) return self(self, l + 1, r - 1);
                return 1 + std::min(self(self, l + 1, r), self(self, l, r - 1));
            };
            for (int l = 0; l <= n; ++l)
                for (int r = l; r <= n && r - l <= 9; ++r) {
                    ID_CHECK(p.get(l, r) == br(br, l, r));
                    ++cases;
                }
        }
    }
    std::cout << "[dp_interval] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[dp_interval] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n1 = TEST_SCALE == 1 ? 120 : 400, n2 = TEST_SCALE == 1 ? 800 : 5000;
    std::vector<ll> a(n2), pre(n2 + 1, 0);
    for (auto& x : a) x = rnd(rng, 1, 1000000000);
    for (int i = 0; i < n2; ++i) pre[i + 1] = pre[i] + a[i];
    auto t0 = std::chrono::steady_clock::now();
    ll chk = 0;
    {
        IntervalDP<ll> dp(n1);
        dp.split([](int) { return 0LL; }, [&](int l, int, int r) { return pre[r] - pre[l]; });
        chk += dp.all();
    }
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    {
        IntervalDP<ll> dp(n2);
        dp.split_knuth([](int) { return 0LL; }, [&](int l, int r) { return pre[r] - pre[l]; });
        chk += dp.all();
    }
    const ll t2 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    {
        IntervalDP<ll> dp(n2);
        dp.run([](int) { return 0LL; }, [&](int l, int r, const IntervalDP<ll>& g) {
            if (r - l == 1) return a[l];
            return std::max(a[l] - g.get(l + 1, r), a[r - 1] - g.get(l, r - 1));  // EDPC L（Deque）の形
        });
        chk += dp.all();
    }
    const ll t3 = elapsed_ms(t0);
    std::cout << "[dp_interval] benchmark: split(n = " << n1 << ") " << t1 << " ms / split_knuth(n = " << n2 << ") " << t2 << " ms / run(n = "
              << n2 << ") " << t3 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace intervaldp_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace intervaldp_demo {

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
#define SHOW(expr, comment) intervaldp_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    intervaldp_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("区間 DP（dp_interval） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  IntervalDP<T> dp(n)   run(base, trans) / split(leaf, cost, is_min) / split_knuth(leaf, w)\n"
              << "  get(l, r) / get_closed(l, r) / all() / split_point(l, r) / merge_order() / size()\n";

    api("【1】 split(leaf, cost)：隣どうしを合体させる最小コスト（EDPC N Slimes）", "void",
        "dp[i][i+1] = leaf(i)、dp[l][r] = min_m dp[l][m] + dp[m][r] + cost(l, m, r)。Slimes は cost = 区間の和", "O(n^3)");
    RUN(const std::vector<long long> a = {10, 20, 30, 40});
    RUN(std::vector<long long> pre(5, 0));
    RUN(for (int i = 0; i < 4; ++i) pre[i + 1] = pre[i] + a[i]);
    RUN(IntervalDP<long long> dp(4));
    RUN(dp.split([](int) { return 0LL; }, [&](int l, int, int r) { return pre[r] - pre[l]; }));
    SHOW(dp.all(), "(10 + 20) → (30 + 30) → (60 + 40) の 30 + 60 + 100");
    SHOW(dp.get(1, 4), "{20, 30, 40} だけなら");
    SHOW(dp.get_closed(1, 3), "閉区間 [1, 3] = 同じ");

    api("【2】 split_point(l, r) / merge_order()", "int / vector<tuple<int, int, int>>", "[l, r) を分けた位置 / 合体の順番 (l, m, r)（子が先）",
        "O(1) / O(n)");
    SHOW(dp.split_point(0, 4), "[0, 3) と [3, 4) に分ける");
    SHOW(dp.merge_order().size(), "合体は n − 1 回");
    {
        std::string s;
        for (auto [l, m, r] : dp.merge_order()) s += "([" + to_s(l) + "," + to_s(m) + ")+[" + to_s(m) + "," + to_s(r) + ")) ";
        note("merge_order() = " + s);
    }

    api("【3】 split_knuth(leaf, w)", "void", "split で cost が w(l, r) だけで決まり、四角不等式・単調性を満たすとき（区間の和など）。結果は split と同じ",
        "O(n^2)");
    RUN(IntervalDP<long long> kn(4));
    RUN(kn.split_knuth([](int) { return 0LL; }, [&](int l, int r) { return pre[r] - pre[l]; }));
    SHOW(kn.all(), "split と同じ 190");

    api("【4】 split(leaf, cost) で cost が m による例：行列の積の順番", "void", "行列 i が d[i] × d[i+1]。掛け算の回数 d[l] d[m] d[r] の和の最小", "O(n^3)");
    RUN(const std::vector<long long> dim = {10, 30, 5, 60});
    RUN(IntervalDP<long long> mc(3));
    RUN(mc.split([](int) { return 0LL; }, [&](int l, int m, int r) { return dim[l] * dim[m] * dim[r]; }));
    SHOW(mc.all(), "(A B) C = 10·30·5 + 10·5·60");
    SHOW(mc.split_point(0, 3), "");

    api("【5】 split(…, false)：最大化", "void", "is_min = false にすると max", "O(n^3)");
    RUN(IntervalDP<long long> mx(4));
    RUN(mx.split([](int) { return 0LL; }, [&](int l, int, int r) { return pre[r] - pre[l]; }, false));
    SHOW(mx.all(), "合体コストの最大（(10 + (20 + (30 + 40))) など）");

    api("【6】 run(base, trans)：一般形", "void",
        "dp[i][i] = base(i)、dp[l][r] = trans(l, r, dp)。trans の中で dp.get(a, b) は [l, r) に含まれる短い区間だけ", "O(n^2 × trans)");
    code("// 回文にするために足す文字の最小の数");
    RUN(const std::string s = "abcab");
    RUN(IntervalDP<int> pal((int)s.size()));
    RUN(pal.run([](int) { return 0; }, [&](int l, int r, const IntervalDP<int>& g) { if (r - l == 1) return 0; if (s[l] == s[r - 1]) return g.get(l + 1, r - 1); return 1 + std::min(g.get(l + 1, r), g.get(l, r - 1)); }));
    SHOW(pal.all(), "\"abcab\" → \"bacbcab\" など 2 文字");
    SHOW(pal.get(0, 3), "\"abc\" は 2 文字");
    note("EDPC L（Deque）: dp[l][r] = max(a[l] − dp[l+1][r], a[r−1] − dp[l][r−1]) も run で書ける");

    api("【7】 円環の区間 DP", "—", "列を 2 周分つないだ長さ 2n で作り、長さ n の区間 get(i, i + n) の最小を取る", "—");
    {
        code("// 円く並んだスライム {40, 10, 20, 30}（40 と 30 も隣）");
        const std::vector<long long> c = {40, 10, 20, 30};
        const int n = 4;
        std::vector<long long> p2(2 * n + 1, 0);
        for (int i = 0; i < 2 * n; ++i) p2[i + 1] = p2[i] + c[i % n];
        IntervalDP<long long> cd(2 * n);
        cd.split([](int) { return 0LL; }, [&](int l, int, int r) { return p2[r] - p2[l]; });
        long long best = cd.get(0, n);
        for (int i = 1; i < n; ++i) best = std::min(best, cd.get(i, i + n));
        note("→ 最小コスト " + to_s(best) + "（並べ方は【1】と同じなので 190）");
    }
}

}  // namespace intervaldp_demo

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
    intervaldp_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    intervaldp_test::benchmark();    // 速度計測
    intervaldp_demo::demo();         // 使い方デモ
    return 0;
}
