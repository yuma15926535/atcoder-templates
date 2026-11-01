// =============================================================================
//  SWAG（Sliding Window Aggregation：キュー・両端キューの中身全体の集約を、逆元なしで O(1)（ならし））   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  キューに push（後ろに入れる）・pop（前から出す）しながら、中身全部を op でまとめた値を答える。
//  op は結合的（(a op b) op c = a op (b op c)）なら何でもよく、交換できなくても・逆元がなくてもよい
//  （min・gcd・行列の積・一次関数の合成など）。スタック 2 本で持ち、片方が空になったら移す
//    例: 長さ K の窓の gcd を全部、窓の中の一次関数をこの順に合成した値、窓の行列の積
//
//  機能一覧（S は値の型、op(a, b) は「a の後ろに b」をまとめる、e は空のときの値）              計算量
//    SWAG<S, Op> q(op, e)               キュー版                                                     O(1)
//      q.push(x) / q.pop()              後ろに入れる / 前から出す                                    ならし O(1)
//      q.fold()                         中身を前から順に op でまとめた値（空なら e）                 O(1)
//      q.front() / q.size() / q.empty() / q.clear()
//    DequeSWAG<S, Op> d(op, e)          両端キュー版
//      d.push_front(x) / d.push_back(x) / d.pop_front() / d.pop_back()                            ならし O(1)
//      d.fold() / d.front() / d.back() / d.size() / d.empty() / d.clear()
//    make_swag(op, e) / make_deque_swag(op, e)   ラムダを渡して作る（型を書かなくてよい）
//
//  注意: 空のときに pop・front すると止める。min / max だけなら SlidingWindowDequeTemplateTest.txt（単調な deque）でもよい
//  速さの目安: push・pop・fold を 10^7 回ずつ（op = 一次関数の合成 mod）約 0.3 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    swag_test::random_test()  … 中身を毎回まとめ直す愚直解との突き合わせ（交換できない op で）
//    swag_test::benchmark()    … 速度計測
//    swag_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <deque>
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
namespace swag_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace swag_internal

template <class S, class Op>
struct SWAG {
    SWAG(Op op, S e) : op_(op), e_(e), back_fold_(e) {}
    void push(const S& x) {
        back_.push_back(x);
        back_fold_ = op_(back_fold_, x);
    }
    void pop() {
        if (empty()) swag_internal::fail("SWAG::pop", "空");
        if (front_.empty()) move_back();
        front_.pop_back();
    }
    S fold() const { return front_.empty() ? back_fold_ : op_(front_.back().second, back_fold_); }
    const S& front() {
        if (empty()) swag_internal::fail("SWAG::front", "空");
        if (front_.empty()) move_back();
        return front_.back().first;
    }
    size_t size() const { return front_.size() + back_.size(); }
    bool empty() const { return size() == 0; }
    void clear() { front_.clear(), back_.clear(), back_fold_ = e_; }

   private:
    Op op_;
    S e_;
    std::vector<std::pair<S, S>> front_;  // (値, その値から前の山の底まで（キューの順）をまとめた値)。末尾がキューの先頭
    std::vector<S> back_;
    S back_fold_;
    void move_back() {
        for (size_t i = back_.size(); i-- > 0;) {
            const S f = front_.empty() ? back_[i] : op_(back_[i], front_.back().second);
            front_.push_back({back_[i], f});
        }
        back_.clear(), back_fold_ = e_;
    }
};

template <class S, class Op>
struct DequeSWAG {
    DequeSWAG(Op op, S e) : op_(op), e_(e) {}
    void push_front(const S& x) { F_.push_back({x, F_.empty() ? x : op_(x, F_.back().second)}); }
    void push_back(const S& x) { B_.push_back({x, B_.empty() ? x : op_(B_.back().second, x)}); }
    void pop_front() {
        if (empty()) swag_internal::fail("DequeSWAG::pop_front", "空");
        if (F_.empty()) rebalance(true);
        F_.pop_back();
    }
    void pop_back() {
        if (empty()) swag_internal::fail("DequeSWAG::pop_back", "空");
        if (B_.empty()) rebalance(false);
        B_.pop_back();
    }
    S fold() const {
        if (F_.empty() && B_.empty()) return e_;
        if (F_.empty()) return B_.back().second;
        if (B_.empty()) return F_.back().second;
        return op_(F_.back().second, B_.back().second);
    }
    const S& front() {
        if (empty()) swag_internal::fail("DequeSWAG::front", "空");
        if (F_.empty()) rebalance(true);
        return F_.back().first;
    }
    const S& back() {
        if (empty()) swag_internal::fail("DequeSWAG::back", "空");
        if (B_.empty()) rebalance(false);
        return B_.back().first;
    }
    size_t size() const { return F_.size() + B_.size(); }
    bool empty() const { return size() == 0; }
    void clear() { F_.clear(), B_.clear(); }

   private:
    Op op_;
    S e_;
    std::vector<std::pair<S, S>> F_, B_;  // F_ の末尾 = 先頭の要素、B_ の末尾 = 最後の要素
    // 空になった側へ、反対側の半分を移す（中身の順番：F_ の末尾 → F_ の先頭 → B_ の先頭 → B_ の末尾）
    void rebalance(bool to_front) {
        std::vector<S> all;  // キューの順
        for (size_t i = F_.size(); i-- > 0;) all.push_back(F_[i].first);
        for (const auto& p : B_) all.push_back(p.first);
        F_.clear(), B_.clear();
        const size_t n = all.size();
        size_t mid = n / 2;  // [0, mid) を前、[mid, n) を後ろへ
        if (to_front && mid == 0) mid = 1;
        if (!to_front && mid == n) mid = n - 1;
        for (size_t i = mid; i-- > 0;) push_front(all[i]);
        for (size_t i = mid; i < n; ++i) push_back(all[i]);
    }
};
template <class S, class Op>
SWAG<S, Op> make_swag(Op op, S e) {
    return SWAG<S, Op>(op, e);
}
template <class S, class Op>
DequeSWAG<S, Op> make_deque_swag(Op op, S e) {
    return DequeSWAG<S, Op>(op, e);
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace swag_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define SW_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_swag] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
using F = std::pair<ll, ll>;  // 一次関数 x → a x + b（mod）。op(f, g) = 先に f、次に g
static const ll MOD = 998244353;
static F comp(const F& f, const F& g) { return {f.first * g.first % MOD, (f.second * g.first + g.second) % MOD}; }

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        auto q = make_swag(comp, F{1, 0});
        auto d = make_deque_swag(comp, F{1, 0});
        std::deque<F> bq, bd;
        auto fold = [](const std::deque<F>& v) {
            F r{1, 0};
            for (const F& f : v) r = comp(r, f);
            return r;
        };
        const int ops = (int)rnd(rng, 1, 400);
        for (int t = 0; t < ops; ++t) {
            const F f{rnd(rng, 0, MOD - 1), rnd(rng, 0, MOD - 1)};
            const int k = (int)rnd(rng, 0, 9);
            if (k < 5 || bq.empty()) q.push(f), bq.push_back(f);
            else {
                SW_CHECK(q.front() == bq.front());
                q.pop(), bq.pop_front();
            }
            const int k2 = (int)rnd(rng, 0, 9);
            if (bd.empty() || k2 < 3) d.push_back(f), bd.push_back(f);
            else if (k2 < 5) d.push_front(f), bd.push_front(f);
            else if (k2 < 7) {
                SW_CHECK(d.front() == bd.front());
                d.pop_front(), bd.pop_front();
            } else {
                SW_CHECK(d.back() == bd.back());
                d.pop_back(), bd.pop_back();
            }
            SW_CHECK(q.size() == bq.size() && q.empty() == bq.empty() && q.fold() == fold(bq));
            SW_CHECK(d.size() == bd.size() && d.empty() == bd.empty() && d.fold() == fold(bd));
            cases += 2;
        }
        q.clear(), d.clear();
        SW_CHECK(q.empty() && d.empty() && q.fold() == (F{1, 0}) && d.fold() == (F{1, 0}));
    }
    // 片側だけ使い続ける（偏った移し替え）
    {
        auto d = make_deque_swag([](ll a, ll b) { return a * 10 + b; }, 0LL);  // 交換できない：数字を並べる
        for (int i = 1; i <= 5; ++i) d.push_back(i);
        d.pop_front(), d.pop_front();
        SW_CHECK(d.fold() == 345);
        for (int i = 0; i < 3; ++i) d.pop_back();
        SW_CHECK(d.empty() && d.fold() == 0);
        d.push_front(7), d.push_front(8);
        SW_CHECK(d.fold() == 87 && d.back() == 7 && d.front() == 8);
        cases += 3;
    }
    std::cout << "[ds_swag] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_swag] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 1000000 : 10000000, K = 1000;
    std::vector<F> fs(1 << 16);
    for (auto& f : fs) f = {rnd(rng, 0, MOD - 1), rnd(rng, 0, MOD - 1)};
    const auto t0 = std::chrono::steady_clock::now();
    auto q = make_swag(comp, F{1, 0});
    ll chk = 0;
    for (int i = 0; i < N; ++i) {
        q.push(fs[i & 65535]);
        if ((int)q.size() > K) q.pop();
        chk += q.fold().second;
    }
    std::cout << "[ds_swag] benchmark: 窓 " << K << " で push・pop・fold を " << N << " 回（一次関数の合成） " << elapsed_ms(t0) << " ms  (check " << chk % 1000
              << ")\n";
}

}  // namespace swag_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace swag_demo {

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
#define SHOW(expr, comment) swag_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    swag_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("SWAG（ds_swag） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  SWAG<S, Op> / make_swag(op, e): push / pop / fold / front / size / empty / clear\n"
              << "  DequeSWAG<S, Op> / make_deque_swag(op, e): push_front / push_back / pop_front / pop_back / fold / front / back / …\n";

    api("【1】 make_swag(op, e) / push(x) / pop() / fold()", "— / void / void / S", "キューの中身を前から op でまとめた値（逆元がない op でよい）", "ならし O(1)");
    RUN(auto q = make_swag([](long long a, long long b) { return std::gcd(a, b); }, 0LL));
    RUN(q.push(12); q.push(18); q.push(30));
    SHOW(q.fold(), "gcd(12, 18, 30)");
    RUN(q.pop());
    SHOW(q.fold(), "gcd(18, 30)");
    SHOW(q.front(), "");
    {
        code("// 長さ 3 の窓の gcd を全部");
        const std::vector<long long> a = {6, 9, 15, 10, 20, 8};
        auto w = make_swag([](long long x, long long y) { return std::gcd(x, y); }, 0LL);
        std::vector<long long> res;
        for (size_t i = 0; i < a.size(); ++i) {
            w.push(a[i]);
            if (w.size() > 3) w.pop();
            if (w.size() == 3) res.push_back(w.fold());
        }
        note("a = " + to_s(a) + " → " + to_s(res));
    }

    api("【2】 交換できない op：一次関数の合成", "—", "op(f, g) = 「先に f、次に g」。窓の中の関数をこの順に当てた結果", "—");
    RUN(using Fn = std::pair<long long, long long>);  // x → a x + b
    RUN(auto fq = make_swag([](Fn f, Fn g) { return Fn{f.first * g.first, f.second * g.first + g.second}; }, Fn{1, 0}));
    RUN(fq.push({2, 1}); fq.push({3, 0}));  // x → 2x + 1 → 3(2x + 1)
    SHOW(fq.fold(), "6x + 3");

    api("【3】 make_deque_swag(op, e)", "—", "両端から入れたり出したりできる版", "ならし O(1)");
    RUN(auto d = make_deque_swag([](std::string a, std::string b) { return a + b; }, std::string()));
    RUN(d.push_back("b"); d.push_back("c"); d.push_front("a"));
    SHOW(d.fold(), "");
    RUN(d.pop_back());
    SHOW(d.fold(), "");
}

}  // namespace swag_demo

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
    swag_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    swag_test::benchmark();    // 速度計測
    swag_demo::demo();         // 使い方デモ
    return 0;
}
