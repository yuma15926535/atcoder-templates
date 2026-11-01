// =============================================================================
//  Disjoint Sparse Table（更新のない配列で、任意のモノイドの区間積を O(1)）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  普通の Sparse Table は区間が重なってもよい演算（min・max・gcd など、同じものを 2 回まとめても変わらない）専用。
//  Disjoint Sparse Table は「真ん中から左右へ伸ばした積」を段ごとに持つので、和・積・mod の積・行列・一次関数の合成など
//  結合的なら何でも O(1) で答えられる（前処理 O(N log N)）
//    例: 「区間の積 mod M を Q 回」（逆元がなくてよい）、「区間の一次関数をこの順に合成した値」
//
//  機能一覧（区間は 0-indexed・半開区間 [l, r)。閉区間は _closed 版）                            計算量
//    make_disjoint_sparse_table(a, op, e)   配列 a から作る（op は結合的、e は空の区間の値）           O(N log N)
//    prod(l, r) / prod_closed(l, r)         op(a[l], …, a[r − 1]) / [l, r]（l = r なら e）             O(1)
//    size()                                 N
//
//  注意: 更新はできない（更新があるならセグメント木）。メモリは N log N 個の S
//  速さの目安: N = 2×10^5 で構築 + prod を 2×10^5 回 約 0.03 秒（op = 和）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    dst_test::random_test()  … 区間を 1 つずつまとめる愚直解との突き合わせ（交換できない op で）
//    dst_test::benchmark()    … 速度計測
//    dst_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace dst_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace dst_internal

template <class S, class Op>
struct DisjointSparseTable {
    DisjointSparseTable(const std::vector<S>& a, Op op, S e) : n_((int)a.size()), op_(op), e_(e) {
        lg_ = 1;
        while ((1 << lg_) < n_) ++lg_;
        const int sz = 1 << lg_;
        t_.assign(lg_ + 1, std::vector<S>(sz, e_));
        t_[0].assign(a.begin(), a.end());
        t_[0].resize(sz, e_);
        for (int h = 1; h <= lg_; ++h) {
            const int half = 1 << (h - 1);
            auto& row = t_[h];
            for (int mid = half; mid < sz; mid += 2 * half) {
                // 左側：row[i] = a[i] … a[mid − 1]、右側：row[i] = a[mid] … a[i]
                row[mid - 1] = t_[0][mid - 1];
                for (int i = mid - 2; i >= mid - half; --i) row[i] = op_(t_[0][i], row[i + 1]);
                row[mid] = t_[0][mid];
                for (int i = mid + 1; i < mid + half; ++i) row[i] = op_(row[i - 1], t_[0][i]);
            }
        }
    }
    int size() const { return n_; }
    S prod(int l, int r) const {
        if (l < 0 || r > n_ || l > r) dst_internal::fail("DisjointSparseTable::prod", "0 ≤ l ≤ r ≤ N の半開区間 [l, r)");
        if (l == r) return e_;
        --r;
        if (l == r) return t_[0][l];
        const int h = 32 - __builtin_clz((unsigned)(l ^ r));  // l と r が分かれる段
        return op_(t_[h][l], t_[h][r]);
    }
    S prod_closed(int l, int r) const { return prod(l, r + 1); }

   private:
    int n_, lg_;
    Op op_;
    S e_;
    std::vector<std::vector<S>> t_;
};
template <class S, class Op>
DisjointSparseTable<S, Op> make_disjoint_sparse_table(const std::vector<S>& a, Op op, S e) {
    return DisjointSparseTable<S, Op>(a, op, e);
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace dst_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define DS_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_disjoint_sparse_table] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
using F = std::pair<ll, ll>;
static const ll MOD = 998244353;
static F comp(const F& f, const F& g) { return {f.first * g.first % MOD, (f.second * g.first + g.second) % MOD}; }

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, rnd(rng, 0, 3) ? 20 : 130);
        std::vector<F> a(n);
        for (auto& f : a) f = {rnd(rng, 0, MOD - 1), rnd(rng, 0, MOD - 1)};
        const auto dst = make_disjoint_sparse_table(a, comp, F{1, 0});
        DS_CHECK(dst.size() == n);
        for (int q = 0; q < 60; ++q) {
            int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
            if (l > r) std::swap(l, r);
            F x{1, 0};
            for (int i = l; i < r; ++i) x = comp(x, a[i]);
            DS_CHECK(dst.prod(l, r) == x);
            if (l < r) DS_CHECK(dst.prod_closed(l, r - 1) == x);
            ++cases;
        }
    }
    std::cout << "[ds_disjoint_sparse_table] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_disjoint_sparse_table] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    std::vector<ll> a(N);
    for (auto& x : a) x = rnd(rng, 0, 1000000000);
    const auto t0 = std::chrono::steady_clock::now();
    const auto dst = make_disjoint_sparse_table(a, [](ll x, ll y) { return x + y; }, 0LL);
    ll chk = 0;
    for (int i = 0; i < N; ++i) {
        int l = (int)rnd(rng, 0, N), r = (int)rnd(rng, 0, N);
        if (l > r) std::swap(l, r);
        chk += dst.prod(l, r);
    }
    std::cout << "[ds_disjoint_sparse_table] benchmark (N = " << N << "): 構築 + prod を " << N << " 回 " << elapsed_ms(t0) << " ms  (check " << chk % 1000
              << ")\n";
}

}  // namespace dst_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace dst_demo {

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
#define SHOW(expr, comment) dst_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    dst_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("Disjoint Sparse Table（ds_disjoint_sparse_table） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  make_disjoint_sparse_table(a, op, e): prod(l, r) / prod_closed(l, r) / size()\n";

    api("【1】 make_disjoint_sparse_table(a, op, e) / prod(l, r)", "— / S", "更新のない配列で、区間の op を O(1)（重なってはいけない演算でもよい）",
        "O(N log N) / O(1)");
    RUN(const long long M = 1000000007);
    RUN(const auto mul = make_disjoint_sparse_table(std::vector<long long>{2, 3, 5, 7, 11}, [&](long long x, long long y) { return x * y % M; }, 1LL));
    SHOW(mul.prod(0, 5), "2 · 3 · 5 · 7 · 11（逆元がなくてもよい）");
    SHOW(mul.prod(1, 3), "3 · 5");
    SHOW(mul.prod_closed(4, 4), "");
    SHOW(mul.prod(2, 2), "空の区間は e");

    api("【2】 交換できない op：文字列をつなぐ", "—", "op(a, b) = a の後ろに b", "O(1)");
    RUN(const auto cat = make_disjoint_sparse_table(std::vector<std::string>{"ab", "c", "de", "f"}, [](const std::string& x, const std::string& y) { return x + y; }, std::string()));
    SHOW(cat.prod(1, 4), "");
    note("min・max・gcd だけなら SparseTableTemplateTest.txt（普通の Sparse Table）でもよい");
}

}  // namespace dst_demo

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
    dst_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    dst_test::benchmark();    // 速度計測
    dst_demo::demo();         // 使い方デモ
    return 0;
}
