// =============================================================================
//  永続セグメント木（一般のモノイド。1 点更新で新しい版を作り、どの版でも区間の集約・二分探索ができる）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  set(v, i, x) は版 v の根から葉までの O(log N) 個の頂点だけをコピーして新しい版を作る。古い版も残る。
//  S・op・e は ACL の segtree と同じ（op は結合的、e は単位元）
//    例: 「時刻 t の時点での区間の和」を時刻ごとに答える（時刻 = 版）、版 r = 先頭 r 個を入れた木 として区間の数を数える
//
//  機能一覧（区間は 0-indexed・半開区間 [l, r)。閉区間は _closed 版。版は int の番号で 0 が最初）   計算量
//    make_persistent_segtree(n, op, e)       全部 e の長さ n（版 0）                                   O(n)
//    make_persistent_segtree(a, op, e)       配列 a（版 0）
//    set(v, i, x) / update(v, i, f)          版 v の a[i] を x / f(a[i]) にした新しい版の番号           O(log n)
//    get(v, i)                               版 v の a[i]                                              O(log n)
//    prod(v, l, r) / prod_closed(v, l, r)    版 v の op(a[l], …, a[r − 1]) / [l, r]                    O(log n)
//    all_prod(v)                             版 v の全体                                               O(1)
//    max_right(v, l, pred)                   pred(prod(v, l, r)) が true の最大の r（pred(e) は true）   O(log n)
//    versions() / size()                     版の数 / n
//
//  注意: 頂点は更新 1 回で約 log2(n) + 1 個増える（n = 2×10^5、更新 2×10^5 回で約 400 万個）
//  速さの目安: n = 2×10^5 で set と prod を 2×10^5 回ずつ（版はばらばら）約 0.4 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    pseg_test::random_test()  … 版ごとに配列を持つ愚直解との突き合わせ（交換できない op で）
//    pseg_test::benchmark()    … 速度計測
//    pseg_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
namespace pseg_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace pseg_internal

template <class S, class Op>
struct PersistentSegTree {
    PersistentSegTree(int n, Op op, S e) : PersistentSegTree(std::vector<S>(n, e), op, e) {}
    PersistentSegTree(const std::vector<S>& a, Op op, S e) : n_((int)a.size()), op_(op), e_(e) {
        nd_.reserve((size_t)2 * std::max(1, n_));
        roots_.push_back(n_ == 0 ? new_node(-1, -1, e_) : build(a, 0, n_));
    }
    int size() const { return n_; }
    int versions() const { return (int)roots_.size(); }
    int set(int v, int i, const S& x) {
        check_v(v, "PersistentSegTree::set"), check_i(i, "PersistentSegTree::set");
        roots_.push_back(set_rec(roots_[v], 0, n_, i, x));
        return (int)roots_.size() - 1;
    }
    template <class F>
    int update(int v, int i, F f) {
        return set(v, i, f(get(v, i)));
    }
    S get(int v, int i) const {
        check_v(v, "PersistentSegTree::get"), check_i(i, "PersistentSegTree::get");
        int cur = roots_[v], lo = 0, hi = n_;
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            if (i < mid) cur = nd_[cur].l, hi = mid;
            else cur = nd_[cur].r, lo = mid;
        }
        return nd_[cur].val;
    }
    S prod(int v, int l, int r) const {
        check_v(v, "PersistentSegTree::prod");
        if (l < 0 || r > n_ || l > r) pseg_internal::fail("PersistentSegTree::prod", "0 ≤ l ≤ r ≤ n の半開区間 [l, r)");
        if (l == r) return e_;
        return prod_rec(roots_[v], 0, n_, l, r);
    }
    S prod_closed(int v, int l, int r) const { return prod(v, l, r + 1); }
    S all_prod(int v) const {
        check_v(v, "PersistentSegTree::all_prod");
        return n_ == 0 ? e_ : nd_[roots_[v]].val;
    }
    template <class P>
    int max_right(int v, int l, P pred) const {
        check_v(v, "PersistentSegTree::max_right");
        if (l < 0 || l > n_) pseg_internal::fail("PersistentSegTree::max_right", "0 ≤ l ≤ n");
        if (!pred(e_)) pseg_internal::fail("PersistentSegTree::max_right", "pred(e) は true にする");
        if (l == n_) return n_;
        S acc = e_;
        const int r = mr_rec(roots_[v], 0, n_, l, pred, acc);
        return r < 0 ? n_ : r;
    }

   private:
    struct Node {
        int l, r;
        S val;
    };
    int n_;
    Op op_;
    S e_;
    std::vector<Node> nd_;
    std::vector<int> roots_;
    int new_node(int l, int r, const S& x) {
        nd_.push_back({l, r, x});
        return (int)nd_.size() - 1;
    }
    int build(const std::vector<S>& a, int lo, int hi) {
        if (hi - lo == 1) return new_node(-1, -1, a[lo]);
        const int mid = (lo + hi) / 2;
        const int l = build(a, lo, mid), r = build(a, mid, hi);
        return new_node(l, r, op_(nd_[l].val, nd_[r].val));
    }
    int set_rec(int node, int lo, int hi, int i, const S& x) {
        if (hi - lo == 1) return new_node(-1, -1, x);
        const int mid = (lo + hi) / 2, l = nd_[node].l, r = nd_[node].r;
        if (i < mid) {
            const int nl = set_rec(l, lo, mid, i, x);
            return new_node(nl, r, op_(nd_[nl].val, nd_[r].val));
        }
        const int nr = set_rec(r, mid, hi, i, x);
        return new_node(l, nr, op_(nd_[l].val, nd_[nr].val));
    }
    S prod_rec(int node, int lo, int hi, int l, int r) const {
        if (l <= lo && hi <= r) return nd_[node].val;
        const int mid = (lo + hi) / 2;
        if (r <= mid) return prod_rec(nd_[node].l, lo, mid, l, r);
        if (mid <= l) return prod_rec(nd_[node].r, mid, hi, l, r);
        return op_(prod_rec(nd_[node].l, lo, mid, l, r), prod_rec(nd_[node].r, mid, hi, l, r));
    }
    // [l, n) を左から見て、pred が false になる最初の位置（無ければ −1）
    template <class P>
    int mr_rec(int node, int lo, int hi, int l, P& pred, S& acc) const {
        if (hi <= l) return -1;
        if (l <= lo) {
            const S t = op_(acc, nd_[node].val);
            if (pred(t)) {
                acc = t;
                return -1;
            }
            if (hi - lo == 1) return lo;
        }
        const int mid = (lo + hi) / 2;
        const int a = mr_rec(nd_[node].l, lo, mid, l, pred, acc);
        if (a >= 0) return a;
        return mr_rec(nd_[node].r, mid, hi, l, pred, acc);
    }
    void check_v(int v, const char* f) const {
        if (v < 0 || v >= (int)roots_.size()) pseg_internal::fail(f, "版の番号が範囲外");
    }
    void check_i(int i, const char* f) const {
        if (i < 0 || i >= n_) pseg_internal::fail(f, "i が範囲外");
    }
};
template <class S, class Op>
PersistentSegTree<S, Op> make_persistent_segtree(int n, Op op, S e) {
    return PersistentSegTree<S, Op>(n, op, e);
}
template <class S, class Op>
PersistentSegTree<S, Op> make_persistent_segtree(const std::vector<S>& a, Op op, S e) {
    return PersistentSegTree<S, Op>(a, op, e);
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace pseg_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define PS_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_persistent_segtree] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 12);
        std::vector<F> a(n);
        for (auto& f : a) f = {rnd(rng, 0, 5), rnd(rng, 0, 5)};
        const bool use_a = rnd(rng, 0, 1);
        auto st = use_a ? make_persistent_segtree(a, comp, F{1, 0}) : make_persistent_segtree(n, comp, F{1, 0});
        std::vector<std::vector<F>> ver{use_a ? a : std::vector<F>(n, F{1, 0})};
        PS_CHECK(st.size() == n);
        for (int i = 0; i < n; ++i) PS_CHECK(st.get(0, i) == ver[0][i]);
        auto fold = [](const std::vector<F>& v, int l, int r) {
            F x{1, 0};
            for (int i = l; i < r; ++i) x = comp(x, v[i]);
            return x;
        };
        for (int t = 0; t < 100; ++t) {
            const int v = (int)rnd(rng, 0, (int)ver.size() - 1);
            const int k = (int)rnd(rng, 0, 3);
            if (n > 0 && k == 0) {
                const int i = (int)rnd(rng, 0, n - 1);
                const F x{rnd(rng, 0, 5), rnd(rng, 0, 5)};
                PS_CHECK(st.set(v, i, x) == (int)ver.size());
                ver.push_back(ver[v]), ver.back()[i] = x;
            } else if (n > 0 && k == 1) {
                const int i = (int)rnd(rng, 0, n - 1);
                st.update(v, i, [](F f) { return F{(f.first + 1) % MOD, f.second}; });
                ver.push_back(ver[v]), ver.back()[i].first = (ver.back()[i].first + 1) % MOD;
            } else if (k == 2) {
                int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
                if (l > r) std::swap(l, r);
                PS_CHECK(st.prod(v, l, r) == fold(ver[v], l, r));
                if (l < r) PS_CHECK(st.prod_closed(v, l, r - 1) == fold(ver[v], l, r));
                PS_CHECK(st.all_prod(v) == fold(ver[v], 0, n));
            } else {
                // max_right：pred = 「合成した関数の b が limit 以下」（単調ではないので和で）→ 和のモノイドで別に確かめる
                const int l = (int)rnd(rng, 0, n);
                const ll lim = rnd(rng, 0, 40);
                // 和の木を同じ版の値から作り直して比べる
                std::vector<ll> s(n);
                for (int i = 0; i < n; ++i) s[i] = ver[v][i].first + ver[v][i].second;
                auto ss = make_persistent_segtree(s, [](ll x, ll y) { return x + y; }, 0LL);
                int exp = l;
                ll acc = 0;
                while (exp < n && acc + s[exp] <= lim) acc += s[exp], ++exp;
                PS_CHECK(ss.max_right(0, l, [&](ll x) { return x <= lim; }) == exp);
            }
            ++cases;
        }
        PS_CHECK(st.versions() == (int)ver.size());
    }
    std::cout << "[ds_persistent_segtree] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_persistent_segtree] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    const auto t0 = std::chrono::steady_clock::now();
    auto st = make_persistent_segtree(N, [](ll a, ll b) { return a + b; }, 0LL);
    ll chk = 0;
    for (int i = 0; i < N; ++i) {
        st.set(st.versions() - 1, (int)rnd(rng, 0, N - 1), i);
        int l = (int)rnd(rng, 0, N), r = (int)rnd(rng, 0, N);
        if (l > r) std::swap(l, r);
        chk += st.prod((int)rnd(rng, 0, st.versions() - 1), l, r);
    }
    std::cout << "[ds_persistent_segtree] benchmark (n = " << N << "): set と prod を " << N << " 回ずつ " << elapsed_ms(t0) << " ms  (check " << chk % 1000
              << ")\n";
}

}  // namespace pseg_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace pseg_demo {

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
#define SHOW(expr, comment) pseg_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    pseg_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("永続セグメント木（ds_persistent_segtree） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  make_persistent_segtree(n or a, op, e): set / update → 新しい版 / get / prod / prod_closed / all_prod / max_right / versions / size\n";

    api("【1】 make_persistent_segtree(a, op, e) / set(v, i, x) / prod(v, l, r)", "— / int（新しい版）/ S",
        "版 v を元に 1 点変えた新しい版を作る。どの版でも区間の集約を聞ける", "O(n) / O(log n) / O(log n)");
    RUN(auto st = make_persistent_segtree(std::vector<long long>{1, 2, 3, 4}, [](long long a, long long b) { return a + b; }, 0LL));
    RUN(int v1 = st.set(0, 0, 10));
    RUN(int v2 = st.set(v1, 3, 40));
    RUN(int v3 = st.set(0, 2, 0));  // 版 0 から分かれる
    SHOW(st.all_prod(0), "1 + 2 + 3 + 4");
    SHOW(st.all_prod(v2), "10 + 2 + 3 + 40");
    SHOW(st.prod(v3, 1, 4), "2 + 0 + 4");
    SHOW(st.prod_closed(v1, 0, 1), "10 + 2");

    api("【2】 update(v, i, f) / get(v, i)", "int / S", "版 v の a[i] を f(a[i]) にした新しい版 / 1 点の値", "O(log n)");
    RUN(int v4 = st.update(v2, 1, [](long long x) { return x * 100; }));
    SHOW(st.get(v4, 1), "");
    SHOW(st.get(v2, 1), "元の版はそのまま");

    api("【3】 max_right(v, l, pred)", "int", "版 v で pred(prod(v, l, r)) が true の最大の r", "O(log n)");
    SHOW(st.max_right(v2, 0, [](long long x) { return x <= 15; }), "10 + 2 + 3 = 15 まで");

    api("【4】 使い方の例：時刻ごとの区間の和", "—", "時刻 t の更新で版 t を作っておけば、あとから「時刻 t のときの [l, r) の和」を聞ける", "—");
    code("std::vector<int> at(T + 1); at[0] = 0;  for (t = 1..T) at[t] = st.set(at[t − 1], i_t, x_t);");
    code("答え = st.prod(at[t], l, r);");
}

}  // namespace pseg_demo

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
    pseg_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pseg_test::benchmark();    // 速度計測
    pseg_demo::demo();         // 使い方デモ
    return 0;
}
