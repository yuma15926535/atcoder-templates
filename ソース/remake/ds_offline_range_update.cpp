// =============================================================================
//  【制限付き最速版】区間更新の後読み（区間代入・区間 chmin / chmax を予約して、最後にまとめて各位置の値にする）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : 更新の途中で値を読めない（更新を全部終えてから、最後にまとめて読む）。区間の和・最小などの区間の質問はない
//    できること           : 区間 [l, r) に「代入（後に来たものが勝つ）」「chmin / chmax（など冪等な演算）」を何回でも予約し、
//                           build() で全位置の最終的な値を出す。どの代入がその位置を決めたかもわかる
//    一般版との比較       : RangeAssignPointGetTemplate.txt・遅延セグ木（lazy_segtree.cpp）は 1 回 O(log N)
//                           → 代入は全体で O((N + Q) α(N))、chmin / chmax は 1 回 O(1) + 最後に O(N log N)
//    制限に合わないとき   : 更新と質問がまざるなら RangeAssignPointGetTemplate.txt / セグメント木/追加/lazy_segtree.cpp
//
//  OfflineRangeAssign : 後ろの代入から順に見て「まだ決まっていない位置」だけを塗る。決まった位置は Union-Find で飛ばすので、
//                       各位置を塗るのは 1 回だけ
//  OfflineRangeIdempotent : 区間を長さ 2^k の 2 つの区間で覆って表に書き込み（Sparse Table を逆向きに使う）、
//                       build() で長い方から半分ずつに配る。op は min・max・gcd・AND・OR のように
//                       結合的・可換・冪等（op(x, x) = x）なもの
//    例: 「区間を色 c で塗る」を Q 回 → 最後の色、「区間の上限を x にする」を Q 回 → 各位置の値
//
//  機能一覧（位置は 0〜N−1。区間は半開区間 [l, r)。閉区間は _closed 版）                          計算量
//    OfflineRangeAssign<T> ra(n, init) / ra(a)   全部 init / 配列 a から                                    O(N)
//      ra.assign(l, r, x) / assign_closed(l, r, x)   [l, r) を x にする（予約。呼んだ順に番号 0, 1, 2, …）    O(1)
//      ra.build()                     全位置の値（後の assign が勝つ）を作って返す。何度でも呼べる              O((N + Q) α(N))
//      ra.get(i)                      build 後：位置 i の値                                                     O(1)
//      ra.last_query(i)               build 後：位置 i を決めた assign の番号（どれも塗っていなければ −1）       O(1)
//      ra.to_vector() / ra.query_count() / ra.size() / ra.clear()（予約を全部消す）
//    OfflineRangeIdempotent<T, Op>    冪等な演算の版（下の make_… で作る）
//      make_offline_range_chmin(a) / make_offline_range_chmax(a)   a[i] = min(a[i], x) / max の版
//      make_offline_range_apply(a, op, e)   一般の冪等な op（e は op(e, x) = x となる値）
//      (n, init) を渡す形もある（例: make_offline_range_chmin(n, init)）
//      oi.apply(l, r, x) / apply_closed(l, r, x)   [l, r) の各 a[i] を op(a[i], x) にする（予約）            O(1)
//      oi.build()                     全位置の値を作って返す（予約は反映して消える。続けて apply してよい）       O(N log N)
//      oi.get(i) / oi.to_vector() / oi.size()
//
//  注意: OfflineRangeIdempotent の表は使った長さの段だけ作る（最大 N log N 個の T。N = 2×10^5・long long で約 30 MB）。
//        メモリが厳しいときの chmin は OfflineRangeAssign で「x の大きい順に assign」して、最後に元の値との min を取ればよい
//        （小さい x が後から塗るので勝つ。ソートの分 O(Q log Q) かかる）
//  速さの目安: N = 10^6・assign 10^6 回 + build 約 0.03 秒（一般版の遅延セグ木（区間代入・1 点取得）は同じ量で 約 0.13 秒）、
//              N = 2×10^5・chmin 10^6 回 + build 約 0.04 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    oru_test::random_test()  … 1 つずつ書き込む愚直解との突き合わせ
//    oru_test::benchmark()    … 速度計測（一般版の遅延セグ木との比較つき）
//    oru_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
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
namespace oru_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
inline int floor_log2(unsigned x) {  // x ≥ 1
    return 31 - __builtin_clz(x);
}
}  // namespace oru_internal

// ---- 区間代入（後に来たものが勝つ）を予約して、最後にまとめて各位置の値にする ----
template <class T>
struct OfflineRangeAssign {
    explicit OfflineRangeAssign(int n = 0, T init = T()) : n_(n), init_(n < 0 ? 0 : n, init) {
        if (n < 0) oru_internal::fail("OfflineRangeAssign", "n は 0 以上");
    }
    explicit OfflineRangeAssign(const std::vector<T>& a) : n_((int)a.size()), init_(a) {}

    void assign(int l, int r, T x) {  // [l, r) を x に（予約）
        if (l < 0 || r > n_ || l > r) oru_internal::fail("OfflineRangeAssign::assign", "0 ≤ l ≤ r ≤ n の半開区間 [l, r)");
        qs_.push_back({l, r, std::move(x)});
    }
    void assign_closed(int l, int r, T x) { assign(l, r + 1, std::move(x)); }  // [l, r] を x に（予約）

    const std::vector<T>& build() {  // 後ろの assign から、まだ決まっていない位置だけ塗る
        val_ = init_;
        who_.assign(n_, -1);
        // Union-Find（union by size + 経路半分化）。各集合は「塗り終わった位置の連続」+「その右のまだの位置」で、根の alive がその右のまだの位置
        std::vector<int> par(n_ + 1, -1), alive(n_ + 1);
        std::iota(alive.begin(), alive.end(), 0);
        auto find = [&](int x) {
            while (par[x] >= 0) {
                if (par[par[x]] >= 0) par[x] = par[par[x]];
                x = par[x];
            }
            return x;
        };
        int rest = n_;
        for (int k = (int)qs_.size() - 1; k >= 0 && rest > 0; --k) {
            const Q& q = qs_[k];
            int i = alive[find(q.l)];
            while (i < q.r) {
                val_[i] = q.x, who_[i] = k, --rest;
                int a = find(i), b = find(i + 1);  // i を右と併合（i は塗り終わり）
                const int nr = alive[b];
                if (par[a] > par[b]) std::swap(a, b);
                par[a] += par[b], par[b] = a, alive[a] = nr;
                i = nr;
            }
        }
        built_ = true;
        return val_;
    }
    const T& get(int i) const {
        check_built("get");
        if (i < 0 || i >= n_) oru_internal::fail("OfflineRangeAssign::get", "位置が範囲外");
        return val_[i];
    }
    int last_query(int i) const {  // 位置 i を決めた assign の番号（呼んだ順に 0, 1, 2, …）。どれも塗っていなければ −1
        check_built("last_query");
        if (i < 0 || i >= n_) oru_internal::fail("OfflineRangeAssign::last_query", "位置が範囲外");
        return who_[i];
    }
    const std::vector<T>& to_vector() const {
        check_built("to_vector");
        return val_;
    }
    int query_count() const { return (int)qs_.size(); }
    int size() const { return n_; }
    void clear() { qs_.clear(), built_ = false; }  // 予約を全部消す（初期値はそのまま）

   private:
    struct Q {
        int l, r;
        T x;
    };
    int n_;
    std::vector<T> init_, val_;
    std::vector<int> who_;
    std::vector<Q> qs_;
    bool built_ = false;
    void check_built(const char* f) const {
        if (!built_) oru_internal::fail(f, "build() の前に呼んだ");
    }
};

// ---- 冪等な演算（min・max・gcd・AND・OR）を区間に予約して、最後にまとめて各位置の値にする ----
template <class T, class Op>
struct OfflineRangeIdempotent {
    OfflineRangeIdempotent(const std::vector<T>& a, Op op, T e) : n_((int)a.size()), a_(a), op_(op), e_(e) {
        tab_.resize(n_ == 0 ? 1 : oru_internal::floor_log2((unsigned)n_) + 1);
    }
    void apply(int l, int r, const T& x) {  // [l, r) の各 a[i] を op(a[i], x) に（予約）
        if (l < 0 || r > n_ || l > r) oru_internal::fail("OfflineRangeIdempotent::apply", "0 ≤ l ≤ r ≤ n の半開区間 [l, r)");
        if (l == r) return;
        const int k = oru_internal::floor_log2((unsigned)(r - l));
        std::vector<T>& t = level(k);
        t[l] = op_(t[l], x);
        const int j = r - (1 << k);
        t[j] = op_(t[j], x);
        built_ = false;
    }
    void apply_closed(int l, int r, const T& x) { apply(l, r + 1, x); }
    const std::vector<T>& build() {  // 長い段から半分ずつ配る。予約は a に反映して消える
        for (int k = (int)tab_.size() - 1; k >= 1; --k) {
            if (tab_[k].empty()) continue;
            std::vector<T>& lo = level(k - 1);
            const std::vector<T>& hi = tab_[k];
            const int h = 1 << (k - 1);
            for (int i = 0; i < (int)hi.size(); ++i) {
                lo[i] = op_(lo[i], hi[i]);
                lo[i + h] = op_(lo[i + h], hi[i]);
            }
            std::vector<T>().swap(tab_[k]);
        }
        if (!tab_.empty() && !tab_[0].empty()) {
            for (int i = 0; i < n_; ++i) a_[i] = op_(a_[i], tab_[0][i]);
            std::vector<T>().swap(tab_[0]);
        }
        built_ = true;
        return a_;
    }
    const T& get(int i) const {
        if (!built_) oru_internal::fail("OfflineRangeIdempotent::get", "build() の前に呼んだ（apply の後は build し直す）");
        if (i < 0 || i >= n_) oru_internal::fail("OfflineRangeIdempotent::get", "位置が範囲外");
        return a_[i];
    }
    const std::vector<T>& to_vector() const {
        if (!built_) oru_internal::fail("OfflineRangeIdempotent::to_vector", "build() の前に呼んだ（apply の後は build し直す）");
        return a_;
    }
    int size() const { return n_; }

   private:
    int n_;
    std::vector<T> a_;
    Op op_;
    T e_;
    std::vector<std::vector<T>> tab_;  // tab_[k][i] = [i, i + 2^k) に予約した値（使った段だけ作る）
    bool built_ = true;
    std::vector<T>& level(int k) {
        if (tab_[k].empty()) tab_[k].assign(n_ - (1 << k) + 1, e_);
        return tab_[k];
    }
};
template <class T, class Op>
OfflineRangeIdempotent<T, Op> make_offline_range_apply(const std::vector<T>& a, Op op, T e) {
    return OfflineRangeIdempotent<T, Op>(a, op, e);
}
template <class T, class Op>
OfflineRangeIdempotent<T, Op> make_offline_range_apply(int n, T init, Op op, T e) {
    return OfflineRangeIdempotent<T, Op>(std::vector<T>(n, init), op, e);
}
template <class T>
auto make_offline_range_chmin(const std::vector<T>& a) {
    auto op = [](const T& x, const T& y) { return y < x ? y : x; };
    return OfflineRangeIdempotent<T, decltype(op)>(a, op, std::numeric_limits<T>::max());
}
template <class T>
auto make_offline_range_chmin(int n, T init) { return make_offline_range_chmin(std::vector<T>(n, init)); }
template <class T>
auto make_offline_range_chmax(const std::vector<T>& a) {
    auto op = [](const T& x, const T& y) { return x < y ? y : x; };
    return OfflineRangeIdempotent<T, decltype(op)>(a, op, std::numeric_limits<T>::lowest());
}
template <class T>
auto make_offline_range_chmax(int n, T init) { return make_offline_range_chmax(std::vector<T>(n, init)); }
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace oru_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define ORU_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_offline_range_update] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 一般版の代わり：区間代入・1 点取得の遅延セグ木（時刻つき・push なし。RangeAssignPointGetTemplate.txt と同じやり方）
struct NaiveAssignSeg {
    int sz;
    std::vector<std::pair<int, ll>> t;  // (時刻, 値)
    explicit NaiveAssignSeg(int n) : sz(1) {
        while (sz < n) sz <<= 1;
        t.assign(2 * sz, {-1, 0});
    }
    void assign(int l, int r, ll x, int time) {
        for (l += sz, r += sz; l < r; l >>= 1, r >>= 1) {
            if (l & 1) t[l++] = {time, x};
            if (r & 1) t[--r] = {time, x};
        }
    }
    ll get(int i, ll init) const {
        std::pair<int, ll> best{-1, init};
        for (i += sz; i >= 1; i >>= 1)
            if (t[i].first > best.first) best = t[i];
        return best.second;
    }
};

static void random_test() {
    std::mt19937_64 rng(20261101);
    long long cases = 0;
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 30);
        std::vector<ll> a(n);
        for (auto& v : a) v = rnd(rng, -20, 20);
        // 代入
        {
            OfflineRangeAssign<ll> ra(a);
            std::vector<ll> b = a;
            std::vector<int> who(n, -1);
            const int q = (int)rnd(rng, 0, 30);
            for (int k = 0; k < q; ++k) {
                int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
                if (l > r) std::swap(l, r);
                const ll x = rnd(rng, -50, 50);
                if (rnd(rng, 0, 1) && l < r) {
                    ra.assign_closed(l, r - 1, x);
                } else {
                    ra.assign(l, r, x);
                }
                for (int i = l; i < r; ++i) b[i] = x, who[i] = k;
                if (rnd(rng, 0, 3) == 0) {  // 途中で build しても続けられる
                    ORU_CHECK(ra.build() == b);
                }
            }
            ORU_CHECK(ra.build() == b);
            ORU_CHECK(ra.to_vector() == b && ra.query_count() == q && ra.size() == n);
            for (int i = 0; i < n; ++i) ORU_CHECK(ra.get(i) == b[i] && ra.last_query(i) == who[i]);
            ra.clear();
            ORU_CHECK(ra.build() == a);
        }
        // chmin / chmax / gcd
        {
            auto mn = make_offline_range_chmin(a);
            auto mx = make_offline_range_chmax(a);
            std::vector<ll> g0(n);
            for (auto& v : g0) v = rnd(rng, 0, 60);
            auto gcd_op = [](ll x, ll y) { return std::gcd(x, y); };
            auto gd = make_offline_range_apply(g0, gcd_op, 0LL);
            std::vector<ll> bmn = a, bmx = a, bg = g0;
            const int q = (int)rnd(rng, 0, 40);
            for (int k = 0; k < q; ++k) {
                int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
                if (l > r) std::swap(l, r);
                const ll x = rnd(rng, -30, 30), y = rnd(rng, 0, 60);
                mn.apply(l, r, x), mx.apply(l, r, x);
                if (l < r) gd.apply_closed(l, r - 1, y);
                for (int i = l; i < r; ++i) bmn[i] = std::min(bmn[i], x), bmx[i] = std::max(bmx[i], x), bg[i] = std::gcd(bg[i], y);
                if (rnd(rng, 0, 4) == 0) ORU_CHECK(mn.build() == bmn);
            }
            ORU_CHECK(mn.build() == bmn && mx.build() == bmx && gd.build() == bg);
            for (int i = 0; i < n; ++i) ORU_CHECK(mn.get(i) == bmn[i] && mx.to_vector()[i] == bmx[i]);
            ORU_CHECK(mn.size() == n);
        }
        ++cases;
    }
    // (n, init) で作る形
    {
        OfflineRangeAssign<int> ra(5, 7);
        ra.assign(1, 3, 2);
        ORU_CHECK(ra.build() == std::vector<int>({7, 2, 2, 7, 7}));
        auto mn = make_offline_range_chmin(5, 10);
        mn.apply(0, 5, 4);
        ORU_CHECK(mn.build() == std::vector<int>(5, 4));
    }
    std::cout << "[ds_offline_range_update] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_offline_range_update] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    {
        const int N = TEST_SCALE == 1 ? 100000 : 1000000, Q = N;
        std::vector<std::array<int, 3>> qs(Q);
        for (auto& q : qs) {
            int l = (int)rnd(rng, 0, N), r = (int)rnd(rng, 0, N);
            if (l > r) std::swap(l, r);
            q = {l, r, (int)rnd(rng, 0, 1000000000)};
        }
        auto t0 = std::chrono::steady_clock::now();
        OfflineRangeAssign<ll> ra(N, 0);
        for (auto& q : qs) ra.assign(q[0], q[1], q[2]);
        const auto& v = ra.build();
        ll chk = 0;
        for (int i = 0; i < N; i += 7) chk += v[i];
        const ll ms1 = elapsed_ms(t0);
        t0 = std::chrono::steady_clock::now();
        NaiveAssignSeg sg(N);
        for (int k = 0; k < Q; ++k) sg.assign(qs[k][0], qs[k][1], qs[k][2], k);
        ll chk2 = 0;
        for (int i = 0; i < N; i += 7) chk2 += sg.get(i, 0);
        const ll ms2 = elapsed_ms(t0);
        ORU_CHECK(chk == chk2);
        std::cout << "[ds_offline_range_update] benchmark: OfflineRangeAssign N = " << N << "・assign " << Q << " 回 + build " << ms1
                  << " ms（一般版の遅延セグ木 " << ms2 << " ms）  (check " << chk % 1000 << ")\n";
    }
    {
        const int N = TEST_SCALE == 1 ? 50000 : 200000, Q = TEST_SCALE == 1 ? 200000 : 1000000;
        auto t0 = std::chrono::steady_clock::now();
        auto mn = make_offline_range_chmin(N, (ll)4e18);
        for (int k = 0; k < Q; ++k) {
            int l = (int)rnd(rng, 0, N), r = (int)rnd(rng, 0, N);
            if (l > r) std::swap(l, r);
            mn.apply(l, r, rnd(rng, 0, 1000000000));
        }
        const auto& v = mn.build();
        ll chk = 0;
        for (int i = 0; i < N; i += 7) chk += v[i];
        std::cout << "[ds_offline_range_update] benchmark: chmin N = " << N << "・apply " << Q << " 回 + build " << elapsed_ms(t0) << " ms  (check "
                  << chk % 1000 << ")\n";
    }
}

}  // namespace oru_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace oru_demo {

using ll = long long;
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
#define SHOW(expr, comment) oru_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    oru_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】区間更新の後読み（ds_offline_range_update） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 更新の途中で値を読めない（全部予約してから最後に build で読む）\n"
              << "■ できること: 区間代入（後勝ち）・区間 chmin / chmax / gcd などを予約 → 全位置の最終的な値\n"
              << "■ API 一覧\n"
              << "  OfflineRangeAssign<T> ra(n, init) / ra(a): assign / assign_closed / build / get / last_query / to_vector / query_count / size / clear\n"
              << "  make_offline_range_chmin(a) / make_offline_range_chmax(a) / make_offline_range_apply(a, op, e): apply / apply_closed / build / get / to_vector / size\n";

    api("【1】 OfflineRangeAssign<T> ra(n, init)  assign(l, r, x) → build()", "const vector<T>&", "[l, r) を x にする予約を全部反映した、全位置の値（後の assign が勝つ）",
        "assign O(1)、build O((N + Q) α(N))");
    RUN(OfflineRangeAssign<int> ra(10, 0));
    RUN(ra.assign(0, 6, 1); ra.assign(3, 9, 2); ra.assign(5, 7, 3));
    SHOW(ra.build(), "[0,6)→1、[3,9)→2、[5,7)→3 の順に塗った結果");
    SHOW(ra.get(4), "位置 4 の値");
    SHOW(ra.last_query(4), "位置 4 を決めたのは 1 番目（0 始まり）の assign");
    SHOW(ra.last_query(9), "どれも塗っていない → −1");
    SHOW(ra.query_count(), "予約した assign の数");

    api("【2】 assign_closed(l, r, x) / build を何度でも / clear()", "void / const vector<T>& / void", "閉区間 [l, r] で予約 / 続けて予約して build し直せる / 予約を全部消す（初期値に戻る）",
        "O(1) / O((N + Q) α(N)) / O(1)");
    RUN(ra.assign_closed(8, 9, 5));
    SHOW(ra.build(), "[8, 9] を 5 に（さっきの予約も残っている）");
    RUN(ra.clear());
    SHOW(ra.build(), "予約を消したので全部初期値 0");

    api("【3】 make_offline_range_chmin(a)  apply(l, r, x) → build()", "const vector<T>&", "[l, r) の各 a[i] を min(a[i], x) にする予約を反映した値",
        "apply O(1)、build O(N log N)");
    RUN(std::vector<ll> a = {5, 8, 3, 9, 7, 6});
    RUN(auto mn = make_offline_range_chmin(a));
    RUN(mn.apply(0, 4, 6); mn.apply(2, 6, 4); mn.apply_closed(5, 5, 1));
    SHOW(mn.build(), "[0,4) の上限 6、[2,6) の上限 4、[5,5] の上限 1");
    SHOW(mn.get(1), "");
    RUN(mn.apply(0, 6, 2));
    SHOW(mn.build(), "続けて apply → build（前の結果に重ねる）");

    api("【4】 make_offline_range_chmax(a) / make_offline_range_apply(a, op, e)", "OfflineRangeIdempotent（auto で受ける）", "max の版 / 一般の冪等な op（gcd・AND・OR など。e は op(e, x) = x となる値）",
        "apply O(1)、build O(N log N)");
    RUN(auto mx = make_offline_range_chmax(4, 0LL));
    RUN(mx.apply(1, 3, 5); mx.apply(0, 2, 2));
    SHOW(mx.build(), "(n, init) で作る形。[1,3) の下限 5、[0,2) の下限 2");
    RUN(auto gd = make_offline_range_apply(std::vector<ll>{12, 18, 24, 30}, [](ll x, ll y) { return std::gcd(x, y); }, 0LL));
    RUN(gd.apply(0, 3, 8); gd.apply(2, 4, 9));
    SHOW(gd.build(), "[0,3) と 8 の gcd、[2,4) と 9 の gcd");

    api("【5】 メモリが厳しいときの chmin", "（使い方の例）", "OfflineRangeAssign で x の大きい順に assign → 最後に元の値との min（小さい x が後から塗るので勝つ）",
        "O(Q log Q + (N + Q) α(N))");
    {
        std::vector<ll> b = {5, 8, 3, 9, 7, 6};
        std::vector<std::array<ll, 3>> qs = {{0, 4, 6}, {2, 6, 4}, {5, 6, 1}};
        std::sort(qs.begin(), qs.end(), [](auto& p, auto& q) { return p[2] > q[2]; });
        OfflineRangeAssign<ll> r2((int)b.size(), (ll)4e18);
        for (auto& q : qs) r2.assign((int)q[0], (int)q[1], q[2]);
        auto v = r2.build();
        for (int i = 0; i < (int)b.size(); ++i) b[i] = std::min(b[i], v[i]);
        note("【3】の最初の 3 回と同じ結果 → " + to_s(b));
    }
}

}  // namespace oru_demo

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
    oru_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    oru_test::benchmark();    // 速度計測
    oru_demo::demo();         // 使い方デモ
    return 0;
}
