// =============================================================================
//  【制限付き最速版】等差数列・多項式のいもす（区間に a, a+d, a+2d, … や多項式を足す予約を O(1) で、最後にまとめて各位置の値にする）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : 足している途中で値を読めない（全部足してから build() で読む）。区間の和などの質問はない（build 後に累積和を取る）
//    できること           : 区間 [l, r) に「同じ値」「等差数列 a, a + d, a + 2d, …」「(i − l) の多項式」を足す予約を何回でもして、
//                           build() で全位置の値を出す（build は何度でも呼べる）
//    一般版との比較       : 遅延セグ木で等差数列を足す（作用を自分で書く）と 1 回 O(log N)
//                           → この版は等差数列 1 回 O(1)、K 次の多項式 1 回 O(K^2)、build O(NK)
//    制限に合わないとき   : 足すのと読むのがまざるなら遅延セグ木（セグメント木/追加/lazy_segtree.cpp に等差数列の作用を書く）
//                           同じ値を足すだけなら DifferenceArrayTwoPointersTemplate.txt（1 階のいもす）でもよい
//
//  しくみ: (i − p) を j 個選ぶ二項係数 C(i − p, j) を位置 p から足すのは、差分を j + 1 回重ねた配列の 1 か所に書くのと同じ。
//          区間の終わりでは、続きの多項式を r の位置から引く。build() で累積和を K + 1 回取る
//    例: 「区間 [l, r) に 1, 2, 3, … を足す」を Q 回 → add_arith(l, r, 1, 1)、「(i − l)^2 を足す」→ add_poly(l, r, {0, 0, 1})
//
//  機能一覧（位置は 0〜N−1。区間は半開区間 [l, r)。閉区間は _closed 版。区間が外にはみ出した部分は捨てる）  計算量
//    PolyImos<T> im(n, K = 1)              長さ n、多項式の次数の上限 K（等差数列なら 1、同じ値だけなら 0）      O(NK)
//      im.add(l, r, x) / add_closed        [l, r) の全部に x を足す                                             O(1)
//      im.add_arith(l, r, a, d) / add_arith_closed   位置 l, l+1, … に a, a + d, a + 2d, … を足す（K ≥ 1）      O(1)
//      im.add_poly(l, r, c) / add_poly_closed        位置 i に c[0] + c[1](i − l) + c[2](i − l)^2 + … を足す（c の長さ ≤ K + 1）  O(K^2)
//      im.build()                          全位置の値を作って返す（予約は残るので、続けて足して何度でも呼べる）   O(NK)
//      im.get(i) / im.to_vector()          build 後：位置 i の値 / 全体                                          O(1)
//      im.clear() / im.size() / im.degree()   予約を全部消す / n / K
//
//  T は long long（ふつう）/ __int128 / modint（T(long long) で作れて +, −, * ができる型）。値は T で計算するので溢れに注意
//  注意: K は 18 まで（多項式を二項係数に直す係数 S(k, j) × j! が long long に収まる範囲）
//  速さの目安: N = 10^6・add_arith 10^6 回 + build 約 0.06 秒、N = 10^6・K = 3 の add_poly 10^6 回 + build 約 0.17 秒（どちらも乱数の生成込み）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    pia_test::random_test()  … 1 つずつ足す愚直解との突き合わせ
//    pia_test::benchmark()    … 速度計測
//    pia_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
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
namespace pia_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace pia_internal

template <class T = long long>
struct PolyImos {
    using Buf = std::array<T, 19>;  // 多項式の係数（次数 18 まで）
    explicit PolyImos(int n = 0, int degree = 1) : n_(n), K_(degree) {
        if (n < 0) pia_internal::fail("PolyImos", "n は 0 以上");
        if (degree < 0 || degree > 18) pia_internal::fail("PolyImos", "次数 K は 0〜18");
        E_.assign(K_ + 1, std::vector<T>(n_, T(0)));
        // sf_[k][j] = S(k, j) × j!（x^k = Σ_j sf_[k][j] C(x, j)）、bin_[k][t] = C(k, t)
        sf_.assign(K_ + 1, std::vector<long long>(K_ + 1, 0));
        bin_.assign(K_ + 1, std::vector<long long>(K_ + 1, 0));
        std::vector<std::vector<long long>> st(K_ + 1, std::vector<long long>(K_ + 1, 0));
        st[0][0] = 1;
        for (int k = 1; k <= K_; ++k)
            for (int j = 1; j <= k; ++j) st[k][j] = st[k - 1][j - 1] + (long long)j * st[k - 1][j];
        for (int k = 0; k <= K_; ++k) {
            long long f = 1;
            for (int j = 0; j <= k; ++j) {
                if (j) f *= j;
                sf_[k][j] = st[k][j] * f;
            }
            bin_[k][0] = 1;
            for (int t = 1; t <= k; ++t) bin_[k][t] = bin_[k - 1][t - 1] + (t <= k - 1 ? bin_[k - 1][t] : 0);
        }
    }

    void add(long long l, long long r, const T& x) {  // [l, r) の全部に x
        if (!clip(l, r)) return;
        put(0, l, x);
        put(0, r, T(0) - x);
    }
    void add_closed(long long l, long long r, const T& x) { add(l, r + 1, x); }
    void add_arith(long long l, long long r, const T& a, const T& d) {  // 位置 l + k に a + k d（l ≤ l + k < r）
        if (K_ < 1) pia_internal::fail("PolyImos::add_arith", "等差数列は次数 K ≥ 1 で作る");
        if (l >= r || r <= 0 || l >= n_) return;
        if (l < 0) {  // 左にはみ出した分は初項をずらす
            add_arith_in(0, r, a + T(-l) * d, d);
            return;
        }
        add_arith_in(l, r, a, d);
    }
    void add_arith_closed(long long l, long long r, const T& a, const T& d) { add_arith(l, r + 1, a, d); }
    void add_poly(long long l, long long r, const std::vector<T>& c) {  // 位置 i に Σ_k c[k] (i − l)^k
        if ((int)c.size() > K_ + 1) pia_internal::fail("PolyImos::add_poly", "係数の数が K + 1 を超える（作るときの K を大きくする）");
        if (l >= r || c.empty()) return;
        long long l2 = l, r2 = r;
        if (!clip(l2, r2)) return;
        const int k = (int)c.size();
        Buf c0, b, c1;
        for (int a = 0; a < k; ++a) c0[a] = c[a];
        if (l2 != l) shift(c0, k, l2 - l, c1), c0 = c1;  // 左にはみ出したら (i − l2) の式に直す
        to_binom(c0, k, b);
        for (int j = 0; j < k; ++j) put(j, l2 + j, b[j]);
        if (r2 < n_) {  // r から先は同じ多項式を引く
            shift(c0, k, r2 - l2, c1);
            to_binom(c1, k, b);
            for (int j = 0; j < k; ++j) put(j, r2 + j, T(0) - b[j]);
        }
    }
    void add_poly_closed(long long l, long long r, const std::vector<T>& c) { add_poly(l, r + 1, c); }

    const std::vector<T>& build() {  // V = P(E0 + P(E1 + P(E2 + …)))（P は累積和）
        val_ = E_[K_];
        for (int m = K_ - 1; m >= 0; --m) {
            prefix(val_);
            const std::vector<T>& e = E_[m];
            for (int i = 0; i < n_; ++i) val_[i] += e[i];
        }
        prefix(val_);
        built_ = true;
        return val_;
    }
    const T& get(int i) const {
        if (!built_) pia_internal::fail("PolyImos::get", "build() の前に呼んだ");
        if (i < 0 || i >= n_) pia_internal::fail("PolyImos::get", "位置が範囲外");
        return val_[i];
    }
    const std::vector<T>& to_vector() const {
        if (!built_) pia_internal::fail("PolyImos::to_vector", "build() の前に呼んだ");
        return val_;
    }
    void clear() {
        for (auto& e : E_) std::fill(e.begin(), e.end(), T(0));
        built_ = false;
    }
    int size() const { return n_; }
    int degree() const { return K_; }

   private:
    int n_, K_;
    std::vector<std::vector<T>> E_;  // E_[j][p] = 「位置 p − j から C(i − (p − j), j) を足す」の係数
    std::vector<T> val_;
    std::vector<std::vector<long long>> sf_, bin_;
    bool built_ = false;

    bool clip(long long& l, long long& r) const {
        if (l < 0) l = 0;
        if (r > n_) r = n_;
        return l < r;
    }
    void put(int j, long long p, const T& v) {
        if (0 <= p && p < n_) E_[j][p] += v;
    }
    void add_arith_in(long long l, long long r, const T& a, const T& d) {
        if (l >= n_) return;
        put(0, l, a);
        put(1, l + 1, d);
        if (r < n_) {
            put(0, r, T(0) - (a + T(r - l) * d));
            put(1, r + 1, T(0) - d);
        }
    }
    void shift(const Buf& c, int k, long long L, Buf& res) const {  // Σ c[a] (x + L)^a を x の式に
        Buf pw;
        pw[0] = T(1);
        for (int t = 1; t < k; ++t) pw[t] = pw[t - 1] * T(L);
        for (int t = 0; t < k; ++t) res[t] = T(0);
        for (int a = 0; a < k; ++a)
            for (int t = 0; t <= a; ++t) res[t] += c[a] * T(bin_[a][t]) * pw[a - t];
    }
    void to_binom(const Buf& c, int k, Buf& b) const {  // Σ c[a] x^a = Σ b[j] C(x, j)
        for (int j = 0; j < k; ++j) b[j] = T(0);
        for (int a = 0; a < k; ++a)
            for (int j = 0; j <= a; ++j) b[j] += c[a] * T(sf_[a][j]);
    }
    void prefix(std::vector<T>& v) const {
        for (int i = 1; i < n_; ++i) v[i] += v[i - 1];
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace pia_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define PIA_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[pre_imos_arith] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    std::mt19937_64 rng(20261101);
    long long cases = 0;
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 25), K = (int)rnd(rng, 0, 4);
        PolyImos<ll> im(n, K);
        std::vector<ll> b(n, 0);
        const int q = (int)rnd(rng, 0, 25);
        for (int t = 0; t < q; ++t) {
            ll l = rnd(rng, -5, n + 5), r = rnd(rng, -5, n + 5);
            if (l > r) std::swap(l, r);
            const int kind = (int)rnd(rng, 0, 5);
            if (kind == 0) {
                const ll x = rnd(rng, -9, 9);
                im.add(l, r, x);
                for (ll i = std::max(l, 0LL); i < std::min(r, (ll)n); ++i) b[i] += x;
            } else if (kind == 1) {
                const ll x = rnd(rng, -9, 9);
                im.add_closed(l, r, x);
                for (ll i = std::max(l, 0LL); i <= std::min(r, (ll)n - 1); ++i) b[i] += x;
            } else if (kind == 2 && K >= 1) {
                const ll a = rnd(rng, -9, 9), d = rnd(rng, -9, 9);
                im.add_arith(l, r, a, d);
                for (ll i = std::max(l, 0LL); i < std::min(r, (ll)n); ++i) b[i] += a + (i - l) * d;
            } else if (kind == 3 && K >= 1) {
                const ll a = rnd(rng, -9, 9), d = rnd(rng, -9, 9);
                im.add_arith_closed(l, r, a, d);
                for (ll i = std::max(l, 0LL); i <= std::min(r, (ll)n - 1); ++i) b[i] += a + (i - l) * d;
            } else {
                std::vector<ll> c((int)rnd(rng, 0, K + 1));
                for (auto& v : c) v = rnd(rng, -5, 5);
                auto f = [&](ll x) {
                    ll s = 0, p = 1;
                    for (ll v : c) s += v * p, p *= x;
                    return s;
                };
                if (kind == 5) {
                    im.add_poly_closed(l, r, c);
                    for (ll i = std::max(l, 0LL); i <= std::min(r, (ll)n - 1); ++i) b[i] += f(i - l);
                } else {
                    im.add_poly(l, r, c);
                    for (ll i = std::max(l, 0LL); i < std::min(r, (ll)n); ++i) b[i] += f(i - l);
                }
            }
            if (rnd(rng, 0, 4) == 0) PIA_CHECK(im.build() == b);
        }
        PIA_CHECK(im.build() == b && im.to_vector() == b);
        for (int i = 0; i < n; ++i) PIA_CHECK(im.get(i) == b[i]);
        PIA_CHECK(im.size() == n && im.degree() == K);
        im.clear();
        PIA_CHECK(im.build() == std::vector<ll>(n, 0));
        ++cases;
    }
    // 大きい次数（K = 18）と符号なしの剰余（2^64）でも合う
    {
        const int n = 40, K = 18;
        PolyImos<unsigned long long> im(n, K);
        std::vector<unsigned long long> b(n, 0), c(K + 1);
        for (int k = 0; k <= K; ++k) c[k] = (unsigned long long)(k * 7 + 3);
        im.add_poly(5, 33, c);
        for (int i = 5; i < 33; ++i) {
            unsigned long long s = 0, p = 1, x = (unsigned long long)(i - 5);
            for (auto v : c) s += v * p, p *= x;
            b[i] += s;
        }
        PIA_CHECK(im.build() == b);
    }
    std::cout << "[pre_imos_arith] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[pre_imos_arith] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 100000 : 1000000, Q = N;
    {
        const auto t0 = std::chrono::steady_clock::now();
        PolyImos<ll> im(N, 1);
        for (int i = 0; i < Q; ++i) {
            ll l = rnd(rng, 0, N), r = rnd(rng, 0, N);
            if (l > r) std::swap(l, r);
            im.add_arith(l, r, rnd(rng, -1000, 1000), rnd(rng, -1000, 1000));
        }
        const auto& v = im.build();
        ll chk = 0;
        for (int i = 0; i < N; i += 7) chk += v[i];
        std::cout << "[pre_imos_arith] benchmark: N = " << N << "・add_arith " << Q << " 回 + build " << elapsed_ms(t0) << " ms  (check " << chk % 1000 << ")\n";
    }
    {
        const auto t0 = std::chrono::steady_clock::now();
        PolyImos<unsigned long long> im(N, 3);  // 3 次式は long long を溢れるので 2^64 の剰余で
        std::vector<unsigned long long> c(4);
        for (int i = 0; i < Q; ++i) {
            ll l = rnd(rng, 0, N), r = rnd(rng, 0, N);
            if (l > r) std::swap(l, r);
            for (auto& x : c) x = (unsigned long long)rnd(rng, 0, 9);
            im.add_poly(l, r, c);
        }
        const auto& v = im.build();
        unsigned long long chk = 0;
        for (int i = 0; i < N; i += 7) chk += v[i];
        std::cout << "[pre_imos_arith] benchmark: N = " << N << "・K = 3 の add_poly " << Q << " 回 + build " << elapsed_ms(t0) << " ms  (check " << chk % 1000
                  << ")\n";
    }
}

}  // namespace pia_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace pia_demo {

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
#define SHOW(expr, comment) pia_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    pia_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】等差数列・多項式のいもす（pre_imos_arith） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 足している途中で値を読めない（全部足してから build で読む）\n"
              << "■ できること: 区間に 同じ値・等差数列・(i − l) の多項式 を O(1)（多項式は O(K^2)）で足す予約 → 全位置の値\n"
              << "■ API 一覧\n"
              << "  PolyImos<T> im(n, K = 1): add / add_closed / add_arith / add_arith_closed / add_poly / add_poly_closed / build / get / to_vector / clear / size / degree\n";

    api("【1】 PolyImos<T> im(n, K)  add_arith(l, r, a, d) → build()", "const vector<T>&", "位置 l, l+1, …, r−1 に a, a+d, a+2d, … を足した結果",
        "add_arith O(1)、build O(NK)");
    RUN(PolyImos<ll> im(10, 1));
    RUN(im.add_arith(2, 6, 1, 1));
    SHOW(im.build(), "[2, 6) に 1, 2, 3, 4");
    RUN(im.add_arith(0, 10, 10, -1));
    SHOW(im.build(), "続けて [0, 10) に 10, 9, …, 1 を足す（予約は残る）");
    SHOW(im.get(3), "");

    api("【2】 add(l, r, x) / add_closed / add_arith_closed", "void", "[l, r) に同じ値 x / 閉区間 [l, r] の版", "O(1)");
    RUN(PolyImos<ll> im2(8, 1));
    RUN(im2.add(1, 4, 5); im2.add_closed(3, 5, 1); im2.add_arith_closed(5, 7, 100, 10));
    SHOW(im2.build(), "[1,4) に 5、[3,5] に 1、[5,7] に 100, 110, 120");
    RUN(im2.add_arith(-2, 3, 0, 1));
    SHOW(im2.build(), "左にはみ出しても OK（位置 −2 から 0, 1, 2, … のつもりで、位置 0 は 2）");

    api("【3】 add_poly(l, r, c) / add_poly_closed", "void", "位置 i に c[0] + c[1](i − l) + c[2](i − l)^2 + … を足す（c の長さ ≤ K + 1）", "O(K^2)");
    RUN(PolyImos<ll> im3(8, 2));
    RUN(im3.add_poly(1, 6, {0, 0, 1}));
    SHOW(im3.build(), "[1, 6) に (i − 1)^2 = 0, 1, 4, 9, 16");
    RUN(im3.add_poly_closed(0, 7, {1, 1}));
    SHOW(im3.build(), "[0, 7] に 1 + i（K 以下の次数なら何でも）");

    api("【4】 clear() / size() / degree()", "void / int / int", "予約を全部消す / 長さ / 次数の上限 K", "O(NK) / O(1) / O(1)");
    RUN(im3.clear());
    SHOW(im3.build(), "全部 0 に戻る");
    SHOW(im3.size(), "");
    SHOW(im3.degree(), "");
    note("類題の形: 「区間 [l, r) に 1, 2, 3, … を足す」を Q 回して最後に全部出力 → add_arith(l, r, 1, 1) を Q 回 → build()");
}

}  // namespace pia_demo

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
    pia_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pia_test::benchmark();    // 速度計測
    pia_demo::demo();         // 使い方デモ
    return 0;
}
