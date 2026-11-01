// =============================================================================
//  【制限付き最速版】単調な CHT（Convex Hull Trick：直線を傾きの順に追加し、x の順に最小値 / 最大値を聞く。ならし O(1)）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : 直線は傾きの順にしか追加できない（最小値なら傾きが減る順、最大値なら増える順。同じ傾きは可）。
//                           速い質問 query_inc / query_dec は x が単調（増える一方 / 減る一方）のときだけ。線分（x の範囲つき）は追加できない
//    できること           : 直線 y = a x + b の追加と、x での最小値（最大値）。x がばらばらなら query(x) で O(log n)
//    一般版との比較       : li_chao_tree.cpp は追加・質問 O(log n)（順番の制約なし・線分も可）
//                           → この版は追加・質問とも ならし O(1)（全体 O(N)）
//    制限に合わないとき   : 傾きや x の順番がばらばら・線分を足すなら li_chao_tree.cpp
//                           （x がばらばらなだけなら query(x)（O(log n)）でもよいが、オーダーは Li Chao と同じ）
//
//  しくみ: 傾きの順に並んだ直線の下側の包絡線（最小値をとる直線の列）を配列で持つ。新しい直線で要らなくなった直線は後ろから捨てる。
//          x が増える順に聞くなら、要らなくなった前の直線を捨てながら前から見る（各直線は 1 回しか捨てないので、ならし O(1)）
//    例: EDPC Z（Frog 3）… dp[i] = min_j (dp[j] + (h_i − h_j)^2) + C。h が増える順なので、傾き −2 h_j が減る順に追加し、x = h_i が増える順に聞く
//
//  機能一覧（T は long long（ふつう）/ __int128 / double などの符号つきの型。IS_MIN = true で最小値、false で最大値）  計算量
//    MonotoneCHTMin<T> ch / MonotoneCHTMax<T> ch      空（MonotoneCHT<T, IS_MIN> の別名。T の既定は long long）
//      ch.add_line(a, b)               直線 y = a x + b を追加（Min は a が前以下、Max は a が前以上。違えば止める）   ならし O(1)
//      ch.query_inc(x)                 x での最小値（Max は最大値）。x は前の query_inc 以上                         ならし O(1)
//      ch.query_dec(x)                 同じく。x は前の query_dec 以下（query_inc とはまぜられない）                ならし O(1)
//      ch.query(x)                     x の順番は自由（二分探索）。query_inc / query_dec も使うなら、前に聞いた x 以上（以下）だけ  O(log n)
//      ch.size() / ch.empty() / ch.clear()   包絡線に残っている直線の数 / 空か / 全部消す
//      MonotoneCHTMin<T>::EMPTY        直線がないときの値（Min は T の最大値、Max は最小値）
//  オーバーフロー注意: a x + b を T で計算する（|a x|, |b| ≤ 4e18 なら long long）。直線を捨てるかの判定は long long なら __int128 で正確、
//                      __int128 / double なら long double
//
//  速さの目安: 直線の追加と query_inc を 10^6 回ずつ 約 0.015 秒（一般版の Li Chao Tree は 約 0.11 秒）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mcht_test::random_test()  … 全部の直線を見る愚直解との突き合わせ
//    mcht_test::benchmark()    … 速度計測（Li Chao Tree との比較つき）
//    mcht_demo::demo()         … 全機能の使い方デモ（日本語の説明付き。EDPC Z の形も）
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
#include <type_traits>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace mcht_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
template <class T>
constexpr T max_value() {
    if constexpr (std::is_same_v<T, __int128>) {
        return (T)(((unsigned __int128)1 << 127) - 1);
    } else {
        return std::numeric_limits<T>::max();
    }
}
template <class T>
constexpr T lowest_value() {
    if constexpr (std::is_same_v<T, __int128>) {
        return -max_value<T>() - 1;
    } else {
        return std::numeric_limits<T>::lowest();
    }
}
// 傾き a1 > a2 > a3 の 3 本で、真ん中の直線 2 が最小値に一度も出てこないか
template <class T>
bool bad(T a1, T b1, T a2, T b2, T a3, T b3) {
    if constexpr (std::is_integral_v<T> && !std::is_same_v<T, __int128> && sizeof(T) <= 8) {
        return ((__int128)b3 - b1) * ((__int128)a1 - a2) <= ((__int128)b2 - b1) * ((__int128)a1 - a3);
    } else {
        return ((long double)b3 - (long double)b1) * ((long double)a1 - (long double)a2) <= ((long double)b2 - (long double)b1) * ((long double)a1 - (long double)a3);
    }
}
}  // namespace mcht_internal

template <class T = long long, bool IS_MIN = true>
struct MonotoneCHT {
    static constexpr T EMPTY = IS_MIN ? mcht_internal::max_value<T>() : mcht_internal::lowest_value<T>();

    void add_line(T a, T b) {  // y = a x + b（中では最小値の形に直して持つ：最大値なら符号を反転）
        if (!IS_MIN) a = -a, b = -b;
        if (has_last_ && a > last_a_) mcht_internal::fail("MonotoneCHT::add_line", "傾きの順番が違う（Min は傾きが減る順、Max は増える順に追加する）");
        has_last_ = true, last_a_ = a;
        if (size() > 0 && ls_.back().a == a) {
            if (ls_.back().b <= b) return;  // 同じ傾きで上にある直線は要らない
            ls_.pop_back();
        }
        while (size() >= 2) {
            const Line& l1 = ls_[ls_.size() - 2];
            const Line& l2 = ls_.back();
            if (!mcht_internal::bad(l1.a, l1.b, l2.a, l2.b, a, b)) break;
            ls_.pop_back();
        }
        ls_.push_back({a, b});
    }
    T query_inc(T x) {  // x は前の query_inc 以上
        set_mode(1, x);
        if (size() == 0) return EMPTY;
        while (size() >= 2 && f(ls_[head_ + 1], x) <= f(ls_[head_], x)) ++head_;
        return out(f(ls_[head_], x));
    }
    T query_dec(T x) {  // x は前の query_dec 以下
        set_mode(2, x);
        if (size() == 0) return EMPTY;
        while (size() >= 2 && f(ls_[ls_.size() - 2], x) <= f(ls_.back(), x)) ls_.pop_back();
        return out(f(ls_.back(), x));
    }
    T query(T x) const {  // x の順番は自由（二分探索）。query_inc / query_dec と一緒に使うなら、それまでに聞いた x より先だけ
        if ((mode_ == 1 && x < last_x_) || (mode_ == 2 && x > last_x_))
            mcht_internal::fail("MonotoneCHT::query", "query_inc（query_dec）で捨てた直線が要る x（前に聞いた x より小さい（大きい））は聞けない");
        if (size() == 0) return EMPTY;
        int lo = head_, hi = (int)ls_.size() - 1;
        while (lo < hi) {
            const int mid = (lo + hi) / 2;
            if (f(ls_[mid], x) <= f(ls_[mid + 1], x)) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        return out(f(ls_[lo], x));
    }
    int size() const { return (int)ls_.size() - head_; }
    bool empty() const { return size() == 0; }
    void clear() { ls_.clear(), head_ = 0, mode_ = 0, has_last_ = false; }

   private:
    struct Line {
        T a, b;
    };
    std::vector<Line> ls_;  // 有効なのは [head_, size)
    int head_ = 0, mode_ = 0;
    bool has_last_ = false;
    T last_a_{}, last_x_{};
    static T f(const Line& l, T x) { return l.a * x + l.b; }
    static T out(T v) { return IS_MIN ? v : -v; }
    void set_mode(int m, T x) {
        if (mode_ != 0 && mode_ != m) mcht_internal::fail("MonotoneCHT::query_inc / query_dec", "query_inc と query_dec はまぜられない");
        if (mode_ == m && (m == 1 ? x < last_x_ : x > last_x_))
            mcht_internal::fail("MonotoneCHT::query_inc / query_dec", "x の順番が違う（query_inc は増える順、query_dec は減る順。ばらばらなら query）");
        mode_ = m, last_x_ = x;
    }
};
template <class T = long long>
using MonotoneCHTMin = MonotoneCHT<T, true>;
template <class T = long long>
using MonotoneCHTMax = MonotoneCHT<T, false>;
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace mcht_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define MCHT_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[dp_monotone_cht] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 一般版：Li Chao Tree（x は [lo, hi) の整数。最小値）
struct SimpleLiChao {
    struct L {
        ll a, b;
        ll f(ll x) const { return a * x + b; }
    };
    ll lo, hi;
    std::vector<L> t;
    std::vector<char> has;
    std::vector<int> lc, rc;
    SimpleLiChao(ll lo_, ll hi_) : lo(lo_), hi(hi_) { newnode(); }
    int newnode() {
        t.push_back({0, 0}), has.push_back(0), lc.push_back(-1), rc.push_back(-1);
        return (int)t.size() - 1;
    }
    void add(L ln) {
        int k = 0;
        ll l = lo, r = hi;
        while (true) {
            if (!has[k]) {
                t[k] = ln, has[k] = 1;
                return;
            }
            const ll m = l + (r - l) / 2;
            const bool lef = ln.f(l) < t[k].f(l), mid = ln.f(m) < t[k].f(m);
            if (mid) std::swap(t[k], ln);
            if (r - l == 1) return;
            if (lef != mid) {
                if (lc[k] < 0) {
                    const int c = newnode();
                    lc[k] = c;
                }
                k = lc[k], r = m;
            } else {
                if (rc[k] < 0) {
                    const int c = newnode();
                    rc[k] = c;
                }
                k = rc[k], l = m;
            }
        }
    }
    ll query(ll x) const {
        ll res = std::numeric_limits<ll>::max(), l = lo, r = hi;
        int k = 0;
        while (k >= 0) {
            if (has[k]) res = std::min(res, t[k].f(x));
            const ll m = l + (r - l) / 2;
            if (x < m) {
                k = lc[k], r = m;
            } else {
                k = rc[k], l = m;
            }
        }
        return res;
    }
};

static void random_test() {
    std::mt19937_64 rng(20261101);
    long long cases = 0;
    for (int it = 0; it < 600 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 25), q = 30;
        const ll A = rnd(rng, 1, 20), B = rnd(rng, 1, 1000);
        std::vector<ll> as(n);
        for (auto& a : as) a = rnd(rng, -A, A);
        std::sort(as.rbegin(), as.rend());  // 減る順（Min 用）
        std::vector<std::pair<ll, ll>> lines;
        MonotoneCHTMin<ll> mn_inc, mn_dec, mn_any;
        MonotoneCHTMax<ll> mx_inc;
        std::vector<ll> xi(q), xd(q);
        for (auto& x : xi) x = rnd(rng, -30, 30);
        xd = xi;
        std::sort(xi.begin(), xi.end());
        std::sort(xd.rbegin(), xd.rend());
        // 追加と質問をまぜる
        int qi = 0;
        for (int i = 0; i <= n; ++i) {
            if (i < n) {
                const ll b = rnd(rng, -B, B);
                mn_inc.add_line(as[i], b), mn_dec.add_line(as[i], b), mn_any.add_line(as[i], b);
                mx_inc.add_line(-as[i], -b);  // 傾きは増える順
                lines.push_back({as[i], b});
            }
            const int k = (int)rnd(rng, 0, 2);
            for (int t = 0; t < k && qi < q; ++t, ++qi) {
                auto best = [&](ll x) {
                    ll r = MonotoneCHTMin<ll>::EMPTY;
                    for (auto [a, b] : lines) r = std::min(r, a * x + b);
                    return r;
                };
                MCHT_CHECK(mn_inc.query_inc(xi[qi]) == best(xi[qi]));
                MCHT_CHECK(mn_dec.query_dec(xd[qi]) == best(xd[qi]));
                const ll xa = xi[qi] + rnd(rng, 0, 5), xb = xd[qi] - rnd(rng, 0, 5);  // 前に聞いた x より先なら query もまぜてよい
                MCHT_CHECK(mn_inc.query(xa) == best(xa) && mn_dec.query(xb) == best(xb));
                const ll xr = rnd(rng, -30, 30);
                MCHT_CHECK(mn_any.query(xr) == best(xr));
                const ll mxv = mx_inc.query_inc(xi[qi]);
                MCHT_CHECK(lines.empty() ? mxv == MonotoneCHTMax<ll>::EMPTY : mxv == -best(xi[qi]));
            }
        }
        MCHT_CHECK(mn_any.empty() == lines.empty() && mn_any.size() <= (int)lines.size());
        mn_any.clear();
        MCHT_CHECK(mn_any.empty() && mn_any.query(0) == MonotoneCHTMin<ll>::EMPTY);
        ++cases;
    }
    // 大きい値（判定は __int128 で正確）と double
    {
        MonotoneCHTMin<ll> ch;  // 判定の積は 10^27 くらいになる（long long では溢れる）
        const ll E18 = 1000000000000000000LL, E9 = 1000000000LL;
        ch.add_line(E9, 0), ch.add_line(0, -E18), ch.add_line(-E9, 0);
        MCHT_CHECK(ch.size() == 3 && ch.query(0) == -E18 && ch.query(2 * E9) == -2 * E18 && ch.query(-2 * E9) == -2 * E18);
        MonotoneCHTMax<double> cd;
        cd.add_line(-0.5, 1.0), cd.add_line(0.25, 0.0), cd.add_line(2.0, -3.0);
        MCHT_CHECK(std::abs(cd.query_inc(0.0) - 1.0) < 1e-12 && std::abs(cd.query_inc(10.0) - 17.0) < 1e-12);
    }
    std::cout << "[dp_monotone_cht] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[dp_monotone_cht] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    {
        const int N = TEST_SCALE == 1 ? 1000000 : 10000000;
        const auto t0 = std::chrono::steady_clock::now();
        MonotoneCHTMin<ll> ch;
        ll chk = 0, a = 1000000000, x = -1000000000;
        for (int i = 0; i < N; ++i) {
            a -= rnd(rng, 0, 100), x += rnd(rng, 0, 100);
            ch.add_line(a, rnd(rng, -1000000000, 1000000000));
            chk += ch.query_inc(x) & 1023;
        }
        std::cout << "[dp_monotone_cht] benchmark: add_line と query_inc を " << N << " 回ずつ " << elapsed_ms(t0) << " ms  (check " << chk % 1000 << ")\n";
    }
    {
        const int N = TEST_SCALE == 1 ? 100000 : 1000000;
        std::vector<std::array<ll, 3>> ops(N);
        ll a = 1000000, x = -100000000;
        for (auto& o : ops) a -= rnd(rng, 0, 2), x += rnd(rng, 0, 200), o = {a, rnd(rng, -1000000000, 1000000000), x};
        auto t0 = std::chrono::steady_clock::now();
        MonotoneCHTMin<ll> ch;
        ll c1 = 0;
        for (auto& o : ops) ch.add_line(o[0], o[1]), c1 += ch.query_inc(o[2]) & 1023;
        const ll ms1 = elapsed_ms(t0);
        t0 = std::chrono::steady_clock::now();
        SimpleLiChao lc(-100000000, 300000000);
        ll c2 = 0;
        for (auto& o : ops) lc.add({o[0], o[1]}), c2 += lc.query(o[2]) & 1023;
        const ll ms2 = elapsed_ms(t0);
        MCHT_CHECK(c1 == c2);
        std::cout << "[dp_monotone_cht] benchmark: 同じ " << N << " 回ずつで MonotoneCHT " << ms1 << " ms / 一般版（Li Chao Tree）" << ms2 << " ms\n";
    }
}

}  // namespace mcht_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace mcht_demo {

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
#define SHOW(expr, comment) mcht_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    mcht_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】単調な CHT（dp_monotone_cht） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 直線は傾きの順にしか追加できない（Min は減る順、Max は増える順）・速い質問は x が単調なときだけ・線分は不可\n"
              << "■ できること: 直線の追加と x での最小値 / 最大値が ならし O(1)（x がばらばらなら query で O(log n)）\n"
              << "■ API 一覧\n"
              << "  MonotoneCHTMin<T> / MonotoneCHTMax<T>: add_line / query_inc / query_dec / query / size / empty / clear / EMPTY\n";

    api("【1】 MonotoneCHTMin<T> ch  add_line(a, b) / query_inc(x)", "void / T", "直線 y = a x + b を追加（傾きは減る順）/ x（増える順）での最小値", "ならし O(1)");
    RUN(MonotoneCHTMin<ll> ch);
    RUN(ch.add_line(2, 0); ch.add_line(0, 3); ch.add_line(-1, 7));
    SHOW(ch.query_inc(-1), "min(−2, 3, 8)");
    SHOW(ch.query_inc(2), "min(4, 3, 5)");
    SHOW(ch.query_inc(10), "min(20, 3, −3)");
    SHOW(ch.size(), "x = 10 より前の直線はもう要らないので捨てた");

    api("【2】 query_dec(x) / query(x)", "T", "x が減る順の質問 / 順番が自由な質問（二分探索）", "ならし O(1) / O(log n)");
    RUN(MonotoneCHTMin<ll> c2);
    RUN(c2.add_line(2, 0); c2.add_line(0, 3); c2.add_line(-1, 7));
    SHOW(c2.query(10), "");
    SHOW(c2.query(-1), "query は順番を気にしない");
    SHOW(c2.query_dec(2), "");
    SHOW(c2.query_dec(-1), "");

    api("【3】 MonotoneCHTMax<T>", "T", "最大値の版（傾きは増える順に追加）。直線がないときは EMPTY", "ならし O(1)");
    RUN(MonotoneCHTMax<ll> mx);
    SHOW(mx.query(0) == MonotoneCHTMax<ll>::EMPTY, "直線がない");
    RUN(mx.add_line(-1, 5); mx.add_line(1, 0); mx.add_line(3, -10));
    SHOW(mx.query_inc(0), "max(5, 0, −10)");
    SHOW(mx.query_inc(4), "max(1, 4, 2)");
    SHOW(mx.query_inc(9), "max(−4, 9, 17)");

    api("【4】 EDPC Z（Frog 3）の形", "（使い方の例）", "dp[i] = min_j (dp[j] + (h_i − h_j)^2) + C = h_i^2 + C + min_j(−2 h_j · h_i + dp[j] + h_j^2)", "全体 O(N)");
    {
        const std::vector<ll> h = {1, 2, 3, 4, 5};
        const ll C = 6;
        const int n = (int)h.size();
        std::vector<ll> dp(n, 0);
        MonotoneCHTMin<ll> c;
        c.add_line(-2 * h[0], dp[0] + h[0] * h[0]);
        for (int i = 1; i < n; ++i) {
            dp[i] = c.query_inc(h[i]) + h[i] * h[i] + C;
            c.add_line(-2 * h[i], dp[i] + h[i] * h[i]);
        }
        note("N = 5, C = 6, h = [1, 2, 3, 4, 5] → dp = " + to_s(dp) + "（答え dp[N−1] = 20 は入力例 1 と同じ）");
    }
}

}  // namespace mcht_demo

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
    mcht_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mcht_test::benchmark();    // 速度計測
    mcht_demo::demo();         // 使い方デモ
    return 0;
}
