// =============================================================================
//  平方分割（√N ずつのブロックに分けて、端は 1 つずつ・中はブロックごとにまとめて処理する型）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  セグメント木に乗らない操作（「区間に足す」と「区間で x 未満の数」を両方、など）を O(√N) か O(√N log N) で。
//  区間 [l, r) を「左の端のブロックの一部」「まるごと入るブロック」「右の端のブロックの一部」に分け、
//  一部は要素を 1 つずつ、まるごとはブロックの情報（遅延・並べ替えた写し など）で処理する
//    例: 区間に足す・区間で x 以上の個数、区間に足す・区間で x 以下の最大、ブロックごとに何かを数え直す
//
//  機能一覧（区間は 0-indexed・半開区間 [l, r)）                                               計算量
//    SqrtBlocks sb(n, B = 0)            n 個を大きさ B のブロックに（B = 0 なら √n くらい）             O(1)
//      sb.block(i) / sb.lo(b) / sb.hi(b) / sb.count()    i のブロック / ブロック b の範囲 [lo, hi) / ブロックの数
//      sb.for_range(l, r, part, whole)  一部 → part(b, l', r')（ブロック b の [l', r') だけ）、まるごと → whole(b)   O(n / B + B)
//    RangeAddCountLess<T> rc(a, B = 0)  区間に足す + 区間で値の大小を数える（ブロックごとに並べた写しと遅延。B = 0 なら √n）
//      rc.add(l, r, x)                  [l, r) の全部に x を足す（端のブロックはマージで並べ直す）         O(n / B + B)
//      rc.count_lt / count_le / count_gt / count_ge(l, r, x)   [l, r) で x 未満 / 以下 / より大きい / 以上の数   O((n / B) log B + B)
//      rc.max_lt(l, r, x) / rc.min_ge(l, r, x)   [l, r) で x 未満の最大 / x 以上の最小（無ければ nullopt）
//      rc.get(i) / rc.set(i, x) / rc.size()                                                        O(1) / O(B)
//
//  注意: B は「一部の処理の重さ」と「ブロックの数」のつり合いで決める。並べ替えを使うなら √(n log n) くらいが速いこともある
//  速さの目安: RangeAddCountLess（n = 10^5）で add と count を 5×10^4 回ずつ 約 0.55 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    sqrtdec_test::random_test()  … 配列をそのまま持つ愚直解との突き合わせ
//    sqrtdec_test::benchmark()    … 速度計測
//    sqrtdec_demo::demo()         … 全機能の使い方デモ（日本語の説明付き。自分で書き換える例も）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
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
namespace sqrtdec_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace sqrtdec_internal

struct SqrtBlocks {
    int n, B, nb;
    explicit SqrtBlocks(int n_ = 0, int B_ = 0) : n(n_) {
        if (n_ < 0) sqrtdec_internal::fail("SqrtBlocks", "n は 0 以上");
        B = B_ > 0 ? B_ : std::max(1, (int)std::sqrt((double)std::max(1, n_)));
        nb = (n + B - 1) / B;
    }
    int block(int i) const { return i / B; }
    int lo(int b) const { return b * B; }
    int hi(int b) const { return std::min(n, (b + 1) * B); }
    int count() const { return nb; }
    template <class Part, class Whole>
    void for_range(int l, int r, Part part, Whole whole) const {
        if (l < 0 || r > n || l > r) sqrtdec_internal::fail("SqrtBlocks::for_range", "0 ≤ l ≤ r ≤ n の半開区間 [l, r)");
        if (l == r) return;
        const int bl = block(l), br = block(r - 1);
        if (bl == br) {
            if (l == lo(bl) && r == hi(bl)) whole(bl);
            else part(bl, l, r);
            return;
        }
        if (l == lo(bl)) whole(bl);
        else part(bl, l, hi(bl));
        for (int b = bl + 1; b < br; ++b) whole(b);
        if (r == hi(br)) whole(br);
        else part(br, lo(br), r);
    }
};

template <class T>
struct RangeAddCountLess {
    explicit RangeAddCountLess(const std::vector<T>& a, int B = 0)
        : sb_((int)a.size(), B > 0 ? B : default_block((int)a.size())), a_(a), lazy_(sb_.count(), T(0)), sorted_(sb_.count()), key_(sb_.count()) {
        for (int b = 0; b < sb_.count(); ++b) {
            for (int i = sb_.lo(b); i < sb_.hi(b); ++i) sorted_[b].push_back({a_[i], i});
            std::sort(sorted_[b].begin(), sorted_[b].end());
            refresh_key(b);
        }
    }
    int size() const { return sb_.n; }
    T get(int i) const {
        check(i, "RangeAddCountLess::get");
        return a_[i] + lazy_[sb_.block(i)];
    }
    void set(int i, T x) {
        check(i, "RangeAddCountLess::set");
        const int b = sb_.block(i);
        auto& s = sorted_[b];
        s.erase(std::find(s.begin(), s.end(), std::make_pair(a_[i], i)));
        a_[i] = x - lazy_[b];
        s.insert(std::lower_bound(s.begin(), s.end(), std::make_pair(a_[i], i)), {a_[i], i});
        refresh_key(b);
    }
    void add(int l, int r, T x) {
        sb_.for_range(
            l, r,
            [&](int b, int pl, int pr) {  // 並べた写しを「足すもの」「足さないもの」に分けてマージ（O(B)）
                auto& s = sorted_[b];
                in_.clear(), out_.clear();
                for (const auto& p : s)
                    if (pl <= p.second && p.second < pr) in_.push_back({p.first + x, p.second});
                    else out_.push_back(p);
                for (int i = pl; i < pr; ++i) a_[i] += x;
                std::merge(in_.begin(), in_.end(), out_.begin(), out_.end(), s.begin());
                refresh_key(b);
            },
            [&](int b) { lazy_[b] += x; });
    }
    int count_lt(int l, int r, T x) const { return count_if_(l, r, x, false); }
    int count_le(int l, int r, T x) const { return count_if_(l, r, x, true); }
    int count_ge(int l, int r, T x) const { return (r - l) - count_lt(l, r, x); }
    int count_gt(int l, int r, T x) const { return (r - l) - count_le(l, r, x); }
    std::optional<T> max_lt(int l, int r, T x) const {
        std::optional<T> best;
        auto upd = [&](T v) {
            if (v < x && (!best || *best < v)) best = v;
        };
        sb_.for_range(
            l, r, [&](int b, int pl, int pr) { for (int i = pl; i < pr; ++i) upd(a_[i] + lazy_[b]); },
            [&](int b) {
                const auto& k = key_[b];
                auto it = std::lower_bound(k.begin(), k.end(), x - lazy_[b]);
                if (it != k.begin()) upd(*std::prev(it) + lazy_[b]);
            });
        return best;
    }
    std::optional<T> min_ge(int l, int r, T x) const {
        std::optional<T> best;
        auto upd = [&](T v) {
            if (!(v < x) && (!best || v < *best)) best = v;
        };
        sb_.for_range(
            l, r, [&](int b, int pl, int pr) { for (int i = pl; i < pr; ++i) upd(a_[i] + lazy_[b]); },
            [&](int b) {
                const auto& k = key_[b];
                auto it = std::lower_bound(k.begin(), k.end(), x - lazy_[b]);
                if (it != k.end()) upd(*it + lazy_[b]);
            });
        return best;
    }

   private:
    SqrtBlocks sb_;
    std::vector<T> a_, lazy_;                        // 本当の値 = a_[i] + lazy_[ブロック]
    std::vector<std::vector<std::pair<T, int>>> sorted_;  // ブロックごとの (a_[i], i) を並べたもの
    std::vector<std::vector<T>> key_;                // その値だけ（二分探索用）
    std::vector<std::pair<T, int>> in_, out_;        // add の作業用
    static int default_block(int n) { return std::max(1, (int)std::sqrt((double)std::max(1, n))); }
    void refresh_key(int b) {
        key_[b].resize(sorted_[b].size());
        for (size_t k = 0; k < sorted_[b].size(); ++k) key_[b][k] = sorted_[b][k].first;
    }
    int count_if_(int l, int r, T x, bool le) const {
        int c = 0;
        sb_.for_range(
            l, r,
            [&](int b, int pl, int pr) {
                for (int i = pl; i < pr; ++i) {
                    const T v = a_[i] + lazy_[b];
                    c += le ? !(x < v) : (v < x);
                }
            },
            [&](int b) {
                const auto& k = key_[b];
                const T y = x - lazy_[b];
                c += (int)((le ? std::upper_bound(k.begin(), k.end(), y) : std::lower_bound(k.begin(), k.end(), y)) - k.begin());
            });
        return c;
    }
    void check(int i, const char* f) const {
        if (i < 0 || i >= sb_.n) sqrtdec_internal::fail(f, "i が範囲外");
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace sqrtdec_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define SQ_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_sqrt_decomposition] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 40), B = (int)rnd(rng, 0, 7);
        // for_range が区間をちょうど 1 回ずつ覆うか
        {
            const SqrtBlocks sb(n, B);
            for (int q = 0; q < 10; ++q) {
                int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
                if (l > r) std::swap(l, r);
                std::vector<int> cover(n, 0);
                sb.for_range(
                    l, r,
                    [&](int b, int pl, int pr) {
                        SQ_CHECK(sb.lo(b) <= pl && pl < pr && pr <= sb.hi(b) && !(pl == sb.lo(b) && pr == sb.hi(b)));
                        for (int i = pl; i < pr; ++i) cover[i]++;
                    },
                    [&](int b) {
                        for (int i = sb.lo(b); i < sb.hi(b); ++i) cover[i]++;
                    });
                for (int i = 0; i < n; ++i) SQ_CHECK(cover[i] == (l <= i && i < r ? 1 : 0));
                ++cases;
            }
        }
        std::vector<ll> a(n);
        for (auto& x : a) x = rnd(rng, -10, 10);
        RangeAddCountLess<ll> rc(a, B);
        for (int t = 0; t < 80; ++t) {
            int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
            if (l > r) std::swap(l, r);
            const ll x = rnd(rng, -15, 15);
            const int k = (int)rnd(rng, 0, 3);
            if (k == 0) {
                rc.add(l, r, x);
                for (int i = l; i < r; ++i) a[i] += x;
            } else if (k == 1 && n > 0) {
                const int i = (int)rnd(rng, 0, n - 1);
                rc.set(i, x), a[i] = x;
            } else {
                int lt = 0, le = 0;
                std::optional<ll> mx, mn;
                for (int i = l; i < r; ++i) {
                    lt += a[i] < x, le += a[i] <= x;
                    if (a[i] < x && (!mx || *mx < a[i])) mx = a[i];
                    if (a[i] >= x && (!mn || a[i] < *mn)) mn = a[i];
                }
                SQ_CHECK(rc.count_lt(l, r, x) == lt && rc.count_le(l, r, x) == le);
                SQ_CHECK(rc.count_ge(l, r, x) == (r - l) - lt && rc.count_gt(l, r, x) == (r - l) - le);
                SQ_CHECK(rc.max_lt(l, r, x) == mx && rc.min_ge(l, r, x) == mn);
            }
            for (int i = 0; i < n; i += 5) SQ_CHECK(rc.get(i) == a[i]);
            ++cases;
        }
    }
    std::cout << "[ds_sqrt_decomposition] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_sqrt_decomposition] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 10000 : 100000;
    std::vector<ll> a(N);
    for (auto& x : a) x = rnd(rng, 0, 1000000000);
    const auto t0 = std::chrono::steady_clock::now();
    RangeAddCountLess<ll> rc(a);
    ll chk = 0;
    for (int i = 0; i < N; ++i) {
        int l = (int)rnd(rng, 0, N), r = (int)rnd(rng, 0, N);
        if (l > r) std::swap(l, r);
        if (i & 1) rc.add(l, r, rnd(rng, -1000, 1000));
        else chk += rc.count_ge(l, r, rnd(rng, 0, 1000000000));
    }
    std::cout << "[ds_sqrt_decomposition] benchmark: RangeAddCountLess（n = " << N << "）で add と count_ge を " << N / 2 << " 回ずつ " << elapsed_ms(t0)
              << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace sqrtdec_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace sqrtdec_demo {

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
#define SHOW(expr, comment) sqrtdec_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    sqrtdec_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("平方分割（ds_sqrt_decomposition） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  SqrtBlocks sb(n, B): block / lo / hi / count / for_range(l, r, part, whole)\n"
              << "  RangeAddCountLess<T> rc(a, B): add / count_lt / count_le / count_gt / count_ge / max_lt / min_ge / get / set / size\n";

    api("【1】 RangeAddCountLess<T> rc(a) / add(l, r, x) / count_ge(l, r, x)", "— / void / int", "区間に足す・区間で x 以上の数（セグ木に乗らない組み合わせ）",
        "O(n log n) / O(√n log n) / O(√n log n)");
    RUN(RangeAddCountLess<long long> rc({5, 1, 4, 1, 5, 9, 2, 6}));
    SHOW(rc.count_ge(0, 8, 5), "5, 5, 9, 6");
    RUN(rc.add(1, 4, 10));  // 1, 4, 1 → 11, 14, 11
    SHOW(rc.count_ge(0, 8, 5), "");
    SHOW(rc.count_lt(0, 4, 12), "5, 11, 11");
    SHOW(rc.max_lt(0, 8, 10), "10 未満の最大");
    SHOW(rc.min_ge(0, 8, 12), "12 以上の最小");

    api("【2】 SqrtBlocks sb(n, B) / for_range(l, r, part, whole)", "—", "自分で平方分割を書くときの枠：区間を「端の一部」と「まるごとのブロック」に分ける",
        "O(n / B + B)");
    {
        code("// 例: 区間に足す・区間の和（自分でブロックごとの和と遅延を持つ）");
        const int n = 10;
        SqrtBlocks sb(n, 3);
        std::vector<long long> a(n, 1), sum(sb.count()), lazy(sb.count(), 0);
        for (int b = 0; b < sb.count(); ++b) sum[b] = sb.hi(b) - sb.lo(b);
        auto add = [&](int l, int r, long long x) {
            sb.for_range(
                l, r, [&](int b, int pl, int pr) { for (int i = pl; i < pr; ++i) a[i] += x, sum[b] += x; },
                [&](int b) { lazy[b] += x; });
        };
        auto get_sum = [&](int l, int r) {
            long long s = 0;
            sb.for_range(
                l, r, [&](int b, int pl, int pr) { for (int i = pl; i < pr; ++i) s += a[i] + lazy[b]; },
                [&](int b) { s += sum[b] + lazy[b] * (sb.hi(b) - sb.lo(b)); });
            return s;
        };
        add(2, 8, 5);
        note("全部 1 の長さ 10 で [2, 8) に 5 を足す → 和 [0, 10) = " + to_s(get_sum(0, 10)) + "、[1, 3) = " + to_s(get_sum(1, 3)));
    }
}

}  // namespace sqrtdec_demo

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
    sqrtdec_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    sqrtdec_test::benchmark();    // 速度計測
    sqrtdec_demo::demo();         // 使い方デモ
    return 0;
}
