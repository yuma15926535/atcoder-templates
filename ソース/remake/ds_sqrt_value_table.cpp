// =============================================================================
//  【制限付き最速版】平方分割の値の表（値の範囲が小さい多重集合：追加・削除 O(1)、k 番目・mex・個数 O(√V)。Mo と組み合わせる用）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : 値は 0〜V−1 の整数だけ（大きい値は座標圧縮しておく）。質問は O(√V) と遅くなる（追加・削除だけが速い）
//    できること           : 多重集合への追加・削除（O(1)）と、k 番目・mex・値の範囲の個数・種類数・最小・最大・前後の値（O(√V)）
//    一般版との比較       : OrderedBITTemplate.txt・MexBITTemplate.txt は追加・削除・質問すべて O(log V)
//                           → Mo（追加・削除が O(N√Q) 回、質問が Q 回）と組み合わせると全体 O(N√Q log V) → O(N√Q + Q√V)
//    制限に合わないとき   : 質問が追加・削除と同じくらい多いなら OrderedBITTemplate.txt（どれも O(log V)）
//
//  しくみ: 値を √V 個ずつのブロックに分け、値ごとの個数・ブロックごとの個数・ブロックごとの種類数を持つ。
//          追加・削除は 3 か所を書き換えるだけ。質問はブロックを飛ばしながら数え、最後のブロックの中だけ 1 つずつ見る
//    例: Mo で「区間の k 番目（中央値）」「区間の mex」「区間で値が [a, b) のものの個数」を答える
//
//  機能一覧（値は 0〜V−1。値の範囲は半開区間 [a, b)。閉区間は _closed 版）                          計算量
//    SqrtValueTable st(V, B = 0)          空の多重集合（B = ブロックの大きさ。0 なら √V くらいの 2 のべき）     O(V)
//      st.add(x) / st.erase(x)            x を 1 個入れる / 1 個消す（無ければ false）                        O(1)
//      st.count(x) / st.contains(x)       x の個数 / 入っているか                                            O(1)
//      st.size() / st.distinct()          全部の個数 / 種類数                                               O(1)
//      st.count_lt(x) / count_le / count_gt / count_ge   x 未満 / 以下 / より大きい / 以上の個数               O(√V)
//      st.count_range(a, b) / count_range_closed(a, b)   値が [a, b) / [a, b] の個数                           O(√V)
//      st.distinct_range(a, b) / distinct_range_closed   値が [a, b) / [a, b] の種類数                         O(√V)
//      st.kth(k)                          小さい方から k 番目（0 始まり・同じ値も数える）。k ≥ size なら −1      O(√V)
//      st.mex()                           入っていない最小の値（全部あれば V）                                O(√V)
//      st.min() / st.max()                最小 / 最大（空なら −1）                                            O(√V)
//      st.next(x) / st.prev(x)            x 以上で最小 / x 以下で最大（無ければ −1）                           O(√V)
//      st.clear() / st.value_range()      全部消す / V                                                       O(V) / O(1)
//
//  速さの目安: N = Q = 2×10^5・値 < 2×10^5 で Mo + 区間の中央値 約 0.28 秒（一般版の Mo + BIT の k 番目なら 約 1.1 秒）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    svt_test::random_test()  … 値ごとの個数をそのまま数える愚直解との突き合わせ（Mo の形も）
//    svt_test::benchmark()    … 速度計測（Mo + BIT との比較つき）
//    svt_demo::demo()         … 全機能の使い方デモ（日本語の説明付き。Mo の形も）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
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
namespace svt_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace svt_internal

struct SqrtValueTable {
    explicit SqrtValueTable(int V = 0, int B = 0) : V_(V) {
        if (V < 0) svt_internal::fail("SqrtValueTable", "V は 0 以上");
        if (B <= 0) B = std::max(1, (int)std::sqrt((double)std::max(V, 1)));
        while ((1 << s_) < B) ++s_;  // ブロックの大きさは 2 のべき（割り算をシフトにする）
        B_ = 1 << s_;
        nb_ = (V_ + B_ - 1) >> s_;
        cnt_.assign(V_, 0), bcnt_.assign(nb_, 0), bdist_.assign(nb_, 0);
    }
    void add(int x) {
        const int b = x >> s_;
        if (cnt_[x]++ == 0) ++bdist_[b], ++distinct_;
        ++bcnt_[b], ++size_;
    }
    bool erase(int x) {
        if (cnt_[x] == 0) return false;
        const int b = x >> s_;
        if (--cnt_[x] == 0) --bdist_[b], --distinct_;
        --bcnt_[b], --size_;
        return true;
    }
    int count(int x) const { return (0 <= x && x < V_) ? cnt_[x] : 0; }
    bool contains(int x) const { return count(x) > 0; }
    long long size() const { return size_; }
    int distinct() const { return distinct_; }
    int value_range() const { return V_; }

    long long count_lt(long long x) const {  // x 未満の個数
        const int r = clamp(x);
        long long s = 0;
        const int bl = r >> s_;
        for (int b = 0; b < bl; ++b) s += bcnt_[b];
        for (int v = bl << s_; v < r; ++v) s += cnt_[v];
        return s;
    }
    long long count_le(long long x) const { return count_lt(x + 1); }
    long long count_gt(long long x) const { return size_ - count_le(x); }
    long long count_ge(long long x) const { return size_ - count_lt(x); }
    long long count_range(long long a, long long b) const { return a >= b ? 0 : count_lt(b) - count_lt(a); }
    long long count_range_closed(long long a, long long b) const { return count_range(a, b + 1); }
    int distinct_range(long long a, long long b) const { return a >= b ? 0 : distinct_lt(b) - distinct_lt(a); }
    int distinct_range_closed(long long a, long long b) const { return distinct_range(a, b + 1); }

    int kth(long long k) const {  // 小さい方から k 番目（0 始まり）。無ければ −1
        if (k < 0 || k >= size_) return -1;
        int b = 0;
        while (k >= bcnt_[b]) k -= bcnt_[b++];
        int v = b << s_;
        while (k >= cnt_[v]) k -= cnt_[v++];
        return v;
    }
    int mex() const {
        int b = 0;
        while (b < nb_ && bdist_[b] == block_len(b)) ++b;
        if (b == nb_) return V_;
        int v = b << s_;
        while (cnt_[v] > 0) ++v;
        return v;
    }
    int min() const { return next(0); }
    int max() const { return prev((long long)V_ - 1); }
    int next(long long x) const {  // x 以上で最小（無ければ −1）
        if (x < 0) x = 0;
        if (x >= V_) return -1;
        int v = (int)x;
        const int end = std::min(V_, ((v >> s_) + 1) << s_);
        for (; v < end; ++v)
            if (cnt_[v]) return v;
        int b = (int)(x >> s_) + 1;
        while (b < nb_ && bcnt_[b] == 0) ++b;
        if (b == nb_) return -1;
        v = b << s_;
        while (cnt_[v] == 0) ++v;
        return v;
    }
    int prev(long long x) const {  // x 以下で最大（無ければ −1）
        if (x >= V_) x = (long long)V_ - 1;
        if (x < 0) return -1;
        int v = (int)x;
        const int beg = (v >> s_) << s_;
        for (; v >= beg; --v)
            if (cnt_[v]) return v;
        int b = (beg >> s_) - 1;
        while (b >= 0 && bcnt_[b] == 0) --b;
        if (b < 0) return -1;
        v = std::min(V_, (b + 1) << s_) - 1;
        while (cnt_[v] == 0) --v;
        return v;
    }
    void clear() {
        std::fill(cnt_.begin(), cnt_.end(), 0), std::fill(bcnt_.begin(), bcnt_.end(), 0), std::fill(bdist_.begin(), bdist_.end(), 0);
        size_ = 0, distinct_ = 0;
    }

   private:
    int V_, s_ = 0, B_ = 1, nb_ = 0, distinct_ = 0;
    long long size_ = 0;
    std::vector<int> cnt_, bcnt_, bdist_;  // 値ごとの個数 / ブロックごとの個数 / ブロックごとの種類数
    int clamp(long long x) const { return x < 0 ? 0 : (x > V_ ? V_ : (int)x); }
    int block_len(int b) const { return std::min(B_, V_ - (b << s_)); }
    int distinct_lt(long long x) const {
        const int r = clamp(x);
        int s = 0;
        const int bl = r >> s_;
        for (int b = 0; b < bl; ++b) s += bdist_[b];
        for (int v = bl << s_; v < r; ++v) s += cnt_[v] > 0;
        return s;
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace svt_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define SVT_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_sqrt_value_table] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// Mo の並び（テスト・速度計測用。本番は MoSortTemplate.txt の run を使ってよい）
static std::vector<int> mo_order(int n, const std::vector<std::pair<int, int>>& qs) {
    const int q = (int)qs.size();
    const int B = std::max(1, (int)(n / std::sqrt((double)std::max(q, 1))));
    std::vector<int> ord(q);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        const int ba = qs[a].first / B, bb = qs[b].first / B;
        if (ba != bb) return ba < bb;
        return (ba & 1) ? qs[a].second > qs[b].second : qs[a].second < qs[b].second;
    });
    return ord;
}

static void random_test() {
    std::mt19937_64 rng(20261101);
    long long cases = 0;
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int V = (int)rnd(rng, 0, 40), B = (int)rnd(rng, 0, 9);
        SqrtValueTable st(V, B);
        std::vector<int> cnt(V, 0);
        for (int t = 0; t < 120; ++t) {
            if (V > 0) {
                const int x = (int)rnd(rng, 0, V - 1);
                if (rnd(rng, 0, 2)) {
                    st.add(x), ++cnt[x];
                } else {
                    SVT_CHECK(st.erase(x) == (cnt[x] > 0));
                    if (cnt[x]) --cnt[x];
                }
            }
            long long sz = 0;
            int ds = 0;
            for (int v = 0; v < V; ++v) sz += cnt[v], ds += cnt[v] > 0;
            SVT_CHECK(st.size() == sz && st.distinct() == ds && st.value_range() == V);
            const ll a = rnd(rng, -3, V + 3), b = rnd(rng, -3, V + 3);
            auto cl = [&](ll lo, ll hi) {  // [lo, hi) の個数・種類数
                long long c = 0;
                int d = 0;
                for (ll v = std::max(lo, 0LL); v < std::min(hi, (ll)V); ++v) c += cnt[v], d += cnt[v] > 0;
                return std::make_pair(c, d);
            };
            SVT_CHECK(st.count(a) == ((0 <= a && a < V) ? cnt[a] : 0) && st.contains(a) == (st.count(a) > 0));
            SVT_CHECK(st.count_lt(a) == cl(-5, a).first && st.count_le(a) == cl(-5, a + 1).first);
            SVT_CHECK(st.count_gt(a) == cl(a + 1, V + 5).first && st.count_ge(a) == cl(a, V + 5).first);
            SVT_CHECK(st.count_range(a, b) == cl(a, b).first && st.count_range_closed(a, b) == cl(a, b + 1).first);
            SVT_CHECK(st.distinct_range(a, b) == cl(a, b).second && st.distinct_range_closed(a, b) == cl(a, b + 1).second);
            // kth
            const ll k = rnd(rng, -1, sz + 1);
            int exp = -1;
            if (0 <= k && k < sz) {
                ll c = 0;
                for (int v = 0; v < V; ++v) {
                    c += cnt[v];
                    if (c > k) {
                        exp = v;
                        break;
                    }
                }
            }
            SVT_CHECK(st.kth(k) == exp);
            int mex = 0;
            while (mex < V && cnt[mex] > 0) ++mex;
            SVT_CHECK(st.mex() == mex);
            int nx = -1, pv = -1;
            for (ll v = std::max(a, 0LL); v < V; ++v)
                if (cnt[v]) {
                    nx = (int)v;
                    break;
                }
            for (ll v = std::min(a, (ll)V - 1); v >= 0; --v)
                if (cnt[v]) {
                    pv = (int)v;
                    break;
                }
            SVT_CHECK(st.next(a) == nx && st.prev(a) == pv);
            SVT_CHECK(st.min() == st.next(0) && st.max() == st.prev(V - 1));
        }
        st.clear();
        SVT_CHECK(st.size() == 0 && st.mex() == 0 && st.min() == -1);
        ++cases;
    }
    // Mo の形：区間の mex と k 番目
    for (int it = 0; it < 30 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 60), V = (int)rnd(rng, 1, 30), q = (int)rnd(rng, 1, 60);
        std::vector<int> a(n);
        for (auto& x : a) x = (int)rnd(rng, 0, V - 1);
        std::vector<std::pair<int, int>> qs(q);
        for (auto& p : qs) {
            int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
            if (l > r) std::swap(l, r);
            p = {l, r};
        }
        SqrtValueTable st(V);
        int L = 0, R = 0;
        for (int id : mo_order(n, qs)) {
            while (R < qs[id].second) st.add(a[R++]);
            while (L > qs[id].first) st.add(a[--L]);
            while (R > qs[id].second) st.erase(a[--R]);
            while (L < qs[id].first) st.erase(a[L++]);
            std::vector<int> b(a.begin() + L, a.begin() + R);
            std::sort(b.begin(), b.end());
            int mex = 0;
            for (int x : b)
                if (x == mex) ++mex;
            SVT_CHECK(st.mex() == mex);
            SVT_CHECK(st.kth((R - L) / 2) == (R > L ? b[(R - L) / 2] : -1));
        }
        ++cases;
    }
    std::cout << "[ds_sqrt_value_table] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_sqrt_value_table] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 30000 : 200000, Q = N, V = N;
    std::vector<int> a(N);
    for (auto& x : a) x = (int)rnd(rng, 0, V - 1);
    std::vector<std::pair<int, int>> qs(Q);
    for (auto& p : qs) {
        int l = (int)rnd(rng, 0, N - 1), r = (int)rnd(rng, 1, N);
        if (l >= r) std::swap(l, r), ++r;
        if (r > N) r = N;
        if (l >= r) l = r - 1;
        p = {l, r};
    }
    const std::vector<int> ord = mo_order(N, qs);
    ll chk1 = 0, chk2 = 0;
    auto t0 = std::chrono::steady_clock::now();
    {
        SqrtValueTable st(V);
        int L = 0, R = 0;
        for (int id : ord) {
            while (R < qs[id].second) st.add(a[R++]);
            while (L > qs[id].first) st.add(a[--L]);
            while (R > qs[id].second) st.erase(a[--R]);
            while (L < qs[id].first) st.erase(a[L++]);
            chk1 += st.kth((R - L) / 2);
        }
    }
    const ll ms1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    {
        // 一般版の形：BIT（追加・削除・k 番目すべて O(log V)）
        std::vector<int> bit(V + 1, 0);
        int LOG = 1;
        while ((1 << LOG) <= V) ++LOG;
        auto add = [&](int x, int d) {
            for (++x; x <= V; x += x & -x) bit[x] += d;
        };
        auto kth = [&](int k) {  // 0 始まり
            int pos = 0;
            for (int pw = 1 << (LOG - 1); pw; pw >>= 1)
                if (pos + pw <= V && bit[pos + pw] <= k) pos += pw, k -= bit[pos];
            return pos;
        };
        int L = 0, R = 0;
        for (int id : ord) {
            while (R < qs[id].second) add(a[R++], 1);
            while (L > qs[id].first) add(a[--L], 1);
            while (R > qs[id].second) add(a[--R], -1);
            while (L < qs[id].first) add(a[L++], -1);
            chk2 += kth((R - L) / 2);
        }
    }
    const ll ms2 = elapsed_ms(t0);
    SVT_CHECK(chk1 == chk2);
    std::cout << "[ds_sqrt_value_table] benchmark: N = Q = " << N << " の Mo + 区間の中央値  SqrtValueTable " << ms1 << " ms / 一般版（BIT）" << ms2
              << " ms  (check " << chk1 % 1000 << ")\n";
}

}  // namespace svt_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace svt_demo {

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
#define SHOW(expr, comment) svt_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    svt_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】平方分割の値の表（ds_sqrt_value_table） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 値は 0〜V−1 の整数だけ・質問は O(√V)（追加・削除だけが O(1)）\n"
              << "■ できること: 追加・削除 O(1)、k 番目・mex・範囲の個数・種類数・最小・最大・前後 O(√V)（Mo と組み合わせる）\n"
              << "■ API 一覧\n"
              << "  SqrtValueTable st(V, B = 0): add / erase / count / contains / size / distinct / count_lt / count_le / count_gt / count_ge /\n"
              << "    count_range / count_range_closed / distinct_range / distinct_range_closed / kth / mex / min / max / next / prev / clear / value_range\n";

    api("【1】 SqrtValueTable st(V)  add(x) / erase(x)", "void / bool", "x を 1 個入れる / 1 個消す（無ければ false）", "O(1)");
    RUN(SqrtValueTable st(20));
    RUN(st.add(3); st.add(3); st.add(0); st.add(1); st.add(7); st.add(12));
    SHOW(st.erase(5), "入っていない → false");
    SHOW(st.size(), "");
    SHOW(st.distinct(), "種類数（3 が 2 個）");
    SHOW(st.count(3), "");
    SHOW(st.contains(4), "");

    api("【2】 count_lt / count_le / count_gt / count_ge / count_range / distinct_range", "long long / int", "値の大小・範囲の個数 / 種類数", "O(√V)");
    SHOW(st.count_lt(3), "3 未満 = {0, 1}");
    SHOW(st.count_le(3), "");
    SHOW(st.count_gt(3), "");
    SHOW(st.count_ge(3), "");
    SHOW(st.count_range(1, 8), "値が [1, 8) = {1, 3, 3, 7}");
    SHOW(st.count_range_closed(1, 7), "");
    SHOW(st.distinct_range(1, 8), "種類は {1, 3, 7}");
    SHOW(st.distinct_range_closed(0, 3), "");

    api("【3】 kth(k) / mex() / min() / max() / next(x) / prev(x)", "int", "k 番目（0 始まり）/ 入っていない最小 / 最小 / 最大 / x 以上で最小 / x 以下で最大", "O(√V)");
    SHOW(st.kth(2), "小さい順 {0, 1, 3, 3, 7, 12} の 2 番目");
    SHOW(st.kth(6), "個数以上 → −1");
    SHOW(st.mex(), "");
    SHOW(st.min(), "");
    SHOW(st.max(), "");
    SHOW(st.next(4), "");
    SHOW(st.prev(11), "");
    SHOW(st.next(13), "無い → −1");
    RUN(st.clear());
    SHOW(st.size(), "clear で空に");

    api("【4】 Mo と組み合わせる", "（使い方の例）", "区間を 1 つずつ伸ばし縮めして add / erase、質問のたびに kth・mex", "全体 O(N√Q + Q√V)");
    {
        const std::vector<int> a = {2, 0, 3, 1, 0, 4, 2, 5};
        const std::vector<std::pair<int, int>> qs = {{0, 4}, {2, 7}, {1, 3}, {4, 8}};
        SqrtValueTable t(6);
        std::vector<int> mex(qs.size()), med(qs.size());
        int L = 0, R = 0;
        for (int id : svt_test::mo_order((int)a.size(), qs)) {
            while (R < qs[id].second) t.add(a[R++]);
            while (L > qs[id].first) t.add(a[--L]);
            while (R > qs[id].second) t.erase(a[--R]);
            while (L < qs[id].first) t.erase(a[L++]);
            mex[id] = t.mex(), med[id] = t.kth((R - L - 1) / 2);
        }
        note("a = [2, 0, 3, 1, 0, 4, 2, 5]、区間 [0,4) [2,7) [1,3) [4,8)");
        note("区間の mex → " + to_s(mex) + "、区間の中央値（小さい方）→ " + to_s(med));
        note("本番では MoSortTemplate.txt の run(add, erase, out) に add / erase と、out で st.mex() などを記録する関数を渡す");
    }
}

}  // namespace svt_demo

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
    svt_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    svt_test::benchmark();    // 速度計測
    svt_demo::demo();         // 使い方デモ
    return 0;
}
