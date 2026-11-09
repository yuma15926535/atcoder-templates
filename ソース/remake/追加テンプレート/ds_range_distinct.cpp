// =============================================================================
//  区間の種類数（区間 [l, r) にある値の種類数・k 回以上 / ちょうど k 回出る値の数・種類ごとの和）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  オフライン（質問を先に全部読める）: 右端の順に見ていき、「各値の最後の出現の位置だけ 1」を BIT で持つ。
//    [l, r) の種類数 = 位置 l〜r−1 の 1 の数。k 回以上なら「後ろから k 番目の出現」に 1 を置く。O((N + Q) log N)
//  オンライン（質問に答えてから次の質問が来る）: prev[i]（同じ値の 1 つ前の位置）を使い、
//    種類数 = #{ i ∈ [l, r) : prev[i] < l } を Wavelet Matrix で数える。1 回 O(log N)
//    例: ABC174 F（Range Set Query）… range_distinct_offline(c, queries)
//
//  機能一覧（区間は 0-indexed・半開区間 [l, r)。閉区間 [l, r] の質問は {l, r + 1} にして渡す）    計算量
//   オフライン（qs = vector<pair<int, int>> の質問 {l, r}。答えは質問の順）
//    range_distinct_offline(a, qs)              [l, r) にある値の種類数                         O((N + Q) log N)
//    range_count_at_least_offline(a, qs, k)     [l, r) に k 回以上出てくる値の種類数（k ≥ 1）   O((N k + Q) log N)
//    range_count_exactly_offline(a, qs, k)      [l, r) にちょうど k 回出てくる値の種類数（k ≥ 1）
//    range_distinct_sum_offline(a, qs, weight)  [l, r) に出てくる値 x（種類ごとに 1 回）の weight(x) の和   O((N + Q) log N)
//   オンライン
//    RangeDistinct rd(a)                        前処理                                         O(N log N)
//    rd.distinct(l, r) / rd.distinct_closed(l, r)   [l, r) / [l, r] の種類数                   O(log N)
//    rd.size()                                  N
//
//  a の値は == と < で比べられる型なら何でもよい（int・long long・string・pair など。中で座標圧縮する）
//  l = r の空の区間は 0。範囲外の質問（l < 0、r > N、l > r）は止める
//  速さの目安（N = Q = 2×10^5）: range_distinct_offline 約 0.06 秒、RangeDistinct（前処理 + 質問）約 0.08 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    rangedistinct_test::random_test()  … 区間を全部数える愚直解との突き合わせ
//    rangedistinct_test::benchmark()    … 速度計測
//    rangedistinct_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace rangedistinct_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
template <class V>
std::vector<int> compress(const std::vector<V>& a) {  // 値を 0〜(種類数 − 1) の番号に（並べ替えた順に番号を振る）
    std::vector<int> ord(a.size()), id(a.size());
    for (size_t i = 0; i < a.size(); ++i) ord[i] = (int)i;
    std::sort(ord.begin(), ord.end(), [&](int x, int y) { return a[x] < a[y]; });
    int k = -1;
    for (size_t i = 0; i < ord.size(); ++i) {
        if (i == 0 || a[ord[i - 1]] < a[ord[i]]) ++k;
        id[ord[i]] = k;
    }
    return id;
}
inline void check_queries(int n, const std::vector<std::pair<int, int>>& qs, const char* f) {
    for (const auto& q : qs)
        if (q.first < 0 || q.second > n || q.first > q.second) fail(f, "質問は 0 ≤ l ≤ r ≤ N の半開区間 [l, r)");
}
// 値 x の後ろから j 番目（j = 1〜K）の出現の位置に coef[j−1] を置いたときの、[l, r) の和（W は int・long long など）
template <class W, class V, class Coef>
std::vector<W> kth_last(const std::vector<V>& a, const std::vector<std::pair<int, int>>& qs, int K, Coef coef, const char* f) {
    const int n = (int)a.size();
    check_queries(n, qs, f);
    const std::vector<int> id = compress(a);
    // 値の番号ごとに、最後の K + 1 回の出現位置を輪のように持つ（o(x, t) = t 回目（0 始まり）の出現位置）
    const size_t R = (size_t)K + 1;
    std::vector<int> ring((size_t)n * R), seen(n, 0);
    auto o = [&](int x, int t) -> int& { return ring[(size_t)x * R + (size_t)t % R]; };
    std::vector<W> bit(n + 1, W(0));
    auto add = [&](int i, W v) {
        for (++i; i <= n; i += i & -i) bit[i] += v;
    };
    auto pre = [&](int i) {  // 位置 0〜i−1 の和
        W s = W(0);
        for (; i > 0; i -= i & -i) s += bit[i];
        return s;
    };
    std::vector<int> ord(qs.size());
    for (size_t i = 0; i < qs.size(); ++i) ord[i] = (int)i;
    std::sort(ord.begin(), ord.end(), [&](int x, int y) { return qs[x].second < qs[y].second; });
    std::vector<W> ans(qs.size(), W(0));
    size_t qi = 0;
    for (int r = 0; r <= n; ++r) {
        while (qi < ord.size() && qs[ord[qi]].second == r) {
            const auto [l, rr] = qs[ord[qi]];
            ans[ord[qi]] = pre(rr) - pre(l);
            ++qi;
        }
        if (r == n) break;
        const int x = id[r];
        o(x, seen[x]) = r;
        const int m = ++seen[x];  // x の出現回数（今まで）
        for (int j = 1; j <= K && m - j >= 0; ++j) {  // m − j 回目の出現は後ろから j − 1 番目 → j 番目に
            const W now = coef(a[r], j), before = j >= 2 ? coef(a[r], j - 1) : W(0);
            if (!(now == before)) add(o(x, m - j), now - before);
        }
        if (m - 1 - K >= 0) {  // 後ろから K 番目だったもの → K + 1 番目（係数 0）
            const W before = coef(a[r], K);
            if (!(before == W(0))) add(o(x, m - 1 - K), W(0) - before);
        }
    }
    return ans;
}
}  // namespace rangedistinct_internal

// ---- オフライン ----
template <class V>
std::vector<int> range_distinct_offline(const std::vector<V>& a, const std::vector<std::pair<int, int>>& qs) {
    return rangedistinct_internal::kth_last<int>(a, qs, 1, [](const V&, int) { return 1; }, "range_distinct_offline");
}
template <class V>
std::vector<int> range_count_at_least_offline(const std::vector<V>& a, const std::vector<std::pair<int, int>>& qs, int k) {
    if (k < 1) rangedistinct_internal::fail("range_count_at_least_offline", "k は 1 以上");
    return rangedistinct_internal::kth_last<int>(a, qs, k, [k](const V&, int j) { return j == k ? 1 : 0; }, "range_count_at_least_offline");
}
template <class V>
std::vector<int> range_count_exactly_offline(const std::vector<V>& a, const std::vector<std::pair<int, int>>& qs, int k) {
    if (k < 1) rangedistinct_internal::fail("range_count_exactly_offline", "k は 1 以上");
    // k 回以上 − (k + 1) 回以上
    return rangedistinct_internal::kth_last<int>(
        a, qs, k + 1, [k](const V&, int j) { return j == k ? 1 : (j == k + 1 ? -1 : 0); }, "range_count_exactly_offline");
}
template <class V, class F>
auto range_distinct_sum_offline(const std::vector<V>& a, const std::vector<std::pair<int, int>>& qs, F weight) {
    using W = decltype(weight(a[0]));
    return rangedistinct_internal::kth_last<W>(a, qs, 1, [&](const V& x, int) { return weight(x); }, "range_distinct_sum_offline");
}

// ---- オンライン（prev[i] < l の数を Wavelet Matrix で数える）----
struct RangeDistinct {
    RangeDistinct() {}
    template <class V>
    explicit RangeDistinct(const std::vector<V>& a) : n_((int)a.size()) {
        const std::vector<int> id = rangedistinct_internal::compress(a);
        std::vector<int> last(n_, -1);
        std::vector<unsigned> cur(n_);  // prev[i] + 1（0〜N）
        for (int i = 0; i < n_; ++i) cur[i] = (unsigned)(last[id[i]] + 1), last[id[i]] = i;
        lg_ = 1;
        while ((1u << lg_) <= (unsigned)n_) ++lg_;
        const size_t blocks = (size_t)n_ / 64 + 1;
        bits_.assign((size_t)lg_ * blocks, 0);
        cnt_.assign((size_t)lg_ * (blocks + 1), 0);
        zeros_.assign(lg_, 0);
        std::vector<unsigned> nxt(n_);
        for (int lv = lg_ - 1; lv >= 0; --lv) {
            uint64_t* b = &bits_[(size_t)lv * blocks];
            int z = 0;
            for (int i = 0; i < n_; ++i)
                if (cur[i] >> lv & 1) b[i >> 6] |= 1ULL << (i & 63);
                else ++z;
            int* c = &cnt_[(size_t)lv * (blocks + 1)];
            for (size_t k = 0; k < blocks; ++k) c[k + 1] = c[k] + __builtin_popcountll(b[k]);
            zeros_[lv] = z;
            int p0 = 0, p1 = z;
            for (int i = 0; i < n_; ++i) (cur[i] >> lv & 1 ? nxt[p1++] : nxt[p0++]) = cur[i];
            std::swap(cur, nxt);
        }
    }
    int size() const { return n_; }
    int distinct(int l, int r) const {
        if (l < 0 || r > n_ || l > r) rangedistinct_internal::fail("RangeDistinct::distinct", "0 ≤ l ≤ r ≤ N の半開区間 [l, r)");
        return count_lt(l, r, (unsigned)l + 1);  // prev[i] + 1 ≤ l ⇔ prev[i] < l
    }
    int distinct_closed(int l, int r) const { return distinct(l, r + 1); }

   private:
    int n_ = 0, lg_ = 1;
    std::vector<uint64_t> bits_;
    std::vector<int> cnt_, zeros_;
    int rank1(int lv, int i) const {  // レベル lv の先頭 i 個の 1 の数
        const size_t blocks = (size_t)n_ / 64 + 1;
        const uint64_t* b = &bits_[(size_t)lv * blocks];
        const int* c = &cnt_[(size_t)lv * (blocks + 1)];
        const int k = i >> 6, s = i & 63;
        return c[k] + (s ? __builtin_popcountll(b[k] & ((1ULL << s) - 1)) : 0);
    }
    int count_lt(int l, int r, unsigned upper) const {  // 位置 [l, r) で値 < upper の数
        if (upper >= (1u << lg_)) return r - l;
        int res = 0;
        for (int lv = lg_ - 1; lv >= 0; --lv) {
            const int l1 = rank1(lv, l), r1 = rank1(lv, r), l0 = l - l1, r0 = r - r1;
            if (upper >> lv & 1) {
                res += r0 - l0;
                l = zeros_[lv] + l1, r = zeros_[lv] + r1;
            } else {
                l = l0, r = r0;
            }
        }
        return res;
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace rangedistinct_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define RD_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_range_distinct] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

template <class V>
static void check_all(const std::vector<V>& a, std::mt19937_64& rng, long long& cases) {
    const int n = (int)a.size();
    std::vector<std::pair<int, int>> qs;
    const int q = (int)rnd(rng, 0, 40);
    for (int i = 0; i < q; ++i) {
        int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
        if (l > r) std::swap(l, r);
        qs.push_back({l, r});
    }
    const RangeDistinct rd(a);
    RD_CHECK(rd.size() == n);
    const auto d1 = range_distinct_offline(a, qs);
    std::vector<std::vector<int>> atl(4), exa(4);
    for (int k = 1; k <= 3; ++k) atl[k] = range_count_at_least_offline(a, qs, k), exa[k] = range_count_exactly_offline(a, qs, k);
    const auto ws = range_distinct_sum_offline(a, qs, [](const V& x) { return (ll)x * 3 + 1; });
    for (int i = 0; i < q; ++i) {
        const auto [l, r] = qs[i];
        std::map<V, int> cnt;
        for (int j = l; j < r; ++j) cnt[a[j]]++;
        RD_CHECK(d1[i] == (int)cnt.size());
        RD_CHECK(rd.distinct(l, r) == (int)cnt.size());
        if (l < r) RD_CHECK(rd.distinct_closed(l, r - 1) == (int)cnt.size());
        for (int k = 1; k <= 3; ++k) {
            int al = 0, ex = 0;
            for (auto& [x, c] : cnt) al += c >= k, ex += c == k;
            RD_CHECK(atl[k][i] == al && exa[k][i] == ex);
        }
        ll s = 0;
        for (auto& [x, c] : cnt) s += (ll)x * 3 + 1;
        RD_CHECK(ws[i] == s);
        cases += 9;
    }
    // オンラインは全区間も
    for (int l = 0; l <= n && n <= 30; ++l) {
        std::map<V, int> cnt;
        for (int r = l; r <= n; ++r) {
            RD_CHECK(rd.distinct(l, r) == (int)cnt.size());
            if (r < n) cnt[a[r]]++;
            ++cases;
        }
    }
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 600 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, rnd(rng, 0, 3) ? 30 : 200);
        const ll V = rnd(rng, 1, 8);
        std::vector<ll> a(n);
        for (auto& x : a) x = rnd(rng, -V, V) * (rnd(rng, 0, 1) ? 1 : 1000000000000LL);
        check_all(a, rng, cases);
    }
    // 長さが 2 の冪の前後（Wavelet Matrix の段の数が変わるところ）
    for (int n : {0, 1, 2, 3, 63, 64, 65, 127, 128, 129, 255, 256, 257}) {
        std::vector<int> a(n);
        for (auto& x : a) x = (int)rnd(rng, 0, n / 3 + 1);
        check_all(a, rng, cases);
    }
    // 文字列の値
    {
        std::vector<std::string> a = {"ab", "c", "ab", "d", "c", "c"};
        const auto r = range_distinct_offline(a, {{0, 6}, {1, 5}, {2, 3}, {3, 3}});
        RD_CHECK(r == std::vector<int>({3, 3, 1, 0}));
        RD_CHECK(RangeDistinct(a).distinct(0, 6) == 3);
        cases += 2;
    }
    std::cout << "[ds_range_distinct] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_range_distinct] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000, Q = N;
    std::vector<int> a(N);
    for (auto& x : a) x = (int)rnd(rng, 1, N);
    std::vector<std::pair<int, int>> qs(Q);
    for (auto& [l, r] : qs) {
        l = (int)rnd(rng, 0, N), r = (int)rnd(rng, 0, N);
        if (l > r) std::swap(l, r);
    }
    ll chk = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int x : range_distinct_offline(a, qs)) chk += x;
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (int x : range_count_exactly_offline(a, qs, 2)) chk += x;
    const ll t2 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    {
        const RangeDistinct rd(a);
        for (auto [l, r] : qs) chk += rd.distinct(l, r);
    }
    const ll t3 = elapsed_ms(t0);
    std::cout << "[ds_range_distinct] benchmark (N = Q = " << N << "): range_distinct_offline " << t1 << " ms / range_count_exactly_offline(k = 2) "
              << t2 << " ms / RangeDistinct（前処理 + 質問） " << t3 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace rangedistinct_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace rangedistinct_demo {

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
#define SHOW(expr, comment) rangedistinct_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    rangedistinct_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("区間の種類数（ds_range_distinct） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  range_distinct_offline / range_count_at_least_offline / range_count_exactly_offline / range_distinct_sum_offline\n"
              << "  RangeDistinct rd(a): distinct / distinct_closed / size\n";

    RUN(const std::vector<int> c = {1, 2, 1, 3, 2, 2, 5});
    RUN(const std::vector<std::pair<int, int>> qs = {{0, 3}, {1, 4}, {2, 3}, {0, 7}, {4, 6}, {3, 3}});

    api("【1】 range_distinct_offline(a, qs)", "vector<int>（質問の順）", "[l, r) にある値の種類数（ABC174 F）。閉区間 [l, r] は {l, r + 1} にして渡す",
        "O((N + Q) log N)");
    SHOW(range_distinct_offline(c, qs), "{1,2,1}→2, {2,1,3}→3, {1}→1, 全体→4, {2,2}→1, 空→0");

    api("【2】 range_count_at_least_offline(a, qs, k) / range_count_exactly_offline(a, qs, k)", "vector<int>",
        "[l, r) に k 回以上 / ちょうど k 回出てくる値の種類数", "O((N k + Q) log N)");
    SHOW(range_count_at_least_offline(c, qs, 2), "2 回以上出る値の数");
    SHOW(range_count_exactly_offline(c, qs, 1), "1 回だけ出る値の数");
    SHOW(range_count_exactly_offline(c, qs, 3), "全体では 2 が 3 回");

    api("【3】 range_distinct_sum_offline(a, qs, weight)", "vector<W>（W は weight の返り値の型）",
        "[l, r) に出てくる値 x（種類ごとに 1 回）の weight(x) の和", "O((N + Q) log N)");
    SHOW(range_distinct_sum_offline(c, qs, [](int x) { return (long long)x; }), "種類ごとの値の和：全体は 1 + 2 + 3 + 5");

    api("【4】 RangeDistinct rd(a) / rd.distinct(l, r) / rd.distinct_closed(l, r)", "— / int / int",
        "オンライン（質問を先読みできないとき）。[l, r) / [l, r] の種類数", "O(N log N) / O(log N)");
    RUN(const RangeDistinct rd(c));
    SHOW(rd.distinct(0, 3), "");
    SHOW(rd.distinct(0, 7), "");
    SHOW(rd.distinct_closed(4, 5), "[4, 5] = {2, 2}");

    api("【5】 値は文字列・pair などでもよい", "—", "中で座標圧縮するので == と < が使える型なら何でも", "—");
    RUN(const std::vector<std::string> s = {"apple", "pen", "apple", "pineapple"});
    SHOW(range_distinct_offline(s, std::vector<std::pair<int, int>>{{0, 3}, {0, 4}}), "");
}

}  // namespace rangedistinct_demo

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
    rangedistinct_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    rangedistinct_test::benchmark();    // 速度計測
    rangedistinct_demo::demo();         // 使い方デモ
    return 0;
}
