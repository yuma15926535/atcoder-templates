// =============================================================================
//  半分全列挙（meet in the middle：N ≤ 40 の部分集合・グリッドの経路を半分ずつ列挙して組み合わせる）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  2^N 通りを全部見るのが間に合わない（N = 40 → 10^12）とき、前半と後半に分けて 2^(N/2) 通りずつ列挙し、
//  ソートしてから二分探索・尺取りで組み合わせる。O(2^(N/2) × N) 程度
//    例: ABC184 F（Programming Contest）… 和が T 以下の最大 → mitm_max_le(A, T)
//        ABC271 F（XOR on Grid Path）  … 経路の XOR が 0 の数 → grid_path_count_xor(A, 0)
//
//  機能一覧（a は vector<long long> など。部分集合は空集合も含む。N は a の長さ）               計算量
//   部分集合の列挙
//    subset_sums_sorted(a)                 全部分集合の和（2^N 個・小さい順。同じ値も残す）           O(2^N)
//    subset_fold(a, op, e)                 f[mask] = mask に入る要素を op でまとめた値（mask の順）  O(2^N)
//   半分全列挙（部分集合の和）
//    mitm_max_le(a, S) / mitm_min_ge(a, S) 和が S 以下の最大 / S 以上の最小（無ければ nullopt）        O(2^(N/2) N)
//    mitm_count_eq(a, S)                   和がちょうど S の部分集合の数                              O(2^(N/2) N)
//    mitm_count_le(a, S) / mitm_count_lt   和が S 以下 / S 未満の部分集合の数
//    mitm_count_range(a, L, R)             和が [L, R) に入る数（閉区間は mitm_count_range_closed）
//    mitm_min_abs_diff(a, S)               |和 − S| の最小
//    mitm_find_eq(a, S)                    和がちょうど S の部分集合を 1 つ（mask。N ≤ 64）。無ければ nullopt
//    mitm_count_xor(a, X)                  XOR がちょうど X の部分集合の数（a は 0 以上）              O(2^(N/2) N)
//    mitm_knapsack(w, v, W)                重さの和 ≤ W で価値の和の最大（N ≤ 40、W が大きい）。W < 0 なら nullopt
//   グリッドの経路（(0, 0) から (H−1, W−1) へ右か下に進む。H + W − 2 ≤ 40 くらい）
//    grid_path_count_xor(A, X)             通るマス（両端を含む）の XOR が X の経路の数                O(2^((H+W)/2) (H+W))
//    grid_path_count_sum(A, S)             通るマスの和が S の経路の数
//    grid_path_count_meet(A, op, e, need)  一般形：前半（真ん中の斜めのマスを含む）の値 L、後半の値 R として
//                                          R == need(L) となる組の数。op は結合的・e は単位元（後半は op(マス, R) の順にまとめる）
//
//  注意
//    和は T（a の要素の型）で計算する。|和| が T に収まるように（long long なら |a_i| ≤ 2×10^17 で N = 40 まで大丈夫）
//    負の値も使える（mitm_count_xor・grid_path_count_xor は 0 以上の整数）
//    数は long long（最大 2^40 くらい）。N が 50 を超えると時間もメモリも足りない
//  速さの目安: N = 40 で mitm_count_le・mitm_max_le・mitm_knapsack それぞれ 約 0.05 秒
//             grid_path_count_xor（20 × 20）約 0.15 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mitm_test::random_test()  … 全部の部分集合・全部の経路を調べる愚直解との突き合わせ
//    mitm_test::benchmark()    … 速度計測
//    mitm_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace mitm_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// 全部分集合の和を小さい順に（マージで O(2^n)）
template <class T>
std::vector<T> sums_sorted(const T* a, int n) {
    if (n > 30) fail("subset_sums_sorted", "要素が多すぎる（2^n 個を作る。n ≤ 30 くらいまで）");
    std::vector<T> v{T(0)}, w, x;
    v.reserve((size_t)1 << n), w.reserve((size_t)1 << n), x.reserve((size_t)1 << n);
    for (int i = 0; i < n; ++i) {
        x.resize(v.size());
        for (size_t j = 0; j < v.size(); ++j) x[j] = v[j] + a[i];
        w.resize(v.size() * 2);
        std::merge(v.begin(), v.end(), x.begin(), x.end(), w.begin());
        std::swap(v, w);
    }
    return v;
}
// 和と mask の組を和の小さい順に（mask は offset ビット目から）
template <class T>
std::vector<std::pair<T, unsigned long long>> sums_with_mask(const T* a, int n, int offset) {
    std::vector<std::pair<T, unsigned long long>> v{{T(0), 0ULL}}, w, x;
    for (int i = 0; i < n; ++i) {
        x.resize(v.size());
        for (size_t j = 0; j < v.size(); ++j) x[j] = {v[j].first + a[i], v[j].second | (1ULL << (i + offset))};
        w.resize(v.size() * 2);
        std::merge(v.begin(), v.end(), x.begin(), x.end(), w.begin(),
                   [](const auto& p, const auto& q) { return p.first < q.first; });
        std::swap(v, w);
    }
    return v;
}
template <class T>
std::pair<std::vector<T>, std::vector<T>> halves(const std::vector<T>& a, const char* func) {
    if (a.size() > 62) fail(func, "要素が多すぎる（半分ずつでも 2^(N/2) 個を作る）");
    const int n = (int)a.size(), h = n / 2;
    return {sums_sorted(a.data(), h), sums_sorted(a.data() + h, n - h)};
}
// L, R（どちらも小さい順）で l + r ≤ S（strict なら <）の組の数
template <class T>
long long count_pair_le(const std::vector<T>& L, const std::vector<T>& R, T S, bool strict) {
    long long c = 0;
    size_t j = R.size();
    for (const T& l : L) {
        while (j > 0 && (strict ? !(l + R[j - 1] < S) : (S < l + R[j - 1]))) --j;
        c += (long long)j;
    }
    return c;
}
}  // namespace mitm_internal

// ---- 部分集合の列挙 ----
template <class T>
std::vector<T> subset_sums_sorted(const std::vector<T>& a) {
    return mitm_internal::sums_sorted(a.data(), (int)a.size());
}
template <class T, class Op>
std::vector<T> subset_fold(const std::vector<T>& a, Op op, T e) {
    const int n = (int)a.size();
    if (n > 30) mitm_internal::fail("subset_fold", "要素が多すぎる（2^n 個を作る。n ≤ 30 くらいまで）");
    std::vector<T> f((size_t)1 << n);
    f[0] = e;
    for (int i = 0; i < n; ++i)
        for (size_t m = 0; m < ((size_t)1 << i); ++m) f[m | ((size_t)1 << i)] = op(f[m], a[i]);
    return f;
}

// ---- 半分全列挙（部分集合の和）----
template <class T>
std::optional<T> mitm_max_le(const std::vector<T>& a, T S) {
    auto [L, R] = mitm_internal::halves(a, "mitm_max_le");
    std::optional<T> best;
    size_t j = R.size();
    for (const T& l : L) {
        while (j > 0 && S < l + R[j - 1]) --j;
        if (j == 0) break;
        const T v = l + R[j - 1];
        if (!best || *best < v) best = v;
    }
    return best;
}
template <class T>
std::optional<T> mitm_min_ge(const std::vector<T>& a, T S) {
    auto [L, R] = mitm_internal::halves(a, "mitm_min_ge");
    std::optional<T> best;
    size_t j = 0;  // R[j] が l + R[j] ≥ S となる最初
    for (size_t i = L.size(); i-- > 0;) {
        const T& l = L[i];
        while (j < R.size() && l + R[j] < S) ++j;
        if (j == R.size()) break;
        const T v = l + R[j];
        if (!best || v < *best) best = v;
    }
    return best;
}
template <class T>
long long mitm_count_le(const std::vector<T>& a, T S) {
    auto [L, R] = mitm_internal::halves(a, "mitm_count_le");
    return mitm_internal::count_pair_le(L, R, S, false);
}
template <class T>
long long mitm_count_lt(const std::vector<T>& a, T S) {
    auto [L, R] = mitm_internal::halves(a, "mitm_count_lt");
    return mitm_internal::count_pair_le(L, R, S, true);
}
template <class T>
long long mitm_count_range(const std::vector<T>& a, T lo, T hi) {  // 和 ∈ [lo, hi)
    if (!(lo < hi)) return 0;
    auto [L, R] = mitm_internal::halves(a, "mitm_count_range");
    return mitm_internal::count_pair_le(L, R, hi, true) - mitm_internal::count_pair_le(L, R, lo, true);
}
template <class T>
long long mitm_count_range_closed(const std::vector<T>& a, T lo, T hi) {  // 和 ∈ [lo, hi]
    if (hi < lo) return 0;
    auto [L, R] = mitm_internal::halves(a, "mitm_count_range_closed");
    return mitm_internal::count_pair_le(L, R, hi, false) - mitm_internal::count_pair_le(L, R, lo, true);
}
template <class T>
long long mitm_count_eq(const std::vector<T>& a, T S) {
    auto [L, R] = mitm_internal::halves(a, "mitm_count_eq");
    long long c = 0;
    size_t i = 0, j = R.size();
    while (i < L.size() && j > 0) {
        const T v = L[i] + R[j - 1];
        if (v < S) ++i;
        else if (S < v) --j;
        else {
            size_t i2 = i, j2 = j;
            while (i2 < L.size() && L[i2] == L[i]) ++i2;
            while (j2 > 0 && R[j2 - 1] == R[j - 1]) --j2;
            c += (long long)(i2 - i) * (long long)(j - j2);
            i = i2, j = j2;
        }
    }
    return c;
}
template <class T>
T mitm_min_abs_diff(const std::vector<T>& a, T S) {
    auto [L, R] = mitm_internal::halves(a, "mitm_min_abs_diff");
    bool first = true;
    T best{};
    size_t j = R.size();  // R[j-1] が l + R[j-1] ≤ S となる最後
    for (const T& l : L) {
        while (j > 0 && S < l + R[j - 1]) --j;
        if (j > 0) {
            const T d = S - (l + R[j - 1]);
            if (first || d < best) best = d, first = false;
        }
        if (j < R.size()) {
            const T d = (l + R[j]) - S;
            if (first || d < best) best = d, first = false;
        }
    }
    return best;
}
template <class T>
std::optional<unsigned long long> mitm_find_eq(const std::vector<T>& a, T S) {
    if (a.size() > 62) mitm_internal::fail("mitm_find_eq", "要素が多すぎる");
    const int n = (int)a.size(), h = n / 2;
    const auto L = mitm_internal::sums_with_mask(a.data(), h, 0);
    const auto R = mitm_internal::sums_with_mask(a.data() + h, n - h, h);
    size_t i = 0, j = R.size();
    while (i < L.size() && j > 0) {
        const T v = L[i].first + R[j - 1].first;
        if (v < S) ++i;
        else if (S < v) --j;
        else return L[i].second | R[j - 1].second;
    }
    return std::nullopt;
}
template <class T>
long long mitm_count_xor(const std::vector<T>& a, T X) {
    if (a.size() > 62) mitm_internal::fail("mitm_count_xor", "要素が多すぎる");
    const int n = (int)a.size(), h = n / 2;
    auto fold = [](const T* p, int m) {
        std::vector<T> f((size_t)1 << m);
        f[0] = T(0);
        for (int i = 0; i < m; ++i)
            for (size_t k = 0; k < ((size_t)1 << i); ++k) f[k | ((size_t)1 << i)] = f[k] ^ p[i];
        std::sort(f.begin(), f.end());
        return f;
    };
    const std::vector<T> L = fold(a.data(), h), R = fold(a.data() + h, n - h);
    long long c = 0;
    for (size_t i = 0; i < L.size();) {
        size_t i2 = i;
        while (i2 < L.size() && L[i2] == L[i]) ++i2;
        const T need = L[i] ^ X;
        const auto r = std::equal_range(R.begin(), R.end(), need);
        c += (long long)(i2 - i) * (long long)(r.second - r.first);
        i = i2;
    }
    return c;
}
// 重さ w・価値 v の品物から、重さの和 ≤ W で価値の和の最大（W < 0 なら nullopt。空集合は価値 0）
template <class T>
std::optional<T> mitm_knapsack(const std::vector<T>& w, const std::vector<T>& v, T W) {
    if (w.size() != v.size()) mitm_internal::fail("mitm_knapsack", "w と v の長さが違う");
    if (w.size() > 62) mitm_internal::fail("mitm_knapsack", "要素が多すぎる");
    if (W < T(0)) return std::nullopt;
    const int n = (int)w.size(), h = n / 2;
    auto items = [&](int s, int m) {  // (重さ, 価値) の全部分集合を重さの小さい順に（マージで O(2^m)）
        std::vector<std::pair<T, T>> f{{T(0), T(0)}}, g, x;
        f.reserve((size_t)1 << m), g.reserve((size_t)1 << m), x.reserve((size_t)1 << m);
        for (int i = 0; i < m; ++i) {
            x.resize(f.size());
            for (size_t k = 0; k < f.size(); ++k) x[k] = {f[k].first + w[s + i], f[k].second + v[s + i]};
            g.resize(f.size() * 2);
            std::merge(f.begin(), f.end(), x.begin(), x.end(), g.begin(),
                       [](const std::pair<T, T>& p, const std::pair<T, T>& q) { return p.first < q.first; });
            std::swap(f, g);
        }
        return f;
    };
    const auto L = items(0, h);
    auto R = items(h, n - h);
    // 重さの順に、価値をそれまでの最大にする（重いのに価値が低いものは使わない）
    for (size_t i = 1; i < R.size(); ++i) R[i].second = std::max(R[i].second, R[i - 1].second);
    std::optional<T> best;
    size_t j = R.size();  // R[j-1] が重さ ≤ W − lw の最後
    for (const auto& [lw, lv] : L) {
        if (W < lw) break;
        const T rest = W - lw;
        while (j > 0 && rest < R[j - 1].first) --j;
        if (j == 0) break;  // 重さが負の品物があるときだけ起きうる
        const T val = lv + R[j - 1].second;
        if (!best || *best < val) best = val;
    }
    return best;
}

// ---- グリッドの経路 ----
// (0, 0) → (H−1, W−1)（右か下）。前半 = (0, 0) から斜め i + j = D のマスまで（そのマスを含む）を op でまとめた値 L、
// 後半 = そのマスの次から (H−1, W−1) まで（そのマスを含まない）の値 R。R == need(L) の組の数
template <class T, class Op, class Need>
long long grid_path_count_meet(const std::vector<std::vector<T>>& A, Op op, T e, Need need) {
    const int H = (int)A.size();
    if (H == 0) return 0;
    const int W = (int)A[0].size();
    if (W == 0) return 0;
    for (const auto& r : A)
        if ((int)r.size() != W) mitm_internal::fail("grid_path_count_meet", "行の長さがそろっていない");
    if (H + W - 2 > 50) mitm_internal::fail("grid_path_count_meet", "グリッドが大きすぎる（H + W − 2 ≤ 50 くらいまで）");
    const int D = (H + W - 2) / 2;
    // 前半：斜めのマス（行 i で決まる）ごとに need(L) を集める
    std::vector<std::vector<T>> needs(H);
    auto dfs1 = [&](auto&& self, int i, int j, T val) -> void {
        if (i + j == D) {
            needs[i].push_back(need(val));
            return;
        }
        if (i + 1 < H) self(self, i + 1, j, op(val, A[i + 1][j]));
        if (j + 1 < W) self(self, i, j + 1, op(val, A[i][j + 1]));
    };
    dfs1(dfs1, 0, 0, op(e, A[0][0]));
    for (auto& v : needs) std::sort(v.begin(), v.end());
    long long cnt = 0;
    auto hit = [&](int i, const T& val) {
        const auto r = std::equal_range(needs[i].begin(), needs[i].end(), val);
        cnt += (long long)(r.second - r.first);
    };
    // 後半：(H−1, W−1) から戻る。斜めのマスに着いたら、そのマスの値は入れずに数える
    if (H - 1 + W - 1 == D) {  // H = W = 1
        hit(H - 1, e);
        return cnt;
    }
    auto dfs2 = [&](auto&& self, int i, int j, T val) -> void {  // (i, j) まで入れた値 val
        if (i > 0) {
            if (i - 1 + j == D) hit(i - 1, val);
            else self(self, i - 1, j, op(A[i - 1][j], val));
        }
        if (j > 0) {
            if (i + j - 1 == D) hit(i, val);
            else self(self, i, j - 1, op(A[i][j - 1], val));
        }
    };
    dfs2(dfs2, H - 1, W - 1, op(A[H - 1][W - 1], e));
    return cnt;
}
template <class T>
long long grid_path_count_xor(const std::vector<std::vector<T>>& A, T X) {
    return grid_path_count_meet(A, [](T a, T b) { return a ^ b; }, T(0), [X](T L) { return L ^ X; });
}
template <class T>
long long grid_path_count_sum(const std::vector<std::vector<T>>& A, T S) {
    return grid_path_count_meet(A, [](T a, T b) { return a + b; }, T(0), [S](T L) { return S - L; });
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace mitm_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define MM_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[search_meet_in_middle] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static void test_subsets(std::mt19937_64& rng, long long& cases) {
    const int n = (int)rnd(rng, 0, 14);
    const ll V = rnd(rng, 0, 1) ? 5 : 1000000000000LL;
    const bool neg = rnd(rng, 0, 2) == 0;
    std::vector<ll> a(n);
    for (auto& x : a) x = rnd(rng, neg ? -V : 0, V);
    std::vector<ll> all;
    for (int m = 0; m < (1 << n); ++m) {
        ll s = 0;
        for (int i = 0; i < n; ++i)
            if (m >> i & 1) s += a[i];
        all.push_back(s);
    }
    std::vector<ll> srt = all;
    std::sort(srt.begin(), srt.end());
    MM_CHECK(subset_sums_sorted(a) == srt);
    MM_CHECK(subset_fold(a, [](ll x, ll y) { return x + y; }, 0LL) == all);
    for (int q = 0; q < 6; ++q) {
        ll S = rnd(rng, 0, 1) && !all.empty() ? all[rng() % all.size()] + rnd(rng, -1, 1) : rnd(rng, -V * 3, V * 3 + 1);
        ll lo = rnd(rng, -V * 2, V * 2), hi = rnd(rng, -V * 2, V * 2);
        std::optional<ll> mx, mn;
        ll ceq = 0, cle = 0, clt = 0, crg = 0, crc = 0, mad = -1;
        for (ll s : all) {
            if (s <= S && (!mx || *mx < s)) mx = s;
            if (s >= S && (!mn || s < *mn)) mn = s;
            ceq += s == S, cle += s <= S, clt += s < S;
            crg += lo <= s && s < hi, crc += lo <= s && s <= hi;
            const ll d = s > S ? s - S : S - s;
            if (mad < 0 || d < mad) mad = d;
        }
        MM_CHECK(mitm_max_le(a, S) == mx);
        MM_CHECK(mitm_min_ge(a, S) == mn);
        MM_CHECK(mitm_count_eq(a, S) == ceq);
        MM_CHECK(mitm_count_le(a, S) == cle);
        MM_CHECK(mitm_count_lt(a, S) == clt);
        MM_CHECK(mitm_count_range(a, lo, hi) == crg);
        MM_CHECK(mitm_count_range_closed(a, lo, hi) == crc);
        MM_CHECK(mitm_min_abs_diff(a, S) == mad);
        const auto f = mitm_find_eq(a, S);
        MM_CHECK(f.has_value() == (ceq > 0));
        if (f) {
            ll s = 0;
            MM_CHECK((*f >> n) == 0);
            for (int i = 0; i < n; ++i)
                if (*f >> i & 1) s += a[i];
            MM_CHECK(s == S);
        }
        cases += 9;
    }
    // XOR
    std::vector<ll> b(n);
    const ll B = rnd(rng, 0, 1) ? 7 : (1LL << 40);
    for (auto& x : b) x = rnd(rng, 0, B);
    for (int q = 0; q < 4; ++q) {
        const ll X = rnd(rng, 0, B);
        ll c = 0;
        for (int m = 0; m < (1 << n); ++m) {
            ll x = 0;
            for (int i = 0; i < n; ++i)
                if (m >> i & 1) x ^= b[i];
            c += x == X;
        }
        MM_CHECK(mitm_count_xor(b, X) == c);
        ++cases;
    }
    // ナップサック
    std::vector<ll> w(n), v(n);
    for (int i = 0; i < n; ++i) w[i] = rnd(rng, 0, V), v[i] = rnd(rng, neg ? -V : 0, V);
    for (int q = 0; q < 4; ++q) {
        const ll W = rnd(rng, -2, V * n / 2 + 1);
        std::optional<ll> best;
        if (W >= 0)
            for (int m = 0; m < (1 << n); ++m) {
                ll sw = 0, sv = 0;
                for (int i = 0; i < n; ++i)
                    if (m >> i & 1) sw += w[i], sv += v[i];
                if (sw <= W && (!best || *best < sv)) best = sv;
            }
        MM_CHECK(mitm_knapsack(w, v, W) == best);
        ++cases;
    }
}

static void test_grid(std::mt19937_64& rng, long long& cases) {
    const int H = (int)rnd(rng, 1, 6), W = (int)rnd(rng, 1, 6);
    const ll V = rnd(rng, 0, 1) ? 3 : 1000;
    std::vector<std::vector<ll>> A(H, std::vector<ll>(W));
    for (auto& r : A)
        for (auto& x : r) x = rnd(rng, 0, V);
    std::vector<ll> xs, ss;  // 全経路の XOR と和
    auto dfs = [&](auto&& self, int i, int j, ll x, ll s) -> void {
        x ^= A[i][j], s += A[i][j];
        if (i == H - 1 && j == W - 1) {
            xs.push_back(x), ss.push_back(s);
            return;
        }
        if (i + 1 < H) self(self, i + 1, j, x, s);
        if (j + 1 < W) self(self, i, j + 1, x, s);
    };
    dfs(dfs, 0, 0, 0, 0);
    for (int q = 0; q < 4; ++q) {
        const ll X = rnd(rng, 0, 1) ? xs[rng() % xs.size()] : rnd(rng, 0, V);
        const ll S = rnd(rng, 0, 1) ? ss[rng() % ss.size()] : rnd(rng, 0, V * (H + W));
        MM_CHECK(grid_path_count_xor(A, X) == (ll)std::count(xs.begin(), xs.end(), X));
        MM_CHECK(grid_path_count_sum(A, S) == (ll)std::count(ss.begin(), ss.end(), S));
        cases += 2;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 400 * TEST_SCALE; ++it) test_subsets(rng, cases);
    for (int it = 0; it < 400 * TEST_SCALE; ++it) test_grid(rng, cases);
    // 端：空・1 個
    MM_CHECK(mitm_max_le(std::vector<ll>{}, 0LL) == std::optional<ll>(0) && !mitm_max_le(std::vector<ll>{}, -1LL));
    MM_CHECK(mitm_count_eq(std::vector<ll>{}, 0LL) == 1 && mitm_count_eq(std::vector<ll>{5}, 5LL) == 1);
    MM_CHECK(grid_path_count_xor(std::vector<std::vector<ll>>{{3}}, 3LL) == 1);
    MM_CHECK(grid_path_count_xor(std::vector<std::vector<ll>>{{3}}, 0LL) == 0);
    cases += 4;
    std::cout << "[search_meet_in_middle] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[search_meet_in_middle] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 34 : 40;
    std::vector<ll> a(N), w(N), v(N);
    for (int i = 0; i < N; ++i) a[i] = rnd(rng, 1, 1000000000), w[i] = rnd(rng, 1, 1000000000000LL), v[i] = rnd(rng, 1, 1000000000000LL);
    ll chk = 0;
    auto t0 = std::chrono::steady_clock::now();
    chk += mitm_count_le(a, 10000000000LL);
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    chk += *mitm_max_le(a, 10000000000LL);
    const ll t2 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    chk += *mitm_knapsack(w, v, 10000000000000LL);
    const ll t3 = elapsed_ms(t0);
    const int G = TEST_SCALE == 1 ? 12 : 20;
    std::vector<std::vector<ll>> A(G, std::vector<ll>(G));
    for (auto& r : A)
        for (auto& x : r) x = rnd(rng, 0, (1LL << 30) - 1);
    t0 = std::chrono::steady_clock::now();
    chk += grid_path_count_xor(A, 0LL);
    const ll t4 = elapsed_ms(t0);
    std::cout << "[search_meet_in_middle] benchmark (N = " << N << "): mitm_count_le " << t1 << " ms / mitm_max_le " << t2
              << " ms / mitm_knapsack " << t3 << " ms / grid_path_count_xor(" << G << " × " << G << ") " << t4 << " ms  (check " << chk % 1000
              << ")\n";
}

}  // namespace mitm_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace mitm_demo {

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
#define SHOW(expr, comment) mitm_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    mitm_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("半分全列挙（search_meet_in_middle） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  subset_sums_sorted / subset_fold\n"
              << "  mitm_max_le / mitm_min_ge / mitm_count_eq / mitm_count_le / mitm_count_lt / mitm_count_range(_closed)\n"
              << "  mitm_min_abs_diff / mitm_find_eq / mitm_count_xor / mitm_knapsack\n"
              << "  grid_path_count_xor / grid_path_count_sum / grid_path_count_meet\n";

    api("【1】 subset_sums_sorted(a) / subset_fold(a, op, e)", "vector<T>",
        "全部分集合の和（小さい順） / mask の順に op でまとめた値（mask の i ビット目が 1 ⇔ a[i] を使う）", "O(2^N)");
    RUN(const std::vector<long long> a = {3, 5, 9});
    SHOW(subset_sums_sorted(a), "空集合の 0 も入る");
    SHOW(subset_fold(a, [](long long x, long long y) { return x + y; }, 0LL), "f[mask]：mask = 0〜7 の順");
    SHOW(subset_fold(a, [](long long x, long long y) { return std::max(x, y); }, 0LL), "op = max でも");

    api("【2】 mitm_max_le(a, S) / mitm_min_ge(a, S)", "optional<T>", "和が S 以下の最大 / S 以上の最小（無ければ nullopt）", "O(2^(N/2) N)");
    RUN(const std::vector<long long> b = {2, 7, 11, 13, 19, 23});
    SHOW(mitm_max_le(b, 30LL), "例 7 + 23 = 30");
    SHOW(mitm_max_le(b, 1LL), "空集合の 0");
    SHOW(mitm_max_le(b, -1LL), "空集合でも超える → nullopt");
    SHOW(mitm_min_ge(b, 76LL), "全部足しても 75 → nullopt");

    api("【3】 mitm_count_eq / mitm_count_le / mitm_count_lt / mitm_count_range(a, L, R)", "long long",
        "和がちょうど S / S 以下 / S 未満 / [L, R) に入る部分集合の数（閉区間は mitm_count_range_closed）", "O(2^(N/2) N)");
    SHOW(mitm_count_eq(b, 20LL), "{7, 13}, {2, 7, 11} → 2 個");
    SHOW(mitm_count_le(b, 10LL), "{}, {2}, {7}, {2, 7}");
    SHOW(mitm_count_range(b, 10LL, 20LL), "和が 10〜19");
    SHOW(mitm_count_range_closed(b, 10LL, 20LL), "和が 10〜20");

    api("【4】 mitm_min_abs_diff(a, S) / mitm_find_eq(a, S)", "T / optional<unsigned long long>",
        "|和 − S| の最小 / 和がちょうど S の部分集合を 1 つ（mask）", "O(2^(N/2) N)");
    SHOW(mitm_min_abs_diff(b, 76LL), "全部足しても 75 → 差 1");
    SHOW(mitm_find_eq(b, 20LL), "mask = 0b1010 → b[1] + b[3] = 7 + 13");
    SHOW(mitm_find_eq(b, 1LL), "作れない");

    api("【5】 mitm_count_xor(a, X)", "long long", "XOR がちょうど X の部分集合の数（a は 0 以上）", "O(2^(N/2) N)");
    SHOW(mitm_count_xor(std::vector<long long>{1, 2, 3, 4}, 0LL), "{}, {1, 2, 3} の 2 個");

    api("【6】 mitm_knapsack(w, v, W)", "optional<T>", "重さの和 ≤ W で価値の和の最大（N ≤ 40、W が大きくて普通の DP ができないとき）",
        "O(2^(N/2) N)");
    SHOW(mitm_knapsack(std::vector<long long>{4, 3, 5, 6}, std::vector<long long>{5, 4, 6, 9}, 10LL), "重さ 4 + 6 → 価値 14");

    api("【7】 grid_path_count_xor(A, X) / grid_path_count_sum(A, S)", "long long",
        "(0, 0) から (H−1, W−1) へ右か下に進む経路で、通るマス（両端を含む）の XOR / 和が X / S のものの数",
        "O(2^((H+W)/2) (H+W))");
    RUN(const std::vector<std::vector<long long>> A = {{1, 5, 3}, {7, 2, 6}, {4, 0, 1}});
    SHOW(grid_path_count_xor(A, 3LL), "XOR が 3 の経路（下右右下・下下右右）");
    SHOW(grid_path_count_sum(A, 16LL), "和が 16 の経路（右右下下）");

    api("【8】 grid_path_count_meet(A, op, e, need)", "long long",
        "一般形。前半（真ん中の斜めのマスまで・そのマスを含む）を op でまとめた L、後半の R として R == need(L) の数", "同上");
    code("// 例: 経路の積が P（mod 素数）→ op = 掛け算、e = 1、need(L) = P × L^{-1}");
    SHOW(grid_path_count_meet(A, [](long long x, long long y) { return x ^ y; }, 0LL, [](long long L) { return L ^ 3LL; }),
         "XOR が 3 の経路（xor 版と同じ）");

    api("【9】 使い方の例: ABC184 F（Programming Contest）", "—", "N ≤ 40 個の A_i から選んだ和が T 以下で最大 → mitm_max_le(A, T)", "O(2^(N/2) N)");
    code("std::vector<long long> A(N);  …  cout << *mitm_max_le(A, T) << '\\n';   // 空集合 0 ≤ T なので必ず値がある");
}

}  // namespace mitm_demo

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
    mitm_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mitm_test::benchmark();    // 速度計測
    mitm_demo::demo();         // 使い方デモ
    return 0;
}
