#ifndef GG_COMBINATORICS
#define GG_COMBINATORICS
// ---- 数え上げ（数学第二世代）  スターリング数・ベル数・分割数・完全順列・全射・写像 12 相・包除原理・二項変換・組分けの列挙・q 二項係数・ヤング盤 ----
// mint の mod は素数（Binomial を使う）。_table は O(N^2) の DP、_row は畳み込みで O(n log n)
namespace math_internal {
template <class mint>
Binomial<mint>& comb_binom() {  // ファイルの中で使い回す階乗の表（必要なら自動で伸びる）
    static Binomial<mint> B(0);
    return B;
}
// T が整数型なら溢れを調べる足し算・掛け算（long long の表で黙って溢れないように）
template <class T>
T cb_add(const T& a, const T& b, const char* f) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail(f, "値が整数型に収まらない（mint で使う）");
        return r;
    } else {
        return a + b;
    }
}
template <class T>
T cb_mul(const T& a, const T& b, const char* f) {
    if constexpr (std::is_integral<T>::value) {
        T r;
        if (__builtin_mul_overflow(a, b, &r)) fail(f, "値が整数型に収まらない（mint で使う）");
        return r;
    } else {
        return a * b;
    }
}
template <class T>
using CbAcc = typename std::conditional<std::is_integral<T>::value, Int128, T>::type;  // 符号が交互の和は Int128 で足してから戻す
template <class T>
T cb_from_acc(const CbAcc<T>& v, const char* f) {
    if constexpr (std::is_integral<T>::value) {
        if (v < (Int128)std::numeric_limits<T>::min() || v > (Int128)std::numeric_limits<T>::max()) fail(f, "値が整数型に収まらない（mint で使う）");
        return (T)v;
    } else {
        return v;
    }
}
inline void cb_check_n(long long N, const char* f) {
    if (N < 0) fail(f, "N ≥ 0 で使う");
}
template <class mint>
std::vector<mint> comb_series_inv(const std::vector<mint>& Q, size_t n) {  // 1/Q の最初の n 項（Q[0] = 1。掛け算だけのニュートン法）
    std::vector<mint> g{mint(1)};
    for (size_t k = 1; k < n;) {
        k *= 2;
        std::vector<mint> q(Q.begin(), Q.begin() + (long long)std::min(Q.size(), k));
        std::vector<mint> t = convolve(q, g);  // g ← g (2 − Q g)
        t.resize(k);
        for (auto& x : t) x = -x;
        t[0] += mint(2);
        g = convolve(g, t);
        g.resize(k);
    }
    g.resize(n);
    return g;
}
template <class mint>
mint comb_falling(long long a, long long n) {  // a (a−1) … (a−n+1)
    mint r = 1;
    for (long long i = 0; i < n; ++i) r *= mint(a - i);
    return r;
}
}  // namespace math_internal

// ---- 第 2 種スターリング数 S(n, k)：n 個の区別できるものを k 個の空でない（区別しない）組に分ける数 ----
template <class mint>
mint stirling2(long long n, long long k) {  // 1 つだけ：Σ_{i=0}^{k} (−1)^{k−i} C(k, i) i^n / k!  O(k log n)
    if (n < 0 || k < 0 || k > n) return mint(0);
    auto& B = math_internal::comb_binom<mint>();
    mint s = 0;
    for (long long i = 0; i <= k; ++i) {
        const mint t = B.C(k, i) * mint(i).pow(n);
        s += ((k - i) % 2 ? -t : t);
    }
    return s * B.inv_fact(k);
}
template <class mint>
std::vector<mint> stirling2_row(int n) {  // S(n, 0), …, S(n, n)（畳み込み 1 回）  O(n log n)
    math_internal::cb_check_n(n, "stirling2_row");
    auto& B = math_internal::comb_binom<mint>();
    std::vector<mint> a((size_t)n + 1), b((size_t)n + 1);
    for (int i = 0; i <= n; ++i) a[(size_t)i] = mint(i).pow(n) * B.inv_fact(i), b[(size_t)i] = (i % 2 ? -B.inv_fact(i) : B.inv_fact(i));
    std::vector<mint> c = convolve(a, b);
    c.resize((size_t)n + 1);
    return c;
}
template <class T>
std::vector<std::vector<T>> stirling2_table(int N) {  // t[n][k]（0 ≤ k ≤ n ≤ N。S(n, k) = k S(n−1, k) + S(n−1, k−1)。整数型は溢れたら止める）  O(N^2)
    math_internal::cb_check_n(N, "stirling2_table");
    std::vector<std::vector<T>> t((size_t)N + 1);
    for (int n = 0; n <= N; ++n) {
        t[(size_t)n].assign((size_t)n + 1, T(0));
        t[(size_t)n][0] = (n == 0 ? T(1) : T(0));
        for (int k = 1; k <= n; ++k)
            t[(size_t)n][(size_t)k] = math_internal::cb_add(k < n ? math_internal::cb_mul(t[(size_t)n - 1][(size_t)k], T(k), "stirling2_table") : T(0), t[(size_t)n - 1][(size_t)k - 1], "stirling2_table");
    }
    return t;
}
// ---- 第 1 種スターリング数（符号なし）c(n, k)：n 個の順列でサイクルがちょうど k 個のものの数（符号つきは (−1)^{n−k} 倍） ----
template <class mint>
std::vector<mint> stirling1_row(int n) {  // c(n, 0), …, c(n, n) = x(x+1)…(x+n−1) の係数（分割統治）  O(n log^2 n)
    math_internal::cb_check_n(n, "stirling1_row");
    std::vector<std::vector<mint>> ps;
    for (int i = 0; i < n; ++i) ps.push_back({mint(i), mint(1)});
    return convolve_many(std::move(ps));
}
template <class T>
std::vector<std::vector<T>> stirling1_table(int N) {  // t[n][k] = c(n, k)（c(n, k) = (n−1) c(n−1, k) + c(n−1, k−1)）  O(N^2)
    math_internal::cb_check_n(N, "stirling1_table");
    std::vector<std::vector<T>> t((size_t)N + 1);
    for (int n = 0; n <= N; ++n) {
        t[(size_t)n].assign((size_t)n + 1, T(0));
        t[(size_t)n][0] = (n == 0 ? T(1) : T(0));
        for (int k = 1; k <= n; ++k)
            t[(size_t)n][(size_t)k] = math_internal::cb_add(k < n ? math_internal::cb_mul(t[(size_t)n - 1][(size_t)k], T(n - 1), "stirling1_table") : T(0), t[(size_t)n - 1][(size_t)k - 1], "stirling1_table");
    }
    return t;
}
// ---- ベル数 B_n（n 個を空でない組にいくつかに分ける数）----
template <class mint>
mint bell(int n) {  // B_n = Σ_k S(n, k)  O(n log n)
    mint s = 0;
    for (const mint& x : stirling2_row<mint>(n)) s += x;
    return s;
}
template <class T>
std::vector<T> bell_table(int N) {  // B_0, …, B_N（ベルの三角形。n − 1 行目の最後が B_n）  O(N^2)
    math_internal::cb_check_n(N, "bell_table");
    std::vector<T> res((size_t)N + 1), row{T(1)};
    res[0] = T(1);
    for (int n = 1; n <= N; ++n) {
        res[(size_t)n] = row.back();
        if (n == N) break;
        std::vector<T> nxt{row.back()};
        for (const T& x : row) nxt.push_back(math_internal::cb_add(nxt.back(), x, "bell_table"));
        row.swap(nxt);
    }
    return res;
}
// ---- 分割数 ----
template <class T>
std::vector<T> partition_table(int N) {  // p(0), …, p(N)：n を正の整数の和に分ける数（順番は区別しない。五角数定理）  O(N √N)
    math_internal::cb_check_n(N, "partition_table");
    std::vector<T> p((size_t)N + 1, T(0));
    p[0] = T(1);
    for (int n = 1; n <= N; ++n) {
        math_internal::CbAcc<T> s = 0;
        for (long long j = 1;; ++j) {
            const long long a = j * (3 * j - 1) / 2, b = j * (3 * j + 1) / 2;
            if (a > n) break;
            const bool plus = j % 2 == 1;
            s = plus ? s + p[(size_t)(n - a)] : s - p[(size_t)(n - a)];
            if (b <= n) s = plus ? s + p[(size_t)(n - b)] : s - p[(size_t)(n - b)];
        }
        p[(size_t)n] = math_internal::cb_from_acc<T>(s, "partition_table");
    }
    return p;
}
template <class mint>
std::vector<mint> partition_row(int N) {  // p(0), …, p(N) を 1/Π(1 − x^k)（五角数の級数の逆数）で  O(N log N)
    math_internal::cb_check_n(N, "partition_row");
    std::vector<mint> Q((size_t)N + 1, mint(0));  // Π(1 − x^k) = Σ_j (−1)^j x^{j(3j−1)/2}（j は ±）
    Q[0] = 1;
    for (long long j = 1;; ++j) {
        const long long a = j * (3 * j - 1) / 2, b = j * (3 * j + 1) / 2;
        if (a > N) break;
        const mint sg = j % 2 ? mint(-1) : mint(1);
        Q[(size_t)a] += sg;
        if (b <= N) Q[(size_t)b] += sg;
    }
    return math_internal::comb_series_inv(Q, (size_t)N + 1);
}
template <class T>
std::vector<std::vector<T>> partition_count_table(int N, int K) {  // t[k][n]：n を k 個以下に分ける数（= k 以下の数の和にする数）  O(N K)
    math_internal::cb_check_n(N, "partition_count_table"), math_internal::cb_check_n(K, "partition_count_table");
    std::vector<std::vector<T>> t((size_t)K + 1, std::vector<T>((size_t)N + 1, T(0)));
    t[0][0] = T(1);
    for (int k = 1; k <= K; ++k)
        for (int n = 0; n <= N; ++n) t[(size_t)k][(size_t)n] = math_internal::cb_add(t[(size_t)k - 1][(size_t)n], n >= k ? t[(size_t)k][(size_t)(n - k)] : T(0), "partition_count_table");
    return t;
}
// ---- 完全順列・全射 ----
template <class T>
std::vector<T> derangement_table(int N) {  // D_0, …, D_N：どの i も i 番目に来ない順列の数（D_n = (n−1)(D_{n−1} + D_{n−2})）
    math_internal::cb_check_n(N, "derangement_table");
    std::vector<T> d((size_t)std::max(N + 1, 2), T(0));
    d[0] = T(1), d[1] = T(0);
    for (int n = 2; n <= N; ++n)
        d[(size_t)n] = math_internal::cb_mul(T(n - 1), math_internal::cb_add(d[(size_t)n - 1], d[(size_t)n - 2], "derangement_table"), "derangement_table");
    d.resize((size_t)N + 1);
    return d;
}
template <class mint>
mint surjections(long long n, long long k) {  // n 個から k 個への全射の数 k! S(n, k)
    if (n < 0 || k < 0 || k > n) return mint(0);
    return math_internal::comb_binom<mint>().fact(k) * stirling2<mint>(n, k);
}

// ---- 写像 12 相：n 個の玉を k 個の箱に入れる方法の数 ----
// balls_distinct / boxes_distinct：区別するか。rule：0 = 何個でもよい、1 = 各箱に 1 個以下、2 = 各箱に 1 個以上
template <class mint>
mint balls_into_boxes(long long n, long long k, bool balls_distinct, bool boxes_distinct, int rule) {
    if (n < 0 || k < 0 || rule < 0 || rule > 2) math_internal::fail("balls_into_boxes", "n, k ≥ 0、rule は 0・1・2");
    auto& B = math_internal::comb_binom<mint>();
    if (rule == 1 && n > k) return mint(0);
    if (rule == 2 && n < k) return mint(0);
    const bool big = n + k > 20000000;  // k が巨大なら階乗の表を作らず、n 個の積で
    if (balls_distinct && boxes_distinct) {
        if (rule == 0) return mint(k).pow(n);
        if (rule == 1) return big ? math_internal::comb_falling<mint>(k, n) : B.P(k, n);
        return surjections<mint>(n, k);
    }
    if (!balls_distinct && boxes_distinct) {
        if (rule == 0) {
            if (n == 0) return mint(1);
            return big ? math_internal::comb_falling<mint>(n + k - 1, n) * B.inv_fact(n) : B.H(k, n);
        }
        if (rule == 1) return big ? math_internal::comb_falling<mint>(k, n) * B.inv_fact(n) : B.C(k, n);
        return n == 0 ? mint(k == 0 ? 1 : 0) : B.C(n - 1, k - 1);
    }
    if (balls_distinct && !boxes_distinct) {
        if (rule == 1) return mint(1);
        if (rule == 2) return stirling2<mint>(n, k);
        // Σ_{i ≤ k} S(n, i) = Σ_j j^n / j! · Σ_{t ≤ k−j} (−1)^t / t!  O(k log n)
        const long long kk = std::min(n, k);
        std::vector<mint> e((size_t)kk + 1);  // e[t] = Σ_{s ≤ t} (−1)^s / s!
        for (long long t = 0; t <= kk; ++t) e[(size_t)t] = (t ? e[(size_t)t - 1] : mint(0)) + (t % 2 ? -B.inv_fact(t) : B.inv_fact(t));
        mint s = 0;
        for (long long j = 0; j <= kk; ++j) s += mint(j).pow(n) * B.inv_fact(j) * e[(size_t)(kk - j)];
        return s;
    }
    // 玉も箱も区別しない：分割数（m を k 個以下に分ける数）
    if (rule == 1) return mint(1);
    const long long m = rule == 2 ? n - k : n;  // 1 個以上なら、先に 1 個ずつ入れて残りを k 個以下に
    if (k == 0) return mint(m == 0 ? 1 : 0);
    if (k >= m) {  // 制限なし：分割数 p(m)
        if (m > 10000000) math_internal::fail("balls_into_boxes", "玉も箱も区別しない場合は n ≤ 10^7 程度");
        return partition_row<mint>((int)m)[(size_t)m];
    }
    if ((Int128)(m + 1) * k > 300000000) math_internal::fail("balls_into_boxes", "玉も箱も区別せず k < n の場合は n·k ≤ 3×10^8 程度");
    std::vector<mint> dp((size_t)m + 1, mint(0));  // 1〜k の数を何個でも使って和を作る
    dp[0] = 1;
    for (long long j = 1; j <= k; ++j)
        for (long long v = j; v <= m; ++v) dp[(size_t)v] += dp[(size_t)(v - j)];
    return dp[(size_t)m];
}

// ---- 包除原理 ----
// Σ_{S ⊆ {0, …, n−1}} (−1)^{|S|} f(S)（n ≤ 30。f は集合 S（ビット）を受け取る。「どの条件も満たさない数」を「条件 S を全部満たす数」から）  O(2^n)
template <class T, class F>
T inclusion_exclusion(int n, F f) {
    if (n < 0 || n > 30) math_internal::fail("inclusion_exclusion", "0 ≤ n ≤ 30 で使う");
    math_internal::CbAcc<T> s = 0;  // 整数型なら Int128 で足して、最後に収まるか確かめる
    for (long long S = 0; S < (1LL << n); ++S) {
        const T v = f(S);
        s = (__builtin_popcountll((unsigned long long)S) % 2) ? s - v : s + v;
    }
    return math_internal::cb_from_acc<T>(s, "inclusion_exclusion");
}
// 対称な包除原理：g(m) = 「n 個の候補のうち、ある m 個だけを使ってよい」ときの数 → n 個を全部使う数 = Σ_{i=0}^{n} (−1)^i C(n, i) g(n − i)  O(n)
template <class mint, class G>
mint inclusion_exclusion_count(int n, G g) {
    math_internal::cb_check_n(n, "inclusion_exclusion_count");
    auto& B = math_internal::comb_binom<mint>();
    mint s = 0;
    for (int i = 0; i <= n; ++i) {
        const mint t = B.C(n, i) * mint(g(n - i));
        s += (i % 2 ? -t : t);
    }
    return s;
}
// ---- 二項変換・「ちょうど k 個」と「k 個を選んで満たす」の変換（mint。畳み込み 1 回）----
namespace math_internal {
template <class mint>
std::vector<mint> cb_egf_mul(const std::vector<mint>& g, bool alt) {  // f(n) = Σ_k C(n, k) g(k) (−1)^{(n−k)·alt}
    const size_t n = g.size();
    if (n == 0) return {};
    auto& B = comb_binom<mint>();
    std::vector<mint> a(n), b(n);
    for (size_t k = 0; k < n; ++k) a[k] = g[k] * B.inv_fact((long long)k), b[k] = (alt && k % 2) ? -B.inv_fact((long long)k) : B.inv_fact((long long)k);
    std::vector<mint> c = convolve(a, b);
    c.resize(n);
    for (size_t k = 0; k < n; ++k) c[k] *= B.fact((long long)k);
    return c;
}
template <class mint>
std::vector<mint> cb_upper(const std::vector<mint>& E, bool alt) {  // A[k] = Σ_{j≥k} C(j, k) E[j] (−1)^{(j−k)·alt}
    const size_t n = E.size();
    if (n == 0) return {};
    auto& B = comb_binom<mint>();
    std::vector<mint> R(n), b(n);
    for (size_t i = 0; i < n; ++i) R[i] = B.fact((long long)(n - 1 - i)) * E[n - 1 - i], b[i] = (alt && i % 2) ? -B.inv_fact((long long)i) : B.inv_fact((long long)i);
    const std::vector<mint> c = convolve(R, b);
    std::vector<mint> A(n);
    for (size_t k = 0; k < n; ++k) A[k] = c[n - 1 - k] * B.inv_fact((long long)k);
    return A;
}
}  // namespace math_internal
// 二項変換 f(n) = Σ_{k=0}^{n} C(n, k) g(k)（n = 0〜|g|−1）  O(n log n)
template <class mint>
std::vector<mint> binomial_transform(const std::vector<mint>& g) { return math_internal::cb_egf_mul(g, false); }
// その逆 g(n) = Σ_{k=0}^{n} (−1)^{n−k} C(n, k) f(k)（「n 個の候補のうち使ってよいもの」f から「n 個全部使う」g へ。包除原理）
template <class mint>
std::vector<mint> binomial_inverse(const std::vector<mint>& f) { return math_internal::cb_egf_mul(f, true); }
// E[j] = ちょうど j 個の条件を満たすものの数 → A[k] = Σ_{j≥k} C(j, k) E[j]（「k 個の条件を選び、それを全部満たす」の和）
template <class mint>
std::vector<mint> at_least_from_exactly(const std::vector<mint>& E) { return math_internal::cb_upper(E, false); }
// その逆：A[k]（k 個の条件を選んで全部満たす、の和。数えやすい）→ E[k] = Σ_{j≥k} (−1)^{j−k} C(j, k) A[j]（ちょうど k 個）
template <class mint>
std::vector<mint> exactly_from_at_least(const std::vector<mint>& A) { return math_internal::cb_upper(A, true); }

// [l, r) の整数で、ds のどれかで割り切れるものの数（ds ≥ 1。包除原理。lcm が範囲を超えた組はまとめて打ち切る）  O(2^|ds|) 以下
inline long long count_divisible_by_any(long long l, long long r, std::vector<long long> ds) {
    for (long long d : ds)
        if (d < 1) math_internal::fail("count_divisible_by_any", "ds ≥ 1 で使う");
    if (r <= l || ds.empty()) return 0;
    // 他の数の倍数になっている数はいらない（その倍数は全部もう数えている）
    std::sort(ds.begin(), ds.end());
    ds.erase(std::unique(ds.begin(), ds.end()), ds.end());
    {
        std::vector<long long> keep;
        for (long long d : ds) {
            bool red = false;
            for (long long e : keep)
                if (d % e == 0) { red = true; break; }
            if (!red) keep.push_back(d);
        }
        ds.swap(keep);
    }
    using U = unsigned long long;
    const U ub = std::max(math_internal::uabs(l), math_internal::uabs(r));  // lcm がこれより大きいと、範囲にある倍数は 0 だけ
    const int n = (int)ds.size();
    const long long zero = (l <= 0 && 0 < r) ? 1 : 0;
    Int128 total = 0;
    auto dfs = [&](auto&& self, int i, U lc, int cnt) -> void {
        if (i == n) {
            if (cnt > 0) {
                const Int128 c = math_internal::ceil_div128(r, (Int128)lc) - math_internal::ceil_div128(l, (Int128)lc);
                total += (cnt % 2 ? c : -c);
            }
            return;
        }
        self(self, i + 1, lc, cnt);
        const U d = (U)ds[(size_t)i], g = math_internal::ugcd(lc, d);
        const UInt128 nl = (UInt128)(lc / g) * d;
        if (nl > ub) {  // この先の組はみな「0 だけ」なので、残りを足し引きすると打ち消し合う（最後の 1 個のときだけ残る）
            if (i + 1 == n) total += ((cnt + 1) % 2 ? zero : -zero);
            return;
        }
        self(self, i + 1, (U)nl, cnt + 1);
    };
    dfs(dfs, 0, 1, 0);
    if (total > (Int128)LLONG_MAX) math_internal::fail("count_divisible_by_any", "答えが long long に収まらない");
    return (long long)total;
}
inline long long count_divisible_by_any_closed(long long L, long long R, const std::vector<long long>& ds) {  // [L, R]
    if (R < L) return 0;
    if (R < LLONG_MAX) return count_divisible_by_any(L, R + 1, ds);
    long long c = count_divisible_by_any(L, R, ds);
    for (long long d : ds)
        if (d >= 1 && R % d == 0) return c + 1;
    return c;
}

// ---- 分け方を全部たどる ----
// {0, …, n−1} を空でない組に分ける方法をすべて f(g, k) に渡す（g[i] = i の組の番号（初めて出る順に 0, 1, …）、k = 組の数）
//   ベル数 B_n 通り（n = 10 で 11 万、12 で 420 万、13 で 2760 万）  O(B_n · n)
template <class F>
void for_each_set_partition(int n, F f) {
    if (n < 0 || n > 15) math_internal::fail("for_each_set_partition", "0 ≤ n ≤ 15 で使う（B_15 は 13 億）");
    std::vector<int> g((size_t)n, 0);
    auto rec = [&](auto&& self, int i, int k) -> void {
        if (i == n) {
            f((const std::vector<int>&)g, k);
            return;
        }
        for (int b = 0; b <= k; ++b) g[(size_t)i] = b, self(self, i + 1, b == k ? k + 1 : k);
    };
    rec(rec, 0, 0);
}
// n を正の整数の和に分ける方法（順番は区別しない）をすべて f(parts) に渡す（parts は大きい順。max_part を渡すと各部分は max_part 以下）
//   分割数 p(n) 通り（n = 50 で 20 万、70 で 400 万）  O(p(n) · n)
template <class F>
void for_each_integer_partition(int n, F f, int max_part = -1) {
    if (n < 0) math_internal::fail("for_each_integer_partition", "n ≥ 0 で使う");
    if (max_part < 0 || max_part > n) max_part = n;
    std::vector<int> parts;
    if (n == 0) {
        f((const std::vector<int>&)parts);
        return;
    }
    struct Frame {
        int rem, x;  // 残り rem を、次は x 以下の部分から試す
    };
    std::vector<Frame> st;  // 再帰しない（部分の数が多くても深くならない）
    st.push_back(Frame{n, std::min(n, max_part)});
    while (!st.empty()) {
        Frame& t = st.back();
        if (t.x == 0) {  // この段は全部試した：ひとつ前の段で選んだ部分を戻す
            st.pop_back();
            if (!parts.empty()) parts.pop_back();
            continue;
        }
        const int x = t.x--, rem = t.rem - x;
        parts.push_back(x);
        if (rem == 0) {
            f((const std::vector<int>&)parts);
            parts.pop_back();
            continue;
        }
        st.push_back({rem, std::min(rem, x)});
    }
}

// ---- q 二項係数・ヤング盤 ----
// q 二項係数 [n, k]_q = Π_{i=0}^{k−1} (1 − q^{n−i}) / (1 − q^{i+1})（q が素数べきなら F_q^n の k 次元部分空間の数）
//   q は何でも（q^d ≡ 1 となる d があれば q-Lucas の定理で）。n ≤ 10^7  O(n)
template <class mint>
mint q_binomial(long long n, long long k, long long q) {
    if (n < 0 || k < 0 || k > n) return mint(0);
    if (n > 10000000) math_internal::fail("q_binomial", "n ≤ 10^7 で使う");
    k = std::min(k, n - k);
    const mint Q(q);
    if (Q == mint(0)) return mint(1);  // [n, k]_0 = 1
    long long d = 0;                   // q の位数（q^d = 1 となる最小の d ≤ n。なければ 0）
    {
        mint t = 1;
        for (long long j = 1; j <= n; ++j)
            if ((t *= Q) == mint(1)) {
                d = j;
                break;
            }
    }
    auto direct = [&](long long a, long long b) {  // [a, b]_q（b ≤ a で分母が 0 にならないとき）
        mint num = 1, den = 1, x = Q.pow(a - b + 1), y = Q;
        for (long long i = 0; i < b; ++i) num *= mint(1) - x, den *= mint(1) - y, x *= Q, y *= Q;
        return num / den;
    };
    if (d == 0) return direct(n, k);
    // q-Lucas：[n, k]_q ≡ C(⌊n/d⌋, ⌊k/d⌋) · [n mod d, k mod d]_q。C は mod が小さいときは Lucas の定理で
    auto& B = math_internal::comb_binom<mint>();
    const long long p = (long long)mint::mod();
    mint c = 1;
    for (long long a = n / d, b = k / d; (a > 0 || b > 0) && c != mint(0); a /= p, b /= p) c *= B.C(a % p, b % p);
    return c == mint(0) ? c : c * direct(n % d, k % d);
}
// 形 shape（各行の長さ。大きい順）の標準ヤング盤の数 = n! / Π（各マスのフックの長さ）（フック長公式。n はマスの数）  O(n)
template <class mint, class V>
mint young_tableaux_count(const V& shape) {
    std::vector<long long> row;
    for (const auto& x : shape) row.push_back(math_internal::to_ll_checked(x, "young_tableaux_count"));
    while (!row.empty() && row.back() == 0) row.pop_back();
    long long n = 0;
    for (size_t i = 0; i < row.size(); ++i) {
        if (row[i] < 0 || (i > 0 && row[i] > row[i - 1])) math_internal::fail("young_tableaux_count", "行の長さは 0 以上で大きい順に並べる");
        if (row[i] > 100000000 - n) math_internal::fail("young_tableaux_count", "マスが多すぎる（10^8 まで）");
        n += row[i];
    }
    std::vector<long long> col(row.empty() ? 0 : (size_t)row[0], 0);  // col[j] = j 列目の長さ
    for (long long r : row)
        for (long long j = 0; j < r; ++j) ++col[(size_t)j];
    auto& B = math_internal::comb_binom<mint>();
    mint prod = 1;
    for (size_t i = 0; i < row.size(); ++i)
        for (long long j = 0; j < row[i]; ++j) prod *= mint((row[i] - j - 1) + (col[(size_t)j] - (long long)i - 1) + 1);
    return B.fact(n) / prod;
}
template <class mint>
mint young_tableaux_count(std::initializer_list<long long> shape) { return young_tableaux_count<mint, std::initializer_list<long long>>(shape); }
#endif  // GG_COMBINATORICS
