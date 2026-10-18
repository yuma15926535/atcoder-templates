#ifndef GG_DIVISOR_TRANSFORM
#define GG_DIVISOR_TRANSFORM
// ---- 約数・倍数の変換（数学第二世代）  約数・倍数のゼータ・メビウス変換、gcd・lcm の畳み込み、gcd ごとの組の数、Σφ・Σμ、ディリクレ積、M の約数の上の変換 ----
// 配列は 1〜N を使う（f[0] は使わない・変えない）。T は mint・long long（溢れたら止める）・double など
namespace math_internal {
template <class T>
inline T dt_add(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail("divisor_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a + b;
    }
}
template <class T>
inline T dt_sub(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_sub_overflow(a, b, &r)) fail("divisor_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a - b;
    }
}
template <class T>
inline T dt_mul(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_mul_overflow(a, b, &r)) fail("divisor_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a * b;
    }
}
struct DtAdd {
    template <class T>
    T operator()(const T& a, const T& b) const { return dt_add(a, b); }
};
struct DtSub {
    template <class T>
    T operator()(const T& a, const T& b) const { return dt_sub(a, b); }
};
inline const std::vector<int>& dt_primes(int N) {  // N 以下の素数（前より大きい N が来たら作り直す）
    static std::vector<int> ps;
    static int built = -1;
    if (built < N) ps = primes_up_to<int>(std::max(N, 2)), built = std::max(N, 2);
    return ps;
}
}  // namespace math_internal

// g[n] = Σ_{d | n} f[d]（n = 1〜N。op を渡すと Σ の代わり）  O(N log log N)
template <class T, class Op = math_internal::DtAdd>
void divisor_zeta(std::vector<T>& f, Op op = Op{}) {
    const int N = (int)f.size() - 1;
    for (int p : math_internal::dt_primes(N)) {
        if (p > N) break;
        for (int i = 1; (long long)i * p <= N; ++i) f[(size_t)i * p] = op(f[(size_t)i * p], f[(size_t)i]);
    }
}
// divisor_zeta の逆（f[n] = Σ_{d | n} μ(n/d) g[d]）
template <class T, class Op = math_internal::DtSub>
void divisor_mobius(std::vector<T>& f, Op inv_op = Op{}) {
    const int N = (int)f.size() - 1;
    for (int p : math_internal::dt_primes(N)) {
        if (p > N) break;
        for (int i = N / p; i >= 1; --i) f[(size_t)i * p] = inv_op(f[(size_t)i * p], f[(size_t)i]);
    }
}
// g[d] = Σ_{d | n, n ≤ N} f[n]（倍数の和）
template <class T, class Op = math_internal::DtAdd>
void multiple_zeta(std::vector<T>& f, Op op = Op{}) {
    const int N = (int)f.size() - 1;
    for (int p : math_internal::dt_primes(N)) {
        if (p > N) break;
        for (int i = N / p; i >= 1; --i) f[(size_t)i] = op(f[(size_t)i], f[(size_t)i * p]);
    }
}
// multiple_zeta の逆
template <class T, class Op = math_internal::DtSub>
void multiple_mobius(std::vector<T>& f, Op inv_op = Op{}) {
    const int N = (int)f.size() - 1;
    for (int p : math_internal::dt_primes(N)) {
        if (p > N) break;
        for (int i = 1; (long long)i * p <= N; ++i) f[(size_t)i] = inv_op(f[(size_t)i], f[(size_t)i * p]);
    }
}
// c[k] = Σ_{gcd(i, j) = k} a[i] b[j]（i, j, k = 1〜N。長さは長い方にそろえる）  O(N log log N)
template <class T>
std::vector<T> gcd_convolution(std::vector<T> a, std::vector<T> b) {
    const size_t n = std::max(a.size(), b.size());
    a.resize(n, T(0)), b.resize(n, T(0));
    multiple_zeta(a), multiple_zeta(b);
    for (size_t i = 1; i < n; ++i) a[i] = math_internal::dt_mul(a[i], b[i]);
    multiple_mobius(a);
    if (n) a[0] = T(0);
    return a;
}
// c[k] = Σ_{lcm(i, j) = k} a[i] b[j]（k ≤ N のものだけ）
template <class T>
std::vector<T> lcm_convolution(std::vector<T> a, std::vector<T> b) {
    const size_t n = std::max(a.size(), b.size());
    a.resize(n, T(0)), b.resize(n, T(0));
    divisor_zeta(a), divisor_zeta(b);
    for (size_t i = 1; i < n; ++i) a[i] = math_internal::dt_mul(a[i], b[i]);
    divisor_mobius(a);
    if (n) a[0] = T(0);
    return a;
}
// cnt[v] = 値 v（1〜N）の個数のとき、gcd がちょうど g になる組 {i < j} の数（g = 1〜N）  O(N log log N)
inline std::vector<long long> count_pairs_by_gcd(const std::vector<long long>& cnt) {
    std::vector<long long> m(cnt);
    if (!m.empty()) m[0] = 0;
    multiple_zeta(m);  // m[d] = d の倍数の個数
    for (size_t d = 1; d < m.size(); ++d) {  // C(m, 2)（偶数の方を先に 2 で割る。溢れたら止める）
        const long long x = m[d], a = x % 2 == 0 ? x / 2 : x, b = x % 2 == 0 ? x - 1 : (x - 1) / 2;
        if (__builtin_mul_overflow(a, b, &m[d])) math_internal::fail("count_pairs_by_gcd", "組の数が long long に収まらない");
    }
    multiple_mobius(m);
    return m;
}
// ディリクレ積 c[n] = Σ_{d | n} a[d] b[n/d]（n = 1〜N。長さは長い方にそろえる）  O(N log N)
template <class T>
std::vector<T> dirichlet_convolution(const std::vector<T>& a, const std::vector<T>& b) {
    const size_t n = std::max(a.size(), b.size());
    std::vector<T> c(n, T(0));
    for (size_t i = 1; i < a.size(); ++i)
        for (size_t j = 1; j < b.size() && i * j < n; ++j) c[i * j] = math_internal::dt_add(c[i * j], math_internal::dt_mul(a[i], b[j]));
    return c;
}
// ---- M の約数の上での変換（M ≤ 10^18。約数の数 D ≤ 103680）----
// DivisorsOf D(M)：D.divisors() は M の約数（昇順）。配列 f は約数の番号（D.index(d)）で引く（長さ D.size()）
struct DivisorsOf {
  private:
    long long M_;
    std::vector<PrimePower> fs_;
    std::vector<long long> ds_;
    std::vector<std::vector<std::pair<int, int>>> steps_;  // 素数ごとに (i, j)：ds_[j] = ds_[i] · p（i の昇順）
    template <class T>
    void check_size(const std::vector<T>& f, const char* func) const {
        if (f.size() != ds_.size()) math_internal::fail(func, "配列の長さを約数の数 size() にする");
    }

  public:
    explicit DivisorsOf(long long M) : M_(M) {
        math_internal::check_positive(M, "DivisorsOf");
        fs_ = ::factorize(M), ds_ = divisors_from_factors(fs_);
        for (const PrimePower& f : fs_) {
            std::vector<std::pair<int, int>> st;
            size_t j = 0;
            for (size_t i = 0; i < ds_.size(); ++i) {
                if ((M / ds_[i]) % f.p != 0) continue;  // ds_[i] · p が M を割り切るものだけ
                const long long t = ds_[i] * f.p;
                while (ds_[j] < t) ++j;
                st.push_back({(int)i, (int)j});
            }
            steps_.push_back(std::move(st));
        }
    }
    long long value() const { return M_; }
    int size() const { return (int)ds_.size(); }                          // 約数の数
    const std::vector<long long>& divisors() const { return ds_; }       // 約数（昇順）
    const std::vector<PrimePower>& factors() const { return fs_; }       // M の素因数分解
    long long operator[](int i) const { return ds_[(size_t)i]; }          // i 番目の約数
    int index(long long d) const {  // 約数 d の番号（d が M の約数でなければ -1）。二分探索（std::lower_bound はデバッグビルドで O(D) になるので自分で）
        size_t lo = 0, hi = ds_.size();
        while (lo < hi) {
            const size_t mid = (lo + hi) / 2;
            if (ds_[mid] < d) lo = mid + 1;
            else hi = mid;
        }
        return lo < ds_.size() && ds_[lo] == d ? (int)lo : -1;
    }
    // g[d] = Σ_{e | d} f[e]（d | M）  O(D · ω(M))
    template <class T, class Op = math_internal::DtAdd>
    void zeta(std::vector<T>& f, Op op = Op{}) const {
        check_size(f, "DivisorsOf::zeta");
        for (const auto& st : steps_)
            for (const auto& ij : st) f[(size_t)ij.second] = op(f[(size_t)ij.second], f[(size_t)ij.first]);
    }
    template <class T, class Op = math_internal::DtSub>
    void mobius(std::vector<T>& f, Op inv_op = Op{}) const {  // zeta の逆
        check_size(f, "DivisorsOf::mobius");
        for (const auto& st : steps_)
            for (auto it = st.rbegin(); it != st.rend(); ++it) f[(size_t)it->second] = inv_op(f[(size_t)it->second], f[(size_t)it->first]);
    }
    // g[d] = Σ_{d | e | M} f[e]（倍数の和）
    template <class T, class Op = math_internal::DtAdd>
    void multiple_zeta(std::vector<T>& f, Op op = Op{}) const {
        check_size(f, "DivisorsOf::multiple_zeta");
        for (const auto& st : steps_)
            for (auto it = st.rbegin(); it != st.rend(); ++it) f[(size_t)it->first] = op(f[(size_t)it->first], f[(size_t)it->second]);
    }
    template <class T, class Op = math_internal::DtSub>
    void multiple_mobius(std::vector<T>& f, Op inv_op = Op{}) const {  // multiple_zeta の逆
        check_size(f, "DivisorsOf::multiple_mobius");
        for (const auto& st : steps_)
            for (const auto& ij : st) f[(size_t)ij.first] = inv_op(f[(size_t)ij.first], f[(size_t)ij.second]);
    }
};
// Σ_{i=1}^{N} φ(i)・Σ μ(i)（N ≤ 10^11 程度。杜教篩：小さいところは表、大きいところは ⌊N/k⌋ ごとに）  O(N^{2/3})
namespace math_internal {
struct DtPrefix {  // ⌊N/i⌋ の形の v での和 S(v)：v ≤ T は表 small、それより大きい v = ⌊N/k⌋ は big[k]
    long long N = 0, T = 0;
    std::vector<long long> small;
    std::vector<Int128> big;
    Int128 operator()(long long v) const { return v <= T ? (Int128)small[(size_t)v] : big[(size_t)(N / v)]; }
};
template <class Base>
void dt_du_sieve(DtPrefix& P, Base base) {  // S(x) = base(x) − Σ_{d=2}^{x} S(⌊x/d⌋)
    const long long N = P.N, T = P.T;
    if (N <= T) return;
    const long long K = N / T;
    P.big.assign((size_t)K + 1, 0);
    for (long long k = K; k >= 1; --k) {
        const long long x = N / k;
        Int128 r = base(x);
        for (long long d = 2; d <= x;) {
            const long long q = x / d, d2 = x / q;
            const Int128 v = q <= T ? (Int128)P.small[(size_t)q] : P.big[(size_t)(k * d)];
            r -= (Int128)(d2 - d + 1) * v;
            d = d2 + 1;
        }
        P.big[(size_t)k] = r;
    }
}
inline long long dt_threshold(long long N) {  // 表を作る大きさ（N^{2/3}。10^7 まで）
    long long t = (long long)std::cbrt((long double)N);
    t = t * t;
    return std::max<long long>(1, std::min<long long>({N, t, 10000000LL}));
}
inline DtPrefix dt_phi_prefix(long long N, const char* func) {  // Σφ を ⌊N/i⌋ の形の値すべてで
    if (N > 10000000000000LL) fail(func, "N ≤ 10^13 で使う（表の大きさ。10^11 で約 1.5 秒）");
    DtPrefix P;
    P.N = std::max(N, 0LL), P.T = dt_threshold(std::max(N, 1LL));
    P.small = phi_table((int)P.T);
    for (size_t i = 1; i < P.small.size(); ++i) P.small[i] += P.small[i - 1];
    dt_du_sieve(P, [](long long x) { return (Int128)x * (x + 1) / 2; });
    return P;
}
}  // namespace math_internal
inline Int128 phi_sum(long long N) {
    if (N <= 0) return 0;
    return math_internal::dt_phi_prefix(N, "phi_sum")(N);
}
// Σ_{i=1}^{N} μ(i)（メルテンス関数。N ≤ 10^11 程度）  O(N^{2/3})
inline long long mertens(long long N) {
    if (N <= 0) return 0;
    if (N > 10000000000000LL) math_internal::fail("mertens", "N ≤ 10^13 で使う（表の大きさ。10^11 で約 1.5 秒）");
    math_internal::DtPrefix P;
    P.N = N, P.T = math_internal::dt_threshold(N);
    const std::vector<int> mu = mu_table((int)P.T);
    P.small.assign(mu.size(), 0);
    for (size_t i = 1; i < P.small.size(); ++i) P.small[i] = P.small[i - 1] + mu[i];
    math_internal::dt_du_sieve(P, [](long long) { return (Int128)1; });
    return (long long)P(N);
}
// 1 ≤ a, b ≤ N で gcd(a, b) = 1 となる組 (a, b) の数（順序あり）= 2 Σφ − 1
inline Int128 count_coprime_pairs(long long N) { return N <= 0 ? 0 : 2 * phi_sum(N) - 1; }
// 長さ len、値が 1〜K の列すべての gcd の和 mod m（ABC162 E）= Σ_{d ≤ K} φ(d) ⌊K/d⌋^len（⌊K/d⌋ が同じ d をまとめる。len = 0 なら 0）
//   K ≤ 10^13 程度  O(K^{2/3} + √K log len)
inline long long sum_gcd_tuples(long long K, long long len, long long m) {
    if (K < 0 || len < 0 || m < 1) math_internal::fail("sum_gcd_tuples", "K, len ≥ 0、m ≥ 1 で使う");
    if (K == 0 || len == 0) return 0;  // 列がない / 空の列の gcd は 0
    const math_internal::DtPrefix P = math_internal::dt_phi_prefix(K, "sum_gcd_tuples");
    long long s = 0;
    Int128 prev = 0;  // Φ(d − 1)
    for (long long d = 1; d <= K;) {
        const long long q = K / d, d2 = K / q;
        const Int128 cur = P(d2);
        const long long blk = (long long)((cur - prev) % m);
        s = (long long)(((Int128)s + (Int128)blk * mod_pow(q, len, m)) % m);
        prev = cur, d = d2 + 1;
    }
    return s;
}
#endif  // GG_DIVISOR_TRANSFORM
