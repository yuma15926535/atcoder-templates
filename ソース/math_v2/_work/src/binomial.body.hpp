#ifndef GG_BINOMIAL
#define GG_BINOMIAL
// ---- 二項係数（数学第二世代）  Binomial<mint>（階乗の表）・n が大きい nCr・Lucas・合成数の mod・厳密な値・パスカルの三角形・Σ_{i≤m} C(n, i)（Mo 用） ----
// Binomial<mint> B(N)：0〜N の階乗・階乗の逆元の表（mint の mod は素数。N を超える n が来たら自動で表を伸ばす）
template <class mint>
struct Binomial {
    explicit Binomial(long long n = 0) { ensure(std::max(n, 1LL)); }
    mint fact(long long n) const {  // n!
        if (n < 0) math_internal::fail("Binomial::fact", "n ≥ 0 で使う");
        ensure(n);
        return f_[(size_t)n];
    }
    mint inv_fact(long long n) const {  // 1 / n!
        if (n < 0) math_internal::fail("Binomial::inv_fact", "n ≥ 0 で使う");
        ensure(n);
        return fi_[(size_t)n];
    }
    mint inv(long long n) const {  // 1 / n（n ≥ 1）
        if (n < 1) math_internal::fail("Binomial::inv", "n ≥ 1 で使う");
        ensure(n);
        return fi_[(size_t)n] * f_[(size_t)n - 1];
    }
    mint C(long long n, long long r) const {  // nCr（r < 0・r > n・n < 0 なら 0）
        if (r < 0 || r > n) return mint(0);
        ensure(n);
        return f_[(size_t)n] * fi_[(size_t)r] * fi_[(size_t)(n - r)];
    }
    mint P(long long n, long long r) const {  // nPr = n! / (n − r)!（範囲外は 0）
        if (r < 0 || r > n) return mint(0);
        ensure(n);
        return f_[(size_t)n] * fi_[(size_t)(n - r)];
    }
    mint H(long long n, long long r) const {  // nHr = (n + r − 1)Cr（n 種類から重複を許して r 個。H(n, 0) = 1）
        if (n < 0 || r < 0) return mint(0);
        if (r == 0) return mint(1);
        return C(n + r - 1, r);
    }
    template <class V>
    mint multinomial(const V& ks) const {  // (Σ k)! / Π k!（負の k があれば 0）
        long long s = 0;
        for (auto kv : ks) {
            const long long k = math_internal::to_ll_checked(kv, "Binomial::multinomial");
            if (k < 0) return mint(0);
            if (__builtin_add_overflow(s, k, &s)) math_internal::fail("Binomial::multinomial", "個数の和が long long に収まらない");
        }
        ensure(s);
        mint r = f_[(size_t)s];
        for (auto k : ks) r *= fi_[(size_t)k];
        return r;
    }
    mint multinomial(std::initializer_list<long long> ks) const { return multinomial<std::initializer_list<long long>>(ks); }
    mint catalan(long long n) const {  // カタラン数 C(2n, n) / (n + 1)
        if (n < 0) return mint(0);
        return C(2 * n, n) - C(2 * n, n + 1);
    }
    mint paths(long long x, long long y) const { return x < 0 || y < 0 ? mint(0) : C(x + y, x); }  // (0, 0) → (x, y) の右・上だけの道
    // (0, 0) → (x, y) の右・上だけの道で、ずっと y ≤ x + k のもの（k ≥ 0。鏡像法：C(x + y, x) − C(x + y, y − k − 1)）
    mint paths_not_crossing(long long x, long long y, long long k) const {
        if (x < 0 || y < 0 || k < 0 || y > x + k) return mint(0);
        return C(x + y, x) - C(x + y, y - k - 1);
    }

  private:
    mutable std::vector<mint> f_, fi_;
    mutable int built_mod_ = 0;
    void ensure(long long n) const {
        if (n >= (long long)f_.size() || built_mod_ != mint::mod()) grow(n);
    }
    void grow(long long n) const {
        const long long M = mint::mod();
        if (built_mod_ != M) {  // はじめて、または DynModInt の mod が変わった
            if (!math_internal::is_prime_u32((unsigned)M)) math_internal::fail("Binomial", "mod が素数でない（合成数の mod は pascal_table か binom_mod_any）");
            f_.clear(), fi_.clear(), built_mod_ = (int)M;
        }
        if (n >= M) math_internal::fail("Binomial", "n ≥ mod では n! が 0 になる（binom_lucas を使う）");
        if (n > 100000000) math_internal::fail("Binomial", "n が大きすぎる（表は 10^8 まで。r が小さいなら binom_large_n）");
        const size_t old = f_.size();
        if ((long long)old > n) return;
        size_t sz = std::max<size_t>((size_t)n + 1, std::min<size_t>(old + old / 2, (size_t)M));  // 1.5 倍ずつ伸ばす
        sz = std::min<size_t>(sz, (size_t)100000001);
        f_.resize(sz), fi_.resize(sz);
        if (old == 0) f_[0] = 1;
        for (size_t i = std::max<size_t>(old, 1); i < sz; ++i) f_[i] = f_[i - 1] * mint((long long)i);
        fi_[sz - 1] = f_[sz - 1].inv();
        for (size_t i = sz - 1; i > old; --i) fi_[i - 1] = fi_[i] * mint((long long)i);
    }
};

// S(n, m) = Σ_{i=0}^{m} C(n, i)（n 個から m 個以下を選ぶ数）を、n・m を 1 ずつ動かしながら O(1) で保つ（Mo のアルゴリズムで多くのクエリに）
//   n ≥ 0。m < 0 なら 0、m ≥ n なら 2^n。mod は 3 以上の素数（n を減らすとき 2 で割る）
template <class mint>
struct BinomialPrefix {
    explicit BinomialPrefix(long long max_n = 0) : B_(std::max(max_n, 1LL)) {
        if (mint::mod() == 2) math_internal::fail("BinomialPrefix", "mod = 2 では 2 で割れない");
    }
    long long n() const { return n_; }
    long long m() const { return m_; }
    mint value() const { return s_; }  // Σ_{i=0}^{m} C(n, i)
    void inc_n() { s_ = s_ + s_ - B_.C(n_, m_), ++n_; }  // S(n+1, m) = 2 S(n, m) − C(n, m)
    void dec_n() {                                       // S(n−1, m) = (S(n, m) + C(n−1, m)) / 2
        if (n_ == 0) math_internal::fail("BinomialPrefix::dec_n", "n ≥ 0 のまま使う");
        --n_, s_ = (s_ + B_.C(n_, m_)) * inv2();
    }
    void inc_m() { ++m_, s_ += B_.C(n_, m_); }
    void dec_m() { s_ -= B_.C(n_, m_), --m_; }
    mint move_to(long long n, long long m) {  // (n, m) まで 1 ずつ動かして値を返す（|Δn| + |Δm| 回。m ≥ n なら値は 2^n なので、m は min(m, n) にして渡すと速い）
        if (n < 0) math_internal::fail("BinomialPrefix::move_to", "n ≥ 0 で使う");
        while (n_ < n) inc_n();
        while (m_ > m) dec_m();
        while (n_ > n) dec_n();
        while (m_ < m) inc_m();
        return s_;
    }

  private:
    Binomial<mint> B_;
    long long n_ = 0, m_ = 0;
    mint s_ = mint(1);  // S(0, 0) = 1
    mint inv2() const {
        static long long cached_mod = -1;
        static mint v;
        if (cached_mod != (long long)mint::mod()) cached_mod = mint::mod(), v = mint(2).inv();
        return v;
    }
};

// n が大きく r が小さい nCr（n ≤ 10^18、r ≤ 10^7 程度、r < mod。mod は素数）  O(r)
template <class mint>
mint binom_large_n(long long n, long long r) {
    if (r < 0 || r > n) return mint(0);
    if (r >= mint::mod()) math_internal::fail("binom_large_n", "r ≥ mod（binom_lucas を使う）");
    mint num = 1, den = 1;
    for (long long i = 0; i < r; ++i) num *= mint(n - i), den *= mint(i + 1);
    return num / den;
}

// nCr mod p（p は素数で 10^7 以下。n, r ≤ 10^18。Lucas の定理）  前処理 O(p)（同じ p なら 1 回）・O(log_p n)
inline long long binom_lucas(long long n, long long r, long long p) {
    if (p < 2 || p > 10000000 || !math_internal::is_prime_u32((unsigned)p)) math_internal::fail("binom_lucas", "p は 10^7 以下の素数");
    if (r < 0 || r > n) return 0;
    static long long cached = 0;
    static std::vector<int> f, fi;
    if (cached != p) {
        f.assign((size_t)p, 1), fi.assign((size_t)p, 1);
        for (long long i = 1; i < p; ++i) f[(size_t)i] = (int)(f[(size_t)i - 1] * i % p);
        fi[(size_t)p - 1] = (int)mod_pow(f[(size_t)p - 1], p - 2, p);
        for (long long i = p - 1; i > 0; --i) fi[(size_t)i - 1] = (int)(fi[(size_t)i] * i % p);
        cached = p;
    }
    long long res = 1;
    while (r > 0) {
        const long long ni = n % p, ri = r % p;
        if (ri > ni) return 0;
        res = res * f[(size_t)ni] % p * fi[(size_t)ri] % p * fi[(size_t)(ni - ri)] % p;
        n /= p, r /= p;
    }
    return res;
}

namespace math_internal {
struct BinomPrimePower {  // nCr mod p^e のための表（p で割れる数を除いた積）
    long long p, pe;
    int e;
    std::vector<int> prod;  // prod[i] = Π_{j ≤ i, p ∤ j} j mod p^e
    long long fact_wo_p(long long n) const {  // n! から p をすべて除いたもの mod p^e
        long long r = 1;
        while (n > 0) {
            r = r * mod_pow(prod[(size_t)pe], n / pe, pe) % pe * prod[(size_t)(n % pe)] % pe;
            n /= p;
        }
        return r;
    }
    long long legendre(long long n) const {  // n! の中の p の個数
        long long c = 0;
        while (n > 0) n /= p, c += n;
        return c;
    }
    long long binom(long long n, long long r) const {
        const long long k = legendre(n) - legendre(r) - legendre(n - r);
        if (k >= e) return 0;
        long long v = fact_wo_p(n) * mod_inv(fact_wo_p(r) * fact_wo_p(n - r) % pe, pe) % pe;
        for (long long i = 0; i < k; ++i) v = v * p % pe;
        return v;
    }
};
}  // namespace math_internal
// nCr mod m（m は合成数でもよい。m の素数べきの部分はそれぞれ 10^7 以下。n, r ≤ 10^18）  前処理 O(Σ p^e)（同じ m なら 1 回）
inline long long binom_mod_any(long long n, long long r, long long m) {
    if (m < 1) math_internal::fail("binom_mod_any", "m ≥ 1 で使う");
    if (r < 0 || r > n) return 0;
    if (m == 1) return 0;
    static long long cached = 0;
    static std::vector<math_internal::BinomPrimePower> parts;
    if (cached != m) {
        parts.clear();
        long long x = m;
        for (long long p = 2; p <= 10000000 && p * p <= x; ++p) {  // 10^7 より大きい素因数は下で止める
            if (x % p) continue;
            math_internal::BinomPrimePower b{p, 1, 0, {}};
            while (x % p == 0) x /= p, b.pe *= p, b.e++;
            parts.push_back(b);
        }
        if (x > 1) parts.push_back(math_internal::BinomPrimePower{x, x, 1, {}});
        for (auto& b : parts) {
            if (b.pe > 10000000) math_internal::fail("binom_mod_any", "m の素数べきの部分が 10^7 を超える");
            b.prod.assign((size_t)b.pe + 1, 1);
            for (long long i = 1; i <= b.pe; ++i) b.prod[(size_t)i] = (int)(i % b.p == 0 ? b.prod[(size_t)i - 1] : b.prod[(size_t)i - 1] * i % b.pe);
        }
        cached = m;
    }
    long long res = 0, M = 1;  // 中国剰余定理で合わせる（素数べきは互いに素）
    for (const auto& b : parts) {
        const long long a = b.binom(n, r);
        const long long t = mod_mul(a - res % b.pe, mod_inv(M % b.pe, b.pe), b.pe);
        res += M * t, M *= b.pe;
    }
    return res;
}

// 厳密な nCr（long long に収まらなければ止める）/ min(nCr, cap)（溢れない）  O(min(r, n − r))（値が大きくなる前に終わる）
namespace math_internal {
inline Int128 binom_sat(long long n, long long r, Int128 cap) {  // min(nCr, cap)（0 ≤ r ≤ n、cap ≤ 2^63）
    r = std::min(r, n - r);
    Int128 v = 1;
    for (long long i = 1; i <= r; ++i) {  // v = C(n − r + i, i)。i が増えると増えるので、cap を超えたら終わり
        v = v * (n - r + i) / i;
        if (v > cap) return cap;
    }
    return std::min(v, cap);  // r = 0 で cap = 0 のとき
}
}  // namespace math_internal
inline long long binom_exact(long long n, long long r) {
    if (r < 0 || r > n) return 0;
    const Int128 L = (Int128)9223372036854775807LL;
    const Int128 v = math_internal::binom_sat(n, r, L + 1);
    if (v > L) math_internal::fail("binom_exact", "nCr が long long に収まらない（binom_cap か mod で）");
    return (long long)v;
}
inline long long binom_cap(long long n, long long r, long long cap) {
    if (cap < 0) math_internal::fail("binom_cap", "cap ≥ 0 で使う");
    if (r < 0 || r > n) return 0;
    return (long long)math_internal::binom_sat(n, r, cap);
}

// パスカルの三角形 t[n][r] = nCr（0 ≤ r ≤ n ≤ N）。T は何でも（mint・DynModInt・double・long long は溢れたら止める）  O(N^2)
template <class T = long long>
std::vector<std::vector<T>> pascal_table(int N) {
    if (N < 0) math_internal::fail("pascal_table", "N ≥ 0 で使う");
    std::vector<std::vector<T>> t((size_t)N + 1);
    for (int n = 0; n <= N; ++n) {
        t[(size_t)n].assign((size_t)n + 1, T(1));
        for (int r = 1; r < n; ++r) {
            if constexpr (std::is_integral_v<T>) {
                if (__builtin_add_overflow(t[(size_t)n - 1][(size_t)r - 1], t[(size_t)n - 1][(size_t)r], &t[(size_t)n][(size_t)r]))
                    math_internal::fail("pascal_table", "値が型に収まらない（long long なら N ≤ 66）");
            } else {
                t[(size_t)n][(size_t)r] = t[(size_t)n - 1][(size_t)r - 1] + t[(size_t)n - 1][(size_t)r];
            }
        }
    }
    return t;
}
#endif  // GG_BINOMIAL
