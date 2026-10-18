#ifndef GG_SIEVE
#define GG_SIEVE
// ---- 篩（数学第二世代）  Sieve S(N)：N 以下の数を何度も素因数分解・素数判定、φ・μ などの表、区間篩 ----
// Sieve は N を超える数も受け付ける（02_factorize の is_prime / factorize に任せる）
namespace math_internal {
inline int checked_table_size(long long N, const char* func) {  // 0〜N の表を作る前に大きさを確かめる
    if (N < 0) fail(func, "N ≥ 0 で使う");
    if (N > 100000000) fail(func, "N が大きすぎる（10^8 まで）");
    return (int)N;
}
}  // namespace math_internal
struct Sieve {
  private:
    int n_;
    std::vector<int> spf_;
    std::vector<long long> primes_;
    static int checked_size(long long N) {  // 表を作る前に大きさを確かめる（N = 0 は 1 として作る）
        if (N < 0) math_internal::fail("Sieve", "N ≥ 0 で使う");
        if (N > 100000000) math_internal::fail("Sieve", "N が大きすぎる（10^8 まで。素数の一覧だけなら primes_up_to）");
        return (int)std::max(N, 1LL);
    }
    std::vector<PrimePower> small_factorize(long long x) const {
        std::vector<PrimePower> fs;
        while (x > 1) {
            const int p = spf_[x];
            int e = 0;
            while (x % p == 0) x /= p, ++e;
            fs.push_back({p, e});
        }
        return fs;
    }

  public:
    explicit Sieve(long long N) : n_(checked_size(N)), spf_((size_t)n_ + 1, 0) {  // 線形篩 O(N)。N ≤ 5×10^7 程度（メモリ 4N バイト）
        spf_[1] = 1;
        for (int i = 2; i <= n_; ++i) {
            if (spf_[i] == 0) spf_[i] = i, primes_.push_back(i);
            for (const long long p : primes_) {
                if (p > spf_[i] || p * i > n_) break;
                spf_[p * i] = (int)p;
            }
        }
    }
    int size() const { return n_; }
    const std::vector<long long>& primes() const { return primes_; }  // N 以下の素数（昇順。long long なので p * p も溢れない）
    int spf(long long x) const {                                        // 最小の素因数（2 ≤ x ≤ N。spf(1) = 1）
        if (x < 1 || x > n_) math_internal::fail("Sieve::spf", "1 ≤ x ≤ N で使う");
        return spf_[x];
    }
    bool is_prime(long long x) const { return x <= n_ ? (x >= 2 && spf_[x] == x) : ::is_prime(x); }
    std::vector<PrimePower> factorize(long long x) const {  // x ≥ 1（N を超えたら 02_factorize の factorize）
        math_internal::check_positive(x, "Sieve::factorize");
        return x <= n_ ? small_factorize(x) : ::factorize(x);
    }
    std::vector<long long> prime_factors(long long x) const {
        std::vector<long long> r;
        for (const PrimePower& f : factorize(x)) r.push_back(f.p);
        return r;
    }
    std::vector<long long> divisors(long long x) const { return divisors_from_factors(factorize(x)); }  // 昇順
    long long count_divisors(long long x) const {
        long long c = 1;
        for (const PrimePower& f : factorize(x)) c *= f.e + 1;
        return c;
    }
    Int128 sum_divisors(long long x) const {
        Int128 s = 1;
        for (const PrimePower& f : factorize(x)) {
            Int128 t = 1, pk = 1;
            for (int k = 0; k < f.e; ++k) pk *= f.p, t += pk;
            s *= t;
        }
        return s;
    }
    long long phi(long long x) const {
        long long r = x;
        for (const PrimePower& f : factorize(x)) r = r / f.p * (f.p - 1);
        return r;
    }
    int mu(long long x) const {
        const std::vector<PrimePower> fs = factorize(x);
        for (const PrimePower& f : fs)
            if (f.e >= 2) return 0;
        return fs.size() % 2 ? -1 : 1;
    }
    int omega(long long x) const { return (int)factorize(x).size(); }  // 異なる素因数の数
    int big_omega(long long x) const {                                  // 素因数の数（重複も数える）
        int c = 0;
        for (const PrimePower& f : factorize(x)) c += f.e;
        return c;
    }
    long long largest_prime_factor(long long x) const {  // x = 1 なら 1
        const std::vector<PrimePower> fs = factorize(x);
        return fs.empty() ? 1 : fs.back().p;
    }
    SquareFreeParts square_free_decompose(long long x) const {  // x = a^2 · b
        long long a = 1, b = 1;
        for (const PrimePower& f : factorize(x)) {
            for (int k = 0; k < f.e / 2; ++k) a *= f.p;
            if (f.e % 2) b *= f.p;
        }
        return {a, b};
    }
    bool is_square_free(long long x) const { return square_free_decompose(x).a == 1; }
    long long prime_count(long long x) const {  // x 以下の素数の数（x ≤ N）
        if (x > n_) math_internal::fail("Sieve::prime_count", "x ≤ N で使う（大きい x は 02_prime_count）");
        int lo = 0, hi = (int)primes_.size();  // primes_[lo..hi) で x を超える最初の位置を二分探索
        while (lo < hi) {
            const int mid = (lo + hi) / 2;
            if (primes_[mid] <= x) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
    long long count_primes(long long l, long long r) const {  // [l, r) の素数の数（r ≤ N + 1）
        if (l >= r) return 0;
        return prime_count(r - 1) - (l >= 1 ? prime_count(l - 1) : 0);
    }
    long long count_primes_closed(long long L, long long R) const {  // [L, R] の素数の数（R ≤ N）
        if (L > R) return 0;
        if (R > n_) math_internal::fail("Sieve::count_primes_closed", "R ≤ N で使う（大きい R は 02_prime_count）");
        return count_primes(L, R + 1);
    }
    long long kth_prime(long long k) const {  // 小さい方から k 番目の素数（0 番目が 2。N 以下の素数の中で）
        if (k < 0 || k >= (long long)primes_.size()) math_internal::fail("Sieve::kth_prime", "k が範囲の外（N を大きくする）");
        return primes_[k];
    }
    template <class V>
    long long lcm_mod(const V& vs, long long m) const {  // 全部の数（≥ 1）の lcm を mod m で（vector<int> なども可）
        if (m < 1) math_internal::fail("Sieve::lcm_mod", "m ≥ 1 で使う");
        std::vector<long long> xs;
        for (const auto& xv : vs) {  // 先に全部を確かめる
            const long long x = math_internal::to_ll_checked(xv, "Sieve::lcm_mod");
            if (x < 0) math_internal::fail("Sieve::lcm_mod", "数は 0 以上で使う");
            xs.push_back(x);
        }
        std::vector<PrimePower> all;
        for (const long long x : xs) {
            if (x == 0) return 0;
            for (const PrimePower& f : factorize(x)) all.push_back(f);
        }
        std::sort(all.begin(), all.end(), [](const PrimePower& a, const PrimePower& b) { return a.p != b.p ? a.p < b.p : a.e > b.e; });
        long long r = 1 % m;
        for (size_t i = 0; i < all.size(); ++i)
            if (i == 0 || all[i].p != all[i - 1].p) r = mod_mul(r, mod_pow(all[i].p, all[i].e, m), m);
        return r;
    }
    long long lcm_mod(std::initializer_list<long long> vs, long long m) const { return lcm_mod<std::initializer_list<long long>>(vs, m); }
};

// ---- 0〜N の表（O(N)。0 ≤ N ≤ 10^8）----
inline std::vector<long long> phi_table(long long N_) {  // phi[x] = φ(x)（phi[0] = 0）
    const int N = math_internal::checked_table_size(N_, "phi_table");
    std::vector<long long> phi(N + 1);
    std::vector<int> ps;
    for (int i = 0; i <= N; ++i) phi[i] = i;
    std::vector<char> comp(N + 1, 0);
    for (int i = 2; i <= N; ++i) {
        if (!comp[i]) ps.push_back(i), phi[i] = i - 1;
        for (int p : ps) {
            if ((long long)p * i > N) break;
            comp[p * i] = 1;
            if (i % p == 0) {
                phi[p * i] = phi[i] * p;
                break;
            }
            phi[p * i] = phi[i] * (p - 1);
        }
    }
    return phi;
}
inline std::vector<int> mu_table(long long N_) {  // mu[x] = μ(x)（mu[0] = 0）
    const int N = math_internal::checked_table_size(N_, "mu_table");
    std::vector<int> mu(N + 1, 0), ps;
    std::vector<char> comp(N + 1, 0);
    if (N >= 1) mu[1] = 1;
    for (int i = 2; i <= N; ++i) {
        if (!comp[i]) ps.push_back(i), mu[i] = -1;
        for (int p : ps) {
            if ((long long)p * i > N) break;
            comp[p * i] = 1;
            if (i % p == 0) {
                mu[p * i] = 0;
                break;
            }
            mu[p * i] = -mu[i];
        }
    }
    return mu;
}
// 以下も線形篩（x = p·i（p は x の最小の素因数）を i から決める。e[x] は最小の素因数の指数）
inline std::vector<int> divisor_count_table(long long N_) {  // d[x] = x の約数の個数（d[0] = 0）  O(N)
    const int N = math_internal::checked_table_size(N_, "divisor_count_table");
    std::vector<int> d(N + 1, 0), ps;
    std::vector<unsigned char> e(N + 1, 0);
    if (N >= 1) d[1] = 1;
    for (int i = 2; i <= N; ++i) {
        if (d[i] == 0) ps.push_back(i), d[i] = 2, e[i] = 1;  // まだ決まっていない = 素数
        for (int p : ps) {
            if ((long long)p * i > N) break;
            const int j = p * i;
            if (i % p == 0) {
                e[j] = (unsigned char)(e[i] + 1), d[j] = d[i] / (e[i] + 1) * (e[i] + 2);
                break;
            }
            e[j] = 1, d[j] = d[i] * 2;
        }
    }
    return d;
}
inline std::vector<long long> divisor_sum_table(long long N_) {  // s[x] = x の約数の和（s[0] = 0）  O(N)
    const int N = math_internal::checked_table_size(N_, "divisor_sum_table");
    std::vector<long long> s(N + 1, 0);
    std::vector<unsigned char> e(N + 1, 0);
    std::vector<int> ps;
    if (N >= 1) s[1] = 1;
    for (int i = 2; i <= N; ++i) {
        if (s[i] == 0) ps.push_back(i), s[i] = i + 1, e[i] = 1;
        for (int p : ps) {
            if ((long long)p * i > N) break;
            const int j = p * i;
            if (i % p == 0) {  // i = p^e · rest → j = p^{e+1} · rest
                long long pk = 1, term = 1, pw = 1;
                for (int t = 0; t < e[i]; ++t) pk *= p;
                for (int t = 0; t <= e[i]; ++t) pw *= p, term += pw;
                e[j] = (unsigned char)(e[i] + 1), s[j] = s[i / pk] * term;
                break;
            }
            e[j] = 1, s[j] = s[i] * (p + 1);
        }
    }
    return s;
}
inline std::vector<int> omega_table(long long N_) {  // w[x] = x の異なる素因数の数（w[0] = w[1] = 0）  O(N)
    const int N = math_internal::checked_table_size(N_, "omega_table");
    std::vector<int> w(N + 1, 0), ps;
    std::vector<char> comp(N + 1, 0);
    for (int i = 2; i <= N; ++i) {
        if (!comp[i]) ps.push_back(i), w[i] = 1;
        for (int p : ps) {
            if ((long long)p * i > N) break;
            const int j = p * i;
            comp[j] = 1;
            if (i % p == 0) {
                w[j] = w[i];
                break;
            }
            w[j] = w[i] + 1;
        }
    }
    return w;
}
// N 以下の素数（エラトステネス。奇数だけ・区間ごと。N < 2^31。10^9 なら 5×10^7 個・約 3 秒）
// 返り値は vector<long long>（p * p も溢れない）。メモリを半分にするなら primes_up_to<int>(N)（N = 10^9 で 400MB → 200MB）
template <class T = long long>
std::vector<T> primes_up_to(long long N) {
    if (N > 2147483647LL) math_internal::fail("primes_up_to", "N < 2^31 で使う");
    std::vector<T> ps;
    if (N < 2) return ps;
    ps.push_back(2);
    const int n = (int)N, sq = (int)isqrt(N);
    std::vector<char> small(sq + 1, 1);
    std::vector<int> base;  // √N 以下の奇素数
    for (int i = 3; i <= sq; i += 2)
        if (small[i]) {
            base.push_back(i);
            for (long long j = (long long)i * i; j <= sq; j += 2 * i) small[j] = 0;
        }
    const int B = 1 << 16;  // 奇数 B 個ずつ
    std::vector<char> blk(B);
    std::vector<long long> nxt(base.size());  // 次に消す奇数の倍数
    for (size_t i = 0; i < base.size(); ++i) nxt[i] = (long long)base[i] * base[i];
    for (long long lo = 3; lo <= n; lo += 2LL * B) {  // [lo, lo + 2B) の奇数
        const long long hi = std::min<long long>(lo + 2LL * B, (long long)n + 1);
        std::fill(blk.begin(), blk.end(), 1);
        for (size_t i = 0; i < base.size(); ++i) {
            long long j = nxt[i];
            const long long step = 2LL * base[i];
            for (; j < hi; j += step) blk[(j - lo) >> 1] = 0;
            nxt[i] = j;
        }
        for (long long x = lo; x < hi; x += 2)
            if (blk[(x - lo) >> 1]) ps.push_back((T)x);
    }
    return ps;
}

// ---- 区間篩（[l, r) の数。r ≤ 10^14 程度、r − l ≤ 10^7 程度。√r 以下の素数で割る）----
inline std::vector<long long> segment_primes(long long l, long long r) {  // [l, r) の素数（昇順）
    std::vector<long long> res;
    l = std::max(l, 2LL);
    if (l >= r) return res;
    const long long sq = isqrt(r - 1);
    std::vector<char> small(sq + 1, 1), seg((size_t)(r - l), 1);
    for (long long p = 2; p <= sq; ++p) {
        if (!small[p]) continue;
        for (long long j = p * p; j <= sq; j += p) small[j] = 0;
        for (long long j = std::max(p * p, ceil_multiple(l, p)); j < r; j += p) seg[(size_t)(j - l)] = 0;
    }
    for (long long x = l; x < r; ++x)
        if (seg[(size_t)(x - l)]) res.push_back(x);
    return res;
}
inline std::vector<long long> segment_primes_closed(long long L, long long R) { return segment_primes(L, R + 1); }
// [l, r) の各数 x（≥ 1）の素因数分解を、x の小さい順に f(x, fs) へ渡す（fs は p の昇順の vector<PrimePower>）
template <class F>
void segment_factorize(long long l, long long r, F f) {
    l = std::max(l, 1LL);
    if (l >= r) return;
    const long long sq = isqrt(r - 1);
    std::vector<char> small(sq + 1, 1);
    std::vector<long long> base;
    for (long long p = 2; p <= sq; ++p) {
        if (!small[p]) continue;
        base.push_back(p);
        for (long long j = p * p; j <= sq; j += p) small[j] = 0;
    }
    const long long B = 1 << 18;
    std::vector<long long> rem, nxt(base.size());
    for (size_t k = 0; k < base.size(); ++k) nxt[k] = ceil_multiple(l, base[k]);  // 次に割る倍数
    std::vector<std::vector<PrimePower>> fs((size_t)std::min(B, r - l));
    for (long long lo = l; lo < r; lo += B) {
        const long long hi = std::min(lo + B, r);
        rem.resize((size_t)(hi - lo));
        for (long long x = lo; x < hi; ++x) rem[(size_t)(x - lo)] = x, fs[(size_t)(x - lo)].clear();
        for (size_t k = 0; k < base.size(); ++k) {
            const long long p = base[k];
            long long j = nxt[k];
            for (; j < hi; j += p) {
                long long& v = rem[(size_t)(j - lo)];
                int e = 0;
                while (v % p == 0) v /= p, ++e;
                fs[(size_t)(j - lo)].push_back({p, e});
            }
            nxt[k] = j;
        }
        for (long long x = lo; x < hi; ++x) {
            const size_t i = (size_t)(x - lo);
            if (rem[i] > 1) fs[i].push_back({rem[i], 1});
            f(x, (const std::vector<PrimePower>&)fs[i]);
        }
    }
}
template <class F>
void segment_factorize_closed(long long L, long long R, F f) { segment_factorize(L, R + 1, f); }
#endif  // GG_SIEVE
