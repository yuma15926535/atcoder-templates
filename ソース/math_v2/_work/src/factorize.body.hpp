#ifndef GG_FACTORIZE
#define GG_FACTORIZE
// ---- 素因数分解（数学第二世代）  1 つの大きな数（n < 2^63）の素数判定・素因数分解・約数・φ・μ など ----
struct PrimePower {  // 素因数分解の 1 項 p^e（素因数分解は p の昇順の vector<PrimePower>）
    long long p;
    int e;
    friend bool operator==(const PrimePower& a, const PrimePower& b) { return a.p == b.p && a.e == b.e; }
    friend bool operator!=(const PrimePower& a, const PrimePower& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const PrimePower& f) { return f.e == 1 ? os << f.p : os << f.p << "^" << f.e; }
};
struct SquareFreeParts {  // n = a^2 · b（b は平方因子を持たない）
    long long a, b;
};

namespace math_internal {
struct Mont64 {  // 64 ビットのモンゴメリ乗算（n は奇数。値は [0, n) で持つ）
    unsigned long long n, ninv, r2;
    explicit Mont64(unsigned long long n_) : n(n_), ninv(n_), r2(0) {
        for (int i = 0; i < 6; ++i) ninv *= 2 - n * ninv;  // n^{-1} mod 2^64
        r2 = (unsigned long long)(((UInt128)0 - (UInt128)n) % n);  // 2^128 mod n
    }
    unsigned long long reduce(UInt128 t) const {  // t < n·2^64 → t·2^{-64} mod n
        const unsigned long long q = (unsigned long long)t * ninv;
        const unsigned long long hi = (unsigned long long)(t >> 64), m = (unsigned long long)(((UInt128)q * n) >> 64);
        return hi >= m ? hi - m : hi - m + n;
    }
    unsigned long long mul(unsigned long long a, unsigned long long b) const { return reduce((UInt128)a * b); }
    unsigned long long to(unsigned long long a) const { return mul(a % n, r2); }
    unsigned long long from(unsigned long long a) const { return reduce(a); }
    unsigned long long pow(unsigned long long a, unsigned long long e) const {  // a はモンゴメリ表現
        unsigned long long r = to(1);
        for (; e; e >>= 1, a = mul(a, a))
            if (e & 1) r = mul(r, a);
        return r;
    }
};
inline bool miller_rabin(unsigned long long n) {  // n は 5 以上の奇数
    unsigned long long d = n - 1;
    int s = 0;
    while (!(d & 1)) d >>= 1, ++s;
    if (n < (1ULL << 32)) {
        for (unsigned long long a : {2ULL, 7ULL, 61ULL}) {
            if (a % n == 0) continue;
            unsigned long long x = 1, b = a % n;
            for (unsigned long long e = d; e; e >>= 1, b = b * b % n)
                if (e & 1) x = x * b % n;
            if (x == 1 || x == n - 1) continue;
            bool comp = true;
            for (int i = 1; i < s && comp; ++i) {
                x = x * x % n;
                if (x == n - 1) comp = false;
            }
            if (comp) return false;
        }
        return true;
    }
    const Mont64 mt(n);
    const unsigned long long one = mt.to(1), mone = mt.to(n - 1);
    for (unsigned long long a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (a % n == 0) continue;
        unsigned long long x = mt.pow(mt.to(a), d);
        if (x == one || x == mone) continue;
        bool comp = true;
        for (int i = 1; i < s && comp; ++i) {
            x = mt.mul(x, x);
            if (x == mone) comp = false;
        }
        if (comp) return false;
    }
    return true;
}
inline unsigned long long find_factor(unsigned long long n) {  // 合成数 n の 1 でも n でもない約数（ポラード・ロー、Brent 版）
    if (n % 2 == 0) return 2;
    const Mont64 mt(n);
    for (unsigned long long c0 = 1;; ++c0) {
        const unsigned long long c = mt.to(c0);
        auto f = [&](unsigned long long x) {
            const unsigned long long y = mt.mul(x, x) + c;
            return y >= n ? y - n : y;
        };
        unsigned long long x = 0, y = mt.to(2), ys = 0, q = mt.to(1), g = 1;
        const unsigned long long M = 128;
        for (unsigned long long r = 1; g == 1; r <<= 1) {
            x = y;
            for (unsigned long long i = 0; i < r; ++i) y = f(y);
            for (unsigned long long k = 0; k < r && g == 1; k += M) {
                ys = y;
                for (unsigned long long i = 0; i < M && i < r - k; ++i) y = f(y), q = mt.mul(q, x > y ? x - y : y - x);
                g = ugcd(q, n);
            }
        }
        if (g == n) {
            do {
                ys = f(ys);
                g = ugcd(x > ys ? x - ys : ys - x, n);
            } while (g == 1);
        }
        if (g != n) return g;
    }
}
inline std::vector<PrimePower> group_primes(std::vector<long long>& ps) {
    std::sort(ps.begin(), ps.end());
    std::vector<PrimePower> fs;
    for (long long p : ps) {
        if (!fs.empty() && fs.back().p == p) ++fs.back().e;
        else fs.push_back({p, 1});
    }
    return fs;
}
inline void check_positive(long long n, const char* func) {
    if (n < 1) fail(func, "n ≥ 1 で使う");
}
}  // namespace math_internal

inline bool is_prime(long long n) {  // 素数か（n < 2^63。決定的ミラー・ラビン）
    if (n < 2) return false;
    for (long long p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        if (n % p == 0) return n == p;
    }
    if (n < 37 * 37) return true;
    return math_internal::miller_rabin((unsigned long long)n);
}
inline std::vector<PrimePower> factorize(long long n) {  // 素因数分解（n ≥ 1。p の昇順。n = 1 なら空。ポラード・ロー）
    math_internal::check_positive(n, "factorize");
    std::vector<long long> ps;
    for (long long p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47})
        while (n % p == 0) ps.push_back(p), n /= p;
    std::vector<unsigned long long> st;
    if (n > 1) st.push_back((unsigned long long)n);
    while (!st.empty()) {
        const unsigned long long m = st.back();
        st.pop_back();
        if (m < 53 * 53 || is_prime((long long)m)) {
            ps.push_back((long long)m);
            continue;
        }
        const unsigned long long d = math_internal::find_factor(m);
        st.push_back(d), st.push_back(m / d);
    }
    return math_internal::group_primes(ps);
}
inline std::vector<PrimePower> factorize_trial(long long n) {  // 試し割りの素因数分解（n ≥ 1。O(√n)。確かめ用）
    math_internal::check_positive(n, "factorize_trial");
    std::vector<PrimePower> fs;
    for (long long p = 2; p <= n / p; p += (p == 2 ? 1 : 2)) {
        if (n % p) continue;
        int e = 0;
        while (n % p == 0) n /= p, ++e;
        fs.push_back({p, e});
    }
    if (n > 1) fs.push_back({n, 1});
    return fs;
}
inline std::vector<long long> prime_factors(long long n) {  // 異なる素因数（昇順）
    std::vector<long long> r;
    for (const PrimePower& f : factorize(n)) r.push_back(f.p);
    return r;
}
template <class F>
void for_each_divisor(const std::vector<PrimePower>& fs, F f) {  // 素因数分解から全部の約数を f(d) に渡す（順番は決まっていない）
    std::vector<long long> ds = {1};
    for (const PrimePower& pe : fs) {
        const size_t sz = ds.size();
        long long pk = 1;
        for (int k = 1; k <= pe.e; ++k) {
            pk *= pe.p;
            for (size_t i = 0; i < sz; ++i) ds.push_back(ds[i] * pk);
        }
    }
    for (long long d : ds) f(d);
}
inline std::vector<long long> divisors_from_factors(const std::vector<PrimePower>& fs) {  // 素因数分解から約数（昇順）
    std::vector<long long> ds;
    for_each_divisor(fs, [&](long long d) { ds.push_back(d); });
    std::sort(ds.begin(), ds.end());
    return ds;
}
inline std::vector<long long> divisors_of(long long n) { return divisors_from_factors(factorize(n)); }  // 約数（昇順。n ≥ 1）
inline long long count_divisors(long long n) {  // 約数の個数
    long long c = 1;
    for (const PrimePower& f : factorize(n)) c *= f.e + 1;
    return c;
}
inline Int128 sum_divisors(long long n) {  // 約数の和（Int128）
    Int128 s = 1;
    for (const PrimePower& f : factorize(n)) {
        Int128 t = 1, pk = 1;
        for (int k = 0; k < f.e; ++k) pk *= f.p, t += pk;
        s *= t;
    }
    return s;
}
inline long long euler_phi(long long n) {  // オイラーの φ(n)（1 以上 n 以下で n と互いに素な数の個数）
    long long r = n;
    for (const PrimePower& f : factorize(n)) r = r / f.p * (f.p - 1);
    return r;
}
inline int mobius(long long n) {  // メビウスの μ(n)（平方因子があれば 0、なければ (-1)^(素因数の数)）
    const std::vector<PrimePower> fs = factorize(n);
    for (const PrimePower& f : fs)
        if (f.e >= 2) return 0;
    return fs.size() % 2 ? -1 : 1;
}
inline int omega(long long n) { return (int)factorize(n).size(); }  // 異なる素因数の個数 ω(n)（12 = 2^2·3 → 2）
inline int big_omega(long long n) {  // 重複を含めた素因数の個数 Ω(n)（12 = 2^2·3 → 3）
    int c = 0;
    for (const PrimePower& f : factorize(n)) c += f.e;
    return c;
}
inline long long largest_prime_factor(long long n) {  // 最大の素因数（n = 1 なら 1）
    const std::vector<PrimePower> fs = factorize(n);
    return fs.empty() ? 1 : fs.back().p;
}
inline SquareFreeParts square_free_decompose(long long n) {  // n = a^2 · b（b は平方因子を持たない）
    long long a = 1, b = 1;
    for (const PrimePower& f : factorize(n)) {
        for (int k = 0; k < f.e / 2; ++k) a *= f.p;
        if (f.e % 2) b *= f.p;
    }
    return {a, b};
}
inline bool is_square_free(long long n) { return square_free_decompose(n).a == 1; }  // 平方因子を持たないか
inline long long radical(long long n) {  // 異なる素因数の積
    long long r = 1;
    for (const PrimePower& f : factorize(n)) r *= f.p;
    return r;
}
inline bool is_prime_power(long long n) {  // n = p^k（p は素数、k ≥ 1）か
    return n >= 2 && factorize(n).size() == 1;
}
inline long long factorial_exponent(long long n, long long p) {  // n! が素数 p で何回割れるか（ルジャンドル。n ≥ 0）
    if (p < 2) math_internal::fail("factorial_exponent", "p ≥ 2 で使う");
    long long c = 0;
    for (long long m = n / p; m > 0; m /= p) c += m;
    return c;
}
inline long long factorial_max_power(long long n, long long k) {  // k^e が n! を割り切る最大の e（k ≥ 2 は合成数でもよい）
    if (k < 2) math_internal::fail("factorial_max_power", "k ≥ 2 で使う");
    long long best = std::numeric_limits<long long>::max();
    for (const PrimePower& f : factorize(k)) best = std::min(best, factorial_exponent(n, f.p) / f.e);
    return best;
}
inline long long min_factorial_multiple(long long k) {  // k が n! を割り切る最小の n ≥ 1（k ≥ 1）
    long long ans = 1;
    for (const PrimePower& f : factorize(k)) {
        long long n = 0, cnt = 0;
        while (cnt < f.e) {
            n += f.p;
            for (long long m = n; m % f.p == 0; m /= f.p) ++cnt;
        }
        ans = std::max(ans, n);
    }
    return ans;
}
inline bool is_perfect_power(long long n) {  // n = a^b（a ≥ 1, b ≥ 2）と書けるか（1 は true。0 以下は false）
    if (n < 1) return false;
    if (n == 1) return true;
    for (int b = 2; b < 63; ++b) {
        const long long a = kth_root_floor(n, b);
        if (a < 2) break;
        if (is_kth_power(n, b)) return true;
    }
    return false;
}
inline long long count_perfect_powers(long long N) {  // 1 以上 N 以下の、a^b（a ≥ 1, b ≥ 2）と書ける数の個数（1 も数える）
    if (N < 1) return 0;
    std::vector<long long> v;  // 平方数でない 3 乗以上の累乗
    for (long long a = 2; a <= N / a / a; ++a) {
        if (is_square(a)) continue;  // a が平方数なら a^b は (√a)^(2b) で、下の数え方に含まれる
        Int128 x = (Int128)a * a * a;
        for (int b = 3; x <= N; ++b, x *= a)
            if (b % 2 == 1) v.push_back((long long)x);  // b が偶数なら平方数
    }
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return isqrt(N) + (long long)v.size();
}
inline long long count_coprime_up_to(long long N, long long M) {  // 1 以上 N 以下で M（≥ 1）と互いに素な数の個数（包除原理）
    if (N < 1) return 0;
    const std::vector<long long> ps = prime_factors(M);
    const int k = (int)ps.size();
    long long r = 0;
    for (int s = 0; s < (1 << k); ++s) {
        Int128 prod = 1;
        int bits = 0;
        for (int i = 0; i < k; ++i)
            if (s >> i & 1) prod *= ps[i], ++bits;
        if (prod > N) continue;
        r += (bits % 2 ? -1 : 1) * (N / (long long)prod);
    }
    return r;
}
template <class V>
long long lcm_mod_of(const V& vs, long long m) {  // 全部の数（≥ 1）の lcm を mod m で（m は何でも。0 があれば 0。vector<int> なども可）
    if (m < 1) math_internal::fail("lcm_mod_of", "m ≥ 1 で使う");
    std::vector<long long> xs;
    for (const auto& xv : vs) {  // 先に全部を確かめる
        const long long x = math_internal::to_ll_checked(xv, "lcm_mod_of");
        if (x < 0) math_internal::fail("lcm_mod_of", "数は 0 以上で使う");
        xs.push_back(x);
    }
    std::vector<PrimePower> best;
    for (const long long x : xs) {
        if (x == 0) return 0;
        for (const PrimePower& f : factorize(x)) best.push_back(f);
    }
    std::sort(best.begin(), best.end(), [](const PrimePower& a, const PrimePower& b) { return a.p != b.p ? a.p < b.p : a.e > b.e; });
    long long r = 1 % m;
    for (size_t i = 0; i < best.size(); ++i)
        if (i == 0 || best[i].p != best[i - 1].p) r = mod_mul(r, mod_pow(best[i].p, best[i].e, m), m);
    return r;
}
inline long long lcm_mod_of(std::initializer_list<long long> vs, long long m) { return lcm_mod_of<std::initializer_list<long long>>(vs, m); }
#endif  // GG_FACTORIZE
