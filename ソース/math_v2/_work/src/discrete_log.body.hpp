#ifndef GG_DISCRETE_LOG
#define GG_DISCRETE_LOG
// ---- 位数・原始根・離散対数・平方根（数学第二世代）  m は 10^18 程度まで（掛け算は 128 ビット） ----
// カーマイケル関数 λ(n)：すべての gcd(a, n) = 1 で a^λ ≡ 1 となる最小の数  O(n の素因数分解)
inline long long carmichael_lambda(long long n) {
    if (n < 1) math_internal::fail("carmichael_lambda", "n ≥ 1 で使う");
    long long l = 1;
    if (n == 1) return 1;
    for (const PrimePower& f : factorize(n)) {
        long long t;
        if (f.p == 2) t = f.e <= 2 ? (1LL << (f.e - 1)) : (1LL << (f.e - 2));
        else t = pow_int(f.p, f.e - 1) * (f.p - 1);
        l = l / (long long)math_internal::ugcd((unsigned long long)l, (unsigned long long)t) * t;
    }
    return l;
}
// 位数：a^k ≡ 1 (mod m) となる最小の k ≥ 1。同じ m で何度も聞くときは MultOrder ord(m); ord(a)（λ(m) の分解を 1 回だけ）
struct MultOrder {
    explicit MultOrder(long long m) : m_(m) {
        if (m < 1) math_internal::fail("MultOrder", "m ≥ 1 で使う");
        lam_ = carmichael_lambda(m);
        if (lam_ > 1)
            for (const PrimePower& f : factorize(lam_)) qs_.push_back(f.p);
    }
    Maybe<long long> operator()(long long a) const {  // gcd(a, m) ≠ 1 なら ok = false
        a = mod_floor(a, m_);
        if (m_ == 1) return {true, 1};
        if (math_internal::ugcd((unsigned long long)a, (unsigned long long)m_) != 1) return {false, 0};
        long long k = lam_;
        for (long long q : qs_)
            while (k % q == 0 && mod_pow(a, k / q, m_) == 1) k /= q;
        return {true, k};
    }
    long long lambda() const { return lam_; }

  private:
    long long m_, lam_;
    std::vector<long long> qs_;
};
inline Maybe<long long> mult_order(long long a, long long m) { return MultOrder(m)(a); }
// 素数 p の原始根（いちばん小さいもの。p ≤ 10^18）
inline long long primitive_root(long long p) {
    if (p < 2 || !is_prime(p)) math_internal::fail("primitive_root", "p は素数");
    if (p == 2) return 1;
    const std::vector<long long> qs = prime_factors(p - 1);
    for (long long g = 2;; ++g) {
        bool ok = true;
        for (long long q : qs)
            if (mod_pow(g, (p - 1) / q, p) == 1) { ok = false; break; }
        if (ok) return g;
    }
}
// x^2 ≡ a (mod p) の x（p は素数。もう 1 つの解は p − x。なければ ok = false。トネリ・シャンクス）
inline Maybe<long long> mod_sqrt(long long a, long long p) {
    if (p < 2 || !is_prime(p)) math_internal::fail("mod_sqrt", "p は素数");
    a = mod_floor(a, p);
    if (a == 0 || p == 2) return {true, a};
    if (mod_pow(a, (p - 1) / 2, p) != 1) return {false, 0};
    long long q = p - 1;
    int s = 0;
    while (q % 2 == 0) q /= 2, ++s;
    long long z = 2;
    while (mod_pow(z, (p - 1) / 2, p) == 1) ++z;
    long long c = mod_pow(z, q, p), x = mod_pow(a, (q + 1) / 2, p), t = mod_pow(a, q, p);
    int m = s;
    while (t != 1) {
        int i = 0;
        for (long long tt = t; tt != 1; tt = mod_mul(tt, tt, p)) ++i;
        long long b = c;
        for (int j = 0; j < m - i - 1; ++j) b = mod_mul(b, b, p);
        x = mod_mul(x, b, p), c = mod_mul(b, b, p), t = mod_mul(t, c, p), m = i;
    }
    return {true, std::min(x, p - x)};
}
// a^x ≡ b (mod m) となる最小の x ≥ 0（m は合成数でもよい。0^0 = 1。なければ ok = false。BSGS：m ≤ 10^13 程度）  O(√m)
inline Maybe<long long> discrete_log(long long a, long long b, long long m) {
    if (m < 1) math_internal::fail("discrete_log", "m ≥ 1 で使う");
    a = mod_floor(a, m), b = mod_floor(b, m);
    if (m == 1) return {true, 0};
    {  // 小さい x はそのまま調べる（a と m が互いに素でないときの最初の部分も）
        long long cur = 1 % m;
        for (long long x = 0; x < 64; ++x) {
            if (cur == b) return {true, x};
            cur = mod_mul(cur, a, m);
        }
    }
    // a^x = a^k · a^{x−k}：gcd(a, m) を m から除いていく
    long long k = 0, c = 1 % m;
    for (long long g; (g = (long long)math_internal::ugcd((unsigned long long)a, (unsigned long long)m)) != 1;) {
        if (b % g != 0) return {false, 0};
        b /= g, m /= g, ++k;
        c = mod_mul(c, a / g, m);
    }
    // c · a^y ≡ b (mod m)、gcd(a, m) = 1：a^y ≡ b c^{−1}
    const long long t = mod_mul(b, mod_inv(c, m), m);
    const long long n = isqrt(m) + 1;
    if (n > 4000000) math_internal::fail("discrete_log", "m が大きすぎる（√m の表を作るので m ≤ 10^13 程度）");
    std::unordered_map<long long, long long> baby;  // a^j·t → j（大きい方の j を残す）
    baby.reserve((size_t)n * 2);
    long long cur = t;
    for (long long j = 0; j < n; ++j) baby[cur] = j, cur = mod_mul(cur, a, m);
    const long long an = mod_pow(a, n, m);
    long long giant = an;
    for (long long i = 1; i <= n; ++i) {  // a^{i n} = a^j t ⇒ y = i n − j
        const auto it = baby.find(giant);
        if (it != baby.end()) return {true, i * n - it->second + k};
        giant = mod_mul(giant, an, m);
    }
    return {false, 0};
}
// a_0^(a_1^(a_2^(…))) mod m（a_i ≥ 0。0^0 = 1。指数は φ(m) で割った余り + φ(m) にして小さくする）
inline long long pow_tower_mod(const std::vector<long long>& as, long long m) {
    if (m < 1) math_internal::fail("pow_tower_mod", "m ≥ 1 で使う");
    for (long long a : as)
        if (a < 0) math_internal::fail("pow_tower_mod", "a_i ≥ 0 で使う");
    if (as.empty()) return 1 % m;
    const size_t n = as.size();
    // cap[i] = min(a_i^(a_{i+1}^…), 2^62)（正確な小さい値）
    const long long CAP = 1LL << 62;
    std::vector<long long> cap(n);
    cap[n - 1] = std::min(as[n - 1], CAP);
    for (size_t i = n - 1; i-- > 0;) cap[i] = pow_cap(as[i], cap[i + 1], CAP);
    // 上から mod を φ にしながら
    std::function<long long(size_t, long long)> rec = [&](size_t i, long long mod) -> long long {
        if (mod == 1) return 0;
        if (i == n - 1) return as[i] % mod;
        if (cap[i + 1] < 64) return mod_pow(as[i], cap[i + 1], mod);  // 指数が小さい（正確に分かる）
        const long long ph = euler_phi(mod);                            // 指数 ≥ 64 ≥ mod の素因数の重複度
        return mod_pow(as[i], rec(i + 1, ph) + ph, mod);
    };
    return rec(0, m);
}
#endif  // GG_DISCRETE_LOG
