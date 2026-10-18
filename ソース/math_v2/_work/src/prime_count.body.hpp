#ifndef GG_PRIME_COUNT
#define GG_PRIME_COUNT
// ---- 素数の個数・和（数学第二世代）  N が大きい（10^11〜10^12）π(N)・素数の和、⌊N/i⌋ ごとの表、P-smooth な数の個数、乗法的関数の和（Min_25） ----
// ⌊N/i⌋ の形の値 v ごとの π(v)（Lucy の方法）。P.count(v) は v ≤ √N か v = ⌊N/i⌋ のとき  O(N^{3/4} / log N)
struct PrimeCountTable {
    explicit PrimeCountTable(long long N) : n_(N) {
        if (N < 0) math_internal::fail("PrimeCountTable", "N ≥ 0 で使う");
        r_ = isqrt(N);
        lo_.assign((size_t)r_ + 1, 0), hi_.assign((size_t)r_ + 1, 0);
        for (long long v = 1; v <= r_; ++v) lo_[(size_t)v] = v - 1;    // lo[v] = π(v) の候補（はじめは 2〜v の個数）
        for (long long i = 1; i <= r_; ++i) hi_[(size_t)i] = N / i - 1;  // hi[i] = π(N/i) の候補
        for (long long p = 2; p <= r_; ++p) {
            if (lo_[(size_t)p] == lo_[(size_t)p - 1]) continue;  // p は合成数
            const long long pc = lo_[(size_t)p - 1], p2 = p * p, lim = std::min(r_, N / p2);
            for (long long i = 1; i <= lim; ++i) {  // N/i ≥ p^2 のもの：p を最小の素因数に持つ数を除く
                const long long d = i * p;
                hi_[(size_t)i] -= (d <= r_ ? hi_[(size_t)d] : lo_[(size_t)(N / d)]) - pc;
            }
            for (long long v = r_; v >= p2; --v) lo_[(size_t)v] -= lo_[(size_t)(v / p)] - pc;
        }
    }
    long long count(long long v) const {  // π(v)（v ≤ √N か v = ⌊N/i⌋）
        if (v < 0) return 0;
        if (v <= r_) return lo_[(size_t)v];
        if (v > n_ || n_ / (n_ / v) != v) math_internal::fail("PrimeCountTable::count", "v は √N 以下か ⌊N/i⌋ の形");
        return hi_[(size_t)(n_ / v)];
    }
    long long operator()(long long v) const { return count(v); }
    long long n() const { return n_; }

  private:
    long long n_, r_;
    std::vector<long long> lo_, hi_;
};
inline long long prime_count(long long N) { return N < 2 ? 0 : PrimeCountTable(N).count(N); }  // N 以下の素数の個数

// N 以下の素数の和（Int128。N ≤ 10^12 程度）  O(N^{3/4} / log N)
inline Int128 prime_sum(long long N) {
    if (N < 2) return 0;
    const long long r = isqrt(N);
    std::vector<Int128> lo((size_t)r + 1, 0), hi((size_t)r + 1, 0);
    auto tri = [](long long v) { return (Int128)v * (v + 1) / 2 - 1; };  // 2 + 3 + … + v
    for (long long v = 1; v <= r; ++v) lo[(size_t)v] = tri(v);
    for (long long i = 1; i <= r; ++i) hi[(size_t)i] = tri(N / i);
    for (long long p = 2; p <= r; ++p) {
        if (lo[(size_t)p] == lo[(size_t)p - 1]) continue;
        const Int128 pc = lo[(size_t)p - 1];
        const long long p2 = p * p, lim = std::min(r, N / p2);
        for (long long i = 1; i <= lim; ++i) {
            const long long d = i * p;
            hi[(size_t)i] -= (Int128)p * ((d <= r ? hi[(size_t)d] : lo[(size_t)(N / d)]) - pc);
        }
        for (long long v = r; v >= p2; --v) lo[(size_t)v] -= (Int128)p * (lo[(size_t)(v / p)] - pc);
    }
    return hi[1];
}

// N 以下の正の整数で、素因数がすべて P 以下のものの個数（1 も数える。P ≤ 100 程度・N ≤ 10^18。素数を 2 組に分けて数え合わせる）
inline long long count_smooth(long long N, long long P) {
    if (N < 1) return 0;
    if (P > 1000) math_internal::fail("count_smooth", "P ≤ 1000 で使う（N が大きいときは P ≤ 100 程度）");
    std::vector<long long> ps;
    for (long long p = 2; p <= P; ++p) {
        bool pr = true;
        for (long long q = 2; q * q <= p; ++q) pr &= p % q != 0;
        if (pr) ps.push_back(p);
    }
    std::vector<long long> g1, g2;  // 素数を交互に 2 組に
    for (size_t i = 0; i < ps.size(); ++i) (i % 2 ? g2 : g1).push_back(ps[i]);
    size_t total = 0;
    auto enumerate = [&](const std::vector<long long>& g) {  // g の素数だけでできた N 以下の数（1 も）
        std::vector<long long> v = {1};
        for (long long p : g) {
            const size_t m = v.size();
            for (size_t i = 0; i < m; ++i)
                for (long long x = v[i]; x <= N / p;) {
                    x *= p, v.push_back(x);
                    if (++total > 40000000) math_internal::fail("count_smooth", "途中で作る数が多すぎる（N か P が大きすぎる。(10^18, 97) くらいまで）");
                }
        }
        std::sort(v.begin(), v.end());
        return v;
    };
    const std::vector<long long> A = enumerate(g1), B = enumerate(g2);
    long long cnt = 0;
    long long j = (long long)B.size() - 1;
    for (long long a : A) {  // a·b ≤ N となる b の数（a が増えると b の上限は減る）
        while (j >= 0 && B[(size_t)j] > N / a) --j;
        if (j < 0) break;
        cnt += j + 1;
    }
    return cnt;
}

// ---- 乗法的関数の和（Min_25 篩）----
namespace math_internal {
template <class T>
constexpr bool pc_is_int = std::is_integral<T>::value || std::is_same<T, Int128>::value;
template <class T>
inline T pc_add(const T& a, const T& b) {
    if constexpr (pc_is_int<T>) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail("multiplicative_sum", "値が整数の型に収まらない（mint か Int128 で）");
        return r;
    } else {
        return a + b;
    }
}
template <class T>
inline T pc_sub(const T& a, const T& b) {
    if constexpr (pc_is_int<T>) {
        T r;
        if (__builtin_sub_overflow(a, b, &r)) fail("multiplicative_sum", "値が整数の型に収まらない（mint か Int128 で）");
        return r;
    } else {
        return a - b;
    }
}
template <class T>
inline T pc_mul(const T& a, const T& b) {
    if constexpr (pc_is_int<T>) {
        T r;
        if (__builtin_mul_overflow(a, b, &r)) fail("multiplicative_sum", "値が整数の型に収まらない（mint か Int128 で）");
        return r;
    } else {
        return a * b;
    }
}
template <class T>
T pc_power_sum(long long v, int k) {  // 1^k + 2^k + … + v^k（k ≤ 3。割り算は先に整数で済ませる）
    if (k == 0) return T(v);
    long long a = v, b = v + 1;
    (a % 2 == 0 ? a : b) /= 2;
    if (k == 1) return pc_mul(T(a), T(b));
    if (k == 3) {
        const T s = pc_mul(T(a), T(b));
        return pc_mul(s, s);
    }
    long long x = v, y = v + 1, z = 2 * v + 1;  // k = 2：v(v+1)(2v+1)/6
    (x % 2 == 0 ? x : y) /= 2;
    (x % 3 == 0 ? x : (y % 3 == 0 ? y : z)) /= 3;
    return pc_mul(pc_mul(T(x), T(y)), T(z));
}
}  // namespace math_internal
// Σ_{i=1}^{N} f(i)：f は乗法的（f(1) = 1）で、素数 p では f(p) = coef[0] + coef[1] p + coef[2] p^2 + coef[3] p^3（3 次まで）、
// 素数の累乗では f(p^e) = fpe(p, e)（e ≥ 1。fpe(p, 1) は f(p) と同じ値にする）。T は mint・long long・Int128（整数は溢れたら止める）
//   例：約数の個数 d … coef {2}、fpe = e + 1 ／ φ … coef {−1, 1}、fpe = p^e − p^(e−1) ／ μ … coef {−1}、fpe = (e == 1 ? −1 : 0)
//   N ≤ 10^12（mint で 10^10 なら約 0.2 秒、10^11 で約 1 秒、10^12 で約 4 秒）  O(N^{3/4} / log N) 程度
template <class T, class Fpe>
T multiplicative_sum(long long N, const std::vector<T>& coef, Fpe fpe) {
    using namespace math_internal;
    if (N < 1) return T(0);
    if (coef.empty() || coef.size() > 4) fail("multiplicative_sum", "coef は 1〜4 個（f(p) は p の 3 次まで）");
    if (N > 1000000000000LL) fail("multiplicative_sum", "N ≤ 10^12 で使う");
    const long long r = isqrt(N);
    std::vector<long long> ps;  // √N 以下の素数
    {
        std::vector<char> comp((size_t)r + 1, 0);
        for (long long p = 2; p <= r; ++p)
            if (!comp[(size_t)p]) {
                ps.push_back(p);
                for (long long q = p * p; q <= r; q += p) comp[(size_t)q] = 1;
            }
    }
    // ⌊N/i⌋ の形の値 v での Σ_{素数 q ≤ v} f(q)：v ≤ r は Glo[v]、v > r は Ghi[N/v]（Lucy の方法を p^k ごとに）
    std::vector<T> Glo((size_t)r + 1, T(0)), Ghi((size_t)r + 1, T(0));
    for (int k = 0; k < (int)coef.size(); ++k) {
        if (coef[(size_t)k] == T(0)) continue;
        std::vector<T> lo((size_t)r + 1), hi((size_t)r + 1);
        lo[0] = T(0);
        for (long long v = 1; v <= r; ++v) lo[(size_t)v] = pc_sub(pc_power_sum<T>(v, k), T(1));
        for (long long i = 1; i <= r; ++i) hi[(size_t)i] = pc_sub(pc_power_sum<T>(N / i, k), T(1));
        for (const long long p : ps) {
            T pw = T(1);
            for (int t = 0; t < k; ++t) pw = pc_mul(pw, T(p));
            const T pc = lo[(size_t)p - 1];  // p より小さい素数の q^k の和
            const long long p2 = p * p, lim = std::min(r, N / p2);
            if (k == 0) {  // p^0 = 1 なので掛け算なし
                for (long long i = 1; i <= lim; ++i) {
                    const long long d = i * p;
                    hi[(size_t)i] = pc_sub(hi[(size_t)i], pc_sub(d <= r ? hi[(size_t)d] : lo[(size_t)(N / d)], pc));
                }
                for (long long v = r; v >= p2; --v) lo[(size_t)v] = pc_sub(lo[(size_t)v], pc_sub(lo[(size_t)(v / p)], pc));
            } else {
                for (long long i = 1; i <= lim; ++i) {
                    const long long d = i * p;
                    hi[(size_t)i] = pc_sub(hi[(size_t)i], pc_mul(pw, pc_sub(d <= r ? hi[(size_t)d] : lo[(size_t)(N / d)], pc)));
                }
                for (long long v = r; v >= p2; --v) lo[(size_t)v] = pc_sub(lo[(size_t)v], pc_mul(pw, pc_sub(lo[(size_t)(v / p)], pc)));
            }
        }
        for (long long v = 1; v <= r; ++v) Glo[(size_t)v] = pc_add(Glo[(size_t)v], pc_mul(coef[(size_t)k], lo[(size_t)v]));
        for (long long i = 1; i <= r; ++i) Ghi[(size_t)i] = pc_add(Ghi[(size_t)i], pc_mul(coef[(size_t)k], hi[(size_t)i]));
    }
    const int P = (int)ps.size();
    // S(v, j) = Σ_{2 ≤ x ≤ v、x の最小の素因数 ≥ ps[j]} f(x)
    auto S = [&](auto&& self, long long v, int j) -> T {
        const long long last = j == 0 ? 1 : ps[(size_t)j - 1];
        if (v <= last) return T(0);
        T res = pc_sub(v <= r ? Glo[(size_t)v] : Ghi[(size_t)(N / v)], Glo[(size_t)last]);  // ps[j] 以上の素数
        for (int k = j; k < P && ps[(size_t)k] * ps[(size_t)k] <= v; ++k) {
            const long long p = ps[(size_t)k];
            long long pe = p;
            T fe = T(fpe(p, 1));
            for (int e = 1; pe <= v / p; ++e, pe *= p) {
                const T fe1 = T(fpe(p, e + 1));
                res = pc_add(res, pc_add(pc_mul(fe, self(self, v / pe, k + 1)), fe1));
                fe = fe1;
            }
        }
        return res;
    };
    return pc_add(S(S, N, 0), T(1));
}
#endif  // GG_PRIME_COUNT
