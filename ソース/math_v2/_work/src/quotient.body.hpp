#ifndef GG_QUOTIENT
#define GG_QUOTIENT
// ---- 商の列挙・床関数の和（数学第二世代）  ⌊N/i⌋ が同じ i をまとめる、Σ⌊N/i⌋、Σ⌊(a i + b)/m⌋・Σ (a i + b) mod m・Σ⌊⌋²、万能ユークリッド ----
struct QuotientBlock {  // l ≤ i < r で ⌊N/i⌋（ceil 版は ⌈N/i⌉）が q
    long long q, l, r;
};
struct QuotientBlock2 {  // l ≤ i < r で ⌊N/i⌋ = q1、⌊M/i⌋ = q2
    long long q1, q2, l, r;
};

// i = 1, 2, …, N を ⌊N/i⌋ が同じ区間に分けて、i の小さい順に f(q, l, r)（区間の数は 2√N 以下）  O(√N)
//   f が bool を返すなら、false を返したところで止める
template <class F>
void for_each_quotient(long long N, F f) {  // （0 ≤ N < 2^63 − 1。区間の右端 r = N + 1 が long long に入るように）
    if (N < 0 || N == 9223372036854775807LL) math_internal::fail("for_each_quotient", "0 ≤ N < LLONG_MAX で使う");
    for (long long l = 1; l <= N;) {
        const long long q = N / l, r = N / q + 1;
        if constexpr (std::is_same<decltype(f(q, l, r)), bool>::value) {
            if (!f(q, l, r)) return;
        } else {
            f(q, l, r);
        }
        l = r;
    }
}
inline std::vector<QuotientBlock> quotient_blocks(long long N) {  // 上を vector で（q は大きい順）
    std::vector<QuotientBlock> v;
    for_each_quotient(N, [&](long long q, long long l, long long r) { v.push_back(QuotientBlock{q, l, r}); });
    return v;
}
inline std::vector<QuotientBlock> quotient_blocks_ceil(long long N) {  // i = 1〜N を ⌈N/i⌉ が同じ区間に分ける（N ≥ 1）
    if (N < 0 || N == 9223372036854775807LL) math_internal::fail("quotient_blocks_ceil", "0 ≤ N < LLONG_MAX で使う");
    std::vector<QuotientBlock> v;
    for (long long l = 1; l <= N;) {
        const long long q = (N - 1) / l + 1;              // ⌈N/l⌉（N ≥ 1。N + l − 1 は溢れることがある）
        const long long r = q == 1 ? N + 1 : (N - 1) / (q - 1) + 1;  // ⌈N/i⌉ = q となる最大の i は ⌊(N−1)/(q−1)⌋
        v.push_back(QuotientBlock{q, l, r});
        l = r;
    }
    return v;
}
inline std::vector<QuotientBlock2> quotient_blocks2(long long N, long long M) {  // i = 1〜min(N, M) を ⌊N/i⌋・⌊M/i⌋ が両方同じ区間に
    if (N < 0 || M < 0 || std::min(N, M) == 9223372036854775807LL) math_internal::fail("quotient_blocks2", "N, M ≥ 0、min(N, M) < LLONG_MAX で使う");
    std::vector<QuotientBlock2> v;
    const long long lim = std::min(N, M);
    for (long long l = 1; l <= lim;) {
        const long long q1 = N / l, q2 = M / l, r = std::min(N / q1, M / q2) + 1;
        v.push_back(QuotientBlock2{q1, q2, l, r});
        l = r;
    }
    return v;
}

// ---- ⌊N/i⌋ の和（N ≥ 0。答えは Int128） ----
inline Int128 sum_floor_div(long long N) {  // Σ_{i=1}^{N} ⌊N/i⌋（= 1〜N の約数の個数の和）  O(√N)
    if (N < 0) math_internal::fail("sum_floor_div", "N ≥ 0 で使う");
    const long long s = isqrt(N);
    Int128 r = 0;
    for (long long i = 1; i <= s; ++i) r += N / i;
    return 2 * r - (Int128)s * s;  // i ≤ √N と ⌊N/i⌋ ≤ √N の対称性
}
inline Int128 sum_i_floor_div(long long N) {  // Σ_{i=1}^{N} i·⌊N/i⌋（= 1〜N の約数の和の和）  O(√N)
    if (N < 0) math_internal::fail("sum_i_floor_div", "N ≥ 0 で使う");
    Int128 r = 0;
    for_each_quotient(N, [&](long long q, long long l, long long R) {
        const long long w = R - l;
        const Int128 t = (Int128)l + R - 1;  // l + … + (R − 1) = w·t / 2（どちらかは偶数。t は 2^63 を超えることがある）
        r += (Int128)q * (w % 2 == 0 ? (Int128)(w / 2) * t : (Int128)w * (t / 2));
    });
    return r;
}
inline Int128 sum_mod_i(long long N) {  // Σ_{i=1}^{N} (N mod i)  O(√N)
    return (Int128)N * N - sum_i_floor_div(N);
}

// ---- Σ_{i=0}^{n−1} ⌊(a i + b) / m⌋（a, b は負でもよい。n ≥ 0、m ≥ 1。答えは Int128。途中で溢れたら止める）  O(log m) ----
namespace math_internal {
inline Int128 qmul(Int128 x, Int128 y) {
    Int128 r;
    if (__builtin_mul_overflow(x, y, &r)) fail("linear_floor_sum", "Int128 に収まらない");
    return r;
}
inline Int128 qadd(Int128 x, Int128 y) {
    Int128 r;
    if (__builtin_add_overflow(x, y, &r)) fail("linear_floor_sum", "Int128 に収まらない");
    return r;
}
inline Int128 tri_count(Int128 n) { return (n % 2 == 0) ? qmul(n / 2, n - 1) : qmul(n, (n - 1) / 2); }  // n(n−1)/2
// a, b ≥ 0（ACL と同じユークリッドの互除法の形）。U は n, m < 2^32 なら unsigned long long（速い）、それ以外は Int128
template <class U>
inline Int128 floor_sum_nonneg(U n, U m, U a, U b) {
    Int128 ans = 0;
    while (true) {
        if constexpr (std::is_same_v<U, unsigned long long>) {  // n, m < 2^32 なら答えは 2^64 未満なので確認しない
            if (a >= m) ans += (Int128)(n * (n - 1) / 2) * (Int128)(a / m), a %= m;
            if (b >= m) ans += (Int128)n * (Int128)(b / m), b %= m;
        } else {
            if (a >= m) ans = qadd(ans, qmul(tri_count((Int128)n), (Int128)(a / m))), a %= m;
            if (b >= m) ans = qadd(ans, qmul((Int128)n, (Int128)(b / m))), b %= m;
        }
        const U y_max = a * n + b;  // a < m、b < m なので m(n+1) 未満（n は増えず、m は減る）
        if (y_max < m) break;
        n = y_max / m, b = y_max % m;
        std::swap(m, a);
    }
    return ans;
}
}  // namespace math_internal
inline Int128 linear_floor_sum(long long n, long long m, long long a, long long b) {
    if (n < 0 || m < 1) math_internal::fail("linear_floor_sum", "n ≥ 0、m ≥ 1 で使う");
    if (n == 0) return 0;
    if (n < (1LL << 32) && m < (1LL << 32)) {  // よくある大きさ：64 ビットで計算（答えは 2^127 未満なので溢れない）
        const long long ka = floor_div(a, m), kb = floor_div(b, m);
        const unsigned long long un = (unsigned long long)n;
        return (Int128)(un * (un - 1) / 2) * ka + (Int128)n * kb +
               math_internal::floor_sum_nonneg<unsigned long long>(un, (unsigned long long)m, (unsigned long long)mod_floor(a, m), (unsigned long long)mod_floor(b, m));
    }
    using math_internal::qadd;
    using math_internal::qmul;
    Int128 A = a, B = b, ans = 0;
    const Int128 N = n, M = m;
    if (A < 0 || A >= M) {  // a = a' + k m（0 ≤ a' < m）にして k·n(n−1)/2 を外に出す
        const Int128 k = math_internal::floor_div128(A, M);
        ans = qadd(ans, qmul(math_internal::tri_count(N), k)), A -= k * M;
    }
    if (B < 0 || B >= M) {
        const Int128 k = math_internal::floor_div128(B, M);
        ans = qadd(ans, qmul(N, k)), B -= k * M;
    }
    return qadd(ans, math_internal::floor_sum_nonneg<Int128>(N, M, A, B));
}
// a x + b y ≤ c となる x, y ≥ 0 の整数の組の数（a, b ≥ 1。c < 0 なら 0）  O(log)
inline Int128 count_lattice_under_line(long long a, long long b, long long c) {
    if (a < 1 || b < 1) math_internal::fail("count_lattice_under_line", "a, b ≥ 1 で使う");
    if (c < 0) return 0;
    // x = 0 は ⌊c/b⌋ + 1 個。x = 1〜q（q = ⌊c/a⌋）は Σ_{i=0}^{q−1} (⌊(c − a − a i)/b⌋ + 1)（q + 1 は溢れることがあるので分ける）
    const long long q = c / a;
    return linear_floor_sum(q, b, -a, c - a) + q + c / b + 1;
}
// Σ_{i=0}^{n−1} ((a i + b) mod m)（mod は 0〜m−1 の方。a, b は負でも。答えは n(m−1) 以下なので Int128）  O(log m)
inline Int128 linear_mod_sum(long long n, long long m, long long a, long long b) {
    if (n < 0 || m < 1) math_internal::fail("linear_mod_sum", "n ≥ 0、m ≥ 1 で使う");
    if (n == 0) return 0;
    const long long a2 = mod_floor(a, m), b2 = mod_floor(b, m);
    // Σ (a2 i + b2) − m Σ ⌊(a2 i + b2)/m⌋。途中は 2^128 を法として計算してよい（答えは 0 以上 2^126 未満）
    const UInt128 N = (UInt128)n, T = (N % 2 == 0) ? (N / 2) * (N - 1) : N * ((N - 1) / 2);
    const UInt128 F = (UInt128)linear_floor_sum(n, m, a2, b2);
    return (Int128)((UInt128)a2 * T + (UInt128)b2 * N - (UInt128)m * F);
}

// ---- 床関数の和の拡張：Σ⌊⌋・Σ i⌊⌋・Σ⌊⌋² を mint で（mod は奇数の素数）----
template <class T>
struct FloorSumExt {
    T sum, sum_i, sum_sq;  // Σ_{i=0}^{n−1} ⌊(a i + b)/m⌋、Σ i⌊…⌋、Σ ⌊…⌋^2
};
namespace math_internal {
template <class T>
T fs_s1(long long n) {  // 0 + 1 + … + n（n ≥ −1）
    long long x = n, y = n + 1;
    (x % 2 == 0 ? x : y) /= 2;
    return T(x) * T(y);
}
template <class T>
T fs_from128(Int128 v) {  // 0 ≤ v < 2^125 を T に（T は long long から作れればよい）
    return T((long long)(v >> 62)) * T(1LL << 62) + T((long long)(v & ((1LL << 62) - 1)));
}
template <class T>
T fs_s2(long long n) {  // 0^2 + 1^2 + … + n^2
    Int128 x = n, y = (Int128)n + 1, z = 2 * (Int128)n + 1;
    (x % 2 == 0 ? x : y) /= 2;
    (x % 3 == 0 ? x : (y % 3 == 0 ? y : z)) /= 3;
    return fs_from128<T>(x) * fs_from128<T>(y) * fs_from128<T>(z);
}
template <class T>
FloorSumExt<T> fs_ext(long long n, long long a, long long b, long long c) {  // x = 0〜n、a, b ≥ 0、c ≥ 1
    if (n < 0) return {T(0), T(0), T(0)};
    const T N(n), S1 = fs_s1<T>(n), S2 = fs_s2<T>(n);
    if (a == 0) {
        const T q(b / c);
        return {q * (N + T(1)), q * S1, q * q * (N + T(1))};
    }
    if (a >= c || b >= c) {
        const FloorSumExt<T> t = fs_ext<T>(n, a % c, b % c, c);
        const T qa(a / c), qb(b / c);
        FloorSumExt<T> r;
        r.sum = t.sum + qa * S1 + qb * (N + T(1));
        r.sum_i = t.sum_i + qa * S2 + qb * S1;
        r.sum_sq = qa * qa * S2 + qb * qb * (N + T(1)) + T(2) * qa * qb * S1 + t.sum_sq + T(2) * qb * t.sum + T(2) * qa * t.sum_i;
        return r;
    }
    const long long m = (long long)(((Int128)a * n + b) / c);  // a, b < c なので n + 1 未満
    const FloorSumExt<T> t = fs_ext<T>(m - 1, c, c - b - 1, a);
    const T M(m);
    FloorSumExt<T> r;
    r.sum = N * M - t.sum;
    r.sum_i = (M * N * (N + T(1)) - t.sum_sq - t.sum) / T(2);
    r.sum_sq = N * M * (M + T(1)) - T(2) * t.sum_i - T(2) * t.sum - r.sum;
    return r;
}
}  // namespace math_internal
// n ≥ 0、m ≥ 1、a, b は負でもよい。T は mint（mod は奇数の素数）  O(log m)
template <class T>
FloorSumExt<T> linear_floor_sum_ext(long long n, long long m, long long a, long long b) {
    if (n < 0 || m < 1) math_internal::fail("linear_floor_sum_ext", "n ≥ 0、m ≥ 1 で使う");
    if (n == 0) return {T(0), T(0), T(0)};
    const long long qa = floor_div(a, m), qb = floor_div(b, m);  // ⌊(a i + b)/m⌋ = qa i + qb + ⌊(a' i + b')/m⌋
    const FloorSumExt<T> t = math_internal::fs_ext<T>(n - 1, mod_floor(a, m), mod_floor(b, m), m);
    const T A(qa), B(qb), N(n), S1 = math_internal::fs_s1<T>(n - 1), S2 = math_internal::fs_s2<T>(n - 1);
    FloorSumExt<T> r;
    r.sum = t.sum + A * S1 + B * N;
    r.sum_i = t.sum_i + A * S2 + B * S1;
    r.sum_sq = A * A * S2 + B * B * N + T(2) * A * B * S1 + t.sum_sq + T(2) * B * t.sum + T(2) * A * t.sum_i;
    return r;
}

// ---- 万能ユークリッド（universal Euclid）----
// 直線 y = ⌊(a x + b)/c⌋ に沿って、x を 1 増やすたびに R、y が 1 増えるたびに U を掛けた積
//   U^{y(0)} · Π_{x=1}^{n} ( U^{y(x) − y(x−1)} · R )（a, b ≥ 0、c ≥ 1、n ≥ 0）を求める。
//   M は結合的な掛け算 operator* と単位元 M()（既定のコンストラクタ）を持つ型（行列・自作の struct など。可換でなくてよい）
//   Σ_x A^x B^{y(x)} のような「床関数を含む和」を、M に必要な量を持たせて求める  O(log max(a, c) × M の掛け算 × log)
namespace math_internal {
template <class M>
M ue_pow(M x, long long e) {  // 最後の 2 乗はしない（使わない値で溢れないように）
    M r;
    while (e > 0) {
        if (e & 1) r = r * x;
        e >>= 1;
        if (e > 0) x = x * x;
    }
    return r;
}
template <class M>
M ue_solve(long long p, long long q, long long r, long long l, const M& U, const M& R) {  // 0 ≤ r < q。y = ⌊(p x + r)/q⌋、x = 1〜l
    if (l == 0) return M();
    if (p >= q) return ue_solve<M>(p % q, q, r, l, U, ue_pow(U, p / q) * R);
    const long long m = (long long)(((Int128)l * p + r) / q);
    if (m == 0) return ue_pow(R, l);
    const long long cnt = l - (long long)(((Int128)q * m - r - 1) / p);
    return ue_pow(R, (q - r - 1) / p) * U * ue_solve<M>(q, p, (q - r - 1) % p, m - 1, R, U) * ue_pow(R, cnt);
}
}  // namespace math_internal
template <class M>
M universal_euclid(long long n, long long a, long long b, long long c, const M& U, const M& R) {
    if (n < 0 || a < 0 || b < 0 || c < 1) math_internal::fail("universal_euclid", "n, a, b ≥ 0、c ≥ 1 で使う");
    return math_internal::ue_pow(U, b / c) * math_internal::ue_solve<M>(a, c, b % c, n, U, R);
}
#endif  // GG_QUOTIENT
