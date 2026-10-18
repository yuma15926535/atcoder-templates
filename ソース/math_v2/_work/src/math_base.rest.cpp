// =============================================================================
//  ランダムテスト
// =============================================================================
namespace mathb_test {

{COMMON_TEST}

static const long long LMAX = std::numeric_limits<long long>::max(), LMIN = std::numeric_limits<long long>::min();
static int sgn128(Int128 v) { return (v > 0) - (v < 0); }
// 愚直な ⌊a/b⌋（引き算で余りを 0 以上にそろえる）
static Int128 bf_floor(Int128 a, Int128 b) {
    if (b < 0) a = -a, b = -b;
    Int128 q = a / b, r = a - q * b;
    while (r < 0) r += b, --q;
    while (r >= b) r -= b, ++q;
    return q;
}
static Int128 bf_ceil(Int128 a, Int128 b) { return -bf_floor(-a, b); }
static unsigned long long bf_gcd(unsigned long long a, unsigned long long b) { return b ? bf_gcd(b, a % b) : a; }
static Int128 bf_gcd128(Int128 a, Int128 b) { return b ? bf_gcd128(b, a % b) : a; }
// y^k を ±2^100 で止める（|y| < 2^64）
static Int128 bf_pow_sat(Int128 y, long long k) {
    const Int128 SAT = (Int128)1 << 100;
    Int128 r = 1;
    const Int128 ay = y < 0 ? -y : y;
    for (long long i = 0; i < k; ++i) {
        if (r == 0) return 0;
        const Int128 ar = r < 0 ? -r : r;
        if (ay > 1 && ar > SAT / ay) return ((r < 0) != (y < 0)) ? -SAT : SAT;
        r *= y;
        if (ay <= 1 && (r == 1 || r == 0)) return (y == -1 && (k - i - 1) % 2) ? -r : r;
        if (ay <= 1 && r == -1) return ((k - i - 1) % 2) ? 1 : -1;
    }
    return r;
}
static long long bf_mod_pow(long long a, unsigned long long e, long long m) {  // 別の書き方（右から左ではなく左から右）
    Int128 r = 1 % m, x = ((Int128)a % m + m) % m;
    for (int bit = 63; bit >= 0; --bit) {
        r = r * r % m;
        if ((e >> bit) & 1) r = r * x % m;
    }
    return (long long)r;
}

static void test_division(std::mt19937_64& rng, int& cases) {
    auto check = [&](long long a, long long b) {
        if (b == 0 || (a == LMIN && b == -1)) return;
        ++cases;
        MATHB_CHECK(floor_div(a, b) == (long long)bf_floor(a, b));
        MATHB_CHECK(ceil_div(a, b) == (long long)bf_ceil(a, b));
        if (b != LMIN) {
            const long long r = mod_floor(a, b), ab = b < 0 ? -b : b;
            MATHB_CHECK(0 <= r && r < ab && ((Int128)a - r) % ab == 0);
        }
        // 四捨五入：-1/2 ≤ a/b - q < 1/2（ちょうど 1/2 のときは大きい方）
        const long long q = round_div(a, b);
        Int128 x = a, y = b;
        if (y < 0) x = -x, y = -y;
        const Int128 diff2 = 2 * (x - (Int128)q * y);  // 2(a/b - q)·y
        MATHB_CHECK(-y <= diff2 && diff2 < y);
    };
    for (long long a = -40; a <= 40; ++a)
        for (long long b = -12; b <= 12; ++b) check(a, b);
    for (int it = 0; it < 200000 * TEST_SCALE; ++it) check(rnd_ll(rng), rnd_ll(rng));
    MATHB_CHECK(floor_div(-7, 2) == -4 && ceil_div(-7, 2) == -3 && mod_floor(-7, 3) == 2 && round_div(-5, 2) == -2 && round_div(5, 2) == 3);
    MATHB_CHECK(floor_div(LMIN, 1) == LMIN && ceil_div(LMAX, -1) == -LMAX && mod_floor(LMIN, 7) == (long long)(((Int128)LMIN % 7 + 7) % 7));
}

static void test_multiples(std::mt19937_64& rng, int& cases) {
    // 小さい範囲は数える
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long l = rnd(rng, -60, 60), r = rnd(rng, -60, 60), d = rnd(rng, 1, 13), a = rnd(rng, -30, 30);
        long long c1 = 0, c2 = 0, c3 = 0, c4 = 0;
        for (long long x = l; x < r; ++x) c1 += (x % d == 0), c3 += (((x - a) % d) == 0);
        for (long long x = l; x <= r; ++x) c2 += (x % d == 0), c4 += (((x - a) % d) == 0);
        MATHB_CHECK(count_multiples(l, r, d) == c1 && count_multiples_closed(l, r, d) == c2);
        MATHB_CHECK(count_congruent(l, r, a, d) == c3 && count_congruent_closed(l, r, a, d) == c4);
    }
    // 大きい範囲は別の式（⌊(r-1)/d⌋ - ⌊(l-1)/d⌋）と性質で
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        long long l = rnd_ll(rng), r = rnd_ll(rng);
        if (l > r) std::swap(l, r);
        while ((Int128)r - l > ((Int128)1 << 62)) l /= 2, r /= 2;  // 個数が long long に収まる範囲
        long long d = rnd_ll(rng);
        if (d <= 0) d = d == LMIN ? 3 : -d + (d == 0);
        MATHB_CHECK((Int128)count_multiples(l, r, d) == bf_floor((Int128)r - 1, d) - bf_floor((Int128)l - 1, d));
        MATHB_CHECK((Int128)count_multiples_closed(l, r, d) == bf_floor(r, d) - bf_floor((Int128)l - 1, d));
        const long long a = rnd_ll(rng);
        MATHB_CHECK((Int128)count_congruent(l, r, a, d) == bf_floor((Int128)r - 1 - a, d) - bf_floor((Int128)l - 1 - a, d));
        MATHB_CHECK((Int128)count_congruent_closed(l, r, a, d) == bf_floor((Int128)r - a, d) - bf_floor((Int128)l - 1 - a, d));
        // 倍数・合同な値の切り上げ・切り捨て
        const long long x = rnd_ll(rng) / 2;
        d = d % (1LL << 61) + 1;  // 答えが long long に収まる範囲
        const long long cm = ceil_multiple(x, d), fm = floor_multiple(x, d);
        MATHB_CHECK(cm % d == 0 && cm >= x && (Int128)cm - d < x && fm % d == 0 && fm <= x && (Int128)fm + d > x);
        const long long cc = ceil_congruent(x, a, d), fc = floor_congruent(x, a, d);
        MATHB_CHECK(((Int128)cc - a) % d == 0 && cc >= x && (Int128)cc - d < x && ((Int128)fc - a) % d == 0 && fc <= x && (Int128)fc + d > x);
    }
}

static void test_roots(std::mt19937_64& rng, int& cases) {
    auto check = [&](long long n, int k) {
        if (n < 0) return;
        ++cases;
        const long long x = kth_root_floor(n, k), c = kth_root_ceil(n, k);
        MATHB_CHECK(x >= 0 && bf_pow_sat(x, k) <= n && bf_pow_sat((Int128)x + 1, k) > n);
        MATHB_CHECK(c >= 0 && bf_pow_sat(c, k) >= n && (c == 0 || bf_pow_sat((Int128)c - 1, k) < n));
        MATHB_CHECK(is_kth_power(n, k) == (bf_pow_sat(x, k) == n));
        if (k == 2) MATHB_CHECK(is_square(n) == (bf_pow_sat(x, 2) == n) && x == isqrt(n) && c == isqrt_ceil(n));
    };
    for (long long n = 0; n <= 3000; ++n)
        for (int k = 1; k <= 13; ++k) check(n, k);
    for (int it = 0; it < 30000 * TEST_SCALE; ++it) {
        const int k = (int)rnd(rng, 1, 66);
        long long n = rnd_ll(rng);
        if (n < 0) n = n == LMIN ? LMAX : -n;
        check(n, k);
        // 累乗数の前後
        const long long base = kth_root_floor(n, k);
        for (Int128 y = std::max(0LL, base - 1); y <= (Int128)base + 1; ++y) {
            const Int128 p = bf_pow_sat(y, k);
            if (p > LMAX) continue;
            check((long long)p, k);
            if (p > 0) check((long long)p - 1, k);
            if (p < LMAX) check((long long)p + 1, k);
        }
    }
    check(LMAX, 2), check(LMAX, 3), check(LMAX, 62), check(LMAX, 63), check(LMAX, 64);
    MATHB_CHECK(kth_root_floor(1000000000000000000LL, 3) == 1000000 && kth_root_floor(999999999999999999LL, 3) == 999999);
    MATHB_CHECK(kth_root_floor(LMAX, 62) == 2 && kth_root_floor(LMAX, 63) == 1 && kth_root_floor(4611686018427387904LL, 62) == 2);
    MATHB_CHECK(!is_square(-4) && is_square(0) && is_square(1) && !is_square(2) && is_kth_power(0, 5) && !is_kth_power(-8, 3));
}

static void test_pow(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 200000 * TEST_SCALE; ++it, ++cases) {
        long long a = rnd_ll(rng), b = rnd_ll(rng);
        MATHB_CHECK(mul_overflows(a, b) == ((Int128)a * b > LMAX || (Int128)a * b < LMIN));
        MATHB_CHECK(add_overflows(a, b) == ((Int128)a + b > LMAX || (Int128)a + b < LMIN));
        a = a < 0 ? (a == LMIN ? LMAX : -a) : a, b = b < 0 ? (b == LMIN ? LMAX : -b) : b;
        long long cap = rnd_ll(rng);
        cap = cap < 0 ? (cap == LMIN ? 0 : -cap) : cap;
        MATHB_CHECK((Int128)mul_cap(a, b, cap) == std::min<Int128>((Int128)a * b, cap));
        MATHB_CHECK((Int128)add_cap(a, b, cap) == std::min<Int128>((Int128)a + b, cap));
        // 累乗
        const long long base = it % 3 == 0 ? rnd(rng, -12, 12) : (it % 3 == 1 ? rnd(rng, -3100000000LL, 3100000000LL) : rnd_ll(rng));
        const long long e = it % 4 == 0 ? rnd(rng, 0, 70) : (it % 4 == 1 ? rnd(rng, 0, 3) : rnd(rng, 0, LMAX));
        const Int128 p = bf_pow_sat(base, std::min<long long>(e, 200));
        const long long ab = base < 0 ? (base == LMIN ? LMAX : -base) : base;
        const Int128 pa = bf_pow_sat(ab, std::min<long long>(e, 200));
        MATHB_CHECK((Int128)pow_cap(ab, e, cap) == std::min<Int128>(pa, cap));
        if (p >= LMIN && p <= LMAX && e <= 200) MATHB_CHECK(pow_int(base, e) == (long long)p);
    }
    MATHB_CHECK(pow_int(-2, 63) == LMIN && pow_int(2, 62) == 4611686018427387904LL && pow_int(0, 0) == 1 && pow_int(-1, LMAX) == -1);
    MATHB_CHECK(pow_cap(0, 0, 5) == 1 && pow_cap(0, 0, 0) == 0 && pow_cap(2, 1000000000000000000LL, LMAX) == LMAX && pow_cap(1, LMAX, 7) == 1);
    MATHB_CHECK(pow_cap(10, 18, LMAX) == 1000000000000000000LL && pow_cap(10, 19, LMAX) == LMAX);
    for (long long n = 1; n <= 5000; ++n)
        for (long long b = 2; b <= 12; ++b) {
            int f = 0, c = 0;
            for (long long p = 1; p * b <= n; p *= b) ++f;
            for (long long p = 1; p < n; p *= b) ++c;
            MATHB_CHECK(ilog_floor(n, b) == f && ilog_ceil(n, b) == c);
            ++cases;
        }
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        long long n = rnd_ll(rng);
        n = n <= 0 ? (n <= -LMAX ? LMAX : -n + 1) : n;
        const long long b = it % 2 ? rnd(rng, 2, 40) : (rnd_ll(rng) & LMAX) | 2;
        const int f = ilog_floor(n, b), c = ilog_ceil(n, b);
        MATHB_CHECK(bf_pow_sat(b, f) <= n && bf_pow_sat(b, f + 1) > n && bf_pow_sat(b, c) >= n && (c == 0 || bf_pow_sat(b, c - 1) < n));
    }
    for (int k = 0; k <= 18; ++k) MATHB_CHECK((Int128)pow10_ll(k) == bf_pow_sat(10, k));
    MATHB_CHECK(ilog_floor(LMAX, 2) == 62 && ilog_floor(999, 10) == 2 && ilog_floor(1000, 10) == 3 && ilog_ceil(1001, 10) == 4);
}

static void test_gcd(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 200000 * TEST_SCALE; ++it, ++cases) {
        long long a = rnd_ll(rng), b = rnd_ll(rng);
        if (it % 5 == 0) b = a / rnd(rng, 3, 100) * rnd(rng, -3, 3);  // 公約数の多い組
        const unsigned long long ua = a < 0 ? 0ULL - (unsigned long long)a : (unsigned long long)a, ub = b < 0 ? 0ULL - (unsigned long long)b : (unsigned long long)b;
        const unsigned long long g = bf_gcd(ua, ub);
        if (g > (unsigned long long)LMAX) continue;  // gcd = 2^63（LLONG_MIN と 0 など）は収まらない
        const ExtGcd e = ext_gcd(a, b);
        MATHB_CHECK((unsigned long long)e.g == g && (Int128)a * e.x + (Int128)b * e.y == e.g);
        if (g) {
            MATHB_CHECK((Int128)std::llabs(e.x) <= std::max<Int128>(1, (Int128)ub / g) && (Int128)std::llabs(e.y) <= std::max<Int128>(1, (Int128)ua / g));
        } else {
            MATHB_CHECK(e.x == 0 && e.y == 0);
        }
        // lcm
        const long long pa = (long long)(ua & (unsigned long long)LMAX), pb = (long long)(ub & (unsigned long long)LMAX);
        long long cap = rnd_ll(rng);
        cap = cap < 0 ? (cap == LMIN ? 0 : -cap) : cap;
        const unsigned long long gg = bf_gcd(pa, pb);
        const Int128 l = (pa == 0 || pb == 0) ? 0 : (Int128)(pa / gg) * pb;
        MATHB_CHECK((Int128)lcm_cap(pa, pb, cap) == std::min<Int128>(l, cap));
    }
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 8);
        std::vector<long long> v(n);
        for (auto& x : v) x = it % 2 ? rnd(rng, -60, 60) : rnd_ll(rng) / (long long)rnd(rng, 1, 1000000);
        unsigned long long g = 0;
        Int128 l = 1;
        bool zero = false;
        const Int128 SAT = (Int128)1 << 100;
        for (long long x : v) {
            const unsigned long long u = x < 0 ? 0ULL - (unsigned long long)x : (unsigned long long)x;
            g = bf_gcd(g, u);
            if (u == 0) zero = true;
            else if (l < SAT) {
                const Int128 t = l / bf_gcd128(l, (Int128)u);
                l = t > SAT / (Int128)u ? SAT : std::min<Int128>(SAT, t * (Int128)u);
            }
        }
        long long cap = rnd(rng, 0, it % 3 ? 1000000 : LMAX);
        MATHB_CHECK((unsigned long long)gcd_of(v) == g);
        MATHB_CHECK((Int128)lcm_of_cap(v, cap) == (zero ? 0 : std::min<Int128>(l, cap)));
        if (it % 2) {  // vector<int> でも
            std::vector<int> vi;
            for (long long x : v) vi.push_back((int)(x % 1000));
            unsigned long long gi = 0;
            for (int x : vi) gi = bf_gcd(gi, (unsigned long long)(x < 0 ? -x : x));
            MATHB_CHECK((unsigned long long)gcd_of(vi) == gi);
        }
        // try_lcm：2 つの数の lcm を 128 ビットで
        const long long a = it % 3 ? rnd(rng, -100000, 100000) : rnd_ll(rng), b = it % 3 ? rnd(rng, -100000, 100000) : rnd_ll(rng);
        const UInt128 ua = a < 0 ? (UInt128)0 - (UInt128)(Int128)a : (UInt128)a, ub = b < 0 ? (UInt128)0 - (UInt128)(Int128)b : (UInt128)b;
        const Maybe<long long> tl = try_lcm(a, b);
        if (ua == 0 || ub == 0) MATHB_CHECK(tl.ok && tl.value == 0);
        else {
            const UInt128 L = ua / bf_gcd128((Int128)ua, (Int128)ub) * ub;
            MATHB_CHECK(tl.ok == (L <= (UInt128)LMAX) && (!tl.ok || (UInt128)tl.value == L));
        }
    }
    MATHB_CHECK(gcd_of({12, -18, 30}) == 6 && lcm_of_cap({4, 6}, 100) == 12 && gcd_of(std::vector<int>{}) == 0 && lcm_of_cap(std::vector<unsigned>{4u, 6u}, 5) == 5);
    MATHB_CHECK(!try_lcm(LMAX, LMAX - 1).ok && try_lcm(LMIN + 1, 1).value == LMAX && try_lcm(-4, 6).value == 12);
    cases += 2;
}

static void test_mod(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        long long m;
        switch (it % 5) {
            case 0: m = rnd(rng, 1, 30); break;
            case 1: m = rnd(rng, 1, 4294967296LL); break;
            case 2: m = rnd(rng, 4294967297LL, LMAX); break;
            case 3: m = it % 2 ? 998244353 : 1000000007; break;
            default: m = LMAX - rnd(rng, 0, 3); break;
        }
        const long long a = rnd_ll(rng), b = rnd_ll(rng);
        MATHB_CHECK(mod_mul(a, b, m) == (long long)((((Int128)a * b) % m + m) % m));
        const unsigned long long e = it % 3 ? (unsigned long long)rnd(rng, 0, 20) : rng() >> 1;
        const long long p = mod_pow(a, (long long)e, m);
        MATHB_CHECK(p == bf_mod_pow(a, e, m));
        const unsigned long long g = bf_gcd((unsigned long long)mod_floor(a, m), (unsigned long long)m);
        const Maybe<long long> inv = try_mod_inv(a, m);
        MATHB_CHECK(inv.ok == (g == 1));
        if (inv.ok) {
            MATHB_CHECK(0 <= inv.value && inv.value < m && mod_mul(a, inv.value, m) == 1 % m && mod_inv(a, m) == inv.value);
            MATHB_CHECK(mod_mul(mod_pow(a, -(long long)e, m), p, m) == 1 % m);
        }
    }
    MATHB_CHECK(mod_pow(0, 0, 7) == 1 && mod_pow(5, 0, 1) == 0 && mod_pow(-2, 3, 5) == 2 && mod_pow(3, -1, 7) == 5);
    MATHB_CHECK(mod_pow(2, LMAX, LMAX) == bf_mod_pow(2, (unsigned long long)LMAX, LMAX) && mod_inv(LMAX - 1, LMAX) == LMAX - 1);
    // 等比数列の和
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long m = it % 3 == 0 ? rnd(rng, 1, 50) : (it % 3 == 1 ? rnd(rng, 1, 4294967296LL) : rnd(rng, 1, LMAX));
        const long long r = rnd_ll(rng), n = rnd(rng, 0, 300);
        Int128 s = 0, pw = 1 % m;
        const Int128 rr = ((Int128)r % m + m) % m;
        for (long long i = 0; i < n; ++i) s = (s + pw) % m, pw = pw * rr % m;
        MATHB_CHECK(geom_sum_mod(r, n, m) == (long long)s);
        // 大きい n：S(2n) = S(n)(1 + r^n)、S(n+1) = 1 + r S(n)
        const long long N = rnd(rng, 0, LMAX / 2 - 1);
        const long long sn = geom_sum_mod(r, N, m), rn = mod_pow(r, N, m);
        MATHB_CHECK(geom_sum_mod(r, 2 * N, m) == mod_mul(sn, (long long)(((Int128)rn + 1) % m), m));
        MATHB_CHECK(geom_sum_mod(r, N + 1, m) == (long long)(((Int128)mod_mul(r, sn, m) + 1) % m));
    }
    MATHB_CHECK(geom_sum_mod(2, 10, 1000000) == 1023 && geom_sum_mod(1, 1000000000000000000LL, 998244353) == 1000000000000000000LL % 998244353);
}

static void test_sums(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long a = it % 2 ? rnd(rng, -1000, 1000) : rnd_ll(rng), d = it % 2 ? rnd(rng, -1000, 1000) : rnd_ll(rng), n = rnd(rng, 0, 60);
        Int128 s = 0;
        for (long long i = 0; i < n; ++i) s += (Int128)a + (Int128)i * d;
        MATHB_CHECK(arith_sum(a, d, n) == s);
        const long long m = it % 3 ? rnd(rng, 1, 1000000007) : rnd(rng, 1, LMAX);
        MATHB_CHECK(arith_sum_mod(a, d, n, m) == (long long)((s % m + m) % m));
        const long long l = it % 2 ? rnd(rng, -100, 100) : rnd_ll(rng) / 4, r = l + rnd(rng, -3, 80);
        Int128 t = 0, tc = 0;
        for (long long x = l; x < r; ++x) t += x;
        for (long long x = l; x <= r; ++x) tc += x;
        MATHB_CHECK(sum_range(l, r) == t && sum_range_closed(l, r) == tc);
        MATHB_CHECK(sum_range_mod(l, r, m) == (long long)((t % m + m) % m) && sum_range_closed_mod(l, r, m) == (long long)((tc % m + m) % m));
    }
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {  // 端の値：sum_range（Int128）と
        long long l = rnd_ll(rng), r = rnd_ll(rng);
        if (l > r) std::swap(l, r);
        const long long m = it % 2 ? rnd(rng, 1, 1000) : rnd(rng, 1, LMAX);
        if ((Int128)r - l < ((Int128)1 << 62)) {
            const Int128 t = sum_range(l, r);
            MATHB_CHECK(sum_range_mod(l, r, m) == (long long)((t % m + m) % m));
        }
        MATHB_CHECK(sum_range_closed_mod(l, r, m) == (long long)(((Int128)sum_range_mod(l, r, m) + mod_floor(r, m)) % m));
    }
    MATHB_CHECK(sum_range_closed_mod(LMIN, LMAX, 1000000007) == (long long)((-((Int128)1 << 63) % 1000000007 + 1000000007) % 1000000007));
    // 大きい n は式どうしで
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 4000000000000000000LL), a = rnd(rng, -1000000000000000000LL, 1000000000000000000LL), d = rnd(rng, -10, 10);
        const Int128 expect = (Int128)n * a + (Int128)d * ((Int128)n * (n - 1) / 2);
        MATHB_CHECK(arith_sum(a, d, n) == expect);
        const long long m = rnd(rng, 1, LMAX);
        MATHB_CHECK(arith_sum_mod(a, d, n, m) == (long long)((expect % m + m) % m));
        const long long dm = rnd_ll(rng), nm = rnd(rng, 0, LMAX);
        // n(n-1)/2 を 2 つに分けて計算した値と比べる
        const Int128 h = nm % 2 == 0 ? (Int128)(nm / 2) % m * ((nm - (nm > 0)) % m) : (Int128)(nm % m) * (((nm - 1) / 2) % m);
        MATHB_CHECK(arith_sum_mod(0, dm, nm, m) == (long long)(((Int128)mod_floor(dm, m) * (h % m)) % m));
    }
    MATHB_CHECK(sum_range(LMIN, LMAX) == -((Int128)1 << 64) + 1 && sum_range_closed(LMIN, LMAX) == -((Int128)1 << 63));
    MATHB_CHECK(sum_range_closed(1, 100) == 5050 && sum_range(5, 5) == 0 && sum_range_closed(3, 2) == 0 && arith_sum(1, 2, 3) == 9);
}

static void test_frac(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 200000 * TEST_SCALE; ++it, ++cases) {
        long long a = rnd_ll(rng), b = rnd_ll(rng), c = rnd_ll(rng), d = rnd_ll(rng);
        if (b == 0) b = 1;
        if (d == 0) d = -1;
        if (it % 4 == 0) {  // 等しい分数を作る
            a /= 4, b /= 4;
            if (b == 0) b = 1;
            const long long k = rnd(rng, 1, 3) * ((rng() & 1) ? 1 : -1);
            c = a * k, d = b * k;
        }
        // 愚直：符号を b·d の符号で合わせて a·d と c·b を比べる
        const int s = sgn128((Int128)b) * sgn128((Int128)d);
        const Int128 l = (Int128)a * d * s, r = (Int128)c * b * s;
        MATHB_CHECK(frac_compare(a, b, c, d) == (l > r) - (l < r));
        MATHB_CHECK(frac_compare(c, d, a, b) == -frac_compare(a, b, c, d));
    }
    MATHB_CHECK(frac_compare(1, 3, 2, 6) == 0 && frac_compare(-1, 2, 1, -2) == 0 && frac_compare(1, -3, 0, 5) == -1 && frac_compare(LMAX, LMAX - 1, LMAX - 1, LMAX - 2) == -1);
}

static void test_misc(int& cases) {
    const Maybe<long long> m1 = try_mod_inv(4, 8), m2 = try_mod_inv(3, 8);
    MATHB_CHECK(!m1.ok && m2.ok && m2.value == 3);
    auto [ok, v] = try_mod_inv(2, 9);
    MATHB_CHECK(ok && v == 5);
    const Maybe<int> none;
    MATHB_CHECK(!none.ok && none.value == 0);
    ++cases;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_division(rng, cases);
    test_multiples(rng, cases);
    test_roots(rng, cases);
    test_pow(rng, cases);
    test_gcd(rng, cases);
    test_mod(rng, cases);
    test_sums(rng, cases);
    test_frac(rng, cases);
    test_misc(cases);
    std::cout << "[math_base] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[math_base] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 2000000;
    std::vector<long long> xs(n), ys(n);
    for (int i = 0; i < n; ++i) xs[i] = rnd_ll(rng), ys[i] = rnd_ll(rng) | 1, ys[i] = ys[i] == -1 ? 3 : ys[i];  // LLONG_MIN / -1 を避ける
    auto t0 = std::chrono::steady_clock::now();
    unsigned long long s1 = 0;
    for (int i = 0; i < n; ++i) s1 += (unsigned long long)(floor_div(xs[i], ys[i]) + mod_floor(xs[i], ys[i]));
    const long long t_div = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s2 = 0;
    for (int i = 0; i < n / 4; ++i) s2 += (unsigned long long)(mod_pow(xs[i], ys[i] & LMAX, 1000000007));
    const long long t_pow32 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s3 = 0;
    for (int i = 0; i < n / 4; ++i) s3 += (unsigned long long)(mod_pow(xs[i], ys[i] & LMAX, 4611686018427387847LL));
    const long long t_pow64 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s4 = 0;
    for (int i = 0; i < n; ++i) s4 += (unsigned long long)(kth_root_floor(xs[i] & LMAX, 3));
    const long long t_cbrt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s5 = 0;
    for (int i = 0; i < n; ++i) s5 += (unsigned long long)(ext_gcd(xs[i], ys[i]).x);
    const long long t_gcd = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s6 = 0;
    for (int i = 0; i < n / 4; ++i) s6 += (unsigned long long)(geom_sum_mod(xs[i], ys[i] & LMAX, 998244353));
    const long long t_geom = ms(t0);
    t0 = std::chrono::steady_clock::now();
    unsigned long long s7 = 0;
    for (int i = 0; i < n; ++i) s7 += (unsigned long long)(pow_cap(xs[i] & 1023, ys[i] & 63, 1000000000000000000LL));
    const long long t_cap = ms(t0);
    std::cout << "[math_base] benchmark: floor_div + mod_floor × " << n << " " << t_div << " ms, kth_root_floor(·, 3) × " << n << " " << t_cbrt
              << " ms, ext_gcd × " << n << " " << t_gcd << " ms, pow_cap × " << n << " " << t_cap << " ms\n"
              << "[math_base] benchmark: mod_pow（m = 10^9+7）× " << n / 4 << " " << t_pow32 << " ms, mod_pow（m ≈ 4.6×10^18）× " << n / 4 << " "
              << t_pow64 << " ms, geom_sum_mod × " << n / 4 << " " << t_geom << " ms  (check "
              << (s1 + s2 + s3 + s4 + s5 + s6 + s7) % 7 << ")\n";
}

}  // namespace mathb_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace mathb_demo {

{COMMON_DEMO}

static void demo() {
    banner("数学基盤（math_base） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  floor_div / ceil_div / mod_floor / round_div      負の数でも正しい切り捨て・切り上げ・余り・四捨五入\n"
              << "  ceil_multiple / floor_multiple / ceil_congruent / floor_congruent   x 以上・以下で最も近い倍数・合同な数\n"
              << "  count_multiples(_closed) / count_congruent(_closed)                  範囲の中の倍数・合同な数の個数\n"
              << "  isqrt / isqrt_ceil / kth_root_floor / kth_root_ceil / is_square / is_kth_power   累乗根（誤差なし）\n"
              << "  pow_int / pow_cap / mul_cap / add_cap / mul_overflows / add_overflows  累乗と溢れ\n"
              << "  ilog_floor / ilog_ceil / pow10_ll                 対数の切り捨て・切り上げ / 10 の累乗\n"
              << "  ext_gcd / lcm_cap / try_lcm / gcd_of / lcm_of_cap  拡張ユークリッド / 溢れない lcm / 全部の gcd・lcm\n"
              << "  mod_mul / mod_pow / mod_inv / try_mod_inv         64 ビットの mod（合成数の mod でも）\n"
              << "  geom_sum_mod / arith_sum / arith_sum_mod / sum_range(_closed)(_mod)  等比・等差数列の和\n"
              << "  frac_compare / Maybe<T>                           分数の大小（誤差なし）/ ないことがある値\n"
              << "  （Int128・isqrt・parse_fixed は 00_int_common。cout << Int128 もそのまま使える）\n";

    api("【1】 floor_div(a, b) / ceil_div(a, b) / mod_floor(a, m) / round_div(a, b)", "long long",
        "⌊a/b⌋・⌈a/b⌉・0 以上の余り・四捨五入。C++ の / と % は 0 に向かって切り捨てるので、負の数ではこちらを使う", "O(1)");
    SHOW(-7 / 2, "C++ の割り算（0 に向かって切り捨て）");
    SHOW(floor_div(-7, 2), "⌊-3.5⌋ = -4");
    SHOW(ceil_div(-7, 2), "⌈-3.5⌉ = -3");
    SHOW((-7 + 3 - 1) / 3, "よくある切り上げ (a + b - 1) / b は負の数で誤る（⌈-7/3⌉ = -2）");
    SHOW(ceil_div(-7, 3), "正しい切り上げ");
    SHOW(-7 % 3, "C++ の余り（負になる）");
    SHOW(mod_floor(-7, 3), "0 以上の余り");
    SHOW(round_div(7, 2), "3.5 → 4");
    SHOW(round_div(-7, 2), "-3.5 → -3（.5 は大きい方へ）");

    api("【2】 count_multiples(l, r, d) / count_multiples_closed(L, R, d) / count_congruent(l, r, a, m)", "long long",
        "[l, r)（_closed なら [L, R]）にある d の倍数・x ≡ a (mod m) の個数。負の範囲でもよい", "O(1)");
    SHOW(count_multiples_closed(1, 100, 7), "1〜100 の 7 の倍数");
    SHOW(count_multiples(-10, 10, 3), "[-10, 10) の 3 の倍数（-9, -6, -3, 0, 3, 6, 9）");
    SHOW(count_congruent_closed(1, 100, 2, 7), "1〜100 で 7 で割って 2 余る数");
    {
        // ABC131 C: A 以上 B 以下で、C でも D でも割り切れない数
        const long long A = 4, B = 9, C = 2, D = 3;
        const long long cd = lcm_cap(C, D, B + 1);  // lcm が B を超えたら、どちらでも割り切れる数はない
        const long long ans = (B - A + 1) - count_multiples_closed(A, B, C) - count_multiples_closed(A, B, D) + count_multiples_closed(A, B, cd);
        note("ABC131 C（C でも D でも割り切れない数）: 全体 − C の倍数 − D の倍数 + lcm(C, D) の倍数");
        std::cout << "    A = 4, B = 9, C = 2, D = 3 → " << ans << "（5, 7）\n";
    }
    SHOW(ceil_multiple(-7, 3), "-7 以上で最小の 3 の倍数");
    SHOW(floor_multiple(-7, 3), "-7 以下で最大の 3 の倍数");
    SHOW(ceil_congruent(10, 2, 7), "10 以上で 7 で割って 2 余る最小の数");
    SHOW(floor_congruent(10, 2, 7), "10 以下で最大");

    api("【3】 isqrt(n) / kth_root_floor(n, k) / kth_root_ceil(n, k) / is_square(n) / is_kth_power(n, k)", "long long / bool",
        "誤差のない ⌊√n⌋・⌊n^(1/k)⌋・⌈n^(1/k)⌉。std::pow / std::cbrt は 10^18 付近で 1 ずれることがある", "O(k)");
    SHOW(kth_root_floor(1000000000000000000LL, 3), "10^18 の立方根");
    SHOW(kth_root_floor(999999999999999999LL, 3), "10^18 − 1 → 999999（cbrt だと 1000000 になりうる）");
    SHOW(kth_root_ceil(10, 3), "x^3 ≥ 10 となる最小の x");
    SHOW(is_square(144), "12^2");
    SHOW(is_kth_power(1024, 5), "4^5");
    {
        // ABC227 C: A ≤ B ≤ C かつ ABC ≤ N となる組の数
        const long long N = 100;
        long long cnt = 0;
        for (long long a = 1; a <= kth_root_floor(N, 3); ++a)
            for (long long b = a; b <= isqrt(N / a); ++b) cnt += N / (a * b) - b + 1;
        note("ABC227 C（ABC ≤ N の数え上げ）: A ≤ ∛N、B ≤ √(N/A) の範囲で C の個数 ⌊N/(AB)⌋ − B + 1 を足す");
        std::cout << "    N = 100 → " << cnt << "\n";
    }

    api("【4】 pow_cap(a, e, cap) / mul_cap / add_cap / pow_int(a, e) / mul_overflows / add_overflows", "long long / bool",
        "溢れずに「cap を超えるか」を調べる（_cap は min(値, cap)）。pow_int は厳密な値（溢れたら止める）", "O(log e) / O(1)");
    SHOW(pow_cap(10, 30, 1000000000000000001LL), "10^30 は 10^18 + 1 で止まる → 10^18 を超える");
    SHOW(pow_int(-3, 5), "厳密な値");
    SHOW(mul_overflows(3000000000LL, 4000000000LL), "1.2×10^19 は long long に収まらない");
    {
        // ABC169 B: A1 × … × AN。10^18 を超えたら -1（0 があれば 0）
        const std::vector<long long> A = {1000000000LL, 1000000000LL, 10};
        const long long LIM = 1000000000000000000LL;
        long long prod = 1;
        bool zero = false;
        for (long long a : A) zero |= (a == 0);
        for (long long a : A) prod = mul_cap(prod, a, LIM + 1);
        note("ABC169 B（積が 10^18 を超えたら -1）: mul_cap(積, a, 10^18 + 1) で溢れずに掛ける。0 があれば先に 0");
        std::cout << "    {10^9, 10^9, 10} → " << (zero ? 0 : (prod > LIM ? -1 : prod)) << "\n";
    }
    {
        // 典型90 020: log2(a) < b·log2(c) ⇔ a < c^b
        const long long a = 4, b = 3, c = 2;
        note("典型90 020（log2 a < b log2 c）: a < c^b を pow_cap(c, b, a + 1) > a で（実数の log を使わない）");
        std::cout << "    a = 4, b = 3, c = 2 → " << (pow_cap(c, b, a + 1) > a ? "Yes" : "No") << "\n";
    }

    api("【5】 ilog_floor(n, b) / ilog_ceil(n, b) / pow10_ll(k)", "int / long long", "b^k ≤ n となる最大の k / b^k ≥ n となる最小の k / 10^k", "O(log n)");
    SHOW(ilog_floor(999, 10) + 1, "10 進の桁数");
    SHOW(ilog_ceil(1000000, 2), "2^k ≥ 10^6 となる最小の k（ダブリングの段数など）");
    SHOW(pow10_ll(18), "10^18");

    api("【6】 ext_gcd(a, b) / lcm_cap(a, b, cap) / try_lcm(a, b) / gcd_of(v) / lcm_of_cap(v, cap)", "ExtGcd{g, x, y} / long long / Maybe",
        "a x + b y = gcd の解（負の数でも）/ 溢れない lcm / 全部の gcd・lcm（gcd・lcm 1 つずつなら std::gcd / std::lcm）", "O(log)");
    RUN(ExtGcd e = ext_gcd(240, 46));
    SHOW(e.g, "gcd(240, 46)");
    SHOW(240 * e.x + 46 * e.y, "240 x + 46 y = 2 を満たす");
    SHOW(lcm_cap(1000000000000LL, 999999999999LL, 1000000000000000001LL), "lcm ≈ 10^24 は 10^18 + 1 で止まる（典型90 038 の Large 判定）");
    SHOW(lcm_cap(1000000000LL, 999999999LL, 1000000000000000001LL), "超えなければ lcm そのもの");
    SHOW(gcd_of({12, -18, 30}), "全部の gcd（絶対値で。vector<int> などでも）");
    SHOW(lcm_of_cap({4, 6, 10}, 1000), "全部の lcm");
    SHOW(try_lcm(1000000000000LL, 999999999999LL).ok, "lcm が long long に収まらない（Maybe）");
    SHOW(try_lcm(-4, 6).value, "収まれば lcm（絶対値で）");
    note("多くの数の lcm を mod で求めるときは、mod を取りながら lcm を取ってはいけない（02_factorize の lcm_mod_of）。");

    api("【7】 mod_pow(a, e, m) / mod_inv(a, m) / try_mod_inv(a, m) / mod_mul(a, b, m)", "long long / Maybe<long long>",
        "64 ビットの mod（m < 2^63）。逆元は合成数の mod でも gcd(a, m) = 1 なら求まる。mod が 2^31 未満の定数なら 00_modint の方が速い", "O(log e)");
    SHOW(mod_pow(2, 100, 1000000007), "2^100 mod (10^9+7)");
    SHOW(mod_pow(3, -1, 10), "3 の逆元 mod 10（3 × 7 = 21 ≡ 1）");
    SHOW(mod_inv(3, 10), "同じ");
    SHOW(try_mod_inv(4, 10).ok, "gcd(4, 10) = 2 なので逆元はない");
    SHOW(mod_mul(4000000000000000000LL, 3, 9000000000000000000LL), "1.2×10^19 mod 9×10^18（128 ビットで計算）");
    {
        // ABC228 E: M^(K^N) mod p（p = 998244353）。指数は p - 1 で割った余りに。ただし M が p の倍数なら 0
        const long long N = 2, K = 2, M = 2, p = 998244353;
        const long long ans = M % p == 0 ? 0 : mod_pow(M, mod_pow(K, N, p - 1), p);
        note("ABC228 E（M^(K^N) mod p）: フェルマーで指数を p − 1 で割った余りにする。M ≡ 0 のときは答えが 0（指数が 0 になっても 1 にしない）");
        std::cout << "    N = 2, K = 2, M = 2 → 2^4 = " << ans << "\n";
    }

    api("【8】 geom_sum_mod(r, n, m)", "long long", "1 + r + r^2 + … + r^(n-1) を mod m で。逆元を使わないので m が合成数・r - 1 が割れないときも正しい", "O(log n)");
    SHOW(geom_sum_mod(3, 4, 7), "ABC293 E: 1 + 3 + 9 + 27 = 40 ≡ 5 (mod 7)");
    {
        // ABC357 D: N を N 回つなげた数 mod 998244353 = N × (1 + 10^d + 10^2d + …)
        const long long N = 5, p = 998244353;
        const int d = ilog_floor(N, 10) + 1;
        note("ABC357 D（N を N 回つなげた数）: N × geom_sum_mod(10^d, N, p)（d は N の桁数。10^d は mod_pow(10, d, p)）");
        std::cout << "    N = 5 → " << mod_mul(N % p, geom_sum_mod(mod_pow(10, d, p), N, p), p) << "（55555）\n";
    }

    api("【9】 arith_sum(a, d, n) / arith_sum_mod(a, d, n, m) / sum_range(l, r) / sum_range_closed(L, R) / sum_range_mod(l, r, m)", "Int128 / long long",
        "等差数列の和（Int128 に収まらなければ止める。mod で欲しいなら _mod）。sum_range は l 以上 r 未満の整数の和", "O(1)");
    SHOW(arith_sum(1, 2, 5), "1 + 3 + 5 + 7 + 9");
    SHOW(sum_range_closed(1, 1000000000000LL), "1 + … + 10^12（Int128。cout でそのまま出せる）");
    SHOW(arith_sum_mod(1, 1, 1000000000000000000LL, 998244353), "1 + … + 10^18 を mod で");
    SHOW(sum_range_mod(100, 1000000000000000000LL, 998244353), "100 + … + (10^18 − 1) を mod で（引数は区間の端）");

    api("【10】 frac_compare(a, b, c, d)", "int（-1 / 0 / 1）", "a/b と c/d を誤差なしで比べる（分母は負でもよい。128 ビットで掛け算）", "O(1)");
    SHOW(frac_compare(1, 3, 333333333, 1000000000), "1/3 > 0.333333333");
    SHOW(frac_compare(2, -4, -1, 2), "-1/2 = -1/2");
    {
        // ABC308 C: 成功率 A/(A+B) の大きい順（同じなら番号の小さい順）
        std::vector<long long> A = {1, 3, 2}, B = {3, 1, 2};
        std::vector<int> idx = {0, 1, 2};
        std::stable_sort(idx.begin(), idx.end(), [&](int i, int j) { return frac_compare(A[i], A[i] + B[i], A[j], A[j] + B[j]) > 0; });
        note("ABC308 C（成功率で並べる）: stable_sort の比較に frac_compare(A_i, A_i + B_i, A_j, A_j + B_j) > 0");
        std::cout << "    → " << idx[0] + 1 << " " << idx[1] + 1 << " " << idx[2] + 1 << "\n";
    }
    note("何度も分数を足したり掛けたりするなら 08_fraction の Fraction。");

    api("【11】 Maybe<T>{ok, value}", "構造体", "ないことがある値。try_mod_inv などが返す。ok を確かめてから value を使う（構造化束縛でも）", "O(1)");
    RUN(auto [ok, inv] = try_mod_inv(7, 12));
    SHOW(ok, "7 と 12 は互いに素");
    SHOW(inv, "7 × 7 = 49 ≡ 1 (mod 12)");
}

}  // namespace mathb_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    mathb_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mathb_test::benchmark();    // 速度計測
    mathb_demo::demo();         // 使い方デモ
    return 0;
}
