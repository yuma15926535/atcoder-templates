// =============================================================================
//  ランダムテスト
// =============================================================================
namespace meq_test {

{COMMON_TEST}

static const long long LMAX = std::numeric_limits<long long>::max();
static const long long LMIN = std::numeric_limits<long long>::min();
static unsigned long long bf_gcd(unsigned long long a, unsigned long long b) {
    while (b) {
        const unsigned long long t = a % b;
        a = b, b = t;
    }
    return a;
}
static unsigned long long uab(long long a) { return a < 0 ? 0ULL - (unsigned long long)a : (unsigned long long)a; }
static Int128 mod128(Int128 a, Int128 m) {
    a %= m;
    return a < 0 ? a + m : a;
}
static long long no_lmin(long long v) { return v == LMIN ? LMIN + 1 : v; }

static void test_linear_congruence(std::mt19937_64& rng, int& cases) {
    for (long long m = 1; m <= 40; ++m)
        for (int it = 0; it < 300; ++it, ++cases) {
            const long long a = it % 7 == 0 ? rnd_ll(rng) : rnd(rng, -100, 100), b = it % 5 == 0 ? rnd_ll(rng) : rnd(rng, -100, 100);
            const ModSolution s = linear_congruence(a, b, m);
            std::vector<long long> sols, expect;
            for (long long x = 0; x < m; ++x)
                if (mod128((Int128)a * x - b, m) == 0) sols.push_back(x);
            MEQ_CHECK(s.ok == !sols.empty() && !s.overflow);
            if (!s.ok) continue;
            MEQ_CHECK(0 <= s.r && s.r < s.m && m % s.m == 0);
            for (long long x = 0; x < m; ++x)
                if (x % s.m == s.r) expect.push_back(x);
            MEQ_CHECK(sols == expect);
        }
    for (int it = 0; it < 300000 * TEST_SCALE; ++it, ++cases) {
        const long long m = it % 3 == 0 ? rnd(rng, 1, LMAX) : (it % 3 == 1 ? rnd(rng, 1, 1000000000) : std::max<long long>(1, (long long)uab(no_lmin(rnd_ll(rng)))));
        long long a = rnd_ll(rng), b = rnd_ll(rng);
        if (it % 4 == 0) a = a / 1000 * 6, b = b / 1000 * 6;  // 公約数を持たせる
        const ModSolution s = linear_congruence(a, b, m);
        const unsigned long long g = bf_gcd((unsigned long long)mod128(a, m), (unsigned long long)m);
        MEQ_CHECK(s.ok == (mod128(b, m) % (Int128)g == 0));
        if (s.ok) MEQ_CHECK(s.m == (long long)((unsigned long long)m / g) && 0 <= s.r && s.r < s.m && mod128((Int128)a * s.r - b, m) == 0);
    }
}

static void test_diophantine(std::mt19937_64& rng, int& cases) {
    // 小さい値：窓の中の解を全部調べる
    for (int it = 0; it < 25000 * TEST_SCALE; ++it, ++cases) {
        const long long a = rnd(rng, -12, 12), b = rnd(rng, -12, 12), c = rnd(rng, -60, 60);
        const Diophantine d = linear_diophantine(a, b, c);
        if (a == 0 && b == 0) {
            MEQ_CHECK(d.ok == (c == 0) && d.dx == 0 && d.dy == 0);
        } else if (b == 0) {
            MEQ_CHECK(d.ok == (c % a == 0));
            if (d.ok) MEQ_CHECK(d.x0 == c / a && d.dx == 0 && d.dy == -1 && d.y(5) == 5 && d.x(7) == c / a);
        } else {
            std::vector<long long> xs;
            for (long long x = -150; x <= 150; ++x)
                if ((c - a * x) % b == 0) xs.push_back(x);
            MEQ_CHECK(d.ok == !xs.empty());
            if (d.ok) {
                MEQ_CHECK(0 <= d.x0 && d.x0 < d.dx && d.dx == std::llabs(b) / (long long)bf_gcd(uab(a), uab(b)));
                std::vector<long long> got;
                for (long long k = -170; k <= 170; ++k) {  // x0 < 12、dx ≥ 1 なので窓 [−150, 150] を覆う
                    const long long x = d.x(k);
                    if (-150 <= x && x <= 150) {
                        got.push_back(x);
                        MEQ_CHECK(d.y(k) == (c - a * x) / b);
                    }
                }
                MEQ_CHECK(got == xs);
            }
        }
        // 範囲の中の解の数・x が最小の解
        const long long xl = rnd(rng, -25, 25), xr = rnd(rng, -25, 25), yl = rnd(rng, -25, 25), yr = rnd(rng, -25, 25);
        long long cnt = 0, cnt_closed = 0;
        for (long long x = xl; x <= xr; ++x)
            for (long long y = yl; y <= yr; ++y)
                if (a * x + b * y == c) cnt_closed++, cnt += (x < xr && y < yr);
        MEQ_CHECK(diophantine_count(a, b, c, xl, xr, yl, yr) == cnt && diophantine_count_closed(a, b, c, xl, xr, yl, yr) == cnt_closed);
        const DiophantinePoint p = diophantine_min_x(a, b, c, xl);
        long long bx = LMAX;
        for (long long x = xl; x < xl + 200 && bx == LMAX; ++x)
            if (b == 0 ? a * x == c : (c - a * x) % b == 0) bx = x;
        MEQ_CHECK(p.ok == (bx != LMAX));
        if (p.ok) MEQ_CHECK(p.x == bx && a * p.x + b * p.y == c && (b != 0 || p.y == 0));
    }
    // 大きい値：128 ビットで検算
    for (int it = 0; it < 150000 * TEST_SCALE; ++it, ++cases) {
        long long a = no_lmin(rnd_ll(rng)), b = no_lmin(rnd_ll(rng)), c = rnd_ll(rng);
        if (it % 3 == 0) {  // 解を持つように作る
            const long long x = rnd(rng, -(1LL << 26), 1LL << 26), y = rnd(rng, -(1LL << 26), 1LL << 26);  // |a x| < 2^61
            a = rnd(rng, -(1LL << 31), 1LL << 31) * (it % 2 ? 1 : 12), b = rnd(rng, -(1LL << 31), 1LL << 31) * (it % 2 ? 1 : 18);
            c = a * x + b * y;
        }
        const unsigned long long g = bf_gcd(uab(a), uab(b));
        const bool solvable = g == 0 ? c == 0 : (Int128)c % (Int128)g == 0;
        // y0 が long long に入るときだけ linear_diophantine を呼ぶ（入らなければ止まる仕様）
        const bool fits = b == 0 ? !(a == -1 && c == LMIN) : (Int128)uab(c) / uab(b) + (Int128)uab(a) / g < (Int128)LMAX;
        if (fits) {
            const Diophantine d = linear_diophantine(a, b, c);
            MEQ_CHECK(d.ok == solvable);
            if (d.ok && b != 0) {
                MEQ_CHECK((Int128)a * d.x0 + (Int128)b * d.y0 == c && 0 <= d.x0 && d.x0 < d.dx && (unsigned long long)d.dx == uab(b) / g);
                MEQ_CHECK((Int128)d.dy == (Int128)(a / (long long)g) * (b < 0 ? -1 : 1));
                const long long k = rnd(rng, -3, 3);
                const Int128 X = (Int128)d.x0 + (Int128)k * d.dx, Y = (Int128)d.y0 - (Int128)k * d.dy;
                const Int128 L = (Int128)LMAX;
                if (-L <= X && X <= L && -L <= Y && Y <= L) MEQ_CHECK((Int128)a * X + (Int128)b * Y == c);  // 128 ビットに収まるときだけ
            }
        }
        // 範囲の中の解の数：x の範囲が狭いときに 1 つずつ調べる
        const long long xl = rnd_ll(rng), w = rnd(rng, 0, 100);
        const long long xr = xl > LMAX - w ? LMAX : xl + w;
        long long yl = rnd_ll(rng), yr = rnd_ll(rng);
        if (yl > yr) std::swap(yl, yr);
        if (it % 4 == 0) yl = LMIN, yr = LMAX;
        long long cnt = 0;
        long long first = LMAX;
        bool first_found = false;
        if (b == 0) continue;  // b = 0 は小さい値のテストで
        Int128 first_y = 0;
        for (long long x = xl; x < xr; ++x) {
            const Int128 rest = (Int128)c - (Int128)a * x;
            if (rest % b != 0) continue;
            const Int128 y = rest / b;
            if (!first_found) first = x, first_y = y, first_found = true;
            if (yl <= y && y < yr) cnt++;
        }
        MEQ_CHECK(diophantine_count(a, b, c, xl, xr, yl, yr) == cnt);
        if (first_found && uab(b) / g <= 100 && first_y <= (Int128)LMAX && first_y >= (Int128)LMIN) {  // 窓の中に最小の解がある（y が収まるとき）
            const DiophantinePoint p = diophantine_min_x(a, b, c, xl);
            MEQ_CHECK(p.ok && p.x == first && (Int128)p.y == first_y);
        }
    }
    // 端
    MEQ_CHECK(diophantine_count_closed(3, 5, 100, 0, LMAX, 0, LMAX) == 7);  // (0,20) (5,17) … (30,2)
    MEQ_CHECK(diophantine_count_closed(1, 1, LMAX - 1, 0, LMAX, 0, LMAX) == LMAX);
    MEQ_CHECK(diophantine_count(0, 0, 0, 0, 3, 0, 4) == 12 && diophantine_count(0, 0, 1, 0, 3, 0, 4) == 0);
    MEQ_CHECK(!diophantine_min_x(2, 4, 5, 0).ok && diophantine_min_x(1, 0, 7, 8).ok == false && diophantine_min_x(1, 0, 7, 7).x == 7);
    {
        const Diophantine d = linear_diophantine(LMAX, LMAX - 1, 1);  // x = 1, y = −1
        MEQ_CHECK(d.ok && d.x0 == 1 && d.y0 == -1 && d.dx == LMAX - 1 && d.dy == LMAX);
    }
}

static void test_crt(std::mt19937_64& rng, int& cases) {
    // 小さい値：lcm まで全部調べる
    for (int it = 0; it < 40000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 4);
        std::vector<long long> rs(n), ms(n);
        long long L = 1;
        for (int i = 0; i < n; ++i) {
            ms[i] = rnd(rng, 1, 30), rs[i] = it % 2 ? rnd(rng, -100, 100) : rnd_ll(rng);
            L = L / std::gcd(L, ms[i]) * ms[i];
        }
        if (L > 3000) continue;  // 全部調べるのは lcm が小さいときだけ
        long long best = -1;
        for (long long x = 0; x < L && best < 0; ++x) {
            bool ok = true;
            for (int i = 0; i < n && ok; ++i) ok = mod128(x - (Int128)rs[i], ms[i]) == 0;
            if (ok) best = x;
        }
        const ModSolution s = crt_all(rs, ms);
        MEQ_CHECK(!s.overflow && s.ok == (best >= 0));
        if (s.ok) MEQ_CHECK(s.r == best && s.m == L);
        if (n == 2) {
            const ModSolution t = crt_merge(rs[0], ms[0], rs[1], ms[1]);
            MEQ_CHECK(t.ok == s.ok && t.r == s.r && t.m == s.m && !t.overflow);
        }
    }
    // 大きい値：答え x を先に決めて作る
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 6);
        std::vector<long long> rs(n), ms(n);
        const Int128 x = it % 3 == 0 ? (Int128)rnd(rng, 0, LMAX) : (Int128)rnd(rng, 0, 1LL << 40);
        Int128 L = 1;
        bool big = false;
        for (int i = 0; i < n; ++i) {
            const int t = (int)rnd(rng, 0, 3);
            ms[i] = t == 0 ? rnd(rng, 1, LMAX) : (t == 1 ? rnd(rng, 1, 1000000000) : (t == 2 ? rnd(rng, 1, 64) : (1LL << rnd(rng, 0, 62))));
            rs[i] = (long long)mod128(x, ms[i]);  // 余りに mod の倍数を足してもよい
            if (it % 4 == 0 && rs[i] <= LMAX - ms[i]) rs[i] += ms[i];
            if (it % 4 == 1) rs[i] -= ms[i];
            if (!big) {
                const unsigned long long g = bf_gcd((unsigned long long)(L % ms[i]), (unsigned long long)ms[i]);
                L = L / g * ms[i];
                if (L > (Int128)LMAX) big = true;
            }
        }
        const ModSolution s = crt_all(rs, ms);
        MEQ_CHECK(s.overflow == big);
        if (!big) {
            MEQ_CHECK(s.ok && s.m == (long long)L && (Int128)s.r == mod128(x, L));
        } else {  // lcm > 2^63 > x なので x が最小の解
            MEQ_CHECK(s.ok && (Int128)s.r == x && s.m == 0);
        }
        // 1 つの余りをずらす → 解がないか、あっても lcm が収まらないときは最小の解が変わる
        if (n >= 2) {
            std::vector<long long> rs2 = rs;
            const int i = (int)rnd(rng, 0, n - 1);
            if (ms[i] == 1) continue;
            rs2[i] = (long long)mod128((Int128)rs2[i] + 1, ms[i]);
            const ModSolution u = crt_all(rs2, ms);
            // 1 つずつ確かめる：u.ok なら u.r はすべての式を満たす
            if (u.ok) {
                for (int j = 0; j < n; ++j) MEQ_CHECK(mod128((Int128)u.r - rs2[j], ms[j]) == 0);
                if (!u.overflow) MEQ_CHECK(u.r < u.m && u.m == s.m);
            }
        }
    }
    // lcm が収まらず、最小の解も long long を超える
    {
        const long long p = 3037000507LL, q = 3037000537LL, r2 = 4611686014132420609LL;  // p, q は素数（p·q > 2^63）
        const Int128 x = (Int128)LMAX + 12345;
        const ModSolution s = crt_merge((long long)(x % p), p, (long long)(x % q), q);
        MEQ_CHECK(s.overflow && !s.ok);
        const ModSolution t = crt_merge(5, p, 5, q);
        MEQ_CHECK(t.overflow && t.ok && t.r == 5);
        const ModSolution u = crt_all({5, 5, 6}, {p, q, 7});  // 5 ≢ 6 (mod 7) なので long long に入る解はない
        MEQ_CHECK(u.overflow && !u.ok);
        MEQ_CHECK(crt_all({0}, {r2}).ok && crt_all({}, {}).ok && crt_all({}, {}).m == 1);
    }
}

static void test_garner(std::mt19937_64& rng, int& cases) {
    const std::vector<long long> primes = {998244353, 1000000007, 1000000009, 167772161, 469762049, 754974721, 2147483647, 1000000000000000003LL, 999999999999999989LL};
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        std::vector<long long> ms;
        Int128 P = 1;
        // 互いに素な mod：素数から選ぶか、小さい数を互いに素になるように選ぶ
        if (it % 2 == 0) {
            std::vector<long long> pool = primes;
            std::shuffle(pool.begin(), pool.end(), rng);
            for (long long p : pool) {
                if (P > ((Int128)1 << 125) / p) break;
                if (ms.size() >= 4 || rng() % 4 == 0) break;
                ms.push_back(p), P *= p;
            }
        } else {
            for (int k = 0; k < 8; ++k) {
                const long long m = rnd(rng, 1, 1000);
                bool ok = true;
                for (long long q : ms) ok &= std::gcd(q, m) == 1;
                if (ok && P <= ((Int128)1 << 110) / m) ms.push_back(m), P *= m;
            }
        }
        const UInt128 ux = ((UInt128)rng() << 64) | (UInt128)rng();
        const Int128 x = (Int128)(ux % (UInt128)P);
        std::vector<long long> rs;
        for (long long m : ms) rs.push_back((long long)(x % m) - (it % 3 == 0 ? m : 0));
        const long long mod = it % 3 == 0 ? 998244353 : rnd(rng, 1, LMAX);
        MEQ_CHECK(garner(rs, ms, mod) == (long long)(x % mod));
    }
    MEQ_CHECK(garner({}, {}, 7) == 0 && garner({3}, {1}, 10) == 0);
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_linear_congruence(rng, cases);
    test_diophantine(rng, cases);
    test_crt(rng, cases);
    test_garner(rng, cases);
    std::cout << "[mod_equation] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[mod_equation] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 1000000;
    std::vector<long long> xs(n), ys(n);
    for (int i = 0; i < n; ++i) xs[i] = rnd(rng, 1, 1000000000000000000LL), ys[i] = rnd(rng, 1, 1000000000000000000LL);
    auto t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (int i = 0; i < n; ++i) {
        const ModSolution s = linear_congruence(xs[i], ys[i], ys[(i + 1) % n]);
        s1 += s.ok ? s.r % 7 : 1;
    }
    const long long t_lc = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for (int i = 0; i < n; ++i) {
        const Diophantine d = linear_diophantine(xs[i] % 1000000000, ys[i] % 1000000000, xs[(i + 1) % n]);
        s2 += d.ok ? d.x0 % 7 : 1;
    }
    const long long t_dio = ms(t0);
    // crt_all：10^6 個の式（mod は 720720 の約数。答えはすべての式を満たす）
    std::vector<long long> divs;
    for (long long d = 1; d <= 720720; ++d)
        if (720720 % d == 0) divs.push_back(d);
    std::vector<long long> rs(n), mods(n);
    const long long x = 123456;
    for (int i = 0; i < n; ++i) mods[i] = divs[rng() % divs.size()], rs[i] = x % mods[i];
    t0 = std::chrono::steady_clock::now();
    const ModSolution c = crt_all(rs, mods);
    long long s3 = c.r;
    // 互いに素な 2 つ（素数の積は 10^18 程度）× 10^6
    for (int i = 0; i < n; ++i) s3 += crt_merge(xs[i], 1000000007, ys[i], 999999937).r % 7;
    const long long t_crt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s4 = 0;
    const std::vector<long long> gm = {998244353, 1000000007, 1000000009, 167772161, 469762049};
    for (int i = 0; i < n / 5; ++i) s4 += garner({xs[i] % gm[0], xs[i] % gm[1], ys[i] % gm[2], ys[i] % gm[3], xs[i] % gm[4]}, gm, 998244353);
    const long long t_gar = ms(t0);
    std::cout << "[mod_equation] benchmark: linear_congruence × 10^6 " << t_lc << " ms, linear_diophantine × 10^6 " << t_dio << " ms, crt_all（10^6 個）+ crt_merge × 10^6 "
              << t_crt << " ms, garner（5 個）× 2×10^5 " << t_gar << " ms  (check " << (s1 + s2 + s3 + s4) % 7 << ")\n";
}

}  // namespace meq_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace meq_demo {

{COMMON_DEMO}

// 使い方関数：周期的に繰り返す 2 つの区間が重なる最初の時刻（ABC193 E の形）
//   t mod T1 ∈ [l1, r1) かつ t mod T2 ∈ [l2, r2) となる最小の t ≥ 0（なければ −1）。(r1 − l1)(r2 − l2) 回の crt_merge
static long long first_common_time(long long T1, long long l1, long long r1, long long T2, long long l2, long long r2) {
    long long best = -1;
    for (long long a = l1; a < r1; ++a)
        for (long long b = l2; b < r2; ++b) {
            const ModSolution s = crt_merge(a, T1, b, T2);
            if (s.ok && (best < 0 || s.r < best)) best = s.r;
        }
    return best;
}

static void demo() {
    banner("合同式・一次不定方程式・中国剰余定理（mod_equation） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  linear_congruence(a, b, m)                    a x ≡ b (mod m) → ModSolution{ok, r, m}（x ≡ r (mod m)）\n"
              << "  linear_diophantine(a, b, c)                   a x + b y = c → Diophantine{ok, x0, y0, dx, dy}\n"
              << "  diophantine_count / _closed / diophantine_min_x   範囲の中の解の数 / x ≥ x_min で最小の解\n"
              << "  crt_merge(r1, m1, r2, m2) / crt_all(rs, ms)   中国剰余定理（互いに素でなくてよい）\n"
              << "  garner(rs, ms, mod)                           互いに素な mod での答えを別の mod で\n";

    api("【1】 linear_congruence(a, b, m)", "ModSolution{ok, r, m}",
        "a x ≡ b (mod m) の解。x ≡ r (mod m')（m' = m / gcd(a, m)、0 ≤ r < m'）。解がなければ ok = false", "O(log m)");
    RUN(ModSolution s = linear_congruence(3, 6, 10));
    SHOW(s.ok, "解がある");
    SHOW(s.r, "x ≡ 2 (mod 10)");
    SHOW(linear_congruence(4, 6, 10).r, "4x ≡ 6 (mod 10) → x ≡ 4 (mod 5)");
    SHOW(linear_congruence(4, 6, 10).m, "解の mod は 10 / gcd(4, 10) = 5");
    SHOW(linear_congruence(2, 1, 4).ok, "2x ≡ 1 (mod 4) は解なし");
    {
        // ABC186 E（Throne）：S + K x ≡ 0 (mod N) となる最小の x ≥ 0
        note("ABC186 E（Throne）: K x ≡ −S (mod N) を解く。r が答え（解がなければ −1）");
        const long long cases[4][3] = {{10, 4, 3}, {1000, 11, 2}, {998244353, 897581057, 595591169}, {10000, 6, 14}};
        for (auto& t : cases) {
            const ModSolution r = linear_congruence(t[2], -t[1], t[0]);
            std::cout << "    N = " << t[0] << ", S = " << t[1] << ", K = " << t[2] << " → " << (r.ok ? r.r : -1) << "\n";
        }
    }

    api("【2】 linear_diophantine(a, b, c)", "Diophantine{ok, x0, y0, dx, dy}",
        "a x + b y = c の整数解すべて：x = x0 + k·dx、y = y0 − k·dy（d.x(k) / d.y(k)）。0 ≤ x0 < dx（x0 は 0 以上で最小の x）", "O(log)");
    RUN(Diophantine d = linear_diophantine(3, 5, 100));
    SHOW(d.x0, "x = 0 + 5k");
    SHOW(d.y0, "y = 20 − 3k");
    SHOW(d.x(2), "k = 2 なら x = 10");
    SHOW(d.y(2), "y = 14（3·10 + 5·14 = 100）");
    SHOW(linear_diophantine(4, 6, 7).ok, "4x + 6y = 7 は解なし（gcd = 2 が 7 を割らない）");

    api("【3】 diophantine_count(a, b, c, xl, xr, yl, yr) / diophantine_count_closed(…) / diophantine_min_x(a, b, c, x_min)", "long long / DiophantinePoint{ok, x, y}",
        "x, y の範囲の中の解の数（[xl, xr) / 閉区間）。x ≥ x_min で x が最小の解", "O(log)");
    SHOW(diophantine_count_closed(3, 5, 100, 0, LLONG_MAX, 0, LLONG_MAX), "x, y ≥ 0 の解：(0,20) (5,17) … (30,2)");
    SHOW(diophantine_count(1, 1, 10, 0, 5, 0, 100), "x ∈ [0, 5) で x + y = 10");
    RUN(DiophantinePoint p = diophantine_min_x(3, 5, 100, 1));
    SHOW(p.x, "x ≥ 1 で最小の x");
    SHOW(p.y, "そのときの y");
    note("例：a 円玉と b 円玉でちょうど c 円にする組み合わせの数 = diophantine_count_closed(a, b, c, 0, LLONG_MAX, 0, LLONG_MAX)");

    api("【4】 crt_merge(r1, m1, r2, m2) / crt_all(rs, ms)", "ModSolution{ok, r, m, overflow}",
        "x ≡ r_i (mod m_i) をすべて満たす x ≡ r (mod lcm)。r は 0 以上で最小の解。mod は互いに素でなくてよい", "O(n log)");
    SHOW(crt_merge(2, 3, 3, 5).r, "x ≡ 2 (mod 3)、x ≡ 3 (mod 5) → 8");
    SHOW(crt_merge(2, 3, 3, 5).m, "lcm = 15");
    SHOW(crt_all({1, 3, 5}, {4, 6, 8}).ok, "互いに素でなくても、公約数の分で矛盾しなければ解がある");
    SHOW(crt_all({1, 3, 5}, {4, 6, 8}).r, "21 ≡ 1 (mod 4)、≡ 3 (mod 6)、≡ 5 (mod 8)");
    SHOW(crt_all({0, 1}, {4, 6}).ok, "0 と 1 は mod 2 で違うので解なし");
    RUN(ModSolution big = crt_merge(5, 3037000507LL, 5, 3037000537LL));
    SHOW(big.overflow, "lcm が long long を超えた");
    SHOW(big.ok, "それでも long long に入る最小の解はわかる");
    SHOW(big.r, "5");
    note("ABC193 E（Train）：使い方関数 first_common_time（下の例）。t mod (2X+2Y) ∈ [X, X+Y) かつ t mod (P+Q) ∈ [P, P+Q)");
    {
        const long long cases[3][4] = {{5, 2, 7, 6}, {1, 1, 3, 1}, {999999999, 1, 1000000000, 1}};
        for (auto& c : cases) {
            const long long X = c[0], Y = c[1], P = c[2], Q = c[3];
            const long long t = first_common_time(2 * X + 2 * Y, X, X + Y, P + Q, P, P + Q);
            std::cout << "    X = " << X << ", Y = " << Y << ", P = " << P << ", Q = " << Q << " → " << (t < 0 ? std::string("infinity") : std::to_string(t)) << "\n";
        }
    }

    api("【5】 garner(rs, ms, mod)", "long long",
        "互いに素な ms での余りから、0 ≤ x < Π ms の x（Π ms は long long を超えてよい）を mod で", "O(n^2)");
    SHOW(garner({2, 3, 2}, {3, 5, 7}, 1000000007), "x = 23（2 mod 3、3 mod 5、2 mod 7）");
    RUN(const long long X = 123456789012345678LL);
    SHOW(garner({X % 998244353, X % 1000000007, X % 1000000009}, {998244353, 1000000007, 1000000009}, 1000000000000LL), "X mod 10^12 = 789012345678");
    note("いくつかの素数 mod で計算した答えから本当の値（の mod）を戻すとき（任意 mod の畳み込みなど）に使う");
}

}  // namespace meq_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    meq_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    meq_test::benchmark();    // 速度計測
    meq_demo::demo();         // 使い方デモ
    return 0;
}
