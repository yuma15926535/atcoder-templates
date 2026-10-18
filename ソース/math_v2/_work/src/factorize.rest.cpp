// =============================================================================
//  ランダムテスト
// =============================================================================
namespace fac_test {

{COMMON_TEST}

// 愚直な篩（N 以下の素数判定）
static std::vector<char> brute_sieve(int N) {
    std::vector<char> ip(N + 1, 1);
    ip[0] = 0;
    if (N >= 1) ip[1] = 0;
    for (int i = 2; (long long)i * i <= N; ++i)
        if (ip[i])
            for (int j = i * i; j <= N; j += i) ip[j] = 0;
    return ip;
}
static long long prod_of(const std::vector<PrimePower>& fs) {
    Int128 r = 1;
    for (const PrimePower& f : fs)
        for (int k = 0; k < f.e; ++k) r *= f.p;
    return (long long)r;
}
static bool well_formed(const std::vector<PrimePower>& fs) {  // p が昇順・重複なし・すべて素数・e ≥ 1
    for (size_t i = 0; i < fs.size(); ++i) {
        if (fs[i].e < 1 || !is_prime(fs[i].p)) return false;
        if (i && fs[i - 1].p >= fs[i].p) return false;
    }
    return true;
}
static long long rnd_prime(std::mt19937_64& rng, long long lo, long long hi) {
    while (true) {
        const long long x = rnd(rng, lo, hi);
        if (is_prime(x)) return x;
    }
}

static void test_is_prime(std::mt19937_64& rng, int& cases) {
    const int N = 2000000;
    const std::vector<char> ip = brute_sieve(N);
    for (int x = -5; x <= N; ++x, ++cases) FAC_CHECK(is_prime(x) == (x >= 0 && ip[x]));
    // 試し割りで確かめられる大きさ
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const long long x = rnd(rng, 1, 10000000000LL);
        bool pr = x >= 2;
        for (long long d = 2; d * d <= x && pr; ++d) pr = x % d != 0;
        FAC_CHECK(is_prime(x) == pr);
    }
    // 強擬素数・カーマイケル数（合成数）と、知られている大きな素数
    for (long long c : {561LL, 41041LL, 825265LL, 321197185LL, 5394826801LL, 232250619601LL, 9746347772161LL, 3215031751LL, 2152302898747LL,
                        3474749660383LL, 341550071728321LL, 3825123056546413051LL, 318665857834031151LL, 7999252175582851LL,
                        55245642489451LL, 1122004669633LL, 25326001LL, 4759123141LL, 1000000016000000063LL})
        FAC_CHECK(!is_prime(c));
    for (long long p : {2305843009213693951LL, 1000000000000000003LL, 1000000000000000009LL, 998244353LL, 1000000007LL,
                        4611686018427387847LL, 9223372036854775783LL, 4294967311LL, 4294967291LL, 2147483647LL})
        FAC_CHECK(is_prime(p));
    FAC_CHECK(!is_prime(9223372036854775807LL) && !is_prime(-7) && !is_prime(0) && !is_prime(1));
}

static void test_factorize(std::mt19937_64& rng, int& cases) {
    // 小さい数は全部、試し割りと
    for (long long n = 1; n <= 200000; ++n, ++cases) FAC_CHECK(factorize(n) == factorize_trial(n));
    for (int it = 0; it < 1000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 1, it % 4 ? 1000000 : 10000000000LL);
        FAC_CHECK(factorize(n) == factorize_trial(n));
    }
    // 大きい数は「作り方の分かっている数」と
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        std::vector<long long> ps;
        Int128 n = 1;
        const int kind = it % 5;
        if (kind == 0) {  // 2 つの大きい素数の積（ポラード・ローの苦手な形）
            const long long p = rnd_prime(rng, 1000000000LL, 3037000000LL), q = rnd_prime(rng, 2, 9223372036854775807LL / p);
            ps = {p, q};
        } else if (kind == 1) {  // 素数の累乗
            const long long p = rnd_prime(rng, 2, 3000000);
            for (Int128 x = p; x * p <= (Int128)9223372036854775807LL && (ps.size() < 2 || rng() % 3); x *= p) ps.push_back(p);
            if (ps.empty()) ps.push_back(p);
        } else {  // 大きさのばらばらな素数をいくつか
            while (true) {
                const int bits = (int)rnd(rng, 2, 40);
                const long long p = rnd_prime(rng, 2, (1LL << bits));
                if (n * p > (Int128)9223372036854775807LL) break;
                n *= p, ps.push_back(p);
                if (rng() % 4 == 0) break;
            }
            if (ps.empty()) ps.push_back(rnd_prime(rng, 2, 1000000000000000000LL));
        }
        n = 1;
        for (long long p : ps) n *= p;
        std::sort(ps.begin(), ps.end());
        std::vector<PrimePower> expect;
        for (long long p : ps) {
            if (!expect.empty() && expect.back().p == p) ++expect.back().e;
            else expect.push_back({p, 1});
        }
        FAC_CHECK(factorize((long long)n) == expect);
    }
    // ランダムな 64 ビットの数：積が戻ること・素数であること・昇順
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 1, 9223372036854775807LL);
        const std::vector<PrimePower> fs = factorize(n);
        FAC_CHECK(prod_of(fs) == n && well_formed(fs));
    }
    FAC_CHECK(factorize(1).empty() && factorize(9223372036854775807LL) == (std::vector<PrimePower>{{7, 2}, {73, 1}, {127, 1}, {337, 1}, {92737, 1}, {649657, 1}}));
    FAC_CHECK(factorize(1000000016000000063LL) == (std::vector<PrimePower>{{1000000007, 1}, {1000000009, 1}}));
}

static void test_divisor_functions(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 + 300 * TEST_SCALE; ++it, ++cases) {
        const long long n = it < 2000 ? it + 1 : rnd(rng, 1, it % 2 ? 1000000 : 10000000000LL);
        std::vector<long long> ds;
        for (long long d = 1; d * d <= n; ++d)
            if (n % d == 0) {
                ds.push_back(d);
                if (d != n / d) ds.push_back(n / d);
            }
        std::sort(ds.begin(), ds.end());
        FAC_CHECK(divisors_of(n) == ds && count_divisors(n) == (long long)ds.size());
        Int128 s = 0;
        for (long long d : ds) s += d;
        FAC_CHECK(sum_divisors(n) == s);
        // φ・μ・最大の素因数・平方因子・根基を試し割りの分解から
        const std::vector<PrimePower> fs = factorize_trial(n);
        long long phi = n, rad = 1, a = 1, b = 1;
        int mu = 1;
        for (const PrimePower& f : fs) {
            phi = phi / f.p * (f.p - 1), rad *= f.p;
            mu = f.e >= 2 ? 0 : -mu;
            for (int k = 0; k < f.e / 2; ++k) a *= f.p;
            if (f.e % 2) b *= f.p;
        }
        FAC_CHECK(euler_phi(n) == phi && mobius(n) == mu && radical(n) == rad && largest_prime_factor(n) == (fs.empty() ? 1 : fs.back().p));
        const SquareFreeParts sp = square_free_decompose(n);
        FAC_CHECK(sp.a == a && sp.b == b && is_square_free(n) == (a == 1) && is_prime_power(n) == (fs.size() == 1));
        FAC_CHECK(prime_factors(n).size() == fs.size() && omega(n) == (int)fs.size());
        {
            int bo = 0;
            for (const PrimePower& f : fs) bo += f.e;
            FAC_CHECK(big_omega(n) == bo);
        }
        // for_each_divisor は同じ集合
        std::vector<long long> d2;
        for_each_divisor(fs, [&](long long d) { d2.push_back(d); });
        std::sort(d2.begin(), d2.end());
        FAC_CHECK(d2 == ds && divisors_from_factors(fs) == ds);
        if (n <= 100000) {  // φ を数える
            long long c = 0;
            for (long long x = 1; x <= n; ++x) c += std::gcd(x, n) == 1;
            FAC_CHECK(c == phi);
        }
    }
}

static void test_factorial_and_powers(std::mt19937_64& rng, int& cases) {
    // ルジャンドル
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 3000), p = rnd_prime(rng, 2, 60);
        long long c = 0;
        for (long long i = 1; i <= n; ++i)
            for (long long m = i; m % p == 0; m /= p) ++c;
        FAC_CHECK(factorial_exponent(n, p) == c);
        const long long k = rnd(rng, 2, 200);
        long long e = 0;
        while (true) {  // k^(e+1) | n! か（素因数ごとに）
            bool ok = true;
            for (const PrimePower& f : factorize_trial(k)) ok &= factorial_exponent(n, f.p) >= (long long)f.e * (e + 1);
            if (!ok) break;
            ++e;
        }
        FAC_CHECK(factorial_max_power(n, k) == e);
    }
    for (long long k = 1; k <= 3000; ++k, ++cases) {  // k | n! となる最小の n
        long long n = 1;
        Int128 f = 1;
        while ((f % k) != 0) ++n, f = f * n % ((Int128)k * 1000000007);
        FAC_CHECK(min_factorial_multiple(k) == n);
    }
    FAC_CHECK(min_factorial_multiple(1000000007LL * 2) == 1000000007LL && min_factorial_multiple(1LL << 40) == 44);
    FAC_CHECK(factorial_exponent(1000000000000000000LL, 2) == 1000000000000000000LL - __builtin_popcountll(1000000000000000000ULL));
    // 累乗数
    {
        const long long LIM = 1000000;
        std::vector<char> pp(LIM + 1, 0);
        pp[1] = 1;
        for (long long a = 2; a * a <= LIM; ++a)
            for (long long x = a * a; x <= LIM; x *= a) pp[x] = 1;
        long long c = 0;
        for (long long x = 1; x <= LIM; ++x) {
            c += pp[x];
            FAC_CHECK(is_perfect_power(x) == (bool)pp[x]);
            if (x % 997 == 0 || x < 3000) FAC_CHECK(count_perfect_powers(x) == c);
            ++cases;
        }
    }
    // 大きい N はメビウスの式 1 − Σ_{k ≥ 2} μ(k)(⌊N^{1/k}⌋ − 1) と
    for (int it = 0; it < 4 * TEST_SCALE; ++it, ++cases) {
        const long long N = it == 0 ? 1000000000000000000LL : rnd(rng, 1, 9223372036854775807LL >> (rng() % 40));
        long long expect = 1;
        for (int k = 2; k < 63; ++k) expect -= mobius(k) * (kth_root_floor(N, k) - 1);
        FAC_CHECK(count_perfect_powers(N) == expect);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい累乗数
        const long long a = rnd(rng, 2, 3000000), b = rnd(rng, 2, 40);
        const long long x = pow_cap(a, b, 9223372036854775807LL);
        if (x == 9223372036854775807LL) continue;
        FAC_CHECK(is_perfect_power(x));
        if (b % 2 == 1 && !is_square(a)) FAC_CHECK(!is_square(x));
    }
}

static void test_misc(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long N = rnd(rng, 0, 800), M = it % 2 ? rnd(rng, 1, 3000) : rnd(rng, 1, 1000000000000000000LL);
        long long c = 0;
        for (long long x = 1; x <= N; ++x) c += std::gcd(x, M) == 1;
        FAC_CHECK(count_coprime_up_to(N, M) == c);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 6);
        std::vector<long long> vs(n);
        Int128 l = 1;
        bool big = false;
        for (auto& v : vs) {
            v = it % 2 ? rnd(rng, 1, 60) : rnd(rng, 1, 1000000000000LL);
            if (big) continue;
            const Int128 t = l / std::gcd((long long)(l % v), v);
            if (t > ((Int128)1 << 100) / v) big = true;
            else l = t * v;
        }
        const long long m = rnd(rng, 1, 1000000007);
        if (!big) FAC_CHECK(lcm_mod_of(vs, m) == (long long)(l % m));
    }
    FAC_CHECK(lcm_mod_of({0, 5}, 7) == 0 && lcm_mod_of({}, 7) == 1 && lcm_mod_of({4, 6}, 1) == 0 && lcm_mod_of(std::vector<int>{4, 6, 10}, 1000) == 60);
    FAC_CHECK(omega(1) == 0 && big_omega(1) == 0 && omega(720) == 3 && big_omega(720) == 7 && big_omega(1LL << 62) == 62);
    std::ostringstream os;
    for (const PrimePower& f : factorize(360)) os << f << ' ';
    FAC_CHECK(os.str() == "2^3 3^2 5 ");
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_is_prime(rng, cases);
    test_factorize(rng, cases);
    test_divisor_functions(rng, cases);
    test_factorial_and_powers(rng, cases);
    test_misc(rng, cases);
    std::cout << "[factorize] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[factorize] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 100000;
    std::vector<long long> xs(n);
    for (auto& x : xs) x = rnd(rng, 1, 1000000000000000000LL) | 1;
    auto t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (long long x : xs) s1 += is_prime(x);
    const long long t_pr = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for (int i = 0; i < 10000; ++i) s2 += (long long)factorize(xs[i]).size();
    const long long t_fac = ms(t0);
    std::vector<long long> semi(1000);
    for (auto& x : semi) x = rnd_prime(rng, 1000000000LL, 3000000000LL) * rnd_prime(rng, 1000000000LL, 3000000000LL);
    t0 = std::chrono::steady_clock::now();
    long long s3 = 0;
    for (long long x : semi) s3 += factorize(x)[0].p;
    const long long t_semi = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const long long s4 = count_perfect_powers(1000000000000000000LL);
    const long long t_pp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const long long s5 = (long long)divisors_of(963761198400LL).size();
    const long long t_div = ms(t0);
    std::cout << "[factorize] benchmark: is_prime（10^18 程度）× " << n << " " << t_pr << " ms, factorize（ランダム 10^18）× 10^4 " << t_fac
              << " ms, factorize（約 30 ビットの素数 2 つの積）× 1000 " << t_semi << " ms\n"
              << "[factorize] benchmark: count_perfect_powers(10^18) " << t_pp << " ms, divisors_of(963761198400)（約数 6720 個）" << t_div
              << " ms  (check " << (s1 + s2 + s3 + s4 + s5) % 7 << ")\n";
}

}  // namespace fac_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace fac_demo {

{COMMON_DEMO}

static void demo() {
    banner("素因数分解（factorize） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧（n < 2^63）\n"
              << "  is_prime(n) / factorize(n) / prime_factors(n)        素数判定 / 素因数分解（{p, e} の昇順）/ 異なる素因数\n"
              << "  divisors_of(n) / divisors_from_factors(fs) / for_each_divisor(fs, f)   約数\n"
              << "  count_divisors / sum_divisors / euler_phi / mobius / omega / big_omega / largest_prime_factor / radical\n"
              << "  square_free_decompose(n) / is_square_free(n) / is_prime_power(n)          n = a²·b、平方因子なし、素数の累乗\n"
              << "  factorial_exponent(n, p) / factorial_max_power(n, k) / min_factorial_multiple(k)   n! の素因数\n"
              << "  is_perfect_power(n) / count_perfect_powers(N)        a^b の形の数\n"
              << "  count_coprime_up_to(N, M) / lcm_mod_of(vs, m)        M と互いに素な数の個数 / 多くの数の lcm を mod で\n"
              << "  （N 以下の数を何度も分解するなら 02_sieve の Sieve。このファイルの関数も使える）\n";

    api("【1】 is_prime(n) / factorize(n)", "bool / vector<PrimePower>",
        "10^18 の数でも速い（ミラー・ラビンとポラード・ロー）。factorize は p の昇順の {p, e}（PrimePower）。cout << で p^e と出る", "O(log n) / 期待 O(n^{1/4})");
    SHOW(is_prime(1000000007), "素数");
    SHOW(is_prime(1000000016000000063LL), "10^9+7 と 10^9+9 の積");
    SHOW(factorize(360), "2^3 · 3^2 · 5");
    SHOW(factorize(1000000016000000063LL), "大きい素数 2 つの積も速い");
    RUN(for (auto [p, e] : factorize(720)) std::cout << "  " << p << "^" << e; std::cout << "\n");
    {
        // ABC169 D: N = p1^e1 ... に、z = p^k（使っていないもの）で割る操作を何回できるか。各 e について 1 + 2 + … + k ≤ e となる最大の k
        const long long N = 24;
        long long ans = 0;
        for (auto [p, e] : factorize(N))
            for (long long k = 1, rest = e; rest >= k; rest -= k, ++k) ++ans;
        note("ABC169 D（Div Game）: 各素因数の指数 e について 1 + 2 + … ≤ e となる回数を足す");
        std::cout << "    N = 24 → " << ans << "\n";
    }
    {
        // ABC284 D: N = p^2 q（N ≤ 9×10^18）の p と q
        const long long N = 2023;
        auto fs = factorize(N);
        const long long p = fs[0].e == 2 ? fs[0].p : fs[1].p, q = fs[0].e == 2 ? fs[1].p : fs[0].p;
        note("ABC284 D（N = p²q）: factorize(N) で分解して指数 2 の方が p");
        std::cout << "    N = 2023 → p = " << p << ", q = " << q << "\n";
    }

    api("【2】 divisors_of(n) / count_divisors(n) / sum_divisors(n)", "vector<long long> / long long / Int128",
        "約数（昇順）・約数の個数・約数の和。素因数分解が分かっているなら divisors_from_factors(fs)", "O(約数の個数 · log)");
    SHOW(divisors_of(36), "昇順");
    SHOW(count_divisors(963761198400LL), "10^12 以下で約数が最も多い数の 1 つ");
    SHOW(sum_divisors(28), "完全数（28 + 28）");
    {
        // ABC190 D: 和が N となる連続する整数の列の数 = 2 × (N の奇数の約数の個数)
        long long N = 12, odd = N;
        while (odd % 2 == 0) odd /= 2;
        note("ABC190 D（和が N の等差数列の数）: 2 × count_divisors(N から 2 を除いた数)");
        std::cout << "    N = 12 → " << 2 * count_divisors(odd) << "\n";
    }

    api("【3】 euler_phi(n) / mobius(n) / omega(n) / big_omega(n) / largest_prime_factor(n) / radical(n) / is_prime_power(n)", "long long / int / bool",
        "φ(n)（n 以下で n と互いに素な数の個数）・μ(n)・最大の素因数・異なる素因数の積・素数の累乗か", "factorize と同じ");
    SHOW(euler_phi(36), "36 と互いに素な 1〜36 の数");
    SHOW(mobius(30), "2·3·5 → (-1)^3");
    SHOW(omega(720), "720 = 2^4·3^2·5：異なる素因数は 3 個");
    SHOW(big_omega(720), "重複を含めると 4 + 2 + 1 = 7 個");
    SHOW(mobius(12), "平方因子があれば 0");
    SHOW(largest_prime_factor(1001), "7·11·13");
    SHOW(is_prime_power(243), "3^5");

    api("【4】 square_free_decompose(n) / is_square_free(n)", "SquareFreeParts{a, b}", "n = a²·b（b は平方因子を持たない）。「掛けて平方数になる組」は b が同じ組", "factorize と同じ");
    RUN(auto [a, b] = square_free_decompose(72));
    SHOW(a, "72 = 6² · 2");
    SHOW(b, "平方因子を除いた部分");
    note("ABC254 D（i·j が平方数になる組）: i·j が平方数 ⇔ i と j の b が同じ（N 以下を全部なら 02_sieve の Sieve::square_free_decompose）");

    api("【5】 factorial_exponent(n, p) / factorial_max_power(n, k) / min_factorial_multiple(k)", "long long",
        "n! が素数 p で何回割れるか（ルジャンドル）/ k^e | n! の最大の e / k | n! となる最小の n", "O(log n) / factorize(k)");
    SHOW(factorial_exponent(100, 5), "100! の末尾の 0 の数（5 の方が少ない）");
    SHOW(factorial_max_power(10, 12), "12^e | 10! の最大の e（2^8 · 3^4 → 4）");
    SHOW(min_factorial_multiple(30), "ABC280 D: 30 | 5!");
    SHOW(min_factorial_multiple(123456789011LL), "大きい素数なら k そのもの");
    note("ABC052 C（N! の約数の個数）: 各素数 p ≤ N について factorial_exponent(N, p) + 1 を掛ける");

    api("【6】 is_perfect_power(n) / count_perfect_powers(N)", "bool / long long", "a^b（a ≥ 1, b ≥ 2）の形か / 1〜N にいくつあるか（1 も数える）", "O(N^{1/3} log)");
    SHOW(is_perfect_power(1000000), "10^6 = 1000^2 = 100^3");
    SHOW(count_perfect_powers(99), "ABC361 F: 1, 4, 8, 9, 16, 25, 27, 32, 36, 49, 64, 81");
    SHOW(count_perfect_powers(1000000000000000000LL), "N = 10^18 でも速い");

    api("【7】 count_coprime_up_to(N, M) / lcm_mod_of(vs, m)", "long long", "1〜N で M と互いに素な数の個数（包除）/ 多くの数の lcm を mod で", "O(2^(M の素因数の数)) / 分解");
    SHOW(count_coprime_up_to(100, 30), "2・3・5 のどれでも割り切れない数");
    SHOW(lcm_mod_of({1000000007, 1000000009, 6}, 998244353), "lcm は約 6×10^18。mod を取りながら lcm を取ると誤る");
    note("N 以下の数をたくさん分解するときは 02_sieve の Sieve S(N)（このファイルの関数も一緒に入っている）。");
}

}  // namespace fac_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    fac_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    fac_test::benchmark();    // 速度計測
    fac_demo::demo();         // 使い方デモ
    return 0;
}
