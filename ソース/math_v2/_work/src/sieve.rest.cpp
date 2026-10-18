// =============================================================================
//  ランダムテスト
// =============================================================================
namespace sieve_test {

{COMMON_TEST}

static std::vector<char> brute_sieve(int N) {
    std::vector<char> ip(N + 1, 1);
    ip[0] = 0;
    if (N >= 1) ip[1] = 0;
    for (int i = 2; (long long)i * i <= N; ++i)
        if (ip[i])
            for (int j = i * i; j <= N; j += i) ip[j] = 0;
    return ip;
}

static void check_members(const Sieve& S, long long x) {
    const std::vector<PrimePower> fs = factorize_trial(x);
    SIEVE_CHECK(S.factorize(x) == fs);
    long long cd = 1, phi = x, a = 1, b = 1, lp = 1;
    Int128 sd = 1;
    int mu = 1, big = 0;
    for (const PrimePower& f : fs) {
        cd *= f.e + 1, phi = phi / f.p * (f.p - 1), big += f.e, lp = f.p;
        Int128 t = 1, pk = 1;
        for (int k = 0; k < f.e; ++k) pk *= f.p, t += pk;
        sd *= t;
        mu = f.e >= 2 ? 0 : -mu;
        for (int k = 0; k < f.e / 2; ++k) a *= f.p;
        if (f.e % 2) b *= f.p;
    }
    SIEVE_CHECK(S.count_divisors(x) == cd && S.sum_divisors(x) == sd && S.phi(x) == phi && S.mu(x) == mu);
    SIEVE_CHECK(S.omega(x) == (int)fs.size() && S.big_omega(x) == big && S.largest_prime_factor(x) == lp);
    const SquareFreeParts sp = S.square_free_decompose(x);
    SIEVE_CHECK(sp.a == a && sp.b == b && S.is_square_free(x) == (a == 1) && S.prime_factors(x).size() == fs.size());
    SIEVE_CHECK(S.is_prime(x) == (fs.size() == 1 && fs[0].e == 1));
}

static void test_sieve(std::mt19937_64& rng, int& cases) {
    for (int N : {1, 2, 3, 4, 10, 97, 100, 1000, 65537, 200000}) {
        const Sieve S(N);
        const std::vector<char> ip = brute_sieve(N);
        std::vector<long long> ps;
        for (int x = 0; x <= N; ++x)
            if (ip[x]) ps.push_back(x);
        SIEVE_CHECK(S.primes() == ps && S.size() == N);
        long long cnt = 0;
        for (int x = 0; x <= N; ++x, ++cases) {
            cnt += ip[x];
            SIEVE_CHECK(S.is_prime(x) == (bool)ip[x] && S.prime_count(x) == cnt);
            if (x >= 2) SIEVE_CHECK(S.spf(x) == factorize_trial(x)[0].p);
            if (x >= 1 && (N <= 1000 || x % 7 == 0)) check_members(S, x);
        }
        for (int k = 0; k < (int)ps.size(); ++k) SIEVE_CHECK(S.kth_prime(k) == ps[k]);
        for (int it = 0; it < 200; ++it, ++cases) {
            const long long l = rnd(rng, -3, N + 1), r = rnd(rng, -3, N + 1);
            long long c = 0;
            for (long long x = std::max(0LL, l); x < r; ++x) c += ip[x];
            SIEVE_CHECK(S.count_primes(l, r) == c && S.count_primes_closed(l, r - 1) == c);
        }
        // N を超える数は 02_factorize に任せる
        for (int it = 0; it < 300; ++it, ++cases) {
            const long long x = rnd(rng, N + 1, it % 2 ? 10000000000LL : 1000000000000000000LL);
            SIEVE_CHECK(S.factorize(x) == factorize(x) && S.is_prime(x) == is_prime(x));
            if (x <= 10000000000LL) check_members(S, x);
        }
    }
    // lcm_mod
    const Sieve S(100000);
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        std::vector<long long> vs((size_t)rnd(rng, 0, 6));
        for (auto& v : vs) v = rnd(rng, 1, it % 2 ? 100000 : 1000000000000LL);
        const long long m = rnd(rng, 1, 1000000007);
        SIEVE_CHECK(S.lcm_mod(vs, m) == lcm_mod_of(vs, m));
    }
    SIEVE_CHECK(S.lcm_mod(std::vector<int>{4, 6, 10}, 1000) == 60 && S.lcm_mod({8, 12}, 1000) == 24 && S.count_primes_closed(2, 2) == 1 && S.count_primes_closed(5, 4) == 0);
    {  // long long の N でも作れる・primes() は long long
        const long long NL = 1000;
        const Sieve S2(NL);
        long long sq = 0;
        for (auto p : S2.primes()) sq += p * p;
        SIEVE_CHECK(S2.size() == 1000 && sq > 0 && std::is_same<std::decay_t<decltype(S2.primes()[0])>, long long>::value);
    }
}

static void test_tables(std::mt19937_64& rng, int& cases) {
    for (int N : {0, 1, 2, 3, 10, 1000, 100000}) {
        const Sieve S(std::max(N, 1));
        const std::vector<long long> phi = phi_table(N), ds = divisor_sum_table(N);
        const std::vector<int> mu = mu_table(N), dc = divisor_count_table(N), om = omega_table(N);
        SIEVE_CHECK((int)phi.size() == N + 1 && (int)mu.size() == N + 1 && (int)ds.size() == N + 1 && (int)dc.size() == N + 1 && (int)om.size() == N + 1);
        SIEVE_CHECK(phi[0] == 0 && mu[0] == 0 && dc[0] == 0 && ds[0] == 0 && om[0] == 0);
        for (int x = 1; x <= N; ++x, ++cases)
            SIEVE_CHECK(phi[x] == S.phi(x) && mu[x] == S.mu(x) && dc[x] == S.count_divisors(x) && ds[x] == (long long)S.sum_divisors(x) && om[x] == S.omega(x));
    }
    // primes_up_to（区間の境目の前後も）
    const int M = 3000000;
    const std::vector<char> ip = brute_sieve(M);
    std::vector<int> all;
    for (int x = 0; x <= M; ++x)
        if (ip[x]) all.push_back(x);
    std::vector<long long> Ns = {0, 1, 2, 3, 4, 5, 9, 25, 131073, 131074, 131075, 131076, 262147, 262149, M};
    for (int it = 0; it < 20; ++it) Ns.push_back(rnd(rng, 0, M));
    for (long long N : Ns) {
        const std::vector<int> ps = primes_up_to<int>(N);
        const std::vector<long long> pl = primes_up_to(N);
        size_t k = 0;
        while (k < all.size() && all[k] <= N) ++k;
        SIEVE_CHECK(ps == std::vector<int>(all.begin(), all.begin() + (long)k) && pl == std::vector<long long>(all.begin(), all.begin() + (long)k));
        ++cases;
    }
    SIEVE_CHECK(primes_up_to(20000000).size() == 1270607 && primes_up_to<int>(20000000).size() == 1270607);
    SIEVE_CHECK(phi_table(3LL).size() == 4 && omega_table(0).size() == 1);
    {  // 大きめの N で、約数を全部たどる方法と比べる
        const int N = 2000000;
        std::vector<int> dc(N + 1, 0), om(N + 1, 0);
        std::vector<long long> ds(N + 1, 0);
        for (int i = 1; i <= N; ++i)
            for (int j = i; j <= N; j += i) ++dc[j], ds[j] += i;
        for (int p = 2; p <= N; ++p)
            if (om[p] == 0 && dc[p] == 2)
                for (int j = p; j <= N; j += p) ++om[j];
        SIEVE_CHECK(divisor_count_table(N) == dc && divisor_sum_table(N) == ds && omega_table(N) == om);
        cases += 3;
    }
}

static void test_segment(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 60 * TEST_SCALE; ++it) {
        long long l;  // √r までの篩を毎回作るので、大きい l（10^14 付近）ははじめの 24 回だけ
        switch (it < 24 ? it % 8 : it % 6) {
            case 0: case 1: case 2: l = rnd(rng, -5, 100); break;
            case 3: case 4: case 5: l = rnd(rng, 1, 1000000000LL); break;
            case 6: l = rnd(rng, 1, 100000000000000LL); break;
            default: l = 100000000000000LL - rnd(rng, 0, 3000); break;
        }
        const long long r = l + rnd(rng, 0, 2000);
        std::vector<long long> expect;
        for (long long x = std::max(l, 0LL); x < r; ++x)
            if (is_prime(x)) expect.push_back(x);
        SIEVE_CHECK(segment_primes(l, r) == expect && segment_primes_closed(l, r - 1) == expect);
        long long next = std::max(l, 1LL);
        segment_factorize(l, r, [&](long long x, const std::vector<PrimePower>& fs) {
            SIEVE_CHECK(x == next && fs == factorize(x));
            ++next, ++cases;
        });
        SIEVE_CHECK(next == std::max(r, std::max(l, 1LL)));
    }
    // ブロック（2^18 個）をまたぐ
    long long next = 1;
    segment_factorize_closed(1, 600000, [&](long long x, const std::vector<PrimePower>& fs) {
        SIEVE_CHECK(x == next && fs == factorize(x));
        ++next, ++cases;
    });
    SIEVE_CHECK(next == 600001);
    next = 99999999700000LL;  // 10^14 の手前（√r ≈ 10^7）
    segment_factorize(99999999700000LL, 100000000000000LL, [&](long long x, const std::vector<PrimePower>& fs) {
        SIEVE_CHECK(x == next++);
        if (x % 97 == 0) SIEVE_CHECK(fs == factorize(x));
        Int128 p = 1;
        for (const PrimePower& f : fs) {
            for (int k = 0; k < f.e; ++k) p *= f.p;
        }
        SIEVE_CHECK(p == x);
        ++cases;
    });
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_sieve(rng, cases);
    test_tables(rng, cases);
    test_segment(rng, cases);
    std::cout << "[sieve] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[sieve] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    const Sieve S(10000000);
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (int x = 1; x <= 10000000; ++x) s1 += (long long)S.factorize(x).size();
    const long long t_fac = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<long long> phi = phi_table(10000000);
    const long long t_phi = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<int> om = omega_table(10000000);
    const long long t_om = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const size_t np = primes_up_to(100000000).size();
    const long long t_pr = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const size_t ns = segment_primes(100000000000000LL, 100000000000000LL + 10000000).size();
    const long long t_seg = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    segment_factorize(100000000000000LL, 100000000000000LL + 1000000, [&](long long, const std::vector<PrimePower>& fs) { s2 += (long long)fs.size(); });
    const long long t_segf = ms(t0);
    std::cout << "[sieve] benchmark: Sieve(10^7) " << t_build << " ms, S.factorize を 1〜10^7 全部 " << t_fac << " ms, phi_table(10^7) " << t_phi
              << " ms, omega_table(10^7) " << t_om << " ms\n"
              << "[sieve] benchmark: primes_up_to(10^8)（" << np << " 個）" << t_pr << " ms, segment_primes（10^14 から幅 10^7。" << ns << " 個）" << t_seg
              << " ms, segment_factorize（10^14 から幅 10^6）" << t_segf << " ms  (check " << (s1 + s2 + phi[12345] + om[30]) % 7 << ")\n";
}

}  // namespace sieve_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace sieve_demo {

{COMMON_DEMO}

static void demo() {
    banner("篩（sieve） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Sieve S(N)                          N 以下の最小の素因数の表（N ≤ 5×10^7 程度）\n"
              << "    is_prime / spf / primes / factorize / prime_factors / divisors / count_divisors / sum_divisors\n"
              << "    phi / mu / omega / big_omega / largest_prime_factor / square_free_decompose / is_square_free\n"
              << "    prime_count(x) / count_primes(l, r) / count_primes_closed(L, R) / kth_prime(k) / lcm_mod(vs, m)   （N を超える数は 02_factorize に任せる）\n"
              << "  phi_table / mu_table / divisor_count_table / divisor_sum_table / omega_table   0〜N の表（N ≤ 10^8）\n"
              << "  primes_up_to(N) / primes_up_to<int>(N)   N 以下の素数 vector<long long> / メモリ半分の vector<int>（N < 2^31）\n"
              << "  segment_primes(l, r) / segment_factorize(l, r, f)（_closed も）   区間篩（r ≤ 10^14）\n";

    api("【1】 Sieve S(N);  S.is_prime(x) / S.factorize(x) / S.divisors(x)", "bool / vector<PrimePower> / vector<long long>",
        "最小の素因数の表を作り、N 以下の数の素数判定・素因数分解を O(log x) で。N を超える数を渡してもよい（ポラード・ロー）", "構築 O(N)、1 回 O(log x)");
    RUN(Sieve S(1000000));
    SHOW(S.is_prime(999983), "100 万以下で最大の素数");
    SHOW(S.factorize(720720), "2^4 · 3^2 · 5 · 7 · 11 · 13");
    SHOW(S.divisors(60), "昇順");
    SHOW(S.phi(36), "φ(36)");
    SHOW(S.prime_count(1000000), "100 万以下の素数の数");
    SHOW(S.count_primes_closed(100, 200), "[100, 200] の素数の数（count_primes(l, r) は [l, r)）");
    SHOW(S.kth_prime(0), "0 番目の素数");
    SHOW(S.factorize(1000000016000000063LL), "N を超えたら 02_factorize（ポラード・ロー）");
    {
        // ABC177 E: すべての組が互いに素（pairwise）/ 全体の gcd が 1（setwise）/ どちらでもない
        const std::vector<int> A = {6, 10, 15};
        std::vector<int> seen(1000001, 0);
        bool pairwise = true;
        for (int a : A)
            for (long long p : S.prime_factors(a)) pairwise &= ++seen[p] == 1;
        long long g = 0;
        for (int a : A) g = std::gcd(g, (long long)a);
        note("ABC177 E（Coprime）: 各素因数が 2 つ以上の数に出てこなければ pairwise。そうでなく全体の gcd が 1 なら setwise");
        std::cout << "    {6, 10, 15} → " << (pairwise ? "pairwise coprime" : g == 1 ? "setwise coprime" : "not coprime") << "\n";
    }
    {
        // ABC254 D: i·j が平方数になる 1 ≤ i, j ≤ N の組の数 = Σ_b (平方因子を除いた部分が b の数)^2
        const int N = 4;
        std::map<long long, long long> cnt;
        for (int i = 1; i <= N; ++i) cnt[S.square_free_decompose(i).b]++;
        long long ans = 0;
        for (auto& kv : cnt) ans += kv.second * kv.second;
        note("ABC254 D（Together Square）: 平方因子を除いた部分 b が同じ組なら積が平方数");
        std::cout << "    N = 4 → " << ans << "\n";
    }

    api("【2】 phi_table(N) / mu_table(N) / divisor_count_table(N) / divisor_sum_table(N) / omega_table(N)", "vector",
        "0〜N の値の表（添字がそのまま x）。たくさんの x の φ・μ などが要るとき", "O(N)（線形篩）");
    RUN(auto om = omega_table(15));
    SHOW(om, "異なる素因数の数");
    {
        // 典型90 030: N 以下で異なる素因数を K 個以上持つ数の個数
        const int N = 15, K = 2;
        int c = 0;
        for (int x = 1; x <= N; ++x) c += om[x] >= K;
        note("典型90 030（K Factors）: omega_table(N) で数える");
        std::cout << "    N = 15, K = 2 → " << c << "（6, 10, 12, 14, 15）\n";
    }
    SHOW(phi_table(10), "φ(0)〜φ(10)");
    SHOW(mu_table(10), "μ(0)〜μ(10)");

    api("【3】 primes_up_to(N) / primes_up_to<int>(N)", "vector<long long> / vector<int>",
        "N 以下の素数（昇順）。Sieve より速く省メモリ（10^8 で約 0.3 秒・素数 576 万個）。long long なので p * p も溢れない。<int> はメモリ半分", "O(N log log N)");
    RUN(auto ps = primes_up_to(30));
    SHOW(ps, "30 以下の素数");
    {
        // ABC383 D: N 以下で約数がちょうど 9 個の数 = p^8 または p^2 q^2（p < q）
        const long long N = 200;
        const std::vector<long long> qs = primes_up_to(isqrt(N));
        long long c = 0;
        for (size_t i = 0; i < qs.size(); ++i) {
            if (pow_cap(qs[i], 8, N + 1) <= N) ++c;
            for (size_t j = i + 1; j < qs.size() && qs[i] * qs[j] <= isqrt(N); ++j) ++c;
        }
        note("ABC383 D（9 Divisors）: p^8 か p²q²。√N 以下の素数で数える");
        std::cout << "    N = 200 → " << c << "（36, 100, 196）\n";
    }

    api("【4】 segment_primes(l, r) / segment_factorize(l, r, f)（_closed も）", "vector<long long> / なし",
        "[l, r) の素数 / 各数の素因数分解（r ≤ 10^14、幅 ≤ 10^7 程度。√r 以下の素数だけで割る）", "O((r − l) log log r + √r)");
    SHOW(segment_primes_closed(1000000000000LL, 1000000000100LL), "10^12 〜 10^12 + 100 の素数");
    {
        // ABC412 E: A_n = lcm(1, …, n)。A_L, …, A_R の値の種類 = 1 + (L, R] の素数の累乗の数
        const long long L = 4, R = 12;
        long long c = 1;
        segment_factorize_closed(L + 1, R, [&](long long, const std::vector<PrimePower>& fs) { c += fs.size() == 1; });
        note("ABC412 E（LCM Sequence）: lcm(1..n) は n が素数の累乗のときだけ変わる。(L, R] の素数の累乗を segment_factorize で数える");
        std::cout << "    L = 4, R = 12 → " << c << "（5, 7, 8, 9, 11 で変わる）\n";
    }
    note("x の小さい順に f(x, fs) が呼ばれる。fs は p の昇順の vector<PrimePower>（1 なら空）。");
}

}  // namespace sieve_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    sieve_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    sieve_test::benchmark();    // 速度計測
    sieve_demo::demo();         // 使い方デモ
    return 0;
}
