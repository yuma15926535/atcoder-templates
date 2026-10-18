// =============================================================================
//  ランダムテスト
// =============================================================================
namespace divt_test {

{COMMON_TEST}

static long long bgcd(long long a, long long b) { return b ? bgcd(b, a % b) : a; }

static void test_transforms(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int N = (int)rnd(rng, 0, it % 2 ? 40 : 300);
        std::vector<long long> f((size_t)N + 1), g((size_t)N + 1);
        for (auto& x : f) x = rnd(rng, -1000, 1000);
        for (auto& x : g) x = rnd(rng, -1000, 1000);
        std::vector<long long> dz((size_t)N + 1, 0), mz((size_t)N + 1, 0);
        dz[0] = f[0], mz[0] = f[0];
        for (int n = 1; n <= N; ++n)
            for (int d = 1; d <= N; ++d) {
                if (n % d == 0) dz[(size_t)n] += f[(size_t)d];
                if (d % n == 0) mz[(size_t)n] += f[(size_t)d];
            }
        std::vector<long long> a = f, b = f;
        divisor_zeta(a), multiple_zeta(b);
        DIVT_CHECK(a == dz && b == mz);
        divisor_mobius(a), multiple_mobius(b);
        DIVT_CHECK(a == f && b == f);
        // gcd・lcm の畳み込み
        if (N <= 80) {
            std::vector<long long> cg((size_t)N + 1, 0), cl((size_t)N + 1, 0);
            for (int i = 1; i <= N; ++i)
                for (int j = 1; j <= N; ++j) {
                    const long long gg = bgcd(i, j), l = (long long)i / gg * j;
                    cg[(size_t)gg] += f[(size_t)i] * g[(size_t)j];
                    if (l <= N) cl[(size_t)l] += f[(size_t)i] * g[(size_t)j];
                }
            DIVT_CHECK(gcd_convolution(f, g) == cg && lcm_convolution(f, g) == cl);
            // gcd ごとの組の数
            std::vector<long long> cnt((size_t)N + 1, 0), e((size_t)N + 1, 0);
            std::vector<long long> vals;
            for (int t = 0; t < (int)rnd(rng, 0, 30) && N >= 1; ++t) {
                const long long v = rnd(rng, 1, N);
                cnt[(size_t)v]++, vals.push_back(v);
            }
            for (size_t i = 0; i < vals.size(); ++i)
                for (size_t j = i + 1; j < vals.size(); ++j) e[(size_t)bgcd(vals[i], vals[j])]++;
            DIVT_CHECK(count_pairs_by_gcd(cnt) == e);
        }
    }
    // 演算を渡す：約数の最大値
    std::vector<long long> f = {0, 5, 1, 7, 2, 3, 9};
    divisor_zeta(f, [](long long x, long long y) { return std::max(x, y); });
    DIVT_CHECK(f == std::vector<long long>({0, 5, 5, 7, 5, 5, 9}));
    // ディリクレ積
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int na = (int)rnd(rng, 0, 120), nb = it % 3 ? na : (int)rnd(rng, 0, 120);
        std::vector<long long> a((size_t)na), b((size_t)nb);
        for (auto& x : a) x = rnd(rng, -1000, 1000);
        for (auto& x : b) x = rnd(rng, -1000, 1000);
        const size_t n = std::max(a.size(), b.size());
        std::vector<long long> c(n, 0);
        for (size_t k = 1; k < n; ++k)
            for (size_t d = 1; d <= k; ++d)
                if (k % d == 0 && d < a.size() && k / d < b.size()) c[k] += a[d] * b[k / d];
        DIVT_CHECK(dirichlet_convolution(a, b) == c);
    }
    {  // μ * 1 = ε、φ * 1 = id
        const int N = 1000;
        const std::vector<int> mu = mu_table(N);
        const std::vector<long long> phi = phi_table(N);
        const std::vector<long long> one((size_t)N + 1, 1), mul(mu.begin(), mu.end());
        const std::vector<long long> e = dirichlet_convolution(mul, one), id = dirichlet_convolution(phi, one);
        for (int k = 1; k <= N; ++k) DIVT_CHECK(e[(size_t)k] == (k == 1) && id[(size_t)k] == k);
    }
    // M の約数の上での変換
    std::vector<long long> Ms = {1, 2, 12, 1024, 720720, 735134400LL, 963761198400LL, 1000000007LL, 999999999999999989LL, 897612484786617600LL};
    for (int it = 0; it < 300 * TEST_SCALE; ++it) Ms.push_back(rnd(rng, 1, it % 3 == 0 ? 1000000000000000000LL : (it % 3 == 1 ? 1000000 : 1000)));
    for (const long long M : Ms) {
        const DivisorsOf D(M);
        const std::vector<long long> ds = divisors_of(M);
        const int K = D.size();
        DIVT_CHECK(D.divisors() == ds && K == (int)ds.size() && D.value() == M && D.factors() == factorize(M));
        for (int i = 0; i < K; ++i) DIVT_CHECK(D.index(ds[(size_t)i]) == i && D[i] == ds[(size_t)i]);
        DIVT_CHECK(D.index(M + 1) == -1 && (M < 3 || M % 3 == 0 || D.index(3) == -1));
        std::vector<long long> fv((size_t)K);
        for (auto& x : fv) x = rnd(rng, -1000000, 1000000);
        std::vector<long long> a = fv, b = fv;
        D.zeta(a), D.multiple_zeta(b);
        if (K <= 1500) {
            for (int i = 0; i < K; ++i, ++cases) {
                long long sa = 0, sb = 0;
                for (int j = 0; j < K; ++j) {
                    if (ds[(size_t)i] % ds[(size_t)j] == 0) sa += fv[(size_t)j];
                    if (ds[(size_t)j] % ds[(size_t)i] == 0) sb += fv[(size_t)j];
                }
                DIVT_CHECK(a[(size_t)i] == sa && b[(size_t)i] == sb);
            }
        } else {  // 大きいときは全部 1 の変換 = 約数の個数 / 倍数の個数
            std::vector<long long> o((size_t)K, 1), o2((size_t)K, 1);
            D.zeta(o), D.multiple_zeta(o2);
            for (int i = 0; i < K; i += 37, ++cases) DIVT_CHECK(o[(size_t)i] == count_divisors(ds[(size_t)i]) && o2[(size_t)i] == count_divisors(M / ds[(size_t)i]));
        }
        D.mobius(a), D.multiple_mobius(b);
        DIVT_CHECK(a == fv && b == fv);
    }
    {  // 演算を渡す：倍数の最小値
        const DivisorsOf D(12);  // 1 2 3 4 6 12
        std::vector<long long> g = {5, 3, 9, 8, 1, 7};
        D.multiple_zeta(g, [](long long x, long long y) { return std::min(x, y); });
        DIVT_CHECK(g == std::vector<long long>({1, 1, 1, 7, 1, 7}));
    }
}

static void test_sums(std::mt19937_64& rng, int& cases) {
    const int M = 2000000;
    const std::vector<long long> phi = phi_table(M);
    const std::vector<int> mu = mu_table(M);
    std::vector<Int128> ps((size_t)M + 1, 0);
    std::vector<long long> ms((size_t)M + 1, 0);
    for (int i = 1; i <= M; ++i) ps[(size_t)i] = ps[(size_t)i - 1] + phi[(size_t)i], ms[(size_t)i] = ms[(size_t)i - 1] + mu[(size_t)i];
    for (int it = 0; it < 60 * TEST_SCALE; ++it, ++cases) {
        const long long N = it < 20 ? it : rnd(rng, 1, M);
        DIVT_CHECK(phi_sum(N) == ps[(size_t)N] && mertens(N) == ms[(size_t)N] && count_coprime_pairs(N) == (N ? 2 * ps[(size_t)N] - 1 : 0));
    }
    for (int it = 0; it < 4; ++it, ++cases) {  // 大きい N：1 つ前との差は φ(N)・μ(N)
        const long long N = rnd(rng, 3000000000LL, it % 2 ? 10000000000LL : 3000000000LL + 100);
        DIVT_CHECK(phi_sum(N) - phi_sum(N - 1) == euler_phi(N) && mertens(N) - mertens(N - 1) == mobius(N));
    }
    // gcd の和：小さい場合を全部数える
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const long long K = rnd(rng, 1, 6), len = rnd(rng, 1, 5), m = rnd(rng, 1, 1000000007);
        long long total = 0;
        std::vector<long long> v((size_t)len, 1);
        while (true) {
            long long gg = 0;
            for (long long x : v) gg = bgcd(gg, x);
            total += gg;
            size_t i = 0;
            while (i < v.size() && v[i] == K) v[i] = 1, ++i;
            if (i == v.size()) break;
            v[i]++;
        }
        DIVT_CHECK(sum_gcd_tuples(K, len, m) == total % m);
    }
    // 大きい K：φ の表で Σ φ(d) ⌊K/d⌋^len を直接
    {
        const std::vector<long long> ph = phi_table(3000000);
        for (int it = 0; it < 6; ++it, ++cases) {
            const long long K = it == 0 ? 3000000 : rnd(rng, 1, 300000), len = it == 0 ? 3 : rnd(rng, 1, 1000000000000000000LL), m = rnd(rng, 1, 1000000000000LL);
            long long e = 0;
            for (long long d = 1; d <= K; ++d) e = (long long)(((Int128)e + (Int128)(ph[(size_t)d] % m) * mod_pow(K / d, len, m)) % m);
            DIVT_CHECK(sum_gcd_tuples(K, len, m) == e);
        }
    }
    DIVT_CHECK(sum_gcd_tuples(5, 0, 7) == 0 && sum_gcd_tuples(0, 3, 7) == 0);
    DIVT_CHECK(count_pairs_by_gcd({0, 3000000000LL})[1] == 4499999998500000000LL);
    cases += 2;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_transforms(rng, cases);
    test_sums(rng, cases);
    std::cout << "[divisor_transform] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[divisor_transform] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 10000000;
    std::vector<long long> f((size_t)N + 1);
    for (auto& x : f) x = rnd(rng, 0, 1000);
    auto t0 = std::chrono::steady_clock::now();
    divisor_zeta(f), divisor_mobius(f), multiple_zeta(f), multiple_mobius(f);
    const long long t_tr = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Int128 s1 = phi_sum(100000000000LL);
    const long long t_phi = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const long long s2 = mertens(100000000000LL);
    const long long t_mu = ms(t0);
    std::cout << "[divisor_transform] benchmark: 約数・倍数のゼータ・メビウス 4 回（N = 10^7）" << t_tr << " ms, phi_sum(10^11) " << t_phi << " ms, mertens(10^11) " << t_mu
              << " ms  (check " << (long long)((s1 + s2 + f[12345]) % 7 + 7) % 7 << ")\n";
}

}  // namespace divt_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace divt_demo {

{COMMON_DEMO}

static void demo() {
    banner("約数・倍数の変換（divisor_transform） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  divisor_zeta / divisor_mobius / multiple_zeta / multiple_mobius（f[1..N] を書きかえる。op を渡せる）\n"
              << "  gcd_convolution / lcm_convolution / count_pairs_by_gcd / dirichlet_convolution\n"
              << "  phi_sum / mertens / count_coprime_pairs / sum_gcd_tuples\n"
              << "  DivisorsOf D(M)：M の約数の上で zeta / mobius / multiple_zeta / multiple_mobius（M ≤ 10^18）\n";

    api("【1】 divisor_zeta(f) / divisor_mobius(f) / multiple_zeta(f) / multiple_mobius(f)", "なし（f を書きかえる）",
        "約数の和 g[n] = Σ_{d|n} f[d] / 倍数の和 g[d] = Σ_{d|n} f[n] / それぞれの逆（f[0] は使わない）", "O(N log log N)");
    RUN(std::vector<long long> f = {0, 1, 1, 1, 1, 1, 1}; divisor_zeta(f));
    SHOW(f, "1〜6 の約数の個数");
    RUN(std::vector<long long> g = {0, 1, 1, 1, 1, 1, 1}; multiple_zeta(g));
    SHOW(g, "6 以下の d の倍数の個数 = ⌊6/d⌋");
    note("「gcd がちょうど g」の数は、まず「g で割り切れる」で数えて multiple_mobius（包除）で戻す");

    api("【2】 gcd_convolution(a, b) / lcm_convolution(a, b) / count_pairs_by_gcd(cnt)", "vector",
        "c[k] = Σ_{gcd(i,j)=k} a[i] b[j] / lcm の版 / gcd がちょうど g の組 {i < j} の数", "O(N log log N)");
    SHOW(gcd_convolution<long long>({0, 1, 1, 1, 1}, {0, 1, 1, 1, 1}), "1〜4 から 2 つ（順序あり）選んだ gcd ごとの数");
    SHOW(count_pairs_by_gcd({0, 0, 2, 1, 1}), "{2, 2, 3, 4} の 2 つ組の gcd：gcd 1 が 3 組、gcd 2 が 3 組");
    {
        // ABC206 E：[L, R] の組 (x, y) で gcd = g ≥ 2 かつ x ≠ g かつ y ≠ g の数
        auto solve = [](long long L, long long R) {
            std::vector<long long> F((size_t)R + 1, 0);
            for (long long d = 1; d <= R; ++d) {
                const long long c = R / d - (L - 1) / d;  // [L, R] の d の倍数
                F[(size_t)d] = c * c;
            }
            multiple_mobius(F);  // gcd がちょうど d の組の数
            long long ans = 0;
            for (long long d = 2; d <= R; ++d) {
                ans += F[(size_t)d];
                if (d >= L) ans -= 2 * (R / d - (L - 1) / d) - 1;  // x = d または y = d
            }
            return ans;
        };
        note("ABC206 E（Divide Both）: 「g の倍数どうしの組」c_g^2 を multiple_mobius して gcd ちょうど g の数に");
        std::cout << "    (3, 7) → " << solve(3, 7) << "、(4, 10) → " << solve(4, 10) << "、(1, 10^6) → " << solve(1, 1000000) << "\n";
    }

    api("【3】 phi_sum(N) / mertens(N) / count_coprime_pairs(N) / sum_gcd_tuples(K, len, m)", "Int128 / long long",
        "Σφ（Int128）/ Σμ（N ≤ 10^11 程度）/ 互いに素な組の数 / 列の gcd の和 mod m", "O(N^{2/3}) / O(K log len)");
    SHOW(phi_sum(10), "1 + 1 + 2 + 2 + 4 + 2 + 6 + 4 + 6 + 4");
    SHOW(phi_sum(10000000000LL), "N = 10^10");
    SHOW(mertens(1000000000LL), "Σμ(i)（i ≤ 10^9）");
    SHOW(count_coprime_pairs(3), "(1,1) (1,2) (2,1) (1,3) (3,1) (2,3) (3,2)");
    note("ABC162 E（Sum of gcd of Tuples (Hard)）: sum_gcd_tuples(K, N, 10^9+7)");
    std::cout << "    (N, K) = (3, 2) → " << sum_gcd_tuples(2, 3, 1000000007) << "、(3, 200) → " << sum_gcd_tuples(200, 3, 1000000007) << "、(10^5, 10^5) → "
              << sum_gcd_tuples(100000, 100000, 1000000007) << "\n";

    api("【4】 dirichlet_convolution(a, b)", "vector", "ディリクレ積 c[n] = Σ_{d|n} a[d] b[n/d]（n = 1〜N。乗法的関数の組み合わせ）", "O(N log N)");
    SHOW(dirichlet_convolution<long long>({0, 1, 1, 1, 1, 1, 1}, {0, 1, 1, 1, 1, 1, 1}), "1 * 1 = 約数の個数");
    SHOW(dirichlet_convolution<long long>({0, 1, 2, 3, 4, 5, 6}, {0, 1, -1, -1, 0, -1, 1}), "id * μ = φ");

    api("【5】 DivisorsOf D(M);  D.zeta(f) / D.mobius(f) / D.multiple_zeta(f) / D.multiple_mobius(f)", "なし（f を書きかえる）",
        "M（≤ 10^18）の約数だけの上の約数の和・倍数の和とその逆。f は約数の番号で引く（長さ D.size()、番号は D.index(d)、値は D[i]）",
        "O(D · ω(M))（D は約数の数）");
    RUN(const DivisorsOf D(12));
    SHOW(D.divisors(), "約数（昇順）");
    SHOW(D.index(6), "6 は 4 番目（0 から）");
    {
        // 1〜M の x で gcd(x, M) = d になる数：まず「d の倍数」M/d 個を数えて multiple_mobius
        const long long M = 1000000000000LL;
        const DivisorsOf E(M);
        std::vector<long long> c((size_t)E.size());
        for (int i = 0; i < E.size(); ++i) c[(size_t)i] = M / E[i];
        E.multiple_mobius(c);  // c[i] = gcd(x, M) = E[i] になる x の数（= φ(M / E[i])）
        Int128 s = 0;
        for (int i = 0; i < E.size(); ++i) s += (Int128)E[i] * c[(size_t)i];
        note("Σ_{x=1}^{M} gcd(x, M)：gcd ごとの個数を「倍数の個数 M/d」の multiple_mobius で");
        std::cout << "    M = 10^12（約数 " << E.size() << " 個）→ " << s << "\n";
    }
    {
        // M の約数 a, b で lcm(a, b) = M の組 (a, b) の数：「d の約数どうし」(約数の数)^2 を mobius
        const long long M = 720720;
        const DivisorsOf E(M);
        std::vector<long long> cnt((size_t)E.size(), 1);
        E.zeta(cnt);  // cnt[i] = E[i] の約数の個数
        for (auto& x : cnt) x *= x;
        E.mobius(cnt);  // cnt[i] = lcm がちょうど E[i] の組の数
        note("lcm(a, b) = M となる約数の組 (a, b)：「どちらも d の約数」の数を mobius で「lcm がちょうど d」に");
        std::cout << "    M = 720720 → " << cnt.back() << "\n";
    }
}

}  // namespace divt_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    divt_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    divt_test::benchmark();    // 速度計測
    divt_demo::demo();         // 使い方デモ
    return 0;
}
