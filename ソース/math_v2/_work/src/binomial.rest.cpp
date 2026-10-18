// =============================================================================
//  ランダムテスト
// =============================================================================
namespace binom_test {

{COMMON_TEST}

static const long long LMAX = std::numeric_limits<long long>::max();
static const int EN = 128;  // C(128, 64) < 2^127 なので 128 ビットで厳密に持てる
static std::vector<std::vector<Int128>> EX;
static Int128 ex(long long n, long long r) { return (r < 0 || r > n || n < 0) ? 0 : EX[(size_t)n][(size_t)r]; }
static long long exmod(long long n, long long r, long long m) { return (long long)(ex(n, r) % m); }

template <class mint>
static void test_binomial(std::mt19937_64& rng, int& cases, int iters) {
    const long long p = mint::mod();
    Binomial<mint> B;  // 小さい表から自動で伸びる
    const long long nmax = std::min<long long>(EN, p - 1);
    for (int it = 0; it < iters; ++it, ++cases) {
        const long long n = rnd(rng, -2, nmax), r = rnd(rng, -3, std::max(0LL, n) + 3);
        BINOM_CHECK(B.C(n, r) == mint(exmod(n, r, p)));
        if (n >= 0) {
            mint pr = (r >= 0 && r <= n) ? mint(1) : mint(0);
            for (long long i = 0; i < r && r <= n; ++i) pr *= mint(n - i);
            BINOM_CHECK(B.P(n, r) == pr);
            BINOM_CHECK(B.fact(n) * B.inv_fact(n) == mint(1));
            if (n >= 1) BINOM_CHECK(B.inv(n) * mint(n) == mint(1));
        }
        // nHr
        const long long hn = rnd(rng, 0, 30), hr = rnd(rng, 0, 30);
        if (hn + hr - 1 <= nmax) BINOM_CHECK(B.H(hn, hr) == mint(hr == 0 ? 1 : exmod(hn + hr - 1, hr, p)));
        // 多項係数 = 二項係数の積
        std::vector<int> ks((size_t)rnd(rng, 0, 4));
        long long s = 0;
        mint e = 1;
        for (int& k : ks) k = (int)rnd(rng, 0, 12), s += k, e *= mint(exmod(s, k, p));
        if (s <= nmax) BINOM_CHECK(B.multinomial(ks) == e);
        // カタラン数
        const long long c = rnd(rng, 0, 40);
        if (2 * c <= nmax) BINOM_CHECK(B.catalan(c) == mint((long long)((ex(2 * c, c) / (c + 1)) % p)));
    }
    if (p > 4) BINOM_CHECK(B.multinomial({2, 2}) == mint(6 % p));
    BINOM_CHECK(B.H(0, 0) == mint(1) && B.H(0, 3) == mint(0) && B.catalan(-1) == mint(0));
    // 鏡像法：DP と比べる
    for (int it = 0; it < iters / 50; ++it, ++cases) {
        const int X = (int)rnd(rng, 0, 15), Y = (int)rnd(rng, 0, 15), K = (int)rnd(rng, 0, 4);
        if (X + Y > nmax) continue;
        std::vector<std::vector<long long>> dp(X + 1, std::vector<long long>(Y + 1, 0));
        for (int x = 0; x <= X; ++x)
            for (int y = 0; y <= Y; ++y) {
                if (y > x + K) continue;
                if (x == 0 && y == 0) dp[0][0] = 1;
                if (x > 0) dp[x][y] += dp[x - 1][y];
                if (y > 0) dp[x][y] += dp[x][y - 1];
            }
        BINOM_CHECK(B.paths_not_crossing(X, Y, K) == mint(dp[X][Y] % p) && B.paths(X, Y) == mint(exmod(X + Y, X, p)));
    }
}

template <class mint>
static void test_large(std::mt19937_64& rng, int& cases, int iters) {
    Binomial<mint> B(1000);
    for (int it = 0; it < iters; ++it, ++cases) {  // 大きい n：パスカルの関係式と逆元
        const long long n = rnd(rng, 1, 2000000), r = rnd(rng, 0, n);
        BINOM_CHECK(B.C(n, r) == B.C(n - 1, r - 1) + B.C(n - 1, r));
        BINOM_CHECK(B.fact(n) == B.fact(n - 1) * mint(n) && B.inv(n) * mint(n) == mint(1));
        const long long rr = rnd(rng, 0, 300);
        BINOM_CHECK(binom_large_n<mint>(n, rr) == B.C(n, rr));
    }
    for (int it = 0; it < iters / 10; ++it, ++cases) {  // n が 10^18 程度：C(n, r)·r = C(n, r − 1)·(n − r + 1)
        const long long n = rnd(rng, 1, 1000000000000000000LL), r = rnd(rng, 1, 200);
        BINOM_CHECK(binom_large_n<mint>(n, r) * mint(r) == binom_large_n<mint>(n, r - 1) * mint(n - r + 1));
    }
}

static void test_dynamic(std::mt19937_64& rng, int& cases) {
    Binomial<DynModInt> B(50);
    const int primes[] = {2, 3, 13, 101, 65537, 1000000007, 998244353, 2147483647};
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const int p = primes[it % 8];
        DynModInt::set_mod(p);  // mod を変えても表は作り直される
        const long long n = rnd(rng, 0, std::min<long long>(EN, p - 1)), r = rnd(rng, -1, n + 1);
        BINOM_CHECK(B.C(n, r).val() == exmod(n, r, p));
    }
}

static void test_lucas_and_any(std::mt19937_64& rng, int& cases) {
    const long long ps[] = {2, 3, 5, 7, 11, 13, 97, 101, 9999991};
    for (long long p : ps)
        for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
            const long long n = rnd(rng, -1, EN), r = rnd(rng, -1, std::max(0LL, n) + 1);
            BINOM_CHECK(binom_lucas(n, r, p) == exmod(n, r, p));
        }
    const long long ms[] = {1, 2, 4, 6, 12, 100, 1000000000, 1 << 20, 4782969, 720720, 9999991, 30030, 2 * 9999991LL, 9765625LL * 1024 * 59049};
    for (long long m : ms) {
        for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
            const long long n = rnd(rng, -1, EN), r = rnd(rng, -1, std::max(0LL, n) + 1);
            BINOM_CHECK(binom_mod_any(n, r, m) == exmod(n, r, m));
        }
        if (m == 1) continue;
        for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {  // 大きい n：パスカルの関係式
            const long long n = rnd(rng, 1, 1000000000000000000LL), r = it % 2 ? rnd(rng, 0, n) : rnd(rng, 0, 100);
            BINOM_CHECK(binom_mod_any(n, r, m) == (binom_mod_any(n - 1, r - 1, m) + binom_mod_any(n - 1, r, m)) % m);
        }
    }
    // 平方因子のない m は、素数ごとの Lucas と中国剰余定理で
    const long long sq[] = {2, 3, 5, 7, 11, 13};
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, LMAX), r = rnd(rng, 0, n);
        const long long got = binom_mod_any(n, r, 30030);
        for (long long p : sq) BINOM_CHECK(got % p == binom_lucas(n, r, p));
    }
}

static void test_exact(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, -2, EN), r = rnd(rng, -2, std::max(0LL, n) + 2);
        const Int128 e = ex(n, r);
        long long cap = rnd_ll(rng);
        cap = cap < 0 ? (cap == std::numeric_limits<long long>::min() ? 0 : -cap) : cap;
        if (it % 3 == 0) cap = (long long)std::min<Int128>(e + rnd(rng, -2, 2), LMAX), cap = std::max(cap, 0LL);
        BINOM_CHECK((Int128)binom_cap(n, r, cap) == std::min<Int128>(e, cap));
        if (e <= LMAX) BINOM_CHECK((Int128)binom_exact(n, r) == e);
    }
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {  // 大きい n で r ≤ 2
        const long long n = rnd(rng, 0, it % 2 ? 3037000499LL : LMAX);
        const Int128 c2 = (Int128)n * (n - 1) / 2;
        BINOM_CHECK(binom_exact(n, 0) == 1 && binom_exact(n, n) == 1 && (n == 0 || binom_exact(n, 1) == n));
        BINOM_CHECK((Int128)binom_cap(n, 2, LMAX) == std::min<Int128>(c2, LMAX) && (Int128)binom_cap(n, n - 2, LMAX) == std::min<Int128>(n >= 2 ? c2 : 0, LMAX));
    }
    BINOM_CHECK(binom_exact(66, 33) == 7219428434016265740LL && binom_cap(67, 33, LMAX) == LMAX && binom_cap(1000000000000000000LL, 500000000000000000LL, 5) == 5);
    // パスカルの三角形
    const auto t = pascal_table(66);
    for (int n = 0; n <= 66; ++n)
        for (int r = 0; r <= n; ++r, ++cases) BINOM_CHECK((Int128)t[(size_t)n][(size_t)r] == ex(n, r));
    const auto tm = pascal_table<ModInt998>(EN);
    const auto td = pascal_table<double>(40);
    for (int n = 0; n <= EN; ++n)
        for (int r = 0; r <= n; ++r) {
            BINOM_CHECK(tm[(size_t)n][(size_t)r] == ModInt998(ex(n, r)));
            if (n <= 40) BINOM_CHECK(td[(size_t)n][(size_t)r] == (double)(long long)ex(n, r));
        }
    DynModInt::set_mod(1000000000);  // 合成数の mod でも足し算だけなので使える
    const auto tdm = pascal_table<DynModInt>(EN);
    for (int n = 0; n <= EN; ++n)
        for (int r = 0; r <= n; ++r) BINOM_CHECK(tdm[(size_t)n][(size_t)r].val() == exmod(n, r, 1000000000));
}

template <class mint>
static void test_prefix(std::mt19937_64& rng, int& cases) {
    const Binomial<mint> B(3000);
    BinomialPrefix<mint> P;
    BINOM_CHECK(P.value() == mint(1) && P.n() == 0 && P.m() == 0);
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        switch (rng() % 5) {
            case 0: P.inc_n(); break;
            case 1: if (P.n() > 0) P.dec_n(); break;
            case 2: P.inc_m(); break;
            case 3: P.dec_m(); break;
            default: P.move_to(rnd(rng, 0, 1500), rnd(rng, -3, 1600)); break;
        }
        if (P.n() > 2500) P.move_to(0, P.m());
        if (P.m() > 2600 || P.m() < -5) P.move_to(P.n(), 0);
        mint s = 0;
        for (long long i = 0; i <= P.m(); ++i) s += B.C(P.n(), i);
        BINOM_CHECK(P.value() == s);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    EX.assign(EN + 1, {});
    for (int n = 0; n <= EN; ++n) {
        EX[(size_t)n].assign((size_t)n + 1, 1);
        for (int r = 1; r < n; ++r) EX[(size_t)n][(size_t)r] = EX[(size_t)n - 1][(size_t)r - 1] + EX[(size_t)n - 1][(size_t)r];
    }
    int cases = 0;
    test_binomial<ModInt998>(rng, cases, 100000 * TEST_SCALE);
    test_binomial<ModInt107>(rng, cases, 100000 * TEST_SCALE);
    test_binomial<ModInt<13>>(rng, cases, 20000 * TEST_SCALE);
    test_binomial<ModInt<2>>(rng, cases, 2000 * TEST_SCALE);
    test_large<ModInt998>(rng, cases, 50000 * TEST_SCALE);
    test_large<ModInt107>(rng, cases, 20000 * TEST_SCALE);
    test_dynamic(rng, cases);
    test_lucas_and_any(rng, cases);
    test_exact(rng, cases);
    test_prefix<ModInt998>(rng, cases);
    test_prefix<ModInt107>(rng, cases);
    std::cout << "[binomial] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[binomial] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 10000000;
    auto t0 = std::chrono::steady_clock::now();
    Binomial<ModInt998> B(N);
    const long long t_build = ms(t0);
    std::vector<int> ns(N), rs(N);
    for (int i = 0; i < N; ++i) ns[i] = (int)rnd(rng, 0, N), rs[i] = (int)rnd(rng, 0, ns[i]);
    t0 = std::chrono::steady_clock::now();
    ModInt998 s1 = 0;
    for (int i = 0; i < N; ++i) s1 += B.C(ns[i], rs[i]);
    const long long t_c = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const ModInt998 s2 = binom_large_n<ModInt998>(1000000000000000000LL, 10000000);
    const long long t_large = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s3 = 0;
    for (int i = 0; i < 1000000; ++i) s3 += binom_lucas((long long)rng() >> 2, (long long)rng() >> 3, 1000003);
    const long long t_lucas = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s4 = 0;
    for (int i = 0; i < 100000; ++i) s4 += binom_mod_any((long long)rng() >> 2, (long long)rng() >> 3, 1000000000);
    const long long t_any = ms(t0);
    t0 = std::chrono::steady_clock::now();
    DynModInt::set_mod(1000000000);
    const auto t = pascal_table<DynModInt>(5000);
    const long long t_pascal = ms(t0);
    std::cout << "[binomial] benchmark: Binomial(10^7) の表 " << t_build << " ms, C × 10^7 " << t_c << " ms, binom_large_n(10^18, 10^7) " << t_large
              << " ms\n[binomial] benchmark: binom_lucas（p = 10^6+3）× 10^6 " << t_lucas << " ms, binom_mod_any（m = 10^9）× 10^5 " << t_any
              << " ms, pascal_table<DynModInt>(5000) " << t_pascal << " ms  (check " << (s1 + s2).val() % 7 + (s3 + s4 + t[5000][2500].val()) % 7 << ")\n";
}

}  // namespace binom_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace binom_demo {

{COMMON_DEMO}

static void demo() {
    banner("二項係数（binomial） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Binomial<mint> B(N)：B.C / B.P / B.H / B.fact / B.inv_fact / B.inv / B.multinomial / B.catalan / B.paths / B.paths_not_crossing\n"
              << "  BinomialPrefix<mint> P           Σ_{i≤m} C(n, i) を n・m を ±1 して O(1)（Mo のアルゴリズムに）\n"
              << "  binom_large_n<mint>(n, r)        n が大きく r が小さい nCr\n"
              << "  binom_lucas(n, r, p)             小さい素数 p での nCr（Lucas）\n"
              << "  binom_mod_any(n, r, m)           合成数の mod での nCr\n"
              << "  binom_exact / binom_cap          厳密な nCr / 上限で止めた nCr\n"
              << "  pascal_table<T>(N)               パスカルの三角形（T は何でも）\n";

    api("【1】 Binomial<mint> B(N)：B.C(n, r) / B.P(n, r) / B.H(n, r) / B.fact(n) / B.inv(n)", "mint",
        "階乗の表を作って nCr などを O(1) で。範囲外（r < 0、r > n）は 0。N を超える n が来たら表を伸ばす（mod は素数）", "前処理 O(N)・O(1)");
    note("using mint = ModInt107;  // 10^9+7（998244353 なら ModInt998）");
    using mint = ModInt107;
    RUN(Binomial<mint> B(200000));
    SHOW(B.C(5, 2), "5C2");
    SHOW(B.C(5, 7), "範囲外は 0");
    SHOW(B.P(5, 2), "5P2 = 20");
    SHOW(B.H(3, 2), "3 種類から重複を許して 2 個 = 4C2 = 6");
    SHOW(B.fact(10), "10!");
    SHOW(B.inv(2), "1/2 mod 10^9+7");
    SHOW(B.C(1000000, 500000), "表は自動で 10^6 まで伸びる");
    {
        // ABC132 D：K 個の青と N−K 個の赤。青をちょうど i 回で回収 = C(N−K+1, i)·C(K−1, i−1)
        note("ABC132 D（Blue and Red Balls）: C(N − K + 1, i) · C(K − 1, i − 1)");
        for (auto [N, K] : std::vector<std::pair<int, int>>{{5, 3}, {2000, 3}}) {
            std::cout << "    N = " << N << ", K = " << K << " →";
            for (int i = 1; i <= K; ++i) std::cout << " " << B.C(N - K + 1, i) * B.C(K - 1, i - 1);
            std::cout << "\n";
        }
        // ABC145 D：(+1, +2) を a 回、(+2, +1) を b 回。a + 2b = X、2a + b = Y
        auto knight = [&](long long X, long long Y) {
            if ((X + Y) % 3 != 0) return mint(0);
            const long long a = (2 * Y - X) / 3, b = (2 * X - Y) / 3;
            return B.C(a + b, a);  // a, b が負なら 0
        };
        note("ABC145 D（Knight）: 動き方の回数 a, b を連立方程式で出して C(a + b, a)");
        std::cout << "    (3, 3) → " << knight(3, 3) << "、(2, 2) → " << knight(2, 2) << "、(999999, 999999) → " << knight(999999, 999999) << "\n";
        // ABC110 D：M = Π p^e。各素数の e 個を N 個の箱に分ける = H(N, e)
        auto factorized_count = [&](long long N, long long M) {
            mint r = 1;
            for (long long p = 2; p * p <= M; ++p) {
                int e = 0;
                while (M % p == 0) M /= p, ++e;
                if (e) r *= B.H(N, e);
            }
            if (M > 1) r *= B.H(N, 1);
            return r;
        };
        note("ABC110 D（Factorization）: 素因数 p^e ごとに H(N, e) を掛ける（素因数分解は 02_factorize でも）");
        std::cout << "    (2, 6) → " << factorized_count(2, 6) << "、(3, 12) → " << factorized_count(3, 12) << "、(100000, 10^9) → " << factorized_count(100000, 1000000000) << "\n";
    }

    api("【2】 B.multinomial(ks) / B.catalan(n) / B.paths(x, y) / B.paths_not_crossing(x, y, k)", "mint",
        "多項係数 (Σk)!/Πk! / カタラン数 / 右・上だけの道の数 / ずっと y ≤ x + k の道の数（鏡像法）", "O(|ks|) / O(1)");
    SHOW(B.multinomial({2, 1, 1}), "AABC の並べ方 = 4! / (2! 1! 1!) = 12");
    SHOW(B.catalan(3), "正しい括弧列（3 組）");
    SHOW(B.paths(2, 3), "5C2 = 10");
    SHOW(B.paths_not_crossing(3, 3, 0), "対角線を越えない（= カタラン数）");
    note("ABC205 E（White and Black Balls）: 白 N・黒 M、どの途中でも 白 ≤ 黒 + K → paths_not_crossing(M, N, K)");
    std::cout << "    (N, M, K) = (2, 3, 1) → " << B.paths_not_crossing(3, 2, 1) << "、(1, 0, 0) → " << B.paths_not_crossing(0, 1, 0)
              << "、(10^6, 10^6, 10^6) → " << B.paths_not_crossing(1000000, 1000000, 1000000) << "\n";
    {
        // ABC167 E：隣が同じ色の組が k 個 → M (M−1)^{N−1−k} C(N−1, k) を k = 0..K で足す（mod 998244353）
        Binomial<ModInt998> B2(1);
        auto solve = [&](long long N, long long M, long long K) {
            ModInt998 s = 0;
            for (long long k = 0; k <= K; ++k) s += ModInt998(M) * ModInt998(M - 1).pow(N - 1 - k) * B2.C(N - 1, k);
            return s;
        };
        note("ABC167 E（Colorful Blocks）: Σ_{k ≤ K} M (M − 1)^{N−1−k} C(N − 1, k)（Binomial<ModInt998>）");
        std::cout << "    (3, 2, 1) → " << solve(3, 2, 1) << "、(100, 100, 0) → " << solve(100, 100, 0) << "、(60522, 114575, 7559) → " << solve(60522, 114575, 7559) << "\n";
    }

    api("【3】 binom_large_n<mint>(n, r)", "mint", "n ≤ 10^18 で r が小さい nCr（r < mod）。分子を r 個掛けて r! で割る", "O(r)");
    SHOW(binom_large_n<mint>(1000000000000000000LL, 3), "10^18 C 3");
    {
        auto solve = [](long long n, long long a, long long b) { return mint(2).pow(n) - 1 - binom_large_n<mint>(n, a) - binom_large_n<mint>(n, b); };
        note("ABC156 D（Bouquet）: 2^n − 1 − C(n, a) − C(n, b)（n ≤ 10^9、a, b ≤ 2×10^5）");
        std::cout << "    (4, 1, 3) → " << solve(4, 1, 3) << "、(10^9, 141421, 173205) → " << solve(1000000000, 141421, 173205) << "\n";
    }

    api("【4】 binom_lucas(n, r, p) / binom_mod_any(n, r, m)", "long long",
        "小さい素数 p での nCr（Lucas）/ 合成数の mod での nCr（素数べきごとに計算して中国剰余定理）。n, r ≤ 10^18", "O(log n)・前処理 O(p) / O(Σ p^e)");
    SHOW(binom_lucas(100, 50, 7), "100 = (202)_7、50 = (101)_7 → C(2,1)·C(0,0)·C(2,1) = 4");
    SHOW(binom_lucas(1000000000000000000LL, 12345, 1000003), "n が 10^18 でも（p = 10^6+3）");
    SHOW(binom_mod_any(10, 5, 12), "252 mod 12 = 0");
    SHOW(binom_mod_any(1000000000000000000LL, 500000000000000000LL, 1000000000), "mod 10^9（2^9·5^9）");
    {
        // ARC117 C：B, W, R を 0, 1, 2。いちばん上 = (−1)^{N−1} Σ C(N−1, i) a_i (mod 3)
        auto top = [](const std::string& s) {
            const int n = (int)s.size();
            long long t = 0;
            for (int i = 0; i < n; ++i) t += binom_lucas(n - 1, i, 3) * (s[(size_t)i] == 'B' ? 0 : s[(size_t)i] == 'W' ? 1 : 2);
            if (n % 2 == 0) t = -t;
            return "BWR"[mod_floor(t, 3)];
        };
        note("ARC117 C（Tricolor Pyramid）: 色を 0, 1, 2 にして (−1)^{N−1} Σ C(N−1, i) a_i mod 3（binom_lucas）");
        std::cout << "    BWR → " << top("BWR") << "、RRBB → " << top("RRBB") << "、BWWRBW → " << top("BWWRBW") << "、WWBRBBWB → " << top("WWBRBBWB")
                  << "、BWBRRBBRWBRBBBRRBWWWR → " << top("BWBRRBBRWBRBBBRRBWWWR") << "\n";
    }

    api("【5】 binom_exact(n, r) / binom_cap(n, r, cap)", "long long", "厳密な nCr（収まらなければ止める）/ min(nCr, cap)（溢れない）", "O(min(r, 64))");
    SHOW(binom_exact(60, 30), "厳密な値");
    SHOW(binom_exact(1000000000, 2), "n が大きくても r が小さければ");
    SHOW(binom_cap(100, 50, 1000000000000000000LL), "10^18 を超えるので 10^18");

    api("【6】 pascal_table<T>(N)", "vector<vector<T>>",
        "t[n][r] = nCr（0 ≤ r ≤ n ≤ N）。足し算だけなので合成数の mod（DynModInt）でも使える。long long は N ≤ 66", "O(N^2)");
    RUN(auto t = pascal_table(60));
    SHOW(t[60][30], "long long で厳密に");
    {
        // ABC202 D：a A 個・b B 個の文字列で辞書順 K 番目（K は 1 から）
        auto kth = [&](int A, int Bn, long long K) {
            std::string s;
            while (A > 0 || Bn > 0) {
                const long long c = A > 0 ? t[(size_t)(A + Bn - 1)][(size_t)(A - 1)] : 0;  // 先頭を a にしたときの数
                if (A > 0 && K <= c) s += 'a', --A;
                else s += 'b', K -= c, --Bn;
            }
            return s;
        };
        note("ABC202 D（aab aba baa）: 先頭を a にしたときの数 C(A+B−1, A−1) と K を比べて 1 文字ずつ決める");
        std::cout << "    (2, 2, 4) → " << kth(2, 2, 4) << "、(30, 30, 118264581564861424) → " << kth(30, 30, 118264581564861424LL) << "\n";
    }
    {
        // ABC425 E：Σ C_i 個を並べる数 = Π C(累積, C_i)（mod M は合成数のこともある）
        auto solve = [](long long M, const std::vector<std::vector<int>>& cases) {
            DynModInt::set_mod((int)M);
            const auto pt = pascal_table<DynModInt>(5000);
            std::string out;
            for (const auto& c : cases) {
                DynModInt r = 1;
                int s = 0;
                for (int x : c) s += x, r *= pt[(size_t)s][(size_t)x];
                out += (out.empty() ? "" : ", ") + std::to_string(r.val());
            }
            return out;
        };
        note("ABC425 E（Count Sequences 2）: mod M は合成数のこともあるので pascal_table<DynModInt>(5000) で Π C(累積, C_i)");
        std::cout << "    M = 10^9: " << solve(1000000000, {{2, 2}, {1, 1, 1, 1, 1}, {1, 2, 3, 4, 5, 6}}) << "\n"
                  << "    M = 998244353: " << solve(998244353, {{1}, {4, 2, 5}, std::vector<int>(10, 500)}) << "\n";
    }

    api("【7】 BinomialPrefix<mint> P(最大の n);  P.move_to(n, m) / P.inc_n() / P.dec_n() / P.inc_m() / P.dec_m() / P.value()", "mint",
        "S(n, m) = Σ_{i=0}^{m} C(n, i)（n 個から m 個以下を選ぶ数）。n・m を 1 動かすのが O(1)。クエリを Mo の順に並べて全体 O(N √Q)。mod は 3 以上の素数",
        "1 回の移動 O(1)");
    RUN(BinomialPrefix<ModInt998> P(100));
    SHOW(P.move_to(10, 3), "C(10,0) + C(10,1) + C(10,2) + C(10,3) = 1 + 10 + 45 + 120");
    SHOW(P.move_to(10, 10), "m ≥ n なら 2^n");
    {
        // 多くのクエリ (n, m)：n を √ ごとのブロックに分け、ブロック内は m の順（Mo のアルゴリズム）
        const std::vector<std::pair<int, int>> qs = {{5, 2}, {100, 50}, {7, 7}, {100000, 49999}, {3, 0}, {99999, 50000}};
        const int Bsz = 300;
        std::vector<int> ord(qs.size());
        std::iota(ord.begin(), ord.end(), 0);
        std::sort(ord.begin(), ord.end(), [&](int a, int b) {
            const int ba = qs[(size_t)a].first / Bsz, bb = qs[(size_t)b].first / Bsz;
            return ba != bb ? ba < bb : (ba % 2 ? qs[(size_t)a].second > qs[(size_t)b].second : qs[(size_t)a].second < qs[(size_t)b].second);
        });
        BinomialPrefix<ModInt998> Q(100000);
        std::vector<ModInt998> ans(qs.size());
        for (int i : ord) ans[(size_t)i] = Q.move_to(qs[(size_t)i].first, qs[(size_t)i].second);
        note("Mo のアルゴリズム：n のブロック順・ブロック内は m の順に move_to すると、移動の合計が O(N √Q)");
        std::cout << "    ";
        for (size_t i = 0; i < qs.size(); ++i) std::cout << "(" << qs[i].first << ", " << qs[i].second << ") → " << ans[i] << (i + 1 < qs.size() ? "、" : "\n");
    }
}

}  // namespace binom_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    binom_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    binom_test::benchmark();    // 速度計測
    binom_demo::demo();         // 使い方デモ
    return 0;
}
