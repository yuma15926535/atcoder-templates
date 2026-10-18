// =============================================================================
//  ランダムテスト
// =============================================================================
namespace dlog_test {

{COMMON_TEST}

static long long bgcd(long long a, long long b) { return b ? bgcd(b, a % b) : a; }
static long long brute_order(long long a, long long m) {  // なければ −1
    if (m == 1) return 1;
    if (bgcd(a % m, m) != 1) return -1;
    long long x = a % m, k = 1;
    while (x != 1 % m) x = x * (a % m) % m, ++k;
    return k;
}
static long long rand_prime(std::mt19937_64& rng, long long lo, long long hi) {
    while (true) {
        const long long p = rnd(rng, lo, hi);
        if (is_prime(p)) return p;
    }
}

static void test_order(std::mt19937_64& rng, int& cases) {
    for (long long n = 1; n <= 300; ++n, ++cases) {  // λ(n)：互いに素な a の位数の lcm
        long long l = 1;
        for (long long a = 1; a <= n; ++a) {
            const long long o = brute_order(a, n);
            if (o > 0) l = l / bgcd(l, o) * o;
        }
        DLOG_CHECK(carmichael_lambda(n) == l);
        const MultOrder ord(n);
        for (int t = 0; t < 10; ++t) {
            const long long a = rnd(rng, -3 * n, 3 * n), o = brute_order(mod_floor(a, n), n);
            const Maybe<long long> r = ord(a);
            DLOG_CHECK(r.ok == (o > 0) && (!r.ok || r.value == o));
        }
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい m：a^k ≡ 1 で、k の約数では 1 にならない
        const long long m = it % 2 ? rnd(rng, 2, 1000000000000000000LL) : rand_prime(rng, 2, 1000000000000000000LL);
        const long long a = rnd(rng, 0, m - 1);
        const Maybe<long long> r = mult_order(a, m);
        DLOG_CHECK(r.ok == (bgcd(a, m) == 1));
        if (!r.ok) continue;
        DLOG_CHECK(mod_pow(a, r.value, m) == 1 && carmichael_lambda(m) % r.value == 0);
        for (long long q : prime_factors(r.value)) DLOG_CHECK(mod_pow(a, r.value / q, m) != 1);
    }
    // 原始根
    for (long long p = 2; p < 3000; ++p) {
        if (!is_prime(p)) continue;
        ++cases;
        long long g = 1;
        while (brute_order(g, p) != p - 1) ++g;
        DLOG_CHECK(primitive_root(p) == g);
    }
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const long long p = rand_prime(rng, 3, 1000000000000000000LL), g = primitive_root(p);
        DLOG_CHECK(mult_order(g, p).value == p - 1);
    }
    DLOG_CHECK(primitive_root(998244353) == 3 && primitive_root(1000000007) == 5);
}

static void test_sqrt(std::mt19937_64& rng, int& cases) {
    for (long long p = 2; p < 600; ++p) {
        if (!is_prime(p)) continue;
        for (long long a = 0; a < p; ++a, ++cases) {
            long long e = -1;
            for (long long x = 0; x < p && e < 0; ++x)
                if (x * x % p == a) e = x;  // 小さい方
            const Maybe<long long> r = mod_sqrt(a, p);
            DLOG_CHECK(r.ok == (e >= 0) && (!r.ok || r.value == e));
        }
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long p = rand_prime(rng, 3, it % 2 ? 1000000000LL : 1000000000000000000LL), a = rnd(rng, 0, p - 1);
        const Maybe<long long> r = mod_sqrt(a, p);
        DLOG_CHECK(r.ok == (a == 0 || mod_pow(a, (p - 1) / 2, p) == 1));
        if (r.ok) DLOG_CHECK(mod_mul(r.value, r.value, p) == a && r.value <= p - r.value);
    }
}

static void test_dlog(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 6000 * TEST_SCALE; ++it, ++cases) {
        const long long m = rnd(rng, 1, it % 3 ? 60 : 600), a = rnd(rng, -5, 3 * m), b = rnd(rng, -5, 3 * m);
        long long e = -1, cur = 1 % m;
        for (long long x = 0; x <= 2 * m + 70 && e < 0; ++x) {  // 周期は m 以下（前の部分は 64 以下）
            if (cur == mod_floor(b, m)) e = x;
            cur = (long long)((Int128)cur * mod_floor(a, m) % m);
        }
        const Maybe<long long> r = discrete_log(a, b, m);
        DLOG_CHECK(r.ok == (e >= 0) && (!r.ok || r.value == e));
    }
    for (int it = 0; it < 30 * TEST_SCALE; ++it, ++cases) {  // 大きい m：答えから作る（√m の表を作るので回数は少なめ）
        const long long m = it % 2 ? rnd(rng, 2, 10000000000LL) : rand_prime(rng, 2, it >= 4 ? 10000000000LL : 1000000000000LL);
        const long long a = rnd(rng, 0, m - 1), x = rnd(rng, 0, m), b = mod_pow(a, x, m);
        const Maybe<long long> r = discrete_log(a, b, m);
        DLOG_CHECK(r.ok && r.value <= x && mod_pow(a, r.value, m) == b);
        if (r.value > 0) DLOG_CHECK(r.value < 64 || mod_pow(a, r.value - 1, m) != b);
    }
}

static void test_tower(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long m = rnd(rng, 1, it % 2 ? 100 : 1000000000000000000LL);
        // 3 段：b^c が long long に入るとき mod_pow(a, b^c, m) と
        const long long a = rnd(rng, 0, 30), b = rnd(rng, 0, 30), c = rnd(rng, 0, 12);
        const long long e = pow_cap(b, c, LLONG_MAX);
        if (e < LLONG_MAX) DLOG_CHECK(pow_tower_mod({a, b, c}, m) == mod_pow(a, e, m));
        // 2 段・1 段
        DLOG_CHECK(pow_tower_mod({a, b}, m) == mod_pow(a, b, m) && pow_tower_mod({a}, m) == a % m);
    }
    DLOG_CHECK(pow_tower_mod({2, 2, 2, 2}, 1000000007) == 65536 && pow_tower_mod({2, 2, 2, 2, 2}, 1000000007) == mod_pow(2, 65536, 1000000007));
    DLOG_CHECK(pow_tower_mod({2, 3, 4}, 1000) == mod_pow(2, 81, 1000) && pow_tower_mod({0, 0}, 7) == 1 && pow_tower_mod({0, 5, 0}, 7) == 0 && pow_tower_mod({}, 7) == 1);
    // 指数が巨大：周期で確かめる（a^E と a^{E + λ}、E ≥ 64）
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long m = rnd(rng, 1, 100000), a = rnd(rng, 0, 1000), b = rnd(rng, 2, 1000), c = rnd(rng, 30, 1000);
        // E = b^c は巨大。E mod φ(m) を直接計算して比べる（E ≥ 64）
        const long long ph = euler_phi(m), r = mod_pow(b, c, ph);
        DLOG_CHECK(pow_tower_mod({a, b, c}, m) == mod_pow(a, r + ph, m));
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_order(rng, cases);
    test_sqrt(rng, cases);
    test_dlog(rng, cases);
    test_tower(rng, cases);
    std::cout << "[discrete_log] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[discrete_log] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const long long P = 9999999999971LL;  // 素数
    auto t0 = std::chrono::steady_clock::now();
    const MultOrder ord(P);
    long long s1 = 0;
    for (int i = 0; i < 200000; ++i) s1 += ord(rnd(rng, 1, P - 1)).value % 7;
    const long long t_ord = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for (int i = 0; i < 100000; ++i) s2 += mod_sqrt(rnd(rng, 0, P - 1), P).value % 7;
    const long long t_sqrt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Maybe<long long> d = discrete_log(3, 123456789, 1000000000039LL);
    const long long t_dlog = ms(t0);
    std::cout << "[discrete_log] benchmark: MultOrder（m ≈ 10^13）× 2×10^5 " << t_ord << " ms, mod_sqrt（p ≈ 10^13）× 10^5 " << t_sqrt << " ms, discrete_log（m ≈ 10^12）"
              << t_dlog << " ms  (check " << (s1 + s2 + d.value) % 7 << ")\n";
}

}  // namespace dlog_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace dlog_demo {

{COMMON_DEMO}

static void demo() {
    banner("位数・原始根・離散対数（discrete_log） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  carmichael_lambda(n) / mult_order(a, m) / MultOrder ord(m)\n"
              << "  primitive_root(p) / mod_sqrt(a, p) / discrete_log(a, b, m) / pow_tower_mod(as, m)\n";

    api("【1】 mult_order(a, m) / MultOrder ord(m); ord(a) / carmichael_lambda(n)", "Maybe<long long> / long long",
        "a^k ≡ 1 (mod m) となる最小の k（gcd(a, m) ≠ 1 なら ok = false）/ λ(n)", "O(log^2)（m の素因数分解）");
    SHOW(mult_order(10, 7).value, "10^6 ≡ 1 (mod 7)");
    SHOW(mult_order(2, 6).ok, "gcd(2, 6) ≠ 1 なので位数はない");
    SHOW(carmichael_lambda(1000000000), "λ(10^9) = lcm(λ(2^9), λ(5^9))");
    {
        // ABC222 G：222…2（k 桁）が K の倍数になる最小の k。2(10^k − 1)/9 ≡ 0 (mod K) → 10^k ≡ 1 (mod 9K / gcd(K, 2))
        auto solve = [](long long K) {
            const long long M = 9 * (K % 2 == 0 ? K / 2 : K);
            const Maybe<long long> o = mult_order(10, M);
            return o.ok ? o.value : -1;
        };
        note("ABC222 G（222）: 10^k ≡ 1 (mod 9K / gcd(K, 2)) の最小の k = mult_order(10, M)");
        std::cout << "    K = 1 → " << solve(1) << "、7 → " << solve(7) << "、10 → " << solve(10) << "、999983 → " << solve(999983) << "\n";
    }
    {
        // ABC335 G：A_i^k ≡ A_j となる k がある ⇔ ord(A_j) | ord(A_i)（P は素数）
        auto solve = [](long long P, const std::vector<long long>& A) {
            const MultOrder ord(P);
            std::vector<long long> o;
            for (long long a : A) o.push_back(ord(a).value);
            long long c = 0;
            for (long long x : o)
                for (long long y : o) c += x % y == 0;
            return c;
        };
        note("ABC335 G（Discrete Logarithm Problems）: MultOrder で各 A_i の位数を出し、位数が割り切る組を数える（例は少ないので全部の組で）");
        std::cout << "    例 1 → " << solve(13, {2, 3, 5}) << "、例 2 → " << solve(2, {1, 1, 1, 1, 1}) << "、例 3 → "
                  << solve(9999999999971LL, {141592653589LL, 793238462643LL, 383279502884LL, 197169399375LL, 105820974944LL, 592307816406LL, 286208998628LL, 34825342117LL, 67982148086LL, 513282306647LL})
                  << "\n";
    }

    api("【2】 primitive_root(p) / mod_sqrt(a, p)", "long long / Maybe<long long>", "素数 p の原始根（最小）/ x^2 ≡ a (mod p) の小さい方の x（なければ ok = false）", "O(log^2 p)");
    SHOW(primitive_root(998244353), "998244353 の原始根");
    SHOW(mod_sqrt(2, 7).value, "3^2 = 9 ≡ 2（もう 1 つは 4）");
    SHOW(mod_sqrt(3, 7).ok, "3 は mod 7 の平方数ではない");

    api("【3】 discrete_log(a, b, m)", "Maybe<long long>", "a^x ≡ b (mod m) となる最小の x ≥ 0（m は合成数でもよい。BSGS）", "O(√m)");
    SHOW(discrete_log(2, 3, 5).value, "2^3 = 8 ≡ 3");
    SHOW(discrete_log(2, 0, 8).value, "2^3 ≡ 0 (mod 8)：互いに素でなくても");
    SHOW(discrete_log(4, 2, 7).value, "4^0 = 1、4^1 = 4、4^2 = 16 ≡ 2");
    SHOW(discrete_log(3, 0, 7).ok, "3^x ≡ 0 にはならない");

    api("【4】 pow_tower_mod(as, m)", "long long", "a_0^(a_1^(a_2^(…))) mod m（指数を φ(m) で割った余り + φ(m) にしていく。a は 0 以上）", "O(段数 × m の素因数分解)");
    SHOW(pow_tower_mod({2, 3, 4}, 1000), "2^(3^4) = 2^81 mod 1000");
    SHOW(pow_tower_mod({3, 3, 3, 3}, 1000000007), "3^(3^(3^3))");
    note("ABC228 E（Integer Sequence Fair）: M^(K^N) mod 998244353 = pow_tower_mod({M, K, N}, 998244353)（M が p の倍数のときも正しい）");
    std::cout << "    (N, K, M) = (2, 2, 2) → " << pow_tower_mod({2, 2, 2}, 998244353) << "、(3, 14, 15926535) → " << pow_tower_mod({15926535, 14, 3}, 998244353) << "\n";
}

}  // namespace dlog_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    dlog_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    dlog_test::benchmark();    // 速度計測
    dlog_demo::demo();         // 使い方デモ
    return 0;
}
