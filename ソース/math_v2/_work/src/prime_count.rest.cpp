// =============================================================================
//  ランダムテスト
// =============================================================================
namespace pcnt_test {

{COMMON_TEST}

struct Mod7 {  // テスト用の小さな mod 型（mint の代わり。10^9+7）
    static constexpr unsigned long long MD = 1000000007;
    long long v;
    Mod7(long long x = 0) : v(x % (long long)MD) {
        if (v < 0) v += MD;
    }
    friend Mod7 operator+(Mod7 a, Mod7 b) {
        a.v += b.v;
        if (a.v >= (long long)MD) a.v -= MD;
        return a;
    }
    friend Mod7 operator-(Mod7 a, Mod7 b) {
        a.v -= b.v;
        if (a.v < 0) a.v += MD;
        return a;
    }
    friend Mod7 operator*(Mod7 a, Mod7 b) {
        a.v = (long long)((unsigned long long)a.v * (unsigned long long)b.v % MD);
        return a;
    }
    friend bool operator==(Mod7 a, Mod7 b) { return a.v == b.v; }
};
static long long ipow(long long p, int e) {
    long long r = 1;
    while (e-- > 0) r *= p;
    return r;
}
static void test_multiplicative(std::mt19937_64& rng, int& cases) {
    // 小さい N：最小の素因数の表から各関数を作って累積和と比べる
    const int M = 200000;
    std::vector<int> spf((size_t)M + 1, 0);
    for (int i = 2; i <= M; ++i)
        if (!spf[(size_t)i])
            for (int j = i; j <= M; j += i)
                if (!spf[(size_t)j]) spf[(size_t)j] = i;
    // 0: d, 1: σ1, 2: φ, 3: μ, 4: μ², 5: σ2, 6: i^3
    const int F = 7;
    std::vector<std::vector<Int128>> pre(F, std::vector<Int128>((size_t)M + 1, 0));
    for (int x = 1; x <= M; ++x) {
        Int128 v[F] = {1, 1, 1, 1, 1, 1, 1};
        for (int y = x; y > 1;) {
            const int p = spf[(size_t)y];
            int e = 0;
            while (y % p == 0) y /= p, ++e;
            const Int128 pe = ipow(p, e);
            Int128 s1 = 0, s2 = 0;
            for (int t = 0; t <= e; ++t) s1 += ipow(p, t), s2 += (Int128)ipow(p, t) * ipow(p, t);
            v[0] *= e + 1, v[1] *= s1, v[2] *= pe - pe / p, v[3] *= e == 1 ? -1 : 0, v[4] *= e == 1 ? 1 : 0, v[5] *= s2, v[6] *= pe * pe * pe;
        }
        for (int f = 0; f < F; ++f) pre[(size_t)f][(size_t)x] = pre[(size_t)f][(size_t)x - 1] + v[f];
    }
    auto sig = [](long long p, int e, int k) {  // 1 + p^k + … + p^{ke}
        Int128 s = 0, t = 1;
        for (int i = 0; i <= e; ++i) s += t, t *= ipow(p, k);
        return s;
    };
    auto calc = [&](int f, long long N) -> Int128 {
        switch (f) {
            case 0: return multiplicative_sum<Int128>(N, {2}, [](long long, int e) { return (Int128)e + 1; });
            case 1: return multiplicative_sum<Int128>(N, {1, 1}, [&](long long p, int e) { return sig(p, e, 1); });
            case 2: return multiplicative_sum<Int128>(N, {-1, 1}, [](long long p, int e) { return (Int128)ipow(p, e) - ipow(p, e - 1); });
            case 3: return multiplicative_sum<Int128>(N, {-1}, [](long long, int e) { return (Int128)(e == 1 ? -1 : 0); });
            case 4: return multiplicative_sum<Int128>(N, {1}, [](long long, int e) { return (Int128)(e == 1); });
            case 5: return multiplicative_sum<Int128>(N, {1, 0, 1}, [&](long long p, int e) { return sig(p, e, 2); });
            default: return multiplicative_sum<Int128>(N, {0, 0, 0, 1}, [](long long p, int e) { const Int128 q = ipow(p, e); return q * q * q; });
        }
    };
    for (int it = 0; it < 200 * TEST_SCALE; ++it, ++cases) {
        const long long N = it < 40 ? it : rnd(rng, 0, it % 3 ? 3000 : M);
        const int f = it % F;
        PCNT_CHECK(calc(f, N) == pre[(size_t)f][(size_t)N]);
        if (f == 0) PCNT_CHECK(multiplicative_sum<long long>(N, {2}, [](long long, int e) { return e + 1LL; }) == (long long)pre[0][(size_t)N]);
        if (f == 2) {  // mod 型でも
            const Mod7 m = multiplicative_sum<Mod7>(N, {Mod7(-1), Mod7(1)}, [](long long p, int e) { return Mod7(ipow(p, e) - ipow(p, e - 1)); });
            PCNT_CHECK(m.v == (long long)(pre[2][(size_t)N] % 1000000007));
        }
    }
    // 大きい N：別の方法で
    const long long NB = TEST_SCALE > 1 ? 10000000000LL : 100000000LL;
    {
        Int128 dsum = 0, s1 = 0;  // Σ d(i) = Σ ⌊N/i⌋、Σ σ1(i) = Σ i ⌊N/i⌋
        for (long long i = 1; i <= NB;) {
            const long long q = NB / i, j = NB / q;
            dsum += (Int128)q * (j - i + 1), s1 += (Int128)q * ((Int128)(i + j) * (j - i + 1) / 2);
            i = j + 1;
        }
        PCNT_CHECK(calc(0, NB) == dsum && calc(1, NB) == s1);
        const Int128 n = NB / 100;  // Σ i^3 = (n(n+1)/2)^2
        PCNT_CHECK(calc(6, NB / 100) == (n * (n + 1) / 2) * (n * (n + 1) / 2));
        cases += 3;
    }
    {  // Σφ・Σμ をメモ化再帰で（Φ(N) = N(N+1)/2 − Σ_{d≥2} Φ(N/d)、M(N) = 1 − Σ_{d≥2} M(N/d)）
        const long long N = TEST_SCALE > 1 ? 1000000000LL : 10000000LL;
        std::map<long long, std::pair<Int128, Int128>> memo;
        auto rec = [&](auto&& self, long long v) -> std::pair<Int128, Int128> {
            if (v <= M) return {pre[2][(size_t)v], pre[3][(size_t)v]};
            auto itr = memo.find(v);
            if (itr != memo.end()) return itr->second;
            Int128 a = (Int128)v * (v + 1) / 2, b = 1;
            for (long long d = 2; d <= v;) {
                const long long q = v / d, d2 = v / q;
                const auto s = self(self, q);
                a -= (Int128)(d2 - d + 1) * s.first, b -= (Int128)(d2 - d + 1) * s.second;
                d = d2 + 1;
            }
            return memo[v] = {a, b};
        };
        const auto ex = rec(rec, N);
        PCNT_CHECK(calc(2, N) == ex.first && calc(3, N) == ex.second);
        cases += 2;
    }
    PCNT_CHECK(calc(4, 100) == 61 && calc(3, 1) == 1 && calc(0, 0) == 0);  // 100 以下の平方因子なしは 61 個
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    const int M = 3000000;
    std::vector<char> comp((size_t)M + 1, 0);
    std::vector<long long> pc((size_t)M + 1, 0);
    std::vector<Int128> ps((size_t)M + 1, 0);
    for (int i = 2; i <= M; ++i) {
        if (!comp[(size_t)i])
            for (long long j = (long long)i * i; j <= M; j += i) comp[(size_t)j] = 1;
        pc[(size_t)i] = pc[(size_t)i - 1] + !comp[(size_t)i];
        ps[(size_t)i] = ps[(size_t)i - 1] + (comp[(size_t)i] ? 0 : i);
    }
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const long long N = it < 50 ? it : rnd(rng, 0, it % 3 ? 100000 : M);
        PCNT_CHECK(prime_count(N) == pc[(size_t)N] && prime_sum(N) == ps[(size_t)N]);
        if (it % 10 == 0) {  // 表のすべての値
            const PrimeCountTable T(N);
            for (long long i = 1; i <= N; i = N / (N / i) + 1) PCNT_CHECK(T.count(N / i) == pc[(size_t)(N / i)]);
            for (long long v = 0; v * v <= N; ++v) PCNT_CHECK(T(v) == pc[(size_t)v]);
        }
    }
    // よく知られた値
    PCNT_CHECK(prime_count(10000000000LL) == 455052511 && prime_count(100000000000LL) == 4118054813LL);
    PCNT_CHECK(prime_sum(2000000) == 142913828922LL);  // 200 万未満の素数の和
    // P-smooth：小さい N は全部調べる
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const long long N = rnd(rng, 0, 20000), P = rnd(rng, 1, 60);
        long long e = 0;
        for (long long x = 1; x <= N; ++x) {
            long long y = x;
            for (long long p = 2; p <= P && y > 1; ++p)
                while (y % p == 0) y /= p;
            e += y == 1;
        }
        PCNT_CHECK(count_smooth(N, P) == e);
    }
    PCNT_CHECK(count_smooth(36, 3) == 14 && count_smooth(1000000000000000000LL, 2) == 60 && count_smooth(LLONG_MAX, 1) == 1);
    test_multiplicative(rng, cases);
    std::cout << "[prime_count] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[prime_count] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    const long long c11 = prime_count(100000000000LL);
    const long long t11 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const long long c12 = prime_count(1000000000000LL);
    const long long t12 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Int128 s11 = prime_sum(100000000000LL);
    const long long ts = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const long long sm = count_smooth(10000000000000000LL, 97);
    const long long tsm = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Mod7 m10 = multiplicative_sum<Mod7>(10000000000LL, {Mod7(-1), Mod7(1)}, [](long long p, int e) { return Mod7(ipow(p, e) - ipow(p, e - 1)); });
    const long long tm10 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Mod7 m11 = multiplicative_sum<Mod7>(100000000000LL, {Mod7(-1), Mod7(1)}, [](long long p, int e) { return Mod7(ipow(p, e) - ipow(p, e - 1)); });
    const long long tm11 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Mod7 m12 = multiplicative_sum<Mod7>(1000000000000LL, {Mod7(2)}, [](long long, int e) { return Mod7(e + 1); });
    const long long tm12 = ms(t0);
    std::cout << "[prime_count] benchmark: prime_count(10^11) " << t11 << " ms, prime_count(10^12) " << t12 << " ms, prime_sum(10^11) " << ts
              << " ms, count_smooth(10^16, 97) " << tsm << " ms\n"
              << "[prime_count] benchmark: multiplicative_sum（Σφ mod p）10^10 " << tm10 << " ms, 10^11 " << tm11 << " ms, （Σd mod p）10^12 " << tm12
              << " ms  (check " << (long long)((c11 + c12 + s11 + sm + m10.v + m11.v + m12.v) % 7) << ")\n";
}

}  // namespace pcnt_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace pcnt_demo {

{COMMON_DEMO}

static void demo() {
    banner("素数の個数・和（prime_count） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  prime_count(N) / prime_sum(N)        N 以下の素数の個数 / 和（N ≤ 10^12 程度）\n"
              << "  PrimeCountTable P(N)：P.count(v)       ⌊N/i⌋ の形の v での π(v)\n"
              << "  count_smooth(N, P)                   素因数がすべて P 以下の数の個数\n"
              << "  multiplicative_sum<T>(N, coef, fpe)  乗法的関数の和 Σ_{i≤N} f(i)（Min_25 篩。N ≤ 10^12）\n";

    api("【1】 prime_count(N) / prime_sum(N)", "long long / Int128", "N 以下の素数の個数 / 和（篩では足りない大きな N で。Lucy の方法）", "O(N^{3/4} / log N)");
    SHOW(prime_count(100), "2, 3, 5, …, 97");
    SHOW(prime_count(100000000000LL), "π(10^11)");
    SHOW(prime_sum(100), "100 以下の素数の和");
    SHOW(prime_sum(10000000000LL), "10^10 以下の素数の和（Int128）");

    api("【2】 PrimeCountTable P(N); P.count(v)", "long long", "⌊N/i⌋ の形のすべての v での π(v)（√N 以下の v も）。乗法的関数の和を求める DP の部品", "前処理 O(N^{3/4} / log N)・O(1)");
    RUN(PrimeCountTable P(100));
    SHOW(P.count(100 / 3), "π(33)");
    SHOW(P(7), "π(7)");
    note("例：N 以下で「素数 p < q の積 p·q」になる数 = Σ_{p ≤ √N} (π(N/p) − π(p))（N/p は ⌊N/i⌋ の形なので表にある）");
    {
        auto semiprimes = [](long long N) {
            const PrimeCountTable T(N);
            long long c = 0;
            for (long long p = 2; p * p <= N; ++p)
                if (T(p) != T(p - 1)) c += T(N / p) - T(p);
            return c;
        };
        std::cout << "    N = 100 → " << semiprimes(100) << "（6, 10, 14, 15, …, 95）、N = 10^12 → " << semiprimes(1000000000000LL) << "\n";
    }

    api("【3】 count_smooth(N, P)", "long long", "N 以下で素因数がすべて P 以下の数の個数（1 も数える。素数を 2 組に分けて数え合わせる）", "O(組の大きさ log)");
    SHOW(count_smooth(36, 3), "ABC300 G 例 1：2^a 3^b ≤ 36 は 14 個");
    SHOW(count_smooth(10000000000000000LL, 97), "ABC300 G 例 2");

    api("【4】 multiplicative_sum<T>(N, coef, fpe)", "T（mint・long long・Int128）",
        "Σ_{i=1}^{N} f(i)（f は乗法的）。素数で f(p) = coef[0] + coef[1] p + …（3 次まで）、素数の累乗で f(p^e) = fpe(p, e)。mod なら T = mint",
        "O(N^{3/4} / log N) 程度（10^11 で約 1 秒）");
    SHOW(multiplicative_sum<long long>(10000000000LL, {2}, [](long long, int e) { return e + 1LL; }), "Σ d(i)（約数の個数）、N = 10^10");
    SHOW(multiplicative_sum<long long>(10000000000LL, {-1}, [](long long, int e) { return e == 1 ? -1LL : 0LL; }), "Σ μ(i)（メルテンス関数）");
    SHOW(multiplicative_sum<long long>(10000000000LL, {1}, [](long long, int e) { return e == 1 ? 1LL : 0LL; }), "平方因子を持たない数の個数（μ(i)^2 の和）");
    SHOW(multiplicative_sum<Int128>(10000000000LL, {-1, 1}, [](long long p, int e) { return (Int128)pow_int(p, e) - pow_int(p, e - 1); }), "Σ φ(i)：f(p) = p − 1 → coef {−1, 1}");
    SHOW(multiplicative_sum<Int128>(1000000000LL, {1, 1}, [](long long p, int e) { return (Int128)(pow_int(p, e + 1) - 1) / (p - 1); }), "Σ σ(i)（約数の和）、N = 10^9");
    note("fpe(p, 1) は f(p)（coef の多項式）と同じ値にする。整数の T は溢れたら止まるので、mod の答えなら 00_modint の mint を T に");
}

}  // namespace pcnt_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    pcnt_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pcnt_test::benchmark();    // 速度計測
    pcnt_demo::demo();         // 使い方デモ
    return 0;
}
