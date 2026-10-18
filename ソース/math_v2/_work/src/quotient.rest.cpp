// =============================================================================
//  ランダムテスト
// =============================================================================
namespace quot_test {

{COMMON_TEST}

static const long long LMAX = std::numeric_limits<long long>::max();

static void check_blocks(long long N, int& cases) {  // 区間がすき間なく並び、各区間で値が一定で、区間が最大
    ++cases;
    const std::vector<QuotientBlock> f = quotient_blocks(N), c = quotient_blocks_ceil(N);
    long long x = 1;
    for (const auto& b : f) {
        QUOT_CHECK(b.l == x && b.l < b.r && N / b.l == b.q && N / (b.r - 1) == b.q && (b.r > N || N / b.r != b.q));
        x = b.r;
    }
    QUOT_CHECK(x == N + 1 || (N == 0 && f.empty()));
    QUOT_CHECK((long long)f.size() <= 2 * isqrt(N) + 1);
    x = 1;
    for (const auto& b : c) {
        auto cd = [&](long long i) { return (N + i - 1) / i; };
        QUOT_CHECK(b.l == x && b.l < b.r && cd(b.l) == b.q && cd(b.r - 1) == b.q && (b.r > N || cd(b.r) != b.q));
        x = b.r;
    }
    QUOT_CHECK(x == N + 1 || (N == 0 && c.empty()));
}

static void test_blocks(std::mt19937_64& rng, int& cases) {
    for (long long N = 0; N <= 2000; ++N) {
        check_blocks(N, cases);
        // 1 つずつ数える
        Int128 s1 = 0, s2 = 0, s3 = 0;
        for (long long i = 1; i <= N; ++i) s1 += N / i, s2 += (Int128)i * (N / i), s3 += N % i;
        QUOT_CHECK(sum_floor_div(N) == s1 && sum_i_floor_div(N) == s2 && sum_mod_i(N) == s3);
        long long cnt = 0;
        for_each_quotient(N, [&](long long q, long long l, long long r) {
            for (long long i = l; i < r; ++i) QUOT_CHECK(N / i == q);
            cnt += r - l;
        });
        QUOT_CHECK(cnt == N);
    }
    for (int it = 0; it < 100 * TEST_SCALE; ++it) check_blocks(rnd(rng, 0, 100000000LL), cases);
    check_blocks(1000000000000LL, cases), check_blocks(999999999999LL, cases);
    {  // N = 10^14（区間は 2×10^7 個）：for_each_quotient で区間がつながっていることと q の確認（一部）
        const long long N = 100000000000000LL;
        long long x = 1, k = 0;
        for_each_quotient(N, [&](long long q, long long l, long long r) {
            if (++k % 1000003 == 0 || l < 1000) QUOT_CHECK(l == x && N / l == q && N / (r - 1) == q && (r > N || N / r != q));
            x = r;
        });
        QUOT_CHECK(x == N + 1 && k <= 2 * isqrt(N));
    }
    // 2 つの商
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long N = rnd(rng, 0, it % 2 ? 300 : 1000000), M = rnd(rng, 0, it % 3 ? 300 : 1000000);
        long long x = 1;
        for (const auto& b : quotient_blocks2(N, M)) {
            QUOT_CHECK(b.l == x && b.l < b.r && N / b.l == b.q1 && M / b.l == b.q2 && N / (b.r - 1) == b.q1 && M / (b.r - 1) == b.q2);
            QUOT_CHECK(b.r > std::min(N, M) || N / b.r != b.q1 || M / b.r != b.q2);
            x = b.r;
        }
        QUOT_CHECK(x == std::min(N, M) + 1 || (std::min(N, M) == 0 && x == 1));
    }
    // 大きい N：1 つ前との差は約数の個数・約数の和
    for (int it = 0; it < 40 * TEST_SCALE; ++it, ++cases) {
        const long long N = rnd(rng, 2, it % 4 ? 1000000000LL : 100000000000LL);
        long long d = 0;
        Int128 sig = 0;
        for (long long i = 1; i * i <= N; ++i)
            if (N % i == 0) d += (i * i == N ? 1 : 2), sig += i + (i * i == N ? 0 : N / i);
        QUOT_CHECK(sum_floor_div(N) - sum_floor_div(N - 1) == d && sum_i_floor_div(N) - sum_i_floor_div(N - 1) == sig);
    }
}

static void test_floor_sum(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 40);
        const long long m = it % 3 == 0 ? rnd(rng, 1, 10) : rnd(rng, 1, it % 2 ? 1000000000LL : LMAX);
        const long long a = it % 4 == 0 ? rnd(rng, -20, 20) : rnd(rng, -1000000000000LL, 1000000000000LL);
        const long long b = it % 5 == 0 ? rnd(rng, -20, 20) : rnd_ll(rng);
        Int128 e = 0;
        for (long long i = 0; i < n; ++i) e += math_internal::floor_div128((Int128)a * i + b, m);
        QUOT_CHECK(linear_floor_sum(n, m, a, b) == e);
    }
    for (int it = 0; it < 50000 * TEST_SCALE; ++it, ++cases) {  // 大きい n：1 つ前との差
        const long long n = rnd(rng, 1, 1000000000000000000LL), m = rnd(rng, 1, 1000000000000000000LL);
        const long long a = rnd(rng, -1000000000LL, 1000000000LL), b = rnd(rng, -1000000000000000000LL, 1000000000000000000LL);
        QUOT_CHECK(linear_floor_sum(n, m, a, b) - linear_floor_sum(n - 1, m, a, b) == math_internal::floor_div128((Int128)a * (n - 1) + b, m));
    }
    QUOT_CHECK(linear_floor_sum(1000000000, 1000000000, 999999999, 999999999) == 499999999500000000LL);
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long a = rnd(rng, 1, 30), b = rnd(rng, 1, 30), c = rnd(rng, -5, 300);
        long long e = 0;
        for (long long x = 0; a * x <= c; ++x)
            for (long long y = 0; a * x + b * y <= c; ++y) ++e;
        QUOT_CHECK(count_lattice_under_line(a, b, c) == e);
    }
    QUOT_CHECK(count_lattice_under_line(1, 1, 1000000000000000000LL) == (Int128)1000000000000000001LL * 1000000000000000002LL / 2);
    QUOT_CHECK(count_lattice_under_line(1, 1, LMAX) == ((Int128)1 << 125) + ((Int128)1 << 62));  // c = LLONG_MAX でも（(c+1)(c+2)/2）
    for (int it = 0; it < 3000; ++it, ++cases) {  // c が LLONG_MAX 近く：x の数が少ないように a を大きくして 1 列ずつ足す。a と b を入れかえても同じ
        const long long c = rnd(rng, LMAX - 1000000, LMAX), a = rnd(rng, c / 1000, c), b = it % 2 ? rnd(rng, 1, 1000) : rnd(rng, 1, LMAX);
        Int128 e = 0;
        for (long long x = 0; x <= c / a; ++x) e += (c - a * x) / b + 1;
        QUOT_CHECK(count_lattice_under_line(a, b, c) == e && count_lattice_under_line(b, a, c) == e);
    }
}

struct ModP {  // テスト用の mod 型（998244353。割り算も）
    static constexpr long long MD = 998244353;
    long long v;
    ModP(long long x = 0) : v(x % MD) {
        if (v < 0) v += MD;
    }
    friend ModP operator+(ModP a, ModP b) { return ModP(a.v + b.v); }
    friend ModP operator-(ModP a, ModP b) { return ModP(a.v - b.v); }
    friend ModP operator*(ModP a, ModP b) { return ModP(a.v * b.v % MD); }
    friend ModP operator/(ModP a, ModP b) { return a * ModP(mod_pow(b.v, MD - 2, MD)); }
    friend bool operator==(ModP a, ModP b) { return a.v == b.v; }
};
static ModP m128(Int128 x) { return ModP((long long)(((x % ModP::MD) + ModP::MD) % ModP::MD)); }
struct UeNode {  // 万能ユークリッドで Σ y・Σ x y・Σ y^2 を求めるための値（x = R の数、y = U の数）
    Int128 x = 0, y = 0, s = 0, sx = 0, sy2 = 0;
    friend UeNode operator*(const UeNode& a, const UeNode& b) {  // a のあとに b（b の各点の x, y は a の分だけずれる）
        UeNode r;
        r.x = a.x + b.x, r.y = a.y + b.y;
        r.s = a.s + b.s + a.y * b.x;
        r.sx = a.sx + b.sx + a.x * b.s + a.y * (b.x * (b.x + 1) / 2) + a.x * a.y * b.x;
        r.sy2 = a.sy2 + b.sy2 + 2 * a.y * b.s + a.y * a.y * b.x;
        return r;
    }
};
struct StrMonoid {  // 可換でない積の確かめ：文字列をつなげる
    std::string s;
    friend StrMonoid operator*(const StrMonoid& a, const StrMonoid& b) { return StrMonoid{a.s + b.s}; }
};

static void test_extensions(std::mt19937_64& rng, int& cases) {
    // for_each_quotient を途中で止める
    for (int it = 0; it < 2000; ++it, ++cases) {
        const long long N = rnd(rng, 0, 100000), stop = rnd(rng, 0, 50);
        long long calls = 0, all = 0;
        for_each_quotient(N, [&](long long, long long, long long) { ++all; });
        for_each_quotient(N, [&](long long, long long, long long) { return ++calls < stop; });
        QUOT_CHECK(calls == std::min(all, std::max(stop, 1LL)));
    }
    // Σ (a i + b) mod m
    for (int it = 0; it < 50000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 60), m = it % 3 == 0 ? rnd(rng, 1, 10) : rnd(rng, 1, it % 2 ? 1000000000LL : LMAX);
        const long long a = it % 4 == 0 ? rnd(rng, -20, 20) : rnd_ll(rng), b = rnd_ll(rng);
        Int128 e = 0;
        for (long long i = 0; i < n; ++i) {
            const Int128 v = ((Int128)a * i + b) % m;
            e += v < 0 ? v + m : v;
        }
        QUOT_CHECK(linear_mod_sum(n, m, a, b) == e);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい n：m ごとの周期で
        const long long m = rnd(rng, 1, 300), n = rnd(rng, 0, LMAX), a = rnd_ll(rng), b = rnd_ll(rng);
        Int128 per = 0, rest = 0;
        for (long long i = 0; i < m; ++i) {
            const long long v = mod_floor((long long)(((Int128)mod_floor(a, m) * i + mod_floor(b, m)) % m), m);
            per += v;
            if (i < n % m) rest += v;
        }
        QUOT_CHECK(linear_mod_sum(n, m, a, b) == per * (n / m) + rest);
    }
    QUOT_CHECK(linear_mod_sum(LMAX, LMAX, LMAX - 1, LMAX - 1) > 0 && linear_mod_sum(0, 5, 1, 1) == 0);
    // Σ⌊⌋・Σ i⌊⌋・Σ⌊⌋^2（mod p）
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 50), m = it % 3 == 0 ? rnd(rng, 1, 10) : rnd(rng, 1, it % 2 ? 1000000000LL : LMAX);
        const long long a = it % 4 == 0 ? rnd(rng, -20, 20) : rnd_ll(rng), b = it % 5 == 0 ? rnd(rng, -20, 20) : rnd_ll(rng);
        Int128 s = 0, si = 0;
        ModP sq = 0;
        for (long long i = 0; i < n; ++i) {
            const Int128 f = math_internal::floor_div128((Int128)a * i + b, m);
            s += f, si += f * i, sq = sq + m128(f) * m128(f);
        }
        const FloorSumExt<ModP> r = linear_floor_sum_ext<ModP>(n, m, a, b);
        QUOT_CHECK(r.sum == m128(s) && r.sum_i == m128(si) && r.sum_sq == sq);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい n：1 つ前との差
        const long long n = rnd(rng, 1, LMAX), m = rnd(rng, 1, LMAX), a = rnd_ll(rng), b = rnd_ll(rng);
        const FloorSumExt<ModP> r = linear_floor_sum_ext<ModP>(n, m, a, b), q = linear_floor_sum_ext<ModP>(n - 1, m, a, b);
        const ModP f = m128(math_internal::floor_div128((Int128)a * (n - 1) + b, m));
        QUOT_CHECK(r.sum - q.sum == f && r.sum_i - q.sum_i == f * ModP(n - 1) && r.sum_sq - q.sum_sq == f * f);
        if (m >= 1000000000LL && a >= -1000000000LL && a <= 1000000000LL) QUOT_CHECK(r.sum == m128(linear_floor_sum(n, m, a, b)));  // Int128 に収まる範囲で
    }
    // 万能ユークリッド：文字列（順番）と Σ の値
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, 30), a = rnd(rng, 0, 20), b = rnd(rng, 0, 40), c = rnd(rng, 1, 20);
        std::string e(std::string((size_t)(b / c), 'U'));
        for (long long x = 1; x <= n; ++x) e += std::string((size_t)((a * x + b) / c - (a * (x - 1) + b) / c), 'U') + "R";
        QUOT_CHECK(universal_euclid(n, a, b, c, StrMonoid{"U"}, StrMonoid{"R"}).s == e);
        UeNode U, R;
        U.y = 1, R.x = 1;
        const UeNode got = universal_euclid(n, a, b, c, U, R);
        Int128 s = 0, sx = 0, sy2 = 0;
        for (long long x = 1; x <= n; ++x) {
            const Int128 y = (a * x + b) / c;
            s += y, sx += x * y, sy2 += y * y;
        }
        QUOT_CHECK(got.x == n && got.y == (a * n + b) / c && got.s == s && got.sx == sx && got.sy2 == sy2);
    }
    {  // 溢れを調べる掛け算のモノイド：使わない最後の 2 乗をしないこと（2^40 は収まるが 2^64 は収まらない）
        struct ProdM {
            long long x = 1;
            ProdM operator*(const ProdM& o) const {
                ProdM r;
                if (__builtin_mul_overflow(x, o.x, &r.x)) math_internal::fail("ProdM", "溢れた");
                return r;
            }
        };
        ProdM U, R;
        U.x = 3, R.x = 2;
        QUOT_CHECK(universal_euclid(40, 0, 0, 1, U, R).x == (1LL << 40) && universal_euclid(10, 1, 0, 1, U, R).x == 60466176LL);  // 6^10 = 3^10 2^10
        cases += 2;
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい値：Σ y は linear_floor_sum と同じ
        const long long n = rnd(rng, 0, 1000000000LL), a = rnd(rng, 0, 1000000000LL), b = rnd(rng, 0, 1000000000LL), c = rnd(rng, 1, 1000000000LL);
        UeNode U, R;
        U.y = 1, R.x = 1;
        QUOT_CHECK(universal_euclid(n, a, b, c, U, R).s == linear_floor_sum(n + 1, c, a, b) - b / c);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_blocks(rng, cases);
    test_floor_sum(rng, cases);
    test_extensions(rng, cases);
    std::cout << "[quotient] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[quotient] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    const Int128 s1 = sum_floor_div(100000000000000LL);
    const long long t_sum = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Int128 s2 = sum_i_floor_div(100000000000000LL);
    const long long t_isum = ms(t0);
    std::vector<std::array<long long, 4>> qs(1000000);
    for (auto& q : qs) q = {rnd(rng, 1, 1000000000), rnd(rng, 1, 1000000000), rnd(rng, 0, 1000000000), rnd(rng, 0, 1000000000)};
    t0 = std::chrono::steady_clock::now();
    Int128 s3 = 0;
    for (const auto& q : qs) s3 += linear_floor_sum(q[0], q[1], q[2], q[3]);
    const long long t_fs = ms(t0);
    std::cout << "[quotient] benchmark: sum_floor_div(10^14) " << t_sum << " ms, sum_i_floor_div(10^14) " << t_isum << " ms, linear_floor_sum × 10^6 " << t_fs
              << " ms  (check " << (long long)((s1 + s2 + s3) % 7) << ")\n";
}

}  // namespace quot_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace quot_demo {

{COMMON_DEMO}

static void demo() {
    banner("商の列挙・床関数の和（quotient） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  for_each_quotient(N, f) / quotient_blocks(N) / quotient_blocks_ceil(N) / quotient_blocks2(N, M)   ⌊N/i⌋ が同じ i をまとめる\n"
              << "  sum_floor_div(N) / sum_i_floor_div(N) / sum_mod_i(N)        Σ⌊N/i⌋ / Σ i⌊N/i⌋ / Σ (N mod i)\n"
              << "  linear_floor_sum(n, m, a, b) / count_lattice_under_line(a, b, c)   Σ⌊(a i + b)/m⌋ / 直線の下の格子点\n";

    api("【1】 for_each_quotient(N, f) / quotient_blocks(N)", "なし（f(q, l, r) を呼ぶ）/ vector<QuotientBlock{q, l, r}>",
        "i = 1〜N を ⌊N/i⌋ = q が同じ区間 [l, r) に分ける（区間は 2√N 個以下）", "O(√N)");
    RUN(std::vector<std::string> v; for_each_quotient(10, [&](long long q, long long l, long long r) { v.push_back("q=" + std::to_string(q) + ":[" + std::to_string(l) + "," + std::to_string(r) + ")"); }));
    SHOW(v, "10/1 = 10、10/2 = 5、10/3 = 3、10/4 = 10/5 = 2、10/6〜10/10 = 1");
    SHOW(quotient_blocks(1000000000000LL).size(), "N = 10^12 でも 2×10^6 個くらい");
    SHOW(quotient_blocks_ceil(10).size(), "⌈10/i⌉：10, 5, 4, 3, 2（i = 5〜9）, 1（i = 10）");
    SHOW(quotient_blocks2(10, 7).size(), "⌊10/i⌋ と ⌊7/i⌋ が両方同じ区間（i ≤ 7）");
    {
        // ABC172 D：Σ_{K ≤ N} K·d(K) = Σ_i i · T(⌊N/i⌋)（T(m) = m(m+1)/2）
        auto solve = [](long long N) {
            Int128 s = 0;
            for_each_quotient(N, [&](long long q, long long l, long long r) { s += sum_range(l, r) * (Int128)(q * (q + 1) / 2); });
            return s;
        };
        note("ABC172 D（Sum of Divisors）: 区間 [l, r) の i の和（sum_range）× T(q)。O(√N)");
        std::cout << "    N = 4 → " << solve(4) << "、N = 100 → " << solve(100) << "、N = 10^7 → " << solve(10000000) << "\n";
    }
    {
        // ABC044 D：b ≤ √n は直接、b > √n は 2 桁で ⌊n/b⌋ = q が一定の区間ごとに b = (n − s)/q + 1
        auto solve = [](long long n, long long s) {
            if (s == n) return n + 1;  // b ≤ n なら 2 桁以上で桁和 < n なので、b = n + 1（1 桁）
            long long best = -1;
            for (long long b = 2; b * b <= n && best < 0; ++b) {
                long long t = 0;
                for (long long x = n; x; x /= b) t += x % b;
                if (t == s) best = b;
            }
            if (best >= 0) return best;
            for_each_quotient(n, [&](long long q, long long l, long long r) {
                if (best >= 0 || n < s || (n - s) % q != 0) return;
                const long long b = (n - s) / q + 1;  // 桁和 = q + (n − q b)
                if (l <= b && b < r && b > isqrt(n) && b >= 2) best = b;  // b·b は溢れることがあるので isqrt
            });
            return best;
        };
        note("ABC044 D（桁和）: b ≤ √n は直接。b > √n は ⌊n/b⌋ = q が同じ区間で n − q(b − 1) = s を解く");
        std::cout << "    (87654, 30) → " << solve(87654, 30) << "、(87654, 138) → " << solve(87654, 138) << "、(87654, 45678) → " << solve(87654, 45678)
                  << "、(31415926535, 1) → " << solve(31415926535LL, 1) << "、(1, 31415926535) → " << solve(1, 31415926535LL) << "\n";
    }

    api("【2】 sum_floor_div(N) / sum_i_floor_div(N) / sum_mod_i(N)", "Int128",
        "Σ⌊N/i⌋（約数の個数の和）/ Σ i⌊N/i⌋（約数の和の和）/ Σ (N mod i)", "O(√N)");
    SHOW(sum_floor_div(3), "ABC230 E: 3 + 1 + 1");
    SHOW(sum_floor_div(10000000000LL), "ABC230 E: N = 10^10");
    SHOW(sum_i_floor_div(4), "σ(1) + σ(2) + σ(3) + σ(4) = 1 + 3 + 4 + 7");
    SHOW(sum_mod_i(10), "0 + 0 + 1 + 2 + 0 + 4 + 3 + 2 + 1 + 0");

    api("【3】 linear_floor_sum(n, m, a, b)", "Int128", "Σ_{i=0}^{n−1} ⌊(a i + b) / m⌋（a, b は負でもよい。ACL の floor_sum と同じで範囲が広い）", "O(log m)");
    SHOW(linear_floor_sum(4, 10, 6, 3), "ACL Practice C 例 1");
    SHOW(linear_floor_sum(6, 5, 4, 3), "例 2");
    SHOW(linear_floor_sum(1, 1, 0, 0), "例 3");
    SHOW(linear_floor_sum(31415, 92653, 58979, 32384), "例 4");
    SHOW(linear_floor_sum(1000000000, 1000000000, 999999999, 999999999), "例 5");
    SHOW(linear_floor_sum(5, 3, -2, 1), "負の係数も：⌊1/3⌋ + ⌊−1/3⌋ + ⌊−3/3⌋ + ⌊−5/3⌋ + ⌊−7/3⌋ = 0 − 1 − 1 − 2 − 3");
    note("引数の順番は ACL の floor_sum(n, m, a, b) と同じ（個数 n、割る数 m、傾き a、切片 b）");
    SHOW(linear_mod_sum(5, 7, 3, 1), "Σ (3i + 1) mod 7 = 1 + 4 + 0 + 3 + 6");
    note("(a i + b) mod m の和は「Σ (a i + b) − m Σ ⌊(a i + b)/m⌋」なので linear_mod_sum で O(log m)");

    api("【4】 count_lattice_under_line(a, b, c)", "Int128", "a x + b y ≤ c となる x, y ≥ 0 の組の数（三角形の中の格子点）", "O(log)");
    SHOW(count_lattice_under_line(2, 3, 6), "(0,0) (0,1) (0,2) (1,0) (1,1) (2,0) (3,0)");
    SHOW(count_lattice_under_line(1, 1, 1000000000000LL), "(10^12 + 1)(10^12 + 2) / 2");

    api("【5】 linear_floor_sum_ext<mint>(n, m, a, b)", "FloorSumExt<mint>{sum, sum_i, sum_sq}",
        "Σ⌊(a i + b)/m⌋・Σ i⌊…⌋・Σ⌊…⌋^2（i = 0〜n−1）を mod で（mod は奇数の素数）", "O(log m)");
    {
        const auto r = linear_floor_sum_ext<quot_test::ModP>(4, 3, 2, 1);  // ⌊1/3⌋, ⌊3/3⌋, ⌊5/3⌋, ⌊7/3⌋ = 0, 1, 1, 2
        std::cout << "  linear_floor_sum_ext<mint>(4, 3, 2, 1) -> {" << r.sum.v << ", " << r.sum_i.v << ", " << r.sum_sq.v
                  << "}   // 値 0, 1, 1, 2：和 4、Σ i·値 = 0 + 1 + 2 + 6 = 9、Σ 値^2 = 6\n";
    }

    api("【6】 universal_euclid(n, a, b, c, U, R)", "M（自作の型）",
        "直線 y = ⌊(a x + b)/c⌋ に沿って、y が 1 増えるたびに U、x が 1 増えるたびに R を掛けた積（x = 1〜n。M は結合的な * と単位元 M()）",
        "O(log × M の掛け算 × log)");
    note("Σ_x A^x B^{⌊(ax+b)/c⌋}（行列）のような床関数を含む和に。M に「U の数・R の数・欲しい和」を持たせて * を書く");
    {
        const auto s = universal_euclid(4, 2, 1, 3, quot_test::StrMonoid{"U"}, quot_test::StrMonoid{"R"}).s;
        std::cout << "  universal_euclid(4, 2, 1, 3, \"U\", \"R\") -> " << s << "   // y(x) = ⌊(2x+1)/3⌋ = 1, 1, 2, 3（x = 1〜4）\n";
    }
}

}  // namespace quot_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    quot_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    quot_test::benchmark();    // 速度計測
    quot_demo::demo();         // 使い方デモ
    return 0;
}
