// =============================================================================
//  ランダムテスト
// =============================================================================
namespace ser_test {

{COMMON_TEST}

template <class mint>
static void test_sum_powers(std::mt19937_64& rng, int& cases, int kmax) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int k = (int)rnd(rng, 0, kmax);
        const long long n = rnd(rng, 0, 150);
        mint e = 0;
        for (long long i = 0; i < n; ++i) e += mint(i).pow(k);  // 0^0 = 1
        SER_CHECK(sum_powers<mint>(n, k) == e);
        const long long l = rnd(rng, 0, 150), r = rnd(rng, 0, 150);
        mint e2 = 0;
        for (long long i = l; i < r; ++i) e2 += mint(i).pow(k);
        SER_CHECK(sum_powers<mint>(l, r, k) == e2 && sum_powers_closed<mint>(l, r - 1, k) == e2);
    }
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {  // 大きい n：1 つ前との差、k = 1, 2 の公式
        const long long n = rnd(rng, 1, 1000000000000000000LL);
        const int k = (int)rnd(rng, 1, kmax);
        SER_CHECK(sum_powers<mint>(n + 1, k) - sum_powers<mint>(n, k) == mint(n).pow(k));
        const mint N(n);
        SER_CHECK(sum_powers<mint>(n, 1) * mint(2) == N * (N - mint(1)));
        SER_CHECK(sum_powers<mint>(n, 2) * mint(6) == (N - mint(1)) * N * (mint(2) * N - mint(1)));
    }
}

template <class mint>
static void test_arith_geom(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const mint r = it % 5 == 0 ? mint(rnd(rng, -2, 2)) : mint((long long)rng());
        const long long n = rnd(rng, 0, 60);
        mint e = 0, p = 1;
        for (long long i = 0; i < n; ++i) e += mint(i) * p, p *= r;
        SER_CHECK(arith_geom_sum(r, n) == e);
        const long long N = rnd(rng, 0, 1000000000000000000LL);
        SER_CHECK(arith_geom_sum(r, N + 1) - arith_geom_sum(r, N) == mint(N) * r.pow(N));
    }
}

template <class mint>
static void test_interpolate(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 5000 * TEST_SCALE; ++it, ++cases) {
        const int d = (int)rnd(rng, 0, 25);
        std::vector<mint> c((size_t)d + 1), ys((size_t)d + 1);
        for (auto& x : c) x = (long long)rng();
        for (int i = 0; i <= d; ++i) ys[(size_t)i] = poly_eval(c, mint(i));
        const long long x = it % 3 == 0 ? rnd(rng, -5, 40) : rnd_ll(rng);
        SER_CHECK(lagrange_interpolate(ys, x) == poly_eval(c, mint(x)));
        // ばらばらの点
        const int n = (int)rnd(rng, 1, 20);
        std::vector<mint> xs, vs;
        while ((int)xs.size() < n) {
            const mint v = it % 2 ? mint(rnd(rng, -30, 30)) : mint((long long)rng());
            if (std::find(xs.begin(), xs.end(), v) == xs.end()) xs.push_back(v), vs.push_back(mint((long long)rng()));
        }
        const std::vector<mint> co = poly_from_points(xs, vs);
        SER_CHECK((int)co.size() == n);
        for (int i = 0; i < n; ++i) SER_CHECK(poly_eval(co, xs[(size_t)i]) == vs[(size_t)i]);
        const mint t = (long long)rng();
        SER_CHECK(lagrange_interpolate_at(xs, vs, t) == poly_eval(co, t));
    }
    SER_CHECK(lagrange_interpolate(std::vector<mint>{}, 5) == mint(0) && lagrange_interpolate(std::vector<mint>{mint(7)}, 1000000000000000000LL) == mint(7));
}

static void test_poly(std::mt19937_64& rng, int& cases) {
    using mint = ModInt998;
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 12), m = (int)rnd(rng, 1, 8);
        std::vector<mint> a((size_t)n), b((size_t)m);
        for (auto& x : a) x = it % 3 ? mint((long long)rng()) : mint(rnd(rng, 0, 2));
        for (auto& x : b) x = it % 3 ? mint((long long)rng()) : mint(rnd(rng, 0, 2));
        if (b.back() == mint(0)) b.back() = 1;
        // 掛け算
        const std::vector<mint> c = poly_mul_naive(a, b);
        const mint x = (long long)rng();
        if (!a.empty()) SER_CHECK(poly_eval(c, x) == poly_eval(a, x) * poly_eval(b, x) && c.size() == a.size() + b.size() - 1);
        // 割り算：a = q b + r、deg r < deg b
        const PolyDivMod<mint> d = poly_divmod(a, b);
        std::vector<mint> back = poly_mul_naive(d.q, b);
        back.resize(std::max(back.size(), a.size()), mint(0));
        for (size_t i = 0; i < d.r.size(); ++i) back[i] += d.r[i];
        while (!back.empty() && back.back() == mint(0)) back.pop_back();
        std::vector<mint> at = a;
        while (!at.empty() && at.back() == mint(0)) at.pop_back();
        SER_CHECK(back == at && d.r.size() < b.size() && (d.r.empty() || d.r.back() != mint(0)));
        // 整数：b の最高次の係数が ±1 か、全部が割り切れるように作る
        std::vector<long long> qi((size_t)rnd(rng, 1, 6)), bi((size_t)m), ri((size_t)(m - 1));
        for (auto& v : qi) v = rnd(rng, -100, 100);
        for (auto& v : bi) v = rnd(rng, -100, 100);
        for (auto& v : ri) v = rnd(rng, -100, 100);
        bi.back() = it % 2 ? 1 : -1;
        if (qi.back() == 0) qi.back() = 3;
        std::vector<long long> ai = poly_mul_naive(qi, bi);
        for (size_t i = 0; i < ri.size(); ++i) ai[i] += ri[i];
        const PolyDivMod<long long> di = poly_divmod(ai, bi);
        while (!ri.empty() && ri.back() == 0) ri.pop_back();
        SER_CHECK(di.q == qi && di.r == ri);
    }
    // double
    const PolyDivMod<double> dd = poly_divmod<double>({-1, 0, 0, 1}, {-1, 1});  // x^3 − 1 = (x − 1)(x^2 + x + 1)
    SER_CHECK(dd.q == std::vector<double>({1, 1, 1}) && dd.r.empty());
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_sum_powers<ModInt998>(rng, cases, 12);
    test_sum_powers<ModInt107>(rng, cases, 12);
    test_sum_powers<ModInt<101>>(rng, cases, 12);
    test_arith_geom<ModInt998>(rng, cases);
    DynModInt::set_mod(1000000000);  // 合成数の mod でも使える
    test_arith_geom<DynModInt>(rng, cases);
    test_interpolate<ModInt998>(rng, cases);
    test_interpolate<ModInt107>(rng, cases);
    test_poly(rng, cases);
    std::cout << "[series] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[series] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    using mint = ModInt107;
    auto t0 = std::chrono::steady_clock::now();
    const mint s1 = sum_powers<mint>(1000000000000000000LL, 2500000);
    const long long t_sp = ms(t0);
    std::vector<mint> ys(1000001);
    for (auto& y : ys) y = (long long)rng();
    t0 = std::chrono::steady_clock::now();
    const mint s2 = lagrange_interpolate(ys, 1000000000000000000LL);
    const long long t_li = ms(t0);
    std::vector<mint> xs(2000), vs(2000);
    for (int i = 0; i < 2000; ++i) xs[(size_t)i] = mint(i * 7 + 3), vs[(size_t)i] = (long long)rng();
    t0 = std::chrono::steady_clock::now();
    const std::vector<mint> co = poly_from_points(xs, vs);
    const long long t_pp = ms(t0);
    std::vector<mint> a(4000), b(2000);
    for (auto& x : a) x = (long long)rng();
    for (auto& x : b) x = (long long)rng();
    t0 = std::chrono::steady_clock::now();
    const PolyDivMod<mint> d = poly_divmod(a, b);
    const long long t_div = ms(t0);
    std::cout << "[series] benchmark: sum_powers(10^18, k = 2.5×10^6) " << t_sp << " ms, lagrange_interpolate（10^6 + 1 点）" << t_li << " ms, poly_from_points（2000 点）" << t_pp
              << " ms, poly_divmod（4000 ÷ 2000）" << t_div << " ms  (check " << (s1 + s2 + co[5] + d.q[3]).val() % 7 << ")\n";
}

}  // namespace ser_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace ser_demo {

{COMMON_DEMO}

static void demo() {
    banner("べき乗和・ラグランジュ補間・小さい多項式（series） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  sum_powers<mint>(n, k) / sum_powers<mint>(l, r, k) / sum_powers_closed<mint>(L, R, k)   べき乗和\n"
              << "  arith_geom_sum(r, n)                       Σ i·r^i\n"
              << "  lagrange_interpolate(ys, x) / lagrange_interpolate_at(xs, ys, x) / poly_from_points(xs, ys)   補間\n"
              << "  poly_eval / poly_mul_naive / poly_divmod   次数の小さい多項式\n";
    note("using mint = ModInt107;");
    using mint = ModInt107;

    api("【1】 sum_powers<mint>(n, k) / sum_powers<mint>(l, r, k) / sum_powers_closed<mint>(L, R, k)", "mint",
        "Σ i^k（[0, n) / [l, r) / [L, R]）。k + 2 点の値を作ってラグランジュ補間", "O(k)");
    SHOW(sum_powers<mint>(11, 2), "0^2 + 1^2 + … + 10^2 = 385");
    SHOW(sum_powers_closed<mint>(1, 1000000000000000000LL, 3), "1^3 + … + (10^18)^3");
    SHOW(sum_powers<mint>(5, 8, 1), "5 + 6 + 7");
    SHOW(sum_powers<mint>(3, 0), "0^0 + 1^0 + 2^0（0^0 = 1）");

    api("【2】 arith_geom_sum(r, n)", "mint", "Σ_{i=0}^{n−1} i·r^i（等差 × 等比。mod は何でもよい）", "O(log n)");
    SHOW(arith_geom_sum(mint(2), 4), "0 + 2 + 8 + 24 = 34");
    SHOW(arith_geom_sum(mint(3), 1000000000000000000LL), "n = 10^18");

    api("【3】 lagrange_interpolate(ys, x)", "mint", "f(0)〜f(d) が ys の d 次以下の多項式の f(x)（x は 10^18 でもよい）", "O(d)");
    SHOW(lagrange_interpolate(std::vector<mint>{1, 3, 7}, 3), "ARC033 D 例 1：f(x) = x^2 + x + 1 → f(3) = 13");
    SHOW(lagrange_interpolate(std::vector<mint>{4, 16, 106, 484, 1624, 4384}, 1000000000), "ARC033 D 例 2");
    {
        // ABC208 F：f(n, m) は n の K + M 次の多項式。n = 0〜K+M の値を DP で作って補間
        auto solve = [](long long N, int M, int K) {
            const int D = K + M + 1;
            std::vector<mint> f((size_t)D);
            for (int n = 1; n < D; ++n) f[(size_t)n] = mint(n).pow(K);  // m = 0
            for (int m = 1; m <= M; ++m)
                for (int n = 1; n < D; ++n) f[(size_t)n] += f[(size_t)n - 1];  // f(n, m) = Σ_{j ≤ n} f(j, m − 1)
            return lagrange_interpolate(f, N);
        };
        note("ABC208 F（Cumulative Sum）: n = 0〜K+M の f(n, M) を累積和で作り、lagrange_interpolate(f, N)");
        std::cout << "    (3, 4, 2) → " << solve(3, 4, 2) << "、(0, 1, 2) → " << solve(0, 1, 2) << "、(10^18, 30, 123456) → " << solve(1000000000000000000LL, 30, 123456) << "\n";
    }

    api("【4】 lagrange_interpolate_at(xs, ys, x) / poly_from_points(xs, ys)", "mint / vector<mint>",
        "ばらばらの点 (xs[i], ys[i]) を通る多項式の f(x) / 係数（低い次数から）", "O(d^2)");
    RUN(std::vector<mint> xs = {1, 2, 4}, ys = {2, 5, 17});
    SHOW(poly_from_points(xs, ys), "f(x) = x^2 + 1 → 係数 [1, 0, 1]");
    SHOW(lagrange_interpolate_at(xs, ys, mint(10)), "f(10) = 101");

    api("【5】 poly_eval(c, x) / poly_mul_naive(a, b) / poly_divmod(a, b)", "T / vector<T> / PolyDivMod{q, r}",
        "係数は低い次数から。T は mint・long long・double。割り算は a = q b + r（整数なら割り切れないと止める）", "O(n) / O(nm)");
    SHOW(poly_eval<long long>({1, 0, 1}, 3), "x^2 + 1 で x = 3");
    SHOW(poly_mul_naive<long long>({2, 1}, {6, 4, 2}), "(x + 2)(2x^2 + 4x + 6)");
    RUN(PolyDivMod<long long> d = poly_divmod<long long>({5, 0, 0, 1}, {1, 1}));
    SHOW(d.q, "(x^3 + 5) ÷ (x + 1) の商 x^2 − x + 1");
    SHOW(d.r, "余り 4");
    note("ABC245 D（Polynomial division）: B = poly_divmod(C, A).q（係数は整数）");
    std::cout << "    例 1 → " << to_s(poly_divmod<long long>({12, 14, 8, 2}, {2, 1}).q) << "、例 2 → " << to_s(poly_divmod<long long>({10000, 0, -1}, {100, 1}).q) << "\n";
}

}  // namespace ser_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    ser_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ser_test::benchmark();    // 速度計測
    ser_demo::demo();         // 使い方デモ
    return 0;
}
