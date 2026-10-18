// =============================================================================
//  ランダムテスト
// =============================================================================
namespace fps_test {

{COMMON_TEST}

template <class mint>
static std::vector<mint> rand_vec(std::mt19937_64& rng, size_t n) {
    std::vector<mint> v(n);
    for (auto& x : v) x = rng() % 5 == 0 ? mint(rnd(rng, 0, 3)) : mint((long long)rng());
    return v;
}
template <class mint>
static std::vector<mint> mul_trunc(const std::vector<mint>& a, const std::vector<mint>& b, size_t n) {
    std::vector<mint> c(n);
    for (size_t i = 0; i < a.size() && i < n; ++i)
        for (size_t j = 0; j < b.size() && i + j < n; ++j) c[i + j] += a[i] * b[j];
    return c;
}
template <class mint>
static std::vector<mint> brute_inv(const std::vector<mint>& f, size_t n) {  // g_k = −(Σ_{i≥1} f_i g_{k−i}) / f_0
    std::vector<mint> g(n);
    const mint i0 = f[0].inv();
    for (size_t k = 0; k < n; ++k) {
        mint s = k == 0 ? mint(1) : mint(0);
        for (size_t i = 1; i <= k && i < f.size(); ++i) s -= f[i] * g[k - i];
        g[k] = s * i0;
    }
    return g;
}
template <class mint>
static std::vector<mint> brute_exp(const std::vector<mint>& f, size_t n) {  // g_k = (1/k) Σ_{i=1}^{k} i f_i g_{k−i}
    std::vector<mint> g(n);
    if (n) g[0] = 1;
    for (size_t k = 1; k < n; ++k) {
        mint s = 0;
        for (size_t i = 1; i <= k && i < f.size(); ++i) s += mint((long long)i) * f[i] * g[k - i];
        g[k] = s / mint((long long)k);
    }
    return g;
}
template <class mint>
static std::vector<mint> brute_log(const std::vector<mint>& f, size_t n) {  // f L' = f'
    std::vector<mint> d(n), L(n);
    for (size_t k = 0; k + 1 < n; ++k) {
        mint s = (k + 1 < f.size()) ? f[k + 1] * mint((long long)(k + 1)) : mint(0);
        for (size_t i = 1; i <= k && i < f.size(); ++i) s -= f[i] * d[k - i];
        d[k] = s;
        L[k + 1] = d[k] / mint((long long)(k + 1));
    }
    return L;
}

template <class mint>
static void test_basic(std::mt19937_64& rng, int& cases, int iters, bool prime_mod) {
    for (int it = 0; it < iters; ++it, ++cases) {
        const size_t n = (size_t)(it % 3 == 0 ? rnd(rng, 0, 8) : rnd(rng, 1, it % 2 ? 80 : 400));
        std::vector<mint> f = rand_vec<mint>(rng, (size_t)rnd(rng, 1, (long long)n + 3));
        if (f[0] == mint(0)) f[0] = 1;
        if (!prime_mod)
            while (!f[0].is_unit()) f[0] += 1;
        FPS_CHECK(fps_inv(f, (int)n) == Fps<mint>(brute_inv(f, n)));
        if (!prime_mod) continue;
        // log・exp
        std::vector<mint> g = f;
        g[0] = 1;
        FPS_CHECK(fps_log(g, (int)n) == Fps<mint>(brute_log(g, n)));
        g[0] = 0;
        FPS_CHECK(fps_exp(g, (int)n) == Fps<mint>(brute_exp(g, n)));
        // 累乗（先頭に 0 があってもよい）
        std::vector<mint> h(f);
        if (it % 3 == 0) h.insert(h.begin(), (size_t)rnd(rng, 1, 4), mint(0));
        const long long k = rnd(rng, 0, 5);
        std::vector<mint> e(n);
        if (n) e[0] = 1;
        for (long long t = 0; t < k; ++t) e = mul_trunc(e, h, n);
        FPS_CHECK(fps_pow(h, k, (int)n) == Fps<mint>(e));
        const long long k1 = rnd(rng, 0, 1000000000000000000LL), k2 = rnd(rng, 0, 1000);
        if (it % 4 == 0) FPS_CHECK(Fps<mint>(mul_trunc<mint>(fps_pow(f, k1, (int)n), fps_pow(f, k2, (int)n), n)) == fps_pow(f, k1 + k2, (int)n));
        // 平方根：g^2 の平方根の 2 乗がもとに戻る（先頭の 0 は偶数個）
        std::vector<mint> sq = mul_trunc(h, h, n + 6);
        const Maybe<Fps<mint>> s = fps_sqrt(sq, (int)n);
        FPS_CHECK(s.ok && Fps<mint>(mul_trunc<mint>(s.value, s.value, n)) == Fps<mint>(sq).pre(n));
        // 微分・積分
        FPS_CHECK(fps_derivative(fps_integral(f)) == Fps<mint>(f));
    }
}

template <class mint>
static void test_poly(std::mt19937_64& rng, int& cases, int iters) {
    for (int it = 0; it < iters; ++it, ++cases) {
        // 割り算
        std::vector<mint> f = rand_vec<mint>(rng, (size_t)rnd(rng, 0, 300)), g = rand_vec<mint>(rng, (size_t)rnd(rng, 1, 150));
        if (g.back() == mint(0)) g.back() = 1;
        const FpsDivMod<mint> d = fps_divmod(f, g);
        Fps<mint> back = Fps<mint>(g) * d.q + d.r;
        while (!back.empty() && back.back() == mint(0)) back.pop_back();
        while (!f.empty() && f.back() == mint(0)) f.pop_back();
        FPS_CHECK(back == Fps<mint>(f) && d.r.size() < g.size() && (d.r.empty() || d.r.back() != mint(0)));
        // テイラーシフト：いくつかの点で f(x + c) と比べる
        const Fps<mint> p(rand_vec<mint>(rng, (size_t)rnd(rng, 0, 200)));
        const mint c = (long long)rng();
        const Fps<mint> sh = taylor_shift(p, c);
        for (int t = 0; t < 3; ++t) {
            const mint x = (long long)rng();
            FPS_CHECK(sh.eval(x) == p.eval(x + c));
        }
        // P/Q の K 次の係数：級数を展開したものと
        std::vector<mint> P = rand_vec<mint>(rng, (size_t)rnd(rng, 0, 20)), Q = rand_vec<mint>(rng, (size_t)rnd(rng, 1, 20));
        if (Q[0] == mint(0)) Q[0] = 1;
        const size_t K = (size_t)rnd(rng, 0, 300);
        const std::vector<mint> ser = mul_trunc(P, brute_inv(Q, K + 1), K + 1);
        FPS_CHECK(fps_coeff_of_fraction(P, Q, (long long)K) == ser[K]);
        if (it % 4 == 0) {  // 次数が大きい（NTT で Q(−x) を作る方）・P の方が長い
            std::vector<mint> P2 = rand_vec<mint>(rng, (size_t)rnd(rng, 0, 300)), Q2 = rand_vec<mint>(rng, (size_t)rnd(rng, 16, 200));
            Q2[0] = 1;
            const size_t K2 = (size_t)rnd(rng, 0, 1500);
            const std::vector<mint> ser2 = mul_trunc(P2, brute_inv(Q2, K2 + 1), K2 + 1);
            FPS_CHECK(fps_coeff_of_fraction(P2, Q2, (long long)K2) == ser2[K2]);
        }
        // Π (1 + a_i x)
        const std::vector<mint> as = rand_vec<mint>(rng, (size_t)rnd(rng, 0, 100));
        std::vector<mint> e = {mint(1)};
        for (const mint& a : as) e = mul_trunc<mint>(e, {mint(1), a}, e.size() + 1);
        FPS_CHECK(prod_linear(as) == Fps<mint>(e));
    }
    // フィボナッチ：x / (1 − x − x^2) の K 次の係数を倍化の公式と
    for (int it = 0; it < 200; ++it, ++cases) {
        const long long K = rnd(rng, 0, 1000000000000000000LL);
        mint a = 0, b = 1;  // (F_k, F_{k+1}) を上の桁から
        for (int bit = 62; bit >= 0; --bit) {
            const mint c2 = a * (b * 2 - a), d2 = a * a + b * b;
            a = c2, b = d2;
            if ((K >> bit) & 1) {
                const mint t = a + b;
                a = b, b = t;
            }
        }
        FPS_CHECK(fps_coeff_of_fraction<mint>({0, 1}, {1, -1, -1}, K) == a);
    }
}

static void test_edges(int& cases) {
    using mint = ModInt998;
    FPS_CHECK(fps_inv<mint>({mint(2)}, 0).empty() && fps_pow<mint>({0, 0, 3}, 2, 4) == Fps<mint>{0, 0, 0, 0} && fps_pow<mint>({0, 0, 3}, 1, 4) == Fps<mint>{0, 0, 3, 0});
    FPS_CHECK(fps_pow<mint>({}, 0, 3) == Fps<mint>{1, 0, 0} && fps_pow<mint>({}, 5, 3) == Fps<mint>{0, 0, 0});
    FPS_CHECK(!fps_sqrt<mint>({0, 1}, 5).ok && !fps_sqrt<mint>({3}, 2).ok && fps_sqrt<mint>({0, 0, 4}, 3).value.pre(2) * fps_sqrt<mint>({0, 0, 4}, 3).value.pre(2) == Fps<mint>{0, 0, 4});  // 3 は平方剰余でない
    FPS_CHECK(fps_coeff_of_fraction<mint>({}, {1}, 5) == mint(0) && fps_coeff_of_fraction<mint>({7}, {1}, 0) == mint(7));
    {  // 小さい素数の mod：項の数 < mod なら exp・pow も使える
        DynModInt::set_mod(7);
        using D = DynModInt;
        const Fps<D> e = fps_exp(Fps<D>{0, 1}, 5);
        FPS_CHECK(e == Fps<D>{1, 1, D(2).inv(), D(6).inv(), D(24).inv()});
        FPS_CHECK(fps_pow(Fps<D>{1, 1}, 3, 5) == Fps<D>{1, 3, 3, 1, 0} && fps_exp(Fps<D>{0, 3}, 6).size() == 6);
        DynModInt::set_mod(998244353);
    }
    cases += 6;
}

template <class mint>
static void test_points(std::mt19937_64& rng, int& cases, int iters) {
    for (int it = 0; it < iters; ++it, ++cases) {
        // 多点評価（小さいときは 1 点ずつ、大きいときは木）
        const int n = (int)rnd(rng, 0, it % 4 ? 60 : 3000), m = (int)rnd(rng, 0, it % 4 ? 60 : 3000);
        Fps<mint> f((size_t)n);
        for (auto& x : f) x = (long long)rng();
        std::vector<mint> xs((size_t)m);
        for (auto& x : xs) x = it % 3 ? mint((long long)rng()) : mint(rnd(rng, 0, 5));  // 同じ点が何度も出てもよい
        const std::vector<mint> ev = multipoint_eval<mint>(f, xs);
        for (int i = 0; i < m; i += std::max(1, m / 50)) FPS_CHECK(ev[(size_t)i] == f.eval(xs[(size_t)i]));
        // 補間（異なる点）
        const int k = (int)rnd(rng, 0, it % 4 ? 70 : 2000);
        std::vector<mint> px, py;
        std::set<long long> used;
        while ((int)px.size() < k) {
            const long long v = it % 2 ? rnd(rng, 0, 1000000) : (long long)(rng() % mint::mod());
            if (used.insert(v).second) px.push_back(mint(v)), py.push_back(mint((long long)rng()));
        }
        const Fps<mint> P = interpolate(px, py);
        FPS_CHECK((int)P.size() == k);
        for (int i = 0; i < k; i += std::max(1, k / 50)) FPS_CHECK(P.eval(px[(size_t)i]) == py[(size_t)i]);
        // べき乗和
        const int na = (int)rnd(rng, 0, it % 4 ? 30 : 1000), L = (int)rnd(rng, 0, 40);
        std::vector<mint> as((size_t)na), ws((size_t)na);
        for (auto& x : as) x = (long long)rng();
        for (auto& x : ws) x = (long long)rng();
        const Fps<mint> p1 = power_sums(as, L), pw = power_sums(as, L, ws);
        FPS_CHECK((int)p1.size() == L && (int)pw.size() == L);
        for (int j = 0; j < L; j += std::max(1, L / 10)) {
            mint e1 = 0, ew = 0;
            for (int i = 0; i < na; ++i) e1 += as[(size_t)i].pow(j), ew += ws[(size_t)i] * as[(size_t)i].pow(j);
            FPS_CHECK(p1[(size_t)j] == e1 && pw[(size_t)j] == ew);
        }
    }
    // 分数の和：Σ P_i/Q_i を 1 つにまとめた N/D が、それぞれの和と同じ級数になる
    for (int it = 0; it < 50 * TEST_SCALE; ++it, ++cases) {
        const int c = (int)rnd(rng, 1, 8), L = 20;
        std::vector<Fps<mint>> nums((size_t)c), dens((size_t)c);
        Fps<mint> e((size_t)L, mint(0));
        for (int i = 0; i < c; ++i) {
            nums[(size_t)i].resize((size_t)rnd(rng, 0, 4)), dens[(size_t)i].resize((size_t)rnd(rng, 1, 4));
            for (auto& x : nums[(size_t)i]) x = (long long)rng();
            for (auto& x : dens[(size_t)i]) x = (long long)rng();
            dens[(size_t)i][0] = rnd(rng, 1, 100);
            e += (nums[(size_t)i] * fps_inv(dens[(size_t)i], L)).pre(L);
        }
        const FpsFraction<mint> fr = fps_sum_of_fractions(nums, dens);
        FPS_CHECK((fr.num * fps_inv(fr.den, L)).pre(L) == e);
    }
    {
        const std::vector<mint> xs = {1, 2, 3};
        FPS_CHECK(interpolate<mint>(xs, {2, 4, 6}) == Fps<mint>({0, 2, 0}) && multipoint_eval<mint>({}, xs) == std::vector<mint>(3, mint(0)));
        FPS_CHECK(interpolate<mint>({}, {}).empty() && power_sums<mint>({}, 3) == Fps<mint>(3, mint(0)));
        cases += 2;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_points<ModInt998>(rng, cases, 60 * TEST_SCALE);
    test_points<ModInt107>(rng, cases, 20 * TEST_SCALE);
    test_basic<ModInt998>(rng, cases, 300 * TEST_SCALE, true);
    test_basic<ModInt107>(rng, cases, 150 * TEST_SCALE, true);  // NTT が使えない素数（3 つの素数の畳み込み）
    DynModInt::set_mod(1000000000);
    test_basic<DynModInt>(rng, cases, 100 * TEST_SCALE, false);  // 合成数の mod でも fps_inv は使える
    test_poly<ModInt998>(rng, cases, 200 * TEST_SCALE);
    test_poly<ModInt107>(rng, cases, 60 * TEST_SCALE);
    test_edges(cases);
    std::cout << "[fps] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[fps] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    using mint = ModInt998;
    const int n = 500000;
    std::vector<mint> f(n);
    for (auto& x : f) x = (long long)rng();
    f[0] = 1;
    auto t0 = std::chrono::steady_clock::now();
    const auto gi = fps_inv(f, n);
    const long long t_inv = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto gl = fps_log(f, n);
    const long long t_log = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto ge = fps_exp(gl, n);
    const long long t_exp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto gp = fps_pow(f, 1000000000000000000LL, n);
    const long long t_pow = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto gs = fps_sqrt(f, n);
    const long long t_sqrt = ms(t0);
    std::vector<mint> Q(100001);
    for (auto& x : Q) x = (long long)rng();
    Q[0] = 1;
    t0 = std::chrono::steady_clock::now();
    const mint cf = fps_coeff_of_fraction<mint>({1}, Q, 1000000000000000000LL);
    const long long t_bm = ms(t0);
    std::vector<mint> xs(100000), ys(100000);
    for (size_t i = 0; i < xs.size(); ++i) xs[i] = mint((long long)i * 7 + 3), ys[i] = (long long)rng();
    const std::vector<mint> f1(f.begin(), f.begin() + 100000);
    t0 = std::chrono::steady_clock::now();
    const auto ev = multipoint_eval(f1, xs);
    const long long t_mp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto ip = interpolate(xs, ys);
    const long long t_ip = ms(t0);
    std::cout << "[fps] benchmark: n = 5×10^5 で inv " << t_inv << " ms, log " << t_log << " ms, exp " << t_exp << " ms, pow(k = 10^18) " << t_pow << " ms, sqrt " << t_sqrt
              << " ms, coeff_of_fraction（次数 10^5、K = 10^18）" << t_bm << " ms\n"
              << "[fps] benchmark: multipoint_eval（次数 10^5、10^5 点）" << t_mp << " ms, interpolate（10^5 点）" << t_ip << " ms  (check "
              << (gi[5] + ge[7] + gp[9] + cf + ev[77] + ip[55] + (gs.ok ? gs.value[3] : mint(0))).val() % 7 << (ge == Fps<mint>(f) ? "" : " exp(log f) ≠ f") << ")\n";
}

}  // namespace fps_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace fps_demo {

{COMMON_DEMO}

static void demo() {
    banner("形式的冪級数（fps） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Fps<mint>：+ − × / f.pre(n) / f.eval(x)\n"
              << "  fps_inv / fps_log / fps_exp / fps_pow / fps_sqrt（最初の n 項）\n"
              << "  fps_divmod / fps_derivative / fps_integral / taylor_shift\n"
              << "  fps_coeff_of_fraction(P, Q, K) / prod_linear(as)\n";
    note("using mint = ModInt998;");
    using mint = ModInt998;

    api("【1】 Fps<mint> f = {…} / f * g / fps_inv(f, n)", "Fps<mint>", "係数は低い次数から。× は畳み込み。1/f の最初の n 項（f[0] ≠ 0）", "O(n log n)");
    RUN(Fps<mint> f = {1, -1, -1});
    SHOW(fps_inv(f, 10), "1/(1 − x − x^2) = フィボナッチ数 1, 1, 2, 3, 5, …");
    SHOW(f * f, "(1 − x − x^2)^2");
    {
        // 分割数 p(n)：1/Π(1 − x^k) の係数
        const int N = 101;
        Fps<mint> e(N);
        e[0] = 1;
        for (long long j = 1;; ++j) {  // Π(1 − x^k) = 1 + Σ_{j≥1} (−1)^j (x^{j(3j−1)/2} + x^{j(3j+1)/2})（五角数定理）
            const long long a = j * (3 * j - 1) / 2, b = j * (3 * j + 1) / 2;
            if (a >= N) break;
            e[a] += (j % 2 ? mint(-1) : mint(1));
            if (b < N) e[b] += (j % 2 ? mint(-1) : mint(1));
        }
        const Fps<mint> p = fps_inv(e, N);
        note("分割数 p(n) = [x^n] 1/Π(1 − x^k)（Π は五角数定理で作る）");
        std::cout << "    p(10) = " << p[10] << "、p(100) = " << p[100] << "\n";
    }

    api("【2】 fps_log(f, n) / fps_exp(f, n) / fps_pow(f, k, n) / fps_sqrt(f, n)", "Fps<mint> / Maybe<Fps<mint>>",
        "log f（f[0] = 1）/ exp f（f[0] = 0）/ f^k（k ≤ 10^18）/ √f（ないときは ok = false）", "O(n log n)");
    {
        // ベル数：exp(e^x − 1) の係数 × n!
        const int N = 11;
        Fps<mint> ex(N);
        mint fact = 1;
        for (int i = 1; i < N; ++i) fact *= mint(i), ex[i] = fact.inv();  // e^x − 1
        const Fps<mint> b = fps_exp(ex, N);
        mint f10 = 1;
        for (int i = 1; i <= 10; ++i) f10 *= mint(i);
        note("ベル数 B_n = n! [x^n] exp(e^x − 1)");
        std::cout << "    B_10 = " << b[10] * f10 << "（115975）\n";
    }
    SHOW(fps_pow<mint>({1, 1}, 5, 7), "(1 + x)^5");
    SHOW(fps_pow<mint>({0, 1, 1}, 1000000000000000000LL, 3), "(x + x^2)^{10^18} は x^{10^18} から始まるので最初の 3 項は 0");
    {
        // カタラン数：C(x) = (1 − √(1 − 4x)) / (2x)
        const Maybe<Fps<mint>> s = fps_sqrt<mint>({1, -4}, 12);
        Fps<mint> c(11);
        for (int i = 0; i < 11; ++i) c[i] = -s.value[i + 1] / 2;
        note("カタラン数 = (1 − √(1 − 4x)) / (2x) の係数（fps_sqrt）");
        std::cout << "    " << to_s(std::vector<mint>(c.begin(), c.end())) << "\n";
    }
    SHOW(fps_sqrt<mint>({0, 1}, 3).ok, "x の平方根は冪級数にならない");

    api("【3】 fps_divmod(f, g) / fps_derivative(f) / fps_integral(f) / taylor_shift(f, c)", "FpsDivMod{q, r} / Fps<mint>",
        "多項式の割り算 f = q g + r / 微分 / 積分 / f(x + c)", "O(n log n)");
    RUN(FpsDivMod<mint> d = fps_divmod<mint>({5, 0, 0, 1}, {1, 1}));
    SHOW(d.q, "(x^3 + 5) ÷ (x + 1) の商 x^2 − x + 1（−1 は mod で 998244352）");
    SHOW(d.r, "余り 4");
    SHOW(fps_derivative<mint>({1, 2, 3}), "1 + 2x + 3x^2 の微分");
    SHOW(taylor_shift<mint>({0, 0, 1}, mint(1)), "x^2 → (x + 1)^2 = 1 + 2x + x^2");

    api("【4】 fps_coeff_of_fraction(P, Q, K) / prod_linear(as)", "mint / Fps<mint>",
        "[x^K] P/Q（線形漸化式の第 K 項。K ≤ 10^18）/ Π (1 + a_i x)", "O(d log d log K) / O(n log^2 n)");
    SHOW(fps_coeff_of_fraction<mint>({0, 1}, {1, -1, -1}, 10), "フィボナッチ F_10 = 55（x/(1 − x − x^2)）");
    SHOW(fps_coeff_of_fraction<mint>({0, 1}, {1, -1, -1}, 1000000000000000000LL), "F_{10^18} mod 998244353");
    SHOW(prod_linear<mint>({1, 2, 3}), "(1 + x)(1 + 2x)(1 + 3x)：基本対称式");
    note("a_i を 1 つ以上選んだ積の和など、k 個選ぶ積の和は prod_linear の k 次の係数");

    api("【5】 multipoint_eval(f, xs) / interpolate(xs, ys)", "vector<mint> / Fps<mint>",
        "f(x_0), …, f(x_{m−1}) をまとめて / 点 (x_i, y_i)（x は異なる）を通る次数 n − 1 以下の多項式", "O(n log^2 n)");
    SHOW(multipoint_eval<mint>({1, 0, 1}, {0, 1, 2, 3}), "1 + x^2 の x = 0, 1, 2, 3 での値");
    SHOW(interpolate<mint>({0, 1, 2}, {1, 2, 5}), "(0,1) (1,2) (2,5) を通る：1 + x^2");
    note("x が 0, 1, …, n−1 の等間隔で 1 点の値だけ欲しいなら 06_series の lagrange_interpolate（O(n)）が速い");

    api("【6】 power_sums(as, n, ws) / fps_sum_of_fractions(nums, dens)", "Fps<mint> / FpsFraction<mint>{num, den}",
        "p_k = Σ w_i a_i^k（k = 0〜n−1。ws を省くと 1）/ Σ P_i/Q_i を 1 つの分数に（分割統治で通分）", "O(m log^2 m + n log n)");
    SHOW(power_sums<mint>({1, 2, 3}, 4), "Σ a^k：3, 6, 14, 36");
    SHOW(power_sums<mint>({2, 3}, 3, {10, 1}), "10·2^k + 3^k：11, 23, 49");
    note("Σ_i a_i^k をすべての k = 0〜n−1 について一度に（1/(1 − a_i x) の和の係数）");
}

}  // namespace fps_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    fps_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    fps_test::benchmark();    // 速度計測
    fps_demo::demo();         // 使い方デモ
    return 0;
}
