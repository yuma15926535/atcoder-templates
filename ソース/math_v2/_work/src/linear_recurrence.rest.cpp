// =============================================================================
//  ランダムテスト
// =============================================================================
namespace linrec_test {

{COMMON_TEST}

template <class mint>
static std::vector<mint> naive_terms(const std::vector<mint>& a, const std::vector<mint>& c, size_t n) {
    std::vector<mint> r(a.begin(), a.begin() + (long long)std::min(a.size(), n));
    r.resize(n);
    for (size_t i = a.size(); i < n; ++i) {
        mint s = 0;
        for (size_t j = 0; j < c.size(); ++j) s += c[j] * r[i - 1 - j];
        r[i] = s;
    }
    return r;
}
template <class mint>
static std::vector<mint> rand_vec(std::mt19937_64& rng, size_t n) {
    std::vector<mint> v(n);
    for (auto& x : v) x = rng() % 4 == 0 ? mint(rnd(rng, -2, 2)) : mint((long long)(rng() >> 2));
    return v;
}

template <class mint>
static void test_kth(std::mt19937_64& rng, int& cases, int iters) {
    for (int it = 0; it < iters; ++it) {
        const size_t d = (size_t)(it % 3 == 0 ? rnd(rng, 0, 6) : rnd(rng, 1, it % 2 ? 40 : 200));
        const std::vector<mint> c = rand_vec<mint>(rng, d), a = naive_terms(rand_vec<mint>(rng, d), c, d + (size_t)rnd(rng, 0, 2));  // 多めに渡してもよい
        const size_t N = d * 3 + 50;
        const std::vector<mint> t = naive_terms(a, c, N);
        for (int q = 0; q < 6; ++q, ++cases) {
            const long long K = rnd(rng, 0, (long long)N - 1);
            LINREC_CHECK(linear_recurrence_kth(a, c, K) == t[(size_t)K]);
        }
        ++cases;
        LINREC_CHECK(linear_recurrence_terms(a, c, (long long)N) == t);
        // K が巨大：t 項ずらした列で同じ値になるか
        const long long K = rnd(rng, 0, 1000000000000000000LL), sh = rnd(rng, 0, (long long)N - (long long)d);
        const std::vector<mint> a2(t.begin() + sh, t.begin() + sh + (long long)d);
        if (K >= sh) LINREC_CHECK(linear_recurrence_kth(a2, c, K - sh) == linear_recurrence_kth(a, c, K));
    }
}

template <class mint>
static void test_bm(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 80 * TEST_SCALE; ++it, ++cases) {
        const size_t d = (size_t)rnd(rng, 0, it % 2 ? 10 : 80);
        const std::vector<mint> c = rand_vec<mint>(rng, d), a = rand_vec<mint>(rng, d);
        const std::vector<mint> s = naive_terms(a, c, 2 * d + (size_t)rnd(rng, 0, 5));
        const std::vector<mint> g = berlekamp_massey(s);
        LINREC_CHECK(g.size() <= d);
        const std::vector<mint> t = naive_terms(a, c, 2 * d + 100);
        LINREC_CHECK(naive_terms(std::vector<mint>(s.begin(), s.begin() + (long long)g.size()), g, t.size()) == t);
        const long long K = rnd(rng, 0, 1000000000000000000LL);
        LINREC_CHECK(guess_kth(s, K) == linear_recurrence_kth(a, c, K));
    }
    // n^3 は (1 − x)^4 の漸化式
    std::vector<mint> cube;
    for (long long n = 0; n < 12; ++n) cube.push_back(mint(n * n * n));
    LINREC_CHECK(berlekamp_massey(cube) == std::vector<mint>({mint(4), mint(-6), mint(4), mint(-1)}));
    LINREC_CHECK(guess_kth(cube, 1000000) == mint(1000000).pow(3));
    LINREC_CHECK(berlekamp_massey(std::vector<mint>(10, mint(0))).empty() && berlekamp_massey(std::vector<mint>()).empty());
    cases += 3;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_kth<ModInt998>(rng, cases, 60 * TEST_SCALE);
    test_kth<ModInt107>(rng, cases, 25 * TEST_SCALE);  // NTT が使えない素数
    DynModInt::set_mod(1000000000);
    test_kth<DynModInt>(rng, cases, 15 * TEST_SCALE);  // 合成数の mod でも第 K 項・最初の n 項は使える
    DynModInt::set_mod(998244353);
    test_kth<DynModInt>(rng, cases, 8 * TEST_SCALE);
    test_bm<ModInt998>(rng, cases);
    test_bm<ModInt107>(rng, cases);
    {  // フィボナッチ数：倍数公式と
        using mint = ModInt998;
        auto fib = [](unsigned long long n) {  // (F_n, F_{n+1})
            mint x = 0, y = 1;
            for (int b = 63; b >= 0; --b) {
                const mint x2 = x * (y * mint(2) - x), y2 = x * x + y * y;
                x = x2, y = y2;
                if (n >> b & 1) {
                    const mint t = x + y;
                    x = y, y = t;
                }
            }
            return x;
        };
        for (int it = 0; it < 1000; ++it, ++cases) {
            const long long K = rnd(rng, 0, 1000000000000000000LL);
            LINREC_CHECK(linear_recurrence_kth<mint>({0, 1}, {1, 1}, K) == fib((unsigned long long)K));
        }
    }
    {  // 長い列（O(n log n) の方）：漸化式どおりの値と
        using mint = ModInt998;
        for (int it = 0; it < 2; ++it, ++cases) {
            const size_t d = 100, N = 310000;
            const std::vector<mint> c = rand_vec<mint>(rng, d), a = rand_vec<mint>(rng, d);
            LINREC_CHECK(linear_recurrence_terms(a, c, (long long)N) == naive_terms(a, c, N));
        }
        const std::vector<ModInt107> c7 = rand_vec<ModInt107>(rng, 120), a7 = rand_vec<ModInt107>(rng, 120);
        LINREC_CHECK(linear_recurrence_terms(a7, c7, 260000) == naive_terms(a7, c7, 260000));
    }
    std::cout << "[linear_recurrence] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[linear_recurrence] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    using mint = ModInt998;
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::vector<mint> c = rand_vec<mint>(rng, 50000), a = rand_vec<mint>(rng, 50000);
    auto t0 = std::chrono::steady_clock::now();
    const mint x = linear_recurrence_kth(a, c, 1000000000000000000LL);
    const long long t_k = ms(t0);
    const std::vector<ModInt107> c7 = rand_vec<ModInt107>(rng, 10000), a7 = rand_vec<ModInt107>(rng, 10000);
    t0 = std::chrono::steady_clock::now();
    const ModInt107 y = linear_recurrence_kth(a7, c7, 1000000000000000000LL);
    const long long t_k7 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<mint> t = linear_recurrence_terms(a, c, 1000000);
    const long long t_t = ms(t0);
    const std::vector<mint> s = linear_recurrence_terms(std::vector<mint>(a.begin(), a.begin() + 2000), std::vector<mint>(c.begin(), c.begin() + 2000), 4000);
    t0 = std::chrono::steady_clock::now();
    const std::vector<mint> g = berlekamp_massey(s);
    const long long t_bm = ms(t0);
    std::cout << "[linear_recurrence] benchmark: kth（d = 5×10^4, K = 10^18）" << t_k << " ms, kth（10^9+7, d = 10^4）" << t_k7 << " ms, terms（d = 5×10^4, n = 10^6）"
              << t_t << " ms, berlekamp_massey（4000 項）" << t_bm << " ms  (check " << (x.val() + y.val() + t.back().val() + g.size()) % 7 << ")\n";
}

}  // namespace linrec_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace linrec_demo {

{COMMON_DEMO}

static void demo() {
    using mint = ModInt998;
    banner("線形漸化式（linear_recurrence） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  linear_recurrence_kth(a, c, K) / linear_recurrence_terms(a, c, n) / berlekamp_massey(s) / guess_kth(s, K)\n"
              << "  a_i = c_1 a_{i−1} + … + c_d a_{i−d}、a = {a_0, …, a_{d−1}}、c = {c_1, …, c_d}（この中では using mint = ModInt998;）\n";

    api("【1】 linear_recurrence_kth(a, c, K)", "mint", "第 K 項（K ≤ 10^18。Bostan–Mori。mod は何でも）", "O(d log d log K)");
    SHOW(linear_recurrence_kth<mint>({0, 1}, {1, 1}, 10), "フィボナッチ数 F_10");
    SHOW(linear_recurrence_kth<mint>({0, 1}, {1, 1}, 1000000000000000000LL), "F_{10^18} mod 998244353");
    SHOW(linear_recurrence_kth<mint>({0, 0, 1}, {1, 1, 1}, 30), "トリボナッチ数");
    SHOW(linear_recurrence_kth<ModInt107>({1, 3}, {4, -1}, 1000000000000000000LL), "a_i = 4a_{i−1} − a_{i−2}（mod 10^9+7）");

    api("【2】 linear_recurrence_terms(a, c, n)", "vector<mint>", "最初の n 項", "O(n log n)（d が小さいときは O(n d)）");
    SHOW(linear_recurrence_terms<mint>({0, 1}, {1, 1}, 12), "F_0〜F_11");
    SHOW(linear_recurrence_terms<mint>({1, 1, 1}, {0, 1, 1}, 12), "パドヴァン数列 a_i = a_{i−2} + a_{i−3}");

    api("【3】 berlekamp_massey(s) / guess_kth(s, K)", "vector<mint> / mint",
        "列をみたす最短の漸化式の c（次数 d には 2d 項以上いる）/ それで第 K 項", "O(N^2) + O(d log d log K)");
    SHOW(berlekamp_massey(std::vector<mint>({0, 1, 1, 2, 3, 5, 8, 13})), "フィボナッチ：c = {1, 1}");
    SHOW(berlekamp_massey(std::vector<mint>({0, 1, 8, 27, 64, 125, 216, 343, 512, 729})), "n^3：c = {4, −6, 4, −1}（(1 − x)^4）");
    {
        // 3 × n の盤面をドミノで敷き詰める数：列ごとの DP で最初の 40 項を出し、guess_kth で n = 10^18
        std::vector<mint> f;
        for (int n = 0; n < 40; ++n) {
            std::vector<mint> dp(8, mint(0));  // 列の境目で、次の列にはみ出しているマス（3 ビット）
            dp[0] = 1;
            for (int col = 0; col < n; ++col) {
                std::vector<mint> nd(8, mint(0));
                for (int s = 0; s < 8; ++s) {
                    if (dp[(size_t)s] == mint(0)) continue;
                    std::function<void(int, int, int)> fill = [&](int r, int cur, int nxt) {  // r 行目から置く。cur は埋まった、nxt ははみ出し
                        if (r == 3) { nd[(size_t)nxt] += dp[(size_t)s]; return; }
                        if (cur >> r & 1) { fill(r + 1, cur, nxt); return; }
                        fill(r + 1, cur | 1 << r, nxt | 1 << r);  // 横に置く
                        if (r + 1 < 3 && !(cur >> (r + 1) & 1)) fill(r + 2, cur | 3 << r, nxt);  // 縦に置く
                    };
                    fill(0, s, 0);
                }
                dp.swap(nd);
            }
            f.push_back(dp[0]);
        }
        note("3 × n をドミノで敷き詰める数：DP で最初の 40 項 → berlekamp_massey → guess_kth で n = 10^18");
        std::cout << "    最初の項 " << to_s(std::vector<mint>(f.begin(), f.begin() + 9)) << " …\n";
        std::cout << "    漸化式 c = " << to_s(berlekamp_massey(f)) << "（a_n = 4a_{n−2} − a_{n−4}）\n";
        std::cout << "    n = 10^18 → " << guess_kth(f, 1000000000000000000LL) << "\n";
    }
    {
        // ABC300 Ex（Fibonacci: Revisited）：Σ_{m ⊆ N} a_m。P/Q の Bostan–Mori で、N のビットが 1 なら偶奇両方の項を足す
        auto solve = [](int K, long long N) {
            std::vector<mint> Q((size_t)K + 1), A((size_t)K, mint(1));  // Q = 1 − x − … − x^K、P = A Q mod x^K
            Q[0] = 1;
            for (int j = 1; j <= K; ++j) Q[(size_t)j] = -1;
            std::vector<mint> P = convolve(A, Q);
            P.resize((size_t)K);
            for (; N > 0; N >>= 1) {
                std::vector<mint> Qm(Q);
                for (size_t i = 1; i < Qm.size(); i += 2) Qm[i] = -Qm[i];
                const std::vector<mint> U = convolve(P, Qm), V = convolve(Q, Qm);
                std::vector<mint> nP((size_t)K, mint(0)), nQ((size_t)K + 1, mint(0));
                for (size_t i = 0; i < U.size(); ++i)
                    if (i % 2 == 0 || (N & 1)) nP[i / 2] += U[i];  // N のビットが 0 なら偶数次だけ、1 なら両方
                for (size_t i = 0; i < V.size(); i += 2) nQ[i / 2] = V[i];
                P.swap(nP), Q.swap(nQ);
            }
            return P[0] / Q[0];
        };
        note("ABC300 Ex（Fibonacci: Revisited）: Bostan–Mori で、N のビットが 1 の段は偶数次と奇数次の項を両方足す（convolve で）");
        std::cout << "    (2, 6) → " << solve(2, 6) << "、(2, 8) → " << solve(2, 8) << "、(1, 123456789) → " << solve(1, 123456789) << "、(300, 20230429) → "
                  << solve(300, 20230429) << "\n    (42923, 999999999558876113) → " << solve(42923, 999999999558876113LL) << "\n";
    }
}

}  // namespace linrec_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    linrec_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    linrec_test::benchmark();    // 速度計測
    linrec_demo::demo();         // 使い方デモ
    return 0;
}
