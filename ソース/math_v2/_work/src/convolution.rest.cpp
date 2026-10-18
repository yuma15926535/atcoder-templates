// =============================================================================
//  ランダムテスト
// =============================================================================
namespace conv_test {

{COMMON_TEST}

template <class mint>
static std::vector<mint> naive(const std::vector<mint>& a, const std::vector<mint>& b) {
    if (a.empty() || b.empty()) return {};
    std::vector<mint> c(a.size() + b.size() - 1);
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j) c[i + j] += a[i] * b[j];
    return c;
}
static size_t rand_len(std::mt19937_64& rng, int it) {
    switch (it % 6) {
        case 0: return (size_t)rnd(rng, 0, 5);
        case 1: return (size_t)rnd(rng, 55, 70);  // 愚直と NTT の境目
        case 2: return (size_t)rnd(rng, 1, 300);
        case 3: return ((size_t)1 << rnd(rng, 0, 9)) + (size_t)rnd(rng, -1, 1);  // 2 の冪の前後
        default: return (size_t)rnd(rng, 1, TEST_SCALE == 1 ? 400 : 1500);  // デバッグビルドでは短く
    }
}

template <class mint>
static void test_mint(std::mt19937_64& rng, int& cases, int iters) {
    for (int it = 0; it < iters; ++it, ++cases) {
        std::vector<mint> a(rand_len(rng, it)), b(rand_len(rng, it + 3));
        for (auto& x : a) x = (long long)rng();
        for (auto& x : b) x = it % 5 == 0 ? mint(mint::mod() - 1) : mint((long long)rng());  // 大きい値の連続
        CONV_CHECK(convolve(a, b) == naive(a, b));
        if (it % 4 == 0) CONV_CHECK(convolve(a, a) == naive(a, a));  // 2 乗（同じ列）
    }
}

static void test_dynamic(std::mt19937_64& rng, int& cases) {
    const long long mods[] = {998244353, 469762049, 1000000007, 1000000000, 1 << 30, 2147483647, 65537, 7, 2, 1, 754974721};
    for (long long m : mods) {
        DynModInt::set_mod(m);
        test_mint<DynModInt>(rng, cases, 150 * TEST_SCALE);
    }
}

static void test_integer(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const size_t n = rand_len(rng, it), m = rand_len(rng, it + 1);
        const int kind = it % 3;
        const long long R = kind == 0 ? 1000000 : (kind == 1 ? 1000000000 : (1LL << 55));
        std::vector<long long> a(n), b(m);
        for (auto& x : a) x = rnd(rng, -R, R);
        for (auto& x : b) x = rnd(rng, -R, R);
        if (it % 7 == 0)
            for (auto& x : a) x = R;  // 端の値をそろえる
        std::vector<Int128> e(n && m ? n + m - 1 : 0, 0);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < m; ++j) e[i + j] += (Int128)a[i] * b[j];
        CONV_CHECK(convolve_i128(a, b) == e);
        if (kind == 0) {
            std::vector<long long> el(e.begin(), e.end());
            CONV_CHECK(convolve_ll(a, b) == el);
        }
        const long long md = it % 2 ? rnd(rng, 1, 2147483647) : 998244353;
        std::vector<long long> em(e.size());
        for (size_t i = 0; i < e.size(); ++i) em[i] = (long long)(((e[i] % md) + md) % md);
        CONV_CHECK(convolve_mod(a, b, md) == em);
    }
}

static void test_others(std::mt19937_64& rng, int& cases) {
    using mint = ModInt998;
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        // 多くの多項式の積
        std::vector<std::vector<mint>> ps((size_t)rnd(rng, 0, 12));
        for (auto& p : ps) {
            p.resize((size_t)rnd(rng, 1, it % 2 ? 5 : 200));
            for (auto& x : p) x = (long long)rng();
        }
        std::vector<mint> e = {mint(1)};
        for (const auto& p : ps) e = naive(e, p);
        CONV_CHECK(convolve_many(ps) == e);
        // 相関
        std::vector<mint> a((size_t)rnd(rng, 1, 300)), b((size_t)rnd(rng, 1, (long long)a.size()));
        for (auto& x : a) x = (long long)rng();
        for (auto& x : b) x = (long long)rng();
        std::vector<mint> c(a.size() - b.size() + 1);
        for (size_t k = 0; k < c.size(); ++k)
            for (size_t i = 0; i < b.size(); ++i) c[k] += a[k + i] * b[i];
        CONV_CHECK(correlate(a, b) == c);
        // ntt → intt で元に戻る
        std::vector<mint> f((size_t)1 << rnd(rng, 0, 12));
        for (auto& x : f) x = (long long)rng();
        std::vector<mint> g = f;
        ntt(g), intt(g);
        CONV_CHECK(g == f);
    }
    CONV_CHECK(convolve_many<mint>({}) == std::vector<mint>{mint(1)} && convolve(std::vector<mint>{}, std::vector<mint>{mint(1)}).empty());
    // 巡回畳み込み
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        std::vector<mint> a((size_t)rnd(rng, 0, 200)), b((size_t)rnd(rng, 0, 200));
        for (auto& x : a) x = (long long)rng();
        for (auto& x : b) x = (long long)rng();
        const size_t n = (size_t)rnd(rng, 1, 250);
        std::vector<mint> c(n, 0);
        for (size_t i = 0; i < a.size(); ++i)
            for (size_t j = 0; j < b.size(); ++j) c[(i + j) % n] += a[i] * b[j];
        CONV_CHECK(convolve_cyclic(a, b, n) == c);
    }
    // オンラインの畳み込み：各時点の c_n が全体の畳み込みと同じ
    for (int it = 0; it < 30 * TEST_SCALE; ++it) {
        const size_t N = (size_t)rnd(rng, 1, it % 3 ? 300 : 3000);
        std::vector<mint> a(N), b(N);
        for (auto& x : a) x = (long long)rng();
        for (auto& x : b) x = (long long)rng();
        const std::vector<mint> full = convolve(a, b);
        RelaxedConvolution<mint> rc;
        for (size_t n = 0; n < N; ++n, ++cases) CONV_CHECK(rc.push(a[n], b[n]) == full[n] && rc.size() == n + 1 && rc.get(n) == full[n]);
    }
    {  // 答えを次の入力に使う：カタラン数 C_{n+1} = Σ C_i C_{n−i}
        RelaxedConvolution<ModInt107> rc;
        std::vector<ModInt107> C = {1};
        for (int n = 0; n < 2000; ++n) C.push_back(rc.push(C[(size_t)n], C[(size_t)n]));
        ModInt107 e = 1;  // C_2000 = (4000)! / (2001! 2000!) を直接
        for (int i = 2001; i <= 4000; ++i) e *= ModInt107(i);
        for (int i = 1; i <= 2000; ++i) e /= ModInt107(i);
        e /= ModInt107(2001);
        CONV_CHECK(C[2000] == e);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_mint<ModInt998>(rng, cases, 600 * TEST_SCALE);
    test_mint<ModInt107>(rng, cases, 400 * TEST_SCALE);
    test_mint<ModInt<167772161>>(rng, cases, 200 * TEST_SCALE);
    test_dynamic(rng, cases);
    test_integer(rng, cases);
    test_others(rng, cases);
    std::cout << "[convolution] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[convolution] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 500000;
    std::vector<ModInt998> a(n), b(n);
    for (auto& x : a) x = (long long)rng();
    for (auto& x : b) x = (long long)rng();
    auto t0 = std::chrono::steady_clock::now();
    const auto c = convolve(a, b);
    const long long t998 = ms(t0);
    std::vector<ModInt107> a7(n), b7(n);
    for (auto& x : a7) x = (long long)rng();
    for (auto& x : b7) x = (long long)rng();
    t0 = std::chrono::steady_clock::now();
    const auto c7 = convolve(a7, b7);
    const long long t107 = ms(t0);
    std::vector<long long> la(200000), lb(200000);
    for (auto& x : la) x = rnd(rng, -1000000000, 1000000000);
    for (auto& x : lb) x = rnd(rng, -1000000000, 1000000000);
    t0 = std::chrono::steady_clock::now();
    const auto cl = convolve_i128(la, lb);
    const long long tll = ms(t0);
    std::vector<std::vector<ModInt998>> ps(100000, std::vector<ModInt998>{1, 1});
    for (auto& p : ps) p[1] = (long long)rng();
    t0 = std::chrono::steady_clock::now();
    const auto pm = convolve_many(ps);
    const long long tmany = ms(t0);
    t0 = std::chrono::steady_clock::now();
    RelaxedConvolution<ModInt998> rc;
    ModInt998 last = 1;
    for (int i = 0; i < 200000; ++i) last = rc.push(last + a[(size_t)i], b[(size_t)i]);
    const long long trc = ms(t0);
    std::cout << "[convolution] benchmark: convolve<ModInt998>（5×10^5 + 5×10^5）" << t998 << " ms, convolve<ModInt107>（同じ。3 つの素数）" << t107
              << " ms, convolve_i128（2×10^5、|値| ≤ 10^9）" << tll << " ms\n"
              << "[convolution] benchmark: convolve_many（(1 + a x) を 10^5 個）" << tmany << " ms, RelaxedConvolution（2×10^5 回 push）" << trc << " ms  (check "
              << (c[12345] + pm[777] + last).val() % 7 + c7[999].val() % 7 + (long long)((cl[4321] % 7 + 7) % 7) << ")\n";
}

}  // namespace conv_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace conv_demo {

{COMMON_DEMO}

static void demo() {
    banner("畳み込み（convolution） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  convolve(a, b)                 mint の列（どの mod でも。NTT 素数でなければ自動で 3 つの素数）\n"
              << "  convolve_ll / convolve_i128    整数のまま厳密に\n"
              << "  convolve_mod(a, b, m)          long long の列を mod m で\n"
              << "  convolve_many(ps) / correlate(a, b) / ntt / intt\n";
    note("using mint = ModInt998;");
    using mint = ModInt998;

    api("【1】 convolve(a, b)", "vector<mint>", "c[k] = Σ_{i+j=k} a[i] b[j]（多項式の積）。mod は何でもよい", "O((n+m) log(n+m))");
    SHOW(convolve<mint>({1, 2, 3, 4}, {5, 6, 7, 8, 9}), "ACL Practice F の例 1");
    SHOW(convolve<mint>({10000000}, {10000000}), "例 2（10^14 mod 998244353）");
    SHOW(convolve<ModInt107>({1000000006, 2}, {3, 4}), "10^9+7 でも（−1·3、−4 + 6、8）");
    note("ATC001 C（高速フーリエ変換）: 定食の値段 k の組の数 = A と B の畳み込み（mod なしなら convolve_ll）");
    SHOW(convolve_ll({0, 1, 2, 3, 4}, {0, 1, 2, 4, 8}), "値段 i の数を i 番目に（0 円は 0 個）→ 1〜8 番目が 0 1 4 11 26 36 40 32");

    api("【2】 convolve_ll(a, b) / convolve_i128(a, b) / convolve_mod(a, b, m)", "vector<long long> / vector<Int128>",
        "整数のまま厳密に（負の数も。必要なだけの素数を使う）/ mod m（1 ≤ m < 2^31）で", "O((n+m) log(n+m))");
    SHOW(convolve_ll({-1, 2}, {3, -4}), "(−1 + 2x)(3 − 4x)");
    SHOW(convolve_i128({1000000000000000000LL, 1}, {1000000000000000000LL, 1}), "10^36 は long long に入らないので Int128 で");
    SHOW(convolve_mod({5, 6}, {7, 8}, 10), "35, 88, 48 を mod 10");

    api("【3】 convolve_many(ps) / correlate(a, b)", "vector<mint>",
        "多くの多項式の積（短いものから掛ける）/ c[k] = Σ_i a[k+i] b[i]（パターン b を a の上でずらして重ねる）", "O(N log^2 N) / O(n log n)");
    SHOW(convolve_many<mint>({{1, 1}, {1, 2}, {1, 3}}), "(1 + x)(1 + 2x)(1 + 3x)");
    SHOW(correlate<mint>({1, 2, 3, 4}, {1, 1}), "隣り合う 2 つの和：3, 5, 7");
    {
        // ABC196 F：T を S のどこに重ねると違う文字が最も少ないか。一致数 = 1 どうしの相関 + 0 どうしの相関
        auto solve = [](const std::string& S, const std::string& T) {
            std::vector<mint> s1(S.size()), s0(S.size()), t1(T.size()), t0(T.size());
            for (size_t i = 0; i < S.size(); ++i) s1[i] = S[i] == '1', s0[i] = S[i] == '0';
            for (size_t i = 0; i < T.size(); ++i) t1[i] = T[i] == '1', t0[i] = T[i] == '0';
            const auto c1 = correlate(s1, t1), c0 = correlate(s0, t0);
            long long best = (long long)T.size();
            for (size_t k = 0; k < c1.size(); ++k) best = std::min(best, (long long)T.size() - (c1[k] + c0[k]).val());
            return best;
        };
        note("ABC196 F（Substring 2）: 一致する文字の数 = correlate(S の 1, T の 1) + correlate(S の 0, T の 0)");
        std::cout << "    例 1 → " << solve("0001", "101") << "、例 2 → " << solve("0101010", "1010101") << "、例 3 → " << solve("10101000010011011110", "0010011111") << "\n";
    }
    {
        // ABC291 G：A を回して Σ (A_i | B_i) を最大に。ビットごとに「両方 0」の数を相関で（A を 2 周並べる）
        auto solve = [](const std::vector<int>& A, const std::vector<int>& B) {
            const size_t n = A.size();
            std::vector<long long> total(n, 0);
            for (int bit = 0; bit < 5; ++bit) {
                std::vector<mint> a(2 * n), b(n);
                for (size_t i = 0; i < 2 * n; ++i) a[i] = !((A[i % n] >> bit) & 1);
                for (size_t i = 0; i < n; ++i) b[i] = !((B[i] >> bit) & 1);
                const auto c = correlate(a, b);
                for (size_t k = 0; k < n; ++k) total[k] += ((long long)n - c[k].val()) << bit;
            }
            return *std::max_element(total.begin(), total.end());
        };
        note("ABC291 G（OR Sum）: ビットごとに「両方 0」の位置の数を correlate で（A は 2 周並べる）");
        std::cout << "    例 1 → " << solve({0, 1, 3}, {0, 2, 3}) << "、例 2 → " << solve({1, 6, 1, 4, 3}, {0, 6, 4, 0, 1}) << "\n";
    }

    api("【4】 ntt(a) / intt(a)", "なし（a を書きかえる）", "長さ 2^h の NTT（mod は 998244353 などの NTT 素数）。結果の順はビットを逆にした順なので、点ごとの積に使う", "O(n log n)");
    RUN(std::vector<mint> f = {1, 2, 0, 0}, g = {3, 4, 0, 0}; ntt(f); ntt(g); for (int i = 0; i < 4; ++i) f[i] *= g[i]; intt(f));
    SHOW(f, "(1 + 2x)(3 + 4x) = 3 + 10x + 8x^2");

    api("【5】 convolve_cyclic(a, b, n) / RelaxedConvolution<mint> rc;  rc.push(a_n, b_n)", "vector<mint> / mint",
        "巡回畳み込み（添字を mod n で）/ オンラインの畳み込み：a_n, b_n を 1 つずつ渡すと c_n = Σ_{i+j=n} a_i b_j が返る（答えを次の入力に使う DP に）",
        "O(n log n) / 全体 O(N log^2 N)");
    SHOW(convolve_cyclic<mint>({1, 2, 3}, {1, 1, 0}, 3), "(1 + 2x + 3x^2)(1 + x) を x^3 = 1 で：{1+3, 2+1, 3+2}");
    {
        RelaxedConvolution<mint> rc;
        std::vector<mint> C = {1};  // カタラン数 C_{n+1} = Σ_{i=0}^{n} C_i C_{n−i}
        for (int n = 0; n < 10; ++n) C.push_back(rc.push(C[(size_t)n], C[(size_t)n]));
        note("f_{n+1} = Σ_{i=0}^{n} f_i g_{n−i} の形：rc.push(f_n, g_n) が f_{n+1} になる（例はカタラン数）");
        std::cout << "    C_0〜C_10 → ";
        for (const mint& x : C) std::cout << x << " ";
        std::cout << "\n";
    }
}

}  // namespace conv_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    conv_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    conv_test::benchmark();    // 速度計測
    conv_demo::demo();         // 使い方デモ
    return 0;
}
