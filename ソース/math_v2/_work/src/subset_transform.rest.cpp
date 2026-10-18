// =============================================================================
//  ランダムテスト
// =============================================================================
namespace subt_test {

{COMMON_TEST}

template <class T, class Gen>
static void test_type(std::mt19937_64& rng, int& cases, int iters, Gen gen) {
    for (int it = 0; it < iters; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 7);
        const size_t N = (size_t)1 << n;
        std::vector<T> a(N), b(N);
        for (auto& x : a) x = gen(rng);
        for (auto& x : b) x = gen(rng);
        // ゼータ・メビウス
        std::vector<T> z = a, zs = a, e(N, T(0)), es(N, T(0));
        for (size_t s = 0; s < N; ++s)
            for (size_t t = 0; t < N; ++t) {
                if ((t & s) == t) e[s] = e[s] + a[t];
                if ((t & s) == s) es[s] = es[s] + a[t];
            }
        subset_zeta(z), superset_zeta(zs);
        SUBT_CHECK(z == e && zs == es);
        subset_mobius(z), superset_mobius(zs);
        SUBT_CHECK(z == a && zs == a);
        // 畳み込み
        std::vector<T> co(N, T(0)), ca(N, T(0)), cx(N, T(0)), cs(N, T(0));
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j) {
                co[i | j] = co[i | j] + a[i] * b[j], ca[i & j] = ca[i & j] + a[i] * b[j], cx[i ^ j] = cx[i ^ j] + a[i] * b[j];
                if ((i & j) == 0) cs[i | j] = cs[i | j] + a[i] * b[j];
            }
        SUBT_CHECK(or_convolution(a, b) == co && and_convolution(a, b) == ca && xor_convolution(a, b) == cx && subset_convolution(a, b) == cs);
    }
}

template <class T>
static std::vector<T> brute_exp(const std::vector<T>& f) {  // g[S] = Σ_{T ∋ S の最小の要素、T ⊆ S} f[T] g[S∖T]
    const size_t N = f.size();
    std::vector<T> g(N, T(0));
    g[0] = T(1);
    for (size_t s = 1; s < N; ++s) {
        const size_t low = s & (~s + 1), rest = s ^ low;
        for (size_t t = rest;; t = (t - 1) & rest) {  // 最小の要素を含む組 low | t
            g[s] = g[s] + f[low | t] * g[rest ^ t];
            if (t == 0) break;
        }
    }
    return g;
}
static void test_exp_log(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const size_t N = (size_t)1 << rnd(rng, 0, 9);
        std::vector<ModInt998> f(N);
        for (auto& x : f) x = (long long)rng();
        f[0] = 0;
        const std::vector<ModInt998> g = subset_exp(f);
        SUBT_CHECK(g == brute_exp(f) && subset_log(g) == f);
        std::vector<long long> fl(N);  // long long でも（小さい値）
        for (auto& x : fl) x = rnd(rng, -3, 3);
        fl[0] = 0;
        const std::vector<long long> gl = subset_exp(fl);
        SUBT_CHECK(gl == brute_exp(fl) && subset_log(gl) == fl);
    }
    // 全部 1 の exp はベル数、2^{C(|S|,2)}（グラフの数）の log は連結なグラフの数
    const int n = 6;
    std::vector<long long> one((size_t)1 << n, 1), graphs((size_t)1 << n);
    one[0] = 0;
    for (size_t s = 0; s < graphs.size(); ++s) {
        const int k = __builtin_popcountll(s);
        graphs[s] = 1LL << (k * (k - 1) / 2);
    }
    const std::vector<long long> bell = subset_exp(one), conn = subset_log(graphs);
    SUBT_CHECK(bell.back() == 203 && conn[1] == 1 && conn[3] == 1 && conn[7] == 4 && conn[15] == 38 && conn[31] == 728 && conn.back() == 26704);
    DynModInt::set_mod(1000000000);  // 合成数の mod でも（割り算を使わない）
    std::vector<DynModInt> d((size_t)1 << 7);
    for (auto& x : d) x = (long long)rng();
    d[0] = 0;
    SUBT_CHECK(subset_log(subset_exp(d)) == d);
    cases += 2;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_type<ModInt998>(rng, cases, 2000 * TEST_SCALE, [](std::mt19937_64& r) { return ModInt998((long long)r()); });
    test_type<long long>(rng, cases, 2000 * TEST_SCALE, [](std::mt19937_64& r) { return rnd(r, -1000, 1000); });
    DynModInt::set_mod(1000000000);  // 合成数・偶数の mod でも（XOR の逆変換は 2^n の逆元が要るので除く）
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {
        const size_t N = (size_t)1 << rnd(rng, 0, 6);
        std::vector<DynModInt> a(N), b(N), co(N), cs(N);
        for (auto& x : a) x = (long long)rng();
        for (auto& x : b) x = (long long)rng();
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j) {
                co[i | j] += a[i] * b[j];
                if ((i & j) == 0) cs[i | j] += a[i] * b[j];
            }
        SUBT_CHECK(or_convolution(a, b) == co && subset_convolution(a, b) == cs);
    }
    // 演算を渡す：部分集合の最大値
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {
        const size_t N = (size_t)1 << rnd(rng, 0, 8);
        std::vector<long long> a(N), e(N, LLONG_MIN);
        for (auto& x : a) x = rnd(rng, -100, 100);
        for (size_t s = 0; s < N; ++s)
            for (size_t t = 0; t < N; ++t)
                if ((t & s) == t) e[s] = std::max(e[s], a[t]);
        subset_zeta(a, [](long long x, long long y) { return std::max(x, y); });
        SUBT_CHECK(a == e);
    }
    // 長さがちがうものは 2 の冪にそろえる
    SUBT_CHECK(xor_convolution<long long>({1, 2, 3}, {4, 5}) == std::vector<long long>({4 + 10, 5 + 8, 12 + 0, 15 + 0}));
    test_exp_log(rng, cases);
    std::cout << "[subset_transform] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[subset_transform] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const size_t N = (size_t)1 << 20;
    std::vector<ModInt998> a(N), b(N);
    for (auto& x : a) x = (long long)rng();
    for (auto& x : b) x = (long long)rng();
    auto t0 = std::chrono::steady_clock::now();
    const auto c1 = xor_convolution(a, b);
    const auto c2 = or_convolution(a, b);
    const long long t1 = ms(t0);
    std::vector<ModInt998> sa(a.begin(), a.begin() + (1 << 18)), sb(b.begin(), b.begin() + (1 << 18));
    t0 = std::chrono::steady_clock::now();
    const auto c3 = subset_convolution(sa, sb);
    const long long t2 = ms(t0);
    std::vector<ModInt998> e(sa);
    e[0] = 0;
    t0 = std::chrono::steady_clock::now();
    const auto ex = subset_exp(e);
    const long long t3 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto lg = subset_log(ex);
    const long long t4 = ms(t0);
    std::cout << "[subset_transform] benchmark: xor_convolution + or_convolution（2^20）" << t1 << " ms, subset_convolution（2^18）" << t2 << " ms, subset_exp（2^18）"
              << t3 << " ms, subset_log（2^18）" << t4 << " ms  (check " << (c1[5] + c2[7] + c3[9]).val() % 7 << (lg == e ? "" : " log(exp) ≠ f") << ")\n";
}

}  // namespace subt_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace subt_demo {

{COMMON_DEMO}

static void demo() {
    banner("部分集合の変換（subset_transform） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  subset_zeta / subset_mobius / superset_zeta / superset_mobius（f を書きかえる。op を渡せる）\n"
              << "  walsh_hadamard / or_convolution / and_convolution / xor_convolution / subset_convolution\n";

    api("【1】 subset_zeta(f) / subset_mobius(f) / superset_zeta(f) / superset_mobius(f)", "なし（f を書きかえる）",
        "f[S] ← Σ_{T ⊆ S} f[T] / その逆 / ⊇ の版。op を渡すと max・min・gcd などでまとめる", "O(n 2^n)");
    RUN(std::vector<long long> f = {1, 2, 3, 4}; subset_zeta(f));
    SHOW(f, "f[3] = f[0] + f[1] + f[2] + f[3]");
    RUN(subset_mobius(f));
    SHOW(f, "元に戻る");
    RUN(std::vector<long long> g = {5, 1, 7, 2}; subset_zeta(g, [](long long x, long long y) { return std::max(x, y); }));
    SHOW(g, "部分集合の最大値");
    note("bit DP：「S の部分集合すべての和」を全 S で求めるときに。superset は「S を含む集合すべての和」");

    api("【2】 or_convolution(a, b) / and_convolution(a, b) / xor_convolution(a, b)", "vector<T>",
        "c[k] = Σ_{i ∘ j = k} a[i] b[j]（∘ は OR・AND・XOR）", "O(n 2^n)");
    SHOW(or_convolution<long long>({1, 1}, {1, 1}), "{0, 1} から 2 つ選んだ OR：0 が 1 通り、1 が 3 通り");
    SHOW(xor_convolution<long long>({0, 1, 1, 0}, {0, 1, 1, 0}), "{1, 2} から 2 つ選んだ XOR");
    {
        // ABC212 H：A の数を 1〜N 個並べて XOR ≠ 0 になる数。XOR = 0 の数をアダマール変換の上で等比数列の和に
        auto solve = [](long long N, const std::vector<int>& A) {
            using mint = ModInt998;
            std::vector<mint> h(1 << 16);
            for (int a : A) h[(size_t)a] = 1;
            walsh_hadamard(h);
            for (mint& x : h) x = (x == mint(1)) ? mint(N) : x * (x.pow(N) - 1) / (x - 1);  // x + x^2 + … + x^N
            walsh_hadamard(h, true);
            const mint K = (long long)A.size();
            const mint all = (K == mint(1)) ? mint(N) : K * (K.pow(N) - 1) / (K - 1);
            return all - h[0];
        };
        note("ABC212 H（Nim Counting）: アダマール変換して点ごとに x + x^2 + … + x^N、戻して XOR = 0 の数を引く");
        std::cout << "    例 1 → " << solve(2, {1, 2}) << "、例 2 → " << solve(100, {3, 5, 7}) << "\n";
    }

    api("【3】 subset_convolution(a, b)", "vector<T>", "c[S] = Σ_{T ⊆ S} a[T] b[S∖T]（集合 S を 2 つに分けるすべての方法）", "O(n^2 2^n)");
    SHOW(subset_convolution<long long>({1, 1, 1, 1}, {1, 1, 1, 1}), "{0, 1} の部分集合 S を 2 つに分ける方法の数 = 2^|S|");
    note("集合を k 個のグループに分ける数え上げ（bit DP の O(3^n) を O(n^2 2^n) に）で使う");

    api("【4】 subset_exp(f) / subset_log(f)", "vector<T>",
        "g[S] = Σ_{S の組分け} Π f[組]（f[∅] = 0）/ その逆（f[∅] = 1）。割り算なしなので long long・合成数の mod でも", "O(n^2 2^n)");
    SHOW(subset_exp<long long>({0, 1, 1, 1, 1, 1, 1, 1}), "どの組も 1 通り：組分けの数 = ベル数（|S| = 3 で 5）");
    {
        // 頂点の集合ごとの「辺の数 e(S)」から、連結な部分グラフ（全頂点を使う辺の選び方）の数：全体 2^{e(S)} の log
        const int n = 4;
        const std::vector<std::pair<int, int>> es = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}};
        std::vector<ModInt998> all((size_t)1 << n);
        for (size_t s = 0; s < all.size(); ++s) {
            int e = 0;
            for (auto [u, v] : es) e += (s >> u & 1) && (s >> v & 1);
            all[s] = ModInt998(2).pow(e);  // S の中の辺を自由に選ぶ
        }
        const std::vector<ModInt998> conn = subset_log(all);
        note("連結の数え上げ：「S の中の辺を自由に選ぶ」2^{e(S)} を subset_log すると「S が連結になる選び方」");
        std::cout << "    4 頂点（辺 01 12 23 30 02）で全頂点が連結になる辺の選び方 → " << conn.back() << "\n";
    }
}

}  // namespace subt_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    subt_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    subt_test::benchmark();    // 速度計測
    subt_demo::demo();         // 使い方デモ
    return 0;
}
