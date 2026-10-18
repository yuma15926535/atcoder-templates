// =============================================================================
//  ランダムテスト
// =============================================================================
namespace comb_test {

{COMMON_TEST}

// 全部の写像 [n] → [k] を並べて、玉・箱の区別に合わせて同じものをまとめて数える
static long long brute_twelve(int n, int k, bool bd, bool xd, int rule) {
    std::set<std::vector<int>> seen;
    std::vector<int> f((size_t)n, 0);
    long long total = 1;
    for (int i = 0; i < n; ++i) total *= k;
    if (k == 0) total = (n == 0 ? 1 : 0);
    for (long long code = 0; code < total; ++code) {
        long long c = code;
        for (int i = 0; i < n; ++i) f[(size_t)i] = (int)(c % k), c /= k;
        std::vector<int> cnt((size_t)k, 0);
        for (int x : f) ++cnt[(size_t)x];
        bool ok = true;
        for (int x : cnt)
            if ((rule == 1 && x > 1) || (rule == 2 && x == 0)) ok = false;
        if (!ok) continue;
        std::vector<int> key;
        if (bd && xd) key = f;
        else if (!bd && xd) key = cnt;
        else if (bd && !xd) {  // 箱を出てきた順に番号を振り直す
            std::vector<int> lab((size_t)k, -1);
            int nx = 0;
            for (int x : f) {
                if (lab[(size_t)x] < 0) lab[(size_t)x] = nx++;
                key.push_back(lab[(size_t)x]);
            }
        } else {
            key = cnt;
            std::sort(key.begin(), key.end());
        }
        seen.insert(key);
    }
    return (long long)seen.size();
}
static std::vector<std::vector<int>> all_perms(int n) {
    std::vector<int> p((size_t)n);
    std::iota(p.begin(), p.end(), 0);
    std::vector<std::vector<int>> res;
    do res.push_back(p);
    while (std::next_permutation(p.begin(), p.end()));
    return res;
}

template <class mint>
static void test_stirling(std::mt19937_64& rng, int& cases) {
    const int N = 70;
    const auto t2 = stirling2_table<mint>(N);
    const auto t1 = stirling1_table<mint>(N);
    for (int n = 0; n <= N; ++n) {
        const auto r2 = stirling2_row<mint>(n);
        const auto r1 = stirling1_row<mint>(n);
        COMB_CHECK(r2.size() == (size_t)n + 1 && r1.size() == (size_t)n + 1);
        for (int k = 0; k <= n; ++k, ++cases) {
            COMB_CHECK(r2[(size_t)k] == t2[(size_t)n][(size_t)k] && stirling2<mint>(n, k) == t2[(size_t)n][(size_t)k]);
            COMB_CHECK(r1[(size_t)k] == t1[(size_t)n][(size_t)k]);
        }
        COMB_CHECK(stirling2<mint>(n, n + 1) == mint(0) && stirling2<mint>(n, -1) == mint(0));
        mint b = 0;
        for (int k = 0; k <= n; ++k) b += t2[(size_t)n][(size_t)k];
        COMB_CHECK(bell<mint>(n) == b && bell_table<mint>(N)[(size_t)n] == b);
        // Σ_k c(n, k) = n!、x^n = Σ_k S(n, k) x(x−1)…(x−k+1)
        mint s1 = 0, f = 1;
        for (int k = 0; k <= n; ++k) s1 += t1[(size_t)n][(size_t)k];
        for (int i = 1; i <= n; ++i) f *= mint(i);
        COMB_CHECK(s1 == f);
        const long long x = rnd(rng, 0, 1000);
        mint lhs = mint(x).pow(n), rhs = 0, ff = 1;
        for (int k = 0; k <= n; ++k) rhs += t2[(size_t)n][(size_t)k] * ff, ff *= mint(x - k);
        COMB_CHECK(lhs == rhs);
    }
    // 大きい n：1 つだけ版と行の版
    for (int it = 0; it < 4 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1000, 3000);
        const auto r2 = stirling2_row<mint>(n);
        for (int t = 0; t < 20; ++t) {
            const int k = (int)rnd(rng, 0, n);
            COMB_CHECK(stirling2<mint>(n, k) == r2[(size_t)k]);
        }
        const auto r1 = stirling1_row<mint>(n);  // c(n, k) = (n−1) c(n−1, k) + c(n−1, k−1) を 1 段だけ確かめる
        const auto q1 = stirling1_row<mint>(n - 1);
        for (int k = 1; k < n; ++k) COMB_CHECK(r1[(size_t)k] == mint(n - 1) * q1[(size_t)k] + q1[(size_t)k - 1]);
    }
}

static void test_small_exact(int& cases) {
    // 順列を全部並べて：サイクルの数・完全順列
    const auto t1 = stirling1_table<long long>(8);
    const auto der = derangement_table<long long>(8);
    for (int n = 0; n <= 8; ++n) {
        std::vector<long long> cyc((size_t)n + 1, 0);
        long long d = 0;
        for (const auto& p : all_perms(n)) {
            std::vector<bool> vis((size_t)n, false);
            int c = 0;
            bool fixed = false;
            for (int i = 0; i < n; ++i) {
                if (p[(size_t)i] == i) fixed = true;
                if (vis[(size_t)i]) continue;
                ++c;
                for (int j = i; !vis[(size_t)j]; j = p[(size_t)j]) vis[(size_t)j] = true;
            }
            ++cyc[(size_t)c];
            if (!fixed) ++d;
        }
        for (int k = 0; k <= n; ++k, ++cases) COMB_CHECK(t1[(size_t)n][(size_t)k] == cyc[(size_t)k]);
        COMB_CHECK(der[(size_t)n] == d);
    }
    // 集合の分割を全部並べて（制限成長列）：S(n, k)・ベル数
    const auto t2 = stirling2_table<long long>(10);
    const auto bt = bell_table<long long>(10);
    for (int n = 0; n <= 10; ++n) {
        std::vector<long long> cnt((size_t)n + 1, 0);
        std::vector<int> a((size_t)n, 0);
        std::function<void(int, int)> rec = [&](int i, int mx) {
            if (i == n) { ++cnt[(size_t)mx]; return; }
            for (int v = 0; v <= mx; ++v) a[(size_t)i] = v, rec(i + 1, std::max(mx, v + 1));
        };
        rec(0, 0);
        long long b = 0;
        for (int k = 0; k <= n; ++k, ++cases) {
            COMB_CHECK(t2[(size_t)n][(size_t)k] == cnt[(size_t)k]);
            b += cnt[(size_t)k];
        }
        COMB_CHECK(bt[(size_t)n] == b);
    }
    COMB_CHECK(bell_table<long long>(25)[25] == 4638590332229999353LL);
    // 分割数：ナップサックの DP と
    const int N = 400;
    const auto p = partition_table<long long>(N);
    const auto pc = partition_count_table<long long>(N, 30);
    for (int k = 0; k <= 30; ++k) {  // k 以下の数を何個でも使うナップサック
        std::vector<long long> q((size_t)N + 1, 0);
        q[0] = 1;
        for (int j = 1; j <= k; ++j)
            for (int v = j; v <= N; ++v) q[(size_t)v] += q[(size_t)(v - j)];
        for (int v = 0; v <= N; ++v, ++cases) COMB_CHECK(pc[(size_t)k][(size_t)v] == q[(size_t)v]);
    }
    {
        std::vector<long long> q((size_t)N + 1, 0);
        q[0] = 1;
        for (int j = 1; j <= N; ++j)
            for (int v = j; v <= N; ++v) q[(size_t)v] += q[(size_t)(v - j)];
        for (int v = 0; v <= N; ++v, ++cases) COMB_CHECK(p[(size_t)v] == q[(size_t)v]);
    }
    COMB_CHECK(partition_table<long long>(100)[100] == 190569292LL);
    COMB_CHECK(derangement_table<long long>(0).size() == 1 && derangement_table<long long>(0)[0] == 1);
}

template <class mint>
static void test_twelve(std::mt19937_64& rng, int& cases) {
    for (int n = 0; n <= 6; ++n)
        for (int k = 0; k <= 6; ++k)
            for (int m = 0; m < 4; ++m)
                for (int rule = 0; rule < 3; ++rule, ++cases) {
                    const bool bd = m & 1, xd = m & 2;
                    COMB_CHECK(balls_into_boxes<mint>(n, k, bd, xd, rule) == mint(brute_twelve(n, k, bd, xd, rule)));
                }
    // 大きい n：Σ_{i ≤ k} S(n, i) と、全射の数（包除原理）と
    for (int it = 0; it < 200 * TEST_SCALE; ++it, ++cases) {
        const long long n = it % 2 ? rnd(rng, 0, 1000000000000000000LL) : rnd(rng, 0, 60), k = rnd(rng, 0, 40);
        mint s = 0;
        for (long long i = 0; i <= k; ++i) s += stirling2<mint>(n, i);
        COMB_CHECK(balls_into_boxes<mint>(n, k, true, false, 0) == s);
        const mint sur = inclusion_exclusion_count<mint>((int)k, [&](int m) { return mint(m).pow(n); });
        COMB_CHECK(surjections<mint>(n, k) == sur && balls_into_boxes<mint>(n, k, true, true, 2) == sur);
    }
    // 玉も箱も区別しない：分割数の表と
    const auto pc = partition_count_table<mint>(300, 300);
    for (int it = 0; it < 300; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 300), k = (int)rnd(rng, 0, 300);
        COMB_CHECK(balls_into_boxes<mint>(n, k, false, false, 0) == pc[(size_t)k][(size_t)n]);
        COMB_CHECK(balls_into_boxes<mint>(n, k, false, false, 2) == (n >= k ? pc[(size_t)k][(size_t)(n - k)] : mint(0)));
    }
}

static void test_inclusion(std::mt19937_64& rng, int& cases) {
    // [1, X] で、どの d_i でも割り切れない数（inclusion_exclusion）
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 8);
        const long long X = rnd(rng, 0, 3000);
        std::vector<long long> d((size_t)n);
        for (auto& x : d) x = rnd(rng, 1, 40);
        const long long got = inclusion_exclusion<long long>(n, [&](long long S) {
            long long l = 1;
            for (int i = 0; i < n; ++i)
                if (S >> i & 1) l = lcm_cap(l, d[(size_t)i], 1000000);
            return X / l;
        });
        long long e = 0;
        for (long long v = 1; v <= X; ++v) {
            bool any = false;
            for (long long x : d) any = any || v % x == 0;
            e += !any;
        }
        COMB_CHECK(got == e);
        // count_divisible_by_any（負の範囲も）
        const long long l = rnd(rng, -2000, 2000), r = rnd(rng, -2000, 2000);
        long long e2 = 0;
        for (long long v = l; v < r; ++v) {
            bool any = false;
            for (long long x : d) any = any || v % x == 0;
            e2 += any;
        }
        COMB_CHECK(count_divisible_by_any(l, r, d) == e2);
        COMB_CHECK(count_divisible_by_any_closed(l, r, d) == e2 + (r >= l && [&] { for (long long x : d) if (r % x == 0) return true; return false; }() ? 1 : 0));
    }
    // 大きい範囲・大きい d：Int128 で 2 つ・3 つの包除を直接
    using math_internal::ceil_div128;
    auto cm = [&](long long l, long long r, Int128 d) { return ceil_div128(r, d) - ceil_div128(l, d); };
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        long long l = rnd_ll(rng), r = rnd_ll(rng);
        if (l > r) std::swap(l, r);
        const long long a = it % 3 ? rnd(rng, 1, 1000000000000LL) : std::max(1LL, std::llabs(rnd_ll(rng) / 2)), b = rnd(rng, 2, 1000000), c = rnd(rng, 2, 1000);
        auto lc = [](Int128 x, Int128 y) {
            Int128 g = x, h = y;
            while (h) { const Int128 t = g % h; g = h, h = t; }
            return x / g * y;
        };
        const Int128 e = cm(l, r, a) + cm(l, r, b) + cm(l, r, c) - cm(l, r, lc(a, b)) - cm(l, r, lc(a, c)) - cm(l, r, lc(b, c)) + cm(l, r, lc(lc(a, b), c));
        if (e > (Int128)LLONG_MAX) continue;
        COMB_CHECK((Int128)count_divisible_by_any(l, r, {a, b, c}) == e);
    }
    COMB_CHECK((Int128)count_divisible_by_any(LLONG_MIN, LLONG_MAX, {4, 6}) == cm(LLONG_MIN, LLONG_MAX, 4) + cm(LLONG_MIN, LLONG_MAX, 6) - cm(LLONG_MIN, LLONG_MAX, 12));
    COMB_CHECK(count_divisible_by_any_closed(LLONG_MAX, LLONG_MAX, {7, 73}) == 1 && count_divisible_by_any_closed(LLONG_MAX, LLONG_MAX, {2}) == 0);
    COMB_CHECK(count_divisible_by_any(0, 1, {1000000007, 998244353}) == 1 && count_divisible_by_any(1, 1, {1}) == 0 && count_divisible_by_any(5, 9, {}) == 0);
    // 素数を 40 個：lcm がすぐ範囲を超えるので速い
    {
        std::vector<long long> ps;
        for (long long v = 2; ps.size() < 40; ++v) {
            bool pr = true;
            for (long long q = 2; q * q <= v; ++q) pr = pr && v % q != 0;
            if (pr) ps.push_back(v);
        }
        const long long X = 100000;
        long long e = 0;
        for (long long v = 1; v <= X; ++v) {
            bool any = false;
            for (long long p : ps) any = any || v % p == 0;
            e += any;
        }
        COMB_CHECK(count_divisible_by_any_closed(1, X, ps) == e);
        COMB_CHECK(count_divisible_by_any_closed(1, X, {2, 4, 6, 3, 9, 2}) == X / 2 + X / 3 - X / 6);
        cases += 2;
    }
}

template <class mint>
static void test_transforms(std::mt19937_64& rng, int& cases) {
    const Binomial<mint> B(500);
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, it % 4 ? 30 : 300);
        std::vector<mint> g((size_t)n);
        for (auto& x : g) x = mint(rnd(rng, 0, 1000000000));
        std::vector<mint> f((size_t)n, 0), gi((size_t)n, 0), A((size_t)n, 0), E((size_t)n, 0);
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b) {
                if (b <= a) f[(size_t)a] += B.C(a, b) * g[(size_t)b], gi[(size_t)a] += ((a - b) % 2 ? -B.C(a, b) : B.C(a, b)) * g[(size_t)b];
                if (b >= a) A[(size_t)a] += B.C(b, a) * g[(size_t)b], E[(size_t)a] += ((b - a) % 2 ? -B.C(b, a) : B.C(b, a)) * g[(size_t)b];
            }
        COMB_CHECK(binomial_transform(g) == f && binomial_inverse(g) == gi && at_least_from_exactly(g) == A && exactly_from_at_least(g) == E);
        COMB_CHECK(binomial_inverse(binomial_transform(g)) == g && exactly_from_at_least(at_least_from_exactly(g)) == g);
    }
    // q 二項係数：q-パスカルの漸化式 [n, k] = [n−1, k−1] + q^k [n−1, k] と比べる（位数の小さい q も）
    std::vector<long long> qs = {0, 1, 2, 3, -1, (long long)mint::mod() - 1, (long long)mint::mod() + 1, 1000000};
    if (mint::mod() == 998244353)
        for (long long d : {2, 4, 7, 8, 14, 17}) qs.push_back((long long)mint(3).pow((998244353 - 1) / d).val());  // 位数 d
    if (mint::mod() == 1000000007)
        for (long long d : {2, 3, 6, 41}) qs.push_back((long long)mint(5).pow((1000000007 - 1) / d).val());
    for (long long q : qs) {
        const int N = 60;
        std::vector<std::vector<mint>> t(N + 1, std::vector<mint>(N + 1, 0));
        const mint Q(q);
        for (int n = 0; n <= N; ++n) {
            t[(size_t)n][0] = 1;
            for (int k = 1; k <= n; ++k) t[(size_t)n][(size_t)k] = t[(size_t)n - 1][(size_t)k - 1] + Q.pow(k) * t[(size_t)n - 1][(size_t)k];
        }
        for (int n = 0; n <= N; ++n)
            for (int k = -1; k <= n + 1; ++k, ++cases) COMB_CHECK(q_binomial<mint>(n, k, q) == (k < 0 || k > n ? mint(0) : t[(size_t)n][(size_t)k]));
    }
}
static void test_q_binomial_small_mod(int& cases) {  // mod が小さい素数（q の位数 d で ⌊n/d⌋ ≥ mod になる）でも
    for (int p : {2, 3, 5, 7, 13}) {
        DynModInt::set_mod(p);
        const int N = 45;
        for (long long q = -2; q <= p + 1; ++q) {
            const DynModInt Q(q);
            std::vector<std::vector<DynModInt>> t(N + 1, std::vector<DynModInt>(N + 1, DynModInt(0)));
            for (int n = 0; n <= N; ++n) {
                t[(size_t)n][0] = 1;
                for (int k = 1; k <= n; ++k) t[(size_t)n][(size_t)k] = t[(size_t)n - 1][(size_t)k - 1] + Q.pow(k) * t[(size_t)n - 1][(size_t)k];
            }
            for (int n = 0; n <= N; ++n)
                for (int k = 0; k <= n; ++k, ++cases) COMB_CHECK(q_binomial<DynModInt>(n, k, q) == t[(size_t)n][(size_t)k]);
        }
    }
    DynModInt::set_mod(998244353);
}

static long long brute_syt(std::vector<int> sh, std::map<std::vector<int>, long long>& memo) {  // 角を 1 つずつ取る
    while (!sh.empty() && sh.back() == 0) sh.pop_back();
    if (sh.empty()) return 1;
    auto it = memo.find(sh);
    if (it != memo.end()) return it->second;
    long long r = 0;
    for (size_t i = 0; i < sh.size(); ++i)
        if (i + 1 == sh.size() || sh[i + 1] < sh[i]) {
            --sh[i];
            r += brute_syt(sh, memo);
            ++sh[i];
        }
    return memo[sh] = r;
}

static void test_enumerate(std::mt19937_64& rng, int& cases) {
    const std::vector<long long> bell_ll = bell_table<long long>(12), part_ll = partition_table<long long>(40);
    for (int n = 0; n <= 10; ++n) {
        long long cnt = 0;
        std::set<std::vector<int>> seen;
        for_each_set_partition(n, [&](const std::vector<int>& g, int k) {
            ++cnt;
            int mx = -1;
            for (int x : g) {  // 初めて出る順の番号
                COMB_CHECK(x <= mx + 1);
                mx = std::max(mx, x);
            }
            COMB_CHECK(mx + 1 == k && (int)g.size() == n);
            if (n <= 7) seen.insert(g);
        });
        COMB_CHECK(cnt == bell_ll[(size_t)n] && (n > 7 || (long long)seen.size() == cnt));
        ++cases;
    }
    const auto pct = partition_count_table<long long>(40, 40);  // pct[k][n]：n を k 個以下に分ける数 = 最大の部分が k 以下
    for (int n = 0; n <= 40; ++n) {
        long long cnt = 0;
        for_each_integer_partition(n, [&](const std::vector<int>& p) {
            ++cnt;
            int s = 0;
            for (size_t i = 0; i < p.size(); ++i) {
                s += p[i];
                COMB_CHECK(p[i] >= 1 && (i == 0 || p[i] <= p[i - 1]));
            }
            COMB_CHECK(s == n);
        });
        COMB_CHECK(cnt == part_ll[(size_t)n]);
        const int mp = (int)rnd(rng, 0, n + 1);
        long long c2 = 0;
        std::vector<int> prev;
        for_each_integer_partition(
            n, [&](const std::vector<int>& p) {
                ++c2;
                COMB_CHECK(p.empty() || p[0] <= mp);
                COMB_CHECK(prev.empty() || p < prev);  // 辞書順で大きい方から
                prev = p;
            },
            mp);
        COMB_CHECK(c2 == (mp == 0 ? (n == 0) : pct[(size_t)std::min(mp, n)][(size_t)n]));
        cases += 2;
    }
    // フック長公式
    std::map<std::vector<int>, long long> memo;
    for (int it = 0; it < 500; ++it, ++cases) {
        std::vector<int> sh;
        int rem = (int)rnd(rng, 0, 14), mx = rem;
        while (rem > 0) {
            const int x = (int)rnd(rng, 1, std::min(rem, mx));
            sh.push_back(x), rem -= x, mx = x;
        }
        COMB_CHECK(young_tableaux_count<ModInt998>(sh) == ModInt998(brute_syt(sh, memo)));
    }
    COMB_CHECK(young_tableaux_count<ModInt998>({3, 2}) == ModInt998(5) && young_tableaux_count<ModInt998>(std::vector<int>{}) == ModInt998(1));
    COMB_CHECK(young_tableaux_count<ModInt998>({5, 5}) == ModInt998(42));  // 2 行 n 列はカタラン数
    {  // 部分の数が多くても再帰しない
        long long calls = 0;
        size_t len = 0;
        for_each_integer_partition(1000000, [&](const std::vector<int>& p) { ++calls, len = p.size(); }, 1);
        COMB_CHECK(calls == 1 && len == 1000000);
        long long c0 = 0;
        for_each_integer_partition(5, [&](const std::vector<int>&) { ++c0; }, 0);
        for_each_integer_partition(0, [&](const std::vector<int>& p) { c0 += 10 + (long long)p.size(); });
        COMB_CHECK(c0 == 10);
        cases += 2;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_transforms<ModInt998>(rng, cases);
    test_transforms<ModInt107>(rng, cases);
    test_enumerate(rng, cases);
    test_q_binomial_small_mod(cases);
    test_stirling<ModInt998>(rng, cases);
    test_stirling<ModInt107>(rng, cases);  // NTT が使えない素数
    test_small_exact(cases);
    test_twelve<ModInt998>(rng, cases);
    test_twelve<ModInt107>(rng, cases);
    test_inclusion(rng, cases);
    {  // partition_row（級数の逆数）と五角数定理の表
        const auto a = partition_row<ModInt998>(5000), b = partition_table<ModInt998>(5000);
        const auto c = partition_row<ModInt107>(3000), d = partition_table<ModInt107>(3000);
        COMB_CHECK(a == b && c == d && partition_row<ModInt998>(0).size() == 1);
        cases += 3;
    }
    {  // 小さい素数の mod（DynModInt）：mod 7 でも n < 7 なら使える
        DynModInt::set_mod(7);
        COMB_CHECK(stirling2<DynModInt>(3, 2) == DynModInt(3) && surjections<DynModInt>(3, 2) == DynModInt(6) && bell<DynModInt>(4) == DynModInt(15 % 7));
        COMB_CHECK(balls_into_boxes<DynModInt>(3, 2, true, false, 0) == DynModInt(4) && surjections<DynModInt>(3, 2000000000) == DynModInt(0));
        DynModInt::set_mod(998244353);
        // k が巨大：階乗の表を作らずに
        using mint = ModInt998;
        COMB_CHECK(balls_into_boxes<mint>(2, 1000000000, true, true, 1) == mint(1000000000) * mint(999999999));
        COMB_CHECK(balls_into_boxes<mint>(3, 1000000000, false, true, 1) == mint(1000000000) * mint(999999999) * mint(999999998) / mint(6));
        COMB_CHECK(balls_into_boxes<mint>(2, 1000000000, false, true, 0) == mint(1000000001) * mint(1000000000) / mint(2));
        COMB_CHECK(balls_into_boxes<mint>(100000, 100000, false, false, 0) == partition_row<mint>(100000)[100000]);
        cases += 4;
    }
    std::cout << "[combinatorics] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[combinatorics] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    using mint = ModInt998;
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    const auto r2 = stirling2_row<mint>(500000);
    const long long t_s2 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto r1 = stirling1_row<mint>(200000);
    const long long t_s1 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto p = partition_row<mint>(1000000);
    const long long t_p = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto t2 = stirling2_table<mint>(3000);
    const long long t_tab = ms(t0);
    t0 = std::chrono::steady_clock::now();
    std::vector<long long> ds;
    for (long long v = 2; ds.size() < 20; ++v) {
        bool pr = true;
        for (long long q = 2; q * q <= v; ++q) pr = pr && v % q != 0;
        if (pr) ds.push_back(v);
    }
    const long long c = count_divisible_by_any(0, 1000000000000000000LL, ds);
    const long long t_ie = ms(t0);
    std::cout << "[combinatorics] benchmark: stirling2_row(5×10^5) " << t_s2 << " ms, stirling1_row(2×10^5) " << t_s1 << " ms, partition_row(10^6) " << t_p
              << " ms, stirling2_table(3000) " << t_tab << " ms, count_divisible_by_any（10^18, 素数 20 個）" << t_ie << " ms  (check "
              << (r2[1234] + r1[777] + p.back() + t2[3000][1500] + mint(c)).val() % 7 << ")\n";
}

}  // namespace comb_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace comb_demo {

{COMMON_DEMO}

static void demo() {
    using mint = ModInt998;
    banner("数え上げ（combinatorics） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  stirling2 / stirling2_row / stirling2_table / stirling1_row / stirling1_table / bell / bell_table\n"
              << "  partition_table / partition_count_table / derangement_table / surjections / balls_into_boxes\n"
              << "  inclusion_exclusion / inclusion_exclusion_count / count_divisible_by_any(_closed)\n"
              << "  （この中では using mint = ModInt998; としている）\n";

    api("【1】 stirling2<mint>(n, k) / stirling2_row<mint>(n) / stirling2_table<T>(N)", "mint / vector<mint> / vector<vector<T>>",
        "第 2 種スターリング数 S(n, k)：n 個の区別できるものを k 個の空でない（区別しない）組に分ける数", "O(k log n) / O(n log n) / O(N^2)");
    SHOW(stirling2<mint>(5, 2), "{1..5} を 2 組に：2^4 − 1 = 15");
    SHOW(stirling2_row<mint>(5), "S(5, 0..5)");
    SHOW(stirling2<mint>(1000000000000000000LL, 3), "n が 10^18 でも k が小さければ");
    SHOW(stirling2_table<long long>(4)[4], "long long の表も（mod なし）");

    api("【2】 stirling1_row<mint>(n) / stirling1_table<T>(N)", "vector<mint> / vector<vector<T>>",
        "第 1 種スターリング数（符号なし）c(n, k)：n 個の順列でサイクルがちょうど k 個のものの数。x(x+1)…(x+n−1) の係数", "O(n log^2 n) / O(N^2)");
    SHOW(stirling1_row<mint>(4), "c(4, 0..4)。合計 24 = 4!");
    note("符号つき s(n, k) は (−1)^{n−k} c(n, k)");

    api("【3】 bell<mint>(n) / bell_table<T>(N)", "mint / vector<T>", "ベル数 B_n：n 個を空でない組にいくつかに分ける数（= Σ_k S(n, k)）", "O(n log n) / O(N^2)");
    SHOW(bell_table<long long>(8), "B_0..B_8");
    SHOW(bell<mint>(100000), "B_100000 mod 998244353");

    api("【4】 partition_row<mint>(N) / partition_table<T>(N) / partition_count_table<T>(N, K)", "vector<mint> / vector<T> / vector<vector<T>>",
        "分割数 p(n)：n を正の整数の和に分ける数（順番は区別しない）/ 同じく long long でも使える版 / t[k][n]：k 個以下に分ける数", "O(N log N) / O(N √N) / O(N K)");
    SHOW(partition_row<mint>(100000)[100000], "p(10^5) mod 998244353");
    SHOW(partition_table<long long>(10), "p(0..10)：5 = 5 = 4+1 = 3+2 = 3+1+1 = 2+2+1 = 2+1+1+1 = 1+1+1+1+1 で 7");
    SHOW(partition_count_table<long long>(10, 3)[3][10], "10 を 3 個以下に分ける");

    api("【5】 derangement_table<T>(N) / surjections<mint>(n, k)", "vector<T> / mint",
        "完全順列（どの i も i 番目に来ない順列）の数 / n 個から k 個への全射（どの行き先にも 1 個以上）の数", "O(N) / O(k log n)");
    SHOW(derangement_table<long long>(6), "D_0..D_6（プレゼント交換で全員が他人のものをもらう）");
    SHOW(surjections<mint>(5, 3), "3! S(5, 3) = 6 × 25");

    api("【6】 balls_into_boxes<mint>(n, k, 玉を区別するか, 箱を区別するか, rule)", "mint",
        "写像 12 相：n 個の玉を k 個の箱に入れる方法の数。rule：0 = 何個でもよい、1 = 各箱 1 個以下、2 = 各箱 1 個以上", "O(k log n) など（玉も箱も区別しないときは分割数）");
    {
        const char* names[4] = {"玉 区別しない・箱 区別しない", "玉 区別する  ・箱 区別しない", "玉 区別しない・箱 区別する  ", "玉 区別する  ・箱 区別する  "};
        note("n = 5、k = 3 の 12 通り：                  何個でも   1 個以下   1 個以上");
        for (int m = 0; m < 4; ++m) {
            std::cout << "    " << names[m];
            for (int rule = 0; rule < 3; ++rule) std::cout << std::setw(11) << balls_into_boxes<mint>(5, 3, m & 1, m & 2, rule).val();
            std::cout << "\n";
        }
        note("（3^5 = 243、S(5,1)+S(5,2)+S(5,3) = 41、p(5 を 3 個以下) = 5、7H5 = C(7, 5) = 21、全射 150、C(4, 2) = 6 …）");
    }
    SHOW(balls_into_boxes<mint>(2, 3, true, true, 1), "「1 個以下」は n ≤ k のときだけ 0 でない：3P2 = 6");

    api("【7】 inclusion_exclusion<T>(n, f)", "T", "Σ_{S ⊆ {0..n−1}} (−1)^{|S|} f(S)：f(S) =「条件 S を全部破る」数 → どの条件も破らない数", "O(2^n × f)");
    {
        // ABC297 F（Minimum Bounding Box 2）：h × w の長方形がちょうど外接長方形になる選び方を「4 辺のどれかが空」の包除で
        auto solve = [](long long H, long long W, long long K) {
            Binomial<mint> B(H * W);
            mint s = 0;
            for (long long h = 1; h <= H; ++h)
                for (long long w = 1; w <= W; ++w) {
                    const mint ways = inclusion_exclusion<mint>(4, [&](long long S) {  // ビット 0,1 = 上・下の行、2,3 = 左・右の列が空
                        const long long hh = std::max(0LL, h - (S & 1) - (S >> 1 & 1)), ww = std::max(0LL, w - (S >> 2 & 1) - (S >> 3 & 1));
                        return B.C(hh * ww, K);
                    });
                    s += ways * mint(h * w) * mint((H - h + 1) * (W - w + 1));
                }
            return s / B.C(H * W, K);
        };
        note("ABC297 F（Minimum Bounding Box 2）: 外接長方形がちょうど h × w になる選び方を「上下左右の 4 辺のどれかが空」の包除で");
        std::cout << "    (2, 2, 2) → " << solve(2, 2, 2) << "、(10, 10, 1) → " << solve(10, 10, 1) << "、(314, 159, 2653) → " << solve(314, 159, 2653) << "\n";
    }

    api("【8】 inclusion_exclusion_count<mint>(n, g)", "mint", "Σ_i (−1)^i C(n, i) g(n − i)：g(m) =「n 個のうち決まった m 個だけを使ってよい」数 → n 個を全部使う数", "O(n × g)");
    SHOW(inclusion_exclusion_count<mint>(3, [](int m) { return mint(m).pow(5); }), "5 個を 3 個へ、3 個とも使う（全射 150）");
    {
        // ABC242 F（Black and White Rooks）：r 行 c 列を全部使って k 個を置く数 f(r, c, k) を行・列の 2 重の包除で
        auto solve = [](int N, int M, int Bk, int Wk) {
            Binomial<mint> B(N * M + N + M);
            auto f = [&](int k) {
                std::vector<std::vector<mint>> t((size_t)N + 1, std::vector<mint>((size_t)M + 1));
                for (int r = 0; r <= N; ++r)
                    for (int c = 0; c <= M; ++c)
                        t[(size_t)r][(size_t)c] = inclusion_exclusion_count<mint>(r, [&](int rr) { return inclusion_exclusion_count<mint>(c, [&](int cc) { return B.C(rr * cc, k); }); });
                return t;
            };
            const auto fb = f(Bk), fw = f(Wk);
            mint s = 0;
            for (int r1 = 0; r1 <= N; ++r1)
                for (int c1 = 0; c1 <= M; ++c1) {
                    if (fb[(size_t)r1][(size_t)c1] == mint(0)) continue;
                    for (int r2 = 0; r1 + r2 <= N; ++r2)
                        for (int c2 = 0; c1 + c2 <= M; ++c2)
                            s += B.C(N, r1) * B.C(N - r1, r2) * B.C(M, c1) * B.C(M - c1, c2) * fb[(size_t)r1][(size_t)c1] * fw[(size_t)r2][(size_t)c2];
                }
            return s;
        };
        note("ABC242 F（Black and White Rooks）: 「r 行 c 列を全部使う」置き方を inclusion_exclusion_count の 2 重で");
        std::cout << "    (2, 2, 1, 1) → " << solve(2, 2, 1, 1) << "、(1, 2, 1, 1) → " << solve(1, 2, 1, 1) << "、(40, 40, 30, 30) → " << solve(40, 40, 30, 30) << "\n";
    }

    api("【9】 count_divisible_by_any(l, r, ds) / _closed(L, R, ds)", "long long", "[l, r)（[L, R]）の整数で ds のどれかで割り切れるものの数（包除原理。lcm が大きすぎる組は打ち切る）", "O(2^|ds|) 以下");
    SHOW(count_divisible_by_any_closed(1, 100, {2, 3, 5}), "1〜100 で 2・3・5 のどれかの倍数");
    SHOW(count_divisible_by_any(-10, 10, {4, 6}), "負の範囲も（0 も数える）");
    {
        // ABC131 C（Anti-Division）：A〜B で C でも D でも割り切れない数
        auto solve = [](long long A, long long B, long long C, long long D) { return (B - A + 1) - count_divisible_by_any_closed(A, B, {C, D}); };
        note("ABC131 C（Anti-Division）: (B − A + 1) − count_divisible_by_any_closed(A, B, {C, D})");
        std::cout << "    (4, 9, 2, 3) → " << solve(4, 9, 2, 3) << "、(10, 40, 6, 8) → " << solve(10, 40, 6, 8) << "、例 3 → "
                  << solve(314159265358979323LL, 846264338327950288LL, 419716939, 937510582) << "\n";
    }

    api("【10】 binomial_transform / binomial_inverse / at_least_from_exactly / exactly_from_at_least（<mint> の vector）", "vector<mint>",
        "f(n) = Σ_k C(n, k) g(k) とその逆 / A[k] = Σ_{j≥k} C(j, k) E[j] とその逆（「k 個を選んで満たす」の和 A から「ちょうど k 個」E へ）", "O(n log n)");
    note("「ちょうど k 個の条件を満たす数」は直接数えにくい → 「k 個の条件を選び、それを全部満たす」組の数 A[k] を数えて exactly_from_at_least");
    {
        // 例：1〜n の順列で、ちょうど k 個の i が p_i = i（不動点）になる数。A[k] = C(n, k) (n − k)!（k 個を選んで固定、残りは自由）
        const int n = 5;
        std::vector<mint> A(n + 1);
        const Binomial<mint> B(n);
        for (int k = 0; k <= n; ++k) A[(size_t)k] = B.C(n, k) * B.fact(n - k);
        std::cout << "    n = 5 の順列の不動点の数ごとの個数 → ";
        for (const mint& x : exactly_from_at_least(A)) std::cout << x << " ";
        std::cout << "（0 個は完全順列 44）\n";
    }
    SHOW(binomial_inverse<ModInt998>({1, 2, 4, 8, 16}), "2^n = Σ C(n, k)·1 なので 1 に戻る");

    api("【11】 for_each_set_partition(n, f) / for_each_integer_partition(n, f, max_part)", "なし（f を呼ぶ）",
        "{0..n−1} の組分けを全部 f(g, k)（g[i] = i の組の番号、k = 組の数）/ n の分割を全部 f(parts)（大きい順）", "O(B_n · n) / O(p(n) · n)");
    {
        long long c = 0;
        for_each_set_partition(4, [&](const std::vector<int>&, int k) { c += k == 2; });
        std::cout << "  for_each_set_partition(4, …)：4 個を 2 組に分ける方法 → " << c << "（S(4, 2) = 7）\n";
        // ABC390 D：袋の石をいくつかの組にまとめ、組ごとの和の XOR として出てくる値の種類（N ≤ 12）
        const std::vector<long long> a = {2, 5, 7};
        std::set<long long> vals;
        std::vector<long long> sum(a.size());
        for_each_set_partition((int)a.size(), [&](const std::vector<int>& g, int k) {
            std::fill(sum.begin(), sum.begin() + k, 0);
            for (size_t i = 0; i < a.size(); ++i) sum[(size_t)g[i]] += a[i];
            long long x = 0;
            for (int j = 0; j < k; ++j) x ^= sum[(size_t)j];
            vals.insert(x);
        });
        note("ABC390 D（Stone XOR）: N ≤ 12 の組分けを全部（B_12 = 420 万通り）たどって XOR の種類を数える");
        std::cout << "    {2, 5, 7} → " << vals.size() << "\n";
        std::cout << "  for_each_integer_partition(5, …) → ";
        for_each_integer_partition(5, [&](const std::vector<int>& p) {
            std::cout << "{";
            for (size_t i = 0; i < p.size(); ++i) std::cout << p[i] << (i + 1 < p.size() ? "," : "");
            std::cout << "} ";
        });
        std::cout << " // 大きい順、p(5) = 7 通り\n";
    }

    api("【12】 q_binomial<mint>(n, k, q) / young_tableaux_count<mint>(shape)", "mint",
        "q 二項係数 [n, k]_q（F_q^n の k 次元部分空間の数。q^d ≡ 1 でも q-Lucas で正しく）/ 標準ヤング盤の数（フック長公式）", "O(n)");
    SHOW(q_binomial<ModInt998>(4, 2, 2), "F_2^4 の 2 次元部分空間の数");
    SHOW(q_binomial<ModInt998>(10, 3, 1), "q = 1 なら普通の C(10, 3)");
    SHOW(young_tableaux_count<ModInt998>({3, 2, 1}), "形 (3, 2, 1) に 1〜6 を行・列とも増えるように入れる数");
    SHOW(young_tableaux_count<ModInt998>({10, 10}), "2 行 10 列 = カタラン数 C_10");
}

}  // namespace comb_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    comb_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    comb_test::benchmark();    // 速度計測
    comb_demo::demo();         // 使い方デモ
    return 0;
}
