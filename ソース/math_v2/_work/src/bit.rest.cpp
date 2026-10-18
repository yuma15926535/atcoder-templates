// =============================================================================
//  ランダムテスト
// =============================================================================
namespace bit_test {

{COMMON_TEST}

static const long long LMAX = std::numeric_limits<long long>::max();

static void test_single(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300000 * TEST_SCALE; ++it, ++cases) {
        const long long x = rnd_ll(rng);
        const unsigned long long u = (unsigned long long)x;
        int cnt = 0, top = -1, low = -1;
        std::vector<int> idx;
        for (int i = 0; i < 64; ++i)
            if ((u >> i) & 1) {
                cnt++, top = i, idx.push_back(i);
                if (low < 0) low = i;
            }
        BIT_CHECK(bit_count(x) == cnt && top_bit(x) == top && low_bit(x) == low && bit_length(x) == top + 1 && bit_indices(x) == idx);
        std::vector<int> got;
        for_each_bit(x, [&](int i) { got.push_back(i); });
        BIT_CHECK(got == idx);
        const int i = (int)rnd(rng, 0, 63);
        BIT_CHECK(has_bit(x, i) == (((u >> i) & 1) == 1));
        const int n = (int)rnd(rng, 0, 63);
        unsigned long long rev = 0;
        std::string s;
        for (int j = n - 1; j >= 0; --j) s += ((u >> j) & 1) ? '1' : '0';
        for (int j = 0; j < n; ++j)
            if ((u >> j) & 1) rev |= 1ULL << (n - 1 - j);
        BIT_CHECK(bit_reverse(x, n) == (long long)rev && mask_to_string(x, n) == s);
        // グレイコード
        const long long g = gray_code(x);
        BIT_CHECK(gray_inverse(g) == x && bit_count(gray_code(x) ^ gray_code((long long)(u + 1))) == 1);
    }
    BIT_CHECK(mask_to_string(-1, 64) == std::string(64, '1') && mask_to_string(5, 0) == "" && bit_count(-1) == 64 && top_bit(-1) == 63);
}

static void test_enumerate(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 11);
        const long long full = (1LL << n) - 1;
        const long long mask = (long long)(rng() & (unsigned long long)full), mask2 = (long long)(rng() & (unsigned long long)full);
        // 部分集合（4 通りのフラグ）
        for (int fl = 0; fl < 4; ++fl) {
            const bool z = fl & 1, self = fl & 2;
            std::vector<long long> expect, got;
            for (long long s = 0; s <= full; ++s)
                if ((s & ~mask) == 0 && (s != 0 || z) && (s != mask || self)) expect.push_back(s);
            for_each_submask(mask, [&](long long s) { got.push_back(s); }, z, self);
            BIT_CHECK(got == expect);
        }
        {
            std::vector<long long> expect, got;
            for (long long s = 0; s <= full; ++s)
                if ((s & mask) == mask) expect.push_back(s);
            for_each_supermask(mask, n, [&](long long s) { got.push_back(s); });
            BIT_CHECK(got == expect);
        }
        // must ⊆ s ⊆ allowed（must = mask & mask2 なら必ず含まれる。そうでないときは空）
        const long long must = it % 3 == 0 ? (mask & mask2) : (long long)(rng() & (unsigned long long)full), allowed = mask;
        const int k = (int)rnd(rng, -1, n + 1);
        std::vector<long long> e1, e2, g1, g2;
        for (long long s = 0; s <= full; ++s)
            if ((s & must) == must && (s & ~allowed) == 0) {
                e1.push_back(s);
                if (bit_count(s) == k) e2.push_back(s);
            }
        for_each_mask_between(must, allowed, [&](long long s) { g1.push_back(s); });
        for_each_mask_between(must, allowed, k, [&](long long s) { g2.push_back(s); });
        BIT_CHECK(g1 == e1 && g2 == e2);
        // k 個の部分集合
        std::vector<long long> e3, g3;
        for (long long s = 0; s <= full; ++s)
            if (bit_count(s) == k) e3.push_back(s);
        for_each_subset_k(n, k, [&](long long s) { g3.push_back(s); });
        BIT_CHECK(g3 == e3);
    }
    // 大きいビット：63 ビット目まで
    long long c = 0, last = 0;
    for_each_subset_k(63, 62, [&](long long s) { c++, last = s; });
    BIT_CHECK(c == 63 && last == LMAX - 1);  // 1〜62 ビット目
    c = 0;
    for_each_subset_k(63, 1, [&](long long s) { c++, last = s; });
    BIT_CHECK(c == 63 && last == (1LL << 62));
    c = 0;
    for_each_mask_between(1LL << 62, LMAX, 2, [&](long long s) { c++, last = s; });
    BIT_CHECK(c == 62 && last == (3LL << 61));
    c = 0;
    for_each_supermask(LMAX - 1, 63, [&](long long s) { c++, last = s; });
    BIT_CHECK(c == 2 && last == LMAX);
    c = 0;
    for_each_submask(1LL << 62 | 1, [&](long long s) { c++, last = s; });
    BIT_CHECK(c == 4 && last == (1LL << 62 | 1));
}

static void test_counting(std::mt19937_64& rng, int& cases) {
    long long x = 0;
    Int128 pc = 0;
    std::vector<long long> cnt(63, 0);
    for (long long n = 0; n <= 300000; ++n, ++cases) {
        x ^= n, pc += bit_count(n);
        for (int k = 0; k < 20; ++k) cnt[k] += (n >> k) & 1;
        BIT_CHECK(xor_up_to(n) == x && popcount_sum_up_to(n) == pc);
        const int k = (int)(n % 23);
        BIT_CHECK(count_bit_up_to(n, k) == cnt[k]);
    }
    BIT_CHECK(xor_up_to(-1) == 0 && popcount_sum_up_to(-5) == 0 && count_bit_up_to(-1, 3) == 0 && count_bit_up_to(LMAX, 63) == 0);
    for (int it = 0; it < 300000 * TEST_SCALE; ++it, ++cases) {  // 大きい n は 1 つ前との差
        const long long n = rnd(rng, 1, LMAX);
        const int k = (int)rnd(rng, 0, 63);
        BIT_CHECK((xor_up_to(n) ^ xor_up_to(n - 1)) == n);
        BIT_CHECK(count_bit_up_to(n, k) - count_bit_up_to(n - 1, k) == (long long)((n >> k) & 1));
        BIT_CHECK(popcount_sum_up_to(n) - popcount_sum_up_to(n - 1) == bit_count(n));
        // 区間の XOR
        const long long l = rnd(rng, 0, 300), r = rnd(rng, 0, 300);
        long long e = 0;
        for (long long v = l; v < r; ++v) e ^= v;
        BIT_CHECK(xor_range(l, r) == e && xor_range_closed(l, r - 1) == e);
        const long long L = rnd(rng, 0, LMAX), R = L + std::min<long long>(LMAX - L, rnd(rng, 0, 50));
        long long e2 = 0;
        for (long long v = L;; ++v) {
            e2 ^= v;
            if (v == R) break;
        }
        BIT_CHECK(xor_range_closed(L, R) == e2);
    }
    for (int k = 0; k < 63; ++k) BIT_CHECK(count_bit_up_to(LMAX, k) == (1LL << 62));
    BIT_CHECK(popcount_sum_up_to(LMAX) == (Int128)63 << 62);
}

static void test_xor_sums(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 30);
        std::vector<long long> a(n);
        for (auto& v : a) v = it % 2 ? rnd(rng, 0, 15) : rnd(rng, 0, LMAX);
        Int128 e1 = 0, e2 = 0;
        for (int i = 0; i < n; ++i) {
            long long cur = 0;
            for (int j = i; j < n; ++j) {
                cur ^= a[j], e2 += cur;
                if (j > i) e1 += a[i] ^ a[j];
            }
        }
        BIT_CHECK(xor_pair_sum(a) == e1 && subarray_xor_sum(a) == e2);
        bool small = true;
        for (long long v : a) small &= v <= INT_MAX;
        if (small) {  // vector<int> でも
            const std::vector<int> ai(a.begin(), a.end());
            BIT_CHECK(xor_pair_sum(ai) == e1 && subarray_xor_sum(ai) == e2);
        }
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_single(rng, cases);
    test_enumerate(rng, cases);
    test_counting(rng, cases);
    test_xor_sums(rng, cases);
    std::cout << "[bit] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[bit] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    const int N = 16;
    for (long long m = 0; m < (1LL << N); ++m) for_each_submask(m, [&](long long s) { s1 += s; });  // 3^16 ≈ 4.3×10^7
    const long long t_sub = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for_each_subset_k(40, 6, [&](long long s) { s2 ^= s; });  // C(40, 6) ≈ 3.8×10^6
    for_each_mask_between(1, (1LL << 40) - 1, 6, [&](long long s) { s2 += s & 7; });
    const long long t_k = ms(t0);
    const int n = 1000000;
    std::vector<long long> a(n);
    for (auto& v : a) v = rnd(rng, 0, (1LL << 60) - 1);
    t0 = std::chrono::steady_clock::now();
    const Int128 s3 = xor_pair_sum(a) + subarray_xor_sum(a);
    const long long t_sum = ms(t0);
    t0 = std::chrono::steady_clock::now();
    Int128 s4 = 0;
    for (int i = 0; i < n; ++i) s4 += popcount_sum_up_to(a[i]) + count_bit_up_to(a[i], i % 60) + xor_up_to(a[i]);
    const long long t_cnt = ms(t0);
    std::cout << "[bit] benchmark: 全 mask の for_each_submask（n = 16、3^16 回）" << t_sub << " ms, for_each_subset_k(40, 6) + for_each_mask_between(k = 6) "
              << t_k << " ms\n"
              << "[bit] benchmark: xor_pair_sum + subarray_xor_sum（n = 10^6）" << t_sum << " ms, popcount_sum_up_to + count_bit_up_to + xor_up_to × 10^6 " << t_cnt
              << " ms  (check " << (long long)((s3 + s4 + s1 + s2) % 7) << ")\n";
}

}  // namespace bit_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace bit_demo {

{COMMON_DEMO}

static void demo() {
    banner("ビット（bit） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  bit_count / top_bit / low_bit / bit_length / has_bit / bit_indices / for_each_bit   1 つの数のビット\n"
              << "  bit_reverse / mask_to_string                                  反転 / \"0101\" の形\n"
              << "  for_each_submask / for_each_supermask / for_each_mask_between / for_each_subset_k   部分集合の列挙（小さい順）\n"
              << "  gray_code / gray_inverse                                      グレイコード\n"
              << "  xor_up_to / xor_range / xor_range_closed / count_bit_up_to / popcount_sum_up_to   0〜n の XOR・ビットの数\n"
              << "  xor_pair_sum / subarray_xor_sum                               XOR の和\n";

    api("【1】 bit_count(x) / top_bit(x) / low_bit(x) / bit_length(x) / has_bit(x, i) / bit_indices(x)", "int / bool / vector<int>",
        "1 の数 / 最上位・最下位の 1 の位置（0 なら −1）/ 桁数 / i ビット目 / 1 の位置の一覧", "O(1)");
    SHOW(bit_count(13), "13 = 1101");
    SHOW(top_bit(13), "最上位は 3 ビット目");
    SHOW(low_bit(12), "12 = 1100 の最下位の 1 は 2 ビット目（値なら 12 & -12 = 4）");
    SHOW(bit_length(0), "0 は 0 桁");
    SHOW(has_bit(13, 1), "1 ビット目は 0");
    SHOW(bit_indices(13), "1 のビットの位置");
    SHOW(mask_to_string(13, 6), "上の桁が左");
    SHOW(bit_reverse(1, 4), "0001 → 1000");

    api("【2】 for_each_submask(mask, f, include_zero, include_self) / for_each_supermask(mask, n, f)", "なし（f(s) を呼ぶ）",
        "mask の部分集合 / mask を含む n ビットの集合を小さい順に。0 や mask 自身を除くこともできる", "O(個数)");
    RUN(std::vector<long long> v; for_each_submask(11, [&](long long s) { v.push_back(s); }));
    SHOW(v, "ABC269 C（Submask）: 11 = 1011 の部分集合を小さい順に（答えそのもの）");
    RUN(v.clear(); for_each_submask(6, [&](long long s) { v.push_back(s); }, false, false));
    SHOW(v, "0 と自分自身を除く（空でない真部分集合）");
    RUN(v.clear(); for_each_supermask(5, 4, [&](long long s) { v.push_back(s); }));
    SHOW(v, "0101 を含む 4 ビットの集合");
    note("部分集合 DP：for (m = 0; m < 2^n; ++m) for_each_submask(m, …) は全部で 3^n 回（n = 16 で約 4×10^7）");

    api("【3】 for_each_mask_between(must, allowed, f) / for_each_mask_between(must, allowed, k, f) / for_each_subset_k(n, k, f)", "なし（f(s) を呼ぶ）",
        "必ず 1 のビット must と使ってよいビット allowed を決めた列挙（1 の数 k も指定できる）/ n 個から k 個選ぶ", "O(個数) / O(個数 × k)");
    RUN(v.clear(); for_each_mask_between(0b0001, 0b1011, [&](long long s) { v.push_back(s); }));
    SHOW(v, "0 ビット目は必ず 1、2 ビット目は必ず 0");
    RUN(v.clear(); for_each_mask_between(0b0001, 0b1111, 2, [&](long long s) { v.push_back(s); }));
    SHOW(v, "さらに 1 の数がちょうど 2");
    RUN(v.clear(); for_each_subset_k(4, 2, [&](long long s) { v.push_back(s); }));
    SHOW(v, "4 個から 2 個選ぶ（C(4, 2) = 6 通り）");

    head("【4】 bit 全探索（使い方関数の形）  ABC197 C（ORXOR）：区切り方を 2^(N−1) 通り全部試す");
    {
        auto solve = [](const std::vector<long long>& a) {
            const int n = (int)a.size();
            long long best = LLONG_MAX;
            for (long long cut = 0; cut < (1LL << (n - 1)); ++cut) {  // i ビット目が 1 なら a[i] と a[i+1] のあいだで区切る
                long long x = 0, o = 0;
                for (int i = 0; i < n; ++i) {
                    o |= a[i];
                    if (i == n - 1 || has_bit(cut, i)) x ^= o, o = 0;
                }
                best = std::min(best, x);
            }
            return best;
        };
        code("for (long long cut = 0; cut < (1LL << (n - 1)); ++cut)  // has_bit(cut, i) なら i と i+1 のあいだで区切る");
        std::cout << "    A = {1, 5, 7} → " << solve({1, 5, 7}) << "、A = {10, 10, 10} → " << solve({10, 10, 10}) << "、A = {1, 3, 3, 1} → " << solve({1, 3, 3, 1}) << "\n";
    }

    api("【5】 gray_code(i) / gray_inverse(g)", "long long", "i 番目のグレイコード（隣どうしは 1 ビットだけ違う）/ その逆", "O(1)");
    RUN(v.clear(); for (int i = 0; i < 8; ++i) v.push_back(gray_code(i)));
    SHOW(v, "0, 1, 3, 2, 6, 7, 5, 4");
    SHOW(gray_inverse(5), "5 は 6 番目");
    note("集合を 1 つずつ足し引きしながら全部の部分集合を回る（差分更新の bit 全探索）ときに使う");

    api("【6】 xor_up_to(n) / xor_range(l, r) / xor_range_closed(L, R) / count_bit_up_to(n, k) / popcount_sum_up_to(n)", "long long / Int128",
        "0 ⊕ … ⊕ n / 区間の XOR / 0〜n で k ビット目が 1 の数 / 0〜n の 1 の数の合計", "O(1) / O(63)");
    SHOW(xor_range_closed(2, 4), "ABC121 D: 2 ⊕ 3 ⊕ 4 = 5");
    SHOW(xor_range_closed(123, 456), "ABC121 D: 435");
    SHOW(xor_range(2, 5), "[2, 5) も同じ");
    SHOW(count_bit_up_to(10, 1), "0〜10 で 1 ビット目が 1：2, 3, 6, 7, 10");
    SHOW(popcount_sum_up_to(7), "0〜7 の 1 の数の合計 = 12");
    {
        // ABC356 D：Σ_{k=0}^{N} popcount(k & M) mod 998244353 = Σ_{M の i ビット目が 1} count_bit_up_to(N, i)
        auto solve = [](long long N, long long M) {
            long long s = 0;
            for_each_bit(M, [&](int i) { s = (s + count_bit_up_to(N, i)) % 998244353; });
            return s;
        };
        note("ABC356 D（Masked Popcount）: M の 1 のビット i について count_bit_up_to(N, i) を足す");
        std::cout << "    N = 4, M = 3 → " << solve(4, 3) << "、N = M = 2^60 − 1 → " << solve((1LL << 60) - 1, (1LL << 60) - 1) << "\n";
    }

    api("【7】 xor_pair_sum(a) / subarray_xor_sum(a)", "Int128",
        "Σ_{i<j} (a_i ⊕ a_j) / すべての区間の XOR の和（a_i ≥ 0。ビットごとに数える）", "O(63 n)");
    SHOW(xor_pair_sum({1, 2, 3}), "ABC147 D: 3 + 2 + 1 = 6（答えは % 1000000007）");
    SHOW((long long)(xor_pair_sum({3, 14, 159, 2653, 58979, 323846, 2643383, 27950288, 419716939, 9375105820LL}) % 1000000007), "ABC147 D の例 3");
    SHOW(subarray_xor_sum({1, 3, 2}), "長さ 1 の区間も含む：1 + 3 + 2 + 2 + 1 + 0 = 9");
    {
        const std::vector<long long> a = {2, 5, 6, 5, 2, 1, 7};
        Int128 s = subarray_xor_sum(a);
        for (long long x : a) s -= x;
        note("ABC365 E（Xor Sigma Problem）: 長さ 2 以上なので subarray_xor_sum(a) − Σ a_i");
        std::cout << "    A = {1, 3, 2} → " << subarray_xor_sum({1, 3, 2}) - 6 << "、A = {2, 5, 6, 5, 2, 1, 7} → " << s << "\n";
    }
    head("【8】 ビットの式の例  ABC238 D（AND and SUM）：x & y = a、x + y = s となる x, y ≥ 0 があるか");
    {
        auto ok = [](long long a, long long s) { return s - 2 * a >= 0 && ((s - 2 * a) & a) == 0; };  // x + y = 2(x & y) + (x ⊕ y)
        code("s − 2a ≥ 0 かつ ((s − 2a) & a) == 0   // x + y = 2(x & y) + (x ⊕ y)");
        std::cout << "    (a, s) = (1, 8) → " << (ok(1, 8) ? "Yes" : "No") << "、(4, 2) → " << (ok(4, 2) ? "Yes" : "No") << "\n";
    }
}

}  // namespace bit_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    bit_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    bit_test::benchmark();    // 速度計測
    bit_demo::demo();         // 使い方デモ
    return 0;
}
