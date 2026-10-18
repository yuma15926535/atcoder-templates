// =============================================================================
//  ランダムテスト
// =============================================================================
namespace perm_test {

{COMMON_TEST}

template <class C>
static long long brute_inv(const C& a) {
    long long s = 0;
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = i + 1; j < a.size(); ++j) s += a[j] < a[i];
    return s;
}
static Int128 rnd128(std::mt19937_64& rng, Int128 hi) {  // [0, hi)
    const UInt128 x = ((UInt128)rng() << 64) | rng();
    return (Int128)(x % (UInt128)hi);
}

static void test_perm(std::mt19937_64& rng, int& cases) {
    using mint = ModInt998;
    for (int n = 0; n <= 7; ++n) {
        std::vector<int> p((size_t)n);
        std::iota(p.begin(), p.end(), 0);
        long long idx = 0;
        const long long off = rnd(rng, -5, 5), mul = rnd(rng, 1, 4);
        do {
            ++cases;
            PERM_CHECK(perm_rank(p) == idx && perm_unrank(n, idx) == p && perm_rank_mod<mint>(p) == mint(idx));
            std::vector<long long> q;  // 値をずらしても同じ番号
            for (int x : p) q.push_back(x * mul + off);
            PERM_CHECK(perm_rank(q) == idx);
            PERM_CHECK(inversion_count(p) == brute_inv(p) && perm_parity(p) == brute_inv(p) % 2 && perm_parity(q) == brute_inv(p) % 2);
            ++idx;
        } while (std::next_permutation(p.begin(), p.end()));
    }
    // n = 33：番号が Int128 いっぱいまで
    Int128 f33 = 1;
    for (int i = 1; i <= 33; ++i) f33 *= i;
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 33);
        Int128 fn = 1;
        for (int i = 1; i <= n; ++i) fn *= i;
        const Int128 k = it == 0 ? fn - 1 : rnd128(rng, fn);
        const std::vector<int> p = perm_unrank(n, k);
        PERM_CHECK(perm_rank(p) == k && perm_rank_mod<mint>(p) == mint((long long)(k % 998244353)));
    }
    // n > 33：前の方は動かない（k < 34!）
    for (int it = 0; it < 30 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 34, 300);
        const Int128 k = rnd128(rng, math_internal::PERM_INF);
        const std::vector<int> p = perm_unrank(n, k);
        for (int i = 0; i < n - 34; ++i) PERM_CHECK(p[(size_t)i] == i);
        PERM_CHECK(p[(size_t)n - 34] == n - 34 + (int)(k / f33));
        const std::vector<int> tail(p.end() - 33, p.end());
        PERM_CHECK(perm_rank(tail) == k % f33);
    }
    // 大きい n の mod 版：O(n^2) の愚直と
    for (int it = 0; it < 20 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 2000);
        std::vector<int> p((size_t)n);
        std::iota(p.begin(), p.end(), 0);
        std::shuffle(p.begin(), p.end(), rng);
        mint r = 0;
        ModInt107 r7 = 0;
        for (int i = 0; i < n; ++i) {
            int less = 0;
            for (int j = i + 1; j < n; ++j) less += p[(size_t)j] < p[(size_t)i];
            r = r * mint(n - i) + mint(less);
            r7 = r7 * ModInt107(n - i) + ModInt107(less);
        }
        PERM_CHECK(perm_rank_mod<mint>(p) == r && perm_rank_mod<ModInt107>(p) == r7);
    }
    // 転倒数：同じ値・文字列・負の値
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 60);
        std::vector<long long> a((size_t)n);
        for (auto& x : a) x = it % 2 ? rnd(rng, -3, 3) : rnd_ll(rng);
        PERM_CHECK(inversion_count(a) == brute_inv(a));
        std::string s;
        for (int i = 0; i < n; ++i) s += (char)('a' + rnd(rng, 0, 4));
        PERM_CHECK(inversion_count(s) == brute_inv(s));
    }
}

static void test_multiset(std::mt19937_64& rng, int& cases) {
    using mint = ModInt998;
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 8);
        std::vector<int> a((size_t)n);
        for (auto& x : a) x = (int)rnd(rng, 0, 3) * 7 - 5;
        std::sort(a.begin(), a.end());
        std::vector<int> orig = a;
        std::shuffle(orig.begin(), orig.end(), rng);
        long long idx = 0;
        do {
            ++cases;
            PERM_CHECK(multiset_perm_rank(a) == idx && multiset_perm_unrank(orig, idx) == a && multiset_perm_rank_mod<mint>(a) == mint(idx));
            ++idx;
        } while (std::next_permutation(a.begin(), a.end()));
    }
    {  // 文字列
        std::string s = "aabbc";
        long long idx = 0;
        do {
            ++cases;
            PERM_CHECK(multiset_perm_rank(s) == idx && multiset_perm_unrank(std::string("cbaba"), idx) == s);
            ++idx;
        } while (std::next_permutation(s.begin(), s.end()));
        PERM_CHECK(idx == 30);
    }
    // 大きい：60 文字（並べ方 ≈ 5.8×10^26）で行って戻る
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        std::string s;
        const int n = (int)rnd(rng, 1, 60), kinds = (int)rnd(rng, 1, 3);
        for (int i = 0; i < n; ++i) s += (char)('a' + rnd(rng, 0, kinds - 1));
        const Int128 r = multiset_perm_rank(s);
        PERM_CHECK(multiset_perm_unrank(s, r) == s && multiset_perm_rank_mod<mint>(s) == mint((long long)(r % 998244353)));
    }
    PERM_CHECK(multiset_perm_unrank(std::string(30, 'a') + std::string(30, 'b'), (Int128)118264581564861423LL) == std::string(30, 'b') + std::string(30, 'a'));
    // mod 版は n が大きくても（Horner の O(N^2) と）
    for (int it = 0; it < 10 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 1500);
        std::vector<int> a((size_t)n);
        for (auto& x : a) x = (int)rnd(rng, 0, 5);
        // 残りの並べ方 M を mint で持ち、M × (小さい値の個数) / 残りの長さ を足す（逆元は mint で）
        std::vector<int> cnt(6, 0);
        for (int x : a) ++cnt[(size_t)x];
        mint M = 1;
        for (int i = 1; i <= n; ++i) M *= mint(i);
        for (int c : cnt)
            for (int i = 1; i <= c; ++i) M /= mint(i);
        mint r = 0;
        for (int i = 0; i < n; ++i) {
            int less = 0;
            for (int v = 0; v < a[(size_t)i]; ++v) less += cnt[(size_t)v];
            r += M * mint(less) / mint(n - i);
            M = M * mint(cnt[(size_t)a[(size_t)i]]) / mint(n - i);
            --cnt[(size_t)a[(size_t)i]];
        }
        PERM_CHECK(multiset_perm_rank_mod<mint>(a) == r);
    }
}

static void test_comb(std::mt19937_64& rng, int& cases) {
    for (int n = 0; n <= 12; ++n)
        for (int r = 0; r <= n; ++r) {
            std::vector<std::vector<int>> all;
            std::vector<int> cur;
            std::function<void(int)> rec = [&](int x) {
                if ((int)cur.size() == r) { all.push_back(cur); return; }
                for (int y = x; y < n; ++y) cur.push_back(y), rec(y + 1), cur.pop_back();
            };
            rec(0);
            for (size_t i = 0; i < all.size(); ++i, ++cases) {
                PERM_CHECK(kth_combination(n, r, (Int128)i) == all[i] && combination_rank(n, all[i]) == (Int128)i);
                std::vector<int> sh = all[i];
                std::shuffle(sh.begin(), sh.end(), rng);
                PERM_CHECK(combination_rank(n, sh) == (Int128)i);
            }
        }
    // C(120, 60) ≈ 10^35：行って戻る
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 120), r = (int)rnd(rng, 0, n);
        const Int128 tot = math_internal::binom128(n, r);
        const Int128 k = it % 5 == 0 ? tot - 1 : rnd128(rng, tot);
        const std::vector<int> c = kth_combination(n, r, k);
        PERM_CHECK((int)c.size() == r && std::is_sorted(c.begin(), c.end()) && combination_rank(n, c) == k);
    }
    // C(n, r) が Int128 を超えても kth_combination は使える（前の方から決まる）
    for (int it = 0; it < 50 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 200, 3000), r = (int)rnd(rng, 100, n - 50);
        const Int128 k = rnd128(rng, math_internal::PERM_INF);
        const std::vector<int> c = kth_combination(n, r, k);
        PERM_CHECK((int)c.size() == r && std::adjacent_find(c.begin(), c.end(), std::greater_equal<int>()) == c.end() && c.back() < n);
        // 後ろの方を除いた部分（右端の 60 個の位置から選ぶ部分）は C(60 + …) で確かめにくいので、先頭が 0, 1, 2, … と始まることだけ
        PERM_CHECK(c[0] == 0);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_perm(rng, cases);
    test_multiset(rng, cases);
    test_comb(rng, cases);
    std::cout << "[permutation] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[permutation] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    using mint = ModInt998;
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 1000000;
    std::vector<int> p((size_t)N);
    std::iota(p.begin(), p.end(), 0);
    std::shuffle(p.begin(), p.end(), rng);
    std::vector<long long> a((size_t)N);
    for (auto& x : a) x = rnd(rng, 0, 1000000000);
    auto t0 = std::chrono::steady_clock::now();
    const long long inv = inversion_count(a);
    const long long t_inv = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const mint r = perm_rank_mod<mint>(p);
    const long long t_rank = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<int> q = perm_unrank(N, math_internal::PERM_INF - 12345);
    const long long t_un = ms(t0);
    std::vector<int> b((size_t)N);
    for (auto& x : b) x = (int)rnd(rng, 0, 9);
    t0 = std::chrono::steady_clock::now();
    const mint mr = multiset_perm_rank_mod<mint>(b);
    const long long t_ms = ms(t0);
    std::cout << "[permutation] benchmark: N = 10^6 で inversion_count " << t_inv << " ms, perm_rank_mod " << t_rank << " ms, perm_unrank " << t_un
              << " ms, multiset_perm_rank_mod " << t_ms << " ms  (check " << (inv + r.val() + q[777777] + mr.val()) % 7 << ")\n";
}

}  // namespace perm_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace perm_demo {

{COMMON_DEMO}

static void demo() {
    using mint = ModInt998;
    banner("順列・組合せの番号（permutation） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  inversion_count / perm_parity / perm_rank / perm_rank_mod / perm_unrank\n"
              << "  multiset_perm_rank / multiset_perm_rank_mod / multiset_perm_unrank / kth_combination / combination_rank\n"
              << "  （番号はすべて 0 から。この中では using mint = ModInt998; としている）\n";

    api("【1】 inversion_count(a) / perm_parity(a)", "long long / int", "転倒数：i < j かつ a_i > a_j の組の数（値は何でも。同じ値の組は数えない。隣どうしの交換でソートする回数）", "O(N log N)");
    SHOW(inversion_count(std::vector<int>{3, 1, 2}), "(3,1)・(3,2)");
    SHOW(inversion_count(std::string("banana")), "文字列も");
    SHOW(perm_parity(std::vector<int>{3, 1, 2}), "置換の偶奇（0 = 偶置換）。転倒数 2 の偶奇。値がすべて異なる列で O(N)");
    SHOW(perm_parity(std::string("dcba")), "転倒数 6 → 偶");
    note("「2 つを入れ替える」操作で並べ替えられるかは偶奇が決め手（15 パズル・3 つの巡回は偶置換、など）");
    {
        // ABC190 F（Shift and Inversions）：先頭の a_0 を末尾に回すと、転倒数は −a_0 + (N − 1 − a_0) 変わる
        auto solve = [](std::vector<long long> a) {
            const long long n = (long long)a.size();
            long long c = inversion_count(a);
            std::vector<long long> res;
            for (long long k = 0; k < n; ++k) res.push_back(c), c += (n - 1 - a[(size_t)k]) - a[(size_t)k];
            return res;
        };
        note("ABC190 F（Shift and Inversions）: 最初だけ inversion_count、あとは先頭を回すごとに差分で");
        std::cout << "    例 1 → " << to_s(solve({0, 1, 2, 3})) << "\n    例 2 → " << to_s(solve({0, 3, 1, 5, 4, 2, 9, 6, 8, 7})) << "\n";
    }

    api("【2】 perm_rank(p) / perm_rank_mod<mint>(p) / perm_unrank(n, k)", "Int128 / mint / vector<int>",
        "値がすべて異なる列が辞書順で何番目か（値は何でも。n ≤ 33。大きいなら mod）/ 0〜n−1 の順列の k 番目", "O(N log N)");
    SHOW(perm_rank(std::vector<int>{0, 1, 2}), "いちばん最初");
    SHOW(perm_rank(std::vector<int>{3, 1, 2}), "1〜3 の順列：123, 132, 213, 231, 312 → 4 番目（0 から）");
    SHOW(perm_unrank(4, 23), "0〜3 の順列の最後（4! − 1）");
    SHOW(perm_rank_mod<mint>(perm_unrank(100, 1000000007)), "n が大きいときは mod で");
    {
        // ABC276 C（Previous Permutation）：perm_unrank(N, perm_rank(P) − 1) に 1 を足す
        auto solve = [](const std::vector<int>& P) {
            std::vector<int> q = perm_unrank((int)P.size(), perm_rank(P) - 1);
            for (int& x : q) ++x;
            return q;
        };
        note("ABC276 C（Previous Permutation）: perm_unrank(N, perm_rank(P) − 1)（N ≤ 33 のとき。std::prev_permutation でもよい）");
        std::cout << "    例 1 → " << to_s(solve({3, 1, 2})) << "\n    例 2 → " << to_s(solve({9, 8, 6, 5, 10, 3, 1, 2, 4, 7})) << "\n";
    }

    api("【3】 multiset_perm_rank(a) / multiset_perm_rank_mod<mint>(a) / multiset_perm_unrank(a, k)", "Int128 / mint / a と同じ型",
        "同じ値があってよい列の並べ替え（同じものは区別しない）で何番目か / k 番目の並び（a の順番は何でもよい）", "O(N log N) / O(N × 値の種類)");
    SHOW(multiset_perm_rank(std::string("abba")), "aabb, abab, abba → 2");
    SHOW(multiset_perm_unrank(std::string("aabb"), 3), "4 番目（0 から 3）");
    SHOW(multiset_perm_rank(std::vector<int>{2, 2, 1, 1, 1}), "最後：5! / (2! 3!) − 1 = 9");
    {
        // ABC202 D（aab aba baa）：a が A 個、b が B 個の文字列で辞書順 K 番目（1 から）
        auto solve = [](int A, int B, long long K) { return multiset_perm_unrank(std::string(A, 'a') + std::string(B, 'b'), K - 1); };
        note("ABC202 D（aab aba baa）: multiset_perm_unrank(\"a…ab…b\", K − 1)");
        std::cout << "    (2, 2, 4) → " << solve(2, 2, 4) << "\n    (30, 30, 118264581564861424) → " << solve(30, 30, 118264581564861424LL) << "\n";
    }

    api("【4】 kth_combination(n, r, k) / combination_rank(n, c)", "vector<int> / Int128",
        "{0..n−1} から r 個選ぶ組（小さい順の列）の辞書順 k 番目 / 組 c が何番目か", "O(n × min(r, 130)) / O(n + r log r)");
    SHOW(kth_combination(5, 3, 0), "最初");
    SHOW(kth_combination(5, 3, 9), "最後（C(5, 3) − 1 = 9）");
    SHOW(combination_rank(5, std::vector<int>{0, 3, 4}), "012, 013, 014, 023, 024, 034 → 5");
}

}  // namespace perm_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    perm_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    perm_test::benchmark();    // 速度計測
    perm_demo::demo();         // 使い方デモ
    return 0;
}
