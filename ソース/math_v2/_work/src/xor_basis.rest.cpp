// =============================================================================
//  ランダムテスト
// =============================================================================
namespace xorb_test {

{COMMON_TEST}

using u64 = unsigned long long;

static std::vector<u64> span_of(const std::vector<u64>& a) {  // 部分集合の XOR を全部（小さい順・重複なし）
    std::set<u64> s{0};
    for (u64 x : a) {
        std::set<u64> t = s;
        for (u64 y : s) t.insert(x ^ y);
        s.swap(t);
    }
    return std::vector<u64>(s.begin(), s.end());
}

static void test_basis(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 800 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 12), bits = (int)rnd(rng, 1, 7);
        std::vector<u64> a((size_t)n);
        for (auto& x : a) x = rng() & ((1ULL << bits) - 1);
        XorBasis B;
        std::vector<u64> pre;
        for (u64 x : a) {
            const size_t before = span_of(pre).size();
            pre.push_back(x);
            XORB_CHECK(B.insert(x) == (span_of(pre).size() > before));
        }
        const std::vector<u64> sp = span_of(a);
        ++cases;
        XORB_CHECK(B.count() == (Int128)sp.size() && (1ULL << B.rank()) == sp.size() && B.inserted() == n);
        XORB_CHECK(B.zero_by_nonempty() == (n > B.rank()));
        for (size_t k = 0; k < sp.size(); ++k) XORB_CHECK(B.kth(k) == sp[k]);
        for (u64 x = 0; x < (1ULL << bits) + 3; ++x, ++cases) {
            u64 mn = ~0ULL, mx = 0;
            for (u64 s : sp) mn = std::min(mn, x ^ s), mx = std::max(mx, x ^ s);
            XORB_CHECK(B.contains(x) == std::binary_search(sp.begin(), sp.end(), x));
            XORB_CHECK(B.min_xor(x) == mn && B.reduce(x) == mn && B.max_xor(x) == mx);
            XORB_CHECK(B.count_le(x) == (Int128)(std::upper_bound(sp.begin(), sp.end(), x) - sp.begin()));
            const Maybe<std::vector<int>> r = B.represent(x);
            XORB_CHECK(r.ok == B.contains(x));
            if (r.ok) {
                u64 y = 0;
                for (size_t i = 0; i < r.value.size(); ++i) {
                    XORB_CHECK(r.value[i] >= 0 && r.value[i] < n && (i == 0 || r.value[i - 1] < r.value[i]));
                    y ^= a[(size_t)r.value[i]];
                }
                XORB_CHECK(y == x);
            }
        }
        // merge：2 つに分けて作っても同じ
        XorBasis C, D;
        for (int i = 0; i < n; ++i) (i % 2 ? C : D).insert(a[(size_t)i]);
        C.merge(D);
        XORB_CHECK(C.basis() == B.basis());
    }
    // 64 ビットいっぱいの値
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        XorBasis B;
        const int n = (int)rnd(rng, 0, 80);
        std::vector<u64> a((size_t)n);
        for (auto& x : a) x = it % 3 ? rng() : (rng() & rng() & rng());
        for (u64 x : a) B.insert(x);
        const std::vector<u64> bs = B.basis();
        XORB_CHECK((int)bs.size() == B.rank());
        for (size_t i = 0; i < bs.size(); ++i) {  // 簡約：先頭ビットは他の基底に出てこない
            const int h = 63 - __builtin_clzll(bs[i]);
            for (size_t j = 0; j < bs.size(); ++j) XORB_CHECK(i == j || !(bs[j] >> h & 1));
        }
        const u64 x = rng();
        const Int128 c = B.count_le(x);
        XORB_CHECK(c >= 1 && B.kth((u64)(c - 1)) <= x && (c == B.count() || B.kth((u64)c) > x));
        const u64 k = B.rank() == 64 ? rng() : rng() % (1ULL << B.rank());
        XORB_CHECK(B.count_le(B.kth(k)) == (Int128)k + 1);
        u64 s = 0;  // 作れる値 → represent
        for (u64 v : a)
            if (rng() & 1) s ^= v;
        const Maybe<std::vector<int>> r = B.represent(s);
        XORB_CHECK(r.ok);
        u64 y = 0;
        for (int i : r.value) y ^= a[(size_t)i];
        XORB_CHECK(y == s && B.contains(s) && B.max_xor(s) >= s);
    }
}

static void test_prefix(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 40 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 40), bits = (int)rnd(rng, 1, 64);
        std::vector<u64> a((size_t)n);
        for (auto& x : a) x = bits == 64 ? rng() : rng() & ((1ULL << bits) - 1);
        PrefixXorBasis P;
        for (int r = 1; r <= n; ++r) {
            P.push_back(a[(size_t)r - 1]);
            XORB_CHECK(P.size() == r);
            for (int l = 0; l <= r; ++l, ++cases) {
                XorBasis B;
                for (int i = l; i < r; ++i) B.insert(a[(size_t)i]);
                const u64 x = rng() % 3 ? (rng() & (bits == 64 ? ~0ULL : (1ULL << bits) - 1)) : B.kth(B.rank() == 64 ? rng() : rng() % (1ULL << B.rank()));
                XORB_CHECK(P.contains(l, x) == B.contains(x) && P.max_xor(l, x) == B.max_xor(x) && P.rank(l) == B.rank());
                XORB_CHECK(P.basis(l).basis() == B.basis());
            }
        }
    }
}

static void test_gf2(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const int m = (int)rnd(rng, 0, 10), n = (int)rnd(rng, 0, 12);
        Gf2System S(m);
        std::vector<std::vector<int>> eq((size_t)n);
        std::vector<int> rhs((size_t)n);
        for (int i = 0; i < n; ++i) {
            eq[(size_t)i].assign((size_t)m, 0);
            for (auto& c : eq[(size_t)i]) c = (int)(rng() % 3 == 0);
            rhs[(size_t)i] = (int)(rng() & 1);
            if (i % 2) S.add_row(eq[(size_t)i], rhs[(size_t)i]);
            else {
                std::vector<int> vars;
                for (int j = 0; j < m; ++j)
                    if (eq[(size_t)i][(size_t)j]) vars.push_back(j), vars.push_back(j), vars.push_back(j);  // 3 回 = 1 回
                S.add(vars, rhs[(size_t)i]);
            }
        }
        auto sat = [&](const std::vector<int>& x, bool homo) {
            for (int i = 0; i < n; ++i) {
                int s = 0;
                for (int j = 0; j < m; ++j) s ^= eq[(size_t)i][(size_t)j] & x[(size_t)j];
                if (s != (homo ? 0 : rhs[(size_t)i])) return false;
            }
            return true;
        };
        long long cnt = 0;
        for (int mask = 0; mask < (1 << m); ++mask) {
            std::vector<int> x((size_t)m);
            for (int j = 0; j < m; ++j) x[(size_t)j] = mask >> j & 1;
            cnt += sat(x, false);
        }
        XORB_CHECK(S.solvable() == (cnt > 0) && S.equations() == n);
        if (cnt > 0) XORB_CHECK(cnt == (1LL << S.free_count()) && S.rank() == m - S.free_count());
        const Maybe<std::vector<int>> x = S.solve();
        XORB_CHECK(x.ok == (cnt > 0) && (!x.ok || sat(x.value, false)));
        const auto K = S.kernel();
        XORB_CHECK((int)K.size() == m - S.rank());
        XorBasis kb;  // 基底が一次独立か（m ≤ 10 なのでビットにして）
        for (const auto& v : K) {
            XORB_CHECK(sat(v, true));
            u64 b = 0;
            for (int j = 0; j < m; ++j) b |= (u64)v[(size_t)j] << j;
            XORB_CHECK(kb.insert(b));
        }
    }
    // 大きい：答えを決めてから式を作る
    for (int it = 0; it < 20 * TEST_SCALE; ++it, ++cases) {
        const int m = (int)rnd(rng, 1, 400), n = (int)rnd(rng, 1, 450);
        std::vector<int> ans((size_t)m);
        for (auto& v : ans) v = (int)(rng() & 1);
        Gf2System S(m);
        std::vector<std::vector<int>> eq;
        std::vector<int> rhs;
        for (int i = 0; i < n; ++i) {
            std::vector<int> c((size_t)m);
            int s = 0;
            for (int j = 0; j < m; ++j) c[(size_t)j] = (int)(rng() & 1), s ^= c[(size_t)j] & ans[(size_t)j];
            S.add_row(c, s), eq.push_back(c), rhs.push_back(s);
        }
        const Maybe<std::vector<int>> x = S.solve();
        XORB_CHECK(x.ok);
        for (int i = 0; i < n; ++i) {
            int s = 0;
            for (int j = 0; j < m; ++j) s ^= eq[(size_t)i][(size_t)j] & x.value[(size_t)j];
            XORB_CHECK(s == rhs[(size_t)i]);
        }
        S.add({0}, 1 - x.value[0]);  // 式を足すと解き直す
        XORB_CHECK(S.equations() == n + 1);
        if (S.solvable()) XORB_CHECK(S.solve().value[0] == 1 - x.value[0]);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_basis(rng, cases);
    test_prefix(rng, cases);
    test_gf2(rng, cases);
    std::cout << "[xor_basis] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[xor_basis] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    auto t0 = std::chrono::steady_clock::now();
    XorBasis B;
    for (int i = 0; i < 1000000; ++i) B.insert(rng() >> 4);
    u64 s1 = 0;
    for (int i = 0; i < 1000000; ++i) s1 ^= B.max_xor(rng()) ^ B.reduce(rng());
    const long long t_b = ms(t0);
    t0 = std::chrono::steady_clock::now();
    PrefixXorBasis P;
    u64 s2 = 0;
    for (int i = 0; i < 400000; ++i) {
        P.push_back(rng() >> 4);
        if (i % 2) s2 += P.contains((int)rnd(rng, 0, i), rng() >> 4);
    }
    const long long t_p = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const int M = 2000;
    Gf2System S(M);
    for (int i = 0; i < M; ++i) {
        std::vector<int> c((size_t)M);
        for (auto& v : c) v = (int)(rng() & 1);
        S.add_row(c, (int)(rng() & 1));
    }
    const int rk = S.rank();
    const long long t_g = ms(t0);
    std::cout << "[xor_basis] benchmark: XorBasis insert × 10^6 + クエリ × 2×10^6 " << t_b << " ms, PrefixXorBasis（4×10^5 個 + クエリ 2×10^5）" << t_p
              << " ms, Gf2System 2000 × 2000 " << t_g << " ms  (check " << (s1 + s2 + (u64)rk) % 7 << ")\n";
}

}  // namespace xorb_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace xorb_demo {

{COMMON_DEMO}

static void demo() {
    banner("XOR の基底（xor_basis） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  XorBasis：insert / contains / max_xor / min_xor / kth / count_le / rank / count / zero_by_nonempty / represent / merge / basis\n"
              << "  PrefixXorBasis：push_back / contains(l, x) / max_xor(l, x) / rank(l) / basis(l)\n"
              << "  Gf2System：add / add_row / rank / solvable / free_count / solve / kernel\n";

    api("【1】 XorBasis B; B.insert(x) / contains / max_xor / min_xor / rank / count", "bool / bool / u64 / u64 / int / Int128",
        "XOR で作れる値（0 個の XOR = 0 も含む）。insert は作れる値が増えたら true", "O(64)");
    RUN(XorBasis B);
    SHOW(B.insert(2), "");
    SHOW(B.insert(21), "");
    SHOW(B.insert(17), "");
    SHOW(B.insert(21), "21 は作れる（増えない）");
    SHOW(B.contains(6), "2 ^ 21 ^ 17 = 6");
    SHOW(B.contains(1), "");
    SHOW(B.max_xor(), "作れる値の最大");
    SHOW(B.max_xor(8), "8 XOR (作れる値) の最大");
    SHOW(B.min_xor(23), "23 XOR (作れる値) の最小");
    SHOW(B.rank(), "");
    SHOW(B.count(), "2^3 = 8 個");
    SHOW(B.zero_by_nonempty(), "21 ^ 21 = 0（1 個以上の XOR で 0 を作れる）");
    {
        // ABC141 F（Xor Sum 3）：全体の XOR T で立っているビットは必ず 1 回だけ足される。残りのビットは 2 倍 → T を消した値の max_xor
        auto solve = [](const std::vector<unsigned long long>& A) {
            unsigned long long T = 0;
            for (auto a : A) T ^= a;
            XorBasis Bs;
            for (auto a : A) Bs.insert(a & ~T);
            return T + 2 * Bs.max_xor();
        };
        note("ABC141 F（Xor Sum 3）: 全体の XOR を T として、各 A_i から T のビットを消した基底の max_xor × 2 + T");
        std::cout << "    例 1 → " << solve({3, 6, 5}) << "、例 2 → " << solve({23, 36, 66, 65}) << "、例 3 → "
                  << solve({1008288677408720767ULL, 539403903321871999ULL, 1044301017184589821ULL, 215886900497862655ULL, 504277496111605629ULL,
                            972104334925272829ULL, 792625803473366909ULL, 972333547668684797ULL, 467386965442856573ULL, 755861732751878143ULL,
                            1151846447448561405ULL, 467257771752201853ULL, 683930041385277311ULL, 432010719984459389ULL, 319104378117934975ULL,
                            611451291444233983ULL, 647509226592964607ULL, 251832107792119421ULL, 827811265410084479ULL, 864032478037725181ULL})
                  << "\n";
    }

    api("【2】 B.kth(k) / B.count_le(x)", "u64 / Int128", "作れる値を小さい順に並べた k 番目（0 から。0 番目は 0）/ x 以下の作れる値の個数", "O(64)");
    SHOW(B.kth(0), "");
    SHOW(B.kth(7), "最大");
    SHOW(B.count_le(18), "0, 2, 4, 6, 17 の 5 個");
    {
        // ABC283 G（Partial Xor Enumeration）：s_L〜s_R（1 から）= kth(L−1)〜kth(R−1)
        auto solve = [](long long L, long long R, const std::vector<unsigned long long>& A) {
            XorBasis Bs;
            for (auto a : A) Bs.insert(a);
            std::vector<unsigned long long> res;
            for (long long k = L - 1; k < R; ++k) res.push_back(Bs.kth((unsigned long long)k));
            return res;
        };
        note("ABC283 G（Partial Xor Enumeration）: kth(L − 1)〜kth(R − 1)");
        std::cout << "    例 1 → " << to_s(solve(1, 8, {2, 21, 17, 21})) << "\n    例 2 → " << to_s(solve(3, 7, {2, 21, 17, 21})) << "\n    例 3 → "
                  << to_s(solve(1, 1, {0, 0, 0, 0, 0})) << "、例 4 → " << to_s(solve(21, 21, {88, 44, 82, 110, 121, 80})) << "\n";
    }

    api("【3】 B.represent(x) / B.basis() / B.merge(C)", "Maybe<vector<int>> / vector<u64> / なし",
        "x を作る要素の番号（insert を呼んだ順に 0 から）/ 基底（簡約済み。大きい順）/ C で作れる値も加える", "O(64)");
    SHOW(B.represent(6).value, "2 ^ 21 ^ 17 → insert の 0, 1, 2 番目");
    SHOW(B.represent(1).ok, "作れない");
    SHOW(B.basis(), "簡約した基底");

    api("【4】 PrefixXorBasis P; P.push_back(a) / P.contains(l, x) / P.max_xor(l, x) / P.rank(l)", "bool / u64 / int",
        "a_0, a_1, … を順に足し、今の末尾 r = size() までの区間 [l, r) の要素で（クエリは r の順に並べて答える）", "O(64)");
    {
        // ABC223 H（Xor Query）：クエリを R の順に並べ、A_R まで足したところで contains(L − 1, X)
        auto solve = [](const std::vector<unsigned long long>& A, const std::vector<std::array<long long, 3>>& Q) {
            std::vector<int> ord(Q.size());
            std::iota(ord.begin(), ord.end(), 0);
            std::sort(ord.begin(), ord.end(), [&](int a, int b) { return Q[(size_t)a][1] < Q[(size_t)b][1]; });
            PrefixXorBasis P;
            std::vector<std::string> ans(Q.size());
            for (int q : ord) {
                while (P.size() < Q[(size_t)q][1]) P.push_back(A[(size_t)P.size()]);
                ans[(size_t)q] = P.contains((int)Q[(size_t)q][0] - 1, (unsigned long long)Q[(size_t)q][2]) ? "Yes" : "No";
            }
            return ans;
        };
        note("ABC223 H（Xor Query）: クエリを R の順に並べて、A_1〜A_R を push_back したところで contains(L − 1, X)");
        std::cout << "    例 1 → " << to_s(solve({3, 1, 4, 1, 5}, {{1, 3, 7}, {2, 5, 7}})) << "\n    例 2 → "
                  << to_s(solve({8, 45, 56, 9, 38, 28, 33, 5, 15, 19}, {{10, 10, 53}, {3, 8, 60}, {1, 10, 29}, {5, 7, 62}, {3, 7, 51}, {8, 8, 52}, {1, 4, 60}, {6, 8, 32}, {4, 8, 58}, {5, 9, 2}}))
                  << "\n";
    }

    api("【5】 Gf2System S(m); S.add(vars, rhs) / S.solve() / S.kernel() / S.free_count()", "なし / Maybe<vector<int>> / vector<vector<int>> / int",
        "mod 2 の連立一次方程式（vars の変数の XOR = rhs）。解の 1 つと、右辺 0 の解の基底", "O(式 × rank × m / 64)");
    RUN(Gf2System S(3));
    RUN(S.add({0, 1}, 1));
    RUN(S.add({1, 2}, 0));
    SHOW(S.solve().value, "x0 ^ x1 = 1、x1 ^ x2 = 0 の解の 1 つ（x2 = 0）");
    SHOW(S.kernel(), "(1, 1, 1) を XOR しても解");
    SHOW(S.free_count(), "解は 2^1 = 2 個");
    RUN(S.add({0, 2}, 0));
    SHOW(S.solvable(), "x0 ^ x2 = 0 を足すと x0 ^ x1 ^ x1 ^ x2 = 1 と矛盾");
}

}  // namespace xorb_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    xorb_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    xorb_test::benchmark();    // 速度計測
    xorb_demo::demo();         // 使い方デモ
    return 0;
}
