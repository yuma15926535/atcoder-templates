// =============================================================================
//  ランダムテスト
// =============================================================================
namespace shash_test {

{COMMON_TEST}

static std::vector<long long> rnd_vals(std::mt19937_64& rng, int n) {
    std::vector<long long> v(n);
    const int t = (int)(rng() % 3);
    const long long range = 1 + (long long)(rng() % 6);
    for (auto& x : v) x = t == 2 ? rnd_ll(rng) : (long long)(rng() % range) - (t == 1 ? 3 : 0);
    return v;
}

static void test_prefix(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 12), m = (int)(rng() % 12);
        const auto a = rnd_vals(rng, n);
        std::vector<long long> b = rng() % 2 ? rnd_vals(rng, m) : a;
        if (rng() % 2) std::shuffle(b.begin(), b.end(), rng);
        const auto sa = prefix_set_hashes(a), sb = prefix_set_hashes(b);
        const auto ma = prefix_multiset_hashes(a), mb = prefix_multiset_hashes(b);
        SHASH_CHECK((int)sa.size() == n + 1 && (int)mb.size() == (int)b.size() + 1);
        for (int x = 0; x <= n; ++x)
            for (int y = 0; y <= (int)b.size(); ++y) {
                const std::set<long long> A(a.begin(), a.begin() + x), B(b.begin(), b.begin() + y);
                SHASH_CHECK((sa[x] == sb[y]) == (A == B));
                std::multiset<long long> MA(a.begin(), a.begin() + x), MB(b.begin(), b.begin() + y);
                SHASH_CHECK((ma[x] == mb[y]) == (MA == MB));
            }
        // 区間
        const MultisetPrefix PA(a), PB(b);
        for (int l1 = 0; l1 <= n; ++l1)
            for (int r1 = l1; r1 <= n; ++r1) {
                const int l2 = (int)(rng() % (b.size() + 1)), r2 = l2 + (int)(rng() % (b.size() - l2 + 1));
                std::vector<long long> x(a.begin() + l1, a.begin() + r1), y(b.begin() + l2, b.begin() + r2);
                std::sort(x.begin(), x.end()), std::sort(y.begin(), y.end());
                SHASH_CHECK((PA.get(l1, r1) == PB.get(l2, r2)) == (x == y));
                const int l3 = (int)(rng() % (n + 1)), r3 = l3 + (int)(rng() % (n - l3 + 1));
                std::vector<long long> z(a.begin() + l3, a.begin() + r3);
                std::sort(z.begin(), z.end());
                SHASH_CHECK(PA.same(l1, r1, l3, r3) == (x == z));
            }
        ++cases;
    }
    // 文字列でも同じ
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const std::string s = rnd_str(rng, (int)(rng() % 20), 3);
        const std::vector<long long> v(s.begin(), s.end());
        SHASH_CHECK(prefix_set_hashes(s) == prefix_set_hashes(v) && prefix_multiset_hashes(s) == prefix_multiset_hashes(v));
        ++cases;
    }
}

static void test_dynamic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        SetHash s1, s2;
        MultisetHash m1, m2;
        std::map<long long, long long> c1, c2;
        const long long range = 1 + (long long)(rng() % 5);
        for (int q = 0; q < 100; ++q) {
            const long long x = (long long)(rng() % range) - 2;
            const int who = (int)(rng() % 2), op = (int)(rng() % 3);
            SetHash& s = who ? s2 : s1;
            MultisetHash& m = who ? m2 : m1;
            auto& c = who ? c2 : c1;
            if (op < 2) {
                s.add(x), m.add(x), ++c[x];
            } else {
                const bool had = c.count(x) && c[x] > 0;
                SHASH_CHECK(s.remove(x) == had);
                if (had) {
                    m.remove(x);
                    if (--c[x] == 0) c.erase(x);
                }
            }
            std::set<long long> k1, k2;
            for (auto [v, k] : c1) k1.insert(v);
            for (auto [v, k] : c2) k2.insert(v);
            SHASH_CHECK((s1 == s2) == (k1 == k2) && (m1 == m2) == (c1 == c2));
            SHASH_CHECK(s1.distinct() == (long long)k1.size() && s1.count(x) == (c1.count(x) ? c1[x] : 0));
            // 多重集合の和
            MultisetHash u = m1;
            u.merge(m2);
            MultisetHash w;
            for (auto [v, k] : c1) w.add(v, k);
            for (auto [v, k] : c2) w.add(v, k);
            SHASH_CHECK(u == w);
        }
        ++cases;
    }
    // 大きな個数・負の個数
    {
        MultisetHash a, b;
        a.add(5, 1000000000000LL), a.add(7, -3), a.add(7, 3);
        for (int i = 0; i < 10; ++i) b.add(5, 100000000000LL);
        SHASH_CHECK(a == b && a.count == 1000000000000LL);
        a.add(5, std::numeric_limits<long long>::min()), a.add(5, std::numeric_limits<long long>::max()), a.add(5, 1);
        SHASH_CHECK(a == b);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_prefix(rng, cases);
    test_dynamic(rng, cases);
    std::cout << "[set_hash] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[set_hash] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 1000000;
    std::vector<long long> a(N);
    for (auto& x : a) x = (long long)(rng() % 1000000000);
    auto t0 = std::chrono::steady_clock::now();
    const auto ps = prefix_set_hashes(a);
    const long long t_set = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const MultisetPrefix mp(a);
    unsigned long long x = ps[N / 2];
    for (int q = 0; q < N; ++q) {
        int l = (int)(rng() % (N + 1)), r = (int)(rng() % (N + 1));
        if (l > r) std::swap(l, r);
        x += mp.get(l, r);
    }
    const long long t_ms = ms(t0);
    t0 = std::chrono::steady_clock::now();
    SetHash sh;
    for (int q = 0; q < N; ++q) {
        if (q % 3 == 2) sh.remove(a[rng() % q]);
        else sh.add(a[q]);
    }
    x += sh.value;
    const long long t_dyn = ms(t0);
    std::cout << "[set_hash] benchmark: prefix_set_hashes（N = 10^6）" << t_set << " ms, MultisetPrefix 作成 + get × 10^6 " << t_ms << " ms, SetHash 追加・削除 10^6 回 "
              << t_dyn << " ms  (check " << x % 7 << ")\n";
}

}  // namespace shash_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace shash_demo {

{COMMON_DEMO}

static void demo() {
    banner("集合・多重集合のハッシュ（set_hash） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  value_hash / prefix_set_hashes / prefix_multiset_hashes / MultisetPrefix / MultisetHash / SetHash\n"
              << "  値ごとに乱数を割り当てて和を取る（mod 2^61 - 1）。等しいか だけを O(1) で比べる（大小は意味がない）\n";

    api("【1】 prefix_set_hashes(a)", "vector<unsigned long long>（長さ N + 1）", "r[i] = a[0, i) に現れる値の集合のハッシュ（何回現れても 1 回）", "O(N)（平均）");
    {
        const std::vector<int> A = {1, 2, 3, 4, 5}, B = {1, 2, 2, 4, 3};
        code("const vector<int> A = {1, 2, 3, 4, 5}, B = {1, 2, 2, 4, 3};   // ABC250 E（Prefix Equality）の入力例 1");
        const auto ha = prefix_set_hashes(A), hb = prefix_set_hashes(B);
        std::string out;
        const int qs[][2] = {{1, 1}, {2, 2}, {2, 3}, {3, 3}, {4, 4}, {4, 5}, {5, 5}};
        for (const auto& q : qs) out += ha[q[0]] == hb[q[1]] ? "Yes " : "No ";
        note("先頭 x 個と先頭 y 個が集合として同じか: ha[x] == hb[y] → " + out);
    }

    api("【2】 prefix_multiset_hashes(a) / MultisetPrefix mp(a)", "vector<unsigned long long> / mp.get(l, r)・mp.same(l1, r1, l2, r2)",
        "多重集合（個数も見る）のハッシュ。区間 [l, r) は差で O(1)。並べ替えて同じ ⇔ 等しい", "O(N) / O(1)");
    {
        const std::vector<int> A = {1, 2, 3, 2, 4}, B = {2, 3, 1, 4, 2};
        code("const vector<int> A = {1, 2, 3, 2, 4}, B = {2, 3, 1, 4, 2};   // ABC367 F（Rearrange Query）の入力例 1");
        const MultisetPrefix PA(A), PB(B);
        std::string out;
        const int qs[][4] = {{1, 3, 1, 3}, {1, 2, 3, 5}, {1, 4, 2, 5}, {1, 5, 1, 5}};
        for (const auto& q : qs) out += PA.get(q[0] - 1, q[1]) == PB.get(q[2] - 1, q[3]) ? "Yes " : "No ";
        note("A の l〜r 番目と B の L〜R 番目: PA.get(l - 1, r) == PB.get(L - 1, R) → " + out);
        SHOW(PA.same(0, 2, 1, 3), "A[0, 2) = {1, 2} と A[1, 3) = {2, 3}");
    }

    api("【3】 MultisetHash h", "h.add(x, k) / remove(x, k) / merge(o) / ==", "足し引きできる多重集合（個数の表は持たない。value と count だけ）", "O(1)");
    {
        RUN(MultisetHash p, q);
        RUN(p.add(3), p.add(5, 2), q.add(5), q.add(3), q.add(5));
        SHOW(p == q, "{3, 5, 5} と {5, 3, 5}");
        RUN(q.remove(5));
        SHOW(p == q, "");
        SHOW(p.count, "要素の数");
    }

    api("【4】 SetHash h", "h.add(x) / remove(x) / count(x) / distinct() / ==", "足し引きできる集合（値ごとの個数を持ち、0 個 ⇔ 1 個以上 が変わるときだけ変わる）", "O(1)（平均）");
    {
        RUN(SetHash s, t);
        RUN(s.add(1), s.add(2), s.add(2), t.add(2), t.add(1));
        SHOW(s == t, "{1, 2} と {2, 1}");
        SHOW(s.remove(2), "1 つ取り除いても 2 はまだある");
        SHOW(s == t, "");
        SHOW(s.remove(9), "入っていない → false");
        SHOW(s.distinct(), "");
    }

    api("【5】 value_hash(x)", "unsigned long long", "値 x に割り当てた乱数（自分で和を取るとき。2 つの値の組 (a, b) は a × (b の上限) + b を 1 つの値にして渡す）", "O(1)");
    SHOW(value_hash(42) == value_hash(42LL), "同じ値には同じ乱数（実行ごとに変わる）");
    note("文字列の集合は 01_rolling_hash の hash_of(s).value を値にして SetHash に入れる");
}

}  // namespace shash_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    shash_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    shash_test::benchmark();    // 速度計測
    shash_demo::demo();         // 使い方デモ
    return 0;
}
