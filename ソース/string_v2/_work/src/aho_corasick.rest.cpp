// =============================================================================
//  ランダムテスト
// =============================================================================
namespace ac_test {

{COMMON_TEST}

static bool ends_with_at(const std::string& t, int end, const std::string& p) {  // t[0, end) の末尾が p か
    return (int)p.size() <= end && t.compare(end - p.size(), p.size(), p) == 0;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const int k = 1 + (int)(rng() % 6);
        std::vector<std::string> ps(k);
        for (auto& p : ps) p = rnd_str_shaped(rng, 1 + (int)(rng() % 5), sg);
        if (k >= 2 && rng() % 4 == 0) ps[1] = ps[0];  // 同じ文字列
        const std::string t = rnd_str_shaped(rng, (int)(rng() % 40), sg);
        AhoCorasick ac;
        for (int i = 0; i < k; ++i) AC_CHECK(ac.add(ps[i]) == i);
        ac.build();
        AC_CHECK(ac.pattern_count() == k);
        // 出現
        std::vector<long long> cnt(k, 0);
        std::vector<int> first(k, -1);
        std::vector<std::pair<int, int>> occ;  // (終わり, 番号) の順に並べる
        long long total = 0;
        for (int e = 1; e <= (int)t.size(); ++e)
            for (int i = 0; i < k; ++i)
                if (ends_with_at(t, e, ps[i])) {
                    ++cnt[i], ++total;
                    if (first[i] < 0) first[i] = e - (int)ps[i].size();
                    occ.push_back({e, i});
                }
        AC_CHECK(ac.count_each(t) == cnt && ac.count_total(t) == total && ac.first_occurrence_each(t) == first && ac.contains_any(t) == (total > 0));
        std::vector<std::pair<int, int>> got;
        ac.for_each_match(t, [&](int st, int id) { got.push_back({st + (int)ps[id].size(), id}); });
        std::vector<std::pair<int, int>> a = got, b = occ;
        std::sort(a.begin(), a.end()), std::sort(b.begin(), b.end());
        AC_CHECK(a == b);
        for (std::size_t q = 1; q < got.size(); ++q) AC_CHECK(got[q - 1].first <= got[q].first);  // 終わる位置の順
        // 状態：t の各接頭辞を読んだあと
        int s = 0;
        for (int e = 1; e <= (int)t.size(); ++e) {
            const int s2 = ac.next(s, t[e - 1]);
            AC_CHECK(s2 == ac.next_id(s, ac.alphabet().id(t[e - 1])));
            s = s2;
            long long mc = 0;
            int lm = 0, dep = 0;
            std::vector<int> ids;
            for (int i = 0; i < k; ++i)
                if (ends_with_at(t, e, ps[i])) ++mc, lm = std::max(lm, (int)ps[i].size()), ids.push_back(i);
            for (int d = 0; d <= e; ++d) {  // t[0, e) の接尾辞で、どれかのパターンの接頭辞であるものの最長
                const std::string suf = t.substr(e - d, d);
                for (const auto& p : ps)
                    if (p.compare(0, d, suf) == 0 && (int)p.size() >= d) dep = std::max(dep, d);
            }
            AC_CHECK(ac.match_count(s) == mc && ac.longest_match(s) == lm && ac.depth(s) == dep);
            auto mi = ac.match_ids(s);
            std::sort(mi.begin(), mi.end());
            AC_CHECK(mi == ids);
        }
        // 包含
        for (int i = 0; i < k; ++i) {
            std::vector<int> ex;
            for (int j = 0; j < k; ++j)
                if (j != i && ps[i].find(ps[j]) != std::string::npos) ex.push_back(j);
            AC_CHECK(ac.contained_patterns(i) == ex);
            AC_CHECK(ac.depth(ac.node_of(i)) == (int)ps[i].size() && ac.pattern_length(i) == (int)ps[i].size());
        }
        // 幅優先の順でリンクの先が先
        std::vector<int> pos(ac.size(), -1);
        const auto& ord = ac.bfs_order();
        AC_CHECK((int)ord.size() == ac.size() && ord[0] == 0);
        for (int q = 0; q < (int)ord.size(); ++q) pos[ord[q]] = q;
        for (int v = 1; v < ac.size(); ++v) AC_CHECK(pos[ac.link(v)] < pos[v] && ac.depth(ac.link(v)) < ac.depth(v));
        ++cases;
    }
}

static void test_ints(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const Alphabet A = Alphabet::ints(-1, 2);
        const int k = 1 + (int)(rng() % 4);
        std::vector<std::vector<int>> ps(k);
        for (auto& p : ps) {
            p.resize(1 + rng() % 4);
            for (auto& x : p) x = (int)(rng() % 3) - 1;
        }
        std::vector<long long> t(rng() % 30);
        for (auto& x : t) x = (long long)(rng() % 3) - 1;
        AhoCorasick ac(A);
        for (auto& p : ps) ac.add(p);
        ac.build();
        std::vector<long long> cnt(k, 0);
        for (int i = 0; i < k; ++i)
            for (int st = 0; st + (int)ps[i].size() <= (int)t.size(); ++st) cnt[i] += std::equal(ps[i].begin(), ps[i].end(), t.begin() + st);
        AC_CHECK(ac.count_each(t) == cnt);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_basic(rng, cases);
    test_ints(rng, cases);
    std::cout << "[aho_corasick] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[aho_corasick] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::string t = rnd_str(rng, 1000000, 4);
    std::vector<std::string> ps;
    long long tot = 0;
    while (tot < 1000000) {
        const int len = 1 + (int)(rng() % 12), st = (int)(rng() % (t.size() - len));
        ps.push_back(t.substr(st, len)), tot += len;
    }
    auto t0 = std::chrono::steady_clock::now();
    AhoCorasick ac;
    for (const auto& p : ps) ac.add(p);
    ac.build();
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto c = ac.count_each(t);
    const long long t_cnt = ms(t0);
    long long x = c[0] + ac.count_total(t);
    t0 = std::chrono::steady_clock::now();
    AhoCorasick big;
    std::vector<std::string> ps2;
    tot = 0;
    while (tot < 1000000) ps2.push_back(rnd_str(rng, 1 + (int)(rng() % 20), 26)), tot += (long long)ps2.back().size();
    for (const auto& p : ps2) big.add(p);
    big.build();
    x += big.count_each(rnd_str(rng, 1000000, 26))[0];
    const long long t_big = ms(t0);
    std::cout << "[aho_corasick] benchmark: パターンの長さの合計 10^6（σ = 4）build " << t_build << " ms, count_each（|t| = 10^6）" << t_cnt
              << " ms, σ = 26 で build + count_each " << t_big << " ms, ノード " << big.size() << "  (check " << x % 7 << ")\n";
}

}  // namespace ac_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace ac_demo {

{COMMON_DEMO}

static void demo() {
    banner("Aho–Corasick（aho_corasick） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  add / build                                          パターンを足す・作る\n"
              << "  count_each / count_total / contains_any / first_occurrence_each / for_each_match   文字列の中の出現\n"
              << "  contained_patterns                                   パターンどうしの包含\n"
              << "  next / next_id / match_count / longest_match / match_ids / depth / link / bfs_order   状態（DP 用）\n";

    RUN(AhoCorasick ac);
    RUN(for (const std::string p : {"i", "s", "a", "is", "missisippi"}) ac.add(p));
    RUN(ac.build());

    api("【1】 count_each(t) / count_total(t) / contains_any(t)", "vector<long long> / long long / bool", "パターンごとの出現回数 / その合計 / どれかが現れるか", "O(|t| + ノード数)");
    SHOW(ac.count_each("missisippi"), "ABC362 G（Count Substring Query）の入力例 1");
    SHOW(ac.count_total("missisippi"), "");
    SHOW(ac.contains_any("xyz"), "");

    api("【2】 first_occurrence_each(t) / for_each_match(t, f)", "vector<int> / なし（f(開始位置, 番号)）", "パターンごとの最初の出現（なければ -1）/ すべての出現（終わる位置の順）", "O(|t| + ノード数) / O(|t| + 出現の数)");
    SHOW(ac.first_occurrence_each("misa"), "");
    {
        std::string out;
        ac.for_each_match("misis", [&](int st, int id) { out += "(" + std::to_string(st) + "," + std::to_string(id) + ")"; });
        code("ac.for_each_match(\"misis\", f) で (開始位置, 番号) -> " + out);
    }

    api("【3】 contained_patterns(id)", "vector<int>", "パターン id の中に部分文字列として現れるほかのパターン（同じ文字列も）。包含関係の DAG・反鎖（ABC354 G は最大流で）", "O(|p| + 個数)");
    {
        AhoCorasick a2;
        for (const std::string p : {"atcoder", "at", "coder", "code"}) a2.add(p);
        a2.build();
        code("{\"atcoder\", \"at\", \"coder\", \"code\"}   // ABC354 G（Select Strings）の入力例 1");
        SHOW(a2.contained_patterns(0), "atcoder の中に at・coder・code");
        SHOW(a2.contained_patterns(2), "coder の中に code");
        note("同じ文字列は互いに含むので、ABC354 G では先に重複を除く（重みの大きい方を残す）");
    }

    api("【4】 状態：next(s, c) / match_count(s) / longest_match(s) / match_ids(s) / depth(s)", "int / long long / int / vector<int> / int",
        "読んだ文字列の「どれかのパターンの接頭辞である最長の接尾辞」が状態。そこで終わるパターンの数・最長・番号", "O(1)（match_ids は個数）");
    {
        int s = 0;
        for (char c : std::string("mis")) s = ac.next(s, c);
        code("\"mis\" を読んだ状態 s");
        SHOW(ac.match_count(s), "\"mis\" の末尾で終わる：is と s");
        SHOW(ac.depth(s), "\"mis\" は missisippi の接頭辞");
        s = ac.next(s, 'i'), s = ac.next(s, 's');
        SHOW(ac.match_ids(s), "\"misis\" の末尾で終わる：is, s（長い順）");
        SHOW(ac.longest_match(s), "");
    }
    {
        note("ABC268 Ex（Taboo）: 状態を進め、match_count(s) > 0 になったらその文字を * にして根（0）に戻す（貪欲）");
        AhoCorasick a3;
        for (const std::string p : {"abcd", "ijk", "ghi"}) a3.add(p);
        a3.build();
        int s = 0, ans = 0;
        for (char c : std::string("abcdefghijklmn")) {
            s = a3.next(s, c);
            if (a3.match_count(s) > 0) ++ans, s = 0;
        }
        note("  入力例 1 → " + to_s(ans));
    }
    {
        note("パターンを含まない長さ L の文字列の数（DP）: dp[状態] を next で進め、match_count(s) > 0 の状態へは行かない");
        AhoCorasick a4(Alphabet::of("ab"));
        a4.add("aa"), a4.add("bab");
        a4.build();
        std::vector<long long> dp(a4.size(), 0);
        dp[0] = 1;
        for (int step = 0; step < 6; ++step) {
            std::vector<long long> nd(a4.size(), 0);
            for (int v = 0; v < a4.size(); ++v)
                for (int c = 0; c < 2; ++c) {
                    const int u = a4.next_id(v, c);
                    if (a4.match_count(u) == 0) nd[u] += dp[v];
                }
            dp = nd;
        }
        long long total = 0;
        for (long long x : dp) total += x;
        note("  {aa, bab} を含まない a・b の長さ 6 の文字列 → " + to_s(total) + " 個（L が 10^18 なら math_v2 05_matrix で遷移行列を累乗）");
    }
}

}  // namespace ac_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    ac_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ac_test::benchmark();    // 速度計測
    ac_demo::demo();         // 使い方デモ
    return 0;
}
