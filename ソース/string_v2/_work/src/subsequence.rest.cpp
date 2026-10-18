// =============================================================================
//  ランダムテスト
// =============================================================================
namespace subq_test {

{COMMON_TEST}

// 全部の部分列（位置の選び方ごと）
static std::vector<std::string> all_subseq(const std::string& s, int min_step = 1) {
    std::vector<std::string> r;
    const int n = (int)s.size();
    for (int mask = 0; mask < (1 << n); ++mask) {
        int last = -1000;
        bool ok = true;
        std::string t;
        for (int i = 0; i < n; ++i)
            if (mask >> i & 1) {
                if (i - last < min_step) ok = false;
                last = i, t += s[i];
            }
        if (ok) r.push_back(t);
    }
    return r;
}

static void test_all(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 12), sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, n, sg), t = rnd_str_shaped(rng, (int)(rng() % 5), sg);
        const auto subs = all_subseq(s);
        std::set<std::string> dist(subs.begin(), subs.end());
        // NextTable
        const NextTable nx(s);
        for (int from = 0; from <= n; ++from) {
            const std::string rest = s.substr(from);
            const bool sub = is_subsequence(t, rest);
            SUBQ_CHECK(nx.is_subsequence(t, from) == sub);
            int k = 0;
            while (k <= (int)t.size() && is_subsequence(t.substr(0, k), rest)) ++k;
            SUBQ_CHECK(nx.match_length(t, from) == k - 1);
            int e = -1;  // 貪欲に当てはめた終わり
            if (sub) {
                int pos = from;
                for (char c : t) pos = (int)s.find(c, pos) + 1;
                e = pos;
            }
            SUBQ_CHECK(nx.match_end(t, from) == e);
            for (char c = 'a'; c < 'a' + sg; ++c) {
                const auto f = s.find(c, from);
                SUBQ_CHECK(nx.next(from, c) == (f == std::string::npos ? -1 : (int)f));
                const auto b = from ? s.rfind(c, from - 1) : std::string::npos;
                SUBQ_CHECK(nx.prev(from, c) == (b == std::string::npos ? -1 : (int)b));
            }
        }
        // 先頭・末尾から何文字
        const auto pre = prefix_match_lengths(s, t), suf = suffix_match_lengths(s, t);
        for (int i = 0; i <= n; ++i) {
            int a = 0, b = 0;
            while (a < (int)t.size() && is_subsequence(t.substr(0, a + 1), s.substr(0, i))) ++a;
            while (b < (int)t.size() && is_subsequence(t.substr(t.size() - b - 1), s.substr(i))) ++b;
            SUBQ_CHECK(pre[i] == a && suf[i] == b);
        }
        // 辞書順最小・最大
        for (int k = 0; k <= n; ++k) {
            std::string mn, mx;
            bool any = false;
            for (const auto& u : subs)
                if ((int)u.size() == k) {
                    if (!any || u < mn) mn = u;
                    if (!any || u > mx) mx = u;
                    any = true;
                }
            SUBQ_CHECK(lex_smallest_subsequence(s, k) == mn && lex_largest_subsequence(s, k) == mx);
        }
        // 部分列にならない最短
        {
            const Alphabet A = Alphabet::of(std::string("abc").substr(0, sg));
            std::string best;
            for (int len = 1; best.empty(); ++len) {
                std::vector<std::string> cand = {""};
                for (int l = 0; l < len; ++l) {
                    std::vector<std::string> nc;
                    for (const auto& w : cand)
                        for (int c = 0; c < sg; ++c) nc.push_back(w + (char)('a' + c));
                    cand = nc;
                }
                for (const auto& w : cand)
                    if (!dist.count(w)) { best = w; break; }
            }
            SUBQ_CHECK(shortest_non_subsequence(s, A) == best);
        }
        // 数
        long long occ = 0;
        for (const auto& u : subs) occ += u == t;
        SUBQ_CHECK(count_subsequence_occurrences<long long>(s, t) == occ);
        SUBQ_CHECK(count_distinct_subsequences<long long>(s) == (long long)dist.size() - 1 && count_distinct_subsequences<long long>(s, true) == (long long)dist.size());
        for (int g = 1; g <= 3; ++g) {
            const auto sg2 = all_subseq(s, g);
            std::set<std::string> d2(sg2.begin(), sg2.end());
            SUBQ_CHECK(count_distinct_subsequences<long long>(s, false, g) == (long long)d2.size() - 1);
        }
        const auto byl = distinct_subsequences_by_length<long long>(s);
        std::vector<long long> ex(n + 1, 0);
        for (const auto& u : dist) ++ex[u.size()];
        SUBQ_CHECK(byl == ex);
        // 数の列
        std::vector<long long> v(n);
        for (int i = 0; i < n; ++i) v[i] = (s[i] - 'a') * 1000000000000LL - 7;
        SUBQ_CHECK(count_distinct_subsequences<long long>(v) == (long long)dist.size() - 1);
        SUBQ_CHECK(lex_smallest_subsequence(v, n / 2).size() == (size_t)(n / 2));
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_all(rng, cases);
    // 途中の数が long long を超えても、答えが収まれば止まらない
    SUBQ_CHECK(count_subsequence_occurrences<long long>(std::string(67, 'a'), std::string(33, 'a') + "b") == 0);
    ++cases;
    std::cout << "[subsequence] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[subsequence] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 1000000;
    const std::string s = rnd_str(rng, n, 26);
    auto t0 = std::chrono::steady_clock::now();
    const NextTable nx(s);
    long long x = 0;
    for (int i = 0; i < 100000; ++i) x += nx.match_length(rnd_str(rng, 20, 26), (int)(rng() % n));
    const long long t_nx = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += (long long)lex_smallest_subsequence(s, n / 2).size() + (long long)shortest_non_subsequence(s).size();
    const long long t_misc = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += count_distinct_subsequences<long long>(rnd_str(rng, 50, 26)) + (long long)(count_subsequence_occurrences<double>(rnd_str(rng, 1000000, 3), "abcab") > 0);
    const long long t_cnt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += distinct_subsequences_by_length<double>(rnd_str(rng, 5000, 3))[100] > 0;
    const long long t_len = ms(t0);
    std::cout << "[subsequence] benchmark: NextTable（N = 10^6）+ match_length × 10^5 " << t_nx << " ms, 辞書順最小の部分列 + 部分列でない最短（N = 10^6）" << t_misc
              << " ms, 部分列の数（N = 10^6）" << t_cnt << " ms, 長さごとの数（N = 5000）" << t_len << " ms  (check " << x % 7 << ")\n";
}

}  // namespace subq_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace subq_demo {

{COMMON_DEMO}

static void demo() {
    banner("部分列（subsequence） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  NextTable nx(s, A)：next / prev / is_subsequence / match_length / match_end      次の出現と部分列の判定（O(1) / O(|t|)）\n"
              << "  prefix_match_lengths / suffix_match_lengths    t の先頭（末尾）が何文字まで入るか\n"
              << "  lex_smallest_subsequence / lex_largest_subsequence   長さ k の辞書順最小・最大の部分列\n"
              << "  shortest_non_subsequence                       部分列にならない最短の文字列\n"
              << "  count_subsequence_occurrences<T> / count_distinct_subsequences<T> / distinct_subsequences_by_length<T>\n";

    api("【1】 NextTable nx(s, A) / nx.next(i, c) / nx.is_subsequence(t, from) / nx.match_length(t, from)", "構造体",
        "次の出現の表（O(1)）。同じ s に何度も「t は部分列か」「先頭が何文字入るか」を聞くとき", "構築 O(Nσ)・判定 O(|t|)");
    RUN(const NextTable nx("abracadabra"));
    SHOW(nx.next(2, 'a'), "2 以上で最初の 'a'");
    SHOW(nx.is_subsequence("aaa"), "");
    SHOW(nx.is_subsequence("rad", 3), "s[3..) = \"acadabra\" には r の後に d がない");
    SHOW(nx.match_length("bard"), "b, a, r まで入る");
    SHOW(nx.match_end("abc"), "\"abc\" を当てはめた最後の文字（c は 4）の次");
    note("1 回だけなら 00_str_base の is_subsequence(t, s)。");

    api("【2】 prefix_match_lengths(s, t) / suffix_match_lengths(s, t)", "vector<int>（長さ |s| + 1）",
        "s の先頭 i 文字に t の先頭が何文字まで入るか / s[i..) に t の末尾が何文字まで入るか", "O(|s|)");
    SHOW(prefix_match_lengths("axbc", "abc"), "");
    SHOW(suffix_match_lengths("axbc", "abc"), "");
    {
        // ABC324 E: S_i + S_j が T を部分列に含む組 (i, j) の数 = pre(S_i) + suf(S_j) ≥ |T|
        const std::vector<std::string> S = {"abba", "bcb", "aaca"};
        const std::string T = "bac";
        std::vector<int> a, b;
        for (const auto& x : S) a.push_back(prefix_match_lengths(x, T).back()), b.push_back(suffix_match_lengths(x, T)[0]);
        std::sort(b.begin(), b.end());
        long long ans = 0;
        for (int x : a) ans += (long long)(b.end() - std::lower_bound(b.begin(), b.end(), (int)T.size() - x));
        note("ABC324 E（Joint Two Strings）: 各 S_i の先頭からの数 a_i と末尾からの数 b_j で a_i + b_j ≥ |T| の組を数える");
        std::cout << "    入力例 1 → " << ans << "\n";
    }

    api("【3】 lex_smallest_subsequence(s, k) / lex_largest_subsequence(s, k)", "列", "長さ k の部分列のうち辞書順最小 / 最大（スタック）", "O(N)");
    SHOW(lex_smallest_subsequence("atcoder", 3), "典型90 006（Smallest Subsequence）の入力例 1");
    SHOW(lex_largest_subsequence("atcoder", 3), "");

    api("【4】 shortest_non_subsequence(s, A)", "string", "A の文字でできる文字列のうち s の部分列にならない最短のもの（その中で辞書順最小）", "O(Nσ)");
    SHOW(shortest_non_subsequence("atcoderregularcontest"), "ARC081 E（Don't Be a Subsequence）の入力例 1");
    SHOW(shortest_non_subsequence("abcabc", Alphabet::of("abc")), "");

    api("【5】 count_subsequence_occurrences<T>(s, t)", "T", "t が s の部分列として現れる数（位置の選び方の数）", "O(|s||t|)");
    SHOW(count_subsequence_occurrences<long long>("chchokudaiokudai", "chokudai"), "ABC211 C の形（入力例は mod 10^9+7）");
    note("答えを mod で出すときは math_v2 の ModInt107 などを T に（典型90 008 AtCounter・ABC211 C chokudai）。");

    api("【6】 count_distinct_subsequences<T>(s, include_empty, min_step) / distinct_subsequences_by_length<T>(s)", "T / vector<T>",
        "異なる部分列の数（同じ文字列は 1 回。min_step は選ぶ位置の間隔の最小）/ 長さごとの数", "O(N) / O(N^2)");
    SHOW(count_distinct_subsequences<long long>("abab"), "a, b, aa, ab, ba, bb, aab, aba, abb, bab, abab");
    SHOW(count_distinct_subsequences<long long>("abc", false, 2), "隣どうしを選ばない（ABC214 F の入力例 1）");
    SHOW(distinct_subsequences_by_length<long long>("abab"), "r[len]");
    note("ABC230 F（Predilection）は累積和の列（最後を除く）の異なる部分列の数 + 1。数の列もそのまま渡せる。");
}

}  // namespace subq_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    subq_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    subq_test::benchmark();    // 速度計測
    subq_demo::demo();         // 使い方デモ
    return 0;
}
