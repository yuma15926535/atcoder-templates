// =============================================================================
//  ランダムテスト
// =============================================================================
namespace lcs_test {

{COMMON_TEST}

static std::vector<std::vector<int>> bf_lcs_table(const std::string& a, const std::string& b) {
    const int n = (int)a.size(), m = (int)b.size();
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= m; ++j) dp[i][j] = a[i - 1] == b[j - 1] ? dp[i - 1][j - 1] + 1 : std::max(dp[i - 1][j], dp[i][j - 1]);
    return dp;
}
static long long bf_edit(const std::string& a, const std::string& b, long long ins, long long del, long long sub) {
    const int n = (int)a.size(), m = (int)b.size();
    std::vector<std::vector<long long>> dp(n + 1, std::vector<long long>(m + 1));
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= m; ++j) {
            if (i == 0) dp[i][j] = j * ins;
            else if (j == 0) dp[i][j] = i * del;
            else dp[i][j] = std::min({dp[i - 1][j] + del, dp[i][j - 1] + ins, dp[i - 1][j - 1] + (a[i - 1] == b[j - 1] ? 0 : sub)});
        }
    return dp[n][m];
}

static void test_lcs(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 4);
        const int n = (int)(rng() % (it % 10 == 0 ? 200 : 15)), m = (int)(rng() % (it % 10 == 0 ? 200 : 15));
        const std::string a = rnd_str_shaped(rng, n, sg), b = rnd_str_shaped(rng, m, sg);
        const auto dp = bf_lcs_table(a, b);
        const int L = dp[n][m];
        LCS_CHECK(lcs_length(a, b) == L);
        const std::string r = lcs_restore(a, b);
        LCS_CHECK((int)r.size() == L && is_subsequence(r, a) && is_subsequence(r, b));
        const std::string sc = shortest_common_supersequence(a, b);
        LCS_CHECK((int)sc.size() == n + m - L && is_subsequence(a, sc) && is_subsequence(b, sc));
        // 数の列・混ぜる
        const std::vector<long long> va(a.begin(), a.end());
        LCS_CHECK(lcs_length(va, b) == L && (int)lcs_restore(va, b).size() == L);
        ++cases;
    }
    // 長い列（64 ビットの境目をまたぐ）
    for (int it = 0; it < 30 * TEST_SCALE; ++it) {
        const std::string a = rnd_str(rng, 64 * (1 + (int)(rng() % 5)) + (int)(rng() % 3) - 1, 2), b = rnd_str(rng, 64 * (1 + (int)(rng() % 5)) + (int)(rng() % 3) - 1, 2);
        LCS_CHECK(lcs_length(a, b) == bf_lcs_table(a, b)[a.size()][b.size()]);
        ++cases;
    }
}

static void test_edit(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string a = rnd_str_shaped(rng, (int)(rng() % 12), sg);
        std::string b = rng() % 2 ? rnd_str_shaped(rng, (int)(rng() % 12), sg) : a;
        for (int t = (int)(rng() % 3); t > 0 && !b.empty(); --t) {  // a を少しだけ変える
            const int p = (int)(rng() % b.size());
            if (rng() % 3 == 0) b.erase(p, 1);
            else if (rng() % 2) b.insert(b.begin() + p, (char)('a' + rng() % sg));
            else b[p] = (char)('a' + rng() % sg);
        }
        const long long ins = 1 + (long long)(rng() % 3), del = 1 + (long long)(rng() % 3), sub = 1 + (long long)(rng() % 5);
        const long long d1 = bf_edit(a, b, 1, 1, 1), dw = bf_edit(a, b, ins, del, sub);
        LCS_CHECK(edit_distance(a, b) == d1 && edit_distance(a, b, ins, del, sub) == dw);
        // 操作の列を当てはめると b になり、費用が最小
        for (int w = 0; w < 2; ++w) {
            const auto ops = w ? edit_script(a, b, ins, del, sub) : edit_script(a, b);
            std::string out;
            long long cost = 0;
            int pa = 0, pb = 0;
            for (const EditOp& op : ops) {
                LCS_CHECK(pa <= (int)a.size() && pb <= (int)b.size());
                if (op.type == '=') {
                    LCS_CHECK(op.i == pa && op.j == pb && pa < (int)a.size() && pb < (int)b.size() && a[pa] == b[pb]);
                    out += a[pa++], ++pb;
                } else if (op.type == 'S') {
                    LCS_CHECK(op.i == pa && op.j == pb && pa < (int)a.size() && pb < (int)b.size() && a[pa] != b[pb]);
                    out += b[pb], ++pa, ++pb, cost += w ? sub : 1;
                } else if (op.type == 'D') {
                    LCS_CHECK(op.i == pa && pa < (int)a.size());
                    ++pa, cost += w ? del : 1;
                } else {
                    LCS_CHECK(op.type == 'I' && op.j == pb && pb < (int)b.size());
                    out += b[pb++], cost += w ? ins : 1;
                }
            }
            LCS_CHECK(out == b && pa == (int)a.size() && cost == (w ? dw : d1));
        }
        for (int k = 0; k <= 6; ++k) LCS_CHECK(edit_distance_bounded(a, b, k) == (d1 <= k ? (int)d1 : -1));
        LCS_CHECK(within_one_edit(a, b) == (d1 <= 1));
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_lcs(rng, cases);
    test_edit(rng, cases);
    // 費用が大きくても溢れない・k が大きくても速い
    {
        const std::string x(20, 'x');
        LCS_CHECK(edit_distance(x, x, 1000000000000000000LL, 1000000000000000000LL, 1) == 0);
        const auto ops = edit_script(x, x, 1000000000000000000LL, 1000000000000000000LL, 1);
        LCS_CHECK(ops.size() == 20 && std::all_of(ops.begin(), ops.end(), [](const EditOp& o) { return o.type == '='; }));
        LCS_CHECK(edit_distance("ab", "ba", 4000000000000000000LL, 4000000000000000000LL, 3) == 6);
        LCS_CHECK(edit_distance_bounded("abc", "abd", 2147483647) == 1 && edit_distance_bounded("abc", "xyz", 100000000) == 3);
        ++cases;
    }
    std::cout << "[lcs_edit] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[lcs_edit] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::string a = rnd_str(rng, 50000, 4), b = rnd_str(rng, 50000, 4);
    auto t0 = std::chrono::steady_clock::now();
    long long x = lcs_length(a, b);
    const long long t_len = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += (long long)lcs_restore(a.substr(0, 30000), b.substr(0, 30000)).size();
    const long long t_res = ms(t0);
    const std::string c = rnd_str(rng, 10000, 26), d = rnd_str(rng, 10000, 26);
    t0 = std::chrono::steady_clock::now();
    x += edit_distance(c, d);
    const long long t_ed = ms(t0);
    const std::string e = rnd_str(rng, 500000, 26);
    std::string f = e;
    for (int i = 0; i < 15; ++i) f[rng() % f.size()] = 'z';
    t0 = std::chrono::steady_clock::now();
    x += edit_distance_bounded(e, f, 20);
    const long long t_bd = ms(t0);
    std::cout << "[lcs_edit] benchmark: lcs_length（5×10^4 × 5×10^4）" << t_len << " ms, lcs_restore（3×10^4 × 3×10^4）" << t_res << " ms, edit_distance（10^4 × 10^4）"
              << t_ed << " ms, edit_distance_bounded（5×10^5、k = 20）" << t_bd << " ms  (check " << x % 7 << ")\n";
}

}  // namespace lcs_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace lcs_demo {

{COMMON_DEMO}

static void demo() {
    banner("LCS・編集距離（lcs_edit） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  lcs_length / lcs_restore / shortest_common_supersequence   最長共通部分列（飛び飛び）・それを含む最短の列\n"
              << "  edit_distance / edit_script / edit_distance_bounded / within_one_edit   編集距離・操作の列・k 以内・1 回以内\n"
              << "  連続する部分の最長共通（最長共通部分文字列）は 03_suffix_array の longest_common_substring\n";

    api("【1】 lcs_length(a, b) / lcs_restore(a, b)", "int / string",
        "最長共通部分列（飛び飛びに選んだ共通の列）の長さと、その 1 つ（ビット並列なので 5×10^4 × 5×10^4 でも速い）", "O(|a||b|/64)");
    SHOW(lcs_length("axyb", "abyxb"), "EDPC F（LCS）の入力例 1 の形");
    SHOW(lcs_restore("axyb", "abyxb"), "");
    SHOW(lcs_restore("abracadabra", "avadakedavra"), "");

    api("【2】 shortest_common_supersequence(a, b)", "string", "a と b の両方を部分列に含む最短の列（長さ |a| + |b| - LCS）", "O(|a||b|/64)");
    SHOW(shortest_common_supersequence("abac", "cab"), "");

    api("【3】 edit_distance(a, b, ins, del, sub)", "long long",
        "a を b にする最小の費用（挿入 ins・削除 del・置換 sub。既定はすべて 1 = レーベンシュタイン距離）", "O(|a||b|)");
    SHOW(edit_distance("kitten", "sitting"), "");
    SHOW(edit_distance("abc", "acb", 1, 1, 100), "置換が高いので削除と挿入");
    {
        // ABC185 E: A と B から要素を取り除いて同じ長さにし、取り除いた数 + 違う位置の数 を最小に → 削除 1・置換 1 の編集距離
        const std::vector<int> A = {1, 2, 1, 3}, B = {1, 3, 1};
        note("ABC185 E（Sequence Matching）: 「A から消す・B から消す（= A に挿入）・違う位置」がすべて費用 1 の編集距離");
        std::cout << "    入力例 1 → " << edit_distance(A, B) << "\n";
    }

    api("【4】 edit_script(a, b, ins, del, sub)", "vector<EditOp{type, i, j}>",
        "最小の費用の操作の列（前から）。'=' そのまま・'S' a[i] を b[j] に・'D' a[i] を消す・'I' b[j] を入れる", "O(|a||b|)");
    {
        const auto ops = edit_script("kitten", "sitting");
        std::string t;
        for (const auto& op : ops) t += op.type;
        std::cout << "  " << pad("edit_script(\"kitten\", \"sitting\")", 34) << " -> " << t << "  // k→s, e→i, g を挿入\n";
    }

    api("【5】 edit_distance_bounded(a, b, k) / within_one_edit(a, b)", "int / bool",
        "費用 1 の編集距離が k 以下ならその値、超えるなら -1（|a|, |b| が大きく k が小さいとき）/ 1 回以下の操作で一致するか", "O((|a| + |b|) k) / O(N)");
    SHOW(edit_distance_bounded("abcdef", "abxdefg", 3), "");
    SHOW(edit_distance_bounded("abcdef", "uvwxyz", 3), "超える → -1");
    SHOW(within_one_edit("atcoder", "atcorder"), "ABC324 C（Error Correction）の判定");
    note("ABC386 F（Operate K）は edit_distance_bounded(S, T, K) >= 0。");
}

}  // namespace lcs_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    lcs_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    lcs_test::benchmark();    // 速度計測
    lcs_demo::demo();         // 使い方デモ
    return 0;
}
