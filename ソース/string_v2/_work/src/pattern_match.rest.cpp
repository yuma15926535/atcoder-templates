// =============================================================================
//  ランダムテスト
// =============================================================================
namespace pm_test {

{COMMON_TEST}

static std::vector<int> bf_find(const std::string& t, const std::string& p) {
    std::vector<int> r;
    for (int i = 0; i + (int)p.size() <= (int)t.size(); ++i)
        if (t.compare(i, p.size(), p) == 0) r.push_back(i);
    return r;
}
static int bf_lcp(const std::string& a, int i, const std::string& b, int j) {
    int k = 0;
    while (i + k < (int)a.size() && j + k < (int)b.size() && a[i + k] == b[j + k]) ++k;
    return k;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 8000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string t = rnd_str_shaped(rng, (int)(rng() % 25), sg), p = rnd_str_shaped(rng, (int)(rng() % 5), sg);
        const int n = (int)t.size(), m = (int)p.size();
        // z / π
        const auto z = z_array(t);
        for (int i = 0; i < n; ++i) PM_CHECK(z[i] == bf_lcp(t, 0, t, i));
        const auto pi = prefix_function(t);
        for (int i = 0; i < n; ++i) {
            int b = 0;
            for (int k = 1; k <= i; ++k)
                if (t.compare(0, k, t, i + 1 - k, k) == 0) b = k;
            PM_CHECK(pi[i] == b);
        }
        // 一致の長さ・出現
        const auto ml = match_lengths(t, p);
        for (int i = 0; i < n; ++i) PM_CHECK(ml[i] == bf_lcp(t, i, p, 0));
        std::vector<int> ex = bf_find(t, p);
        if (m == 0) {
            ex.clear();
            for (int i = 0; i <= n; ++i) ex.push_back(i);
        }
        PM_CHECK(find_all(t, p) == ex && count_occurrences(t, p) == (long long)ex.size());
        PM_CHECK(find_last(t, p) == (ex.empty() ? -1 : ex.back()));
        const int from = (int)(rng() % (n + 1));
        int ff = -1;
        for (int x : ex)
            if (x >= from) { ff = x; break; }
        PM_CHECK(find_first(t, p, from) == ff && find_first(t, p) == (ex.empty() ? -1 : ex[0]));
        if (m > 0) {
            long long c = 0;
            std::string rep;
            for (int i = 0; i < n;) {
                if (t.compare(i, m, p) == 0 && i + m <= n) ++c, rep += "XY", i += m;
                else rep += t[i++];
            }
            PM_CHECK(count_non_overlapping(t, p) == c && replace_all(t, p, "XY") == rep);
        }
        // 周期・境界
        std::vector<int> bd, pd;
        for (int k = n - 1; k >= 1; --k)
            if (t.compare(0, k, t, n - k, k) == 0) bd.push_back(k);
        for (int q = 1; q <= n; ++q) {
            bool ok = true;
            for (int i = 0; i + q < n; ++i) ok &= t[i] == t[i + q];
            if (ok) pd.push_back(q);
        }
        PM_CHECK(all_borders(t) == bd && all_periods(t) == pd && min_period(t) == (n ? pd[0] : 0));
        int root = n;
        for (int q : pd)
            if (n % q == 0) { root = q; break; }
        PM_CHECK(repetition_root(t) == root);
        const auto pmp = prefix_min_periods(t);
        PM_CHECK((int)pmp.size() == n + 1 && pmp[0] == 0);
        for (int i = 1; i <= n; i += 3) PM_CHECK(pmp[i] == min_period(t.substr(0, i)));
        // 回転
        std::string u = t;
        if (n && rng() % 2) std::rotate(u.begin(), u.begin() + rng() % n, u.end());
        else if (n && rng() % 2) u[rng() % n] = (char)('a' + rng() % 3);
        std::vector<int> ro;
        for (int k = 0; k < n; ++k)
            if (t.substr(k) + t.substr(0, k) == u) ro.push_back(k);
        if (n == 0) ro = {0};
        PM_CHECK(rotation_offsets(t, u) == ro && is_rotation(t, u) == !ro.empty());
        // 無限に繰り返した中の出現（pat は長くてもよい）
        const std::string q = rnd_str_shaped(rng, (int)(rng() % 12), sg);
        if (n > 0) {
            std::string big;
            while (big.size() < (size_t)n + q.size()) big += t;
            std::vector<int> cy;
            for (int i = 0; i < n; ++i)
                if (big.compare(i, q.size(), q) == 0) cy.push_back(i);
            PM_CHECK(find_all_cyclic(t, q) == cy);
        }
        // 重なり
        int ov = 0;
        for (int k = 1; k <= std::min(n, (int)q.size()); ++k)
            if (t.compare(n - k, k, q, 0, k) == 0) ov = k;
        PM_CHECK(max_overlap(t, q) == ov);
        ++cases;
    }
}

static void test_automaton(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string p = rnd_str_shaped(rng, (int)(rng() % 6), sg), t = rnd_str(rng, (int)(rng() % 30), sg);
        const KmpAutomaton KA(p);
        PM_CHECK(KA.size() == (int)p.size() + 1);
        int st = 0;
        long long acc = 0;
        for (int i = 0; i < (int)t.size(); ++i) {
            st = KA.next(st, t[i]);
            acc += KA.accepting(st);
            // 状態 = 読んだ文字列の接尾辞で p の接頭辞であるものの最長（p の全体まで）
            int best = 0;
            for (int k = 1; k <= std::min((int)p.size(), i + 1); ++k)
                if (t.compare(i + 1 - k, k, p, 0, k) == 0) best = k;
            PM_CHECK(st == best);
        }
        if (!p.empty()) PM_CHECK(acc == count_occurrences(t, p));
        ++cases;
    }
    // 括弧の文字で
    const KmpAutomaton KB("()", Alphabet::of("()"));
    PM_CHECK(KB.next(KB.next(0, '('), ')') == 2 && KB.accepting(2));
}

static void test_values(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 20), m = 1 + (int)(rng() % 3);
        std::vector<long long> t(n), p(m);
        const long long big = 1000000000000000000LL;
        for (auto& x : t) x = (rng() % 2 ? big : -big) * (long long)(rng() % 2);
        for (auto& x : p) x = (rng() % 2 ? big : -big) * (long long)(rng() % 2);
        std::vector<int> ex;
        for (int i = 0; i + m <= n; ++i)
            if (std::equal(p.begin(), p.end(), t.begin() + i)) ex.push_back(i);
        PM_CHECK(find_all(t, p) == ex);
        ++cases;
    }
    // string と vector<int> を混ぜる・リテラル
    PM_CHECK(find_all(std::vector<int>{'a', 'b', 'a', 'b'}, "ab") == std::vector<int>({0, 2}));
    PM_CHECK(count_occurrences("aaaa", "aa") == 3 && count_non_overlapping("aaaa", "aa") == 2 && replace_all("aaaa", "aa", "b") == "bb");
    PM_CHECK(replace_all(std::vector<int>{1, 2, 1, 2}, std::vector<int>{1, 2}, std::vector<int>{9}) == std::vector<int>({9, 9}));
}

// 最終照合で足したもの：置き換えの繰り返し・区間の出現
static std::string bf_rewrite(std::string s, const std::string& f, const std::string& t) {
    while (true) {
        const std::size_t i = s.find(f);
        if (i == std::string::npos) return s;
        s.replace(i, f.size(), t);
    }
}
static void test_final_additions(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 30), sg);
        const std::string f = rnd_str_shaped(rng, 1 + (int)(rng() % 3), sg);
        std::string t = rnd_str(rng, (int)(rng() % (f.size() + 1)), sg);  // |t| ≤ |f|
        if (t.find(f) != std::string::npos) t = "";
        if (t.size() == f.size() && rng() % 2) t = std::string(f.rbegin(), f.rend());  // 同じ長さ（ab → ba など）
        if (t.find(f) != std::string::npos) t = "";
        PM_CHECK(rewrite_until_stable(s, f, t) == bf_rewrite(s, f, t));
        PM_CHECK(erase_until_stable(s, f) == bf_rewrite(s, f, ""));
        // 区間の出現
        const std::string p = rnd_str_shaped(rng, (int)(rng() % 3), sg);
        const Occurrences occ(s, p);
        const int n = (int)s.size();
        PM_CHECK(occ.positions() == find_all(s, p) && occ.size() == (int)find_all(s, p).size());
        for (int q = 0; q < 20; ++q) {
            int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
            if (l > r) std::swap(l, r);
            long long c = 0;
            for (int i = l; i + (int)p.size() <= r; ++i) c += s.compare(i, p.size(), p) == 0;
            PM_CHECK(occ.count(l, r) == c && occ.exists(l, r) == (c > 0));
            int nx = -1, pv = -1;
            for (int x : find_all(s, p)) {
                if (x >= l && nx < 0) nx = x;
                if (x < l) pv = x;
            }
            PM_CHECK(occ.next(l) == nx && occ.prev(l) == pv);
        }
        ++cases;
    }
    // 数の列・1 文字の pat・伸びる置き換え
    PM_CHECK(rewrite_until_stable(std::vector<int>{1, 2, 1, 2, 3}, std::vector<int>{1, 2}, std::vector<int>{2}) == std::vector<int>({2, 2, 3}));
    PM_CHECK(erase_until_stable(std::string("AABCBCC"), "ABC") == "C" && rewrite_until_stable(std::string("WWWWA"), "WA", "AC") == "ACCCC");
    PM_CHECK(find_all("abca", 'a') == std::vector<int>({0, 3}) && replace_all("a,b", ',', '.') == "a.b");
    PM_CHECK(rewrite_until_stable(std::string("xyz"), "a", "aa") == "xyz");  // 置き換えが起きなければ to が from を含んでもよい
    ++cases;
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_basic(rng, cases);
    test_automaton(rng, cases);
    test_values(rng, cases);
    test_final_additions(rng, cases);
    std::cout << "[pattern_match] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[pattern_match] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 10000000;
    const std::string s = rnd_str_shaped(rng, n, 2), a(n, 'a');
    auto t0 = std::chrono::steady_clock::now();
    const auto z = z_array(s);
    const auto z2 = z_array(a);
    const long long t_z = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto pi = prefix_function(s);
    const long long t_pi = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto occ = find_all(s, s.substr(1234567, 20));
    const long long t_find = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto ml = match_lengths(s, s.substr(0, 100000));
    const long long t_ml = ms(t0);
    const std::string p = rnd_str(rng, 100000, 26);
    t0 = std::chrono::steady_clock::now();
    const KmpAutomaton KA(p);
    int st = 0;
    for (int i = 0; i < n; ++i) st = KA.next_id(st, (int)(rng() % 26));
    const long long t_ka = ms(t0);
    std::string nest = std::string(500000, '(') + "xx" + std::string(500000, ')'), wa = std::string(1000000, 'W') + "A";
    t0 = std::chrono::steady_clock::now();
    const auto r1 = rewrite_until_stable(nest, "(xx)", "xx"), r2 = rewrite_until_stable(wa, "WA", "AC");
    const long long t_rw = ms(t0);
    std::cout << "[pattern_match] benchmark: z_array（N = 10^7 × 2）" << t_z << " ms, prefix_function（N = 10^7）" << t_pi << " ms, find_all（N = 10^7）"
              << t_find << " ms, match_lengths（N = 10^7）" << t_ml << " ms, KmpAutomaton（M = 10^5, σ = 26）+ 10^7 回の遷移 " << t_ka
              << " ms, rewrite_until_stable（入れ子 10^6・WA→AC 10^6）" << t_rw << " ms  (check "
              << (z[5] + z2[3] + pi[7] + (long long)occ.size() + ml[9] + st + (long long)r1.size() + (long long)r2.size()) % 7 << ")\n";
}

}  // namespace pm_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace pm_demo {

{COMMON_DEMO}

static void demo() {
    banner("パターン検索（pattern_match） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  z_array / prefix_function / match_lengths              Z 配列 / KMP の境界 / 各位置で pat と何文字一致するか\n"
              << "  find_all / count_occurrences / find_first / find_last  出現位置・回数（重なりも）\n"
              << "  count_non_overlapping / replace_all                    重ならない出現の数 / 全部置き換え（1 回）\n"
              << "  rewrite_until_stable / erase_until_stable              置き換え・消すを現れなくなるまで繰り返す\n"
              << "  Occurrences                                            区間 [l, r) に収まる出現の数・次の出現\n"
              << "  all_borders / min_period / all_periods / repetition_root / prefix_min_periods   境界・周期\n"
              << "  rotation_offsets / is_rotation / find_all_cyclic        回転・無限に繰り返した中の出現\n"
              << "  max_overlap                                             a の後ろと b の前の重なり\n"
              << "  KmpAutomaton                                            パターンを含まない文字列の DP\n";

    api("【1】 z_array(s) / prefix_function(s)", "vector<int>",
        "z[i] = s と s[i..) の最長共通接頭辞（z[0] = |s|。ACL の z_algorithm と同じ）/ π[i] = s[0, i] の最長の真の境界", "O(N)");
    SHOW(z_array("aabaab"), "");
    SHOW(prefix_function("aabaab"), "");

    api("【2】 match_lengths(text, pat)", "vector<int>（長さ |text|）",
        "各 i で text[i..) と pat の最長共通接頭辞の長さ（Z を text と pat に）", "O(N + M)");
    SHOW(match_lengths("abcab", "ab"), "");
    {
        // ABC257 G: T を S の接頭辞をつなげて作る最小の個数（位置 i から i + 1〜i + m[i] へ進める → 貪欲）
        const std::string S = "aba", T = "ababaab";
        const auto m = match_lengths(T, S);
        int cnt = 0, cur_end = 0, far = 0;
        bool ok = true;
        for (int i = 0; i < (int)T.size() && ok; ++i) {
            far = std::max(far, i + m[i]);
            if (i == cur_end) {
                if (far == i) ok = false;
                ++cnt, cur_end = far;
            }
        }
        note("ABC257 G（Prefix Concatenation）: m = match_lengths(T, S) で「i から i + m[i] まで進める」→ 貪欲に遠くへ");
        std::cout << "    S = aba, T = ababaab → " << (ok ? cnt : -1) << "\n";
    }

    api("【3】 find_all / count_occurrences / find_first / find_last / count_non_overlapping / replace_all", "vector<int> / long long / int / string",
        "pat の出現（重なりも）。重ならない数・全部置き換えは左から貪欲に。数の列にも（ABC150 F のように差や XOR の列にして）", "O(N + M)");
    SHOW(find_all("abababa", "aba"), "重なりも");
    SHOW(count_occurrences("abababa", "aba"), "");
    SHOW(count_non_overlapping("abababa", "aba"), "左から重ならないように");
    SHOW(find_first("abababa", "aba", 1), "1 以上で最初");
    SHOW(find_last("abababa", "aba"), "");
    SHOW(find_first("abc", "x"), "ない → -1");
    SHOW(replace_all("abababa", "aba", "X"), "");
    SHOW(find_all(std::vector<int>{3, 1, 4, 1, 4}, std::vector<int>{1, 4}), "数の列");
    note("ABC343 G のように「a が b に含まれるか」は find_first(b, a) >= 0。1 文字の pat も渡せる（find_all(s, 'a')）。");

    api("【3b】 rewrite_until_stable(text, from, to) / erase_until_stable(text, pat)", "string（vector なら vector）",
        "いちばん左の from を to にすることを、from が現れなくなるまで繰り返した結果（replace_all をループで回すと O(N^2) になる形も速い）",
        "O(N + 置き換えの回数 × (|from| + |to|))");
    SHOW(rewrite_until_stable("WWWWA", "WA", "AC"), "ABC394 C（Debug）：WA → AC をなくなるまで");
    SHOW(erase_until_stable("BAABCBCCABCAC", "ABC"), "ABC328 D（Take ABC）の入力例 1：ABC を消し続ける");
    SHOW(rewrite_until_stable("((xx))(xx)", "(xx)", "xx"), "ABC454 D：(xx) → xx を繰り返した形どうしを比べる");
    note("|to| < |from| なら必ず終わる。to が from を含むと終わらないので止まる。");

    api("【3c】 Occurrences occ(text, pat)：count(l, r) / exists(l, r) / next(i) / prev(i) / positions()", "構造体",
        "区間 [l, r) に収まる（開始が l 以上・終わりが r 以下）出現の数・あるか / i 以上で最初・i より前で最後の出現の開始位置", "前処理 O(N + M)・O(log N)");
    RUN(const Occurrences occ("abcabcab", "ab"));
    SHOW(occ.count(1, 8), "[1, 8) に収まる \"ab\"（3 と 6）");
    SHOW(occ.exists(1, 4), "\"bca\" には \"ab\" がない");
    SHOW(occ.next(1), "1 以上で最初の出現");
    note("ABC477 C（区間に T が現れるか）。部分文字列 s[l, r) が全体に何回現れるかは 03_suffix_array の count_occurrences(l, r)（意味が違う）。");

    api("【4】 all_borders / min_period / all_periods / repetition_root / prefix_min_periods", "vector<int> / int",
        "境界（接頭辞 = 接尾辞）・周期（s[i] = s[i + p]）・s = t^k となる最短の t・各接頭辞の最小の周期", "O(N)");
    SHOW(all_borders("abaabaab"), "大きい順");
    SHOW(min_period("abaabaab"), "割り切れなくてよい");
    SHOW(all_periods("abaabaab"), "");
    SHOW(repetition_root("abcabcabc"), "\"abc\" × 3");
    SHOW(repetition_root("abcabca"), "繰り返しでない → 長さそのもの");
    SHOW(prefix_min_periods("abab"), "r[i] = s[0, i) の最小の周期");

    api("【5】 rotation_offsets(a, b) / is_rotation(a, b) / find_all_cyclic(text, pat)", "vector<int> / bool / vector<int>",
        "a を左に k 回まわす（a[k..] + a[..k]）と b になる k / text を無限に繰り返した中で pat が始まる位置（|pat| > |text| でもよい）", "O(N + M)");
    SHOW(rotation_offsets("abcab", "cabab"), "");
    SHOW(rotation_offsets("abab", "baba"), "周期があると複数");
    SHOW(is_rotation("abc", "acb"), "");
    SHOW(find_all_cyclic("ab", "babab"), "\"ababab…\" の 1 文字目から");
    note("ABC135 F（Strings of Eternity）: find_all_cyclic(s, t) で i → i + |t| (mod |s|) の辺を作り、");
    note("graph_grid_v2 の 06_functional_graph で閉路と最長の道を調べる。");

    api("【6】 max_overlap(a, b)", "int", "a の接尾辞と b の接頭辞が一致する最大の長さ（a と b を重ねて短くつなぐ）", "O(N + M)");
    SHOW(max_overlap("abcab", "cabd"), "\"cab\" が重なる");
    SHOW(max_overlap("abc", "xyz"), "");
    note("ABC343 G（Compress Strings）: 含まれる文字列を除き、max_overlap を重みにして bit DP。");

    api("【7】 KmpAutomaton KA(pat, A) / KA.next(state, c) / KA.accepting(state)", "構造体",
        "状態 = 読んだ文字列の末尾が pat の先頭何文字と一致しているか。pat を含まない文字列の数え上げ DP に", "構築 O(Mσ)・遷移 O(1)");
    {
        // 長さ 6 の a, b からなる文字列で "aba" を含まないものの数（DP）
        const KmpAutomaton KA("aba", Alphabet::of("ab"));
        std::vector<long long> dp(KA.size(), 0);
        dp[0] = 1;
        for (int step = 0; step < 6; ++step) {
            std::vector<long long> nd(KA.size(), 0);
            for (int st = 0; st < KA.size(); ++st)
                for (char c : {'a', 'b'}) {
                    const int to = KA.next(st, c);
                    if (!KA.accepting(to)) nd[to] += dp[st];
                }
            dp = nd;
        }
        note("長さ 6 の {a, b} の文字列で \"aba\" を含まないもの: dp[状態] を 6 回遷移（accepting に入る遷移を捨てる）");
        std::cout << "    → " << std::accumulate(dp.begin(), dp.end(), 0LL) << "（全 64 個から）\n";
        note("長さが 10^18 なら、遷移を行列にして math_v2 の 05_matrix（Matrix<mint>::pow）。パターンが多いなら 02_aho_corasick。");
    }
}

}  // namespace pm_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    pm_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pm_test::benchmark();    // 速度計測
    pm_demo::demo();         // 使い方デモ
    return 0;
}
