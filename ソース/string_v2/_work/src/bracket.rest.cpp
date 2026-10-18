// =============================================================================
//  ランダムテスト
// =============================================================================
namespace brk_test {

{COMMON_TEST}

static int kind_of(char c, const std::string& pairs) {  // 0 = 括弧以外、k + 1 = 開き、-(k + 1) = 閉じ
    const auto p = pairs.find(c);
    if (p == std::string::npos) return 0;
    return p % 2 == 0 ? (int)(p / 2) + 1 : -(int)(p / 2) - 1;
}
// 愚直：括弧だけを並べ、隣り合う同じ種類の開き・閉じを消すことを繰り返す
static std::vector<int> bf_match(const std::string& s, const std::string& pairs) {
    std::vector<int> m(s.size(), -1), idx;
    for (int i = 0; i < (int)s.size(); ++i)
        if (kind_of(s[i], pairs) != 0) idx.push_back(i);
    for (bool changed = true; changed;) {
        changed = false;
        for (std::size_t k = 0; k + 1 < idx.size(); ++k) {
            const int a = kind_of(s[idx[k]], pairs), b = kind_of(s[idx[k + 1]], pairs);
            if (a > 0 && b == -a) {
                m[idx[k]] = idx[k + 1], m[idx[k + 1]] = idx[k];
                idx.erase(idx.begin() + k, idx.begin() + k + 2);
                changed = true;
                break;
            }
        }
    }
    return m;
}
static bool bf_balanced(const std::string& s, const std::string& pairs) {
    const auto m = bf_match(s, pairs);
    for (int i = 0; i < (int)s.size(); ++i)
        if (kind_of(s[i], pairs) != 0 && m[i] < 0) return false;
    return true;
}
// 足す括弧の最小の数（区間 DP。何種類でも）
static int bf_min_insert(const std::string& s, const std::string& pairs) {
    const int n = (int)s.size();
    std::vector<std::vector<int>> f(n + 1, std::vector<int>(n + 1, 0));  // f[l][r] = s[l, r)
    for (int len = 1; len <= n; ++len)
        for (int l = 0; l + len <= n; ++l) {
            const int r = l + len, k = kind_of(s[l], pairs);
            if (k == 0) {
                f[l][r] = f[l + 1][r];
                continue;
            }
            int best = 1 + f[l + 1][r];
            if (k > 0)
                for (int m = l + 1; m < r; ++m)
                    if (kind_of(s[m], pairs) == -k) best = std::min(best, f[l + 1][m] + f[m + 1][r]);
            f[l][r] = best;
        }
    return f[0][n];
}
static std::string rnd_brackets(std::mt19937_64& rng, int n, const std::string& pairs, bool letters) {
    std::string s;
    const int t = (int)(rng() % 3);
    if (t == 0) {  // 一様
        for (int i = 0; i < n; ++i) s += (letters && rng() % 4 == 0) ? 'x' : pairs[rng() % pairs.size()];
    } else {  // 正しい括弧列を作ってから少し壊す
        std::string st;
        for (int i = 0; i < n; ++i) {
            if (letters && rng() % 5 == 0) {
                s += 'x';
            } else if (!st.empty() && (rng() % 2 || (int)st.size() >= n - i)) {
                s += pairs[pairs.find(st.back()) + 1], st.pop_back();
            } else {
                const char o = pairs[2 * (rng() % (pairs.size() / 2))];
                s += o, st += o;
            }
        }
        for (int k = (t == 1 ? 0 : 1 + (int)(rng() % 2)); k > 0 && n > 0; --k) s[rng() % n] = pairs[rng() % pairs.size()];
    }
    return s;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    const std::string kinds[] = {"()", "[]", "()[]", "()[]<>", ")("};
    for (int it = 0; it < 4000 * TEST_SCALE; ++it) {
        const std::string pairs = kinds[rng() % 5];
        const bool letters = rng() % 2;
        const int n = (int)(rng() % (it % 10 == 0 ? 40 : 14));
        const std::string s = rnd_brackets(rng, n, pairs, letters);
        const auto m = bf_match(s, pairs);
        BRK_CHECK(match_brackets(s, pairs) == m);
        BRK_CHECK(is_balanced(s, pairs) == bf_balanced(s, pairs));
        // 深さ・親
        const auto dep = bracket_depths(s, pairs), par = bracket_parents(s, pairs);
        for (int i = 0; i < n; ++i) {
            int d = 0, p = -1;
            for (int l = 0; l < n; ++l)
                if (m[l] > l && l < i && i < m[l]) ++d, p = std::max(p, l);
            BRK_CHECK(dep[i] == d && par[i] == p);
        }
        // 組を消す：繰り返し「開き + 括弧以外 + 閉じ」を消す（ABC307 D の操作）
        {
            std::string t = s;
            for (bool changed = true; changed;) {
                changed = false;
                for (int l = 0; l < (int)t.size() && !changed; ++l) {
                    const int k = kind_of(t[l], pairs);
                    if (k <= 0) continue;
                    int r = l + 1;
                    while (r < (int)t.size() && kind_of(t[r], pairs) == 0) ++r;
                    if (r < (int)t.size() && kind_of(t[r], pairs) == -k) t.erase(l, r - l + 1), changed = true;
                }
            }
            BRK_CHECK(erase_matched_pairs(s, pairs) == t);
        }
        // 部分文字列
        long long cnt = 0;
        StrSpan best{0, 0};
        const StaticBracket sb(s, pairs);
        for (int l = 0; l <= n; ++l)
            for (int r = l; r <= n; ++r) {
                const bool ok = bf_balanced(s.substr(l, r - l), pairs);
                BRK_CHECK(sb.is_balanced(l, r) == ok);
                if (ok && r > l) {
                    ++cnt;
                    if (r - l > best.len()) best = {l, r};
                }
            }
        for (int i = 0; i < n; ++i) BRK_CHECK(sb.partner(i) == m[i]);
        BRK_CHECK(count_balanced_substrings(s, pairs) == cnt);
        BRK_CHECK(longest_balanced_substring(s, pairs) == best);
        // 1 種類
        if (pairs.size() == 2) {
            const BracketBalance b = balance_of(s, pairs);
            long long cur = 0, mn = 0;
            for (char c : s) {
                const int k = kind_of(c, pairs);
                cur += k > 0 ? 1 : (k < 0 ? -1 : 0);
                mn = std::min(mn, cur);
            }
            BRK_CHECK(b.need_open == -mn && b.need_close == cur - mn && b.balanced() == bf_balanced(s, pairs));
            BRK_CHECK(min_insertions_to_balance(s, pairs) == bf_min_insert(s, pairs));
            const std::string fx = fix_brackets(s, pairs);
            BRK_CHECK(bf_balanced(fx, pairs) && (long long)fx.size() == n + min_insertions_to_balance(s, pairs) && is_subsequence(s, fx));
            // 連結
            const int cut = (int)(rng() % (n + 1));
            BRK_CHECK(balance_of(s.substr(0, cut), pairs) + balance_of(s.substr(cut), pairs) == b);
            BracketBalance acc;
            for (char c : s) acc += balance_of(c, pairs);
            BRK_CHECK(acc == b);
        } else {
            BRK_CHECK(bf_min_insert(s, pairs) >= 0);
        }
        ++cases;
    }
}

static void test_concat(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int k = 1 + (int)(rng() % 6);
        std::vector<std::string> vs(k);
        for (auto& v : vs) v = rnd_brackets(rng, 1 + (int)(rng() % 6), "()", rng() % 4 == 0);
        std::vector<int> perm(k);
        std::iota(perm.begin(), perm.end(), 0);
        bool exists = false;
        do {
            std::string t;
            for (int i : perm) t += vs[i];
            if (bf_balanced(t, "()")) exists = true;
        } while (!exists && std::next_permutation(perm.begin(), perm.end()));
        const BracketOrder r = bracket_concat_order(vs);
        BRK_CHECK(r.ok == exists);
        if (r.ok) {
            BRK_CHECK((int)r.order.size() == k);
            std::vector<int> seen(k, 0);
            std::string t;
            for (int i : r.order) {
                BRK_CHECK(i >= 0 && i < k && !seen[i]);
                seen[i] = 1, t += vs[i];
            }
            BRK_CHECK(bf_balanced(t, "()"));
        } else {
            BRK_CHECK(r.order.empty());
        }
        ++cases;
    }
}

static void test_segtree(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const std::string pairs = rng() % 2 ? "()" : "<>";
        const int n = 1 + (int)(rng() % 40);
        std::string s = rnd_brackets(rng, n, pairs, rng() % 2);
        BracketSegTree bt(s, pairs);
        for (int q = 0; q < 200; ++q) {
            const int op = (int)(rng() % 4);
            if (op == 0) {
                const int i = (int)(rng() % n);
                const char c = rng() % 5 == 0 ? 'x' : pairs[rng() % 2];
                s[i] = c, bt.set(i, c);
            } else if (op == 1) {
                const int i = (int)(rng() % n);
                if (s[i] == pairs[0]) s[i] = pairs[1];
                else if (s[i] == pairs[1]) s[i] = pairs[0];
                bt.flip(i);
            } else {
                int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
                if (l > r) std::swap(l, r);
                const std::string sub = s.substr(l, r - l);
                BRK_CHECK(bt.is_balanced(l, r) == bf_balanced(sub, pairs) && bt.balance(l, r) == balance_of(sub, pairs));
            }
            const int i = (int)(rng() % n);
            BRK_CHECK(bt.get(i) == s[i] && bt.all_balanced() == bf_balanced(s, pairs) && bt.to_string() == s);
        }
        ++cases;
    }
}

static void test_enumerate(int& cases) {
    for (const std::string pairs : {"()", ")(", "[]"}) {
        for (int n = 0; n <= 14; ++n) {
            std::vector<std::string> got, want;
            for_each_balanced(n, [&](const std::string& s) { got.push_back(s); }, pairs);
            for (int mask = 0; mask < (1 << n); ++mask) {
                std::string s(n, pairs[0]);
                for (int i = 0; i < n; ++i)
                    if (mask >> (n - 1 - i) & 1) s[i] = pairs[1];  // mask の昇順 = 開き < 閉じ の辞書順
                if (bf_balanced(s, pairs)) want.push_back(s);
            }
            BRK_CHECK(got == want);
            ++cases;
        }
    }
}

static void test_final_additions(std::mt19937_64& rng, int& cases) {  // count_balanced_completions
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 13);
        std::string s(n, ' ');
        for (char& c : s) c = "()?"[rng() % 3];
        std::vector<int> q;
        for (int i = 0; i < n; ++i)
            if (s[i] == '?') q.push_back(i);
        long long cnt = 0;
        for (int mask = 0; mask < (1 << q.size()); ++mask) {
            std::string t = s;
            for (std::size_t k = 0; k < q.size(); ++k) t[q[k]] = (mask >> k & 1) ? '(' : ')';
            cnt += is_balanced(t);
        }
        BRK_CHECK(count_balanced_completions<long long>(s) == cnt && count_balanced_completions<double>(s) == (double)cnt);
        std::string u = s;
        for (char& c : u) c = c == '(' ? 'A' : c == ')' ? 'B' : '*';
        BRK_CHECK(count_balanced_completions<long long>(u, '*', "AB") == cnt);
        ++cases;
    }
    // 途中だけ溢れても答えが収まれば止まらない：? が 70 個で答えは C(70, 35)/36（カタラン数）
    BRK_CHECK(count_balanced_completions<long long>(std::string(70, '?')) == 3116285494907301262LL);
}
static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_basic(rng, cases);
    test_concat(rng, cases);
    test_segtree(rng, cases);
    test_final_additions(rng, cases);
    test_enumerate(cases);
    std::cout << "[bracket] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[bracket] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 1000000;
    const std::string s = rnd_brackets(rng, N, "()[]", false);
    auto t0 = std::chrono::steady_clock::now();
    long long x = is_balanced(s, "()[]");
    x += match_brackets(s, "()[]")[N / 2];
    x += count_balanced_substrings(s, "()[]");
    const long long t_basic = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const StaticBracket sb(s, "()[]");
    for (int q = 0; q < N; ++q) {
        int l = (int)(rng() % (N + 1)), r = (int)(rng() % (N + 1));
        if (l > r) std::swap(l, r);
        x += sb.is_balanced(l, r);
    }
    const long long t_static = ms(t0);
    const int M = 500000;
    BracketSegTree bt(rnd_brackets(rng, M, "()", false));
    t0 = std::chrono::steady_clock::now();
    for (int q = 0; q < M; ++q) {
        if (q % 2) bt.flip((int)(rng() % M));
        else {
            int l = (int)(rng() % (M + 1)), r = (int)(rng() % (M + 1));
            if (l > r) std::swap(l, r);
            x += bt.is_balanced(l, r);
        }
    }
    const long long t_seg = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long cnt = 0;
    for_each_balanced(24, [&](const std::string& t) { cnt += t[1] == '('; });
    const long long t_enum = ms(t0);
    std::cout << "[bracket] benchmark: N = 10^6 で is_balanced + match_brackets + count_balanced_substrings " << t_basic << " ms, StaticBracket 作成 + 区間 10^6 回 "
              << t_static << " ms, BracketSegTree（N = Q = 5×10^5）" << t_seg << " ms, for_each_balanced(24)（208012 個）" << t_enum << " ms  (check " << (x + cnt) % 7 << ")\n";
}

}  // namespace brk_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace brk_demo {

{COMMON_DEMO}

static void demo() {
    banner("括弧列（bracket） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  is_balanced / match_brackets / bracket_depths / bracket_parents / erase_matched_pairs   判定・対応・深さ・入れ子の木・組を消す\n"
              << "  BracketBalance / balance_of / bracket_concat_order / min_insertions_to_balance / fix_brackets   足りない括弧の数（1 種類）\n"
              << "  count_balanced_substrings / longest_balanced_substring   正しい部分文字列の数・最長\n"
              << "  StaticBracket（区間の判定 O(1)）/ BracketSegTree（書き換えと区間の判定。1 種類）/ for_each_balanced（全列挙）\n"
              << "  count_balanced_completions<T>（? を決めて正しい括弧列にする方法の数）\n"
              << "  括弧の種類は pairs（既定 \"()\"。\"()[]<>\" なら 3 種類）。括弧以外の文字は無視する\n";

    api("【1】 is_balanced(s, pairs)", "bool", "正しい括弧列か（隣り合う同じ種類の開き・閉じを消して空になるか）", "O(N)");
    SHOW(is_balanced("(()())"), "");
    SHOW(is_balanced("())("), "");
    SHOW(is_balanced("a(b)c"), "括弧以外は無視");
    SHOW(is_balanced("([])<>()", "()[]<>"), "ABC394 D（Colorful Bracket Sequence）の入力例 1");
    SHOW(is_balanced("([<)]>", "()[]<>"), "入力例 2：種類が交差");

    api("【2】 match_brackets(s, pairs) / bracket_depths(s, pairs) / bracket_parents(s, pairs)", "vector<int>",
        "対応する括弧の位置（なければ -1）/ 各位置を囲む組の数 / 囲む一番内側の組の開き括弧の位置（なければ -1）", "O(N)");
    SHOW(match_brackets("(()))("), "最後の 2 つは対応しない");
    SHOW(bracket_depths("(a(b))c"), "");
    SHOW(bracket_parents("(a(b))c"), "");

    api("【3】 erase_matched_pairs(s, pairs)", "string", "対応する括弧の組を中身ごと消した残り", "O(N)");
    SHOW(erase_matched_pairs("a(b(d))c"), "ABC307 D（Mismatched Parentheses）の入力例 1");
    SHOW(erase_matched_pairs("a(b)("), "入力例 2");

    api("【4】 BracketBalance{need_open, need_close} / balance_of(s, pairs)（1 種類）", "BracketBalance",
        "対応しない閉じ括弧の数（前に足す開き）と開き括弧の数（後ろに足す閉じ）。a + b でつなげたもの。balanced() で正しいか", "O(N) / 連結 O(1)");
    SHOW(balance_of("))(()("), "{前に 2 個, 後ろに 2 個}");
    SHOW(balance_of(")(") + balance_of("()"), "");
    SHOW((balance_of("((") + balance_of(")") + balance_of(")")).balanced(), "");
    {
        // ABC428 C：1 文字足す・最後の 1 文字を消す を繰り返して、毎回 正しい括弧列か
        note("ABC428 C（Brackets Stack Query）: 「先頭からの BracketBalance」を積み、足すときは 1 文字ぶんを + して push、消すときは pop");
        std::vector<BracketBalance> st = {BracketBalance{}};
        std::string out;
        const std::vector<std::string> qs = {"1 (", "2", "1 (", "1 )", "2", "1 (", "1 )", "1 )"};
        for (const auto& q : qs) {
            if (q[0] == '1') st.push_back(st.back() + balance_of(q[2]));
            else st.pop_back();
            out += st.back().balanced() ? "Yes " : "No ";
        }
        note("入力例 1 → " + out);
    }

    api("【5】 bracket_concat_order(vs, pairs)（1 種類）", "BracketOrder{ok, order}",
        "全部の文字列を並べ替えてつなげ、正しい括弧列にする順（開きが多いものを必要な開きの少ない順 → 閉じが多いものを残す開きの多い順）", "O(Σ|vs| + K log K)");
    {
        const std::vector<std::string> vs = {")", "(()"};
        const BracketOrder r = bracket_concat_order(vs);
        code("bracket_concat_order({\")\", \"(()\"})   // ABC167 F（Bracket Sequencing）の入力例 1");
        std::cout << "    -> ok = " << to_s(r.ok) << ", order = " << to_s(r.order) << "\n";
        SHOW(bracket_concat_order(std::vector<std::string>{")(", "()"}).ok, "入力例 2");
    }

    api("【6】 min_insertions_to_balance(s, pairs) / fix_brackets(s, pairs)（1 種類）", "long long / string",
        "正しくするのに足す括弧の最小の数 / 前に開き・後ろに閉じを足して正しくしたもの", "O(N)");
    SHOW(min_insertions_to_balance("())("), "");
    SHOW(fix_brackets("())("), "");

    api("【7】 count_balanced_substrings(s, pairs) / longest_balanced_substring(s, pairs)", "long long / StrSpan",
        "正しい括弧列である空でない部分文字列の数 / 最長のもの（同じ長さなら左）", "O(N)");
    SHOW(count_balanced_substrings("()()(()"), "() が 3 つ・()() が 1 つ");
    SHOW(longest_balanced_substring(")()())"), "");

    api("【8】 StaticBracket sb(s, pairs)", "sb.is_balanced(l, r) -> bool / sb.partner(i) -> int",
        "書き換えのない文字列で、区間 [l, r) が正しい括弧列か O(1)（何種類でも）。L 文字目〜R 文字目は (L - 1, R)", "前処理 O(N)・質問 O(1)");
    {
        RUN(const StaticBracket sb("(())("));
        SHOW(sb.is_balanced(0, 4), "ABC223 F の入力例 1：1〜4 文字目");
        SHOW(sb.is_balanced(0, 2), "1〜2 文字目");
        SHOW(sb.partner(1), "");
    }

    api("【9】 BracketSegTree bt(s, pairs)（1 種類）", "set(i, c) / flip(i) / get(i) / is_balanced(l, r) / balance(l, r) / all_balanced()",
        "1 文字の書き換え・開き閉じの反転と、区間が正しい括弧列か（セグメント木）", "O(log N)（all_balanced は O(1)）");
    {
        RUN(BracketSegTree bt("(())("));
        SHOW(bt.is_balanced(0, 4), "ABC223 F の入力例 2：判定");
        note("1 4 の入れ替え: a = bt.get(0), b = bt.get(3) を読んで bt.set(0, b), bt.set(3, a)");
        const char a = bt.get(0), b = bt.get(3);
        bt.set(0, b), bt.set(3, a);
        SHOW(bt.to_string(), "");
        SHOW(bt.is_balanced(0, 4), "");
        RUN(bt.flip(4));
        SHOW(bt.balance(0, 5), "");
    }

    api("【10】 for_each_balanced(n, f, pairs)（1 種類）", "なし（f(const string&) を呼ぶ）", "長さ n の正しい括弧列を辞書順（開き < 閉じ）に全部。n が奇数なら呼ばない", "O(個数 × n)");
    {
        std::vector<std::string> all;
        for_each_balanced(6, [&](const std::string& t) { all.push_back(t); });
        note("典型90 002（Encyclopedia of Parentheses）の形。n = 6 → " + join(all));
    }

    head("問題での使い方：ABC350 F（Transpose）");
    note("組の中身を反転して大文字小文字を入れ替える → m = match_brackets(S)、d = bracket_depths(S) を作り、");
    note("i = 0 から向き +1 で歩く。括弧なら i = m[i] に飛んで向きを反転、文字なら d[i] が奇数なら入れ替えて出力。");
    {
        const std::string S = "((XYZ)n(X(y)Z))";
        const auto m = match_brackets(S);
        const auto d = bracket_depths(S);
        std::string out;
        for (int i = 0, dir = 1; 0 <= i && i < (int)S.size(); i += dir) {
            if (m[i] >= 0) i = m[i], dir = -dir;
            else out += d[i] % 2 ? swap_case(std::string(1, S[i]))[0] : S[i];
        }
        note("入力例 2 \"" + S + "\" → " + out);
    }

    api("【11】 count_balanced_completions<T>(s, wild = '?', pairs)（1 種類）", "T（mint・long long。答えが溢れたら止める）",
        "wild をそれぞれ開き・閉じに決めて正しい括弧列にする方法の数（dp[開き − 閉じ]）", "O(N^2)");
    SHOW(count_balanced_completions<long long>(std::string("(??") + "?(?"), "ABC312 D の入力例 1（2 つに分けて書いたのは ? の 3 連続の警告を避けるため）");
    SHOW(count_balanced_completions<long long>(std::string(10, '?')), "カタラン数 C_5");
    note("答えを mod で出すときは count_balanced_completions<mint>(S)。");

    head("【型】 途中で閉じが開きより多くならないか（区間・書き換えあり。ABC473 F）");
    note("「正しい括弧列」ではなく「どの接頭辞でも 閉じ ≤ 開き（最後に開きが余ってよい）」は BracketSegTree の balance(l, r).need_open == 0。");
    note("括弧でない 2 文字も pairs に渡せる：BracketSegTree bt(S, \"AB\")（A が開き、B が閉じ）。");
    {
        BracketSegTree bt("AABAB", "AB");
        std::cout << "    S = AABAB：balance(0, 5).need_open = " << bt.balance(0, 5).need_open << "（途中で B が多くならない）、balance(2, 5).need_open = "
                  << bt.balance(2, 5).need_open << "（\"BAB\" は先頭の B が多い）\n";
    }
}

}  // namespace brk_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    brk_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    brk_test::benchmark();    // 速度計測
    brk_demo::demo();         // 使い方デモ
    return 0;
}
