// =============================================================================
//  ランダムテスト
// =============================================================================
namespace pal_test {

{COMMON_TEST}

static bool bf_pal(const std::string& s, int l, int r) {
    for (int i = l, j = r - 1; i < j; ++i, --j)
        if (s[i] != s[j]) return false;
    return true;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 16), sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, n, sg);
        const Palindromes P(s);
        PAL_CHECK(P.size() == n);
        long long cnt = 0;
        StrSpan best{0, 0};
        for (int l = 0; l <= n; ++l)
            for (int r = l; r <= n; ++r) {
                const bool p = bf_pal(s, l, r);
                PAL_CHECK(P.is_palindrome(l, r) == p);
                if (p && r > l) {
                    ++cnt;
                    if (r - l > best.len()) best = {l, r};
                }
            }
        PAL_CHECK(P.count_all() == cnt && P.longest() == best);
        const auto L = center_palindrome_lengths(s);
        PAL_CHECK((int)L.size() == std::max(0, 2 * n - 1));
        for (int c = 0; c < 2 * n - 1; ++c) {
            int mx = 0;
            for (int l = 0; l <= n; ++l) {
                const int r = c + 1 - l;
                if (r >= l && r <= n && r - l >= mx && bf_pal(s, l, r) && (c % 2 == 0 ? r > l : true)) mx = std::max(mx, r - l);
            }
            PAL_CHECK(L[c] == mx && P.length_at_center(c) == mx && P.span_at_center(c).len() == mx);
        }
        for (int i = 0; i < n; ++i) PAL_CHECK(2 * P.radius_odd(i) - 1 == L[2 * i]);
        for (int i = 0; i <= n; ++i) PAL_CHECK(2 * P.radius_even(i) == (i == 0 || i == n ? 0 : L[2 * i - 1]));
        int pre = 0, suf = 0;
        for (int k = 1; k <= n; ++k) {
            if (bf_pal(s, 0, k)) pre = k;
            if (bf_pal(s, n - k, n)) suf = k;
        }
        PAL_CHECK(P.longest_prefix_palindrome() == pre && P.longest_suffix_palindrome() == suf);
        const auto st = P.longest_starting_at(), en = P.longest_ending_at();
        for (int i = 0; i < n; ++i) {
            int a = 0, b = 0;
            for (int k = 1; i + k <= n; ++k)
                if (bf_pal(s, i, i + k)) a = k;
            for (int k = 1; k <= i + 1; ++k)
                if (bf_pal(s, i + 1 - k, i + 1)) b = k;
            PAL_CHECK(st[i] == a && en[i] == b);
        }
        // 回文にする
        const std::string wp = min_palindrome_with_prefix(s), ws = min_palindrome_with_suffix(s);
        PAL_CHECK(is_palindrome(wp) && wp.compare(0, n, s) == 0 && (int)wp.size() == 2 * n - suf);
        PAL_CHECK(is_palindrome(ws) && ws.compare(ws.size() - n, n, s) == 0 && (int)ws.size() == 2 * n - pre);
        int ch = 0;
        for (int i = 0; i < n / 2; ++i) ch += s[i] != s[n - 1 - i];
        PAL_CHECK(min_changes_to_palindrome(s) == ch);
        std::map<char, int> mc;
        for (char x : s) ++mc[x];
        int odd = 0;
        for (auto& [k, v] : mc) odd += v & 1;
        PAL_CHECK(can_rearrange_to_palindrome(s) == (odd <= 1));
        // 回文分割（O(N^3) の DP）
        std::vector<int> dp(n + 1, 1 << 30);
        dp[0] = 0;
        for (int r = 1; r <= n; ++r)
            for (int l = 0; l < r; ++l)
                if (bf_pal(s, l, r)) dp[r] = std::min(dp[r], dp[l] + 1);
        const auto parts = palindrome_partition_min(s);
        PAL_CHECK((int)parts.size() == dp[n]);
        int pos = 0;
        for (auto sp : parts) {
            PAL_CHECK(sp.l == pos && sp.r > sp.l && bf_pal(s, sp.l, sp.r));
            pos = sp.r;
        }
        PAL_CHECK(pos == n);
        ++cases;
    }
}

static void test_subseq(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 13), sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, n, sg);
        int best = 0;
        long long cnt = 0;
        for (int mask = 1; mask < (1 << n); ++mask) {
            std::string t;
            for (int i = 0; i < n; ++i)
                if (mask >> i & 1) t += s[i];
            if (is_palindrome(t)) ++cnt, best = std::max(best, (int)t.size());
        }
        PAL_CHECK(longest_palindromic_subsequence_length(s) == best);
        const std::string lps = longest_palindromic_subsequence(s);
        PAL_CHECK((int)lps.size() == best && is_palindrome(lps) && is_subsequence(lps, s));
        PAL_CHECK(count_palindromic_subsequences<long long>(s) == cnt && count_palindromic_subsequences<double>(s) == (double)cnt);
        // 数の列でも
        const std::vector<long long> v(s.begin(), s.end());
        PAL_CHECK(longest_palindromic_subsequence_length(v) == best && Palindromes(v).count_all() == Palindromes(s).count_all());
        ++cases;
    }
}


// 長さ n の全部の文字列（A の番号の列）を辞書順で
static std::vector<std::vector<int>> all_id_strings(int n, int sg) {
    std::vector<std::vector<int>> r(1);
    for (int k = 0; k < n; ++k) {
        std::vector<std::vector<int>> nr;
        for (const auto& v : r)
            for (int c = 0; c < sg; ++c) {
                auto w = v;
                w.push_back(c);
                nr.push_back(w);
            }
        r = nr;
    }
    return r;
}
static void test_final_additions(std::mt19937_64& rng, int& cases) {  // count_palindromes_leq
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3), n = (int)(rng() % 7);
        std::string letters = std::string("abc").substr(0, sg);
        std::shuffle(letters.begin(), letters.end(), rng);
        const Alphabet A = Alphabet::of(letters);  // 独自の順も
        std::string s(n, ' ');
        for (char& c : s) c = letters[rng() % sg];
        const std::vector<int> sid = to_ids(s, A);
        long long cnt = 0;
        for (const auto& v : all_id_strings(n, sg)) cnt += std::equal(v.begin(), v.end(), v.rbegin()) && v <= sid;
        PAL_CHECK(count_palindromes_leq<long long>(s, A) == cnt && count_palindromes_leq<double>(s, A) == (double)cnt);
        ++cases;
    }
    PAL_CHECK(count_palindromes_leq<long long>(std::string(26, 'z')) == 2481152873203736576LL);  // 26^13
    PAL_CHECK(count_palindromes_leq<long long>("AXA", Alphabet::upper()) == 24 && count_palindromes_leq<long long>("ABCZAZ", Alphabet::upper()) == 29);  // ABC242 E
}
static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_basic(rng, cases);
    test_subseq(rng, cases);
    test_final_additions(rng, cases);
    // 答えが long long に収まれば、途中の和が大きくても止まらない
    {
        const std::string s = "b" + std::string(15, 'a') + "b" + std::string(20, 'a') + "b" + std::string(7, 'a') + "b" + std::string(8, 'a') + "b" + std::string(12, 'a');
        PAL_CHECK(count_palindromic_subsequences<long long>(s) == 6172198434589173569LL);
        ++cases;
    }
    std::cout << "[palindrome] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[palindrome] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 10000000;
    const std::string s = rnd_str_shaped(rng, n, 2), a(n, 'a');
    auto t0 = std::chrono::steady_clock::now();
    const Palindromes P(s), PA(a);
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = 0;
    for (int i = 0; i < 10000000; ++i) {
        int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
        if (l > r) std::swap(l, r);
        x += P.is_palindrome(l, r);
    }
    const long long t_q = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += P.count_all() + PA.longest().len() + P.longest_starting_at()[3];
    const long long t_misc = ms(t0);
    const std::string m = rnd_str(rng, 5000, 3);
    t0 = std::chrono::steady_clock::now();
    x += (long long)longest_palindromic_subsequence(m).size() + (long long)palindrome_partition_min(m).size();
    const long long t_n2 = ms(t0);
    std::cout << "[palindrome] benchmark: Palindromes 構築（N = 10^7 × 2）" << t_build << " ms, is_palindrome × 10^7 " << t_q
              << " ms, count_all + longest + longest_starting_at（N = 10^7）" << t_misc << " ms, 最長回文部分列 + 回文分割（N = 5000）" << t_n2
              << " ms  (check " << x % 7 << ")\n";
}

}  // namespace pal_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace pal_demo {

{COMMON_DEMO}

static void demo() {
    banner("回文（palindrome） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Palindromes P(s)：is_palindrome(l, r) / radius_odd / radius_even / length_at_center / span_at_center / longest / count_all\n"
              << "                    longest_prefix_palindrome / longest_suffix_palindrome / longest_starting_at / longest_ending_at\n"
              << "  center_palindrome_lengths(s)                         Manacher の生の値（中心 2N-1 個）\n"
              << "  min_palindrome_with_prefix / min_palindrome_with_suffix   s で始まる / 終わる最短の回文\n"
              << "  min_changes_to_palindrome / can_rearrange_to_palindrome  書き換え・並べ替えて回文に\n"
              << "  longest_palindromic_subsequence(_length) / palindrome_partition_min / count_palindromic_subsequences<T>\n"
              << "  count_palindromes_leq<T>(s, A)                       同じ長さで s 以下の回文の数\n";

    api("【1】 Palindromes P(s) / P.is_palindrome(l, r)", "bool", "s[l, r) が回文か（Manacher で前処理）。何度でも O(1)", "構築 O(N)・判定 O(1)");
    RUN(const std::string S = "abacabad");
    RUN(const Palindromes P(S));
    SHOW(P.is_palindrome(0, 7), "\"abacaba\"");
    SHOW(P.is_palindrome(1, 4), "\"bac\"");
    note("問題の「L 文字目から R 文字目」は P.is_palindrome(L - 1, R)。全体だけなら is_palindrome(s)（00_str_base）。");

    api("【2】 P.radius_odd(i) / P.radius_even(i) / P.length_at_center(c) / center_palindrome_lengths(s)", "int / vector<int>",
        "s[i] が中心の最長の回文の半径 r（長さ 2r - 1）/ s[i-1] と s[i] の間が中心の半径 r（長さ 2r）/ 中心 c = 0〜2N-2 ごとの長さ", "O(1) / O(N)");
    SHOW(P.radius_odd(3), "\"abacaba\" は半径 4");
    SHOW(P.radius_even(3), "\"a|c\" の間は回文なし");
    SHOW(center_palindrome_lengths("abba"), "c = 0〜6（奇数の c が文字の間）");

    api("【3】 P.longest() / P.count_all() / P.longest_prefix_palindrome() / P.longest_suffix_palindrome()", "StrSpan / long long / int",
        "最長の回文（同じ長さなら左）・回文である部分文字列の数（位置が違えば別）・最長の回文の接頭辞 / 接尾辞", "O(N)");
    SHOW(substr_of(S, P.longest()), "");
    SHOW(P.count_all(), "");
    SHOW(P.longest_prefix_palindrome(), "\"abacaba\"");
    SHOW(P.longest_suffix_palindrome(), "\"d\"");

    api("【4】 P.longest_starting_at() / P.longest_ending_at()", "vector<int>",
        "各位置から始まる / で終わる最長の回文の長さ（回文で区切る DP の前処理などに）", "O(N)");
    SHOW(P.longest_starting_at(), "");
    SHOW(P.longest_ending_at(), "");

    api("【5】 min_palindrome_with_prefix(s) / min_palindrome_with_suffix(s)", "string",
        "s で始まる最短の回文（後ろに足す）/ s で終わる最短の回文（前に足す）", "O(N)");
    SHOW(min_palindrome_with_prefix("ABC"), "ABC398 F（ABCBA）の入力例 1");
    SHOW(min_palindrome_with_prefix("abcba"), "もう回文");
    SHOW(min_palindrome_with_suffix("aab"), "前に足す");

    api("【6】 min_changes_to_palindrome(s) / can_rearrange_to_palindrome(s)", "int / bool",
        "1 文字ずつ書き換えて回文にする最小の回数 / 並べ替えて回文にできるか（奇数個の文字が 1 種類以下）", "O(N)");
    SHOW(min_changes_to_palindrome("abcd"), "a ↔ d、b ↔ c");
    SHOW(can_rearrange_to_palindrome("aabbc"), "abcba");
    SHOW(can_rearrange_to_palindrome("aabc"), "");

    api("【7】 longest_palindromic_subsequence(s) / _length(s) / count_palindromic_subsequences<T>(s)", "string / int / T",
        "最長回文部分列（飛び飛び）を 1 つ・その長さ・回文である部分列の数（位置が違えば別。T は mint や long long）", "O(N^2)");
    SHOW(longest_palindromic_subsequence("character"), "");
    SHOW(longest_palindromic_subsequence_length("character"), "");
    SHOW(count_palindromic_subsequences<long long>("aab"), "a, a, b, aa");
    note("答えを mod で出すときは math_v2 の using mint = ModInt998; で count_palindromic_subsequences<mint>(s)。");
    note("long long で溢れるときはメッセージを出して止まる。");

    api("【8】 palindrome_partition_min(s)", "vector<StrSpan>", "回文に分ける最小の個数の分け方（.size() が個数）", "O(N^2)");
    RUN(const auto parts = palindrome_partition_min("abaccd"));
    SHOW(parts.size(), "aba | cc | d");
    SHOW(substr_of("abaccd", parts[1]), "");
    note("N が大きい（10^5 以上）ときや、異なる回文の種類・出現回数は 04_palindromic_tree。");

    api("【9】 count_palindromes_leq<T>(s, A)", "T（mint・long long。long long で溢れたら止める）",
        "A の文字からなる長さ |s| の回文のうち、辞書順（A の番号の順）で s 以下のものの数（前半を σ 進数で読む）", "O(N)");
    SHOW(count_palindromes_leq<long long>("AXA", Alphabet::upper()), "ABC242 E の入力例 1（AAA〜AXA の 24 個）");
    SHOW(count_palindromes_leq<long long>("ABCZAZ", Alphabet::upper()), "");
    note("答えを mod で出すときは count_palindromes_leq<mint>(S, Alphabet::upper())。同じ長さの文字列全体での順位は 05_lex_order の lex_rank_same_length<T>。");
    note("k 文字まで書き換えてよい回文（中心ごと）は 01_rolling_hash の center_palindrome_lengths_with_changes（ABC468 D）。");
}

}  // namespace pal_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    pal_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pal_test::benchmark();    // 速度計測
    pal_demo::demo();         // 使い方デモ
    return 0;
}
