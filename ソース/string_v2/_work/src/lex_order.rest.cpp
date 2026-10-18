// =============================================================================
//  ランダムテスト
// =============================================================================
namespace lex_test {

{COMMON_TEST}

static void test_concat(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int k = (int)(rng() % 6);
        std::vector<std::string> vs(k);
        for (auto& s : vs) s = rnd_str_shaped(rng, (int)(rng() % 4), 1 + (int)(rng() % 2));
        std::vector<int> p(k);
        std::iota(p.begin(), p.end(), 0);
        std::string best;
        bool first = true;
        do {
            std::string t;
            for (int i : p) t += vs[i];
            if (first || t < best) best = t, first = false;
        } while (std::next_permutation(p.begin(), p.end()));
        LEX_CHECK(min_concat(vs) == best);
        const auto ord = concat_order(vs);
        std::string t;
        for (int i : ord) t += vs[i];
        LEX_CHECK(t == best);
        ++cases;
    }
}

static void test_order(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        std::string ord = "abcdefghijklmnopqrstuvwxyz";
        std::shuffle(ord.begin(), ord.end(), rng);
        const LexLess lt = lex_less_by(ord);
        std::vector<int> pos(256);
        for (int i = 0; i < 26; ++i) pos[(unsigned char)ord[i]] = i;
        for (int t = 0; t < 20; ++t) {
            const std::string a = rnd_str(rng, (int)(rng() % 4), 3), b = rnd_str(rng, (int)(rng() % 4), 3);
            std::vector<int> x, y;
            for (char c : a) x.push_back(pos[(unsigned char)c]);
            for (char c : b) y.push_back(pos[(unsigned char)c]);
            LEX_CHECK(lt(a, b) == (x < y) && lt.compare(a, b) == (x < y ? -1 : (x == y ? 0 : 1)));
            ++cases;
        }
    }
}

static void test_rotation(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 5000 * TEST_SCALE; ++it) {
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 14), 1 + (int)(rng() % 3));
        const int n = (int)s.size();
        int mn = 0, mx = 0;
        std::string bmn = s, bmx = s;
        for (int k = 0; k < n; ++k) {
            const std::string r = s.substr(k) + s.substr(0, k);
            if (r < bmn) bmn = r, mn = k;
            if (r > bmx) bmx = r, mx = k;
        }
        LEX_CHECK(min_rotation(s) == mn && max_rotation(s) == mx);
        LEX_CHECK(min_rotation(std::vector<int>(s.begin(), s.end())) == mn);
        ++cases;
    }
}

static void test_names(std::mt19937_64& rng, int& cases) {
    // 列名：1 から順に作ったものと
    std::string cur;
    for (long long n = 1; n <= 20000; ++n) {
        // cur を「次」に（z の次は aa）
        int i = (int)cur.size() - 1;
        while (i >= 0 && cur[i] == 'z') cur[i--] = 'a';
        if (i < 0) cur.insert(cur.begin(), 'a');
        else ++cur[i];
        LEX_CHECK(bijective_to_string(n) == cur && bijective_from_string(cur) == n);
        ++cases;
    }
    LEX_CHECK(bijective_to_string(0).empty() && bijective_to_string(28, Alphabet::upper()) == "AB");
    for (int it = 0; it < 1000 * TEST_SCALE; ++it) {
        const long long n = rnd_ll(rng) & std::numeric_limits<long long>::max();
        LEX_CHECK(bijective_from_string(bijective_to_string(n)) == n);
        const Alphabet A = Alphabet::of("xy");
        LEX_CHECK(bijective_from_string(bijective_to_string(n, A), A) == n);
    }
    // 辞書の番号：全部の文字列を並べたものと
    for (int sg = 1; sg <= 3; ++sg)
        for (int L = 1; L <= 6; ++L) {
            std::string chars = std::string("abc").substr(0, sg);
            const Alphabet A = Alphabet::of(chars);
            std::vector<std::string> all;
            std::vector<std::string> layer = {""};
            for (int len = 1; len <= L; ++len) {
                std::vector<std::string> nl;
                for (const auto& w : layer)
                    for (char c : chars) nl.push_back(w + c);
                all.insert(all.end(), nl.begin(), nl.end());
                layer = nl;
            }
            std::sort(all.begin(), all.end());
            for (long long k = 0; k < (long long)all.size(); ++k) {
                LEX_CHECK(dict_rank(all[k], L, A) == k && dict_unrank(k, L, A) == all[k]);
                ++cases;
            }
            // 同じ長さの次
            for (const auto& w : layer) {
                auto it2 = std::upper_bound(layer.begin(), layer.end(), w);
                const std::string nx = it2 == layer.end() ? "" : *it2;
                LEX_CHECK(next_string(w, A) == nx);
            }
        }
    // 大きな max_len
    LEX_CHECK(dict_rank("aaa", 1000000000000LL, Alphabet::of("a")) == 2 && dict_rank("aab", 3) == 3 && dict_unrank(0, 1000000000000000000LL) == "a");
    LEX_CHECK(dict_rank("ab", 40, Alphabet::of("ab")) == (1LL << 39));  // "a" と、"aa" で始まる 2^39 - 1 個
    LEX_CHECK(dict_unrank(1LL << 39, 40, Alphabet::of("ab")) == "ab");
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
static void test_final_additions(std::mt19937_64& rng, int& cases) {  // lex_rank_same_length
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3), n = (int)(rng() % 7);
        std::string letters = std::string("xyz").substr(0, sg);
        std::shuffle(letters.begin(), letters.end(), rng);
        const Alphabet A = Alphabet::of(letters);
        std::string s(n, ' ');
        for (char& c : s) c = letters[rng() % sg];
        const std::vector<int> sid = to_ids(s, A);
        long long cnt = 0;
        for (const auto& v : all_id_strings(n, sg)) cnt += v < sid;
        LEX_CHECK(lex_rank_same_length<long long>(s, A) == cnt && lex_rank_same_length<double>(s, A) == (double)cnt);
        ++cases;
    }
    LEX_CHECK(lex_rank_same_length<long long>(std::vector<int>{1, 0, 1}, Alphabet::ints(2)) == 5);
}
static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_concat(rng, cases);
    test_order(rng, cases);
    test_rotation(rng, cases);
    test_names(rng, cases);
    test_final_additions(rng, cases);
    std::cout << "[lex_order] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[lex_order] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 10000000;
    const std::string s = rnd_str_shaped(rng, n, 2), a(n, 'a');
    auto t0 = std::chrono::steady_clock::now();
    const long long r = min_rotation(s) + max_rotation(s) + min_rotation(a);
    const long long t_rot = ms(t0);
    std::vector<std::string> vs(200000);
    for (auto& x : vs) x = rnd_str(rng, 1 + (int)(rng() % 10), 3);
    t0 = std::chrono::steady_clock::now();
    const auto ord = concat_order(vs);
    const long long t_cat = ms(t0);
    std::vector<std::string> ws(200000);
    for (auto& x : ws) x = rnd_str(rng, 10, 26);
    t0 = std::chrono::steady_clock::now();
    std::sort(ws.begin(), ws.end(), lex_less_by("qwertyuiopasdfghjklzxcvbnm"));
    const long long t_sort = ms(t0);
    std::cout << "[lex_order] benchmark: min_rotation + max_rotation + min_rotation（N = 10^7）" << t_rot << " ms, concat_order（2×10^5 個）" << t_cat
              << " ms, lex_less_by で sort（2×10^5 個）" << t_sort << " ms  (check " << (r + ord[5] + ws[3][0]) % 7 << ")\n";
}

}  // namespace lex_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace lex_demo {

{COMMON_DEMO}

static void demo() {
    banner("辞書順（lex_order） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  concat_order / min_concat                 つなげて辞書順最小（a + b < b + a の順）\n"
              << "  lex_less_by(order) / LexLess::compare    独自の文字の順の比較（sort に渡す）\n"
              << "  min_rotation / max_rotation              辞書順最小・最大の回転の開始位置\n"
              << "  bijective_to_string / bijective_from_string   列名（a, …, z, aa, …）と番号\n"
              << "  dict_rank / dict_unrank                  長さ max_len 以下の全部の文字列での辞書順の番号\n"
              << "  next_string                              同じ長さで辞書順の次\n"
              << "  lex_rank_same_length<T>(s, A)            同じ長さの文字列の中での辞書順の番号（mod でも）\n";

    api("【1】 concat_order(vs) / min_concat(vs)", "vector<int> / 列",
        "全部をつなげて辞書順最小になる並べ方（a + b < b + a で並べる。単に辞書順に並べるのとは違う）", "O(N log N × 長さ)");
    RUN(const std::vector<std::string> vs = {"b", "ba", "bab"});
    SHOW(min_concat(vs), "辞書順の b, ba, bab → \"bbabab\" より小さい");
    SHOW(concat_order(vs), "");
    note("ABC225 F（String Cards）: concat_order で並べてから、後ろから「K 枚選ぶ」DP（dp[i][j] = min(dp[i+1][j], S_i + dp[i+1][j-1])）。");

    api("【2】 lex_less_by(order) / LexLess::compare(a, b)", "比較関数",
        "order に並べた文字の順での辞書順（order にない文字は止まる）。sort にそのまま渡せる", "O(長さ)");
    {
        std::vector<std::string> names = {"abx", "bzz", "bzy", "caa"};
        std::sort(names.begin(), names.end(), lex_less_by("bacdefghijklmnopqrstuvwxzy"));
        note("ABC219 C（Neo-lexicographic Ordering）の入力例 1：sort(names, lex_less_by(X))");
        std::cout << "    → " << join(names) << "\n";
    }
    SHOW(lex_less_by("cba").compare("ab", "ba"), "c < b < a なので \"ab\" > \"ba\"");

    api("【3】 min_rotation(s) / max_rotation(s)", "int",
        "辞書順最小 / 最大の回転 s[k..] + s[..k] の k（同じ回転が何度も出るなら最小の k）", "O(N)");
    RUN(const std::string S = "cabca");
    SHOW(min_rotation(S), "\"acabc\"");
    SHOW(max_rotation(S), "\"cabca\"");
    {
        // ABC223 B: 回転の中で辞書順最小と最大
        const std::string T = "aaba";
        const int a = min_rotation(T), b = max_rotation(T);
        note("ABC223 B（String Shifting）: min_rotation と max_rotation の位置で回したもの");
        std::cout << "    T = aaba → " << T.substr(a) + T.substr(0, a) << " / " << T.substr(b) + T.substr(0, b) << "\n";
    }

    api("【4】 bijective_to_string(n, A) / bijective_from_string(s, A)", "string / long long",
        "列名（全単射の記数法）：1 → a、26 → z、27 → aa、28 → ab。0 は空。A で文字を変えられる", "O(log n)");
    SHOW(bijective_to_string(2), "ABC171 C（One Quadrillion and One Dalmatians）の入力例 1");
    SHOW(bijective_to_string(27), "");
    SHOW(bijective_to_string(123456789), "ABC171 C の入力例 3");
    SHOW(bijective_from_string("aa"), "ABC285 C");
    SHOW(bijective_from_string("BRUTMHYHIIZP", Alphabet::upper()), "ABC285 C の入力例 3");

    api("【5】 dict_rank(s, max_len, A) / dict_unrank(k, max_len, A)", "long long / string",
        "長さ 1〜max_len の全部の文字列を辞書順に並べたときの番号（0 から）とその逆。a, aa, aaa, …, ab, … の順", "O(max_len)");
    SHOW(dict_rank("a", 3), "最初");
    SHOW(dict_rank("ab", 3), "a, aa, aaa, …, aaz, ab");
    SHOW(dict_unrank(28, 3), "");
    SHOW(dict_rank("ba", 2, Alphabet::of("ab")), "a, aa, ab, b, ba");

    api("【6】 next_string(s, A)", "string", "同じ長さで辞書順の次（\"az\" → \"ba\"。最後の \"zz…\" なら空）", "O(N)");
    SHOW(next_string("az"), "");
    SHOW(next_string("zz"), "ない → 空");
    SHOW(next_string("011", Alphabet::of("01")), "2 進数の +1");

    api("【7】 lex_rank_same_length<T>(s, A)", "T（mint・long long。long long で溢れたら止める）",
        "A の文字からなる長さ |s| の文字列で s より小さいものの数（s を σ 進数で読んだ値。N = 10^6 でも mint で）", "O(N)");
    SHOW(lex_rank_same_length<long long>("ba"), "aa〜az の 26 個");
    SHOW(lex_rank_same_length<long long>(std::vector<int>{1, 0, 1}, Alphabet::ints(2)), "2 進数の 101");
    note("長さ max_len 以下の全部の文字列での順位は dict_rank。同じ文字の個数が決まった並べ替えの中での順位・k 番目は math_v2 04_permutation の multiset_perm_rank / multiset_perm_unrank（ABC202 D・ABC215 C）。");
}

}  // namespace lex_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    lex_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    lex_test::benchmark();    // 速度計測
    lex_demo::demo();         // 使い方デモ
    return 0;
}
