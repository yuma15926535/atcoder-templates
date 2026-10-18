// =============================================================================
//  ランダムテスト
// =============================================================================
namespace sam_test {

{COMMON_TEST}

// 愚直：文字列の組（つながらない）の全部の部分文字列
struct Brute {
    std::vector<std::string> ss;
    std::map<std::string, long long> occ;      // 出現回数
    std::map<std::string, long long> first;    // 最初の開始位置（通し番号）
    explicit Brute(const std::vector<std::string>& v) : ss(v) {
        long long base = 0;
        for (const auto& s : ss) {
            for (int l = 0; l < (int)s.size(); ++l)
                for (int r = l + 1; r <= (int)s.size(); ++r) {
                    const std::string t = s.substr(l, r - l);
                    ++occ[t];
                    if (!first.count(t)) first[t] = base + l;
                    else first[t] = std::min(first[t], base + l);
                }
            base += (long long)s.size();
        }
    }
};

static void check(std::mt19937_64& rng, const SuffixAutomaton& sam, const Brute& B, int sg) {
    long long total_len = 0, dup = 0, chars = 0;
    for (const auto& s : B.ss) chars += (long long)s.size();
    for (auto& [t, c] : B.occ) total_len += (long long)t.size(), dup += c;
    SAM_CHECK(sam.count_distinct_substrings() == (long long)B.occ.size() && sam.total_length_distinct() == total_len && sam.count_substrings_with_dup() == dup);
    // 部分文字列か・回数・最初
    for (int it = 0; it < 30; ++it) {
        std::string t = rnd_str(rng, (int)(rng() % 6), sg + 1);
        if (!B.occ.empty() && rng() % 2) {
            auto itr = B.occ.begin();
            std::advance(itr, rng() % B.occ.size());
            t = itr->first;
        }
        const bool in = B.occ.count(t) > 0;
        SAM_CHECK(sam.contains(t) == (in || t.empty()));
        SAM_CHECK(sam.count_occurrences(t) == (t.empty() ? chars + (long long)std::max<std::size_t>(B.ss.size(), 1) : (in ? B.occ.at(t) : 0)));
        SAM_CHECK(sam.first_occurrence(t) == (t.empty() ? 0 : (in ? B.first.at(t) : -1)));
        // 最長共通部分文字列（t の中で）
        int best = 0, bl = 0;
        for (int l = 0; l < (int)t.size(); ++l)
            for (int r = l + 1; r <= (int)t.size(); ++r)
                if (B.occ.count(t.substr(l, r - l)) && r - l > best) best = r - l, bl = l;
        const StrSpan sp = sam.longest_common_substring(t);
        SAM_CHECK(sp.len() == best && (best == 0 || sp.l == bl));
    }
    // k 番目
    std::vector<std::string> dist, dups;
    for (auto& [t, c] : B.occ) {
        dist.push_back(t);
        for (long long k = 0; k < c; ++k) dups.push_back(t);
    }
    for (long long k = 0; k < (long long)dist.size(); ++k) SAM_CHECK(sam.kth_substring(k) == dist[k]);
    for (long long k = 0; k < (long long)dups.size(); ++k) SAM_CHECK(sam.kth_substring(k, true) == dups[k]);
    // 状態
    const auto ord = sam.order_by_len();
    SAM_CHECK((int)ord.size() == sam.size() && ord[0] == 0 && sam.link(0) == -1);
    for (int k = 1; k < sam.size(); ++k) SAM_CHECK(sam.len(ord[k - 1]) <= sam.len(ord[k]));
    for (int v = 1; v < sam.size(); ++v) SAM_CHECK(sam.len(sam.link(v)) < sam.len(v));
    // 状態ごとの回数：その状態の最長の文字列の回数
    for (auto& [t, c] : B.occ) {
        int v = 0;
        for (char ch : t) v = sam.next(v, ch);
        SAM_CHECK(v > 0 && sam.endpos_count(v) == c && sam.first_end(v) == B.first.at(t) + (long long)t.size() - 1);
        SAM_CHECK(sam.len(sam.link(v)) < (int)t.size() && (int)t.size() <= sam.len(v));
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 25), sg);
        const SuffixAutomaton sam(s);
        SAM_CHECK(s.size() < 2 || sam.size() <= 2 * (int)s.size() - 1);
        check(rng, sam, Brute({s}), sg);
        // 1 文字ずつ足しても同じ
        SuffixAutomaton s2;
        for (char c : s) s2.extend(c);
        SAM_CHECK(s2.size() == sam.size() && s2.count_distinct_substrings() == sam.count_distinct_substrings());
        ++cases;
    }
    // 複数の文字列
    for (int it = 0; it < 800 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        std::vector<std::string> v(1 + rng() % 4);
        for (auto& s : v) s = rnd_str_shaped(rng, (int)(rng() % 10), sg);
        SuffixAutomaton sam;
        for (const auto& s : v) sam.add_string(s);
        check(rng, sam, Brute(v), sg);
        ++cases;
    }
    // extend・add_string と問い合わせを混ぜる（数え直し・last・続きの文字列）
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        SuffixAutomaton sam;
        std::vector<std::string> v;
        for (int op = 0; op < 15; ++op) {
            if (v.empty() || rng() % 4 == 0) {
                v.push_back(rnd_str(rng, (int)(rng() % 4), sg));
                sam.add_string(v.back());
            } else {
                const char c = (char)('a' + rng() % sg);
                v.back() += c;
                sam.extend(c);
            }
            const Brute B(v);
            long long total_len = 0;
            for (auto& [t, c] : B.occ) total_len += (long long)t.size();
            SAM_CHECK(sam.count_distinct_substrings() == (long long)B.occ.size() && sam.total_length_distinct() == total_len);
            int w = 0;
            for (char ch : v.back()) w = sam.next(w, ch);
            SAM_CHECK(sam.last() == w);
            const std::string t = rnd_str(rng, (int)(rng() % 3), sg);
            long long chars = 0;
            for (auto& s : v) chars += (long long)s.size();
            SAM_CHECK(sam.count_occurrences(t) == (t.empty() ? chars + (long long)v.size() : (B.occ.count(t) ? B.occ.at(t) : 0)));
        }
        ++cases;
    }
    // 数の列
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        std::vector<int> a(rng() % 20);
        for (auto& x : a) x = (int)(rng() % 3);
        const SuffixAutomaton sam(a, Alphabet::ints(3));
        std::set<std::vector<int>> st;
        for (int l = 0; l < (int)a.size(); ++l)
            for (int r = l + 1; r <= (int)a.size(); ++r) st.insert(std::vector<int>(a.begin() + l, a.begin() + r));
        SAM_CHECK(sam.count_distinct_substrings() == (long long)st.size());
        long long k = 0;
        for (const auto& t : st) SAM_CHECK(sam.kth_substring_ids(k++) == t);
        ++cases;
    }
    std::cout << "[suffix_automaton] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[suffix_automaton] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::string s = rnd_str(rng, 500000, 26);
    auto t0 = std::chrono::steady_clock::now();
    const SuffixAutomaton sam(s);
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = sam.count_distinct_substrings();
    for (int q = 0; q < 100000; ++q) {
        const int l = (int)(rng() % 499990);
        x += sam.count_occurrences(s.substr(l, 1 + rng() % 5));
    }
    const long long t_q = ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (int q = 0; q < 100; ++q) x += (long long)sam.kth_substring((long long)(rng() % (unsigned long long)sam.count_substrings_with_dup()), true).size();
    const long long t_k = ms(t0);
    const std::string ab = rnd_str(rng, 1000000, 2);
    t0 = std::chrono::steady_clock::now();
    const SuffixAutomaton s2(ab, Alphabet::of("ab"));
    x += s2.size();
    const long long t_ab = ms(t0);
    std::cout << "[suffix_automaton] benchmark: 作成（N = 5×10^5, σ = 26）" << t_build << " ms, count_occurrences × 10^5（初回の数え込み）" << t_q
              << " ms, kth_substring(重複あり。答えの長さは平均 N / 3) × 100 " << t_k << " ms, 作成（N = 10^6, σ = 2）" << t_ab << " ms  (check " << x % 7 << ")\n";
}

}  // namespace sam_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace sam_demo {

{COMMON_DEMO}

static void demo() {
    banner("接尾辞オートマトン（suffix_automaton） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  contains / count_occurrences / first_occurrence / count_distinct_substrings / total_length_distinct\n"
              << "  kth_substring / kth_substring_ids / count_substrings_with_dup / longest_common_substring\n"
              << "  extend / add_string / size / len / link / next / next_id / endpos_count / first_end / order_by_len / last   状態（DP 用）\n";

    RUN(const SuffixAutomaton sam(std::string("missisippi")));

    api("【1】 contains(t) / count_occurrences(t) / first_occurrence(t)", "bool / long long / long long", "部分文字列か / 出現回数（重なりも）/ 最初の開始位置（なければ -1）", "O(|t|)");
    SHOW(sam.count_occurrences("i"), "ABC362 G の入力例 1：i が 4 回");
    SHOW(sam.count_occurrences("is"), "");
    SHOW(sam.contains("sip"), "");
    SHOW(sam.first_occurrence("ss"), "");

    api("【2】 count_distinct_substrings() / total_length_distinct() / count_substrings_with_dup()", "long long",
        "異なる部分文字列の数 / その長さの合計（この 2 つは足すたびに更新するので 1 文字ごとに聞いてもよい）/ 位置が違えば別として数えた数", "O(1) / O(1) / O(状態数 × σ)");
    SHOW(sam.count_distinct_substrings(), "");
    SHOW(sam.total_length_distinct(), "");
    SHOW(sam.count_substrings_with_dup(), "10 × 11 / 2");

    api("【3】 kth_substring(k, with_dup)", "string（kth_substring_ids は番号の列）", "辞書順 k 番目（0 から）。with_dup = true なら同じ文字列も出現の数だけ数える", "O(|答え| σ)");
    SHOW(sam.kth_substring(0), "");
    SHOW(sam.kth_substring(5), "i, ip, ipp, ippi, is, isi, ...");
    SHOW(sam.kth_substring(5, true), "i が 4 回、ip、ipp …");

    api("【4】 longest_common_substring(t)", "StrSpan（t での位置）", "t の部分文字列で、こちらにも現れる最長のもの（同じ長さなら左）", "O(|t|)");
    SHOW(sam.longest_common_substring("mississippi"), "\"missis\"");

    api("【5】 extend(c) / add_string(t)", "なし", "1 文字ずつ足す（最後に足した文字列の続き。オンライン）/ 別の文字列として足す（文字列の組の部分文字列をまとめて数える）", "O(σ)（ならし）");
    {
        RUN(SuffixAutomaton g);
        RUN(g.add_string("abab"), g.add_string("cab"));
        SHOW(g.count_distinct_substrings(), "abab と cab の部分文字列の種類");
        SHOW(g.count_occurrences("ab"), "2 つの文字列で合わせて 3 回");
        SHOW(g.first_occurrence("ca"), "通し番号（abab の 4 文字の後の 0 文字目）");
        SHOW(g.kth_substring(13, true), "ABC280 Ex の入力例 1 の 14 番目（2 番目の文字列の 1 文字目。どの文字列かは first_occurrence の通し番号から）");
        note("kth_substring は答えの長さに比例するので、答えが長くなる問い合わせを多数（ABC280 Ex の本番）なら 03_suffix_array の kth_substring_with_dup（StrSpan で返す。O(log N)）");
    }

    api("【6】 状態：len(v) / link(v) / next(v, c) / endpos_count(v) / order_by_len()", "int / int / int / long long / vector<int>",
        "状態 v は「出現位置の集合が同じ部分文字列」のまとまり。最長の長さ len、最短は len(link) + 1。order_by_len の逆順にたどると部分木の和（DP 用）",
        "O(1)（order_by_len は O(状態数 + N)）");
    {
        int v = 0;
        for (char c : std::string("ss")) v = sam.next(v, c);
        SHOW(sam.len(v), "\"ss\" の状態の最長（\"miss\"）");
        SHOW(sam.len(sam.link(v)), "最短は len(link) + 1 = 2");
        SHOW(sam.endpos_count(v), "");
        SHOW(sam.size(), "状態の数（2N 以下）");
    }
}

}  // namespace sam_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    sam_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    sam_test::benchmark();    // 速度計測
    sam_demo::demo();         // 使い方デモ
    return 0;
}
