// =============================================================================
//  ランダムテスト
// =============================================================================
namespace sa_test {

{COMMON_TEST}

using VL = std::vector<long long>;
// 愚直：文字列の組（of_strings と同じ並び：各文字列のあとに区切り）
struct Brute {
    std::vector<VL> strs;
    VL all;             // 値（区切りは使わない）
    std::vector<int> owner, lim;
    std::vector<VL> keys;  // 位置 p の接尾辞の並べ替え用の列（区切りは文字より小さく互いに違う）
    int n = 0;
    bool is_multi = false;
    explicit Brute(const std::vector<VL>& vs, bool multi) : strs(vs), is_multi(multi) {
        const int K = (int)vs.size();
        std::vector<long long> code;  // 並べ替え用の値
        for (int t = 0; t < K; ++t) {
            const int st = (int)all.size();
            for (long long x : vs[t]) all.push_back(x), owner.push_back(t), lim.push_back(st + (int)vs[t].size()), code.push_back(x);
            if (multi) all.push_back(0), owner.push_back(-1), lim.push_back((int)all.size() - 1), code.push_back(std::numeric_limits<long long>::min() + t);
        }
        n = (int)all.size();
        for (int p = 0; p < n; ++p) keys.push_back(VL(code.begin() + p, code.end()));
    }
    std::vector<int> sa() const {
        std::vector<int> r(n);
        std::iota(r.begin(), r.end(), 0);
        std::sort(r.begin(), r.end(), [&](int a, int b) { return keys[a] < keys[b]; });
        return r;
    }
    VL sub(int l, int r) const { return VL(all.begin() + l, all.begin() + r); }
    int lcp(int i, int j) const {
        if (i == n || j == n) return 0;
        int k = 0;
        while (i + k < lim[i] && j + k < lim[j] && all[i + k] == all[j + k]) ++k;
        return k;
    }
};

static void check_sa(std::mt19937_64& rng, const SuffixArray& S, const Brute& B, int& cases) {
    const int n = B.n;
    SA_CHECK(S.size() == n && S.sa == B.sa());
    for (int k = 0; k < n; ++k) SA_CHECK(S.rank[S.sa[k]] == k && S.length_at(k) == B.lim[k] - k);
    for (int k = 0; k + 1 < n; ++k) SA_CHECK(S.lcp_adj[k] == B.lcp(S.sa[k], S.sa[k + 1]));
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j) SA_CHECK(S.lcp(i, j) == B.lcp(i, j));
    // 部分文字列の一覧
    std::vector<std::pair<int, int>> subs;  // (l, r)
    for (int l = 0; l < n; ++l)
        for (int r = l + 1; r <= B.lim[l]; ++r) subs.push_back({l, r});
    std::vector<VL> multi;
    for (auto [l, r] : subs) multi.push_back(B.sub(l, r));
    std::sort(multi.begin(), multi.end());
    std::vector<VL> dist = multi;
    dist.erase(std::unique(dist.begin(), dist.end()), dist.end());
    SA_CHECK(S.count_distinct_substrings() == (long long)dist.size() && S.count_substrings_with_dup() == (long long)multi.size());
    for (long long k = 0; k < (long long)dist.size(); ++k) {
        const StrSpan sp = S.kth_substring(k);
        SA_CHECK(B.sub(sp.l, sp.r) == dist[k]);
    }
    for (long long k = 0; k < (long long)multi.size(); ++k) {
        const StrSpan sp = S.kth_substring_with_dup(k);
        SA_CHECK(B.sub(sp.l, sp.r) == multi[k]);
    }
    const auto byl = S.distinct_substrings_by_length();
    std::vector<long long> exl(byl.size(), 0);
    for (const auto& v : dist) {
        SA_CHECK(v.size() < byl.size());
        ++exl[v.size()];
    }
    SA_CHECK(byl == exl);
    for (auto [l, r] : subs) {
        const VL v = B.sub(l, r);
        long long occ = 0;
        for (auto [l2, r2] : subs)
            if (r2 - l2 == r - l && B.sub(l2, r2) == v) ++occ;
        SA_CHECK(S.count_occurrences(l, r) == occ);
        SA_CHECK(S.substring_rank(l, r) == (long long)(std::lower_bound(dist.begin(), dist.end(), v) - dist.begin()));
        ++cases;
    }
    for (int l = 0; l <= n; ++l) SA_CHECK(S.count_occurrences(l, l) == n + (B.is_multi ? 0 : 1));  // 空の部分文字列
    // 比較
    for (int t = 0; t < 60 && !subs.empty(); ++t) {
        auto [l1, r1] = subs[rng() % subs.size()];
        auto [l2, r2] = subs[rng() % subs.size()];
        if (rng() % 4 == 0) r1 = l1;
        const VL a = B.sub(l1, r1), b = B.sub(l2, r2);
        SA_CHECK(S.compare(l1, r1, l2, r2) == (a < b ? -1 : (a == b ? 0 : 1)));
    }
    // パターン
    for (int t = 0; t < 30; ++t) {
        VL pat;
        if (!subs.empty() && rng() % 2) {
            auto [l, r] = subs[rng() % subs.size()];
            pat = B.sub(l, r);
            if (rng() % 3 == 0 && !pat.empty() && pat.back() > -1000000 && pat.back() < 1000000) pat.back() += (rng() % 2 ? 1 : -1);
        } else {
            for (int k = (int)(rng() % 4); k > 0; --k) pat.push_back(B.all.empty() ? 0 : B.all[rng() % B.all.size()]);
        }
        std::vector<int> ex;
        for (int p = 0; p < n; ++p)
            if (p + (int)pat.size() <= B.lim[p] && B.sub(p, p + (int)pat.size()) == pat) ex.push_back(p);
        if (pat.empty()) {  // 空なら置ける位置すべて（1 つの文字列なら 0〜N、of_strings なら 0〜N-1 = 各文字列の 0〜長さ）
            ex.clear();
            for (int p = 0; p < n + (B.is_multi ? 0 : 1); ++p) ex.push_back(p);
        }
        SA_CHECK(S.find_all(pat) == ex && S.count(pat) == (long long)ex.size());
    }
    // 繰り返し
    int best = 0, best_nov = 0;
    for (auto [l, r] : subs)
        for (int q = l + 1; q < n; ++q)
            if (q + (r - l) <= B.lim[q] && B.sub(q, q + (r - l)) == B.sub(l, r)) {
                best = std::max(best, r - l);
                if (q >= r) best_nov = std::max(best_nov, r - l);
            }
    const CommonSubstring c = S.longest_repeated(), d = S.longest_repeated_non_overlapping();
    SA_CHECK(c.len == best && d.len == best_nov);
    if (best) {
        SA_CHECK(c.pos_a < c.pos_b && B.sub(c.pos_a, c.pos_a + best) == B.sub(c.pos_b, c.pos_b + best) && c.pos_b + best <= B.lim[c.pos_b]);
        VL smallest;  // 最長のうち辞書順最小
        bool any = false;
        for (auto [l, r] : subs)
            if (r - l == best && S.count_occurrences(l, r) >= 2 && (!any || B.sub(l, r) < smallest)) smallest = B.sub(l, r), any = true;
        SA_CHECK(B.sub(c.pos_a, c.pos_a + best) == smallest);
    } else {
        SA_CHECK(c.pos_a == -1 && c.pos_b == -1);
    }
    if (best_nov) SA_CHECK(d.pos_a + best_nov <= d.pos_b && B.sub(d.pos_a, d.pos_a + best_nov) == B.sub(d.pos_b, d.pos_b + best_nov));
    // LCP の和
    const auto sw = S.sum_lcp_with_all();
    long long pairs = 0;
    for (int i = 0; i < n; ++i) {
        long long t = 0;
        for (int j = 0; j < n; ++j) t += i == j ? B.lim[i] - i : B.lcp(i, j);
        for (int j = i + 1; j < n; ++j) pairs += B.lcp(i, j);
        SA_CHECK(sw[i] == t);
    }
    SA_CHECK(S.sum_lcp_all_pairs() == pairs);
}

static void test_single(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 13), sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, n, sg);
        const SuffixArray S(s);
        check_sa(rng, S, Brute({VL(s.begin(), s.end())}, false), cases);
        SA_CHECK(make_suffix_array(s) == S.sa && make_lcp_array(s, S.sa) == S.lcp_adj);
        // 数の列（大きな値・負の値）
        VL a(n);
        const long long vals[] = {std::numeric_limits<long long>::min(), -5, 0, 1000000000000000000LL, std::numeric_limits<long long>::max()};
        for (auto& x : a) x = vals[rng() % (2 + sg)];
        check_sa(rng, SuffixArray(a), Brute({a}, false), cases);
    }
}

static void test_multi(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int K = 1 + (int)(rng() % 3), sg = 1 + (int)(rng() % 3);
        std::vector<std::string> vs(K);
        std::vector<VL> vv(K);
        for (int t = 0; t < K; ++t) vs[t] = rnd_str_shaped(rng, (int)(rng() % 6), sg), vv[t] = VL(vs[t].begin(), vs[t].end());
        const SuffixArray S = SuffixArray::of_strings(vs);
        const Brute B(vv, true);
        check_sa(rng, S, B, cases);
        for (int p = 0; p < B.n; ++p) {
            SA_CHECK(S.owner(p) == B.owner[p]);
            if (B.owner[p] >= 0) SA_CHECK(S.position(S.owner(p), S.offset(p)) == p);
        }
        // 2 つ・K 個の最長共通部分文字列
        if (K >= 2) {
            const std::string& a = vs[0];
            const std::string& b = vs[1];
            int best = 0;
            std::set<std::string> alls;
            for (int i = 0; i < (int)a.size(); ++i)
                for (int j = 0; j < (int)b.size(); ++j) {
                    int k = 0;
                    while (i + k < (int)a.size() && j + k < (int)b.size() && a[i + k] == b[j + k]) ++k;
                    if (k > best) best = k, alls.clear();
                    if (k == best && k > 0) alls.insert(a.substr(i, k));
                }
            const CommonSubstring c = longest_common_substring(a, b);
            SA_CHECK(c.len == best);
            if (best) SA_CHECK(a.substr(c.pos_a, best) == b.substr(c.pos_b, best));
            else SA_CHECK(c.pos_a == -1);
            std::vector<std::string> got;
            for (StrSpan sp : longest_common_substrings_all(a, b)) got.push_back(substr_of(a, sp));
            SA_CHECK(got == std::vector<std::string>(alls.begin(), alls.end()));
            // K 個
            int bk = 0;
            for (int i = 0; i < (int)a.size(); ++i)
                for (int L = 1; i + L <= (int)a.size(); ++L) {
                    bool ok = true;
                    for (int t = 1; t < K; ++t) ok &= vs[t].find(a.substr(i, L)) != std::string::npos;
                    if (ok) bk = std::max(bk, L);
                }
            const StrSpan sk = longest_common_substring_k(vs);
            SA_CHECK(sk.len() == bk);
            for (int t = 1; t < K && bk; ++t) SA_CHECK(vs[t].find(substr_of(a, sk)) != std::string::npos);
        }
    }
}

static void test_sais(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 40 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % (it % 8 == 0 ? 3000 : 400)), sg = 1 + (int)(rng() % 4);
        const std::string s = rnd_str_shaped(rng, n, sg);
        std::vector<int> ex(n);
        std::iota(ex.begin(), ex.end(), 0);
        std::sort(ex.begin(), ex.end(), [&](int a, int b) { return s.compare(a, std::string::npos, s, b, std::string::npos) < 0; });
        const auto sa = make_suffix_array(s);
        SA_CHECK(sa == ex);
        const auto lcp = make_lcp_array(s, sa);
        for (int k = 0; k + 1 < n; ++k) {
            int h = 0;
            while (sa[k] + h < n && sa[k + 1] + h < n && s[sa[k] + h] == s[sa[k + 1] + h]) ++h;
            SA_CHECK(lcp[k] == h);
        }
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_single(rng, cases);
    test_multi(rng, cases);
    test_sais(rng, cases);
    std::cout << "[suffix_array] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[suffix_array] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 1000000;
    const std::string r26 = rnd_str(rng, n, 26), r2 = rnd_str_shaped(rng, n, 2), aa(n, 'a');
    auto t0 = std::chrono::steady_clock::now();
    const auto s1 = make_suffix_array(r26);
    const auto s2 = make_suffix_array(r2);
    const auto s3 = make_suffix_array(aa);
    const long long t_sa = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto l1 = make_lcp_array(r26, s1);
    const long long t_lcp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const SuffixArray S(r2);
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = 0;
    for (int i = 0; i < 1000000; ++i) x += S.lcp((int)(rng() % n), (int)(rng() % n));
    const long long t_q = ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < 100000; ++i) x += S.kth_substring_with_dup((long long)(rng() % (unsigned long long)S.count_substrings_with_dup())).len();
    const long long t_k = ms(t0);
    std::vector<long long> v(n);
    for (auto& y : v) y = rnd_ll(rng);
    t0 = std::chrono::steady_clock::now();
    const auto s4 = make_suffix_array(v);
    const long long t_v = ms(t0);
    std::cout << "[suffix_array] benchmark: make_suffix_array（N = 10^6。σ = 26 / 2 / 全部 a）合計 " << t_sa << " ms, make_lcp_array " << t_lcp
              << " ms, SuffixArray 構築 " << t_build << " ms, lcp(i, j) × 10^6（表の構築込み）" << t_q << " ms, kth_substring_with_dup × 10^5（初回の準備込み）"
              << t_k << " ms, long long の列（N = 10^6）" << t_v << " ms  (check " << (s1[7] + s2[5] + s3[3] + l1[9] + x + s4[1]) % 7 << ")\n";
}

}  // namespace sa_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace sa_demo {

{COMMON_DEMO}

static void demo() {
    banner("接尾辞配列（suffix_array） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  make_suffix_array / make_lcp_array                   接尾辞配列と隣どうしの LCP（ACL の suffix_array / lcp_array と同じ）\n"
              << "  SuffixArray SA(s)：sa / rank / lcp_adj / lcp(i, j) / compare / find_range / count / find_all / count_occurrences\n"
              << "     count_distinct_substrings / distinct_substrings_by_length / kth_substring / substring_rank / kth_substring_with_dup\n"
              << "     longest_repeated / longest_repeated_non_overlapping / sum_lcp_with_all / sum_lcp_all_pairs\n"
              << "  SuffixArray::of_strings(vs)：複数の文字列（owner / offset / position）\n"
              << "  longest_common_substring / longest_common_substrings_all / longest_common_substring_k   最長共通部分文字列（連続）\n";

    api("【1】 make_suffix_array(s) / make_lcp_array(s, sa)", "vector<int>",
        "接尾辞を辞書順に並べた開始位置（SA-IS）と、隣どうしの最長共通接頭辞（長さ N - 1）。数の列にも", "O(N)");
    RUN(const std::string S = "banana");
    SHOW(make_suffix_array(S), "a, ana, anana, banana, na, nana");
    SHOW(make_lcp_array(S, make_suffix_array(S)), "");

    api("【2】 SuffixArray SA(s) / SA.lcp(i, j) / SA.compare(l1, r1, l2, r2)", "int",
        "任意の 2 つの接尾辞の LCP・部分文字列の辞書順を O(1) で（初回に疎テーブルを作る）", "O(1)");
    RUN(const SuffixArray SA(S));
    SHOW(SA.sa, "");
    SHOW(SA.rank, "rank[p] = 位置 p の順位");
    SHOW(SA.lcp(1, 3), "\"anana\" と \"ana\"");
    SHOW(SA.compare(0, 2, 2, 4), "\"ba\" と \"na\"");

    api("【3】 SA.count(pat) / SA.find_all(pat) / SA.count_occurrences(l, r)", "long long / vector<int>",
        "パターンの出現回数・位置（何度でも O(M log N)）/ 部分文字列 s[l, r) の出現回数（O(log N)）", "—");
    SHOW(SA.count("ana"), "重なりも");
    SHOW(SA.find_all("an"), "");
    SHOW(SA.count_occurrences(1, 3), "\"an\"");
    note("ABC362 G（Count Substring Query）: クエリごとに SA.count(T_i)。");

    api("【4】 count_distinct_substrings / distinct_substrings_by_length / kth_substring(k) / substring_rank(l, r)", "long long / vector / StrSpan",
        "異なる部分文字列の数・長さごとの数・辞書順 k 番目（0 から）・s[l, r) が何番目か", "O(N) / O(log N)");
    SHOW(SA.count_distinct_substrings(), "a, an, ana, … の 15 個");
    SHOW(SA.distinct_substrings_by_length(), "r[len]");
    SHOW(substr_of(S, SA.kth_substring(0)), "辞書順で最初");
    SHOW(substr_of(S, SA.kth_substring(4)), "a, an, ana, anan, anana");
    SHOW(SA.substring_rank(1, 4), "\"ana\" より小さい異なる部分文字列の数");

    api("【5】 kth_substring_with_dup(k) / count_substrings_with_dup()", "StrSpan / long long",
        "位置が違えば別として数えたときの辞書順 k 番目（同じ文字列が出現回数だけ並ぶ）。複数の文字列でも（ABC280 Ex）", "O(log N)（初回 O(N)）");
    SHOW(SA.count_substrings_with_dup(), "N(N+1)/2");
    SHOW(substr_of(S, SA.kth_substring_with_dup(2)), "a, a, a, an, …");
    SHOW(substr_of(S, SA.kth_substring_with_dup(3)), "");

    api("【6】 longest_repeated() / longest_repeated_non_overlapping()", "CommonSubstring{len, pos_a, pos_b}",
        "2 か所以上に現れる最長の部分文字列（重なってよい / 重ならない）。ないときは len = 0", "O(N) / O(N log N)");
    SHOW(SA.longest_repeated().len, "\"ana\"（位置 1 と 3。重なる）");
    SHOW(SA.longest_repeated_non_overlapping().len, "\"an\" か \"na\"");
    note("ABC141 E（Who Says a Pun?）は longest_repeated_non_overlapping().len。");

    api("【7】 sum_lcp_with_all() / sum_lcp_all_pairs()", "vector<long long> / long long",
        "各 i で Σ_j lcp(i, j)（j = i も含む）/ すべての組の LCP の和", "O(N)");
    SHOW(SA.sum_lcp_with_all(), "ABC213 F（Common Prefixes）の答えそのもの");
    SHOW(SA.sum_lcp_all_pairs(), "");

    api("【8】 SuffixArray::of_strings(vs) / owner(p) / offset(p) / position(t, off)", "SuffixArray",
        "複数の文字列をまとめて（区切りはどの文字より小さい。部分文字列は 1 つの文字列の中だけ）。sa には区切りの位置（owner = -1）も入る", "O(全体)");
    RUN(const SuffixArray M = SuffixArray::of_strings(std::vector<std::string>{"ab", "b"}));
    SHOW(M.sa, "位置 2 と 4 は区切り");
    SHOW(M.owner(3), "位置 3 は 1 番目の文字列");
    SHOW(M.count("b"), "両方の文字列での出現");
    SHOW(M.count_distinct_substrings(), "a, ab, b（文字列をまたがない）");

    api("【9】 longest_common_substring(a, b) / longest_common_substrings_all(a, b) / longest_common_substring_k(vs)", "CommonSubstring / vector<StrSpan> / StrSpan",
        "最長共通部分文字列（連続する部分）。2 つ・すべて・K 個。飛び飛びの「部分列」の LCS は 06_lcs_edit", "O(全体)");
    RUN(const CommonSubstring c = longest_common_substring("xabcde", "abcyabc"));
    SHOW(c.len, "");
    SHOW(c.pos_a, "a での位置");
    SHOW(c.pos_b, "b での位置");
    SHOW(longest_common_substrings_all("abxcd", "cdyab").size(), "\"ab\" と \"cd\"");
    SHOW(substr_of("abcd", longest_common_substring_k(std::vector<std::string>{"abcd", "zbcx", "bcq"})), "3 つすべてに");
}

}  // namespace sa_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    sa_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    sa_test::benchmark();    // 速度計測
    sa_demo::demo();         // 使い方デモ
    return 0;
}
