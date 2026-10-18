// =============================================================================
//  ランダムテスト
// =============================================================================
namespace trie_test {

{COMMON_TEST}

using Ids = std::vector<int>;
static int lcp_ids(const Ids& a, const Ids& b) {
    std::size_t k = 0;
    while (k < a.size() && k < b.size() && a[k] == b[k]) ++k;
    return (int)k;
}
static bool is_prefix(const Ids& p, const Ids& s) { return p.size() <= s.size() && std::equal(p.begin(), p.end(), s.begin()); }

// 愚直：番号の列 → 数 の map（map の順 = 番号の辞書順）
template <class Make>
static void run_ops(std::mt19937_64& rng, const Alphabet& A, int sg, Make make, int& cases) {
    Trie TR(A);
    std::map<Ids, long long> M;
    std::vector<Ids> pool;  // 使った文字列（何度も使う）
    auto rnd_ids = [&]() {
        if (!pool.empty() && rng() % 2) {
            Ids s = pool[rng() % pool.size()];
            if (rng() % 3 == 0 && !s.empty()) s.resize(rng() % s.size());  // 接頭辞
            else if (rng() % 3 == 0) s.push_back((int)(rng() % sg));
            return s;
        }
        Ids s(rng() % 6);
        for (auto& c : s) c = (int)(rng() % sg);
        pool.push_back(s);
        return s;
    };
    for (int q = 0; q < 300; ++q) {
        const int op = (int)(rng() % 10);
        const Ids s = rnd_ids();
        const auto S = make(s);  // 文字列 / 数の列
        if (op <= 2) {
            const long long cnt = rng() % 4 == 0 ? 0 : 1 + (long long)(rng() % 3);
            const int v = TR.insert(S, cnt);
            M[s] += cnt;
            TRIE_CHECK(v == TR.find(S) && TR.depth(v) == (int)s.size());
        } else if (op <= 4) {
            const long long cnt = (long long)(rng() % 3);
            const long long have = M.count(s) ? M[s] : 0;
            TRIE_CHECK(TR.erase(S, cnt) == (have >= cnt));
            if (have >= cnt) M[s] -= cnt;
        } else if (op == 5 && rng() % 3 == 0) {
            long long removed = 0;
            for (auto& [t, c] : M)
                if (is_prefix(s, t)) removed += c, c = 0;
            TRIE_CHECK(TR.erase_prefix(S) == removed);
        }
        // 数える
        long long total = 0, cnt_s = 0, pre = 0, pre_of = 0, rk = 0, sum_with = 0;
        int best = 0, best_ex = 0;
        bool have_s = false;
        for (const auto& [t, c] : M) {
            if (c == 0) continue;
            total += c;
            if (t == s) cnt_s = c, have_s = true;
            if (is_prefix(s, t)) pre += c;
            if (is_prefix(t, s)) pre_of += c;
            if (t < s) rk += c;
            sum_with += c * lcp_ids(s, t);
            best = std::max(best, lcp_ids(s, t));
            if (t != s || c >= 2) best_ex = std::max(best_ex, lcp_ids(s, t));
        }
        TRIE_CHECK(TR.size() == total && TR.empty() == (total == 0));
        TRIE_CHECK(TR.count(S) == cnt_s && TR.contains(S) == have_s);
        TRIE_CHECK(TR.count_prefix(S) == pre && TR.count_prefixes_of(S) == pre_of);
        TRIE_CHECK(TR.rank(S) == rk && TR.sum_lcp_with(S) == sum_with);
        TRIE_CHECK(TR.max_lcp(S) == best);
        if (have_s) TRIE_CHECK(TR.max_lcp(S, true) == best_ex);
        if (q % 10 == 0) {
            // 辞書順・組の LCP の和
            std::vector<Ids> all;
            for (const auto& [t, c] : M)
                for (long long k = 0; k < c; ++k) all.push_back(t);
            TRIE_CHECK(TR.words_ids() == all);
            for (long long k = 0; k < (long long)all.size(); ++k) TRIE_CHECK(TR.kth_ids(k) == all[k]);
            long long pairs = 0;
            for (std::size_t i = 0; i < all.size(); ++i)
                for (std::size_t j = i + 1; j < all.size(); ++j) pairs += lcp_ids(all[i], all[j]);
            TRIE_CHECK(TR.sum_lcp_pairs() == pairs);
            // ノード：行きがけ順は辞書順で、親が先
            const auto pre_order = TR.preorder();
            TRIE_CHECK((int)pre_order.size() == TR.node_count() && pre_order[0] == TR.root());
            std::vector<int> seen(TR.node_count(), 0);
            Ids prev_ids;
            for (std::size_t k = 0; k < pre_order.size(); ++k) {
                const int v = pre_order[k];
                TRIE_CHECK(!seen[v] && (v == 0 || seen[TR.parent(v)]));
                seen[v] = 1;
                Ids ids;
                for (int u = v; u > 0; u = TR.parent(u)) ids.push_back(TR.label(u));
                std::reverse(ids.begin(), ids.end());
                TRIE_CHECK((int)ids.size() == TR.depth(v) && (k == 0 || prev_ids < ids));
                prev_ids = ids;
                long long pc = 0, ec = 0;
                for (const auto& [t, c] : M) pc += is_prefix(ids, t) ? c : 0, ec += t == ids ? c : 0;
                TRIE_CHECK(TR.pass_count(v) == pc && TR.end_count(v) == ec);
                if (v > 0) TRIE_CHECK(TR.child_by_id(TR.parent(v), TR.label(v)) == v);
            }
        }
        // path
        const auto p = TR.path(S);
        TRIE_CHECK(!p.empty() && p[0] == 0);
        for (std::size_t k = 1; k < p.size(); ++k) TRIE_CHECK(TR.parent(p[k]) == p[k - 1] && TR.label(p[k]) == s[k - 1]);
        TRIE_CHECK(p.size() == s.size() + 1 || TR.find(S) == -1);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    for (int it = 0; it < 60 * TEST_SCALE; ++it) {
        const int t = it % 4;
        if (t == 0) {
            const Alphabet A = Alphabet::lower();
            const int sg = 1 + (int)(rng() % 3);
            run_ops(rng, A, sg, [&](const Ids& s) { return from_ids(s, A); }, cases);
        } else if (t == 1) {
            const Alphabet A = Alphabet::of("zyx");  // 並べた順が辞書順
            run_ops(rng, A, 3, [&](const Ids& s) { return from_ids(s, A); }, cases);
        } else if (t == 2) {
            const Alphabet A = Alphabet::ints(-5, -2);  // 数の列
            run_ops(rng, A, 3, [&](const Ids& s) {
                std::vector<long long> v;
                for (int c : s) v.push_back(c - 5);
                return v;
            }, cases);
        } else {
            const Alphabet A = Alphabet::ints(2);
            run_ops(rng, A, 2, [&](const Ids& s) { return s; }, cases);
        }
    }
    // 文字列で返すもの
    {
        Trie TR;
        for (const std::string s : {"ab", "a", "abc", "b", "ab"}) TR.insert(s);
        TRIE_CHECK((TR.words() == std::vector<std::string>{"a", "ab", "ab", "abc", "b"}) && TR.kth(3) == "abc" && TR.string_of(TR.find("abc")) == "abc");
        TRIE_CHECK(TR.child(TR.root(), 'a') == TR.find("a") && TR.child(TR.root(), 'z') == -1);
        TR.clear();
        TRIE_CHECK(TR.size() == 0 && TR.node_count() == 1);
        TR.insert("a", 3100000000LL);  // 数が大きくても C(数, 2) を溢れずに
        TRIE_CHECK(TR.sum_lcp_pairs() == 4804999998450000000LL);
        ++cases;
    }
    std::cout << "[trie] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[trie] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    std::vector<std::string> vs;
    long long total = 0;
    while (total < 1000000) {
        vs.push_back(rnd_str(rng, 1 + (int)(rng() % 20), 26));
        total += (long long)vs.back().size();
    }
    auto t0 = std::chrono::steady_clock::now();
    Trie TR;
    for (const auto& s : vs) TR.insert(s);
    const long long t_ins = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = 0;
    for (const auto& s : vs) x += TR.count(s) + TR.count_prefix(s.substr(0, 3)) + TR.max_lcp(s, true);
    const long long t_q = ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < 100000; ++k) x += (long long)TR.kth(rng() % TR.size()).size() + TR.rank(vs[rng() % vs.size()]);
    const long long t_k = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += TR.sum_lcp_pairs();
    const long long t_p = ms(t0);
    std::cout << "[trie] benchmark: 文字の合計 10^6（" << vs.size() << " 個）insert " << t_ins << " ms, count + count_prefix + max_lcp " << t_q
              << " ms, kth + rank ×10^5 " << t_k << " ms, sum_lcp_pairs " << t_p << " ms, ノード " << TR.node_count() << "  (check " << x % 7 << ")\n";
}

}  // namespace trie_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace trie_demo {

{COMMON_DEMO}

static void demo() {
    banner("Trie（trie） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  insert / erase / erase_prefix / clear                      追加（cnt 個）・削除・接頭辞で全部消す\n"
              << "  count / contains / count_prefix / count_prefixes_of / size   数・接頭辞の数・接頭辞である数\n"
              << "  kth / kth_ids / rank / words / words_ids                     辞書順 k 番目・順位・全部\n"
              << "  max_lcp / sum_lcp_with / sum_lcp_pairs                       LCP の最大・和\n"
              << "  find / path / child / child_by_id / parent / depth / label / pass_count / end_count / string_of / preorder   ノード\n"
              << "  Trie TR(A) の A で文字の種類（既定は英小文字。Alphabet::upper()・of(\"01\")・ints(K) など。番号の順が辞書順）\n";

    RUN(Trie TR);
    RUN(for (const std::string s : {"apple", "app", "apply", "banana", "app"}) TR.insert(s));

    api("【1】 insert(s, cnt) / erase(s, cnt) / erase_prefix(p)", "int（ノード）/ bool / long long",
        "cnt 個入れる（0 個でもノードは作る）/ cnt 個消す（足りなければ何もせず false）/ p で始まる文字列を全部消す（消した数）", "O(|s|)");
    SHOW(TR.size(), "app は 2 個");
    SHOW(TR.erase("ap"), "入っていない → false");

    api("【2】 count(s) / contains(s) / count_prefix(p) / count_prefixes_of(s)", "long long / bool",
        "s の数 / 入っているか / p で始まる文字列の数（p 自身も）/ s の接頭辞である文字列の数（s 自身も）", "O(|s|)");
    SHOW(TR.count("app"), "");
    SHOW(TR.contains("appl"), "接頭辞だけでは入っていない");
    SHOW(TR.count_prefix("app"), "app, app, apple, apply");
    SHOW(TR.count_prefixes_of("apples"), "app, app, apple");

    api("【3】 kth(k) / rank(s) / words()", "string / long long / vector<string>", "辞書順 k 番目（0 から。同じものも数える）/ s より小さい数 / 全部を辞書順に", "O(|答え| σ) / O(|s| σ)");
    SHOW(TR.kth(2), "app, app, apple, ...");
    SHOW(TR.rank("apply"), "");
    SHOW(TR.words(), "");

    api("【4】 max_lcp(s, exclude_self) / sum_lcp_with(s) / sum_lcp_pairs()", "int / long long / long long",
        "入っている文字列との LCP の最大（exclude_self = true で s を 1 つ除く）/ 全部との LCP の和 / 全部の組（i < j）の LCP の和", "O(|s|) / O(ノード数)");
    SHOW(TR.max_lcp("applause"), "appl");
    SHOW(TR.max_lcp("banana", true), "自分を除くと共通なし");
    SHOW(TR.sum_lcp_with("apx"), "2 + 2 + 2 + 2 + 0");
    SHOW(TR.sum_lcp_pairs(), "");

    api("【5】 ノード：find(s) / path(s) / child(v, c) / parent(v) / depth(v) / label(v) / pass_count(v) / end_count(v) / string_of(v) / preorder()",
        "int / vector<int> / ...", "ノード 0 が根。find はなければ -1（数が 0 でもノードがあれば返す）。path は根から s をたどったノードの列。preorder は辞書順の行きがけ順（逆にすれば子が先。DP 用）", "O(|s|) / O(1)");
    {
        RUN(const int v = TR.find("app"));
        SHOW(TR.pass_count(v), "");
        SHOW(TR.end_count(v), "");
        SHOW(TR.depth(v), "");
        SHOW(TR.string_of(TR.child(v, 'l')), "");
        SHOW(TR.path("apex"), "a, ap まで（e がない）");
    }
    SHOW(TR.erase_prefix("app"), "app で始まる 4 個を消す");
    SHOW(TR.words(), "");

    head("問題での使い方");
    {
        Trie K;
        const std::vector<std::string> S = {"abc", "abb", "aac"};
        for (const auto& s : S) K.insert(s);
        std::vector<int> ans;
        for (const auto& s : S) ans.push_back(K.max_lcp(s, true));
        note("ABC287 E（Karuta）: 全部 insert して max_lcp(S_i, true)。入力例 1 → " + join(ans));
    }
    {
        Trie K;
        for (const std::string s : {"ab", "abc", "arc"}) K.insert(s);
        note("ABC353 E（Yet Another Sigma Problem）: 全部 insert して sum_lcp_pairs()。入力例 1 → " + to_s(K.sum_lcp_pairs()));
    }
    {
        note("ABC268 G（Random Student ID）: 自分の接頭辞の数 P = count_prefixes_of(S_i) - 1、自分で始まる数 Q = count_prefix(S_i) - 1 として");
        note("  期待値 = 1 + P + (N - 1 - P - Q) / 2（mod は math_v2 の ModInt998 で）");
        Trie K;
        const std::vector<std::string> S = {"a", "aa", "ab"};
        for (const auto& s : S) K.insert(s);
        std::vector<std::string> ans;
        for (const auto& s : S) {
            const long long P = K.count_prefixes_of(s) - 1, Q = K.count_prefix(s) - 1, N = (long long)S.size();
            const long long twice = 2 + 2 * P + (N - 1 - P - Q);  // 2 倍した値
            ans.push_back(twice % 2 ? to_s(twice) + "/2" : to_s(twice / 2));
        }
        note("  入力例 1 → " + join(ans) + "（5/2 は mod 998244353 で 499122179）");
    }
    {
        note("ABC377 G（Edit to Match）: ノードごとに best[v] = min(|S_j| - depth(v)) を持ち、S_k の path の上で |S_k| - depth(v) + best[v] の最小");
        Trie K;
        std::vector<long long> best;
        std::vector<long long> ans;
        for (const std::string s : {"snuke", "snuki", "snuuk"}) {
            const long long n = (long long)s.size();
            long long a = n;
            for (int v : K.path(s))
                if (v < (int)best.size()) a = std::min(a, n - K.depth(v) + best[v]);
            ans.push_back(a);
            K.insert(s);
            best.resize(K.node_count(), (long long)4e18);
            for (int v : K.path(s)) best[v] = std::min(best[v], n - K.depth(v));
        }
        note("  入力例 1 → " + join(ans));
    }
    {
        note("ABC403 E（Forbidden Prefix）: X は insert(X, 0) でノードを作って印を付け、erase_prefix(X)。Y は path に印がなければ insert。答えは size()");
        Trie K;
        std::vector<char> banned;
        std::vector<long long> ans;
        const std::vector<std::pair<int, std::string>> qs = {{1, "at"}, {2, "watcoder"}, {2, "atcoder"}, {1, "wa"}};
        for (const auto& [t, s] : qs) {
            if (t == 1) {
                const int v = K.insert(s, 0);
                banned.resize(K.node_count(), 0);
                banned[v] = 1;
                K.erase_prefix(s);
            } else {
                bool ok = true;
                for (int v : K.path(s)) ok = ok && !(v < (int)banned.size() && banned[v]);
                if (ok) K.insert(s);
            }
            ans.push_back(K.size());
        }
        note("  入力例 1 → " + join(ans));
    }
}

}  // namespace trie_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    trie_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    trie_test::benchmark();    // 速度計測
    trie_demo::demo();         // 使い方デモ
    return 0;
}
