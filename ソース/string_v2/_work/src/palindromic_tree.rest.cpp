// =============================================================================
//  ランダムテスト
// =============================================================================
namespace pt_test {

{COMMON_TEST}

struct Mod {  // テスト用の mod 998244353
    unsigned v = 0;
    Mod() = default;
    Mod(long long x) : v((unsigned)(((x % 998244353) + 998244353) % 998244353)) {}
    Mod operator+(const Mod& o) const { return Mod((long long)v + o.v); }
    bool operator==(const Mod& o) const { return v == o.v; }
};
static bool is_pal(const std::string& s) { return std::equal(s.begin(), s.end(), s.rbegin()); }

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 30), sg);
        const int n = (int)s.size();
        const PalindromicTree pt(s);
        // 異なる回文と出現回数
        std::map<std::string, long long> occ;
        std::map<std::string, int> first_end;
        for (int r = 1; r <= n; ++r)
            for (int l = 0; l < r; ++l) {
                const std::string t = s.substr(l, r - l);
                if (!is_pal(t)) continue;
                ++occ[t];
                if (!first_end.count(t)) first_end[t] = r - 1;
            }
        PT_CHECK(pt.count_distinct() == (int)occ.size() && pt.size() == (int)occ.size() + 2 && pt.length() == n);
        long long total = 0;
        for (auto& [t, c] : occ) total += c;
        PT_CHECK(pt.count_total() == total);
        const auto oc = pt.occurrences();
        const auto pals = pt.palindromes();
        std::set<std::string> seen;
        for (int v = 2; v < pt.size(); ++v) {
            const StrSpan sp = pt.span(v);
            const std::string t = s.substr(sp.l, sp.len());
            PT_CHECK(pals[v - 2] == sp && is_pal(t) && (int)t.size() == pt.len(v) && occ.count(t) && !seen.count(t));
            seen.insert(t);
            PT_CHECK(oc[v] == occ[t] && pt.first_end(v) == first_end[t]);
            // link：最長の真の回文接尾辞 / parent：両端を除いたもの / next：両側に足したもの
            int best = 0;
            for (int k = 1; k < (int)t.size(); ++k)
                if (is_pal(t.substr(t.size() - k))) best = k;
            PT_CHECK(pt.len(pt.link(v)) == best);
            const int pv = pt.parent(v);
            PT_CHECK((int)t.size() <= 2 ? pv == 1 : (pt.len(pv) == (int)t.size() - 2 && s.substr(pt.span(pv).l, pt.len(pv)) == t.substr(1, t.size() - 2)));
            if ((int)t.size() >= 3) PT_CHECK(pt.next(pv, t[0]) == v);
        }
        // 各接頭辞の最長の回文接尾辞
        for (int i = 0; i <= n; ++i) {
            int best = 0;
            for (int k = 1; k <= i; ++k)
                if (is_pal(s.substr(i - k, k))) best = k;
            PT_CHECK(pt.len(pt.node_of_prefix(i)) == best);
        }
        // 分割
        std::vector<int> dp(n + 1, 1 << 29);
        std::vector<long long> ways(n + 1, 0);
        dp[0] = 0, ways[0] = 1;
        for (int i = 1; i <= n; ++i)
            for (int j = 0; j < i; ++j)
                if (is_pal(s.substr(j, i - j))) dp[i] = std::min(dp[i], dp[j] + 1), ways[i] += ways[j];
        PT_CHECK(pt.min_partition() == dp[n] && pt.count_partitions<long long>() == ways[n] && pt.count_partitions<Mod>() == Mod(ways[n]));
        const auto mpp = pt.min_partition_prefix();
        const auto cnp = pt.count_partitions_prefix<long long>();
        PT_CHECK((int)mpp.size() == n + 1 && (int)cnp.size() == n + 1);
        for (int i = 0; i <= n; ++i) PT_CHECK(mpp[i] == dp[i] && cnp[i] == ways[i]);
        const auto parts = pt.min_partition_spans();
        int pos = 0;
        for (const StrSpan& sp : parts) {
            PT_CHECK(sp.l == pos && is_pal(s.substr(sp.l, sp.len())));
            pos = sp.r;
        }
        PT_CHECK(pos == n && (int)parts.size() == dp[n]);
        ++cases;
    }
    // 1 文字ずつ足しても同じ・数の列
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        std::vector<int> a(rng() % 25);
        for (auto& x : a) x = (int)(rng() % 3);
        PalindromicTree pt(Alphabet::ints(3));
        for (int x : a) pt.add(x);
        std::set<std::vector<int>> st;
        for (int l = 0; l < (int)a.size(); ++l)
            for (int r = l + 1; r <= (int)a.size(); ++r) {
                std::vector<int> t(a.begin() + l, a.begin() + r);
                if (std::equal(t.begin(), t.end(), t.rbegin())) st.insert(t);
            }
        PT_CHECK(pt.count_distinct() == (int)st.size());
        ++cases;
    }
    // 問い合わせと add を混ぜる（出現回数の作り直し）
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        PalindromicTree pt;
        std::string s;
        for (int k = 0; k < 15; ++k) {
            const char c = (char)('a' + rng() % sg);
            s += c, pt.add(c);
            const auto& oc = pt.occurrences();
            long long total = 0;
            for (int v = 2; v < pt.size(); ++v) {
                const std::string t = s.substr(pt.span(v).l, pt.len(v));
                long long c2 = 0;
                for (std::size_t l = 0; l + t.size() <= s.size(); ++l) c2 += s.compare(l, t.size(), t) == 0;
                PT_CHECK(oc[v] == c2);
                total += c2;
            }
            PT_CHECK(pt.count_total() == total && pt.min_partition_prefix().back() == pt.min_partition());
        }
        ++cases;
    }
    std::cout << "[palindromic_tree] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[palindromic_tree] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::string s = rnd_str(rng, 1000000, 2);
    auto t0 = std::chrono::steady_clock::now();
    const PalindromicTree pt(s, Alphabet::of("ab"));
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = pt.count_total() + pt.count_distinct();
    const long long t_cnt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += pt.min_partition();
    const long long t_part = ms(t0);
    const std::string a = std::string(1000000, 'a');
    t0 = std::chrono::steady_clock::now();
    const PalindromicTree pa(a);
    x += pa.count_partitions<double>() > 0;
    const long long t_aa = ms(t0);
    std::cout << "[palindromic_tree] benchmark: 作成（N = 10^6, σ = 2）" << t_build << " ms, count_total " << t_cnt << " ms, min_partition " << t_part
              << " ms, 全部 a（10^6）の作成 + count_partitions " << t_aa << " ms  (check " << x % 7 << ")\n";
}

}  // namespace pt_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace pt_demo {

{COMMON_DEMO}

static void demo() {
    banner("回文木（palindromic_tree） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  add / count_distinct / count_total / occurrences / palindromes / span\n"
              << "  min_partition / min_partition_spans / min_partition_prefix / count_partitions<T> / count_partitions_prefix<T>\n"
              << "  size / length / len / link / parent / next / node_of_prefix / first_end   ノード（0 は長さ -1 の根、1 は空文字列）\n";

    RUN(const PalindromicTree pt(std::string("ababb")));

    api("【1】 count_distinct() / count_total() / occurrences() / palindromes()", "int / long long / vector<long long> / vector<StrSpan>",
        "異なる回文の数 / 回文の部分文字列の数（位置が違えば別）/ ノードごとの出現回数 / それぞれの最初の出現", "O(N)");
    SHOW(pt.count_distinct(), "a, b, aba, bab, bb");
    SHOW(pt.count_total(), "");
    SHOW(pt.palindromes(), "");
    SHOW(pt.occurrences(), "ノード 0・1 は 0");

    api("【2】 ノード：len(v) / link(v) / parent(v) / next(v, c) / node_of_prefix(i)", "int",
        "長さ / 最長の真の回文接尾辞 / 両端を除いた回文 / 両側に c を足した回文（なければ -1）/ s[0, i) の最長の回文接尾辞", "O(1)");
    {
        RUN(const int v = pt.node_of_prefix(4));
        SHOW(pt.len(v), "\"abab\" の最長の回文接尾辞は bab");
        SHOW(pt.len(pt.link(v)), "bab の最長の真の回文接尾辞は b");
        SHOW(pt.len(pt.parent(v)), "bab の両端を除くと a");
        SHOW(pt.next(pt.parent(v), 'b') == v, "");
    }

    api("【3】 min_partition() / min_partition_spans() / count_partitions<T>() / min_partition_prefix() / count_partitions_prefix<T>()",
        "int / vector<StrSpan> / T / vector<int> / vector<T>",
        "回文に分ける最小の個数 / その分け方 / 分け方の数（N = 10^6 でも速い。N ≤ 5000 なら 04_palindrome の palindrome_partition_min でもよい）", "O(N log N)");
    SHOW(pt.min_partition(), "aba + bb");
    SHOW(pt.min_partition_spans(), "");
    SHOW(pt.count_partitions<long long>(), "");
    SHOW(pt.min_partition_prefix(), "各接頭辞 s[0, i) の最小の個数（count_partitions_prefix<T>() は各接頭辞の分け方の数）");

    head("問題での使い方：ABC237 Ex（Hakata）");
    note("回文どうしが含まれない最大の選び方 → 回文のノードを頂点に、v → link(v)・v → parent(v) をたどれる関係（推移閉包）で");
    note("「どれも含まない」最大の集合（反鎖）= 頂点の数 − 二部グラフの最大マッチング（graph_grid_v2 08_bipartite_matching）。入力例:");
    for (const std::string s : {"ababb", "xyz", "xxxxxxxxxx"}) {
        const PalindromicTree t(s);
        const int m = t.size();
        std::vector<std::vector<char>> reach(m, std::vector<char>(m, 0));
        for (int v = 2; v < m; ++v) {  // v から link・parent をたどって着く回文（ノード 2 以上）
            std::vector<int> st = {v};
            while (!st.empty()) {
                const int u = st.back();
                st.pop_back();
                for (int w : {t.link(u), t.parent(u)})
                    if (w >= 2 && !reach[v][w]) reach[v][w] = 1, st.push_back(w);
            }
        }
        // 二部マッチング（左 v → 右 w：v が w を含む）
        std::vector<int> mt(m, -1);
        int match = 0;
        for (int v = 2; v < m; ++v) {
            std::vector<char> used(m, 0);
            auto dfs = [&](auto&& self, int x) -> bool {
                for (int w = 2; w < m; ++w)
                    if (reach[x][w] && !used[w]) {
                        used[w] = 1;
                        if (mt[w] < 0 || self(self, mt[w])) return mt[w] = x, true;
                    }
                return false;
            };
            match += dfs(dfs, v);
        }
        note("  " + s + " → " + to_s(t.count_distinct() - match));
    }
}

}  // namespace pt_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    pt_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    pt_test::benchmark();    // 速度計測
    pt_demo::demo();         // 使い方デモ
    return 0;
}
