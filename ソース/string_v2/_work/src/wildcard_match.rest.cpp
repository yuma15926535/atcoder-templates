// =============================================================================
//  ランダムテスト
// =============================================================================
namespace wc_test {

{COMMON_TEST}

template <class V>
static std::vector<int> bf_mismatch(const V& t, const V& p, bool has_wild, long long wild = 0) {  // 愚直（has_wild = false なら wild なし）
    const int n = (int)t.size(), m = (int)p.size();
    std::vector<int> r;
    if (m > n) return r;
    for (int i = 0; i + m <= n; ++i) {
        int c = 0;
        for (int j = 0; j < m; ++j) {
            const long long a = str_internal::code_of(t[i + j]), b = str_internal::code_of(p[j]);
            if (!(has_wild && (a == wild || b == wild)) && a != b) ++c;
        }
        r.push_back(c);
    }
    return r;
}

static void test_string(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        std::string t = rnd_str_shaped(rng, (int)(rng() % 40), sg), p = rnd_str_shaped(rng, (int)(rng() % 8), sg);
        const int wt = (int)(rng() % 3);  // wild を入れる場所（0 なし・1 pat だけ・2 両方）
        if (wt >= 1)
            for (auto& c : p)
                if (rng() % 4 == 0) c = '?';
        if (wt == 2)
            for (auto& c : t)
                if (rng() % 6 == 0) c = '?';
        // 不一致の数
        const auto bf = bf_mismatch(t, p, true, '?');
        WC_CHECK(mismatch_counts(t, p, '?') == bf);
        const auto bf0 = bf_mismatch(t, p, false);
        WC_CHECK(mismatch_counts(t, p) == bf0);
        const int k = (int)(rng() % 3);
        std::vector<int> ek, ew;
        for (int i = 0; i < (int)bf.size(); ++i) {
            if (bf[i] <= k) ek.push_back(i);
            if (bf[i] == 0) ew.push_back(i);
        }
        WC_CHECK(match_with_k_mismatches(t, p, k, '?') == ek);
        WC_CHECK(wildcard_match(t, p) == ew && wildcard_match(t, p, '?') == ew);
        ++cases;
    }
}

static void test_values(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1000 * TEST_SCALE; ++it) {
        std::vector<long long> t(rng() % 30), p(rng() % 6);
        const bool big = rng() % 2;
        for (auto& x : t) x = big ? rnd_ll(rng) % 5 * 1000000000000LL : (long long)(rng() % 4) - 1;
        for (auto& x : p) x = big ? rnd_ll(rng) % 5 * 1000000000000LL : (long long)(rng() % 4) - 1;
        if (!big && !p.empty() && rng() % 2) p[rng() % p.size()] = -1;  // -1 を wild に
        WC_CHECK(mismatch_counts(t, p, -1) == bf_mismatch(t, p, true, -1) && mismatch_counts(t, p) == bf_mismatch(t, p, false));
        std::vector<int> ew;
        const auto bf = bf_mismatch(t, p, true, -1);
        for (int i = 0; i < (int)bf.size(); ++i)
            if (bf[i] == 0) ew.push_back(i);
        WC_CHECK(wildcard_match(t, p, -1) == ew);
        ++cases;
    }
    // 端の値（LLONG_MIN・LLONG_MAX）・128 以上の char を wild に・pat にしかない値
    const long long LO = std::numeric_limits<long long>::min(), HI = std::numeric_limits<long long>::max();
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        std::vector<long long> t(rng() % 20), p(rng() % 5);
        const long long pool[5] = {LO, HI, 0, -1, (long long)(rng() % 1000) + 5};
        for (auto& x : t) x = pool[rng() % 4];
        for (auto& x : p) x = pool[rng() % 5];
        const long long w = pool[rng() % 4];
        WC_CHECK(mismatch_counts(t, p) == bf_mismatch(t, p, false) && mismatch_counts(t, p, w) == bf_mismatch(t, p, true, w));
        std::vector<int> ew;
        const auto bf = bf_mismatch(t, p, true, w);
        for (int i = 0; i < (int)bf.size(); ++i)
            if (bf[i] == 0) ew.push_back(i);
        WC_CHECK(wildcard_match(t, p, w) == ew && match_with_k_mismatches(t, p, 0, w) == ew);
        std::string a = rnd_str(rng, (int)(rng() % 20), 2), b = rnd_str(rng, (int)(rng() % 5), 3);
        for (auto& c : a)
            if (rng() % 4 == 0) c = '\xff';
        for (auto& c : b)
            if (rng() % 3 == 0) c = '\xff';
        const auto bfs = bf_mismatch(a, b, true, 255);
        std::vector<int> es;
        for (int i = 0; i < (int)bfs.size(); ++i)
            if (bfs[i] == 0) es.push_back(i);
        WC_CHECK(mismatch_counts(a, b, '\xff') == bfs && wildcard_match(a, b, '\xff') == es);
        ++cases;
    }
    // 長い列（NTT の大きさが変わる所）
    for (int it = 0; it < 20 * TEST_SCALE; ++it) {
        const std::string t = rnd_str(rng, 500 + (int)(rng() % 3000), 2), p = rnd_str(rng, 1 + (int)(rng() % 300), 2);
        WC_CHECK(mismatch_counts(t, p) == bf_mismatch(t, p, false));
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_string(rng, cases);
    test_values(rng, cases);
    std::cout << "[wildcard_match] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[wildcard_match] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::string s = rnd_str(rng, 1000000, 2, '0'), t = rnd_str(rng, 500000, 2, '0');
    auto t0 = std::chrono::steady_clock::now();
    const auto mc = mismatch_counts(s, t);
    const long long t_mc = ms(t0);
    std::string a = rnd_str(rng, 500000, 26), b = rnd_str(rng, 200000, 26);
    for (auto& c : b)
        if (rng() % 3 == 0) c = '?';
    t0 = std::chrono::steady_clock::now();
    const auto wm = wildcard_match(a, b);
    const long long t_wm = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto mc26 = mismatch_counts(a, b, '?');
    const long long t_mc26 = ms(t0);
    std::cout << "[wildcard_match] benchmark: mismatch_counts（σ = 2、10^6 と 5×10^5）" << t_mc << " ms, wildcard_match（σ = 26、5×10^5 と 2×10^5）" << t_wm
              << " ms, mismatch_counts（同じ大きさ・σ = 26）" << t_mc26 << " ms  (check " << (mc[0] + wm.size() + mc26[0]) % 7 << ")\n";
}

}  // namespace wc_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace wc_demo {

{COMMON_DEMO}

static void demo() {
    banner("ワイルドカード・不一致の数（wildcard_match） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  mismatch_counts / match_with_k_mismatches / wildcard_match\n"
              << "  ずらし i は pat を text[i, i + M) に重ねたもの（0〜N - M）\n";

    api("【1】 mismatch_counts(text, pat) / mismatch_counts(text, pat, wild)", "vector<int>（長さ N - M + 1）",
        "ずらしごとに文字が違う数 / wild つき：両方とも wild でない位置で違う数（σ' は text と pat の両方に現れる値の数）", "O(σ' N log N)");
    SHOW(mismatch_counts("0001", "101"), "ABC196 F（Substring 2）の入力例 1 → 最小は 1");
    SHOW(mismatch_counts("abcab", "a?b", '?'), "? はどの文字とも一致");
    {
        const auto mc = mismatch_counts("10101000010011011110", "0010011111");
        SHOW(*std::min_element(mc.begin(), mc.end()), "ABC196 F の入力例 3");
    }

    api("【2】 match_with_k_mismatches(text, pat, k) / (text, pat, k, wild)", "vector<int>", "不一致が k 個以下のずらしの一覧", "同上");
    SHOW(match_with_k_mismatches("abcabd", "abd", 1), "");

    api("【3】 wildcard_match(text, pat) / wildcard_match(text, pat, wild)", "vector<int>",
        "wild（文字列で省くと '?'）を任意の 1 文字として一致するずらしの一覧（両方の wild が使える。乱数を使う。違うのに一致とする確率は 1 か所あたり約 10^-17）", "O(N log N)");
    SHOW(wildcard_match("abcab?b", "a?b"), "text の ? も任意の 1 文字");
    SHOW(wildcard_match(std::vector<int>{3, 1, 4, 1, 5}, std::vector<int>{0, 1}, 0), "数の列では wild の値（ここでは 0）を必ず渡す（省くとコンパイルエラー）");
}

}  // namespace wc_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    wc_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    wc_test::benchmark();    // 速度計測
    wc_demo::demo();         // 使い方デモ
    return 0;
}
