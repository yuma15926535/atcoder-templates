// =============================================================================
//  ランダムテスト
// =============================================================================
namespace lyn_test {

{COMMON_TEST}

using VL = std::vector<long long>;
static bool bf_lyndon(const VL& w) {
    if (w.empty()) return false;
    for (std::size_t k = 1; k < w.size(); ++k)
        if (!(w < VL(w.begin() + k, w.end()))) return false;
    return true;
}
static int bf_period(const VL& w) {
    const int m = (int)w.size();
    for (int p = 1; p <= m; ++p) {
        bool ok = true;
        for (int i = 0; i + p < m && ok; ++i) ok = w[i] == w[i + p];
        if (ok) return p;
    }
    return m;
}

template <class S>
static void check(const S& s) {
    const VL v = str_internal::codes_of(s);
    const int n = (int)v.size();
    // 分解
    const auto f = lyndon_factorization(s);
    int pos = 0;
    VL prev;
    for (std::size_t k = 0; k < f.size(); ++k) {
        LYN_CHECK(f[k].l == pos && f[k].len() > 0);
        const VL w(v.begin() + f[k].l, v.begin() + f[k].r);
        LYN_CHECK(bf_lyndon(w) && (k == 0 || !(prev < w)));
        prev = w, pos = f[k].r;
    }
    LYN_CHECK(pos == n && is_lyndon(s) == bf_lyndon(v));
    // 各接頭辞の最小の接尾辞
    const auto ms = min_suffix_of_prefixes(s);
    LYN_CHECK((int)ms.size() == n + 1 && ms[0] == 0);
    for (int i = 1; i <= n; ++i) {
        int best = 0;
        for (int j = 1; j < i; ++j)
            if (VL(v.begin() + j, v.begin() + i) < VL(v.begin() + best, v.begin() + i)) best = j;
        LYN_CHECK(ms[i] == best);
    }
    // Lyndon 配列
    const auto lam = lyndon_array(s);
    for (int i = 0; i < n; ++i) {
        int best = 0;
        for (int l = 1; i + l <= n; ++l)
            if (bf_lyndon(VL(v.begin() + i, v.begin() + i + l))) best = l;
        LYN_CHECK(lam[i] == best);
    }
    // runs と平方
    std::vector<Repetition> runs;
    for (int l = 0; l < n; ++l)
        for (int r = l + 2; r <= n; ++r) {
            const int p = bf_period(VL(v.begin() + l, v.begin() + r));
            if (2 * p > r - l) continue;
            if (l > 0 && v[l - 1] == v[l - 1 + p]) continue;
            if (r < n && v[r] == v[r - p]) continue;
            runs.push_back({l, r, p});
        }
    LYN_CHECK(maximal_repetitions(s) == runs);
    long long sq = 0;
    for (int i = 0; i < n; ++i)
        for (int L = 1; i + 2 * L <= n; ++L) sq += std::equal(v.begin() + i, v.begin() + i + L, v.begin() + i + L);
    LYN_CHECK(count_square_substrings(s) == sq);
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 26), 1 + (int)(rng() % 3));
        check(s);
        ++cases;
    }
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        std::vector<long long> a(rng() % 20);
        for (auto& x : a) x = rng() % 3 == 0 ? rnd_ll(rng) % 3 : (long long)(rng() % 2) - 1;
        check(a);
        std::string b = rnd_str(rng, (int)(rng() % 15), 2);
        for (auto& c : b) c = (char)(c == 'a' ? 0x80 : 0x7f);  // 128 以上の char
        check(b);
        ++cases;
    }
    std::cout << "[lyndon] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[lyndon] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const std::string s = rnd_str(rng, 1000000, 2);
    auto t0 = std::chrono::steady_clock::now();
    long long x = (long long)lyndon_factorization(s).size() + min_suffix_of_prefixes(s)[777];
    const long long t_f = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += lyndon_array(s)[5];
    const long long t_a = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += (long long)maximal_repetitions(s).size();
    const long long t_r = ms(t0);
    std::string fib = "a", fb = "ab";
    while (fb.size() < 1000000) {
        std::string c = fb + fib;
        fib = fb, fb = c;
    }
    t0 = std::chrono::steady_clock::now();
    x += count_square_substrings(fb);
    const long long t_fib = ms(t0);
    std::cout << "[lyndon] benchmark: N = 10^6 で Lyndon 分解 + 各接頭辞の最小の接尾辞 " << t_f << " ms, lyndon_array " << t_a << " ms, maximal_repetitions " << t_r
              << " ms, フィボナッチ文字列（" << fb.size() << "）の平方の数 " << t_fib << " ms  (check " << x % 7 << ")\n";
}

}  // namespace lyn_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace lyn_demo {

{COMMON_DEMO}

static void demo() {
    banner("Lyndon・runs（lyndon） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  lyndon_factorization / is_lyndon / min_suffix_of_prefixes / lyndon_array   Lyndon 語\n"
              << "  maximal_repetitions / count_square_substrings                           繰り返し（runs）\n";

    api("【1】 lyndon_factorization(s) / is_lyndon(s)", "vector<StrSpan> / bool", "s を、辞書順で増えない（l1 ≥ l2 ≥ …）Lyndon 語の並びに一意に分ける", "O(N)");
    SHOW(lyndon_factorization("banana"), "b | an | an | a");
    SHOW(is_lyndon("aabab"), "");
    SHOW(is_lyndon("abab"), "ab の 2 回は Lyndon 語でない");
    note("最後の因子は s の辞書順最小の接尾辞。最小の回転は 05_lex_order の min_rotation");

    api("【2】 min_suffix_of_prefixes(s)", "vector<int>（長さ N + 1）", "r[i] = s[0, i) の辞書順最小の接尾辞の開始位置", "O(N)");
    SHOW(min_suffix_of_prefixes("abaab"), "");

    api("【3】 lyndon_array(s)", "vector<int>", "λ[i] = i から始まる最長の Lyndon 語の長さ（= 次に小さい接尾辞までの距離）", "O(N)");
    SHOW(lyndon_array("abaab"), "");

    api("【4】 maximal_repetitions(s) / count_square_substrings(s)", "vector<Repetition{l, r, period}> / long long",
        "極大な繰り返し（runs）の一覧 / xx の形の部分文字列の数（位置が違えば別）", "O(N log N)");
    SHOW(maximal_repetitions("aabaabab"), "aa・aabaaba（周期 3）・aa・abab");
    SHOW(count_square_substrings("aabaabab"), "");
    SHOW(count_square_substrings("aaaa"), "aa を 3 か所・aaaa を 1 か所");
}

}  // namespace lyn_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    lyn_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    lyn_test::benchmark();    // 速度計測
    lyn_demo::demo();         // 使い方デモ
    return 0;
}
