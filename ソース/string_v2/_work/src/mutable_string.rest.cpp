// =============================================================================
//  ランダムテスト
// =============================================================================
namespace ms_test {

{COMMON_TEST}

// 愚直：今の列を持ち、毎回その場で計算する
template <class V>
static void check_all(std::mt19937_64& rng, const MutableString& ms, const V& cur, const Alphabet& A) {
    const int n = (int)cur.size();
    MS_CHECK(ms.size() == n);
    std::vector<long long> vals(cur.begin(), cur.end());
    for (auto& x : vals) x = str_internal::code_of((typename V::value_type)x);
    MS_CHECK(ms.values() == vals);
    for (int t = 0; t < 30; ++t) {
        int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
        if (l > r) std::swap(l, r);
        const V sub(cur.begin() + l, cur.begin() + r);
        V rev(sub.rbegin(), sub.rend());
        MS_CHECK(ms.hash(l, r) == hash_of(sub) && ms.reversed_hash(l, r) == hash_of(rev));
        MS_CHECK(ms.is_palindrome(l, r) == (sub == rev));
        int l2 = (int)(rng() % (n + 1)), r2 = (int)(rng() % (n + 1));
        if (l2 > r2) std::swap(l2, r2);
        const V sub2(cur.begin() + l2, cur.begin() + r2);
        MS_CHECK(ms.same(l, r, l2, r2) == (sub == sub2));
        MS_CHECK(ms.compare(l, r, l2, r2) == str_internal::lex_compare(sub, sub2));
        int k = 0;
        while (l + k < n && l2 + k < n && cur[l + k] == cur[l2 + k]) ++k;
        MS_CHECK(ms.lcp(l, l2) == k);
        // 文字の数・並び
        std::vector<long long> cnt(A.size(), 0);
        for (auto c : sub) ++cnt[A.id(c)];
        MS_CHECK(ms.counts(l, r) == cnt);
        const int c = (int)(rng() % A.size());
        MS_CHECK(ms.count(l, r, (typename V::value_type)A.value(c)) == cnt[c]);
        bool sorted = true;
        for (int i = l; i + 1 < r; ++i) sorted &= A.id(cur[i]) <= A.id(cur[i + 1]);
        MS_CHECK(ms.is_sorted(l, r) == sorted);
        if (r > l) {
            std::vector<int> ids;
            for (auto x : sub) ids.push_back(A.id(x));
            std::sort(ids.begin(), ids.end());
            const int kk = (int)(rng() % (r - l));
            MS_CHECK(ms.kth_char(l, r, kk) == A.value(ids[kk]));
        }
    }
}

static void test_string(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const Alphabet A = rng() % 2 ? Alphabet::lower() : Alphabet::of("cba");
        std::string s = rnd_str_shaped(rng, (int)(rng() % 30), sg);
        MutableString ms(s, A);
        const bool lazy_first = rng() % 2;  // 表を作る前に書き換える場合も
        for (int q = 0; q < 60; ++q) {
            if (!s.empty() && (q < 10 ? lazy_first : rng() % 2)) {
                const int i = (int)(rng() % s.size());
                const char c = (char)('a' + rng() % sg);
                s[i] = c, ms.set(i, c);
            }
            if (q >= 10 || !lazy_first) check_all(rng, ms, s, A);
        }
        MS_CHECK(ms.to_string() == s);
        ++cases;
    }
}

static void test_values(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 200 * TEST_SCALE; ++it) {
        const Alphabet A = Alphabet::ints(-2, 2);
        std::vector<long long> v(rng() % 25);
        for (auto& x : v) x = (long long)(rng() % 4) - 2;
        MutableString ms(v, A);
        for (int q = 0; q < 40; ++q) {
            if (!v.empty()) {
                const int i = (int)(rng() % v.size());
                v[i] = (long long)(rng() % 4) - 2, ms.set(i, v[i]);
            }
            check_all(rng, ms, v, A);
        }
        ++cases;
    }
    // 大きな値（ハッシュだけ使う。Alphabet の範囲は確かめない）
    for (int it = 0; it < 200 * TEST_SCALE; ++it) {
        std::vector<long long> v(rng() % 25);
        for (auto& x : v) x = rnd_ll(rng);
        MutableString ms(v);
        for (int q = 0; q < 40; ++q) {
            if (!v.empty() && rng() % 2) {
                const int i = (int)(rng() % v.size());
                v[i] = rng() % 2 ? rnd_ll(rng) : v[rng() % v.size()], ms.set(i, v[i]);
            }
            int l = (int)(rng() % (v.size() + 1)), r = (int)(rng() % (v.size() + 1));
            if (l > r) std::swap(l, r);
            const std::vector<long long> sub(v.begin() + l, v.begin() + r);
            MS_CHECK(ms.hash(l, r) == hash_of(sub) && ms.is_palindrome(l, r) == std::equal(sub.begin(), sub.end(), sub.rbegin()));
        }
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_string(rng, cases);
    test_values(rng, cases);
    {  // RollingHash と同じ値
        const std::string s = "abracadabra";
        MutableString ms(s);
        const RollingHash H(s);
        MS_CHECK(ms.hash(2, 9) == H.get(2, 9) && ms.hash(0, 11) == H.all());
        ++cases;
    }
    std::cout << "[mutable_string] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[mutable_string] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 1000000, Q = 1000000;
    const std::string s = rnd_str(rng, N, 2);
    auto t0 = std::chrono::steady_clock::now();
    MutableString m(s);
    const long long t_build = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long x = 0;
    for (int q = 0; q < Q; ++q) {
        const int i = (int)(rng() % N);
        if (q % 2) m.set(i, (char)('a' + rng() % 2));
        else x += m.is_palindrome(i, std::min(N, i + (int)(rng() % 20)));
    }
    const long long t_pal = ms(t0);
    t0 = std::chrono::steady_clock::now();
    for (int q = 0; q < 100000; ++q) x += m.lcp((int)(rng() % N), (int)(rng() % N));
    const long long t_lcp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    MutableString m2(rnd_str(rng, 100000, 26));
    for (int q = 0; q < 100000; ++q) {
        const int i = (int)(rng() % 100000);
        if (q % 2) m2.set(i, (char)('a' + rng() % 26));
        else x += m2.is_sorted(i, std::min(100000, i + 3)) + m2.count(0, i, 'c');
    }
    const long long t_cnt = ms(t0);
    std::cout << "[mutable_string] benchmark: 作成（N = 10^6）" << t_build << " ms, set + is_palindrome 10^6 回 " << t_pal << " ms, lcp × 10^5 " << t_lcp
              << " ms, set + is_sorted + count 10^5 回（N = 10^5、表の作成込み）" << t_cnt << " ms  (check " << x % 7 << ")\n";
}

}  // namespace ms_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace ms_demo {

{COMMON_DEMO}

static void demo() {
    banner("書き換えのある文字列（mutable_string） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  set / get / hash / reversed_hash / same / is_palindrome / lcp / compare   書き換え・ハッシュ・比較・回文\n"
              << "  count / counts / is_sorted / kth_char                                    文字の数・並んでいるか・k 番目（表は初めて使うときに作る）\n";

    RUN(MutableString ms(std::string("abcbacb")));

    api("【1】 set(i, c) / get(i) / is_palindrome(l, r)", "なし / long long / bool", "1 文字の書き換えと、区間 [l, r) が回文か（L 文字目〜R 文字目は (L - 1, R)）", "O(log N)");
    {
        note("ABC331 F（Palindrome Query）の入力例 1:");
        std::string out;
        out += ms.is_palindrome(0, 5) ? "Yes " : "No ";
        out += ms.is_palindrome(3, 7) ? "Yes " : "No ";
        out += ms.is_palindrome(1, 2) ? "Yes " : "No ";
        ms.set(4, 'c');
        out += ms.is_palindrome(0, 5) ? "Yes " : "No ";
        out += ms.is_palindrome(3, 7) ? "Yes " : "No ";
        ms.set(3, 'c');
        out += ms.is_palindrome(2, 6) ? "Yes " : "No ";
        note("  → " + out);
        SHOW(ms.to_string(), "");
        SHOW((char)ms.get(4), "");
    }

    api("【2】 hash(l, r) / reversed_hash(l, r) / same(l1, r1, l2, r2)", "StrHash / StrHash / bool",
        "区間のハッシュ（01_rolling_hash の RollingHash::get・hash_of と同じ値）/ 反転したもののハッシュ / 区間どうしが等しいか", "O(log N)");
    SHOW(ms.same(1, 3, 4, 6), "\"bc\" と \"cc\"");
    SHOW(ms.hash(0, 3) == hash_of("abc"), "ほかの文字列のハッシュとも比べられる");
    SHOW(ms.reversed_hash(0, 3) == hash_of("cba"), "");

    api("【3】 lcp(i, j) / compare(l1, r1, l2, r2)", "int / int（-1 / 0 / 1）", "書き換えたあとの最長共通接頭辞・辞書順（値の大小で）", "O(log^2 N)");
    SHOW(ms.lcp(2, 3), "");
    SHOW(ms.compare(0, 2, 5, 7), "\"ab\" と \"cb\"");

    api("【4】 count(l, r, c) / counts(l, r) / is_sorted(l, r) / kth_char(l, r, k)", "int / vector<int> / bool / long long",
        "区間の文字 c の数・番号ごとの数・広義単調増加か（Alphabet の番号の順）・並べ替えた k 番目（0 から）", "O(log N)・O(σ log N)");
    {
        note("ABC285 F（Substring of Sorted String）: s[l, r) が並べ替えた S の部分文字列 ⇔ is_sorted(l, r) かつ");
        note("  両端の文字より真ん中の文字（s[l] < c < s[r-1]）がすべて区間の中にある（count(l, r, c) == count(0, N, c)）");
        MutableString t(std::string("abcdcf"));
        auto ok = [&](int l, int r) {
            if (!t.is_sorted(l, r)) return false;
            const int a = t.alphabet().id(t.get(l)), b = t.alphabet().id(t.get(r - 1));
            for (int c = a + 1; c < b; ++c)
                if (t.count(l, r, t.alphabet().ch(c)) != t.count(0, t.size(), t.alphabet().ch(c))) return false;
            return true;
        };
        std::string out;
        out += ok(0, 3) ? "Yes " : "No ";
        out += ok(1, 6) ? "Yes " : "No ";
        t.set(4, 'e');
        out += ok(1, 6) ? "Yes " : "No ";
        note("  入力例 1 → " + out);
        SHOW(t.counts(0, 6), "a〜z の数");
        SHOW((char)t.kth_char(0, 6, 3), "abcdef の 3 番目");
    }
}

}  // namespace ms_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    ms_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    ms_test::benchmark();    // 速度計測
    ms_demo::demo();         // 使い方デモ
    return 0;
}
