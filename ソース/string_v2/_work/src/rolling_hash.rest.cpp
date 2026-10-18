// =============================================================================
//  ランダムテスト
// =============================================================================
namespace rh_test {

{COMMON_TEST}

static int bf_lcp(const std::string& a, int i, const std::string& b, int j) {
    int k = 0;
    while (i + k < (int)a.size() && j + k < (int)b.size() && a[i + k] == b[j + k]) ++k;
    return k;
}
static int bf_cmp(const std::string& a, const std::string& b) { return a < b ? -1 : (a == b ? 0 : 1); }

static void test_1d(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int sg = 1 + (int)(rng() % 3);
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 14), sg), t = rnd_str_shaped(rng, (int)(rng() % 14), sg);
        const RollingHash H(s), T(t);
        const int n = (int)s.size(), m = (int)t.size();
        RH_CHECK(H.size() == n && H.all() == hash_of(s) && T.all() == hash_of(t));
        for (int l1 = 0; l1 <= n; ++l1)
            for (int r1 = l1; r1 <= n; ++r1) {
                const std::string a = s.substr(l1, r1 - l1);
                RH_CHECK(H.get(l1, r1) == hash_of(a) && H.get_len(l1, r1 - l1) == H.get(l1, r1));
                for (int mid = l1; mid <= r1; ++mid) RH_CHECK(H.get(l1, mid) + H.get(mid, r1) == H.get(l1, r1));
                for (int l2 = 0; l2 <= m; ++l2)
                    for (int r2 = l2; r2 <= m; ++r2) {
                        const std::string b = t.substr(l2, r2 - l2);
                        RH_CHECK((H.get(l1, r1) == T.get(l2, r2)) == (a == b));
                        RH_CHECK(H.compare(l1, r1, T, l2, r2) == bf_cmp(a, b));
                        ++cases;
                    }
                for (int l2 = 0; l2 <= n; ++l2) {
                    const int r2 = l2 + (int)(rng() % (n - l2 + 1));
                    const std::string b = s.substr(l2, r2 - l2);
                    RH_CHECK(H.same(l1, r1, l2, r2) == (a == b) && H.compare(l1, r1, l2, r2) == bf_cmp(a, b));
                }
            }
        for (int i = 0; i <= n; ++i) {
            for (int j = 0; j <= n; ++j) {
                RH_CHECK(H.lcp(i, j) == bf_lcp(s, i, s, j));
                int k = 0;
                while (k < i && k < j && s[i - 1 - k] == s[j - 1 - k]) ++k;
                RH_CHECK(H.common_suffix(i, j) == k);
            }
            for (int j = 0; j <= m; ++j) {
                RH_CHECK(H.lcp(i, T, j) == bf_lcp(s, i, t, j));
                int k = 0;
                while (k < i && k < j && s[i - 1 - k] == t[j - 1 - k]) ++k;
                RH_CHECK(H.common_suffix(i, T, j) == k);
            }
        }
        // 出現位置
        std::vector<int> ex;
        for (int i = 0; i + m <= n; ++i)
            if (s.compare(i, m, t) == 0) ex.push_back(i);
        RH_CHECK(H.find_all(t) == ex && H.find_all(hash_of(t)) == ex && H.count(t) == (long long)ex.size());
        // 連結・繰り返し・1 文字
        RH_CHECK(hash_of(s) + hash_of(t) == hash_of(s + t) && hash_of(s + t) == hash_of(s) + T.all());
        StrHash acc = empty_hash();
        for (char c : s) acc += hash_of_char(c);
        RH_CHECK(acc == hash_of(s));
        const long long k = rng() % 6;
        RH_CHECK(hash_of(s).repeat(k) == hash_of(repeat_string(s, k)));
        RH_CHECK(hash_of(s).repeat(0) == empty_hash());
        // 途中に挿入した文字列のハッシュ
        const int p = (int)(rng() % (n + 1));
        RH_CHECK(H.prefix(p) + T.all() + H.suffix(p) == hash_of(s.substr(0, p) + t + s.substr(p)));
    }
    // 巨大な繰り返し（長さだけでなく値も：(a^k)^2 = a^(2k)）
    const StrHash a = hash_of("abc");
    const StrHash big = a.repeat(1000000000000LL);
    RH_CHECK(big.len == 3000000000000LL && big + big == a.repeat(2000000000000LL) && big.repeat(3) == a.repeat(3000000000000LL));
}

static void test_values(std::mt19937_64& rng, int& cases) {
    // 数の列：大きな値・負の値・2^61 - 1 の倍数ずれも区別する
    const long long M = (long long)((1ULL << 61) - 1);
    const std::vector<long long> sp = {0, 1, -1, M, M + 1, 2 * M, -M, (1LL << 60) - 1, 1LL << 60, std::numeric_limits<long long>::max(), std::numeric_limits<long long>::min()};
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int n = 1 + (int)(rng() % 6);
        std::vector<long long> a(n), b(n);
        for (auto& x : a) x = rng() % 2 ? sp[rng() % sp.size()] : rnd_ll(rng);
        b = a;
        if (rng() % 2) b[rng() % n] = rng() % 2 ? sp[rng() % sp.size()] : rnd_ll(rng);
        RH_CHECK((hash_of(a) == hash_of(b)) == (a == b));
        const RollingHash A(a), B(b);
        int k = 0;
        while (k < n && a[k] == b[k]) ++k;
        RH_CHECK(A.lcp(0, B, 0) == k);
        RH_CHECK(A.compare(0, n, B, 0, n) == (a < b ? -1 : (a == b ? 0 : 1)));
        ++cases;
    }
    // string と vector<int> は同じ値なら同じハッシュ
    const std::string s = "hello";
    RH_CHECK(hash_of(s) == hash_of(std::vector<int>(s.begin(), s.end())) && hash_of("hello") == hash_of(s));
    // 128 以上の char
    std::string hi = "a";
    hi += (char)200;
    RH_CHECK(RollingHash(hi).compare(0, 1, 1, 2) == -1);
}

static void test_distinct(std::mt19937_64& rng, int& cases) {
    // 異なる部分文字列の数をハッシュと set<string> で（衝突しないことの確認）
    for (int it = 0; it < 6 * TEST_SCALE; ++it) {
        const std::string s = rnd_str_shaped(rng, 250, 2 + (int)(rng() % 2));
        const RollingHash H(s);
        std::unordered_set<StrHash> hs;
        std::set<std::string> ss;
        for (int l = 0; l < (int)s.size(); ++l)
            for (int r = l + 1; r <= (int)s.size(); ++r) hs.insert(H.get(l, r)), ss.insert(s.substr(l, r - l));
        RH_CHECK(hs.size() == ss.size());
        cases += (int)ss.size();
    }
    std::unordered_map<StrHash, int> mp;
    mp[hash_of("ab")] = 3;
    RH_CHECK(mp[hash_of("a") + hash_of("b")] == 3);
}

static void test_2d(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int H = (int)(rng() % 6), W = H ? 1 + (int)(rng() % 6) : 0, sg = 1 + (int)(rng() % 2);
        std::vector<std::string> g(H);
        for (auto& row : g) row = rnd_str(rng, W, sg);
        const RollingHash2D G(g);
        RH_CHECK(G.rows() == H && G.cols() == W);
        auto sub = [&](int i1, int j1, int i2, int j2) {
            std::vector<std::string> r;
            for (int i = i1; i < i2; ++i) r.push_back(g[i].substr(j1, j2 - j1));
            return r;
        };
        for (int t = 0; t < 30; ++t) {
            const int i1 = (int)(rng() % (H + 1)), i2 = i1 + (int)(rng() % (H - i1 + 1));
            const int j1 = (int)(rng() % (W + 1)), j2 = j1 + (int)(rng() % (W - j1 + 1));
            const int k1 = (int)(rng() % (H - (i2 - i1) + 1)), l1 = (int)(rng() % (W - (j2 - j1) + 1));
            const auto a = sub(i1, j1, i2, j2), b = sub(k1, l1, k1 + (i2 - i1), l1 + (j2 - j1));
            RH_CHECK((G.get(i1, j1, i2, j2) == G.get(k1, l1, k1 + i2 - i1, l1 + j2 - j1)) == (a == b || i1 == i2 || j1 == j2));
            ++cases;
        }
        // 模様探し
        const int h = (int)(rng() % 3) + 1, w = (int)(rng() % 3) + 1;
        std::vector<std::string> pat(h);
        for (auto& row : pat) row = rnd_str(rng, w, sg);
        std::vector<std::pair<int, int>> ex;
        for (int i = 0; i + h <= H; ++i)
            for (int j = 0; j + w <= W; ++j)
                if (sub(i, j, i + h, j + w) == pat) ex.push_back({i, j});
        RH_CHECK(G.find_all(pat) == ex);
        std::vector<std::vector<int>> gi(H, std::vector<int>(W));
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) gi[i][j] = g[i][j];
        RH_CHECK(RollingHash2D(gi).all() == G.all());
    }
}

static void test_changes(std::mt19937_64& rng, int& cases) {  // center_palindrome_lengths_with_changes
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 25), k = (int)(rng() % 4);
        const std::string s = rnd_str_shaped(rng, n, 1 + (int)(rng() % 3));
        const auto L = center_palindrome_lengths_with_changes(s, k);
        RH_CHECK((int)L.size() == std::max(0, 2 * n - 1));
        for (int c = 0; c < 2 * n - 1; ++c) {
            int best = 0;
            for (int len = c % 2 == 0 ? 1 : 2; len <= n; len += 2) {  // 中心 c の長さ len の区間
                const int l = (c + 1 - len) / 2;
                if (l < 0 || l + len > n) break;
                int mis = 0;
                for (int t = 0; t < len / 2; ++t) mis += s[l + t] != s[l + len - 1 - t];
                if (mis <= k) best = len;
                else break;
            }
            RH_CHECK(L[c] == best);
        }
        if (k == 0) {  // 04_palindrome の center_palindrome_lengths と同じ（愚直に）
            for (int c = 0; c < 2 * n - 1; ++c) {
                int len = c % 2 == 0 ? 1 : 0, l = c % 2 == 0 ? c / 2 - 1 : (c - 1) / 2, r = c % 2 == 0 ? c / 2 + 1 : (c + 1) / 2;
                while (l >= 0 && r < n && s[l] == s[r]) --l, ++r, len += 2;
                RH_CHECK(L[c] == len);
            }
        }
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_1d(rng, cases);
    test_values(rng, cases);
    test_distinct(rng, cases);
    test_2d(rng, cases);
    test_changes(rng, cases);
    std::cout << "[rolling_hash] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[rolling_hash] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 10000000;
    const std::string s = rnd_str(rng, n, 26);
    auto t0 = std::chrono::steady_clock::now();
    const RollingHash H(s);
    const long long t_build = ms(t0);
    std::vector<int> ls(n), rs(n);
    for (int i = 0; i < n; ++i) {
        ls[i] = (int)(rng() % (n + 1)), rs[i] = (int)(rng() % (n + 1));
        if (ls[i] > rs[i]) std::swap(ls[i], rs[i]);
    }
    t0 = std::chrono::steady_clock::now();
    unsigned long long x = 0;
    for (int i = 0; i < n; ++i) x += H.get(ls[i], rs[i]).value;
    const long long t_get = ms(t0);
    const std::string s2 = rnd_str_shaped(rng, 1000000, 2);
    const RollingHash H2(s2);
    t0 = std::chrono::steady_clock::now();
    long long y = 0;
    for (int i = 0; i < 1000000; ++i) y += H2.lcp((int)(rng() % 1000001), (int)(rng() % 1000001));
    const long long t_lcp = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const auto occ = H.find_all(s.substr(5000000, 3));
    const long long t_find = ms(t0);
    std::vector<std::string> g(2000, std::string(2000, 'a'));
    for (auto& row : g)
        for (char& c : row) c = (char)('a' + rng() % 2);
    t0 = std::chrono::steady_clock::now();
    const RollingHash2D G(g);
    const auto f2 = G.find_all(std::vector<std::string>{"ab", "ba"});
    const long long t_2d = ms(t0);
    std::cout << "[rolling_hash] benchmark: 構築（N = 10^7）" << t_build << " ms, get × 10^7 " << t_get << " ms, lcp × 10^6（N = 10^6）" << t_lcp
              << " ms, find_all（N = 10^7）" << t_find << " ms, 2D 構築 + find_all（2000×2000）" << t_2d << " ms  (check "
              << (x + (unsigned long long)y + occ.size() + f2.size()) % 7 << ")\n";
}

}  // namespace rh_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace rh_demo {

{COMMON_DEMO}

static void demo() {
    banner("ローリングハッシュ（rolling_hash） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  RollingHash H(s)：get / get_len / prefix / suffix / all / same / lcp / common_suffix / compare / find_all / count\n"
              << "  StrHash{value, len}：+（連結）/ repeat(k) / == / unordered_map のキーに\n"
              << "  hash_of(s) / hash_of_char(c) / empty_hash()\n"
              << "  RollingHash2D G(grid)：get(i1, j1, i2, j2) / find_all(pat) / all\n"
              << "  center_palindrome_lengths_with_changes(s, k)   k 文字まで書き換えてよい回文（中心ごとの最長）\n"
              << "  基数は実行ごとの乱数で全部の文字列で共通なので、別の文字列の部分文字列どうしも比べられる\n";

    api("【1】 RollingHash H(s) / H.get(l, r) / H.same(l1, r1, l2, r2)", "StrHash / bool",
        "s[l, r) のハッシュ。同じ長さで値が同じなら（ほぼ確実に）同じ文字列。何度でも O(1)", "構築 O(N)・取得 O(1)");
    RUN(const std::string S = "abcabcab");
    RUN(const RollingHash H(S));
    SHOW(H.same(0, 3, 3, 6), "\"abc\" と \"abc\"");
    SHOW(H.same(0, 3, 5, 8), "\"abc\" と \"cab\"");
    SHOW(H.get(1, 4) == hash_of("bca"), "別に作ったハッシュとも比べられる");
    note("問題の「L 文字目から R 文字目」は H.get(L - 1, R)。");

    api("【2】 H.lcp(i, j) / H.common_suffix(i, j) / H.compare(l1, r1, l2, r2)", "int",
        "s[i..) と s[j..) の最長共通接頭辞・s[0, i) と s[0, j) の最長共通接尾辞・部分文字列の辞書順（-1 / 0 / 1）。二分探索", "O(log N)");
    SHOW(H.lcp(0, 3), "\"abcabcab\" と \"abcab\" → 5");
    SHOW(H.common_suffix(8, 5), "\"abcabcab\" と \"abcab\" の共通の末尾");
    SHOW(H.compare(1, 3, 0, 3), "\"bc\" と \"abc\" → 1（大きい）");
    RUN(const RollingHash T("abd"));
    SHOW(H.lcp(0, T, 0), "別の文字列 \"abd\" とも");
    SHOW(H.compare(0, 3, T, 0, 3), "\"abc\" < \"abd\"");
    note("O(1) の LCP・比較が要るなら 03_suffix_array の SuffixArray。");

    api("【3】 StrHash の +（連結）/ repeat(k) / hash_of(s) / hash_of_char(c)", "StrHash",
        "ハッシュどうしをつなげる（文字列を作らない）。途中に挿入・削除・入れ替えた文字列と比べるときに", "O(1) / O(log k)");
    SHOW(H.get(0, 2) + H.get(5, 8) == hash_of("abcab"), "\"ab\" + \"cab\"");
    SHOW(hash_of("ab").repeat(3) == hash_of("ababab"), "3 回つなげる");
    SHOW(hash_of("x").repeat(1000000000000LL).len, "長さ 10^12 の文字列のハッシュも O(log)");
    {
        // ABC284 F: T（長さ 2N）= S[0, i) + rev(S) + S[i, N) となる S と i
        const std::string Tt = "abcbac";
        const int N = 3;
        const RollingHash F(Tt), R(reversed(Tt));
        int ans = -1;
        for (int i = 0; i <= N && ans < 0; ++i)  // T[0, i) + T[i + N, 2N) と、T[i, i + N) を反転したもの
            if (F.prefix(i) + F.suffix(i + N) == R.get(N - i, 2 * N - i)) ans = i;
        note("ABC284 F（ABCBAC）: F.prefix(i) + F.suffix(i + N) と、反転した T のハッシュ R.get(N - i, 2N - i) を比べる");
        std::cout << "    T = abcbac → i = " << ans << "（S = " << Tt.substr(0, ans) + Tt.substr(ans + N) << "）\n";
    }

    api("【4】 unordered_map<StrHash, V> / H.find_all(pat) / H.count(pat)", "—",
        "ハッシュを map のキーに（同じ部分文字列をまとめる）。find_all は pat の出現位置（重なりも）", "O(1) / O(N)");
    {
        // ABC141 E: 重ならずに 2 回以上現れる部分文字列の最長（長さで二分探索）
        const std::string S2 = "ababa";
        const RollingHash G(S2);
        const int n = (int)S2.size();
        auto ok = [&](int len) {
            std::unordered_map<StrHash, int> first;  // ハッシュ → 最初の開始位置
            for (int i = 0; i + len <= n; ++i) {
                auto [it, ins] = first.emplace(G.get(i, i + len), i);
                if (!ins && i - it->second >= len) return true;
            }
            return false;
        };
        int lo = 0, hi = n / 2;
        while (lo < hi) {
            const int mid = (lo + hi + 1) / 2;
            if (ok(mid)) lo = mid;
            else hi = mid - 1;
        }
        note("ABC141 E（Who Says a Pun?）: 長さ len で二分探索し、unordered_map<StrHash, int> に最初の位置を覚える");
        std::cout << "    S = ababa → " << lo << "\n";
    }
    SHOW(H.find_all("ab"), "\"ab\" の出現位置");
    SHOW(H.count("abcab"), "重なりも数える");

    api("【5】 数の列にも使える", "—", "vector<int>・vector<long long> の部分列の一致（10^18 級・負の値も区別する）", "—");
    RUN(const std::vector<long long> A = {5, -1, 1000000000000000000LL, 5, -1});
    RUN(const RollingHash HA(A));
    SHOW(HA.same(0, 2, 3, 5), "{5, -1} と {5, -1}");

    api("【6】 RollingHash2D G(grid) / G.get(i1, j1, i2, j2) / G.find_all(pat)", "GridHash / vector<pair<int,int>>",
        "グリッドの長方形 [i1, i2) × [j1, j2) のハッシュ。大きなグリッドの中で模様を探す（小さいなら graph_grid_v2 の find_subgrid）", "構築 O(HW)・取得 O(1)");
    RUN(const std::vector<std::string> g = {"#..#", ".##.", "#..#"});
    RUN(const RollingHash2D G(g));
    SHOW(G.get(0, 0, 1, 2) == G.get(2, 2, 3, 4), "\"#.\" と \".#\" は違う");
    SHOW(G.get(0, 0, 1, 4) == G.get(2, 0, 3, 4), "1 行目と 3 行目は同じ");
    SHOW(G.find_all(std::vector<std::string>{"#.", ".#"}).size(), "2×2 の模様の数");
    SHOW(G.find_all(std::vector<std::string>{"#.", ".#"})[0].second, "最初の位置の列");

    api("【7】 center_palindrome_lengths_with_changes(s, k)", "vector<int>（長さ 2N-1）",
        "中心 c ごとに、食い違う組が k 個以下の最長の長さ（偶数 c は s[c/2] が中心の奇数長、奇数 c は間が中心の偶数長）。中心から LCP で k + 1 回広げる", "O(N k log N)");
    SHOW(center_palindrome_lengths_with_changes("abcda", 1), "1 文字まで書き換えてよい");
    {
        // ABC468 D 型：1 文字まで書き換えて回文にできる部分文字列の数 = 中心ごとに (L + 1) / 2 を足す（同じ偶奇の長さは全部できる）
        const auto L = center_palindrome_lengths_with_changes(std::string("abca"), 1);
        long long cnt = 0;
        for (int c = 0; c < (int)L.size(); ++c) cnt += c % 2 == 0 ? (L[c] + 1) / 2 : L[c] / 2;
        note("1 文字まで書き換えて回文になる部分文字列の数：中心ごとに長さ L 以下で同じ偶奇のもの（奇数中心は (L+1)/2 個、偶数中心は L/2 個）");
        std::cout << "    S = abca → " << cnt << "（a, b, c, a, ab, bc, ca, abc, bca, abca の 10 個）\n";
    }
}

}  // namespace rh_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    rh_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    rh_test::benchmark();    // 速度計測
    rh_demo::demo();         // 使い方デモ
    return 0;
}
