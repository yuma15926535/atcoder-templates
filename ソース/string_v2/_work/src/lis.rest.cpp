// =============================================================================
//  ランダムテスト
// =============================================================================
namespace lis_test {

{COMMON_TEST}

// テスト用の mod 998244353（数える関数に整数以外の型を渡す確認）
struct Mod {
    unsigned v = 0;
    Mod() = default;
    Mod(long long x) : v((unsigned)(((x % 998244353) + 998244353) % 998244353)) {}
    Mod operator+(const Mod& o) const { return Mod((long long)v + o.v); }
    Mod operator-(const Mod& o) const { return Mod((long long)v - o.v); }
    bool operator==(const Mod& o) const { return v == o.v; }
};

static bool less_by(long long x, long long y, bool strict) { return strict ? x < y : x <= y; }

// O(N^2) の DP：e[i] = i で終わる最長、s[i] = i から始まる最長、ce / cs = その数（mod）
struct Bf {
    std::vector<int> e, s;
    std::vector<Mod> ce, cs;
    int L = 0;
    Mod total;
};
static Bf bf_dp(const std::vector<long long>& v, bool strict) {
    const int n = (int)v.size();
    Bf b;
    b.e.assign(n, 1), b.s.assign(n, 1), b.ce.assign(n, Mod(1)), b.cs.assign(n, Mod(1));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j)
            if (less_by(v[j], v[i], strict)) {
                if (b.e[j] + 1 > b.e[i]) b.e[i] = b.e[j] + 1, b.ce[i] = b.ce[j];
                else if (b.e[j] + 1 == b.e[i]) b.ce[i] = b.ce[i] + b.ce[j];
            }
    for (int i = n - 1; i >= 0; --i)
        for (int j = i + 1; j < n; ++j)
            if (less_by(v[i], v[j], strict)) {
                if (b.s[j] + 1 > b.s[i]) b.s[i] = b.s[j] + 1, b.cs[i] = b.cs[j];
                else if (b.s[j] + 1 == b.s[i]) b.cs[i] = b.cs[i] + b.cs[j];
            }
    for (int i = 0; i < n; ++i) b.L = std::max(b.L, b.e[i]);
    b.total = n == 0 ? Mod(1) : Mod(0);
    for (int i = 0; i < n; ++i)
        if (b.e[i] == b.L) b.total = b.total + b.ce[i];
    return b;
}
// 増加部分列に分ける最小の数（バックトラック。n ≤ 9）
static int bf_cover(const std::vector<long long>& v, bool strict, bool increasing) {
    const int n = (int)v.size();
    for (int k = 0; k <= n; ++k) {
        std::vector<long long> last;
        bool ok = false;
        auto dfs = [&](auto&& self, int i) -> void {
            if (ok) return;
            if (i == n) {
                ok = true;
                return;
            }
            for (std::size_t c = 0; c < last.size(); ++c) {
                const bool fit = increasing ? less_by(last[c], v[i], strict) : less_by(v[i], last[c], strict);
                if (fit) {
                    const long long keep = last[c];
                    last[c] = v[i];
                    self(self, i + 1);
                    last[c] = keep;
                }
            }
            if ((int)last.size() < k) {
                last.push_back(v[i]);
                self(self, i + 1);
                last.pop_back();
            }
        };
        dfs(dfs, 0);
        if (ok) return k;
    }
    return n;
}

static std::vector<long long> rnd_values(std::mt19937_64& rng, int n) {
    std::vector<long long> v(n);
    const int t = (int)(rng() % 4);
    const long long range = 1 + (long long)(rng() % (t == 0 ? 3 : 2 * n + 2));
    for (int i = 0; i < n; ++i) {
        if (t == 3) v[i] = rnd_ll(rng);
        else v[i] = (long long)(rng() % range) - (t == 2 ? range / 2 : 0);
    }
    if (rng() % 8 == 0) std::sort(v.begin(), v.end());
    if (rng() % 8 == 0) std::sort(v.rbegin(), v.rend());
    return v;
}

static void check_restore(const std::vector<long long>& v, const std::vector<int>& r, int L, bool strict) {
    LIS_CHECK((int)r.size() == L);
    for (std::size_t k = 0; k < r.size(); ++k) {
        LIS_CHECK(r[k] >= 0 && r[k] < (int)v.size());
        if (k) LIS_CHECK(r[k - 1] < r[k] && less_by(v[r[k - 1]], v[r[k]], strict));
    }
}

static void test_small_exhaustive(std::mt19937_64& rng, int& cases) {
    // 全部の部分集合を調べる：LIS の長さ・数・入りうる位置・必ず入る位置
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 13);
        const auto v = rnd_values(rng, n);
        const bool strict = rng() % 2;
        int L = 0;
        long long count = 0;
        std::vector<int> through(n, 0);
        for (int mask = 0; mask < (1 << n); ++mask) {
            int prev = -1, len = 0;
            bool ok = true;
            for (int i = 0; i < n && ok; ++i)
                if (mask >> i & 1) {
                    if (prev >= 0 && !less_by(v[prev], v[i], strict)) ok = false;
                    prev = i, ++len;
                }
            if (!ok) continue;
            if (len > L) L = len, count = 0, std::fill(through.begin(), through.end(), 0);
            if (len == L) {
                ++count;
                for (int i = 0; i < n; ++i)
                    if (mask >> i & 1) ++through[i];
            }
        }
        LIS_CHECK(lis_length(v, strict) == L);
        check_restore(v, lis_restore(v, strict), L, strict);
        LIS_CHECK(count_lis<long long>(v, strict) == count);
        LIS_CHECK(count_lis<Mod>(v, strict) == Mod(count));
        const LisMembers mb = lis_members(v, strict);
        LIS_CHECK(mb.length == L && (int)mb.in_some.size() == n && (int)mb.in_all.size() == n);
        for (int i = 0; i < n; ++i) LIS_CHECK(mb.in_some[i] == (through[i] > 0) && mb.in_all[i] == (L > 0 && through[i] == count));
        if (n <= 9) {
            LIS_CHECK(min_increasing_cover(v, strict) == bf_cover(v, strict, true));
            LIS_CHECK(min_decreasing_cover(v, strict) == bf_cover(v, strict, false));
        }
        ++cases;
    }
}

static void test_dp(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % (it % 20 == 0 ? 600 : 60));
        const auto v = rnd_values(rng, n);
        for (int st = 0; st < 2; ++st) {
            const bool strict = st == 1;
            const Bf b = bf_dp(v, strict);
            LIS_CHECK(lis_length(v, strict) == b.L);
            LIS_CHECK(lis_ending_at(v, strict) == b.e && lis_starting_at(v, strict) == b.s);
            check_restore(v, lis_restore(v, strict), b.L, strict);
            LIS_CHECK(count_lis<Mod>(v, strict) == b.total);
            // 減少：大小を逆にした LIS
            std::vector<long long> w(v);
            for (auto& x : w) x = ~x;
            LIS_CHECK(lds_length(v, strict) == bf_dp(w, strict).L);
            LIS_CHECK(min_increasing_cover(v, strict) == lds_length(v, !strict) && min_decreasing_cover(v, strict) == lis_length(v, !strict));
            const LisMembers mb = lis_members(v, strict);
            for (int i = 0; i < n; ++i) {
                const bool some = b.e[i] + b.s[i] - 1 == b.L;
                // すべての LIS に入る ⇔ i を通る LIS の数（ce[i] × cs[i]）が全体の数と同じ（mod で比べる）
                const Mod through((long long)((unsigned long long)b.ce[i].v * b.cs[i].v % 998244353ULL));
                LIS_CHECK(mb.in_some[i] == some);
                LIS_CHECK(mb.in_all[i] == (some && through == b.total));
            }
        }
        // 文字列でも同じ
        const std::string s = rnd_str_shaped(rng, (int)(rng() % 40), 1 + (int)(rng() % 26));
        const std::vector<long long> sv(s.begin(), s.end());
        LIS_CHECK(lis_length(s) == lis_length(sv) && lis_length(s, false) == lis_length(sv, false) && lis_restore(s) == lis_restore(sv));
        ++cases;
    }
}

static void test_2d(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)(rng() % 60);
        const long long range = 1 + (long long)(rng() % 12);
        std::vector<std::pair<long long, int>> ps(n);
        for (auto& p : ps) p = {(long long)(rng() % range) - 3, (int)(rng() % range)};
        // 愚直：x で並べて O(N^2)
        std::vector<int> ord(n);
        for (int i = 0; i < n; ++i) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](int i, int j) { return ps[i].first < ps[j].first; });
        std::vector<int> dp(n, 1);
        int L = 0;
        for (int a = 0; a < n; ++a) {
            for (int b = 0; b < a; ++b)
                if (ps[ord[b]].first < ps[ord[a]].first && ps[ord[b]].second < ps[ord[a]].second) dp[a] = std::max(dp[a], dp[b] + 1);
            L = std::max(L, dp[a]);
        }
        const auto r = lis_2d(ps);
        LIS_CHECK((int)r.size() == L);
        for (std::size_t k = 1; k < r.size(); ++k) LIS_CHECK(ps[r[k - 1]].first < ps[r[k]].first && ps[r[k - 1]].second < ps[r[k]].second);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261012);
    int cases = 0;
    test_small_exhaustive(rng, cases);
    test_dp(rng, cases);
    test_2d(rng, cases);
    // 途中の数が long long を超えても、答えが収まれば止まらない
    {
        std::vector<long long> v;
        for (int k = 0; k <= 62; ++k) v.push_back(1000 + 2 * k + 1), v.push_back(1000 + 2 * k);
        for (int k = 0; k < 64; ++k) v.push_back(k);
        LIS_CHECK(lis_length(v) == 64 && count_lis<long long>(v) == 1);
        ++cases;
    }
    // lis_2d は < だけで比べる（double・unsigned long long）
    {
        LIS_CHECK(lis_2d(std::vector<std::pair<double, double>>{{0, 0.1}, {1, 0.5}}).size() == 2);
        LIS_CHECK(lis_2d(std::vector<std::pair<long long, unsigned long long>>{{0, 1}, {1, (1ULL << 63) + 5}}).size() == 2);
        ++cases;
    }
    std::cout << "[lis] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[lis] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 1000000;
    std::vector<int> a(N);
    for (auto& x : a) x = (int)(rng() % 1000000000);
    auto t0 = std::chrono::steady_clock::now();
    long long x = lis_length(a);
    const long long t_len = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += (long long)lis_restore(a).size();
    const long long t_res = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const LisMembers mb = lis_members(a);
    x += mb.in_some[N / 2];
    const long long t_mem = ms(t0);
    t0 = std::chrono::steady_clock::now();
    x += count_lis<Mod>(a).v;
    const long long t_cnt = ms(t0);
    std::vector<std::pair<int, int>> ps(N);
    for (auto& p : ps) p = {(int)(rng() % 1000000), (int)(rng() % 1000000)};
    t0 = std::chrono::steady_clock::now();
    x += (long long)lis_2d(ps).size();
    const long long t_2d = ms(t0);
    std::cout << "[lis] benchmark: N = 10^6  lis_length " << t_len << " ms, lis_restore " << t_res << " ms, lis_members " << t_mem << " ms, count_lis "
              << t_cnt << " ms, lis_2d " << t_2d << " ms  (check " << x % 7 << ")\n";
}

}  // namespace lis_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace lis_demo {

{COMMON_DEMO}

static void demo() {
    banner("最長増加部分列（lis） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  lis_length / lis_restore / lis_ending_at / lis_starting_at   長さ・1 つ（添字）・各位置で終わる / 始まる長さ\n"
              << "  lis_members / count_lis<T>                                 入りうる・必ず入る位置 / LIS の数\n"
              << "  lds_length / min_increasing_cover / min_decreasing_cover    減少・増加列に分ける最小の数（Dilworth）\n"
              << "  lis_2d                                                     (x, y) が両方とも増える最長の列（入れ子の箱）\n"
              << "  strict = true（既定）で狭義（<）、false で広義（≤）。列は vector<整数> でも string でもよい\n";

    const std::vector<int> a = {3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5};
    code("const vector<int> a = {3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5};");

    api("【1】 lis_length(a, strict) / lis_restore(a, strict)", "int / vector<int>（添字。昇順）", "最長増加部分列の長さと、その 1 つ", "O(N log N)");
    SHOW(lis_length(a), "3, 4, 5, 6 など（狭義）");
    SHOW(lis_length(a, false), "1, 1, 2, 3, 5 など（広義 = 減らない）");
    SHOW(lis_restore(a), "a[3], a[6], a[9], a[10] = 1, 2, 3, 5（いくつかあるうちの 1 つ）");
    SHOW(lis_length(std::string("abracadabra")), "文字列でも（a < b < c < d < r）");

    api("【2】 lis_ending_at(a, strict) / lis_starting_at(a, strict)", "vector<int>", "各位置で終わる / 始まる最長の増加部分列の長さ（両方を足して 1 を引くと、その位置を通る最長）", "O(N log N)");
    SHOW(lis_ending_at(a), "");
    SHOW(lis_starting_at(a), "");

    api("【3】 lis_members(a, strict)", "LisMembers{length, in_some, in_all}", "各位置が「どれかの LIS に入る」「すべての LIS に入る」か", "O(N log N)");
    {
        const LisMembers mb = lis_members(a);
        SHOW(mb.in_some, "ABC354 F（Useless for LIS）はこれが true の位置");
        SHOW(mb.in_all, "");
        const std::vector<int> c = {2, 5, 3, 4, 3, 4};
        const LisMembers m2 = lis_members(c);
        std::vector<int> idx;
        for (int i = 0; i < (int)c.size(); ++i)
            if (m2.in_some[i]) idx.push_back(i + 1);
        note("ABC354 F の入力例 2 の 1 つ目（2 5 3 4 3 4）→ " + to_s((int)idx.size()) + " 個: " + join(idx));
    }

    api("【4】 count_lis<T>(a, strict)", "T", "LIS の数（添字の列として。同じ値の並びでも位置が違えば別）。T は ModInt998・long long（溢れたら止める）など", "O(N log N)");
    SHOW(count_lis<long long>(a), "");
    SHOW(count_lis<long long>(std::vector<int>{2, 1, 4, 3, 6, 5}), "2 × 2 × 2");

    api("【5】 lds_length(a, strict) / min_increasing_cover(a, strict) / min_decreasing_cover(a, strict)", "int",
        "最長減少部分列 / 増加（減少）部分列に分ける最小の数（Dilworth：狭義増加に分ける数 = 最長の「減らない」部分列）", "O(N log N)");
    SHOW(lds_length(a), "9, 6, 5, 3 など");
    SHOW(min_increasing_cover(std::vector<int>{2, 1, 4, 5, 3}), "ABC134 E（Sequence Decomposing）の入力例 1");
    SHOW(min_increasing_cover(std::vector<int>{0, 0, 0, 0}), "入力例 2：同じ値は同じ色にできない");

    api("【6】 lis_2d(ps)", "vector<int>（添字。x の小さい順）", "(x, y) が両方とも狭義に増える最長の列（x の昇順・同じ x は y の降順に並べて y の LIS）", "O(N log N)");
    {
        const std::vector<std::pair<int, int>> box = {{8, 8}, {5, 3}, {2, 2}, {4, 2}, {2, 1}};
        code("const vector<pair<int, int>> box = {{8, 8}, {5, 3}, {2, 2}, {4, 2}, {2, 1}};   // ABC038 D（プレゼント）の入力例 4");
        SHOW(lis_2d(box).size(), "入れ子にできる箱の最大の数");
        SHOW(lis_2d(box), "(2, 1) → (4, 2) → (5, 3) → (8, 8)");
    }

    head("問題での使い方：ABC360 G（Suitable Edit for LIS）");
    note("1 か所を好きな整数に変えて LIS を最大に → e = lis_ending_at(A)、s = lis_starting_at(A) を作り、");
    note("j < k、k - j ≥ 2、A[k] - A[j] ≥ 2 の組で e[j] + 1 + s[k] の最大（間の 1 つを A[j] + 1 に変える）。両端の外も番兵で。");
    note("本番は j を値の BIT（最大値）に入れて O(N log N)。ここでは入力例を O(N^2) で:");
    {
        for (const auto& A : {std::vector<long long>{3, 2, 2, 4}, std::vector<long long>{4, 5, 3, 6, 7}}) {
            const int n = (int)A.size();
            const auto e = lis_ending_at(A), s = lis_starting_at(A);
            int best = 0;
            for (int x = 0; x < n; ++x) {  // x を変える。左は j < x、右は k > x
                int left = 0, right = 0;
                for (int j = 0; j < x; ++j) left = std::max(left, e[j]);
                for (int k = x + 1; k < n; ++k) right = std::max(right, s[k]);
                best = std::max({best, left + 1, 1 + right});
                for (int j = 0; j < x; ++j)
                    for (int k = x + 1; k < n; ++k)
                        if (A[k] - A[j] >= 2) best = std::max(best, e[j] + 1 + s[k]);
            }
            note("A = " + to_s(A) + " → " + to_s(best));
        }
    }
}

}  // namespace lis_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    lis_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    lis_test::benchmark();    // 速度計測
    lis_demo::demo();         // 使い方デモ
    return 0;
}
