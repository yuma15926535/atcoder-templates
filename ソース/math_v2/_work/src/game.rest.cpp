// =============================================================================
//  ランダムテスト
// =============================================================================
namespace game_test {

{COMMON_TEST}

// 山の組を状態にして、全部の手を試す（勝ち負けだけ）
static bool brute_win(std::vector<int> piles, const std::vector<int>& S, bool misere, std::map<std::vector<int>, bool>& memo) {
    std::sort(piles.begin(), piles.end());
    auto it = memo.find(piles);
    if (it != memo.end()) return it->second;
    bool any_move = false, win = false;
    for (size_t i = 0; i < piles.size() && !win; ++i)
        for (int s : S) {
            if (s > piles[i]) continue;
            any_move = true;
            std::vector<int> q = piles;
            q[i] -= s;
            if (!brute_win(q, S, misere, memo)) { win = true; break; }
        }
    if (!any_move) win = misere;  // 動けない：普通は負け、逆形は（最後に取ったのが相手なので）勝ち
    return memo[piles] = win;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 10);
        std::vector<long long> v((size_t)n);
        for (auto& x : v) x = rnd(rng, -2, 12);
        int e = 0;
        while (std::find(v.begin(), v.end(), (long long)e) != v.end()) ++e;
        GAME_CHECK(mex(v) == e);
        std::vector<unsigned> u;
        for (long long x : v)
            if (x >= 0) u.push_back((unsigned)x);
        GAME_CHECK(mex(u) == e);
    }
    // Nim・逆形 Nim：全部の手を試す
    std::map<std::vector<int>, bool> m1, m2;
    const std::vector<int> all{1, 2, 3, 4, 5, 6, 7};
    for (int it = 0; it < 3000; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 4);
        std::vector<int> p((size_t)n);
        for (auto& x : p) x = (int)rnd(rng, 0, 7);
        GAME_CHECK(nim_first_wins(p) == brute_win(p, all, false, m1));
        GAME_CHECK(misere_nim_first_wins(p) == brute_win(p, all, true, m2));
    }
}

static void test_dag(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 40);
        std::vector<int> perm((size_t)n);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        std::vector<std::vector<int>> mv((size_t)n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j)
                if (rng() % 4 == 0) mv[(size_t)perm[(size_t)i]].push_back(perm[(size_t)j]);  // 番号の大きい方から小さい方へ（並べ替えた番号で）
        for (auto& v : mv)
            if (!v.empty() && rng() % 3 == 0) v.push_back(v[0]);  // 同じ行き先が 2 回
        const std::vector<int> g = grundy_dag(n, mv);
        std::vector<int> e((size_t)n, -1);  // 小さい方から順に
        for (int i = 0; i < n; ++i) {
            const int v = perm[(size_t)i];
            std::vector<int> ws;
            for (int w : mv[(size_t)v]) ws.push_back(e[(size_t)w]);
            e[(size_t)v] = mex(ws);
        }
        GAME_CHECK(g == e);
    }
    // 長い 1 本道（再帰しないので深くても大丈夫）
    const int N = 1000000;
    std::vector<std::vector<int>> mv((size_t)N);
    for (int i = 1; i < N; ++i) mv[(size_t)i].push_back(i - 1);
    const std::vector<int> g = grundy_dag(N, mv);
    GAME_CHECK(g[0] == 0 && g[1] == 1 && g[N - 1] == (N - 1) % 2);
    ++cases;
}

static void test_multi_and_outcome(std::mt19937_64& rng, int& cases) {
    // 1 手で複数のゲームに分かれる
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 30);
        std::vector<std::vector<std::vector<int>>> mv((size_t)n);
        std::vector<std::vector<int>> single((size_t)n);
        bool all_single = it % 3 == 0;
        for (int v = 1; v < n; ++v) {  // 小さい番号へだけ
            const int k = (int)rnd(rng, 0, 4);
            for (int t = 0; t < k; ++t) {
                std::vector<int> to;
                const int len = all_single ? 1 : (int)rnd(rng, 0, 3);
                for (int j = 0; j < len; ++j) to.push_back((int)rnd(rng, 0, v - 1));
                mv[(size_t)v].push_back(to);
                if (all_single) single[(size_t)v].push_back(to[0]);
            }
        }
        std::vector<int> e((size_t)n);
        for (int v = 0; v < n; ++v) {
            std::vector<int> vals;
            for (const auto& to : mv[(size_t)v]) {
                int x = 0;
                for (int w : to) x ^= e[(size_t)w];
                vals.push_back(x);
            }
            e[(size_t)v] = mex(vals);
        }
        GAME_CHECK(grundy_dag_multi(n, mv) == e && (!all_single || grundy_dag(n, single) == e));
    }
    // 後退解析：手数ごとに決める愚直解と（閉路あり）
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 25);
        std::vector<std::vector<int>> mv((size_t)n);
        for (int v = 0; v < n; ++v)
            for (int w = 0; w < n; ++w)
                if (rng() % (it % 2 ? 4 : 8) == 0) mv[(size_t)v].push_back(w);
        std::vector<int> res((size_t)n, GAME_DRAW), st((size_t)n, -1);
        for (int v = 0; v < n; ++v)
            if (mv[(size_t)v].empty()) res[(size_t)v] = GAME_LOSE, st[(size_t)v] = 0;
        for (int k = 1; k <= 2 * n + 2; ++k) {  // 手数 k で決まるもの
            std::vector<int> nr = res, ns = st;
            for (int u = 0; u < n; ++u) {
                if (res[(size_t)u] != GAME_DRAW) continue;
                bool win = false, all_win = true;
                for (int w : mv[(size_t)u]) {
                    win |= res[(size_t)w] == GAME_LOSE && st[(size_t)w] == k - 1;
                    all_win &= res[(size_t)w] == GAME_WIN;
                }
                if (win) nr[(size_t)u] = GAME_WIN, ns[(size_t)u] = k;
                else if (all_win) nr[(size_t)u] = GAME_LOSE, ns[(size_t)u] = k;
            }
            res = nr, st = ns;
        }
        const GameOutcome o = game_outcome(n, mv);
        GAME_CHECK(o.result == res && o.steps == st);
    }
}

static void test_take(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 400 * TEST_SCALE; ++it, ++cases) {
        const long long L = rnd(rng, 1, 20), R = rnd(rng, L, 40);
        std::vector<long long> S;
        for (long long s = L; s <= R; ++s) S.push_back(s);
        const std::vector<int> g = subtraction_game_grundy(S, 300);
        for (int x = 0; x <= 300; ++x) GAME_CHECK(range_take_grundy(x, L, R) == g[(size_t)x]);
        const long long X = rnd(rng, 0, 1000000000000000000LL);
        const SubtractionGame G(S);
        GAME_CHECK(G.grundy(X) == range_take_grundy(X, L, R) && G.period() == L + R);
    }
    GAME_CHECK(range_take_grundy(LLONG_MAX, LLONG_MAX, LLONG_MAX) == 1 && range_take_grundy(5, 1, 1) == 1);
    // 集合 S：メモ化探索の Grundy 数と、周期
    for (int it = 0; it < 400 * TEST_SCALE; ++it, ++cases) {
        const int k = (int)rnd(rng, 1, 5);
        std::vector<long long> S;
        for (int i = 0; i < k; ++i) S.push_back(rnd(rng, 1, it % 2 ? 10 : 40));
        const int N = 3000;
        const std::vector<int> g = subtraction_game_grundy(S, N);
        for (int x = 0; x <= N; ++x) {  // 定義どおり
            std::vector<int> ws;
            for (long long s : S)
                if (s <= x) ws.push_back(g[(size_t)(x - s)]);
            GAME_CHECK(g[(size_t)x] == mex(ws));
        }
        const SubtractionGame G(S);
        for (int x = 0; x <= N; ++x) GAME_CHECK(G.grundy(x) == g[(size_t)x]);
        const long long P = G.preperiod(), T = G.period();
        GAME_CHECK(T >= 1 && P >= 0);
        if (P + 2 * T <= N) {
            for (long long x = P; x + T <= N; ++x) GAME_CHECK(g[(size_t)x] == g[(size_t)(x + T)]);
            GAME_CHECK(P == 0 || g[(size_t)P - 1] != g[(size_t)(P - 1 + T)]);  // 周期はできるだけ前から
        }
        const long long X = rnd(rng, P + T, 1000000000000000000LL);
        GAME_CHECK(G.grundy(X) == G.grundy(X - T));
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_basic(rng, cases);
    test_dag(rng, cases);
    test_multi_and_outcome(rng, cases);
    test_take(rng, cases);
    std::cout << "[game] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[game] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 200000;  // 約数へ動けるゲーム（ABC368 F の形）
    std::vector<std::vector<int>> mv((size_t)N + 1);
    for (int d = 1; d <= N; ++d)
        for (int m = 2 * d; m <= N; m += d) mv[(size_t)m].push_back(d);
    auto t0 = std::chrono::steady_clock::now();
    const std::vector<int> g = grundy_dag(N + 1, mv);
    const long long t_dag = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<int> h = subtraction_game_grundy({1, 3, 7, 19, 31, 47, 53, 61, 97, 101}, 3000000);
    const long long t_sub = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const SubtractionGame G({2, 5, 11, 23, 37, 89, 123});
    const long long t_per = ms(t0);
    std::cout << "[game] benchmark: grundy_dag（2×10^5 状態, 約 2.2×10^6 辺）" << t_dag << " ms, subtraction_game_grundy（3×10^6, |S| = 10）" << t_sub
              << " ms, SubtractionGame の周期 " << t_per << " ms（周期 " << G.period() << "）  (check " << (g[N] + h.back() + G.grundy(1000000000000LL)) % 7 << ")\n";
}

}  // namespace game_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace game_demo {

{COMMON_DEMO}

static void demo() {
    banner("ゲーム（game） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  mex / grundy_dag / nim_first_wins / misere_nim_first_wins / range_take_grundy / subtraction_game_grundy / SubtractionGame\n"
              << "  先手必勝 ⇔ 並べたゲームの Grundy 数の XOR ≠ 0\n";

    api("【1】 mex(vs)", "int", "vs に含まれない最小の 0 以上の整数（Grundy 数 = 行き先の Grundy 数の mex）", "O(|vs|)");
    SHOW(mex(std::vector<int>{0, 1, 3}), "");
    SHOW(mex(std::vector<int>{1, 2}), "");
    SHOW(mex(std::vector<int>{}), "行き先がない（負けの状態）→ 0");
    {
        // ABC206 F（Interval Game 2）：g[l][r] = mex{ g[l][L_i] ^ g[R_i][r] : [L_i, R_i) ⊆ [l, r) }
        auto solve = [](const std::vector<std::pair<int, int>>& iv) {
            std::vector<std::vector<int>> g(101, std::vector<int>(101, 0));
            for (int len = 1; len <= 99; ++len)
                for (int l = 1; l + len <= 100; ++l) {
                    const int r = l + len;
                    std::vector<int> ws;
                    for (const auto& p : iv)
                        if (l <= p.first && p.second <= r) ws.push_back(g[(size_t)l][(size_t)p.first] ^ g[(size_t)p.second][(size_t)r]);
                    g[(size_t)l][(size_t)r] = mex(ws);
                }
            return std::string(g[1][100] ? "Alice" : "Bob");
        };
        note("ABC206 F（Interval Game 2）: 区間 [l, r) の Grundy 数 = mex（選んだ区間の左右の Grundy 数の XOR）");
        std::cout << "    例 → " << solve({{53, 98}, {8, 43}, {12, 53}}) << ", "
                  << solve({{4, 7}, {5, 7}, {3, 7}, {4, 5}, {5, 8}, {6, 9}, {4, 8}, {5, 10}, {1, 9}, {5, 10}}) << ", " << solve({{58, 98}, {11, 29}}) << ", "
                  << solve({{79, 83}, {44, 83}, {38, 74}, {49, 88}, {18, 45}, {64, 99}}) << ", " << solve({{5, 9}}) << "\n";
    }

    api("【2】 grundy_dag(n, moves)", "vector<int>", "状態 0〜n−1、moves[v] = v から 1 手で行ける状態。各状態の Grundy 数（閉路があれば止める。再帰しない）", "O(n + 辺の数)");
    RUN(std::vector<std::vector<int>> mv = {{}, {0}, {0, 1}, {1, 2}});
    SHOW(grundy_dag(4, mv), "0 → 行き先なし、1 → {0}、2 → {0, 1}、3 → {1, 2}");
    {
        // ABC368 F（Dividing Game）：A_i を真の約数にする。状態 = 値、行き先 = 真の約数
        auto solve = [](const std::vector<int>& A) {
            const int M = *std::max_element(A.begin(), A.end());
            std::vector<std::vector<int>> moves((size_t)M + 1);
            for (int d = 1; d <= M; ++d)
                for (int m = 2 * d; m <= M; m += d) moves[(size_t)m].push_back(d);
            const std::vector<int> g = grundy_dag(M + 1, moves);
            int x = 0;
            for (int a : A) x ^= g[(size_t)a];
            return std::string(x ? "Anna" : "Bruno");
        };
        note("ABC368 F（Dividing Game）: 状態 = 値、行き先 = 真の約数 で grundy_dag、XOR");
        std::cout << "    例 1 → " << solve({2, 3, 4}) << "、例 2 → " << solve({2, 3, 4, 6}) << "\n";
    }

    api("【3】 nim_first_wins(piles) / misere_nim_first_wins(piles)", "bool", "Nim（最後に取った方が勝ち）/ 逆形 Nim（最後に取った方が負け）で先手必勝か", "O(N)");
    SHOW(nim_first_wins(std::vector<int>{1, 2, 3}), "1 ^ 2 ^ 3 = 0 → 後手必勝");
    SHOW(nim_first_wins(std::vector<int>{3, 4, 5}), "");
    SHOW(misere_nim_first_wins(std::vector<int>{1, 1}), "全部 1：1 の山が偶数個なら先手必勝");
    SHOW(misere_nim_first_wins(std::vector<int>{1, 2, 3}), "2 以上がある：普通の Nim と同じ");

    api("【4】 range_take_grundy(x, L, R)", "long long", "1 回に L〜R 個取れる石取りの x 個の山の Grundy 数 ⌊(x mod (L + R)) / L⌋", "O(1)");
    SHOW(range_take_grundy(10, 1, 2), "10 mod 3 = 1");
    SHOW(range_take_grundy(7, 3, 5), "7 mod 8 = 7 → ⌊7/3⌋");
    {
        // ABC297 G（Constrained Nim 2）
        auto solve = [](long long L, long long R, const std::vector<long long>& A) {
            long long x = 0;
            for (long long a : A) x ^= range_take_grundy(a, L, R);
            return std::string(x ? "First" : "Second");
        };
        note("ABC297 G（Constrained Nim 2）: 各山の range_take_grundy の XOR");
        std::cout << "    例 1 → " << solve(1, 2, {2, 3, 3}) << "、例 2 → " << solve(1, 1, {3, 1, 4, 1, 5}) << "、例 3 → " << solve(3, 14, {10, 20, 30, 40, 50, 60, 70}) << "\n";
    }

    api("【5】 subtraction_game_grundy(S, N) / SubtractionGame G(S)", "vector<int> / 構造体",
        "1 回に取れる数が集合 S の石取りの Grundy 数 g(0)〜g(N) / x が巨大なとき：周期を見つけて G.grundy(x)", "O(N |S|) / O(周期が見つかるまで × |S|)");
    SHOW(subtraction_game_grundy({1, 3, 4}, 14), "取れる数 {1, 3, 4}");
    RUN(SubtractionGame G({1, 3, 4}));
    SHOW(G.period(), "周期 7（0 1 0 1 2 3 2 のくり返し）");
    SHOW(G.preperiod(), "最初から周期的");
    SHOW(G.grundy(1000000000000000000LL), "10^18 mod 7 = 1 → g(1)");
    SHOW(subtraction_game_grundy({2, 4, 7}, 16), "取れる数 {2, 4, 7}");
    RUN(SubtractionGame H({2, 4, 7}));
    SHOW(H.preperiod(), "g(8) から周期的になる（途中から周期的な例）");
    SHOW(H.period(), "1 0 2 のくり返し");

    api("【6】 grundy_dag_multi(n, moves)", "vector<int>",
        "1 手で「いくつかのゲームに分かれる」とき：moves[v] は手の一覧、1 つの手は行き先の状態の列（その手の値 = 行き先の Grundy 数の XOR）", "O(n + 行き先の総数)");
    {
        // グランディのゲーム：山を大きさの違う 2 つの山に分ける
        const int N = 16;
        std::vector<std::vector<std::vector<int>>> moves2(N);
        for (int x = 0; x < N; ++x)
            for (int a = 1; 2 * a < x; ++a) moves2[(size_t)x].push_back({a, x - a});
        const std::vector<int> g = grundy_dag_multi(N, moves2);
        std::cout << "  グランディのゲーム（山を大きさの違う 2 つに分ける）g(0)〜g(15) -> ";
        for (int x : g) std::cout << x << " ";
        std::cout << "\n";
    }

    api("【7】 game_outcome(n, moves)", "GameOutcome{result, steps}",
        "閉路があってもよいゲームの後退解析。result[v] = GAME_WIN / GAME_LOSE / GAME_DRAW（v で手番の人から見て）、steps[v] = 決着までの手数（勝つ方は最短、負ける方は最長）",
        "O(n + 辺の数)");
    {
        // ABC209 E（Shiritori）：単語 = 状態、末尾 3 文字 = 次の単語の先頭 3 文字なら行ける。高橋くんが単語 i を言うと青木くんの手番
        const std::vector<std::string> w = {"abcd", "bcda", "ada"};
        const int n = (int)w.size();
        std::vector<std::vector<int>> moves2(n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (w[(size_t)i].substr(w[(size_t)i].size() - 3) == w[(size_t)j].substr(0, 3)) moves2[(size_t)i].push_back(j);
        const GameOutcome o = game_outcome(n, moves2);
        note("ABC209 E（Shiritori）: 単語ごとの状態で後退解析（実際は「先頭 3 文字」の状態にまとめて辺を減らす）");
        std::cout << "    ";
        for (int i = 0; i < n; ++i) std::cout << (o.result[(size_t)i] == GAME_WIN ? "Aoki" : o.result[(size_t)i] == GAME_LOSE ? "Takahashi" : "Draw") << " ";
        std::cout << "（steps: " << o.steps[0] << " " << o.steps[1] << " " << o.steps[2] << "）\n";
    }
}

}  // namespace game_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    game_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    game_test::benchmark();    // 速度計測
    game_demo::demo();         // 使い方デモ
    return 0;
}
