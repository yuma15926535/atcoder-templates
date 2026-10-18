#ifndef GG_GAME
#define GG_GAME
// ---- ゲーム（数学第二世代）  mex・Grundy 数（DAG・石取り・分かれる手）・Nim・閉路のあるゲームの後退解析 ----
// 「動けなくなった方が負け」の 2 人ゲームで、いくつかのゲームを並べたもの（どれか 1 つで 1 手）の先手必勝 ⇔ 各 Grundy 数の XOR ≠ 0
template <class T>
int mex(const std::vector<T>& vs) {  // vs に含まれない最小の 0 以上の整数（負の値は無視）  O(|vs|)
    std::vector<char> seen(vs.size() + 1, 0);
    for (const T& v : vs) {
        if constexpr (std::is_signed<T>::value)
            if (v < 0) continue;
        if ((unsigned long long)v <= vs.size()) seen[(size_t)v] = 1;
    }
    int m = 0;
    while (seen[(size_t)m]) ++m;
    return m;
}
// 状態 0〜n−1、moves[v] = v から 1 手で行ける状態。各状態の Grundy 数（閉路があれば止める）  O(n + 辺の数)
inline std::vector<int> grundy_dag(int n, const std::vector<std::vector<int>>& moves) {
    if (n < 0 || (int)moves.size() != n) math_internal::fail("grundy_dag", "moves の大きさは n");
    std::vector<int> g((size_t)n, -1), st;
    std::vector<char> state((size_t)n, 0);  // 0 = まだ、1 = 調べている途中、2 = 決まった
    std::vector<size_t> it((size_t)n, 0);
    std::vector<char> seen;
    for (int s = 0; s < n; ++s) {
        if (state[(size_t)s]) continue;
        st.push_back(s), state[(size_t)s] = 1;
        while (!st.empty()) {  // 行き先を全部決めてから自分を決める（再帰を使わない深さ優先）
            const int v = st.back();
            const std::vector<int>& mv = moves[(size_t)v];
            if (it[(size_t)v] < mv.size()) {
                const int w = mv[it[(size_t)v]++];
                if (w < 0 || w >= n) math_internal::fail("grundy_dag", "行き先の番号は 0〜n−1");
                if (state[(size_t)w] == 1) math_internal::fail("grundy_dag", "閉路がある（同じ状態に戻れる）");
                if (state[(size_t)w] == 0) st.push_back(w), state[(size_t)w] = 1;
                continue;
            }
            seen.assign(mv.size() + 1, 0);
            for (int w : mv)
                if ((size_t)g[(size_t)w] <= mv.size()) seen[(size_t)g[(size_t)w]] = 1;
            int m = 0;
            while (seen[(size_t)m]) ++m;
            g[(size_t)v] = m, state[(size_t)v] = 2, st.pop_back();
        }
    }
    return g;
}
// Nim：好きな山から 1 個以上取る。最後に取った方が勝ち ⇔ XOR ≠ 0 なら先手必勝
template <class T>
bool nim_first_wins(const std::vector<T>& piles) {
    unsigned long long x = 0;
    for (const T& a : piles) x ^= (unsigned long long)a;
    return x != 0;
}
// 逆形 Nim（最後に取った方が負け）：全部の山が 1 以下なら「1 の山が偶数個」、そうでなければ XOR ≠ 0 なら先手必勝
template <class T>
bool misere_nim_first_wins(const std::vector<T>& piles) {
    unsigned long long x = 0;
    bool big = false;
    long long ones = 0;
    for (const T& a : piles) x ^= (unsigned long long)a, big = big || a > 1, ones += a == 1;
    return big ? x != 0 : ones % 2 == 0;
}
// 1 回に L〜R 個取れる石取り（L ≥ 1）の、x 個の山の Grundy 数 ⌊(x mod (L + R)) / L⌋  O(1)
inline long long range_take_grundy(long long x, long long L, long long R) {
    if (x < 0 || L < 1 || R < L) math_internal::fail("range_take_grundy", "x ≥ 0、1 ≤ L ≤ R で使う");
    return (long long)(((UInt128)x % ((UInt128)L + (UInt128)R)) / (UInt128)L);
}
// 1 回に取れる数が集合 S（S の要素は 1 以上）の石取りの Grundy 数 g(0), …, g(N)  O(N |S|)
inline std::vector<int> subtraction_game_grundy(const std::vector<long long>& S, int N) {
    for (long long s : S)
        if (s < 1) math_internal::fail("subtraction_game_grundy", "取れる数は 1 以上");
    if (N < 0) return {};
    std::vector<int> g((size_t)N + 1, 0);
    std::vector<int> mark(S.size() + 2, -1);
    for (int x = 0; x <= N; ++x) {
        for (long long s : S)
            if (s <= x && (size_t)g[(size_t)(x - s)] <= S.size()) mark[(size_t)g[(size_t)(x - s)]] = x;
        int m = 0;
        while (mark[(size_t)m] == x) ++m;
        g[(size_t)x] = m;
    }
    return g;
}
// 同じ石取りで x が巨大（10^18）なとき：Grundy 数の列はいつか周期的になるので、周期を見つけて使う
//   SubtractionGame G(S); G.grundy(x) / G.preperiod()（周期が始まる位置）/ G.period()
//   周期は limit までに見つからなければ止める。O(limit × |S|)
struct SubtractionGame {
    explicit SubtractionGame(const std::vector<long long>& S, int limit = 2000000) {
        if (S.empty()) {
            pre_ = 0, per_ = 1, g_ = {0};
            return;
        }
        long long smax = 0;
        for (long long s : S) smax = std::max(smax, s);
        if (smax > limit) math_internal::fail("SubtractionGame", "取れる数が limit より大きい");
        const int w = (int)smax;  // 直前 w 個の Grundy 数が同じなら、その先は全部同じ
        const unsigned long long B = 0x9E3779B97F4A7C15ULL;
        unsigned long long Bw = 1;
        for (int i = 0; i < w; ++i) Bw *= B;
        std::unordered_map<unsigned long long, int> first;  // 窓のハッシュ → 位置
        std::vector<int> mark(S.size() + 2, -1);
        unsigned long long h = 0;
        for (int x = 0; x <= limit; ++x) {
            if (x >= w) {  // 窓 g[x−w, x) を調べる
                auto itr = first.find(h);
                if (itr != first.end() && std::equal(g_.begin() + (itr->second - w), g_.begin() + itr->second, g_.begin() + (x - w))) {
                    per_ = x - itr->second, pre_ = itr->second - w;
                    while (pre_ > 0 && g_[(size_t)pre_ - 1] == g_[(size_t)(pre_ - 1 + per_)]) --pre_;  // 周期をもっと前から
                    return;
                }
                first.emplace(h, x);
            }
            for (long long s : S)
                if (s <= x && (size_t)g_[(size_t)(x - s)] <= S.size()) mark[(size_t)g_[(size_t)(x - s)]] = x;
            int m = 0;
            while (mark[(size_t)m] == x) ++m;
            g_.push_back(m);
            h = h * B + (unsigned long long)(m + 1);
            if (x >= w) h -= Bw * (unsigned long long)(g_[(size_t)(x - w)] + 1);
        }
        math_internal::fail("SubtractionGame", "limit までに周期が見つからない");
    }
    int grundy(long long x) const {
        if (x < 0) math_internal::fail("SubtractionGame::grundy", "x ≥ 0 で使う");
        if (x < (long long)g_.size()) return g_[(size_t)x];
        return g_[(size_t)(pre_ + (x - pre_) % per_)];
    }
    long long preperiod() const { return pre_; }
    long long period() const { return per_; }

  private:
    std::vector<int> g_;
    long long pre_ = 0, per_ = 1;
};

// ---- 1 手で「いくつかのゲームに分かれる」ときの Grundy 数 ----
// moves[v] = v の手の一覧。1 つの手は「行き先の状態の列」で、その手の値は行き先の Grundy 数の XOR
//   （山を 2 つに分ける：{a, b}、山を取り除くだけ：{}（値 0））。閉路があれば止める  O(n + 行き先の総数)
inline std::vector<int> grundy_dag_multi(int n, const std::vector<std::vector<std::vector<int>>>& moves) {
    if (n < 0 || (int)moves.size() != n) math_internal::fail("grundy_dag_multi", "moves の大きさは n");
    std::vector<int> g((size_t)n, -1), st;
    std::vector<char> state((size_t)n, 0);  // 0 = まだ、1 = 調べている途中、2 = 決まった
    std::vector<size_t> it_move((size_t)n, 0), it_pos((size_t)n, 0);
    std::vector<char> seen;
    for (int s = 0; s < n; ++s) {
        if (state[(size_t)s]) continue;
        st.push_back(s), state[(size_t)s] = 1;
        while (!st.empty()) {
            const int v = st.back();
            const auto& mv = moves[(size_t)v];
            bool pushed = false;
            while (it_move[(size_t)v] < mv.size()) {  // 次に調べる行き先
                const auto& to = mv[it_move[(size_t)v]];
                if (it_pos[(size_t)v] == to.size()) {
                    ++it_move[(size_t)v], it_pos[(size_t)v] = 0;
                    continue;
                }
                const int w = to[it_pos[(size_t)v]++];
                if (w < 0 || w >= n) math_internal::fail("grundy_dag_multi", "行き先の番号は 0〜n−1");
                if (state[(size_t)w] == 1) math_internal::fail("grundy_dag_multi", "閉路がある（同じ状態に戻れる）");
                if (state[(size_t)w] == 0) {
                    st.push_back(w), state[(size_t)w] = 1, pushed = true;
                    break;
                }
            }
            if (pushed) continue;
            seen.assign(mv.size() + 1, 0);
            for (const auto& to : mv) {
                int x = 0;
                for (int w : to) x ^= g[(size_t)w];
                if ((size_t)x <= mv.size()) seen[(size_t)x] = 1;
            }
            int m = 0;
            while (seen[(size_t)m]) ++m;
            g[(size_t)v] = m, state[(size_t)v] = 2, st.pop_back();
        }
    }
    return g;
}

// ---- 閉路があってもよいゲームの勝ち負け（後退解析）----
// 状態 0〜n−1、moves[v] = v から 1 手で行ける状態。動けない状態は手番の人の負け。
// result[v] = GAME_WIN（手番の人が勝つ）/ GAME_LOSE / GAME_DRAW（どちらも勝てず続く）、
// steps[v] = 決着までの手数（勝つ方は最短、負ける方は最長で。引き分けは −1）  O(n + 辺の数)
enum GameResultKind : int { GAME_LOSE = 0, GAME_WIN = 1, GAME_DRAW = 2 };
struct GameOutcome {
    std::vector<int> result, steps;
};
inline GameOutcome game_outcome(int n, const std::vector<std::vector<int>>& moves) {
    if (n < 0 || (int)moves.size() != n) math_internal::fail("game_outcome", "moves の大きさは n");
    std::vector<int> cnt((size_t)n), head((size_t)n + 1, 0), radj;  // 逆向きの辺（CSR）
    for (int v = 0; v < n; ++v) {
        cnt[(size_t)v] = (int)moves[(size_t)v].size();
        for (int w : moves[(size_t)v]) {
            if (w < 0 || w >= n) math_internal::fail("game_outcome", "行き先の番号は 0〜n−1");
            ++head[(size_t)w + 1];
        }
    }
    for (int v = 0; v < n; ++v) head[(size_t)v + 1] += head[(size_t)v];
    radj.resize((size_t)head[(size_t)n]);
    {
        std::vector<int> pos(head.begin(), head.end() - 1);
        for (int v = 0; v < n; ++v)
            for (int w : moves[(size_t)v]) radj[(size_t)pos[(size_t)w]++] = v;
    }
    GameOutcome r;
    r.result.assign((size_t)n, GAME_DRAW), r.steps.assign((size_t)n, -1);
    std::vector<int> q;
    q.reserve((size_t)n);
    for (int v = 0; v < n; ++v)
        if (cnt[(size_t)v] == 0) r.result[(size_t)v] = GAME_LOSE, r.steps[(size_t)v] = 0, q.push_back(v);
    for (size_t h = 0; h < q.size(); ++h) {  // 手数の小さい順に決まる
        const int v = q[h];
        for (int i = head[(size_t)v]; i < head[(size_t)v + 1]; ++i) {
            const int u = radj[(size_t)i];
            if (r.result[(size_t)u] != GAME_DRAW) continue;
            if (r.result[(size_t)v] == GAME_LOSE) {  // 負けの状態に行ければ勝ち
                r.result[(size_t)u] = GAME_WIN, r.steps[(size_t)u] = r.steps[(size_t)v] + 1, q.push_back(u);
            } else if (--cnt[(size_t)u] == 0) {  // どこへ行っても相手の勝ち
                r.result[(size_t)u] = GAME_LOSE, r.steps[(size_t)u] = r.steps[(size_t)v] + 1, q.push_back(u);
            }
        }
    }
    return r;
}
#endif  // GG_GAME
