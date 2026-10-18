#ifndef GG_LCS_EDIT
#define GG_LCS_EDIT
// ---- 最長共通部分列（LCS）と編集距離（文字列第二世代） ----
// LCS は飛び飛びに選んだ「部分列」。連続する「部分文字列」の最長共通は 03_suffix_array の longest_common_substring
namespace str_internal {
// LCS のビット並列（Hyyrö）：b の位置をビットにして a を 1 文字ずつ読む。rows に各行の V を残す（復元用）
struct LcsBits {
    int n = 0, m = 0, W = 0;
    std::vector<unsigned long long> V;  // 最後の行（keep_rows なら (n+1) 行ぶん）
    template <class A, class B>
    LcsBits(const A& a0, const B& b0, bool keep_rows, const char* func = "lcs_restore") {
        const auto a = codes_of(a0), b = codes_of(b0);
        n = (int)a.size(), m = (int)b.size(), W = (m + 63) / 64;
        if (keep_rows && (long long)(n + 1) * W > 30000000LL) fail(func, "|a| × |b| が大きすぎる（約 1.9×10^9 まで）。長さだけなら lcs_length");
        // b の値ごとの位置のビット
        std::vector<long long> vals(b);
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
        std::vector<unsigned long long> M(vals.size() * (std::size_t)W, 0);
        auto vid = [&](long long x) {
            const int k = upper_index(vals, x) - 1;
            return (k >= 0 && vals[k] == x) ? k : -1;
        };
        for (int j = 0; j < m; ++j) M[(std::size_t)vid(b[j]) * W + j / 64] |= 1ULL << (j % 64);
        const std::size_t rows = keep_rows ? (std::size_t)n + 1 : 1;
        V.assign(rows * W, ~0ULL);
        std::vector<unsigned long long> cur(W, ~0ULL);
        for (int i = 0; i < n; ++i) {
            const int k = vid(a[i]);
            if (k >= 0) {
                const unsigned long long* mk = M.data() + (std::size_t)k * W;
                unsigned long long carry = 0;
                for (int w = 0; w < W; ++w) {
                    // V' = (V + U) | (V - U)、U = V & M（U は V の一部なので V - U = V & ~M。足し算だけ繰り上がりを運ぶ）
                    const unsigned long long v = cur[w], u = v & mk[w];
                    const unsigned long long s1 = v + u, c1 = s1 < v;
                    const unsigned long long s2 = s1 + carry, c2 = s2 < s1;
                    carry = c1 | c2;
                    cur[w] = s2 | (v & ~mk[w]);
                }
            }
            if (keep_rows) std::copy(cur.begin(), cur.end(), V.begin() + (std::size_t)(i + 1) * W);
        }
        if (!keep_rows) std::copy(cur.begin(), cur.end(), V.begin());
    }
    // LCS(a[0, i), b[0, j)) = j - （行 i の下 j ビットの 1 の数）
    int value(int i, int j) const {
        const unsigned long long* r = V.data() + (std::size_t)i * W;
        int ones = 0;
        for (int w = 0; w < j / 64; ++w) ones += __builtin_popcountll(r[w]);
        if (j % 64) ones += __builtin_popcountll(r[j / 64] & ((1ULL << (j % 64)) - 1));
        return j - ones;
    }
};
// LCS の対応する位置の組（a の位置, b の位置）を前から
template <class A, class B>
std::vector<std::pair<int, int>> lcs_pairs(const A& a0, const B& b0, const char* func) {
    const auto a = codes_of(a0), b = codes_of(b0);
    const LcsBits L(a, b, true, func);
    std::vector<std::pair<int, int>> r;
    int i = L.n, j = L.m, cur = L.value(i, j);
    while (i > 0 && j > 0 && cur > 0) {
        if (a[i - 1] == b[j - 1]) {
            r.push_back({i - 1, j - 1}), --i, --j, --cur;
        } else if (L.value(i - 1, j) == cur) {
            --i;
        } else {
            --j;
        }
    }
    std::reverse(r.begin(), r.end());
    return r;
}
}  // namespace str_internal

template <class A, class B>
int lcs_length(const A& a, const B& b) {  // 最長共通部分列（飛び飛び）の長さ。ビット並列  O(|a||b|/64)
    const str_internal::LcsBits L(a, b, false);
    return L.value(0, L.m);
}
template <class A, class B>
str_internal::seq_t<A> lcs_restore(const A& a0, const B& b) {  // 最長共通部分列を 1 つ（a の型で）。|a| × |b| が約 1.9×10^9 まで（両方 4×10^4 くらい）  O(|a||b|/64)
    const auto& a = str_internal::seq(a0);
    str_internal::seq_t<A> r;
    for (auto [i, j] : str_internal::lcs_pairs(a, b, "lcs_restore")) r.push_back(a[i]), (void)j;
    return r;
}
template <class A, class B>
str_internal::seq_t<A> shortest_common_supersequence(const A& a0, const B& b0) {  // a と b の両方を部分列に含む最短の列（長さ |a| + |b| - LCS）
    const auto& a = str_internal::seq(a0);
    const auto& b = str_internal::seq(b0);
    str_internal::seq_t<A> r;
    int pa = 0, pb = 0;
    using E = typename str_internal::seq_t<A>::value_type;
    auto conv = [&](long long v) {  // b の値を a の要素の型に（収まらなければ止める）
        const E e = (E)v;
        if (str_internal::code_of(e) != v) str_internal::fail("shortest_common_supersequence", "b の値 " + std::to_string(v) + " が a の要素の型に収まらない（a を数の列にする）");
        return e;
    };
    for (auto [i, j] : str_internal::lcs_pairs(a, b, "shortest_common_supersequence")) {
        while (pa < i) r.push_back(a[pa++]);
        while (pb < j) r.push_back(conv(str_internal::code_of(b[pb++])));
        r.push_back(a[i]), ++pa, ++pb;
    }
    while (pa < (int)a.size()) r.push_back(a[pa++]);
    while (pb < (int)b.size()) r.push_back(conv(str_internal::code_of(b[pb++])));
    return r;
}

// ---- 編集距離（a を b にする。挿入 ins・削除 del・置換 sub の費用（0 以上）。既定はすべて 1） ----
namespace str_internal {
constexpr long long ed_inf = std::numeric_limits<long long>::max();
inline long long ed_add(long long x, long long y) { return x > ed_inf - y ? ed_inf : x + y; }  // 0 以上どうしの和（上で止める）
inline void ed_check(const char* func, long long ins, long long del, long long sub) {
    if (ins < 0 || del < 0 || sub < 0) fail(func, "費用 ins・del・sub は 0 以上で使う");
}
}  // namespace str_internal
template <class A, class B>
long long edit_distance(const A& a0, const B& b0, long long ins = 1, long long del = 1, long long sub = 1) {  // O(|a||b|) 時間・O(|b|) メモリ
    using str_internal::ed_add;
    str_internal::ed_check("edit_distance", ins, del, sub);
    const auto a = str_internal::codes_of(a0), b = str_internal::codes_of(b0);
    const int n = (int)a.size(), m = (int)b.size();
    std::vector<long long> dp(m + 1);
    dp[0] = 0;
    for (int j = 1; j <= m; ++j) dp[j] = ed_add(dp[j - 1], ins);
    for (int i = 1; i <= n; ++i) {
        long long diag = dp[0];
        dp[0] = ed_add(dp[0], del);
        for (int j = 1; j <= m; ++j) {
            const long long up = dp[j];
            dp[j] = std::min({ed_add(up, del), ed_add(dp[j - 1], ins), a[i - 1] == b[j - 1] ? diag : ed_add(diag, sub)});
            diag = up;
        }
    }
    if (dp[m] == str_internal::ed_inf) str_internal::fail("edit_distance", "答えが long long に収まらない");
    return dp[m];
}
struct EditOp {
    char type;  // '=' そのまま（a[i] == b[j]）/ 'S' a[i] を b[j] に置換 / 'D' a[i] を削除 / 'I' b[j] を挿入（a[i] の前に）
    int i, j;   // a・b での位置（D の j・I の i は「今いる位置」）
};
template <class A, class B>
std::vector<EditOp> edit_script(const A& a0, const B& b0, long long ins = 1, long long del = 1, long long sub = 1) {  // 最小の費用の操作の列（前から）  O(|a||b|) 時間・メモリ
    using str_internal::ed_add;
    str_internal::ed_check("edit_script", ins, del, sub);
    const auto a = str_internal::codes_of(a0), b = str_internal::codes_of(b0);
    const int n = (int)a.size(), m = (int)b.size();
    if ((long long)(n + 1) * (m + 1) > 200000000LL) str_internal::fail("edit_script", "|a| × |b| が 2×10^8 を超える（距離だけなら edit_distance）");
    std::vector<unsigned char> from((std::size_t)(n + 1) * (m + 1));  // 0 = 対角, 1 = 上（削除）, 2 = 左（挿入）
    std::vector<long long> prev(m + 1), cur(m + 1);
    prev[0] = 0, from[0] = 2;
    for (int j = 1; j <= m; ++j) prev[j] = ed_add(prev[j - 1], ins), from[j] = 2;
    for (int i = 1; i <= n; ++i) {
        cur[0] = ed_add(prev[0], del), from[(std::size_t)i * (m + 1)] = 1;
        for (int j = 1; j <= m; ++j) {
            long long best = a[i - 1] == b[j - 1] ? prev[j - 1] : ed_add(prev[j - 1], sub);
            unsigned char f = 0;
            if (ed_add(prev[j], del) < best) best = ed_add(prev[j], del), f = 1;
            if (ed_add(cur[j - 1], ins) < best) best = ed_add(cur[j - 1], ins), f = 2;
            cur[j] = best, from[(std::size_t)i * (m + 1) + j] = f;
        }
        std::swap(prev, cur);
    }
    if (prev[m] == str_internal::ed_inf) str_internal::fail("edit_script", "費用の合計が long long に収まらない");
    std::vector<EditOp> ops;
    for (int i = n, j = m; i > 0 || j > 0;) {
        const unsigned char f = from[(std::size_t)i * (m + 1) + j];
        if (i > 0 && j > 0 && f == 0) ops.push_back({a[i - 1] == b[j - 1] ? '=' : 'S', i - 1, j - 1}), --i, --j;
        else if (i > 0 && (j == 0 || f == 1)) ops.push_back({'D', i - 1, j}), --i;
        else ops.push_back({'I', i, j - 1}), --j;
    }
    std::reverse(ops.begin(), ops.end());
    return ops;
}
template <class A, class B>
int edit_distance_bounded(const A& a0, const B& b0, int k) {  // すべて費用 1 の編集距離。k を超えるなら -1  O((|a| + |b|) k)
    const auto a = str_internal::codes_of(a0), b = str_internal::codes_of(b0);
    const int n = (int)a.size(), m = (int)b.size();
    if (k < 0) str_internal::fail("edit_distance_bounded", "k ≥ 0 で使う");
    if (std::abs(n - m) > k) return -1;
    k = std::min(k, std::max(n, m));  // 距離は max(|a|, |b|) 以下なので、それより大きい k は同じ
    const int W = 2 * k + 1, over = k + 1;  // over = 「k を超えた」
    std::vector<int> prev(W, over), cur(W, over);  // 列 j = i + (d - k)
    for (int d = 0; d < W; ++d) {
        const int j = d - k;
        if (j >= 0 && j <= m) prev[d] = j;
    }
    for (int i = 1; i <= n; ++i) {
        for (int d = 0; d < W; ++d) {
            const int j = i + d - k;
            int v = over;
            if (j >= 0 && j <= m) {
                if (j == 0) v = i;
                else {
                    v = prev[d] + (a[i - 1] == b[j - 1] ? 0 : 1);         // (i-1, j-1) は前の行の同じ d
                    if (d + 1 < W) v = std::min(v, prev[d + 1] + 1);      // (i-1, j)：削除
                    if (d > 0) v = std::min(v, cur[d - 1] + 1);           // (i, j-1)：挿入
                }
            }
            cur[d] = std::min(v, over);
        }
        std::swap(prev, cur);
    }
    const int r = prev[m - n + k];
    return r <= k ? r : -1;
}
template <class A, class B>
bool within_one_edit(const A& a0, const B& b0) {  // 挿入・削除・置換を 1 回以下で一致するか  O(N)
    const auto a = str_internal::codes_of(a0), b = str_internal::codes_of(b0);
    const int n = (int)a.size(), m = (int)b.size();
    if (std::abs(n - m) > 1) return false;
    int i = 0;
    while (i < n && i < m && a[i] == b[i]) ++i;
    if (n == m) {
        int j = n - 1;
        while (j > i && a[j] == b[j]) --j;
        return j <= i;
    }
    const auto& L = n > m ? a : b;
    const auto& S = n > m ? b : a;
    for (int j = i; j < (int)S.size(); ++j)
        if (L[j + 1] != S[j]) return false;
    return true;
}
#endif  // GG_LCS_EDIT
