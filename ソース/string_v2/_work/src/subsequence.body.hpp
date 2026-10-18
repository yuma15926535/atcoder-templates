#ifndef GG_SUBSEQUENCE
#define GG_SUBSEQUENCE
// ---- 部分列（文字列第二世代）  次の出現の表・部分列の判定・何文字まで入るか・辞書順最小の部分列・部分列でない最短・部分列の数 ----
// 「部分列」は飛び飛びに選んだもの（順番はそのまま）。連続する「部分文字列」は 02_pattern_match・03_suffix_array
// 数える関数は型 T を受け取る（math_v2 の ModInt998 など。long long なら溢れたら止める）

// ---- 次 / 前の出現を O(1) で（中は 00_str_base の CharPrefix。メモリ (N+1)σ 個の int） ----
struct NextTable {
    NextTable() = default;
    template <class Seq>
    explicit NextTable(const Seq& s, const Alphabet& A = Alphabet::lower()) : cp_(s, A) {}
    int size() const { return cp_.size(); }
    template <class C>
    int next(int i, C c) const { return cp_.next(i, c); }  // i 以上で最初の c（なければ -1。0 ≤ i ≤ N）  O(1)
    template <class C>
    int prev(int i, C c) const { return cp_.prev(i, c); }  // i より前で最後の c（なければ -1）  O(1)
    template <class T>
    int match_length(const T& t0, int from = 0) const {  // s[from..) に t の先頭が何文字まで部分列として入るか（貪欲）  O(|t|)
        const auto& t = str_internal::seq(t0);
        int pos = from, k = 0;
        if (from < 0 || from > size()) str_internal::fail("NextTable::match_length", "from は 0〜N");
        for (; k < (int)t.size(); ++k) {
            const int q = cp_.next(pos, t[k]);
            if (q < 0) break;
            pos = q + 1;
        }
        return k;
    }
    template <class T>
    int match_end(const T& t0, int from = 0) const {  // s[from..) に t を貪欲に当てはめたときの、最後の文字の次の位置（入らなければ -1）  O(|t|)
        const auto& t = str_internal::seq(t0);
        if (from < 0 || from > size()) str_internal::fail("NextTable::match_end", "from は 0〜N");
        int pos = from;
        for (std::size_t k = 0; k < t.size(); ++k) {
            const int q = cp_.next(pos, t[k]);
            if (q < 0) return -1;
            pos = q + 1;
        }
        return pos;
    }
    template <class T>
    bool is_subsequence(const T& t, int from = 0) const { return match_end(t, from) >= 0; }  // t が s[from..) の部分列か  O(|t|)
    const CharPrefix& char_prefix() const { return cp_; }  // 区間の文字の数なども使える

   private:
    CharPrefix cp_;
};

// ---- 先頭から・末尾から何文字まで入るか（ABC324 E） ----
template <class S, class T>
std::vector<int> prefix_match_lengths(const S& s0, const T& t0) {  // r[i] = t の先頭何文字までが s[0, i) の部分列になるか（長さ |s| + 1）  O(|s|)
    const auto& s = str_internal::seq(s0);
    const auto& t = str_internal::seq(t0);
    std::vector<int> r(s.size() + 1, 0);
    int k = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (k < (int)t.size() && str_internal::code_of(s[i]) == str_internal::code_of(t[k])) ++k;
        r[i + 1] = k;
    }
    return r;
}
template <class S, class T>
std::vector<int> suffix_match_lengths(const S& s0, const T& t0) {  // r[i] = t の末尾何文字までが s[i, N) の部分列になるか（長さ |s| + 1）  O(|s|)
    const auto& s = str_internal::seq(s0);
    const auto& t = str_internal::seq(t0);
    const int n = (int)s.size(), m = (int)t.size();
    std::vector<int> r(n + 1, 0);
    int k = 0;
    for (int i = n - 1; i >= 0; --i) {
        if (k < m && str_internal::code_of(s[i]) == str_internal::code_of(t[m - 1 - k])) ++k;
        r[i] = k;
    }
    return r;
}

// ---- 辞書順最小 / 最大の部分列（長さ k） ----
namespace str_internal {
template <class Seq>
seq_t<Seq> extreme_subsequence(const Seq& s0, int k, bool want_min) {
    const auto& s = seq(s0);
    const int n = (int)s.size();
    if (k < 0 || k > n) fail(want_min ? "lex_smallest_subsequence" : "lex_largest_subsequence", "0 ≤ k ≤ N で使う");
    seq_t<Seq> st;
    for (int i = 0; i < n; ++i) {
        const long long x = code_of(s[i]);
        while (!st.empty() && (want_min ? code_of(st.back()) > x : code_of(st.back()) < x) && (int)st.size() - 1 + (n - i) >= k) st.pop_back();
        if ((int)st.size() < k) st.push_back(s[i]);
    }
    return st;
}
}  // namespace str_internal
template <class Seq>
str_internal::seq_t<Seq> lex_smallest_subsequence(const Seq& s, int k) { return str_internal::extreme_subsequence(s, k, true); }  // 典型90 006  O(N)
template <class Seq>
str_internal::seq_t<Seq> lex_largest_subsequence(const Seq& s, int k) { return str_internal::extreme_subsequence(s, k, false); }

// ---- 部分列にならない最短の文字列（その中で辞書順最小。ARC081 E） ----
template <class Seq>
std::vector<int> shortest_non_subsequence_ids(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 番号の列で返す  O(Nσ)
    const auto& s = str_internal::seq(s0);
    const int n = (int)s.size(), sg = A.size();
    std::vector<int> id(n);
    for (int i = 0; i < n; ++i) id[i] = A.id(s[i]);
    // dp[i] = s[i..) の部分列でない最短の長さ。nx[i][c] = i 以上で最初の c（なければ n）
    std::vector<int> dp(n + 2, 0), choice(n + 1, 0), last(sg, n);
    dp[n] = 1;  // 空の残りには長さ 1 の文字列が入らない
    std::vector<int> nxt((std::size_t)(n + 1) * sg);
    for (int c = 0; c < sg; ++c) nxt[(std::size_t)n * sg + c] = n;
    for (int i = n - 1; i >= 0; --i) {
        last[id[i]] = i;
        for (int c = 0; c < sg; ++c) nxt[(std::size_t)i * sg + c] = last[c];
    }
    for (int i = n - 1; i >= 0; --i) {
        int best = 1 << 30, bc = 0;
        for (int c = 0; c < sg; ++c) {
            const int q = nxt[(std::size_t)i * sg + c];
            const int v = q == n ? 1 : 1 + dp[q + 1];
            if (v < best) best = v, bc = c;
        }
        dp[i] = best, choice[i] = bc;
    }
    std::vector<int> r;
    for (int i = 0;;) {
        if (i >= n) {
            r.push_back(0);
            break;
        }
        const int c = choice[i], q = nxt[(std::size_t)i * sg + c];
        r.push_back(c);
        if (q == n) break;
        i = q + 1;
    }
    return r;
}
template <class Seq>
std::string shortest_non_subsequence(const Seq& s, const Alphabet& A = Alphabet::lower()) { return from_ids(shortest_non_subsequence_ids(s, A), A); }  // 文字列で返す

// ---- 部分列の数 ----
template <class T, class S, class P>
T count_subsequence_occurrences(const S& s0, const P& t0) {  // t が s の部分列として現れる数（位置の選び方の数）  O(|s||t|)
    const auto& s = str_internal::seq(s0);
    const auto t = str_internal::codes_of(t0);
    const int m = (int)t.size();
    std::vector<T> dp(m + 1, T(0));  // dp[j] = t の先頭 j 文字の現れ方の数（整数型で溢れたら over[j] に印。答えの dp[m] のときだけ止める）
    std::vector<char> over(m + 1, 0);
    dp[0] = T(1);
    for (std::size_t i = 0; i < s.size(); ++i) {
        const long long x = str_internal::code_of(s[i]);
        for (int j = m; j >= 1; --j)
            if (t[j - 1] == x) {
                bool o = over[j];
                str_internal::add_count(dp[j], o, dp[j - 1], over[j - 1]);
                over[j] = o;
            }
    }
    if (over[m]) str_internal::count_overflow("count_subsequence_occurrences");
    return dp[m];
}
template <class T, class Seq>
T count_distinct_subsequences(const Seq& s0, bool include_empty = false, int min_step = 1) {  // 異なる部分列の数（選ぶ位置の間を min_step 以上あける。ABC214 F は 2）  O(N log N)（値の圧縮）
    if (min_step < 1) str_internal::fail("count_distinct_subsequences", "min_step ≥ 1 で使う");
    const CompressedValues cv = compress_values(s0);
    const int n = (int)cv.ids.size(), g = min_step;
    std::vector<int> last(cv.size(), -1);
    std::vector<T> f(n, T(0)), P(n + 1, T(0));  // f[i] = 貪欲に左から当てはめると i で終わる異なる部分列の数、P = f の累積和
    T total = include_empty ? T(1) : T(0);
    for (int i = 0; i < n; ++i) {
        const int prev = last[cv.ids[i]];
        last[cv.ids[i]] = i;
        T v = prev < 0 ? T(1) : T(0);
        const int lo = prev < 0 ? 0 : std::max(0, prev - g + 1), hi = i - g;  // 1 つ前に選んだ位置 j の範囲 [lo, hi]
        if (lo <= hi) v = str_internal::checked_add(v, P[hi + 1] - P[lo], "count_distinct_subsequences");
        f[i] = v;
        P[i + 1] = str_internal::checked_add(P[i], v, "count_distinct_subsequences");
        total = str_internal::checked_add(total, v, "count_distinct_subsequences");
    }
    return total;
}
template <class T, class Seq>
std::vector<T> distinct_subsequences_by_length(const Seq& s0) {  // r[len] = 長さ len の異なる部分列の数（len = 0〜N。r[0] = 1）  O(N^2)
    const CompressedValues cv = compress_values(s0);
    const int n = (int)cv.ids.size();
    std::vector<int> prev(n), last(cv.size(), -1);
    for (int i = 0; i < n; ++i) prev[i] = last[cv.ids[i]], last[cv.ids[i]] = i;
    std::vector<T> r(n + 1, T(0)), f(n, T(0)), P(n + 1, T(0));
    r[0] = T(1);
    for (int i = 0; i < n; ++i) f[i] = prev[i] < 0 ? T(1) : T(0);  // 長さ 1
    for (int len = 1; len <= n; ++len) {
        T sum = T(0);
        for (int i = 0; i < n; ++i) sum = str_internal::checked_add(sum, f[i], "distinct_subsequences_by_length");
        r[len] = sum;
        if (len == n) break;
        P[0] = T(0);
        for (int i = 0; i < n; ++i) P[i + 1] = str_internal::checked_add(P[i], f[i], "distinct_subsequences_by_length");
        for (int i = 0; i < n; ++i) {  // 次の長さ：j ∈ [max(0, prev), i - 1]
            const int lo = std::max(0, prev[i]), hi = i - 1;
            f[i] = lo <= hi ? P[hi + 1] - P[lo] : T(0);
        }
    }
    return r;
}
#endif  // GG_SUBSEQUENCE
