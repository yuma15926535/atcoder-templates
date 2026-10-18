#ifndef GG_LYNDON
#define GG_LYNDON
// ---- Lyndon 分解・Lyndon 配列・極大な繰り返し（runs）（文字列第二世代） ----
// Lyndon 語 = 自分のどの真の接尾辞よりも辞書順で小さい（空でない）文字列。辞書順は std::string と同じ（接頭辞は小さい）
// Lyndon 配列と runs は中で 03_suffix_array（同梱済み）を使う。数の列は値の大小で比べる（unsigned long long は 2^63 未満）
namespace str_internal {
inline std::vector<int> lyndon_from_rank(const std::vector<int>& rank) {  // λ[i] = （i より後で接尾辞が小さい最初の位置）− i
    const int n = (int)rank.size();
    std::vector<int> lam(n), st;
    st.reserve(n);
    for (int i = n - 1; i >= 0; --i) {
        while (!st.empty() && rank[st.back()] > rank[i]) st.pop_back();
        lam[i] = (st.empty() ? n : st.back()) - i;
        st.push_back(i);
    }
    return lam;
}
inline std::vector<int> rank_of(const std::vector<int>& sa) {  // 接尾辞配列の逆
    std::vector<int> r(sa.size());
    for (std::size_t k = 0; k < sa.size(); ++k) r[sa[k]] = (int)k;
    return r;
}
}  // namespace str_internal

template <class Seq>
std::vector<StrSpan> lyndon_factorization(const Seq& s0) {  // s = l1 l2 … lk（どれも Lyndon 語で l1 ≥ l2 ≥ … ≥ lk）。Duval  O(N)
    const auto s = str_internal::codes_of(s0);
    const int n = (int)s.size();
    std::vector<StrSpan> r;
    for (int i = 0; i < n;) {
        int j = i + 1, k = i;
        while (j < n && s[k] <= s[j]) {
            if (s[k] < s[j]) k = i;
            else ++k;
            ++j;
        }
        while (i <= k) r.push_back({i, i + j - k}), i += j - k;
    }
    return r;
}
template <class Seq>
bool is_lyndon(const Seq& s) {  // Lyndon 語か（空なら false）  O(N)
    const auto f = lyndon_factorization(s);
    return f.size() == 1;
}
template <class Seq>
std::vector<int> min_suffix_of_prefixes(const Seq& s0) {  // r[i] = s[0, i) の辞書順最小の接尾辞の開始位置（i = 1〜N。r[0] = 0）  O(N)
    const auto s = str_internal::codes_of(s0);
    const int n = (int)s.size();
    std::vector<int> ans(n + 1, 0);
    int done = 0;
    auto set = [&](int len, int v) {
        if (len > done) ans[len] = v, done = len;
    };
    for (int i = 0; i < n;) {  // Duval の途中の形 s[i, j) = w^q w'（w は Lyndon 語、w' は w の真の接頭辞）から読む
        int j = i + 1, k = i;
        set(j, i);
        while (j < n && s[k] <= s[j]) {
            if (s[k] < s[j]) k = i;
            else ++k;
            ++j;
            const int wl = j - k, q = (j - i) / wl, rem = (j - i) % wl;
            set(j, rem == 0 ? j - wl : ans[i + rem] + q * wl);
        }
        while (i <= k) i += j - k;
    }
    return ans;
}
template <class Seq>
std::vector<int> lyndon_array(const Seq& s) {  // λ[i] = i から始まる最長の Lyndon 語の長さ  O(N)（接尾辞配列。数の列は O(N log N)）
    return str_internal::lyndon_from_rank(str_internal::rank_of(make_suffix_array(s)));
}

// ---- 極大な繰り返し（runs）：s[l, r) が周期 p（最小）をもち 2p ≤ r − l で、左右に伸ばせない ----
struct Repetition {
    int l = 0, r = 0, period = 0;
    int len() const { return r - l; }  // 長さ
    friend bool operator==(const Repetition& a, const Repetition& b) { return a.l == b.l && a.r == b.r && a.period == b.period; }
    friend bool operator<(const Repetition& a, const Repetition& b) { return a.l != b.l ? a.l < b.l : (a.r != b.r ? a.r < b.r : a.period < b.period); }
    friend std::ostream& operator<<(std::ostream& os, const Repetition& x) { return os << "{[" << x.l << ", " << x.r << "), p = " << x.period << "}"; }  // デバッグ用
};
namespace str_internal {
template <class S>
std::vector<Repetition> runs_impl(const S& f, const S& rv, const S& inv) {  // f と、その反転 rv、文字の順を逆にした inv から
    const int n = (int)f.size();
    std::vector<Repetition> res;
    if (n < 2) return res;
    const SuffixArray F(f), R(rv);  // 前向き・後ろ向きの LCP（疎テーブル）
    const std::vector<int> rank_inv = rank_of(make_suffix_array(inv));
    for (int o = 0; o < 2; ++o) {  // 2 つの順それぞれの Lyndon 根から、左右に伸ばす
        const auto lam = lyndon_from_rank(o == 0 ? F.rank : rank_inv);
        for (int i = 0; i < n; ++i) {
            const int j = i + lam[i], p = lam[i];
            const int r = j + (j < n ? F.lcp(i, j) : 0);
            const int l = i - (i > 0 ? R.lcp(n - i, n - j) : 0);
            if (r - l >= 2 * p) res.push_back({l, r, p});
        }
    }
    std::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    return res;
}
}  // namespace str_internal
template <class Seq>
std::vector<Repetition> maximal_repetitions(const Seq& s0) {  // runs の一覧（l, r の昇順。個数は N 未満）  O(N log N)（接尾辞配列 3 本。N = 10^6 で約 0.8 秒・約 220MB）
    if constexpr (std::is_same<str_internal::elem_t<Seq>, char>::value) {  // 文字列は char のまま（接尾辞配列の速い道）
        const auto& t = str_internal::seq(s0);
        const std::string f(t.begin(), t.end()), rv(f.rbegin(), f.rend());
        std::string inv(f);
        for (auto& c : inv) c = (char)(255 - (unsigned char)c);  // 文字の順を逆に
        return str_internal::runs_impl(f, rv, inv);
    } else {
        const auto s = str_internal::codes_of(s0);
        std::vector<long long> rv(s.rbegin(), s.rend()), inv(s);
        for (auto& x : inv) x = ~x;  // 値の順を逆に
        return str_internal::runs_impl(s, rv, inv);
    }
}
template <class Seq>
long long count_square_substrings(const Seq& s) {  // xx の形の部分文字列（空でない x。位置が違えば別）の数  O(N log N)
    long long c = 0;
    for (const Repetition& x : maximal_repetitions(s))
        for (long long q = x.period; 2 * q <= x.len(); q += x.period) c += x.len() - 2 * q + 1;
    return c;
}
#endif  // GG_LYNDON
