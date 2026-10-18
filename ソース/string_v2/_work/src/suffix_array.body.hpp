#ifndef GG_SUFFIX_ARRAY
#define GG_SUFFIX_ARRAY
// ---- 接尾辞配列（文字列第二世代）  SA-IS・LCP（Kasai）・O(1) の LCP と比較・出現・異なる部分文字列・k 番目・複数の文字列 ----
// ACL の suffix_array / lcp_array と名前を変えている（make_suffix_array / make_lcp_array。using namespace atcoder; と一緒でも衝突しない）
namespace str_internal {
inline std::vector<int> sa_naive(const std::vector<int>& s) {
    const int n = (int)s.size();
    std::vector<int> sa(n);
    std::iota(sa.begin(), sa.end(), 0);
    std::sort(sa.begin(), sa.end(), [&](int l, int r) {
        if (l == r) return false;
        while (l < n && r < n) {
            if (s[l] != s[r]) return s[l] < s[r];
            ++l, ++r;
        }
        return l == n;
    });
    return sa;
}
inline std::vector<int> sa_doubling(const std::vector<int>& s) {
    const int n = (int)s.size();
    std::vector<int> sa(n), rnk = s, tmp(n);
    std::iota(sa.begin(), sa.end(), 0);
    for (int k = 1; k < n; k *= 2) {
        auto cmp = [&](int x, int y) {
            if (rnk[x] != rnk[y]) return rnk[x] < rnk[y];
            const int rx = x + k < n ? rnk[x + k] : -1, ry = y + k < n ? rnk[y + k] : -1;
            return rx < ry;
        };
        std::sort(sa.begin(), sa.end(), cmp);
        tmp[sa[0]] = 0;
        for (int i = 1; i < n; ++i) tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
        std::swap(tmp, rnk);
    }
    return sa;
}
// SA-IS（s の値は 0〜upper）
inline std::vector<int> sa_is(const std::vector<int>& s, int upper) {
    const int n = (int)s.size();
    if (n == 0) return {};
    if (n == 1) return {0};
    if (n == 2) return s[0] < s[1] ? std::vector<int>{0, 1} : std::vector<int>{1, 0};
    if (n < 10) return sa_naive(s);
    if (n < 40) return sa_doubling(s);
    std::vector<int> sa(n);
    std::vector<bool> ls(n);  // true = S 型（vector<bool> の方が速い：induce で飛び飛びに読むので小さい方がよい）
    for (int i = n - 2; i >= 0; --i) ls[i] = (s[i] == s[i + 1]) ? ls[i + 1] : (s[i] < s[i + 1]);
    std::vector<int> sum_l(upper + 1), sum_s(upper + 1);
    for (int i = 0; i < n; ++i) {
        if (!ls[i]) ++sum_s[s[i]];
        else ++sum_l[s[i] + 1];
    }
    for (int i = 0; i <= upper; ++i) {
        sum_s[i] += sum_l[i];
        if (i < upper) sum_l[i + 1] += sum_s[i];
    }
    std::vector<int> buf(upper + 1);
    auto induce = [&](const std::vector<int>& lms) {
        std::fill(sa.begin(), sa.end(), -1);
        std::copy(sum_s.begin(), sum_s.end(), buf.begin());
        for (int d : lms)
            if (d != n) sa[buf[s[d]]++] = d;
        std::copy(sum_l.begin(), sum_l.end(), buf.begin());
        sa[buf[s[n - 1]]++] = n - 1;
        for (int i = 0; i < n; ++i) {
            const int v = sa[i];
            if (v >= 1 && !ls[v - 1]) sa[buf[s[v - 1]]++] = v - 1;
        }
        std::copy(sum_l.begin(), sum_l.end(), buf.begin());
        for (int i = n - 1; i >= 0; --i) {
            const int v = sa[i];
            if (v >= 1 && ls[v - 1]) sa[--buf[s[v - 1] + 1]] = v - 1;
        }
    };
    std::vector<int> lms_map(n + 1, -1), lms;
    int m = 0;
    for (int i = 1; i < n; ++i)
        if (!ls[i - 1] && ls[i]) lms_map[i] = m++;
    lms.reserve(m);
    for (int i = 1; i < n; ++i)
        if (!ls[i - 1] && ls[i]) lms.push_back(i);
    induce(lms);
    if (m) {
        std::vector<int> sorted_lms;
        sorted_lms.reserve(m);
        for (int v : sa)
            if (lms_map[v] != -1) sorted_lms.push_back(v);
        std::vector<int> rec_s(m);
        int rec_upper = 0;
        rec_s[lms_map[sorted_lms[0]]] = 0;
        for (int i = 1; i < m; ++i) {
            int l = sorted_lms[i - 1], r = sorted_lms[i];
            const int end_l = (lms_map[l] + 1 < m) ? lms[lms_map[l] + 1] : n;
            const int end_r = (lms_map[r] + 1 < m) ? lms[lms_map[r] + 1] : n;
            bool same = true;
            if (end_l - l != end_r - r) {
                same = false;
            } else {
                while (l < end_l) {
                    if (s[l] != s[r]) break;
                    ++l, ++r;
                }
                if (l == n || s[l] != s[r]) same = false;
            }
            if (!same) ++rec_upper;
            rec_s[lms_map[sorted_lms[i]]] = rec_upper;
        }
        const std::vector<int> rec_sa = sa_is(rec_s, rec_upper);
        for (int i = 0; i < m; ++i) sorted_lms[i] = lms[rec_sa[i]];
        induce(sorted_lms);
    }
    return sa;
}
inline std::vector<int> lcp_kasai(const std::vector<int>& s, const std::vector<int>& sa) {
    const int n = (int)s.size();
    if (n == 0) return {};
    std::vector<int> rnk(n), lcp(n - 1);
    for (int i = 0; i < n; ++i) rnk[sa[i]] = i;
    int h = 0;
    for (int i = 0; i < n; ++i) {
        if (h > 0) --h;
        if (rnk[i] == 0) continue;
        const int j = sa[rnk[i] - 1];
        while (j + h < n && i + h < n && s[j + h] == s[i + h]) ++h;
        lcp[rnk[i] - 1] = h;
    }
    return lcp;
}
// 列 → 0〜upper の整数の列（char は 0〜255。数は圧縮）
template <class Seq>
std::vector<int> sa_ints(const Seq& s0, int& upper) {
    const auto& s = seq(s0);
    std::vector<int> r(s.size());
    if constexpr (std::is_same<elem_t<Seq>, char>::value || std::is_same<elem_t<Seq>, signed char>::value || std::is_same<elem_t<Seq>, unsigned char>::value) {
        for (std::size_t i = 0; i < s.size(); ++i) r[i] = (int)(unsigned char)s[i];
        upper = 255;
    } else {
        const CompressedValues c = compress_values(s);
        r = c.ids;
        upper = std::max(0, c.size() - 1);
    }
    return r;
}
}  // namespace str_internal

template <class Seq>
std::vector<int> make_suffix_array(const Seq& s) {  // 接尾辞を辞書順に並べた開始位置（SA-IS）  O(N)（数の列は圧縮で O(N log N)）
    int upper = 0;
    const auto v = str_internal::sa_ints(s, upper);
    return str_internal::sa_is(v, upper);
}
template <class Seq>
std::vector<int> make_lcp_array(const Seq& s, const std::vector<int>& sa) {  // lcp[i] = sa[i] と sa[i+1] から始まる接尾辞の LCP（長さ N - 1）  O(N)
    int upper = 0;
    return str_internal::lcp_kasai(str_internal::sa_ints(s, upper), sa);
}

// 2 か所の同じ部分文字列（長さ len と 2 つの開始位置。ないときは len = 0, pos = -1）
//   longest_repeated 系：同じ文字列の中の 2 か所（pos_a < pos_b）/ longest_common_substring：a での位置 pos_a と b での位置 pos_b
struct CommonSubstring {
    int len = 0;
    int pos_a = -1, pos_b = -1;
};

struct SuffixArray {
    std::vector<int> sa;       // sa[k] = 辞書順で k 番目の接尾辞の開始位置
    std::vector<int> rank;     // rank[p] = 位置 p から始まる接尾辞の順位（sa の逆）
    std::vector<int> lcp_adj;  // lcp_adj[k] = sa[k] と sa[k+1] の接尾辞の LCP（長さ N - 1）

    SuffixArray() = default;
    template <class Seq>
    explicit SuffixArray(const Seq& s0) {
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        val_.resize(n_);
        for (int i = 0; i < n_; ++i) val_[i] = str_internal::code_of(s[i]);
        int upper = 0;
        const auto v = str_internal::sa_ints(s, upper);
        build(v, upper);
    }
    // 複数の文字列をまとめて（間に区切りを入れる。区切りはどの文字より小さく、互いに違う）。sa には区切りの位置（owner = -1）も入る
    template <class Seq>
    static SuffixArray of_strings(const std::vector<Seq>& vs) {
        SuffixArray R;
        const int K = (int)vs.size();
        long long total = K;
        for (const auto& x : vs) total += (long long)str_internal::seq(x).size();
        if (total > 2000000000LL) str_internal::fail("SuffixArray::of_strings", "長さの合計が大きすぎる");
        R.n_ = (int)total, R.multi_ = true;
        R.val_.resize(R.n_), R.lim_.resize(R.n_), R.owner_.resize(R.n_), R.start_.resize(K + 1);
        std::vector<long long> all;
        for (const auto& x : vs) {
            const auto& s = str_internal::seq(x);
            for (std::size_t i = 0; i < s.size(); ++i) all.push_back(str_internal::code_of(s[i]));
        }
        const CompressedValues cv = compress_values(all);
        std::vector<int> v(R.n_);
        int p = 0, q = 0;
        for (int t = 0; t < K; ++t) {
            const int len = (int)str_internal::seq(vs[t]).size();
            R.start_[t] = p;
            for (int i = 0; i < len; ++i, ++p, ++q) v[p] = K + cv.ids[q], R.val_[p] = all[q], R.lim_[p] = R.start_[t] + len, R.owner_[p] = t;
            v[p] = t, R.val_[p] = 0, R.lim_[p] = p, R.owner_[p] = -1, ++p;  // 区切り
        }
        R.start_[K] = p;
        R.build(v, K + std::max(0, cv.size() - 1));
        return R;
    }

    int size() const { return n_; }  // 全体の長さ（of_strings なら区切りも含む）
    int length_at(int p) const {  // 位置 p から始まる接尾辞の（自分の文字列の中での）長さ（0 ≤ p ≤ N。区切りと p = N は 0）
        check_pos("SuffixArray::length_at", p);
        return p == n_ ? 0 : lim(p) - p;
    }
    int owner(int p) const {  // 位置 p がどの文字列か（区切りは -1。0 ≤ p < N）
        check_index("SuffixArray::owner", p);
        return multi_ ? owner_[p] : 0;
    }
    int offset(int p) const {  // 位置 p がその文字列の何文字目か（区切りの位置では止める。owner(p) ≥ 0 を先に確かめる）
        check_index("SuffixArray::offset", p);
        if (!multi_) return p;
        if (owner_[p] < 0) str_internal::fail("SuffixArray::offset", "位置 " + std::to_string(p) + " は区切り（owner = -1）");
        return p - start_[owner_[p]];
    }
    int position(int t, int off) const {  // t 番目の文字列の off 文字目の位置（0 ≤ off ≤ その長さ。off = 長さ は直後の区切り）
        const int K = multi_ ? (int)start_.size() - 1 : 1;
        if (t < 0 || t >= K) str_internal::fail("SuffixArray::position", "文字列の番号 " + std::to_string(t) + " は 0〜" + std::to_string(K - 1) + " の外");
        const int len = multi_ ? start_[t + 1] - start_[t] - 1 : n_;
        if (off < 0 || off > len) str_internal::fail("SuffixArray::position", "off = " + std::to_string(off) + " は 0〜" + std::to_string(len) + " の外");
        return multi_ ? start_[t] + off : off;
    }
    long long value(int p) const {  // 位置 p の値（char は 0〜255。区切りは 0）
        check_index("SuffixArray::value", p);
        return val_[p];
    }

    // 位置 i, j から始まる接尾辞の LCP（0 ≤ i, j ≤ N）  O(1)（初回に疎テーブルを O(N log N) で作る）
    int lcp(int i, int j) const {
        check_pos("SuffixArray::lcp", i), check_pos("SuffixArray::lcp", j);
        if (i == n_ || j == n_) return 0;
        if (i == j) return length_at(i);
        int a = rank[i], b = rank[j];
        if (a > b) std::swap(a, b);
        return range_min(a, b);
    }
    // 部分文字列 s[l1, r1) と s[l2, r2) の辞書順（-1 / 0 / 1）  O(1)
    int compare(int l1, int r1, int l2, int r2) const {
        check_sub("SuffixArray::compare", l1, r1), check_sub("SuffixArray::compare", l2, r2);
        const int m1 = r1 - l1, m2 = r2 - l2, m = std::min(m1, m2);
        if (m == 0 || lcp(l1, l2) >= m) return m1 == m2 ? 0 : (m1 < m2 ? -1 : 1);
        return rank[l1] < rank[l2] ? -1 : 1;
    }
    // pat で始まる接尾辞の SA 上の区間 [l, r)（なければ空。pat が空なら全部 [0, N)）  O(M log N)
    template <class P>
    StrSpan find_range(const P& pat0) const {
        const auto& pat = str_internal::seq(pat0);
        std::vector<long long> p(pat.size());
        for (std::size_t i = 0; i < pat.size(); ++i) p[i] = str_internal::code_of(pat[i]);
        int lo = 0, hi = n_;  // 最初の「接尾辞 ≥ pat」
        while (lo < hi) {
            const int mid = (lo + hi) >> 1;
            if (cmp_pat(sa[mid], p) < 0) lo = mid + 1;
            else hi = mid;
        }
        int lo2 = lo, hi2 = n_;  // 最初の「接尾辞 > pat（pat で始まらない）」
        while (lo2 < hi2) {
            const int mid = (lo2 + hi2) >> 1;
            if (cmp_pat(sa[mid], p) <= 0) lo2 = mid + 1;
            else hi2 = mid;
        }
        return {lo, lo2};
    }
    template <class P>
    long long count(const P& pat) const {  // pat の出現回数（pat が空なら置ける位置の数：N + 1。of_strings では各文字列の長さ + 1 の和 = N）  O(M log N)
        if (str_internal::seq(pat).size() == 0) return empty_count();
        return find_range(pat).len();
    }
    template <class P>
    std::vector<int> find_all(const P& pat) const {  // pat の出現位置（昇順。pat が空なら 0〜N（of_strings では 0〜N-1））  O(M log N + 個数 log)
        if (str_internal::seq(pat).size() == 0) {
            std::vector<int> v((std::size_t)empty_count());
            for (std::size_t i = 0; i < v.size(); ++i) v[i] = (int)i;
            return v;
        }
        const StrSpan r = find_range(pat);
        std::vector<int> v(sa.begin() + r.l, sa.begin() + r.r);
        std::sort(v.begin(), v.end());
        return v;
    }
    long long count_occurrences(int l, int r) const {  // 部分文字列 s[l, r) の出現回数（自分も含む）  O(log N)
        check_sub("SuffixArray::count_occurrences", l, r);
        if (l == r) return empty_count();
        const StrSpan g = occurrence_range(l, r);
        return g.len();
    }
    StrSpan occurrence_range(int l, int r) const {  // s[l, r)（空でない）で始まる接尾辞の SA 上の区間  O(log N)
        check_sub("SuffixArray::occurrence_range", l, r);
        if (l == r) str_internal::fail("SuffixArray::occurrence_range", "空の部分文字列");
        const int m = r - l, k = rank[l];
        ensure_table();
        int lo = 0, hi = k;  // 最小の a で、a..k の LCP がすべて m 以上
        while (lo < hi) {
            const int mid = (lo + hi) >> 1;
            if (range_min(mid, k) >= m) hi = mid;
            else lo = mid + 1;
        }
        int lo2 = k, hi2 = n_ - 1;  // 最大の b で、k..b の LCP がすべて m 以上
        while (lo2 < hi2) {
            const int mid = (lo2 + hi2 + 1) >> 1;
            if (range_min(k, mid) >= m) lo2 = mid;
            else hi2 = mid - 1;
        }
        return {lo, lo2 + 1};
    }

    // ---- 異なる部分文字列 ----
    long long count_distinct_substrings() const {  // 異なる部分文字列（空でない）の数  O(1)（初回 O(N)）
        ensure_newcnt();
        return newcnt_.back();
    }
    std::vector<long long> distinct_substrings_by_length() const {  // r[len] = 長さ len の異なる部分文字列の数（len = 0〜最長）  O(N)
        int mx = 0;
        for (int p = 0; p < n_; ++p) mx = std::max(mx, length_at(p));
        std::vector<long long> d(mx + 2, 0);
        for (int k = 0; k < n_; ++k) {
            const int lo = (k ? lcp_adj[k - 1] : 0) + 1, hi = length_at(sa[k]);
            if (lo <= hi) ++d[lo], --d[hi + 1];
        }
        std::vector<long long> r(mx + 1, 0);
        for (int t = 1; t <= mx; ++t) d[t] += d[t - 1], r[t] = d[t];
        return r;
    }
    StrSpan kth_substring(long long k) const {  // 異なる部分文字列の辞書順 k 番目（0 から）  O(log N)
        ensure_newcnt();
        if (k < 0 || k >= newcnt_.back()) str_internal::fail("SuffixArray::kth_substring", "k は 0〜（異なる部分文字列の数 - 1）");
        const int idx = str_internal::upper_index(newcnt_, k) - 1;
        const int len = (idx ? lcp_adj[idx - 1] : 0) + (int)(k - newcnt_[idx]) + 1;
        return {sa[idx], sa[idx] + len};
    }
    long long substring_rank(int l, int r) const {  // 異なる部分文字列のうち s[l, r)（空でない）より小さいものの数  O(log N)
        const StrSpan g = occurrence_range(l, r);
        ensure_newcnt();
        const int k0 = g.l;
        return newcnt_[k0] + (r - l - 1 - (k0 ? lcp_adj[k0 - 1] : 0));
    }
    long long count_substrings_with_dup() const {  // 部分文字列（位置が違えば別）の数 = Σ 長さ  O(1)（初回 O(N)）
        ensure_blocks();
        return bsum_.back();
    }
    StrSpan kth_substring_with_dup(long long k) const {  // 位置が違えば別として数えたときの辞書順 k 番目（0 から）  O(log N)（初回 O(N)）
        ensure_blocks();
        if (k < 0 || k >= bsum_.back()) str_internal::fail("SuffixArray::kth_substring_with_dup", "k は 0〜（部分文字列の数 - 1）");
        const int b = str_internal::upper_index(bsum_, k) - 1;
        const Block& B = blocks_[b];
        const int len = B.t_lo + (int)((k - bsum_[b]) / B.mult);
        return {sa[B.k], sa[B.k] + len};
    }

    // ---- 繰り返し・LCP の和 ----
    CommonSubstring longest_repeated() const {  // 2 か所以上に現れる最長の部分文字列（重なってよい。同じ長さなら辞書順最小）  O(N)
        CommonSubstring c;
        for (int k = 0; k + 1 < n_; ++k)
            if (lcp_adj[k] > c.len) c.len = lcp_adj[k], c.pos_a = std::min(sa[k], sa[k + 1]), c.pos_b = std::max(sa[k], sa[k + 1]);
        return c;
    }
    CommonSubstring longest_repeated_non_overlapping() const {  // 重ならずに 2 か所に現れる最長の部分文字列  O(N log N)
        int mx = 0;
        for (int x : lcp_adj) mx = std::max(mx, x);
        auto test = [&](int L, CommonSubstring* out) {
            for (int k = 0; k < n_;) {  // lcp_adj ≥ L でつながる SA の区間ごとに、最小と最大の位置
                int j = k, mn = sa[k], mxp = sa[k];
                while (j + 1 < n_ && lcp_adj[j] >= L) ++j, mn = std::min(mn, sa[j]), mxp = std::max(mxp, sa[j]);
                if (j > k && mxp - mn >= L) {
                    if (out) *out = {L, mn, mxp};
                    return true;
                }
                k = j + 1;
            }
            return false;
        };
        int lo = 0, hi = mx;
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (test(mid, nullptr)) lo = mid;
            else hi = mid - 1;
        }
        CommonSubstring c;
        if (lo > 0) test(lo, &c);
        return c;
    }
    std::vector<long long> sum_lcp_with_all() const {  // r[i] = Σ_j lcp(i, j)（j はすべての位置。j = i も含む）  O(N)
        std::vector<long long> left(n_, 0), right(n_, 0), r(n_, 0);
        sweep(left, false), sweep(right, true);
        for (int k = 0; k < n_; ++k) r[sa[k]] = left[k] + right[k] + length_at(sa[k]);
        return r;
    }
    long long sum_lcp_all_pairs() const {  // Σ_{i < j} lcp(i, j)（すべての位置の組）  O(N)
        std::vector<long long> left(n_, 0);
        sweep(left, false);
        long long t = 0;
        for (long long x : left) t += x;
        return t;
    }

   private:
    struct Block {
        int k, t_lo, t_hi;
        long long mult;
    };
    int lim(int p) const { return multi_ ? lim_[p] : n_; }
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    void check_index(const char* func, int i) const {
        if (i < 0 || i >= n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N-1（N = " + std::to_string(n_) + "）の外");
    }
    long long empty_count() const { return multi_ ? n_ : (long long)n_ + 1; }  // 空の列を置ける位置の数
    void check_sub(const char* func, int l, int r) const {
        if (l < 0 || l > n_ || r < l || (l < n_ && r > lim(l)) || (l == n_ && r != l))
            str_internal::fail(func, "部分文字列 [" + std::to_string(l) + ", " + std::to_string(r) + ") が範囲外（of_strings では 1 つの文字列の中で）");
    }
    void build(const std::vector<int>& v, int upper) {
        sa = str_internal::sa_is(v, upper);
        rank.assign(n_, 0);
        for (int k = 0; k < n_; ++k) rank[sa[k]] = k;
        lcp_adj = str_internal::lcp_kasai(v, sa);
    }
    // 接尾辞 p と pat：p が pat より小さい（pat で始まらず小さい）なら -1、pat で始まるなら 0、大きいなら 1
    int cmp_pat(int p, const std::vector<long long>& pat) const {
        const int e = lim(p);
        for (std::size_t t = 0; t < pat.size(); ++t) {
            if (p + (int)t >= e) return -1;
            const long long x = val_[p + t];
            if (x != pat[t]) return x < pat[t] ? -1 : 1;
        }
        return 0;
    }
    void ensure_table() const {
        if (built_ || n_ <= 1) return;
        built_ = true;
        const int m = n_ - 1;
        int lg = 1;
        while ((1 << lg) <= m) ++lg;
        tab_.assign(lg, std::vector<int>());
        tab_[0] = lcp_adj;
        for (int j = 1; j < lg; ++j) {
            const int len = m - (1 << j) + 1;
            if (len <= 0) break;
            tab_[j].resize(len);
            for (int i = 0; i < len; ++i) tab_[j][i] = std::min(tab_[j - 1][i], tab_[j - 1][i + (1 << (j - 1))]);
        }
    }
    int range_min(int a, int b) const {  // min lcp_adj[a, b)（a < b）
        ensure_table();
        const int j = 31 - __builtin_clz((unsigned)(b - a));
        return std::min(tab_[j][a], tab_[j][b - (1 << j)]);
    }
    void ensure_newcnt() const {
        if (!newcnt_.empty()) return;
        newcnt_.assign(n_ + 1, 0);
        for (int k = 0; k < n_; ++k) newcnt_[k + 1] = newcnt_[k] + std::max(0, length_at(sa[k]) - (k ? lcp_adj[k - 1] : 0));
    }
    void ensure_blocks() const {
        if (!bsum_.empty()) return;
        // k を大きい方から。st は lcp_adj[k..] の「その位置までの最小」が下がる位置（値は狭義に減る）
        std::vector<std::pair<int, int>> st;  // (値, 位置)
        std::vector<Block> tmp;
        for (int k = n_ - 1; k >= 0; --k) {
            if (k <= n_ - 2) {
                while (!st.empty() && st.back().first >= lcp_adj[k]) st.pop_back();
                st.push_back({lcp_adj[k], k});
            }
            const int L = length_at(sa[k]), lp = k ? lcp_adj[k - 1] : 0;
            // t ∈ (v_0, L]：自分だけ（mult 1）。t ∈ (v_{i+1}, v_i]：mult = j_{i+1} - k + 1（なければ n - k）
            const int top = (int)st.size() - 1;
            if (k <= n_ - 2) {
                const int v0 = st[top].first;
                if (L > std::max(lp, v0)) tmp.push_back({k, std::max(lp, v0) + 1, L, 1});
                for (int i = top; i >= 0 && st[i].first > lp; --i) {
                    const int vi = st[i].first, vnext = i > 0 ? st[i - 1].first : 0;
                    const long long mult = (i > 0 ? st[i - 1].second : n_ - 1) - k + 1;
                    const int lo = std::max(lp, vnext) + 1;
                    if (lo <= vi) tmp.push_back({k, lo, vi, mult});
                }
            } else if (L > lp) {
                tmp.push_back({k, lp + 1, L, 1});
            }
        }
        // tmp は k の降順、同じ k の中では t の降順 → 反転すると k の昇順・t の昇順
        blocks_.assign(tmp.rbegin(), tmp.rend());
        bsum_.assign(blocks_.size() + 1, 0);
        for (std::size_t b = 0; b < blocks_.size(); ++b) bsum_[b + 1] = bsum_[b] + blocks_[b].mult * (blocks_[b].t_hi - blocks_[b].t_lo + 1);
    }
    void sweep(std::vector<long long>& out, bool from_right) const {  // out[k] = Σ_{k' < k}（右からなら k' > k）min lcp
        std::vector<std::pair<int, long long>> st;  // (値, 個数)
        long long cur = 0;
        for (int step = 1; step < n_; ++step) {
            const int k = from_right ? n_ - 1 - step : step;
            const int v = from_right ? lcp_adj[k] : lcp_adj[k - 1];
            long long c = 1;
            while (!st.empty() && st.back().first >= v) cur -= (long long)st.back().first * st.back().second, c += st.back().second, st.pop_back();
            st.push_back({v, c}), cur += (long long)v * c;
            out[k] = cur;
        }
    }
    int n_ = 0;
    bool multi_ = false;
    std::vector<long long> val_;
    std::vector<int> lim_, owner_, start_;
    mutable bool built_ = false;
    mutable std::vector<std::vector<int>> tab_;
    mutable std::vector<long long> newcnt_, bsum_;
    mutable std::vector<Block> blocks_;
};

// ---- 2 つ・K 個の文字列の最長共通部分文字列（連続する部分。飛び飛びの「部分列」は 06_lcs_edit） ----
template <class A, class B>
CommonSubstring longest_common_substring(const A& a, const B& b) {  // 長さと、a での位置 pos_a・b での位置 pos_b（なければ len = 0, -1）  O(N + M)
    using V = std::vector<long long>;
    const SuffixArray S = SuffixArray::of_strings(std::vector<V>{str_internal::codes_of(a), str_internal::codes_of(b)});
    CommonSubstring c;
    for (int k = 0; k + 1 < S.size(); ++k) {
        const int p = S.sa[k], q = S.sa[k + 1];
        if (S.owner(p) < 0 || S.owner(q) < 0 || S.owner(p) == S.owner(q)) continue;
        if (S.lcp_adj[k] > c.len) {
            c.len = S.lcp_adj[k];
            c.pos_a = S.offset(S.owner(p) == 0 ? p : q), c.pos_b = S.offset(S.owner(p) == 0 ? q : p);
        }
    }
    return c;
}
template <class A, class B>
std::vector<StrSpan> longest_common_substrings_all(const A& a, const B& b) {  // 最長共通部分文字列をすべて（異なるもの・辞書順。a での位置）  O(N + M)
    using V = std::vector<long long>;
    const SuffixArray S = SuffixArray::of_strings(std::vector<V>{str_internal::codes_of(a), str_internal::codes_of(b)});
    int best = 0;
    for (int k = 0; k + 1 < S.size(); ++k) {
        const int p = S.sa[k], q = S.sa[k + 1];
        if (S.owner(p) >= 0 && S.owner(q) >= 0 && S.owner(p) != S.owner(q)) best = std::max(best, S.lcp_adj[k]);
    }
    std::vector<StrSpan> r;
    if (best == 0) return r;
    int last = -1;  // 直前に入れたものの位置（全体の位置）
    for (int k = 0; k + 1 < S.size(); ++k) {
        const int p = S.sa[k], q = S.sa[k + 1];
        if (S.owner(p) < 0 || S.owner(q) < 0 || S.owner(p) == S.owner(q) || S.lcp_adj[k] != best) continue;
        const int pa = S.owner(p) == 0 ? p : q;
        if (last >= 0 && S.lcp(last, pa) >= best) continue;
        last = pa;
        r.push_back({S.offset(pa), S.offset(pa) + best});
    }
    return r;
}
template <class Seq>
StrSpan longest_common_substring_k(const std::vector<Seq>& vs) {  // すべての文字列に現れる最長の部分文字列（vs[0] での位置。なければ {0, 0}）  O(全体の長さ)
    const int K = (int)vs.size();
    if (K == 0) return {0, 0};
    if (K == 1) return {0, (int)str_internal::seq(vs[0]).size()};
    const SuffixArray S = SuffixArray::of_strings(vs);
    const int n = S.size();
    std::vector<int> cnt(K, 0);
    int kinds = 0, best = 0, best_pos = -1, last0 = -1;  // last0 = 窓の中で最後の「0 番目の文字列」の SA の位置
    std::vector<int> dq(n + 1);  // lcp_adj の位置の単調な列（最小を先頭に）
    int qh = 0, qt = 0;
    for (int l = 0, r = 0; r < n; ++r) {
        const int o = S.owner(S.sa[r]);
        if (o >= 0 && cnt[o]++ == 0) ++kinds;
        if (o == 0) last0 = r;
        if (r > 0) {
            while (qt > qh && S.lcp_adj[dq[qt - 1]] >= S.lcp_adj[r - 1]) --qt;
            dq[qt++] = r - 1;
        }
        while (kinds == K) {
            if (r > l) {
                while (qh < qt && dq[qh] < l) ++qh;
                const int v = S.lcp_adj[dq[qh]];
                if (v > best) best = v, best_pos = S.offset(S.sa[last0]);
            }
            const int o2 = S.owner(S.sa[l]);
            if (o2 >= 0 && --cnt[o2] == 0) --kinds;
            ++l;
        }
    }
    if (best == 0) return {0, 0};
    return {best_pos, best_pos + best};
}
#endif  // GG_SUFFIX_ARRAY
