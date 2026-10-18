#ifndef GG_PATTERN_MATCH
#define GG_PATTERN_MATCH
// ---- パターン検索（文字列第二世代）  Z・KMP・出現・周期・境界・回転・重なり・置き換えの繰り返し・区間の出現・KMP オートマトン ----
// text・pat は string でも vector<整数> でもよい（混ぜてもよい。値で比べる）。文字列リテラルも渡せる
namespace str_internal {
// 列（string・vector<整数>）をそのまま読む（コピーしない）。要素は code_of で比べる
template <class S>
std::vector<int> z_codes(const S& s) {
    const int n = (int)s.size();
    std::vector<int> z(n);
    if (n == 0) return z;
    z[0] = n;
    for (int i = 1, l = 0, r = 0; i < n; ++i) {
        int k = i < r ? std::min(r - i, z[i - l]) : 0;
        while (i + k < n && s[k] == s[i + k]) ++k;
        z[i] = k;
        if (i + k > r) l = i, r = i + k;
    }
    return z;
}
template <class S>
std::vector<int> pi_codes(const S& s) {
    const int n = (int)s.size();
    std::vector<int> pi(n);
    for (int i = 1; i < n; ++i) {
        int k = pi[i - 1];
        while (k > 0 && !(s[i] == s[k])) k = pi[k - 1];
        pi[i] = k + (s[i] == s[k]);
    }
    return pi;
}
// text[from..) の中の pat の出現開始位置（text での位置。重なりも）。pat が空なら from〜|text|
template <class T, class P>
std::vector<int> kmp_find(const T& t, const P& p, int from = 0) {
    const int n = (int)t.size(), m = (int)p.size();
    std::vector<int> r;
    if (m == 0) {
        for (int i = from; i <= n; ++i) r.push_back(i);
        return r;
    }
    const std::vector<int> pi = pi_codes(p);
    for (int i = from, k = 0; i < n; ++i) {
        const long long c = code_of(t[i]);
        while (k > 0 && (k == m || c != code_of(p[k]))) k = pi[k - 1];
        if (c == code_of(p[k])) ++k;
        if (k == m) r.push_back(i - m + 1);
    }
    return r;
}
}  // namespace str_internal

template <class Seq>
std::vector<int> z_array(const Seq& s) { return str_internal::z_codes(str_internal::seq(s)); }  // z[i] = s と s[i..) の最長共通接頭辞（z[0] = |s|）  O(N)
template <class Seq>
std::vector<int> prefix_function(const Seq& s) { return str_internal::pi_codes(str_internal::seq(s)); }  // π[i] = s[0, i] の最長の真の境界の長さ  O(N)
template <class T, class P>
std::vector<int> match_lengths(const T& text, const P& pat) {  // 各 i で text[i..) と pat の最長共通接頭辞の長さ（長さ |text|）  O(N + M)
    const auto& t = str_internal::seq(text);
    const auto& p = str_internal::seq(pat);
    const int n = (int)t.size(), m = (int)p.size();
    const std::vector<int> zp = str_internal::z_codes(p);
    std::vector<int> res(n);
    for (int i = 0, l = 0, r = 0; i < n; ++i) {  // [l, r) は text[l, r) == pat[0, r - l) となる右端が最大の区間
        int k = 0;
        if (i < r) k = std::min(r - i, zp[i - l]);
        if (i >= r || zp[i - l] >= r - i) {
            while (k < m && i + k < n && str_internal::code_of(t[i + k]) == str_internal::code_of(p[k])) ++k;
        }
        res[i] = k;
        if (i + k > r) l = i, r = i + k;
    }
    return res;
}
template <class T, class P>
std::vector<int> find_all(const T& text, const P& pat) {  // pat が現れる開始位置（重なりも。昇順。pat が空なら 0〜|text|）  O(N + M)
    return str_internal::kmp_find(str_internal::seq(text), str_internal::seq(pat));
}
template <class T, class P>
long long count_occurrences(const T& text, const P& pat) { return (long long)find_all(text, pat).size(); }  // 重なりも数える
template <class T, class P>
int find_first(const T& text, const P& pat, int from = 0) {  // from 以上で最初に現れる位置（なければ -1）  O(N + M)
    const auto& t = str_internal::seq(text);
    if (from < 0 || from > (int)t.size()) str_internal::fail("find_first", "from は 0〜|text| で指定する");
    const auto r = str_internal::kmp_find(t, str_internal::seq(pat), from);
    return r.empty() ? -1 : r[0];
}
template <class T, class P>
int find_last(const T& text, const P& pat) {  // 最後に現れる位置（なければ -1）
    const auto r = find_all(text, pat);
    return r.empty() ? -1 : r.back();
}
template <class T, class P>
long long count_non_overlapping(const T& text, const P& pat) {  // 左から重ならないように取った数（pat は空でない）
    const auto& p = str_internal::seq(pat);
    if (p.empty()) str_internal::fail("count_non_overlapping", "pat が空");
    long long c = 0, last = -1;
    for (int i : str_internal::kmp_find(str_internal::seq(text), p))
        if (i >= last) ++c, last = i + (long long)p.size();
    return c;
}
template <class T, class F, class U>
str_internal::seq_t<T> replace_all(const T& text0, const F& from, const U& to0) {  // from を左から重ならないように全部 to に（from は空でない）  O(N + M)
    const auto& text = str_internal::seq(text0);
    const auto& to = str_internal::seq(to0);
    const auto& f = str_internal::seq(from);
    if (f.empty()) str_internal::fail("replace_all", "from が空");
    str_internal::seq_t<T> r;
    std::size_t pos = 0;
    for (int i : str_internal::kmp_find(text, f)) {
        if ((std::size_t)i < pos) continue;
        r.insert(r.end(), text.begin() + pos, text.begin() + i);
        for (std::size_t j = 0; j < to.size(); ++j) r.push_back((typename str_internal::seq_t<T>::value_type)str_internal::code_of(to[j]));
        pos = (std::size_t)i + f.size();
    }
    r.insert(r.end(), text.begin() + pos, text.end());
    return r;
}

// ---- 周期・境界 ----
template <class Seq>
std::vector<int> all_borders(const Seq& s) {  // 0 < k < |s| で s[0, k) == s[|s|-k, |s|) となる k（大きい順）  O(N)
    const auto pi = prefix_function(s);
    std::vector<int> r;
    if (pi.empty()) return r;
    for (int k = pi.back(); k > 0; k = pi[k - 1]) r.push_back(k);
    return r;
}
template <class Seq>
int min_period(const Seq& s) {  // 最小の周期 p（s[i] == s[i + p]。割り切れなくてよい。空なら 0）  O(N)
    const auto pi = prefix_function(s);
    return pi.empty() ? 0 : (int)pi.size() - pi.back();
}
template <class Seq>
std::vector<int> all_periods(const Seq& s) {  // 1〜|s| のすべての周期（小さい順。|s| 自身も含む）  O(N)
    const int n = (int)str_internal::seq(s).size();
    std::vector<int> r;
    for (int k : all_borders(s)) r.push_back(n - k);
    if (n > 0) r.push_back(n);
    return r;
}
template <class Seq>
int repetition_root(const Seq& s) {  // s = t^k（t を k 回）となる最短の t の長さ（空なら 0）  O(N)
    const int n = (int)str_internal::seq(s).size(), p = min_period(s);
    return (n > 0 && n % p == 0) ? p : n;
}
template <class Seq>
std::vector<int> prefix_min_periods(const Seq& s) {  // r[i] = s[0, i) の最小の周期（長さ |s| + 1。r[0] = 0）  O(N)
    const auto pi = prefix_function(s);
    std::vector<int> r(pi.size() + 1, 0);
    for (std::size_t i = 0; i < pi.size(); ++i) r[i + 1] = (int)(i + 1) - pi[i];
    return r;
}

// ---- 回転・重なり ----
template <class A, class B>
std::vector<int> rotation_offsets(const A& a, const B& b) {  // a を左に k 回まわす（a[k..] + a[..k]）と b になる k（0 ≤ k < |a|。昇順。a, b がともに空なら {0}）  O(N)
    const auto x = str_internal::codes_of(a), y = str_internal::codes_of(b);
    std::vector<int> r;
    if (x.size() != y.size()) return r;
    if (x.empty()) return {0};
    std::vector<long long> xx(x);
    xx.insert(xx.end(), x.begin(), x.end() - 1);
    return str_internal::kmp_find(xx, y);
}
template <class A, class B>
bool is_rotation(const A& a, const B& b) { return !rotation_offsets(a, b).empty(); }  // 回転して一致するか
template <class T, class P>
std::vector<int> find_all_cyclic(const T& text, const P& pat) {  // text を無限に繰り返した列で pat が始まる位置 i（0 ≤ i < |text|。|pat| > |text| でもよい）  O(N + M)
    const auto t = str_internal::codes_of(text), p = str_internal::codes_of(pat);
    const std::size_t n = t.size(), m = p.size();
    if (n == 0) return {};
    std::vector<long long> tt;
    tt.reserve(n + m);
    while (tt.size() < n + (m > 0 ? m - 1 : 0)) tt.insert(tt.end(), t.begin(), t.end());
    std::vector<int> r;
    for (int i : str_internal::kmp_find(tt, p))
        if ((std::size_t)i < n) r.push_back(i);
    return r;
}
template <class A, class B>
int max_overlap(const A& a, const B& b) {  // a の接尾辞と b の接頭辞が一致する最大の長さ（a + b を短くつなげるとき）  O(N + M)
    const auto x = str_internal::codes_of(a), y = str_internal::codes_of(b);
    const int m = (int)y.size();
    if (m == 0) return 0;
    const std::vector<int> pi = str_internal::pi_codes(y);
    int k = 0;
    for (long long c : x) {
        if (k == m) k = pi[k - 1];
        while (k > 0 && c != y[k]) k = pi[k - 1];
        if (c == y[k]) ++k;
    }
    return k;
}

// ---- 置き換え・消すを、現れなくなるまで繰り返す（いちばん左の出現から） ----
//   スタックに KMP の状態を積み、置き換えた to は読み直す。O(N + 置き換えの回数 × (|from| + |to|))
//   |to| < |from| なら必ず終わる。置き換えが起きて to が from を含むと終わらないので止める（伸び続けるものは長さ 2×10^7 で止める）
//   |to| ≥ |from| で置き換えの回数が多いと遅い（ba → ab は Θ(N^2) 回）。to の値は text の要素の型に収まること
template <class T, class F, class U>
str_internal::seq_t<T> rewrite_until_stable(const T& text0, const F& from0, const U& to0) {  // from を to に（ABC394 C：WA → AC）
    using R = str_internal::seq_t<T>;
    using E = typename R::value_type;
    const auto& text = str_internal::seq(text0);
    const std::vector<long long> f = str_internal::codes_of(from0), to = str_internal::codes_of(to0);
    const int m = (int)f.size();
    if (m == 0) str_internal::fail("rewrite_until_stable", "from が空");
    const bool loops = !str_internal::kmp_find(to, f).empty();  // 1 回でも置き換えると終わらない
    const std::vector<int> pi = str_internal::pi_codes(f);
    std::vector<E> out;
    std::vector<int> st;          // st[k] = out[0, k] を読んだあとの KMP の状態
    std::vector<long long> back;  // 読み直す値（後ろから取り出す）
    out.reserve(text.size()), st.reserve(text.size());
    std::size_t idx = 0;
    while (idx < text.size() || !back.empty()) {
        E e;
        long long c;
        if (!back.empty()) c = back.back(), back.pop_back(), e = (E)c;
        else e = text[idx++], c = str_internal::code_of(e);
        int q = st.empty() ? 0 : st.back();
        while (q > 0 && f[q] != c) q = pi[q - 1];
        if (f[q] == c) ++q;
        out.push_back(e), st.push_back(q);
        if (q == m) {  // いちばん左の出現（スタックの中には出現がない）
            if (loops) str_internal::fail("rewrite_until_stable", "to が from を含む（いつまでも終わらない）");
            out.resize(out.size() - m), st.resize(st.size() - m);
            back.insert(back.end(), to.rbegin(), to.rend());
            if (out.size() + back.size() + (text.size() - idx) > 20000000) str_internal::fail("rewrite_until_stable", "長さが 2×10^7 を超えた（終わらない置き換え？）");
        }
    }
    return R(out.begin(), out.end());
}
template <class T, class P>
str_internal::seq_t<T> erase_until_stable(const T& text, const P& pat) {  // pat を消すことを現れなくなるまで（ABC328 D：ABC を消す）  O(N + M)
    return rewrite_until_stable(text, pat, std::vector<long long>());
}

// ---- 出現の一覧から、区間 [l, r) に収まる出現の数・i 以降で最初の出現を O(log N) で（ABC477 C） ----
struct Occurrences {
    Occurrences() = default;
    template <class T, class P>
    Occurrences(const T& text, const P& pat) : pos_(find_all(text, pat)), n_((int)str_internal::seq(text).size()), m_((int)str_internal::seq(pat).size()) {}
    const std::vector<int>& positions() const { return pos_; }  // 開始位置（昇順）
    int size() const { return (int)pos_.size(); }               // 出現の数
    long long count(int l, int r) const {                        // [l, r) に収まる（開始が l 以上・終わりが r 以下）出現の数  O(log N)
        str_internal::check_range("Occurrences::count", l, r, n_);
        if (r - l < m_) return 0;
        return idx(r - m_ + 1) - idx(l);
    }
    bool exists(int l, int r) const { return count(l, r) > 0; }  // [l, r) に収まる出現があるか
    int next(int i) const {                                       // i 以上で最初の出現の開始位置（なければ -1。0 ≤ i ≤ N）
        check_pos("Occurrences::next", i);
        const int k = idx(i);
        return k < size() ? pos_[k] : -1;
    }
    int prev(int i) const {  // i より前（開始位置 < i）で最後の出現の開始位置（なければ -1。0 ≤ i ≤ N）
        check_pos("Occurrences::prev", i);
        const int k = idx(i);
        return k > 0 ? pos_[k - 1] : -1;
    }

   private:
    int idx(int x) const { return str_internal::lower_index(pos_.data(), 0, size(), x); }  // 開始位置が x 以上の最初の番号
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    std::vector<int> pos_;
    int n_ = 0, m_ = 0;
};

// ---- KMP オートマトン：状態 = pat の先頭何文字まで一致しているか（0〜|pat|）。|pat| に着いたら pat が現れた ----
struct KmpAutomaton {
    KmpAutomaton() = default;
    template <class P>
    explicit KmpAutomaton(const P& pat0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        const auto& pat = str_internal::seq(pat0);
        m_ = (int)pat.size();
        sg_ = A_.size();
        std::vector<int> p(m_);
        for (int i = 0; i < m_; ++i) p[i] = A_.id(pat[i]);
        nx_.assign((std::size_t)(m_ + 1) * sg_, 0);
        int fail = 0;  // 状態 q の「1 つ短い一致」
        for (int q = 0; q <= m_; ++q) {
            for (int c = 0; c < sg_; ++c) nx_[(std::size_t)q * sg_ + c] = (q < m_ && p[q] == c) ? q + 1 : (q == 0 ? 0 : nx_[(std::size_t)fail * sg_ + c]);
            if (q > 0 && q < m_) fail = nx_[(std::size_t)fail * sg_ + p[q]];
        }
    }
    int size() const { return m_ + 1; }  // 状態の数（0〜|pat|）
    template <class C>
    int next(int state, C c) const { return next_id(state, A_.id(c)); }  // 状態 state で文字 c を読んだあとの状態  O(1)
    int next_id(int state, int id) const {
        if (state < 0 || state > m_ || id < 0 || id >= sg_) str_internal::fail("KmpAutomaton::next", "状態か文字の番号が範囲外");
        return nx_[(std::size_t)state * sg_ + id];
    }
    bool accepting(int state) const { return state == m_; }  // pat がちょうど現れたところか
    const Alphabet& alphabet() const { return A_; }

   private:
    Alphabet A_;
    int m_ = 0, sg_ = 0;
    std::vector<int> nx_;
};
#endif  // GG_PATTERN_MATCH
