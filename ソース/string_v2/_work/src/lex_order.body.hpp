#ifndef GG_LEX_ORDER
#define GG_LEX_ORDER
// ---- 辞書順（文字列第二世代）  つなげて最小・独自の文字の順・最小の回転・列名（a, …, z, aa）・辞書の番号・同じ長さでの番号・次の文字列 ----
// 辞書順は「先頭から比べ、最初に違う文字で決まる。片方がもう片方の接頭辞なら短い方が小さい」（std::string と同じ）
namespace str_internal {
__extension__ typedef __int128 lex_i128;
// a + b と b + a の辞書順（-1 / 0 / 1）を、つなげた列を作らずに
template <class S>
int concat_cmp(const S& a, const S& b) {
    const std::size_t n = a.size(), m = b.size(), t = n + m;
    for (std::size_t i = 0; i < t; ++i) {
        const long long x = code_of(i < n ? a[i] : b[i - n]), y = code_of(i < m ? b[i] : a[i - m]);
        if (x != y) return x < y ? -1 : 1;
    }
    return 0;
}
inline lex_i128 lex_big() { return (lex_i128)1000000000000000000LL * 10; }  // 10^19（long long の最大より大きい）
inline lex_i128 lex_cap(lex_i128 x) { return x > lex_big() ? lex_big() : x; }
// 長さ 0〜t の文字列の数 1 + σ + … + σ^t（10^19 で止める）
inline lex_i128 lex_count_upto(long long sigma, long long t) {
    if (t < 0) return 0;
    if (sigma == 1) return lex_cap((lex_i128)t + 1);
    lex_i128 sum = 0, p = 1;
    for (long long u = 0; u <= t; ++u) {
        sum = lex_cap(sum + p);
        p = lex_cap(p * sigma);
        if (sum >= lex_big()) break;
    }
    return sum;
}
}  // namespace str_internal

// ---- つなげて辞書順最小 ----
template <class S>
std::vector<int> concat_order(const std::vector<S>& vs) {  // 全部をつなげて辞書順最小になる並べ方（添字。a + b < b + a の順。同じなら元の順。空の列は最後）  O(N log N × 最長の長さ)
    std::vector<int> idx, empty;
    for (int i = 0; i < (int)vs.size(); ++i) (str_internal::seq(vs[i]).size() ? idx : empty).push_back(i);  // 空の列は順序を壊すので外す
    std::stable_sort(idx.begin(), idx.end(), [&](int i, int j) { return str_internal::concat_cmp(vs[i], vs[j]) < 0; });
    idx.insert(idx.end(), empty.begin(), empty.end());
    return idx;
}
template <class S>
S min_concat(const std::vector<S>& vs) {  // 全部をつなげた辞書順最小の列
    S r;
    for (int i : concat_order(vs)) r.insert(r.end(), vs[i].begin(), vs[i].end());
    return r;
}

// ---- 独自の文字の順での辞書順（order に並べた順。order にない文字は止める） ----
struct LexLess {
    explicit LexLess(const Alphabet& A) : A_(A) {}
    template <class X, class Y>
    int compare(const X& a0, const Y& b0) const {  // -1 / 0 / 1
        const auto& a = str_internal::seq(a0);
        const auto& b = str_internal::seq(b0);
        const std::size_t k = std::min(a.size(), b.size());
        for (std::size_t i = 0; i < k; ++i) {
            const int x = A_.id(a[i]), y = A_.id(b[i]);
            if (x != y) return x < y ? -1 : 1;
        }
        return a.size() == b.size() ? 0 : (a.size() < b.size() ? -1 : 1);
    }
    template <class X, class Y>
    bool operator()(const X& a, const Y& b) const { return compare(a, b) < 0; }  // sort の比較関数に

   private:
    Alphabet A_;
};
inline LexLess lex_less_by(const std::string& order) { return LexLess(Alphabet::of(order)); }  // sort(v.begin(), v.end(), lex_less_by("bac…"))
inline LexLess lex_less_by(const Alphabet& A) { return LexLess(A); }

// ---- 回転 ----
namespace str_internal {
template <class Seq>
int extreme_rotation(const Seq& s0, bool want_min) {
    const auto s = codes_of(s0);
    const int n = (int)s.size();
    if (n == 0) return 0;
    int i = 0, j = 1, k = 0;
    while (i < n && j < n && k < n) {
        const int ia = i + k < n ? i + k : i + k - n, ib = j + k < n ? j + k : j + k - n;
        const long long a = s[ia], b = s[ib];
        if (a == b) {
            ++k;
            continue;
        }
        if ((a > b) == want_min) i += k + 1;
        else j += k + 1;
        if (i == j) ++j;
        k = 0;
    }
    int r = std::min(i, j);
    // 周期があれば最小の位置に（s = t^m なら r mod |t|）
    std::vector<int> pi(n);
    for (int x = 1; x < n; ++x) {
        int q = pi[x - 1];
        while (q > 0 && s[x] != s[q]) q = pi[q - 1];
        pi[x] = q + (s[x] == s[q]);
    }
    const int p = n - pi[n - 1];
    if (n % p == 0) r %= p;
    return r;
}
}  // namespace str_internal
template <class Seq>
int min_rotation(const Seq& s) { return str_internal::extreme_rotation(s, true); }  // 辞書順最小の回転 s[k..] + s[..k] の k（同じなら最小の k）  O(N)
template <class Seq>
int max_rotation(const Seq& s) { return str_internal::extreme_rotation(s, false); }  // 辞書順最大の回転の k  O(N)

// ---- 列名（a, b, …, z, aa, ab, …）：全単射の記数法 ----
inline std::string bijective_to_string(long long n, const Alphabet& A = Alphabet::lower()) {  // 1 → "a"、26 → "z"、27 → "aa"（0 なら空）。ABC171 C
    if (n < 0) str_internal::fail("bijective_to_string", "n ≥ 0 で使う");
    const long long k = A.size();
    std::string r;
    while (n > 0) {
        --n;
        r += A.ch((int)(n % k));
        n /= k;
    }
    std::reverse(r.begin(), r.end());
    return r;
}
template <class Seq>
long long bijective_from_string(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // "a" → 1、"z" → 26、"aa" → 27（ABC285 C）。long long に収まらなければ止める
    const auto& s = str_internal::seq(s0);
    const long long k = A.size();
    str_internal::lex_i128 v = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        v = v * k + (A.id(s[i]) + 1);
        if (v > (str_internal::lex_i128)std::numeric_limits<long long>::max()) str_internal::fail("bijective_from_string", "long long に収まらない");
    }
    return (long long)v;
}

// ---- 辞書の番号：長さ 1〜max_len のすべての文字列（A の文字）を辞書順に並べたときの番号（0 から） ----
template <class Seq>
long long dict_rank(const Seq& s0, long long max_len, const Alphabet& A = Alphabet::lower()) {  // s より小さい文字列の数  O(|s| + max_len)
    const auto& s = str_internal::seq(s0);
    const long long m = (long long)s.size(), sg = A.size();
    if (m < 1 || m > max_len) str_internal::fail("dict_rank", "1 ≤ |s| ≤ max_len で使う");
    str_internal::lex_i128 r = m - 1;  // s の真の接頭辞（空でない）
    for (long long d = 0; d < m; ++d) {
        const int c = A.id(s[d]);
        if (c) r += (str_internal::lex_i128)c * str_internal::lex_count_upto(sg, max_len - d - 1);
        if (r > (str_internal::lex_i128)std::numeric_limits<long long>::max()) str_internal::fail("dict_rank", "番号が long long に収まらない");
    }
    return (long long)r;
}
inline std::string dict_unrank(long long k, long long max_len, const Alphabet& A = Alphabet::lower()) {  // 番号 k（0 から）の文字列  O(max_len)
    if (k < 0 || max_len < 1) str_internal::fail("dict_unrank", "k ≥ 0, max_len ≥ 1 で使う");
    const long long sg = A.size();
    if (str_internal::lex_count_upto(sg, max_len) - 1 <= k) str_internal::fail("dict_unrank", "k が文字列の数以上");
    std::string r;
    str_internal::lex_i128 x = k;
    for (long long d = 0; d < max_len; ++d) {
        const str_internal::lex_i128 block = str_internal::lex_count_upto(sg, max_len - d - 1);  // r + c で始まる文字列の数
        const long long c = (long long)(x / block);
        r += A.ch((int)c);
        x -= (str_internal::lex_i128)c * block;
        if (x == 0) break;
        --x;  // r + c そのものを飛ばす
    }
    return r;
}
template <class Seq>
str_internal::seq_t<Seq> next_string(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // 同じ長さで辞書順の次（最後なら空）  O(N)
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    for (int i = (int)s.size() - 1; i >= 0; --i) {
        const int c = A.id(s[i]);
        if (c + 1 < A.size()) {
            s[i] = (typename str_internal::seq_t<Seq>::value_type)A.value(c + 1);
            for (int j = i + 1; j < (int)s.size(); ++j) s[j] = (typename str_internal::seq_t<Seq>::value_type)A.value(0);
            return s;
        }
    }
    return {};
}
template <class T, class Seq>
T lex_rank_same_length(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // A の文字からなる長さ |s| の文字列で s より小さいものの数（= s を σ 進数で読んだ値）  O(N)
    const auto& s = str_internal::seq(s0);
    T r = T(0);
    for (std::size_t i = 0; i < s.size(); ++i)
        r = str_internal::checked_add(str_internal::checked_mul(r, T(A.size()), "lex_rank_same_length"), T(A.id(s[i])), "lex_rank_same_length");
    return r;
}
#endif  // GG_LEX_ORDER
