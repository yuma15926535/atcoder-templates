#ifndef GG_PALINDROME
#define GG_PALINDROME
// ---- 回文（文字列第二世代）  Manacher・区間の回文判定 O(1)・最長の回文・回文の個数・回文にする・回文部分列・回文分割・s 以下の回文の数 ----
// 中心の番号 c（0 ≤ c < 2N - 1）：c = 2i は s[i] が中心（奇数長）、c = 2i - 1 は s[i-1] と s[i] の間が中心（偶数長）
namespace str_internal {
inline std::vector<int> manacher_codes(const std::vector<long long>& s) {  // L[c] = 中心 c の最長の回文の長さ
    const int n = (int)s.size();
    if (n == 0) return {};
    std::vector<int> d1(n), d2(n), L(2 * n - 1);
    for (int i = 0, l = 0, r = -1; i < n; ++i) {
        int k = (i > r) ? 1 : std::min(d1[l + r - i], r - i + 1);
        while (i - k >= 0 && i + k < n && s[i - k] == s[i + k]) ++k;
        d1[i] = k--;
        if (i + k > r) l = i - k, r = i + k;
    }
    for (int i = 0, l = 0, r = -1; i < n; ++i) {
        int k = (i > r) ? 0 : std::min(d2[l + r - i + 1], r - i + 1);
        while (i - k - 1 >= 0 && i + k < n && s[i - k - 1] == s[i + k]) ++k;
        d2[i] = k--;
        if (i + k > r) l = i - k - 1, r = i + k;
    }
    for (int i = 0; i < n; ++i) L[2 * i] = 2 * d1[i] - 1;
    for (int i = 1; i < n; ++i) L[2 * i - 1] = 2 * d2[i];
    return L;
}
}  // namespace str_internal

template <class Seq>
std::vector<int> center_palindrome_lengths(const Seq& s) {  // 中心 c（0〜2N-2）ごとの最長の回文の長さ（Manacher）  O(N)
    return str_internal::manacher_codes(str_internal::codes_of(s));
}

struct Palindromes {
    Palindromes() = default;
    template <class Seq>
    explicit Palindromes(const Seq& s) : n_((int)str_internal::seq(s).size()), L_(center_palindrome_lengths(s)) {}
    int size() const { return n_; }
    bool is_palindrome(int l, int r) const {  // s[l, r) が回文か（空・1 文字は回文）  O(1)
        str_internal::check_range("Palindromes::is_palindrome", l, r, n_);
        return r - l <= 1 || L_[l + r - 1] >= r - l;
    }
    int radius_odd(int i) const {  // s[i] が中心の最長の回文の半径 r（長さ 2r - 1）
        if (i < 0 || i >= n_) str_internal::fail("Palindromes::radius_odd", "i は 0〜N-1");
        return (L_[2 * i] + 1) / 2;
    }
    int radius_even(int i) const {  // s[i-1] と s[i] の間が中心の最長の回文の半径 r（長さ 2r。i = 0, N なら 0）
        if (i < 0 || i > n_) str_internal::fail("Palindromes::radius_even", "i は 0〜N");
        return (i == 0 || i == n_) ? 0 : L_[2 * i - 1] / 2;
    }
    int length_at_center(int c) const {  // 中心 c の最長の回文の長さ
        if (c < 0 || c >= 2 * n_ - 1) str_internal::fail("Palindromes::length_at_center", "c は 0〜2N-2");
        return L_[c];
    }
    StrSpan span_at_center(int c) const {  // 中心 c の最長の回文の位置
        const int len = length_at_center(c), l = (c + 1 - len) / 2;
        return {l, l + len};
    }
    StrSpan longest() const {  // 最長の回文の部分文字列（同じ長さなら左）  O(N)
        StrSpan best{0, 0};
        for (int c = 0; c < 2 * n_ - 1; ++c)
            if (L_[c] > best.len() || (L_[c] == best.len() && (c + 1 - L_[c]) / 2 < best.l)) best = span_at_center(c);
        return best;
    }
    long long count_all() const {  // 回文である部分文字列（空でない）の数。位置が違えば別  O(N)
        long long t = 0;
        for (int c = 0; c < 2 * n_ - 1; ++c) t += (L_[c] + (c % 2 == 0)) / 2;
        return t;
    }
    int longest_prefix_palindrome() const {  // s[0, r) が回文となる最大の r  O(N)
        for (int r = n_; r > 0; --r)
            if (is_palindrome(0, r)) return r;
        return 0;
    }
    int longest_suffix_palindrome() const {  // s[N - k, N) が回文となる最大の k  O(N)
        for (int l = 0; l < n_; ++l)
            if (is_palindrome(l, n_)) return n_ - l;
        return 0;
    }
    std::vector<int> longest_starting_at() const {  // r[i] = s[i, i + len) が回文となる最大の len  O(N)
        std::vector<int> best(n_, -1), r(n_);
        for (int c = 0; c < 2 * n_ - 1; ++c) {
            const int l = (c + 1 - L_[c]) / 2;
            best[l] = std::max(best[l], c);
        }
        for (int i = 0, mx = -1; i < n_; ++i) mx = std::max(mx, best[i]), r[i] = mx - 2 * i + 1;
        return r;
    }
    std::vector<int> longest_ending_at() const {  // r[j] = s[j - len + 1, j + 1) が回文となる最大の len（j が最後の文字）  O(N)
        std::vector<int> best(n_, 2 * n_), r(n_);
        for (int c = 0; c < 2 * n_ - 1; ++c) {
            const int e = (c + 1 - L_[c]) / 2 + L_[c] - 1;  // 最後の文字の位置
            best[e] = std::min(best[e], c);
        }
        for (int j = n_ - 1, mn = 2 * n_; j >= 0; --j) mn = std::min(mn, best[j]), r[j] = 2 * j - mn + 1;
        return r;
    }

   private:
    int n_ = 0;
    std::vector<int> L_;
};

// ---- 回文にする・回文の個数 ----
template <class Seq>
str_internal::seq_t<Seq> min_palindrome_with_prefix(const Seq& s0) {  // s で始まる最短の回文（後ろに足す。ABC398 F）  O(N)
    auto s = str_internal::seq_t<Seq>(str_internal::seq(s0));
    const int n = (int)s.size(), k = Palindromes(s).longest_suffix_palindrome();
    for (int i = n - k - 1; i >= 0; --i) s.push_back(s[i]);
    return s;
}
template <class Seq>
str_internal::seq_t<Seq> min_palindrome_with_suffix(const Seq& s0) {  // s で終わる最短の回文（前に足す）  O(N)
    const auto& s = str_internal::seq(s0);
    const int n = (int)s.size(), k = Palindromes(s).longest_prefix_palindrome();
    str_internal::seq_t<Seq> r;
    for (int i = n - 1; i >= k; --i) r.push_back(s[i]);
    r.insert(r.end(), s.begin(), s.end());
    return r;
}
template <class Seq>
int min_changes_to_palindrome(const Seq& s0) {  // 1 文字ずつ書き換えて回文にする最小の回数  O(N)
    const auto& s = str_internal::seq(s0);
    int c = 0;
    for (std::size_t i = 0, n = s.size(); i < n / 2; ++i) c += !(s[i] == s[n - 1 - i]);
    return c;
}
template <class Seq>
bool can_rearrange_to_palindrome(const Seq& s0) {  // 並べ替えて回文にできるか（奇数個の文字が 1 種類以下）  O(N log N)
    auto c = str_internal::codes_of(s0);
    std::sort(c.begin(), c.end());
    int odd = 0;
    for (std::size_t i = 0, j; i < c.size(); i = j) {
        for (j = i; j < c.size() && c[j] == c[i];) ++j;
        odd += (int)((j - i) & 1);
    }
    return odd <= 1;
}
template <class Seq>
int longest_palindromic_subsequence_length(const Seq& s0) {  // 最長回文部分列（飛び飛び）の長さ  O(N^2) 時間・O(N) メモリ
    const auto s = str_internal::codes_of(s0);
    const int n = (int)s.size();
    std::vector<int> nxt(n + 1, 0), cur(n + 1, 0);  // nxt[j] = LPS(s[i+1..j])
    for (int i = n - 1; i >= 0; --i) {
        cur[i] = 1;
        for (int j = i + 1; j < n; ++j) cur[j] = s[i] == s[j] ? (j - 1 >= i + 1 ? nxt[j - 1] : 0) + 2 : std::max(nxt[j], cur[j - 1]);
        std::swap(cur, nxt);
    }
    return n ? nxt[n - 1] : 0;
}
template <class Seq>
str_internal::seq_t<Seq> longest_palindromic_subsequence(const Seq& s0) {  // 最長回文部分列（飛び飛び）を 1 つ  O(N^2) 時間・2N^2 バイト（N ≤ 10^4）
    const auto& s = str_internal::seq(s0);
    const auto c = str_internal::codes_of(s0);
    const int n = (int)c.size();
    if (n > 10000) str_internal::fail("longest_palindromic_subsequence", "N ≤ 10^4 で使う（長さだけなら longest_palindromic_subsequence_length）");
    std::vector<unsigned short> dp((std::size_t)n * n, 0);  // dp[i * n + j] = LPS(s[i..j])（i ≤ j）
    auto at = [&](int i, int j) -> unsigned short& { return dp[(std::size_t)i * n + j]; };
    for (int i = n - 1; i >= 0; --i) {
        at(i, i) = 1;
        for (int j = i + 1; j < n; ++j) at(i, j) = c[i] == c[j] ? (unsigned short)((j - 1 >= i + 1 ? at(i + 1, j - 1) : 0) + 2) : std::max(at(i + 1, j), at(i, j - 1));
    }
    str_internal::seq_t<Seq> left, right;
    int i = 0, j = n - 1;
    while (i <= j) {
        if (i == j) {
            left.push_back(s[i]);
            break;
        }
        if (c[i] == c[j]) left.push_back(s[i]), right.push_back(s[j]), ++i, --j;
        else if (at(i + 1, j) >= at(i, j - 1)) ++i;
        else --j;
    }
    left.insert(left.end(), right.rbegin(), right.rend());
    return left;
}
template <class Seq>
std::vector<StrSpan> palindrome_partition_min(const Seq& s) {  // 回文に分ける最小の個数で分けたもの（.size() が個数）  O(N^2)
    const Palindromes P(s);
    const int n = P.size();
    std::vector<int> dp(n + 1, 1 << 30), from(n + 1, -1);
    dp[0] = 0;
    for (int r = 1; r <= n; ++r)
        for (int l = 0; l < r; ++l)
            if (dp[l] + 1 < dp[r] && P.is_palindrome(l, r)) dp[r] = dp[l] + 1, from[r] = l;
    std::vector<StrSpan> parts;
    for (int r = n; r > 0; r = from[r]) parts.push_back({from[r], r});
    std::reverse(parts.begin(), parts.end());
    return parts;
}
template <class T, class Seq>
T count_palindromic_subsequences(const Seq& s0) {  // 回文である部分列（空でない。位置が違えば別）の数  O(N^2) 時間・O(N) メモリ
    const auto c = str_internal::codes_of(s0);
    const int n = (int)c.size();
    std::vector<T> nxt(n + 1, T(0)), cur(n + 1, T(0));  // f(i, j) = s[i..j] の中の数。nxt は i + 1 の行、添字は j + 1
    for (int i = n - 1; i >= 0; --i) {
        cur[i] = T(0), cur[i + 1] = T(1);  // f(i, i - 1) = 0, f(i, i) = 1
        for (int j = i + 1; j < n; ++j) {
            const T A = nxt[j + 1], B = cur[j], C = nxt[j];  // f(i+1, j), f(i, j-1), f(i+1, j-1)
            cur[j + 1] = c[i] == c[j] ? str_internal::checked_add(str_internal::checked_add(A, B, "count_palindromic_subsequences"), T(1), "count_palindromic_subsequences")
                                      : str_internal::checked_add(A, B - C, "count_palindromic_subsequences");  // B ≥ C なので先に引く（途中で溢れない）
        }
        std::swap(cur, nxt);
    }
    return n ? nxt[n] : T(0);
}
template <class T, class Seq>
T count_palindromes_leq(const Seq& s0, const Alphabet& A = Alphabet::lower()) {  // A の文字からなる長さ |s| の回文で s 以下（Alphabet の順）のものの数（ABC242 E）  O(N)
    const auto& s = str_internal::seq(s0);
    const int n = (int)s.size(), h = (n + 1) / 2;
    const char* func = "count_palindromes_leq";
    T r = T(0);
    for (int i = 0; i < h; ++i) r = str_internal::checked_add(str_internal::checked_mul(r, T(A.size()), func), T(A.id(s[i])), func);  // 前半が小さいもの
    int cmp = 0;  // 前半をそのまま折り返した回文と s の比較
    for (int i = h; i < n && cmp == 0; ++i) {
        const int a = A.id(s[n - 1 - i]), b = A.id(s[i]);
        cmp = a < b ? -1 : (a > b ? 1 : 0);
    }
    if (cmp <= 0) r = str_internal::checked_add(r, T(1), func);
    return r;
}
#endif  // GG_PALINDROME
