#ifndef GG_ROLLING_HASH
#define GG_ROLLING_HASH
// ---- ローリングハッシュ（文字列第二世代）  mod 2^61 - 1。基数は実行ごとの乱数で、すべての文字列で共通 ----
// 別の文字列の部分文字列どうし・つなげたものも比べられる。累乗の表は基数を変えずに伸ばす
// 1 回の比較で衝突する確率は約 N / 2^61（N は長さ）。k 文字まで書き換えてよい回文の長さ（center_palindrome_lengths_with_changes）もここ
namespace str_internal {
__extension__ typedef unsigned __int128 rh_u128;
struct RHCore {
    static constexpr unsigned long long MOD = (1ULL << 61) - 1;
    static unsigned long long mul(unsigned long long a, unsigned long long b) {
        const rh_u128 t = (rh_u128)a * b;
        const unsigned long long r = (unsigned long long)(t >> 61) + ((unsigned long long)t & MOD);
        return r >= MOD ? r - MOD : r;
    }
    static unsigned long long add(unsigned long long a, unsigned long long b) {
        const unsigned long long r = a + b;
        return r >= MOD ? r - MOD : r;
    }
    static unsigned long long sub(unsigned long long a, unsigned long long b) { return a >= b ? a - b : a + MOD - b; }
    static unsigned long long splitmix(unsigned long long x) {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }
    inline static unsigned long long seed = 0, base1 = 0, base2 = 0;
    inline static std::vector<unsigned long long> pw1, pw2;
    static void init() {
        if (base1) return;
        std::random_device rd;
        seed = splitmix(((unsigned long long)rd() << 32) ^ rd() ^ (unsigned long long)std::chrono::steady_clock::now().time_since_epoch().count());
        const unsigned long long lo = 1ULL << 40, w = MOD - 2 * lo;
        base1 = lo + splitmix(seed ^ 1) % w;
        do base2 = lo + splitmix(seed ^ base1) % w;
        while (base2 == base1);
        pw1.assign(1, 1), pw2.assign(1, 1);
    }
    static void ensure(std::size_t n) {  // pw1 の大きさを n + 1 以上に（基数は変えない）
        init();
        while (pw1.size() <= n) pw1.push_back(mul(pw1.back(), base1));
    }
    static void ensure2(std::size_t n) {  // pw2 も（2 次元用）
        ensure(n);
        while (pw2.size() <= n) pw2.push_back(mul(pw2.back(), base2));
    }
    static unsigned long long pow_mod(unsigned long long b, unsigned long long e) {
        unsigned long long r = 1;
        for (; e; e >>= 1, b = mul(b, b))
            if (e & 1) r = mul(r, b);
        return r;
    }
    static unsigned long long pow1(long long n) {  // base1^n（n が表の外なら繰り返し二乗）
        init();
        return (unsigned long long)n < pw1.size() ? pw1[(std::size_t)n] : pow_mod(base1, (unsigned long long)n);
    }
    static unsigned long long pow2(long long n) {
        init();
        return (unsigned long long)n < pw2.size() ? pw2[(std::size_t)n] : pow_mod(base2, (unsigned long long)n);
    }
    // 要素の値 → [1, MOD)。0〜2^60-1 はそのまま +1、それ以外は乱数で [2^60+1, MOD) に散らす
    static unsigned long long value(long long v) {
        if ((unsigned long long)v < (1ULL << 60)) return (unsigned long long)v + 1;
        init();
        return (1ULL << 60) + 1 + splitmix((unsigned long long)v ^ seed) % (MOD - (1ULL << 60) - 1);
    }
};
}  // namespace str_internal

// ---- 部分文字列のハッシュ値（長さも持つ。+ でつなげられる） ----
struct StrHash {
    unsigned long long value = 0;
    long long len = 0;
    friend StrHash operator+(const StrHash& a, const StrHash& b) {  // a のあとに b をつなげたもののハッシュ  O(1)（b.len が表の外なら O(log)）
        using C = str_internal::RHCore;
        return {C::add(C::mul(a.value, C::pow1(b.len)), b.value), a.len + b.len};
    }
    StrHash& operator+=(const StrHash& b) { return *this = *this + b; }
    StrHash repeat(long long k) const {  // k 回つなげたもの  O(log k)
        using C = str_internal::RHCore;
        if (k < 0 || (len > 0 && k > 4000000000000000000LL / len)) str_internal::fail("StrHash::repeat", "k ≥ 0 かつ長さ × k ≤ 4×10^18 で使う");
        // sum = 1 + P + … + P^(k-1)（P = base^len）を倍々で
        unsigned long long sum = 0, cur_sum = 1, cur_p = C::pow1(len);  // cur_sum = 1 + P + … + P^(2^b - 1)、cur_p = P^(2^b)
        for (long long e = k; e > 0; e >>= 1) {
            if (e & 1) sum = C::add(C::mul(sum, cur_p), cur_sum);
            cur_sum = C::add(cur_sum, C::mul(cur_sum, cur_p)), cur_p = C::mul(cur_p, cur_p);
        }
        return {C::mul(value, sum), len * k};
    }
    friend bool operator==(const StrHash& a, const StrHash& b) { return a.value == b.value && a.len == b.len; }
    friend bool operator!=(const StrHash& a, const StrHash& b) { return !(a == b); }
    friend bool operator<(const StrHash& a, const StrHash& b) { return a.len != b.len ? a.len < b.len : a.value < b.value; }  // map・sort 用（辞書順ではない）
    friend std::ostream& operator<<(std::ostream& os, const StrHash& h) { return os << "#" << h.value << "(len " << h.len << ")"; }  // デバッグ用
};
namespace std {
template <>
struct hash<StrHash> {  // unordered_map<StrHash, V> / unordered_set<StrHash> にそのまま使える
    size_t operator()(const StrHash& h) const noexcept { return (size_t)(h.value ^ ((unsigned long long)h.len * 0x9E3779B97F4A7C15ULL)); }
};
}  // namespace std

template <class Seq>
StrHash hash_of(const Seq& s0) {  // 列全体のハッシュ  O(N)
    using C = str_internal::RHCore;
    const auto& s = str_internal::seq(s0);
    C::init();
    unsigned long long h = 0;
    for (std::size_t i = 0; i < s.size(); ++i) h = C::add(C::mul(h, C::base1), C::value(str_internal::code_of(s[i])));
    return {h, (long long)s.size()};
}
template <class T>
StrHash hash_of_char(T c) { return {str_internal::RHCore::value(str_internal::code_of(c)), 1}; }  // 1 文字（1 要素）のハッシュ
inline StrHash empty_hash() { return {0, 0}; }

// ---- 列の部分文字列のハッシュを O(1) で ----
struct RollingHash {
    RollingHash() : h_(1, 0) { str_internal::RHCore::ensure(0); }  // 空の文字列と同じ
    template <class Seq>
    explicit RollingHash(const Seq& s0) {
        using C = str_internal::RHCore;
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        C::ensure((std::size_t)n_);
        h_.resize(n_ + 1);
        h_[0] = 0;
        // 比べるときに使う元の値は数の列のときだけ持つ（char の値 0〜255 はハッシュから戻せる：1 文字のハッシュ = 値 + 1）
        constexpr bool narrow = std::is_same<str_internal::elem_t<Seq>, char>::value || std::is_same<str_internal::elem_t<Seq>, signed char>::value ||
                                std::is_same<str_internal::elem_t<Seq>, unsigned char>::value;
        if (!narrow) c64_.resize(n_);
        const unsigned long long b = C::base1;
        for (int i = 0; i < n_; ++i) {
            const long long c = str_internal::code_of(s[i]);
            if (!narrow) c64_[i] = c;
            h_[i + 1] = C::add(C::mul(h_[i], b), C::value(c));
        }
    }
    int size() const { return n_; }
    long long code(int i) const {  // i 文字目の値（char は 0〜255）
        if (i < 0 || i >= n_) str_internal::fail("RollingHash::code", "位置 " + std::to_string(i) + " は 0〜N-1 の外");
        return code_at(i);
    }
    StrHash get(int l, int r) const {  // s[l, r)  O(1)
        str_internal::check_range("RollingHash::get", l, r, n_);
        return {raw(l, r), (long long)(r - l)};
    }
    StrHash get_len(int l, int len) const { return get(l, l + len); }  // s[l, l + len)
    StrHash prefix(int r) const { return get(0, r); }                  // s[0, r)
    StrHash suffix(int l) const { return get(l, n_); }                 // s[l, N)
    StrHash all() const { return get(0, n_); }
    bool same(int l1, int r1, int l2, int r2) const { return get(l1, r1) == get(l2, r2); }  // s[l1, r1) == s[l2, r2) か  O(1)
    int lcp(int i, int j) const { return lcp(i, *this, j); }  // s[i..) と s[j..) の最長共通接頭辞の長さ  O(log N)
    int lcp(int i, const RollingHash& t, int j) const {       // s[i..) と t[j..) の最長共通接頭辞の長さ  O(log N)
        check_pos("RollingHash::lcp", i), t.check_pos("RollingHash::lcp", j);
        int lo = 0, hi = std::min(n_ - i, t.n_ - j);
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (raw(i, i + mid) == t.raw(j, j + mid)) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    }
    int common_suffix(int i, int j) const { return common_suffix(i, *this, j); }  // s[0, i) と s[0, j) の最長共通接尾辞の長さ  O(log N)
    int common_suffix(int i, const RollingHash& t, int j) const {
        check_pos("RollingHash::common_suffix", i), t.check_pos("RollingHash::common_suffix", j);
        int lo = 0, hi = std::min(i, j);
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (raw(i - mid, i) == t.raw(j - mid, j)) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    }
    int compare(int l1, int r1, int l2, int r2) const { return compare(l1, r1, *this, l2, r2); }  // s[l1, r1) と s[l2, r2) の辞書順（-1 / 0 / 1）  O(log N)
    int compare(int l1, int r1, const RollingHash& t, int l2, int r2) const {                      // s[l1, r1) と t[l2, r2)
        str_internal::check_range("RollingHash::compare", l1, r1, n_), str_internal::check_range("RollingHash::compare", l2, r2, t.n_);
        const int m = std::min(r1 - l1, r2 - l2);
        int lo = 0, hi = m;
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (raw(l1, l1 + mid) == t.raw(l2, l2 + mid)) lo = mid;
            else hi = mid - 1;
        }
        if (lo == m) return (r1 - l1) == (r2 - l2) ? 0 : ((r1 - l1) < (r2 - l2) ? -1 : 1);
        return code_at(l1 + lo) < t.code_at(l2 + lo) ? -1 : 1;
    }
    std::vector<int> find_all(const StrHash& p) const {  // ハッシュが p と同じ部分文字列の開始位置（昇順）  O(N)
        std::vector<int> r;
        if (p.len > n_) return r;
        const int m = (int)p.len;
        for (int i = 0; i + m <= n_; ++i)
            if (raw(i, i + m) == p.value) r.push_back(i);
        return r;
    }
    template <class Seq>
    std::vector<int> find_all(const Seq& pat) const { return find_all(hash_of(pat)); }  // pat が現れる開始位置（重なりも）
    template <class Seq>
    long long count(const Seq& pat) const { return (long long)find_all(pat).size(); }

   private:
    void check_pos(const char* func, int i) const {
        if (i < 0 || i > n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜N（= " + std::to_string(n_) + "）の外");
    }
    unsigned long long raw(int l, int r) const {  // 範囲の確認なしのハッシュ値
        using C = str_internal::RHCore;
        return C::sub(h_[r], C::mul(h_[l], C::pw1[(std::size_t)(r - l)]));
    }
    long long code_at(int i) const { return c64_.empty() ? (long long)raw(i, i + 1) - 1 : c64_[i]; }
    int n_ = 0;
    std::vector<unsigned long long> h_;  // h_[i] = s[0, i) のハッシュ
    std::vector<long long> c64_;         // 元の値（数の列のとき。char の列なら空）
};

// ---- k 文字まで書き換えてよい回文：中心 c（0〜2N-2。偶数 c は s[c/2] が中心の奇数長、奇数 c は間が中心の偶数長）ごとの最長の長さ ----
//   中心から広げ、食い違いを k 回まで飛ばす（ハッシュの LCP を k + 1 回）。k = 0 なら 04_palindrome の center_palindrome_lengths と同じ（ABC468 D）
//   中心 c で書き換え k 回以内の回文は、長さ L[c] 以下で同じ偶奇のもの全部  O(N k log N)（k が小さいとき用。N = 2×10^5 なら k は 10 程度まで）
template <class Seq>
std::vector<int> center_palindrome_lengths_with_changes(const Seq& s0, int k) {
    if (k < 0) str_internal::fail("center_palindrome_lengths_with_changes", "k は 0 以上");
    const auto& s = str_internal::seq(s0);
    const int n = (int)s.size();
    if (n == 0) return {};
    const RollingHash F(s), B(reversed(s));
    std::vector<int> res(2 * n - 1);
    for (int c = 0; c < 2 * n - 1; ++c) {
        int l = c % 2 == 0 ? c / 2 - 1 : (c - 1) / 2, r = c % 2 == 0 ? c / 2 + 1 : (c + 1) / 2;  // 次に比べる左右の位置
        int len = c % 2 == 0 ? 1 : 0;
        for (int used = 0;; ++used) {
            if (l < 0 || r >= n) break;
            const int e = F.lcp(r, B, n - 1 - l);  // s[r..) と s[..l] を逆から読んだものの一致
            l -= e, r += e, len += 2 * e;
            if (l < 0 || r >= n || used == k) break;
            --l, ++r, len += 2;  // 食い違いを 1 つ書き換える
        }
        res[c] = len;
    }
    return res;
}

// ---- 2 次元のハッシュ（グリッドの長方形 [i1, i2) × [j1, j2) を O(1) で） ----
struct GridHash {
    unsigned long long value = 0;
    int h = 0, w = 0;
    friend bool operator==(const GridHash& a, const GridHash& b) { return a.value == b.value && a.h == b.h && a.w == b.w; }
    friend bool operator!=(const GridHash& a, const GridHash& b) { return !(a == b); }
    friend bool operator<(const GridHash& a, const GridHash& b) { return std::tie(a.h, a.w, a.value) < std::tie(b.h, b.w, b.value); }
};
namespace std {
template <>
struct hash<GridHash> {
    size_t operator()(const GridHash& g) const noexcept { return (size_t)(g.value ^ ((unsigned long long)g.h * 0x9E3779B97F4A7C15ULL) ^ ((unsigned long long)g.w << 32)); }
};
}  // namespace std
struct RollingHash2D {
    RollingHash2D() : P_(1, 0) { str_internal::RHCore::ensure2(0); }  // 0 × 0 のグリッドと同じ
    template <class Grid>
    explicit RollingHash2D(const Grid& g) {  // vector<string> でも vector<vector<int>> でも（行の長さはそろえる）
        using C = str_internal::RHCore;
        H_ = (int)g.size(), W_ = H_ ? (int)g[0].size() : 0;
        C::ensure2((std::size_t)std::max(H_, W_));
        P_.assign((std::size_t)(H_ + 1) * (W_ + 1), 0);
        const unsigned long long b12 = C::mul(C::base1, C::base2);
        for (int i = 0; i < H_; ++i) {
            if ((int)g[i].size() != W_) str_internal::fail("RollingHash2D", "行の長さがそろっていない");
            for (int j = 0; j < W_; ++j) {
                unsigned long long v = C::add(C::mul(at(i, j + 1), C::base2), C::mul(at(i + 1, j), C::base1));
                v = C::sub(v, C::mul(at(i, j), b12));
                at(i + 1, j + 1) = C::add(v, C::value(str_internal::code_of(g[i][j])));
            }
        }
    }
    int rows() const { return H_; }
    int cols() const { return W_; }
    GridHash get(int i1, int j1, int i2, int j2) const {  // 行 [i1, i2) × 列 [j1, j2)  O(1)
        str_internal::check_range("RollingHash2D::get", i1, i2, H_), str_internal::check_range("RollingHash2D::get", j1, j2, W_);
        using C = str_internal::RHCore;
        const unsigned long long ph = C::pw2[(std::size_t)(i2 - i1)], pw = C::pw1[(std::size_t)(j2 - j1)];
        unsigned long long v = C::sub(cat(i2, j2), C::mul(cat(i1, j2), ph));
        v = C::sub(v, C::mul(cat(i2, j1), pw));
        v = C::add(v, C::mul(C::mul(cat(i1, j1), ph), pw));
        return {v, i2 - i1, j2 - j1};
    }
    GridHash all() const { return get(0, 0, H_, W_); }
    template <class Grid>
    std::vector<std::pair<int, int>> find_all(const Grid& pat) const {  // pat と一致する左上の位置 (i, j) の一覧（行の順・列の順）  O(HW)
        std::vector<std::pair<int, int>> r;
        const RollingHash2D P(pat);
        const int h = P.rows(), w = P.cols();
        if (h > H_ || w > W_) return r;
        const GridHash want = P.all();
        for (int i = 0; i + h <= H_; ++i)
            for (int j = 0; j + w <= W_; ++j)
                if (get(i, j, i + h, j + w) == want) r.push_back({i, j});
        return r;
    }

   private:
    unsigned long long& at(int i, int j) { return P_[(std::size_t)i * (W_ + 1) + j]; }
    unsigned long long cat(int i, int j) const { return P_[(std::size_t)i * (W_ + 1) + j]; }
    int H_ = 0, W_ = 0;
    std::vector<unsigned long long> P_;  // P_[i][j] = 左上 [0, i) × [0, j) のハッシュ
};
#endif  // GG_ROLLING_HASH
