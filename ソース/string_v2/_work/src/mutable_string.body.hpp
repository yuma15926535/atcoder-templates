#ifndef GG_MUTABLE_STRING
#define GG_MUTABLE_STRING
// ---- 1 文字の書き換えがある文字列（文字列第二世代） ----
// ハッシュは 01_rolling_hash と同じ基数（RollingHash・hash_of で作ったほかの文字列とも比べられる）。フェニック木で O(log N)
// 文字の数・並び（count・counts・is_sorted・kth_char）の表は初めて使うときに作る（メモリ (N+1)σ 個の int。そのときに A の範囲を確かめる）
struct MutableString {
    MutableString() = default;
    template <class Seq>
    explicit MutableString(const Seq& s0, const Alphabet& A = Alphabet::lower()) : A_(A) {
        using C = str_internal::RHCore;
        const auto& s = str_internal::seq(s0);
        n_ = (int)s.size();
        C::ensure((std::size_t)n_ + 1);
        const unsigned long long inv = C::pow_mod(C::base1, C::MOD - 2);
        ipw_.assign(n_ + 1, 1);
        for (int i = 1; i <= n_; ++i) ipw_[i] = C::mul(ipw_[i - 1], inv);
        val_.resize(n_);
        fw_.assign(n_ + 1, 0), bw_.assign(n_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            val_[i] = str_internal::code_of(s[i]);
            const unsigned long long h = C::value(val_[i]);
            fw_[i + 1] = C::mul(h, ipw_[i]), bw_[i + 1] = C::mul(h, C::pw1[i]);
        }
        for (int i = 1; i <= n_; ++i) {  // フェニック木を O(N) で作る
            const int j = i + (i & -i);
            if (j <= n_) fw_[j] = C::add(fw_[j], fw_[i]), bw_[j] = C::add(bw_[j], bw_[i]);
        }
    }
    int size() const { return n_; }
    long long get(int i) const {  // i 文字目の値（char は 0〜255）  O(1)
        check_pos(i, "MutableString::get");
        return val_[i];
    }
    template <class T>
    void set(int i, T c) {  // i 文字目を c に  O(log N)
        using C = str_internal::RHCore;
        check_pos(i, "MutableString::set");
        const long long nv = str_internal::code_of(c), ov = val_[i];
        if (nv == ov) return;
        if (counts_built_) {
            const int oc = A_.id(ov), nc = A_.id(nv);
            bit_add(cnt_, oc, i, -1), bit_add(cnt_, nc, i, +1);
        }
        val_[i] = nv;
        if (desc_built_) {  // 前後との「下がる所」を付け直す
            for (int k = std::max(0, i - 1); k <= std::min(n_ - 2, i); ++k) {
                const int d = A_.id(val_[k]) > A_.id(val_[k + 1]) ? 1 : 0;
                if (d != desc_flag_[k]) bit_add_desc(k, d - desc_flag_[k]), desc_flag_[k] = (char)d;
            }
        }
        const unsigned long long diff = C::sub(C::value(nv), C::value(ov));
        for (int k = i + 1; k <= n_; k += k & -k) fw_[k] = C::add(fw_[k], C::mul(diff, ipw_[i])), bw_[k] = C::add(bw_[k], C::mul(diff, C::pw1[i]));
    }
    StrHash hash(int l, int r) const {  // s[l, r) のハッシュ（RollingHash::get と同じ値）  O(log N)
        str_internal::check_range("MutableString::hash", l, r, n_);
        return {raw(l, r), (long long)(r - l)};
    }
    StrHash reversed_hash(int l, int r) const {  // s[l, r) を反転した文字列のハッシュ  O(log N)
        str_internal::check_range("MutableString::reversed_hash", l, r, n_);
        return {raw_rev(l, r), (long long)(r - l)};
    }
    bool same(int l1, int r1, int l2, int r2) const {  // s[l1, r1) == s[l2, r2) か  O(log N)
        str_internal::check_range("MutableString::same", l1, r1, n_), str_internal::check_range("MutableString::same", l2, r2, n_);
        return r1 - l1 == r2 - l2 && hash(l1, r1) == hash(l2, r2);
    }
    bool is_palindrome(int l, int r) const {  // s[l, r) が回文か（ABC331 F）  O(log N)
        str_internal::check_range("MutableString::is_palindrome", l, r, n_);
        return raw(l, r) == raw_rev(l, r);
    }
    int lcp(int i, int j) const {  // s[i..) と s[j..) の最長共通接頭辞の長さ  O(log^2 N)
        if (i < 0 || i > n_ || j < 0 || j > n_) str_internal::fail("MutableString::lcp", "位置は 0〜N");
        int lo = 0, hi = std::min(n_ - i, n_ - j);
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (raw(i, i + mid) == raw(j, j + mid)) lo = mid;
            else hi = mid - 1;
        }
        return lo;
    }
    int compare(int l1, int r1, int l2, int r2) const {  // s[l1, r1) と s[l2, r2) の辞書順（値の大小で。-1 / 0 / 1）  O(log^2 N)
        str_internal::check_range("MutableString::compare", l1, r1, n_), str_internal::check_range("MutableString::compare", l2, r2, n_);
        const int m = std::min(r1 - l1, r2 - l2);
        int lo = 0, hi = m;
        while (lo < hi) {
            const int mid = (lo + hi + 1) >> 1;
            if (raw(l1, l1 + mid) == raw(l2, l2 + mid)) lo = mid;
            else hi = mid - 1;
        }
        if (lo == m) return (r1 - l1) == (r2 - l2) ? 0 : ((r1 - l1) < (r2 - l2) ? -1 : 1);
        return val_[l1 + lo] < val_[l2 + lo] ? -1 : 1;
    }

    // ---- 文字の数・並び（Alphabet の番号の順） ----
    template <class T>
    long long count(int l, int r, T c) const {  // s[l, r) にある c の数  O(log N)
        str_internal::check_range("MutableString::count", l, r, n_);
        ensure_counts();
        const int id = A_.id(c);
        return bit_sum(cnt_, id, r) - bit_sum(cnt_, id, l);
    }
    std::vector<long long> counts(int l, int r) const {  // s[l, r) の番号ごとの数（長さ σ）  O(σ log N)
        str_internal::check_range("MutableString::counts", l, r, n_);
        ensure_counts();
        std::vector<long long> res(A_.size());
        for (int c = 0; c < A_.size(); ++c) res[c] = bit_sum(cnt_, c, r) - bit_sum(cnt_, c, l);
        return res;
    }
    bool is_sorted(int l, int r) const {  // s[l, r) が広義単調増加（Alphabet の番号の順）か（ABC285 F）  O(log N)
        str_internal::check_range("MutableString::is_sorted", l, r, n_);
        ensure_desc();
        if (r - l <= 1) return true;
        return desc_sum(r - 1) - desc_sum(l) == 0;  // 下がる所 k（s[k] > s[k+1]）が [l, r-1) にない
    }
    long long kth_char(int l, int r, int k) const {  // s[l, r) を並べ替えた k 番目（0 から）の値  O(σ log N)
        str_internal::check_range("MutableString::kth_char", l, r, n_);
        if (k < 0 || k >= r - l) str_internal::fail("MutableString::kth_char", "k は 0〜(r - l - 1)");
        ensure_counts();
        for (int c = 0; c < A_.size(); ++c) {
            const int x = bit_sum(cnt_, c, r) - bit_sum(cnt_, c, l);
            if (k < x) return A_.value(c);
            k -= x;
        }
        return -1;  // ここには来ない
    }
    std::vector<long long> values() const { return val_; }  // 今の値の列
    std::string to_string() const {                         // 今の文字列（値が 0〜255 のとき）
        std::string r(n_, '\0');
        for (int i = 0; i < n_; ++i) {
            if (val_[i] < 0 || val_[i] > 255) str_internal::fail("MutableString::to_string", "値が文字（0〜255）でない");
            r[i] = (char)(unsigned char)val_[i];
        }
        return r;
    }
    const Alphabet& alphabet() const { return A_; }

   private:
    unsigned long long raw(int l, int r) const {  // Σ v_i B^(r-1-i)
        using C = str_internal::RHCore;
        if (l >= r) return 0;
        return C::mul(C::sub(fsum(fw_, r), fsum(fw_, l)), C::pw1[r - 1]);
    }
    unsigned long long raw_rev(int l, int r) const {  // Σ v_i B^(i-l)
        using C = str_internal::RHCore;
        if (l >= r) return 0;
        return C::mul(C::sub(fsum(bw_, r), fsum(bw_, l)), ipw_[l]);
    }
    static unsigned long long fsum(const std::vector<unsigned long long>& t, int i) {  // t の [0, i) の和
        unsigned long long s = 0;
        for (; i > 0; i -= i & -i) s = str_internal::RHCore::add(s, t[i]);
        return s;
    }
    void check_pos(int i, const char* func) const {
        if (i < 0 || i >= n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜" + std::to_string(n_ - 1) + " の外");
    }
    // 文字ごとのフェニック木（cnt_[c * (N+1) + i]）
    void ensure_counts() const {
        if (counts_built_) return;
        counts_built_ = true;
        const std::size_t W = (std::size_t)n_ + 1;
        if ((long long)A_.size() * (long long)W > 200000000LL) str_internal::fail("MutableString", "σ × (N+1) が 2×10^8（約 800MB）を超える（文字の数の表を作れない）");
        cnt_.assign((std::size_t)A_.size() * W, 0);
        for (int i = 0; i < n_; ++i) ++cnt_[(std::size_t)A_.id(val_[i]) * W + i + 1];
        for (int c = 0; c < A_.size(); ++c) {
            int* t = cnt_.data() + (std::size_t)c * W;
            for (int i = 1; i <= n_; ++i) {
                const int j = i + (i & -i);
                if (j <= n_) t[j] += t[i];
            }
        }
    }
    int bit_sum(const std::vector<int>& t, int c, int i) const {
        if (n_ == 0) return 0;
        const int* b = t.data() + (std::size_t)c * ((std::size_t)n_ + 1);
        int s = 0;
        for (; i > 0; i -= i & -i) s += b[i];
        return s;
    }
    void bit_add(std::vector<int>& t, int c, int i, int d) {
        int* b = t.data() + (std::size_t)c * ((std::size_t)n_ + 1);
        for (int k = i + 1; k <= n_; k += k & -k) b[k] += d;
    }
    // 下がる所（s[k] > s[k+1]）のフェニック木
    void ensure_desc() const {
        if (desc_built_) return;
        desc_built_ = true;
        desc_.assign(n_ + 1, 0), desc_flag_.assign(std::max(0, n_ - 1), 0);
        for (int k = 0; k + 1 < n_; ++k) desc_flag_[k] = A_.id(val_[k]) > A_.id(val_[k + 1]) ? 1 : 0, desc_[k + 1] = desc_flag_[k];
        for (int i = 1; i <= n_; ++i) {
            const int j = i + (i & -i);
            if (j <= n_) desc_[j] += desc_[i];
        }
    }
    int desc_sum(int i) const {
        int s = 0;
        for (; i > 0; i -= i & -i) s += desc_[i];
        return s;
    }
    void bit_add_desc(int k, int d) {
        for (int x = k + 1; x <= n_; x += x & -x) desc_[x] += d;
    }
    Alphabet A_;
    int n_ = 0;
    std::vector<long long> val_;
    std::vector<unsigned long long> ipw_ = std::vector<unsigned long long>(1, 1), fw_ = std::vector<unsigned long long>(1, 0), bw_ = std::vector<unsigned long long>(1, 0);
    mutable std::vector<int> cnt_, desc_;
    mutable std::vector<char> desc_flag_;
    mutable bool counts_built_ = false, desc_built_ = false;
};
#endif  // GG_MUTABLE_STRING
