#ifndef GG_SET_HASH
#define GG_SET_HASH
// ---- 集合・多重集合のハッシュ（文字列第二世代） ----
// 値ごとに乱数（実行ごとに変わる）を割り当て、集合は「入っている値の乱数の和」、多重集合は「個数 × 乱数の和」を mod 2^61 - 1 で持つ
// 並べ替えて同じか・集合として同じかを O(1) で比べられる（違うのに等しくなる確率は 1 回あたり約 2^-61）。表（値 → 乱数）は持たない
namespace str_internal {
__extension__ typedef unsigned __int128 sh_u128;
struct SetHashCore {
    static constexpr unsigned long long MOD = (1ULL << 61) - 1;
    static unsigned long long add(unsigned long long a, unsigned long long b) {
        const unsigned long long r = a + b;
        return r >= MOD ? r - MOD : r;
    }
    static unsigned long long sub(unsigned long long a, unsigned long long b) { return a >= b ? a - b : a + MOD - b; }
    static unsigned long long mul(unsigned long long a, unsigned long long b) {
        const sh_u128 t = (sh_u128)a * b;
        const unsigned long long r = (unsigned long long)(t >> 61) + ((unsigned long long)t & MOD);
        return r >= MOD ? r - MOD : r;
    }
    static unsigned long long mix(unsigned long long x) {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    }
    inline static unsigned long long seed = 0;
    static unsigned long long of(long long v) {  // 値 → [1, MOD) の乱数
        if (seed == 0) {
            std::random_device rd;
            seed = mix(((unsigned long long)rd() << 32) ^ rd() ^ (unsigned long long)std::chrono::steady_clock::now().time_since_epoch().count()) | 1;
        }
        const unsigned long long h = mix(mix((unsigned long long)v ^ seed) + seed) % (MOD - 1);
        return h + 1;
    }
    static unsigned long long times(unsigned long long h, long long k) {  // h × k（k は負でもよい）
        const unsigned long long kk = k >= 0 ? (unsigned long long)k % MOD : MOD - (unsigned long long)(-(k + 1)) % MOD - 1;
        return mul(h, kk);
    }
};
struct SetHashKey {  // unordered_map の鍵のハッシュ（実行ごとの乱数。std::hash<long long> は値そのままで、狙った入力で遅くなる）
    std::size_t operator()(long long v) const {
        static const unsigned long long salt = SetHashCore::mix((unsigned long long)std::chrono::steady_clock::now().time_since_epoch().count() ^ ((unsigned long long)std::random_device()() << 32));
        return (std::size_t)SetHashCore::mix((unsigned long long)v + salt);
    }
};
}  // namespace str_internal

template <class T>
unsigned long long value_hash(T x) { return str_internal::SetHashCore::of(str_internal::code_of(x)); }  // 値 x の乱数（[1, 2^61 - 1)。同じ実行の中では同じ値に同じ乱数）

// ---- 先頭からの集合・多重集合（区間の多重集合は差で O(1)） ----
template <class Seq>
std::vector<unsigned long long> prefix_set_hashes(const Seq& a0) {  // r[i] = a[0, i) に現れる値の「集合」のハッシュ（重複は 1 回。長さ N + 1。ABC250 E）  O(N log N)
    using C = str_internal::SetHashCore;
    const auto& a = str_internal::seq(a0);
    const int n = (int)a.size();
    std::vector<std::pair<long long, int>> vi(n);  // 並べて、それぞれの値が最初に現れる位置を見つける（ハッシュ表を使わない）
    for (int i = 0; i < n; ++i) vi[i] = {str_internal::code_of(a[i]), i};
    std::sort(vi.begin(), vi.end());
    std::vector<unsigned long long> r(n + 1, 0);
    for (int k = 0; k < n; ++k)
        if (k == 0 || vi[k].first != vi[k - 1].first) r[vi[k].second + 1] = C::of(vi[k].first);
    for (int i = 0; i < n; ++i) r[i + 1] = C::add(r[i], r[i + 1]);
    return r;
}
template <class Seq>
std::vector<unsigned long long> prefix_multiset_hashes(const Seq& a0) {  // r[i] = a[0, i) の「多重集合」のハッシュ（長さ N + 1）  O(N)
    using C = str_internal::SetHashCore;
    const auto& a = str_internal::seq(a0);
    std::vector<unsigned long long> r(a.size() + 1, 0);
    for (std::size_t i = 0; i < a.size(); ++i) r[i + 1] = C::add(r[i], C::of(str_internal::code_of(a[i])));
    return r;
}
struct MultisetPrefix {  // 区間 [l, r) の多重集合のハッシュを O(1) で（ABC367 F）
    MultisetPrefix() = default;
    template <class Seq>
    explicit MultisetPrefix(const Seq& a) : p_(prefix_multiset_hashes(a)) {}
    int size() const { return (int)p_.size() - 1; }
    unsigned long long get(int l, int r) const {  // a[l, r) の多重集合のハッシュ（並べ替えて同じ ⇔ 等しい）  O(1)
        str_internal::check_range("MultisetPrefix::get", l, r, size());
        return str_internal::SetHashCore::sub(p_[r], p_[l]);
    }
    bool same(int l1, int r1, int l2, int r2) const {  // 並べ替えて同じか  O(1)
        str_internal::check_range("MultisetPrefix::same", l1, r1, size()), str_internal::check_range("MultisetPrefix::same", l2, r2, size());
        return r1 - l1 == r2 - l2 && get(l1, r1) == get(l2, r2);
    }

   private:
    std::vector<unsigned long long> p_ = std::vector<unsigned long long>(1, 0);
};

// ---- 追加・削除できる集合・多重集合 ----
struct MultisetHash {  // 多重集合（値の個数を持たない。足す・引く・合わせるだけ。個数 k は |k| < 2^61 − 1）
    unsigned long long value = 0;  // ハッシュ（空なら 0）
    long long count = 0;           // 要素の数
    template <class T>
    void add(T x, long long k = 1) {  // x を k 個足す（k が負なら取り除く）  O(1)
        using C = str_internal::SetHashCore;
        value = C::add(value, C::times(C::of(str_internal::code_of(x)), k)), count += k;
    }
    template <class T>
    void remove(T x, long long k = 1) {  // x を k 個取り除く（入っているかは確かめない）  O(1)
        using C = str_internal::SetHashCore;
        value = C::sub(value, C::times(C::of(str_internal::code_of(x)), k)), count -= k;
    }
    void merge(const MultisetHash& o) { value = str_internal::SetHashCore::add(value, o.value), count += o.count; }  // o を全部足す
    friend bool operator==(const MultisetHash& a, const MultisetHash& b) { return a.value == b.value && a.count == b.count; }
    friend bool operator!=(const MultisetHash& a, const MultisetHash& b) { return !(a == b); }
};
struct SetHash {  // 集合（値ごとの個数を持ち、0 個 ⇔ 1 個以上 が変わったときだけハッシュが変わる）
    unsigned long long value = 0;  // 入っている異なる値の集合のハッシュ
    template <class T>
    void add(T x) {  // x を 1 つ足す（すでにあれば集合は変わらない）  O(1)（平均）
        const long long v = str_internal::code_of(x);
        if (cnt_[v]++ == 0) value = str_internal::SetHashCore::add(value, str_internal::SetHashCore::of(v)), ++distinct_;
    }
    template <class T>
    bool remove(T x) {  // x を 1 つ取り除く（入っていなければ何もせず false）  O(1)（平均）
        const long long v = str_internal::code_of(x);
        auto it = cnt_.find(v);
        if (it == cnt_.end()) return false;
        if (--it->second == 0) cnt_.erase(it), value = str_internal::SetHashCore::sub(value, str_internal::SetHashCore::of(v)), --distinct_;
        return true;
    }
    template <class T>
    long long count(T x) const {  // x がいくつ足されているか
        auto it = cnt_.find(str_internal::code_of(x));
        return it == cnt_.end() ? 0 : it->second;
    }
    long long distinct() const { return distinct_; }  // 異なる値の数
    friend bool operator==(const SetHash& a, const SetHash& b) { return a.value == b.value && a.distinct_ == b.distinct_; }
    friend bool operator!=(const SetHash& a, const SetHash& b) { return !(a == b); }

   private:
    std::unordered_map<long long, long long, str_internal::SetHashKey> cnt_;
    long long distinct_ = 0;
};
#endif  // GG_SET_HASH
