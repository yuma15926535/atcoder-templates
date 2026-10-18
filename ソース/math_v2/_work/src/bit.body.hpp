#ifndef GG_BIT
#define GG_BIT
// ---- ビット（数学第二世代）  ビットの数・位置、部分集合の列挙、グレイコード、0〜n の XOR・ビットの数、XOR の和 ----
//   x は 64 ビットとして扱う（負の数は 2 の補数のビット）。位置 i は 0〜63（0 が最下位）
inline int bit_count(long long x) { return __builtin_popcountll((unsigned long long)x); }                // 1 の数
inline int top_bit(long long x) { return x == 0 ? -1 : 63 - __builtin_clzll((unsigned long long)x); }   // 最上位の 1 の位置（0 なら −1）
inline int low_bit(long long x) { return x == 0 ? -1 : __builtin_ctzll((unsigned long long)x); }       // 最下位の 1 の位置（0 なら −1。値がほしいなら x & -x）
inline int bit_length(long long x) { return top_bit(x) + 1; }                                          // 2 進数の桁数（0 なら 0）
inline bool has_bit(long long x, int i) {                                                              // i ビット目が 1 か
    if (i < 0 || i > 63) math_internal::fail("has_bit", "i は 0〜63");
    return ((unsigned long long)x >> i) & 1;
}
inline std::vector<int> bit_indices(long long x) {  // 1 のビットの位置（小さい順）
    std::vector<int> r;
    for (unsigned long long u = (unsigned long long)x; u; u &= u - 1) r.push_back(__builtin_ctzll(u));
    return r;
}
template <class F>
void for_each_bit(long long x, F f) {  // 1 のビットの位置 i を小さい順に f(i)
    for (unsigned long long u = (unsigned long long)x; u; u &= u - 1) f(__builtin_ctzll(u));
}
inline long long bit_reverse(long long x, int n) {  // 下 n ビットを左右反転（n は 0〜63。上のビットは 0 になる）
    if (n < 0 || n > 63) math_internal::fail("bit_reverse", "n は 0〜63");
    unsigned long long r = 0;
    for (int i = 0; i < n; ++i) r |= (((unsigned long long)x >> i) & 1) << (n - 1 - i);
    return (long long)r;
}
inline std::string mask_to_string(long long mask, int n) {  // 下 n ビットを "0101" の形に（上の桁が左。n は 0〜64）
    if (n < 0 || n > 64) math_internal::fail("mask_to_string", "n は 0〜64");
    std::string s((size_t)n, '0');
    for (int i = 0; i < n; ++i)
        if (((unsigned long long)mask >> i) & 1) s[(size_t)(n - 1 - i)] = '1';
    return s;
}

// ---- 部分集合の列挙（小さい順） ----
namespace math_internal {
inline void check_width(int n, const char* func) {
    if (n < 0 || n > 63) fail(func, "n は 0〜63");
}
}  // namespace math_internal
template <class F>
void for_each_submask(long long mask, F f, bool include_zero = true, bool include_self = true) {  // s ⊆ mask を小さい順に f(s)（3^n 通りの部分集合 DP にも）
    if (mask < 0) math_internal::fail("for_each_submask", "mask ≥ 0 で使う");
    const unsigned long long m = (unsigned long long)mask;
    unsigned long long s = 0;
    do {
        if ((s != 0 || include_zero) && (s != m || include_self)) f((long long)s);
        s = (s - m) & m;  // 次に大きい部分集合
    } while (s != 0);
}
template <class F>
void for_each_supermask(long long mask, int n, F f) {  // mask ⊆ s ⊆ {0, …, n−1} を小さい順に f(s)
    math_internal::check_width(n, "for_each_supermask");
    if (mask < 0 || (unsigned long long)mask >> n) math_internal::fail("for_each_supermask", "mask は 0 以上 2^n 未満");
    const unsigned long long full = (1ULL << n) - 1;
    for (unsigned long long s = (unsigned long long)mask;; s = (s + 1) | (unsigned long long)mask) {
        f((long long)s);
        if (s == full) break;
    }
}
template <class F>
void for_each_mask_between(long long must, long long allowed, F f) {  // must ⊆ s ⊆ allowed を小さい順に f(s)（必ず 1 のビットと使ってよいビット。must ⊄ allowed なら何もしない）
    if (must < 0 || allowed < 0) math_internal::fail("for_each_mask_between", "must, allowed ≥ 0 で使う");
    if (must & ~allowed) return;
    for_each_submask(allowed & ~must, [&](long long s) { f(must | s); });
}
template <class F>
void for_each_mask_between(long long must, long long allowed, int k, F f) {  // さらに 1 の数がちょうど k のものだけ（小さい順）
    if (must < 0 || allowed < 0) math_internal::fail("for_each_mask_between", "must, allowed ≥ 0 で使う");
    if (must & ~allowed) return;
    const std::vector<int> pos = bit_indices(allowed & ~must);  // 自由なビット
    const int c = (int)pos.size(), t = k - bit_count(must);
    if (t < 0 || t > c) return;
    if (t == 0) { f(must); return; }
    // 自由なビットの番号の上でゴスパーの方法（allowed ≥ 0 なので c ≤ 63）
    const unsigned long long end = 1ULL << c;
    for (unsigned long long x = (1ULL << t) - 1; x < end;) {
        unsigned long long s = (unsigned long long)must;
        for (unsigned long long u = x; u; u &= u - 1) s |= 1ULL << pos[(size_t)__builtin_ctzll(u)];
        f((long long)s);
        const unsigned long long lo = x & (0 - x), y = x + lo;
        x = (((x ^ y) >> 2) / lo) | y;
    }
}
template <class F>
void for_each_subset_k(int n, int k, F f) {  // {0, …, n−1} の k 個の部分集合（1 の数が k の n ビットの数）を小さい順に f(s)
    math_internal::check_width(n, "for_each_subset_k");
    if (k < 0 || k > n) return;
    if (k == 0) { f(0LL); return; }
    const unsigned long long end = 1ULL << n;
    for (unsigned long long x = (1ULL << k) - 1; x < end;) {
        f((long long)x);
        const unsigned long long lo = x & (0 - x), y = x + lo;
        x = (((x ^ y) >> 2) / lo) | y;
    }
}

// ---- グレイコード（隣どうしが 1 ビットだけ違う並び） ----
inline long long gray_code(long long i) { return i ^ (long long)((unsigned long long)i >> 1); }  // i 番目のグレイコード
inline long long gray_inverse(long long g) {                                                       // gray_code の逆（何番目か）
    unsigned long long u = (unsigned long long)g;
    for (int s = 1; s < 64; s <<= 1) u ^= u >> s;
    return (long long)u;
}

// ---- 0〜n の XOR・ビットの数（n ≥ 0） ----
inline long long xor_up_to(long long n) {  // 0 ⊕ 1 ⊕ … ⊕ n（n < 0 なら 0）
    if (n < 0) return 0;
    switch (n & 3) {
        case 0: return n;
        case 1: return 1;
        case 2: return n + 1;
        default: return 0;
    }
}
inline long long xor_range(long long l, long long r) {  // l ⊕ … ⊕ (r − 1)（0 ≤ l。空なら 0）
    if (l < 0) math_internal::fail("xor_range", "l ≥ 0 で使う");
    return r <= l ? 0 : xor_up_to(r - 1) ^ xor_up_to(l - 1);
}
inline long long xor_range_closed(long long L, long long R) {  // L ⊕ … ⊕ R（0 ≤ L。ABC121 D）
    if (L < 0) math_internal::fail("xor_range_closed", "L ≥ 0 で使う");
    return R < L ? 0 : xor_up_to(R) ^ xor_up_to(L - 1);
}
inline long long count_bit_up_to(long long n, int k) {  // 0〜n のうち k ビット目が 1 の数（k は 0〜63）
    if (k < 0 || k > 63) math_internal::fail("count_bit_up_to", "k は 0〜63");
    if (n < 0 || k == 63) return 0;
    const unsigned long long m = (unsigned long long)n + 1, p = 1ULL << k;  // 0〜n は m 個。周期 2p のうち後半 p 個が 1
    const unsigned long long full = (m >> (k + 1)) << k, rest = m & (2 * p - 1);
    return (long long)(full + (rest > p ? rest - p : 0));
}
inline Int128 popcount_sum_up_to(long long n) {  // Σ_{x=0}^{n} bit_count(x)（n < 0 なら 0）
    if (n < 0) return 0;
    const unsigned long long m = (unsigned long long)n + 1;  // 0〜m−1 を、m の上の桁から 1 のビットごとに区切って数える
    UInt128 s = 0;
    unsigned long long ones = 0;  // ここより上で 1 のビットの数
    for (unsigned long long u = m; u; ++ones) {
        const int k = 63 - __builtin_clzll(u);
        u ^= 1ULL << k;
        s += ((UInt128)(unsigned long long)k << k >> 1) + ((UInt128)ones << k);  // 下 k ビットが全部の形（k·2^(k−1) 個）＋ 上の 1 の分
    }
    return (Int128)s;
}

// ---- XOR の和（a_i ≥ 0） ----
// （a は vector<long long>・vector<int> など。{1, 2, 3} とも書ける）
template <class V>
Int128 xor_pair_sum(const V& a) {  // Σ_{i<j} (a_i ⊕ a_j)（ABC147 D）  O(63 n)
    long long cnt[63] = {};
    for (const auto& x : a) {
        const long long v = math_internal::to_ll_checked(x, "xor_pair_sum");
        if (v < 0) math_internal::fail("xor_pair_sum", "a_i ≥ 0 で使う");
        for (unsigned long long u = (unsigned long long)v; u; u &= u - 1) cnt[__builtin_ctzll(u)]++;
    }
    const long long n = (long long)a.size();
    Int128 s = 0;
    for (int k = 0; k < 63; ++k) s += (Int128)cnt[k] * (n - cnt[k]) << k;
    return s;
}
inline Int128 xor_pair_sum(std::initializer_list<long long> a) { return xor_pair_sum<std::initializer_list<long long>>(a); }
template <class V>
Int128 subarray_xor_sum(const V& a) {  // すべての区間 [l, r)（長さ 1 以上）の XOR の和（ABC365 E）  O(63 n)
    for (const auto& x : a)
        if (math_internal::to_ll_checked(x, "subarray_xor_sum") < 0) math_internal::fail("subarray_xor_sum", "a_i ≥ 0 で使う");
    Int128 s = 0;
    for (int k = 0; k < 63; ++k) {
        long long c[2] = {1, 0};  // 累積 XOR の k ビット目が 0 / 1 の個数（はじめの 0 を含む）
        int cur = 0;
        for (const auto& x : a) cur ^= (int)(((long long)x >> k) & 1), c[cur]++;
        s += (Int128)c[0] * c[1] << k;
    }
    return s;
}
inline Int128 subarray_xor_sum(std::initializer_list<long long> a) { return subarray_xor_sum<std::initializer_list<long long>>(a); }
#endif  // GG_BIT
