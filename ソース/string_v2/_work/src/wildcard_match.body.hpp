#ifndef GG_WILDCARD_MATCH
#define GG_WILDCARD_MATCH
// ---- ワイルドカードつきの一致・ずらしごとの不一致の数（畳み込み。文字列第二世代） ----
// 中で NTT（mod 998244353）を使う（math_v2 には依存しない）。text の長さ N（2^23 まで）、pat の長さ M。ずらし i は 0〜N-M
// wild（文字列の wildcard_match では省くと '?'）は text・pat のどちらにあっても任意の 1 文字と一致する
namespace str_internal {
struct WcNtt {  // mod 998244353 の NTT（4 基数。結果はビットを逆にした順で、点ごとの積と ntt_inv の組で使う）
    static constexpr unsigned MOD = 998244353;
    static unsigned mul(unsigned a, unsigned b) { return (unsigned)((unsigned long long)a * b % MOD); }
    static unsigned add(unsigned a, unsigned b) { return a + b >= MOD ? a + b - MOD : a + b; }
    static unsigned sub(unsigned a, unsigned b) { return a >= b ? a - b : a + MOD - b; }
    static int size_for(long long need, const char* func) {  // need 以上の最小の 2 冪（2^23 まで）
        if (need > (1LL << 23)) fail(func, "text の長さ N が 2^23 を超える（NTT の長さの上限）");
        int sz = 1;
        while (sz < need) sz <<= 1;
        return sz;
    }
    static unsigned pw(unsigned long long a, unsigned long long e) {
        unsigned long long r = 1;
        a %= MOD;
        for (; e; e >>= 1, a = a * a % MOD)
            if (e & 1) r = r * a % MOD;
        return (unsigned)r;
    }
    struct Table {
        unsigned root[24], iroot[24], rate2[23], irate2[23], rate3[22], irate3[22];
        Table() {
            const int r = 23;  // 998244353 - 1 = 119 × 2^23
            root[r] = pw(3, (MOD - 1) >> r), iroot[r] = pw(root[r], MOD - 2);
            for (int i = r - 1; i >= 0; --i) root[i] = mul(root[i + 1], root[i + 1]), iroot[i] = mul(iroot[i + 1], iroot[i + 1]);
            unsigned prod = 1, iprod = 1;
            for (int i = 0; i + 2 <= r; ++i) {
                rate2[i] = mul(root[i + 2], prod), irate2[i] = mul(iroot[i + 2], iprod);
                prod = mul(prod, iroot[i + 2]), iprod = mul(iprod, root[i + 2]);
            }
            prod = 1, iprod = 1;
            for (int i = 0; i + 3 <= r; ++i) {
                rate3[i] = mul(root[i + 3], prod), irate3[i] = mul(iroot[i + 3], iprod);
                prod = mul(prod, iroot[i + 3]), iprod = mul(iprod, root[i + 3]);
            }
        }
    };
    static const Table& table() {
        static const Table t;
        return t;
    }
    static void ntt(std::vector<unsigned>& a) {  // 長さ 2^h
        const int n = (int)a.size();
        int h = 0;
        while ((1 << h) < n) ++h;
        const Table& t = table();
        for (int len = 0; len < h;) {
            if (h - len == 1) {
                const int p = 1 << (h - len - 1);
                unsigned rot = 1;
                for (int s = 0; s < (1 << len); ++s) {
                    const int off = s << (h - len);
                    for (int i = 0; i < p; ++i) {
                        const unsigned l = a[i + off], r = mul(a[i + off + p], rot);
                        a[i + off] = add(l, r), a[i + off + p] = sub(l, r);
                    }
                    if (s + 1 != (1 << len)) rot = mul(rot, t.rate2[__builtin_ctz(~(unsigned)s)]);
                }
                ++len;
            } else {
                const int p = 1 << (h - len - 2);
                unsigned rot = 1;
                const unsigned imag = t.root[2];
                for (int s = 0; s < (1 << len); ++s) {
                    const unsigned rot2 = mul(rot, rot), rot3 = mul(rot2, rot);
                    const int off = s << (h - len);
                    for (int i = 0; i < p; ++i) {
                        const unsigned a0 = a[i + off], a1 = mul(a[i + off + p], rot), a2 = mul(a[i + off + 2 * p], rot2), a3 = mul(a[i + off + 3 * p], rot3);
                        const unsigned x = mul(sub(a1, a3), imag), s02 = add(a0, a2), d02 = sub(a0, a2), s13 = add(a1, a3);
                        a[i + off] = add(s02, s13);
                        a[i + off + p] = sub(s02, s13);
                        a[i + off + 2 * p] = add(d02, x);
                        a[i + off + 3 * p] = sub(d02, x);
                    }
                    if (s + 1 != (1 << len)) rot = mul(rot, t.rate3[__builtin_ctz(~(unsigned)s)]);
                }
                len += 2;
            }
        }
    }
    static void ntt_inv(std::vector<unsigned>& a) {  // ntt の逆（1/n 倍まで）
        const int n = (int)a.size();
        int h = 0;
        while ((1 << h) < n) ++h;
        const Table& t = table();
        for (int len = h; len;) {
            if (len == 1) {
                const int p = 1 << (h - len);
                unsigned irot = 1;
                for (int s = 0; s < (1 << (len - 1)); ++s) {
                    const int off = s << (h - len + 1);
                    for (int i = 0; i < p; ++i) {
                        const unsigned l = a[i + off], r = a[i + off + p];
                        a[i + off] = add(l, r), a[i + off + p] = mul(sub(l, r), irot);
                    }
                    if (s + 1 != (1 << (len - 1))) irot = mul(irot, t.irate2[__builtin_ctz(~(unsigned)s)]);
                }
                --len;
            } else {
                const int p = 1 << (h - len);
                unsigned irot = 1;
                const unsigned iimag = t.iroot[2];
                for (int s = 0; s < (1 << (len - 2)); ++s) {
                    const unsigned irot2 = mul(irot, irot), irot3 = mul(irot2, irot);
                    const int off = s << (h - len + 2);
                    for (int i = 0; i < p; ++i) {
                        const unsigned a0 = a[i + off], a1 = a[i + off + p], a2 = a[i + off + 2 * p], a3 = a[i + off + 3 * p];
                        const unsigned x = mul(sub(a2, a3), iimag);
                        a[i + off] = add(add(a0, a1), add(a2, a3));
                        a[i + off + p] = mul(add(sub(a0, a1), x), irot);
                        a[i + off + 2 * p] = mul(sub(add(a0, a1), add(a2, a3)), irot2);
                        a[i + off + 3 * p] = mul(sub(sub(a0, a1), x), irot3);
                    }
                    if (s + 1 != (1 << (len - 2))) irot = mul(irot, t.irate3[__builtin_ctz(~(unsigned)s)]);
                }
                len -= 2;
            }
        }
        const unsigned inv_n = pw(n, MOD - 2);
        for (auto& x : a) x = mul(x, inv_n);
    }
};
// 本体（has_wild = false なら wild は使わない）
template <class T, class P>
std::vector<int> mismatch_counts_impl(const T& text0, const P& pat0, bool has_wild, long long wild) {
    using W = WcNtt;
    const auto t = codes_of(text0), p = codes_of(pat0);
    const int n = (int)t.size(), m = (int)p.size();
    if (m > n) return {};
    if (m == 0) return std::vector<int>(n + 1, 0);
    const int sz = W::size_for(n, "mismatch_counts");  // 読むのは添字 m-1〜n-1 だけなので、巡回の長さは N 以上で足りる
    auto is_w = [&](long long x) { return has_wild && x == wild; };
    std::vector<long long> vals, tv;  // 両方に現れる値（片方にしかない値は一致しない）
    for (long long x : p)
        if (!is_w(x)) vals.push_back(x);
    for (long long x : t)
        if (!is_w(x)) tv.push_back(x);
    std::sort(vals.begin(), vals.end()), std::sort(tv.begin(), tv.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    vals.erase(std::remove_if(vals.begin(), vals.end(), [&](long long x) { return !std::binary_search(tv.begin(), tv.end(), x); }), vals.end());
    std::vector<unsigned> acc(sz, 0), fa(sz), fb(sz);
    auto add_product = [&](auto ft, auto fp, int sign) {  // acc += sign × F(text の印) × F(pat を反転した印)
        std::fill(fa.begin(), fa.end(), 0u), std::fill(fb.begin(), fb.end(), 0u);
        for (int i = 0; i < n; ++i) fa[i] = ft(t[i]) ? 1u : 0u;
        for (int j = 0; j < m; ++j) fb[m - 1 - j] = fp(p[j]) ? 1u : 0u;
        W::ntt(fa), W::ntt(fb);
        for (int k = 0; k < sz; ++k) {
            const unsigned v = (unsigned)((unsigned long long)fa[k] * fb[k] % W::MOD);
            acc[k] = sign > 0 ? (acc[k] + v >= W::MOD ? acc[k] + v - W::MOD : acc[k] + v) : (acc[k] >= v ? acc[k] - v : acc[k] + W::MOD - v);
        }
    };
    // 不一致 = （両方とも wild でない組の数）−（同じ値の組の数）
    if (has_wild) add_product([&](long long x) { return !is_w(x); }, [&](long long x) { return !is_w(x); }, +1);
    for (long long c : vals) add_product([&](long long x) { return x == c; }, [&](long long x) { return x == c; }, -1);
    W::ntt_inv(acc);
    std::vector<int> r(n - m + 1);
    for (int i = 0; i <= n - m; ++i) {
        const long long both = has_wild ? 0 : m;  // wild がなければ組はいつも m 個
        long long v = (long long)acc[i + m - 1];
        if (v > (long long)W::MOD / 2) v -= W::MOD;  // 負の値（−一致の数）を戻す
        r[i] = (int)(both + v);
    }
    return r;
}
template <class T, class P>
std::vector<int> wildcard_match_impl(const T& text0, const P& pat0, long long wild) {
    using W = WcNtt;
    const auto t = codes_of(text0), p = codes_of(pat0);
    const int n = (int)t.size(), m = (int)p.size();
    std::vector<int> res;
    if (m > n) return res;
    if (m == 0) {
        for (int i = 0; i <= n; ++i) res.push_back(i);
        return res;
    }
    const int sz = W::size_for(n, "wildcard_match");
    std::vector<char> ok(n - m + 1, 1);
    static std::mt19937_64 rng((unsigned long long)std::chrono::steady_clock::now().time_since_epoch().count() ^ 0x5DEECE66DULL);
    for (int round = 0; round < 2; ++round) {
        const unsigned long long salt = rng();
        auto rv = [&](long long x) -> unsigned long long {  // 値 → [1, MOD) の乱数（wild は 0）
            if (x == wild) return 0;
            unsigned long long z = (unsigned long long)x ^ salt;
            z += 0x9E3779B97F4A7C15ULL, z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL, z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL, z ^= z >> 31;
            return 1 + z % (W::MOD - 1);
        };
        // Σ_j P_j T_{i+j} (P_j - T_{i+j})^2 = Σ P^3 T − 2 P^2 T^2 + P T^3
        std::vector<unsigned> t1(sz, 0), t2(sz, 0), t3(sz, 0), p1(sz, 0), p2(sz, 0), p3(sz, 0);
        for (int i = 0; i < n; ++i) {
            const unsigned long long v = rv(t[i]);
            t1[i] = (unsigned)v, t2[i] = (unsigned)(v * v % W::MOD), t3[i] = (unsigned)((unsigned long long)t2[i] * v % W::MOD);
        }
        for (int j = 0; j < m; ++j) {
            const unsigned long long v = rv(p[j]);
            const int k = m - 1 - j;
            p1[k] = (unsigned)v, p2[k] = (unsigned)(v * v % W::MOD), p3[k] = (unsigned)((unsigned long long)p2[k] * v % W::MOD);
        }
        W::ntt(t1), W::ntt(t2), W::ntt(t3), W::ntt(p1), W::ntt(p2), W::ntt(p3);
        for (int k = 0; k < sz; ++k) {
            const unsigned long long a = (unsigned long long)p3[k] * t1[k] % W::MOD, b = (unsigned long long)p2[k] * t2[k] % W::MOD,
                                     c = (unsigned long long)p1[k] * t3[k] % W::MOD;
            t1[k] = (unsigned)((a + c + 2ULL * (W::MOD - b)) % W::MOD);
        }
        W::ntt_inv(t1);
        for (int i = 0; i <= n - m; ++i)
            if (t1[i + m - 1] != 0) ok[i] = 0;
    }
    for (int i = 0; i <= n - m; ++i)
        if (ok[i]) res.push_back(i);
    return res;
}
}  // namespace str_internal

// ずらし i ごとの不一致の数（r[i] = pat[j] ≠ text[i + j] の j の数。長さ N − M + 1）  O(σ' N log N)（σ' は text と pat の両方に現れる値の数）
template <class T, class P>
std::vector<int> mismatch_counts(const T& text, const P& pat) { return str_internal::mismatch_counts_impl(text, pat, false, 0); }
// wild つき：両方とも wild でない位置で値が違う数（wild は '?' のような文字か、数の列なら値）
template <class T, class P, class W>
std::vector<int> mismatch_counts(const T& text, const P& pat, W wild) { return str_internal::mismatch_counts_impl(text, pat, true, str_internal::code_of(wild)); }
template <class T, class P>
std::vector<int> match_with_k_mismatches(const T& text, const P& pat, int k) {  // 不一致が k 以下のずらし i の一覧（wild なし）
    std::vector<int> r;
    const auto mc = mismatch_counts(text, pat);
    for (int i = 0; i < (int)mc.size(); ++i)
        if (mc[i] <= k) r.push_back(i);
    return r;
}
template <class T, class P, class W>
std::vector<int> match_with_k_mismatches(const T& text, const P& pat, int k, W wild) {  // wild つき
    std::vector<int> r;
    const auto mc = mismatch_counts(text, pat, wild);
    for (int i = 0; i < (int)mc.size(); ++i)
        if (mc[i] <= k) r.push_back(i);
    return r;
}
// wild を任意の 1 文字として一致するずらし i の一覧（乱数を使う。違うのに一致とする確率は 1 か所あたり約 10^-17）  O(N log N)
template <class T, class P, class W>
std::vector<int> wildcard_match(const T& text, const P& pat, W wild) { return str_internal::wildcard_match_impl(text, pat, str_internal::code_of(wild)); }
template <class T, class P>
std::vector<int> wildcard_match(const T& text, const P& pat) {  // wild = '?'（文字列だけ。数の列では wild の値を渡す）
    static_assert(std::is_same<str_internal::elem_t<T>, char>::value && std::is_same<str_internal::elem_t<P>, char>::value,
                  "数の列では wild の値を渡す：wildcard_match(text, pat, wild)");
    return str_internal::wildcard_match_impl(text, pat, '?');
}
#endif  // GG_WILDCARD_MATCH
