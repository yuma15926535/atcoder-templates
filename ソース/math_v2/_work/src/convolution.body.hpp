#ifndef GG_CONVOLUTION
#define GG_CONVOLUTION
// ---- 畳み込み（数学第二世代）  mint の列（NTT・任意の mod）、整数のまま（long long・Int128）、多くの多項式の積、相関、巡回、オンライン ----
namespace math_internal {
template <class mint>
struct NttTable {  // mod ごとの 1 の冪根の表（ModInt は 1 回、DynModInt は mod が変わったら作り直す）
    int mod = 0, rank2 = 0;
    std::vector<mint> root, iroot, rate2, irate2, rate3, irate3;
};
inline unsigned primitive_root_u32(unsigned p) {  // 素数 p の原始根（いちばん小さいもの）
    if (p == 2) return 1;
    std::vector<unsigned> qs;
    unsigned x = p - 1;
    for (unsigned q = 2; (unsigned long long)q * q <= x; ++q)
        if (x % q == 0) {
            qs.push_back(q);
            while (x % q == 0) x /= q;
        }
    if (x > 1) qs.push_back(x);
    for (unsigned g = 2;; ++g) {
        bool ok = true;
        for (unsigned q : qs) {
            unsigned long long r = 1, b = g, e = (p - 1) / q;
            for (; e; e >>= 1, b = b * b % p)
                if (e & 1) r = r * b % p;
            if (r == 1) { ok = false; break; }
        }
        if (ok) return g;
    }
}
template <class mint>
const NttTable<mint>& ntt_table() {
    static NttTable<mint> t;
    if (t.mod == mint::mod()) return t;
    const unsigned p = (unsigned)mint::mod();
    if (!is_prime_u32(p)) fail("ntt", "mod が素数でない（convolve なら自動で 3 つの素数を使う）");
    t = NttTable<mint>();
    t.mod = (int)p;
    t.rank2 = __builtin_ctz(p - 1);
    const int r = t.rank2;
    t.root.assign(r + 1, mint(1)), t.iroot.assign(r + 1, mint(1));
    t.rate2.assign(std::max(0, r - 1), mint(1)), t.irate2.assign(std::max(0, r - 1), mint(1));
    t.rate3.assign(std::max(0, r - 2), mint(1)), t.irate3.assign(std::max(0, r - 2), mint(1));
    const mint g = mint((long long)primitive_root_u32(p));
    t.root[r] = g.pow((long long)((p - 1) >> r));
    t.iroot[r] = t.root[r].inv();
    for (int i = r - 1; i >= 0; --i) t.root[i] = t.root[i + 1] * t.root[i + 1], t.iroot[i] = t.iroot[i + 1] * t.iroot[i + 1];
    mint prod = 1, iprod = 1;
    for (int i = 0; i + 2 <= r; ++i) {
        t.rate2[i] = t.root[i + 2] * prod, t.irate2[i] = t.iroot[i + 2] * iprod;
        prod *= t.iroot[i + 2], iprod *= t.root[i + 2];
    }
    prod = 1, iprod = 1;
    for (int i = 0; i + 3 <= r; ++i) {
        t.rate3[i] = t.root[i + 3] * prod, t.irate3[i] = t.iroot[i + 3] * iprod;
        prod *= t.iroot[i + 3], iprod *= t.root[i + 3];
    }
    return t;
}
template <class mint>
bool ntt_usable(size_t len) {  // mod が素数で、長さ len（2 の冪）の NTT ができるか
    const unsigned p = (unsigned)mint::mod();
    if (p < 3 || !is_prime_u32(p)) return false;
    return (size_t)1 << __builtin_ctz(p - 1) >= len;
}
inline int ceil_log2(size_t n) {
    int h = 0;
    while (((size_t)1 << h) < n) ++h;
    return h;
}
}  // namespace math_internal

// 長さ 2^h の NTT（結果はビットを逆にした順。intt と組で、点ごとの積をとるのに使う。mod は素数で 2^h | mod − 1）
template <class mint>
void ntt(std::vector<mint>& a) {
    const int n = (int)a.size(), h = math_internal::ceil_log2(a.size());
    if (n == 0) return;
    if ((1 << h) != n) math_internal::fail("ntt", "長さは 2 の冪");
    const auto& t = math_internal::ntt_table<mint>();
    if (h > t.rank2) math_internal::fail("ntt", "この mod では長さが大きすぎる");
    int len = 0;
    while (len < h) {
        if (h - len == 1) {
            const int p = 1 << (h - len - 1);
            mint rot = 1;
            for (int s = 0; s < (1 << len); ++s) {
                const int off = s << (h - len);
                for (int i = 0; i < p; ++i) {
                    const mint l = a[i + off], r = a[i + off + p] * rot;
                    a[i + off] = l + r, a[i + off + p] = l - r;
                }
                if (s + 1 != (1 << len)) rot *= t.rate2[__builtin_ctz(~(unsigned)s)];
            }
            ++len;
        } else {
            const int p = 1 << (h - len - 2);
            mint rot = 1;
            const mint imag = t.root[2];
            for (int s = 0; s < (1 << len); ++s) {
                const mint rot2 = rot * rot, rot3 = rot2 * rot;
                const int off = s << (h - len);
                for (int i = 0; i < p; ++i) {
                    const mint a0 = a[i + off], a1 = a[i + off + p] * rot, a2 = a[i + off + 2 * p] * rot2, a3 = a[i + off + 3 * p] * rot3;
                    const mint x = (a1 - a3) * imag;
                    a[i + off] = a0 + a2 + a1 + a3;
                    a[i + off + p] = a0 + a2 - a1 - a3;
                    a[i + off + 2 * p] = a0 - a2 + x;
                    a[i + off + 3 * p] = a0 - a2 - x;
                }
                if (s + 1 != (1 << len)) rot *= t.rate3[__builtin_ctz(~(unsigned)s)];
            }
            len += 2;
        }
    }
}
template <class mint>
void intt(std::vector<mint>& a, bool divide = true) {  // ntt の逆（divide なら 1/n 倍まで）
    const int n = (int)a.size(), h = math_internal::ceil_log2(a.size());
    if (n == 0) return;
    if ((1 << h) != n) math_internal::fail("intt", "長さは 2 の冪");
    const auto& t = math_internal::ntt_table<mint>();
    if (h > t.rank2) math_internal::fail("intt", "この mod では長さが大きすぎる");
    int len = h;
    while (len) {
        if (len == 1) {
            const int p = 1 << (h - len);
            mint irot = 1;
            for (int s = 0; s < (1 << (len - 1)); ++s) {
                const int off = s << (h - len + 1);
                for (int i = 0; i < p; ++i) {
                    const mint l = a[i + off], r = a[i + off + p];
                    a[i + off] = l + r, a[i + off + p] = (l - r) * irot;
                }
                if (s + 1 != (1 << (len - 1))) irot *= t.irate2[__builtin_ctz(~(unsigned)s)];
            }
            --len;
        } else {
            const int p = 1 << (h - len);
            mint irot = 1;
            const mint iimag = t.iroot[2];
            for (int s = 0; s < (1 << (len - 2)); ++s) {
                const mint irot2 = irot * irot, irot3 = irot2 * irot;
                const int off = s << (h - len + 2);
                for (int i = 0; i < p; ++i) {
                    const mint a0 = a[i + off], a1 = a[i + off + p], a2 = a[i + off + 2 * p], a3 = a[i + off + 3 * p];
                    const mint x = (a2 - a3) * iimag;
                    a[i + off] = a0 + a1 + a2 + a3;
                    a[i + off + p] = (a0 - a1 + x) * irot;
                    a[i + off + 2 * p] = (a0 + a1 - a2 - a3) * irot2;
                    a[i + off + 3 * p] = (a0 - a1 - x) * irot3;
                }
                if (s + 1 != (1 << (len - 2))) irot *= t.irate3[__builtin_ctz(~(unsigned)s)];
            }
            len -= 2;
        }
    }
    if (divide) {
        const mint inv_n = mint(n).inv();
        for (mint& x : a) x *= inv_n;
    }
}

namespace math_internal {
template <class mint>
std::vector<mint> conv_naive(const std::vector<mint>& a, const std::vector<mint>& b) {
    std::vector<mint> c(a.size() + b.size() - 1);
    if (a.size() < b.size()) {
        for (size_t j = 0; j < b.size(); ++j)
            for (size_t i = 0; i < a.size(); ++i) c[i + j] += a[i] * b[j];
    } else {
        for (size_t i = 0; i < a.size(); ++i)
            for (size_t j = 0; j < b.size(); ++j) c[i + j] += a[i] * b[j];
    }
    return c;
}
template <class mint>
std::vector<mint> conv_ntt(std::vector<mint> a, std::vector<mint> b) {  // mod が NTT に使える素数のとき
    const size_t z = a.size() + b.size() - 1, n = (size_t)1 << ceil_log2(z);
    const bool same = a == b;
    const mint inv_n = mint((long long)n).inv();  // intt の 1/n は点ごとの積に混ぜる
    a.resize(n), ntt(a);
    if (same) {
        for (size_t i = 0; i < n; ++i) a[i] *= a[i] * inv_n;
    } else {
        b.resize(n), ntt(b);
        for (size_t i = 0; i < n; ++i) a[i] *= b[i] * inv_n;
    }
    intt(a, false);
    a.resize(z);
    return a;
}
// 3〜5 個の NTT 素数（2^23 以上の 2 の冪で割れる）での畳み込みの余り
constexpr unsigned NTT_P0 = 754974721, NTT_P1 = 167772161, NTT_P2 = 469762049, NTT_P3 = 998244353, NTT_P4 = 2013265921;
template <unsigned P>
std::vector<unsigned> conv_residue(const std::vector<long long>& a, const std::vector<long long>& b) {
    using M = ModInt<P>;
    std::vector<M> x(a.size()), y(b.size());
    for (size_t i = 0; i < a.size(); ++i) x[i] = M(a[i]);
    for (size_t i = 0; i < b.size(); ++i) y[i] = M(b[i]);
    const std::vector<M> z = (std::min(a.size(), b.size()) <= 60) ? conv_naive(x, y) : conv_ntt(std::move(x), std::move(y));
    std::vector<unsigned> r(z.size());
    for (size_t i = 0; i < z.size(); ++i) r[i] = (unsigned)z[i].val();
    return r;
}
// 素数 ps[0..k) での余り rs から、−M/2 < x < M/2（M = Π ps）の x を 2^128 で割った余りとして戻す（つり合った混合基数のガーナー）
struct GarnerBalanced {
    int k;
    long long p[5], inv[5], pm[5][5];  // inv[j] = (p_0 … p_{j−1})^{−1} mod p_j、pm[i][j] = p_0 … p_{i−1} mod p_j
    UInt128 base[5];                   // p_0 … p_{j−1}（2^128 で割った余り）
    GarnerBalanced(const unsigned* ps, int k_) : k(k_) {
        for (int j = 0; j < k; ++j) p[j] = ps[j];
        for (int j = 0; j < k; ++j) {
            long long pr = 1;
            for (int i = 0; i <= j; ++i) pm[i][j] = pr, pr = pr * p[i] % p[j];
            inv[j] = mod_inv(pm[j][j], p[j]);
        }
        UInt128 b = 1;
        for (int j = 0; j < k; ++j) base[j] = b, b *= (UInt128)p[j];
    }
    Int128 operator()(const unsigned* rs) const {
        long long t[5];
        UInt128 x = 0;
        for (int j = 0; j < k; ++j) {
            long long cur = 0;  // すでに決めた部分 Σ_{i<j} t_i (p_0 … p_{i−1}) を p_j で割った余り
            for (int i = 0; i < j; ++i) cur = (cur + (t[i] % p[j] + p[j]) % p[j] * pm[i][j]) % p[j];
            long long v = ((long long)rs[j] - cur) % p[j];
            if (v < 0) v += p[j];
            v = v * inv[j] % p[j];
            if (v > p[j] / 2) v -= p[j];  // (−p/2, p/2] に
            t[j] = v;
            x += (UInt128)(Int128)v * base[j];
        }
        return (Int128)x;
    }
};
}  // namespace math_internal

// mint の列の畳み込み c[k] = Σ_{i+j=k} a[i] b[j]。mod が 998244353 などの NTT 素数なら NTT、それ以外（10^9+7・合成数）は
// 3 つの素数で計算してガーナーで戻す（自動）。短い方が 60 以下なら愚直  O((n + m) log(n + m))
template <class mint>
std::vector<mint> convolve(const std::vector<mint>& a, const std::vector<mint>& b) {
    if (a.empty() || b.empty()) return {};
    if (std::min(a.size(), b.size()) <= 60) return math_internal::conv_naive(a, b);
    const size_t z = a.size() + b.size() - 1;
    if (math_internal::ntt_usable<mint>((size_t)1 << math_internal::ceil_log2(z))) return math_internal::conv_ntt(a, b);
    if (z > (1u << 24)) math_internal::fail("convolve", "長さが 2^24 を超える");
    using namespace math_internal;
    std::vector<long long> x(a.size()), y(b.size());
    for (size_t i = 0; i < a.size(); ++i) x[i] = a[i].val();
    for (size_t i = 0; i < b.size(); ++i) y[i] = b[i].val();
    const std::vector<unsigned> r0 = conv_residue<NTT_P0>(x, y), r1 = conv_residue<NTT_P1>(x, y), r2 = conv_residue<NTT_P2>(x, y);
    // 値は 0 以上 min(n, m)·(mod − 1)^2 < P0 P1 P2 なので、ふつうのガーナーで mod の値に
    const long long m = mint::mod();
    const long long i01 = mod_inv(NTT_P0, NTT_P1), i012 = mod_inv((long long)NTT_P0 * NTT_P1 % NTT_P2, NTT_P2);
    const long long p01_m = (long long)NTT_P0 * NTT_P1 % m;
    std::vector<mint> c(z);
    for (size_t i = 0; i < z; ++i) {
        const long long t0 = r0[i];
        const long long t1 = (long long)(r1[i] + NTT_P1 - t0 % NTT_P1) % NTT_P1 * i01 % NTT_P1;
        const long long cur2 = (t0 + (long long)NTT_P0 % NTT_P2 * t1) % NTT_P2;
        const long long t2 = ((long long)r2[i] - cur2 % NTT_P2 + NTT_P2) % NTT_P2 * i012 % NTT_P2;
        c[i] = mint((t0 + (Int128)NTT_P0 * t1 + (Int128)p01_m * t2) % m);
    }
    return c;
}
namespace math_internal {
// 整数のまま畳み込む（結果の絶対値が 2^126 未満と分かっているときだけ。足りるだけの素数を使う）
inline std::vector<Int128> convolve_exact(const std::vector<long long>& a, const std::vector<long long>& b, const char* func) {
    if (a.empty() || b.empty()) return {};
    const size_t z = a.size() + b.size() - 1;
    Int128 ma = 0, mb = 0;
    for (long long v : a) ma = std::max<Int128>(ma, v < 0 ? -(Int128)v : (Int128)v);
    for (long long v : b) mb = std::max<Int128>(mb, v < 0 ? -(Int128)v : (Int128)v);
    // 係数の大きさの上限 ≤ min(n, m)·max|a|·max|b|（long double で見積もる）
    const long double bound = (long double)std::min(a.size(), b.size()) * (long double)ma * (long double)mb;
    if (bound >= 8.0e37L) fail(func, "結果が Int128 に収まらないかもしれない");
    if (std::min(a.size(), b.size()) <= 60) {  // 愚直（128 ビットで。上の見積もりで溢れない）
        std::vector<Int128> c(z, 0);
        for (size_t i = 0; i < a.size(); ++i)
            for (size_t j = 0; j < b.size(); ++j) c[i + j] += (Int128)a[i] * b[j];
        return c;
    }
    const int k = bound < 4.0e24L ? 3 : 5;  // 3 つの素数の積 ≈ 5.9×10^25、5 つ ≈ 1.2×10^44（2^126 ≈ 8.5×10^37 まで）
    if (z > (1u << (k == 3 ? 24 : 23))) fail(func, "長さが大きすぎる");
    const unsigned ps[5] = {NTT_P0, NTT_P1, NTT_P2, NTT_P3, NTT_P4};
    std::vector<unsigned> rs[5];
    rs[0] = conv_residue<NTT_P0>(a, b), rs[1] = conv_residue<NTT_P1>(a, b), rs[2] = conv_residue<NTT_P2>(a, b);
    if (k == 5) rs[3] = conv_residue<NTT_P3>(a, b), rs[4] = conv_residue<NTT_P4>(a, b);
    std::vector<Int128> c(z);
    if (k == 3) {  // 3 つの素数：ふつうのガーナーで [0, M) にしてから M/2 を超えたら M を引く
        const long long P0 = NTT_P0, P1 = NTT_P1, P2 = NTT_P2;
        const long long i01 = mod_inv(P0, P1), i012 = mod_inv(P0 * P1 % P2, P2), p0_2 = P0 % P2;
        const Int128 M = (Int128)P0 * P1 * P2, P01 = (Int128)P0 * P1;
        for (size_t i = 0; i < z; ++i) {
            const long long t0 = rs[0][i];
            const long long t1 = (rs[1][i] - t0 % P1 + P1) % P1 * i01 % P1;
            const long long t2 = ((long long)rs[2][i] - (t0 + p0_2 * t1) % P2 + P2) % P2 * i012 % P2;
            Int128 x = t0 + (Int128)P0 * t1 + P01 * t2;
            if (x > M / 2) x -= M;
            c[i] = x;
        }
        return c;
    }
    const GarnerBalanced g(ps, k);
    for (size_t i = 0; i < z; ++i) {
        unsigned r[5];
        for (int j = 0; j < k; ++j) r[j] = rs[j][i];
        c[i] = g(r);
    }
    return c;
}
}  // namespace math_internal
// 整数の列の畳み込み（厳密。負の数も。結果が long long に収まらなければ止める）
inline std::vector<long long> convolve_ll(const std::vector<long long>& a, const std::vector<long long>& b) {
    const std::vector<Int128> c = math_internal::convolve_exact(a, b, "convolve_ll");
    std::vector<long long> r(c.size());
    for (size_t i = 0; i < c.size(); ++i) {
        if (c[i] > (Int128)9223372036854775807LL || c[i] < -(Int128)9223372036854775807LL - 1) math_internal::fail("convolve_ll", "結果が long long に収まらない（convolve_i128 を使う）");
        r[i] = (long long)c[i];
    }
    return r;
}
// 整数の列の畳み込みを Int128 で（結果の絶対値は 2^126 程度まで）
inline std::vector<Int128> convolve_i128(const std::vector<long long>& a, const std::vector<long long>& b) { return math_internal::convolve_exact(a, b, "convolve_i128"); }
// long long の列を mod m（1 ≤ m < 2^31）で畳み込む
inline std::vector<long long> convolve_mod(const std::vector<long long>& a, const std::vector<long long>& b, long long m) {
    if (m < 1 || m >= (1LL << 31)) math_internal::fail("convolve_mod", "1 ≤ m < 2^31 で使う");
    using M = DynModIntT<-1>;  // 利用者の DynModInt とは別の mod を持つ
    const long long keep = M::mod();
    M::set_mod(m);
    std::vector<M> x(a.size()), y(b.size());
    for (size_t i = 0; i < a.size(); ++i) x[i] = M(a[i]);
    for (size_t i = 0; i < b.size(); ++i) y[i] = M(b[i]);
    const std::vector<M> z = convolve(x, y);
    std::vector<long long> r(z.size());
    for (size_t i = 0; i < z.size(); ++i) r[i] = z[i].val();
    M::set_mod(keep);
    return r;
}
// 多くの多項式の積（短いものから 2 つずつ掛ける）  O(N log^2 N)
template <class mint>
std::vector<mint> convolve_many(std::vector<std::vector<mint>> ps) {
    if (ps.empty()) return {mint(1)};
    for (const auto& p : ps)
        if (p.empty()) return {};
    while (ps.size() > 1) {  // 短い順に並べ、隣どうしを掛けて半分にする（短いものから掛けていく）
        std::sort(ps.begin(), ps.end(), [](const std::vector<mint>& x, const std::vector<mint>& y) { return x.size() < y.size(); });
        std::vector<std::vector<mint>> nx;
        nx.reserve(ps.size() / 2 + 1);
        for (size_t i = 0; i + 1 < ps.size(); i += 2) nx.push_back(convolve(ps[i], ps[i + 1]));
        if (ps.size() % 2) nx.push_back(std::move(ps.back()));
        ps.swap(nx);
    }
    return ps[0];
}
// 相関 c[k] = Σ_{i<|b|} a[k+i]·b[i]（k = 0〜|a|−|b|。a の中で b をずらして重ねた積の和。|a| ≥ |b| ≥ 1）
template <class mint>
std::vector<mint> correlate(const std::vector<mint>& a, const std::vector<mint>& b) {
    if (b.empty() || a.size() < b.size()) math_internal::fail("correlate", "|a| ≥ |b| ≥ 1 で使う");
    const std::vector<mint> rb(b.rbegin(), b.rend());
    const std::vector<mint> c = convolve(a, rb);
    return std::vector<mint>(c.begin() + (long)(b.size() - 1), c.begin() + (long)a.size());
}
// 巡回畳み込み c[k] = Σ_{(i+j) mod n = k} a[i] b[j]（長さ n。a, b の添字は mod n で見る）
template <class mint>
std::vector<mint> convolve_cyclic(const std::vector<mint>& a, const std::vector<mint>& b, size_t n) {
    if (n == 0) math_internal::fail("convolve_cyclic", "n ≥ 1 で使う");
    std::vector<mint> c(n, mint(0));
    const std::vector<mint> d = convolve(a, b);
    for (size_t k = 0; k < d.size(); ++k) c[k % n] += d[k];
    return c;
}
// オンラインの畳み込み：a_n, b_n を 1 つずつ渡すと、その時点で c_n = Σ_{i+j=n} a_i b_j を返す（後の値を知らずに）
//   f_{n+1} = Σ_{i=0}^{n} f_i g_{n−i} のように、答えを次の入力に使う DP に。全体で O(N log^2 N)
template <class mint>
struct RelaxedConvolution {
    mint push(const mint& x, const mint& y) {  // a_n = x、b_n = y を加えて c_n を返す（n = これまでに渡した数）
        const size_t n = a_.size();
        a_.push_back(x), b_.push_back(y);
        if (c_.size() < 2 * n + 2) c_.resize(2 * n + 2, mint(0));
        c_[n] += n == 0 ? a_[0] * b_[0] : a_[n] * b_[0] + a_[0] * b_[n];
        // 添字が 1 以上の部分を大きさ s の正方形に分けて、使う値がそろった時点で掛ける（[s, 2s) × [m s, (m+1) s) と、m ≥ 2 ならその逆）
        for (size_t s = 1; 2 * s <= n + 1 && (n + 1) % s == 0; s *= 2) {
            const size_t m = (n + 1) / s - 1;
            add_block(a_, s, b_, m * s, s, (m + 1) * s);
            if (m >= 2) add_block(b_, s, a_, m * s, s, (m + 1) * s);
        }
        return c_[n];
    }
    size_t size() const { return a_.size(); }  // これまでに渡した数
    mint get(size_t n) const {                // 確定した c_n（n < size()）
        if (n >= a_.size()) math_internal::fail("RelaxedConvolution::get", "n < size() で使う");
        return c_[n];
    }

  private:
    std::vector<mint> a_, b_, c_;
    void add_block(const std::vector<mint>& p, size_t pl, const std::vector<mint>& q, size_t ql, size_t s, size_t out) {  // p[pl, pl+s) × q[ql, ql+s) を c[out + …] に
        if (s <= 64) {
            for (size_t i = 0; i < s; ++i)
                for (size_t j = 0; j < s; ++j) c_[out + i + j] += p[pl + i] * q[ql + j];
            return;
        }
        const std::vector<mint> u(p.begin() + (long)pl, p.begin() + (long)(pl + s)), v(q.begin() + (long)ql, q.begin() + (long)(ql + s));
        const std::vector<mint> w = convolve(u, v);
        for (size_t k = 0; k < w.size(); ++k) c_[out + k] += w[k];
    }
};
#endif  // GG_CONVOLUTION
