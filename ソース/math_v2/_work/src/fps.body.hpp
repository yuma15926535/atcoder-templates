#ifndef GG_FPS
#define GG_FPS
// ---- 形式的冪級数（数学第二世代）  Fps<mint>：逆数・log・exp・累乗・平方根・割り算・テイラーシフト・P/Q の K 次の係数・多点評価・補間・べき乗和 ----
// Fps<mint> は vector<mint>（低い次数から）に + − × を足したもの。n を渡す関数は「最初の n 項（x^0〜x^{n−1}）」を返す
template <class mint>
struct Fps : std::vector<mint> {
    using std::vector<mint>::vector;
    Fps() = default;
    Fps(const std::vector<mint>& v) : std::vector<mint>(v) {}
    Fps(std::vector<mint>&& v) : std::vector<mint>(std::move(v)) {}
    Fps pre(size_t n) const {  // 最初の n 項（足りなければ 0 を足す）
        Fps r(this->begin(), this->begin() + (long)std::min(n, this->size()));
        r.resize(n);
        return r;
    }
    mint eval(const mint& x) const {  // f(x)
        mint r = 0;
        for (size_t i = this->size(); i > 0; --i) r = r * x + (*this)[i - 1];
        return r;
    }
    Fps& operator+=(const Fps& g) {
        if (g.size() > this->size()) this->resize(g.size());
        for (size_t i = 0; i < g.size(); ++i) (*this)[i] += g[i];
        return *this;
    }
    Fps& operator-=(const Fps& g) {
        if (g.size() > this->size()) this->resize(g.size());
        for (size_t i = 0; i < g.size(); ++i) (*this)[i] -= g[i];
        return *this;
    }
    Fps& operator*=(const Fps& g) { return *this = Fps(convolve<mint>(*this, g)); }
    Fps& operator*=(const mint& s) {
        for (mint& x : *this) x *= s;
        return *this;
    }
    friend Fps operator+(Fps f, const Fps& g) { return f += g; }
    friend Fps operator-(Fps f, const Fps& g) { return f -= g; }
    friend Fps operator*(const Fps& f, const Fps& g) { return Fps(convolve<mint>(f, g)); }
    friend Fps operator*(Fps f, const mint& s) { return f *= s; }
    friend Fps operator*(const mint& s, Fps f) { return f *= s; }
    Fps operator-() const {
        Fps r(*this);
        for (mint& x : r) x = -x;
        return r;
    }
    friend std::ostream& operator<<(std::ostream& os, const Fps& f) {  // 係数を空白区切りで（そのまま答えの出力に使える）
        for (size_t i = 0; i < f.size(); ++i) os << (i ? " " : "") << f[i];
        return os;
    }
};

namespace math_internal {
template <class mint>
void fps_need_prime(long long n, const char* func) {  // 1〜n の逆数を使う関数（mod は n より大きい素数）
    if (!is_prime_u32((unsigned)mint::mod()) || n >= mint::mod()) fail(func, "mod が項の数より大きい素数のときに使う");
}
template <class mint>
std::vector<mint> fps_inverses(int n) {  // 1/0（使わない）, 1/1, …, 1/n
    std::vector<mint> iv((size_t)n + 1, mint(1));
    const long long p = mint::mod();
    for (int i = 2; i <= n; ++i) iv[(size_t)i] = -iv[(size_t)(p % i)] * mint(p / i);
    return iv;
}
}  // namespace math_internal

template <class mint>
Fps<mint> fps_derivative(const std::vector<mint>& f) {  // f'
    Fps<mint> r(f.size() > 1 ? f.size() - 1 : 0);
    for (size_t i = 1; i < f.size(); ++i) r[i - 1] = f[i] * mint((long long)i);
    return r;
}
template <class mint>
Fps<mint> fps_integral(const std::vector<mint>& f) {  // ∫f（定数項 0）
    math_internal::fps_need_prime<mint>((long long)f.size(), "fps_integral");
    const std::vector<mint> iv = math_internal::fps_inverses<mint>((int)f.size());
    Fps<mint> r(f.size() + 1);
    for (size_t i = 0; i < f.size(); ++i) r[i + 1] = f[i] * iv[i + 1];
    return r;
}
// 1/f の最初の n 項（f[0] ≠ 0。f[0] が mod と互いに素なら合成数の mod でも）  O(n log n)
template <class mint>
Fps<mint> fps_inv(const std::vector<mint>& f, int n) {
    if (n < 0) math_internal::fail("fps_inv", "n ≥ 0 で使う");
    if (n == 0) return {};
    if (f.empty() || f[0] == mint(0)) math_internal::fail("fps_inv", "f[0] が 0");
    Fps<mint> g{f[0].inv()};
    const bool fast = math_internal::ntt_usable<mint>((size_t)1 << math_internal::ceil_log2((size_t)2 * n));
    for (int m = 1; m < n; m *= 2) {
        if (fast) {  // NTT 5 回：g ← g − g (f g − 1)
            std::vector<mint> F((size_t)2 * m), G((size_t)2 * m);
            for (int i = 0; i < 2 * m && i < (int)f.size(); ++i) F[i] = f[i];
            for (int i = 0; i < m; ++i) G[i] = g[i];
            ntt(F), ntt(G);
            for (int i = 0; i < 2 * m; ++i) F[i] *= G[i];
            intt(F);
            for (int i = 0; i < m; ++i) F[i] = 0;  // f g の最初の m 項は 1, 0, …（使わない）
            ntt(F);
            for (int i = 0; i < 2 * m; ++i) F[i] *= G[i];
            intt(F);
            g.resize((size_t)2 * m);
            for (int i = m; i < 2 * m; ++i) g[i] = -F[i];
        } else {  // どの mod でも：g ← g (2 − f g)
            Fps<mint> h = Fps<mint>(f).pre((size_t)2 * m) * g;
            h.resize((size_t)2 * m);
            for (mint& x : h) x = -x;
            h[0] += 2;
            g = (g * h).pre((size_t)2 * m);
        }
    }
    return g.pre((size_t)n);
}
// log f の最初の n 項（f[0] = 1。mod は n より大きい素数）  O(n log n)
template <class mint>
Fps<mint> fps_log(const std::vector<mint>& f, int n) {
    if (n < 0) math_internal::fail("fps_log", "n ≥ 0 で使う");
    if (n == 0) return {};
    if (f.empty() || f[0] != mint(1)) math_internal::fail("fps_log", "f[0] = 1 で使う");
    Fps<mint> d = fps_derivative(Fps<mint>(f).pre((size_t)n)) * fps_inv(f, n);
    d.resize((size_t)n - 1);
    return fps_integral(d).pre((size_t)n);
}
// exp f の最初の n 項（f[0] = 0。mod は n より大きい素数）  O(n log n)
template <class mint>
Fps<mint> fps_exp(const std::vector<mint>& f, int n) {
    if (n < 0) math_internal::fail("fps_exp", "n ≥ 0 で使う");
    if (n == 0) return {};
    if (!f.empty() && f[0] != mint(0)) math_internal::fail("fps_exp", "f[0] = 0 で使う");
    math_internal::fps_need_prime<mint>(n, "fps_exp");
    Fps<mint> g{mint(1)};
    for (int m = 1; m < n; m *= 2) {  // g ← g (1 − log g + f)（最後の段は n 項まで：mod が小さい素数でも 1/n までしか使わない）
        const int len = std::min(2 * m, n);
        Fps<mint> h = fps_log(g, len);
        for (int i = 0; i < len; ++i) h[i] = (i < (int)f.size() ? f[i] : mint(0)) - h[i];
        h[0] += 1;
        g = (g * h).pre((size_t)len);
    }
    return g.pre((size_t)n);
}
// f^k の最初の n 項（k ≥ 0。k は 10^18 でもよい。f[0] = 0 でもよい。mod は n より大きい素数）  O(n log n)
template <class mint>
Fps<mint> fps_pow(const std::vector<mint>& f, long long k, int n) {
    if (n < 0 || k < 0) math_internal::fail("fps_pow", "n, k ≥ 0 で使う");
    Fps<mint> r((size_t)n);
    if (n == 0) return r;
    if (k == 0) {
        r[0] = 1;
        return r;
    }
    size_t d = 0;
    while (d < f.size() && f[d] == mint(0)) ++d;
    if (d == f.size() || (Int128)d * k >= n) return r;  // 全部 0
    const size_t shift = (size_t)d * (size_t)k, m = (size_t)n - shift;
    const mint c = f[d], ic = c.inv();
    Fps<mint> g(f.begin() + (long)d, f.begin() + (long)std::min(f.size(), d + m));
    for (mint& x : g) x *= ic;  // g = f / (c x^d)、g[0] = 1
    Fps<mint> lg = fps_log(g, (int)m);
    for (mint& x : lg) x *= mint(k);
    const Fps<mint> e = fps_exp(lg, (int)m);
    const mint ck = c.pow(k);
    for (size_t i = 0; i < m; ++i) r[shift + i] = e[i] * ck;
    return r;
}
namespace math_internal {
template <class mint>
Maybe<mint> fps_mod_sqrt(const mint& a) {  // mod は奇素数。トネリ・シャンクス
    const long long p = mint::mod();
    if (a == mint(0)) return {true, mint(0)};
    if (p == 2) return {true, a};
    if (a.pow((p - 1) / 2) != mint(1)) return {false, mint(0)};
    long long q = p - 1;
    int s = 0;
    while (q % 2 == 0) q /= 2, ++s;
    mint z = 2;
    while (z.pow((p - 1) / 2) == mint(1)) z += 1;
    mint c = z.pow(q), x = a.pow((q + 1) / 2), t = a.pow(q);
    int m = s;
    while (t != mint(1)) {
        int i = 0;
        mint tt = t;
        while (tt != mint(1)) tt *= tt, ++i;
        mint b = c;
        for (int j = 0; j < m - i - 1; ++j) b *= b;
        x *= b, c = b * b, t *= c, m = i;
    }
    return {true, x};
}
}  // namespace math_internal
// √f の最初の n 項（なければ ok = false。f[0] = 0 でもよい。mod は奇素数）  O(n log n)
template <class mint>
Maybe<Fps<mint>> fps_sqrt(const std::vector<mint>& f, int n) {
    if (n < 0) math_internal::fail("fps_sqrt", "n ≥ 0 で使う");
    if (mint::mod() == 2 || !math_internal::is_prime_u32((unsigned)mint::mod())) math_internal::fail("fps_sqrt", "mod が奇素数のときに使う");
    Fps<mint> r((size_t)n);
    if (n == 0) return {true, r};
    size_t d = 0;
    while (d < f.size() && f[d] == mint(0)) ++d;
    if (d == f.size()) return {true, r};  // f = 0
    if (d % 2 == 1) return {false, Fps<mint>()};
    const Maybe<mint> s0 = math_internal::fps_mod_sqrt(f[d]);
    if (!s0.ok) return {false, Fps<mint>()};
    if (d / 2 >= (size_t)n) return {true, r};  // 最初の n 項は 0
    const size_t m = (size_t)n - d / 2;
    const Fps<mint> g(f.begin() + (long)d, f.begin() + (long)std::min(f.size(), d + m));
    Fps<mint> h{s0.value};
    const mint inv2 = mint(2).inv();
    for (size_t k = 1; k < m; k *= 2) {  // h ← (h + g/h) / 2
        Fps<mint> t = (g.pre(2 * k) * fps_inv(h, (int)(2 * k))).pre(2 * k);
        h.resize(2 * k);
        for (size_t i = 0; i < 2 * k; ++i) h[i] = (h[i] + t[i]) * inv2;
    }
    for (size_t i = 0; i < m; ++i) r[d / 2 + i] = h[i];
    return {true, r};
}
template <class mint>
struct FpsDivMod {  // f = q g + r（deg r < deg g）
    Fps<mint> q, r;
};
// 多項式の割り算（g の最高次の係数で割る。末尾の 0 は除く）  O(n log n)
template <class mint>
FpsDivMod<mint> fps_divmod(std::vector<mint> f, std::vector<mint> g) {
    while (!f.empty() && f.back() == mint(0)) f.pop_back();
    while (!g.empty() && g.back() == mint(0)) g.pop_back();
    if (g.empty()) math_internal::fail("fps_divmod", "0 で割った");
    FpsDivMod<mint> res;
    if (f.size() < g.size()) {
        res.r = f;
        return res;
    }
    const size_t k = f.size() - g.size() + 1;
    Fps<mint> rf(f.rbegin(), f.rend()), rg(g.rbegin(), g.rend());
    Fps<mint> q = (rf.pre(k) * fps_inv(rg, (int)k)).pre(k);
    std::reverse(q.begin(), q.end());
    Fps<mint> r = Fps<mint>(f) - Fps<mint>(g) * q;
    r.resize(g.size() - 1);
    while (!r.empty() && r.back() == mint(0)) r.pop_back();
    while (!q.empty() && q.back() == mint(0)) q.pop_back();
    res.q = q, res.r = r;
    return res;
}
// f(x + c) の係数（mod は f の長さより大きい素数）  O(n log n)
template <class mint>
Fps<mint> taylor_shift(const std::vector<mint>& f, const mint& c) {
    const int n = (int)f.size();
    if (n == 0) return {};
    math_internal::fps_need_prime<mint>(n, "taylor_shift");
    std::vector<mint> fact((size_t)n, mint(1)), ifact((size_t)n, mint(1));
    for (int i = 1; i < n; ++i) fact[i] = fact[i - 1] * mint(i);
    ifact[n - 1] = fact[n - 1].inv();
    for (int i = n - 1; i > 0; --i) ifact[i - 1] = ifact[i] * mint(i);
    std::vector<mint> a((size_t)n), b((size_t)n);
    mint pw = 1;
    for (int i = 0; i < n; ++i) a[i] = f[i] * fact[i], b[(size_t)(n - 1 - i)] = pw * ifact[i], pw *= c;
    const std::vector<mint> conv = convolve(a, b);
    Fps<mint> r((size_t)n);
    for (int j = 0; j < n; ++j) r[j] = conv[(size_t)(n - 1 + j)] * ifact[j];
    return r;
}
// [x^K] P(x)/Q(x)（Q[0] ≠ 0。K ≤ 10^18。Bostan–Mori）  O(d log d log K)（d = 次数）
template <class mint>
mint fps_coeff_of_fraction(std::vector<mint> P, std::vector<mint> Q, long long K) {
    if (K < 0) return mint(0);
    if (Q.empty() || Q[0] == mint(0)) math_internal::fail("fps_coeff_of_fraction", "Q[0] が 0");
    while (K > 0) {
        if (P.empty()) return mint(0);
        const size_t up = P.size() + Q.size() - 1, vq = 2 * Q.size() - 1;
        const size_t n = (size_t)1 << math_internal::ceil_log2(std::max(up, vq));
        if (Q.size() >= 16 && math_internal::ntt_usable<mint>(n)) {  // NTT の値の並びでは j と j^1 が ±x の組（Q(−x) の値は隣と入れ替えるだけ）
            P.resize(n), Q.resize(n);
            ntt(P), ntt(Q);
            for (size_t j = 0; j < n; j += 2) {
                const mint q0 = Q[j], q1 = Q[j + 1];
                P[j] *= q1, P[j + 1] *= q0;
                Q[j] = Q[j + 1] = q0 * q1;
            }
            intt(P), intt(Q);
            P.resize(up), Q.resize(vq);
        } else {
            std::vector<mint> Qm(Q);
            for (size_t i = 1; i < Qm.size(); i += 2) Qm[i] = -Qm[i];  // Q(−x)
            P = convolve(P, Qm), Q = convolve(Q, Qm);
        }
        std::vector<mint> np, nq;
        for (size_t i = (size_t)(K & 1); i < P.size(); i += 2) np.push_back(P[i]);
        for (size_t i = 0; i < Q.size(); i += 2) nq.push_back(Q[i]);
        P.swap(np), Q.swap(nq), K >>= 1;
    }
    return P.empty() ? mint(0) : P[0] / Q[0];
}
// Π (1 + a_i x)（分割統治）  O(n log^2 n)
template <class mint>
Fps<mint> prod_linear(const std::vector<mint>& as) {
    std::vector<std::vector<mint>> ps;
    ps.reserve(as.size());
    for (const mint& a : as) ps.push_back({mint(1), a});
    return Fps<mint>(convolve_many(std::move(ps)));
}

// ---- 多点評価・補間・分数の和 ----
namespace math_internal {
template <class mint>
struct FpsProductTree {  // 区間 [l, r) の点の Π (x − x_i) を持つ木（葉の近くは 32 点ずつまとめる）
    static constexpr int LEAF = 32;
    const std::vector<mint>& xs;
    std::vector<Fps<mint>> t;
    explicit FpsProductTree(const std::vector<mint>& x) : xs(x), t(4 * std::max<size_t>(x.size(), 1)) { build(1, 0, (int)x.size()); }
    void build(int k, int l, int r) {
        if (r - l <= LEAF) {
            Fps<mint> p = {mint(1)};
            for (int i = l; i < r; ++i) {  // p ← p (x − x_i)
                p.push_back(mint(0));
                for (size_t j = p.size() - 1; j > 0; --j) p[j] = p[j - 1] - p[j] * xs[(size_t)i];
                p[0] = -p[0] * xs[(size_t)i];
            }
            t[(size_t)k] = p;
            return;
        }
        const int m = (l + r) / 2;
        build(2 * k, l, m), build(2 * k + 1, m, r);
        t[(size_t)k] = t[(size_t)(2 * k)] * t[(size_t)(2 * k + 1)];
    }
    void eval(int k, int l, int r, Fps<mint> f, std::vector<mint>& out) const {  // f の次数は r − l 未満にしてから
        if (r - l <= LEAF) {
            for (int i = l; i < r; ++i) out[(size_t)i] = f.eval(xs[(size_t)i]);
            return;
        }
        const int m = (l + r) / 2;
        eval(2 * k, l, m, fps_divmod(f, t[(size_t)(2 * k)]).r, out);
        eval(2 * k + 1, m, r, fps_divmod(f, t[(size_t)(2 * k + 1)]).r, out);
    }
    Fps<mint> combine(int k, int l, int r, const std::vector<mint>& c) const {  // Σ c_i Π_{j ≠ i} (x − x_j)
        if (r - l <= LEAF) {
            Fps<mint> res((size_t)(r - l), mint(0));
            for (int i = l; i < r; ++i) {  // t[k] / (x − x_i) を組み立て除法で
                const Fps<mint>& p = t[(size_t)k];
                mint carry = 0;
                for (size_t j = p.size() - 1; j > 0; --j) {
                    carry = p[j] + carry * xs[(size_t)i];
                    res[j - 1] += carry * c[(size_t)i];
                }
            }
            return res;
        }
        const int m = (l + r) / 2;
        return combine(2 * k, l, m, c) * t[(size_t)(2 * k + 1)] + combine(2 * k + 1, m, r, c) * t[(size_t)(2 * k)];
    }
};
}  // namespace math_internal
// 多点評価：f(x_0), f(x_1), …, f(x_{m−1})  O(n log n + m log^2 m)（n = |f|, m = |xs|）
template <class mint>
std::vector<mint> multipoint_eval(const std::vector<mint>& f, const std::vector<mint>& xs) {
    std::vector<mint> out(xs.size());
    if (xs.empty()) return out;
    if ((long long)f.size() * (long long)xs.size() <= 1000000 || std::min(f.size(), xs.size()) <= 64) {  // 小さいとき・片方が短いときは 1 点ずつ
        const Fps<mint> F(f);
        for (size_t i = 0; i < xs.size(); ++i) out[i] = F.eval(xs[i]);
        return out;
    }
    const math_internal::FpsProductTree<mint> T(xs);
    T.eval(1, 0, (int)xs.size(), fps_divmod(f, T.t[1]).r, out);
    return out;
}
// 補間：x_i（すべて異なる）で y_i をとる次数 n − 1 以下の多項式（係数は低い次数から。長さ n）  O(n log^2 n)
template <class mint>
Fps<mint> interpolate(const std::vector<mint>& xs, const std::vector<mint>& ys) {
    if (xs.size() != ys.size()) math_internal::fail("interpolate", "xs と ys の長さが違う");
    const int n = (int)xs.size();
    if (n == 0) return {};
    const math_internal::FpsProductTree<mint> T(xs);
    const std::vector<mint> d = multipoint_eval<mint>(fps_derivative(T.t[1]), xs);  // M'(x_i) = Π_{j ≠ i} (x_i − x_j)
    std::vector<mint> c((size_t)n);
    for (int i = 0; i < n; ++i) {
        if (d[(size_t)i] == mint(0)) math_internal::fail("interpolate", "x が重複している");
        c[(size_t)i] = ys[(size_t)i] / d[(size_t)i];
    }
    Fps<mint> r = T.combine(1, 0, n, c);
    r.resize((size_t)n);
    return r;
}
template <class mint>
struct FpsFraction {  // num / den
    Fps<mint> num, den;
};
// 分数の和 Σ nums[i] / dens[i] を 1 つの分数に（分割統治。通分するだけで約分はしない）  O(N log^2 N)（N = 次数の合計）
template <class mint>
FpsFraction<mint> fps_sum_of_fractions(const std::vector<Fps<mint>>& nums, const std::vector<Fps<mint>>& dens) {
    if (nums.size() != dens.size()) math_internal::fail("fps_sum_of_fractions", "nums と dens の長さが違う");
    if (nums.empty()) return {Fps<mint>{}, Fps<mint>{mint(1)}};
    auto rec = [&](auto&& self, size_t l, size_t r) -> FpsFraction<mint> {
        if (r - l == 1) return {nums[l], dens[l]};
        const size_t m = (l + r) / 2;
        const FpsFraction<mint> a = self(self, l, m), b = self(self, m, r);
        return {a.num * b.den + b.num * a.den, a.den * b.den};
    };
    return rec(rec, 0, nums.size());
}
// べき乗和 p_k = Σ_i w_i a_i^k（k = 0〜n−1。ws を省くと w_i = 1）= Σ w_i / (1 − a_i x) の係数  O(m log^2 m + n log n)
template <class mint>
Fps<mint> power_sums(const std::vector<mint>& as, int n, const std::vector<mint>& ws = {}) {
    if (!ws.empty() && ws.size() != as.size()) math_internal::fail("power_sums", "ws と as の長さが違う");
    if (n <= 0) return {};
    if (as.empty()) return Fps<mint>((size_t)n, mint(0));
    std::vector<Fps<mint>> nums(as.size()), dens(as.size());
    for (size_t i = 0; i < as.size(); ++i) nums[i] = {ws.empty() ? mint(1) : ws[i]}, dens[i] = {mint(1), -as[i]};
    const FpsFraction<mint> fr = fps_sum_of_fractions(nums, dens);
    return (fr.num.pre((size_t)n) * fps_inv(fr.den, n)).pre((size_t)n);
}
#endif  // GG_FPS
