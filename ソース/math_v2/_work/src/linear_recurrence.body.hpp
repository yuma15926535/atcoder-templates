#ifndef GG_LINEAR_RECURRENCE
#define GG_LINEAR_RECURRENCE
// ---- 線形漸化式（数学第二世代）  a_i = c_1 a_{i−1} + … + c_d a_{i−d} の第 K 項（Bostan–Mori）・最初の数項から漸化式を求める（Berlekamp–Massey） ----
// c は {c_1, c_2, …, c_d}（c[0] が 1 つ前の項の係数）。a は最初の d 項 a_0〜a_{d−1}
namespace math_internal {
template <class mint>
mint linrec_bostan_mori(std::vector<mint> P, std::vector<mint> Q, unsigned long long K) {  // [x^K] P/Q（deg P < deg Q = d、Q[0] = 1）
    const size_t d = Q.size() - 1;
    P.resize(d);
    while (K > 0) {
        const size_t n = (size_t)1 << ceil_log2(2 * d + 1);
        if (d >= 16 && ntt_usable<mint>(n)) {  // NTT の値の並びでは j と j^1 が ±x の組（Q(−x) の値は隣と入れ替えるだけ）
            P.resize(n), Q.resize(n);
            ntt(P), ntt(Q);
            for (size_t j = 0; j < n; j += 2) {
                const mint q0 = Q[j], q1 = Q[j + 1];
                P[j] *= q1, P[j + 1] *= q0;
                Q[j] = Q[j + 1] = q0 * q1;
            }
            intt(P), intt(Q);
        } else {
            std::vector<mint> Qm(Q);
            for (size_t i = 1; i < Qm.size(); i += 2) Qm[i] = -Qm[i];  // Q(−x)
            P = d == 0 ? std::vector<mint>() : convolve(P, Qm), Q = convolve(Q, Qm);
        }
        // P(x) Q(−x) の K の偶奇と同じ次数の項、Q(x) Q(−x) の偶数次の項
        std::vector<mint> np, nq;
        np.reserve(d), nq.reserve(d + 1);
        for (size_t i = (size_t)(K & 1); i < P.size() && np.size() < d; i += 2) np.push_back(P[i]);
        for (size_t i = 0; i < Q.size() && nq.size() < d + 1; i += 2) nq.push_back(Q[i]);
        np.resize(d), nq.resize(d + 1);
        P.swap(np), Q.swap(nq), K >>= 1;
    }
    return d == 0 ? mint(0) : P[0];
}
template <class mint>
std::vector<mint> linrec_series_inv(const std::vector<mint>& Q, size_t n) {  // 1/Q の最初の n 項（Q[0] = 1。掛け算だけのニュートン法なので mod は何でも）
    std::vector<mint> g{mint(1)};
    for (size_t k = 1; k < n;) {
        k *= 2;
        std::vector<mint> q(Q.begin(), Q.begin() + (long long)std::min(Q.size(), k));
        std::vector<mint> t = convolve(q, g);  // g ← g (2 − Q g)
        t.resize(k);
        for (auto& x : t) x = -x;
        t[0] += mint(2);
        g = convolve(g, t);
        g.resize(k);
    }
    g.resize(n);
    return g;
}
template <class mint>
void linrec_check(const std::vector<mint>& a, const std::vector<mint>& c, const char* func) {
    if (a.size() < c.size()) fail(func, "a は c と同じ数（d 個）以上の項が必要");
}
}  // namespace math_internal

// 第 K 項 a_K（a は最初の d 項以上、c = {c_1, …, c_d}。K ≤ 10^18。mod は何でも）  O(d log d log K)
template <class mint>
mint linear_recurrence_kth(const std::vector<mint>& a, const std::vector<mint>& c, long long K) {
    math_internal::linrec_check(a, c, "linear_recurrence_kth");
    if (K < 0) math_internal::fail("linear_recurrence_kth", "K ≥ 0 で使う");
    if (K < (long long)a.size()) return a[(size_t)K];
    const size_t d = c.size();
    if (d == 0) return mint(0);
    std::vector<mint> Q(d + 1), A(a.begin(), a.begin() + (long long)d);  // Q = 1 − c_1 x − … − c_d x^d、P = A Q mod x^d
    Q[0] = 1;
    for (size_t j = 0; j < d; ++j) Q[j + 1] = -c[j];
    std::vector<mint> P = convolve(A, Q);
    P.resize(d);
    return math_internal::linrec_bostan_mori(P, Q, (unsigned long long)K);
}
// 最初の n 項 a_0〜a_{n−1}  O(n log n)（d が小さければ O(n d)）
template <class mint>
std::vector<mint> linear_recurrence_terms(const std::vector<mint>& a, const std::vector<mint>& c, long long n) {
    math_internal::linrec_check(a, c, "linear_recurrence_terms");
    if (n < 0) math_internal::fail("linear_recurrence_terms", "n ≥ 0 で使う");
    const size_t d = c.size(), N = (size_t)n;
    std::vector<mint> res(a.begin(), a.begin() + (long long)std::min(a.size(), N));
    if (res.size() == N) return res;
    if (d == 0) {
        res.resize(N);
        return res;
    }
    if ((double)d * (double)N <= 3e7 || d <= 64) {  // 漸化式どおりに
        res.resize(N);
        for (size_t i = a.size(); i < N; ++i) {
            mint s = 0;
            for (size_t j = 0; j < d; ++j) s += c[j] * res[i - 1 - j];
            res[i] = s;
        }
        return res;
    }
    std::vector<mint> Q(d + 1), A(a.begin(), a.begin() + (long long)d);
    Q[0] = 1;
    for (size_t j = 0; j < d; ++j) Q[j + 1] = -c[j];
    std::vector<mint> P = convolve(A, Q);
    P.resize(d);
    std::vector<mint> r = convolve(P, math_internal::linrec_series_inv(Q, N));
    r.resize(N);
    return r;
}
// 列 s をみたす最短の線形漸化式の c = {c_1, …, c_d}（mod は素数。2d 項以上あれば正しい漸化式が出る）  O(N^2)
template <class mint>
std::vector<mint> berlekamp_massey(const std::vector<mint>& s) {
    std::vector<mint> C{mint(1)}, B{mint(1)};
    mint b = 1;
    size_t L = 0, m = 1;
    for (size_t n = 0; n < s.size(); ++n) {
        mint dlt = s[n];
        for (size_t i = 1; i <= L; ++i) dlt += C[i] * s[n - i];
        if (dlt == mint(0)) {
            ++m;
            continue;
        }
        const mint coef = dlt / b;
        if (2 * L <= n) {
            const std::vector<mint> T = C;
            if (C.size() < B.size() + m) C.resize(B.size() + m, mint(0));
            for (size_t i = 0; i < B.size(); ++i) C[i + m] -= coef * B[i];
            L = n + 1 - L, B = T, b = dlt, m = 1;
        } else {
            if (C.size() < B.size() + m) C.resize(B.size() + m, mint(0));
            for (size_t i = 0; i < B.size(); ++i) C[i + m] -= coef * B[i];
            ++m;
        }
    }
    C.resize(L + 1, mint(0));
    std::vector<mint> c(L);
    for (size_t j = 1; j <= L; ++j) c[j - 1] = -C[j];
    return c;
}
// 最初の数項 s から漸化式を推測して第 K 項（漸化式の次数 d に対して 2d 項以上いる。K ≤ 10^18）
template <class mint>
mint guess_kth(const std::vector<mint>& s, long long K) {
    if (K < 0) math_internal::fail("guess_kth", "K ≥ 0 で使う");
    if (K < (long long)s.size()) return s[(size_t)K];
    return linear_recurrence_kth(s, berlekamp_massey(s), K);
}
#endif  // GG_LINEAR_RECURRENCE
