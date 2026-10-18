#ifndef GG_SERIES
#define GG_SERIES
// ---- べき乗和・ラグランジュ補間・次数の小さい多項式（数学第二世代） ----
namespace math_internal {
template <class mint>
void check_prime_mod_over(long long need, const char* func) {  // mod が素数で need より大きいか（階乗の逆元を使うため）
    if (!is_prime_u32((unsigned)mint::mod())) fail(func, "mod が素数でない");
    if (need >= mint::mod()) fail(func, "点の数が mod 以上");
}
template <class mint>
std::vector<mint> inv_factorials(int n) {  // 1/0!〜1/n!
    std::vector<mint> f((size_t)n + 1, mint(1));
    for (int i = 1; i <= n; ++i) f[(size_t)i] = f[(size_t)i - 1] * mint(i);
    f[(size_t)n] = f[(size_t)n].inv();
    for (int i = n; i > 0; --i) f[(size_t)i - 1] = f[(size_t)i] * mint(i);
    return f;
}
}  // namespace math_internal

// ---- ラグランジュ補間 ----
// f(0), f(1), …, f(d) が ys のとき、d 次以下の多項式 f の f(x)（x は 10^18 でも負でもよい。mod は d より大きい素数）  O(d)
template <class mint>
mint lagrange_interpolate(const std::vector<mint>& ys, long long x) {
    const int n = (int)ys.size();  // 点の数 = d + 1
    if (n == 0) return mint(0);
    math_internal::check_prime_mod_over<mint>(n - 1, "lagrange_interpolate");
    if (0 <= x && x < n) return ys[(size_t)x];
    const mint X(x);
    std::vector<mint> pre((size_t)n + 1, mint(1)), suf((size_t)n + 1, mint(1));  // Π (x − j) の左右からの積
    for (int i = 0; i < n; ++i) pre[(size_t)i + 1] = pre[(size_t)i] * (X - mint(i));
    for (int i = n - 1; i >= 0; --i) suf[(size_t)i] = suf[(size_t)i + 1] * (X - mint(i));
    const std::vector<mint> fi = math_internal::inv_factorials<mint>(n - 1);
    mint r = 0;
    for (int i = 0; i < n; ++i) {
        mint t = ys[(size_t)i] * pre[(size_t)i] * suf[(size_t)i + 1] * fi[(size_t)i] * fi[(size_t)(n - 1 - i)];
        r += ((n - 1 - i) % 2 == 0) ? t : -t;
    }
    return r;
}
// 点 (xs[i], ys[i])（xs は互いに違う）を通る多項式の x での値（mod は素数）  O(d^2)
template <class mint>
mint lagrange_interpolate_at(const std::vector<mint>& xs, const std::vector<mint>& ys, mint x) {
    if (xs.size() != ys.size()) math_internal::fail("lagrange_interpolate_at", "xs と ys の長さが違う");
    const size_t n = xs.size();
    mint r = 0;
    for (size_t i = 0; i < n; ++i) {
        mint num = 1, den = 1;
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            if (xs[i] == xs[j]) math_internal::fail("lagrange_interpolate_at", "xs に同じ値がある");
            num *= x - xs[j], den *= xs[i] - xs[j];
        }
        r += ys[i] * num / den;
    }
    return r;
}
// 点 (xs[i], ys[i]) を通る多項式の係数（低い次数から。長さ n）  O(n^2)
template <class mint>
std::vector<mint> poly_from_points(const std::vector<mint>& xs, const std::vector<mint>& ys) {
    if (xs.size() != ys.size()) math_internal::fail("poly_from_points", "xs と ys の長さが違う");
    const size_t n = xs.size();
    std::vector<mint> all(n + 1, mint(0)), res(n, mint(0));  // all = Π (x − xs[j])
    all[0] = 1;
    for (size_t j = 0; j < n; ++j) {  // all ← all × (x − xs[j])（上の次数から更新）
        for (size_t k = j + 1; k > 0; --k) all[k] = all[k - 1] - xs[j] * all[k];
        all[0] = -xs[j] * all[0];
    }
    for (size_t i = 0; i < n; ++i) {
        mint den = 1;
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            if (xs[i] == xs[j]) math_internal::fail("poly_from_points", "xs に同じ値がある");
            den *= xs[i] - xs[j];
        }
        const mint c = ys[i] / den;
        // q = all / (x − xs[i])（組立除法）を c 倍して足す
        mint carry = 0;
        for (size_t k = n; k > 0; --k) {
            carry = all[k] + carry * xs[i];  // q の k−1 次の係数
            res[k - 1] += c * carry;
        }
    }
    return res;
}

// ---- べき乗和 ----
// Σ_{i=0}^{n−1} i^k（0^0 = 1。n ≤ 10^18、k ≤ 10^7 程度。mod は k + 1 より大きい素数）  O(k)
template <class mint>
mint sum_powers(long long n, int k) {
    if (n < 0 || k < 0) math_internal::fail("sum_powers", "n, k ≥ 0 で使う");
    // f(x) = Σ_{i<x} i^k は k + 1 次の多項式。f(0)〜f(k+1) を作って補間
    const int m = k + 2;
    if ((long long)m - 1 >= mint::mod()) math_internal::fail("sum_powers", "k + 1 が mod 以上");
    std::vector<mint> pw((size_t)m, mint(0));  // pw[i] = i^k（線形篩で：素数だけ累乗）
    std::vector<int> spf((size_t)m, 0), primes;
    if (m > 0) pw[0] = (k == 0 ? mint(1) : mint(0));
    if (m > 1) pw[1] = 1;
    for (int i = 2; i < m; ++i) {
        if (spf[(size_t)i] == 0) spf[(size_t)i] = i, primes.push_back(i), pw[(size_t)i] = mint(i).pow(k);
        for (int p : primes) {
            if (p > spf[(size_t)i] || (long long)p * i >= m) break;
            spf[(size_t)(p * i)] = p, pw[(size_t)(p * i)] = pw[(size_t)p] * pw[(size_t)i];
        }
    }
    std::vector<mint> ys((size_t)m, mint(0));
    for (int i = 1; i < m; ++i) ys[(size_t)i] = ys[(size_t)i - 1] + pw[(size_t)i - 1];
    return lagrange_interpolate(ys, n);
}
template <class mint>
mint sum_powers(long long l, long long r, int k) {  // Σ_{i=l}^{r−1} i^k（0 ≤ l。空なら 0）
    if (l < 0) math_internal::fail("sum_powers", "l ≥ 0 で使う");
    if (r <= l) return mint(0);
    return sum_powers<mint>(r, k) - sum_powers<mint>(l, k);
}
template <class mint>
mint sum_powers_closed(long long L, long long R, int k) {  // Σ_{i=L}^{R} i^k（0 ≤ L ≤ R < 2^63 − 1）
    if (L < 0) math_internal::fail("sum_powers_closed", "L ≥ 0 で使う");
    if (R < L) return mint(0);
    return sum_powers<mint>(L, R + 1, k);
}
// Σ_{i=0}^{n−1} i·r^i（mod は何でもよい）  O(log n)
template <class mint>
mint arith_geom_sum(mint r, long long n) {
    if (n < 0) math_internal::fail("arith_geom_sum", "n ≥ 0 で使う");
    // (len, S0 = Σ r^i, S1 = Σ i r^i, P = r^len) を上の桁から 2 倍・+1 で作る
    mint s0 = 0, s1 = 0, p = 1;
    long long len = 0;
    for (int b = 62; b >= 0; --b) {
        // 2 倍：[0, 2len) = [0, len) + len + [0, len)
        s1 = s1 + p * (s1 + mint(len) * s0), s0 = s0 * (mint(1) + p), p = p * p, len *= 2;
        if ((n >> b) & 1) s1 += mint(len) * p, s0 += p, p *= r, ++len;  // 1 項足す（i = len の項 len·r^len）
    }
    return s1;
}

// ---- 次数の小さい多項式（係数は低い次数から。T は mint・long long・double など） ----
template <class T>
T poly_eval(const std::vector<T>& c, const T& x) {  // c(x)（ホーナー法）  O(n)
    T r = T(0);
    for (size_t i = c.size(); i > 0; --i) r = r * x + c[i - 1];
    return r;
}
template <class T>
std::vector<T> poly_mul_naive(const std::vector<T>& a, const std::vector<T>& b) {  // a × b（長さ |a| + |b| − 1）  O(|a||b|)
    if (a.empty() || b.empty()) return {};
    std::vector<T> c(a.size() + b.size() - 1, T(0));
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j) c[i + j] = c[i + j] + a[i] * b[j];
    return c;
}
template <class T>
struct PolyDivMod {  // a = q × b + r（deg r < deg b）
    std::vector<T> q, r;
};
// a ÷ b（b の最高次の係数で割る。T が整数なら割り切れないと止める。商・余りの末尾の 0 は除く）  O(|a||b|)
template <class T>
PolyDivMod<T> poly_divmod(std::vector<T> a, std::vector<T> b) {
    auto trim = [](std::vector<T>& v) {
        while (!v.empty() && v.back() == T(0)) v.pop_back();
    };
    trim(a), trim(b);
    if (b.empty()) math_internal::fail("poly_divmod", "0 で割った");
    PolyDivMod<T> res;
    if (a.size() < b.size()) {
        res.r = a;
        return res;
    }
    const size_t n = a.size(), m = b.size();
    res.q.assign(n - m + 1, T(0));
    const T lc = b.back();
    for (size_t i = n - m + 1; i > 0; --i) {
        const size_t top = i - 1 + m - 1;  // 今の最高次
        T c;
        if constexpr (std::is_integral_v<T>) {
            if (a[top] % lc != 0) math_internal::fail("poly_divmod", "整数の多項式で割り切れない（最高次の係数で割れない）");
            c = a[top] / lc;
        } else {
            c = a[top] / lc;
        }
        res.q[i - 1] = c;
        for (size_t j = 0; j < m; ++j) a[i - 1 + j] = a[i - 1 + j] - c * b[j];
    }
    a.resize(m - 1);
    trim(a), trim(res.q);
    res.r = a;
    return res;
}
#endif  // GG_SERIES
