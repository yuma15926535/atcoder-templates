#ifndef GG_SUBSET_TRANSFORM
#define GG_SUBSET_TRANSFORM
// ---- 部分集合の変換（数学第二世代）  高速ゼータ・メビウス変換、OR・AND・XOR の畳み込み、部分集合の畳み込み・exp・log ----
// 配列の長さは 2^n（i 番目が集合 i）。T は mint・long long（溢れたら止める）・double など
namespace math_internal {
template <class T>
inline T st_add(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_add_overflow(a, b, &r)) fail("subset_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a + b;
    }
}
template <class T>
inline T st_sub(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_sub_overflow(a, b, &r)) fail("subset_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a - b;
    }
}
template <class T>
inline T st_mul(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T r;
        if (__builtin_mul_overflow(a, b, &r)) fail("subset_transform", "値が整数の型に収まらない（mint で）");
        return r;
    } else {
        return a * b;
    }
}
inline int st_log2(size_t n, const char* func) {  // n = 2^k の k（2 の冪でなければ止める）
    if (n == 0 || (n & (n - 1))) fail(func, "長さは 2 の冪（2^n）");
    return __builtin_ctzll((unsigned long long)n);
}
struct StAdd {
    template <class T>
    T operator()(const T& a, const T& b) const { return st_add(a, b); }
};
struct StSub {
    template <class T>
    T operator()(const T& a, const T& b) const { return st_sub(a, b); }
};
}  // namespace math_internal

// f[S] ← Σ_{T ⊆ S} f[T]（op を渡すと Σ の代わりに op でまとめる：max・min・gcd・OR など）  O(n 2^n)
template <class T, class Op = math_internal::StAdd>
void subset_zeta(std::vector<T>& f, Op op = Op{}) {
    const size_t n = f.size();
    math_internal::st_log2(n, "subset_zeta");
    for (size_t w = 1; w < n; w <<= 1)
        for (size_t s = 0; s < n; ++s)
            if (s & w) f[s] = op(f[s], f[s ^ w]);
}
// subset_zeta の逆（f[S] ← Σ_{T ⊆ S} (−1)^{|S∖T|} f[T]。inv_op は op の逆）
template <class T, class Op = math_internal::StSub>
void subset_mobius(std::vector<T>& f, Op inv_op = Op{}) {
    const size_t n = f.size();
    math_internal::st_log2(n, "subset_mobius");
    for (size_t w = 1; w < n; w <<= 1)
        for (size_t s = 0; s < n; ++s)
            if (s & w) f[s] = inv_op(f[s], f[s ^ w]);
}
// f[S] ← Σ_{T ⊇ S} f[T]
template <class T, class Op = math_internal::StAdd>
void superset_zeta(std::vector<T>& f, Op op = Op{}) {
    const size_t n = f.size();
    math_internal::st_log2(n, "superset_zeta");
    for (size_t w = 1; w < n; w <<= 1)
        for (size_t s = 0; s < n; ++s)
            if (!(s & w)) f[s] = op(f[s], f[s | w]);
}
template <class T, class Op = math_internal::StSub>
void superset_mobius(std::vector<T>& f, Op inv_op = Op{}) {
    const size_t n = f.size();
    math_internal::st_log2(n, "superset_mobius");
    for (size_t w = 1; w < n; w <<= 1)
        for (size_t s = 0; s < n; ++s)
            if (!(s & w)) f[s] = inv_op(f[s], f[s | w]);
}
// アダマール変換（XOR の畳み込み用。inverse なら 2^n で割る：T が整数なら割り切れることを確かめる）
template <class T>
void walsh_hadamard(std::vector<T>& f, bool inverse = false) {
    const size_t n = f.size();
    math_internal::st_log2(n, "walsh_hadamard");
    for (size_t w = 1; w < n; w <<= 1)
        for (size_t s = 0; s < n; ++s)
            if (!(s & w)) {
                const T a = f[s], b = f[s | w];
                f[s] = math_internal::st_add(a, b), f[s | w] = math_internal::st_sub(a, b);
            }
    if (inverse) {
        if constexpr (std::is_integral_v<T>) {
            for (T& x : f) {
                if (x % (T)n != 0) math_internal::fail("walsh_hadamard", "2^n で割り切れない");
                x /= (T)n;
            }
        } else {
            const T inv = T(1) / T((long long)n);
            for (T& x : f) x = x * inv;
        }
    }
}
namespace math_internal {
template <class T>
void st_pad(std::vector<T>& a, std::vector<T>& b) {  // 同じ 2 の冪の長さにそろえる
    size_t n = 1;
    while (n < std::max(a.size(), b.size())) n <<= 1;
    a.resize(n, T(0)), b.resize(n, T(0));
}
}  // namespace math_internal
// c[k] = Σ_{i | j = k} a[i] b[j]（長さは 2 の冪にそろえる）  O(n 2^n)
template <class T>
std::vector<T> or_convolution(std::vector<T> a, std::vector<T> b) {
    if (a.empty() || b.empty()) return {};
    math_internal::st_pad(a, b);
    subset_zeta(a), subset_zeta(b);
    for (size_t i = 0; i < a.size(); ++i) a[i] = math_internal::st_mul(a[i], b[i]);
    subset_mobius(a);
    return a;
}
// c[k] = Σ_{i & j = k} a[i] b[j]
template <class T>
std::vector<T> and_convolution(std::vector<T> a, std::vector<T> b) {
    if (a.empty() || b.empty()) return {};
    math_internal::st_pad(a, b);
    superset_zeta(a), superset_zeta(b);
    for (size_t i = 0; i < a.size(); ++i) a[i] = math_internal::st_mul(a[i], b[i]);
    superset_mobius(a);
    return a;
}
// c[k] = Σ_{i ^ j = k} a[i] b[j]
template <class T>
std::vector<T> xor_convolution(std::vector<T> a, std::vector<T> b) {
    if (a.empty() || b.empty()) return {};
    math_internal::st_pad(a, b);
    walsh_hadamard(a), walsh_hadamard(b);
    for (size_t i = 0; i < a.size(); ++i) a[i] = math_internal::st_mul(a[i], b[i]);
    walsh_hadamard(a, true);
    return a;
}
// 部分集合の畳み込み c[S] = Σ_{T ⊆ S} a[T] b[S∖T]（ランク付きのゼータ変換）  O(n^2 2^n)
template <class T>
std::vector<T> subset_convolution(std::vector<T> a, std::vector<T> b) {
    if (a.empty() || b.empty()) return {};
    math_internal::st_pad(a, b);
    const size_t N = a.size();
    const int n = math_internal::st_log2(N, "subset_convolution");
    std::vector<std::vector<T>> ra((size_t)n + 1, std::vector<T>(N, T(0))), rb = ra, rc = ra;
    for (size_t s = 0; s < N; ++s) {
        const int k = __builtin_popcountll((unsigned long long)s);
        ra[(size_t)k][s] = a[s], rb[(size_t)k][s] = b[s];
    }
    for (int k = 0; k <= n; ++k) subset_zeta(ra[(size_t)k]), subset_zeta(rb[(size_t)k]);
    for (int i = 0; i <= n; ++i)
        for (int j = 0; i + j <= n; ++j)
            for (size_t s = 0; s < N; ++s) rc[(size_t)(i + j)][s] = math_internal::st_add(rc[(size_t)(i + j)][s], math_internal::st_mul(ra[(size_t)i][s], rb[(size_t)j][s]));
    for (int k = 0; k <= n; ++k) subset_mobius(rc[(size_t)k]);
    std::vector<T> c(N, T(0));
    for (size_t s = 0; s < N; ++s) c[s] = rc[(size_t)__builtin_popcountll((unsigned long long)s)][s];
    return c;
}
namespace math_internal {
template <class T>
std::vector<T> st_subset_inv(const std::vector<T>& a) {  // subset_convolution の逆元（a[∅] = 1。割り算はいらない）  O(n^2 2^n)
    const size_t N = a.size();
    const int n = st_log2(N, "subset_log");
    std::vector<std::vector<T>> r((size_t)n + 1, std::vector<T>(N, T(0)));
    for (size_t s = 0; s < N; ++s) r[(size_t)__builtin_popcountll((unsigned long long)s)][s] = a[s];
    for (int k = 0; k <= n; ++k) subset_zeta(r[(size_t)k]);
    std::vector<T> q((size_t)n + 1);
    for (size_t s = 0; s < N; ++s) {  // 各 s で多項式 Σ_k r[k][s] z^k（定数項 1）の逆数
        q[0] = T(1);
        for (int m = 1; m <= n; ++m) {
            T acc = T(0);
            for (int k = 1; k <= m; ++k) acc = st_add(acc, st_mul(r[(size_t)k][s], q[(size_t)(m - k)]));
            q[(size_t)m] = st_sub(T(0), acc);
        }
        for (int m = 0; m <= n; ++m) r[(size_t)m][s] = q[(size_t)m];
    }
    for (int k = 0; k <= n; ++k) subset_mobius(r[(size_t)k]);
    std::vector<T> c(N);
    for (size_t s = 0; s < N; ++s) c[s] = r[(size_t)__builtin_popcountll((unsigned long long)s)][s];
    return c;
}
}  // namespace math_internal
// 集合のべき級数の exp：g[S] = Σ_{S の組分け {B_1, …, B_k}} f[B_1] … f[B_k]（f[∅] = 0、g[∅] = 1。長さ 2^n）
//   「つながった部分の数え方」f から「全体の数え方」g へ。割り算を使わないので long long・合成数の mod でも  O(n^2 2^n)
template <class T>
std::vector<T> subset_exp(const std::vector<T>& f) {
    const int n = math_internal::st_log2(f.size(), "subset_exp");
    if (!(f[0] == T(0))) math_internal::fail("subset_exp", "f[∅] = 0 で使う");
    std::vector<T> g = {T(1)};
    for (int k = 0; k < n; ++k) {  // 要素 k を含む集合：k を含む組 {k} ∪ T と残りの組分け
        const std::vector<T> h(f.begin() + (1LL << k), f.begin() + (2LL << k));
        const std::vector<T> hi = subset_convolution(h, g);
        g.insert(g.end(), hi.begin(), hi.end());
    }
    return g;
}
// subset_exp の逆：f[∅] = 1 の f から g（f = subset_exp(g)、g[∅] = 0）。「全体」から「つながったもの」を数える  O(n^2 2^n)
template <class T>
std::vector<T> subset_log(const std::vector<T>& f) {
    const int n = math_internal::st_log2(f.size(), "subset_log");
    if (!(f[0] == T(1))) math_internal::fail("subset_log", "f[∅] = 1 で使う");
    std::vector<T> g = {T(0)};
    for (int k = 0; k < n; ++k) {  // f の上半分 = g の上半分 * f の下半分 なので、f の下半分の逆元を掛ける
        const std::vector<T> lo(f.begin(), f.begin() + (1LL << k)), hi(f.begin() + (1LL << k), f.begin() + (2LL << k));
        const std::vector<T> gh = subset_convolution(hi, math_internal::st_subset_inv(lo));
        g.insert(g.end(), gh.begin(), gh.end());
    }
    return g;
}
#endif  // GG_SUBSET_TRANSFORM
