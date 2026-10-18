#ifndef GG_LINEAR_ALGEBRA
#define GG_LINEAR_ALGEBRA
// ---- 線形代数（数学第二世代）  体（mint・double・Fraction）の上のガウスの消去法・行列式・逆行列・連立一次方程式・行列木定理・特性多項式 ----
// T は mint（mod は素数）・double・long double・08_fraction の Fraction など「割り算ができる」型。double は絶対値最大のピボット（|x| < 1e-12 は 0）
namespace math_internal {
template <class T>
bool la_zero(const T& x) {
    if constexpr (std::is_floating_point<T>::value) return std::abs(x) < (T)1e-12;
    else return x == T(0);
}
template <class T>
int la_pivot(const Matrix<T>& a, int col, int from) {  // col 列の from 行目以降でピボットにする行（なければ −1）
    static_assert(!std::is_integral<T>::value, "整数型の行列では割り算が切り捨てになる。mint・double を使う（合成数の mod の行列式は det_mod_any）");
    int p = -1;
    for (int r = from; r < a.rows(); ++r) {
        if (la_zero(a(r, col))) continue;
        if constexpr (std::is_floating_point<T>::value) {
            if (p < 0 || std::abs(a(r, col)) > std::abs(a(p, col))) p = r;
        } else {
            return r;
        }
    }
    return p;
}
template <class T>
void la_swap_rows(Matrix<T>& a, int r1, int r2) {
    if (r1 == r2) return;
    for (int j = 0; j < a.cols(); ++j) std::swap(a(r1, j), a(r2, j));
}
// 行の基本変形で簡約階段形にする（左の cols_used 列だけでピボットを選ぶ）。ピボットの列を返す
template <class T>
std::vector<int> la_rref(Matrix<T>& a, int cols_used) {
    const int n = a.rows(), m = a.cols();
    std::vector<int> piv;
    int r = 0;
    for (int c = 0; c < cols_used && r < n; ++c) {
        const int p = la_pivot(a, c, r);
        if (p < 0) continue;
        la_swap_rows(a, p, r);
        const T iv = T(1) / a(r, c);
        for (int j = c; j < m; ++j) a(r, j) = a(r, j) * iv;
        T* ar = a[r];
        for (int i = 0; i < n; ++i) {
            if (i == r) continue;
            const T f = a(i, c);
            if (la_zero(f)) continue;
            T* ai = a[i];
            for (int j = c; j < m; ++j) ai[j] = ai[j] - f * ar[j];
        }
        piv.push_back(c), ++r;
    }
    return piv;
}
}  // namespace math_internal

// 階数  O(n m min(n, m))
template <class T>
int mat_rank(Matrix<T> a) { return (int)math_internal::la_rref(a, a.cols()).size(); }
// 行列式（正方行列）  O(n^3)
template <class T>
T mat_det(Matrix<T> a) {
    if (a.rows() != a.cols()) math_internal::fail("mat_det", "正方行列で使う");
    const int n = a.rows();
    T d = T(1);
    for (int c = 0; c < n; ++c) {
        const int p = math_internal::la_pivot(a, c, c);
        if (p < 0) return T(0);
        if (p != c) math_internal::la_swap_rows(a, p, c), d = -d;
        d = d * a(c, c);
        const T iv = T(1) / a(c, c);
        const T* ac = a[c];
        for (int r = c + 1; r < n; ++r) {
            const T f = a(r, c) * iv;
            if (math_internal::la_zero(f)) continue;
            T* ar = a[r];
            for (int j = c; j < n; ++j) ar[j] = ar[j] - f * ac[j];
        }
    }
    return d;
}
// 逆行列（なければ ok = false）  O(n^3)
template <class T>
Maybe<Matrix<T>> mat_inverse(const Matrix<T>& a) {
    if (a.rows() != a.cols()) math_internal::fail("mat_inverse", "正方行列で使う");
    const int n = a.rows();
    Matrix<T> b(n, 2 * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) b(i, j) = a(i, j);
        b(i, n + i) = T(1);
    }
    if ((int)math_internal::la_rref(b, n).size() < n) return {false, Matrix<T>()};
    Matrix<T> r(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) r(i, j) = b(i, n + j);
    return {true, r};
}
// A x = b の解：ok（解があるか）、x（解の 1 つ。自由な変数は 0）、kernel（A y = 0 の解の基底。解は x + kernel の一次結合で全部）
template <class T>
struct LinearSolution {
    bool ok = false;
    std::vector<T> x;
    std::vector<std::vector<T>> kernel;
};
template <class T>
LinearSolution<T> solve_linear(const Matrix<T>& A, const std::vector<T>& b) {  // A は n × m、b の長さは n  O(n m min(n, m))
    const int n = A.rows(), m = A.cols();
    if ((int)b.size() != n) math_internal::fail("solve_linear", "b の長さが A の行数と違う");
    Matrix<T> a(n, m + 1);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) a(i, j) = A(i, j);
        a(i, m) = b[(size_t)i];
    }
    const std::vector<int> piv = math_internal::la_rref(a, m);
    const int rk = (int)piv.size();
    LinearSolution<T> res;
    for (int i = rk; i < n; ++i)
        if (!math_internal::la_zero(a(i, m))) return res;  // 0 = (0 でない数) の行
    res.ok = true;
    res.x.assign((size_t)m, T(0));
    for (int i = 0; i < rk; ++i) res.x[(size_t)piv[(size_t)i]] = a(i, m);
    std::vector<bool> is_piv((size_t)m, false);
    for (int c : piv) is_piv[(size_t)c] = true;
    for (int f = 0; f < m; ++f) {
        if (is_piv[(size_t)f]) continue;
        std::vector<T> v((size_t)m, T(0));
        v[(size_t)f] = T(1);
        for (int i = 0; i < rk; ++i) v[(size_t)piv[(size_t)i]] = -a(i, f);
        res.kernel.push_back(v);
    }
    return res;
}

// 行列式を合成数の mod で（m ≥ 1。m ≤ 10^18 程度。割り算を使わず互除法で消す）  O(n^3 log m)
inline long long det_mod_any(Matrix<long long> a, long long m) {
    if (a.rows() != a.cols()) math_internal::fail("det_mod_any", "正方行列で使う");
    if (m < 1) math_internal::fail("det_mod_any", "m ≥ 1 で使う");
    const int n = a.rows();
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) a(i, j) = mod_floor(a(i, j), m);
    long long d = 1 % m;
    for (int c = 0; c < n; ++c) {
        for (int r = c + 1; r < n; ++r) {
            while (a(r, c) != 0) {  // 行 c と行 r で互除法（c 列目を小さくしていく）
                const long long q = a(c, c) / a(r, c);
                if (q != 0)
                    for (int j = c; j < n; ++j) a(c, j) = mod_floor(a(c, j) - mod_mul(q, a(r, j), m), m);
                math_internal::la_swap_rows(a, c, r), d = (m - d) % m;
            }
        }
        d = mod_mul(d, a(c, c), m);
        if (d == 0) return 0;
    }
    return d;
}

// ---- 行列木定理 ----
// 無向グラフ（頂点 0〜n−1、辺 us[i]–vs[i]。多重辺は別々に数える。自己ループは無視）の全域木の数。T は mint か double  O(n^3)
template <class T>
T count_spanning_trees(int n, const std::vector<int>& us, const std::vector<int>& vs) {
    if (n < 1 || us.size() != vs.size()) math_internal::fail("count_spanning_trees", "n ≥ 1、us と vs は同じ長さ");
    Matrix<T> L(n - 1, n - 1);  // 頂点 n−1 の行と列を除いたラプラシアン
    for (size_t i = 0; i < us.size(); ++i) {
        const int u = us[i], v = vs[i];
        if (u < 0 || u >= n || v < 0 || v >= n) math_internal::fail("count_spanning_trees", "頂点の番号は 0〜n−1");
        if (u == v) continue;
        if (u < n - 1) L(u, u) = L(u, u) + T(1);
        if (v < n - 1) L(v, v) = L(v, v) + T(1);
        if (u < n - 1 && v < n - 1) L(u, v) = L(u, v) - T(1), L(v, u) = L(v, u) - T(1);
    }
    return mat_det(L);
}
// 有向グラフ（辺 us[i] → vs[i]）で root から全頂点へ向かう全域木（有向木）の数  O(n^3)
template <class T>
T count_arborescences(int n, const std::vector<int>& us, const std::vector<int>& vs, int root) {
    if (n < 1 || us.size() != vs.size() || root < 0 || root >= n) math_internal::fail("count_arborescences", "n ≥ 1、0 ≤ root < n、us と vs は同じ長さ");
    auto id = [&](int v) { return v < root ? v : v - 1; };  // root を除いた番号
    Matrix<T> L(n - 1, n - 1);                              // 入次数 − 隣接行列
    for (size_t i = 0; i < us.size(); ++i) {
        const int u = us[i], v = vs[i];
        if (u < 0 || u >= n || v < 0 || v >= n) math_internal::fail("count_arborescences", "頂点の番号は 0〜n−1");
        if (u == v || v == root) continue;
        L(id(v), id(v)) = L(id(v), id(v)) + T(1);
        if (u != root) L(id(u), id(v)) = L(id(u), id(v)) - T(1);
    }
    return mat_det(L);
}

// ---- 特性多項式・det(A + xB) ----
// det(x I − A) の係数（低い次数から。長さ n + 1、最高次は 1）。ヘッセンベルグ形にしてから  O(n^3)
template <class T>
std::vector<T> char_poly(Matrix<T> a) {
    if (a.rows() != a.cols()) math_internal::fail("char_poly", "正方行列で使う");
    const int n = a.rows();
    for (int i = 0; i + 2 < n; ++i) {  // 相似変換で上ヘッセンベルグ形に（i+2 行目以降の i 列目を 0 に）
        const int p = math_internal::la_pivot(a, i, i + 1);
        if (p < 0) continue;
        if (p != i + 1) {
            math_internal::la_swap_rows(a, p, i + 1);
            for (int k = 0; k < n; ++k) std::swap(a(k, p), a(k, i + 1));
        }
        const T iv = T(1) / a(i + 1, i);
        std::vector<T> c((size_t)n, T(0));  // 行 j −= c_j 行 (i+1) を全部してから、列 (i+1) += Σ c_j 列 j（E A E^{−1}）
        bool any = false;
        for (int j = i + 2; j < n; ++j) {
            c[(size_t)j] = a(j, i) * iv;
            if (math_internal::la_zero(c[(size_t)j])) continue;
            any = true;
            T* aj = a[j];
            const T* ap = a[i + 1];
            for (int k = i; k < n; ++k) aj[k] = aj[k] - c[(size_t)j] * ap[k];  // 行 (i+1) の i より左は 0
        }
        if (!any) continue;
        for (int k = 0; k < n; ++k) {
            T* ak = a[k];
            T s = ak[i + 1];
            for (int j = i + 2; j < n; ++j) s = s + c[(size_t)j] * ak[j];
            ak[i + 1] = s;
        }
    }
    // p_{i+1}(x) = (x − h_ii) p_i(x) − Σ_{j<i} h_ji (h_{j+1,j} … h_{i,i−1}) p_j(x)
    std::vector<std::vector<T>> p((size_t)n + 1);
    p[0].assign(1, T(1));
    for (int i = 0; i < n; ++i) {
        std::vector<T>& q = p[(size_t)i + 1];
        q.assign((size_t)i + 2, T(0));
        for (int k = 0; k <= i; ++k) q[(size_t)k + 1] = q[(size_t)k + 1] + p[(size_t)i][(size_t)k], q[(size_t)k] = q[(size_t)k] - a(i, i) * p[(size_t)i][(size_t)k];
        T prod = T(1);
        for (int j = i - 1; j >= 0; --j) {
            prod = prod * a(j + 1, j);
            const T c = a(j, i) * prod;
            if (math_internal::la_zero(c)) continue;
            for (int k = 0; k <= j; ++k) q[(size_t)k] = q[(size_t)k] - c * p[(size_t)j][(size_t)k];
        }
    }
    return p[(size_t)n];
}
// det(A + x B) の係数（低い次数から。長さ n + 1）。T は mint（double では誤差が大きいことがある）  O(n^3)
template <class T>
std::vector<T> det_affine(Matrix<T> A, Matrix<T> B) {
    if (A.rows() != A.cols() || B.rows() != A.rows() || B.cols() != A.cols()) math_internal::fail("det_affine", "同じ大きさの正方行列で使う");
    const int n = A.rows();
    T adj = T(1);    // もとの det = adj × （今の det）/ x^shift
    int shift = 0;
    for (int p = 0; p < n; ++p) {
        int piv = math_internal::la_pivot(B, p, p);
        while (piv < 0) {  // B の p 列目（p 行目以降）が全部 0：列の変形で上も 0 にしてから、その列に x を掛ける
            if (++shift > n) return std::vector<T>((size_t)n + 1, T(0));
            for (int r = 0; r < p; ++r) {  // B の r 列目は単位ベクトル e_r になっている
                const T v = B(r, p);
                if (math_internal::la_zero(v)) continue;
                B(r, p) = T(0);
                for (int i = 0; i < n; ++i) A(i, p) = A(i, p) - v * A(i, r);
            }
            for (int i = 0; i < n; ++i) std::swap(A(i, p), B(i, p));  // 列 p は定数だけ → x 倍すると B 側へ
            piv = math_internal::la_pivot(B, p, p);
        }
        if (piv != p) math_internal::la_swap_rows(A, piv, p), math_internal::la_swap_rows(B, piv, p), adj = -adj;
        const T v = B(p, p), iv = T(1) / v;
        adj = adj * v;
        for (int j = 0; j < n; ++j) A(p, j) = A(p, j) * iv, B(p, j) = B(p, j) * iv;
        const T *ap = A[p], *bp = B[p];
        for (int r = 0; r < n; ++r) {
            if (r == p) continue;
            const T f = B(r, p);
            if (math_internal::la_zero(f)) continue;
            T *ar = A[r], *br = B[r];
            for (int j = 0; j < n; ++j) ar[j] = ar[j] - f * ap[j];
            for (int j = p; j < n; ++j) br[j] = br[j] - f * bp[j];  // B の p より左の列は単位ベクトルなので B(p, j) = 0
        }
    }
    // B = I：det(A + x I) = det(x I − (−A))
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) A(i, j) = -A(i, j);
    std::vector<T> cp = char_poly(A), res((size_t)n + 1, T(0));
    for (int k = shift; k <= n; ++k) res[(size_t)(k - shift)] = cp[(size_t)k] * adj;
    return res;
}
#endif  // GG_LINEAR_ALGEBRA
