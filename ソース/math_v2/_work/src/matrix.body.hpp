#ifndef GG_MATRIX
#define GG_MATRIX
// ---- 行列（数学第二世代）  Matrix<T>（mint・long long・double・半環）・MatrixFixed<T, K>（小さい行列を速く）・漸化式の行列 ----
namespace math_internal {
template <class T, class = void>
struct has_semiring_one : std::false_type {};
template <class T>
struct has_semiring_one<T, std::void_t<decltype(T::one())>> : std::true_type {};
template <class T>
T sr_zero() {  // 足し算の単位元（普通の数は 0、半環は T::zero()）
    if constexpr (has_semiring_one<T>::value) return T::zero();
    else return T(0);
}
template <class T>
T sr_one() {  // 掛け算の単位元（普通の数は 1、半環は T::one()）
    if constexpr (has_semiring_one<T>::value) return T::one();
    else return T(1);
}
template <class T>
inline void mat_fma(T& acc, const T& a, const T& b) {  // acc += a * b（整数の型は溢れたら止める）
    if constexpr (std::is_integral_v<T>) {
        T p;
        if (__builtin_mul_overflow(a, b, &p) || __builtin_add_overflow(acc, p, &acc)) fail("Matrix", "値が整数の型に収まらない（数え上げは mint で）");
    } else {
        acc = acc + a * b;
    }
}
template <class T>
inline T mat_sub(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T s;
        if (__builtin_sub_overflow(a, b, &s)) fail("Matrix", "値が整数の型に収まらない（数え上げは mint で）");
        return s;
    } else {
        return a - b;
    }
}
template <class T>
inline T mat_mul(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T s;
        if (__builtin_mul_overflow(a, b, &s)) fail("Matrix", "値が整数の型に収まらない（数え上げは mint で）");
        return s;
    } else {
        return a * b;
    }
}
template <class T>
inline T mat_add(const T& a, const T& b) {
    if constexpr (std::is_integral_v<T>) {
        T s;
        if (__builtin_add_overflow(a, b, &s)) fail("Matrix", "値が整数の型に収まらない（数え上げは mint で）");
        return s;
    } else {
        return a + b;
    }
}
}  // namespace math_internal

// ---- 半環（行列の足し算・掛け算を別のものに変える） ----
template <class T>
struct MinPlus {  // 足し算 = min、掛け算 = +（ちょうど K 歩の最短路）。道がないのは INF
    static constexpr T INF = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();
    T v;
    MinPlus() : v(INF) {}
    MinPlus(T x) : v(x) {}
    static MinPlus zero() { return MinPlus(); }
    static MinPlus one() { return MinPlus(T(0)); }
    bool is_inf() const { return v == INF; }
    friend MinPlus operator+(const MinPlus& a, const MinPlus& b) { return MinPlus(a.v < b.v ? a.v : b.v); }
    friend MinPlus operator*(const MinPlus& a, const MinPlus& b) {
        if (a.v == INF || b.v == INF) return MinPlus();
        if constexpr (std::is_integral_v<T>) {
            T s;
            if (__builtin_add_overflow(a.v, b.v, &s) || s == INF) math_internal::fail("MinPlus", "長さが型に収まらない");
            return MinPlus(s);
        } else {
            return MinPlus(a.v + b.v);
        }
    }
    friend bool operator==(const MinPlus& a, const MinPlus& b) { return a.v == b.v; }
    friend bool operator!=(const MinPlus& a, const MinPlus& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const MinPlus& a) { return a.v == INF ? os << "INF" : os << a.v; }
};
template <class T>
struct MaxPlus {  // 足し算 = max、掛け算 = +（ちょうど K 歩の最長路）。道がないのは −INF
    static constexpr T NEG_INF = std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest();
    T v;
    MaxPlus() : v(NEG_INF) {}
    MaxPlus(T x) : v(x) {}
    static MaxPlus zero() { return MaxPlus(); }
    static MaxPlus one() { return MaxPlus(T(0)); }
    bool is_inf() const { return v == NEG_INF; }
    friend MaxPlus operator+(const MaxPlus& a, const MaxPlus& b) { return MaxPlus(a.v > b.v ? a.v : b.v); }
    friend MaxPlus operator*(const MaxPlus& a, const MaxPlus& b) {
        if (a.v == NEG_INF || b.v == NEG_INF) return MaxPlus();
        if constexpr (std::is_integral_v<T>) {
            T s;
            if (__builtin_add_overflow(a.v, b.v, &s) || s == NEG_INF) math_internal::fail("MaxPlus", "長さが型に収まらない");
            return MaxPlus(s);
        } else {
            return MaxPlus(a.v + b.v);
        }
    }
    friend bool operator==(const MaxPlus& a, const MaxPlus& b) { return a.v == b.v; }
    friend bool operator!=(const MaxPlus& a, const MaxPlus& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const MaxPlus& a) { return a.v == NEG_INF ? os << "-INF" : os << a.v; }
};
template <class T>
struct MinMax {  // 足し算 = min、掛け算 = max（ちょうど K 歩で、通る辺の重みの最大を最小に。ボトルネック）。道がないのは INF
    static constexpr T INF = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();
    static constexpr T NEG_INF = std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest();
    T v;
    MinMax() : v(INF) {}
    MinMax(T x) : v(x) {}
    static MinMax zero() { return MinMax(); }         // 道がない
    static MinMax one() { return MinMax(NEG_INF); }   // 0 歩（辺を通らない）
    bool is_inf() const { return v == INF; }
    friend MinMax operator+(const MinMax& a, const MinMax& b) { return MinMax(a.v < b.v ? a.v : b.v); }
    friend MinMax operator*(const MinMax& a, const MinMax& b) { return MinMax(a.v > b.v ? a.v : b.v); }
    friend bool operator==(const MinMax& a, const MinMax& b) { return a.v == b.v; }
    friend bool operator!=(const MinMax& a, const MinMax& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const MinMax& a) { return a.v == INF ? os << "INF" : (a.v == NEG_INF ? os << "-INF" : os << a.v); }
};
template <class T>
struct MaxMin {  // 足し算 = max、掛け算 = min（ちょうど K 歩で、通る辺の重みの最小を最大に。最大容量路）。道がないのは −INF
    static constexpr T INF = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max();
    static constexpr T NEG_INF = std::numeric_limits<T>::has_infinity ? -std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest();
    T v;
    MaxMin() : v(NEG_INF) {}
    MaxMin(T x) : v(x) {}
    static MaxMin zero() { return MaxMin(); }      // 道がない
    static MaxMin one() { return MaxMin(INF); }    // 0 歩
    bool is_inf() const { return v == NEG_INF; }
    friend MaxMin operator+(const MaxMin& a, const MaxMin& b) { return MaxMin(a.v > b.v ? a.v : b.v); }
    friend MaxMin operator*(const MaxMin& a, const MaxMin& b) { return MaxMin(a.v < b.v ? a.v : b.v); }
    friend bool operator==(const MaxMin& a, const MaxMin& b) { return a.v == b.v; }
    friend bool operator!=(const MaxMin& a, const MaxMin& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const MaxMin& a) { return a.v == NEG_INF ? os << "-INF" : (a.v == INF ? os << "INF" : os << a.v); }
};
struct XorAnd {  // 足し算 = XOR、掛け算 = AND（ビットごとの漸化式。ABC009 D）
    unsigned long long v;
    XorAnd(unsigned long long x = 0) : v(x) {}
    static XorAnd zero() { return XorAnd(0); }
    static XorAnd one() { return XorAnd(~0ULL); }
    friend XorAnd operator+(const XorAnd& a, const XorAnd& b) { return XorAnd(a.v ^ b.v); }
    friend XorAnd operator*(const XorAnd& a, const XorAnd& b) { return XorAnd(a.v & b.v); }
    friend bool operator==(const XorAnd& a, const XorAnd& b) { return a.v == b.v; }
    friend bool operator!=(const XorAnd& a, const XorAnd& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const XorAnd& a) { return os << a.v; }
};
struct OrAnd {  // 足し算 = OR、掛け算 = AND（ちょうど K 歩で行けるか）
    bool v;
    OrAnd(bool x = false) : v(x) {}
    static OrAnd zero() { return OrAnd(false); }
    static OrAnd one() { return OrAnd(true); }
    friend OrAnd operator+(const OrAnd& a, const OrAnd& b) { return OrAnd(a.v || b.v); }
    friend OrAnd operator*(const OrAnd& a, const OrAnd& b) { return OrAnd(a.v && b.v); }
    friend bool operator==(const OrAnd& a, const OrAnd& b) { return a.v == b.v; }
    friend bool operator!=(const OrAnd& a, const OrAnd& b) { return a.v != b.v; }
    friend std::ostream& operator<<(std::ostream& os, const OrAnd& a) { return os << (a.v ? 1 : 0); }
};

// ---- 行列 Matrix<T>（n 行 m 列。T は mint・long long・double・上の半環） ----
template <class T>
struct Matrix {
    Matrix() : n_(0), m_(0) {}
    Matrix(int n, int m) : n_(n), m_(m) {  // すべて 0（半環なら zero）
        if (n < 0 || m < 0) math_internal::fail("Matrix", "大きさは 0 以上");
        a_.assign((size_t)n * (size_t)m, math_internal::sr_zero<T>());
    }
    explicit Matrix(int n) : Matrix(n, n) {}
    Matrix(const std::vector<std::vector<T>>& v) : Matrix((int)v.size(), v.empty() ? 0 : (int)v[0].size()) {
        for (int i = 0; i < n_; ++i) {
            if ((int)v[(size_t)i].size() != m_) math_internal::fail("Matrix", "行の長さがそろっていない");
            for (int j = 0; j < m_; ++j) (*this)(i, j) = v[(size_t)i][(size_t)j];
        }
    }
    Matrix(std::initializer_list<std::initializer_list<T>> v) : Matrix(std::vector<std::vector<T>>(v.begin(), v.end())) {}
    static Matrix identity(int n) {  // 単位行列
        Matrix r(n, n);
        for (int i = 0; i < n; ++i) r(i, i) = math_internal::sr_one<T>();
        return r;
    }
    static Matrix zero(int n, int m) { return Matrix(n, m); }
    int rows() const { return n_; }
    int cols() const { return m_; }
    T& operator()(int i, int j) { return a_[(size_t)i * (size_t)m_ + (size_t)j]; }
    const T& operator()(int i, int j) const { return a_[(size_t)i * (size_t)m_ + (size_t)j]; }
    T* operator[](int i) { return a_.data() + (size_t)i * (size_t)m_; }  // A[i][j] とも書ける
    const T* operator[](int i) const { return a_.data() + (size_t)i * (size_t)m_; }

    Matrix& operator+=(const Matrix& b) {
        same_size(b, "Matrix::+");
        for (size_t i = 0; i < a_.size(); ++i) a_[i] = math_internal::mat_add(a_[i], b.a_[i]);
        return *this;
    }
    Matrix& operator-=(const Matrix& b) {  // 半環では使えない
        same_size(b, "Matrix::-");
        for (size_t i = 0; i < a_.size(); ++i) a_[i] = math_internal::mat_sub(a_[i], b.a_[i]);
        return *this;
    }
    Matrix& operator*=(const Matrix& b) { return *this = *this * b; }
    Matrix& operator*=(const T& s) {  // スカラー倍
        for (T& x : a_) x = math_internal::mat_mul(x, s);
        return *this;
    }
    friend Matrix operator+(Matrix a, const Matrix& b) { return a += b; }
    friend Matrix operator-(Matrix a, const Matrix& b) { return a -= b; }
    friend Matrix operator*(Matrix a, const T& s) { return a *= s; }
    friend Matrix operator*(const T& s, Matrix a) { return a *= s; }
    friend Matrix operator*(const Matrix& a, const Matrix& b) {  // O(n m k)
        if (a.m_ != b.n_) math_internal::fail("Matrix::*", "左の列数と右の行数が違う");
        const T z = math_internal::sr_zero<T>();
        Matrix c(a.n_, b.m_);
        for (int i = 0; i < a.n_; ++i) {
            T* ci = c[i];
            for (int k = 0; k < a.m_; ++k) {
                const T aik = a(i, k);
                if (aik == z) continue;
                const T* bk = b[k];
                for (int j = 0; j < b.m_; ++j) math_internal::mat_fma(ci[j], aik, bk[j]);
            }
        }
        return c;
    }
    Matrix pow(long long e) const {  // A^e（e ≥ 0。e ≤ 10^18 でよい）  O(n^3 log e)
        if (n_ != m_) math_internal::fail("Matrix::pow", "正方行列で使う");
        if (e < 0) math_internal::fail("Matrix::pow", "e ≥ 0 で使う");
        Matrix r = identity(n_), b = *this;
        for (; e > 0; e >>= 1) {
            if (e & 1) r = r * b;
            if (e > 1) b = b * b;
        }
        return r;
    }
    std::vector<T> apply(const std::vector<T>& v) const {  // A v（行列 × 縦ベクトル）  O(n m)
        if ((int)v.size() != m_) math_internal::fail("Matrix::apply", "ベクトルの長さが列数と違う");
        std::vector<T> r((size_t)n_, math_internal::sr_zero<T>());
        for (int i = 0; i < n_; ++i)
            for (int j = 0; j < m_; ++j) math_internal::mat_fma(r[(size_t)i], (*this)(i, j), v[(size_t)j]);
        return r;
    }
    Matrix transpose() const {
        Matrix r(m_, n_);
        for (int i = 0; i < n_; ++i)
            for (int j = 0; j < m_; ++j) r(j, i) = (*this)(i, j);
        return r;
    }
    T trace() const {  // 対角成分の和
        if (n_ != m_) math_internal::fail("Matrix::trace", "正方行列で使う");
        T s = math_internal::sr_zero<T>();
        for (int i = 0; i < n_; ++i) s = math_internal::mat_add(s, (*this)(i, i));
        return s;
    }
    friend bool operator==(const Matrix& a, const Matrix& b) { return a.n_ == b.n_ && a.m_ == b.m_ && a.a_ == b.a_; }
    friend bool operator!=(const Matrix& a, const Matrix& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const Matrix& a) {  // 1 行ずつ（行の中は空白区切り）
        for (int i = 0; i < a.n_; ++i) {
            for (int j = 0; j < a.m_; ++j) os << (j ? " " : "") << a(i, j);
            if (i + 1 < a.n_) os << '\n';
        }
        return os;
    }

  private:
    int n_, m_;
    std::vector<T> a_;
    void same_size(const Matrix& b, const char* f) const {
        if (n_ != b.n_ || m_ != b.m_) math_internal::fail(f, "大きさが違う");
    }
};

// ---- 線形漸化式の行列 ----
// a_i = c_0 a_{i−1} + c_1 a_{i−2} + … + c_{k−1} a_{i−k} の行列（1 行目が c、その下が 1 つずらす形）
//   (a_{i}, a_{i−1}, …, a_{i−k+1})^T = M (a_{i−1}, …, a_{i−k})^T
template <class T>
Matrix<T> companion_matrix(const std::vector<T>& c) {
    const int k = (int)c.size();
    if (k == 0) math_internal::fail("companion_matrix", "c が空");
    Matrix<T> m(k, k);
    for (int j = 0; j < k; ++j) m(0, j) = c[(size_t)j];
    for (int i = 1; i < k; ++i) m(i, i - 1) = math_internal::sr_one<T>();
    return m;
}
// 上の漸化式で a_0〜a_{k−1} が a のとき a_n（n は 0 から数える。n ≤ 10^18。半環でもよい）  O(k^3 log n)
template <class T>
T recurrence_nth_matrix(const std::vector<T>& c, const std::vector<T>& a, long long n) {
    const int k = (int)c.size();
    if ((int)a.size() != k) math_internal::fail("recurrence_nth_matrix", "c と a の長さが違う");
    if (n < 0) math_internal::fail("recurrence_nth_matrix", "n ≥ 0 で使う");
    if (n < k) return a[(size_t)n];
    std::vector<T> s(a.rbegin(), a.rend());  // (a_{k−1}, …, a_0)
    return companion_matrix(c).pow(n - k + 1).apply(s)[0];
}

// ---- 小さい行列を速く MatrixFixed<T, K>（K×K。大きさがコンパイル時に決まるので展開される） ----
template <class T, int K>
struct MatrixFixed {
    std::array<std::array<T, K>, K> a;
    MatrixFixed() {
        for (auto& row : a) row.fill(math_internal::sr_zero<T>());
    }
    MatrixFixed(std::initializer_list<std::initializer_list<T>> v) : MatrixFixed() {
        if ((int)v.size() != K) math_internal::fail("MatrixFixed", "行の数が K でない");
        int i = 0;
        for (const auto& row : v) {
            if ((int)row.size() != K) math_internal::fail("MatrixFixed", "列の数が K でない");
            int j = 0;
            for (const T& x : row) a[(size_t)i][(size_t)j++] = x;
            ++i;
        }
    }
    static MatrixFixed identity() {
        MatrixFixed r;
        for (int i = 0; i < K; ++i) r.a[(size_t)i][(size_t)i] = math_internal::sr_one<T>();
        return r;
    }
    T& operator()(int i, int j) { return a[(size_t)i][(size_t)j]; }
    const T& operator()(int i, int j) const { return a[(size_t)i][(size_t)j]; }
    std::array<T, K>& operator[](int i) { return a[(size_t)i]; }
    const std::array<T, K>& operator[](int i) const { return a[(size_t)i]; }
    friend MatrixFixed operator+(MatrixFixed x, const MatrixFixed& y) {
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) x.a[i][j] = math_internal::mat_add(x.a[i][j], y.a[i][j]);
        return x;
    }
    friend MatrixFixed operator-(MatrixFixed x, const MatrixFixed& y) {  // 半環では使えない
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) x.a[i][j] = math_internal::mat_sub(x.a[i][j], y.a[i][j]);
        return x;
    }
    friend MatrixFixed operator*(const MatrixFixed& x, const MatrixFixed& y) {
        MatrixFixed r;
        for (int i = 0; i < K; ++i)
            for (int k = 0; k < K; ++k)
                for (int j = 0; j < K; ++j) math_internal::mat_fma(r.a[i][j], x.a[i][k], y.a[k][j]);
        return r;
    }
    friend MatrixFixed operator*(MatrixFixed x, const T& s) {  // スカラー倍
        for (auto& row : x.a)
            for (T& v : row) v = math_internal::mat_mul(v, s);
        return x;
    }
    friend MatrixFixed operator*(const T& s, const MatrixFixed& x) { return x * s; }
    MatrixFixed& operator*=(const T& s) { return *this = *this * s; }
    MatrixFixed& operator+=(const MatrixFixed& y) { return *this = *this + y; }
    MatrixFixed& operator-=(const MatrixFixed& y) { return *this = *this - y; }
    MatrixFixed& operator*=(const MatrixFixed& y) { return *this = *this * y; }
    MatrixFixed pow(long long e) const {
        if (e < 0) math_internal::fail("MatrixFixed::pow", "e ≥ 0 で使う");
        MatrixFixed r = identity(), b = *this;
        for (; e > 0; e >>= 1) {
            if (e & 1) r = r * b;
            if (e > 1) b = b * b;
        }
        return r;
    }
    std::array<T, K> apply(const std::array<T, K>& v) const {
        std::array<T, K> r;
        r.fill(math_internal::sr_zero<T>());
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) math_internal::mat_fma(r[(size_t)i], a[i][j], v[(size_t)j]);
        return r;
    }
    MatrixFixed transpose() const {
        MatrixFixed r;
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) r.a[j][i] = a[i][j];
        return r;
    }
    T trace() const {
        T s = math_internal::sr_zero<T>();
        for (int i = 0; i < K; ++i) s = math_internal::mat_add(s, a[i][i]);
        return s;
    }
    T det() const {  // 行列式（T は mint（mod は素数）か double）  O(K^3)
        static_assert(!std::is_integral_v<T>, "det は mint か double で使う");
        MatrixFixed b = *this;
        T d = T(1);
        for (int c = 0; c < K; ++c) {
            const int p = pivot(b, c);
            if (p < 0) return T(0);
            if (p != c) std::swap(b.a[p], b.a[c]), d = -d;
            d = d * b.a[c][c];
            const T iv = T(1) / b.a[c][c];
            for (int r = c + 1; r < K; ++r) {
                const T f = b.a[r][c] * iv;
                for (int j = c; j < K; ++j) b.a[r][j] = b.a[r][j] - f * b.a[c][j];
            }
        }
        return d;
    }
    Maybe<MatrixFixed> try_inv() const {  // 逆行列（なければ ok = false。T は mint（mod は素数）か double）
        static_assert(!std::is_integral_v<T>, "inv は mint か double で使う");
        MatrixFixed b = *this, r = identity();
        for (int c = 0; c < K; ++c) {
            const int p = pivot(b, c);
            if (p < 0) return Maybe<MatrixFixed>{false, MatrixFixed()};
            std::swap(b.a[p], b.a[c]), std::swap(r.a[p], r.a[c]);
            const T iv = T(1) / b.a[c][c];
            for (int j = 0; j < K; ++j) b.a[c][j] = b.a[c][j] * iv, r.a[c][j] = r.a[c][j] * iv;
            for (int i = 0; i < K; ++i) {
                if (i == c) continue;
                const T f = b.a[i][c];
                for (int j = 0; j < K; ++j) b.a[i][j] = b.a[i][j] - f * b.a[c][j], r.a[i][j] = r.a[i][j] - f * r.a[c][j];
            }
        }
        return Maybe<MatrixFixed>{true, r};
    }
    MatrixFixed inv() const {  // 逆行列（なければ止める）
        const Maybe<MatrixFixed> r = try_inv();
        if (!r.ok) math_internal::fail("MatrixFixed::inv", "逆行列がない（行列式が 0）");
        return r.value;
    }
    friend bool operator==(const MatrixFixed& x, const MatrixFixed& y) { return x.a == y.a; }
    friend bool operator!=(const MatrixFixed& x, const MatrixFixed& y) { return !(x == y); }
    friend std::ostream& operator<<(std::ostream& os, const MatrixFixed& x) {
        for (int i = 0; i < K; ++i) {
            for (int j = 0; j < K; ++j) os << (j ? " " : "") << x.a[i][j];
            if (i + 1 < K) os << '\n';
        }
        return os;
    }

  private:
    static int pivot(const MatrixFixed& b, int c) {  // c 列目で c 行目以降のピボット（double は絶対値が最大のもの。なければ −1）
        if constexpr (std::is_floating_point_v<T>) {
            int p = c;
            for (int r = c + 1; r < K; ++r)
                if (std::abs(b.a[r][c]) > std::abs(b.a[p][c])) p = r;
            return std::abs(b.a[p][c]) < 1e-12 ? -1 : p;
        } else {
            for (int r = c; r < K; ++r)
                if (b.a[r][c] != T(0)) return r;
            return -1;
        }
    }
};
#endif  // GG_MATRIX
