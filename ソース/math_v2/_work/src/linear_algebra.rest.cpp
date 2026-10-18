// =============================================================================
//  ランダムテスト
// =============================================================================
namespace linalg_test {

{COMMON_TEST}

using mint = ModInt998;
struct Q {  // テスト用の小さい有理数（「割り算ができる」型なら何でも使えることの確認。Fraction の代わり）
    Int128 p = 0, q = 1;
    Q() {}
    Q(long long v) : p(v), q(1) {}
    Q(Int128 a, Int128 b) {
        if (b < 0) a = -a, b = -b;
        Int128 g = a < 0 ? -a : a, h = b;
        while (h) {
            const Int128 t = g % h;
            g = h, h = t;
        }
        p = a / g, q = b / g;
    }
    friend Q operator+(const Q& a, const Q& b) { return Q(a.p * b.q + b.p * a.q, a.q * b.q); }
    friend Q operator-(const Q& a, const Q& b) { return Q(a.p * b.q - b.p * a.q, a.q * b.q); }
    friend Q operator*(const Q& a, const Q& b) { return Q(a.p * b.p, a.q * b.q); }
    friend Q operator/(const Q& a, const Q& b) { return Q(a.p * b.q, a.q * b.p); }
    Q operator-() const { return Q(-p, q); }
    friend bool operator==(const Q& a, const Q& b) { return a.p == b.p && a.q == b.q; }
    friend bool operator!=(const Q& a, const Q& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const Q& a) { return os << a.p << "/" << a.q; }
};

template <class T>
static Matrix<T> rand_mat(std::mt19937_64& rng, int n, int m, int lo, int hi) {
    Matrix<T> a(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j) a(i, j) = T((long long)rnd(rng, lo, hi));
    return a;
}
template <class T>
static T brute_det(const Matrix<T>& a) {  // 置換の和
    const int n = a.rows();
    std::vector<int> p((size_t)n);
    std::iota(p.begin(), p.end(), 0);
    T s = T(0);
    do {
        T t = T(1);
        for (int i = 0; i < n; ++i) t = t * a(i, p[(size_t)i]);
        int inv = 0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) inv += p[(size_t)i] > p[(size_t)j];
        s = inv % 2 ? s - t : s + t;
    } while (std::next_permutation(p.begin(), p.end()));
    return s;
}
static int brute_rank(const Matrix<Q>& a) {  // 0 でない小行列式の最大の大きさ
    const int n = a.rows(), m = a.cols();
    int best = 0;
    for (int R = 1; R < (1 << n); ++R)
        for (int C = 1; C < (1 << m); ++C) {
            const int k = __builtin_popcount((unsigned)R);
            if (k != __builtin_popcount((unsigned)C) || k <= best) continue;
            Matrix<Q> s(k, k);
            int ii = 0;
            for (int i = 0; i < n; ++i) {
                if (!(R >> i & 1)) continue;
                int jj = 0;
                for (int j = 0; j < m; ++j)
                    if (C >> j & 1) s(ii, jj++) = a(i, j);
                ++ii;
            }
            if (brute_det(s) != Q(0)) best = k;
        }
    return best;
}
template <class T>
static std::vector<T> mul_vec(const Matrix<T>& a, const std::vector<T>& x) {
    std::vector<T> r((size_t)a.rows(), T(0));
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < a.cols(); ++j) r[(size_t)i] = r[(size_t)i] + a(i, j) * x[(size_t)j];
    return r;
}

template <class T>
static void check_solve(const Matrix<T>& A, const std::vector<T>& b, int& cases) {
    ++cases;
    const LinearSolution<T> s = solve_linear(A, b);
    const int rk = mat_rank(A);
    Matrix<T> Ab(A.rows(), A.cols() + 1);
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) Ab(i, j) = A(i, j);
        Ab(i, A.cols()) = b[(size_t)i];
    }
    LINALG_CHECK(s.ok == (mat_rank(Ab) == rk));
    if (!s.ok) return;
    LINALG_CHECK(mul_vec(A, s.x) == b && (int)s.kernel.size() == A.cols() - rk);
    for (const auto& k : s.kernel) LINALG_CHECK(mul_vec(A, k) == std::vector<T>((size_t)A.rows(), T(0)));
}

static void test_exact(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 5), m = (int)rnd(rng, 0, 5), r = (int)rnd(rng, 0, 3);
        // 階数の低い行列も作る（U V）
        Matrix<Q> A = it % 2 ? rand_mat<Q>(rng, n, m, -2, 2) : rand_mat<Q>(rng, n, r, -2, 2) * rand_mat<Q>(rng, r, m, -2, 2);
        Matrix<mint> Am(n, m);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) Am(i, j) = mint((long long)A(i, j).p);
        const int rk = brute_rank(A);
        ++cases;
        LINALG_CHECK(mat_rank(A) == rk && mat_rank(Am) == rk);
        std::vector<Q> b((size_t)n);
        for (auto& x : b) x = Q(rnd(rng, -3, 3));
        if (it % 3 == 0) b = mul_vec(A, std::vector<Q>((size_t)m, Q(1)));  // 解があるもの
        check_solve(A, b, cases);
        std::vector<mint> bm((size_t)n);
        for (int i = 0; i < n; ++i) bm[(size_t)i] = mint((long long)b[(size_t)i].p);
        check_solve(Am, bm, cases);
        if (n == m) {
            const Q d = brute_det(A);
            LINALG_CHECK(mat_det(A) == d && mat_det(Am) == mint((long long)d.p));
            const Maybe<Matrix<Q>> iv = mat_inverse(A);
            LINALG_CHECK(iv.ok == (d != Q(0)));
            if (iv.ok) {
                Matrix<Q> I(n, n);
                for (int i = 0; i < n; ++i) I(i, i) = Q(1);
                LINALG_CHECK(A * iv.value == I && iv.value * A == I);
            }
        }
    }
    // 大きい mint の行列
    for (int it = 0; it < 30 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 60);
        const Matrix<mint> A = rand_mat<mint>(rng, n, n, 0, 998244352), B = rand_mat<mint>(rng, n, n, 0, 998244352);
        LINALG_CHECK(mat_det(A * B) == mat_det(A) * mat_det(B));
        const Maybe<Matrix<mint>> iv = mat_inverse(A);
        LINALG_CHECK(iv.ok && A * iv.value == Matrix<mint>::identity(n));
        std::vector<mint> b((size_t)n);
        for (auto& x : b) x = mint(rnd(rng, 0, 998244352));
        const LinearSolution<mint> s = solve_linear(A, b);
        LINALG_CHECK(s.ok && mul_vec(A, s.x) == b && s.kernel.empty());
        const int r = (int)rnd(rng, 0, n);
        const Matrix<mint> L = rand_mat<mint>(rng, n, r, 0, 998244352) * rand_mat<mint>(rng, r, n + 3, 0, 998244352);
        LINALG_CHECK(mat_rank(L) == r);
        check_solve(L, mul_vec(L, std::vector<mint>((size_t)n + 3, mint(7))), cases);
    }
}

static void test_double(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 30);
        Matrix<double> A(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) A(i, j) = (double)rnd(rng, -1000, 1000) / 100;
        std::vector<double> x((size_t)n), b;
        for (auto& v : x) v = (double)rnd(rng, -1000, 1000) / 10;
        b = mul_vec(A, x);
        const LinearSolution<double> s = solve_linear(A, b);
        LINALG_CHECK(s.ok);
        for (int i = 0; i < n; ++i) LINALG_CHECK(std::abs(s.x[(size_t)i] - x[(size_t)i]) < 1e-6);
        const Maybe<Matrix<double>> iv = mat_inverse(A);
        LINALG_CHECK(iv.ok);
        const Matrix<double> P = A * iv.value;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) LINALG_CHECK(std::abs(P(i, j) - (i == j ? 1.0 : 0.0)) < 1e-8);
        if (n <= 6) {
            const double d = brute_det(A), e = mat_det(A);
            LINALG_CHECK(std::abs(d - e) <= 1e-9 * std::max(1.0, std::abs(d)));
        }
    }
    // 特異な行列（2 行が同じ）
    Matrix<double> S = {{1, 2, 3}, {2, 4, 6}, {1, 0, 1}};
    LINALG_CHECK(mat_rank(S) == 2 && !mat_inverse(S).ok && std::abs(mat_det(S)) < 1e-12);
    ++cases;
}

static void test_det_mod_any(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 6);
        const long long m = it % 3 == 0 ? rnd(rng, 1, 30) : it % 3 == 1 ? rnd(rng, 1, 1000000000000LL) : rnd(rng, 1, 1000000000000000000LL);
        Matrix<long long> A(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) A(i, j) = it % 2 ? rnd(rng, -5, 5) : rnd_ll(rng);
        // 置換の和を mod m で
        std::vector<int> p((size_t)n);
        std::iota(p.begin(), p.end(), 0);
        long long s = 0;
        do {
            long long t = 1 % m;
            for (int i = 0; i < n; ++i) t = mod_mul(t, mod_floor(A(i, p[(size_t)i]), m), m);
            int inv = 0;
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j) inv += p[(size_t)i] > p[(size_t)j];
            s = (long long)(((Int128)s + (inv % 2 ? (Int128)m - t : (Int128)t)) % m);
        } while (std::next_permutation(p.begin(), p.end()));
        LINALG_CHECK(det_mod_any(A, m) == s);
    }
    for (int it = 0; it < 20 * TEST_SCALE; ++it, ++cases) {  // 大きい n：素数の mod は mat_det と
        const int n = (int)rnd(rng, 1, 80);
        Matrix<long long> A(n, n);
        Matrix<mint> B(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) A(i, j) = rnd_ll(rng), B(i, j) = mint(A(i, j));
        LINALG_CHECK(det_mod_any(A, 998244353) == (long long)mat_det(B).val());
    }
}

static void test_trees(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 6), E = (int)rnd(rng, 0, 9);
        std::vector<int> us((size_t)E), vs((size_t)E);
        for (int i = 0; i < E; ++i) us[(size_t)i] = (int)rnd(rng, 0, n - 1), vs[(size_t)i] = (int)rnd(rng, 0, n - 1);
        // 無向：n−1 本の辺の組を全部試す
        long long cnt = 0;
        for (int S = 0; S < (1 << E); ++S) {
            if (__builtin_popcount((unsigned)S) != n - 1) continue;
            std::vector<int> par((size_t)n);
            std::iota(par.begin(), par.end(), 0);
            std::function<int(int)> find = [&](int x) { return par[(size_t)x] == x ? x : par[(size_t)x] = find(par[(size_t)x]); };
            bool ok = true;
            for (int i = 0; i < E && ok; ++i)
                if (S >> i & 1) {
                    const int a = find(us[(size_t)i]), b = find(vs[(size_t)i]);
                    if (a == b) ok = false;
                    else par[(size_t)a] = b;
                }
            cnt += ok;
        }
        LINALG_CHECK(count_spanning_trees<mint>(n, us, vs) == mint(cnt) && count_spanning_trees<Q>(n, us, vs) == Q(cnt));
        // 有向：root 以外の各頂点に入ってくる辺を 1 本ずつ選び、親をたどって root に着くか
        const int root = (int)rnd(rng, 0, n - 1);
        long long acnt = 0;
        std::vector<int> choice((size_t)n, -1);
        std::function<void(int)> rec = [&](int v) {
            if (v == n) {
                for (int s = 0; s < n; ++s) {
                    int x = s, steps = 0;
                    while (x != root && steps <= n) x = us[(size_t)choice[(size_t)x]], ++steps;
                    if (x != root) return;
                }
                ++acnt;
                return;
            }
            if (v == root) { rec(v + 1); return; }
            for (int i = 0; i < E; ++i)
                if (vs[(size_t)i] == v && us[(size_t)i] != v) choice[(size_t)v] = i, rec(v + 1);
        };
        rec(0);
        LINALG_CHECK(count_arborescences<mint>(n, us, vs, root) == mint(acnt));
    }
    // 完全グラフ K_n の全域木は n^{n−2}（ケイリーの公式）
    for (int n = 2; n <= 40; ++n, ++cases) {
        std::vector<int> us, vs;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) us.push_back(i), vs.push_back(j);
        LINALG_CHECK(count_spanning_trees<mint>(n, us, vs) == mint(n).pow(n - 2));
        if (n <= 12) {  // double でも
            const double d = count_spanning_trees<double>(n, us, vs), e = std::pow((double)n, n - 2);
            LINALG_CHECK(std::abs(d - e) <= 1e-9 * e);
        }
    }
}

static void test_poly(std::mt19937_64& rng, int& cases) {
    auto eval = [](const std::vector<mint>& p, mint x) {
        mint s = 0;
        for (size_t i = p.size(); i-- > 0;) s = s * x + p[i];
        return s;
    };
    for (int it = 0; it < 400 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, it % 4 ? 8 : 40);
        const int kind = (int)rnd(rng, 0, 3);
        auto gen = [&](int k) {  // 0 の多い行列・階数の低い行列・普通の行列
            if (k == 0) {
                Matrix<mint> a(n, n);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        if (rnd(rng, 0, 3) == 0) a(i, j) = mint(rnd(rng, -3, 3));
                return a;
            }
            if (k == 1) {
                const int r = (int)rnd(rng, 0, n);
                return rand_mat<mint>(rng, n, r, -2, 2) * rand_mat<mint>(rng, r, n, -2, 2);
            }
            return rand_mat<mint>(rng, n, n, 0, 998244352);
        };
        const Matrix<mint> A = gen(kind), B = gen((int)rnd(rng, 0, 3));
        const std::vector<mint> cp = char_poly(A), da = det_affine(A, B);
        LINALG_CHECK((int)cp.size() == n + 1 && cp[(size_t)n] == mint(1) && (int)da.size() == n + 1);
        for (int t = 0; t < 3; ++t) {
            const mint x = mint(rnd(rng, 0, 998244352));
            Matrix<mint> M = Matrix<mint>::identity(n) * x - A, N = A + B * x;
            LINALG_CHECK(eval(cp, x) == mat_det(M) && eval(da, x) == mat_det(N));
        }
    }
    // B = 0 なら定数 det(A)、B = I なら char_poly(−A)
    const Matrix<mint> A = rand_mat<mint>(rng, 5, 5, -3, 3);
    std::vector<mint> c0 = det_affine(A, Matrix<mint>(5, 5));
    LINALG_CHECK(c0[0] == mat_det(A) && c0[1] == mint(0) && c0[5] == mint(0));
    LINALG_CHECK(det_affine(A, Matrix<mint>::identity(5)) == char_poly(A * mint(-1)));
    // double の特性多項式（小さい行列）：トレースと行列式
    Matrix<double> D = {{2, 1, 0}, {1, 3, 1}, {0, 1, 4}};
    const std::vector<double> cd = char_poly(D);
    LINALG_CHECK(std::abs(cd[2] + 9) < 1e-9 && std::abs(cd[0] + 18) < 1e-9 && std::abs(cd[3] - 1) < 1e-12);
    cases += 3;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_exact(rng, cases);
    test_double(rng, cases);
    test_det_mod_any(rng, cases);
    test_trees(rng, cases);
    test_poly(rng, cases);
    std::cout << "[linear_algebra] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[linear_algebra] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = 500;
    const Matrix<mint> A = rand_mat<mint>(rng, N, N, 0, 998244352), B = rand_mat<mint>(rng, N, N, 0, 998244352);
    auto t0 = std::chrono::steady_clock::now();
    const mint d = mat_det(A);
    const long long t_det = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const Maybe<Matrix<mint>> iv = mat_inverse(A);
    const long long t_inv = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::vector<mint> da = det_affine(A, B);
    const long long t_aff = ms(t0);
    Matrix<double> D(N, N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) D(i, j) = (double)rnd(rng, -1000, 1000);
    t0 = std::chrono::steady_clock::now();
    const LinearSolution<double> s = solve_linear(D, std::vector<double>((size_t)N, 1.0));
    const long long t_sol = ms(t0);
    Matrix<long long> L(200, 200);
    for (int i = 0; i < 200; ++i)
        for (int j = 0; j < 200; ++j) L(i, j) = rnd_ll(rng);
    t0 = std::chrono::steady_clock::now();
    const long long dm = det_mod_any(L, 1000000000000000000LL);
    const long long t_dm = ms(t0);
    std::cout << "[linear_algebra] benchmark: n = 500 の mint で mat_det " << t_det << " ms, mat_inverse " << t_inv << " ms, det_affine " << t_aff
              << " ms, double の solve_linear " << t_sol << " ms, det_mod_any（n = 200, m = 10^18）" << t_dm << " ms  (check "
              << (d.val() + iv.value(3, 4).val() + da[7].val() + (long long)s.ok + dm % 7) % 7 << ")\n";
}

}  // namespace linalg_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace linalg_demo {

{COMMON_DEMO}

static void demo() {
    using mint = ModInt998;
    banner("線形代数（linear_algebra） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  mat_rank / mat_det / mat_inverse / solve_linear（LinearSolution{ok, x, kernel}）/ det_mod_any\n"
              << "  count_spanning_trees / count_arborescences / char_poly / det_affine\n"
              << "  （この中では using mint = ModInt998; としている）\n";

    api("【1】 mat_rank(A) / mat_det(A) / mat_inverse(A)", "int / T / Maybe<Matrix<T>>", "階数 / 行列式 / 逆行列（ok が false なら逆行列なし）。T は mint・double・Fraction", "O(n^3)");
    RUN(Matrix<mint> A = {{1, 2}, {3, 4}});
    SHOW(mat_det(A), "1·4 − 2·3 = −2（mod 998244353）");
    SHOW(mat_inverse(A).value * A == Matrix<mint>::identity(2), "A^{−1} A = I");
    SHOW(mat_inverse(Matrix<double>({{2, 1}, {1, 1}})).value(0, 1), "[[2, 1], [1, 1]] の逆行列は [[1, −1], [−1, 2]]（cout << で行ごとに出せる）");
    RUN(Matrix<double> D = {{1, 2, 3}, {2, 4, 6}, {1, 0, 1}});
    SHOW(mat_rank(D), "2 行目は 1 行目の 2 倍");
    SHOW(mat_inverse(D).ok, "正則でないので逆行列なし");

    api("【2】 solve_linear(A, b)", "LinearSolution<T>{ok, x, kernel}", "A x = b の解の 1 つ x と、A y = 0 の解の基底 kernel（解は x + Σ c_i kernel_i で全部）", "O(n m min(n, m))");
    RUN(Matrix<double> P = {{2, 1}, {1, 3}});
    SHOW(solve_linear(P, {5.0, 10.0}).x, "2x + y = 5、x + 3y = 10 → (1, 3)");
    RUN(Matrix<mint> R = {{1, 1, 1}, {0, 1, 2}});
    RUN(LinearSolution<mint> s = solve_linear(R, {mint(6), mint(8)}));
    SHOW(s.ok, "解がある");
    SHOW(s.x, "解の 1 つ（自由な変数 z = 0）");
    SHOW(s.kernel, "(1, −2, 1) の倍数を足しても解（−2 は mod で 998244351）");
    SHOW(solve_linear(R, {mint(1), mint(2)}).ok, "これも解がある");
    SHOW(solve_linear(Matrix<mint>({{1, 1}, {1, 1}}), {mint(1), mint(2)}).ok, "x + y = 1 と x + y = 2 は同時に成り立たない");

    api("【3】 det_mod_any(A, m)", "long long", "Matrix<long long> の行列式を mod m で（m は合成数でもよい。割り算を使わず互除法で）", "O(n^3 log m)");
    SHOW(det_mod_any(Matrix<long long>({{2, 1}, {1, 3}}), 4), "det = 5 → 5 mod 4");
    SHOW(det_mod_any(Matrix<long long>({{6, 4}, {4, 6}}), 1000000000), "det = 20");

    api("【4】 count_spanning_trees<T>(n, us, vs) / count_arborescences<T>(n, us, vs, root)", "T",
        "無向グラフの全域木の数（行列木定理。多重辺は別々に数える）/ 有向グラフで root から全頂点へ向かう有向木の数", "O(n^3)");
    SHOW(count_spanning_trees<mint>(4, {0, 1, 2, 3, 0}, {1, 2, 3, 0, 2}), "四角形 + 対角線 1 本");
    SHOW(count_spanning_trees<mint>(5, {0, 0, 0, 0, 1, 1, 1, 2, 2, 3}, {1, 2, 3, 4, 2, 3, 4, 3, 4, 4}), "K_5：5^3 = 125");
    SHOW(count_arborescences<mint>(3, {0, 0, 1, 2}, {1, 2, 2, 1}, 0), "0→1, 0→2, 1→2, 2→1：0 から全部へ");
    note("graph_grid_v2 の辺（G.es）から us・vs を作って渡せばよい");
    {
        // ABC253 Ex（We Love Forest）：S ごとの全域木の数 → S を木に分ける数 → 辺 K 本の森 × K! / M^K
        auto solve = [](int N, const std::vector<std::pair<int, int>>& E) {
            const int M = (int)E.size(), FULL = (1 << N) - 1;
            std::vector<mint> tree((size_t)FULL + 1, mint(0));
            for (int S = 1; S <= FULL; ++S) {
                std::vector<int> id((size_t)N, -1), us, vs;
                int k = 0;
                for (int v = 0; v < N; ++v)
                    if (S >> v & 1) id[(size_t)v] = k++;
                for (const auto& e : E)
                    if (id[(size_t)e.first] >= 0 && id[(size_t)e.second] >= 0) us.push_back(id[(size_t)e.first]), vs.push_back(id[(size_t)e.second]);
                tree[(size_t)S] = count_spanning_trees<mint>(k, us, vs);
            }
            // g[c][S]：S を c 本の木に分ける数（S のいちばん下のビットを含む木を選ぶ）
            std::vector<std::vector<mint>> g((size_t)N + 1, std::vector<mint>((size_t)FULL + 1, mint(0)));
            g[0][0] = 1;
            for (int S = 1; S <= FULL; ++S) {
                const int low = S & -S;
                for (int T = S; T > 0; T = (T - 1) & S) {
                    if (!(T & low)) continue;
                    for (int c = 1; c <= N; ++c) g[(size_t)c][(size_t)S] += tree[(size_t)T] * g[(size_t)c - 1][(size_t)(S ^ T)];
                }
            }
            std::vector<mint> res;
            mint fact = 1;
            for (int K = 1; K < N; ++K) fact *= mint(K), res.push_back(g[(size_t)(N - K)][(size_t)FULL] * fact / mint(M).pow(K));
            return res;
        };
        note("ABC253 Ex（We Love Forest）: 頂点の集合 S ごとに count_spanning_trees → 森の数の DP → × K! / M^K");
        std::cout << "    例 1 → " << to_s(solve(3, {{0, 1}, {1, 2}})) << "\n    例 2 → " << to_s(solve(4, {{0, 1}, {0, 1}, {0, 3}, {1, 2}, {1, 3}})) << "\n";
    }

    api("【5】 char_poly(A) / det_affine(A, B)", "vector<T>", "det(x I − A) の係数 / det(A + x B) の係数（低い次数から。長さ n + 1。B は正則でなくてよい）", "O(n^3)");
    SHOW(char_poly(Matrix<mint>({{2, 1}, {1, 2}})), "x^2 − 4x + 3 = (x − 1)(x − 3)");
    SHOW(det_affine(Matrix<mint>({{1, 0}, {0, 1}}), Matrix<mint>({{1, 0}, {0, 0}})), "(1 + x)·1 = 1 + x");
    {
        // ABC323 G（Inversion of Tree）：転倒の辺の重みを x にしたラプラシアン A + x B の行列式の係数
        auto solve = [](const std::vector<int>& perm) {
            const int N = (int)perm.size();
            Matrix<mint> L0(N - 1, N - 1), L1(N - 1, N - 1);  // 頂点 N−1 の行と列を除く。L0 + x L1
            for (int i = 0; i < N; ++i)
                for (int j = i + 1; j < N; ++j) {
                    Matrix<mint>& L = perm[(size_t)i] > perm[(size_t)j] ? L1 : L0;
                    for (int v : {i, j})
                        if (v < N - 1) L(v, v) += mint(1);
                    if (j < N - 1) L(i, j) -= mint(1), L(j, i) -= mint(1);
                }
            return det_affine(L0, L1);  // x^K の係数 = 転倒の辺がちょうど K 本の全域木の数
        };
        note("ABC323 G（Inversion of Tree）: 転倒の辺だけ重み x のラプラシアンで det_affine（x^K の係数が答え）");
        std::cout << "    例 1 → " << to_s(solve({1, 3, 2})) << "\n    例 2 → " << to_s(solve({3, 1, 4, 10, 8, 6, 9, 2, 7, 5})) << "\n";
    }
}

}  // namespace linalg_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    linalg_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    linalg_test::benchmark();    // 速度計測
    linalg_demo::demo();         // 使い方デモ
    return 0;
}
