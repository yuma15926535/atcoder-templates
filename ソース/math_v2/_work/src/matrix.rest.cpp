// =============================================================================
//  ランダムテスト
// =============================================================================
namespace mat_test {

{COMMON_TEST}

using mint = ModInt998;
template <class T>
static Matrix<T> naive_mul(const Matrix<T>& a, const Matrix<T>& b) {  // 愚直（i, j, k の順）
    Matrix<T> c(a.rows(), b.cols());
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < b.cols(); ++j) {
            T s = math_internal::sr_zero<T>();
            for (int k = 0; k < a.cols(); ++k) s = s + a(i, k) * b(k, j);
            c(i, j) = s;
        }
    return c;
}
static Matrix<mint> rand_mat(std::mt19937_64& rng, int n, int m) {
    Matrix<mint> a(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j) a(i, j) = rng() % 4 == 0 ? mint(0) : mint((long long)rng());
    return a;
}

static void test_mint(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 0, 6), m = (int)rnd(rng, 0, 6), k = (int)rnd(rng, 0, 6);
        const Matrix<mint> A = rand_mat(rng, n, m), B = rand_mat(rng, m, k), C = rand_mat(rng, n, m);
        MAT_CHECK(A * B == naive_mul(A, B));
        MAT_CHECK((A * B).transpose() == B.transpose() * A.transpose());
        const mint s = (long long)rng();
        Matrix<mint> D = A + C, E = A - C, F = A * s, G = s * A;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) MAT_CHECK(D(i, j) == A(i, j) + C(i, j) && E[i][j] == A[i][j] - C[i][j] && F(i, j) == A(i, j) * s && G(i, j) == F(i, j));
        // apply は列ベクトルとの掛け算
        std::vector<mint> v((size_t)m);
        Matrix<mint> V(m, 1);
        for (int j = 0; j < m; ++j) v[(size_t)j] = (long long)rng(), V(j, 0) = v[(size_t)j];
        const std::vector<mint> av = A.apply(v);
        const Matrix<mint> AV = A * V;
        for (int i = 0; i < n; ++i) MAT_CHECK(av[(size_t)i] == AV(i, 0));
        // 累乗
        const Matrix<mint> S = rand_mat(rng, n, n);
        Matrix<mint> P = Matrix<mint>::identity(n);
        const int e = (int)rnd(rng, 0, 12);
        for (int t = 0; t < e; ++t) P = P * S;
        MAT_CHECK(S.pow(e) == P);
        const long long e1 = rnd(rng, 0, 1000000000000000000LL), e2 = rnd(rng, 0, 1000000000000000000LL);
        if (it % 8 == 0) MAT_CHECK(S.pow(e1) * S.pow(e2) == S.pow(e1 + e2));
        mint tr = 0;
        for (int i = 0; i < n; ++i) tr += S(i, i);
        MAT_CHECK(S.trace() == tr);
    }
    MAT_CHECK(Matrix<mint>::identity(3) * Matrix<mint>({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}) == Matrix<mint>({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}));
}

static void test_integer_and_real(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 20000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 5), m = (int)rnd(rng, 1, 5), k = (int)rnd(rng, 1, 5);
        Matrix<long long> A(n, m), B(m, k);
        Matrix<double> Ad(n, m), Bd(m, k);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) A(i, j) = rnd(rng, -1000000, 1000000), Ad(i, j) = (double)A(i, j);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < k; ++j) B(i, j) = rnd(rng, -1000000, 1000000), Bd(i, j) = (double)B(i, j);
        const Matrix<long long> C = A * B;
        const Matrix<double> Cd = Ad * Bd;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < k; ++j) {
                Int128 s = 0;
                for (int t = 0; t < m; ++t) s += (Int128)A(i, t) * B(t, j);
                MAT_CHECK((Int128)C(i, j) == s && std::abs(Cd(i, j) - (double)s) <= 1e-6 * (1 + std::abs((double)s)));
            }
    }
    // フィボナッチ：long long で F_90 まで厳密に
    const Matrix<long long> F = Matrix<long long>({{1, 1}, {1, 0}}).pow(90);
    MAT_CHECK(F(0, 1) == 2880067194370816120LL);
}

template <class S>
static void test_semiring_paths(std::mt19937_64& rng, int& cases, bool minimize) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 6), K = (int)rnd(rng, 0, 12);
        Matrix<S> A(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (rng() % 3) A(i, j) = S(rnd(rng, -20, 100));
        // DP：ちょうど t 歩
        const long long NONE = std::numeric_limits<long long>::min();
        std::vector<std::vector<long long>> d((size_t)n, std::vector<long long>((size_t)n, NONE));
        for (int i = 0; i < n; ++i) d[(size_t)i][(size_t)i] = 0;
        for (int t = 0; t < K; ++t) {
            std::vector<std::vector<long long>> nd((size_t)n, std::vector<long long>((size_t)n, NONE));
            for (int i = 0; i < n; ++i)
                for (int u = 0; u < n; ++u) {
                    if (d[(size_t)i][(size_t)u] == NONE) continue;
                    for (int v = 0; v < n; ++v) {
                        if (A(u, v).is_inf()) continue;
                        const long long c = d[(size_t)i][(size_t)u] + A(u, v).v;
                        long long& x = nd[(size_t)i][(size_t)v];
                        if (x == NONE || (minimize ? c < x : c > x)) x = c;
                    }
                }
            d = nd;
        }
        const Matrix<S> P = A.pow(K);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) MAT_CHECK(d[(size_t)i][(size_t)j] == NONE ? P(i, j).is_inf() : (!P(i, j).is_inf() && P(i, j).v == d[(size_t)i][(size_t)j]));
    }
}

template <class S>
static void test_bottleneck(std::mt19937_64& rng, int& cases, bool min_of_max) {  // MinMax / MaxMin
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = (int)rnd(rng, 1, 6), K = (int)rnd(rng, 0, 12);
        Matrix<S> A(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (rng() % 3) A(i, j) = S(rnd(rng, -50, 50));
        const long long NONE = LLONG_MIN, START = min_of_max ? LLONG_MIN + 1 : LLONG_MAX;  // 0 歩は「辺なし」
        std::vector<std::vector<long long>> d((size_t)n, std::vector<long long>((size_t)n, NONE));
        for (int i = 0; i < n; ++i) d[(size_t)i][(size_t)i] = START;
        for (int t = 0; t < K; ++t) {
            std::vector<std::vector<long long>> nd((size_t)n, std::vector<long long>((size_t)n, NONE));
            for (int i = 0; i < n; ++i)
                for (int u = 0; u < n; ++u) {
                    if (d[(size_t)i][(size_t)u] == NONE) continue;
                    for (int v = 0; v < n; ++v) {
                        if (A(u, v).is_inf()) continue;
                        const long long w = A(u, v).v, c = min_of_max ? std::max(d[(size_t)i][(size_t)u], w) : std::min(d[(size_t)i][(size_t)u], w);
                        long long& x = nd[(size_t)i][(size_t)v];
                        if (x == NONE || (min_of_max ? c < x : c > x)) x = c;
                    }
                }
            d = nd;
        }
        const Matrix<S> P = A.pow(K);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                if (d[(size_t)i][(size_t)j] == NONE) MAT_CHECK(P(i, j).is_inf());
                else if (d[(size_t)i][(size_t)j] == START) MAT_CHECK(P(i, j) == S::one());
                else MAT_CHECK(!P(i, j).is_inf() && P(i, j).v == d[(size_t)i][(size_t)j]);
            }
    }
}

static void test_bool_and_xor(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {  // OrAnd：ちょうど K 歩で行けるか
        const int n = (int)rnd(rng, 1, 7), K = (int)rnd(rng, 0, 15);
        Matrix<OrAnd> A(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) A(i, j) = rng() % 4 == 0;
        std::vector<std::vector<bool>> r((size_t)n, std::vector<bool>((size_t)n, false));
        for (int i = 0; i < n; ++i) r[(size_t)i][(size_t)i] = true;
        for (int t = 0; t < K; ++t) {
            std::vector<std::vector<bool>> nr((size_t)n, std::vector<bool>((size_t)n, false));
            for (int i = 0; i < n; ++i)
                for (int u = 0; u < n; ++u)
                    if (r[(size_t)i][(size_t)u])
                        for (int v = 0; v < n; ++v)
                            if (A(u, v).v) nr[(size_t)i][(size_t)v] = true;
            r = nr;
        }
        const Matrix<OrAnd> P = A.pow(K);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) MAT_CHECK(P(i, j).v == r[(size_t)i][(size_t)j]);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // XorAnd と mint の漸化式：直接計算と比べる
        const int k = (int)rnd(rng, 1, 6);
        std::vector<XorAnd> c((size_t)k), a((size_t)k);
        std::vector<mint> cm((size_t)k), am((size_t)k);
        for (int i = 0; i < k; ++i) c[(size_t)i] = rng(), a[(size_t)i] = rng(), cm[(size_t)i] = (long long)rng(), am[(size_t)i] = (long long)rng();
        std::vector<XorAnd> seq(a);
        std::vector<mint> sm(am);
        const int N = 150;
        while ((int)seq.size() <= N) {
            XorAnd x = 0;
            mint y = 0;
            const size_t i = seq.size();
            for (int j = 0; j < k; ++j) x = x + c[(size_t)j] * seq[i - 1 - (size_t)j], y += cm[(size_t)j] * sm[i - 1 - (size_t)j];
            seq.push_back(x), sm.push_back(y);
        }
        const long long n = rnd(rng, 0, N);
        MAT_CHECK(recurrence_nth_matrix(c, a, n) == seq[(size_t)n] && recurrence_nth_matrix(cm, am, n) == sm[(size_t)n]);
    }
}

template <int K>
static void test_fixed(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 5000 * TEST_SCALE; ++it, ++cases) {
        MatrixFixed<mint, K> A, B;
        Matrix<mint> a(K, K), b(K, K);
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) {
                const long long x = it % 3 == 0 ? rnd(rng, 0, 2) : (long long)rng(), y = (long long)rng();
                A(i, j) = x, a(i, j) = x, B[i][j] = y, b(i, j) = y;
            }
        const MatrixFixed<mint, K> C = A * B, D = A + B, E = A - B, P = A.pow(it % 50);
        const Matrix<mint> c = a * b, p = a.pow(it % 50);
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) MAT_CHECK(C(i, j) == c(i, j) && D(i, j) == a(i, j) + b(i, j) && E(i, j) == a(i, j) - b(i, j) && P(i, j) == p(i, j));
        // 行列式：置換で展開
        std::vector<int> perm(K);
        std::iota(perm.begin(), perm.end(), 0);
        mint det = 0;
        do {
            int inv = 0;
            for (int i = 0; i < K; ++i)
                for (int j = i + 1; j < K; ++j) inv += perm[(size_t)i] > perm[(size_t)j];
            mint t = inv % 2 ? mint(-1) : mint(1);
            for (int i = 0; i < K; ++i) t *= A(i, perm[(size_t)i]);
            det += t;
        } while (std::next_permutation(perm.begin(), perm.end()));
        MAT_CHECK(A.det() == det && (A * B).det() == A.det() * B.det());
        const Maybe<MatrixFixed<mint, K>> iv = A.try_inv();
        MAT_CHECK(iv.ok == (det != mint(0)));
        if (iv.ok) MAT_CHECK(A * iv.value == (MatrixFixed<mint, K>::identity()) && iv.value * A == (MatrixFixed<mint, K>::identity()));
        std::array<mint, K> v;
        for (auto& x : v) x = (long long)rng();
        const std::array<mint, K> av = A.apply(v);
        for (int i = 0; i < K; ++i) {
            mint s = 0;
            for (int j = 0; j < K; ++j) s += A(i, j) * v[(size_t)j];
            MAT_CHECK(av[(size_t)i] == s);
        }
        MAT_CHECK(A.transpose().det() == det && A.trace() == a.trace());
        const mint sc = (long long)rng();
        MatrixFixed<mint, K> S1 = A * sc, S2 = sc * A, S3 = A;
        S3 *= sc;
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) MAT_CHECK(S1(i, j) == A(i, j) * sc && S2 == S1 && S3 == S1);
    }
    // double
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        MatrixFixed<double, K> A;
        for (int i = 0; i < K; ++i)
            for (int j = 0; j < K; ++j) A(i, j) = (double)rnd(rng, -100, 100) / 7;
        const double d = A.det();
        const Maybe<MatrixFixed<double, K>> iv = A.try_inv();
        if (std::abs(d) > 1e-6) {
            MAT_CHECK(iv.ok);
            const MatrixFixed<double, K> I = A * iv.value;
            for (int i = 0; i < K; ++i)
                for (int j = 0; j < K; ++j) MAT_CHECK(std::abs(I(i, j) - (i == j ? 1.0 : 0.0)) < 1e-8);
            MAT_CHECK(std::abs(iv.value.det() * d - 1) < 1e-8);
        }
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_mint(rng, cases);
    test_integer_and_real(rng, cases);
    test_semiring_paths<MinPlus<long long>>(rng, cases, true);
    test_semiring_paths<MaxPlus<long long>>(rng, cases, false);
    test_bottleneck<MinMax<long long>>(rng, cases, true);
    test_bottleneck<MaxMin<long long>>(rng, cases, false);
    test_bool_and_xor(rng, cases);
    test_fixed<1>(rng, cases);
    test_fixed<2>(rng, cases);
    test_fixed<3>(rng, cases);
    test_fixed<4>(rng, cases);
    std::cout << "[matrix] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[matrix] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    using M107 = ModInt107;
    Matrix<M107> A(300, 300), B(300, 300);
    for (int i = 0; i < 300; ++i)
        for (int j = 0; j < 300; ++j) A(i, j) = (long long)rng(), B(i, j) = (long long)rng();
    auto t0 = std::chrono::steady_clock::now();
    const Matrix<M107> C = A * B;
    const long long t_mul = ms(t0);
    Matrix<M107> G(50, 50);
    for (int i = 0; i < 50; ++i)
        for (int j = 0; j < 50; ++j) G(i, j) = (long long)(rng() % 2);
    t0 = std::chrono::steady_clock::now();
    const Matrix<M107> GP = G.pow(1000000000000000000LL);
    const long long t_pow = ms(t0);
    t0 = std::chrono::steady_clock::now();
    MatrixFixed<M107, 2> F = {{1, 1}, {1, 0}};
    M107 s = 0;
    for (int i = 0; i < 100000; ++i) s += F.pow((long long)(rng() >> 4))(0, 1);
    const long long t_fixed = ms(t0);
    Matrix<MinPlus<long long>> D(100, 100);
    for (int i = 0; i < 100; ++i)
        for (int j = 0; j < 100; ++j)
            if (rng() % 2) D(i, j) = rnd(rng, 1, 1000000);
    t0 = std::chrono::steady_clock::now();
    const Matrix<MinPlus<long long>> DP = D.pow(1000);
    const long long t_mp = ms(t0);
    std::cout << "[matrix] benchmark: Matrix<mint> 300×300 の掛け算 " << t_mul << " ms, 50×50 の 10^18 乗 " << t_pow << " ms, MatrixFixed<mint, 2> の累乗 × 10^5 " << t_fixed
              << " ms, Matrix<MinPlus> 100×100 の 1000 乗 " << t_mp << " ms  (check " << (C(1, 2) + GP(3, 4) + s).val() % 7 + DP(0, 0).v % 7 << ")\n";
}

}  // namespace mat_test


// デモで行列を 1 行に表示する（行を " / " でつなぐ。テンプレートの cout << A は 1 行ずつ改行）
template <class T>
static std::string to_s(const Matrix<T>& a) {
    std::ostringstream os;
    os << "[";
    for (int i = 0; i < a.rows(); ++i) {
        os << (i ? " / " : "");
        for (int j = 0; j < a.cols(); ++j) os << (j ? " " : "") << a(i, j);
    }
    return os.str() + "]";
}
template <class T, int K>
static std::string to_s(const MatrixFixed<T, K>& a) {
    std::ostringstream os;
    os << "[";
    for (int i = 0; i < K; ++i) {
        os << (i ? " / " : "");
        for (int j = 0; j < K; ++j) os << (j ? " " : "") << a(i, j);
    }
    return os.str() + "]";
}

// =============================================================================
//  使い方デモ
// =============================================================================
namespace mat_demo {

{COMMON_DEMO}

static void demo() {
    banner("行列（matrix） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Matrix<T>：A(i, j) / A[i][j] / + - * / A.pow(e) / A.apply(v) / A.transpose() / A.trace() / identity(n)\n"
              << "  半環：MinPlus<T> / MaxPlus<T> / MinMax<T> / MaxMin<T> / XorAnd / OrAnd（Matrix<MinPlus<long long>> のように）\n"
              << "  companion_matrix(c) / recurrence_nth_matrix(c, a, n)    線形漸化式\n"
              << "  MatrixFixed<T, K>：K×K の小さい行列（det / inv も）\n";

    api("【1】 Matrix<T> A = {{…}, {…}} / A * B / A.pow(e) / A.apply(v)", "Matrix<T> / vector<T>",
        "行列の計算。T は mint・long long（溢れたら止める）・double など", "掛け算 O(n^3)・累乗 O(n^3 log e)");
    note("using mint = ModInt107;");
    using mint = ModInt107;
    RUN(Matrix<mint> A = {{1, 1}, {1, 0}});
    SHOW(A.pow(10), "フィボナッチ：右上が F_10 = 55");
    SHOW(A.pow(1000000000000000000LL)(0, 1), "F_{10^18} mod 10^9+7");
    SHOW(A * A, "掛け算");
    SHOW((A.apply({mint(1), mint(0)})), "A × (1, 0)^T");
    SHOW(Matrix<long long>({{1, 2, 3}, {4, 5, 6}}).transpose(), "転置（2×3 → 3×2）");
    SHOW(Matrix<long long>({{1, 2}, {3, 4}}).trace(), "対角成分の和");
    {
        // EDPC R：隣接行列の K 乗の全成分の和 = 長さ K の道の数
        auto walks = [](const std::vector<std::vector<int>>& adj, long long K) {
            const int n = (int)adj.size();
            Matrix<mint> G(n, n);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) G(i, j) = adj[(size_t)i][(size_t)j];
            const Matrix<mint> P = G.pow(K);
            mint s = 0;
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) s += P(i, j);
            return s;
        };
        note("EDPC R（Walk）: 隣接行列の K 乗の全成分の和（グラフは graph_grid_v2 で読んで隣接行列にしてもよい）");
        std::cout << "    例 1 → " << walks({{0, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {1, 0, 0, 0}}, 2) << "、例 2 → " << walks({{0, 1, 0}, {1, 0, 1}, {0, 0, 0}}, 3)
                  << "、例 3 → " << walks({{0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0}}, 2)
                  << "、例 4 → " << walks({{0}}, 1) << "\n";
    }
    {
        // ABC293 E：[[A, 1], [0, 1]]^X の右上 = Σ_{i<X} A^i（mod M は合成数のこともある → DynModInt）
        auto solve = [](long long a, long long X, long long M) {
            DynModInt::set_mod(M);
            Matrix<DynModInt> T = {{DynModInt(a), DynModInt(1)}, {DynModInt(0), DynModInt(1)}};
            return T.pow(X)(0, 1);
        };
        note("ABC293 E（Geometric Sequence Sum）: [[A, 1], [0, 1]]^X の右上（Matrix<DynModInt>。00_math_base の geom_sum_mod でも）");
        std::cout << "    (3, 4, 7) → " << solve(3, 4, 7) << "、(8, 10, 9) → " << solve(8, 10, 9) << "、(10^9, 10^12, 998244353) → " << solve(1000000000, 1000000000000LL, 998244353) << "\n";
    }

    api("【2】 companion_matrix(c) / recurrence_nth_matrix(c, a, n)", "Matrix<T> / T",
        "a_i = c_0 a_{i−1} + … + c_{k−1} a_{i−k}（a_0〜a_{k−1} が a）の a_n。n は 0 から数える", "O(k^3 log n)");
    SHOW(recurrence_nth_matrix<mint>({1, 1}, {0, 1}, 90), "フィボナッチ F_90 mod 10^9+7");
    SHOW(recurrence_nth_matrix<mint>({1, 1, 1}, {0, 0, 1}, 1000000000000000000LL), "トリボナッチの 10^18 項目");
    SHOW(companion_matrix<long long>({2, 3}), "a_i = 2 a_{i−1} + 3 a_{i−2} の行列");
    note("k が大きい（数百以上）ときは 06_linear_recurrence（きたまさ法・Bostan–Mori）のほうが速い");

    api("【3】 半環：Matrix<MinPlus<long long>> / Matrix<OrAnd> / Matrix<XorAnd>", "Matrix<半環>",
        "+ と × を min と +（最短路）、OR と AND（行けるか）、XOR と AND（ビットの漸化式）にした行列の累乗", "O(n^3 log K)");
    RUN(Matrix<MinPlus<long long>> D(3, 3); D(0, 1) = 5; D(1, 2) = 1; D(2, 0) = 2; D(1, 0) = 10);
    SHOW(D.pow(4)(0, 1), "0 → 1 にちょうど 4 歩で行く最短：0→1→2→0→1 = 5+1+2+5");
    SHOW(D.pow(2)(0, 2), "ちょうど 2 歩：0→1→2 = 6");
    SHOW(D.pow(2)(0, 0), "ちょうど 2 歩で 0 → 0：0→1→0 = 15");
    SHOW(D.pow(1)(0, 2), "1 歩では行けない");
    RUN(Matrix<MinMax<int>> W(3, 3); W(0, 1) = 7; W(1, 2) = 3; W(0, 2) = 9; W(2, 2) = 1);
    SHOW(W.pow(2)(0, 2), "MinMax：ちょうど 2 歩で 0 → 2、通る辺の最大を最小に（0→1→2 なら max(7, 3) = 7）");
    RUN(Matrix<MaxMin<int>> Cap(3, 3); Cap(0, 1) = 7; Cap(1, 2) = 3; Cap(0, 2) = 9; Cap(2, 2) = 1);
    SHOW(Cap.pow(2)(0, 2), "MaxMin：通る辺の最小を最大に（0→1→2 なら 3、0→2→2 なら 1）");
    RUN(Matrix<OrAnd> R = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}});
    SHOW(R.pow(5)(0, 2), "ちょうど 5 歩で 0 → 2 に行けるか（3 歩で一周）");
    {
        // ABC009 D：A_{N+K} = XOR_j (C_j AND A_{N+K−j}) の M 番目（1 から数える）
        auto solve = [](const std::vector<unsigned long long>& as, const std::vector<unsigned long long>& cs, long long M) {
            std::vector<XorAnd> a(as.begin(), as.end()), c(cs.begin(), cs.end());
            return recurrence_nth_matrix(c, a, M - 1).v;
        };
        note("ABC009 D（漸化式）: XorAnd で recurrence_nth_matrix(C, A, M − 1)");
        std::cout << "    例 1 → " << solve({10, 20, 30}, {7, 19, 13}, 5) << "、例 2 → "
                  << solve({2345678901ULL, 1001001001ULL, 3333333333ULL, 3141592653ULL, 1234567890ULL}, {2147483648ULL, 2147483647ULL, 4294967295ULL, 4294967294ULL, 3434343434ULL}, 100) << "\n";
    }

    api("【4】 MatrixFixed<T, K>：M.pow(e) / M.det() / M.inv() / M.try_inv()", "MatrixFixed / T / Maybe",
        "K×K の小さい行列（大きさがコンパイル時に決まるので速い）。行列式・逆行列は T が mint か double", "O(K^3)");
    RUN(MatrixFixed<mint, 2> F = {{1, 1}, {1, 0}});
    SHOW(F.pow(10)(0, 1), "F_10");
    SHOW(F.det(), "1·0 − 1·1 = −1（mod 10^9+7）");
    SHOW(F.inv(), "逆行列 [[0, 1], [1, −1]]");
    RUN(MatrixFixed<double, 3> R3 = {{2, 0, 0}, {0, 3, 0}, {1, 0, 4}});
    SHOW(R3.det(), "2·3·4");
    SHOW(R3.try_inv().ok, "逆行列がある");
    SHOW((MatrixFixed<double, 2>{{1, 2}, {2, 4}}.try_inv().ok), "行列式 0 なので逆行列なし");
}

}  // namespace mat_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    mat_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mat_test::benchmark();    // 速度計測
    mat_demo::demo();         // 使い方デモ
    return 0;
}
