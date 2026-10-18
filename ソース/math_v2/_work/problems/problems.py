# 代表的な問題を「Main_base ＋ テンプレート本体だけ」で解いて、例を確かめる（作業用）
PROBLEMS = []
def P(name, files, solv, samples):
    PROBLEMS.append((name, files, solv, samples))

P('ABC186_E', ['03_mod_equation'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        ll N, S, K; cin >> N >> S >> K;
        ModSolution r = linear_congruence(K, -S, N);   // K x ≡ −S (mod N)
        cout << (r.ok ? r.r : -1) << '\n';
    }
}''', [("4\n10 4 3\n1000 11 2\n998244353 897581057 595591169\n10000 6 14\n", "2\n-1\n249561088\n3571\n")])

P('ABC193_E', ['03_mod_equation'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        ll X, Y, P, Q; cin >> X >> Y >> P >> Q;
        ll best = -1;
        for (ll a = X; a < X + Y; ++a)
            for (ll b = P; b < P + Q; ++b) {
                ModSolution s = crt_merge(a, 2 * X + 2 * Y, b, P + Q);
                if (s.ok && (best < 0 || s.r < best)) best = s.r;
            }
        if (best < 0) cout << "infinity\n"; else cout << best << '\n';
    }
}''', [("3\n5 2 7 6\n1 1 3 1\n999999999 1 1000000000 1\n", "20\ninfinity\n1000000000999999999\n")])

P('ABC156_D', ['04_binomial'], r'''
void solv(){
    using mint = ModInt107;
    ll n, a, b; cin >> n >> a >> b;
    cout << mint(2).pow(n) - 1 - binom_large_n<mint>(n, a) - binom_large_n<mint>(n, b) << '\n';
}''', [("4 1 3\n", "7\n"), ("1000000000 141421 173205\n", "34076506\n")])

P('ABC205_E', ['04_binomial'], r'''
void solv(){
    using mint = ModInt107;
    ll N, M, K; cin >> N >> M >> K;
    Binomial<mint> B(N + M);
    cout << B.paths_not_crossing(M, N, K) << '\n';   // 黒 M 個を x、白 N 個を y にして y ≤ x + K
}''', [("2 3 1\n", "9\n"), ("1 0 0\n", "0\n"), ("1000000 1000000 1000000\n", "192151600\n")])

P('ABC425_E', ['04_binomial'], r'''
void solv(){
    int T; ll M; cin >> T >> M;
    DynModInt::set_mod(M);
    auto C = pascal_table<DynModInt>(5000);
    while (T--) {
        int N; cin >> N;
        DynModInt r = 1; int s = 0;
        for (int i = 0; i < N; ++i) { int c; cin >> c; s += c; r *= C[s][c]; }
        cout << r << '\n';
    }
}''', [("3 1000000000\n2\n2 2\n5\n1 1 1 1 1\n6\n1 2 3 4 5 6\n", "6\n120\n230379200\n"),
       ("3 998244353\n1\n1\n3\n4 2 5\n10\n500 500 500 500 500 500 500 500 500 500\n", "1\n6930\n261233246\n")])

P('ARC117_C', ['04_binomial'], r'''
void solv(){
    int N; string s; cin >> N >> s;
    ll t = 0;
    for (int i = 0; i < N; ++i) t += binom_lucas(N - 1, i, 3) * (s[i] == 'B' ? 0 : s[i] == 'W' ? 1 : 2);
    if (N % 2 == 0) t = -t;
    cout << "BWR"[mod_floor(t, 3)] << '\n';
}''', [("3\nBWR\n", "W\n"), ("4\nRRBB\n", "W\n"), ("6\nBWWRBW\n", "B\n"), ("8\nWWBRBBWB\n", "R\n"), ("21\nBWBRRBBRWBRBBBRRBWWWR\n", "B\n")])

P('EDPC_R', ['05_matrix'], r'''
void solv(){
    using mint = ModInt107;
    int N; ll K; cin >> N >> K;
    Matrix<mint> A(N, N);
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) { int x; cin >> x; A(i, j) = x; }
    Matrix<mint> P = A.pow(K);
    mint s = 0;
    for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) s += P(i, j);
    cout << s << '\n';
}''', [("4 2\n0 1 0 0\n0 0 1 1\n0 0 0 1\n1 0 0 0\n", "6\n"), ("3 3\n0 1 0\n1 0 1\n0 0 0\n", "3\n"),
       ("6 2\n0 0 0 0 0 0\n0 0 1 0 0 0\n0 0 0 0 0 0\n0 0 0 0 1 0\n0 0 0 0 0 1\n0 0 0 0 0 0\n", "1\n"), ("1 1\n0\n", "0\n")])

P('ABC293_E', ['05_matrix'], r'''
void solv(){
    ll A, X, M; cin >> A >> X >> M;
    DynModInt::set_mod(M);
    Matrix<DynModInt> T = {{DynModInt(A), DynModInt(1)}, {DynModInt(0), DynModInt(1)}};
    cout << T.pow(X)(0, 1) << '\n';
}''', [("3 4 7\n", "5\n"), ("8 10 9\n", "0\n"), ("1000000000 1000000000000 998244353\n", "919667211\n")])

P('ABC009_D', ['05_matrix'], r'''
void solv(){
    int K; ll M; cin >> K >> M;
    vector<XorAnd> a(K), c(K);
    for (auto& x : a) { unsigned long long v; cin >> v; x = v; }
    for (auto& x : c) { unsigned long long v; cin >> v; x = v; }
    cout << recurrence_nth_matrix(c, a, M - 1) << '\n';
}''', [("3 5\n10 20 30\n7 19 13\n", "16\n"),
       ("5 100\n2345678901 1001001001 3333333333 3141592653 1234567890\n2147483648 2147483647 4294967295 4294967294 3434343434\n", "1067078691\n")])

P('ABC230_E', ['06_quotient'], r'''
void solv(){
    ll N; cin >> N;
    cout << sum_floor_div(N) << '\n';
}''', [("3\n", "5\n"), ("10000000000\n", "231802823220\n")])

P('PRACTICE2_C', ['06_quotient'], r'''
void solv(){
    int T; cin >> T;
    while (T--) { ll N, M, A, B; cin >> N >> M >> A >> B; cout << linear_floor_sum(N, M, A, B) << '\n'; }
}''', [("5\n4 10 6 3\n6 5 4 3\n1 1 0 0\n31415 92653 58979 32384\n1000000000 1000000000 999999999 999999999\n", "3\n13\n0\n314095480\n499999999500000000\n")])

P('ABC172_D', ['06_quotient'], r'''
void solv(){
    ll N; cin >> N;
    Int128 s = 0;
    for_each_quotient(N, [&](ll q, ll l, ll r) { s += sum_range(l, r) * (Int128)(q * (q + 1) / 2); });
    cout << s << '\n';
}''', [("4\n", "23\n"), ("100\n", "26879\n"), ("10000000\n", "838627288460105\n")])

P('ABC208_F', ['06_series'], r'''
void solv(){
    using mint = ModInt107;
    ll N; int M, K; cin >> N >> M >> K;
    const int D = K + M + 1;
    vector<mint> f(D);
    for (int n = 1; n < D; ++n) f[n] = mint(n).pow(K);
    for (int m = 1; m <= M; ++m) for (int n = 1; n < D; ++n) f[n] += f[n - 1];
    cout << lagrange_interpolate(f, N) << '\n';
}''', [("3 4 2\n", "35\n"), ("0 1 2\n", "0\n"), ("1000000000000000000 30 123456\n", "297085514\n")])

P('ABC245_D', ['06_series'], r'''
void solv(){
    int N, M; cin >> N >> M;
    vector<ll> A(N + 1), C(N + M + 1);
    cin >> A >> C;
    auto d = poly_divmod(C, A);
    for (int i = 0; i <= M; ++i) cout << d.q[i] << (i == M ? '\n' : ' ');
}''', [("1 2\n2 1\n12 14 8 2\n", "6 4 2\n"), ("1 1\n100 1\n10000 0 -1\n", "100 -1\n")])

P('ARC033_D', ['06_series'], r'''
void solv(){
    using mint = ModInt107;
    int N; cin >> N;
    vector<mint> y(N + 1);
    for (auto& v : y) { ll a; cin >> a; v = a; }
    ll T; cin >> T;
    cout << lagrange_interpolate(y, T) << '\n';
}''', [("2\n1 3 7\n3\n", "13\n"), ("5\n4 16 106 484 1624 4384\n1000000000\n", "999984471\n")])

P('ABC308_C', ['08_fraction'], r'''
void solv(){
    int N; cin >> N;
    vector<Fraction> r(N);
    for (auto& x : r) { ll a, b; cin >> a >> b; x = Fraction(a, a + b); }
    vector<int> id(N); iota(nall(id), 0);
    stable_sort(nall(id), [&](int i, int j) { return r[i] > r[j]; });
    for (int i = 0; i < N; ++i) cout << id[i] + 1 << (i + 1 == N ? '\n' : ' ');
}''', [("3\n1 3\n3 1\n2 2\n", "2 3 1\n"), ("2\n1 3\n2 6\n", "1 2\n"),
       ("4\n999999999 1000000000\n333333333 999999999\n1000000000 999999997\n999999998 1000000000\n", "3 1 4 2\n")])

P('ABC269_C', ['09_bit'], r'''
void solv(){
    ll N; cin >> N;
    for_each_submask(N, [&](ll s) { cout << s << '\n'; });
}''', [("11\n", "0\n1\n2\n3\n8\n9\n10\n11\n"), ("0\n", "0\n"),
       ("576461302059761664\n", "0\n524288\n549755813888\n549756338176\n576460752303423488\n576460752303947776\n576461302059237376\n576461302059761664\n")])

P('ABC356_D', ['09_bit'], r'''
void solv(){
    ll N, M; cin >> N >> M;
    ll s = 0;
    for_each_bit(M, [&](int i) { s = (s + count_bit_up_to(N, i)) % 998244353; });
    cout << s << '\n';
}''', [("4 3\n", "4\n"), ("0 0\n", "0\n"), ("1152921504606846975 1152921504606846975\n", "499791890\n")])

P('ABC121_D', ['09_bit'], r'''
void solv(){
    ll A, B; cin >> A >> B;
    cout << xor_range_closed(A, B) << '\n';
}''', [("2 4\n", "5\n"), ("123 456\n", "435\n"), ("123456789012 123456789012\n", "123456789012\n")])

P('ABC365_E', ['09_bit'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> a(N); cin >> a;
    Int128 s = subarray_xor_sum(a);
    for (ll x : a) s -= x;
    cout << s << '\n';
}''', [("3\n1 3 2\n", "3\n"), ("7\n2 5 6 5 2 1 7\n", "83\n")])

P('ABC363_D', ['01_digits'], r'''
void solv(){
    ll N; cin >> N;
    cout << kth_palindrome(N - 1) << '\n';
}''', [("46\n", "363\n"), ("1\n", "0\n"), ("1000000000000000000\n", "90000000000000000000000000000000009\n")])

P('ABC084_D', ['02_sieve'], r'''
void solv(){
    const int MX = 100000;
    Sieve S(MX);
    vector<int> pre(MX + 2, 0);
    for (int x = 1; x <= MX; ++x) pre[x + 1] = pre[x] + (x % 2 == 1 && S.is_prime(x) && S.is_prime((x + 1) / 2));
    int Q; cin >> Q;
    while (Q--) { int l, r; cin >> l >> r; cout << pre[r + 1] - pre[l] << '\n'; }
}''', [("1\n3 7\n", "2\n"), ("4\n13 13\n7 11\n7 11\n2017 2017\n", "1\n0\n0\n1\n"), ("6\n1 53\n13 91\n37 55\n19 51\n73 91\n13 49\n", "4\n4\n1\n1\n1\n2\n")])

P('ABC169_D', ['02_factorize'], r'''
void solv(){
    ll N; cin >> N;
    ll ans = 0;
    if (N > 1) for (const PrimePower& f : factorize(N)) { int e = f.e; for (int k = 1; e >= k; ++k) e -= k, ++ans; }
    cout << ans << '\n';
}''', [("24\n", "3\n"), ("1\n", "0\n"), ("64\n", "3\n"), ("1000000007\n", "1\n"), ("997764507000\n", "7\n")])

P('ABC110_D', ['04_binomial', '02_factorize'], r'''
void solv(){
    using mint = ModInt107;
    ll N, M; cin >> N >> M;
    Binomial<mint> B(N + 40);
    mint r = 1;
    if (M > 1) for (const PrimePower& f : factorize(M)) r *= B.H(N, f.e);
    cout << r << '\n';
}''', [("2 6\n", "4\n"), ("3 12\n", "18\n"), ("100000 1000000000\n", "957870001\n")])

# ---------------- Phase 4（☆） ----------------
P('PRACTICE2_F', ['07_convolution'], r'''
void solv(){
    using mint = ModInt998;
    int N, M; cin >> N >> M;
    vector<mint> a(N), b(M);
    for (auto& x : a) { ll v; cin >> v; x = v; }
    for (auto& x : b) { ll v; cin >> v; x = v; }
    vector<mint> c = convolve(a, b);
    for (size_t i = 0; i < c.size(); ++i) cout << c[i] << (i + 1 < c.size() ? ' ' : '\n');
}''', [("4 5\n1 2 3 4\n5 6 7 8 9\n", "5 16 34 60 70 70 59 36\n"), ("1 1\n10000000\n10000000\n", "871938225\n")])

P('ATC001_C', ['07_convolution'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> A(N + 1, 0), B(N + 1, 0);
    for (int i = 1; i <= N; ++i) { ll a, b; cin >> a >> b; A[i] = a, B[i] = b; }
    vector<ll> C = convolve_ll(A, B);
    for (int k = 1; k <= 2 * N; ++k) cout << C[k] << '\n';
}''', [("4\n1 1\n2 2\n3 4\n4 8\n", "0\n1\n4\n11\n26\n36\n40\n32\n")])

P('ABC196_F', ['07_convolution'], r'''
void solv(){
    using mint = ModInt998;
    string S, T; cin >> S >> T;
    vector<mint> s1(S.size()), s0(S.size()), t1(T.size()), t0(T.size());
    for (size_t i = 0; i < S.size(); ++i) s1[i] = S[i] == '1', s0[i] = S[i] == '0';
    for (size_t i = 0; i < T.size(); ++i) t1[i] = T[i] == '1', t0[i] = T[i] == '0';
    auto c1 = correlate(s1, t1), c0 = correlate(s0, t0);
    ll best = (ll)T.size();
    for (size_t k = 0; k < c1.size(); ++k) chmin(best, (ll)T.size() - (ll)(c1[k] + c0[k]).val());
    cout << best << '\n';
}''', [("0001\n101\n", "1\n"), ("0101010\n1010101\n", "7\n"), ("10101000010011011110\n0010011111\n", "1\n")])

P('ABC230_Ex', ['07_fps'], r'''
void solv(){
    // F = x (MSET(F + G) − 1)。F(x^k)（k ≥ 2）は前半の係数だけで決まるので、ニュートン法で係数を倍々に
    using mint = ModInt998;
    int W, K; cin >> W >> K;
    vector<mint> G(W + 1, 0);
    rep(i, K) { int w; cin >> w; G[w] = 1; }
    const int n = W + 1;
    Fps<mint> F = {0};
    for (int m = 1; m < n;) {
        const int len = min(2 * m, n);
        Fps<mint> S(len, 0);  // Σ_{k≥2} (F(x^k) + G(x^k)) / k ＋ G
        for (int k = 2; k < len; ++k) {
            const mint ik = mint(k).inv();
            for (int i = 1; (ll)i * k < len; ++i) S[i * k] += ((i < (int)F.size() ? F[i] : mint(0)) + G[i]) * ik;
        }
        for (int i = 0; i < len; ++i) S[i] += G[i];
        const Fps<mint> B = fps_exp(S, len);  // F = x (B e^F − 1) を解く
        Fps<mint> eF = fps_exp(F.pre(len), len);
        Fps<mint> BeF = (B * eF).pre(len);
        Fps<mint> phi(len, 0), dphi(len, 0);  // Φ = F − x (B e^F − 1)、Φ' = 1 − x B e^F
        for (int i = 0; i < len; ++i) {
            const mint xb = i ? BeF[i - 1] : mint(0);
            phi[i] = (i < (int)F.size() ? F[i] : mint(0)) - xb + (i == 1 ? mint(1) : mint(0));
            dphi[i] = (i == 0 ? mint(1) : mint(0)) - xb;
        }
        const Fps<mint> corr = (phi * fps_inv(dphi, len)).pre(len);
        F.resize(len);
        for (int i = 0; i < len; ++i) F[i] -= corr[i];
        m = len;
    }
    for (int w = 2; w <= W; ++w) cout << F[w] << '\n';
}''', [("4 1\n1\n", "1\n2\n4\n"), ("10 10\n1 2 3 4 5 6 7 8 9 10\n", "1\n3\n7\n18\n45\n121\n325\n904\n2546\n")])

P('ABC212_H', ['07_subset_transform'], r'''
void solv(){
    using mint = ModInt998;
    ll N; int K; cin >> N >> K;
    vector<mint> h(1 << 16, 0);
    rep(i, K) { int a; cin >> a; h[a] = 1; }
    walsh_hadamard(h);
    for (mint& x : h) x = (x == mint(1)) ? mint(N) : x * (x.pow(N) - 1) / (x - 1);
    walsh_hadamard(h, true);
    const mint k = K;
    const mint all = (k == mint(1)) ? mint(N) : k * (k.pow(N) - 1) / (k - 1);
    cout << all - h[0] << '\n';
}''', [("2 2\n1 2\n", "4\n"), ("100 3\n3 5 7\n", "112184936\n")])

P('ABC206_E', ['02_divisor_transform'], r'''
void solv(){
    ll L, R; cin >> L >> R;
    // gcd(x, y) = g ≥ 2 の組を数え、g = x か g = y のものを引く
    const int N = (int)R;
    vector<ll> f(N + 1, 0);  // f[g] = gcd がちょうど g の (x, y) の数（x, y ∈ [L, R]）
    vector<ll> cnt(N + 1, 0);
    for (ll g = 1; g <= N; ++g) { ll c = R / g - (L - 1) / g; f[g] = c * c; }
    multiple_mobius(f);
    ll ans = 0;
    for (ll g = max(2LL, 1LL); g <= N; ++g) {
        ans += f[g];
        if (g >= L) ans -= 2 * (R / g) - 1;  // x = g か y = g（g 自身の倍数）
    }
    cout << ans << '\n';
}''', [("3 7\n", "2\n"), ("4 10\n", "12\n"), ("1 1000000\n", "392047955148\n")])

P('ABC162_E', ['02_divisor_transform'], r'''
void solv(){
    ll N, K; cin >> N >> K;
    cout << sum_gcd_tuples(K, N, 1000000007) << '\n';
}''', [("3 2\n", "9\n"), ("3 200\n", "10813692\n"), ("100000 100000\n", "742202979\n")])

P('ABC300_G', ['02_prime_count'], r'''
void solv(){
    ll N, P; cin >> N >> P;
    cout << count_smooth(N, P) << '\n';
}''', [("36 3\n", "14\n"), ("10000000000000000 97\n", "2345134674\n")])

P('ABC222_G', ['03_discrete_log'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        ll K; cin >> K;
        const ll M = 9 * (K % 2 == 0 ? K / 2 : K);
        const Maybe<ll> o = mult_order(10, M);
        cout << (o.ok ? o.value : -1) << '\n';
    }
}''', [("4\n1\n7\n10\n999983\n", "1\n6\n-1\n999982\n")])

P('ABC335_G', ['03_discrete_log'], r'''
void solv(){
    int N; ll P; cin >> N >> P;
    const MultOrder ord(P);
    map<ll, ll> cnt;
    rep(i, N) { ll a; cin >> a; cnt[ord(a).value]++; }
    ll ans = 0;
    for (auto& [x, cx] : cnt) for (auto& [y, cy] : cnt) if (x % y == 0) ans += cx * cy;
    cout << ans << '\n';
}''', [("3 13\n2 3 5\n", "5\n"), ("5 2\n1 1 1 1 1\n", "25\n"),
       ("10 9999999999971\n141592653589 793238462643 383279502884 197169399375 105820974944 592307816406 286208998628 034825342117 067982148086 513282306647\n", "63\n")])

P('ABC228_E', ['03_discrete_log'], r'''
void solv(){
    ll N, K, M; cin >> N >> K >> M;
    cout << pow_tower_mod({M, K, N}, 998244353) << '\n';
}''', [("2 2 2\n", "16\n"), ("3 14 15926535\n", "109718301\n")])

P('ABC131_C', ['04_combinatorics'], r'''
void solv(){
    ll A, B, C, D; cin >> A >> B >> C >> D;
    cout << (B - A + 1) - count_divisible_by_any_closed(A, B, {C, D}) << '\n';
}''', [("4 9 2 3\n", "2\n"), ("10 40 6 8\n", "23\n"), ("314159265358979323 846264338327950288 419716939 937510582\n", "532105071133627368\n")])

P('ABC242_F', ['04_combinatorics'], r'''
void solv(){
    using mint = ModInt998;
    int N, M, Bk, Wk; cin >> N >> M >> Bk >> Wk;
    Binomial<mint> B(N * M + N + M);
    auto f = [&](int k) {
        vector<vector<mint>> t(N + 1, vector<mint>(M + 1));
        erep(r, N) erep(c, M) t[r][c] = inclusion_exclusion_count<mint>((int)r, [&](int rr) { return inclusion_exclusion_count<mint>((int)c, [&](int cc) { return B.C(rr * cc, k); }); });
        return t;
    };
    auto fb = f(Bk), fw = f(Wk);
    mint s = 0;
    erep(r1, N) erep(c1, M) {
        if (fb[r1][c1] == mint(0)) continue;
        for (int r2 = 0; r1 + r2 <= N; ++r2) for (int c2 = 0; c1 + c2 <= M; ++c2)
            s += B.C(N, r1) * B.C(N - r1, r2) * B.C(M, c1) * B.C(M - c1, c2) * fb[r1][c1] * fw[r2][c2];
    }
    cout << s << '\n';
}''', [("2 2 1 1\n", "4\n"), ("1 2 1 1\n", "0\n"), ("40 40 30 30\n", "467620384\n")])

P('ABC297_F', ['04_combinatorics'], r'''
void solv(){
    using mint = ModInt998;
    ll H, W, K; cin >> H >> W >> K;
    Binomial<mint> B(H * W);
    mint s = 0;
    for (ll h = 1; h <= H; ++h) for (ll w = 1; w <= W; ++w) {
        const mint ways = inclusion_exclusion<mint>(4, [&](ll S) {
            const ll hh = max(0LL, h - (S & 1) - (S >> 1 & 1)), ww = max(0LL, w - (S >> 2 & 1) - (S >> 3 & 1));
            return B.C(hh * ww, K);
        });
        s += ways * mint(h * w) * mint((H - h + 1) * (W - w + 1));
    }
    cout << s / B.C(H * W, K) << '\n';
}''', [("2 2 2\n", "665496238\n"), ("10 10 1\n", "1\n"), ("314 159 2653\n", "639716353\n")])

P('ABC190_F', ['04_permutation'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> a(N); cin >> a;
    ll c = inversion_count(a);
    rep(k, N) { cout << c << '\n'; c += (N - 1 - a[k]) - a[k]; }
}''', [("4\n0 1 2 3\n", "0\n3\n4\n3\n"), ("10\n0 3 1 5 4 2 9 6 8 7\n", "9\n18\n21\n28\n27\n28\n33\n24\n21\n14\n")])

P('ABC276_C', ['04_permutation'], r'''
void solv(){
    int N; cin >> N;
    vector<int> P(N); cin >> P;
    vector<int> q = perm_unrank(N, perm_rank(P) - 1);
    rep(i, N) cout << q[i] + 1 << (i + 1 < N ? ' ' : '\n');
}''', [("3\n3 1 2\n", "2 3 1\n"), ("10\n9 8 6 5 10 3 1 2 4 7\n", "9 8 6 5 10 2 7 4 3 1\n")])

P('ABC202_D', ['04_permutation'], r'''
void solv(){
    int A, B; ll K; cin >> A >> B >> K;
    cout << multiset_perm_unrank(string(A, 'a') + string(B, 'b'), K - 1) << '\n';
}''', [("2 2 4\n", "baab\n"), ("30 30 118264581564861424\n", "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n")])

P('ABC253_Ex', ['05_linear_algebra'], r'''
void solv(){
    using mint = ModInt998;
    int N, M; cin >> N >> M;
    vector<pair<int, int>> E(M);
    for (auto& [u, v] : E) { cin >> u >> v; --u, --v; }
    const int FULL = (1 << N) - 1;
    vector<mint> tree(FULL + 1, 0);
    for (int S = 1; S <= FULL; ++S) {
        vector<int> id(N, -1), us, vs; int k = 0;
        rep(v, N) if (S >> v & 1) id[v] = k++;
        for (auto& [u, v] : E) if (id[u] >= 0 && id[v] >= 0) us.push_back(id[u]), vs.push_back(id[v]);
        tree[S] = count_spanning_trees<mint>(k, us, vs);
    }
    vector<vector<mint>> g(N + 1, vector<mint>(FULL + 1, 0));
    g[0][0] = 1;
    for (int S = 1; S <= FULL; ++S) {
        const int low = S & -S;
        for (int T = S; T > 0; T = (T - 1) & S) if (T & low) erep(c, 1, N) g[c][S] += tree[T] * g[c - 1][S ^ T];
    }
    mint fact = 1;
    erep(K, 1, N - 1) { fact *= mint(K); cout << g[N - K][FULL] * fact / mint(M).pow(K) << '\n'; }
}''', [("3 2\n1 2\n2 3\n", "1\n499122177\n"), ("4 5\n1 2\n1 2\n1 4\n2 3\n2 4\n", "1\n758665709\n918384805\n")])

P('ABC323_G', ['05_linear_algebra'], r'''
void solv(){
    using mint = ModInt998;
    int N; cin >> N;
    vector<int> P(N); cin >> P;
    Matrix<mint> L0(N - 1, N - 1), L1(N - 1, N - 1);
    rep(i, N) rep(j, i + 1, N) {
        Matrix<mint>& L = P[i] > P[j] ? L1 : L0;
        for (int v : {(int)i, (int)j}) if (v < N - 1) L(v, v) += mint(1);
        if (j < N - 1) L(i, j) -= mint(1), L(j, i) -= mint(1);
    }
    vector<mint> d = det_affine(L0, L1);
    rep(k, N) cout << d[k] << (k + 1 < N ? ' ' : '\n');
}''', [("3\n1 3 2\n", "1 2 0\n"), ("10\n3 1 4 10 8 6 9 2 7 5\n", "294448 2989776 12112684 25422152 30002820 20184912 7484084 1397576 108908 2640\n")])

P('ABC141_F', ['05_xor_basis'], r'''
void solv(){
    int N; cin >> N;
    vector<ull> A(N); cin >> A;
    ull T = 0; for (ull a : A) T ^= a;
    XorBasis B; for (ull a : A) B.insert(a & ~T);
    cout << T + 2 * B.max_xor() << '\n';
}''', [("3\n3 6 5\n", "12\n"), ("4\n23 36 66 65\n", "188\n"),
       ("20\n1008288677408720767 539403903321871999 1044301017184589821 215886900497862655 504277496111605629 972104334925272829 792625803473366909 972333547668684797 467386965442856573 755861732751878143 1151846447448561405 467257771752201853 683930041385277311 432010719984459389 319104378117934975 611451291444233983 647509226592964607 251832107792119421 827811265410084479 864032478037725181\n", "2012721721873704572\n")])

P('ABC223_H', ['05_xor_basis'], r'''
void solv(){
    int N, Q; cin >> N >> Q;
    vector<ull> A(N); cin >> A;
    vector<array<ll, 3>> qs(Q);
    for (auto& q : qs) cin >> q[0] >> q[1] >> q[2];
    vector<int> ord(Q); iota(nall(ord), 0);
    sort(nall(ord), [&](int a, int b) { return qs[a][1] < qs[b][1]; });
    PrefixXorBasis P; vector<bool> ans(Q);
    for (int q : ord) {
        while (P.size() < qs[q][1]) P.push_back(A[P.size()]);
        ans[q] = P.contains((int)qs[q][0] - 1, (ull)qs[q][2]);
    }
    rep(i, Q) { if (ans[i]) YES; else NO; }
}''', [("5 2\n3 1 4 1 5\n1 3 7\n2 5 7\n", "Yes\nNo\n"),
       ("10 10\n8 45 56 9 38 28 33 5 15 19\n10 10 53\n3 8 60\n1 10 29\n5 7 62\n3 7 51\n8 8 52\n1 4 60\n6 8 32\n4 8 58\n5 9 2\n", "No\nNo\nYes\nNo\nYes\nNo\nNo\nNo\nYes\nYes\n")])

P('ABC283_G', ['05_xor_basis'], r'''
void solv(){
    int N; ll L, R; cin >> N >> L >> R;
    XorBasis B; rep(i, N) { ull a; cin >> a; B.insert(a); }
    for (ll k = L - 1; k < R; ++k) cout << B.kth((ull)k) << (k + 1 < R ? ' ' : '\n');
}''', [("4 1 8\n2 21 17 21\n", "0 2 4 6 17 19 21 23\n"), ("4 3 7\n2 21 17 21\n", "4 6 17 19 21\n"), ("5 1 1\n0 0 0 0 0\n", "0\n"), ("6 21 21\n88 44 82 110 121 80\n", "41\n")])

P('ABC300_Ex', ['06_linear_recurrence'], r'''
void solv(){
    using mint = ModInt998;
    int K; ll N; cin >> K >> N;
    vector<mint> Q(K + 1), A(K, mint(1));
    Q[0] = 1; erep(j, 1, K) Q[j] = -1;
    vector<mint> P = convolve(A, Q); P.resize(K);
    for (; N > 0; N >>= 1) {
        vector<mint> Qm(Q);
        for (size_t i = 1; i < Qm.size(); i += 2) Qm[i] = -Qm[i];
        const vector<mint> U = convolve(P, Qm), V = convolve(Q, Qm);
        vector<mint> nP(K, 0), nQ(K + 1, 0);
        rep(i, sizei(U)) if (i % 2 == 0 || (N & 1)) nP[i / 2] += U[i];
        for (int i = 0; i < sizei(V); i += 2) nQ[i / 2] = V[i];
        P.swap(nP), Q.swap(nQ);
    }
    cout << P[0] / Q[0] << '\n';
}''', [("2 6\n", "21\n"), ("2 8\n", "35\n"), ("1 123456789\n", "65536\n"), ("300 20230429\n", "125461938\n"), ("42923 999999999558876113\n", "300300300\n")])

P('ABC169_B', ['08_bigint'], r'''
void solv(){
    int N; cin >> N;
    vector<BigInt> A(N); cin >> A;
    BigInt p = 1;
    for (auto& a : A) if (a.is_zero()) { cout << 0 << '\n'; return; }
    for (auto& a : A) { p *= a; if (p > BigInt(1000000000000000000LL)) { cout << -1 << '\n'; return; } }
    cout << p << '\n';
}''', [("2\n1000000000 1000000000\n", "1000000000000000000\n"), ("3\n101 9901 999999000001\n", "-1\n"),
       ("31\n4 1 5 9 2 6 5 3 5 8 9 7 9 3 2 3 8 4 6 2 6 4 3 3 8 3 2 7 9 5 0\n", "0\n")])

P('ABC206_F', ['10_game'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        int N; cin >> N;
        vector<pair<int, int>> iv(N); for (auto& [l, r] : iv) cin >> l >> r;
        vector<vector<int>> g(101, vector<int>(101, 0));
        for (int len = 1; len <= 99; ++len) for (int l = 1; l + len <= 100; ++l) {
            const int r = l + len; vector<int> ws;
            for (auto& [a, b] : iv) if (l <= a && b <= r) ws.push_back(g[l][a] ^ g[b][r]);
            g[l][r] = mex(ws);
        }
        cout << (g[1][100] ? "Alice" : "Bob") << '\n';
    }
}''', [("5\n3\n53 98\n8 43\n12 53\n10\n4 7\n5 7\n3 7\n4 5\n5 8\n6 9\n4 8\n5 10\n1 9\n5 10\n2\n58 98\n11 29\n6\n79 83\n44 83\n38 74\n49 88\n18 45\n64 99\n1\n5 9\n",
        "Bob\nAlice\nBob\nAlice\nAlice\n")])

P('ABC297_G', ['10_game'], r'''
void solv(){
    int N; ll L, R; cin >> N >> L >> R;
    ll x = 0; rep(i, N) { ll a; cin >> a; x ^= range_take_grundy(a, L, R); }
    cout << (x ? "First" : "Second") << '\n';
}''', [("3 1 2\n2 3 3\n", "First\n"), ("5 1 1\n3 1 4 1 5\n", "Second\n"), ("7 3 14\n10 20 30 40 50 60 70\n", "First\n")])

P('ABC368_F', ['10_game'], r'''
void solv(){
    int N; cin >> N;
    vector<int> A(N); cin >> A;
    const int M = 100000;
    vector<vector<int>> moves(M + 1);
    for (int d = 1; d <= M; ++d) for (int m = 2 * d; m <= M; m += d) moves[m].push_back(d);
    const vector<int> g = grundy_dag(M + 1, moves);
    int x = 0; for (int a : A) x ^= g[a];
    cout << (x ? "Anna" : "Bruno") << '\n';
}''', [("3\n2 3 4\n", "Anna\n"), ("4\n2 3 4 6\n", "Bruno\n")])

# ================= 最終照合（P1）で解いた問題 =================
# 数学テンプレート（第二世代）の最終照合 P1：Main_base ＋ テンプレート本体だけで解いて、問題ページのサンプルを確かめる

# ---------------- 00_int_common / 00_math_base ----------------
P('TYP020', ['00_math_base'], r'''
void solv(){
    ll a, b, c; cin >> a >> b >> c;
    // log2 a < b log2 c ⇔ a < c^b（a ≤ 9×10^18 なので a + 1 で止める）
    if (pow_cap(c, b, a + 1) > a) YES; else NO;
}''', [("4 3 2\n", "Yes\n"), ("16 3 2\n", "No\n"), ("8 3 2\n", "No\n")])

P('ABC169_C', ['00_int_common'], r'''
void solv(){
    ll A; cin >> A;
    const ll B100 = read_fixed(2);   // B × 100（誤差なし）
    cout << A * B100 / 100 << '\n';
}''', [("198 1.10\n", "217\n"), ("1 0.01\n", "0\n"), ("1000000000000000 9.99\n", "9990000000000000\n")])

P('TYP038', ['00_math_base'], r'''
void solv(){
    ll A, B; cin >> A >> B;
    const ll LIM = 1000000000000000000LL;
    const ll l = lcm_cap(A, B, LIM + 1);
    if (l > LIM) cout << "Large\n"; else cout << l << '\n';
}''', [("4 6\n", "12\n"), ("1000000000000000000 3\n", "Large\n"), ("1000000000000000000 1\n", "1000000000000000000\n")])

P('TYP022', ['00_math_base'], r'''
void solv(){
    ll A, B, C; cin >> A >> B >> C;
    const ll g = gcd_of({A, B, C});
    cout << (A / g - 1) + (B / g - 1) + (C / g - 1) << '\n';
}''', [("2 2 3\n", "4\n"), ("2 2 4\n", "1\n"),
       ("1000000000000000000 999999999999999999 999999999999999998\n", "2999999999999999994\n")])

P('ABC227_C', ['00_math_base'], r'''
void solv(){
    ll N; cin >> N;
    ll cnt = 0;
    const ll A3 = kth_root_floor(N, 3);
    for (ll a = 1; a <= A3; ++a) {
        const ll B2 = isqrt(N / a);   // a b^2 ≤ N
        for (ll b = a; b <= B2; ++b) cnt += N / (a * b) - b + 1;
    }
    cout << cnt << '\n';
}''', [("4\n", "5\n"), ("100\n", "323\n"), ("100000000000\n", "5745290566750\n")])

P('ARC059_C', ['00_math_base'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> a(N); cin >> a;
    const ll y = round_div(accumulate(nall(a), 0LL), N);   // 平均にいちばん近い整数
    ll cost = 0;
    for (ll x : a) cost += (x - y) * (x - y);
    cout << cost << '\n';
}''', [("2\n4 8\n", "8\n"), ("3\n1 1 3\n", "3\n"), ("3\n4 2 5\n", "5\n"), ("4\n-100 -100 -100 -100\n", "0\n")])

# ---------------- 00_modint ----------------
P('TYP069', ['00_modint'], r'''
void solv(){
    using mint = ModInt107;
    ll N, K; cin >> N >> K;
    // N = 1 を分けないと (K − 2).pow(−1) が逆元になる
    const mint ans = N == 1 ? mint(K) : mint(K) * (K - 1) * mint(K - 2).pow(N - 2);
    cout << ans << '\n';
}''', [("2 3\n", "6\n"), ("10 2\n", "0\n"), ("2021 617\n", "53731843\n")])

P('ABC280_E', ['00_modint'], r'''
void solv(){
    using mint = ModInt998;
    int N, P; cin >> N >> P;
    const mint p = mint(P) / 100;
    vector<mint> dp(N + 1, 0);
    erep(i, 1, N) dp[i] = 1 + p * (i >= 2 ? dp[i - 2] : mint(0)) + (1 - p) * dp[i - 1];
    cout << dp[N] << '\n';
}''', [("3 10\n", "229596204\n"), ("5 100\n", "3\n"), ("280 59\n", "567484387\n")])

P('ABC263_E', ['00_modint'], r'''
void solv(){
    using mint = ModInt998;
    int N; cin >> N;
    vector<ll> A(N, 0);
    rep(i, 1, N) cin >> A[i];   // A[1..N−1]
    // dp[i] = マス i から N までの回数の期待値。dp[i] = 1 + (dp[i] + dp[i+1] + … + dp[i+A]) / (A + 1)
    //  → dp[i] = (A + 1 + dp[i+1] + … + dp[i+A]) / A
    vector<mint> dp(N + 2, 0), suf(N + 2, 0);   // suf[i] = dp[i] + dp[i+1] + … + dp[N]
    for (int i = N - 1; i >= 1; --i) {
        dp[i] = (mint(A[i] + 1) + suf[i + 1] - suf[i + A[i] + 1]) / A[i];
        suf[i] = suf[i + 1] + dp[i];
    }
    cout << dp[1] << '\n';
}''', [("3\n1 1\n", "4\n"), ("5\n3 1 2 1\n", "332748122\n")])

# ---------------- 01_digits ----------------
P('ABC083_B', ['01_digits'], r'''
void solv(){
    ll N, A, B; cin >> N >> A >> B;
    ll s = 0;
    erep(x, 1, N) {
        const ll d = digit_sum(x);
        if (A <= d && d <= B) s += x;
    }
    cout << s << '\n';
}''', [("20 2 5\n", "84\n"), ("10 1 2\n", "13\n"), ("100 4 16\n", "4554\n")])

P('ABC105_C', ['01_digits'], r'''
void solv(){
    ll N; cin >> N;
    cout << to_base(N, -2) << '\n';
}''', [("-9\n", "1011\n"), ("123456789\n", "11000101011001101110100010101\n"), ("0\n", "0\n")])

P('ABC186_C', ['01_digits'], r'''
void solv(){
    ll N; cin >> N;
    ll c = 0;
    erep(x, 1, N) {
        const vector<int> d10 = digits_of(x), d8 = digits_of(x, 8);
        if (count(nall(d10), 7) == 0 && count(nall(d8), 7) == 0) ++c;
    }
    cout << c << '\n';
}''', [("20\n", "17\n"), ("100000\n", "30555\n")])

P('ABC336_C', ['01_digits'], r'''
void solv(){
    ll N; cin >> N;
    vector<int> d = digits_of(N - 1, 5);   // N − 1 の 5 進数の各桁を 2 倍
    for (int& x : d) x *= 2;
    cout << from_digits(d) << '\n';
}''', [("8\n", "24\n"), ("133\n", "2024\n"), ("31415926535\n", "2006628868244228\n")])

P('TYP067', ['01_digits'], r'''
void solv(){
    string s; int K; cin >> s >> K;
    rep(k, K) {
        s = to_base(from_base(s, 8), 9);
        for (char& c : s) if (c == '8') c = '5';
    }
    cout << s << '\n';
}''', [("21 1\n", "15\n"), ("1330 1\n", "555\n"), ("2311640221315 15\n", "474547\n")])

P('ABC192_D', ['01_digits'], r'''
void solv(){
    string X; ll M; cin >> X >> M;
    if (sizei(X) == 1) { cout << ((X[0] - '0') <= M ? 1 : 0) << '\n'; return; }
    ll dmax = 0;
    for (char c : X) chmax(dmax, (ll)(c - '0'));
    ll lo = dmax, hi = M + 1;   // lo：満たす最大の進数（dmax は「満たさない」の番兵）
    while (hi - lo > 1) {
        const ll mid = lo + (hi - lo) / 2;
        (base_value_leq(X, mid, M) ? lo : hi) = mid;
    }
    cout << lo - dmax << '\n';
}''', [("22\n10\n", "2\n"), ("999\n1500\n", "3\n"),
       ("100000000000000000000000000000000000000000000000000000000000\n1000000000000000000\n", "1\n")])

P('TYP082', ['01_digits'], r'''
void solv(){
    ll L, R; cin >> L >> R;
    const ll MOD = 1000000007;
    ll ans = 0;
    for (const DigitLengthBlock& b : digit_length_blocks(L, R + 1))
        ans = (ans + mod_mul(b.len, arith_sum_mod(b.lo, 1, b.hi - b.lo, MOD), MOD)) % MOD;
    cout << ans << '\n';
}''', [("3 5\n", "12\n"), ("98 100\n", "694\n"), ("1001 869120\n", "59367733\n"),
       ("381453331666495446 746254773042091083\n", "584127830\n")])

P('ABC357_D', ['01_digits'], r'''
void solv(){
    ll N; cin >> N;
    cout << repeat_concat_mod(N, N, 998244353) << '\n';
}''', [("5\n", "55555\n"), ("9\n", "1755646\n"), ("10000000000\n", "468086693\n")])

P('ARC141_A', ['01_digits', '02_factorize'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        string s; cin >> s;
        const int L = sizei(s);
        const ll N = substr_value(s, 0, L);
        ll best = L >= 3 ? pow10_ll(L - 1) - 1 : 0;   // L − 1 桁の 99…9（2 桁以上なら周期的）
        for (ll d : divisors_of(L)) {
            if (d == L) continue;
            const ll rep1 = (pow10_ll(L) - 1) / (pow10_ll((int)d) - 1);   // 100…0100…01（d 桁ごとに 1）
            ll P = substr_value(s, 0, (int)d);   // 先頭 d 桁をくり返す
            if (P * rep1 > N) --P;
            if (P >= pow10_ll((int)d - 1)) chmax(best, P * rep1);
        }
        cout << best << '\n';
    }
}''', [("3\n1412\n23\n498650499498649123\n", "1313\n22\n498650498650498650\n")])

# ---------------- 02_sieve ----------------
P('ABC300_D', ['02_sieve'], r'''
void solv(){
    ll N; cin >> N;
    const vector<ll> ps = primes_up_to(isqrt(N / 12));   // c ≤ √(N / (2^2·3))
    const int n = sizei(ps);
    ll cnt = 0;
    for (int i = 0; i < n; ++i) {
        const ll a = ps[i];
        if (pow_cap(a, 5, N + 1) > N) break;   // a^2 b c^2 > a^5
        for (int j = i + 1; j < n; ++j) {
            const ll b = ps[j];
            if (mul_cap(a * a, pow_cap(b, 3, N + 1), N + 1) > N) break;   // c > b
            const ll cmax = isqrt(N / (a * a * b));
            const ll k = upper_bound(nall(ps), cmax) - ps.begin();   // c は ps[j+1..k)
            if (k <= j + 1) break;
            cnt += k - (j + 1);
        }
    }
    cout << cnt << '\n';
}''', [("1000\n", "3\n"), ("1000000000000\n", "2817785\n")])

P('ABC383_D', ['02_sieve'], r'''
void solv(){
    ll N; cin >> N;
    const ll R = isqrt(N);
    const auto ps = primes_up_to(R);
    ll c = upper_bound(nall(ps), kth_root_floor(N, 8)) - ps.begin();   // p^8 ≤ N
    for (size_t i = 0; i < ps.size(); ++i)   // p^2 q^2 ≤ N ⇔ p q ≤ ⌊√N⌋
        for (size_t j = i + 1; j < ps.size() && (ll)ps[i] * ps[j] <= R; ++j) ++c;
    cout << c << '\n';
}''', [("200\n", "3\n"), ("4000000000000\n", "407073\n")])

P('ABC400_E', ['02_sieve'], r'''
void solv(){
    const int M = 1000000;
    const vector<int> om = omega_table(M);
    vector<ll> good;   // 異なる素因数がちょうど 2 個の a（400 number = a^2）
    erep(a, 1, M) if (om[a] == 2) good.push_back(a);
    int Q; cin >> Q;
    while (Q--) {
        ll A; cin >> A;
        const ll a = *prev(upper_bound(nall(good), isqrt(A)));
        cout << a * a << '\n';
    }
}''', [("5\n404\n36\n60\n1000000000000\n123456789\n", "400\n36\n36\n1000000000000\n123454321\n")])

P('ABC177_E', ['02_sieve'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> A(N); cin >> A;
    const Sieve S(1000000);
    vector<int> seen(1000001, 0);
    bool pairwise = true;
    for (ll a : A) for (ll p : S.prime_factors(a)) if (++seen[p] >= 2) pairwise = false;
    if (pairwise) cout << "pairwise coprime\n";
    else if (gcd_of(A) == 1) cout << "setwise coprime\n";
    else cout << "not coprime\n";
}''', [("3\n3 4 5\n", "pairwise coprime\n"), ("3\n6 10 15\n", "setwise coprime\n"), ("3\n6 10 16\n", "not coprime\n")])

# ================= 最終照合（P2）で解いた問題 =================
# 数学テンプレート（第二世代）の最終照合 P2：Main_base ＋ テンプレート本体だけで解いて、問題ページの例を確かめる

# ---------------- 素数・約数 ----------------
P('TYP030', ['02_sieve'], r'''
void solv(){
    int N, K; cin >> N >> K;
    const vector<int> om = omega_table(N);   // 異なる素因数の数の表
    int c = 0;
    erep(x, 2, N) c += om[x] >= K;
    cout << c << '\n';
}''', [("15 2\n", "5\n"), ("30 2\n", "13\n"), ("200 4\n", "0\n"), ("869120 1\n", "869119\n"), ("10000000 3\n", "6798027\n")])

P('ABC412_E', ['02_sieve'], r'''
void solv(){
    ll L, R; cin >> L >> R;
    // lcm(1..n) は n が素数の累乗のときだけ変わる → 1 + ((L, R] の素数の累乗の数)
    ll c = 1;
    segment_factorize_closed(L + 1, R, [&](ll, const vector<PrimePower>& fs) { c += fs.size() == 1; });
    cout << c << '\n';
}''', [("4 12\n", "6\n"), ("123456789 123456789\n", "1\n"), ("99999990000000 100000000000000\n", "310458\n")])

P('ABC180_C', ['02_factorize'], r'''
void solv(){
    ll N; cin >> N;
    for (ll d : divisors_of(N)) cout << d << '\n';
}''', [("6\n", "1\n2\n3\n6\n"),
       ("720\n", "1\n2\n3\n4\n5\n6\n8\n9\n10\n12\n15\n16\n18\n20\n24\n30\n36\n40\n45\n48\n60\n72\n80\n90\n120\n144\n180\n240\n360\n720\n"),
       ("1000000007\n", "1\n1000000007\n")])

P('ABC190_D', ['02_factorize'], r'''
void solv(){
    ll N; cin >> N;
    ll odd = N;
    while (odd % 2 == 0) odd /= 2;
    cout << 2 * count_divisors(odd) << '\n';   // 2 × (N の奇数の約数の個数)
}''', [("12\n", "4\n"), ("1\n", "2\n"), ("963761198400\n", "1920\n")])

P('ABC215_D', ['02_sieve'], r'''
void solv(){
    int N, M; cin >> N >> M;
    const int MX = 100000;
    Sieve S(MX);
    vector<char> badp(MX + 1, 0), bad(M + 1, 0);
    rep(i, N) { ll a; cin >> a; for (ll p : S.prime_factors(a)) badp[p] = 1; }
    erep(p, 2, MX) if (badp[p]) for (ll m = p; m <= M; m += p) bad[m] = 1;
    vector<int> ans;
    erep(k, 1, M) if (!bad[k]) ans.push_back((int)k);
    cout << ans.size() << '\n';
    for (int k : ans) cout << k << '\n';
}''', [("3 12\n6 1 5\n", "3\n1\n7\n11\n")])

P('TYP075', ['02_factorize'], r'''
void solv(){
    ll N; cin >> N;
    int c = 0;
    for (const PrimePower& f : factorize(N)) c += f.e;   // 素因数の数（重複も）
    cout << ilog_ceil(c, 2) << '\n';                      // 1 回で個数が最大 2 倍 → ⌈log2 c⌉
}''', [("42\n", "2\n"), ("48\n", "3\n"), ("54\n", "2\n"), ("53\n", "0\n")])

P('TYP085', ['02_factorize'], r'''
void solv(){
    ll K; cin >> K;
    const vector<ll> d = divisors_of(K);
    const ll a_max = kth_root_floor(K, 3);
    ll c = 0;
    for (size_t i = 0; i < d.size() && d[i] <= a_max; ++i) {
        const ll a = d[i], rest = K / a, b_max = isqrt(rest);
        for (size_t j = i; j < d.size() && d[j] <= b_max; ++j) c += rest % d[j] == 0;   // b ≤ c ⇔ b ≤ √(K/a)
    }
    cout << c << '\n';
}''', [("42\n", "5\n"), ("7\n", "1\n"), ("192\n", "16\n")])

P('ABC052_C', ['02_sieve', '00_modint'], r'''
void solv(){
    using mint = ModInt107;
    ll N; cin >> N;
    mint r = 1;
    for (int p : primes_up_to(N)) r *= factorial_exponent(N, p) + 1;
    cout << r << '\n';
}''', [("3\n", "4\n"), ("6\n", "30\n"), ("1000\n", "972926972\n")])

P('ABC280_D', ['02_factorize'], r'''
void solv(){
    ll K; cin >> K;
    cout << min_factorial_multiple(K) << '\n';
}''', [("30\n", "5\n"), ("123456789011\n", "123456789011\n"), ("280\n", "7\n")])

P('ABC254_D', ['02_sieve'], r'''
void solv(){
    int N; cin >> N;
    Sieve S(N);
    vector<ll> cnt(N + 1, 0);
    erep(i, 1, N) cnt[S.square_free_decompose(i).b]++;   // i·j が平方数 ⇔ 平方因子を除いた部分 b が同じ
    ll ans = 0;
    for (ll c : cnt) ans += c * c;
    cout << ans << '\n';
}''', [("4\n", "6\n"), ("254\n", "896\n")])

P('ABC284_D', ['02_factorize'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        ll N; cin >> N;
        const vector<PrimePower> fs = factorize(N);
        const ll p = fs[0].e == 2 ? fs[0].p : fs[1].p, q = fs[0].e == 2 ? fs[1].p : fs[0].p;
        cout << p << ' ' << q << '\n';
    }
}''', [("3\n2023\n63\n1059872604593911\n", "17 7\n3 7\n104149 97711\n")])

P('ABC361_F', ['02_factorize'], r'''
void solv(){
    ll N; cin >> N;
    cout << count_perfect_powers(N) << '\n';   // 1 も数える
}''', [("99\n", "12\n"), ("1000000000000000000\n", "1001003332\n")])

# ---------------- 数え上げ ----------------
P('ABC132_D', ['04_binomial'], r'''
void solv(){
    using mint = ModInt107;
    ll N, K; cin >> N >> K;
    Binomial<mint> B(N + 1);
    erep(i, 1, K) cout << B.C(N - K + 1, i) * B.C(K - 1, i - 1) << '\n';
}''', [("5 3\n", "3\n6\n1\n"), ("2000 3\n", "1998\n3990006\n327341989\n")])

P('ABC145_D', ['04_binomial'], r'''
void solv(){
    using mint = ModInt107;
    ll X, Y; cin >> X >> Y;
    if ((X + Y) % 3 != 0) { cout << 0 << '\n'; return; }
    const ll a = (2 * Y - X) / 3, b = (2 * X - Y) / 3;   // (+1, +2) を a 回、(+2, +1) を b 回
    Binomial<mint> B(X + Y);
    cout << B.paths(a, b) << '\n';                        // C(a + b, a)。負なら 0
}''', [("3 3\n", "2\n"), ("2 2\n", "0\n"), ("999999 999999\n", "151840682\n")])

P('ABC167_E', ['04_binomial'], r'''
void solv(){
    using mint = ModInt998;
    ll N, M, K; cin >> N >> M >> K;
    Binomial<mint> B(N);
    mint s = 0;
    erep(k, 0, K) s += mint(M) * mint(M - 1).pow(N - 1 - k) * B.C(N - 1, k);
    cout << s << '\n';
}''', [("3 2 1\n", "6\n"), ("100 100 0\n", "73074801\n"), ("60522 114575 7559\n", "479519525\n")])

P('ABC262_E', ['04_binomial'], r'''
void solv(){
    using mint = ModInt998;
    int N, M, K; cin >> N >> M >> K;
    vector<int> deg(N, 0);
    rep(i, M) { int u, v; cin >> u >> v; deg[u - 1]++, deg[v - 1]++; }
    int D = 0;
    for (int d : deg) D += d % 2;   // 赤どうしでない辺の数の偶奇 = 赤い頂点の次数の和の偶奇 → 奇数次数の赤が偶数個
    Binomial<mint> B(N);
    mint s = 0;
    for (int i = 0; i <= K; i += 2) s += B.C(D, i) * B.C(N - D, K - i);
    cout << s << '\n';
}''', [("4 4 2\n1 2\n1 3\n2 3\n3 4\n", "2\n"), ("10 10 3\n1 2\n2 4\n1 5\n3 6\n3 9\n4 10\n7 8\n9 10\n5 9\n3 4\n", "64\n")])

P('ARC175_A', ['00_modint'], r'''
void solv(){
    using mint = ModInt998;
    int N; cin >> N;
    vector<int> P(N); cin >> P;
    string S; cin >> S;
    // 全員が左を取る（最初の人が L）か、全員が右を取る（最初の人が R）のどちらかしかない
    auto count = [&](char hand, int dir) {
        vector<char> done(N, 0);
        mint r = 1;
        for (int k = 0; k < N; ++k) {
            const int i = P[k] - 1, nb = (i + dir + N) % N;   // 利き手の側のスプーンを先に取る相手
            if (done[nb]) { if (S[i] == '?') r *= 2; }         // 片方しか残っていない → どちらの利き手でもよい
            else if (S[i] != '?' && S[i] != hand) r = 0;       // 両方ある → 利き手が決まる
            done[i] = 1;
        }
        return r;
    };
    cout << count('L', 1) + count('R', -1) << '\n';
}''', [("3\n1 2 3\nL??\n", "2\n"), ("3\n1 3 2\nR?L\n", "0\n"), ("12\n6 2 9 3 1 4 11 5 12 10 7 8\n????????????\n", "160\n")])

# ---------------- 桁・和 ----------------
P('ABC044_D', ['06_quotient', '01_digits'], r'''
void solv(){
    ll n, s; cin >> n >> s;
    if (s == n) { cout << n + 1 << '\n'; return; }   // b ≤ n なら 2 桁以上で桁和 < n
    for (ll b = 2; b * b <= n; ++b)
        if (digit_sum(n, b) == s) { cout << b << '\n'; return; }
    ll best = -1;   // b > √n は 2 桁：n = q b + r、q + r = s → b = (n − s)/q + 1（i = b の小さい順に回る）
    if (n > s)
        for_each_quotient(n, [&](ll q, ll l, ll r) {
            if (best >= 0 || (n - s) % q != 0) return;
            const ll b = (n - s) / q + 1;
            if (l <= b && b < r && b > isqrt(n) && digit_sum(n, b) == s) best = b;
        });
    cout << best << '\n';
}''', [("87654\n30\n", "10\n"), ("87654\n138\n", "100\n"), ("87654\n45678\n", "-1\n"), ("31415926535\n1\n", "31415926535\n"), ("1\n31415926535\n", "-1\n")])

P('ABC372_G', ['08_fraction', '06_quotient'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        int N; cin >> N;
        vector<ll> A(N), B(N), C(N);
        ll xmax = LLONG_MAX;
        rep(i, N) {
            cin >> A[i] >> B[i] >> C[i];
            C[i] -= 1;                                     // A x + B y ≤ C − 1
            chmin(xmax, floor_div(C[i] - B[i], A[i]));     // y = 1 が取れる x の上限
        }
        if (xmax < 1) { cout << 0 << '\n'; continue; }
        // 直線 y = (C_i − A_i x) / B_i の下側の包絡線（x ∈ [1, xmax] ではどれも 1 以上）
        vector<int> id(N); iota(nall(id), 0);
        sort(nall(id), [&](int i, int j) {
            const Fraction si(A[i], B[i]), sj(A[j], B[j]);   // 傾きは −A/B。A/B の小さい順 = 左から
            if (si != sj) return si < sj;
            return Fraction(C[i], B[i]) < Fraction(C[j], B[j]);
        });
        auto cross = [&](int i, int j) { return Fraction(C[i] * B[j] - C[j] * B[i], A[i] * B[j] - A[j] * B[i]); };
        vector<int> h;
        for (int i : id) {
            if (!h.empty() && Fraction(A[h.back()], B[h.back()]) == Fraction(A[i], B[i])) continue;
            while (h.size() >= 2 && cross(h[h.size() - 2], i) <= cross(h[h.size() - 2], h.back())) h.pop_back();
            h.push_back(i);
        }
        Int128 ans = 0;
        ll lo = 1;
        for (size_t k = 0; k < h.size() && lo <= xmax; ++k) {
            ll hi = xmax;
            if (k + 1 < h.size()) chmin(hi, cross(h[k], h[k + 1]).floor());
            if (hi < lo) continue;
            const int i = h[k];
            ans += linear_floor_sum(hi - lo + 1, B[i], -A[i], C[i] - A[i] * lo);   // Σ_{x=lo}^{hi} ⌊(C − A x)/B⌋
            lo = hi + 1;
        }
        cout << ans << '\n';
    }
}''', [("2\n2\n1 1 4\n1 2 5\n1\n1 1 2\n", "2\n0\n"),
       ("3\n7\n138 16011 918976\n5478 7748 499926\n5234 17727 748589\n1157 10511 643136\n31200 3005 721285\n28839 14469 798851\n1933 5378 864127\n9\n17775 1665 386430\n37001 863 922418\n9756 4182 746671\n12379 9106 807578\n3984 4049 640539\n25333 9869 780810\n20372 7000 688738\n16107 11974 827227\n10779 10531 770510\n5\n4916 14132 460944\n11856 45422 610561\n56014 18216 825793\n10363 6220 945356\n37418 33866 851593\n",
        "660\n995\n140\n")])

P('ABC114_C', ['01_digits'], r'''
void solv(){
    ll N; cin >> N;
    ll c = 0;
    for (ll x : numbers_from_digits({3, 5, 7}, num_digits(N))) {
        if (x > N) break;   // 小さい順に並んでいる
        const vector<int> ds = digits_of(x);
        if (count(nall(ds), 3) && count(nall(ds), 5) && count(nall(ds), 7)) ++c;
    }
    cout << c << '\n';
}''', [("575\n", "4\n"), ("3600\n", "13\n"), ("999999999\n", "26484\n")])

P('ABC161_D', [], r'''
void solv(){
    // テンプレートに合う関数がないので BFS（小さい順に、末尾に ±1 の数字をつなげる）
    int K; cin >> K;
    deque<ll> q;
    erep(d, 1, 9) q.push_back(d);
    ll x = 0;
    rep(k, K) {
        x = q.front(); q.pop_front();
        const ll d = x % 10;
        for (ll e = max(0LL, d - 1); e <= min(9LL, d + 1); ++e) q.push_back(x * 10 + e);
    }
    cout << x << '\n';
}''', [("15\n", "23\n"), ("1\n", "1\n"), ("13\n", "21\n"), ("100000\n", "3234566667\n")])

P('ABC448_E', ['01_digits'], r'''
void solv(){
    ll K, M; cin >> K >> M;
    // ⌊N/M⌋ mod 10007 = ⌊(N mod 10007 M) / M⌋
    const ll mm = M * 10007;
    ll r = 0;
    rep(i, K) {
        ll c, l; cin >> c >> l;
        r = (mod_mul(r, pow10_mod(l, mm), mm) + repeat_concat_mod(c, l, mm)) % mm;   // 後ろに c を l 個つなげる
    }
    cout << r / M << '\n';
}''', [("6 7\n3 1\n1 1\n6 1\n2 2\n7 2\n6 2\n", "3797\n"), ("1 1\n1 1\n", "1\n"),
       ("10 9999\n9 419921892\n9 923650333\n6 476449815\n1 8837775\n2 141135534\n5 462618481\n3 202652735\n0 771538044\n4 321458589\n0 570032864\n", "8437\n")])

# ---------------- ビット ----------------
P('ABC147_D', ['09_bit'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> a(N); cin >> a;
    cout << xor_pair_sum(a) % 1000000007 << '\n';
}''', [("3\n1 2 3\n", "6\n"), ("10\n3 1 4 1 5 9 2 6 5 3\n", "237\n"),
       ("10\n3 14 159 2653 58979 323846 2643383 27950288 419716939 9375105820\n", "103715602\n")])

P('ABC197_C', ['09_bit'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> a(N); cin >> a;
    ll best = LLONG_MAX;
    for_each_submask((1LL << (N - 1)) - 1, [&](ll cut) {   // has_bit(cut, i) なら a[i] と a[i+1] のあいだで区切る
        ll x = 0, o = 0;
        rep(i, N) {
            o |= a[i];
            if (i == N - 1 || has_bit(cut, (int)i)) x ^= o, o = 0;
        }
        chmin(best, x);
    });
    cout << best << '\n';
}''', [("3\n1 5 7\n", "2\n"), ("3\n10 10 10\n", "0\n"), ("4\n1 3 3 1\n", "0\n")])

P('ABC238_D', ['09_bit'], r'''
void solv(){
    // x + y = 2 (x AND y) + (x XOR y)。x XOR y = s − 2a が 0 以上で a と共通のビットがなければよい
    int T; cin >> T;
    while (T--) {
        ll a, s; cin >> a >> s;
        const ll d = s - 2 * a;
        bool ok = d >= 0;
        if (ok) for_each_bit(a, [&](int i) { ok &= !has_bit(d, i); });
        if (ok) YES; else NO;
    }
}''', [("2\n1 8\n4 2\n", "Yes\nNo\n"),
       ("4\n201408139683277485 381410962404666524\n360288799186493714 788806911317182736\n18999951915747344 451273909320288229\n962424162689761932 1097438793187620758\n",
        "No\nYes\nYes\nNo\n")])

# ---------------- lcm ----------------
P('ABC445_E', ['02_sieve', '00_modint'], r'''
void solv(){
    using mint = ModInt998;
    Sieve S(10000000);
    int T; cin >> T;
    while (T--) {
        int N; cin >> N;
        vector<array<ll, 3>> fs;   // (p, e, 何番目の数か)
        rep(i, N) {
            ll a; cin >> a;
            for (const PrimePower& f : S.factorize(a)) fs.push_back({f.p, (ll)f.e, (ll)i});
        }
        sort(nall(fs));
        mint total = 1;
        vector<mint> adj(N, mint(1));
        for (size_t l = 0, r; l < fs.size(); l = r) {
            for (r = l; r < fs.size() && fs[r][0] == fs[l][0]; ++r) {}
            const ll p = fs[l][0], mx = fs[r - 1][1], sec = r - l >= 2 ? fs[r - 2][1] : 0;
            total *= mint(p).pow(mx);
            if (sec < mx) adj[fs[r - 1][2]] *= mint(p).pow(sec - mx);   // 最大の指数を持つ数が 1 つだけ → それを除くと 2 番目に下がる
        }
        rep(i, N) cout << total * adj[i] << (i + 1 < N ? ' ' : '\n');
    }
}''', [("3\n5\n9 12 25 8 15\n8\n592 26 167 912 171 321 327 651\n10\n9566582 2785103 1924635 359502 6912831 4893928 2820155 5071742 1836019 9037883\n",
        "600 1800 360 900 1800\n447582749 506009091 523568328 588652065 196217355 471970745 921220985 738745385\n747032704 756838459 344127037 146466685 159487731 555429485 826726159 884617928 322846201 298477407\n")])

# ================= 最終照合で足した機能を使う問題 =================
# ---------------- 最終照合で足した機能（公式の例） ----------------
P('ABC390_D', ['04_combinatorics'], r'''
void solv(){
    int N; cin >> N;
    vector<ll> A(N); cin >> A;
    vector<ll> xs, sum(N);
    for_each_set_partition(N, [&](const vector<int>& g, int k) {   // B_12 = 420 万通り
        fill(sum.begin(), sum.begin() + k, 0);
        rep(i, N) sum[g[i]] += A[i];
        ll x = 0; rep(j, k) x ^= sum[j];
        xs.push_back(x);
    });
    sort(nall(xs));
    cout << unique(nall(xs)) - xs.begin() << '\n';
}''', [("3\n2 5 7\n", "3\n"), ("2\n100000000000000000 100000000000000000\n", "2\n"), ("6\n71 74 45 34 31 60\n", "84\n")])

P('ABC209_E', ['10_game'], r'''
void solv(){
    int N; cin >> N;
    vector<string> s(N); cin >> s;
    map<string, int> id;   // 3 文字の状態
    auto get = [&](const string& t) { auto it = id.find(t); if (it != id.end()) return it->second; int k = id.size(); id[t] = k; return k; };
    vector<pair<int, int>> e(N);
    rep(i, N) e[i] = {get(s[i].substr(0, 3)), get(s[i].substr(s[i].size() - 3))};
    vector<vector<int>> mv(id.size());
    for (auto [a, b] : e) mv[a].push_back(b);   // 先頭 3 文字の状態から末尾 3 文字の状態へ
    const GameOutcome o = game_outcome(id.size(), mv);
    rep(i, N) {   // 高橋くんが単語 i を言った後、青木くんの手番
        const int r = o.result[e[i].second];
        cout << (r == GAME_WIN ? "Aoki" : r == GAME_LOSE ? "Takahashi" : "Draw") << '\n';
    }
}''', [("3\nabcd\nbcda\nada\n", "Aoki\nTakahashi\nDraw\n"), ("1\nABC\n", "Draw\n"),
       ("5\neaaaabaa\neaaaacaa\ndaaaaaaa\neaaaadaa\ndaaaafaa\n", "Takahashi\nTakahashi\nTakahashi\nAoki\nTakahashi\n")])

P('ABC172_D_min25', ['02_prime_count'], r'''
void solv(){
    ll N; cin >> N;   // Σ K·d(K)：f(K) = K d(K) は乗法的で f(p) = 2p、f(p^e) = p^e (e + 1)
    cout << multiplicative_sum<ll>(N, {0, 2}, [](ll p, int e) { return pow_int(p, e) * (e + 1); }) << '\n';
}''', [("4\n", "23\n"), ("100\n", "26879\n"), ("10000000\n", "838627288460105\n")])

# ---------------- 最終照合で足した機能（自作の例。期待値は Python の愚直解で作った） ----------------
P('SELF_multipoint', ['07_fps'], r'''
void solv(){
    using mint = ModInt998;
    int N, M; cin >> N >> M;
    vector<mint> c(N), p(M);
    for (auto& x : c) { ll v; cin >> v; x = v; }
    for (auto& x : p) { ll v; cin >> v; x = v; }
    for (const mint& v : multipoint_eval(c, p)) cout << v << ' ';
    cout << '\n';
}''', [("5 7\n1 2 3 4 5\n5 6 7 8 9 0 998244352\n", '3711 7465 13539 22737 35983 1 3\n')])

P('SELF_interpolate', ['07_fps'], r'''
void solv(){
    using mint = ModInt998;
    int N; cin >> N;
    vector<mint> x(N), y(N);
    for (auto& v : x) { ll t; cin >> t; v = t; }
    for (auto& v : y) { ll t; cin >> t; v = t; }
    cout << interpolate(x, y) << '\n';
}''', [("4\n1 5 2 9\n3 1 4 1\n", '802160641 133693444 142606335 918028289\n')])

P('SELF_set_exp', ['07_subset_transform'], r'''
void solv(){
    using mint = ModInt998;
    int N; cin >> N;
    vector<mint> b(1 << N);
    for (auto& v : b) { ll t; cin >> t; v = t; }
    for (const mint& v : subset_exp(b)) cout << v << ' ';
    cout << '\n';
}''', [("3\n0 1 2 3 4 5 6 7\n", '1 1 2 5 4 9 14 43\n')])

P('SELF_connected', ['07_subset_transform'], r'''
void solv(){
    using mint = ModInt998;
    int T; cin >> T;
    while (T--) {
        int n, m; cin >> n >> m;
        vector<pair<int, int>> es(m); for (auto& [u, v] : es) cin >> u >> v;
        vector<mint> all(1 << n);
        rep(s, 1 << n) { int e = 0; for (auto [u, v] : es) e += (s >> u & 1) && (s >> v & 1); all[s] = mint(2).pow(e); }
        cout << subset_log(all).back() << '\n';   // 全頂点が連結になる辺の選び方
    }
}''', [('3\n4 5\n0 1\n1 2\n2 3\n3 0\n0 2\n3 1\n0 1\n5 10\n0 1\n0 2\n0 3\n0 4\n1 2\n1 3\n1 4\n2 3\n2 4\n3 4\n', '14\n0\n728\n')])

P('SELF_totient_sum', ['02_prime_count', '00_modint'], r'''
void solv(){
    using mint = ModInt998;
    ll N; cin >> N;
    cout << multiplicative_sum<mint>(N, {-1, 1}, [](ll p, int e) { return mint(p).pow(e) - mint(p).pow(e - 1); }) << '\n';
}''', [("10\n", "32\n"), ("100000\n", "44917695\n"), ("10000000000\n", "866849765\n")])

P('SELF_floor_sum_ext', ['06_quotient', '00_modint'], r'''
void solv(){
    using mint = ModInt998;
    int T; cin >> T;
    while (T--) {
        ll n, m, a, b; cin >> n >> m >> a >> b;
        const auto r = linear_floor_sum_ext<mint>(n, m, a, b);
        cout << r.sum << ' ' << r.sum_i << ' ' << r.sum_sq << '\n';
    }
}''', [('6\n4 10 6 3\n6 5 4 3\n1 1 0 0\n31415 92653 58979 32384\n1000 7 -5 12345\n100000 999999937 123456789 -987654321\n', '3 8 5\n13 46 39\n0 0 0\n314095480 694303041 674212984\n1406357 642951786 23868045\n617129629 599495294 552365150\n')])

P('SELF_linear_mod_sum', ['06_quotient'], r'''
void solv(){
    int T; cin >> T;
    while (T--) { ll n, m, a, b; cin >> n >> m >> a >> b; cout << linear_mod_sum(n, m, a, b) << '\n'; }
}''', [('6\n4 10 6 3\n6 5 4 3\n1 1 0 0\n31415 92653 58979 32384\n1000 7 -5 12345\n100000 999999937 123456789 -987654321\n', '18\n13\n0\n1455586415\n3001\n49416607616627\n')])

P('SELF_relaxed_catalan', ['07_convolution'], r'''
void solv(){
    using mint = ModInt998;
    int Q; cin >> Q;
    RelaxedConvolution<mint> rc;
    vector<mint> C = {1};
    while (Q--) {
        int n; cin >> n;
        while ((int)C.size() <= n) { const mint c = C.back(); C.push_back(rc.push(c, c)); }   // C_{k+1} = Σ C_i C_{k−i}
        cout << C[n] << '\n';
    }
}''', [("2\n10\n199999\n", '16796\n566295761\n')])

P('SELF_gcd_sum', ['02_divisor_transform'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        ll M; cin >> M;
        const DivisorsOf D(M);
        vector<ll> c(D.size());
        rep(i, D.size()) c[i] = M / D[i];   // d の倍数の個数
        D.multiple_mobius(c);               // gcd(x, M) = d の個数
        Int128 s = 0; rep(i, D.size()) s += (Int128)D[i] * c[i];
        cout << s << '\n';
    }
}''', [('7\n1\n6\n12\n720720\n1000000000000\n999999999989\n963761198400\n', '1\n15\n40\n61916400\n74200000000000\n1999999999977\n1853278146720000\n')])

P('SELF_q_binomial', ['04_combinatorics'], r'''
void solv(){
    using mint = ModInt998;
    int T; cin >> T;
    while (T--) { ll n, k, q; cin >> n >> k >> q; cout << q_binomial<mint>(n, k, q) << '\n'; }
}''', [('6\n4 2 2\n10 3 1\n30 15 3\n50 20 998244352\n60 7 0\n40 13 779057549\n', '35\n120\n817230929\n3268760\n1\n0\n')])

P('SELF_binom_prefix_mo', ['04_binomial'], r'''
void solv(){
    using mint = ModInt998;
    int Q; cin >> Q;
    vector<ll> n(Q), m(Q); rep(i, Q) cin >> n[i] >> m[i];
    vector<int> ord(Q); iota(nall(ord), 0);
    const int B = 300;
    sort(nall(ord), [&](int a, int b) { return n[a] / B != n[b] / B ? n[a] < n[b] : ((n[a] / B) % 2 ? m[a] > m[b] : m[a] < m[b]); });
    BinomialPrefix<mint> P(100000);
    vector<mint> ans(Q);
    for (int i : ord) ans[i] = P.move_to(n[i], m[i]);
    rep(i, Q) cout << ans[i] << '\n';
}''', [('8\n5 2\n100 50\n7 7\n100000 49999\n3 0\n99999 50000\n0 0\n20 -1\n', '16\n41441083\n128\n413175703\n1\n988764691\n1\n0\n')])

P('SELF_big_mul', ['01_digits', '07_convolution'], r'''
void solv(){
    int T; cin >> T;
    while (T--) {
        string a, b; cin >> a >> b;
        bool neg = false;
        if (a[0] == '-') neg = !neg, a = a.substr(1);
        if (b[0] == '-') neg = !neg, b = b.substr(1);
        vector<ll> x(a.rbegin(), a.rend()), y(b.rbegin(), b.rend());   // 下の桁から
        for (auto& v : x) v -= '0';
        for (auto& v : y) v -= '0';
        const string s = carry_to_string(convolve_ll(x, y));
        cout << (neg && s != "0" ? "-" : "") << s << '\n';
    }
}''', [('5\n0 0\n12 -34\n-99999999999999999999 99999999999999999999\n888635009363308772389590261263111804134211405022045064349843702920313005978970156998061486981438338347608813582822017835137570863827918927183799571006424549234097440696641880988643944993252504042052940806197660519026522591697946710778467916664723035798888265339643789912985122130532920872122203515185 -5791572010985296464704062486791447554855212774242921280027526934055334555611155636965156875542968830575318093059131006498609683231669862267293035312455773611688481699699972618357987656249527506389205594896592448943614286293925767905175898235987966342\n100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001\n', '0\n-408\n-9999999999999999999800000000000000000001\n-5146593648210195940296941371817942141949234207766308273703129674324840083110755134974856553088534056983753980964217497739636108403884197870436877639336939772189939946595581012800266857305129460844282902095340607221059834080315215629969154322673586278049604986750328495500091117963145681427770851704747112271361464559010735039899278704872618832080238080524690027330702107464661028829678588846990861375742643525040597853858867610932671179622058775579877941488443488443819377341433206987481763773891136371851261096473497913968265784662860382351865903270\n1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000\n')])

P('SELF_next_palindrome', ['01_digits'], r'''
void solv(){
    int T; cin >> T;
    while (T--) { string s; cin >> s; cout << ceil_palindrome_str(s) << '\n'; }
}''', [('10\n0\n9\n10\n99\n123456\n808\n999999\n1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000\n999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999\n12345678901234567890123456789\n', '0\n9\n11\n99\n124421\n808\n999999\n1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001\n999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999\n12345678901234643210987654321\n')])

P('SELF_perm_parity', ['04_permutation'], r'''
void solv(){
    int T; cin >> T;
    while (T--) { int n; cin >> n; vector<int> p(n); cin >> p; cout << (perm_parity(p) ? "Odd" : "Even") << '\n'; }
}''', [('5\n1\n1\n2\n2 1\n3\n3 1 2\n4\n4 3 2 1\n5\n2 3 4 5 1\n', 'Even\nOdd\nEven\nEven\nEven\n')])
