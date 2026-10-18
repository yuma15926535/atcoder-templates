// =============================================================================
//  ランダムテスト
// =============================================================================
namespace mint_test {

{COMMON_TEST}

static long long refmod(Int128 x, long long m) {
    x %= m;
    if (x < 0) x += m;
    return (long long)x;
}
static long long ref_pow(long long a, unsigned long long e, long long m) {
    Int128 r = 1 % m, x = refmod(a, m);
    for (; e; e >>= 1, x = x * x % m)
        if (e & 1) r = r * x % m;
    return (long long)r;
}
static bool brute_is_prime(long long n) {
    if (n < 2) return false;
    for (long long d = 2; d * d <= n; ++d)
        if (n % d == 0) return false;
    return true;
}
static std::string rnd_digits(std::mt19937_64& rng, int len) {
    std::string s;
    for (int i = 0; i < len; ++i) s += char('0' + rng() % 10);
    return s;
}

// M の全機能を、128 ビットの愚直計算と突き合わせる（M::mod() は呼ぶ前に決まっていること）
template <class M>
static void test_ops(std::mt19937_64& rng, int& cases, int iters) {
    const long long m = M::mod();
    for (int it = 0; it < iters; ++it, ++cases) {
        const long long x = rnd_ll(rng), y = it % 7 == 0 ? x / 4 * (long long)rnd(rng, -2, 2) : rnd_ll(rng);  // 7 回に 1 回は x と関係のある値
        const M a = x, b = y;
        const long long rx = refmod(x, m), ry = refmod(y, m);
        MINT_CHECK(a.val() == rx && b.val() == ry && 0 <= a.val() && a.val() < m);
        MINT_CHECK((a + b).val() == refmod((Int128)x + y, m) && (a - b).val() == refmod((Int128)x - y, m));
        MINT_CHECK((a * b).val() == refmod((Int128)rx * ry, m) && (-a).val() == refmod(-(Int128)x, m) && (+a) == a);
        MINT_CHECK((a * 3 + 1).val() == refmod((Int128)rx * 3 + 1, m) && (2 - a).val() == refmod(2 - (Int128)x, m) && (a == b) == (rx == ry));
        MINT_CHECK((a != b) == (rx != ry) && (a < b) == (rx < ry) && (a <= b) == (rx <= ry) && (a > b) == (rx > ry) && (a >= b) == (rx >= ry));
        MINT_CHECK((bool)a == (rx != 0) && (a == x) && (b == y));
        M c = a;
        c += b, c -= a, c *= a, c -= b * a;
        MINT_CHECK(c == 0);
        c = a;
        MINT_CHECK((c++).val() == rx && c.val() == refmod((Int128)x + 1, m) && (--c).val() == rx && (++c).val() == refmod((Int128)x + 1, m) && (c--) == a + 1 && c == a);
        // 別の型から
        const unsigned long long ux = (unsigned long long)x;
        const Int128 big = (Int128)x * (long long)(rng() >> 1) + (Int128)y;
        const UInt128 ubig = ((UInt128)rng() << 64) | rng();
        MINT_CHECK(M(ux).val() == (long long)(ux % (unsigned long long)m) && M(big).val() == refmod(big, m) && M(ubig).val() == (long long)(ubig % (UInt128)m));
        MINT_CHECK(M((int)x).val() == refmod((int)x, m) && M((unsigned)x).val() == (long long)((unsigned)x % (unsigned long long)m) && M(true).val() == 1 % m && M((char)7).val() == 7 % m);
        MINT_CHECK(M::raw((int)rx) == a);
        // 累乗・逆元
        const unsigned long long e = it % 3 ? (unsigned long long)rnd(rng, 0, 30) : rng() >> 1;
        MINT_CHECK(a.pow((long long)e).val() == ref_pow(rx, e, m));
        const bool unit = std::gcd(rx, m) == 1;
        MINT_CHECK(a.is_unit() == unit);
        const Maybe<M> inv = a.try_inv();
        MINT_CHECK(inv.ok == unit);
        if (unit) {
            MINT_CHECK((a * inv.value).val() == 1 % m && a.inv() == inv.value && (b / a * a) == b && (b / a).val() == refmod((Int128)ry * inv.value.val(), m));
            MINT_CHECK(a.pow(-(long long)e) * a.pow((long long)e) == M(1) && a.pow(-1) == inv.value);
        }
        // 入出力
        std::ostringstream os;
        os << a << ' ' << -b;
        MINT_CHECK(os.str() == std::to_string(rx) + " " + std::to_string(refmod(-(Int128)y, m)));
        const std::string d = rnd_digits(rng, (int)rnd(rng, 1, 36));
        Int128 dv = 0;
        for (char ch : d) dv = dv * 10 + (ch - '0');
        std::istringstream is(d + " -" + d + " +" + d);
        M r1, r2, r3;
        is >> r1 >> r2 >> r3;
        MINT_CHECK(r1.val() == refmod(dv, m) && r2.val() == refmod(-dv, m) && r3 == r1);
    }
    // 長い計算の連鎖（遅延簡約の値がずれないか）
    {
        M acc = 1, acc2 = 0;
        long long ra = 1 % m, rb = 0;
        for (int i = 0; i < 20000; ++i, ++cases) {
            const long long v = rnd_ll(rng);
            const long long rv = refmod(v, m);
            switch (rng() % 6) {
                case 0: acc *= v, ra = refmod((Int128)ra * rv, m); break;
                case 1: acc += v, ra = refmod((Int128)ra + rv, m); break;
                case 2: acc -= v, ra = refmod((Int128)ra - rv, m); break;
                case 3: acc = -acc, ra = refmod(-(Int128)ra, m); break;
                case 4: acc2 += acc * acc, rb = refmod((Int128)rb + (Int128)ra * ra, m); break;
                default: acc2 = acc2 * acc2 - acc, rb = refmod((Int128)rb * rb - ra, m); break;
            }
            MINT_CHECK(acc.val() == ra && acc2.val() == rb);
        }
    }
    // ハッシュ・map・sort
    {
        std::unordered_set<M> hs;
        std::map<M, int> mp;
        std::vector<long long> vals;
        for (int i = 0; i < 2000; ++i) {
            const long long v = rnd(rng, -100000, 100000);
            hs.insert(M(v)), mp[M(v)]++, vals.push_back(refmod(v, m));
        }
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
        MINT_CHECK(hs.size() == vals.size() && mp.size() == vals.size() && mp.begin()->first.val() == vals[0]);
        ++cases;
    }
}

template <unsigned MOD>
static void test_static(std::mt19937_64& rng, int& cases) {
    static_assert(ModInt<MOD>::mod() == (int)MOD, "mod");
    MINT_CHECK(ModInt<MOD>::IS_PRIME == brute_is_prime(MOD));
    test_ops<ModInt<MOD>>(rng, cases, 3000 * TEST_SCALE);
}

static void test_dynamic(std::mt19937_64& rng, int& cases) {
    std::vector<long long> mods = {998244353, 1000000007, 1, 2, 3, 4, 6, 7, 2147483647, 2147483646, 1073741823, 1073741824, 1073741825,
                                   1073741789, 1000000000, 999999999, 65536, 167772161};
    for (int i = 0; i < 60; ++i) mods.push_back(rnd(rng, 1, i % 3 == 0 ? 100 : 2147483647));
    for (long long m : mods) {
        DynModInt::set_mod(m);
        MINT_CHECK(DynModInt::mod() == m && DynModInt::mod_is_prime() == brute_is_prime(m));
        test_ops<DynModInt>(rng, cases, 500 * TEST_SCALE);
    }
    // 2 つの mod を同時に
    DynModInt::set_mod(1000000007);
    DynModIntT<1>::set_mod(998244353);
    for (int it = 0; it < 2000 * TEST_SCALE; ++it, ++cases) {
        const long long x = rnd_ll(rng), y = rnd_ll(rng);
        const DynModInt a = x, b = y;
        const DynModIntT<1> c = x, d = y;
        MINT_CHECK((a * b).val() == refmod((Int128)refmod(x, 1000000007) * refmod(y, 1000000007), 1000000007));
        MINT_CHECK((c * d).val() == refmod((Int128)refmod(x, 998244353) * refmod(y, 998244353), 998244353));
        MINT_CHECK((a * b).val() == (ModInt107(x) * ModInt107(y)).val() && (c + d).val() == (ModInt998(x) + ModInt998(y)).val());
    }
    DynModInt::set_mod(998244353);
}

template <class M>
static void test_reconstruct(std::mt19937_64& rng, int& cases) {
    const long long m = M::mod(), B = isqrt(m / 2);
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const long long lim = it % 2 ? std::min<long long>(B, 1000) : B;
        long long p = rnd(rng, -lim, lim), q = rnd(rng, 1, std::max<long long>(1, lim));
        const long long g = std::gcd(std::llabs(p), q);
        p /= g, q /= g;
        if (std::gcd(q, m) != 1) continue;
        const RationalGuess r = rational_reconstruct(M(p) / M(q));
        MINT_CHECK(r.ok && r.num == p && r.den == q);
        // 小さい分数でないかもしれない値：見つかったら条件を満たすこと
        const M z = rnd_ll(rng);
        const RationalGuess s = rational_reconstruct(z);
        if (s.ok) MINT_CHECK(s.den >= 1 && s.den <= B && std::llabs(s.num) <= B && M(s.num) == z * M(s.den));
    }
    std::ostringstream os;
    os << rational_reconstruct(M(1) / M(3)) << ' ' << rational_reconstruct(M(5));
    MINT_CHECK(os.str() == "1/3 5");
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_static<998244353>(rng, cases);
    test_static<1000000007>(rng, cases);
    test_static<167772161>(rng, cases);
    test_static<2147483647>(rng, cases);   // 素数・2^30 以上（モンゴメリ・遅延簡約なし）
    test_static<1073741823>(rng, cases);   // 合成数・2^30 未満（遅延簡約）
    test_static<1073741825>(rng, cases);   // 合成数・2^30 以上
    test_static<999999999>(rng, cases);
    test_static<1000000000>(rng, cases);   // 偶数
    test_static<65536>(rng, cases);
    test_static<7>(rng, cases);
    test_static<3>(rng, cases);
    test_static<2>(rng, cases);
    test_static<1>(rng, cases);
    test_dynamic(rng, cases);
    test_reconstruct<ModInt998>(rng, cases);
    test_reconstruct<ModInt107>(rng, cases);
    DynModInt::set_mod(1000000000);
    test_reconstruct<DynModInt>(rng, cases);
    DynModInt::set_mod(998244353);
    std::cout << "[modint] random_test: OK (" << cases << " cases)\n";
}

template <class M>
static void bench_one(const char* name) {
    std::mt19937_64 rng(1);
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    unsigned long long chk = 0;
    auto t0 = std::chrono::steady_clock::now();
    {
        M x = 12345, y = 678901, c = 3;
        for (int i = 0; i < 20000000; ++i) x = x * y + c;
        chk += (unsigned long long)x.val();
    }
    const long long t_chain = ms(t0);
    const int n = 1 << 20;
    std::vector<M> a(n), b(n);
    for (int i = 0; i < n; ++i) a[i] = (long long)(rng() >> 2), b[i] = (long long)(rng() >> 2);
    t0 = std::chrono::steady_clock::now();
    {
        M s = 0;
        for (int r = 0; r < 16; ++r)
            for (int i = 0; i < n; ++i) s += a[i] * b[i];
        chk += (unsigned long long)s.val();
    }
    const long long t_dot = ms(t0);
    t0 = std::chrono::steady_clock::now();
    {
        const int N = 10000000;
        std::vector<M> f(N + 1), g(N + 1);
        f[0] = 1;
        for (int i = 1; i <= N; ++i) f[i] = f[i - 1] * i;
        g[N] = f[N].inv();
        for (int i = N; i > 0; --i) g[i - 1] = g[i] * i;
        chk += (unsigned long long)g[7].val();
    }
    const long long t_fact = ms(t0);
    t0 = std::chrono::steady_clock::now();
    {
        M s = 0;
        for (int i = 1; i < 1000000; ++i) s += a[i].pow(1000000000 + i) + a[i].inv();
        chk += (unsigned long long)s.val();
    }
    const long long t_pow = ms(t0);
    std::cout << "[modint] benchmark: " << name << "  x = x*y+c × 2×10^7 " << t_chain << " ms, 内積 1.6×10^7 回 " << t_dot
              << " ms, 階乗と逆元の表 N = 10^7 " << t_fact << " ms, pow + inv × 10^6 " << t_pow << " ms  (check " << chk % 7 << ")\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[modint] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    bench_one<ModInt998>("ModInt998    ");
    bench_one<ModInt107>("ModInt107    ");
    DynModInt::set_mod(998244353);
    bench_one<DynModInt>("DynModInt(998244353)");
    DynModInt::set_mod(2147483647);
    bench_one<DynModInt>("DynModInt(2^31-1)");
    DynModInt::set_mod(998244353);
}

}  // namespace mint_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace mint_demo {

{COMMON_DEMO}

static void demo() {
    banner("modint 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧（ModInt と DynModInt で同じ）\n"
              << "  using mint = ModInt998;  / ModInt107 / ModInt<M> / DynModInt（mint::set_mod(M)）   使う側で 1 行書く\n"
              << "  mint x = 整数 / mint::raw(x)            作る（負の数・long long・Int128 から）\n"
              << "  + - * / ++ -- 比較                      int と同じように（整数と混ざってもよい）\n"
              << "  x.val() / mint::mod()                   値 / mod\n"
              << "  x.pow(e) / x.inv() / x.try_inv() / x.is_unit()   累乗（負の指数も）/ 逆元 / 止まらない逆元 / 逆元があるか\n"
              << "  cout << x / cin >> x / unordered_map<mint, V>    入出力 / ハッシュ\n"
              << "  rational_reconstruct(x)                 mint の値を小さい分数に戻す（デバッグ用）\n";

    api("【1】 using mint = ModInt998;  mint a = x;  + - * /  cout << a", "mint",
        "998244353 で割った余りで計算する整数。int・long long と同じように書ける。val() なしで cout できる", "O(1)（/ は O(log mod)）");
    using mint = ModInt998;
    RUN(mint a = -3; mint b = 10);
    SHOW(a, "負の数は mod を足した値（998244350）");
    SHOW(a * b + 5, "整数と混ざってもよい");
    SHOW(b / 4, "割り算は逆元を掛ける（10 / 4 = 5/2）");
    SHOW(b / 4 * 4, "掛け戻すと 10");
    SHOW(mint(1000000000000000000LL), "long long から（10^18 mod 998244353）");
    SHOW(mint::mod(), "mod");
    RUN(mint c = 5; c *= c; c -= 1; ++c);
    SHOW(c, "複合代入・++ も");
    SHOW(c.val() + 1LL, "値を整数として使うとき val()");

    api("【2】 x.pow(e) / x.inv()", "mint", "累乗（e は long long。負なら逆元の累乗）/ 逆元", "O(log e) / O(log mod)");
    SHOW(mint(2).pow(10), "1024");
    SHOW(mint(2).pow(-1), "2 の逆元（499122177）");
    SHOW(mint(3).inv() * 3, "逆元を掛けると 1");
    {
        // 典型90 069: N 個のブロックを K 色で、距離 2 以内の色がすべて違う塗り方。K × (K-1) × (K-2)^(N-2)
        using mint107 = ModInt107;
        const long long N = 1000000000000000000LL, K = 3;
        const mint107 ans = N == 1 ? mint107(K) : mint107(K) * (K - 1) * mint107(K - 2).pow(N - 2);
        note("典型90 069（N ≤ 10^18 の累乗）: K × (K − 1) × (K − 2)^(N − 2) を mod 10^9+7 で（ModInt107）");
        std::cout << "    N = 10^18, K = 3 → " << ans << "\n";
        note("注意：pow の指数が負だと逆元の累乗になる。N − 2 のような式が負になりうるとき（N = 1）は、上のように場合分けする");
    }

    api("【3】 DynModInt::set_mod(M)", "なし", "mod が入力で決まるとき。素数でなくてよい・偶数でもよい（逆元は mod と互いに素なときだけ）", "O(1)");
    {
        using dmint = DynModInt;
        RUN(dmint::set_mod(100); dmint x = 7; dmint y = 15);
        SHOW(x * y, "105 mod 100");
        SHOW(x.inv(), "7 × 43 = 301 ≡ 1 (mod 100)");
        SHOW(y.is_unit(), "15 と 100 は互いに素でない");
        SHOW(y.try_inv().ok, "逆元はない（inv() を呼ぶと止まる）");
        SHOW(dmint::mod_is_prime(), "100 は素数でない");
        dmint::set_mod(998244353);
        note("ABC425 E のように mod が素数でないときは、割り算を使わない形（パスカルの三角形など。04_binomial）で数える。");
        note("2 つの mod を同時に使うなら DynModIntT<1>::set_mod(M2) で別の型を作る。");
    }

    api("【4】 期待値・確率を mod で", "mint", "確率 p/100 は mint(p) / 100。分数のまま計算した値と同じものが mod で求まる", "-");
    {
        // ABC280 E: 体力 N。確率 P% で 2、それ以外で 1 のダメージ。体力が 0 以下になるまでの攻撃回数の期待値
        const int N = 3, P = 10;
        const mint p = mint(P) / 100;
        std::vector<mint> dp(N + 1);
        for (int i = 1; i <= N; ++i) dp[i] = 1 + p * (i >= 2 ? dp[i - 2] : mint(0)) + (1 - p) * dp[i - 1];
        note("ABC280 E: dp[i] = 1 + p·dp[i−2] + (1 − p)·dp[i−1]（p = P/100 を mint で）");
        std::cout << "    N = 3, P = 10 → " << dp[N] << "\n";
        SHOW(rational_reconstruct(dp[N]), "分数に戻すと 281/100 = 2.81（手計算と合うか確かめられる）");
    }

    api("【5】 rational_reconstruct(x)", "RationalGuess{ok, num, den}", "mint の値を |分子|, 分母 ≤ √(mod/2) の分数に戻す（デバッグ用。見つからなければ ok = false）", "O(log mod)");
    SHOW(rational_reconstruct(mint(1) / 3), "1/3");
    SHOW(rational_reconstruct(mint(-7) / 12), "-7/12");
    SHOW(rational_reconstruct(mint(123456789)), "小さい分数がない値");

    api("【6】 比較・ハッシュ・入力", "-", "比較は 0〜mod-1 の値で（sort・map に使える）。unordered_map<mint, V> もそのまま。cin は長い数も読める", "-");
    {
        RUN(std::unordered_map<mint, int> cnt; cnt[mint(5)]++; cnt[mint(998244358)]++);
        SHOW(cnt[mint(5)], "998244358 ≡ 5 なので同じキー");
        std::istringstream is("123456789012345678901234567890");
        mint big;
        is >> big;
        SHOW(big, "30 桁の数を読んで mod を取る（long long に入らない入力）");
        SHOW(mint(3) < mint(-1), "値での比較（-1 は 998244352）");
        SHOW(ModInt998::IS_PRIME, "mod が素数か（コンパイル時）");
        SHOW(mint::raw(5), "0 ≤ x < mod と分かっているとき（剰余を取らないので少し速い）");
    }
    note("組合せ（nCr）は 04_binomial の Binomial<mint>、行列は 05_matrix の Matrix<mint>、畳み込みは 07_convolution。");
    note("mod が 2^31 以上（10^18 など）のときは modint ではなく 00_math_base の mod_mul / mod_pow を使う。");
}

}  // namespace mint_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    mint_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mint_test::benchmark();    // 速度計測
    mint_demo::demo();         // 使い方デモ
    return 0;
}
