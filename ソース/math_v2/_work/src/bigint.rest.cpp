// =============================================================================
//  ランダムテスト
// =============================================================================
namespace big_test {

{COMMON_TEST}

static Int128 rnd128(std::mt19937_64& rng) {  // 大きさがばらばらの Int128
    const int bits = (int)rnd(rng, 0, 126);
    UInt128 u = ((UInt128)rng() << 64) | rng();
    u = bits == 0 ? 0 : u >> (128 - bits);
    return (rng() & 1) ? -(Int128)u : (Int128)u;
}
static std::string rand_digits(std::mt19937_64& rng, int len) {  // len 桁（先頭は 0 でない。符号もランダム）
    std::string s = (rng() & 1) ? "-" : "";
    s += (char)('1' + rng() % 9);
    for (int i = 1; i < len; ++i) s += (char)('0' + (rng() % 5 == 0 ? (rng() % 2) * 9 : rng() % 10));  // 0 と 9 を多めに
    return s;
}
static const long long PRIMES[3] = {998244353, 1000000007, 4611686018427387847LL};

static void check_div(const BigInt& a, const BigInt& b) {
    const auto qr = divmod(a, b);
    const BigInt& q = qr.first;
    const BigInt& r = qr.second;
    BIG_CHECK(q * b + r == a && abs(r) < abs(b));
    BIG_CHECK(r.is_zero() || r.sign() == a.sign());
    BIG_CHECK(a / b == q && a % b == r);
}
static void check_mul(const BigInt& a, const BigInt& b) {
    const BigInt c = a * b;
    for (long long p : PRIMES) BIG_CHECK(c.mod_ll(p) == (long long)((Int128)a.mod_ll(p) * b.mod_ll(p) % p));
    BIG_CHECK(c.sign() == a.sign() * b.sign());
}

static void test_small(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const Int128 x = rnd128(rng), y = rnd128(rng);
        const BigInt a = x, b = y;
        BIG_CHECK(a.to_string() == int128_to_string(x) && BigInt(int128_to_string(x)) == a && a.to_i128() == x);
        BIG_CHECK((a < b) == (x < y) && (a == b) == (x == y) && (a >= b) == (x >= y) && compare(a, b) == (x < y ? -1 : x > y ? 1 : 0));
        // 足し算・引き算（Int128 で溢れない範囲だけ比べる：|x|, |y| < 2^126）
        BIG_CHECK((a + b).to_i128() == x + y && (a - b).to_i128() == x - y && (-a).to_i128() == -x);
        if (y != 0) BIG_CHECK((a / b).to_i128() == x / y && (a % b).to_i128() == x % y);
        const Int128 xs = x >> 64, ys = y >> 64;  // 掛け算は小さくして
        BIG_CHECK((BigInt(xs) * BigInt(ys)).to_i128() == xs * ys);
        const long long l = (long long)(x >> 64);
        BIG_CHECK(BigInt(l).fits_ll() && BigInt(l).to_ll() == l);
        BIG_CHECK(a.fits_ll() == (x >= LLONG_MIN && x <= LLONG_MAX));
        const long long m = rnd(rng, 1, 1000000000000000000LL);
        BIG_CHECK(a.mod_ll(m) == (long long)(((x % m) + m) % m));
        BIG_CHECK(a.num_digits() == (int)int128_to_string(x < 0 ? -x : x).size());
        if (x >= 0) {  // 平方根
            const BigInt s = isqrt(a);
            BIG_CHECK(s * s <= a && (s + 1) * (s + 1) > a);
        }
    }
    // 境目
    BIG_CHECK(BigInt(LLONG_MIN).to_ll() == LLONG_MIN && BigInt(LLONG_MAX).to_ll() == LLONG_MAX);
    BIG_CHECK(!(BigInt(LLONG_MAX) + 1).fits_ll() && !(BigInt(LLONG_MIN) - 1).fits_ll());
    const Int128 I128MIN = (Int128)((UInt128)1 << 127), I128MAX = (Int128)(((UInt128)1 << 127) - 1);
    BIG_CHECK(BigInt(I128MIN).to_i128() == I128MIN && BigInt(I128MAX).to_i128() == I128MAX && !(BigInt(I128MAX) + 1).fits_i128() && !(BigInt(I128MIN) - 1).fits_i128());
    BIG_CHECK(BigInt(I128MIN).to_string() == "-170141183460469231731687303715884105728");
    BIG_CHECK(BigInt((UInt128)-1).to_string() == "340282366920938463463374607431768211455" && BigInt(18446744073709551615ULL).to_string() == "18446744073709551615");
    BIG_CHECK(BigInt("-0") == BigInt(0) && BigInt("-0").sign() == 0 && BigInt("+000123") == 123 && BigInt("0000").to_string() == "0");
    BIG_CHECK((BigInt(-7) / 2) == -3 && (BigInt(-7) % 2) == -1 && (BigInt(7) / -2) == -3 && (BigInt(7) % -2) == 1);
    BIG_CHECK(BigInt(5) - 5 == 0 && (BigInt(5) - 5).sign() == 0 && -(BigInt(0)) == 0);
    BigInt t = 41;
    BIG_CHECK((++t) == 42 && (--t) == 41 && t.is_even() == false);
    cases += 10;
}

static void test_big(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int la = (int)rnd(rng, 1, it % 10 == 0 ? 2000 : 120), lb = (int)rnd(rng, 1, it % 3 == 0 ? la + 5 : 60);
        const BigInt a(rand_digits(rng, la)), b(rand_digits(rng, lb));
        BIG_CHECK(a.num_digits() == la && BigInt(a.to_string()) == a);
        check_mul(a, b);
        check_div(a, b);
        check_div(a * b + b / 3, b);
        BIG_CHECK((a * b) / b == a && (a * b) % b == 0 && (a + b) - b == a);
        // 9 で割った余りと各桁の和
        BIG_CHECK(a.digit_sum() % 9 == abs(a).mod_ll(9));
    }
    // 長い掛け算（カラツバ・10^3 進の方も）
    for (int len : {500, 3000, 20000, 100000, 120000}) {
        for (int t = 0; t < 2; ++t, ++cases) {
            const BigInt a(rand_digits(rng, len)), b(rand_digits(rng, len - t * 77));
            check_mul(a, b);
            if (len <= 20000) check_div(a * b + a, b);
        }
    }
    // 長さが大きく違う掛け算（短い方で区切る）
    for (int t = 0; t < 6; ++t, ++cases) {
        const BigInt a(rand_digits(rng, 50000 + t * 777)), b(rand_digits(rng, 300 + t * 401));
        check_mul(a, b);
        check_mul(b, a);
        BIG_CHECK(a * b == b * a && (a * b) / b == a);
    }
    // 大きい平方根・int との gcd・後置の ++ −−
    for (int len : {40, 41, 100, 1000, 5000, 20000}) {
        const BigInt n = abs(BigInt(rand_digits(rng, len)));
        const BigInt s = isqrt(n);
        BIG_CHECK(s * s <= n && (s + 1) * (s + 1) > n && isqrt(s * s) == s && isqrt(s * s - 1) == s - 1);
        ++cases;
    }
    {
        BigInt x = 5;
        BIG_CHECK(x++ == 5 && x == 6 && x-- == 6 && x == 5);
        BIG_CHECK(gcd(BigInt(12), 18) == 6 && gcd(30LL, BigInt(-45)) == 15 && gcd(BigInt(7), 0) == 7);
        cases += 2;
    }
    {  // 全部 9（繰り上がりがいちばん多い）
        const BigInt a(std::string(150000, '9'));
        const BigInt c = a * a;  // (10^n − 1)^2 = 10^{2n} − 2·10^n + 1
        BIG_CHECK(c == pow(BigInt(10), 300000) - 2 * pow(BigInt(10), 150000) + 1);
        ++cases;
    }
    // gcd・isqrt・pow
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {
        const BigInt g(rand_digits(rng, (int)rnd(rng, 1, 40))), a(rand_digits(rng, (int)rnd(rng, 1, 60))), b(rand_digits(rng, (int)rnd(rng, 1, 60)));
        const BigInt h = gcd(a * g, b * g);
        BIG_CHECK(h.sign() > 0 && (a * g) % h == 0 && (b * g) % h == 0 && h % abs(g) == 0);
        BIG_CHECK(gcd((a * g) / h, (b * g) / h) == 1);
        const BigInt l = lcm(a * g, b * g);
        BIG_CHECK(l * h == abs(a * g * b * g) && (a.is_zero() || b.is_zero() || g.is_zero() || (l % (a * g) == 0 && l % (b * g) == 0)));
        const BigInt n = abs(a) * abs(b) + abs(g);
        const BigInt s = isqrt(n);
        BIG_CHECK(s * s <= n && (s + 1) * (s + 1) > n && isqrt(s * s) == s);
        const long long e = rnd(rng, 0, 30);
        BigInt p = 1;
        for (int i = 0; i < e; ++i) p *= g;
        BIG_CHECK(pow(g, e) == p);
    }
    BIG_CHECK(gcd(BigInt(0), BigInt(0)) == 0 && gcd(BigInt(0), BigInt(-12)) == 12 && pow(BigInt(0), 0) == 1);
    BIG_CHECK(lcm(BigInt(4), 6) == 12 && lcm(-4LL, BigInt(6)) == 12 && lcm(BigInt(0), 5) == 0 && lcm(BigInt(1000000007), BigInt(998244353)) == BigInt(1000000007LL * 998244353LL));
    // よく知られた値
    BIG_CHECK(pow(BigInt(2), 1000).digit_sum() == 1366);  // Project Euler 16
    BigInt f = 1;
    for (int i = 1; i <= 100; ++i) f *= i;
    BIG_CHECK(f.digit_sum() == 648 && f.num_digits() == 158);  // Project Euler 20
    BIG_CHECK(isqrt(pow(BigInt(10), 100)) == pow(BigInt(10), 50) && isqrt(pow(BigInt(10), 101)).to_string().substr(0, 10) == "3162277660");
    cases += 4;
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_small(rng, cases);
    test_big(rng, cases);
    std::cout << "[bigint] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[bigint] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const BigInt a(rand_digits(rng, 100000)), b(rand_digits(rng, 100000)), c(rand_digits(rng, 50000));
    auto t0 = std::chrono::steady_clock::now();
    const BigInt p = a * b;
    const long long t_mul = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const BigInt q = p / c;
    const long long t_div = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const BigInt sq = isqrt(abs(p));
    const long long t_sqrt = ms(t0);
    t0 = std::chrono::steady_clock::now();
    BigInt f = 1;
    for (int i = 1; i <= 10000; ++i) f *= i;
    const long long t_fact = ms(t0);
    t0 = std::chrono::steady_clock::now();
    const std::string s = f.to_string();
    const BigInt g(s);
    const long long t_str = ms(t0);
    std::cout << "[bigint] benchmark: 10^5 桁 × 10^5 桁 " << t_mul << " ms, 2×10^5 桁 ÷ 5×10^4 桁 " << t_div << " ms, isqrt（2×10^5 桁）" << t_sqrt << " ms, 10000!（35660 桁）" << t_fact
              << " ms, 文字列との変換 " << t_str << " ms  (check " << (p.mod_ll(7) + q.mod_ll(7) + g.mod_ll(7) + sq.mod_ll(7)) % 7 << ")\n";
}

}  // namespace big_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace big_demo {

{COMMON_DEMO}

static void demo() {
    banner("多倍長整数（bigint） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  BigInt(long long / Int128 / \"文字列\")、+ − × ÷ %、比較、divmod / pow / isqrt / gcd / abs\n"
              << "  to_string / fits_ll / to_ll / to_i128 / to_ld / num_digits / digit_sum / mod_ll / sign / is_zero / is_even、cin >> / cout <<\n";

    api("【1】 BigInt a(\"…\") / BigInt b = 10 / + − × ÷ %", "BigInt", "符号つきの多倍長整数。int・long long と混ぜて書ける（÷ は 0 の方向に切り捨て）", "× O(n^1.58)、÷ O(n m)");
    RUN(BigInt a("123456789012345678901234567890"));
    SHOW(a * a, "");
    SHOW(a * 2 + 1, "");
    SHOW(a / 7, "");
    SHOW(a % 7, "");
    SHOW(BigInt(-7) / 2, "C++ と同じく 0 の方向へ");
    SHOW(BigInt(-7) % 2, "余りの符号は割られる数と同じ");
    SHOW(divmod(a, BigInt(1000000007)).second, "{商, 余り}");
    SHOW(a > BigInt(LLONG_MAX), "比較");

    api("【2】 pow(a, e) / isqrt(a) / gcd(a, b) / lcm(a, b) / abs(a)", "BigInt", "累乗 / ⌊√a⌋ / 最大公約数・最小公倍数（0 以上）/ 絶対値", "pow O(M(n) log e)");
    SHOW(pow(BigInt(2), 100), "2^100");
    SHOW(isqrt(BigInt("1000000000000000000000000000000")), "⌊√10^30⌋ = 10^15");
    SHOW(isqrt(pow(BigInt(10), 41)), "√10^41 = 316227766016837933199.88…");
    SHOW(gcd(pow(BigInt(6), 30), pow(BigInt(10), 20)), "2^20");
    {
        BigInt l = 1;
        for (int i = 1; i <= 100; ++i) l = lcm(l, i);
        std::cout << "  lcm(1, 2, …, 100)（BigInt で順に）-> " << l << "\n";
    }
    SHOW(abs(BigInt(-5)), "");

    api("【3】 to_string / fits_ll / to_ll / to_i128 / to_ld", "string / bool / long long / Int128 / long double", "変換（to_ll・to_i128 は収まらなければ止める）", "O(n)");
    SHOW(BigInt(LLONG_MAX).fits_ll(), "");
    SHOW((BigInt(LLONG_MAX) + 1).fits_ll(), "2^63 は long long に入らない");
    SHOW((BigInt(LLONG_MAX) + 1).to_i128(), "Int128 には入る");
    SHOW(pow(BigInt(3), 100).to_ld(), "おおよその値");

    api("【4】 num_digits / digit_sum / mod_ll(m)", "int / long long / long long", "10 進の桁数 / 各桁の和 / a mod m（[0, m)）", "O(n)");
    SHOW(pow(BigInt(2), 1000).num_digits(), "2^1000 は 302 桁");
    SHOW(pow(BigInt(2), 1000).digit_sum(), "各桁の和（Project Euler 16）");
    {
        BigInt f = 1;
        for (int i = 1; i <= 100; ++i) f *= i;
        note("100! = " + f.to_string().substr(0, 30) + "…（" + std::to_string(f.num_digits()) + " 桁、各桁の和 " + std::to_string(f.digit_sum()) + "）");
    }
    SHOW(BigInt(std::string(1000, '9')).mod_ll(998244353), "1000 桁の数 mod 998244353（文字列で来る巨大な N の mod）");
    note("比べるだけなら桁数 → 文字列の比較でもよい。答えが Int128（±1.7×10^38）に収まるなら Int128 の方が速い");
}

}  // namespace big_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    big_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    big_test::benchmark();    // 速度計測
    big_demo::demo();         // 使い方デモ
    return 0;
}
