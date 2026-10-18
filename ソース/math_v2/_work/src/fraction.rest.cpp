// =============================================================================
//  ランダムテスト
// =============================================================================
namespace frac_test {

{COMMON_TEST}

static const long long LMAX = std::numeric_limits<long long>::max();
static Int128 g128(Int128 a, Int128 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        const Int128 t = a % b;
        a = b, b = t;
    }
    return a;
}
static long long nz(std::mt19937_64& rng, long long lo, long long hi) {  // 0 でない乱数
    long long v;
    do v = rnd(rng, lo, hi);
    while (v == 0);
    return v;
}
template <class F>
static bool same_value(const F& x, Int128 n, Int128 d) {  // x == n/d（d ≠ 0）か（既約・分母 > 0 も確かめる）
    if (x.den <= 0 || g128(x.num, x.den) != 1) return false;
    return (Int128)x.num * d == n * (Int128)x.den;
}

static void test_basic(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 50000 * TEST_SCALE; ++it, ++cases) {
        const long long R = it % 3 == 0 ? 20 : (it % 3 == 1 ? 1000000000 : LMAX);
        const long long a = rnd(rng, -R, R), b = nz(rng, -R, R);
        const Fraction x(a, b);
        FRAC_CHECK(same_value(x, a, b));
        // 四則演算（値が long long に収まる大きさで）
        const long long S = 1000000000;
        const long long p = rnd(rng, -S, S), q = nz(rng, -S, S), r = rnd(rng, -S, S), s = nz(rng, -S, S);
        const Fraction u(p, q), v(r, s);
        FRAC_CHECK(same_value(u + v, (Int128)p * s + (Int128)r * q, (Int128)q * s));
        FRAC_CHECK(same_value(u - v, (Int128)p * s - (Int128)r * q, (Int128)q * s));
        FRAC_CHECK(same_value(u * v, (Int128)p * r, (Int128)q * s));
        if (r != 0) FRAC_CHECK(same_value(u / v, (Int128)p * s, (Int128)q * r));
        FRAC_CHECK(same_value(-u, -p, q) && same_value(u * 3 + 1, (Int128)3 * p + q, q) && same_value(2 - u, (Int128)2 * q - p, q));
        // 比較（大きい値でも）
        const long long c = rnd(rng, -R, R), d = nz(rng, -R, R);
        const Fraction y(c, d);
        const Int128 lhs = (Int128)a * d * (b < 0 ? -1 : 1) * (d < 0 ? -1 : 1), rhs = (Int128)c * b * (b < 0 ? -1 : 1) * (d < 0 ? -1 : 1);
        const int e = lhs < rhs ? -1 : (lhs > rhs ? 1 : 0);
        FRAC_CHECK(compare(x, y) == e && (x < y) == (e < 0) && (x == y) == (e == 0) && (x >= y) == (e >= 0));
        // 切り捨て・切り上げ・四捨五入
        FRAC_CHECK((Int128)x.floor() == math_internal::floor_div128(a, b) && (Int128)x.ceil() == math_internal::ceil_div128(a, b));
        const Int128 fl = math_internal::floor_div128(a, b), two_rem = 2 * ((Int128)x.num - fl * x.den);  // 2 × (x − ⌊x⌋) × den
        Int128 rd = two_rem < x.den ? fl : (two_rem > x.den ? fl + 1 : (x.num >= 0 ? fl + 1 : fl));
        FRAC_CHECK((Int128)x.round() == rd);
        FRAC_CHECK(x.sign() == (a == 0 ? 0 : ((a < 0) != (b < 0) ? -1 : 1)) && x.abs() >= Fraction(0) && (x.is_integer() == (a % b == 0)));
        if (a != 0) FRAC_CHECK(x * x.reciprocal() == Fraction(1));
        FRAC_CHECK(std::fabs((double)x.to_real() - (double)a / (double)b) <= 1e-9 * (1 + std::fabs((double)a / (double)b)));
    }
}

static void test_text(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long a = rnd(rng, -LMAX, LMAX), b = nz(rng, 1, it % 2 ? 1000 : LMAX);
        const Fraction x(a, b);
        // 出して読み戻す
        std::ostringstream os;
        os << x;
        std::istringstream is(os.str());
        Fraction y;
        is >> y;
        FRAC_CHECK(x == y && (x.den != 1 || os.str().find('/') == std::string::npos));
        // 小数 k 桁（切り捨て）：愚直な割り算
        const int k = (int)rnd(rng, 0, 25);
        const Int128 A = a < 0 ? -(Int128)a : a;
        Int128 rem = A % b;
        std::string digits;
        for (int i = 0; i < k; ++i) rem *= 10, digits += (char)('0' + (int)(rem / b)), rem %= b;
        const std::string body = int128_to_string(A / b) + (k ? "." + digits : "");
        const bool zero = A / b == 0 && digits.find_first_not_of('0') == std::string::npos;
        const std::string e = (a < 0 && !zero ? "-" : "") + body;
        FRAC_CHECK(Fraction(a, b).to_decimal(k) == e);
        // 小数を読む：to_decimal の文字列はその値（切り捨て）を表す
        Int128 p10 = 1;
        for (int i = 0; i < k; ++i) p10 *= 10;
        if (k > 0 && k <= 18 && (A / b + 1) * p10 < (Int128)LMAX) {  // 読んだ値が long long に収まるとき
            const Fraction t = Fraction::parse(e);
            FRAC_CHECK((Int128)t.num * p10 == (parse_int128(int128_to_string(A / b) + digits) * (a < 0 && !zero ? -1 : 1)) * (Int128)t.den);
        }
    }
    FRAC_CHECK(Fraction::parse("-1.25") == Fraction(-5, 4) && Fraction::parse("3/-6") == Fraction(-1, 2) && Fraction::parse(".5") == Fraction(1, 2) && Fraction::parse("7") == Fraction(7));
    FRAC_CHECK(Fraction(-1, 3).to_decimal(5) == "-0.33333" && Fraction(-1, 3000000).to_decimal(5) == "0.00000" && Fraction(7, 2).to_decimal(0) == "3");
    // ハッシュ
    std::unordered_map<Fraction, int> mp;
    for (int i = 1; i <= 100; ++i) mp[Fraction(i, 2 * i)]++, mp[Fraction(-i, 3 * i)]++;
    FRAC_CHECK(mp.size() == 2 && mp[Fraction(1, 2)] == 100);
}

static void test_128(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 15000 * TEST_SCALE; ++it, ++cases) {
        const long long R = 1LL << 60;
        const long long a = rnd(rng, -R, R), b = nz(rng, -R, R), c = rnd(rng, -R, R), d = nz(rng, -R, R);
        const Fraction128 x(a, b), y(c, d);
        FRAC_CHECK(same_value(x, a, b) && x + y == Fraction128((Int128)a * d + (Int128)c * b, (Int128)b * d) && x * y == Fraction128((Int128)a * c, (Int128)b * d));
        FRAC_CHECK(x - y == Fraction128((Int128)a * d - (Int128)c * b, (Int128)b * d) && (c == 0 || x / y == Fraction128((Int128)a * d, (Int128)b * c)));
        const Int128 lhs = (Int128)a * d * (b < 0 ? -1 : 1) * (d < 0 ? -1 : 1), rhs = (Int128)c * b * (b < 0 ? -1 : 1) * (d < 0 ? -1 : 1);
        const int e = lhs < rhs ? -1 : (lhs > rhs ? 1 : 0);
        FRAC_CHECK(compare(x, y) == e && compare(y, x) == -e);
        // 大きい分子・分母：比較を long double と（十分離れているとき）
        const Int128 P60 = (Int128)1 << 60;
        const Fraction128 big1 = Fraction128((Int128)rnd(rng, -LMAX, LMAX) * P60 + rnd_ll(rng), (Int128)rnd(rng, 1, LMAX) * P60 + 1);
        const Fraction128 big2 = Fraction128((Int128)rnd(rng, -LMAX, LMAX) * P60 + rnd_ll(rng), (Int128)rnd(rng, 1, LMAX) * P60 + 3);
        const long double r1 = big1.to_real(), r2 = big2.to_real();
        if (std::fabs((double)(r1 - r2)) > 1e-6 * (std::fabs((double)r1) + std::fabs((double)r2) + 1)) FRAC_CHECK((big1 < big2) == (r1 < r2));
        FRAC_CHECK(big1.to_decimal(30).size() > 30 && Fraction128::parse(big1.to_decimal(0)) == Fraction128(big1.num / big1.den));
        FRAC_CHECK(compare(big1, big1) == 0 && compare(big1 + Fraction128(1, 3), big1) == 1 && compare(big1, big1 + Fraction128(1, 3)) == -1);
    }
}

static void test_edge128(int& cases) {  // 途中の値が −2^127 になる足し算、分母が 2^125 を超える小数表示
    const Int128 P126 = (Int128)1 << 126;
    const Fraction128 x(-P126 - 1, 2), y(-P126 + 1, 2), z(P126 - 1, 2);
    FRAC_CHECK(x + y == Fraction128(-P126) && x - z == Fraction128(-P126));
    const Int128 d = (Int128)3 << 125;  // (d − 1)/d = 0.999…
    FRAC_CHECK(Fraction128(d - 1, d).to_decimal(3) == "0.999" && Fraction128(-(d - 1), d).to_decimal(5) == "-0.99999" && Fraction128(1, d).to_decimal(3) == "0.000");
    cases += 3;
}

static void test_approx(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 15000 * TEST_SCALE; ++it, ++cases) {
        const long long q = nz(rng, -1000000, 1000000), p = it % 2 ? rnd(rng, -1000000, 1000000) : rnd(rng, -50, 50);
        const long long M = rnd(rng, 1, 60);
        // 愚直：分母 k = 1〜M、分子は p k / q の前後
        Int128 bh = 0, bk = 0;
        for (long long k = 1; k <= M; ++k) {
            const Int128 f = math_internal::floor_div128((Int128)p * k, q);
            for (Int128 h = f; h <= f + 1; ++h) {
                // |p/q − h/k| = |p k − h q| / (|q| k)。前の候補 bh/bk と比べる
                const Int128 dn = (Int128)p * k - h * q, dist = dn < 0 ? -dn : dn;
                if (bk == 0) {
                    bh = h, bk = k;
                    continue;
                }
                const Int128 bdn = (Int128)p * bk - bh * q, bdist = bdn < 0 ? -bdn : bdn;
                if (dist * bk < bdist * k) bh = h, bk = k;
            }
        }
        const Fraction got = best_rational_approx(p, q, M);
        FRAC_CHECK(got == Fraction(bh, bk));
    }
    FRAC_CHECK(best_rational_approx(314159265358979LL, 100000000000000LL, 1000) == Fraction(355, 113) && best_rational_approx(3, 4, 2) == Fraction(1));
    // 分母が max_den 以下で p/q 以下の最大・以上の最小
    for (int it = 0; it < 15000 * TEST_SCALE; ++it, ++cases) {
        const long long q = nz(rng, -100000, 100000), p = it % 2 ? rnd(rng, -1000000, 1000000) : rnd(rng, -50, 50);
        const long long M = rnd(rng, 1, 60);
        const Fraction x(p, q);
        Fraction lo, hi;
        bool has = false;
        for (long long k = 1; k <= M; ++k) {
            const Fraction f(math_internal::floor_div128((Int128)p * k, q), k), c(math_internal::ceil_div128((Int128)p * k, q), k);
            if (!has || lo < f) lo = f;
            if (!has || c < hi) hi = c;
            has = true;
        }
        FRAC_CHECK(rational_floor(p, q, M) == lo && rational_ceil(p, q, M) == hi && lo <= x && x <= hi);
    }
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {  // 大きい値：隣り合う（右 − 左 = 1/(分母の積)）こと、分母の和が max_den を超えること
        const long long q = nz(rng, -LLONG_MAX, LLONG_MAX), p = rnd(rng, -LLONG_MAX, LLONG_MAX), M = rnd(rng, 1, it % 2 ? 1000000 : LLONG_MAX);  // Fraction の分子は −2^63 を使わない
        const Fraction l = rational_floor(p, q, M), r = rational_ceil(p, q, M);
        const Fraction128 x(p, q);
        FRAC_CHECK(l.den <= M && r.den <= M && Fraction128(l.num, l.den) <= x && x <= Fraction128(r.num, r.den));
        if (!(l == r)) FRAC_CHECK((Int128)r.num * l.den - (Int128)l.num * r.den == 1 && (Int128)l.den + r.den > M);
    }
    // いちばん簡単な分数（開区間）
    for (int it = 0; it < 15000 * TEST_SCALE; ++it, ++cases) {
        Fraction x(rnd(rng, -30, 30), rnd(rng, 1, 30)), y(rnd(rng, -30, 30), rnd(rng, 1, 30));
        if (x == y) continue;
        if (y < x) std::swap(x, y);
        Fraction best;
        bool found = false;
        for (long long k = 1; k <= 61 && !found; ++k) {  // 分母の小さい順、分子は絶対値の小さい順
            for (long long a = 0; a <= 2000 && !found; ++a)
                for (long long s : {1, -1}) {
                    const Fraction f(s * a, k);
                    if (x < f && f < y) {
                        best = f, found = true;
                        break;
                    }
                }
        }
        FRAC_CHECK(found && simplest_between(x, y) == best);
    }
    FRAC_CHECK(simplest_between(Fraction(1, 3), Fraction(1, 2)) == Fraction(2, 5) && simplest_between(Fraction(-7, 2), Fraction(-3)) == Fraction(-10, 3));
    FRAC_CHECK(simplest_between(Fraction(999999999999999999LL, 1000000000000000000LL), Fraction(1)).den == 1000000000000000001LL);
    // 連分数
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long q = nz(rng, -LLONG_MAX, LLONG_MAX), p = rnd(rng, -LLONG_MAX, LLONG_MAX);
        const std::vector<long long> cf = continued_fraction(p, q);
        for (size_t i = 1; i < cf.size(); ++i) FRAC_CHECK(cf[i] >= 1 && (i + 1 < cf.size() || cf[i] >= 2));
        FRAC_CHECK(from_continued_fraction(cf) == Fraction(p, q));
    }
    FRAC_CHECK(continued_fraction(355, 113) == std::vector<long long>({3, 7, 16}) && continued_fraction(-7, 3) == std::vector<long long>({-3, 1, 2}) && from_continued_fraction({1, 2, 2, 2}) == Fraction(17, 12));
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_basic(rng, cases);
    test_text(rng, cases);
    test_128(rng, cases);
    test_edge128(cases);
    test_approx(rng, cases);
    std::cout << "[fraction] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[fraction] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 1000000;
    std::vector<Fraction> v(n);
    for (auto& x : v) x = Fraction(rnd(rng, 0, 1000000000), rnd(rng, 1, 1000000000));
    auto t0 = std::chrono::steady_clock::now();
    std::sort(v.begin(), v.end());
    const long long t_sort = ms(t0);
    t0 = std::chrono::steady_clock::now();
    Fraction s = 0;
    for (int i = 0; i < n; ++i) s = Fraction(v[(size_t)i].num % 1000 + 1, v[(size_t)i].den % 1000 + 1) * Fraction(i % 7 + 1, i % 5 + 1) + Fraction(1, i % 9 + 1);
    const long long t_ops = ms(t0);
    std::vector<Fraction128> w(n / 5);
    for (auto& x : w) x = Fraction128(((Int128)rnd(rng, 0, LMAX) << 40) + 1, ((Int128)rnd(rng, 1, LMAX) << 40) + 1);
    t0 = std::chrono::steady_clock::now();
    std::sort(w.begin(), w.end());
    const long long t_sort128 = ms(t0);
    std::cout << "[fraction] benchmark: Fraction 10^6 個の sort " << t_sort << " ms, 四則演算 × 10^6 " << t_ops << " ms, Fraction128（分子・分母 ≈ 2^100）2×10^5 個の sort "
              << t_sort128 << " ms  (check " << (long long)(s.num % 7 + w[5].num % 7) << ")\n";
}

}  // namespace frac_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace frac_demo {

{COMMON_DEMO}

static void demo() {
    banner("分数（fraction） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  Fraction x(n, d) / x.num / x.den / + − × ÷ / == < など / compare(x, y)\n"
              << "  x.floor() / x.ceil() / x.round() / x.sign() / x.abs() / x.reciprocal() / x.is_integer() / x.to_real() / x.to_decimal(k)\n"
              << "  cout << x / cin >> x / Fraction::parse(s) / unordered_map<Fraction, V>\n"
              << "  best_rational_approx(p, q, max_den)      分母が小さい近い分数\n"
              << "  Fraction128                              分子・分母が Int128\n";

    api("【1】 Fraction x(n, d) / + − × ÷ / 比較", "Fraction",
        "いつも既約で分母 > 0。途中は 128 ビットで計算し、結果が long long に収まらなければ止める", "O(log)");
    RUN(Fraction a(1, 3), b(-2, 4));
    SHOW(b, "約分され、符号は分子に");
    SHOW(a + b, "1/3 − 1/2");
    SHOW(a * b, "掛け算");
    SHOW(a / b, "割り算");
    SHOW(a * 3, "整数と混ぜても");
    SHOW(a < b, "誤差なしの比較");
    SHOW(compare(Fraction(2, 6), a), "等しいので 0");
    SHOW(a.num, "分子");
    SHOW(a.den, "分母");
    {
        // ABC308 C：A/(A+B) の大きい順、同じなら番号の小さい順（stable_sort + 誤差なしの比較）
        auto solve = [](const std::vector<std::pair<long long, long long>>& ab) {
            std::vector<int> id(ab.size());
            std::iota(id.begin(), id.end(), 0);
            std::stable_sort(id.begin(), id.end(), [&](int i, int j) {
                return Fraction(ab[(size_t)i].first, ab[(size_t)i].first + ab[(size_t)i].second) > Fraction(ab[(size_t)j].first, ab[(size_t)j].first + ab[(size_t)j].second);
            });
            std::string s;
            for (int i : id) s += (s.empty() ? "" : " ") + std::to_string(i + 1);
            return s;
        };
        note("ABC308 C（Standings）: 成功率 A/(A+B) を Fraction で比べて stable_sort");
        std::cout << "    例 1 → " << solve({{1, 3}, {3, 1}, {2, 2}}) << "、例 2 → " << solve({{1, 3}, {2, 6}})
                  << "、例 3 → " << solve({{999999999, 1000000000}, {333333333, 999999999}, {1000000000, 999999997}, {999999998, 1000000000}}) << "\n";
    }

    api("【2】 x.floor() / x.ceil() / x.round() / x.sign() / x.abs() / x.reciprocal() / x.to_real() / x.to_decimal(k)", "long long / Fraction / long double / string",
        "切り捨て・切り上げ・四捨五入（.5 は 0 から遠い方）/ 符号 / 絶対値 / 逆数 / 実数 / 小数第 k 位まで（切り捨て）", "O(1) / O(k)");
    RUN(Fraction c(-7, 2));
    SHOW(c.floor(), "⌊−3.5⌋");
    SHOW(c.ceil(), "⌈−3.5⌉");
    SHOW(c.round(), "−3.5 の四捨五入（0 から遠い方）");
    SHOW(c.abs(), "絶対値");
    SHOW(c.reciprocal(), "逆数");
    SHOW(Fraction(1, 7).to_decimal(20), "1/7 を小数第 20 位まで（誤差なし）");
    note("出力の桁を決めて実数で出すなら cout << fixed << setprecision(15) << (double)x.to_real()");

    api("【3】 cout << x / cin >> x / Fraction::parse(s) / unordered_map<Fraction, V>", "", "\"p/q\"（整数なら \"p\"）で出す / \"p/q\"・\"p\"・\"-1.25\" を読む / ハッシュ", "O(桁数)");
    SHOW(Fraction::parse("-1.25"), "小数を読む");
    SHOW(Fraction::parse("6/-4"), "分数を読む");
    RUN(std::unordered_map<Fraction, int> cnt; cnt[Fraction(1, 2)]++; cnt[Fraction(2, 4)]++);
    SHOW(cnt[Fraction(1, 2)], "同じ値は同じキー（傾きを数えるときなど）");

    api("【4】 best_rational_approx(p, q, max_den)", "Fraction", "p/q にいちばん近い、分母が max_den 以下の分数（同じ近さなら分母の小さい方）", "O(log)");
    SHOW(best_rational_approx(314159265358979LL, 100000000000000LL, 1000), "π ≈ 355/113");
    SHOW(best_rational_approx(314159265358979LL, 100000000000000LL, 100), "分母 100 以下なら 311/99");
    SHOW(best_rational_approx(1, 3, 2), "1/3 に近い分母 2 以下：1/2 と 0 は 1/6 と 1/3 離れている");

    api("【5】 Fraction128", "Fraction128", "分子・分母が Int128（10^38 まで）。比較は掛け算をしない方法で誤差も溢れもない", "O(log)");
    RUN(Fraction128 h(1000000000000000000LL, 3));
    SHOW(h * h, "10^36 / 9（long long には入らない）");
    SHOW(h * h > Fraction128((Int128)111111111111111111LL * 1000000000000000000LL), "比較");

    api("【6】 rational_floor(p, q, N) / rational_ceil(p, q, N) / simplest_between(x, y) / continued_fraction(p, q) / from_continued_fraction(cf)",
        "Fraction / vector<long long>", "分母 N 以下で p/q 以下の最大・以上の最小 / x < r < y で分母最小の r / 連分数展開とその逆", "O(log)");
    SHOW(rational_floor(314159265358979LL, 100000000000000LL, 100), "π 以下で分母 100 以下の最大");
    SHOW(rational_ceil(314159265358979LL, 100000000000000LL, 100), "π 以上で分母 100 以下の最小");
    SHOW(simplest_between(Fraction(1, 3), Fraction(1, 2)), "1/3 と 1/2 のあいだで分母最小");
    SHOW(simplest_between(Fraction(-7, 2), Fraction(-3)), "負でも");
    SHOW(continued_fraction(355, 113), "355/113 = 3 + 1/(7 + 1/16)");
    SHOW(from_continued_fraction({1, 2, 2, 2}), "√2 の近似 17/12");
    note("区間 (x, y) に入る分数の分母の最小値は simplest_between(x, y).den（それより小さい分母の分数は入らない）");
}

}  // namespace frac_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    frac_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    frac_test::benchmark();    // 速度計測
    frac_demo::demo();         // 使い方デモ
    return 0;
}
