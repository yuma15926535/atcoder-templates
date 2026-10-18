// =============================================================================
//  ランダムテスト
// =============================================================================
namespace dig_test {

{COMMON_TEST}

static const long long LMAX = std::numeric_limits<long long>::max();
static std::string brute_base(unsigned long long u, unsigned long long b) {  // 正の進数（桁は 0-9A-Z。36 まで）
    static const char* D = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    if (u == 0) return "0";
    std::string s;
    while (u) s = D[u % b] + s, u /= b;
    return s;
}
static bool brute_pal(const std::string& s) { return std::string(s.rbegin(), s.rend()) == s; }
static Int128 parse128(const std::string& s) {  // 10 進（先頭の - / + 可）
    size_t i = (s[0] == '-' || s[0] == '+') ? 1 : 0;
    Int128 v = 0;
    for (; i < s.size(); ++i) v = v * 10 + (s[i] - '0');
    return s[0] == '-' ? -v : v;
}

static void test_digits_and_base(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd_ll(rng);
        const long long base = it % 3 == 0 ? 10 : (it % 3 == 1 ? rnd(rng, 2, 36) : rnd(rng, 2, LMAX));
        const unsigned long long u = n < 0 ? 0ULL - (unsigned long long)n : (unsigned long long)n;
        std::vector<int> ds;
        unsigned long long t = u, s = 0;
        do ds.insert(ds.begin(), (int)(t % (unsigned long long)base)), s += t % (unsigned long long)base, t /= (unsigned long long)base;
        while (t);
        DIG_CHECK(digit_sum(n, base) == (long long)s && num_digits(n, base) == (int)ds.size());
        if (base <= 2147483647LL) {
            DIG_CHECK(digits_of(n, base) == ds);
            if (u <= (unsigned long long)LMAX) DIG_CHECK(from_digits(ds, base) == (long long)u);
        }
        if (base <= 36) {
            const std::string b = brute_base(u, (unsigned long long)base);
            DIG_CHECK(to_base(n, (int)base) == (n < 0 ? "-" : "") + b);
            DIG_CHECK(from_base(to_base(n, (int)base), (int)base) == n);
            std::string low = b;
            for (char& c : low) c = (char)std::tolower(c);
            if (n >= 0) DIG_CHECK(from_base(low, (int)base) == n);
        }
        // 負の進数：値・桁の範囲・先頭が 0 でないこと
        const int nb = -(int)rnd(rng, 2, 36);
        const std::string ns = to_base(n, nb);
        Int128 v = 0;
        for (char c : ns) {
            const int d = c <= '9' ? c - '0' : c - 'A' + 10;
            DIG_CHECK(0 <= d && d < -nb);
            v = v * nb + d;
        }
        DIG_CHECK(v == n && (ns == "0" || ns[0] != '0') && from_base(ns, nb) == n);
    }
    DIG_CHECK(to_base(-9, -2) == "1011" && to_base(0, -2) == "0" && to_base(255, 16) == "FF" && to_base(-255, 16) == "-FF");
    DIG_CHECK(to_base(std::numeric_limits<long long>::min(), 2) == "-1" + std::string(63, '0') && from_base("-1" + std::string(63, '0'), 2) == std::numeric_limits<long long>::min());
    DIG_CHECK(digits_of(0) == std::vector<int>{0} && num_digits(0) == 1 && digit_sum(-123) == 6);
    // from_base_cap / base_value_leq
    for (int it = 0; it < 100000 * TEST_SCALE; ++it, ++cases) {
        const int len = (int)rnd(rng, 1, 25);
        const long long base = it % 2 ? rnd(rng, 2, 36) : rnd(rng, 10, it % 4 == 0 ? LMAX : 100000);
        std::string s;
        for (int i = 0; i < len; ++i) s += "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[rng() % std::min<long long>(base, it % 3 ? 10 : 36)];
        const Int128 SAT = (Int128)1 << 100;
        Int128 v = 0;
        for (char c : s) {
            const int d = c <= '9' ? c - '0' : c - 'A' + 10;
            v = v > SAT / base ? SAT : std::min<Int128>(v * base + d, SAT);  // 愚直解も 2^100 で止める
        }
        long long cap = rnd_ll(rng);
        cap = cap < 0 ? (cap == std::numeric_limits<long long>::min() ? 0 : -cap) : cap;
        DIG_CHECK((Int128)from_base_cap(s, base, cap) == std::min<Int128>(v, cap));
        long long lim = it % 2 ? rnd_ll(rng) : (v <= LMAX ? (long long)v + rnd(rng, -1, 1) : LMAX);
        DIG_CHECK(base_value_leq(s, base, lim) == (v <= lim));
    }
}

static void test_palindrome(std::mt19937_64& rng, int& cases) {
    const long long M = 2000000;
    std::vector<long long> pals;
    for (long long x = 0; x <= M; ++x)
        if (brute_pal(std::to_string(x))) pals.push_back(x);
    long long cnt = 0;  // [0, x] の回文数の数
    size_t k = 0;
    for (long long x = 0; x <= M; ++x, ++cases) {
        const bool p = k < pals.size() && pals[k] == x;
        cnt += p;
        DIG_CHECK(is_palindrome_num(x) == p && count_palindromes_up_to(x) == cnt && count_palindromes_below(x) == cnt - p);
        if (x < 1900000) DIG_CHECK(ceil_palindrome(x) == pals[k]);  // k はまだ x 以上で最初の回文数の番号
        if (p) {
            DIG_CHECK(palindrome_rank(x) == (long long)k && kth_palindrome((long long)k) == std::to_string(x));
            ++k;
        }
        DIG_CHECK(floor_palindrome(x) == pals[k - 1]);  // k − 1 は x 以下で最後の回文数の番号
        if (x < 3000) DIG_CHECK(is_palindrome_num(x, 2) == brute_pal(brute_base((unsigned long long)x, 2)));
    }
    // 大きい数：kth と rank、ceil / floor のあいだに回文数がないこと
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long kk = rnd(rng, 0, it % 2 ? 1000000000LL : 1900000000LL);
        const std::string s = kth_palindrome(kk);
        DIG_CHECK(brute_pal(s) && (s.size() == 1 || s[0] != '0'));
        DIG_CHECK(palindrome_rank(std::stoll(s)) == kk);
        const long long n = rnd(rng, 0, it % 3 ? 1000000000000000LL : 9000000000000000000LL);
        const long long c = ceil_palindrome(n), f = floor_palindrome(n);
        DIG_CHECK(c >= n && f <= n && is_palindrome_num(c) && is_palindrome_num(f));
        DIG_CHECK(count_palindromes_below(c) == count_palindromes_below(n) && count_palindromes_up_to(f) == count_palindromes_up_to(n));
        const long long half = rnd(rng, 1, 999999999);
        std::string hs = std::to_string(half);
        DIG_CHECK(std::to_string(make_palindrome(half, true)) == hs + std::string(hs.rbegin() + 1, hs.rend()) && std::to_string(make_palindrome(half, false)) == hs + std::string(hs.rbegin(), hs.rend()));
    }
    // 長い回文数（文字列）：桁数ごとの個数と合うこと
    long long base = 0;
    for (int len = 1; len <= 36; ++len) {
        const long long c = count_palindromes_of_length(len);
        DIG_CHECK(kth_palindrome(base).size() == (size_t)len && kth_palindrome(base + c - 1) == std::string((size_t)len, '9'));
        if (__builtin_add_overflow(base, c, &base)) break;
        ++cases;
    }
    DIG_CHECK(kth_palindrome(45) == "363" && ceil_palindrome(LMAX - 1000000000000LL) == 9223371037301733229LL && floor_palindrome(LMAX) == 9223372036302733229LL);
    // 文字列の回文数：小さい数は全部、long long の範囲は long long 版と、長い数は性質で
    {
        size_t kk = 0;
        Int128 c = 0;
        for (long long x = 0; x <= 300000; ++x, ++cases) {
            const std::string s = std::to_string(x);
            while (pals[kk] < x) ++kk;
            const bool p = pals[kk] == x;
            c += p;
            DIG_CHECK(ceil_palindrome_str(s) == std::to_string(pals[kk]) && floor_palindrome_str(s) == std::to_string(pals[kk - !p]));
            DIG_CHECK(count_palindromes_up_to_str(s) == c && count_palindromes_below_str(s) == c - p && is_palindrome_str(s) == p);
            if (x % 1000 == 7) DIG_CHECK(ceil_palindrome_str("00" + s) == ceil_palindrome_str(s) && count_palindromes_up_to_str("+0" + s) == c);
        }
    }
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long n = rnd(rng, 0, it % 2 ? 1000000000000LL : 9000000000000000000LL);
        const std::string s = std::to_string(n);
        DIG_CHECK(ceil_palindrome_str(s) == std::to_string(ceil_palindrome(n)) && floor_palindrome_str(s) == std::to_string(floor_palindrome(n)));
        DIG_CHECK(count_palindromes_up_to_str(s) == count_palindromes_up_to(n) && count_palindromes_below_str(s) == count_palindromes_below(n));
        // 長い数（0・1・9 が多いと境目を通る）
        const int len = (int)rnd(rng, 1, 75);
        std::string t(1, char('1' + rng() % 9));
        const char* pool = it % 3 == 0 ? "019" : (it % 3 == 1 ? "0123456789" : "09");
        const int psz = (int)std::char_traits<char>::length(pool);
        while ((int)t.size() < len) t += pool[rng() % psz];
        if (it % 5 == 0) t = make_palindrome_str(t.substr(0, (len + 1) / 2), len % 2) + "";  // 回文そのもの
        if (it % 7 == 0) t = std::string((size_t)len, '9');
        if (it % 11 == 0) t = "1" + std::string((size_t)len - 1, '0');
        const std::string ce = ceil_palindrome_str(t), fl = floor_palindrome_str(t);
        DIG_CHECK(brute_pal(ce) && brute_pal(fl) && str_num_compare(fl, t) <= 0 && str_num_compare(t, ce) <= 0);
        DIG_CHECK((ce.size() == 1 || ce[0] != '0') && (fl.size() == 1 || fl[0] != '0'));
        if (t.size() <= 74) {
            const Int128 up = count_palindromes_up_to_str(t), below = count_palindromes_below_str(t);
            DIG_CHECK(count_palindromes_below_str(ce) == below && count_palindromes_up_to_str(fl) == up);  // あいだに回文数がない
            DIG_CHECK(up - below == (brute_pal(t) ? 1 : 0) && count_palindromes_up_to_str(ce) == below + 1);
        }
    }
    DIG_CHECK(ceil_palindrome_str("100000000000000000000000") == "100000000000000000000001" && floor_palindrome_str("100000000000000000000000") == "99999999999999999999999");
    DIG_CHECK(ceil_palindrome_str("123456789012345678901234567890") == "123456789012346643210987654321" && ceil_palindrome_str("199") == "202");
    DIG_CHECK(count_palindromes_up_to_str(kth_palindrome(1000000000000000000LL)) == (Int128)1000000000000000001LL && make_palindrome_str("123", true) == "12321");
    DIG_CHECK(count_palindromes_up_to_str(std::string(75, '9')) > 0 && floor_palindrome_str("0") == "0" && ceil_palindrome_str("000") == "0");
}

static void test_concat_and_strings(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {
        const long long a = rnd(rng, 0, it % 2 ? 1000 : LMAX), b = rnd(rng, 0, it % 3 ? 1000 : LMAX), m = rnd(rng, 1, it % 2 ? 1000000007 : LMAX);
        const std::string s = std::to_string(a) + std::to_string(b);
        DIG_CHECK(concat_value(a, b) == parse128(s));
        DIG_CHECK(concat_mod(a, b, m) == (long long)(parse128(s) % m) && str_mod(s, m) == (long long)(parse128(s) % m));
        const long long k = rnd(rng, 0, 12), x = rnd(rng, 0, it % 2 ? 99 : 999999);
        long long r = 0, rep = 0;
        const long long p = pow10_mod(num_digits(x), m);
        for (long long i = 0; i < k; ++i) r = (long long)(((Int128)r * p + x) % m), rep = (long long)(((Int128)rep * 10 + 1) % m);
        DIG_CHECK(repeat_concat_mod(x, k, m) == r && repunit_mod(k, m) == rep);
        {
            std::string rs;
            for (long long i = 0; i < k; ++i) rs += std::to_string(x);
            if (rs.size() <= 38) DIG_CHECK(repeat_concat_value(x, k) == (rs.empty() ? 0 : parse128(rs)));
        }
        // base 進数の文字列の mod（36 進 20 桁まで）
        {
            const int bs = (int)rnd(rng, 2, 36), len = (int)rnd(rng, 1, 20);
            std::string t;
            Int128 val = 0;
            for (int i = 0; i < len; ++i) {
                const int d = (int)(rng() % bs);
                t += (d < 10 ? char('0' + d) : char((rng() % 2 ? 'A' : 'a') + d - 10));
                val = val * bs + d;
            }
            const bool neg = rng() % 2;
            DIG_CHECK(str_mod((neg ? "-" : "") + t, m, bs) == (long long)(((neg ? -val : val) % m + m) % m));
            const long long y = rnd_ll(rng);
            DIG_CHECK(str_mod(to_base(y, bs), m, bs) == (long long)(((Int128)y % m + m) % m));
        }
        // 繰り上がり
        {
            const int bs = it % 2 ? 10 : (int)rnd(rng, 2, 36), len = (int)rnd(rng, 0, 8);
            std::vector<long long> co(len);
            Int128 val = 0, pw = 1;
            for (int i = 0; i < len; ++i) co[i] = rnd(rng, it % 3 == 0 ? -1000 : 0, it % 4 == 0 ? 1000000000000LL : 50), val += co[i] * pw, pw *= bs;
            if (val >= 0) {
                std::string expect;
                for (Int128 v = val; v > 0; v /= bs) expect += "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[(int)(v % bs)];
                if (expect.empty()) expect = "0";
                std::reverse(expect.begin(), expect.end());
                DIG_CHECK(carry_to_string(co, bs) == expect);
            }
            // 桁の畳み込み = 掛け算
            const long long A = rnd(rng, 0, LMAX), B = rnd(rng, 0, LMAX);
            std::vector<int> da = digits_of(A), db = digits_of(B);
            std::reverse(da.begin(), da.end()), std::reverse(db.begin(), db.end());
            std::vector<long long> cv(da.size() + db.size() - 1, 0);
            for (size_t i = 0; i < da.size(); ++i)
                for (size_t j = 0; j < db.size(); ++j) cv[i + j] += da[i] * db[j];
            std::ostringstream os;
            os << (Int128)A * B;
            DIG_CHECK(carry_to_string(cv) == os.str());
        }
        // 文字列の比較・mod（36 桁まで。符号・先頭の 0）
        auto rnds = [&]() {
            std::string t = std::to_string(rnd(rng, 0, 9));
            const int len = (int)rnd(rng, 0, 35);
            for (int i = 0; i < len; ++i) t += char('0' + rng() % 10);
            if (rng() % 4 == 0) t = "00" + t;
            if (rng() % 2) t = "-" + t;
            return t;
        };
        const std::string u = rnds(), v = rnds();
        const Int128 iu = parse128(u), iv = parse128(v);
        DIG_CHECK(str_num_compare(u, v) == (iu > iv) - (iu < iv) && str_num_compare(u, u) == 0);
        DIG_CHECK(str_mod(u, m) == (long long)((iu % m + m) % m));
        // 部分文字列
        std::string w;
        for (int i = 0; i < 40; ++i) w += char('0' + rng() % 10);
        const DigitPrefix P(w, m);
        const int l = (int)rnd(rng, 0, 40), rr = (int)rnd(rng, l, std::min(40, l + 18));
        DIG_CHECK(substr_value(w, l, rr) == (rr > l ? std::stoll(w.substr(l, rr - l)) : 0) && P.value(l, rr) == substr_value(w, l, rr) % m);
        const int l2 = (int)rnd(rng, 0, 40), r2 = (int)rnd(rng, l2, 40);
        DIG_CHECK(P.value(l2, r2) == str_mod(r2 > l2 ? w.substr(l2, r2 - l2) : "0", m));
    }
    DIG_CHECK(str_num_compare("-0", "0") == 0 && str_num_compare("007", "7") == 0 && str_num_compare("-10", "-9") == -1);
    DIG_CHECK(repeat_concat_value(9, 38) == parse128(std::string(38, '9')) && repeat_concat_value(0, LMAX) == 0 && repeat_concat_value(123, 0) == 0);
    DIG_CHECK(carry_to_string({12, 3}) == "42" && carry_to_string(std::vector<int>{}) == "0" && carry_to_string({0, 0, 0}) == "0" && carry_to_string({-1, 1}) == "9");
}

static void test_counting(std::mt19937_64& rng, int& cases) {
    Int128 td = 0, dsum_total = 0;
    Int128 cnt[10] = {0};
    for (long long x = 1; x <= 300000; ++x, ++cases) {
        const std::string s = std::to_string(x);
        td += (Int128)s.size();
        for (char c : s) cnt[c - '0']++, dsum_total += c - '0';
        DIG_CHECK(total_digits_up_to(x) == td && digit_sum_total_up_to(x) == dsum_total);
        if (x % 13 == 0 || x < 2000)
            for (int d = 0; d <= 9; ++d) DIG_CHECK(count_digit_up_to(x, d) == cnt[d]);
    }
    for (int it = 0; it < 30000 * TEST_SCALE; ++it, ++cases) {  // 大きい N は 1 つ前との差
        const long long N = rnd(rng, 2, it % 2 ? 1000000000000LL : LMAX);
        const std::string s = std::to_string(N);
        DIG_CHECK(total_digits_up_to(N) - total_digits_up_to(N - 1) == (Int128)s.size());
        long long dsum = 0;
        for (char c : s) dsum += c - '0';
        DIG_CHECK(digit_sum_total_up_to(N) - digit_sum_total_up_to(N - 1) == dsum);
        const int d = (int)rnd(rng, 0, 9);
        DIG_CHECK(count_digit_up_to(N, d) - count_digit_up_to(N - 1, d) == (Int128)std::count(s.begin(), s.end(), char('0' + d)));
        // 桁数ごとの区間
        const long long l = rnd(rng, 0, N), r = rnd(rng, l, LMAX);
        long long x = l;
        for (const DigitLengthBlock& b : digit_length_blocks(l, r)) {
            DIG_CHECK(b.lo == x && b.lo < b.hi && num_digits(b.lo) == b.len && num_digits(b.hi - 1) == b.len);
            x = b.hi;
        }
        DIG_CHECK(x == std::max(l, r) || (l == r && x == l));
    }
    DIG_CHECK(total_digits_up_to(LMAX) - total_digits_up_to(LMAX - 1) == 19);
    // 使える数字が決まった数
    for (int it = 0; it < 200; ++it, ++cases) {
        std::vector<int> ds;
        for (int d = 0; d <= 9; ++d)
            if (rng() % 3 == 0) ds.push_back(d);
        const int L = (int)rnd(rng, 1, 5);
        std::vector<long long> expect;
        for (long long y = 0; y < pow10_ll(L); ++y) {
            bool ok = true;
            for (char c : std::to_string(y)) ok &= std::find(ds.begin(), ds.end(), c - '0') != ds.end();
            if (ok) expect.push_back(y);
        }
        DIG_CHECK(numbers_from_digits(ds, L) == expect);
    }
}

static void random_test() {
    std::mt19937_64 rng(20261010);
    int cases = 0;
    test_digits_and_base(rng, cases);
    test_palindrome(rng, cases);
    test_concat_and_strings(rng, cases);
    test_counting(rng, cases);
    std::cout << "[digits] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[digits] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int n = 1000000;
    std::vector<long long> xs(n);
    for (auto& x : xs) x = rnd(rng, 0, LMAX);
    auto t0 = std::chrono::steady_clock::now();
    long long s1 = 0;
    for (long long x : xs) s1 += digit_sum(x) + num_digits(x);
    const long long t_ds = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s2 = 0;
    for (int i = 0; i < n / 4; ++i) s2 += (long long)to_base(xs[i], -2).size() + from_base(to_base(xs[i], 7), 7) % 3;
    const long long t_base = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s3 = 0;
    for (int i = 0; i < n / 4; ++i) s3 += ceil_palindrome(xs[i] / 2) % 7 + count_palindromes_up_to(xs[i]) % 7 + (long long)kth_palindrome(xs[i] / 10).size();
    const long long t_pal = ms(t0);
    t0 = std::chrono::steady_clock::now();
    long long s4 = 0;
    for (int i = 0; i < n / 4; ++i) s4 += (long long)(count_digit_up_to(xs[i], i % 10) % 7) + repeat_concat_mod(xs[i] % 1000000, xs[i], 998244353);
    const long long t_cnt = ms(t0);
    std::string big(1000000, '7');
    t0 = std::chrono::steady_clock::now();
    const DigitPrefix P(big, 998244353);
    long long s5 = str_mod(big, 1000000007);
    for (int i = 0; i < n; ++i) s5 += P.value((int)(xs[i] % 500000), (int)(500000 + xs[i] % 500000));
    const long long t_str = ms(t0);
    std::cout << "[digits] benchmark: digit_sum + num_digits × 10^6 " << t_ds << " ms, to_base(−2) + to_base / from_base(7) × 2.5×10^5 " << t_base
              << " ms, 回文数 3 種 × 2.5×10^5 " << t_pal << " ms\n"
              << "[digits] benchmark: count_digit_up_to + repeat_concat_mod × 2.5×10^5 " << t_cnt << " ms, 10^6 桁の str_mod + DigitPrefix の前処理と 10^6 回の value " << t_str
              << " ms  (check " << (s1 + s2 + s3 + s4 + s5) % 7 << ")\n";
}

}  // namespace dig_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace dig_demo {

{COMMON_DEMO}

static void demo() {
    banner("桁・進数・回文数（digits） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  digits_of / digit_sum / num_digits / from_digits              桁（base を指定できる）\n"
              << "  to_base / from_base（±2〜36 進）/ from_base_cap / base_value_leq（36 を超える進数も）\n"
              << "  is_palindrome_num / kth_palindrome / palindrome_rank / count_palindromes_up_to / _below / _of_length\n"
              << "  ceil_palindrome / floor_palindrome / make_palindrome           回文数\n"
              << "  ceil_palindrome_str / floor_palindrome_str / count_palindromes_up_to_str / _below_str / make_palindrome_str / is_palindrome_str   大きな回文数\n"
              << "  concat_value / repeat_concat_value / concat_mod / repeat_concat_mod / repunit_mod / pow10_mod   数をつなげる\n"
              << "  str_mod / str_num_compare / substr_value / DigitPrefix / carry_to_string   大きな数の文字列\n"
              << "  digit_length_blocks / total_digits_up_to / count_digit_up_to / digit_sum_total_up_to   桁の数え上げ\n"
              << "  numbers_from_digits                                            使える数字が決まった数の列挙\n";

    api("【1】 digits_of(n, base) / digit_sum(n, base) / num_digits(n, base) / from_digits(ds, base)", "vector<int> / long long / int",
        "各桁（上の桁から）/ 桁和 / 桁数。base は 2 以上（digits_of は 2^31 − 1 まで）。負の数は絶対値で", "O(桁数)");
    SHOW(digits_of(2026), "上の桁から");
    SHOW(digits_of(10, 2), "2 進数の桁");
    SHOW(digit_sum(999), "桁和");
    SHOW(num_digits(1000000), "桁数");
    SHOW(from_digits({1, 0, 1}, 2), "桁の並びから（2 進の 101 = 5）");
    {
        // ABC083 B: 1〜N で桁和が A 以上 B 以下の数の和
        const long long N = 20, A = 2, B = 5;
        long long s = 0;
        for (long long x = 1; x <= N; ++x)
            if (A <= digit_sum(x) && digit_sum(x) <= B) s += x;
        note("ABC083 B（Some Sums）: digit_sum で桁和を調べる");
        std::cout << "    N = 20, A = 2, B = 5 → " << s << "\n";
    }

    api("【2】 to_base(n, base) / from_base(s, base)", "string / long long", "2〜36 進と −2〜−36 進の文字列との変換（10 以上の桁は A〜Z。0・負の数も）", "O(桁数)");
    SHOW(to_base(-9, -2), "ABC105 C: −2 進数");
    SHOW(to_base(255, 16), "16 進数");
    SHOW(from_base("ff", 16), "小文字も読める");
    SHOW(to_base(-10, 3), "正の進数の負の数は - が付く");
    {
        // 典型90 067: 8 進数の N を 9 進数にして 8 を 5 に置き換える、を K 回
        std::string s = "21";
        const int K = 1;
        for (int k = 0; k < K; ++k) {
            s = to_base(from_base(s, 8), 9);
            for (char& c : s)
                if (c == '8') c = '5';
        }
        note("典型90 067（Base 8 to 9）: from_base(s, 8) → to_base(·, 9) → 8 を 5 に");
        std::cout << "    N = 21, K = 1 → " << s << "\n";
    }
    {
        // ABC336 C: 各桁が偶数の数の N 番目（0 が 1 番目）= N − 1 の 5 進数の各桁を 2 倍
        const long long N = 8;
        std::string s = to_base(N - 1, 5);
        for (char& c : s) c = char('0' + (c - '0') * 2);
        note("ABC336 C（Even Digits）: N − 1 を 5 進数にして各桁を 2 倍");
        std::cout << "    N = 8 → " << s << "\n";
    }
    {
        // ABC186 C: N 以下で 10 進でも 8 進でも 7 を含まない数
        const long long N = 20;
        long long c = 0;
        for (long long x = 1; x <= N; ++x) {
            const std::string a = to_base(x, 10), b = to_base(x, 8);
            c += a.find('7') == std::string::npos && b.find('7') == std::string::npos;
        }
        note("ABC186 C（Unlucky 7）: to_base(x, 8) で 8 進数の文字列にして 7 を探す");
        std::cout << "    N = 20 → " << c << "\n";
    }

    api("【3】 from_base_cap(s, base, cap) / base_value_leq(s, base, limit)", "long long / bool",
        "base は 2 以上いくつでも（10^18 でも）。値が cap・limit を超えるかを溢れずに調べる", "O(|s|)");
    SHOW(from_base_cap("99", 1000000000000LL, 1000000000000000000LL), "99 × 10^12 + 9 は上限 10^18 で止まる");
    SHOW(base_value_leq("22", 4, 10), "22 (4 進) = 10 ≤ 10");
    {
        // ABC192 D: X を n 進数（n > X の最大の桁）として読んだ値が M 以下になる n の個数（二分探索）
        const std::string X = "22";
        const long long M = 10;
        long long dmax = 0;
        for (char c : X) dmax = std::max<long long>(dmax, c - '0');
        long long ans;
        if (X.size() == 1) {
            ans = (X[0] - '0') <= M ? 1 : 0;
        } else {
            long long lo = dmax, hi = M + 1;  // lo は条件を満たす最大の進数（dmax は「満たさない」の番兵）
            while (hi - lo > 1) {
                const long long mid = lo + (hi - lo) / 2;
                (base_value_leq(X, mid, M) ? lo : hi) = mid;
            }
            ans = lo - dmax;
        }
        note("ABC192 D（Base n）: base_value_leq(X, n, M) で二分探索（X が 1 桁なら別扱い）");
        std::cout << "    X = 22, M = 10 → " << ans << "（3 進と 4 進）\n";
    }

    api("【4】 kth_palindrome(k) / palindrome_rank(n) / count_palindromes_up_to(n) / ceil_palindrome(n) / floor_palindrome(n)", "string / long long",
        "回文数（10 進）。k 番目は 0 から数える（0 番目が 0）。kth_palindrome は 10^18 を超えるので文字列", "O(桁数)");
    SHOW(kth_palindrome(45), "ABC363 D: 46 番目（N = 46 なら kth_palindrome(N − 1)）");
    SHOW(kth_palindrome(1000000000000000000LL), "とても大きい k");
    SHOW(palindrome_rank(363), "363 は 45 番目");
    SHOW(count_palindromes_up_to(1000), "0〜1000 の回文数");
    SHOW(ceil_palindrome(1234), "1234 以上で最小");
    SHOW(floor_palindrome(1000), "1000 以下で最大");
    SHOW(make_palindrome(123, true), "左半分から（真ん中を共有）");
    SHOW(is_palindrome_num(9, 2), "2 進数で 1001");
    note("long long を超える数は文字列版（_str）。count_palindromes_up_to_str は Int128 を返す（75 桁まで）");
    SHOW(ceil_palindrome_str("123456789012345678901234567890"), "30 桁");
    SHOW(floor_palindrome_str("1000000000000000000000"), "10^21 以下で最大");
    SHOW(count_palindromes_up_to_str("100000000000000000000000000000"), "0〜10^29 の回文数");
    SHOW(make_palindrome_str("abc", false), "文字列でも作れる");
    SHOW(is_palindrome_str("racecar"), "文字列が回文か");

    api("【5】 concat_value(a, b) / concat_mod(a, b, m) / repeat_concat_mod(x, k, m) / repunit_mod(k, m)", "Int128 / long long",
        "数をつなげた値 / その mod / x を k 回つなげた数の mod（k ≤ 10^18）/ 11…1 の mod", "O(1) / O(log k)");
    SHOW(concat_value(123, 45), "12345");
    SHOW(repeat_concat_value(2026, 3), "202620262026");
    SHOW(concat_mod(1000000000LL, 1000000000LL, 998244353), "10^9 の後ろに 10^9（19 桁）を mod で");
    SHOW(repeat_concat_mod(5, 5, 998244353), "ABC357 D: 5 を 5 回 → 55555");
    SHOW(repeat_concat_mod(10, 3, 1000000007), "101010");
    SHOW(repunit_mod(1000000000000000000LL, 998244353), "1 を 10^18 個並べた数 mod p");

    api("【6】 str_mod(s, m, base) / str_num_compare(a, b) / substr_value(s, l, r) / DigitPrefix P(s, m) / carry_to_string(c, base)", "long long / int / string",
        "long long に入らない数（文字列）の mod・大小 / 部分文字列の値（P.value(l, r) は mod で O(1)）/ 係数の繰り上がり", "O(|s|) / O(1)");
    SHOW(str_mod("123456789012345678901234567890", 9), "桁和と同じ（9 で割った余り）");
    SHOW(str_num_compare("100000000000000000000", "-5"), "符号も見る");
    SHOW(str_num_compare("0012", "12"), "先頭の 0 は無視");
    SHOW(substr_value("20261010", 4, 8), "部分文字列 [4, 8)");
    RUN(DigitPrefix P("31415926535", 1000000007));
    SHOW(P.value(2, 7), "\"41592\"");
    SHOW(str_mod("ZZ", 1000, 36), "36 進数の ZZ = 1295 の mod");
    note("carry_to_string：桁ごとの係数（下の桁から）を繰り上げて文字列に。大きな数の掛け算を畳み込みでするときの最後の一歩");
    SHOW(carry_to_string({18, 27, 28, 13, 4}), "123 × 456：桁 {3, 2, 1} と {6, 5, 4} の畳み込み → 56088");
    SHOW(carry_to_string({1, 1, 0, 1}, 2), "2 進数の 1011");

    api("【7】 digit_length_blocks(l, r) / total_digits_up_to(N) / count_digit_up_to(N, d) / digit_sum_total_up_to(N)", "vector / Int128",
        "[l, r) を桁数ごとの区間に分ける / 1〜N の桁数の合計・数字 d の個数・桁和の合計", "O(桁数)");
    SHOW(total_digits_up_to(100), "\"123…100\" の長さ（9 + 180 + 3）");
    SHOW(count_digit_up_to(100, 0), "1〜100 の 0 の個数（10, 20, …, 100）");
    SHOW(digit_sum_total_up_to(10), "1 + 2 + … + 9 + 1");
    {
        // 典型90 082: Σ_{x=L}^{R} x × (x の桁数) mod 10^9+7
        const long long L = 98, R = 100, MOD = 1000000007;
        long long ans = 0;
        for (const DigitLengthBlock& b : digit_length_blocks(L, R + 1))
            ans = (ans + mod_mul(b.len, arith_sum_mod(b.lo, 1, b.hi - b.lo, MOD), MOD)) % MOD;
        note("典型90 082（Counting Numbers）: 桁数ごとの区間に分けて、桁数 × 区間の和（arith_sum_mod）");
        std::cout << "    L = 98, R = 100 → " << ans << "（98·2 + 99·2 + 100·3）\n";
    }

    api("【8】 numbers_from_digits(ds, max_len)", "vector<long long>", "数字 ds だけでできた 1〜max_len 桁の数（小さい順。0 が ds にあれば 0 も）", "O(個数)");
    SHOW(numbers_from_digits({1, 3}, 2), "1, 3, 11, 13, 31, 33");
    {
        // ABC114 C: N 以下で、数字が 3・5・7 だけで、3・5・7 をすべて含む数（七五三数）
        const long long N = 575;
        long long c = 0;
        for (long long x : numbers_from_digits({3, 5, 7}, num_digits(N))) {
            const std::string s = std::to_string(x);
            if (x <= N && s.find('3') != std::string::npos && s.find('5') != std::string::npos && s.find('7') != std::string::npos) ++c;
        }
        note("ABC114 C（755）: numbers_from_digits({3, 5, 7}, 桁数) から条件を満たすものを数える");
        std::cout << "    N = 575 → " << c << "\n";
    }
}

}  // namespace dig_demo

{SETUP_CONSOLE}

int main() {
    setup_console();
    dig_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    dig_test::benchmark();    // 速度計測
    dig_demo::demo();         // 使い方デモ
    return 0;
}
