#ifndef GG_DIGITS
#define GG_DIGITS
// ---- 桁・進数・回文数（数学第二世代）  桁和・n 進数・回文数（文字列の大きなものも）・数をつなげる・大きな数の文字列・繰り上がり・桁の数え上げ ----
namespace math_internal {
inline int digit_value(char c) {  // '0'〜'9'・'A'〜'Z'・'a'〜'z' の値（それ以外は -1）
    if ('0' <= c && c <= '9') return c - '0';
    if ('A' <= c && c <= 'Z') return c - 'A' + 10;
    if ('a' <= c && c <= 'z') return c - 'a' + 10;
    return -1;
}
inline void check_base(long long base, const char* func) {
    if (base < 2) fail(func, "base ≥ 2 で使う");
}
inline long long checked_mod(long long m, const char* func) {  // m ≥ 1 を確かめて返す（% m の前に）
    if (m < 1) fail(func, "m ≥ 1 で使う");
    return m;
}
inline Int128 pow10_128(int k) {
    Int128 r = 1;
    while (k-- > 0) r *= 10;
    return r;
}
inline std::string normalize_nonneg_str(const std::string& s, const char* func) {  // 先頭の + と余分な 0 を取る（0 以上の 10 進でなければ止める）
    size_t i = (!s.empty() && s[0] == '+') ? 1 : 0;
    if (i == s.size()) fail(func, "数字がない");
    for (size_t j = i; j < s.size(); ++j)
        if (s[j] < '0' || s[j] > '9') fail(func, "0 以上の 10 進の数として読めない");
    while (i + 1 < s.size() && s[i] == '0') ++i;
    return s.substr(i);
}
inline std::string str_inc(std::string s) {  // 0 以上の 10 進の文字列 + 1
    int i = (int)s.size() - 1;
    while (i >= 0 && s[i] == '9') s[i--] = '0';
    if (i < 0) s.insert(s.begin(), '1');
    else ++s[i];
    return s;
}
inline std::string str_dec(std::string s) {  // 1 以上の 10 進の文字列 − 1（先頭の 0 を取る）
    int i = (int)s.size() - 1;
    while (s[i] == '0') s[i--] = '9';
    --s[i];
    if (s.size() > 1 && s[0] == '0') s.erase(s.begin());
    return s;
}
inline Int128 str_to_128(const std::string& s) {  // 38 桁までの 0 以上の 10 進
    Int128 v = 0;
    for (char c : s) v = v * 10 + (c - '0');
    return v;
}
}  // namespace math_internal

// ---- 桁 ----
inline std::vector<int> digits_of(long long n, long long base = 10) {  // |n| の各桁（上の桁から。0 なら {0}。base ≤ 2^31 − 1）
    math_internal::check_base(base, "digits_of");
    if (base > 2147483647LL) math_internal::fail("digits_of", "桁を int で返すので base ≤ 2^31 − 1 で使う");
    unsigned long long u = math_internal::uabs(n);
    std::vector<int> ds;
    do ds.push_back((int)(u % (unsigned long long)base)), u /= (unsigned long long)base;
    while (u);
    std::reverse(ds.begin(), ds.end());
    return ds;
}
inline long long digit_sum(long long n, long long base = 10) {  // |n| の桁和
    math_internal::check_base(base, "digit_sum");
    unsigned long long u = math_internal::uabs(n), s = 0;
    for (; u; u /= (unsigned long long)base) s += u % (unsigned long long)base;
    return (long long)s;
}
inline int num_digits(long long n, long long base = 10) {  // |n| の桁数（0 は 1 桁）
    math_internal::check_base(base, "num_digits");
    unsigned long long u = math_internal::uabs(n);
    int d = 1;
    while (u >= (unsigned long long)base) u /= (unsigned long long)base, ++d;
    return d;
}
inline long long from_digits(const std::vector<int>& ds, long long base = 10) {  // 上の桁からの並びを数に（溢れたら止める）
    math_internal::check_base(base, "from_digits");
    long long v = 0;
    for (int d : ds) {
        if (d < 0 || d >= base) math_internal::fail("from_digits", "桁が 0 以上 base 未満でない");
        if (__builtin_mul_overflow(v, base, &v) || __builtin_add_overflow(v, (long long)d, &v)) math_internal::fail("from_digits", "long long に収まらない");
    }
    return v;
}

// ---- n 進数（2〜36 進と −2〜−36 進。10 以上の桁は A〜Z） ----
inline std::string to_base(long long n, int base) {
    if (!((2 <= base && base <= 36) || (-36 <= base && base <= -2))) math_internal::fail("to_base", "base は 2〜36 か −36〜−2");
    static const char* D = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    if (n == 0) return "0";
    std::string s;
    if (base > 0) {
        unsigned long long u = math_internal::uabs(n);
        for (; u; u /= (unsigned long long)base) s += D[u % (unsigned long long)base];
        if (n < 0) s += '-';
    } else {  // 負の進数：余りを 0〜|base|-1 にそろえる
        Int128 x = n;
        while (x != 0) {
            Int128 r = x % base;
            x /= base;
            if (r < 0) r -= base, x += 1;
            s += D[(int)r];
        }
    }
    std::reverse(s.begin(), s.end());
    return s;
}
inline long long from_base(const std::string& s, int base) {  // to_base の逆（小文字も可。正の進数なら先頭の - も可。溢れたら止める）
    if (!((2 <= base && base <= 36) || (-36 <= base && base <= -2))) math_internal::fail("from_base", "base は 2〜36 か −36〜−2");
    size_t i = 0;
    const bool neg = base > 0 && !s.empty() && s[0] == '-';
    if (neg || (!s.empty() && s[0] == '+')) ++i;
    if (i == s.size()) math_internal::fail("from_base", "数字がない");
    Int128 v = 0;
    const Int128 lim = (Int128)9223372036854775807LL + 1;
    for (; i < s.size(); ++i) {
        const int d = math_internal::digit_value(s[i]);
        if (d < 0 || d >= (base < 0 ? -base : base)) math_internal::fail("from_base", "その進数の桁でない文字がある");
        v = v * base + d;
        if (v > lim || v < -lim) math_internal::fail("from_base", "long long に収まらない");
    }
    if (neg) v = -v;
    if (v > (Int128)9223372036854775807LL || v < -lim) math_internal::fail("from_base", "long long に収まらない");
    return (long long)v;
}
// s（0〜9・A〜Z）を base 進数（2 以上。36 を超えてもよい）として読んだ値と cap の小さい方（溢れない。0 以上の数として読む）
inline long long from_base_cap(const std::string& s, long long base, long long cap) {
    math_internal::check_base(base, "from_base_cap");
    if (cap < 0) math_internal::fail("from_base_cap", "cap ≥ 0 で使う");
    if (s.empty()) math_internal::fail("from_base_cap", "数字がない");
    Int128 v = 0;
    for (char c : s) {
        const int d = math_internal::digit_value(c);
        if (d < 0 || d >= base) math_internal::fail("from_base_cap", "その進数の桁でない文字がある");
        if (v < cap) v = std::min<Int128>(v * base + d, cap);  // cap に届いたらそのまま（桁の確認だけ続ける）
    }
    return (long long)v;
}
inline bool base_value_leq(const std::string& s, long long base, long long limit) {  // s を base 進数で読んだ値 ≤ limit か（溢れない。ABC192 D）
    math_internal::check_base(base, "base_value_leq");
    if (s.empty()) math_internal::fail("base_value_leq", "数字がない");
    Int128 v = 0;
    bool over = limit < 0;
    for (char c : s) {
        const int d = math_internal::digit_value(c);
        if (d < 0 || d >= base) math_internal::fail("base_value_leq", "その進数の桁でない文字がある");
        if (!over) v = v * base + d, over = v > limit;
    }
    return !over;
}

// ---- 回文数（0 以上の整数。k 番目は 0 から数える：0 番目が 0） ----
inline bool is_palindrome_num(long long n, long long base = 10) {  // n ≥ 0 が base 進数で回文か
    if (n < 0) return false;
    const std::vector<int> d = digits_of(n, base);
    for (size_t i = 0, j = d.size() - 1; i < j; ++i, --j)
        if (d[i] != d[j]) return false;
    return true;
}
inline long long count_palindromes_of_length(int len) {  // len 桁の回文数の個数（1 桁は 0〜9 の 10 個。len ≤ 38 まで）
    if (len < 1) return 0;
    if (len == 1) return 10;
    return 9 * pow10_ll((len + 1) / 2 - 1);
}
inline std::string kth_palindrome(long long k) {  // 小さい方から k 番目の回文数（0 番目は 0）。10^18 を超えるので文字列
    if (k < 0) math_internal::fail("kth_palindrome", "k ≥ 0 で使う");
    if (k < 10) return std::string(1, char('0' + k));
    k -= 10;
    for (int len = 2;; ++len) {
        const long long c = count_palindromes_of_length(len);
        if (k < c) {
            const std::string half = std::to_string(pow10_ll((len + 1) / 2 - 1) + k);
            std::string rev(half.rbegin() + (len % 2), half.rend());
            return half + rev;
        }
        k -= c;
    }
}
inline long long make_palindrome(long long half, bool odd) {  // half（≥ 0）を左半分にした回文数（odd なら真ん中の桁を共有：123 → 12321 / 123321）
    if (half < 0) math_internal::fail("make_palindrome", "half ≥ 0 で使う");
    Int128 v = half;
    for (long long t = odd ? half / 10 : half; t > 0; t /= 10) {
        v = v * 10 + t % 10;
        if (v > (Int128)9223372036854775807LL) math_internal::fail("make_palindrome", "long long に収まらない");
    }
    return (long long)v;
}
inline long long count_palindromes_up_to(long long n) {  // 0 以上 n 以下の回文数の個数（n < 0 なら 0）
    if (n < 0) return 0;
    const int len = num_digits(n);
    long long c = 0;
    for (int l = 1; l < len; ++l) c += count_palindromes_of_length(l);
    if (len == 1) return n + 1;
    const long long hdiv = pow10_ll(len / 2), half = n / hdiv, lo = pow10_ll((len + 1) / 2 - 1);
    c += half - lo;  // 左半分が half より小さいもの
    if (make_palindrome(half, len % 2) <= n) ++c;
    return c;
}
inline long long count_palindromes_below(long long n) { return n <= 0 ? 0 : count_palindromes_up_to(n - 1); }  // [0, n) の回文数の個数
inline long long palindrome_rank(long long n) {  // 回文数 n が小さい方から何番目か（0 から数える。kth_palindrome の逆）
    if (!is_palindrome_num(n)) math_internal::fail("palindrome_rank", "n が回文数でない");
    return count_palindromes_below(n);
}
inline long long ceil_palindrome(long long n) {  // n 以上で最小の回文数（long long に収まらなければ止める）
    if (n <= 0) return 0;
    const int len = num_digits(n);
    if (len == 1) return n;
    const long long half = n / pow10_ll(len / 2);
    const long long p = make_palindrome(half, len % 2);
    if (p >= n) return p;
    if (half + 1 == pow10_ll(num_digits(half))) {  // 99…9 → 次の桁数の最小 10…01
        if (len >= 19) math_internal::fail("ceil_palindrome", "long long に収まらない");
        return pow10_ll(len) + 1;
    }
    return make_palindrome(half + 1, len % 2);
}
inline long long floor_palindrome(long long n) {  // n 以下で最大の回文数（n ≥ 0）
    if (n < 0) math_internal::fail("floor_palindrome", "n ≥ 0 で使う");
    const int len = num_digits(n);
    if (len == 1) return n;
    const long long half = n / pow10_ll(len / 2);
    const long long p = make_palindrome(half, len % 2);
    if (p <= n) return p;
    if (half == pow10_ll(num_digits(half) - 1)) return pow10_ll(len - 1) - 1;  // 10…0 → 前の桁数の最大 99…9
    return make_palindrome(half - 1, len % 2);
}

// ---- 大きな回文数（10 進の文字列。0 以上。先頭の 0 は無視。long long を超える数に） ----
inline bool is_palindrome_str(const std::string& s) {  // 文字列 s が回文か（数でなくてもよい）
    for (size_t i = 0, j = s.size(); i + 1 < j; ++i, --j)
        if (s[i] != s[j - 1]) return false;
    return true;
}
inline std::string make_palindrome_str(const std::string& half, bool odd) {  // half を左半分にした回文（odd なら真ん中の文字を共有："123" → "12321" / "123321"）
    if (half.empty()) math_internal::fail("make_palindrome_str", "half が空");
    std::string r = half;
    r.append(half.rbegin() + (odd ? 1 : 0), half.rend());
    return r;
}
inline std::string ceil_palindrome_str(const std::string& n) {  // n 以上で最小の回文数
    const std::string s = math_internal::normalize_nonneg_str(n, "ceil_palindrome_str");
    const size_t len = s.size();
    if (len == 1) return s;
    const std::string half = s.substr(0, (len + 1) / 2);
    const std::string p = make_palindrome_str(half, len % 2);
    if (p >= s) return p;  // 同じ長さなので文字列の比較 = 数の比較
    if (half.find_first_not_of('9') == std::string::npos) return "1" + std::string(len - 1, '0') + "1";  // 99…9 → 10…01
    return make_palindrome_str(math_internal::str_inc(half), len % 2);
}
inline std::string floor_palindrome_str(const std::string& n) {  // n 以下で最大の回文数
    const std::string s = math_internal::normalize_nonneg_str(n, "floor_palindrome_str");
    const size_t len = s.size();
    if (len == 1) return s;
    const std::string half = s.substr(0, (len + 1) / 2);
    const std::string p = make_palindrome_str(half, len % 2);
    if (p <= s) return p;
    if (half == "1" + std::string(half.size() - 1, '0')) return std::string(len - 1, '9');  // 10…0 → 99…9
    return make_palindrome_str(math_internal::str_dec(half), len % 2);
}
inline Int128 count_palindromes_up_to_str(const std::string& n) {  // [0, n] の回文数の個数（n は 75 桁まで）
    const std::string s = math_internal::normalize_nonneg_str(n, "count_palindromes_up_to_str");
    const int len = (int)s.size();
    if (len > 75) math_internal::fail("count_palindromes_up_to_str", "Int128 に収まるよう 75 桁まで");
    if (len == 1) return s[0] - '0' + 1;
    Int128 c = 10;
    for (int l = 2; l < len; ++l) c += 9 * math_internal::pow10_128((l + 1) / 2 - 1);
    const std::string half = s.substr(0, (len + 1) / 2);
    c += math_internal::str_to_128(half) - math_internal::pow10_128((len + 1) / 2 - 1);  // 左半分が half より小さいもの
    if (make_palindrome_str(half, len % 2) <= s) ++c;
    return c;
}
inline Int128 count_palindromes_below_str(const std::string& n) {  // [0, n) の回文数の個数（n は 75 桁まで）
    const std::string s = math_internal::normalize_nonneg_str(n, "count_palindromes_below_str");
    return count_palindromes_up_to_str(s) - (is_palindrome_str(s) ? 1 : 0);
}

// ---- 数をつなげる（a, b ≥ 0） ----
inline Int128 concat_value(long long a, long long b) {  // a の後ろに b をつなげた数（"12" + "345" → 12345）
    if (a < 0 || b < 0) math_internal::fail("concat_value", "a, b ≥ 0 で使う");
    return (Int128)a * math_internal::pow10_128(num_digits(b)) + b;
}
inline Int128 repeat_concat_value(long long x, long long k) {  // x（≥ 0）を k 回つなげた数（Int128 に収まらなければ止める。k = 0 なら 0）
    if (x < 0 || k < 0) math_internal::fail("repeat_concat_value", "x, k ≥ 0 で使う");
    if (x == 0) return 0;
    const Int128 p = math_internal::pow10_128(num_digits(x)), LIM = (Int128)(~(UInt128)0 >> 1);
    Int128 v = 0;
    for (long long i = 0; i < k; ++i) {
        if (v > (LIM - x) / p) math_internal::fail("repeat_concat_value", "Int128 に収まらない");
        v = v * p + x;
    }
    return v;
}
inline long long pow10_mod(long long k, long long m) { return mod_pow(10, k, m); }  // 10^k mod m
inline long long concat_mod(long long a, long long b, long long m) {  // (a の後ろに b をつなげた数) mod m
    if (a < 0 || b < 0) math_internal::fail("concat_mod", "a, b ≥ 0 で使う");
    math_internal::checked_mod(m, "concat_mod");
    return (long long)(((Int128)mod_mul(a, pow10_mod(num_digits(b), m), m) + b % m) % m);
}
inline long long repeat_concat_mod(long long x, long long k, long long m) {  // x（≥ 0）を k 回つなげた数 mod m（ABC357 D）
    if (x < 0 || k < 0) math_internal::fail("repeat_concat_mod", "x, k ≥ 0 で使う");
    math_internal::checked_mod(m, "repeat_concat_mod");
    return mod_mul(x % m, geom_sum_mod(pow10_mod(num_digits(x), m), k, m), m);
}
inline long long repunit_mod(long long k, long long m) { return geom_sum_mod(10, k, m); }  // 1 を k 個並べた数 mod m

// ---- 大きな数の文字列（10 進。先頭の - と 0 があってもよい） ----
inline long long str_mod(const std::string& s, long long m, int base = 10) {  // 値 mod m（[0, m)）。base 進数（2〜36。10 以上の桁は A〜Z / a〜z）としても読める
    if (m < 1) math_internal::fail("str_mod", "m ≥ 1 で使う");
    if (base < 2 || base > 36) math_internal::fail("str_mod", "base は 2〜36");
    size_t i = (!s.empty() && (s[0] == '-' || s[0] == '+')) ? 1 : 0;
    if (i == s.size()) math_internal::fail("str_mod", "数字がない");
    long long r = 0;
    for (; i < s.size(); ++i) {
        const int d = math_internal::digit_value(s[i]);
        if (d < 0 || d >= base) math_internal::fail("str_mod", "数として読めない");
        r = (long long)(((Int128)r * base + d) % m);
    }
    return s[0] == '-' ? (m - r) % m : r;
}
// 繰り上がり：c[i] が base^i の係数（下の桁から。畳み込みの結果など。負もよいが全体は 0 以上）の数を base 進の文字列（上の桁から）に
template <class V>
std::string carry_to_string(const V& c, int base = 10) {
    if (base < 2 || base > 36) math_internal::fail("carry_to_string", "base は 2〜36");
    static const char* D = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    std::string s;
    Int128 carry = 0;
    auto push = [&](Int128 v) {  // v を足して 1 桁出す（余りは 0〜base−1、商は切り捨て）
        Int128 r = v % base, q = v / base;
        if (r < 0) r += base, --q;
        s += D[(int)r];
        carry = q;
    };
    for (const auto& x : c) push(carry + (Int128)x);
    while (carry > 0) push(carry);
    if (carry < 0) math_internal::fail("carry_to_string", "値が負");
    while (s.size() > 1 && s.back() == '0') s.pop_back();
    if (s.empty()) s = "0";
    std::reverse(s.begin(), s.end());
    return s;
}
inline std::string carry_to_string(std::initializer_list<long long> c, int base = 10) { return carry_to_string<std::initializer_list<long long>>(c, base); }
inline int str_num_compare(const std::string& a, const std::string& b) {  // 10 進の整数の文字列の大小（-1 / 0 / 1）
    auto parse = [](const std::string& s, bool& neg) {
        size_t i = (!s.empty() && (s[0] == '-' || s[0] == '+')) ? 1 : 0;
        neg = !s.empty() && s[0] == '-';
        while (i + 1 < s.size() && s[i] == '0') ++i;
        std::string d = s.substr(i);
        if (d.empty()) math_internal::fail("str_num_compare", "数字がない");
        for (char c : d)
            if (c < '0' || c > '9') math_internal::fail("str_num_compare", "数として読めない");
        if (d == "0") neg = false;
        return d;
    };
    bool na, nb;
    const std::string x = parse(a, na), y = parse(b, nb);
    if (na != nb) return na ? -1 : 1;
    int c = x.size() != y.size() ? (x.size() < y.size() ? -1 : 1) : (x < y ? -1 : (x > y ? 1 : 0));
    return na ? -c : c;
}
inline long long substr_value(const std::string& s, int l, int r) {  // 数字の文字列 s の [l, r) の値（l ≥ r なら 0。溢れたら止める）
    if (l < r && (l < 0 || r > (int)s.size())) math_internal::fail("substr_value", "範囲が文字列の外");
    long long v = 0;
    for (int i = l; i < r; ++i) {
        if (s[i] < '0' || s[i] > '9') math_internal::fail("substr_value", "数字でない文字がある");
        if (__builtin_mul_overflow(v, 10LL, &v) || __builtin_add_overflow(v, (long long)(s[i] - '0'), &v)) math_internal::fail("substr_value", "long long に収まらない");
    }
    return v;
}
struct DigitPrefix {  // 数字の文字列 s の部分文字列 [l, r) の値を mod m で O(1)
    long long m;
    std::vector<long long> pre, pw;  // pre[i] = s[0, i) の値 mod m、pw[i] = 10^i mod m
    DigitPrefix(const std::string& s, long long m_) : m(math_internal::checked_mod(m_, "DigitPrefix")), pre(s.size() + 1, 0), pw(s.size() + 1, 1 % m) {
        for (size_t i = 0; i < s.size(); ++i) {
            if (s[i] < '0' || s[i] > '9') math_internal::fail("DigitPrefix", "数字でない文字がある");
            pre[i + 1] = (long long)(((Int128)pre[i] * 10 + (s[i] - '0')) % m);
            pw[i + 1] = (long long)((Int128)pw[i] * 10 % m);
        }
    }
    long long value(int l, int r) const {  // s[l, r) の値 mod m（l ≥ r なら 0。0 ≤ l、r ≤ |s|）
        if (l >= r) return 0;
        if (l < 0 || r >= (int)pre.size()) math_internal::fail("DigitPrefix::value", "範囲が文字列の外");
        return (long long)((((Int128)pre[r] - (Int128)pre[l] * pw[r - l]) % m + m) % m);
    }
};

// ---- 桁の数え上げ（1〜N） ----
struct DigitLengthBlock {  // 桁数が len の数の範囲 [lo, hi)
    int len;
    long long lo, hi;
};
inline std::vector<DigitLengthBlock> digit_length_blocks(long long l, long long r) {  // [l, r)（0 ≤ l）を桁数ごとに分ける
    std::vector<DigitLengthBlock> res;
    if (l < 0) math_internal::fail("digit_length_blocks", "l ≥ 0 で使う");
    for (long long x = l; x < r;) {
        const int len = num_digits(x);
        const long long end = len >= 19 ? r : std::min(r, pow10_ll(len));
        res.push_back({len, x, end});
        x = end;
    }
    return res;
}
inline Int128 total_digits_up_to(long long N) {  // 1〜N の桁数の合計（"123…N" の長さ）
    Int128 s = 0;
    if (N < 1) return 0;
    for (const DigitLengthBlock& b : digit_length_blocks(1, N + (N < std::numeric_limits<long long>::max()))) s += (Int128)b.len * (b.hi - b.lo);
    if (N == std::numeric_limits<long long>::max()) s += 19;
    return s;
}
inline Int128 count_digit_up_to(long long N, int d) {  // 1〜N を書いたとき数字 d（0〜9）が出る回数（先頭の 0 は数えない）
    if (d < 0 || d > 9) math_internal::fail("count_digit_up_to", "d は 0〜9");
    if (N < 1) return 0;
    Int128 c = 0;
    for (Int128 p = 1; p <= N; p *= 10) {  // 各位について
        const Int128 high = N / (p * 10), cur = N / p % 10, low = N % p;
        if (d == 0) {
            if (high == 0) continue;  // 先頭の位
            c += (high - 1) * p + (cur > 0 ? p : low + 1);
        } else {
            c += high * p + (cur > d ? p : (cur == d ? low + 1 : 0));
        }
    }
    return c;
}
inline Int128 digit_sum_total_up_to(long long N) {  // 1〜N の桁和の合計
    Int128 s = 0;
    for (int d = 1; d <= 9; ++d) s += count_digit_up_to(N, d) * d;
    return s;
}

// ---- 使える数字が決まった数の列挙 ----
// 数字の集合 ds（0〜9）だけを使う 1〜max_len 桁の数（先頭は 0 以外。0 が ds にあれば 0 も）を小さい順に。max_len ≤ 18
inline std::vector<long long> numbers_from_digits(std::vector<int> ds, int max_len) {
    if (max_len < 1 || max_len > 18) math_internal::fail("numbers_from_digits", "1 ≤ max_len ≤ 18 で使う");
    for (int d : ds)
        if (d < 0 || d > 9) math_internal::fail("numbers_from_digits", "数字は 0〜9");
    std::sort(ds.begin(), ds.end());
    ds.erase(std::unique(ds.begin(), ds.end()), ds.end());
    std::vector<long long> res, cur;
    if (!ds.empty() && ds[0] == 0) res.push_back(0);
    for (int d : ds)
        if (d > 0) cur.push_back(d);
    for (int len = 1; len <= max_len && !cur.empty(); ++len) {
        res.insert(res.end(), cur.begin(), cur.end());
        if (len == max_len) break;
        if (cur.size() * ds.size() > 50000000) math_internal::fail("numbers_from_digits", "数が多すぎる（5×10^7 個まで）");
        std::vector<long long> nxt;
        nxt.reserve(cur.size() * ds.size());
        for (long long x : cur)
            for (int d : ds) nxt.push_back(x * 10 + d);
        cur.swap(nxt);
    }
    return res;
}
#endif  // GG_DIGITS
