#ifndef GG_FRACTION
#define GG_FRACTION
// ---- 分数（数学第二世代）  FractionT<T>：いつも既約・分母 > 0。Fraction = long long、Fraction128 = Int128 ----
namespace math_internal {
inline Int128 frac_abs(Int128 a) { return a < 0 ? -a : a; }
inline int frac_ctz(UInt128 x) {  // x ≠ 0
    const unsigned long long lo = (unsigned long long)x;
    return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((unsigned long long)(x >> 64));
}
inline UInt128 frac_uabs(Int128 a) { return a < 0 ? (UInt128)0 - (UInt128)a : (UInt128)a; }  // −2^127 でも溢れない
inline Int128 frac_gcd(Int128 sa, Int128 sb) {  // gcd(|a|, |b|)（128 ビットは二進 GCD、64 ビットに入ったら速い方で。結果 < 2^127 のときに使う）
    UInt128 a = frac_uabs(sa), b = frac_uabs(sb);
    if (a == 0) return (Int128)b;
    if (b == 0) return (Int128)a;
    const int za = frac_ctz(a), zb = frac_ctz(b), k = za < zb ? za : zb;
    a >>= za, b >>= zb;
    while (true) {  // a, b は奇数
        if ((a >> 64) == 0 && (b >> 64) == 0) return (Int128)((UInt128)ugcd((unsigned long long)a, (unsigned long long)b) << k);
        if (a > b) std::swap(a, b);
        b -= a;
        if (b == 0) return (Int128)(a << k);
        b >>= frac_ctz(b);
    }
}
inline Int128 frac_mul(Int128 a, Int128 b) {
    Int128 r;
    if (__builtin_mul_overflow(a, b, &r)) fail("Fraction", "途中の値が 128 ビットに収まらない");
    return r;
}
inline Int128 frac_add(Int128 a, Int128 b) {
    Int128 r;
    if (__builtin_add_overflow(a, b, &r)) fail("Fraction", "途中の値が 128 ビットに収まらない");
    return r;
}
// a/b と c/d（b, d > 0）の大小（-1 / 0 / 1）。掛け算をしない（連分数のように整数部分を比べていく）
inline int frac_cmp_euclid(Int128 a, Int128 b, Int128 c, Int128 d) {
    int sign = 1;
    while (true) {
        Int128 q1 = a / b, r1 = a % b, q2 = c / d, r2 = c % d;  // 掛け算をしないので溢れない
        if (r1 < 0) q1 -= 1, r1 += b;                          // 床関数の商と 0 ≤ r1 < b に
        if (r2 < 0) q2 -= 1, r2 += d;
        if (q1 != q2) return q1 < q2 ? -sign : sign;
        if (r1 == 0 || r2 == 0) return r1 == r2 ? 0 : (r1 == 0 ? -sign : sign);
        a = b, b = r1, c = d, d = r2;  // r1/b と r2/d の比較は逆数 b/r1 と d/r2 の比較（向きが逆）
        sign = -sign;
    }
}
}  // namespace math_internal

template <class T>
struct FractionT {
    static_assert(std::is_same_v<T, long long> || std::is_same_v<T, Int128>, "FractionT は long long か Int128");
    T num, den;  // いつも既約で den > 0（直接書きかえない）

    FractionT() : num(0), den(1) {}
    template <class I, std::enable_if_t<std::is_integral_v<I> || std::is_same_v<I, Int128>, int> = 0>
    FractionT(I n) : num(0), den(1) {  // 整数から
        *this = make(to_i128(n), 1);
    }
    template <class I, class J, std::enable_if_t<(std::is_integral_v<I> || std::is_same_v<I, Int128>) && (std::is_integral_v<J> || std::is_same_v<J, Int128>), int> = 0>
    FractionT(I n, J d) : num(0), den(1) {  // n / d（d ≠ 0。約分する）
        if (d == 0) math_internal::fail("Fraction", "分母が 0");
        *this = make(to_i128(n), to_i128(d));
    }

    friend FractionT operator+(const FractionT& x, const FractionT& y) {
        using namespace math_internal;
        const Int128 g = frac_gcd(x.den, y.den);
        const Int128 t = frac_add(frac_mul(x.num, y.den / g), frac_mul(y.num, x.den / g));
        if (g == 1) return reduced(t, frac_mul(x.den, y.den));  // 分母が互いに素なら約分は要らない
        const Int128 g2 = frac_gcd(t, g);                      // gcd(t, 分母) = gcd(t, g)
        return reduced(t / g2, frac_mul(x.den / g, y.den / g2));
    }
    friend FractionT operator-(const FractionT& x, const FractionT& y) { return x + (-y); }
    friend FractionT operator*(const FractionT& x, const FractionT& y) {
        using namespace math_internal;
        const Int128 g1 = frac_gcd(x.num, y.den), g2 = frac_gcd(y.num, x.den);
        if (g1 == 0 || g2 == 0) return FractionT();  // どちらかが 0
        return reduced(frac_mul(x.num / g1, y.num / g2), frac_mul(x.den / g2, y.den / g1));  // 先に約分したので既約
    }
    friend FractionT operator/(const FractionT& x, const FractionT& y) { return x * y.reciprocal(); }
    FractionT operator-() const {
        FractionT r;
        r.num = -num, r.den = den;  // num は T の最小値にならないようにしてある
        return r;
    }
    FractionT operator+() const { return *this; }
    FractionT& operator+=(const FractionT& y) { return *this = *this + y; }
    FractionT& operator-=(const FractionT& y) { return *this = *this - y; }
    FractionT& operator*=(const FractionT& y) { return *this = *this * y; }
    FractionT& operator/=(const FractionT& y) { return *this = *this / y; }

    friend int compare(const FractionT& x, const FractionT& y) {  // −1 / 0 / 1
        if constexpr (std::is_same_v<T, long long>) {
            const Int128 l = (Int128)x.num * y.den, r = (Int128)y.num * x.den;
            return l < r ? -1 : (l > r ? 1 : 0);
        } else {
            const Int128 L = (Int128)1 << 63;  // 4 つとも 2^63 未満なら掛け算で
            if (-L < x.num && x.num < L && x.den < L && -L < y.num && y.num < L && y.den < L) {
                const Int128 l = x.num * y.den, r = y.num * x.den;
                return l < r ? -1 : (l > r ? 1 : 0);
            }
            return math_internal::frac_cmp_euclid(x.num, x.den, y.num, y.den);
        }
    }
    friend bool operator==(const FractionT& x, const FractionT& y) { return x.num == y.num && x.den == y.den; }
    friend bool operator!=(const FractionT& x, const FractionT& y) { return !(x == y); }
    friend bool operator<(const FractionT& x, const FractionT& y) { return compare(x, y) < 0; }
    friend bool operator>(const FractionT& x, const FractionT& y) { return compare(x, y) > 0; }
    friend bool operator<=(const FractionT& x, const FractionT& y) { return compare(x, y) <= 0; }
    friend bool operator>=(const FractionT& x, const FractionT& y) { return compare(x, y) >= 0; }

    T floor() const { return (T)math_internal::floor_div128(num, den); }  // ⌊x⌋
    T ceil() const { return (T)math_internal::ceil_div128(num, den); }    // ⌈x⌉
    T round() const {                                                    // 四捨五入（ちょうど .5 は 0 から遠い方）
        const T q = num / den, r = num % den;                             // 0 に向かって切り捨て
        const T ar = r < 0 ? -r : r;
        if (ar >= den - ar) return r < 0 ? q - 1 : q + 1;
        return q;
    }
    int sign() const { return num > 0 ? 1 : (num < 0 ? -1 : 0); }
    FractionT abs() const { return num < 0 ? -*this : *this; }
    FractionT reciprocal() const {  // 1/x（0 なら止める）
        if (num == 0) math_internal::fail("Fraction::reciprocal", "0 の逆数");
        FractionT r;
        r.num = num < 0 ? -den : den, r.den = num < 0 ? -num : num;
        return r;
    }
    bool is_integer() const { return den == 1; }
    long double to_real() const { return (long double)num / (long double)den; }
    std::string to_decimal(int k) const {  // 小数第 k 位まで（それより下は切り捨て。k = 0 なら整数部分だけ）
        if (k < 0) math_internal::fail("Fraction::to_decimal", "k ≥ 0 で使う");
        const Int128 a = math_internal::frac_abs(num);
        const Int128 ip = a / den;
        Int128 r = a % den;
        std::string digits;
        const UInt128 D = (UInt128)den;
        for (int i = 0; i < k; ++i) {  // 10 r を den で割る（r を 10 回足しながら引く。10 r が 128 ビットを超えても大丈夫）
            UInt128 acc = 0, R = (UInt128)r;
            int dg = 0;
            for (int t = 0; t < 10; ++t) {
                acc += R;  // acc < den、R < den なので 2^128 未満
                if (acc >= D) acc -= D, ++dg;
            }
            digits += (char)('0' + dg), r = (Int128)acc;
        }
        const bool zero = ip == 0 && digits.find_first_not_of('0') == std::string::npos;
        return (num < 0 && !zero ? "-" : "") + int128_to_string(ip) + (k ? "." + digits : "");
    }
    friend std::ostream& operator<<(std::ostream& os, const FractionT& x) {  // "p/q"（整数なら "p"）
        os << int128_to_string(x.num);
        if (x.den != 1) os << '/' << int128_to_string(x.den);
        return os;
    }
    friend std::istream& operator>>(std::istream& is, FractionT& x) {  // "p/q"・"p"・"-1.25" を読む
        std::string s;
        if (!(is >> s)) return is;
        x = parse(s);
        return is;
    }
    static FractionT parse(const std::string& s) {  // "p/q"・"p"・"-1.25"・".5" を読む
        using namespace math_internal;
        const size_t slash = s.find('/');
        if (slash != std::string::npos) return FractionT(parse_int128(s.substr(0, slash)), parse_int128(s.substr(slash + 1)));
        const size_t dot = s.find('.');
        if (dot == std::string::npos) return FractionT(parse_int128(s));
        const size_t st = (!s.empty() && (s[0] == '-' || s[0] == '+')) ? 1 : 0;
        Int128 n = 0, d = 1;
        int cnt = 0;  // 数字の数
        for (size_t i = st; i < s.size(); ++i) {
            if (i == dot) continue;
            if (s[i] < '0' || s[i] > '9') fail("Fraction::parse", "数の形が正しくない");
            n = frac_add(frac_mul(n, 10), s[i] - '0'), ++cnt;
            if (i > dot) d = frac_mul(d, 10);
        }
        if (cnt == 0) fail("Fraction::parse", "数字がない");
        return FractionT(st == 1 && s[0] == '-' ? -n : n, d);
    }

  private:
    template <class I>
    static Int128 to_i128(I v) {  // 整数を Int128 に（2^127 以上の UInt128 は止める）
        if constexpr (sizeof(I) == 16 && !std::is_same_v<I, Int128>) {
            if (v > (I)((UInt128)1 << 127) - 1) math_internal::fail("Fraction", "値が型に収まらない");
        }
        return (Int128)v;
    }
    static FractionT reduced(Int128 n, Int128 d) {  // 既約で d > 0 と分かっているもの（T に入るかだけ確かめる）
        if constexpr (std::is_same_v<T, long long>) {
            const Int128 L = (Int128)9223372036854775807LL;
            if (n > L || n < -L || d > L) math_internal::fail("Fraction", "値が long long に収まらない（Fraction128 を使う）");
        } else {
            const Int128 MIN128 = (Int128)((UInt128)1 << 127);
            if (n == MIN128) math_internal::fail("Fraction", "値が型に収まらない");
        }
        FractionT r;
        r.num = (T)n, r.den = (T)d;
        return r;
    }
    static FractionT make(Int128 n, Int128 d) {  // 約分して分母を正に。T に入らなければ止める
        using namespace math_internal;
        if (d == 0) fail("Fraction", "分母が 0");
        const Int128 MIN128 = (Int128)((UInt128)1 << 127);
        if (n == MIN128 || d == MIN128) fail("Fraction", "値が型に収まらない");
        if (d < 0) n = -n, d = -d;
        const Int128 g = frac_gcd(n, d);
        if (g > 1) n /= g, d /= g;
        if (n == 0) d = 1;
        if constexpr (std::is_same_v<T, long long>) {
            const Int128 L = (Int128)9223372036854775807LL;
            if (n > L || n < -L || d > L) fail("Fraction", "値が long long に収まらない（Fraction128 を使う）");
        }
        FractionT r;
        r.num = (T)n, r.den = (T)d;
        return r;
    }
};
using Fraction = FractionT<long long>;
using Fraction128 = FractionT<Int128>;
namespace std {
template <class T>
struct hash<FractionT<T>> {  // unordered_map<Fraction, V> でそのまま使える
    size_t operator()(const FractionT<T>& x) const {
        const unsigned long long a = (unsigned long long)(UInt128)x.num, b = (unsigned long long)((UInt128)x.num >> 64);
        const unsigned long long c = (unsigned long long)(UInt128)x.den, d = (unsigned long long)((UInt128)x.den >> 64);
        unsigned long long h = a * 0x9E3779B97F4A7C15ULL ^ (b + 0x632BE59BD9B4E019ULL);
        h = (h ^ (h >> 29)) * 0xBF58476D1CE4E5B9ULL ^ c;
        h = (h ^ (h >> 31)) * 0x94D049BB133111EBULL ^ (d + 0x7F4A7C15ULL);
        return (size_t)(h ^ (h >> 32));
    }
};
}  // namespace std

// p/q（q ≠ 0）にいちばん近い、分母が max_den 以下の分数（同じ近さなら分母の小さい方、分母も同じなら値の小さい方）  O(log)
inline Fraction best_rational_approx(long long p, long long q, long long max_den) {
    if (q == 0 || max_den < 1) math_internal::fail("best_rational_approx", "q ≠ 0、max_den ≥ 1 で使う");
    Int128 P = p, Q = q;
    if (Q < 0) P = -P, Q = -Q;
    {
        const Int128 g = math_internal::frac_gcd(P, Q);
        P /= g, Q /= g;
        if (Q <= max_den) return Fraction(P, Q);
    }
    // 連分数の近似分数 h/k を分母が max_den を超える手前まで
    Int128 h0 = 0, h1 = 1, k0 = 1, k1 = 0, a = P, b = Q;
    while (b != 0) {
        const Int128 t = math_internal::floor_div128(a, b);
        const Int128 k2 = t * k1 + k0;
        if (k2 > max_den) break;
        const Int128 h2 = t * h1 + h0;
        h0 = h1, h1 = h2, k0 = k1, k1 = k2;
        const Int128 r = a - t * b;
        a = b, b = r;
    }
    // 候補：最後の近似分数 h1/k1 と、中間近似分数 (h0 + m h1)/(k0 + m k1)（m は分母が収まる最大）
    const Int128 m = (max_den - k0) / k1;
    const Int128 hs = h0 + m * h1, ks = k0 + m * k1;
    const Int128 d1 = math_internal::frac_abs(P * k1 - h1 * Q), d2 = math_internal::frac_abs(P * ks - hs * Q);  // 距離 × Q × 分母
    const Int128 l = d1 * ks, r = d2 * k1;  // d1/k1 と d2/ks を比べる
    if (l < r || (l == r && k1 <= ks)) return Fraction(h1, k1);
    return Fraction(hs, ks);
}
namespace math_internal {
inline void frac_bounded(long long p, long long q, long long max_den, const char* func, Int128& lp, Int128& lq, Int128& rp, Int128& rq) {
    if (q == 0 || max_den < 1) fail(func, "q ≠ 0、max_den ≥ 1 で使う");
    Int128 P = p, Q = q;
    if (Q < 0) P = -P, Q = -Q;
    const Int128 g = frac_gcd(P, Q);
    P /= g, Q /= g;
    if (Q <= max_den) {  // x そのもの
        lp = rp = P, lq = rq = Q;
        return;
    }
    // L < x < R を隣り合う分数（Stern–Brocot 木）のまま、分母 max_den 以下で x に近づける
    lp = floor_div128(P, Q), lq = 1, rp = lp + 1, rq = 1;
    while (true) {
        Int128 k1 = (P * lq - lp * Q) / (rp * Q - P * rq);  // L + k R ≤ x となる最大の k
        k1 = std::min(k1, (max_den - lq) / rq);
        lp += k1 * rp, lq += k1 * rq;
        Int128 k2 = (rp * Q - P * rq) / (P * lq - lp * Q);  // R + k L ≥ x となる最大の k
        k2 = std::min(k2, (max_den - rq) / lq);
        rp += k2 * lp, rq += k2 * lq;
        if (k1 == 0 && k2 == 0) return;
    }
}
}  // namespace math_internal
// p/q 以下で最大の、分母が max_den 以下の分数 / p/q 以上で最小のもの（p/q の分母が max_den 以下なら p/q）  O(log)
inline Fraction rational_floor(long long p, long long q, long long max_den) {
    Int128 lp, lq, rp, rq;
    math_internal::frac_bounded(p, q, max_den, "rational_floor", lp, lq, rp, rq);
    return Fraction(lp, lq);
}
inline Fraction rational_ceil(long long p, long long q, long long max_den) {
    Int128 lp, lq, rp, rq;
    math_internal::frac_bounded(p, q, max_den, "rational_ceil", lp, lq, rp, rq);
    return Fraction(rp, rq);
}
// x < p/q < y となる分数のうち、分母が最小のもの（分母が同じなら |分子| が最小。開区間。Stern–Brocot 木でいちばん上）  O(log)
inline Fraction simplest_between(const Fraction& x, const Fraction& y) {
    if (!(x < y)) math_internal::fail("simplest_between", "x < y で使う");
    if (x.sign() < 0 && y.sign() > 0) return Fraction(0);
    if (y.sign() <= 0) return -simplest_between(-y, -x);
    // 0 ≤ a/b < c/d の開区間で（分母最小 = 分子も最小）
    auto rec = [](auto&& self, Int128 a, Int128 b, Int128 c, Int128 d) -> std::pair<Int128, Int128> {
        const Int128 k = a / b;
        if ((k + 1) * d < c) return {k + 1, 1};  // 整数 k + 1 が入る
        a -= k * b, c -= k * d;                  // 0 ≤ a/b < c/d ≤ 1 に
        if (a == 0) {                            // (0, c/d) で分母最小：1/q（q = ⌊d/c⌋ + 1）
            const Int128 q = d / c + 1;
            return {k * q + 1, q};
        }
        const std::pair<Int128, Int128> r = self(self, d, c, b, a);  // 逆数をとる：d/c < P/Q < b/a
        return {k * r.first + r.second, r.first};
    };
    const std::pair<Int128, Int128> r = rec(rec, x.num, x.den, y.num, y.den);
    return Fraction(r.first, r.second);
}
// 連分数展開 p/q = a0 + 1/(a1 + 1/(a2 + …))（q ≠ 0。a0 は負もある。最後の項は 2 以上（p/q が整数なら [a0] だけ））  O(log)
inline std::vector<long long> continued_fraction(long long p, long long q) {
    if (q == 0) math_internal::fail("continued_fraction", "q ≠ 0 で使う");
    Int128 a = p, b = q;
    if (b < 0) a = -a, b = -b;
    std::vector<long long> cf;
    while (b != 0) {
        const Int128 t = math_internal::floor_div128(a, b);
        if (t > (Int128)LLONG_MAX) math_internal::fail("continued_fraction", "項が long long に収まらない（p/q = 2^63）");
        cf.push_back((long long)t);
        const Int128 r = a - t * b;
        a = b, b = r;
    }
    return cf;
}
// 連分数 [a0; a1, a2, …] の値（a1 以降は 1 以上。溢れたら止める）  O(長さ)
template <class V>
Fraction from_continued_fraction(const V& cfv) {
    std::vector<long long> cf;
    for (const auto& x : cfv) cf.push_back(math_internal::to_ll_checked(x, "from_continued_fraction"));
    if (cf.empty()) math_internal::fail("from_continued_fraction", "空の連分数");
    Fraction x(cf.back());
    for (size_t i = cf.size() - 1; i > 0; --i) {
        if (cf[i] < 1) math_internal::fail("from_continued_fraction", "a1 以降は 1 以上で使う");
        x = Fraction(cf[i - 1]) + x.reciprocal();
    }
    return x;
}
inline Fraction from_continued_fraction(std::initializer_list<long long> cf) { return from_continued_fraction<std::initializer_list<long long>>(cf); }
#endif  // GG_FRACTION
