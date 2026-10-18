#ifndef GG_BIGINT
#define GG_BIGINT
// ---- 多倍長整数 BigInt（数学第二世代）  符号つき。10^9 進で持つ。× は桁が多いとカラツバ、÷ は筆算（Knuth） ----
// 割り算は C++ の整数と同じ（0 の方向に切り捨て。% の符号は割られる数と同じ）。0 で割ると止める
struct BigInt {
    static constexpr unsigned BASE = 1000000000U;
    static constexpr int BASE_DIGITS = 9;

    BigInt() {}
    template <class T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    BigInt(T v) {  // int・long long・unsigned など（Int128 も）
        if constexpr (std::is_signed<T>::value) {
            neg_ = v < 0;
            set_mag(neg_ ? (UInt128)0 - (UInt128)(Int128)v : (UInt128)(Int128)v);
        } else {
            set_mag((UInt128)v);
        }
    }
    BigInt(Int128 v) : neg_(v < 0) { set_mag(v < 0 ? (UInt128)0 - (UInt128)v : (UInt128)v); }
    BigInt(UInt128 v) { set_mag(v); }
    explicit BigInt(const std::string& s) { parse(s); }
    explicit BigInt(const char* s) { parse(std::string(s)); }

    // ---- 変換 ----
    std::string to_string() const {
        if (d_.empty()) return "0";
        std::string s = neg_ ? "-" : "";
        s += std::to_string(d_.back());
        char buf[16];
        for (size_t i = d_.size() - 1; i-- > 0;) {
            std::snprintf(buf, sizeof(buf), "%09u", d_[i]);
            s += buf;
        }
        return s;
    }
    bool fits_ll() const { return cmp_mag_u128((UInt128)9223372036854775807ULL + (neg_ ? 1 : 0)) <= 0; }
    long long to_ll() const {  // long long に収まらなければ止める
        if (!fits_ll()) math_internal::fail("BigInt::to_ll", "long long に収まらない");
        const UInt128 m = mag_u128();
        return neg_ ? (long long)((UInt128)0 - m) : (long long)m;
    }
    bool fits_i128() const { return cmp_mag_u128(((UInt128)1 << 127) - (neg_ ? 0 : 1)) <= 0; }
    Int128 to_i128() const {
        if (!fits_i128()) math_internal::fail("BigInt::to_i128", "Int128 に収まらない");
        const UInt128 m = mag_u128();
        return neg_ ? (Int128)((UInt128)0 - m) : (Int128)m;
    }
    long double to_ld() const {  // おおよその値
        long double r = 0;
        for (size_t i = d_.size(); i-- > 0;) r = r * BASE + d_[i];
        return neg_ ? -r : r;
    }

    // ---- 性質 ----
    bool is_zero() const { return d_.empty(); }
    int sign() const { return d_.empty() ? 0 : (neg_ ? -1 : 1); }
    bool is_even() const { return d_.empty() || d_[0] % 2 == 0; }
    int num_digits() const {  // 10 進の桁数（0 は 1 桁。符号は数えない）
        if (d_.empty()) return 1;
        int k = 0;
        for (unsigned t = d_.back(); t; t /= 10) ++k;
        return k + BASE_DIGITS * (int)(d_.size() - 1);
    }
    long long digit_sum() const {  // 10 進の各桁の和（符号は無視）
        long long s = 0;
        for (unsigned x : d_)
            for (; x; x /= 10) s += x % 10;
        return s;
    }
    long long mod_ll(long long m) const {  // この数を m（≥ 1）で割った余り（[0, m)。負の数も）  O(桁数)
        if (m < 1) math_internal::fail("BigInt::mod_ll", "m ≥ 1 で使う");
        UInt128 r = 0;
        for (size_t i = d_.size(); i-- > 0;) r = (r * BASE + d_[i]) % (UInt128)m;
        return neg_ && r ? m - (long long)r : (long long)r;
    }

    // ---- 比較 ----
    friend int compare(const BigInt& a, const BigInt& b) {  // −1 / 0 / 1
        if (a.neg_ != b.neg_) return a.neg_ ? -1 : 1;
        const int c = cmp_mag(a.d_, b.d_);
        return a.neg_ ? -c : c;
    }
    friend bool operator==(const BigInt& a, const BigInt& b) { return a.neg_ == b.neg_ && a.d_ == b.d_; }
    friend bool operator!=(const BigInt& a, const BigInt& b) { return !(a == b); }
    friend bool operator<(const BigInt& a, const BigInt& b) { return compare(a, b) < 0; }
    friend bool operator>(const BigInt& a, const BigInt& b) { return compare(a, b) > 0; }
    friend bool operator<=(const BigInt& a, const BigInt& b) { return compare(a, b) <= 0; }
    friend bool operator>=(const BigInt& a, const BigInt& b) { return compare(a, b) >= 0; }

    // ---- 四則演算 ----
    BigInt operator-() const {
        BigInt r = *this;
        if (!r.d_.empty()) r.neg_ = !r.neg_;
        return r;
    }
    BigInt operator+() const { return *this; }
    friend BigInt operator+(const BigInt& a, const BigInt& b) {
        if (a.neg_ == b.neg_) return make(add_mag(a.d_, b.d_), a.neg_);
        const int c = cmp_mag(a.d_, b.d_);
        if (c == 0) return BigInt();
        return c > 0 ? make(sub_mag(a.d_, b.d_), a.neg_) : make(sub_mag(b.d_, a.d_), b.neg_);
    }
    friend BigInt operator-(const BigInt& a, const BigInt& b) { return a + (-b); }
    friend BigInt operator*(const BigInt& a, const BigInt& b) { return make(mul_mag(a.d_, b.d_), a.neg_ != b.neg_); }
    friend BigInt operator/(const BigInt& a, const BigInt& b) { return divmod(a, b).first; }
    friend BigInt operator%(const BigInt& a, const BigInt& b) { return divmod(a, b).second; }
    BigInt& operator+=(const BigInt& b) { return *this = *this + b; }
    BigInt& operator-=(const BigInt& b) { return *this = *this - b; }
    BigInt& operator*=(const BigInt& b) { return *this = *this * b; }
    BigInt& operator/=(const BigInt& b) { return *this = *this / b; }
    BigInt& operator%=(const BigInt& b) { return *this = *this % b; }
    BigInt& operator++() { return *this += 1; }
    BigInt& operator--() { return *this -= 1; }
    BigInt operator++(int) {
        BigInt t = *this;
        *this += 1;
        return t;
    }
    BigInt operator--(int) {
        BigInt t = *this;
        *this -= 1;
        return t;
    }
    // 商と余り（0 の方向に切り捨て）  O(桁数 × 割る数の桁数)
    friend std::pair<BigInt, BigInt> divmod(const BigInt& a, const BigInt& b) {
        if (b.d_.empty()) math_internal::fail("BigInt::/", "0 で割った");
        std::vector<unsigned> q, r;
        div_mag(a.d_, b.d_, q, r);
        return {make(std::move(q), a.neg_ != b.neg_), make(std::move(r), a.neg_)};
    }

    // ⌊√n⌋（n ≥ 0）：上の半分の桁の平方根を初めの値にしたニュートン法（割り算 3 回ほど）  O(n^2)
    friend BigInt isqrt(const BigInt& n) {
        if (n.neg_) math_internal::fail("isqrt(BigInt)", "n ≥ 0 で使う");
        if (n.d_.size() <= 4) {  // 10^36 未満は UInt128 で
            const UInt128 v = n.mag_u128();
            UInt128 r = (UInt128)std::sqrt((long double)v);
            while (r > 0 && r * r > v) --r;
            while ((r + 1) * (r + 1) <= v) ++r;
            return BigInt(r);
        }
        const size_t s = n.d_.size() / 4;  // 下の 2s 桁（10^9 進）を落とした数の平方根 + 1 を BASE^s 倍すると √n 以上
        BigInt hi;
        hi.d_.assign(n.d_.begin() + (long long)(2 * s), n.d_.end());
        BigInt x = isqrt(hi) + 1;
        x.d_.insert(x.d_.begin(), s, 0U);
        while (true) {  // 上から近づく
            BigInt y = (x + n / x) / 2;
            if (y >= x) return x;
            x = std::move(y);
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const BigInt& a) { return os << a.to_string(); }
    friend std::istream& operator>>(std::istream& is, BigInt& a) {
        std::string s;
        if (is >> s) a = BigInt(s);
        return is;
    }

  private:
    std::vector<unsigned> d_;  // 10^9 進の桁（下の桁から。上の 0 は持たない。0 は空）
    bool neg_ = false;         // 0 のときは false

    void set_mag(UInt128 m) {
        d_.clear();
        for (; m; m /= BASE) d_.push_back((unsigned)(m % BASE));
        if (d_.empty()) neg_ = false;
    }
    UInt128 mag_u128() const {
        UInt128 r = 0;
        for (size_t i = d_.size(); i-- > 0;) r = r * BASE + d_[i];
        return r;
    }
    int cmp_mag_u128(UInt128 x) const {  // |this| と x の比較
        if (d_.size() > 5) return 1;  // 10^45 > 2^128
        BigInt t((UInt128)x);
        return cmp_mag(d_, t.d_);
    }
    void parse(const std::string& s) {
        size_t i = 0;
        bool ng = false;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) ng = s[i] == '-', ++i;
        if (i == s.size()) math_internal::fail("BigInt", "数字がない");
        for (size_t j = i; j < s.size(); ++j)
            if (s[j] < '0' || s[j] > '9') math_internal::fail("BigInt", "数として読めない文字がある");
        d_.clear();
        for (size_t end = s.size(); end > i;) {
            const size_t beg = end >= i + BASE_DIGITS ? end - BASE_DIGITS : i;
            unsigned x = 0;
            for (size_t j = beg; j < end; ++j) x = x * 10 + (unsigned)(s[j] - '0');
            d_.push_back(x), end = beg;
        }
        trim(d_);
        neg_ = ng && !d_.empty();
    }
    static void trim(std::vector<unsigned>& a) {
        while (!a.empty() && a.back() == 0) a.pop_back();
    }
    static BigInt make(std::vector<unsigned> m, bool ng) {
        trim(m);
        BigInt r;
        r.neg_ = ng && !m.empty();
        r.d_ = std::move(m);
        return r;
    }
    static int cmp_mag(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {
        if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
        for (size_t i = a.size(); i-- > 0;)
            if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
        return 0;
    }
    static std::vector<unsigned> add_mag(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {
        std::vector<unsigned> r(std::max(a.size(), b.size()) + 1, 0);
        unsigned c = 0;
        for (size_t i = 0; i + 1 < r.size(); ++i) {
            unsigned s = c + (i < a.size() ? a[i] : 0) + (i < b.size() ? b[i] : 0);
            c = s >= BASE;
            r[i] = c ? s - BASE : s;
        }
        r.back() = c;
        return r;
    }
    static std::vector<unsigned> sub_mag(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {  // |a| ≥ |b|
        std::vector<unsigned> r(a.size());
        long long br = 0;
        for (size_t i = 0; i < a.size(); ++i) {
            long long t = (long long)a[i] - br - (i < b.size() ? b[i] : 0);
            br = t < 0;
            r[i] = (unsigned)(br ? t + BASE : t);
        }
        return r;
    }
    // ---- 掛け算：短い方が 40 桁（10^9 進）以下なら筆算、それより長ければ 10^6 進にしてカラツバ ----
    static std::vector<unsigned> mul_mag(const std::vector<unsigned>& a, const std::vector<unsigned>& b) {
        if (a.empty() || b.empty()) return {};
        if (std::min(a.size(), b.size()) <= 40) {
            std::vector<unsigned> r(a.size() + b.size(), 0);
            for (size_t i = 0; i < a.size(); ++i) {
                unsigned long long c = 0;
                for (size_t j = 0; j < b.size() || c; ++j) {
                    const unsigned long long cur = r[i + j] + c + (j < b.size() ? (unsigned long long)a[i] * b[j] : 0);
                    r[i + j] = (unsigned)(cur % BASE), c = cur / BASE;
                }
            }
            return r;
        }
        if (a.size() < b.size()) return mul_mag(b, a);
        if (a.size() >= 2 * b.size()) {  // 長さが大きく違う：長い方を短い方の長さで区切って掛け、ずらして足す
            std::vector<unsigned> r(a.size() + b.size() + 1, 0);
            for (size_t off = 0; off < a.size(); off += b.size()) {
                std::vector<unsigned> ch(a.begin() + (long long)off, a.begin() + (long long)std::min(off + b.size(), a.size()));
                trim(ch);
                std::vector<unsigned> pr = mul_mag(ch, b);
                trim(pr);  // カラツバの結果は後ろに 0 が付いていることがある
                unsigned c = 0;
                size_t i = 0;
                for (; i < pr.size() || c; ++i) {
                    const unsigned t = r[off + i] + c + (i < pr.size() ? pr[i] : 0);
                    c = t >= BASE;
                    r[off + i] = c ? t - BASE : t;
                }
            }
            return r;
        }
        // 10^k 進（k = 6、桁がとても多ければ 3）に直してカラツバ。途中の値が long long に収まるように k を選ぶ
        //   長さ n（2 の冪）のとき途中の値は n^2 (10^k)^2 / 32 以下：k = 6 なら n ≤ 2^14、k = 3 なら n ≤ 2^23
        int k = 6;
        size_t n = pow2_len(std::max(a.size(), b.size()), k);
        if (n > ((size_t)1 << 14)) k = 3, n = pow2_len(std::max(a.size(), b.size()), k);
        if (n > ((size_t)1 << 23)) math_internal::fail("BigInt::*", "桁が多すぎる（2×10^7 桁程度まで）");
        std::vector<long long> x = to_small(a, k), y = to_small(b, k);
        x.resize(n), y.resize(n);
        const std::vector<long long> z = karatsuba(x, y);
        const long long S = k == 6 ? 1000000 : 1000;
        std::vector<long long> ds;  // 10^k 進で繰り上げ
        long long c = 0;
        for (size_t i = 0; i < z.size() || c; ++i) {
            const long long cur = (i < z.size() ? z[i] : 0) + c;
            ds.push_back(cur % S), c = cur / S;
        }
        const size_t g = (size_t)(18 / k);  // 10^k 進 g 桁 = 10^18 未満 = 10^9 進 2 桁
        std::vector<unsigned> r;
        for (size_t i = 0; i < ds.size(); i += g) {
            unsigned long long v = 0;
            for (size_t t = g; t-- > 0;) v = v * (unsigned long long)S + (i + t < ds.size() ? (unsigned long long)ds[i + t] : 0);
            r.push_back((unsigned)(v % BASE)), r.push_back((unsigned)(v / BASE));
        }
        return r;
    }
    static size_t pow2_len(size_t limbs, int k) {  // 10^9 進 limbs 桁を 10^k 進にした長さ以上の 2 の冪
        const size_t len = (limbs + 1) / 2 * (size_t)(18 / k);
        size_t n = 1;
        while (n < len) n *= 2;
        return n;
    }
    static std::vector<long long> to_small(const std::vector<unsigned>& a, int k) {  // 10^9 進 2 桁（< 10^18）→ 10^k 進 18/k 桁
        const unsigned long long S = k == 6 ? 1000000 : 1000;
        std::vector<long long> r;
        for (size_t i = 0; i < a.size(); i += 2) {
            unsigned long long v = a[i] + (i + 1 < a.size() ? (unsigned long long)a[i + 1] * BASE : 0);
            for (int t = 0; t < 18 / k; ++t) r.push_back((long long)(v % S)), v /= S;
        }
        return r;
    }
    static std::vector<long long> karatsuba(const std::vector<long long>& a, const std::vector<long long>& b) {  // 長さは同じ 2 の冪。結果は長さ 2n
        const size_t n = a.size();
        std::vector<long long> r(2 * n, 0);
        if (n <= 32) {
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j) r[i + j] += a[i] * b[j];
            return r;
        }
        const size_t k = n / 2;
        const std::vector<long long> a1(a.begin(), a.begin() + (long long)k), a2(a.begin() + (long long)k, a.end());
        const std::vector<long long> b1(b.begin(), b.begin() + (long long)k), b2(b.begin() + (long long)k, b.end());
        const std::vector<long long> a1b1 = karatsuba(a1, b1), a2b2 = karatsuba(a2, b2);
        std::vector<long long> s1(k), s2(k);
        for (size_t i = 0; i < k; ++i) s1[i] = a1[i] + a2[i], s2[i] = b1[i] + b2[i];
        std::vector<long long> m = karatsuba(s1, s2);
        for (size_t i = 0; i < 2 * k; ++i) m[i] -= a1b1[i] + a2b2[i];
        for (size_t i = 0; i < 2 * k; ++i) r[i] += a1b1[i], r[i + k] += m[i], r[i + n] += a2b2[i];
        return r;
    }
    // ---- 割り算（Knuth のアルゴリズム D）----
    static void div_mag(const std::vector<unsigned>& a, const std::vector<unsigned>& b, std::vector<unsigned>& q, std::vector<unsigned>& r) {
        if (cmp_mag(a, b) < 0) {
            q.clear(), r = a;
            return;
        }
        if (b.size() == 1) {  // 1 桁で割る
            q.assign(a.size(), 0);
            unsigned long long rem = 0;
            for (size_t i = a.size(); i-- > 0;) {
                const unsigned long long cur = rem * BASE + a[i];
                q[i] = (unsigned)(cur / b[0]), rem = cur % b[0];
            }
            r.clear();
            if (rem) r.push_back((unsigned)rem);
            trim(q);
            return;
        }
        const size_t n = b.size(), m = a.size() - b.size();
        const unsigned long long norm = BASE / ((unsigned long long)b.back() + 1);  // 割る数の最上位を BASE/2 以上に
        std::vector<unsigned> u(a.size() + 1, 0), v(n, 0);
        {
            unsigned long long c = 0;
            for (size_t i = 0; i < a.size(); ++i) {
                const unsigned long long cur = a[i] * norm + c;
                u[i] = (unsigned)(cur % BASE), c = cur / BASE;
            }
            u[a.size()] = (unsigned)c, c = 0;
            for (size_t i = 0; i < n; ++i) {
                const unsigned long long cur = b[i] * norm + c;
                v[i] = (unsigned)(cur % BASE), c = cur / BASE;
            }
        }
        q.assign(m + 1, 0);
        for (size_t j = m + 1; j-- > 0;) {
            const unsigned long long num = (unsigned long long)u[j + n] * BASE + u[j + n - 1];
            unsigned long long qh = num / v[n - 1], rh = num % v[n - 1];
            while (qh >= BASE || qh * v[n - 2] > rh * BASE + u[j + n - 2]) {
                --qh, rh += v[n - 1];
                if (rh >= BASE) break;
            }
            // u[j..j+n] −= qh × v
            long long br = 0;
            unsigned long long c = 0;
            for (size_t i = 0; i < n; ++i) {
                const unsigned long long p = qh * v[i] + c;
                c = p / BASE;
                long long t = (long long)u[i + j] - (long long)(p % BASE) - br;
                br = t < 0;
                u[i + j] = (unsigned)(br ? t + BASE : t);
            }
            long long top = (long long)u[j + n] - (long long)c - br;
            if (top < 0) {  // 引きすぎ：1 回分戻す
                --qh;
                unsigned long long cc = 0;
                for (size_t i = 0; i < n; ++i) {
                    const unsigned long long s = (unsigned long long)u[i + j] + v[i] + cc;
                    u[i + j] = (unsigned)(s % BASE), cc = s / BASE;
                }
                top += (long long)cc;
            }
            u[j + n] = (unsigned)top;
            q[j] = (unsigned)qh;
        }
        trim(q);
        r.assign(n, 0);  // 余り = u の下 n 桁 / norm
        unsigned long long rem = 0;
        for (size_t i = n; i-- > 0;) {
            const unsigned long long cur = rem * BASE + u[i];
            r[i] = (unsigned)(cur / norm), rem = cur % norm;
        }
        trim(r);
    }
};
inline BigInt abs(const BigInt& a) { return a.sign() < 0 ? -a : a; }
inline BigInt pow(BigInt a, long long e) {  // a^e（e ≥ 0。0^0 = 1）
    if (e < 0) math_internal::fail("pow(BigInt)", "e ≥ 0 で使う");
    BigInt r = 1;
    for (; e > 0; e >>= 1) {
        if (e & 1) r *= a;
        if (e > 1) a *= a;
    }
    return r;
}
inline BigInt gcd(BigInt a, BigInt b) {  // 最大公約数（0 以上。gcd(0, 0) = 0）
    a = abs(a), b = abs(b);
    while (!b.is_zero()) {
        BigInt t = a % b;
        a = std::move(b), b = std::move(t);
    }
    return a;
}
template <class T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
inline BigInt gcd(const BigInt& a, T b) { return gcd(a, BigInt(b)); }
template <class T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
inline BigInt gcd(T a, const BigInt& b) { return gcd(BigInt(a), b); }
inline BigInt lcm(const BigInt& a, const BigInt& b) {  // 最小公倍数（0 以上。どちらかが 0 なら 0）
    if (a.is_zero() || b.is_zero()) return BigInt(0);
    return abs(a) / gcd(a, b) * abs(b);
}
template <class T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
inline BigInt lcm(const BigInt& a, T b) { return lcm(a, BigInt(b)); }
template <class T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
inline BigInt lcm(T a, const BigInt& b) { return lcm(BigInt(a), b); }
#endif  // GG_BIGINT
