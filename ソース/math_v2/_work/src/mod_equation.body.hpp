#ifndef GG_MOD_EQUATION
#define GG_MOD_EQUATION
// ---- 合同式・一次不定方程式・中国剰余定理（数学第二世代） ----
struct ModSolution {  // 解は x ≡ r (mod m)（0 ≤ r < m）。ok が false なら解なし
    bool ok = false;
    long long r = 0, m = 1;
    bool overflow = false;  // crt で lcm が long long を超えた。そのとき ok は「long long に入る解があるか」、r はその最小。m は 0（使わない）
};
struct Diophantine {  // a x + b y = c の整数解すべて：x = x0 + k·dx、y = y0 − k·dy（k は整数）。ok が false なら解なし
    bool ok = false;
    long long x0 = 0, y0 = 0, dx = 0, dy = 0;  // dx ≥ 0。b ≠ 0 なら 0 ≤ x0 < dx（x0 は 0 以上で最小の x）
    long long x(long long k) const { return check(x0 + (Int128)k * dx); }  // k 番目の解の x（収まらなければ止める）
    long long y(long long k) const { return check(y0 - (Int128)k * dy); }

  private:
    static long long check(Int128 v) {
        if (v > (Int128)9223372036854775807LL || v < -(Int128)9223372036854775807LL - 1) math_internal::fail("Diophantine", "long long に収まらない");
        return (long long)v;
    }
};
struct DiophantinePoint {  // a x + b y = c の 1 つの解（ok が false なら条件を満たす解なし）
    bool ok = false;
    long long x = 0, y = 0;
};

namespace math_internal {
inline bool fits_ll(Int128 v) { return -(Int128)9223372036854775807LL - 1 <= v && v <= (Int128)9223372036854775807LL; }
// x ≡ r1 (mod m1) と x ≡ r2 (mod m2)（0 ≤ r1 < m1 ≤ 2^63 − 1、0 ≤ r2 < m2 ≤ 2^63 − 1）を合わせる。r < m = lcm < 2^126
struct Crt128 {
    bool ok;
    Int128 r, m;
};
inline Crt128 crt_merge128(long long r1, long long m1, long long r2, long long m2) {
    if (m1 < m2) std::swap(r1, r2), std::swap(m1, m2);
    if (m1 % m2 == 0) return Crt128{r1 % m2 == r2, r1, m1};  // m2 | m1（よくある形を速く）
    const ExtGcd e = ext_gcd(m1, m2);  // m1·e.x + m2·e.y = g なので (m1 / g)·e.x ≡ 1 (mod m2 / g)
    const long long g = e.g, diff = r2 - r1;  // |diff| < 2^63
    if (diff % g != 0) return Crt128{false, 0, 1};
    const long long mm = m2 / g;  // r1 + m1·t ≡ r2 (mod m2) ⇔ (m1 / g)·t ≡ diff / g (mod mm)
    const long long t = (long long)((Int128)mod_floor(diff / g, mm) * mod_floor(e.x, mm) % mm);
    return Crt128{true, (Int128)r1 + (Int128)m1 * t, (Int128)m1 * mm};
}
// 以下 crt_all の本体：(r, m) に (r2, m2) を足していく。lcm が long long を超えたら「r がそのまま解か」だけを見る
struct CrtState {
    bool ok = true, overflow = false;
    long long r = 0, m = 1;
    void add(long long r2, long long m2, const char* func) {
        if (m2 < 1) fail(func, "mod は 1 以上で使う");
        if (!ok) return;
        r2 = mod_floor(r2, m2);
        if (overflow) {  // 解は r + t·(lcm)（lcm > 2^63）。r が合わなければ次の解はもう long long に入らない
            if (r % m2 != r2) ok = false;
            return;
        }
        const Crt128 c = crt_merge128(r, m, r2, m2);
        if (!c.ok) { ok = false; return; }
        if (c.m > (Int128)9223372036854775807LL) {
            overflow = true;
            if (c.r > (Int128)9223372036854775807LL) { ok = false; return; }
            r = (long long)c.r, m = 0;
            return;
        }
        r = (long long)c.r, m = (long long)c.m;
    }
    ModSolution result() const {
        if (!ok) return ModSolution{false, 0, overflow ? 0 : 1, overflow};
        return ModSolution{true, r, overflow ? 0 : m, overflow};
    }
};
}  // namespace math_internal

// ---- 一次合同式 a x ≡ b (mod m) ----
inline ModSolution linear_congruence(long long a, long long b, long long m) {  // 解を x ≡ r (mod m / g) の形で（g = gcd(a, m)。a, b は負でも）
    if (m < 1) math_internal::fail("linear_congruence", "m ≥ 1 で使う");
    a = mod_floor(a, m), b = mod_floor(b, m);
    const ExtGcd e = ext_gcd(a, m);  // a·e.x + m·e.y = g なので (a / g)·e.x ≡ 1 (mod m / g)。a = 0 なら g = m
    const long long g = e.g;
    if (b % g != 0) return ModSolution{false, 0, 1, false};
    const long long m2 = m / g;
    return ModSolution{true, (long long)((Int128)(b / g) * mod_floor(e.x, m2) % m2), m2, false};
}

// ---- 一次不定方程式 a x + b y = c ----
namespace math_internal {
struct Dio128 {  // Diophantine と同じ形（128 ビット。y0 が long long に入らないときにも使える）
    bool ok = false;
    Int128 x0 = 0, y0 = 0, dx = 0, dy = 0;
};
inline Dio128 dio_solve(long long a, long long b, long long c, const char* func) {
    const long long LMIN = -9223372036854775807LL - 1;
    if (a == LMIN || b == LMIN) fail(func, "a, b に LLONG_MIN は使えない");
    Dio128 d;
    if (a == 0 && b == 0) {
        d.ok = c == 0;
        return d;
    }
    if (b == 0) {
        if ((Int128)c % a != 0) return d;
        d.ok = true, d.x0 = (Int128)c / a, d.y0 = 0, d.dx = 0, d.dy = -1;
        return d;
    }
    const ExtGcd e = ext_gcd(a, b);  // a·e.x + b·e.y = g
    const long long g = e.g;
    if (c % g != 0) return d;
    d.ok = true;
    d.dx = (b < 0 ? -(Int128)b : (Int128)b) / g;
    d.dy = (Int128)(a / g) * (b < 0 ? -1 : 1);
    d.x0 = (Int128)e.x * (c / g) % d.dx;  // |e.x| ≤ |b| / g なので収まる
    if (d.x0 < 0) d.x0 += d.dx;
    d.y0 = ((Int128)c - (Int128)a * d.x0) / b;  // 割り切れる
    return d;
}
}  // namespace math_internal
// 解を x = x0 + k·dx、y = y0 − k·dy で（g = gcd(a, b)、dx = |b| / g、dy = a / g × (b の符号)）。0 ≤ x0 < dx
//   b = 0 なら x は 1 つに決まり y は何でもよい（dx = 0、dy = −1、y0 = 0）。a = 0 なら x が何でもよい（dx = 1、dy = 0）
//   a = b = 0 なら c = 0 のとき x, y は何でもよい（ok = true、dx = dy = 0）。a, b に LLONG_MIN は使えない
inline Diophantine linear_diophantine(long long a, long long b, long long c) {
    const math_internal::Dio128 d = math_internal::dio_solve(a, b, c, "linear_diophantine");
    if (!math_internal::fits_ll(d.x0) || !math_internal::fits_ll(d.y0)) math_internal::fail("linear_diophantine", "x0 / y0 が long long に収まらない");
    return Diophantine{d.ok, (long long)d.x0, (long long)d.y0, (long long)d.dx, (long long)d.dy};
}
namespace math_internal {
// a x + b y = c、xl ≤ x < xr、yl ≤ y < yr の解の数（2^100 を超えたら 2^100）
inline Int128 dio_count(long long a, long long b, long long c, Int128 xl, Int128 xr, Int128 yl, Int128 yr, const char* func) {
    const Dio128 d = dio_solve(a, b, c, func);
    if (xl >= xr || yl >= yr || !d.ok) return 0;
    const Int128 INF = (Int128)1 << 100;
    if (d.dx == 0 && d.dy == 0) {  // a = b = 0, c = 0：どの組も解
        const Int128 w = xr - xl, h = yr - yl;
        return w > INF / h ? INF : w * h;
    }
    Int128 lo = -INF, hi = INF;  // k の範囲 [lo, hi)
    if (d.dx == 0) {
        if (!(xl <= d.x0 && d.x0 < xr)) return 0;
    } else {  // xl ≤ x0 + k dx < xr
        lo = std::max(lo, ceil_div128(xl - d.x0, d.dx));
        hi = std::min(hi, ceil_div128(xr - d.x0, d.dx));
    }
    if (d.dy == 0) {
        if (!(yl <= d.y0 && d.y0 < yr)) return 0;
    } else if (d.dy > 0) {  // yl ≤ y0 − k dy < yr ⇔ (y0 − yr) / dy < k ≤ (y0 − yl) / dy
        lo = std::max(lo, floor_div128(d.y0 - yr, d.dy) + 1);
        hi = std::min(hi, floor_div128(d.y0 - yl, d.dy) + 1);
    } else {  // y = y0 + k e（e = −dy > 0）
        lo = std::max(lo, ceil_div128(yl - d.y0, -d.dy));
        hi = std::min(hi, ceil_div128(yr - d.y0, -d.dy));
    }
    return hi > lo ? hi - lo : 0;
}
}  // namespace math_internal
inline long long diophantine_count(long long a, long long b, long long c, long long xl, long long xr, long long yl, long long yr) {  // xl ≤ x < xr、yl ≤ y < yr の解の数
    const Int128 n = math_internal::dio_count(a, b, c, xl, xr, yl, yr, "diophantine_count");
    if (n > (Int128)9223372036854775807LL) math_internal::fail("diophantine_count", "個数が long long に収まらない");
    return (long long)n;
}
inline long long diophantine_count_closed(long long a, long long b, long long c, long long XL, long long XR, long long YL, long long YR) {  // XL ≤ x ≤ XR、YL ≤ y ≤ YR
    const Int128 n = math_internal::dio_count(a, b, c, XL, (Int128)XR + 1, YL, (Int128)YR + 1, "diophantine_count_closed");
    if (n > (Int128)9223372036854775807LL) math_internal::fail("diophantine_count_closed", "個数が long long に収まらない");
    return (long long)n;
}
// x ≥ x_min で x が最小の解（b = 0 なら y = 0 を返す。x が long long に入らなければ ok = false。y が入らなければ止める）
inline DiophantinePoint diophantine_min_x(long long a, long long b, long long c, long long x_min) {
    const math_internal::Dio128 d = math_internal::dio_solve(a, b, c, "diophantine_min_x");
    if (!d.ok) return DiophantinePoint{};
    Int128 x = x_min, y = 0;  // a = b = c = 0 ならこのまま
    if (d.dx == 0 && d.dy != 0) {  // b = 0：x は 1 つ
        if (d.x0 < x_min) return DiophantinePoint{};
        x = d.x0;
    } else if (d.dx != 0) {
        const Int128 k = math_internal::ceil_div128((Int128)x_min - d.x0, d.dx);
        x = d.x0 + k * d.dx, y = d.y0 - k * d.dy;
    }
    if (x > (Int128)9223372036854775807LL) return DiophantinePoint{};  // long long に入る解がない
    if (!math_internal::fits_ll(y)) math_internal::fail("diophantine_min_x", "y が long long に収まらない");
    return DiophantinePoint{true, (long long)x, (long long)y};
}

// ---- 中国剰余定理（mod は互いに素でなくてよい。答えは 0 ≤ r < lcm で最小） ----
inline ModSolution crt_merge(long long r1, long long m1, long long r2, long long m2) {  // x ≡ r1 (mod m1) かつ x ≡ r2 (mod m2)
    math_internal::CrtState s;
    s.add(r1, m1, "crt_merge");
    s.add(r2, m2, "crt_merge");
    return s.result();
}
inline ModSolution crt_all(const std::vector<long long>& rs, const std::vector<long long>& ms) {  // x ≡ rs[i] (mod ms[i]) をすべて（空なら x ≡ 0 (mod 1)）
    if (rs.size() != ms.size()) math_internal::fail("crt_all", "rs と ms の長さが違う");
    math_internal::CrtState s;
    for (size_t i = 0; i < rs.size(); ++i) s.add(rs[i], ms[i], "crt_all");
    return s.result();
}
// x ≡ rs[i] (mod ms[i])（ms は互いに素）を満たす 0 ≤ x < Π ms の x を mod で（Π ms が long long を超えてよい）  O(n^2)
inline long long garner(const std::vector<long long>& rs, const std::vector<long long>& ms, long long mod) {
    if (rs.size() != ms.size()) math_internal::fail("garner", "rs と ms の長さが違う");
    if (mod < 1) math_internal::fail("garner", "mod ≥ 1 で使う");
    for (long long m : ms)
        if (m < 1) math_internal::fail("garner", "mod は 1 以上で使う");
    const size_t n = ms.size();
    std::vector<long long> coef(n + 1, 1), cons(n + 1, 0), mm(ms);
    mm.push_back(mod);
    for (size_t i = 0; i < n; ++i) {
        const Maybe<long long> inv = try_mod_inv(coef[i], ms[i]);
        if (!inv.ok) math_internal::fail("garner", "ms が互いに素でない");
        const long long t = mod_mul(mod_floor(rs[i], ms[i]) - cons[i], inv.value, ms[i]);
        for (size_t j = i + 1; j <= n; ++j) {
            cons[j] = (long long)(((Int128)cons[j] + (Int128)t * coef[j]) % mm[j]);
            coef[j] = (long long)((Int128)coef[j] * ms[i] % mm[j]);
        }
    }
    return cons[n];
}
#endif  // GG_MOD_EQUATION
