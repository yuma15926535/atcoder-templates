#ifndef GG_PERMUTATION
#define GG_PERMUTATION
// ---- 順列・組合せの番号（数学第二世代）  辞書順で何番目か・k 番目の並び（番号は 0 から）・転倒数・置換の偶奇 ----
// 列は vector でも string でもよい。厳密な番号は Int128（収まらなければ止める。そのときは _mod）
namespace math_internal {
struct PermFenwick {  // 0/1 の出現数を数える BIT
    explicit PermFenwick(int n) : n_(n), t_((size_t)n + 1, 0) {}
    void add(int i, int v) {
        for (++i; i <= n_; i += i & -i) t_[(size_t)i] += v;
    }
    int prefix(int i) const {  // [0, i) の和
        int s = 0;
        for (; i > 0; i -= i & -i) s += t_[(size_t)i];
        return s;
    }
    int kth(int k) const {  // 和が k + 1 になる最小の位置（k 番目の 1。0 から）
        int pos = 0, lg = 1;
        while (lg * 2 <= n_) lg *= 2;
        for (; lg > 0; lg /= 2)
            if (pos + lg <= n_ && t_[(size_t)(pos + lg)] <= k) pos += lg, k -= t_[(size_t)pos];
        return pos;
    }

  private:
    int n_;
    std::vector<int> t_;
};
inline Int128 gcd128(Int128 a, Int128 b) {  // a, b ≥ 0
    while (b) {
        const Int128 t = a % b;
        a = b, b = t;
    }
    return a;
}
const Int128 PERM_INF = (Int128)(((UInt128)1 << 127) - 1);  // Int128 の最大値（「とても大きい」の印にも使う）
// M × a / b（割り切れると分かっているとき）。PERM_INF を超えたら PERM_INF
inline Int128 mul_div_exact(Int128 M, Int128 a, Int128 b) {
    const Int128 g = gcd128(M, b);
    M /= g, b /= g;
    a /= b;  // gcd(M/g, b/g) = 1 なので b/g は a を割り切る
    if (a != 0 && M > PERM_INF / a) return PERM_INF;
    return M * a;
}
inline Int128 binom128(long long n, long long r) {  // min(nCr, PERM_INF)（PERM_INF は nCr になれないので「超えた」の印）
    if (r < 0 || r > n) return 0;
    r = std::min(r, n - r);
    Int128 v = 1;
    for (long long i = 1; i <= r; ++i) {  // v = C(n − r + i, i)：増えていくので PERM_INF になったら終わり
        v = mul_div_exact(v, n - r + i, i);
        if (v == PERM_INF) return v;
    }
    return v;
}
template <class C>
std::vector<int> perm_compress(const C& a, int& kinds) {  // 値を 0, 1, … に（同じ値は同じ番号）
    using T = typename std::decay<decltype(a[0])>::type;
    const size_t n = a.size();
    std::vector<int> r(n);
    kinds = 0;
    if (n == 0) return r;
    if constexpr (std::is_integral<T>::value) {  // 値の幅が狭い（1〜n など）なら数えるだけ  O(N)
        const auto mm = std::minmax_element(a.begin(), a.end());
        const T lo = *mm.first, hi = *mm.second;
        if ((UInt128)((Int128)hi - (Int128)lo) < (UInt128)2 * n + 64) {
            std::vector<int> id((size_t)((Int128)hi - (Int128)lo) + 1, 0);
            for (size_t i = 0; i < n; ++i) id[(size_t)((Int128)a[i] - (Int128)lo)] = 1;
            for (int& x : id) x = x ? kinds++ : 0;
            for (size_t i = 0; i < n; ++i) r[i] = id[(size_t)((Int128)a[i] - (Int128)lo)];
            return r;
        }
    }
    std::vector<std::pair<T, int>> v(n);
    for (size_t i = 0; i < n; ++i) v[i] = {a[i], (int)i};
    std::sort(v.begin(), v.end());
    for (size_t i = 0; i < n; ++i) {
        if (i && v[i - 1].first < v[i].first) ++kinds;
        r[(size_t)v[i].second] = kinds;
    }
    ++kinds;
    return r;
}
}  // namespace math_internal

// 転倒数：i < j かつ a_i > a_j の組の数（値は比べられれば何でも。同じ値は数えない）  O(N log N)
template <class C>
long long inversion_count(const C& a) {
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(a, k);
    math_internal::PermFenwick bit(k);
    long long s = 0;
    for (size_t i = 0; i < c.size(); ++i) {
        s += (long long)i - bit.prefix(c[i] + 1);
        bit.add(c[i], 1);
    }
    return s;
}
// 置換の偶奇：値がすべて異なる列 a をソートする置換が偶置換なら 0、奇置換なら 1（= 転倒数の偶奇 = (N − サイクル数) mod 2）  O(N log N)（値が 0〜N 程度なら O(N)）
template <class C>
int perm_parity(const C& a) {
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(a, k);
    if (k != (int)c.size()) math_internal::fail("perm_parity", "値がすべて異なる列で使う（同じ値があると偶奇が決まらない）");
    std::vector<char> seen(c.size(), 0);
    int cycles = 0;
    for (size_t i = 0; i < c.size(); ++i)
        if (!seen[i]) {
            ++cycles;
            for (size_t j = i; !seen[j]; j = (size_t)c[j]) seen[j] = 1;
        }
    return ((int)c.size() - cycles) % 2;
}

// ---- 順列（値がすべて異なる列）----
// p が「同じ値の並べ替え全体」の辞書順で何番目か（0 から）。値は何でもよい（1〜n でも 0〜n−1 でも）。n ≤ 33  O(N log N)
template <class C>
Int128 perm_rank(const C& p) {
    const int n = (int)p.size();
    if (n > 33) math_internal::fail("perm_rank", "n ≤ 33 で使う（番号が Int128 に収まらない。大きいときは perm_rank_mod）");
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(p, k);
    if (k != n) math_internal::fail("perm_rank", "値がすべて異なる列で使う（同じ値があるなら multiset_perm_rank）");
    math_internal::PermFenwick bit(n);
    Int128 r = 0;
    for (int i = 0; i < n; ++i) {
        const int less = c[(size_t)i] - bit.prefix(c[(size_t)i]);  // まだ使っていない、より小さい値の数
        r = r * (n - i) + less;                                     // ホーナー法：Σ less_i (n−1−i)!
        bit.add(c[(size_t)i], 1);
    }
    return r;
}
template <class mint, class C>
mint perm_rank_mod(const C& p) {  // perm_rank を mod で（n はいくつでも）  O(N log N)
    const int n = (int)p.size();
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(p, k);
    if (k != n) math_internal::fail("perm_rank_mod", "値がすべて異なる列で使う（同じ値があるなら multiset_perm_rank_mod）");
    math_internal::PermFenwick bit(n);
    mint r = 0;
    for (int i = 0; i < n; ++i) {
        r = r * mint(n - i) + mint(c[(size_t)i] - bit.prefix(c[(size_t)i]));
        bit.add(c[(size_t)i], 1);
    }
    return r;
}
// 0〜n−1 の順列で辞書順 k 番目（0 から）のもの（0 ≤ k < n!）  O(N log N)
inline std::vector<int> perm_unrank(int n, Int128 k) {
    if (n < 0 || k < 0) math_internal::fail("perm_unrank", "n ≥ 0、k ≥ 0 で使う");
    std::vector<Int128> f((size_t)n + 1, 1);  // f[m] = m!（34 以上は PERM_INF）
    for (int m = 1; m <= n; ++m) f[(size_t)m] = m >= 34 ? math_internal::PERM_INF : f[(size_t)m - 1] * m;
    if (n <= 33 && k >= f[(size_t)n]) math_internal::fail("perm_unrank", "k < n! で使う");
    math_internal::PermFenwick bit(n);
    for (int i = 0; i < n; ++i) bit.add(i, 1);
    std::vector<int> p((size_t)n);
    for (int i = 0; i < n; ++i) {
        const Int128 fm = f[(size_t)(n - 1 - i)];
        const int d = fm == math_internal::PERM_INF ? 0 : (int)(k / fm);
        if (fm != math_internal::PERM_INF) k %= fm;
        p[(size_t)i] = bit.kth(d), bit.add(p[(size_t)i], -1);
    }
    return p;
}

// ---- 重複のある並び（同じ値があってよい。"aabbc" の並べ替えなど）----
namespace math_internal {
inline Int128 multiset_count(const std::vector<int>& cnt, const char* func) {  // 並べ方の数 (Σ c)! / Π c!（Int128 に収まらなければ止める）
    Int128 M = 1;
    long long s = 0;
    for (int c : cnt)
        for (int j = 1; j <= c; ++j) {
            ++s;
            M = mul_div_exact(M, s, j);  // 1 個ずつ足す：M × s / (その値の個数)
            if (M == PERM_INF) fail(func, "並べ方の数が Int128 に収まらない（_mod を使う）");
        }
    return M;
}
}  // namespace math_internal
// a が「a の並べ替え（同じものは区別しない）」全体の辞書順で何番目か（0 から）  O(N log N)
template <class C>
Int128 multiset_perm_rank(const C& a) {
    const int n = (int)a.size();
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(a, k);
    std::vector<int> cnt((size_t)k, 0);
    for (int x : c) ++cnt[(size_t)x];
    Int128 M = math_internal::multiset_count(cnt, "multiset_perm_rank");  // 残りの並べ方の数
    math_internal::PermFenwick bit(k);
    for (int v = 0; v < k; ++v) bit.add(v, cnt[(size_t)v]);
    Int128 r = 0;
    for (int i = 0; i < n; ++i) {
        const long long s = n - i;
        const int u = c[(size_t)i], less = bit.prefix(u);  // 先頭に u より小さい値を置く並べ方：M × less / s
        r += math_internal::mul_div_exact(M, less, s);
        M = math_internal::mul_div_exact(M, cnt[(size_t)u], s);
        --cnt[(size_t)u], bit.add(u, -1);
    }
    return r;
}
template <class mint, class C>
mint multiset_perm_rank_mod(const C& a) {  // multiset_perm_rank を mod で（mod は素数で n より大きい）  O(N log N)
    const int n = (int)a.size();
    int k = 0;
    const std::vector<int> c = math_internal::perm_compress(a, k);
    std::vector<int> cnt((size_t)k, 0);
    for (int x : c) ++cnt[(size_t)x];
    Binomial<mint> B(n);
    mint M = B.fact(n);
    for (int x : cnt) M *= B.inv_fact(x);
    math_internal::PermFenwick bit(k);
    for (int v = 0; v < k; ++v) bit.add(v, cnt[(size_t)v]);
    mint r = 0;
    for (int i = 0; i < n; ++i) {
        const int u = c[(size_t)i];
        const mint is = B.inv(n - i);
        r += M * mint(bit.prefix(u)) * is;
        M *= mint(cnt[(size_t)u]) * is;
        --cnt[(size_t)u], bit.add(u, -1);
    }
    return r;
}
// a の並べ替え（同じものは区別しない）で辞書順 k 番目（0 から）のもの（a の順番は何でもよい）  O(N × 値の種類)
template <class C>
C multiset_perm_unrank(C a, Int128 k) {
    std::sort(a.begin(), a.end());
    const int n = (int)a.size();
    int kinds = 0;
    const std::vector<int> c = math_internal::perm_compress(a, kinds);
    std::vector<int> cnt((size_t)kinds, 0), first((size_t)kinds, 0);
    for (int i = n - 1; i >= 0; --i) ++cnt[(size_t)c[(size_t)i]], first[(size_t)c[(size_t)i]] = i;
    Int128 M = math_internal::multiset_count(cnt, "multiset_perm_unrank");
    if (k < 0 || k >= M) math_internal::fail("multiset_perm_unrank", "0 ≤ k < 並べ方の数 で使う");
    C res = a;
    for (int i = 0; i < n; ++i) {
        const long long s = n - i;
        for (int v = 0; v < kinds; ++v) {
            if (cnt[(size_t)v] == 0) continue;
            const Int128 blk = math_internal::mul_div_exact(M, cnt[(size_t)v], s);  // 先頭が v の並べ方
            if (k < blk) {
                res[(size_t)i] = a[(size_t)first[(size_t)v]], M = blk, --cnt[(size_t)v];
                break;
            }
            k -= blk;
        }
    }
    return res;
}

// ---- 組合せ（{0, …, n−1} から r 個。小さい順に並べた列の辞書順）----
// k 番目（0 から）の r 個の組（小さい順）。0 ≤ k < C(n, r)  O(n × min(r, 130))
inline std::vector<int> kth_combination(int n, int r, Int128 k) {
    if (n < 0 || r < 0 || r > n) math_internal::fail("kth_combination", "0 ≤ r ≤ n で使う");
    const Int128 tot = math_internal::binom128(n, r);
    if (k < 0 || (tot != math_internal::PERM_INF && k >= tot)) math_internal::fail("kth_combination", "0 ≤ k < C(n, r) で使う");
    std::vector<int> res;
    for (int x = 0; x < n && (int)res.size() < r; ++x) {
        const int need = r - (int)res.size();
        const Int128 c = math_internal::binom128(n - 1 - x, need - 1);  // x を選ぶ組の数（PERM_INF は「k より大きい」）
        if (c == math_internal::PERM_INF || k < c) res.push_back(x);
        else k -= c;
    }
    return res;
}
// c（{0, …, n−1} の異なる r 個。順番は何でもよい）が辞書順で何番目か（0 から。C(n, r) が Int128 に収まること）  O(n + r log r)
template <class C>
Int128 combination_rank(int n, const C& c) {
    std::vector<long long> v(c.begin(), c.end());
    std::sort(v.begin(), v.end());
    const int r = (int)v.size();
    for (int i = 0; i < r; ++i)
        if (v[(size_t)i] < 0 || v[(size_t)i] >= n || (i && v[(size_t)i] == v[(size_t)i - 1])) math_internal::fail("combination_rank", "0〜n−1 の異なる値で使う");
    if (math_internal::binom128(n, r) == math_internal::PERM_INF) math_internal::fail("combination_rank", "C(n, r) が Int128 に収まらない");
    Int128 res = 0;
    long long prev = -1;
    for (int i = 0; i < r; ++i) {  // prev+1 〜 v_i − 1 を i 番目に置く組：Σ_x C(n−1−x, r−i−1) = C(n−prev−1, r−i) − C(n−v_i, r−i)
        res += math_internal::binom128(n - prev - 1, r - i) - math_internal::binom128(n - v[(size_t)i], r - i);
        prev = v[(size_t)i];
    }
    return res;
}
#endif  // GG_PERMUTATION
