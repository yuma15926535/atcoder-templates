// =============================================================================
//  ① Segment Tree  (1点更新 / 区間モノイド積)   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  値の型      : long long 固定
//  区間の流儀  : 半開区間 [l, r)  / 0-indexed
//  依存        : 標準ライブラリのみ (C++17)。bits/stdc++.h も使わないので MSVC でも通る。
//
//  計算量
//    構築            O(n)
//    point_set       O(log n)
//    point_get       O(1)
//    range_query     O(log n)
//    max_right/left  O(log n)
//
//  モノイドの切り替え
//    下の「★ 使用するモノイドをここで選ぶ ★」ブロックで、
//    有効にしたい 1 行のコメントを外し、他を // でコメントアウトするだけ。
// =============================================================================

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

using ll = long long;

// =============================================================================
//  モノイド定義
//    S      : 扱う値の型 (long long 固定)
//    e()    : 単位元
//    op(a,b): 二項演算 (結合則を満たすこと)
// =============================================================================

// 区間和
struct MonoidSum {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return a + b; }
};

// 区間最小値
struct MonoidMin {
    using S = ll;
    static S e() { return std::numeric_limits<ll>::max(); }
    static S op(S a, S b) { return a < b ? a : b; }
};

// 区間最大値
struct MonoidMax {
    using S = ll;
    static S e() { return std::numeric_limits<ll>::min(); }
    static S op(S a, S b) { return a > b ? a : b; }
};

// 区間 GCD  (注: 負数を渡すと絶対値の gcd になる。LLONG_MIN は使用不可)
struct MonoidGcd {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return std::gcd(a, b); }
};

// 区間 LCM  (注: 非負整数を想定。0 は吸収元。オーバーフローに注意)
struct MonoidLcm {
    using S = ll;
    static S e() { return 1; }
    static S op(S a, S b) {
        if (a == 0 || b == 0) return 0;
        S g = std::gcd(a, b);
        return (a / g) * b;  // a*b/g だとオーバーフローしやすいので先に割る
    }
};

// 区間 XOR
struct MonoidXor {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return a ^ b; }
};

// 区間 AND
struct MonoidAnd {
    using S = ll;
    static S e() { return ~0LL; }  // 全ビット 1
    static S op(S a, S b) { return a & b; }
};

// 区間 OR
struct MonoidOr {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return a | b; }
};

// =============================================================================
//  Segment Tree 本体 (非再帰)
// =============================================================================
template <class M>
struct SegTree {
    using S = typename M::S;

    SegTree() : SegTree(0) {}
    explicit SegTree(int n) : SegTree(std::vector<S>(n, M::e())) {}
    explicit SegTree(const std::vector<S>& v) : _n((int)v.size()) {
        _log = 0;
        while ((1 << _log) < _n) ++_log;
        _size = 1 << _log;
        d.assign(2 * (size_t)_size, M::e());
        for (int i = 0; i < _n; i++) d[_size + i] = v[i];
        for (int i = _size - 1; i >= 1; i--) pull(i);
    }

    int size() const { return _n; }

    // p 番目を x に書き換える
    void point_set(int p, S x) {
        assert(0 <= p && p < _n);
        p += _size;
        d[p] = x;
        for (int i = 1; i <= _log; i++) pull(p >> i);
    }

    // p 番目の値を取得
    S point_get(int p) const {
        assert(0 <= p && p < _n);
        return d[p + _size];
    }

    // [l, r) のモノイド積
    S range_query(int l, int r) const {
        assert(0 <= l && l <= r && r <= _n);
        S sml = M::e(), smr = M::e();
        l += _size;
        r += _size;
        while (l < r) {
            if (l & 1) sml = M::op(sml, d[l++]);
            if (r & 1) smr = M::op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }
        return M::op(sml, smr);
    }

    // 全体のモノイド積 (= range_query(0, n) だが O(1))
    S all_query() const { return d[1]; }

    // f(range_query(l, r)) が true となる最大の r を返す。
    //   f(e()) == true が前提。f は「単調」(true が続いた後 false になる) であること。
    //   例: l から右へ和を伸ばして総和が X 以下に収まる最大の r
    template <class F>
    int max_right(int l, F f) const {
        assert(0 <= l && l <= _n);
        assert(f(M::e()));
        if (l == _n) return _n;
        l += _size;
        S sm = M::e();
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(M::op(sm, d[l]))) {
                while (l < _size) {
                    l = 2 * l;
                    if (f(M::op(sm, d[l]))) {
                        sm = M::op(sm, d[l]);
                        l++;
                    }
                }
                return l - _size;
            }
            sm = M::op(sm, d[l]);
            l++;
        } while ((l & -l) != l);
        return _n;
    }

    // f(range_query(l, r)) が true となる最小の l を返す (= 左へ最大限伸ばす)。
    //   f(e()) == true が前提。ACL の min_left と同じもの。
    template <class F>
    int max_left(int r, F f) const {
        assert(0 <= r && r <= _n);
        assert(f(M::e()));
        if (r == 0) return 0;
        r += _size;
        S sm = M::e();
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(M::op(d[r], sm))) {
                while (r < _size) {
                    r = 2 * r + 1;
                    if (f(M::op(d[r], sm))) {
                        sm = M::op(d[r], sm);
                        r--;
                    }
                }
                return r + 1 - _size;
            }
            sm = M::op(d[r], sm);
        } while ((r & -r) != r);
        return 0;
    }

   private:
    int _n = 0, _size = 1, _log = 0;
    std::vector<S> d;
    void pull(int k) { d[k] = M::op(d[2 * k], d[2 * k + 1]); }
};

// =============================================================================
//  ★ 使用するモノイドをここで選ぶ ★  (1 つだけ有効にする)
// =============================================================================
using M = MonoidSum;
// using M = MonoidMin;
// using M = MonoidMax;
// using M = MonoidGcd;
// using M = MonoidLcm;
// using M = MonoidXor;
// using M = MonoidAnd;
// using M = MonoidOr;

using Seg = SegTree<M>;

// =============================================================================
//  ここから下は「ランダムテスト」と「使用例デモ」。提出時は丸ごと消してよい。
// =============================================================================

namespace tester {

// 非可換モノイド (左右の順序が保たれているかの検証用): 「最後の非ゼロ値」
struct MonoidLastNonZero {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return b != 0 ? b : a; }
};

static ll rnd_range(std::mt19937_64& rng, ll lo, ll hi) {  // [lo, hi]
    return lo + (ll)(rng() % (unsigned long long)(hi - lo + 1));
}

// LCM テスト用: 720720 (= 2^4*3^2*5*7*11*13) の約数だけを使えば lcm は 720720 以下で収まる
static std::vector<ll> divisors_of_720720() {
    std::vector<ll> d;
    for (ll x = 1; x <= 720720; x++)
        if (720720 % x == 0) d.push_back(x);
    return d;
}

// 素朴な区間積
template <class Mo>
static typename Mo::S naive_prod(const std::vector<ll>& a, int l, int r) {
    typename Mo::S s = Mo::e();
    for (int i = l; i < r; i++) s = Mo::op(s, a[i]);
    return s;
}

// 1 つのモノイドについて、point_set / point_get / range_query / all_query を検証
template <class Mo>
static void test_one_monoid(const char* name, std::mt19937_64& rng, int n, int q, ll lo, ll hi,
                            const std::vector<ll>& pool = std::vector<ll>()) {
    auto gen = [&]() -> ll { return pool.empty() ? rnd_range(rng, lo, hi) : pool[rng() % pool.size()]; };
    std::vector<ll> a(n);
    for (int i = 0; i < n; i++) a[i] = gen();
    SegTree<Mo> seg(a);

    // 構築直後の全区間チェック
    for (int l = 0; l <= n; l++)
        for (int r = l; r <= n; r++)
            assert(seg.range_query(l, r) == naive_prod<Mo>(a, l, r));

    for (int t = 0; t < q; t++) {
        int type = (int)(rng() % 4);
        if (type == 0 && n > 0) {  // point_set
            int p = (int)(rng() % n);
            ll x = gen();
            a[p] = x;
            seg.point_set(p, x);
        } else if (type == 1 && n > 0) {  // point_get
            int p = (int)(rng() % n);
            assert(seg.point_get(p) == a[p]);
        } else if (type == 2) {  // range_query
            int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
            if (l > r) std::swap(l, r);
            assert(seg.range_query(l, r) == naive_prod<Mo>(a, l, r));
        } else {  // all_query
            assert(seg.all_query() == naive_prod<Mo>(a, 0, n));
        }
    }
    (void)name;
}

// max_right / max_left を素朴解と突き合わせる
template <class Mo, class F>
static void check_binsearch(const SegTree<Mo>& seg, const std::vector<ll>& a, F f) {
    int n = (int)a.size();
    for (int l = 0; l <= n; l++) {
        int expect = l;
        typename Mo::S sm = Mo::e();
        while (expect < n) {
            typename Mo::S nx = Mo::op(sm, a[expect]);
            if (!f(nx)) break;
            sm = nx;
            expect++;
        }
        assert(seg.max_right(l, f) == expect);
    }
    for (int r = 0; r <= n; r++) {
        int expect = r;
        typename Mo::S sm = Mo::e();
        while (expect > 0) {
            typename Mo::S nx = Mo::op(a[expect - 1], sm);
            if (!f(nx)) break;
            sm = nx;
            expect--;
        }
        assert(seg.max_left(r, f) == expect);
    }
}

static void test_binsearch(std::mt19937_64& rng) {
    for (int iter = 0; iter < 1500; iter++) {
        int n = (int)(rng() % 40);
        // --- Sum: 「和が X 以下」 (非負値なので単調) ---
        {
            std::vector<ll> a(n);
            for (auto& x : a) x = rnd_range(rng, 0, 20);
            SegTree<MonoidSum> seg(a);
            ll X = rnd_range(rng, 0, 200);
            check_binsearch(seg, a, [X](ll v) { return v <= X; });
        }
        // --- Min: 「最小値が X 以上」 ---
        {
            std::vector<ll> a(n);
            for (auto& x : a) x = rnd_range(rng, -20, 20);
            SegTree<MonoidMin> seg(a);
            ll X = rnd_range(rng, -25, 25);
            check_binsearch(seg, a, [X](ll v) { return v >= X; });
        }
        // --- Max: 「最大値が X 以下」 ---
        {
            std::vector<ll> a(n);
            for (auto& x : a) x = rnd_range(rng, -20, 20);
            SegTree<MonoidMax> seg(a);
            ll X = rnd_range(rng, -25, 25);
            check_binsearch(seg, a, [X](ll v) { return v <= X; });
        }
    }
}

static void random_test() {
    std::mt19937_64 rng(20260916);
    const std::vector<ll> lcm_pool = divisors_of_720720();

    // 小さい n を全網羅的に、全モノイドで (複数シード分繰り返す)
    for (int rep = 0; rep < 5; rep++) {
        for (int n = 0; n <= 17; n++) {
            test_one_monoid<MonoidSum>("Sum", rng, n, 300, -1000000000LL, 1000000000LL);
            test_one_monoid<MonoidMin>("Min", rng, n, 300, -1000000000LL, 1000000000LL);
            test_one_monoid<MonoidMax>("Max", rng, n, 300, -1000000000LL, 1000000000LL);
            test_one_monoid<MonoidGcd>("Gcd", rng, n, 300, 0, 720720);
            test_one_monoid<MonoidLcm>("Lcm", rng, n, 300, 0, 0, lcm_pool);
            test_one_monoid<MonoidXor>("Xor", rng, n, 300, 0, (1LL << 40));
            test_one_monoid<MonoidAnd>("And", rng, n, 300, -(1LL << 40), (1LL << 40));
            test_one_monoid<MonoidOr>("Or", rng, n, 300, -(1LL << 40), (1LL << 40));
            test_one_monoid<MonoidLastNonZero>("NonComm", rng, n, 300, -3, 3);
        }
    }
    // やや大きい n でも (2 のべき乗ちょうど / ±1 も含める)
    for (int rep = 0; rep < 6; rep++) {
        static const int ns[] = {63, 64, 65, 127, 128, 255};
        int n = ns[rep];
        test_one_monoid<MonoidSum>("Sum", rng, n, 4000, -1000000000000LL, 1000000000000LL);
        test_one_monoid<MonoidMin>("Min", rng, n, 4000, -1000000000000LL, 1000000000000LL);
        test_one_monoid<MonoidXor>("Xor", rng, n, 4000, 0, (1LL << 62));
        test_one_monoid<MonoidAnd>("And", rng, n, 4000, -(1LL << 62), (1LL << 62));
    }
    // 二分探索
    test_binsearch(rng);

    // 境界値 (単位元まわり)
    {
        std::vector<ll> a = {std::numeric_limits<ll>::max(), std::numeric_limits<ll>::min(), 0};
        SegTree<MonoidMin> smin(a);
        assert(smin.range_query(0, 3) == std::numeric_limits<ll>::min());
        SegTree<MonoidMax> smax(a);
        assert(smax.range_query(0, 3) == std::numeric_limits<ll>::max());
        SegTree<MonoidMin> s1(a);
        assert(s1.range_query(1, 1) == MonoidMin::e());  // 空区間は単位元
    }
    std::cout << "[segtree] random test: OK" << std::endl;
}

}  // namespace tester

// =============================================================================
//  使用例デモ : 全機能を 1 箇所で使う。実戦ではここを見れば思い出せる。
// =============================================================================
static void demo() {
    std::cout << "\n===== ① Segment Tree 使用例デモ =====\n";

    // 【構築1】初期値の配列から構築する  O(n)
    std::vector<ll> a = {5, 3, 8, 1, 9, 2, 7, 4};
    Seg seg(a);

    // 【構築2】要素数だけ指定する。全要素が単位元 e() で初期化される。
    Seg seg_empty(8);
    std::cout << "size()            = " << seg.size() << "\n";
    std::cout << "空木の all_query() = " << seg_empty.all_query() << "   (= 単位元 e())\n";

    // 【point_set】1 点更新: a[3] を 10 にする
    seg.point_set(3, 10);

    // 【point_get】1 点取得
    std::cout << "point_get(3)      = " << seg.point_get(3) << "\n";

    // 【range_query】半開区間 [2, 6) のモノイド積（今は M = MonoidSum なので区間和）
    std::cout << "range_query(2, 6) = " << seg.range_query(2, 6) << "\n";

    // 【range_query の空区間】l == r のときは単位元が返る
    std::cout << "range_query(4, 4) = " << seg.range_query(4, 4) << "   (空区間 → e())\n";

    // 【all_query】全体の積。range_query(0, n) と同じだが O(1)
    std::cout << "all_query()       = " << seg.all_query() << "\n";

    // 【max_right】左端 l を固定し、条件を満たす最大の r を返す（f(e()) は true 必須）
    //   例: 区間和が 20 以下に収まる最大の r
    {
        SegTree<MonoidSum> s(a);
        int r = s.max_right(1, [](ll v) { return v <= 20; });
        std::cout << "max_right(1, sum<=20) = " << r
                  << "   → [1, " << r << ") の和 = " << s.range_query(1, r) << "\n";
    }

    // 【max_left】右端 r を固定し、条件を満たす最小の l を返す（左へ最大限のばす）
    {
        SegTree<MonoidSum> s(a);
        int l = s.max_left(7, [](ll v) { return v <= 20; });
        std::cout << "max_left(7, sum<=20)  = " << l
                  << "   → [" << l << ", 7) の和 = " << s.range_query(l, 7) << "\n";
    }

    // 【モノイドの切り替え】使いたいものを using M = ...; で選ぶだけ
    {
        std::vector<ll> b = {12, 18, 24, 6, 30};
        std::cout << "--- 同じ配列 {12,18,24,6,30} の [0,5) を各モノイドで ---\n";
        std::cout << "  Sum = " << SegTree<MonoidSum>(b).range_query(0, 5) << "\n";
        std::cout << "  Min = " << SegTree<MonoidMin>(b).range_query(0, 5) << "\n";
        std::cout << "  Max = " << SegTree<MonoidMax>(b).range_query(0, 5) << "\n";
        std::cout << "  Gcd = " << SegTree<MonoidGcd>(b).range_query(0, 5) << "\n";
        std::cout << "  Lcm = " << SegTree<MonoidLcm>(b).range_query(0, 5) << "\n";
        std::cout << "  Xor = " << SegTree<MonoidXor>(b).range_query(0, 5) << "\n";
        std::cout << "  And = " << SegTree<MonoidAnd>(b).range_query(0, 5) << "\n";
        std::cout << "  Or  = " << SegTree<MonoidOr>(b).range_query(0, 5) << "\n";
    }
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    tester::random_test();
    demo();
    return 0;
}
