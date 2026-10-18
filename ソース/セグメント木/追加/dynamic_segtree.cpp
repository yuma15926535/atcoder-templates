// =============================================================================
//  ④ Dynamic Segment Tree (動的セグメント木・一般モノイド版)
//                                          --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  値の型      : long long 固定
//  座標        : long long の全範囲 [LLONG_MIN, LLONG_MAX] をそのまま扱える。
//                座標圧縮も範囲指定も不要。必要なノードだけ動的に作る。
//  区間の流儀  : 半開区間 [l, r)
//                ※ 半開なので range_query で LLONG_MAX 自身は含められない。
//                   その 1 点が要るときは point_get を使うこと。
//  依存        : 標準ライブラリのみ (C++17)
//
//  一度も point_set していない座標の値は「単位元 e()」として扱われる。
//    Sum なら 0 / Min なら +INF / Max なら -INF / And なら全ビット 1 …
//
//  実装: 葉だけが実データを持つ「パス圧縮済み二分トライ (基数木)」。
//        単一の子しか持たない連鎖を作らないので、
//          ノード数 <= 2 * (point_set した相異なる座標数) - 1
//        となり、64 段そのまま掘る実装に比べてメモリが 30 分の 1 以下になる。
//        左右の順序は保たれるので非可換モノイドでも正しく動く。
//
// -----------------------------------------------------------------------------
//  ★ 設計の要点 (ここが「動的」セグメント木として効く部分)
//
//  (1) 取得 (range_query / point_get / all_query) はノードを一切生成せず、
//      木も書き換えない (すべて const)。遅延伝播が無いので push そのものが存在しない。
//      → 取得を何百万回やっても node_count は 1 も増えない。
//
//  (2) point_get は range_query 経由ではなく専用の 1 本道降下 (ループ)。
//      根から各段で必要な子 1 つだけを辿る。未生成の枝に出たら単位元 e() を返す。
//
//  (3) range_query は子へ降りる前に「その子の担当区間と交差するか」を
//      中点 mid だけで判定してから再帰する。無駄な再帰呼び出しをしない。
//
//  (4) ノードが作られるのは point_set のときだけ。しかも新しい座標 1 個につき
//      高々 2 個 (葉 1 + 分岐点 1)。既存座標の上書きは 0 個。
//
// -----------------------------------------------------------------------------
//  計算量 (木の深さは高々 D = 64)
//    point_set   : O(D)。新規座標 1 個につき生成ノードは高々 2 個
//    point_get   : O(D)。1 本道。生成ノード 0 個
//    range_query : O(D)。訪れるノードは高々 2D 個。生成ノード 0 個
//    all_query   : O(1)。生成ノード 0 個
//
//  メモリ : 1 ノード 32 バイト、ノード数は「相異なる座標数 * 2 - 1」以下。
//           実測: 64bit 全域ランダムな座標に point_set 2*10^5 回 → 40 万ノード / 約 12 MB
//           (※ 遅延なしの ④ はこの通り軽い。区間更新が要る ⑤⑥ は桁違いに重いので注意)
// =============================================================================

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <vector>

using ll = long long;
using ull = unsigned long long;

// =============================================================================
//  モノイド定義 (① segtree.cpp と同じもの)
// =============================================================================
struct MonoidSum {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return a + b; }
};
struct MonoidMin {
    using S = ll;
    static S e() { return std::numeric_limits<ll>::max(); }
    static S op(S a, S b) { return a < b ? a : b; }
};
struct MonoidMax {
    using S = ll;
    static S e() { return std::numeric_limits<ll>::min(); }
    static S op(S a, S b) { return a > b ? a : b; }
};
struct MonoidGcd {  // 注: 負数は絶対値の gcd になる
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return std::gcd(a, b); }
};
struct MonoidLcm {  // 注: 非負整数を想定。0 は吸収元。オーバーフロー注意
    using S = ll;
    static S e() { return 1; }
    static S op(S a, S b) {
        if (a == 0 || b == 0) return 0;
        S g = std::gcd(a, b);
        return (a / g) * b;
    }
};
struct MonoidXor {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return a ^ b; }
};
struct MonoidAnd {
    using S = ll;
    static S e() { return ~0LL; }
    static S op(S a, S b) { return a & b; }
};
struct MonoidOr {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return a | b; }
};

// =============================================================================
//  Dynamic Segment Tree 本体
// =============================================================================
template <class M>
struct DynamicSegTree {
    using S = typename M::S;

    DynamicSegTree() = default;
    void reserve(size_t n) { t.reserve(n); }  // 事前確保しておくと速い (任意)

    // a[key] = x
    void point_set(ll key, S x) {
        int nr = set_rec(root, enc(key), x);
        root = nr;
    }
    // a[key] (未設定なら e())
    // 専用の 1 本道降下。const / ノード生成 0 個。
    S point_get(ll key) const {
        ull k = enc(key);
        int p = root;
        while (p != -1) {
            const Node& nd = t[p];
            if (k < nd.lo || nd.hi < k) return M::e();  // その枝には無い
            if (nd.lo == nd.hi) return nd.prod;         // 葉 = その座標の値
            ull mid = nd.lo + (nd.hi - nd.lo) / 2;
            p = (k <= mid) ? nd.lc : nd.rc;
        }
        return M::e();  // 未設定の座標は単位元
    }
    // [l, r) のモノイド積 (l >= r なら e())。const / ノード生成 0 個。
    S range_query(ll l, ll r) const {
        if (l >= r || root == -1) return M::e();
        ull ql = enc(l), qr = enc(r - 1);
        if (t[root].hi < ql || qr < t[root].lo) return M::e();
        return query_rec(root, ql, qr);
    }
    // 全体のモノイド積 O(1)
    S all_query() const { return root == -1 ? M::e() : t[root].prod; }

    size_t node_count() const { return t.size(); }

   private:
    struct Node {
        S prod;      // 部分木のモノイド積
        ull lo, hi;  // このノードが担当する整列区間 (閉区間, 変換後の符号なし座標)
        int lc, rc;  // 子 (-1 = なし)。lo == hi の葉のみが実データを持つ
    };
    std::vector<Node> t;
    int root = -1;

    // long long の順序を符号なし 64bit の順序に移す
    static ull enc(ll x) { return (ull)x ^ (1ULL << 63); }

    static int top_bit(ull x) {  // x != 0 の最上位ビット位置
        int r = 0;
        if (x >> 32) { x >>= 32; r += 32; }
        if (x >> 16) { x >>= 16; r += 16; }
        if (x >> 8) { x >>= 8; r += 8; }
        if (x >> 4) { x >>= 4; r += 4; }
        if (x >> 2) { x >>= 2; r += 2; }
        if (x >> 1) { r += 1; }
        return r;
    }
    // a, b を両方含む最小の整列区間
    static void cover(ull a, ull b, ull& lo, ull& hi) {
        if (a == b) { lo = hi = a; return; }
        int h = top_bit(a ^ b);
        ull mask = (h == 63) ? ~0ULL : ((1ULL << (h + 1)) - 1);
        lo = a & ~mask;
        hi = a | mask;
    }

    int new_leaf(ull k, S x) {
        t.push_back(Node{x, k, k, -1, -1});
        return (int)t.size() - 1;
    }
    int new_internal(ull lo, ull hi) {
        t.push_back(Node{M::e(), lo, hi, -1, -1});
        return (int)t.size() - 1;
    }
    void pull(int p) {
        S a = (t[p].lc == -1) ? M::e() : t[t[p].lc].prod;
        S b = (t[p].rc == -1) ? M::e() : t[t[p].rc].prod;
        t[p].prod = M::op(a, b);
    }

    int set_rec(int p, ull k, S x) {
        if (p == -1) return new_leaf(k, x);
        if (t[p].lo <= k && k <= t[p].hi) {
            if (t[p].lo == t[p].hi) {  // 同じ座標の葉
                t[p].prod = x;
                return p;
            }
            ull mid = t[p].lo + (t[p].hi - t[p].lo) / 2;
            if (k <= mid) {
                int c = set_rec(t[p].lc, k, x);
                t[p].lc = c;
            } else {
                int c = set_rec(t[p].rc, k, x);
                t[p].rc = c;
            }
            pull(p);
            return p;
        }
        // p の担当区間の外 → 両方を含む最小の整列区間を作って左右に振り分ける
        ull a = std::min(t[p].lo, k), b = std::max(t[p].hi, k);
        ull nlo, nhi;
        cover(a, b, nlo, nhi);
        int leaf = new_leaf(k, x);
        int q = new_internal(nlo, nhi);
        ull mid = nlo + (nhi - nlo) / 2;
        if (k <= mid) { t[q].lc = leaf; t[q].rc = p; }
        else { t[q].lc = p; t[q].rc = leaf; }
        pull(q);
        return q;
    }

    // 前提: p != -1 かつ [t[p].lo, t[p].hi] と [ql, qr] は交差している。
    // 子へ降りる前に中点 mid だけで交差判定するので、無駄な再帰呼び出しをしない。
    S query_rec(int p, ull ql, ull qr) const {
        const Node& nd = t[p];
        if (ql <= nd.lo && nd.hi <= qr) return nd.prod;  // 丸ごと入る
        // ここへ来た時点で葉ではない (葉は lo == hi なので上で必ず返る)
        ull mid = nd.lo + (nd.hi - nd.lo) / 2;
        int lc = nd.lc, rc = nd.rc;
        bool hitL = (lc != -1) && !(qr < nd.lo || mid < ql);
        bool hitR = (rc != -1) && !(qr < mid + 1 || nd.hi < ql);
        if (hitL && hitR) return M::op(query_rec(lc, ql, qr), query_rec(rc, ql, qr));
        if (hitL) return query_rec(lc, ql, qr);
        if (hitR) return query_rec(rc, ql, qr);
        return M::e();
    }
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

using DSeg = DynamicSegTree<M>;

// =============================================================================
//  ここから下は「ランダムテスト」と「使用例デモ」。提出時は丸ごと消してよい。
// =============================================================================

namespace tester {

// 非可換モノイド (左右の順序が保たれているかの検証用): 「最後の非ゼロ値」
// パス圧縮しても葉の左右順序が壊れていないことを、これで確認する。
struct MonoidLastNonZero {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return b != 0 ? b : a; }
};

static ll rnd_range(std::mt19937_64& rng, ll lo, ll hi) {
    return lo + (ll)(rng() % (ull)(hi - lo + 1));
}

// LCM テスト用: 720720 (= 2^4*3^2*5*7*11*13) の約数だけを使えば lcm は 720720 以下で収まる
static std::vector<ll> divisors_of_720720() {
    std::vector<ll> d;
    for (ll x = 1; x <= 720720; x++)
        if (720720 % x == 0) d.push_back(x);
    return d;
}

// 素朴解: std::map (キー昇順) を走査
template <class Mo>
static typename Mo::S naive(const std::map<ll, ll>& mp, ll l, ll r) {
    typename Mo::S s = Mo::e();
    if (l >= r) return s;
    for (auto it = mp.lower_bound(l); it != mp.end() && it->first < r; ++it) s = Mo::op(s, it->second);
    return s;
}

// 座標候補: 極端な値 + 小さい値 をたくさん混ぜる
static std::vector<ll> make_coords(std::mt19937_64& rng, int cnt) {
    std::vector<ll> c;
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();
    c.push_back(LMIN);
    c.push_back(LMIN + 1);
    c.push_back(LMIN / 2);
    c.push_back(-(1LL << 62));
    c.push_back(-(1LL << 31));
    c.push_back(-2);
    c.push_back(-1);
    c.push_back(0);
    c.push_back(1);
    c.push_back(2);
    c.push_back(1LL << 31);
    c.push_back(1LL << 62);
    c.push_back(LMAX / 2);
    c.push_back(LMAX - 1);
    c.push_back(LMAX);
    while ((int)c.size() < cnt) {
        int kind = (int)(rng() % 3);
        if (kind == 0) c.push_back(rnd_range(rng, -8, 8));
        else if (kind == 1) c.push_back(rnd_range(rng, LMIN / 2, LMAX / 2));
        else c.push_back((ll)rng());  // 完全ランダム (全域)
    }
    return c;
}

template <class Mo>
static void test_one_monoid(std::mt19937_64& rng, int ncoord, int q, ll vlo, ll vhi,
                            const std::vector<ll>& pool = std::vector<ll>()) {
    auto gen = [&]() -> ll { return pool.empty() ? rnd_range(rng, vlo, vhi) : pool[rng() % pool.size()]; };
    std::vector<ll> coords = make_coords(rng, ncoord);
    DynamicSegTree<Mo> seg;
    std::map<ll, ll> mp;

    for (int t = 0; t < q; t++) {
        int type = (int)(rng() % 4);
        if (type == 0) {  // point_set
            ll k = coords[rng() % coords.size()];
            ll x = gen();
            seg.point_set(k, x);
            mp[k] = x;
        } else if (type == 1) {  // point_get
            ll k = coords[rng() % coords.size()];
            typename Mo::S expect = mp.count(k) ? (typename Mo::S)mp[k] : Mo::e();
            assert(seg.point_get(k) == expect);
        } else if (type == 2) {  // range_query
            ll l = coords[rng() % coords.size()], r = coords[rng() % coords.size()];
            if (l > r) std::swap(l, r);
            assert(seg.range_query(l, r) == naive<Mo>(mp, l, r));
            assert(seg.range_query(r, l) == Mo::e() || l == r);  // l >= r は空
        } else {  // all_query
            typename Mo::S s = Mo::e();
            for (auto& kv : mp) s = Mo::op(s, kv.second);
            assert(seg.all_query() == s);
        }
    }
    // 最後に全ペアで突き合わせ
    std::vector<ll> cs = coords;
    std::sort(cs.begin(), cs.end());
    cs.erase(std::unique(cs.begin(), cs.end()), cs.end());
    for (size_t i = 0; i < cs.size(); i++)
        for (size_t j = i; j < cs.size(); j++)
            assert(seg.range_query(cs[i], cs[j]) == naive<Mo>(mp, cs[i], cs[j]));
}

// -----------------------------------------------------------------------------
//  境界値テスト: LLONG_MIN / LLONG_MAX 近傍と符号のまたぎ
// -----------------------------------------------------------------------------
static void test_boundary() {
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();
    std::vector<ll> keys = {LMIN, LMIN + 1, LMIN + 2, -3, -2, -1, 0, 1, 2, 3, LMAX - 2, LMAX - 1, LMAX};

    DynamicSegTree<MonoidSum> seg;
    for (size_t i = 0; i < keys.size(); i++) seg.point_set(keys[i], (ll)i + 1);
    for (size_t i = 0; i < keys.size(); i++) assert(seg.point_get(keys[i]) == (ll)i + 1);
    assert(seg.point_get(12345) == 0);  // 未設定は単位元

    assert(seg.all_query() == 91);                 // 1..13 の総和
    assert(seg.range_query(LMIN, LMAX) == 78);     // 半開なので LMAX は入らない
    assert(seg.range_query(LMIN, LMIN + 1) == 1);  // 1 要素
    assert(seg.range_query(LMAX - 1, LMAX) == 12);
    assert(seg.range_query(-3, 4) == 4 + 5 + 6 + 7 + 8 + 9 + 10);  // 符号をまたぐ
    assert(seg.range_query(5, 5) == 0);                            // 空区間
    assert(seg.range_query(10, 3) == 0);                           // 逆転区間

    // min / max モノイドでも境界を確認
    DynamicSegTree<MonoidMin> smin;
    DynamicSegTree<MonoidMax> smax;
    for (size_t i = 0; i < keys.size(); i++) {
        smin.point_set(keys[i], (ll)i + 1);
        smax.point_set(keys[i], (ll)i + 1);
    }
    assert(smin.range_query(LMIN, LMAX) == 1);
    assert(smax.range_query(LMIN, LMAX) == 12);
    assert(smin.range_query(4, 5) == MonoidMin::e());  // 未設定だけの区間 → 単位元
    assert(smax.range_query(4, 5) == MonoidMax::e());
}

// ★ 大量の取得でノードが 1 つも増えないこと (全域ランダム座標)
static void test_query_allocates_nothing(std::mt19937_64& rng) {
    DynamicSegTree<MonoidSum> seg;
    for (int i = 0; i < 2000; i++) seg.point_set((ll)rng(), rnd_range(rng, -1000, 1000));
    size_t before = seg.node_count();
    ll acc = 0;
    for (int i = 0; i < 100000; i++) {
        ll l = (ll)rng(), r = (ll)rng();
        if (l > r) std::swap(l, r);
        acc ^= seg.range_query(l, r);
        acc ^= seg.point_get((ll)rng());
        acc ^= seg.all_query();
    }
    (void)acc;
    assert(seg.node_count() == before);
    // 既存座標の上書きでもノードは増えない
    std::vector<ll> ks;
    for (int i = 0; i < 100; i++) ks.push_back((ll)rng());
    for (ll k : ks) seg.point_set(k, 1);
    size_t mid = seg.node_count();
    for (ll k : ks) seg.point_set(k, 2);
    assert(seg.node_count() == mid);
}

// -----------------------------------------------------------------------------
//  ノード数ベンチ (A: 更新中心 / B: 取得中心 / C: 取得のみ)
// -----------------------------------------------------------------------------
static void bench_nodes() {
    std::cout << "\n----- ノード数ベンチ (座標は 64bit 全域ランダム) -----\n";
    const int Q = 200000;
    {  // A. 更新中心
        std::mt19937_64 rng(1);
        DynamicSegTree<MonoidSum> seg;
        for (int i = 0; i < Q; i++) seg.point_set((ll)rng(), (ll)(rng() % 1000));
        std::cout << "A point_set " << Q << " 回      : node_count = " << seg.node_count()
                  << "  (1 更新あたり " << (double)seg.node_count() / Q << " ノード, 約 "
                  << (seg.node_count() * 32 / 1048576) << " MB)\n";
    }
    {  // B. 少数更新 + 大量取得
        std::mt19937_64 rng(2);
        DynamicSegTree<MonoidSum> seg;
        for (int i = 0; i < 1000; i++) seg.point_set((ll)rng(), (ll)(rng() % 1000));
        size_t before = seg.node_count();
        ll acc = 0;
        for (int i = 0; i < Q; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            acc ^= seg.range_query(l, r) ^ seg.point_get((ll)rng());
        }
        (void)acc;
        std::cout << "B 更新 1000 + 取得 " << Q << " 回 : node_count = " << before << " → "
                  << seg.node_count() << "  (増加 " << seg.node_count() - before << ")\n";
    }
    {  // C. 取得のみ
        std::mt19937_64 rng(3);
        DynamicSegTree<MonoidSum> seg;
        for (int i = 0; i < 50; i++) seg.point_set((ll)rng(), (ll)(rng() % 1000));
        // (遅延伝播が無いのでこの木は元々 query で増えようがない)
        size_t before = seg.node_count();
        ll acc = 0;
        for (int i = 0; i < Q; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            acc ^= seg.range_query(l, r) ^ seg.point_get((ll)rng()) ^ seg.all_query();
        }
        (void)acc;
        std::cout << "C 取得のみ " << Q << " 回      : node_count = " << before << " → "
                  << seg.node_count() << "  (増加 " << seg.node_count() - before << ")\n";
    }
}

static void random_test() {
    std::mt19937_64 rng(20260916ULL);
    const std::vector<ll> lcm_pool = divisors_of_720720();
    for (int rep = 0; rep < 12; rep++) {
        int nc = (int)rnd_range(rng, 16, 40);
        test_one_monoid<MonoidSum>(rng, nc, 800, -1000000000LL, 1000000000LL);
        test_one_monoid<MonoidMin>(rng, nc, 800, -1000000000LL, 1000000000LL);
        test_one_monoid<MonoidMax>(rng, nc, 800, -1000000000LL, 1000000000LL);
        test_one_monoid<MonoidGcd>(rng, nc, 800, 0, 720720);
        test_one_monoid<MonoidLcm>(rng, nc, 800, 0, 0, lcm_pool);
        test_one_monoid<MonoidXor>(rng, nc, 800, 0, (1LL << 40));
        test_one_monoid<MonoidAnd>(rng, nc, 800, -(1LL << 40), (1LL << 40));
        test_one_monoid<MonoidOr>(rng, nc, 800, -(1LL << 40), (1LL << 40));
        test_one_monoid<MonoidLastNonZero>(rng, nc, 800, -3, 3);
    }

    // 連続した座標 0..N-1 だけを使い、静的な素朴配列と突き合わせる
    {
        const int N = 300;
        std::vector<ll> a(N, 0);
        DynamicSegTree<MonoidSum> seg;
        for (int t = 0; t < 20000; t++) {
            int type = (int)(rng() % 2);
            if (type == 0) {
                int p = (int)(rng() % N);
                ll x = rnd_range(rng, -1000000, 1000000);
                a[p] = x;
                seg.point_set(p, x);
            } else {
                int l = (int)(rng() % (N + 1)), r = (int)(rng() % (N + 1));
                if (l > r) std::swap(l, r);
                ll s = 0;
                for (int i = l; i < r; i++) s += a[i];
                assert(seg.range_query(l, r) == s);
            }
        }
        // ノード数は「相異なる座標数 * 2 - 1」以下
        assert(seg.node_count() <= (size_t)(2 * N));
    }

    test_boundary();
    test_query_allocates_nothing(rng);
    std::cout << "[dynamic_segtree] random test: OK" << std::endl;
}

}  // namespace tester

// =============================================================================
//  使用例デモ : 全機能を 1 箇所で使う。実戦ではここを見れば思い出せる。
// =============================================================================
static void demo() {
    std::cout << "\n===== ④ Dynamic Segment Tree 使用例デモ =====\n";
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();

    // 【構築】サイズも座標範囲も指定しない。そのまま使える。
    DSeg seg;
    // 【reserve】使う座標数が読めるなら事前確保すると速い (任意)
    seg.reserve(1 << 10);

    // 【point_set】long long のどんな座標でも OK
    seg.point_set(0, 10);
    seg.point_set(5, 20);
    seg.point_set(1000000000000LL, 30);   // 10^12
    seg.point_set(-1000000000000LL, 40);  // -10^12
    seg.point_set(LMIN, 50);              // 最小値
    seg.point_set(LMAX, 60);              // 最大値

    // 【point_get】未設定の座標は単位元 e() が返る
    std::cout << "point_get(5)        = " << seg.point_get(5) << "\n";
    std::cout << "point_get(LLONG_MAX)= " << seg.point_get(LMAX) << "\n";
    std::cout << "point_get(777)      = " << seg.point_get(777) << "   (未設定 → 単位元 e())\n";

    // 【point_set の上書き】同じ座標に再代入
    seg.point_set(5, 21);
    std::cout << "上書き後 point_get(5) = " << seg.point_get(5) << "\n";

    // 【range_query】半開区間 [l, r)
    std::cout << "range_query(0, 6)                 = " << seg.range_query(0, 6) << "\n";
    std::cout << "range_query(-10^12, 10^12+1)      = "
              << seg.range_query(-1000000000000LL, 1000000000001LL) << "\n";
    std::cout << "range_query(LLONG_MIN, LLONG_MAX) = " << seg.range_query(LMIN, LMAX)
              << "   (半開なので LLONG_MAX 自身は含まれない)\n";

    // 【空区間】l >= r なら単位元
    std::cout << "range_query(10, 10)               = " << seg.range_query(10, 10) << "   (空区間)\n";

    // 【all_query】全座標のモノイド積を O(1) で。LLONG_MAX も含まれる。
    std::cout << "all_query()                       = " << seg.all_query() << "\n";

    // 【node_count】生成されたノード数 (メモリの目安)。取得では絶対に増えない。
    size_t before = seg.node_count();
    ll acc = 0;
    for (int i = 0; i < 10000; i++) acc ^= seg.range_query(-i, i + 1) ^ seg.point_get(i);
    (void)acc;
    std::cout << "node_count()                      = " << before << " → " << seg.node_count()
              << "  (座標 6 個ぶん / 取得 1 万回してもノードは増えない)\n";

    // 【モノイドの切り替え】using M = ...; を差し替えるだけ
    {
        std::cout << "--- 座標 {-10^18, 0, 10^18} に {12, 18, 30} を置いて各モノイドで ---\n";
        auto build = [](auto tag) {
            using Mo = decltype(tag);
            DynamicSegTree<Mo> s;
            s.point_set(-1000000000000000000LL, 12);
            s.point_set(0, 18);
            s.point_set(1000000000000000000LL, 30);
            return s.all_query();
        };
        std::cout << "  Sum = " << build(MonoidSum{}) << "\n";
        std::cout << "  Min = " << build(MonoidMin{}) << "\n";
        std::cout << "  Max = " << build(MonoidMax{}) << "\n";
        std::cout << "  Gcd = " << build(MonoidGcd{}) << "\n";
        std::cout << "  Lcm = " << build(MonoidLcm{}) << "\n";
        std::cout << "  Xor = " << build(MonoidXor{}) << "\n";
        std::cout << "  And = " << build(MonoidAnd{}) << "\n";
        std::cout << "  Or  = " << build(MonoidOr{}) << "\n";
    }
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    tester::random_test();
    tester::bench_nodes();
    demo();
    return 0;
}
