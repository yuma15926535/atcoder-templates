// =============================================================================
//  ② Lazy Segment Tree (区間作用 / 区間取得)   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  値の型      : long long 固定
//  区間の流儀  : 半開区間 [l, r)  / 0-indexed
//  依存        : 標準ライブラリのみ (C++17)
//
//  作用は「アフィン変換 x -> b*x + c」ひとつに統一し、
//      range_add(l, r, v)    = range_affine(l, r, 1, v)
//      range_assign(l, r, v) = range_affine(l, r, 0, v)
//  として実装している。
//  各ノードは (sum, min, max, len) を持つので、1 本で和・最小・最大が全部取れる。
//  b < 0 のときは min と max が入れ替わることに注意 (mapping 内で処理済み)。
//
//  計算量: 構築 O(n) / 各クエリ O(log n)
//
//  オーバーフロー注意: sum は long long。|値| * n が 9.2e18 を超えないこと。
//                      アフィンを繰り返すと値が指数的に増えるので注意。
// =============================================================================

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

using ll = long long;

struct LazySegTree {
    static constexpr ll INF = std::numeric_limits<ll>::max();  // 空区間の min として返る値

    // ---- ノードが持つ情報 ----
    struct S {
        ll sum, mn, mx, len;
    };
    // ---- 作用素: x -> b*x + c ----
    struct F {
        ll b, c;
    };

    static S e() { return S{0, INF, -INF, 0}; }  // 単位元 (len = 0)
    static S op(const S& x, const S& y) {
        return S{x.sum + y.sum, std::min(x.mn, y.mn), std::max(x.mx, y.mx), x.len + y.len};
    }
    static F id() { return F{1, 0}; }
    static S mapping(const F& f, const S& x) {
        if (x.len == 0) return x;  // 単位元には作用させない
        ll lo = f.b * x.mn + f.c, hi = f.b * x.mx + f.c;
        if (f.b < 0) std::swap(lo, hi);  // 負の傾きだと min/max が入れ替わる
        return S{f.b * x.sum + f.c * x.len, lo, hi, x.len};
    }
    // f(g(x)) : g を先に、f を後に適用
    static F composition(const F& f, const F& g) { return F{f.b * g.b, f.b * g.c + f.c}; }

    // -------------------------------------------------------------------------
    LazySegTree() : LazySegTree(0) {}
    explicit LazySegTree(int n) : LazySegTree(std::vector<ll>(n, 0)) {}
    explicit LazySegTree(const std::vector<ll>& v) : _n((int)v.size()) {
        _log = 0;
        while ((1 << _log) < _n) ++_log;
        _size = 1 << _log;
        d.assign(2 * (size_t)_size, e());
        lz.assign((size_t)_size, id());
        for (int i = 0; i < _n; i++) d[_size + i] = S{v[i], v[i], v[i], 1};
        for (int i = _size - 1; i >= 1; i--) update(i);
    }

    int size() const { return _n; }

    // ---- 点操作 ----
    void point_set(int p, ll x) {  // a[p] = x
        assert(0 <= p && p < _n);
        p += _size;
        for (int i = _log; i >= 1; i--) push(p >> i);
        d[p] = S{x, x, x, 1};
        for (int i = 1; i <= _log; i++) update(p >> i);
    }
    ll point_get(int p) {  // a[p]
        assert(0 <= p && p < _n);
        p += _size;
        for (int i = _log; i >= 1; i--) push(p >> i);
        return d[p].sum;  // len == 1 なので sum がそのまま値
    }

    // ---- 区間更新 ----
    void range_affine(int l, int r, ll b, ll c) { apply_range(l, r, F{b, c}); }  // a[i] = b*a[i]+c
    void range_add(int l, int r, ll v) { apply_range(l, r, F{1, v}); }           // a[i] += v
    void range_assign(int l, int r, ll v) { apply_range(l, r, F{0, v}); }        // a[i] = v

    // ---- 区間取得 ----
    S range_query(int l, int r) {  // sum/min/max/len をまとめて 1 回で取る
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return e();
        l += _size;
        r += _size;
        for (int i = _log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        S sml = e(), smr = e();
        while (l < r) {
            if (l & 1) sml = op(sml, d[l++]);
            if (r & 1) smr = op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }
        return op(sml, smr);
    }
    ll range_sum(int l, int r) { return range_query(l, r).sum; }  // 空区間は 0
    ll range_min(int l, int r) { return range_query(l, r).mn; }   // 空区間は  INF
    ll range_max(int l, int r) { return range_query(l, r).mx; }   // 空区間は -INF
    S all_query() const { return d[1]; }

   private:
    int _n = 0, _size = 1, _log = 0;
    std::vector<S> d;
    std::vector<F> lz;

    void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
    void all_apply(int k, const F& f) {
        d[k] = mapping(f, d[k]);
        if (k < _size) lz[k] = composition(f, lz[k]);
    }
    void push(int k) {
        all_apply(2 * k, lz[k]);
        all_apply(2 * k + 1, lz[k]);
        lz[k] = id();
    }
    void apply_range(int l, int r, const F& f) {
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
        l += _size;
        r += _size;
        for (int i = _log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        {
            int l2 = l, r2 = r;
            while (l < r) {
                if (l & 1) all_apply(l++, f);
                if (r & 1) all_apply(--r, f);
                l >>= 1;
                r >>= 1;
            }
            l = l2;
            r = r2;
        }
        for (int i = 1; i <= _log; i++) {
            if (((l >> i) << i) != l) update(l >> i);
            if (((r >> i) << i) != r) update((r - 1) >> i);
        }
    }
};

// =============================================================================
//  ここから下は「ランダムテスト」と「使用例デモ」。提出時は丸ごと消してよい。
// =============================================================================

namespace tester {

static ll rnd_range(std::mt19937_64& rng, ll lo, ll hi) {
    return lo + (ll)(rng() % (unsigned long long)(hi - lo + 1));
}
static ll absll(ll x) { return x < 0 ? -x : x; }

static void random_test() {
    std::mt19937_64 rng(20260916ULL);
    const ll VAL = 1000;

    for (int rep = 0; rep < 400; rep++) {
        int n = (int)rnd_range(rng, 1, 40);
        std::vector<ll> a(n);
        for (auto& x : a) x = rnd_range(rng, -VAL, VAL);
        LazySegTree seg(a);

        for (int t = 0; t < 400; t++) {
            // --- 値が発散しないようにガード (アフィンを繰り返すと指数的に増えるため) ---
            ll big = 0;
            for (ll x : a) big = std::max(big, absll(x));
            if (big > 1000000000000LL) {
                ll v = rnd_range(rng, -VAL, VAL);
                for (auto& x : a) x = v;
                seg.range_assign(0, n, v);
            }

            int l = (int)rnd_range(rng, 0, n), r = (int)rnd_range(rng, 0, n);
            if (l > r) std::swap(l, r);
            int p = (int)rnd_range(rng, 0, n - 1);
            int type = (int)rnd_range(rng, 0, 7);

            if (type == 0) {  // range_add
                ll v = rnd_range(rng, -VAL, VAL);
                for (int i = l; i < r; i++) a[i] += v;
                seg.range_add(l, r, v);
            } else if (type == 1) {  // range_assign
                ll v = rnd_range(rng, -VAL, VAL);
                for (int i = l; i < r; i++) a[i] = v;
                seg.range_assign(l, r, v);
            } else if (type == 2) {  // range_affine
                ll b = rnd_range(rng, -3, 3), c = rnd_range(rng, -VAL, VAL);
                for (int i = l; i < r; i++) a[i] = b * a[i] + c;
                seg.range_affine(l, r, b, c);
            } else if (type == 3) {  // point_set
                ll x = rnd_range(rng, -VAL, VAL);
                a[p] = x;
                seg.point_set(p, x);
            } else if (type == 4) {  // point_get
                assert(seg.point_get(p) == a[p]);
            } else if (type == 5) {  // range_sum
                ll s = 0;
                for (int i = l; i < r; i++) s += a[i];
                assert(seg.range_sum(l, r) == s);
            } else if (type == 6) {  // range_min
                ll m = LazySegTree::INF;
                for (int i = l; i < r; i++) m = std::min(m, a[i]);
                assert(seg.range_min(l, r) == m);
            } else {  // range_max
                ll m = -LazySegTree::INF;
                for (int i = l; i < r; i++) m = std::max(m, a[i]);
                assert(seg.range_max(l, r) == m);
            }
        }
        // 最後に全区間を突き合わせ
        for (int l = 0; l <= n; l++) {
            for (int r = l; r <= n; r++) {
                ll s = 0, mn = LazySegTree::INF, mx = -LazySegTree::INF;
                for (int i = l; i < r; i++) {
                    s += a[i];
                    mn = std::min(mn, a[i]);
                    mx = std::max(mx, a[i]);
                }
                LazySegTree::S q = seg.range_query(l, r);
                assert(q.sum == s && q.mn == mn && q.mx == mx && q.len == r - l);
            }
        }
    }

    // 2 のべき乗ちょうど / ±1 のサイズも確認
    for (int n : {1, 2, 3, 7, 8, 9, 63, 64, 65}) {
        std::vector<ll> a(n);
        for (auto& x : a) x = rnd_range(rng, -VAL, VAL);
        LazySegTree seg(a);
        for (int t = 0; t < 2000; t++) {
            int l = (int)rnd_range(rng, 0, n), r = (int)rnd_range(rng, 0, n);
            if (l > r) std::swap(l, r);
            ll b = rnd_range(rng, -1, 2), c = rnd_range(rng, -VAL, VAL);
            for (int i = l; i < r; i++) a[i] = b * a[i] + c;
            seg.range_affine(l, r, b, c);
            int l2 = (int)rnd_range(rng, 0, n), r2 = (int)rnd_range(rng, 0, n);
            if (l2 > r2) std::swap(l2, r2);
            ll s = 0, mn = LazySegTree::INF, mx = -LazySegTree::INF;
            for (int i = l2; i < r2; i++) {
                s += a[i];
                mn = std::min(mn, a[i]);
                mx = std::max(mx, a[i]);
            }
            assert(seg.range_sum(l2, r2) == s);
            assert(seg.range_min(l2, r2) == mn);
            assert(seg.range_max(l2, r2) == mx);
        }
    }

    std::cout << "[lazy_segtree] random test: OK" << std::endl;
}

}  // namespace tester

// =============================================================================
//  使用例デモ : 全機能を 1 箇所で使う。実戦ではここを見れば思い出せる。
// =============================================================================
static void dump(const char* label, LazySegTree& seg) {
    std::cout << label << " : ";
    for (int i = 0; i < seg.size(); i++) std::cout << seg.point_get(i) << (i + 1 == seg.size() ? '\n' : ' ');
}

static void demo() {
    std::cout << "\n===== ② Lazy Segment Tree 使用例デモ =====\n";

    // 【構築1】初期値の配列から  O(n)
    std::vector<ll> a = {5, 3, 8, 1, 9, 2, 7, 4};
    LazySegTree seg(a);
    dump("初期状態        ", seg);

    // 【構築2】要素数だけ指定 (全要素 0 で初期化)
    LazySegTree seg_zero(8);
    std::cout << "LazySegTree(8) の全体和 = " << seg_zero.all_query().sum << "\n";

    // 【range_add】[1, 4) に +100
    seg.range_add(1, 4, 100);
    dump("range_add(1,4,+100)", seg);

    // 【range_assign】[5, 8) を 7 で塗りつぶす
    seg.range_assign(5, 8, 7);
    dump("range_assign(5,8,7)", seg);

    // 【range_affine】[0, 4) に x -> 2x + 1
    seg.range_affine(0, 4, 2, 1);
    dump("range_affine(0,4,b=2,c=1)", seg);

    // 【range_affine で符号反転】x -> -1*x + 0  (min と max が入れ替わるケース)
    seg.range_affine(4, 6, -1, 0);
    dump("range_affine(4,6,b=-1,c=0)", seg);

    // 【point_set】1 点代入
    seg.point_set(2, -50);
    // 【point_get】1 点取得
    std::cout << "point_get(2)       = " << seg.point_get(2) << "\n";

    // 【range_sum / range_min / range_max】区間取得
    std::cout << "range_sum(0, 8)    = " << seg.range_sum(0, 8) << "\n";
    std::cout << "range_min(0, 8)    = " << seg.range_min(0, 8) << "\n";
    std::cout << "range_max(0, 8)    = " << seg.range_max(0, 8) << "\n";
    std::cout << "range_sum(2, 5)    = " << seg.range_sum(2, 5) << "\n";

    // 【range_query】sum/min/max/len を 1 回の走査でまとめて取る (3 回呼ぶより速い)
    LazySegTree::S q = seg.range_query(1, 6);
    std::cout << "range_query(1,6)   = {sum=" << q.sum << ", min=" << q.mn << ", max=" << q.mx
              << ", len=" << q.len << "}\n";

    // 【all_query】全体を O(1) で取得 (遅延は根には残らないので push 不要)
    std::cout << "all_query()        = {sum=" << seg.all_query().sum << ", min=" << seg.all_query().mn
              << ", max=" << seg.all_query().mx << "}\n";

    // 【空区間】l == r のときは sum=0, min=INF, max=-INF
    std::cout << "range_sum(3,3)     = " << seg.range_sum(3, 3) << "  (空区間)\n";
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    tester::random_test();
    demo();
    return 0;
}
