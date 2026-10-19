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
//  二分探索 (セグ木上の二分探索。ACL の max_right / min_left と同じ仕様):
//      max_right(l, g) : g(prod[l, r)) が true になる最大の r   (g は S を受け取る)
//      min_left(r, g)  : g(prod[l, r)) が true になる最小の l
//      find_first_geq(l, x) など : 「a[i] >= x となる最初の i」のような値の検索 (中身は上の 2 つ)
//
//  計算量: 構築 O(n) / 各クエリ・二分探索 O(log n)
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

    // ---- 二分探索 (ACL の max_right / min_left と同じ仕様) ----
    //  g は「const S& を受け取って bool を返す関数」。g(e()) == true であること。
    //  g は単調であること (区間を広げると true → false の向きにしか変わらない)。
    //  例: [&](const LazySegTree::S& s) { return s.mx < x; }       (区間の最大値 < x)
    //      [&](const LazySegTree::S& s) { return s.sum <= k; }      (区間の和 <= k。値が非負のとき単調)
    //
    //  max_right(l, g) : g(op(a[l], ..., a[r-1])) == true となる最大の r  (l <= r <= n)
    //                    r == n なら [l, n) 全体で true、そうでなければ a[r] を足すと初めて false
    template <class G>
    int max_right(int l, G g) {  // g にはラムダをそのまま渡せる
        assert(0 <= l && l <= _n);
        assert(g(e()));
        if (l == _n) return _n;
        l += _size;
        for (int i = _log; i >= 1; i--) push(l >> i);
        S sm = e();
        do {
            while (l % 2 == 0) l >>= 1;
            if (!g(op(sm, d[l]))) {
                while (l < _size) {
                    push(l);
                    l = 2 * l;
                    if (g(op(sm, d[l]))) {
                        sm = op(sm, d[l]);
                        l++;
                    }
                }
                return l - _size;
            }
            sm = op(sm, d[l]);
            l++;
        } while ((l & -l) != l);
        return _n;
    }
    //  min_left(r, g) : g(op(a[l], ..., a[r-1])) == true となる最小の l  (0 <= l <= r)
    //                   l == 0 なら [0, r) 全体で true、そうでなければ a[l-1] を足すと初めて false
    template <class G>
    int min_left(int r, G g) {
        assert(0 <= r && r <= _n);
        assert(g(e()));
        if (r == 0) return 0;
        r += _size;
        for (int i = _log; i >= 1; i--) push((r - 1) >> i);
        S sm = e();
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!g(op(d[r], sm))) {
                while (r < _size) {
                    push(r);
                    r = 2 * r + 1;
                    if (g(op(d[r], sm))) {
                        sm = op(d[r], sm);
                        r--;
                    }
                }
                return r + 1 - _size;
            }
            sm = op(d[r], sm);
        } while ((r & -r) != r);
        return 0;
    }

    // ---- 値で探す (max_right / min_left の よく使う形) ----
    //  x は ±INF 以外 (ふつうの値) を渡す。より大きい / より小さい は x+1 / x-1 を渡す。
    //  [l, n) で a[i] >= x となる最小の i。なければ n
    int find_first_geq(int l, ll x) { return max_right(l, [x](const S& s) { return s.mx < x; }); }
    //  [l, n) で a[i] <= x となる最小の i。なければ n
    int find_first_leq(int l, ll x) { return max_right(l, [x](const S& s) { return s.mn > x; }); }
    //  [0, r) で a[i] >= x となる最大の i。なければ -1
    int find_last_geq(int r, ll x) { return min_left(r, [x](const S& s) { return s.mx < x; }) - 1; }
    //  [0, r) で a[i] <= x となる最大の i。なければ -1
    int find_last_leq(int r, ll x) { return min_left(r, [x](const S& s) { return s.mn > x; }) - 1; }

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

// ---- 二分探索 (max_right / min_left / find_*) のランダムテスト ----
// 素朴解: 1 つずつ足していき、g が false になる手前で止める (g は単調な述語)
static LazySegTree::S leaf(ll v) { return LazySegTree::S{v, v, v, 1}; }
template <class G>
static int naive_max_right(const std::vector<ll>& a, int l, G g) {
    LazySegTree::S sm = LazySegTree::e();
    int r = l;
    while (r < (int)a.size()) {
        LazySegTree::S nx = LazySegTree::op(sm, leaf(a[r]));
        if (!g(nx)) break;
        sm = nx;
        r++;
    }
    return r;
}
template <class G>
static int naive_min_left(const std::vector<ll>& a, int r, G g) {
    LazySegTree::S sm = LazySegTree::e();
    int l = r;
    while (l > 0) {
        LazySegTree::S nx = LazySegTree::op(leaf(a[l - 1]), sm);
        if (!g(nx)) break;
        sm = nx;
        l--;
    }
    return l;
}

// 区間更新をランダムに 1 回 (nonneg = true のときは値が負にならない更新だけ)
static void random_update(std::mt19937_64& rng, std::vector<ll>& a, LazySegTree& seg, bool nonneg) {
    const ll VAL = 1000;
    int n = (int)a.size();
    ll big = 0;
    for (ll x : a) big = std::max(big, absll(x));
    if (big > 1000000000000LL) {  // 値が発散しないようにガード
        ll v = rnd_range(rng, nonneg ? 0 : -VAL, VAL);
        for (auto& x : a) x = v;
        seg.range_assign(0, n, v);
    }
    int l = (int)rnd_range(rng, 0, n), r = (int)rnd_range(rng, 0, n);
    if (l > r) std::swap(l, r);
    int type = (int)rnd_range(rng, 0, 3);
    if (type == 0) {
        ll v = rnd_range(rng, nonneg ? 0 : -VAL, VAL);
        for (int i = l; i < r; i++) a[i] += v;
        seg.range_add(l, r, v);
    } else if (type == 1) {
        ll v = rnd_range(rng, nonneg ? 0 : -VAL, VAL);
        for (int i = l; i < r; i++) a[i] = v;
        seg.range_assign(l, r, v);
    } else if (type == 2) {
        ll b = rnd_range(rng, nonneg ? 0 : -3, 3), c = rnd_range(rng, nonneg ? 0 : -VAL, VAL);
        for (int i = l; i < r; i++) a[i] = b * a[i] + c;
        seg.range_affine(l, r, b, c);
    } else {
        int p = (int)rnd_range(rng, 0, n - 1);
        ll x = rnd_range(rng, nonneg ? 0 : -VAL, VAL);
        a[p] = x;
        seg.point_set(p, x);
    }
}

// 述語を 1 つランダムに選び、位置 pos で max_right (right = true) か min_left を 1 回だけ突き合わせる
static void check_bsearch_one(std::mt19937_64& rng, const std::vector<ll>& a, LazySegTree& seg, bool nonneg, bool right, int pos) {
    using S = LazySegTree::S;
    const int n = (int)a.size();
    ll lo = a[0], hi = a[0];
    for (ll v : a) lo = std::min(lo, v), hi = std::max(hi, v);
    // x は「配列の値ちょうど」「その ±1」「範囲の中」「範囲の外」が出るように選ぶ
    ll x;
    int tx = (int)rnd_range(rng, 0, 3);
    if (tx == 0) x = a[rnd_range(rng, 0, n - 1)];
    else if (tx == 1) x = a[rnd_range(rng, 0, n - 1)] + rnd_range(rng, -1, 1);
    else if (tx == 2) x = rnd_range(rng, lo - 2, hi + 2);
    else x = rnd_range(rng, 0, 1) ? lo - 5 : hi + 5;
    const ll k = rnd_range(rng, 0, n + 1);
    const ll s0 = rnd_range(rng, 0, 3000);
    auto g_mx = [&](const S& s) { return s.mx < x; };
    auto g_mn = [&](const S& s) { return s.mn > x; };
    auto g_len = [&](const S& s) { return s.len <= k; };
    auto g_mix = [&](const S& s) { return s.mx < x && s.len <= k; };
    auto g_sum = [&](const S& s) { return s.sum <= s0; };  // 値が非負のときだけ単調
    const int type = (int)rnd_range(rng, 0, nonneg ? 6 : 5);
    if (right) {
        const int l = pos;
        if (type == 0) assert(seg.max_right(l, g_mx) == naive_max_right(a, l, g_mx));
        if (type == 1) assert(seg.max_right(l, g_mn) == naive_max_right(a, l, g_mn));
        if (type == 2) assert(seg.max_right(l, g_len) == naive_max_right(a, l, g_len));
        if (type == 3) assert(seg.max_right(l, g_mix) == naive_max_right(a, l, g_mix));
        if (type == 4) assert(seg.find_first_geq(l, x) == naive_max_right(a, l, g_mx));
        if (type == 5) assert(seg.find_first_leq(l, x) == naive_max_right(a, l, g_mn));
        if (type == 6) assert(seg.max_right(l, g_sum) == naive_max_right(a, l, g_sum));
    } else {
        const int r = pos;
        if (type == 0) assert(seg.min_left(r, g_mx) == naive_min_left(a, r, g_mx));
        if (type == 1) assert(seg.min_left(r, g_mn) == naive_min_left(a, r, g_mn));
        if (type == 2) assert(seg.min_left(r, g_len) == naive_min_left(a, r, g_len));
        if (type == 3) assert(seg.min_left(r, g_mix) == naive_min_left(a, r, g_mix));
        if (type == 4) assert(seg.find_last_geq(r, x) == naive_min_left(a, r, g_mx) - 1);
        if (type == 5) assert(seg.find_last_leq(r, x) == naive_min_left(a, r, g_mn) - 1);
        if (type == 6) assert(seg.min_left(r, g_sum) == naive_min_left(a, r, g_sum));
    }
}

// 更新の直後に 1 回だけ (遅延が残った状態で) 探索する、を繰り返したあと、全位置をばらばらの順で突き合わせる
static long long check_round(std::mt19937_64& rng, std::vector<ll>& a, LazySegTree& seg, bool nonneg, int steps) {
    const int n = (int)a.size();
    long long checks = 0;
    for (int t = 0; t < steps; t++) {
        random_update(rng, a, seg, nonneg);
        bool right = rnd_range(rng, 0, 1) == 1;
        check_bsearch_one(rng, a, seg, nonneg, right, (int)rnd_range(rng, 0, n));
        checks++;
    }
    std::vector<int> order(n + 1);
    for (int i = 0; i <= n; i++) order[i] = i;
    std::shuffle(order.begin(), order.end(), rng);
    for (int p : order) check_bsearch_one(rng, a, seg, nonneg, true, p), checks++;
    std::shuffle(order.begin(), order.end(), rng);
    for (int p : order) check_bsearch_one(rng, a, seg, nonneg, false, p), checks++;
    // 二分探索の途中の push で壊れていないか、全要素を確認
    for (int i = 0; i < n; i++) assert(seg.point_get(i) == a[i]);
    return checks;
}

static void random_test_bsearch() {
    std::mt19937_64 rng(20261019ULL);
    const ll VAL = 1000;
    long long checks = 0;

    for (int rep = 0; rep < 400; rep++) {
        const bool nonneg = rep % 2 == 1;  // 半分は値を非負に保って「和」の二分探索も試す
        int n = (int)rnd_range(rng, 1, 40);
        std::vector<ll> a(n);
        // 値の種類を少なくして同じ値・境界ちょうどが出やすくする
        ll range = rnd_range(rng, 0, 1) ? 3 : VAL;
        for (auto& x : a) x = rnd_range(rng, nonneg ? 0 : -range, range);
        LazySegTree seg(a);
        for (int round = 0; round < 10; round++) checks += check_round(rng, a, seg, nonneg, 40);
    }

    // 2 のべき乗ちょうど / ±1 のサイズ (末尾の詰め物の葉をまたぐ探索)
    for (int n : {1, 2, 3, 7, 8, 9, 15, 16, 17, 63, 64, 65}) {
        for (int rep = 0; rep < 20; rep++) {
            const bool nonneg = rep % 2 == 1;
            std::vector<ll> a(n);
            for (auto& x : a) x = rnd_range(rng, nonneg ? 0 : -VAL, VAL);
            LazySegTree seg(a);
            for (int round = 0; round < 5; round++) checks += check_round(rng, a, seg, nonneg, 20);
        }
    }

    // 全体で true になる述語 (r == n / l == 0 を返す) と、先頭ですぐ false になる述語
    {
        LazySegTree seg(std::vector<ll>{5, 1, 7});
        auto all_true = [](const LazySegTree::S&) { return true; };
        auto len0 = [](const LazySegTree::S& s) { return s.len == 0; };
        for (int l = 0; l <= 3; l++) assert(seg.max_right(l, all_true) == 3 && seg.max_right(l, len0) == l);
        for (int r = 0; r <= 3; r++) assert(seg.min_left(r, all_true) == 0 && seg.min_left(r, len0) == r);
        assert(seg.find_first_geq(0, 100) == 3 && seg.find_last_leq(3, -100) == -1);
    }

    std::cout << "[lazy_segtree] random test (max_right / min_left / find_*): OK (" << checks << " 回)" << std::endl;
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

    // -------------------------------------------------------------------------
    //  二分探索 (max_right / min_left)
    // -------------------------------------------------------------------------
    std::cout << "\n----- 二分探索 (max_right / min_left / find_*) -----\n";
    using S = LazySegTree::S;  // ラムダの引数の型 (sum, mn, mx, len を持つ)
    LazySegTree t(std::vector<ll>{3, 1, 4, 1, 5, 9, 2, 6});
    dump("配列 t             ", t);

    // 【max_right(l, g)】l から右へ伸ばし、g([l, r) の S) が true のままでいられる最大の r を返す
    std::cout << "max_right(0, mx<5)    = " << t.max_right(0, [](const S& s) { return s.mx < 5; })
              << "  ([0, r) の最大値が 5 未満のままでいられる最大の r = 最初に 5 以上が現れる位置)\n";
    std::cout << "max_right(2, sum<=10) = " << t.max_right(2, [](const S& s) { return s.sum <= 10; })
              << "  ([2, r) の和が 10 以下のままでいられる最大の r。値が非負なので単調)\n";
    std::cout << "max_right(1, len<=3)  = " << t.max_right(1, [](const S& s) { return s.len <= 3; })
              << "  (len も使える。[1, r) の長さが 3 以下 → r = 4)\n";
    std::cout << "max_right(0, mx<100)  = " << t.max_right(0, [](const S& s) { return s.mx < 100; })
              << "  (最後まで true なら n を返す)\n";

    // 【min_left(r, g)】r から左へ伸ばし、g([l, r) の S) が true のままでいられる最小の l を返す
    std::cout << "min_left(8, mn>1)     = " << t.min_left(8, [](const S& s) { return s.mn > 1; })
              << "  ([l, 8) の最小値が 1 より大きいままでいられる最小の l = 最後の 1 以下の位置 + 1)\n";
    std::cout << "min_left(5, sum<=6)   = " << t.min_left(5, [](const S& s) { return s.sum <= 6; })
              << "  ([l, 5) の和が 6 以下のままでいられる最小の l)\n";
    std::cout << "min_left(3, mn>-100)  = " << t.min_left(3, [](const S& s) { return s.mn > -100; })
              << "  (最後まで true なら 0 を返す)\n";

    // 【find_first_geq / find_first_leq / find_last_geq / find_last_leq】値で位置を探す
    std::cout << "find_first_geq(0, 5)  = " << t.find_first_geq(0, 5) << "  ([0, n) で 5 以上になる最初の位置)\n";
    std::cout << "find_first_geq(6, 7)  = " << t.find_first_geq(6, 7) << "  (見つからなければ n)\n";
    std::cout << "find_first_leq(2, 1)  = " << t.find_first_leq(2, 1) << "  ([2, n) で 1 以下になる最初の位置)\n";
    std::cout << "find_last_geq(8, 6)   = " << t.find_last_geq(8, 6) << "  ([0, 8) で 6 以上になる最後の位置)\n";
    std::cout << "find_last_leq(5, 1)   = " << t.find_last_leq(5, 1) << "  ([0, 5) で 1 以下になる最後の位置)\n";
    std::cout << "find_last_leq(3, 0)   = " << t.find_last_leq(3, 0) << "  (見つからなければ -1)\n";
    std::cout << "find_first_geq(0, 4+1)= " << t.find_first_geq(0, 4 + 1) << "  (「より大きい」は x+1 を渡す: 4 より大きい最初の位置)\n";

    // 【使い方の例 1】昇順の列で lower_bound の代わりに使い、値が [3, 5] の位置だけに +1 (ABC389 F の形)
    LazySegTree u(std::vector<ll>{1, 2, 3, 4, 5, 6});
    int lo = u.find_first_geq(0, 3), hi = u.find_first_geq(0, 5 + 1);  // 値 >= 3 の先頭, 値 > 5 の先頭
    u.range_add(lo, hi, 1);  // 昇順は崩れない
    std::cout << "例1: 値が [3,5] の位置 [" << lo << ", " << hi << ") に +1 ";
    dump("", u);

    // 【使い方の例 2】0/1 の列で「k 番目 (1 始まり) の 1 の位置」= 和が k 未満でいられる最大の r
    LazySegTree w(std::vector<ll>{0, 1, 1, 0, 1, 0, 1});
    ll k = 3;
    int pos = w.max_right(0, [&](const S& s) { return s.sum < k; });
    std::cout << "例2: {0,1,1,0,1,0,1} の 3 番目の 1 の位置 = " << pos << "  (k 番目がなければ n)\n";
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    tester::random_test();
    tester::random_test_bsearch();
    demo();
    return 0;
}
