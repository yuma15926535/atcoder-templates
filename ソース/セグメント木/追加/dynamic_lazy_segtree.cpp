// =============================================================================
//  ⑤ Dynamic Lazy Segment Tree (動的遅延伝播セグメント木)  --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  値の型      : long long 固定
//  座標        : long long の全範囲 [LLONG_MIN, LLONG_MAX]。範囲指定も座標圧縮も不要。
//                一度も触っていない座標の値は 0。
//  区間の流儀  : 半開区間 [l, r)
//                ※ 半開なので区間更新/取得で LLONG_MAX 自身は含められない。
//                   その 1 点は point_set / point_add / point_get を使うこと (全域 OK)。
//  依存        : 標準ライブラリのみ (C++17)
//
//  更新はすべて「アフィン変換 x -> b*x + c」に統一:
//      range_add(l, r, v)    = range_affine(l, r, 1, v)
//      range_assign(l, r, v) = range_affine(l, r, 0, v)
//  各ノードは (sum, min, max) を持つので 1 本で和・最小・最大が取れる。
//  b < 0 のときは min/max が入れ替わる (apply_node / query で処理済み)。
//
// -----------------------------------------------------------------------------
//  ★ 設計の要点 (ここが「動的」セグメント木として効く部分)
//
//  (1) 未生成の子は「一様値スロット」で持つ
//      各ノードは子のインデックス lc / rc と、その子が未生成 (-1) のときの
//      一様値 uniL / uniR を持つ。「その子の区間は全部この値」という意味。
//      → 区間更新で丸ごと覆われた部分木は、ノードを作らずに一様値を書き換えるだけ。
//         (range_assign / range_add / range_affine はどれも一様区間を一様のまま保つ)
//      → push は子を実体化しない。実体化するのは「分割が必要になった枝」だけ。
//         つまり「左だけ必要なら左だけ生成」。両方生成は絶対にしない。
//
//  (2) 取得 (query / point_get) はノードを一切生成せず、木も書き換えない (const)
//      遅延を実ノードへ押し下げる代わりに、合成した (b, c) を再帰の引数として
//      持ち回る。未生成の子に当たったら「一様値に (b, c) を適用した値」を返す。
//
//  (3) point_get は range_query 経由ではなく専用の 1 本道降下
//      根から目的の葉まで、各段で必要な子 1 つだけを見る。O(木の深さ)。
//
// -----------------------------------------------------------------------------
//  計算量 (木の深さを D とする。既定の TOP = 63 なので D = 64)
//    point_get                        : O(D)。1 本道。生成ノード 0 個。
//    range_query / range_sum/min/max  : O(D)。訪れるノードは高々 2D 個。生成ノード 0 個。
//    all_query                        : O(1)
//    point_set / point_add            : O(D)
//    range_add / range_assign / range_affine
//                                     : O(D)。押し下げは 1 段 O(1)。
//
//  メモリ  : 1 ノード 64 バイト。
//            ノードが作られるのは「区間更新の左右の境界が通る経路」だけで、
//            丸ごと覆われる部分は一様値スロットで済むため 0 個。
//            → 区間更新 1 回あたり最悪 2D ≒ 126 ノード、実際にはもっと少ない。
//            取得系は 0 個。
//            座標圧縮できるなら ② lazy_segtree.cpp のほうが圧倒的に速くて軽い。
//
//  オーバーフロー注意 (long long のまま扱うので仕様上避けられないもの)
//    ・sum : 長さ 10^18 級の区間に range_add / range_assign すると当然あふれる。
//            全域 (2^64 要素) を 1 で埋めれば sum は意味を持たない。
//    ・mul_len(値, 長さ) は符号なしで掛けてから long long に戻している。
//      未定義動作は起こさないが、値としてはラップアラウンドする。
//    ・アフィンを何度も重ねると b が指数的に大きくなるので、
//      |b| >= 2 のアフィンを繰り返す使い方では値の発散に注意。
// =============================================================================

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

using ll = long long;
using ull = unsigned long long;

struct DynamicLazySegTree {
    // ★ 扱う座標範囲: [-2^TOP, 2^TOP)。既定の TOP = 63 で long long 全域。
    //    (座標の絶対値の上限が確実に分かっている場合だけ下げてよい。
    //     下げるとノード数も実行時間もほぼ (TOP+1)/64 倍になる。既定は全域のまま)
    static constexpr int TOP = 63;

    static constexpr ll INF = std::numeric_limits<ll>::max();  // 空区間の min として返る値

    struct Agg {
        ll sum, mn, mx;
    };

    DynamicLazySegTree() = default;
    void reserve(size_t n) { t.reserve(n); }  // 事前確保しておくと速い (任意)

    // =========================================================================
    //  区間更新 : [l, r)
    // =========================================================================
    void range_affine(ll l, ll r, ll b, ll c) {  // a[i] = b*a[i] + c
        if (l >= r) return;
        if (b == 1 && c == 0) return;  // 恒等変換
        ull ql = enc(l), qr = enc(r - 1);
        assert(in_range(ql) && in_range(qr));
        for (int s = 0; s < 2; s++) {
            Ref res = apply_rec(root_ref(s), root_lo(s), TOP, ql, qr, b, c);
            store_root(s, res);
        }
    }
    void range_add(ll l, ll r, ll v) { range_affine(l, r, 1, v); }     // a[i] += v
    void range_assign(ll l, ll r, ll v) { range_affine(l, r, 0, v); }  // a[i] = v

    // =========================================================================
    //  点操作 (LLONG_MAX / LLONG_MIN も含め全座標で使える)
    // =========================================================================
    void point_set(ll key, ll x) { point_affine(key, 0, x); }  // a[key] = x
    void point_add(ll key, ll v) { point_affine(key, 1, v); }  // a[key] += v

    // a[key] を 1 本道で降りて返す。ノードは生成しないし木も書き換えない。
    ll point_get(ll key) const {
        ull k = enc(key);
        assert(in_range(k));
        int s = side_of(k);
        ull lo = root_lo(s);
        int level = TOP;
        Ref cur = root_ref(s);
        ll b = 1, c = 0;  // 祖先から持ち越した遅延
        while (true) {
            if (cur.idx == -1) return b * cur.uni + c;  // 未生成 = 一様
            const Node& nd = t[cur.idx];
            if (level == 0) return b * nd.sum + c;  // 葉 (長さ 1 なので sum が値)
            ll nb = b * nd.b, nc = b * nd.c + c;    // 自分の遅延を合成してから降りる
            int cl = level - 1;
            ull mid = lo + (1ULL << cl);
            if (k < mid) {
                cur = Ref{nd.lc, nd.uniL};
            } else {
                cur = Ref{nd.rc, nd.uniR};
                lo = mid;
            }
            b = nb;
            c = nc;
            level = cl;
        }
    }

    // =========================================================================
    //  区間取得 : [l, r)   ---- すべて const / ノード生成 0 個
    // =========================================================================
    Agg range_query(ll l, ll r) const {  // sum/min/max をまとめて 1 回で取る
        if (l >= r) return Agg{0, INF, -INF};
        ull ql = enc(l), qr = enc(r - 1);
        assert(in_range(ql) && in_range(qr));
        bool hit0 = (ql < HALF);   // 左の根と交差するか
        bool hit1 = (qr >= HALF);  // 右の根と交差するか
        if (hit0 && hit1) {
            Agg a = query_rec(root_ref(0), root_lo(0), TOP, ql, qr, 1, 0);
            Agg b = query_rec(root_ref(1), root_lo(1), TOP, ql, qr, 1, 0);
            return merge(a, b);
        }
        if (hit0) return query_rec(root_ref(0), root_lo(0), TOP, ql, qr, 1, 0);
        return query_rec(root_ref(1), root_lo(1), TOP, ql, qr, 1, 0);
    }
    ll range_sum(ll l, ll r) const { return range_query(l, r).sum; }  // 空区間は 0
    ll range_min(ll l, ll r) const { return range_query(l, r).mn; }   // 空区間は  INF
    ll range_max(ll l, ll r) const { return range_query(l, r).mx; }   // 空区間は -INF

    // 全座標 (LLONG_MAX 含む) の集約。長さ 2^64 なので sum はまず溢れる点に注意。
    Agg all_query() const {
        return merge(agg_of(root_ref(0), TOP), agg_of(root_ref(1), TOP));
    }

    size_t node_count() const { return t.size(); }

   private:
    static constexpr ull HALF = 1ULL << 63;
    static constexpr ull BASE = HALF - (1ULL << TOP);  // 左の根の開始位置 (TOP=63 なら 0)
    static bool in_range(ull u) { return (((u - BASE) >> TOP) >> 1) == 0; }
    static ull enc(ll x) { return (ull)x ^ HALF; }
    static int side_of(ull k) { return k < HALF ? 0 : 1; }
    static ull root_lo(int s) { return s == 0 ? BASE : HALF; }
    // 値 * 長さ を符号なしで計算 (あふれてもラップするだけで未定義動作にはならない)
    static ll mul_len(ll v, ull len) { return (ll)((ull)v * len); }
    static Agg merge(const Agg& a, const Agg& b) {
        return Agg{a.sum + b.sum, std::min(a.mn, b.mn), std::max(a.mx, b.mx)};
    }

    struct Node {
        ll sum, mn, mx;  // この部分木の集約 (自分の遅延は反映済み)
        ll b, c;         // 子に未反映の遅延 x -> b*x + c
        ll uniL;         // lc == -1 のとき、左の子の区間はすべてこの値
        ll uniR;         // rc == -1 のとき、右の子の区間はすべてこの値
        int lc, rc;
    };
    std::vector<Node> t;
    int rootIdx[2] = {-1, -1};
    ll rootUni[2] = {0, 0};  // 根が未生成のときの一様値 (初期状態は全部 0)

    // 子への参照: idx == -1 なら「その区間は全部 uni」
    struct Ref {
        int idx;
        ll uni;
    };
    Ref root_ref(int s) const { return Ref{rootIdx[s], rootUni[s]}; }
    void store_root(int s, Ref r) {
        rootIdx[s] = r.idx;
        if (r.idx == -1) rootUni[s] = r.uni;
    }

    static Agg uniform_agg(ll v, int level) {  // 「v が 2^level 個」
        return Agg{mul_len(v, 1ULL << level), v, v};
    }
    Agg agg_of(Ref r, int level) const {
        if (r.idx == -1) return uniform_agg(r.uni, level);
        const Node& n = t[r.idx];
        return Agg{n.sum, n.mn, n.mx};
    }

    int new_uniform(ll v, int level) {  // 「v が 2^level 個」のノードを 1 つ作る
        assert(level >= 1);
        t.push_back(Node{mul_len(v, 1ULL << level), v, v, 1, 0, v, v, -1, -1});
        return (int)t.size() - 1;
    }

    void apply_node(int p, ll b, ll c, int level) {
        t[p].sum = b * t[p].sum + mul_len(c, 1ULL << level);
        ll lo = b * t[p].mn + c, hi = b * t[p].mx + c;
        if (b < 0) std::swap(lo, hi);  // 負の傾きだと min/max が入れ替わる
        t[p].mn = lo;
        t[p].mx = hi;
        t[p].b = b * t[p].b;
        t[p].c = b * t[p].c + c;
    }

    // 遅延を子へ降ろす。★ ここでノードは一切生成しない ★
    //   ・実体のある子 : apply_node
    //   ・未生成の子   : 一様値を b*uni + c に更新するだけ
    void push(int p, int level) {
        if (t[p].b == 1 && t[p].c == 0) return;  // 押し下げるものが無い
        ll b = t[p].b, c = t[p].c;
        int cl = level - 1;
        if (t[p].lc != -1) apply_node(t[p].lc, b, c, cl);
        else t[p].uniL = b * t[p].uniL + c;
        if (t[p].rc != -1) apply_node(t[p].rc, b, c, cl);
        else t[p].uniR = b * t[p].uniR + c;
        t[p].b = 1;
        t[p].c = 0;
    }

    void pull(int p, int level) {
        int cl = level - 1;
        Agg a = agg_of(Ref{t[p].lc, t[p].uniL}, cl);
        Agg d = agg_of(Ref{t[p].rc, t[p].uniR}, cl);
        t[p].sum = a.sum + d.sum;
        t[p].mn = std::min(a.mn, d.mn);
        t[p].mx = std::max(a.mx, d.mx);
    }

    void store_left(int p, Ref r) {
        t[p].lc = r.idx;
        if (r.idx == -1) t[p].uniL = r.uni;
    }
    void store_right(int p, Ref r) {
        t[p].rc = r.idx;
        if (r.idx == -1) t[p].uniR = r.uni;
    }

    void point_affine(ll key, ll b, ll c) {
        ull k = enc(key);
        assert(in_range(k));
        int s = side_of(k);
        Ref res = apply_rec(root_ref(s), root_lo(s), TOP, k, k, b, c);
        store_root(s, res);
    }

    Ref apply_rec(Ref cur, ull lo, int level, ull ql, ull qr, ll b, ll c) {
        ull hi = lo + ((1ULL << level) - 1);
        if (qr < lo || hi < ql) return cur;  // 交差なし
        if (cur.idx == -1) {
            ll nv = b * cur.uni + c;
            if (nv == cur.uni) return cur;  // 何も変わらない (例: +0, 同じ値の assign)
            if (ql <= lo && hi <= qr) return Ref{-1, nv};  // 一様のまま。ノード生成なし
        } else if (ql <= lo && hi <= qr) {
            apply_node(cur.idx, b, c, level);
            return cur;
        }
        int p = cur.idx;
        if (p == -1) p = new_uniform(cur.uni, level);  // 分割が要るのでここだけ実体化
        push(p, level);
        int cl = level - 1;
        ull mid = lo + (1ULL << cl);
        Ref L = apply_rec(Ref{t[p].lc, t[p].uniL}, lo, cl, ql, qr, b, c);
        store_left(p, L);
        Ref R = apply_rec(Ref{t[p].rc, t[p].uniR}, mid, cl, ql, qr, b, c);
        store_right(p, R);
        pull(p, level);
        return Ref{p, 0};
    }

    // 遅延 (b, c) を引数で持ち回るので木を書き換えない。
    // 前提: [lo, hi] と [ql, qr] は交差する。
    Agg query_rec(Ref cur, ull lo, int level, ull ql, ull qr, ll b, ll c) const {
        ull hi = lo + ((1ULL << level) - 1);
        if (cur.idx == -1) {  // 一様なので交差長ぶん
            ll v = b * cur.uni + c;
            ull ilen = std::min(hi, qr) - std::max(lo, ql) + 1;
            return Agg{mul_len(v, ilen), v, v};
        }
        const Node& nd = t[cur.idx];
        if (ql <= lo && hi <= qr) {  // 丸ごと入る
            ll s = b * nd.sum + mul_len(c, 1ULL << level);
            ll mn = b * nd.mn + c, mx = b * nd.mx + c;
            if (b < 0) std::swap(mn, mx);
            return Agg{s, mn, mx};
        }
        ll nb = b * nd.b, nc = b * nd.c + c;  // 自分の遅延を合成してから子へ
        int cl = level - 1;
        ull mid = lo + (1ULL << cl);
        Ref lref{nd.lc, nd.uniL}, rref{nd.rc, nd.uniR};
        bool hitL = !(qr < lo || mid - 1 < ql);
        bool hitR = !(qr < mid || hi < ql);
        if (hitL && hitR) {
            Agg a = query_rec(lref, lo, cl, ql, qr, nb, nc);
            Agg d = query_rec(rref, mid, cl, ql, qr, nb, nc);
            return merge(a, d);
        }
        if (hitL) return query_rec(lref, lo, cl, ql, qr, nb, nc);
        return query_rec(rref, mid, cl, ql, qr, nb, nc);
    }
};

// =============================================================================
//  ここから下は「ランダムテスト」と「使用例デモ」。提出時は丸ごと消してよい。
// =============================================================================

namespace tester {

static ll rnd_range(std::mt19937_64& rng, ll lo, ll hi) {
    return lo + (ll)(rng() % (ull)(hi - lo + 1));
}
static ll absll(ll x) { return x < 0 ? -x : x; }

// -----------------------------------------------------------------------------
//  窓テスト: 座標 base + i (i = 0..W-1) を素朴配列と突き合わせる。
//  更新は内側 [PAD, W-PAD) のみ。外側 PAD は常に 0 なので「未設定 = 0」も検証できる。
//  取得のたびに node_count が増えていないことも確認する。
// -----------------------------------------------------------------------------
static void test_window(std::mt19937_64& rng, ll base, int W, int iters) {
    const int PAD = 4;
    assert(W > 2 * PAD);
    std::vector<ll> a(W, 0);
    DynamicLazySegTree seg;
    const ll VAL = 1000;

    for (int t = 0; t < iters; t++) {
        // 値が発散しないようにガード
        ll big = 0;
        for (ll x : a) big = std::max(big, absll(x));
        if (big > 1000000000000LL) {
            ll v = rnd_range(rng, -VAL, VAL);
            for (int i = PAD; i < W - PAD; i++) a[i] = v;
            seg.range_assign(base + PAD, base + W - PAD, v);
        }

        int il = (int)rnd_range(rng, PAD, W - PAD), ir = (int)rnd_range(rng, PAD, W - PAD);
        if (il > ir) std::swap(il, ir);
        int p = (int)rnd_range(rng, PAD, W - PAD - 1);
        int type = (int)rnd_range(rng, 0, 7);

        if (type == 0) {  // range_add
            ll v = rnd_range(rng, -VAL, VAL);
            for (int i = il; i < ir; i++) a[i] += v;
            seg.range_add(base + il, base + ir, v);
        } else if (type == 1) {  // range_assign
            ll v = rnd_range(rng, -VAL, VAL);
            for (int i = il; i < ir; i++) a[i] = v;
            seg.range_assign(base + il, base + ir, v);
        } else if (type == 2) {  // range_affine
            ll b = rnd_range(rng, -3, 3), c = rnd_range(rng, -VAL, VAL);
            for (int i = il; i < ir; i++) a[i] = b * a[i] + c;
            seg.range_affine(base + il, base + ir, b, c);
        } else if (type == 3) {  // point_set
            ll x = rnd_range(rng, -VAL, VAL);
            a[p] = x;
            seg.point_set(base + p, x);
        } else if (type == 4) {  // point_add
            ll v = rnd_range(rng, -VAL, VAL);
            a[p] += v;
            seg.point_add(base + p, v);
        } else {  // 取得 (クエリ範囲は外側 PAD も含めてよい)
            size_t before = seg.node_count();
            int ql = (int)rnd_range(rng, 0, W), qr = (int)rnd_range(rng, 0, W);
            if (ql > qr) std::swap(ql, qr);
            ll s = 0, mn = DynamicLazySegTree::INF, mx = -DynamicLazySegTree::INF;
            for (int i = ql; i < qr; i++) {
                s += a[i];
                mn = std::min(mn, a[i]);
                mx = std::max(mx, a[i]);
            }
            DynamicLazySegTree::Agg q = seg.range_query(base + ql, base + qr);
            assert(q.sum == s && q.mn == mn && q.mx == mx);
            assert(seg.range_sum(base + ql, base + qr) == s);
            assert(seg.range_min(base + ql, base + qr) == mn);
            assert(seg.range_max(base + ql, base + qr) == mx);
            assert(seg.node_count() == before);  // ★ 取得でノードが増えないこと
        }
        size_t before = seg.node_count();
        int gp = (int)rnd_range(rng, 0, W - 1);
        assert(seg.point_get(base + gp) == a[gp]);
        assert(seg.node_count() == before);  // ★ point_get でノードが増えないこと
    }

    // 最後に全区間を突き合わせ
    size_t before = seg.node_count();
    for (int l = 0; l <= W; l++) {
        for (int r = l; r <= W; r++) {
            ll s = 0, mn = DynamicLazySegTree::INF, mx = -DynamicLazySegTree::INF;
            for (int i = l; i < r; i++) {
                s += a[i];
                mn = std::min(mn, a[i]);
                mx = std::max(mx, a[i]);
            }
            DynamicLazySegTree::Agg q = seg.range_query(base + l, base + r);
            assert(q.sum == s && q.mn == mn && q.mx == mx);
        }
    }
    for (int i = 0; i < W; i++) assert(seg.point_get(base + i) == a[i]);
    assert(seg.node_count() == before);  // ★ 全区間なめてもノードは増えない
}

// -----------------------------------------------------------------------------
//  操作の組み合わせ (add→assign, assign→affine, affine→add, point_set→... など)
// -----------------------------------------------------------------------------
static void test_op_pairs(std::mt19937_64& rng) {
    enum Op { ADD, ASSIGN, AFFINE, PSET, PADD, OPCNT };
    for (int o1 = 0; o1 < OPCNT; o1++) {
        for (int o2 = 0; o2 < OPCNT; o2++) {
            for (int rep = 0; rep < 40; rep++) {
                int n = (int)rnd_range(rng, 1, 12);
                ll base = rnd_range(rng, -50, 50);
                std::vector<ll> a(n);
                DynamicLazySegTree seg;
                for (int i = 0; i < n; i++) {
                    a[i] = rnd_range(rng, -5, 5);
                    seg.point_set(base + i, a[i]);
                }
                int ops[2] = {o1, o2};
                for (int k = 0; k < 2; k++) {
                    int l = (int)rnd_range(rng, 0, n), r = (int)rnd_range(rng, 0, n);
                    if (l > r) std::swap(l, r);
                    ll x = rnd_range(rng, -6, 6);
                    if (ops[k] == ADD) {
                        for (int i = l; i < r; i++) a[i] += x;
                        seg.range_add(base + l, base + r, x);
                    } else if (ops[k] == ASSIGN) {
                        for (int i = l; i < r; i++) a[i] = x;
                        seg.range_assign(base + l, base + r, x);
                    } else if (ops[k] == AFFINE) {
                        ll b = rnd_range(rng, -3, 3);
                        for (int i = l; i < r; i++) a[i] = b * a[i] + x;
                        seg.range_affine(base + l, base + r, b, x);
                    } else if (ops[k] == PSET) {
                        int p = (int)rnd_range(rng, 0, n - 1);
                        a[p] = x;
                        seg.point_set(base + p, x);
                    } else {
                        int p = (int)rnd_range(rng, 0, n - 1);
                        a[p] += x;
                        seg.point_add(base + p, x);
                    }
                    for (int i = 0; i < n; i++) assert(seg.point_get(base + i) == a[i]);
                    for (int l2 = 0; l2 <= n; l2++) {
                        for (int r2 = l2; r2 <= n; r2++) {
                            ll s = 0, mn = DynamicLazySegTree::INF, mx = -DynamicLazySegTree::INF;
                            for (int i = l2; i < r2; i++) {
                                s += a[i];
                                mn = std::min(mn, a[i]);
                                mx = std::max(mx, a[i]);
                            }
                            DynamicLazySegTree::Agg q = seg.range_query(base + l2, base + r2);
                            assert(q.sum == s && q.mn == mn && q.mx == mx);
                        }
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
//  境界値テスト: LLONG_MIN / LLONG_MAX 近傍と符号のまたぎ
// -----------------------------------------------------------------------------
static void test_boundary() {
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();
    std::vector<ll> keys = {LMIN, LMIN + 1, LMIN + 2, -3, -2, -1, 0, 1, 2, 3, LMAX - 2, LMAX - 1, LMAX};

    DynamicLazySegTree seg;
    for (size_t i = 0; i < keys.size(); i++) seg.point_set(keys[i], (ll)i + 1);
    for (size_t i = 0; i < keys.size(); i++) assert(seg.point_get(keys[i]) == (ll)i + 1);

    // 半開区間なので [LMIN, LMAX) は LMAX 以外の全部
    assert(seg.range_sum(LMIN, LMAX) == 78);
    assert(seg.range_max(LMIN, LMAX) == 12);
    assert(seg.point_get(LMAX) == 13);

    // 1 要素だけの区間
    assert(seg.range_sum(LMIN, LMIN + 1) == 1);
    assert(seg.range_sum(LMAX - 1, LMAX) == 12);
    assert(seg.range_sum(0, 1) == 7);

    // 符号をまたぐ区間
    assert(seg.range_sum(-3, 4) == 4 + 5 + 6 + 7 + 8 + 9 + 10);
    assert(seg.range_min(-3, 4) == 4 && seg.range_max(-3, 4) == 10);

    // 空区間 / 逆転区間
    assert(seg.range_sum(5, 5) == 0);
    assert(seg.range_min(5, 5) == DynamicLazySegTree::INF);
    assert(seg.range_max(5, 5) == -DynamicLazySegTree::INF);
    assert(seg.range_sum(10, 3) == 0);

    // 境界をまたぐ更新
    seg.range_add(LMIN, LMAX, 1000);
    for (size_t i = 0; i < keys.size(); i++) {
        ll expect = (keys[i] == LMAX) ? 13 : (ll)i + 1 + 1000;
        assert(seg.point_get(keys[i]) == expect);
    }
    seg.range_affine(LMIN + 1, LMAX - 1, -1, 0);
    assert(seg.point_get(LMIN) == 1001);
    assert(seg.point_get(LMIN + 1) == -1002);
    assert(seg.point_get(LMAX - 2) == -1011);
    assert(seg.point_get(LMAX - 1) == 1012);
    assert(seg.point_get(LMAX) == 13);
    // 未設定座標 (値 0) も -1 倍されているので 0 のまま
    assert(seg.point_get(12345) == -1000);

    // 点操作は LLONG_MAX でも使える
    seg.point_add(LMAX, 7);
    assert(seg.point_get(LMAX) == 20);
}

// 巨大区間に対する更新が「未生成の座標」にも効いていることを確認
static void test_huge() {
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();
    DynamicLazySegTree seg;

    seg.range_assign(LMIN, LMAX, 7);
    assert(seg.point_get(0) == 7);
    assert(seg.point_get(LMIN) == 7);
    assert(seg.point_get(LMAX - 1) == 7);
    assert(seg.point_get(LMAX) == 0);  // 半開なので含まれない
    assert(seg.range_sum(0, 10) == 70);
    assert(seg.range_min(-5, 5) == 7 && seg.range_max(-5, 5) == 7);

    seg.range_add(-1000000000000000000LL, 1000000000000000000LL, 3);
    assert(seg.point_get(0) == 10);
    assert(seg.point_get(-999999999999999999LL) == 10);
    assert(seg.point_get(1000000000000000000LL) == 7);  // 半開なので含まれない
    assert(seg.range_sum(-5, 5) == 100);

    seg.range_affine(-100, 100, -2, 1);  // x -> -2x + 1
    assert(seg.point_get(0) == -19);
    assert(seg.range_min(-100, 100) == -19 && seg.range_max(-100, 100) == -19);
    assert(seg.range_min(-101, 101) == -19 && seg.range_max(-101, 101) == 10);
    assert(seg.range_sum(-100, 100) == -19 * 200);

    seg.point_set(LMAX, 12345);
    assert(seg.point_get(LMAX) == 12345);
    seg.point_add(LMAX, 5);
    assert(seg.point_get(LMAX) == 12350);
    seg.point_set(LMIN, -1);
    assert(seg.point_get(LMIN) == -1);
    assert(seg.range_min(LMIN, LMIN + 1) == -1);

    assert(seg.range_sum(5, 5) == 0);
    assert(seg.range_min(5, 5) == DynamicLazySegTree::INF);
    assert(seg.range_max(5, 5) == -DynamicLazySegTree::INF);
}

// ★ 大量の取得でノードが 1 つも増えないこと (全域ランダム座標)
static void test_query_allocates_nothing(std::mt19937_64& rng) {
    DynamicLazySegTree seg;
    for (int i = 0; i < 300; i++) {
        ll l = (ll)rng(), r = (ll)rng();
        if (l > r) std::swap(l, r);
        int ty = (int)(rng() % 4);
        if (ty == 0) seg.range_add(l, r, rnd_range(rng, -1000, 1000));
        else if (ty == 1) seg.range_assign(l, r, rnd_range(rng, -1000, 1000));
        else if (ty == 2) seg.range_affine(l, r, rnd_range(rng, -2, 2), rnd_range(rng, -100, 100));
        else seg.point_set((ll)rng(), rnd_range(rng, -1000, 1000));
    }
    size_t before = seg.node_count();
    ll acc = 0;
    for (int i = 0; i < 50000; i++) {
        ll l = (ll)rng(), r = (ll)rng();
        if (l > r) std::swap(l, r);
        acc ^= seg.range_sum(l, r) ^ seg.range_min(l, r) ^ seg.range_max(l, r);
        acc ^= seg.point_get((ll)rng());
        acc ^= seg.all_query().sum;
    }
    (void)acc;
    assert(seg.node_count() == before);
}

static void random_test() {
    std::mt19937_64 rng(20260917ULL);
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();

    std::vector<ll> bases = {
        0,                  // 0 起点
        -20,                // 負→正 をまたぐ (符号ビットの境界)
        LMIN,               // 最小値ぴったり
        LMIN + 3,           //
        LMAX - 60,          // 最大値ぎりぎり
        -(1LL << 62),       //
        (1LL << 62),        //
        -(1LL << 40) - 7,   //
        123456789012345LL,  //
    };
    for (ll base : bases) test_window(rng, base, 40, 1200);

    for (int rep = 0; rep < 12; rep++) {
        ll base = (ll)rng();
        if (base > LMAX - 100) base = LMAX - 100;
        if (base < LMIN + 100) base = LMIN + 100;
        test_window(rng, base, 30, 800);
    }

    test_op_pairs(rng);
    test_boundary();
    test_huge();
    test_query_allocates_nothing(rng);
    std::cout << "[dynamic_lazy_segtree] random test: OK" << std::endl;
}

// -----------------------------------------------------------------------------
//  ノード数ベンチ (A: 更新中心 / B: 取得中心 / C: 取得のみ)
// -----------------------------------------------------------------------------
static void bench_nodes() {
    std::cout << "\n----- ノード数ベンチ (座標は 64bit 全域ランダム) -----\n";
    const int Q = 200000;
    const int QA = 50000;  // A はノードを大量に作るので控えめに
    {                      // A. 更新中心
        std::mt19937_64 rng(1);
        DynamicLazySegTree seg;
        for (int i = 0; i < QA; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            int ty = (int)(rng() % 4);
            if (ty == 0) seg.range_add(l, r, (ll)(rng() % 1000));
            else if (ty == 1) seg.range_assign(l, r, (ll)(rng() % 1000));
            else if (ty == 2) seg.range_affine(l, r, 1, (ll)(rng() % 10));
            else seg.point_set((ll)rng(), (ll)(rng() % 1000));
        }
        std::cout << "A 更新 " << QA << " 回           : node_count = " << seg.node_count()
                  << "  (1 更新あたり " << (double)seg.node_count() / QA << " ノード, 約 "
                  << (seg.node_count() * 64 / 1048576) << " MB)\n";
    }
    {  // B. 少数更新 + 大量取得
        std::mt19937_64 rng(2);
        DynamicLazySegTree seg;
        for (int i = 0; i < 1000; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            seg.range_add(l, r, (ll)(rng() % 1000));
        }
        size_t before = seg.node_count();
        ll acc = 0;
        for (int i = 0; i < Q; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            acc ^= seg.range_sum(l, r) ^ seg.range_min(l, r) ^ seg.point_get((ll)rng());
        }
        (void)acc;
        std::cout << "B 更新 1000 + 取得 " << Q << " 回 : node_count = " << before << " → "
                  << seg.node_count() << "  (増加 " << seg.node_count() - before << ")\n";
    }
    {  // C. 取得のみ (少数の区間更新で遅延を残した状態から)
        std::mt19937_64 rng(3);
        DynamicLazySegTree seg;
        for (int i = 0; i < 50; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            seg.range_add(l, r, (ll)(rng() % 1000));
        }
        size_t before = seg.node_count();
        ll acc = 0;
        for (int i = 0; i < Q; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            acc ^= seg.range_sum(l, r) ^ seg.range_min(l, r) ^ seg.range_max(l, r);
            acc ^= seg.point_get((ll)rng());
        }
        (void)acc;
        std::cout << "C 取得のみ " << Q << " 回      : node_count = " << before << " → "
                  << seg.node_count() << "  (増加 " << seg.node_count() - before << ")\n";
    }
}

}  // namespace tester

// =============================================================================
//  使用例デモ : 全機能を 1 箇所で使う。実戦ではここを見れば思い出せる。
// =============================================================================
static void demo() {
    std::cout << "\n===== ⑤ Dynamic Lazy Segment Tree 使用例デモ =====\n";
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();

    // 【構築】サイズも座標範囲も指定しない。未設定の座標はすべて 0。
    DynamicLazySegTree seg;
    // 【reserve】事前確保しておくと速い (任意)
    seg.reserve(1 << 12);

    // 【point_set】1 点代入
    seg.point_set(0, 10);
    seg.point_set(1, 20);
    seg.point_set(1000000000000LL, 30);  // 10^12 のような遠い座標も OK

    // 【point_add】1 点加算
    seg.point_add(1, 5);

    // 【point_get】未設定の座標は 0。専用の 1 本道降下なのでノードは増えない
    std::cout << "point_get(1)      = " << seg.point_get(1) << "\n";
    std::cout << "point_get(999)    = " << seg.point_get(999) << "   (未設定 → 0)\n";

    // 【range_add】[0, 5) に +100
    seg.range_add(0, 5, 100);
    std::cout << "range_add(0,5,+100) 後 point_get(3) = " << seg.point_get(3) << "\n";

    // 【range_assign】[2, 4) を 7 で塗りつぶす
    seg.range_assign(2, 4, 7);

    // 【range_affine】[0, 5) に x -> 2x + 1
    seg.range_affine(0, 5, 2, 1);

    // 【range_affine で符号反転】x -> -x  (min と max が入れ替わるケース)
    seg.range_affine(3, 5, -1, 0);

    std::cout << "座標 0..5 の値     : ";
    for (ll i = 0; i <= 5; i++) std::cout << seg.point_get(i) << (i == 5 ? '\n' : ' ');

    // 【range_sum / range_min / range_max】半開区間 [l, r)
    std::cout << "range_sum(0, 5)   = " << seg.range_sum(0, 5) << "\n";
    std::cout << "range_min(0, 5)   = " << seg.range_min(0, 5) << "\n";
    std::cout << "range_max(0, 5)   = " << seg.range_max(0, 5) << "\n";

    // 【range_query】sum/min/max を 1 回の走査でまとめて取る (3 回呼ぶより速い)
    DynamicLazySegTree::Agg q = seg.range_query(0, 10);
    std::cout << "range_query(0,10) = {sum=" << q.sum << ", min=" << q.mn << ", max=" << q.mx << "}\n";

    // 【巨大区間への更新】未生成の座標にもちゃんと効く (ノードはほとんど作られない)
    seg.range_assign(-1000000000000000000LL, 1000000000000000000LL, 4);
    std::cout << "巨大 range_assign 後 point_get(-777777) = " << seg.point_get(-777777) << "\n";
    std::cout << "                     range_sum(0, 100)  = " << seg.range_sum(0, 100) << "\n";

    // 【LLONG_MAX の扱い】区間は半開なので LLONG_MAX は含まれない。点操作なら OK。
    seg.point_set(LMAX, 999);
    std::cout << "point_get(LLONG_MAX)        = " << seg.point_get(LMAX) << "\n";
    seg.point_set(LMIN, -999);
    std::cout << "point_get(LLONG_MIN)        = " << seg.point_get(LMIN) << "\n";
    std::cout << "range_min(LLONG_MIN, 0)     = " << seg.range_min(LMIN, 0) << "\n";

    // 【all_query】全座標の集約 (長さ 2^64 なので sum は溢れる点に注意)
    DynamicLazySegTree::Agg all = seg.all_query();
    std::cout << "all_query()                 = {min=" << all.mn << ", max=" << all.mx << "}\n";

    // 【空区間】l >= r なら sum=0, min=INF, max=-INF
    std::cout << "range_sum(7, 7)             = " << seg.range_sum(7, 7) << "   (空区間)\n";

    // 【node_count】取得ではノードが増えないことを確認
    size_t before = seg.node_count();
    ll acc = 0;
    for (int i = 0; i < 10000; i++) acc ^= seg.range_sum(-i, i + 1) ^ seg.point_get(i);
    (void)acc;
    std::cout << "node_count()                = " << before << " → " << seg.node_count()
              << "  (取得 1 万回してもノードは増えない)\n";
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    tester::random_test();
    tester::bench_nodes();
    demo();
    return 0;
}
