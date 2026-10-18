// =============================================================================
//  ③ Segment Tree Beats   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  値の型      : long long 固定
//  区間の流儀  : 半開区間 [l, r)  / 0-indexed
//  依存        : 標準ライブラリのみ (C++17)
//
//  更新 : point_set / range_add / range_chmin / range_chmax
//  取得 : point_get / range_sum / range_min / range_max
//
//  仕組み: 各ノードに「最大値・狭義2番目の最大値・最大値の個数」
//          「最小値・狭義2番目の最小値・最小値の個数」「総和」「加算の遅延」を持つ。
//          chmin(x) は「max_v <= x なら何もしない」「2nd_max < x < max_v なら
//          最大値だけを一括で x に落とす」で再帰を打ち切る。これで
//          ならし計算量 O(n log^2 n) になる。
//
//  計算量: 構築 O(n) / 各クエリ ならし O(log^2 n)
//  オーバーフロー注意: sum は long long。|値| * n が 9.2e18 を超えないこと。
// =============================================================================

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

using ll = long long;

struct SegTreeBeats {
    static constexpr ll INF = std::numeric_limits<ll>::max();

    SegTreeBeats() : SegTreeBeats(0) {}
    explicit SegTreeBeats(int n) : SegTreeBeats(std::vector<ll>(n, 0)) {}
    explicit SegTreeBeats(const std::vector<ll>& a) : _n((int)a.size()) {
        n0 = 1;
        while (n0 < _n) n0 <<= 1;
        size_t m = 2 * (size_t)n0;
        sum_.assign(m, 0);
        max_v.assign(m, 0);
        smax_v.assign(m, 0);
        min_v.assign(m, 0);
        smin_v.assign(m, 0);
        max_c.assign(m, 0);
        min_c.assign(m, 0);
        len_.assign(m, 0);
        ladd.assign(m, 0);
        for (int i = 0; i < n0; i++) {
            int k = n0 - 1 + i;
            if (i < _n) {
                sum_[k] = max_v[k] = min_v[k] = a[i];
                smax_v[k] = -INF;
                smin_v[k] = INF;
                max_c[k] = min_c[k] = 1;
                len_[k] = 1;
            } else {  // 使わない葉 (長さ 0 の番兵)
                sum_[k] = 0;
                max_v[k] = smax_v[k] = -INF;
                min_v[k] = smin_v[k] = INF;
                max_c[k] = min_c[k] = 0;
                len_[k] = 0;
            }
        }
        for (int k = n0 - 2; k >= 0; k--) {
            len_[k] = len_[2 * k + 1] + len_[2 * k + 2];
            pull(k);
        }
    }

    int size() const { return _n; }

    // ---- 更新 ----
    void point_set(int p, ll x) {  // a[p] = x
        assert(0 <= p && p < _n);
        _point_set(p, x, 0, 0, n0);
    }
    void range_add(int l, int r, ll x) {  // a[i] += x
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
        _range_add(x, l, r, 0, 0, n0);
    }
    void range_chmin(int l, int r, ll x) {  // a[i] = min(a[i], x)
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
        _range_chmin(x, l, r, 0, 0, n0);
    }
    void range_chmax(int l, int r, ll x) {  // a[i] = max(a[i], x)
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
        _range_chmax(x, l, r, 0, 0, n0);
    }

    // ---- 取得 ----
    ll point_get(int p) {  // a[p]
        assert(0 <= p && p < _n);
        return _query_max(p, p + 1, 0, 0, n0);
    }
    ll range_sum(int l, int r) {  // 空区間は 0
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return 0;
        return _query_sum(l, r, 0, 0, n0);
    }
    ll range_min(int l, int r) {  // 空区間は  INF
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return INF;
        return _query_min(l, r, 0, 0, n0);
    }
    ll range_max(int l, int r) {  // 空区間は -INF
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return -INF;
        return _query_max(l, r, 0, 0, n0);
    }

   private:
    int _n = 0, n0 = 1;
    std::vector<ll> sum_, max_v, smax_v, min_v, smin_v, max_c, min_c, len_, ladd;

    // 最大値だけを x に落とす (前提: smax_v[k] < x < max_v[k])
    void update_node_max(int k, ll x) {
        sum_[k] += (x - max_v[k]) * max_c[k];
        if (max_v[k] == min_v[k]) {
            max_v[k] = min_v[k] = x;
        } else if (max_v[k] == smin_v[k]) {
            max_v[k] = smin_v[k] = x;
        } else {
            max_v[k] = x;
        }
    }
    // 最小値だけを x に上げる (前提: min_v[k] < x < smin_v[k])
    void update_node_min(int k, ll x) {
        sum_[k] += (x - min_v[k]) * min_c[k];
        if (max_v[k] == min_v[k]) {
            max_v[k] = min_v[k] = x;
        } else if (smax_v[k] == min_v[k]) {
            min_v[k] = smax_v[k] = x;
        } else {
            min_v[k] = x;
        }
    }
    void addall(int k, ll x) {
        max_v[k] += x;
        if (smax_v[k] != -INF) smax_v[k] += x;
        min_v[k] += x;
        if (smin_v[k] != INF) smin_v[k] += x;
        sum_[k] += x * len_[k];
        ladd[k] += x;
    }
    void push(int k) {
        if (k >= n0 - 1) return;  // 葉
        if (ladd[k] != 0) {
            addall(2 * k + 1, ladd[k]);
            addall(2 * k + 2, ladd[k]);
            ladd[k] = 0;
        }
        if (max_v[k] < max_v[2 * k + 1]) update_node_max(2 * k + 1, max_v[k]);
        if (min_v[2 * k + 1] < min_v[k]) update_node_min(2 * k + 1, min_v[k]);
        if (max_v[k] < max_v[2 * k + 2]) update_node_max(2 * k + 2, max_v[k]);
        if (min_v[2 * k + 2] < min_v[k]) update_node_min(2 * k + 2, min_v[k]);
    }
    void pull(int k) {
        int a = 2 * k + 1, b = 2 * k + 2;
        sum_[k] = sum_[a] + sum_[b];
        if (max_v[a] < max_v[b]) {
            max_v[k] = max_v[b];
            max_c[k] = max_c[b];
            smax_v[k] = std::max(max_v[a], smax_v[b]);
        } else if (max_v[a] > max_v[b]) {
            max_v[k] = max_v[a];
            max_c[k] = max_c[a];
            smax_v[k] = std::max(smax_v[a], max_v[b]);
        } else {
            max_v[k] = max_v[a];
            max_c[k] = max_c[a] + max_c[b];
            smax_v[k] = std::max(smax_v[a], smax_v[b]);
        }
        if (min_v[a] < min_v[b]) {
            min_v[k] = min_v[a];
            min_c[k] = min_c[a];
            smin_v[k] = std::min(smin_v[a], min_v[b]);
        } else if (min_v[a] > min_v[b]) {
            min_v[k] = min_v[b];
            min_c[k] = min_c[b];
            smin_v[k] = std::min(min_v[a], smin_v[b]);
        } else {
            min_v[k] = min_v[a];
            min_c[k] = min_c[a] + min_c[b];
            smin_v[k] = std::min(smin_v[a], smin_v[b]);
        }
    }

    void _range_chmin(ll x, int a, int b, int k, int l, int r) {
        if (b <= l || r <= a || max_v[k] <= x) return;         // 触る必要なし
        if (a <= l && r <= b && smax_v[k] < x) {               // 最大値だけ潰せば済む
            update_node_max(k, x);
            return;
        }
        push(k);
        int m = (l + r) >> 1;
        _range_chmin(x, a, b, 2 * k + 1, l, m);
        _range_chmin(x, a, b, 2 * k + 2, m, r);
        pull(k);
    }
    void _range_chmax(ll x, int a, int b, int k, int l, int r) {
        if (b <= l || r <= a || x <= min_v[k]) return;
        if (a <= l && r <= b && x < smin_v[k]) {
            update_node_min(k, x);
            return;
        }
        push(k);
        int m = (l + r) >> 1;
        _range_chmax(x, a, b, 2 * k + 1, l, m);
        _range_chmax(x, a, b, 2 * k + 2, m, r);
        pull(k);
    }
    void _range_add(ll x, int a, int b, int k, int l, int r) {
        if (b <= l || r <= a) return;
        if (a <= l && r <= b) {
            addall(k, x);
            return;
        }
        push(k);
        int m = (l + r) >> 1;
        _range_add(x, a, b, 2 * k + 1, l, m);
        _range_add(x, a, b, 2 * k + 2, m, r);
        pull(k);
    }
    void _point_set(int p, ll x, int k, int l, int r) {
        if (r - l == 1) {
            sum_[k] = max_v[k] = min_v[k] = x;
            smax_v[k] = -INF;
            smin_v[k] = INF;
            max_c[k] = min_c[k] = 1;
            ladd[k] = 0;
            return;
        }
        push(k);
        int m = (l + r) >> 1;
        if (p < m) _point_set(p, x, 2 * k + 1, l, m);
        else _point_set(p, x, 2 * k + 2, m, r);
        pull(k);
    }
    ll _query_sum(int a, int b, int k, int l, int r) {
        if (b <= l || r <= a) return 0;
        if (a <= l && r <= b) return sum_[k];
        push(k);
        int m = (l + r) >> 1;
        return _query_sum(a, b, 2 * k + 1, l, m) + _query_sum(a, b, 2 * k + 2, m, r);
    }
    ll _query_max(int a, int b, int k, int l, int r) {
        if (b <= l || r <= a) return -INF;
        if (a <= l && r <= b) return max_v[k];
        push(k);
        int m = (l + r) >> 1;
        return std::max(_query_max(a, b, 2 * k + 1, l, m), _query_max(a, b, 2 * k + 2, m, r));
    }
    ll _query_min(int a, int b, int k, int l, int r) {
        if (b <= l || r <= a) return INF;
        if (a <= l && r <= b) return min_v[k];
        push(k);
        int m = (l + r) >> 1;
        return std::min(_query_min(a, b, 2 * k + 1, l, m), _query_min(a, b, 2 * k + 2, m, r));
    }
};

// =============================================================================
//  ここから下は「ランダムテスト」と「使用例デモ」。提出時は丸ごと消してよい。
// =============================================================================

namespace tester {

static ll rnd_range(std::mt19937_64& rng, ll lo, ll hi) {
    return lo + (ll)(rng() % (unsigned long long)(hi - lo + 1));
}

static void random_test() {
    std::mt19937_64 rng(20260916ULL);
    const ll VAL = 500;

    for (int rep = 0; rep < 400; rep++) {
        int n = (int)rnd_range(rng, 1, 40);
        std::vector<ll> a(n);
        for (auto& x : a) x = rnd_range(rng, -VAL, VAL);
        SegTreeBeats seg(a);

        for (int t = 0; t < 500; t++) {
            int l = (int)rnd_range(rng, 0, n), r = (int)rnd_range(rng, 0, n);
            if (l > r) std::swap(l, r);
            int p = (int)rnd_range(rng, 0, n - 1);
            int type = (int)rnd_range(rng, 0, 7);

            if (type == 0) {  // point_set
                ll x = rnd_range(rng, -VAL, VAL);
                a[p] = x;
                seg.point_set(p, x);
            } else if (type == 1) {  // range_add
                ll x = rnd_range(rng, -VAL / 2, VAL / 2);
                for (int i = l; i < r; i++) a[i] += x;
                seg.range_add(l, r, x);
            } else if (type == 2) {  // range_chmin
                ll x = rnd_range(rng, -VAL, VAL);
                for (int i = l; i < r; i++) a[i] = std::min(a[i], x);
                seg.range_chmin(l, r, x);
            } else if (type == 3) {  // range_chmax
                ll x = rnd_range(rng, -VAL, VAL);
                for (int i = l; i < r; i++) a[i] = std::max(a[i], x);
                seg.range_chmax(l, r, x);
            } else if (type == 4) {  // point_get
                assert(seg.point_get(p) == a[p]);
            } else if (type == 5) {  // range_sum
                ll s = 0;
                for (int i = l; i < r; i++) s += a[i];
                assert(seg.range_sum(l, r) == s);
            } else if (type == 6) {  // range_min
                ll m = SegTreeBeats::INF;
                for (int i = l; i < r; i++) m = std::min(m, a[i]);
                assert(seg.range_min(l, r) == m);
            } else {  // range_max
                ll m = -SegTreeBeats::INF;
                for (int i = l; i < r; i++) m = std::max(m, a[i]);
                assert(seg.range_max(l, r) == m);
            }
        }
        // 最後に全区間を突き合わせ
        for (int l = 0; l <= n; l++) {
            for (int r = l; r <= n; r++) {
                ll s = 0, mn = SegTreeBeats::INF, mx = -SegTreeBeats::INF;
                for (int i = l; i < r; i++) {
                    s += a[i];
                    mn = std::min(mn, a[i]);
                    mx = std::max(mx, a[i]);
                }
                assert(seg.range_sum(l, r) == s);
                assert(seg.range_min(l, r) == mn);
                assert(seg.range_max(l, r) == mx);
            }
        }
    }

    // 値の重複が非常に多いケース (beats の「2 番目の最大値 == 更新値」の分岐を踏みやすい)
    // + 2 のべき乗境界。最大値の個数が壊れる系のバグはここで落ちる。
    for (int n : {1, 2, 3, 4, 7, 8, 9, 16, 31, 32, 33, 64}) {
        std::vector<ll> a(n);
        for (auto& x : a) x = rnd_range(rng, 0, 2);
        SegTreeBeats seg(a);
        for (int t = 0; t < 4000; t++) {
            int l = (int)rnd_range(rng, 0, n), r = (int)rnd_range(rng, 0, n);
            if (l > r) std::swap(l, r);
            int type = (int)rnd_range(rng, 0, 3);
            ll x = rnd_range(rng, -2, 4);
            if (type == 0) {
                for (int i = l; i < r; i++) a[i] = std::min(a[i], x);
                seg.range_chmin(l, r, x);
            } else if (type == 1) {
                for (int i = l; i < r; i++) a[i] = std::max(a[i], x);
                seg.range_chmax(l, r, x);
            } else if (type == 2) {
                ll v = rnd_range(rng, -1, 1);
                for (int i = l; i < r; i++) a[i] += v;
                seg.range_add(l, r, v);
            } else {
                int p = (int)rnd_range(rng, 0, n - 1);
                a[p] = x;
                seg.point_set(p, x);
            }
            for (int rep2 = 0; rep2 < 3; rep2++) {
                int l2 = (int)rnd_range(rng, 0, n), r2 = (int)rnd_range(rng, 0, n);
                if (l2 > r2) std::swap(l2, r2);
                ll s = 0, mn = SegTreeBeats::INF, mx = -SegTreeBeats::INF;
                for (int i = l2; i < r2; i++) {
                    s += a[i];
                    mn = std::min(mn, a[i]);
                    mx = std::max(mx, a[i]);
                }
                assert(seg.range_sum(l2, r2) == s);
                assert(seg.range_min(l2, r2) == mn);
                assert(seg.range_max(l2, r2) == mx);
            }
            for (int i = 0; i < n; i++) assert(seg.point_get(i) == a[i]);
        }
        // 最後に全区間を突き合わせ
        for (int l = 0; l <= n; l++) {
            for (int r = l; r <= n; r++) {
                ll s = 0, mn = SegTreeBeats::INF, mx = -SegTreeBeats::INF;
                for (int i = l; i < r; i++) {
                    s += a[i];
                    mn = std::min(mn, a[i]);
                    mx = std::max(mx, a[i]);
                }
                assert(seg.range_sum(l, r) == s);
                assert(seg.range_min(l, r) == mn);
                assert(seg.range_max(l, r) == mx);
            }
        }
    }

    std::cout << "[segtree_beats] random test: OK" << std::endl;
}

}  // namespace tester

// =============================================================================
//  使用例デモ : 全機能を 1 箇所で使う。実戦ではここを見れば思い出せる。
// =============================================================================
static void dump(const char* label, SegTreeBeats& seg) {
    std::cout << label << " : ";
    for (int i = 0; i < seg.size(); i++) std::cout << seg.point_get(i) << (i + 1 == seg.size() ? '\n' : ' ');
}

static void demo() {
    std::cout << "\n===== ③ Segment Tree Beats 使用例デモ =====\n";

    // 【構築1】初期値の配列から  O(n)
    std::vector<ll> a = {5, 3, 8, 1, 9, 2, 7, 4};
    SegTreeBeats seg(a);
    dump("初期状態                ", seg);

    // 【構築2】要素数だけ指定 (全要素 0)
    SegTreeBeats seg_zero(8);
    std::cout << "SegTreeBeats(8) の全体和 = " << seg_zero.range_sum(0, 8) << "\n";

    // 【range_chmin】[0, 8) の各要素を min(a[i], 6) にする
    seg.range_chmin(0, 8, 6);
    dump("range_chmin(0,8,6)      ", seg);

    // 【range_chmax】[2, 6) の各要素を max(a[i], 4) にする
    seg.range_chmax(2, 6, 4);
    dump("range_chmax(2,6,4)      ", seg);

    // 【range_add】[1, 5) に +10
    seg.range_add(1, 5, 10);
    dump("range_add(1,5,+10)      ", seg);

    // 【point_set】1 点代入
    seg.point_set(7, 100);
    dump("point_set(7,100)        ", seg);

    // 【point_get】1 点取得
    std::cout << "point_get(7)  = " << seg.point_get(7) << "\n";

    // 【range_sum / range_min / range_max】区間取得
    std::cout << "range_sum(0,8) = " << seg.range_sum(0, 8) << "\n";
    std::cout << "range_min(0,8) = " << seg.range_min(0, 8) << "\n";
    std::cout << "range_max(0,8) = " << seg.range_max(0, 8) << "\n";
    std::cout << "range_sum(2,5) = " << seg.range_sum(2, 5) << "\n";
    std::cout << "range_min(2,5) = " << seg.range_min(2, 5) << "\n";
    std::cout << "range_max(2,5) = " << seg.range_max(2, 5) << "\n";

    // 【chmin + chmax の合わせ技】区間を [3, 5] にクランプする
    seg.range_chmax(0, 8, 3);
    seg.range_chmin(0, 8, 5);
    dump("[3,5] にクランプ        ", seg);

    // 【空区間】l == r のときは sum=0, min=INF, max=-INF
    std::cout << "range_sum(4,4) = " << seg.range_sum(4, 4) << "  (空区間)\n";
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    tester::random_test();
    demo();
    return 0;
}
