// =============================================================================
//  【ごり押し】最小・最大の位置のセグ木（区間の最小・最大と、それをとる最も左・最も右の位置・個数を 1 回で。1 点更新版・区間更新版）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ ごり押しテンプレート（考察を飛ばして、直さずにそのまま呼べる形）
//    省ける考察     : 位置を取るために「最大値を聞いてから二分探索」の 2 手に分ける・(値, 位置) の pair を自分で載せる
//    向いている問題 : 「1 か所だけ変えた N 通りを全部試す（差分だけ更新 → 全体を 1 回聞く → 戻す）」、得票の 1 位（同点は番号が小さい方）、
//                     累積和の列に区間加算して max − min を見る、など（類題：ABC329 D、ABC483 C・D）
//    切り替え       : ArgMinMaxSegTree<MAX_ONLY> / <MIN_ONLY> / <BOTH>（省略すると BOTH）。問題を読んで要る側だけにすると速い
//                     （中に持つものも片側だけになる）。返す struct はどれでも同じ 8 項目で、使わない側は空の値が入る
//    今あるもので書くと : lazy_segtree_bs.cpp の all_prod().mx → find_first_geq(0, mx) の 2 手、RangeMaxSegTreeBinary.txt の first_ge
//    2 番目に大きい・小さい値も要るとき : ds_segtree_argminmax2.cpp（ArgMinMax2SegTree。この 8 項目 + 2 番目の 8 項目）
//
//  持つ情報（区間 [l, r) について。個数 = その値をとる位置の数）
//    min / first_min / last_min / min_count   … 最小値・それをとる最も左の位置・最も右の位置・個数
//    max / first_max / last_max / max_count   … 最大値・それをとる最も左の位置・最も右の位置・個数
//    例 [3, 1, 4, 1, 5] の全体：min 1・first_min 1・last_min 3・min_count 2 / max 5・first_max 4・last_max 4・max_count 1
//    空区間と、使わない側（MAX_ONLY のときの最小側など）：min = INF（LLONG_MAX）・max = −INF・位置 −1・個数 0
//
//  機能一覧（値は long long。位置は 0〜N−1。区間は半開区間 [l, r)）                                        計算量
//    ArgMinMaxSegTree<MODE> seg(a)          1 点更新版。a は vector<long long> でも vector<int> でもよい        O(N)
//    ArgMinMaxSegTree<MODE> seg(n)          長さ n・全部 0                                                     O(N)
//    ArgMinMaxLazySegTree<MODE> seg(a) / seg(n)   区間更新版（下の「区間更新版だけ」の関数も使える）           O(N)
//    共通
//      seg.set(i, x) / seg.add(i, x)        a[i] = x / a[i] += x                                               O(log N)
//      seg.get(i)                           a[i]                                       1 点更新版 O(1)・区間更新版 O(log N)
//      seg.prod(l, r)                       8 項目をまとめて struct で返す（auto s = seg.prod(l, r); s.first_max）   O(log N)
//      seg.all_prod()                       全体の 8 項目                                                      O(1)
//      seg.min(l, r) / seg.max(l, r) / seg.first_min(l, r) / seg.last_min(l, r) / seg.min_count(l, r)
//        / seg.first_max(l, r) / seg.last_max(l, r) / seg.max_count(l, r)     区間の 1 項目だけ                 O(log N)
//      seg.min() / seg.max() / … / seg.max_count()    全体の 1 項目だけ                                        O(1)
//      seg.size() / seg.dump()（a を標準エラーに表示）
//    区間更新版だけ
//      seg.range_add(l, r, x)               [l, r) の各値に x を足す                                           O(log N)
//      seg.assign(l, r, x)                  [l, r) の各値を x にする                                           O(log N)
//      seg.range_affine(l, r, a, b)         [l, r) の各値を a·x + b にする（a < 0 は BOTH のときだけ。
//                                           符号が反転すると最小と最大が入れ替わる。MAX_ONLY / MIN_ONLY で a < 0 は assert で止まる）  O(log N)
//    使わない側の 1 項目の関数（MAX_ONLY での seg.min() など）を呼ぶとコンパイルエラー
//
//  速さの違い（どれを呼ぶか迷ったら）
//    1 項目の関数（区間）: BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。
//                          同じ区間で 2 つ以上の値が要るなら prod を 1 回
//    1 項目の関数（全体）: どの切り替えでも all_prod と同じ速さ（O(1)）
//    区間更新版では 1 項目の関数と prod の差は小さい（遅延の伝搬はどちらも同じだけする）
//    1 点更新しかしないなら 1 点更新版の方が速い（区間更新版の約 1.2〜1.3 倍）
//
//  速さの目安（N = Q = 2×10^5、更新と区間の質問が半分ずつ、prod で 1 項目を読む。計測した環境は AtCoder より遅め）
//    1 点更新版（1 点更新 + 区間）: MAX_ONLY 約 0.06 秒・BOTH 約 0.10 秒（最大値だけのセグ木は約 0.025 秒）
//    区間更新版（区間加算 + 区間）: MAX_ONLY 約 0.15 秒・BOTH 約 0.19 秒（lazy_segtree_bs.cpp は約 0.18 秒）
//    10^6 だと 1 点更新版 BOTH 約 0.65 秒、区間更新版 BOTH 約 1.5 秒（lazy_segtree_bs.cpp は約 1.3 秒）なので、10^6 級で両方要るときは時間に注意
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（ds_segtree_argminmax2.cpp と両方貼っても衝突しない）
//    amm_test::random_test()  … 3 つの切り替え × 1 点更新版・区間更新版を、毎回数え直す愚直解と突き合わせる
//    amm_test::benchmark()    … 速度計測（切り替え・1 項目の関数・最大値だけのセグ木との比較）
//    amm_demo::demo()         … 全機能の使い方デモ（日本語の説明付き。ABC483 C・D の例つき）
// =============================================================================
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef ARGMINMAX_MODE_DEFINED
#define ARGMINMAX_MODE_DEFINED  // ds_segtree_argminmax2.cpp と両方貼っても二重定義にならないための目印
// 最大だけ・最小だけ・両方（テンプレート引数で指定する。省略すると BOTH）
enum ArgMinMaxMode : int { MAX_ONLY = 1, MIN_ONLY = 2, BOTH = 3 };
#endif

namespace amm_internal {
using ll = long long;
constexpr ll INF = std::numeric_limits<ll>::max();
// 片側の 4 項目（値・最も左・最も右・個数。c = 0 なら空）
struct Side {
    ll v;
    int f, l, c;
};
template <bool MX>
constexpr Side empty_side() { return Side{MX ? -INF : INF, -1, -1, 0}; }
template <bool MX>
inline Side merge(const Side& L, const Side& R) {  // L が左、R が右の区間
    if (MX ? L.v > R.v : L.v < R.v) return L;
    if (MX ? R.v > L.v : R.v < L.v) return R;
    return Side{L.v, L.c ? L.f : R.f, R.c ? R.l : L.l, L.c + R.c};
}
inline void affine_side(Side& s, ll a, ll b) {
    if (s.c) s.v = a * s.v + b;
}
// 返す struct（どの切り替えでも同じ 8 項目。使わない側は空の値）
struct Info {
    ll min, max;
    int first_min, last_min, min_count, first_max, last_max, max_count;
    friend std::ostream& operator<<(std::ostream& os, const Info& s) {
        return os << "{min " << s.min << " first_min " << s.first_min << " last_min " << s.last_min << " min_count " << s.min_count
                  << " / max " << s.max << " first_max " << s.first_max << " last_max " << s.last_max << " max_count " << s.max_count << "}";
    }
};
// 中に持つノード（切り替えで形が変わる：使わない側は持たない）
template <int MODE>
struct Core {
    struct NMax { Side mx; };
    struct NMin { Side mn; };
    struct NBoth { Side mn, mx; };
    using Node = std::conditional_t<MODE == MAX_ONLY, NMax, std::conditional_t<MODE == MIN_ONLY, NMin, NBoth>>;
    static constexpr bool HAS_MIN = (MODE & MIN_ONLY) != 0;
    static constexpr bool HAS_MAX = (MODE & MAX_ONLY) != 0;
    static Node e() {
        Node n{};
        if constexpr (HAS_MIN) n.mn = empty_side<false>();
        if constexpr (HAS_MAX) n.mx = empty_side<true>();
        return n;
    }
    static Node leaf(int i, ll v) {
        Node n{};
        if constexpr (HAS_MIN) n.mn = Side{v, i, i, 1};
        if constexpr (HAS_MAX) n.mx = Side{v, i, i, 1};
        return n;
    }
    static Node op(const Node& x, const Node& y) {
        Node n{};
        if constexpr (HAS_MIN) n.mn = merge<false>(x.mn, y.mn);
        if constexpr (HAS_MAX) n.mx = merge<true>(x.mx, y.mx);
        return n;
    }
    template <bool MX>
    static const Side& side(const Node& n) {
        if constexpr (MX) return n.mx;
        else return n.mn;
    }
    static const Side& any_side(const Node& n) {
        if constexpr (HAS_MAX) return n.mx;
        else return n.mn;
    }
    static Info to_info(const Node& n) {
        Info s{INF, -INF, -1, -1, 0, -1, -1, 0};
        if constexpr (HAS_MIN) s.min = n.mn.v, s.first_min = n.mn.f, s.last_min = n.mn.l, s.min_count = n.mn.c;
        if constexpr (HAS_MAX) s.max = n.mx.v, s.first_max = n.mx.f, s.last_max = n.mx.l, s.max_count = n.mx.c;
        return s;
    }
    // 作用 a·x + b（a = 0 は代入）。lo・len はノードが受け持つ区間（代入で使う）
    static void apply(Node& n, ll a, ll b, int lo, int len) {
        if (any_side(n).c == 0) return;
        if (a == 0) {
            if constexpr (HAS_MIN) n.mn = Side{b, lo, lo + len - 1, len};
            if constexpr (HAS_MAX) n.mx = Side{b, lo, lo + len - 1, len};
        } else if (a > 0) {
            if constexpr (HAS_MIN) affine_side(n.mn, a, b);
            if constexpr (HAS_MAX) affine_side(n.mx, a, b);
        } else {
            if constexpr (MODE == BOTH) {
                std::swap(n.mn, n.mx);  // 符号が反転すると最小と最大が入れ替わる
                affine_side(n.mn, a, b), affine_side(n.mx, a, b);
            }
        }
    }
};
}  // namespace amm_internal

// 1 点更新版
template <int MODE = BOTH>
struct ArgMinMaxSegTree {
    static_assert(MODE == MAX_ONLY || MODE == MIN_ONLY || MODE == BOTH, "MODE は MAX_ONLY / MIN_ONLY / BOTH");
    using S = amm_internal::Info;
    static constexpr long long INF = amm_internal::INF;

    ArgMinMaxSegTree() : ArgMinMaxSegTree(0) {}
    explicit ArgMinMaxSegTree(int n) : ArgMinMaxSegTree(std::vector<long long>(n, 0)) {}
    template <class V>
    explicit ArgMinMaxSegTree(const std::vector<V>& a) : n_((int)a.size()) {
        sz_ = 1;
        while (sz_ < n_) sz_ <<= 1;
        d_.assign(2 * (size_t)sz_, C::e());
        for (int i = 0; i < n_; ++i) d_[sz_ + i] = C::leaf(i, (long long)a[i]);
        for (int k = sz_ - 1; k >= 1; --k) d_[k] = C::op(d_[2 * k], d_[2 * k + 1]);
    }
    int size() const { return n_; }

    // a[i] = x
    void set(int i, long long x) {
        assert(0 <= i && i < n_);
        int k = i + sz_;
        d_[k] = C::leaf(i, x);
        for (k >>= 1; k; k >>= 1) d_[k] = C::op(d_[2 * k], d_[2 * k + 1]);
    }
    // a[i] += x
    void add(int i, long long x) { set(i, get(i) + x); }
    // a[i]
    long long get(int i) const {
        assert(0 <= i && i < n_);
        return C::any_side(d_[i + sz_]).v;
    }

    // [l, r) の 8 項目をまとめて返す。BOTH で片側の値だけ要るときは 1 項目の関数の方が速い。
    // MAX_ONLY / MIN_ONLY では 1 項目の関数と同じ速さ。2 つ以上の値が要るならこちらを 1 回
    S prod(int l, int r) const {
        assert(0 <= l && l <= r && r <= n_);
        Node sl = C::e(), sr = C::e();
        for (l += sz_, r += sz_; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sl = C::op(sl, d_[l++]);
            if (r & 1) sr = C::op(d_[--r], sr);
        }
        return C::to_info(C::op(sl, sr));
    }
    // 全体の 8 項目。1 項目の関数の全体版と同じ速さ（O(1)）
    S all_prod() const { return C::to_info(d_[1]); }

    // ---- 1 項目の関数（(l, r) で区間 [l, r)・引数なしで全体）。区間版は BOTH のときだけ prod より速い。全体版は all_prod と同じ速さ
    // [l, r) の最小値（空なら INF）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long min(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).v; }
    // 全体の最小値（空なら INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.v; }
    // [l, r) の最小値をとる最も左の位置（空なら −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_min(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).f; }
    // 全体の最小値をとる最も左の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.f; }
    // [l, r) の最小値をとる最も右の位置（空なら −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_min(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).l; }
    // 全体の最小値をとる最も右の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.l; }
    // [l, r) の最小値の個数（空なら 0）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int min_count(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).c; }
    // 全体の最小値の個数（空なら 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int min_count() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.c; }
    // [l, r) の最大値（空なら −INF）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long max(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).v; }
    // 全体の最大値（空なら −INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.v; }
    // [l, r) の最大値をとる最も左の位置（空なら −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_max(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).f; }
    // 全体の最大値をとる最も左の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.f; }
    // [l, r) の最大値をとる最も右の位置（空なら −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_max(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).l; }
    // 全体の最大値をとる最も右の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.l; }
    // [l, r) の最大値の個数（空なら 0）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int max_count(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).c; }
    // 全体の最大値の個数（空なら 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int max_count() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.c; }

    // a を標準エラーに表示（デバッグ用）
    void dump() const {
        std::cerr << "a = [";
        for (int i = 0; i < n_; ++i) std::cerr << (i ? ", " : "") << get(i);
        std::cerr << "]\n";
    }

   private:
    using C = amm_internal::Core<MODE>;
    using Node = typename C::Node;
    using Side = amm_internal::Side;
    int n_ = 0, sz_ = 1;
    std::vector<Node> d_;
    template <bool MX>
    Side side_prod(int l, int r) const {
        assert(0 <= l && l <= r && r <= n_);
        Side sl = amm_internal::empty_side<MX>(), sr = sl;
        for (l += sz_, r += sz_; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sl = amm_internal::merge<MX>(sl, C::template side<MX>(d_[l++]));
            if (r & 1) sr = amm_internal::merge<MX>(C::template side<MX>(d_[--r]), sr);
        }
        return amm_internal::merge<MX>(sl, sr);
    }
};

// 区間更新版（作用は a·x + b。加算 a = 1・代入 a = 0・符号反転 a = −1 はその特別な場合）
template <int MODE = BOTH>
struct ArgMinMaxLazySegTree {
    static_assert(MODE == MAX_ONLY || MODE == MIN_ONLY || MODE == BOTH, "MODE は MAX_ONLY / MIN_ONLY / BOTH");
    using S = amm_internal::Info;
    static constexpr long long INF = amm_internal::INF;

    ArgMinMaxLazySegTree() : ArgMinMaxLazySegTree(0) {}
    explicit ArgMinMaxLazySegTree(int n) : ArgMinMaxLazySegTree(std::vector<long long>(n, 0)) {}
    template <class V>
    explicit ArgMinMaxLazySegTree(const std::vector<V>& a) : n_((int)a.size()) {
        lg_ = 0;
        while ((1 << lg_) < n_) ++lg_;
        sz_ = 1 << lg_;
        d_.assign(2 * (size_t)sz_, C::e());
        lz_.assign((size_t)sz_, F{1, 0});
        for (int i = 0; i < n_; ++i) d_[sz_ + i] = C::leaf(i, (long long)a[i]);
        for (int k = sz_ - 1; k >= 1; --k) update(k);
    }
    int size() const { return n_; }

    // a[i] = x
    void set(int i, long long x) {
        assert(0 <= i && i < n_);
        int k = i + sz_;
        for (int h = lg_; h >= 1; --h) push(k >> h);
        d_[k] = C::leaf(i, x);
        for (int h = 1; h <= lg_; ++h) update(k >> h);
    }
    // a[i] += x
    void add(int i, long long x) { set(i, get(i) + x); }
    // a[i]
    long long get(int i) {
        assert(0 <= i && i < n_);
        int k = i + sz_;
        for (int h = lg_; h >= 1; --h) push(k >> h);
        return C::any_side(d_[k]).v;
    }

    // [l, r) の各値に x を足す
    void range_add(int l, int r, long long x) { apply_range(l, r, F{1, x}); }
    // [l, r) の各値を x にする
    void assign(int l, int r, long long x) { apply_range(l, r, F{0, x}); }
    // [l, r) の各値を a·x + b にする（a < 0 は BOTH のときだけ）
    void range_affine(int l, int r, long long a, long long b) {
        assert((a >= 0 || MODE == BOTH) && "負の数を掛けるのは BOTH のときだけ（片側だけでは掛けたあとの最大・最小が分からない）");
        apply_range(l, r, F{a, b});
    }

    // [l, r) の 8 項目をまとめて返す。BOTH で片側の値だけ要るときは 1 項目の関数の方が速い。
    // MAX_ONLY / MIN_ONLY では 1 項目の関数と同じ速さ。2 つ以上の値が要るならこちらを 1 回
    S prod(int l, int r) {
        assert(0 <= l && l <= r && r <= n_);
        if (l == r) return C::to_info(C::e());
        push_range(l + sz_, r + sz_);
        Node sl = C::e(), sr = C::e();
        for (l += sz_, r += sz_; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sl = C::op(sl, d_[l++]);
            if (r & 1) sr = C::op(d_[--r], sr);
        }
        return C::to_info(C::op(sl, sr));
    }
    // 全体の 8 項目。1 項目の関数の全体版と同じ速さ（O(1)）
    S all_prod() const { return C::to_info(d_[1]); }

    // ---- 1 項目の関数（(l, r) で区間 [l, r)・引数なしで全体）。区間版は BOTH のときだけ prod より速い。全体版は all_prod と同じ速さ
    // [l, r) の最小値（空なら INF）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long min(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).v; }
    // 全体の最小値（空なら INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.v; }
    // [l, r) の最小値をとる最も左の位置（空なら −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_min(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).f; }
    // 全体の最小値をとる最も左の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.f; }
    // [l, r) の最小値をとる最も右の位置（空なら −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_min(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).l; }
    // 全体の最小値をとる最も右の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.l; }
    // [l, r) の最小値の個数（空なら 0）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int min_count(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).c; }
    // 全体の最小値の個数（空なら 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int min_count() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.c; }
    // [l, r) の最大値（空なら −INF）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long max(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).v; }
    // 全体の最大値（空なら −INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.v; }
    // [l, r) の最大値をとる最も左の位置（空なら −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_max(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).f; }
    // 全体の最大値をとる最も左の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.f; }
    // [l, r) の最大値をとる最も右の位置（空なら −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_max(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).l; }
    // 全体の最大値をとる最も右の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.l; }
    // [l, r) の最大値の個数（空なら 0）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int max_count(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).c; }
    // 全体の最大値の個数（空なら 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int max_count() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.c; }

    // a を標準エラーに表示（デバッグ用）
    void dump() {
        std::cerr << "a = [";
        for (int i = 0; i < n_; ++i) std::cerr << (i ? ", " : "") << get(i);
        std::cerr << "]\n";
    }

   private:
    using C = amm_internal::Core<MODE>;
    using Node = typename C::Node;
    using Side = amm_internal::Side;
    struct F {
        long long a, b;
    };
    int n_ = 0, lg_ = 0, sz_ = 1;
    std::vector<Node> d_;
    std::vector<F> lz_;
    void update(int k) { d_[k] = C::op(d_[2 * k], d_[2 * k + 1]); }
    void all_apply(int k, const F& f) {
        if (f.a == 0) {  // 代入：ノード k が受け持つ区間 [lo, lo + len)
            const int dep = 31 - __builtin_clz((unsigned)k), len = sz_ >> dep;
            C::apply(d_[k], 0, f.b, (k - (1 << dep)) * len, len);
        } else {
            C::apply(d_[k], f.a, f.b, 0, 0);
        }
        if (k < sz_) lz_[k] = F{f.a * lz_[k].a, f.a * lz_[k].b + f.b};
    }
    void push(int k) {
        if (lz_[k].a == 1 && lz_[k].b == 0) return;
        all_apply(2 * k, lz_[k]), all_apply(2 * k + 1, lz_[k]);
        lz_[k] = F{1, 0};
    }
    void push_range(int l, int r) {  // l, r は葉の番号
        for (int h = lg_; h >= 1; --h) {
            if (((l >> h) << h) != l) push(l >> h);
            if (((r >> h) << h) != r) push((r - 1) >> h);
        }
    }
    void apply_range(int l, int r, F f) {
        assert(0 <= l && l <= r && r <= n_);
        if (l == r) return;
        l += sz_, r += sz_;
        push_range(l, r);
        for (int a = l, b = r; a < b; a >>= 1, b >>= 1) {
            if (a & 1) all_apply(a++, f);
            if (b & 1) all_apply(--b, f);
        }
        for (int h = 1; h <= lg_; ++h) {
            if (((l >> h) << h) != l) update(l >> h);
            if (((r >> h) << h) != r) update((r - 1) >> h);
        }
    }
    template <bool MX>
    Side side_prod(int l, int r) {
        assert(0 <= l && l <= r && r <= n_);
        Side sl = amm_internal::empty_side<MX>(), sr = sl;
        if (l == r) return sl;
        push_range(l + sz_, r + sz_);
        for (l += sz_, r += sz_; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sl = amm_internal::merge<MX>(sl, C::template side<MX>(d_[l++]));
            if (r & 1) sr = amm_internal::merge<MX>(C::template side<MX>(d_[--r]), sr);
        }
        return amm_internal::merge<MX>(sl, sr);
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace amm_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define AMM_CHECK(...)                                                                                  \
    do {                                                                                                \
        if (!(__VA_ARGS__)) {                                                                           \
            std::cerr << "[ds_segtree_argminmax] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                               \
        }                                                                                               \
    } while (0)

using ll = long long;

// [lo, hi] の一様乱数
[[maybe_unused]] static ll rnd(std::mt19937_64& rng, ll lo, ll hi) {
    const unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return (ll)((unsigned long long)lo + (w == 0 ? rng() : rng() % w));
}
[[maybe_unused]] static long long elapsed_ms(std::chrono::steady_clock::time_point a) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
}

// 愚直：[l, r) の 8 項目を数え直す（使わない側は空の値）
static amm_internal::Info brute(const std::vector<ll>& a, int l, int r, int mode) {
    const ll INF = amm_internal::INF;
    amm_internal::Info s{INF, -INF, -1, -1, 0, -1, -1, 0};
    for (int i = l; i < r; ++i) {
        if (mode & MIN_ONLY) {
            if (a[i] < s.min) s.min = a[i], s.first_min = s.last_min = i, s.min_count = 1;
            else if (a[i] == s.min) s.last_min = i, ++s.min_count;
        }
        if (mode & MAX_ONLY) {
            if (a[i] > s.max) s.max = a[i], s.first_max = s.last_max = i, s.max_count = 1;
            else if (a[i] == s.max) s.last_max = i, ++s.max_count;
        }
    }
    return s;
}
static bool same(const amm_internal::Info& x, const amm_internal::Info& y) {
    return x.min == y.min && x.first_min == y.first_min && x.last_min == y.last_min && x.min_count == y.min_count && x.max == y.max &&
           x.first_max == y.first_max && x.last_max == y.last_max && x.max_count == y.max_count;
}

// 1 項目の関数（区間版・全体版）が prod / all_prod と同じか
template <int MODE, class Seg>
static void check_single(Seg& seg, int l, int r, const amm_internal::Info& want) {
    const int n = seg.size();
    const amm_internal::Info all = seg.all_prod();
    if constexpr ((MODE & MIN_ONLY) != 0) {
        AMM_CHECK(seg.min(l, r) == want.min && seg.first_min(l, r) == want.first_min && seg.last_min(l, r) == want.last_min && seg.min_count(l, r) == want.min_count);
        AMM_CHECK(seg.min() == all.min && seg.first_min() == all.first_min && seg.last_min() == all.last_min && seg.min_count() == all.min_count);
    }
    if constexpr ((MODE & MAX_ONLY) != 0) {
        AMM_CHECK(seg.max(l, r) == want.max && seg.first_max(l, r) == want.first_max && seg.last_max(l, r) == want.last_max && seg.max_count(l, r) == want.max_count);
        AMM_CHECK(seg.max() == all.max && seg.first_max() == all.first_max && seg.last_max() == all.last_max && seg.max_count() == all.max_count);
    }
    (void)n;
}

template <int MODE, bool LAZY>
static long long random_test_one(std::mt19937_64& rng) {
    long long cases = 0;
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 14);
        const ll V = rnd(rng, 1, 4);  // 値の範囲を狭くして同じ値をたくさん作る
        std::vector<ll> a(n);
        for (auto& x : a) x = rnd(rng, -V, V);
        std::conditional_t<LAZY, ArgMinMaxLazySegTree<MODE>, ArgMinMaxSegTree<MODE>> seg(a);
        for (int q = 0; q < 40; ++q) {
            const int type = (int)rnd(rng, 0, LAZY ? 7 : 3);
            int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
            if (l > r) std::swap(l, r);
            if (n > 0 && type == 0) {
                const int i = (int)rnd(rng, 0, n - 1);
                const ll x = rnd(rng, -V, V);
                seg.set(i, x), a[i] = x;
            } else if (n > 0 && type == 1) {
                const int i = (int)rnd(rng, 0, n - 1);
                const ll x = rnd(rng, -2, 2);
                seg.add(i, x), a[i] += x;
            } else if constexpr (LAZY) {
                if (type == 4) {
                    const ll x = rnd(rng, -2, 2);
                    seg.range_add(l, r, x);
                    for (int i = l; i < r; ++i) a[i] += x;
                } else if (type == 5) {
                    const ll x = rnd(rng, -V, V);
                    seg.assign(l, r, x);
                    for (int i = l; i < r; ++i) a[i] = x;
                } else if (type == 6) {
                    ll ca = rnd(rng, (MODE == BOTH) ? -1 : 0, 2);
                    const ll cb = rnd(rng, -2, 2);
                    seg.range_affine(l, r, ca, cb);
                    for (int i = l; i < r; ++i) a[i] = ca * a[i] + cb;
                }
            }
            const auto want = brute(a, l, r, MODE);
            AMM_CHECK(same(seg.prod(l, r), want));
            AMM_CHECK(same(seg.all_prod(), brute(a, 0, n, MODE)));
            check_single<MODE>(seg, l, r, want);
            if (n > 0) {
                const int i = (int)rnd(rng, 0, n - 1);
                AMM_CHECK(seg.get(i) == a[i]);
            }
            ++cases;
        }
        AMM_CHECK(seg.size() == n);
    }
    return cases;
}

static void random_test() {
    std::mt19937_64 rng(20261109);
    long long cases = 0;
    cases += random_test_one<MAX_ONLY, false>(rng);
    cases += random_test_one<MIN_ONLY, false>(rng);
    cases += random_test_one<BOTH, false>(rng);
    cases += random_test_one<MAX_ONLY, true>(rng);
    cases += random_test_one<MIN_ONLY, true>(rng);
    cases += random_test_one<BOTH, true>(rng);
    // 端：長さ 0・長さ 1・全部同じ値・大きい値
    {
        ArgMinMaxSegTree<BOTH> e(0);
        AMM_CHECK(e.all_prod().max == -amm_internal::INF && e.all_prod().first_max == -1 && e.all_prod().min_count == 0);
        ArgMinMaxLazySegTree<BOTH> one(std::vector<int>{7});
        one.range_affine(0, 1, -1, 0);
        AMM_CHECK(one.max() == -7 && one.min() == -7 && one.first_max() == 0 && one.max_count() == 1);
        ArgMinMaxLazySegTree<MAX_ONLY> same_v(std::vector<ll>(10, 5));
        same_v.assign(2, 7, 9);
        AMM_CHECK(same_v.max() == 9 && same_v.first_max() == 2 && same_v.last_max() == 6 && same_v.max_count() == 5);
        AMM_CHECK(same_v.all_prod().min == amm_internal::INF && same_v.all_prod().first_min == -1);  // 使わない側は空の値
        ArgMinMaxSegTree<MIN_ONLY> big(std::vector<ll>{(ll)4e18, (ll)-4e18, (ll)-4e18});
        AMM_CHECK(big.min() == (ll)-4e18 && big.first_min() == 1 && big.last_min() == 2);
        cases += 4;
    }
    std::cout << "[ds_segtree_argminmax] random_test: OK (" << cases << " cases)\n";
}

// 比べる相手：最大値だけのセグ木（RangeMaxQuerySegTree.txt と同じ形）
struct PlainMaxSeg {
    int sz;
    std::vector<ll> d;
    explicit PlainMaxSeg(const std::vector<ll>& a) {
        sz = 1;
        while (sz < (int)a.size()) sz <<= 1;
        d.assign(2 * sz, std::numeric_limits<ll>::min());
        for (size_t i = 0; i < a.size(); ++i) d[sz + i] = a[i];
        for (int k = sz - 1; k >= 1; --k) d[k] = std::max(d[2 * k], d[2 * k + 1]);
    }
    void set(int i, ll x) {
        int k = i + sz;
        d[k] = x;
        for (k >>= 1; k; k >>= 1) d[k] = std::max(d[2 * k], d[2 * k + 1]);
    }
    ll prod(int l, int r) const {
        ll s = std::numeric_limits<ll>::min();
        for (l += sz, r += sz; l < r; l >>= 1, r >>= 1) {
            if (l & 1) s = std::max(s, d[l++]);
            if (r & 1) s = std::max(s, d[--r]);
        }
        return s;
    }
};

struct Query {
    int t, l, r;
    ll x;
};
static std::vector<Query> make_queries(std::mt19937_64& rng, int n, int q, bool range_update) {
    std::vector<Query> qs(q);
    for (auto& e : qs) {
        int l = (int)rnd(rng, 0, n - 1), r = (int)rnd(rng, 0, n - 1);
        if (l > r) std::swap(l, r);
        e = Query{(int)rnd(rng, 0, 1), l, r + 1, rnd(rng, range_update ? -3 : 0, range_update ? 3 : 999)};
    }
    return qs;
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_segtree_argminmax] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000, Q = N;
    std::vector<ll> a(N);
    for (auto& x : a) x = rnd(rng, 0, 999);
    const auto pq = make_queries(rng, N, Q, false), rq = make_queries(rng, N, Q, true);
    std::cout << "[ds_segtree_argminmax] benchmark: N = Q = " << N << "（更新と区間の質問が半分ずつ）\n";
    auto line = [](const std::string& name, long long ms, ll chk) {
        std::cout << "  " << name << " " << ms << " ms  (check " << chk % 1000 << ")\n";
    };
    {
        auto t0 = std::chrono::steady_clock::now();
        PlainMaxSeg s(a);
        ll chk = 0;
        for (auto& e : pq) e.t ? s.set(e.l, e.x) : (void)(chk += s.prod(e.l, e.r));
        line("比べる相手：最大値だけのセグ木（1 点更新 + 区間の最大値）       ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMaxSegTree<MAX_ONLY> s(a);
        ll chk = 0;
        for (auto& e : pq) e.t ? s.set(e.l, e.x) : (void)(chk += s.prod(e.l, e.r).first_max);
        line("1 点更新版 MAX_ONLY（1 点更新 + prod）                           ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMaxSegTree<BOTH> s(a);
        ll chk = 0;
        for (auto& e : pq) e.t ? s.set(e.l, e.x) : (void)(chk += s.prod(e.l, e.r).first_max);
        line("1 点更新版 BOTH（1 点更新 + prod）                               ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMaxSegTree<BOTH> s(a);
        ll chk = 0;
        for (auto& e : pq) e.t ? s.set(e.l, e.x) : (void)(chk += s.first_max(e.l, e.r));
        line("1 点更新版 BOTH（1 点更新 + 1 項目の関数 first_max(l, r)）      ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMaxLazySegTree<BOTH> s(a);
        ll chk = 0;
        for (auto& e : pq) e.t ? s.set(e.l, e.x) : (void)(chk += s.prod(e.l, e.r).first_max);
        line("区間更新版 BOTH で 1 点更新だけ（1 点更新 + prod）             ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMaxLazySegTree<MAX_ONLY> s(a);
        ll chk = 0;
        for (auto& e : rq) e.t ? s.range_add(e.l, e.r, e.x) : (void)(chk += s.prod(e.l, e.r).first_max);
        line("区間更新版 MAX_ONLY（区間加算 + prod）                         ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMaxLazySegTree<BOTH> s(a);
        ll chk = 0;
        for (auto& e : rq) e.t ? s.range_add(e.l, e.r, e.x) : (void)(chk += s.prod(e.l, e.r).first_max);
        line("区間更新版 BOTH（区間加算 + prod）                             ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMaxLazySegTree<BOTH> s(a);
        ll chk = 0;
        int k = 0;
        for (auto& e : rq) {
            if (e.t) {
                ++k;
                if (k % 3 == 0) s.range_add(e.l, e.r, e.x);
                else if (k % 3 == 1) s.assign(e.l, e.r, e.x + 500);
                else s.range_affine(e.l, e.r, -1, 0);
            } else {
                chk += s.prod(e.l, e.r).first_max;
            }
        }
        line("区間更新版 BOTH（加算・代入・符号反転 + prod）                 ", elapsed_ms(t0), chk);
    }
}

}  // namespace amm_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace amm_demo {

template <class T>
static std::string to_s(const std::vector<T>& v);
template <class T>
static std::string to_s(const T& v) {  // << で出せるもの（整数・struct など）
    std::ostringstream os;
    os << v;
    return os.str();
}
template <class T>
static std::string to_s(const std::vector<T>& v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "]";
}
[[maybe_unused]] static std::string pad(const std::string& s, size_t w) {
    size_t cols = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = (unsigned char)s[i];
        if ((c & 0xC0) == 0x80) continue;
        cols += (c >= 0xE0 ? 2 : 1);
    }
    return cols >= w ? s : s + std::string(w - cols, ' ');
}
[[maybe_unused]] static void banner(const std::string& s) {
    std::cout << "\n=====================================================================\n"
              << "  " << s << "\n"
              << "=====================================================================\n";
}
[[maybe_unused]] static void api(const std::string& usage, const std::string& ret, const std::string& what, const std::string& cost) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << usage << "\n"
              << "    返り値 : " << ret << "\n"
              << "    内容   : " << what << "\n"
              << "    計算量 : " << cost << "\n";
}
[[maybe_unused]] static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
[[maybe_unused]] static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
[[maybe_unused]] static void stmt(const std::string& s) { code(s + ";"); }
template <class V>
[[maybe_unused]] static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 38) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 66) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) amm_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    amm_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    using ll = long long;
    banner("【ごり押し】最小・最大の位置のセグ木（ds_segtree_argminmax） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 区間 [l, r) の最小・最大と、それをとる最も左・最も右の位置・個数（8 項目）を 1 回で取る\n"
              << "■ 切り替え: <MAX_ONLY> / <MIN_ONLY> / <BOTH>（省略すると BOTH）。要る側だけにすると速い\n"
              << "■ API 一覧\n"
              << "  ArgMinMaxSegTree<MODE> seg(a) / seg(n): set / add / get / prod / all_prod / 1 項目の関数 8 つ（(l, r) で区間・() で全体）/ size / dump\n"
              << "  ArgMinMaxLazySegTree<MODE> seg(a) / seg(n): 上と同じ + range_add / assign / range_affine\n";

    api("【1】 ArgMinMaxSegTree<BOTH> seg(a)  prod(l, r) / all_prod()", "struct（min / first_min / last_min / min_count / max / first_max / last_max / max_count）",
        "区間 [l, r) の 8 項目をまとめて（auto s = seg.prod(l, r); で受けて s.first_max などで読む）", "O(log N) / O(1)");
    RUN(ArgMinMaxSegTree<BOTH> seg(std::vector<ll>{3, 1, 4, 1, 5, 9, 2, 6, 5}));
    SHOW(seg.all_prod(), "全体。最小 1 は位置 1 と 3（2 個）、最大 9 は位置 5");
    RUN(auto s = seg.prod(0, 5));
    SHOW(s.max, "[3, 1, 4, 1, 5] の最大");
    SHOW(s.first_min, "最小 1 の最も左の位置");
    SHOW(s.last_min, "最小 1 の最も右の位置");
    SHOW(s.min_count, "最小 1 の個数");
    SHOW(seg.prod(2, 2), "空区間：min = INF・max = −INF・位置 −1・個数 0");

    api("【2】 1 項目の関数  min / max / first_min / last_min / min_count / first_max / last_max / max_count", "long long（値）/ int（位置・個数）",
        "(l, r) で区間 [l, r)、引数なしで全体。BOTH のときだけ区間版は prod より速い（片側だけ計算する）", "区間 O(log N) / 全体 O(1)");
    SHOW(seg.first_max(), "全体の最大 9 の位置");
    SHOW(seg.max(6, 9), "[2, 6, 5] の最大");
    SHOW(seg.last_max(0, 9), "最大 9 の最も右の位置");
    SHOW(seg.min_count(), "全体の最小 1 の個数");
    note("同じ区間で 2 つ以上の値が要るなら prod を 1 回呼ぶ方が速い（1 項目の関数は呼ぶたびに木をたどる）");

    api("【3】 set(i, x) / add(i, x) / get(i)", "void / void / long long", "a[i] = x / a[i] += x / a[i]", "O(log N) / O(log N) / 1 点更新版 O(1)");
    RUN(seg.set(1, 9));
    SHOW(seg.first_max(), "最大 9 が位置 1 と 5 の 2 か所 → 最も左は 1");
    SHOW(seg.max_count(), "");
    RUN(seg.add(5, -8));
    SHOW(seg.get(5), "9 − 8");
    SHOW(seg.all_prod(), "");

    api("【4】 切り替え  ArgMinMaxSegTree<MAX_ONLY> / <MIN_ONLY>", "（BOTH と同じ struct）",
        "要る側だけ持つので速い。struct の使わない側は空の値。使わない側の 1 項目の関数はコンパイルエラー", "（BOTH より速い）");
    RUN(ArgMinMaxSegTree<MAX_ONLY> mx(std::vector<int>{2, 7, 7, 1}));
    SHOW(mx.first_max(), "最大 7 の最も左");
    SHOW(mx.all_prod(), "最小側は空の値（min = INF・位置 −1・個数 0）");
    note("mx.min() と書くとコンパイルエラー（static_assert：最小側は MIN_ONLY / BOTH で使う）");

    api("【5】 ArgMinMaxLazySegTree<BOTH>  range_add(l, r, x) / assign(l, r, x) / range_affine(l, r, a, b)", "void",
        "[l, r) に x を足す / x にする / a·x + b にする（a < 0 は BOTH だけ。符号が反転すると最小と最大が入れ替わる）", "O(log N)");
    RUN(ArgMinMaxLazySegTree<BOTH> lz(std::vector<ll>{5, 1, 4, 1, 5}));
    SHOW(lz.all_prod(), "");
    RUN(lz.range_add(1, 3, 10));
    SHOW(lz.all_prod(), "[5, 11, 14, 1, 5]");
    RUN(lz.assign(0, 2, 14));
    SHOW(lz.all_prod(), "[14, 14, 14, 1, 5] → 最大 14 が 3 個（位置 0〜2）");
    RUN(lz.range_affine(0, 5, -1, 0));
    SHOW(lz.all_prod(), "符号を反転：[-14, -14, -14, -1, -5] → 最小と最大が入れ替わる");
    SHOW(lz.get(3), "");
    note("MAX_ONLY / MIN_ONLY で負の数を掛けると assert で止まる（片側だけでは掛けたあとの最大・最小が分からない）");

    api("【6】 例：ABC483 C（候補 k が辞退したときの当選者を k = 1〜N で全部）", "（使い方の例）",
        "第一希望が k の票だけ第二希望へ移す → 全体の最大の最も左 → 戻す。各票が動くのは全体で 1 回ずつ", "全体 O((N + M) log N)");
    {
        const int n = 5;
        const std::vector<std::pair<int, int>> votes = {{0, 2}, {0, 2}, {1, 0}, {1, 0}, {2, 3}, {3, 1}};  // 入力例 1（0 始まり）
        std::vector<ll> cnt(n, 0);
        std::vector<std::vector<int>> by(n);
        for (int i = 0; i < (int)votes.size(); ++i) cnt[votes[i].first]++, by[votes[i].first].push_back(i);
        ArgMinMaxSegTree<MAX_ONLY> st(cnt);
        std::vector<int> ans;
        for (int k = 0; k < n; ++k) {
            for (int i : by[k]) st.add(votes[i].first, -1), st.add(votes[i].second, +1);
            ans.push_back(st.first_max() + 1);
            for (int i : by[k]) st.add(votes[i].first, +1), st.add(votes[i].second, -1);
        }
        note("入力例 1 の答え（1 始まり）→ " + to_s(ans) + "（正解は [3, 1, 1, 2, 1]）");
    }

    api("【7】 例：ABC483 D（A_i を 1 つだけ符号反転して、通る整数座標の数を最大に）", "（使い方の例）",
        "位置の累積和 P_0〜P_N を載せ、A_i の反転 = P_{i+1} 以降に −2A_i を足す → max − min + 1 → 戻す", "全体 O(N log N)");
    {
        const std::vector<ll> A = {-2, 3, -2};  // 入力例 1
        const int n = (int)A.size();
        std::vector<ll> P(n + 1, 0);
        for (int i = 0; i < n; ++i) P[i + 1] = P[i] + A[i];
        ArgMinMaxLazySegTree<BOTH> st(P);
        ll best = st.max() - st.min() + 1;
        for (int i = 0; i < n; ++i) {
            st.range_add(i + 1, n + 1, -2 * A[i]);
            best = std::max(best, st.max() - st.min() + 1);
            st.range_add(i + 1, n + 1, 2 * A[i]);
        }
        note("入力例 1 の答え → " + to_s(best) + "（正解は 8）");
    }
}

}  // namespace amm_demo

// =============================================================================
//  Windows のコンソールで日本語が文字化けしないようにする（Windows 以外では何もしない）
// =============================================================================
#ifdef _WIN32
extern "C" __declspec(dllimport) int __stdcall SetConsoleOutputCP(unsigned int);
#endif
static void setup_console() {
#ifdef _WIN32
    const char* probe = "あ";  // UTF-8 なら 0xE3 0x81 0x82、CP932 なら 0x82 0xA0
    SetConsoleOutputCP((unsigned char)probe[0] == 0xE3 ? 65001u : 932u);
#endif
}

int main() {
    setup_console();
    amm_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    amm_test::benchmark();    // 速度計測
    amm_demo::demo();         // 使い方デモ
    return 0;
}
