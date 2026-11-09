// =============================================================================
//  【ごり押し】2 番目まで持つ最小・最大の位置のセグ木（最小・最大と 2 番目の値・位置・個数の 16 項目を 1 回で。i を除いた最大・重複込みの 2 番目も。1 点更新版・区間更新版）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ ごり押しテンプレート（考察を飛ばして、直さずにそのまま呼べる形）
//    省ける考察     : (最大, 個数, 2 番目, 個数) を合わせるモノイドを自分で書く・「i を除いた最大」を前後の累積から作る
//    向いている問題 : 区間の 2 番目に大きい値とその個数（1 点更新あり。ABC343 F）、a[i] を除いた最大値を i ごとに（ABC134 C）、
//                     1 位が抜けたときの次点、など
//    切り替え       : ArgMinMax2SegTree<MAX_ONLY> / <MIN_ONLY> / <BOTH>（省略すると BOTH）。問題を読んで要る側だけにすると速い
//                     （中に持つものも片側だけになる）。返す struct はどれでも同じ 16 項目で、使わない側は空の値が入る
//    今あるもので書くと : static_seq_toolkit.cpp の second_max / second_max_distinct（更新なし）、MonoidPrefixTemplateTest.txt の except(i)（更新なし）
//    2 番目が要らないとき : ds_segtree_argminmax.cpp（ArgMinMaxSegTree。8 項目だけなので約 1.5〜2 倍速い）
//
//  持つ情報（区間 [l, r) について。個数 = その値をとる位置の数）
//    min  / first_min  / last_min  / min_count    … 最小値・それをとる最も左の位置・最も右の位置・個数
//    min2 / first_min2 / last_min2 / min2_count   … 2 番目に小さい値（同じ値は 1 つと数える：最小より大きい値の中で最小）とその位置・個数
//    max  / first_max  / last_max  / max_count    … 最大値・位置・個数
//    max2 / first_max2 / last_max2 / max2_count   … 2 番目に大きい値（同じ値は 1 つと数える：最大より小さい値の中で最大）とその位置・個数
//    例 [3, 1, 4, 1, 5] の全体：min 1（位置 1・3、2 個）・min2 3（位置 0、1 個）/ max 5（位置 4、1 個）・max2 4（位置 2、1 個）
//    空区間と、使わない側（MAX_ONLY のときの最小側など）：min・min2 = INF（LLONG_MAX）・max・max2 = −INF・位置 −1・個数 0
//    2 番目が無いとき（区間の値が全部同じ・長さ 1）：2 番目の 4 項目だけ空区間と同じ値
//
//  機能一覧（値は long long。位置は 0〜N−1。区間は半開区間 [l, r)）                                        計算量
//    ArgMinMax2SegTree<MODE> seg(a)         1 点更新版。a は vector<long long> でも vector<int> でもよい        O(N)
//    ArgMinMax2SegTree<MODE> seg(n)         長さ n・全部 0                                                     O(N)
//    ArgMinMax2LazySegTree<MODE> seg(a) / seg(n)   区間更新版（下の「区間更新版だけ」の関数も使える）          O(N)
//    共通
//      seg.set(i, x) / seg.add(i, x)        a[i] = x / a[i] += x                                               O(log N)
//      seg.get(i)                           a[i]                                       1 点更新版 O(1)・区間更新版 O(log N)
//      seg.prod(l, r)                       16 項目をまとめて struct で返す（auto s = seg.prod(l, r); s.max2_count）  O(log N)
//      seg.all_prod()                       全体の 16 項目                                                     O(1)
//      seg.min(l, r) / seg.first_min(l, r) / … / seg.max2_count(l, r)   区間の 1 項目だけ（16 個の名前）           O(log N)
//      seg.min() / … / seg.max2_count()     全体の 1 項目だけ                                                  O(1)
//      seg.second_max(l, r) / seg.second_min(l, r)   重複込みの 2 番目（[5, 5, 3] なら 5。max2 は 3）。
//                                           要素が 1 個以下なら −INF / INF。引数なしで全体                       O(log N)・全体 O(1)
//      seg.max_except(i) / seg.min_except(i)          a[i] を除いた全体の最大・最小（ほかに要素が無ければ −INF / INF）  O(1)
//      seg.max_except(l, r, i) / seg.min_except(l, r, i)   [l, r) から a[i] を除いた最大・最小（l ≤ i < r）       O(log N)
//      seg.size() / seg.dump()（a を標準エラーに表示）
//    区間更新版だけ
//      seg.range_add(l, r, x)               [l, r) の各値に x を足す                                           O(log N)
//      seg.assign(l, r, x)                  [l, r) の各値を x にする（その区間の 2 番目は「無い」になる）      O(log N)
//      seg.range_affine(l, r, a, b)         [l, r) の各値を a·x + b にする（a < 0 は BOTH のときだけ。
//                                           符号が反転すると最小側と最大側が入れ替わる。MAX_ONLY / MIN_ONLY で a < 0 は assert で止まる）  O(log N)
//    使わない側の関数（MAX_ONLY での seg.min2() や seg.min_except(i) など）を呼ぶとコンパイルエラー
//
//  速さの違い（どれを呼ぶか迷ったら）
//    1 番目の 4 項目の関数（区間）: 2 番目を計算しないので、どの切り替えでも prod より速い
//    2 番目の 4 項目・second_max・max_except の関数（区間）: BOTH のときだけ prod より速い（片側だけ計算する）
//    全体の関数（引数なし・max_except(i)）: どの切り替えでも all_prod と同じ速さ（O(1)）
//    同じ区間で 2 つ以上の値が要るなら prod を 1 回。区間更新版では 1 項目の関数と prod の差は小さい（遅延の伝搬はどちらも同じだけする）
//    1 点更新しかしないなら 1 点更新版の方が速い（区間更新版の約 1.1〜1.3 倍）
//
//  速さの目安（N = Q = 2×10^5、更新と区間の質問が半分ずつ、prod で 1 項目を読む。計測した環境は AtCoder より遅め）
//    1 点更新版（1 点更新 + 区間）: MAX_ONLY 約 0.10 秒・BOTH 約 0.22 秒（ds_segtree_argminmax.cpp は約 0.06 秒・0.10 秒）
//    区間更新版（区間加算 + 区間）: MAX_ONLY 約 0.23 秒・BOTH 約 0.40 秒（ds_segtree_argminmax.cpp は約 0.15 秒・0.19 秒）
//    10^6 だと 1 点更新版 BOTH 約 1.4 秒・区間更新版 BOTH 約 3.0 秒なので、10^6 級で BOTH が要るときは時間に注意
//    （片側だけで済むなら MAX_ONLY / MIN_ONLY にする。2 番目が要らないなら ds_segtree_argminmax.cpp にする）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい（ds_segtree_argminmax.cpp と両方貼っても衝突しない）
//    amm2_test::random_test()  … 3 つの切り替え × 1 点更新版・区間更新版を、毎回数え直す愚直解と突き合わせる
//    amm2_test::benchmark()    … 速度計測（切り替え・1 項目の関数・2 番目が無い版との比較）
//    amm2_demo::demo()         … 全機能の使い方デモ（日本語の説明付き。ABC343 F・ABC134 C の例つき）
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
#define ARGMINMAX_MODE_DEFINED  // ds_segtree_argminmax.cpp と両方貼っても二重定義にならないための目印
// 最大だけ・最小だけ・両方（テンプレート引数で指定する。省略すると BOTH）
enum ArgMinMaxMode : int { MAX_ONLY = 1, MIN_ONLY = 2, BOTH = 3 };
#endif

namespace amm2_internal {
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
template <bool MX>
inline Side pick(const Side& X, const Side& Y) {  // 空を飛ばして合わせる（X が左、Y が右）
    if (X.c == 0) return Y;
    if (Y.c == 0) return X;
    return merge<MX>(X, Y);
}
// 片側の 8 項目（a = 1 番目、b = 2 番目（同じ値は 1 つと数える））
struct Side2 {
    Side a, b;
};
template <bool MX>
constexpr Side2 empty_side2() { return Side2{empty_side<MX>(), empty_side<MX>()}; }
template <bool MX>
inline Side2 merge2(const Side2& L, const Side2& R) {
    if (L.a.c == 0) return R;
    if (R.a.c == 0) return L;
    if (MX ? L.a.v > R.a.v : L.a.v < R.a.v) return Side2{L.a, pick<MX>(L.b, R.a)};
    if (MX ? R.a.v > L.a.v : R.a.v < L.a.v) return Side2{R.a, pick<MX>(L.a, R.b)};
    return Side2{Side{L.a.v, L.a.f, R.a.l, L.a.c + R.a.c}, pick<MX>(L.b, R.b)};
}
inline void affine_side(Side& s, ll a, ll b) {
    if (s.c) s.v = a * s.v + b;
}
// 重複込みの 2 番目（1 番目が 2 個以上あれば 1 番目、そうでなければ 2 番目）
inline ll second_dup(const Side2& s) { return s.a.c >= 2 ? s.a.v : s.b.v; }
// 位置 i を除いた 1 番目（1 番目が i の 1 個だけなら 2 番目）
inline ll except(const Side2& s, int i) { return (s.a.c == 1 && s.a.f == i) ? s.b.v : s.a.v; }
// 返す struct（どの切り替えでも同じ 16 項目。使わない側は空の値）
struct Info {
    ll min, min2, max, max2;
    int first_min, last_min, min_count, first_min2, last_min2, min2_count;
    int first_max, last_max, max_count, first_max2, last_max2, max2_count;
    friend std::ostream& operator<<(std::ostream& os, const Info& s) {
        return os << "{min " << s.min << " (" << s.first_min << ", " << s.last_min << ", " << s.min_count << ")"
                  << " min2 " << s.min2 << " (" << s.first_min2 << ", " << s.last_min2 << ", " << s.min2_count << ")"
                  << " / max " << s.max << " (" << s.first_max << ", " << s.last_max << ", " << s.max_count << ")"
                  << " max2 " << s.max2 << " (" << s.first_max2 << ", " << s.last_max2 << ", " << s.max2_count << ")}";
    }
};
// 中に持つノード（切り替えで形が変わる：使わない側は持たない）
template <int MODE>
struct Core {
    struct NMax { Side2 mx; };
    struct NMin { Side2 mn; };
    struct NBoth { Side2 mn, mx; };
    using Node = std::conditional_t<MODE == MAX_ONLY, NMax, std::conditional_t<MODE == MIN_ONLY, NMin, NBoth>>;
    static constexpr bool HAS_MIN = (MODE & MIN_ONLY) != 0;
    static constexpr bool HAS_MAX = (MODE & MAX_ONLY) != 0;
    static Node e() {
        Node n;
        if constexpr (HAS_MIN) n.mn = empty_side2<false>();
        if constexpr (HAS_MAX) n.mx = empty_side2<true>();
        return n;
    }
    static Node leaf(int i, ll v) {
        Node n;
        if constexpr (HAS_MIN) n.mn = Side2{Side{v, i, i, 1}, empty_side<false>()};
        if constexpr (HAS_MAX) n.mx = Side2{Side{v, i, i, 1}, empty_side<true>()};
        return n;
    }
    static Node op(const Node& x, const Node& y) {
        Node n;
        if constexpr (HAS_MIN) n.mn = merge2<false>(x.mn, y.mn);
        if constexpr (HAS_MAX) n.mx = merge2<true>(x.mx, y.mx);
        return n;
    }
    template <bool MX>
    static const Side2& side(const Node& n) {
        if constexpr (MX) return n.mx;
        else return n.mn;
    }
    static const Side2& any_side(const Node& n) {
        if constexpr (HAS_MAX) return n.mx;
        else return n.mn;
    }
    static Info to_info(const Node& n) {
        Info s{INF, INF, -INF, -INF, -1, -1, 0, -1, -1, 0, -1, -1, 0, -1, -1, 0};
        if constexpr (HAS_MIN) {
            s.min = n.mn.a.v, s.first_min = n.mn.a.f, s.last_min = n.mn.a.l, s.min_count = n.mn.a.c;
            s.min2 = n.mn.b.v, s.first_min2 = n.mn.b.f, s.last_min2 = n.mn.b.l, s.min2_count = n.mn.b.c;
        }
        if constexpr (HAS_MAX) {
            s.max = n.mx.a.v, s.first_max = n.mx.a.f, s.last_max = n.mx.a.l, s.max_count = n.mx.a.c;
            s.max2 = n.mx.b.v, s.first_max2 = n.mx.b.f, s.last_max2 = n.mx.b.l, s.max2_count = n.mx.b.c;
        }
        return s;
    }
    // 作用 a·x + b（a = 0 は代入）。lo・len はノードが受け持つ区間（代入で使う）
    static void apply(Node& n, ll a, ll b, int lo, int len) {
        if (any_side(n).a.c == 0) return;
        if (a == 0) {  // 代入すると全部同じ値になるので 2 番目は無い
            if constexpr (HAS_MIN) n.mn = Side2{Side{b, lo, lo + len - 1, len}, empty_side<false>()};
            if constexpr (HAS_MAX) n.mx = Side2{Side{b, lo, lo + len - 1, len}, empty_side<true>()};
        } else if (a > 0) {
            if constexpr (HAS_MIN) affine_side(n.mn.a, a, b), affine_side(n.mn.b, a, b);
            if constexpr (HAS_MAX) affine_side(n.mx.a, a, b), affine_side(n.mx.b, a, b);
        } else {
            if constexpr (MODE == BOTH) {
                std::swap(n.mn, n.mx);  // 符号が反転すると最小側と最大側が入れ替わる
                affine_side(n.mn.a, a, b), affine_side(n.mn.b, a, b), affine_side(n.mx.a, a, b), affine_side(n.mx.b, a, b);
                if (n.mn.b.c == 0) n.mn.b = empty_side<false>();  // 空の 2 番目の値（±INF）の向きを直す
                if (n.mx.b.c == 0) n.mx.b = empty_side<true>();
            }
        }
    }
};
}  // namespace amm2_internal

// 1 点更新版
template <int MODE = BOTH>
struct ArgMinMax2SegTree {
    static_assert(MODE == MAX_ONLY || MODE == MIN_ONLY || MODE == BOTH, "MODE は MAX_ONLY / MIN_ONLY / BOTH");
    using S = amm2_internal::Info;
    static constexpr long long INF = amm2_internal::INF;

    ArgMinMax2SegTree() : ArgMinMax2SegTree(0) {}
    explicit ArgMinMax2SegTree(int n) : ArgMinMax2SegTree(std::vector<long long>(n, 0)) {}
    template <class V>
    explicit ArgMinMax2SegTree(const std::vector<V>& a) : n_((int)a.size()) {
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
        return C::any_side(d_[i + sz_]).a.v;
    }

    // [l, r) の 16 項目をまとめて返す。1 項目だけ要るときは 1 項目の関数の方が速い（1 番目の項目はどの切り替えでも、2 番目の項目は BOTH のとき）。
    // 2 つ以上の値が要るならこちらを 1 回
    S prod(int l, int r) const {
        assert(0 <= l && l <= r && r <= n_);
        Node sl = C::e(), sr = C::e();
        for (l += sz_, r += sz_; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sl = C::op(sl, d_[l++]);
            if (r & 1) sr = C::op(d_[--r], sr);
        }
        return C::to_info(C::op(sl, sr));
    }
    // 全体の 16 項目。1 項目の関数の全体版と同じ速さ（O(1)）
    S all_prod() const { return C::to_info(d_[1]); }

    // ---- 1 番目の 4 項目の関数（(l, r) で区間 [l, r)・引数なしで全体）。区間版は 2 番目を計算しないので、どの切り替えでも prod より速い
    // [l, r) の最小値（空なら INF）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long min(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).v; }
    // 全体の最小値（空なら INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.a.v; }
    // [l, r) の最小値をとる最も左の位置（空なら −1）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_min(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).f; }
    // 全体の最小値をとる最も左の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.a.f; }
    // [l, r) の最小値をとる最も右の位置（空なら −1）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_min(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).l; }
    // 全体の最小値をとる最も右の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.a.l; }
    // [l, r) の最小値の個数（空なら 0）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int min_count(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).c; }
    // 全体の最小値の個数（空なら 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int min_count() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.a.c; }
    // [l, r) の最大値（空なら −INF）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long max(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).v; }
    // 全体の最大値（空なら −INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.a.v; }
    // [l, r) の最大値をとる最も左の位置（空なら −1）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_max(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).f; }
    // 全体の最大値をとる最も左の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.a.f; }
    // [l, r) の最大値をとる最も右の位置（空なら −1）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_max(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).l; }
    // 全体の最大値をとる最も右の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.a.l; }
    // [l, r) の最大値の個数（空なら 0）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int max_count(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).c; }
    // 全体の最大値の個数（空なら 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int max_count() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.a.c; }

    // ---- 2 番目の 4 項目の関数。区間版は BOTH のときだけ prod より速い（片側だけ計算する）
    // [l, r) の2 番目に小さい値（同じ値は 1 つと数える。無ければ INF）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long min2(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side2_prod<false>(l, r).b.v; }
    // 全体の2 番目に小さい値（同じ値は 1 つと数える。無ければ INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long min2() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.b.v; }
    // [l, r) の2 番目に小さい値をとる最も左の位置（無ければ −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_min2(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side2_prod<false>(l, r).b.f; }
    // 全体の2 番目に小さい値をとる最も左の位置（無ければ −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_min2() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.b.f; }
    // [l, r) の2 番目に小さい値をとる最も右の位置（無ければ −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_min2(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side2_prod<false>(l, r).b.l; }
    // 全体の2 番目に小さい値をとる最も右の位置（無ければ −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_min2() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.b.l; }
    // [l, r) の2 番目に小さい値の個数（無ければ 0）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int min2_count(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side2_prod<false>(l, r).b.c; }
    // 全体の2 番目に小さい値の個数（無ければ 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int min2_count() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.b.c; }
    // [l, r) の2 番目に大きい値（同じ値は 1 つと数える。無ければ −INF）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long max2(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side2_prod<true>(l, r).b.v; }
    // 全体の2 番目に大きい値（同じ値は 1 つと数える。無ければ −INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long max2() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.b.v; }
    // [l, r) の2 番目に大きい値をとる最も左の位置（無ければ −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_max2(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side2_prod<true>(l, r).b.f; }
    // 全体の2 番目に大きい値をとる最も左の位置（無ければ −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_max2() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.b.f; }
    // [l, r) の2 番目に大きい値をとる最も右の位置（無ければ −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_max2(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side2_prod<true>(l, r).b.l; }
    // 全体の2 番目に大きい値をとる最も右の位置（無ければ −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_max2() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.b.l; }
    // [l, r) の2 番目に大きい値の個数（無ければ 0）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int max2_count(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side2_prod<true>(l, r).b.c; }
    // 全体の2 番目に大きい値の個数（無ければ 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int max2_count() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.b.c; }

    // ---- 重複込みの 2 番目・i を除いた最大値 / 最小値（値だけを返す）
    // [l, r) で重複込みの 2 番目に小さい値（[1, 1, 3] なら 1。要素が 1 個以下なら INF）。速さは min2(l, r) と同じ
    long long second_min(int l, int r) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return amm2_internal::second_dup(side2_prod<false>(l, r)); }
    // 全体で重複込みの 2 番目に小さい値（[1, 1, 3] なら 1。要素が 1 個以下なら INF）。O(1)
    long long second_min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return amm2_internal::second_dup(d_[1].mn); }
    // a[i] を除いた全体の最小値（ほかに要素が無ければ INF）。O(1)
    long long min_except(int i) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); assert(0 <= i && i < n_); return amm2_internal::except(d_[1].mn, i); }
    // [l, r) から a[i] を除いた最小値（l ≤ i < r。ほかに要素が無ければ INF）。速さは min2(l, r) と同じ
    long long min_except(int l, int r, int i) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); assert(l <= i && i < r); return amm2_internal::except(side2_prod<false>(l, r), i); }
    // [l, r) で重複込みの 2 番目に大きい値（[5, 5, 3] なら 5。要素が 1 個以下なら −INF）。速さは max2(l, r) と同じ
    long long second_max(int l, int r) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return amm2_internal::second_dup(side2_prod<true>(l, r)); }
    // 全体で重複込みの 2 番目に大きい値（[5, 5, 3] なら 5。要素が 1 個以下なら −INF）。O(1)
    long long second_max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return amm2_internal::second_dup(d_[1].mx); }
    // a[i] を除いた全体の最大値（ほかに要素が無ければ −INF）。O(1)
    long long max_except(int i) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); assert(0 <= i && i < n_); return amm2_internal::except(d_[1].mx, i); }
    // [l, r) から a[i] を除いた最大値（l ≤ i < r。ほかに要素が無ければ −INF）。速さは max2(l, r) と同じ
    long long max_except(int l, int r, int i) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); assert(l <= i && i < r); return amm2_internal::except(side2_prod<true>(l, r), i); }

    // a を標準エラーに表示（デバッグ用）
    void dump() const {
        std::cerr << "a = [";
        for (int i = 0; i < n_; ++i) std::cerr << (i ? ", " : "") << get(i);
        std::cerr << "]\n";
    }

   private:
    using C = amm2_internal::Core<MODE>;
    using Node = typename C::Node;
    using Side = amm2_internal::Side;
    using Side2 = amm2_internal::Side2;
    int n_ = 0, sz_ = 1;
    std::vector<Node> d_;
    template <bool MX>
    Side side_prod(int l, int r) const {  // 1 番目だけ
        assert(0 <= l && l <= r && r <= n_);
        Side sl = amm2_internal::empty_side<MX>(), sr = sl;
        for (l += sz_, r += sz_; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sl = amm2_internal::merge<MX>(sl, C::template side<MX>(d_[l++]).a);
            if (r & 1) sr = amm2_internal::merge<MX>(C::template side<MX>(d_[--r]).a, sr);
        }
        return amm2_internal::merge<MX>(sl, sr);
    }
    template <bool MX>
    Side2 side2_prod(int l, int r) const {  // 片側の 1 番目と 2 番目
        assert(0 <= l && l <= r && r <= n_);
        Side2 sl = amm2_internal::empty_side2<MX>(), sr = sl;
        for (l += sz_, r += sz_; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sl = amm2_internal::merge2<MX>(sl, C::template side<MX>(d_[l++]));
            if (r & 1) sr = amm2_internal::merge2<MX>(C::template side<MX>(d_[--r]), sr);
        }
        return amm2_internal::merge2<MX>(sl, sr);
    }
};

// 区間更新版（作用は a·x + b。加算 a = 1・代入 a = 0・符号反転 a = −1 はその特別な場合）
template <int MODE = BOTH>
struct ArgMinMax2LazySegTree {
    static_assert(MODE == MAX_ONLY || MODE == MIN_ONLY || MODE == BOTH, "MODE は MAX_ONLY / MIN_ONLY / BOTH");
    using S = amm2_internal::Info;
    static constexpr long long INF = amm2_internal::INF;

    ArgMinMax2LazySegTree() : ArgMinMax2LazySegTree(0) {}
    explicit ArgMinMax2LazySegTree(int n) : ArgMinMax2LazySegTree(std::vector<long long>(n, 0)) {}
    template <class V>
    explicit ArgMinMax2LazySegTree(const std::vector<V>& a) : n_((int)a.size()) {
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
        return C::any_side(d_[k]).a.v;
    }

    // [l, r) の各値に x を足す
    void range_add(int l, int r, long long x) { apply_range(l, r, F{1, x}); }
    // [l, r) の各値を x にする（その区間の 2 番目は「無い」になる）
    void assign(int l, int r, long long x) { apply_range(l, r, F{0, x}); }
    // [l, r) の各値を a·x + b にする（a < 0 は BOTH のときだけ）
    void range_affine(int l, int r, long long a, long long b) {
        assert((a >= 0 || MODE == BOTH) && "負の数を掛けるのは BOTH のときだけ（片側だけでは掛けたあとの最大・最小が分からない）");
        apply_range(l, r, F{a, b});
    }

    // [l, r) の 16 項目をまとめて返す。1 項目だけ要るときは 1 項目の関数の方が速い（1 番目の項目はどの切り替えでも、2 番目の項目は BOTH のとき）。
    // 2 つ以上の値が要るならこちらを 1 回
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
    // 全体の 16 項目。1 項目の関数の全体版と同じ速さ（O(1)）
    S all_prod() const { return C::to_info(d_[1]); }

    // ---- 1 番目の 4 項目の関数（(l, r) で区間 [l, r)・引数なしで全体）。区間版は 2 番目を計算しないので、どの切り替えでも prod より速い
    // [l, r) の最小値（空なら INF）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long min(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).v; }
    // 全体の最小値（空なら INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.a.v; }
    // [l, r) の最小値をとる最も左の位置（空なら −1）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_min(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).f; }
    // 全体の最小値をとる最も左の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.a.f; }
    // [l, r) の最小値をとる最も右の位置（空なら −1）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_min(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).l; }
    // 全体の最小値をとる最も右の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.a.l; }
    // [l, r) の最小値の個数（空なら 0）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int min_count(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side_prod<false>(l, r).c; }
    // 全体の最小値の個数（空なら 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int min_count() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.a.c; }
    // [l, r) の最大値（空なら −INF）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long max(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).v; }
    // 全体の最大値（空なら −INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.a.v; }
    // [l, r) の最大値をとる最も左の位置（空なら −1）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_max(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).f; }
    // 全体の最大値をとる最も左の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.a.f; }
    // [l, r) の最大値をとる最も右の位置（空なら −1）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_max(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).l; }
    // 全体の最大値をとる最も右の位置（空なら −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.a.l; }
    // [l, r) の最大値の個数（空なら 0）。2 番目を計算しないので、どの切り替えでも prod より速い（BOTH なら片側だけなのでさらに速い）。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int max_count(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side_prod<true>(l, r).c; }
    // 全体の最大値の個数（空なら 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int max_count() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.a.c; }

    // ---- 2 番目の 4 項目の関数。区間版は BOTH のときだけ prod より速い（片側だけ計算する）
    // [l, r) の2 番目に小さい値（同じ値は 1 つと数える。無ければ INF）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long min2(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side2_prod<false>(l, r).b.v; }
    // 全体の2 番目に小さい値（同じ値は 1 つと数える。無ければ INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long min2() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.b.v; }
    // [l, r) の2 番目に小さい値をとる最も左の位置（無ければ −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_min2(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side2_prod<false>(l, r).b.f; }
    // 全体の2 番目に小さい値をとる最も左の位置（無ければ −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_min2() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.b.f; }
    // [l, r) の2 番目に小さい値をとる最も右の位置（無ければ −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_min2(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side2_prod<false>(l, r).b.l; }
    // 全体の2 番目に小さい値をとる最も右の位置（無ければ −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_min2() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.b.l; }
    // [l, r) の2 番目に小さい値の個数（無ければ 0）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int min2_count(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return side2_prod<false>(l, r).b.c; }
    // 全体の2 番目に小さい値の個数（無ければ 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int min2_count() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return d_[1].mn.b.c; }
    // [l, r) の2 番目に大きい値（同じ値は 1 つと数える。無ければ −INF）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    long long max2(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side2_prod<true>(l, r).b.v; }
    // 全体の2 番目に大きい値（同じ値は 1 つと数える。無ければ −INF）。どの切り替えでも all_prod と同じ速さ（O(1)）
    long long max2() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.b.v; }
    // [l, r) の2 番目に大きい値をとる最も左の位置（無ければ −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int first_max2(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side2_prod<true>(l, r).b.f; }
    // 全体の2 番目に大きい値をとる最も左の位置（無ければ −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int first_max2() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.b.f; }
    // [l, r) の2 番目に大きい値をとる最も右の位置（無ければ −1）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int last_max2(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side2_prod<true>(l, r).b.l; }
    // 全体の2 番目に大きい値をとる最も右の位置（無ければ −1）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int last_max2() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.b.l; }
    // [l, r) の2 番目に大きい値の個数（無ければ 0）。BOTH のときだけ prod より速い（片側だけ計算する）。MAX_ONLY / MIN_ONLY では prod と同じ速さ。同じ区間で 2 つ以上の値が要るなら prod を 1 回
    int max2_count(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return side2_prod<true>(l, r).b.c; }
    // 全体の2 番目に大きい値の個数（無ければ 0）。どの切り替えでも all_prod と同じ速さ（O(1)）
    int max2_count() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return d_[1].mx.b.c; }

    // ---- 重複込みの 2 番目・i を除いた最大値 / 最小値（値だけを返す）
    // [l, r) で重複込みの 2 番目に小さい値（[1, 1, 3] なら 1。要素が 1 個以下なら INF）。速さは min2(l, r) と同じ
    long long second_min(int l, int r) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return amm2_internal::second_dup(side2_prod<false>(l, r)); }
    // 全体で重複込みの 2 番目に小さい値（[1, 1, 3] なら 1。要素が 1 個以下なら INF）。O(1)
    long long second_min() const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); return amm2_internal::second_dup(d_[1].mn); }
    // a[i] を除いた全体の最小値（ほかに要素が無ければ INF）。O(1)
    long long min_except(int i) const { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); assert(0 <= i && i < n_); return amm2_internal::except(d_[1].mn, i); }
    // [l, r) から a[i] を除いた最小値（l ≤ i < r。ほかに要素が無ければ INF）。速さは min2(l, r) と同じ
    long long min_except(int l, int r, int i) { static_assert(C::HAS_MIN, "最小側は MIN_ONLY / BOTH で使う"); assert(l <= i && i < r); return amm2_internal::except(side2_prod<false>(l, r), i); }
    // [l, r) で重複込みの 2 番目に大きい値（[5, 5, 3] なら 5。要素が 1 個以下なら −INF）。速さは max2(l, r) と同じ
    long long second_max(int l, int r) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return amm2_internal::second_dup(side2_prod<true>(l, r)); }
    // 全体で重複込みの 2 番目に大きい値（[5, 5, 3] なら 5。要素が 1 個以下なら −INF）。O(1)
    long long second_max() const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); return amm2_internal::second_dup(d_[1].mx); }
    // a[i] を除いた全体の最大値（ほかに要素が無ければ −INF）。O(1)
    long long max_except(int i) const { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); assert(0 <= i && i < n_); return amm2_internal::except(d_[1].mx, i); }
    // [l, r) から a[i] を除いた最大値（l ≤ i < r。ほかに要素が無ければ −INF）。速さは max2(l, r) と同じ
    long long max_except(int l, int r, int i) { static_assert(C::HAS_MAX, "最大側は MAX_ONLY / BOTH で使う"); assert(l <= i && i < r); return amm2_internal::except(side2_prod<true>(l, r), i); }

    // a を標準エラーに表示（デバッグ用）
    void dump() {
        std::cerr << "a = [";
        for (int i = 0; i < n_; ++i) std::cerr << (i ? ", " : "") << get(i);
        std::cerr << "]\n";
    }

   private:
    using C = amm2_internal::Core<MODE>;
    using Node = typename C::Node;
    using Side = amm2_internal::Side;
    using Side2 = amm2_internal::Side2;
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
    Side side_prod(int l, int r) {  // 1 番目だけ
        assert(0 <= l && l <= r && r <= n_);
        Side sl = amm2_internal::empty_side<MX>(), sr = sl;
        if (l == r) return sl;
        push_range(l + sz_, r + sz_);
        for (l += sz_, r += sz_; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sl = amm2_internal::merge<MX>(sl, C::template side<MX>(d_[l++]).a);
            if (r & 1) sr = amm2_internal::merge<MX>(C::template side<MX>(d_[--r]).a, sr);
        }
        return amm2_internal::merge<MX>(sl, sr);
    }
    template <bool MX>
    Side2 side2_prod(int l, int r) {  // 片側の 1 番目と 2 番目
        assert(0 <= l && l <= r && r <= n_);
        Side2 sl = amm2_internal::empty_side2<MX>(), sr = sl;
        if (l == r) return sl;
        push_range(l + sz_, r + sz_);
        for (l += sz_, r += sz_; l < r; l >>= 1, r >>= 1) {
            if (l & 1) sl = amm2_internal::merge2<MX>(sl, C::template side<MX>(d_[l++]));
            if (r & 1) sr = amm2_internal::merge2<MX>(C::template side<MX>(d_[--r]), sr);
        }
        return amm2_internal::merge2<MX>(sl, sr);
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace amm2_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define AMM2_CHECK(...)                                                                                  \
    do {                                                                                                 \
        if (!(__VA_ARGS__)) {                                                                            \
            std::cerr << "[ds_segtree_argminmax2] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                                \
        }                                                                                                \
    } while (0)

using ll = long long;
using Info = amm2_internal::Info;
constexpr ll INF = amm2_internal::INF;

// [lo, hi] の一様乱数
[[maybe_unused]] static ll rnd(std::mt19937_64& rng, ll lo, ll hi) {
    const unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return (ll)((unsigned long long)lo + (w == 0 ? rng() : rng() % w));
}
[[maybe_unused]] static long long elapsed_ms(std::chrono::steady_clock::time_point a) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
}

// 愚直：[l, r) の 16 項目を数え直す（使わない側は空の値）
static Info brute(const std::vector<ll>& a, int l, int r, int mode) {
    Info s{INF, INF, -INF, -INF, -1, -1, 0, -1, -1, 0, -1, -1, 0, -1, -1, 0};
    auto scan = [&](ll v, int& f, int& la, int& c) {  // 値 v の位置と個数
        for (int i = l; i < r; ++i)
            if (a[i] == v) {
                if (c == 0) f = i;
                la = i, ++c;
            }
    };
    if (l == r) return s;
    if (mode & MIN_ONLY) {
        s.min = *std::min_element(a.begin() + l, a.begin() + r);
        scan(s.min, s.first_min, s.last_min, s.min_count);
        for (int i = l; i < r; ++i)
            if (a[i] > s.min && a[i] < s.min2) s.min2 = a[i];
        if (s.min2 != INF) scan(s.min2, s.first_min2, s.last_min2, s.min2_count);
    }
    if (mode & MAX_ONLY) {
        s.max = *std::max_element(a.begin() + l, a.begin() + r);
        scan(s.max, s.first_max, s.last_max, s.max_count);
        for (int i = l; i < r; ++i)
            if (a[i] < s.max && a[i] > s.max2) s.max2 = a[i];
        if (s.max2 != -INF) scan(s.max2, s.first_max2, s.last_max2, s.max2_count);
    }
    return s;
}
static bool same(const Info& x, const Info& y) {
    return x.min == y.min && x.first_min == y.first_min && x.last_min == y.last_min && x.min_count == y.min_count && x.min2 == y.min2 &&
           x.first_min2 == y.first_min2 && x.last_min2 == y.last_min2 && x.min2_count == y.min2_count && x.max == y.max &&
           x.first_max == y.first_max && x.last_max == y.last_max && x.max_count == y.max_count && x.max2 == y.max2 &&
           x.first_max2 == y.first_max2 && x.last_max2 == y.last_max2 && x.max2_count == y.max2_count;
}
// 愚直：重複込みの 2 番目・i を除いた最大 / 最小
static ll brute_second(const std::vector<ll>& a, int l, int r, bool mx) {
    std::vector<ll> v(a.begin() + l, a.begin() + r);
    if (v.size() < 2) return mx ? -INF : INF;
    std::sort(v.begin(), v.end());
    return mx ? v[v.size() - 2] : v[1];
}
static ll brute_except(const std::vector<ll>& a, int l, int r, int i, bool mx) {
    ll s = mx ? -INF : INF;
    for (int j = l; j < r; ++j)
        if (j != i) s = mx ? std::max(s, a[j]) : std::min(s, a[j]);
    return s;
}

// 1 項目の関数（区間版・全体版）と second / except が愚直と同じか
template <int MODE, class Seg>
static void check_single(Seg& seg, const std::vector<ll>& a, int l, int r, const Info& want) {
    const int n = seg.size();
    const Info all = seg.all_prod();
    std::mt19937_64 rng((unsigned)(l * 1000 + r));
    if constexpr ((MODE & MIN_ONLY) != 0) {
        AMM2_CHECK(seg.min(l, r) == want.min && seg.first_min(l, r) == want.first_min && seg.last_min(l, r) == want.last_min && seg.min_count(l, r) == want.min_count);
        AMM2_CHECK(seg.min2(l, r) == want.min2 && seg.first_min2(l, r) == want.first_min2 && seg.last_min2(l, r) == want.last_min2 && seg.min2_count(l, r) == want.min2_count);
        AMM2_CHECK(seg.min() == all.min && seg.first_min() == all.first_min && seg.last_min() == all.last_min && seg.min_count() == all.min_count);
        AMM2_CHECK(seg.min2() == all.min2 && seg.first_min2() == all.first_min2 && seg.last_min2() == all.last_min2 && seg.min2_count() == all.min2_count);
        AMM2_CHECK(seg.second_min(l, r) == brute_second(a, l, r, false) && seg.second_min() == brute_second(a, 0, n, false));
        if (l < r) {
            const int i = (int)rnd(rng, l, r - 1);
            AMM2_CHECK(seg.min_except(l, r, i) == brute_except(a, l, r, i, false));
        }
        if (n > 0) {
            const int i = (int)rnd(rng, 0, n - 1);
            AMM2_CHECK(seg.min_except(i) == brute_except(a, 0, n, i, false));
        }
    }
    if constexpr ((MODE & MAX_ONLY) != 0) {
        AMM2_CHECK(seg.max(l, r) == want.max && seg.first_max(l, r) == want.first_max && seg.last_max(l, r) == want.last_max && seg.max_count(l, r) == want.max_count);
        AMM2_CHECK(seg.max2(l, r) == want.max2 && seg.first_max2(l, r) == want.first_max2 && seg.last_max2(l, r) == want.last_max2 && seg.max2_count(l, r) == want.max2_count);
        AMM2_CHECK(seg.max() == all.max && seg.first_max() == all.first_max && seg.last_max() == all.last_max && seg.max_count() == all.max_count);
        AMM2_CHECK(seg.max2() == all.max2 && seg.first_max2() == all.first_max2 && seg.last_max2() == all.last_max2 && seg.max2_count() == all.max2_count);
        AMM2_CHECK(seg.second_max(l, r) == brute_second(a, l, r, true) && seg.second_max() == brute_second(a, 0, n, true));
        if (l < r) {
            const int i = (int)rnd(rng, l, r - 1);
            AMM2_CHECK(seg.max_except(l, r, i) == brute_except(a, l, r, i, true));
        }
        if (n > 0) {
            const int i = (int)rnd(rng, 0, n - 1);
            AMM2_CHECK(seg.max_except(i) == brute_except(a, 0, n, i, true));
        }
    }
}

template <int MODE, bool LAZY>
static long long random_test_one(std::mt19937_64& rng) {
    long long cases = 0;
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 14);
        const ll V = rnd(rng, 1, 4);  // 値の範囲を狭くして同じ値をたくさん作る
        std::vector<ll> a(n);
        for (auto& x : a) x = rnd(rng, -V, V);
        std::conditional_t<LAZY, ArgMinMax2LazySegTree<MODE>, ArgMinMax2SegTree<MODE>> seg(a);
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
                    const ll ca = rnd(rng, (MODE == BOTH) ? -1 : 0, 2), cb = rnd(rng, -2, 2);
                    seg.range_affine(l, r, ca, cb);
                    for (int i = l; i < r; ++i) a[i] = ca * a[i] + cb;
                }
            }
            const auto want = brute(a, l, r, MODE);
            AMM2_CHECK(same(seg.prod(l, r), want));
            AMM2_CHECK(same(seg.all_prod(), brute(a, 0, n, MODE)));
            check_single<MODE>(seg, a, l, r, want);
            if (n > 0) {
                const int i = (int)rnd(rng, 0, n - 1);
                AMM2_CHECK(seg.get(i) == a[i]);
            }
            ++cases;
        }
        AMM2_CHECK(seg.size() == n);
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
        ArgMinMax2SegTree<BOTH> e(0);
        AMM2_CHECK(e.all_prod().max2 == -INF && e.all_prod().first_max2 == -1 && e.all_prod().min2_count == 0 && e.second_max() == -INF);
        ArgMinMax2LazySegTree<BOTH> one(std::vector<int>{7});
        one.range_affine(0, 1, -1, 0);
        AMM2_CHECK(one.max() == -7 && one.max2() == -INF && one.min2() == INF && one.max_except(0) == -INF && one.second_min() == INF);
        ArgMinMax2LazySegTree<MAX_ONLY> same_v(std::vector<ll>(10, 5));
        AMM2_CHECK(same_v.max2() == -INF && same_v.second_max() == 5 && same_v.max_except(3) == 5);
        same_v.assign(2, 7, 9);
        AMM2_CHECK(same_v.max() == 9 && same_v.max2() == 5 && same_v.first_max2() == 0 && same_v.last_max2() == 9 && same_v.max2_count() == 5);
        AMM2_CHECK(same_v.all_prod().min2 == INF && same_v.all_prod().first_min2 == -1);  // 使わない側は空の値
        ArgMinMax2SegTree<MIN_ONLY> big(std::vector<ll>{(ll)4e18, (ll)-4e18, (ll)-4e18});
        AMM2_CHECK(big.min() == (ll)-4e18 && big.min2() == (ll)4e18 && big.min_except(1) == (ll)-4e18 && big.second_min() == (ll)-4e18);
        cases += 4;
    }
    std::cout << "[ds_segtree_argminmax2] random_test: OK (" << cases << " cases)\n";
}

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
    if (TEST_SCALE == 1) std::cout << "[ds_segtree_argminmax2] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000, Q = N;
    std::vector<ll> a(N);
    for (auto& x : a) x = rnd(rng, 0, 999);
    const auto pq = make_queries(rng, N, Q, false), rq = make_queries(rng, N, Q, true);
    std::cout << "[ds_segtree_argminmax2] benchmark: N = Q = " << N << "（更新と区間の質問が半分ずつ）\n";
    auto line = [](const std::string& name, long long ms, ll chk) {
        std::cout << "  " << name << " " << ms << " ms  (check " << chk % 1000 << ")\n";
    };
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMax2SegTree<MAX_ONLY> s(a);
        ll chk = 0;
        for (auto& e : pq) e.t ? s.set(e.l, e.x) : (void)(chk += s.prod(e.l, e.r).max2_count);
        line("1 点更新版 MAX_ONLY（1 点更新 + prod）                           ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMax2SegTree<MAX_ONLY> s(a);
        ll chk = 0;
        for (auto& e : pq) e.t ? s.set(e.l, e.x) : (void)(chk += s.max2_count(e.l, e.r));
        line("1 点更新版 MAX_ONLY（1 点更新 + max2_count(l, r)）              ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMax2SegTree<BOTH> s(a);
        ll chk = 0;
        for (auto& e : pq) e.t ? s.set(e.l, e.x) : (void)(chk += s.prod(e.l, e.r).max2_count);
        line("1 点更新版 BOTH（1 点更新 + prod）                               ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMax2SegTree<BOTH> s(a);
        ll chk = 0;
        for (auto& e : pq) e.t ? s.set(e.l, e.x) : (void)(chk += s.max2_count(e.l, e.r));
        line("1 点更新版 BOTH（1 点更新 + max2_count(l, r)）                  ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMax2LazySegTree<BOTH> s(a);
        ll chk = 0;
        for (auto& e : pq) e.t ? s.set(e.l, e.x) : (void)(chk += s.prod(e.l, e.r).max2_count);
        line("区間更新版 BOTH で 1 点更新だけ（1 点更新 + prod）             ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMax2LazySegTree<MAX_ONLY> s(a);
        ll chk = 0;
        for (auto& e : rq) e.t ? s.range_add(e.l, e.r, e.x) : (void)(chk += s.prod(e.l, e.r).max2_count);
        line("区間更新版 MAX_ONLY（区間加算 + prod）                         ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMax2LazySegTree<BOTH> s(a);
        ll chk = 0;
        for (auto& e : rq) e.t ? s.range_add(e.l, e.r, e.x) : (void)(chk += s.prod(e.l, e.r).max2_count);
        line("区間更新版 BOTH（区間加算 + prod）                             ", elapsed_ms(t0), chk);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        ArgMinMax2LazySegTree<BOTH> s(a);
        ll chk = 0;
        int k = 0;
        for (auto& e : rq) {
            if (e.t) {
                ++k;
                if (k % 3 == 0) s.range_add(e.l, e.r, e.x);
                else if (k % 3 == 1) s.assign(e.l, e.r, e.x + 500);
                else s.range_affine(e.l, e.r, -1, 0);
            } else {
                chk += s.prod(e.l, e.r).max2_count;
            }
        }
        line("区間更新版 BOTH（加算・代入・符号反転 + prod）                 ", elapsed_ms(t0), chk);
    }
}

}  // namespace amm2_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace amm2_demo {

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
#define SHOW(expr, comment) amm2_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    amm2_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    using ll = long long;
    banner("【ごり押し】2 番目まで持つ最小・最大の位置のセグ木（ds_segtree_argminmax2） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 区間 [l, r) の最小・最大と 2 番目（同じ値は 1 つと数える）の値・最も左・最も右の位置・個数（16 項目）を 1 回で取る\n"
              << "■ 切り替え: <MAX_ONLY> / <MIN_ONLY> / <BOTH>（省略すると BOTH）。要る側だけにすると速い\n"
              << "■ API 一覧\n"
              << "  ArgMinMax2SegTree<MODE> seg(a) / seg(n): set / add / get / prod / all_prod / 1 項目の関数 16 個（(l, r) で区間・() で全体）\n"
              << "      / second_max / second_min（重複込みの 2 番目）/ max_except / min_except（i を除いた最大・最小）/ size / dump\n"
              << "  ArgMinMax2LazySegTree<MODE> seg(a) / seg(n): 上と同じ + range_add / assign / range_affine\n";

    api("【1】 ArgMinMax2SegTree<BOTH> seg(a)  prod(l, r) / all_prod()", "struct（min / first_min / last_min / min_count / min2 / … / max2_count の 16 項目）",
        "区間 [l, r) の 16 項目をまとめて（auto s = seg.prod(l, r); で受けて s.max2 などで読む）。表示は 値 (最も左, 最も右, 個数)", "O(log N) / O(1)");
    RUN(ArgMinMax2SegTree<BOTH> seg(std::vector<ll>{3, 1, 4, 1, 5, 9, 2, 6, 5}));
    SHOW(seg.all_prod(), "全体");
    RUN(auto s = seg.prod(0, 5));
    SHOW(s.max2, "[3, 1, 4, 1, 5] の 2 番目に大きい値");
    SHOW(s.first_max2, "その位置");
    SHOW(s.min2, "2 番目に小さい値（1 は 2 個あるが、同じ値は 1 つと数える）");
    SHOW(s.min2_count, "");
    SHOW(seg.prod(4, 4), "空区間");
    SHOW(seg.prod(3, 4).max2, "長さ 1 の区間には 2 番目が無い → −INF");

    api("【2】 1 項目の関数（16 個）  min / first_min / … / max2_count", "long long（値）/ int（位置・個数）",
        "(l, r) で区間 [l, r)、引数なしで全体。1 番目の 4 項目は 2 番目を計算しないので prod より速い", "区間 O(log N) / 全体 O(1)");
    SHOW(seg.max2(), "全体の 2 番目に大きい値");
    SHOW(seg.first_max2(), "6 の位置");
    SHOW(seg.max2(5, 9), "[9, 2, 6, 5] の 2 番目に大きい値");
    SHOW(seg.max(5, 9), "1 番目（2 番目を計算しない）");
    SHOW(seg.min2_count(), "全体の 2 番目に小さい値 2 の個数");

    api("【3】 second_max(l, r) / second_min(l, r)（引数なしで全体）", "long long",
        "重複込みの 2 番目（[5, 5, 3] なら 5。max2 は 3）。要素が 1 個以下なら −INF / INF", "O(log N)・全体 O(1)");
    SHOW(seg.second_min(), "1 が 2 個あるので重複込みの 2 番目に小さい値は 1");
    SHOW(seg.min2(), "同じ値は 1 つと数える 2 番目は 2");
    SHOW(seg.second_max(4, 9), "[5, 9, 2, 6, 5] → 6");

    api("【4】 max_except(i) / min_except(i)、max_except(l, r, i) / min_except(l, r, i)", "long long",
        "a[i] を除いた全体（[l, r)）の最大・最小。ほかに要素が無ければ −INF / INF", "全体 O(1)・区間 O(log N)");
    SHOW(seg.max_except(5), "9 を除くと 6");
    SHOW(seg.max_except(0), "3 を除いても 9");
    SHOW(seg.min_except(1), "1 は位置 3 にもあるので 1");
    SHOW(seg.max_except(4, 7, 5), "[5, 9, 2] から 9 を除くと 5");

    api("【5】 set(i, x) / add(i, x) / get(i)", "void / void / long long", "a[i] = x / a[i] += x / a[i]", "O(log N) / O(log N) / 1 点更新版 O(1)");
    RUN(seg.set(7, 9));
    SHOW(seg.max_count(), "9 が 2 個");
    SHOW(seg.max2(), "2 番目は 5");
    SHOW(seg.max_except(5), "9 は位置 7 にもあるので 9");
    RUN(seg.add(7, -8));
    SHOW(seg.get(7), "");

    api("【6】 切り替え  ArgMinMax2SegTree<MAX_ONLY> / <MIN_ONLY>", "（BOTH と同じ struct）",
        "要る側だけ持つので速い。struct の使わない側は空の値。使わない側の関数はコンパイルエラー", "（BOTH より速い）");
    RUN(ArgMinMax2SegTree<MAX_ONLY> mx(std::vector<int>{2, 7, 7, 1}));
    SHOW(mx.max2(), "");
    SHOW(mx.all_prod(), "最小側は空の値");
    note("mx.min2() や mx.min_except(0) と書くとコンパイルエラー（static_assert：最小側は MIN_ONLY / BOTH で使う）");

    api("【7】 ArgMinMax2LazySegTree<BOTH>  range_add(l, r, x) / assign(l, r, x) / range_affine(l, r, a, b)", "void",
        "[l, r) に x を足す / x にする / a·x + b にする（a < 0 は BOTH だけ。符号が反転すると最小側と最大側が入れ替わる）", "O(log N)");
    RUN(ArgMinMax2LazySegTree<BOTH> lz(std::vector<ll>{5, 1, 4, 1, 5}));
    SHOW(lz.all_prod(), "");
    RUN(lz.range_add(1, 3, 10));
    SHOW(lz.all_prod(), "[5, 11, 14, 1, 5]");
    RUN(lz.assign(0, 5, 7));
    SHOW(lz.all_prod(), "全部 7 → 2 番目は無い");
    RUN(lz.range_affine(3, 5, -1, 0));
    SHOW(lz.all_prod(), "[7, 7, 7, -7, -7]");
    SHOW(lz.get(4), "");
    note("MAX_ONLY / MIN_ONLY で負の数を掛けると assert で止まる（片側だけでは掛けたあとの最大・最小が分からない）");

    api("【8】 例：ABC343 F（1 点更新しながら、区間の 2 番目に大きい値の個数を答える）", "（使い方の例）",
        "同じ値は 1 つと数える 2 番目なので max2_count(l, r) そのもの（2 番目が無ければ 0）", "全体 O((N + Q) log N)");
    {
        ArgMinMax2SegTree<MAX_ONLY> st(std::vector<int>{3, 3, 1, 4, 5});  // 入力例 1
        std::vector<int> ans;
        ans.push_back(st.max2_count(0, 3));  // 2 1 3
        ans.push_back(st.max2_count(4, 5));  // 2 5 5
        st.set(2, 3);                        // 1 3 3
        ans.push_back(st.max2_count(1, 4));  // 2 2 4
        note("入力例 1 の答え → " + to_s(ans) + "（正解は [1, 0, 2]）");
    }

    api("【9】 例：ABC134 C（A_i を除いた最大値を i ごとに）", "（使い方の例）", "max_except(i) を i = 0〜N−1 で呼ぶだけ", "全体 O(N)");
    {
        ArgMinMax2SegTree<MAX_ONLY> st(std::vector<int>{1, 4, 3});  // 入力例 1
        std::vector<ll> ans;
        for (int i = 0; i < st.size(); ++i) ans.push_back(st.max_except(i));
        note("入力例 1 の答え → " + to_s(ans) + "（正解は [4, 3, 4]）");
        ArgMinMax2SegTree<MAX_ONLY> st2(std::vector<int>{5, 5});  // 入力例 2
        note("入力例 2 の答え → " + to_s(std::vector<ll>{st2.max_except(0), st2.max_except(1)}) + "（正解は [5, 5]）");
    }
}

}  // namespace amm2_demo

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
    amm2_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    amm2_test::benchmark();    // 速度計測
    amm2_demo::demo();         // 使い方デモ
    return 0;
}
