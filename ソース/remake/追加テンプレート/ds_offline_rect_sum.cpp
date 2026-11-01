// =============================================================================
//  長方形の和（座標が大きい・オフライン）：点に足す + 長方形の和、長方形に足す + 1 点の値   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  座標が 10^9 のように大きくて 2 次元 BIT の表を作れないとき、使う座標を先に全部読んでおいて（オフライン）、
//  BIT の各頂点に「そこに来る y の一覧」だけを持つ（圧縮した 2 次元 BIT）。足す・聞くを好きな順に混ぜてよい。
//  点が最初から全部あって長方形の和を聞くだけなら、x の順に掃く static_rectangle_sum の方が速い
//    例: 「点 (x, y) に w を足す」「[x1, x2) × [y1, y2) の和」がまざった Q 個の質問（座標 ≤ 10^9）
//
//  機能一覧（座標は long long。長方形は半開区間 [x1, x2) × [y1, y2)）                            計算量
//    OfflineBIT2D<T> bit(points)          足すかもしれない点 (x, y) の一覧を先に渡す                          O(N log^2 N)
//      bit.add(x, y, w)                   点 (x, y) に w を足す（points にある点だけ）                        O(log^2 N)
//      bit.sum(x1, y1, x2, y2)            [x1, x2) × [y1, y2) の点の重みの和（座標はどこでもよい）             O(log^2 N)
//      bit.prefix_sum(x, y)               x' < x かつ y' < y の和
//    RectAddPointGet<T> ra(rects)         足すかもしれない長方形 {x1, y1, x2, y2} の一覧を先に渡す            O(N log^2 N)
//      ra.add(x1, y1, x2, y2, w)          長方形に w を足す（rects にあるものだけ）                            O(log^2 N)
//      ra.get(x, y)                       点 (x, y) の値（覆っている長方形の w の和）                         O(log^2 N)
//    static_rectangle_sum(points, queries)   重み付きの点 {x, y, w} と長方形 {x1, y1, x2, y2} の一覧 → 長方形ごとの和  O((N + Q) log N)
//    point_add_rect_sum_offline(ops)      足す {0, x, y, w, 0}・聞く {1, x1, y1, x2, y2} の列を全部渡す → 聞いた順の答え
//                                         （時刻で分割統治。OfflineBIT2D より速い）                          O((N + Q) log^2 N)
//
//  注意: points にない点に add すると止める。重みの和は T（long long など）
//  速さの目安（N = Q = 2×10^5、座標がばらばら）: OfflineBIT2D 約 1.6 秒、point_add_rect_sum_offline 約 0.7 秒、static_rectangle_sum 約 0.3 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    rectsum_test::random_test()  … 点を全部見る愚直解との突き合わせ
//    rectsum_test::benchmark()    … 速度計測
//    rectsum_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace rectsum_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace rectsum_internal

template <class T = long long>
struct OfflineBIT2D {
    OfflineBIT2D() {}
    explicit OfflineBIT2D(const std::vector<std::pair<long long, long long>>& points) {
        for (const auto& p : points) xs_.push_back(p.first), yv_.push_back(p.second);
        std::sort(xs_.begin(), xs_.end());
        xs_.erase(std::unique(xs_.begin(), xs_.end()), xs_.end());
        std::sort(yv_.begin(), yv_.end());
        yv_.erase(std::unique(yv_.begin(), yv_.end()), yv_.end());
        const int n = (int)xs_.size();
        // BIT の頂点 i（1〜n）に来る y の番号を、1 本の配列 ys_[st_[i] .. st_[i+1]) にまとめる（4 バイトの int で持つ）
        std::vector<std::vector<int>> tmp(n + 1);
        for (const auto& p : points) {
            const int yr = yid(p.second);
            for (int i = xid(p.first) + 1; i <= n; i += i & -i) tmp[i].push_back(yr);
        }
        st_.assign(n + 2, 0);
        for (int i = 1; i <= n; ++i) {
            auto& v = tmp[i];
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
            st_[i + 1] = st_[i] + (int)v.size();
        }
        ys_.resize(st_[n + 1]);
        for (int i = 1; i <= n; ++i) std::copy(tmp[i].begin(), tmp[i].end(), ys_.begin() + st_[i]);
        bit_.assign(ys_.size() + n + 1, T(0));  // 頂点 i の BIT は bit_[st_[i] + i .. ]（1 始まりの分だけずらす）
    }
    void add(long long x, long long y, T w) {
        const int xi = xid(x), n = (int)xs_.size(), yr = yid(y);
        if (xi >= n || xs_[xi] != x || yr >= (int)yv_.size() || yv_[yr] != y) bad("OfflineBIT2D::add");
        for (int i = xi + 1; i <= n; i += i & -i) {
            const int* b = ys_.data() + st_[i];
            const int len = st_[i + 1] - st_[i];
            const int yi = (int)(std::lower_bound(b, b + len, yr) - b);
            if (yi >= len || b[yi] != yr) bad("OfflineBIT2D::add");
            T* t = bit_.data() + st_[i] + i;
            for (int j = yi + 1; j <= len; j += j & -j) t[j] += w;
        }
    }
    T prefix_sum(long long x, long long y) const { return col(xid(x), 0, yid(y)); }  // x' < x、y' < y
    T sum(long long x1, long long y1, long long x2, long long y2) const {
        if (x1 >= x2 || y1 >= y2) return T(0);
        const int a = yid(y1), b = yid(y2);
        return col(xid(x2), a, b) - col(xid(x1), a, b);
    }

   private:
    std::vector<long long> xs_, yv_;
    std::vector<int> ys_, st_;
    std::vector<T> bit_;
    int xid(long long x) const { return (int)(std::lower_bound(xs_.begin(), xs_.end(), x) - xs_.begin()); }  // x 未満の数
    int yid(long long y) const { return (int)(std::lower_bound(yv_.begin(), yv_.end(), y) - yv_.begin()); }  // y 未満の数
    [[noreturn]] static void bad(const char* f) { rectsum_internal::fail(f, "先に渡した点ではない"); }
    // x の番号が k 未満の点のうち、y の番号が [ya, yb) の和
    T col(int k, int ya, int yb) const {
        T s = T(0);
        for (int i = k; i > 0; i -= i & -i) {
            const int* b = ys_.data() + st_[i];
            const int len = st_[i + 1] - st_[i];
            const T* t = bit_.data() + st_[i] + i;
            const int hi = (int)(std::lower_bound(b, b + len, yb) - b);
            const int lo = ya == 0 ? 0 : (int)(std::lower_bound(b, b + hi, ya) - b);
            for (int j = hi; j > 0; j -= j & -j) s += t[j];
            for (int j = lo; j > 0; j -= j & -j) s -= t[j];
        }
        return s;
    }
};

template <class T = long long>
struct RectAddPointGet {
    explicit RectAddPointGet(const std::vector<std::array<long long, 4>>& rects) {
        std::vector<std::pair<long long, long long>> pts;
        for (const auto& r : rects) pts.push_back({r[0], r[1]}), pts.push_back({r[2], r[1]}), pts.push_back({r[0], r[3]}), pts.push_back({r[2], r[3]});
        bit_ = OfflineBIT2D<T>(pts);
    }
    void add(long long x1, long long y1, long long x2, long long y2, T w) {
        if (x1 >= x2 || y1 >= y2) return;
        bit_.add(x1, y1, w), bit_.add(x2, y1, T(0) - w), bit_.add(x1, y2, T(0) - w), bit_.add(x2, y2, w);
    }
    T get(long long x, long long y) const { return bit_.prefix_sum(x + 1, y + 1); }  // 角が (x, y) 以下の和

   private:
    OfflineBIT2D<T> bit_;
};

// points = {x, y, w}、queries = {x1, y1, x2, y2} → 長方形ごとの和
template <class T = long long>
std::vector<T> static_rectangle_sum(const std::vector<std::array<long long, 3>>& points, const std::vector<std::array<long long, 4>>& queries) {
    std::vector<long long> ys;
    for (const auto& p : points) ys.push_back(p[1]);
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    const int m = (int)ys.size();
    std::vector<T> bit(m + 1, T(0)), ans(queries.size(), T(0));
    // 質問を 4 つの角 F(X, Y) = Σ_{x < X, y < Y} に分けて、X の小さい順に掃く
    struct Ev {
        long long X, Y;
        int qi, sign;
    };
    std::vector<Ev> ev;
    for (size_t i = 0; i < queries.size(); ++i) {
        const auto& q = queries[i];
        if (q[0] >= q[2] || q[1] >= q[3]) continue;
        ev.push_back({q[2], q[3], (int)i, 1}), ev.push_back({q[0], q[3], (int)i, -1});
        ev.push_back({q[2], q[1], (int)i, -1}), ev.push_back({q[0], q[1], (int)i, 1});
    }
    std::sort(ev.begin(), ev.end(), [](const Ev& a, const Ev& b) { return a.X < b.X; });
    std::vector<int> ord(points.size());
    for (size_t i = 0; i < points.size(); ++i) ord[i] = (int)i;
    std::sort(ord.begin(), ord.end(), [&](int a, int b) { return points[a][0] < points[b][0]; });
    size_t pi = 0;
    for (const Ev& e : ev) {
        while (pi < ord.size() && points[ord[pi]][0] < e.X) {
            const auto& p = points[ord[pi++]];
            for (int j = (int)(std::lower_bound(ys.begin(), ys.end(), p[1]) - ys.begin()) + 1; j <= m; j += j & -j) bit[j] += T(p[2]);
        }
        T s = T(0);
        for (int j = (int)(std::lower_bound(ys.begin(), ys.end(), e.Y) - ys.begin()); j > 0; j -= j & -j) s += bit[j];
        ans[e.qi] += e.sign > 0 ? s : T(0) - s;
    }
    return ans;
}
// 足す・聞くの列を先に全部読んで答える（時刻で分割統治 = CDQ）。ops[i] = {0, x, y, w, 0}（点に足す）か
// {1, x1, y1, x2, y2}（[x1, x2) × [y1, y2) の和を聞く）。答えは「聞く」の順に
template <class T = long long>
std::vector<T> point_add_rect_sum_offline(const std::vector<std::array<long long, 5>>& ops) {
    std::vector<long long> ys;
    for (const auto& o : ops)
        if (o[0] == 0) ys.push_back(o[2]);
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    const int m = (int)ys.size();
    auto yid = [&](long long y) { return (int)(std::lower_bound(ys.begin(), ys.end(), y) - ys.begin()); };
    struct Ev {
        long long x;
        int kind;  // 0 = 聞く角（x' < X の点を数える）、1 = 足す
        int y;     // 足す：y の番号、角：y の番号の上限（それ未満を数える）
        int q, sg;
        T w;
    };
    std::vector<Ev> E;  // 時刻の順
    int nq = 0;
    for (const auto& o : ops) {
        if (o[0] == 0) {
            E.push_back({o[1], 1, yid(o[2]), -1, 0, T(o[3])});
        } else {
            if (o[1] < o[3] && o[2] < o[4]) {
                const int ya = yid(o[2]), yb = yid(o[4]);
                E.push_back({o[3], 0, yb, nq, 1, T(0)}), E.push_back({o[3], 0, ya, nq, -1, T(0)});
                E.push_back({o[1], 0, yb, nq, -1, T(0)}), E.push_back({o[1], 0, ya, nq, 1, T(0)});
            }
            ++nq;
        }
    }
    std::vector<T> ans(nq, T(0)), bit(m + 1, T(0));
    std::vector<Ev> buf(E.size());
    auto key_less = [](const Ev& p, const Ev& q) { return p.x != q.x ? p.x < q.x : p.kind < q.kind; };
    // 時刻で半分に分け、左の「足す」を右の「角」に効かせる。各段で x の順にマージする（マージソートと同じ形）
    auto rec = [&](auto&& self, int l, int r) -> void {
        if (r - l <= 1) return;
        const int mid = (l + r) / 2;
        self(self, l, mid), self(self, mid, r);
        int i = l;
        for (int j = mid; j < r; ++j) {
            if (E[j].kind != 0) continue;
            while (i < mid && key_less(E[i], E[j])) {
                if (E[i].kind == 1)
                    for (int k = E[i].y + 1; k <= m; k += k & -k) bit[k] += E[i].w;
                ++i;
            }
            T sm = T(0);
            for (int k = E[j].y; k > 0; k -= k & -k) sm += bit[k];
            ans[E[j].q] += E[j].sg > 0 ? sm : T(0) - sm;
        }
        for (int t = l; t < i; ++t)
            if (E[t].kind == 1)
                for (int k = E[t].y + 1; k <= m; k += k & -k) bit[k] = T(0);
        std::merge(E.begin() + l, E.begin() + mid, E.begin() + mid, E.begin() + r, buf.begin() + l, key_less);
        std::copy(buf.begin() + l, buf.begin() + r, E.begin() + l);
    };
    rec(rec, 0, (int)E.size());
    return ans;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace rectsum_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define RS_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_offline_rect_sum] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                      \
        }                                                                                      \
    } while (0)

// [lo, hi] の一様乱数
[[maybe_unused]] static long long rnd(std::mt19937_64& rng, long long lo, long long hi) {
    const unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return (long long)((unsigned long long)lo + (w == 0 ? rng() : rng() % w));
}
[[maybe_unused]] static long long elapsed_ms(std::chrono::steady_clock::time_point a) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
}

using ll = long long;

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const ll C = rnd(rng, 0, 1) ? 8 : 1000000000;
        const int n = (int)rnd(rng, 0, 30);
        std::vector<std::pair<ll, ll>> pts(n);
        for (auto& p : pts) p = {rnd(rng, -C, C), rnd(rng, -C, C)};
        OfflineBIT2D<ll> bit(pts);
        std::vector<std::array<ll, 3>> added;  // 足した点（愚直）
        for (int t = 0; t < 60; ++t) {
            if (n > 0 && rnd(rng, 0, 1)) {
                const auto p = pts[rng() % n];
                const ll w = rnd(rng, -5, 5);
                bit.add(p.first, p.second, w);
                added.push_back({p.first, p.second, w});
            } else {
                ll x1 = rnd(rng, -C - 1, C + 1), x2 = rnd(rng, -C - 1, C + 1), y1 = rnd(rng, -C - 1, C + 1), y2 = rnd(rng, -C - 1, C + 1);
                if (n > 0 && rnd(rng, 0, 1)) x1 = pts[rng() % n].first, y2 = pts[rng() % n].second;  // ちょうど点の座標も
                ll s = 0;
                for (const auto& a : added)
                    if (x1 <= a[0] && a[0] < x2 && y1 <= a[1] && a[1] < y2) s += a[2];
                RS_CHECK(bit.sum(x1, y1, x2, y2) == s);
                ++cases;
            }
        }
        // point_add_rect_sum_offline：同じ操作の列で
        {
            std::vector<std::array<ll, 5>> ops;
            std::vector<ll> exp;
            std::vector<std::array<ll, 3>> sofar;
            for (int t = 0; t < 40; ++t) {
                if (rnd(rng, 0, 1)) {
                    const std::array<ll, 3> a{rnd(rng, -C, C), rnd(rng, -C, C), rnd(rng, -5, 5)};
                    ops.push_back({0, a[0], a[1], a[2], 0}), sofar.push_back(a);
                } else {
                    const ll x1 = rnd(rng, -C - 1, C + 1), x2 = rnd(rng, -C - 1, C + 1), y1 = rnd(rng, -C - 1, C + 1), y2 = rnd(rng, -C - 1, C + 1);
                    ops.push_back({1, x1, y1, x2, y2});
                    ll sm = 0;
                    for (const auto& a : sofar)
                        if (x1 <= a[0] && a[0] < x2 && y1 <= a[1] && a[1] < y2) sm += a[2];
                    exp.push_back(sm);
                }
            }
            RS_CHECK(point_add_rect_sum_offline<ll>(ops) == exp);
            ++cases;
        }
        // static_rectangle_sum
        std::vector<std::array<ll, 4>> qs(rnd(rng, 0, 20));
        for (auto& q : qs) q = {rnd(rng, -C - 1, C + 1), rnd(rng, -C - 1, C + 1), rnd(rng, -C - 1, C + 1), rnd(rng, -C - 1, C + 1)};
        const auto res = static_rectangle_sum<ll>(added, qs);
        for (size_t i = 0; i < qs.size(); ++i) {
            ll s = 0;
            for (const auto& a : added)
                if (qs[i][0] <= a[0] && a[0] < qs[i][2] && qs[i][1] <= a[1] && a[1] < qs[i][3]) s += a[2];
            RS_CHECK(res[i] == s);
            ++cases;
        }
        // RectAddPointGet
        std::vector<std::array<ll, 4>> rects(rnd(rng, 0, 10));
        for (auto& r : rects) {
            r = {rnd(rng, -C, C), rnd(rng, -C, C), rnd(rng, -C, C), rnd(rng, -C, C)};
            if (r[0] > r[2]) std::swap(r[0], r[2]);
            if (r[1] > r[3]) std::swap(r[1], r[3]);
        }
        RectAddPointGet<ll> ra(rects);
        std::vector<std::pair<std::array<ll, 4>, ll>> done;
        for (int t = 0; t < 30; ++t) {
            if (!rects.empty() && rnd(rng, 0, 1)) {
                const auto r = rects[rng() % rects.size()];
                const ll w = rnd(rng, -5, 5);
                ra.add(r[0], r[1], r[2], r[3], w);
                done.push_back({r, w});
            } else {
                const ll x = rnd(rng, -C, C), y = rnd(rng, -C, C);
                ll s = 0;
                for (const auto& [r, w] : done)
                    if (r[0] <= x && x < r[2] && r[1] <= y && y < r[3]) s += w;
                RS_CHECK(ra.get(x, y) == s);
                ++cases;
            }
        }
    }
    std::cout << "[ds_offline_rect_sum] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_offline_rect_sum] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    const ll C = 1000000000;
    std::vector<std::pair<ll, ll>> pts(N);
    for (auto& p : pts) p = {rnd(rng, 0, C), rnd(rng, 0, C)};
    auto t0 = std::chrono::steady_clock::now();
    OfflineBIT2D<ll> bit(pts);
    ll chk = 0;
    for (int i = 0; i < N; ++i) {
        bit.add(pts[i].first, pts[i].second, i);
        ll a = rnd(rng, 0, C), b = rnd(rng, 0, C), c = rnd(rng, 0, C), d = rnd(rng, 0, C);
        chk += bit.sum(std::min(a, b), std::min(c, d), std::max(a, b), std::max(c, d));
    }
    const ll t1 = elapsed_ms(t0);
    std::vector<std::array<ll, 3>> p3(N);
    std::vector<std::array<ll, 4>> q4(N);
    for (int i = 0; i < N; ++i) {
        p3[i] = {pts[i].first, pts[i].second, i};
        ll a = rnd(rng, 0, C), b = rnd(rng, 0, C), c = rnd(rng, 0, C), d = rnd(rng, 0, C);
        q4[i] = {std::min(a, b), std::min(c, d), std::max(a, b), std::max(c, d)};
    }
    t0 = std::chrono::steady_clock::now();
    for (ll v : static_rectangle_sum<ll>(p3, q4)) chk += v;
    const ll t2 = elapsed_ms(t0);
    std::vector<std::array<ll, 5>> ops;
    for (int i = 0; i < N; ++i) ops.push_back({0, p3[i][0], p3[i][1], p3[i][2], 0}), ops.push_back({1, q4[i][0], q4[i][1], q4[i][2], q4[i][3]});
    t0 = std::chrono::steady_clock::now();
    for (ll v : point_add_rect_sum_offline<ll>(ops)) chk += v;
    const ll t3 = elapsed_ms(t0);
    std::cout << "[ds_offline_rect_sum] benchmark (N = Q = " << N << "): OfflineBIT2D（構築 + add + sum） " << t1 << " ms / static_rectangle_sum " << t2
              << " ms / point_add_rect_sum_offline " << t3 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace rectsum_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace rectsum_demo {

[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(const std::string& s) { return "\"" + s + "\""; }
[[maybe_unused]] static std::string to_s(const char* s) { return to_s(std::string(s)); }
[[maybe_unused]] static std::string to_s(char c) { return std::string("'") + c + "'"; }
[[maybe_unused]] static std::string to_s(double v) {
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}
template <class T>
static std::string to_s(const T& v);
template <class A, class B>
static std::string to_s(const std::pair<A, B>& p);
template <class T>
static std::string to_s(const std::vector<T>& v);
template <class T>
static std::string to_s(const std::optional<T>& v);
template <class T>
static std::string to_s(const T& v) {  // << で出せるもの（整数など）
    std::ostringstream os;
    os << v;
    return os.str();
}
template <class A, class B>
static std::string to_s(const std::pair<A, B>& p) {
    return "(" + to_s(p.first) + ", " + to_s(p.second) + ")";
}
template <class T>
static std::string to_s(const std::vector<T>& v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "]";
}
template <class T>
static std::string to_s(const std::optional<T>& v) {
    return v ? to_s(*v) : std::string("nullopt");
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
#define SHOW(expr, comment) rectsum_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    rectsum_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("長方形の和（ds_offline_rect_sum） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  OfflineBIT2D<T> bit(points): add(x, y, w) / sum(x1, y1, x2, y2) / prefix_sum(x, y)\n"
              << "  RectAddPointGet<T> ra(rects): add(x1, y1, x2, y2, w) / get(x, y)\n"
              << "  static_rectangle_sum(points, queries)\n";

    api("【1】 OfflineBIT2D<> bit(points) / add(x, y, w) / sum(x1, y1, x2, y2)", "— / void / T",
        "足す点を先に全部渡す（座標は大きくてよい）。足すのと聞くのは好きな順に", "O(N log^2 N) / O(log^2 N) / O(log^2 N)");
    RUN(OfflineBIT2D<> bit({{1, 1}, {1000000000, 5}, {3, 1000000000}}));
    RUN(bit.add(1, 1, 10));
    SHOW(bit.sum(0, 0, 2, 2), "(1, 1) だけ");
    RUN(bit.add(1000000000, 5, 7); bit.add(3, 1000000000, 100));
    SHOW(bit.sum(0, 0, 1000000001, 1000000001), "全部");
    SHOW(bit.sum(2, 0, 1000000001, 6), "x ≥ 2 かつ y < 6");

    api("【2】 RectAddPointGet<> ra(rects) / add(x1, y1, x2, y2, w) / get(x, y)", "— / void / T", "長方形に足して、1 点の値を聞く（足す長方形は先に全部渡す）",
        "O(N log^2 N) / O(log^2 N) / O(log^2 N)");
    RUN(RectAddPointGet<> ra({{0, 0, 10, 10}, {5, 5, 20, 20}}));
    RUN(ra.add(0, 0, 10, 10, 1); ra.add(5, 5, 20, 20, 2));
    SHOW(ra.get(7, 7), "両方に入る");
    SHOW(ra.get(15, 15), "2 つ目だけ");
    SHOW(ra.get(10, 3), "x = 10 は 1 つ目に入らない（半開区間）");

    api("【3】 static_rectangle_sum(points, queries)", "vector<T>", "点が最初から全部あるとき、長方形ごとの和（x の順に掃く）", "O((N + Q) log N)");
    SHOW(static_rectangle_sum<long long>({{1, 1, 5}, {2, 3, 7}, {4, 4, 1}}, {{0, 0, 3, 4}, {2, 2, 5, 5}}), "");

    api("【4】 point_add_rect_sum_offline(ops)", "vector<T>", "足す {0, x, y, w, 0}・聞く {1, x1, y1, x2, y2} の列を全部渡すと、聞いた順の答え（OfflineBIT2D より速い）",
        "O((N + Q) log^2 N)");
    SHOW(point_add_rect_sum_offline<long long>({{0, 1, 1, 5, 0}, {1, 0, 0, 3, 3}, {0, 2, 2, 7, 0}, {1, 0, 0, 3, 3}}), "2 回目は 5 + 7");
}

}  // namespace rectsum_demo

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
    rectsum_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    rectsum_test::benchmark();    // 速度計測
    rectsum_demo::demo();         // 使い方デモ
    return 0;
}
