// =============================================================================
//  走査線（平面走査）の型・縦線と横線の交差数・最も多く重なる点・長方形の和集合の周の長さ   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  縦の線（x = 一定）を左から右へ動かし、「始まる・終わる・聞く」の出来事（イベント）を x の順に処理する。
//  同じ x のイベントの順番（追加が先か削除が先か）で、端点で触れるものを数えるかが決まるので、それを order で指定する型にした
//    例: 「横線と縦線が何か所で交わるか」「長方形が最も多く重なる点」「長方形をいくつか置いたときの外周の長さ」
//
//  機能一覧（座標は long long。|座標| ≤ 10^9 くらいまで）                                        計算量
//    SweepLine<X, D> sw                  イベントを x の順に処理する型（D はイベントに付ける情報）
//    sw.add(x, order, data)              イベントを足す。同じ x なら order の小さい順、それも同じなら足した順    O(1)
//    sw.run(f)                           f(const Event& e)（e.x, e.order, e.data）を順に呼ぶ                     O(K log K)
//    sw.run_grouped(f)                   同じ x のイベントをまとめて f(x, b, e)（[b, e) がそのイベント。b->data など）
//    sw.events()                         並べたあとのイベントの一覧
//    count_hv_intersections(H, V)        横線 H {y, x1, x2} と縦線 V {x, y1, y2} の交わる組の数（端で触れるのも数える） O((H + V) log)
//    hv_intersection_counts(H, V)        {per_h, per_v}：それぞれの線が何本と交わるか                             O((H + V) log)
//    max_overlap_intervals(iv)           半開区間 [l, r) が最も多く重なる点 {count, x}（x は一番左。count = 0 なら x = 0） O(N log N)
//    max_overlap_rects(rs)               半開の長方形 [x1, x2) × [y1, y2) が最も多く重なる点 {count, x, y}            O(N log N)
//    union_perimeter(rs)                 長方形の和集合の周の長さ（穴の周も含む。面積は geometry_v2/08_rect.cpp の union_area）  O(N log N)
//
//  注意: 長方形はすべて Rect {x1, y1, x2, y2}（x1 < x2, y1 < y2。空のものは無視）。格子点の閉区間 [L, R] は x2 = R + 1 で半開にする
//        横線・縦線は端を含む（x1 > x2 でも入れ替えて扱う）。長さ 0 の線（点）もよい
//  速さの目安: count_hv_intersections（各 2×10^5 本）約 0.2 秒 / max_overlap_rects（2×10^5 個）約 0.35 秒 / union_perimeter（2×10^5 個）約 0.9 秒
//              （走査を縦横 2 回。同じ計測環境で geometry_v2/08_rect.cpp の union_area 2×10^5 個が 約 0.3 秒）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    sweep_test::random_test()  … 全部の組・全部のマスを調べる愚直解との突き合わせ
//    sweep_test::benchmark()    … 速度計測
//    sweep_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
template <class X, class D = int>
struct SweepLine {
    struct Event {
        X x;
        int order;
        D data;
    };
    void add(const X& x, int order, const D& data) {
        ev_.push_back({x, order, data});
        sorted_ = false;
    }
    size_t size() const { return ev_.size(); }
    template <class F>
    void run(F f) {
        prepare();
        for (const Event& e : ev_) f(e);
    }
    template <class F>
    void run_grouped(F f) {
        prepare();
        for (size_t i = 0; i < ev_.size();) {
            size_t j = i;
            while (j < ev_.size() && !(ev_[i].x < ev_[j].x)) ++j;
            f(ev_[i].x, ev_.cbegin() + i, ev_.cbegin() + j);
            i = j;
        }
    }
    const std::vector<Event>& events() {
        prepare();
        return ev_;
    }

   private:
    std::vector<Event> ev_;
    bool sorted_ = true;
    void prepare() {
        if (sorted_) return;
        std::stable_sort(ev_.begin(), ev_.end(), [](const Event& a, const Event& b) { return a.x < b.x || (!(b.x < a.x) && a.order < b.order); });
        sorted_ = true;
    }
};

struct HSeg {
    long long y, x1, x2;
};
struct VSeg {
    long long x, y1, y2;
};
struct Rect {
    long long x1, y1, x2, y2;
};
struct HVCounts {
    std::vector<long long> per_h, per_v;
};
struct Overlap1D {
    int count;
    long long x;
};
struct Overlap2D {
    int count;
    long long x, y;
};

namespace sweep_internal {
struct BIT {
    std::vector<long long> t;
    explicit BIT(int n) : t(n + 1, 0) {}
    void add(int i, long long v) {
        for (++i; i < (int)t.size(); i += i & -i) t[i] += v;
    }
    long long sum(int r) const {  // [0, r)
        long long s = 0;
        for (; r > 0; r -= r & -r) s += t[r];
        return s;
    }
};
// 横線を BIT に置き、縦線ごとに [y1, y2] の本数を数える
inline std::vector<long long> per_vertical(const std::vector<HSeg>& H, const std::vector<VSeg>& V) {
    std::vector<long long> ys;
    ys.reserve(H.size());
    for (const auto& h : H) ys.push_back(h.y);
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    auto id = [&](long long y) { return (int)(std::lower_bound(ys.begin(), ys.end(), y) - ys.begin()); };
    SweepLine<long long, int> sw;  // order：0 横線が始まる / 1 縦線 / 2 横線が終わる（端で触れるのも数える）
    for (int i = 0; i < (int)H.size(); ++i) {
        sw.add(std::min(H[i].x1, H[i].x2), 0, i);
        sw.add(std::max(H[i].x1, H[i].x2), 2, i);
    }
    for (int i = 0; i < (int)V.size(); ++i) sw.add(V[i].x, 1, i);
    std::vector<int> hid(H.size());
    for (size_t i = 0; i < H.size(); ++i) hid[i] = id(H[i].y);
    std::vector<std::pair<int, int>> vid(V.size());  // 縦線の y の範囲 → BIT の [l, r)
    for (size_t i = 0; i < V.size(); ++i)
        vid[i] = {id(std::min(V[i].y1, V[i].y2)), (int)(std::upper_bound(ys.begin(), ys.end(), std::max(V[i].y1, V[i].y2)) - ys.begin())};
    BIT bit((int)ys.size());
    std::vector<long long> res(V.size(), 0);
    sw.run([&](const auto& e) {
        if (e.order == 1) {
            const auto [l, r] = vid[e.data];
            res[e.data] = l < r ? bit.sum(r) - bit.sum(l) : 0;
        } else {
            bit.add(hid[e.data], e.order == 0 ? 1 : -1);
        }
    });
    return res;
}
// 区間に足す・全体の最大とその位置（遅延を下ろさない形）。節点は「左の子 = k + 1、右の子 = k + 2 × (左の長さ)」に並べて 2n − 1 個
// （部分木が連続に並ぶので、深いところの読み書きが近くにまとまって速い）
struct AddMax {
    struct Node {
        int mx = 0, lz = 0;
    };
    int n;
    std::vector<Node> t;
    explicit AddMax(int n_) : n(n_), t(std::max(1, 2 * n_ - 1)) {}
    void add(int l, int r, int v) {
        if (l < r) add(0, 0, n, l, r, v);
    }
    void add(int k, int nl, int nr, int l, int r, int v) {
        if (l <= nl && nr <= r) {
            t[k].mx += v, t[k].lz += v;
            return;
        }
        const int m = (nl + nr) / 2, a = k + 1, b = k + 2 * (m - nl);
        if (l < m) add(a, nl, m, l, r, v);
        if (m < r) add(b, m, nr, l, r, v);
        t[k].mx = std::max(t[a].mx, t[b].mx) + t[k].lz;
    }
    int top() const { return t[0].mx; }
    int argmax() const {
        int k = 0, nl = 0, nr = n;
        while (nr - nl > 1) {
            const int m = (nl + nr) / 2, a = k + 1, b = k + 2 * (m - nl);
            if (t[a].mx >= t[b].mx) k = a, nr = m;
            else k = b, nl = m;
        }
        return nl;
    }
};
// 覆われている長さ（並べ方は AddMax と同じ）
struct CoverLen {
    struct Node {
        long long width = 0, len = 0;
        int cnt = 0;
    };
    int n;
    std::vector<Node> t;
    explicit CoverLen(const std::vector<long long>& ys) : n((int)ys.size() - 1), t(std::max(1, 2 * ((int)ys.size() - 1) - 1)) {
        if (n > 0) build(0, 0, n, ys);
    }
    void build(int k, int nl, int nr, const std::vector<long long>& ys) {
        t[k].width = ys[nr] - ys[nl];
        if (nr - nl == 1) return;
        const int m = (nl + nr) / 2;
        build(k + 1, nl, m, ys), build(k + 2 * (m - nl), m, nr, ys);
    }
    void add(int l, int r, int v) {
        if (l < r) add(0, 0, n, l, r, v);
    }
    void add(int k, int nl, int nr, int l, int r, int v) {
        const int m = (nl + nr) / 2, a = k + 1, b = k + 2 * (m - nl);
        if (l <= nl && nr <= r) {
            t[k].cnt += v;
        } else {
            if (l < m) add(a, nl, m, l, r, v);
            if (m < r) add(b, m, nr, l, r, v);
        }
        if (t[k].cnt > 0) t[k].len = t[k].width;
        else if (nr - nl == 1) t[k].len = 0;
        else t[k].len = t[a].len + t[b].len;
    }
    long long covered() const { return n > 0 ? t[0].len : 0; }
};
// x 方向に走査して、縦の辺の長さの和（|覆われた長さの変化| の和。同じ x では追加が先）
inline long long vertical_edges(const std::vector<Rect>& rs) {
    std::vector<long long> ys;
    for (const auto& r : rs) ys.push_back(r.y1), ys.push_back(r.y2);
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    if (ys.size() < 2) return 0;
    auto id = [&](long long y) { return (int)(std::lower_bound(ys.begin(), ys.end(), y) - ys.begin()); };
    std::vector<std::pair<int, int>> yid(rs.size());
    for (size_t i = 0; i < rs.size(); ++i) yid[i] = {id(rs[i].y1), id(rs[i].y2)};
    SweepLine<long long, int> sw;
    for (int i = 0; i < (int)rs.size(); ++i) sw.add(rs[i].x1, 0, i), sw.add(rs[i].x2, 1, i);
    CoverLen cl(ys);
    long long total = 0;
    sw.run([&](const auto& e) {
        const long long before = cl.covered();
        cl.add(yid[e.data].first, yid[e.data].second, e.order == 0 ? 1 : -1);
        const long long after = cl.covered();
        total += after > before ? after - before : before - after;
    });
    return total;
}
inline std::vector<Rect> valid_rects(const std::vector<Rect>& rs) {
    std::vector<Rect> v;
    for (const auto& r : rs)
        if (r.x1 < r.x2 && r.y1 < r.y2) v.push_back(r);
    return v;
}
}  // namespace sweep_internal

inline HVCounts hv_intersection_counts(const std::vector<HSeg>& H, const std::vector<VSeg>& V) {
    HVCounts res;
    res.per_v = sweep_internal::per_vertical(H, V);
    // 縦横を入れ替えて、横線ごとの数
    std::vector<HSeg> h2;
    std::vector<VSeg> v2;
    for (const auto& v : V) h2.push_back({v.x, v.y1, v.y2});
    for (const auto& h : H) v2.push_back({h.y, h.x1, h.x2});
    res.per_h = sweep_internal::per_vertical(h2, v2);
    return res;
}
inline long long count_hv_intersections(const std::vector<HSeg>& H, const std::vector<VSeg>& V) {
    long long s = 0;
    for (long long c : sweep_internal::per_vertical(H, V)) s += c;
    return s;
}
inline Overlap1D max_overlap_intervals(const std::vector<std::pair<long long, long long>>& iv) {
    SweepLine<long long, int> sw;  // 同じ x では終わり（−1）が先：[l, r) なので r では重ならない
    for (const auto& [l, r] : iv)
        if (l < r) sw.add(l, 1, +1), sw.add(r, 0, -1);
    Overlap1D best{0, 0};
    int cur = 0;
    sw.run_grouped([&](long long x, auto b, auto e) {
        for (auto it = b; it != e; ++it) cur += it->data;
        if (cur > best.count) best = {cur, x};
    });
    return best;
}
inline Overlap2D max_overlap_rects(const std::vector<Rect>& rects) {
    const auto rs = sweep_internal::valid_rects(rects);
    Overlap2D best{0, 0, 0};
    if (rs.empty()) return best;
    std::vector<long long> ys;
    for (const auto& r : rs) ys.push_back(r.y1), ys.push_back(r.y2);
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
    auto id = [&](long long y) { return (int)(std::lower_bound(ys.begin(), ys.end(), y) - ys.begin()); };
    SweepLine<long long, int> sw;  // 同じ x では削除（x2 = x）が先、追加（x1 = x）が後。まとめて処理してから最大を見る
    for (int i = 0; i < (int)rs.size(); ++i) sw.add(rs[i].x1, 1, i), sw.add(rs[i].x2, 0, i);
    std::vector<std::pair<int, int>> yid(rs.size());
    for (size_t i = 0; i < rs.size(); ++i) yid[i] = {id(rs[i].y1), id(rs[i].y2)};
    sweep_internal::AddMax seg((int)ys.size() - 1);
    sw.run_grouped([&](long long x, auto b, auto e) {
        for (auto it = b; it != e; ++it) seg.add(yid[it->data].first, yid[it->data].second, it->order == 1 ? 1 : -1);
        if (seg.top() > best.count) best = {seg.top(), x, ys[seg.argmax()]};
    });
    return best;
}
inline long long union_perimeter(const std::vector<Rect>& rects) {
    const auto rs = sweep_internal::valid_rects(rects);
    std::vector<Rect> tr;
    for (const auto& r : rs) tr.push_back({r.y1, r.x1, r.y2, r.x2});  // x と y を入れ替え
    return sweep_internal::vertical_edges(rs) + sweep_internal::vertical_edges(tr);
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace sweep_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define SW_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[geo_sweep_line] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const ll C = rnd(rng, 1, 12);
        // 横線・縦線
        {
            const int nh = (int)rnd(rng, 0, 10), nv = (int)rnd(rng, 0, 10);
            std::vector<HSeg> H(nh);
            std::vector<VSeg> V(nv);
            for (auto& h : H) h = {rnd(rng, -C, C), rnd(rng, -C, C), rnd(rng, -C, C)};
            for (auto& v : V) v = {rnd(rng, -C, C), rnd(rng, -C, C), rnd(rng, -C, C)};
            std::vector<ll> ph(nh, 0), pv(nv, 0);
            ll tot = 0;
            for (int i = 0; i < nh; ++i)
                for (int j = 0; j < nv; ++j) {
                    const bool hit = std::min(H[i].x1, H[i].x2) <= V[j].x && V[j].x <= std::max(H[i].x1, H[i].x2) && std::min(V[j].y1, V[j].y2) <= H[i].y &&
                                     H[i].y <= std::max(V[j].y1, V[j].y2);
                    if (hit) ++ph[i], ++pv[j], ++tot;
                }
            SW_CHECK(count_hv_intersections(H, V) == tot);
            const auto hc = hv_intersection_counts(H, V);
            SW_CHECK(hc.per_h == ph && hc.per_v == pv);
            cases += 2;
        }
        // 区間
        {
            const int n = (int)rnd(rng, 0, 8);
            std::vector<std::pair<ll, ll>> iv(n);
            for (auto& [l, r] : iv) l = rnd(rng, -C, C), r = rnd(rng, -C, C + 2);
            int bc = 0;
            ll bx = 0;
            for (ll x = -C; x <= C + 2; ++x) {
                int c = 0;
                for (auto [l, r] : iv) c += l <= x && x < r;
                if (c > bc) bc = c, bx = x;
            }
            const auto o = max_overlap_intervals(iv);
            SW_CHECK(o.count == bc && o.x == bx);
            ++cases;
        }
        // 長方形：重なりの最大と周の長さ（0 ≤ 座標 ≤ C のマスで愚直）
        {
            const int n = (int)rnd(rng, 0, 7);
            std::vector<Rect> rs(n);
            for (auto& r : rs) r = {rnd(rng, 0, C), rnd(rng, 0, C), rnd(rng, 0, C), rnd(rng, 0, C)};
            std::vector<std::vector<int>> cov(C + 2, std::vector<int>(C + 2, 0));
            int bc = 0;
            for (ll x = 0; x <= C; ++x)
                for (ll y = 0; y <= C; ++y) {
                    for (const auto& r : rs) cov[x][y] += r.x1 <= x && x < r.x2 && r.y1 <= y && y < r.y2;
                    bc = std::max(bc, cov[x][y]);
                }
            const auto o = max_overlap_rects(rs);
            SW_CHECK(o.count == bc);
            if (bc > 0) {
                int c = 0;
                for (const auto& r : rs) c += r.x1 <= o.x && o.x < r.x2 && r.y1 <= o.y && o.y < r.y2;
                SW_CHECK(c == bc);
            }
            // マス (x, y) = [x, x + 1) × [y, y + 1) が覆われているか。隣と違う辺を数える
            auto in = [&](ll x, ll y) { return x >= 0 && y >= 0 && x <= C && y <= C && cov[x][y] > 0; };
            ll per = 0;
            for (ll x = -1; x <= C; ++x)
                for (ll y = -1; y <= C; ++y) per += (in(x, y) != in(x + 1, y)) + (in(x, y) != in(x, y + 1));
            SW_CHECK(union_perimeter(rs) == per);
            cases += 2;
        }
        // SweepLine の順番：x、order、足した順
        {
            SweepLine<int, int> sw;
            const int k = (int)rnd(rng, 0, 12);
            std::vector<std::array<int, 3>> ev(k);
            for (int i = 0; i < k; ++i) ev[i] = {(int)rnd(rng, 0, 3), (int)rnd(rng, 0, 2), i}, sw.add(ev[i][0], ev[i][1], i);
            std::sort(ev.begin(), ev.end());
            std::vector<int> got;
            sw.run([&](const auto& e) { got.push_back(e.data); });
            std::vector<int> want;
            for (auto& e : ev) want.push_back(e[2]);
            SW_CHECK(got == want);
            int groups = 0, total = 0;
            sw.run_grouped([&](int x, auto b, auto e) {
                ++groups;
                for (auto p = b; p != e; ++p) {
                    SW_CHECK(p->x == x);
                    ++total;
                }
            });
            int distinct = 0;
            for (int i = 0; i < k; ++i) distinct += i == 0 || ev[i][0] != ev[i - 1][0];
            SW_CHECK(groups == distinct && total == k);
            cases += 2;
        }
    }
    std::cout << "[geo_sweep_line] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[geo_sweep_line] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 3000 : 200000;  // デバッグビルドは lower_bound の検査（毎回 O(N)）が重いので小さく
    const ll C = 1000000000;
    {
        std::vector<HSeg> H(N);
        std::vector<VSeg> V(N);
        for (auto& h : H) {
            const ll x = rnd(rng, 0, C);
            h = {rnd(rng, 0, C), x, x + rnd(rng, 0, C / 100)};
        }
        for (auto& v : V) {
            const ll y = rnd(rng, 0, C);
            v = {rnd(rng, 0, C), y, y + rnd(rng, 0, C / 100)};
        }
        const auto t0 = std::chrono::steady_clock::now();
        const ll c = count_hv_intersections(H, V);
        std::cout << "[geo_sweep_line] benchmark: count_hv_intersections（各 " << N << " 本） " << elapsed_ms(t0) << " ms  (交点 " << c << ")\n";
    }
    std::vector<Rect> rs(N);
    for (auto& r : rs) {
        const ll x = rnd(rng, 0, C), y = rnd(rng, 0, C);
        r = {x, y, x + rnd(rng, 1, C / 50), y + rnd(rng, 1, C / 50)};
    }
    {
        const auto t0 = std::chrono::steady_clock::now();
        const auto o = max_overlap_rects(rs);
        std::cout << "[geo_sweep_line] benchmark: max_overlap_rects（" << N << " 個） " << elapsed_ms(t0) << " ms  (最大 " << o.count << ")\n";
    }
    {
        const auto t0 = std::chrono::steady_clock::now();
        const ll p = union_perimeter(rs);
        std::cout << "[geo_sweep_line] benchmark: union_perimeter（" << N << " 個） " << elapsed_ms(t0) << " ms  (check " << p % 1000 << ")\n";
    }
}

}  // namespace sweep_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace sweep_demo {

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
#define SHOW(expr, comment) sweep_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    sweep_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("走査線（geo_sweep_line） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  SweepLine<X, D>: add(x, order, data) / run(f) / run_grouped(f) / events()\n"
              << "  count_hv_intersections(H, V) / hv_intersection_counts(H, V) / max_overlap_intervals(iv) / max_overlap_rects(rs) / union_perimeter(rs)\n";

    api("【1】 SweepLine<X, D>：add(x, order, data) → run(f) / run_grouped(f)", "void", "x の順（同じ x なら order の小さい順）にイベントを処理する", "O(K log K)");
    code("// 人の出入り：同じ時刻なら「出る（order 0）」を先にすると、ちょうど入れ替わる人は重ならない");
    RUN(SweepLine<int, int> sw);
    RUN(sw.add(1, 1, +1); sw.add(4, 0, -1); sw.add(4, 1, +1); sw.add(6, 0, -1); sw.add(2, 1, +1); sw.add(5, 0, -1));
    RUN(int cur = 0, best = 0);
    RUN(sw.run([&](const auto& e) { cur += e.data, best = std::max(best, cur); }));
    SHOW(best, "[1, 4)・[4, 6)・[2, 5) が同時にいる最大の人数");
    RUN(std::vector<std::string> log);
    RUN(sw.run_grouped([&](int x, auto b, auto e) { log.push_back(std::to_string(x) + ":" + std::to_string(e - b) + "件"); }));
    SHOW(log, "同じ x のイベントをまとめて");

    api("【2】 count_hv_intersections(H, V) / hv_intersection_counts(H, V)", "long long / HVCounts", "横線 {y, x1, x2} と縦線 {x, y1, y2} の交わる組の数 / 線ごとの数", "O((H + V) log)");
    RUN(const std::vector<HSeg> H = {{0, 0, 10}, {5, 2, 4}});
    RUN(const std::vector<VSeg> V = {{3, -1, 6}, {10, 0, 0}, {7, 1, 9}});
    SHOW(count_hv_intersections(H, V), "(3, 0)・(3, 5)・(10, 0)（端で触れる）");
    SHOW(hv_intersection_counts(H, V).per_h, "横線ごと");
    SHOW(hv_intersection_counts(H, V).per_v, "縦線ごと");

    api("【3】 max_overlap_intervals(iv) / max_overlap_rects(rs)", "Overlap1D {count, x} / Overlap2D {count, x, y}", "半開の区間・長方形が最も多く重なる点", "O(N log N)");
    SHOW(max_overlap_intervals({{1, 4}, {4, 6}, {2, 5}}).count, "x = 4 で [4, 6) と [2, 5)（[1, 4) は 4 を含まない）");
    RUN(const std::vector<Rect> rs = {{0, 0, 4, 4}, {2, 2, 6, 6}, {3, 0, 5, 3}});
    RUN(const auto o = max_overlap_rects(rs));
    SHOW(o.count, "");
    SHOW((std::pair<long long, long long>{o.x, o.y}), "3 つが重なる点の 1 つ");

    api("【4】 union_perimeter(rs)", "long long", "長方形の和集合の周の長さ（穴の周も含む）", "O(N log N)");
    SHOW(union_perimeter({{0, 0, 2, 2}, {1, 1, 3, 3}}), "2 × 2 を 2 つずらして重ねた形 = 12");
    SHOW(union_perimeter({{0, 0, 3, 1}, {0, 2, 3, 3}, {0, 0, 1, 3}, {2, 0, 3, 3}}), "真ん中に 1 × 1 の穴：外 12 + 穴 4");
}

}  // namespace sweep_demo

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
    sweep_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    sweep_test::benchmark();    // 速度計測
    sweep_demo::demo();         // 使い方デモ
    return 0;
}
