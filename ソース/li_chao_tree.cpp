// =============================================================================
//  Li Chao Tree（直線・線分を追加し、点 x での最小値 / 最大値を答える）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  直線 y = a x + b を追加していき、ある x での最小値（または最大値）を答える。CHT（Convex Hull Trick）の代わりに使える
//  傾きの順番・x の順番・追加と質問の順番に制約がない。x の範囲を限った「線分」も追加できる
//    例: dp[i] = min_j (dp[j] + (h[i] - h[j])^2) → j ごとに 直線 a = -2 h[j], b = dp[j] + h[j]^2 を足し、x = h[i] で query
//
//  2 種類（どちらも IS_MIN = true で最小値、false で最大値）
//    LiChaoTree<T, IS_MIN>        … 質問する x を先に全部渡す（座標圧縮）。速くてメモリも少ない。x は整数でも実数でもよい
//    DynamicLiChaoTree<T, IS_MIN> … x の範囲 [x_lo, x_hi)（long long の整数）だけ渡す。x を先読みできないとき
//    別名 : LiChaoTreeMin<T> / LiChaoTreeMax<T> / DynamicLiChaoTreeMin<T> / DynamicLiChaoTreeMax<T>（T の既定は long long）
//
//  機能一覧（x の範囲は半開区間 [l, r)。閉区間は _closed 版）         計算量（n = x の種類数、W = x_hi - x_lo）
//    LiChaoTreeMin<> lc(xs)              xs = 質問する x の一覧（重複・順不同でよい。vector<int> なども可）  O(n log n)
//    DynamicLiChaoTreeMin<> dl(lo, hi)   x は [lo, hi) の整数                                                 O(1)
//    add_line(a, b)                      直線 y = a x + b を追加                                              O(log n) / O(log W)
//    add_segment(l, r, a, b)             x ∈ [l, r) だけにある線分を追加（l >= r なら何もしない）              O(log² n) / O(log² W)
//    add_segment_closed(l, r, a, b)      x ∈ [l, r] だけにある線分を追加                                      O(log² n) / O(log² W)
//    query(x)                            x での最小値（Max 版は最大値）。線がないときは EMPTY                  O(log n) / O(log W)
//    query_line(x)                       x で最小（最大）になる直線 {a, b}。線がないときは {0, EMPTY}          O(log n) / O(log W)
//    clear()                             線を全部消す（x の一覧・範囲はそのまま）                             O(n) / O(1)
//    EMPTY                               線がないときの値（Min 版は INF = T の最大値、Max 版は -INF）
//    LiChaoTree だけ        : size()（x の種類数）/ coord(i)（i 番目に小さい x）/ query_index(i)（= query(coord(i))）
//    DynamicLiChaoTree だけ : node_count()（作った頂点の数。add_line で最大 1、add_segment で最大 約 4 log W 増える）
//
//  値の型 T : long long（ふつう）/ __int128（a x + b が long long を超えるとき）/ double（x や係数が実数のとき）
//  オーバーフロー注意 : a * x + b を T で計算する。質問しうる x すべてで |a x| ≤ 4e18 かつ |b| ≤ 4e18 なら long long で足りる
//                       （例: |a|, |x| ≤ 1e9 かつ |b| ≤ 4e18）。超えるなら T = __int128 にする
//  query(x) の x は、LiChaoTree では xs に入れた値、DynamicLiChaoTree では [x_lo, x_hi) の値（外ならメッセージを出して止める）
//  速さの目安（N = 2×10^5、x ∈ [-1e9, 1e9]、追加と質問を N 回ずつ）: 直線だけならどちらも 0.1〜0.2 秒
//    線分は LiChaoTree が 0.2〜0.3 秒。DynamicLiChaoTree は 0.7〜0.8 秒・頂点 約 700 万個（約 170 MB）になるので、
//    線分が多くて x を先読みできるなら LiChaoTree を使う
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    lichao_test::random_test()  … 素朴な全探索と全機能を突き合わせるランダムテスト
//    lichao_test::benchmark()    … 速度計測
//    lichao_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace lichao_internal {
template <class T>
constexpr T inf() {
    if constexpr (std::numeric_limits<T>::is_specialized) return std::numeric_limits<T>::max();
    else return (T)(~(unsigned __int128)0 >> 1);  // __int128 を gnu 拡張なし（-std=c++17）で使うとき
}
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace lichao_internal

// 直線 y = a x + b（query_line の返り値）
template <class T>
struct LiChaoLine {
    T a, b;
    T operator()(T x) const { return a * x + b; }
    bool operator==(const LiChaoLine& o) const { return a == o.a && b == o.b; }
};

// ---- 質問する x を先に全部渡す版 ----
template <class T = long long, bool IS_MIN = true>
struct LiChaoTree {
    using Line = LiChaoLine<T>;
    static constexpr T INF = lichao_internal::inf<T>();
    static constexpr T EMPTY = IS_MIN ? INF : -INF;  // 線がないときに query が返す値

    LiChaoTree() : LiChaoTree(std::vector<T>{}) {}
    // xs = 質問する x の一覧（重複・順不同でよい）
    template <class U>
    explicit LiChaoTree(const std::vector<U>& xs_) {
        xs.reserve(xs_.size());
        for (const U& x : xs_) xs.push_back((T)x);
        std::sort(xs.begin(), xs.end());
        xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
        n = (int)xs.size();
        sz = 1;
        while (sz < n) sz <<= 1;
        xs.resize(sz, n ? xs[n - 1] : T(0));  // 末尾は最後の x で埋める（同じ x が並んでも正しく動く）
        ln.assign(2 * (size_t)sz, empty_line());
    }

    int size() const { return n; }                                     // x の種類数
    T coord(int i) const { return check_index(i, "coord"), xs[i]; }    // i 番目に小さい x

    void add_line(T a, T b) { insert(make(a, b), 1, 0, sz); }                                      // y = a x + b
    void add_segment(T l, T r, T a, T b) { insert_range(make(a, b), lower(l), lower(r)); }         // x ∈ [l, r)
    void add_segment_closed(T l, T r, T a, T b) { insert_range(make(a, b), lower(l), upper(r)); }  // x ∈ [l, r]

    // x での最小値（Max 版は最大値）。線がないときは EMPTY
    T query(T x) const { return query_index(index_of(x, "query")); }
    T query_index(int i) const {
        check_index(i, "query_index");
        const T x = xs[i];
        T res = INF;
        for (int k = i + sz; k >= 1; k >>= 1) {
            const T v = ln[k](x);
            if (v < res) res = v;
        }
        return IS_MIN ? res : -res;
    }
    // x で最小（最大）になる直線。線がないときは {0, EMPTY}
    Line query_line(T x) const {
        const int i = index_of(x, "query_line");
        const T xv = xs[i];
        Line best = empty_line();
        T res = INF;
        for (int k = i + sz; k >= 1; k >>= 1) {
            const T v = ln[k](xv);
            if (v < res) res = v, best = ln[k];
        }
        return IS_MIN ? best : Line{-best.a, -best.b};
    }
    void clear() { std::fill(ln.begin(), ln.end(), empty_line()); }

   private:
    int n = 0, sz = 1;
    std::vector<T> xs;     // 質問する x（昇順・重複なし。長さ sz まで最後の値で埋める）
    std::vector<Line> ln;  // ln[k] = 頂点 k の直線（Max 版は符号を反転して最小値として持つ）

    static Line empty_line() { return Line{T(0), INF}; }
    static Line make(T a, T b) { return IS_MIN ? Line{a, b} : Line{-a, -b}; }
    // xs[i] >= x（upper は xs[i] > x）となる最小の i（なければ n）
    // 長さ sz（2 のべき乗）を分岐なしで二分探索（std::lower_bound は _GLIBCXX_DEBUG だと 1 回 O(n) になるため自前）
    int lower(T x) const {
        const T* p = xs.data();
        for (int len = sz; len > 1; len >>= 1) p += (p[(len >> 1) - 1] < x) ? (len >> 1) : 0;
        const int i = (int)(p - xs.data()) + (*p < x);
        return i < n ? i : n;
    }
    int upper(T x) const {
        const T* p = xs.data();
        for (int len = sz; len > 1; len >>= 1) p += (p[(len >> 1) - 1] <= x) ? (len >> 1) : 0;
        const int i = (int)(p - xs.data()) + (*p <= x);
        return i < n ? i : n;
    }
    int index_of(T x, const char* func) const {
        const int i = lower(x);
        if (i == n || xs[i] != x) lichao_internal::fail(func, "x が最初に渡した x の一覧にない");
        return i;
    }
    void check_index(int i, const char* func) const {
        if (i < 0 || i >= n) lichao_internal::fail(func, "i が範囲 [0, size()) の外");
    }
    // 頂点 k（添字 [l, r) の x を担当）より下に直線 f を入れる
    void insert(Line f, int k, int l, int r) {
        while (true) {
            Line& g = ln[k];
            const T xl = xs[l], xr = xs[r - 1];
            const T fl = f(xl), gl = g(xl), fr = f(xr), gr = g(xr);
            if (fl >= gl && fr >= gr) return;              // f はこの範囲のどこでも g 以上 → 要らない
            if (fl <= gl && fr <= gr) return void(g = f);  // f はこの範囲のどこでも g 以下 → g が要らない
            // 範囲の中で交わる（このとき r - l >= 2）→ 真ん中で小さいほうを頂点に残し、負けたほうを片側の子へ
            const int m = (l + r) >> 1;
            const bool swapped = f(xs[m]) < g(xs[m]);
            if (swapped) std::swap(f, g);
            if ((fl < gl) != swapped) k = 2 * k, r = m;  // 負けたほうが小さくなりうるのは左の子 [l, m) だけ
            else k = 2 * k + 1, l = m;                   // 右の子 [m, r) だけ
        }
    }
    // 添字 [l, r) の x だけに直線 f を入れる（区間を O(log n) 個の頂点に分ける）
    void insert_range(const Line& f, int l, int r) {
        if (l >= r) return;
        int lo = l + sz, hi = r + sz, w = 1;
        while (lo < hi) {
            if (lo & 1) insert(f, lo, lo * w - sz, (lo + 1) * w - sz), lo++;
            if (hi & 1) --hi, insert(f, hi, hi * w - sz, (hi + 1) * w - sz);
            lo >>= 1, hi >>= 1, w <<= 1;
        }
    }
};

// ---- x の範囲 [x_lo, x_hi)（整数）だけ渡す版（使う頂点だけ作る） ----
template <class T = long long, bool IS_MIN = true>
struct DynamicLiChaoTree {
    using Line = LiChaoLine<T>;
    static constexpr T INF = lichao_internal::inf<T>();
    static constexpr T EMPTY = IS_MIN ? INF : -INF;  // 線がないときに query が返す値

    // x は [x_lo, x_hi) の整数。reserve = 頂点数の見込み（省略可）
    DynamicLiChaoTree(long long x_lo, long long x_hi, int reserve = 0) : lo(x_lo), hi(x_hi) {
        if (!(x_lo < x_hi)) lichao_internal::fail("DynamicLiChaoTree", "x_lo < x_hi でない");
        nodes.reserve(reserve > 1 ? reserve : 1);
        nodes.push_back(Node{empty_line(), {-1, -1}});
    }

    void add_line(T a, T b) { insert(make(a, b), 0, lo, hi); }  // y = a x + b
    void add_segment(long long l, long long r, T a, T b) {      // x ∈ [l, r)（[x_lo, x_hi) の外は無視）
        l = std::max(l, lo), r = std::min(r, hi);
        if (l < r) insert_range(make(a, b), 0, lo, hi, l, r);
    }
    void add_segment_closed(long long l, long long r, T a, T b) { add_segment(l, r >= hi ? hi : r + 1, a, b); }  // x ∈ [l, r]

    // x での最小値（Max 版は最大値）。線がないときは EMPTY
    T query(long long x) const {
        check_x(x, "query");
        T res = INF;
        long long l = lo, r = hi;
        for (int k = 0; k >= 0;) {
            const T v = nodes[k].f(x);
            if (v < res) res = v;
            const long long m = l + (r - l) / 2;
            if (x < m) k = nodes[k].ch[0], r = m;
            else k = nodes[k].ch[1], l = m;
        }
        return IS_MIN ? res : -res;
    }
    // x で最小（最大）になる直線。線がないときは {0, EMPTY}
    Line query_line(long long x) const {
        check_x(x, "query_line");
        Line best = empty_line();
        T res = INF;
        long long l = lo, r = hi;
        for (int k = 0; k >= 0;) {
            const T v = nodes[k].f(x);
            if (v < res) res = v, best = nodes[k].f;
            const long long m = l + (r - l) / 2;
            if (x < m) k = nodes[k].ch[0], r = m;
            else k = nodes[k].ch[1], l = m;
        }
        return IS_MIN ? best : Line{-best.a, -best.b};
    }
    void clear() { nodes.assign(1, Node{empty_line(), {-1, -1}}); }
    int node_count() const { return (int)nodes.size(); }

   private:
    struct Node {
        Line f;     // この頂点の直線（Max 版は符号を反転して最小値として持つ）
        int ch[2];  // 子（ない子は -1）
    };
    long long lo, hi;
    std::vector<Node> nodes;  // nodes[0] が根（担当 [lo, hi)）

    static Line empty_line() { return Line{T(0), INF}; }
    static Line make(T a, T b) { return IS_MIN ? Line{a, b} : Line{-a, -b}; }
    void check_x(long long x, const char* func) const {
        if (x < lo || x >= hi) lichao_internal::fail(func, "x が範囲 [x_lo, x_hi) の外");
    }
    // 頂点 k（担当 [l, r)）より下に直線 f を入れる。子がなければ作ってそこに置く
    void insert(Line f, int k, long long l, long long r) {
        while (true) {
            Line& g = nodes[k].f;
            const T fl = f(l), gl = g(l), fr = f(r - 1), gr = g(r - 1);
            if (fl >= gl && fr >= gr) return;              // f はこの範囲のどこでも g 以上 → 要らない
            if (fl <= gl && fr <= gr) return void(g = f);  // f はこの範囲のどこでも g 以下 → g が要らない
            // 範囲の中で交わる（このとき r - l >= 2）→ 真ん中で小さいほうを頂点に残し、負けたほうを片側の子へ
            const long long m = l + (r - l) / 2;
            const bool swapped = f(m) < g(m);
            if (swapped) std::swap(f, g);
            int side;
            if ((fl < gl) != swapped) side = 0, r = m;  // 左の子 [l, m) だけ
            else side = 1, l = m;                       // 右の子 [m, r) だけ
            const int c = nodes[k].ch[side];
            if (c < 0) {
                nodes[k].ch[side] = (int)nodes.size();
                nodes.push_back(Node{f, {-1, -1}});
                return;
            }
            k = c;
        }
    }
    // [sl, sr) を頂点に分けて入れる（頂点 k の担当 [l, r) は [sl, sr) と交わる）
    void insert_range(const Line& f, int k, long long l, long long r, long long sl, long long sr) {
        if (sl <= l && r <= sr) return insert(f, k, l, r);
        const long long m = l + (r - l) / 2;
        if (sl < m) insert_range(f, child(k, 0), l, m, sl, sr);
        if (m < sr) insert_range(f, child(k, 1), m, r, sl, sr);
    }
    int child(int k, int side) {
        if (nodes[k].ch[side] < 0) {
            nodes[k].ch[side] = (int)nodes.size();
            nodes.push_back(Node{empty_line(), {-1, -1}});
        }
        return nodes[k].ch[side];
    }
};

template <class T = long long>
using LiChaoTreeMin = LiChaoTree<T, true>;
template <class T = long long>
using LiChaoTreeMax = LiChaoTree<T, false>;
template <class T = long long>
using DynamicLiChaoTreeMin = DynamicLiChaoTree<T, true>;
template <class T = long long>
using DynamicLiChaoTreeMax = DynamicLiChaoTree<T, false>;
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace lichao_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define LC_CHECK(...)                                                                                \
    do {                                                                                             \
        if (!(__VA_ARGS__)) {                                                                        \
            std::cerr << "[li_chao_tree] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                            \
        }                                                                                            \
    } while (0)

using i128 = __int128;

// [lo, hi] の一様乱数（long long の端まで使ってよい）
static long long rnd(std::mt19937_64& rng, long long lo, long long hi) {
    const unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return (long long)((unsigned long long)lo + (w == 0 ? rng() : rng() % w));
}
static i128 rnd128(std::mt19937_64& rng, i128 lim) {  // [-lim, lim]（lim < 2^100 くらいまで）
    const i128 v = ((i128)(rng() >> 1) << 63 | (i128)(rng() >> 1)) % (2 * lim + 1);
    return v - lim;
}

// 素朴解：追加した直線・線分を全部持ち、x ごとに全部見る（値は W で計算）
template <class T, class W>
struct Naive {
    struct Seg {
        T a, b, l, r;
        int kind;  // 0: 直線, 1: [l, r), 2: [l, r]
    };
    bool is_min;
    std::vector<Seg> v;
    bool covers(const Seg& s, T x) const {
        if (s.kind == 0) return true;
        if (s.kind == 1) return s.l <= x && x < s.r;
        return s.l <= x && x <= s.r;
    }
    // x で有効な線があれば最良値を out に入れて true
    bool best(T x, W& out) const {
        bool found = false;
        for (const Seg& s : v) {
            if (!covers(s, x)) continue;
            const W y = (W)s.a * (W)x + (W)s.b;
            if (!found || (is_min ? y < out : y > out)) out = y, found = true;
        }
        return found;
    }
    // 直線 {a, b} が x で有効な線の中にあるか
    bool has(T x, T a, T b) const {
        for (const Seg& s : v)
            if (covers(s, x) && s.a == a && s.b == b) return true;
        return false;
    }
};

// query / query_line / query_index を素朴解と突き合わせる
template <class Tree, class T, class W>
static void check_point(const Tree& tr, const Naive<T, W>& nv, T x, int& cases) {
    W b{};
    const bool f = nv.best(x, b);
    const T q = tr.query(x);
    const LiChaoLine<T> L = tr.query_line(x);
    if (!f) {
        LC_CHECK(q == Tree::EMPTY);
        LC_CHECK(L.a == T(0) && L.b == Tree::EMPTY);
    } else {
        LC_CHECK((W)q == b);
        LC_CHECK((W)L.a * (W)x + (W)L.b == b);
        LC_CHECK(nv.has(x, L.a, L.b));
    }
    ++cases;
}

// LiChaoTree：x の一覧を作り、追加と質問をランダムに混ぜる
template <class T, bool IS_MIN, class W, class GX, class GA, class GB>
static void test_static(std::mt19937_64& rng, int n_max, int ops, GX gen_x, GA gen_a, GB gen_b, int& cases) {
    const int n0 = (int)rnd(rng, 0, n_max);
    std::vector<T> xs(n0);
    for (auto& x : xs) x = gen_x();
    LiChaoTree<T, IS_MIN> lc(xs);
    std::vector<T> ux = xs;
    std::sort(ux.begin(), ux.end());
    ux.erase(std::unique(ux.begin(), ux.end()), ux.end());
    const int n = (int)ux.size();
    LC_CHECK(lc.size() == n);
    for (int i = 0; i < n; ++i) LC_CHECK(lc.coord(i) == ux[i]);
    Naive<T, W> nv{IS_MIN, {}};
    // 線分の端は「一覧の x ちょうど」と「一覧にない値」の両方を出す
    auto gen_end = [&]() -> T { return (n > 0 && rnd(rng, 0, 1)) ? ux[rnd(rng, 0, n - 1)] : gen_x(); };
    for (int t = 0; t < ops; ++t) {
        const int type = (int)rnd(rng, 0, 11);
        if (type <= 2) {
            const T a = gen_a(), b = gen_b();
            lc.add_line(a, b);
            nv.v.push_back({a, b, T(0), T(0), 0});
        } else if (type <= 5) {
            T l = gen_end(), r = gen_end();
            if (rnd(rng, 0, 4) && r < l) std::swap(l, r);  // ときどき l > r（何もしない）も試す
            const T a = gen_a(), b = gen_b();
            const bool closed = type == 5;
            if (closed) lc.add_segment_closed(l, r, a, b);
            else lc.add_segment(l, r, a, b);
            nv.v.push_back({a, b, l, r, closed ? 2 : 1});
        } else if (type == 6 && rnd(rng, 0, 30) == 0) {
            lc.clear();
            nv.v.clear();
        } else if (n > 0) {
            const int i = (int)rnd(rng, 0, n - 1);
            check_point(lc, nv, ux[i], cases);
            LC_CHECK(lc.query_index(i) == lc.query(ux[i]));
        }
    }
    for (int i = 0; i < n; ++i) {
        check_point(lc, nv, ux[i], cases);
        LC_CHECK(lc.query_index(i) == lc.query(ux[i]));
    }
}

// DynamicLiChaoTree：範囲 [lo, hi) で追加と質問をランダムに混ぜる（all = true なら最後に全 x を確認）
template <class T, bool IS_MIN, class W, class GA, class GB>
static void test_dynamic(std::mt19937_64& rng, long long lo, long long hi, int ops, bool all, GA gen_a, GB gen_b, int& cases) {
    DynamicLiChaoTree<T, IS_MIN> dl(lo, hi);
    Naive<T, W> nv{IS_MIN, {}};
    auto gen_x = [&]() { return rnd(rng, lo, hi - 1); };
    // 線分の端：範囲の中・端ちょうど・少し外・long long の端
    auto gen_end = [&]() -> long long {
        const int t = (int)rnd(rng, 0, 9);
        if (t <= 5) return gen_x();
        if (t == 6) return rnd(rng, 0, 1) ? lo : hi;
        if (t == 7) return rnd(rng, 0, 1) ? lo - (lo > std::numeric_limits<long long>::min() + 3 ? rnd(rng, 1, 3) : 0)
                                          : hi + (hi < std::numeric_limits<long long>::max() - 3 ? rnd(rng, 1, 3) : 0);
        if (t == 8) return rnd(rng, 0, 1) ? std::numeric_limits<long long>::min() : std::numeric_limits<long long>::max();
        return rnd(rng, 0, 1) ? hi - 1 : lo + 1;
    };
    for (int t = 0; t < ops; ++t) {
        const int type = (int)rnd(rng, 0, 11);
        if (type <= 2) {
            const T a = gen_a(), b = gen_b();
            dl.add_line(a, b);
            nv.v.push_back({a, b, T(0), T(0), 0});
        } else if (type <= 5) {
            long long l = gen_end(), r = gen_end();
            if (rnd(rng, 0, 4) && r < l) std::swap(l, r);
            const T a = gen_a(), b = gen_b();
            const bool closed = type == 5;
            if (closed) dl.add_segment_closed(l, r, a, b);
            else dl.add_segment(l, r, a, b);
            nv.v.push_back({a, b, (T)l, (T)r, closed ? 2 : 1});
        } else if (type == 6 && rnd(rng, 0, 30) == 0) {
            dl.clear();
            nv.v.clear();
            LC_CHECK(dl.node_count() == 1);
        } else {
            const long long x = rnd(rng, 0, 3) == 0 ? (rnd(rng, 0, 1) ? lo : hi - 1) : gen_x();
            check_point(dl, nv, (T)x, cases);
        }
    }
    if (all)
        for (long long x = lo; x < hi; ++x) check_point(dl, nv, (T)x, cases);
}

// DP の使い方の確認（EDPC Z Frog 3）：O(N^2) の DP と一致するか
static void test_frog3(std::mt19937_64& rng, int& cases) {
    for (int it = 0; it < 200 * TEST_SCALE; ++it) {
        const int N = (int)rnd(rng, 2, 60);
        const long long C = rnd(rng, 1, 1000000000000LL);
        std::vector<long long> h(N);
        long long cur = 0;
        for (auto& x : h) x = (cur += rnd(rng, 1, 20000));  // 狭義単調増加（最大 1.2e6）
        std::vector<long long> naive(N, 0), d1(N, 0), d2(N, 0);
        for (int i = 1; i < N; ++i) {
            naive[i] = std::numeric_limits<long long>::max();
            for (int j = 0; j < i; ++j) naive[i] = std::min(naive[i], naive[j] + (h[i] - h[j]) * (h[i] - h[j]) + C);
        }
        LiChaoTreeMin<> f(h);
        DynamicLiChaoTreeMin<> g(1, cur + 1);
        for (int i = 0; i < N; ++i) {
            if (i > 0) d1[i] = f.query(h[i]) + h[i] * h[i] + C, d2[i] = g.query(h[i]) + h[i] * h[i] + C;
            f.add_line(-2 * h[i], d1[i] + h[i] * h[i]);
            g.add_line(-2 * h[i], d2[i] + h[i] * h[i]);
        }
        LC_CHECK(d1 == naive && d2 == naive);
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261019);
    int cases = 0;
    using ll = long long;

    // ---- LiChaoTree（long long・小さい値。同じ値・境界ちょうどが出やすい）----
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const ll X = rnd(rng, 1, 60);
        auto gx = [&]() { return (ll)rnd(rng, -X, X); };
        auto ga = [&]() { return (ll)rnd(rng, -20, 20); };
        auto gb = [&]() { return (ll)rnd(rng, -1000, 1000); };
        test_static<ll, true, i128>(rng, 70, 200, gx, ga, gb, cases);
        test_static<ll, false, i128>(rng, 70, 200, gx, ga, gb, cases);
    }
    // ---- LiChaoTree（long long・大きい値：|a|, |x| ≤ 1e9、|b| ≤ 4e18）----
    for (int it = 0; it < 60 * TEST_SCALE; ++it) {
        auto gx = [&]() { return (ll)rnd(rng, -1000000000, 1000000000); };
        auto ga = [&]() { return (ll)rnd(rng, -1000000000, 1000000000); };
        auto gb = [&]() { return (ll)rnd(rng, -4000000000000000000LL, 4000000000000000000LL); };
        test_static<ll, true, i128>(rng, 300, 600, gx, ga, gb, cases);
        test_static<ll, false, i128>(rng, 300, 600, gx, ga, gb, cases);
    }
    // ---- LiChaoTree（__int128：|a|, |x| ≤ 1e18、|b| ≤ 1e36。座標は vector<long long> で渡す）----
    for (int it = 0; it < 30 * TEST_SCALE; ++it) {
        const ll LIM = 1000000000000000000LL;
        auto ga = [&]() { return (i128)rnd(rng, -LIM, LIM); };
        auto gb = [&]() { return rnd128(rng, (i128)LIM * LIM); };
        auto gx = [&]() { return (i128)rnd(rng, -LIM, LIM); };
        test_static<i128, true, i128>(rng, 200, 400, gx, ga, gb, cases);
        test_static<i128, false, i128>(rng, 200, 400, gx, ga, gb, cases);
        std::vector<ll> xs = {-LIM, 0, LIM};
        LiChaoTreeMin<i128> t(xs);  // vector<long long> から作れる
        t.add_line((i128)LIM, (i128)LIM);
        LC_CHECK(t.query(LIM) == (i128)LIM * LIM + LIM && t.query(-LIM) == -(i128)LIM * LIM + LIM);
    }
    // ---- LiChaoTree（double・整数の値なので誤差なしで比べられる）----
    for (int it = 0; it < 60 * TEST_SCALE; ++it) {
        auto gx = [&]() { return (double)rnd(rng, -40, 40) / 4; };
        auto ga = [&]() { return (double)rnd(rng, -20, 20); };
        auto gb = [&]() { return (double)rnd(rng, -1000, 1000); };
        test_static<double, true, double>(rng, 60, 200, gx, ga, gb, cases);
        test_static<double, false, double>(rng, 60, 200, gx, ga, gb, cases);
    }
    // ---- LiChaoTree：x が 0 個・1 個、2 のべき乗ちょうど / ±1 ----
    {
        LiChaoTreeMin<> e(std::vector<ll>{});
        e.add_line(1, 2), e.add_segment(0, 5, 1, 1), e.add_segment_closed(0, 5, 1, 1);  // 何も起きない
        LC_CHECK(e.size() == 0);
        LiChaoTreeMax<> one(std::vector<int>{7, 7, 7});
        LC_CHECK(one.size() == 1 && one.query(7) == LiChaoTreeMax<>::EMPTY);
        one.add_segment(7, 7, 1, 0);  // 空の区間 [7, 7)
        LC_CHECK(one.query(7) == LiChaoTreeMax<>::EMPTY);
        one.add_segment_closed(7, 7, 1, 0);  // [7, 7] は 7 を含む
        LC_CHECK(one.query(7) == 7);
        ++cases;
        for (int n : {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 64, 65}) {
            for (int rep = 0; rep < 4 * TEST_SCALE; ++rep) {
                std::vector<ll> xs(n);
                for (int i = 0; i < n; ++i) xs[i] = 3 * i - n;
                std::shuffle(xs.begin(), xs.end(), rng);
                LiChaoTreeMin<> lc(xs);
                Naive<ll, i128> nv{true, {}};
                for (int t = 0; t < 3 * n + 5; ++t) {
                    const ll a = rnd(rng, -10, 10), b = rnd(rng, -100, 100);
                    if (rnd(rng, 0, 1)) {
                        lc.add_line(a, b), nv.v.push_back({a, b, 0, 0, 0});
                    } else {
                        ll l = rnd(rng, -n - 2, 2 * n + 2), r = rnd(rng, -n - 2, 2 * n + 2);
                        if (l > r) std::swap(l, r);
                        lc.add_segment(l, r, a, b), nv.v.push_back({a, b, l, r, 1});
                    }
                }
                for (int i = 0; i < n; ++i) check_point(lc, nv, lc.coord(i), cases);
            }
        }
    }

    // ---- DynamicLiChaoTree（小さい範囲：最後に全 x を確認）----
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const ll lo = rnd(rng, -40, 40), hi = lo + rnd(rng, 1, 70);
        auto ga = [&]() { return (ll)rnd(rng, -20, 20); };
        auto gb = [&]() { return (ll)rnd(rng, -1000, 1000); };
        test_dynamic<ll, true, i128>(rng, lo, hi, 200, true, ga, gb, cases);
        test_dynamic<ll, false, i128>(rng, lo, hi, 200, true, ga, gb, cases);
    }
    // ---- DynamicLiChaoTree（x ∈ [-1e9, 1e9]、|a| ≤ 1e9、|b| ≤ 4e18）----
    for (int it = 0; it < 60 * TEST_SCALE; ++it) {
        auto ga = [&]() { return (ll)rnd(rng, -1000000000, 1000000000); };
        auto gb = [&]() { return (ll)rnd(rng, -4000000000000000000LL, 4000000000000000000LL); };
        test_dynamic<ll, true, i128>(rng, -1000000000, 1000000001, 600, false, ga, gb, cases);
        test_dynamic<ll, false, i128>(rng, -1000000000, 1000000001, 600, false, ga, gb, cases);
    }
    // ---- DynamicLiChaoTree（範囲の幅が long long の最大値ちょうど：[-2^62, 2^62 - 1)）----
    for (int it = 0; it < 30 * TEST_SCALE; ++it) {
        const ll lo = -(1LL << 62), hi = (1LL << 62) - 1;
        auto ga = [&]() { return (ll)rnd(rng, -1, 1); };
        auto gb = [&]() { return (ll)rnd(rng, -1000000000000000000LL, 1000000000000000000LL); };
        test_dynamic<ll, true, i128>(rng, lo, hi, 400, false, ga, gb, cases);
        test_dynamic<ll, false, i128>(rng, lo, hi, 400, false, ga, gb, cases);
        auto ga2 = [&]() { return (i128)rnd(rng, -1000000000000000000LL, 1000000000000000000LL); };
        auto gb2 = [&]() { return rnd128(rng, (i128)1000000000000000000LL * 1000000000000000000LL); };
        test_dynamic<i128, true, i128>(rng, lo, hi, 400, false, ga2, gb2, cases);
        test_dynamic<i128, false, i128>(rng, lo, hi, 400, false, ga2, gb2, cases);
    }
    // ---- DynamicLiChaoTree：幅 1 の範囲 ----
    {
        DynamicLiChaoTreeMin<> d(5, 6);
        LC_CHECK(d.query(5) == DynamicLiChaoTreeMin<>::EMPTY);
        d.add_segment(0, 5, 1, 0);  // 範囲 [5, 6) と交わらない
        LC_CHECK(d.query(5) == DynamicLiChaoTreeMin<>::EMPTY);
        d.add_segment_closed(0, 5, 1, 0), d.add_line(2, -1);
        LC_CHECK(d.query(5) == 5 && d.query_line(5) == (LiChaoLine<ll>{1, 0}));
        ++cases;
    }

    test_frog3(rng, cases);
    std::cout << "[li_chao_tree] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[li_chao_tree] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    auto ms = [](std::chrono::steady_clock::time_point a) {
        return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    };
    const int N = TEST_SCALE == 1 ? 20000 : 200000;  // デバッグビルドでは 1/10
    const long long X = 1000000000;
    std::vector<long long> xs(N), qa(N), qb(N), ql(N), qr(N), qx(N);
    // 放物線 y = -x^2 の接線（どの直線も下側の包絡線に残る、Li Chao Tree にとって重い入力）
    for (int i = 0; i < N; ++i) {
        const long long t = rnd(rng, -X, X);
        xs[i] = rnd(rng, -X, X), qa[i] = -2 * t, qb[i] = t * t;
        ql[i] = rnd(rng, -X, X), qr[i] = rnd(rng, -X, X);
        if (ql[i] > qr[i]) std::swap(ql[i], qr[i]);
    }
    for (int i = 0; i < N; ++i) qx[i] = xs[rng() % N];
    unsigned long long chk = 0;  // 計算を消されないための合計（あふれてよいので unsigned）

    auto t0 = std::chrono::steady_clock::now();
    {
        LiChaoTreeMin<> lc(xs);
        for (int i = 0; i < N; ++i) lc.add_line(qa[i], qb[i]), chk += (unsigned long long)lc.query(qx[i]);
    }
    const long long t1 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    {
        LiChaoTreeMin<> lc(xs);
        for (int i = 0; i < N; ++i) lc.add_segment(ql[i], qr[i], qa[i], qb[i]), chk += (unsigned long long)lc.query(qx[i]);
    }
    const long long t2 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    {
        DynamicLiChaoTreeMin<> dl(-X, X + 1);
        for (int i = 0; i < N; ++i) dl.add_line(qa[i], qb[i]), chk += (unsigned long long)dl.query(qx[i]);
    }
    const long long t3 = ms(t0);
    t0 = std::chrono::steady_clock::now();
    int nodes = 0;
    {
        DynamicLiChaoTreeMin<> dl(-X, X + 1);
        for (int i = 0; i < N; ++i) dl.add_segment(ql[i], qr[i], qa[i], qb[i]), chk += (unsigned long long)dl.query(qx[i]);
        nodes = dl.node_count();
    }
    const long long t4 = ms(t0);
    std::cout << "[li_chao_tree] benchmark (N = " << N << "、x ∈ [-1e9, 1e9]、追加と質問を N 回ずつ): LiChaoTree 直線 " << t1 << " ms / 線分 " << t2
              << " ms, DynamicLiChaoTree 直線 " << t3 << " ms / 線分 " << t4 << " ms（頂点 " << nodes << " 個）  (check " << chk % 7 << ")\n";
}

}  // namespace lichao_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace lichao_demo {

[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(double v) {
    std::ostringstream os;
    os << std::setprecision(12) << v;
    return os.str();
}
template <class T>
[[maybe_unused]] static std::string to_s(const T& v) {  // << で出せるもの（整数など）
    std::ostringstream os;
    os << v;
    return os.str();
}
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "]";
}
template <class T>
[[maybe_unused]] static std::string to_s(const LiChaoLine<T>& f) {
    return "{a = " + to_s(f.a) + ", b = " + to_s(f.b) + "}";
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
#define SHOW(expr, comment) lichao_demo::show(#expr, (expr), comment)
#define RUN(...)                    \
    lichao_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("Li Chao Tree（li_chao_tree） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  LiChaoTreeMin<T> / LiChaoTreeMax<T>                … 質問する x を先に全部渡す版（T の既定は long long）\n"
              << "  DynamicLiChaoTreeMin<T> / DynamicLiChaoTreeMax<T>  … x の範囲 [x_lo, x_hi) だけ渡す版\n"
              << "  add_line / add_segment / add_segment_closed / query / query_line / clear / EMPTY\n"
              << "  LiChaoTree だけ: size / coord / query_index      DynamicLiChaoTree だけ: node_count\n";

    api("【1】 LiChaoTreeMin<> lc(xs)", "—", "質問する x の一覧 xs を渡して作る（重複・順不同でよい）。x での最小値を答える", "O(n log n)");
    RUN(const std::vector<long long> xs = {5, -3, 2, 0, 2});
    RUN(LiChaoTreeMin<> lc(xs));
    SHOW(lc.size(), "x の種類数（2 の重複は 1 つ）");
    SHOW(lc.query(0), "線がまだない → EMPTY（= INF = long long の最大値）");
    SHOW(lc.query(0) == lc.EMPTY, "線があるかどうかは EMPTY と比べる");

    api("【2】 add_line(a, b) / query(x)", "void / T", "直線 y = a x + b を追加 / x での最小値（x は xs に入れた値）", "O(log n) / O(log n)");
    RUN(lc.add_line(2, 1));   // y = 2x + 1
    RUN(lc.add_line(-1, 4));  // y = -x + 4
    RUN(lc.add_line(0, 3));   // y = 3
    SHOW(lc.query(-3), "min(2·(-3)+1, 3+4, 3) = -5");
    SHOW(lc.query(0), "min(1, 4, 3) = 1");
    SHOW(lc.query(2), "min(5, 2, 3) = 2");
    SHOW(lc.query(5), "min(11, -1, 3) = -1");

    api("【3】 query_line(x)", "LiChaoLine<T>（.a と .b を持つ）", "x で最小になる直線。線がないときは {0, EMPTY}", "O(log n)");
    SHOW(lc.query_line(5), "x = 5 では y = -x + 4 が最小");
    SHOW(lc.query_line(-3).a, "傾きだけ取り出す（y = 2x + 1 の 2）");

    api("【4】 add_segment(l, r, a, b) / add_segment_closed(l, r, a, b)", "void",
        "x ∈ [l, r)（_closed は [l, r]）だけにある線分を追加。l, r は xs にない値でもよい。l >= r なら何もしない", "O(log² n)");
    RUN(lc.add_segment(-1, 2, 0, -10));  // x ∈ [-1, 2) で y = -10（xs のうち 0 だけが入る）
    SHOW(lc.query(0), "x = 0 は線分の中 → -10");
    SHOW(lc.query(2), "x = 2 は [-1, 2) の外（右端は含まない）→ 2 のまま");
    RUN(lc.add_segment_closed(2, 5, 1, -20));  // x ∈ [2, 5] で y = x - 20
    SHOW(lc.query(2), "2 - 20 = -18");
    SHOW(lc.query(5), "5 - 20 = -15（_closed なので右端の 5 も含む）");
    SHOW(lc.query(-3), "-3 は線分の外 → -5 のまま");

    api("【5】 size() / coord(i) / query_index(i)（LiChaoTree だけ）", "int / T / T",
        "x の種類数 / i 番目に小さい x / query(coord(i)) と同じ（x を探す二分探索を省ける）", "O(1) / O(1) / O(log n)");
    {
        code("for (int i = 0; i < lc.size(); i++)  → lc.coord(i) : lc.query_index(i)");
        std::string s;
        for (int i = 0; i < lc.size(); i++) s += to_s(lc.coord(i)) + " : " + to_s(lc.query_index(i)) + "   ";
        note("→ " + s);
    }

    api("【6】 clear()", "void", "線を全部消す（x の一覧はそのまま）", "O(n)");
    RUN(lc.clear());
    SHOW(lc.query(5) == lc.EMPTY, "線がなくなった");

    api("【7】 LiChaoTreeMax<> lx(xs)", "—", "最大値を答える版（使い方は Min 版と同じ）。線がないときの EMPTY は -INF", "Min 版と同じ");
    RUN(LiChaoTreeMax<> lx(xs));
    SHOW(lx.query(0), "線がまだない → EMPTY（= -INF）");
    RUN(lx.add_line(2, 1), lx.add_line(-1, 4), lx.add_line(0, 3));
    SHOW(lx.query(-3), "max(-5, 7, 3) = 7");
    SHOW(lx.query(0), "max(1, 4, 3) = 4");
    SHOW(lx.query(5), "max(11, -1, 3) = 11");
    SHOW(lx.query_line(0), "x = 0 では y = -x + 4 が最大");

    api("【8】 DynamicLiChaoTreeMin<> dl(x_lo, x_hi)", "—",
        "x の範囲 [x_lo, x_hi)（整数）だけ渡して作る。x を先読みできない（質問が来るまで x がわからない）ときに使う", "O(1)（線を足すたびに頂点が増える）");
    RUN(DynamicLiChaoTreeMin<> dl(-1000000000, 1000000001));  // x ∈ [-1e9, 1e9]
    RUN(dl.add_line(3, 0), dl.add_line(-2, 10));
    SHOW(dl.query(1), "min(3, 8) = 3（範囲の中ならどの x でも質問できる）");
    SHOW(dl.query(1000000000), "min(3e9, -2e9 + 10)");
    RUN(dl.add_segment_closed(-5, 5, 0, -100));  // x ∈ [-5, 5] で y = -100
    SHOW(dl.query(5), "線分 [-5, 5] の中 → -100");
    SHOW(dl.query(6), "線分の外 → min(18, -2) = -2");
    SHOW(dl.query_line(6), "x = 6 で最小の直線 y = -2x + 10");
    SHOW(dl.node_count(), "作った頂点の数（メモリの目安）");
    RUN(dl.clear());
    SHOW(dl.query(0) == dl.EMPTY, "線を全部消した");
    note("最大値は DynamicLiChaoTreeMax<> を使う（使い方は同じ）");

    api("【9】 使い方の例: DP の高速化（EDPC Z Frog 3）", "—",
        "dp[i] = min_{j<i} (dp[j] + (h[i] - h[j])²) + C。展開すると dp[j] + h[j]² - 2 h[j]·h[i] + h[i]² + C なので、"
        "j ごとに 直線 a = -2 h[j], b = dp[j] + h[j]² を足し、x = h[i] で query する",
        "O(N log N)");
    {
        code("N = 5, C = 6, h = {1, 2, 3, 4, 5}   // EDPC Z の入力例 1（答え 20）");
        code("LiChaoTreeMin<> f(h);");
        code("for (int i = 0; i < N; i++) {");
        code("    if (i > 0) dp[i] = f.query(h[i]) + h[i] * h[i] + C;");
        code("    f.add_line(-2 * h[i], dp[i] + h[i] * h[i]);");
        code("}");
        const long long C = 6;
        const std::vector<long long> h = {1, 2, 3, 4, 5};
        const int N = (int)h.size();
        LiChaoTreeMin<> f(h);
        std::vector<long long> dp(N, 0);
        for (int i = 0; i < N; i++) {
            if (i > 0) dp[i] = f.query(h[i]) + h[i] * h[i] + C;
            f.add_line(-2 * h[i], dp[i] + h[i] * h[i]);
        }
        note("→ dp = " + to_s(dp) + "、答え dp[N - 1] = " + to_s(dp.back()));
    }

    api("【10】 値が long long を超えるとき: LiChaoTreeMin<__int128>", "__int128",
        "a x + b を __int128 で計算する。xs は vector<long long> のまま渡せる（cout できないので表示するときは変換する）", "long long 版と同じ（少し遅い）");
    RUN(LiChaoTreeMin<__int128> big(std::vector<long long>{-1000000000000000000LL, 1000000000000000000LL}));
    RUN(big.add_line(4000000000000000000LL, 5));  // 4e18 × 1e18 = 4e36
    SHOW((double)big.query(1000000000000000000LL), "4e36 + 5（double にして表示）");
    SHOW((double)big.query(-1000000000000000000LL), "-4e36 + 5");
    note("x や係数が実数なら LiChaoTreeMin<double>（xs も double の一覧で渡す）");
}

}  // namespace lichao_demo

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
    lichao_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    lichao_test::benchmark();    // 速度計測
    lichao_demo::demo();         // 使い方デモ
    return 0;
}
