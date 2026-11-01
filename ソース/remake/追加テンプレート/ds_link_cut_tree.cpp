// =============================================================================
//  Link-Cut Tree（辺の追加・削除がある森で、パスの集約・連結判定・LCA を O(log N)（ならし））   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  森（木の集まり）を「よく使うパス」ごとにスプレー木で持つ。根を好きな頂点に変えられる（evert）ので、
//  辺 u–v をつなぐ・切る、u–v のパス上の頂点の値をまとめる、u と v が同じ木か、などをすべて O(log N)（ならし）で答える
//    例: 「辺の追加・削除と、パス上の和・最大を答えるクエリ」、動的な最小全域木（辺を頂点にして最大の辺を探す）
//
//  機能一覧（頂点は 0〜N−1。値は頂点に持つ。S・op・e はモノイド（交換できなくてよい。パスは u から v の順））     計算量（ならし）
//    make_link_cut_tree(n, op, e)          全部の値が e、辺なしの森                                         O(N)
//    make_link_cut_tree(a, op, e)          頂点 i の値が a[i]
//    link(u, v)                            辺 u–v を足す（u と v は別の木であること）                         O(log N)
//    cut(u, v)                             辺 u–v を消す（その辺があること）                                 O(log N)
//    connected(u, v)                       同じ木か                                                         O(log N)
//    prod(u, v)                            パス u → v の頂点の値を順に op でまとめた値（同じ木であること）      O(log N)
//    set(v, x) / get(v)                    頂点 v の値を x に / 頂点 v の値                                  O(log N)
//    evert(v) / root(v)                    v を木の根にする / v の木の今の根                                  O(log N)
//    lca(u, v, r)                          根を r にしたときの u と v の LCA（同じ木であること）               O(log N)
//    path_length(u, v)                     パス u–v の辺の数                                                 O(log N)
//    has_edge(u, v)                        辺 u–v があるか                                                   O(log N)
//
//  注意: link でつなぐと閉路ができる場合・cut で無い辺を指定した場合は止める。辺に値を持たせたいときは辺を頂点にする
//  速さの目安: N = 2×10^5 で link を N 回、そのあと cut・link・prod を 2×10^5 回ずつ 約 0.4 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    lct_test::random_test()  … 森を隣接リストで持ち、パスを幅優先探索で探す愚直解との突き合わせ
//    lct_test::benchmark()    … 速度計測
//    lct_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace lct_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace lct_internal

template <class S, class Op>
struct LinkCutTree {
    LinkCutTree(const std::vector<S>& a, Op op, S e) : n_((int)a.size()), op_(op), e_(e), nd_(a.size() + 1) {
        nd_[0].val = nd_[0].sum = nd_[0].rsum = e_;
        for (int i = 0; i < n_; ++i) nd_[i + 1].val = nd_[i + 1].sum = nd_[i + 1].rsum = a[i];
    }
    int size() const { return n_; }
    void link(int u, int v) {
        chk(u, "link"), chk(v, "link");
        if (connected(u, v)) lct_internal::fail("LinkCutTree::link", "u と v はもう同じ木（閉路ができる）");
        ++u, ++v;
        evert_(u);
        nd_[u].p = v;
    }
    void cut(int u, int v) {
        chk(u, "cut"), chk(v, "cut");
        ++u, ++v;
        evert_(u);
        access(v);
        splay(v);
        if (nd_[v].ch[0] != u || nd_[u].ch[1] != 0) lct_internal::fail("LinkCutTree::cut", "辺 u–v が無い");
        nd_[v].ch[0] = 0, nd_[u].p = 0;
        update(v);
    }
    bool has_edge(int u, int v) {
        chk(u, "has_edge"), chk(v, "has_edge");
        if (u == v || !connected(u, v)) return false;
        ++u, ++v;
        evert_(u);
        access(v);
        splay(v);
        return nd_[v].ch[0] == u && nd_[u].ch[1] == 0;
    }
    bool connected(int u, int v) {
        chk(u, "connected"), chk(v, "connected");
        return root_(u + 1) == root_(v + 1);
    }
    S prod(int u, int v) {
        chk(u, "prod"), chk(v, "prod");
        if (!connected(u, v)) lct_internal::fail("LinkCutTree::prod", "u と v が別の木");
        ++u, ++v;
        evert_(u);
        access(v);
        splay(v);
        return nd_[v].sum;
    }
    int path_length(int u, int v) {
        chk(u, "path_length"), chk(v, "path_length");
        if (!connected(u, v)) lct_internal::fail("LinkCutTree::path_length", "u と v が別の木");
        ++u, ++v;
        evert_(u);
        access(v);
        splay(v);
        return nd_[v].sz - 1;
    }
    void set(int v, const S& x) {
        chk(v, "set");
        ++v;
        access(v);
        splay(v);
        nd_[v].val = x;
        update(v);
    }
    S get(int v) const { return nd_[v + 1].val; }
    void evert(int v) {
        chk(v, "evert");
        evert_(v + 1);
    }
    int root(int v) {
        chk(v, "root");
        return root_(v + 1) - 1;
    }
    int lca(int u, int v, int r) {
        chk(u, "lca"), chk(v, "lca"), chk(r, "lca");
        if (!connected(u, v) || !connected(u, r)) lct_internal::fail("LinkCutTree::lca", "u, v, r が同じ木にない");
        evert_(r + 1);
        access(u + 1);
        return access(v + 1) - 1;
    }

   private:
    struct Node {
        int ch[2] = {0, 0}, p = 0, sz = 1;
        bool rev = false;
        S val{}, sum{}, rsum{};  // sum = 左 → 右の順、rsum = 右 → 左の順
    };
    int n_;
    Op op_;
    S e_;
    std::vector<Node> nd_;
    void chk(int v, const char* f) const {
        if (v < 0 || v >= n_) lct_internal::fail(f, "頂点番号が範囲外");
    }
    bool is_root(int x) const {
        const int p = nd_[x].p;
        return p == 0 || (nd_[p].ch[0] != x && nd_[p].ch[1] != x);
    }
    void toggle(int x) {
        if (!x) return;
        std::swap(nd_[x].ch[0], nd_[x].ch[1]);
        std::swap(nd_[x].sum, nd_[x].rsum);
        nd_[x].rev ^= 1;
    }
    void push(int x) {
        if (nd_[x].rev) toggle(nd_[x].ch[0]), toggle(nd_[x].ch[1]), nd_[x].rev = false;
    }
    void update(int x) {
        const int l = nd_[x].ch[0], r = nd_[x].ch[1];
        nd_[x].sz = 1 + (l ? nd_[l].sz : 0) + (r ? nd_[r].sz : 0);
        nd_[x].sum = op_(op_(l ? nd_[l].sum : e_, nd_[x].val), r ? nd_[r].sum : e_);
        nd_[x].rsum = op_(op_(r ? nd_[r].rsum : e_, nd_[x].val), l ? nd_[l].rsum : e_);
    }
    void rotate(int x) {
        const int p = nd_[x].p, g = nd_[p].p;
        const int d = nd_[p].ch[1] == x;
        const int c = nd_[x].ch[d ^ 1];
        if (!is_root(p)) nd_[g].ch[nd_[g].ch[1] == p] = x;
        nd_[x].p = g;
        nd_[x].ch[d ^ 1] = p, nd_[p].p = x;
        nd_[p].ch[d] = c;
        if (c) nd_[c].p = p;
        update(p), update(x);
    }
    std::vector<int> stk_;
    void splay(int x) {
        stk_.clear();
        for (int y = x;; y = nd_[y].p) {
            stk_.push_back(y);
            if (is_root(y)) break;
        }
        for (size_t k = stk_.size(); k-- > 0;) push(stk_[k]);
        while (!is_root(x)) {
            const int p = nd_[x].p;
            if (!is_root(p)) rotate((nd_[nd_[p].p].ch[1] == p) == (nd_[p].ch[1] == x) ? p : x);
            rotate(x);
        }
    }
    int access(int x) {  // 根から x までを 1 本のパスに。最後にくっつけた頂点を返す（LCA に使う）
        int last = 0;
        for (int y = x; y; y = nd_[y].p) {
            splay(y);
            nd_[y].ch[1] = last;
            update(y);
            last = y;
        }
        splay(x);
        return last;
    }
    void evert_(int x) {
        access(x);
        toggle(x);
        push(x);
    }
    int root_(int x) {
        access(x);
        while (true) {
            push(x);
            if (!nd_[x].ch[0]) break;
            x = nd_[x].ch[0];
        }
        splay(x);
        return x;
    }
};
template <class S, class Op>
LinkCutTree<S, Op> make_link_cut_tree(int n, Op op, S e) {
    return LinkCutTree<S, Op>(std::vector<S>(n, e), op, e);
}
template <class S, class Op>
LinkCutTree<S, Op> make_link_cut_tree(const std::vector<S>& a, Op op, S e) {
    return LinkCutTree<S, Op>(a, op, e);
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace lct_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define LC_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_link_cut_tree] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
using F = std::pair<ll, ll>;
static const ll MOD = 998244353;
static F comp(const F& f, const F& g) { return {f.first * g.first % MOD, (f.second * g.first + g.second) % MOD}; }

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 200 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 15);
        std::vector<F> val(n);
        for (auto& f : val) f = {rnd(rng, 1, 9), rnd(rng, 0, 9)};
        auto lct = make_link_cut_tree(val, comp, F{1, 0});
        std::vector<std::vector<int>> adj(n);
        auto path = [&](int s, int t) {  // s から t への頂点の列（無ければ空）
            std::vector<int> par(n, -2);
            std::queue<int> q;
            q.push(s), par[s] = -1;
            while (!q.empty()) {
                const int x = q.front();
                q.pop();
                for (int y : adj[x])
                    if (par[y] == -2) par[y] = x, q.push(y);
            }
            std::vector<int> p;
            if (par[t] == -2) return p;
            for (int x = t; x != -1; x = par[x]) p.push_back(x);
            std::reverse(p.begin(), p.end());
            return p;
        };
        for (int t = 0; t < 120; ++t) {
            const int u = (int)rnd(rng, 0, n - 1), v = (int)rnd(rng, 0, n - 1);
            const int k = (int)rnd(rng, 0, 7);
            const auto pth = path(u, v);
            const bool conn = !pth.empty();
            LC_CHECK(lct.connected(u, v) == conn);
            if (k <= 1 && !conn) {
                lct.link(u, v);
                adj[u].push_back(v), adj[v].push_back(u);
            } else if (k == 2) {
                // 辺をランダムに 1 本消す
                std::vector<std::pair<int, int>> es;
                for (int a = 0; a < n; ++a)
                    for (int b : adj[a])
                        if (a < b) es.push_back({a, b});
                if (!es.empty()) {
                    auto [a, b] = es[rng() % es.size()];
                    LC_CHECK(lct.has_edge(a, b) && lct.has_edge(b, a));
                    if (rnd(rng, 0, 1)) lct.cut(a, b);
                    else lct.cut(b, a);
                    adj[a].erase(std::find(adj[a].begin(), adj[a].end(), b));
                    adj[b].erase(std::find(adj[b].begin(), adj[b].end(), a));
                }
            } else if (k == 3) {
                const F x{rnd(rng, 1, 9), rnd(rng, 0, 9)};
                lct.set(u, x), val[u] = x;
                LC_CHECK(lct.get(u) == x);
            } else if (conn) {
                F s{1, 0};
                for (int x : pth) s = comp(s, val[x]);
                LC_CHECK(lct.prod(u, v) == s);
                LC_CHECK(lct.path_length(u, v) == (int)pth.size() - 1);
                LC_CHECK(lct.has_edge(u, v) == (pth.size() == 2));
                // LCA：根 r からのパスで、u と v の共通の最後の頂点
                const int r = pth[rng() % pth.size()];  // 同じ木の頂点
                const auto pu = path(r, u), pv = path(r, v);
                int w = r;
                for (size_t i = 0; i < std::min(pu.size(), pv.size()) && pu[i] == pv[i]; ++i) w = pu[i];
                LC_CHECK(lct.lca(u, v, r) == w);
                LC_CHECK(lct.root(u) == r && lct.root(v) == r);  // lca で r が根になっている
            }
            ++cases;
        }
    }
    std::cout << "[ds_link_cut_tree] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_link_cut_tree] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    std::vector<ll> a(N);
    for (auto& x : a) x = rnd(rng, 0, 1000000000);
    const auto t0 = std::chrono::steady_clock::now();
    auto lct = make_link_cut_tree(a, [](ll x, ll y) { return x + y; }, 0LL);
    std::vector<int> par(N, -1);
    ll chk = 0;
    for (int i = 1; i < N; ++i) par[i] = (int)rnd(rng, 0, i - 1), lct.link(i, par[i]);
    for (int i = 0; i < N; ++i) {
        const int v = (int)rnd(rng, 1, N - 1);
        // 辺をつなぎ替える（v を切って、v の部分木に入らない頂点 0 につなぐ）
        lct.cut(v, par[v]);
        lct.link(v, 0);
        par[v] = 0;
        chk += lct.prod((int)rnd(rng, 0, N - 1), (int)rnd(rng, 0, N - 1));
    }
    std::cout << "[ds_link_cut_tree] benchmark (N = " << N << "): link " << N - 1 << " 回の後、cut・link・prod を " << N << " 回ずつ " << elapsed_ms(t0)
              << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace lct_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace lct_demo {

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
#define SHOW(expr, comment) lct_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    lct_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("Link-Cut Tree（ds_link_cut_tree） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  make_link_cut_tree(n or a, op, e): link / cut / connected / has_edge / prod / path_length / set / get / evert / root / lca\n";

    api("【1】 make_link_cut_tree(a, op, e) / link(u, v) / prod(u, v)", "— / void / S", "辺をつなぎながら、パスの頂点の値をまとめる", "O(log N)（ならし）");
    RUN(auto lct = make_link_cut_tree(std::vector<long long>{1, 2, 3, 4, 5}, [](long long x, long long y) { return x + y; }, 0LL));
    RUN(lct.link(0, 1); lct.link(1, 2); lct.link(1, 3));  // 0-1, 1-2, 1-3
    SHOW(lct.prod(2, 3), "2 → 1 → 3 で 3 + 2 + 4");
    SHOW(lct.connected(0, 4), "4 はまだつながっていない");
    RUN(lct.link(4, 3));
    SHOW(lct.prod(0, 4), "0 → 1 → 3 → 4");
    SHOW(lct.path_length(0, 4), "");

    api("【2】 cut(u, v) / has_edge(u, v)", "void / bool", "辺を消す / 辺があるか", "O(log N)");
    RUN(lct.cut(1, 3));
    SHOW(lct.connected(0, 4), "");
    SHOW(lct.has_edge(3, 4), "");

    api("【3】 set(v, x) / lca(u, v, r) / root(v)", "void / int / int", "頂点の値を変える / 根を r にしたときの LCA / 今の根", "O(log N)");
    RUN(lct.set(1, 100));
    SHOW(lct.prod(0, 2), "1 + 100 + 3");
    SHOW(lct.lca(0, 2, 2), "根 2 で 0 と 2 の LCA は 2");
    SHOW(lct.lca(0, 2, 1), "根 1 なら 1");
}

}  // namespace lct_demo

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
    lct_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    lct_test::benchmark();    // 速度計測
    lct_demo::demo();         // 使い方デモ
    return 0;
}
