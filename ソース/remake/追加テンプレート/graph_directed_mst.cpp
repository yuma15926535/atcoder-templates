// =============================================================================
//  有向グラフの最小全域木（最小有向全域木 / Chu-Liu・Edmonds：根から全頂点へ向かう最小コストの木）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  各頂点（根以外）に入ってくる一番安い辺を選び、閉路ができたら閉路を 1 つの頂点に縮めて（入る辺のコストから
//  閉路の中の辺のコストを引いて）続ける。ヒープの併合（左偏木）と巻き戻せる Union-Find で O(M log N)。使った辺も復元する
//    例: 「根から全部の頂点へ情報を届ける一方通行の通信路を最小コストで選ぶ」、有向グラフで「全員にたどり着ける」最小の辺の組
//
//  機能一覧（頂点は 0〜N−1。辺は {from, to, cost}（重みは負でもよい））                              計算量
//    directed_mst(n, root, edges)        DirectedMSTResult {ok, cost, parent, edge_id}                        O(M log N)
//                                          ok        … 全頂点に根から届くか（届かなければ false）
//                                          cost      … 最小コスト
//                                          parent[v] … v の親（根は −1）、edge_id[v] … v に入る辺の edges での番号（根は −1）
//
//  注意: 自己ループの辺は使わない。多重辺はあってよい。無向グラフの最小全域木なら 03_kruskal.cpp
//  速さの目安: N = 10^5、M = 2×10^5 で 約 0.04 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    dmst_test::random_test()  … 各頂点に入る辺の選び方を全部試す愚直解との突き合わせ
//    dmst_test::benchmark()    … 速度計測
//    dmst_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace dmst_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
struct E {
    int a, b, id;
    long long w;
};
// 左偏木（遅延で全体から引く）
struct Heap {
    struct N {
        E key;
        long long delta;
        int l, r, rk;
    };
    std::vector<N> nd;
    void prop(int x) {
        if (nd[x].delta == 0) return;
        nd[x].key.w += nd[x].delta;
        if (nd[x].l >= 0) nd[nd[x].l].delta += nd[x].delta;
        if (nd[x].r >= 0) nd[nd[x].r].delta += nd[x].delta;
        nd[x].delta = 0;
    }
    int rank(int x) const { return x < 0 ? 0 : nd[x].rk; }
    int merge(int a, int b) {  // 右の背骨の長さは O(log) なので再帰でよい
        if (a < 0 || b < 0) return a < 0 ? b : a;
        prop(a), prop(b);
        if (nd[a].key.w > nd[b].key.w) std::swap(a, b);
        const int r = merge(nd[a].r, b);
        nd[a].r = r;
        if (rank(nd[a].l) < rank(nd[a].r)) std::swap(nd[a].l, nd[a].r);
        nd[a].rk = rank(nd[a].r) + 1;
        return a;
    }
    int make(const E& e) {
        nd.push_back({e, 0, -1, -1, 1});
        return (int)nd.size() - 1;
    }
    E top(int x) {
        prop(x);
        return nd[x].key;
    }
    int pop(int x) {
        prop(x);
        return merge(nd[x].l, nd[x].r);
    }
};
// 巻き戻せる Union-Find
struct RUF {
    std::vector<int> e;
    std::vector<std::pair<int, int>> st;
    explicit RUF(int n) : e(n, -1) {}
    int find(int x) const { return e[x] < 0 ? x : find(e[x]); }
    int time() const { return (int)st.size(); }
    void rollback(int t) {
        for (int i = time(); i-- > t;) e[st[i].first] = st[i].second;
        st.resize(t);
    }
    bool join(int a, int b) {
        a = find(a), b = find(b);
        if (a == b) return false;
        if (e[a] > e[b]) std::swap(a, b);
        st.push_back({a, e[a]}), st.push_back({b, e[b]});
        e[a] += e[b], e[b] = a;
        return true;
    }
};
}  // namespace dmst_internal

struct DirectedMSTResult {
    bool ok = false;
    long long cost = 0;
    std::vector<int> parent, edge_id;
};

inline DirectedMSTResult directed_mst(int n, int root, const std::vector<std::array<long long, 3>>& edges) {
    using namespace dmst_internal;
    if (root < 0 || root >= n) fail("directed_mst", "root が範囲外");
    Heap H;
    RUF uf(n);
    std::vector<int> heap(n, -1);
    for (int i = 0; i < (int)edges.size(); ++i) {
        const int a = (int)edges[i][0], b = (int)edges[i][1];
        if (a < 0 || a >= n || b < 0 || b >= n) fail("directed_mst", "頂点番号が範囲外");
        if (a == b) continue;
        heap[b] = H.merge(heap[b], H.make({a, b, i, edges[i][2]}));
    }
    DirectedMSTResult res;
    std::vector<int> seen(n, -1), path(n);
    seen[root] = root;
    std::vector<E> Q(n), in(n, {-1, -1, -1, 0});
    std::deque<std::tuple<int, int, std::vector<E>>> cycs;
    long long cost = 0;
    for (int s = 0; s < n; ++s) {
        int u = s, qi = 0, w;
        while (seen[u] < 0) {
            if (heap[u] < 0) return res;  // 届かない
            const E e = H.top(heap[u]);
            H.nd[heap[u]].delta -= e.w;
            heap[u] = H.pop(heap[u]);
            Q[qi] = e, path[qi++] = u, seen[u] = s;
            cost += e.w, u = uf.find(e.a);
            if (seen[u] == s) {  // 閉路を縮める
                int cyc = -1;
                const int end = qi, t = uf.time();
                do cyc = H.merge(cyc, heap[w = path[--qi]]);
                while (uf.join(u, w));
                u = uf.find(u), heap[u] = cyc, seen[u] = -1;
                cycs.push_front({u, t, std::vector<E>(Q.begin() + qi, Q.begin() + end)});
            }
        }
        for (int i = 0; i < qi; ++i) in[uf.find(Q[i].b)] = Q[i];
    }
    for (auto& [u, t, comp] : cycs) {  // 縮めた閉路を戻して、使う辺を決める
        uf.rollback(t);
        const E in_edge = in[u];
        for (const E& e : comp) in[uf.find(e.b)] = e;
        in[uf.find(in_edge.b)] = in_edge;
    }
    res.ok = true, res.cost = 0;
    res.parent.assign(n, -1), res.edge_id.assign(n, -1);
    for (int v = 0; v < n; ++v)
        if (v != root) {
            res.parent[v] = in[v].a, res.edge_id[v] = in[v].id;
            res.cost += edges[in[v].id][2];
        }
    (void)cost;
    return res;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace dmst_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define DM_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_directed_mst] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 6), m = (int)rnd(rng, 0, 12);
        std::vector<std::array<ll, 3>> es(m);
        for (auto& e : es) e = {rnd(rng, 0, n - 1), rnd(rng, 0, n - 1), rnd(rng, -5, 20)};
        const int root = (int)rnd(rng, 0, n - 1);
        const auto r = directed_mst(n, root, es);
        // 愚直：根以外の各頂点に入る辺を 1 本ずつ選ぶ組み合わせを全部
        std::vector<std::vector<int>> inc(n);
        for (int i = 0; i < m; ++i)
            if (es[i][0] != es[i][1]) inc[es[i][1]].push_back(i);
        bool found = false;
        ll best = 0;
        std::vector<int> pick(n, -1);
        auto rec = [&](auto&& self, int v) -> void {
            if (v == n) {
                // 親をたどると根に着くか
                for (int x = 0; x < n; ++x) {
                    int y = x, steps = 0;
                    while (y != root && steps <= n) y = (int)es[pick[y]][0], ++steps;
                    if (y != root) return;
                }
                ll c = 0;
                for (int x = 0; x < n; ++x)
                    if (x != root) c += es[pick[x]][2];
                if (!found || c < best) best = c, found = true;
                return;
            }
            if (v == root) return self(self, v + 1);
            for (int i : inc[v]) pick[v] = i, self(self, v + 1);
        };
        rec(rec, 0);
        DM_CHECK(r.ok == found);
        if (found) {
            DM_CHECK(r.cost == best);
            // 返した辺が本当に木で、そのコスト
            ll c = 0;
            for (int x = 0; x < n; ++x) {
                if (x == root) {
                    DM_CHECK(r.parent[x] == -1 && r.edge_id[x] == -1);
                    continue;
                }
                const auto& e = es[r.edge_id[x]];
                DM_CHECK(e[1] == x && e[0] == r.parent[x]);
                c += e[2];
                int y = x, steps = 0;
                while (y != root && steps <= n) y = r.parent[y], ++steps;
                DM_CHECK(y == root);
            }
            DM_CHECK(c == best);
        }
        ++cases;
    }
    std::cout << "[graph_directed_mst] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_directed_mst] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 10000 : 100000, M = 2 * N;
    std::vector<std::array<ll, 3>> es;
    for (int i = 1; i < N; ++i) es.push_back({rnd(rng, 0, i - 1), i, rnd(rng, 1, 1000000000)});  // 必ず届く
    while ((int)es.size() < M) es.push_back({rnd(rng, 0, N - 1), rnd(rng, 0, N - 1), rnd(rng, 1, 1000000000)});
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = directed_mst(N, 0, es);
    std::cout << "[graph_directed_mst] benchmark (N = " << N << ", M = " << M << "): " << elapsed_ms(t0) << " ms  (check " << r.cost % 1000 << ")\n";
}

}  // namespace dmst_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace dmst_demo {

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
#define SHOW(expr, comment) dmst_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    dmst_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("最小有向全域木（graph_directed_mst） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  directed_mst(n, root, edges{from, to, cost}) → DirectedMSTResult {ok, cost, parent, edge_id}\n";

    api("【1】 directed_mst(n, root, edges)", "DirectedMSTResult", "根から全頂点へ向かう木で、辺のコストの和が最小のもの", "O(M log N)");
    RUN(const std::vector<std::array<long long, 3>> es = {{0, 1, 10}, {0, 2, 3}, {2, 1, 4}, {1, 3, 2}, {3, 2, 1}, {2, 3, 8}});
    RUN(const auto r = directed_mst(4, 0, es));
    SHOW(r.ok, "");
    SHOW(r.cost, "0 → 2（3）、2 → 1（4）、1 → 3（2）");
    SHOW(r.parent, "各頂点の親（根は −1）");
    SHOW(r.edge_id, "使った辺の番号");
    SHOW(directed_mst(3, 0, {{1, 2, 1}}).ok, "根から届かない頂点がある");
}

}  // namespace dmst_demo

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
    dmst_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    dmst_test::benchmark();    // 速度計測
    dmst_demo::demo();         // 使い方デモ
    return 0;
}
