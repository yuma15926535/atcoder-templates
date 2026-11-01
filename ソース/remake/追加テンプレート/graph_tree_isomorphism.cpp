// =============================================================================
//  木の同型判定（根付き木・根なし木が同じ形か、部分木ごとの形の番号）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  根付き木の形に番号をつける（AHU）：葉から順に「子の番号を並べ替えた列」を辞書に入れ、同じ列なら同じ番号。
//  ハッシュではないので衝突しない。根なし木は中心（直径の真ん中。1 つか 2 つ）を根にした番号の小さい方を使う
//    例: 「2 つの木は頂点の番号を付け替えて同じになるか」「部分木の形の種類数」「同じ形の部分木の組の数」
//
//  機能一覧（g は木の隣接リスト vector<vector<int>>。頂点は 0〜N−1）                                計算量
//    TreeCanon tc                        形の辞書（いくつもの木で番号を共通にする）
//    tc.subtree_ids(g, root)             各頂点 v の部分木の形の番号（同じ番号 ⇔ 同じ形）                     O(N log N)
//    tc.rooted(g, root)                  根付き木の形の番号                                                 O(N log N)
//    tc.unrooted(g)                      根なし木の形の番号（中心を根にした番号の小さい方）                     O(N log N)
//    tree_centers(g)                     木の中心（1 つか 2 つ）                                            O(N)
//    is_isomorphic(g1, g2)               根なし木として同じ形か                                             O(N log N)
//    is_isomorphic_rooted(g1, r1, g2, r2)   根付き木として同じ形か
//
//  注意: 再帰を使わない。違う TreeCanon で作った番号どうしは比べられない
//  速さの目安: N = 2×10^5 で unrooted 約 0.06 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    treeiso_test::random_test()  … 頂点の付け替えを全部試す愚直解・括弧の文字列との突き合わせ
//    treeiso_test::benchmark()    … 速度計測
//    treeiso_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace treeiso_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
inline void order_from(const std::vector<std::vector<int>>& g, int root, std::vector<int>& ord, std::vector<int>& par, const char* f) {
    const int n = (int)g.size();
    if (root < 0 || root >= n) fail(f, "root が範囲外");
    ord.clear(), par.assign(n, -2);
    std::vector<int> st{root};
    par[root] = -1;
    while (!st.empty()) {
        const int v = st.back();
        st.pop_back();
        ord.push_back(v);
        for (int u : g[v]) {
            if (u == par[v]) continue;
            if (par[u] != -2) fail(f, "g は木にする（閉路がある）");
            par[u] = v, st.push_back(u);
        }
    }
    if ((int)ord.size() != n) fail(f, "g は連結な木にする");
}
}  // namespace treeiso_internal

inline std::vector<int> tree_centers(const std::vector<std::vector<int>>& g) {
    const int n = (int)g.size();
    if (n == 0) return {};
    if (n == 1) return {0};
    // 葉を外側から削る
    std::vector<int> deg(n), cur;
    for (int v = 0; v < n; ++v) {
        deg[v] = (int)g[v].size();
        if (deg[v] <= 1) cur.push_back(v);
    }
    int remain = n;
    while (remain > 2) {
        std::vector<int> nxt;
        remain -= (int)cur.size();
        for (int v : cur)
            for (int u : g[v])
                if (--deg[u] == 1) nxt.push_back(u);
        cur.swap(nxt);
    }
    std::sort(cur.begin(), cur.end());
    return cur;
}

struct TreeCanon {
    std::vector<int> subtree_ids(const std::vector<std::vector<int>>& g, int root) {
        std::vector<int> ord, par;
        treeiso_internal::order_from(g, root, ord, par, "TreeCanon::subtree_ids");
        const int n = (int)g.size();
        std::vector<int> id(n);
        std::vector<std::vector<int>> ch(n);
        for (int i = n - 1; i >= 0; --i) {
            const int v = ord[i];
            auto& c = ch[v];
            std::sort(c.begin(), c.end());
            auto it = dict_.find(c);
            if (it == dict_.end()) it = dict_.emplace(c, (int)dict_.size()).first;
            id[v] = it->second;
            if (par[v] >= 0) ch[par[v]].push_back(id[v]);
            std::vector<int>().swap(c);
        }
        return id;
    }
    int rooted(const std::vector<std::vector<int>>& g, int root) { return subtree_ids(g, root)[root]; }
    int unrooted(const std::vector<std::vector<int>>& g) {
        if (g.empty()) return -1;
        const auto c = tree_centers(g);
        int best = rooted(g, c[0]);
        if (c.size() == 2) best = std::min(best, rooted(g, c[1]));
        return best;
    }
    int kinds() const { return (int)dict_.size(); }

   private:
    std::map<std::vector<int>, int> dict_;
};

inline bool is_isomorphic(const std::vector<std::vector<int>>& g1, const std::vector<std::vector<int>>& g2) {
    if (g1.size() != g2.size()) return false;
    TreeCanon tc;
    return tc.unrooted(g1) == tc.unrooted(g2);
}
inline bool is_isomorphic_rooted(const std::vector<std::vector<int>>& g1, int r1, const std::vector<std::vector<int>>& g2, int r2) {
    if (g1.size() != g2.size()) return false;
    TreeCanon tc;
    return tc.rooted(g1, r1) == tc.rooted(g2, r2);
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace treeiso_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define TI_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[graph_tree_isomorphism] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
using G = std::vector<std::vector<int>>;

static G random_tree(std::mt19937_64& rng, int n) {
    G g(n);
    for (int i = 1; i < n; ++i) {
        const int p = (int)rnd(rng, 0, i - 1);
        g[i].push_back(p), g[p].push_back(i);
    }
    // 番号を混ぜる
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    G h(n);
    for (int v = 0; v < n; ++v)
        for (int u : g[v]) h[perm[v]].push_back(perm[u]);
    return h;
}
// 括弧の文字列（AHU の文字列版）
static std::string canon_str(const G& g, int v, int p) {
    std::vector<std::string> cs;
    for (int u : g[v])
        if (u != p) cs.push_back(canon_str(g, u, v));
    std::sort(cs.begin(), cs.end());
    std::string s = "(";
    for (auto& c : cs) s += c;
    return s + ")";
}
static bool brute_iso(const G& a, const G& b) {
    const int n = (int)a.size();
    if ((int)b.size() != n) return false;
    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 0);
    std::set<std::pair<int, int>> eb;
    for (int v = 0; v < n; ++v)
        for (int u : b[v]) eb.insert({v, u});
    do {
        bool ok = true;
        for (int v = 0; v < n && ok; ++v)
            for (int u : a[v])
                if (!eb.count({p[v], p[u]})) ok = false;
        if (ok) return true;
    } while (std::next_permutation(p.begin(), p.end()));
    return false;
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 7);
        const G a = random_tree(rng, n), b = rnd(rng, 0, 1) ? random_tree(rng, n) : [&] {
            // a の番号を付け替えたもの（必ず同じ形）
            std::vector<int> perm(n);
            std::iota(perm.begin(), perm.end(), 0);
            std::shuffle(perm.begin(), perm.end(), rng);
            G h(n);
            for (int v = 0; v < n; ++v)
                for (int u : a[v]) h[perm[v]].push_back(perm[u]);
            return h;
        }();
        TI_CHECK(is_isomorphic(a, b) == brute_iso(a, b));
        const int r1 = (int)rnd(rng, 0, n - 1), r2 = (int)rnd(rng, 0, n - 1);
        TI_CHECK(is_isomorphic_rooted(a, r1, b, r2) == (canon_str(a, r1, -1) == canon_str(b, r2, -1)));
        // 部分木の番号：同じ番号 ⇔ 同じ括弧の文字列
        TreeCanon tc;
        const auto ia = tc.subtree_ids(a, r1), ib = tc.subtree_ids(b, r2);
        std::vector<int> pa(n, -1), pb(n, -1);
        {
            std::vector<int> ord;
            treeiso_internal::order_from(a, r1, ord, pa, "t");
            treeiso_internal::order_from(b, r2, ord, pb, "t");
        }
        std::vector<std::string> sa(n), sb(n);
        for (int v = 0; v < n; ++v) sa[v] = canon_str(a, v, pa[v]), sb[v] = canon_str(b, v, pb[v]);
        for (int x = 0; x < n; ++x)
            for (int y = 0; y < n; ++y) TI_CHECK((ia[x] == ib[y]) == (sa[x] == sb[y]) && (ia[x] == ia[y]) == (sa[x] == sa[y]));
        // 中心：直径の真ん中
        const auto c = tree_centers(a);
        auto ecc = [&](int s) {
            std::vector<int> d(n, -1), q{s};
            d[s] = 0;
            for (size_t k = 0; k < q.size(); ++k)
                for (int u : a[q[k]])
                    if (d[u] < 0) d[u] = d[q[k]] + 1, q.push_back(u);
            return *std::max_element(d.begin(), d.end());
        };
        int mn = n;
        for (int v = 0; v < n; ++v) mn = std::min(mn, ecc(v));
        std::vector<int> ec;
        for (int v = 0; v < n; ++v)
            if (ecc(v) == mn) ec.push_back(v);
        TI_CHECK(c == ec);
        cases += 4;
    }
    std::cout << "[graph_tree_isomorphism] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[graph_tree_isomorphism] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    const G a = random_tree(rng, N);
    const auto t0 = std::chrono::steady_clock::now();
    TreeCanon tc;
    const int x = tc.unrooted(a);
    std::cout << "[graph_tree_isomorphism] benchmark: unrooted（N = " << N << "） " << elapsed_ms(t0) << " ms  (check " << x % 1000 << ", 形の種類 " << tc.kinds()
              << ")\n";
}

}  // namespace treeiso_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace treeiso_demo {

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
#define SHOW(expr, comment) treeiso_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    treeiso_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("木の同型判定（graph_tree_isomorphism） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  TreeCanon tc: subtree_ids(g, root) / rooted(g, root) / unrooted(g) / kinds()\n"
              << "  tree_centers(g) / is_isomorphic(g1, g2) / is_isomorphic_rooted(g1, r1, g2, r2)\n";

    RUN(const std::vector<std::vector<int>> a = {{1, 2}, {0, 3}, {0}, {1}});      // 道 3-1-0-2
    RUN(const std::vector<std::vector<int>> b = {{3}, {2, 3}, {1}, {0, 1}});      // 道 0-3-1-2
    RUN(const std::vector<std::vector<int>> s = {{1, 2, 3}, {0}, {0}, {0}});      // 星

    api("【1】 is_isomorphic(g1, g2) / is_isomorphic_rooted(g1, r1, g2, r2)", "bool", "根なし木として / 根付き木として同じ形か", "O(N log N)");
    SHOW(is_isomorphic(a, b), "どちらも長さ 4 の道");
    SHOW(is_isomorphic(a, s), "道と星");
    SHOW(is_isomorphic_rooted(a, 0, b, 3), "どちらも道の 2 番目を根に");
    SHOW(is_isomorphic_rooted(a, 0, b, 0), "道の端を根にしたものとは違う");

    api("【2】 TreeCanon tc / subtree_ids(g, root)", "vector<int>", "各頂点の部分木の形の番号（同じ番号 ⇔ 同じ形）", "O(N log N)");
    RUN(TreeCanon tc);
    SHOW(tc.subtree_ids(s, 0), "葉はみな同じ番号");
    SHOW(tc.subtree_ids(a, 0), "同じ tc なので番号を比べられる");
    SHOW(tc.kinds(), "今までに出た形の種類");

    api("【3】 tree_centers(g) / unrooted(g)", "vector<int> / int", "木の中心（1 つか 2 つ）/ 根なし木の形の番号", "O(N) / O(N log N)");
    SHOW(tree_centers(a), "長さ 4 の道の真ん中 2 つ");
    SHOW(tc.unrooted(a) == tc.unrooted(b), "");
}

}  // namespace treeiso_demo

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
    treeiso_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    treeiso_test::benchmark();    // 速度計測
    treeiso_demo::demo();         // 使い方デモ
    return 0;
}
