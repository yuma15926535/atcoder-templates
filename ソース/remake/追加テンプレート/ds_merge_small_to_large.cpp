// =============================================================================
//  マージテク（小さい方を大きい方に混ぜる Union-Find・DSU on tree）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  集合を併合するとき、要素の少ない方を多い方へ 1 つずつ移すと、1 つの要素が移るたびに入っている集合の大きさが 2 倍以上になる
//  ので、各要素が移るのは高々 log2 N 回。全体で O(N log N) 回の移動（set なら O(N log^2 N)）。
//  木の上では「重い子の情報はそのまま使い、軽い子の部分木だけ足し直す」（DSU on tree / sack）
//    例: ABC183 F（Confluence）… 成分ごとに map<クラス, 人数> を持って併合する
//        「部分木ごとの色の種類数」「部分木で一番多い色」を全部の頂点について
//
//  機能一覧                                                                                  計算量
//    SmallToLargeUF<C> uf(init)          頂点 i の成分の入れ物を init[i]（set・map・vector など）で始める     O(N)
//      uf.unite(x, y)                    併合。入れ物は小さい方を大きい方へ big.insert(small.begin(), small.end())   全体 O(N log N) 回の挿入
//      uf.unite(x, y, merge)             merge(big, small) で混ぜ方を決める（vector・自分で作った型・個数の足し算など）
//      uf.find(x) / uf.same(x, y) / uf.size(x)     代表 / 同じ成分か / 成分の頂点の数                     ならし O(α(N))
//      uf.data(x)                        x の成分の入れ物（参照）                                          O(α(N))
//    dsu_on_tree(g, root, add, remove, query)
//        各頂点 v について、v の部分木の頂点を全部 add した状態で query(v) を呼ぶ（それ以外の頂点は入っていない）。
//        add(u) / remove(u) は全体で O(N log N) 回。再帰を使わない                          O(N log N × add)
//
//  注意: merge の中で small から big へ移す（small は使い終わったら空にしてよい）。unite の返り値は併合したら true
//  速さの目安: N = 2×10^5 で SmallToLargeUF<set<int>> をランダムに全部併合 約 0.03 秒、dsu_on_tree（色の種類数）約 0.07 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    stl_test::random_test()  … 集合を毎回作り直す・部分木を全部見る愚直解との突き合わせ
//    stl_test::benchmark()    … 速度計測
//    stl_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
namespace stl_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace stl_internal

template <class C>
struct SmallToLargeUF {
    explicit SmallToLargeUF(std::vector<C> init) : par_(init.size()), sz_(init.size(), 1), data_(std::move(init)) {
        std::iota(par_.begin(), par_.end(), 0);
    }
    int find(int x) {
        check(x, "SmallToLargeUF::find");
        while (par_[x] != x) x = par_[x] = par_[par_[x]];
        return x;
    }
    bool same(int x, int y) { return find(x) == find(y); }
    int size(int x) { return sz_[find(x)]; }
    C& data(int x) { return data_[find(x)]; }
    bool unite(int x, int y) {
        return unite(x, y, [](C& big, C& small) { big.insert(small.begin(), small.end()); });
    }
    template <class Merge>
    bool unite(int x, int y, Merge merge) {
        x = find(x), y = find(y);
        if (x == y) return false;
        if (sz_[x] < sz_[y]) std::swap(x, y);  // 頂点の多い方を代表に
        par_[y] = x, sz_[x] += sz_[y];
        if (data_[x].size() < data_[y].size()) std::swap(data_[x], data_[y]);  // 入れ物は中身の多い方へ
        merge(data_[x], data_[y]);
        C().swap(data_[y]);
        return true;
    }

   private:
    std::vector<int> par_, sz_;
    std::vector<C> data_;
    void check(int x, const char* f) const {
        if (x < 0 || x >= (int)par_.size()) stl_internal::fail(f, "頂点番号が範囲外");
    }
};

template <class Add, class Remove, class Query>
void dsu_on_tree(const std::vector<std::vector<int>>& g, int root, Add add, Remove remove, Query query) {
    const int n = (int)g.size();
    if (n == 0) return;
    if (root < 0 || root >= n) stl_internal::fail("dsu_on_tree", "root が範囲外");
    // 行きがけ順・親・部分木の大きさ・重い子・オイラーツアーの順番
    std::vector<int> par(n, -1), ord, sz(n, 1), heavy(n, -1), tin(n), tout(n), ver(n);
    {
        std::vector<int> st{root};
        par[root] = root;
        while (!st.empty()) {
            const int v = st.back();
            st.pop_back();
            ord.push_back(v);
            for (int c : g[v])
                if (par[c] == -1) par[c] = v, st.push_back(c);
                else if (c != par[v]) stl_internal::fail("dsu_on_tree", "g は木にする");
        }
        if ((int)ord.size() != n) stl_internal::fail("dsu_on_tree", "g は連結な木にする");
        for (int i = n - 1; i > 0; --i) sz[par[ord[i]]] += sz[ord[i]];
        for (int v : ord)
            for (int c : g[v])
                if (c != par[v] && (heavy[v] < 0 || sz[c] > sz[heavy[v]])) heavy[v] = c;
        // 重い子を先に並べたオイラーツアー（部分木 = ver[tin, tout)）
        int t = 0;
        std::vector<std::pair<int, int>> s2{{root, 0}};
        while (!s2.empty()) {
            auto& [v, k] = s2.back();
            if (k == 0) tin[v] = t, ver[t++] = v;
            // 0 番目に重い子、そのあと軽い子
            int next = -1;
            if (k == 0 && heavy[v] >= 0) next = heavy[v], k = 1;
            else {
                if (k == 0) k = 1;
                while (k - 1 < (int)g[v].size()) {
                    const int c = g[v][k - 1];
                    ++k;
                    if (c != par[v] && c != heavy[v]) {
                        next = c;
                        break;
                    }
                }
            }
            if (next >= 0) s2.push_back({next, 0});
            else tout[v] = t, s2.pop_back();
        }
    }
    // 本体（再帰の代わりにフレームを積む）。phase 0: 軽い子を順に（残さない）、1: 重い子（残す）、2: まとめ
    struct Frame {
        int v;
        bool keep;
        int phase, idx;
    };
    std::vector<Frame> st{{root, true, 0, 0}};
    while (!st.empty()) {
        Frame& f = st.back();
        const int v = f.v;
        if (f.phase == 0) {
            bool pushed = false;
            while (f.idx < (int)g[v].size()) {
                const int c = g[v][f.idx++];
                if (c != par[v] && c != heavy[v]) {
                    st.push_back({c, false, 0, 0});
                    pushed = true;
                    break;
                }
            }
            if (pushed) continue;
            f.phase = 1;
            if (heavy[v] >= 0) {
                st.push_back({heavy[v], true, 0, 0});
                continue;
            }
        }
        if (f.phase == 1) {
            // 軽い子の部分木と v を足す
            for (int c : g[v])
                if (c != par[v] && c != heavy[v])
                    for (int i = tin[c]; i < tout[c]; ++i) add(ver[i]);
            add(v);
            query(v);
            if (!f.keep)
                for (int i = tin[v]; i < tout[v]; ++i) remove(ver[i]);
            st.pop_back();
        }
    }
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace stl_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define ST_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_merge_small_to_large] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static std::vector<std::vector<int>> random_tree(std::mt19937_64& rng, int n, int shape) {
    std::vector<std::vector<int>> g(n);
    for (int i = 1; i < n; ++i) {
        const int p = shape == 1 ? i - 1 : (shape == 2 ? 0 : (int)rnd(rng, 0, i - 1));
        g[i].push_back(p), g[p].push_back(i);
    }
    return g;
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 30);
        std::vector<std::set<int>> init(n);
        std::vector<std::map<int, int>> cnt(n);
        std::vector<int> cls(n);
        for (int i = 0; i < n; ++i) init[i] = {(int)rnd(rng, 0, 9)}, cls[i] = (int)rnd(rng, 0, 4), cnt[i][cls[i]] = 1;
        SmallToLargeUF<std::set<int>> uf(init);
        SmallToLargeUF<std::map<int, int>> uc(cnt);
        std::vector<int> comp(n);
        std::iota(comp.begin(), comp.end(), 0);
        for (int t = 0; t < 60; ++t) {
            const int x = (int)rnd(rng, 0, n - 1), y = (int)rnd(rng, 0, n - 1);
            const bool merged = comp[x] != comp[y];
            ST_CHECK(uf.unite(x, y) == merged);
            uc.unite(x, y, [](std::map<int, int>& big, std::map<int, int>& small) {
                for (auto [k, c] : small) big[k] += c;
            });
            if (merged) {
                const int cx = comp[x], cy = comp[y];
                for (auto& c : comp)
                    if (c == cy) c = cx;
            }
            const int q = (int)rnd(rng, 0, n - 1);
            std::set<int> es;
            std::map<int, int> ec;
            int sz = 0;
            for (int i = 0; i < n; ++i)
                if (comp[i] == comp[q]) es.insert(*init[i].begin()), ec[cls[i]]++, ++sz;
            ST_CHECK(uf.data(q) == es && uc.data(q) == ec && uf.size(q) == sz);
            const int r = (int)rnd(rng, 0, n - 1);
            ST_CHECK(uf.same(q, r) == (comp[q] == comp[r]));
            cases += 3;
        }
    }
    // dsu_on_tree：部分木の色の種類数・一番多い色の数
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 40);
        const auto g = random_tree(rng, n, (int)rnd(rng, 0, 2));
        std::vector<int> col(n);
        for (auto& c : col) c = (int)rnd(rng, 0, 5);
        const int root = (int)rnd(rng, 0, n - 1);
        std::vector<int> cnt(6, 0), distinct(n), mx(n);
        int kinds = 0, adds = 0;
        std::vector<char> in(n, 0);
        dsu_on_tree(
            g, root,
            [&](int u) {
                ST_CHECK(!in[u]);
                in[u] = 1, ++adds;
                if (cnt[col[u]]++ == 0) ++kinds;
            },
            [&](int u) {
                ST_CHECK(in[u]);
                in[u] = 0;
                if (--cnt[col[u]] == 0) --kinds;
            },
            [&](int v) { distinct[v] = kinds, mx[v] = *std::max_element(cnt.begin(), cnt.end()); });
        // 愚直：根からの親を作り、各頂点の部分木を全部見る
        std::vector<int> par(n, -2), ord;
        std::vector<int> st{root};
        par[root] = -1;
        while (!st.empty()) {
            const int v = st.back();
            st.pop_back();
            ord.push_back(v);
            for (int c : g[v])
                if (par[c] == -2) par[c] = v, st.push_back(c);
        }
        for (int v = 0; v < n; ++v) {
            std::vector<int> c2(6, 0);
            for (int u = 0; u < n; ++u) {
                int w = u;
                while (w != -1 && w != v) w = par[w];
                if (w == v) c2[col[u]]++;
            }
            int k = 0;
            for (int x : c2) k += x > 0;
            ST_CHECK(distinct[v] == k && mx[v] == *std::max_element(c2.begin(), c2.end()));
            ++cases;
        }
        int lg = 1;
        while ((1 << lg) < n) ++lg;
        ST_CHECK(adds <= n * (lg + 1));  // O(N log N) 回
    }
    std::cout << "[ds_merge_small_to_large] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_merge_small_to_large] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    std::vector<std::set<int>> init(N);
    for (int i = 0; i < N; ++i) init[i] = {i};
    auto t0 = std::chrono::steady_clock::now();
    SmallToLargeUF<std::set<int>> uf(init);
    for (int i = 1; i < N; ++i) uf.unite((int)rnd(rng, 0, i - 1), i);
    const ll t1 = elapsed_ms(t0);
    const auto g = random_tree(rng, N, 0);
    std::vector<int> col(N), cnt(N + 1, 0);
    for (auto& c : col) c = (int)rnd(rng, 0, N);
    ll kinds = 0, chk = (ll)uf.data(0).size();
    t0 = std::chrono::steady_clock::now();
    dsu_on_tree(
        g, 0, [&](int u) { kinds += cnt[col[u]]++ == 0; }, [&](int u) { kinds -= --cnt[col[u]] == 0; }, [&](int) { chk += kinds; });
    const ll t2 = elapsed_ms(t0);
    std::cout << "[ds_merge_small_to_large] benchmark (N = " << N << "): SmallToLargeUF<set<int>> 全部併合 " << t1 << " ms / dsu_on_tree（色の種類数） " << t2
              << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace stl_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace stl_demo {

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
#define SHOW(expr, comment) stl_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    stl_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("マージテク（ds_merge_small_to_large） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  SmallToLargeUF<C> uf(init): unite(x, y) / unite(x, y, merge) / find / same / size / data\n"
              << "  dsu_on_tree(g, root, add, remove, query)\n";

    api("【1】 SmallToLargeUF<C> uf(init) / unite(x, y) / data(x)", "— / bool / C&", "成分ごとに入れ物を持ち、併合で小さい方を大きい方へ移す", "全体 O(N log N) 回の挿入");
    RUN(SmallToLargeUF<std::set<int>> uf({{5}, {3}, {5}, {8}}));
    RUN(uf.unite(0, 1); uf.unite(2, 3));
    SHOW(uf.data(0).size(), "{5, 3}");
    RUN(uf.unite(1, 3));
    SHOW(uf.data(2).size(), "{3, 5, 8}（5 は 1 つ）");
    SHOW(uf.size(2), "頂点は 4 つ");

    api("【2】 unite(x, y, merge)：ABC183 F（Confluence）の形", "bool", "map<クラス, 人数> を持って、併合では人数を足す。成分の中のクラス y の人数を答える",
        "全体 O(N log^2 N)");
    {
        code("SmallToLargeUF<std::map<int, int>> uf(init);   // init[i] = {{c_i, 1}}");
        code("uf.unite(a, b, [](auto& big, auto& small) { for (auto [k, c] : small) big[k] += c; });");
        const std::vector<int> c = {1, 2, 3, 2, 1};
        std::vector<std::map<int, int>> init(5);
        for (int i = 0; i < 5; ++i) init[i][c[i]] = 1;
        SmallToLargeUF<std::map<int, int>> mc(init);
        auto merge = [](std::map<int, int>& big, std::map<int, int>& small) {
            for (auto [k, x] : small) big[k] += x;
        };
        mc.unite(0, 1, merge), mc.unite(3, 4, merge), mc.unite(1, 3, merge);
        note("クラス c = {1, 2, 3, 2, 1}、{0, 1, 3, 4} を併合 → 頂点 0 の成分のクラス 2 の人数 = " + to_s(mc.data(0)[2]) + "、クラス 3 = " + to_s(mc.data(0)[3]));
    }

    api("【3】 dsu_on_tree(g, root, add, remove, query)", "void", "各頂点 v で、v の部分木の頂点だけが add された状態で query(v) を呼ぶ", "O(N log N × add)");
    {
        code("// 部分木ごとの色の種類数");
        const std::vector<std::vector<int>> g = {{1, 2}, {0, 3, 4}, {0}, {1}, {1}};
        const std::vector<int> col = {1, 2, 1, 3, 2};
        std::vector<int> cnt(4, 0), ans(5);
        int kinds = 0;
        dsu_on_tree(
            g, 0, [&](int u) { kinds += cnt[col[u]]++ == 0; }, [&](int u) { kinds -= --cnt[col[u]] == 0; }, [&](int v) { ans[v] = kinds; });
        note("木 0-1, 0-2, 1-3, 1-4、色 {1, 2, 1, 3, 2} → 部分木の色の種類数 " + to_s(ans));
    }
}

}  // namespace stl_demo

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
    stl_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    stl_test::benchmark();    // 速度計測
    stl_demo::demo();         // 使い方デモ
    return 0;
}
