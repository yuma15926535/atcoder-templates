// =============================================================================
//  巻き戻せる Union-Find  RollbackUnionFind                          --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 03 Union-Find・最小全域木
//  依存       : なし（単独で貼る）
//  由来       : 9/23 作成の dynamic_connectivity_offline.cpp に入っている RollbackUnionFind をそのまま取り出し、
//               第二世代の形（ファイル名・#ifndef ガード）にした。API は同じ
//  頂点番号   : 0-indexed（0 ~ n-1）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る（<cassert> を使う）
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ふつうの UnionFind と違い、unite を「後から取り消せる」。経路圧縮をしないので find は O(log n)。
//  使いどころ: DFS で木（操作の履歴の木など）をたどり、戻るときに取り消す / 試しに併合して元に戻す /
//              辺の追加と削除が混ざる問題（オフラインなら dynamic_connectivity_offline.cpp）
//  ※ dynamic_connectivity_offline.cpp にも同じ RollbackUnionFind が入っている。両方を貼るときはこちらを貼らない
//
//  機能一覧                                                                       計算量
//    RollbackUnionFind uf(n)            頂点数 n で作る                              O(n)
//    uf.unite(x, y)        -> bool      併合。新しく併合したら true（false でも 1 回と数える）  O(log n)
//    uf.same(x, y) / uf.find(x) / uf.size(x)                                         O(log n)
//    uf.group_count()                   グループの個数                               O(1)
//    uf.undo()                          直前の unite 1 回分を取り消す                O(1)
//    uf.snapshot()         -> int       今の状態の番号                               O(1)
//    uf.rollback(s)                     snapshot() が s を返した状態まで戻す         O(戻す回数)
//    uf.groups()                        全グループ（最小要素の昇順・各グループ内も昇順）  O(n log n)
//    uf.init(n) / uf.n()                作り直す / 頂点数
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    ruf_test::random_test() / ruf_test::benchmark() / ruf_demo::demo()
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <cassert>

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_ROLLBACK_UNION_FIND
#define GG_ROLLBACK_UNION_FIND
// ---- 巻き戻せる Union-Find（9/23 版 dynamic_connectivity_offline.cpp の RollbackUnionFind をそのまま採用） ----
// 経路圧縮をしない代わりに、unite を 1 回ずつ取り消せる。find は O(log n)。必要なヘッダ: <vector>, <utility>, <cassert>
// ※ dynamic_connectivity_offline.cpp にも同じ RollbackUnionFind が入っているので、両方を貼るときはこちらを貼らない
struct RollbackUnionFind {
  public:
    RollbackUnionFind() : RollbackUnionFind(0) {}
    explicit RollbackUnionFind(int n) { init(n); }

    // 頂点数 n で作り直す（履歴も消える）  O(n)
    void init(int n) {
        n_ = n;
        cnt_ = n;
        par_.assign(n, -1);
        hist_.clear();
    }

    // 頂点数
    int n() const { return n_; }

    // x の根（代表元）  O(log n)
    int find(int x) const {
        while (par_[x] >= 0) x = par_[x];
        return x;
    }

    // 併合。新しく併合したら true。false の場合も「1 回の unite」として履歴に積む  O(log n)
    bool unite(int x, int y) {
        x = find(x);
        y = find(y);
        if (x == y) {
            hist_.emplace_back(-1, 0);
            return false;
        }
        if (par_[x] > par_[y]) std::swap(x, y);  // サイズの大きい方 x を根にする
        hist_.emplace_back(y, par_[y]);
        par_[x] += par_[y];
        par_[y] = x;
        --cnt_;
        return true;
    }

    // x と y が同じグループなら true  O(log n)
    bool same(int x, int y) const { return find(x) == find(y); }

    // x のグループの要素数  O(log n)
    int size(int x) const { return -par_[find(x)]; }

    // グループの個数  O(1)
    int group_count() const { return cnt_; }

    // 直前の unite を 1 回分取り消す（false だった unite も 1 回と数える）  O(1)
    void undo() {
        assert(!hist_.empty());
        const int y = hist_.back().first, py = hist_.back().second;
        hist_.pop_back();
        if (y < 0) return;
        const int x = par_[y];
        par_[x] -= py;
        par_[y] = py;
        ++cnt_;
    }

    // 今の状態を表す番号（= 取り消していない unite の回数）。rollback(s) でこの状態に戻れる  O(1)
    int snapshot() const { return (int)hist_.size(); }

    // snapshot() が s を返した状態まで戻す  O(戻す回数)
    void rollback(int s) {
        while ((int)hist_.size() > s) undo();
    }

    // 全グループ（最小要素の昇順・各グループ内も昇順）  O(n log n)
    std::vector<std::vector<int>> groups() const {
        std::vector<int> id(n_, -1);
        std::vector<std::vector<int>> res;
        for (int i = 0; i < n_; ++i) {
            const int r = find(i);
            if (id[r] < 0) {
                id[r] = (int)res.size();
                res.emplace_back();
            }
            res[id[r]].push_back(i);
        }
        return res;
    }

  private:
    int n_ = 0;
    int cnt_ = 0;
    std::vector<int> par_;                   // 根なら -(サイズ)、それ以外は親
    std::vector<std::pair<int, int>> hist_;  // (吸収された根, その元の par_ 値)。false の unite は (-1, 0)
};
#endif  // GG_ROLLBACK_UNION_FIND
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "rollback_union_find";
// ---- テスト用の小道具 ----
#define CHECK(...)                                                                                \
    do {                                                                                          \
        if (!(__VA_ARGS__)) {                                                                     \
            std::cerr << "[" << TEST_TAG << "] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                         \
        }                                                                                         \
    } while (0)
// テスト量の倍率。デバッグビルド（最適化なし・範囲チェックあり）は遅いので自動で 1/8 にする
#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif
[[maybe_unused]] static int rnd(std::mt19937_64& rng, int lo, int hi) { return lo + (int)(rng() % (unsigned long long)(hi - lo + 1)); }
[[maybe_unused]] static long long ms_since(std::chrono::steady_clock::time_point t0) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
struct CinFromString {  // std::cin を文字列に差し替える
    std::istringstream ss;
    std::streambuf* old;
    explicit CinFromString(const std::string& s) : ss(s), old(std::cin.rdbuf(ss.rdbuf())) {}
    ~CinFromString() { std::cin.rdbuf(old); }
};

namespace ruf_test {

// 愚直: 今までの unite の列（取り消したものは除く）からラベルを毎回作り直す
struct Naive {
    int n;
    std::vector<std::pair<int, int>> ops;
    std::vector<int> label() const {
        std::vector<int> lab(n);
        std::iota(lab.begin(), lab.end(), 0);
        for (bool ch = true; ch;) {
            ch = false;
            for (auto [a, b] : ops) {
                const int m = std::min(lab[a], lab[b]);
                if (lab[a] != m || lab[b] != m) lab[a] = lab[b] = m, ch = true;
            }
        }
        return lab;
    }
};

static void check(const RollbackUnionFind& uf, const Naive& nv) {
    const std::vector<int> lab = nv.label();
    const int n = nv.n;
    std::vector<int> cnt(n, 0);
    int groups = 0;
    for (int v = 0; v < n; ++v) groups += (lab[v] == v), ++cnt[lab[v]];
    CHECK(uf.group_count() == groups && uf.n() == n);
    for (int a = 0; a < n; ++a) {
        CHECK(uf.size(a) == cnt[lab[a]] && lab[uf.find(a)] == lab[a]);
        for (int b = 0; b < n; ++b) CHECK(uf.same(a, b) == (lab[a] == lab[b]));
    }
    std::vector<std::vector<int>> g;
    std::vector<int> id(n, -1);
    for (int v = 0; v < n; ++v) {
        if (id[lab[v]] < 0) id[lab[v]] = (int)g.size(), g.emplace_back();
        g[id[lab[v]]].push_back(v);
    }
    CHECK(uf.groups() == g);
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 12), q = rnd(rng, 1, 40);
        RollbackUnionFind uf(n);
        Naive nv{n, {}};
        std::vector<std::pair<int, size_t>> snaps;  // (snapshot の値, そのときの ops の長さ)
        for (int k = 0; k < q; ++k) {
            const int t = rnd(rng, 0, 9);
            if (t <= 5) {
                const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
                const bool was = uf.same(a, b);
                CHECK(uf.unite(a, b) == !was);
                nv.ops.emplace_back(a, b);
            } else if (t <= 7) {
                if (nv.ops.empty()) continue;
                uf.undo();
                nv.ops.pop_back();
            } else if (t == 8) {
                CHECK(uf.snapshot() == (int)nv.ops.size());
                snaps.emplace_back(uf.snapshot(), nv.ops.size());
            } else {
                // 今の長さ以下の snapshot に戻る
                std::vector<std::pair<int, size_t>> ok;
                for (auto& s : snaps)
                    if (s.second <= nv.ops.size()) ok.push_back(s);
                if (ok.empty()) continue;
                const auto s = ok[rng() % ok.size()];
                uf.rollback(s.first);
                nv.ops.resize(s.second);
            }
            check(uf, nv);
        }
        if (it % 50 == 0) {
            uf.init(n);
            nv.ops.clear();
            check(uf, nv);
        }
    }
    std::cout << "[rollback_union_find] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[rollback_union_find] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {200000, 1000000}) {
        if (TEST_SCALE == 1 && n > 200000) break;
        RollbackUnionFind uf(n);
        auto t0 = std::chrono::steady_clock::now();
        long long chk = 0;
        for (int k = 0; k < 4 * n; ++k) {
            const int t = (int)(rng() % 4);
            const int a = (int)(rng() % n), b = (int)(rng() % n);
            if (t <= 1) uf.unite(a, b);
            else if (t == 2) chk += uf.same(a, b);
            else if (uf.snapshot() > 0) uf.undo();
        }
        std::cout << "[rollback_union_find] benchmark: n = " << n << ", " << 4 * n << " 回（unite : same : undo = 2 : 1 : 1）" << ms_since(t0) << " ms  (check "
                  << chk + uf.group_count() << ")\n";
    }
}

}  // namespace ruf_test


namespace ruf_demo {
// ---- デモ用の小道具（見方: コード -> 返り値 // 意味） ----
template <class A, class B>
[[maybe_unused]] static std::string to_s(const std::pair<A, B>& p);
template <class T, size_t N>
[[maybe_unused]] static std::string to_s(const std::array<T, N>& a);
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v);
#ifdef GG_GRAPH_BASE
[[maybe_unused]] static std::string to_s(const Edge& e) { return "{to=" + std::to_string(e.to) + ", id=" + std::to_string(e.id) + ", w=" + std::to_string(e.w) + "}"; }
[[maybe_unused]] static std::string to_s(const RawEdge& e) { return "{" + std::to_string(e.from) + "-" + std::to_string(e.to) + ", w=" + std::to_string(e.w) + "}"; }
#endif
[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(char c) { return std::string("'") + c + "'"; }
[[maybe_unused]] static std::string to_s(const std::string& s) { return "\"" + s + "\""; }
[[maybe_unused]] static std::string to_s(const char* s) { return "\"" + std::string(s) + "\""; }
template <class T, std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool> && !std::is_same_v<T, char>, int> = 0>
[[maybe_unused]] static std::string to_s(T v) {
    if (std::is_same_v<T, long long> && (long long)v == (1LL << 60)) return "DIST_INF";
    return std::to_string(v);
}
template <class A, class B>
[[maybe_unused]] static std::string to_s(const std::pair<A, B>& p) { return "(" + to_s(p.first) + ", " + to_s(p.second) + ")"; }
template <class T, size_t N>
[[maybe_unused]] static std::string to_s(const std::array<T, N>& a) {
    std::string s = "(";
    for (size_t i = 0; i < N; ++i) s += (i ? ", " : "") + to_s(a[i]);
    return s + ")";
}
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
template <class T, class = decltype(std::declval<const T&>().to_vector())>
[[maybe_unused]] static std::string to_s(const T& v) { return to_s(v.to_vector()); }  // IntSpan など
[[maybe_unused]] static std::string pad(const std::string& s, size_t w) {  // 全角を 2 文字分として幅をそろえる
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
[[maybe_unused]] static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
[[maybe_unused]] static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
[[maybe_unused]] static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
[[maybe_unused]] static void stmt(const std::string& s) {
    const bool block = !s.empty() && s.back() == '}' &&
                       (s.compare(0, 4, "for ") == 0 || s.compare(0, 3, "if ") == 0 || s.compare(0, 6, "while ") == 0);
    code(block ? s : s + ";");
}
template <class V>
[[maybe_unused]] static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 36) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 56) + "  // " + comment;
    std::cout << line << "\n";
}
struct Row { const char* use; const char* ret; const char* what; };
template <size_t N>
[[maybe_unused]] static void cheat_sheet(const std::string& title, const Row (&rows)[N], size_t w1 = 46, size_t w2 = 22) {
    std::cout << "\n■ " << title << "\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, w1) << pad(ret, w2) << r.what << "\n";
    }
}
#define SHOW(expr, comment) show(#expr, (expr), comment)
#define RUN(...)             \
    stmt(#__VA_ARGS__);      \
    __VA_ARGS__

static void demo() {
    banner("巻き戻せる Union-Find 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"RollbackUnionFind uf(n);", "", "頂点数 n"},
        {"uf.unite(x, y)", "bool", "併合（新しく併合したら true）"},
        {"uf.same(x, y) / uf.find(x) / uf.size(x)", "", "同じグループか / 根 / 要素数（O(log n)）"},
        {"uf.group_count()", "int", "グループの個数"},
        {"uf.undo()", "", "直前の unite を 1 回取り消す（false だった unite も 1 回）"},
        {"int s = uf.snapshot(); ... uf.rollback(s);", "", "s の状態まで戻す"},
        {"uf.groups()", "vector<vector<int>>", "全グループ"},
    };
    cheat_sheet("API 一覧", rows, 44, 24);

    api("【1】 unite と undo", "bool / なし", "unite の履歴を積んでおき、新しいものから 1 回ずつ取り消せる", "unite O(log n) / undo O(1)");
    RUN(RollbackUnionFind uf(5));
    SHOW(uf.unite(0, 1), "");
    SHOW(uf.unite(1, 2), "");
    SHOW(uf.unite(0, 2), "すでに同じ（false）。これも 1 回と数える");
    SHOW(uf.groups(), "");
    RUN(uf.undo());
    SHOW(uf.same(0, 2), "false だった unite を取り消しただけなので、まだつながっている");
    RUN(uf.undo());
    SHOW(uf.same(0, 2), "1-2 の unite を取り消した");
    SHOW(uf.groups(), "");

    api("【2】 snapshot と rollback", "int / なし", "今の状態の番号を覚えておき、まとめてそこまで戻す", "O(戻す回数)");
    RUN(int s = uf.snapshot());
    RUN(uf.unite(2, 3); uf.unite(3, 4); uf.unite(0, 4));
    SHOW(uf.group_count(), "全部つながった");
    RUN(uf.rollback(s));
    SHOW(uf.group_count(), "s の時点に戻った（{0, 1}, {2}, {3}, {4}）");

    head("典型的な使い方");
    code("// 操作の履歴が木になっている問題（バージョン管理・オフラインのクエリ）: 木を DFS して、");
    code("//   入るときに unite、出るときに undo する。各頂点で「その時点の状態」の質問に答えられる");
    code("// 「この辺を足したら矛盾するか」を試す: s = uf.snapshot(); 試しに unite; 調べる; uf.rollback(s);");
    note("辺の削除も混ざる連結性の問題（全部の入力を読んでから答えてよい）は dynamic_connectivity_offline.cpp。");
    note("削除しかないなら、逆順に見て「追加だけ」にすれば 03_union_find でよい（巻き戻し不要）。");
    note("二部グラフ判定を巻き戻したいときは、頂点 v を v と v + n の 2 つにして unite(a, b + n), unite(a + n, b) とし、same(v, v + n) なら矛盾。");
}

}  // namespace ruf_demo

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
    ruf_test::random_test();
    ruf_test::benchmark();
    ruf_demo::demo();
    return 0;
}
