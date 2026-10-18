// =============================================================================
//  部分永続 Union-Find（過去の時刻を問い合わせられる）  PersistentUnionFind   --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 03 Union-Find・最小全域木
//  依存       : なし（単独で貼れる）
//  頂点番号   : 0-indexed（0 ~ n-1）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  unite のたびに時刻が 1 進む（k 回目の unite の直後が時刻 k）。過去の時刻 t の状態を、あとからいつでも問い合わせられる:
//    - 「u と v が初めてつながったのはいつか」（when）… 並列二分探索やクエリの先読みをしなくてよい
//    - 「時刻 t の時点で同じ成分か・成分の大きさ・成分の数」
//  変更（unite）は今の時刻にだけできる（過去を書き換えることはできない = 部分永続）
//
//  機能一覧                                                                       計算量
//    PersistentUnionFind uf(n)      頂点数 n、時刻 0                                O(n)
//    uf.unite(x, y)        -> bool  併合して時刻を 1 進める。新しく併合したら true   O(log n)
//    uf.now()              -> int   今の時刻（unite を呼んだ回数。併合しなかった回も数える）  O(1)
//    uf.find(x, t = 今)    -> int   時刻 t での根                                   O(log n)
//    uf.same(x, y, t = 今) -> bool  時刻 t で同じ成分か                             O(log n)
//    uf.size(x, t = 今)    -> int   時刻 t での x の成分の大きさ                     O(log n)
//    uf.group_count(t = 今) -> int  時刻 t での成分の数                              O(1)
//    uf.when(x, y)         -> int   初めて同じ成分になった時刻（x == y なら 0、今もつながっていなければ -1）  O(log n)
//    uf.n() / uf.init(n)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    puf_test::random_test() / puf_test::benchmark() / puf_demo::demo()
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
#include <climits>

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_PERSISTENT_UNION_FIND
#define GG_PERSISTENT_UNION_FIND
// ---- 部分永続 Union-Find（グラフ・グリッド第二世代）  過去の時刻の状態を問い合わせられる Union-Find ----
// 時刻 = それまでに unite を呼んだ回数（最初は 0。k 回目の unite の直後が時刻 k。併合しなかった unite も 1 回と数える）。
// 変更（unite）は今の時刻にだけでき、問い合わせ（find / same / size / group_count）は 0 ~ now() のどの時刻でもできる。
// 経路圧縮をしない（union by size のみ）ので find は O(log n)。必要なヘッダ: <vector>, <utility>, <climits>
struct PersistentUnionFind {
    PersistentUnionFind() : PersistentUnionFind(0) {}
    explicit PersistentUnionFind(int n) { init(n); }
    void init(int n) {
        n_ = n, now_ = 0;
        par_.resize(n), time_.assign(n, INT_MAX), hist_.assign(n, {});
        for (int i = 0; i < n; ++i) par_[i] = i, hist_[i].push_back({0, 1});
        cnt_.assign(1, n);
    }
    int n() const { return n_; }
    int now() const { return now_; }  // 今の時刻（unite を呼んだ回数）
    // x と y を併合する（時刻が 1 進む）。新しく併合したら true  O(log n)
    bool unite(int x, int y) {
        ++now_;
        x = find(x), y = find(y);
        if (x == y) {
            cnt_.push_back(cnt_.back());
            return false;
        }
        int sx = hist_[x].back().second, sy = hist_[y].back().second;
        if (sx < sy) std::swap(x, y), std::swap(sx, sy);
        par_[y] = x, time_[y] = now_;
        hist_[x].push_back({now_, sx + sy});
        cnt_.push_back(cnt_.back() - 1);
        return true;
    }
    // 時刻 t での x の根（t を省くと今）  O(log n)
    int find(int x, int t = INT_MAX) const {
        while (par_[x] != x && time_[x] <= t) x = par_[x];
        return x;
    }
    bool same(int x, int y, int t = INT_MAX) const { return find(x, t) == find(y, t); }
    // 時刻 t での x の成分の大きさ  O(log n)
    int size(int x, int t = INT_MAX) const {
        const std::vector<std::pair<int, int>>& h = hist_[find(x, t)];
        int lo = 0, hi = (int)h.size();  // h[k].first <= t となる最後の k（h[0].first = 0）
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            (h[mid].first <= t ? lo : hi) = mid;
        }
        return h[lo].second;
    }
    // 時刻 t での成分の個数  O(1)
    int group_count(int t = INT_MAX) const { return cnt_[t < now_ ? t : now_]; }
    // x と y が初めて同じ成分になった時刻（最初から同じ（x == y）なら 0、今もつながっていなければ -1）  O(log n)
    int when(int x, int y) const {
        if (!same(x, y)) return -1;
        int t = 0;
        while (x != y) {  // つながった時刻の早いほうを親へ進める（親へ進むほど時刻は大きくなる）
            if (time_[x] < time_[y]) t = time_[x], x = par_[x];
            else t = time_[y], y = par_[y];
        }
        return t;
    }

  private:
    int n_ = 0, now_ = 0;
    std::vector<int> par_, time_;                        // 親、親につながった時刻（根は INT_MAX）
    std::vector<std::vector<std::pair<int, int>>> hist_;  // 根だった間の (時刻, 成分の大きさ) の履歴
    std::vector<int> cnt_;                               // cnt_[t] = 時刻 t での成分の個数
};
#endif  // GG_PERSISTENT_UNION_FIND
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "persistent_union_find";
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

namespace puf_test {

struct Naive {  // 時刻ごとに成分の番号を全部持つ
    std::vector<std::vector<int>> comp;  // comp[t][v]
    explicit Naive(int n) {
        comp.push_back(std::vector<int>(n));
        for (int v = 0; v < n; ++v) comp[0][v] = v;
    }
    void unite(int x, int y) {
        std::vector<int> c = comp.back();
        const int a = c[x], b = c[y];
        for (int& z : c)
            if (z == b) z = a;
        comp.push_back(c);
    }
};

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 10), q = rnd(rng, 0, 25);
        PersistentUnionFind uf(n);
        Naive nv(n);
        CHECK(uf.n() == n && uf.now() == 0 && uf.group_count() == n);
        for (int k = 0; k < q; ++k) {
            const int x = rnd(rng, 0, n - 1), y = rnd(rng, 0, n - 1);
            const bool was = nv.comp.back()[x] == nv.comp.back()[y];
            CHECK(uf.unite(x, y) == !was);
            nv.unite(x, y);
            CHECK(uf.now() == k + 1);
        }
        for (int t = 0; t <= q + 2; ++t) {
            const std::vector<int>& c = nv.comp[std::min(t, q)];
            int groups = 0;
            for (int v = 0; v < n; ++v) {
                int sz = 0;
                for (int u = 0; u < n; ++u) sz += (c[u] == c[v]);
                CHECK(uf.size(v, t) == sz);
                bool first = true;
                for (int u = 0; u < v; ++u)
                    if (c[u] == c[v]) first = false;
                groups += first;
                for (int u = 0; u < n; ++u) CHECK(uf.same(u, v, t) == (c[u] == c[v]));
                CHECK(c[uf.find(v, t)] == c[v]);
                CHECK(uf.find(uf.find(v, t), t) == uf.find(v, t));
            }
            CHECK(uf.group_count(t) == groups);
        }
        for (int u = 0; u < n; ++u)
            for (int v = 0; v < n; ++v) {
                int w = -1;
                for (int t = 0; t <= q; ++t)
                    if (nv.comp[t][u] == nv.comp[t][v]) {
                        w = t;
                        break;
                    }
                CHECK(uf.when(u, v) == w);
            }
        // 今の時刻の問い合わせ（t を省く）
        for (int v = 0; v < n; ++v) CHECK(uf.size(v) == uf.size(v, q) && uf.find(v) == uf.find(v, q));
    }
    std::cout << "[persistent_union_find] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[persistent_union_find] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    const int n = 200000, q = 200000;
    auto t0 = std::chrono::steady_clock::now();
    PersistentUnionFind uf(n);
    for (int k = 0; k < q; ++k) uf.unite((int)(rng() % n), (int)(rng() % n));
    const long long t1 = ms_since(t0);
    t0 = std::chrono::steady_clock::now();
    long long chk = 0;
    for (int k = 0; k < q; ++k) {
        const int u = (int)(rng() % n), v = (int)(rng() % n), t = (int)(rng() % (q + 1));
        chk += uf.when(u, v) + uf.size(u, t) + uf.same(u, v, t);
    }
    std::cout << "[persistent_union_find] benchmark: N = " << n << "  unite " << q << " 回 " << t1 << " ms, when + size + same を " << q << " 回 " << ms_since(t0)
              << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace puf_test


namespace puf_demo {
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
    banner("部分永続 Union-Find 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"PersistentUnionFind uf(n);", "", "時刻 0（全部バラバラ）"},
        {"uf.unite(x, y)", "bool", "併合して時刻を 1 進める"},
        {"uf.now()", "int", "今の時刻 = unite を呼んだ回数"},
        {"uf.same(x, y, t) / uf.size(x, t) / uf.find(x, t)", "", "時刻 t の状態（t を省くと今）"},
        {"uf.group_count(t)", "int", "時刻 t の成分の数"},
        {"uf.when(x, y)", "int", "初めてつながった時刻（つながっていなければ -1）"},
    };
    cheat_sheet("API 一覧", rows, 50, 10);

    api("【1】 時刻を進めながら併合する", "bool", "k 回目の unite の直後が時刻 k", "O(log n)");
    RUN(PersistentUnionFind uf(5));
    SHOW(uf.unite(0, 1), "時刻 1");
    SHOW(uf.unite(2, 3), "時刻 2");
    SHOW(uf.unite(1, 0), "時刻 3（もう同じ。時刻は進む）");
    SHOW(uf.unite(1, 3), "時刻 4");
    SHOW(uf.now(), "");

    api("【2】 過去の時刻の状態", "bool / int", "same(x, y, t) / size(x, t) / group_count(t)", "O(log n)");
    SHOW(uf.same(0, 3, 3), "時刻 3 ではまだ別");
    SHOW(uf.same(0, 3, 4), "時刻 4 でつながる");
    SHOW(uf.size(0, 1), "時刻 1 の {0, 1}");
    SHOW(uf.size(0), "今は {0, 1, 2, 3}");
    SHOW(uf.group_count(0), "");
    SHOW(uf.group_count(2), "");
    SHOW(uf.group_count(), "");

    api("【3】 いつ初めてつながったか", "int", "x == y なら 0、今もつながっていなければ -1", "O(log n)");
    SHOW(uf.when(0, 3), "");
    SHOW(uf.when(2, 3), "");
    SHOW(uf.when(4, 0), "");
    SHOW(uf.when(4, 4), "");

    head("【4】 例: 年が新しい道から順に unite（ABC040 D 道路の老朽化対策）");
    note("道を「年が新しい順」に並べて unite していく。「年 w より新しい道だけ使うと、v から何都市行けるか」は、");
    note("年 w より新しい道の本数を t として uf.size(v, t)（t は二分探索で求める）。クエリを先読みして並べ替えなくてよい。");
    head("【5】 例: 辺を順に足したとき、u と v が初めてつながる辺は何本目か");
    note("uf.when(u, v)。最小全域木の順に足せば「u - v のパスの最大の辺をなるべく小さく」した値もわかる（最大値の最小化）。");

    head("注意");
    note("経路圧縮をしないので find は O(log n)（今の状態だけなら 03_union_find のほうが速い）。");
    note("取り消し（undo）がしたいときは 03_rollback_union_find。");
}

}  // namespace puf_demo

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
    puf_test::random_test();
    puf_test::benchmark();
    puf_demo::demo();
    return 0;
}
