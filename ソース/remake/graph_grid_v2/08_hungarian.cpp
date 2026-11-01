// =============================================================================
//  ハンガリアン法（割り当て問題）  hungarian / hungarian_max           --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 08 フロー
//  依存       : なし（単独で貼る）
//  入力       : a[i][j] = 行 i（人）に列 j（仕事）を割り当てるコスト（n 行 m 列。0-indexed）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  各行に違う列を 1 つずつ割り当て、コストの合計を最小（最大）にする。n ≤ m なら全部の行に、n > m なら全部の列に割り当てる。
//  全部の組にコストがある（密な）割り当てなら、最小費用流より速くて書くのも楽。
//  割り当ててはいけない組には大きいコスト（10^15 など。n 個足してもあふれない値）を入れる。
//
//  機能一覧                                                                       計算量
//    hungarian(a)        最小コストの割り当て                                       O(min(n,m)^2 max(n,m))
//    hungarian_max(a)    最大（得点の合計が最大）の割り当て                          同上
//    結果 AssignmentResult r: r.cost / r.row_to_col[i] / r.col_to_row[j]（なければ -1）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    hg_test::random_test() / hg_test::benchmark() / hg_demo::demo()
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

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_HUNGARIAN
#define GG_HUNGARIAN
// ---- ハンガリアン法（グラフ・グリッド第二世代）  割り当て問題: n 行 m 列のコストで、各行に違う列を 1 つずつ割り当てる最小 ----
// n ≤ m なら全部の行に列を割り当てる。n > m なら全部の列に行を割り当てる（割り当てられない行は -1）。
// 割り当ててはいけない組には大きいコスト（10^15 など。足してもあふれない値）を入れる
struct AssignmentResult {
    long long cost = 0;             // 割り当てたコストの合計
    std::vector<int> row_to_col;    // 行 i に割り当てた列（なければ -1）
    std::vector<int> col_to_row;    // 列 j に割り当てた行（なければ -1）
};
namespace gg_internal {
// n ≤ m の場合（行をすべて割り当てる）。O(n^2 m)
inline std::vector<int> hungarian_rows(const std::vector<std::vector<long long>>& a, int n, int m) {
    const long long INF = (long long)4e18;
    std::vector<long long> u(n + 1, 0), v(m + 1, 0), minv(m + 1);
    std::vector<int> p(m + 1, 0), way(m + 1, 0);
    std::vector<char> used(m + 1);
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::fill(minv.begin(), minv.end(), INF), std::fill(used.begin(), used.end(), 0);
        do {
            used[j0] = 1;
            const int i0 = p[j0];
            long long delta = INF;
            int j1 = 0;
            const std::vector<long long>& row = a[i0 - 1];
            for (int j = 1; j <= m; ++j) {
                if (used[j]) continue;
                const long long cur = row[j - 1] - u[i0] - v[j];
                if (cur < minv[j]) minv[j] = cur, way[j] = j0;
                if (minv[j] < delta) delta = minv[j], j1 = j;
            }
            for (int j = 0; j <= m; ++j) {
                if (used[j]) u[p[j]] += delta, v[j] -= delta;
                else minv[j] -= delta;
            }
            j0 = j1;
        } while (p[j0] != 0);
        do {
            const int j1 = way[j0];
            p[j0] = p[j1], j0 = j1;
        } while (j0 != 0);
    }
    std::vector<int> r(n, -1);
    for (int j = 1; j <= m; ++j)
        if (p[j] != 0) r[p[j] - 1] = j - 1;
    return r;
}
}  // namespace gg_internal
// 最小コストの割り当て。a は n 行 m 列。O(min(n, m)^2 max(n, m))
inline AssignmentResult hungarian(const std::vector<std::vector<long long>>& a) {
    const int n = (int)a.size(), m = n ? (int)a[0].size() : 0;
    AssignmentResult res;
    res.row_to_col.assign(n, -1), res.col_to_row.assign(m, -1);
    if (n == 0 || m == 0) return res;
    if (n <= m) {
        res.row_to_col = gg_internal::hungarian_rows(a, n, m);
    } else {
        std::vector<std::vector<long long>> b(m, std::vector<long long>(n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j) b[j][i] = a[i][j];
        res.col_to_row = gg_internal::hungarian_rows(b, m, n);
        for (int j = 0; j < m; ++j) res.row_to_col[res.col_to_row[j]] = j;
    }
    for (int i = 0; i < n; ++i)
        if (res.row_to_col[i] >= 0) res.col_to_row[res.row_to_col[i]] = i, res.cost += a[i][res.row_to_col[i]];
    return res;
}
// 最大（得点の合計が最大）の割り当て
inline AssignmentResult hungarian_max(const std::vector<std::vector<long long>>& a) {
    std::vector<std::vector<long long>> b = a;
    for (auto& row : b)
        for (auto& x : row) x = -x;
    AssignmentResult r = hungarian(b);
    r.cost = -r.cost;
    return r;
}
#endif  // GG_HUNGARIAN
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "hungarian";
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

namespace hg_test {

// 愚直: 小さい側の各要素に、大きい側の異なる要素を全通り割り当てる
static long long brute(const std::vector<std::vector<long long>>& a, bool mx) {
    const int n = (int)a.size(), m = n ? (int)a[0].size() : 0;
    const int k = std::min(n, m);
    long long best = mx ? -(1LL << 62) : (1LL << 62);
    std::vector<int> used(std::max(n, m), 0);
    std::function<void(int, long long)> rec = [&](int i, long long c) {
        if (i == k) {
            best = mx ? std::max(best, c) : std::min(best, c);
            return;
        }
        for (int j = 0; j < std::max(n, m); ++j) {
            if (used[j]) continue;
            used[j] = 1;
            rec(i + 1, c + (n <= m ? a[i][j] : a[j][i]));
            used[j] = 0;
        }
    };
    rec(0, 0);
    return k == 0 ? 0 : best;
}

static void random_test() {
    std::mt19937_64 rng(20261007);
    int cases = 0;
    for (int it = 0; it < 3000 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 1, 6), m = rnd(rng, 1, 6);
        std::vector<std::vector<long long>> a(n, std::vector<long long>(m));
        for (auto& row : a)
            for (auto& x : row) x = (it % 4 == 0) ? (long long)(rng() % 2000000000000ULL) - 1000000000000LL : rnd(rng, -10, 10);
        for (bool mx : {false, true}) {
            const AssignmentResult r = mx ? hungarian_max(a) : hungarian(a);
            CHECK(r.cost == brute(a, mx));
            CHECK((int)r.row_to_col.size() == n && (int)r.col_to_row.size() == m);
            long long c = 0;
            int cnt = 0;
            for (int i = 0; i < n; ++i)
                if (r.row_to_col[i] >= 0) {
                    CHECK(r.col_to_row[r.row_to_col[i]] == i);
                    c += a[i][r.row_to_col[i]], ++cnt;
                }
            for (int j = 0; j < m; ++j)
                if (r.col_to_row[j] >= 0) CHECK(r.row_to_col[r.col_to_row[j]] == j);
            CHECK(c == r.cost && cnt == std::min(n, m));
        }
    }
    CHECK(hungarian({}).cost == 0 && hungarian({}).row_to_col.empty());
    std::cout << "[hungarian] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[hungarian] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {300, 1000}) {
        if (TEST_SCALE == 1 && n > 300) break;
        std::vector<std::vector<long long>> a(n, std::vector<long long>(n));
        for (auto& row : a)
            for (auto& x : row) x = (long long)(rng() % 1000000000);
        auto t0 = std::chrono::steady_clock::now();
        const AssignmentResult r = hungarian(a);
        std::cout << "[hungarian] benchmark: " << n << " × " << n << "  " << ms_since(t0) << " ms  (cost " << r.cost % 1000 << ")\n";
    }
}

}  // namespace hg_test


namespace hg_demo {
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
    banner("ハンガリアン法（割り当て問題） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"auto r = hungarian(a);", "AssignmentResult", "コスト最小の割り当て（a は n 行 m 列）"},
        {"auto r = hungarian_max(a);", "AssignmentResult", "得点最大の割り当て"},
        {"r.cost / r.row_to_col[i] / r.col_to_row[j]", "", "合計 / 行 i の列 / 列 j の行（なければ -1）"},
    };
    cheat_sheet("API 一覧", rows, 44, 18);

    head("デモで使うコスト: 3 人 × 3 つの仕事");
    RUN(std::vector<std::vector<long long>> a = {{4, 1, 3}, {2, 0, 5}, {3, 2, 2}});

    api("【1】 hungarian(a)", "AssignmentResult", "各人に違う仕事を 1 つずつ。コストの合計の最小", "O(n^3)");
    RUN(auto r = hungarian(a));
    SHOW(r.cost, "1 + 2 + 2");
    SHOW(r.row_to_col, "人 i の仕事");
    SHOW(r.col_to_row, "仕事 j の人");

    api("【2】 hungarian_max(a)", "AssignmentResult", "得点の合計の最大", "O(n^3)");
    SHOW(hungarian_max(a).cost, "4 + 5 + 2（人 0 → 仕事 0、人 1 → 仕事 2、人 2 → 仕事 1）");

    api("【3】 行と列の数が違うとき", "AssignmentResult", "n ≤ m なら全部の行に、n > m なら全部の列に割り当てる（余りは -1）", "O(min^2 max)");
    RUN(std::vector<std::vector<long long>> b = {{5, 1, 9, 2}, {8, 7, 1, 6}});
    RUN(auto rb = hungarian(b));
    SHOW(rb.cost, "人 0 → 仕事 1、人 1 → 仕事 2");
    SHOW(rb.col_to_row, "使われない仕事は -1");

    head("典型的な使い方");
    note("人と仕事の割り当て（全部の組にコスト）: hungarian(a)。N ≤ 500 程度なら一瞬、1000 でも 1 秒かからない");
    note("割り当ててはいけない組                : 大きいコスト（10^15 など）を入れ、結果のコストが大きすぎたら「不可能」");
    note("一部だけ割り当てる・容量がある        : 08_min_cost_flow");
}

}  // namespace hg_demo

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
    hg_test::random_test();
    hg_test::benchmark();
    hg_demo::demo();
    return 0;
}
