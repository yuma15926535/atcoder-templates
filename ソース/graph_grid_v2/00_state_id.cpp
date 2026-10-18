// =============================================================================
//  状態の番号付け StateId（(頂点, 状態) ↔ 連番）                  --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 00 基盤（状態付き探索の土台。ほかのファイルには同梱しないので、使うときに一緒に貼る）
//  依存       : なし
//  番号       : (a, b, c, d) → ((a * B + b) * C + c) * D + d。a が同じ状態は番号が連続する
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  状態付き探索の組み方（専用テンプレートは作らない）
//    1. 状態を (頂点, 状態) の組にして、StateId S(N, K) で番号を付ける（3 つ・4 つの組も可）
//    2. 状態の数 S.size() が 10^7 程度以下か確かめる
//    3. 移動のコストで探索の本体を選ぶ：全部 1 → bfs_core、0 か 1 → zero_one_bfs_core、0 以上 → dijkstra_core
//         auto r = dijkstra_core(S.size(), S.id(s, 0), [&](int x, auto&& push) {
//             auto [v, k] = S.get(x);
//             for (auto& e : G[v]) { push(S.id(e.to, k), e.w); if (k < K) push(S.id(e.to, k + 1), 0); }
//         });
//    4. 答えは終点のすべての状態の最小：for (int x : S.ids_of(t)) ans = min(ans, r.dist[x]);
//
//  機能一覧                                                                    計算量
//    StateId S(A, B) / S(A, B, C) / S(A, B, C, D)   各成分の個数を指定して作る    O(1)
//    S.size()                            状態の数 A*B*C*D                         O(1)
//    S.id(a, b) / S.id(a, b, c) / S.id(a, b, c, d)   状態の番号                    O(1)
//    S.get(x)            -> pair         (a, b)。auto [v, k] = S.get(x);          O(1)
//    S.get3(x) / S.get4(x) -> array      (a, b, c) / (a, b, c, d)                 O(1)
//    S.first(x)                          1 つ目の成分（頂点・マスの番号）だけ      O(1)
//    S.in(a, b, ...)                     各成分が範囲内か                          O(1)
//    S.ids_of(a)         -> vector<int>  1 つ目の成分が a の状態すべて             O(B*C*D)
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    sid_test::random_test()  … 全機能のランダムテスト
//    sid_demo::demo()         … 使い方デモ（状態付きダイクストラの組み方の例つき）
// =============================================================================
#include <algorithm>
#include <array>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_STATE_ID
#define GG_STATE_ID
// ---- 状態の番号付け（グラフ・グリッド第二世代）  (a, b[, c[, d]]) ↔ 連番 0 ~ size()-1 ----
struct StateId {
    int A = 0, B = 1, C = 1, D = 1;  // 各成分の個数（a は 0 ~ A-1、b は 0 ~ B-1、...）
    int cd = 1, bcd = 1;             // C*D, B*C*D
    StateId() {}
    StateId(int a, int b, int c = 1, int d = 1) : A(a), B(b), C(c), D(d), cd(c * d), bcd(b * c * d) {}
    int size() const { return A * bcd; }  // 状態の数（int に収まる前提。目安 10^7 程度まで）
    int id(int a, int b, int c = 0, int d = 0) const { return a * bcd + b * cd + c * D + d; }
    bool in(int a, int b, int c = 0, int d = 0) const {  // 各成分が範囲内か
        return (unsigned)a < (unsigned)A && (unsigned)b < (unsigned)B && (unsigned)c < (unsigned)C && (unsigned)d < (unsigned)D;
    }
    int first(int x) const { return x / bcd; }  // 1 つ目の成分（頂点・マスの番号）だけ取り出す
    std::pair<int, int> get(int x) const { return {x / bcd, x % bcd / cd}; }  // (a, b)。auto [v, k] = S.get(x);
    std::array<int, 3> get3(int x) const { return {x / bcd, x % bcd / cd, x % cd / D}; }  // (a, b, c)
    std::array<int, 4> get4(int x) const { return {x / bcd, x % bcd / cd, x % cd / D, x % D}; }  // (a, b, c, d)
    std::vector<int> ids_of(int a) const {  // 1 つ目の成分が a の状態すべて（終点のすべての状態の最小を取るときなど）
        std::vector<int> r(bcd);
        for (int k = 0; k < bcd; ++k) r[k] = a * bcd + k;
        return r;
    }
};
#endif  // GG_STATE_ID
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト
// =============================================================================
namespace sid_test {

#define SID_CHECK(...)                                                                           \
    do {                                                                                         \
        if (!(__VA_ARGS__)) {                                                                    \
            std::cerr << "[state_id] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                        \
        }                                                                                        \
    } while (0)

static int rnd(std::mt19937_64& rng, int lo, int hi) { return lo + (int)(rng() % (unsigned long long)(hi - lo + 1)); }

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 3000; ++it, ++cases) {
        const int dims = rnd(rng, 2, 4);
        const int A = rnd(rng, 1, 6), B = rnd(rng, 1, 6), C = dims >= 3 ? rnd(rng, 1, 5) : 1, D = dims >= 4 ? rnd(rng, 1, 4) : 1;
        const StateId S = (dims == 2) ? StateId(A, B) : (dims == 3) ? StateId(A, B, C) : StateId(A, B, C, D);
        SID_CHECK(S.size() == A * B * C * D);
        std::vector<int> seen(S.size(), 0);
        int expect = 0;  // 辞書順に並べたときの番号と一致する
        for (int a = 0; a < A; ++a)
            for (int b = 0; b < B; ++b)
                for (int c = 0; c < C; ++c)
                    for (int d = 0; d < D; ++d, ++expect) {
                        const int x = (dims == 2) ? S.id(a, b) : (dims == 3) ? S.id(a, b, c) : S.id(a, b, c, d);
                        SID_CHECK(x == expect && S.id(a, b, c, d) == x);
                        SID_CHECK(0 <= x && x < S.size() && !seen[x]);
                        seen[x] = 1;
                        SID_CHECK(S.first(x) == a);
                        SID_CHECK(S.get(x) == std::make_pair(a, b));
                        SID_CHECK(S.get3(x) == (std::array<int, 3>{a, b, c}));
                        SID_CHECK(S.get4(x) == (std::array<int, 4>{a, b, c, d}));
                        SID_CHECK(S.in(a, b, c, d));
                    }
        for (int k = 0; k < 30; ++k) {
            const int a = rnd(rng, -2, A + 1), b = rnd(rng, -2, B + 1), c = rnd(rng, -2, C + 1), d = rnd(rng, -2, D + 1);
            SID_CHECK(S.in(a, b, c, d) == (0 <= a && a < A && 0 <= b && b < B && 0 <= c && c < C && 0 <= d && d < D));
        }
        for (int a = 0; a < A; ++a) {
            const std::vector<int> ids = S.ids_of(a);
            SID_CHECK((int)ids.size() == B * C * D);
            for (int x : ids) SID_CHECK(S.first(x) == a);
        }
    }
    const StateId Z;
    SID_CHECK(Z.size() == 0);
    std::cout << "[state_id] random_test: OK (" << cases << " cases)\n";
}

}  // namespace sid_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace sid_demo {

static std::string to_s(bool b) { return b ? "true" : "false"; }
template <class T, std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool>, int> = 0>
static std::string to_s(T v) { return std::to_string(v); }
static std::string to_s(const std::pair<int, int>& p) { return "(" + std::to_string(p.first) + ", " + std::to_string(p.second) + ")"; }
template <size_t N>
static std::string to_s(const std::array<int, N>& a) {
    std::string s = "(";
    for (size_t i = 0; i < N; ++i) s += (i ? ", " : "") + std::to_string(a[i]);
    return s + ")";
}
static std::string to_s(const std::vector<int>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + std::to_string(v[i]);
    return s + "}";
}
static std::string pad(const std::string& s, size_t w) {
    size_t cols = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = (unsigned char)s[i];
        if ((c & 0xC0) == 0x80) continue;
        cols += (c >= 0xE0 ? 2 : 1);
    }
    return cols >= w ? s : s + std::string(w - cols, ' ');
}
static void banner(const std::string& s) {
    std::cout << "\n=====================================================================\n"
              << "  " << s << "\n"
              << "=====================================================================\n";
}
static void api(const std::string& usage, const std::string& ret, const std::string& what, const std::string& cost) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << usage << "\n"
              << "    返り値 : " << ret << "\n"
              << "    内容   : " << what << "\n"
              << "    計算量 : " << cost << "\n";
}
static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
static void stmt(const std::string& s) { code(s + ";"); }
template <class V>
static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 28) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 46) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) sid_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    sid_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("StateId 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  StateId S(A, B[, C[, D]]);             各成分の個数を指定（a は 0 ~ A-1 ...）\n"
              << "  S.size()                               状態の数\n"
              << "  S.id(a, b[, c[, d]])                   状態の番号\n"
              << "  S.get(x) / S.get3(x) / S.get4(x)       番号 → (a, b) / (a, b, c) / (a, b, c, d)\n"
              << "  S.first(x)                             1 つ目の成分だけ\n"
              << "  S.in(a, b[, c[, d]])                   範囲内か\n"
              << "  S.ids_of(a)                            1 つ目の成分が a の状態すべて\n";

    api("【1】 StateId S(N, K);  S.id(v, k);  S.get(x);  S.size()", "int / pair<int,int>",
        "(頂点 v, 状態 k) と番号 x = v * K + k を行き来する。auto [v, k] = S.get(x); で受け取る", "O(1)");
    RUN(StateId S(4, 3));
    SHOW(S.size(), "4 頂点 × 状態 3 種類 = 12 状態");
    SHOW(S.id(2, 1), "(頂点 2, 状態 1) = 2 * 3 + 1");
    SHOW(S.get(7), "番号 7 は (頂点 2, 状態 1)");
    SHOW(S.first(7), "頂点だけなら first");
    SHOW(S.in(2, 3), "状態 3 は範囲外（0 ~ 2）");

    api("【2】 StateId S(H*W, 4);  S(A, B, C);  S.get3(x)  /  S(A, B, C, D);  S.get4(x)", "array<int, 3> / array<int, 4>",
        "3 つ・4 つの組にも番号を付けられる。例: (マス, 向き)、(マス, 向き, 残り回数)", "O(1)");
    RUN(StateId T(6, 4, 3));
    SHOW(T.id(5, 2, 1), "(マス 5, 向き 2, 残り 1)");
    SHOW(T.get3(T.id(5, 2, 1)), "auto [cell, dir, rest] = T.get3(x);");
    RUN(StateId U(2, 3, 4, 5));
    SHOW(U.get4(U.id(1, 2, 3, 4)), "4 つの組");

    api("【3】 S.ids_of(a)", "vector<int>", "1 つ目の成分が a の状態すべて（終点のどの状態でもよいときに最小を取る）", "O(B*C*D)");
    SHOW(S.ids_of(3), "頂点 3 の状態 0, 1, 2");

    head("【4】 例: クーポン問題（K 枚まで辺のコストを 0 にできる。頂点 0 → 3 の最小コスト）");
    note("状態 = (頂点, 使ったクーポンの枚数)。辺を普通に通る（コスト w）か、クーポンを使う（コスト 0、枚数 +1）。");
    note("本番では 02_dijkstra の dijkstra_core を使う。ここでは同じことを priority_queue で書いて見せる。");
    code("// 本番の書き方:");
    code("//   StateId S(N, K + 1);");
    code("//   auto r = dijkstra_core(S.size(), S.id(0, 0), [&](int x, auto&& push) {");
    code("//       auto [v, k] = S.get(x);");
    code("//       for (auto& e : G[v]) {");
    code("//           push(S.id(e.to, k), e.w);                 // 普通に通る");
    code("//           if (k < K) push(S.id(e.to, k + 1), 0);    // クーポンを使う");
    code("//       }");
    code("//   });");
    code("//   long long ans = DIST_INF;  for (int x : S.ids_of(N - 1)) ans = min(ans, r.dist[x]);");
    {
        // 頂点 0-1 (10), 1-3 (10), 0-2 (3), 2-3 (30)。K = 1
        const int N = 4, K = 1;
        std::vector<std::vector<std::pair<int, long long>>> G(N);
        auto add = [&](int u, int v, long long w) { G[u].push_back({v, w}), G[v].push_back({u, w}); };
        add(0, 1, 10), add(1, 3, 10), add(0, 2, 3), add(2, 3, 30);
        const StateId SK(N, K + 1);
        std::vector<long long> dist(SK.size(), (long long)1e18);
        std::priority_queue<std::pair<long long, int>, std::vector<std::pair<long long, int>>, std::greater<>> pq;
        dist[SK.id(0, 0)] = 0, pq.push({0, SK.id(0, 0)});
        while (!pq.empty()) {
            const long long d = pq.top().first;
            const int x = pq.top().second;
            pq.pop();
            if (d != dist[x]) continue;
            const auto [v, k] = SK.get(x);
            auto push = [&](int y, long long w) {
                if (d + w < dist[y]) dist[y] = d + w, pq.push({dist[y], y});
            };
            for (auto [to, w] : G[v]) {
                push(SK.id(to, k), w);
                if (k < K) push(SK.id(to, k + 1), 0);
            }
        }
        long long ans = (long long)1e18;
        for (int x : SK.ids_of(N - 1)) ans = std::min(ans, dist[x]);
        note("グラフ: 0-1 (10), 1-3 (10), 0-2 (3), 2-3 (30)、クーポン 1 枚");
        note("クーポンなしなら 20、2-3 にクーポンを使うと 3 + 0 = 3");
        std::cout << "    答え = " << ans << "\n";
    }

    head("状態の例（どれも StateId で番号を付けて、探索の本体に隣接を渡す）");
    note("鍵を持っているか / スイッチの状態     : StateId S(N, 2)");
    note("残りの回数（クーポン・壁を壊す回数）  : StateId S(N, K + 1)");
    note("向き（グリッド）                      : StateId S(H * W, 4)。向きを変えるのにコスト 1 なら 0-1 BFS");
    note("時刻 mod K / 偶奇                     : StateId S(N, K) / StateId S(N, 2)");
    note("訪れた頂点の集合（bit）               : StateId S(N, 1 << K)（K が小さいとき）");
    note("類題: ABC132 E（3 歩ごと = 歩数 mod 3）、ABC246 E（斜めの向き）、ABC395 E（辺の向きの反転 = 偶奇）");
}

}  // namespace sid_demo


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
    sid_test::random_test();
    sid_demo::demo();
    return 0;
}
