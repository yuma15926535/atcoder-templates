// =============================================================================
//  隣接交換の最小回数（目的の並びにする・並べ替える・回文にする・好きな 2 つの交換）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  隣どうしの入れ替えで a を b にする最小回数は、「a の各要素が b の何番目に行くか」の列の転倒数。
//  同じ値が何個もあるときは、同じ値どうしは順番を変えずに対応させるのが最小（a の k 個目の x → b の k 個目の x）
//    例: 「文字列 S を T にする隣接交換の最小回数」「0 / 1 の列を 0 → 1 の順に並べる」「回文にする最小回数」
//
//  機能一覧（列は vector<T>（T は == と < で比べられる型）か string）                            計算量
//    min_adjacent_swaps(a, b)            隣接交換で a を b にする最小回数（並べ替えて同じでなければ −1）     O(N log N)
//    min_adjacent_swaps_sort(a)          a を小さい順（同じ値は元の順）に並べる最小回数（= 転倒数）         O(N log N)
//    target_positions(a, b)              a[i] が行く b の位置 p[i]（min_adjacent_swaps は p の転倒数）       O(N log N)
//    min_adjacent_swaps_palindrome(s)    隣接交換で s を回文にする最小回数（できなければ −1）               O(N log N)
//    min_swaps_any(a, b)                 好きな 2 つを入れ替える操作で a を b にする最小回数（a は相異なる値。−1 は無理）  O(N log N)
//
//  注意: 回数は long long（最大 N(N − 1)/2）。0 / 1 の列を 1 を K 個つなげる回数は gather_cost.cpp の gather_consecutive
//  速さの目安: N = 2×10^5 で min_adjacent_swaps 約 0.03 秒、min_adjacent_swaps_palindrome 約 0.02 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    adjswap_test::random_test()  … 幅優先探索の愚直解との突き合わせ
//    adjswap_test::benchmark()    … 速度計測
//    adjswap_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace adjswap_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
inline long long inversions(const std::vector<int>& p) {  // p は 0〜n−1 の値（重複なし）
    const int n = (int)p.size();
    std::vector<int> bit(n + 1, 0);
    long long inv = 0;
    for (int i = n - 1; i >= 0; --i) {
        for (int k = p[i]; k > 0; k -= k & -k) inv += bit[k];  // 右にある p[i] 未満の数
        for (int k = p[i] + 1; k <= n; k += k & -k) ++bit[k];
    }
    return inv;
}
// 安定な並べ替えの順番（値が同じなら元の順）
template <class V>
std::vector<int> stable_order(const V& a) {
    std::vector<int> id(a.size());
    for (size_t i = 0; i < a.size(); ++i) id[i] = (int)i;
    std::stable_sort(id.begin(), id.end(), [&](int x, int y) { return a[x] < a[y]; });
    return id;
}
}  // namespace adjswap_internal

// a[i] が行く b の位置（同じ値は順番どおり）。並べ替えて同じでなければ空
template <class V>
std::vector<int> target_positions(const V& a, const V& b) {
    if (a.size() != b.size()) return {};
    const auto oa = adjswap_internal::stable_order(a), ob = adjswap_internal::stable_order(b);
    std::vector<int> p(a.size());
    for (size_t k = 0; k < a.size(); ++k) {
        if (a[oa[k]] < b[ob[k]] || b[ob[k]] < a[oa[k]]) return {};
        p[oa[k]] = ob[k];  // 同じ値の中で k 番目どうし
    }
    return p;
}
template <class V>
long long min_adjacent_swaps(const V& a, const V& b) {
    if (a.size() != b.size()) return -1;
    const auto p = target_positions(a, b);
    if (p.size() != a.size()) return -1;
    return adjswap_internal::inversions(p);
}
template <class V>
long long min_adjacent_swaps_sort(const V& a) {
    const auto o = adjswap_internal::stable_order(a);
    std::vector<int> p(a.size());
    for (size_t k = 0; k < o.size(); ++k) p[o[k]] = (int)k;
    return adjswap_internal::inversions(p);
}
// 左から見て、まだ使っていない文字を、同じ文字の一番右のものと組にする（左は左端から、右は右端から詰める）
template <class V>
long long min_adjacent_swaps_palindrome(const V& s) {
    const int n = (int)s.size();
    const auto o = adjswap_internal::stable_order(s);
    std::vector<int> grp(n);  // 値の番号
    int g = -1;
    for (int k = 0; k < n; ++k) {
        if (k == 0 || s[o[k - 1]] < s[o[k]]) ++g;
        grp[o[k]] = g;
    }
    std::vector<std::deque<int>> pos(g + 1);
    for (int i = 0; i < n; ++i) pos[grp[i]].push_back(i);
    int odd = 0;
    for (const auto& d : pos) odd += (int)(d.size() % 2);
    if (odd > n % 2) return -1;
    std::vector<int> t(n, -1);
    std::vector<char> used(n, 0);
    int l = 0, r = n - 1;
    for (int i = 0; i < n; ++i) {
        if (used[i]) continue;
        auto& d = pos[grp[i]];
        if (d.size() == 1) {  // 真ん中に行く文字
            t[i] = n / 2, used[i] = 1;
            d.pop_front();
            continue;
        }
        const int j = d.back();
        d.pop_front(), d.pop_back();
        t[i] = l++, t[j] = r--;
        used[i] = used[j] = 1;
    }
    return adjswap_internal::inversions(t);
}
template <class V>
long long min_swaps_any(const V& a, const V& b) {
    const auto p = target_positions(a, b);
    if (p.size() != a.size()) return -1;
    {  // 相異なるか
        const auto o = adjswap_internal::stable_order(a);
        for (size_t q = 0; q + 1 < o.size(); ++q)
            if (!(a[o[q]] < a[o[q + 1]])) adjswap_internal::fail("min_swaps_any", "a の値は相異なるようにする");
    }
    const int n = (int)p.size();
    std::vector<char> seen(n, 0);
    long long c = 0;
    for (int i = 0; i < n; ++i) {
        if (seen[i]) continue;
        ++c;
        for (int x = i; !seen[x]; x = p[x]) seen[x] = 1;
    }
    return n - c;  // n − (巡回の数)
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace adjswap_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define AS_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[gather_adjacent_swaps] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

// 幅優先探索：a から（隣接 / 任意の）交換で各並びへの最小回数
static std::map<std::string, int> bfs(const std::string& a, bool adjacent) {
    std::map<std::string, int> d;
    std::queue<std::string> q;
    d[a] = 0, q.push(a);
    while (!q.empty()) {
        const std::string s = q.front();
        q.pop();
        for (size_t i = 0; i < s.size(); ++i)
            for (size_t j = i + 1; j < s.size(); ++j) {
                if (adjacent && j != i + 1) break;
                std::string t = s;
                std::swap(t[i], t[j]);
                if (!d.count(t)) d[t] = d[s] + 1, q.push(t);
            }
    }
    return d;
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 7);
        const int sigma = (int)rnd(rng, 1, 4);
        std::string a(n, 'a');
        for (auto& c : a) c = (char)('a' + rnd(rng, 0, sigma - 1));
        const auto d = bfs(a, true);
        // a → b（並べ替えた b と、関係ない b）
        std::string b = a;
        std::shuffle(b.begin(), b.end(), rng);
        if (rnd(rng, 0, 4) == 0 && n > 0) b[0] = (char)('a' + rnd(rng, 0, sigma));
        AS_CHECK(min_adjacent_swaps(a, b) == (d.count(b) ? d.at(b) : -1));
        std::string s = a;
        std::sort(s.begin(), s.end());
        AS_CHECK(min_adjacent_swaps_sort(a) == d.at(s));
        // 回文
        int best = -1;
        for (const auto& [t, c] : d) {
            std::string rv(t.rbegin(), t.rend());
            if (rv == t && (best < 0 || c < best)) best = c;
        }
        AS_CHECK(min_adjacent_swaps_palindrome(a) == best);
        // 好きな 2 つの交換（相異なる値）
        std::string u;
        for (int i = 0; i < n; ++i) u += (char)('a' + i);
        std::shuffle(u.begin(), u.end(), rng);
        std::string v = u;
        std::shuffle(v.begin(), v.end(), rng);
        AS_CHECK(min_swaps_any(u, v) == bfs(u, false).at(v));
        // vector<int> でも
        std::vector<int> va(a.begin(), a.end()), vb(b.begin(), b.end());
        AS_CHECK(min_adjacent_swaps(va, vb) == min_adjacent_swaps(a, b));
        cases += 5;
    }
    std::cout << "[gather_adjacent_swaps] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[gather_adjacent_swaps] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    std::string a(N, 'a');
    for (auto& c : a) c = (char)('a' + rnd(rng, 0, 25));
    std::string b = a;
    std::shuffle(b.begin(), b.end(), rng);
    std::string p = a.substr(0, N / 2);
    p += std::string(p.rbegin(), p.rend());
    std::shuffle(p.begin(), p.end(), rng);
    auto t0 = std::chrono::steady_clock::now();
    ll chk = min_adjacent_swaps(a, b);
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    chk += min_adjacent_swaps_palindrome(p);
    const ll t2 = elapsed_ms(t0);
    std::cout << "[gather_adjacent_swaps] benchmark (N = " << N << "): min_adjacent_swaps " << t1 << " ms / min_adjacent_swaps_palindrome " << t2
              << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace adjswap_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace adjswap_demo {

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
#define SHOW(expr, comment) adjswap_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    adjswap_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("隣接交換の最小回数（gather_adjacent_swaps） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  min_adjacent_swaps(a, b) / min_adjacent_swaps_sort(a) / target_positions(a, b)\n"
              << "  min_adjacent_swaps_palindrome(s) / min_swaps_any(a, b)\n";

    api("【1】 min_adjacent_swaps(a, b) / target_positions(a, b)", "long long / vector<int>",
        "隣接交換で a を b にする最小回数（できなければ −1）/ a[i] が行く b の位置（同じ文字は順番どおり）", "O(N log N)");
    SHOW(min_adjacent_swaps(std::string("abcab"), std::string("ababc")), "");
    SHOW(target_positions(std::string("abcab"), std::string("ababc")), "この列の転倒数が答え");
    SHOW(min_adjacent_swaps(std::string("abc"), std::string("abd")), "文字が違う → −1");

    api("【2】 min_adjacent_swaps_sort(a)", "long long", "小さい順に並べる最小回数（転倒数）。0 / 1 の列なら「1 の後ろにある 0」の数の和", "O(N log N)");
    SHOW(min_adjacent_swaps_sort(std::vector<int>{1, 0, 1, 0, 0}), "1 の後ろの 0：3 + 2");
    SHOW(min_adjacent_swaps_sort(std::vector<int>{3, 1, 2}), "");

    api("【3】 min_adjacent_swaps_palindrome(s)", "long long", "隣接交換で回文にする最小回数（できなければ −1）", "O(N log N)");
    SHOW(min_adjacent_swaps_palindrome(std::string("aabb")), "abba にする");
    SHOW(min_adjacent_swaps_palindrome(std::string("ntiin")), "nitin にする");
    SHOW(min_adjacent_swaps_palindrome(std::string("abc")), "回文にできない");

    api("【4】 min_swaps_any(a, b)", "long long", "好きな 2 つを入れ替える操作で a を b にする最小回数 = N − (巡回の数)（a は相異なる値）", "O(N log N)");
    SHOW(min_swaps_any(std::vector<int>{2, 0, 1, 3}, std::vector<int>{0, 1, 2, 3}), "2 → 0 → 1 の巡回 1 つ（長さ 3）で 2 回");
}

}  // namespace adjswap_demo

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
    adjswap_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    adjswap_test::benchmark();    // 速度計測
    adjswap_demo::demo();         // 使い方デモ
    return 0;
}
