// =============================================================================
//  全探索の列挙（直積・重複組合せ・組合せ・順列 nPk・和の分け方・整数の分割）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  「長さ N の列を全部試す」ときの多重ループ・再帰を、関数 1 つで書けるようにする。どれも辞書順に f(列) を呼ぶ。
//  f が bool を返すときは、false を返したらそこで止める（見つかったら終わり、に使う）
//    例: ABC165 C（Many Requirements）… 1 ≤ A_1 ≤ … ≤ A_N ≤ M を全部 → for_each_multichoose(N, 1, M, f)
//
//  機能一覧（f は const vector<int>& を受け取る。返り値は呼んだ回数）                          個数
//    for_each_product(sizes, f)              a[i] ∈ [0, sizes[i]) の組を全部（直積）                 Π sizes[i]
//    for_each_product(n, k, f)               a[i] ∈ [0, k) の長さ n の列を全部（k 進数の数え上げ）    k^n
//    for_each_multichoose(n, lo, hi, f)      lo ≤ a[0] ≤ a[1] ≤ … ≤ a[n−1] ≤ hi（重複組合せ）        C(hi − lo + n, n)
//    for_each_combination(m, k, f)           0〜m−1 から k 個を小さい順に（a[0] < a[1] < …）         C(m, k)
//    for_each_permutation_k(n, k, f)         0〜n−1 から違う k 個を並べる（順番も区別）              n! / (n − k)!
//    for_each_composition(total, parts, f)   a[0] + … + a[parts−1] = total（a[i] ≥ 0。玉を箱に分ける）  C(total + parts − 1, parts − 1)
//    for_each_partition(n, f)                n を大きい順の正の整数の和に（4 = 3 + 1 = 2 + 2 など）     分割数 p(n)
//
//  注意: 個数がすぐ大きくなる。10^7 くらいまでにする（f の中身が軽ければ 10^8 も可）。n = 0 の列（空の列）も 1 回呼ぶ
//  速さの目安: for_each_product(10, 5, f)（約 10^7 個）約 0.015 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    enumerate_test::random_test()  … 再帰で作った愚直な列挙・個数の公式との突き合わせ
//    enumerate_test::benchmark()    … 速度計測
//    enumerate_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace enumerate_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// f を呼び、続けるなら true（f が void なら常に true）
template <class F>
bool call(F& f, const std::vector<int>& v) {
    if constexpr (std::is_same_v<decltype(f(v)), bool>) return f(v);
    else {
        f(v);
        return true;
    }
}
}  // namespace enumerate_internal

template <class F>
long long for_each_product(const std::vector<int>& sizes, F f) {
    const int n = (int)sizes.size();
    for (int s : sizes)
        if (s <= 0) return 0;
    std::vector<int> a(n, 0);
    long long cnt = 0;
    while (true) {
        ++cnt;
        if (!enumerate_internal::call(f, a)) return cnt;
        int i = n - 1;
        while (i >= 0 && a[i] + 1 == sizes[i]) a[i] = 0, --i;
        if (i < 0) return cnt;
        ++a[i];
    }
}
template <class F>
long long for_each_product(int n, int k, F f) {
    if (n < 0) enumerate_internal::fail("for_each_product", "n は 0 以上");
    return for_each_product(std::vector<int>(n, k), f);
}
template <class F>
long long for_each_multichoose(int n, int lo, int hi, F f) {
    if (n < 0) enumerate_internal::fail("for_each_multichoose", "n は 0 以上");
    if (n > 0 && lo > hi) return 0;
    std::vector<int> a(n, lo);
    long long cnt = 0;
    while (true) {
        ++cnt;
        if (!enumerate_internal::call(f, a)) return cnt;
        int i = n - 1;
        while (i >= 0 && a[i] == hi) --i;
        if (i < 0) return cnt;
        const int v = a[i] + 1;
        for (int j = i; j < n; ++j) a[j] = v;
    }
}
template <class F>
long long for_each_combination(int m, int k, F f) {
    if (k < 0 || m < 0) enumerate_internal::fail("for_each_combination", "m, k は 0 以上");
    if (k > m) return 0;
    std::vector<int> a(k);
    for (int i = 0; i < k; ++i) a[i] = i;
    long long cnt = 0;
    while (true) {
        ++cnt;
        if (!enumerate_internal::call(f, a)) return cnt;
        int i = k - 1;
        while (i >= 0 && a[i] == m - k + i) --i;
        if (i < 0) return cnt;
        ++a[i];
        for (int j = i + 1; j < k; ++j) a[j] = a[j - 1] + 1;
    }
}
template <class F>
long long for_each_permutation_k(int n, int k, F f) {
    if (k < 0 || n < 0) enumerate_internal::fail("for_each_permutation_k", "n, k は 0 以上");
    if (k > n) return 0;
    std::vector<int> a(k);
    std::vector<char> used(n, 0);
    long long cnt = 0;
    bool stop = false;
    auto rec = [&](auto&& self, int i) -> void {
        if (stop) return;
        if (i == k) {
            ++cnt;
            if (!enumerate_internal::call(f, a)) stop = true;
            return;
        }
        for (int v = 0; v < n && !stop; ++v) {
            if (used[v]) continue;
            used[v] = 1, a[i] = v;
            self(self, i + 1);
            used[v] = 0;
        }
    };
    rec(rec, 0);
    return cnt;
}
template <class F>
long long for_each_composition(int total, int parts, F f) {
    if (total < 0 || parts < 0) enumerate_internal::fail("for_each_composition", "total, parts は 0 以上");
    if (parts == 0) {
        if (total != 0) return 0;
        std::vector<int> e;
        enumerate_internal::call(f, e);
        return 1;
    }
    std::vector<int> a(parts, 0);
    a[parts - 1] = total;  // 辞書順の最初は (0, …, 0, total)
    long long cnt = 0;
    while (true) {
        ++cnt;
        if (!enumerate_internal::call(f, a)) return cnt;
        // 最後の要素以外で、右に残っている分から 1 つ持ってこられる一番右の位置を増やす
        int i = parts - 2;
        while (i >= 0) {
            int rest = 0;
            for (int j = i + 1; j < parts; ++j) rest += a[j];
            if (rest > 0) break;
            --i;
        }
        if (i < 0) return cnt;
        int rest = -1;
        for (int j = i + 1; j < parts; ++j) rest += a[j];
        ++a[i];
        for (int j = i + 1; j < parts; ++j) a[j] = 0;
        a[parts - 1] = rest;
    }
}
template <class F>
long long for_each_partition(int n, F f) {
    if (n < 0) enumerate_internal::fail("for_each_partition", "n は 0 以上");
    std::vector<int> a;
    long long cnt = 0;
    bool stop = false;
    // 辞書順（大きい部分から）：a[0] ≥ a[1] ≥ …
    auto rec = [&](auto&& self, int rest, int mx) -> void {
        if (stop) return;
        if (rest == 0) {
            ++cnt;
            if (!enumerate_internal::call(f, a)) stop = true;
            return;
        }
        for (int v = 1; v <= std::min(rest, mx) && !stop; ++v) {
            a.push_back(v);
            self(self, rest - v, v);
            a.pop_back();
        }
    };
    rec(rec, n, n);
    return cnt;
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace enumerate_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define EN_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[search_enumerate] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
using VV = std::vector<std::vector<int>>;

template <class G>
static VV collect(G g) {
    VV out;
    g([&](const std::vector<int>& v) { out.push_back(v); });
    return out;
}
// 長さ n、各値 [0, k) の全部の列（再帰）から条件に合うものを辞書順で
template <class P>
static VV brute(int n, int lo, int hi, P ok) {
    VV out;
    std::vector<int> a(n);
    auto rec = [&](auto&& self, int i) -> void {
        if (i == n) {
            if (ok(a)) out.push_back(a);
            return;
        }
        for (int v = lo; v <= hi; ++v) a[i] = v, self(self, i + 1);
    };
    rec(rec, 0);
    return out;
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        // 直積
        {
            const int n = (int)rnd(rng, 0, 4);
            std::vector<int> sz(n);
            for (auto& s : sz) s = (int)rnd(rng, 0, 4);
            int mx = 0;
            for (int s : sz) mx = std::max(mx, s);
            const auto got = collect([&](auto f) { for_each_product(sz, f); });
            const auto exp = brute(n, 0, std::max(0, mx - 1), [&](const std::vector<int>& a) {
                for (int i = 0; i < n; ++i)
                    if (a[i] >= sz[i]) return false;
                return true;
            });
            bool any0 = false;
            for (int s : sz) any0 = any0 || s == 0;
            EN_CHECK(got == (any0 ? VV{} : exp));
            const int k = (int)rnd(rng, 1, 3);
            EN_CHECK(collect([&](auto f) { for_each_product(n, k, f); }) == brute(n, 0, k - 1, [](const std::vector<int>&) { return true; }));
            cases += 2;
        }
        // 重複組合せ・組合せ・nPk
        {
            const int n = (int)rnd(rng, 0, 4), lo = (int)rnd(rng, -2, 2), hi = lo + (int)rnd(rng, -1, 3);
            EN_CHECK(collect([&](auto f) { for_each_multichoose(n, lo, hi, f); }) ==
                     brute(n, lo, std::max(lo, hi), [&](const std::vector<int>& a) {
                         if (hi < lo && n > 0) return false;
                         for (int i = 0; i + 1 < n; ++i)
                             if (a[i] > a[i + 1]) return false;
                         for (int x : a)
                             if (x > hi) return false;
                         return true;
                     }));
            const int m = (int)rnd(rng, 0, 6), k = (int)rnd(rng, 0, 4);
            EN_CHECK(collect([&](auto f) { for_each_combination(m, k, f); }) ==
                     (k > m ? VV{} : brute(k, 0, std::max(0, m - 1), [&](const std::vector<int>& a) {
                         for (int i = 0; i + 1 < k; ++i)
                             if (a[i] >= a[i + 1]) return false;
                         return m > 0 || k == 0;
                     })));
            EN_CHECK(collect([&](auto f) { for_each_permutation_k(m, k, f); }) ==
                     (k > m ? VV{} : brute(k, 0, std::max(0, m - 1), [&](const std::vector<int>& a) {
                         for (int i = 0; i < k; ++i)
                             for (int j = i + 1; j < k; ++j)
                                 if (a[i] == a[j]) return false;
                         return m > 0 || k == 0;
                     })));
            cases += 3;
        }
        // 和の分け方・整数の分割
        {
            const int total = (int)rnd(rng, 0, 5), parts = (int)rnd(rng, 0, 4);
            EN_CHECK(collect([&](auto f) { for_each_composition(total, parts, f); }) == brute(parts, 0, total, [&](const std::vector<int>& a) {
                         int s = 0;
                         for (int x : a) s += x;
                         return s == total;
                     }));
            const int n = (int)rnd(rng, 0, 9);
            VV exp;  // 正の整数の並べ方（2^(n−1) 通り。区切りの位置を bit で）のうち、大きい順のもの
            if (n == 0) exp.push_back({});
            for (int m = 0; n > 0 && m < (1 << (n - 1)); ++m) {
                std::vector<int> a{1};
                for (int i = 0; i < n - 1; ++i) {
                    if (m >> i & 1) a.push_back(1);
                    else ++a.back();
                }
                if (std::is_sorted(a.rbegin(), a.rend())) exp.push_back(a);
            }
            std::sort(exp.begin(), exp.end(), [](const std::vector<int>& x, const std::vector<int>& y) {
                return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
            });
            EN_CHECK(collect([&](auto f) { for_each_partition(n, f); }) == exp);
            cases += 2;
        }
    }
    // 途中で止める（bool を返す）
    {
        int seen = 0;
        const ll c = for_each_product(3, 10, [&](const std::vector<int>& a) {
            ++seen;
            return !(a[0] == 1 && a[1] == 2 && a[2] == 3);
        });
        EN_CHECK(c == 124 && seen == 124);
        EN_CHECK(for_each_partition(10, [](const std::vector<int>&) {}) == 42);
        EN_CHECK(for_each_composition(4, 3, [](const std::vector<int>&) {}) == 15);
        cases += 3;
    }
    std::cout << "[search_enumerate] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[search_enumerate] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = TEST_SCALE == 1 ? 8 : 10;
    ll chk = 0;
    auto t0 = std::chrono::steady_clock::now();
    const ll c1 = for_each_product(n, 5, [&](const std::vector<int>& a) { chk += a[0] + a[n - 1]; });
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    const ll c2 = for_each_multichoose(10, 1, TEST_SCALE == 1 ? 8 : 15, [&](const std::vector<int>& a) { chk += a[5]; });
    const ll t2 = elapsed_ms(t0);
    std::cout << "[search_enumerate] benchmark: for_each_product(" << n << ", 5)（" << c1 << " 個） " << t1 << " ms / for_each_multichoose（" << c2
              << " 個） " << t2 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace enumerate_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace enumerate_demo {

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
#define SHOW(expr, comment) enumerate_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    enumerate_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static std::string list_all(const std::function<void(std::function<void(const std::vector<int>&)>)>& g) {
    std::string s;
    g([&](const std::vector<int>& v) { s += to_s(v) + " "; });
    return s;
}

static void demo() {
    banner("全探索の列挙（search_enumerate） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  for_each_product(sizes, f) / for_each_product(n, k, f) / for_each_multichoose / for_each_combination\n"
              << "  for_each_permutation_k / for_each_composition / for_each_partition   （f が false を返すと止まる。返り値は呼んだ回数）\n";

    api("【1】 for_each_product(sizes, f) / for_each_product(n, k, f)", "long long（呼んだ回数）", "a[i] ∈ [0, sizes[i]) の全部 / [0, k)^n の全部（辞書順）",
        "O(個数 × n)");
    note(list_all([](auto f) { for_each_product(std::vector<int>{2, 3}, f); }));
    SHOW(for_each_product(3, 4, [](const std::vector<int>&) {}), "4^3");

    api("【2】 for_each_multichoose(n, lo, hi, f)", "long long", "lo ≤ a[0] ≤ … ≤ a[n−1] ≤ hi を全部（ABC165 C の形）", "O(個数 × n)");
    note(list_all([](auto f) { for_each_multichoose(2, 1, 3, f); }));
    {
        code("// 1 ≤ A1 ≤ A2 ≤ A3 ≤ 4 で、条件 A[b] − A[a] == c のとき d 点。点の和の最大（自作の例）");
        const std::vector<std::vector<int>> req = {{0, 2, 3, 100}, {0, 1, 1, 10}, {1, 2, 2, 5}};  // {a, b, c, d}
        long long best = 0;
        for_each_multichoose(3, 1, 4, [&](const std::vector<int>& A) {
            long long s = 0;
            for (const auto& r : req)
                if (A[r[1]] - A[r[0]] == r[2]) s += r[3];
            best = std::max(best, s);
        });
        note("→ 最大 " + to_s(best) + "（例: A = {1, 2, 4} で 100 + 10 + 5）");
    }

    api("【3】 for_each_combination(m, k, f) / for_each_permutation_k(n, k, f)", "long long", "0〜m−1 から k 個（小さい順）/ 違う k 個を並べる（順番も区別）",
        "O(個数 × k)");
    note(list_all([](auto f) { for_each_combination(4, 2, f); }));
    SHOW(for_each_permutation_k(5, 3, [](const std::vector<int>&) {}), "5 × 4 × 3");

    api("【4】 for_each_composition(total, parts, f) / for_each_partition(n, f)", "long long",
        "a[0] + … + a[parts−1] = total（0 以上）/ n を大きい順の正の整数の和に", "O(個数 × 長さ)");
    note(list_all([](auto f) { for_each_composition(2, 3, f); }));
    note(list_all([](auto f) { for_each_partition(5, f); }));

    api("【5】 f が bool を返すと途中で止められる", "—", "false を返したところで終わる（返り値はそこまでに呼んだ回数）", "—");
    SHOW(for_each_product(2, 10, [](const std::vector<int>& a) { return a[0] * 10 + a[1] != 37; }), "37 で止まる（00 から 38 個目）");
}

}  // namespace enumerate_demo

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
    enumerate_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    enumerate_test::benchmark();    // 速度計測
    enumerate_demo::demo();         // 使い方デモ
    return 0;
}
