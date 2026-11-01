// =============================================================================
//  線形計画（単体法：一次不等式 A x ≤ b、x ≥ 0 のもとで c·x を最大化）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  実行可能領域（凸多面体）の頂点を、目的関数がよくなる向きにたどる。b に負があれば、まず実行可能な頂点を探す（2 段階法）。
//  循環しないように、同じなら番号の小さい変数を選ぶ（Bland の規則）。実数（double / long double）で計算する
//    例: 「材料の制約のもとで利益を最大に」「最小費用流・最大流を LP で書いたもの」、双対（min b·y、Aᵀy ≥ c、y ≥ 0）も同じ値
//
//  機能一覧（m 本の制約・n 個の変数）                                                           計算量
//    SimplexLP<T> lp(A, b, c)            A は m × n、b は長さ m、c は長さ n（T は double か long double）    O(m n)
//    lp.solve()                          LPResult<T> {status, value, x}                                    多くは速い（最悪は指数）
//                                          status：LPStatus::Optimal / Infeasible（解なし）/ Unbounded（いくらでも大きい）
//    最小化は c を −1 倍、等式は ≤ と ≥（−1 倍して ≤）の 2 本、x が負でもよい変数は x = x⁺ − x⁻ に分ける
//
//  注意: 誤差の許し幅 eps = 1e-9（long double なら小さくしてもよい）。答えが整数のはずなら四捨五入する
//        m, n ≤ 数百くらい。整数計画（x が整数）は解けない
//  速さの目安: m = n = 300 のランダムな問題（係数が正）で 約 0.005 秒。ピボット 1 回が O(m n)、回数は問題しだい
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    simplex_test::random_test()  … 頂点を全部調べる愚直解（変数 3 個まで）との突き合わせ
//    simplex_test::benchmark()    … 速度計測
//    simplex_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
enum class LPStatus { Optimal, Infeasible, Unbounded };
template <class T>
struct LPResult {
    LPStatus status;
    T value;
    std::vector<T> x;
};

template <class T = double>
struct SimplexLP {
    SimplexLP(const std::vector<std::vector<T>>& A, const std::vector<T>& b, const std::vector<T>& c, T eps = (T)1e-9)
        : m_((int)b.size()), n_((int)c.size()), eps_(eps), N_(n_ + 1), B_(m_), D_(m_ + 2, std::vector<T>(n_ + 2, T(0))) {
        if ((int)A.size() != m_) fail("A の行の数と b の長さが違う");
        for (int i = 0; i < m_; ++i) {
            if ((int)A[i].size() != n_) fail("A の列の数と c の長さが違う");
            for (int j = 0; j < n_; ++j) D_[i][j] = A[i][j];
            B_[i] = n_ + i, D_[i][n_] = T(-1), D_[i][n_ + 1] = b[i];
        }
        for (int j = 0; j < n_; ++j) N_[j] = j, D_[m_][j] = -c[j];
        N_[n_] = -1, D_[m_ + 1][n_] = T(1);
    }
    LPResult<T> solve() {
        LPResult<T> res{LPStatus::Optimal, T(0), std::vector<T>(n_, T(0))};
        if (m_ > 0) {
            int r = 0;
            for (int i = 1; i < m_; ++i)
                if (D_[i][n_ + 1] < D_[r][n_ + 1]) r = i;
            if (D_[r][n_ + 1] < -eps_) {  // b に負がある → まず実行可能な頂点を探す
                pivot(r, n_);
                if (!simplex(2) || D_[m_ + 1][n_ + 1] < -eps_) {
                    res.status = LPStatus::Infeasible;
                    return res;
                }
                for (int i = 0; i < m_; ++i)
                    if (B_[i] == -1) {
                        int s = 0;
                        for (int j = 1; j <= n_; ++j) better(D_[i], j, s);
                        pivot(i, s);
                    }
            }
        }
        if (!simplex(1)) {
            res.status = LPStatus::Unbounded;
            return res;
        }
        for (int i = 0; i < m_; ++i)
            if (B_[i] < n_) res.x[B_[i]] = D_[i][n_ + 1];
        res.value = D_[m_][n_ + 1];
        return res;
    }

   private:
    int m_, n_;
    T eps_;
    std::vector<int> N_, B_;
    std::vector<std::vector<T>> D_;
    [[noreturn]] static void fail(const char* msg) {
        std::cerr << "SimplexLP: " << msg << "\n";
        std::abort();
    }
    void better(const std::vector<T>& X, int j, int& s) const {  // 値が小さい方、同じなら番号の小さい方（Bland）
        if (s == -1 || std::make_pair(X[j], N_[j]) < std::make_pair(X[s], N_[s])) s = j;
    }
    void pivot(int r, int s) {
        const T inv = T(1) / D_[r][s];
        for (int i = 0; i < m_ + 2; ++i)
            if (i != r && std::abs(D_[i][s]) > eps_) {
                const T inv2 = D_[i][s] * inv;
                for (int j = 0; j < n_ + 2; ++j) D_[i][j] -= D_[r][j] * inv2;
                D_[i][s] = D_[r][s] * inv2;
            }
        for (int j = 0; j < n_ + 2; ++j)
            if (j != s) D_[r][j] *= inv;
        for (int i = 0; i < m_ + 2; ++i)
            if (i != r) D_[i][s] *= -inv;
        D_[r][s] = inv;
        std::swap(B_[r], N_[s]);
    }
    bool simplex(int phase) {
        const int x = m_ + phase - 1;
        while (true) {
            int s = -1;
            for (int j = 0; j <= n_; ++j)
                if (N_[j] != -phase) better(D_[x], j, s);
            if (s == -1 || D_[x][s] >= -eps_) return true;
            int r = -1;
            for (int i = 0; i < m_; ++i) {
                if (D_[i][s] <= eps_) continue;
                if (r == -1 || std::make_pair(D_[i][n_ + 1] / D_[i][s], B_[i]) < std::make_pair(D_[r][n_ + 1] / D_[r][s], B_[r])) r = i;
            }
            if (r == -1) return false;
            pivot(r, s);
        }
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace simplex_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define SX_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[math_simplex] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
using VD = std::vector<double>;

// 愚直：n 本の等式（制約 + x_j ≥ 0 から選ぶ）を解いた点のうち、実行可能で c·x が最大
static std::optional<double> brute(const std::vector<VD>& A, const VD& b, const VD& c) {
    const int m = (int)b.size(), n = (int)c.size();
    std::vector<VD> rows;
    VD rhs;
    for (int i = 0; i < m; ++i) rows.push_back(A[i]), rhs.push_back(b[i]);
    for (int j = 0; j < n; ++j) {
        VD e(n, 0);
        e[j] = -1;
        rows.push_back(e), rhs.push_back(0);  // −x_j ≤ 0
    }
    const int R = (int)rows.size();
    std::optional<double> best;
    std::vector<int> pick(n);
    auto rec = [&](auto&& self, int k, int start) -> void {
        if (k == n) {
            std::vector<VD> M(n, VD(n + 1));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) M[i][j] = rows[pick[i]][j];
                M[i][n] = rhs[pick[i]];
            }
            for (int col = 0; col < n; ++col) {  // ガウスの消去
                int pv = -1;
                for (int i = col; i < n; ++i)
                    if (std::abs(M[i][col]) > 1e-9 && (pv < 0 || std::abs(M[i][col]) > std::abs(M[pv][col]))) pv = i;
                if (pv < 0) return;
                std::swap(M[col], M[pv]);
                for (int i = 0; i < n; ++i)
                    if (i != col) {
                        const double f = M[i][col] / M[col][col];
                        for (int j = col; j <= n; ++j) M[i][j] -= f * M[col][j];
                    }
            }
            VD x(n);
            for (int i = 0; i < n; ++i) x[i] = M[i][n] / M[i][i];
            for (int i = 0; i < R; ++i) {
                double s = 0;
                for (int j = 0; j < n; ++j) s += rows[i][j] * x[j];
                if (s > rhs[i] + 1e-7) return;
            }
            double v = 0;
            for (int j = 0; j < n; ++j) v += c[j] * x[j];
            if (!best || v > *best) best = v;
            return;
        }
        for (int i = start; i < R; ++i) pick[k] = i, self(self, k + 1, i + 1);
    };
    rec(rec, 0, 0);
    return best;
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 3), m = (int)rnd(rng, 0, 5);
        std::vector<VD> A(m, VD(n));
        VD b(m), c(n);
        for (auto& r : A)
            for (auto& v : r) v = (double)rnd(rng, -5, 5);
        for (auto& v : b) v = (double)rnd(rng, -5, 10);
        for (auto& v : c) v = (double)rnd(rng, -5, 5);
        // 必ず有界にする：各 x_j ≤ 10 を足す
        for (int j = 0; j < n; ++j) {
            VD e(n, 0);
            e[j] = 1;
            A.push_back(e), b.push_back((double)rnd(rng, 0, 10));
        }
        const auto r = SimplexLP<double>(A, b, c).solve();
        const auto br = brute(A, b, c);
        SX_CHECK((r.status == LPStatus::Optimal) == br.has_value());
        SX_CHECK(r.status != LPStatus::Unbounded);
        if (br) {
            SX_CHECK(std::abs(r.value - *br) < 1e-6);
            // x が本当に実行可能で、値が c·x
            double v = 0;
            for (int j = 0; j < n; ++j) {
                SX_CHECK(r.x[j] >= -1e-7);
                v += c[j] * r.x[j];
            }
            for (size_t i = 0; i < A.size(); ++i) {
                double s = 0;
                for (int j = 0; j < n; ++j) s += A[i][j] * r.x[j];
                SX_CHECK(s <= b[i] + 1e-6);
            }
            SX_CHECK(std::abs(v - r.value) < 1e-6);
        }
        ++cases;
    }
    // 有界でない・解なしの例
    {
        const auto u = SimplexLP<double>({{1, -1}}, {1}, {1, 0}).solve();  // x − y ≤ 1 で x を最大に → y を増やせばいくらでも
        SX_CHECK(u.status == LPStatus::Unbounded);
        const auto f = SimplexLP<double>({{1, 1}, {-1, -1}}, {1, -3}, {1, 1}).solve();  // x + y ≤ 1 かつ x + y ≥ 3
        SX_CHECK(f.status == LPStatus::Infeasible);
        const auto z = SimplexLP<long double>({}, {}, {-1, -2}).solve();  // 制約なし、係数が負 → 0
        SX_CHECK(z.status == LPStatus::Optimal && std::abs((double)z.value) < 1e-12);
        cases += 3;
    }
    std::cout << "[math_simplex] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[math_simplex] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int n = TEST_SCALE == 1 ? 60 : 300, m = n;
    std::vector<VD> A(m, VD(n));
    VD b(m), c(n);
    for (auto& r : A)
        for (auto& v : r) v = (double)rnd(rng, 1, 100);
    for (auto& v : b) v = (double)rnd(rng, 1000, 10000);
    for (auto& v : c) v = (double)rnd(rng, 1, 100);
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = SimplexLP<double>(A, b, c).solve();
    std::cout << "[math_simplex] benchmark (m = n = " << n << "): " << elapsed_ms(t0) << " ms  (値 " << std::fixed << std::setprecision(3) << r.value << ")\n"
              << std::defaultfloat;
}

}  // namespace simplex_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace simplex_demo {

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
#define SHOW(expr, comment) simplex_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    simplex_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static std::string st(LPStatus s) { return s == LPStatus::Optimal ? "Optimal" : (s == LPStatus::Infeasible ? "Infeasible" : "Unbounded"); }

static void demo() {
    banner("線形計画（math_simplex） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  SimplexLP<T> lp(A, b, c) → lp.solve() → LPResult {status, value, x}（A x ≤ b、x ≥ 0 で c·x を最大）\n";

    api("【1】 SimplexLP<double>(A, b, c).solve()", "LPResult<double>", "A x ≤ b、x ≥ 0 のもとで c·x を最大化", "多くは速い");
    code("// 製品 1 は材料 A を 1・B を 3、製品 2 は A を 2・B を 1 使う。A は 14、B は 18 まで。利益 3x + 2y を最大に");
    RUN(const auto r = SimplexLP<double>({{1, 2}, {3, 1}}, {14, 18}, {3, 2}).solve());
    SHOW(st(r.status), "");
    SHOW(r.value, "x = 4.4, y = 4.8 で 3 × 4.4 + 2 × 4.8 = 22.8");
    SHOW(r.x, "");

    api("【2】 解なし・有界でない", "LPStatus", "Infeasible（条件を全部満たす x が無い）/ Unbounded（いくらでも大きくできる）", "—");
    SHOW(st(SimplexLP<double>({{1, 1}, {-1, -1}}, {1, -3}, {1, 1}).solve().status), "x + y ≤ 1 と x + y ≥ 3");
    SHOW(st(SimplexLP<double>({{1, -1}}, {1}, {1, 0}).solve().status), "x − y ≤ 1 で x を最大に");

    api("【3】 最小化・等式・負になれる変数", "—", "c を −1 倍して最大化 → 値を −1 倍 / 等式は 2 本の不等式 / x = x⁺ − x⁻", "—");
    code("// x + y = 4、x ≥ 1（−x ≤ −1）で 2x + 3y を最小に → 最大化 −2x − 3y");
    RUN(const auto q = SimplexLP<double>({{1, 1}, {-1, -1}, {-1, 0}}, {4, -4, -1}, {-2, -3}).solve());
    SHOW(-q.value, "x = 4, y = 0 で 8");
}

}  // namespace simplex_demo

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
    simplex_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    simplex_test::benchmark();    // 速度計測
    simplex_demo::demo();         // 使い方デモ
    return 0;
}
