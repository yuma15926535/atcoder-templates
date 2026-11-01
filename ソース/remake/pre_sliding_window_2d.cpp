// =============================================================================
//  【制限付き最速版】2 次元の固定サイズの窓の最小・最大（h × w の窓をずらした全部の位置の答えを O(HW)）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : 窓の大きさ h × w が固定。好きな長方形の最小を聞くことはできない。値の更新もない
//    できること           : 全部の位置 (i, j) について、[i, i + h) × [j, j + w) の最小・最大（または好きな「良い方」の比較）を一度に出す。
//                           1 次元（長さ w の窓）の版もある
//    一般版との比較       : 2 次元の区間最小を一般に答えるなら 2 次元 Sparse Table（前処理 O(HW log H log W)・メモリも同じ）か
//                           2 次元セグ木（1 回 O(log H log W)）が要る（今のテンプレートには無い。愚直だと O(HW hw)）
//                           → この版は全部の窓で O(HW)・メモリ O(HW)
//    制限に合わないとき   : 窓の大きさが毎回違うなら 2 次元 Sparse Table を書く。1 次元で窓が左右に伸び縮みするなら SlidingWindowDequeTemplateTest.txt
//
//  しくみ: 行ごとに長さ w の窓の最小を出し（O(HW)）、その結果を列ごとに長さ h の窓で同じようにする（O(HW)）。
//          1 次元の窓は「長さ w のブロックごとの前からの累積と後ろからの累積」の 2 つを比べるだけ（van Herk / Gil-Werman 法。分岐が少なく速い）
//    例: 「h × w の範囲の最大値の最小」「すべての h × w の正方形で max − min ≤ K か」
//
//  機能一覧（a は vector<vector<T>>。H = a の行数、W = 列数。結果は (H − h + 1) × (W − w + 1)）           計算量
//    sliding_window_min_2d(a, h, w)       res[i][j] = [i, i + h) × [j, j + w) の最小                              O(HW)
//    sliding_window_max_2d(a, h, w)       同じく最大                                                             O(HW)
//    sliding_window_2d(a, h, w, better)   better(x, y) = 「x が y より良い」（最小なら x < y）で一番良い値（同じ良さならどれか）  O(HW)
//    sliding_window_min_1d(v, w) / sliding_window_max_1d(v, w) / sliding_window_1d(v, w, better)
//                                         res[i] = [i, i + w) の最小 / 最大 / 一番良い値（長さ N − w + 1）            O(N)
//  h > H や w > W なら空を返す。h・w が 0 以下、行の長さがそろっていないときは止める
//  速さの目安: H = W = 3000・窓 50 × 70 の min で 約 0.07 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    sw2_test::random_test()  … 窓の中を全部見る愚直解との突き合わせ
//    sw2_test::benchmark()    … 速度計測
//    sw2_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace sw2_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
// src[0..n) の長さ w の窓を dst[0..n − w] に書く（van Herk / Gil-Werman 法：長さ w のブロックごとの前からと後ろからの累積。
// 分岐がほとんどないので、単調な両端キューより速い）。pre・suf は作業用（n 個以上）
template <class T, class Better>
inline void slide(const T* src, int n, int w, T* dst, Better& better, T* pre, T* suf) {
    for (int b = 0; b < n; b += w) {
        const int e = std::min(n, b + w);
        pre[b] = src[b];
        for (int i = b + 1; i < e; ++i) pre[i] = better(src[i], pre[i - 1]) ? src[i] : pre[i - 1];
        suf[e - 1] = src[e - 1];
        for (int i = e - 2; i >= b; --i) suf[i] = better(src[i], suf[i + 1]) ? src[i] : suf[i + 1];
    }
    for (int i = 0; i + w <= n; ++i) dst[i] = better(pre[i + w - 1], suf[i]) ? pre[i + w - 1] : suf[i];
}
}  // namespace sw2_internal

template <class T, class Better>
std::vector<T> sliding_window_1d(const std::vector<T>& v, int w, Better better) {
    if (w <= 0) sw2_internal::fail("sliding_window_1d", "窓の長さ w は 1 以上");
    const int n = (int)v.size();
    if (w > n) return {};
    std::vector<T> res(n - w + 1), pre(n), suf(n);
    sw2_internal::slide(v.data(), n, w, res.data(), better, pre.data(), suf.data());
    return res;
}
template <class T>
std::vector<T> sliding_window_min_1d(const std::vector<T>& v, int w) { return sliding_window_1d(v, w, std::less<T>()); }
template <class T>
std::vector<T> sliding_window_max_1d(const std::vector<T>& v, int w) { return sliding_window_1d(v, w, std::greater<T>()); }

template <class T, class Better>
std::vector<std::vector<T>> sliding_window_2d(const std::vector<std::vector<T>>& a, int h, int w, Better better) {
    if (h <= 0 || w <= 0) sw2_internal::fail("sliding_window_2d", "窓の大きさ h, w は 1 以上");
    const int H = (int)a.size();
    const int W = H ? (int)a[0].size() : 0;
    for (int i = 0; i < H; ++i)
        if ((int)a[i].size() != W) sw2_internal::fail("sliding_window_2d", "行の長さがそろっていない");
    if (h > H || w > W) return {};
    const int W2 = W - w + 1, H2 = H - h + 1;
    std::vector<T> row((size_t)H * W2);  // 行ごとの結果（H × W2、1 本の配列）
    std::vector<T> pre(std::max(H, W)), suf(std::max(H, W));
    for (int i = 0; i < H; ++i) sw2_internal::slide(a[i].data(), W, w, row.data() + (long long)i * W2, better, pre.data(), suf.data());
    std::vector<std::vector<T>> res(H2, std::vector<T>(W2));
    // 列方向：BJ 列ずつ作業用の配列に縦に並べ直してから窓を取る（飛び飛びに読むよりキャッシュに乗って速い）
    const int BJ = 16;
    std::vector<T> col((long long)BJ * H), out((long long)BJ * H2);
    for (int j0 = 0; j0 < W2; j0 += BJ) {
        const int nb = std::min(BJ, W2 - j0);
        for (int i = 0; i < H; ++i) {
            const T* src = row.data() + (long long)i * W2 + j0;
            for (int jj = 0; jj < nb; ++jj) col[(long long)jj * H + i] = src[jj];
        }
        for (int jj = 0; jj < nb; ++jj) sw2_internal::slide(col.data() + (long long)jj * H, H, h, out.data() + (long long)jj * H2, better, pre.data(), suf.data());
        for (int i = 0; i < H2; ++i) {
            T* dst = res[i].data() + j0;
            for (int jj = 0; jj < nb; ++jj) dst[jj] = out[(long long)jj * H2 + i];
        }
    }
    return res;
}
template <class T>
std::vector<std::vector<T>> sliding_window_min_2d(const std::vector<std::vector<T>>& a, int h, int w) { return sliding_window_2d(a, h, w, std::less<T>()); }
template <class T>
std::vector<std::vector<T>> sliding_window_max_2d(const std::vector<std::vector<T>>& a, int h, int w) {
    return sliding_window_2d(a, h, w, std::greater<T>());
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace sw2_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define SW2_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[pre_sliding_window_2d] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    std::mt19937_64 rng(20261101);
    long long cases = 0;
    for (int it = 0; it < 600 * TEST_SCALE; ++it) {
        const int H = (int)rnd(rng, 1, 9), W = (int)rnd(rng, 1, 9);
        std::vector<std::vector<int>> a(H, std::vector<int>(W));
        for (auto& r : a)
            for (auto& x : r) x = (int)rnd(rng, -5, 5);
        const int h = (int)rnd(rng, 1, H + 1), w = (int)rnd(rng, 1, W + 1);
        auto mn = sliding_window_min_2d(a, h, w), mx = sliding_window_max_2d(a, h, w);
        auto ab = sliding_window_2d(a, h, w, [](int x, int y) { return std::abs(x) < std::abs(y); });  // 絶対値が小さい方（同じ良さならどれか）
        if (h > H || w > W) {
            SW2_CHECK(mn.empty() && mx.empty() && ab.empty());
        } else {
            SW2_CHECK((int)mn.size() == H - h + 1 && (int)mn[0].size() == W - w + 1);
            for (int i = 0; i + h <= H; ++i)
                for (int j = 0; j + w <= W; ++j) {
                    int lo = 1 << 30, hi = -(1 << 30), ba = 1 << 30;
                    for (int x = i; x < i + h; ++x)
                        for (int y = j; y < j + w; ++y) {
                            lo = std::min(lo, a[x][y]), hi = std::max(hi, a[x][y]);
                            ba = std::min(ba, std::abs(a[x][y]));
                        }
                    SW2_CHECK(mn[i][j] == lo && mx[i][j] == hi && std::abs(ab[i][j]) == ba);
                }
        }
        // 1 次元
        const int n = (int)rnd(rng, 0, 12), w1 = (int)rnd(rng, 1, 13);
        std::vector<ll> v(n);
        for (auto& x : v) x = rnd(rng, -9, 9);
        auto m1 = sliding_window_min_1d(v, w1), M1 = sliding_window_max_1d(v, w1);
        if (w1 > n) {
            SW2_CHECK(m1.empty() && M1.empty());
        } else {
            for (int i = 0; i + w1 <= n; ++i) {
                SW2_CHECK(m1[i] == *std::min_element(v.begin() + i, v.begin() + i + w1));
                SW2_CHECK(M1[i] == *std::max_element(v.begin() + i, v.begin() + i + w1));
            }
        }
        ++cases;
    }
    std::cout << "[pre_sliding_window_2d] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[pre_sliding_window_2d] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int H = TEST_SCALE == 1 ? 1000 : 3000, W = H;
    std::vector<std::vector<int>> a(H, std::vector<int>(W));
    for (auto& r : a)
        for (auto& x : r) x = (int)rnd(rng, 0, 1000000000);
    const auto t0 = std::chrono::steady_clock::now();
    auto res = sliding_window_min_2d(a, 50, 70);
    ll chk = 0;
    for (auto& r : res) chk += r[r.size() / 2];
    std::cout << "[pre_sliding_window_2d] benchmark: H = W = " << H << "・窓 50 × 70 の min " << elapsed_ms(t0) << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace sw2_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace sw2_demo {

using ll = long long;
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
#define SHOW(expr, comment) sw2_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    sw2_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】2 次元の固定サイズの窓の最小・最大（pre_sliding_window_2d） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 窓の大きさ h × w が固定（好きな長方形は聞けない）・更新なし\n"
              << "■ できること: 全部の位置の h × w の窓の最小・最大（好きな比較でも）を O(HW) で一度に。1 次元版も\n"
              << "■ API 一覧\n"
              << "  sliding_window_min_2d / sliding_window_max_2d / sliding_window_2d(a, h, w, better)\n"
              << "  sliding_window_min_1d / sliding_window_max_1d / sliding_window_1d(v, w, better)\n";

    RUN(std::vector<std::vector<int>> a = {{3, 1, 4, 1}, {5, 9, 2, 6}, {5, 3, 5, 8}});
    api("【1】 sliding_window_min_2d(a, h, w) / sliding_window_max_2d(a, h, w)", "vector<vector<T>>（(H − h + 1) × (W − w + 1)）",
        "res[i][j] = [i, i + h) × [j, j + w) の最小 / 最大", "O(HW)");
    SHOW(sliding_window_min_2d(a, 2, 2), "2 × 2 の窓の最小（3 行 4 列 → 2 行 3 列）");
    SHOW(sliding_window_max_2d(a, 2, 3), "2 × 3 の窓の最大");
    SHOW(sliding_window_min_2d(a, 4, 1), "窓が大きすぎる → 空");

    api("【2】 sliding_window_2d(a, h, w, better)", "vector<vector<T>>", "better(x, y) = x が y より良い（最小なら x < y）の意味で一番良い値", "O(HW)");
    RUN(auto far5 = [](int x, int y) { return std::abs(x - 5) > std::abs(y - 5); });
    SHOW(sliding_window_2d(a, 2, 2, far5), "5 から一番遠い値（同じ良さが複数あるときはそのどれか）");

    api("【3】 sliding_window_min_1d(v, w) / sliding_window_max_1d(v, w) / sliding_window_1d(v, w, better)", "vector<T>（長さ N − w + 1）",
        "res[i] = [i, i + w) の最小 / 最大 / 一番良い値", "O(N)");
    RUN(std::vector<ll> v = {4, 2, 12, 3, 8, 5, 1});
    SHOW(sliding_window_min_1d(v, 3), "長さ 3 の窓の最小");
    SHOW(sliding_window_max_1d(v, 3), "長さ 3 の窓の最大");
    {
        code("// 使い方の例：すべての 2 × 2 の窓で max − min ≤ 4 か");
        auto mn = sliding_window_min_2d(a, 2, 2), mx = sliding_window_max_2d(a, 2, 2);
        bool ok = true;
        for (size_t i = 0; i < mn.size(); ++i)
            for (size_t j = 0; j < mn[0].size(); ++j) ok &= mx[i][j] - mn[i][j] <= 4;
        note("答え → " + to_s(ok) + "（左上の窓 {3,1,5,9} は 9 − 1 = 8）");
    }
}

}  // namespace sw2_demo

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
    sw2_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    sw2_test::benchmark();    // 速度計測
    sw2_demo::demo();         // 使い方デモ
    return 0;
}
