// =============================================================================
//  【制限付き最速版】巻き戻し配列（直前の変更から順に取り消せる配列。読み書き・取り消しすべて O(1)）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : 戻れるのは「今から履歴を逆にたどった時点」だけ。取り消した変更をやり直す・別の版に飛ぶ・
//                           版を枝分かれさせることはできない（永続配列ではない）
//    できること           : 配列の読み書き、直前の変更の取り消し（undo）、ある時点を覚えて（snapshot）そこまで戻す（rollback）
//    一般版との比較       : 追加テンプレート/ds_persistent_array.cpp の PersistentArray は読み書き O(log N)（どの版にも飛べる）
//                           → この版は読む・書く・undo すべて O(1)、rollback は戻す変更の数だけ
//    制限に合わないとき   : 過去の好きな版を読む・版を分けるなら ds_persistent_array.cpp（スタック・キューだけなら同じファイルの
//                           PersistentStack / PersistentQueue が O(1)）
//
//  しくみ: 書くたびに (位置, 前の値) を履歴に積み、undo で 1 つ取り出して元に戻す
//    例: DFS で子へ進むときに書き換え、戻るときに rollback（木の根からのパスの状態）、巻き戻し Mo（ds_mo_extended.cpp の
//        mo_rollback の snapshot / rollback）、「直前の操作を取り消す」クエリ
//
//  機能一覧（位置は 0〜N−1）                                                                       計算量
//    RollbackArray<T> ra(n, init) / ra(a)   全部 init / 配列 a から                                     O(N)
//      ra.get(i) / ra[i]                    a[i]                                                       O(1)
//      ra.set(i, x)                         a[i] = x（前の値を履歴に積む）                              O(1)
//      ra.add(i, x)                         a[i] += x（set(i, a[i] + x) と同じ）                         O(1)
//      ra.undo()                            最後の set / add を 1 つ取り消す（履歴が空なら false）       O(1)
//      ra.snapshot()                        今の時点の番号（= 履歴の長さ）                               O(1)
//      ra.rollback(t)                       snapshot() が t を返した時点まで戻す                          O(戻す回数)
//      ra.history_size() / ra.size()        履歴の長さ / N                                               O(1)
//      ra.data()                            今の配列全体（const vector<T>&）                             O(1)
//      ra.clear_history()                   今の値はそのままで履歴だけ消す（もう戻れない）                 O(1)
//
//  速さの目安: N = 2×10^5 で set 10^7 回・undo 5×10^6 回 約 0.05 秒（ds_persistent_array.cpp の目安では PersistentArray の set と get 2×10^5 回ずつで 約 0.3 秒）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    rba_test::random_test()  … 時点ごとに配列を丸ごとコピーして持つ愚直解との突き合わせ
//    rba_test::benchmark()    … 速度計測
//    rba_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
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
namespace rba_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace rba_internal

template <class T>
struct RollbackArray {
    explicit RollbackArray(int n = 0, T init = T()) : a_(n < 0 ? 0 : n, init) {
        if (n < 0) rba_internal::fail("RollbackArray", "n は 0 以上");
    }
    explicit RollbackArray(const std::vector<T>& a) : a_(a) {}

    const T& get(int i) const { return a_[i]; }
    const T& operator[](int i) const { return a_[i]; }
    void set(int i, const T& x) {
        hist_.emplace_back(i, a_[i]);
        a_[i] = x;
    }
    void add(int i, const T& x) { set(i, a_[i] + x); }
    bool undo() {  // 最後の変更を 1 つ取り消す
        if (hist_.empty()) return false;
        a_[hist_.back().first] = hist_.back().second;
        hist_.pop_back();
        return true;
    }
    int snapshot() const { return (int)hist_.size(); }
    void rollback(int t) {  // snapshot() == t の時点まで戻す
        if (t < 0 || t > (int)hist_.size()) rba_internal::fail("RollbackArray::rollback", "その時点はもう無い（戻した先より新しい時点には戻れない）");
        while ((int)hist_.size() > t) undo();
    }
    int history_size() const { return (int)hist_.size(); }
    int size() const { return (int)a_.size(); }
    const std::vector<T>& data() const { return a_; }
    void clear_history() { hist_.clear(); }

   private:
    std::vector<T> a_;
    std::vector<std::pair<int, T>> hist_;  // (位置, 書く前の値)
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace rba_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define RBA_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_rollback_array] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 15);
        std::vector<ll> init(n);
        for (auto& v : init) v = rnd(rng, -9, 9);
        RollbackArray<ll> ra = rnd(rng, 0, 1) ? RollbackArray<ll>(init) : RollbackArray<ll>(n, init[0]);
        if (ra.data() != init) init = ra.data();
        std::vector<std::vector<ll>> states = {init};  // states[t] = 時点 t の配列
        for (int t = 0; t < 150; ++t) {
            const int k = (int)rnd(rng, 0, 9);
            const int i = (int)rnd(rng, 0, n - 1);
            const ll x = rnd(rng, -9, 9);
            if (k <= 3) {
                std::vector<ll> s = states.back();
                s[i] = x;
                ra.set(i, x), states.push_back(s);
            } else if (k <= 5) {
                std::vector<ll> s = states.back();
                s[i] += x;
                ra.add(i, x), states.push_back(s);
            } else if (k <= 7) {
                const bool ok = states.size() > 1;
                RBA_CHECK(ra.undo() == ok);
                if (ok) states.pop_back();
            } else if (k == 8) {
                const int to = (int)rnd(rng, 0, (int)states.size() - 1);
                ra.rollback(to), states.resize(to + 1);
            } else if (rnd(rng, 0, 3) == 0) {
                ra.clear_history();
                std::vector<ll> s = states.back();
                states = {s};
            }
            RBA_CHECK(ra.snapshot() == (int)states.size() - 1 && ra.history_size() == ra.snapshot());
            RBA_CHECK(ra.data() == states.back() && ra.size() == n);
            for (int j = 0; j < n; ++j) RBA_CHECK(ra.get(j) == states.back()[j] && ra[j] == ra.get(j));
        }
        ++cases;
    }
    std::cout << "[ds_rollback_array] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_rollback_array] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = 200000, R = TEST_SCALE == 1 ? 1000000 : 10000000;
    std::vector<int> pos(1 << 16);
    for (auto& p : pos) p = (int)rnd(rng, 0, N - 1);
    const auto t0 = std::chrono::steady_clock::now();
    RollbackArray<ll> ra(N, 0);
    ll chk = 0;
    for (int r = 0; r < R; ++r) {
        const int p = pos[r & 0xFFFF];
        ra.set(p, r);
        if (r & 1) ra.undo();
        chk += ra[p];
        if ((r & 0xFFFFF) == 0xFFFFF) ra.rollback(0);
    }
    std::cout << "[ds_rollback_array] benchmark: N = " << N << " で set " << R << " 回・undo " << R / 2 << " 回 " << elapsed_ms(t0) << " ms  (check " << chk % 1000
              << ")\n";
}

}  // namespace rba_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace rba_demo {

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
#define SHOW(expr, comment) rba_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    rba_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】巻き戻し配列（ds_rollback_array） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 戻れるのは履歴を逆にたどった時点だけ（やり直し・別の版へ飛ぶ・枝分かれはできない）\n"
              << "■ できること: 読む・書く・直前の取り消し（undo）が O(1)、覚えた時点まで戻す（rollback）\n"
              << "■ API 一覧\n"
              << "  RollbackArray<T> ra(n, init) / ra(a): get / [i] / set / add / undo / snapshot / rollback / history_size / size / data / clear_history\n";

    api("【1】 RollbackArray<T> ra(n, init)  get / set / add", "const T& / void / void", "a[i] を読む / a[i] = x / a[i] += x（どちらも前の値を履歴に積む）", "O(1)");
    RUN(RollbackArray<int> ra(5, 0));
    RUN(ra.set(1, 10); ra.set(3, 7); ra.add(1, 5));
    SHOW(ra.data(), "");
    SHOW(ra.get(1), "");
    SHOW(ra[3], "");

    api("【2】 undo()", "bool", "最後の set / add を 1 つ取り消す（履歴が空なら false）", "O(1)");
    SHOW(ra.undo(), "add(1, 5) を取り消す");
    SHOW(ra.data(), "");
    SHOW(ra.history_size(), "残りの履歴は 2 つ");

    api("【3】 snapshot() / rollback(t)", "int / void", "今の時点の番号を覚える / その時点まで戻す", "O(1) / O(戻す回数)");
    RUN(int t = ra.snapshot());
    RUN(ra.set(0, 4); ra.set(4, 9); ra.set(0, 8));
    SHOW(ra.data(), "");
    RUN(ra.rollback(t));
    SHOW(ra.data(), "覚えた時点に戻る");
    RUN(ra.clear_history());
    SHOW(ra.undo(), "clear_history で履歴を消すと、もう戻れない");

    api("【4】 DFS で使う", "（使い方の例）", "子へ進む前に snapshot、戻ったら rollback（根から今の頂点までの状態を持てる）", "全体 O(N)");
    {
        // 木の各頂点で「根からその頂点までに出てきた色の種類数」を数える
        const std::vector<std::vector<int>> ch = {{1, 2}, {3, 4}, {5}, {}, {}, {}};
        const std::vector<int> color = {0, 1, 0, 2, 1, 1};
        RollbackArray<int> cnt(3, 0);
        int kinds = 0;
        std::vector<int> ans(6);
        auto dfs = [&](auto&& self, int v) -> void {
            const int snap = cnt.snapshot(), k = kinds;
            if (cnt[color[v]] == 0) ++kinds;
            cnt.add(color[v], 1);
            ans[v] = kinds;
            for (int c : ch[v]) self(self, c);
            cnt.rollback(snap), kinds = k;
        };
        dfs(dfs, 0);
        note("各頂点で根からの色の種類数 → " + to_s(ans));
    }
}

}  // namespace rba_demo

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
    rba_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    rba_test::benchmark();    // 速度計測
    rba_demo::demo();         // 使い方デモ
    return 0;
}
