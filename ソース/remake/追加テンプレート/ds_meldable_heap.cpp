// =============================================================================
//  マージ可能なヒープ（左偏木：2 つのヒープの併合 O(log N)・全体への加算（遅延））   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  priority_queue は 2 つをまとめるのに小さい方を全部入れ直す必要がある（マージテクで O(N log^2 N)）。
//  左偏木（leftist heap）は 2 つのヒープを O(log N) で併合でき、ヒープ全体に x を足す（遅延）こともできる
//    例: Union-Find と組み合わせて「成分ごとの最小の辺」を持つ（有向最小全域木 Chu-Liu/Edmonds など）、
//        木の上で子のヒープを親に併合しながら上がる（スロープトリックの木版など）
//
//  機能一覧（ヒープは int の番号（根）で表す。−1 が空。T は値の型、IS_MIN = false なら最大ヒープ）     計算量
//    MeldableHeap<T, IS_MIN> H          頂点をためておく入れ物（いくつものヒープを持てる）
//    H.push(h, x)                       ヒープ h に x を入れた新しい根                                      O(log N)
//    H.meld(a, b)                       ヒープ a と b を併合した根（a, b はもう使わない）                     O(log N)
//    H.top(h) / H.pop(h)                最小（最大）の値 / それを取り除いた新しい根                         O(1) / O(log N)
//    H.add_all(h, x)                    ヒープ h の全部の値に x を足す（遅延）                              O(1)
//    H.size(h) / H.empty(h)             要素の数 / 空か（h == −1）                                          O(1)
//    H.to_vector(h)                     中身（順不同）                                                     O(サイズ)
//
//  注意: 根の番号は push・meld・pop のたびに変わるので、返り値で上書きする（h = H.push(h, x)）
//        空のヒープに top・pop すると止める
//  速さの目安: push と pop を 10^6 回ずつ 約 0.9 秒（ふつうのヒープなら priority_queue の方が速い）、10^5 個のヒープを全部併合 約 0.005 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    meld_test::random_test()  … multiset を使う愚直解との突き合わせ
//    meld_test::benchmark()    … 速度計測
//    meld_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace meld_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace meld_internal

template <class T, bool IS_MIN = true>
struct MeldableHeap {
    int push(int h, const T& x) {
        nd_.push_back({x, T(0), -1, -1, 1, 1});
        return meld(h, (int)nd_.size() - 1);
    }
    int meld(int a, int b) {
        // 右の背骨に沿って併合（深さは O(log N) なので再帰でもよいが、ここでは積んでから戻る）
        std::vector<int>& st = stack_;
        st.clear();
        while (a != -1 && b != -1) {
            push_down(a), push_down(b);
            if (better(nd_[b].val, nd_[a].val)) std::swap(a, b);
            st.push_back(a);
            a = nd_[a].r;
        }
        int cur = a == -1 ? b : a;
        for (size_t k = st.size(); k-- > 0;) {
            const int v = st[k];
            nd_[v].r = cur;
            int l = nd_[v].l;
            if (rank(l) < rank(cur)) std::swap(nd_[v].l, nd_[v].r);
            nd_[v].rk = rank(nd_[v].r) + 1;
            nd_[v].sz = 1 + size(nd_[v].l) + size(nd_[v].r);
            cur = v;
        }
        return cur;
    }
    const T& top(int h) const {
        if (h == -1) meld_internal::fail("MeldableHeap::top", "空");
        return nd_[h].val;
    }
    int pop(int h) {
        if (h == -1) meld_internal::fail("MeldableHeap::pop", "空");
        push_down(h);
        return meld(nd_[h].l, nd_[h].r);
    }
    void add_all(int h, const T& x) {
        if (h == -1) return;
        nd_[h].val += x, nd_[h].lazy += x;
    }
    int size(int h) const { return h == -1 ? 0 : nd_[h].sz; }
    bool empty(int h) const { return h == -1; }
    std::vector<T> to_vector(int h) {
        std::vector<T> out;
        std::vector<int> st;
        if (h != -1) st.push_back(h);
        while (!st.empty()) {
            const int v = st.back();
            st.pop_back();
            push_down(v);
            out.push_back(nd_[v].val);
            if (nd_[v].l != -1) st.push_back(nd_[v].l);
            if (nd_[v].r != -1) st.push_back(nd_[v].r);
        }
        return out;
    }

   private:
    struct Node {
        T val, lazy;
        int l, r, rk, sz;
    };
    std::vector<Node> nd_;
    std::vector<int> stack_;
    static bool better(const T& a, const T& b) { return IS_MIN ? a < b : b < a; }
    int rank(int v) const { return v == -1 ? 0 : nd_[v].rk; }
    void push_down(int v) {
        if (nd_[v].lazy == T(0)) return;
        for (int c : {nd_[v].l, nd_[v].r})
            if (c != -1) nd_[c].val += nd_[v].lazy, nd_[c].lazy += nd_[v].lazy;
        nd_[v].lazy = T(0);
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace meld_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define MH_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_meldable_heap] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

template <bool IS_MIN>
static void test_kind(std::mt19937_64& rng, long long& cases) {
    MeldableHeap<ll, IS_MIN> H;
    const int K = (int)rnd(rng, 1, 6);
    std::vector<int> h(K, -1);
    std::vector<std::multiset<ll>> b(K);
    for (int t = 0; t < 300; ++t) {
        const int i = (int)rnd(rng, 0, K - 1), k = (int)rnd(rng, 0, 9);
        if (k < 4) {
            const ll x = rnd(rng, -20, 20);
            h[i] = H.push(h[i], x), b[i].insert(x);
        } else if (k < 6 && !b[i].empty()) {
            const ll top = IS_MIN ? *b[i].begin() : *b[i].rbegin();
            MH_CHECK(H.top(h[i]) == top);
            h[i] = H.pop(h[i]);
            b[i].erase(b[i].find(top));
        } else if (k == 6) {
            const int j = (int)rnd(rng, 0, K - 1);
            if (i != j) {
                h[i] = H.meld(h[i], h[j]), h[j] = -1;
                b[i].insert(b[j].begin(), b[j].end()), b[j].clear();
            }
        } else if (k == 7) {
            const ll x = rnd(rng, -5, 5);
            H.add_all(h[i], x);
            std::multiset<ll> nb;
            for (ll v : b[i]) nb.insert(v + x);
            b[i] = nb;
        }
        MH_CHECK(H.size(h[i]) == (int)b[i].size() && H.empty(h[i]) == b[i].empty());
        if (!b[i].empty()) MH_CHECK(H.top(h[i]) == (IS_MIN ? *b[i].begin() : *b[i].rbegin()));
        if (rnd(rng, 0, 30) == 0) {
            auto v = H.to_vector(h[i]);
            std::sort(v.begin(), v.end());
            MH_CHECK(v == std::vector<ll>(b[i].begin(), b[i].end()));
        }
        ++cases;
    }
}

static void random_test() {
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 300 * TEST_SCALE; ++it) test_kind<true>(rng, cases), test_kind<false>(rng, cases);
    std::cout << "[ds_meldable_heap] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_meldable_heap] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 100000 : 1000000, M = TEST_SCALE == 1 ? 10000 : 100000;
    auto t0 = std::chrono::steady_clock::now();
    MeldableHeap<ll> H;
    int h = -1;
    ll chk = 0;
    for (int i = 0; i < N; ++i) h = H.push(h, rnd(rng, 0, 1000000000));
    for (int i = 0; i < N; ++i) chk += H.top(h), h = H.pop(h);
    const ll t1 = elapsed_ms(t0);
    t0 = std::chrono::steady_clock::now();
    MeldableHeap<ll> G;
    std::vector<int> hs(M, -1);
    for (int i = 0; i < M; ++i) hs[i] = G.push(-1, rnd(rng, 0, 1000000000));
    for (int step = 1; step < M; step *= 2)
        for (int i = 0; i + step < M; i += 2 * step) hs[i] = G.meld(hs[i], hs[i + step]);
    chk += G.top(hs[0]);
    const ll t2 = elapsed_ms(t0);
    std::cout << "[ds_meldable_heap] benchmark: push と pop を " << N << " 回ずつ " << t1 << " ms / " << M << " 個のヒープを全部併合 " << t2 << " ms  (check "
              << chk % 1000 << ")\n";
}

}  // namespace meld_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace meld_demo {

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
#define SHOW(expr, comment) meld_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    meld_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("マージ可能なヒープ（ds_meldable_heap） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  MeldableHeap<T, IS_MIN> H: push(h, x) / meld(a, b) / pop(h) → 新しい根 / top(h) / add_all(h, x) / size(h) / empty(h) / to_vector(h)\n";

    api("【1】 push(h, x) / top(h) / pop(h)", "int / const T& / int", "ヒープは根の番号（−1 が空）。返り値で上書きする", "O(log N) / O(1) / O(log N)");
    RUN(MeldableHeap<long long> H);
    RUN(int a = -1, b = -1);
    RUN(a = H.push(a, 5); a = H.push(a, 2); a = H.push(a, 8));
    RUN(b = H.push(b, 3); b = H.push(b, 7));
    SHOW(H.top(a), "");
    SHOW(H.top(b), "");

    api("【2】 meld(a, b) / add_all(h, x)", "int / void", "2 つのヒープを併合 / 全部に x を足す（遅延）", "O(log N) / O(1)");
    RUN(H.add_all(b, -2));  // b = {1, 5}
    RUN(int c = H.meld(a, b));
    SHOW(H.size(c), "");
    SHOW(H.top(c), "b の 3 − 2 = 1");
    RUN(c = H.pop(c));
    SHOW(H.top(c), "");

    api("【3】 MeldableHeap<T, false>：最大ヒープ", "—", "IS_MIN = false で大きい順", "—");
    RUN(MeldableHeap<int, false> M);
    RUN(int m = M.push(M.push(M.push(-1, 4), 9), 1));
    SHOW(M.top(m), "");
}

}  // namespace meld_demo

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
    meld_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    meld_test::benchmark();    // 速度計測
    meld_demo::demo();         // 使い方デモ
    return 0;
}
