// =============================================================================
//  【制限付き最速版】追加だけ / 削除だけの mex（入っていない最小の非負整数を、追加・削除ともならし O(1)）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : 追加と削除をまぜられない（AddOnlyMex は追加だけ、EraseOnlyMex は削除だけ）。
//                           mex の上限 cap を最初に決める（cap 以上の値は mex に関係しないので数えない）。k 番目の欠番などはない
//    できること           : 多重集合への追加（または削除）と、mex・個数・含むかの取得
//    一般版との比較       : MexCounterTemplateTest.txt は追加・削除 O(log N)（mex O(1)）、MexBITTemplate.txt は O(log MAX)
//                           → この版は追加・削除・mex すべて ならし O(1)（全体で O(N + cap)）
//    制限に合わないとき   : 追加と削除がまざるなら MexCounterTemplateTest.txt / 330E.txt、欠番の k 番目なども要るなら MexBITTemplate.txt
//
//  しくみ: 追加だけなら mex は増える一方なので、mex の位置を右へ進めるだけ。削除だけなら mex は減る一方で、
//          個数が 0 になった値 x と比べて小さい方にするだけ
//    例: 「数を 1 つずつ追加して毎回 mex を出力」、オフラインで追加だけの問題を逆から見て「削除だけ」にする
//
//  機能一覧（値は long long で渡してよい。負の値・cap 以上の値は個数を数えず size だけ増減する）           計算量
//    AddOnlyMex mx(cap)                   空の集合。cap = mex の上限（入れる回数の上限 N でよい：mex ≤ 入れた個数）  O(cap)
//      mx.insert(x)                       x を 1 個入れる                                                      ならし O(1)
//      mx.mex()                           入っていない最小の非負整数（cap 以上なら cap）                       O(1)
//    EraseOnlyMex mx(a, cap = a.size())   a の要素を全部入れた集合（cap = mex の上限。省略なら a の長さ）       O(N + cap)
//      mx.erase(x)                        x を 1 個消す（入っていなければ何もせず false）                       O(1)
//      mx.mex()                           入っていない最小の非負整数（cap 以上なら cap）                       O(1)
//    共通: count(x)（0 ≤ x < cap の個数。それ以外は 0）/ contains(x) / size()（入っている個数）/ capacity()
//    EraseOnlyMex の erase で負・cap 以上の値を渡すと、そういう値がまだ残っていれば 1 個減らして true（どの値かは区別しない）
//
//  速さの目安: AddOnlyMex に 10^7 回 insert + 毎回 mex 約 0.11 秒。10^6 回で比べると 一般版の形（欠番を set で持つ MexCounter 型）約 0.55 秒、
//              この版 約 0.07 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    mmx_test::random_test()  … 毎回 0 から数える愚直解との突き合わせ
//    mmx_test::benchmark()    … 速度計測（一般版との比較つき）
//    mmx_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace mmx_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace mmx_internal

// ---- 追加だけの mex ----
struct AddOnlyMex {
    explicit AddOnlyMex(int cap = 0) : cap_(cap), cnt_(cap < 0 ? 0 : cap, 0) {
        if (cap < 0) mmx_internal::fail("AddOnlyMex", "cap は 0 以上");
    }
    void insert(long long x) {
        ++size_;
        if (x < 0 || x >= cap_) return;
        if (cnt_[x]++ == 0 && x == mex_)
            while (mex_ < cap_ && cnt_[mex_] > 0) ++mex_;
    }
    int mex() const { return mex_; }
    int count(long long x) const { return (0 <= x && x < cap_) ? cnt_[x] : 0; }
    bool contains(long long x) const { return count(x) > 0; }
    long long size() const { return size_; }
    int capacity() const { return cap_; }

   private:
    int cap_, mex_ = 0;
    long long size_ = 0;
    std::vector<int> cnt_;
};

// ---- 削除だけの mex ----
struct EraseOnlyMex {
    template <class T>
    explicit EraseOnlyMex(const std::vector<T>& a, int cap = -1) : cap_(cap < 0 ? (int)a.size() : cap), size_((long long)a.size()) {
        cnt_.assign(cap_, 0);
        for (const T& v : a) {
            const long long x = (long long)v;
            if (0 <= x && x < cap_) {
                ++cnt_[x];
            } else {
                ++others_;
            }
        }
        mex_ = 0;
        while (mex_ < cap_ && cnt_[mex_] > 0) ++mex_;
    }
    bool erase(long long x) {  // 1 個消す（入っていなければ false。負・cap 以上の値はまとめて数えているので、そういう値が残っていれば true）
        if (x < 0 || x >= cap_) {
            if (others_ == 0) return false;
            --others_, --size_;
            return true;
        }
        if (cnt_[x] == 0) return false;
        --size_;
        if (--cnt_[x] == 0 && x < mex_) mex_ = (int)x;
        return true;
    }
    int mex() const { return mex_; }
    int count(long long x) const { return (0 <= x && x < cap_) ? cnt_[x] : 0; }
    bool contains(long long x) const { return count(x) > 0; }
    long long size() const { return size_; }
    int capacity() const { return cap_; }

   private:
    int cap_, mex_ = 0;
    long long size_, others_ = 0;  // others_ = 負・cap 以上の値の個数（どの値かは区別しない）
    std::vector<int> cnt_;
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace mmx_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define MMX_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_monotone_mex] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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

static int brute_mex(const std::vector<int>& cnt, int cap) {
    int m = 0;
    while (m < cap && m < (int)cnt.size() && cnt[m] > 0) ++m;
    return std::min(m, cap);
}

static void random_test() {
    std::mt19937_64 rng(20261101);
    long long cases = 0;
    for (int it = 0; it < 1000 * TEST_SCALE; ++it) {
        const int cap = (int)rnd(rng, 0, 20);
        // 追加だけ
        {
            AddOnlyMex mx(cap);
            std::vector<int> cnt(64, 0);
            long long sz = 0;
            for (int t = 0; t < 40; ++t) {
                const ll x = rnd(rng, -3, 25);
                mx.insert(x), ++sz;
                if (0 <= x) ++cnt[x];
                MMX_CHECK(mx.mex() == brute_mex(cnt, cap));
                const ll y = rnd(rng, -3, 25);
                MMX_CHECK(mx.count(y) == ((0 <= y && y < cap) ? cnt[y] : 0) && mx.contains(y) == (mx.count(y) > 0));
                MMX_CHECK(mx.size() == sz && mx.capacity() == cap);
            }
        }
        // 削除だけ
        {
            const int n = (int)rnd(rng, 0, 30);
            std::vector<ll> a(n);
            std::vector<int> cnt(64, 0);
            for (auto& v : a) v = rnd(rng, -2, 25), cnt[std::max(v, 0LL)] += (v >= 0);
            const int c2 = rnd(rng, 0, 1) ? n : cap;
            EraseOnlyMex mx(a, rnd(rng, 0, 1) ? -1 : c2);
            const int cp = mx.capacity();
            MMX_CHECK(cp == (cp == n ? n : c2));
            long long sz = n, others = 0;
            for (ll v : a) others += (v < 0 || v >= cp);
            MMX_CHECK(mx.mex() == brute_mex(cnt, cp));
            for (int t = 0; t < 40; ++t) {
                const ll x = rnd(rng, -3, 25);
                const bool out = (x < 0 || x >= cp);
                const bool has = out ? others > 0 : cnt[x] > 0;
                MMX_CHECK(mx.erase(x) == has);
                if (has) {
                    --sz;
                    if (out) {
                        --others;
                    } else {
                        --cnt[x];
                    }
                }
                MMX_CHECK(mx.mex() == brute_mex(cnt, cp) && mx.size() == sz);
                const ll y = rnd(rng, -3, 25);
                MMX_CHECK(mx.count(y) == ((0 <= y && y < cp) ? cnt[y] : 0) && mx.contains(y) == (mx.count(y) > 0));
            }
        }
        ++cases;
    }
    std::cout << "[ds_monotone_mex] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_monotone_mex] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    {
        const int N = TEST_SCALE == 1 ? 1000000 : 10000000;
        std::vector<int> xs(N);
        for (auto& x : xs) x = (int)rnd(rng, 0, N - 1);
        const auto t0 = std::chrono::steady_clock::now();
        AddOnlyMex mx(N);
        ll chk = 0;
        for (int x : xs) mx.insert(x), chk += mx.mex();
        std::cout << "[ds_monotone_mex] benchmark: AddOnlyMex insert + mex を " << N << " 回 " << elapsed_ms(t0) << " ms  (check " << chk % 1000 << ")\n";
    }
    {
        // 一般版の形（欠番を std::set で持つ。MexCounterTemplateTest.txt と同じやり方）
        const int N = TEST_SCALE == 1 ? 100000 : 1000000;
        std::vector<int> xs(N);
        for (auto& x : xs) x = (int)rnd(rng, 0, N - 1);
        auto t0 = std::chrono::steady_clock::now();
        std::vector<int> cnt(N + 1, 0);
        std::set<int> missing;
        for (int i = 0; i <= N; ++i) missing.insert(missing.end(), i);
        ll chk = 0;
        for (int x : xs) {
            if (cnt[x]++ == 0) missing.erase(x);
            chk += *missing.begin();
        }
        const ll ms_set = elapsed_ms(t0);
        t0 = std::chrono::steady_clock::now();
        AddOnlyMex mx(N);
        ll chk2 = 0;
        for (int x : xs) mx.insert(x), chk2 += mx.mex();
        MMX_CHECK(chk == chk2);
        std::cout << "[ds_monotone_mex] benchmark: 同じ " << N << " 回で 一般版（set で欠番）" << ms_set << " ms / AddOnlyMex " << elapsed_ms(t0) << " ms\n";
    }
}

}  // namespace mmx_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace mmx_demo {

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
#define SHOW(expr, comment) mmx_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    mmx_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】追加だけ / 削除だけの mex（ds_monotone_mex） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 追加と削除をまぜられない・mex の上限 cap を最初に決める\n"
              << "■ できること: 追加（または削除）・mex・個数・含むか が すべて（ならし）O(1)\n"
              << "■ API 一覧\n"
              << "  AddOnlyMex mx(cap): insert / mex / count / contains / size / capacity\n"
              << "  EraseOnlyMex mx(a, cap = a.size()): erase / mex / count / contains / size / capacity\n";

    api("【1】 AddOnlyMex mx(cap)  insert(x) / mex()", "void / int", "x を 1 個入れる / 入っていない最小の非負整数（cap 以上なら cap）", "ならし O(1) / O(1)");
    RUN(AddOnlyMex mx(10));
    RUN(mx.insert(1); mx.insert(0); mx.insert(3));
    SHOW(mx.mex(), "{0, 1, 3} → 2");
    RUN(mx.insert(2));
    SHOW(mx.mex(), "2 が入って 0〜3 がそろった → 4");
    RUN(mx.insert(-5); mx.insert(100));
    SHOW(mx.mex(), "負の値・cap 以上の値は mex に関係しない");
    SHOW(mx.count(3), "");
    SHOW(mx.contains(100), "cap 以上の値は数えていない（false）");
    SHOW(mx.size(), "入れた個数（数えない値も含む）");

    api("【2】 EraseOnlyMex mx(a, cap)  erase(x) / mex()", "bool / int", "x を 1 個消す（無ければ false）/ 入っていない最小の非負整数", "O(1) / O(1)");
    RUN(EraseOnlyMex em(std::vector<int>{0, 1, 1, 2, 3, 5}));
    SHOW(em.mex(), "{0, 1, 1, 2, 3, 5} → 4（cap は a の長さ 6）");
    SHOW(em.erase(1), "1 を 1 個消す（まだ 1 個ある）");
    SHOW(em.mex(), "");
    SHOW(em.erase(1), "");
    SHOW(em.mex(), "1 が無くなった → 1");
    SHOW(em.erase(7), "入っていない → false");
    SHOW(em.capacity(), "");

    api("【3】 追加だけの問題を逆から見る", "（使い方の例）", "オフラインで「全部入れてから 1 つずつ消す」に直せば EraseOnlyMex でも解ける", "全体 O(N + cap)");
    {
        std::vector<int> add = {2, 0, 1, 4, 3};
        EraseOnlyMex e(add, 5);
        std::vector<int> ans(add.size());
        for (int i = (int)add.size() - 1; i >= 0; --i) ans[i] = e.mex(), e.erase(add[i]);
        note("add = [2, 0, 1, 4, 3] を前から入れたときの各時点の mex → " + to_s(ans));
    }
}

}  // namespace mmx_demo

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
    mmx_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    mmx_test::benchmark();    // 速度計測
    mmx_demo::demo();         // 使い方デモ
    return 0;
}
