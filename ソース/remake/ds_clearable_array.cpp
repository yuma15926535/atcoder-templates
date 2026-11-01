// =============================================================================
//  【制限付き最速版】O(1) で全部消せる配列・集合（世代番号つき配列・訪問済みの印・sparse set）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  ★ 制限付き最速版（機能を絞る代わりに、一般版よりオーダーが速い）
//    制限（できないこと） : キー（添字）は 0〜n−1 の整数だけ。順序（最小・lower_bound・小さい順の列挙）はない。範囲外は確かめない
//    できること           : 配列の読み書き・訪問済みの印・整数の集合（追加・削除・含むか・中身の列挙）を、全消去 O(1) 付きで
//    一般版との比較       : vector を毎回 assign / fill し直すと全消去 O(n)、std::set は追加・削除 O(log N)
//                           → この版は全消去 O(1)、集合の操作 O(1)、中身の列挙 O(要素数)（bool 配列だと O(n)）
//                           （HashSetTemplate.txt も平均 O(1) でオーダーは同じ。こちらはハッシュがない分の定数倍が速い）
//    制限に合わないとき   : キーが大きい・負の値なら HashSetTemplate.txt / HashMapTemplate.txt、順序が要るなら std::set や OrderedBITTemplate.txt
//
//  しくみ: 各要素に「書いたときの世代」を持たせ、clear は世代を 1 つ進めるだけ（古い世代のものは初期値とみなす）。
//          SparseSet は「中身の配列 dense」と「各値が dense の何番目か pos」を持つ（clear は個数を 0 にするだけ）
//    例: スタート地点を変えて狭い範囲の BFS・DFS を何度も回す、質問ごとに訪問済み・作業用の配列を作り直す（N × Q が大きいとき）
//
//  機能一覧（添字・値は 0〜n−1）                                                                 計算量
//    ClearableArray<T> a(n, def = T())    全部 def の配列                                             O(n)
//      a.get(i) / a.set(i, x)             読む（clear 後に書いていなければ def）/ 書く                   O(1)
//      a.ref(i)                           書き換え用の参照（a.ref(i) += x など）                        O(1)
//      a.is_set(i)                        clear の後に set / ref したか                                  O(1)
//      a.clear() / a.size()               全部 def に戻す / n                                            O(1)
//    VisitedStamp vis(n)                  訪問済みの印だけ（bool 配列の代わり）
//      vis.visit(i)                       初めてなら印を付けて true、もう付いていれば false              O(1)
//      vis.visited(i) / vis.unvisit(i)    印が付いているか / 印を外す                                   O(1)
//      vis.clear() / vis.size()           全部の印を外す / n                                             O(1)
//    SparseSet st(n)                      0〜n−1 の整数の集合
//      st.insert(x) / st.erase(x)         入れる / 消す（変わったら true）                              O(1)
//      st.contains(x) / st.size() / st.empty()                                                         O(1)
//      st[k] / for (int x : st)           中身（並びは入れた順。erase で最後の要素がその場所に移る）    O(1) / O(要素数)
//      st.pop_back() / st.to_vector()     最後の要素を取り出す / 中身の vector                           O(1) / O(要素数)
//      st.clear() / st.capacity()         全部消す / n                                                   O(1)
//
//  速さの目安: n = 10^6 で「VisitedStamp の clear + 100 個に印」を 10^6 回 約 0.13 秒
//              （vector<char> を毎回 fill すると 10^4 回で 約 0.16 秒 → 10^6 回なら 約 16 秒）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    clr_test::random_test()  … vector と std::set を使う愚直解との突き合わせ
//    clr_test::benchmark()    … 速度計測（毎回 fill する場合との比較つき）
//    clr_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
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
// ---- 全部消すのが O(1) の配列（世代番号つき） ----
template <class T>
struct ClearableArray {
    explicit ClearableArray(int n = 0, T def = T()) : def_(def), val_(n, def), st_(n, 0) {}
    T get(int i) const { return st_[i] == cur_ ? val_[i] : def_; }
    void set(int i, const T& x) { val_[i] = x, st_[i] = cur_; }
    T& ref(int i) {  // 書き換え用の参照（clear 後に初めてなら def にしてから返す）
        if (st_[i] != cur_) st_[i] = cur_, val_[i] = def_;
        return val_[i];
    }
    bool is_set(int i) const { return st_[i] == cur_; }
    void clear() {  // 世代を進めるだけ（約 40 億回に 1 回だけ全部書き直す）
        if (++cur_ == 0) std::fill(st_.begin(), st_.end(), 0u), cur_ = 1;
    }
    int size() const { return (int)val_.size(); }

   private:
    T def_;
    std::vector<T> val_;
    std::vector<unsigned> st_;
    unsigned cur_ = 1;
};

// ---- 訪問済みの印（bool 配列の代わり。全部外すのが O(1)） ----
struct VisitedStamp {
    explicit VisitedStamp(int n = 0) : st_(n, 0) {}
    bool visit(int i) {  // 初めてなら印を付けて true
        if (st_[i] == cur_) return false;
        st_[i] = cur_;
        return true;
    }
    bool visited(int i) const { return st_[i] == cur_; }
    void unvisit(int i) { st_[i] = 0; }  // 世代は 1 以上なので 0 は「印なし」
    void clear() {
        if (++cur_ == 0) std::fill(st_.begin(), st_.end(), 0u), cur_ = 1;
    }
    int size() const { return (int)st_.size(); }

   private:
    std::vector<unsigned> st_;
    unsigned cur_ = 1;
};

// ---- 0〜n−1 の整数の集合（追加・削除・含むか・全消去が O(1)、中身を O(要素数) で列挙） ----
struct SparseSet {
    explicit SparseSet(int n = 0) : dense_(n), pos_(n, 0) {}
    bool contains(int x) const {
        const int p = pos_[x];
        return p < sz_ && dense_[p] == x;
    }
    bool insert(int x) {
        if (contains(x)) return false;
        dense_[sz_] = x, pos_[x] = sz_++;
        return true;
    }
    bool erase(int x) {  // 最後の要素を x の場所へ移す
        if (!contains(x)) return false;
        const int p = pos_[x], y = dense_[--sz_];
        dense_[p] = y, pos_[y] = p;
        return true;
    }
    int pop_back() { return dense_[--sz_]; }  // 最後の要素を取り出す（空なら使わない）
    void clear() { sz_ = 0; }
    int size() const { return sz_; }
    bool empty() const { return sz_ == 0; }
    int operator[](int k) const { return dense_[k]; }
    const int* begin() const { return dense_.data(); }
    const int* end() const { return dense_.data() + sz_; }
    std::vector<int> to_vector() const { return std::vector<int>(begin(), end()); }
    int capacity() const { return (int)dense_.size(); }

   private:
    std::vector<int> dense_, pos_;
    int sz_ = 0;
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace clr_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define CLR_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_clearable_array] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    for (int it = 0; it < 500 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 20);
        ClearableArray<ll> ca(n, -7);
        std::vector<ll> ba(n, -7);
        std::vector<char> bset(n, 0);
        VisitedStamp vs(n);
        std::vector<char> bv(n, 0);
        SparseSet ss(n);
        std::set<int> bs;
        for (int t = 0; t < 200; ++t) {
            const int k = (int)rnd(rng, 0, 9), i = (int)rnd(rng, 0, n - 1);
            const ll x = rnd(rng, -9, 9);
            if (k == 0) {
                ca.set(i, x), ba[i] = x, bset[i] = 1;
            } else if (k == 1) {
                ca.ref(i) += x, ba[i] += x, bset[i] = 1;
            } else if (k == 2 && rnd(rng, 0, 5) == 0) {
                ca.clear(), std::fill(ba.begin(), ba.end(), -7), std::fill(bset.begin(), bset.end(), 0);
                vs.clear(), std::fill(bv.begin(), bv.end(), 0);
                ss.clear(), bs.clear();
            } else if (k == 3) {
                CLR_CHECK(vs.visit(i) == !bv[i]);
                bv[i] = 1;
            } else if (k == 4) {
                vs.unvisit(i), bv[i] = 0;
            } else if (k == 5) {
                CLR_CHECK(ss.insert(i) == !bs.count(i));
                bs.insert(i);
            } else if (k == 6) {
                CLR_CHECK(ss.erase(i) == (bs.count(i) > 0));
                bs.erase(i);
            } else if (k == 7 && !ss.empty()) {
                const int y = ss.pop_back();
                CLR_CHECK(bs.count(y) > 0);
                bs.erase(y);
            }
            for (int j = 0; j < n; ++j) {
                CLR_CHECK(ca.get(j) == ba[j] && ca.is_set(j) == (bool)bset[j]);
                CLR_CHECK(vs.visited(j) == (bool)bv[j]);
                CLR_CHECK(ss.contains(j) == (bs.count(j) > 0));
            }
            std::vector<int> got = ss.to_vector(), got2;
            for (int v : ss) got2.push_back(v);
            for (int q = 0; q < ss.size(); ++q) CLR_CHECK(ss[q] == got[q]);
            std::sort(got.begin(), got.end());
            CLR_CHECK(got == std::vector<int>(bs.begin(), bs.end()) && got2.size() == got.size());
            CLR_CHECK(ss.size() == (int)bs.size() && ss.empty() == bs.empty());
        }
        CLR_CHECK(ca.size() == n && vs.size() == n && ss.capacity() == n);
        ++cases;
    }
    // clear の直後は初期値
    {
        ClearableArray<int> c(2, 5);
        c.set(0, 1);
        c.clear();
        CLR_CHECK(c.get(0) == 5 && !c.is_set(0));
    }
    std::cout << "[ds_clearable_array] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_clearable_array] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = 1000000;
    std::vector<int> pts(100);
    for (auto& p : pts) p = (int)rnd(rng, 0, N - 1);
    {
        const int R = TEST_SCALE == 1 ? 100000 : 1000000;
        const auto t0 = std::chrono::steady_clock::now();
        VisitedStamp vis(N);
        ll chk = 0;
        for (int r = 0; r < R; ++r) {
            vis.clear();
            for (int p : pts) chk += vis.visit((p + r) % N);
        }
        std::cout << "[ds_clearable_array] benchmark: n = " << N << " で VisitedStamp の clear + 100 個に印 を " << R << " 回 " << elapsed_ms(t0) << " ms  (check "
                  << chk % 1000 << ")\n";
    }
    {
        const int R = TEST_SCALE == 1 ? 1000 : 10000;
        const auto t0 = std::chrono::steady_clock::now();
        std::vector<char> vis(N);
        ll chk = 0;
        for (int r = 0; r < R; ++r) {
            std::fill(vis.begin(), vis.end(), 0);
            for (int p : pts) {
                const int q = (p + r) % N;
                chk += !vis[q], vis[q] = 1;
            }
        }
        std::cout << "[ds_clearable_array] benchmark: （比較）vector<char> を毎回 fill して同じことを " << R << " 回 " << elapsed_ms(t0) << " ms  (check "
                  << chk % 1000 << ")\n";
    }
}

}  // namespace clr_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace clr_demo {

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
#define SHOW(expr, comment) clr_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    clr_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("【制限付き最速版】O(1) で全部消せる配列・集合（ds_clearable_array） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ 制限: 添字・値は 0〜n−1 の整数だけ・順序はない・範囲外は確かめない\n"
              << "■ できること: 配列 / 訪問済みの印 / 整数の集合 を、全消去 O(1) 付きで（集合の操作も O(1)）\n"
              << "■ API 一覧\n"
              << "  ClearableArray<T> a(n, def): get / set / ref / is_set / clear / size\n"
              << "  VisitedStamp vis(n): visit / visited / unvisit / clear / size\n"
              << "  SparseSet st(n): insert / erase / contains / size / empty / [k] / for / pop_back / to_vector / clear / capacity\n";

    api("【1】 ClearableArray<T> a(n, def)  get / set / ref / is_set / clear", "T / void / T& / bool / void", "読む / 書く / 書き換え用の参照 / clear 後に書いたか / 全部 def に戻す",
        "すべて O(1)");
    RUN(ClearableArray<int> a(5, -1));
    RUN(a.set(2, 10); a.ref(3) += 4);
    SHOW(a.get(2), "");
    SHOW(a.get(3), "ref は def（−1）から足す → 3");
    SHOW(a.get(0), "書いていない → def");
    SHOW(a.is_set(3), "");
    RUN(a.clear());
    SHOW(a.get(2), "clear は O(1) で全部 def に戻る");
    SHOW(a.is_set(3), "");

    api("【2】 VisitedStamp vis(n)  visit / visited / unvisit / clear", "bool / bool / void / void", "初めてなら印を付けて true / 印があるか / 印を外す / 全部外す",
        "すべて O(1)");
    RUN(VisitedStamp vis(6));
    SHOW(vis.visit(4), "初めて → true");
    SHOW(vis.visit(4), "2 回目 → false");
    SHOW(vis.visited(4), "");
    RUN(vis.unvisit(4));
    SHOW(vis.visited(4), "");
    RUN(vis.visit(1); vis.clear());
    SHOW(vis.visited(1), "clear で全部外れる");
    {
        code("// 使い方の例：各頂点から距離 1 以内の頂点数を、訪問済みを使い回して数える（毎回 fill しない）");
        const std::vector<std::vector<int>> g = {{1, 2}, {0, 2}, {0, 1, 3}, {2}};
        VisitedStamp v((int)g.size());
        std::vector<int> cnt;
        for (int s = 0; s < (int)g.size(); ++s) {
            v.clear();
            int c = v.visit(s);
            for (int t : g[s]) c += v.visit(t);
            cnt.push_back(c);
        }
        note("答え → " + to_s(cnt));
    }

    api("【3】 SparseSet st(n)  insert / erase / contains / 列挙 / pop_back / clear", "bool / bool / bool / … / int / void",
        "入れる / 消す（変わったら true）/ 含むか / 中身を入れた順に（erase で並びが変わる）/ 最後の要素を取り出す / 全部消す", "O(1)（列挙は O(要素数)）");
    RUN(SparseSet st(10));
    RUN(st.insert(7); st.insert(2); st.insert(5); st.insert(9));
    SHOW(st.to_vector(), "入れた順");
    SHOW(st.erase(2), "2 を消す → 最後の 9 が 2 の場所へ");
    SHOW(st.to_vector(), "");
    SHOW(st.contains(5), "");
    SHOW(st.insert(5), "もう入っている → false");
    SHOW(st[0], "0 番目");
    SHOW(st.pop_back(), "最後の要素を取り出す");
    SHOW(st.size(), "");
    RUN(st.clear());
    SHOW(st.empty(), "clear は O(1)");
    SHOW(st.capacity(), "");
}

}  // namespace clr_demo

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
    clr_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    clr_test::benchmark();    // 速度計測
    clr_demo::demo();         // 使い方デモ
    return 0;
}
