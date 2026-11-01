// =============================================================================
//  永続配列・永続スタック・永続キュー（過去のどの版にも戻れる。更新すると新しい版の番号が返る）   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  「元に戻す」「ある時点を保存して、あとで呼び出す」「版 v から分かれて別々に更新する」を、
//  毎回コピーせずにできる。更新は新しい版を作り、古い版はそのまま残る
//    例: ABC273 E（Notebook）… 末尾に追加・末尾を消す・保存・呼び出し → PersistentStack
//
//  機能一覧（版は int の番号。0 が最初の版）                                                     計算量
//    PersistentArray<T> pa(a)           配列 a を版 0 にする（長さは変わらない）                        O(N)
//      pa.set(v, i, x)                  版 v の a[i] を x にした新しい版の番号                          O(log N)
//      pa.get(v, i)                     版 v の a[i]                                                   O(log N)
//      pa.to_vector(v) / pa.size()      版 v の配列全体 / 長さ N                                        O(N) / O(1)
//    PersistentStack<T> ps              版 0 = 空のスタック
//      ps.push(v, x) / ps.pop(v)        版 v に x を積んだ版 / 版 v の一番上を除いた版                  O(1)
//      ps.top(v) / ps.size(v) / ps.empty(v)                                                           O(1)
//    PersistentQueue<T> pq              版 0 = 空のキュー
//      pq.push(v, x) / pq.pop(v)        版 v の後ろに x を入れた版 / 版 v の先頭を除いた版              O(1) / O(1)
//      pq.front(v) / pq.back(v) / pq.size(v) / pq.empty(v)                                            O(log N) / O(1)
//
//  注意: 版の番号は更新のたびに 1 つずつ増える（pop も新しい番号を返す）。空のときに pop・top・front すると止める
//  速さの目安: PersistentArray（N = 2×10^5）で set と get を 2×10^5 回ずつ（版はばらばら）約 0.3 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    persist_test::random_test()  … 版ごとに丸ごとコピーして持つ愚直解との突き合わせ
//    persist_test::benchmark()    … 速度計測
//    persist_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace persist_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace persist_internal

template <class T>
struct PersistentArray {
    PersistentArray() {}
    explicit PersistentArray(const std::vector<T>& a) : n_((int)a.size()) {
        sz_ = 1;
        while (sz_ < n_) sz_ <<= 1;
        nd_.reserve((size_t)2 * sz_);
        roots_.push_back(build(a, 0, sz_));
    }
    int size() const { return n_; }
    int versions() const { return (int)roots_.size(); }
    int set(int v, int i, const T& x) {
        check(v, i, "PersistentArray::set");
        roots_.push_back(set_rec(roots_[v], 0, sz_, i, x));
        return (int)roots_.size() - 1;
    }
    T get(int v, int i) const {
        check(v, i, "PersistentArray::get");
        int cur = roots_[v], lo = 0, hi = sz_;
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            if (i < mid) cur = nd_[cur].l, hi = mid;
            else cur = nd_[cur].r, lo = mid;
        }
        return nd_[cur].val;
    }
    std::vector<T> to_vector(int v) const {
        std::vector<T> out(n_);
        for (int i = 0; i < n_; ++i) out[i] = get(v, i);
        return out;
    }

   private:
    struct Node {
        int l, r;
        T val;
    };
    int n_ = 0, sz_ = 1;
    std::vector<Node> nd_;
    std::vector<int> roots_;
    int new_node(int l, int r, const T& x) {
        nd_.push_back({l, r, x});
        return (int)nd_.size() - 1;
    }
    int set_rec(int node, int lo, int hi, int i, const T& x) {  // 根から葉までの道だけコピーする
        if (hi - lo == 1) return new_node(-1, -1, x);
        const int mid = (lo + hi) / 2, l = nd_[node].l, r = nd_[node].r;
        if (i < mid) {
            const int nl = set_rec(l, lo, mid, i, x);
            return new_node(nl, r, T{});
        }
        const int nr = set_rec(r, mid, hi, i, x);
        return new_node(l, nr, T{});
    }
    int build(const std::vector<T>& a, int lo, int hi) {
        if (hi - lo == 1) return new_node(-1, -1, lo < n_ ? a[lo] : T{});
        const int mid = (lo + hi) / 2;
        const int l = build(a, lo, mid), r = build(a, mid, hi);
        return new_node(l, r, T{});
    }
    void check(int v, int i, const char* f) const {
        if (v < 0 || v >= (int)roots_.size()) persist_internal::fail(f, "版の番号が範囲外");
        if (i < 0 || i >= n_) persist_internal::fail(f, "i が範囲外");
    }
};

template <class T>
struct PersistentStack {
    PersistentStack() { node_.push_back({T{}, -1, 0}); }  // 版 0 = 空（番兵の頂点 0）
    int push(int v, const T& x) {
        check(v, "PersistentStack::push");
        const int h = ver_[v];
        node_.push_back({x, h, node_[h].size + 1});
        ver_.push_back((int)node_.size() - 1);
        return (int)ver_.size() - 1;
    }
    int pop(int v) {
        check(v, "PersistentStack::pop");
        if (empty(v)) persist_internal::fail("PersistentStack::pop", "空");
        ver_.push_back(node_[ver_[v]].parent);
        return (int)ver_.size() - 1;
    }
    const T& top(int v) const {
        check(v, "PersistentStack::top");
        if (empty(v)) persist_internal::fail("PersistentStack::top", "空");
        return node_[ver_[v]].x;
    }
    int size(int v) const {
        check(v, "PersistentStack::size");
        return node_[ver_[v]].size;
    }
    bool empty(int v) const { return size(v) == 0; }
    int versions() const { return (int)ver_.size(); }

   private:
    struct Node {
        T x;
        int parent, size;
    };
    std::vector<Node> node_;
    std::vector<int> ver_{0};  // 版 → 一番上の頂点
    void check(int v, const char* f) const {
        if (v < 0 || v >= (int)ver_.size()) persist_internal::fail(f, "版の番号が範囲外");
    }
};

// 入れた要素を木（親 = 1 つ前に入れた要素）にし、版 = (最後に入れた頂点, 長さ)。先頭は「深さ = 最後の深さ − 長さ + 1」の祖先
template <class T>
struct PersistentQueue {
    PersistentQueue() { node_.push_back({T{}, 0, 0, 0}); }  // 頂点 0 = 根（深さ 0）
    int push(int v, const T& x) {
        check(v, "PersistentQueue::push");
        const int p = ver_[v].first;
        const int pj = node_[p].jump, pjj = node_[pj].jump;
        const int jump = (node_[p].depth - node_[pj].depth == node_[pj].depth - node_[pjj].depth) ? pjj : p;
        node_.push_back({x, p, jump, node_[p].depth + 1});
        ver_.push_back({(int)node_.size() - 1, ver_[v].second + 1});
        return (int)ver_.size() - 1;
    }
    int pop(int v) {
        check(v, "PersistentQueue::pop");
        if (empty(v)) persist_internal::fail("PersistentQueue::pop", "空");
        ver_.push_back({ver_[v].first, ver_[v].second - 1});
        return (int)ver_.size() - 1;
    }
    const T& front(int v) const {
        check(v, "PersistentQueue::front");
        if (empty(v)) persist_internal::fail("PersistentQueue::front", "空");
        const int last = ver_[v].first;
        return node_[ancestor(last, node_[last].depth - ver_[v].second + 1)].x;
    }
    const T& back(int v) const {
        check(v, "PersistentQueue::back");
        if (empty(v)) persist_internal::fail("PersistentQueue::back", "空");
        return node_[ver_[v].first].x;
    }
    int size(int v) const {
        check(v, "PersistentQueue::size");
        return ver_[v].second;
    }
    bool empty(int v) const { return size(v) == 0; }
    int versions() const { return (int)ver_.size(); }

   private:
    struct Node {
        T x;
        int parent, jump, depth;
    };
    std::vector<Node> node_;
    std::vector<std::pair<int, int>> ver_{{0, 0}};
    int ancestor(int u, int d) const {  // 深さ d の祖先（ジャンプポインタで O(log N)）
        while (node_[u].depth > d) u = node_[node_[u].jump].depth >= d ? node_[u].jump : node_[u].parent;
        return u;
    }
    void check(int v, const char* f) const {
        if (v < 0 || v >= (int)ver_.size()) persist_internal::fail(f, "版の番号が範囲外");
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace persist_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define PA_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_persistent_array] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    std::mt19937_64 rng(20261031);
    long long cases = 0;
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        // 配列
        const int n = (int)rnd(rng, 1, 13);
        std::vector<ll> a(n);
        for (auto& x : a) x = rnd(rng, -9, 9);
        PersistentArray<ll> pa(a);
        std::vector<std::vector<ll>> va{a};
        // スタック・キュー
        PersistentStack<ll> ps;
        PersistentQueue<ll> pq;
        std::vector<std::vector<ll>> vs{{}};
        std::vector<std::deque<ll>> vq{{}};
        for (int t = 0; t < 200; ++t) {
            {
                const int v = (int)rnd(rng, 0, (int)va.size() - 1);
                if (rnd(rng, 0, 1)) {
                    const int i = (int)rnd(rng, 0, n - 1);
                    const ll x = rnd(rng, -9, 9);
                    PA_CHECK(pa.set(v, i, x) == (int)va.size());
                    va.push_back(va[v]), va.back()[i] = x;
                } else {
                    const int i = (int)rnd(rng, 0, n - 1);
                    PA_CHECK(pa.get(v, i) == va[v][i]);
                }
            }
            {
                const int v = (int)rnd(rng, 0, (int)vs.size() - 1);
                if (vs[v].empty() || rnd(rng, 0, 2)) {
                    const ll x = rnd(rng, -99, 99);
                    PA_CHECK(ps.push(v, x) == (int)vs.size());
                    vs.push_back(vs[v]), vs.back().push_back(x);
                } else {
                    PA_CHECK(ps.top(v) == vs[v].back());
                    PA_CHECK(ps.pop(v) == (int)vs.size());
                    vs.push_back(vs[v]), vs.back().pop_back();
                }
                const int w = (int)vs.size() - 1;
                PA_CHECK(ps.size(w) == (int)vs[w].size() && ps.empty(w) == vs[w].empty());
            }
            {
                const int v = (int)rnd(rng, 0, (int)vq.size() - 1);
                if (vq[v].empty() || rnd(rng, 0, 2)) {
                    const ll x = rnd(rng, -99, 99);
                    PA_CHECK(pq.push(v, x) == (int)vq.size());
                    vq.push_back(vq[v]), vq.back().push_back(x);
                } else {
                    PA_CHECK(pq.front(v) == vq[v].front() && pq.back(v) == vq[v].back());
                    PA_CHECK(pq.pop(v) == (int)vq.size());
                    vq.push_back(vq[v]), vq.back().pop_front();
                }
                const int w = (int)vq.size() - 1;
                PA_CHECK(pq.size(w) == (int)vq[w].size());
                if (!vq[w].empty()) PA_CHECK(pq.front(w) == vq[w].front() && pq.back(w) == vq[w].back());
            }
            cases += 3;
        }
        for (int v = 0; v < (int)va.size(); v += 7) PA_CHECK(pa.to_vector(v) == va[v]);
    }
    // 長いキュー（ジャンプポインタの深さ）
    {
        PersistentQueue<int> q;
        int v = 0;
        for (int i = 0; i < 100000; ++i) v = q.push(v, i);
        for (int i = 0; i < 99990; ++i) v = q.pop(v);
        PA_CHECK(q.front(v) == 99990 && q.back(v) == 99999 && q.size(v) == 10);
        ++cases;
    }
    std::cout << "[ds_persistent_array] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_persistent_array] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 20000 : 200000;
    std::vector<ll> a(N, 0);
    const auto t0 = std::chrono::steady_clock::now();
    PersistentArray<ll> pa(a);
    ll chk = 0;
    for (int i = 0; i < N; ++i) {
        const int v = (int)rnd(rng, 0, pa.versions() - 1);
        pa.set(v, (int)rnd(rng, 0, N - 1), i);
        chk += pa.get((int)rnd(rng, 0, pa.versions() - 1), (int)rnd(rng, 0, N - 1));
    }
    std::cout << "[ds_persistent_array] benchmark: PersistentArray（N = " << N << "）で set と get を " << N << " 回ずつ " << elapsed_ms(t0) << " ms  (check "
              << chk % 1000 << ")\n";
}

}  // namespace persist_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace persist_demo {

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
#define SHOW(expr, comment) persist_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    persist_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("永続配列・永続スタック・永続キュー（ds_persistent_array） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  PersistentArray<T> pa(a): set(v, i, x) → 新しい版 / get(v, i) / to_vector(v) / size / versions\n"
              << "  PersistentStack<T>: push(v, x) / pop(v) → 新しい版 / top(v) / size(v) / empty(v)\n"
              << "  PersistentQueue<T>: push(v, x) / pop(v) → 新しい版 / front(v) / back(v) / size(v) / empty(v)\n";

    api("【1】 PersistentArray<T> pa(a) / set(v, i, x) / get(v, i)", "— / int（新しい版）/ T", "版 v を元に 1 か所変えた新しい版を作る。古い版も読める",
        "O(N) / O(log N) / O(log N)");
    RUN(PersistentArray<int> pa(std::vector<int>{1, 2, 3}));
    RUN(int v1 = pa.set(0, 1, 20));
    RUN(int v2 = pa.set(v1, 2, 30));
    RUN(int v3 = pa.set(0, 0, 100));  // 版 0 から分かれる
    SHOW(pa.to_vector(0), "版 0 は変わらない");
    SHOW(pa.to_vector(v2), "");
    SHOW(pa.to_vector(v3), "");

    api("【2】 PersistentStack<T>：ABC273 E（Notebook）の形", "—",
        "ADD x = push、DELETE = pop（空なら何もしない）、SAVE y = page[y] = 今の版、LOAD z = 今の版 = page[z]、答え = top（空なら −1）", "1 回 O(1)");
    {
        PersistentStack<int> ps;
        std::vector<int> page(10, 0);
        int cur = 0;
        std::vector<int> out;
        const std::vector<std::pair<std::string, int>> ops = {{"ADD", 3}, {"SAVE", 2}, {"ADD", 4}, {"SAVE", 3}, {"LOAD", 2}, {"DELETE", 0}, {"DELETE", 0},
                                                              {"ADD", 10}, {"LOAD", 3}, {"ADD", 5}};
        for (auto [op, x] : ops) {
            if (op == "ADD") cur = ps.push(cur, x);
            else if (op == "DELETE") {
                if (!ps.empty(cur)) cur = ps.pop(cur);
            } else if (op == "SAVE") page[x] = cur;
            else cur = page[x];
            out.push_back(ps.empty(cur) ? -1 : ps.top(cur));
        }
        note("操作 ADD 3, SAVE 2, ADD 4, SAVE 3, LOAD 2, DELETE, DELETE, ADD 10, LOAD 3, ADD 5");
        note("→ 毎回の末尾: " + to_s(out));
    }

    api("【3】 PersistentQueue<T> pq / push(v, x) / pop(v) / front(v)", "int / int / const T&", "キューを版ごとに残す。どの版からでも続けられる",
        "O(1) / O(1) / O(log N)");
    RUN(PersistentQueue<int> pq);
    RUN(int a = pq.push(0, 1); int b = pq.push(a, 2); int c = pq.push(b, 3));
    RUN(int d = pq.pop(c));
    RUN(int e = pq.push(b, 9));  // 版 b（1, 2）から分かれる
    SHOW(pq.front(d), "版 d = (2, 3)");
    SHOW(pq.size(d), "");
    SHOW(pq.back(e), "版 e = (1, 2, 9)");
    SHOW(pq.front(e), "");
}

}  // namespace persist_demo

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
    persist_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    persist_test::benchmark();    // 速度計測
    persist_demo::demo();         // 使い方デモ
    return 0;
}
