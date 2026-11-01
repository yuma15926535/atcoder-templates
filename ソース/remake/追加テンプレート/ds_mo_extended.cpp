// =============================================================================
//  Mo の拡張（更新つき Mo・木のパスの Mo・巻き戻し Mo（足すだけでよい））   --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）
//  文字コード     : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  普通の Mo（MoSortTemplate.txt）は「区間を 1 つ伸ばす・縮める」で区間の質問をまとめて答える。ここではその拡張 3 つ
//    更新つき Mo  : 1 点の更新がまざるとき、時刻を 3 つ目の軸にする（ブロック N^(2/3)）。O(N^(5/3))
//    木のパスの Mo: オイラーツアー（入る・出るで 2 回並べる）で、パス u–v を列の区間に直す（LCA は別に入れる）
//    巻き戻し Mo  : 消す操作が難しいとき（最大値・最頻値の回数など）、足すだけで済むように左端はブロックの端から毎回足して戻す
//
//  機能一覧（区間は 0-indexed・半開区間 [l, r)。out(i) は i 番目の質問のときに呼ばれるので、そこで答えを記録する）
//    mo_with_updates(n, queries, num_updates, add, erase, toggle, out)
//        queries[i] = {l, r, t}（t = その質問より前にある更新の数）。toggle(t, l, r) は t 番目の更新を「入れる / 戻す」
//        （値の入れ替えで書くと同じ関数でよい。今の区間 [l, r) に位置が入っていれば erase → 値を変える → add）   O(N^(5/3))
//    mo_tree_path(g, root, queries, add, erase, out)
//        queries[i] = {u, v}。パス u–v の頂点（両端を含む）がちょうど add された状態で out(i) を呼ぶ             O(N √Q)
//    mo_rollback(n, queries, add, snapshot, rollback, reset, out)
//        add(i) で位置 i を入れるだけ。snapshot() で今を覚え、rollback() で覚えたところへ戻す。reset() で空にする   O(N √Q)
//
//  注意: 関数の中で使う配列などは自分で持つ（下のデモのような形）。ブロックの大きさは自動で決める
//  速さの目安: N = 10^5 で mo_with_updates（質問・更新 5×10^4 回ずつ）約 0.4 秒、mo_rollback（質問 10^5 回・最頻値の回数）約 0.1 秒
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    moext_test::random_test()  … 区間・パスを毎回数え直す愚直解との突き合わせ
//    moext_test::benchmark()    … 速度計測
//    moext_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
namespace moext_internal {
[[noreturn]] inline void fail(const char* func, const char* msg) {
    std::cerr << func << ": " << msg << "\n";
    std::abort();
}
}  // namespace moext_internal

template <class Add, class Erase, class Toggle, class Out>
void mo_with_updates(int n, const std::vector<std::array<int, 3>>& queries, int num_updates, Add add, Erase erase, Toggle toggle, Out out) {
    const int q = (int)queries.size();
    for (const auto& x : queries)
        if (x[0] < 0 || x[1] > n || x[0] > x[1] || x[2] < 0 || x[2] > num_updates)
            moext_internal::fail("mo_with_updates", "質問は 0 ≤ l ≤ r ≤ n、0 ≤ t ≤ 更新の数");
    const int B = std::max(1, (int)std::cbrt((double)std::max(1, n) * std::max(1, n)));  // n^(2/3)
    std::vector<int> ord(q);
    for (int i = 0; i < q; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        const int la = queries[a][0] / B, lb = queries[b][0] / B;
        if (la != lb) return la < lb;
        const int ra = queries[a][1] / B, rb = queries[b][1] / B;
        if (ra != rb) return (la & 1) ? ra > rb : ra < rb;
        return ((la + ra) & 1) ? queries[a][2] > queries[b][2] : queries[a][2] < queries[b][2];
    });
    int l = 0, r = 0, t = 0;
    for (int i : ord) {
        const int ql = queries[i][0], qr = queries[i][1], qt = queries[i][2];
        while (l > ql) add(--l);
        while (r < qr) add(r++);
        while (l < ql) erase(l++);
        while (r > qr) erase(--r);
        while (t < qt) toggle(t++, l, r);
        while (t > qt) toggle(--t, l, r);
        out(i);
    }
}

template <class Add, class Erase, class Out>
void mo_tree_path(const std::vector<std::vector<int>>& g, int root, const std::vector<std::pair<int, int>>& queries, Add add, Erase erase, Out out) {
    const int n = (int)g.size();
    if (n == 0) return;
    // オイラーツアー（入る・出るで 2 回）と LCA（ダブリング）
    std::vector<int> tin(n), tout(n), euler(2 * n), depth(n, 0);
    int LOG = 1;
    while ((1 << LOG) < n) ++LOG;
    std::vector<std::vector<int>> up(LOG, std::vector<int>(n, root));
    {
        std::vector<std::pair<int, int>> st{{root, 0}};
        std::vector<char> seen(n, 0);
        seen[root] = 1;
        int t = 0;
        while (!st.empty()) {
            auto& [v, k] = st.back();
            if (k == 0) tin[v] = t, euler[t++] = v;
            if (k < (int)g[v].size()) {
                const int c = g[v][k++];
                if (!seen[c]) seen[c] = 1, up[0][c] = v, depth[c] = depth[v] + 1, st.push_back({c, 0});
            } else {
                tout[v] = t, euler[t++] = v;
                st.pop_back();
            }
        }
        if (t != 2 * n) moext_internal::fail("mo_tree_path", "g は連結な木にする");
        for (int k = 1; k < LOG; ++k)
            for (int v = 0; v < n; ++v) up[k][v] = up[k - 1][up[k - 1][v]];
    }
    auto lca = [&](int a, int b) {
        if (depth[a] < depth[b]) std::swap(a, b);
        for (int k = LOG - 1; k >= 0; --k)
            if (depth[a] - (1 << k) >= depth[b]) a = up[k][a];
        if (a == b) return a;
        for (int k = LOG - 1; k >= 0; --k)
            if (up[k][a] != up[k][b]) a = up[k][a], b = up[k][b];
        return up[0][a];
    };
    // パス u–v を列の区間 [L, R) に（u が v の祖先なら [tin u, tin v]、そうでなければ [tout u, tin v] + LCA）
    const int q = (int)queries.size();
    std::vector<std::array<int, 3>> rg(q);  // {L, R, lca か −1}
    for (int i = 0; i < q; ++i) {
        int u = queries[i].first, v = queries[i].second;
        if (u < 0 || u >= n || v < 0 || v >= n) moext_internal::fail("mo_tree_path", "頂点番号が範囲外");
        if (tin[u] > tin[v]) std::swap(u, v);
        const int w = lca(u, v);
        if (w == u) rg[i] = {tin[u], tin[v] + 1, -1};
        else rg[i] = {tout[u], tin[v] + 1, w};
    }
    const int m = 2 * n;
    const int B = std::max(1, (int)(m / std::sqrt((double)std::max(1, q))));
    std::vector<int> ord(q);
    for (int i = 0; i < q; ++i) ord[i] = i;
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        const int ba = rg[a][0] / B, bb = rg[b][0] / B;
        if (ba != bb) return ba < bb;
        return (ba & 1) ? rg[a][1] > rg[b][1] : rg[a][1] < rg[b][1];
    });
    std::vector<char> in(n, 0);
    auto flip = [&](int pos) {
        const int v = euler[pos];
        if (in[v]) erase(v);
        else add(v);
        in[v] ^= 1;
    };
    int l = 0, r = 0;
    for (int i : ord) {
        while (l > rg[i][0]) flip(--l);
        while (r < rg[i][1]) flip(r++);
        while (l < rg[i][0]) flip(l++);
        while (r > rg[i][1]) flip(--r);
        if (rg[i][2] >= 0) add(rg[i][2]);
        out(i);
        if (rg[i][2] >= 0) erase(rg[i][2]);
    }
}

template <class Add, class Snapshot, class Rollback, class Reset, class Out>
void mo_rollback(int n, const std::vector<std::pair<int, int>>& queries, Add add, Snapshot snapshot, Rollback rollback, Reset reset, Out out) {
    const int q = (int)queries.size();
    for (const auto& x : queries)
        if (x.first < 0 || x.second > n || x.first > x.second) moext_internal::fail("mo_rollback", "質問は 0 ≤ l ≤ r ≤ n の半開区間");
    const int B = std::max(1, (int)(n / std::sqrt((double)std::max(1, q))));
    std::vector<int> ord;
    // 短い区間（同じブロックの中）はそのたびに作る
    for (int i = 0; i < q; ++i) {
        const auto [l, r] = queries[i];
        if (r - l <= B) {
            reset();
            for (int k = l; k < r; ++k) add(k);
            out(i);
        } else {
            ord.push_back(i);
        }
    }
    std::sort(ord.begin(), ord.end(), [&](int a, int b) {
        const int ba = queries[a].first / B, bb = queries[b].first / B;
        return ba != bb ? ba < bb : queries[a].second < queries[b].second;
    });
    int cur_block = -1, r = 0, bend = 0;
    for (int i : ord) {
        const auto [ql, qr] = queries[i];
        const int b = ql / B;
        if (b != cur_block) {
            cur_block = b, bend = std::min(n, (b + 1) * B);
            reset();
            r = bend;
        }
        while (r < qr) add(r++);  // 右はブロックの端から伸ばすだけ
        snapshot();
        for (int k = bend - 1; k >= ql; --k) add(k);  // 左は毎回足して、あとで戻す
        out(i);
        rollback();
    }
}
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト・速度計測
// =============================================================================
namespace moext_test {

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif

#define MO_CHECK(...)                                                                          \
    do {                                                                                       \
        if (!(__VA_ARGS__)) {                                                                  \
            std::cerr << "[ds_mo_extended] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
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
    // 更新つき：区間の種類数（1 点の書き換えがまざる）
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 25), V = (int)rnd(rng, 1, 6);
        std::vector<int> a(n), init;
        for (auto& x : a) x = (int)rnd(rng, 0, V);
        init = a;
        std::vector<std::pair<int, int>> upd;  // {位置, 新しい値}
        std::vector<std::array<int, 3>> qs;
        std::vector<int> exp;
        std::vector<int> cur = a;
        for (int t = 0; t < 40; ++t) {
            if (rnd(rng, 0, 2) == 0) {
                const int p = (int)rnd(rng, 0, n - 1), x = (int)rnd(rng, 0, V);
                upd.push_back({p, x}), cur[p] = x;
            } else {
                int l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
                if (l > r) std::swap(l, r);
                qs.push_back({l, r, (int)upd.size()});
                std::vector<int> seen(V + 1, 0);
                int d = 0;
                for (int i = l; i < r; ++i) d += !seen[cur[i]]++;
                exp.push_back(d);
            }
        }
        std::vector<int> cnt(V + 1, 0), got(qs.size(), -1);
        int kinds = 0;
        a = init;
        auto upd2 = upd;  // 入れ替えで持つ（toggle すると元に戻る）
        mo_with_updates(
            n, qs, (int)upd2.size(), [&](int i) { kinds += cnt[a[i]]++ == 0; }, [&](int i) { kinds -= --cnt[a[i]] == 0; },
            [&](int t, int l, int r) {
                const int p = upd2[t].first;
                const bool inside = l <= p && p < r;
                if (inside) kinds -= --cnt[a[p]] == 0;
                std::swap(a[p], upd2[t].second);
                if (inside) kinds += cnt[a[p]]++ == 0;
            },
            [&](int i) { got[i] = kinds; });
        MO_CHECK(got == exp);
        cases += (long long)qs.size();
    }
    // 木のパス：パス上の色の種類数
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 1, 25), V = (int)rnd(rng, 1, 5);
        std::vector<std::vector<int>> g(n);
        std::vector<int> par(n, -1);
        for (int i = 1; i < n; ++i) {
            const int p = (int)rnd(rng, 0, i - 1);
            g[i].push_back(p), g[p].push_back(i), par[i] = p;
        }
        std::vector<int> col(n);
        for (auto& c : col) c = (int)rnd(rng, 0, V);
        std::vector<std::pair<int, int>> qs(rnd(rng, 1, 30));
        for (auto& [u, v] : qs) u = (int)rnd(rng, 0, n - 1), v = (int)rnd(rng, 0, n - 1);
        std::vector<int> cnt(V + 1, 0), got(qs.size());
        int kinds = 0;
        mo_tree_path(
            g, 0, qs, [&](int v) { kinds += cnt[col[v]]++ == 0; }, [&](int v) { kinds -= --cnt[col[v]] == 0; }, [&](int i) { got[i] = kinds; });
        for (size_t i = 0; i < qs.size(); ++i) {
            // 愚直：根からの祖先の列で LCA を探してパスを集める
            auto anc = [&](int x) {
                std::vector<int> v;
                for (; x != -1; x = par[x]) v.push_back(x);
                return v;
            };
            const auto au = anc(qs[i].first), av = anc(qs[i].second);
            std::vector<char> onu(n, 0);
            for (int x : au) onu[x] = 1;
            int w = -1;
            for (int x : av)
                if (onu[x]) {
                    w = x;
                    break;
                }
            std::vector<int> seen(V + 1, 0);
            int d = 0;
            for (int x : au) {
                d += !seen[col[x]]++;
                if (x == w) break;
            }
            for (int x : av) {
                if (x == w) break;
                d += !seen[col[x]]++;
            }
            MO_CHECK(got[i] == d);
            ++cases;
        }
    }
    // 巻き戻し：区間の最頻値の回数（足すだけ）
    for (int it = 0; it < 300 * TEST_SCALE; ++it) {
        const int n = (int)rnd(rng, 0, 60), V = (int)rnd(rng, 1, 6);
        std::vector<int> a(n);
        for (auto& x : a) x = (int)rnd(rng, 0, V);
        std::vector<std::pair<int, int>> qs(rnd(rng, 1, 40));
        for (auto& [l, r] : qs) {
            l = (int)rnd(rng, 0, n), r = (int)rnd(rng, 0, n);
            if (l > r) std::swap(l, r);
        }
        std::vector<int> cnt(V + 1, 0), got(qs.size(), -1);
        int best = 0;
        std::vector<std::pair<int, int>> hist;  // (値, 足す前の best)
        size_t mark = 0;
        mo_rollback(
            n, qs,
            [&](int i) {
                hist.push_back({a[i], best});
                best = std::max(best, ++cnt[a[i]]);
            },
            [&]() { mark = hist.size(); },
            [&]() {
                while (hist.size() > mark) {
                    --cnt[hist.back().first];
                    best = hist.back().second;
                    hist.pop_back();
                }
            },
            [&]() {
                std::fill(cnt.begin(), cnt.end(), 0);
                best = 0, hist.clear(), mark = 0;
            },
            [&](int i) { got[i] = best; });
        for (size_t i = 0; i < qs.size(); ++i) {
            std::vector<int> c(V + 1, 0);
            int b = 0;
            for (int k = qs[i].first; k < qs[i].second; ++k) b = std::max(b, ++c[a[k]]);
            MO_CHECK(got[i] == b);
            ++cases;
        }
    }
    std::cout << "[ds_mo_extended] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    std::mt19937_64 rng(1);
    if (TEST_SCALE == 1) std::cout << "[ds_mo_extended] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    const int N = TEST_SCALE == 1 ? 10000 : 100000;
    std::vector<int> a(N);
    for (auto& x : a) x = (int)rnd(rng, 0, N - 1);
    // 更新つき
    std::vector<std::pair<int, int>> upd;
    std::vector<std::array<int, 3>> qs;
    for (int i = 0; i < N; ++i) {
        if (i & 1) upd.push_back({(int)rnd(rng, 0, N - 1), (int)rnd(rng, 0, N - 1)});
        else {
            int l = (int)rnd(rng, 0, N), r = (int)rnd(rng, 0, N);
            if (l > r) std::swap(l, r);
            qs.push_back({l, r, (int)upd.size()});
        }
    }
    std::vector<int> cnt(N, 0);
    ll kinds = 0, chk = 0;
    auto t0 = std::chrono::steady_clock::now();
    mo_with_updates(
        N, qs, (int)upd.size(), [&](int i) { kinds += cnt[a[i]]++ == 0; }, [&](int i) { kinds -= --cnt[a[i]] == 0; },
        [&](int t, int l, int r) {
            const int p = upd[t].first;
            const bool in = l <= p && p < r;
            if (in) kinds -= --cnt[a[p]] == 0;
            std::swap(a[p], upd[t].second);
            if (in) kinds += cnt[a[p]]++ == 0;
        },
        [&](int) { chk += kinds; });
    const ll t1 = elapsed_ms(t0);
    // 巻き戻し
    std::vector<std::pair<int, int>> q2(N);
    for (auto& [l, r] : q2) {
        l = (int)rnd(rng, 0, N), r = (int)rnd(rng, 0, N);
        if (l > r) std::swap(l, r);
    }
    std::fill(cnt.begin(), cnt.end(), 0);
    int best = 0;
    std::vector<std::pair<int, int>> hist;
    size_t mark = 0;
    t0 = std::chrono::steady_clock::now();
    mo_rollback(
        N, q2,
        [&](int i) {
            hist.push_back({a[i], best});
            best = std::max(best, ++cnt[a[i]]);
        },
        [&]() { mark = hist.size(); },
        [&]() {
            while (hist.size() > mark) --cnt[hist.back().first], best = hist.back().second, hist.pop_back();
        },
        [&]() {
            for (auto& h : hist) --cnt[h.first];
            best = 0, hist.clear(), mark = 0;
        },
        [&](int) { chk += best; });
    const ll t2 = elapsed_ms(t0);
    std::cout << "[ds_mo_extended] benchmark (N = " << N << "): mo_with_updates（質問・更新 " << N / 2 << " 回ずつ） " << t1 << " ms / mo_rollback（質問 " << N
              << " 回） " << t2 << " ms  (check " << chk % 1000 << ")\n";
}

}  // namespace moext_test

// =============================================================================
//  使い方デモ
// =============================================================================
namespace moext_demo {

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
#define SHOW(expr, comment) moext_demo::show(#expr, (expr), comment)
#define RUN(...)                  \
    moext_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void demo() {
    banner("Mo の拡張（ds_mo_extended） 使い方デモ    見方:  コード  ->  返り値    // 意味");
    std::cout << "\n■ API 一覧\n"
              << "  mo_with_updates(n, queries{l, r, t}, num_updates, add, erase, toggle, out)\n"
              << "  mo_tree_path(g, root, queries{u, v}, add, erase, out)\n"
              << "  mo_rollback(n, queries{l, r}, add, snapshot, rollback, reset, out)\n";

    api("【1】 mo_with_updates", "void", "更新がまざる区間の種類数。toggle(t, l, r) は値の入れ替えで書く（同じ関数で入れる・戻すになる）", "O(N^(5/3))");
    {
        std::vector<int> a = {1, 2, 1, 3, 2};
        std::vector<std::pair<int, int>> upd = {{2, 5}};               // 更新 0：a[2] = 5
        const std::vector<std::array<int, 3>> qs = {{0, 5, 0}, {0, 5, 1}, {1, 3, 1}};  // {l, r, その前の更新の数}
        std::vector<int> cnt(6, 0), ans(3);
        int kinds = 0;
        mo_with_updates(
            5, qs, 1, [&](int i) { kinds += cnt[a[i]]++ == 0; }, [&](int i) { kinds -= --cnt[a[i]] == 0; },
            [&](int t, int l, int r) {
                const int p = upd[t].first;
                const bool in = l <= p && p < r;
                if (in) kinds -= --cnt[a[p]] == 0;
                std::swap(a[p], upd[t].second);
                if (in) kinds += cnt[a[p]]++ == 0;
            },
            [&](int i) { ans[i] = kinds; });
        note("a = {1, 2, 1, 3, 2}、更新 a[2] = 5 → 質問 [0,5) 更新前 / [0,5) 更新後 / [1,3) 更新後 の種類数 = " + to_s(ans));
    }

    api("【2】 mo_tree_path", "void", "パス u–v（両端を含む）の頂点だけが add された状態で out(i)", "O(N √Q)");
    {
        const std::vector<std::vector<int>> g = {{1, 2}, {0, 3}, {0}, {1}};  // 0-1, 0-2, 1-3
        const std::vector<int> col = {1, 2, 1, 3};
        std::vector<int> cnt(4, 0), ans(2);
        int kinds = 0;
        mo_tree_path(
            g, 0, {{3, 2}, {1, 3}}, [&](int v) { kinds += cnt[col[v]]++ == 0; }, [&](int v) { kinds -= --cnt[col[v]] == 0; }, [&](int i) { ans[i] = kinds; });
        note("木 0-1, 0-2, 1-3、色 {1, 2, 1, 3} → パス 3–2（3,1,0,2）、1–3 の色の種類数 = " + to_s(ans));
    }

    api("【3】 mo_rollback", "void", "消すのが難しいとき（区間の最頻値の回数など）。add だけ書き、snapshot / rollback は「足した履歴を戻す」で書く",
        "O(N √Q)");
    code("// add(i): 履歴に (a[i], 今の best) を積んで cnt[a[i]]++ / snapshot(): mark = 履歴の長さ / rollback(): mark まで戻す / reset(): 空に");
    {
        const std::vector<int> a = {1, 2, 2, 3, 3, 3, 1};
        std::vector<int> cnt(4, 0), ans(3);
        int best = 0;
        std::vector<std::pair<int, int>> hist;
        size_t mark = 0;
        mo_rollback(
            7, {{0, 3}, {2, 7}, {0, 7}},
            [&](int i) {
                hist.push_back({a[i], best});
                best = std::max(best, ++cnt[a[i]]);
            },
            [&]() { mark = hist.size(); },
            [&]() {
                while (hist.size() > mark) --cnt[hist.back().first], best = hist.back().second, hist.pop_back();
            },
            [&]() {
                std::fill(cnt.begin(), cnt.end(), 0);
                best = 0, hist.clear(), mark = 0;
            },
            [&](int i) { ans[i] = best; });
        note("a = {1, 2, 2, 3, 3, 3, 1} → [0,3)、[2,7)、[0,7) の最頻値の回数 = " + to_s(ans));
    }
}

}  // namespace moext_demo

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
    moext_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    moext_test::benchmark();    // 速度計測
    moext_demo::demo();         // 使い方デモ
    return 0;
}
