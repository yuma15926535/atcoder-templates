// =============================================================================
//  グリッドの変換（回転・反転・切り取り・形の一致・点の回転）  grid_transform  --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 10 グリッドの変換
//  依存       : なし（ほかのテンプレートと一緒に貼っても衝突しない）
//  対象       : vector<string>（Grid なら g.s）と vector<vector<T>> の両方。結果は新しいグリッド（元は変えない）
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  回転の向きは見た目のとおり（時計回り = 右に 90 度）。H×W を 90 度回すと W×H になる。
//
//  機能一覧                                                                       計算量
//    rotate_cw(a) / rotate_ccw(a) / rotate180(a)    時計回り 90 度 / 反時計回り 90 度 / 180 度   O(HW)
//    rotate_k(a, k) / all_rotations(a)              時計回りに k 回（負も可）/ 4 通り全部          O(HW)
//    transpose(a) / flip_lr(a) / flip_ud(a)         転置 / 左右反転 / 上下反転                    O(HW)
//    shift_torus(a, di, dj)                         下に di・右に dj ずらす（はみ出しは反対側へ。ABC300 B）  O(HW)
//    trim(a, empty)                                 empty 以外を含む最小の長方形に切り取る          O(HW)
//    match_rotation(a, b, empty)                    回転と平行移動で b に重なる最小の回数 k（なければ -1。ABC218 C）  O(HW)
//    match_symmetry(a, b, empty)                    裏返しも許す（4 ~ 7 は左右反転してから回転）       O(HW)
//    find_subgrid(a, pat [, wildcard])              pat と一致する左上の位置 (i, j) の一覧（wildcard は何でも可）  O(HW hw)
//    rotate_point_cw / _ccw(i, j, ci, cj) / rotate_point(i, j, ci, cj, k)   点を中心 (ci, cj) のまわりに回す（ABC275 C）
//    manhattan_to_chebyshev(x, y) / chebyshev_to_manhattan(u, v)            45 度回転（マンハッタン距離 ↔ チェビシェフ距離）
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    gt_test::random_test() / gt_test::benchmark() / gt_demo::demo()
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_GRID_TRANSFORM
#define GG_GRID_TRANSFORM
// ---- グリッドの変換（グラフ・グリッド第二世代）  vector<string> と vector<vector<T>> のどちらにも使える ----
// 回転は見た目のとおり（時計回り = 右に 90 度）。H×W のグリッドを 90 度回すと W×H になる
namespace gg_internal {
template <class G>
G make_grid(int h, int w) {
    using Row = typename G::value_type;
    return G(h, Row(w, typename Row::value_type()));
}
template <class G>
int width_of(const G& a) { return a.empty() ? 0 : (int)a[0].size(); }
}  // namespace gg_internal
template <class G>
G rotate_cw(const G& a) {  // 時計回りに 90 度: (i, j) → (j, H-1-i)
    const int H = (int)a.size(), W = gg_internal::width_of(a);
    G r = gg_internal::make_grid<G>(W, H);
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) r[j][H - 1 - i] = a[i][j];
    return r;
}
template <class G>
G rotate_ccw(const G& a) {  // 反時計回りに 90 度: (i, j) → (W-1-j, i)
    const int H = (int)a.size(), W = gg_internal::width_of(a);
    G r = gg_internal::make_grid<G>(W, H);
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) r[W - 1 - j][i] = a[i][j];
    return r;
}
template <class G>
G rotate180(const G& a) {  // 180 度: (i, j) → (H-1-i, W-1-j)
    const int H = (int)a.size(), W = gg_internal::width_of(a);
    G r = gg_internal::make_grid<G>(H, W);
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) r[H - 1 - i][W - 1 - j] = a[i][j];
    return r;
}
template <class G>
G rotate_k(const G& a, int k) {  // 時計回りに 90 度を k 回（負なら反時計回り）
    k = ((k % 4) + 4) % 4;
    return k == 0 ? a : k == 1 ? rotate_cw(a) : k == 2 ? rotate180(a) : rotate_ccw(a);
}
template <class G>
G transpose(const G& a) {  // 転置: (i, j) → (j, i)
    const int H = (int)a.size(), W = gg_internal::width_of(a);
    G r = gg_internal::make_grid<G>(W, H);
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j) r[j][i] = a[i][j];
    return r;
}
template <class G>
G flip_lr(const G& a) {  // 左右反転: (i, j) → (i, W-1-j)
    G r = a;
    for (auto& row : r) std::reverse(row.begin(), row.end());
    return r;
}
template <class G>
G flip_ud(const G& a) {  // 上下反転: (i, j) → (H-1-i, j)
    G r = a;
    std::reverse(r.begin(), r.end());
    return r;
}
template <class G>
std::vector<G> all_rotations(const G& a) {  // {そのまま, 時計回り 90, 180, 反時計回り 90}
    return {a, rotate_cw(a), rotate180(a), rotate_ccw(a)};
}
// トーラスとしてずらす: r[(i + di) mod H][(j + dj) mod W] = a[i][j]（下に di、右に dj。負も可。ABC300 B）
template <class G>
G shift_torus(const G& a, long long di, long long dj) {
    const int H = (int)a.size(), W = gg_internal::width_of(a);
    if (H == 0 || W == 0) return a;
    const int si = (int)((di % H + H) % H), sj = (int)((dj % W + W) % W);
    G r = gg_internal::make_grid<G>(H, W);
    for (int i = 0; i < H; ++i) {
        const int ni = (i + si) % H;
        for (int j = 0; j < W; ++j) r[ni][(j + sj) % W] = a[i][j];
    }
    return r;
}
// 余白の切り取り: empty 以外の値があるマスをすべて含む最小の長方形（全部 empty なら空のグリッド）
template <class G, class T>
G trim(const G& a, const T& empty) {
    const int H = (int)a.size(), W = gg_internal::width_of(a);
    int t = H, b = -1, l = W, rr = -1;
    for (int i = 0; i < H; ++i)
        for (int j = 0; j < W; ++j)
            if (!(a[i][j] == empty)) t = std::min(t, i), b = std::max(b, i), l = std::min(l, j), rr = std::max(rr, j);
    if (b < 0) return G();
    G r = gg_internal::make_grid<G>(b - t + 1, rr - l + 1);
    for (int i = t; i <= b; ++i)
        for (int j = l; j <= rr; ++j) r[i - t][j - l] = a[i][j];
    return r;
}
// 形の一致（平行移動は無視）: a を時計回りに k 回回すと b と同じ形になる最小の k（0 ~ 3）。なければ -1。ABC218 C
template <class G, class T>
int match_rotation(const G& a, const G& b, const T& empty) {
    const G tb = trim(b, empty);
    G x = a;
    for (int k = 0; k < 4; ++k, x = rotate_cw(x))
        if (trim(x, empty) == tb) return k;
    return -1;
}
// 裏返しも許す形の一致: 0 ~ 3 は match_rotation と同じ、4 ~ 7 は左右反転してから時計回りに (値 - 4) 回。なければ -1
template <class G, class T>
int match_symmetry(const G& a, const G& b, const T& empty) {
    const int k = match_rotation(a, b, empty);
    if (k >= 0) return k;
    const int f = match_rotation(flip_lr(a), b, empty);
    return f >= 0 ? 4 + f : -1;
}
// 部分一致: pat の左上を a の (i, j) に置くと全部のマスが一致する (i, j) の一覧（行優先）。O(HW hw)
template <class G>
std::vector<std::pair<int, int>> find_subgrid(const G& a, const G& pat) {
    const int H = (int)a.size(), W = gg_internal::width_of(a), h = (int)pat.size(), w = gg_internal::width_of(pat);
    std::vector<std::pair<int, int>> r;
    for (int i = 0; i + h <= H; ++i)
        for (int j = 0; j + w <= W; ++j) {
            bool ok = true;
            for (int x = 0; x < h && ok; ++x)
                for (int y = 0; y < w; ++y)
                    if (!(a[i + x][j + y] == pat[x][y])) {
                        ok = false;
                        break;
                    }
            if (ok) r.emplace_back(i, j);
        }
    return r;
}
// pat のうち wildcard のマスは何とでも一致する版（例: '?' や '.' を「何でもよい」にする）
template <class G, class T>
std::vector<std::pair<int, int>> find_subgrid(const G& a, const G& pat, const T& wildcard) {
    const int H = (int)a.size(), W = gg_internal::width_of(a), h = (int)pat.size(), w = gg_internal::width_of(pat);
    std::vector<std::pair<int, int>> r;
    for (int i = 0; i + h <= H; ++i)
        for (int j = 0; j + w <= W; ++j) {
            bool ok = true;
            for (int x = 0; x < h && ok; ++x)
                for (int y = 0; y < w; ++y)
                    if (!(pat[x][y] == wildcard) && !(a[i + x][j + y] == pat[x][y])) {
                        ok = false;
                        break;
                    }
            if (ok) r.emplace_back(i, j);
        }
    return r;
}
// 点の 90 度回転（グリッドの座標 (i, j)。中心 (ci, cj) のまわり、見た目で時計回り / 反時計回り）。
// 中心が半整数（マスの角など）のときは、座標を全部 2 倍してから使う
inline std::pair<long long, long long> rotate_point_cw(long long i, long long j, long long ci, long long cj) {
    return {ci + (j - cj), cj - (i - ci)};
}
inline std::pair<long long, long long> rotate_point_ccw(long long i, long long j, long long ci, long long cj) {
    return {ci - (j - cj), cj + (i - ci)};
}
inline std::pair<long long, long long> rotate_point(long long i, long long j, long long ci, long long cj, int k) {  // 時計回りに k 回
    k = ((k % 4) + 4) % 4;
    for (int t = 0; t < k; ++t) {
        const auto p = rotate_point_cw(i, j, ci, cj);
        i = p.first, j = p.second;
    }
    return {i, j};
}
// 45 度回転: マンハッタン距離 |dx| + |dy| が、変換後のチェビシェフ距離 max(|du|, |dv|) になる
inline std::pair<long long, long long> manhattan_to_chebyshev(long long x, long long y) { return {x + y, x - y}; }
inline std::pair<long long, long long> chebyshev_to_manhattan(long long u, long long v) { return {(u + v) / 2, (u - v) / 2}; }  // u, v は同じ偶奇
#endif  // GG_GRID_TRANSFORM
// ======================= ここまでテンプレート本体 =======================


static const char* TEST_TAG = "grid_transform";
// ---- テスト用の小道具 ----
#define CHECK(...)                                                                                \
    do {                                                                                          \
        if (!(__VA_ARGS__)) {                                                                     \
            std::cerr << "[" << TEST_TAG << "] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                         \
        }                                                                                         \
    } while (0)
// テスト量の倍率。デバッグビルド（最適化なし・範囲チェックあり）は遅いので自動で 1/8 にする
#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
[[maybe_unused]] static const int TEST_SCALE = 1;
#else
[[maybe_unused]] static const int TEST_SCALE = 8;
#endif
[[maybe_unused]] static int rnd(std::mt19937_64& rng, int lo, int hi) { return lo + (int)(rng() % (unsigned long long)(hi - lo + 1)); }
[[maybe_unused]] static long long ms_since(std::chrono::steady_clock::time_point t0) {
    return (long long)std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
struct CinFromString {  // std::cin を文字列に差し替える
    std::istringstream ss;
    std::streambuf* old;
    explicit CinFromString(const std::string& s) : ss(s), old(std::cin.rdbuf(ss.rdbuf())) {}
    ~CinFromString() { std::cin.rdbuf(old); }
};

namespace gt_test {

using VS = std::vector<std::string>;
using VV = std::vector<std::vector<int>>;

static VS random_vs(std::mt19937_64& rng, int H, int W, const char* chars) {
    VS a(H, std::string(W, '.'));
    const int k = (int)std::string(chars).size();
    for (auto& row : a)
        for (char& c : row) c = chars[rng() % k];
    return a;
}
static VS brute_trim(const VS& a, char e) {
    int t = 1 << 30, b = -1, l = 1 << 30, r = -1;
    for (int i = 0; i < (int)a.size(); ++i)
        for (int j = 0; j < (int)a[i].size(); ++j)
            if (a[i][j] != e) t = std::min(t, i), b = std::max(b, i), l = std::min(l, j), r = std::max(r, j);
    if (b < 0) return {};
    VS res;
    for (int i = t; i <= b; ++i) res.push_back(a[i].substr(l, r - l + 1));
    return res;
}
// 平行移動で重なるか（集合として比べる）
static bool same_shape(const VS& a, const VS& b, char e) {
    std::vector<std::tuple<int, int, char>> pa, pb;
    for (int i = 0; i < (int)a.size(); ++i)
        for (int j = 0; j < (int)a[i].size(); ++j)
            if (a[i][j] != e) pa.emplace_back(i, j, a[i][j]);
    for (int i = 0; i < (int)b.size(); ++i)
        for (int j = 0; j < (int)b[i].size(); ++j)
            if (b[i][j] != e) pb.emplace_back(i, j, b[i][j]);
    if (pa.size() != pb.size()) return false;
    if (pa.empty()) return true;
    const int di = std::get<0>(pb[0]) - std::get<0>(pa[0]), dj = std::get<1>(pb[0]) - std::get<1>(pa[0]);
    for (size_t k = 0; k < pa.size(); ++k)
        if (std::get<0>(pa[k]) + di != std::get<0>(pb[k]) || std::get<1>(pa[k]) + dj != std::get<1>(pb[k]) || std::get<2>(pa[k]) != std::get<2>(pb[k])) return false;
    return true;
}

static void random_test() {
    std::mt19937_64 rng(20261006);
    int cases = 0;
    for (int it = 0; it < 4000 * TEST_SCALE; ++it, ++cases) {
        const int H = rnd(rng, 1, 6), W = rnd(rng, 1, 6);
        const VS a = random_vs(rng, H, W, (it % 2) ? ".#" : "..#ab");
        VV v(H, std::vector<int>(W));
        for (auto& row : v)
            for (int& x : row) x = rnd(rng, -2, 2);
        // 定義どおりか
        const VS cw = rotate_cw(a), ccw = rotate_ccw(a), r180 = rotate180(a), tr = transpose(a), lr = flip_lr(a), ud = flip_ud(a);
        CHECK((int)cw.size() == W && (int)cw[0].size() == H && (int)ccw.size() == W && (int)tr.size() == W && (int)r180.size() == H);
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) {
                CHECK(cw[j][H - 1 - i] == a[i][j] && ccw[W - 1 - j][i] == a[i][j] && r180[H - 1 - i][W - 1 - j] == a[i][j]);
                CHECK(tr[j][i] == a[i][j] && lr[i][W - 1 - j] == a[i][j] && ud[H - 1 - i][j] == a[i][j]);
            }
        // 関係
        CHECK(rotate_cw(rotate_cw(rotate_cw(rotate_cw(a)))) == a && rotate_ccw(cw) == a && rotate_cw(rotate_cw(a)) == r180);
        CHECK(cw == flip_lr(tr) && ccw == flip_ud(tr) && transpose(tr) == a && flip_lr(lr) == a);
        const std::vector<VS> rots = all_rotations(a);
        CHECK(rots.size() == 4 && rots[0] == a && rots[1] == cw && rots[2] == r180 && rots[3] == ccw);
        for (int k = -6; k <= 6; ++k) CHECK(rotate_k(a, k) == rots[((k % 4) + 4) % 4]);
        // vector<vector<int>> でも同じ
        {
            const VV vcw = rotate_cw(v);
            for (int i = 0; i < H; ++i)
                for (int j = 0; j < W; ++j) CHECK(vcw[j][H - 1 - i] == v[i][j]);
            CHECK(rotate_ccw(vcw) == v && transpose(transpose(v)) == v && flip_ud(flip_ud(v)) == v && rotate180(rotate180(v)) == v);
            const VV vt = trim(v, 0);
            int cnt = 0, cnt2 = 0;
            for (auto& row : v)
                for (int x : row) cnt += (x != 0);
            for (auto& row : vt)
                for (int x : row) cnt2 += (x != 0);
            CHECK(cnt == cnt2);
        }
        // shift_torus
        {
            const long long di = rnd(rng, -20, 20), dj = rnd(rng, -20, 20);
            const VS sh = shift_torus(a, di, dj);
            for (int i = 0; i < H; ++i)
                for (int j = 0; j < W; ++j) CHECK(sh[(int)(((i + di) % H + H) % H)][(int)(((j + dj) % W + W) % W)] == a[i][j]);
            CHECK(shift_torus(a, H, -W) == a && shift_torus(sh, -di, -dj) == a);
        }
        // trim
        CHECK(trim(a, '.') == brute_trim(a, '.'));
        // match_rotation / match_symmetry: b は a を回転・反転して余白をつけたもの、ときどき 1 マス変える
        {
            const int k = rnd(rng, 0, 3);
            const bool fl = rng() % 3 == 0;
            VS b = rotate_k(fl ? flip_lr(a) : a, k);
            const int pt = rnd(rng, 0, 2), pl = rnd(rng, 0, 2), pb = rnd(rng, 0, 2), pr = rnd(rng, 0, 2);
            const int bw = (int)b[0].size() + pl + pr;
            VS c(pt, std::string(bw, '.'));
            for (auto& row : b) c.push_back(std::string(pl, '.') + row + std::string(pr, '.'));
            for (int t = 0; t < pb; ++t) c.push_back(std::string(bw, '.'));
            if (rng() % 4 == 0) {
                const int i = rnd(rng, 0, (int)c.size() - 1), j = rnd(rng, 0, bw - 1);
                c[i][j] = (c[i][j] == '.') ? '#' : '.';
            }
            int exr = -1, exs = -1;
            for (int t = 0; t < 4 && exr < 0; ++t)
                if (same_shape(rotate_k(a, t), c, '.')) exr = t;
            exs = exr;
            for (int t = 0; t < 4 && exs < 0; ++t)
                if (same_shape(rotate_k(flip_lr(a), t), c, '.')) exs = 4 + t;
            CHECK(match_rotation(a, c, '.') == exr && match_symmetry(a, c, '.') == exs);
        }
        // find_subgrid
        {
            const int h = rnd(rng, 1, H), w = rnd(rng, 1, W);
            VS pat = random_vs(rng, h, w, ".#");
            if (rng() % 2) {  // a の一部を切り出す
                const int i0 = rnd(rng, 0, H - h), j0 = rnd(rng, 0, W - w);
                for (int x = 0; x < h; ++x) pat[x] = a[i0 + x].substr(j0, w);
            }
            std::vector<std::pair<int, int>> ex, exw;
            for (int i = 0; i + h <= H; ++i)
                for (int j = 0; j + w <= W; ++j) {
                    bool ok = true, okw = true;
                    for (int x = 0; x < h; ++x)
                        for (int y = 0; y < w; ++y) ok &= (a[i + x][j + y] == pat[x][y]), okw &= (pat[x][y] == '.' || a[i + x][j + y] == pat[x][y]);
                    if (ok) ex.emplace_back(i, j);
                    if (okw) exw.emplace_back(i, j);
                }
            CHECK(find_subgrid(a, pat) == ex && find_subgrid(a, pat, '.') == exw);
        }
        // 点の回転: 正方形 N×N のグリッドの回転と一致する（中心は ((N-1)/2, (N-1)/2)。座標を 2 倍して整数にする）
        {
            const int N = std::min(H, W);
            VS sq(N);
            for (int i = 0; i < N; ++i) sq[i] = a[i].substr(0, N);
            const VS sc = rotate_cw(sq), scc = rotate_ccw(sq);
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    const auto p = rotate_point_cw(2 * i, 2 * j, N - 1, N - 1);
                    CHECK(p.first % 2 == 0 && p.second % 2 == 0 && sc[p.first / 2][p.second / 2] == sq[i][j]);
                    const auto q = rotate_point_ccw(2 * i, 2 * j, N - 1, N - 1);
                    CHECK(scc[q.first / 2][q.second / 2] == sq[i][j]);
                }
            const long long x = rnd(rng, -50, 50), y = rnd(rng, -50, 50), ci = rnd(rng, -9, 9), cj = rnd(rng, -9, 9);
            const auto p1 = rotate_point_cw(x, y, ci, cj);
            CHECK(rotate_point_ccw(p1.first, p1.second, ci, cj) == std::make_pair(x, y));
            for (int k = -5; k <= 5; ++k) {
                std::pair<long long, long long> e{x, y};
                for (int t = 0; t < ((k % 4) + 4) % 4; ++t) e = rotate_point_cw(e.first, e.second, ci, cj);
                CHECK(rotate_point(x, y, ci, cj, k) == e);
            }
            // 見た目で時計回り: 中心の右のマスは下へ、上のマスは右へ
            CHECK(rotate_point_cw(ci, cj + 1, ci, cj) == std::make_pair(ci + 1, cj) && rotate_point_cw(ci - 1, cj, ci, cj) == std::make_pair(ci, cj + 1));
        }
        // 45 度回転
        {
            const long long x1 = rnd(rng, -100, 100), y1 = rnd(rng, -100, 100), x2 = rnd(rng, -100, 100), y2 = rnd(rng, -100, 100);
            const auto a1 = manhattan_to_chebyshev(x1, y1), a2 = manhattan_to_chebyshev(x2, y2);
            CHECK(std::abs(x1 - x2) + std::abs(y1 - y2) == std::max(std::abs(a1.first - a2.first), std::abs(a1.second - a2.second)));
            CHECK(chebyshev_to_manhattan(a1.first, a1.second) == std::make_pair(x1, y1));
        }
    }
    {
        const VS e;
        CHECK(rotate_cw(e).empty() && trim(e, '.').empty() && trim(VS{"..", ".."}, '.').empty() && match_rotation(VS{".."}, VS{"."}, '.') == 0);
        CHECK(find_subgrid(VS{"ab"}, VS{"abc"}).empty());
    }
    std::cout << "[grid_transform] random_test: OK (" << cases << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[grid_transform] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {1000, 2000}) {
        if (TEST_SCALE == 1 && n > 1000) break;
        const VS a = random_vs(rng, n, n, ".#");
        VV v(n, std::vector<int>(n));
        for (auto& row : v)
            for (int& x : row) x = (int)(rng() % 10);
        auto t0 = std::chrono::steady_clock::now();
        const VS r = rotate_cw(a);
        const long long t1 = ms_since(t0);
        t0 = std::chrono::steady_clock::now();
        const VV rv = rotate_cw(v);
        const long long t2 = ms_since(t0);
        t0 = std::chrono::steady_clock::now();
        const int k = match_symmetry(a, rotate_ccw(flip_lr(a)), '.');
        std::cout << "[grid_transform] benchmark: " << n << " × " << n << "  rotate_cw(string) " << t1 << " ms, rotate_cw(int) " << t2 << " ms, match_symmetry " << ms_since(t0)
                  << " ms  (k = " << k << ", check " << (r[0][0] + rv[0][0]) << ")\n";
    }
}

}  // namespace gt_test


namespace gt_demo {
// ---- デモ用の小道具（見方: コード -> 返り値 // 意味） ----
template <class A, class B>
[[maybe_unused]] static std::string to_s(const std::pair<A, B>& p);
template <class T, size_t N>
[[maybe_unused]] static std::string to_s(const std::array<T, N>& a);
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v);
#ifdef GG_GRAPH_BASE
[[maybe_unused]] static std::string to_s(const Edge& e) { return "{to=" + std::to_string(e.to) + ", id=" + std::to_string(e.id) + ", w=" + std::to_string(e.w) + "}"; }
[[maybe_unused]] static std::string to_s(const RawEdge& e) { return "{" + std::to_string(e.from) + "-" + std::to_string(e.to) + ", w=" + std::to_string(e.w) + "}"; }
#endif
[[maybe_unused]] static std::string to_s(bool b) { return b ? "true" : "false"; }
[[maybe_unused]] static std::string to_s(char c) { return std::string("'") + c + "'"; }
[[maybe_unused]] static std::string to_s(const std::string& s) { return "\"" + s + "\""; }
[[maybe_unused]] static std::string to_s(const char* s) { return "\"" + std::string(s) + "\""; }
template <class T, std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool> && !std::is_same_v<T, char>, int> = 0>
[[maybe_unused]] static std::string to_s(T v) {
    if (std::is_same_v<T, long long> && (long long)v == (1LL << 60)) return "DIST_INF";
    return std::to_string(v);
}
template <class A, class B>
[[maybe_unused]] static std::string to_s(const std::pair<A, B>& p) { return "(" + to_s(p.first) + ", " + to_s(p.second) + ")"; }
template <class T, size_t N>
[[maybe_unused]] static std::string to_s(const std::array<T, N>& a) {
    std::string s = "(";
    for (size_t i = 0; i < N; ++i) s += (i ? ", " : "") + to_s(a[i]);
    return s + ")";
}
template <class T>
[[maybe_unused]] static std::string to_s(const std::vector<T>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
template <class T, class = decltype(std::declval<const T&>().to_vector())>
[[maybe_unused]] static std::string to_s(const T& v) { return to_s(v.to_vector()); }  // IntSpan など
[[maybe_unused]] static std::string pad(const std::string& s, size_t w) {  // 全角を 2 文字分として幅をそろえる
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
[[maybe_unused]] static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
[[maybe_unused]] static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
[[maybe_unused]] static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
[[maybe_unused]] static void stmt(const std::string& s) {
    const bool block = !s.empty() && s.back() == '}' &&
                       (s.compare(0, 4, "for ") == 0 || s.compare(0, 3, "if ") == 0 || s.compare(0, 6, "while ") == 0);
    code(block ? s : s + ";");
}
template <class V>
[[maybe_unused]] static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 36) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 56) + "  // " + comment;
    std::cout << line << "\n";
}
struct Row { const char* use; const char* ret; const char* what; };
template <size_t N>
[[maybe_unused]] static void cheat_sheet(const std::string& title, const Row (&rows)[N], size_t w1 = 46, size_t w2 = 22) {
    std::cout << "\n■ " << title << "\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, w1) << pad(ret, w2) << r.what << "\n";
    }
}
#define SHOW(expr, comment) show(#expr, (expr), comment)
#define RUN(...)             \
    stmt(#__VA_ARGS__);      \
    __VA_ARGS__

static void print_vs(const std::vector<std::string>& a) {
    for (const auto& row : a) std::cout << "      " << row << "\n";
}

static void demo() {
    banner("グリッドの変換 使い方デモ    見方:  コード  ->  返り値    // 意味");
    const Row rows[] = {
        {"rotate_cw(a) / rotate_ccw(a) / rotate180(a)", "同じ型", "時計回り 90 度 / 反時計回り 90 度 / 180 度"},
        {"rotate_k(a, k) / all_rotations(a)", "同じ型 / vector", "時計回りに k 回 / 4 通り"},
        {"transpose(a) / flip_lr(a) / flip_ud(a)", "同じ型", "転置 / 左右反転 / 上下反転"},
        {"shift_torus(a, di, dj)", "同じ型", "下に di、右に dj ずらす（はみ出しは反対側）"},
        {"trim(a, '.')", "同じ型", "'.' 以外を含む最小の長方形"},
        {"match_rotation(a, b, '.')", "int", "回転 + 平行移動で重なる最小の k（なければ -1）"},
        {"match_symmetry(a, b, '.')", "int", "裏返しも許す（4 ~ 7 は左右反転してから回転）"},
        {"find_subgrid(a, pat [, '?'])", "vector<pair>", "pat と一致する左上 (i, j) の一覧"},
        {"rotate_point_cw(i, j, ci, cj) など", "pair<ll, ll>", "点を (ci, cj) のまわりに 90 度回す"},
        {"manhattan_to_chebyshev(x, y)", "pair<ll, ll>", "(x + y, x - y)。|dx| + |dy| = max(|du|, |dv|)"},
    };
    cheat_sheet("API 一覧", rows, 46, 20);

    head("デモで使うグリッド（3 × 4）");
    RUN(std::vector<std::string> a = {"#...", "##..", "#.#."});
    print_vs(a);

    api("【1】 rotate_cw(a) / rotate_ccw(a) / rotate180(a) / rotate_k(a, k)", "同じ型（H×W → W×H）", "見た目のとおりに回す。(i, j) は時計回りで (j, H-1-i) へ", "O(HW)");
    RUN(auto r = rotate_cw(a));
    print_vs(r);
    SHOW(rotate_ccw(a), "反時計回り");
    SHOW(rotate180(a), "");
    SHOW(rotate_k(a, -1) == rotate_ccw(a), "負の k は反時計回り");
    note("Grid g なら rotate_cw(g.s)。Grid に戻すなら Grid(rotate_cw(g.s))。");

    api("【2】 transpose(a) / flip_lr(a) / flip_ud(a)", "同じ型", "転置（行と列を入れ替える）/ 左右反転 / 上下反転", "O(HW)");
    SHOW(transpose(a), "");
    SHOW(flip_lr(a), "");
    SHOW(flip_ud(a), "");

    api("【3】 vector<vector<T>> にも使える", "同じ型", "数値のグリッドも同じ関数で回す", "O(HW)");
    RUN(std::vector<std::vector<int>> m = {{1, 2, 3}, {4, 5, 6}});
    SHOW(rotate_cw(m), "");
    SHOW(transpose(m), "");

    api("【4】 shift_torus(a, di, dj)", "同じ型", "下に di、右に dj ずらす（はみ出した分は反対側から出てくる）。ABC300 B は全部の (di, dj) を試す", "O(HW)");
    SHOW(shift_torus(a, 1, 1), "");

    api("【5】 trim(a, empty) / match_rotation(a, b, empty) / match_symmetry(a, b, empty)", "同じ型 / int",
        "形の比較: 余白（empty）を切り取って、回転（と裏返し）で重なるか。平行移動は自由", "O(HW)");
    RUN(std::vector<std::string> b = {".....", ".###.", "..#..", ".#...", "....."});
    print_vs(b);
    SHOW(trim(b, '.'), "'.' だけの行・列を切り取る");
    SHOW(trim(a, '.'), "a の右端の列も切り取られる");
    SHOW(match_rotation(a, b, '.'), "a を時計回りに 1 回まわすと b と同じ形（重ならなければ -1）");
    RUN(std::vector<std::string> c = {"..#", ".##", "#.#"});
    SHOW(match_rotation(a, c, '.'), "回転だけでは重ならない");
    SHOW(match_symmetry(a, c, '.'), "裏返すと重なる（4 ~ 7 は左右反転してから (値 - 4) 回）");
    note("ABC218 C は match_rotation(S, T, '.') != -1 で Yes。");

    api("【6】 find_subgrid(a, pat) / find_subgrid(a, pat, wildcard)", "vector<pair<int,int>>", "pat を左上 (i, j) に置くと全部一致する位置。wildcard のマスは何でもよい", "O(HW hw)");
    SHOW(find_subgrid(a, std::vector<std::string>{"#", "#"}), "縦に # が 2 つ並ぶ位置");
    SHOW(find_subgrid(a, std::vector<std::string>{"#?", "??"}, '?'), "'?' は何でもよい");

    api("【7】 rotate_point_cw(i, j, ci, cj) / rotate_point_ccw / rotate_point(i, j, ci, cj, k)", "pair<long long, long long>",
        "点 (i, j) を中心 (ci, cj) のまわりに見た目で 90 度回す（グリッドの座標: i は下向き、j は右向き）", "O(1)");
    SHOW(rotate_point_cw(0, 1, 0, 0), "中心の右 → 中心の下");
    SHOW(rotate_point_cw(-1, 0, 0, 0), "中心の上 → 中心の右");
    SHOW(rotate_point(2, 3, 1, 1, 2), "180 度");
    note("中心がマスの角（半整数）のときは座標を全部 2 倍して使う。正方形 ABC275 C は 2 頂点から残りを回転で作る。");

    api("【8】 manhattan_to_chebyshev(x, y) / chebyshev_to_manhattan(u, v)", "pair<long long, long long>",
        "45 度回転 (x + y, x - y)。マンハッタン距離が max(|u の差|, |v の差|) になり、最遠点・範囲の計算が軸ごとに分かれる", "O(1)");
    SHOW(manhattan_to_chebyshev(3, 1), "");
    SHOW(chebyshev_to_manhattan(4, 2), "元に戻す（u, v は同じ偶奇）");
    note("「マンハッタン距離の最大」は u, v それぞれの (最大 - 最小) の大きいほう（ABC178 E）。");
}

}  // namespace gt_demo

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
    gt_test::random_test();
    gt_test::benchmark();
    gt_demo::demo();
    return 0;
}
