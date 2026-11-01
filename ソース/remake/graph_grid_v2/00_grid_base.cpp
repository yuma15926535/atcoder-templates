// =============================================================================
//  グリッド基盤（Grid / CharSet / 隣接 / 移動リスト / 入力）     --- グラフ・グリッド第二世代
// -----------------------------------------------------------------------------
//  分類       : 00 基盤（09 グリッド探索・状態付き探索はこれを同梱している）
//  依存       : なし
//  座標       : マス (i, j) は 0-indexed（i = 行 0 ~ H-1、j = 列 0 ~ W-1）。番号 id = i * W + j
//  依存ライブラリ : 標準ライブラリのみ（C++17 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  方向の番号（全テンプレート共通）
//    0 上 (U, ^)   1 右 (R, >)   2 下 (D, v)   3 左 (L, <)      ← ここまで 4 方向
//    4 左上        5 右上        6 右下        7 左下           ← 8 方向
//    Grid::DI[d], Grid::DJ[d] が方向 d の (行, 列) の変化量。Grid::DIR_CHAR[d] = "URDL"[d]
//
//  高速化のポイント
//    - マスは番号 id（int 1 つ）で扱う。探索の距離配列などを 1 次元にでき、キャッシュに乗りやすい
//    - for_each4(v, f) は割り算 1 回と比較 4 回だけ。f はテンプレートなのでインライン展開される
//    - 壁の判定は CharSet（256 個の bool の表）で O(1)
//
//  機能一覧                                                                       計算量
//    Grid g(H, W, fill = '.') / Grid g(rows)   H×W のグリッド（rows は vector<string>）  O(HW)
//    read_grid(H, W) / read_grid()            H 行の文字列を読む / 「H W」から読む       O(HW)
//    read_matrix<T>(H, W)                     数値のグリッド（vector<vector<T>>）を読む   O(HW)
//    g.H, g.W, g.s[i][j], g(i, j), g[i][j], g.at(id)   大きさ・文字の読み書き           O(1)
//    g.size() / g.in(i, j) / g.id(i, j) / g.pos(id) / g.row(id) / g.col(id)            O(1)
//    g.count(c) / g.find(c) / g.find_all(c) / g.mask(c)   文字（の集合）の個数・位置     O(HW)
//    g.neighbor(id, d)                        方向 d の隣（範囲外なら -1）               O(1)
//    g.for_each4(id, f) / g.for_each8(id, f)  隣のマスごとに f(u) か f(u, d)              O(1)
//    g.for_each4(i, j, f) / g.for_each8(i, j, f)   (i, j) 版。f(ni, nj) か f(ni, nj, d)  O(1)
//    g.for_each_offsets(id, offs, f)          任意の移動リスト。f(u) か f(u, k)           O(|offs|)
//    g.for_each4_torus / for_each8_torus / for_each_offsets_torus   トーラス版          O(1)
//    g.slide(id, d, walls)                    壁の手前まで滑ったマス                     O(距離)
//    g.wrap(i, j) / g.wrap_id(i, j)           トーラスの座標の正規化                     O(1)
//    g.padded(c = '#')                        外周を 1 マスずつ c で囲んだコピー          O(HW)
//    g.border_cells() / g.is_border(id)       外周のマス                                  O(HW) / O(1)
//    g.path_to_dirs(path, torus = false)      経路（番号の列）→ "URDL" の文字列          O(長さ)
//    Grid::dir_of(c) / reverse_dir(d) / turn_right(d) / turn_left(d)   方向の変換        O(1)
//    Grid::offsets_manhattan(D) / offsets_chebyshev(D) / offsets_euclid_sq(d2) / offsets_knight()
//                                             移動リストを作る                            O(個数)
//    CharSet walls("#X") / CharSet('#') / ~CharSet(".")   文字の集合（補集合も作れる）   O(1) で判定
//    cout << g                                各行を改行区切りで出力
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    grb_test::random_test()  … 愚直な実装と全機能を突き合わせるランダムテスト
//    grb_test::benchmark()    … 速度計測
//    grb_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// ======================= ここからテンプレート本体（コピー範囲） =======================
#ifndef GG_GRID_BASE
#define GG_GRID_BASE
// ---- グリッド基盤（グラフ・グリッド第二世代）  Grid / CharSet / 隣接 / 移動 / 入力 ----
struct CharSet {  // 文字の集合（壁の文字などを O(1) で判定）。'#' / "#X" / std::string から暗黙に作れる
    std::array<bool, 256> f{};
    CharSet() {}
    CharSet(char c) { f[(unsigned char)c] = true; }
    CharSet(const char* s) { for (; *s; ++s) f[(unsigned char)*s] = true; }
    CharSet(const std::string& s) { for (char c : s) f[(unsigned char)c] = true; }
    bool operator()(char c) const { return f[(unsigned char)c]; }
    CharSet operator~() const {  // 補集合（例: ~CharSet(".") =「'.' 以外すべて」）
        CharSet r;
        for (int i = 0; i < 256; ++i) r.f[i] = !f[i];
        return r;
    }
};

namespace gg_internal {
template <class F>
inline void call_cell(F& f, int u, int k) {  // f(u) と f(u, 方向 or 番号) の両方を受け付ける
    if constexpr (std::is_invocable_v<F&, int, int>) f(u, k);
    else f(u);
}
template <class F>
inline void call_ij(F& f, int i, int j, int k) {  // f(i, j) と f(i, j, 方向 or 番号) の両方を受け付ける
    if constexpr (std::is_invocable_v<F&, int, int, int>) f(i, j, k);
    else f(i, j);
}
}  // namespace gg_internal

struct Grid {
    int H = 0, W = 0;
    std::vector<std::string> s;  // s[i][j] = マス (i, j) の文字
    // 方向: 0 上, 1 右, 2 下, 3 左（ここまで 4 方向）, 4 左上, 5 右上, 6 右下, 7 左下（8 方向）
    static constexpr int DI[8] = {-1, 0, 1, 0, -1, -1, 1, 1};
    static constexpr int DJ[8] = {0, 1, 0, -1, -1, 1, 1, -1};
    static constexpr char DIR_CHAR[5] = "URDL";  // DIR_CHAR[d] = 方向 d の文字

    Grid() {}
    Grid(int h, int w, char fill = '.') : H(h), W(w), s(h, std::string(w, fill)) {}
    explicit Grid(std::vector<std::string> rows) : H((int)rows.size()), W(rows.empty() ? 0 : (int)rows[0].size()), s(std::move(rows)) {}

    // ---- 座標と番号（番号 id = i * W + j） ----
    int size() const { return H * W; }
    bool in(int i, int j) const { return (unsigned)i < (unsigned)H && (unsigned)j < (unsigned)W; }
    bool in(std::pair<int, int> p) const { return in(p.first, p.second); }
    int id(int i, int j) const { return i * W + j; }
    int id(std::pair<int, int> p) const { return p.first * W + p.second; }
    std::pair<int, int> pos(int v) const { return {v / W, v % W}; }
    int row(int v) const { return v / W; }
    int col(int v) const { return v % W; }
    char& operator()(int i, int j) { return s[i][j]; }
    char operator()(int i, int j) const { return s[i][j]; }
    char& at(int v) { return s[v / W][v % W]; }
    char at(int v) const { return s[v / W][v % W]; }
    std::string& operator[](int i) { return s[i]; }
    const std::string& operator[](int i) const { return s[i]; }

    // ---- 文字の検索 ----
    int count(const CharSet& c) const {
        int r = 0;
        for (const std::string& t : s)
            for (char x : t) r += c(x);
        return r;
    }
    int find(const CharSet& c) const {  // 行優先で最初のマスの番号（なければ -1）
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j)
                if (c(s[i][j])) return i * W + j;
        return -1;
    }
    std::vector<int> find_all(const CharSet& c) const {  // 当てはまるマスの番号すべて（行優先）
        std::vector<int> r;
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j)
                if (c(s[i][j])) r.push_back(i * W + j);
        return r;
    }
    std::vector<char> mask(const CharSet& c) const {  // m[id] = 1（当てはまる） / 0
        std::vector<char> m(size());
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) m[i * W + j] = c(s[i][j]);
        return m;
    }

    // ---- 隣接の列挙（f は f(u) か f(u, 方向) / (i, j) 版は f(ni, nj) か f(ni, nj, 方向)） ----
    int neighbor(int v, int d) const {  // 方向 d の隣のマス（範囲外なら -1）
        const int i = v / W + DI[d], j = v % W + DJ[d];
        return in(i, j) ? i * W + j : -1;
    }
    template <class F>
    void for_each4(int v, F&& f) const {
        const int i = v / W, j = v - i * W;
        if (i > 0) gg_internal::call_cell(f, v - W, 0);
        if (j + 1 < W) gg_internal::call_cell(f, v + 1, 1);
        if (i + 1 < H) gg_internal::call_cell(f, v + W, 2);
        if (j > 0) gg_internal::call_cell(f, v - 1, 3);
    }
    template <class F>
    void for_each8(int v, F&& f) const {
        const int i = v / W, j = v - i * W;
        for (int d = 0; d < 8; ++d) {
            const int ni = i + DI[d], nj = j + DJ[d];
            if (in(ni, nj)) gg_internal::call_cell(f, ni * W + nj, d);
        }
    }
    template <class F>
    void for_each4(int i, int j, F&& f) const {
        for (int d = 0; d < 4; ++d) {
            const int ni = i + DI[d], nj = j + DJ[d];
            if (in(ni, nj)) gg_internal::call_ij(f, ni, nj, d);
        }
    }
    template <class F>
    void for_each8(int i, int j, F&& f) const {
        for (int d = 0; d < 8; ++d) {
            const int ni = i + DI[d], nj = j + DJ[d];
            if (in(ni, nj)) gg_internal::call_ij(f, ni, nj, d);
        }
    }
    // 任意の移動リスト（桂馬跳び・決まった距離へのジャンプなど）。f(u) か f(u, k)（k = offs の何番目か）
    template <class F>
    void for_each_offsets(int v, const std::vector<std::pair<int, int>>& offs, F&& f) const {
        const int i = v / W, j = v - i * W;
        for (int k = 0; k < (int)offs.size(); ++k) {
            const int ni = i + offs[k].first, nj = j + offs[k].second;
            if (in(ni, nj)) gg_internal::call_cell(f, ni * W + nj, k);
        }
    }
    template <class F>
    void for_each_offsets(int i, int j, const std::vector<std::pair<int, int>>& offs, F&& f) const {
        for (int k = 0; k < (int)offs.size(); ++k) {
            const int ni = i + offs[k].first, nj = j + offs[k].second;
            if (in(ni, nj)) gg_internal::call_ij(f, ni, nj, k);
        }
    }
    // トーラス（上下・左右の端がつながっている）版
    template <class F>
    void for_each4_torus(int v, F&& f) const {
        const int i = v / W, j = v - i * W;
        for (int d = 0; d < 4; ++d) gg_internal::call_cell(f, wrap_id(i + DI[d], j + DJ[d]), d);
    }
    template <class F>
    void for_each8_torus(int v, F&& f) const {
        const int i = v / W, j = v - i * W;
        for (int d = 0; d < 8; ++d) gg_internal::call_cell(f, wrap_id(i + DI[d], j + DJ[d]), d);
    }
    template <class F>
    void for_each_offsets_torus(int v, const std::vector<std::pair<int, int>>& offs, F&& f) const {
        const int i = v / W, j = v - i * W;
        for (int k = 0; k < (int)offs.size(); ++k) gg_internal::call_cell(f, wrap_id(i + offs[k].first, j + offs[k].second), k);
    }
    // 方向 d に、壁かグリッドの外の手前まで滑ったマス（すぐ前が壁なら v のまま）
    int slide(int v, int d, const CharSet& walls) const {
        int i = v / W, j = v % W;
        while (true) {
            const int ni = i + DI[d], nj = j + DJ[d];
            if (!in(ni, nj) || walls(s[ni][nj])) break;
            i = ni, j = nj;
        }
        return i * W + j;
    }
    // トーラスの座標の正規化（負の数・H 以上でもよい）
    std::pair<int, int> wrap(int i, int j) const { return {(i % H + H) % H, (j % W + W) % W}; }
    int wrap_id(int i, int j) const { return ((i % H + H) % H) * W + (j % W + W) % W; }

    // ---- 形を変えたコピー・外周 ----
    Grid padded(char c = '#') const {  // 外周を c で 1 マスずつ囲んだ (H+2)×(W+2)。元の (i, j) は (i+1, j+1) になる
        Grid g(H + 2, W + 2, c);
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) g.s[i + 1][j + 1] = s[i][j];
        return g;
    }
    bool is_border(int v) const {
        const int i = v / W, j = v % W;
        return i == 0 || j == 0 || i == H - 1 || j == W - 1;
    }
    std::vector<int> border_cells() const {  // 外周のマス（行優先・重複なし）。外から塗る多始点 BFS の始点に
        std::vector<int> r;
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j)
                if (i == 0 || j == 0 || i == H - 1 || j == W - 1) r.push_back(i * W + j);
        return r;
    }

    // ---- 方向 ----
    // 経路（マスの番号の列）を "URDL" の文字列にする。同じ行・列の何マス先への移動も向きの文字になる。
    // 斜めなどは '?'。torus = true なら端をまたぐ 1 マスの移動も正しい向きにする
    std::string path_to_dirs(const std::vector<int>& path, bool torus = false) const {
        std::string r;
        for (size_t k = 0; k + 1 < path.size(); ++k) {
            int di = path[k + 1] / W - path[k] / W, dj = path[k + 1] % W - path[k] % W;
            if (torus) {
                if (H > 2 && di == H - 1) di = -1;
                else if (H > 2 && di == -(H - 1)) di = 1;
                if (W > 2 && dj == W - 1) dj = -1;
                else if (W > 2 && dj == -(W - 1)) dj = 1;
            }
            const int si = (di > 0) - (di < 0), sj = (dj > 0) - (dj < 0);
            r += (si == -1 && sj == 0) ? 'U' : (si == 0 && sj == 1) ? 'R' : (si == 1 && sj == 0) ? 'D' : (si == 0 && sj == -1) ? 'L' : '?';
        }
        return r;
    }
    static int dir_of(char c) {  // 'U','R','D','L'（小文字・矢印 ^ > v < も可）→ 0 ~ 3。それ以外は -1
        switch (c) {
            case 'U': case 'u': case '^': return 0;
            case 'R': case 'r': case '>': return 1;
            case 'D': case 'd': case 'v': return 2;
            case 'L': case 'l': case '<': return 3;
            default: return -1;
        }
    }
    static int reverse_dir(int d) { return d < 4 ? (d + 2) & 3 : ((d - 2) & 3) + 4; }  // 逆向き（8 方向にも対応）
    static int turn_right(int d) { return (d + 1) & 3; }                                // 時計回りに 90 度
    static int turn_left(int d) { return (d + 3) & 3; }                                 // 反時計回りに 90 度

    // ---- 移動リストを作る（for_each_offsets に渡す） ----
    static std::vector<std::pair<int, int>> offsets_manhattan(int D, bool with_center = false) {  // |di| + |dj| <= D
        std::vector<std::pair<int, int>> r;
        for (int di = -D; di <= D; ++di) {
            const int rem = D - (di < 0 ? -di : di);
            for (int dj = -rem; dj <= rem; ++dj)
                if (with_center || di != 0 || dj != 0) r.emplace_back(di, dj);
        }
        return r;
    }
    static std::vector<std::pair<int, int>> offsets_chebyshev(int D, bool with_center = false) {  // max(|di|, |dj|) <= D
        std::vector<std::pair<int, int>> r;
        for (int di = -D; di <= D; ++di)
            for (int dj = -D; dj <= D; ++dj)
                if (with_center || di != 0 || dj != 0) r.emplace_back(di, dj);
        return r;
    }
    static std::vector<std::pair<int, int>> offsets_euclid_sq(long long d2) {  // di^2 + dj^2 == d2 ちょうど（ABC272 D 型）
        std::vector<std::pair<int, int>> r;
        if (d2 < 0) return r;
        auto isqrt = [](long long x) {
            long long y = (long long)std::sqrt((double)x);
            while (y * y > x) --y;
            while ((y + 1) * (y + 1) <= x) ++y;
            return y;
        };
        const long long R = isqrt(d2);
        for (long long di = -R; di <= R; ++di) {
            const long long rem = d2 - di * di, dj = isqrt(rem);
            if (dj * dj != rem) continue;
            if (dj == 0) r.emplace_back((int)di, 0);
            else r.emplace_back((int)di, (int)-dj), r.emplace_back((int)di, (int)dj);
        }
        return r;
    }
    static std::vector<std::pair<int, int>> offsets_knight() {  // 桂馬（チェスのナイト）の 8 方向
        return {{-2, -1}, {-2, 1}, {-1, -2}, {-1, 2}, {1, -2}, {1, 2}, {2, -1}, {2, 1}};
    }
};

// ---- 入力 ----
inline Grid read_grid(int H, int W) {  // H 行の文字列を読む
    Grid g(H, W);
    for (std::string& r : g.s) std::cin >> r;
    return g;
}
inline Grid read_grid() {  // 「H W」に続いて H 行を読む
    int H, W;
    std::cin >> H >> W;
    return read_grid(H, W);
}
template <class T>
inline std::vector<std::vector<T>> read_matrix(int H, int W) {  // 数値のグリッド（空白区切り）。形は Grid g(H, W) で扱える
    std::vector<std::vector<T>> a(H, std::vector<T>(W));
    for (auto& r : a)
        for (T& x : r) std::cin >> x;
    return a;
}
inline std::ostream& operator<<(std::ostream& os, const Grid& g) {  // 各行を改行区切りで出力（最後の行の後ろには改行なし）
    for (int i = 0; i < g.H; ++i) os << (i ? "\n" : "") << g.s[i];
    return os;
}
#endif  // GG_GRID_BASE
// ======================= ここまでテンプレート本体 =======================


struct CinFromString {
    std::istringstream ss;
    std::streambuf* old;
    explicit CinFromString(const std::string& s) : ss(s), old(std::cin.rdbuf(ss.rdbuf())) {}
    ~CinFromString() { std::cin.rdbuf(old); }
};

// =============================================================================
//  ランダムテスト
// =============================================================================
namespace grb_test {

#define GRB_CHECK(...)                                                                           \
    do {                                                                                     \
        if (!(__VA_ARGS__)) {                                                                \
            std::cerr << "[grid_base] TEST FAILED  line " << __LINE__ << ": " #__VA_ARGS__ "\n"; \
            std::exit(1);                                                                    \
        }                                                                                    \
    } while (0)

#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
static const int TEST_SCALE = 1;
#else
static const int TEST_SCALE = 8;
#endif

static int rnd(std::mt19937_64& rng, int lo, int hi) { return lo + (int)(rng() % (unsigned long long)(hi - lo + 1)); }
static int sgn(int x) { return (x > 0) - (x < 0); }

static void random_test() {
    std::mt19937_64 rng(20261006);
    const std::string alpha = ".#SGx";
    const int iters = 1500 * TEST_SCALE;
    for (int it = 0; it < iters; ++it) {
        const int H = rnd(rng, 1, 7), W = rnd(rng, 1, 7);
        std::vector<std::string> rows(H, std::string(W, '.'));
        for (auto& r : rows)
            for (char& c : r) c = alpha[rnd(rng, 0, (int)alpha.size() - 1)];
        const Grid g(rows);
        GRB_CHECK(g.H == H && g.W == W && g.size() == H * W && g.s == rows);
        // 1. 座標と番号
        for (int i = -2; i <= H + 1; ++i)
            for (int j = -2; j <= W + 1; ++j) GRB_CHECK(g.in(i, j) == (0 <= i && i < H && 0 <= j && j < W) && g.in({i, j}) == g.in(i, j));
        for (int v = 0; v < H * W; ++v) {
            const auto [i, j] = g.pos(v);
            GRB_CHECK(g.id(i, j) == v && g.id({i, j}) == v && g.row(v) == i && g.col(v) == j);
            GRB_CHECK(g.at(v) == rows[i][j] && g(i, j) == rows[i][j] && g[i][j] == rows[i][j]);
        }
        // 2. 隣接（番号版・(i, j) 版・方向付き・方向なし）
        for (int v = 0; v < H * W; ++v) {
            const int i = v / W, j = v % W;
            for (int K : {4, 8}) {
                std::vector<std::pair<int, int>> ex, got1, got2, got3, got4;
                for (int d = 0; d < K; ++d) {
                    const int ni = i + Grid::DI[d], nj = j + Grid::DJ[d];
                    if (0 <= ni && ni < H && 0 <= nj && nj < W) ex.emplace_back(ni * W + nj, d);
                    GRB_CHECK(g.neighbor(v, d) == ((0 <= ni && ni < H && 0 <= nj && nj < W) ? ni * W + nj : -1));
                }
                auto f_ud = [&](int u, int d) { got1.emplace_back(u, d); };
                auto f_u = [&](int u) { got2.emplace_back(u, -1); };
                auto f_ijd = [&](int a, int b, int d) { got3.emplace_back(a * W + b, d); };
                auto f_ij = [&](int a, int b) { got4.emplace_back(a * W + b, -1); };
                if (K == 4) g.for_each4(v, f_ud), g.for_each4(v, f_u), g.for_each4(i, j, f_ijd), g.for_each4(i, j, f_ij);
                else g.for_each8(v, f_ud), g.for_each8(v, f_u), g.for_each8(i, j, f_ijd), g.for_each8(i, j, f_ij);
                GRB_CHECK(got1 == ex && got3 == ex && got2.size() == ex.size() && got4.size() == ex.size());
                for (size_t k = 0; k < ex.size(); ++k) GRB_CHECK(got2[k].first == ex[k].first && got4[k].first == ex[k].first);
                // トーラス版
                std::vector<std::pair<int, int>> tex, tgot;
                for (int d = 0; d < K; ++d) {
                    const int ni = ((i + Grid::DI[d]) % H + H) % H, nj = ((j + Grid::DJ[d]) % W + W) % W;
                    tex.emplace_back(ni * W + nj, d);
                }
                if (K == 4) g.for_each4_torus(v, f_ud);
                else g.for_each8_torus(v, f_ud);
                tgot.assign(got1.begin() + (long)ex.size(), got1.end());
                GRB_CHECK(tgot == tex);
            }
        }
        // 3. 移動リスト
        {
            const int D = rnd(rng, 0, 3);
            const auto om = Grid::offsets_manhattan(D), omc = Grid::offsets_manhattan(D, true);
            const auto oc = Grid::offsets_chebyshev(D), occ = Grid::offsets_chebyshev(D, true);
            std::set<std::pair<int, int>> sm, smc, sc, scc;
            for (int a = -D; a <= D; ++a)
                for (int b = -D; b <= D; ++b) {
                    if (std::abs(a) + std::abs(b) <= D) {
                        smc.insert({a, b});
                        if (a || b) sm.insert({a, b});
                    }
                    scc.insert({a, b});
                    if (a || b) sc.insert({a, b});
                }
            GRB_CHECK(std::set<std::pair<int, int>>(om.begin(), om.end()) == sm && om.size() == sm.size());
            GRB_CHECK(std::set<std::pair<int, int>>(omc.begin(), omc.end()) == smc && omc.size() == smc.size());
            GRB_CHECK(std::set<std::pair<int, int>>(oc.begin(), oc.end()) == sc && oc.size() == sc.size());
            GRB_CHECK(std::set<std::pair<int, int>>(occ.begin(), occ.end()) == scc && occ.size() == scc.size());
            const long long d2 = rnd(rng, 0, 60);
            const auto oe = Grid::offsets_euclid_sq(d2);
            std::set<std::pair<int, int>> se;
            for (int a = -9; a <= 9; ++a)
                for (int b = -9; b <= 9; ++b)
                    if ((long long)a * a + (long long)b * b == d2) se.insert({a, b});
            GRB_CHECK(std::set<std::pair<int, int>>(oe.begin(), oe.end()) == se && oe.size() == se.size());
            GRB_CHECK(Grid::offsets_knight().size() == 8 && Grid::offsets_euclid_sq(5) == Grid::offsets_knight());
            const auto& offs = (it & 1) ? om : Grid::offsets_knight();
            for (int v = 0; v < H * W; ++v) {
                const int i = v / W, j = v % W;
                std::vector<std::pair<int, int>> ex, tex, got, got2, tgot;
                for (int k = 0; k < (int)offs.size(); ++k) {
                    const int ni = i + offs[k].first, nj = j + offs[k].second;
                    if (g.in(ni, nj)) ex.emplace_back(ni * W + nj, k);
                    tex.emplace_back(((ni % H + H) % H) * W + (nj % W + W) % W, k);
                }
                g.for_each_offsets(v, offs, [&](int u, int k) { got.emplace_back(u, k); });
                g.for_each_offsets(i, j, offs, [&](int a, int b, int k) { got2.emplace_back(a * W + b, k); });
                g.for_each_offsets_torus(v, offs, [&](int u, int k) { tgot.emplace_back(u, k); });
                GRB_CHECK(got == ex && got2 == ex && tgot == tex);
                int cnt = 0;
                g.for_each_offsets(v, offs, [&](int) { ++cnt; });
                GRB_CHECK(cnt == (int)ex.size());
            }
        }
        // 4. 文字の検索・CharSet
        {
            const CharSet cs = (it % 3 == 0) ? CharSet('#') : (it % 3 == 1) ? CharSet("#x") : CharSet(std::string("SG"));
            const CharSet ncs = ~cs;
            std::vector<int> ex;
            for (int v = 0; v < H * W; ++v) {
                GRB_CHECK(ncs(g.at(v)) == !cs(g.at(v)));
                if (cs(g.at(v))) ex.push_back(v);
            }
            GRB_CHECK(g.find_all(cs) == ex && g.count(cs) == (int)ex.size() && g.find(cs) == (ex.empty() ? -1 : ex[0]));
            const auto m = g.mask(cs);
            for (int v = 0; v < H * W; ++v) GRB_CHECK(m[v] == (char)cs(g.at(v)));
            GRB_CHECK(g.count('.') + g.count("#SGx") == H * W);
        }
        // 5. slide
        for (int v = 0; v < H * W; ++v)
            for (int d = 0; d < 8; ++d) {
                int i = v / W, j = v % W;
                while (g.in(i + Grid::DI[d], j + Grid::DJ[d]) && g(i + Grid::DI[d], j + Grid::DJ[d]) != '#') i += Grid::DI[d], j += Grid::DJ[d];
                GRB_CHECK(g.slide(v, d, '#') == i * W + j);
            }
        // 6. wrap
        for (int k = 0; k < 20; ++k) {
            const int i = rnd(rng, -30, 30), j = rnd(rng, -30, 30);
            const int ei = ((i % H) + H) % H, ej = ((j % W) + W) % W;
            GRB_CHECK(g.wrap(i, j) == std::make_pair(ei, ej) && g.wrap_id(i, j) == ei * W + ej);
        }
        // 7. padded / 外周
        {
            const Grid p = g.padded('@');
            GRB_CHECK(p.H == H + 2 && p.W == W + 2);
            for (int i = 0; i < H + 2; ++i)
                for (int j = 0; j < W + 2; ++j)
                    GRB_CHECK(p(i, j) == ((i == 0 || j == 0 || i == H + 1 || j == W + 1) ? '@' : g(i - 1, j - 1)));
            std::vector<int> b;
            for (int v = 0; v < H * W; ++v) {
                const bool e = (v / W == 0 || v / W == H - 1 || v % W == 0 || v % W == W - 1);
                GRB_CHECK(g.is_border(v) == e);
                if (e) b.push_back(v);
            }
            GRB_CHECK(g.border_cells() == b);
        }
        // 8. path_to_dirs / 方向の変換
        {
            std::vector<int> path{rnd(rng, 0, H * W - 1)};
            std::string ex;
            for (int k = 0; k < 10; ++k) {
                const int v = path.back(), d = rnd(rng, 0, 3);
                const int u = g.neighbor(v, d);
                if (u == -1) continue;
                path.push_back(u), ex += Grid::DIR_CHAR[d];
            }
            GRB_CHECK(g.path_to_dirs(path) == ex && g.path_to_dirs(path, true) == ex);
            std::vector<int> tpath{rnd(rng, 0, H * W - 1)};
            std::string tex;
            for (int k = 0; k < 10; ++k) {
                const int v = tpath.back(), d = rnd(rng, 0, 3);
                const int u = g.wrap_id(v / W + Grid::DI[d], v % W + Grid::DJ[d]);
                if ((d % 2 == 0 && H <= 2) || (d % 2 == 1 && W <= 2)) continue;  // 2 以下だと向きが決まらない
                tpath.push_back(u), tex += Grid::DIR_CHAR[d];
            }
            GRB_CHECK(g.path_to_dirs(tpath, true) == tex);
            const int a = rnd(rng, 0, H * W - 1), b2 = rnd(rng, 0, H * W - 1);
            const int di = b2 / W - a / W, dj = b2 % W - a % W;
            const char ec = (di < 0 && dj == 0) ? 'U' : (di == 0 && dj > 0) ? 'R' : (di > 0 && dj == 0) ? 'D' : (di == 0 && dj < 0) ? 'L' : '?';
            if (a != b2) GRB_CHECK(g.path_to_dirs({a, b2}) == std::string(1, ec));
            GRB_CHECK(sgn(0) == 0);
        }
    }
    // 方向の表
    for (int d = 0; d < 4; ++d) {
        GRB_CHECK(Grid::dir_of(Grid::DIR_CHAR[d]) == d && Grid::dir_of((char)(Grid::DIR_CHAR[d] - 'A' + 'a')) == d);
        GRB_CHECK(Grid::DI[Grid::reverse_dir(d)] == -Grid::DI[d] && Grid::DJ[Grid::reverse_dir(d)] == -Grid::DJ[d]);
        GRB_CHECK(Grid::turn_left(Grid::turn_right(d)) == d);
        // 時計回り: (di, dj) → (dj, -di)
        GRB_CHECK(Grid::DI[Grid::turn_right(d)] == Grid::DJ[d] && Grid::DJ[Grid::turn_right(d)] == -Grid::DI[d]);
    }
    for (int d = 4; d < 8; ++d) GRB_CHECK(Grid::DI[Grid::reverse_dir(d)] == -Grid::DI[d] && Grid::DJ[Grid::reverse_dir(d)] == -Grid::DJ[d]);
    GRB_CHECK(Grid::dir_of('^') == 0 && Grid::dir_of('>') == 1 && Grid::dir_of('v') == 2 && Grid::dir_of('<') == 3 && Grid::dir_of('#') == -1);
    // 入力・出力
    {
        CinFromString cf("3 4\n#..S\n.##.\nG...\n");
        const Grid g = read_grid();
        GRB_CHECK(g.H == 3 && g.W == 4 && g.s[0] == "#..S" && g.find('S') == 3 && g.find('G') == 8);
        std::ostringstream os;
        os << g;
        GRB_CHECK(os.str() == "#..S\n.##.\nG...");
        CinFromString cf2("1 2 3\n4 5 6\n");
        const auto a = read_matrix<long long>(2, 3);
        GRB_CHECK(a.size() == 2 && a[1][2] == 6 && a[0][0] == 1);
        const Grid e(2, 3);
        GRB_CHECK(e.count('.') == 6 && e.find('#') == -1);
        const Grid empty;
        GRB_CHECK(empty.size() == 0 && empty.find_all('.').empty());
    }
    std::cout << "[grid_base] random_test: OK (" << iters << " cases)\n";
}

static void benchmark() {
    if (TEST_SCALE == 1) std::cout << "[grid_base] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    for (int n : {1000, 2000}) {
        if (TEST_SCALE == 1 && n > 1000) break;
        Grid g(n, n);
        for (auto& r : g.s)
            for (char& c : r) c = (rng() % 4 == 0) ? '#' : '.';
        g.s[0][0] = '.';
        auto t0 = std::chrono::steady_clock::now();
        // for_each4 だけで書いた BFS（09_grid_bfs と同じ形）
        const std::vector<char> wall = g.mask('#');
        std::vector<int> dist(g.size(), -1), q;
        q.reserve(g.size());
        dist[0] = 0, q.push_back(0);
        for (size_t h = 0; h < q.size(); ++h) {
            const int v = q[h];
            g.for_each4(v, [&](int u) {
                if (!wall[u] && dist[u] < 0) dist[u] = dist[v] + 1, q.push_back(u);
            });
        }
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        t0 = std::chrono::steady_clock::now();
        long long sum = 0;
        for (int v = 0; v < g.size(); ++v) g.for_each8(v, [&](int u, int d) { sum += u ^ d; });
        const double ms8 = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "[grid_base] benchmark: " << n << "x" << n << "  for_each4 での BFS " << (long long)ms << " ms（到達 " << q.size()
                  << " マス）, 全マス for_each8 " << (long long)ms8 << " ms  (check " << sum % 1000 << ")\n";
    }
}

}  // namespace grb_test


// =============================================================================
//  使い方デモ
// =============================================================================
namespace grb_demo {

static std::string to_s(bool b) { return b ? "true" : "false"; }
template <class T, std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool> && !std::is_same_v<T, char>, int> = 0>
static std::string to_s(T v) { return std::to_string(v); }
static std::string to_s(char c) { return std::string("'") + c + "'"; }
static std::string to_s(const std::string& s) { return "\"" + s + "\""; }
static std::string to_s(const std::pair<int, int>& p) { return "(" + std::to_string(p.first) + ", " + std::to_string(p.second) + ")"; }
static std::string to_s(const std::vector<int>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + std::to_string(v[i]);
    return s + "}";
}
static std::string to_s(const std::vector<std::pair<int, int>>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
static std::string pad(const std::string& s, size_t w) {
    size_t cols = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = (unsigned char)s[i];
        if ((c & 0xC0) == 0x80) continue;
        cols += (c >= 0xE0 ? 2 : 1);
    }
    return cols >= w ? s : s + std::string(w - cols, ' ');
}
static void banner(const std::string& s) {
    std::cout << "\n=====================================================================\n"
              << "  " << s << "\n"
              << "=====================================================================\n";
}
static void api(const std::string& usage, const std::string& ret, const std::string& what, const std::string& cost) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << usage << "\n"
              << "    返り値 : " << ret << "\n"
              << "    内容   : " << what << "\n"
              << "    計算量 : " << cost << "\n";
}
static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
static void stmt(const std::string& s) {
    const bool block = !s.empty() && s.back() == '}' &&
                       (s.compare(0, 4, "for ") == 0 || s.compare(0, 3, "if ") == 0 || s.compare(0, 6, "while ") == 0);
    code(block ? s : s + ";");
}
template <class V>
static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 34) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 52) + "  // " + comment;
    std::cout << line << "\n";
}
#define SHOW(expr, comment) grb_demo::show(#expr, (expr), comment)
#define RUN(...)                   \
    grb_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__
static void print_grid(const Grid& g) {
    for (int i = 0; i < g.H; ++i) std::cout << "      " << g.s[i] << "\n";
}

static void cheat_sheet() {
    struct Row { const char* use; const char* ret; const char* what; };
    const Row rows[] = {
        {"Grid g = read_grid(H, W);", "", "H 行の文字列を読む（read_grid() は H W から）"},
        {"Grid g(H, W, '.') / Grid g(rows)", "", "大きさと文字で作る / vector<string> から作る"},
        {"g.H, g.W, g.s[i][j], g(i, j), g.at(id)", "", "大きさ・マスの文字（読み書きできる）"},
        {"g.in(i, j) / g.id(i, j) / g.pos(id)", "bool / int / pair", "範囲内か / 番号 / 座標"},
        {"g.find(c) / g.find_all(c) / g.count(c)", "int / vector / int", "文字（の集合）の位置・個数"},
        {"g.mask(c)", "vector<char>", "m[id] = 1（文字が当てはまる）"},
        {"g.for_each4(id, f) / for_each8", "", "隣のマスごとに f(u) か f(u, 方向)"},
        {"g.for_each4(i, j, f) / for_each8", "", "(i, j) 版。f(ni, nj) か f(ni, nj, 方向)"},
        {"g.neighbor(id, d)", "int", "方向 d の隣（範囲外は -1）"},
        {"g.for_each_offsets(id, offs, f)", "", "任意の移動リスト。f(u) か f(u, k)"},
        {"Grid::offsets_manhattan(D) など", "vector<pair>", "移動リストを作る（chebyshev / euclid_sq / knight）"},
        {"g.for_each4_torus(id, f) など", "", "トーラス版（8・offsets もある）"},
        {"g.slide(id, d, walls)", "int", "壁の手前まで滑ったマス"},
        {"g.wrap(i, j) / g.wrap_id(i, j)", "pair / int", "トーラスの座標の正規化"},
        {"g.padded('#')", "Grid", "外周を 1 マス囲んだコピー"},
        {"g.border_cells() / g.is_border(id)", "vector / bool", "外周のマス"},
        {"g.path_to_dirs(path)", "string", "経路 → \"URDL\" の文字列"},
        {"Grid::dir_of(c) / reverse_dir / turn_right / turn_left", "int", "方向の変換"},
        {"CharSet walls(\"#X\") / ~CharSet(\".\")", "", "文字の集合 / 補集合"},
        {"read_matrix<T>(H, W)", "vector<vector<T>>", "数値のグリッドを読む"},
    };
    std::cout << "\n■ API 一覧（マスは (i, j) か番号 id = i * W + j。方向 0 上 1 右 2 下 3 左 4~7 斜め）\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, 50) << pad(ret, 22) << r.what << "\n";
    }
}

static void demo() {
    banner("グリッド基盤 使い方デモ    見方:  コード  ->  返り値    // 意味");
    cheat_sheet();

    api("【1】 read_grid(H, W) / read_grid()", "Grid", "H 行の文字列を読む。read_grid() は先頭の「H W」も読む", "O(HW)");
    CinFromString cf("3 5\nS.#..\n..#.G\n.....\n");
    code("// 入力:  3 5 / S.#.. / ..#.G / .....");
    RUN(Grid g = read_grid());
    print_grid(g);
    SHOW(g.H, "行数");
    SHOW(g.W, "列数");
    SHOW(g(0, 2), "マス (0, 2) の文字（g.s[0][2] / g[0][2] でも同じ）");
    note("数値のグリッドは read_matrix<long long>(H, W) で読み、形だけ Grid g(H, W) で扱う（隣接の列挙などが使える）。");

    api("【2】 g.id(i, j) / g.pos(id) / g.in(i, j) / g.at(id) / g.size()", "", "マスの番号 id = i * W + j との変換。探索の距離配列は 1 次元にできる", "O(1)");
    SHOW(g.id(1, 4), "(1, 4) の番号 = 1 * 5 + 4");
    SHOW(g.pos(9), "番号 9 の座標（auto [i, j] = g.pos(v); で受け取る）");
    SHOW(g.in(3, 0), "範囲外");
    SHOW(g.at(9), "番号でマスの文字");
    SHOW(g.size(), "マスの数 H * W");

    api("【3】 g.find(c) / g.find_all(c) / g.count(c) / g.mask(c)", "int / vector<int> / int / vector<char>",
        "文字（または文字の集合）の最初の位置・全部の位置・個数・表。c は '#'、\"#X\"、CharSet どれでもよい", "O(HW)");
    SHOW(g.find('S'), "スタートの番号");
    SHOW(g.find('G'), "ゴールの番号");
    SHOW(g.find_all('#'), "壁のマス全部");
    SHOW(g.count("#G"), "'#' か 'G' のマスの数");
    SHOW(g.find('X'), "ない文字は -1");

    api("【4】 g.for_each4(id, f) / g.for_each8(id, f)", "なし（f を呼ぶ）",
        "隣のマス（範囲内のものだけ）ごとに f を呼ぶ。f(u) でも f(u, 方向) でもよい。順番は方向の番号順", "O(1)");
    RUN(std::vector<int> nb);
    RUN(g.for_each4(g.id(1, 1), [&](int u) { nb.push_back(u); }));
    SHOW(nb, "(1, 1) の上・右・下・左");
    RUN(std::vector<std::pair<int, int>> nd);
    RUN(g.for_each8(0, [&](int u, int d) { nd.push_back({u, d}); }));
    SHOW(nd, "(0, 0) の 8 近傍 {番号, 方向}（範囲外は呼ばれない）");
    note("BFS の中では: g.for_each4(v, [&](int u) { if (g.at(u) != '#' && dist[u] == -1) { ... } });");

    api("【5】 g.for_each4(i, j, f) / g.for_each8(i, j, f)", "なし", "座標版。f(ni, nj) か f(ni, nj, 方向)", "O(1)");
    RUN(std::vector<std::pair<int, int>> nij);
    RUN(g.for_each4(2, 4, [&](int ni, int nj) { nij.push_back({ni, nj}); }));
    SHOW(nij, "(2, 4) の隣（右と下は範囲外）");

    api("【6】 g.neighbor(id, d) / Grid::DI[d], Grid::DJ[d] / Grid::DIR_CHAR[d]", "int", "方向 d の隣（範囲外は -1）。方向の表", "O(1)");
    SHOW(g.neighbor(0, 2), "(0, 0) の下 = (1, 0)");
    SHOW(g.neighbor(0, 0), "(0, 0) の上は範囲外");
    SHOW(Grid::DIR_CHAR[3], "方向 3 の文字");

    api("【7】 Grid::dir_of(c) / reverse_dir(d) / turn_right(d) / turn_left(d)", "int",
        "'U','R','D','L'（小文字・矢印 ^ > v < も可）→ 方向番号、逆向き、右・左に 90 度回す", "O(1)");
    SHOW(Grid::dir_of('L'), "'L' は 3");
    SHOW(Grid::dir_of('v'), "矢印 'v' は下 = 2");
    SHOW(Grid::reverse_dir(1), "右の逆は左");
    SHOW(Grid::turn_right(0), "上を右に 90 度回すと右");
    SHOW(Grid::turn_left(0), "上を左に 90 度回すと左");

    api("【8】 g.for_each_offsets(id, offs, f) と Grid::offsets_xxx", "なし / vector<pair<int,int>>",
        "任意の移動のしかた（桂馬・距離 D 以内へのワープ・決まった距離へのジャンプ）で隣を列挙。f(u) か f(u, k)",
        "O(移動の種類数)");
    SHOW(Grid::offsets_knight(), "桂馬（ナイト）の 8 方向");
    SHOW(Grid::offsets_manhattan(1), "マンハッタン距離 1 以内（中心なし）= 4 方向");
    SHOW(Grid::offsets_chebyshev(1).size(), "チェビシェフ距離 1 以内 = 8 方向（ABC176 D は 2 以内）");
    SHOW(Grid::offsets_euclid_sq(5), "di^2 + dj^2 = 5 ちょうど（ABC272 D 型）");
    RUN(std::vector<int> kn);
    RUN(g.for_each_offsets(0, Grid::offsets_knight(), [&](int u) { kn.push_back(u); }));
    SHOW(kn, "(0, 0) から桂馬で行けるマス");
    note("offsets_manhattan(D, true) / offsets_chebyshev(D, true) で中心 (0, 0) も含められる。");

    api("【9】 g.slide(id, d, walls)", "int", "方向 d に、壁かグリッドの外の手前まで滑ったマス（氷の上を滑る問題: ABC311 D）", "O(滑った距離)");
    SHOW(g.slide(0, 1, '#'), "(0, 0) から右へ滑ると (0, 1) で止まる（(0, 2) が壁）");
    SHOW(g.slide(0, 2, '#'), "下へ滑ると (2, 0)");

    api("【10】 g.wrap(i, j) / g.wrap_id(i, j) / g.for_each4_torus(id, f)", "pair / int",
        "トーラス（上下・左右の端がつながる）の座標を 0 ~ H-1, 0 ~ W-1 に直す。負の数でもよい", "O(1)");
    SHOW(g.wrap(-1, 7), "(-1, 7) は (2, 2)");
    RUN(std::vector<int> tn);
    RUN(g.for_each4_torus(0, [&](int u) { tn.push_back(u); }));
    SHOW(tn, "(0, 0) のトーラスでの上・右・下・左");

    api("【11】 g.padded(c = '#')", "Grid", "外周を 1 マスずつ c で囲んだコピー。元の (i, j) は (i+1, j+1) になる", "O(HW)");
    RUN(Grid p = g.padded());
    print_grid(p);

    api("【12】 g.border_cells() / g.is_border(id)", "vector<int> / bool",
        "外周のマス。「外とつながっていない領域」は外周の '.' から多始点 BFS で塗って、残ったところ", "O(HW) / O(1)");
    SHOW(g.border_cells().size(), "3×5 の外周は 12 マス");
    SHOW(g.is_border(g.id(1, 2)), "(1, 2) は外周ではない");

    api("【13】 g.path_to_dirs(path, torus = false)", "string",
        "経路（マスの番号の列。r.path(t) の結果）を \"URDL\" の文字列にする。同じ行・列の何マス先でも向きの文字になる", "O(経路の長さ)");
    RUN(std::vector<int> path = {g.id(0, 0), g.id(1, 0), g.id(1, 1), g.id(2, 1), g.id(2, 4)});
    SHOW(g.path_to_dirs(path), "下・右・下・右（3 マス先も R）");

    api("【14】 CharSet", "", "文字の集合。壁などを O(1) で判定する。'#'、\"#X\"、string から暗黙に作れ、~ で補集合", "判定 O(1)");
    RUN(CharSet walls("#X"));
    SHOW(walls('#'), "'#' は含まれる");
    SHOW(walls('.'), "'.' は含まれない");
    RUN(CharSet not_dot = ~CharSet('.'));
    SHOW(not_dot('S'), "「'.' 以外」の集合に 'S' は含まれる");

    api("【15】 cout << g", "", "各行を改行区切りで出力（最後に改行なし）。答えのグリッドを出力するときに", "O(HW)");
    code("std::cout << g << '\\n';");
    std::cout << g << "\n";

    head("典型的な使い方");
    note("迷路の最短路     : 09_grid_bfs の grid_bfs(g, g.find('S'), \"#\")  →  r.dist[g.find('G')]");
    note("壁を壊す・ワープ : 09_grid_01bfs（隣接処理に for_each_offsets を使う）");
    note("向きを状態に持つ : StateId S(g.size(), 4); と for_each4(v, [&](int u, int d) {...}) の方向 d");
    note("経路の出力       : g.path_to_dirs(r.path(t))");
}

}  // namespace grb_demo


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
    grb_test::random_test();
    grb_test::benchmark();
    grb_demo::demo();
    return 0;
}
