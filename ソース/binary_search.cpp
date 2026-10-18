#pragma region template
#include <bits/stdc++.h>
#include <atcoder/all>
#include <random>
using namespace std;
using namespace atcoder;
using ll = long long;
using ull = unsigned long long;
#define GET_REP(_1,_2,_3,NAME,...) NAME
#define rep3(i,s,n) for(ll i=(s);i<(n);i++)      // [s,n)
#define rep2(i,n) rep3(i,0,n)
#define rep(...) GET_REP(__VA_ARGS__, rep3, rep2)(__VA_ARGS__)
#define erep3(i,s,n) for(ll i=(s);i<=(n);i++)    // [s,n]
#define erep2(i,n) erep3(i,0,n)
#define erep(...) GET_REP(__VA_ARGS__, erep3, erep2)(__VA_ARGS__)
#define rrep3(i,s,n) for(ll i=(s)-1;i>=(n);i--)  // s-1,s-2,...,n
#define rrep2(i,n) rrep3(i,n,0)
#define rrep(...) GET_REP(__VA_ARGS__, rrep3, rrep2)(__VA_ARGS__)
#define nall(a) a.begin(),a.end()
#define rall(a) a.rbegin(),a.rend()
#define YES cout<<"Yes\n"
#define NO cout<<"No\n"
#define YN {cout<<"Yes\n";}else{cout<<"No\n";}
#define BUFF cin.tie(nullptr);ios::sync_with_stdio(false)
#define next_p(v) next_permutation(v.begin(),v.end())
template<class T> using PQ = priority_queue<T, vector<T>>;//大きい順
template<class T> using PQ_g = priority_queue<T, vector<T>, greater<T>>;//小さい順
template<class T> int sizei(const T& x){return (int)x.size();}
template<class T>bool chmax(T& a,const T& b){return a<b?a=b,1:0;}
template<class T>bool chmin(T& a,const T& b){return a>b?a=b,1:0;}
template<class T>int sgn(T x){return (x>0)-(x<0);}
inline int popcnt(long long x){return __builtin_popcountll(x);}
const ll LINF = (1LL<<60);
const int INF = (1<<30);
template<class A, class B> istream& operator>>(istream& is, pair<A, B>& p) { is >> p.first >> p.second; return is; }
//template<class A, class B> ostream& operator<<(ostream& os, const pair<A, B>& p) { os << p.first << " " << p.second; return os; }
template<class T> istream& operator>>(istream& is, vector<T>& v) { for(auto& x : v) is >> x; return is; }
//template<class T> ostream& operator<<(ostream& os, const vector<T>& v) {for(int i = 0; i < (int)v.size(); i++) os << v[i] << (i + 1 == (int)v.size() ? "" : " ");return os;}
//cout<<FIXED(15)<<(double)1/7<<'\n';
#define FIXED(n) fixed << setprecision(n)
#pragma endregion
#pragma region GraphIO
struct GraphIO {
    // 非重みグラフ入力
    static vector<vector<int>>
    read_unweighted(
        int n,
        int m,
        vector<int>* indeg = nullptr,   // 必要なら渡す
        vector<int>* outdeg = nullptr,  // 必要なら渡す
        bool directed = false,
        int input_index = 1
    ){
        vector<vector<int>> adj(n);
        if(indeg)  indeg->assign(n, 0);
        if(outdeg) outdeg->assign(n, 0);
        for(int i = 0; i < m; i++){
            int u, v;
            cin >> u >> v;
            u -= input_index;
            v -= input_index;
            adj[u].push_back(v);
            if(outdeg) (*outdeg)[u]++;
            if(indeg)  (*indeg)[v]++;
            if(!directed){
                adj[v].push_back(u);
                if(outdeg) (*outdeg)[v]++;
                if(indeg)  (*indeg)[u]++;
            }
        }
        return adj; // RVOでコピーほぼなし
    }
    // 重み付きグラフ入力
    static vector<vector<pair<int,long long>>>
    read_weighted(
        int n,
        int m,
        vector<int>* indeg = nullptr,
        vector<int>* outdeg = nullptr,
        bool directed = false,
        int input_index = 1
    ){
        vector<vector<pair<int,long long>>> adj(n);
        if(indeg)  indeg->assign(n, 0);
        if(outdeg) outdeg->assign(n, 0);
        for(int i = 0; i < m; i++){
            int u, v;
            long long w;
            cin >> u >> v >> w;
            u -= input_index;
            v -= input_index;
            adj[u].push_back({v, w});
            if(outdeg) (*outdeg)[u]++;
            if(indeg)  (*indeg)[v]++;
            if(!directed){
                adj[v].push_back({u, w});
                if(outdeg) (*outdeg)[v]++;
                if(indeg)  (*indeg)[u]++;
            }
        }
        return adj;
    }
};
#pragma endregion
#pragma region debug
// ===== debug color =====
#define DBG_RESET "\033[0m"
#define DBG_CYAN "\033[36m"
#define DBG_COLOR DBG_CYAN
template<class T,class=void>struct is_iter:std::false_type{};
template<class T>struct is_iter<T,std::void_t<decltype(begin(std::declval<T>())),decltype(end(std::declval<T>()))>>:std::true_type{};
template<class T>constexpr bool is_iter_v=is_iter<T>::value&&!std::is_same_v<T,std::string>;
template<class A,class B>
ostream&operator<<(ostream&os,const pair<A,B>&p){return os<<p.first<<","<<p.second;}
template<size_t I,class...Ts>
void _tuple_out(ostream&os,const tuple<Ts...>&t){
    if constexpr(I<sizeof...(Ts)){
        if(I) os<<",";
        os<<get<I>(t);
        _tuple_out<I+1>(os,t);
    }
}
template<class...Ts>
ostream&operator<<(ostream&os,const tuple<Ts...>&t){
    _tuple_out<0>(os,t);
    return os;
}
template<class T>
std::enable_if_t<is_iter_v<T>,ostream&>
operator<<(ostream&os,const T&v){
    using V = std::decay_t<decltype(*begin(v))>;
    if constexpr(is_iter_v<V>){
        bool first_row=true;
        os<<"\n";
        for(auto&&row:v){
            if(!first_row) os<<"\n";
            first_row=false;
            bool f=0;
            for(auto&&x:row){
                if(f) os<<" ";
                os<<x;
                f=1;
            }
        }
    }
    else{
        os<<"[";
        bool f=0;
        for(auto&&x:v){
            if(f) os<<", ";
            os<<x;
            f=1;
        }
        os<<"]";
    }
    return os;
}
void _dbg(const char*){
    cerr<<DBG_COLOR<<"========================================\n"<<DBG_RESET;
}
template<class T,class...U>
void _dbg(const char*names,T&&x,U&&...y){
    const char*comma=strchr(names,',');
    if(!comma){
        cerr<<DBG_COLOR<<names<<" = "<<x<<"\n"<<DBG_RESET;
    }else{
        cerr<<DBG_COLOR;
        cerr.write(names,comma-names)<<" = "<<x<<" | ";
        _dbg(comma+1,y...);
    }
}
// run.shを引数ありで実行するとdebug有効、なしだと無効
#ifdef LOCAL
#define debug(...) _dbg(#__VA_ARGS__ , ##__VA_ARGS__)
#else
#define debug(...)
#endif
#pragma endregion

#pragma region BinarySearch
// =============================================================================
//  BinarySearch : vector<long long> 用「全部入り」二分探索テンプレート（AtCoder 用）
// -----------------------------------------------------------------------------
//  旧 3 本（BinarySearch.txt / BinarySearchTemplateTest.txt / SortedArrayHelperTemplateTest.txt）
//  を 1 本にまとめたもの。実戦では「#pragma region BinarySearch ～ endregion」だけをコピーする。
//
//  ■ 使い方は 2 通り
//   (1) 配列の二分探索 : BinarySearch bs(A, false);   // 第 2 引数 = A がソート済みなら true
//                        bs.count_le(x)  bs.min_ge(x)  bs.sum_lt(x)  bs.kth(k) ...
//   (2) 答えで二分探索 : first_true(lo, hi, f)  last_true  meguru  first_true_real ...（自由関数）
//
//  ■ 共通の約束
//   ・index は 0-indexed。index の区間は半開区間 [l, r)
//   ・値の区間 [L, R] は閉区間（L <= a_i <= R）。L > R なら空。  ※ L <= a_i < R は [L, R - 1]
//   ・「該当なし」の返り値 : index を返す関数は n か -1（各関数に記載）、
//                            値を返す関数は std::nullopt（→ .value_or(既定値) で好きな値に）
//   ・和（sum_*）は mod 2^64 で計算する。「答えが long long に収まれば」途中であふれても正しい
//   ・第 2 引数 true（ソート済み）で渡したときは assert で昇順か確かめる（O(n)）
//   ・GCC / Clang 用（__int128 と __builtin_ctz を使う）。C++17 以上
//
//  ■ 高速化
//   ・探索の本体は静的 B+ 木（S+ tree）。8 要素 = 64 バイト（キャッシュ 1 行）ずつのノードに分けた
//     9 分木を上から下る。1 段で 8 個を同時に比較するので、段数は普通の二分探索の約 1/3
//       AVX-512 が使えるビルド : 1 段 = 比較 1 命令     （-march=native など）
//       AVX2 が使えるビルド    : 1 段 = 比較 2 命令     （-mavx2 など）
//       どちらも無いビルド     : 1 段 = 分岐なしの比較 8 回
//   ・upper_bound(x) は lower_bound(x + 1) で求める（整数なので同じ）→ 探索ルーチンは 1 本だけ
//   ・lower_bound_many : 多数のクエリをまとめて同時に下り、メモリ待ちを重ねて隠す
//   ・ソート : n >= 2048 は 11bit × 最大 6 パスの基数ソート（全要素で同じ桁のパスは飛ばす）。
//             元からソート済みならソートしない
//   ・累積和・重複除去の配列は、初めて使ったときに作る（使わなければ時間もメモリも 0）
//   ・答えの二分探索は中点を unsigned で計算 → long long 全域 [-2^63, 2^63) であふれない
//   ・実数の二分探索は double のビット列を整数とみなして二分探索する
//     → 最大 64 回の評価で「隣り合う 2 つの double」まで厳密に求まる（回数も誤差も気にしなくてよい）
//
//  ■ 計算量・メモリ（n = 要素数）
//     構築 O(n log n)（ソート済みなら O(n)）。メモリは約 9n バイト（vector<long long> の約 1.13 倍）
//     和を使うと +8n、相異なる値を使うと +8×(種類数) バイト。1 回の探索 O(log n)
//
//  ■ 計測（n = 2*10^5 のランダム配列。lower_bound 1 回あたり。GCC 13 -O2、作成時のマシン）
//       std::lower_bound 約 140 ns → BinarySearch  SIMD なし 約 50 ns / AVX2 約 25 ns / AVX-512 約 16 ns
//       lower_bound_many（まとめて）約 13～33 ns。構築（ソート込み）は std::sort 単体の約 1/3 の時間
//       kth_pair_prod（n = 2*10^5、ABC155 D 型）約 55 ms
//     自分の環境では bs_test::benchmark() で測れる
//
//  ■ 旧テンプレートからの置き換え
//     sorted(A) + 各関数           → BinarySearch bs(A, false);
//     count_less / count_lt(x)      → bs.count_lt(x)      count_leq / count_le → bs.count_le(x)
//     count_geq / count_ge(x)       → bs.count_ge(x)      count_greater / count_gt → bs.count_gt(x)
//     count_range(A, L, R)          → bs.count_range(L, R)（閉区間のまま）
//     range_count(A, L, R)          → bs.range_index(L, R)  ※返り値は半開区間 [l, r)（個数 = r - l）
//     exists(A, x)                  → bs.exists(x)
//     closest(A, x)（差）           → bs.min_dist(x)      値がほしいときは bs.nearest(x)
//     max_right_leq(A, i, x)        → j = bs.upper_bound(x) - 1;  j > i ? j : -1
//     max_len_diff_leq(A, K)        → bs.max_window_count(K)
//     count_pair_leq(A, K)          → bs.count_pairs_sum_le(K)
//     kth_value(A, k)（1-indexed）  → bs.kth(k - 1)（0-indexed）
//     BinarySearchArray の first_true(l, r, f) / last_true / exists（f は index を受け取る）
//                                   → 自由関数 first_true(l, r, f) / last_true(l, r, f) / first_true(l, r, f) < r
//
// -----------------------------------------------------------------------------
//  【API 一覧】                                   x, L, R, K : long long   i, l, r, k : int
// -----------------------------------------------------------------------------
//  ◆ 構築
//    BinarySearch bs(A, is_sorted)   A をコピー（is_sorted = false ならソート）  O(n log n) / O(n)
//    BinarySearch bs(std::move(A), is_sorted)   コピーせずに受け取る
//    bs.build(A, is_sorted)          作り直す                                  同上
//  ◆ 基本                                                                        すべて O(1)
//    size() empty() bs[i] kth(k) kth_largest(k) min() max()   k は 0-indexed
//    begin() end()（for (long long v : bs) が使える） data()（const long long*）
//    to_vector()（ソート済みの vector を作って返す O(n)）
//  ◆ 位置（index）                                                               O(log n)
//    lower_bound(x)    a_i >= x となる最小の i（無ければ n）  = count_lt(x)
//    upper_bound(x)    a_i >  x となる最小の i（無ければ n）  = count_le(x)
//    equal_range(x)    x の範囲 [lower_bound(x), upper_bound(x))
//    range_index(L, R) L <= a_i <= R となる i の範囲 [l, r)（L > R なら空）
//    first_index(x) / last_index(x)   x が最初 / 最後に現れる i（無ければ -1）
//  ◆ 個数                                                                        O(log n)
//    count_lt(x) count_le(x) count_gt(x) count_ge(x) count(x)（== x）
//    count_range(L, R)（L <= a_i <= R）  exists(x)  exists_range(L, R)
//  ◆ 前後の値（std::optional<long long>。無ければ nullopt）                     O(log n)
//    max_lt(x)（x 未満の最大） max_le(x)（x 以下の最大）
//    min_gt(x)（x より大きい最小） min_ge(x)（x 以上の最小）
//    nearest(x)（x に最も近い値。同じ距離なら小さい方）
//    min_dist(x) → long long  min |a_i - x|（空なら LLONG_MAX）
//  ◆ 和（初回だけ累積和を作る O(n)）                                             O(log n)
//    prefix_sum(i)（a_0 + … + a_{i-1}） sum_index(l, r)（a_l + … + a_{r-1}） sum_all()
//    sum_lt(x) sum_le(x) sum_gt(x) sum_ge(x) sum_range(L, R)（その条件を満たす a_i の和）
//    sum_smallest(k) / sum_largest(k)（小さい / 大きい方から k 個の和）
//    sum_abs_diff(x)（Σ|a_i - x|）
//  ◆ 相異なる値（初回だけ重複除去 O(n)）                                        O(log n)
//    distinct_count() distinct_count_range(L, R)  distinct_values()（const vector&）
//    distinct_rank(x)（x 未満の相異なる値の個数 = 座標圧縮後の番号） distinct_kth(k)
//  ◆ 条件で探す（f は値 a_i を受け取る。a の昇順に単調であること）      O(log n) 回 f を呼ぶ
//    first_index_if(f)       f(a_i) が true の最小の i（F…F T…T。無ければ n）
//    last_index_if(f)        f(a_i) が true の最大の i（T…T F…F。無ければ -1）
//    first_index_if(l, r, f) / last_index_if(l, r, f)   [l, r) 版（無ければ r / l - 1）
//  ◆ まとめて探索                                                    O(q log n)（1 回ずつより速い）
//    lower_bound_many(xs) / upper_bound_many(xs)   → vector<int>
//  ◆ ペア（i < j の組を数える。和・差・積はあふれないように比較するので、値は long long 全域でよい）
//    count_pairs_sum_le/lt/ge/gt(K)    a_i + a_j ⋈ K の組の数             O(n)
//    count_pairs_diff_le/lt/ge/gt(K)   |a_i - a_j| ⋈ K の組の数           O(n)
//    count_pairs_prod_le/lt/ge/gt(K)   a_i * a_j ⋈ K の組の数（負・0 も可）O(n)
//    max_window_count(W)   幅 W の値の区間 [v, v + W] に入る要素数の最大  O(n)
//    kth_pair_sum(k) / kth_pair_diff(k) / kth_pair_prod(k)
//                          その値を小さい順に並べたときの k 番目（0-indexed）
//                          それぞれ O(n log V)（V = 値の幅）   ※答えが long long に収まること
//  ◆ 2 つの配列（i は bs、j は b から 1 個ずつ。n × m 組）
//    count_sum_le/lt/ge/gt(b, K)   a_i + b_j ⋈ K の組の数                 O(n + m)
//    kth_sum(b, k)                 a_i + b_j の k 番目に小さい値（0-indexed）  O((n + m) log V)
//
//  ◆ 自由関数 : 答えで二分探索（f は long long を受け取り bool を返す）
//    first_true(lo, hi, f)   [lo, hi) で F…F T…T のとき、true の最小（無ければ hi）
//    last_true(lo, hi, f)    [lo, hi) で T…T F…F のとき、true の最大（無ければ lo - 1）
//    meguru(ok, ng, f)       めぐる式。f(ok) = true, f(ng) = false とみなし（呼ばない）境界の ok を返す
//                            ok < ng でも ok > ng でもよい
//    first_true_gallop(lo, hi, f) / last_true_gallop(lo, hi, f)
//                            結果は first_true / last_true と同じ。lo から 1, 2, 4, … と幅を広げて探す
//                            → 評価回数 O(log(答え - lo))。f を lo + 2*(答え - lo) + 2 までしか呼ばない
//                              （大きい値で f があふれる、hi が分からない、答えが lo に近いときに）
//    評価回数 : 最大 ceil(log2(hi - lo + 1)) 回（long long 全域でも 64 回）
//  ◆ 自由関数 : 実数（double）で二分探索（f は double を受け取る）              最大 64 回
//    first_true_real(lo, hi, f)   F…F T…T。f(lo) = false, f(hi) = true とみなす。true の最小の x
//    last_true_real(lo, hi, f)    T…T F…F。f(lo) = true, f(hi) = false とみなす。true の最大の x
//    meguru_real(ok, ng, f)       めぐる式（ok < ng でも ok > ng でもよい）
//    ※ lo, hi に ±INFINITY を渡してもよい
//  ◆ 自由関数 : 凸関数の最小・最大（整数は傾きの符号で二分探索）
//    argmin_convex(lo, hi, f)     [lo, hi) で下に凸な f の最小値をとる最小の x   2 log(hi - lo) 回
//    argmax_concave(lo, hi, f)    [lo, hi) で上に凸な f の最大値をとる最小の x
//    argmin_convex_real(lo, hi, f, iter = 100) / argmax_concave_real(…)   [lo, hi] の黄金分割探索
//  ◆ 自由関数 : 並列二分探索
//    parallel_binary_search(Q, T, reset, apply, check) → vector<int>
//      イベント 0..T-1 を順に適用していくとき、各クエリ q について
//      「最初の t 個を適用した状態で check(q) が true」となる最小の t（0..T）を返す。無ければ T + 1
//      reset() : 状態を「何も適用していない」に戻す / apply(t) : イベント t を適用 / check(q) : bool
//      計算量 O((T + Q) log T) 回の apply / check
// =============================================================================
#if defined(__AVX512F__) || defined(__AVX2__)
#include <immintrin.h>
#endif

namespace bs_internal {
using i128 = __int128;
using u64 = unsigned long long;

// 64 バイト境界に揃えて確保するアロケータ（B+ 木の 1 ノード = キャッシュ 1 行にするため）
template <class T>
struct AlignedAlloc {
    using value_type = T;
    AlignedAlloc() = default;
    template <class U>
    AlignedAlloc(const AlignedAlloc<U>&) {}
    T* allocate(std::size_t n) { return static_cast<T*>(::operator new(n * sizeof(T), std::align_val_t(64))); }
    void deallocate(T* p, std::size_t) { ::operator delete(p, std::align_val_t(64)); }
    template <class U>
    bool operator==(const AlignedAlloc<U>&) const { return true; }
    template <class U>
    bool operator!=(const AlignedAlloc<U>&) const { return false; }
};

// 昇順ソート。n < 2048 は std::sort、それ以上は 11bit × 6 パスの LSD 基数ソート。
// 符号ビットを反転すると「unsigned の大小 = signed の大小」になる。全要素で同じ桁のパスは飛ばす
// （値が 0 ～ 10^9 程度なら 3 パスで終わる）
inline void sort_ll(std::vector<long long>& v) {
    const std::size_t n = v.size();
    if (n < 2048) {
        std::sort(v.begin(), v.end());
        return;
    }
    constexpr int R = 11, K = 1 << R, P = 6;
    constexpr u64 FLIP = 1ULL << 63;
    std::vector<unsigned> cnt((std::size_t)P * K, 0);
    std::vector<u64> buf(n);
    u64* src = reinterpret_cast<u64*>(v.data());
    u64* dst = buf.data();
    for (std::size_t i = 0; i < n; ++i) {
        const u64 x = src[i] ^ FLIP;
        for (int d = 0; d < P; ++d) ++cnt[(std::size_t)d * K + ((x >> (R * d)) & (K - 1))];
    }
    bool in_buf = false;
    for (int d = 0; d < P; ++d) {
        unsigned* c = cnt.data() + (std::size_t)d * K;
        const int sh = R * d;
        if (c[((src[0] ^ FLIP) >> sh) & (K - 1)] == n) continue;  // この桁は全要素同じ
        unsigned s = 0;
        for (int b = 0; b < K; ++b) {
            const unsigned t = c[b];
            c[b] = s;
            s += t;
        }
        for (std::size_t i = 0; i < n; ++i) {
            const u64 x = src[i];
            dst[c[((x ^ FLIP) >> sh) & (K - 1)]++] = x;
        }
        std::swap(src, dst);
        in_buf = !in_buf;
    }
    if (in_buf) std::memcpy(v.data(), src, n * sizeof(u64));
}

// double ↔ 順序を保つ 64bit 整数（-inf < … < -0.0 < +0.0 < … < +inf が連続した整数になる）
inline u64 dkey(double d) {
    u64 u;
    std::memcpy(&u, &d, sizeof(u));
    return (u >> 63) ? ~u : (u | (1ULL << 63));
}
inline double dval(u64 k) {
    const u64 u = (k >> 63) ? (k & ~(1ULL << 63)) : ~k;
    double d;
    std::memcpy(&d, &u, sizeof(d));
    return d;
}
}  // namespace bs_internal

// =============================================================================
//  答えで二分探索（整数）
// =============================================================================

// [lo, hi) で f が F…F T…T のとき、f(x) = true となる最小の x。無ければ hi
// long long 全域で使える（中点を unsigned で計算するのであふれない）。f は [lo, hi) の中でしか呼ばない
template <class F>
long long first_true(long long lo, long long hi, F f) {
    using u64 = unsigned long long;
    if (lo >= hi) return hi;
    u64 len = (u64)hi - (u64)lo;  // 答えは [base, base + len]（base + len = hi は「無し」）
    long long base = lo;
    while (len > 0) {
        const u64 half = len >> 1;
        const long long mid = (long long)((u64)base + half);
        if (f(mid)) {
            len = half;
        } else {
            base = mid + 1;
            len -= half + 1;
        }
    }
    return base;
}

// [lo, hi) で f が T…T F…F のとき、f(x) = true となる最大の x。無ければ lo - 1（lo > LLONG_MIN のこと）
template <class F>
long long last_true(long long lo, long long hi, F f) {
    return first_true(lo, hi, [&](long long x) { return !f(x); }) - 1;
}

// めぐる式二分探索。f(ok) = true, f(ng) = false とみなし（この 2 点では f を呼ばない）、
// |ok - ng| = 1 になるまで縮めたときの ok を返す。ok < ng でも ok > ng でもよい
template <class F>
long long meguru(long long ok, long long ng, F f) {
    if (ok < ng) return last_true(ok + 1, ng, f);  // (ok, ng) は T…T F…F。無ければ ok
    if (ok > ng) return first_true(ng + 1, ok, f);  // (ng, ok) は F…F T…T。無ければ ok
    return ok;
}

// first_true と同じ結果を、lo から幅 1, 2, 4, 8, … と広げながら探す（ギャロップ探索）。
// 評価回数 O(log(答え - lo))、f を呼ぶのは lo + 2*(答え - lo) + 2 以下の x だけ
template <class F>
long long first_true_gallop(long long lo, long long hi, F f) {
    using u64 = unsigned long long;
    if (lo >= hi) return hi;
    const u64 len = (u64)hi - (u64)lo;
    u64 l = 0, step = 1;  // オフセット [0, l) は false と確定
    while (true) {
        if (step == 0 || step - 1 >= len - l) return first_true((long long)((u64)lo + l), hi, f);
        const u64 p = l + step - 1;
        const long long x = (long long)((u64)lo + p);
        if (f(x)) return first_true((long long)((u64)lo + l), x, f);  // 答えは [l, p]
        l = p + 1;
        step <<= 1;
    }
}

// last_true と同じ結果を、lo 側からギャロップで探す（T の区間が短いときに速い）。無ければ lo - 1
template <class F>
long long last_true_gallop(long long lo, long long hi, F f) {
    return first_true_gallop(lo, hi, [&](long long x) { return !f(x); }) - 1;
}

// =============================================================================
//  答えで二分探索（実数）: double のビット列で二分探索 → 最大 64 回で隣り合う double まで厳密
// =============================================================================

// めぐる式。f(ok) = true, f(ng) = false とみなし、ok と ng が隣り合う double になったときの ok を返す
template <class F>
double meguru_real(double ok, double ng, F f) {
    using u64 = unsigned long long;
    u64 a = bs_internal::dkey(ok), b = bs_internal::dkey(ng);
    while (true) {
        const u64 d = a < b ? b - a : a - b;
        if (d <= 1) break;
        const u64 m = a < b ? a + (d >> 1) : a - (d >> 1);
        if (f(bs_internal::dval(m))) a = m;
        else b = m;
    }
    return bs_internal::dval(a);
}
// F…F T…T。f(lo) = false, f(hi) = true とみなし、f(x) = true となる最小の x
template <class F>
double first_true_real(double lo, double hi, F f) {
    return meguru_real(hi, lo, f);
}
// T…T F…F。f(lo) = true, f(hi) = false とみなし、f(x) = true となる最大の x
template <class F>
double last_true_real(double lo, double hi, F f) {
    return meguru_real(lo, hi, f);
}

// =============================================================================
//  凸関数の最小・最大
// =============================================================================

// [lo, hi)（lo < hi）で下に凸な f の最小値をとる x のうち最小のもの。
// 傾き f(x + 1) - f(x) >= 0 となる最初の x を二分探索する（平らな部分があってもよい）
template <class F>
long long argmin_convex(long long lo, long long hi, F f) {
    assert(lo < hi);
    return first_true(lo, hi - 1, [&](long long x) { return !(f(x + 1) < f(x)); });
}
// [lo, hi)（lo < hi）で上に凸な f の最大値をとる x のうち最小のもの
template <class F>
long long argmax_concave(long long lo, long long hi, F f) {
    assert(lo < hi);
    return first_true(lo, hi - 1, [&](long long x) { return !(f(x) < f(x + 1)); });
}
// [lo, hi] で下に凸（狭義の単峰）な f の最小点。黄金分割探索（1 回ごとに f を 1 回呼ぶ）
// iter = 100 で区間は 0.618^100 ≈ 1e-21 倍になる
template <class F>
double argmin_convex_real(double lo, double hi, F f, int iter = 100) {
    const double r = (3.0 - std::sqrt(5.0)) / 2.0;  // 0.3819…
    double a = lo, b = hi;
    double c = a + r * (b - a), d = b - r * (b - a);
    auto fc = f(c), fd = f(d);
    for (int i = 0; i < iter; ++i) {
        if (fc < fd) {
            b = d, d = c, fd = fc;
            c = a + r * (b - a), fc = f(c);
        } else {
            a = c, c = d, fc = fd;
            d = b - r * (b - a), fd = f(d);
        }
    }
    return (a + b) / 2;
}
// [lo, hi] で上に凸（狭義の単峰）な f の最大点
template <class F>
double argmax_concave_real(double lo, double hi, F f, int iter = 100) {
    return argmin_convex_real(lo, hi, [&](double x) { return -f(x); }, iter);
}

// =============================================================================
//  並列二分探索
// =============================================================================
// イベント 0..T-1 を順に適用していくとき、各クエリ q (0..Q-1) について
// 「最初の t 個のイベントを適用した状態で check(q) が true」となる最小の t (0..T) を返す。
// T 個すべて適用しても false なら T + 1。check(q) は t について単調（一度 true になったらずっと true）
//   reset()  : 状態を「何も適用していない」に戻す
//   apply(t) : イベント t を適用する
//   check(q) : 今の状態でクエリ q が条件を満たすか
// 計算量 : ceil(log2(T + 2)) 周 × (apply 最大 T 回 + check Q 回)
template <class Reset, class Apply, class Check>
std::vector<int> parallel_binary_search(int Q, int T, Reset reset, Apply apply, Check check) {
    std::vector<int> ok(Q, T + 1), ng(Q, -1), start(T + 2), order(Q);
    while (true) {
        std::fill(start.begin(), start.end(), 0);
        int last = -1;
        for (int q = 0; q < Q; ++q)
            if (ok[q] - ng[q] > 1) {
                const int mid = (ok[q] + ng[q]) >> 1;
                ++start[mid + 1];
                last = std::max(last, mid);
            }
        if (last < 0) break;
        for (int t = 0; t <= T; ++t) start[t + 1] += start[t];
        {
            std::vector<int> pos(start.begin(), start.end() - 1);
            for (int q = 0; q < Q; ++q)
                if (ok[q] - ng[q] > 1) order[pos[(ok[q] + ng[q]) >> 1]++] = q;
        }
        reset();
        for (int t = 0; t <= last; ++t) {
            for (int i = start[t]; i < start[t + 1]; ++i) {
                const int q = order[i];
                if (check(q)) ok[q] = t;
                else ng[q] = t;
            }
            if (t < last) apply(t);
        }
    }
    return ok;
}

// =============================================================================
//  BinarySearch 本体
// =============================================================================
struct BinarySearch {
    using i128 = bs_internal::i128;
    using u64 = bs_internal::u64;

    // ------------------------------------------------------------------ 構築
    BinarySearch() { build(std::vector<long long>(), true); }
    BinarySearch(std::vector<long long> v, bool is_sorted) { build(std::move(v), is_sorted); }

    // is_sorted = true  : v はソート済み（assert で確認するだけ。O(n)）
    // is_sorted = false : ソートする（元からソート済みならしない）
    void build(std::vector<long long> v, bool is_sorted) {
        if (is_sorted) {
            assert(std::is_sorted(v.begin(), v.end()) && "BinarySearch: is_sorted = true なのに昇順でない");
        } else if (!std::is_sorted(v.begin(), v.end())) {
            bs_internal::sort_ll(v);
        }
        assert(v.size() < (std::size_t)INT_MAX - 64);
        n_ = (int)v.size();
        // 各段のノード数。葉は n/8 + 1 個にして、a[n] に必ず番兵 LLONG_MAX が入るようにする
        std::size_t cnt[MAXH];
        int H = 0;
        cnt[H++] = (std::size_t)n_ / B + 1;
        while (cnt[H - 1] > 1) {
            cnt[H] = (cnt[H - 1] + B) / (B + 1);
            ++H;
        }
        H_ = H;
        std::size_t tot = 0;  // 根を先頭に、葉を最後に並べる
        for (int h = H - 1; h >= 0; --h) {
            off_[h] = tot * B;
            tot += cnt[h];
        }
        tr_.assign(tot * B, LLONG_MAX);
        long long* leaf = tr_.data() + off_[0];
        std::copy(v.begin(), v.end(), leaf);
        // 段 h のノード j のキー c = 子 (j*9 + c + 1) の部分木の最小値（= その部分木の最初の葉の先頭）
        std::size_t span = B;  // 段 h-1 の 1 ノードの下にある葉の要素数 = 8 * 9^(h-1)
        for (int h = 1; h < H; ++h) {
            long long* node = tr_.data() + off_[h];
            for (std::size_t j = 0; j < cnt[h]; ++j)
                for (int c = 0; c < B; ++c) {
                    const std::size_t child = j * (B + 1) + c + 1;
                    if (child >= cnt[h - 1]) break;
                    node[j * B + c] = leaf[child * span];
                }
            span *= B + 1;
        }
        pre_.clear();
        pre_.shrink_to_fit();
        uniq_.clear();
        uniq_.shrink_to_fit();
        uniq_built_ = false;
    }

    // ------------------------------------------------------------------ 基本
    int size() const { return n_; }
    bool empty() const { return n_ == 0; }
    long long operator[](int i) const { return A()[i]; }
    long long kth(int k) const {  // k 番目に小さい値（0-indexed）
        assert(0 <= k && k < n_);
        return A()[k];
    }
    long long kth_largest(int k) const {  // k 番目に大きい値（0-indexed）
        assert(0 <= k && k < n_);
        return A()[n_ - 1 - k];
    }
    long long min() const {
        assert(n_ > 0);
        return A()[0];
    }
    long long max() const {
        assert(n_ > 0);
        return A()[n_ - 1];
    }
    const long long* data() const { return A(); }
    const long long* begin() const { return A(); }
    const long long* end() const { return A() + n_; }
    std::vector<long long> to_vector() const { return std::vector<long long>(begin(), end()); }

    // ------------------------------------------------------------------ 位置（index）
    // a_i >= x となる最小の i（無ければ n）
    int lower_bound(long long x) const {
        const long long* t = tr_.data();
        std::size_t k = 0;
        for (int h = H_ - 1; h > 0; --h) k = k * (B + 1) + node_count(t + off_[h] + k * B, x);
        return (int)(k * B + node_count(t + off_[0] + k * B, x));
    }
    // a_i > x となる最小の i（無ければ n）
    int upper_bound(long long x) const { return x == LLONG_MAX ? n_ : lower_bound(x + 1); }
    std::pair<int, int> equal_range(long long x) const { return {lower_bound(x), upper_bound(x)}; }
    // L <= a_i <= R となる i の範囲 [l, r)。L > R なら空（l == r）
    std::pair<int, int> range_index(long long L, long long R) const {
        const int l = lower_bound(L);
        return {l, L > R ? l : upper_bound(R)};
    }
    int first_index(long long x) const {  // x が最初に現れる i（無ければ -1）
        const int i = lower_bound(x);
        return (i < n_ && A()[i] == x) ? i : -1;
    }
    int last_index(long long x) const {  // x が最後に現れる i（無ければ -1）
        const int i = upper_bound(x) - 1;
        return (i >= 0 && A()[i] == x) ? i : -1;
    }

    // ------------------------------------------------------------------ 個数
    int count_lt(long long x) const { return lower_bound(x); }
    int count_le(long long x) const { return upper_bound(x); }
    int count_gt(long long x) const { return n_ - upper_bound(x); }
    int count_ge(long long x) const { return n_ - lower_bound(x); }
    int count(long long x) const { return upper_bound(x) - lower_bound(x); }
    int count_range(long long L, long long R) const { return L > R ? 0 : upper_bound(R) - lower_bound(L); }
    bool exists(long long x) const {
        const int i = lower_bound(x);
        return i < n_ && A()[i] == x;
    }
    bool exists_range(long long L, long long R) const {
        if (L > R) return false;
        const int i = lower_bound(L);
        return i < n_ && A()[i] <= R;
    }

    // ------------------------------------------------------------------ 前後の値
    std::optional<long long> max_lt(long long x) const {
        const int i = lower_bound(x);
        if (i == 0) return std::nullopt;
        return A()[i - 1];
    }
    std::optional<long long> max_le(long long x) const {
        const int i = upper_bound(x);
        if (i == 0) return std::nullopt;
        return A()[i - 1];
    }
    std::optional<long long> min_gt(long long x) const {
        const int i = upper_bound(x);
        if (i == n_) return std::nullopt;
        return A()[i];
    }
    std::optional<long long> min_ge(long long x) const {
        const int i = lower_bound(x);
        if (i == n_) return std::nullopt;
        return A()[i];
    }
    // x に最も近い値（距離が同じなら小さい方）
    std::optional<long long> nearest(long long x) const {
        if (n_ == 0) return std::nullopt;
        const int i = lower_bound(x);
        if (i == n_) return A()[n_ - 1];
        if (i == 0) return A()[0];
        const u64 dl = (u64)x - (u64)A()[i - 1], dr = (u64)A()[i] - (u64)x;
        return dl <= dr ? A()[i - 1] : A()[i];
    }
    // min |a_i - x|。空なら LLONG_MAX（差が 2^63 以上のときも LLONG_MAX）
    long long min_dist(long long x) const {
        if (n_ == 0) return LLONG_MAX;
        const int i = lower_bound(x);
        u64 d = ~0ULL;
        if (i < n_) d = (u64)A()[i] - (u64)x;
        if (i > 0) d = std::min(d, (u64)x - (u64)A()[i - 1]);
        return d > (u64)LLONG_MAX ? LLONG_MAX : (long long)d;
    }

    // ------------------------------------------------------------------ 和（mod 2^64 で計算）
    long long prefix_sum(int i) const {  // a_0 + … + a_{i-1}（0 <= i <= n）
        ensure_pre();
        return (long long)pre_[i];
    }
    long long sum_index(int l, int r) const {  // a_l + … + a_{r-1}（0 <= l <= r <= n）
        ensure_pre();
        return (long long)(pre_[r] - pre_[l]);
    }
    long long sum_all() const {
        ensure_pre();
        return (long long)pre_[n_];
    }
    long long sum_lt(long long x) const { return sum_index(0, lower_bound(x)); }
    long long sum_le(long long x) const { return sum_index(0, upper_bound(x)); }
    long long sum_gt(long long x) const { return sum_index(upper_bound(x), n_); }
    long long sum_ge(long long x) const { return sum_index(lower_bound(x), n_); }
    long long sum_range(long long L, long long R) const {  // L <= a_i <= R の a_i の和
        if (L > R) return 0;
        return sum_index(lower_bound(L), upper_bound(R));
    }
    long long sum_smallest(int k) const {  // 小さい方から k 個の和（0 <= k <= n）
        assert(0 <= k && k <= n_);
        return sum_index(0, k);
    }
    long long sum_largest(int k) const {  // 大きい方から k 個の和（0 <= k <= n）
        assert(0 <= k && k <= n_);
        return sum_index(n_ - k, n_);
    }
    long long sum_abs_diff(long long x) const {  // Σ|a_i - x|
        ensure_pre();
        const int i = lower_bound(x);
        const u64 X = (u64)x;
        const u64 left = X * (u64)i - pre_[i];
        const u64 right = (pre_[n_] - pre_[i]) - X * (u64)(n_ - i);
        return (long long)(left + right);
    }

    // ------------------------------------------------------------------ 相異なる値
    int distinct_count() const {
        ensure_uniq();
        return (int)uniq_.size();
    }
    int distinct_count_range(long long L, long long R) const {  // L <= v <= R の相異なる値の個数
        if (L > R) return 0;
        ensure_uniq();
        const int m = (int)uniq_.size();
        const int r = R == LLONG_MAX ? m : plain_lb(uniq_.data(), m, R + 1);
        return r - plain_lb(uniq_.data(), m, L);
    }
    int distinct_rank(long long x) const {  // x 未満の相異なる値の個数（x があれば座標圧縮後の番号）
        ensure_uniq();
        return plain_lb(uniq_.data(), (int)uniq_.size(), x);
    }
    long long distinct_kth(int k) const {  // k 番目に小さい相異なる値（0-indexed）
        ensure_uniq();
        assert(0 <= k && k < (int)uniq_.size());
        return uniq_[k];
    }
    const std::vector<long long>& distinct_values() const {
        ensure_uniq();
        return uniq_;
    }

    // ------------------------------------------------------------------ 条件で探す
    // [l, r) で f(a_i) が F…F T…T のとき、true となる最小の i（無ければ r）
    template <class F>
    int first_index_if(int l, int r, F f) const {
        const long long* a = A();
        int lo = l, len = r - l;
        while (len > 0) {
            const int half = len >> 1;
            if (f(a[lo + half])) {
                len = half;
            } else {
                lo += half + 1;
                len -= half + 1;
            }
        }
        return lo;
    }
    template <class F>
    int first_index_if(F f) const { return first_index_if(0, n_, f); }
    // [l, r) で f(a_i) が T…T F…F のとき、true となる最大の i（無ければ l - 1）
    template <class F>
    int last_index_if(int l, int r, F f) const {
        return first_index_if(l, r, [&](long long v) { return !f(v); }) - 1;
    }
    template <class F>
    int last_index_if(F f) const { return last_index_if(0, n_, f); }

    // ------------------------------------------------------------------ まとめて探索
    std::vector<int> lower_bound_many(const std::vector<long long>& xs) const { return many<false>(xs); }
    std::vector<int> upper_bound_many(const std::vector<long long>& xs) const { return many<true>(xs); }

    // ------------------------------------------------------------------ ペア（i < j）
    long long count_pairs_sum_le(long long K) const { return pairs_sum_le((i128)K); }
    long long count_pairs_sum_lt(long long K) const { return pairs_sum_le((i128)K - 1); }
    long long count_pairs_sum_ge(long long K) const { return total_pairs() - pairs_sum_le((i128)K - 1); }
    long long count_pairs_sum_gt(long long K) const { return total_pairs() - pairs_sum_le((i128)K); }
    long long count_pairs_diff_le(long long K) const { return pairs_diff_le((i128)K); }
    long long count_pairs_diff_lt(long long K) const { return pairs_diff_le((i128)K - 1); }
    long long count_pairs_diff_ge(long long K) const { return total_pairs() - pairs_diff_le((i128)K - 1); }
    long long count_pairs_diff_gt(long long K) const { return total_pairs() - pairs_diff_le((i128)K); }
    long long count_pairs_prod_le(long long K) const { return pairs_prod_le((i128)K); }
    long long count_pairs_prod_lt(long long K) const { return pairs_prod_le((i128)K - 1); }
    long long count_pairs_prod_ge(long long K) const { return total_pairs() - pairs_prod_le((i128)K - 1); }
    long long count_pairs_prod_gt(long long K) const { return total_pairs() - pairs_prod_le((i128)K); }

    // 値の区間 [v, v + W] に入る要素数の最大（W < 0 なら 0）
    int max_window_count(long long W) const {
        if (W < 0) return 0;
        const long long* a = A();
        int best = 0, j = 0;
        for (int i = 0; i < n_; ++i) {
            const long long t = clamp_ll((i128)a[i] + W);  // a_j <= t
            while (j < n_ && a[j] <= t) ++j;
            best = std::max(best, j - i);
        }
        return best;
    }
    // a_i + a_j（i < j）を小さい順に並べたときの k 番目（0-indexed）。和が long long に収まること
    long long kth_pair_sum(long long k) const {
        assert(n_ >= 2 && 0 <= k && k < total_pairs());
        const long long* a = A();
        return first_true(a[0] + a[1], a[n_ - 2] + a[n_ - 1], [&](long long s) { return count_pairs_sum_le(s) > k; });
    }
    // |a_i - a_j|（i < j）の k 番目に小さい値（0-indexed）。差が long long に収まること
    long long kth_pair_diff(long long k) const {
        assert(n_ >= 2 && 0 <= k && k < total_pairs());
        const long long* a = A();
        return first_true(0, a[n_ - 1] - a[0], [&](long long d) { return count_pairs_diff_le(d) > k; });
    }
    // a_i * a_j（i < j）の k 番目に小さい値（0-indexed）。積が long long に収まること
    long long kth_pair_prod(long long k) const {
        assert(n_ >= 2 && 0 <= k && k < total_pairs());
        const long long* a = A();
        // 積の最小・最大は、index {0, 1, n-2, n-1} の中の 2 つの積のどれか
        const int id[4] = {0, 1, n_ - 2, n_ - 1};
        i128 lo = 0, hi = 0;
        bool first = true;
        for (int p = 0; p < 4; ++p)
            for (int q = p + 1; q < 4; ++q) {
                if (id[p] == id[q]) continue;
                const i128 v = (i128)a[id[p]] * a[id[q]];
                if (first || v < lo) lo = v;
                if (first || v > hi) hi = v;
                first = false;
            }
        assert(LLONG_MIN <= lo && hi <= LLONG_MAX);
        return first_true((long long)lo, (long long)hi, [&](long long p) { return count_pairs_prod_le(p) > k; });
    }

    // ------------------------------------------------------------------ 2 つの配列（n × m 組）
    long long count_sum_le(const BinarySearch& b, long long K) const { return cross_sum_le(b, (i128)K); }
    long long count_sum_lt(const BinarySearch& b, long long K) const { return cross_sum_le(b, (i128)K - 1); }
    long long count_sum_ge(const BinarySearch& b, long long K) const {
        return (long long)n_ * b.n_ - cross_sum_le(b, (i128)K - 1);
    }
    long long count_sum_gt(const BinarySearch& b, long long K) const {
        return (long long)n_ * b.n_ - cross_sum_le(b, (i128)K);
    }
    // a_i + b_j の k 番目に小さい値（0-indexed）。和が long long に収まること
    long long kth_sum(const BinarySearch& b, long long k) const {
        assert(0 <= k && k < (long long)n_ * b.n_);
        return first_true(A()[0] + b.A()[0], A()[n_ - 1] + b.A()[b.n_ - 1],
                          [&](long long s) { return count_sum_le(b, s) > k; });
    }

  private:
    static constexpr int B = 8;      // 1 ノードのキー数（8 × 8 バイト = 64 バイト）
    static constexpr int MAXH = 16;  // 段数の上限（n < 2^31 なら 10 段以下）
#if defined(__AVX512F__) || defined(__AVX2__)
    static constexpr int G = 32;  // lower_bound_many で同時に下るクエリ数
#else
    static constexpr int G = 8;
#endif

    int n_ = 0, H_ = 1;
    std::array<std::size_t, MAXH> off_{};  // 段 h の先頭要素の位置（0 = 葉）
    std::vector<long long, bs_internal::AlignedAlloc<long long>> tr_;
    mutable std::vector<u64> pre_;  // 累積和（mod 2^64）。空なら未構築
    mutable std::vector<long long> uniq_;
    mutable bool uniq_built_ = false;

    const long long* A() const { return tr_.data() + off_[0]; }  // 葉の段 = ソート済み配列（a[n] は番兵）

    // ノード p（昇順の 8 キー）のうち x 未満のキーの個数。未満のキーは必ず先頭に並ぶ
    static inline int node_count(const long long* p, long long x) {
#if defined(__AVX512F__)
        const unsigned m = _mm512_cmpgt_epi64_mask(_mm512_set1_epi64(x), _mm512_load_si512((const void*)p));
        return __builtin_ctz(m + 1);
#elif defined(__AVX2__)
        const __m256i xv = _mm256_set1_epi64x(x);
        const unsigned m0 = (unsigned)_mm256_movemask_pd(
            _mm256_castsi256_pd(_mm256_cmpgt_epi64(xv, _mm256_load_si256((const __m256i*)p))));
        const unsigned m1 = (unsigned)_mm256_movemask_pd(
            _mm256_castsi256_pd(_mm256_cmpgt_epi64(xv, _mm256_load_si256((const __m256i*)(p + 4)))));
        return __builtin_ctz((m0 | (m1 << 4)) + 1);
#else
        return (p[0] < x) + (p[1] < x) + (p[2] < x) + (p[3] < x) + (p[4] < x) + (p[5] < x) + (p[6] < x) +
               (p[7] < x);
#endif
    }

    // 分岐なし二分探索（重複除去した配列用）
    static int plain_lb(const long long* a, int n, long long x) {
        if (n == 0) return 0;
        const long long* base = a;
        int len = n;
        while (len > 1) {
            const int half = len >> 1;
            base += (base[half - 1] < x) * half;
            len -= half;
        }
        return (int)(base - a) + (*base < x);
    }

    template <bool UP>
    std::vector<int> many(const std::vector<long long>& xs) const {
        const int q = (int)xs.size();
        std::vector<int> res(q);
        const long long* t = tr_.data();
        int i = 0;
        for (; i + G <= q; i += G) {
            long long key[G];
            std::size_t k[G];
            for (int g = 0; g < G; ++g) {
                const long long x = xs[i + g];
                key[g] = UP ? x + (x != LLONG_MAX) : x;
                k[g] = 0;
            }
            for (int h = H_ - 1; h > 0; --h) {
                const long long* base = t + off_[h];
                for (int g = 0; g < G; ++g) k[g] = k[g] * (B + 1) + node_count(base + k[g] * B, key[g]);
            }
            const long long* base = t + off_[0];
            for (int g = 0; g < G; ++g) res[i + g] = (int)(k[g] * B + node_count(base + k[g] * B, key[g]));
            if (UP)
                for (int g = 0; g < G; ++g)
                    if (xs[i + g] == LLONG_MAX) res[i + g] = n_;
        }
        for (; i < q; ++i) res[i] = UP ? upper_bound(xs[i]) : lower_bound(xs[i]);
        return res;
    }

    void ensure_pre() const {
        if (!pre_.empty()) return;
        pre_.resize((std::size_t)n_ + 1);
        const long long* a = A();
        u64 s = 0;
        pre_[0] = 0;
        for (int i = 0; i < n_; ++i) pre_[i + 1] = (s += (u64)a[i]);
    }
    void ensure_uniq() const {
        if (uniq_built_) return;
        uniq_.clear();
        const long long* a = A();
        for (int i = 0; i < n_; ++i)
            if (i == 0 || a[i] != a[i - 1]) uniq_.push_back(a[i]);
        uniq_built_ = true;
    }

    long long total_pairs() const { return (long long)n_ * (n_ - 1) / 2; }
    // 尺取り法の比較を long long 1 回で済ませるため、しきい値を long long の範囲に丸める
    static long long clamp_ll(i128 v) { return v > LLONG_MAX ? LLONG_MAX : v < LLONG_MIN ? LLONG_MIN : (long long)v; }
    long long pairs_sum_le(i128 K) const {  // a_i + a_j <= K（i < j）
        const long long* a = A();
        long long res = 0;
        int j = n_ - 1;
        for (int i = 0; i < j; ++i) {
            const i128 t = K - a[i];  // a_j <= t
            if (t < LLONG_MIN) break;  // この i 以降は相手がいない
            const long long tt = clamp_ll(t);
            while (j > i && a[j] > tt) --j;
            res += j - i;
        }
        return res;
    }
    long long pairs_diff_le(i128 K) const {  // a_j - a_i <= K（i < j）
        if (K < 0) return 0;
        const long long* a = A();
        long long res = 0;
        int i = 0;
        for (int j = 0; j < n_; ++j) {
            const i128 t = (i128)a[j] - K;  // a_i >= t
            if (t > LLONG_MAX) {
                i = j;
            } else {
                const long long tt = clamp_ll(t);
                while (a[i] < tt) ++i;
            }
            res += j - i;
        }
        return res;
    }
    long long pairs_prod_le(i128 K) const {  // a_i * a_j <= K（i < j）
        // 負 [0, ng)・0 [ng, ps)・正 [ps, n) に分け、組の種類ごとに尺取り法で数える（割り算なし・O(n)）
        const long long* a = A();
        const int n = n_, ng = lower_bound(0), ps = upper_bound(0);
        long long res = 0;
        {  // 正 × 正 : i が小さい方、j が大きい方
            int j = n - 1;
            for (int i = ps; i < j; ++i) {
                while (j > i && (i128)a[i] * a[j] > K) --j;
                res += j - i;
            }
        }
        {  // 負 × 負 : 絶対値の小さい方 i（右から）と大きい方 j（左から）
            int j = 0;
            for (int i = ng - 1; i > j; --i) {
                while (j < i && (i128)a[i] * a[j] > K) ++j;
                res += i - j;
            }
        }
        if (K >= 0) {  // 負 × 正（積 < 0）と 0 を含む組はすべて K 以下
            const long long z = ps - ng;
            res += (long long)ng * (n - ps) + z * (n - z) + z * (z - 1) / 2;
        } else {  // 負 × 正 : a_i が 0 に近づくほど条件を満たす a_j（正の側の右端部分）が減る
            int j = ps;
            for (int i = 0; i < ng; ++i) {
                while (j < n && (i128)a[i] * a[j] > K) ++j;
                res += n - j;
            }
        }
        return res;
    }
    long long cross_sum_le(const BinarySearch& o, i128 K) const {  // a_i + b_j <= K
        const long long* a = A();
        const long long* b = o.A();
        long long res = 0;
        int j = o.n_ - 1;
        for (int i = 0; i < n_; ++i) {
            const i128 t = K - a[i];  // b_j <= t
            if (t < LLONG_MIN) break;
            const long long tt = clamp_ll(t);
            while (j >= 0 && b[j] > tt) --j;
            if (j < 0) break;
            res += j + 1;
        }
        return res;
    }
};
#pragma endregion

#pragma region BinarySearch_test
// =============================================================================
//  ランダムテスト : 全機能を愚直解（全探索）と突き合わせる。失敗したら行番号を出して終了
// =============================================================================
namespace bs_test {

// NDEBUG（Release ビルド）でも無効にならないチェック
#define BS_CHECK(cond)                                                                        \
    do {                                                                                      \
        ++bs_test::g_checks;                                                                  \
        if (!(cond)) {                                                                        \
            std::cerr << "[BinarySearch] TEST FAILED  line " << __LINE__ << ": " #cond "\n"; \
            std::exit(1);                                                                     \
        }                                                                                     \
    } while (0)

// テスト量の倍率。最適化なし / デバッグビルドは遅いので自動で小さくする
#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
static const int TEST_SCALE = 1;
#else
static const int TEST_SCALE = 4;
#endif

using i128 = __int128;
using u64 = unsigned long long;
static long long g_checks = 0;
static std::mt19937_64 rng(20260926);

static long long rnd(long long lo, long long hi) {  // [lo, hi]
    const u64 w = (u64)hi - (u64)lo + 1;
    return w == 0 ? (long long)rng() : (long long)((u64)lo + rng() % w);
}
static long long sat_add(long long a, long long d) {  // 飽和加算
    const i128 v = (i128)a + d;
    return v > LLONG_MAX ? LLONG_MAX : v < LLONG_MIN ? LLONG_MIN : (long long)v;
}
static long long wrap(i128 v) { return (long long)(u64)v; }  // mod 2^64

// 値の分布（mode）
//   0: [-5, 5]  1: [0, 30]  2: [-1e9, 1e9]  3: [-1e18, 1e18]  4: 64bit 全域  5: 両端の値を多めに混ぜる
static const int MODES = 6;
static long long gen(int mode) {
    switch (mode) {
        case 0: return rnd(-5, 5);
        case 1: return rnd(0, 30);
        case 2: return rnd(-1000000000, 1000000000);
        case 3: return rnd(-1000000000000000000LL, 1000000000000000000LL);
        case 4: return (long long)rng();
        default: {
            static const long long sp[] = {LLONG_MIN, LLONG_MIN + 1, LLONG_MIN + 2, -2, -1, 0, 1, 2, LLONG_MAX - 2, LLONG_MAX - 1, LLONG_MAX};
            return rng() % 3 ? sp[rng() % 11] : (long long)rng();
        }
    }
}
// クエリ用の値: 分布どおり / 要素 ±1 / 両端
static long long gen_query(int mode, const std::vector<long long>& s) {
    const int t = (int)(rng() % 8);
    if (t < 3 || s.empty()) return gen(mode);
    if (t < 6) return sat_add(s[rng() % s.size()], rnd(-1, 1));
    static const long long sp[] = {LLONG_MIN, LLONG_MIN + 1, -1, 0, 1, LLONG_MAX - 1, LLONG_MAX};
    return sp[rng() % 7];
}

// ------------------------------------------------------------------ 愚直解
struct Naive {
    std::vector<long long> s;  // ソート済み
    int n() const { return (int)s.size(); }
    int lt(long long x) const { int c = 0; for (long long v : s) c += v < x; return c; }
    int le(long long x) const { int c = 0; for (long long v : s) c += v <= x; return c; }
    int eq(long long x) const { int c = 0; for (long long v : s) c += v == x; return c; }
    int in(long long L, long long R) const { int c = 0; for (long long v : s) c += (L <= v && v <= R); return c; }
    i128 sum_if(long long L, long long R) const { i128 t = 0; for (long long v : s) if (L <= v && v <= R) t += v; return t; }
};

// ------------------------------------------------------------------ 1 つの配列に対する全チェック
static void check_array(const BinarySearch& bs, const std::vector<long long>& sorted_v, int mode, int nq, bool pairs) {
    Naive nv{sorted_v};
    const int n = nv.n();
    // 基本
    BS_CHECK(bs.size() == n && bs.empty() == (n == 0));
    BS_CHECK(bs.to_vector() == sorted_v);
    for (int i = 0; i < n; ++i) BS_CHECK(bs[i] == sorted_v[i] && bs.kth(i) == sorted_v[i] && bs.kth_largest(i) == sorted_v[n - 1 - i]);
    if (n) BS_CHECK(bs.min() == sorted_v[0] && bs.max() == sorted_v[n - 1]);
    { int k = 0; for (long long v : bs) BS_CHECK(v == sorted_v[k++]); BS_CHECK(k == n); }
    BS_CHECK(bs.end() - bs.begin() == n && bs.data() == bs.begin());

    std::vector<long long> uq = sorted_v;
    uq.erase(std::unique(uq.begin(), uq.end()), uq.end());
    BS_CHECK(bs.distinct_count() == (int)uq.size() && bs.distinct_values() == uq);
    for (int k = 0; k < (int)uq.size(); ++k) BS_CHECK(bs.distinct_kth(k) == uq[k]);

    // 和（mod 2^64 で比較。答えが long long に収まれば完全一致を意味する）
    {
        i128 t = 0;
        BS_CHECK(bs.prefix_sum(0) == 0);
        for (int i = 0; i < n; ++i) { t += sorted_v[i]; BS_CHECK(bs.prefix_sum(i + 1) == wrap(t)); }
        BS_CHECK(bs.sum_all() == wrap(t));
        for (int k = 0; k <= n; ++k) {
            i128 a = 0, b = 0;
            for (int i = 0; i < k; ++i) a += sorted_v[i], b += sorted_v[n - 1 - i];
            BS_CHECK(bs.sum_smallest(k) == wrap(a) && bs.sum_largest(k) == wrap(b));
            if (k > 40) break;
        }
        for (int it = 0; it < 30; ++it) {
            int l = (int)rnd(0, n), r = (int)rnd(0, n);
            if (l > r) std::swap(l, r);
            i128 u = 0;
            for (int i = l; i < r; ++i) u += sorted_v[i];
            BS_CHECK(bs.sum_index(l, r) == wrap(u));
        }
    }

    // 値 x に対するクエリ
    for (int it = 0; it < nq; ++it) {
        const long long x = gen_query(mode, sorted_v);
        const int lt = nv.lt(x), le = nv.le(x), eq = nv.eq(x);
        BS_CHECK(bs.lower_bound(x) == lt);
        BS_CHECK(bs.upper_bound(x) == le);
        BS_CHECK(bs.equal_range(x) == std::make_pair(lt, le));
        BS_CHECK(bs.count_lt(x) == lt && bs.count_le(x) == le && bs.count_gt(x) == n - le && bs.count_ge(x) == n - lt);
        BS_CHECK(bs.count(x) == eq && bs.exists(x) == (eq > 0));
        BS_CHECK(bs.first_index(x) == (eq ? lt : -1) && bs.last_index(x) == (eq ? le - 1 : -1));
        // 前後の値
        {
            std::optional<long long> mlt, mle, mgt, mge;
            for (long long v : sorted_v) {
                if (v < x) mlt = v;
                if (v <= x) mle = v;
                if (v > x && !mgt) mgt = v;
                if (v >= x && !mge) mge = v;
            }
            BS_CHECK(bs.max_lt(x) == mlt && bs.max_le(x) == mle && bs.min_gt(x) == mgt && bs.min_ge(x) == mge);
            std::optional<long long> nr;
            i128 bd = -1;
            for (long long v : sorted_v) {
                const i128 d = v > x ? (i128)v - x : (i128)x - v;
                if (bd < 0 || d < bd || (d == bd && v < *nr)) bd = d, nr = v;
            }
            BS_CHECK(bs.nearest(x) == nr);
            const long long md = n == 0 ? LLONG_MAX : (bd > LLONG_MAX ? LLONG_MAX : (long long)bd);
            BS_CHECK(bs.min_dist(x) == md);
        }
        // 和
        {
            i128 slt = 0, sle = 0, sgt = 0, sge = 0, sad = 0;
            for (long long v : sorted_v) {
                if (v < x) slt += v;
                if (v <= x) sle += v;
                if (v > x) sgt += v;
                if (v >= x) sge += v;
                sad += v > x ? (i128)v - x : (i128)x - v;
            }
            BS_CHECK(bs.sum_lt(x) == wrap(slt) && bs.sum_le(x) == wrap(sle) && bs.sum_gt(x) == wrap(sgt) && bs.sum_ge(x) == wrap(sge));
            BS_CHECK(bs.sum_abs_diff(x) == wrap(sad));
        }
        // 相異なる値
        {
            int r = 0;
            for (long long v : uq) r += v < x;
            BS_CHECK(bs.distinct_rank(x) == r);
        }
        // 条件で探す（しきい値の述語）
        {
            BS_CHECK(bs.first_index_if([&](long long v) { return v >= x; }) == lt);
            BS_CHECK(bs.first_index_if([&](long long v) { return v > x; }) == le);
            BS_CHECK(bs.last_index_if([&](long long v) { return v < x; }) == lt - 1);
            BS_CHECK(bs.last_index_if([&](long long v) { return v <= x; }) == le - 1);
            int l = (int)rnd(0, n), r = (int)rnd(0, n);
            if (l > r) std::swap(l, r);
            int f = r, g = l - 1;
            for (int i = l; i < r; ++i) {
                if (sorted_v[i] >= x && f == r) f = i;
                if (sorted_v[i] <= x) g = i;
            }
            int calls = 0;
            BS_CHECK(bs.first_index_if(l, r, [&](long long v) { ++calls; return v >= x; }) == f);
            BS_CHECK(calls <= 32);
            BS_CHECK(bs.last_index_if(l, r, [&](long long v) { return v <= x; }) == g);
        }
    }
    // 値の区間 [L, R]
    for (int it = 0; it < nq; ++it) {
        long long L = gen_query(mode, sorted_v), R = gen_query(mode, sorted_v);
        if (rng() % 4 == 0) std::swap(L, R);  // L > R も試す
        const int c = nv.in(L, R);
        BS_CHECK(bs.count_range(L, R) == c && bs.exists_range(L, R) == (c > 0));
        const auto pr = bs.range_index(L, R);
        BS_CHECK(pr.second - pr.first == c && pr.first == nv.lt(L));
        if (c) BS_CHECK(sorted_v[pr.first] >= L && sorted_v[pr.second - 1] <= R);
        BS_CHECK(bs.sum_range(L, R) == wrap(nv.sum_if(L, R)));
        int dc = 0;
        for (long long v : uq) dc += (L <= v && v <= R);
        BS_CHECK(bs.distinct_count_range(L, R) == dc);
    }
    // まとめて探索（G の倍数と端数の両方を試す）
    {
        const int q = (int)rnd(0, 100);
        std::vector<long long> xs(q);
        for (auto& x : xs) x = gen_query(mode, sorted_v);
        const auto lo = bs.lower_bound_many(xs), up = bs.upper_bound_many(xs);
        BS_CHECK((int)lo.size() == q && (int)up.size() == q);
        for (int i = 0; i < q; ++i) BS_CHECK(lo[i] == nv.lt(xs[i]) && up[i] == nv.le(xs[i]));
    }
    if (!pairs) return;

    // ペア（i < j）: 全組を列挙して比較
    std::vector<i128> ps, pd, pp;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            ps.push_back((i128)sorted_v[i] + sorted_v[j]);
            pd.push_back((i128)sorted_v[j] - sorted_v[i]);
            pp.push_back((i128)sorted_v[i] * sorted_v[j]);
        }
    std::sort(ps.begin(), ps.end());
    std::sort(pd.begin(), pd.end());
    std::sort(pp.begin(), pp.end());
    auto cnt_le = [](const std::vector<i128>& v, i128 K) { return (long long)(std::upper_bound(v.begin(), v.end(), K) - v.begin()); };
    const long long tot = (long long)ps.size();
    auto pick = [&](const std::vector<i128>& v) -> long long {  // 境界付近の K を作る
        const int t = (int)(rng() % 4);
        if (t == 0 || v.empty()) return gen(mode);
        const i128 w = v[rng() % v.size()] + rnd(-1, 1);
        if (t == 3) { static const long long sp[] = {LLONG_MIN, -1, 0, 1, LLONG_MAX}; return sp[rng() % 5]; }
        return w > LLONG_MAX ? LLONG_MAX : w < LLONG_MIN ? LLONG_MIN : (long long)w;
    };
    for (int it = 0; it < 40; ++it) {
        long long K = pick(ps);
        BS_CHECK(bs.count_pairs_sum_le(K) == cnt_le(ps, K) && bs.count_pairs_sum_lt(K) == cnt_le(ps, (i128)K - 1));
        BS_CHECK(bs.count_pairs_sum_gt(K) == tot - cnt_le(ps, K) && bs.count_pairs_sum_ge(K) == tot - cnt_le(ps, (i128)K - 1));
        K = pick(pd);
        BS_CHECK(bs.count_pairs_diff_le(K) == cnt_le(pd, K) && bs.count_pairs_diff_lt(K) == cnt_le(pd, (i128)K - 1));
        BS_CHECK(bs.count_pairs_diff_gt(K) == tot - cnt_le(pd, K) && bs.count_pairs_diff_ge(K) == tot - cnt_le(pd, (i128)K - 1));
        K = pick(pp);
        BS_CHECK(bs.count_pairs_prod_le(K) == cnt_le(pp, K) && bs.count_pairs_prod_lt(K) == cnt_le(pp, (i128)K - 1));
        BS_CHECK(bs.count_pairs_prod_gt(K) == tot - cnt_le(pp, K) && bs.count_pairs_prod_ge(K) == tot - cnt_le(pp, (i128)K - 1));
        // 幅 W の窓
        const long long W = rng() % 3 ? pick(pd) : rnd(-2, 3);
        int best = 0;
        for (int i = 0; i < n; ++i) {
            int c = 0;
            for (int j = i; j < n; ++j) c += ((i128)sorted_v[j] - sorted_v[i] <= W);
            best = std::max(best, c);
        }
        BS_CHECK(bs.max_window_count(W) == (W < 0 ? 0 : best));
    }
    // k 番目（値が long long に収まるときだけ）
    auto fits = [](const std::vector<i128>& v) { return v.empty() || (v.front() >= LLONG_MIN && v.back() <= LLONG_MAX); };
    for (long long k = 0; k < tot; ++k) {
        if (tot > 60 && rng() % 8) continue;
        if (fits(ps)) BS_CHECK(bs.kth_pair_sum(k) == ps[k]);
        if (fits(pd)) BS_CHECK(bs.kth_pair_diff(k) == pd[k]);
        if (fits(pp)) BS_CHECK(bs.kth_pair_prod(k) == pp[k]);
    }
    // 2 つの配列
    {
        const int m = (int)rnd(0, 12);
        std::vector<long long> w(m);
        for (auto& x : w) x = gen(mode);
        BinarySearch b(w, false);
        std::sort(w.begin(), w.end());
        std::vector<i128> cs;
        for (long long x : sorted_v)
            for (long long y : w) cs.push_back((i128)x + y);
        std::sort(cs.begin(), cs.end());
        const long long nm = (long long)cs.size();
        for (int it = 0; it < 20; ++it) {
            const long long K = pick(cs);
            BS_CHECK(bs.count_sum_le(b, K) == cnt_le(cs, K) && bs.count_sum_lt(b, K) == cnt_le(cs, (i128)K - 1));
            BS_CHECK(bs.count_sum_gt(b, K) == nm - cnt_le(cs, K) && bs.count_sum_ge(b, K) == nm - cnt_le(cs, (i128)K - 1));
        }
        if (fits(cs))
            for (long long k = 0; k < nm; ++k)
                if (nm <= 40 || rng() % 8 == 0) BS_CHECK(bs.kth_sum(b, k) == cs[k]);
    }
}

// ------------------------------------------------------------------ 配列系のテスト
static void test_arrays() {
    // 小さい n（ノード境界 8, 72 をまたぐ）を大量に。ソート済みで渡す / 未ソートで渡すの両方
    for (int it = 0; it < 1500 * TEST_SCALE; ++it) {
        const int mode = it % MODES;
        const int n = (int)(it % 3 == 0 ? rnd(0, 80) : rnd(0, 20));
        std::vector<long long> v(n);
        for (auto& x : v) x = gen(mode);
        std::vector<long long> s = v;
        std::sort(s.begin(), s.end());
        const bool pre_sorted = it & 1;
        BinarySearch bs(pre_sorted ? s : v, pre_sorted);
        check_array(bs, s, mode, 40, true);
    }
    // 中くらいの n（段数 2～4 の木）: ペア以外を全部
    for (int it = 0; it < 60 * TEST_SCALE; ++it) {
        const int mode = it % MODES;
        static const int ns[] = {71, 72, 73, 647, 648, 649, 656, 657, 1000, 2047, 2048, 2049, 3000, 5831, 5832, 5833};
        const int n = it % 3 ? ns[rng() % 16] : (int)rnd(81, 6000);
        std::vector<long long> v(n);
        for (auto& x : v) x = gen(mode);
        std::vector<long long> s = v;
        std::sort(s.begin(), s.end());
        BinarySearch bs(v, false);
        check_array(bs, s, mode, 60, false);
    }
    // 大きい n: std::lower_bound / std::sort と比較（基数ソートの全パスと多段の木）
    for (int it = 0; it < 6 * TEST_SCALE; ++it) {
        const int mode = it % MODES;
        static const int ns[] = {52487, 52488, 52489, 52496, 200000, 472392, 472400};
        const int n = ns[it % 7];
        std::vector<long long> v(n);
        for (auto& x : v) x = gen(mode);
        std::vector<long long> s = v;
        std::sort(s.begin(), s.end());
        BinarySearch bs(v, false);
        BS_CHECK(std::equal(bs.begin(), bs.end(), s.begin()) && bs.size() == n);
        std::vector<long long> xs(20000);
        for (auto& x : xs) x = gen_query(mode, s);
        const auto lo = bs.lower_bound_many(xs), up = bs.upper_bound_many(xs);
        for (int i = 0; i < (int)xs.size(); ++i) {
            const int a = (int)(std::lower_bound(s.begin(), s.end(), xs[i]) - s.begin());
            const int b = (int)(std::upper_bound(s.begin(), s.end(), xs[i]) - s.begin());
            BS_CHECK(bs.lower_bound(xs[i]) == a && bs.upper_bound(xs[i]) == b && lo[i] == a && up[i] == b);
        }
    }
    // 基数ソート単体（ソート済み・逆順・全部同じ・2 値 なども）
    for (int it = 0; it < 40 * TEST_SCALE; ++it) {
        const int n = (int)rnd(2048, 30000);
        std::vector<long long> v(n);
        const int kind = it % 6;
        for (int i = 0; i < n; ++i) {
            if (kind == 0) v[i] = (long long)rng();
            else if (kind == 1) v[i] = rnd(-3, 3);
            else if (kind == 2) v[i] = 7;
            else if (kind == 3) v[i] = rng() & 1 ? LLONG_MIN : LLONG_MAX;
            else if (kind == 4) v[i] = n - i;  // 逆順
            else v[i] = (long long)(rng() >> rnd(0, 63)) * (rng() & 1 ? 1 : -1);
        }
        std::vector<long long> s = v;
        std::sort(s.begin(), s.end());
        bs_internal::sort_ll(v);
        BS_CHECK(v == s);
    }
    // コピー・作り直し（遅延構築した累積和 / 重複除去が正しく作り直されるか）
    {
        std::vector<long long> v = {5, 1, 4, 1, 5, 9, 2, 6};
        BinarySearch* p = new BinarySearch(v, false);
        BS_CHECK(p->sum_all() == 33 && p->distinct_count() == 6);
        BinarySearch c = *p;
        delete p;
        BS_CHECK(c.sum_all() == 33 && c.count_le(5) == 6 && c.distinct_rank(5) == 3);
        c.build({3, 3, 3}, true);
        BS_CHECK(c.sum_all() == 9 && c.distinct_count() == 1 && c.count(3) == 3 && c.size() == 3);
        BinarySearch e;
        BS_CHECK(e.size() == 0 && e.lower_bound(0) == 0 && e.sum_all() == 0 && !e.nearest(0) && e.min_dist(5) == LLONG_MAX);
        BinarySearch m(std::move(v), false);
        BS_CHECK(m.kth(0) == 1 && m.kth(7) == 9);
    }
}

// ------------------------------------------------------------------ 答えの二分探索（整数）
static void test_answer_search() {
    auto gen_range = [&](long long& lo, long long& hi) {
        const int t = (int)(rng() % 5);
        if (t == 0) lo = rnd(-20, 20), hi = rnd(-20, 20);
        else if (t == 1) lo = LLONG_MIN, hi = LLONG_MAX;
        else if (t == 2) lo = (long long)rng(), hi = (long long)rng();
        else if (t == 3) lo = LLONG_MIN + rnd(0, 3), hi = rnd(-5, 5);
        else lo = rnd(-5, 5), hi = LLONG_MAX - rnd(0, 3);
        if (lo > hi) std::swap(lo, hi);
    };
    for (int it = 0; it < 30000 * TEST_SCALE; ++it) {
        long long lo, hi;
        gen_range(lo, hi);
        // しきい値 t（範囲の外側・端・内側）
        long long t;
        const int k = (int)(rng() % 6);
        if (k == 0) t = lo;
        else if (k == 1) t = hi;
        else if (k == 2) t = sat_add(lo, 1);
        else if (k == 3) t = sat_add(hi, -1);
        else if (k == 4) t = (long long)rng();
        else t = rnd(lo, hi);
        const u64 len = (u64)hi - (u64)lo;
        int max_calls = 0;
        for (u64 w = len; w; w >>= 1) ++max_calls;  // ceil(log2(len + 1))
        // first_true : f(x) = (x >= t)
        {
            int calls = 0;
            bool inside = true;
            const long long r = first_true(lo, hi, [&](long long x) { ++calls; inside &= (lo <= x && x < hi); return x >= t; });
            const long long expect = t <= lo ? lo : (t >= hi ? hi : t);
            BS_CHECK(r == expect && inside && calls <= max_calls);
            // ギャロップ
            calls = 0;
            u64 far = 0;
            const long long g = first_true_gallop(lo, hi, [&](long long x) { ++calls; inside &= (lo <= x && x < hi); far = std::max(far, (u64)x - (u64)lo); return x >= t; });
            BS_CHECK(g == expect && inside);
            if (expect < hi) {
                const u64 d = (u64)expect - (u64)lo;
                BS_CHECK((i128)far <= 2 * (i128)d + 2);
                int lg = 0;
                for (u64 w = d + 1; w; w >>= 1) ++lg;
                BS_CHECK(calls <= 2 * lg + 2);
            } else {
                BS_CHECK(calls <= 2 * 64 + 2);
            }
        }
        // last_true : f(x) = (x <= t)（lo > LLONG_MIN のとき）
        if (lo > LLONG_MIN) {
            bool inside = true;
            const long long r = last_true(lo, hi, [&](long long x) { inside &= (lo <= x && x < hi); return x <= t; });
            const long long expect = lo >= hi ? lo - 1 : (t < lo ? lo - 1 : (t >= hi - 1 ? hi - 1 : t));
            BS_CHECK(r == expect && inside);
            BS_CHECK(last_true_gallop(lo, hi, [&](long long x) { return x <= t; }) == expect);
        }
        // めぐる式（両方向）。ok と ng では呼ばない
        if (lo < hi) {
            const long long tt = (long long)((u64)lo + (u64)rng() % len);  // [lo, hi - 1]
            bool strict = true;
            // ok = lo < ng = hi : f(x) = (x <= tt) → tt
            long long r = meguru(lo, hi, [&](long long x) { strict &= (lo < x && x < hi); return x <= tt; });
            BS_CHECK(r == tt && strict);
            // ok = hi > ng = lo : f(x) = (x >= tt + 1) → tt + 1
            r = meguru(hi, lo, [&](long long x) { strict &= (lo < x && x < hi); return x >= tt + 1; });
            BS_CHECK(r == tt + 1 && strict);
        }
    }
    // 空区間・1 点
    BS_CHECK(first_true(5, 5, [](long long) { return true; }) == 5);
    BS_CHECK(first_true(7, 3, [](long long) { return true; }) == 3);
    BS_CHECK(last_true(5, 5, [](long long) { return true; }) == 4);
    BS_CHECK(meguru(3, 4, [](long long) { return true; }) == 3 && meguru(4, 3, [](long long) { return true; }) == 4);
    BS_CHECK(first_true(LLONG_MAX - 1, LLONG_MAX, [](long long) { return false; }) == LLONG_MAX);
    BS_CHECK(first_true(LLONG_MIN, LLONG_MAX, [](long long x) { return x >= LLONG_MAX - 1; }) == LLONG_MAX - 1);
    BS_CHECK(first_true(LLONG_MIN, LLONG_MAX, [](long long) { return true; }) == LLONG_MIN);
}

// ------------------------------------------------------------------ 答えの二分探索（実数）
static double rnd_double() {
    const int t = (int)(rng() % 6);
    double d;
    if (t == 0) {  // 任意のビット列（NaN / inf 以外）
        do {
            const u64 u = rng();
            std::memcpy(&d, &u, 8);
        } while (!std::isfinite(d));
    } else if (t == 1) {
        d = std::ldexp((double)(rng() >> 11), -(int)rnd(0, 1100));  // 非正規化数まで
    } else if (t == 2) {
        d = (double)rnd(-1000000, 1000000) / 1000.0;
    } else if (t == 3) {
        static const double sp[] = {0.0, -0.0, 1.0, -1.0, 1e-308, 5e-324, -5e-324, 1e308, -1e308, 0.1};
        d = sp[rng() % 10];
    } else {
        d = std::ldexp((double)(rng() >> 11), (int)rnd(-60, 60) - 53);
    }
    return rng() & 1 ? -d : d;
}
static void test_real_search() {
    const double INF = std::numeric_limits<double>::infinity();
    for (int it = 0; it < 20000 * TEST_SCALE; ++it) {
        const double t = rnd_double();
        // lo < t <= hi
        double lo = rng() % 4 == 0 ? -INF : std::nextafter(t, -INF) - std::fabs(rnd_double());
        double hi = rng() % 4 == 0 ? INF : t + std::fabs(rnd_double());
        if (!(lo < t)) lo = -INF;
        int calls = 0;
        const double r = first_true_real(lo, hi, [&](double x) { ++calls; return x >= t; });
        BS_CHECK(r == t && calls <= 64);
        BS_CHECK(!(std::nextafter(r, -INF) >= t) || std::nextafter(r, -INF) == t);  // 1 つ手前は false（±0 は同じ値）
        // lo <= t < hi
        double lo2 = rng() % 4 == 0 ? -INF : t - std::fabs(rnd_double());
        double hi2 = rng() % 4 == 0 ? INF : std::nextafter(t, INF) + std::fabs(rnd_double());
        if (!(t < hi2)) hi2 = INF;
        calls = 0;
        const double r2 = last_true_real(lo2, hi2, [&](double x) { ++calls; return x <= t; });
        BS_CHECK(r2 == t && calls <= 64);
        BS_CHECK(meguru_real(hi, lo, [&](double x) { return x >= t; }) == t);
    }
    // 典型: sqrt(2)（隣り合う double まで厳密）
    const double s = first_true_real(0.0, 2.0, [](double x) { return x * x >= 2.0; });
    BS_CHECK(s * s >= 2.0 && std::nextafter(s, 0.0) * std::nextafter(s, 0.0) < 2.0);
}

// ------------------------------------------------------------------ 凸関数
static void test_convex() {
    for (int it = 0; it < 3000 * TEST_SCALE; ++it) {
        const long long lo = rnd(-50, 50), hi = lo + rnd(1, 80);
        // f(x) = Σ w_i |x - c_i| + a (x - c)^2（下に凸。平らな部分もある）
        const int m = (int)rnd(0, 4);
        std::vector<long long> w(m), c(m);
        for (int i = 0; i < m; ++i) w[i] = rnd(0, 3), c[i] = rnd(-60, 60);
        const long long a = rnd(0, 2), cc = rnd(-60, 60), b0 = rnd(-5, 5);
        auto f = [&](long long x) {
            long long v = a * (x - cc) * (x - cc) + b0 * 0;
            for (int i = 0; i < m; ++i) v += w[i] * (x > c[i] ? x - c[i] : c[i] - x);
            return v;
        };
        long long best = lo;
        for (long long x = lo; x < hi; ++x)
            if (f(x) < f(best)) best = x;
        BS_CHECK(argmin_convex(lo, hi, f) == best);
        BS_CHECK(argmax_concave(lo, hi, [&](long long x) { return -f(x); }) == best);
        BS_CHECK(argmin_convex(lo, hi, [&](long long x) { return (double)f(x); }) == best);
    }
    for (int it = 0; it < 2000 * TEST_SCALE; ++it) {
        const double c = (double)rnd(-1000000, 1000000) / 7.0, w = (double)rnd(1, 1000) / 10.0;
        const double lo = c - (double)rnd(1, 1000000) / 3.0, hi = c + (double)rnd(1, 1000000) / 3.0;
        const double x = argmin_convex_real(lo, hi, [&](double y) { return w * (y - c) * (y - c) + 3.0; });
        BS_CHECK(std::fabs(x - c) <= 1e-5 * std::max(1.0, std::fabs(c)));
        const double y = argmax_concave_real(lo, hi, [&](double z) { return -std::fabs(z - c) * w; });
        BS_CHECK(std::fabs(y - c) <= 1e-7 * std::max(1.0, std::fabs(c)) + 1e-9 * (hi - lo));
    }
}

// ------------------------------------------------------------------ 並列二分探索
static void test_parallel() {
    for (int it = 0; it < 400 * TEST_SCALE; ++it) {
        const int T = (int)rnd(0, 60), M = (int)rnd(1, 10), Q = (int)rnd(0, 50);
        std::vector<int> ep(T);
        std::vector<long long> ew(T);
        for (int t = 0; t < T; ++t) ep[t] = (int)rnd(0, M - 1), ew[t] = rnd(0, 5);
        std::vector<int> ql(Q), qr(Q);
        std::vector<long long> qs(Q);
        for (int q = 0; q < Q; ++q) {
            ql[q] = (int)rnd(0, M - 1), qr[q] = (int)rnd(ql[q] + 1, M);
            qs[q] = rnd(-1, 40);
        }
        // イベント t : 位置 ep[t] に ew[t] を足す。クエリ q : [ql, qr) の和が qs 以上になる最初の時刻
        std::vector<long long> cur(M);
        int applied = 0, resets = 0;
        const auto ans = parallel_binary_search(
            Q, T, [&]() { std::fill(cur.begin(), cur.end(), 0); applied = 0; ++resets; },
            [&](int t) { BS_CHECK(t == applied); cur[ep[t]] += ew[t]; ++applied; },
            [&](int q) { long long s = 0; for (int i = ql[q]; i < qr[q]; ++i) s += cur[i]; return s >= qs[q]; });
        for (int q = 0; q < Q; ++q) {
            std::vector<long long> b(M);
            int expect = T + 1;
            for (int t = 0; t <= T; ++t) {
                long long s = 0;
                for (int i = ql[q]; i < qr[q]; ++i) s += b[i];
                if (s >= qs[q]) { expect = t; break; }
                if (t < T) b[ep[t]] += ew[t];
            }
            BS_CHECK(ans[q] == expect);
        }
        int lg = 0;
        for (int w = T + 1; w; w >>= 1) ++lg;
        BS_CHECK(resets <= lg + 1);
    }
}

static void run_all() {
    const auto t0 = std::chrono::steady_clock::now();
    test_arrays();
    test_answer_search();
    test_real_search();
    test_convex();
    test_parallel();
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
#if defined(__AVX512F__)
    const char* path = "AVX-512";
#elif defined(__AVX2__)
    const char* path = "AVX2";
#else
    const char* path = "SIMD なし";
#endif
    std::cout << "[BinarySearch] random test: OK  (" << g_checks << " checks, " << std::fixed << std::setprecision(2) << sec
              << " s, 探索の実装 = " << path << ")\n";
    std::cout.unsetf(std::ios::fixed);
}

// ------------------------------------------------------------------ 速度計測
static void benchmark() {
    std::mt19937_64 r(1);
    const int ns[] = {200000, 1000000};
    for (int n : ns) {
        std::vector<long long> v(n);
        for (auto& x : v) x = (long long)(r() % 1000000000);
        const int Q = 2000000;
        std::vector<long long> xs(Q);
        for (auto& x : xs) x = (long long)(r() % 1000000000);
        auto ms = [](std::chrono::steady_clock::time_point a) { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count(); };
        auto t = std::chrono::steady_clock::now();
        std::vector<long long> s = v;
        std::sort(s.begin(), s.end());
        const double t_sort = ms(t);
        t = std::chrono::steady_clock::now();
        BinarySearch bs(v, false);
        const double t_build = ms(t);
        long long c1 = 0, c2 = 0, c3 = 0;
        t = std::chrono::steady_clock::now();
        for (int i = 0; i < Q; ++i) c1 += std::lower_bound(s.begin(), s.end(), xs[i]) - s.begin();
        const double t_std = ms(t);
        t = std::chrono::steady_clock::now();
        for (int i = 0; i < Q; ++i) c2 += bs.lower_bound(xs[i]);
        const double t_bs = ms(t);
        t = std::chrono::steady_clock::now();
        for (int x : bs.lower_bound_many(xs)) c3 += x;
        const double t_many = ms(t);
        BS_CHECK(c1 == c2 && c2 == c3);
        std::cout << std::fixed << std::setprecision(1) << "[BinarySearch] benchmark n = " << n << ", " << Q << " 回:"
                  << "  std::sort " << t_sort << " ms / 構築(ソート込み) " << t_build << " ms"
                  << "  |  std::lower_bound " << t_std * 1e6 / Q << " ns/回  bs.lower_bound " << t_bs * 1e6 / Q
                  << " ns/回  lower_bound_many " << t_many * 1e6 / Q << " ns/回\n";
        std::cout.unsetf(std::ios::fixed);
    }
}

}  // namespace bs_test
#pragma endregion

#pragma region BinarySearch_demo
// =============================================================================
//  使い方デモ : 全 API を 1 つずつ実行して「コード -> 返り値 // 意味」を表示する
// =============================================================================
namespace bs_demo {

static std::string to_s(bool b) { return b ? "true" : "false"; }
static std::string to_s(int v) { return std::to_string(v); }
static std::string to_s(long long v) { return std::to_string(v); }
static std::string to_s(double v) {
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}
static std::string to_s(const std::optional<long long>& v) { return v ? std::to_string(*v) : std::string("nullopt"); }
static std::string to_s(const std::pair<int, int>& p) { return "(" + std::to_string(p.first) + ", " + std::to_string(p.second) + ")"; }
template <class T>
static std::string to_s(const std::vector<T>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
static std::string pad(const std::string& s, size_t w) {  // 表示幅（全角 = 2）で揃える
    size_t width = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = (unsigned char)s[i];
        if ((c & 0xC0) == 0x80) continue;
        width += c >= 0xE0 ? 2 : 1;
    }
    return width >= w ? s : s + std::string(w - width, ' ');
}
static void banner(const std::string& s) {
    std::cout << "\n=====================================================================\n"
              << "  " << s << "\n"
              << "=====================================================================\n";
}
static void head(const std::string& s) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
static void note(const std::string& s) { std::cout << "    " << s << "\n"; }
static void code(const std::string& s) { std::cout << "  " << s << "\n"; }
// 文の末尾に ; を補う（if / for / while で } で終わる文には付けない）
static void stmt(const std::string& s) {
    const bool block = !s.empty() && s.back() == '}' &&
                       (s.compare(0, 3, "if ") == 0 || s.compare(0, 4, "for ") == 0 || s.compare(0, 6, "while ") == 0);
    code(block ? s : s + ";");
}
template <class V>
static void show(const char* expr, const V& v, const std::string& comment) {
    std::string line = "  " + pad(expr, 38) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 56) + "  // " + comment;
    std::cout << line << "\n";
}
// SHOW(式, "説明") : 式・返り値・説明を 1 行で表示
#define SHOW(expr, comment) bs_demo::show(#expr, (expr), comment)
// RUN(文) : 文を表示してから実行（宣言にも使える）
#define RUN(...)                 \
    bs_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void cheat_sheet() {
    head("API 一覧（詳しい説明はファイル先頭のコメント、使い方は下の【1】～【16】）");
    const char* rows[] = {
        "BinarySearch bs(A, is_sorted)      構築（false ならソートする）",
        "size empty [i] kth kth_largest min max begin end data to_vector",
        "lower_bound upper_bound equal_range range_index first_index last_index",
        "count_lt count_le count_gt count_ge count count_range exists exists_range",
        "max_lt max_le min_gt min_ge nearest（optional）  min_dist",
        "prefix_sum sum_index sum_all sum_lt/le/gt/ge sum_range sum_smallest sum_largest sum_abs_diff",
        "distinct_count distinct_count_range distinct_rank distinct_kth distinct_values",
        "first_index_if last_index_if（f は値を受け取る）",
        "lower_bound_many upper_bound_many",
        "count_pairs_{sum,diff,prod}_{le,lt,ge,gt}  max_window_count  kth_pair_{sum,diff,prod}",
        "count_sum_{le,lt,ge,gt}(b, K)  kth_sum(b, k)",
        "first_true last_true meguru first_true_gallop last_true_gallop（自由関数・整数）",
        "first_true_real last_true_real meguru_real（自由関数・実数）",
        "argmin_convex argmax_concave argmin_convex_real argmax_concave_real",
        "parallel_binary_search(Q, T, reset, apply, check)",
    };
    for (const char* r : rows) note(r);
}

static void demo() {
    banner("BinarySearch 使い方デモ    見方:  コード  ->  返り値    // 意味");
    cheat_sheet();

    // ------------------------------------------------------------------ 【1】
    head("【1】 構築  BinarySearch bs(A, is_sorted)");
    note("第 2 引数: A がソート済みなら true（assert で確認するだけ）、未ソートなら false（中でソート）");
    RUN(std::vector<long long> A = {5, 1, 4, 2, 3, 3, 9, 7});
    RUN(BinarySearch bs(A, false));
    SHOW(bs.to_vector(), "中身はソート済みの配列（A 自体は変わらない）");
    RUN(BinarySearch bs2({1, 2, 2, 8}, true));
    SHOW(bs2.size(), "ソート済みの配列はそのまま使う");
    note("大きな配列は BinarySearch bs(std::move(A), false); でコピーを省ける（A は空になる）");
    note("bs.build(新しい配列, is_sorted) で作り直せる。BinarySearch e; は空の配列");

    // ------------------------------------------------------------------ 【2】
    head("【2】 基本（O(1)）  ※ k は 0-indexed");
    SHOW(bs.size(), "要素数");
    SHOW(bs.empty(), "空か");
    SHOW(bs[2], "ソート後の 2 番目（0-indexed）");
    SHOW(bs.kth(0), "最小");
    SHOW(bs.kth(3), "小さい方から 4 番目 = kth(3)");
    SHOW(bs.kth_largest(0), "最大");
    SHOW(bs.kth_largest(1), "大きい方から 2 番目");
    SHOW(bs.min(), "最小値");
    SHOW(bs.max(), "最大値");
    RUN(long long total = 0; for (long long v : bs) total += v);
    SHOW(total, "範囲 for で全要素を回せる");

    // ------------------------------------------------------------------ 【3】
    head("【3】 位置（index）  配列 = {1, 2, 3, 3, 4, 5, 7, 9}");
    SHOW(bs.lower_bound(3), "3 以上が始まる位置（= 3 未満の個数）");
    SHOW(bs.upper_bound(3), "3 より大きいが始まる位置（= 3 以下の個数）");
    SHOW(bs.lower_bound(100), "無ければ n");
    SHOW(bs.equal_range(3), "3 がある範囲 [2, 4)");
    SHOW(bs.range_index(2, 5), "2 <= a_i <= 5 の index 範囲 [1, 6)（個数 = 5）");
    SHOW(bs.range_index(6, 6), "無ければ l == r");
    SHOW(bs.first_index(3), "3 が最初に出る位置");
    SHOW(bs.last_index(3), "3 が最後に出る位置");
    SHOW(bs.first_index(6), "無ければ -1");
    note("「x 未満で最大の要素の index」= lower_bound(x) - 1、「x 以下で最大」= upper_bound(x) - 1（無ければ -1）");

    // ------------------------------------------------------------------ 【4】
    head("【4】 個数（O(log n)）");
    SHOW(bs.count_lt(4), "4 未満の個数");
    SHOW(bs.count_le(4), "4 以下の個数");
    SHOW(bs.count_gt(4), "4 より大きい個数");
    SHOW(bs.count_ge(4), "4 以上の個数");
    SHOW(bs.count(3), "3 の個数");
    SHOW(bs.count_range(2, 5), "2 <= a_i <= 5 の個数（閉区間）");
    SHOW(bs.count_range(2, 5 - 1), "2 <= a_i < 5 は右端を -1 する");
    SHOW(bs.count_range(5, 2), "L > R なら 0");
    SHOW(bs.exists(7), "7 はある");
    SHOW(bs.exists(6), "6 はない");
    SHOW(bs.exists_range(6, 6), "[6, 6] に要素はない");
    SHOW(bs.exists_range(6, 8), "[6, 8] には 7 がある");

    // ------------------------------------------------------------------ 【5】
    head("【5】 前後の値（std::optional<long long>。無ければ nullopt）");
    SHOW(bs.max_lt(4), "4 未満の最大");
    SHOW(bs.max_le(4), "4 以下の最大");
    SHOW(bs.min_gt(4), "4 より大きい最小");
    SHOW(bs.min_ge(6), "6 以上の最小");
    SHOW(bs.max_lt(1), "1 未満は無い");
    SHOW(bs.max_lt(1).value_or(-1), "無いときの値は value_or で決める");
    SHOW(*bs.min_ge(6), "あると分かっていれば * で取り出す");
    SHOW(bs.nearest(6), "6 に最も近い値（5 と 7 は同じ距離 → 小さい方）");
    SHOW(bs.nearest(8), "8 に最も近い値（7 と 9 → 7）");
    SHOW(bs.min_dist(6), "min |a_i - 6|（旧 closest）");
    RUN(if (auto v = bs.min_gt(9)) { std::cout << "  あった: " << *v << "\n"; } else { std::cout << "  9 より大きい値は無い\n"; });

    // ------------------------------------------------------------------ 【6】
    head("【6】 和（初回だけ累積和 O(n)、以降 O(log n)）");
    SHOW(bs.sum_all(), "全体の和");
    SHOW(bs.prefix_sum(3), "a_0 + a_1 + a_2（ソート後）");
    SHOW(bs.sum_index(2, 5), "a_2 + a_3 + a_4");
    SHOW(bs.sum_lt(4), "4 未満の要素の和 = 1+2+3+3");
    SHOW(bs.sum_le(4), "4 以下の要素の和");
    SHOW(bs.sum_gt(4), "4 より大きい要素の和 = 5+7+9");
    SHOW(bs.sum_ge(4), "4 以上の要素の和");
    SHOW(bs.sum_range(3, 5), "3 <= a_i <= 5 の和 = 3+3+4+5");
    SHOW(bs.sum_smallest(3), "小さい方から 3 個の和");
    SHOW(bs.sum_largest(2), "大きい方から 2 個の和");
    SHOW(bs.sum_abs_diff(6), "Σ|a_i - 6|（全員を 6 に揃えるコスト）");
    SHOW(bs.sum_abs_diff(bs.kth(bs.size() / 2)), "中央値に揃えると最小");
    note("和は mod 2^64 で計算する。答えが long long に収まれば、途中であふれても正しい");

    // ------------------------------------------------------------------ 【7】
    head("【7】 相異なる値（初回だけ重複除去 O(n)）");
    SHOW(bs.distinct_count(), "相異なる値の個数（3 が 2 個）");
    SHOW(bs.distinct_values(), "相異なる値の一覧");
    SHOW(bs.distinct_rank(4), "4 未満の相異なる値の個数 = 座標圧縮後の 4 の番号");
    SHOW(bs.distinct_rank(6), "無い値でも「6 未満の相異なる値の個数」");
    SHOW(bs.distinct_kth(4), "4 番目に小さい相異なる値（0-indexed）");
    SHOW(bs.distinct_count_range(3, 7), "3 <= v <= 7 の相異なる値 {3,4,5,7}");

    // ------------------------------------------------------------------ 【8】
    head("【8】 条件で探す（f は値 a_i を受け取る。昇順に並べて F…F T…T / T…T F…F）");
    SHOW(bs.first_index_if([](long long v) { return v * v >= 20; }), "v^2 >= 20 となる最初の index（v = 5）");
    SHOW(bs.last_index_if([](long long v) { return v <= 4; }), "v <= 4 となる最後の index（T…T F…F）");
    SHOW(bs.first_index_if([](long long v) { return v > 100; }), "無ければ n");
    SHOW(bs.first_index_if(0, 3, [](long long v) { return v >= 3; }), "[0, 3) だけで探す");
    SHOW(bs.last_index_if(4, 8, [](long long v) { return v < 4; }), "[4, 8) に無ければ l - 1 = 3");
    note("index を受け取る述語で探したいときは自由関数 first_true(0, n, f) を使う（【12】）");

    // ------------------------------------------------------------------ 【9】
    head("【9】 まとめて探索（クエリが多いとき。1 回ずつより速い）");
    RUN(std::vector<long long> qs = {0, 3, 6, 10});
    SHOW(bs.lower_bound_many(qs), "それぞれの lower_bound（= 未満の個数）");
    SHOW(bs.upper_bound_many(qs), "それぞれの upper_bound（= 以下の個数）");

    // ------------------------------------------------------------------ 【10】
    head("【10】 ペア（i < j の組）  配列 = {1, 2, 3, 3, 4, 5, 7, 9}");
    SHOW(bs.count_pairs_sum_le(6), "a_i + a_j <= 6 の組の数（旧 count_pair_leq）");
    SHOW(bs.count_pairs_sum_lt(6), "a_i + a_j < 6");
    SHOW(bs.count_pairs_sum_ge(10), "a_i + a_j >= 10");
    SHOW(bs.count_pairs_sum_gt(10), "a_i + a_j > 10");
    SHOW(bs.count_pairs_diff_le(1), "|a_i - a_j| <= 1 の組の数");
    SHOW(bs.count_pairs_diff_lt(1), "|a_i - a_j| < 1 = 同じ値の組の数");
    SHOW(bs.count_pairs_diff_ge(5), "|a_i - a_j| >= 5");
    SHOW(bs.count_pairs_prod_le(6), "a_i * a_j <= 6（負の数・0 があってもよい）");
    SHOW(bs.count_pairs_prod_gt(20), "a_i * a_j > 20");
    SHOW(bs.max_window_count(2), "幅 2 の値の区間 [v, v+2] に入る最大個数（旧 max_len_diff_leq）");
    SHOW(bs.kth_pair_sum(0), "和の最小（1 + 2）");
    SHOW(bs.kth_pair_sum(5), "和を小さい順に並べた 5 番目（0-indexed）");
    SHOW(bs.kth_pair_diff(0), "差の最小（3 と 3）");
    SHOW(bs.kth_pair_prod(27), "積の最大（7 * 9。組は全部で 28 個）");
    RUN(BinarySearch neg({-3, -1, 0, 2, 5}, true));
    SHOW(neg.kth_pair_prod(0), "負の数を含む積の最小（-3 * 5）");
    SHOW(neg.kth_pair_prod(4), "積の 4 番目（ABC155 D Pairs 型）");

    // ------------------------------------------------------------------ 【11】
    head("【11】 2 つの配列（bs から 1 個、b から 1 個。n × m 組）");
    RUN(BinarySearch b({10, 20, 30}, true));
    SHOW(bs.count_sum_le(b, 13), "a_i + b_j <= 13 の組の数");
    SHOW(bs.count_sum_ge(b, 35), "a_i + b_j >= 35 の組の数");
    SHOW(bs.kth_sum(b, 0), "a_i + b_j の最小");
    SHOW(bs.kth_sum(b, 10), "a_i + b_j の 10 番目（0-indexed）");

    // ------------------------------------------------------------------ 【12】
    head("【12】 答えで二分探索（自由関数・整数）  f は long long を受け取る");
    note("first_true(lo, hi, f) : [lo, hi) で F…F T…T のとき true の最小（無ければ hi）");
    SHOW(first_true(0, 100, [](long long x) { return x * x >= 50; }), "x^2 >= 50 となる最小の x");
    SHOW(first_true(0, 5, [](long long x) { return x > 10; }), "無ければ hi");
    note("last_true(lo, hi, f)  : [lo, hi) で T…T F…F のとき true の最大（無ければ lo - 1）");
    SHOW(last_true(0, 2000000000, [](long long x) { return x * x <= 1000000000LL; }), "floor(sqrt(10^9))");
    note("meguru(ok, ng, f) : f(ok) = true, f(ng) = false とみなし（呼ばない）境界の ok を返す");
    SHOW(meguru(0, 2000000000, [](long long x) { return x * x <= 1000000000LL; }), "ok < ng（最大化）");
    SHOW(meguru(2000000000, 0, [](long long x) { return x * x >= 1000000000LL; }), "ok > ng（最小化）");
    note("long long 全域 [LLONG_MIN, LLONG_MAX) でもあふれない（中点を unsigned で計算）");
    SHOW(first_true(LLONG_MIN, LLONG_MAX, [](long long x) { return x >= -123; }), "全域で探しても最大 64 回");
    note("配列の index で探す（旧 BinarySearchArray の first_true / exists）");
    RUN(std::vector<long long> a = {1, 2, 2, 4, 6, 6, 6, 9});
    SHOW(first_true(0, (long long)a.size(), [&](long long i) { return a[i] >= 4; }), "a[i] >= 4 の最初の i");
    SHOW(first_true(2, 7, [&](long long i) { return a[i] >= 6; }), "[2, 7) で探す");
    SHOW(last_true(0, (long long)a.size(), [&](long long i) { return a[i] <= 6; }), "a[i] <= 6 の最後の i");
    SHOW(first_true(0, 3, [&](long long i) { return a[i] >= 4; }) < 3, "exists 相当: [0, 3) に true があるか");
    note("典型: 「1 日の上限 cap で 2 日以内に終わるか」を f にして、最小の cap を求める");
    RUN(std::vector<long long> work = {7, 2, 5, 10, 8});
    RUN(auto can = [&](long long cap) { long long days = 1, cur = 0; for (long long w : work) { if (w > cap) return false; if (cur + w > cap) ++days, cur = 0; cur += w; } return days <= 2; });
    SHOW(first_true(1, 1000, can), "仕事を順番どおり 2 日に分けるときの 1 日の量の最大の最小値");

    note("典型: 「x 以下の個数 >= k」で k 番目の値を求める（N × M の九九表で k 番目に小さい値）");
    RUN(const long long N = 3000, M = 2000, k = 1234567);
    RUN(auto cnt_le = [&](long long x) { long long c = 0; for (long long i = 1; i <= N; ++i) c += std::min(M, x / i); return c; });
    SHOW(first_true(1, N * M + 1, [&](long long x) { return cnt_le(x) >= k; }), "k 番目（1-indexed）に小さい積 i * j");

    // ------------------------------------------------------------------ 【13】
    head("【13】 ギャロップ探索 first_true_gallop / last_true_gallop（結果は first_true と同じ）");
    note("lo から 1, 2, 4, 8, … と幅を広げる。評価回数 O(log(答え - lo))、");
    note("f を呼ぶのは lo + 2*(答え - lo) + 2 以下だけ → 大きい x で f があふれる / hi が分からないとき");
    RUN(int calls = 0);
    SHOW(first_true_gallop(0, LLONG_MAX, [&](long long x) { ++calls; return x * x * x >= 1000; }), "x^3 >= 1000 の最小（x^3 があふれる所は呼ばない）");
    SHOW(calls, "f を呼んだ回数");
    SHOW(last_true_gallop(0, LLONG_MAX, [](long long x) { return x * x <= 30; }), "x^2 <= 30 の最大");

    // ------------------------------------------------------------------ 【14】
    head("【14】 実数で二分探索（double）  最大 64 回で隣り合う double まで厳密。回数を決めなくてよい");
    SHOW(first_true_real(0.0, 2.0, [](double x) { return x * x >= 2.0; }), "sqrt(2)：x^2 >= 2 の最小（f(lo)=false, f(hi)=true とみなす）");
    SHOW(last_true_real(0.0, 1e18, [](double x) { return x * x * x <= 10.0; }), "cbrt(10)：x^3 <= 10 の最大（f(lo)=true, f(hi)=false）");
    SHOW(meguru_real(1e9, 0.0, [](double x) { return x * std::log(x) >= 100.0; }), "めぐる式（ok = 1e9, ng = 0）");
    note("典型: 比の最大化「Σv / Σw >= x にできるか」⇔「v_i - x*w_i の大きい方から K 個の和 >= 0」");
    RUN(std::vector<double> v = {3, 5, 2}, w = {2, 4, 1});
    RUN(auto ratio_ok = [&](double x) { std::vector<double> d; for (size_t i = 0; i < v.size(); ++i) d.push_back(v[i] - x * w[i]); std::sort(d.rbegin(), d.rend()); return d[0] + d[1] >= 0; });
    SHOW(last_true_real(0.0, 100.0, ratio_ok), "2 個選んで Σv/Σw を最大化 → (3+2)/(2+1) = 5/3");
    note("±INFINITY を端に渡してもよい（例: first_true_real(-INFINITY, INFINITY, f)）");

    // ------------------------------------------------------------------ 【15】
    head("【15】 凸関数の最小・最大");
    SHOW(argmin_convex(-100, 101, [](long long x) { return (x - 7) * (x - 7) + 3; }), "[-100, 101) で (x-7)^2+3 を最小にする x");
    SHOW(argmin_convex(0, 100, [](long long x) { return std::abs(x - 10) + std::abs(x - 20); }), "平らな部分があれば最小の x");
    SHOW(argmax_concave(0, 50, [](long long x) { return x * (40 - x); }), "x(40 - x) を最大にする x");
    note("整数版は f(x+1) - f(x) の符号で二分探索するので、三分探索より評価が少なく必ず正しい");
    SHOW(argmin_convex_real(-10.0, 10.0, [](double x) { return (x - 1.5) * (x - 1.5); }), "実数版（黄金分割探索、100 回）");
    SHOW(argmax_concave_real(0.0, 3.0, [](double x) { return std::sin(x); }), "sin の最大点 π/2");

    // ------------------------------------------------------------------ 【16】
    head("【16】 並列二分探索 parallel_binary_search(Q, T, reset, apply, check)");
    note("T 個のイベントを順に適用するとき、クエリごとに「check が true になる最初の時刻」を一斉に二分探索する");
    note("返り値: 最初の t 個を適用した状態で true になる最小の t（0..T）。最後まで false なら T + 1");
    note("例: 辺を順に追加（イベント）し、頂点 u, v がつながる最初の時刻（クエリ）を求める");
    RUN(std::vector<std::pair<int, int>> edges = {{0, 1}, {2, 3}, {1, 2}, {3, 4}});
    RUN(std::vector<std::pair<int, int>> query = {{0, 3}, {2, 3}, {0, 4}, {4, 5}});
    RUN(atcoder::dsu uf(6));
    RUN(auto when = parallel_binary_search((int)query.size(), (int)edges.size(), [&]() { uf = atcoder::dsu(6); }, [&](int t) { uf.merge(edges[t].first, edges[t].second); }, [&](int q) { return uf.same(query[q].first, query[q].second); }));
    SHOW(when, "0-3 は 3 本目、2-3 は 2 本目、0-4 は 4 本目で連結。4-5 は連結しない（T + 1）");
    note("計算量: ceil(log2(T + 2)) 周 × (apply 最大 T 回 + check Q 回)");

    head("おわり（実戦ではファイル先頭の「#pragma region BinarySearch」だけをコピーする）");
}

}  // namespace bs_demo
#pragma endregion

void solv(){
    bs_test::run_all();    // ランダムテスト（全機能を愚直解と突き合わせる。失敗したら行番号を表示して終了）
    bs_test::benchmark();  // 速度計測（std::lower_bound / std::sort との比較）
    bs_demo::demo();       // 使い方デモ（全機能の「コード -> 返り値 // 意味」を表示）
    return;
}

int main(){
    BUFF;
    int T=1;
    //cin>>T;
    while(T--)solv();
    return 0;
}