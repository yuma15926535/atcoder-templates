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

// ============================================================================
//  StaticSeq : 静的な vector<long long> 専用 「全部入り」高速クエリテンプレート
//  ---------------------------------------------------------------------------
//  使い方 :   StaticSeq S(vector<long long>);   → あとは S.xxx(l, r, ...) を呼ぶだけ
//
//  ★ 共通の約束
//     ・配列の index は 0..N-1。
//     ・区間はすべて 半開区間 [l, r) で、l < r を仮定する。
//       (例外 : distinct_count は l>=r で 0、mode は l>=r で (0,0) を返す)
//     ・prefix_sum(r) / distinct_prefix(i) / palindrome_radius_even(i) /
//       first_sum_ge(l,x) / last_sum_ge(r,x) など一部は 0..N の位置を取る
//       (各関数のコメントを参照)。
//     ・範囲外の index / 区間に対するチェックは入っていない (assert もしない)。
//       assert があるのは「全要素が 0 以上」(XOR 系・線形基底) と
//       「modulus > 0」(剰余クエリ) の 2 つだけ。
//
//  【設計方針】
//   ・「静的（構築後に値が変わらない）」ことをフルに使って、動的構造では
//     出せない速度・機能（区間 mex / 区間最頻値 / 区間 distinct / 区間線形基底 /
//     部分列の O(1) 一致判定 など）まで提供する。
//   ・重い内部構造は すべて遅延構築。その機能を初めて呼んだ時に一度だけ作られる。
//     → 使わない機能の前処理コスト・メモリは一切かからない。
//   ・前処理は最大でも O(N sqrt N)（区間最頻値のみ）、他は O(N log N) 以下。
//     クエリは O(1) / O(log N) / O(sqrt N) / O(N)（下の API 一覧に明記）。
//   ・細部の高速化 :
//       - O(N) 構築 / O(1) クエリ の 64bit マスク RMQ
//       - rank 併載の簡潔ビットベクトル (全レベルを 1 本の連続配列に平坦化)
//       - 座標圧縮した Wavelet Matrix (レベル数 = ceil(log2(ユニーク値数)))
//       - predecessor / successor / 最頻値判定 / Σ|a_i-x| などを「1 回の下降」で処理
//       - しきい値検索 (next_ge 等) と累積和検索はブロック Sparse Table を 1 回下降
//       - 分岐の少ない二分探索、除算を使わない binary GCD、64x64→128 の乗算
//       - 座標圧縮はソート 1 回で rank / 出現位置 / CSR を同時生成
//
//  ★【遅延構築 : どの関数が何を作るか】初回呼び出し時に一度だけ構築される
//     build_rmin / build_rmax   … min max argmin argmax range_width minmax
//                                 next_ge/gt/le/lt prev_ge/gt/le/lt
//     build_prefix_agg          … prefix_/suffix_{min,max,argmin,argmax,gcd,and,or}
//     Sparse Table (gcd/and/or) … gcd(l,r) bit_and(l,r) bit_or(l,r)  ※演算ごとに別
//     build_compress            … 値と順位の変換、出現位置、全体情報の大半
//     build_wm                  … count_* kth median quantile prev/next_value nearest
//     build_wm_sum              … sum_lt/le/ge/gt sum_range sum_smallest/largest
//                                 sum_abs(l,r,x) median_cost
//     build_occ / build_wm_prev … prev/next/first_occurrence distinct_count
//     build_mex                 … mex(l,r) mex()
//     build_mode                … mode(l,r)
//     build_near                … prev/next_smaller/greater 系 cartesian_parent
//                                 minimum_containing_interval largest_rectangle
//     build_lift                … jump_next_greater / jump_next_smaller
//     build_run                 … is_non_decreasing 系 run_left/right/length
//                                 is_arithmetic is_constant is_sorted
//     build_psrmq / build_sst   … first/last_sum_ge/le max/min_subarray_sum
//                                 max/min_prefix_sum max/min_suffix_sum
//     build_wm_bits             … count_xor_* kth_xor max_xor min_xor (要 0 以上)
//     build_basis               … basis_* (要 0 以上)
//     build_suffix_array        … equal compare lcp_suffix lcp is_periodic
//                                 suffix_array suffix_rank distinct_subarray_count
//     build_manacher            … is_palindrome longest_palindrome
//                                 palindrome_radius_odd / _even
//     build_z / build_prefix_function / build_hash / build_moment /
//     build_sign / build_bitcnt / build_signchange … 対応する関数
//     build_all()               … 上記すべてを一度に構築する。
//                                 ※通常の競プロコードで呼ぶ必要は無い。
//                                   使う機能だけが自動で構築されるので、
//                                   build_all() は「まとめて先に作りたい」時だけ。
//
//  【前処理コストとメモリの目安 (N = 2*10^5, 値 < 2^30)】※呼んだ機能のぶんだけ
//     コンストラクタ                 O(N)          〜 3 MB
//     RMQ (min/max)                 O(N)          〜 4 MB   クエリ O(1)
//     Sparse Table (gcd/and/or)     O(N log N)    各 29 MB  クエリ O(1)
//     座標圧縮 + 出現位置            O(N log N)    〜 5 MB
//     Wavelet Matrix (個数/順位)     O(N log N)    〜 3 MB   クエリ O(log N)
//     Wavelet Matrix の区間和        O(N log N)    〜 30 MB
//     区間 distinct                  O(N log N)    〜 3 MB   クエリ O(log N)
//     区間 mex (永続セグ木)          O(N log N)    〜 48 MB  クエリ O(log N)
//     区間最頻値 (ブロック表)        O(N sqrt N)   〜 4 MB   クエリ O(sqrt N)
//     区間線形基底                   O(N log A)    〜 72 MB  クエリ O(log A)
//     Suffix Array + LCP            O(N log N)    〜 5 MB   クエリ O(1)
//     ぜんぶ (build_all)            0.5 秒程度     〜 350 MB
//   ※メモリは配列サイズから計算した目安。Wavelet Matrix の区間和・区間 mex・
//     区間線形基底は特に大きいので、使わないなら呼ばないこと。
//
//  【main() に入っているもの】
//     sademo::demo()   … 全機能の使い方を日本語解説つきで表示する動作確認関数
//     satest::run(seed)… 全機能を素朴解 (brute force) と突き合わせるランダムテスト
//   実戦で使うときは StaticSeq 本体 (先頭 〜 "StaticSeq ここまで" のコメント) だけを
//   コピーしてください。
//
//  ★【注意】
//   ・sum 系の戻り値は long long。オーバーフローが心配なときは sum128 / sum_sq
//     (__int128 を返す) を使ってください。
//     Wavelet Matrix の sum 系 (sum_lt など) も long long で持っています。
//   ・XOR 系クエリ・線形基底は「全要素が 0 以上」のときだけ使えます
//     (build_wm_bits / build_basis 内の assert(allNonNeg) で保護)。
//   ・NONE (= LLONG_MIN) は prev_value / next_value / nearest / majority /
//     second_min / second_max などの「該当なし」を表す番兵です。
//   ・hash_value は Rolling Hash なので理論上は衝突しえます。確実な一致判定が
//     必要なら equal() (Suffix Array + LCP, 決定的) を使ってください。
//   ・コピーしても安全です (内部にポインタを保持しません)。
//
//  【意図的に入れなかった機能】(オンラインで高速に answer できないため)
//   ・区間の転倒数 inversion_count(l,r) / 区間の大小ペア数 count_increasing_pairs(l,r)
//       … オンラインでは O(N sqrt N) メモリ級の構造が必要。全体の転倒数のみ提供。
//   ・区間の Σ freq_x^2 / 「ちょうど 1 回出現する値の個数」
//       … Mo's algorithm (オフライン) 専用の量。
//   ・区間内の最小差 min|a_i - a_j| / 和が x になるペア数 / 積 product(l,r)
//       … 高速化手段が無い、または long long でオーバーフローする。
//   ・LCM(l,r) … 即座にオーバーフローするため。GCD のみ提供。
//   ・任意 mod の区間積 / 素数判定カウント … 前処理が mod や値域に依存し汎用性が低い。
//   ・Mo's algorithm 等のオフライン処理 … 「配列そのものの機能」ではないため分離推奨。
// ============================================================================
//  【API 一覧】  N = 要素数, K = 相異なる値の個数, A = 値のビット幅
//    ・右端の計算量は「必要な内部構造が構築済みの状態での 1 回あたり」。
//      初回呼び出しのときだけ上の遅延構築コストが追加でかかります。
//    ・O(N) / O(N log N) と書いたものは呼ぶたびに毎回走査します
//      (結果をキャッシュしないので、ループ内で何度も呼ばないこと)。
// ----------------------------------------------------------------------------
//  ■ 基本         size() empty() get(i) [i] front() back() raw() diff(i)         O(1)
//                 get(l,r) … 部分配列を vector で返す                         O(r-l)
//  ■ 累積和       sum(l,r) prefix_sum(r) bit_xor(l,r)                            O(1)
//                 sum128 sum_sq sum_ia sum_i2a sum_sq_distance(l,r,x) ※__int128 O(1)
//                 count_positive/negative/zero/nonzero/odd/even(l,r)             O(1)
//                 sum_positive/negative/abs(l,r) sum_popcount count_bit(l,r,b)   O(1)
//  ■ 統計         mean variance stddev variance_index covariance_index
//                 correlation_index fit_line          ※long double で返す       O(1)
//  ■ min/max      min max argmin argmax minmax argminmax min_with_index
//                 max_with_index range_width                                     O(1)
//                 min_count max_count first_min last_min first_max last_max  O(log N)
//  ■ 区間集約     gcd(l,r) bit_and(l,r) bit_or(l,r)                              O(1)
//  ■ prefix/suffix prefix_/suffix_{min,max,argmin,argmax,gcd,and,or}(i)          O(1)
//  ■ 値の個数/和  count_lt/le/ge/gt/eq(l,r,x) count_range exists_range       O(log N)
//                 sum_lt/le/ge/gt(l,r,x) sum_range(l,r,x,y)                 O(log N)
//  ■ 順位統計     kth kth_largest median quantile sum_smallest sum_largest  O(log N)
//                 smallest_k(l,r,k) largest_k(l,r,k)  O(k + (取り出した異なる値の数)*log N)
//                 second_min second_max second_min/max_distinct             O(log N)
//  ■ 前後の値     prev_value prev_value_strict next_value next_value_strict
//                 max_le min_ge nearest sum_abs(l,r,x) median_cost majority O(log N)
//  ■ 静的専用     distinct_count distinct_prefix mex(l,r) mex()             O(log N)
//                 mode(l,r)                                                O(sqrt N)
//  ■ 出現位置     count(x) count(l,r,x) exists first_pos last_pos
//                 first_equal last_equal value_index rank_lb rank_ub        O(log N)
//                 prev/next/first_occurrence prev_distance next_distance        O(1)
//                 positions(x)                                   O(log N + 出現回数)
//                 first_duplicate()                                             O(N)
//  ■ 全体情報     unique_size compressed decompress is_all_distinct
//                 has_duplicate unique_values()                                 O(1)
//                 kth_global median_global lower/upper_bound_global          O(log N)
//                 max_frequency argmax_frequency frequency_table
//                 equal_pair_count                                              O(K)
//                 sorted() inversion_count() lis_length()                 O(N log N)
//  ■ 大小関係     prev/next_smaller(_equal) prev/next_greater(_equal)            O(1)
//                 cartesian_parent minimum_containing_interval                   O(1)
//                 jump_next_greater/smaller(i,k)                             O(log N)
//                 largest_rectangle()                                           O(N)
//  ■ 単調性/run   is_non_decreasing/non_increasing/strictly_increasing/
//                 strictly_decreasing/constant/arithmetic(l,r) is_sorted()       O(1)
//                 next_break_non_decreasing run_left run_right run_length       O(1)
//                 is_local_min/max count_sign_changes                            O(1)
//                 runs() peaks() valleys() difference() second_difference()     O(N)
//  ■ しきい値検索 next_ge/gt/le/lt(i,x)  prev_ge/gt/le/lt(i,x)               O(log N)
//  ■ 累積和検索   first_sum_ge/le(l,x) last_sum_ge/le(r,x)                   O(log N)
//                 max_subarray_sum min_subarray_sum                        O(log N)
//                 max/min_prefix_sum max/min_suffix_sum                        O(1)
//  ■ XOR          count_xor_lt/le/ge/gt(l,r,x,bound) kth_xor max_xor min_xor O(log A)
//  ■ 線形基底     basis_rank basis_max_xor basis_contains                    O(log A)
//                 basis_kth_xor                                           O(log^2 A)
//  ■ 部分列       equal(l1,r1,l2,r2) compare(...) lcp_suffix lcp compare_suffix
//                 suffix_rank is_periodic hash_value(l,r) suffix_array()
//                 min_suffix max_suffix                                         O(1)
//                 distinct_subarray_count()                                     O(N)
//  ■ 回文/文字列  is_palindrome palindrome_radius_odd/even                       O(1)
//                 z_array() prefix_function()  ※参照を返す                      O(1)
//                 longest_palindrome()                                          O(N)
//  ■ 剰余         count_mod sum_mod count_divisible sum_divisible           O(log N)
//                 ※ modulus を前回と変えるたびに O(N log N) で作り直す
//  ■ まとめ構築   build_all()  … すべての内部構造をまとめて構築 (通常は不要)
// ============================================================================
using i128 = __int128_t;
using u128 = __uint128_t;

// __int128 を 10 進文字列にする (表示用)。O(桁数)。
// LLONG_MIN 相当の最小値でも安全に扱えるよう、符号反転は (-(x+1))+1 で行う。
static string i128str(i128 x) {
    if (x == 0) return "0";
    bool neg = (x < 0);
    u128 v = neg ? (u128)(-(x + 1)) + 1 : (u128)x;
    string s;
    while (v) { s += char('0' + (int)(v % 10)); v /= 10; }
    if (neg) s += '-';
    reverse(s.begin(), s.end());
    return s;
}

// gcd(|x|, |y|) を返す (戻り値は常に 0 以上)。除算を使わない binary GCD。
//   gcd(0,0) = 0、gcd(0,y) = |y|。x = LLONG_MIN でも未定義動作にならない。
//   計算量 O(log(min(|x|,|y|)))
static inline ll ll_gcd(ll x, ll y) {
    ull a = (x < 0) ? (ull)(-(x + 1)) + 1 : (ull)x;
    ull b = (y < 0) ? (ull)(-(y + 1)) + 1 : (ull)y;
    if (a == b) return (ll)a;                  // 同値 (0,0 を含む) は即決
    if (a == 1 || b == 1) return 1;            // 片方が 1 なら即決 (Sparse Table の上位段は大半がこれ)
    if (!a) return (ll)b;
    if (!b) return (ll)a;
    int sh = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do {
        b >>= __builtin_ctzll(b);
        if (a > b) { ull t = a; a = b; b = t; }
        b -= a;
    } while (b);
    return (ll)(a << sh);
}

// 分岐の少ない二分探索 (cmov 化されやすい形)。計算量 O(log n)。
//   ★前提 : a[0..n) は昇順ソート済みであること。n == 0 のときは 0 を返す。
//   lb_idx : a[j] >= x となる最小の j (std::lower_bound と同じ)
//   ub_idx : a[j] >  x となる最小の j (std::upper_bound と同じ)
template <class T>
static inline int lb_idx(const T *a, int n, T x) {
    int base = 0;
    while (n > 1) {
        int half = n >> 1;
        base += (a[base + half - 1] < x) ? half : 0;
        n -= half;
    }
    return base + ((n > 0 && a[base] < x) ? 1 : 0);
}
template <class T>
static inline int ub_idx(const T *a, int n, T x) {
    int base = 0;
    while (n > 1) {
        int half = n >> 1;
        base += (a[base + half - 1] <= x) ? half : 0;
        n -= half;
    }
    return base + ((n > 0 && a[base] <= x) ? 1 : 0);
}

// ============================================================================
//  簡潔ビットベクトルのブロック : 64bit の本体 b と、その手前までの 1 の個数 c を
//  1 つの構造体に同居させ、1 回のキャッシュライン読み込みで rank が求まるようにする。
//  (pad は 16 バイト境界に揃えるための詰め物で、値は使用しない)
// ============================================================================
struct BitBlk { ull b; unsigned c; unsigned pad; };

// 上の BitBlk を使った素直なビットベクトル : 構築 O(N) / rank O(1)。
//   set() でビットを立てたあと必ず build() を呼んでから rank1/rank0 を使う。
//   rank1(i) = [0,i) に含まれる 1 の個数、rank0(i) = [0,i) に含まれる 0 の個数。
//   ※ Wavelet Matrix は全レベルを 1 本の配列に平坦化した独自の持ち方をしているため、
//     この構造体は現在どこからも使われていない (参考実装として残してある)。
struct BitVec {
    int n = 0;
    vector<BitBlk> v;
    void init(int n_) { n = n_; v.assign((size_t)(n >> 6) + 2, BitBlk{0ULL, 0u, 0u}); }
    inline void set(int i) { v[i >> 6].b |= 1ULL << (i & 63); }
    inline int get(int i) const { return (int)((v[i >> 6].b >> (i & 63)) & 1ULL); }
    void build() {
        unsigned s = 0;
        for (auto &blk : v) { unsigned p = (unsigned)__builtin_popcountll(blk.b); blk.c = s; s += p; }
    }
    inline int rank1(int i) const {
        const BitBlk &blk = v[i >> 6];
        return (int)(blk.c + (unsigned)__builtin_popcountll(blk.b & ((1ULL << (i & 63)) - 1)));
    }
    inline int rank0(int i) const { return i - rank1(i); }
};

// ============================================================================
//  RMQ (64bit マスク + ブロック Sparse Table 方式)
//    build   : O(N) 時間 / O(N) メモリ
//    query   : O(1)      … [l,r) の最良位置
//    findFirst / findLast : O(log N) … しきい値を満たす最初/最後の位置
//    MAXI=false → 最小値の位置, MAXI=true → 最大値の位置
//    ★同値のときは常に「最も左」の位置を返す。
//    ★値配列は保持せず、build / query / findFirst / findLast に毎回ポインタを渡す。
//      (StaticSeq をコピーしても内部ポインタが壊れないようにするため)
// ============================================================================
template <class T, bool MAXI>
struct RMQ64 {
    int n = 0, nb = 0, K = 1;
    vector<ull> mask;
    vector<int> sp;

    static inline bool sBetter(const T *v, int i, int j) { return MAXI ? (v[i] > v[j]) : (v[i] < v[j]); }
    static inline int pick(const T *v, int i, int j) {
        if (sBetter(v, i, j)) return i;
        if (sBetter(v, j, i)) return j;
        return i < j ? i : j;              // 同値は左優先
    }
    // [r-sz+1, r] の最良位置 (1 <= sz <= 64, r-sz+1 >= 0 が前提)。mask だけで O(1)。
    inline int smallIdx(int r, int sz) const {
        ull m = mask[r];
        if (sz < 64) m &= (1ULL << sz) - 1;
        return r - (63 - __builtin_clzll(m));
    }
    void build(const T *v, int n_) {
        n = n_;
        mask.assign(n, 0ULL);
        ull cur = 0;
        for (int i = 0; i < n; i++) {
            cur = (cur << 1) | 1ULL;
            ull t = cur & ~1ULL;
            while (t) {
                int k = __builtin_ctzll(t);
                if (sBetter(v, i, i - k)) { cur ^= 1ULL << k; t ^= 1ULL << k; }
                else break;
            }
            mask[i] = cur;
        }
        nb = (n + 63) >> 6;
        K = 1; while ((1 << K) <= nb) K++;
        sp.assign((size_t)K * nb, 0);
        for (int b = 0; b < nb; b++) {
            int r = min((b << 6) + 63, n - 1);
            sp[b] = smallIdx(r, r - (b << 6) + 1);
        }
        for (int k = 1; k < K; k++)
            for (int b = 0; b + (1 << k) <= nb; b++)
                sp[(size_t)k * nb + b] = pick(v, sp[(size_t)(k - 1) * nb + b],
                                                 sp[(size_t)(k - 1) * nb + b + (1 << (k - 1))]);
    }
    // [l,r) 内の最良位置を返す。★l < r が前提 (空区間は不可)。O(1)。
    int query(const T *v, int l, int r) const {
        if (r - l <= 64) return smallIdx(r - 1, r - l);
        int res = pick(v, smallIdx(l + 63, 64), smallIdx(r - 1, 64));
        int lb = (l >> 6) + 1, rb = (r >> 6);
        if (lb < rb) {
            int k = 31 - __builtin_clz((unsigned)(rb - lb));
            res = pick(v, res, pick(v, sp[(size_t)k * nb + lb], sp[(size_t)k * nb + rb - (1 << k)]));
        }
        return res;
    }
    // 条件判定 : MAXI なら v[j] >= x (strict なら v[j] > x)
    //            MAXI でなければ v[j] <= x (strict なら v[j] < x)
    inline bool ok(T val, T x, bool strict) const {
        if (MAXI) return strict ? (val > x) : (val >= x);
        else      return strict ? (val < x) : (val <= x);
    }
    // j >= i で ok(v[j], x, strict) を満たす最小の j。無ければ n を返す。O(log N)。
    //   i >= n のときは n を返す。ブロック内は mask、ブロック間は Sparse Table を
    //   1 回下降して探すので、O(1) RMQ を log 回呼ぶより速い。
    int findFirst(const T *v, int i, T x, bool strict) const {
        if (i >= n) return n;
        int b = i >> 6;
        int be = min(n - 1, (b << 6) + 63);
        if (ok(v[smallIdx(be, be - i + 1)], x, strict)) {
            int lo = i, hi = be;
            while (lo < hi) { int m = (lo + hi) >> 1; if (ok(v[smallIdx(m, m - i + 1)], x, strict)) hi = m; else lo = m + 1; }
            return lo;
        }
        int c = b + 1;
        if (c >= nb) return n;
        for (int k = K - 1; k >= 0; k--)
            if (c + (1 << k) <= nb && !ok(v[sp[(size_t)k * nb + c]], x, strict)) c += (1 << k);
        if (c >= nb || !ok(v[sp[c]], x, strict)) return n;
        int st = c << 6, en = min(n - 1, st + 63);
        int lo = st, hi = en;
        while (lo < hi) { int m = (lo + hi) >> 1; if (ok(v[smallIdx(m, m - st + 1)], x, strict)) hi = m; else lo = m + 1; }
        return lo;
    }
    // j <= i で ok(v[j], x, strict) を満たす最大の j。無ければ -1 を返す。O(log N)。
    //   i < 0 のときは -1 を返す。
    int findLast(const T *v, int i, T x, bool strict) const {
        if (i < 0) return -1;
        int b = i >> 6, st = b << 6;
        if (ok(v[smallIdx(i, i - st + 1)], x, strict)) {
            int lo = st, hi = i;
            while (lo < hi) { int m = (lo + hi + 1) >> 1; if (ok(v[smallIdx(i, i - m + 1)], x, strict)) lo = m; else hi = m - 1; }
            return lo;
        }
        if (b == 0) return -1;
        int cend = b;                           // 探索範囲は ブロック [0, cend)
        for (int k = K - 1; k >= 0; k--)
            if (cend - (1 << k) >= 0 && !ok(v[sp[(size_t)k * nb + cend - (1 << k)]], x, strict)) cend -= (1 << k);
        if (cend <= 0) return -1;
        int blk = cend - 1;
        if (!ok(v[sp[blk]], x, strict)) return -1;
        int bst = blk << 6, ben = min(n - 1, bst + 63);
        int lo = bst, hi = ben;
        while (lo < hi) { int m = (lo + hi + 1) >> 1; if (ok(v[smallIdx(ben, ben - m + 1)], x, strict)) lo = m; else hi = m - 1; }
        return lo;
    }
};

// ============================================================================
//  冪等演算用 Sparse Table : 構築 O(N log N) 時間 / O(N log N) メモリ、クエリ O(1)
//    ★重なりを許して 2 区間を合成するので、F::op は冪等 (op(x,x) == x) かつ
//      結合的・可換である必要がある (gcd / AND / OR は満たす。和は不可)。
//    使い方は build(a) 1 回、または level0(n) で領域を得て直接書き込んでから finish()。
// ============================================================================
template <class F>
struct SparseOp {
    int n = 0, K = 1;
    vector<ll> t;
    // レベル 0 (元配列) を直接書き込むための領域を返す。書き終えたら finish() を呼ぶ。
    // 一時 vector を作らずに済むので、絶対値を取るなどの前処理があるときに使う。
    ll *level0(int n_) {
        n = n_; K = 1; while ((1 << K) <= n) K++;
        t.resize((size_t)K * n);
        return t.data();
    }
    void finish() {
        for (int k = 1; k < K; k++) {
            const ll *p = t.data() + (size_t)(k - 1) * n;
            ll *q = t.data() + (size_t)k * n;
            int half = 1 << (k - 1), lim = n - (1 << k);
            for (int i = 0; i <= lim; i++) q[i] = F::op(p[i], p[i + half]);
        }
    }
    void build(const vector<ll> &a) {
        ll *p = level0((int)a.size());
        copy(a.begin(), a.end(), p);
        finish();
    }
    ll query(int l, int r) const {             // ★[l,r) の集約。l < r が前提。O(1)
        int k = 31 - __builtin_clz((unsigned)(r - l));
        const ll *p = t.data() + (size_t)k * n;
        return F::op(p[l], p[r - (1 << k)]);
    }
};
struct OpGcd { static ll op(ll x, ll y) { return ll_gcd(x, y); } };   // 絶対値の gcd
struct OpAnd { static ll op(ll x, ll y) { return x & y; } };          // ビット AND
struct OpOr  { static ll op(ll x, ll y) { return x | y; } };          // ビット OR

// ============================================================================
//  Wavelet Matrix (静的な整数列に対する 2 次元 (位置, 値) クエリの土台)
//    ★前提 : build に渡す keys は 0 <= key < 2^lg の非負整数列であること。
//             (StaticSeq では「座標圧縮した順位」か「非負の生の値」を渡している)
//    ★区間はすべて 半開区間 [l, r)。位置は 0..n-1。
//    build     : O(N * lg) 時間 / ビットベクトルは O(N * lg / 8) バイト程度
//    buildSum  : O(N * lg) 時間 / (lg+1) * (N+1) 個の long long を追加で使う
//                (sum 系クエリを使うときだけ呼べばよい)
//    各クエリ  : O(lg)  ※dfsSmall / dfsLarge だけは列挙する値の個数に比例
//    全レベルのビットベクトルを 1 本の連続配列に平坦化し、
//    l と r の rank をまとめて求めることでメモリアクセスを減らしている。
// ============================================================================
struct WaveletMatrix {
    int n = 0, lg = 0, nb = 0;                // lg = ビット数, nb = 1 レベルあたりのブロック数
    vector<BitBlk> bit;                       // 全レベル分を平坦化 (レベル t は bit[t*nb ..])
    vector<int> mid;                          // mid[t] = レベル t で bit が 0 の要素数
    vector<vector<ll>> S;                     // S[t][i] : レベル t 時点の並びでの重み累積和
    bool hasSum = false;                      // buildSum 済みか (sum 系はこれが true のとき有効)

    inline const BitBlk *lvl(int t) const { return bit.data() + (size_t)t * nb; }  // レベル t の先頭
    // レベル B の [0,i) に含まれる 1 の個数。O(1)。
    static inline int rk1(const BitBlk *B, int i) {
        const BitBlk &x = B[i >> 6];
        return (int)(x.c + (unsigned)__builtin_popcountll(x.b & ((1ULL << (i & 63)) - 1)));
    }
    // [0,l) と [0,r) の「0 の個数」を同時に求める (l0, r0 に格納)。O(1)。
    //   分岐を入れずに 2 つのブロックを読むことで、2 本のメモリアクセスを
    //   CPU が並行して発行できるようにしている (同一ブロックなら 2 本目は L1 ヒット)
    static inline void rk0pair(const BitBlk *B, int l, int r, int &l0, int &r0) {
        const BitBlk &X = B[l >> 6];
        const BitBlk &Y = B[r >> 6];
        int a1 = (int)(X.c + (unsigned)__builtin_popcountll(X.b & ((1ULL << (l & 63)) - 1)));
        int b1 = (int)(Y.c + (unsigned)__builtin_popcountll(Y.b & ((1ULL << (r & 63)) - 1)));
        l0 = l - a1; r0 = r - b1;
    }
    // レベル t の位置 i のビット (0 or 1)。O(1)。
    inline int getBit(int t, int i) const {
        const BitBlk &x = lvl(t)[i >> 6];
        return (int)((x.b >> (i & 63)) & 1ULL);
    }

    // keys 列 cur (0 <= key < 2^lg_) から本体を構築する。lg_ は 1 以上に丸められる。
    // O(N * lg) 時間。cur は値渡しなので呼び出し側の配列は変化しない。
    void build(vector<ll> cur, int lg_) {
        n = (int)cur.size(); lg = max(lg_, 1);
        nb = (n >> 6) + 2;
        bit.assign((size_t)lg * nb, BitBlk{0ULL, 0u, 0u});
        mid.assign(lg, 0);
        vector<ll> nxt(n);
        for (int t = 0; t < lg; t++) {
            int b = lg - 1 - t;
            BitBlk *B = bit.data() + (size_t)t * nb;
            int z = 0;
            const ll *c = cur.data();
            for (int i = 0; i < n; i++) z += (int)(~(c[i] >> b) & 1);
            mid[t] = z;
            int p0 = 0, p1 = z;
            ll *d = nxt.data();
            for (int i = 0; i < n; i++) {
                ll val = c[i];
                if ((val >> b) & 1) { B[i >> 6].b |= 1ULL << (i & 63); d[p1++] = val; }
                else                { d[p0++] = val; }
            }
            unsigned s = 0;
            for (int i = 0; i < nb; i++) { unsigned p = (unsigned)__builtin_popcountll(B[i].b); B[i].c = s; s += p; }
            swap(cur, nxt);
        }
    }
    // sum 系クエリ用の累積和を追加構築する。★build() の後に呼ぶこと。
    //   cur : 元の並びでの重み列 (StaticSeq では元の値 a をそのまま渡している)
    //   O(N * lg) 時間、(lg+1)*(N+1) 個の long long を確保する。
    void buildSum(vector<ll> cur) {
        S.assign(lg + 1, vector<ll>());
        vector<ll> nxt(n);
        for (int t = 0; t < lg; t++) {
            auto &s = S[t];
            s.resize(n + 1); s[0] = 0;
            const ll *c = cur.data();
            for (int i = 0; i < n; i++) s[i + 1] = s[i] + c[i];
            int p0 = 0, p1 = mid[t];
            ll *d = nxt.data();
            const BitBlk *B = lvl(t);
            for (int i = 0; i < n; i++) {
                if ((B[i >> 6].b >> (i & 63)) & 1ULL) d[p1++] = c[i]; else d[p0++] = c[i];
            }
            swap(cur, nxt);
        }
        S[lg].resize(n + 1); S[lg][0] = 0;
        for (int i = 0; i < n; i++) S[lg][i + 1] = S[lg][i] + cur[i];
        hasSum = true;
    }
    // [l,r) で key < x を満たす要素の個数。O(lg)。
    //   累積和配列 S に触れない軽量版 (buildSum していなくても使える)。
    //   x <= 0 なら 0、x >= 2^lg なら r-l を返す。
    ll cntLt(int l, int r, ll x) const {
        if (x <= 0) return 0;
        if (lg < 63 && x >= (1LL << lg)) return r - l;
        ll c = 0;
        for (int t = 0; t < lg; t++) {
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            if ((x >> (lg - 1 - t)) & 1) {
                c += r0 - l0;
                l = mid[t] + (l - l0); r = mid[t] + (r - r0);
            } else { l = l0; r = r0; }
            if (l >= r) break;
        }
        return c;
    }
    // レベル t0 から下だけを見て key < x となる個数。O(lg - t0)。
    //   x の上位ビット (レベル t0 より上) は無視され、下位ビットだけが使われる。
    //   cntRange の内部専用。
    ll cntLtFrom(int t0, int l, int r, ll x) const {
        ll c = 0;
        for (int t = t0; t < lg; t++) {
            if (l >= r) break;
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            if ((x >> (lg - 1 - t)) & 1) { c += r0 - l0; l = mid[t] + (l - l0); r = mid[t] + (r - r0); }
            else { l = l0; r = r0; }
        }
        return c;
    }
    // [l,r) で lo <= key < hi を満たす要素の個数。O(lg)。
    //   lo は 0 未満なら 0 に、hi は 2^lg より大きければ 2^lg に丸める。lo >= hi なら 0。
    //   lo と hi の上位ビットが一致している間は木の下降を 1 本にまとめている。
    ll cntRange(int l, int r, ll lo, ll hi) const {
        if (lo < 0) lo = 0;
        ll top = (lg < 63) ? (1LL << lg) : LLONG_MAX;
        if (hi > top) hi = top;
        if (lo >= hi) return 0;
        if (lo <= 0) return cntLt(l, r, hi);
        if (hi >= top) return (ll)(r - l) - cntLt(l, r, lo);
        int t = 0;
        while (t < lg && (((lo ^ hi) >> (lg - 1 - t)) & 1) == 0) {   // 共通接頭辞
            if (l >= r) return 0;
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            if ((lo >> (lg - 1 - t)) & 1) { l = mid[t] + (l - l0); r = mid[t] + (r - r0); }
            else { l = l0; r = r0; }
            t++;
        }
        if (l >= r) return 0;
        int l0, r0; rk0pair(lvl(t), l, r, l0, r0);          // ここで lo は bit 0 側, hi は bit 1 側
        int ol = mid[t] + (l - l0), orr = mid[t] + (r - r0);
        ll res = (ll)(r0 - l0) - cntLtFrom(t + 1, l0, r0, lo);
        res += cntLtFrom(t + 1, ol, orr, hi);
        return res;
    }
    // [l,r) で key < x を満たす要素の (個数, 重みの和) を 1 回の下降で返す。O(lg)。
    //   ★和が意味を持つのは buildSum 済みのときだけ (未構築なら和は 0 になる)。
    //   和は long long なのでオーバーフローには注意。
    pair<ll, ll> cntSumLt(int l, int r, ll x) const {
        if (x <= 0) return {0, 0};
        if (lg < 63 && x >= (1LL << lg)) return {r - l, hasSum ? S[0][r] - S[0][l] : 0};
        ll c = 0, s = 0;
        for (int t = 0; t < lg; t++) {
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            if ((x >> (lg - 1 - t)) & 1) {
                c += r0 - l0;
                if (hasSum) { const ll *p = S[t + 1].data(); s += p[r0] - p[l0]; }
                l = mid[t] + (l - l0); r = mid[t] + (r - r0);
            } else { l = l0; r = r0; }
            if (l >= r) break;
        }
        return {c, s};
    }
    // [l,r) で k 番目に小さい key。★k は 0-indexed で 0 <= k < r-l が前提。O(lg)。
    ll kth(int l, int r, ll k) const {
        ll res = 0;
        for (int t = 0; t < lg; t++) {
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            int c = r0 - l0;
            if (k < c) { l = l0; r = r0; }
            else { k -= c; res |= 1LL << (lg - 1 - t); l = mid[t] + (l - l0); r = mid[t] + (r - r0); }
        }
        return res;
    }
    // k 番目に小さい key (0-indexed) と、その key の [l,r) 内での出現個数を
    // 1 回の下降で同時に返す。★0 <= k < r-l が前提。O(lg)。
    void kthCount(int l, int r, ll k, ll &key, ll &cnt) const {
        ll res = 0;
        for (int t = 0; t < lg; t++) {
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            int c = r0 - l0;
            if (k < c) { l = l0; r = r0; }
            else { k -= c; res |= 1LL << (lg - 1 - t); l = mid[t] + (l - l0); r = mid[t] + (r - r0); }
        }
        key = res; cnt = r - l;
    }
    struct SmallRes { ll sum; ll rest; ll key; };
    // 小さい方から k 個の重み和を求めるための下降。★buildSum 済みであること。O(lg)。
    //   sum  : 「key が完全に k 個の中に収まったぶん」の重みの和
    //   rest : sum に含めきれなかった残りの個数 (すべて key と等しい値)
    //   key  : その境界にあたる key
    //   → 呼び出し側で sum + rest * (key に対応する値) とすれば答えになる。
    SmallRes smallestK(int l, int r, ll k) const {
        ll s = 0, res = 0;
        for (int t = 0; t < lg; t++) {
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            int c = r0 - l0;
            if (k < c) { l = l0; r = r0; }
            else {
                k -= c;
                const ll *p = S[t + 1].data(); s += p[r0] - p[l0];
                res |= 1LL << (lg - 1 - t);
                l = mid[t] + (l - l0); r = mid[t] + (r - r0);
            }
        }
        return SmallRes{s, k, res};
    }
    // 部分木 (レベル t, 範囲 [L,R), 確定済み上位ビット pre) に含まれる
    // 最大 key / 最小 key を返す。★[L,R) が空でないことが前提。O(lg - t)。
    ll maxIn(int t, int L, int R, ll pre) const {
        for (; t < lg; t++) {
            int l0, r0; rk0pair(lvl(t), L, R, l0, r0);
            if ((R - L) - (r0 - l0) > 0) { pre |= 1LL << (lg - 1 - t); L = mid[t] + (L - l0); R = mid[t] + (R - r0); }
            else { L = l0; R = r0; }
        }
        return pre;
    }
    ll minIn(int t, int L, int R, ll pre) const {
        for (; t < lg; t++) {
            int l0, r0; rk0pair(lvl(t), L, R, l0, r0);
            if (r0 > l0) { L = l0; R = r0; }
            else { pre |= 1LL << (lg - 1 - t); L = mid[t] + (L - l0); R = mid[t] + (R - r0); }
        }
        return pre;
    }
    // key = X の経路を 1 回だけ下降して
    //   pred = [l,r) 内で key <= X を満たす最大の key (X 自身を認めるかは includeEqual)
    //   succ = [l,r) 内で key >= X を満たす最小の key
    // を同時に得る。該当が無ければ -1 を返す。O(lg)。
    //   nearest のように両方必要なときに下降が 1 回で済む。
    //   ★0 <= X < 2^lg が前提 (呼び出し側で保証すること)。
    void predSucc(int l, int r, ll X, bool includeEqual, ll &pred, ll &succ) const {
        pred = -1; succ = -1;
        int cl = l, cr = r;
        int lL[64], lR[64], lT[64]; ll lP[64]; int lsp = 0;
        int rL[64], rR[64], rT[64]; ll rP[64]; int rsp = 0;
        for (int t = 0; t < lg && cl < cr; t++) {
            int l0, r0; rk0pair(lvl(t), cl, cr, l0, r0);
            int b = lg - 1 - t;
            ll pre = (X >> (b + 1)) << (b + 1);
            int ol = mid[t] + (cl - l0), orr = mid[t] + (cr - r0);
            if ((X >> b) & 1) {
                if (r0 > l0) { lL[lsp] = l0; lR[lsp] = r0; lP[lsp] = pre; lT[lsp] = t + 1; lsp++; }
                cl = ol; cr = orr;
            } else {
                if (orr > ol) { rL[rsp] = ol; rR[rsp] = orr; rP[rsp] = pre | (1LL << b); rT[rsp] = t + 1; rsp++; }
                cl = l0; cr = r0;
            }
        }
        bool has = (cl < cr);
        if (has && includeEqual) pred = X;
        else if (lsp > 0) { --lsp; pred = maxIn(lT[lsp], lL[lsp], lR[lsp], lP[lsp]); }
        if (has) succ = X;
        else if (rsp > 0) { --rsp; succ = minIn(rT[rsp], rL[rsp], rR[rsp], rP[rsp]); }
    }
    // [l,r) 内で key <= X を満たす最大の key。無ければ -1。O(lg)。
    //   X < 0 なら -1、X が 2^lg-1 より大きければ 2^lg-1 に丸める (= 区間最大値)。
    //   ※「個数を数えてから k 番目を取る」のではなく 1 回の下降で求める。
    ll prevKey(int l, int r, ll X) const {
        if (X < 0) return -1;
        ll top = (lg < 63) ? ((1LL << lg) - 1) : LLONG_MAX;
        if (X > top) X = top;
        int cl[64], cr[64], ct[64]; ll cp[64];
        int sp = 0;
        ll res = 0;
        for (int t = 0; t < lg; t++) {
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            int b = lg - 1 - t;
            if ((X >> b) & 1) {
                if (r0 > l0) { cl[sp] = l0; cr[sp] = r0; cp[sp] = res; ct[sp] = t + 1; sp++; }
                res |= 1LL << b;
                l = mid[t] + (l - l0); r = mid[t] + (r - r0);
            } else { l = l0; r = r0; }
            if (l >= r) break;
        }
        if (l < r) return res;
        if (sp == 0) return -1;
        sp--;
        return maxIn(ct[sp], cl[sp], cr[sp], cp[sp]);
    }
    // [l,r) 内で key >= X を満たす最小の key。無ければ -1。O(lg)。
    //   X < 0 は 0 に丸め、X >= 2^lg なら -1 を返す。
    //   ※「個数を数えてから k 番目を取る」のではなく 1 回の下降で求める。
    ll nextKey(int l, int r, ll X) const {
        if (X < 0) X = 0;
        if (lg < 63 && X >= (1LL << lg)) return -1;
        int cl[64], cr[64], ct[64]; ll cp[64];
        int sp = 0;
        ll res = 0;
        for (int t = 0; t < lg; t++) {
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            int b = lg - 1 - t;
            int ol = mid[t] + (l - l0), orr = mid[t] + (r - r0);
            if ((X >> b) & 1) {
                res |= 1LL << b;
                l = ol; r = orr;
            } else {
                if (orr > ol) { cl[sp] = ol; cr[sp] = orr; cp[sp] = res | (1LL << b); ct[sp] = t + 1; sp++; }
                l = l0; r = r0;
            }
            if (l >= r) break;
        }
        if (l < r) return res;
        if (sp == 0) return -1;
        sp--;
        return minIn(ct[sp], cl[sp], cr[sp], cp[sp]);
    }
    // 小さい方 (dfsSmall) / 大きい方 (dfsLarge) から need 個ぶんを
    // (key, 個数) の形で out に積む。dfsSmall は key 昇順、dfsLarge は key 降順。
    //   need は参照で減っていき、0 になった時点で打ち切る。
    //   ★呼び出し時は t = 0, pre = 0, out は空にしておくこと。
    //   計算量は訪れるノード数に比例し、O((取り出す相異なる値の個数) * lg)。
    void dfsSmall(int t, int l, int r, ll pre, ll &need, vector<pair<ll, int>> &out) const {
        if (l >= r || need <= 0) return;
        if (t == lg) { int c = (int)std::min<ll>(need, (ll)(r - l)); out.push_back({pre, c}); need -= c; return; }
        int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
        dfsSmall(t + 1, l0, r0, pre, need, out);
        if (need > 0) dfsSmall(t + 1, mid[t] + (l - l0), mid[t] + (r - r0), pre | (1LL << (lg - 1 - t)), need, out);
    }
    void dfsLarge(int t, int l, int r, ll pre, ll &need, vector<pair<ll, int>> &out) const {
        if (l >= r || need <= 0) return;
        if (t == lg) { int c = (int)std::min<ll>(need, (ll)(r - l)); out.push_back({pre, c}); need -= c; return; }
        int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
        dfsLarge(t + 1, mid[t] + (l - l0), mid[t] + (r - r0), pre | (1LL << (lg - 1 - t)), need, out);
        if (need > 0) dfsLarge(t + 1, l0, r0, pre, need, out);
    }
    // [l,r) で (key ^ x) < bound を満たす個数。O(lg)。
    //   ★x, bound はともに 0 <= 値 < 2^lg に正規化済みであることが前提
    //     (StaticSeq::count_xor_lt が上位ビットを処理してから呼んでいる)。
    ll cntXorLt(int l, int r, ull x, ull bound) const {
        if (bound == 0) return 0;
        if (lg < 63 && bound >= (1ULL << lg)) return r - l;
        ll c = 0;
        for (int t = 0; t < lg; t++) {
            int b = lg - 1 - t;
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            int c0 = r0 - l0, c1 = (r - l) - c0;
            int xb = (int)((x >> b) & 1);
            int czero = xb ? c1 : c0;                  // (key^x) の bit b が 0 の個数
            int ol = mid[t] + (l - l0), orr = mid[t] + (r - r0);
            if ((bound >> b) & 1) {
                c += czero;
                if (xb) { l = l0; r = r0; }            // 進む先は key bit = 1-xb
                else    { l = ol; r = orr; }
            } else {
                if (xb) { l = ol; r = orr; }
                else    { l = l0; r = r0; }
            }
            if (l >= r) break;
        }
        return c;
    }
    // [l,r) の要素について (key ^ x) を昇順に並べたときの k 番目 (0-indexed) の値。O(lg)。
    //   ★0 <= k < r-l、かつ 0 <= x < 2^lg が前提。
    ull kthXor(int l, int r, ull x, ll k) const {
        ull res = 0;
        for (int t = 0; t < lg; t++) {
            int b = lg - 1 - t;
            int l0, r0; rk0pair(lvl(t), l, r, l0, r0);
            int c0 = r0 - l0, c1 = (r - l) - c0;
            int xb = (int)((x >> b) & 1);
            int czero = xb ? c1 : c0;
            int ol = mid[t] + (l - l0), orr = mid[t] + (r - r0);
            if (k < czero) {
                if (xb) { l = ol; r = orr; }
                else    { l = l0; r = r0; }
            } else {
                k -= czero; res |= 1ULL << b;
                if (xb) { l = l0; r = r0; }
                else    { l = ol; r = orr; }
            }
        }
        return res;
    }
};

// ============================================================================
//  永続セグメント木 (区間 mex 専用)
//    値 v (0 <= v < m) をキーとして「その値が最後に現れた位置」を持ち、
//    区間 min を取れるようにしたもの。位置 r まで挿入した状態を version r とする。
//    push  : 1 回あたり O(log m) ノードを新規作成 (時間も O(log m))
//    findFirstLess : O(log m)
//    ★init の reserveNodes は「作られるノード数の見積もり」。足りなくても
//      動作は正しいが再確保が起きる。ノードは push 1 回につき O(log m) 個増える。
//    左右の子と min を 1 つの構造体にまとめ、下降 1 段あたりの読み込みを減らしている。
// ============================================================================
struct PersistentMinSeg {
    struct PNode { int lc, rc, mn; };          // 左右の子の添字と、部分木の min
    int m = 0;
    vector<PNode> t;
    vector<int> root;
    // 値域 [0, m_) で初期化する。version 0 は「まだ何も挿入していない状態」。
    void init(int m_, int reserveNodes) {
        m = m_;
        t.clear(); t.reserve(reserveNodes);
        t.push_back(PNode{0, 0, -1});             // node 0 : 全体が -1 (未出現) を表す番人
        root.clear(); root.reserve(m_ + 2); root.push_back(0);
    }
    // 直前 version から「値 p の最終出現位置を val にした」新しい部分木を作る (内部用)
    int upd(int node, int nl, int nr, int p, int val) {
        int cur = (int)t.size();
        PNode nd = t[node];
        t.push_back(nd);
        if (nr - nl == 1) { t[cur].mn = val; return cur; }
        int nm = (nl + nr) >> 1;
        if (p < nm) { int c = upd(nd.lc, nl, nm, p, val); t[cur].lc = c; t[cur].mn = min(t[c].mn, t[nd.rc].mn); }
        else        { int c = upd(nd.rc, nm, nr, p, val); t[cur].rc = c; t[cur].mn = min(t[nd.lc].mn, t[c].mn); }
        return cur;
    }
    // 最新 version に「値 p が位置 val に現れた」ことを反映した新しい version を追加する。
    void push(int p, int val) { root.push_back(upd(root.back(), 0, m, p, val)); }
    // version ver において「最後の出現位置 < l」を満たす最小の値を返す。O(log m)。
    //   = [l, ver) に含まれない最小の値 (区間 mex そのもの)。
    int findFirstLess(int ver, int l) const {
        const PNode *p = t.data();
        int node = root[ver], nl = 0, nr = m;
        while (nr - nl > 1) {
            int nm = (nl + nr) >> 1;
            const PNode &nd = p[node];
            if (p[nd.lc].mn < l) { node = nd.lc; nr = nm; }
            else                 { node = nd.rc; nl = nm; }
        }
        return nl;
    }
};

// ============================================================================
//  区間の最大部分和 / 最小部分和 用の静的セグメント木
//    ★累積和配列 P[0..N] に対して構築する (元の配列に対してではない)。
//    区間 [l, r] の P について best = max(P[j]-P[i]) (i<j) を持たせることで、
//    元配列の [l, r) 内の連続部分列の和の最大値 / 最小値が求まる。
//    build : O(N) 時間 / ノード数は 2 * (N 以上の 2 冪) 個
//    query : O(log N)   ※マージが可換でないので左右を分けて集約している
//    ★番人に LLONG_MAX/4, LLONG_MIN/4 を使っているので、
//      累積和の絶対値がこの程度を超えるほど大きい入力では使えない。
// ============================================================================
struct SubSegTree {
    struct Node { ll mn, mx, best, worst; };    // best = max(P[j]-P[i]) (i<j), worst = min(...)
    int sz = 1;
    vector<Node> t;
    static inline Node merge(const Node &a, const Node &b) {
        Node c;
        c.mn = min(a.mn, b.mn);
        c.mx = max(a.mx, b.mx);
        c.best  = max(max(a.best, b.best), b.mx - a.mn);
        c.worst = min(min(a.worst, b.worst), b.mn - a.mx);
        return c;
    }
    // p は累積和配列 (長さ N+1) を想定。O(N)。
    void build(const vector<ll> &p) {
        int m = (int)p.size();
        sz = 1; while (sz < m) sz <<= 1;
        t.assign(2 * sz, Node{LLONG_MAX / 4, LLONG_MIN / 4, LLONG_MIN / 4, LLONG_MAX / 4});
        for (int i = 0; i < m; i++) t[sz + i] = Node{p[i], p[i], LLONG_MIN / 4, LLONG_MAX / 4};
        for (int i = sz - 1; i >= 1; i--) t[i] = merge(t[2 * i], t[2 * i + 1]);
    }
    // ★累積和配列の添字での [l, r) を集約する。l < r が前提。O(log N)。
    Node query(int l, int r) const {
        Node resL{LLONG_MAX / 4, LLONG_MIN / 4, LLONG_MIN / 4, LLONG_MAX / 4};
        Node resR = resL;
        int lo = l + sz, hi = r + sz;
        while (lo < hi) {
            if (lo & 1) { resL = merge(resL, t[lo]); lo++; }
            if (hi & 1) { --hi; resR = merge(t[hi], resR); }
            lo >>= 1; hi >>= 1;
        }
        return merge(resL, resR);
    }
};

// ############################################################################
// #                              StaticSeq                                   #
// #   index は 0..n-1、区間は半開区間 [l, r)。                                 #
// #   コンストラクタで作るのは下の 3 本の配列だけ (O(N))。                     #
// #   それ以外の内部構造は、対応する関数を最初に呼んだときに構築される。       #
// ############################################################################
struct StaticSeq {
    // ------------------------------------------------------------------ 基本
    int n = 0;
    vector<ll> a;           // 元の配列 (構築後は変更しない)
    vector<ll> ps;          // 累積和   ps[i] = a[0]+...+a[i-1]      (長さ n+1)
    vector<ll> pxr;         // 累積 XOR pxr[i] = a[0]^...^a[i-1]     (長さ n+1)
    bool allNonNeg = true;  // 全要素が 0 以上か (XOR 系 / 線形基底の前提チェック用)

    StaticSeq() {}
    explicit StaticSeq(vector<ll> v) { init(std::move(v)); }

    // 配列を受け取り直す (既に構築済みの内部構造は破棄されない点に注意)。O(N)。
    void init(vector<ll> v) {
        a = std::move(v);
        n = (int)a.size();
        ps.assign(n + 1, 0);
        pxr.assign(n + 1, 0);
        allNonNeg = true;
        for (int i = 0; i < n; i++) {
            ps[i + 1]  = ps[i] + a[i];
            pxr[i + 1] = pxr[i] ^ a[i];
            if (a[i] < 0) allNonNeg = false;
        }
    }

    inline int size()  const { return n; }                           // 要素数            O(1)
    inline bool empty() const { return n == 0; }                     // 空配列か          O(1)
    inline ll get(int i) const { return a[i]; }                      // a[i] (0<=i<n)     O(1)
    inline ll operator[](int i) const { return a[i]; }               // a[i] と同じ       O(1)
    inline ll front() const { return a[0]; }                         // a[0]              O(1)
    inline ll back()  const { return a[n - 1]; }                     // a[n-1]            O(1)
    // 部分配列 a[l..r-1] のコピーを返す。★コピーするので O(r-l)。
    vector<ll> get(int l, int r) const { return vector<ll>(a.begin() + l, a.begin() + r); }
    const vector<ll>& raw() const { return a; }                      // 元配列への参照    O(1)
    // 階差 a[i+1]-a[i]。★0 <= i < n-1 が前提。                                   O(1)
    inline ll diff(int i) const { return a[i + 1] - a[i]; }

    // ================================================================ 累積和系
    //   この 3 つは prefix 配列の差なので、l == r (空区間) でも 0 を返して安全。
    //   0 <= l <= r <= n を満たしていればよい。戻り値は long long。
    inline ll prefix_sum(int r) const { return ps[r]; }              // a[0..r-1] の総和  O(1)
    inline ll sum(int l, int r) const { return ps[r] - ps[l]; }      // a[l..r-1] の総和  O(1)
    inline ll bit_xor(int l, int r) const { return pxr[r] ^ pxr[l]; }// a[l..r-1] の XOR  O(1)

    // ---- 高次モーメント (初回呼び出し時に O(N) で構築, __int128 で保持) ----
    //   sum128 / sum_sq / sum_ia / sum_i2a / sum_sq_distance と、
    //   それらを使う統計量 (mean / variance / fit_line など) が使う。
    //   __int128 の配列を 4 本持つので、メモリは 16 * 4 * (N+1) バイト。
    bool b_mom = false;
    vector<i128> m_s, m_sq, m_ia, m_i2a;
    void build_moment() {
        if (b_mom) return;
        b_mom = true;
        m_s.assign(n + 1, 0); m_sq.assign(n + 1, 0);
        m_ia.assign(n + 1, 0); m_i2a.assign(n + 1, 0);
        const ll *p = a.data();
        for (int i = 0; i < n; i++) {
            ll x = p[i], ii = i;                 // (i128)ll * ll は 64x64->128 の 1 命令になる
            m_s[i + 1]   = m_s[i]   + x;
            m_sq[i + 1]  = m_sq[i]  + (i128)x * x;
            m_ia[i + 1]  = m_ia[i]  + (i128)ii * x;
            m_i2a[i + 1] = m_i2a[i] + (i128)(ii * ii) * x;
        }
    }
    // すべて [l,r) について、__int128 で返す (long long の桁あふれを避けたいとき用)。
    // 0 <= l <= r <= n であればよく、l == r なら 0。いずれも O(1)。
    i128 sum128(int l, int r) { build_moment(); return m_s[r]  - m_s[l];  }   // Σ a_i
    i128 sum_sq(int l, int r) { build_moment(); return m_sq[r] - m_sq[l]; }   // Σ a_i^2
    i128 sum_ia(int l, int r) { build_moment(); return m_ia[r] - m_ia[l]; }   // Σ i * a_i
    i128 sum_i2a(int l, int r){ build_moment(); return m_i2a[r]- m_i2a[l];}   // Σ i^2 * a_i
    // Σ (a_i - x)^2  (x との二乗距離の総和)                                   O(1)
    i128 sum_sq_distance(int l, int r, ll x) {
        build_moment();
        i128 s = sum128(l, r), q = sum_sq(l, r);
        return q - (i128)2 * x * s + (i128)(r - l) * x * x;
    }

    // ---- 符号 / 偶奇 / popcount 別の累積 (初回呼び出し時に O(N) で構築) ----
    //   count_positive / count_negative / count_zero / count_odd /
    //   sum_positive / sum_negative / sum_abs(l,r) / sum_popcount が使う。
    bool b_sign = false;
    vector<int> c_pos, c_neg, c_zero, c_odd;
    vector<ll>  s_pos, s_neg, s_abs, s_pop;
    void build_sign() {
        if (b_sign) return;
        b_sign = true;
        c_pos.assign(n + 1, 0); c_neg.assign(n + 1, 0); c_zero.assign(n + 1, 0);
        c_odd.assign(n + 1, 0);
        s_pos.assign(n + 1, 0); s_neg.assign(n + 1, 0);
        s_abs.assign(n + 1, 0); s_pop.assign(n + 1, 0);
        const ll *pa = a.data();
        for (int i = 0; i < n; i++) {
            ll x = pa[i];
            c_pos[i + 1]  = c_pos[i]  + (x > 0);
            c_neg[i + 1]  = c_neg[i]  + (x < 0);
            c_zero[i + 1] = c_zero[i] + (x == 0);
            c_odd[i + 1]  = c_odd[i]  + (int)(((ull)x) & 1ULL);
            s_pos[i + 1]  = s_pos[i]  + (x > 0 ? x : 0);
            s_neg[i + 1]  = s_neg[i]  + (x < 0 ? x : 0);
            s_abs[i + 1]  = s_abs[i]  + (x < 0 ? -x : x);
            s_pop[i + 1]  = s_pop[i]  + __builtin_popcountll((ull)x);
        }
    }
    // いずれも [l,r) について O(1)。0 <= l <= r <= n であればよい (l == r なら 0)。
    ll count_positive(int l, int r) { build_sign(); return c_pos[r]  - c_pos[l];  }   // a_i > 0 の個数
    ll count_negative(int l, int r) { build_sign(); return c_neg[r]  - c_neg[l];  }   // a_i < 0 の個数
    ll count_zero(int l, int r)     { build_sign(); return c_zero[r] - c_zero[l]; }   // a_i == 0 の個数
    ll count_nonzero(int l, int r)  { return (r - l) - count_zero(l, r); }            // a_i != 0 の個数
    // 偶奇は 2 の補数表現の最下位ビットで判定するので、負数もそのまま扱える
    ll count_odd(int l, int r)      { build_sign(); return c_odd[r]  - c_odd[l];  }   // 奇数の個数
    ll count_even(int l, int r)     { return (r - l) - count_odd(l, r); }             // 偶数の個数
    ll sum_positive(int l, int r)   { build_sign(); return s_pos[r]  - s_pos[l];  }   // 正の要素だけの和
    ll sum_negative(int l, int r)   { build_sign(); return s_neg[r]  - s_neg[l];  }   // 負の要素だけの和 (<=0)
    ll sum_abs(int l, int r)        { build_sign(); return s_abs[r]  - s_abs[l];  }   // Σ|a_i|
    // popcount は 2 の補数表現に対して取るので、負数は 1 が多く数えられる
    ll sum_popcount(int l, int r)   { build_sign(); return s_pop[r]  - s_pop[l];  }   // Σ popcount(a_i)

    // ---- bit ごとの 1 の個数 (初回呼び出し時に構築) ----
    //   全要素の OR を取り、1 度でも立っているビットについてだけ累積和を持つので、
    //   構築は O(N * (使われているビット数))、メモリも同じオーダー。
    //   負数があると上位ビットも「使われている」と判定される点に注意。
    bool b_bitcnt = false;
    vector<int> bitSlot;                 // bit -> bitPre の行番号 (-1 なら全要素でそのビットが 0)
    vector<vector<int>> bitPre;          // bitPre[行][i] = [0,i) でそのビットが立つ個数
    int bitRows = 0;                     // 行数 (現状ほかから参照していない)
    void build_bitcnt() {
        if (b_bitcnt) return;
        b_bitcnt = true;
        ull orAll = 0;
        for (int i = 0; i < n; i++) orAll |= (ull)a[i];
        bitSlot.assign(64, -1);
        int cnt = 0;
        for (int b = 0; b < 64; b++) if ((orAll >> b) & 1) bitSlot[b] = cnt++;
        bitRows = cnt;
        bitPre.assign(cnt, vector<int>());
        const ll *pa = a.data();
        for (int b = 0; b < 64; b++) {
            int s = bitSlot[b]; if (s < 0) continue;
            // 1 行ずつ「確保 → すぐ書き込み」でキャッシュに載っているうちに埋める
            bitPre[s].resize(n + 1);
            int *v = bitPre[s].data();
            v[0] = 0;
            for (int i = 0; i < n; i++) v[i + 1] = v[i] + (int)(((ull)pa[i] >> b) & 1ULL);
        }
    }
    // [l,r) で bit b (0 <= b < 64) が立っている要素数。O(1)。
    //   全要素でそのビットが 0 なら即 0 を返す。0 <= l <= r <= n であればよい。
    ll count_bit(int l, int r, int b) {
        build_bitcnt();
        int s = bitSlot[b]; if (s < 0) return 0;
        const int *v = bitPre[s].data();
        return v[r] - v[l];
    }

    // ================================================================ 統計量
    //   すべて上のモーメント配列から O(1) で求める。戻り値は long double。
    //   ★(r-l) で割るので l < r が前提 (空区間では 0 除算になる)。
    long double mean(int l, int r) { return (long double)sum128(l, r) / (long double)(r - l); }
    // 母分散 (n-1 ではなく n で割る)。E[X^2]-E[X]^2 で計算するため、
    // 値が大きく分散が小さい入力では桁落ちしうる。
    long double variance(int l, int r) {
        long double m = mean(l, r);
        return (long double)sum_sq(l, r) / (long double)(r - l) - m * m;
    }
    // 標準偏差。桁落ちで variance が 0 以下になった場合は 0 を返す。
    long double stddev(int l, int r) { long double v = variance(l, r); return v > 0 ? sqrtl(v) : 0.0L; }
    // index そのものの和 : idx_sum = Σ_{i=l}^{r-1} i , idx_sum2 = Σ_{i=l}^{r-1} i^2
    // 閉じた式で求める内部用ヘルパ。O(1)。
    static i128 idx_sum(ll l, ll r)  { return (i128)(l + r - 1) * (r - l) / 2; }
    static i128 idx_sum2(ll l, ll r) {
        auto f = [](i128 m) -> i128 { return m < 0 ? (i128)0 : m * (m + 1) * (2 * m + 1) / 6; };
        return f((i128)r - 1) - f((i128)l - 1);
    }
    // index i と値 a_i を 2 変数と見なしたときの共分散。★l < r が前提。O(1)。
    long double covariance_index(int l, int r) {
        long double cnt = (long double)(r - l);
        long double si  = (long double)idx_sum(l, r);
        long double sa  = (long double)sum128(l, r);
        long double sia = (long double)sum_ia(l, r);
        return sia / cnt - (si / cnt) * (sa / cnt);
    }
    // index i 自身の母分散 (値には依存しない)。★l < r が前提。O(1)。
    long double variance_index(int l, int r) {
        long double cnt = (long double)(r - l);
        long double si  = (long double)idx_sum(l, r);
        long double sii = (long double)idx_sum2(l, r);
        return sii / cnt - (si / cnt) * (si / cnt);
    }
    // 最小二乗直線 y = slope * i + intercept の (slope, intercept) を返す。O(1)。
    //   ★l < r が前提。r-l == 1 のように index の分散が 0 のときは slope = 0 とする。
    pair<long double, long double> fit_line(int l, int r) {
        long double vi = variance_index(l, r);
        long double slope = (vi > 0) ? covariance_index(l, r) / vi : 0.0L;
        long double mi = (long double)idx_sum(l, r) / (long double)(r - l);
        return {slope, mean(l, r) - slope * mi};
    }
    // index と値の相関係数。★l < r が前提。O(1)。
    //   どちらかの分散が 0 以下 (定数列など) のときは 0 を返す。
    long double correlation_index(int l, int r) {
        long double vi = variance_index(l, r), va = variance(l, r);
        if (vi <= 0 || va <= 0) return 0.0L;
        return covariance_index(l, r) / sqrtl(vi * va);
    }

    // ============================================================ RMQ (min/max)
    //   初回呼び出し時に O(N) で構築、以降クエリは O(1)。min 用と max 用は別々に構築。
    //   ★区間は [l,r) で l < r が前提 (空区間は不可)。
    //   ★最小値 / 最大値が複数あるときは「最も左の位置」を返す。
    bool b_rmin = false, b_rmax = false;
    RMQ64<ll, false> rmin;
    RMQ64<ll, true>  rmax;
    void build_rmin() { if (b_rmin) return;
 b_rmin = true; rmin.build(a.data(), n); }
    void build_rmax() { if (b_rmax) return;
 b_rmax = true; rmax.build(a.data(), n); }

    int  argmin(int l, int r) { build_rmin(); return rmin.query(a.data(), l, r); }   // 最小値の位置(最左)
    int  argmax(int l, int r) { build_rmax(); return rmax.query(a.data(), l, r); }   // 最大値の位置(最左)
    ll   min(int l, int r)    { return a[argmin(l, r)]; }                  // 最小値
    ll   max(int l, int r)    { return a[argmax(l, r)]; }                  // 最大値
    pair<ll, ll> minmax(int l, int r) { return {min(l, r), max(l, r)}; }             // (最小値, 最大値)
    pair<int, int> argminmax(int l, int r) { return {argmin(l, r), argmax(l, r)}; }  // (最小位置, 最大位置)
    pair<ll, int> min_with_index(int l, int r) { int i = argmin(l, r); return {a[i], i}; }  // (値, 位置)
    pair<ll, int> max_with_index(int l, int r) { int i = argmax(l, r); return {a[i], i}; }  // (値, 位置)
    ll   range_width(int l, int r) { return max(l, r) - min(l, r); }       // max - min

    // ============================================== prefix / suffix の集約
    //   初回呼び出し時に O(N) で 14 本の配列をまとめて構築、以降クエリは O(1)。
    //   ★prefix_xxx(i) は a[0..i]、suffix_xxx(i) は a[i..n-1] (どちらも i を含む)。
    //   ★0 <= i < n が前提。
    //   ★同値のときの位置は prefix 側が最左、suffix 側も最左 (= 最も小さい index)。
    //   ★gcd は絶対値の gcd なので常に 0 以上。and / or は元の値のまま。
    bool b_ps = false;
    vector<ll> pmin, pmax, smin, smax, pgcd, sgcd, pand, sand_, por, sor;
    vector<int> pamin, pamax, samin, samax;
    void build_prefix_agg() {
        if (b_ps) return;
        b_ps = true;
        pmin.resize(n); pmax.resize(n); smin.resize(n); smax.resize(n);
        pamin.resize(n); pamax.resize(n); samin.resize(n); samax.resize(n);
        pgcd.resize(n); sgcd.resize(n); pand.resize(n); sand_.resize(n);
        por.resize(n);  sor.resize(n);
        for (int i = 0; i < n; i++) {
            if (i == 0) {
                pmin[i] = pmax[i] = a[i]; pamin[i] = pamax[i] = 0;
                pgcd[i] = a[i] < 0 ? -a[i] : a[i]; pand[i] = a[i]; por[i] = a[i];
            } else {
                if (a[i] < pmin[i - 1]) { pmin[i] = a[i]; pamin[i] = i; }
                else { pmin[i] = pmin[i - 1]; pamin[i] = pamin[i - 1]; }
                if (a[i] > pmax[i - 1]) { pmax[i] = a[i]; pamax[i] = i; }
                else { pmax[i] = pmax[i - 1]; pamax[i] = pamax[i - 1]; }
                pgcd[i] = ll_gcd(pgcd[i - 1], a[i]);
                pand[i] = pand[i - 1] & a[i];
                por[i]  = por[i - 1] | a[i];
            }
        }
        for (int i = n - 1; i >= 0; i--) {
            if (i == n - 1) {
                smin[i] = smax[i] = a[i]; samin[i] = samax[i] = i;
                sgcd[i] = a[i] < 0 ? -a[i] : a[i]; sand_[i] = a[i]; sor[i] = a[i];
            } else {
                if (a[i] <= smin[i + 1]) { smin[i] = a[i]; samin[i] = i; }
                else { smin[i] = smin[i + 1]; samin[i] = samin[i + 1]; }
                if (a[i] >= smax[i + 1]) { smax[i] = a[i]; samax[i] = i; }
                else { smax[i] = smax[i + 1]; samax[i] = samax[i + 1]; }
                sgcd[i] = ll_gcd(sgcd[i + 1], a[i]);
                sand_[i] = sand_[i + 1] & a[i];
                sor[i]  = sor[i + 1] | a[i];
            }
        }
    }
    ll prefix_min(int i) { build_prefix_agg(); return pmin[i]; }   // a[0..i] の最小値
    ll prefix_max(int i) { build_prefix_agg(); return pmax[i]; }   // a[0..i] の最大値
    ll suffix_min(int i) { build_prefix_agg(); return smin[i]; }   // a[i..n-1] の最小値
    ll suffix_max(int i) { build_prefix_agg(); return smax[i]; }   // a[i..n-1] の最大値
    int prefix_argmin(int i) { build_prefix_agg(); return pamin[i]; }   // a[0..i] の最小位置
    int prefix_argmax(int i) { build_prefix_agg(); return pamax[i]; }   // a[0..i] の最大位置
    int suffix_argmin(int i) { build_prefix_agg(); return samin[i]; }   // a[i..n-1] の最小位置
    int suffix_argmax(int i) { build_prefix_agg(); return samax[i]; }   // a[i..n-1] の最大位置
    ll prefix_gcd(int i) { build_prefix_agg(); return pgcd[i]; }   // a[0..i] の gcd (絶対値)
    ll suffix_gcd(int i) { build_prefix_agg(); return sgcd[i]; }   // a[i..n-1] の gcd (絶対値)
    ll prefix_and(int i) { build_prefix_agg(); return pand[i]; }   // a[0..i] の AND
    ll suffix_and(int i) { build_prefix_agg(); return sand_[i]; }  // a[i..n-1] の AND
    ll prefix_or(int i)  { build_prefix_agg(); return por[i]; }    // a[0..i] の OR
    ll suffix_or(int i)  { build_prefix_agg(); return sor[i]; }    // a[i..n-1] の OR

    // =============================== 区間 GCD / AND / OR  (Sparse Table)
    //   初回呼び出し時に O(N log N) で構築 (演算ごとに別のテーブル)、クエリは O(1)。
    //   ★[l,r) で l < r が前提。
    //   ★gcd は |a_i| の gcd を返す (常に 0 以上)。全要素 0 の区間なら 0。
    //   ★メモリは 1 テーブルあたり long long を N * ceil(log2 N) 個ぶん使うので、
    //     3 種すべて呼ぶとかなり大きくなる。XOR は prefix XOR の bit_xor を使うこと。
    bool b_gcd = false, b_and = false, b_or = false;
    SparseOp<OpGcd> stGcd;
    SparseOp<OpAnd> stAnd;
    SparseOp<OpOr>  stOr;
    ll gcd(int l, int r) {
        if (!b_gcd) {
            b_gcd = true;
            ll *p = stGcd.level0(n);
            for (int i = 0; i < n; i++) p[i] = a[i] < 0 ? -a[i] : a[i];
            stGcd.finish();
        }
        return stGcd.query(l, r);
    }
    ll bit_and(int l, int r) { if (!b_and) { b_and = true; stAnd.build(a); } return stAnd.query(l, r); }
    ll bit_or (int l, int r) { if (!b_or)  { b_or  = true; stOr .build(a); } return stOr .query(l, r); }

    // ================================ 座標圧縮 + 値ごとの出現位置リスト (CSR)
    //   初回呼び出し時に O(N log N) で構築。ここで作る 5 本の配列が、
    //   Wavelet Matrix / 出現位置系 / 区間 mex / 区間最頻値 / Suffix Array の土台になる。
    //   K = 相異なる値の個数。
    bool b_cmp = false;
    vector<ll>  vals;        // 昇順ソート済みのユニーク値 (長さ K)
    vector<int> rnk;         // rnk[i] = a[i] の順位 (0..K-1)。vals[rnk[i]] == a[i]
    vector<int> occSt;       // CSR : 値 v の出現位置は occPos[occSt[v] .. occSt[v+1]) (長さ K+1)
    vector<int> occPos;      // 値ごとにまとめた出現位置。同じ値の中では昇順
    vector<int> occIdx;      // occIdx[i] = i が a[i] の何番目 (0-indexed) の出現か
    void build_compress() {
        if (b_cmp) return;
        b_cmp = true;
        // (値, 位置) を 1 回ソートするだけで vals / rnk / occPos / occIdx / occSt が同時に得られる
        // (元実装の「ソート + N 回の lower_bound + カウントソート」より高速)
        vector<pair<ll, int>> pv(n);
        for (int i = 0; i < n; i++) pv[i] = {a[i], i};
        sort(pv.begin(), pv.end());
        vals.clear(); vals.reserve(n);
        occSt.clear(); occSt.reserve(n + 1);
        rnk.resize(n); occPos.resize(n); occIdx.resize(n);
        int cur = -1, gs = 0;
        for (int j = 0; j < n; j++) {
            if (j == 0 || pv[j].first != pv[j - 1].first) { vals.push_back(pv[j].first); occSt.push_back(j); cur++; gs = j; }
            int i = pv[j].second;
            rnk[i] = cur; occPos[j] = i; occIdx[i] = j - gs;
        }
        occSt.push_back(n);
        vals.shrink_to_fit();
    }
    int unique_size() { build_compress(); return (int)vals.size(); }      // 相異なる値の個数 K   O(1)
    const vector<ll>& unique_values() { build_compress(); return vals; }  // 昇順ユニーク値への参照 O(1)
    int compressed(int i) { build_compress(); return rnk[i]; }            // a[i] の順位 (0<=i<n)  O(1)
    ll  decompress(int v) { build_compress(); return vals[v]; }           // 順位 -> 値 (0<=v<K)   O(1)
    // rank_lb(x) : x 未満のユニーク値の個数 (= x 以上で最小のユニーク値の順位)。O(log K)
    // rank_ub(x) : x 以下のユニーク値の個数 (= x より大きい最小のユニーク値の順位)。O(log K)
    //   どちらも x が配列に無くても使える (0..K の値を返す)。
    int rank_lb(ll x) { build_compress(); return lb_idx(vals.data(), (int)vals.size(), x); }
    int rank_ub(ll x) { build_compress(); return ub_idx(vals.data(), (int)vals.size(), x); }
    int value_index(ll x) {                     // x が存在すればその順位, 無ければ -1  O(log K)
        build_compress();
        int p = lb_idx(vals.data(), (int)vals.size(), x);
        return (p < (int)vals.size() && vals[p] == x) ? p : -1;
    }
    // --- 内部用 : 順位 v (0<=v<K) をそのまま受け取る出現リスト操作。O(log(その値の出現数))。
    //     値から順位を引く二分探索を省けるところで使う。build_compress 済みが前提。 ---
    // [l,r) に含まれる「順位 v の値」の個数
    inline ll count_rank(int l, int r, int v) const {
        const int *b = occPos.data() + occSt[v];
        int m = occSt[v + 1] - occSt[v];
        return (ll)(lb_idx(b, m, r) - lb_idx(b, m, l));
    }
    // [l,r) で「順位 v の値」が最初に現れる位置 (無ければ -1)
    inline int first_eq_rank(int l, int r, int v) const {
        const int *b = occPos.data() + occSt[v];
        int m = occSt[v + 1] - occSt[v];
        int p = lb_idx(b, m, l);
        return (p < m && b[p] < r) ? b[p] : -1;
    }
    // [l,r) で「順位 v の値」が最後に現れる位置 (無ければ -1)
    inline int last_eq_rank(int l, int r, int v) const {
        const int *b = occPos.data() + occSt[v];
        int m = occSt[v + 1] - occSt[v];
        int p = lb_idx(b, m, r);
        return (p > 0 && b[p - 1] >= l) ? b[p - 1] : -1;
    }
    // 値 x の出現位置すべてを昇順で返す。無ければ空。O(log K + 出現回数)。
    //   ★毎回 vector を作って返すので、ループ内で何度も呼ばないこと。
    vector<int> positions(ll x) {
        build_compress();
        int v = value_index(x); if (v < 0) return {};
        return vector<int>(occPos.begin() + occSt[v], occPos.begin() + occSt[v + 1]);
    }
    ll count(ll x) {                            // 配列全体での値 x の出現回数  O(log K)
        build_compress();
        int v = value_index(x); return v < 0 ? 0 : occSt[v + 1] - occSt[v];
    }
    // 配列全体で値 x が最初 / 最後に現れる位置 (無ければ -1)。O(log K)。
    int first_pos(ll x) { build_compress(); int v = value_index(x); return v < 0 ? -1 : occPos[occSt[v]]; }
    int last_pos(ll x)  { build_compress(); int v = value_index(x); return v < 0 ? -1 : occPos[occSt[v + 1] - 1]; }
    // [l,r) に含まれる値 x の個数                                    O(log K + log(出現数))
    ll count(int l, int r, ll x) {
        build_compress();
        int v = value_index(x); if (v < 0) return 0;
        return count_rank(l, r, v);
    }
    // [l,r) に値 x が 1 つでもあるか。個数は数えず最初の出現だけ調べる。
    //                                                              O(log K + log(出現数))
    bool exists(int l, int r, ll x) {
        build_compress();
        int v = value_index(x); if (v < 0) return false;
        return first_eq_rank(l, r, v) >= 0;
    }
    // [l,r) で最初に x が現れる位置 (無ければ -1)                    O(log K + log(出現数))
    int first_equal(int l, int r, ll x) {
        build_compress();
        int v = value_index(x); if (v < 0) return -1;
        return first_eq_rank(l, r, v);
    }
    // [l,r) で最後に x が現れる位置 (無ければ -1)                    O(log K + log(出現数))
    int last_equal(int l, int r, ll x) {
        build_compress();
        int v = value_index(x); if (v < 0) return -1;
        return last_eq_rank(l, r, v);
    }
    // 配列全体をソートしたときの k 番目 (0-indexed) の値。★0 <= k < n が前提。O(log K)。
    ll kth_global(ll k) {
        build_compress();
        int lo = 0, hi = (int)vals.size() - 1;
        while (lo < hi) { int m = (lo + hi) / 2; if ((ll)occSt[m + 1] > k) hi = m; else lo = m + 1; }
        return vals[lo];
    }
    // 配列全体の中央値。要素数が偶数のときは小さい方を返す。★n >= 1 が前提。O(log K)。
    ll median_global() { return kth_global((n - 1) / 2); }
    // 配列全体で x 以上の最小の要素 (無ければ LLONG_MAX)                        O(log K)
    ll lower_bound_global(ll x) {
        build_compress();
        int p = rank_lb(x); return p < (int)vals.size() ? vals[p] : LLONG_MAX;
    }
    // 配列全体で x より大きい最小の要素 (無ければ LLONG_MAX)                   O(log K)
    ll upper_bound_global(ll x) {
        build_compress();
        int p = rank_ub(x); return p < (int)vals.size() ? vals[p] : LLONG_MAX;
    }
    bool is_all_distinct() { build_compress(); return (int)vals.size() == n; }  // 全要素が相異なるか O(1)
    // 最も多い値の出現回数。★呼ぶたびに全ユニーク値を走査するので O(K)。
    ll max_frequency() {
        build_compress(); ll best = 0;
        for (int v = 0; v + 1 < (int)occSt.size(); v++) best = std::max<ll>(best, occSt[v + 1] - occSt[v]);
        return best;
    }
    // 最も多く出現する値そのもの。同数なら最も小さい値。★呼ぶたびに O(K)。
    ll argmax_frequency() {
        build_compress(); ll best = -1; ll bv = 0;
        for (int v = 0; v + 1 < (int)occSt.size(); v++) {
            ll c = occSt[v + 1] - occSt[v];
            if (c > best) { best = c; bv = vals[v]; }
        }
        return bv;
    }
    // (値, 出現回数) を値の昇順で全部返す。★呼ぶたびに長さ K の vector を作るので O(K)。
    vector<pair<ll, int>> frequency_table() {
        build_compress();
        int K = (int)vals.size();
        vector<pair<ll, int>> res; res.reserve(K);
        for (int v = 0; v < K; v++) res.push_back({vals[v], occSt[v + 1] - occSt[v]});
        return res;
    }
    // ソート済みのコピーを返す。★毎回ソートし直すので O(N log N)。
    vector<ll> sorted() { vector<ll> s = a; sort(s.begin(), s.end()); return s; }
    // 区間 [l,r) 内の最小値 / 最大値についての個数と出現位置。★l < r が前提。
    //   argmin/argmax (O(1)) で順位を直接得てから出現リストを二分探索するので O(log N)。
    //   first_min / first_max は O(1) (argmin / argmax そのもの)。
    ll min_count(int l, int r) { build_compress(); return count_rank(l, r, rnk[argmin(l, r)]); }  // 最小値の個数
    ll max_count(int l, int r) { build_compress(); return count_rank(l, r, rnk[argmax(l, r)]); }  // 最大値の個数
    int first_min(int l, int r) { return argmin(l, r); }   // 最小値が最初に現れる位置
    int first_max(int l, int r) { return argmax(l, r); }   // 最大値が最初に現れる位置
    int last_min(int l, int r)  { build_compress(); return last_eq_rank(l, r, rnk[argmin(l, r)]); }  // 最後の位置
    int last_max(int l, int r)  { build_compress(); return last_eq_rank(l, r, rnk[argmax(l, r)]); }  // 最後の位置

    // ====================================================== Wavelet Matrix 本体
    //   座標圧縮した順位 (0..K-1) を並べた Wavelet Matrix。レベル数 lg = ceil(log2 K)。
    //   build_wm     : 初回呼び出し時に O(N log K)。個数系 / 順位統計 / 前後の値が使う。
    //   build_wm_sum : さらに重み (元の値) の累積和を足す。O(N log K) 時間だが
    //                  long long を (lg+1)*(N+1) 個使うのでメモリが大きい。
    //                  sum 系 (sum_lt など) を初めて呼んだときだけ構築される。
    //   ★以下の区間クエリはすべて [l,r) で l < r が前提 (count_range などの一部を除く)。
    bool b_wm = false, b_wmsum = false;
    WaveletMatrix wm;
    vector<pair<ll, int>> kvBuf;                 // smallest_k / largest_k 用の作業バッファ
    void build_wm() {
        if (b_wm) return;
        b_wm = true;
        build_compress();
        int K = (int)vals.size();
        int lg = 1; while ((1 << lg) < K) lg++;
        vector<ll> keys(n);
        for (int i = 0; i < n; i++) keys[i] = rnk[i];
        wm.build(std::move(keys), lg);
    }
    void build_wm_sum() { if (b_wmsum) return;
 b_wmsum = true; build_wm(); wm.buildSum(a); }

    // ---- 個数系 : [l,r) の中で値の条件を満たす要素数。x は配列に無い値でもよい。O(log N) ----
    ll count_lt(int l, int r, ll x) { build_wm(); return wm.cntLt(l, r, rank_lb(x)); }  // a_i <  x
    ll count_le(int l, int r, ll x) { build_wm(); return wm.cntLt(l, r, rank_ub(x)); }  // a_i <= x
    ll count_ge(int l, int r, ll x) { return (r - l) - count_lt(l, r, x); }             // a_i >= x
    ll count_gt(int l, int r, ll x) { return (r - l) - count_le(l, r, x); }             // a_i >  x
    ll count_eq(int l, int r, ll x) { return count(l, r, x); }                          // a_i == x
    // x <= a_i < y の個数 (右端は含まない)。x >= y なら 0。
    //   2 回数えるのではなく、x と y の上位ビットが一致する間は下降を共有する。
    ll count_range(int l, int r, ll x, ll y) {
        if (x >= y) return 0;
        build_wm();
        return wm.cntRange(l, r, rank_lb(x), rank_lb(y));
    }
    // ---- 和系 : 条件を満たす要素の総和。★戻り値は long long なので桁あふれに注意。
    //      初回呼び出しで Wavelet Matrix の累積和 (大きめのメモリ) が構築される。O(log N) ----
    ll sum_lt(int l, int r, ll x) { build_wm_sum(); return wm.cntSumLt(l, r, rank_lb(x)).second; }  // a_i <  x
    ll sum_le(int l, int r, ll x) { build_wm_sum(); return wm.cntSumLt(l, r, rank_ub(x)).second; }  // a_i <= x
    ll sum_ge(int l, int r, ll x) { return sum(l, r) - sum_lt(l, r, x); }                           // a_i >= x
    ll sum_gt(int l, int r, ll x) { return sum(l, r) - sum_le(l, r, x); }                           // a_i >  x
    ll sum_range(int l, int r, ll x, ll y) {      // x <= a_i < y の総和 (x >= y なら 0)
        if (x >= y) return 0;
        return sum_lt(l, r, y) - sum_lt(l, r, x);
    }
    // ---- 順位統計 : ★k はすべて 0-indexed。0 <= k < r-l が前提。O(log N) ----
    ll kth(int l, int r, ll k) { build_wm(); return vals[wm.kth(l, r, k)]; }  // k 番目に小さい値
    ll kth_largest(int l, int r, ll k) { return kth(l, r, (r - l) - 1 - k); } // k 番目に大きい値
    ll median(int l, int r) { return kth(l, r, (r - l - 1) / 2); }            // 中央値(偶数長は小さい方)
    ll quantile(int l, int r, ll k) { return kth(l, r, k); }                  // kth の別名
    // 小さい方 / 大きい方から k 個の和。k <= 0 なら 0、k >= r-l なら区間の総和。O(log N)。
    ll sum_smallest(int l, int r, ll k) {
        if (k <= 0) return 0;
        if (k >= r - l) return sum(l, r);
        build_wm_sum();
        auto res = wm.smallestK(l, r, k);
        return res.sum + res.rest * vals[res.key];
    }
    ll sum_largest(int l, int r, ll k) {
        if (k <= 0) return 0;
        if (k >= r - l) return sum(l, r);
        return sum(l, r) - sum_smallest(l, r, (r - l) - k);
    }
    // 小さい方から k 個を昇順に列挙する (重複はそのまま含む)。
    //   k は min(k, r-l) に丸められ、k <= 0 なら空 vector。
    //   kth を k 回呼ぶのではなく木を 1 回たどるので
    //   O(k + (取り出した相異なる値の個数) * log N)。
    vector<ll> smallest_k(int l, int r, int k) {
        k = std::min(k, r - l);
        vector<ll> res;
        if (k <= 0) return res;
        build_wm();
        kvBuf.clear();
        ll need = k;
        wm.dfsSmall(0, l, r, 0, need, kvBuf);     // 木を 1 回たどって (値,個数) を昇順に集める
        res.reserve(k);
        for (const auto &pr : kvBuf) { ll v = vals[(size_t)pr.first]; for (int i = 0; i < pr.second; i++) res.push_back(v); }
        return res;
    }
    // 大きい方から k 個を降順に列挙する。計算量・丸め方は smallest_k と同じ。
    vector<ll> largest_k(int l, int r, int k) {
        k = std::min(k, r - l);
        vector<ll> res;
        if (k <= 0) return res;
        build_wm();
        kvBuf.clear();
        ll need = k;
        wm.dfsLarge(0, l, r, 0, need, kvBuf);
        res.reserve(k);
        for (const auto &pr : kvBuf) { ll v = vals[(size_t)pr.first]; for (int i = 0; i < pr.second; i++) res.push_back(v); }
        return res;
    }
    // ---- predecessor / successor ----
    //   ★該当する値が無いときは NONE (= LLONG_MIN) を返す。
    //   「個数を数えてから k 番目を取る」のではなく Wavelet Matrix を 1 回下降する。O(log N)。
    static constexpr ll NONE = LLONG_MIN;         // 該当なしを表す番兵
    ll prev_value(int l, int r, ll x) {           // [l,r) 内で x 以下の最大値
        build_wm();
        int ub = rank_ub(x); if (ub == 0) return NONE;
        ll k = wm.prevKey(l, r, ub - 1);
        return k < 0 ? NONE : vals[(size_t)k];
    }
    ll prev_value_strict(int l, int r, ll x) {    // [l,r) 内で x 未満の最大値
        build_wm();
        int lb = rank_lb(x); if (lb == 0) return NONE;
        ll k = wm.prevKey(l, r, lb - 1);
        return k < 0 ? NONE : vals[(size_t)k];
    }
    ll next_value(int l, int r, ll x) {           // [l,r) 内で x 以上の最小値
        build_wm();
        int lb = rank_lb(x); if (lb >= (int)vals.size()) return NONE;
        ll k = wm.nextKey(l, r, lb);
        return k < 0 ? NONE : vals[(size_t)k];
    }
    ll next_value_strict(int l, int r, ll x) {    // [l,r) 内で x より大きい最小値
        build_wm();
        int ub = rank_ub(x); if (ub >= (int)vals.size()) return NONE;
        ll k = wm.nextKey(l, r, ub);
        return k < 0 ? NONE : vals[(size_t)k];
    }
    ll max_le(int l, int r, ll x) { return prev_value(l, r, x); }   // prev_value の別名
    ll min_ge(int l, int r, ll x) { return next_value(l, r, x); }   // next_value の別名
    // [l,r) 内で x に最も近い値。★距離が同じなら小さい方 (x 以下の方) を返す。
    //   pred と succ を 1 回の下降でまとめて取得する。O(log N)。区間が空でなければ NONE は返らない。
    ll nearest(int l, int r, ll x) {
        build_wm();
        int K = (int)vals.size();
        int lb = lb_idx(vals.data(), K, x);
        ll pk, sk;
        if (lb >= K) { pk = wm.prevKey(l, r, K - 1); sk = -1; }   // x が全要素より大きい場合
        else wm.predSucc(l, r, lb, vals[lb] == x, pk, sk);
        ll p = (pk < 0) ? NONE : vals[(size_t)pk];
        ll q = (sk < 0) ? NONE : vals[(size_t)sk];
        if (p == NONE) return q;
        if (q == NONE) return p;
        return (x - p <= q - x) ? p : q;
    }
    // 2 番目に小さい / 大きい値。重複もそのまま数えるので [3,3,5] の second_min は 3。
    //   要素数が 1 以下のときは NONE。O(log N)。
    ll second_min(int l, int r) { return (r - l >= 2) ? kth(l, r, 1) : NONE; }
    ll second_max(int l, int r) { return (r - l >= 2) ? kth(l, r, r - l - 2) : NONE; }
    // 「相異なる値」として 2 番目に小さい / 大きい値。区間が 1 種類の値だけなら NONE。O(log N)。
    ll second_min_distinct(int l, int r) { return next_value_strict(l, r, min(l, r)); }
    ll second_max_distinct(int l, int r) { return prev_value_strict(l, r, max(l, r)); }
    // --- 内部用 : 圧縮後の順位 rk を直接受け取り、1 回の下降で (個数,総和) を得て
    //     Σ|a_i - x| を計算する。★build_wm_sum 済みで呼ぶこと。O(log N)。 ---
    inline ll sum_abs_rank(int l, int r, ll x, int rk) {
        auto cs = wm.cntSumLt(l, r, rk);
        ll c = cs.first, s = cs.second, tot = sum(l, r), len = r - l;
        return (x * c - s) + ((tot - s) - x * (len - c));
    }
    // [l,r) について Σ|a_i - x| を返す。★戻り値は long long。O(log N)。
    ll sum_abs(int l, int r, ll x) {
        build_wm_sum();
        return sum_abs_rank(l, r, x, rank_lb(x));
    }
    // min_x Σ|a_i - x| (最小値は中央値で達成される)。★l < r が前提。O(log N)。
    ll median_cost(int l, int r) {
        build_wm_sum();
        ll rk = wm.kth(l, r, (r - l - 1) / 2);    // 中央値の順位をそのまま使い、二分探索を省く
        return sum_abs_rank(l, r, vals[(size_t)rk], (int)rk);
    }
    // [l,r) で過半数 (出現回数 * 2 > 区間長) を占める値。無ければ NONE。★l < r が前提。O(log N)。
    //   過半数を占める値は必ず中央値になるので、中央値とその出現数を 1 回の下降で求めて判定する。
    ll majority(int l, int r) {
        build_wm();
        ll key, cnt;
        wm.kthCount(l, r, (r - l - 1) / 2, key, cnt);   // 中央値とその出現数を 1 回の下降で
        return (cnt * 2 > r - l) ? vals[(size_t)key] : NONE;
    }

    // ============================== 出現位置 (prev / next) と 区間 distinct 数
    //   build_occ     : 初回呼び出し時に O(N) で prevOcc / nextOcc / firstOcc を作る
    //                   (座標圧縮が未構築ならそれも構築される)。
    //   build_wm_prev : prevOcc を並べた Wavelet Matrix。distinct_count 専用。O(N log N)。
    bool b_occ = false, b_wmprev = false;
    vector<int> prevOcc, nextOcc, firstOcc;   // 同じ値の直前 / 直後 / 最初の出現位置
    WaveletMatrix wmPrev;                     // prevOcc[i]+1 (0..n) を並べたもの
    void build_occ() {
        if (b_occ) return;
        b_occ = true;
        build_compress();
        prevOcc.resize(n); nextOcc.resize(n); firstOcc.resize(n);
        for (int i = 0; i < n; i++) {
            int v = rnk[i], base = occSt[v], c = occSt[v + 1] - occSt[v], k = occIdx[i];
            prevOcc[i]  = (k > 0)     ? occPos[base + k - 1] : -1;
            nextOcc[i]  = (k + 1 < c) ? occPos[base + k + 1] : n;
            firstOcc[i] = occPos[base];
        }
    }
    // すべて 0 <= i < n が前提、クエリは O(1)。
    int prev_occurrence(int i) { build_occ(); return prevOcc[i]; }   // a[i] と同じ値の直前の位置 (無ければ -1)
    int next_occurrence(int i) { build_occ(); return nextOcc[i]; }   // a[i] と同じ値の直後の位置 (無ければ n)
    int first_occurrence(int i){ build_occ(); return firstOcc[i]; }  // a[i] と同じ値が最初に現れる位置
    int prev_distance(int i) { build_occ(); return prevOcc[i] < 0 ? -1 : i - prevOcc[i]; }  // 直前までの距離(無ければ -1)
    int next_distance(int i) { build_occ(); return nextOcc[i] >= n ? -1 : nextOcc[i] - i; } // 直後までの距離(無ければ -1)
    bool has_duplicate() { return !is_all_distinct(); }              // 重複があるか         O(1)
    void build_wm_prev() {
        if (b_wmprev) return;
        b_wmprev = true;
        build_occ();
        int lg = 1; while ((1 << lg) <= n) lg++;
        vector<ll> keys(n);
        for (int i = 0; i < n; i++) keys[i] = prevOcc[i] + 1;   // 0..n
        wmPrev.build(std::move(keys), lg);
    }
    // [l,r) に含まれる相異なる値の個数。O(log N)。
    //   ★l >= r (空区間) のときは 0 を返すので安全。
    //   「その値が区間内で最初に現れる位置か = prevOcc[i] < l か」を数えている。
    ll distinct_count(int l, int r) {
        if (l >= r) return 0;
        build_wm_prev();
        return wmPrev.cntLt(l, r, (ll)l + 1);
    }
    // a[0..i-1] に含まれる相異なる値の個数 (0 <= i <= n)。O(log N)。
    ll distinct_prefix(int i) { return distinct_count(0, i); }

    // ============================================================== 区間 mex
    //   初回呼び出し時に永続セグメント木を O(N log N) で構築する。
    //   ★ノードを N * O(log N) 個作るのでメモリは大きめ。
    //   ★値域 [0, n] の外にある要素 (負数や n より大きい値) は mex に影響しないので
    //     木には入れず、version だけ進める。
    bool b_mex = false;
    PersistentMinSeg pmex;
    void build_mex() {
        if (b_mex) return;
        b_mex = true;
        int m = n + 1;                                  // mex は高々 n
        int lgm = 1; while ((1 << lgm) < m) lgm++;
        pmex.init(m, 2 + n * (lgm + 2));
        for (int i = 0; i < n; i++) {
            if (a[i] >= 0 && a[i] < m) pmex.push((int)a[i], i);
            else pmex.root.push_back(pmex.root.back());
        }
    }
    // [l,r) に含まれない最小の非負整数 (区間 mex)。O(log N)。
    //   ★0 <= l <= r <= n。戻り値は 0..(r-l) の範囲。
    ll mex(int l, int r) { build_mex(); return pmex.findFirstLess(r, l); }
    ll mex() { return mex(0, n); }              // 配列全体の mex。O(log N)。

    // =========================================================== 区間 最頻値
    //   ブロック分割 + 出現位置リストによる古典的な静的 range mode。
    //   初回呼び出し時の構築 O(N sqrt N) 時間 / O(N) メモリ、クエリ O(sqrt N)。
    //   ★このテンプレートの中で唯一 O(sqrt N) のクエリ。多用すると重い。
    //   ★最頻値が複数あるときにどれを返すかは実装依存 (ブロックの取り方に依存)。
    //     出現回数のほうは常に正しい最大値になる。
    bool b_mode = false;
    int modeB = 1, modeNB = 0;       // modeB = ブロック幅 (≒ sqrt(N)), modeNB = ブロック数
    vector<int> modeVal, modeCnt;    // ブロック bi..bj をちょうど覆う区間の (最頻値の順位, 出現数)
    vector<int> mPos, mLo, mHi;      // i ごとの occPos 上の位置 / 所属グループの範囲 (連続アクセス用)
    struct ModeCnt { int stamp, c; };  // stamp = このカウンタを最後に使ったブロック番号
    void build_mode() {
        if (b_mode) return;
        b_mode = true;
        build_compress();
        if (n == 0) return;
        modeB = std::max(1, (int)(sqrt((double)n) + 1));
        modeNB = (n + modeB - 1) / modeB;
        modeVal.assign((size_t)modeNB * modeNB, 0);
        modeCnt.assign((size_t)modeNB * modeNB, 0);
        // 走査のたびにカウンタを 0 に戻す代わりに「どのブロック起点か」を同居させ、
        // リセット用の全走査 (毎回 O(n) のランダム書き込み) を丸ごと無くす
        vector<ModeCnt> cnt(vals.size(), ModeCnt{-1, 0});
        const int *rk = rnk.data();
        for (int bi = 0; bi < modeNB; bi++) {
            int st = bi * modeB, best = 0, bestv = 0, bj = bi;
            int mark = std::min(n, (bj + 1) * modeB) - 1;
            for (int j = st; j < n; j++) {
                int v = rk[j];
                ModeCnt &e = cnt[v];
                if (e.stamp != bi) { e.stamp = bi; e.c = 0; }
                if (++e.c > best) { best = e.c; bestv = v; }
                if (j == mark) {
                    modeVal[(size_t)bi * modeNB + bj] = bestv;
                    modeCnt[(size_t)bi * modeNB + bj] = best;
                    bj++;
                    mark = std::min(n, (bj + 1) * modeB) - 1;
                }
            }
        }
        mPos.resize(n); mLo.resize(n); mHi.resize(n);
        for (int i = 0; i < n; i++) {
            int v = rnk[i], lo = occSt[v];
            mLo[i] = lo; mHi[i] = occSt[v + 1]; mPos[i] = lo + occIdx[i];
        }
    }
    // [l,r) の (最頻値, その出現回数) を返す。O(sqrt N)。
    //   ★l >= r のときは (0, 0) を返す (空区間でも落ちない)。
    //   ★最頻値が複数あるときはそのうちの 1 つ。出現回数は常に正しい。
    //   内側のブロック範囲は前計算表から取り、はみ出した両端だけ
    //   出現位置リストを使って「あと 1 個増やせるか」を繰り返し確認する。
    pair<ll, ll> mode(int l, int r) {
        if (l >= r) return {0, 0};
        build_mode();
        int bl = (l + modeB - 1) / modeB, br = r / modeB;
        const int *rk = rnk.data(), *op = occPos.data();
        const int *mp = mPos.data(), *mlo = mLo.data(), *mhi = mHi.data();
        int best = 0, bestv = rk[l];
        // 大半の p では best が更新されないので、更新が起きたときだけ rk[p] を読む
        if (bl >= br) {
            for (int p = l; p < r; p++) {
                int k = mp[p], hi = mhi[p];
                if (k + best < hi && op[k + best] < r) {
                    do { best++; } while (k + best < hi && op[k + best] < r);
                    bestv = rk[p];
                }
            }
        } else {
            bestv = modeVal[(size_t)bl * modeNB + (br - 1)];
            best  = modeCnt[(size_t)bl * modeNB + (br - 1)];
            int le = bl * modeB, rs = br * modeB;
            for (int p = l; p < le; p++) {
                int k = mp[p], hi = mhi[p];
                if (k + best < hi && op[k + best] < r) {
                    do { best++; } while (k + best < hi && op[k + best] < r);
                    bestv = rk[p];
                }
            }
            for (int p = rs; p < r; p++) {
                int k = mp[p], lo = mlo[p];
                if (k - best >= lo && op[k - best] >= l) {
                    do { best++; } while (k - best >= lo && op[k - best] >= l);
                    bestv = rk[p];
                }
            }
        }
        return {vals[bestv], best};
    }

    // ======================================== 直近のより小さい / 大きい要素
    //   初回呼び出し時に単調スタックで 8 種類 + Cartesian Tree の親を O(N) で作る。
    //   以降のクエリはすべて O(1) (largest_rectangle だけは毎回 O(N) の走査)。
    //   ★見つからないとき、左向きは -1、右向きは n を返す。
    //   ★0 <= i < n が前提。
    bool b_near = false;
    vector<int> pSm, nSm, pSmE, nSmE, pGr, nGr, pGrE, nGrE, cartPar;
    void build_near() {
        if (b_near) return;
        b_near = true;
        pSm.assign(n, -1); nSm.assign(n, n); pSmE.assign(n, -1); nSmE.assign(n, n);
        pGr.assign(n, -1); nGr.assign(n, n); pGrE.assign(n, -1); nGrE.assign(n, n);
        cartPar.assign(n, -1);
        vector<int> st;
        st.clear(); for (int i = 0; i < n; i++) { while (!st.empty() && a[st.back()] >= a[i]) st.pop_back(); pSm[i]  = st.empty() ? -1 : st.back(); st.push_back(i); }
        st.clear(); for (int i = n - 1; i >= 0; i--) { while (!st.empty() && a[st.back()] >= a[i]) st.pop_back(); nSm[i]  = st.empty() ?  n : st.back(); st.push_back(i); }
        st.clear(); for (int i = 0; i < n; i++) { while (!st.empty() && a[st.back()] >  a[i]) st.pop_back(); pSmE[i] = st.empty() ? -1 : st.back(); st.push_back(i); }
        st.clear(); for (int i = n - 1; i >= 0; i--) { while (!st.empty() && a[st.back()] >  a[i]) st.pop_back(); nSmE[i] = st.empty() ?  n : st.back(); st.push_back(i); }
        st.clear(); for (int i = 0; i < n; i++) { while (!st.empty() && a[st.back()] <= a[i]) st.pop_back(); pGr[i]  = st.empty() ? -1 : st.back(); st.push_back(i); }
        st.clear(); for (int i = n - 1; i >= 0; i--) { while (!st.empty() && a[st.back()] <= a[i]) st.pop_back(); nGr[i]  = st.empty() ?  n : st.back(); st.push_back(i); }
        st.clear(); for (int i = 0; i < n; i++) { while (!st.empty() && a[st.back()] <  a[i]) st.pop_back(); pGrE[i] = st.empty() ? -1 : st.back(); st.push_back(i); }
        st.clear(); for (int i = n - 1; i >= 0; i--) { while (!st.empty() && a[st.back()] <  a[i]) st.pop_back(); nGrE[i] = st.empty() ?  n : st.back(); st.push_back(i); }
        // 最小値 Cartesian Tree の親 : 左右の「自分より小さい直近の要素」のうち
        // 値が大きい方を親とする (同値なら左側を親にする)。根は -1。
        for (int i = 0; i < n; i++) {
            int L = pSm[i], R = nSm[i];
            if (L < 0 && R >= n) cartPar[i] = -1;
            else if (L < 0) cartPar[i] = R;
            else if (R >= n) cartPar[i] = L;
            else cartPar[i] = (a[L] >= a[R]) ? L : R;
        }
    }
    int prev_smaller(int i)       { build_near(); return pSm[i];  }  // a[j] <  a[i] となる最大の j<i
    int next_smaller(int i)       { build_near(); return nSm[i];  }  // a[j] <  a[i] となる最小の j>i
    int prev_smaller_equal(int i) { build_near(); return pSmE[i]; }  // a[j] <= a[i] となる最大の j<i
    int next_smaller_equal(int i) { build_near(); return nSmE[i]; }  // a[j] <= a[i] となる最小の j>i
    int prev_greater(int i)       { build_near(); return pGr[i];  }  // a[j] >  a[i] となる最大の j<i
    int next_greater(int i)       { build_near(); return nGr[i];  }  // a[j] >  a[i] となる最小の j>i
    int prev_greater_equal(int i) { build_near(); return pGrE[i]; }  // a[j] >= a[i] となる最大の j<i
    int next_greater_equal(int i) { build_near(); return nGrE[i]; }  // a[j] >= a[i] となる最小の j>i
    int cartesian_parent(int i)   { build_near(); return cartPar[i]; }  // 最小値 Cartesian Tree の親(根は -1)
    // a[i] が最小値でいられる最大の区間 [l, r)。O(1)。
    pair<int, int> minimum_containing_interval(int i) { build_near(); return {pSm[i] + 1, nSm[i]}; }
    // ヒストグラム最大長方形 max(a[i] * 幅)。★非負の値を想定 (負があると意味を持たない)。
    //   ★呼ぶたびに全要素を走査するので O(N)。戻り値は __int128。
    i128 largest_rectangle() {
        build_near(); i128 best = 0;
        for (int i = 0; i < n; i++) best = std::max(best, (i128)a[i] * (nSm[i] - pSm[i] - 1));
        return best;
    }

    // ================================= next_greater / next_smaller の binary lifting
    //   初回呼び出し時に O(N log N) で 2 本のダブリング表を作る (int を 2*(N+1)*log N 個)。
    //   以降 jump_next_greater / jump_next_smaller が O(log N)。
    bool b_lift = false;
    int liftK = 1;
    vector<int> liftG, liftS;
    void build_lift() {
        if (b_lift) return;
        b_lift = true;
        build_near();
        liftK = 1; while ((1 << liftK) <= n) liftK++;
        liftG.assign((size_t)liftK * (n + 1), n);
        liftS.assign((size_t)liftK * (n + 1), n);
        for (int i = 0; i < n; i++) { liftG[i] = nGr[i]; liftS[i] = nSm[i]; }
        liftG[n] = n; liftS[n] = n;
        for (int k = 1; k < liftK; k++)
            for (int i = 0; i <= n; i++) {
                liftG[(size_t)k * (n + 1) + i] = liftG[(size_t)(k - 1) * (n + 1) + liftG[(size_t)(k - 1) * (n + 1) + i]];
                liftS[(size_t)k * (n + 1) + i] = liftS[(size_t)(k - 1) * (n + 1) + liftS[(size_t)(k - 1) * (n + 1) + i]];
            }
    }
    // next_greater を k 回たどった位置。途中で n に到達したら n を返す。O(log N)。
    //   ★0 <= i <= n、k >= 0。k が大きすぎる場合も n を返す。
    int jump_next_greater(int i, ll k) {
        build_lift();
        for (int b = 0; b < liftK && i < n; b++) if ((k >> b) & 1) i = liftG[(size_t)b * (n + 1) + i];
        return (k >> liftK) ? n : i;
    }
    // next_smaller を k 回たどった位置。仕様は jump_next_greater と同じ。O(log N)。
    int jump_next_smaller(int i, ll k) {
        build_lift();
        for (int b = 0; b < liftK && i < n; b++) if ((k >> b) & 1) i = liftS[(size_t)b * (n + 1) + i];
        return (k >> liftK) ? n : i;
    }

    // ============================ 単調性 / 等値 run / 等差
    //   初回呼び出し時に O(N) で「その位置から右にどこまで条件が続くか」を作る。
    //   以降クエリはすべて O(1)。
    //   ★区間は [l,r)。長さ 0 / 1 の区間は常に true (is_arithmetic は長さ 2 以下も true)。
    bool b_run = false;
    vector<int> ndR, niR, siR, sdR, eqL, eqR, dEqR;
    void build_run() {
        if (b_run) return;
        b_run = true;
        ndR.assign(n, 0); niR.assign(n, 0); siR.assign(n, 0); sdR.assign(n, 0);
        eqL.assign(n, 0); eqR.assign(n, 0);
        for (int i = n - 1; i >= 0; i--) {
            if (i == n - 1) { ndR[i] = niR[i] = siR[i] = sdR[i] = eqR[i] = i; }
            else {
                ndR[i] = (a[i] <= a[i + 1]) ? ndR[i + 1] : i;
                niR[i] = (a[i] >= a[i + 1]) ? niR[i + 1] : i;
                siR[i] = (a[i] <  a[i + 1]) ? siR[i + 1] : i;
                sdR[i] = (a[i] >  a[i + 1]) ? sdR[i + 1] : i;
                eqR[i] = (a[i] == a[i + 1]) ? eqR[i + 1] : i;
            }
        }
        for (int i = 0; i < n; i++) eqL[i] = (i == 0 || a[i] != a[i - 1]) ? i : eqL[i - 1];
        int m = std::max(0, n - 1);
        dEqR.assign(std::max(1, m), 0);
        for (int i = m - 1; i >= 0; i--) {
            if (i == m - 1) dEqR[i] = i;
            else dEqR[i] = (a[i + 1] - a[i] == a[i + 2] - a[i + 1]) ? dEqR[i + 1] : i;
        }
    }
    // [l,r) が各種の単調性を満たすか。長さ 1 以下は常に true。すべて O(1)。
    bool is_non_decreasing(int l, int r) { if (r - l <= 1) return true; build_run(); return ndR[l] >= r - 1; }
    bool is_non_increasing(int l, int r) { if (r - l <= 1) return true; build_run(); return niR[l] >= r - 1; }
    bool is_strictly_increasing(int l, int r) { if (r - l <= 1) return true; build_run(); return siR[l] >= r - 1; }
    bool is_strictly_decreasing(int l, int r) { if (r - l <= 1) return true; build_run(); return sdR[l] >= r - 1; }
    bool is_constant(int l, int r) { if (r - l <= 1) return true; build_run(); return eqR[l] >= r - 1; }
    // [l,r) が等差数列か。★長さ 2 以下は常に true。O(1)。
    bool is_arithmetic(int l, int r) { if (r - l <= 2) return true; build_run(); return dEqR[l] >= r - 2; }
    bool is_sorted() { return is_non_decreasing(0, n); }   // 配列全体が非減少か  O(1)
    // a[i..返り値] が非減少となる最大の index (i 自身から右にどこまで非減少か)。O(1)。
    int next_break_non_decreasing(int i) { build_run(); return ndR[i]; }
    // 位置 i を含む「同じ値が連続する区間」の左端 / 右端 / 長さ。0 <= i < n。すべて O(1)。
    int run_left(int i)   { build_run(); return eqL[i]; }
    int run_right(int i)  { build_run(); return eqR[i]; }
    int run_length(int i) { build_run(); return eqR[i] - eqL[i] + 1; }
    // ランレングス圧縮 (値, 連続長) を左から順に返す。★呼ぶたびに全走査するので O(N)。
    vector<pair<ll, int>> runs() {
        vector<pair<ll, int>> res;
        const ll *p = a.data();
        for (int i = 0; i < n;) { int j = i + 1; while (j < n && p[j] == p[i]) j++; res.push_back({p[i], j - i}); i = j; }
        return res;
    }
    // 位置 i が極小 / 極大か。★両隣より「真に」小さい / 大きいときだけ true。
    //   端 (i == 0 / n-1) は外側に隣が無いものとして判定する。O(1)。
    bool is_local_min(int i) const { return (i == 0 || a[i] < a[i - 1]) && (i == n - 1 || a[i] < a[i + 1]); }
    bool is_local_max(int i) const { return (i == 0 || a[i] > a[i - 1]) && (i == n - 1 || a[i] > a[i + 1]); }
    // 極小 / 極大の位置を昇順に全部返す。★呼ぶたびに全走査するので O(N)。
    vector<int> valleys() {
        int c = 0; for (int i = 0; i < n; i++) c += is_local_min(i);
        vector<int> r; r.reserve(c);
        for (int i = 0; i < n; i++) if (is_local_min(i)) r.push_back(i);
        return r;
    }
    vector<int> peaks() {
        int c = 0; for (int i = 0; i < n; i++) c += is_local_max(i);
        vector<int> r; r.reserve(c);
        for (int i = 0; i < n; i++) if (is_local_max(i)) r.push_back(i);
        return r;
    }
    // 符号の変化回数 : 隣り合う 2 要素が (正, 負) または (負, 正) になっている箇所の数。
    //   ★0 は符号の変化として数えない。初回呼び出し時に O(N) で構築、クエリは O(1)。
    bool b_sc = false;
    vector<int> scPre;
    void build_signchange() {
        if (b_sc) return;
        b_sc = true;
        scPre.assign(std::max(1, n), 0);
        for (int i = 0; i + 1 < n; i++)
            scPre[i + 1] = scPre[i] + (((a[i] > 0 && a[i + 1] < 0) || (a[i] < 0 && a[i + 1] > 0)) ? 1 : 0);
    }
    // [l,r) の内部 (a[l..r-1] の隣接ペア) での符号変化回数。長さ 1 以下なら 0。O(1)。
    ll count_sign_changes(int l, int r) {
        if (r - l <= 1) return 0;
        build_signchange();
        return scPre[r - 1] - scPre[l];
    }

    // ================================= しきい値を満たす最初 / 最後の位置
    //   ★next_* は「見つからなければ n」、prev_* は「見つからなければ -1」を返す。
    //   ★i は範囲外でもよい (next_* に i >= n を渡すと n、prev_* に i < 0 を渡すと -1)。
    //   ★区間ではなく「i 以降 / i 以前の配列全体」を探す点に注意 (右端 / 左端の指定は無い)。
    //   内部では RMQ のブロック Sparse Table を 1 回下降するので O(log N)。
    //   初回呼び出しで最大値 RMQ (ge/gt) または最小値 RMQ (le/lt) が構築される。
    int next_ge(int i, ll x) {                  // j>=i で a[j] >= x となる最小の j
        build_rmax(); return rmax.findFirst(a.data(), i, x, false);
    }
    int next_gt(int i, ll x) { build_rmax(); return rmax.findFirst(a.data(), i, x, true);  }  // a[j] >  x
    int next_le(int i, ll x) { build_rmin(); return rmin.findFirst(a.data(), i, x, false); }  // a[j] <= x
    int next_lt(int i, ll x) { build_rmin(); return rmin.findFirst(a.data(), i, x, true);  }  // a[j] <  x
    int prev_ge(int i, ll x) {                  // j<=i で a[j] >= x となる最大の j
        build_rmax(); return rmax.findLast(a.data(), i, x, false);
    }
    int prev_gt(int i, ll x) { build_rmax(); return rmax.findLast(a.data(), i, x, true);  }   // a[j] >  x
    int prev_le(int i, ll x) { build_rmin(); return rmin.findLast(a.data(), i, x, false); }   // a[j] <= x
    int prev_lt(int i, ll x) { build_rmin(); return rmin.findLast(a.data(), i, x, true);  }   // a[j] <  x

    // ======================== 累積和に対する検索 / 区間の最大部分和
    //   累積和配列 ps[0..n] に対して RMQ と専用セグメント木を作る。
    //   build_psrmq : O(N) 、first/last_sum_* と max/min_prefix/suffix_sum が使う
    //   build_sst   : O(N) 、max/min_subarray_sum が使う
    //   ★元の配列に負数があっても正しく動く (累積和の単調性は仮定していない)。
    bool b_psr = false, b_sst = false;
    RMQ64<ll, false> psMin;                  // 累積和の最小位置を返す RMQ
    RMQ64<ll, true>  psMax;                  // 累積和の最大位置を返す RMQ
    SubSegTree   sst;                        // 最大 / 最小部分和用
    void build_psrmq() { if (b_psr) return;
 b_psr = true; psMin.build(ps.data(), n + 1); psMax.build(ps.data(), n + 1); }
    void build_sst()   { if (b_sst) return;
 b_sst = true; sst.build(ps); }
    // ★以下 4 つは 0 <= l <= n / 0 <= r <= n を取り、返る位置も 0..n の範囲。
    //   条件を満たす位置が無ければ -1。いずれも O(log N)。
    //   r == l のとき sum(l,r) == 0 なので、x <= 0 なら first_sum_ge(l,x) == l。
    // sum(l,r) >= x となる最小の r
    int first_sum_ge(int l, ll x) {
        build_psrmq();
        int j = psMax.findFirst(ps.data(), l, ps[l] + x, false);
        return j > n ? -1 : j;
    }
    // sum(l,r) <= x となる最小の r
    int first_sum_le(int l, ll x) {
        build_psrmq();
        int j = psMin.findFirst(ps.data(), l, ps[l] + x, false);
        return j > n ? -1 : j;
    }
    // sum(l,r) >= x となる最大の l
    int last_sum_ge(int r, ll x) {
        build_psrmq();
        return psMin.findLast(ps.data(), r, ps[r] - x, false);
    }
    // sum(l,r) <= x となる最大の l
    int last_sum_le(int r, ll x) {
        build_psrmq();
        return psMax.findLast(ps.data(), r, ps[r] - x, false);
    }
    // [l,r) に完全に含まれる連続部分列の和の最大 / 最小。★空の部分列は許さない
    //   (必ず 1 要素以上を選ぶ)。l < r が前提。O(log N)。
    ll max_subarray_sum(int l, int r) { build_sst(); return sst.query(l, r + 1).best;  }
    ll min_subarray_sum(int l, int r) { build_sst(); return sst.query(l, r + 1).worst; }
    // 以下 4 つは ★l < r が前提、すべて O(1)。
    //   max/min_prefix_sum : 左端を l に固定した和 sum(l,j) (l < j <= r) の最大 / 最小
    //   max/min_suffix_sum : 右端を r に固定した和 sum(i,r) (l <= i < r) の最大 / 最小
    ll max_prefix_sum(int l, int r) { build_psrmq(); return ps[psMax.query(ps.data(), l + 1, r + 1)] - ps[l]; }
    ll max_suffix_sum(int l, int r) { build_psrmq(); return ps[r] - ps[psMin.query(ps.data(), l, r)]; }
    ll min_prefix_sum(int l, int r) { build_psrmq(); return ps[psMin.query(ps.data(), l + 1, r + 1)] - ps[l]; }
    ll min_suffix_sum(int l, int r) { build_psrmq(); return ps[r] - ps[psMax.query(ps.data(), l, r)]; }

    // ============================================ XOR 系
    //   ★前提 : 全要素が 0 以上であること (build_wm_bits の assert で保護)。
    //   座標圧縮ではなく「生の値のビット」を並べた Wavelet Matrix を作る。
    //   A = 最大値のビット幅 (bitLg)。初回構築 O(N log A)、クエリ O(log A)。
    //   ★x に bitLg より上位のビットがあっても正しく扱える
    //     (上位ビットは全要素で 0 なので、定数シフトとして先に処理している)。
    bool b_wmb = false;
    int bitLg = 1;                            // 扱うビット数 (最大値のビット幅)
    WaveletMatrix wmb;
    void build_wm_bits() {
        if (b_wmb) return;
        b_wmb = true;
        assert(allNonNeg && "XOR 系クエリは全要素が 0 以上のときのみ使えます");
        ll mx = 0; for (int i = 0; i < n; i++) mx = std::max(mx, a[i]);
        bitLg = 1; while (bitLg < 63 && (mx >> bitLg)) bitLg++;
        wmb.build(a, bitLg);
    }
    // #{ i in [l,r) : (a_i ^ x) < bound }。★x >= 0 が前提。O(log A)。
    //   bound <= 0 なら 0。
    ll count_xor_lt(int l, int r, ll x, ll bound) {
        if (bound <= 0) return 0;
        build_wm_bits();
        ull mask = (bitLg >= 64) ? ~0ULL : ((1ULL << bitLg) - 1);
        ull xu = (ull)x, base = xu & ~mask, xl = xu & mask, b = (ull)bound;
        if (b <= base) return 0;
        if (b >= base + mask + 1) return r - l;
        return wmb.cntXorLt(l, r, xl, b - base);
    }
    // ★count_xor_le / gt は bound+1 を使うので、bound に LLONG_MAX を渡さないこと。
    ll count_xor_le(int l, int r, ll x, ll bound) { return count_xor_lt(l, r, x, bound + 1); }
    ll count_xor_ge(int l, int r, ll x, ll bound) { return (r - l) - count_xor_lt(l, r, x, bound); }
    ll count_xor_gt(int l, int r, ll x, ll bound) { return (r - l) - count_xor_le(l, r, x, bound); }
    // [l,r) の要素について (a_i ^ x) を昇順に並べたときの k 番目 (0-indexed) の値。
    //   ★0 <= k < r-l、x >= 0 が前提。O(log A)。
    ll kth_xor(int l, int r, ll x, ll k) {
        build_wm_bits();
        ull mask = (bitLg >= 64) ? ~0ULL : ((1ULL << bitLg) - 1);
        ull xu = (ull)x, base = xu & ~mask, xl = xu & mask;
        return (ll)(base + wmb.kthXor(l, r, xl, k));
    }
    ll max_xor(int l, int r, ll x) { return kth_xor(l, r, x, (r - l) - 1); }  // max (a_i ^ x)
    ll min_xor(int l, int r, ll x) { return kth_xor(l, r, x, 0); }            // min (a_i ^ x)

    // ================================ 区間の線形基底 (XOR 基底)
    //   ★前提 : 全要素が 0 以上であること (build_basis の assert で保護)。
    //   「prefix [0,r) の基底を、各ベクトルが由来する位置つきで持つ」方式。
    //   クエリでは位置が l 以上の基底ベクトルだけを使えば [l,r) の基底になる。
    //   初回構築 O(N log A)、メモリは (N+1) * A 個ぶんの (long long + int) で
    //   このテンプレートの中では特に大きい。クエリは O(log A) (kth だけ O(log^2 A))。
    bool b_basis = false;
    int basLg = 1;                            // 扱うビット数 (最大値のビット幅) = A
    vector<ll> basV;                          // basV[r*basLg + b] : prefix [0,r) の bit b の基底
    vector<int> basP;                         // basP[r*basLg + b] : そのベクトルが由来する位置 (無ければ -1)
    void build_basis() {
        if (b_basis) return;
        b_basis = true;
        assert(allNonNeg && "線形基底は全要素が 0 以上のときのみ使えます");
        ll mx = 0; for (int i = 0; i < n; i++) mx = std::max(mx, a[i]);
        basLg = 1; while (basLg < 63 && (mx >> basLg)) basLg++;
        // 行 1..n は直前の行のコピーで必ず上書きされるので、行 0 だけ初期化すればよい
        basV.assign((size_t)(n + 1) * basLg, 0);
        basP.resize((size_t)(n + 1) * basLg);
        std::fill(basP.begin(), basP.begin() + basLg, -1);
        for (int r = 1; r <= n; r++) {
            ll *pv = &basV[(size_t)r * basLg]; int *pp = &basP[(size_t)r * basLg];
            const ll *qv = &basV[(size_t)(r - 1) * basLg]; const int *qp = &basP[(size_t)(r - 1) * basLg];
            for (int b = 0; b < basLg; b++) { pv[b] = qv[b]; pp[b] = qp[b]; }
            ll v = a[r - 1]; int p = r - 1;
            for (int b = basLg - 1; b >= 0 && v; b--) {
                if (!((v >> b) & 1)) continue;
                if (pv[b] == 0) { pv[b] = v; pp[b] = p; v = 0; break; }
                if (pp[b] < p) { std::swap(pv[b], v); std::swap(pp[b], p); }
                v ^= pv[b];
            }
        }
    }
    // [l,r) の要素が張る XOR 空間の次元 (作れる値は 2^次元 通り)。O(log A)。
    int basis_rank(int l, int r) {
        build_basis(); int c = 0;
        const int *pp = &basP[(size_t)r * basLg];
        for (int b = 0; b < basLg; b++) if (pp[b] >= l) c++;
        return c;
    }
    // [l,r) の要素から選ぶ部分集合の XOR の最大値。空集合も選べるので常に 0 以上。O(log A)。
    ll basis_max_xor(int l, int r) {
        build_basis(); ll res = 0;
        const ll *pv = &basV[(size_t)r * basLg]; const int *pp = &basP[(size_t)r * basLg];
        for (int b = basLg - 1; b >= 0; b--) if (pp[b] >= l && !((res >> b) & 1)) res ^= pv[b];
        return res;
    }
    // x を [l,r) の要素の部分集合 XOR で作れるか。x < 0 なら false。O(log A)。
    bool basis_contains(int l, int r, ll x) {
        if (x < 0) return false;
        build_basis();
        const ll *pv = &basV[(size_t)r * basLg]; const int *pp = &basP[(size_t)r * basLg];
        for (int b = basLg - 1; b >= 0; b--) if ((x >> b) & 1) {
            if (pp[b] < l) return false;
            x ^= pv[b];
        }
        return x == 0;
    }
    // 部分集合 XOR で作れる値を昇順に並べたときの k 番目 (0-indexed, 空集合の 0 を含む)。
    //   ★作れる値が 2^rank 通りしかないので、k がそれ以上なら -1 を返す。
    //   基底を簡約するため、この関数だけ O(log^2 A) かかる (一時 vector も 2 本作る)。
    ll basis_kth_xor(int l, int r, ll k) {
        build_basis();
        const ll *pv = &basV[(size_t)r * basLg]; const int *pp = &basP[(size_t)r * basLg];
        vector<ll> t(basLg, 0);
        for (int b = 0; b < basLg; b++) if (pp[b] >= l) t[b] = pv[b];
        for (int b = 0; b < basLg; b++) if (t[b])
            for (int b2 = b + 1; b2 < basLg; b2++) if (t[b2] && ((t[b2] >> b) & 1)) t[b2] ^= t[b];
        vector<ll> bs;
        for (int b = 0; b < basLg; b++) if (t[b]) bs.push_back(t[b]);
        int R = (int)bs.size();
        if (R < 62 && k >= (1LL << R)) return -1;
        ll res = 0;
        for (int i = 0; i < R; i++) if ((k >> i) & 1) res ^= bs[i];
        return res;
    }

    // ====================== Suffix Array + LCP
    //   配列を「文字列」と見なして接尾辞配列と LCP 配列を作る (座標圧縮した順位を使う)。
    //   初回呼び出し時に O(N log N) で構築、以降 lcp_suffix / equal / compare などが O(1)。
    //   ★ハッシュと違って決定的なので、部分列の一致判定はこちらが確実。
    bool b_sa = false;
    //   sa[k]     = 辞書順 k 番目の接尾辞の開始位置
    //   saRank[i] = 接尾辞 a[i..] の辞書順の順位
    //   lcpArr[k] = LCP(a[sa[k]..], a[sa[k+1]..])   (長さ n-1)
    vector<int> sa, saRank, lcpArr;
    RMQ64<int, false> lcpRMQ;        // LCP 配列は int のまま扱う (ll への複製をやめてメモリ半減)
    // ダブリング + 基数ソートによる接尾辞配列。s は 0..upper-1 の値の列。O(N log N)。
    static vector<int> build_sa_doubling(const vector<int> &s, int upper) {
        int n = (int)s.size();
        vector<int> sa(n), rk(s), tmp(n), sa2(n);
        vector<int> cnt(std::max(upper, n) + 2, 0);
        for (int i = 0; i < n; i++) cnt[rk[i]]++;
        for (int i = 1; i < upper; i++) cnt[i] += cnt[i - 1];
        for (int i = n - 1; i >= 0; i--) sa[--cnt[rk[i]]] = i;
        int m = upper;
        for (int k = 1; k < n; k <<= 1) {
            int p = 0;
            for (int i = n - k; i < n; i++) sa2[p++] = i;
            for (int i = 0; i < n; i++) if (sa[i] >= k) sa2[p++] = sa[i] - k;
            fill(cnt.begin(), cnt.begin() + m + 1, 0);
            for (int i = 0; i < n; i++) cnt[rk[i]]++;
            for (int i = 1; i < m; i++) cnt[i] += cnt[i - 1];
            for (int i = n - 1; i >= 0; i--) sa[--cnt[rk[sa2[i]]]] = sa2[i];
            tmp[sa[0]] = 0;
            for (int i = 1; i < n; i++) {
                int x = sa[i - 1], y = sa[i];
                bool same = (rk[x] == rk[y]) &&
                            ((x + k < n ? rk[x + k] : -1) == (y + k < n ? rk[y + k] : -1));
                tmp[y] = tmp[x] + (same ? 0 : 1);
            }
            rk = tmp;
            m = rk[sa[n - 1]] + 1;
            if (m == n) break;
        }
        return sa;
    }
    void build_suffix_array() {
        if (b_sa) return;
        b_sa = true;
        if (n == 0) return;
        build_compress();
        sa = build_sa_doubling(rnk, (int)vals.size());
        saRank.assign(n, 0);
        for (int i = 0; i < n; i++) saRank[sa[i]] = i;
        lcpArr.assign(std::max(1, n - 1), 0);
        int h = 0;
        for (int i = 0; i < n; i++) {
            if (saRank[i] > 0) {
                int j = sa[saRank[i] - 1];
                while (i + h < n && j + h < n && a[i + h] == a[j + h]) h++;
                lcpArr[saRank[i] - 1] = h;
                if (h) h--;
            } else h = 0;
        }
        if (n >= 2) lcpRMQ.build(lcpArr.data(), n - 1);
    }
    // 接尾辞配列そのものへの参照 (辞書順に並べた開始位置)。★n >= 1 が前提。O(1)。
    const vector<int>& suffix_array() { build_suffix_array(); return sa; }
    // 接尾辞 a[i..] の辞書順の順位 (0..n-1)。0 <= i < n。O(1)。
    int suffix_rank(int i) { build_suffix_array(); return saRank[i]; }
    // 接尾辞 a[i..] と a[j..] の最長共通接頭辞長。0 <= i,j < n。O(1)。
    //   i == j のときは n-i (接尾辞の長さそのもの) を返す。
    int lcp_suffix(int i, int j) {
        if (i == j) return n - i;
        build_suffix_array();
        int x = saRank[i], y = saRank[j];
        if (x > y) std::swap(x, y);
        return lcpArr[lcpRMQ.query(lcpArr.data(), x, y)];
    }
    // lcp_suffix に上限 maxlen を付けたもの。O(1)。
    int lcp(int i, int j, int maxlen) { return std::min(lcp_suffix(i, j), maxlen); }
    // 接尾辞 a[i..] と a[j..] の辞書順比較。i<j なら -1、同じなら 0、i>j なら 1。O(1)。
    int compare_suffix(int i, int j) {
        if (i == j) return 0;
        build_suffix_array();
        return saRank[i] < saRank[j] ? -1 : 1;
    }
    // a[l1,r1) と a[l2,r2) が列として等しいか。長さが違えば false。O(1)。
    //   ★ハッシュではなく LCP による判定なので衝突は起こらない (決定的)。
    bool equal(int l1, int r1, int l2, int r2) {
        if (r1 - l1 != r2 - l2) return false;
        if (l1 == l2) return true;
        return lcp_suffix(l1, l2) >= r1 - l1;
    }
    // a[l1,r1) と a[l2,r2) の辞書順比較。小さければ -1、等しければ 0、大きければ 1。O(1)。
    //   共通接頭辞のあとに短い方が尽きた場合は「短い方が小さい」。
    int compare(int l1, int r1, int l2, int r2) {
        int n1 = r1 - l1, n2 = r2 - l2;
        int k = std::min(lcp_suffix(l1, l2), std::min(n1, n2));
        if (k < n1 && k < n2) return a[l1 + k] < a[l2 + k] ? -1 : 1;
        if (n1 == n2) return 0;
        return n1 < n2 ? -1 : 1;
    }
    // a[l,r) が周期 p を持つか (すべての i で a[i] == a[i+p] が成り立つか)。O(1)。
    //   ★p <= 0 は false、p >= 区間長 は true (自明な周期) を返す。
    bool is_periodic(int l, int r, int p) {
        if (p <= 0) return false;
        if (p >= r - l) return true;
        return lcp_suffix(l, l + p) >= (r - l) - p;
    }
    // 相異なる連続部分列 (空でないもの) の個数。n(n+1)/2 から LCP の総和を引く。
    //   ★呼ぶたびに LCP 配列を全走査するので O(N)。戻り値は __int128。
    i128 distinct_subarray_count() {
        build_suffix_array();
        i128 tot = (i128)n * (n + 1) / 2;
        for (int i = 0; i + 1 < n; i++) tot -= lcpArr[i];
        return tot;
    }
    int min_suffix() { build_suffix_array(); return sa[0]; }      // 辞書順最小の接尾辞の開始位置 O(1)
    int max_suffix() { build_suffix_array(); return sa[n - 1]; }  // 辞書順最大の接尾辞の開始位置 O(1)

    // ========================================= Manacher (回文判定 / 最長回文)
    //   区切り文字を挟んだ長さ 2n+1 の配列上で、各中心の回文半径を求める。
    //   初回呼び出し時に O(N) で構築、is_palindrome は O(1)。
    bool b_man = false;
    vector<int> manR;                             // manR[c]-1 = 中心 c の最大半径 (拡張配列上)
    void build_manacher() {
        if (b_man) return;
        b_man = true;
        int m = 2 * n + 1;
        vector<ll> t(m, LLONG_MIN);
        for (int i = 0; i < n; i++) t[2 * i + 1] = a[i];
        manR.assign(m, 0);
        int i = 0, j = 0;
        while (i < m) {
            while (i - j >= 0 && i + j < m && t[i - j] == t[i + j]) j++;
            manR[i] = j;
            int k = 1;
            while (i - k >= 0 && k + manR[i - k] < j) { manR[i + k] = manR[i - k]; k++; }
            i += k; j -= k;
        }
    }
    // a[l,r) が回文か。長さ 1 以下は常に true。O(1)。
    bool is_palindrome(int l, int r) {
        if (r - l <= 1) return true;
        build_manacher();
        return manR[l + r] >= r - l + 1;
    }
    // 最長の回文区間 [l,r) を返す。同じ長さが複数あるときは最初に見つかったもの。
    //   n == 0 なら (0,0)。★呼ぶたびに全中心を走査するので O(N)。
    pair<int, int> longest_palindrome() {
        build_manacher();
        int bl = 0, br = std::min(n, 1), bestLen = std::min(n, 1);
        for (int c = 0; c < 2 * n + 1; c++) {
            int d = manR[c] - 1;
            if (d > bestLen) { bestLen = d; bl = (c - d) / 2; br = (c + d) / 2; }
        }
        return {bl, br};
    }
    // palindrome_radius_odd(i)  : a[i] を中心とする回文の半径 k
    //                               (a[i-k .. i+k] が回文となる最大の k)。0 <= i < n。
    // palindrome_radius_even(i) : a[i-1] と a[i] の間を中心とする回文の半径 k
    //                               (a[i-k .. i+k-1] が回文となる最大の k)。★0 <= i <= n。
    // どちらも O(1)。
    int palindrome_radius_odd(int i)  { build_manacher(); return (manR[2 * i + 1] - 1 - 1) / 2; }
    int palindrome_radius_even(int i) { build_manacher(); return (manR[2 * i] - 1) / 2; }

    // ================================================ Z algorithm / prefix function
    //   どちらも初回呼び出し時に O(N) で構築し、以降は参照を返すだけ (O(1))。
    //   ★n == 0 のときも長さ 1 の配列を確保するので、空配列で呼んでも落ちない
    //     (ただし中身は意味を持たない)。
    bool b_z = false, b_pf = false;
    //   zArr[i]  = a と a[i..] の最長共通接頭辞長 (zArr[0] = n)
    //   pfArr[i] = a[0..i] の最長の「接頭辞かつ接尾辞」の長さ (KMP の失敗関数)
    vector<int> zArr, pfArr;
    void build_z() {
        if (b_z) return;
        b_z = true;
        zArr.assign(std::max(1, n), 0);
        if (n == 0) return;
        zArr[0] = n;
        int l = 0, r = 0;
        for (int i = 1; i < n; i++) {
            int z = 0;
            if (i < r) z = std::min(r - i, zArr[i - l]);
            while (i + z < n && a[z] == a[i + z]) z++;
            zArr[i] = z;
            if (i + z > r) { l = i; r = i + z; }
        }
    }
    const vector<int>& z_array() { build_z(); return zArr; }          // Z 配列への参照  O(1)
    void build_prefix_function() {
        if (b_pf) return;
        b_pf = true;
        pfArr.assign(std::max(1, n), 0);
        for (int i = 1; i < n; i++) {
            int j = pfArr[i - 1];
            while (j > 0 && a[i] != a[j]) j = pfArr[j - 1];
            if (a[i] == a[j]) j++;
            pfArr[i] = j;
        }
    }
    const vector<int>& prefix_function() { build_prefix_function(); return pfArr; }  // KMP の失敗関数 O(1)

    // ================================================ Rolling Hash (mod 2^61-1)
    //   初回呼び出し時に O(N) で構築、hash_value は O(1)。
    //   ★確率的な判定なので理論上は衝突しうる。確実に一致を判定したいときは
    //     equal() (Suffix Array + LCP) を使うこと。
    //   ★base / seed は固定シードの乱数から決めるので、同じ入力なら毎回同じ値になる。
    static constexpr ull HMOD = (1ULL << 61) - 1;
    static ull hmul(ull x, ull y) {
        u128 t = (u128)x * y;
        ull lo = (ull)(t & HMOD), hi = (ull)(t >> 61);
        ull s = lo + hi;
        return s >= HMOD ? s - HMOD : s;
    }
    static ull hadd(ull x, ull y) { ull s = x + y; return s >= HMOD ? s - HMOD : s; }
    static ull hsub(ull x, ull y) { return hadd(x, HMOD - y); }
    // 値を撹拌して [0,HMOD) に写すための関数 (splitmix64)。
    // これにより 0 が並ぶ列でも「長さの違う区間」が同じハッシュになりにくい。
    static ull mix64(ull x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    bool b_hash = false;
    ull hbase = 0, hseed = 0;                     // 基数と値の撹拌用シード (固定シードで決定的)
    vector<ull> hh, hpw;                          // hh = 累積ハッシュ, hpw = 基数の冪
    void build_hash() {
        if (b_hash) return;
        b_hash = true;
        mt19937_64 rg(998244353ULL ^ 0x9e3779b97f4a7c15ULL);
        hbase = rg() % (HMOD - 4) + 2;
        hseed = rg();
        hh.assign(n + 1, 0); hpw.assign(n + 1, 1);
        for (int i = 0; i < n; i++) {
            // 値を乱択ハッシュで [0,HMOD) に写す (0 が並ぶ列でも長さ違いが衝突しない)
            ull t = mix64((ull)a[i] ^ hseed) % HMOD;
            hh[i + 1]  = hadd(hmul(hh[i], hbase), t);
            hpw[i + 1] = hmul(hpw[i], hbase);
        }
    }
    // a[l,r) のハッシュ値。0 <= l <= r <= n (l == r でも呼べる)。O(1)。
    ull hash_value(int l, int r) { build_hash(); return hsub(hh[r], hmul(hh[l], hpw[r - l])); }

    // ============================================================ 全体の統計量
    //   ★どれも結果をキャッシュしないので、呼ぶたびに計算し直す。
    // 最長増加部分列の長さ。strict=true なら狭義単調増加、false なら非減少。
    //   ★呼ぶたびに O(N log N)。
    ll lis_length(bool strict = true) {
        vector<ll> d;
        for (int i = 0; i < n; i++) {
            auto it = strict ? lower_bound(d.begin(), d.end(), a[i]) : upper_bound(d.begin(), d.end(), a[i]);
            if (it == d.end()) d.push_back(a[i]); else *it = a[i];
        }
        return (ll)d.size();
    }
    // 配列全体の転倒数 (i < j かつ a[i] > a[j] となるペア数)。
    //   ★呼ぶたびに BIT を作り直すので O(N log N)。区間版は提供していない。
    ll inversion_count() {
        build_compress();
        int K = (int)vals.size();
        vector<int> bt(K + 1, 0);
        ll res = 0;
        for (int i = n - 1; i >= 0; i--) {
            for (int x = rnk[i]; x > 0; x -= x & -x) res += bt[x];
            for (int x = rnk[i] + 1; x <= K; x += x & -x) bt[x]++;
        }
        return res;
    }
    // 配列全体で同じ値になるペアの数 Σ C(freq, 2)。★呼ぶたびに O(K)。
    ll equal_pair_count() {
        build_compress(); ll res = 0;
        for (int v = 0; v + 1 < (int)occSt.size(); v++) { ll c = occSt[v + 1] - occSt[v]; res += c * (c - 1) / 2; }
        return res;
    }

    // ================================= 剰余に関するクエリ
    //   ★直前に使った modulus を 1 つだけキャッシュする。m を変えるたびに
    //     O(N log N) で作り直すので、複数の m を交互に使うと毎回作り直しになる。
    //   ★m > 0 が前提 (assert あり)。負の要素も数学的な剰余 (0..m-1) に正規化する。
    ll modM = 0;                                  // 現在キャッシュしている modulus
    vector<ll> modVals;                           // 実際に現れる剰余 (昇順)
    vector<vector<int>> modPos;                   // 剰余ごとの出現位置 (昇順)
    vector<vector<ll>>  modPre;                   // 剰余ごとの値の累積和
    void build_mod(ll m) {
        assert(m > 0);
        if (modM == m) return;
        modM = m;
        vector<ll> res(n);
        for (int i = 0; i < n; i++) { ll t = a[i] % m; if (t < 0) t += m; res[i] = t; }
        modVals = res;
        sort(modVals.begin(), modVals.end());
        modVals.erase(unique(modVals.begin(), modVals.end()), modVals.end());
        int K = (int)modVals.size();
        modPos.assign(K, {}); modPre.assign(K, {});
        for (int i = 0; i < n; i++) {
            int g = (int)(lower_bound(modVals.begin(), modVals.end(), res[i]) - modVals.begin());
            modPos[g].push_back(i);
        }
        for (int g = 0; g < K; g++) {
            modPre[g].assign(modPos[g].size() + 1, 0);
            for (size_t j = 0; j < modPos[g].size(); j++) modPre[g][j + 1] = modPre[g][j] + a[modPos[g][j]];
        }
    }
    // 剰余 rem が modVals の何番目かを返す (存在しなければ -1)。内部用。
    int mod_group(ll m, ll rem) {
        build_mod(m);
        ll t = rem % m; if (t < 0) t += m;
        int p = lb_idx(modVals.data(), (int)modVals.size(), t);
        return (p < (int)modVals.size() && modVals[p] == t) ? p : -1;
    }
    // [l,r) で a_i ≡ rem (mod m) を満たす個数。O(log N) (m を変えた直後は O(N log N))。
    ll count_mod(int l, int r, ll m, ll rem) {
        int g = mod_group(m, rem); if (g < 0) return 0;
        const vector<int> &P = modPos[g];
        int sz = (int)P.size();
        return (ll)(lb_idx(P.data(), sz, r) - lb_idx(P.data(), sz, l));
    }
    // [l,r) で a_i ≡ rem (mod m) を満たす要素の総和。O(log N)。
    ll sum_mod(int l, int r, ll m, ll rem) {
        int g = mod_group(m, rem); if (g < 0) return 0;
        const vector<int> &P = modPos[g];
        int sz = (int)P.size();
        int i1 = lb_idx(P.data(), sz, l), i2 = lb_idx(P.data(), sz, r);
        return modPre[g][i2] - modPre[g][i1];
    }
    ll count_divisible(int l, int r, ll d) { return count_mod(l, r, d, 0); }  // d の倍数の個数
    ll sum_divisible(int l, int r, ll d)   { return sum_mod(l, r, d, 0); }    // d の倍数の総和

    // ============================================================== 小物
    // 階差数列 b_i = a[i+1]-a[i] (長さ n-1) を新しい vector で返す。★毎回 O(N)。
    //   StaticSeq D(S.difference()); のようにすれば階差数列にも全機能を使える。
    vector<ll> difference() const {
        vector<ll> d(std::max(0, n - 1));
        for (int i = 0; i + 1 < n; i++) d[i] = a[i + 1] - a[i];
        return d;
    }
    // 2 階差分 c_i = a[i+2]-2a[i+1]+a[i] (長さ n-2) を返す。★毎回 O(N)。
    vector<ll> second_difference() const {
        vector<ll> d(std::max(0, n - 2));
        for (int i = 0; i + 2 < n; i++) d[i] = a[i + 2] - 2 * a[i + 1] + a[i];
        return d;
    }
    // [l,r) に 値が [x,y) の範囲にある要素が存在するか。x >= y なら false。O(log N)。
    bool exists_range(int l, int r, ll x, ll y) { return count_range(l, r, x, y) > 0; }
    // 「同じ値がそれ以前にも出ている」最初の位置 (重複が無ければ -1)。★毎回 O(N)。
    int first_duplicate() {
        build_occ();
        for (int i = 0; i < n; i++) if (prevOcc[i] >= 0) return i;
        return -1;
    }

    // =============================================== 全部まとめて構築したい場合
    //   ★通常のコンテスト用コードで呼ぶ必要はない。使った機能だけが自動で構築される。
    //     「入力直後にまとめて作っておきたい」「計測を安定させたい」ときだけ使う。
    //   ★呼ぶと Wavelet Matrix の累積和・区間 mex・区間線形基底など、
    //     大きい構造もすべて作られるのでメモリを大量に消費する。
    //   ★XOR 系と線形基底は全要素が 0 以上のときだけ構築される。
    void build_all() {
        build_moment(); build_sign(); build_bitcnt();
        build_rmin(); build_rmax(); build_prefix_agg();
        if (n > 0) { gcd(0, n); bit_and(0, n); bit_or(0, n); }
        build_compress(); build_wm(); build_wm_sum(); build_occ(); build_wm_prev();
        build_mex(); build_mode(); build_near(); build_lift(); build_run();
        build_signchange(); build_psrmq(); build_sst();
        build_suffix_array(); build_manacher(); build_z(); build_prefix_function(); build_hash();
        if (allNonNeg) { build_wm_bits(); build_basis(); }
    }
};
// ============================ StaticSeq ここまで ============================
//  （実戦ではここまでをコピーすれば OK。以下はテストとデモです）
// ===========================================================================

// ############################################################################
// #            ランダムテスト : 全機能を素朴解 (brute force) と突き合わせる     #
// #  ・satest::run(seed) で、いろいろな形の配列 × ランダムな区間について       #
// #    StaticSeq の結果と「その場で書いた素朴な実装」の結果を全件比較する。     #
// #  ・検査した件数は g_checks、不一致件数は g_fail に入る。                    #
// #  ・実戦では不要なので、テンプレートをコピーするときは含めなくてよい。       #
// ############################################################################
namespace satest {

static long long g_checks = 0;
static long long g_fail   = 0;
static string    g_ctx;

// 不一致が出たときに、失敗した項目名・行番号・そのときの配列と区間を表示する
static void fail(const char *msg, int line) {
    if (g_fail < 15) fprintf(stderr, "[NG] %s   (line %d)  ctx=%s\n", msg, line, g_ctx.c_str());
    g_fail++;
}
#define CK(cond, msg) do { g_checks++; if (!(cond)) fail(msg, __LINE__); } while (0)

static mt19937_64 rng(20240918ULL);
static ll rnd(ll lo, ll hi) { return lo + (ll)(rng() % (ull)(hi - lo + 1)); }

// 失敗時の表示用に配列を文字列化する
static string vecstr(const vector<ll> &a) {
    string s = "[";
    for (size_t i = 0; i < a.size(); i++) { if (i) s += ","; s += to_string(a[i]); }
    return s + "]";
}
// 浮動小数の比較 (相対誤差で判定する)
static bool close_ld(long double x, long double y, long double eps = 1e-6L) {
    long double d = fabsl(x - y), m = std::max((long double)1.0L, std::max(fabsl(x), fabsl(y)));
    return d / m <= eps;
}

// ---------------------------------------------------------------- 素朴解たち
static ll bf_gcd(const vector<ll> &a, int l, int r) { ll g = 0; for (int i = l; i < r; i++) g = ll_gcd(g, a[i]); return g; }

// 配列 a について、全体に関するクエリを 1 通り、区間クエリを iters 回ぶん検査する
static void test_array(vector<ll> a, int iters) {
    int n = (int)a.size();
    StaticSeq S(a);
    g_ctx = "n=" + to_string(n) + " a=" + (n <= 40 ? vecstr(a) : string("(long)"));

    // ======================================================= 全体に関するもの
    {
        // 座標圧縮 / 出現位置
        vector<ll> sorted_a = a; sort(sorted_a.begin(), sorted_a.end());
        CK(S.sorted() == sorted_a, "sorted");
        vector<ll> uq = sorted_a; uq.erase(unique(uq.begin(), uq.end()), uq.end());
        CK(S.unique_values() == uq, "unique_values");
        CK(S.unique_size() == (int)uq.size(), "unique_size");
        CK(S.is_all_distinct() == ((int)uq.size() == n), "is_all_distinct");
        CK(S.has_duplicate() == ((int)uq.size() != n), "has_duplicate");
        {
            map<ll, int> f; for (ll x : a) f[x]++;
            ll mx = 0, av = 0;
            for (auto &kv : f) if (kv.second > mx) { mx = kv.second; av = kv.first; }
            CK(S.max_frequency() == mx, "max_frequency");
            CK(S.count(S.argmax_frequency()) == mx, "argmax_frequency");
            (void)av;
            ll ep = 0; for (auto &kv : f) ep += (ll)kv.second * (kv.second - 1) / 2;
            CK(S.equal_pair_count() == ep, "equal_pair_count");
            vector<pair<ll,int>> ft; for (auto &kv : f) ft.push_back({kv.first, kv.second});
            CK(S.frequency_table() == ft, "frequency_table");
        }
        for (int t = 0; t < 8; t++) {
            ll x = a[rnd(0, n - 1)];
            if (t & 1) x += rnd(-2, 2);
            ll c = 0; vector<int> ps; for (int i = 0; i < n; i++) if (a[i] == x) { c++; ps.push_back(i); }
            CK(S.count(x) == c, "count(x)");
            CK(S.positions(x) == ps, "positions(x)");
            CK(S.first_pos(x) == (ps.empty() ? -1 : ps.front()), "first_pos");
            CK(S.last_pos(x) == (ps.empty() ? -1 : ps.back()), "last_pos");
            ll lb = LLONG_MAX, ub = LLONG_MAX;
            for (ll v : sorted_a) { if (v >= x && lb == LLONG_MAX) lb = v; if (v > x && ub == LLONG_MAX) ub = v; }
            CK(S.lower_bound_global(x) == lb, "lower_bound_global");
            CK(S.upper_bound_global(x) == ub, "upper_bound_global");
        }
        for (int t = 0; t < 5; t++) { ll k = rnd(0, n - 1); CK(S.kth_global(k) == sorted_a[k], "kth_global"); }
        CK(S.median_global() == sorted_a[(n - 1) / 2], "median_global");
        // 転倒数 / LIS
        if (n <= 400) {
            ll inv = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (a[i] > a[j]) inv++;
            CK(S.inversion_count() == inv, "inversion_count");
            vector<int> dp(n, 1), dp2(n, 1); int best = 0, best2 = 0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    if (a[j] <  a[i]) dp[i]  = std::max(dp[i],  dp[j] + 1);
                    if (a[j] <= a[i]) dp2[i] = std::max(dp2[i], dp2[j] + 1);
                }
                best = std::max(best, dp[i]); best2 = std::max(best2, dp2[i]);
            }
            CK(S.lis_length(true) == best, "lis_length strict");
            CK(S.lis_length(false) == best2, "lis_length nondec");
        }
        // 全体 mex
        {
            set<ll> st(a.begin(), a.end()); ll m = 0; while (st.count(m)) m++;
            CK(S.mex() == m, "global mex");
        }
        // 単調性
        CK(S.is_sorted() == is_sorted(a.begin(), a.end()), "is_sorted");
        // run length encoding
        {
            vector<pair<ll,int>> rl;
            for (int i = 0; i < n;) { int j = i; while (j < n && a[j] == a[i]) j++; rl.push_back({a[i], j - i}); i = j; }
            CK(S.runs() == rl, "runs");
        }
        // 最大長方形
        {
            i128 best = 0;
            for (int i = 0; i < n; i++) { ll mn = a[i]; for (int j = i; j < n; j++) { mn = std::min(mn, a[j]); best = std::max(best, (i128)mn * (j - i + 1)); } }
            CK(S.largest_rectangle() == best, "largest_rectangle");
        }
        // 相異なる連続部分列の数
        if (n <= 40) {
            set<vector<ll>> st;
            for (int i = 0; i < n; i++) { vector<ll> cur; for (int j = i; j < n; j++) { cur.push_back(a[j]); st.insert(cur); } }
            CK(S.distinct_subarray_count() == (i128)st.size(), "distinct_subarray_count");
        }
        // Z / prefix function
        if (n > 0) {
            vector<int> z(n, 0); z[0] = n;
            for (int i = 1; i < n; i++) { int k = 0; while (i + k < n && a[k] == a[i + k]) k++; z[i] = k; }
            CK(S.z_array() == z, "z_array");
            vector<int> pf(n, 0);
            for (int i = 1; i < n; i++) { int j = pf[i - 1]; while (j > 0 && a[i] != a[j]) j = pf[j - 1]; if (a[i] == a[j]) j++; pf[i] = j; }
            CK(S.prefix_function() == pf, "prefix_function");
        }
        // suffix array (辞書順)
        if (n > 0) {
            vector<int> idx(n); iota(idx.begin(), idx.end(), 0);
            sort(idx.begin(), idx.end(), [&](int x, int y) {
                int L = std::min(n - x, n - y);
                for (int k = 0; k < L; k++) if (a[x + k] != a[y + k]) return a[x + k] < a[y + k];
                return n - x < n - y;
            });
            CK(S.suffix_array() == idx, "suffix_array");
            for (int t = 0; t < 10; t++) {
                int i = (int)rnd(0, n - 1), j = (int)rnd(0, n - 1);
                int k = 0; while (i + k < n && j + k < n && a[i + k] == a[j + k]) k++;
                    CK(S.lcp_suffix(i, j) == k, "lcp_suffix");
                CK(S.lcp(i, j, 3) == std::min(k, 3), "lcp(maxlen)");
                int c = 0;
                if (i != j) c = (S.suffix_rank(i) < S.suffix_rank(j)) ? -1 : 1;
                int bc = 0;
                if (i != j) { if (k < n - i && k < n - j) bc = a[i + k] < a[j + k] ? -1 : 1; else bc = (n - i) < (n - j) ? -1 : 1; }
                CK(c == bc, "compare_suffix");
            }
            CK(S.min_suffix() == idx[0], "min_suffix");
            CK(S.max_suffix() == idx[n - 1], "max_suffix");
        }
        // 最長回文
        if (n > 0) {
            if (n <= 80) {
                int bestLen = 1;
                for (int i = 0; i < n; i++) for (int j = i + 1; j <= n; j++) {
                    bool pal = true; for (int k = 0; k < (j - i) / 2; k++) if (a[i + k] != a[j - 1 - k]) { pal = false; break; }
                    if (pal && j - i > bestLen) bestLen = j - i;
                }
                auto pr = S.longest_palindrome();
                CK(pr.second - pr.first == bestLen, "longest_palindrome len");
                CK(S.is_palindrome(pr.first, pr.second), "longest_palindrome is pal");
            }
            for (int i = 0; i < n; i++) {
                int k = 0; while (i - k - 1 >= 0 && i + k + 1 < n && a[i - k - 1] == a[i + k + 1]) k++;
                CK(S.palindrome_radius_odd(i) == k, "palindrome_radius_odd");
            }
            for (int i = 0; i <= n; i++) {
                int k = 0; while (i - k - 1 >= 0 && i + k < n && a[i - k - 1] == a[i + k]) k++;
                CK(S.palindrome_radius_even(i) == k, "palindrome_radius_even");
            }
        }
    }

    // ----------------------------------------------------- 基本アクセス / 小物
    {
        CK(S.size() == n && !S.empty(), "size/empty");
        CK(S.front() == a[0] && S.back() == a[n - 1], "front/back");
        CK(S.raw() == a, "raw");
        for (int t = 0; t < 5; t++) {
            int i = (int)rnd(0, n - 1);
            CK(S.get(i) == a[i] && S[i] == a[i], "get/[]");
            if (i + 1 < n) CK(S.diff(i) == a[i + 1] - a[i], "diff");
            CK(S.is_local_min(i) == ((i == 0 || a[i] < a[i - 1]) && (i == n - 1 || a[i] < a[i + 1])), "is_local_min");
            CK(S.is_local_max(i) == ((i == 0 || a[i] > a[i - 1]) && (i == n - 1 || a[i] > a[i + 1])), "is_local_max");
            int j = i; while (j + 1 < n && a[j] <= a[j + 1]) j++;
            CK(S.next_break_non_decreasing(i) == j, "next_break_non_decreasing");
            CK(S.decompress(S.compressed(i)) == a[i], "compress/decompress");
            ll x = a[i] + rnd(-1, 1);
            CK(S.rank_lb(x) == (int)(lower_bound(S.unique_values().begin(), S.unique_values().end(), x) - S.unique_values().begin()), "rank_lb");
            CK(S.rank_ub(x) == (int)(upper_bound(S.unique_values().begin(), S.unique_values().end(), x) - S.unique_values().begin()), "rank_ub");
            bool ex = false; for (ll v : a) if (v == x) ex = true;
            CK((S.value_index(x) >= 0) == ex, "value_index");
        }
        { vector<int> pk, vl; for (int i = 0; i < n; i++) { if (S.is_local_max(i)) pk.push_back(i); if (S.is_local_min(i)) vl.push_back(i); }
          CK(S.peaks() == pk && S.valleys() == vl, "peaks/valleys"); }
        { vector<ll> d1, d2;
          for (int i = 0; i + 1 < n; i++) d1.push_back(a[i + 1] - a[i]);
          for (int i = 0; i + 2 < n; i++) d2.push_back(a[i + 2] - 2 * a[i + 1] + a[i]);
          CK(S.difference() == d1 && S.second_difference() == d2, "difference"); }
        { int fd = -1;
          for (int i = 0; i < n && fd < 0; i++) for (int j = 0; j < i; j++) if (a[j] == a[i]) { fd = i; break; }
          CK(S.first_duplicate() == fd, "first_duplicate"); }
    }

    // ============================================================ 区間クエリ
    for (int it = 0; it < iters; it++) {
        int l = (int)rnd(0, n - 1), r = (int)rnd(0, n - 1);
        if (l > r) std::swap(l, r);
        r++;                                                   // [l,r), l<r
        int len = r - l;
        vector<ll> b(a.begin() + l, a.begin() + r);
        vector<ll> sb = b; sort(sb.begin(), sb.end());
        g_ctx = "n=" + to_string(n) + " [" + to_string(l) + "," + to_string(r) + ") a=" + (n <= 40 ? vecstr(a) : string("(long)"));

        // ---- 累積和系 ----
        {
            i128 s = 0, sq = 0, sia = 0, si2a = 0;
            ll xr = 0, sp = 0, sn = 0, sab = 0, spop = 0;
            ll cp = 0, cn = 0, cz = 0, co = 0;
            for (int i = l; i < r; i++) {
                i128 x = a[i];
                s += x; sq += x * x; sia += (i128)i * x; si2a += (i128)i * i * x;
                xr ^= a[i];
                if (a[i] > 0) { cp++; sp += a[i]; }
                if (a[i] < 0) { cn++; sn += a[i]; }
                if (a[i] == 0) cz++;
                if (((ull)a[i]) & 1ULL) co++;
                sab += (a[i] < 0 ? -a[i] : a[i]);
                spop += __builtin_popcountll((ull)a[i]);
            }
            CK(S.sum(l, r) == (ll)s, "sum");
            CK(S.prefix_sum(r) - S.prefix_sum(l) == (ll)s, "prefix_sum");
            CK(S.sum128(l, r) == s, "sum128");
            CK(S.sum_sq(l, r) == sq, "sum_sq");
            CK(S.sum_ia(l, r) == sia, "sum_ia");
            CK(S.sum_i2a(l, r) == si2a, "sum_i2a");
            CK(S.bit_xor(l, r) == xr, "bit_xor");
            CK(S.count_positive(l, r) == cp, "count_positive");
            CK(S.count_negative(l, r) == cn, "count_negative");
            CK(S.count_zero(l, r) == cz, "count_zero");
            CK(S.count_nonzero(l, r) == len - cz, "count_nonzero");
            CK(S.count_odd(l, r) == co, "count_odd");
            CK(S.count_even(l, r) == len - co, "count_even");
            CK(S.sum_positive(l, r) == sp, "sum_positive");
            CK(S.sum_negative(l, r) == sn, "sum_negative");
            CK(S.sum_abs(l, r) == sab, "sum_abs(l,r)");
            CK(S.sum_popcount(l, r) == spop, "sum_popcount");
            ll x = rnd(-100, 100);
            i128 sd = 0; for (int i = l; i < r; i++) sd += ((i128)a[i] - x) * ((i128)a[i] - x);
            CK(S.sum_sq_distance(l, r, x) == sd, "sum_sq_distance");
            for (int t = 0; t < 4; t++) {
                int bt = (int)rnd(0, 63); ll c = 0;
                for (int i = l; i < r; i++) c += (ll)(((ull)a[i] >> bt) & 1ULL);
                CK(S.count_bit(l, r, bt) == c, "count_bit");
            }
        }
        // ---- 統計量 ----
        {
            long double mean = 0; for (ll v : b) mean += (long double)v; mean /= len;
            CK(close_ld(S.mean(l, r), mean), "mean");
            long double var = 0; for (ll v : b) var += ((long double)v - mean) * ((long double)v - mean);
            var /= len;
            CK(close_ld(S.variance(l, r), var, 1e-5L), "variance");
            CK(close_ld(S.stddev(l, r), sqrtl(std::max((long double)0.0L, var)), 1e-5L), "stddev");
            long double mi = 0; for (int i = l; i < r; i++) mi += i; mi /= len;
            long double vi = 0, cv = 0;
            for (int i = l; i < r; i++) { vi += ((long double)i - mi) * ((long double)i - mi); cv += ((long double)i - mi) * ((long double)a[i] - mean); }
            vi /= len; cv /= len;
            CK(close_ld(S.variance_index(l, r), vi, 1e-5L), "variance_index");
            CK(close_ld(S.covariance_index(l, r), cv, 1e-5L), "covariance_index");
            auto fl = S.fit_line(l, r);
            long double slope = (vi > 0) ? cv / vi : 0.0L;
            CK(close_ld(fl.first, slope, 1e-5L), "fit_line slope");
            CK(close_ld(fl.second, mean - slope * mi, 1e-5L), "fit_line intercept");
        }
        // ---- min / max ----
        {
            ll mn = *min_element(b.begin(), b.end()), mx = *max_element(b.begin(), b.end());
            int ami = l, ama = l;
            for (int i = l; i < r; i++) { if (a[i] < a[ami]) ami = i; if (a[i] > a[ama]) ama = i; }
            CK(S.min(l, r) == mn && S.max(l, r) == mx, "min/max");
            CK(S.argmin(l, r) == ami && S.argmax(l, r) == ama, "argmin/argmax");
            CK(S.argminmax(l, r) == make_pair(ami, ama), "argminmax");
            CK(S.get(l, r) == b, "get(l,r)");
            CK(S.range_width(l, r) == mx - mn, "range_width");
            CK(S.minmax(l, r) == make_pair(mn, mx), "minmax");
            CK(S.min_with_index(l, r) == make_pair(mn, ami), "min_with_index");
            CK(S.max_with_index(l, r) == make_pair(mx, ama), "max_with_index");
            ll cmn = count(b.begin(), b.end(), mn), cmx = count(b.begin(), b.end(), mx);
            CK(S.min_count(l, r) == cmn, "min_count");
            CK(S.max_count(l, r) == cmx, "max_count");
            int lmi = l, lma = l;
            for (int i = l; i < r; i++) { if (a[i] == mn) lmi = i; if (a[i] == mx) lma = i; }
            CK(S.first_min(l, r) == ami && S.last_min(l, r) == lmi, "first/last_min");
            CK(S.first_max(l, r) == ama && S.last_max(l, r) == lma, "first/last_max");
        }
        // ---- gcd / and / or ----
        {
            ll g = bf_gcd(a, l, r), an = -1LL, orr = 0;
            for (int i = l; i < r; i++) { an &= a[i]; orr |= a[i]; }
            CK(S.gcd(l, r) == g, "gcd");
            CK(S.bit_and(l, r) == an, "bit_and");
            CK(S.bit_or(l, r) == orr, "bit_or");
        }
        // ---- prefix / suffix 集約 ----
        {
            int i = (int)rnd(0, n - 1);
            ll mn = a[0], mx = a[0], g = 0, an = -1LL, orr = 0; int ami = 0, ama = 0;
            for (int j = 0; j <= i; j++) { if (a[j] < mn) { mn = a[j]; ami = j; } if (a[j] > mx) { mx = a[j]; ama = j; } g = ll_gcd(g, a[j]); an &= a[j]; orr |= a[j]; }
            CK(S.prefix_min(i) == mn && S.prefix_max(i) == mx, "prefix_min/max");
            CK(S.prefix_argmin(i) == ami && S.prefix_argmax(i) == ama, "prefix_argmin/max");
            CK(S.prefix_gcd(i) == g && S.prefix_and(i) == an && S.prefix_or(i) == orr, "prefix_gcd/and/or");
            mn = a[i]; mx = a[i]; g = 0; an = -1LL; orr = 0; ami = i; ama = i;
            for (int j = i; j < n; j++) { if (a[j] < mn) { mn = a[j]; ami = j; } if (a[j] > mx) { mx = a[j]; ama = j; } g = ll_gcd(g, a[j]); an &= a[j]; orr |= a[j]; }
            CK(S.suffix_min(i) == mn && S.suffix_max(i) == mx, "suffix_min/max");
            CK(S.suffix_argmin(i) == ami && S.suffix_argmax(i) == ama, "suffix_argmin/max");
            CK(S.suffix_gcd(i) == g && S.suffix_and(i) == an && S.suffix_or(i) == orr, "suffix_gcd/and/or");
        }
        // ---- count / sum with value constraint ----
        for (int t = 0; t < 5; t++) {
            ll x = sb[rnd(0, len - 1)] + rnd(-2, 2);
            ll y = x + rnd(0, 20);
            ll clt = 0, cle = 0, ceq = 0, crg = 0, slt = 0, sle = 0, srg = 0;
            for (ll v : b) {
                if (v <  x) { clt++; slt += v; }
                if (v <= x) { cle++; sle += v; }
                if (v == x) ceq++;
                if (x <= v && v < y) { crg++; srg += v; }
            }
            CK(S.count_lt(l, r, x) == clt, "count_lt");
            CK(S.count_le(l, r, x) == cle, "count_le");
            CK(S.count_ge(l, r, x) == len - clt, "count_ge");
            CK(S.count_gt(l, r, x) == len - cle, "count_gt");
            CK(S.count_eq(l, r, x) == ceq, "count_eq");
            CK(S.count_range(l, r, x, y) == crg, "count_range");
            CK(S.exists_range(l, r, x, y) == (crg > 0), "exists_range");
            CK(S.sum_lt(l, r, x) == slt, "sum_lt");
            CK(S.sum_le(l, r, x) == sle, "sum_le");
            CK(S.sum_ge(l, r, x) == S.sum(l, r) - slt, "sum_ge");
            CK(S.sum_gt(l, r, x) == S.sum(l, r) - sle, "sum_gt");
            CK(S.sum_range(l, r, x, y) == srg, "sum_range");
            // predecessor / successor
            ll pv = StaticSeq::NONE, pvs = StaticSeq::NONE, nv = StaticSeq::NONE, nvs = StaticSeq::NONE;
            for (ll v : b) {
                if (v <= x && (pv  == StaticSeq::NONE || v > pv))  pv  = v;
                if (v <  x && (pvs == StaticSeq::NONE || v > pvs)) pvs = v;
                if (v >= x && (nv  == StaticSeq::NONE || v < nv))  nv  = v;
                if (v >  x && (nvs == StaticSeq::NONE || v < nvs)) nvs = v;
            }
            CK(S.prev_value(l, r, x) == pv, "prev_value");
            CK(S.prev_value_strict(l, r, x) == pvs, "prev_value_strict");
            CK(S.next_value(l, r, x) == nv, "next_value");
            CK(S.next_value_strict(l, r, x) == nvs, "next_value_strict");
            CK(S.max_le(l, r, x) == pv && S.min_ge(l, r, x) == nv, "max_le/min_ge");
            ll near = sb[0];
            for (ll v : b) { i128 d1 = (i128)v - x, d2 = (i128)near - x; if (d1 < 0) d1 = -d1; if (d2 < 0) d2 = -d2; if (d1 < d2 || (d1 == d2 && v < near)) near = v; }
            CK(S.nearest(l, r, x) == near, "nearest");
            // Σ|a_i - x|
            i128 sa2 = 0; for (ll v : b) sa2 += (i128)((v > x) ? v - x : x - v);
            CK(S.sum_abs(l, r, x) == (ll)sa2, "sum_abs(l,r,x)");
        }
        // ---- 順位統計 ----
        {
            for (int t = 0; t < 5; t++) {
                ll k = rnd(0, len - 1);
                CK(S.kth(l, r, k) == sb[k], "kth");
                CK(S.kth_largest(l, r, k) == sb[len - 1 - k], "kth_largest");
                CK(S.quantile(l, r, k) == sb[k], "quantile");
                ll ss = 0; for (ll i = 0; i <= k; i++) ss += sb[i];
                CK(S.sum_smallest(l, r, k + 1) == ss, "sum_smallest");
                ll sl = 0; for (ll i = 0; i <= k; i++) sl += sb[len - 1 - i];
                CK(S.sum_largest(l, r, k + 1) == sl, "sum_largest");
            }
            CK(S.sum_smallest(l, r, 0) == 0 && S.sum_largest(l, r, 0) == 0, "sum_smallest/largest 0");
            CK(S.sum_smallest(l, r, len) == S.sum(l, r), "sum_smallest all");
            CK(S.median(l, r) == sb[(len - 1) / 2], "median");
            int kk = (int)std::min<ll>(len, 5);
            vector<ll> smk(sb.begin(), sb.begin() + kk);
            CK(S.smallest_k(l, r, kk) == smk, "smallest_k");
            vector<ll> lgk; for (int i = 0; i < kk; i++) lgk.push_back(sb[len - 1 - i]);
            CK(S.largest_k(l, r, kk) == lgk, "largest_k");
            CK(S.second_min(l, r) == (len >= 2 ? sb[1] : StaticSeq::NONE), "second_min");
            CK(S.second_max(l, r) == (len >= 2 ? sb[len - 2] : StaticSeq::NONE), "second_max");
            ll smd = StaticSeq::NONE, smx = StaticSeq::NONE;
            for (ll v : sb) if (v > sb[0]) { smd = v; break; }
            for (int i = len - 1; i >= 0; i--) if (sb[i] < sb[len - 1]) { smx = sb[i]; break; }
            CK(S.second_min_distinct(l, r) == smd, "second_min_distinct");
            CK(S.second_max_distinct(l, r) == smx, "second_max_distinct");
            // median_cost : min_x Σ|a_i-x| は中央値で達成される
            ll mc = 0; { ll md = sb[(len - 1) / 2]; for (ll w : sb) mc += (w > md ? w - md : md - w); }
            if (len <= 60) { ll mc2 = LLONG_MAX; for (ll v : sb) { ll c = 0; for (ll w : sb) c += (w > v ? w - v : v - w); mc2 = std::min(mc2, c); } CK(mc == mc2, "median_cost(brute min)"); }
            CK(S.median_cost(l, r) == mc, "median_cost");
            // majority
            ll maj = StaticSeq::NONE;
            { map<ll,int> f; for (ll v : b) f[v]++; for (auto &kv : f) if (kv.second * 2 > len) maj = kv.first; }
            CK(S.majority(l, r) == maj, "majority");
        }
        // ---- distinct / mex / mode ----
        {
            set<ll> st(b.begin(), b.end());
            CK(S.distinct_count(l, r) == (ll)st.size(), "distinct_count");
            ll m = 0; while (st.count(m)) m++;
            CK(S.mex(l, r) == m, "mex(l,r)");
            map<ll,int> f; for (ll v : b) f[v]++;
            int best = 0; for (auto &kv : f) best = std::max(best, kv.second);
            auto md = S.mode(l, r);
            CK(md.second == best && f[md.first] == best, "mode");
        }
        // ---- 単調性 / run / 等差 ----
        {
            bool nd = true, ni = true, si = true, sd = true, cst = true, ar = true;
            for (int i = l; i + 1 < r; i++) {
                if (a[i] >  a[i + 1]) nd = false;
                if (a[i] <  a[i + 1]) ni = false;
                if (a[i] >= a[i + 1]) si = false;
                if (a[i] <= a[i + 1]) sd = false;
                if (a[i] != a[i + 1]) cst = false;
            }
            for (int i = l; i + 2 < r; i++) if (a[i + 1] - a[i] != a[i + 2] - a[i + 1]) ar = false;
            CK(S.is_non_decreasing(l, r) == nd, "is_non_decreasing");
            CK(S.is_non_increasing(l, r) == ni, "is_non_increasing");
            CK(S.is_strictly_increasing(l, r) == si, "is_strictly_increasing");
            CK(S.is_strictly_decreasing(l, r) == sd, "is_strictly_decreasing");
            CK(S.is_constant(l, r) == cst, "is_constant");
            CK(S.is_arithmetic(l, r) == ar, "is_arithmetic");
            ll sc = 0;
            for (int i = l; i + 1 < r; i++) if ((a[i] > 0 && a[i + 1] < 0) || (a[i] < 0 && a[i + 1] > 0)) sc++;
            CK(S.count_sign_changes(l, r) == sc, "count_sign_changes");
            int i = (int)rnd(0, n - 1);
            int rl = i; while (rl > 0 && a[rl - 1] == a[i]) rl--;
            int rr = i; while (rr + 1 < n && a[rr + 1] == a[i]) rr++;
            CK(S.run_left(i) == rl && S.run_right(i) == rr && S.run_length(i) == rr - rl + 1, "run_*");
        }
        // ---- nearest smaller / greater ----
        {
            int i = (int)rnd(0, n - 1);
            int ps = -1, ns = n, pse = -1, nse = n, pg = -1, ng = n, pge = -1, nge = n;
            for (int j = i - 1; j >= 0; j--) { if (a[j] <  a[i]) { ps  = j; break; } }
            for (int j = i - 1; j >= 0; j--) { if (a[j] <= a[i]) { pse = j; break; } }
            for (int j = i - 1; j >= 0; j--) { if (a[j] >  a[i]) { pg  = j; break; } }
            for (int j = i - 1; j >= 0; j--) { if (a[j] >= a[i]) { pge = j; break; } }
            for (int j = i + 1; j <  n; j++) { if (a[j] <  a[i]) { ns  = j; break; } }
            for (int j = i + 1; j <  n; j++) { if (a[j] <= a[i]) { nse = j; break; } }
            for (int j = i + 1; j <  n; j++) { if (a[j] >  a[i]) { ng  = j; break; } }
            for (int j = i + 1; j <  n; j++) { if (a[j] >= a[i]) { nge = j; break; } }
            CK(S.prev_smaller(i) == ps && S.next_smaller(i) == ns, "prev/next_smaller");
            CK(S.prev_smaller_equal(i) == pse && S.next_smaller_equal(i) == nse, "prev/next_smaller_equal");
            CK(S.prev_greater(i) == pg && S.next_greater(i) == ng, "prev/next_greater");
            CK(S.prev_greater_equal(i) == pge && S.next_greater_equal(i) == nge, "prev/next_greater_equal");
            CK(S.minimum_containing_interval(i) == make_pair(ps + 1, ns), "minimum_containing_interval");
            int par;
            if (ps < 0 && ns >= n) par = -1;
            else if (ps < 0) par = ns;
            else if (ns >= n) par = ps;
            else par = (a[ps] >= a[ns]) ? ps : ns;
            CK(S.cartesian_parent(i) == par, "cartesian_parent");
            // binary lifting
            ll kk = rnd(0, 5);
            int cur = i; for (ll t = 0; t < kk && cur < n; t++) { int nx = n; for (int j = cur + 1; j < n; j++) if (a[j] > a[cur]) { nx = j; break; } cur = nx; }
            CK(S.jump_next_greater(i, kk) == cur, "jump_next_greater");
            cur = i; for (ll t = 0; t < kk && cur < n; t++) { int nx = n; for (int j = cur + 1; j < n; j++) if (a[j] < a[cur]) { nx = j; break; } cur = nx; }
            CK(S.jump_next_smaller(i, kk) == cur, "jump_next_smaller");
        }
        // ---- 出現位置 ----
        {
            int i = (int)rnd(0, n - 1);
            int pv = -1, nx = n, fo = i;
            for (int j = i - 1; j >= 0; j--) if (a[j] == a[i]) { pv = j; break; }
            for (int j = i + 1; j < n; j++) if (a[j] == a[i]) { nx = j; break; }
            for (int j = 0; j < n; j++) if (a[j] == a[i]) { fo = j; break; }
            CK(S.prev_occurrence(i) == pv && S.next_occurrence(i) == nx && S.first_occurrence(i) == fo, "prev/next/first_occurrence");
            CK(S.prev_distance(i) == (pv < 0 ? -1 : i - pv), "prev_distance");
            CK(S.next_distance(i) == (nx >= n ? -1 : nx - i), "next_distance");
            ll x = a[rnd(0, n - 1)];
            int fe = -1, le = -1;
            for (int j = l; j < r; j++) if (a[j] == x) { if (fe < 0) fe = j; le = j; }
            CK(S.first_equal(l, r, x) == fe, "first_equal");
            CK(S.last_equal(l, r, x) == le, "last_equal");
            CK(S.exists(l, r, x) == (fe >= 0), "exists");
            ll c = 0; for (int j = l; j < r; j++) if (a[j] == x) c++;
            CK(S.count(l, r, x) == c, "count(l,r,x)");
            set<ll> st(a.begin(), a.begin() + r);
            CK(S.distinct_prefix(r) == (ll)st.size(), "distinct_prefix");
        }
        // ---- 条件を満たす最初 / 最後の位置 ----
        for (int t = 0; t < 3; t++) {
            int i = (int)rnd(0, n - 1);
            ll x = a[rnd(0, n - 1)] + rnd(-2, 2);
            int e1 = n, e2 = n, e3 = n, e4 = n;
            for (int j = i; j < n; j++) { if (a[j] >= x) { e1 = j; break; } }
            for (int j = i; j < n; j++) { if (a[j] >  x) { e2 = j; break; } }
            for (int j = i; j < n; j++) { if (a[j] <= x) { e3 = j; break; } }
            for (int j = i; j < n; j++) { if (a[j] <  x) { e4 = j; break; } }
            CK(S.next_ge(i, x) == e1 && S.next_gt(i, x) == e2, "next_ge/gt");
            CK(S.next_le(i, x) == e3 && S.next_lt(i, x) == e4, "next_le/lt");
            int f1 = -1, f2 = -1, f3 = -1, f4 = -1;
            for (int j = i; j >= 0; j--) { if (a[j] >= x) { f1 = j; break; } }
            for (int j = i; j >= 0; j--) { if (a[j] >  x) { f2 = j; break; } }
            for (int j = i; j >= 0; j--) { if (a[j] <= x) { f3 = j; break; } }
            for (int j = i; j >= 0; j--) { if (a[j] <  x) { f4 = j; break; } }
            CK(S.prev_ge(i, x) == f1 && S.prev_gt(i, x) == f2, "prev_ge/gt");
            CK(S.prev_le(i, x) == f3 && S.prev_lt(i, x) == f4, "prev_le/lt");
        }
        // ---- 累積和の検索 / 最大部分和 ----
        {
            ll best = LLONG_MIN, worst = LLONG_MAX;
            { ll cb = 0, cw = 0;                       // Kadane (実装とは別アルゴリズム)
              for (int i = l; i < r; i++) {
                  cb = std::max(a[i], cb + a[i]); best  = std::max(best, cb);
                  cw = std::min(a[i], cw + a[i]); worst = std::min(worst, cw);
              } }
            if (len <= 60) {
                ll b2 = LLONG_MIN, w2 = LLONG_MAX;
                for (int i = l; i < r; i++) { ll s = 0; for (int j = i; j < r; j++) { s += a[j]; b2 = std::max(b2, s); w2 = std::min(w2, s); } }
                CK(b2 == best && w2 == worst, "kadane sanity");
            }
            CK(S.max_subarray_sum(l, r) == best, "max_subarray_sum");
            CK(S.min_subarray_sum(l, r) == worst, "min_subarray_sum");
            ll mp = LLONG_MIN, mnp = LLONG_MAX, ms = LLONG_MIN, mns = LLONG_MAX, s = 0;
            for (int j = l; j < r; j++) { s += a[j]; mp = std::max(mp, s); mnp = std::min(mnp, s); }
            s = 0;
            for (int j = r - 1; j >= l; j--) { s += a[j]; ms = std::max(ms, s); mns = std::min(mns, s); }
            CK(S.max_prefix_sum(l, r) == mp && S.min_prefix_sum(l, r) == mnp, "max/min_prefix_sum");
            CK(S.max_suffix_sum(l, r) == ms && S.min_suffix_sum(l, r) == mns, "max/min_suffix_sum");
            ll x = rnd(-3 * (len + 1), 3 * (len + 1)) * (ll)(1 + rnd(0, 3));
            int e = -1; { ll c = 0; for (int j = l; j <= n; j++) { if (j > l) c += a[j - 1]; if (c >= x) { e = j; break; } } }
            CK(S.first_sum_ge(l, x) == e, "first_sum_ge");
            e = -1; { ll c = 0; for (int j = l; j <= n; j++) { if (j > l) c += a[j - 1]; if (c <= x) { e = j; break; } } }
            CK(S.first_sum_le(l, x) == e, "first_sum_le");
            e = -1; { for (int j = r; j >= 0; j--) { if (S.sum(j, r) >= x) { e = j; break; } } }
            CK(S.last_sum_ge(r, x) == e, "last_sum_ge");
            e = -1; { for (int j = r; j >= 0; j--) { if (S.sum(j, r) <= x) { e = j; break; } } }
            CK(S.last_sum_le(r, x) == e, "last_sum_le");
        }
        // ---- 部分列の一致 / 比較 / 周期 ----
        {
            int l2 = (int)rnd(0, n - 1), r2 = (int)rnd(0, n - 1);
            if (l2 > r2) std::swap(l2, r2);
            r2++;
            bool eq = (r - l == r2 - l2);
            if (eq) for (int k = 0; k < r - l; k++) if (a[l + k] != a[l2 + k]) { eq = false; break; }
            CK(S.equal(l, r, l2, r2) == eq, "equal");
            int cmp = 0;
            { int L = std::min(r - l, r2 - l2), k = 0;
              while (k < L && a[l + k] == a[l2 + k]) k++;
              if (k < L) cmp = a[l + k] < a[l2 + k] ? -1 : 1;
              else if (r - l != r2 - l2) cmp = (r - l < r2 - l2) ? -1 : 1; }
            CK(S.compare(l, r, l2, r2) == cmp, "compare");
            CK((S.hash_value(l, r) == S.hash_value(l2, r2)) == eq, "hash_value");
            int p = (int)rnd(1, len);
            bool per = true; for (int i = l; i + p < r; i++) if (a[i] != a[i + p]) { per = false; break; }
            CK(S.is_periodic(l, r, p) == per, "is_periodic");
            bool pal = true; for (int k = 0; k < len / 2; k++) if (a[l + k] != a[r - 1 - k]) { pal = false; break; }
            CK(S.is_palindrome(l, r) == pal, "is_palindrome(l,r)");
        }
        // ---- 剰余 ----
        {
            ll m = rnd(1, 7), rem = rnd(0, m - 1);
            ll c = 0, s = 0;
            for (int i = l; i < r; i++) { ll t = a[i] % m; if (t < 0) t += m; if (t == rem) { c++; s += a[i]; } }
            CK(S.count_mod(l, r, m, rem) == c, "count_mod");
            CK(S.sum_mod(l, r, m, rem) == s, "sum_mod");
            c = 0; s = 0;
            for (int i = l; i < r; i++) if (a[i] % m == 0) { c++; s += a[i]; }
            CK(S.count_divisible(l, r, m) == c, "count_divisible");
            CK(S.sum_divisible(l, r, m) == s, "sum_divisible");
        }
        // ---- XOR 系 / 線形基底 (全要素が 0 以上のときのみ) ----
        if (S.allNonNeg) {
            ll x = b[rnd(0, len - 1)] ^ (ll)rnd(0, 15);
            if (rnd(0, 2) == 0) x = rnd(0, 63);
            vector<ll> xs; for (ll v : b) xs.push_back(v ^ x);
            sort(xs.begin(), xs.end());
            ll bound = xs[rnd(0, len - 1)] + rnd(-2, 2);
            if (bound < 0) bound = 0;
            ll clt = 0, cle = 0;
            for (ll v : b) { if ((v ^ x) <  bound) clt++; if ((v ^ x) <= bound) cle++; }
            CK(S.count_xor_lt(l, r, x, bound) == clt, "count_xor_lt");
            CK(S.count_xor_le(l, r, x, bound) == cle, "count_xor_le");
            CK(S.count_xor_ge(l, r, x, bound) == len - clt, "count_xor_ge");
            CK(S.count_xor_gt(l, r, x, bound) == len - cle, "count_xor_gt");
            ll k = rnd(0, len - 1);
            CK(S.kth_xor(l, r, x, k) == xs[k], "kth_xor");
            CK(S.max_xor(l, r, x) == xs[len - 1], "max_xor");
            CK(S.min_xor(l, r, x) == xs[0], "min_xor");
            if (len <= 16) {
                set<ll> span; span.insert(0);
                for (ll v : b) { set<ll> nw = span; for (ll w : span) nw.insert(w ^ v); span = nw; }
                vector<ll> sp(span.begin(), span.end());
                int rk = 0; while ((1 << rk) < (int)sp.size()) rk++;
                CK(S.basis_rank(l, r) == rk, "basis_rank");
                CK(S.basis_max_xor(l, r) == sp.back(), "basis_max_xor");
                for (int t = 0; t < 6; t++) {
                    ll q = (t & 1) ? sp[rnd(0, (ll)sp.size() - 1)] ^ (ll)rnd(0, 3) : rnd(0, 80);
                    CK(S.basis_contains(l, r, q) == (span.count(q) > 0), "basis_contains");
                }
                for (int t = 0; t < 4; t++) {
                    ll kk = rnd(0, (ll)sp.size() + 1);
                    CK(S.basis_kth_xor(l, r, kk) == (kk < (ll)sp.size() ? sp[kk] : -1), "basis_kth_xor");
                }
            }
        }
    }
}

static void run(ull seed = 20240918ULL) {
    rng.seed(seed); g_checks = 0; g_fail = 0;
    struct Conf { int n; ll lo, hi; };
    vector<Conf> confs = {
        {1, 0, 0}, {1, -5, 5}, {2, 0, 1}, {3, 0, 2}, {5, 0, 1}, {5, -3, 3},
        {8, 0, 3}, {8, -8, 8}, {10, 0, 2}, {12, 0, 5}, {15, -4, 4}, {16, 0, 1},
        {20, 0, 6}, {25, -10, 10}, {30, 0, 3}, {33, 0, 100}, {40, -20, 20},
        {50, 0, 4}, {64, 0, 7}, {65, 0, 2}, {70, -50, 50}, {100, 0, 9},
        {128, 0, 3}, {129, -30, 30}, {150, 0, 20}, {200, 0, 1}, {200, -1000, 1000},
        {257, 0, 60}, {300, -5, 5}, {300, 0, 1000000},
        {500, 0, 8}, {512, 0, 300}, {513, -100, 100}, {700, 0, 2}, {1000, 0, 50},
    };
    for (size_t ci = 0; ci < confs.size(); ci++) {
        for (int rep = 0; rep < 10; rep++) {
            int n = confs[ci].n;
            vector<ll> a(n);
            for (int i = 0; i < n; i++) a[i] = rnd(confs[ci].lo, confs[ci].hi);
            test_array(a, n <= 20 ? 60 : (n <= 300 ? 40 : 20));
        }
    }
    // 特殊形状
    {
        for (int n : {1, 2, 7, 33, 100}) {
            vector<ll> inc(n), dec(n), same(n), alt(n), ar(n);
            for (int i = 0; i < n; i++) { inc[i] = i; dec[i] = n - i; same[i] = 7; alt[i] = (i & 1) ? -1 : 1; ar[i] = 3LL * i - 5; }
            test_array(inc, 20); test_array(dec, 20); test_array(same, 20);
            test_array(alt, 20); test_array(ar, 20);
            vector<ll> pal(n); for (int i = 0; i < n; i++) pal[i] = std::min(i, n - 1 - i) % 3;
            test_array(pal, 20);
        }
    }
    // build_all() を先に呼んだ場合と遅延構築の場合で結果が一致するか
    {
        for (int rep = 0; rep < 20; rep++) {
            int n = (int)rnd(1, 80);
            vector<ll> a(n);
            for (int i = 0; i < n; i++) a[i] = rnd(0, 12);
            StaticSeq A(a), B(a);
            A.build_all();
            g_ctx = "build_all n=" + to_string(n);
            for (int t = 0; t < 30; t++) {
                int l = (int)rnd(0, n - 1), r = (int)rnd(0, n - 1);
                if (l > r) std::swap(l, r);
                r++;
                ll x = rnd(0, 12), k = rnd(0, r - l - 1);
                CK(A.min(l,r)==B.min(l,r) && A.max(l,r)==B.max(l,r), "build_all min/max");
                CK(A.gcd(l,r)==B.gcd(l,r) && A.bit_and(l,r)==B.bit_and(l,r) && A.bit_or(l,r)==B.bit_or(l,r), "build_all gcd/and/or");
                CK(A.kth(l,r,k)==B.kth(l,r,k) && A.count_lt(l,r,x)==B.count_lt(l,r,x), "build_all kth/count");
                CK(A.sum_lt(l,r,x)==B.sum_lt(l,r,x) && A.distinct_count(l,r)==B.distinct_count(l,r), "build_all sum/distinct");
                CK(A.mex(l,r)==B.mex(l,r) && A.mode(l,r)==B.mode(l,r), "build_all mex/mode");
                CK(A.max_xor(l,r,x)==B.max_xor(l,r,x) && A.basis_max_xor(l,r)==B.basis_max_xor(l,r), "build_all xor/basis");
                CK(A.max_subarray_sum(l,r)==B.max_subarray_sum(l,r), "build_all maxsub");
                CK(A.is_palindrome(l,r)==B.is_palindrome(l,r) && A.hash_value(l,r)==B.hash_value(l,r), "build_all pal/hash");
            }
        }
    }
    // 大きい値 (オーバーフロー境界の確認)
    {
        for (int rep = 0; rep < 4; rep++) {
            int n = 60; vector<ll> a(n);
            for (int i = 0; i < n; i++) a[i] = rnd(-(ll)2e12, (ll)2e12);
            test_array(a, 15);
            for (int i = 0; i < n; i++) a[i] = rnd(0, (ll)4e12);
            test_array(a, 15);
        }
    }
}

} // namespace satest

// ############################################################################
// #        demo() : 全機能の使い方 (日本語解説つき) — 見れば思い出せる用       #
// #  ・小さな配列に対して各 API を実際に呼び、「式・結果・意味」を並べて表示する。#
// #  ・実戦では不要なので、テンプレートをコピーするときは含めなくてよい。       #
// ############################################################################
namespace sademo {

// 表示用の出力演算子 (__int128 / pair / vector を見やすく出す)
static ostream &operator<<(ostream &os, i128 x) { return os << i128str(x); }
template <class A, class B>
static ostream &operator<<(ostream &os, const pair<A, B> &p) { return os << "(" << p.first << ", " << p.second << ")"; }
template <class T>
static ostream &operator<<(ostream &os, const vector<T> &v) {
    os << "[";
    for (size_t i = 0; i < v.size(); i++) { if (i) os << ", "; os << v[i]; }
    return os << "]";
}
template <class T> static string to_s(const T &x) { ostringstream os; os << x; return os.str(); }
static string to_s(bool b) { return b ? "true" : "false"; }
// long long は NONE (= LLONG_MIN) だけ「該当なし」と表示する
static string to_s(ll x) { return (x == StaticSeq::NONE) ? string("NONE(該当なし)") : to_string(x); }

// 「式そのもの」と「結果」と「説明」を並べて表示する
//   expr はマクロ内で 1 度だけ評価され、#expr で式の文字列がそのまま表示される
#define SHOW(note, expr) do { \
    string v = sademo::to_s(expr); \
    printf("    %-44s = %s\n", #expr, v.c_str()); \
    printf("        %s\n", note); \
} while (0)
#define NOTE(s) printf("    %s\n", s)
static void sec(const char *t) { printf("\n============================================================\n %s\n============================================================\n", t); }

static void demo() {
    printf("################################################################\n");
    printf("#  StaticSeq デモ : 全機能の使い方と意味 (半開区間 [l,r) 固定)  #\n");
    printf("################################################################\n");

    vector<ll> a = {3, 1, 4, 1, 5, 9, 2, 6, 5, 3};
    StaticSeq S(a);
    printf("\n  対象配列 a = %s   (index:   0  1  2  3  4  5  6  7  8  9)\n", to_s(a).c_str());
    NOTE("コンストラクタ : StaticSeq S(vector<long long>);  構築は O(N)。");
    NOTE("重い構造は「その機能を最初に呼んだ時」に自動で作られます (遅延構築)。");

    // ------------------------------------------------------------------ 基本
    sec("1. 基本アクセス");
    SHOW("要素数", S.size());
    SHOW("空かどうか", S.empty());
    SHOW("a[4] を取得 (O(1))", S.get(4));
    SHOW("operator[] でも同じ", S[4]);
    SHOW("先頭要素", S.front());
    SHOW("末尾要素", S.back());
    SHOW("部分配列を vector で取り出す (O(r-l))", S.get(2, 6));
    SHOW("a[i+1]-a[i] (階差)", S.diff(2));

    // -------------------------------------------------------------- 累積和系
    sec("2. 累積和系 (すべて O(1))");
    SHOW("[2,6) の総和 = 4+1+5+9", S.sum(2, 6));
    SHOW("[0,r) の総和", S.prefix_sum(6));
    SHOW("総和を __int128 で (桁あふれ対策)", S.sum128(0, 10));
    SHOW("Σ a_i^2  (__int128)", S.sum_sq(2, 6));
    SHOW("Σ i*a_i  (__int128) 重み付き累積和", S.sum_ia(2, 6));
    SHOW("Σ i^2*a_i (__int128)", S.sum_i2a(2, 6));
    SHOW("Σ (a_i - x)^2   … x=4 との二乗距離", S.sum_sq_distance(2, 6, 4));
    SHOW("[2,6) の XOR", S.bit_xor(2, 6));
    SHOW("正の要素の個数", S.count_positive(0, 10));
    SHOW("負の要素の個数", S.count_negative(0, 10));
    SHOW("0 の個数", S.count_zero(0, 10));
    SHOW("0 でない要素の個数", S.count_nonzero(0, 10));
    SHOW("奇数の個数", S.count_odd(0, 10));
    SHOW("偶数の個数", S.count_even(0, 10));
    SHOW("正の要素だけの和", S.sum_positive(0, 10));
    SHOW("負の要素だけの和", S.sum_negative(0, 10));
    SHOW("Σ |a_i|", S.sum_abs(0, 10));
    SHOW("Σ popcount(a_i)", S.sum_popcount(0, 10));
    SHOW("bit 0 が立っている要素の個数", S.count_bit(0, 10, 0));
    SHOW("bit 2 が立っている要素の個数", S.count_bit(0, 10, 2));

    // ------------------------------------------------------------ 統計量
    sec("3. 区間の統計量 (すべて O(1))");
    printf("    (小数は long double で返ります)\n");
    SHOW("平均", (double)S.mean(0, 10));
    SHOW("分散 (母分散)", (double)S.variance(0, 10));
    SHOW("標準偏差", (double)S.stddev(0, 10));
    SHOW("index の分散", (double)S.variance_index(0, 10));
    SHOW("index と値の共分散", (double)S.covariance_index(0, 10));
    SHOW("index との相関係数", (double)S.correlation_index(0, 10));
    NOTE("S.fit_line(l,r) → 最小二乗直線 y = slope*i + intercept の (slope, intercept)");
    { auto f = S.fit_line(0, 10); printf("        S.fit_line(0,10) = (slope=%.6f, intercept=%.6f)\n", (double)f.first, (double)f.second); }

    // ------------------------------------------------------------ min / max
    sec("4. 区間 min / max  (Sparse Table 相当, 構築 O(N) / クエリ O(1))");
    SHOW("最小値", S.min(2, 8));
    SHOW("最大値", S.max(2, 8));
    SHOW("最小値の位置 (同値なら最も左)", S.argmin(2, 8));
    SHOW("最大値の位置 (同値なら最も左)", S.argmax(2, 8));
    SHOW("(min, max) を同時に", S.minmax(2, 8));
    SHOW("(argmin, argmax) を同時に", S.argminmax(2, 8));
    SHOW("(最小値, その位置)", S.min_with_index(2, 8));
    SHOW("(最大値, その位置)", S.max_with_index(2, 8));
    SHOW("max - min", S.range_width(2, 8));
    SHOW("最小値の個数 (O(log N))", S.min_count(0, 10));
    SHOW("最大値の個数 (O(log N))", S.max_count(0, 10));
    SHOW("最小値が最初に現れる位置", S.first_min(0, 10));
    SHOW("最小値が最後に現れる位置", S.last_min(0, 10));
    SHOW("最大値が最初に現れる位置", S.first_max(0, 10));
    SHOW("最大値が最後に現れる位置", S.last_max(0, 10));

    // ------------------------------------------------------- gcd / bit 集約
    sec("5. 区間 GCD / AND / OR  (Sparse Table, クエリ O(1))");
    SHOW("区間 GCD (負の値は絶対値で扱う)", S.gcd(0, 10));
    SHOW("区間 AND", S.bit_and(0, 10));
    SHOW("区間 OR", S.bit_or(0, 10));
    NOTE("XOR は prefix XOR なので S.bit_xor(l,r) (2 章) を使ってください。");

    // ------------------------------------------------- prefix / suffix 集約
    sec("6. prefix / suffix の集約 (すべて O(1))");
    SHOW("a[0..4] の最小値", S.prefix_min(4));
    SHOW("a[0..4] の最大値", S.prefix_max(4));
    SHOW("a[0..4] の最小値の位置", S.prefix_argmin(4));
    SHOW("a[0..4] の最大値の位置", S.prefix_argmax(4));
    SHOW("a[4..n-1] の最小値", S.suffix_min(4));
    SHOW("a[4..n-1] の最大値", S.suffix_max(4));
    SHOW("a[4..n-1] の最小値の位置", S.suffix_argmin(4));
    SHOW("a[4..n-1] の最大値の位置", S.suffix_argmax(4));
    SHOW("a[0..4] の GCD", S.prefix_gcd(4));
    SHOW("a[4..n-1] の GCD", S.suffix_gcd(4));
    SHOW("a[0..4] の AND", S.prefix_and(4));
    SHOW("a[4..n-1] の AND", S.suffix_and(4));
    SHOW("a[0..4] の OR", S.prefix_or(4));
    SHOW("a[4..n-1] の OR", S.suffix_or(4));

    // -------------------------------------------------------- 個数 / 和 (値条件)
    sec("7. 値の条件つき 個数 / 総和  (Wavelet Matrix, クエリ O(log N))");
    SHOW("a_i <  4 の個数", S.count_lt(0, 10, 4));
    SHOW("a_i <= 4 の個数", S.count_le(0, 10, 4));
    SHOW("a_i >= 4 の個数", S.count_ge(0, 10, 4));
    SHOW("a_i >  4 の個数", S.count_gt(0, 10, 4));
    SHOW("a_i == 3 の個数", S.count_eq(0, 10, 3));
    SHOW("3 <= a_i < 6 の個数 (矩形カウント)", S.count_range(0, 10, 3, 6));
    SHOW("a_i <  4 の総和", S.sum_lt(0, 10, 4));
    SHOW("a_i <= 4 の総和", S.sum_le(0, 10, 4));
    SHOW("a_i >= 4 の総和", S.sum_ge(0, 10, 4));
    SHOW("a_i >  4 の総和", S.sum_gt(0, 10, 4));
    SHOW("3 <= a_i < 6 の総和 (矩形和)", S.sum_range(0, 10, 3, 6));

    // --------------------------------------------------------------- 順位統計
    sec("8. 区間の順位統計 (Wavelet Matrix, クエリ O(log N))");
    SHOW("0-indexed で 3 番目に小さい値", S.kth(0, 10, 3));
    SHOW("0-indexed で 0 番目に大きい値 = 最大", S.kth_largest(0, 10, 0));
    SHOW("中央値 (偶数長なら小さい方)", S.median(0, 10));
    SHOW("quantile は kth の別名", S.quantile(0, 10, 5));
    SHOW("小さい方から 3 個の和", S.sum_smallest(0, 10, 3));
    SHOW("大きい方から 3 個の和", S.sum_largest(0, 10, 3));
    SHOW("小さい方から 4 個を昇順で列挙 (O(k log N))", S.smallest_k(0, 10, 4));
    SHOW("大きい方から 4 個を降順で列挙", S.largest_k(0, 10, 4));
    SHOW("2 番目に小さい値 (重複込み)", S.second_min(0, 10));
    SHOW("2 番目に大きい値 (重複込み)", S.second_max(0, 10));
    SHOW("2 番目に小さい「相異なる」値", S.second_min_distinct(0, 10));
    SHOW("2 番目に大きい「相異なる」値", S.second_max_distinct(0, 10));

    // ------------------------------------------------ predecessor / successor
    sec("9. 区間 predecessor / successor (クエリ O(log N))");
    SHOW("x=4 以下で最大の値", S.prev_value(0, 10, 4));
    SHOW("x=4 未満で最大の値", S.prev_value_strict(0, 10, 4));
    SHOW("x=4 以上で最小の値", S.next_value(0, 10, 4));
    SHOW("x=4 より大きい最小の値", S.next_value_strict(0, 10, 4));
    SHOW("max_le は prev_value の別名", S.max_le(0, 10, 4));
    SHOW("min_ge は next_value の別名", S.min_ge(0, 10, 4));
    SHOW("x=4 に最も近い値 (同距離なら小さい方)", S.nearest(0, 10, 4));
    SHOW("Σ|a_i - 4|", S.sum_abs(0, 10, 4));
    SHOW("min_x Σ|a_i - x| (中央値コスト)", S.median_cost(0, 10));
    SHOW("過半数を占める値 (無ければ NONE)", S.majority(0, 10));
    {
        vector<ll> m = {7, 7, 3, 7, 7, 1, 7};
        StaticSeq M(m);
        printf("    別の配列 m = %s\n", to_s(m).c_str());
        SHOW("過半数 (5 個中 4 個が 7)", M.majority(0, 5));
        SHOW("過半数なし の例", M.majority(2, 4));
    }

    // ------------------------------------------------------- distinct / mex / mode
    sec("10. 静的配列ならではの強力クエリ");
    SHOW("区間内の相異なる値の個数 (O(log N))", S.distinct_count(0, 10));
    SHOW("[3,8) の相異なる値の個数", S.distinct_count(3, 8));
    SHOW("prefix の distinct 数", S.distinct_prefix(5));
    {
        vector<ll> f = {0, 1, 3, 0, 2, 1, 4, 0};
        StaticSeq F(f);
        printf("    別の配列 f = %s   (index: 0 1 2 3 4 5 6 7)\n", to_s(f).c_str());
        SHOW("区間 mex : 含まれない最小の非負整数 (O(log N))", F.mex(0, 8));
        SHOW("[0,3) = {0,1,3} の mex", F.mex(0, 3));
        SHOW("[2,6) = {3,0,2,1} の mex", F.mex(2, 6));
        SHOW("配列全体の mex", F.mex());
    }
    NOTE("S.mode(l,r) → (最頻値, その出現回数)。前処理 O(N sqrt N) / クエリ O(sqrt N) / 同点は任意の 1 つ");
    SHOW("区間の最頻値", S.mode(0, 10));
    SHOW("[3,9) の最頻値", S.mode(3, 9));

    // ------------------------------------------------------------ 出現位置
    sec("11. 値と出現位置");
    SHOW("値 3 の全体での出現回数 (O(log N))", S.count(3));
    SHOW("値 3 の出現位置すべて", S.positions(3));
    SHOW("値 3 が最初に現れる位置", S.first_pos(3));
    SHOW("値 3 が最後に現れる位置", S.last_pos(3));
    SHOW("[2,9) に含まれる値 5 の個数", S.count(2, 9, 5));
    SHOW("[2,9) に値 5 があるか", S.exists(2, 9, 5));
    SHOW("[2,9) で値 5 が最初に現れる位置 (無ければ -1)", S.first_equal(2, 9, 5));
    SHOW("[2,9) で値 5 が最後に現れる位置", S.last_equal(2, 9, 5));
    SHOW("a[8] と同じ値の直前の位置 (無ければ -1)", S.prev_occurrence(8));
    SHOW("a[1] と同じ値の直後の位置 (無ければ n)", S.next_occurrence(1));
    SHOW("a[8] と同じ値が最初に現れる位置", S.first_occurrence(8));
    SHOW("直前の同値までの距離 (無ければ -1)", S.prev_distance(8));
    SHOW("直後の同値までの距離 (無ければ -1)", S.next_distance(1));

    // --------------------------------------------------------- 座標圧縮 / 全体
    sec("12. 座標圧縮 と 配列全体の情報");
    SHOW("ソート済みコピー", S.sorted());
    SHOW("相異なる値 (昇順)", S.unique_values());
    SHOW("相異なる値の個数", S.unique_size());
    SHOW("a[5] の順位 (座標圧縮後の値)", S.compressed(5));
    SHOW("順位 → 元の値", S.decompress(S.compressed(5)));
    SHOW("値 4 未満のユニーク値の個数 (= 4 の下位順位)", S.rank_lb(4));
    SHOW("値 4 以下のユニーク値の個数", S.rank_ub(4));
    SHOW("値 4 の順位 (存在しなければ -1)", S.value_index(4));
    SHOW("全体をソートしたときの 0-indexed 3 番目", S.kth_global(3));
    SHOW("全体の中央値", S.median_global());
    SHOW("4 以上で最小の要素 (無ければ LLONG_MAX)", S.lower_bound_global(4));
    SHOW("4 より大きい最小の要素", S.upper_bound_global(4));
    SHOW("全要素が相異なるか", S.is_all_distinct());
    SHOW("重複があるか", S.has_duplicate());
    SHOW("最大出現回数", S.max_frequency());
    SHOW("最大出現回数をとる値", S.argmax_frequency());
    SHOW("(値, 出現回数) の一覧", S.frequency_table());
    SHOW("同じ値のペア数 Σ C(freq,2)", S.equal_pair_count());
    SHOW("転倒数 (O(N log N))", S.inversion_count());
    SHOW("最長増加部分列の長さ (狭義単調増加)", S.lis_length(true));
    SHOW("最長非減少部分列の長さ", S.lis_length(false));
    SHOW("2 回目の出現が最初に起きる位置 (重複なしなら -1)", S.first_duplicate());
    SHOW("値が [3,6) の要素が [2,8) に存在するか", S.exists_range(2, 8, 3, 6));
    SHOW("階差数列 b_i = a_{i+1}-a_i", S.difference());
    SHOW("2 階差分 c_i = a_{i+2}-2a_{i+1}+a_i", S.second_difference());
    NOTE("→ StaticSeq D(S.difference()); のように、階差数列に対しても全機能が使えます。");

    // ---------------------------------------------------- 直近の大小関係
    sec("13. 直近のより小さい / 大きい要素 (単調スタック前計算, クエリ O(1))");
    NOTE("見つからない場合 : 左向きは -1, 右向きは n を返します。");
    SHOW("a[j] <  a[5] となる最大の j<5", S.prev_smaller(5));
    SHOW("a[j] <  a[5] となる最小の j>5", S.next_smaller(5));
    SHOW("a[j] <= a[5] となる最大の j<5", S.prev_smaller_equal(5));
    SHOW("a[j] <= a[5] となる最小の j>5", S.next_smaller_equal(5));
    SHOW("a[j] >  a[2] となる最大の j<2", S.prev_greater(2));
    SHOW("a[j] >  a[2] となる最小の j>2", S.next_greater(2));
    SHOW("a[j] >= a[2] となる最大の j<2", S.prev_greater_equal(2));
    SHOW("a[j] >= a[2] となる最小の j>2", S.next_greater_equal(2));
    SHOW("最小値 Cartesian Tree の親 (根は -1)", S.cartesian_parent(5));
    SHOW("a[5] が最小値になる最大の区間 [l,r)", S.minimum_containing_interval(5));
    SHOW("ヒストグラム最大長方形 (非負前提, __int128)", S.largest_rectangle());
    SHOW("next_greater を 2 回たどった位置 (binary lifting)", S.jump_next_greater(0, 2));
    SHOW("next_smaller を 2 回たどった位置", S.jump_next_smaller(5, 2));

    // ------------------------------------------- しきい値で最初 / 最後を探す
    sec("14. しきい値検索 : i 以降 / i 以前で最初に条件を満たす位置 (O(log N))");
    SHOW("j>=2 で a[j] >= 6 となる最小の j (無ければ n)", S.next_ge(2, 6));
    SHOW("j>=2 で a[j] >  6 となる最小の j", S.next_gt(2, 6));
    SHOW("j>=2 で a[j] <= 2 となる最小の j", S.next_le(2, 2));
    SHOW("j>=2 で a[j] <  2 となる最小の j", S.next_lt(2, 2));
    SHOW("j<=7 で a[j] >= 6 となる最大の j (無ければ -1)", S.prev_ge(7, 6));
    SHOW("j<=7 で a[j] >  6 となる最大の j", S.prev_gt(7, 6));
    SHOW("j<=7 で a[j] <= 2 となる最大の j", S.prev_le(7, 2));
    SHOW("j<=7 で a[j] <  2 となる最大の j", S.prev_lt(7, 2));

    // ------------------------------------------------ 累積和検索 / 最大部分和
    sec("15. 累積和に対する検索 と 区間の最大部分和 (O(log N))");
    SHOW("sum(2,r) >= 10 となる最小の r (無ければ -1)", S.first_sum_ge(2, 10));
    SHOW("sum(2,r) <= 5  となる最小の r", S.first_sum_le(2, 5));
    SHOW("sum(l,8) >= 10 となる最大の l", S.last_sum_ge(8, 10));
    SHOW("sum(l,8) <= 5  となる最大の l", S.last_sum_le(8, 5));
    SHOW("区間内の連続部分列の和の最大値 (非空)", S.max_subarray_sum(0, 10));
    SHOW("区間内の連続部分列の和の最小値 (非空)", S.min_subarray_sum(0, 10));
    SHOW("左端を l に固定した和の最大値", S.max_prefix_sum(2, 8));
    SHOW("左端を l に固定した和の最小値", S.min_prefix_sum(2, 8));
    SHOW("右端を r に固定した和の最大値", S.max_suffix_sum(2, 8));
    SHOW("右端を r に固定した和の最小値", S.min_suffix_sum(2, 8));

    // ------------------------------------------------------------- 単調性 / run
    sec("16. 単調性 / 等値 run / 等差 (すべて O(1))");
    {
        vector<ll> c = {2, 2, 5, 7, 7, 9, 4, 4, 4, 1};
        StaticSeq T(c);
        printf("    別の配列 c = %s\n", to_s(c).c_str());
        SHOW("[0,6) は非減少か", T.is_non_decreasing(0, 6));
        SHOW("[5,10) は非増加か", T.is_non_increasing(5, 10));
        SHOW("[2,4) は狭義単調増加か", T.is_strictly_increasing(2, 4));
        SHOW("[5,7) は狭義単調減少か", T.is_strictly_decreasing(5, 7));
        SHOW("[6,9) は全部同じ値か", T.is_constant(6, 9));
        SHOW("[0,10) 全体がソート済みか", T.is_sorted());
        SHOW("i=0 から右にどこまで非減少か (その最大 index)", T.next_break_non_decreasing(0));
        SHOW("i=7 を含む等値 run の左端", T.run_left(7));
        SHOW("i=7 を含む等値 run の右端", T.run_right(7));
        SHOW("i=7 を含む等値 run の長さ", T.run_length(7));
        SHOW("ランレングス圧縮 (値, 長さ)", T.runs());
        SHOW("i=3 は極大か", T.is_local_max(3));
        SHOW("i=9 は極小か", T.is_local_min(9));
        SHOW("極大の位置一覧", T.peaks());
        SHOW("極小の位置一覧", T.valleys());
        vector<ll> d = {5, 8, 11, 14, 2, 9};
        StaticSeq U(d);
        printf("    別の配列 d = %s\n", to_s(d).c_str());
        SHOW("[0,4) は等差数列か", U.is_arithmetic(0, 4));
        SHOW("[0,5) は等差数列か", U.is_arithmetic(0, 5));
        vector<ll> e = {3, -1, 4, -1, -5, 9};
        StaticSeq V(e);
        printf("    別の配列 e = %s\n", to_s(e).c_str());
        SHOW("隣接要素の符号が反転する回数", V.count_sign_changes(0, 6));
    }

    // ------------------------------------------------------- 部分列の比較
    sec("17. 部分列の一致 / 辞書順比較 / 周期 / 回文 (Suffix Array + Manacher)");
    {
        vector<ll> c = {1, 2, 1, 2, 1, 2, 3, 1, 2, 1};
        StaticSeq T(c);
        printf("    別の配列 c = %s\n", to_s(c).c_str());
        SHOW("c[0,3) と c[2,5) は等しいか (O(1))", T.equal(0, 3, 2, 5));
        SHOW("c[0,3) と c[6,9) の辞書順比較 (-1/0/1)", T.compare(0, 3, 6, 9));
        SHOW("接尾辞 c[0..] と c[2..] の最長共通接頭辞 (O(1))", T.lcp_suffix(0, 2));
        SHOW("上限つき lcp", T.lcp(0, 2, 3));
        SHOW("接尾辞 c[3..] の辞書順の順位", T.suffix_rank(3));
        SHOW("接尾辞配列", T.suffix_array());
        SHOW("辞書順最小の接尾辞の開始位置", T.min_suffix());
        SHOW("辞書順最大の接尾辞の開始位置", T.max_suffix());
        SHOW("相異なる連続部分列の個数 (__int128)", T.distinct_subarray_count());
        SHOW("c[0,6) は周期 2 を持つか (O(1))", T.is_periodic(0, 6, 2));
        SHOW("c[0,6) は周期 3 を持つか", T.is_periodic(0, 6, 3));
        SHOW("c[6,10) は回文か (O(1))", T.is_palindrome(6, 10));
        SHOW("c[2,5) は回文か", T.is_palindrome(2, 5));
        SHOW("最長回文区間 [l,r)", T.longest_palindrome());
        SHOW("i=4 を中心とする奇数長回文の半径", T.palindrome_radius_odd(4));
        SHOW("i=5 の直前を中心とする偶数長回文の半径", T.palindrome_radius_even(5));
        SHOW("Z algorithm", T.z_array());
        SHOW("prefix function (KMP)", T.prefix_function());
        NOTE("S.hash_value(l,r) : 区間の Rolling Hash (mod 2^61-1, O(1))。");
        printf("        T.hash_value(0,3)==T.hash_value(2,5) → %s\n",
               to_s(T.hash_value(0, 3) == T.hash_value(2, 5)).c_str());
    }

    // ------------------------------------------------------------- XOR / 基底
    sec("18. XOR 系クエリ と 区間の線形基底 (要素が 0 以上のときのみ)");
    SHOW("(a_i ^ 3) < 5 を満たす個数 (O(log A))", S.count_xor_lt(0, 10, 3, 5));
    SHOW("(a_i ^ 3) <= 5 を満たす個数", S.count_xor_le(0, 10, 3, 5));
    SHOW("(a_i ^ 3) >= 5 を満たす個数", S.count_xor_ge(0, 10, 3, 5));
    SHOW("(a_i ^ 3) >  5 を満たす個数", S.count_xor_gt(0, 10, 3, 5));
    SHOW("(a_i ^ 3) の昇順 2 番目の値", S.kth_xor(0, 10, 3, 2));
    SHOW("max (a_i ^ 3)", S.max_xor(0, 10, 3));
    SHOW("min (a_i ^ 3)", S.min_xor(0, 10, 3));
    NOTE("以下は「区間の要素から選ぶ部分集合の XOR」に関するクエリです。");
    SHOW("区間の線形基底の次元 (張る空間は 2^次元 個)", S.basis_rank(0, 10));
    SHOW("部分集合 XOR の最大値", S.basis_max_xor(0, 10));
    SHOW("7 を部分集合 XOR で作れるか", S.basis_contains(0, 10, 7));
    SHOW("部分集合 XOR の昇順 3 番目 (0-indexed, 空集合=0 込み)", S.basis_kth_xor(0, 10, 3));

    // ------------------------------------------------------------------ 剰余
    sec("19. 剰余に関するクエリ (modulus ごとに O(N) でキャッシュ, クエリ O(log N))");
    SHOW("a_i ≡ 1 (mod 3) の個数", S.count_mod(0, 10, 3, 1));
    SHOW("a_i ≡ 1 (mod 3) の総和", S.sum_mod(0, 10, 3, 1));
    SHOW("3 の倍数の個数", S.count_divisible(0, 10, 3));
    SHOW("3 の倍数の総和", S.sum_divisible(0, 10, 3));

    // ------------------------------------------------------------------ 補足
    sec("20. 補足");
    NOTE("・S.build_all() を呼ぶと全構造を先に作れます (入力直後にまとめて作りたい時用)。");
    NOTE("・NONE (= LLONG_MIN) は prev_value / next_value / majority などの「該当なし」。");
    NOTE("・XOR 系 / 線形基底は全要素が 0 以上のときのみ (assert で守っています)。");
    NOTE("・sum 系の戻り値は long long。桁あふれが心配なら sum128 / sum_sq を使ってください。");
    NOTE("・区間はすべて半開区間 [l,r) で、l < r を仮定しています (distinct_count のみ l>=r で 0)。");
    printf("\n");
}

} // namespace sademo

void solv(){
    string mode = "demo";

    bool doDemo = (mode.empty() || mode == "demo");
    int seeds = 100;
    if (!mode.empty() && mode != "demo") seeds = 1;

    if (doDemo) sademo::demo();

    if (seeds > 0) {
        printf("\n============================================================\n");
        printf(" ランダムテスト : 全機能を素朴解 (brute force) と突き合わせます\n");
        printf("============================================================\n");
        long long tot = 0, bad = 0;
        auto st = chrono::steady_clock::now();
        for (int s = 1; s <= seeds; s++) {
            satest::run(1000003ULL * (ull)s + 7);
            tot += satest::g_checks; bad += satest::g_fail;
            printf("  seed %2d : 検査 %8lld 件 / 不一致 %lld 件\n", s, satest::g_checks, satest::g_fail);
            fflush(stdout);
        }
        double sec = chrono::duration<double>(chrono::steady_clock::now() - st).count();
        printf("  ------------------------------------------------------------\n");
        printf("  合計 : 検査 %lld 件 / 不一致 %lld 件   (%.2f 秒)\n", tot, bad, sec);
        printf("  %s\n", bad == 0 ? "=> 全機能が素朴解と完全に一致しました。" : "=> 不一致があります！");
        if (bad) return;
    }
    return;
}

int main(){
    BUFF;
    int T=1;
    //cin>>T;
    while(T--)solv();
    return 0;
}