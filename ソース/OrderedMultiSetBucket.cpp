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
// 【設計見直しメモ】ベースを B+Tree → 平方分割 (バケット) + Fenwick Tree に変更
// ----------------------------------------------------------------------------
// ■ なぜ変えたのか
//   旧版は B+Tree (B=64) だった。B+Tree は「根から葉まで 3〜4 段、各段で
//   ノードを読んで二分探索して子へ降りる」ため、1 クエリあたり
//   「前のロード結果が出るまで次のロードを発行できない」依存ロードの鎖が
//   3〜4 本つながる。各段でキャッシュミスすると、その待ち時間が丸ごと直列に積まれる。
//
//   一方、平方分割は「上段 1 枚 + 下段 1 枚」の 2 段しかない。
//     ・上段 (バケット代表値の配列 bmax / Fenwick) は要素数 N/BK 個ぶんしかなく、
//       N=1e6, BK=128 なら 8000 要素 = 64KB 程度。L1/L2 に丸ごと乗る。
//     ・下段 (バケット本体) だけがコールドなアクセスになる。
//   つまり「コールドな探索は 1 回だけ」になり、依存ロードの鎖が短くなる。
//
// ■ 実測 (同一マシン・同一ワークロード、単位 ms。3 回実行の中央値。合計が小さいほど速い)
//   ワークロード: N 回 insert → count_lt / kth / contains / sum_lt / next を各 N 回
//                 → erase_one を N/2 回 → insert と erase_one を混在させて N 回
//
//     N       旧 B+Tree   新バケット    倍率
//     2e5        247         200       1.24x
//     5e5        741         602       1.23x
//     1e6       1783        1395       1.28x
//     2e6       4991        3281       1.52x
//
//   操作別 (N=1e6 / N=2e6)
//     insert              219 / 630ms  →  192 / 454ms   (1.14x / 1.39x)
//     count_lt (= rank)   254 / 809ms  →  192 / 452ms   (1.32x / 1.79x)
//     contains            221 / 623ms  →  183 / 424ms   (1.21x / 1.47x)
//     next (= prev)       222 / 573ms  →  176 / 417ms   (1.26x / 1.37x)
//     kth                 115 / 382ms  →   84 / 210ms   (1.37x / 1.82x)
//     sum_lt              270 / 784ms  →  226 / 531ms   (1.19x / 1.48x)
//     erase_one           170 / 477ms  →  111 / 265ms   (1.53x / 1.80x)
//     insert/erase 混在   296 / 732ms  →  233 / 535ms   (1.27x / 1.37x)
//   N が大きいほど差が開く (B+Tree は段数が増えて依存ロードの鎖が伸びるため)。
//   min()/max() は O(log N) → O(1)、erase_all はバケットごと落とせるので更に速い。
//   メモリも減る (N=1e6 で概算 30MB 超 → 16MB 程度)。
//
// ■ 「B+Tree をフラット化すれば勝てるのでは？」も検証した
//   旧版の遅さの一因は、各ノードが std::vector を 4 本持っていて
//   「ノード → vector のヒープ領域」という余計な間接参照が段ごとに入ることだった。
//   そこでノードを固定長配列にしてプールに詰めた「フラット B+Tree」を作り、
//   静的構造 (更新なし・理想的な配置) でクエリだけ比較した:
//     N=1e6  flatB+Tree<32>  803ms
//            flatB+Tree<64>  770ms
//            flatB+Tree<128> 803ms
//            bucket<128>     643ms   ← バケットの方が速い
//   段数そのものが減らないため、フラット化しても平方分割には届かなかった。
//   よってベースの変更が正しい判断と結論した。
//
// ■ Fenwick は「遅延再構築するキャッシュ」にしている
//   バケット数が変わる操作 (分割 / 併合 / 空バケット削除) が起きると Fenwick 全体の
//   添字がずれるため、素朴には O(nb) で作り直すことになる。分割は BK 回の挿入ごとに
//   起きるので、これが insert コストの 3 割を占めていた。
//   そこで「Fenwick は bk / bsum から復元できるキャッシュ」と割り切り、
//   構造が変わったら fw_dirty を立てるだけにして、実際に Fenwick を読むクエリの
//   直前でまとめて作り直すようにした。連続した挿入・削除の間は 1 度も作り直さない。
//   これで insert は旧 B+Tree より速くなった (219ms → 192ms at N=1e6)。
//
// ■ どこで損をしているか (正直な注意点)
//   ・insert の伸びは他より小さい (1.14x)。バケット内 memmove (平均 BK/2 要素) が
//     支配的で、そこは B+Tree の葉 (平均 B/2=32 要素) より移動量が多いため。
//     なお素朴に作ると「分割のたびに Fenwick を O(nb) で作り直す」コストで
//     insert が旧版より 1 割遅くなった。そこで Fenwick を遅延再構築 (下記) にしている。
//   ・BK を小さくするとクエリは少し速くなるが、バケット数 nb=N/BK が増えて
//     分割のたびの O(nb) 作業 (配列の中央挿入と Fenwick 再構築) が重くなる。
//     実測では BK=128〜192 がほぼ同等 (差 1% 以内)。既定は 128 にしている。
//   ・N が 4e6 を大きく超えると上段が L2 に乗らなくなるため、この優位は縮む。
//     その規模なら値を座標圧縮して Fenwick に載せる方が速い (静的な問題なら特に)。
// ============================================================================

/**
 * @brief OrderedMultiSet (平方分割 + Fenwick Tree ベースの多重順序集合)
 *
 * 【データ構造の概要】
 * 全要素を昇順に並べたまま、長さ BK 程度の「バケット」に区切って保持する。
 * バケット単位の要素数・総和を Fenwick Tree で管理することで、
 * 「k 番目の要素」「x 未満の個数」「x 未満の総和」などを O(log N) で処理する。
 *
 *   bk[0]              bk[1]              bk[2]        ...     (各バケット内は昇順)
 *   [1,3,3,7]          [8,8,12]           [15,20,20]
 *   bmax = [7,          12,                20]                 ← 代表値の連続配列
 *   Fenwick(個数) = [4, 3, 3]                                   ← 個数の累積
 *   Fenwick(総和) = [14, 28, 55]                                ← 総和の累積
 *
 * 探索は 2 段で終わる:
 *   1. bmax を二分探索して「x が入り得るバケット」を 1 つに絞る (上段・キャッシュに乗る)
 *   2. そのバケット内を二分探索する (下段・ここだけがコールド)
 *
 * 【テンプレート引数】
 *   BK : 1 バケットの目標長 (既定 128)。旧版の B に相当するチューニング用の値だが、
 *        意味は「ノードの最大キー数」から「バケットの目標長」に変わっている。
 *        - 大きいほどバケット数が減り、分割コストと上段の探索が軽くなる。
 *        - 小さいほどバケット内の memmove と二分探索が軽くなる。
 *        - 実測では 128〜192 がほぼ同等 (差 1% 以内)。既定の 128 のままで概ね最速。
 *        - 1 以上なら動作する (小さい値は分割/併合を多発させるストレステスト用)。
 *   T  : 要素の型 (既定 long long)。 operator< で全順序が入っていればよい。
 *   S  : 総和の型 (既定 long long)。 T の総和が long long に収まらない場合は
 *        __int128 などを指定する。例: OrderedMultiSet<128, long long, __int128>
 *   ※ 旧版と同じ書き方 OrderedMultiSet<> / OrderedMultiSet<4> もそのまま使える。
 *
 * 【内部の不変条件】(validate() が全部検査する)
 *   ・空でない限り、どのバケットも空ではなく、長さは 2*BK 以下
 *   ・バケット内は昇順で、bk[i] の全要素 <= bk[i+1] の全要素
 *   ・bmax[i] == bk[i].back()   (bmax は非減少)
 *   ・bsum[i] == bk[i] の総和
 *   ・Fenwick の i 番目の差分が bk[i] の個数 / 総和と一致
 *   ・n_ / sum_ が全体の個数 / 総和と一致
 *   ・長さ 2*BK を超えたら半分に分割し、隣と合わせて BK 以下になったら併合する
 *     (2*BK で分割・BK 以下で併合、というヒステリシスにより
 *      分割と併合が交互に起きて O(nb) を連発することはない)
 *
 * 【計算量】 N = size(), nb = バケット数 ≒ N/BK
 *   insert(x)                : O(BK + log N)         memmove が支配的
 *   erase_one(x)             : O(BK + log N)
 *   erase_at(k)              : O(BK + log N)
 *   erase_all(x)             : O(BK + log N) / 複数バケットに跨るときは O(k/BK + nb)
 *   erase_range(lo,hi)       : 同上 (要素を1つずつ消すより大幅に速い)
 *   erase_index_range(a,b)   : 同上
 *   size / empty / sum       : O(1)
 *   min / max / front / back : O(1)                  ← 旧版は O(log N)
 *   pop_min / pop_max        : O(BK + log N)
 *   rank / count_lt / count_le / count_gt / count_ge / count : O(log N)
 *   count_range / equal_range_index : O(log N)
 *   contains(x)              : O(log N)
 *   sum_lt / sum_le / sum_gt / sum_ge / sum_range : O(BK/2 + log N)
 *   sum_first_k / sum_last_k / sum_index_range    : O(BK/2 + log N)
 *   max_prefix_count(budget) : O(BK + log N)         (非負要素前提)
 *   kth / kth_smallest / kth_largest / at / operator[] : O(log N)
 *   kth_in_range(lo,hi,k)    : O(log N)
 *   next / next_eq / prev / prev_eq : O(log N)
 *   lower_bound / upper_bound / iterator_at (イテレータ) : O(log N)
 *   vec()                    : O(N)
 *   vec_range / vec_index_range / for_each_range : O(log N + 出力数)
 *   build(vector)            : O(N log N) (ソート込み) / ソート済みなら O(N)
 *   insert_all(vector)       : O(N + M) (マージして作り直す)
 *   merge_from(other)        : O(N + M)
 *   split_ge(x)              : O(BK + nb)
 *   clear / コピー            : O(N)
 *   validate()               : O(N + nb log nb)      デバッグ用
 *
 * 【提供関数一覧】
 * - 構築      : OrderedMultiSet(), OrderedMultiSet(vector), build(vector), clear(),
 *               reserve(n), コピー / ムーブ / swap
 * - 挿入      : insert, insert_all, merge_from
 * - 削除      : erase_one, erase_all, erase_at, erase_range, erase_index_range,
 *               pop_min, pop_max, split_ge
 * - 状態      : size, empty, sum, bucket_count
 * - 極値      : min, max, front, back            (front/back は optional で例外を投げない)
 * - 順位      : rank, lower_bound_index, upper_bound_index, count_lt, count_le,
 *               count_gt, count_ge, count, count_range, contains, equal_range_index
 * - 条件付き和: sum_lt, sum_le, sum_gt, sum_ge, sum_range
 * - 順位での和: sum_first_k, sum_last_k, sum_index_range, max_prefix_count
 * - インデックス: kth, kth_smallest, kth_largest, kth_in_range, at, operator[], kth_opt
 * - 前後探索  : next, next_eq, prev, prev_eq
 * - 走査      : vec, vec_range, vec_index_range, for_each_range,
 *               begin, end, lower_bound, upper_bound, iterator_at
 * - 検証      : validate()
 *
 * 【注意点】
 * - min()/max()/kth() は空集合・範囲外で std::runtime_error を投げる。
 *   例外を使いたくない場合は front()/back()/kth_opt() を使う。
 * - sum 系の戻り値は S 型。既定では long long なのでオーバーフローに注意
 *   (要素 1e9 が 1e6 個なら 1e15 で収まるが、1e18 級の値を大量に入れるなら S=__int128)。
 * - max_prefix_count() だけは「全要素が 0 以上」を前提とする (単調性が必要)。
 * - 本構造体は「位置順」ではなく「値順」の集合である。位置指定の挿入・削除が必要なら
 *   DynamicSequence (Implicit Treap) 側を使うこと。
 */
template <int BK = 128, class T = long long, class S = long long>
struct OrderedMultiSet {
    static_assert(BK >= 1, "BK は 1 以上を指定してください");

    // ------------------------------------------------------------------
    // 内部状態
    // ------------------------------------------------------------------
    // bk : 全要素を昇順に並べたものをバケットに区切って保持する。
    //      各バケットは昇順で、bk[i] の全要素 <= bk[i+1] の全要素。
    //      どのバケットも空ではなく、長さは 2*BK 以下。
    vector<vector<T>> bk;
    // bmax[i] == bk[i].back()。連続配列なので二分探索がキャッシュに優しい。
    // 「x が入り得るバケット」はこの配列だけで特定できる。
    vector<T> bmax;
    // bsum[i] == bk[i] の総和。sum_* の「後ろから引く」高速化と分割時に使う。
    vector<S> bsum;
    // バケット単位の Fenwick Tree (1-indexed)。fc は個数、fs は総和。
    //
    // これは bk / bsum から導出できる「キャッシュ」であり、いつでも O(nb) で
    // 作り直せる。バケット数が変わる操作 (分割/併合/削除) のたびに作り直すと
    // 挿入コストの 3 割ほどを食っていたので、fw_dirty を立てておいて
    // 実際に Fenwick が必要になったクエリの直前でまとめて作り直す
    // (遅延再構築)。連続した挿入の途中では 1 度も作り直さずに済む。
    // const なクエリからも作り直せるように mutable にしている。
    mutable vector<long long> fc;
    mutable vector<S> fs;
    mutable int LOG = 0;          // Fenwick 上の二分探索で使う最上位ビット
    mutable bool fw_dirty = true; // true なら fc/fs は未構築
    long long n_ = 0;      // 全体の要素数
    S sum_ = S(0);         // 全体の総和

    int bucket_count() const { return (int)bk.size(); }
    int nb() const { return (int)bk.size(); }

    // ------------------------------------------------------------------
    // ノード内探索の薄いラッパ (メンバ名 lower_bound/upper_bound と衝突するため
    // クラス内では必ず std:: を明示する必要がある。その煩雑さを隠すためのもの)
    // ------------------------------------------------------------------
    static int lower_idx(const vector<T>& v, const T& x) {
        return (int)(std::lower_bound(v.begin(), v.end(), x) - v.begin());
    }
    static int upper_idx(const vector<T>& v, const T& x) {
        return (int)(std::upper_bound(v.begin(), v.end(), x) - v.begin());
    }
    // x が入り得る最初のバケット (bmax[i] >= x となる最小の i)。
    // 返り値が nb() なら「x は全要素より大きい」。
    int find_ge(const T& x) const {
        return (int)(std::lower_bound(bmax.begin(), bmax.end(), x) - bmax.begin());
    }
    // x より大きい要素を含み得る最初のバケット (bmax[i] > x となる最小の i)。
    int find_gt(const T& x) const {
        return (int)(std::upper_bound(bmax.begin(), bmax.end(), x) - bmax.begin());
    }

    // ------------------------------------------------------------------
    // Fenwick Tree (バケット単位)
    // ------------------------------------------------------------------
    // 各バケットの個数・総和から O(nb) で作り直す。
    // バケットの個数が変わる操作 (分割 / 併合 / 削除) の後に呼ぶ。
    void fw_build() const {
        const int m = nb();
        fc.assign(m + 1, 0);
        fs.assign(m + 1, 0);
        for (int i = 0; i < m; i++) { fc[i + 1] = (long long)bk[i].size(); fs[i + 1] = bsum[i]; }
        for (int i = 1; i <= m; i++) {
            const int j = i + (i & -i);
            if (j <= m) { fc[j] += fc[i]; fs[j] += fs[i]; }
        }
        LOG = 0;
        while ((1 << (LOG + 1)) <= m) LOG++;
        fw_dirty = false;
    }
    // Fenwick が必要になった時点で作り直す (遅延再構築)
    void ensure_fw() const { if (fw_dirty) fw_build(); }
    // バケット i の集約だけが dc / ds だけ変わったときの更新。
    // まだ dirty なら触らない (どうせ作り直すので無駄)
    void fw_touch(int i, long long dc, S ds) { if (!fw_dirty) fw_add(i, dc, ds); }
    // バケット i (0-indexed) の個数を dc、総和を ds だけ増やす。 O(log nb)
    void fw_add(int i, long long dc, S ds) {
        for (int p = i + 1; p <= nb(); p += p & -p) { fc[p] += dc; fs[p] += ds; }
    }
    // バケット [0, i) に入っている要素数。 O(log nb)
    long long fw_cnt(int i) const {
        ensure_fw();
        long long r = 0;
        for (int p = i; p > 0; p -= p & -p) r += fc[p];
        return r;
    }
    // バケット [0, i) に入っている要素の総和。 O(log nb)
    S fw_sum(int i) const {
        ensure_fw();
        S r = S(0);
        for (int p = i; p > 0; p -= p & -p) r += fs[p];
        return r;
    }
    // 「小さい方から k 番目 (0-indexed)」が入っているバケット番号を返し、
    // k をそのバケット内のインデックスに書き換える。 O(log nb)
    // 前提: 0 <= k < n_
    int fw_locate(long long& k) const {
        ensure_fw();
        int pos = 0;
        const int m = nb();
        for (int pw = LOG; pw >= 0; pw--) {
            const int np = pos + (1 << pw);
            if (np <= m && fc[np] <= k) { pos = np; k -= fc[np]; }
        }
        return pos;
    }
    // 「累積和が budget を超えない最大のバケット境界」を探し、
    // budget を残額に、cnt をそこまでの個数に書き換えてバケット番号を返す。
    // 前提: 全要素が 0 以上 (累積和が単調非減少であること)。 O(log nb)
    int fw_locate_sum(S& budget, long long& cnt) const {
        ensure_fw();
        int pos = 0;
        const int m = nb();
        cnt = 0;
        for (int pw = LOG; pw >= 0; pw--) {
            const int np = pos + (1 << pw);
            if (np <= m && !(fs[np] > budget)) { pos = np; budget -= fs[np]; cnt += fc[np]; }
        }
        return pos;
    }

    // ------------------------------------------------------------------
    // バケットの分割 / 併合 / 削除
    // ------------------------------------------------------------------
    // バケット i が 2*BK を超えたので半分に割る。 O(BK + nb)
    void split(int i) {
        auto& v = bk[i];
        const int half = (int)v.size() / 2;
        vector<T> right(v.begin() + half, v.end());
        right.reserve(2 * BK + 1);
        v.erase(v.begin() + half, v.end());
        S s0 = S(0), s1 = S(0);
        for (const T& e : v) s0 += S(e);
        for (const T& e : right) s1 += S(e);
        bk.insert(bk.begin() + i + 1, std::move(right));
        bmax.insert(bmax.begin() + i + 1, bk[i + 1].back());
        bsum.insert(bsum.begin() + i + 1, s1);
        bmax[i] = bk[i].back();
        bsum[i] = s0;
        fw_dirty = true;
    }
    // バケット i と i+1 を 1 つにまとめる。 O(BK + nb)
    void merge_next(int i) {
        auto& a = bk[i];
        a.insert(a.end(), bk[i + 1].begin(), bk[i + 1].end());
        bsum[i] += bsum[i + 1];
        bmax[i] = a.back();
        bk.erase(bk.begin() + i + 1);
        bmax.erase(bmax.begin() + i + 1);
        bsum.erase(bsum.begin() + i + 1);
        fw_dirty = true;
    }
    // 空になったバケット i を取り除く。 O(nb)
    void drop(int i) {
        bk.erase(bk.begin() + i);
        bmax.erase(bmax.begin() + i);
        bsum.erase(bsum.begin() + i);
        fw_dirty = true;
    }
    // バケット i が痩せてきたら隣と併合する。
    // 「2*BK で分割・合計 BK 以下で併合」というヒステリシスがあるので、
    // 分割と併合が交互に起きて O(nb) を連発することはない。
    void try_merge(int i) {
        if ((int)bk[i].size() * 2 > BK) return;
        if (i + 1 < nb() && (int)(bk[i].size() + bk[i + 1].size()) <= BK) { merge_next(i); return; }
        if (i > 0 && (int)(bk[i - 1].size() + bk[i].size()) <= BK) { merge_next(i - 1); return; }
    }

    // ソート済み列から全体を作り直す。 O(N)
    void rebuild_from_sorted(const vector<T>& v) {
        bk.clear(); bmax.clear(); bsum.clear();
        n_ = (long long)v.size();
        sum_ = S(0);
        fw_dirty = true;
        if (v.empty()) return;
        // 各バケットが BK 以下になるように均等割りする (分割の余地を残す)
        const size_t g = (v.size() + BK - 1) / BK;
        size_t pos = 0;
        bk.reserve(g); bmax.reserve(g); bsum.reserve(g);
        for (size_t j = 0; j < g; j++) {
            const size_t len = v.size() / g + (j < v.size() % g ? 1 : 0);
            vector<T> cur(v.begin() + pos, v.begin() + pos + len);
            cur.reserve(2 * BK + 1);
            S s = S(0);
            for (const T& e : cur) s += S(e);
            sum_ += s;
            bmax.push_back(cur.back());
            bsum.push_back(s);
            bk.push_back(std::move(cur));
            pos += len;
        }
        fw_dirty = true;
    }

    // ==================================================================
    // 構築 / 破棄 / コピー / ムーブ
    // ==================================================================
    OrderedMultiSet() = default;
    // vector から一括構築する。1 要素ずつ insert するより大幅に速い。 O(N log N)
    explicit OrderedMultiSet(vector<T> v) { build(std::move(v)); }
    // 内部が vector だけなので、コピーとデストラクタは既定実装で安全。
    // (旧版はノードを生ポインタで持っていたため deep copy を手書きする必要があった)
    OrderedMultiSet(const OrderedMultiSet&) = default;
    OrderedMultiSet& operator=(const OrderedMultiSet&) = default;
    // ムーブだけは既定実装に任せられない。
    // vector はムーブされて空になるが、n_ / sum_ はスカラなのでコピーされてしまい、
    // 「バケットは空なのに size() が 0 でない」という不整合が残る。
    // そのため明示的に実装し、ムーブ元を clear() で確実に空にしている。
    OrderedMultiSet(OrderedMultiSet&& o) noexcept
        : bk(std::move(o.bk)), bmax(std::move(o.bmax)), bsum(std::move(o.bsum)),
          fc(std::move(o.fc)), fs(std::move(o.fs)), LOG(o.LOG), fw_dirty(o.fw_dirty),
          n_(o.n_), sum_(o.sum_) {
        o.clear();
    }
    OrderedMultiSet& operator=(OrderedMultiSet&& o) noexcept {
        if (this != &o) {
            bk = std::move(o.bk); bmax = std::move(o.bmax); bsum = std::move(o.bsum);
            fc = std::move(o.fc); fs = std::move(o.fs);
            n_ = o.n_; sum_ = o.sum_; LOG = o.LOG; fw_dirty = o.fw_dirty;
            o.clear();
        }
        return *this;
    }

    void swap(OrderedMultiSet& o) noexcept {
        bk.swap(o.bk); bmax.swap(o.bmax); bsum.swap(o.bsum);
        fc.swap(o.fc); fs.swap(o.fs);
        std::swap(n_, o.n_); std::swap(sum_, o.sum_); std::swap(LOG, o.LOG);
        std::swap(fw_dirty, o.fw_dirty);
    }

    // vector から一括構築する (既存の要素はすべて破棄される)。
    // O(N log N) (ソート込み。ソート済みを渡せば O(N))
    void build(vector<T> v) {
        if (!std::is_sorted(v.begin(), v.end())) std::sort(v.begin(), v.end());
        rebuild_from_sorted(v);
    }
    // 全要素を削除する。 O(N)
    void clear() {
        bk.clear(); bmax.clear(); bsum.clear(); fc.clear(); fs.clear();
        n_ = 0; sum_ = S(0); LOG = 0; fw_dirty = true;
    }
    // 要素数の目安が分かっているならバケット配列を事前確保できる。
    void reserve(long long expected_n) {
        const size_t g = (size_t)(expected_n / BK + 2);
        bk.reserve(g); bmax.reserve(g); bsum.reserve(g);
    }

    // ==================================================================
    // 挿入
    // ==================================================================
    // 値 x を挿入する (重複可)。 O(BK + log N)
    void insert(const T& x) {
        if (bk.empty()) {
            bk.emplace_back();
            bk[0].reserve(2 * BK + 1);
            bk[0].push_back(x);
            bmax.push_back(x);
            bsum.push_back(S(x));
            n_ = 1; sum_ = S(x);
            fw_dirty = true;
            return;
        }
        int i = find_ge(x);
        if (i == nb()) i = nb() - 1;          // 全要素より大きいなら末尾のバケットへ
        auto& v = bk[i];
        v.insert(v.begin() + upper_idx(v, x), x);
        bsum[i] += S(x);
        bmax[i] = v.back();
        fw_touch(i, 1, S(x));
        ++n_; sum_ += S(x);
        if ((int)v.size() > 2 * BK) split(i);
    }

    // 列をまとめて挿入する。1 要素ずつ insert するより速い。 O(N + M)
    void insert_all(vector<T> v) {
        if (v.empty()) return;
        if (!std::is_sorted(v.begin(), v.end())) std::sort(v.begin(), v.end());
        if (bk.empty()) { rebuild_from_sorted(v); return; }
        vector<T> cur = vec();
        vector<T> merged(cur.size() + v.size());
        std::merge(cur.begin(), cur.end(), v.begin(), v.end(), merged.begin());
        rebuild_from_sorted(merged);
    }
    // other の全要素を自分に取り込む (other は空になる)。 O(N + M)
    void merge_from(OrderedMultiSet& other) {
        if (other.empty()) return;
        if (empty()) { swap(other); return; }
        vector<T> a = vec(), b = other.vec();
        vector<T> merged(a.size() + b.size());
        std::merge(a.begin(), a.end(), b.begin(), b.end(), merged.begin());
        other.clear();
        rebuild_from_sorted(merged);
    }

    // ==================================================================
    // 削除
    // ==================================================================
    // 値 x を 1 つだけ削除する。成功したら true。 O(BK + log N)
    bool erase_one(const T& x) {
        if (bk.empty()) return false;
        const int i = find_ge(x);
        if (i == nb()) return false;
        auto& v = bk[i];
        const int p = lower_idx(v, x);
        if (p == (int)v.size() || x < v[p]) return false;
        v.erase(v.begin() + p);
        bsum[i] -= S(x);
        fw_touch(i, -1, -S(x));
        --n_; sum_ -= S(x);
        if (v.empty()) drop(i);
        else { bmax[i] = v.back(); try_merge(i); }
        return true;
    }

    // ソート順で半開区間 [a, b) 番目 (0-indexed) の要素をまとめて削除し、個数を返す。
    // バケットを丸ごと落とせるので、1 要素ずつ削除するより大幅に速い。
    // 1 バケットに収まる場合は O(BK + log N)、跨る場合は O((b-a)/BK + BK + nb)。
    long long erase_index_range(long long a, long long b) {
        if (a < 0) a = 0;
        if (b > n_) b = n_;
        if (a >= b) return 0;
        const long long removed = b - a;

        long long ka = a;
        const int ia = fw_locate(ka);          // ia バケットの ka 番目が先頭
        long long kb = b - 1;
        const int ib = fw_locate(kb);          // ib バケットの kb 番目が末尾 (含む)

        if (ia == ib) {
            // 1 バケット内で収まるケース。Fenwick は差分更新で済む。
            auto& v = bk[ia];
            S ds = S(0);
            for (long long j = ka; j <= kb; j++) ds += S(v[(size_t)j]);
            v.erase(v.begin() + ka, v.begin() + kb + 1);
            bsum[ia] -= ds;
            fw_touch(ia, -removed, -ds);
            n_ -= removed; sum_ -= ds;
            if (v.empty()) drop(ia);
            else { bmax[ia] = v.back(); try_merge(ia); }
            return removed;
        }

        // 複数バケットに跨るケース。
        // ia の末尾側を切り、ib の先頭側を切り、その間のバケットは丸ごと捨てる。
        S ds = S(0);
        {
            auto& v = bk[ia];
            for (size_t j = (size_t)ka; j < v.size(); j++) ds += S(v[j]);
            v.erase(v.begin() + ka, v.end());
        }
        for (int j = ia + 1; j < ib; j++) ds += bsum[j];
        {
            auto& v = bk[ib];
            for (long long j = 0; j <= kb; j++) ds += S(v[(size_t)j]);
            v.erase(v.begin(), v.begin() + kb + 1);
        }
        // 空になったバケットを取り除きつつ、間のバケットを消す
        bk.erase(bk.begin() + ia + 1, bk.begin() + ib);
        bmax.erase(bmax.begin() + ia + 1, bmax.begin() + ib);
        bsum.erase(bsum.begin() + ia + 1, bsum.begin() + ib);
        int left = ia, right = ia + 1;         // 残った両端のバケット位置
        if (bk[right].empty()) { bk.erase(bk.begin() + right); bmax.erase(bmax.begin() + right); bsum.erase(bsum.begin() + right); }
        else { bsum[right] = S(0); for (const T& e : bk[right]) bsum[right] += S(e); bmax[right] = bk[right].back(); }
        if (bk[left].empty()) { bk.erase(bk.begin() + left); bmax.erase(bmax.begin() + left); bsum.erase(bsum.begin() + left); left--; }
        else { bsum[left] = S(0); for (const T& e : bk[left]) bsum[left] += S(e); bmax[left] = bk[left].back(); }
        n_ -= removed; sum_ -= ds;
        fw_dirty = true;
        if (left >= 0 && left < nb()) try_merge(left);
        else if (nb() > 0) try_merge(std::min(std::max(left, 0), nb() - 1));
        return removed;
    }

    // 値 x をすべて削除し、削除した個数を返す。 O(BK + log N) 〜 O(k/BK + nb)
    long long erase_all(const T& x) {
        const long long lo = count_lt(x);
        const long long hi = count_le(x);
        return erase_index_range(lo, hi);
    }
    // 値が半開区間 [lo, hi) に入る要素をすべて削除し、個数を返す。
    long long erase_range(const T& lo, const T& hi) {
        if (!(lo < hi)) return 0;
        return erase_index_range(count_lt(lo), count_lt(hi));
    }
    // 小さい方から k 番目 (0-indexed) の要素を削除する。
    // 範囲外なら false を返して何もしない。 O(BK + log N)
    bool erase_at(long long k) {
        if (k < 0 || k >= n_) return false;
        erase_index_range(k, k + 1);
        return true;
    }
    // 最小値 / 最大値を取り出して削除する。空なら nullopt。 O(BK + log N)
    optional<T> pop_min() {
        if (bk.empty()) return nullopt;
        const T v = bk.front().front();
        erase_index_range(0, 1);
        return v;
    }
    optional<T> pop_max() {
        if (bk.empty()) return nullopt;
        const T v = bk.back().back();
        erase_index_range(n_ - 1, n_);
        return v;
    }

    // 値が x 以上の要素をすべて切り離し、新しい OrderedMultiSet として返す。
    // 自分側には x 未満の要素だけが残る。 O(BK + nb)
    // (バケット列をそのまま付け替えるだけなので要素数に比例しない)
    OrderedMultiSet split_ge(const T& x) {
        OrderedMultiSet res;
        if (bk.empty()) return res;
        const int i = find_ge(x);
        if (i == nb()) return res;             // 全要素が x 未満
        // バケット i を x の位置で 2 つに分ける
        auto& v = bk[i];
        const int p = lower_idx(v, x);
        vector<T> tail(v.begin() + p, v.end());
        v.erase(v.begin() + p, v.end());
        // i+1 以降のバケットを res へ移す
        for (int j = i + 1; j < nb(); j++) {
            res.bk.push_back(std::move(bk[j]));
            res.bmax.push_back(bmax[j]);
            res.bsum.push_back(bsum[j]);
        }
        bk.resize(i + 1); bmax.resize(i + 1); bsum.resize(i + 1);
        if (!tail.empty()) {
            S st = S(0);
            for (const T& e : tail) st += S(e);
            const T tmax = tail.back();
            tail.reserve(2 * BK + 1);
            res.bk.insert(res.bk.begin(), std::move(tail));
            res.bmax.insert(res.bmax.begin(), tmax);
            res.bsum.insert(res.bsum.begin(), st);
        }
        if (bk[i].empty()) { bk.pop_back(); bmax.pop_back(); bsum.pop_back(); }
        else {
            bsum[i] = S(0);
            for (const T& e : bk[i]) bsum[i] += S(e);
            bmax[i] = bk[i].back();
        }
        // 個数と総和を数え直す (バケット数に比例、要素数には比例しない)
        n_ = 0; sum_ = S(0);
        for (int j = 0; j < nb(); j++) { n_ += (long long)bk[j].size(); sum_ += bsum[j]; }
        res.n_ = 0; res.sum_ = S(0);
        for (int j = 0; j < res.nb(); j++) { res.n_ += (long long)res.bk[j].size(); res.sum_ += res.bsum[j]; }
        fw_dirty = true;
        res.fw_dirty = true;
        return res;
    }

    // ==================================================================
    // 状態取得
    // ==================================================================
    long long size() const { return n_; }          // O(1)
    bool empty() const { return n_ == 0; }         // O(1)
    S sum() const { return sum_; }                 // O(1)

    // ==================================================================
    // 極値 (バケットの両端そのものなので O(1))
    // ==================================================================
    T min() const {
        if (bk.empty()) throw runtime_error("OrderedMultiSet::min() : Empty set");
        return bk.front().front();
    }
    T max() const {
        if (bk.empty()) throw runtime_error("OrderedMultiSet::max() : Empty set");
        return bk.back().back();
    }
    optional<T> front() const { if (bk.empty()) return nullopt; return bk.front().front(); }
    optional<T> back() const { if (bk.empty()) return nullopt; return bk.back().back(); }

    // ==================================================================
    // 順位・カウント
    // ==================================================================
    // 値 x 未満の要素数 (= x を挿入するときの 0-indexed 位置)。 O(log N)
    long long count_lt(const T& x) const {
        if (bk.empty()) return 0;
        const int i = find_ge(x);
        if (i == nb()) return n_;
        return fw_cnt(i) + (long long)lower_idx(bk[i], x);
    }
    // 値 x 以下の要素数。 O(log N)
    long long count_le(const T& x) const {
        if (bk.empty()) return 0;
        const int i = find_gt(x);
        if (i == nb()) return n_;
        return fw_cnt(i) + (long long)upper_idx(bk[i], x);
    }
    long long rank(const T& x) const { return count_lt(x); }
    long long lower_bound_index(const T& x) const { return count_lt(x); }
    long long upper_bound_index(const T& x) const { return count_le(x); }
    long long count_gt(const T& x) const { return n_ - count_le(x); }
    long long count_ge(const T& x) const { return n_ - count_lt(x); }
    long long count(const T& x) const { return count_le(x) - count_lt(x); }
    // 値が半開区間 [lo, hi) に入る要素数。 lo >= hi なら 0。 O(log N)
    long long count_range(const T& lo, const T& hi) const {
        if (!(lo < hi)) return 0;
        return count_lt(hi) - count_lt(lo);
    }
    // x と等しい要素が占めるインデックス範囲 [first, last)。 O(log N)
    pair<long long, long long> equal_range_index(const T& x) const {
        return {count_lt(x), count_le(x)};
    }
    // 値 x が含まれるか。 O(log N)
    bool contains(const T& x) const {
        if (bk.empty()) return false;
        const int i = find_ge(x);
        if (i == nb()) return false;
        const auto& v = bk[i];
        const int p = lower_idx(v, x);
        return p != (int)v.size() && !(x < v[p]);
    }

    // ==================================================================
    // 条件付き総和
    // ==================================================================
    // バケット i の先頭 p 個の総和。p が大きいときは「総和から後ろを引く」。 O(BK/2)
    S partial_sum(int i, int p) const {
        const auto& v = bk[i];
        const int m = (int)v.size();
        if (p * 2 <= m) {
            S r = S(0);
            for (int j = 0; j < p; j++) r += S(v[j]);
            return r;
        }
        S r = bsum[i];
        for (int j = p; j < m; j++) r -= S(v[j]);
        return r;
    }
    // 値 x 未満の要素の総和。 O(BK/2 + log N)
    S sum_lt(const T& x) const {
        if (bk.empty()) return S(0);
        const int i = find_ge(x);
        if (i == nb()) return sum_;
        return fw_sum(i) + partial_sum(i, lower_idx(bk[i], x));
    }
    // 値 x 以下の要素の総和。 O(BK/2 + log N)
    S sum_le(const T& x) const {
        if (bk.empty()) return S(0);
        const int i = find_gt(x);
        if (i == nb()) return sum_;
        return fw_sum(i) + partial_sum(i, upper_idx(bk[i], x));
    }
    S sum_gt(const T& x) const { return sum_ - sum_le(x); }
    S sum_ge(const T& x) const { return sum_ - sum_lt(x); }
    // 値が半開区間 [l, r) に入る要素の総和。 l >= r なら 0。 O(BK + log N)
    S sum_range(const T& l, const T& r) const {
        if (!(l < r)) return S(0);
        return sum_lt(r) - sum_lt(l);
    }

    // ==================================================================
    // 順位 (インデックス) による総和
    // ==================================================================
    // 小さい方から k 個の総和。k<=0 なら 0、k>=size() なら全体。 O(BK/2 + log N)
    // 例: 「一番安い k 個を買うときの合計金額」がそのまま求まる。
    S sum_first_k(long long k) const {
        if (bk.empty() || k <= 0) return S(0);
        if (k >= n_) return sum_;
        long long r = k;
        const int i = fw_locate(r);            // r はバケット内インデックス
        return fw_sum(i) + partial_sum(i, (int)r);
    }
    // 大きい方から k 個の総和。 O(BK/2 + log N)
    S sum_last_k(long long k) const {
        if (k <= 0) return S(0);
        if (k >= n_) return sum_;
        return sum_ - sum_first_k(n_ - k);
    }
    // ソート順で半開区間 [a, b) 番目の要素の総和。 O(BK + log N)
    S sum_index_range(long long a, long long b) const {
        if (a < 0) a = 0;
        if (b > n_) b = n_;
        if (a >= b) return S(0);
        return sum_first_k(b) - sum_first_k(a);
    }
    // 小さい方から貪欲に取ったとき、総和が budget 以下に収まる最大個数。
    // 【前提】すべての要素が 0 以上であること (負の値があると単調性が崩れる)。
    // O(BK + log N)
    long long max_prefix_count(S budget) const {
        if (bk.empty() || budget < S(0)) return 0;
        if (!(sum_ > budget)) return n_;
        S rem = budget;
        long long cnt = 0;
        const int i = fw_locate_sum(rem, cnt);
        if (i >= nb()) return cnt;
        for (const T& e : bk[i]) {
            if (S(e) > rem) break;
            rem -= S(e);
            ++cnt;
        }
        return cnt;
    }

    // ==================================================================
    // インデックスアクセス
    // ==================================================================
    // 小さい方から k 番目 (0-indexed)。範囲外なら runtime_error。 O(log N)
    T kth(long long k) const {
        if (k < 0 || k >= n_) throw runtime_error("OrderedMultiSet::kth() : index out of bounds");
        long long r = k;
        const int i = fw_locate(r);
        return bk[i][(size_t)r];
    }
    T kth_smallest(long long k) const { return kth(k); }
    T kth_largest(long long k) const { return kth(n_ - 1 - k); }
    T at(long long k) const { return kth(k); }
    T operator[](long long k) const { return kth(k); }
    // 例外を投げない kth。範囲外なら nullopt。 O(log N)
    optional<T> kth_opt(long long k) const {
        if (k < 0 || k >= n_) return nullopt;
        return kth(k);
    }
    // 値が半開区間 [lo, hi) に入る要素のうち、小さい方から k 番目 (0-indexed)。
    // 足りなければ nullopt。 O(log N)
    optional<T> kth_in_range(const T& lo, const T& hi, long long k) const {
        if (!(lo < hi) || k < 0) return nullopt;
        const long long a = count_lt(lo), b = count_lt(hi);
        if (k >= b - a) return nullopt;
        return kth(a + k);
    }

    // ==================================================================
    // 前後要素探索 (いずれも 1 バケットを見るだけで済む)
    // ==================================================================
    // x より大きい (strict >) 最小の要素。 O(log N)
    optional<T> next(const T& x) const {
        if (bk.empty()) return nullopt;
        const int i = find_gt(x);              // bmax[i] > x なので必ずこのバケットにある
        if (i == nb()) return nullopt;
        return bk[i][(size_t)upper_idx(bk[i], x)];
    }
    // x 以上 (>=) の最小の要素。 O(log N)
    optional<T> next_eq(const T& x) const {
        if (bk.empty()) return nullopt;
        const int i = find_ge(x);              // bmax[i] >= x なので必ずこのバケットにある
        if (i == nb()) return nullopt;
        return bk[i][(size_t)lower_idx(bk[i], x)];
    }
    // x より小さい (strict <) 最大の要素。 O(log N)
    optional<T> prev(const T& x) const {
        if (bk.empty()) return nullopt;
        const int i = find_ge(x);
        if (i == nb()) return bk.back().back();          // x は全要素より大きい
        const int p = lower_idx(bk[i], x);
        if (p > 0) return bk[i][(size_t)(p - 1)];
        if (i > 0) return bk[i - 1].back();
        return nullopt;
    }
    // x 以下 (<=) の最大の要素。 O(log N)
    optional<T> prev_eq(const T& x) const {
        if (bk.empty()) return nullopt;
        const int i = find_gt(x);
        if (i == nb()) return bk.back().back();
        const int p = upper_idx(bk[i], x);
        if (p > 0) return bk[i][(size_t)(p - 1)];
        if (i > 0) return bk[i - 1].back();
        return nullopt;
    }

    // ==================================================================
    // 走査
    // ==================================================================
    // 全要素を昇順に並べた vector を返す。 O(N)
    vector<T> vec() const {
        vector<T> res;
        res.reserve((size_t)n_);
        for (const auto& v : bk) res.insert(res.end(), v.begin(), v.end());
        return res;
    }
    // 値が半開区間 [lo, hi) に入る要素を昇順に列挙する。 O(log N + 出力数)
    vector<T> vec_range(const T& lo, const T& hi) const {
        vector<T> res;
        if (bk.empty() || !(lo < hi)) return res;
        int i = find_ge(lo);
        int p = (i == nb()) ? 0 : lower_idx(bk[i], lo);
        for (; i < nb(); i++, p = 0) {
            const auto& v = bk[i];
            for (int j = p; j < (int)v.size(); j++) {
                if (!(v[j] < hi)) return res;
                res.push_back(v[j]);
            }
        }
        return res;
    }
    // ソート順で半開区間 [a, b) 番目の要素を昇順に列挙する。 O(log N + (b-a))
    vector<T> vec_index_range(long long a, long long b) const {
        vector<T> res;
        if (a < 0) a = 0;
        if (b > n_) b = n_;
        if (a >= b) return res;
        res.reserve((size_t)(b - a));
        long long r = a;
        int i = fw_locate(r);
        long long need = b - a;
        for (int p = (int)r; i < nb(); i++, p = 0) {
            const auto& v = bk[i];
            for (int j = p; j < (int)v.size() && need > 0; j++, --need) res.push_back(v[j]);
            if (need == 0) break;
        }
        return res;
    }
    // 値が [lo, hi) の要素に対して f を順に適用する (vector を作らないので速い)。
    // O(log N + 該当数)
    template <class F>
    void for_each_range(const T& lo, const T& hi, F f) const {
        if (bk.empty() || !(lo < hi)) return;
        int i = find_ge(lo);
        int p = (i == nb()) ? 0 : lower_idx(bk[i], lo);
        for (; i < nb(); i++, p = 0) {
            const auto& v = bk[i];
            for (int j = p; j < (int)v.size(); j++) {
                if (!(v[j] < hi)) return;
                f(v[j]);
            }
        }
    }

    // 昇順の前方向イテレータ (範囲 for に使える)。
    // 注意: 集合を変更するとイテレータは無効になる。
    struct const_iterator {
        const OrderedMultiSet* owner = nullptr;
        int bi = 0;           // バケット番号 (owner->nb() なら end)
        int pos = 0;          // バケット内の位置
        using iterator_category = std::forward_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = const T*;
        using reference = const T&;
        const T& operator*() const { return owner->bk[bi][pos]; }
        const T* operator->() const { return &owner->bk[bi][pos]; }
        const_iterator& operator++() {
            if (++pos >= (int)owner->bk[bi].size()) { ++bi; pos = 0; }
            return *this;
        }
        const_iterator operator++(int) { const_iterator t = *this; ++(*this); return t; }
        bool operator==(const const_iterator& o) const { return bi == o.bi && pos == o.pos; }
        bool operator!=(const const_iterator& o) const { return !(*this == o); }
    };
    const_iterator begin() const { return const_iterator{this, 0, 0}; }
    const_iterator end() const { return const_iterator{this, nb(), 0}; }
    // x 以上の最初の要素を指すイテレータ (std::multiset::lower_bound 相当)。 O(log N)
    const_iterator lower_bound(const T& x) const {
        if (bk.empty()) return end();
        const int i = find_ge(x);
        if (i == nb()) return end();
        return const_iterator{this, i, lower_idx(bk[i], x)};
    }
    // x より大きい最初の要素を指すイテレータ。 O(log N)
    const_iterator upper_bound(const T& x) const {
        if (bk.empty()) return end();
        const int i = find_gt(x);
        if (i == nb()) return end();
        return const_iterator{this, i, upper_idx(bk[i], x)};
    }
    // k 番目 (0-indexed) の要素を指すイテレータ。k==size() なら end()。 O(log N)
    const_iterator iterator_at(long long k) const {
        if (k < 0) k = 0;
        if (k >= n_) return end();
        long long r = k;
        const int i = fw_locate(r);
        return const_iterator{this, i, (int)r};
    }

    // ==================================================================
    // 検証 (デバッグ用)
    // ==================================================================
    // 不変条件をすべて検査する。壊れていたら assert で落ちる。 O(N + nb log nb)
    bool validate() const {
        ensure_fw();          // 遅延再構築を確定させてから検査する
        assert((bk.empty() == (n_ == 0)) && "バケット数と要素数の整合性が壊れています");
        assert(bmax.size() == bk.size() && bsum.size() == bk.size() && "補助配列の長さが不正です");
        // 空のときは Fenwick を未構築 (長さ0) でも良いことにしている
        // (既定構築直後と clear() 直後がこれに該当する。nb()==0 なので誰も参照しない)
        assert(fc.size() == bk.size() + 1 && fs.size() == bk.size() + 1 && "Fenwick の長さが不正です");
        long long total = 0;
        S tsum = S(0);
        for (int i = 0; i < nb(); i++) {
            const auto& v = bk[i];
            assert(!v.empty() && "空のバケットが残っています");
            assert((int)v.size() <= 2 * BK && "バケットが 2*BK を超えています");
            assert(std::is_sorted(v.begin(), v.end()) && "バケット内が昇順ではありません");
            assert(!(bmax[i] < v.back()) && !(v.back() < bmax[i]) && "bmax がバケット末尾と一致しません");
            S s = S(0);
            for (const T& e : v) s += S(e);
            assert(bsum[i] == s && "bsum がバケットの総和と一致しません");
            if (i > 0) assert(!(v.front() < bk[i - 1].back()) && "バケット間の順序が壊れています");
            // Fenwick の差分がバケットの個数 / 総和と一致すること
            assert(fw_cnt(i + 1) - fw_cnt(i) == (long long)v.size() && "Fenwick(個数) が一致しません");
            assert(fw_sum(i + 1) - fw_sum(i) == s && "Fenwick(総和) が一致しません");
            total += (long long)v.size();
            tsum += s;
        }
        assert(total == n_ && "要素数が一致しません");
        assert(tsum == sum_ && "総和が一致しません");
        int lg = 0;
        while ((1 << (lg + 1)) <= nb()) lg++;
        assert(lg == LOG && "LOG が不正です");
        return true;
    }
};

// ============================================================================
// ここから下はテスト / 動作確認用のコードです。
// AtCoder に提出するときは、この領域と solv() 内の呼び出しを削除してください。
// ============================================================================

// ---------------------------------------------------------------------------
// [検証ヘルパ] multiset を正解として、全クエリを突き合わせる
// ---------------------------------------------------------------------------
template <int B>
void verify_all_against_multiset(OrderedMultiSet<B>& bst, const multiset<long long>& ms,
                                 mt19937_64& rng, long long val_range) {
    const vector<long long> v(ms.begin(), ms.end());

    // --- 基本状態 ---
    assert(bst.size() == (long long)ms.size());
    assert(bst.empty() == ms.empty());
    {
        long long s = 0;
        for (long long x : v) s += x;
        assert(bst.sum() == s);
    }

    // --- 全要素列挙 ---
    assert(bst.vec() == v);
    {   // range-for (イテレータ) でも同じ列が得られること
        vector<long long> it_v;
        for (long long x : bst) it_v.push_back(x);
        assert(it_v == v);
    }

    if (ms.empty()) {
        assert(!bst.front().has_value() && !bst.back().has_value());
        return;
    }

    // --- 極値 ---
    assert(bst.min() == v.front());
    assert(bst.max() == v.back());
    assert(bst.front().value() == v.front());
    assert(bst.back().value() == v.back());

    // ある値 x について、値に関する全クエリを愚直計算と突き合わせる
    auto check_one_value = [&](long long x) {
        const long long lt = (long long)(std::lower_bound(v.begin(), v.end(), x) - v.begin());
        const long long le = (long long)(std::upper_bound(v.begin(), v.end(), x) - v.begin());
        assert(bst.rank(x) == lt);
        assert(bst.count_lt(x) == lt);
        assert(bst.count_le(x) == le);
        assert(bst.count_gt(x) == (long long)v.size() - le);
        assert(bst.count_ge(x) == (long long)v.size() - lt);
        assert(bst.count(x) == le - lt);
        assert(bst.contains(x) == (le - lt > 0));
        assert(bst.equal_range_index(x) == make_pair(lt, le));

        // 条件付き総和
        long long s_lt = 0, s_le = 0, s_all = 0;
        for (long long e : v) {
            s_all += e;
            if (e < x) s_lt += e;
            if (e <= x) s_le += e;
        }
        assert(bst.sum_lt(x) == s_lt);
        assert(bst.sum_le(x) == s_le);
        assert(bst.sum_gt(x) == s_all - s_le);
        assert(bst.sum_ge(x) == s_all - s_lt);

        // 前後要素探索
        auto it_next = ms.upper_bound(x);
        if (it_next == ms.end()) assert(!bst.next(x).has_value());
        else assert(bst.next(x).value() == *it_next);

        auto it_next_eq = ms.lower_bound(x);
        if (it_next_eq == ms.end()) assert(!bst.next_eq(x).has_value());
        else assert(bst.next_eq(x).value() == *it_next_eq);

        if (it_next_eq == ms.begin()) assert(!bst.prev(x).has_value());
        else assert(bst.prev(x).value() == *std::prev(it_next_eq));

        if (it_next == ms.begin()) assert(!bst.prev_eq(x).has_value());
        else assert(bst.prev_eq(x).value() == *std::prev(it_next));

        // イテレータ版 lower_bound / upper_bound
        auto lb = bst.lower_bound(x);
        if (it_next_eq == ms.end()) assert(lb == bst.end());
        else assert(*lb == *it_next_eq);
        auto ub = bst.upper_bound(x);
        if (it_next == ms.end()) assert(ub == bst.end());
        else assert(*ub == *it_next);
    };

    // 値域からランダムに取った x と、実際に存在する要素の x の両方で検証する
    // (存在する値を狙わないと count/contains の「ヒット」経路が検査されにくい)
    check_one_value((long long)(rng() % (unsigned long long)(val_range + 4)) - 2);
    check_one_value(v[rng() % v.size()]);

    {   // 値の区間 [l, r)
        long long l = v[rng() % v.size()];
        long long r = v[rng() % v.size()];
        if (l > r) std::swap(l, r);
        long long s = 0, c = 0;
        for (long long e : v) if (l <= e && e < r) { s += e; c++; }
        assert(bst.sum_range(l, r) == s);
        assert(bst.count_range(l, r) == c);
        assert(bst.vec_range(l, r) == vector<long long>(v.begin() + (std::lower_bound(v.begin(), v.end(), l) - v.begin()),
                                                        v.begin() + (std::lower_bound(v.begin(), v.end(), r) - v.begin())));
        // 左端 >= 右端 (逆転区間・空区間) は 0 / 空を返す
        assert(bst.sum_range(r, l) == 0);
        assert(bst.count_range(r, l) == 0);
        assert(bst.vec_range(r, l).empty());
    }

    {   // インデックスアクセス
        const long long k = (long long)(rng() % v.size());
        assert(bst.kth(k) == v[k]);
        assert(bst.kth_smallest(k) == v[k]);
        assert(bst.at(k) == v[k]);
        assert(bst[k] == v[k]);
        assert(bst.kth_largest(k) == v[v.size() - 1 - k]);
        assert(bst.kth_opt(k).value() == v[k]);
        assert(!bst.kth_opt(-1).has_value());
        assert(!bst.kth_opt((long long)v.size()).has_value());
    }

    {   // 順位 (インデックス) による総和
        long long a = (long long)(rng() % (v.size() + 1));
        long long b = (long long)(rng() % (v.size() + 1));
        if (a > b) std::swap(a, b);
        auto prefix = [&](long long k) { long long s = 0; for (long long i = 0; i < k; i++) s += v[i]; return s; };
        const long long pa = prefix(a), pb = prefix(b);
        assert(bst.sum_first_k(a) == pa);
        assert(bst.sum_first_k(-5) == 0);
        assert(bst.sum_first_k((long long)v.size() + 5) == bst.sum());
        long long tail = 0;
        for (long long i = (long long)v.size() - a; i < (long long)v.size(); i++) tail += v[i];
        assert(bst.sum_last_k(a) == tail);
        assert(bst.sum_index_range(a, b) == pb - pa);
        assert(bst.sum_index_range(b, a) == 0);                       // 逆転区間
        assert(bst.sum_index_range(-10, (long long)v.size() + 10) == bst.sum());  // クランプ
    }
}

// ---------------------------------------------------------------------------
// [テスト1] 標準ランダムテスト (multiset と全面突き合わせ)
// ---------------------------------------------------------------------------
template <int B>
void run_tests(int num_queries, int val_range, const string& test_name) {
    OrderedMultiSet<B> bst;
    multiset<long long> ms;
    long long current_sum = 0;
    mt19937_64 rng(1337);

    cout << "Starting " << test_name << " (BK=" << B << ", " << num_queries << " queries)..." << endl;

    for (int i = 0; i < num_queries; ++i) {
        int op = rng() % 100;
        long long val = rng() % val_range;

        if (op < 45) { // insert
            bst.insert(val);
            ms.insert(val);
            current_sum += val;
        }
        else if (op < 60) { // erase_one
            bool b_erased = bst.erase_one(val);
            auto it = ms.find(val);
            bool ms_erased = (it != ms.end());
            assert(b_erased == ms_erased);
            if (ms_erased) {
                ms.erase(it);
                current_sum -= val;
            }
        }
        else if (op < 70) { // erase_all
            long long b_erased = bst.erase_all(val);
            long long ms_erased = ms.count(val);
            assert(b_erased == ms_erased);
            if (ms_erased > 0) {
                ms.erase(val);
                current_sum -= val * ms_erased;
            }
        }
        else if (op < 74 && !ms.empty()) { // erase_at (新規API)
            long long k = rng() % ms.size();
            auto it = ms.begin(); advance(it, k);
            assert(bst.erase_at(k));
            current_sum -= *it;
            ms.erase(it);
        }
        else if (op < 76 && !ms.empty()) { // pop_min / pop_max (新規API)
            if (rng() % 2) {
                auto v = bst.pop_min();
                assert(v.has_value() && v.value() == *ms.begin());
                current_sum -= v.value();
                ms.erase(ms.begin());
            } else {
                auto v = bst.pop_max();
                assert(v.has_value() && v.value() == *ms.rbegin());
                current_sum -= v.value();
                ms.erase(std::prev(ms.end()));
            }
        }

        assert(bst.size() == (long long)ms.size());
        assert(bst.empty() == ms.empty());
        assert(bst.sum() == current_sum);

        // 重いクエリ検証は 10 回に 1 回程度
        if (rng() % 10 == 0) verify_all_against_multiset(bst, ms, rng, val_range);
        // 不変条件の検査は 200 回に 1 回 (O(N) なので)
        if (i % 200 == 0) assert(bst.validate());
    }

    assert(bst.validate());
    vector<long long> bst_vec = bst.vec();
    vector<long long> ms_vec(ms.begin(), ms.end());
    assert(bst_vec == ms_vec);

    cout << test_name << " passed successfully!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト2] 小さい BK でバケットの分割 / 併合 / 削除を多発させる全探索的ストレス
// ---------------------------------------------------------------------------
void stress_exhaustive_test() {
    cout << "Running Exhaustive Small Stress Test..." << endl;
    for (int seed = 0; seed < 300; ++seed) {
        OrderedMultiSet<4> bst;
        multiset<long long> ms;
        mt19937 rng(seed);

        for (int i = 0; i < 500; ++i) {
            int val = rng() % 10;
            if (rng() % 2) {
                bst.insert(val);
                ms.insert(val);
            } else if (!ms.empty()) {
                // 削除対象をランダム選択
                auto it = ms.begin();
                advance(it, rng() % ms.size());
                bst.erase_one(*it);
                ms.erase(it);
            }

            // 検証
            assert(bst.size() == (long long)ms.size());
            assert(bst.validate());
            if (!ms.empty()) {
                assert(bst.min() == *ms.begin());
                assert(bst.max() == *ms.rbegin());

                // 全要素走査の一致
                vector<long long> bst_v = bst.vec();
                vector<long long> ms_v(ms.begin(), ms.end());
                assert(bst_v == ms_v);
            }
        }
    }
    cout << "Exhaustive stress test passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト3] 各種の極端な入力パターン
// ---------------------------------------------------------------------------
void same_value_test() {
    cout << "Running Same Value Test..." << endl;

    OrderedMultiSet<64> bst;
    const int N = 100000;
    for (int i = 0; i < N; i++) bst.insert(7);

    assert(bst.size() == N);
    assert(bst.min() == 7);
    assert(bst.max() == 7);
    assert(bst.count(7) == N);
    assert(bst.count_lt(7) == 0);
    assert(bst.count_le(7) == N);
    assert(bst.sum() == 7LL * N);
    assert(bst.validate());

    // 全要素が同じ値のときの前後探索
    // (同じ値が何バケットにも跨るので、1 バケットだけ見て答えが出ることの確認)
    assert(!bst.next(7).has_value());
    assert(bst.next_eq(7).value() == 7);
    assert(!bst.prev(7).has_value());
    assert(bst.prev_eq(7).value() == 7);
    assert(bst.next(6).value() == 7);
    assert(bst.prev(8).value() == 7);

    for (int i = 0; i < N; i++) assert(bst.kth(i) == 7);
    assert(bst.sum_first_k(10) == 70);

    assert(bst.erase_all(7) == N);
    assert(bst.empty());
    assert(bst.validate());

    cout << "Same Value Test passed!" << endl;
}

void ascending_test() {
    cout << "Running Ascending Test..." << endl;

    OrderedMultiSet<64> bst;
    const int N = 200000;
    for (int i = 0; i < N; i++) bst.insert(i);

    assert(bst.size() == N);
    for (int i = 0; i < N; i++) assert(bst.kth(i) == i);
    assert(bst.min() == 0);
    assert(bst.max() == N - 1);
    for (int i = 0; i < N; i++) assert(bst.count(i) == 1);
    assert(bst.validate());
    // 先頭 k 個の総和 = 0+1+...+(k-1)
    for (int k = 0; k <= 1000; k++) assert(bst.sum_first_k(k) == (long long)k * (k - 1) / 2);

    cout << "Ascending Test passed!" << endl;
}

void descending_test() {
    cout << "Running Descending Test..." << endl;

    OrderedMultiSet<64> bst;
    const int N = 200000;
    for (int i = N - 1; i >= 0; i--) bst.insert(i);

    assert(bst.size() == N);
    for (int i = 0; i < N; i++) assert(bst.kth(i) == i);
    assert(bst.min() == 0);
    assert(bst.max() == N - 1);
    assert(bst.validate());

    cout << "Descending Test passed!" << endl;
}

void erase_all_elements_test() {
    cout << "Running Erase All Elements Test..." << endl;

    OrderedMultiSet<64> bst;
    const int N = 100000;
    for (int i = 0; i < N; i++) bst.insert(i);
    for (int i = 0; i < N; i++) assert(bst.erase_one(i));
    assert(bst.empty());
    assert(bst.validate());

    // 逆順に消す場合も確認 (バケットの併合方向が変わる)
    for (int i = 0; i < N; i++) bst.insert(i);
    for (int i = N - 1; i >= 0; i--) assert(bst.erase_one(i));
    assert(bst.empty());
    assert(bst.validate());

    cout << "Erase All Elements Test passed!" << endl;
}

void rank_kth_exhaustive_test() {
    cout << "Running Rank/Kth Exhaustive Test..." << endl;

    OrderedMultiSet<4> bst;
    multiset<long long> ms;
    mt19937 rng(0);

    for (int step = 0; step < 15000; step++) {
        int op = rng() % 3;
        int x = rng() % 10;

        if (op == 0) { bst.insert(x); ms.insert(x); }
        else if (op == 1) { bst.erase_one(x); auto it = ms.find(x); if (it != ms.end()) ms.erase(it); }
        else { bst.erase_all(x); ms.erase(x); }

        vector<long long> v(ms.begin(), ms.end());

        // rank / count 全探索
        for (int q = -2; q <= 12; q++) {
            assert(bst.count(q) == (long long)ms.count(q));
            assert(bst.count_lt(q) == distance(ms.begin(), ms.lower_bound(q)));
            assert(bst.count_le(q) == distance(ms.begin(), ms.upper_bound(q)));
            assert(bst.count_gt(q) == distance(ms.upper_bound(q), ms.end()));
            assert(bst.count_ge(q) == distance(ms.lower_bound(q), ms.end()));
            assert(bst.rank(q) == distance(ms.begin(), ms.lower_bound(q)));
            assert(bst.contains(q) == (ms.count(q) > 0));
        }

        // kth 全探索
        for (size_t k = 0; k < v.size(); k++) assert(bst.kth(k) == v[k]);
    }

    cout << "Rank/Kth Exhaustive Test passed!" << endl;
}

void kth_smallest_largest_test() {
    cout << "Running kth_smallest / kth_largest Test..." << endl;

    OrderedMultiSet<64> bst;
    multiset<long long> ms;
    mt19937 rng(12345);

    for (int step = 0; step < 30000; step++) {
        int op = rng() % 3;
        int x = rng() % 50;

        if (op == 0) { bst.insert(x); ms.insert(x); }
        else if (op == 1) { bst.erase_one(x); auto it = ms.find(x); if (it != ms.end()) ms.erase(it); }
        else { bst.erase_all(x); ms.erase(x); }

        vector<long long> v(ms.begin(), ms.end());

        for (size_t k = 0; k < v.size(); k++) {
            assert(bst.kth(k) == v[k]);
            assert(bst.kth_smallest(k) == v[k]);
            assert(bst.kth_largest(k) == v[v.size() - 1 - k]);
            assert(bst.kth_smallest(k) == bst.kth(k));
            assert(bst.kth_largest(k) == bst.kth(v.size() - 1 - k));
        }
    }

    cout << "kth_smallest / kth_largest Test passed!" << endl;
}

void kth_boundary_test() {
    cout << "Running kth Boundary Test..." << endl;

    OrderedMultiSet<64> bst;
    for (int i = 0; i < 1000; i++) bst.insert(i);

    assert(bst.kth_smallest(0) == 0);
    assert(bst.kth_smallest(999) == 999);
    assert(bst.kth_largest(0) == 999);
    assert(bst.kth_largest(999) == 0);

    // 範囲外は例外
    bool thrown = false;
    try { bst.kth(1000); } catch (const runtime_error&) { thrown = true; }
    assert(thrown);
    thrown = false;
    try { bst.kth(-1); } catch (const runtime_error&) { thrown = true; }
    assert(thrown);
    // 例外を投げない版
    assert(!bst.kth_opt(1000).has_value());

    cout << "kth Boundary Test passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト4] 総和系関数のランダムテスト (元からあるテスト)
// ---------------------------------------------------------------------------
void random_test() {
    cout << "--- 総和系関数のランダムテストを開始します ---" << endl;

    // バケットの分割・併合を頻発させるため、あえて BK=4 と小さくしてテスト
    OrderedMultiSet<4> st;
    multiset<long long> ms;
    mt19937_64 rng(42);

    for (int i = 0; i < 2000; ++i) {
        int op = rng() % 6;
        if (op == 0) {
            long long x = rng() % 100;
            st.insert(x);
            ms.insert(x);
        } else if (op == 1) {
            if (!ms.empty()) {
                auto it = ms.begin();
                advance(it, rng() % ms.size());
                long long x = *it;
                st.erase_one(x);
                ms.erase(ms.find(x));
            }
        } else {
            long long x = rng() % 100;
            long long naive_sum_lt = 0, naive_sum_le = 0, naive_sum_gt = 0, naive_sum_ge = 0;
            for (long long v : ms) {
                if (v < x) naive_sum_lt += v;
                if (v <= x) naive_sum_le += v;
                if (v > x) naive_sum_gt += v;
                if (v >= x) naive_sum_ge += v;
            }
            assert(st.sum_lt(x) == naive_sum_lt);
            assert(st.sum_le(x) == naive_sum_le);
            assert(st.sum_gt(x) == naive_sum_gt);
            assert(st.sum_ge(x) == naive_sum_ge);

            long long l = rng() % 100, r = rng() % 100;
            if (l > r) std::swap(l, r);
            long long naive_sum_range = 0;
            for (long long v : ms) if (v >= l && v < r) naive_sum_range += v;
            assert(st.sum_range(l, r) == naive_sum_range);
        }
    }
    cout << "ランダムテストが正常に完了しました（すべてのクエリで一致）。" << endl;
}

// ---------------------------------------------------------------------------
// [テスト5・新規] 負数・大きな値を含むランダムテスト
//   旧テストは値域が 0 以上に限られていたため、負数での総和・ランクを検証する。
// ---------------------------------------------------------------------------
template <int B>
void negative_value_test(const string& name, int steps) {
    cout << "Running Negative Value Test (" << name << ")..." << endl;
    OrderedMultiSet<B> bst;
    multiset<long long> ms;
    mt19937_64 rng(20240919);

    for (int i = 0; i < steps; i++) {
        long long x = (long long)(rng() % 41) - 20;          // -20 .. 20
        if (rng() % 100 < 15) x *= 1000000000LL;             // 時々 1e10 級の絶対値
        int op = rng() % 10;
        if (op < 5) { bst.insert(x); ms.insert(x); }
        else if (op < 8) { bool b = bst.erase_one(x); auto it = ms.find(x); assert(b == (it != ms.end())); if (it != ms.end()) ms.erase(it); }
        else { long long c = bst.erase_all(x); assert(c == (long long)ms.count(x)); ms.erase(x); }

        if (i % 7 == 0) verify_all_against_multiset(bst, ms, rng, 41);
        if (i % 200 == 0) assert(bst.validate());
    }
    assert(bst.validate());
    cout << "Negative Value Test (" << name << ") passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト6・新規] 空集合と境界の挙動
// ---------------------------------------------------------------------------
void empty_set_test() {
    cout << "Running Empty Set Test..." << endl;
    OrderedMultiSet<> bst;

    assert(bst.size() == 0);
    assert(bst.empty());
    assert(bst.sum() == 0);
    assert(bst.vec().empty());
    assert(bst.vec_range(-100, 100).empty());
    assert(bst.begin() == bst.end());
    assert(bst.lower_bound(0) == bst.end());
    assert(bst.upper_bound(0) == bst.end());
    assert(bst.count(5) == 0);
    assert(bst.count_lt(5) == 0 && bst.count_le(5) == 0);
    assert(bst.count_gt(5) == 0 && bst.count_ge(5) == 0);
    assert(bst.count_range(-10, 10) == 0);
    assert(bst.rank(5) == 0);
    assert(!bst.contains(5));
    assert(bst.sum_lt(5) == 0 && bst.sum_le(5) == 0 && bst.sum_gt(5) == 0 && bst.sum_ge(5) == 0);
    assert(bst.sum_range(0, 10) == 0);
    assert(bst.sum_first_k(3) == 0 && bst.sum_last_k(3) == 0 && bst.sum_index_range(0, 3) == 0);
    assert(bst.max_prefix_count(100) == 0);
    assert(!bst.next(0).has_value() && !bst.next_eq(0).has_value());
    assert(!bst.prev(0).has_value() && !bst.prev_eq(0).has_value());
    assert(!bst.front().has_value() && !bst.back().has_value());
    assert(!bst.pop_min().has_value() && !bst.pop_max().has_value());
    assert(!bst.erase_one(3));
    assert(bst.erase_all(3) == 0);
    assert(!bst.erase_at(0));
    assert(bst.erase_range(-5, 5) == 0);
    assert(bst.validate());

    bool thrown = false;
    try { bst.min(); } catch (const runtime_error&) { thrown = true; }
    assert(thrown);
    thrown = false;
    try { bst.max(); } catch (const runtime_error&) { thrown = true; }
    assert(thrown);
    thrown = false;
    try { bst.kth(0); } catch (const runtime_error&) { thrown = true; }
    assert(thrown);

    // 1 要素だけの状態
    bst.insert(42);
    assert(bst.size() == 1 && bst.min() == 42 && bst.max() == 42 && bst.sum() == 42);
    assert(bst.kth(0) == 42);
    assert(bst.prev_eq(42).value() == 42 && !bst.prev(42).has_value());
    assert(bst.next_eq(42).value() == 42 && !bst.next(42).has_value());
    assert(bst.validate());
    assert(bst.erase_one(42));
    assert(bst.empty() && bst.validate());

    // clear() の確認
    for (int i = 0; i < 1000; i++) bst.insert(i);
    bst.clear();
    assert(bst.empty() && bst.sum() == 0 && bst.validate());

    cout << "Empty Set Test passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト7・新規] 一括構築 build(vector) の検証
// ---------------------------------------------------------------------------
template <int B>
void build_test(const string& name) {
    cout << "Running Build Test (" << name << ")..." << endl;
    mt19937_64 rng(777);
    for (int n : {0, 1, 2, 3, B - 1, B, B + 1, 2 * B, 5 * B, 1000, 5000}) {
        if (n < 0) continue;
        vector<long long> v(n);
        for (auto& x : v) x = (long long)(rng() % 100) - 50;

        OrderedMultiSet<B> a(v);                       // 一括構築
        OrderedMultiSet<B> b;
        for (long long x : v) b.insert(x);             // 逐次挿入

        vector<long long> sorted_v = v;
        std::sort(sorted_v.begin(), sorted_v.end());

        assert(a.validate());
        assert(a.vec() == sorted_v);
        assert(a.vec() == b.vec());
        assert(a.size() == b.size());
        assert(a.sum() == b.sum());

        // 構築した木に対して挿入・削除を続けても壊れないこと
        for (int i = 0; i < 200; i++) {
            long long x = (long long)(rng() % 100) - 50;
            if (rng() % 2) { a.insert(x); b.insert(x); }
            else { assert(a.erase_one(x) == b.erase_one(x)); }
            assert(a.vec() == b.vec());
        }
        assert(a.validate());
    }
    cout << "Build Test (" << name << ") passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト8・新規] コピー / ムーブ / swap の検証
//   旧版はコピーコンストラクタが無く、コピーすると二重解放していた。
// ---------------------------------------------------------------------------
void copy_move_test() {
    cout << "Running Copy/Move Test..." << endl;
    mt19937_64 rng(31415);

    OrderedMultiSet<8> a;
    for (int i = 0; i < 2000; i++) a.insert((long long)(rng() % 500) - 250);
    const vector<long long> snapshot = a.vec();

    // コピーコンストラクタ: 独立していること
    OrderedMultiSet<8> b = a;
    assert(b.vec() == snapshot);
    assert(b.validate());
    for (int i = 0; i < 300; i++) b.insert(10000 + i);
    assert(a.vec() == snapshot);                 // 元は変化しない
    assert(b.size() == a.size() + 300);
    assert(a.validate() && b.validate());

    // コピー代入
    OrderedMultiSet<8> c;
    c.insert(-99999);
    c = a;
    assert(c.vec() == snapshot && c.validate());
    c = c;                                       // 自己代入
    assert(c.vec() == snapshot && c.validate());

    // ムーブ
    OrderedMultiSet<8> d = std::move(c);
    assert(d.vec() == snapshot && d.validate());
    assert(c.empty());                           // ムーブ元は空
    OrderedMultiSet<8> e;
    e = std::move(d);
    assert(e.vec() == snapshot && e.validate());

    // swap
    OrderedMultiSet<8> f;
    f.insert(1); f.insert(2);
    e.swap(f);
    assert(e.vec() == vector<long long>({1, 2}));
    assert(f.vec() == snapshot);
    assert(e.validate() && f.validate());

    cout << "Copy/Move Test passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト9・新規] 新規 API のランダムテスト
//   erase_at / erase_range / count_range / sum_first_k / sum_last_k /
//   sum_index_range / max_prefix_count / vec_range / イテレータ
// ---------------------------------------------------------------------------
template <int B>
void new_api_random_test(const string& name, int steps) {
    cout << "Running New API Random Test (" << name << ")..." << endl;
    OrderedMultiSet<B> bst;
    multiset<long long> ms;
    mt19937_64 rng(98765);

    for (int step = 0; step < steps; step++) {
        vector<long long> v(ms.begin(), ms.end());
        int op = rng() % 10;

        if (op < 5 || ms.empty()) {                      // 挿入 (非負値: max_prefix_count のため)
            long long x = (long long)(rng() % 60);
            bst.insert(x); ms.insert(x);
        } else if (op < 7) {                             // erase_at
            long long k = rng() % ms.size();
            auto it = ms.begin(); advance(it, k);
            assert(bst.erase_at(k));
            ms.erase(it);
        } else if (op < 8) {                             // pop_min / pop_max
            if (rng() % 2) { assert(bst.pop_min().value() == *ms.begin()); ms.erase(ms.begin()); }
            else { assert(bst.pop_max().value() == *ms.rbegin()); ms.erase(std::prev(ms.end())); }
        } else if (op < 9) {                             // erase_range
            long long l = rng() % 60, r = rng() % 60;
            if (l > r) std::swap(l, r);
            long long expect = 0;
            for (long long e : v) if (l <= e && e < r) expect++;
            assert(bst.erase_range(l, r) == expect);
            ms.erase(ms.lower_bound(l), ms.lower_bound(r));
        } else {                                         // erase_one
            long long x = (long long)(rng() % 60);
            auto it = ms.find(x);
            assert(bst.erase_one(x) == (it != ms.end()));
            if (it != ms.end()) ms.erase(it);
        }

        v.assign(ms.begin(), ms.end());
        assert(bst.vec() == v);

        // max_prefix_count: 非負列に対する貪欲個数
        {
            long long budget = (long long)(rng() % 400);
            long long expect = 0, acc = 0;
            for (long long e : v) { if (acc + e > budget) break; acc += e; expect++; }
            assert(bst.max_prefix_count(budget) == expect);
        }
        if (step % 5 == 0) verify_all_against_multiset(bst, ms, rng, 60);
        if (step % 200 == 0) assert(bst.validate());
    }
    assert(bst.validate());
    cout << "New API Random Test (" << name << ") passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト10・新規] T / S を変えたインスタンス化の確認
// ---------------------------------------------------------------------------
void custom_type_test() {
    cout << "Running Custom Type Test..." << endl;

    // T=int, S=long long : int の総和が 32bit を超えても正しい
    OrderedMultiSet<8, int, long long> a;
    long long expect = 0;
    for (int i = 0; i < 20000; i++) { a.insert(2000000000 - i); expect += 2000000000 - i; }
    assert(a.sum() == expect);
    assert(a.size() == 20000);
    assert(a.min() == 2000000000 - 19999);
    assert(a.max() == 2000000000);
    assert(a.kth(0) == 2000000000 - 19999);
    assert(a.validate());

    // T=double
    OrderedMultiSet<16, double, double> b;
    for (int i = 0; i < 100; i++) b.insert(i * 0.5);
    assert(b.size() == 100);
    assert(b.kth(3) == 1.5);
    assert(b.count_lt(2.0) == 4);
    assert(b.validate());

    cout << "Custom Type Test passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト11・新規] 複数の B について不変条件を毎回検査するストレス
// ---------------------------------------------------------------------------
template <int B>
void validate_stress_test(int steps) {
    OrderedMultiSet<B> bst;
    multiset<long long> ms;
    mt19937_64 rng(555 + B);
    for (int i = 0; i < steps; i++) {
        long long x = (long long)(rng() % 12) - 3;
        int op = rng() % 10;
        if (op < 5) { bst.insert(x); ms.insert(x); }
        else if (op < 8) { auto it = ms.find(x); assert(bst.erase_one(x) == (it != ms.end())); if (it != ms.end()) ms.erase(it); }
        else { assert(bst.erase_all(x) == (long long)ms.count(x)); ms.erase(x); }
        assert(bst.validate());
        assert(bst.vec() == vector<long long>(ms.begin(), ms.end()));
    }
}
void validate_stress_all() {
    cout << "Running Validate Stress Test (B=4,5,6,7,8,16)..." << endl;
    validate_stress_test<4>(8000);
    validate_stress_test<5>(8000);
    validate_stress_test<6>(8000);
    validate_stress_test<7>(8000);
    validate_stress_test<8>(8000);
    validate_stress_test<16>(8000);
    cout << "Validate Stress Test passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト12・新規] 大規模テスト兼 簡易ベンチマーク
//   ・N = 1e6 規模で壊れないこと / 実行時間の目安を確認する
//   ・正解との突き合わせは O(N) で可能なものだけ行う
// ---------------------------------------------------------------------------
void large_scale_test() {
    cout << "Running Large Scale Test (N = 1e6)..." << endl;
    const int N = 1000000;
    mt19937_64 rng(2026);

    // --- 一括構築 vs 逐次挿入の時間比較 ---
    vector<long long> a(N);
    for (auto& x : a) x = (long long)(rng() % 1000000000LL);

    auto t0 = chrono::steady_clock::now();
    OrderedMultiSet<64> bulk(a);
    auto t1 = chrono::steady_clock::now();

    OrderedMultiSet<64> seq;
    for (long long x : a) seq.insert(x);
    auto t2 = chrono::steady_clock::now();

    vector<long long> sorted_a = a;
    std::sort(sorted_a.begin(), sorted_a.end());
    assert(bulk.vec() == sorted_a);
    assert(seq.vec() == sorted_a);
    assert(bulk.sum() == seq.sum());
    assert(bulk.validate() && seq.validate());

    cout << "  build(vector) : "
         << chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << " ms" << endl;
    cout << "  insert x N    : "
         << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() << " ms" << endl;

    // --- クエリ時間の目安 ---
    auto t3 = chrono::steady_clock::now();
    long long checksum = 0;
    for (int i = 0; i < 1000000; i++) {
        long long x = (long long)(rng() % 1000000000LL);
        checksum += seq.count_lt(x);
    }
    auto t4 = chrono::steady_clock::now();
    for (int i = 0; i < 1000000; i++) checksum += seq.kth(rng() % N);
    auto t5 = chrono::steady_clock::now();
    for (int i = 0; i < 1000000; i++) checksum += seq.contains((long long)(rng() % 1000000000LL));
    auto t6 = chrono::steady_clock::now();
    cout << "  count_lt x 1e6: "
         << chrono::duration_cast<chrono::milliseconds>(t4 - t3).count() << " ms" << endl;
    cout << "  kth      x 1e6: "
         << chrono::duration_cast<chrono::milliseconds>(t5 - t4).count() << " ms" << endl;
    cout << "  contains x 1e6: "
         << chrono::duration_cast<chrono::milliseconds>(t6 - t5).count() << " ms" << endl;
    assert(checksum != 0);   // 最適化で消えないようにするためのダミー検査

    // --- 削除しきっても壊れないこと ---
    // (erase_all で一度に大量に消えるケースを通す)
    OrderedMultiSet<64> dup;
    for (int i = 0; i < 300000; i++) dup.insert(i % 100);
    long long total = 0;
    for (int v = 0; v < 100; v++) total += dup.erase_all(v);
    assert(total == 300000);
    assert(dup.empty() && dup.validate());

    cout << "Large Scale Test passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト13・新規] erase_index_range の徹底検証
//   バケットを丸ごと落とす経路があり、実装が最も込み入っているので単独で潰す。
//   ・1 バケットに収まるケース / 複数バケットに跨るケース
//   ・両端のバケットが空になるケース
//   ・範囲が集合の外にはみ出すケース (クランプ)
// ---------------------------------------------------------------------------
template <int BK>
void erase_index_range_test(const string& name) {
    cout << "Running erase_index_range Test (" << name << ")..." << endl;
    mt19937_64 rng(4649);
    for (int iter = 0; iter < 2000; iter++) {
        const int n = (int)(rng() % 400);
        vector<long long> base(n);
        for (auto& x : base) x = (long long)(rng() % 60) - 30;
        std::sort(base.begin(), base.end());
        OrderedMultiSet<BK> st(base);
        assert(st.validate());

        // わざと範囲外も混ぜる
        long long a = (long long)(rng() % (n + 5)) - 2;
        long long b = (long long)(rng() % (n + 5)) - 2;
        if (a > b) std::swap(a, b);

        const long long removed = st.erase_index_range(a, b);
        const long long ca = a < 0 ? 0 : a;
        const long long cb = b > n ? n : b;
        if (ca >= cb) {
            assert(removed == 0);
            assert(st.vec() == base);
        } else {
            assert(removed == cb - ca);
            vector<long long> expect(base.begin(), base.begin() + ca);
            expect.insert(expect.end(), base.begin() + cb, base.end());
            assert(st.vec() == expect);
            long long s = 0;
            for (long long e : expect) s += e;
            assert(st.sum() == s);
            assert(st.size() == (long long)expect.size());
        }
        assert(st.validate());
    }
    cout << "erase_index_range Test (" << name << ") passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト14・新規] erase_all / erase_range が複数バケットを跨いでも正しいこと
//   同じ値を大量に入れて、バケットを何個も跨がせた状態で消す。
// ---------------------------------------------------------------------------
template <int BK>
void bulk_erase_test(const string& name) {
    cout << "Running Bulk Erase Test (" << name << ")..." << endl;
    mt19937_64 rng(1220);
    for (int iter = 0; iter < 400; iter++) {
        // 値域を狭くして 1 つの値が何バケットにも跨るようにする
        const int n = 1 + (int)(rng() % (10 * BK + 20));
        vector<long long> base(n);
        for (auto& x : base) x = (long long)(rng() % 5);
        OrderedMultiSet<BK> st(base);
        multiset<long long> ms(base.begin(), base.end());
        assert(st.validate());

        if (rng() % 2) {                                   // erase_all
            const long long v = (long long)(rng() % 7) - 1;
            const long long c = st.erase_all(v);
            assert(c == (long long)ms.count(v));
            ms.erase(v);
        } else {                                           // erase_range
            long long lo = (long long)(rng() % 7) - 1, hi = (long long)(rng() % 7) - 1;
            if (lo > hi) std::swap(lo, hi);
            long long expect = 0;
            for (long long e : ms) if (lo <= e && e < hi) expect++;
            assert(st.erase_range(lo, hi) == expect);
            ms.erase(ms.lower_bound(lo), ms.lower_bound(hi));
        }
        assert(st.vec() == vector<long long>(ms.begin(), ms.end()));
        assert(st.validate());

        // 全部消しきっても壊れないこと
        for (long long v = -1; v <= 6; v++) st.erase_all(v);
        assert(st.empty() && st.validate());
    }
    cout << "Bulk Erase Test (" << name << ") passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト15・新規] insert_all / merge_from / split_ge
// ---------------------------------------------------------------------------
template <int BK>
void bulk_insert_split_test(const string& name) {
    cout << "Running insert_all / merge_from / split_ge Test (" << name << ")..." << endl;
    mt19937_64 rng(31337);
    for (int iter = 0; iter < 1500; iter++) {
        const int n = (int)(rng() % 300);
        const int m = (int)(rng() % 300);
        vector<long long> a(n), b(m);
        for (auto& x : a) x = (long long)(rng() % 80) - 40;
        for (auto& x : b) x = (long long)(rng() % 80) - 40;

        // --- insert_all ---
        {
            OrderedMultiSet<BK> st(a);
            st.insert_all(b);
            vector<long long> expect = a;
            expect.insert(expect.end(), b.begin(), b.end());
            std::sort(expect.begin(), expect.end());
            assert(st.vec() == expect);
            assert(st.size() == (long long)expect.size());
            long long s = 0; for (long long e : expect) s += e;
            assert(st.sum() == s);
            assert(st.validate());
            // 追加後も通常の操作が使えること
            for (int t = 0; t < 20; t++) {
                const long long x = (long long)(rng() % 80) - 40;
                st.insert(x);
                st.erase_one(x);
            }
            assert(st.vec() == expect);
            assert(st.validate());
        }

        // --- merge_from ---
        {
            OrderedMultiSet<BK> p(a), q(b);
            p.merge_from(q);
            vector<long long> expect = a;
            expect.insert(expect.end(), b.begin(), b.end());
            std::sort(expect.begin(), expect.end());
            assert(p.vec() == expect);
            assert(q.empty() && q.validate());       // 取り込まれた側は空になる
            assert(p.validate());
        }

        // --- split_ge ---
        {
            OrderedMultiSet<BK> st(a);
            const long long x = (long long)(rng() % 90) - 45;
            OrderedMultiSet<BK> hi = st.split_ge(x);
            vector<long long> sorted_a = a;
            std::sort(sorted_a.begin(), sorted_a.end());
            const size_t p = (size_t)(std::lower_bound(sorted_a.begin(), sorted_a.end(), x) - sorted_a.begin());
            vector<long long> lo_expect(sorted_a.begin(), sorted_a.begin() + p);
            vector<long long> hi_expect(sorted_a.begin() + p, sorted_a.end());
            assert(st.vec() == lo_expect);
            assert(hi.vec() == hi_expect);
            assert(st.size() == (long long)lo_expect.size());
            assert(hi.size() == (long long)hi_expect.size());
            long long sl = 0, sh = 0;
            for (long long e : lo_expect) sl += e;
            for (long long e : hi_expect) sh += e;
            assert(st.sum() == sl && hi.sum() == sh);
            assert(st.validate() && hi.validate());
            // 切り離した後も両方が普通に使えること
            for (int t = 0; t < 20; t++) {
                const long long y = (long long)(rng() % 80) - 40;
                st.insert(y); hi.insert(y);
            }
            assert(st.validate() && hi.validate());
            // merge_from で元に戻せること
            OrderedMultiSet<BK> back = st;
            back.merge_from(hi);
            assert(back.validate());
            assert(back.size() == st.size() + (long long)hi_expect.size() + 20);
        }
    }
    cout << "insert_all / merge_from / split_ge Test (" << name << ") passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト16・新規] 走査系の新 API
//   vec_index_range / for_each_range / kth_in_range / iterator_at /
//   lower_bound_index / upper_bound_index
// ---------------------------------------------------------------------------
template <int BK>
void traversal_api_test(const string& name) {
    cout << "Running Traversal API Test (" << name << ")..." << endl;
    mt19937_64 rng(7777);
    for (int iter = 0; iter < 1500; iter++) {
        const int n = (int)(rng() % 300);
        vector<long long> v(n);
        for (auto& x : v) x = (long long)(rng() % 50) - 25;
        std::sort(v.begin(), v.end());
        OrderedMultiSet<BK> st(v);

        // vec_index_range (範囲外はクランプ)
        {
            long long a = (long long)(rng() % (n + 5)) - 2;
            long long b = (long long)(rng() % (n + 5)) - 2;
            if (a > b) std::swap(a, b);
            const long long ca = std::max(0LL, a), cb = std::min<long long>(n, b);
            vector<long long> expect;
            if (ca < cb) expect.assign(v.begin() + ca, v.begin() + cb);
            assert(st.vec_index_range(a, b) == expect);
        }
        // for_each_range は vec_range と同じ結果になること
        {
            long long lo = (long long)(rng() % 60) - 30, hi = (long long)(rng() % 60) - 30;
            if (lo > hi) std::swap(lo, hi);
            vector<long long> got;
            st.for_each_range(lo, hi, [&](long long e) { got.push_back(e); });
            assert(got == st.vec_range(lo, hi));
            vector<long long> expect;
            for (long long e : v) if (lo <= e && e < hi) expect.push_back(e);
            assert(got == expect);
        }
        // kth_in_range
        {
            long long lo = (long long)(rng() % 60) - 30, hi = (long long)(rng() % 60) - 30;
            if (lo > hi) std::swap(lo, hi);
            vector<long long> seg;
            for (long long e : v) if (lo <= e && e < hi) seg.push_back(e);
            for (long long k = -1; k <= (long long)seg.size(); k++) {
                auto r = st.kth_in_range(lo, hi, k);
                if (k < 0 || k >= (long long)seg.size()) assert(!r.has_value());
                else assert(r.value() == seg[(size_t)k]);
            }
            assert(!st.kth_in_range(hi, lo, 0).has_value());   // 逆転区間
        }
        // iterator_at / lower_bound_index / upper_bound_index
        {
            for (long long k = 0; k < n; k++) assert(*st.iterator_at(k) == v[(size_t)k]);
            assert(st.iterator_at(n) == st.end());
            assert(st.iterator_at(-5) == (n == 0 ? st.end() : st.begin()));
            const long long x = (long long)(rng() % 60) - 30;
            assert(st.lower_bound_index(x) == (long long)(std::lower_bound(v.begin(), v.end(), x) - v.begin()));
            assert(st.upper_bound_index(x) == (long long)(std::upper_bound(v.begin(), v.end(), x) - v.begin()));
            // iterator_at(k) から end() まで辿ると残り全部が出てくること
            if (n > 0) {
                const long long k = (long long)(rng() % n);
                vector<long long> got;
                for (auto it = st.iterator_at(k); it != st.end(); ++it) got.push_back(*it);
                assert(got == vector<long long>(v.begin() + k, v.end()));
            }
        }
        // bucket_count は要素数に対して妥当な範囲であること
        {
            assert(st.bucket_count() >= (n > 0 ? 1 : 0));
            assert(st.bucket_count() <= n + 1);
        }
    }
    cout << "Traversal API Test (" << name << ") passed!" << endl;
}

// ---------------------------------------------------------------------------
// [テスト17・新規] バケット数が暴れないこと (大量削除しても痩せたバケットが残らない)
//   併合 (try_merge) が効いているかの確認。
// ---------------------------------------------------------------------------
void bucket_count_health_test() {
    cout << "Running Bucket Count Health Test..." << endl;
    constexpr int B = 128;
    OrderedMultiSet<B> st;
    mt19937_64 rng(20260920);
    const int N = 300000;
    vector<long long> a(N);
    for (auto& x : a) x = (long long)(rng() % 1000000000LL);
    for (long long x : a) st.insert(x);
    const int nb_full = st.bucket_count();
    // バケット数は N/BK 前後に収まっているはず
    assert(nb_full <= 2 * (N / B) + 4);
    // 9 割削除してもバケット数が残り要素数に見合う程度まで縮むこと
    for (int i = 0; i < N * 9 / 10; i++) assert(st.erase_one(a[i]));
    const long long rest = st.size();
    assert(rest == N - N * 9 / 10);
    assert(st.bucket_count() <= 2 * (rest / B) + 4);
    assert(st.validate());
    cout << "  要素 " << N << " → " << rest << " で バケット数 "
         << nb_full << " → " << st.bucket_count() << endl;
    cout << "Bucket Count Health Test passed!" << endl;
}

// ---------------------------------------------------------------------------
// [動作確認] 使い方の説明を兼ねたデモ
// ---------------------------------------------------------------------------
void test() {
    // 1. 構造体の宣言（既定 BK=128, T=long long, S=long long）
    //    BK は「1 バケットの目標長」。128〜192 なら実測差 1% 以内なので、
    //    普段は既定のままでよい。
    OrderedMultiSet<> st;

    cout << "--- 初期状態の確認 ---" << endl;
    cout << "空ですか？: " << (st.empty() ? "はい" : "いいえ") << endl;
    cout << "サイズ: " << st.size() << endl;
    cout << "要素の総和: " << st.sum() << endl;

    cout << "\n--- 要素の挿入 (重複要素を含む) ---" << endl;
    st.insert(10);
    st.insert(20);
    st.insert(20); // 重複
    st.insert(20); // 重複
    st.insert(30);
    st.insert(40);
    st.insert(40); // 重複

    auto print_vec = [&]() {
        cout << "現在の要素一覧 [vec()]: ";
        for (long long x : st.vec()) cout << x << " ";
        cout << endl;
    };
    print_vec();

    cout << "\n--- 基本状態の再確認 ---" << endl;
    cout << "空ですか？: " << (st.empty() ? "はい" : "いいえ") << endl;
    cout << "サイズ: " << st.size() << endl;
    cout << "要素の総和: " << st.sum() << endl;
    // min()/max() はバケットの両端そのものなので O(1) (旧 B+Tree 版は O(log N))
    cout << "最小値 [min()]: " << st.min() << endl;
    cout << "最大値 [max()]: " << st.max() << endl;
    // front()/back() は例外を投げない版 (空なら nullopt)
    cout << "最小値 [front()]: " << st.front().value() << endl;
    cout << "最大値 [back()]: " << st.back().value() << endl;

    cout << "\n--- カウント・順位系クエリ (ターゲット: 20) ---" << endl;
    cout << "20 を含んでいますか？ [contains()]: " << (st.contains(20) ? "はい" : "いいえ") << endl;
    cout << "50 を含んでいますか？ [contains()]: " << (st.contains(50) ? "はい" : "いいえ") << endl;
    cout << "20 の個数 [count()]: " << st.count(20) << endl;
    cout << "20 未満の位置 [rank()]: " << st.rank(20) << endl;
    cout << "20 未満の個数 [count_lt()]: " << st.count_lt(20) << endl;
    cout << "20 以下の個数 [count_le()]: " << st.count_le(20) << endl;
    cout << "20 より大きい個数 [count_gt()]: " << st.count_gt(20) << endl;
    cout << "20 以上の個数 [count_ge()]: " << st.count_ge(20) << endl;

    cout << "\n--- 【新規】値の区間クエリ [lo, hi) ---" << endl;
    // count_range は「値が [lo,hi) に入る要素数」。 rank の差分なので O(log N)。
    cout << "値が [20, 40) に入る個数 [count_range(20,40)]: " << st.count_range(20, 40)
         << " (想定: 20*3 + 30 = 4 個)" << endl;
    cout << "値が [20, 40) の要素一覧 [vec_range(20,40)]: ";
    for (long long x : st.vec_range(20, 40)) cout << x << " ";
    cout << endl;
    // equal_range_index は「x が占めるインデックス範囲 [first,last)」。std::equal_range のインデックス版。
    auto [er_l, er_r] = st.equal_range_index(20);
    cout << "20 が占めるインデックス範囲 [equal_range_index(20)]: ["
         << er_l << ", " << er_r << ")" << endl;

    cout << "\n--- 条件付き総和・区間和クエリ (ターゲット: 20) ---" << endl;
    cout << "20 未満 (strict <) の要素の総和 [sum_lt(20)]: " << st.sum_lt(20) << " (想定: 10)" << endl;
    cout << "20 以下 (<=) の要素の総和 [sum_le(20)]: " << st.sum_le(20) << " (想定: 10 + 20*3 = 70)" << endl;
    cout << "20 より大きい (strict >) 要素の総和 [sum_gt(20)]: " << st.sum_gt(20) << " (想定: 30 + 40*2 = 110)" << endl;
    cout << "20 以上 (>=) の要素の総和 [sum_ge(20)]: " << st.sum_ge(20) << " (想定: 20*3 + 30 + 40*2 = 170)" << endl;

    cout << "\n--- 区間和クエリ [l, r) ---" << endl;
    cout << "[20, 40) の範囲の要素の総和 [sum_range(20, 40)]: " << st.sum_range(20, 40) << " (想定: 20*3 + 30 = 90)" << endl;
    cout << "[15, 25) の範囲の要素の総和 [sum_range(15, 25)]: " << st.sum_range(15, 25) << " (想定: 20*3 = 60)" << endl;
    cout << "逆転区間 [40, 20) の総和 [sum_range(40, 20)]: " << st.sum_range(40, 20) << " (想定: 0)" << endl;

    cout << "\n--- 【新規】順位 (インデックス) を指定した総和 ---" << endl;
    // 現在の列は 10 20 20 20 30 40 40
    cout << "小さい方から 3 個の総和 [sum_first_k(3)]: " << st.sum_first_k(3) << " (想定: 10+20+20 = 50)" << endl;
    cout << "大きい方から 2 個の総和 [sum_last_k(2)]: " << st.sum_last_k(2) << " (想定: 40+40 = 80)" << endl;
    cout << "ソート順 [2, 5) 番目の総和 [sum_index_range(2,5)]: " << st.sum_index_range(2, 5) << " (想定: 20+20+30 = 70)" << endl;
    // max_prefix_count: 「安い方から買えるだけ買う」用途。非負要素であることが前提。
    cout << "予算 60 で小さい方から何個取れるか [max_prefix_count(60)]: " << st.max_prefix_count(60)
         << " (想定: 10+20+20=50 まで → 3 個)" << endl;

    cout << "\n--- インデックスアクセス系クエリ (0-indexed) ---" << endl;
    cout << "小さい方から 2 番目の要素 [kth(2)]: " << st.kth(2) << endl;
    cout << "小さい方から 3 番目の要素 [kth_smallest(3)]: " << st.kth_smallest(3) << endl;
    cout << "大きい方から 0 番目の要素 [kth_largest(0)]: " << st.kth_largest(0) << endl;
    cout << "大きい方から 1 番目の要素 [kth_largest(1)]: " << st.kth_largest(1) << endl;
    cout << "at(1) / operator[](1): " << st.at(1) << " / " << st[1] << endl;
    // 範囲外は kth() だと例外、kth_opt() だと nullopt
    cout << "範囲外 kth_opt(100) は値を持つか: " << (st.kth_opt(100).has_value() ? "はい" : "いいえ") << endl;

    cout << "\n--- 前後要素探索クエリ (optionalを返す関数) ---" << endl;
    auto print_optional = [](string name, optional<long long> opt) {
        cout << name << ": ";
        if (opt.has_value()) cout << opt.value() << endl;
        else cout << "nullopt (存在しません)" << endl;
    };

    print_optional("20 より大きい最小の要素 [next(20)]", st.next(20));
    print_optional("20 以上の最小の要素 [next_eq(20)]", st.next_eq(20));
    print_optional("20 より小さい最大の要素 [prev(20)]", st.prev(20));
    print_optional("20 以下の最大の要素 [prev_eq(20)]", st.prev_eq(20));
    print_optional("40 より大きい最小の要素 [next(40)]", st.next(40));
    print_optional("10 より小さい最大の要素 [prev(10)]", st.prev(10));

    cout << "\n--- 【新規】イテレータによる走査 ---" << endl;
    // begin()/end() は昇順の前方向イテレータ。範囲 for がそのまま使える。
    // 注意: 集合を変更するとイテレータは無効になる。
    cout << "範囲for: ";
    for (long long x : st) cout << x << " ";
    cout << endl;
    // lower_bound / upper_bound は std::multiset と同じ意味 (イテレータを返す)
    cout << "lower_bound(20) が指す値: " << *st.lower_bound(20) << endl;
    cout << "upper_bound(20) が指す値: " << *st.upper_bound(20) << endl;
    cout << "lower_bound(20) から末尾まで: ";
    for (auto it = st.lower_bound(20); it != st.end(); ++it) cout << *it << " ";
    cout << endl;

    cout << "\n--- 削除系関数の違いの検証 ---" << endl;
    print_vec();

    bool erased_one = st.erase_one(20);
    cout << "20 を1つ削除しましたか？ [erase_one(20)]: " << (erased_one ? "はい" : "いいえ") << endl;
    print_vec();
    cout << "現在の要素の総和: " << st.sum() << endl;

    long long erased_all_count = st.erase_all(20);
    cout << "20 をすべて削除しました。削除個数 [erase_all(20)]: " << erased_all_count << " 個" << endl;
    print_vec();
    cout << "現在の要素の総和: " << st.sum() << endl;

    bool erased_none = st.erase_one(99);
    cout << "存在しない 99 を1つ削除しようとしました [erase_one(99)]: " << (erased_none ? "成功" : "失敗") << endl;

    cout << "\n--- 【新規】インデックス指定削除・最小/最大の取り出し・値域削除 ---" << endl;
    print_vec();
    // erase_at(k): 小さい方から k 番目 (0-indexed) を削除
    st.erase_at(0);
    cout << "erase_at(0) 後: ";
    print_vec();
    // pop_min / pop_max: 取り出して削除 (空なら nullopt)
    print_optional("pop_min()", st.pop_min());
    print_optional("pop_max()", st.pop_max());
    cout << "pop 後: ";
    print_vec();

    cout << "\n--- 【新規】vector からの一括構築 ---" << endl;
    // N 回 insert するより高速 (ソート済みなら O(N))
    OrderedMultiSet<> st2(vector<long long>{5, 1, 4, 1, 5, 9, 2, 6});
    cout << "build 結果: ";
    for (long long x : st2.vec()) cout << x << " ";
    cout << endl;
    cout << "サイズ: " << st2.size() << " / 総和: " << st2.sum() << endl;
    // erase_range: 値が [lo,hi) の要素をすべて削除
    cout << "値が [2,6) の要素を削除 [erase_range(2,6)]: " << st2.erase_range(2, 6) << " 個削除" << endl;
    cout << "削除後: ";
    for (long long x : st2.vec()) cout << x << " ";
    cout << endl;

    cout << "\n--- 【新規】コピーは deep copy (元と独立) ---" << endl;
    OrderedMultiSet<> st3 = st2;   // コピーコンストラクタ
    st3.insert(100);
    cout << "コピー元: ";
    for (long long x : st2.vec()) cout << x << " ";
    cout << endl;
    cout << "コピー先(100を追加): ";
    for (long long x : st3.vec()) cout << x << " ";
    cout << endl;

    cout << "\n===== ここから 平方分割ベースで新しく入った機能 =====" << endl;

    cout << "\n--- 内部のバケット構成 ---" << endl;
    // 全要素を昇順に並べたものを長さ BK 程度のバケットに区切って持っている。
    // バケット数は要素数 / BK 程度。デバッグ時の目安に使える。
    OrderedMultiSet<4> small(vector<long long>{1, 3, 3, 7, 8, 8, 12, 15, 20, 20});
    cout << "BK=4 で 10 要素を入れたときのバケット数 [bucket_count()]: "
         << small.bucket_count() << endl;
    cout << "中身: ";
    for (int i = 0; i < small.bucket_count(); i++) {
        cout << "[";
        for (size_t j = 0; j < small.bk[i].size(); j++) cout << (j ? "," : "") << small.bk[i][j];
        cout << "]";
    }
    cout << endl;

    cout << "\n--- 【新規】insert_all / merge_from : 列や集合をまとめて取り込む ---" << endl;
    // 1 要素ずつ insert すると O(M log N) だが、まとめて入れるとマージして
    // 作り直すだけなので O(N + M)。
    OrderedMultiSet<> bulk(vector<long long>{1, 5, 9});
    bulk.insert_all({2, 2, 7});
    cout << "insert_all({2,2,7}) 後: ";
    for (long long x : bulk.vec()) cout << x << " ";
    cout << endl;

    OrderedMultiSet<> other2(vector<long long>{0, 100});
    bulk.merge_from(other2);      // other2 は空になる
    cout << "merge_from({0,100}) 後: ";
    for (long long x : bulk.vec()) cout << x << " ";
    cout << " / 取り込まれた側のサイズ: " << other2.size() << endl;

    cout << "\n--- 【新規】split_ge : 値 x 以上を丸ごと切り離す ---" << endl;
    // バケット列を付け替えるだけなので O(BK + バケット数)。要素数に比例しない。
    OrderedMultiSet<> hi = bulk.split_ge(7);
    cout << "split_ge(7) 後の自分 (7 未満): ";
    for (long long x : bulk.vec()) cout << x << " ";
    cout << endl;
    cout << "split_ge(7) の戻り値 (7 以上): ";
    for (long long x : hi.vec()) cout << x << " ";
    cout << endl;

    cout << "\n--- 【新規】erase_index_range : 順位で指定してまとめて削除 ---" << endl;
    // [a,b) 番目をまとめて削除する。間のバケットを丸ごと落とせるので、
    // 1 要素ずつ消すより大幅に速い。erase_all / erase_range もこれを使っている。
    OrderedMultiSet<> del(vector<long long>{10, 20, 30, 40, 50, 60});
    cout << "erase_index_range(1, 4) の削除個数: " << del.erase_index_range(1, 4) << endl;
    cout << "削除後: ";
    for (long long x : del.vec()) cout << x << " ";
    cout << endl;

    cout << "\n--- 【新規】vec_index_range / for_each_range / kth_in_range ---" << endl;
    OrderedMultiSet<> q(vector<long long>{1, 2, 2, 3, 5, 8, 13, 21});
    cout << "ソート順 [2,5) 番目 [vec_index_range(2,5)]: ";
    for (long long x : q.vec_index_range(2, 5)) cout << x << " ";
    cout << endl;
    // for_each_range は vector を作らずに値域 [lo,hi) を走査する
    cout << "値が [2,9) の要素 [for_each_range]: ";
    q.for_each_range(2, 9, [](long long e) { cout << e << " "; });
    cout << endl;
    // kth_in_range は「値が [lo,hi) の中で k 番目」
    cout << "値が [2,9) の中で 0 番目 [kth_in_range(2,9,0)]: " << q.kth_in_range(2, 9, 0).value() << endl;
    cout << "値が [2,9) の中で 3 番目 [kth_in_range(2,9,3)]: " << q.kth_in_range(2, 9, 3).value() << endl;
    cout << "値が [2,9) の中で 99 番目は存在するか: "
         << (q.kth_in_range(2, 9, 99).has_value() ? "はい" : "いいえ") << endl;

    cout << "\n--- 【新規】iterator_at / lower_bound_index / upper_bound_index ---" << endl;
    // iterator_at(k) は「k 番目の要素を指すイテレータ」。そこから前方へ走査できる。
    cout << "3 番目から末尾まで [iterator_at(3)]: ";
    for (auto it = q.iterator_at(3); it != q.end(); ++it) cout << *it << " ";
    cout << endl;
    // rank の別名。std::lower_bound / upper_bound の「インデックス版」であることを明示した名前。
    cout << "lower_bound_index(2) = " << q.lower_bound_index(2)
         << " / upper_bound_index(2) = " << q.upper_bound_index(2) << endl;

    cout << "\n--- 【新規】reserve : バケット配列の事前確保 ---" << endl;
    OrderedMultiSet<> pre;
    pre.reserve(1000000);       // 100 万要素ぶんのバケット配列を先に確保
    cout << "reserve(1000000) 済み" << endl;

    cout << "\n--- validate(): 不変条件チェック (デバッグ用) ---" << endl;
    // バケットの昇順性・bmax/bsum の一致・Fenwick の整合性・個数と総和 をすべて検査する
    cout << "validate() = " << (st3.validate() ? "true" : "false") << endl;
    cout << "validate() = " << (q.validate() ? "true" : "false") << endl;

    return;
}

void solv() {
    // ---- ランダムテスト -------------------------------------------------
    run_tests<64>(400000, 3000, "Standard Random Test");
    run_tests<4>(400000, 5, "Small Exhaustive Stress Test");
    stress_exhaustive_test();
    run_tests<4>(200000, 5, "B=4");
    run_tests<5>(200000, 7, "B=5 (奇数)");
    run_tests<8>(200000, 5, "B=8");
    run_tests<16>(200000, 5, "B=16");
    run_tests<64>(200000, 5, "B=64");
    validate_stress_all();
    negative_value_test<4>("B=4", 30000);
    negative_value_test<64>("B=64", 30000);
    new_api_random_test<4>("B=4", 15000);
    new_api_random_test<32>("B=32", 15000);
    build_test<4>("B=4");
    build_test<5>("B=5");
    build_test<64>("B=64");
    copy_move_test();
    custom_type_test();
    empty_set_test();
    same_value_test();
    ascending_test();
    descending_test();
    erase_all_elements_test();
    rank_kth_exhaustive_test();
    kth_smallest_largest_test();
    kth_boundary_test();
    random_test();
    // ---- バケット版で新しく入った API のテスト ----
    erase_index_range_test<4>("BK=4");
    erase_index_range_test<128>("BK=128");
    bulk_erase_test<4>("BK=4");
    bulk_erase_test<32>("BK=32");
    bulk_insert_split_test<4>("BK=4");
    bulk_insert_split_test<128>("BK=128");
    traversal_api_test<4>("BK=4");
    traversal_api_test<128>("BK=128");
    bucket_count_health_test();
    large_scale_test();
    cout << "\n===== ALL TESTS PASSED =====\n" << endl;

    // ---- 動作確認 -------------------------------------------------------
    test();
    return;
}

int main() {
    BUFF;
    int T = 1;
    //cin>>T;
    while (T--) solv();
    return 0;
}
