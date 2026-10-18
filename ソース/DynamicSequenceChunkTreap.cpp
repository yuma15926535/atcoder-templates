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
// 【設計見直しメモ】ベースを Implicit Treap → チャンク化 Treap (rope) に変更
// ----------------------------------------------------------------------------
// ■ 旧版の何が遅かったのか
//   旧版は「1 要素 = 1 節点」の Implicit Treap だった。この形だと
//     ・N=1e6 なら節点が 1e6 個。1 節点 336 byte (全部入り) で 300MB 超の領域を
//       ランダムに歩くことになり、木を 1 回降りるだけで 20〜40 回のキャッシュミス。
//     ・しかも木の降下は「今読んだ節点を見ないと次にどこを読むか決まらない」
//       依存ロードの鎖なので、ミスの待ち時間がそのまま直列に積み上がる。
//     ・1 要素ごとに pull() (集約の再計算) が走る。全部入りだと 1 回の pull に
//       gcd 2 回 + Kadane 12 要素の合成が入るので、これも段数ぶん掛かる。
//   つまり「木の段数 = log2(1e6) ≒ 20」がそのままコストの係数になっていた。
//
// ■ 新版: 1 節点に「最大 2*CS 要素の配列」を持たせる
//   節点を 1 要素ではなく CS 要素程度の配列 (チャンク) にすると、
//     ・節点数が N から N/CS に減る (CS=128 なら 1e6 → 約 7800 個)。
//     ・木の段数が log2(N) から log2(N/CS) に減る (20 → 13)。
//     ・チャンク内は連続メモリなので、同じ 1 回のキャッシュミスで
//       64 byte = 8 要素をまとめて持ってこられる (ハードウェアプリフェッチも効く)。
//   treap のままなので split / merge / 遅延伝播はすべてそのまま使え、
//   reverse を含めて計算量は O(log(N/CS)) を保てる。
//   「1 操作あたりに触る節点数 (= キャッシュミス数) を CS 分の 1 にする」のが狙い。
//
// ■ 捨てた案: 平方分割 (ブロックリスト)
//   最初に平方分割 (ブロック列 + ブロック長の Fenwick) を試した。
//   get / insert / sum は確かに最速だったが、reverse がブロック列そのものの
//   並べ替えになるため O(√N) になり、N が大きいと致命的に遅い。
//     N=1e6 (同一ワークロード, ms)   reverse
//       旧 treap                      5052
//       平方分割 (BS=450)            24788   ← 約 5 倍遅い
//       チャンク化 treap (CS=128)     3200
//   reverse は「機能を落とさない」条件に含まれるので平方分割は不採用とした。
//   チャンク化 treap は「平方分割のキャッシュ効率」と「treap の O(log N) reverse」を
//   両取りする形になっている。
//
// ■ 実測 (同一マシン・同一ワークロード、単位 ms。合計が小さいほど速い)
//   ワークロード: build → insert/erase/get を各 N 回 → sum / min+max を各 N 回
//                 → add / assign / reverse(ランダム長) / reverse(長さ64) を各 N 回
//                 → 長さ 64 の add+sum を N 回 → to_vector
//   (どの実装でも最後の sink 値が完全に一致することを毎回確認している)
//
//   ・AGG_BASIC (sum/min/max のみ, CS=128)
//       N       旧 treap    新チャンク    倍率
//       1e5       1700         810       2.10x
//       2e5       4608        1676       2.75x
//       5e5      15134        5297       2.86x
//       1e6      36049       12470       2.89x
//
//   ・AGG_ALL (全部入り, CS=32)
//       N       旧 treap    新チャンク    倍率
//       1e5       3224        1505       2.14x
//       2e5       7173        4199       1.71x
//       5e5      25151       12002       2.10x
//       1e6      58199       32071       1.81x
//
//   操作別 (N=1e6)                   AGG_BASIC            AGG_ALL
//     build (vector から構築)       87 →    5  (17x)    241 →   17  (14x)
//     insert                      2626 →  413  (6.4x)  5951 → 1584  (3.8x)
//     erase                       2696 →  637  (4.2x)  5895 → 1812  (3.3x)
//     get (1点取得)               1385 →  202  (6.9x)  1968 →  591  (3.3x)
//     sum(l,r)                    2565 →  659  (3.9x)  3647 → 1480  (2.5x)
//     min(l,r)+max(l,r)           3304 → 1111  (3.0x)  4011 → 2046  (2.0x)
//     add(l,r,x)                  4785 → 1407  (3.4x)  9681 → 3704  (2.6x)
//     assign(l,r,x)               4364 → 1188  (3.7x)  5995 → 3150  (1.9x)
//     reverse(ランダム長)         5052 → 3200  (1.6x)  7132 → 8025  (0.89x)
//     reverse(長さ64)             4492 → 2504  (1.8x)  5893 → 6795  (0.87x)
//     長さ64の add+sum            4601 → 1143  (4.0x)  7673 → 2862  (2.7x)
//   メモリも減る (節点数が 1/CS になるので、N=1e6 全部入りで概算 340MB → 40MB)。
//
// ■ どこで損をしているか (正直な注意点)
//   ・reverse だけは「AGG に AGG_GCD を含めると旧版と同等〜1 割遅い」。
//     reverse は境界のチャンクを 2 回割る必要があり、割った直後に
//     両半分の集約を作り直す (= 最大 2*CS 要素を舐め直す) ためで、
//     このとき gcd (ユークリッドの互除法) が 1 要素ごとに走るのが効く。
//     旧版は 1 節点 1 要素なので「割る」コストが無い代わりに段数で損をしており、
//     ちょうど釣り合ってしまう。
//     gcd を使わないなら AGG から AGG_GCD を外すだけで reverse も速くなる:
//       N=5e5, AGG_BITWISE|AGG_KADANE    旧 treap → 新チャンク
//         合計                            19178  →  7931   (2.42x)
//         reverse(ランダム長)              2605  →  2374   (1.10x)
//         reverse(長さ64)                  2298  →   725   (3.17x)
//   ・1 点の set / erase はチャンクの集約を作り直すので O(CS) が付く。
//     旧版の O(log N) より「計算量」としては悪いが、CS 要素の連続アクセスは
//     キャッシュミス 1 回の中で終わるため実測では常に速い (上表の通り)。
//   ・N が 1e4 を下回るような小さい列では、旧版との差はほとんど出ない
//     (もともと全体が L2 に乗るのでキャッシュミスが問題にならないため)。
//
// ■ この設計で効かせている工夫 (実測で効果を確認したもの)
//   (1) 遅延タグを 2 本に分ける
//       ct = 子へまだ配っていない遅延 / at = 自分の配列にまだ適用していない遅延。
//       1 本しか持たないと push() のたびに自分の配列 CS 要素を実体化することになり
//       O(CS) が全操作に付く。2 本に分けると push() は O(1) のままで、
//       配列の実体化 (O(CS)) は「実際に配列の中身に触るとき」だけに限定できる。
//       区間 add / assign / sum / min / max は「タグを見たまま」答えられるので、
//       実体化を最後まで遅らせられる。
//   (2) 区間更新は split/merge ではなく「降下してタグを付ける」
//       split/merge で挟むと節点が増えてチャンクが細切れになる。降下方式なら
//       節点数が変わらず、木を 1 回降りるだけで済む (add が 3.4 倍速くなった)。
//   (3) 短い並べ替えは「値を読んで書き戻す」
//       reverse / rotate は長さ 2*CS 以下なら split/merge を使わず、
//       値を取り出して並べ替えて同じ位置に書き戻す。木の形が一切変わらないので
//       チャンクが細切れにならず、再構築も誘発しない
//       (長さ 64 の reverse が AGG_BITWISE|AGG_KADANE で 3.2 倍速くなった)。
//   (4) 1 要素の挿入では集約を O(1) で更新する
//       sum/min/max/gcd/xor/and/or は「要素の集合」だけで決まるので、
//       1 要素増えただけなら配列を舐め直す必要がない (Kadane 系だけ順序に依存)。
//       これで AGG_BASIC の insert が 104ms → 65ms (N=2e5, 20 万回) になった。
//   (5) split で生まれる新しい節点は「一番外側で merge する」
//       分割点がチャンクの内部に落ちると、後半を入れる節点が新しく必要になる。
//       この節点の優先度は新しい乱数なので、再帰の途中でそのままぶら下げると
//       親より優先度が高くなり得て heap 性質 (親 >= 子) が壊れ、木が偏る。
//       そこで新節点は再帰から外へ持ち上げ、最後に merge で正しい深さへ挿し込む。
//       (実装当初はここを間違えており、validate() の heap 検査で検出した。
//        この検査はテストに残してある。)
//   (6) チャンクが痩せたら丸ごと組み直す
//       split を繰り返すとチャンクが細切れになるので、
//       「節点数 * CS > 3 * 要素数 + 8 * CS」になったら O(N) で組み直す。
//       組み直しは O(N) だが、そこまで細切れにするには Ω(N/CS) 回の
//       split が要るので償却 O(1)。この再構築を外すと reverse が 1.5 倍遅くなった。
//
// ■ CS (チャンク長) の既定値について
//   1 要素あたりの集約コストが重いほど、チャンクを短くした方が速い。
//   N=2e5 での合計時間 (小さいほど速い):
//     AGG             CS=32   CS=64   CS=128   既定
//     AGG_BASIC        1514    1381     1334    128
//     AGG_BITWISE      1920    1615     1797     64
//     AGG_KADANE       2170    1706     2015     64
//     AGG_GCD          2260    2451     4227     32
//     AGG_ALL          3042    3201     4807     32
//   この表のとおり「gcd あり → 32 / その他の集約あり → 64 / sum,min,max のみ → 128」
//   を既定にしている。第 3 テンプレート引数で明示的に変更できる
//   (例: DynamicSequence<long long, dseq::AGG_BASIC, 256>)。
//   CS=1 にすると「1 要素 = 1 節点」となり旧版と同じ形になる (テストで使っている)。
// ============================================================================

// ============================================================================
// DynamicSequence<T, AGG, CS>  --  AtCoder用 "全部入り" 動的配列テンプレート
// ----------------------------------------------------------------------------
// 実装方式 : チャンク化 Implicit Treap (rope)
//            乱択優先度つき2分探索木の各節点が「最大 2*CS 要素の配列」を持つ。
//            split / merge を基本操作とし, その組み合わせで全機能を実現する。
//
//   木のイメージ (CS=4 のとき。in-order に読むと元の列になる)
//
//                          [ 9, 2, 7 ]              ← 根の節点 (長さ 1..2*CS)
//                         /           .
//       左部分木 [ 4, 1, 8, 3 ]     [ 5, 6 ] 右部分木
//
//   in-order (左部分木 → 自分の配列 → 右部分木) に並べたものが実際の列:
//       4 1 8 3 | 9 2 7 | 5 6
//
//   ・各節点は「自分の配列だけの集約 blk」と「部分木全体の集約 sub」を持つ。
//   ・遅延は「子へ配る分 ct」と「自分の配列に未適用の分 at」に分けて持つ。
//     こうすると push() が O(1) のままで, 配列の実体化 (O(CS)) は
//     実際に配列の中身に触るときだけになる。
//
// 区間は常に半開区間 [l, r), 添字はすべて 0-based で統一する。
//
// ---------------------------------------------------------------------------
// 【テンプレート引数】
//
//   T   : 要素の型 (long long / int / __int128 など)
//   AGG : 持つ集約の種類 (ビット OR)。 ★高速化のために重要★
//   CS  : 1 節点が持つ配列の基準長 (省略時は AGG から自動で決まる)
//
//   DynamicSequence<long long>                    ... 全部入り (既定, 旧版と同じ)
//   DynamicSequence<long long, dseq::AGG_BASIC>   ... sum/min/max のみ (最速)
//   DynamicSequence<long long, dseq::AGG_KADANE>  ... sum/min/max + Kadane系
//   DynamicSequence<long long, dseq::AGG_BITWISE|dseq::AGG_GCD> ... など自由に OR
//   DynamicSequence<long long, dseq::AGG_BASIC, 256> ... CS も明示する場合
//
//   pull() (節点1個の集約再計算) は split/merge のたびに O(log(N/CS)) 回走る。
//   全部入りだと 1 回の pull で gcd (ユークリッドの互除法!) 2 回と
//   Kadane系 12 要素の合成が入るため, insert/erase/reverse などの
//   構造操作が数倍遅くなる。使わない集約を AGG から外すと,
//   ・pull の計算が消える
//   ・Node のメモリも減る (不要なフィールドは [[no_unique_address]] で消える)
//   ・CS の既定値が大きくなり, 節点数が更に減る
//   ので, 実測で構造操作が 3〜6 倍速くなる。
//   「sum しか使わない」問題では AGG_BASIC を強く推奨。
//
//   なお AGG から外した集約のクエリを呼ぶと static_assert でコンパイルエラーに
//   なるので, 事故る心配はない。
//
// ---------------------------------------------------------------------------
// 【API一覧と計算量】 (N = size(), 特記なき限り expected computational cost)
//   L = log(N/CS) と書く。 CS は定数なので O(L + CS) は O(log N) と同じオーダだが,
//   「木を降りる回数 (= キャッシュミス) は L 回だけ」「残りは連続アクセスの CS 要素」
//   という内訳を明示するためにこの表記にしている。
//
//  size / empty                          O(1)
//  get(i) / operator[](i) / at(i) / front / back  O(L)
//  set(i,x)                              O(L + CS)
//  push_back/pop_back/push_front/pop_front O(L + CS)
//  insert(i,x) / erase(i) / erase(l,r)   O(L + CS)
//  insert(i, vector) / append / prepend  O(len + L)
//  to_vector() / to_vector(l,r)          O(N) / O(r-l + L)
//  split(k) -> (左, 右)                   O(L + CS)
//  merge(other) (自分の末尾に連結)         O(L)      ※ static merge(a,b) も提供
//  swap(other)                           O(1)
//  clear()                               O(N/CS)
//  resize(n, val)                        O(|差分| + L)
//  sub(l,r) (コピーを取り出す)              O(r-l)
//  reverse(l,r)                          O(L + CS)  (長さ<=2*CS なら書き戻し方式)
//  reverse()  (全体)                     O(1)       (根にタグを付けるだけ)
//  rotate(l,r,k)                         O(L + CS)  (長さ<=2*CS なら書き戻し方式)
//  shift_left/shift_right(l,r,k)         O(L + CS)  (rotateの別名)
//  move(l,r,pos)                         O(L + CS)
//  swap_range(l1,r1,l2,r2)               O(L + CS)
//  extract(l,r) -> 新しい DynamicSequence  O(L + CS) (元から切り離す)
//  splice(pos,other,l,r)                 O(L + CS)
//  concat(a,b,...)                       O(L) per merge
//  set_range(l,r,vector)                 O(len + L)
//
//  --- Aggregate (sum/min/max は affine lazy と両立) ---
//  sum/min/max(l,r), 引数なし版          O(L + CS) / 全区間は O(1)
//  minmax(l,r)                           O(L + CS)
//  argmin(l,r) / argmax(l,r)             O(L + CS)
//  gcd_range/xor_sum/and_sum/or_sum(l,r) O(L + CS)  ただし後述の制約あり
//  max_prefix/suffix/min_prefix/suffix(l,r) O(L + CS) 制約あり (Kadane系)
//  max_subarray_sum/min_subarray_sum(l,r) O(L + CS) 制約あり (Kadane系)
//
//  --- Range Lazy ---
//  add(l,r,x) / assign(l,r,x) / multiply(l,r,x) / affine(l,r,a,b)  O(L + CS)
//    -> 内部では 1 種類の affine lazy (has_assign, mul, add) に統一して合成する。
//    -> split/merge を使わず降下してタグを付けるので節点が増えない。
//
//  --- 二分探索 ---
//  max_right(l, pred)                    O(L + CS)
//  min_left(r, pred)                     O(L + CS)
//  lower_bound_prefix(x)                 O(L + CS)
//  lower_bound_sorted(l,r,x) / upper_bound_sorted(l,r,x)  O(L log N)
//
//  --- 便利機能 (O(log N) は不可能, 明記の上で許容) ---
//  sort/stable_sort(l,r)      O(len log len)
//  unique(l,r)                O(len)   (sizeが変化する)
//  erase_if/filter(l,r,pred)  O(len)   (sizeが変化する)
//  partition/stable_partition(l,r,pred) O(len)
//  map(l,r,f)                 O(len)
//  kth_smallest(l,r,k) / kth_largest(l,r,k) / median(l,r) O(len)
//  count_value(l,r,x) / count_if(l,r,pred) / mode(l,r)    O(len) / O(len log len)
//  is_sorted(l,r)             O(len)
//  mex(l,r)                   O(len log len)
//  inversion(l,r)             O(len log len)  (クエリ毎に再計算, 差分維持はしない)
//
//  --- debug / 内部指標 ---
//  dump()/print()/check()/validate()   validate は集約値の整合性まで検査する
//  reserve_pool(n) : 節点プールの事前確保 (再確保コストを消せる)
//  chunk_size()    : CS の値                            ← 新規
//  node_count()    : 節点数 (≒ size()/CS)  O(1)         ← 新規
//  height()        : 木の高さ (≒ log2(節点数))  O(節点数) ← 新規
//    -> 「node_count() が size() の割に多い」「height() が極端に大きい」
//       ときだけ性能を疑えばよい。通常は気にしなくてよい。
//
// ---------------------------------------------------------------------------
// 【非可換/構造化 Aggregate と Lazy の両立性について (重要な設計方針)】
//
//  sum, min, max は affine変換 x -> a*x+b (add/assign/multiply/affine 全て
//  この形に統一される) のもとで O(1) で正しく更新できる。 (min/maxはa<0のとき
//  min/maxが入れ替わることに注意して実装している)
//
//  一方で以下の aggregate は変換の種類によって維持できる/できないが分かれるので,
//  3 種類のビット (BAD_GCD / BAD_BIT / BAD_KAD) に分離して,
//  維持できるものは維持するようにしている。
//
//    変換                | sum/min/max | gcd        | xor/and/or | Kadane系
//    --------------------+-------------+------------+------------+----------
//    assign (定数で埋める) | OK         | OK (再計算) | OK (再計算) | OK (再計算)
//    multiply (x -> a*x) | OK          | OK (|a|倍) | ダメ        | OK (a倍)
//    add (x -> x+b)      | OK          | ダメ       | ダメ        | ダメ
//    affine (a*x+b)      | OK          | ダメ       | ダメ        | ダメ
//
//  ・gcd は gcd(a*x1,...,a*xn) = |a| * gcd(x1,...,xn) なので定数倍に耐える。
//  ・Kadane系は「部分和 S が a 倍されるだけ」なので定数倍に耐える
//    (a<0 のとき max と min が入れ替わる)。
//  ・一方 add は部分和が S + k*b (k は長さ) に変わるため, 「どの長さが最適か」が
//    b によって非線形に変化して O(1) 更新できない。gcd も gcd(a+d,b+d) が
//    gcd(a,b) から決まらないのでダメ。xor/and/or は定数倍でも維持できない。
//
//  維持できない変換が掛かった区間に対してその集約を問い合わせた場合は
//  assert で検出する (黙って誤った値を返すことを避ける)。
//  add と gcd/xor/and/or/Kadane 集約を同時に必要とする問題では,
//  本テンプレートではなく問題特化の Segment Tree を検討すること。
//
// ---------------------------------------------------------------------------
// 【内部の不変条件】(validate() が全部検査する)
//   ・空でない限り, どの節点の配列も空ではなく, 長さは 2*CS 以下
//   ・heap 性質: 親の優先度 >= 子の優先度   ← 新版で追加した検査
//   ・cnt (部分木の要素数) が左部分木 + 自分の配列 + 右部分木 と一致
//   ・nn  (部分木の節点数) が 1 + 左 + 右 と一致  ← 新版で追加した検査
//   ・blk (自分の配列の集約) が配列から作り直した値と一致
//   ・sub (部分木の集約) が 左.sub + blk + 右.sub と一致
//   ・上記は遅延タグを実体化してから比較する (壊れている集約は bad ビットで除外)
//
// ---------------------------------------------------------------------------
// 【実装対象外の機能とその理由】
//
//  - value順 order statistic (値順の lower_bound(x) 等)
//      : 位置順の treap に加えて値順の平衡木を同時維持する必要があり,
//        2つの木を同期させる実装コストが非常に高い。
//        -> 値順の集合が必要なら OrderedMultiSet を併用すること。
//  - range kth_smallest / median を O(log N) や O(log^2 N) で : 一般には困難。
//    Wavelet Tree 等は静的構造で, insert/erase/split/merge に対応する
//    動的 Wavelet Tree は非常に重く実用性が低い。
//    -> 本テンプレートでは O(len) の "遅いが便利" 版のみ提供する。
//  - mode(l,r) (最頻値), count_value 全般高速化 : 一般値域での動的 range mode は
//    既知の効率的アルゴリズムが無く実用的でない。O(len) fallback のみ提供。
//  - mex(l,r) : 同上, 一般値域では O(len log len) fallback のみ。
//  - inversion(l,r) の "維持" (insert/erase/reverseと連動した差分更新)
//      : 一般には困難。本テンプレートはクエリ毎に O(len log len) で再計算する。
//  - range chmin/chmax (Segment Tree Beats 相当)
//      : Segment Tree Beats の O(log^2 N) 償却解析は「静的な完全二分木」を
//        前提とするポテンシャル論法に依存しており, treap の split/merge による
//        構造変化がポテンシャルを何度もリセットし得るため, 償却計算量の保証が
//        成立しない。一般には困難と分類し, 本テンプレートには含めない。
//  - 任意 predicate による O(log N) filter / arbitrary range sort
//      : 一般には不可能 (要素の再配置がΩ(len)必要なため)。
//        "ライブラリとして便利" という位置づけで O(len log len) 版のみ提供する。
//  - Weighted Sequence (value,weight) の専用 kth_by_weight() 等
//      : weight を要素に持たせ sum(l,r) と max_right() を組み合わせれば実現できる
//        ため, 専用APIとしては用意しない。
//  - Pair Aggregate (equal pair count 等)
//      : T をカスタム構造体にして operator+ 等を定義すれば汎用 sum で実現できる。
//
// ---------------------------------------------------------------------------
// 【T に要求される演算】
//   常に必要                : 既定構築, コピー, operator+ / += , operator< ,
//                             T(int) (サイズを掛けるため), operator* (assign/affine)
//   AGG_GCD を使う場合      : 単項 - , std::gcd が使える整数型
//   AGG_BITWISE を使う場合  : ^ & | (整数型)
//   AGG_KADANE を使う場合   : 数値型として自然に比較・加算できること
//   -> long long / int / __int128 などが想定。modint は AGG_BASIC 相当でのみ利用可
//      (大小比較が意味を持たないため min/max も使わないこと)。
//
// 【乱数・優先度】
//   mt19937_64 を chrono seed で初期化する (再現性が必要な場合は
//   rng の初期化を固定値に変更すること)。
//   優先度は 64bit のため衝突確率は無視できるほど小さい。
//
// 【メモリ】
//   全 DynamicSequence<T,AGG,CS> インスタンスは static な節点プールを共有する。
//   extract/splice は節点の所有権を付け替えるだけで, 要素をコピーしない。
//   解放された節点は free list で再利用する (配列の確保済み容量も再利用される)。
//   要素数が分かっているなら reserve_pool(n) を呼ぶと再確保が起きない
//   (n には「節点数の目安 = 要素数 / CS」を渡せば十分)。
//
// 【コピー】
//   節点はプール上の index で共有されるため, 単純なコピーは危険
//   (二重解放/意図しない共有) となる。よって copy constructor / copy
//   assignment operator は削除し, move のみ許可する。深いコピーが必要な
//   場合は明示的に clone() を呼ぶこと (O(N))。
//
// 【境界条件】
//   - 空配列, 長さ1, l==r (空区間, 何もしない), l==0, r==size, 全区間,
//     k==0, kが区間長以上 (mod処理), 負のrotate量, 先頭/末尾へのinsert/erase
//     について, 下記のとおり動作を定義し, assert で範囲外を検出する。
//   - チャンク境界 (長さ CS, 2*CS の前後) は CS=1..5 に落としたランダムテストで
//     総当たり的に検査している (下記テストA〜E)。
// ============================================================================
// 集約の種類を選ぶビットフラグ (AGG テンプレート引数に OR で渡す)
namespace dseq {
    constexpr int AGG_BASIC   = 0;   // sum / min / max のみ (最速・最小メモリ)
    constexpr int AGG_BITWISE = 1;   // xor_sum / and_sum / or_sum
    constexpr int AGG_GCD     = 2;   // gcd_range
    constexpr int AGG_KADANE  = 4;   // max/min_prefix, max/min_suffix, max/min_subarray_sum
    constexpr int AGG_ALL     = AGG_BITWISE | AGG_GCD | AGG_KADANE;
}

// CS は 1 節点が持つ配列の基準長 (実際の長さは 1 以上 2*CS 以下)。
// 1 要素あたりの集約コストが重いほど小さい方が速いので、既定値を AGG から決める。
//   gcd あり      -> 32   (gcd の計算が最も重い)
//   その他の集約  -> 64
//   sum/min/max   -> 128  (1 要素の処理が軽く、節点を減らす効果が勝つ)
template <class T, int AGG = dseq::AGG_ALL,
          int CS = ((AGG & dseq::AGG_GCD) ? 32 : (AGG == dseq::AGG_BASIC ? 128 : 64))>
struct DynamicSequence {
    static_assert(CS >= 1, "CS は 1 以上を指定してください");
    static constexpr bool HAS_BIT = (AGG & dseq::AGG_BITWISE) != 0;
    static constexpr bool HAS_GCD = (AGG & dseq::AGG_GCD) != 0;
    static constexpr bool HAS_KAD = (AGG & dseq::AGG_KADANE) != 0;

    // 「この集約が壊れている」ことを表すビット
    static constexpr unsigned char BAD_GCD = 1;   // gcd_range が無効
    static constexpr unsigned char BAD_BIT = 2;   // xor_sum / and_sum / or_sum が無効
    static constexpr unsigned char BAD_KAD = 4;   // Kadane系が無効

    // ------------------------------------------------------------------
    // 集約フィールドの有無を AGG で切り替えるためのサブ構造体。
    // 使わない集約は空構造体になり、[[no_unique_address]] でサイズが消える。
    // ------------------------------------------------------------------
    struct BitAgg { T xorv{}, andv{}, orv{}; };
    struct GcdAgg { T gcdv{}; };
    struct KadAgg { T pre_max{}, suf_max{}, pre_min{}, suf_min{}, best_max{}, best_min{}; };
    struct BitNone {};
    struct GcdNone {};
    struct KadNone {};
    using BitOpt = std::conditional_t<HAS_BIT, BitAgg, BitNone>;
    using GcdOpt = std::conditional_t<HAS_GCD, GcdAgg, GcdNone>;
    using KadOpt = std::conditional_t<HAS_KAD, KadAgg, KadNone>;

    // ------------------------------------------------------------------
    // Agg : 「ある連続した区間」の集約。
    // 部分木の集約にも、節点が持つ配列の集約にも、クエリの途中結果にも使う。
    // has=false は空区間 (単位元) を表す。
    // ------------------------------------------------------------------
    struct Agg {
        bool has = false;
        unsigned char bad = 0;
        T sum{}, mn{}, mx{};
        [[no_unique_address]] BitOpt bit{};
        [[no_unique_address]] GcdOpt gc{};
        [[no_unique_address]] KadOpt kd{};
    };

    // ------------------------------------------------------------------
    // Lazy : 区間に掛ける変換。
    //   has_assign なら「全部 assign_val にする」
    //   そうでなければ x -> mul*x + add
    //   rev なら並びを反転
    // (has_assign が true の間は mul=1, add=0 に保たれる)
    // ------------------------------------------------------------------
    struct Lazy {
        bool has_assign = false;
        bool rev = false;
        T assign_val{};
        T mul = T(1), add = T(0);
        bool empty() const { return !has_assign && !rev && mul == T(1) && add == T(0); }
    };

    // ------------------------------------------------------------------
    // Node : treap の節点。1 節点が「最大 2*CS 要素の配列」を持つ。
    //
    //   ・a          : この節点が受け持つ要素 (長さ 1 以上 2*CS 以下)
    //   ・blk        : a だけの集約
    //   ・sub        : 部分木 (左部分木 → a → 右部分木) 全体の集約
    //   ・ct         : 子へまだ配っていない遅延
    //   ・at         : a にまだ適用していない遅延
    //
    // 遅延を ct と at に分けているのがこの実装の肝。
    // こうすると push() が O(1) のままで、配列の実体化 (O(CS)) は
    // 「実際に配列の中身に触るとき」だけで済む。
    // ------------------------------------------------------------------
    struct Node {
        int ch[2] = {-1, -1};
        int cnt = 0;          // 部分木の要素数
        int nn = 1;           // 部分木の節点数 (再構築の判定に使う)
        unsigned long long pri = 0;
        Agg blk, sub;
        Lazy ct, at;
        vector<T> a;
    };

    static inline vector<Node> pool;
    static inline vector<int> free_list;
    static inline mt19937_64 rng{(unsigned long long)chrono::steady_clock::now().time_since_epoch().count()};

    int root = -1;

    // ==================================================================
    // 小物
    // ==================================================================
    static int cnt(int t) { return t == -1 ? 0 : pool[t].cnt; }
    static int nnode(int t) { return t == -1 ? 0 : pool[t].nn; }
    static T abs_val(const T& v) { return v < T(0) ? -v : v; }
    // メンバ関数 min/max/sort/swap/move 等が std:: の同名を隠すため、
    // クラス内では vmin/vmax を使い、他の標準関数は必ず std:: を明示する。
    template <class U> static U vmin(U a, U b) { return a < b ? a : b; }
    template <class U> static U vmax(U a, U b) { return a < b ? b : a; }

    static void reserve_pool(size_t n) { pool.reserve(n); }

    // ==================================================================
    // 集約の合成
    // ==================================================================
    // 1 要素の集約
    static Agg agg_single(const T& v) {
        Agg r;
        r.has = true; r.bad = 0;
        r.sum = v; r.mn = v; r.mx = v;
        if constexpr (HAS_GCD) r.gc.gcdv = abs_val(v);
        if constexpr (HAS_BIT) { r.bit.xorv = v; r.bit.andv = v; r.bit.orv = v; }
        if constexpr (HAS_KAD) {
            r.kd.pre_max = v; r.kd.suf_max = v; r.kd.best_max = v;
            r.kd.pre_min = v; r.kd.suf_min = v; r.kd.best_min = v;
        }
        return r;
    }
    // 区間 a の直後に区間 b をつないだものの集約
    static Agg agg_merge(const Agg& a, const Agg& b) {
        if (!a.has) return b;
        if (!b.has) return a;
        Agg r;
        r.has = true;
        r.bad = (unsigned char)(a.bad | b.bad);
        r.sum = a.sum + b.sum;
        r.mn = (a.mn < b.mn) ? a.mn : b.mn;
        r.mx = (a.mx < b.mx) ? b.mx : a.mx;
        if constexpr (HAS_GCD) {
            // どちらかが 1 なら結果は 1 (gcd は重いので早期に打ち切る)
            r.gc.gcdv = (a.gc.gcdv == T(1) || b.gc.gcdv == T(1)) ? T(1) : std::gcd(a.gc.gcdv, b.gc.gcdv);
        }
        if constexpr (HAS_BIT) {
            r.bit.xorv = a.bit.xorv ^ b.bit.xorv;
            r.bit.andv = a.bit.andv & b.bit.andv;
            r.bit.orv  = a.bit.orv  | b.bit.orv;
        }
        if constexpr (HAS_KAD) {
            // prefix は「a の中で完結」か「a 全部 + b の prefix」
            r.kd.pre_max = vmax(a.kd.pre_max, a.sum + b.kd.pre_max);
            r.kd.suf_max = vmax(b.kd.suf_max, b.sum + a.kd.suf_max);
            r.kd.best_max = vmax(vmax(a.kd.best_max, b.kd.best_max), a.kd.suf_max + b.kd.pre_max);
            r.kd.pre_min = vmin(a.kd.pre_min, a.sum + b.kd.pre_min);
            r.kd.suf_min = vmin(b.kd.suf_min, b.sum + a.kd.suf_min);
            r.kd.best_min = vmin(vmin(a.kd.best_min, b.kd.best_min), a.kd.suf_min + b.kd.pre_min);
        }
        return r;
    }
    // 配列 a の [lo,hi) の集約を 1 パスで計算する O(hi-lo)
    static Agg agg_range(const vector<T>& a, int lo, int hi) {
        Agg r;
        if (lo >= hi) return r;
        r.has = true; r.bad = 0;
        T s = a[lo], mn = a[lo], mx = a[lo];
        T g{}, xo{}, an{}, orr{};
        T pM{}, sM{}, bM{}, pm{}, sm{}, bm{};
        if constexpr (HAS_GCD) g = abs_val(a[lo]);
        if constexpr (HAS_BIT) { xo = a[lo]; an = a[lo]; orr = a[lo]; }
        if constexpr (HAS_KAD) { pM = sM = bM = a[lo]; pm = sm = bm = a[lo]; }
        for (int j = lo + 1; j < hi; j++) {
            const T v = a[j];
            if constexpr (HAS_KAD) {
                pM = vmax(pM, s + v);
                sM = vmax(v, sM + v);
                bM = vmax(bM, sM);
                pm = vmin(pm, s + v);
                sm = vmin(v, sm + v);
                bm = vmin(bm, sm);
            }
            s += v;
            if (v < mn) mn = v;
            if (mx < v) mx = v;
            if constexpr (HAS_GCD) { if (!(g == T(1))) g = std::gcd(g, abs_val(v)); }   // 1 になったら以降は 1
            if constexpr (HAS_BIT) { xo ^= v; an &= v; orr |= v; }
        }
        r.sum = s; r.mn = mn; r.mx = mx;
        if constexpr (HAS_GCD) r.gc.gcdv = g;
        if constexpr (HAS_BIT) { r.bit.xorv = xo; r.bit.andv = an; r.bit.orv = orr; }
        if constexpr (HAS_KAD) {
            r.kd.pre_max = pM; r.kd.suf_max = sM; r.kd.best_max = bM;
            r.kd.pre_min = pm; r.kd.suf_min = sm; r.kd.best_min = bm;
        }
        return r;
    }

    // ==================================================================
    // 遅延の合成と適用
    // ==================================================================
    // dst を適用した後に src を適用した結果を dst にまとめる
    static void compose(Lazy& dst, const Lazy& src) {
        if (src.has_assign) {
            dst.has_assign = true;
            dst.assign_val = src.assign_val;
            dst.mul = T(1); dst.add = T(0);
        } else if (!(src.mul == T(1)) || !(src.add == T(0))) {
            if (dst.has_assign) dst.assign_val = src.mul * dst.assign_val + src.add;
            else { dst.mul = src.mul * dst.mul; dst.add = src.mul * dst.add + src.add; }
        }
        dst.rev ^= src.rev;
    }
    // 長さ len の区間の集約 g に変換 z を適用する O(1)
    static void apply_agg(Agg& g, const Lazy& z, long long len) {
        if (!g.has) return;
        if (z.has_assign) {
            const T v = z.assign_val;
            g.sum = v * T(len); g.mn = v; g.mx = v;
            g.bad = 0;                                   // assign なら全ての集約を再計算できる
            if constexpr (HAS_GCD) g.gc.gcdv = abs_val(v);
            if constexpr (HAS_BIT) {
                g.bit.xorv = (len % 2 == 1) ? v : T(0);  // 偶数個の xor は 0
                g.bit.andv = v; g.bit.orv = v;
            }
            if constexpr (HAS_KAD) {
                // 全要素が v なので部分和は k*v。v>=0 なら長いほど大きい。
                if (!(v < T(0))) {
                    g.kd.pre_max = g.kd.suf_max = g.kd.best_max = v * T(len);
                    g.kd.pre_min = g.kd.suf_min = g.kd.best_min = v;
                } else {
                    g.kd.pre_max = g.kd.suf_max = g.kd.best_max = v;
                    g.kd.pre_min = g.kd.suf_min = g.kd.best_min = v * T(len);
                }
            }
        } else if (!(z.mul == T(1)) || !(z.add == T(0))) {
            const T a = z.mul, b = z.add;
            g.sum = a * g.sum + b * T(len);
            const bool neg = (a < T(0));                 // a<0 なら min と max が入れ替わる
            if (!neg) { g.mn = a * g.mn + b; g.mx = a * g.mx + b; }
            else { const T om = g.mn, oM = g.mx; g.mn = a * oM + b; g.mx = a * om + b; }
            const bool pure_mul = (b == T(0));           // 定数倍なら gcd と Kadane は維持できる
            // 「持っていない集約」の bad は立てない (立てると bad==0 を使った高速路が塞がる)
            if constexpr (HAS_BIT) g.bad |= BAD_BIT;     // xor/and/or は定数倍でも維持できない
            if constexpr (HAS_GCD) { if (!pure_mul) g.bad |= BAD_GCD; }
            if constexpr (HAS_KAD) { if (!pure_mul) g.bad |= BAD_KAD; }
            if constexpr (HAS_GCD) { if (pure_mul) g.gc.gcdv = abs_val(a) * g.gc.gcdv; }
            if constexpr (HAS_KAD) {
                if (pure_mul) {
                    if (!neg) {
                        g.kd.pre_max = a * g.kd.pre_max; g.kd.suf_max = a * g.kd.suf_max; g.kd.best_max = a * g.kd.best_max;
                        g.kd.pre_min = a * g.kd.pre_min; g.kd.suf_min = a * g.kd.suf_min; g.kd.best_min = a * g.kd.best_min;
                    } else {
                        const T pM = g.kd.pre_max, pm = g.kd.pre_min;
                        const T sM = g.kd.suf_max, sm = g.kd.suf_min;
                        const T bM = g.kd.best_max, bm = g.kd.best_min;
                        g.kd.pre_max = a * pm; g.kd.pre_min = a * pM;
                        g.kd.suf_max = a * sm; g.kd.suf_min = a * sM;
                        g.kd.best_max = a * bm; g.kd.best_min = a * bM;
                    }
                }
            }
        }
        if (z.rev) {
            if constexpr (HAS_KAD) {
                std::swap(g.kd.pre_max, g.kd.suf_max);   // 反転すると prefix と suffix が入れ替わる
                std::swap(g.kd.pre_min, g.kd.suf_min);
            }
        }
    }
    // 部分木 t 全体に変換 z を適用する O(1)
    static void apply(int t, const Lazy& z) {
        if (t == -1 || z.empty()) return;
        Node& nd = pool[t];
        apply_agg(nd.sub, z, nd.cnt);
        apply_agg(nd.blk, z, (long long)nd.a.size());
        if (z.rev) std::swap(nd.ch[0], nd.ch[1]);
        compose(nd.ct, z);
        compose(nd.at, z);
    }
    // 子へ伝播 O(1)
    static void push(int t) {
        Node& nd = pool[t];
        if (nd.ct.empty()) return;
        const Lazy z = nd.ct;
        nd.ct = Lazy();
        apply(nd.ch[0], z);
        apply(nd.ch[1], z);
    }
    // 自分の配列を実体化する O(CS) (配列の中身に直接触る直前に呼ぶ)
    static void mat(int t) {
        Node& nd = pool[t];
        if (nd.at.empty()) return;
        const Lazy z = nd.at;
        nd.at = Lazy();
        if (z.has_assign) { for (auto& e : nd.a) e = z.assign_val; }
        else if (!(z.mul == T(1)) || !(z.add == T(0))) { for (auto& e : nd.a) e = z.mul * e + z.add; }
        if (z.rev) std::reverse(nd.a.begin(), nd.a.end());
        // blk は apply_agg で更新済みなので作り直さなくてよい
    }
    // 配列を書き換えた後に呼ぶ (blk を作り直す) O(CS)
    static void pull_blk(int t) {
        Node& nd = pool[t];
        nd.blk = agg_range(nd.a, 0, (int)nd.a.size());
    }
    // 子と自分の配列から部分木の情報を作り直す O(1)
    static void pull(int t) {
        Node& nd = pool[t];
        const int l = nd.ch[0], r = nd.ch[1];
        nd.cnt = (int)nd.a.size() + cnt(l) + cnt(r);
        nd.nn = 1 + nnode(l) + nnode(r);
        Agg g = (l == -1) ? Agg{} : pool[l].sub;
        g = agg_merge(g, nd.blk);
        if (r != -1) g = agg_merge(g, pool[r].sub);
        nd.sub = g;
    }

    // ==================================================================
    // 節点の確保 / 解放
    // ==================================================================
    static int alloc_empty() {
        int idx;
        if (!free_list.empty()) { idx = free_list.back(); free_list.pop_back(); }
        else { pool.emplace_back(); idx = (int)pool.size() - 1; }
        Node& nd = pool[idx];
        nd.ch[0] = nd.ch[1] = -1;
        nd.cnt = 0; nd.nn = 1;
        nd.pri = rng();
        nd.ct = Lazy(); nd.at = Lazy();
        nd.a.clear();
        return idx;
    }
    // 配列を詰め終わった葉の情報を確定させる
    static void finish_leaf(int t) {
        pull_blk(t);
        Node& nd = pool[t];
        nd.cnt = (int)nd.a.size();
        nd.nn = 1;
        nd.sub = nd.blk;
    }
    static void free_subtree(int t) {
        if (t == -1) return;
        free_subtree(pool[t].ch[0]);
        free_subtree(pool[t].ch[1]);
        free_list.push_back(t);
    }

    // ==================================================================
    // split / merge
    // ==================================================================
    // split の本体。
    //
    // 【なぜ新節点を外に出すのか】
    // 分割点がチャンクの内部に落ちると「後半のチャンク」を入れる新しい節点が生まれる。
    // この節点の優先度は新しい乱数なので、木の途中にそのままぶら下げると
    // 親より高い優先度になって heap 性質 (親 >= 子) が壊れ、平衡が保証されなくなる。
    // そこで新節点は u_out として外に出し、再帰から戻った一番外側で merge する。
    // 新節点の要素は必ず「右側の列の先頭」に来るので merge(u, R) でよく、
    // merge が正しい位置 (優先度に従った深さ) に挿し込んでくれる。
    // 再帰の中で扱う節点は全て元からある節点なので、優先度は常に
    //   (戻り値の根) <= (その部分木の元の根) <= (呼び出し元の節点)
    // を満たし、heap 性質が保たれる。
    static pair<int,int> split_raw(int t, int k, int& u_out) {
        if (t == -1) return {-1, -1};
        push(t);
        const int lc = pool[t].ch[0];
        const int lcnt = cnt(lc);
        if (k <= lcnt) {
            auto [x, y] = split_raw(lc, k, u_out);
            pool[t].ch[0] = y; pull(t);
            return {x, t};
        }
        const int sz = (int)pool[t].a.size();
        if (k >= lcnt + sz) {
            auto [x, y] = split_raw(pool[t].ch[1], k - lcnt - sz, u_out);
            pool[t].ch[1] = x; pull(t);
            return {t, y};
        }
        // 分割点が自分の配列の内部 -> 配列を割り、後半を新しい節点にして外へ渡す
        mat(t);
        const int off = k - lcnt;
        const int u = alloc_empty();     // pool が再確保され得るので、以降は必ず添字で触ること
        pool[u].a.assign(pool[t].a.begin() + off, pool[t].a.end());
        pool[t].a.resize(off);
        finish_leaf(u);
        pull_blk(t);
        const int rc = pool[t].ch[1];
        pool[t].ch[1] = -1;
        pull(t);
        u_out = u;
        return {t, rc};
    }
    // 先頭 k 要素と残りに分割する。 O(log(N/CS) + CS)
    static pair<int,int> split(int t, int k) {
        int u = -1;
        auto [x, y] = split_raw(t, k, u);
        // 新節点は右側の列の先頭に来る。merge が cnt / 集約も経路上で直してくれる。
        if (u != -1) y = merge(u, y);
        return {x, y};
    }
    // a (前) と b (後) を連結する。 O(log(N/CS))
    static int merge(int a, int b) {
        if (a == -1) return b;
        if (b == -1) return a;
        if (pool[a].pri > pool[b].pri) {
            push(a);
            pool[a].ch[1] = merge(pool[a].ch[1], b);
            pull(a);
            return a;
        }
        push(b);
        pool[b].ch[0] = merge(a, pool[b].ch[0]);
        pull(b);
        return b;
    }

    // ==================================================================
    // 構築
    // ==================================================================
    // vector から O(N) で構築する。
    // CS 個ずつのチャンクに切り、そのチャンク列に対して
    // 「ランダム優先度のデカルト木」を単調スタックで組む。
    static int build(const vector<T>& v) { return build(v, 0, (int)v.size()); }
    static int build(const vector<T>& v, int l, int r) {
        if (l >= r) return -1;
        const int n = r - l;
        const int g = (n + CS - 1) / CS;
        vector<int> ids;
        ids.reserve(g);
        int pos = l;
        for (int j = 0; j < g; j++) {
            const int len = n / g + (j < n % g ? 1 : 0);
            const int id = alloc_empty();
            pool[id].a.assign(v.begin() + pos, v.begin() + pos + len);
            pos += len;
            finish_leaf(id);
            ids.push_back(id);
        }
        // 右スパインを単調スタックで管理してデカルト木を作る (O(チャンク数))
        vector<int> st;
        st.reserve(g);
        for (int i = 0; i < g; i++) {
            const int cur = ids[i];
            int last = -1;
            while (!st.empty() && pool[st.back()].pri < pool[cur].pri) {
                last = st.back(); st.pop_back();
                pull(last);
            }
            pool[cur].ch[0] = last;
            if (!st.empty()) pool[st.back()].ch[1] = cur;
            st.push_back(cur);
        }
        while (st.size() > 1) { const int u = st.back(); st.pop_back(); pull(u); }
        pull(st[0]);
        return st[0];
    }
    // 分割でチャンクが細切れになりすぎたら丸ごと組み直す (償却 O(1))
    void maybe_rebuild() {
        if (root == -1) return;
        if ((long long)pool[root].nn * CS <= 3LL * pool[root].cnt + 8LL * CS) return;
        vector<T> v = to_vector();
        free_subtree(root);
        root = build(v);
    }

    // ==================================================================
    // コンストラクタ / デストラクタ / move
    // ==================================================================
    DynamicSequence() = default;
    explicit DynamicSequence(const vector<T>& a) { root = build(a); }
    explicit DynamicSequence(int n, const T& val) { root = build(vector<T>(n, val)); }
    DynamicSequence(const DynamicSequence&) = delete;
    DynamicSequence& operator=(const DynamicSequence&) = delete;
    DynamicSequence(DynamicSequence&& other) noexcept : root(other.root) { other.root = -1; }
    DynamicSequence& operator=(DynamicSequence&& other) noexcept {
        if (this != &other) { free_subtree(root); root = other.root; other.root = -1; }
        return *this;
    }
    ~DynamicSequence() { free_subtree(root); }

    // 深いコピー。 O(N)
    DynamicSequence clone() const {
        DynamicSequence res;
        res.root = clone_subtree(root);
        return res;
    }
    static int clone_subtree(int t) {
        if (t == -1) return -1;
        push(t); mat(t);
        const int idx = alloc_empty();
        pool[idx].a = pool[t].a;
        pool[idx].pri = pool[t].pri;
        const int lc = clone_subtree(pool[t].ch[0]);
        const int rc = clone_subtree(pool[t].ch[1]);
        pool[idx].ch[0] = lc;
        pool[idx].ch[1] = rc;
        pull_blk(idx);
        pull(idx);
        return idx;
    }

    void clear() { free_subtree(root); root = -1; }
    void reset(const vector<T>& a) { clear(); root = build(a); }
    void swap(DynamicSequence& other) noexcept { std::swap(root, other.root); }

    // ==================================================================
    // 基本情報
    // ==================================================================
    int size() const { return cnt(root); }
    bool empty() const { return root == -1; }

    // ==================================================================
    // 要素の取得・更新
    // ==================================================================
    // O(log(N/CS))
    T get(int i) const {
        assert(0 <= i && i < size());
        int t = root;
        while (true) {
            push(t);
            const int lcnt = cnt(pool[t].ch[0]);
            const int sz = (int)pool[t].a.size();
            if (i < lcnt) { t = pool[t].ch[0]; continue; }
            if (i >= lcnt + sz) { i -= lcnt + sz; t = pool[t].ch[1]; continue; }
            const int off = i - lcnt;
            const Lazy& z = pool[t].at;
            if (z.has_assign) return z.assign_val;
            const int p = z.rev ? sz - 1 - off : off;
            return z.mul * pool[t].a[p] + z.add;
        }
    }
    T operator[](int i) const { return get(i); }
    T at(int i) const { return get(i); }
    T front() const { assert(!empty()); return get(0); }
    T back() const { assert(!empty()); return get(size() - 1); }

    // O(log(N/CS) + CS)
    void set(int i, const T& x) {
        assert(0 <= i && i < size());
        set_impl(root, i, x);
    }
    static void set_impl(int t, int i, const T& x) {
        push(t);
        const int lcnt = cnt(pool[t].ch[0]);
        const int sz = (int)pool[t].a.size();
        if (i < lcnt) set_impl(pool[t].ch[0], i, x);
        else if (i >= lcnt + sz) set_impl(pool[t].ch[1], i - lcnt - sz, x);
        else {
            mat(t);
            pool[t].a[i - lcnt] = x;
            pull_blk(t);
        }
        pull(t);
    }

    // ==================================================================
    // 挿入・削除
    // ==================================================================
    // i の位置に x を挿入する。 O(log(N/CS) + CS)
    void insert(int i, const T& x) {
        assert(0 <= i && i <= size());
        if (root == -1) {
            const int id = alloc_empty();
            pool[id].a.push_back(x);
            finish_leaf(id);
            root = id;
            return;
        }
        // 木を 1 回降りて、該当チャンクの配列に直接挿し込む (split/merge を使わない)。
        long long need = -1;
        ins_impl(root, i, x, 0, need);
        if (need >= 0) {
            // チャンクが 2*CS を超えたので、その中央で構造的に分割する。
            // ここで split を使うのは、新しい節点を優先度に従った正しい深さへ
            // 挿し込めるのが merge (= split の後処理) だけだから。
            // 発生頻度は「CS 回の挿入につき 1 回」程度なので、
            // 1 要素あたりの追加コストは O(log(N/CS) / CS) で無視できる。
            auto [a, b] = split(root, (int)need);
            root = merge(a, b);
        }
        maybe_rebuild();
    }
    // base は t の部分木の先頭が列全体の何番目かを表す (分割位置を絶対位置で返すため)
    static void ins_impl(int t, int i, const T& x, int base, long long& need) {
        push(t);
        const int lcnt = cnt(pool[t].ch[0]);
        const int sz = (int)pool[t].a.size();
        if (i < lcnt) { ins_impl(pool[t].ch[0], i, x, base, need); pull(t); return; }
        if (i > lcnt + sz) { ins_impl(pool[t].ch[1], i - lcnt - sz, x, base + lcnt + sz, need); pull(t); return; }
        mat(t);
        pool[t].a.insert(pool[t].a.begin() + (i - lcnt), x);
        // sum/min/max/gcd/xor/and/or は「要素の集合」だけで決まるので、
        // 1 要素増やすだけなら O(1) で更新できる (配列を舐め直さなくてよい)。
        // Kadane 系だけは並び順に依存するので作り直す。
        if constexpr (HAS_KAD) pull_blk(t);
        else if (pool[t].blk.bad) pull_blk(t);   // 壊れた集約は作り直せるなら作り直す
        else pool[t].blk = agg_merge(pool[t].blk, agg_single(x));
        pull(t);
        // 2*CS を超えたら呼び出し元に「ここで割って」と伝える (一時的に 2*CS+1 になる)
        if ((int)pool[t].a.size() > 2 * CS) need = (long long)base + lcnt + (int)pool[t].a.size() / 2;
    }
    // i の位置に列 v をまとめて挿入する。 O(|v| + log(N/CS))
    void insert(int i, const vector<T>& v) {
        assert(0 <= i && i <= size());
        if (v.empty()) return;
        auto [a, b] = split(root, i);
        root = merge(merge(a, build(v)), b);
        maybe_rebuild();
    }
    void append(const vector<T>& v) { insert(size(), v); }
    void prepend(const vector<T>& v) { insert(0, v); }

    void erase(int i) {
        assert(0 <= i && i < size());
        root = ers_impl(root, i);
        maybe_rebuild();   // 削除でチャンクが痩せ続けたときに組み直す (償却 O(1))
    }
    static int ers_impl(int t, int i) {
        push(t);
        const int lcnt = cnt(pool[t].ch[0]);
        const int sz = (int)pool[t].a.size();
        if (i < lcnt) { pool[t].ch[0] = ers_impl(pool[t].ch[0], i); pull(t); return t; }
        if (i >= lcnt + sz) { pool[t].ch[1] = ers_impl(pool[t].ch[1], i - lcnt - sz); pull(t); return t; }
        mat(t);
        pool[t].a.erase(pool[t].a.begin() + (i - lcnt));
        if (pool[t].a.empty()) {
            const int res = merge(pool[t].ch[0], pool[t].ch[1]);
            free_list.push_back(t);
            return res;
        }
        pull_blk(t);
        pull(t);
        return t;
    }
    // [l,r) を削除。 O(log(N/CS) + CS)
    void erase(int l, int r) {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return;
        auto [a, bc] = split(root, l);
        auto [b, c] = split(bc, r - l);
        free_subtree(b);
        root = merge(a, c);
        maybe_rebuild();
    }

    void push_back(const T& x) { insert(size(), x); }
    void pop_back() { assert(!empty()); erase(size() - 1); }
    void push_front(const T& x) { insert(0, x); }
    void pop_front() { assert(!empty()); erase(0); }

    void resize(int n, const T& val = T()) {
        assert(n >= 0);
        const int cur = size();
        if (n < cur) erase(n, cur);
        else if (n > cur) append(vector<T>(n - cur, val));
    }

    // ==================================================================
    // to_vector
    // ==================================================================
    vector<T> to_vector() const { return to_vector(0, size()); }
    vector<T> to_vector(int l, int r) const {
        assert(0 <= l && l <= r && r <= size());
        vector<T> res;
        res.reserve(r - l);
        collect(root, l, r, res);
        return res;
    }
    static void collect(int t, int l, int r, vector<T>& res) {
        if (t == -1 || l >= r) return;
        push(t);
        const int lcnt = cnt(pool[t].ch[0]);
        const int sz = (int)pool[t].a.size();
        if (l < lcnt) collect(pool[t].ch[0], l, vmin(r, lcnt), res);
        const int a0 = vmax(l, lcnt) - lcnt, a1 = vmin(r, lcnt + sz) - lcnt;
        if (a0 < a1) {
            mat(t);
            res.insert(res.end(), pool[t].a.begin() + a0, pool[t].a.begin() + a1);
        }
        if (r > lcnt + sz) collect(pool[t].ch[1], vmax(0, l - lcnt - sz), r - lcnt - sz, res);
    }

    // ==================================================================
    // split / merge (公開API)
    // ==================================================================
    pair<DynamicSequence, DynamicSequence> split(int k) {
        assert(0 <= k && k <= size());
        auto [a, b] = split(root, k);
        root = -1;
        DynamicSequence L, R;
        L.root = a; R.root = b;
        L.maybe_rebuild(); R.maybe_rebuild();
        return {std::move(L), std::move(R)};
    }
    void merge(DynamicSequence& other) {
        root = merge(root, other.root);
        other.root = -1;
    }
    static DynamicSequence merge(DynamicSequence& a, DynamicSequence& b) {
        DynamicSequence res;
        res.root = merge(a.root, b.root);
        a.root = -1; b.root = -1;
        return res;
    }
    static DynamicSequence concat(vector<DynamicSequence>& seqs) {
        DynamicSequence res;
        for (auto& s : seqs) res.merge(s);
        return res;
    }

    // ==================================================================
    // 構造操作
    // ==================================================================
    // 区間 [l,r) の値を v で上書きする (木の形は一切変えない)。
    // 触ったチャンクだけ実体化して書き込むので O(r-l + (かかるチャンク数)*CS)。
    static void overwrite_rec(int t, int l, int r, const vector<T>& v, int& pos) {
        if (t == -1 || l >= r) return;
        push(t);
        const int lcnt = cnt(pool[t].ch[0]);
        const int sz = (int)pool[t].a.size();
        if (l < lcnt) overwrite_rec(pool[t].ch[0], l, vmin(r, lcnt), v, pos);
        const int a0 = vmax(l, lcnt) - lcnt, a1 = vmin(r, lcnt + sz) - lcnt;
        if (a0 < a1) {
            mat(t);
            for (int j = a0; j < a1; j++) pool[t].a[j] = v[pos++];
            pull_blk(t);
        }
        if (r > lcnt + sz) overwrite_rec(pool[t].ch[1], vmax(0, l - lcnt - sz), r - lcnt - sz, v, pos);
        pull(t);
    }
    // 「短い区間の並べ替え」は split/merge を使わずに値の書き戻しで済ませる。
    // split/merge は節点を増やすので、短い区間に対して何度も呼ぶと
    // 木が細切れになって再構築が頻発する。書き戻しなら節点数が変わらない。
    static constexpr int SMALL_PERM = 2 * CS;   // これ以下の長さなら書き戻し方式
    template <class F>
    void permute_small(int l, int r, F f) {
        vector<T> v = to_vector(l, r);
        f(v);
        int pos = 0;
        overwrite_rec(root, l, r, v, pos);
    }
    void reverse(int l, int r) {
        assert(0 <= l && l <= r && r <= size());
        if (r - l <= 1) return;
        if (l == 0 && r == size()) { reverse(); return; }   // 全体反転は根にタグを付けるだけ O(1)
        if (r - l <= SMALL_PERM) {                          // 短い区間は値を反転して書き戻す
            permute_small(l, r, [](vector<T>& v) { std::reverse(v.begin(), v.end()); });
            return;
        }
        auto [a, bc] = split(root, l);
        auto [b, c] = split(bc, r - l);
        Lazy z; z.rev = true;
        apply(b, z);
        root = merge(merge(a, b), c);
        maybe_rebuild();
    }
    void reverse() { if (root != -1) { Lazy z; z.rev = true; apply(root, z); } }

    void rotate(int l, int r, long long k) {
        assert(0 <= l && l <= r && r <= size());
        const int len = r - l;
        if (len <= 1) return;
        k %= len; if (k < 0) k += len;
        if (k == 0) return;
        if (len <= SMALL_PERM) {                            // 短い区間は値を回して書き戻す
            permute_small(l, r, [&](vector<T>& v) { std::rotate(v.begin(), v.begin() + (len - (int)k), v.end()); });
            return;
        }
        auto [a, bc] = split(root, l);
        auto [b, c] = split(bc, len);
        auto [b1, b2] = split(b, len - (int)k);
        root = merge(merge(a, merge(b2, b1)), c);
        maybe_rebuild();
    }
    void shift_right(int l, int r, long long k) { rotate(l, r, k); }
    void shift_left(int l, int r, long long k) { rotate(l, r, -k); }

    void move(int l, int r, int pos) {
        assert(0 <= l && l <= r && r <= size());
        const int len = r - l;
        assert(0 <= pos && pos <= size() - len);
        if (len == 0) return;
        auto [a, bc] = split(root, l);
        auto [b, c] = split(bc, len);
        const int reduced = merge(a, c);
        auto [r1, r2] = split(reduced, pos);
        root = merge(merge(r1, b), r2);
        maybe_rebuild();
    }

    void swap_range(int l1, int r1, int l2, int r2) {
        assert(0 <= l1 && l1 <= r1);
        assert(0 <= l2 && l2 <= r2);
        if (l1 > l2) { std::swap(l1, l2); std::swap(r1, r2); }
        assert(r1 <= l2 && "swap_range: 区間が重複しています");
        assert(r2 <= size());
        if (r1 - l1 == 0 && r2 - l2 == 0) return;
        auto [a, rest1] = split(root, l1);
        auto [seg1, rest2] = split(rest1, r1 - l1);
        auto [b, rest3] = split(rest2, l2 - r1);
        auto [seg2, c] = split(rest3, r2 - l2);
        root = merge(merge(merge(merge(a, seg2), b), seg1), c);
        maybe_rebuild();
    }

    DynamicSequence extract(int l, int r) {
        assert(0 <= l && l <= r && r <= size());
        auto [a, bc] = split(root, l);
        auto [b, c] = split(bc, r - l);
        root = merge(a, c);
        DynamicSequence res;
        res.root = b;
        maybe_rebuild();
        res.maybe_rebuild();
        return res;
    }
    void splice(int pos, DynamicSequence& other, int l, int r) {
        assert(0 <= pos && pos <= size());
        assert(0 <= l && l <= r && r <= other.size());
        auto [oa, obc] = split(other.root, l);
        auto [seg, oc] = split(obc, r - l);
        other.root = merge(oa, oc);
        auto [a, b] = split(root, pos);
        root = merge(merge(a, seg), b);
        maybe_rebuild();
        other.maybe_rebuild();
    }
    DynamicSequence sub(int l, int r) const {
        assert(0 <= l && l <= r && r <= size());
        DynamicSequence res;
        res.root = build(to_vector(l, r));
        return res;
    }
    void set_range(int l, int r, const vector<T>& v) {
        assert(0 <= l && l <= r && r <= size());
        rebuild_range(l, r, v);
    }

    // ==================================================================
    // Aggregate : sum / min / max (専用の軽量クエリ)
    // ==================================================================
    static T sum_rec(int t, int l, int r) {
        if (t == -1 || l >= r) return T(0);
        if (l <= 0 && r >= pool[t].cnt) return pool[t].sub.sum;
        push(t);
        const int lcnt = cnt(pool[t].ch[0]);
        const int sz = (int)pool[t].a.size();
        T res = T(0);
        if (l < lcnt) res += sum_rec(pool[t].ch[0], l, vmin(r, lcnt));
        const int a0 = vmax(l, lcnt) - lcnt, a1 = vmin(r, lcnt + sz) - lcnt;
        if (a0 < a1) {
            if (a0 == 0 && a1 == sz) res += pool[t].blk.sum;
            else {
                const Lazy& z = pool[t].at;
                if (z.has_assign) res += z.assign_val * T(a1 - a0);
                else {
                    const int s = z.rev ? sz - a1 : a0, e = z.rev ? sz - a0 : a1;
                    T p = T(0);
                    for (int j = s; j < e; j++) p += pool[t].a[j];
                    res += z.mul * p + z.add * T(a1 - a0);
                }
            }
        }
        if (r > lcnt + sz) res += sum_rec(pool[t].ch[1], vmax(0, l - lcnt - sz), r - lcnt - sz);
        return res;
    }
    static void minmax_rec(int t, int l, int r, T& mn, T& mx, bool& has) {
        if (t == -1 || l >= r) return;
        auto take = [&](const T& a, const T& b) {
            if (!has) { mn = a; mx = b; has = true; }
            else { if (a < mn) mn = a; if (mx < b) mx = b; }
        };
        if (l <= 0 && r >= pool[t].cnt) { take(pool[t].sub.mn, pool[t].sub.mx); return; }
        push(t);
        const int lcnt = cnt(pool[t].ch[0]);
        const int sz = (int)pool[t].a.size();
        if (l < lcnt) minmax_rec(pool[t].ch[0], l, vmin(r, lcnt), mn, mx, has);
        const int a0 = vmax(l, lcnt) - lcnt, a1 = vmin(r, lcnt + sz) - lcnt;
        if (a0 < a1) {
            if (a0 == 0 && a1 == sz) take(pool[t].blk.mn, pool[t].blk.mx);
            else {
                const Lazy& z = pool[t].at;
                if (z.has_assign) take(z.assign_val, z.assign_val);
                else {
                    const int s = z.rev ? sz - a1 : a0, e = z.rev ? sz - a0 : a1;
                    T lo = pool[t].a[s], hi = pool[t].a[s];
                    for (int j = s + 1; j < e; j++) { const T v = pool[t].a[j]; if (v < lo) lo = v; if (hi < v) hi = v; }
                    const T p = z.mul * lo + z.add, q = z.mul * hi + z.add;
                    take(vmin(p, q), vmax(p, q));
                }
            }
        }
        if (r > lcnt + sz) minmax_rec(pool[t].ch[1], vmax(0, l - lcnt - sz), r - lcnt - sz, mn, mx, has);
    }

    T sum(int l, int r) const { assert(0 <= l && l <= r && r <= size()); return sum_rec(root, l, r); }
    T sum() const { return root == -1 ? T(0) : pool[root].sub.sum; }
    T min(int l, int r) const {
        assert(0 <= l && l < r && r <= size());
        T mn{}, mx{}; bool has = false;
        minmax_rec(root, l, r, mn, mx, has);
        return mn;
    }
    T max(int l, int r) const {
        assert(0 <= l && l < r && r <= size());
        T mn{}, mx{}; bool has = false;
        minmax_rec(root, l, r, mn, mx, has);
        return mx;
    }
    T min() const { assert(!empty()); return pool[root].sub.mn; }
    T max() const { assert(!empty()); return pool[root].sub.mx; }
    pair<T,T> minmax(int l, int r) const {
        assert(0 <= l && l < r && r <= size());
        T mn{}, mx{}; bool has = false;
        minmax_rec(root, l, r, mn, mx, has);
        return {mn, mx};
    }

    // [l,r) 内で最小値/最大値をとる最も左の位置。 O(log(N/CS) + CS)
    int argmin(int l, int r) const {
        assert(0 <= l && l < r && r <= size());
        const T target = min(l, r);
        const int res = arg_rec(root, l, r, target, true);
        assert(res != -1);
        return res;
    }
    int argmax(int l, int r) const {
        assert(0 <= l && l < r && r <= size());
        const T target = max(l, r);
        const int res = arg_rec(root, l, r, target, false);
        assert(res != -1);
        return res;
    }
    static int arg_rec(int t, int l, int r, const T& target, bool find_min) {
        if (t == -1 || l >= r) return -1;
        if (l <= 0 && r >= pool[t].cnt) {
            if (find_min) { if (target < pool[t].sub.mn) return -1; }
            else { if (pool[t].sub.mx < target) return -1; }
        }
        push(t);
        const int lcnt = cnt(pool[t].ch[0]);
        const int sz = (int)pool[t].a.size();
        if (l < lcnt) {
            const int res = arg_rec(pool[t].ch[0], l, vmin(r, lcnt), target, find_min);
            if (res != -1) return res;
        }
        const int a0 = vmax(l, lcnt) - lcnt, a1 = vmin(r, lcnt + sz) - lcnt;
        if (a0 < a1) {
            bool skip = false;
            if (a0 == 0 && a1 == sz) {
                if (find_min) skip = (target < pool[t].blk.mn);
                else skip = (pool[t].blk.mx < target);
            }
            if (!skip) {
                mat(t);
                for (int j = a0; j < a1; j++) {
                    const T& v = pool[t].a[j];
                    if (!(v < target) && !(target < v)) return lcnt + j;
                }
            }
        }
        if (r > lcnt + sz) {
            const int res = arg_rec(pool[t].ch[1], vmax(0, l - lcnt - sz), r - lcnt - sz, target, find_min);
            if (res != -1) return res + lcnt + sz;
        }
        return -1;
    }

    // ==================================================================
    // Aggregate : gcd / xor / and / or / Kadane系
    // ==================================================================
    static Agg query(int t, int l, int r) {
        if (t == -1 || l >= r) return Agg{};
        if (l <= 0 && r >= pool[t].cnt) return pool[t].sub;
        push(t);
        const int lcnt = cnt(pool[t].ch[0]);
        const int sz = (int)pool[t].a.size();
        Agg res;
        if (l < lcnt) res = query(pool[t].ch[0], l, vmin(r, lcnt));
        const int a0 = vmax(l, lcnt) - lcnt, a1 = vmin(r, lcnt + sz) - lcnt;
        if (a0 < a1) {
            if (a0 == 0 && a1 == sz) res = agg_merge(res, pool[t].blk);
            else { mat(t); res = agg_merge(res, agg_range(pool[t].a, a0, a1)); }
        }
        if (r > lcnt + sz) res = agg_merge(res, query(pool[t].ch[1], vmax(0, l - lcnt - sz), r - lcnt - sz));
        return res;
    }

    T gcd_range(int l, int r) const {
        static_assert(HAS_GCD, "gcd_range を使うには AGG に dseq::AGG_GCD を含めてください");
        assert(0 <= l && l < r && r <= size());
        Agg a = query(root, l, r);
        assert(!(a.bad & BAD_GCD) && "gcd_range: この区間には add/affine が適用されており gcd を維持できません");
        return a.gc.gcdv;
    }
    T xor_sum(int l, int r) const {
        static_assert(HAS_BIT, "xor_sum を使うには AGG に dseq::AGG_BITWISE を含めてください");
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return T(0);
        Agg a = query(root, l, r);
        assert(!(a.bad & BAD_BIT) && "xor_sum: この区間には add/multiply/affine が適用されており維持できません");
        return a.bit.xorv;
    }
    T and_sum(int l, int r) const {
        static_assert(HAS_BIT, "and_sum を使うには AGG に dseq::AGG_BITWISE を含めてください");
        assert(0 <= l && l < r && r <= size());
        Agg a = query(root, l, r);
        assert(!(a.bad & BAD_BIT) && "and_sum: この区間には add/multiply/affine が適用されており維持できません");
        return a.bit.andv;
    }
    T or_sum(int l, int r) const {
        static_assert(HAS_BIT, "or_sum を使うには AGG に dseq::AGG_BITWISE を含めてください");
        assert(0 <= l && l < r && r <= size());
        Agg a = query(root, l, r);
        assert(!(a.bad & BAD_BIT) && "or_sum: この区間には add/multiply/affine が適用されており維持できません");
        return a.bit.orv;
    }
    T max_prefix(int l, int r) const { return kadane_query(l, r).kd.pre_max; }
    T max_suffix(int l, int r) const { return kadane_query(l, r).kd.suf_max; }
    T min_prefix(int l, int r) const { return kadane_query(l, r).kd.pre_min; }
    T min_suffix(int l, int r) const { return kadane_query(l, r).kd.suf_min; }
    T max_subarray_sum(int l, int r) const { return kadane_query(l, r).kd.best_max; }
    T min_subarray_sum(int l, int r) const { return kadane_query(l, r).kd.best_min; }
    Agg kadane_query(int l, int r) const {
        static_assert(HAS_KAD, "Kadane系を使うには AGG に dseq::AGG_KADANE を含めてください");
        assert(0 <= l && l < r && r <= size());
        Agg a = query(root, l, r);
        assert(!(a.bad & BAD_KAD) && "Kadane系: この区間には add/affine が適用されており維持できません");
        return a;
    }

    // ==================================================================
    // Range Lazy : add / assign / multiply / affine
    // ==================================================================
    // 区間更新は「降下して完全に覆う部分木・配列にタグを付ける」方式。
    // split/merge を使わないので節点が増えず、木を 1 回降りるだけで済む。
    static void update_rec(int t, int l, int r, const Lazy& z) {
        if (t == -1 || l >= r) return;
        if (l <= 0 && r >= pool[t].cnt) { apply(t, z); return; }
        push(t);
        const int lcnt = cnt(pool[t].ch[0]);
        const int sz = (int)pool[t].a.size();
        if (l < lcnt) update_rec(pool[t].ch[0], l, vmin(r, lcnt), z);
        const int a0 = vmax(l, lcnt) - lcnt, a1 = vmin(r, lcnt + sz) - lcnt;
        if (a0 < a1) {
            if (a0 == 0 && a1 == sz) {
                // 自分の配列を丸ごと更新 (実体化は不要)
                apply_agg(pool[t].blk, z, sz);
                compose(pool[t].at, z);
            } else {
                mat(t);
                auto& a = pool[t].a;
                if (z.has_assign) { for (int j = a0; j < a1; j++) a[j] = z.assign_val; }
                else { for (int j = a0; j < a1; j++) a[j] = z.mul * a[j] + z.add; }
                pull_blk(t);
            }
        }
        if (r > lcnt + sz) update_rec(pool[t].ch[1], vmax(0, l - lcnt - sz), r - lcnt - sz, z);
        pull(t);
    }
    void add(int l, int r, const T& x) { range_affine(l, r, T(1), x); }
    void assign(int l, int r, const T& x) { range_assign(l, r, x); }
    void multiply(int l, int r, const T& x) { range_affine(l, r, x, T(0)); }
    void affine(int l, int r, const T& a, const T& b) { range_affine(l, r, a, b); }
    void range_assign(int l, int r, const T& x) {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return;
        Lazy z; z.has_assign = true; z.assign_val = x;
        update_rec(root, l, r, z);
    }
    void range_affine(int l, int r, const T& mulv, const T& addv) {
        assert(0 <= l && l <= r && r <= size());
        if (l == r) return;
        if (mulv == T(0)) { range_assign(l, r, addv); return; }   // 0 倍は実質 assign
        if (mulv == T(1) && addv == T(0)) return;                  // 恒等変換
        Lazy z; z.mul = mulv; z.add = addv;
        update_rec(root, l, r, z);
    }

    // ==================================================================
    // 二分探索
    // ==================================================================
    template <class F>
    int max_right(int l, F pred) const {
        assert(0 <= l && l <= size());
        assert(pred(T(0)) && "max_right: pred(0) が true である必要があります");
        if (l == size()) return l;
        T acc = T(0);
        return max_right_rec(root, l, 0, pred, acc);
    }
    template <class F>
    static int max_right_rec(int t, int l, int offset, F& pred, T& acc) {
        if (t == -1) return offset;
        if (l <= offset) {
            const T cand = acc + pool[t].sub.sum;
            if (pred(cand)) { acc = cand; return offset + pool[t].cnt; }
        }
        push(t);
        const int lc = pool[t].ch[0], rc = pool[t].ch[1];
        const int lcnt = cnt(lc);
        const int sz = (int)pool[t].a.size();
        const int mid = offset + lcnt;                // 自分の配列の先頭位置
        if (l < mid) {
            const int res = max_right_rec(lc, l, offset, pred, acc);
            if (res != mid) return res;
        }
        if (l < mid + sz) {
            // 自分の配列を 1 要素ずつ試す (丸ごと入るならまとめて判定)
            const int from = vmax(l, mid) - mid;
            if (from == 0) {
                const T cand = acc + pool[t].blk.sum;
                if (pred(cand)) { acc = cand; goto go_right; }
            }
            mat(t);
            for (int j = from; j < sz; j++) {
                const T cand = acc + pool[t].a[j];
                if (!pred(cand)) return mid + j;
                acc = cand;
            }
        }
        go_right:
        return max_right_rec(rc, l, mid + sz, pred, acc);
    }

    template <class F>
    int min_left(int r, F pred) const {
        assert(0 <= r && r <= size());
        assert(pred(T(0)) && "min_left: pred(0) が true である必要があります");
        if (r == 0) return 0;
        T acc = T(0);
        return min_left_rec(root, r, 0, pred, acc);
    }
    template <class F>
    static int min_left_rec(int t, int r, int offset, F& pred, T& acc) {
        if (t == -1) return offset;
        if (offset + pool[t].cnt <= r) {
            const T cand = acc + pool[t].sub.sum;
            if (pred(cand)) { acc = cand; return offset; }
        }
        push(t);
        const int lc = pool[t].ch[0], rc = pool[t].ch[1];
        const int lcnt = cnt(lc);
        const int sz = (int)pool[t].a.size();
        const int mid = offset + lcnt;
        if (mid + sz < r) {
            const int res = min_left_rec(rc, r, mid + sz, pred, acc);
            if (res != mid + sz) return res;
        }
        if (mid < r) {
            const int to = vmin(r, mid + sz) - mid;   // 配列の [0,to) が対象
            if (to == sz) {
                const T cand = acc + pool[t].blk.sum;
                if (pred(cand)) { acc = cand; goto go_left; }
            }
            mat(t);
            for (int j = to - 1; j >= 0; j--) {
                const T cand = acc + pool[t].a[j];
                if (!pred(cand)) return mid + j + 1;
                acc = cand;
            }
        }
        go_left:
        return min_left_rec(lc, r, offset, pred, acc);
    }

    int lower_bound_prefix(const T& x) const {
        if (!(x > T(0))) return 0;
        const int r = max_right(0, [&](const T& s) { return s < x; });
        return std::min(r + 1, size());
    }
    int lower_bound_sorted(int l, int r, const T& x) const {
        assert(0 <= l && l <= r && r <= size());
        int lo = l, hi = r;
        while (lo < hi) { const int mid = lo + (hi - lo) / 2; if (get(mid) < x) lo = mid + 1; else hi = mid; }
        return lo;
    }
    int upper_bound_sorted(int l, int r, const T& x) const {
        assert(0 <= l && l <= r && r <= size());
        int lo = l, hi = r;
        while (lo < hi) { const int mid = lo + (hi - lo) / 2; if (x < get(mid)) hi = mid; else lo = mid + 1; }
        return lo;
    }

    // ==================================================================
    // 便利機能 (O(log N) は不可能。明示的に O(len [log len]) とする)
    // ==================================================================
    void sort(int l, int r) {
        assert(0<=l&&l<=r&&r<=size());
        if (r-l<=1) return;
        vector<T> v = to_vector(l, r);
        std::sort(v.begin(), v.end());
        rebuild_range(l, r, v);
    }
    void stable_sort(int l, int r) {
        assert(0<=l&&l<=r&&r<=size());
        if (r-l<=1) return;
        vector<T> v = to_vector(l, r);
        std::stable_sort(v.begin(), v.end());
        rebuild_range(l, r, v);
    }
    void unique(int l, int r) {
        assert(0<=l&&l<=r&&r<=size());
        vector<T> v = to_vector(l, r);
        v.erase(std::unique(v.begin(), v.end()), v.end());
        rebuild_range(l, r, v);
    }
    template <class Pred>
    void erase_if(int l, int r, Pred pred) {
        assert(0<=l&&l<=r&&r<=size());
        vector<T> v = to_vector(l, r);
        vector<T> w; w.reserve(v.size());
        for (auto& e : v) if (!pred(e)) w.push_back(e);
        rebuild_range(l, r, w);
    }
    template <class Pred>
    void filter(int l, int r, Pred pred) { erase_if(l, r, [&](const T& e){ return !pred(e); }); }
    template <class Pred>
    void partition(int l, int r, Pred pred) {
        assert(0<=l&&l<=r&&r<=size());
        vector<T> v = to_vector(l, r);
        std::partition(v.begin(), v.end(), pred);
        rebuild_range(l, r, v);
    }
    template <class Pred>
    void stable_partition(int l, int r, Pred pred) {
        assert(0<=l&&l<=r&&r<=size());
        vector<T> v = to_vector(l, r);
        std::stable_partition(v.begin(), v.end(), pred);
        rebuild_range(l, r, v);
    }
    template <class F>
    void map(int l, int r, F f) {
        assert(0<=l&&l<=r&&r<=size());
        vector<T> v = to_vector(l, r);
        for (auto& e : v) e = f(e);
        rebuild_range(l, r, v);
    }
    void rebuild_range(int l, int r, const vector<T>& v) {
        auto [a, bc] = split(root, l);
        auto [b, c] = split(bc, r - l);
        free_subtree(b);
        root = merge(merge(a, build(v)), c);
        maybe_rebuild();
    }

    T kth_smallest(int l, int r, int k) const {
        assert(0<=l&&l<r&&r<=size());
        assert(0<=k&&k<r-l);
        vector<T> v = to_vector(l, r);
        std::nth_element(v.begin(), v.begin()+k, v.end());
        return v[k];
    }
    T kth_largest(int l, int r, int k) const {
        assert(0<=l&&l<r&&r<=size());
        assert(0<=k&&k<r-l);
        return kth_smallest(l, r, (r-l-1) - k);
    }
    T median(int l, int r) const {
        assert(0<=l&&l<r&&r<=size());
        return kth_smallest(l, r, (r-l-1)/2);
    }
    long long count_value(int l, int r, const T& x) const {
        assert(0<=l&&l<=r&&r<=size());
        vector<T> v = to_vector(l, r);
        return std::count(v.begin(), v.end(), x);
    }
    template <class Pred>
    long long count_if(int l, int r, Pred pred) const {
        assert(0<=l&&l<=r&&r<=size());
        vector<T> v = to_vector(l, r);
        return std::count_if(v.begin(), v.end(), pred);
    }
    bool is_sorted(int l, int r) const {
        assert(0<=l&&l<=r&&r<=size());
        vector<T> v = to_vector(l, r);
        return std::is_sorted(v.begin(), v.end());
    }
    T mode(int l, int r) const {
        assert(0<=l&&l<r&&r<=size());
        vector<T> v = to_vector(l, r);
        std::sort(v.begin(), v.end());
        T best = v[0]; int bestc=1, cur=1;
        for (size_t i=1;i<v.size();i++){
            if (!(v[i-1]<v[i])) cur++; else cur=1;
            if (cur>bestc){bestc=cur;best=v[i];}
        }
        return best;
    }
    long long mex(int l, int r) const {
        assert(0<=l&&l<=r&&r<=size());
        vector<T> v = to_vector(l, r);
        std::set<T> s(v.begin(), v.end());
        long long m = 0;
        while (s.count((T)m)) m++;
        return m;
    }
    long long inversion(int l, int r) const {
        assert(0<=l&&l<=r&&r<=size());
        vector<T> v = to_vector(l, r);
        vector<T> uniq = v; std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        const int n = (int)v.size();
        vector<int> bit(n+2, 0);
        auto add_bit=[&](int i,int d){ for(i++;i<=n;i+=i&(-i)) bit[i]+=d; };
        auto sum_bit=[&](int i){ int s=0; for(;i>0;i-=i&(-i)) s+=bit[i]; return s; };
        long long inv=0;
        for (int i=n-1;i>=0;i--){
            const int id = (int)(std::lower_bound(uniq.begin(),uniq.end(),v[i]) - uniq.begin());
            inv += sum_bit(id);
            add_bit(id,1);
        }
        return inv;
    }

    // ==================================================================
    // Debug
    // ==================================================================
    // 木の状態を覗くための補助 (性能の確認・デバッグ用)
    static constexpr int chunk_size() { return CS; }      // 1 節点が持つ配列の基準長
    int node_count() const { return nnode(root); }        // 節点数 (≒ size()/CS) O(1)
    int height() const { return height_rec(root); }       // 木の高さ O(節点数)
    static int height_rec(int t) {
        if (t == -1) return 0;
        return 1 + vmax(height_rec(pool[t].ch[0]), height_rec(pool[t].ch[1]));
    }
    void print() const { auto v = to_vector(); for (auto& e : v) cerr << e << ' '; cerr << '\n'; }
    void dump() const { dump_rec(root, 0); }
    static void dump_rec(int t, int depth) {
        if (t == -1) return;
        push(t);
        dump_rec(pool[t].ch[1], depth+1);
        mat(t);
        cerr << string(depth*2, ' ') << "[";
        for (size_t i = 0; i < pool[t].a.size(); i++) cerr << (i ? "," : "") << pool[t].a[i];
        cerr << "] (cnt=" << pool[t].cnt << ")\n";
        dump_rec(pool[t].ch[0], depth+1);
    }
    // 構造検証: 部分木サイズ / 節点数 / heap 性質 / チャンク長 / 集約の一致 を確認する。
    // 内部で push/mat するため lazy は確定するが、値は変わらない。
    bool validate() const { return check_rec(root) >= 0; }
    bool check() const { return validate(); }
    static int check_rec(int t) {
        if (t == -1) return 0;
        push(t);
        const int lc = pool[t].ch[0], rc = pool[t].ch[1];
        if (lc != -1) assert(pool[t].pri >= pool[lc].pri && "heap property violated");
        if (rc != -1) assert(pool[t].pri >= pool[rc].pri && "heap property violated");
        assert(!pool[t].a.empty() && "空のチャンクが残っています");
        assert((int)pool[t].a.size() <= 2 * CS && "チャンクが 2*CS を超えています");
        const int lsz = check_rec(lc);
        const int rsz = check_rec(rc);
        const int expect = lsz + (int)pool[t].a.size() + rsz;
        assert(pool[t].cnt == expect && "subtree size mismatch");
        assert(pool[t].nn == 1 + nnode(lc) + nnode(rc) && "node count mismatch");

        // 集約を作り直して一致を確認する (壊れている集約は bad ビットが立っているので比較しない)
        const Agg snap_blk = pool[t].blk;
        const Agg snap_sub = pool[t].sub;
        mat(t);
        pull_blk(t);
        pull(t);
        cmp_agg(pool[t].blk, snap_blk, "blk");
        cmp_agg(pool[t].sub, snap_sub, "sub");
        return expect;
    }
    static void cmp_agg(const Agg& a, const Agg& b, const char* what) {
        (void)what;
        assert(a.has == b.has);
        if (!a.has) return;
        assert(!(a.sum < b.sum) && !(b.sum < a.sum) && "sum mismatch");
        assert(!(a.mn < b.mn) && !(b.mn < a.mn) && "min mismatch");
        assert(!(a.mx < b.mx) && !(b.mx < a.mx) && "max mismatch");
        if constexpr (HAS_GCD) { if (!(b.bad & BAD_GCD)) assert(a.gc.gcdv == b.gc.gcdv && "gcd mismatch"); }
        if constexpr (HAS_BIT) {
            if (!(b.bad & BAD_BIT)) {
                assert(a.bit.xorv == b.bit.xorv && "xor mismatch");
                assert(a.bit.andv == b.bit.andv && "and mismatch");
                assert(a.bit.orv == b.bit.orv && "or mismatch");
            }
        }
        if constexpr (HAS_KAD) {
            if (!(b.bad & BAD_KAD)) {
                assert(!(a.kd.pre_max < b.kd.pre_max) && !(b.kd.pre_max < a.kd.pre_max) && "pre_max mismatch");
                assert(!(a.kd.suf_max < b.kd.suf_max) && !(b.kd.suf_max < a.kd.suf_max) && "suf_max mismatch");
                assert(!(a.kd.pre_min < b.kd.pre_min) && !(b.kd.pre_min < a.kd.pre_min) && "pre_min mismatch");
                assert(!(a.kd.suf_min < b.kd.suf_min) && !(b.kd.suf_min < a.kd.suf_min) && "suf_min mismatch");
                assert(!(a.kd.best_max < b.kd.best_max) && !(b.kd.best_max < a.kd.best_max) && "best_max mismatch");
                assert(!(a.kd.best_min < b.kd.best_min) && !(b.kd.best_min < a.kd.best_min) && "best_min mismatch");
            }
        }
    }
};

// ============================================================================
// ここから下はテスト / 動作確認用のコードです。
// AtCoder に提出するときは、この領域と solv() 内の呼び出しを削除してください。
// ----------------------------------------------------------------------------
// 【テストの構成】
//   テスト1〜9 : 機能のランダムテスト (vector を正解として突き合わせる)。
//                旧版から引き継いでいるので「機能が落ちていないこと」の証明になる。
//   テストA〜E : チャンク化 treap 固有のテスト (今回の作り直しで追加)。
//                CS を 1〜16 まで落として、チャンクの分割・消滅・併合や
//                書き戻し方式と split/merge 方式の切り替え境界を総当たりする。
//                A: チャンク境界の総合ストレス (validate() を毎回呼ぶ)
//                B: 木の健全性 (節点数・平均チャンク長・高さが理論どおりか)
//                C: 遅延タグ ct / at の実体化まわり
//                D: 短い並べ替えの分岐境界 (長さ 2*CS の前後)
//                E: チャンクの分割と消滅 (同じ位置に挿し続ける / 消し続ける)
//   validate() は部分木サイズ・節点数・heap 性質・チャンク長・集約の整合性を
//   すべて検査するので、これらのテストは内部不変条件の検査も兼ねている。
// ============================================================================

// ---------------------------------------------------------------------------
// [テスト1] 基本操作の総合ランダムストレス
//   vector を正解として、構造操作・区間更新・区間取得をすべて突き合わせる。
//   旧版のテストに以下を追加した:
//     ・multiply / affine を混ぜる (旧版は add / assign のみ)
//     ・shift_left / shift_right / reverse() (全体) / set_range / insert(vector)
//     ・argmin / argmax / minmax / sum() / min() / max() (引数なし版)
//     ・to_vector(l,r) / front / back / at / lower_bound_sorted
//     ・validate() を毎回呼ぶ (集約値の整合性まで検査する)
// ---------------------------------------------------------------------------
void stress_test(int OPS, int lo, int hi, unsigned seed) {
    mt19937 rng(seed);
    vector<long long> ref(lo);
    for (auto& x : ref) x = (int)(rng() % 20) - 10;
    DynamicSequence<long long> ds(ref);

    long long size_sum = 0;
    for (int op = 0; op < OPS; op++) {
        // 列が空になったり大きくなりすぎたりしないよう長さを [lo,hi] 程度に保つ
        // (放っておくとランダムウォークで長さ 0 に張り付き、テストが空回りする)
        if ((int)ref.size() < lo) {
            vector<long long> v(lo - (int)ref.size() + (int)(rng()%4));
            for (auto& e : v) e = (int)(rng()%20)-10;
            ref.insert(ref.end(), v.begin(), v.end());
            ds.append(v);
        } else if ((int)ref.size() > hi) {
            const int cut = (int)ref.size() - (lo + hi) / 2;
            ref.erase(ref.begin(), ref.begin()+cut);
            ds.erase(0, cut);
        }
        size_sum += (long long)ref.size();

        const int type = rng() % 24;
        const int n = (int)ref.size();
        if (type == 0) { // push_back
            long long x = (int)(rng()%20)-10;
            ref.push_back(x); ds.push_back(x);
        } else if (type == 1 && n>0) { // pop_back
            ref.pop_back(); ds.pop_back();
        } else if (type == 2) { // push_front
            long long x = (int)(rng()%20)-10;
            ref.insert(ref.begin(), x); ds.push_front(x);
        } else if (type == 3 && n>0) { // pop_front
            ref.erase(ref.begin()); ds.pop_front();
        } else if (type == 4) { // insert
            int i = rng()%(n+1); long long x=(int)(rng()%20)-10;
            ref.insert(ref.begin()+i, x); ds.insert(i,x);
        } else if (type == 5 && n>0) { // erase single
            int i = rng()%n;
            ref.erase(ref.begin()+i); ds.erase(i);
        } else if (type == 6 && n>0) { // erase range
            int l = rng()%n, r = l + rng()%(n-l+1); if(l>r) std::swap(l,r); if(r>n) r=n;
            ref.erase(ref.begin()+l, ref.begin()+r); ds.erase(l,r);
        } else if (type == 7 && n>0) { // set
            int i = rng()%n; long long x=(int)(rng()%20)-10;
            ref[i]=x; ds.set(i,x);
        } else if (type == 8 && n>0) { // reverse
            int l = rng()%n, r = l + rng()%(n-l+1);
            std::reverse(ref.begin()+l, ref.begin()+r); ds.reverse(l,r);
        } else if (type == 9 && n>0) { // rotate
            int l = rng()%n, r = l + rng()%(n-l+1);
            if (r>l) {
                long long k = (long long)(rng()%(2*(r-l)+1)) - (r-l);
                long long kk = ((k % (r-l)) + (r-l)) % (r-l);
                // 右にkk回rotate = std::rotateで先頭を (len-kk) にする
                std::rotate(ref.begin()+l, ref.begin()+ (l + ( (r-l) - (int)kk) % (r-l) ), ref.begin()+r);
                ds.rotate(l,r,k);
            }
        } else if (type == 10 && n>0) { // move
            int l = rng()%n, r = l + rng()%(n-l+1);
            int len = r-l;
            int pos = rng()%(n-len+1);
            vector<long long> seg(ref.begin()+l, ref.begin()+r);
            ref.erase(ref.begin()+l, ref.begin()+r);
            ref.insert(ref.begin()+pos, seg.begin(), seg.end());
            ds.move(l,r,pos);
        } else if (type == 11 && n>0) { // swap_range
            int a1 = rng()%n, b1 = a1 + rng()%(n-a1+1);
            if (b1>=n) continue;
            int a2 = b1 + rng()%(n-b1+1), b2v = a2 + (a2<n? rng()%(n-a2+1):0);
            if (a2>n) continue;
            if (b2v>n) b2v=n;
            if (!(a1<=b1 && b1<=a2 && a2<=b2v)) continue;
            vector<long long> s1(ref.begin()+a1, ref.begin()+b1);
            vector<long long> s2(ref.begin()+a2, ref.begin()+b2v);
            vector<long long> nref;
            nref.insert(nref.end(), ref.begin(), ref.begin()+a1);
            nref.insert(nref.end(), s2.begin(), s2.end());
            nref.insert(nref.end(), ref.begin()+b1, ref.begin()+a2);
            nref.insert(nref.end(), s1.begin(), s1.end());
            nref.insert(nref.end(), ref.begin()+b2v, ref.end());
            ref = nref;
            ds.swap_range(a1,b1,a2,b2v);
        } else if (type == 12 && n>0) { // add
            int l = rng()%n, r = l + rng()%(n-l+1);
            if (r>l) { long long x=(int)(rng()%10)-5; for(int i=l;i<r;i++) ref[i]+=x; ds.add(l,r,x); }
        } else if (type == 13 && n>0) { // assign
            int l = rng()%n, r = l + rng()%(n-l+1);
            if (r>l) { long long x=(int)(rng()%10)-5; for(int i=l;i<r;i++) ref[i]=x; ds.assign(l,r,x); }
        } else if (type == 14 && n>0) { // multiply (負数・0 を含む)
            int l = rng()%n, r = l + rng()%(n-l+1);
            if (r>l) { long long x=(int)(rng()%5)-2; for(int i=l;i<r;i++) ref[i]*=x; ds.multiply(l,r,x); }
        } else if (type == 15 && n>0) { // affine
            int l = rng()%n, r = l + rng()%(n-l+1);
            if (r>l) {
                long long a=(int)(rng()%5)-2, b=(int)(rng()%7)-3;
                for(int i=l;i<r;i++) ref[i]=a*ref[i]+b;
                ds.affine(l,r,a,b);
            }
        } else if (type == 16 && n>0) { // shift_left / shift_right
            int l = rng()%n, r = l + rng()%(n-l+1);
            if (r>l) {
                int len = r-l;
                long long k = rng()%(2*len+1);
                long long kk = k % len;
                if (rng()%2) {  // 右シフト
                    std::rotate(ref.begin()+l, ref.begin()+(l+(len-(int)kk)%len), ref.begin()+r);
                    ds.shift_right(l,r,k);
                } else {        // 左シフト
                    std::rotate(ref.begin()+l, ref.begin()+(l+(int)kk), ref.begin()+r);
                    ds.shift_left(l,r,k);
                }
            }
        } else if (type == 17) { // reverse() 全体
            std::reverse(ref.begin(), ref.end()); ds.reverse();
        } else if (type == 18) { // insert(pos, vector) / append / prepend
            int m = rng()%4;
            vector<long long> v(m);
            for (auto& e : v) e = (int)(rng()%20)-10;
            int mode = rng()%3;
            if (mode==0) { int i = rng()%(n+1); ref.insert(ref.begin()+i, v.begin(), v.end()); ds.insert(i, v); }
            else if (mode==1) { ref.insert(ref.end(), v.begin(), v.end()); ds.append(v); }
            else { ref.insert(ref.begin(), v.begin(), v.end()); ds.prepend(v); }
        } else if (type == 19 && n>0) { // set_range (長さが変わる置き換え)
            int l = rng()%n, r = l + rng()%(n-l+1);
            int m = rng()%4;
            vector<long long> v(m);
            for (auto& e : v) e = (int)(rng()%20)-10;
            vector<long long> nref(ref.begin(), ref.begin()+l);
            nref.insert(nref.end(), v.begin(), v.end());
            nref.insert(nref.end(), ref.begin()+r, ref.end());
            ref = nref;
            ds.set_range(l, r, v);
        } else if (type == 20 && n>0) { // extract して別の位置に splice し戻す
            int l = rng()%n, r = l + rng()%(n-l+1);
            DynamicSequence<long long> part = ds.extract(l, r);
            vector<long long> seg(ref.begin()+l, ref.begin()+r);
            ref.erase(ref.begin()+l, ref.begin()+r);
            assert(part.to_vector() == seg);
            assert((int)ref.size() == ds.size());
            int pos = rng()%(ds.size()+1);
            ds.splice(pos, part, 0, part.size());
            ref.insert(ref.begin()+pos, seg.begin(), seg.end());
            assert(part.empty());
        } else if (type == 21) { // split して merge し戻す
            int k = rng()%(n+1);
            auto [L, R] = ds.split(k);
            assert(L.size() == k && R.size() == n-k);
            assert(ds.empty());
            L.merge(R);
            assert(R.empty());
            ds = std::move(L);
        } else if (type == 22) { // clone が独立していること
            DynamicSequence<long long> cp = ds.clone();
            assert(cp.to_vector() == ref);
            if (!cp.empty()) cp.add(0, cp.size(), 1000);   // コピーだけ変更
            assert(ds.to_vector() == ref);
        } else if (type == 23 && n>0) { // sub は元を壊さない
            int l = rng()%n, r = l + rng()%(n-l+1);
            DynamicSequence<long long> s = ds.sub(l, r);
            assert(s.to_vector() == vector<long long>(ref.begin()+l, ref.begin()+r));
            assert(ds.to_vector() == ref);
        }

        // ---- 毎回の整合性チェック ----
        assert((int)ref.size() == ds.size());
        assert(ref.empty() == ds.empty());
        auto v = ds.to_vector();
        if (v != ref) {
            cerr << "MISMATCH at op="<<op<<" type="<<type<<"\n";
            cerr << "ref: "; for(auto e:ref) cerr<<e<<' '; cerr<<"\n";
            cerr << "got: "; for(auto e:v) cerr<<e<<' '; cerr<<"\n";
        }
        assert(v == ref);
        assert(ds.validate());
        if (!ref.empty()) {
            // 引数なし版 (全区間) は O(1) で返るパス
            long long total = 0, gmn = ref[0], gmx = ref[0];
            for (long long e : ref) { total += e; gmn = std::min(gmn, e); gmx = std::max(gmx, e); }
            assert(ds.sum() == total);
            assert(ds.min() == gmn);
            assert(ds.max() == gmx);
            assert(ds.sum(0, ds.size()) == total);
            assert(ds.front() == ref.front());
            assert(ds.back() == ref.back());

            int l = rng()%ref.size(), r = l+rng()%(ref.size()-l+1);
            if (r>l) {
                long long s=0; for(int i=l;i<r;i++) s+=ref[i];
                assert(ds.sum(l,r) == s);
                long long mn=ref[l], mx=ref[l];
                int amn=l, amx=l;
                for(int i=l;i<r;i++){
                    if (ref[i]<mn){mn=ref[i];amn=i;}
                    if (ref[i]>mx){mx=ref[i];amx=i;}
                }
                assert(ds.min(l,r)==mn);
                assert(ds.max(l,r)==mx);
                assert(ds.minmax(l,r)==make_pair(mn,mx));
                // argmin/argmax は「最も左の」位置を返す
                assert(ds.argmin(l,r)==amn);
                assert(ds.argmax(l,r)==amx);
                // 区間の to_vector
                assert(ds.to_vector(l,r) == vector<long long>(ref.begin()+l, ref.begin()+r));
                // ランダムな1点の get / at / operator[]
                int i = l + rng()%(r-l);
                assert(ds.get(i)==ref[i] && ds.at(i)==ref[i] && ds[i]==ref[i]);
            }
            assert(ds.sum(3 % (ds.size()+1), 3 % (ds.size()+1)) == 0);  // 空区間は 0
        } else {
            assert(ds.sum() == 0);
            assert(ds.sum(0,0) == 0);
        }
    }
    assert(ds.to_vector() == ref);
    cerr << "[DynamicSequence] stress_test OK (" << OPS << " ops, 長さ帯 [" << lo << "," << hi
         << "], 平均長 " << (double)size_sum/OPS << ")\n";
}

// ---------------------------------------------------------------------------
// [テスト2] gcd/xor/and/or/Kadane系(非線形aggregate) の検証。
//   これらは add/affine と両立しないので、assign と純粋な構造操作のみを行う。
//   (旧版のテストのまま + validate を追加)
// ---------------------------------------------------------------------------
void stress_test_nonlinear() {
    mt19937 rng(999331);
    int n0 = 1 + rng() % 6;
    vector<long long> ref(n0);
    for (auto& x : ref) x = 1 + (int)(rng() % 12);
    DynamicSequence<long long> ds(ref);

    const int OPS = 400000;
    for (int op = 0; op < OPS; op++) {
        // 長さが 0 に張り付かないように下限を保つ
        while ((int)ref.size() < 4) { long long x = 1+(int)(rng()%12); ref.push_back(x); ds.push_back(x); }
        if ((int)ref.size() > 60) { const int cut = (int)ref.size()-40; ref.erase(ref.begin(), ref.begin()+cut); ds.erase(0, cut); }
        const int type = rng() % 9;
        const int n = (int)ref.size();
        if (type == 0) { long long x=1+(int)(rng()%12); ref.push_back(x); ds.push_back(x); }
        else if (type == 1 && n>0) { ref.pop_back(); ds.pop_back(); }
        else if (type == 2) { int i=rng()%(n+1); long long x=1+(int)(rng()%12); ref.insert(ref.begin()+i,x); ds.insert(i,x); }
        else if (type == 3 && n>0) { int i=rng()%n; ref.erase(ref.begin()+i); ds.erase(i); }
        else if (type == 4 && n>0) { int l=rng()%n,r=l+rng()%(n-l+1); std::reverse(ref.begin()+l,ref.begin()+r); ds.reverse(l,r); }
        else if (type == 5 && n>0) {
            int l=rng()%n,r=l+rng()%(n-l+1);
            if (r>l) { long long k=(long long)(rng()%(2*(r-l)+1))-(r-l); long long kk=((k%(r-l))+(r-l))%(r-l);
                std::rotate(ref.begin()+l, ref.begin()+(l+((r-l)-(int)kk)%(r-l)), ref.begin()+r); ds.rotate(l,r,k); }
        }
        else if (type == 6 && n>0) {
            int l=rng()%n,r=l+rng()%(n-l+1); int len=r-l; int pos=rng()%(n-len+1);
            vector<long long> seg(ref.begin()+l, ref.begin()+r);
            ref.erase(ref.begin()+l, ref.begin()+r); ref.insert(ref.begin()+pos, seg.begin(), seg.end());
            ds.move(l,r,pos);
        }
        else if (type == 7 && n>0) { int l=rng()%n,r=l+rng()%(n-l+1); if(r>l){ long long x=1+(int)(rng()%12); for(int i=l;i<r;i++) ref[i]=x; ds.assign(l,r,x);} }
        else if (type == 8 && n>0) { // swap_range
            int a1=rng()%n,b1=a1+rng()%(n-a1+1); if(b1>=n) continue;
            int a2=b1+rng()%(n-b1+1); if(a2>n) continue; int b2v=a2+rng()%(n-a2+1); if(b2v>n) continue;
            if(!(a1<=b1&&b1<=a2&&a2<=b2v)) continue;
            vector<long long> s1(ref.begin()+a1, ref.begin()+b1), s2(ref.begin()+a2, ref.begin()+b2v);
            vector<long long> nref(ref.begin(), ref.begin()+a1);
            nref.insert(nref.end(), s2.begin(), s2.end());
            nref.insert(nref.end(), ref.begin()+b1, ref.begin()+a2);
            nref.insert(nref.end(), s1.begin(), s1.end());
            nref.insert(nref.end(), ref.begin()+b2v, ref.end());
            ref = nref; ds.swap_range(a1,b1,a2,b2v);
        }

        assert((int)ref.size() == ds.size());
        assert(ds.to_vector() == ref);
        assert(ds.validate());
        if (!ref.empty()) {
            int l = rng()%ref.size(), r = l + 1 + rng()%(ref.size()-l);
            long long g=0,x=0,a=~0LL,o=0;
            for (int i=l;i<r;i++){ g=std::gcd(g,ref[i]); x^=ref[i]; a&=ref[i]; o|=ref[i]; }
            assert(ds.gcd_range(l,r)==g);
            assert(ds.xor_sum(l,r)==x);
            assert(ds.and_sum(l,r)==a);
            assert(ds.or_sum(l,r)==o);
            assert(ds.xor_sum(l,l)==0);   // 空区間は 0
            long long best=-1e18, cur=0, bestmn=1e18, curmn=0;
            long long pre=-1e18, prem=1e18;
            long long suf=-1e18, sufm=1e18;
            for (int i=l;i<r;i++){ cur=std::max((long long)ref[i],cur+ref[i]); best=std::max(best,cur);
                                    curmn=std::min((long long)ref[i],curmn+ref[i]); bestmn=std::min(bestmn,curmn); }
            long long s=0;
            for (int i=l;i<r;i++){ s+=ref[i]; pre=std::max(pre,s); prem=std::min(prem,s); }
            s=0;
            for (int i=r-1;i>=l;i--){ s+=ref[i]; suf=std::max(suf,s); sufm=std::min(sufm,s); }
            assert(ds.max_subarray_sum(l,r)==best);
            assert(ds.min_subarray_sum(l,r)==bestmn);
            assert(ds.max_prefix(l,r)==pre);
            assert(ds.min_prefix(l,r)==prem);
            assert(ds.max_suffix(l,r)==suf);
            assert(ds.min_suffix(l,r)==sufm);
        }
    }
    cerr << "[DynamicSequence] stress_test_nonlinear OK (" << OPS << " ops)\n";
}

// ---------------------------------------------------------------------------
// [テスト3・新規] multiply (定数倍) と gcd / Kadane系 の両立を検証する。
//   旧版は「multiply したら gcd も Kadane も無効」としていたが、実際には
//   定数倍なら維持できる。その新しい挙動が正しいことを確認する。
//   (xor/and/or は定数倍でも維持できないのでここでは検証しない)
// ---------------------------------------------------------------------------
void stress_test_multiply_agg() {
    mt19937 rng(20260919);
    vector<long long> ref(1 + rng()%5);
    for (auto& x : ref) x = 1 + (int)(rng()%12);
    DynamicSequence<long long> ds(ref);

    const int OPS = 300000;
    for (int op = 0; op < OPS; op++) {
        while ((int)ref.size() < 4) { long long x = 1+(int)(rng()%12); ref.push_back(x); ds.push_back(x); }
        if ((int)ref.size() > 50) { const int cut = (int)ref.size()-35; ref.erase(ref.begin(), ref.begin()+cut); ds.erase(0, cut); }
        const int n = (int)ref.size();
        const int type = rng() % 7;
        if (type == 0) { long long x = 1+(int)(rng()%12); ref.push_back(x); ds.push_back(x); }
        else if (type == 1 && n>0) { ref.pop_back(); ds.pop_back(); }
        else if (type == 2) { int i=rng()%(n+1); long long x=1+(int)(rng()%12); ref.insert(ref.begin()+i,x); ds.insert(i,x); }
        else if (type == 3 && n>0) { int i=rng()%n; ref.erase(ref.begin()+i); ds.erase(i); }
        else if (type == 4 && n>0) { int l=rng()%n,r=l+rng()%(n-l+1); std::reverse(ref.begin()+l,ref.begin()+r); ds.reverse(l,r); }
        else if (type == 5 && n>0) { // multiply (定数倍)
            int l=rng()%n,r=l+rng()%(n-l+1);
            if (r>l) {
                long long a = (long long)(rng()%5) - 2;      // -2,-1,0,1,2
                for (int i=l;i<r;i++) ref[i]*=a;
                ds.multiply(l,r,a);
            }
        }
        else if (type == 6 && n>0) { // assign (集約を復活させる)
            int l=rng()%n,r=l+rng()%(n-l+1);
            if (r>l) { long long x=1+(int)(rng()%12); for(int i=l;i<r;i++) ref[i]=x; ds.assign(l,r,x); }
        }

        // 値が膨らみすぎないように、絶対値が大きくなったら全体を assign でリセットする
        bool big = false;
        for (long long e : ref) if (e > 1000000000LL || e < -1000000000LL) big = true;
        if (big && !ref.empty()) {
            long long x = 1 + (int)(rng()%12);
            for (auto& e : ref) e = x;
            ds.assign(0, ds.size(), x);
        }

        assert(ds.to_vector() == ref);
        assert(ds.validate());
        if (!ref.empty()) {
            int l = rng()%ref.size(), r = l + 1 + rng()%(ref.size()-l);
            long long g = 0;
            for (int i=l;i<r;i++) g = std::gcd(g, ref[i]);
            assert(ds.gcd_range(l,r) == g);      // multiply 後も gcd が維持されている

            long long best=LLONG_MIN, cur=0, bestmn=LLONG_MAX, curmn=0;
            for (int i=l;i<r;i++){ cur=std::max(ref[i],cur+ref[i]); best=std::max(best,cur);
                                   curmn=std::min(ref[i],curmn+ref[i]); bestmn=std::min(bestmn,curmn); }
            assert(ds.max_subarray_sum(l,r)==best);   // multiply 後も Kadane が維持されている
            assert(ds.min_subarray_sum(l,r)==bestmn);
            long long s=0, pre=LLONG_MIN, prem=LLONG_MAX;
            for (int i=l;i<r;i++){ s+=ref[i]; pre=std::max(pre,s); prem=std::min(prem,s); }
            assert(ds.max_prefix(l,r)==pre);
            assert(ds.min_prefix(l,r)==prem);
            s=0; long long suf=LLONG_MIN, sufm=LLONG_MAX;
            for (int i=r-1;i>=l;i--){ s+=ref[i]; suf=std::max(suf,s); sufm=std::min(sufm,s); }
            assert(ds.max_suffix(l,r)==suf);
            assert(ds.min_suffix(l,r)==sufm);
        }
    }
    cerr << "[DynamicSequence] stress_test_multiply_agg OK (" << OPS << " ops)\n";
}

// ---------------------------------------------------------------------------
// [テスト4] extract / splice / concat / max_right / min_left / 便利機能 の検証
//   (旧版のテスト + 追加分)
// ---------------------------------------------------------------------------
void stress_test_extra() {
    mt19937 rng(555001);

    // extract/splice/concat
    {
        vector<long long> a = {1,2,3,4,5,6,7,8};
        DynamicSequence<long long> ds(a);
        auto sub = ds.extract(2,5); // {3,4,5}
        assert(sub.to_vector() == vector<long long>({3,4,5}));
        assert(ds.to_vector() == vector<long long>({1,2,6,7,8}));
        ds.splice(1, sub, 0, 2); // insert {3,4} at pos1
        assert(ds.to_vector() == vector<long long>({1,3,4,2,6,7,8}));
        assert(sub.to_vector() == vector<long long>({5}));

        vector<DynamicSequence<long long>> parts;
        parts.push_back(DynamicSequence<long long>(vector<long long>{1,2}));
        parts.push_back(DynamicSequence<long long>(vector<long long>{3,4}));
        parts.push_back(DynamicSequence<long long>(vector<long long>{5,6}));
        auto cc = DynamicSequence<long long>::concat(parts);
        assert(cc.to_vector() == vector<long long>({1,2,3,4,5,6}));
    }

    // max_right / min_left / lower_bound_prefix (非負列を仮定)
    {
        vector<long long> a = {2,1,4,0,3,5,2,1};
        DynamicSequence<long long> ds(a);
        for (int l = 0; l <= (int)a.size(); l++) {
            for (long long lim = 0; lim <= 20; lim++) {
                int expect = l;
                long long acc = 0;
                for (int i = l; i < (int)a.size(); i++) {
                    if (acc + a[i] > lim) break;
                    acc += a[i]; expect = i+1;
                }
                int got = ds.max_right(l, [&](long long s){ return s <= lim; });
                assert(got == expect);
            }
        }
        for (int x = 0; x <= 20; x++) {
            int expect;
            if (x <= 0) expect = 0;
            else {
                long long acc=0; expect=(int)a.size();
                for (int i=0;i<(int)a.size();i++){ acc+=a[i]; if (acc>=x){expect=i+1;break;} }
            }
            int got = ds.lower_bound_prefix(x);
            assert(got == expect);
        }
        // min_left
        for (int r = 0; r <= (int)a.size(); r++) {
            for (long long lim = 0; lim <= 20; lim++) {
                int expect = r;
                long long acc = 0;
                for (int i = r-1; i >= 0; i--) {
                    if (acc + a[i] > lim) break;
                    acc += a[i]; expect = i;
                }
                int got = ds.min_left(r, [&](long long s){ return s <= lim; });
                assert(got == expect);
            }
        }
    }

    // 負数を含む affine / multiply (min/maxの入れ替わり確認)
    {
        vector<long long> a = {3,-1,4,-1,5,-9,2,-6};
        DynamicSequence<long long> ds(a);
        vector<long long> ref = a;
        auto check_all = [&]() {
            assert(ds.to_vector() == ref);
            for (int i=0;i<(int)ref.size();i++) {
                for (int j=i+1;j<=(int)ref.size();j++) {
                    long long mn=ref[i], mx=ref[i], s=0;
                    for (int t=i;t<j;t++){ mn=std::min(mn,ref[t]); mx=std::max(mx,ref[t]); s+=ref[t]; }
                    assert(ds.min(i,j)==mn); assert(ds.max(i,j)==mx); assert(ds.sum(i,j)==s);
                }
            }
        };
        check_all();
        for (int i=0;i<(int)ref.size();i++) ref[i] = ref[i]*(-2)+3;
        ds.affine(0, ds.size(), -2, 3);
        check_all();
        for (int i=2;i<6;i++) ref[i] = ref[i]*(-1);
        ds.multiply(2,6,-1);
        check_all();
        // multiply(.., 0) は実質 assign なので gcd も Kadane も有効に戻る
        for (int i=0;i<3;i++) ref[i] = 0;
        ds.multiply(0,3,0);
        check_all();
        assert(ds.gcd_range(0,3) == 0);
    }

    // sort/stable_sort/unique/erase_if/kth_smallest/median/mode/mex/inversion
    {
        for (int iter = 0; iter < 9000; iter++) {
            int n = 1 + rng()%15;
            vector<long long> a(n);
            for (auto& e : a) e = (int)(rng()%10);
            DynamicSequence<long long> ds(a);

            int l = rng()%n, r = l + rng()%(n-l+1);
            if (r<=l) continue;

            // kth_smallest / kth_largest / median
            int k = rng()%(r-l);
            vector<long long> seg(a.begin()+l, a.begin()+r);
            vector<long long> sorted_seg = seg; std::sort(sorted_seg.begin(), sorted_seg.end());
            assert(ds.kth_smallest(l,r,k) == sorted_seg[k]);
            assert(ds.kth_largest(l,r,k) == sorted_seg[(r-l-1)-k]);
            assert(ds.median(l,r) == sorted_seg[(r-l-1)/2]);

            // mex
            {
                std::set<long long> s(seg.begin(), seg.end());
                long long mex=0; while (s.count(mex)) mex++;
                assert(ds.mex(l,r) == mex);
            }
            // inversion
            {
                long long inv=0;
                for (int i=0;i<(int)seg.size();i++) for (int j=i+1;j<(int)seg.size();j++) if (seg[i]>seg[j]) inv++;
                assert(ds.inversion(l,r) == inv);
            }
            // count_value / count_if / is_sorted
            {
                long long x = rng()%10;
                assert(ds.count_value(l,r,x) == std::count(seg.begin(), seg.end(), x));
                assert(ds.count_if(l,r,[](long long e){ return e%2==0; })
                       == std::count_if(seg.begin(), seg.end(), [](long long e){ return e%2==0; }));
                assert(ds.is_sorted(l,r) == std::is_sorted(seg.begin(), seg.end()));
            }
            // mode (同率なら値が小さい方)
            {
                std::map<long long,int> cnt;
                for (long long e : seg) cnt[e]++;
                int bc = 0; long long bv = 0;
                for (auto& [val,c] : cnt) if (c > bc) { bc = c; bv = val; }
                assert(ds.mode(l,r) == bv);
            }
            // sort(l,r) then compare
            {
                DynamicSequence<long long> ds2(a);
                ds2.sort(l,r);
                vector<long long> expect = a;
                std::sort(expect.begin()+l, expect.begin()+r);
                assert(ds2.to_vector() == expect);
                assert(ds2.is_sorted(l,r));
                assert(ds2.validate());
                // ソート済みなので lower_bound_sorted / upper_bound_sorted が使える
                for (long long x = -1; x <= 10; x++) {
                    assert(ds2.lower_bound_sorted(l,r,x)
                           == l + (int)(std::lower_bound(expect.begin()+l, expect.begin()+r, x) - (expect.begin()+l)));
                    assert(ds2.upper_bound_sorted(l,r,x)
                           == l + (int)(std::upper_bound(expect.begin()+l, expect.begin()+r, x) - (expect.begin()+l)));
                }
            }
            // stable_sort / unique / partition / stable_partition / map
            {
                DynamicSequence<long long> d2(a);
                d2.stable_sort(l,r);
                vector<long long> e2 = a; std::stable_sort(e2.begin()+l, e2.begin()+r);
                assert(d2.to_vector() == e2);

                d2.unique(l, r);   // 連続する重複を削除 (ソート済みなので全重複が消える)
                vector<long long> e3(e2.begin()+l, e2.begin()+r);
                e3.erase(std::unique(e3.begin(), e3.end()), e3.end());
                vector<long long> expect3(e2.begin(), e2.begin()+l);
                expect3.insert(expect3.end(), e3.begin(), e3.end());
                expect3.insert(expect3.end(), e2.begin()+r, e2.end());
                assert(d2.to_vector() == expect3);
                assert(d2.validate());

                DynamicSequence<long long> d3(a);
                auto pred = [](long long e){ return e < 5; };
                d3.partition(l,r,pred);
                vector<long long> e4 = a; std::partition(e4.begin()+l, e4.begin()+r, pred);
                assert(d3.to_vector() == e4);

                DynamicSequence<long long> d4(a);
                d4.stable_partition(l,r,pred);
                vector<long long> e5 = a; std::stable_partition(e5.begin()+l, e5.begin()+r, pred);
                assert(d4.to_vector() == e5);

                DynamicSequence<long long> d5(a);
                d5.map(l,r,[](long long e){ return e*3+1; });
                vector<long long> e6 = a; for (int i=l;i<r;i++) e6[i] = e6[i]*3+1;
                assert(d5.to_vector() == e6);
            }
            // erase_if(l,r, even) / filter
            {
                DynamicSequence<long long> ds3(a);
                ds3.erase_if(l,r,[](long long v){ return v%2==0; });
                vector<long long> expect;
                expect.insert(expect.end(), a.begin(), a.begin()+l);
                for (int i=l;i<r;i++) if (a[i]%2!=0) expect.push_back(a[i]);
                expect.insert(expect.end(), a.begin()+r, a.end());
                assert(ds3.to_vector() == expect);

                DynamicSequence<long long> ds4(a);
                ds4.filter(l,r,[](long long v){ return v%2!=0; });   // 奇数だけ残す = 同じ結果
                assert(ds4.to_vector() == expect);
            }
        }
    }
    cerr << "[DynamicSequence] stress_test_extra OK\n";
}

// ---------------------------------------------------------------------------
// [テスト5・新規] max_right / min_left のランダムストレス
//   ・l/r を動かし、しきい値もランダムに振る
//   ・O(log N) 実装の min_left を愚直計算と突き合わせる (旧版は O(log^2 N) 実装で
//     「sum(l,r) を使った二分探索」だったため、降下型の実装は新規に検証が必要)
// ---------------------------------------------------------------------------
void stress_test_binary_search() {
    mt19937 rng(31337);
    for (int iter = 0; iter < 15000; iter++) {
        const int n = 1 + rng() % 40;
        vector<long long> a(n);
        for (auto& e : a) e = rng() % 12;            // 非負 -> 累積和が単調
        DynamicSequence<long long> ds(a);

        // max_right : sum(l,x) <= lim を満たす最大の x
        for (int trial = 0; trial < 6; trial++) {
            const int l = rng() % (n + 1);
            const long long lim = rng() % (12LL * n + 5);
            int expect = l; long long acc = 0;
            for (int i = l; i < n; i++) { if (acc + a[i] > lim) break; acc += a[i]; expect = i + 1; }
            assert(ds.max_right(l, [&](long long s){ return s <= lim; }) == expect);
        }
        // min_left : sum(x,r) <= lim を満たす最小の x
        for (int trial = 0; trial < 6; trial++) {
            const int r = rng() % (n + 1);
            const long long lim = rng() % (12LL * n + 5);
            int expect = r; long long acc = 0;
            for (int i = r - 1; i >= 0; i--) { if (acc + a[i] > lim) break; acc += a[i]; expect = i; }
            assert(ds.min_left(r, [&](long long s){ return s <= lim; }) == expect);
        }
        // lower_bound_prefix
        for (long long x = -2; x <= 12LL * n + 3; x++) {
            int expect;
            if (x <= 0) expect = 0;
            else { long long acc = 0; expect = n; for (int i=0;i<n;i++){ acc+=a[i]; if (acc>=x){expect=i+1;break;} } }
            assert(ds.lower_bound_prefix(x) == expect);
        }
        // lazy (add/assign) を掛けた後でも二分探索が動くこと
        if (n >= 2) {
            const int l = rng()%n, r = l + 1 + rng()%(n-l);
            const long long x = rng()%5;
            for (int i=l;i<r;i++) a[i] += x;
            ds.add(l,r,x);
            const int rr = rng() % (n + 1);
            const long long lim = rng() % (20LL * n + 5);
            int expect = rr; long long acc = 0;
            for (int i = rr - 1; i >= 0; i--) { if (acc + a[i] > lim) break; acc += a[i]; expect = i; }
            assert(ds.min_left(rr, [&](long long s){ return s <= lim; }) == expect);
            assert(ds.validate());
        }
    }
    cerr << "[DynamicSequence] stress_test_binary_search OK\n";
}

// ---------------------------------------------------------------------------
// [テスト6・新規] AGG (集約の取捨選択) を変えても共通機能が同じ結果になること
//   sum/min/max/構造操作/lazy だけを使って、どの AGG でも一致することを確認する。
// ---------------------------------------------------------------------------
template <int AGG>
void stress_test_agg(const string& name, int ops) {
    mt19937 rng(24680);
    vector<long long> ref(1 + rng()%5);
    for (auto& x : ref) x = (int)(rng()%20)-10;
    DynamicSequence<long long, AGG> ds(ref);

    for (int op = 0; op < ops; op++) {
        while ((int)ref.size() < 4) { long long x = (int)(rng()%20)-10; ref.push_back(x); ds.push_back(x); }
        if ((int)ref.size() > 60) { const int cut = (int)ref.size()-40; ref.erase(ref.begin(), ref.begin()+cut); ds.erase(0, cut); }
        const int n = (int)ref.size();
        const int type = rng() % 9;
        if (type == 0) { long long x=(int)(rng()%20)-10; ref.push_back(x); ds.push_back(x); }
        else if (type == 1 && n>0) { ref.pop_back(); ds.pop_back(); }
        else if (type == 2) { int i=rng()%(n+1); long long x=(int)(rng()%20)-10; ref.insert(ref.begin()+i,x); ds.insert(i,x); }
        else if (type == 3 && n>0) { int i=rng()%n; ref.erase(ref.begin()+i); ds.erase(i); }
        else if (type == 4 && n>0) { int l=rng()%n,r=l+rng()%(n-l+1); std::reverse(ref.begin()+l,ref.begin()+r); ds.reverse(l,r); }
        else if (type == 5 && n>0) { int l=rng()%n,r=l+rng()%(n-l+1); if(r>l){ long long x=(int)(rng()%10)-5; for(int i=l;i<r;i++) ref[i]+=x; ds.add(l,r,x);} }
        else if (type == 6 && n>0) { int l=rng()%n,r=l+rng()%(n-l+1); if(r>l){ long long x=(int)(rng()%10)-5; for(int i=l;i<r;i++) ref[i]=x; ds.assign(l,r,x);} }
        else if (type == 7 && n>0) { int l=rng()%n,r=l+rng()%(n-l+1); if(r>l){ long long a=(int)(rng()%5)-2,b=(int)(rng()%7)-3; for(int i=l;i<r;i++) ref[i]=a*ref[i]+b; ds.affine(l,r,a,b);} }
        else if (type == 8 && n>1) { int l=rng()%n,r=l+1+rng()%(n-l); long long k=rng()%(2*(r-l)+1);
            long long kk=k%(r-l); std::rotate(ref.begin()+l, ref.begin()+(l+((r-l)-(int)kk)%(r-l)), ref.begin()+r); ds.rotate(l,r,k); }

        assert(ds.to_vector() == ref);
        assert(ds.validate());
        if (!ref.empty()) {
            int l = rng()%ref.size(), r = l + 1 + rng()%(ref.size()-l);
            long long s=0, mn=ref[l], mx=ref[l];
            for (int i=l;i<r;i++){ s+=ref[i]; mn=std::min(mn,ref[i]); mx=std::max(mx,ref[i]); }
            assert(ds.sum(l,r)==s);
            assert(ds.min(l,r)==mn);
            assert(ds.max(l,r)==mx);
            assert(ds.minmax(l,r)==make_pair(mn,mx));
            assert(ds.max_right(0,[&](long long v){ return v <= LLONG_MAX/4; }) == ds.size());
        }
    }
    cerr << "[DynamicSequence] stress_test_agg<" << name << "> OK (" << ops << " ops)\n";
}

// AGG_KADANE だけを有効にしたインスタンスで Kadane 系が正しいことを確認する
void stress_test_agg_kadane_only() {
    mt19937 rng(1357);
    vector<long long> ref(1 + rng()%5);
    for (auto& x : ref) x = (int)(rng()%20)-10;
    DynamicSequence<long long, dseq::AGG_KADANE> ds(ref);
    for (int op = 0; op < 200000; op++) {
        while ((int)ref.size() < 4) { long long x = (int)(rng()%20)-10; ref.push_back(x); ds.push_back(x); }
        if ((int)ref.size() > 60) { const int cut = (int)ref.size()-40; ref.erase(ref.begin(), ref.begin()+cut); ds.erase(0, cut); }
        const int n = (int)ref.size();
        const int type = rng() % 5;
        if (type == 0) { long long x=(int)(rng()%20)-10; ref.push_back(x); ds.push_back(x); }
        else if (type == 1 && n>0) { ref.pop_back(); ds.pop_back(); }
        else if (type == 2) { int i=rng()%(n+1); long long x=(int)(rng()%20)-10; ref.insert(ref.begin()+i,x); ds.insert(i,x); }
        else if (type == 3 && n>0) { int i=rng()%n; ref.erase(ref.begin()+i); ds.erase(i); }
        else if (type == 4 && n>0) { int l=rng()%n,r=l+rng()%(n-l+1); std::reverse(ref.begin()+l,ref.begin()+r); ds.reverse(l,r); }
        assert(ds.to_vector() == ref);
        if (!ref.empty()) {
            int l = rng()%ref.size(), r = l + 1 + rng()%(ref.size()-l);
            long long best=LLONG_MIN, cur=0;
            for (int i=l;i<r;i++){ cur=std::max(ref[i],cur+ref[i]); best=std::max(best,cur); }
            assert(ds.max_subarray_sum(l,r)==best);
        }
    }
    cerr << "[DynamicSequence] stress_test_agg_kadane_only OK\n";
}

// ---------------------------------------------------------------------------
// [テスト7・新規] 新規 API と境界条件
// ---------------------------------------------------------------------------
void test_new_api_and_edge() {
    using DS = DynamicSequence<long long>;

    // --- 空列 ---
    {
        DS ds;
        assert(ds.size() == 0 && ds.empty());
        assert(ds.sum() == 0);
        assert(ds.sum(0,0) == 0);
        assert(ds.to_vector().empty());
        assert(ds.validate());
        assert(ds.max_right(0, [](long long){ return true; }) == 0);
        assert(ds.min_left(0, [](long long){ return true; }) == 0);
        assert(ds.lower_bound_prefix(0) == 0);
        assert(ds.lower_bound_prefix(5) == 0);   // 見つからない -> size()==0
        ds.erase(0,0);            // 空区間の削除は何もしない
        ds.reverse(0,0);
        ds.reverse();             // 空でも安全
        ds.add(0,0,5);
        ds.assign(0,0,5);
        ds.rotate(0,0,3);
        assert(ds.validate());
        // 1要素
        ds.push_back(7);
        assert(ds.size()==1 && ds.sum()==7 && ds.min()==7 && ds.max()==7);
        assert(ds.min(0,1)==7 && ds.max(0,1)==7 && ds.sum(0,1)==7);
        assert(ds.argmin(0,1)==0 && ds.argmax(0,1)==0);
        assert(ds.gcd_range(0,1)==7 && ds.xor_sum(0,1)==7 && ds.and_sum(0,1)==7 && ds.or_sum(0,1)==7);
        assert(ds.max_subarray_sum(0,1)==7 && ds.min_subarray_sum(0,1)==7);
        ds.reverse(0,1); ds.rotate(0,1,5);   // 長さ1なら何も起きない
        assert(ds.to_vector() == vector<long long>({7}));
        assert(ds.validate());
    }

    // --- clear / reset / swap / resize / sub / DynamicSequence(n, val) ---
    {
        DS a({1,2,3,4,5});
        DS b({9,8});
        a.swap(b);                        // O(1)
        assert(a.to_vector() == vector<long long>({9,8}));
        assert(b.to_vector() == vector<long long>({1,2,3,4,5}));
        assert(a.validate() && b.validate());

        b.reset({7,7,7});
        assert(b.to_vector() == vector<long long>({7,7,7}));
        b.clear();
        assert(b.empty() && b.validate());

        DS c(4, 3LL);                     // 4 個の 3
        assert(c.to_vector() == vector<long long>({3,3,3,3}));
        c.resize(6, 1LL);
        assert(c.to_vector() == vector<long long>({3,3,3,3,1,1}));
        c.resize(2);
        assert(c.to_vector() == vector<long long>({3,3}));
        c.resize(0);
        assert(c.empty());
        assert(c.validate());

        DS d({1,2,3,4,5,6});
        DS e = d.sub(1,4);                // コピーなので d は変わらない
        assert(e.to_vector() == vector<long long>({2,3,4}));
        assert(d.to_vector() == vector<long long>({1,2,3,4,5,6}));
        e.set(0, 99);
        assert(d.to_vector() == vector<long long>({1,2,3,4,5,6}));
        assert(d.validate() && e.validate());
    }

    // --- insert(vector) / append / prepend / set_range ---
    {
        DS ds({1,2,3});
        ds.insert(1, vector<long long>{10,11});
        assert(ds.to_vector() == vector<long long>({1,10,11,2,3}));
        ds.append({4,5});
        assert(ds.to_vector() == vector<long long>({1,10,11,2,3,4,5}));
        ds.prepend({0});
        assert(ds.to_vector() == vector<long long>({0,1,10,11,2,3,4,5}));
        ds.set_range(1, 4, {100});     // [1,4) を {100} に置き換え (長さが変わる)
        assert(ds.to_vector() == vector<long long>({0,100,2,3,4,5}));
        ds.insert(0, vector<long long>{});   // 空の挿入は何もしない
        assert(ds.to_vector() == vector<long long>({0,100,2,3,4,5}));
        assert(ds.validate());
    }

    // --- clone の独立性 (lazy が残った状態でも正しくコピーされること) ---
    {
        DS ds({1,2,3,4,5});
        ds.add(0, 3, 10);        // lazy を残す
        ds.reverse(1, 4);
        DS cp = ds.clone();
        assert(cp.to_vector() == ds.to_vector());
        cp.add(0, cp.size(), 100);
        assert(cp.to_vector() != ds.to_vector());
        assert(ds.validate() && cp.validate());
    }

    // --- argmin / argmax は「最も左の」位置 ---
    {
        DS ds({3,1,4,1,5,1,9});
        assert(ds.argmin(0,7) == 1);      // 1 は index 1,3,5 にあるが最左
        assert(ds.argmax(0,7) == 6);
        assert(ds.argmin(2,6) == 3);
        assert(ds.argmax(2,6) == 4);
    }

    // --- 例外的な rotate / move / swap_range の境界 ---
    {
        DS ds({1,2,3,4,5});
        ds.rotate(0,5,0);                             // k==0 は何もしない
        assert(ds.to_vector() == vector<long long>({1,2,3,4,5}));
        ds.rotate(0,5,5);                             // 区間長と同じなら元通り
        assert(ds.to_vector() == vector<long long>({1,2,3,4,5}));
        ds.rotate(0,5,-7);                            // 負でも区間長で mod
        assert(ds.to_vector() == vector<long long>({3,4,5,1,2}));
        ds.rotate(0,5,7);
        assert(ds.to_vector() == vector<long long>({1,2,3,4,5}));
        ds.move(0,0,0);                               // 空区間の move は何もしない
        assert(ds.to_vector() == vector<long long>({1,2,3,4,5}));
        ds.swap_range(0,2,2,5);                       // 隣接区間 (長さが違ってもよい)
        assert(ds.to_vector() == vector<long long>({3,4,5,1,2}));
        ds.swap_range(3,5,0,3);                       // 引数の順序が逆でも正規化される
        assert(ds.to_vector() == vector<long long>({1,2,3,4,5}));
        assert(ds.validate());
    }

    // --- reserve_pool (プール事前確保) が動くこと ---
    {
        DS::reserve_pool(1 << 16);
        DS ds(vector<long long>(1000, 1LL));
        assert(ds.sum() == 1000);
        assert(ds.validate());
    }

    cerr << "[DynamicSequence] test_new_api_and_edge OK\n";
}

// ---------------------------------------------------------------------------
// [テスト8・新規] 大規模テスト兼 簡易ベンチマーク
//   ・AGG を絞るとどれだけ速くなるかを実測して表示する
// ---------------------------------------------------------------------------
template <int AGG>
void bench_one(const string& name, int N) {
    mt19937_64 rng(777);
    auto MS = [](chrono::steady_clock::time_point s, chrono::steady_clock::time_point e){
        return (long long)chrono::duration_cast<chrono::milliseconds>(e - s).count(); };

    vector<long long> a(N);
    for (auto& x : a) x = (long long)(rng() % 1000000) - 500000;

    auto t0 = chrono::steady_clock::now();
    DynamicSequence<long long, AGG> ds(a);
    auto t1 = chrono::steady_clock::now();
    for (int i = 0; i < N; i++) ds.insert(rng() % (ds.size() + 1), (long long)(rng() % 1000));
    auto t2 = chrono::steady_clock::now();
    long long sink = 0;
    for (int i = 0; i < N; i++) { int l = rng()%ds.size(), r = l+1+rng()%(ds.size()-l); sink += ds.sum(l,r); }
    auto t3 = chrono::steady_clock::now();
    for (int i = 0; i < N; i++) { int l = rng()%ds.size(), r = l+1+rng()%(ds.size()-l); ds.reverse(l,r); }
    auto t4 = chrono::steady_clock::now();
    for (int i = 0; i < N; i++) { int r = rng()%(ds.size()+1); sink += ds.min_left(r, [&](long long s){ return s < 1000; }); }
    auto t5 = chrono::steady_clock::now();

    cerr << "  " << name
         << " : sizeof(Node)=" << sizeof(typename DynamicSequence<long long, AGG>::Node)
         << "  build=" << MS(t0,t1) << "ms"
         << "  insert=" << MS(t1,t2) << "ms"
         << "  sum=" << MS(t2,t3) << "ms"
         << "  reverse=" << MS(t3,t4) << "ms"
         << "  min_left=" << MS(t4,t5) << "ms\n";
    assert(sink != 0);
    assert((int)ds.size() == 2*N);
}

void large_scale_test() {
    const int N = 200000;
    cerr << "[DynamicSequence] 大規模テスト / ベンチマーク (N=" << N << ", 各操作 " << N << " 回)\n";
    bench_one<dseq::AGG_ALL>("AGG_ALL   ", N);
    bench_one<dseq::AGG_KADANE>("AGG_KADANE", N);
    bench_one<dseq::AGG_BASIC>("AGG_BASIC ", N);

    // 大きい N でも壊れないことの確認 (validate は O(N) なので最後に1回だけ)
    vector<long long> a(N);
    mt19937_64 rng(1);
    for (auto& x : a) x = (long long)(rng() % 1000) - 500;
    DynamicSequence<long long, dseq::AGG_BASIC> ds(a);
    for (int i = 0; i < 100000; i++) {
        const int l = rng()%ds.size(), r = l+1+rng()%(ds.size()-l);
        const int t = rng()%4;
        if (t==0) ds.reverse(l,r);
        else if (t==1) ds.rotate(l,r,(long long)(rng()%1000)-500);
        else if (t==2) ds.add(l,r,1);
        else ds.assign(l,r,3);
    }
    assert(ds.validate());
    long long s = 0;
    for (long long x : ds.to_vector()) s += x;
    assert(ds.sum() == s);
    cerr << "[DynamicSequence] large_scale_test OK\n";
}

// ---------------------------------------------------------------------------
// [動作確認] 使い方の説明を兼ねたデモ
// ---------------------------------------------------------------------------
void test_DynamicSequence() {
    using DS = DynamicSequence<long long>;

    cerr << "========== DynamicSequence 動作確認 ==========\n";

    // ================================================================
    // 集約の取捨選択 (AGG) について
    // ================================================================
    {
        cerr << "\n[集約の取捨選択 (AGG テンプレート引数)]\n";

        // DynamicSequence<T> は既定で「全部入り」。
        // 使わない集約を外すと pull() の計算とノードのメモリが減り、
        // insert/erase/reverse などの構造操作が数倍速くなる。
        //
        //   DynamicSequence<long long>                   ... 全部入り
        //   DynamicSequence<long long, dseq::AGG_BASIC>  ... sum/min/max のみ
        //   DynamicSequence<long long, dseq::AGG_KADANE> ... sum/min/max + Kadane系
        //   DynamicSequence<long long, dseq::AGG_GCD|dseq::AGG_BITWISE> ... など
        //
        // 外した集約のクエリを呼ぶとコンパイルエラーになるので事故らない。
        cerr << "sizeof(Node) 全部入り  : " << sizeof(DS::Node) << " byte\n";
        cerr << "sizeof(Node) AGG_BASIC : "
             << sizeof(DynamicSequence<long long, dseq::AGG_BASIC>::Node) << " byte\n";

        DynamicSequence<long long, dseq::AGG_BASIC> fast({1,2,3,4,5});
        cerr << "AGG_BASIC でも sum/min/max は使える: sum(0,5)=" << fast.sum(0,5)
             << " min(0,5)=" << fast.min(0,5) << " max(0,5)=" << fast.max(0,5) << '\n';
    }

    // ================================================================
    // 基本操作
    // ================================================================
    {
        cerr << "\n[基本操作]\n";

        DS ds({1, 2, 3, 4, 5});

        cerr << "初期列: ";
        ds.print();

        cerr << "size() = " << ds.size() << '\n';
        cerr << "empty() = " << (ds.empty() ? "true" : "false") << '\n';
        cerr << "ds[2] = " << ds[2] << " (at(2) も同じ: " << ds.at(2) << ")\n";
        cerr << "front() = " << ds.front() << '\n';
        cerr << "back() = " << ds.back() << '\n';

        ds.set(2, 10);
        cerr << "set(2, 10): ";
        ds.print();

        // i == size() なら末尾への挿入
        ds.insert(2, 7);
        cerr << "insert(2, 7): ";
        ds.print();

        ds.erase(3);
        cerr << "erase(3): ";
        ds.print();

        // [l,r) の区間削除
        ds.erase(1, 3);
        cerr << "erase(1, 3): ";
        ds.print();

        ds.push_front(0);
        ds.push_back(6);
        cerr << "push_front(0), push_back(6): ";
        ds.print();

        ds.pop_front();
        ds.pop_back();
        cerr << "pop_front(), pop_back(): ";
        ds.print();

        // 【新規】vector をまとめて挿入 (1要素ずつより速い)
        ds.insert(1, vector<long long>{100, 200});
        cerr << "insert(1, {100,200}): ";
        ds.print();

        ds.append({7, 8});
        ds.prepend({-1});
        cerr << "append({7,8}), prepend({-1}): ";
        ds.print();

        // 【新規】長さを変える
        ds.resize(4);
        cerr << "resize(4): ";
        ds.print();
        ds.resize(6, 9);
        cerr << "resize(6, 9): ";
        ds.print();
    }

    // ================================================================
    // split / merge / concat / clone / swap
    // ================================================================
    {
        cerr << "\n[split / merge / concat / clone / swap]\n";

        DS ds({1, 2, 3, 4, 5});

        // split(k):
        // 先頭k個 / 残り に分割する。
        // 元のds自身は空になる。
        auto [L, R] = ds.split(2);

        cerr << "split(2):\n";
        cerr << "左: "; L.print();
        cerr << "右: "; R.print();

        // Lの末尾にRを連結。
        // Rは空になる。
        L.merge(R);

        cerr << "左.merge(右): ";
        L.print();

        DS A({10, 20});
        DS B({30, 40});
        DS C({50, 60});

        vector<DS> seqs;
        seqs.push_back(std::move(A));
        seqs.push_back(std::move(B));
        seqs.push_back(std::move(C));

        // seqs内のSequenceを順番に連結。
        // 各要素はmergeによって空になる。
        DS all = DS::concat(seqs);

        cerr << "concat: ";
        all.print();

        // cloneは深いコピー。
        // 元とコピーを独立して変更できる。
        DS copied = all.clone();
        copied.set(0, 999);

        cerr << "clone後の元: ";
        all.print();

        cerr << "cloneを変更: ";
        copied.print();

        // 【新規】swap は O(1) (根の index を交換するだけ)
        all.swap(copied);
        cerr << "swap後の all: ";
        all.print();
        cerr << "swap後の copied: ";
        copied.print();
    }

    // ================================================================
    // 構造操作
    // ================================================================
    {
        cerr << "\n[reverse / rotate / shift / move / swap_range / sub]\n";

        DS ds({1, 2, 3, 4, 5, 6, 7});

        // [l,r) を反転。
        ds.reverse(1, 6);
        cerr << "reverse(1, 6): ";
        ds.print();

        // 【新規】引数なしの reverse() は全体反転 (O(1)。根にフラグを立てるだけ)
        ds.reverse();
        cerr << "reverse() 全体: ";
        ds.print();
        ds.reverse();

        // [l,r) を右にk個回転。
        // kは負でもよく、負なら左回転。
        ds.rotate(1, 6, 2);
        cerr << "rotate(1, 6, +2): ";
        ds.print();

        ds.rotate(1, 6, -2);
        cerr << "rotate(1, 6, -2): ";
        ds.print();

        // shift_right / shift_left はrotateの別名的なAPI。
        ds.shift_right(0, 3, 1);
        cerr << "shift_right(0, 3, 1): ";
        ds.print();

        ds.shift_left(0, 3, 1);
        cerr << "shift_left(0, 3, 1): ";
        ds.print();

        // move(l,r,pos):
        // [l,r) を一度削除した「後」の列に対して、
        // pos番目の位置へ挿入する。
        //
        // つまりposは「元の列での位置」ではない。
        // 移動対象を除いた後の列での位置。
        ds.move(2, 5, ds.size() - 3);

        cerr << "move(2, 5, 取り出した区間を末尾へ): ";
        ds.print();

        // swap_range:
        // [l1,r1) と [l2,r2) を交換。
        // 2区間は重複禁止。
        // 長さが異なっていてもよい。
        // 隣接区間はOK。
        ds.swap_range(0, 2, 3, 5);

        cerr << "swap_range(0, 2, 3, 5): ";
        ds.print();

        // 【新規】sub(l,r) は「元を壊さずに [l,r) のコピーを取る」
        DS part = ds.sub(1, 4);
        cerr << "sub(1,4) の中身: ";
        part.print();
        cerr << "sub しても元は無傷: ";
        ds.print();
    }

    // ================================================================
    // extract / splice / set_range
    // ================================================================
    {
        cerr << "\n[extract / splice / set_range]\n";

        DS ds({1, 2, 3, 4, 5, 6});

        // [l,r) を元Sequenceから切り離して、
        // 新しいDynamicSequenceとして返す。
        DS part = ds.extract(2, 5);

        cerr << "extract(2, 5):\n";
        cerr << "元: ";
        ds.print();
        cerr << "取り出した区間: ";
        part.print();

        DS other({10, 20, 30});

        // splice(pos, other, l, r):
        // otherの[l,r)を切り取り、
        // 自分のposへ挿入する。
        //
        // 注意:
        // ・other側の区間は削除される
        // ・posは「自分側」の挿入位置
        ds.splice(1, other, 1, 3);

        cerr << "splice(pos=1, other[1,3)):\n";
        cerr << "自分: ";
        ds.print();
        cerr << "other: ";
        other.print();

        // 【新規】set_range(l, r, v) は [l,r) を列 v で置き換える。
        // 長さが変わってもよい。
        ds.set_range(0, 2, {77, 88, 99});
        cerr << "set_range(0, 2, {77,88,99}): ";
        ds.print();
    }

    // ================================================================
    // Aggregate
    // ================================================================
    {
        cerr << "\n[Aggregate]\n";

        DS ds({3, -2, 5, -1, 4, -6, 2});

        cerr << "列: ";
        ds.print();

        // 基本Aggregate。
        // 任意の値に対して使用可能。
        // 引数なし版 (全区間) は O(1) で返る。
        cerr << "sum() = " << ds.sum() << " / min() = " << ds.min() << " / max() = " << ds.max() << '\n';

        cerr << "sum(1, 6) = " << ds.sum(1, 6) << '\n';
        cerr << "min(1, 6) = " << ds.min(1, 6) << '\n';
        cerr << "max(1, 6) = " << ds.max(1, 6) << '\n';

        // 【新規】minmax は1回の走査で両方返す
        auto [mn, mx] = ds.minmax(1, 6);
        cerr << "minmax(1, 6) = (" << mn << ", " << mx << ")\n";

        // 【新規】argmin / argmax は「最小/最大をとる最も左の位置」を O(log N) で返す
        cerr << "argmin(0, 7) = " << ds.argmin(0, 7) << " (値 " << ds[ds.argmin(0,7)] << ")\n";
        cerr << "argmax(0, 7) = " << ds.argmax(0, 7) << " (値 " << ds[ds.argmax(0,7)] << ")\n";

        // 以下のprefix/suffix/subarray系 (Kadane系) は
        // 「add / affine をその区間に適用した後」は無効になる。
        //
        // assign と multiply なら維持できるので使用可能。
        // (multiply は「部分和が一律 a 倍される」だけなので維持できる)

        cerr << "max_prefix(0, 7) = " << ds.max_prefix(0, 7) << '\n';
        cerr << "max_suffix(0, 7) = " << ds.max_suffix(0, 7) << '\n';
        cerr << "min_prefix(0, 7) = " << ds.min_prefix(0, 7) << '\n';
        cerr << "min_suffix(0, 7) = " << ds.min_suffix(0, 7) << '\n';
        cerr << "max_subarray_sum(0, 7) = " << ds.max_subarray_sum(0, 7) << '\n';
        cerr << "min_subarray_sum(0, 7) = " << ds.min_subarray_sum(0, 7) << '\n';
    }

    // ================================================================
    // gcd / xor / and / or
    // ================================================================
    {
        cerr << "\n[gcd / xor / and / or]\n";

        DS ds({12, 18, 6, 10});

        cerr << "列: ";
        ds.print();

        // 維持できる/できないの一覧:
        //   assign   : gcd / xor / and / or / Kadane すべて維持できる (再計算可能)
        //   multiply : gcd と Kadane は維持できる。xor / and / or は維持できない
        //   add      : どれも維持できない
        //   affine   : どれも維持できない
        //
        // 維持できない操作をした区間に問い合わせると assert で落ちる
        // (黙って間違った値を返さない)。
        //
        // gcd は負数も絶対値で扱い、結果は常に非負。

        cerr << "gcd_range(0, 4) = " << ds.gcd_range(0, 4) << '\n';
        cerr << "xor_sum(0, 4) = " << ds.xor_sum(0, 4) << '\n';
        cerr << "and_sum(0, 4) = " << ds.and_sum(0, 4) << '\n';
        cerr << "or_sum(0, 4) = " << ds.or_sum(0, 4) << '\n';

        // multiply の後も gcd は使える (|a| 倍になる)
        ds.multiply(0, 4, 3);
        cerr << "multiply(0, 4, 3) 後: ";
        ds.print();
        cerr << "multiply後のgcd_range = " << ds.gcd_range(0, 4) << " (元の gcd の 3 倍)\n";

        // assignなら全ての非線形Aggregateを復活できる。
        ds.assign(1, 3, 8);
        cerr << "assign(1, 3, 8) 後: ";
        ds.print();
        cerr << "assign後のgcd_range = " << ds.gcd_range(0, 4) << '\n';
        cerr << "assign後のxor_sum  = " << ds.xor_sum(0, 4) << '\n';
    }

    // ================================================================
    // Range Lazy
    // ================================================================
    {
        cerr << "\n[Range Lazy]\n";

        DS ds({1, 2, 3, 4, 5});

        // add:
        // x -> x + b
        //
        // sum/min/max は維持可能。
        // gcd/xor/and/or/Kadane系は無効になる。
        ds.add(1, 4, 10);

        cerr << "add(1, 4, 10): ";
        ds.print();

        // assign:
        // 区間全体をxにする。
        //
        // 全ての集約を再計算できるので、
        // assign後はgcd/xor/and/or/Kadane系も使用可能。
        ds.assign(0, 2, 7);

        cerr << "assign(0, 2, 7): ";
        ds.print();

        // multiply:
        // x -> ax
        //
        // aが0なら実質assign (集約もすべて復活する)。
        // gcdとKadane系は維持される。xor/and/orは無効になる。
        ds.multiply(1, 4, 2);

        cerr << "multiply(1, 4, 2): ";
        ds.print();

        // affine:
        // x -> ax + b
        //
        // sum/min/maxは更新可能。
        // b != 0 なので gcd/xor/and/or/Kadane系は無効になる。
        //
        // a<0 の場合はmin/maxやprefix/suffix等の
        // 最大・最小関係が入れ替わる (実装済み)。
        ds.affine(0, 5, 2, 1);

        cerr << "affine(0, 5, x -> 2x+1): ";
        ds.print();
        cerr << "affine後も sum/min/max は使える: sum=" << ds.sum()
             << " min=" << ds.min() << " max=" << ds.max() << '\n';
    }

    // ================================================================
    // 二分探索
    // ================================================================
    {
        cerr << "\n[二分探索]\n";

        // max_right / min_left は「累積sumの単調性」が必要。
        //
        // 例えば全要素が非負なら、
        // sum(l,r) は r を増やすほど単調非減少になるので、
        // 累積和による境界探索が可能。
        DS ds({2, 3, 1, 4, 2});

        cerr << "列: ";
        ds.print();

        // max_right(l, pred):
        //
        // pred(sum(l,x)) が
        // true -> false
        // と1回だけ変化することを仮定。
        //
        // pred(0) == true が必須。
        int r = ds.max_right(0, [&](long long s) { return s < 6; });
        cerr << "max_right(0, sum < 6) = " << r << '\n';

        // min_left(r, pred):
        //
        // pred(sum(x,r)) が
        // xを左へ広げるにつれて
        // true -> false と境界を持つことを仮定。
        //
        // pred(0) == true が必須。
        //
        // 【高速化】旧版は sum(l,r) を使った二分探索で O(log^2 N) だったが、
        // 現在は max_right と対称な降下型実装で O(log N)。
        int l = ds.min_left(5, [&](long long s) { return s < 8; });
        cerr << "min_left(5, sum < 8) = " << l << '\n';

        // lower_bound_prefix(x):
        //
        // prefix sum が単調非減少であることを仮定。
        // つまり基本的には「全要素が非負」。
        //
        // prefix sum >= x となる最初の位置を返す。
        // 見つからなければsize()。
        cerr << "lower_bound_prefix(7) = " << ds.lower_bound_prefix(7) << '\n';

        // 【新規】lower_bound_sorted / upper_bound_sorted:
        // 区間 [l,r) がソート済みである前提で、値による二分探索を行う。
        // get(i) が O(log N) なので O(log^2 N)。
        DS sorted_ds({1, 3, 3, 5, 7, 9});
        cerr << "lower_bound_sorted(0,6,3) = " << sorted_ds.lower_bound_sorted(0, 6, 3) << '\n';
        cerr << "upper_bound_sorted(0,6,3) = " << sorted_ds.upper_bound_sorted(0, 6, 3) << '\n';
        cerr << "lower_bound_sorted(0,6,4) = " << sorted_ds.lower_bound_sorted(0, 6, 4) << '\n';
    }

    // ================================================================
    // sort / stable_sort / unique / erase_if / filter
    // partition / stable_partition / map
    // ================================================================
    {
        cerr << "\n[便利機能: 並べ替え・削除・変換]\n";

        DS ds({5, 2, 4, 2, 3, 1, 2});

        // [l,r)をソート。
        // O(len log len)
        ds.sort(0, ds.size());

        cerr << "sort: ";
        ds.print();

        // 【新規】is_sorted で確認できる (O(len))
        cerr << "is_sorted(0, size()) = " << (ds.is_sorted(0, ds.size()) ? "true" : "false") << '\n';

        // stable_sort。
        // 値だけを保持するこのSequenceではsortとの差は
        // 同値要素の元の順序だが、Tが単純な値なら通常は同じ。
        ds.stable_sort(0, ds.size());

        cerr << "stable_sort: ";
        ds.print();

        // unique:
        // std::uniqueと同じく「連続する重複」を削除。
        // したがって、重複をすべて消すには先にsortする。
        ds.unique(0, ds.size());

        cerr << "unique: ";
        ds.print();

        // erase_if:
        // pred(e)==true の要素を削除。
        DS a({1, 2, 3, 4, 5, 6});

        a.erase_if(0, a.size(), [](long long x) {
            return x % 2 == 0;
        });

        cerr << "erase_if(偶数を削除): ";
        a.print();

        // filter:
        // pred(e)==true の要素だけを残す。
        DS b({1, 2, 3, 4, 5, 6});

        b.filter(0, b.size(), [](long long x) {
            return x % 2 == 1;
        });

        cerr << "filter(奇数だけ残す): ";
        b.print();

        // partition:
        // predを満たす要素を前半へ移動。
        // 順序は保証されない。
        DS c({1, 2, 3, 4, 5, 6});

        c.partition(0, c.size(), [](long long x) {
            return x <= 3;
        });

        cerr << "partition(x<=3): ";
        c.print();

        // stable_partition:
        // predを満たす要素を前半へ。
        // それぞれのグループ内の相対順序を維持。
        DS d({1, 2, 3, 4, 5, 6});

        d.stable_partition(0, d.size(), [](long long x) {
            return x <= 3;
        });

        cerr << "stable_partition(x<=3): ";
        d.print();

        // map:
        // [l,r)の各要素をf(e)に置換。
        DS e({1, 2, 3, 4, 5});

        e.map(1, 4, [](long long x) {
            return x * 10;
        });

        cerr << "map(1, 4, x*10): ";
        e.print();
    }

    // ================================================================
    // kth / median / count / mode / mex / inversion
    // ================================================================
    {
        cerr << "\n[統計・順序系]\n";

        DS ds({3, 1, 4, 1, 5, 9, 2, 6});

        cerr << "列: ";
        ds.print();

        // kth_smallest:
        // [l,r)を一時vector化してnth_elementする。
        // kは0-based。 O(len)。
        cerr << "kth_smallest(0, 8, 3) = " << ds.kth_smallest(0, 8, 3) << '\n';

        // 【新規】kth_largest は大きい方から k 番目 (0-based)。
        cerr << "kth_largest(0, 8, 0) = " << ds.kth_largest(0, 8, 0) << " (最大値)\n";
        cerr << "kth_largest(0, 8, 1) = " << ds.kth_largest(0, 8, 1) << '\n';

        // median:
        // 下側中央値。
        // 長さが偶数なら小さい方の中央値。
        cerr << "median(0, 8) = " << ds.median(0, 8) << '\n';

        // count_value:
        // 区間内のxの個数を数える。 O(len)。
        cerr << "count_value(0, 8, 1) = " << ds.count_value(0, 8, 1) << '\n';

        // 【新規】count_if: 条件を満たす要素数。 O(len)。
        cerr << "count_if(0, 8, 偶数) = "
             << ds.count_if(0, 8, [](long long x){ return x % 2 == 0; }) << '\n';

        // mode:
        // 最頻値。同率なら値が小さい方。
        cerr << "mode(0, 8) = " << ds.mode(0, 8) << '\n';

        // mex:
        // Tが非負整数として扱えることを仮定。
        // 0以上の値について最小の「存在しない値」を返す。
        DS mex_ds({0, 1, 2, 4, 5});
        cerr << "mex = " << mex_ds.mex(0, mex_ds.size()) << '\n';

        // inversion:
        // [l,r)内の転倒数。
        // 内部で座標圧縮 + BITを使う。 O(len log len)。
        DS inv_ds({3, 1, 2, 5, 4});
        cerr << "inversion = " << inv_ds.inversion(0, inv_ds.size()) << '\n';
    }

    // ================================================================
    // to_vector
    // ================================================================
    {
        cerr << "\n[to_vector]\n";

        DS ds({10, 20, 30, 40, 50});

        // 全要素をvectorへ。
        auto all = ds.to_vector();

        // [l,r)だけvectorへ。
        auto part = ds.to_vector(1, 4);

        cerr << "to_vector(): ";
        for (auto x : all) cerr << x << ' ';
        cerr << '\n';

        cerr << "to_vector(1, 4): ";
        for (auto x : part) cerr << x << ' ';
        cerr << '\n';
    }

    // ================================================================
    // Debug / validate
    // ================================================================
    {
        cerr << "\n[Debug]\n";

        DS ds({5, 1, 4, 2, 3});

        ds.reverse(1, 4);
        ds.add(0, 3, 10);
        ds.insert(2, 99);

        cerr << "print(): ";
        ds.print();

        // 内部構造 (チャンク化treap) を表示。
        // 1 行が 1 節点 = 1 チャンクで、[ ] の中がその節点の持つ配列。
        // 競プロ中のバグ調査用。
        cerr << "dump():\n";
        ds.dump();

        // 【強化】validate() は subtree size / heap property に加えて
        // 集約値 (sum/min/max/gcd/xor/and/or/Kadane) が子から作り直した値と
        // 一致するかまで検査する。trueなら内部構造が正常。
        cerr << "validate() = " << (ds.validate() ? "true" : "false") << '\n';

        // check()はvalidate()の別名。
        cerr << "check() = " << (ds.check() ? "true" : "false") << '\n';

        // 【新規】要素数の目安が分かっているならプールを事前確保できる
        DS::reserve_pool(1 << 20);
        cerr << "reserve_pool(1<<20) 済み\n";

        // 【新規】チャンク化treapの内部指標。
        //   chunk_size() : 1 節点が持つ配列の基準長 CS (実際の長さは 1..2*CS)
        //   node_count() : 節点数。おおよそ size()/CS になる
        //   height()     : 木の高さ。おおよそ log2(node_count()) の定数倍
        // 「node_count() が size() の割に多い」「height() が極端に大きい」
        // ときだけ性能を疑えばよい。
        cerr << "\n[内部指標]\n";
        DS big(vector<long long>(100000, 1));
        cerr << "CS (chunk_size) = " << DS::chunk_size() << '\n';
        cerr << "size = " << big.size()
             << " / node_count = " << big.node_count()
             << " / height = " << big.height() << '\n';
        big.reverse(3, 100000 - 3);          // split/merge が走る操作
        cerr << "reverse 後: node_count = " << big.node_count()
             << " / height = " << big.height() << '\n';
        cerr << "1 節点あたりの平均要素数 = "
             << (double)big.size() / big.node_count() << " (CS 前後なら健全)\n";
    }

    cerr << "\n========== 動作確認終了 ==========\n";
}

// ---------------------------------------------------------------------------
// [テストA] チャンク境界の総当たりストレス (CS を極端に小さくする)
//   CS を 1〜5 にすると「チャンクの分割・消滅・併合」が毎回のように起きるので、
//   通常の CS=32/64/128 では滅多に通らない経路を短時間で踏み抜ける。
//   ・毎回 validate() を呼ぶので、部分木サイズ / 節点数 / heap 性質 /
//     チャンク長 (1 以上 2*CS 以下) / 集約の整合性がすべて検査される。
//   ・split/merge を使う操作 (reverse/rotate/move/swap_range/extract/splice) を
//     厚めに混ぜて、「新しいチャンク節点が正しい優先度の位置に入るか」を突く。
// ---------------------------------------------------------------------------
template <int CS>
void chunk_boundary_test(const string& name, int OPS, unsigned seed) {
    using DS = DynamicSequence<long long, dseq::AGG_ALL, CS>;
    mt19937 rng(seed);
    vector<long long> ref(12);
    for (auto& x : ref) x = (int)(rng() % 20) - 10;
    DS ds(ref);

    long long node_peak = 0, len_sum = 0;
    for (int op = 0; op < OPS; op++) {
        // 長さを [8,60] 程度に保つ (空・巨大に張り付かせない)
        if ((int)ref.size() < 8) {
            const int add = 8 - (int)ref.size() + (int)(rng() % 5);
            for (int i = 0; i < add; i++) {
                const int pos = (int)(rng() % (ref.size() + 1));
                const long long v = (int)(rng() % 20) - 10;
                ref.insert(ref.begin() + pos, v);
                ds.insert(pos, v);
            }
        }
        const int n = (int)ref.size();
        len_sum += n;
        const int kind = (int)(rng() % 14);
        if (kind == 0) {                       // insert
            const int pos = (int)(rng() % (n + 1));
            const long long v = (int)(rng() % 40) - 20;
            ref.insert(ref.begin() + pos, v);
            ds.insert(pos, v);
        } else if (kind == 1) {                // erase (1点)
            const int pos = (int)(rng() % n);
            ref.erase(ref.begin() + pos);
            ds.erase(pos);
        } else if (kind == 2) {                // insert(vector) : 一気に増やす
            const int pos = (int)(rng() % (n + 1));
            const int len = 1 + (int)(rng() % (3 * CS + 2));
            vector<long long> v(len);
            for (auto& e : v) e = (int)(rng() % 40) - 20;
            ref.insert(ref.begin() + pos, v.begin(), v.end());
            ds.insert(pos, v);
        } else if (kind == 3) {                // erase(l,r) : 一気に減らす
            const int l = (int)(rng() % n), r = l + 1 + (int)(rng() % (n - l));
            ref.erase(ref.begin() + l, ref.begin() + r);
            ds.erase(l, r);
        } else if (kind == 4) {                // reverse (短い区間 = 書き戻し経路)
            const int l = (int)(rng() % n);
            const int r = l + 1 + (int)(rng() % std::min<int>(n - l, 2 * CS));
            std::reverse(ref.begin() + l, ref.begin() + r);
            ds.reverse(l, r);
        } else if (kind == 5) {                // reverse (長い区間 = split/merge 経路)
            const int l = (int)(rng() % n), r = l + 1 + (int)(rng() % (n - l));
            std::reverse(ref.begin() + l, ref.begin() + r);
            ds.reverse(l, r);
        } else if (kind == 6) {                // rotate
            const int l = (int)(rng() % n), r = l + 1 + (int)(rng() % (n - l));
            const int k = (int)(rng() % (r - l));
            std::rotate(ref.begin() + l, ref.begin() + (r - k), ref.begin() + r);
            ds.rotate(l, r, k);
        } else if (kind == 7) {                // move
            const int l = (int)(rng() % n), r = l + 1 + (int)(rng() % (n - l));
            const int len = r - l;
            const int pos = (int)(rng() % (n - len + 1));
            vector<long long> seg(ref.begin() + l, ref.begin() + r);
            ref.erase(ref.begin() + l, ref.begin() + r);
            ref.insert(ref.begin() + pos, seg.begin(), seg.end());
            ds.move(l, r, pos);
        } else if (kind == 8) {                // swap_range
            if (n >= 4) {
                const int l1 = (int)(rng() % (n - 3));
                const int r1 = l1 + 1 + (int)(rng() % std::min(3, n - l1 - 2));
                const int l2 = r1 + (int)(rng() % (n - r1));
                const int r2 = l2 + 1 + (int)(rng() % (n - l2));
                vector<long long> s1(ref.begin() + l1, ref.begin() + r1);
                vector<long long> s2(ref.begin() + l2, ref.begin() + r2);
                vector<long long> nxt(ref.begin(), ref.begin() + l1);
                nxt.insert(nxt.end(), s2.begin(), s2.end());
                nxt.insert(nxt.end(), ref.begin() + r1, ref.begin() + l2);
                nxt.insert(nxt.end(), s1.begin(), s1.end());
                nxt.insert(nxt.end(), ref.begin() + r2, ref.end());
                ref = nxt;
                ds.swap_range(l1, r1, l2, r2);
            }
        } else if (kind == 9) {                // extract + splice (別の木との行き来)
            const int l = (int)(rng() % n), r = l + 1 + (int)(rng() % (n - l));
            DS taken = ds.extract(l, r);
            assert(taken.validate());
            vector<long long> seg(ref.begin() + l, ref.begin() + r);
            ref.erase(ref.begin() + l, ref.begin() + r);
            const int pos = (int)(rng() % (ref.size() + 1));
            ds.splice(pos, taken, 0, (int)taken.size());
            ref.insert(ref.begin() + pos, seg.begin(), seg.end());
            assert(taken.size() == 0);
        } else if (kind == 10) {               // add / assign / affine
            const int l = (int)(rng() % n), r = l + 1 + (int)(rng() % (n - l));
            const int w = (int)(rng() % 3);
            if (w == 0) { const long long x = (int)(rng()%9)-4; for (int i=l;i<r;i++) ref[i]+=x; ds.add(l,r,x); }
            else if (w == 1) { const long long x = (int)(rng()%9)-4; for (int i=l;i<r;i++) ref[i]=x; ds.assign(l,r,x); }
            else { const long long a=(int)(rng()%5)-2, b=(int)(rng()%7)-3; for (int i=l;i<r;i++) ref[i]=a*ref[i]+b; ds.affine(l,r,a,b); }
        } else if (kind == 11) {               // set / get
            const int pos = (int)(rng() % n);
            const long long v = (int)(rng() % 40) - 20;
            ref[pos] = v; ds.set(pos, v);
            assert(ds.get(pos) == v);
        } else if (kind == 12) {               // sort / unique (rebuild_range 経路)
            const int l = (int)(rng() % n), r = l + 1 + (int)(rng() % (n - l));
            if (rng() % 2) { std::sort(ref.begin()+l, ref.begin()+r); ds.sort(l, r); }
            else {
                vector<long long> v(ref.begin()+l, ref.begin()+r);
                v.erase(std::unique(v.begin(), v.end()), v.end());
                ref.erase(ref.begin()+l, ref.begin()+r);
                ref.insert(ref.begin()+l, v.begin(), v.end());
                ds.unique(l, r);
            }
        } else {                               // 区間クエリ (集約の整合性)
            const int l = (int)(rng() % n), r = l + 1 + (int)(rng() % (n - l));
            long long s = 0, mn = LLONG_MAX, mx = LLONG_MIN;
            for (int i = l; i < r; i++) { s += ref[i]; mn = std::min(mn, ref[i]); mx = std::max(mx, ref[i]); }
            assert(ds.sum(l, r) == s);
            assert(ds.min(l, r) == mn);
            assert(ds.max(l, r) == mx);
        }
        assert(ds.validate());
        assert((int)ds.size() == (int)ref.size());
        assert(ds.to_vector() == ref);
        node_peak = std::max(node_peak, (long long)ds.node_count());
    }
    cerr << "[DynamicSequence] chunk_boundary_test<CS=" << CS << "> OK (" << name << ", "
         << OPS << " ops, 平均長 " << (double)len_sum / OPS << ", 節点数の最大 " << node_peak << ")\n";
}

// ---------------------------------------------------------------------------
// [テストB] 木の健全性 (節点数と高さが理論どおりに収まっているか)
//   チャンク化 treap の速さは
//     ・節点数 ≒ N/CS  (チャンクがちゃんと詰まっている)
//     ・高さ   = O(log(節点数))  (heap 性質が守られていて偏っていない)
//   の 2 つが前提。split/merge を多用しても崩れないことを数値で確認する。
//   ※ heap 性質そのものは validate() 内で assert している。
// ---------------------------------------------------------------------------
void chunk_health_test() {
    constexpr int CS = 32;
    using DS = DynamicSequence<long long, dseq::AGG_BASIC, CS>;
    mt19937 rng(20240919);
    const int N = 200000;
    vector<long long> init(N);
    for (auto& x : init) x = (int)(rng() % 1000);
    DS ds(init);

    auto report = [&](const char* phase) {
        const int nodes = ds.node_count();
        const int height = ds.height();
        const double fill = (double)ds.size() / std::max(1, nodes);
        // 期待値: 平均チャンク長は CS/3 以上 (再構築のしきい値)、高さは 4*log2(節点数)+8 以内
        const double limit = 4.0 * std::log2(std::max(2, nodes)) + 8.0;
        cerr << "    " << phase << ": 節点数=" << nodes << " 平均チャンク長=" << fill
             << " 高さ=" << height << " (目安 " << (int)limit << " 以下)\n";
        assert(fill >= (double)CS / 3.0 - 1.0 && "チャンクが痩せすぎ (再構築が効いていない)");
        assert((double)height <= limit && "木が偏りすぎ (heap 性質が壊れている疑い)");
        assert(ds.validate());
    };
    report("構築直後   ");

    for (int i = 0; i < N; i++) ds.insert((int)(rng() % (ds.size() + 1)), (long long)(rng() % 1000));
    report("挿入 20 万 ");
    for (int i = 0; i < N; i++) ds.erase((int)(rng() % ds.size()));
    report("削除 20 万 ");
    for (int i = 0; i < N / 2; i++) {
        const int l = (int)(rng() % ds.size());
        const int r = l + 1 + (int)(rng() % (ds.size() - l));
        ds.reverse(l, r);
    }
    report("反転 10 万 ");
    for (int i = 0; i < N / 2; i++) {
        const int l = (int)(rng() % (ds.size() - 8));
        ds.reverse(l, l + 8);            // 短い反転 (書き戻し経路): 節点が増えないはず
    }
    report("短い反転   ");
    cerr << "[DynamicSequence] chunk_health_test OK\n";
}

// ---------------------------------------------------------------------------
// [テストC] 遅延タグ (ct / at) の実体化まわり
//   この実装は遅延を「子へ配る分 ct」と「自分の配列に未適用の分 at」に分ける。
//   そのため「配列を実体化せずに at を見ながら答える」読み出し経路
//   (get / sum / min / max の部分区間) が独自に存在する。
//   ここを反転タグ付きの状態で総当たりし、実体化する経路と一致することを確認する。
// ---------------------------------------------------------------------------
template <int CS>
void lazy_materialize_test(unsigned seed) {
    using DS = DynamicSequence<long long, dseq::AGG_ALL, CS>;
    mt19937 rng(seed);
    for (int iter = 0; iter < 300; iter++) {
        const int n = 1 + (int)(rng() % (6 * CS));
        vector<long long> ref(n);
        for (auto& x : ref) x = (int)(rng() % 100) - 50;
        DS ds(ref);
        for (int step = 0; step < 12; step++) {
            const int l = (int)(rng() % n), r = l + 1 + (int)(rng() % (n - l));
            const int w = (int)(rng() % 4);
            if (w == 0) { const long long x=(int)(rng()%11)-5; for(int i=l;i<r;i++) ref[i]+=x; ds.add(l,r,x); }
            else if (w == 1) { const long long x=(int)(rng()%11)-5; for(int i=l;i<r;i++) ref[i]=x; ds.assign(l,r,x); }
            else if (w == 2) { const long long x=(int)(rng()%5)-2; for(int i=l;i<r;i++) ref[i]*=x; ds.multiply(l,r,x); }
            else { std::reverse(ref.begin()+l, ref.begin()+r); ds.reverse(l,r); }

            // 実体化されていない可能性のある経路を先に叩く
            for (int q = 0; q < 6; q++) {
                const int a = (int)(rng() % n), b = a + 1 + (int)(rng() % (n - a));
                long long s = 0, mn = LLONG_MAX, mx = LLONG_MIN;
                for (int i = a; i < b; i++) { s += ref[i]; mn = std::min(mn, ref[i]); mx = std::max(mx, ref[i]); }
                assert(ds.sum(a, b) == s);
                assert(ds.min(a, b) == mn);
                assert(ds.max(a, b) == mx);
                const int p = (int)(rng() % n);
                assert(ds.get(p) == ref[p]);
            }
            // そのあとで実体化する経路と突き合わせる
            assert(ds.to_vector() == ref);
            assert(ds.validate());
        }
    }
    cerr << "[DynamicSequence] lazy_materialize_test<CS=" << CS << "> OK\n";
}

// ---------------------------------------------------------------------------
// [テストD] 短い並べ替えの境界 (書き戻し経路 ⇔ split/merge 経路の切り替え)
//   reverse / rotate は長さ 2*CS 以下なら「値を書き戻す」方式、
//   それより長ければ split/merge 方式に分岐する。
//   境界のちょうど前後で結果が一致することを全長さについて確認する。
// ---------------------------------------------------------------------------
template <int CS>
void small_permutation_boundary_test() {
    using DS = DynamicSequence<long long, dseq::AGG_ALL, CS>;
    mt19937 rng(98765);
    const int n = 6 * CS + 5;
    vector<long long> base(n);
    for (auto& x : base) x = (int)(rng() % 200) - 100;
    // 長さ 1 .. 2*CS+3 まで全部試す (2*CS が分岐の境界)
    for (int len = 1; len <= 2 * CS + 3; len++) {
        for (int l : {0, 1, CS, CS + 1, 2 * CS, n - len}) {
            if (l < 0 || l + len > n) continue;
            {   // reverse
                vector<long long> ref = base;
                DS ds(base);
                std::reverse(ref.begin() + l, ref.begin() + l + len);
                ds.reverse(l, l + len);
                assert(ds.to_vector() == ref);
                assert(ds.validate());
            }
            for (int k : {0, 1, len / 2, len - 1}) {   // rotate
                vector<long long> ref = base;
                DS ds(base);
                std::rotate(ref.begin() + l, ref.begin() + l + len - k, ref.begin() + l + len);
                ds.rotate(l, l + len, k);
                assert(ds.to_vector() == ref);
                assert(ds.validate());
            }
        }
    }
    // 全体反転 (根にタグを付けるだけの経路)
    {
        vector<long long> ref = base;
        DS ds(base);
        std::reverse(ref.begin(), ref.end());
        ds.reverse(0, n);
        assert(ds.to_vector() == ref);
        ds.reverse();
        std::reverse(ref.begin(), ref.end());
        assert(ds.to_vector() == ref);
        assert(ds.validate());
    }
    cerr << "[DynamicSequence] small_permutation_boundary_test<CS=" << CS << "> OK\n";
}

// ---------------------------------------------------------------------------
// [テストE] チャンク長の自動調整 (分割と消滅)
//   ・同じ位置に挿し続ける  -> チャンクが 2*CS を超えて分割されるのを繰り返す
//   ・同じ位置から消し続ける -> チャンクが空になって節点が消えるのを繰り返す
//   どちらも「毎回 validate()」で長さ 1..2*CS の不変条件を検査する。
// ---------------------------------------------------------------------------
template <int CS>
void chunk_split_merge_test() {
    using DS = DynamicSequence<long long, dseq::AGG_BASIC, CS>;
    // (1) 先頭に挿し続ける
    {
        DS ds;
        vector<long long> ref;
        for (int i = 0; i < 10 * CS; i++) {
            ds.insert(0, i); ref.insert(ref.begin(), i);
            assert(ds.validate());
            assert(ds.to_vector() == ref);
        }
        // (2) 先頭から消し続ける
        while (!ref.empty()) {
            ds.erase(0); ref.erase(ref.begin());
            assert(ds.validate());
            assert(ds.to_vector() == ref);
        }
        assert(ds.empty());
    }
    // (3) 真ん中に挿し続ける (1 つのチャンクが何度も分割される)
    {
        DS ds; vector<long long> ref;
        for (int i = 0; i < 10 * CS; i++) {
            const int pos = (int)ref.size() / 2;
            ds.insert(pos, i); ref.insert(ref.begin() + pos, i);
            assert(ds.validate());
        }
        assert(ds.to_vector() == ref);
        // (4) 真ん中から消し続ける
        while (!ref.empty()) {
            const int pos = (int)ref.size() / 2;
            ds.erase(pos); ref.erase(ref.begin() + pos);
            assert(ds.validate());
        }
    }
    // (5) 末尾に積んでから全部消す (push_back / pop_back)
    {
        DS ds; vector<long long> ref;
        for (int i = 0; i < 10 * CS; i++) { ds.push_back(i); ref.push_back(i); }
        assert(ds.to_vector() == ref);
        assert(ds.validate());
        while (!ref.empty()) { ds.pop_back(); ref.pop_back(); assert(ds.validate()); }
    }
    cerr << "[DynamicSequence] chunk_split_merge_test<CS=" << CS << "> OK\n";
}
void solv(){
    // ---- ランダムテスト -------------------------------------------------
    stress_test(800000, 6, 60, 12345);      // 小さい木: split/merge の境界を潰す
    stress_test(200000, 150, 400, 777);     // 深い木: 再帰の深いパスを通す
    stress_test_nonlinear();
    stress_test_multiply_agg();
    stress_test_extra();
    stress_test_binary_search();
    stress_test_agg<dseq::AGG_BASIC>("AGG_BASIC", 400000);
    stress_test_agg<dseq::AGG_GCD>("AGG_GCD", 400000);
    stress_test_agg<dseq::AGG_BITWISE>("AGG_BITWISE", 400000);
    stress_test_agg<dseq::AGG_KADANE>("AGG_KADANE", 400000);
    stress_test_agg<dseq::AGG_ALL>("AGG_ALL", 400000);
    stress_test_agg_kadane_only();
    test_new_api_and_edge();
    chunk_boundary_test<1>("CS=1 (1要素=1節点: 通常のtreapと同じ形)", 60000, 111);
    chunk_boundary_test<2>("CS=2", 60000, 222);
    chunk_boundary_test<3>("CS=3", 60000, 333);
    chunk_boundary_test<5>("CS=5", 60000, 555);
    chunk_boundary_test<16>("CS=16", 60000, 1616);
    chunk_split_merge_test<1>();
    chunk_split_merge_test<3>();
    chunk_split_merge_test<8>();
    lazy_materialize_test<2>(4242);
    lazy_materialize_test<8>(4343);
    lazy_materialize_test<32>(4444);
    small_permutation_boundary_test<2>();
    small_permutation_boundary_test<8>();
    chunk_health_test();
    large_scale_test();
    cerr << "\n===== ALL TESTS PASSED =====\n\n";

    // ---- 動作確認 -------------------------------------------------------
    test_DynamicSequence();
    return;
}

int main(){
    BUFF;
    int T=1;
    //cin>>T;
    while(T--)solv();
    return 0;
}
