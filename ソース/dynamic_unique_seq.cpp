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

/* =====================================================================================
 *  DynUniqueSeq : 「全要素が相異なる long long」を要素に持つ動的配列
 *                 (AtCoder 用テンプレート / C++17 以降)
 * =====================================================================================
 *
 *  ■ このテンプレートは何か
 *      「位置を持つ動的配列」＋「全要素が相異なることを利用した値の検索・順序検索」
 *      を 1 つのデータ構造にまとめたもの。内部は 3 層構造:
 *
 *        (1) 位置木  : implicit treap (遅延 reverse 付き)。配列としての順序と
 *                      区間集約 (sum / min / max) を担当。
 *        (2) 値木    : 値をキーにした treap (順序統計木)。rank / k 番目 /
 *                      prev / next / 値の総和 / 最小ギャップ を担当。
 *        (3) ハッシュ: 値 -> ノード番号 の O(1) 逆引き。index_of を高速化。
 *
 *      (1) と (2) は「同じノード番号」を共有している。したがって
 *          値 --(3)--> ノード --(1)の親を辿る--> 位置
 *          位置 --(1)を降りる--> ノード --(2)の親を辿る--> 順位
 *      が両方向 O(log N) で行える。ここがこのテンプレートの核心。
 *
 *  ■ 設計上のポイント (機能候補 #19 の「reverse すると 値->index が壊れる」問題)
 *      値 -> index を「保存」せず「その都度 親ポインタを辿って計算」する方式にした。
 *      これにより区間 reverse を遅延伝播で O(log N) のまま行っても値→位置の
 *      整合性は自動的に保たれる (位置は木の形からしか作られないため)。
 *      さらに、区間クエリと index 計算では遅延フラグを「下ろさず」に
 *      反転パリティを引数で持ち回るため、読み取り系は木への書き込みが一切発生しない
 *      (キャッシュに優しく、const で呼べる)。
 *
 *  ■ 前提・約束 (重要)
 *      ・要素の値はすべて相異なる long long。重複を入れた時の動作は未定義。
 *        (デバッグ時は  #define DUS_DEBUG  すると assert で検出できる)
 *      ・位置は 0-indexed、区間はすべて半開区間 [l, r)。
 *      ・「見つからない」ときの返り値:
 *            位置(index) は -1
 *            値は「min を探す系」なら INF、「max を探す系」なら NINF
 *      ・range_sum のオーバーフローは面倒を見ない (必要なら sum を __int128 に)。
 *
 *  ■ 命名規則 (速度のため 3 種類用意してある)
 *        xxx(...)              … 値だけ返す           (速い)
 *        xxx_with_index(...)   … {値, 位置} を返す     (位置計算の分 +O(log N))
 *        xxx_index(...)        … 位置だけ返す
 *      値だけで済む場面で位置計算のコストを払わなくて済むようにしてある。
 *
 * -------------------------------------------------------------------------------------
 *  ■ 計算量一覧 (N = 要素数。treap なので期待計算量)
 *
 *   [基本]
 *     size(), empty(), total_sum()                             O(1)
 *     get(i) / operator[](i) / front() / back()                O(log N)
 *     to_vector()                                              O(N)
 *     build(vector) / コンストラクタ                            O(N log N)  (ソート分)
 *
 *   [挿入・削除・変更]
 *     insert(i,x), erase(i), push_back/front, pop_back/front   O(log N)
 *     erase_value(x)                                           O(log N)
 *     set(i,x)                                                 O(log N)
 *     set_value(old,new)                                       O(log N)
 *
 *   [値 -> 位置]
 *     contains(x), count(x)                                    O(1)   期待
 *     index_of(x)                                              O(log N)
 *     contains_in_range(l,r,x), range_index_of(l,r,x)          O(log N)
 *
 *   [区間集約]
 *     range_sum(l,r)                                           O(log N)
 *     range_info(l,r)  … sum/min/max を 1 回の降下でまとめて     O(log N)
 *     range_min/max(l,r)                       (値のみ)         O(log N)
 *     range_min/max_with_index(l,r)                            O(log N)
 *     range_min/max_index(l,r)                                 O(log N)
 *     range_second_min/max(l,r)                                O(log N)  (min を 3 回)
 *     range_max_min_diff(l,r)                                  O(log N)
 *     get_min(), get_max()                     (値のみ)         O(1)
 *     get_min_with_index(), min_index() 等                     O(log N)
 *
 *   [値の順序 (set 相当)]
 *     rank_of_value(x)   (= x 未満の個数)                       O(log N)
 *     count_less/le/greater/ge(x), count_value_range(lo,hi)    O(log N)
 *     sum_less(x), sum_value_range(lo,hi), sum_of_k_smallest(k) O(log N)
 *     kth_smallest(k), kth_largest(k)                          O(log N)
 *     kth_smallest_index(k) 等                                 O(log N)
 *     prev_value(x), next_value(x)                             O(log N)
 *     lower_bound_value(x), upper_bound_value(x)               O(log N)
 *     rank_of_index(i)                                         O(log N)
 *     prev_index_by_value(i), next_index_by_value(i)           O(log N)
 *     minimum_gap(), minimum_gap_pair(), minimum_gap_indices() O(1) / O(log N)
 *
 *   [並べ替え]
 *     swap_at(i,j)                                             O(log N)
 *     swap_values(x,y)                                         O(log N)
 *     move_element(from,to)                                    O(log N)
 *     reverse_range(l,r)                                       O(log N)
 *     rotate_range(l,m,r)                                      O(log N)
 *     move_range(l,r,pos)                                      O(log N)
 *
 *   [まとめて挿入・削除・取り出し]
 *     erase_range(l,r)                                         O(log N + k log N)
 *     insert_range(pos, vector)                                O(k log N)
 *     range_to_vector(l,r)                                     O(log N + k)
 *
 *   [値の昇順と位置の一括取得]
 *     value_order()   … 値の昇順に並べたときの元の位置の列        O(N)
 *     sorted_values() … 値を昇順に並べた配列                     O(N)
 *
 *   [区間の k 番目 (k が小さいとき専用)]
 *     range_k_smallest(l,r,k), range_k_largest(l,r,k)          O(k (log N + log k))
 *     range_kth_smallest(l,r,k), range_kth_largest(l,r,k)      O(k (log N + log k))
 *
 *   [その他]
 *     mex(from=0)  … from 以上で配列に無い最小の整数             O(log^2 N)
 *     max_right(l,f), min_left(r,f) … ACL 風の区間二分探索       O(log N)
 *
 * -------------------------------------------------------------------------------------
 *  ■ 機能候補.txt に無いが追加した機能 (どれも他の操作の定数倍を悪化させない)
 *      range_info(l,r)          sum/min/max を 1 回の降下でまとめて取得 (別々に呼ぶ約 2 倍速)
 *      sum_less(x) / sum_le(x) / sum_value_range(lo,hi)   値で絞った総和
 *      sum_of_k_smallest(k) / sum_of_k_largest(k)         小さい/大きい方 k 個の総和
 *      count_value_range(lo,hi)                           値が [lo,hi] の個数
 *      swap_values(x,y)         値 2 つの位置を交換 (位置を求めずに済むので速い)
 *      erase_range(l,r) / insert_range(pos,vec) / range_to_vector(l,r)   まとめて操作
 *      value_order() / sorted_values()                    値の昇順の一括取得 (O(N))
 *      range_k_smallest/largest, range_kth_smallest/largest(l,r,k)  下記参照
 *      mex(from=0)              from 以上で配列に無い最小の整数 (O(log^2 N))
 *      max_right(l,f) / min_left(r,f)                     ACL 風の区間二分探索 (O(log N))
 *
 *      「追加を検討したが入れなかったもの」
 *        ・区間の最大部分和 (Kadane 集約): ノードに prefix/suffix/best を持たせる必要があり
 *          PN が 64 -> 96 byte に膨らんで全操作が遅くなる。値が相異なる話とも無関係なので不採用。
 *        ・全要素への一律加算 (offset を 1 個持てば O(1) で可能): 全 API に加減算が入り
 *          ハッシュ・値木のキーの扱いも複雑になる。必要になったときに足すのが良い。
 *        ・2 つの DynUniqueSeq の連結: 位置木は O(log N) だが値木の併合が O(n log(m/n)) で
 *          実装量に見合わないため不採用。
 *
 * -------------------------------------------------------------------------------------
 *  ■ 「入れなかった機能」と理由  (機能候補 #21〜#25)
 *      range_count_less(l,r,x) / range_prev(l,r,x) / range_next(l,r,x) /
 *      range_lower_bound(l,r,x) / 完全な range_kth_smallest(l,r,k)
 *
 *      これらは「位置の区間」×「値の順序」の 2 次元クエリで、動的配列 (任意位置の
 *      挿入削除あり) の上で実現するには "平衡二分探索木の各ノードに平衡二分探索木を
 *      持たせる" 構造が必要になり、
 *          時間 O(log^2 N) 〜 O(log^3 N)、空間 O(N log N)、定数倍も非常に大きい
 *      となって、他の全操作の定数倍まで巻き添えで重くなる。
 *      「速さ優先」という方針に反するため今回は不採用とした。
 *      (必要になったら Wavelet Matrix / マージソート木を別テンプレートとして使うのが良い)
 *
 *      ただし #22/#23 の range_kth_smallest / range_kth_largest は
 *      「k が小さいとき専用」なら区間 min を k 回取り直すだけで O(k (log N + log k)) で
 *      できるので、そちらは採用した (range_k_smallest / range_kth_smallest)。
 *
 *      なお 機能候補 #26 の range_index_of(l,r,x) / contains_in_range(l,r,x) は
 *      「値が一意」なので index_of(x) してから範囲判定するだけで O(log N)。採用済み。
 *      機能候補 #8 range_size(l,r) = r-l、#33 maximum_gap = max-min は自明なので不採用。
 *      機能候補 #27 の prev_index(i)/next_index(i) は i-1 / i+1 そのものなので不採用
 *      (値としての前後 prev_index_by_value / next_index_by_value は採用済み)。
 *
 * -------------------------------------------------------------------------------------
 *  ■ 速度のための実装上の工夫
 *      ・ノードはポインタでなく vector 上の int 添字。PN / VN ともちょうど 64 byte
 *        (= 1 キャッシュライン) に収まるようメンバを配置。
 *      ・区間 min/max は {値, ノード番号} を集約に持つ。値を集約に持つことで
 *        pull() 中に他ノードの val を読みに行くランダムアクセスが発生しない。
 *      ・0 番を番兵ノード (sz=0, sum=0, mnv=INF, mxv=NINF) にして null 判定の分岐を削減。
 *      ・読み取り系 (区間クエリ・位置計算・値検索) は遅延フラグを下ろさず
 *        反転パリティを引数で持ち回る -> 書き込み 0 回、const 可。
 *      ・乱数は xorshift32。値→ノードのハッシュは splitmix64 + 線形探索 +
 *        backward-shift 削除 (tombstone なし)。ハッシュ種は起動ごとにランダム化 (hack 対策)。
 *      ・set / set_value / swap_at は木の形を変えないので split/merge を使わず
 *        「親を辿って pull し直す」だけ。sz/sum/min/max は反転で不変なので
 *        遅延フラグが残っていても正しい。
 *      ・erase は split/merge を使わず該当ノードを直接外す。特に erase_value(x) は
 *        「値 -> ノード」で直接消せるので位置を求める必要すら無い。
 *      ・reverse を一度も使っていない間は index 計算を親 1 往復に短縮 (has_rev)。
 *      ・sum/min/max が同時に欲しいときは range_info(l,r) が木を 1 回降りるだけで済む。
 *
 * -------------------------------------------------------------------------------------
 *  ■ 正しさの担保 (このファイルの下半分のテストで確認済み)
 *      ・全 API を愚直実装 (vector を舐めるだけ) と突き合わせるランダムテストを 5 パターン。
 *        1 回の実行で約 300 万項目を比較する。
 *      ・毎回、木の内部不変条件を全数検査する:
 *          親ポインタの整合性 / treap の heap 条件 / sz・sum・min・max の再計算一致 /
 *          値木が昇順であること / 最小ギャップの再計算一致 /
 *          ハッシュ表が全要素を正しいノードに対応づけていること /
 *          全要素について index_of が正しい位置を返すこと
 *      ・「わざとバグを入れた 25 種類の改変版」がすべてこのテストで検出されることを確認済み
 *        (遅延伝播の抜け、親ポインタの張り忘れ、集約の取りこぼし、ハッシュ削除の誤り、
 *         区間二分探索の降下ミス、mex の境界ずれ 等)
 *      ・AddressSanitizer / UndefinedBehaviorSanitizer でも警告なし。
 *
 * -------------------------------------------------------------------------------------
 *  ■ 使い方 (最小例)
 *
 *      DynUniqueSeq s({50, 10, 80, 30});   // または s.build(vector)
 *      s.insert(2, 70);                    // [50,10,70,80,30]
 *      auto [mn, pos] = s.range_min_with_index(1, 4);  // 区間 [1,4) の最小値と位置
 *      int  p  = s.index_of(80);           // 値 80 が今どこにいるか
 *      long long v = s.kth_smallest(1);    // 2 番目に小さい値
 *      s.reverse_range(0, 3);              // 区間反転しても index_of は正しいまま
 *
 *  ■ このファイルの実行方法
 *
 *      g++ -std=gnu++20 -O2 -o dus dyn_unique_seq.cpp
 *      ./dus            … デモ -> ランダムテスト -> ベンチマーク を全部実行
 *      ./dus demo       … 動作確認デモだけ (出力を見れば使い方を思い出せる)
 *      ./dus test       … ランダムテストだけ
 *      ./dus test 12345 … シード指定でランダムテスト (失敗再現用)
 *      ./dus bench      … ベンチマークだけ
 *
 *      AtCoder へ提出するときは下の
 *          ===== テンプレート本体 ここから =====  〜  ===== ここまで =====
 *      の部分だけをコピーすればよい。
 * ===================================================================================== */
// #define DUS_DEBUG   // 有効にすると「値の重複」「範囲外」を assert で検出する(遅くなる)

// ===================== テンプレート本体 ここから =====================

struct DynUniqueSeq {
    using ll  = long long;
    using ull = unsigned long long;
    static constexpr ll INF  = std::numeric_limits<ll>::max();   // 「min が無い」を表す
    static constexpr ll NINF = std::numeric_limits<ll>::min();   // 「max が無い」を表す

    // ---- 位置木 (implicit treap) のノード : 64 byte ----
    struct PN {
        int l, r, par, sz;      // 左右の子 / 親 / 部分木サイズ
        unsigned pri;           // treap 優先度
        int mni, mxi;           // 部分木の最小値 / 最大値を持つ "ノード番号"
        unsigned char rev;      // 遅延反転フラグ
        ll val, sum, mnv, mxv;  // 自身の値 / 部分木の和 / 最小値 / 最大値
    };
    // ---- 値木 (treap, key = 値) のノード : 64 byte ----
    struct VN {
        int l, r, par, sz;
        unsigned pri;
        ll key, minv, maxv, vsum;  // 値 / 部分木の最小値 / 最大値 / 総和
        ull gap;                   // 部分木内で「値が隣り合う 2 要素」の差の最小値
    };

    std::vector<PN> P;
    std::vector<VN> V;
    std::vector<int> freelist;
    std::vector<int> pbuf;     // 根からノードまでの経路を一時的に入れる作業領域
    int proot = 0, vroot = 0;
    bool has_rev = false;      // 一度でも reverse したか (未使用なら index 計算を 1 往復に短縮)
    unsigned rngv;

    // ---- 値 -> ノード番号 のハッシュ表 (open addressing, hv[i]==0 が空) ----
    std::vector<ll> hk;
    std::vector<int> hv;
    size_t hmask = 0, hcnt = 0;
    ull hseed;

    // ---------------------------------------------------------------- 乱数 / ハッシュ
    inline unsigned next_rand() {
        unsigned x = rngv; x ^= x << 13; x ^= x >> 17; x ^= x << 5; return rngv = x;
    }
    static inline ull splitmix(ull x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    inline size_t hslot(ll x) const {
        size_t i = (size_t)(splitmix((ull)x ^ hseed) & hmask);
        while (hv[i] && hk[i] != x) i = (i + 1) & hmask;
        return i;
    }
    void hinit(size_t cap) {
        size_t s = 8; while (s < cap) s <<= 1;
        hk.assign(s, 0); hv.assign(s, 0); hmask = s - 1; hcnt = 0;
    }
    void hgrow() {
        std::vector<ll> ok; std::vector<int> ov;
        ok.swap(hk); ov.swap(hv);
        size_t ns = (ov.size()) * 2;
        hk.assign(ns, 0); hv.assign(ns, 0); hmask = ns - 1; hcnt = 0;
        for (size_t i = 0; i < ov.size(); i++) if (ov[i]) hput(ok[i], ov[i]);
    }
    inline void hput(ll x, int id) {          // 追加 or 上書き
        if ((hcnt + 1) * 2 > hmask + 1) hgrow();
        size_t i = hslot(x);
        if (!hv[i]) { hk[i] = x; hcnt++; }
        hv[i] = id;
    }
    inline int hget(ll x) const { return hv[hslot(x)]; }
    void hdel(ll x) {                          // backward shift 削除 (tombstone 無し)
        size_t i = hslot(x);
        if (!hv[i]) return;
        hcnt--;
        size_t j = i;
        while (true) {
            j = (j + 1) & hmask;
            if (!hv[j]) break;
            size_t k = (size_t)(splitmix((ull)hk[j] ^ hseed) & hmask);
            bool keep = (i <= j) ? (i < k && k <= j) : (i < k || k <= j);
            if (keep) continue;
            hk[i] = hk[j]; hv[i] = hv[j]; i = j;
        }
        hv[i] = 0;
    }

    // ---------------------------------------------------------------- 構築 / 破棄
    void init_sentinel() {
        P.assign(1, PN{});  V.assign(1, VN{});
        P[0].l = P[0].r = P[0].par = 0; P[0].sz = 0; P[0].pri = 0;
        P[0].mni = P[0].mxi = 0; P[0].rev = 0;
        P[0].val = 0; P[0].sum = 0; P[0].mnv = INF; P[0].mxv = NINF;
        V[0].l = V[0].r = V[0].par = 0; V[0].sz = 0; V[0].pri = 0;
        V[0].key = 0; V[0].minv = INF; V[0].maxv = NINF; V[0].vsum = 0;
        V[0].gap = ~0ULL;
    }
    DynUniqueSeq() {
        ull t = (ull)std::chrono::steady_clock::now().time_since_epoch().count();
        rngv = (unsigned)(splitmix(t) | 1u);
        hseed = splitmix(t ^ 0x5deece66dULL);
        init_sentinel(); hinit(16);
    }
    explicit DynUniqueSeq(const std::vector<ll>& a) : DynUniqueSeq() { build(a); }

    void clear() {
        init_sentinel(); freelist.clear(); proot = vroot = 0; has_rev = false; hinit(16);
    }

    int new_node(ll x) {
        int id;
        if (!freelist.empty()) { id = freelist.back(); freelist.pop_back(); }
        else { id = (int)P.size(); P.emplace_back(); V.emplace_back(); }
        PN& p = P[id];
        p.l = p.r = p.par = 0; p.sz = 1; p.pri = next_rand();
        p.mni = p.mxi = id; p.rev = 0;
        p.val = p.sum = p.mnv = p.mxv = x;
        VN& v = V[id];
        v.l = v.r = v.par = 0; v.sz = 1; v.pri = next_rand();
        v.key = v.minv = v.maxv = v.vsum = x; v.gap = ~0ULL;
        return id;
    }
    inline void free_node(int id) { freelist.push_back(id); }

    // ---------------------------------------------------------------- 位置木の基本操作
    inline void pull(int t) {
        const int l = P[t].l, r = P[t].r;
        P[t].sz  = P[l].sz + P[r].sz + 1;
        P[t].sum = P[l].sum + P[t].val + P[r].sum;
        ll v = P[t].val; int i = t;
        if (P[l].mnv < v) { v = P[l].mnv; i = P[l].mni; }
        if (P[r].mnv < v) { v = P[r].mnv; i = P[r].mni; }
        P[t].mnv = v; P[t].mni = i;
        v = P[t].val; i = t;
        if (P[l].mxv > v) { v = P[l].mxv; i = P[l].mxi; }
        if (P[r].mxv > v) { v = P[r].mxv; i = P[r].mxi; }
        P[t].mxv = v; P[t].mxi = i;
    }
    inline void push(int t) {
        if (P[t].rev) {
            const int l = P[t].l, r = P[t].r;
            P[t].l = r; P[t].r = l;
            if (l) P[l].rev ^= 1;
            if (r) P[r].rev ^= 1;
            P[t].rev = 0;
        }
    }
    // sz/sum/min/max は反転で不変なので、遅延フラグが残っていても親を辿って直して良い
    inline void update_up(int t) { while (t) { pull(t); t = P[t].par; } }

    void split(int t, int k, int& a, int& b) {      // [0,k) と [k,sz)
        if (!t) { a = b = 0; return; }
        push(t);
        if (P[P[t].l].sz >= k) {
            split(P[t].l, k, a, P[t].l);
            if (P[t].l) P[P[t].l].par = t;
            pull(t); b = t;
        } else {
            split(P[t].r, k - P[P[t].l].sz - 1, P[t].r, b);
            if (P[t].r) P[P[t].r].par = t;
            pull(t); a = t;
        }
    }
    inline void split_root(int t, int k, int& a, int& b) {
        split(t, k, a, b);
        if (a) P[a].par = 0;
        if (b) P[b].par = 0;
    }
    int merge(int a, int b) {
        if (!a || !b) return a | b;
        if (P[a].pri > P[b].pri) {
            push(a); int c = merge(P[a].r, b); P[a].r = c; P[c].par = a; pull(a); return a;
        } else {
            push(b); int c = merge(a, P[b].l); P[b].l = c; P[c].par = b; pull(b); return b;
        }
    }

    // 根 -> u の経路の遅延反転を下ろす (reverse を使っていなければ何もしない)
    void push_path(int u) {
        if (!has_rev) return;
        pbuf.clear();
        for (int v = u; v; v = P[v].par) pbuf.push_back(v);
        for (int i = (int)pbuf.size() - 1; i >= 0; i--) push(pbuf[i]);
    }
    // ノード u を木から直接外す (split/merge で位置から消すより速い)
    void erase_node(int u) {
        push_path(u);
        const int t = merge(P[u].l, P[u].r);
        const int p = P[u].par;
        if (t) P[t].par = p;
        if (!p) { proot = t; if (t) P[t].par = 0; }
        else {
            if (P[p].l == u) P[p].l = t; else P[p].r = t;
            update_up(p);
        }
        P[u].l = P[u].r = P[u].par = 0;
        v_erase(u); hdel(P[u].val); free_node(u);
    }

    // ---------------------------------------------------------------- 値木の基本操作
    inline void pull_v(int t) {
        const int l = V[t].l, r = V[t].r;
        V[t].sz   = V[l].sz + V[r].sz + 1;
        V[t].vsum = V[l].vsum + V[t].key + V[r].vsum;
        V[t].minv = l ? V[l].minv : V[t].key;
        V[t].maxv = r ? V[r].maxv : V[t].key;
        ull g = V[l].gap < V[r].gap ? V[l].gap : V[r].gap;
        if (l) { ull d = (ull)V[t].key - (ull)V[l].maxv; if (d < g) g = d; }
        if (r) { ull d = (ull)V[r].minv - (ull)V[t].key; if (d < g) g = d; }
        V[t].gap = g;
    }
    inline void update_up_v(int t) { while (t) { pull_v(t); t = V[t].par; } }

    void v_split(int t, ll k, int& a, int& b) {     // key<k を a、key>=k を b
        if (!t) { a = b = 0; return; }
        if (V[t].key < k) {
            v_split(V[t].r, k, V[t].r, b);
            if (V[t].r) V[V[t].r].par = t;
            pull_v(t); a = t;
        } else {
            v_split(V[t].l, k, a, V[t].l);
            if (V[t].l) V[V[t].l].par = t;
            pull_v(t); b = t;
        }
    }
    int v_merge(int a, int b) {
        if (!a || !b) return a | b;
        if (V[a].pri > V[b].pri) {
            int c = v_merge(V[a].r, b); V[a].r = c; V[c].par = a; pull_v(a); return a;
        } else {
            int c = v_merge(a, V[b].l); V[b].l = c; V[c].par = b; pull_v(b); return b;
        }
    }
    inline void v_isolate(int u) {                  // 単独ノードとして初期化
        V[u].l = V[u].r = V[u].par = 0; pull_v(u);
    }
    void v_insert(int u) {
        int a, b;
        v_split(vroot, V[u].key, a, b);
        if (a) V[a].par = 0;
        if (b) V[b].par = 0;
        vroot = v_merge(v_merge(a, u), b);
        V[vroot].par = 0;
    }
    void v_erase(int u) {
        int t = v_merge(V[u].l, V[u].r);
        int p = V[u].par;
        if (t) V[t].par = p;
        if (!p) { vroot = t; if (t) V[t].par = 0; }
        else {
            if (V[p].l == u) V[p].l = t; else V[p].r = t;
            update_up_v(p);
        }
        V[u].l = V[u].r = V[u].par = 0;
    }

    // ---------------------------------------------------------------- 位置 <-> ノード
    // 位置 i のノード番号。遅延フラグは下ろさず反転パリティ f を持ち回る (書き込み 0)
    int node_at(int i) const {
        if (i < 0 || i >= P[proot].sz) return 0;
        int t = proot, f = P[proot].rev;
        while (true) {
            const int L = f ? P[t].r : P[t].l;
            const int ls = P[L].sz;
            if (i < ls)       { t = L; f ^= P[L].rev; }
            else if (i == ls) { return t; }
            else {
                i -= ls + 1;
                const int R = f ? P[t].l : P[t].r;
                t = R; f ^= P[R].rev;
            }
        }
    }
    // ノード u の現在の位置。親を 2 回辿る (1 回目で反転パリティの総和 F を求める)
    int index_of_node(int u) const {
        if (!has_rev) {                     // reverse 未使用なら 1 往復で済む (こちらが普通)
            int idx = P[P[u].l].sz;
            for (int c = u, p = P[u].par; p; c = p, p = P[p].par)
                if (c == P[p].r) idx += P[P[p].l].sz + 1;
            return idx;
        }
        int F = 0;
        for (int v = u; v; v = P[v].par) F ^= P[v].rev;
        int idx = P[F ? P[u].r : P[u].l].sz;
        int h = P[u].rev;                       // h = rev[c..u] の xor
        for (int c = u, p = P[u].par; p; c = p, p = P[p].par) {
            const int fp = F ^ h;               // = rev[root..p] の xor
            const int lc = fp ? P[p].r : P[p].l;
            if (c != lc) idx += P[lc].sz + 1;
            h ^= P[p].rev;
        }
        return idx;
    }

    // ---------------------------------------------------------------- 区間クエリ (読み取りのみ)
    ll q_sum(int t, int l, int r, int f) const {
        if (l <= 0 && r >= P[t].sz) return P[t].sum;
        const int L = f ? P[t].r : P[t].l, R = f ? P[t].l : P[t].r;
        const int ls = P[L].sz;
        ll res = 0;
        if (l < ls)               res += q_sum(L, l, r < ls ? r : ls, f ^ P[L].rev);
        if (l <= ls && ls < r)    res += P[t].val;
        if (r > ls + 1)           res += q_sum(R, l - ls - 1 > 0 ? l - ls - 1 : 0, r - ls - 1, f ^ P[R].rev);
        return res;
    }
    std::pair<ll,int> q_min(int t, int l, int r, int f) const {
        if (l <= 0 && r >= P[t].sz) return {P[t].mnv, P[t].mni};
        const int L = f ? P[t].r : P[t].l, R = f ? P[t].l : P[t].r;
        const int ls = P[L].sz;
        std::pair<ll,int> res{INF, 0};
        if (l < ls) res = q_min(L, l, r < ls ? r : ls, f ^ P[L].rev);
        if (l <= ls && ls < r && P[t].val < res.first) res = {P[t].val, t};
        if (r > ls + 1) {
            auto s = q_min(R, l - ls - 1 > 0 ? l - ls - 1 : 0, r - ls - 1, f ^ P[R].rev);
            if (s.first < res.first) res = s;
        }
        return res;
    }
    std::pair<ll,int> q_max(int t, int l, int r, int f) const {
        if (l <= 0 && r >= P[t].sz) return {P[t].mxv, P[t].mxi};
        const int L = f ? P[t].r : P[t].l, R = f ? P[t].l : P[t].r;
        const int ls = P[L].sz;
        std::pair<ll,int> res{NINF, 0};
        if (l < ls) res = q_max(L, l, r < ls ? r : ls, f ^ P[L].rev);
        if (l <= ls && ls < r && P[t].val > res.first) res = {P[t].val, t};
        if (r > ls + 1) {
            auto s = q_max(R, l - ls - 1 > 0 ? l - ls - 1 : 0, r - ls - 1, f ^ P[R].rev);
            if (s.first > res.first) res = s;
        }
        return res;
    }
    // 区間の要約。空区間は {len=0, sum=0, mn=INF, mx=NINF}
    struct RangeInfo { int len; ll sum, mn, mx; };
    inline void ri_add_node(RangeInfo& w, int t) const {         // 部分木 t をまるごと足す
        w.len += P[t].sz; w.sum += P[t].sum;
        if (P[t].mnv < w.mn) w.mn = P[t].mnv;
        if (P[t].mxv > w.mx) w.mx = P[t].mxv;
    }
    inline void ri_add_val(RangeInfo& w, ll v) const {           // 1 要素だけ足す
        w.len += 1; w.sum += v;
        if (v < w.mn) w.mn = v;
        if (v > w.mx) w.mx = v;
    }
    RangeInfo q_all(int t, int l, int r, int f) const {
        if (l <= 0 && r >= P[t].sz) return {P[t].sz, P[t].sum, P[t].mnv, P[t].mxv};
        const int L = f ? P[t].r : P[t].l, R = f ? P[t].l : P[t].r;
        const int ls = P[L].sz;
        RangeInfo res{0, 0, INF, NINF};
        if (l < ls) res = q_all(L, l, r < ls ? r : ls, f ^ P[L].rev);
        if (l <= ls && ls < r) ri_add_val(res, P[t].val);
        if (r > ls + 1) {
            RangeInfo s = q_all(R, l - ls - 1 > 0 ? l - ls - 1 : 0, r - ls - 1, f ^ P[R].rev);
            res.len += s.len; res.sum += s.sum;
            if (s.mn < res.mn) res.mn = s.mn;
            if (s.mx > res.mx) res.mx = s.mx;
        }
        return res;
    }
    inline bool clip(int& l, int& r) const {        // 区間を [0,N) に丸める。空なら false
        if (l < 0) l = 0;
        const int n = P[proot].sz;
        if (r > n) r = n;
        return l < r;
    }

    // ================================================================
    //  ここから公開 API
    // ================================================================

    // ---------------- 基本 ----------------
    int  size()  const { return P[proot].sz; }
    bool empty() const { return P[proot].sz == 0; }
    ll   total_sum() const { return P[proot].sum; }              // 全体の総和 O(1)

    ll get(int i) const {                                        // 位置 i の値
#ifdef DUS_DEBUG
        assert(0 <= i && i < size());
#endif
        return P[node_at(i)].val;
    }
    ll operator[](int i) const { return get(i); }
    ll front() const { return get(0); }
    ll back()  const { return get(size() - 1); }

    // 部分木 root を中順に走査して out に値を並べる (遅延フラグは下ろさない)
    void inorder_values(int root, std::vector<ll>& out) const {
        std::vector<std::pair<int,int>> st; st.reserve(64);
        int t = root, f = P[root].rev;
        while (t || !st.empty()) {
            while (t) {
                st.push_back({t, f});
                const int L = f ? P[t].r : P[t].l;
                f ^= P[L].rev; t = L;
            }
            auto [u, fu] = st.back(); st.pop_back();
            out.push_back(P[u].val);
            const int R = fu ? P[u].l : P[u].r;
            t = R; f = fu ^ P[R].rev;
        }
    }
    std::vector<ll> to_vector() const {                          // 配列に書き出す O(N)
        std::vector<ll> out; out.reserve(size());
        inorder_values(proot, out);
        return out;
    }

    void build(const std::vector<ll>& a) {                       // O(N log N)
        clear();
        const int n = (int)a.size();
        if (!n) return;
        P.reserve(n + 1); V.reserve(n + 1);
        std::vector<int> ids(n);
        for (int i = 0; i < n; i++) ids[i] = new_node(a[i]);
        {   // 位置木: 中順が a の順になる Cartesian tree を O(N) で
            std::vector<int> stk; stk.reserve(n);
            for (int i = 0; i < n; i++) {
                int u = ids[i], last = 0;
                while (!stk.empty() && P[stk.back()].pri < P[u].pri) { last = stk.back(); stk.pop_back(); }
                P[u].l = last; if (last) P[last].par = u;
                if (!stk.empty()) { P[stk.back()].r = u; P[u].par = stk.back(); }
                stk.push_back(u);
            }
            proot = stk[0]; P[proot].par = 0;
            pull_all_p(proot);
        }
        {   // 値木: 値でソートしてから同様に
            std::vector<int> ord = ids;
            std::sort(ord.begin(), ord.end(), [&](int x, int y) { return V[x].key < V[y].key; });
#ifdef DUS_DEBUG
            for (int i = 1; i < n; i++) assert(V[ord[i-1]].key != V[ord[i]].key);
#endif
            std::vector<int> stk; stk.reserve(n);
            for (int i = 0; i < n; i++) {
                int u = ord[i], last = 0;
                while (!stk.empty() && V[stk.back()].pri < V[u].pri) { last = stk.back(); stk.pop_back(); }
                V[u].l = last; if (last) V[last].par = u;
                if (!stk.empty()) { V[stk.back()].r = u; V[u].par = stk.back(); }
                stk.push_back(u);
            }
            vroot = stk[0]; V[vroot].par = 0;
            pull_all_v(vroot);
        }
        hinit((size_t)n * 2 + 8);
        for (int i = 0; i < n; i++) hput(a[i], ids[i]);
    }
    void pull_all_p(int t) { if (!t) return; pull_all_p(P[t].l); pull_all_p(P[t].r); pull(t); }
    void pull_all_v(int t) { if (!t) return; pull_all_v(V[t].l); pull_all_v(V[t].r); pull_v(t); }

    // ---------------- 挿入・削除・変更 ----------------
    void insert(int i, ll x) {                                   // 位置 i に値 x を挿入
#ifdef DUS_DEBUG
        assert(0 <= i && i <= size());
        assert(!contains(x) && "DynUniqueSeq: 値が重複しています");
#endif
        int u = new_node(x);
        int a, b;
        split_root(proot, i, a, b);
        proot = merge(merge(a, u), b);
        P[proot].par = 0;
        v_insert(u);
        hput(x, u);
    }
    void push_back(ll x)  { insert(size(), x); }
    void push_front(ll x) { insert(0, x); }

    void erase(int i) {                                          // 位置 i を削除
#ifdef DUS_DEBUG
        assert(0 <= i && i < size());
#endif
        int u = node_at(i);
        if (u) erase_node(u);
    }
    void pop_back()  { erase(size() - 1); }
    void pop_front() { erase(0); }

    bool erase_value(ll x) {                                     // 値 x を削除 (無ければ false)
        int u = hget(x);
        if (!u) return false;
        erase_node(u);                   // 位置を求める必要すら無い
        return true;
    }

    void assign_node_value(int u, ll x) {                        // 内部: ノードの値を x に
        ll old = P[u].val;
        if (old == x) return;
        v_erase(u); hdel(old);
        P[u].val = x; V[u].key = x; v_isolate(u);
        v_insert(u); hput(x, u);
        update_up(u);                                            // 木の形は変わらない
    }
    void set(int i, ll x) {                                      // 位置 i の値を x にする
#ifdef DUS_DEBUG
        assert(0 <= i && i < size());
        assert(get(i) == x || !contains(x));
#endif
        assign_node_value(node_at(i), x);
    }
    bool set_value(ll oldv, ll newv) {                           // 値 oldv を newv に変える
        int u = hget(oldv);
        if (!u) return false;
#ifdef DUS_DEBUG
        assert(oldv == newv || !contains(newv));
#endif
        assign_node_value(u, newv);
        return true;
    }

    // ---------------- 値 -> 位置 ----------------
    bool contains(ll x) const { return hget(x) != 0; }
    int  count(ll x)    const { return hget(x) ? 1 : 0; }
    int  index_of(ll x) const { int u = hget(x); return u ? index_of_node(u) : -1; }
    bool contains_in_range(int l, int r, ll x) const {
        int i = index_of(x); return i >= 0 && l <= i && i < r;
    }
    int range_index_of(int l, int r, ll x) const {
        int i = index_of(x); return (i >= 0 && l <= i && i < r) ? i : -1;
    }

    // ---------------- 区間集約 ----------------
    ll range_sum(int l, int r) const { return clip(l, r) ? q_sum(proot, l, r, P[proot].rev) : 0; }
    ll range_min(int l, int r) const { return clip(l, r) ? q_min(proot, l, r, P[proot].rev).first : INF; }
    ll range_max(int l, int r) const { return clip(l, r) ? q_max(proot, l, r, P[proot].rev).first : NINF; }

    std::pair<ll,int> range_min_with_index(int l, int r) const {
        if (!clip(l, r)) return {INF, -1};
        auto [v, u] = q_min(proot, l, r, P[proot].rev);
        return {v, index_of_node(u)};
    }
    std::pair<ll,int> range_max_with_index(int l, int r) const {
        if (!clip(l, r)) return {NINF, -1};
        auto [v, u] = q_max(proot, l, r, P[proot].rev);
        return {v, index_of_node(u)};
    }
    int range_min_index(int l, int r) const { return range_min_with_index(l, r).second; }
    int range_max_index(int l, int r) const { return range_max_with_index(l, r).second; }

    // 区間 2 番目に小さい値 (要素 1 個以下なら INF)。最小値の位置で区間を割るだけ
    ll range_second_min(int l, int r) const {
        if (!clip(l, r) || r - l < 2) return INF;
        auto [v, p] = range_min_with_index(l, r);
        ll a = (p > l)     ? range_min(l, p)     : INF;
        ll b = (p + 1 < r) ? range_min(p + 1, r) : INF;
        return a < b ? a : b;
    }
    ll range_second_max(int l, int r) const {
        if (!clip(l, r) || r - l < 2) return NINF;
        auto [v, p] = range_max_with_index(l, r);
        ll a = (p > l)     ? range_max(l, p)     : NINF;
        ll b = (p + 1 < r) ? range_max(p + 1, r) : NINF;
        return a > b ? a : b;
    }
    std::pair<ll,int> range_second_min_with_index(int l, int r) const {
        ll v = range_second_min(l, r);
        return {v, v == INF ? -1 : index_of(v)};
    }
    std::pair<ll,int> range_second_max_with_index(int l, int r) const {
        ll v = range_second_max(l, r);
        return {v, v == NINF ? -1 : index_of(v)};
    }
    // sum / min / max を「木を 1 回降りるだけ」でまとめて取る (別々に呼ぶより約 3 倍速い)
    // 位置が欲しいときは値が一意なので index_of(値) で O(log N) で引ける
    RangeInfo range_info(int l, int r) const {
        if (!clip(l, r)) return {0, 0, INF, NINF};
        return q_all(proot, l, r, P[proot].rev);
    }
    ll range_max_min_diff(int l, int r) const {                  // max - min
        if (!clip(l, r)) return 0;
        RangeInfo w = q_all(proot, l, r, P[proot].rev);
        return w.mx - w.mn;
    }

    ll get_min() const { return P[proot].mnv; }                  // 全体最小 O(1)
    ll get_max() const { return P[proot].mxv; }                  // 全体最大 O(1)
    std::pair<ll,int> get_min_with_index() const {
        if (empty()) return {INF, -1};
        return {P[proot].mnv, index_of_node(P[proot].mni)};
    }
    std::pair<ll,int> get_max_with_index() const {
        if (empty()) return {NINF, -1};
        return {P[proot].mxv, index_of_node(P[proot].mxi)};
    }
    int min_index() const { return get_min_with_index().second; }
    int max_index() const { return get_max_with_index().second; }

    // ---------------- 値の順序 (set / 順序統計木 相当) ----------------
    int rank_of_value(ll x) const {                              // x 未満の要素数 (= x の順位)
        int c = 0, t = vroot;
        while (t) {
            if (V[t].key < x) { c += V[V[t].l].sz + 1; t = V[t].r; }
            else t = V[t].l;
        }
        return c;
    }
    int count_less(ll x)    const { return rank_of_value(x); }
    int count_le(ll x)      const { return rank_of_value(x) + (contains(x) ? 1 : 0); }
    int count_greater(ll x) const { return size() - count_le(x); }
    int count_ge(ll x)      const { return size() - count_less(x); }
    int count_value_range(ll lo, ll hi) const { return count_le(hi) - count_less(lo); }  // [lo,hi]

    ll sum_less(ll x) const {                                    // 値が x 未満の要素の総和
        ll s = 0; int t = vroot;
        while (t) {
            if (V[t].key < x) { s += V[V[t].l].vsum + V[t].key; t = V[t].r; }
            else t = V[t].l;
        }
        return s;
    }
    ll sum_le(ll x) const { return sum_less(x) + (contains(x) ? x : 0); }
    ll sum_value_range(ll lo, ll hi) const { return sum_le(hi) - sum_less(lo); }          // [lo,hi]
    ll sum_of_k_smallest(int k) const {                          // 小さい方から k 個の総和
        if (k <= 0) return 0;
        if (k >= size()) return total_sum();
        ll s = 0; int t = vroot;
        while (t && k > 0) {
            int ls = V[V[t].l].sz;
            if (k <= ls) t = V[t].l;
            else { s += V[V[t].l].vsum + V[t].key; k -= ls + 1; t = V[t].r; }
        }
        return s;
    }
    ll sum_of_k_largest(int k) const {
        if (k <= 0) return 0;
        if (k >= size()) return total_sum();
        return total_sum() - sum_of_k_smallest(size() - k);
    }

    int v_kth_node(int k) const {                                // 内部: k 番目に小さいノード
        if (k < 0 || k >= V[vroot].sz) return 0;
        int t = vroot;
        while (true) {
            int ls = V[V[t].l].sz;
            if (k < ls) t = V[t].l;
            else if (k == ls) return t;
            else { k -= ls + 1; t = V[t].r; }
        }
    }
    ll kth_smallest(int k) const { int u = v_kth_node(k); return u ? V[u].key : INF; }
    ll kth_largest(int k)  const { int u = v_kth_node(size() - 1 - k); return u ? V[u].key : NINF; }
    std::pair<ll,int> kth_smallest_with_index(int k) const {
        int u = v_kth_node(k); return u ? std::pair<ll,int>{V[u].key, index_of_node(u)} : std::pair<ll,int>{INF, -1};
    }
    std::pair<ll,int> kth_largest_with_index(int k) const {
        int u = v_kth_node(size() - 1 - k); return u ? std::pair<ll,int>{V[u].key, index_of_node(u)} : std::pair<ll,int>{NINF, -1};
    }
    int kth_smallest_index(int k) const { int u = v_kth_node(k); return u ? index_of_node(u) : -1; }
    int kth_largest_index(int k)  const { int u = v_kth_node(size() - 1 - k); return u ? index_of_node(u) : -1; }

    int rank_of_index(int i) const {                             // 位置 i の値は何番目に小さいか
        int u = node_at(i);
        if (!u) return -1;
        int r = V[V[u].l].sz;
        for (int c = u, p = V[u].par; p; c = p, p = V[p].par)
            if (c == V[p].r) r += V[V[p].l].sz + 1;
        return r;
    }

    int v_prev_node(ll x) const { int t = vroot, res = 0; while (t) { if (V[t].key <  x) { res = t; t = V[t].r; } else t = V[t].l; } return res; }
    int v_next_node(ll x) const { int t = vroot, res = 0; while (t) { if (V[t].key >  x) { res = t; t = V[t].l; } else t = V[t].r; } return res; }
    int v_lb_node(ll x)   const { int t = vroot, res = 0; while (t) { if (V[t].key >= x) { res = t; t = V[t].l; } else t = V[t].r; } return res; }

    ll prev_value(ll x)        const { int u = v_prev_node(x); return u ? V[u].key : NINF; }  // x 未満の最大値
    ll next_value(ll x)        const { int u = v_next_node(x); return u ? V[u].key : INF;  }  // x より大きい最小値
    ll lower_bound_value(ll x) const { int u = v_lb_node(x);   return u ? V[u].key : INF;  }  // x 以上の最小値
    ll upper_bound_value(ll x) const { return next_value(x); }                                // x より大きい最小値
    std::pair<ll,int> prev_value_with_index(ll x) const { int u = v_prev_node(x); return u ? std::pair<ll,int>{V[u].key, index_of_node(u)} : std::pair<ll,int>{NINF, -1}; }
    std::pair<ll,int> next_value_with_index(ll x) const { int u = v_next_node(x); return u ? std::pair<ll,int>{V[u].key, index_of_node(u)} : std::pair<ll,int>{INF, -1}; }
    std::pair<ll,int> lower_bound_value_with_index(ll x) const { int u = v_lb_node(x); return u ? std::pair<ll,int>{V[u].key, index_of_node(u)} : std::pair<ll,int>{INF, -1}; }
    std::pair<ll,int> upper_bound_value_with_index(ll x) const { return next_value_with_index(x); }
    int prev_value_index(ll x) const { int u = v_prev_node(x); return u ? index_of_node(u) : -1; }
    int next_value_index(ll x) const { int u = v_next_node(x); return u ? index_of_node(u) : -1; }

    // 位置 i の値の「値として 1 つ下/上」の要素の位置
    int prev_index_by_value(int i) const { int u = node_at(i); return u ? prev_value_index(P[u].val) : -1; }
    int next_index_by_value(int i) const { int u = node_at(i); return u ? next_value_index(P[u].val) : -1; }

    // 値を昇順に並べたとき隣り合う 2 値の差の最小値 (要素 2 個未満なら INF)
    ll minimum_gap() const { return V[vroot].sz < 2 ? INF : (ll)V[vroot].gap; }
    std::pair<ll,ll> minimum_gap_pair() const {                  // {小さい方の値, 大きい方の値}
        if (V[vroot].sz < 2) return {INF, INF};
        const ull g = V[vroot].gap;
        int t = vroot;
        while (true) {
            const int l = V[t].l, r = V[t].r;
            if (l && V[l].gap == g) { t = l; continue; }
            if (l && (ull)((ull)V[t].key - (ull)V[l].maxv) == g) return {V[l].maxv, V[t].key};
            if (r && (ull)((ull)V[r].minv - (ull)V[t].key) == g) return {V[t].key, V[r].minv};
            t = r;
        }
    }
    std::pair<int,int> minimum_gap_indices() const {             // 上の 2 値の「位置」
        auto [a, b] = minimum_gap_pair();
        if (a == INF) return {-1, -1};
        return {index_of(a), index_of(b)};
    }

    // ---------------- 並べ替え ----------------
    void swap_node_values(int u, int v) {                        // 内部: 2 ノードの値を交換
        ll a = P[u].val, b = P[v].val;
        v_erase(u); v_erase(v);
        P[u].val = b; V[u].key = b; v_isolate(u);
        P[v].val = a; V[v].key = a; v_isolate(v);
        v_insert(u); v_insert(v);
        hput(b, u); hput(a, v);
        update_up(u); update_up(v);
    }
    void swap_at(int i, int j) {                                 // 位置 i と j の要素を交換
        if (i == j) return;
        int u = node_at(i), v = node_at(j);
        if (!u || !v || u == v) return;
        swap_node_values(u, v);
    }
    // 値 x と値 y の「位置」を交換する。位置を求めずに済むので swap_at より速い
    bool swap_values(ll x, ll y) {
        int u = hget(x), v = hget(y);
        if (!u || !v) return false;
        if (u != v) swap_node_values(u, v);
        return true;
    }
    // 位置 from の要素を抜き取り、残った配列の位置 to に挿し込む
    void move_element(int from, int to) {
        if (from < 0 || from >= size()) return;
        const int n2 = size() - 1;
        if (to < 0) to = 0;
        if (to > n2) to = n2;
        if (from == to) return;
        int a, b, c;
        split_root(proot, from, a, b);
        split_root(b, 1, b, c);
        int rest = merge(a, c); if (rest) P[rest].par = 0;
        int x, y;
        split_root(rest, to, x, y);
        proot = merge(merge(x, b), y);
        P[proot].par = 0;
    }
    void reverse_range(int l, int r) {                           // [l,r) を反転
        if (!clip(l, r) || r - l < 2) return;
        has_rev = true;
        int a, b, c;
        split_root(proot, l, a, b);
        split_root(b, r - l, b, c);
        P[b].rev ^= 1;
        proot = merge(merge(a, b), c);
        P[proot].par = 0;
    }
    void rotate_range(int l, int m, int r) {                     // [l,r) を m が先頭に来るよう回転
        if (l < 0) l = 0;
        if (r > size()) r = size();
        if (l >= m || m >= r) return;
        int a, b, c, d;
        split_root(proot, l, a, b);
        split_root(b, m - l, b, c);
        split_root(c, r - m, c, d);
        proot = merge(merge(merge(a, c), b), d);
        P[proot].par = 0;
    }
    // 区間 [l,r) を抜き取り、残った配列の位置 pos に挿し込む
    void move_range(int l, int r, int pos) {
        if (l < 0) l = 0;
        if (r > size()) r = size();
        if (l >= r) return;
        const int n2 = size() - (r - l);
        if (pos < 0) pos = 0;
        if (pos > n2) pos = n2;
        int a, b, c;
        split_root(proot, l, a, b);
        split_root(b, r - l, b, c);
        int rest = merge(a, c); if (rest) P[rest].par = 0;
        int x, y;
        split_root(rest, pos, x, y);
        proot = merge(merge(x, b), y);
        P[proot].par = 0;
    }

    // ---------------- まとめて挿入・削除・取り出し ----------------
    void free_subtree(int t) {                                   // 内部: 部分木の全ノードを解放
        if (!t) return;
        free_subtree(P[t].l); free_subtree(P[t].r);
        v_erase(t); hdel(P[t].val);
        P[t].l = P[t].r = P[t].par = 0;
        free_node(t);
    }
    // 区間 [l,r) をまとめて削除。O(log N + (r-l) log N)
    // 1 個ずつ erase するより位置木の操作が 1 回で済む分ずっと速い
    void erase_range(int l, int r) {
        if (l < 0) l = 0;
        if (r > size()) r = size();
        if (l >= r) return;
        int a, b, c;
        split_root(proot, l, a, b);
        split_root(b, r - l, b, c);
        free_subtree(b);
        proot = merge(a, c);
        if (proot) P[proot].par = 0;
    }
    // 位置 pos に配列 xs をまとめて挿入。O(k log N) (位置木の構築は O(k))
    void insert_range(int pos, const std::vector<ll>& xs) {
        const int k = (int)xs.size();
        if (k == 0) return;
        if (pos < 0) pos = 0;
        if (pos > size()) pos = size();
        std::vector<int> ids(k);
        for (int i = 0; i < k; i++) ids[i] = new_node(xs[i]);
        std::vector<int> stk; stk.reserve(k);
        for (int i = 0; i < k; i++) {
            int u = ids[i], last = 0;
            while (!stk.empty() && P[stk.back()].pri < P[u].pri) { last = stk.back(); stk.pop_back(); }
            P[u].l = last; if (last) P[last].par = u;
            if (!stk.empty()) { P[stk.back()].r = u; P[u].par = stk.back(); }
            stk.push_back(u);
        }
        int mid = stk[0]; P[mid].par = 0;
        pull_all_p(mid);
        for (int i = 0; i < k; i++) { v_insert(ids[i]); hput(xs[i], ids[i]); }
        int a, b;
        split_root(proot, pos, a, b);
        proot = merge(merge(a, mid), b);
        P[proot].par = 0;
    }
    // 区間 [l,r) の値を配列で取り出す。O(log N + (r-l))
    std::vector<ll> range_to_vector(int l, int r) {
        std::vector<ll> res;
        if (!clip(l, r)) return res;
        int a, b, c;
        split_root(proot, l, a, b);
        split_root(b, r - l, b, c);
        res.reserve(r - l);
        inorder_values(b, res);
        proot = merge(merge(a, b), c);
        P[proot].par = 0;
        return res;
    }

    // ---------------- 値の昇順 <-> 位置 の一括取得 ----------------
    // 値を昇順に並べたときの「元の位置」の列。O(N)  (機能候補 #34 value_order)
    std::vector<int> value_order() const {
        std::vector<int> res; res.reserve(size());
        if (!proot) return res;
        std::vector<int> pos(P.size(), 0);
        {   // 位置木を中順に走査して各ノードの位置を記録
            std::vector<std::pair<int,int>> st; st.reserve(64);
            int t = proot, f = P[proot].rev, idx = 0;
            while (t || !st.empty()) {
                while (t) { st.push_back({t, f}); const int L = f ? P[t].r : P[t].l; f ^= P[L].rev; t = L; }
                auto [u, fu] = st.back(); st.pop_back();
                pos[u] = idx++;
                const int R = fu ? P[u].l : P[u].r;
                t = R; f = fu ^ P[R].rev;
            }
        }
        {   // 値木を中順に走査すれば値の昇順になる
            std::vector<int> st; st.reserve(64);
            int t = vroot;
            while (t || !st.empty()) {
                while (t) { st.push_back(t); t = V[t].l; }
                int u = st.back(); st.pop_back();
                res.push_back(pos[u]);
                t = V[u].r;
            }
        }
        return res;
    }
    std::vector<ll> sorted_values() const {                      // 値を昇順に並べた配列。O(N)
        std::vector<ll> res; res.reserve(size());
        std::vector<int> st; st.reserve(64);
        int t = vroot;
        while (t || !st.empty()) {
            while (t) { st.push_back(t); t = V[t].l; }
            int u = st.back(); st.pop_back();
            res.push_back(V[u].key);
            t = V[u].r;
        }
        return res;
    }

    // ---------------- 区間の k 番目 (k が小さいとき専用) ----------------
    //  機能候補 #22/#23 の実用版。完全な range_kth_smallest は O(log^2〜log^3 N) の
    //  入れ子構造が必要で重いが、「k 個取り出す」だけなら
    //      区間最小を取る -> その位置で区間を 2 つに割る -> 優先度付きキューに戻す
    //  を k 回繰り返すだけでよく、O(k (log N + log k))。k が小さい用途ならこちらが速い。
    //  返り値は {値, 位置} の昇順 (k > 区間長 なら区間長で打ち切り)
    std::vector<std::pair<ll,int>> range_k_smallest(int l, int r, int k) const {
        std::vector<std::pair<ll,int>> res;
        if (!clip(l, r) || k <= 0) return res;
        if (k > r - l) k = r - l;
        res.reserve(k);
        struct Item { ll v; int node, lo, hi; };
        auto worse = [](const Item& a, const Item& b) { return a.v > b.v; };
        std::priority_queue<Item, std::vector<Item>, decltype(worse)> pq(worse);
        { auto [v, u] = q_min(proot, l, r, P[proot].rev); pq.push({v, u, l, r}); }
        while ((int)res.size() < k) {
            Item it = pq.top(); pq.pop();
            const int p = index_of_node(it.node);
            res.push_back({it.v, p});
            if (it.lo < p)     { auto [v, u] = q_min(proot, it.lo, p, P[proot].rev); pq.push({v, u, it.lo, p}); }
            if (p + 1 < it.hi) { auto [v, u] = q_min(proot, p + 1, it.hi, P[proot].rev); pq.push({v, u, p + 1, it.hi}); }
        }
        return res;
    }
    std::vector<std::pair<ll,int>> range_k_largest(int l, int r, int k) const {
        std::vector<std::pair<ll,int>> res;
        if (!clip(l, r) || k <= 0) return res;
        if (k > r - l) k = r - l;
        res.reserve(k);
        struct Item { ll v; int node, lo, hi; };
        auto worse = [](const Item& a, const Item& b) { return a.v < b.v; };
        std::priority_queue<Item, std::vector<Item>, decltype(worse)> pq(worse);
        { auto [v, u] = q_max(proot, l, r, P[proot].rev); pq.push({v, u, l, r}); }
        while ((int)res.size() < k) {
            Item it = pq.top(); pq.pop();
            const int p = index_of_node(it.node);
            res.push_back({it.v, p});
            if (it.lo < p)     { auto [v, u] = q_max(proot, it.lo, p, P[proot].rev); pq.push({v, u, it.lo, p}); }
            if (p + 1 < it.hi) { auto [v, u] = q_max(proot, p + 1, it.hi, P[proot].rev); pq.push({v, u, p + 1, it.hi}); }
        }
        return res;
    }
    // 区間 [l,r) で k 番目 (0-indexed) に小さい値。O(k log N) なので k が小さいときだけ使う
    ll range_kth_smallest(int l, int r, int k) const {
        auto v = range_k_smallest(l, r, k + 1);
        return (int)v.size() == k + 1 ? v[k].first : INF;
    }
    ll range_kth_largest(int l, int r, int k) const {
        auto v = range_k_largest(l, r, k + 1);
        return (int)v.size() == k + 1 ? v[k].first : NINF;
    }
    std::pair<ll,int> range_kth_smallest_with_index(int l, int r, int k) const {
        auto v = range_k_smallest(l, r, k + 1);
        return (int)v.size() == k + 1 ? v[k] : std::pair<ll,int>{INF, -1};
    }
    std::pair<ll,int> range_kth_largest_with_index(int l, int r, int k) const {
        auto v = range_k_largest(l, r, k + 1);
        return (int)v.size() == k + 1 ? v[k] : std::pair<ll,int>{NINF, -1};
    }

    // ---------------- mex ----------------
    //  from 以上で「配列に無い」最小の整数。値が相異なるので
    //      「[from, from+x) に x 個ある」 <=> 「from..from+x-1 が全部ある」
    //  が成り立ち、この条件は x について単調なので二分探索できる。O(log^2 N)
    ll mex(ll from = 0) const {
        const int base = rank_of_value(from);
        int lo = 0, hi = size() - base;
        while (lo < hi) {
            const int mid = lo + (hi - lo + 1) / 2;
            if (rank_of_value(from + mid) - base == mid) lo = mid; else hi = mid - 1;
        }
        return from + lo;
    }

    // ---------------- 区間二分探索 (ACL の max_right / min_left 相当) ----------------
    //  f は RangeInfo を受けて bool を返す述語で、「区間を伸ばすと true -> false に
    //  一度だけ変わる」単調性が必要。空区間では true であること。
    //  例)  f = [&](RangeInfo w){ return w.mn >= v; }        … 最小値が v 以上
    //       f = [&](RangeInfo w){ return w.len == 0 || w.mx - w.mn <= d; }
    //  いずれも O(log N)。内部で split/merge するので const ではない。
    template <class F> int mr_desc(int t, RangeInfo& acc, F& f) {
        if (!t) return 0;
        { RangeInfo w = acc; ri_add_node(w, t); if (f(w)) { acc = w; return P[t].sz; } }
        push(t);
        const int L = P[t].l, R = P[t].r;
        const int c = mr_desc(L, acc, f);
        if (c < P[L].sz) return c;
        { RangeInfo w = acc; ri_add_val(w, P[t].val); if (!f(w)) return c; acc = w; }
        return c + 1 + mr_desc(R, acc, f);
    }
    template <class F> int ml_desc(int t, RangeInfo& acc, F& f) {
        if (!t) return 0;
        { RangeInfo w = acc; ri_add_node(w, t); if (f(w)) { acc = w; return P[t].sz; } }
        push(t);
        const int L = P[t].l, R = P[t].r;
        const int c = ml_desc(R, acc, f);
        if (c < P[R].sz) return c;
        { RangeInfo w = acc; ri_add_val(w, P[t].val); if (!f(w)) return c; acc = w; }
        return c + 1 + ml_desc(L, acc, f);
    }
    // f(range_info(l, r)) が true となる最大の r
    template <class F> int max_right(int l, F f) {
        if (l < 0) l = 0;
        const int n = size();
        if (l >= n) return n;
        int a, b;
        split_root(proot, l, a, b);
        RangeInfo acc{0, 0, INF, NINF};
        const int c = mr_desc(b, acc, f);
        proot = merge(a, b);
        if (proot) P[proot].par = 0;
        return l + c;
    }
    // f(range_info(l, r)) が true となる最小の l
    template <class F> int min_left(int r, F f) {
        const int n = size();
        if (r > n) r = n;
        if (r <= 0) return 0;
        int a, b;
        split_root(proot, r, a, b);
        RangeInfo acc{0, 0, INF, NINF};
        const int c = ml_desc(a, acc, f);
        proot = merge(a, b);
        if (proot) P[proot].par = 0;
        return r - c;
    }
};

// ===================== テンプレート本体 ここまで =====================

/* =====================================================================================
 *  ここから下は「検証用コード」。AtCoder に提出するときは不要。
 *      1) Naive        : 愚直実装 (vector を舐めるだけ) ― 正解の基準
 *      2) validate     : 木の内部不変条件 (親ポインタ/heap条件/集約値/ハッシュ表) の全数検査
 *      3) random_test  : 全 API を愚直実装と突き合わせるランダムテスト
 *      4) demo         : 日本語解説付きの動作確認 (出力を読めば使い方を思い出せる)
 *      5) bench        : N = 2*10^5 規模の速度計測
 * ===================================================================================== */

using ll  = long long;
using ull = unsigned long long;
static constexpr ll TINF  = std::numeric_limits<ll>::max();
static constexpr ll TNINF = std::numeric_limits<ll>::min();

// ---------------------------------------------------------------------------- 1) 愚直実装
struct Naive {
    std::vector<ll> a;
    int size() const { return (int)a.size(); }
    std::vector<ll> sorted() const { std::vector<ll> b = a; std::sort(b.begin(), b.end()); return b; }

    int index_of(ll x) const {
        for (int i = 0; i < size(); i++) if (a[i] == x) return i;
        return -1;
    }
    ll range_sum(int l, int r) const { ll s = 0; for (int i = l; i < r; i++) s += a[i]; return s; }
    ll range_min(int l, int r) const { ll v = TINF;  for (int i = l; i < r; i++) v = std::min(v, a[i]); return v; }
    ll range_max(int l, int r) const { ll v = TNINF; for (int i = l; i < r; i++) v = std::max(v, a[i]); return v; }
    int range_min_index(int l, int r) const {
        int p = -1; for (int i = l; i < r; i++) if (p < 0 || a[i] < a[p]) p = i; return p;
    }
    int range_max_index(int l, int r) const {
        int p = -1; for (int i = l; i < r; i++) if (p < 0 || a[i] > a[p]) p = i; return p;
    }
    ll range_second_min(int l, int r) const {
        if (r - l < 2) return TINF;
        std::vector<ll> b(a.begin() + l, a.begin() + r); std::sort(b.begin(), b.end()); return b[1];
    }
    ll range_second_max(int l, int r) const {
        if (r - l < 2) return TNINF;
        std::vector<ll> b(a.begin() + l, a.begin() + r); std::sort(b.rbegin(), b.rend()); return b[1];
    }
    int rank_of_value(ll x) const { int c = 0; for (ll v : a) if (v < x) c++; return c; }
    ll  sum_less(ll x) const { ll s = 0; for (ll v : a) if (v < x) s += v; return s; }
    ll  kth_smallest(int k) const { if (k < 0 || k >= size()) return TINF;  auto b = sorted(); return b[k]; }
    ll  kth_largest(int k)  const { if (k < 0 || k >= size()) return TNINF; auto b = sorted(); return b[size() - 1 - k]; }
    ll  sum_of_k_smallest(int k) const {
        if (k <= 0) return 0;
        auto b = sorted(); ll s = 0;
        for (int i = 0; i < std::min(k, size()); i++) s += b[i];
        return s;
    }
    ll  sum_of_k_largest(int k) const {
        if (k <= 0) return 0;
        auto b = sorted(); ll s = 0;
        for (int i = 0; i < std::min(k, size()); i++) s += b[size() - 1 - i];
        return s;
    }
    int rank_of_index(int i) const { return rank_of_value(a[i]); }
    ll prev_value(ll x) const { ll v = TNINF; for (ll t : a) if (t <  x) v = std::max(v, t); return v; }
    ll next_value(ll x) const { ll v = TINF;  for (ll t : a) if (t >  x) v = std::min(v, t); return v; }
    ll lower_bound_value(ll x) const { ll v = TINF; for (ll t : a) if (t >= x) v = std::min(v, t); return v; }
    std::vector<int> value_order() const {
        std::vector<int> id(size()); std::iota(id.begin(), id.end(), 0);
        std::sort(id.begin(), id.end(), [&](int x, int y) { return a[x] < a[y]; });
        return id;
    }
    std::vector<std::pair<ll,int>> range_k_smallest(int l, int r, int k) const {
        std::vector<std::pair<ll,int>> v;
        if (k <= 0) return v;
        for (int i = l; i < r; i++) v.push_back({a[i], i});
        std::sort(v.begin(), v.end());
        if ((int)v.size() > k) v.resize(k);
        return v;
    }
    std::vector<std::pair<ll,int>> range_k_largest(int l, int r, int k) const {
        std::vector<std::pair<ll,int>> v;
        if (k <= 0) return v;
        for (int i = l; i < r; i++) v.push_back({a[i], i});
        std::sort(v.begin(), v.end(), [](const std::pair<ll,int>& x, const std::pair<ll,int>& y) { return x.first > y.first; });
        if ((int)v.size() > k) v.resize(k);
        return v;
    }
    ll mex(ll from) const {
        std::set<ll> s(a.begin(), a.end());
        ll x = from;
        while (s.count(x)) x++;
        return x;
    }
    // f は (len, mn, mx) を受け取る述語。愚直に全部試す
    template <class F> int max_right(int l, F f) const {
        int best = l;
        for (int r = l; r <= size(); r++) {
            int len = r - l;
            ll mn = TINF, mx = TNINF;
            for (int i = l; i < r; i++) { mn = std::min(mn, a[i]); mx = std::max(mx, a[i]); }
            if (f(len, mn, mx)) best = r; else break;
        }
        return best;
    }
    template <class F> int min_left(int r, F f) const {
        int best = r;
        for (int l = r; l >= 0; l--) {
            int len = r - l;
            ll mn = TINF, mx = TNINF;
            for (int i = l; i < r; i++) { mn = std::min(mn, a[i]); mx = std::max(mx, a[i]); }
            if (f(len, mn, mx)) best = l; else break;
        }
        return best;
    }
    ll minimum_gap() const {
        if (size() < 2) return TINF;
        auto b = sorted(); ull g = ~0ULL;
        for (int i = 1; i < size(); i++) g = std::min(g, (ull)b[i] - (ull)b[i - 1]);
        return (ll)g;
    }
};

// pair / vector を cout に出すための補助 (テスト用)
static std::ostream& operator<<(std::ostream& os, const std::pair<ll,int>& p) {
    return os << "{" << p.first << "," << p.second << "}";
}
static std::ostream& operator<<(std::ostream& os, const std::vector<int>& v) {
    os << "[";
    for (size_t i = 0; i < v.size(); i++) { if (i) os << ","; os << v[i]; }
    return os << "]";
}
static std::ostream& operator<<(std::ostream& os, const std::vector<std::pair<ll,int>>& v) {
    os << "[";
    for (size_t i = 0; i < v.size(); i++) { if (i) os << ","; os << "{" << v[i].first << "," << v[i].second << "}"; }
    return os << "]";
}
static std::ostream& operator<<(std::ostream& os, const std::vector<ll>& v) {
    os << "[";
    for (size_t i = 0; i < v.size(); i++) { if (i) os << ","; os << v[i]; }
    return os << "]";
}

// ---------------------------------------------------------------------------- 2) 内部検査
struct Validator {
    const DynUniqueSeq& s;
    bool ok = true;
    std::string msg;
    explicit Validator(const DynUniqueSeq& s_) : s(s_) {}
    void fail(const std::string& m) { if (ok) { ok = false; msg = m; } }

    struct Agg { int sz; ll sum, mn, mx; };
    Agg checkP(int t, int f, std::vector<ll>& inorder) {
        if (!t) return {0, 0, TINF, TNINF};
        const int pl = s.P[t].l, pr = s.P[t].r;
        if (pl && s.P[pl].par != t) fail("位置木: 親ポインタ不整合(左)");
        if (pr && s.P[pr].par != t) fail("位置木: 親ポインタ不整合(右)");
        if (pl && s.P[pl].pri > s.P[t].pri) fail("位置木: heap 条件違反(左)");
        if (pr && s.P[pr].pri > s.P[t].pri) fail("位置木: heap 条件違反(右)");
        const int L = f ? pr : pl, R = f ? pl : pr;
        Agg x = checkP(L, f ^ s.P[L].rev, inorder);
        inorder.push_back(s.P[t].val);
        Agg y = checkP(R, f ^ s.P[R].rev, inorder);
        Agg c;
        c.sz  = x.sz + y.sz + 1;
        c.sum = x.sum + s.P[t].val + y.sum;
        c.mn  = std::min(s.P[t].val, std::min(x.mn, y.mn));
        c.mx  = std::max(s.P[t].val, std::max(x.mx, y.mx));
        if (s.P[t].sz  != c.sz)  fail("位置木: sz 不整合");
        if (s.P[t].sum != c.sum) fail("位置木: sum 不整合");
        if (s.P[t].mnv != c.mn)  fail("位置木: min 値 不整合");
        if (s.P[t].mxv != c.mx)  fail("位置木: max 値 不整合");
        if (s.P[s.P[t].mni].val != c.mn) fail("位置木: min のノード番号 不整合");
        if (s.P[s.P[t].mxi].val != c.mx) fail("位置木: max のノード番号 不整合");
        return c;
    }
    std::vector<ll> checkV(int t) {
        if (!t) return {};
        const int vl = s.V[t].l, vr = s.V[t].r;
        if (vl && s.V[vl].par != t) fail("値木: 親ポインタ不整合(左)");
        if (vr && s.V[vr].par != t) fail("値木: 親ポインタ不整合(右)");
        if (vl && s.V[vl].pri > s.V[t].pri) fail("値木: heap 条件違反(左)");
        if (vr && s.V[vr].pri > s.V[t].pri) fail("値木: heap 条件違反(右)");
        if (s.V[t].key != s.P[t].val) fail("値木: key と 位置木の val が不一致");
        std::vector<ll> x = checkV(vl), y = checkV(vr);
        std::vector<ll> c = x;
        c.push_back(s.V[t].key);
        for (ll v : y) c.push_back(v);
        for (size_t i = 1; i < c.size(); i++) if (c[i - 1] >= c[i]) fail("値木: 昇順になっていない");
        ll sum = 0; for (ll v : c) sum += v;
        ull g = ~0ULL;
        for (size_t i = 1; i < c.size(); i++) g = std::min(g, (ull)c[i] - (ull)c[i - 1]);
        if (s.V[t].sz   != (int)c.size()) fail("値木: sz 不整合");
        if (s.V[t].vsum != sum)           fail("値木: 総和 不整合");
        if (s.V[t].minv != c.front())     fail("値木: minv 不整合");
        if (s.V[t].maxv != c.back())      fail("値木: maxv 不整合");
        if (s.V[t].gap  != g)             fail("値木: 最小ギャップ 不整合");
        return c;
    }
    bool run(const std::vector<ll>& expect, std::string& why) {
        if (s.size() != (int)expect.size()) fail("size() が期待値と違う");
        if (s.proot && s.P[s.proot].par != 0) fail("位置木: 根の親が 0 でない");
        if (s.vroot && s.V[s.vroot].par != 0) fail("値木: 根の親が 0 でない");
        std::vector<ll> inorder;
        checkP(s.proot, s.P[s.proot].rev, inorder);
        if (inorder != expect) fail("位置木の中順走査が期待する配列と違う");
        std::vector<ll> vs = checkV(s.vroot);
        std::vector<ll> es = expect; std::sort(es.begin(), es.end());
        if (vs != es) fail("値木の中順走査がソート済み配列と違う");
        if (s.hcnt != expect.size()) fail("ハッシュ表の要素数が違う");
        for (size_t i = 0; i < expect.size(); i++) {
            int u = s.hget(expect[i]);
            if (!u) fail("ハッシュ表に値が見つからない");
            else if (s.P[u].val != expect[i]) fail("ハッシュ表が別の値のノードを指している");
            else if (s.index_of_node(u) != (int)i) fail("index_of_node が正しい位置を返さない");
        }
        why = msg;
        return ok;
    }
};

static bool validate(const DynUniqueSeq& s, const std::vector<ll>& expect, std::string& why) {
    Validator v(s);
    return v.run(expect, why);
}

// ---------------------------------------------------------------------------- 3) ランダムテスト
struct Tester {
    std::mt19937_64 rng;
    long long ng = 0;         // 失敗数
    int checks = 0;           // 比較した項目数
    std::string ctx;

    explicit Tester(ull seed) : rng(seed) {}
    ll rnd(ll lo, ll hi) { return lo + (ll)(rng() % (ull)(hi - lo + 1)); }

    template <class T>
    void eq(const T& got, const T& want, const char* name) {
        checks++;
        if (!(got == want)) {
            ng++;
            if (ng <= 5) {
                std::cout << "  [NG] " << name << " : got=" << got << " want=" << want
                          << "   (" << ctx << ")\n";
            }
        }
    }

    // --- 全読み取り API を愚直実装と突き合わせる ---
    void check_queries(DynUniqueSeq& s, const Naive& t) {
        const int n = t.size();
        eq(s.size(), n, "size");
        eq(s.empty(), n == 0, "empty");
        eq(s.to_vector(), t.a, "to_vector");
        eq(s.total_sum(), t.range_sum(0, n), "total_sum");
        eq(s.get_min(), n ? t.range_min(0, n) : TINF, "get_min");
        eq(s.get_max(), n ? t.range_max(0, n) : TNINF, "get_max");
        eq(s.min_index(), n ? t.range_min_index(0, n) : -1, "min_index");
        eq(s.max_index(), n ? t.range_max_index(0, n) : -1, "max_index");
        eq(s.minimum_gap(), t.minimum_gap(), "minimum_gap");
        if (n >= 2) {   // ギャップの組は「差が最小ギャップで、間に値が無い」ことで検証
            auto [p, q] = s.minimum_gap_pair();
            auto [pi, qi] = s.minimum_gap_indices();
            eq(q - p, t.minimum_gap(), "minimum_gap_pair の差");
            eq(t.index_of(p), pi, "minimum_gap_indices[0]");
            eq(t.index_of(q), qi, "minimum_gap_indices[1]");
            eq(t.rank_of_value(q) - t.rank_of_value(p), 1, "minimum_gap_pair が値として隣接");
        }
        eq(s.sorted_values(), t.sorted(), "sorted_values");
        eq(s.value_order(), t.value_order(), "value_order");
        for (int rep = 0; rep < 2; rep++) {        // mex (存在する値の近傍を狙う)
            ll from = (n > 0 && (rng() & 1)) ? t.a[rnd(0, n - 1)] : rnd(vlo, vhi);
            eq(s.mex(from), t.mex(from), "mex");
        }
        eq(s.mex(), t.mex(0), "mex(0)");
        if (n > 0) {
            eq(s.front(), t.a[0], "front");
            eq(s.back(), t.a[n - 1], "back");
            for (int rep = 0; rep < 3; rep++) {
                int i = (int)rnd(0, n - 1);
                eq(s.get(i), t.a[i], "get");
                eq(s[i], t.a[i], "operator[]");
                eq(s.rank_of_index(i), t.rank_of_index(i), "rank_of_index");
                ll pv = t.prev_value(t.a[i]), nv = t.next_value(t.a[i]);
                eq(s.prev_index_by_value(i), pv == TNINF ? -1 : t.index_of(pv), "prev_index_by_value");
                eq(s.next_index_by_value(i), nv == TINF  ? -1 : t.index_of(nv), "next_index_by_value");
            }
            for (int rep = 0; rep < 4; rep++) {     // 区間クエリ
                int l = (int)rnd(0, n - 1), r = (int)rnd(0, n);
                if (l > r) std::swap(l, r);
                if (l == r) r = std::min(n, r + 1);
                eq(s.range_sum(l, r), t.range_sum(l, r), "range_sum");
                eq(s.range_min(l, r), t.range_min(l, r), "range_min");
                eq(s.range_max(l, r), t.range_max(l, r), "range_max");
                eq(s.range_min_index(l, r), t.range_min_index(l, r), "range_min_index");
                eq(s.range_max_index(l, r), t.range_max_index(l, r), "range_max_index");
                eq(s.range_min_with_index(l, r),
                   std::pair<ll,int>{t.range_min(l, r), t.range_min_index(l, r)}, "range_min_with_index");
                eq(s.range_max_with_index(l, r),
                   std::pair<ll,int>{t.range_max(l, r), t.range_max_index(l, r)}, "range_max_with_index");
                eq(s.range_second_min(l, r), t.range_second_min(l, r), "range_second_min");
                eq(s.range_second_max(l, r), t.range_second_max(l, r), "range_second_max");
                ll s2 = t.range_second_min(l, r);
                eq(s.range_second_min_with_index(l, r),
                   std::pair<ll,int>{s2, s2 == TINF ? -1 : t.index_of(s2)}, "range_second_min_with_index");
                ll s3 = t.range_second_max(l, r);
                eq(s.range_second_max_with_index(l, r),
                   std::pair<ll,int>{s3, s3 == TNINF ? -1 : t.index_of(s3)}, "range_second_max_with_index");
                eq(s.range_max_min_diff(l, r), t.range_max(l, r) - t.range_min(l, r), "range_max_min_diff");
                auto ri = s.range_info(l, r);
                eq(ri.len, r - l, "range_info.len");
                eq(ri.sum, t.range_sum(l, r), "range_info.sum");
                eq(ri.mn, t.range_min(l, r), "range_info.mn");
                eq(ri.mx, t.range_max(l, r), "range_info.mx");
                eq(s.range_to_vector(l, r), std::vector<ll>(t.a.begin() + l, t.a.begin() + r), "range_to_vector");
                int kk = (int)rnd(0, std::min(5, r - l + 1));
                eq(s.range_k_smallest(l, r, kk), t.range_k_smallest(l, r, kk), "range_k_smallest");
                eq(s.range_k_largest(l, r, kk), t.range_k_largest(l, r, kk), "range_k_largest");
                {
                    int k2 = (int)rnd(0, r - l);
                    auto ks = t.range_k_smallest(l, r, k2 + 1);
                    auto kl = t.range_k_largest(l, r, k2 + 1);
                    bool in = (int)ks.size() == k2 + 1;
                    eq(s.range_kth_smallest(l, r, k2), in ? ks[k2].first : TINF, "range_kth_smallest");
                    eq(s.range_kth_largest(l, r, k2), in ? kl[k2].first : TNINF, "range_kth_largest");
                    eq(s.range_kth_smallest_with_index(l, r, k2), in ? ks[k2] : std::pair<ll,int>{TINF, -1}, "range_kth_smallest_with_index");
                    eq(s.range_kth_largest_with_index(l, r, k2),  in ? kl[k2] : std::pair<ll,int>{TNINF, -1}, "range_kth_largest_with_index");
                }
            }
            for (int rep = 0; rep < 3; rep++) {     // k 番目
                int k = (int)rnd(-1, n);
                bool in = (0 <= k && k < n);
                eq(s.kth_smallest(k), t.kth_smallest(k), "kth_smallest");
                eq(s.kth_largest(k), t.kth_largest(k), "kth_largest");
                eq(s.kth_smallest_index(k), in ? t.index_of(t.kth_smallest(k)) : -1, "kth_smallest_index");
                eq(s.kth_largest_index(k),  in ? t.index_of(t.kth_largest(k))  : -1, "kth_largest_index");
                eq(s.kth_smallest_with_index(k),
                   std::pair<ll,int>{t.kth_smallest(k), in ? t.index_of(t.kth_smallest(k)) : -1}, "kth_smallest_with_index");
                eq(s.kth_largest_with_index(k),
                   std::pair<ll,int>{t.kth_largest(k), in ? t.index_of(t.kth_largest(k)) : -1}, "kth_largest_with_index");
                eq(s.sum_of_k_smallest(k), t.sum_of_k_smallest(k), "sum_of_k_smallest");
                eq(s.sum_of_k_largest(k),  t.sum_of_k_largest(k),  "sum_of_k_largest");
            }
        }
        for (int rep = 0; rep < 3; rep++) {         // 区間二分探索 (3 種類の単調な述語で)
            const int l = (int)rnd(0, n), r = (int)rnd(0, n);
            const ll v = (n > 0 && (rng() & 1)) ? t.a[rnd(0, n - 1)] : rnd(vlo, vhi);
            const ll d = rnd(0, 30);
            using RI = DynUniqueSeq::RangeInfo;
            auto f1 = [&](RI w) { return w.mn >= v; };
            auto g1 = [&](int, ll mn, ll) { return mn >= v; };
            auto f2 = [&](RI w) { return w.mx <= v; };
            auto g2 = [&](int, ll, ll mx) { return mx <= v; };
            auto f3 = [&](RI w) { return w.len == 0 || w.mx - w.mn <= d; };
            auto g3 = [&](int len, ll mn, ll mx) { return len == 0 || mx - mn <= d; };
            eq(s.max_right(l, f1), t.max_right(l, g1), "max_right(min>=v)");
            eq(s.max_right(l, f2), t.max_right(l, g2), "max_right(max<=v)");
            eq(s.max_right(l, f3), t.max_right(l, g3), "max_right(max-min<=d)");
            eq(s.min_left(r, f1), t.min_left(r, g1), "min_left(min>=v)");
            eq(s.min_left(r, f2), t.min_left(r, g2), "min_left(max<=v)");
            eq(s.min_left(r, f3), t.min_left(r, g3), "min_left(max-min<=d)");
        }
        for (int rep = 0; rep < 5; rep++) {         // 値に関するクエリ (存在する値・しない値の両方)
            ll x;
            if (n > 0 && (rng() & 1)) x = t.a[rnd(0, n - 1)] + (ll)rnd(-1, 1);
            else x = rnd(vlo, vhi);
            eq(s.contains(x), t.index_of(x) >= 0, "contains");
            eq(s.count(x), t.index_of(x) >= 0 ? 1 : 0, "count");
            eq(s.index_of(x), t.index_of(x), "index_of");
            eq(s.rank_of_value(x), t.rank_of_value(x), "rank_of_value");
            eq(s.count_less(x), t.rank_of_value(x), "count_less");
            eq(s.count_le(x), t.rank_of_value(x) + (t.index_of(x) >= 0 ? 1 : 0), "count_le");
            eq(s.count_greater(x), n - (t.rank_of_value(x) + (t.index_of(x) >= 0 ? 1 : 0)), "count_greater");
            eq(s.count_ge(x), n - t.rank_of_value(x), "count_ge");
            eq(s.sum_less(x), t.sum_less(x), "sum_less");
            eq(s.sum_le(x), t.sum_less(x) + (t.index_of(x) >= 0 ? x : 0), "sum_le");
            eq(s.prev_value(x), t.prev_value(x), "prev_value");
            eq(s.next_value(x), t.next_value(x), "next_value");
            eq(s.lower_bound_value(x), t.lower_bound_value(x), "lower_bound_value");
            eq(s.upper_bound_value(x), t.next_value(x), "upper_bound_value");
            ll pv = t.prev_value(x), nv = t.next_value(x), lb = t.lower_bound_value(x);
            eq(s.prev_value_index(x), pv == TNINF ? -1 : t.index_of(pv), "prev_value_index");
            eq(s.next_value_index(x), nv == TINF  ? -1 : t.index_of(nv), "next_value_index");
            eq(s.prev_value_with_index(x), std::pair<ll,int>{pv, pv == TNINF ? -1 : t.index_of(pv)}, "prev_value_with_index");
            eq(s.next_value_with_index(x), std::pair<ll,int>{nv, nv == TINF  ? -1 : t.index_of(nv)}, "next_value_with_index");
            eq(s.lower_bound_value_with_index(x), std::pair<ll,int>{lb, lb == TINF ? -1 : t.index_of(lb)}, "lower_bound_value_with_index");
            ll lo = x, hi = x + rnd(0, 10);
            int cnt = 0; ll sm = 0;
            for (ll v : t.a) if (lo <= v && v <= hi) { cnt++; sm += v; }
            eq(s.count_value_range(lo, hi), cnt, "count_value_range");
            eq(s.sum_value_range(lo, hi), sm, "sum_value_range");
            if (n > 0) {
                int l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n);
                int gi = t.index_of(x);
                bool in = (gi >= 0 && l <= gi && gi < r);
                eq(s.contains_in_range(l, r, x), in, "contains_in_range");
                eq(s.range_index_of(l, r, x), in ? gi : -1, "range_index_of");
            }
        }
    }

    ll vlo = -20, vhi = 20;

    ll fresh_value(const std::set<ll>& used) {   // まだ使われていない値を引く
        for (int i = 0; i < 200; i++) {
            ll x = rnd(vlo, vhi);
            if (!used.count(x)) return x;
        }
        return TINF;   // 引けなかった
    }

    // --- 1 セット分のランダムテスト ---
    void run(int iters, int maxn, ll lo, ll hi, bool do_validate, const char* label) {
        vlo = lo; vhi = hi;
        DynUniqueSeq s;
        Naive t;
        std::set<ll> used;
        // 初期配列を build で作る
        {
            int n0 = (int)rnd(0, std::min(maxn, 10));
            std::vector<ll> init;
            for (int i = 0; i < n0; i++) {
                ll x = fresh_value(used);
                if (x == TINF) break;
                used.insert(x); init.push_back(x);
            }
            s.build(init); t.a = init;
        }
        for (int it = 0; it < iters; it++) {
            const int n = t.size();
            int op = (int)rnd(0, 16);
            if (n == 0) op = 0;
            if (n >= maxn && op == 0) op = 1;
            ctx = std::string(label) + " iter=" + std::to_string(it) + " op=" + std::to_string(op)
                + " n=" + std::to_string(n);
            switch (op) {
            case 0: {   // insert
                ll x = fresh_value(used); if (x == TINF) break;
                int i = (int)rnd(0, n);
                s.insert(i, x); t.a.insert(t.a.begin() + i, x); used.insert(x);
                break;
            }
            case 1: {   // erase
                int i = (int)rnd(0, n - 1);
                used.erase(t.a[i]);
                s.erase(i); t.a.erase(t.a.begin() + i);
                break;
            }
            case 2: {   // erase_value (存在する値 / しない値 の両方を試す)
                ll x = (rng() & 1) ? t.a[rnd(0, n - 1)] : rnd(vlo, vhi);
                int i = t.index_of(x);
                eq(s.erase_value(x), i >= 0, "erase_value の戻り値");
                if (i >= 0) { used.erase(x); t.a.erase(t.a.begin() + i); }
                break;
            }
            case 3: {   // set
                ll x = fresh_value(used); if (x == TINF) break;
                int i = (int)rnd(0, n - 1);
                used.erase(t.a[i]); used.insert(x);
                s.set(i, x); t.a[i] = x;
                break;
            }
            case 4: {   // set_value
                ll nx = fresh_value(used); if (nx == TINF) break;
                ll ox = (rng() & 1) ? t.a[rnd(0, n - 1)] : rnd(vlo, vhi);
                int i = t.index_of(ox);
                eq(s.set_value(ox, nx), i >= 0, "set_value の戻り値");
                if (i >= 0) { used.erase(ox); used.insert(nx); t.a[i] = nx; }
                break;
            }
            case 5: {   // swap_at
                int i = (int)rnd(0, n - 1), j = (int)rnd(0, n - 1);
                s.swap_at(i, j); std::swap(t.a[i], t.a[j]);
                break;
            }
            case 6: {   // move_element
                int f = (int)rnd(0, n - 1), to = (int)rnd(0, n - 1);
                s.move_element(f, to);
                ll v = t.a[f]; t.a.erase(t.a.begin() + f); t.a.insert(t.a.begin() + to, v);
                break;
            }
            case 7: {   // reverse_range
                int l = (int)rnd(0, n - 1), r = (int)rnd(0, n);
                if (l > r) std::swap(l, r);
                s.reverse_range(l, r);
                std::reverse(t.a.begin() + l, t.a.begin() + r);
                break;
            }
            case 8: {   // rotate_range
                if (n < 3) break;
                int l = (int)rnd(0, n - 3), r = (int)rnd(l + 2, n), m = (int)rnd(l + 1, r - 1);
                s.rotate_range(l, m, r);
                std::rotate(t.a.begin() + l, t.a.begin() + m, t.a.begin() + r);
                break;
            }
            case 9: {   // move_range
                int l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n);
                int pos = (int)rnd(0, n - (r - l));
                s.move_range(l, r, pos);
                std::vector<ll> mid(t.a.begin() + l, t.a.begin() + r);
                t.a.erase(t.a.begin() + l, t.a.begin() + r);
                t.a.insert(t.a.begin() + pos, mid.begin(), mid.end());
                break;
            }
            case 10: {  // push_back / push_front
                if (n >= maxn) break;
                ll x = fresh_value(used); if (x == TINF) break;
                used.insert(x);
                if (rng() & 1) { s.push_back(x);  t.a.push_back(x); }
                else           { s.push_front(x); t.a.insert(t.a.begin(), x); }
                break;
            }
            case 11: {  // pop_back / pop_front
                if (rng() & 1) { used.erase(t.a.back());  s.pop_back();  t.a.pop_back(); }
                else           { used.erase(t.a.front()); s.pop_front(); t.a.erase(t.a.begin()); }
                break;
            }
            case 12: {  // 全消し -> build で作り直す
                if ((rng() % 20) != 0) break;
                std::vector<ll> init;
                used.clear();
                int n0 = (int)rnd(0, std::min(maxn, 12));
                for (int i = 0; i < n0; i++) {
                    ll x = fresh_value(used); if (x == TINF) break;
                    used.insert(x); init.push_back(x);
                }
                s.build(init); t.a = init;
                break;
            }
            case 13: {  // erase_range
                int l = (int)rnd(0, n - 1), r = (int)rnd(l + 1, n);
                for (int i = l; i < r; i++) used.erase(t.a[i]);
                s.erase_range(l, r);
                t.a.erase(t.a.begin() + l, t.a.begin() + r);
                break;
            }
            case 14: {  // insert_range
                int k = (int)rnd(1, 5);
                std::vector<ll> xs;
                for (int i = 0; i < k && (int)t.a.size() + (int)xs.size() < maxn; i++) {
                    ll x = fresh_value(used); if (x == TINF) break;
                    used.insert(x); xs.push_back(x);
                }
                if (xs.empty()) break;
                int pos = (int)rnd(0, n);
                s.insert_range(pos, xs);
                t.a.insert(t.a.begin() + pos, xs.begin(), xs.end());
                break;
            }
            case 15: {  // swap_values (存在する値 / しない値 の両方を試す)
                ll x = (rng() & 1) ? t.a[rnd(0, n - 1)] : rnd(vlo, vhi);
                ll y = (rng() & 1) ? t.a[rnd(0, n - 1)] : rnd(vlo, vhi);
                int i = t.index_of(x), j = t.index_of(y);
                eq(s.swap_values(x, y), i >= 0 && j >= 0, "swap_values の戻り値");
                if (i >= 0 && j >= 0) std::swap(t.a[i], t.a[j]);
                break;
            }
            default: {  // clear してから少しずつ積む
                if ((rng() % 30) != 0) break;
                s.clear(); t.a.clear(); used.clear();
                break;
            }
            }
            check_queries(s, t);
            if (do_validate && ((it & 7) == 0)) {
                std::string why;
                if (!validate(s, t.a, why)) { ng++; if (ng <= 5) std::cout << "  [NG] 内部検査: " << why << "  (" << ctx << ")\n"; }
            }
            if (ng > 0) return;   // 失敗したら以降は無意味なので止める
        }
        std::string why;
        if (do_validate && !validate(s, t.a, why)) { ng++; std::cout << "  [NG] 内部検査: " << why << "\n"; }
    }
};

static bool run_all_tests(ull seed) {
    std::cout << "=== ランダムテスト (seed = " << seed << ") ===\n";
    std::cout << "    愚直実装との突き合わせ + 木の内部不変条件の検査を行う\n";
    struct Case { const char* label; int iters; int maxn; ll lo, hi; bool validate; };
    const Case cases[] = {
        {"[1] 値が密 (-20..20)・小さい配列",      4000,  25,  -20,        20,        true },
        {"[2] 値が密 (-6..6)・極小配列",          4000,  8,   -6,         6,         true },
        {"[3] 値が疎 (±10^6)・中くらいの配列",    3000,  60,  -1000000,   1000000,   true },
        {"[4] 値が巨大 (±10^15)",                 2000,  120, -1000000000000000LL, 1000000000000000LL, true },
        {"[5] 大きい配列 (最大 800 要素)",        800,   800,  -1000000000, 1000000000, false},
    };
    long long total_checks = 0, total_ng = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (const Case& c : cases) {
        Tester ts(seed ^ (ull)(c.iters * 1000003LL + c.maxn));
        ts.run(c.iters, c.maxn, c.lo, c.hi, c.validate, c.label);
        total_checks += ts.checks; total_ng += ts.ng;
        std::cout << "  " << (ts.ng == 0 ? "OK  " : "NG  ") << c.label
                  << "   比較 " << ts.checks << " 項目"
                  << (ts.ng ? ("  失敗 " + std::to_string(ts.ng)) : std::string()) << "\n";
    }
    // 境界ケース: 空・1 要素
    {
        DynUniqueSeq s;
        bool ok = true;
        ok &= (s.size() == 0) && s.empty();
        ok &= (s.get_min() == TINF) && (s.get_max() == TNINF);
        ok &= (s.index_of(5) == -1) && !s.contains(5);
        ok &= (s.range_sum(0, 0) == 0) && (s.range_min(0, 0) == TINF) && (s.range_max(0, 0) == TNINF);
        ok &= (s.kth_smallest(0) == TINF) && (s.minimum_gap() == TINF);
        { auto ri = s.range_info(0, 0); ok &= (ri.sum == 0 && ri.mn == TINF && ri.mx == TNINF); }
        ok &= (s.range_min_with_index(0, 0) == std::pair<ll,int>{TINF, -1});
        s.push_back(7);
        ok &= (s.size() == 1) && (s.front() == 7) && (s.back() == 7) && (s.index_of(7) == 0);
        ok &= (s.range_second_min(0, 1) == TINF) && (s.minimum_gap() == TINF);
        ok &= (s.rank_of_value(7) == 0) && (s.rank_of_value(8) == 1) && (s.rank_of_index(0) == 0);
        s.reverse_range(0, 1); s.rotate_range(0, 0, 1); s.swap_at(0, 0);
        ok &= (s.to_vector() == std::vector<ll>{7});
        s.pop_back();
        ok &= s.empty() && (s.to_vector() == std::vector<ll>{});
        std::cout << "  " << (ok ? "OK  " : "NG  ") << "[6] 境界ケース (空 / 1 要素)\n";
        if (!ok) total_ng++;
        total_checks += 15;
    }
    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "--- 合計 " << total_checks << " 項目を比較、失敗 " << total_ng
              << " 件 (" << std::fixed << std::setprecision(0) << ms << " ms)\n";
    std::cout << (total_ng == 0 ? "*** 全テスト成功 ***\n" : "*** テスト失敗 ***\n") << "\n";
    return total_ng == 0;
}

// ---------------------------------------------------------------------------- 4) 動作確認デモ
static void show(const char* label, const DynUniqueSeq& s) {
    std::cout << std::left << std::setw(34) << label << " -> " << s.to_vector() << "\n";
}
template <class T>
static void line(const char* expr, const T& value, const char* comment) {
    std::cout << "  " << std::left << std::setw(36) << expr << " = "
              << std::setw(14) << value << "  // " << comment << "\n";
}

static void demo() {
    std::cout << "============================================================\n";
    std::cout << " DynUniqueSeq 動作確認デモ\n";
    std::cout << "  ・位置は 0-indexed、区間はすべて半開区間 [l, r)\n";
    std::cout << "  ・要素の値はすべて相異なる long long であること\n";
    std::cout << "  ・見つからないとき: 位置は -1 / 値は INF(min系) または -INF(max系)\n";
    std::cout << "============================================================\n\n";

    std::vector<ll> init = {50, 10, 80, 30, 70, 20};
    DynUniqueSeq s(init);
    std::cout << "[0] 構築  DynUniqueSeq s({50,10,80,30,70,20});\n";
    show("    s.to_vector()", s);
    std::cout << "           位置:      0   1   2   3   4   5\n\n";

    std::cout << "[1] 配列としての基本\n";
    line("s.size()", s.size(), "要素数");
    line("s.get(2)  /  s[2]", s.get(2), "位置 2 の値");
    line("s.front()", s.front(), "先頭");
    line("s.back()", s.back(), "末尾");
    line("s.total_sum()", s.total_sum(), "全体の総和 (O(1))");
    std::cout << "\n";

    std::cout << "[2] 区間集約  (区間 [1,5) = 10,80,30,70 が対象)\n";
    line("s.range_sum(1,5)", s.range_sum(1, 5), "区間の総和");
    line("s.range_min(1,5)", s.range_min(1, 5), "区間の最小値 (値だけ・速い)");
    line("s.range_max(1,5)", s.range_max(1, 5), "区間の最大値");
    {
        auto [v, i] = s.range_min_with_index(1, 5);
        line("s.range_min_with_index(1,5)", "{" + std::to_string(v) + "," + std::to_string(i) + "}", "{最小値, その位置}");
        auto [w, j] = s.range_max_with_index(1, 5);
        line("s.range_max_with_index(1,5)", "{" + std::to_string(w) + "," + std::to_string(j) + "}", "{最大値, その位置}");
    }
    line("s.range_min_index(1,5)", s.range_min_index(1, 5), "最小値の位置だけ");
    line("s.range_second_min(1,5)", s.range_second_min(1, 5), "区間で 2 番目に小さい値");
    line("s.range_second_max(1,5)", s.range_second_max(1, 5), "区間で 2 番目に大きい値");
    line("s.range_max_min_diff(1,5)", s.range_max_min_diff(1, 5), "区間の max - min");
    {
        auto ri = s.range_info(1, 5);
        std::cout << "  " << std::left << std::setw(36) << "s.range_info(1,5)" << " = "
                  << "{sum=" << ri.sum << ", mn=" << ri.mn << ", mx=" << ri.mx << "}"
                  << "   // 3 つまとめて 1 回の降下で取得 (別々に呼ぶより速い)\n";
    }
    line("s.get_min()", s.get_min(), "全体の最小値 (O(1))");
    line("s.max_index()", s.max_index(), "全体の最大値がある位置");
    std::cout << "\n";

    std::cout << "[3] 値 -> 位置 の逆引き  (値が全部違うので一意に決まる)\n";
    line("s.contains(80)", s.contains(80), "値 80 はあるか (1=ある)");
    line("s.contains(99)", s.contains(99), "値 99 はあるか (0=無い)");
    line("s.index_of(80)", s.index_of(80), "値 80 は今どこ?");
    line("s.index_of(99)", s.index_of(99), "無い値は -1");
    line("s.contains_in_range(0,3,80)", s.contains_in_range(0, 3, 80), "[0,3) に 80 はあるか");
    line("s.range_index_of(3,6,80)", s.range_index_of(3, 6, 80), "[3,6) に 80 は無いので -1");
    std::cout << "\n";

    std::cout << "[4] 値の順序 (set / 順序統計木 としての使い方)\n";
    std::cout << "      値を昇順に並べると  10 < 20 < 30 < 50 < 70 < 80\n";
    line("s.rank_of_value(30)", s.rank_of_value(30), "30 未満の要素数 (= 30 の順位)");
    line("s.rank_of_value(35)", s.rank_of_value(35), "存在しない値でも使える");
    line("s.kth_smallest(0)", s.kth_smallest(0), "0-indexed で k 番目に小さい値");
    line("s.kth_smallest(2)", s.kth_smallest(2), "3 番目に小さい値");
    line("s.kth_smallest_index(2)", s.kth_smallest_index(2), "その値がいる位置");
    line("s.kth_largest(0)", s.kth_largest(0), "k 番目に大きい値");
    line("s.rank_of_index(0)", s.rank_of_index(0), "位置 0 の値は何番目に小さいか");
    line("s.prev_value(50)", s.prev_value(50), "50 未満で最大の値");
    line("s.next_value(50)", s.next_value(50), "50 より大きく最小の値");
    line("s.lower_bound_value(35)", s.lower_bound_value(35), "35 以上で最小の値");
    line("s.upper_bound_value(30)", s.upper_bound_value(30), "30 より大きく最小の値");
    line("s.prev_value_index(50)", s.prev_value_index(50), "その値がいる位置");
    line("s.prev_index_by_value(0)", s.prev_index_by_value(0), "位置0の値(50)の1つ下の値がいる位置");
    line("s.count_less(50)", s.count_less(50), "50 未満の個数");
    line("s.count_ge(50)", s.count_ge(50), "50 以上の個数");
    line("s.count_value_range(20,50)", s.count_value_range(20, 50), "値が [20,50] の個数");
    line("s.sum_less(50)", s.sum_less(50), "値が 50 未満の要素の総和");
    line("s.sum_value_range(20,50)", s.sum_value_range(20, 50), "値が [20,50] の総和");
    line("s.sum_of_k_smallest(3)", s.sum_of_k_smallest(3), "小さい方から 3 個の総和");
    line("s.minimum_gap()", s.minimum_gap(), "値として隣り合う 2 要素の差の最小");
    {
        auto [p, q] = s.minimum_gap_pair();
        auto [pi, qi] = s.minimum_gap_indices();
        line("s.minimum_gap_pair()", "{" + std::to_string(p) + "," + std::to_string(q) + "}", "その 2 値");
        line("s.minimum_gap_indices()", "{" + std::to_string(pi) + "," + std::to_string(qi) + "}", "その 2 値の位置");
    }
    std::cout << "\n";

    std::cout << "[5] 更新 (挿入・削除・値の変更)\n";
    show("    最初", s);
    s.insert(2, 60);         show("    s.insert(2, 60)", s);
    s.push_back(5);          show("    s.push_back(5)", s);
    s.erase(0);              show("    s.erase(0)", s);
    s.erase_value(80);       show("    s.erase_value(80)", s);
    s.set(1, 15);            show("    s.set(1, 15)", s);
    s.set_value(70, 75);     show("    s.set_value(70, 75)", s);
    std::cout << "      ※ insert/set は「その値がまだ存在しない」ことが前提\n";
    std::cout << "      ※ erase_value / set_value は値が無ければ false を返す\n\n";

    std::cout << "[6] 並べ替え (位置の入れ替え系)\n";
    show("    現在", s);
    s.swap_at(0, 4);                 show("    s.swap_at(0, 4)", s);
    s.reverse_range(1, 5);           show("    s.reverse_range(1, 5)", s);
    s.rotate_range(0, 2, 5);         show("    s.rotate_range(0, 2, 5)", s);
    s.move_element(0, 3);            show("    s.move_element(0, 3)", s);
    s.move_range(0, 2, 3);           show("    s.move_range(0, 2, 3)", s);
    std::cout << "      ※ move_element(from,to) / move_range(l,r,pos) の to,pos は\n";
    std::cout << "         「抜き取った後の配列」での位置\n\n";

    std::cout << "[7] まとめて操作 / 一括取得\n";
    {
        DynUniqueSeq w({10, 20, 30, 40, 50, 60});
        show("    w = ", w);
        line("w.range_to_vector(1,4)", "[20,30,40]", "区間を配列で取り出す O(log N + k)");
        w.erase_range(1, 4);           show("    w.erase_range(1, 4)", w);
        w.insert_range(1, {15, 25});   show("    w.insert_range(1, {15,25})", w);
        w.swap_values(10, 60);         show("    w.swap_values(10, 60)", w);
        std::cout << "      ※ swap_values は位置を求めずに 2 値の位置を交換する\n";
        std::cout << "      w.sorted_values()                  = " << w.sorted_values() << "   // 値を昇順に O(N)\n";
        std::cout << "      w.value_order()                    = " << w.value_order() << "   // 値の昇順に並べたときの元の位置\n";
    }
    std::cout << "\n";

    std::cout << "[8] 区間の k 番目 (k が小さいとき専用 / O(k log N))\n";
    {
        DynUniqueSeq w({8, 3, 10, 5, 1, 7, 2});
        show("    w = ", w);
        std::cout << "      w.range_k_smallest(0,7,3)          = " << w.range_k_smallest(0, 7, 3)
                  << "   // 小さい方から 3 個 {値,位置}\n";
        std::cout << "      w.range_k_largest(0,7,2)           = " << w.range_k_largest(0, 7, 2)
                  << "   // 大きい方から 2 個\n";
        line("w.range_kth_smallest(1,6,2)", w.range_kth_smallest(1, 6, 2), "[1,6) で 3 番目に小さい値");
        std::cout << "      ※ 完全な区間 k 番目 (任意の k で O(log N)) は入れていない。理由は冒頭のコメント参照\n";
    }
    std::cout << "\n";

    std::cout << "[9] mex と 区間二分探索\n";
    {
        DynUniqueSeq w({3, 0, 1, 5, 2, 9});
        show("    w = ", w);
        line("w.mex()", w.mex(), "0 以上で配列に無い最小の整数");
        line("w.mex(5)", w.mex(5), "5 以上で配列に無い最小の整数");
        using RI = DynUniqueSeq::RangeInfo;
        ll v = 1;
        int r1 = w.max_right(1, [&](RI x) { return x.mn >= v; });
        std::cout << "      w.max_right(1, [](w){return w.mn>=1;}) = " << r1
                  << "        // 位置 1 から右へ「最小値が 1 以上」を保てる最大の r\n";
        ll d = 4;
        int r2 = w.max_right(0, [&](RI x) { return x.len == 0 || x.mx - x.mn <= d; });
        std::cout << "      w.max_right(0, max-min<=4)          = " << r2
                  << "        // 位置 0 から「max-min<=4」を保てる最大の r\n";
        int l1 = w.min_left(6, [&](RI x) { return x.len == 0 || x.mx - x.mn <= d; });
        std::cout << "      w.min_left(6, max-min<=4)           = " << l1
                  << "        // 右端 6 から左へ伸ばせる最小の l\n";
        std::cout << "      ※ f は「区間を伸ばすと true->false に一度だけ変わる」述語であること\n";
    }
    std::cout << "\n";

    std::cout << "[10] reverse しても 値->位置 は正しいまま (このテンプレートの要)\n";
    {
        DynUniqueSeq u({10, 20, 30, 40, 50});
        show("    u = ", u);
        std::cout << "      u.index_of(20) = " << u.index_of(20) << "\n";
        u.reverse_range(0, 5);
        show("    u.reverse_range(0,5)", u);
        std::cout << "      u.index_of(20) = " << u.index_of(20) << "   <- 追従している\n";
        std::cout << "      u.range_min_with_index(0,3) = {" << u.range_min_with_index(0, 3).first
                  << "," << u.range_min_with_index(0, 3).second << "}\n";
    }
    std::cout << "\n";
}

// ---------------------------------------------------------------------------- 5) ベンチマーク
static void bench() {
    std::cout << "=== ベンチマーク (g++ -O2 想定) ===\n";
    const int N = 200000, Q = 200000;
    std::mt19937_64 rng(20260920);

    std::vector<ll> a(N);
    {   // 相異なる値を用意
        for (int i = 0; i < N; i++) a[i] = i;
        std::shuffle(a.begin(), a.end(), rng);
        for (int i = 0; i < N; i++) a[i] = a[i] * 1000003LL - 500000000LL;
    }
    auto tm = [&](const char* name, auto f) {
        auto t0 = std::chrono::steady_clock::now();
        auto r = f();
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "  " << std::left << std::setw(46) << name << std::right
                  << std::fixed << std::setprecision(1) << std::setw(8) << ms << " ms"
                  << "   (checksum " << r << ")\n";
    };

    DynUniqueSeq s;
    tm("build  N=200000", [&] { s.build(a); return s.total_sum(); });

    tm("range_min のみ  20万回", [&] {
        ull acc = 0;
        for (int q = 0; q < Q; q++) {
            int l = (int)(rng() % N), r = (int)(rng() % N);
            if (l > r) std::swap(l, r);
            acc += (ull)s.range_min(l, r + 1);
        }
        return acc;
    });
    tm("range_min_with_index (位置付き)  20万回", [&] {
        ull acc = 0;
        for (int q = 0; q < Q; q++) {
            int l = (int)(rng() % N), r = (int)(rng() % N);
            if (l > r) std::swap(l, r);
            acc += (ull)s.range_min_with_index(l, r + 1).second;
        }
        return acc;
    });
    tm("range_info (sum/min/max まとめて)  20万回", [&] {
        ull acc = 0;
        for (int q = 0; q < Q; q++) {
            int l = (int)(rng() % N), r = (int)(rng() % N);
            if (l > r) std::swap(l, r);
            auto w = s.range_info(l, r + 1);
            acc += (ull)w.sum ^ (ull)w.mn ^ (ull)w.mx;
        }
        return acc;
    });
    tm("range_min + range_max + range_sum 別々  20万回", [&] {
        ull acc = 0;
        for (int q = 0; q < Q; q++) {
            int l = (int)(rng() % N), r = (int)(rng() % N);
            if (l > r) std::swap(l, r);
            r++;
            acc += (ull)s.range_min(l, r) ^ (ull)s.range_max(l, r) ^ (ull)s.range_sum(l, r);
        }
        return acc;
    });
    tm("index_of (値->位置)  20万回", [&] {
        ull acc = 0;
        for (int q = 0; q < Q; q++) acc += (ull)s.index_of(a[rng() % N]);
        return acc;
    });
    tm("kth_smallest + rank_of_value  20万回", [&] {
        ull acc = 0;
        for (int q = 0; q < Q; q++) {
            acc += (ull)s.kth_smallest((int)(rng() % N));
            acc += (ull)s.rank_of_value((ll)(rng() % 1000000000ULL) - 500000000LL);
        }
        return acc;
    });
    tm("prev_value / next_value  20万回", [&] {
        ull acc = 0;
        for (int q = 0; q < Q; q++) {
            ll x = (ll)(rng() % 1000000000ULL) - 500000000LL;
            acc += (ull)(s.prev_value(x) / 2) + (ull)(s.next_value(x) / 2);
        }
        return acc;
    });
    tm("insert + erase  各 10万回", [&] {
        ull acc = 0;
        for (int q = 0; q < 100000; q++) {
            ll x = -1000000000000LL - q;
            s.insert((int)(rng() % (s.size() + 1)), x);
            acc += (ull)s.erase_value(x);
        }
        return acc;
    });
    tm("range_kth_smallest (k=4)  20万回", [&] {
        ull acc = 0;
        for (int q = 0; q < Q; q++) {
            int l = (int)(rng() % N), r = (int)(rng() % N);
            if (l > r) std::swap(l, r);
            acc += (ull)s.range_kth_smallest(l, r + 1, 3);
        }
        return acc;
    });
    tm("mex  20万回", [&] {
        ull acc = 0;
        for (int q = 0; q < Q; q++) acc += (ull)s.mex((ll)(rng() % 1000000000ULL) - 500000000LL);
        return acc;
    });
    tm("max_right (min>=v)  20万回", [&] {
        ull acc = 0;
        using RI = DynUniqueSeq::RangeInfo;
        for (int q = 0; q < Q; q++) {
            ll v = (ll)(rng() % 1000000000ULL) - 500000000LL;
            acc += (ull)s.max_right((int)(rng() % N), [&](RI w) { return w.mn >= v; });
        }
        return acc;
    });
    tm("reverse_range  20万回", [&] {
        for (int q = 0; q < Q; q++) {
            int l = (int)(rng() % N), r = (int)(rng() % N);
            if (l > r) std::swap(l, r);
            s.reverse_range(l, r + 1);
        }
        return s.total_sum();
    });
    tm("reverse 後の index_of  20万回", [&] {
        ull acc = 0;
        for (int q = 0; q < Q; q++) acc += (ull)s.index_of(a[rng() % N]);
        return acc;
    });
    std::cout << "\n";
}

void solv(){
    std::string mode = "";
    ull seed = 20260920ULL;

    if (mode == "demo") { demo(); return; }
    if (mode == "test") { run_all_tests(seed);return; }
    if (mode == "bench") { bench(); return; }
    demo();
    bool ok = run_all_tests(seed);
    bench();
    return;
}

int main(){
    BUFF;
    int T=1;
    //cin>>T;
    while(T--)solv();
    return 0;
}