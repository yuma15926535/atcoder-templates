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

// =============================================================================
//  辺の削除もできる Union-Find（オンライン版 / Online Dynamic Connectivity）  --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  辺の追加・削除をしながら、その場で「u と v は連結か」「x の連結成分の大きさ」などに答える。
//
//  ★ どれを使うか
//    - 前の答えで次の入力が決まる（強制オンライン）など、その場で答える必要がある → このファイル
//    - 全部の入力を読んでから答えてよい → dynamic_connectivity_offline.cpp（オフライン版。削除が多いと数倍速い）
//    - 削除しかない → 逆順に見て「追加だけ」にし、普通の UnionFind
//
//  仕組み（Holm-de Lichtenberg-Thorup 法）
//    - 辺を「全域森に使う木辺」と「それ以外の非木辺」に分け、全域森を Euler Tour Tree で持つ
//    - 辺にレベル（0 ~ log2 n）を付け、レベル i 以上の木辺からなる森 F_i を各レベルで持つ
//    - 木辺を消したら、小さい側の連結成分から「反対側につながる非木辺」を探して代わりにする。
//      調べて使えなかった辺はレベルを上げる（各辺のレベルは log2 n 回までしか上がらない）
//
//  頂点番号   : 0-indexed（0 ~ n-1）
//  依存       : 標準ライブラリのみ（C++14 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  高速化のポイント
//    - Euler Tour Tree は親ポインタ付き splay 木。番兵ノード 0 で分岐を減らし、ノードは配列で管理
//    - link は splay 2 回 + 列の回転 1 回だけで済むよう、u の直後に v の木を差し込む
//    - 木辺を消したとき、代わりの辺がすぐ見つかれば木辺のレベル上げ（重い処理）を省略する
//    - 多重辺の削除では非木辺を優先して消す（木辺を消すより圧倒的に軽い）
//    - 辺の検索は自作のハッシュ表（開番地法）。std::unordered_map より速い
//    - 上位レベルは必要になってから作る。頂点の値の合計は set_value を使うまで計算しない
//
//  機能一覧                                                                      計算量
//    OnlineDynamicConnectivity dc(n)   頂点数 n・辺なしで作る                      O(n)
//    dc.add_edge(u, v)     -> bool     辺を追加。別々の連結成分がつながったら true   ならし O(log n)
//    dc.remove_edge(u, v)  -> bool     辺を 1 本削除。連結成分が分かれたら true       ならし O(log^2 n)
//    dc.same(u, v)         -> bool     u と v が連結なら true                       ならし O(log n)
//    dc.size(x)            -> int      x の連結成分の頂点数                         ならし O(log n)
//    dc.group_count()      -> int      連結成分の個数                               O(1)
//    ---- 追加機能 ----
//    dc.members(x)         -> vector<int>          x の連結成分の頂点（順不同）   O(成分の大きさ)
//    dc.groups()           -> vector<vector<int>>  全ての連結成分                 O(n log n)
//    dc.has_edge(u, v)     -> bool     辺 (u, v) が今あるか                         O(1)
//    dc.edge_count()       -> int      今ある辺の本数                               O(1)
//    dc.edges()            -> vector<pair<int,int>>  今ある辺の一覧               O(m log m)
//    dc.set_value(v, x) / dc.add_value(v, x) / dc.value(v)   頂点の値（初期値 0）  ならし O(log n)
//    dc.component_sum(v)   -> long long  v の連結成分の値の合計                     ならし O(log n)
//    dc.init(n) / dc.n()               作り直す / 頂点数
//    dc.level_count()                  （デバッグ用）内部のレベル数
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    dc_test::random_test()  … 愚直解（毎回グラフから連結成分を求め直す）と全機能を突き合わせる
//    dc_test::benchmark()    … 速度計測
//    dc_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================

// ======================= ここからテンプレート本体（コピー範囲） =======================

// -----------------------------------------------------------------------------
//  OnlineDynamicConnectivity: 辺の追加・削除をしながら、その場で連結性を答えられる Union-Find
//    Holm-de Lichtenberg-Thorup (HDT) 法。各レベルの全域森を Euler Tour Tree（splay 木）で持つ。
//    add_edge: ならし O(log n)   remove_edge: ならし O(log^2 n)   same / size: ならし O(log n)
// -----------------------------------------------------------------------------
struct OnlineDynamicConnectivity {
  public:
    OnlineDynamicConnectivity() : OnlineDynamicConnectivity(0) {}
    explicit OnlineDynamicConnectivity(int n) { init(n); }

    // 頂点数 n・辺なしで作り直す（頂点の値も 0 に戻る）  O(n)
    void init(int n) {
        n_ = n;
        cnt_ = n;
        m_ = 0;
        lv_.clear();
        lv_.reserve(64);  // レベル数は log2(n) + 1 以下なので再確保は起きない
        head_.clear();
        enode_.clear();
        eu_.clear(), ev_.clear(), elev_.clear(), echain_.clear(), efree_.clear(), snext_.clear(), sprev_.clear();
        etree_.clear();
        emap_.clear();
        ensure_level(0);
    }

    // 頂点数
    int n() const { return n_; }
    // 連結成分の個数  O(1)
    int group_count() const { return cnt_; }
    // 今ある辺の本数（多重辺・自己ループも 1 本ずつ数える）  O(1)
    int edge_count() const { return m_; }

    // 辺 (u, v) を追加する。別々の連結成分がつながったら true  ならし O(log n)
    bool add_edge(int u, int v) {
        assert(0 <= u && u < n_ && 0 <= v && v < n_);
        const int e = new_edge(u, v);
        if (u == v) {  // 自己ループは連結性に関係しないので辺として数えるだけ
            elev_[e] = -1;
            etree_[e] = 0;
            return false;
        }
        if (!lv_[0].connected(u + 1, v + 1)) {
            tree_link(e, 0);
            --cnt_;
            return true;
        }
        nt_insert(e, 0);
        return false;
    }

    // 辺 (u, v) を 1 本削除する（存在する辺であること）。連結成分が分かれたら true  ならし O(log^2 n)
    bool remove_edge(int u, int v) {
        assert(0 <= u && u < n_ && 0 <= v && v < n_);
        const int e = take_edge(u, v);
        if (elev_[e] < 0 || !etree_[e]) {  // 自己ループ / 木に使っていない辺: 消すだけ
            if (elev_[e] >= 0) nt_erase(e);
            free_edge(e);
            return false;
        }
        const int l = elev_[e];
        const int eu = eu_[e], ev = ev_[e];
        int su[64], sv[64];  // レベルごとの、切った後の eu 側・ev 側の頂点数
        for (int i = 0; i <= l; ++i) {
            const std::pair<int, int> s = lv_[i].cut(enode_[i][e]);
            su[i] = s.first, sv[i] = s.second;
        }
        free_edge(e);
        for (int i = l; i >= 0; --i)  // 高いレベルから代わりの辺を探す（小さい側から調べる）
            if (su[i] <= sv[i] ? replace(eu, i) : replace(ev, i)) return false;
        ++cnt_;
        return true;
    }

    // 辺 (u, v) が 1 本以上あるか  O(1)
    bool has_edge(int u, int v) const { return emap_.find(key(u, v)) >= 0; }

    // u と v が連結なら true  ならし O(log n)
    bool same(int u, int v) { return u == v || lv_[0].connected(u + 1, v + 1); }

    // x の連結成分の頂点数  ならし O(log n)
    int size(int x) { return lv_[0].tree_size(x + 1); }

    // x の連結成分の頂点一覧（順不同）  O(成分の大きさ)
    std::vector<int> members(int x) {
        std::vector<int> res;
        lv_[0].collect(x + 1, res);
        return res;
    }

    // 全ての連結成分（最小要素の昇順・各成分内も昇順）  O(n log n)
    std::vector<std::vector<int>> groups() {
        std::vector<char> seen(n_, 0);
        std::vector<std::vector<int>> res;
        for (int v = 0; v < n_; ++v) {
            if (seen[v]) continue;
            std::vector<int> g = members(v);
            std::sort(g.begin(), g.end());
            for (int w : g) seen[w] = 1;
            res.push_back(std::move(g));
        }
        return res;
    }

    // 今ある辺の一覧 (u, v)（u <= v、昇順。多重辺はその本数だけ並ぶ）  O(m log m)
    std::vector<std::pair<int, int>> edges() const {
        std::vector<std::pair<int, int>> res;
        for (size_t i = 0; i < emap_.keys.size(); ++i)
            if (emap_.keys[i] != EdgeMap::empty())
                for (int e = emap_.vals[i]; e >= 0; e = echain_[e]) res.emplace_back(eu_[e], ev_[e]);
        for (auto& p : res)
            if (p.first > p.second) std::swap(p.first, p.second);
        std::sort(res.begin(), res.end());
        return res;
    }

    // 頂点 v の値を x にする / x を足す / 取得する（初期値は 0）  ならし O(log n)
    //   （値を一度も設定しなければ、合計を管理する分の時間はかからない）
    void set_value(int v, long long x) {
        if (!lv_[0].with_sum) lv_[0].enable_sum();
        lv_[0].set_val(v + 1, x);
    }
    void add_value(int v, long long x) { set_value(v, value(v) + x); }
    long long value(int v) const { return lv_[0].with_sum ? lv_[0].val[v + 1] : 0; }
    // v の連結成分に含まれる頂点の値の合計  ならし O(log n)
    long long component_sum(int v) { return lv_[0].with_sum ? lv_[0].tree_sum(v + 1) : 0; }

    // （デバッグ用）内部で使っているレベルの数。常に floor(log2(n)) + 1 以下
    int level_count() const { return (int)lv_.size(); }

  private:
    // ----- Euler Tour Tree（splay 木）: 1 つのレベルの全域森を持つ -----
    //   ノード 0 は番兵（空）、1 ~ n は頂点ノード、n+1 以降は辺ノード（有向辺 2 本を連番の組で確保）
    struct Ett {
        struct Node {
            int ch[2], p;
            int sz;              // 部分木に含まれる頂点ノードの数
            int eid;             // 辺ノードなら辺番号
            unsigned char own;   // bit0: このレベルの非木辺を持つ頂点, bit1: レベルがちょうどここの木辺
            unsigned char agg;   // 部分木の own の OR
        };
        int n = 0;
        bool with_sum = false;  // 頂点の値の合計を持つか（レベル 0 で set_value が呼ばれてから有効）
        std::vector<Node> t;
        std::vector<long long> val, sum;
        std::vector<int> free_pairs;

        static Node make(int sz, int eid) {
            Node x;
            x.ch[0] = x.ch[1] = x.p = 0;
            x.sz = sz;
            x.eid = eid;
            x.own = x.agg = 0;
            return x;
        }
        void init(int num) {
            n = num;
            with_sum = false;
            t.assign(n + 1, make(1, -1));
            t[0].sz = 0;
            val.clear(), sum.clear();
            free_pairs.clear();
        }
        // 頂点の値の合計を使い始める（それまでは値が全部 0 なので合計も全部 0 で正しい）
        void enable_sum() {
            with_sum = true;
            val.assign(t.size(), 0);
            sum.assign(t.size(), 0);
        }
        void update(int x) {
            Node& a = t[x];
            const Node& L = t[a.ch[0]];
            const Node& R = t[a.ch[1]];
            a.sz = L.sz + R.sz + (x <= n ? 1 : 0);
            a.agg = (unsigned char)(a.own | L.agg | R.agg);
            if (with_sum) sum[x] = sum[a.ch[0]] + sum[a.ch[1]] + val[x];
        }
        void rotate(int x) {
            const int y = t[x].p, z = t[y].p;
            const int dx = (t[y].ch[1] == x);
            const int b = t[x].ch[dx ^ 1];
            t[y].ch[dx] = b;
            if (b) t[b].p = y;
            t[x].ch[dx ^ 1] = y;
            t[y].p = x;
            t[x].p = z;
            if (z) t[z].ch[t[z].ch[1] == y] = x;
            update(y);
        }
        void splay(int x) {
            while (t[x].p) {
                const int y = t[x].p, z = t[y].p;
                if (z) rotate(((t[z].ch[0] == y) == (t[y].ch[0] == x)) ? y : x);
                rotate(x);
            }
            update(x);
        }
        // 根 a の木の後ろに根 b の木をつなげ、新しい根を返す（どちらかが 0 でもよい）
        int join(int a, int b) {
            if (!a) return b;
            if (!b) return a;
            int x = a;
            while (t[x].ch[1]) x = t[x].ch[1];
            splay(x);
            t[x].ch[1] = b;
            t[b].p = x;
            update(x);
            return x;
        }
        // x が先頭に来るよう Euler Tour を回転し、根を返す
        int reroot(int x) {
            splay(x);
            const int L = t[x].ch[0];
            if (!L) return x;
            t[x].ch[0] = 0;
            t[L].p = 0;
            update(x);
            return join(x, L);
        }
        bool connected(int a, int b) {
            if (a == b) return true;
            splay(a);
            splay(b);
            return t[a].p != 0;
        }
        int tree_size(int a) {
            splay(a);
            return t[a].sz;
        }
        long long tree_sum(int a) {
            splay(a);
            return sum[a];
        }
        void set_val(int a, long long x) {
            splay(a);
            val[a] = x;
            update(a);
        }
        int new_pair(int e) {
            int a;
            if (!free_pairs.empty()) {
                a = free_pairs.back();
                free_pairs.pop_back();
            } else {
                a = (int)t.size();
                t.resize(a + 2);
                if (with_sum) val.resize(a + 2), sum.resize(a + 2);
            }
            t[a] = t[a + 1] = make(0, e);
            if (with_sum) val[a] = val[a + 1] = sum[a] = sum[a + 1] = 0;
            return a;
        }
        // 頂点 u, v（ノード番号、別々の木）をつなぐ。辺ノードの組の先頭番号 a を返す
        //   u の木の列 [L, u, R] の u の直後に「a, (v から始まる v の木の列), a+1」を差し込む
        int link(int u, int v, int e) {
            const int a = new_pair(e), b = a + 1;
            const int rv = reroot(v);
            splay(u);
            const int R = t[u].ch[1];
            t[u].ch[1] = 0;
            update(u);
            t[a].ch[0] = u, t[u].p = a;
            t[a].ch[1] = rv, t[rv].p = a;
            update(a);
            t[b].ch[0] = a, t[a].p = b;
            t[b].ch[1] = R;
            if (R) t[R].p = b;
            update(b);
            return a;
        }
        // 辺ノードの組 a（a: u->v, a+1: v->u）を切る。切った後の (u 側の頂点数, v 側の頂点数) を返す
        std::pair<int, int> cut(int a) {
            const int b = a + 1;
            splay(a);
            int x = b;
            while (t[x].p != a) x = t[x].p;
            int first = a, second = b;
            if (t[a].ch[0] == x) {  // b が a より前にある
                first = b, second = a;
                splay(first);
            }
            const int L = t[first].ch[0], R = t[first].ch[1];
            t[first].ch[0] = t[first].ch[1] = 0;
            t[L].p = t[R].p = 0;
            splay(second);
            const int M = t[second].ch[0], R2 = t[second].ch[1];
            t[second].ch[0] = t[second].ch[1] = 0;
            t[M].p = t[R2].p = 0;
            const int rest = join(L, R2);
            free_pairs.push_back(a);
            // 列が [.., u->v, M, v->u, ..] なら M が v 側、[.., v->u, M, u->v, ..] なら M が u 側
            const int sm = t[M].sz, sr = t[rest].sz;
            return first == a ? std::make_pair(sr, sm) : std::make_pair(sm, sr);
        }
        void set_own(int x, int bit, bool on) {
            splay(x);
            if (on) t[x].own |= (unsigned char)bit;
            else t[x].own &= (unsigned char)~bit;
            update(x);
        }
        // x と同じ木で、own に bit が立っているノードを 1 つ返す（なければ 0）
        int find_flag(int x, int bit) {
            splay(x);
            if (!(t[x].agg & bit)) return 0;
            while (!(t[x].own & bit)) x = (t[t[x].ch[0]].agg & bit) ? t[x].ch[0] : t[x].ch[1];
            splay(x);
            return x;
        }
        // a と同じ木の頂点（0-indexed）を全部 out に入れる
        void collect(int a, std::vector<int>& out) {
            splay(a);
            out.reserve(out.size() + t[a].sz);
            std::vector<int> st(1, a);
            while (!st.empty()) {
                const int x = st.back();
                st.pop_back();
                if (x <= n) out.push_back(x - 1);
                if (t[x].ch[0]) st.push_back(t[x].ch[0]);
                if (t[x].ch[1]) st.push_back(t[x].ch[1]);
            }
        }
    };

    int n_ = 0, cnt_ = 0, m_ = 0;
    std::vector<Ett> lv_;                    // lv_[i] = レベル i 以上の木辺からなる全域森
    std::vector<std::vector<int>> head_;     // head_[i][v] = 頂点 v のレベル i の非木辺リストの先頭スロット
    std::vector<std::vector<int>> enode_;    // enode_[i][e] = 木辺 e のレベル i の ETT での辺ノード
    std::vector<int> eu_, ev_, elev_;        // 辺の端点・レベル（自己ループは -1）
    std::vector<char> etree_;                // 木辺なら 1
    std::vector<int> echain_, efree_;        // 同じ端点の辺のリスト / 再利用する辺番号
    std::vector<int> snext_, sprev_;         // 非木辺リスト（スロット = 2 * 辺番号 + 端点側）
    // (u, v) -> その端点を持つ辺の番号（多重辺は echain_ でつながる）のハッシュ表。
    //   開番地法（線形探索）、削除は後ろの要素を詰める方式。std::unordered_map より速い
    struct EdgeMap {
        static unsigned long long empty() { return ~0ULL; }  // 空き番地の印
        std::vector<unsigned long long> keys;
        std::vector<int> vals;
        int used = 0;
        void clear() {
            keys.assign(16, empty());
            vals.assign(16, -1);
            used = 0;
        }
        static size_t hash(unsigned long long x) {
            x ^= x >> 31;
            x *= 0x9e3779b97f4a7c15ULL;
            x ^= x >> 29;
            return (size_t)x;
        }
        size_t slot(unsigned long long k) const {
            const size_t m = keys.size() - 1;
            size_t i = hash(k) & m;
            while (keys[i] != empty() && keys[i] != k) i = (i + 1) & m;
            return i;
        }
        // k の値の位置（なければ -1）
        int find(unsigned long long k) const {
            const size_t i = slot(k);
            return keys[i] == k ? (int)i : -1;
        }
        // k の位置を返す。なければ値 -1 で作る
        size_t insert(unsigned long long k) {
            if (2 * (size_t)(used + 1) > keys.size()) {
                std::vector<unsigned long long> ok;
                std::vector<int> ov;
                ok.swap(keys), ov.swap(vals);
                keys.assign(ok.size() * 2, empty());
                vals.assign(ok.size() * 2, -1);
                for (size_t i = 0; i < ok.size(); ++i)
                    if (ok[i] != empty()) {
                        const size_t j = slot(ok[i]);
                        keys[j] = ok[i], vals[j] = ov[i];
                    }
            }
            const size_t i = slot(k);
            if (keys[i] == empty()) keys[i] = k, vals[i] = -1, ++used;
            return i;
        }
        void erase_at(size_t i) {
            const size_t m = keys.size() - 1;
            for (size_t j = i;;) {
                j = (j + 1) & m;
                if (keys[j] == empty()) break;
                const size_t h = hash(keys[j]) & m;
                if (i <= j ? (i < h && h <= j) : (i < h || h <= j)) continue;  // j はそのままでよい
                keys[i] = keys[j], vals[i] = vals[j], i = j;
            }
            keys[i] = empty(), vals[i] = -1, --used;
        }
    };
    EdgeMap emap_;

    static unsigned long long key(int u, int v) {
        if (u > v) std::swap(u, v);
        return (unsigned long long)u << 32 | (unsigned)v;
    }
    void ensure_level(int i) {
        while ((int)lv_.size() <= i) {
            lv_.emplace_back();
            lv_.back().init(n_);
            head_.emplace_back(n_, -1);
            enode_.emplace_back(eu_.size(), 0);
        }
    }
    int new_edge(int u, int v) {
        int e;
        if (!efree_.empty()) {
            e = efree_.back();
            efree_.pop_back();
        } else {
            e = (int)eu_.size();
            eu_.push_back(0), ev_.push_back(0), elev_.push_back(0), etree_.push_back(0), echain_.push_back(-1);
            snext_.push_back(-1), snext_.push_back(-1), sprev_.push_back(-1), sprev_.push_back(-1);
            for (auto& en : enode_) en.push_back(0);
        }
        eu_[e] = u, ev_[e] = v;
        const size_t pos = emap_.insert(key(u, v));  // 同じ端点の辺があれば連鎖の先頭に入れる
        echain_[e] = emap_.vals[pos];
        emap_.vals[pos] = e;
        ++m_;
        return e;
    }
    // (u, v) の辺を 1 本取り出す（木に使っていない辺を優先）。存在しなければ assert
    int take_edge(int u, int v) {
        const int pos = emap_.find(key(u, v));
        assert(pos >= 0);  // 存在しない辺の削除
        int& head = emap_.vals[pos];
        int best = head, bp = -1;
        for (int prev = -1, cur = head; cur >= 0; prev = cur, cur = echain_[cur])
            if (!etree_[cur]) {
                best = cur, bp = prev;
                break;
            }
        if (bp < 0) head = echain_[best];
        else echain_[bp] = echain_[best];
        if (head < 0) emap_.erase_at(pos);
        --m_;
        return best;
    }
    void free_edge(int e) { efree_.push_back(e); }

    // 非木辺 e をレベル i のリストに入れる
    void nt_insert(int e, int i) {
        ensure_level(i);
        elev_[e] = i;
        etree_[e] = 0;
        std::vector<int>& h = head_[i];
        for (int side = 0; side < 2; ++side) {
            const int v = side ? ev_[e] : eu_[e], s = 2 * e + side, old = h[v];
            snext_[s] = old;
            sprev_[s] = -1;
            if (old >= 0) sprev_[old] = s;
            else lv_[i].set_own(v + 1, 1, true);
            h[v] = s;
        }
    }
    // 非木辺 e をリストから外す
    void nt_erase(int e) {
        const int i = elev_[e];
        std::vector<int>& h = head_[i];
        for (int side = 0; side < 2; ++side) {
            const int v = side ? ev_[e] : eu_[e], s = 2 * e + side, pv = sprev_[s], nx = snext_[s];
            if (pv >= 0) snext_[pv] = nx;
            else h[v] = nx;
            if (nx >= 0) sprev_[nx] = pv;
            if (pv < 0 && nx < 0) lv_[i].set_own(v + 1, 1, false);
        }
    }
    // 辺 e をレベル lvl の木辺にする（レベル 0 ~ lvl の森に link）
    void tree_link(int e, int lvl) {
        etree_[e] = 1;
        elev_[e] = lvl;
        for (int i = 0; i <= lvl; ++i) enode_[i][e] = lv_[i].link(eu_[e] + 1, ev_[e] + 1, e);
        lv_[lvl].set_own(enode_[lvl][e], 2, true);
    }
    // レベル i の木辺 e をレベル i + 1 に上げる
    void push_tree_edge(int e, int i) {
        ensure_level(i + 1);
        lv_[i].set_own(enode_[i][e], 2, false);
        elev_[e] = i + 1;
        enode_[i + 1][e] = lv_[i + 1].link(eu_[e] + 1, ev_[e] + 1, e);
        lv_[i + 1].set_own(enode_[i + 1][e], 2, true);
    }
    // 木辺を切った直後に、レベル i で「u を含む小さい側」から反対側へつながる代わりの辺を探す。
    // 見つかればつないで true
    bool replace(int u, int i) {
        bool pushed = false;
        for (int x; (x = lv_[i].find_flag(u + 1, 1)) != 0;) {  // u 側でレベル i の非木辺を持つ頂点
            const int w = x - 1;
            for (int s = head_[i][w]; s >= 0;) {
                const int nx = snext_[s], e = s >> 1;
                const int y = (eu_[e] == w) ? ev_[e] : eu_[e];
                if (!lv_[i].connected(w + 1, y + 1)) {  // 反対側につながる辺 → これを木辺にして終わり
                    nt_erase(e);
                    tree_link(e, i);
                    return true;
                }
                if (!pushed) {  // 非木辺のレベルを上げる前に、u 側のレベル i の木辺を全部 i + 1 に上げる
                    for (int z; (z = lv_[i].find_flag(u + 1, 2)) != 0;) push_tree_edge(lv_[i].t[z].eid, i);
                    pushed = true;
                }
                nt_erase(e);  // u 側の中で閉じている辺 → レベルを上げる
                nt_insert(e, i + 1);
                s = nx;
            }
        }
        return false;
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト: 愚直解（今ある辺から毎回 BFS で連結成分を求め直す）と全機能を突き合わせる
// =============================================================================
namespace dc_test {

// NDEBUG（VS の Release ビルド）でも無効にならないチェック
#define DC_CHECK(cond)                                                                                        \
    do {                                                                                                      \
        if (!(cond)) {                                                                                        \
            std::cerr << "[OnlineDynamicConnectivity] TEST FAILED  line " << __LINE__ << ": " << #cond << "\n"; \
            std::exit(1);                                                                                     \
        }                                                                                                     \
    } while (0)

// テスト量の倍率。デバッグビルド（最適化なし・範囲チェックあり）は非常に遅いので自動で 1/8 にする。
// もっと念入りに検査したいときはこの値を大きくする
#if defined(_DEBUG) || defined(_GLIBCXX_DEBUG) || (defined(__GNUC__) && !defined(__OPTIMIZE__))
static const int TEST_SCALE = 1;
#else
static const int TEST_SCALE = 8;
#endif

static int rnd(std::mt19937_64& rng, int lo, int hi) {  // [lo, hi]
    return lo + (int)(rng() % (unsigned long long)(hi - lo + 1));
}
static int floor_log2(int n) {
    int k = 0;
    while ((2 << k) <= n) ++k;
    return k;
}

// 愚直解: 今ある辺（多重辺・自己ループ込み）と頂点の値をそのまま持ち、連結成分は BFS で求め直す
struct Naive {
    int n = 0;
    std::vector<std::pair<int, int>> alive;  // 今ある辺（向きは登録したまま）
    std::vector<long long> val;
    std::vector<int> lab;  // 連結成分番号（最小頂点の順に 0, 1, 2, ...）
    bool dirty = true;

    void init(int n_) {
        n = n_;
        alive.clear();
        val.assign(n, 0);
        dirty = true;
    }
    const std::vector<int>& labels() {
        if (dirty) {
            std::vector<std::vector<int>> g(n);
            for (const auto& e : alive) g[e.first].push_back(e.second), g[e.second].push_back(e.first);
            lab.assign(n, -1);
            std::vector<int> q;
            int k = 0;
            for (int s = 0; s < n; ++s) {
                if (lab[s] >= 0) continue;
                lab[s] = k;
                q.assign(1, s);
                for (size_t h = 0; h < q.size(); ++h)
                    for (int w : g[q[h]])
                        if (lab[w] < 0) lab[w] = k, q.push_back(w);
                ++k;
            }
            dirty = false;
        }
        return lab;
    }
    int group_count() {
        int k = 0;
        for (int x : labels()) k = std::max(k, x + 1);
        return k;
    }
    std::vector<int> group_of(int x) {
        std::vector<int> res;
        for (int v = 0; v < n; ++v)
            if (labels()[v] == labels()[x]) res.push_back(v);
        return res;
    }
    long long sum_of(int x) {
        long long s = 0;
        for (int v = 0; v < n; ++v)
            if (labels()[v] == labels()[x]) s += val[v];
        return s;
    }
    std::vector<std::vector<int>> groups() {
        std::vector<std::vector<int>> g(group_count());
        for (int v = 0; v < n; ++v) g[labels()[v]].push_back(v);
        return g;
    }
    int count_edge(int u, int v) const {
        int c = 0;
        for (const auto& e : alive) c += (e.first == u && e.second == v) || (e.first == v && e.second == u);
        return c;
    }
    std::vector<std::pair<int, int>> sorted_edges() const {
        std::vector<std::pair<int, int>> res = alive;
        for (auto& p : res)
            if (p.first > p.second) std::swap(p.first, p.second);
        std::sort(res.begin(), res.end());
        return res;
    }
};

// 全機能の全体チェック。full = true なら全頂点の members と全ペアの same も見る
static void check_all(OnlineDynamicConnectivity& dc, Naive& nv, std::mt19937_64& rng, bool full) {
    const int n = nv.n;
    const std::vector<int> lab = nv.labels();
    DC_CHECK(dc.n() == n);
    DC_CHECK(dc.group_count() == nv.group_count());
    DC_CHECK(dc.edge_count() == (int)nv.alive.size());
    DC_CHECK(dc.level_count() <= floor_log2(std::max(n, 1)) + 1);
    std::vector<int> sz(n, 0);
    std::vector<long long> sm(n, 0);
    for (int v = 0; v < n; ++v) ++sz[lab[v]], sm[lab[v]] += nv.val[v];
    for (int v = 0; v < n; ++v) {
        DC_CHECK(dc.size(v) == sz[lab[v]]);
        DC_CHECK(dc.value(v) == nv.val[v]);
        DC_CHECK(dc.component_sum(v) == sm[lab[v]]);
    }
    DC_CHECK(dc.groups() == nv.groups());
    DC_CHECK(dc.edges() == nv.sorted_edges());
    if (full) {
        for (int x = 0; x < n; ++x) {
            std::vector<int> m = dc.members(x);
            std::sort(m.begin(), m.end());
            DC_CHECK(m == nv.group_of(x));
            for (int y = 0; y < n; ++y) DC_CHECK(dc.same(x, y) == (lab[x] == lab[y]));
        }
    } else if (n > 0) {
        for (int k = 0; k < 2 * n; ++k) {
            const int x = rnd(rng, 0, n - 1), y = rnd(rng, 0, n - 1);
            DC_CHECK(dc.same(x, y) == (lab[x] == lab[y]));
        }
    }
}

// ランダム操作列を 1 本流す。add_pct = 辺の追加の割合（残りの半分ずつが削除と質問）
static long long run_case(std::mt19937_64& rng, int n, int q, int add_pct) {
    OnlineDynamicConnectivity dc(n);
    Naive nv;
    nv.init(n);
    const bool small = (n <= 40);
    const int check_every = std::max(8, n);
    for (int t = 0; t < q; ++t) {
        if (n == 0) {
            DC_CHECK(dc.group_count() == 0 && dc.edge_count() == 0 && dc.groups().empty() && dc.edges().empty());
            break;
        }
        const int r = rnd(rng, 0, 99);
        const int rest = 100 - add_pct;
        if (r < add_pct) {  // 追加（たまに自己ループ・多重辺）
            int u = rnd(rng, 0, n - 1), v = rnd(rng, 0, n - 1);
            if (!nv.alive.empty() && rnd(rng, 0, 7) == 0) {
                const auto e = nv.alive[rnd(rng, 0, (int)nv.alive.size() - 1)];
                u = e.second, v = e.first;
            }
            const bool expect = nv.labels()[u] != nv.labels()[v];
            DC_CHECK(dc.add_edge(u, v) == expect);
            nv.alive.emplace_back(u, v);
            nv.dirty = true;
        } else if (r < add_pct + rest / 2) {  // 削除（端点の順番もランダム）
            if (nv.alive.empty()) continue;
            const int k = rnd(rng, 0, (int)nv.alive.size() - 1);
            std::swap(nv.alive[k], nv.alive.back());
            const std::pair<int, int> e = nv.alive.back();
            nv.alive.pop_back();
            nv.dirty = true;
            const bool expect = nv.labels()[e.first] != nv.labels()[e.second];
            const bool got = (rng() & 1) ? dc.remove_edge(e.first, e.second) : dc.remove_edge(e.second, e.first);
            DC_CHECK(got == expect);
        } else {  // 質問
            const int x = rnd(rng, 0, n - 1), y = rnd(rng, 0, n - 1);
            const std::vector<int>& lab = nv.labels();
            const int kind = rnd(rng, 0, 11);
            if (kind <= 2) {
                DC_CHECK(dc.same(x, y) == (lab[x] == lab[y]));
            } else if (kind <= 4) {
                DC_CHECK(dc.size(x) == (int)nv.group_of(x).size());
            } else if (kind == 5) {
                DC_CHECK(dc.group_count() == nv.group_count());
            } else if (kind == 6) {
                std::vector<int> m = dc.members(x);
                std::sort(m.begin(), m.end());
                DC_CHECK(m == nv.group_of(x));
            } else if (kind == 7) {
                DC_CHECK(dc.has_edge(x, y) == (nv.count_edge(x, y) > 0));
                if (!nv.alive.empty()) {
                    const auto e = nv.alive[rnd(rng, 0, (int)nv.alive.size() - 1)];
                    DC_CHECK(dc.has_edge(e.first, e.second) && dc.has_edge(e.second, e.first));
                }
            } else if (kind == 8) {  // 値の設定
                const long long w = (long long)rnd(rng, -1000000000, 1000000000) * rnd(rng, 1, 1000);
                if (rng() & 1) dc.set_value(x, w), nv.val[x] = w;
                else dc.add_value(x, w), nv.val[x] += w;
            } else if (kind == 9) {
                DC_CHECK(dc.value(x) == nv.val[x]);
                DC_CHECK(dc.component_sum(x) == nv.sum_of(x));
            } else if (kind == 10) {
                DC_CHECK(dc.edge_count() == (int)nv.alive.size());
                if (rnd(rng, 0, 9) == 0) DC_CHECK(dc.edges() == nv.sorted_edges());
            } else if (rnd(rng, 0, 29) == 0) {  // まれに作り直し
                dc.init(n);
                nv.init(n);
            }
        }
        if (t % check_every == 0) check_all(dc, nv, rng, small);
    }
    check_all(dc, nv, rng, small);
    return q;
}

// 完全グラフなど辺の多いグラフを作ってから、全部消し切る（代わりの辺探しとレベル上げを大量に起こす）
static void dense_then_delete(std::mt19937_64& rng, int n, int order) {
    OnlineDynamicConnectivity dc(n);
    Naive nv;
    nv.init(n);
    std::vector<std::pair<int, int>> es;
    for (int u = 0; u < n; ++u)
        for (int v = u + 1; v < n; ++v) es.emplace_back(u, v);
    std::shuffle(es.begin(), es.end(), rng);
    for (const auto& e : es) {
        DC_CHECK(dc.add_edge(e.first, e.second) == (nv.labels()[e.first] != nv.labels()[e.second]));
        nv.alive.push_back(e);
        nv.dirty = true;
    }
    if (order == 1) std::reverse(es.begin(), es.end());  // 追加の逆順（木辺が最後まで残りやすい）
    if (order == 2) std::shuffle(es.begin(), es.end(), rng);
    if (order == 3) std::sort(es.begin(), es.end());  // 頂点 0 の辺（木辺になりやすい）から消す
    for (size_t k = 0; k < es.size(); ++k) {
        const auto e = es[k];
        for (size_t j = 0; j < nv.alive.size(); ++j)
            if (nv.alive[j] == e) {
                std::swap(nv.alive[j], nv.alive.back());
                nv.alive.pop_back();
                break;
            }
        nv.dirty = true;
        DC_CHECK(dc.remove_edge(e.second, e.first) == (nv.labels()[e.first] != nv.labels()[e.second]));
        DC_CHECK(dc.group_count() == nv.group_count());
        if (k % 7 == 0) check_all(dc, nv, rng, n <= 16);
    }
    DC_CHECK(dc.group_count() == n && dc.edge_count() == 0);
}

// 形の決まったケース
static void structured_cases(std::mt19937_64& rng) {
    {  // 閉路: 1 本目を消しても分かれない、2 本目で分かれる
        OnlineDynamicConnectivity dc(5);
        for (int i = 0; i < 5; ++i) DC_CHECK(dc.add_edge(i, (i + 1) % 5) == (i < 4));
        DC_CHECK(!dc.remove_edge(2, 1) && dc.same(1, 2) && dc.group_count() == 1);
        DC_CHECK(dc.remove_edge(3, 4) && !dc.same(3, 4) && dc.group_count() == 2);
        DC_CHECK(dc.size(3) == 2 && dc.size(4) == 3);
    }
    {  // 多重辺・自己ループ
        OnlineDynamicConnectivity dc(3);
        DC_CHECK(dc.add_edge(0, 1) && !dc.add_edge(1, 0) && !dc.add_edge(2, 2) && dc.edge_count() == 3);
        DC_CHECK(!dc.remove_edge(0, 1) && dc.same(0, 1) && dc.has_edge(1, 0));
        DC_CHECK(dc.remove_edge(1, 0) && !dc.same(0, 1) && !dc.has_edge(0, 1));
        DC_CHECK(!dc.remove_edge(2, 2) && dc.edge_count() == 0 && dc.group_count() == 3);
    }
    {  // 同じ辺の追加・削除を何度も
        OnlineDynamicConnectivity dc(4);
        dc.add_edge(0, 1), dc.add_edge(1, 2), dc.add_edge(2, 3), dc.add_edge(3, 0);
        for (int k = 0; k < 1000; ++k) {
            DC_CHECK(!dc.remove_edge(k % 4, (k + 1) % 4));
            DC_CHECK(dc.group_count() == 1);
            DC_CHECK(!dc.add_edge((k + 1) % 4, k % 4));
        }
    }
    for (int n : {2, 3, 4, 5, 8, 13, 20, 33})
        for (int order = 0; order < 4; ++order) dense_then_delete(rng, n, order);
    {  // グリッド: 全部つないでからランダムに消す
        const int H = 20, W = 20, n = H * W;
        OnlineDynamicConnectivity dc(n);
        Naive nv;
        nv.init(n);
        for (int i = 0; i < H; ++i)
            for (int j = 0; j < W; ++j) {
                if (i + 1 < H) dc.add_edge(i * W + j, (i + 1) * W + j), nv.alive.emplace_back(i * W + j, (i + 1) * W + j);
                if (j + 1 < W) dc.add_edge(i * W + j, i * W + j + 1), nv.alive.emplace_back(i * W + j, i * W + j + 1);
            }
        nv.dirty = true;
        std::shuffle(nv.alive.begin(), nv.alive.end(), rng);
        while (!nv.alive.empty()) {
            const auto e = nv.alive.back();
            nv.alive.pop_back();
            nv.dirty = true;
            DC_CHECK(dc.remove_edge(e.first, e.second) == (nv.labels()[e.first] != nv.labels()[e.second]));
            if (nv.alive.size() % 50 == 0) check_all(dc, nv, rng, false);
        }
    }
    {  // n = 0, 1
        OnlineDynamicConnectivity a, b(1);
        DC_CHECK(a.n() == 0 && a.group_count() == 0 && a.groups().empty());
        DC_CHECK(b.same(0, 0) && b.size(0) == 1 && !b.add_edge(0, 0) && !b.remove_edge(0, 0) && b.group_count() == 1);
    }
}

static void random_test() {
    std::mt19937_64 rng(20260923);
    long long ops = 0;
    int cases = 0;
    const int pcts[] = {40, 50, 60, 70};
    for (int it = 0; it < 300 * TEST_SCALE; ++it, ++cases) {  // 小さい n（衝突・代わりの辺探しが多い）
        const int n = rnd(rng, 0, 10), q = rnd(rng, 0, 150);
        ops += run_case(rng, n, q, pcts[it % 4]);
    }
    for (int it = 0; it < 25 * TEST_SCALE; ++it, ++cases) {  // 中くらい
        const int n = rnd(rng, 11, 100);
        ops += run_case(rng, n, 10 * n, pcts[it % 4]);
    }
    for (int it = 0; it < TEST_SCALE / 2; ++it, ++cases) {  // 大きめ（辺が多く、レベルが上がりやすい）
        const int n = rnd(rng, 300, 1000);
        ops += run_case(rng, n, 8 * n, 60 + 10 * (it % 2));
    }
    structured_cases(rng);
    std::cout << "[OnlineDynamicConnectivity] random test: OK  (" << cases << " cases, " << ops
              << " ops + dense / grid / structured cases)\n";
}

// 速度計測: 辺の追加 40%・削除 20%・質問 same 40% をランダムに合計 q 個
static void benchmark() {
    if (TEST_SCALE == 1)
        std::cout << "[OnlineDynamicConnectivity] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    const int n = 100000, q = TEST_SCALE == 1 ? 50000 : 300000;
    OnlineDynamicConnectivity dc(n);
    std::vector<std::pair<int, int>> alive;
    long long chk = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < q; ++i) {
        const int type = rnd(rng, 0, 9);
        if (type < 4 || alive.empty()) {
            const int u = rnd(rng, 0, n - 1), v = rnd(rng, 0, n - 1);
            chk += dc.add_edge(u, v);
            alive.emplace_back(u, v);
        } else if (type < 6) {
            const int k = rnd(rng, 0, (int)alive.size() - 1);
            std::swap(alive[k], alive.back());
            chk += dc.remove_edge(alive.back().first, alive.back().second);
            alive.pop_back();
        } else {
            chk += dc.same(rnd(rng, 0, n - 1), rnd(rng, 0, n - 1));
        }
    }
    const auto t1 = std::chrono::steady_clock::now();
    std::cout << "[OnlineDynamicConnectivity] benchmark: n = " << n << ", " << q
              << " ops (add 40% / remove 20% / same 40%)  "
              << (long long)std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms  (checksum " << chk
              << ")\n";
    // 辺の多いグラフで、木辺の削除（代わりの辺探し）を多く起こす
    const int n2 = TEST_SCALE == 1 ? 10000 : 50000;
    OnlineDynamicConnectivity dc2(n2);
    std::vector<std::pair<int, int>> es;
    for (int i = 0; i < 3 * n2; ++i) es.emplace_back(rnd(rng, 0, n2 - 1), rnd(rng, 0, n2 - 1));
    const auto t2 = std::chrono::steady_clock::now();
    for (const auto& e : es) dc2.add_edge(e.first, e.second);
    std::shuffle(es.begin(), es.end(), rng);
    for (const auto& e : es) chk += dc2.remove_edge(e.first, e.second);
    const auto t3 = std::chrono::steady_clock::now();
    std::cout << "[OnlineDynamicConnectivity] benchmark: n = " << n2 << ", 辺 " << 3 * n2
              << " 本を追加してから全部削除  " << (long long)std::chrono::duration<double, std::milli>(t3 - t2).count()
              << " ms  (levels " << dc2.level_count() << ")\n";
}

}  // namespace dc_test


// =============================================================================
//  使い方デモ: 全 API を 1 つずつ「返り値・内容・計算量・注意点」つきで実行して見せる。
//    出力の見方:   実際に書くコード  ->  返り値    // 何を返しているか（日本語）
//                  [状態] の行 … その時点の連結成分と辺
// =============================================================================
namespace dc_demo {

static std::string to_s(bool b) { return b ? "true" : "false"; }
static std::string to_s(int v) { return std::to_string(v); }
static std::string to_s(long long v) { return std::to_string(v); }
static std::string to_s(const std::string& s) { return s; }
static std::string to_s(const std::pair<int, int>& p) { return std::to_string(p.first) + "-" + std::to_string(p.second); }
template <class V>
static std::string to_s(const std::vector<V>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + to_s(v[i]);
    return s + "}";
}
static std::string pad(const std::string& s, size_t w) { return s.size() >= w ? s : s + std::string(w - s.size(), ' '); }

static void banner(const std::string& s) {
    std::cout << "\n=====================================================================\n"
              << "  " << s << "\n"
              << "=====================================================================\n";
}
// API の見出し: 使い方 / 返り値 / 内容 / 計算量
static void api(const std::string& usage, const std::string& ret, const std::string& what, const std::string& cost) {
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << usage << "\n"
              << "    返り値 : " << ret << "\n"
              << "    内容   : " << what << "\n"
              << "    計算量 : " << cost << "\n";
}
static void head(const std::string& s) {  // 典型パターンなどの見出し
    std::cout << "\n--------------------------------------------------------------------\n"
              << "■ " << s << "\n";
}
static void note(const std::string& s) { std::cout << "    " << s << "\n"; }  // 補足説明
static void code(const std::string& s) { std::cout << "  " << s << "\n"; }    // 実行したコード
// 文の末尾に ; を補う（for / if / while で } で終わる文には付けない）
static void stmt(const std::string& s) {
    const bool block = !s.empty() && s.back() == '}' &&
                       (s.compare(0, 4, "for ") == 0 || s.compare(0, 3, "if ") == 0 || s.compare(0, 6, "while ") == 0);
    code(block ? s : s + ";");
}
// 今の連結成分と辺を 1 行で表示
static void state(OnlineDynamicConnectivity& dc) {
    std::string s;
    for (const auto& g : dc.groups()) {
        s += "{";
        for (size_t i = 0; i < g.size(); ++i) s += (i ? "," : "") + std::to_string(g[i]);
        s += "} ";
    }
    std::string e;
    for (const auto& p : dc.edges()) e += (e.empty() ? "" : " ") + to_s(p);
    std::cout << "    [状態] 成分 " << s << " / 辺 " << (e.empty() ? "なし" : e) << "\n";
}
template <class V>
static void show(const char* expr, const V& v, const std::string& comment) {  // コード -> 返り値 // 説明
    std::string line = "  " + pad(expr, 30) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 44) + "  // " + comment;
    std::cout << line << "\n";
}
// SHOW(式, "説明") : 式そのものと返り値と説明を表示
#define SHOW(expr, comment) dc_demo::show(#expr, (expr), comment)
// RUN(文) : 文を表示してから実行（宣言にも使える）
#define RUN(...)                  \
    dc_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void cheat_sheet() {
    struct Row {
        const char* use;
        const char* ret;
        const char* what;
    };
    const Row rows[] = {
        {"OnlineDynamicConnectivity dc(n);", "", "頂点 0 ~ n-1、辺なしで作る  O(n)"},
        {"dc.add_edge(u, v)", "bool", "辺を追加。別々の成分がつながったら true  ならし O(log n)"},
        {"dc.remove_edge(u, v)", "bool", "辺を 1 本削除。成分が分かれたら true  ならし O(log^2 n)"},
        {"dc.same(u, v)", "bool", "u と v が連結なら true  ならし O(log n)"},
        {"dc.size(x)", "int", "x の連結成分の頂点数  ならし O(log n)"},
        {"dc.group_count()", "int", "連結成分の個数  O(1)"},
        {"dc.members(x)", "vector<int>", "x の連結成分の頂点（順不同）  O(成分の大きさ)"},
        {"dc.groups()", "vector<vector<int>>", "全ての連結成分（最小要素の昇順・中も昇順）  O(n log n)"},
        {"dc.has_edge(u, v)", "bool", "辺 (u, v) が今あるか  O(1)"},
        {"dc.edge_count()", "int", "今ある辺の本数  O(1)"},
        {"dc.edges()", "vector<pair<int,int>>", "今ある辺の一覧（u <= v、昇順）  O(m log m)"},
        {"dc.set_value(v, x)", "", "頂点 v の値を x にする（初期値 0）  ならし O(log n)"},
        {"dc.add_value(v, x)", "", "頂点 v の値に x を足す  ならし O(log n)"},
        {"dc.value(v)", "long long", "頂点 v の値  O(1)"},
        {"dc.component_sum(v)", "long long", "v の連結成分の値の合計  ならし O(log n)"},
        {"dc.init(n) / dc.n()", "", "作り直す / 頂点数"},
        {"dc.level_count()", "int", "（デバッグ用）内部のレベル数"},
    };
    std::cout << "\n■ API 一覧（頂点番号は 0 ~ n-1。詳しくは下の【1】~【15】）\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, 34) << pad(ret, 27) << r.what << "\n";
    }
}

static void demo() {
    banner("辺の削除もできる Union-Find（オンライン版）使い方デモ    見方:  コード  ->  返り値    // 意味");
    cheat_sheet();

    api("【1】 OnlineDynamicConnectivity dc(n);   /   OnlineDynamicConnectivity dc;（空。後から init(n)）",
        "（コンストラクタ）", "頂点 0 ~ n-1、辺なしで作る。最初は全員が別々の連結成分", "O(n)");
    RUN(OnlineDynamicConnectivity dc(6));
    state(dc);
    note("1-indexed の入力は dc.add_edge(a - 1, b - 1) のように 0-indexed に直す");

    api("【2】 dc.add_edge(u, v)", "bool … 別々の連結成分がつながったら true / もともと連結なら false",
        "辺 (u, v) を追加する。多重辺・自己ループも可（u と v の順番はどちらでもよい）", "ならし O(log n)");
    SHOW(dc.add_edge(0, 1), "{0} と {1} がつながった");
    SHOW(dc.add_edge(1, 2), "");
    SHOW(dc.add_edge(2, 0), "0-1-2 は既に連結 → false（閉路ができる）");
    SHOW(dc.add_edge(3, 4), "");
    state(dc);
    SHOW(dc.add_edge(4, 3), "多重辺も追加できる（2 本目）");
    SHOW(dc.add_edge(5, 5), "自己ループも追加できる（連結性には影響しない）");
    state(dc);

    api("【3】 dc.remove_edge(u, v)", "bool … 連結成分が分かれたら true / まだ連結なら false",
        "辺 (u, v) を 1 本削除する。存在する辺であること（u と v の順番はどちらでもよい）",
        "ならし O(log^2 n)（木辺でない辺ならならし O(log n)）");
    note("※ 存在しない辺を消すと assert で止まる。心配なら先に has_edge(u, v) で確認する");
    SHOW(dc.remove_edge(1, 0), "閉路の辺なので 0-2-1 と遠回りでまだつながっている → false");
    state(dc);
    SHOW(dc.remove_edge(2, 1), "これで 1 が切り離される → true");
    state(dc);
    SHOW(dc.remove_edge(3, 4), "多重辺の 1 本だけ消えた。もう 1 本あるのでまだ連結 → false");
    SHOW(dc.remove_edge(5, 5), "自己ループの削除 → false");
    state(dc);

    api("【4】 dc.same(u, v)", "bool", "u と v が連結なら true", "ならし O(log n)");
    SHOW(dc.same(0, 2), "0-2 の辺がある");
    SHOW(dc.same(0, 1), "1 は切り離された");
    SHOW(dc.same(3, 4), "");
    SHOW(dc.same(5, 5), "自分自身とは常に true");

    api("【5】 dc.size(x)", "int", "x の連結成分の頂点数", "ならし O(log n)");
    SHOW(dc.size(0), "{0,2}");
    SHOW(dc.size(1), "{1}");
    SHOW(dc.size(4), "{3,4}");

    api("【6】 dc.group_count()", "int", "連結成分の個数。add_edge / remove_edge が true のたびに 1 減る / 増える", "O(1)");
    SHOW(dc.group_count(), "{0,2} {1} {3,4} {5}");

    api("【7】 dc.members(x)", "vector<int> … x の連結成分の頂点の一覧", "※ 並び順は不定。昇順が必要なら sort する",
        "O(x の連結成分の大きさ)");
    SHOW(dc.members(2), "");
    RUN(std::vector<int> mem = dc.members(3));
    RUN(std::sort(mem.begin(), mem.end()));
    SHOW(mem, "sort 後");

    api("【8】 dc.groups()", "vector<vector<int>> … 全ての連結成分", "連結成分は最小要素の昇順、各成分の中も昇順に並ぶ",
        "O(n log n)  ※ ループの中で何度も呼ばない");
    SHOW(dc.groups(), "");

    api("【9】 dc.has_edge(u, v)", "bool", "辺 (u, v) が今 1 本以上あるか（u と v の順番はどちらでもよい）", "O(1)");
    SHOW(dc.has_edge(2, 0), "");
    SHOW(dc.has_edge(0, 1), "消したので false");
    SHOW(dc.has_edge(4, 3), "多重辺の 1 本が残っている");

    api("【10】 dc.edge_count()", "int", "今ある辺の本数（多重辺・自己ループも 1 本ずつ数える）", "O(1)");
    SHOW(dc.edge_count(), "0-2 と 3-4 の 2 本");

    api("【11】 dc.edges()", "vector<pair<int,int>> … 今ある辺の一覧", "各辺は (小さい頂点, 大きい頂点)。辺の一覧は昇順。多重辺はその本数だけ並ぶ",
        "O(m log m)");
    SHOW(dc.edges(), "");

    api("【12】 dc.set_value(v, x) / dc.add_value(v, x) / dc.value(v)", "set / add: なし、value: long long",
        "頂点 v の値を x にする / x を足す / 取得する。初期値は全部 0", "set / add: ならし O(log n)、value: O(1)");
    note("値を一度も設定しなければ、合計を管理する分の時間はかからない");
    RUN(dc.set_value(0, 10));
    RUN(dc.set_value(2, 20));
    RUN(dc.set_value(3, 5));
    RUN(dc.add_value(3, 100));
    SHOW(dc.value(3), "5 + 100");
    SHOW(dc.value(1), "設定していない頂点は 0");

    api("【13】 dc.component_sum(v)", "long long", "v の連結成分に含まれる頂点の値の合計（辺の追加・削除に合わせて変わる）",
        "ならし O(log n)");
    SHOW(dc.component_sum(0), "{0,2} → 10 + 20");
    SHOW(dc.component_sum(4), "{3,4} → 105 + 0");
    SHOW(dc.add_edge(2, 3), "{0,2} と {3,4} をつなぐ");
    SHOW(dc.component_sum(0), "{0,2,3,4} → 135");
    SHOW(dc.remove_edge(0, 2), "0 を切り離す");
    SHOW(dc.component_sum(4), "{2,3,4} → 125");
    SHOW(dc.component_sum(0), "{0} → 10");

    api("【14】 dc.n() / dc.init(n)", "n: int / init: なし", "頂点数 / 頂点数 n・辺なしで作り直す（値も 0 に戻る）", "O(1) / O(n)");
    SHOW(dc.n(), "");
    RUN(dc.init(3));
    state(dc);
    SHOW(dc.component_sum(0), "値も 0 に戻っている");

    api("【15】 dc.level_count()", "int", "（デバッグ用）内部で使っているレベルの数。常に floor(log2(n)) + 1 以下", "O(1)");
    SHOW(dc.level_count(), "");

    head("使い分けのまとめ: オンライン版は「その場で答えが必要なとき」に使う");
    note("全部の入力を読んでから答えてよい → オフライン版（dynamic_connectivity_offline.cpp）の方が速い（削除が多いと数倍）");
    note("削除しかない → 逆順に見て「追加だけ」にし、普通の UnionFind");

    // ================================================================== 典型パターン
    head("典型1: 強制オンライン … 入力の u, v が「前の答え（Yes = 1, No = 0）」と XOR されている問題");
    note("前の答えが分からないと次の入力が読めないので、オフライン版は使えない");
    note("クエリ: 「1 u v: 辺を追加」「2 u v: 辺を削除」「3 u v: 連結なら Yes」、u, v は 0-indexed で last と XOR");
    note("提出用のコード:");
    code("int N, Q;");
    code("std::cin >> N >> Q;");
    code("OnlineDynamicConnectivity dc(N);");
    code("int last = 0;");
    code("for (int i = 0; i < Q; ++i) {");
    code("    int t, u, v;");
    code("    std::cin >> t >> u >> v;");
    code("    u ^= last, v ^= last;");
    code("    if (t == 1) dc.add_edge(u, v);");
    code("    else if (t == 2) dc.remove_edge(u, v);");
    code("    else {");
    code("        last = dc.same(u, v) ? 1 : 0;");
    code("        std::cout << (last ? \"Yes\" : \"No\") << \"\\n\";");
    code("    }");
    code("}");
    note("↓ 同じ処理を solve(入力の文字列) にして実行した結果（答えを空白区切りで返す）");
    auto solve = [](const std::string& input) {
        std::istringstream in(input);
        int N, Q;
        in >> N >> Q;
        OnlineDynamicConnectivity d(N);
        int last = 0;
        std::string out;
        for (int i = 0; i < Q; ++i) {
            int t, u, v;
            in >> t >> u >> v;
            u ^= last, v ^= last;
            if (t == 1) d.add_edge(u, v);
            else if (t == 2) d.remove_edge(u, v);
            else {
                last = d.same(u, v) ? 1 : 0;
                out += (out.empty() ? "" : " ") + std::string(last ? "Yes" : "No");
            }
        }
        return out;
    };
    SHOW(solve("3 5\n1 0 1\n3 0 1\n2 1 0\n3 1 0\n3 0 2\n"), "");
    note("↑ 2 番目のクエリの答えが Yes なので、3 番目は 1 0 → 0 1 に戻して「0-1 を削除」と読む");

    head("典型2: 橋の判定 … 辺を 1 本消して true が返れば、その辺は橋（多重辺がない場合）。調べたら戻す");
    {
        RUN(OnlineDynamicConnectivity g(5));
        RUN(g.add_edge(0, 1));
        RUN(g.add_edge(1, 2));
        RUN(g.add_edge(2, 0));
        RUN(g.add_edge(2, 3));
        RUN(g.add_edge(3, 4));
        RUN(auto is_bridge = [&](int u, int v) { const bool b = g.remove_edge(u, v); g.add_edge(u, v); return b; });
        SHOW(is_bridge(0, 1), "三角形 0-1-2 の辺は橋ではない");
        SHOW(is_bridge(2, 3), "2-3 を消すと分かれる → 橋");
        SHOW(is_bridge(3, 4), "");
        SHOW(g.edge_count(), "調べた後は元に戻っている");
    }

    head("典型3: 頂点の値の更新と、連結成分ごとの合計（辺の追加・削除と混ざっていてもよい）");
    {
        RUN(OnlineDynamicConnectivity g(4));
        RUN(for (int v = 0; v < 4; ++v) g.set_value(v, v + 1));
        RUN(g.add_edge(0, 1));
        RUN(g.add_edge(1, 2));
        SHOW(g.component_sum(2), "1 + 2 + 3");
        RUN(g.add_value(0, 100));
        SHOW(g.component_sum(2), "頂点 0 に 100 足した");
        SHOW(g.remove_edge(0, 1), "");
        SHOW(g.component_sum(2), "0 が抜けたので 2 + 3");
    }
}

}

void solv(){
    dc_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    dc_test::benchmark();    // 速度計測
    dc_demo::demo();         // 使い方デモ
    return;
}

int main(){
    BUFF;
    int T=1;
    //cin>>T;
    while(T--)solv();
    return 0;
}