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
// Pre-published templates: https://github.com/yuma15926535/atcoder-templates
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
//  Dynamic Segment Tree Beats (動的セグメント木 Beats)
//  （作り直し版：元ファイル セグメント木/追加/dynamic_segtree_beats.cpp）
// -----------------------------------------------------------------------------
//  説明
//    - （詳しい説明はこの下の元のコメントを参照）
//  ファイル構成
//    [テンプレート本体]  ← 「ここから」〜「ここまで」をコピーして Main_base.cpp に貼り付ける
//    それ以降            … テスト・使い方デモ（貼り付け不要）
// =============================================================================
// ======================= ここからテンプレート本体（コピー範囲） =======================
struct DynamicSegTreeBeats {
    // ★ 扱う座標範囲: [-2^TOP, 2^TOP)。既定の TOP = 63 で long long 全域。
    //    (座標の絶対値の上限が確実に分かっている場合だけ下げてよい。
    //     下げるとノード数も実行時間もほぼ (TOP+1)/64 倍になる。既定は全域のまま)
    static constexpr int TOP = 63;

    static constexpr ll INF = std::numeric_limits<ll>::max();  // 空区間の min として返る値

    struct Agg {
        ll sum, mn, mx;
    };

    DynamicSegTreeBeats() = default;
    void reserve(size_t n) { t.reserve(n); }  // 事前確保しておくと速い (任意)

    // =========================================================================
    //  点操作 (LLONG_MAX / LLONG_MIN も含め全座標で使える)
    // =========================================================================
    void set(ll key, ll x) {  // a[key] = x
        ull k = enc(key);
        assert(in_range(k));
        int s = side_of(k);
        Ref r = set_rec(root_ref(s), root_lo(s), TOP, k, x);
        store_root(s, r);
    }

    // a[key] を 1 本道で降りて返す。ノードは生成しないし木も書き換えない。
    ll get(ll key) const {
        ull k = enc(key);
        assert(in_range(k));
        int s = side_of(k);
        ull lo = root_lo(s);
        int level = TOP;
        Ref cur = root_ref(s);
        St st = st_of(cur, level);  // cur の「見かけ上の状態」
        while (true) {
            if (cur.idx == -1) return st.mn1;  // 未生成 = 一様なのでその値
            if (level == 0) return st.mn1;     // 葉 (長さ 1 なので mn1 == mx1 == 値)
            int cl = level - 1;
            ull mid = lo + (1ULL << cl);
            Ref nxt;
            if (k < mid) {
                nxt = Ref{t[cur.idx].lc, t[cur.idx].uniL};
            } else {
                nxt = Ref{t[cur.idx].rc, t[cur.idx].uniR};
                lo = mid;
            }
            st = st_child(st, nxt, cl);  // 仮想 push (実木は触らない)
            cur = nxt;
            level = cl;
        }
    }

    // =========================================================================
    //  区間更新 : [l, r)
    // =========================================================================
    void range_add(ll l, ll r, ll x) {  // a[i] += x
        if (l >= r || x == 0) return;
        ull ql = enc(l), qr = enc(r - 1);
        assert(in_range(ql) && in_range(qr));
        for (int s = 0; s < 2; s++) {
            Ref res = add_rec(root_ref(s), root_lo(s), TOP, ql, qr, x);
            store_root(s, res);
        }
    }
    void range_chmin(ll l, ll r, ll x) {  // a[i] = min(a[i], x)
        if (l >= r) return;
        ull ql = enc(l), qr = enc(r - 1);
        assert(in_range(ql) && in_range(qr));
        for (int s = 0; s < 2; s++) {
            Ref res = chmin_rec(root_ref(s), root_lo(s), TOP, ql, qr, x);
            store_root(s, res);
        }
    }
    void range_chmax(ll l, ll r, ll x) {  // a[i] = max(a[i], x)
        if (l >= r) return;
        ull ql = enc(l), qr = enc(r - 1);
        assert(in_range(ql) && in_range(qr));
        for (int s = 0; s < 2; s++) {
            Ref res = chmax_rec(root_ref(s), root_lo(s), TOP, ql, qr, x);
            store_root(s, res);
        }
    }

    // =========================================================================
    //  区間取得 : [l, r)   ---- すべて const / ノード生成 0 個
    // =========================================================================
    Agg prod(ll l, ll r) const {  // sum/min/max をまとめて 1 回で取る
        if (l >= r) return Agg{0, INF, -INF};
        ull ql = enc(l), qr = enc(r - 1);
        assert(in_range(ql) && in_range(qr));
        bool hit0 = (ql < HALF);   // 左の根と交差するか
        bool hit1 = (qr >= HALF);  // 右の根と交差するか
        if (hit0 && hit1) {
            Agg a = query_side(0, ql, qr);
            Agg b = query_side(1, ql, qr);
            return merge(a, b);
        }
        return hit0 ? query_side(0, ql, qr) : query_side(1, ql, qr);
    }
    ll range_sum(ll l, ll r) const { return prod(l, r).sum; }  // 空区間は 0
    ll range_min(ll l, ll r) const { return prod(l, r).mn; }   // 空区間は  INF
    ll range_max(ll l, ll r) const { return prod(l, r).mx; }   // 空区間は -INF

    // 全座標 (LLONG_MAX 含む) の集約。長さ 2^64 なので sum はまず溢れる点に注意。
    Agg all_prod() const {
        St a = st_of(root_ref(0), TOP), b = st_of(root_ref(1), TOP);
        return merge(Agg{a.sum, a.mn1, a.mx1}, Agg{b.sum, b.mn1, b.mx1});
    }

    size_t node_count() const { return t.size(); }

   private:
    static constexpr ull HALF = 1ULL << 63;
    static constexpr ull BASE = HALF - (1ULL << TOP);  // 左の根の開始位置 (TOP=63 なら 0)
    static bool in_range(ull u) { return (((u - BASE) >> TOP) >> 1) == 0; }
    static ull enc(ll x) { return (ull)x ^ HALF; }
    static int side_of(ull k) { return k < HALF ? 0 : 1; }
    static ull root_lo(int s) { return s == 0 ? BASE : HALF; }
    // 値 * 長さ を符号なしで計算 (あふれてもラップするだけで未定義動作にはならない)
    static ll mul_len(ll v, ull len) { return (ll)((ull)v * len); }
    // sum の加算・掛け算も符号なしで計算 (sum は mod 2^64 でラップする約束。未定義動作にしない)
    static ll wadd(ll a, ll b) { return (ll)((ull)a + (ull)b); }
    static ll wmul(ll a, ll b) { return (ll)((ull)a * (ull)b); }
    static Agg merge(const Agg& a, const Agg& b) {
        return Agg{wadd(a.sum, b.sum), std::min(a.mn, b.mn), std::max(a.mx, b.mx)};
    }

    struct Node {
        ll sum;
        ll mx1, mx2;  // 最大値 / 狭義 2 番目の最大値 (無ければ -INF)
        ull mxc;      // 最大値の個数
        ll mn1, mn2;  // 最小値 / 狭義 2 番目の最小値 (無ければ INF)
        ull mnc;      // 最小値の個数
        ll add;       // 子に未反映の加算
        ll uniL;      // lc == -1 のとき、左の子の区間はすべてこの値
        ll uniR;      // rc == -1 のとき、右の子の区間はすべてこの値
        int lc, rc;
    };
    std::vector<Node> t;
    int rootIdx[2] = {-1, -1};
    ll rootUni[2] = {0, 0};  // 根が未生成のときの一様値 (初期状態は全部 0)

    // 子への参照: idx == -1 なら「その区間は全部 uni」
    struct Ref {
        int idx;
        ll uni;
    };
    Ref root_ref(int s) const { return Ref{rootIdx[s], rootUni[s]}; }
    void store_root(int s, Ref r) {
        rootIdx[s] = r.idx;
        if (r.idx == -1) rootUni[s] = r.uni;
    }

    // -------------------------------------------------------------------------
    //  Beats のノード状態 (ローカルコピー用)。実ノードと同じ情報を持つ。
    //  取得系はこれを使って「仮想的に push した状態」を計算する。
    // -------------------------------------------------------------------------
    struct St {
        ll sum;
        ll mx1, mx2;
        ull mxc;
        ll mn1, mn2;
        ull mnc;
        ll add;
    };
    static St st_uniform(ll v, int level) {  // 「v が 2^level 個」
        ull len = 1ULL << level;
        return St{mul_len(v, len), v, -INF, len, v, INF, len, 0};
    }
    St st_of(Ref r, int level) const {
        if (r.idx == -1) return st_uniform(r.uni, level);
        const Node& n = t[r.idx];
        return St{n.sum, n.mx1, n.mx2, n.mxc, n.mn1, n.mn2, n.mnc, n.add};
    }
    // 以下 3 つは実ノード用の node_* と完全に同じ式 (コピー上で動かすだけ)
    static void st_add(St& s, ll x, int level) {
        s.mx1 += x;
        if (s.mx2 != -INF) s.mx2 += x;
        s.mn1 += x;
        if (s.mn2 != INF) s.mn2 += x;
        s.sum = wadd(s.sum, mul_len(x, 1ULL << level));
        s.add += x;
    }
    static void st_chmin(St& s, ll x) {  // 前提: mx2 < x < mx1
        s.sum = wadd(s.sum, mul_len(x - s.mx1, s.mxc));
        if (s.mx1 == s.mn1) { s.mx1 = s.mn1 = x; }
        else if (s.mx1 == s.mn2) { s.mx1 = s.mn2 = x; }
        else { s.mx1 = x; }
    }
    static void st_chmax(St& s, ll x) {  // 前提: mn1 < x < mn2
        s.sum = wadd(s.sum, mul_len(x - s.mn1, s.mnc));
        if (s.mn1 == s.mx1) { s.mn1 = s.mx1 = x; }
        else if (s.mn1 == s.mx2) { s.mn1 = s.mx2 = x; }
        else { s.mn1 = x; }
    }
    // 親の見かけ上の状態 parent から、子 (childRef) の「push 後の見かけ上の状態」を作る。
    // 実際の push(親) が子に対して行う操作とまったく同じ手順なので、
    // Beats の不変条件 (落とすべき値 > 2 番目の最大値) もそのまま成り立つ。
    St st_child(const St& parent, Ref childRef, int childLevel) const {
        St c = st_of(childRef, childLevel);
        if (parent.add != 0) st_add(c, parent.add, childLevel);
        if (parent.mx1 < c.mx1) st_chmin(c, parent.mx1);
        if (c.mn1 < parent.mn1) st_chmax(c, parent.mn1);
        return c;
    }

    // -------------------------------------------------------------------------
    //  取得の再帰 (const / ノード生成なし)
    //  st は cur の「見かけ上の状態」。呼び出し側が用意する。
    // -------------------------------------------------------------------------
    Agg query_side(int s, ull ql, ull qr) const {
        Ref r = root_ref(s);
        return query_rec(r, st_of(r, TOP), root_lo(s), TOP, ql, qr);
    }
    Agg query_rec(Ref cur, const St& st, ull lo, int level, ull ql, ull qr) const {
        ull hi = lo + ((1ULL << level) - 1);
        if (ql <= lo && hi <= qr) return Agg{st.sum, st.mn1, st.mx1};  // 丸ごと入る
        if (cur.idx == -1) {                                          // 一様なので交差長ぶん
            ull ilen = std::min(hi, qr) - std::max(lo, ql) + 1;
            ll v = st.mn1;  // 一様 (mn1 == mx1)
            return Agg{mul_len(v, ilen), v, v};
        }
        int cl = level - 1;
        ull mid = lo + (1ULL << cl);
        Ref lref{t[cur.idx].lc, t[cur.idx].uniL};
        Ref rref{t[cur.idx].rc, t[cur.idx].uniR};
        bool hitL = !(qr < lo || mid - 1 < ql);
        bool hitR = !(qr < mid || hi < ql);
        if (hitL && hitR) {
            Agg a = query_rec(lref, st_child(st, lref, cl), lo, cl, ql, qr);
            Agg b = query_rec(rref, st_child(st, rref, cl), mid, cl, ql, qr);
            return merge(a, b);
        }
        if (hitL) return query_rec(lref, st_child(st, lref, cl), lo, cl, ql, qr);
        return query_rec(rref, st_child(st, rref, cl), mid, cl, ql, qr);
    }

    // -------------------------------------------------------------------------
    //  更新側 (実ノードを書き換える)
    // -------------------------------------------------------------------------
    int new_uniform(ll v, int level) {  // 「v が 2^level 個」のノードを 1 つ作る
        assert(level >= 1);
        ull len = 1ULL << level;
        t.push_back(Node{mul_len(v, len), v, -INF, len, v, INF, len, 0, v, v, -1, -1});
        return (int)t.size() - 1;
    }
    void node_add(int p, ll x, int level) {
        t[p].mx1 += x;
        if (t[p].mx2 != -INF) t[p].mx2 += x;
        t[p].mn1 += x;
        if (t[p].mn2 != INF) t[p].mn2 += x;
        t[p].sum = wadd(t[p].sum, mul_len(x, 1ULL << level));
        t[p].add += x;
    }
    void node_chmin(int p, ll x) {  // 前提: mx2 < x < mx1
        t[p].sum = wadd(t[p].sum, mul_len(x - t[p].mx1, t[p].mxc));
        if (t[p].mx1 == t[p].mn1) { t[p].mx1 = t[p].mn1 = x; }
        else if (t[p].mx1 == t[p].mn2) { t[p].mx1 = t[p].mn2 = x; }
        else { t[p].mx1 = x; }
    }
    void node_chmax(int p, ll x) {  // 前提: mn1 < x < mn2
        t[p].sum = wadd(t[p].sum, mul_len(x - t[p].mn1, t[p].mnc));
        if (t[p].mn1 == t[p].mx1) { t[p].mn1 = t[p].mx1 = x; }
        else if (t[p].mn1 == t[p].mx2) { t[p].mn1 = t[p].mx2 = x; }
        else { t[p].mn1 = x; }
    }

    // 遅延を子へ降ろす。★ ここでノードは一切生成しない ★
    //   ・実体のある子 : そのまま加算 + クランプ
    //   ・未生成の子   : 一様値を clamp(uni + add, mn1, mx1) に更新するだけ
    void push(int p, int level) {
        int cl = level - 1;
        ll ad = t[p].add;
        ll pmx = t[p].mx1, pmn = t[p].mn1;
        int lc = t[p].lc, rc = t[p].rc;
        if (lc != -1) {
            if (ad != 0) node_add(lc, ad, cl);
            if (pmx < t[lc].mx1) node_chmin(lc, pmx);
            if (t[lc].mn1 < pmn) node_chmax(lc, pmn);
        } else {
            ll u = t[p].uniL + ad;
            if (u > pmx) u = pmx;
            if (u < pmn) u = pmn;
            t[p].uniL = u;
        }
        if (rc != -1) {
            if (ad != 0) node_add(rc, ad, cl);
            if (pmx < t[rc].mx1) node_chmin(rc, pmx);
            if (t[rc].mn1 < pmn) node_chmax(rc, pmn);
        } else {
            ll u = t[p].uniR + ad;
            if (u > pmx) u = pmx;
            if (u < pmn) u = pmn;
            t[p].uniR = u;
        }
        t[p].add = 0;
    }

    void pull(int p, int level) {
        int cl = level - 1;
        St L = st_of(Ref{t[p].lc, t[p].uniL}, cl);
        St R = st_of(Ref{t[p].rc, t[p].uniR}, cl);
        t[p].sum = wadd(L.sum, R.sum);
        if (L.mx1 == R.mx1) {
            t[p].mx1 = L.mx1; t[p].mxc = L.mxc + R.mxc; t[p].mx2 = std::max(L.mx2, R.mx2);
        } else if (L.mx1 > R.mx1) {
            t[p].mx1 = L.mx1; t[p].mxc = L.mxc; t[p].mx2 = std::max(L.mx2, R.mx1);
        } else {
            t[p].mx1 = R.mx1; t[p].mxc = R.mxc; t[p].mx2 = std::max(R.mx2, L.mx1);
        }
        if (L.mn1 == R.mn1) {
            t[p].mn1 = L.mn1; t[p].mnc = L.mnc + R.mnc; t[p].mn2 = std::min(L.mn2, R.mn2);
        } else if (L.mn1 < R.mn1) {
            t[p].mn1 = L.mn1; t[p].mnc = L.mnc; t[p].mn2 = std::min(L.mn2, R.mn1);
        } else {
            t[p].mn1 = R.mn1; t[p].mnc = R.mnc; t[p].mn2 = std::min(R.mn2, L.mn1);
        }
    }

    // 子の結果を書き戻すヘルパ (vector の再確保に備えて必ず呼び出し後に代入する)
    void store_left(int p, Ref r) {
        t[p].lc = r.idx;
        if (r.idx == -1) t[p].uniL = r.uni;
    }
    void store_right(int p, Ref r) {
        t[p].rc = r.idx;
        if (r.idx == -1) t[p].uniR = r.uni;
    }

    Ref add_rec(Ref cur, ull lo, int level, ull ql, ull qr, ll x) {
        ull hi = lo + ((1ULL << level) - 1);
        if (qr < lo || hi < ql) return cur;  // 交差なし
        if (ql <= lo && hi <= qr) {          // 丸ごと入る
            if (cur.idx == -1) return Ref{-1, cur.uni + x};  // 一様のまま。ノード生成なし
            node_add(cur.idx, x, level);
            return cur;
        }
        int p = cur.idx;
        if (p == -1) p = new_uniform(cur.uni, level);  // 分割が要るのでここだけ実体化
        push(p, level);
        int cl = level - 1;
        ull mid = lo + (1ULL << cl);
        Ref L = add_rec(Ref{t[p].lc, t[p].uniL}, lo, cl, ql, qr, x);
        store_left(p, L);
        Ref R = add_rec(Ref{t[p].rc, t[p].uniR}, mid, cl, ql, qr, x);
        store_right(p, R);
        pull(p, level);
        return Ref{p, 0};
    }

    Ref chmin_rec(Ref cur, ull lo, int level, ull ql, ull qr, ll x) {
        ull hi = lo + ((1ULL << level) - 1);
        if (qr < lo || hi < ql) return cur;
        bool covered = (ql <= lo && hi <= qr);
        if (cur.idx == -1) {
            if (cur.uni <= x) return cur;             // 変化なし。ノード生成なし
            if (covered) return Ref{-1, x};           // 一様のまま。ノード生成なし
        } else {
            if (t[cur.idx].mx1 <= x) return cur;      // 触る必要なし
            if (covered && t[cur.idx].mx2 < x) {      // 最大値だけ潰せば済む (Beats の打ち切り)
                node_chmin(cur.idx, x);
                return cur;
            }
        }
        int p = cur.idx;
        if (p == -1) p = new_uniform(cur.uni, level);
        push(p, level);
        int cl = level - 1;
        ull mid = lo + (1ULL << cl);
        Ref L = chmin_rec(Ref{t[p].lc, t[p].uniL}, lo, cl, ql, qr, x);
        store_left(p, L);
        Ref R = chmin_rec(Ref{t[p].rc, t[p].uniR}, mid, cl, ql, qr, x);
        store_right(p, R);
        pull(p, level);
        return Ref{p, 0};
    }

    Ref chmax_rec(Ref cur, ull lo, int level, ull ql, ull qr, ll x) {
        ull hi = lo + ((1ULL << level) - 1);
        if (qr < lo || hi < ql) return cur;
        bool covered = (ql <= lo && hi <= qr);
        if (cur.idx == -1) {
            if (x <= cur.uni) return cur;
            if (covered) return Ref{-1, x};
        } else {
            if (x <= t[cur.idx].mn1) return cur;
            if (covered && x < t[cur.idx].mn2) {
                node_chmax(cur.idx, x);
                return cur;
            }
        }
        int p = cur.idx;
        if (p == -1) p = new_uniform(cur.uni, level);
        push(p, level);
        int cl = level - 1;
        ull mid = lo + (1ULL << cl);
        Ref L = chmax_rec(Ref{t[p].lc, t[p].uniL}, lo, cl, ql, qr, x);
        store_left(p, L);
        Ref R = chmax_rec(Ref{t[p].rc, t[p].uniR}, mid, cl, ql, qr, x);
        store_right(p, R);
        pull(p, level);
        return Ref{p, 0};
    }

    // 1 点代入。必要な側の子だけを辿る (もう片方は一様値スロットのまま)
    Ref set_rec(Ref cur, ull lo, int level, ull k, ll x) {
        if (level == 0) {
            assert(cur.idx == -1);  // 葉にノードは作らない
            return Ref{-1, x};
        }
        if (cur.idx == -1 && cur.uni == x) return cur;  // 変化なし。ノード生成なし
        int p = cur.idx;
        if (p == -1) p = new_uniform(cur.uni, level);
        push(p, level);
        int cl = level - 1;
        ull mid = lo + (1ULL << cl);
        if (k < mid) {
            Ref L = set_rec(Ref{t[p].lc, t[p].uniL}, lo, cl, k, x);
            store_left(p, L);
        } else {
            Ref R = set_rec(Ref{t[p].rc, t[p].uniR}, mid, cl, k, x);
            store_right(p, R);
        }
        pull(p, level);
        return Ref{p, 0};
    }
};
// ======================= ここまでテンプレート本体 =======================

// =============================================================================
//  ここから下は「ランダムテスト」と「使用例デモ」。提出時は丸ごと消してよい。
// =============================================================================

namespace tester {

static ll rnd_range(std::mt19937_64& rng, ll lo, ll hi) {
    // 幅は符号なしで計算（全範囲 [LLONG_MIN, LLONG_MAX] でも溢れない）
    unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return (ll)((unsigned long long)lo + (w ? rng() % w : rng()));
}

// -----------------------------------------------------------------------------
//  窓テスト: 座標 base + i (i = 0..W-1) を素朴配列と突き合わせる。
//  更新は内側 [PAD, W-PAD) のみ。外側 PAD は常に 0 なので「未設定 = 0」も検証できる。
//  取得のたびに node_count が増えていないことも確認する。
// -----------------------------------------------------------------------------
static void test_window(std::mt19937_64& rng, ll base, int W, int iters, ll VAL) {
    const int PAD = 4;
    assert(W > 2 * PAD);
    std::vector<ll> a(W, 0);
    DynamicSegTreeBeats seg;

    for (int t = 0; t < iters; t++) {
        int il = (int)rnd_range(rng, PAD, W - PAD), ir = (int)rnd_range(rng, PAD, W - PAD);
        if (il > ir) std::swap(il, ir);
        int p = (int)rnd_range(rng, PAD, W - PAD - 1);
        int type = (int)rnd_range(rng, 0, 7);
        ll x = rnd_range(rng, -VAL, VAL);

        if (type == 0) {  // set
            a[p] = x;
            seg.set(base + p, x);
        } else if (type == 1) {  // range_add
            ll v = rnd_range(rng, -VAL / 2, VAL / 2);
            for (int i = il; i < ir; i++) a[i] += v;
            seg.range_add(base + il, base + ir, v);
        } else if (type == 2 || type == 3) {  // range_chmin
            for (int i = il; i < ir; i++) a[i] = std::min(a[i], x);
            seg.range_chmin(base + il, base + ir, x);
        } else if (type == 4 || type == 5) {  // range_chmax
            for (int i = il; i < ir; i++) a[i] = std::max(a[i], x);
            seg.range_chmax(base + il, base + ir, x);
        } else {  // 取得 (クエリ範囲は外側 PAD も含めてよい)
            size_t before = seg.node_count();
            int ql = (int)rnd_range(rng, 0, W), qr = (int)rnd_range(rng, 0, W);
            if (ql > qr) std::swap(ql, qr);
            ll s = 0, mn = DynamicSegTreeBeats::INF, mx = -DynamicSegTreeBeats::INF;
            for (int i = ql; i < qr; i++) {
                s += a[i];
                mn = std::min(mn, a[i]);
                mx = std::max(mx, a[i]);
            }
            DynamicSegTreeBeats::Agg q = seg.prod(base + ql, base + qr);
            assert(q.sum == s && q.mn == mn && q.mx == mx);
            assert(seg.range_sum(base + ql, base + qr) == s);
            assert(seg.range_min(base + ql, base + qr) == mn);
            assert(seg.range_max(base + ql, base + qr) == mx);
            assert(seg.node_count() == before);  // ★ 取得でノードが増えないこと
        }
        size_t before = seg.node_count();
        int gp = (int)rnd_range(rng, 0, W - 1);
        assert(seg.get(base + gp) == a[gp]);
        assert(seg.node_count() == before);  // ★ get でノードが増えないこと
    }

    // 最後に全区間を突き合わせ
    size_t before = seg.node_count();
    for (int l = 0; l <= W; l++) {
        for (int r = l; r <= W; r++) {
            ll s = 0, mn = DynamicSegTreeBeats::INF, mx = -DynamicSegTreeBeats::INF;
            for (int i = l; i < r; i++) {
                s += a[i];
                mn = std::min(mn, a[i]);
                mx = std::max(mx, a[i]);
            }
            DynamicSegTreeBeats::Agg q = seg.prod(base + l, base + r);
            assert(q.sum == s && q.mn == mn && q.mx == mx);
        }
    }
    for (int i = 0; i < W; i++) assert(seg.get(base + i) == a[i]);
    assert(seg.node_count() == before);  // ★ 全区間なめてもノードは増えない
}

// -----------------------------------------------------------------------------
//  操作の組み合わせ (add→chmin, chmin→add, chmax→chmin, set→... など) を
//  重点的に確認する。1 回ごとに全要素と全区間を突き合わせる。
// -----------------------------------------------------------------------------
static void test_op_pairs(std::mt19937_64& rng) {
    enum Op { ADD, CHMIN, CHMAX, PSET, OPCNT };
    const char* names[] = {"add", "chmin", "chmax", "set"};
    (void)names;
    for (int o1 = 0; o1 < OPCNT; o1++) {
        for (int o2 = 0; o2 < OPCNT; o2++) {
            for (int rep = 0; rep < 60; rep++) {
                int n = (int)rnd_range(rng, 1, 12);
                ll base = rnd_range(rng, -50, 50);
                std::vector<ll> a(n);
                DynamicSegTreeBeats seg;
                for (int i = 0; i < n; i++) {
                    a[i] = rnd_range(rng, -4, 4);
                    seg.set(base + i, a[i]);
                }
                int ops[2] = {o1, o2};
                for (int k = 0; k < 2; k++) {
                    int l = (int)rnd_range(rng, 0, n), r = (int)rnd_range(rng, 0, n);
                    if (l > r) std::swap(l, r);
                    ll x = rnd_range(rng, -5, 5);
                    if (ops[k] == ADD) {
                        for (int i = l; i < r; i++) a[i] += x;
                        seg.range_add(base + l, base + r, x);
                    } else if (ops[k] == CHMIN) {
                        for (int i = l; i < r; i++) a[i] = std::min(a[i], x);
                        seg.range_chmin(base + l, base + r, x);
                    } else if (ops[k] == CHMAX) {
                        for (int i = l; i < r; i++) a[i] = std::max(a[i], x);
                        seg.range_chmax(base + l, base + r, x);
                    } else {
                        int p = (int)rnd_range(rng, 0, n - 1);
                        a[p] = x;
                        seg.set(base + p, x);
                    }
                    // 毎回フル検証
                    for (int i = 0; i < n; i++) assert(seg.get(base + i) == a[i]);
                    for (int l2 = 0; l2 <= n; l2++) {
                        for (int r2 = l2; r2 <= n; r2++) {
                            ll s = 0, mn = DynamicSegTreeBeats::INF, mx = -DynamicSegTreeBeats::INF;
                            for (int i = l2; i < r2; i++) {
                                s += a[i];
                                mn = std::min(mn, a[i]);
                                mx = std::max(mx, a[i]);
                            }
                            DynamicSegTreeBeats::Agg q = seg.prod(base + l2, base + r2);
                            assert(q.sum == s && q.mn == mn && q.mx == mx);
                        }
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
//  境界値テスト: LLONG_MIN / LLONG_MAX 近傍と符号のまたぎ
// -----------------------------------------------------------------------------
static void test_boundary() {
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();
    std::vector<ll> keys = {LMIN, LMIN + 1, LMIN + 2, -3, -2, -1, 0, 1, 2, 3, LMAX - 2, LMAX - 1, LMAX};

    DynamicSegTreeBeats seg;
    for (size_t i = 0; i < keys.size(); i++) seg.set(keys[i], (ll)i + 1);
    for (size_t i = 0; i < keys.size(); i++) assert(seg.get(keys[i]) == (ll)i + 1);

    // 半開区間なので [LMIN, LMAX) は LMAX 以外の全部
    assert(seg.range_sum(LMIN, LMAX) == 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10 + 11 + 12);
    assert(seg.range_max(LMIN, LMAX) == 12);
    assert(seg.get(LMAX) == 13);

    // 1 要素だけの区間
    assert(seg.range_sum(LMIN, LMIN + 1) == 1);
    assert(seg.range_sum(LMAX - 1, LMAX) == 12);
    assert(seg.range_sum(0, 1) == 7);

    // 符号をまたぐ区間
    assert(seg.range_sum(-3, 4) == 4 + 5 + 6 + 7 + 8 + 9 + 10);
    assert(seg.range_min(-3, 4) == 4 && seg.range_max(-3, 4) == 10);

    // 空区間 / 逆転区間
    assert(seg.range_sum(5, 5) == 0);
    assert(seg.range_min(5, 5) == DynamicSegTreeBeats::INF);
    assert(seg.range_max(5, 5) == -DynamicSegTreeBeats::INF);
    assert(seg.range_sum(10, 3) == 0);

    // 境界をまたぐ chmin / chmax / add
    seg.range_chmin(LMIN, LMAX, 6);
    for (size_t i = 0; i < keys.size(); i++) {
        ll expect = (keys[i] == LMAX) ? 13 : std::min((ll)i + 1, 6LL);
        assert(seg.get(keys[i]) == expect);
    }
    seg.range_chmax(LMIN, LMAX, 3);
    for (size_t i = 0; i < keys.size(); i++) {
        ll expect = (keys[i] == LMAX) ? 13 : std::max(std::min((ll)i + 1, 6LL), 3LL);
        assert(seg.get(keys[i]) == expect);
    }
    seg.range_add(LMIN + 1, LMAX - 1, 100);
    assert(seg.get(LMIN) == 3);
    assert(seg.get(LMIN + 1) == 103);
    assert(seg.get(LMAX - 2) == 106);
    assert(seg.get(LMAX - 1) == 6);
    assert(seg.get(LMAX) == 13);
}

// 巨大区間に対する更新が「未生成の座標」にも効いていることを確認
static void test_huge() {
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();
    DynamicSegTreeBeats seg;

    // 初期値は全部 0。chmin(-5) で LLONG_MAX 以外が -5 になる
    seg.range_chmin(LMIN, LMAX, -5);
    assert(seg.get(0) == -5);
    assert(seg.get(LMIN) == -5);
    assert(seg.get(LMAX - 1) == -5);
    assert(seg.get(LMAX) == 0);  // 半開なので含まれない
    assert(seg.range_sum(0, 10) == -50);

    seg.range_chmax(-100, 100, 3);
    assert(seg.get(0) == 3);
    assert(seg.get(-101) == -5);
    assert(seg.range_sum(-100, 100) == 3 * 200);
    assert(seg.range_min(-101, 101) == -5 && seg.range_max(-101, 101) == 3);

    seg.range_add(-10, 10, 100);
    assert(seg.get(0) == 103);
    assert(seg.get(50) == 3);
    assert(seg.range_sum(-10, 10) == 103 * 20);

    seg.range_chmin(-1000000000000000000LL, 1000000000000000000LL, 50);
    assert(seg.get(0) == 50);
    assert(seg.get(-50) == 3);
    assert(seg.get(-10000000) == -5);

    // 点操作は LLONG_MAX でも使える
    seg.set(LMAX, 777);
    assert(seg.get(LMAX) == 777);
    seg.set(LMIN, -777);
    assert(seg.get(LMIN) == -777);
    assert(seg.range_min(LMIN, LMIN + 1) == -777);

    // 空区間
    assert(seg.range_sum(5, 5) == 0);
    assert(seg.range_min(5, 5) == DynamicSegTreeBeats::INF);
    assert(seg.range_max(5, 5) == -DynamicSegTreeBeats::INF);
}

// ★ 最重要: 大量の取得でノードが 1 つも増えないこと (全域ランダム座標)
static void test_query_allocates_nothing(std::mt19937_64& rng) {
    DynamicSegTreeBeats seg;
    for (int i = 0; i < 300; i++) {
        ll l = (ll)rng(), r = (ll)rng();
        if (l > r) std::swap(l, r);
        int ty = (int)(rng() % 4);
        if (ty == 0) seg.range_add(l, r, rnd_range(rng, -1000, 1000));
        else if (ty == 1) seg.range_chmin(l, r, rnd_range(rng, -1000, 1000));
        else if (ty == 2) seg.range_chmax(l, r, rnd_range(rng, -1000, 1000));
        else seg.set((ll)rng(), rnd_range(rng, -1000, 1000));
    }
    size_t before = seg.node_count();
    ll acc = 0;
    for (int i = 0; i < 50000; i++) {
        ll l = (ll)rng(), r = (ll)rng();
        if (l > r) std::swap(l, r);
        acc ^= seg.range_sum(l, r) ^ seg.range_min(l, r) ^ seg.range_max(l, r);
        acc ^= seg.get((ll)rng());
        acc ^= seg.all_prod().sum;
    }
    (void)acc;
    assert(seg.node_count() == before);
}

static void random_test() {
    std::mt19937_64 rng(20260917ULL);
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();

    std::vector<ll> bases = {
        0, -20, LMIN, LMIN + 3, LMAX - 60, -(1LL << 62), (1LL << 62), -(1LL << 40) - 7, 987654321098765LL,
    };
    for (ll base : bases) {
        test_window(rng, base, 40, 1200, 1000);  // 値がばらけるケース
        test_window(rng, base, 24, 1200, 3);     // 値の重複が多いケース (Beats の分岐)
    }
    for (int rep = 0; rep < 8; rep++) {
        ll base = (ll)rng();
        if (base > LMAX - 100) base = LMAX - 100;
        if (base < LMIN + 100) base = LMIN + 100;
        test_window(rng, base, 30, 800, 5);
    }

    test_op_pairs(rng);
    test_boundary();
    test_huge();
    test_query_allocates_nothing(rng);
    std::cout << "[DynamicSegmentTreeBeats] random test: OK" << std::endl;
}

// -----------------------------------------------------------------------------
//  ノード数ベンチ (A: 更新中心 / B: 取得中心 / C: 取得のみ)
// -----------------------------------------------------------------------------
static void bench_nodes() {
    std::cout << "\n----- ノード数ベンチ (座標は 64bit 全域ランダム) -----\n";
    const int Q = 200000;
    const int QA = 50000;  // A はノードを大量に作るので控えめに (1 操作あたりの平均を出す)
    {  // A. 更新中心
        std::mt19937_64 rng(1);
        DynamicSegTreeBeats seg;
        for (int i = 0; i < QA; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            int ty = (int)(rng() % 4);
            if (ty == 0) seg.range_add(l, r, (ll)(rng() % 1000));
            else if (ty == 1) seg.range_chmin(l, r, (ll)(rng() % 1000));
            else if (ty == 2) seg.range_chmax(l, r, (ll)(rng() % 1000) - 500);
            else seg.set((ll)rng(), (ll)(rng() % 1000));
        }
        std::cout << "A 更新 " << QA << " 回           : node_count = " << seg.node_count()
                  << "  (1 更新あたり " << (double)seg.node_count() / QA << " ノード, 約 "
                  << (seg.node_count() * 88 / 1048576) << " MB)\n";
    }
    {  // B. 少数更新 + 大量取得
        std::mt19937_64 rng(2);
        DynamicSegTreeBeats seg;
        for (int i = 0; i < 1000; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            seg.range_add(l, r, (ll)(rng() % 1000));
        }
        size_t before = seg.node_count();
        ll acc = 0;
        for (int i = 0; i < Q; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            acc ^= seg.range_sum(l, r) ^ seg.range_min(l, r) ^ seg.get((ll)rng());
        }
        (void)acc;
        std::cout << "B 更新 1000 + 取得 " << Q << " 回 : node_count = " << before << " → "
                  << seg.node_count() << "  (増加 " << seg.node_count() - before << ")\n";
    }
    {  // C. 取得のみ (最重要) ... 少数の区間更新で遅延を残した状態から
        std::mt19937_64 rng(3);
        DynamicSegTreeBeats seg;
        for (int i = 0; i < 50; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            seg.range_add(l, r, (ll)(rng() % 1000));
        }
        size_t before = seg.node_count();
        ll acc = 0;
        for (int i = 0; i < Q; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            acc ^= seg.range_sum(l, r) ^ seg.range_min(l, r) ^ seg.range_max(l, r);
            acc ^= seg.get((ll)rng());
        }
        (void)acc;
        std::cout << "C 取得のみ " << Q << " 回      : node_count = " << before << " → "
                  << seg.node_count() << "  (増加 " << seg.node_count() - before << ")\n";
    }
}

}  // namespace tester

// =============================================================================
//  使用例デモ : 全機能を 1 箇所で使う。実戦ではここを見れば思い出せる。
// =============================================================================
static void demo() {
    std::cout << "\n===== ⑥ Dynamic Segment Tree Beats 使用例デモ =====\n";
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();

    // 【構築】サイズも座標範囲も指定しない。未設定の座標はすべて 0。
    DynamicSegTreeBeats seg;
    // 【reserve】事前確保しておくと速い (任意)
    seg.reserve(1 << 12);

    // 【set】1 点代入。遠い座標も OK
    for (ll i = 0; i < 8; i++) seg.set(i, (i * 37) % 11);
    seg.set(1000000000000LL, 100);

    // 【get】未設定の座標は 0。専用の 1 本道降下なのでノードは増えない
    std::cout << "座標 0..7 の値      : ";
    for (ll i = 0; i < 8; i++) std::cout << seg.get(i) << (i == 7 ? '\n' : ' ');
    std::cout << "get(999)      = " << seg.get(999) << "   (未設定 → 0)\n";

    // 【range_chmin】[0, 8) を min(a[i], 6) に
    seg.range_chmin(0, 8, 6);
    std::cout << "range_chmin(0,8,6)  : ";
    for (ll i = 0; i < 8; i++) std::cout << seg.get(i) << (i == 7 ? '\n' : ' ');

    // 【range_chmax】[2, 6) を max(a[i], 4) に
    seg.range_chmax(2, 6, 4);
    std::cout << "range_chmax(2,6,4)  : ";
    for (ll i = 0; i < 8; i++) std::cout << seg.get(i) << (i == 7 ? '\n' : ' ');

    // 【range_add】[1, 5) に +10
    seg.range_add(1, 5, 10);
    std::cout << "range_add(1,5,+10)  : ";
    for (ll i = 0; i < 8; i++) std::cout << seg.get(i) << (i == 7 ? '\n' : ' ');

    // 【range_sum / range_min / range_max】半開区間 [l, r)
    std::cout << "range_sum(0, 8)     = " << seg.range_sum(0, 8) << "\n";
    std::cout << "range_min(0, 8)     = " << seg.range_min(0, 8) << "\n";
    std::cout << "range_max(0, 8)     = " << seg.range_max(0, 8) << "\n";

    // 【prod】sum/min/max を 1 回の走査でまとめて取る (3 回呼ぶより速い)
    DynamicSegTreeBeats::Agg q = seg.prod(2, 6);
    std::cout << "prod(2,6)    = {sum=" << q.sum << ", min=" << q.mn << ", max=" << q.mx << "}\n";

    // 【chmin + chmax の合わせ技】区間を [8, 12] にクランプ
    seg.range_chmax(0, 8, 8);
    seg.range_chmin(0, 8, 12);
    std::cout << "[8,12] にクランプ   : ";
    for (ll i = 0; i < 8; i++) std::cout << seg.get(i) << (i == 7 ? '\n' : ' ');

    // 【巨大区間への更新】未生成の座標にもちゃんと効く (ノードはほとんど作られない)
    seg.range_chmax(-1000000000000000000LL, 1000000000000000000LL, 1);
    std::cout << "巨大 range_chmax(...,1) 後 get(-777777) = " << seg.get(-777777) << "\n";

    // 【LLONG_MAX の扱い】区間は半開なので含まれない。点操作なら OK。
    seg.set(LMAX, 55);
    std::cout << "get(LLONG_MAX)  = " << seg.get(LMAX) << "\n";
    seg.set(LMIN, -55);
    std::cout << "get(LLONG_MIN)  = " << seg.get(LMIN) << "\n";
    std::cout << "range_min(LLONG_MIN,0)= " << seg.range_min(LMIN, 0) << "\n";

    // 【all_prod】全座標の集約 (長さ 2^64 なので sum は溢れる点に注意)
    DynamicSegTreeBeats::Agg all = seg.all_prod();
    std::cout << "all_prod()           = {min=" << all.mn << ", max=" << all.mx << "}\n";

    // 【空区間】l >= r なら sum=0, min=INF, max=-INF
    std::cout << "range_sum(3, 3)       = " << seg.range_sum(3, 3) << "   (空区間)\n";

    // 【node_count】取得ではノードが増えないことを確認
    size_t before = seg.node_count();
    ll acc = 0;
    for (int i = 0; i < 10000; i++) acc ^= seg.range_sum(-i, i + 1) ^ seg.get(i);
    (void)acc;
    std::cout << "node_count()          = " << before << " → " << seg.node_count()
              << "  (取得 1 万回してもノードは増えない)\n";
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    tester::random_test();
    tester::bench_nodes();
    demo();
    return 0;
}
