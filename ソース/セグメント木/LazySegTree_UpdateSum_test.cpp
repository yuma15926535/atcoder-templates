#pragma region template
#include <bits/stdc++.h>
#include <atcoder/all>
#include <random>
using namespace std;
using namespace atcoder;
using ll = long long;
using ull = unsigned long long;
#define rep(i,s,n) for(ll i=s;i<n;i++)//i=s,s+1...n-1 ノーマルfor
#define erep(i,s,n) for(ll i=s;i<=n;i++)//i=s,s+1...n イコールfor
#define rrep(i,s,n) for(ll i = (s)-1; i>=n;i--)//s-1スタートでnまで落ちる
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
//#define popcount __builtin_popcountll
const ll LINF = (1LL<<60);
const int INF = (1<<30);
template<class A, class B> istream& operator>>(istream& is, pair<A, B>& p) { is >> p.first >> p.second; return is; }
//template<class A, class B> ostream& operator<<(ostream& os, const pair<A, B>& p) { os << p.first << " " << p.second; return os; }
template<class T> istream& operator>>(istream& is, vector<T>& v) { for(auto& x : v) is >> x; return is; }
//template<class T> ostream& operator<<(ostream& os, const vector<T>& v) {for(int i = 0; i < (int)v.size(); i++) os << v[i] << (i + 1 == (int)v.size() ? "" : " ");return os;}
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

/*
=====================================================
★ ノード情報（② 区間代入 + 区間和）
=====================================================
*/

struct S {
    ll sum;
    int len;
};

/*
=====================================================
★ 更新情報
=====================================================
*/

using F = ll;

/*
=====================================================
★ 演算定義
=====================================================
*/

// ノード結合
S combine(S a, S b) {
    return {a.sum + b.sum, a.len + b.len};
}

// 更新適用
S apply(F f, S x) {
    if (f != LINF) x.sum = f * x.len;
    return x;
}

// 更新の合成
F compose(F f, F g) {
    return (f != LINF ? f : g);
}

S idS() { return {0, 0}; }
F idF() { return LINF; }

/*
=====================================================
Lazy SegTree 本体
=====================================================
*/

struct LazySegTree {
    int n;
    vector<S> seg;
    vector<F> lazy;

    LazySegTree(const vector<ll>& v) {
        int sz = v.size();
        n = 1;
        while (n < sz) n <<= 1;

        seg.assign(2 * n, idS());
        lazy.assign(2 * n, idF());

        for (int i = 0; i < sz; i++)
            seg[n + i] = {v[i], 1};
        for (int i = sz; i < n; i++)
            seg[n + i] = {0, 1};

        for (int i = n - 1; i > 0; i--)
            seg[i] = combine(seg[i << 1], seg[i << 1 | 1]);
    }

    void push(int k) {
        if (lazy[k] == idF()) return;
        apply_node(k << 1, lazy[k]);
        apply_node(k << 1 | 1, lazy[k]);
        lazy[k] = idF();
    }

    void apply_node(int k, F f) {
        seg[k] = apply(f, seg[k]);
        if (k < n)
            lazy[k] = compose(f, lazy[k]);
    }

    void range_apply(int l, int r, F f, int k, int a, int b) {
        if (r <= a || b <= l) return;
        if (l <= a && b <= r) {
            apply_node(k, f);
            return;
        }
        push(k);
        int mid = (a + b) / 2;
        range_apply(l, r, f, k << 1, a, mid);
        range_apply(l, r, f, k << 1 | 1, mid, b);
        seg[k] = combine(seg[k << 1], seg[k << 1 | 1]);
    }

    void range_apply(int l, int r, F f) {
        range_apply(l, r, f, 1, 0, n);
    }

    S range(int l, int r, int k, int a, int b) {
        if (r <= a || b <= l) return idS();
        if (l <= a && b <= r) return seg[k];
        push(k);
        int mid = (a + b) / 2;
        return combine(
            range(l, r, k << 1, a, mid),
            range(l, r, k << 1 | 1, mid, b)
        );
    }

    ll range_sum(int l, int r) {
        return range(l, r, 1, 0, n).sum;
    }

    // ========================
    // prefix探索
    // ========================
    int ge(ll x, int k, ll acc, int a, int b) {
        if (b - a == 1) return a;
        push(k);
        int mid = (a + b) / 2;
        if (acc + seg[k << 1].sum >= x)
            return ge(x, k << 1, acc, a, mid);
        else
            return ge(x, k << 1 | 1, acc + seg[k << 1].sum, mid, b);
    }

    int ge(ll x) {
        if (seg[1].sum < x) return n;
        return ge(x, 1, 0, 0, n);
    }

    int gt(ll x) { return ge(x + 1); }
    int le(ll x) { return ge(x + 1) - 1; }
    int lt(ll x) { return ge(x) - 1; }
};


// 検証用の愚直な実装
struct NaiveLazySegTree_UpdateSum {
    vector<ll> data;
    int n;

    NaiveLazySegTree_UpdateSum(const vector<ll>& v) : data(v) {
        n = data.size();
    }

    void range_apply(int l, int r, ll val) {
        for (int i = l; i < r; ++i) {
            data[i] = val; // 値を代入
        }
    }

    ll range_sum(int l, int r) {
        ll res = 0;
        for (int i = l; i < r; ++i) {
            res += data[i];
        }
        return res;
    }

    // 愚直な ge 関数
    int ge(ll x) {
        ll current_sum = 0;
        for (int i = 0; i < n; ++i) {
            current_sum += data[i];
            if (current_sum >= x) {
                return i;
            }
        }
        return n; // そのようなプレフィックス和が見つからない場合
    }

    int gt(ll x){ return ge(x+1); }
    int le(ll x){
        return ge(x + 1) - 1;
    }
    int lt(ll x){
        return ge(x) - 1;
    }
};

void solv(){
    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

    const int NUM_TESTS = 1000;
    const int MAX_N = 1000;
    const int MAX_VAL = 100000;
    const int NUM_QUERIES = 10000;
    const int MAX_UPDATE_VAL = 100000; // range_apply (代入) の値

    cout << "LazySegTree_UpdateSum のテストを実行中..." << endl;

    for (int t = 0; t < NUM_TESTS; ++t) {
        int n = uniform_int_distribution<int>(1, MAX_N)(rng);
        vector<ll> initial_values(n);
        for (int i = 0; i < n; ++i) {
            initial_values[i] = uniform_int_distribution<ll>(-MAX_VAL, MAX_VAL)(rng);
        }

        LazySegTree seg_tree(initial_values);
        NaiveLazySegTree_UpdateSum naive_tree(initial_values);

        for (int q = 0; q < NUM_QUERIES; ++q) {
            int type = uniform_int_distribution<int>(0, 1)(rng); // 0: apply, 1: sum query, 2: ge, 3: gt/le/lt
            int l, r;
            ll val_or_x;

            switch (type) {
                case 0: { // range_apply (代入)
                    l = uniform_int_distribution<int>(0, n - 1)(rng);
                    r = uniform_int_distribution<int>(l + 1, n)(rng);
                    val_or_x = uniform_int_distribution<ll>(-MAX_UPDATE_VAL, MAX_UPDATE_VAL)(rng);
                    seg_tree.range_apply(l, r, val_or_x);
                    naive_tree.range_apply(l, r, val_or_x);
                    break;
                }
                case 1: { // range_sum
                    l = uniform_int_distribution<int>(0, n - 1)(rng);
                    r = uniform_int_distribution<int>(l + 1, n)(rng);
                    ll seg_res = seg_tree.range_sum(l, r);
                    ll naive_res = naive_tree.range_sum(l, r);
                    if (seg_res != naive_res) {
                        cout << "テスト " << t << ", クエリ " << q << " が失敗しました (range_sum)!" << endl;
                        cout << "範囲 [" << l << ", " << r << ")" << endl;
                        cout << "SegTree 結果: " << seg_res << endl;
                        cout << "Naive 結果: " << naive_res << endl;
                        return;
                    }
                    break;
                }
            }
        }
    }

    cout << "LazySegTree_UpdateSum の全てのテストが合格しました！" << endl;
    return;
}

int main(){
    BUFF;
    int T=1;
    //cin>>T;
    while(T--)solv();
    return 0;
}