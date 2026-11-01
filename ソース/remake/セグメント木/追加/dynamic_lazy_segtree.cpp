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
//  Dynamic Lazy Segment Tree (動的遅延伝播セグメント木)
//  （作り直し版：元ファイル セグメント木/追加/dynamic_lazy_segtree.cpp）
// -----------------------------------------------------------------------------
//  説明
//    - （詳しい説明はこの下の元のコメントを参照）
//  ファイル構成
//    [テンプレート本体]  ← 「ここから」〜「ここまで」をコピーして Main_base.cpp に貼り付ける
//    それ以降            … テスト・使い方デモ（貼り付け不要）
// =============================================================================
// ======================= ここからテンプレート本体（コピー範囲） =======================
struct DynamicLazySegTree {
    // ★ 扱う座標範囲: [-2^TOP, 2^TOP)。既定の TOP = 63 で long long 全域。
    //    (座標の絶対値の上限が確実に分かっている場合だけ下げてよい。
    //     下げるとノード数も実行時間もほぼ (TOP+1)/64 倍になる。既定は全域のまま)
    static constexpr int TOP = 63;

    static constexpr ll INF = std::numeric_limits<ll>::max();  // 空区間の min として返る値

    struct Agg {
        ll sum, mn, mx;
    };

    DynamicLazySegTree() = default;
    void reserve(size_t n) { t.reserve(n); }  // 事前確保しておくと速い (任意)

    // =========================================================================
    //  区間更新 : [l, r)
    // =========================================================================
    void range_affine(ll l, ll r, ll b, ll c) {  // a[i] = b*a[i] + c
        if (l >= r) return;
        if (b == 1 && c == 0) return;  // 恒等変換
        ull ql = enc(l), qr = enc(r - 1);
        assert(in_range(ql) && in_range(qr));
        for (int s = 0; s < 2; s++) {
            Ref res = apply_rec(root_ref(s), root_lo(s), TOP, ql, qr, b, c);
            store_root(s, res);
        }
    }
    void range_add(ll l, ll r, ll v) { range_affine(l, r, 1, v); }     // a[i] += v
    void range_assign(ll l, ll r, ll v) { range_affine(l, r, 0, v); }  // a[i] = v

    // =========================================================================
    //  点操作 (LLONG_MAX / LLONG_MIN も含め全座標で使える)
    // =========================================================================
    void set(ll key, ll x) { point_affine(key, 0, x); }  // a[key] = x
    void point_add(ll key, ll v) { point_affine(key, 1, v); }  // a[key] += v

    // a[key] を 1 本道で降りて返す。ノードは生成しないし木も書き換えない。
    ll get(ll key) const {
        ull k = enc(key);
        assert(in_range(k));
        int s = side_of(k);
        ull lo = root_lo(s);
        int level = TOP;
        Ref cur = root_ref(s);
        ll b = 1, c = 0;  // 祖先から持ち越した遅延
        while (true) {
            if (cur.idx == -1) return b * cur.uni + c;  // 未生成 = 一様
            const Node& nd = t[cur.idx];
            if (level == 0) return b * nd.sum + c;  // 葉 (長さ 1 なので sum が値)
            ll nb = b * nd.b, nc = b * nd.c + c;    // 自分の遅延を合成してから降りる
            int cl = level - 1;
            ull mid = lo + (1ULL << cl);
            if (k < mid) {
                cur = Ref{nd.lc, nd.uniL};
            } else {
                cur = Ref{nd.rc, nd.uniR};
                lo = mid;
            }
            b = nb;
            c = nc;
            level = cl;
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
            Agg a = query_rec(root_ref(0), root_lo(0), TOP, ql, qr, 1, 0);
            Agg b = query_rec(root_ref(1), root_lo(1), TOP, ql, qr, 1, 0);
            return merge(a, b);
        }
        if (hit0) return query_rec(root_ref(0), root_lo(0), TOP, ql, qr, 1, 0);
        return query_rec(root_ref(1), root_lo(1), TOP, ql, qr, 1, 0);
    }
    ll range_sum(ll l, ll r) const { return prod(l, r).sum; }  // 空区間は 0
    ll range_min(ll l, ll r) const { return prod(l, r).mn; }   // 空区間は  INF
    ll range_max(ll l, ll r) const { return prod(l, r).mx; }   // 空区間は -INF

    // 全座標 (LLONG_MAX 含む) の集約。長さ 2^64 なので sum はまず溢れる点に注意。
    Agg all_prod() const {
        return merge(agg_of(root_ref(0), TOP), agg_of(root_ref(1), TOP));
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
        ll sum, mn, mx;  // この部分木の集約 (自分の遅延は反映済み)
        ll b, c;         // 子に未反映の遅延 x -> b*x + c
        ll uniL;         // lc == -1 のとき、左の子の区間はすべてこの値
        ll uniR;         // rc == -1 のとき、右の子の区間はすべてこの値
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

    static Agg uniform_agg(ll v, int level) {  // 「v が 2^level 個」
        return Agg{mul_len(v, 1ULL << level), v, v};
    }
    Agg agg_of(Ref r, int level) const {
        if (r.idx == -1) return uniform_agg(r.uni, level);
        const Node& n = t[r.idx];
        return Agg{n.sum, n.mn, n.mx};
    }

    int new_uniform(ll v, int level) {  // 「v が 2^level 個」のノードを 1 つ作る
        assert(level >= 1);
        t.push_back(Node{mul_len(v, 1ULL << level), v, v, 1, 0, v, v, -1, -1});
        return (int)t.size() - 1;
    }

    void apply_node(int p, ll b, ll c, int level) {
        t[p].sum = wadd(wmul(b, t[p].sum), mul_len(c, 1ULL << level));
        ll lo = b * t[p].mn + c, hi = b * t[p].mx + c;
        if (b < 0) std::swap(lo, hi);  // 負の傾きだと min/max が入れ替わる
        t[p].mn = lo;
        t[p].mx = hi;
        t[p].b = b * t[p].b;
        t[p].c = b * t[p].c + c;
    }

    // 遅延を子へ降ろす。★ ここでノードは一切生成しない ★
    //   ・実体のある子 : apply_node
    //   ・未生成の子   : 一様値を b*uni + c に更新するだけ
    void push(int p, int level) {
        if (t[p].b == 1 && t[p].c == 0) return;  // 押し下げるものが無い
        ll b = t[p].b, c = t[p].c;
        int cl = level - 1;
        if (t[p].lc != -1) apply_node(t[p].lc, b, c, cl);
        else t[p].uniL = b * t[p].uniL + c;
        if (t[p].rc != -1) apply_node(t[p].rc, b, c, cl);
        else t[p].uniR = b * t[p].uniR + c;
        t[p].b = 1;
        t[p].c = 0;
    }

    void pull(int p, int level) {
        int cl = level - 1;
        Agg a = agg_of(Ref{t[p].lc, t[p].uniL}, cl);
        Agg d = agg_of(Ref{t[p].rc, t[p].uniR}, cl);
        t[p].sum = wadd(a.sum, d.sum);
        t[p].mn = std::min(a.mn, d.mn);
        t[p].mx = std::max(a.mx, d.mx);
    }

    void store_left(int p, Ref r) {
        t[p].lc = r.idx;
        if (r.idx == -1) t[p].uniL = r.uni;
    }
    void store_right(int p, Ref r) {
        t[p].rc = r.idx;
        if (r.idx == -1) t[p].uniR = r.uni;
    }

    void point_affine(ll key, ll b, ll c) {
        ull k = enc(key);
        assert(in_range(k));
        int s = side_of(k);
        Ref res = apply_rec(root_ref(s), root_lo(s), TOP, k, k, b, c);
        store_root(s, res);
    }

    Ref apply_rec(Ref cur, ull lo, int level, ull ql, ull qr, ll b, ll c) {
        ull hi = lo + ((1ULL << level) - 1);
        if (qr < lo || hi < ql) return cur;  // 交差なし
        if (cur.idx == -1) {
            ll nv = b * cur.uni + c;
            if (nv == cur.uni) return cur;  // 何も変わらない (例: +0, 同じ値の assign)
            if (ql <= lo && hi <= qr) return Ref{-1, nv};  // 一様のまま。ノード生成なし
        } else if (ql <= lo && hi <= qr) {
            apply_node(cur.idx, b, c, level);
            return cur;
        }
        int p = cur.idx;
        if (p == -1) p = new_uniform(cur.uni, level);  // 分割が要るのでここだけ実体化
        push(p, level);
        int cl = level - 1;
        ull mid = lo + (1ULL << cl);
        Ref L = apply_rec(Ref{t[p].lc, t[p].uniL}, lo, cl, ql, qr, b, c);
        store_left(p, L);
        Ref R = apply_rec(Ref{t[p].rc, t[p].uniR}, mid, cl, ql, qr, b, c);
        store_right(p, R);
        pull(p, level);
        return Ref{p, 0};
    }

    // 遅延 (b, c) を引数で持ち回るので木を書き換えない。
    // 前提: [lo, hi] と [ql, qr] は交差する。
    Agg query_rec(Ref cur, ull lo, int level, ull ql, ull qr, ll b, ll c) const {
        ull hi = lo + ((1ULL << level) - 1);
        if (cur.idx == -1) {  // 一様なので交差長ぶん
            ll v = b * cur.uni + c;
            ull ilen = std::min(hi, qr) - std::max(lo, ql) + 1;
            return Agg{mul_len(v, ilen), v, v};
        }
        const Node& nd = t[cur.idx];
        if (ql <= lo && hi <= qr) {  // 丸ごと入る
            ll s = wadd(wmul(b, nd.sum), mul_len(c, 1ULL << level));
            ll mn = b * nd.mn + c, mx = b * nd.mx + c;
            if (b < 0) std::swap(mn, mx);
            return Agg{s, mn, mx};
        }
        ll nb = b * nd.b, nc = b * nd.c + c;  // 自分の遅延を合成してから子へ
        int cl = level - 1;
        ull mid = lo + (1ULL << cl);
        Ref lref{nd.lc, nd.uniL}, rref{nd.rc, nd.uniR};
        bool hitL = !(qr < lo || mid - 1 < ql);
        bool hitR = !(qr < mid || hi < ql);
        if (hitL && hitR) {
            Agg a = query_rec(lref, lo, cl, ql, qr, nb, nc);
            Agg d = query_rec(rref, mid, cl, ql, qr, nb, nc);
            return merge(a, d);
        }
        if (hitL) return query_rec(lref, lo, cl, ql, qr, nb, nc);
        return query_rec(rref, mid, cl, ql, qr, nb, nc);
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
static ll absll(ll x) { return x < 0 ? -x : x; }

// -----------------------------------------------------------------------------
//  窓テスト: 座標 base + i (i = 0..W-1) を素朴配列と突き合わせる。
//  更新は内側 [PAD, W-PAD) のみ。外側 PAD は常に 0 なので「未設定 = 0」も検証できる。
//  取得のたびに node_count が増えていないことも確認する。
// -----------------------------------------------------------------------------
static void test_window(std::mt19937_64& rng, ll base, int W, int iters) {
    const int PAD = 4;
    assert(W > 2 * PAD);
    std::vector<ll> a(W, 0);
    DynamicLazySegTree seg;
    const ll VAL = 1000;

    for (int t = 0; t < iters; t++) {
        // 値が発散しないようにガード
        ll big = 0;
        for (ll x : a) big = std::max(big, absll(x));
        if (big > 1000000000000LL) {
            ll v = rnd_range(rng, -VAL, VAL);
            for (int i = PAD; i < W - PAD; i++) a[i] = v;
            seg.range_assign(base + PAD, base + W - PAD, v);
        }

        int il = (int)rnd_range(rng, PAD, W - PAD), ir = (int)rnd_range(rng, PAD, W - PAD);
        if (il > ir) std::swap(il, ir);
        int p = (int)rnd_range(rng, PAD, W - PAD - 1);
        int type = (int)rnd_range(rng, 0, 7);

        if (type == 0) {  // range_add
            ll v = rnd_range(rng, -VAL, VAL);
            for (int i = il; i < ir; i++) a[i] += v;
            seg.range_add(base + il, base + ir, v);
        } else if (type == 1) {  // range_assign
            ll v = rnd_range(rng, -VAL, VAL);
            for (int i = il; i < ir; i++) a[i] = v;
            seg.range_assign(base + il, base + ir, v);
        } else if (type == 2) {  // range_affine
            ll b = rnd_range(rng, -3, 3), c = rnd_range(rng, -VAL, VAL);
            for (int i = il; i < ir; i++) a[i] = b * a[i] + c;
            seg.range_affine(base + il, base + ir, b, c);
        } else if (type == 3) {  // set
            ll x = rnd_range(rng, -VAL, VAL);
            a[p] = x;
            seg.set(base + p, x);
        } else if (type == 4) {  // point_add
            ll v = rnd_range(rng, -VAL, VAL);
            a[p] += v;
            seg.point_add(base + p, v);
        } else {  // 取得 (クエリ範囲は外側 PAD も含めてよい)
            size_t before = seg.node_count();
            int ql = (int)rnd_range(rng, 0, W), qr = (int)rnd_range(rng, 0, W);
            if (ql > qr) std::swap(ql, qr);
            ll s = 0, mn = DynamicLazySegTree::INF, mx = -DynamicLazySegTree::INF;
            for (int i = ql; i < qr; i++) {
                s += a[i];
                mn = std::min(mn, a[i]);
                mx = std::max(mx, a[i]);
            }
            DynamicLazySegTree::Agg q = seg.prod(base + ql, base + qr);
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
            ll s = 0, mn = DynamicLazySegTree::INF, mx = -DynamicLazySegTree::INF;
            for (int i = l; i < r; i++) {
                s += a[i];
                mn = std::min(mn, a[i]);
                mx = std::max(mx, a[i]);
            }
            DynamicLazySegTree::Agg q = seg.prod(base + l, base + r);
            assert(q.sum == s && q.mn == mn && q.mx == mx);
        }
    }
    for (int i = 0; i < W; i++) assert(seg.get(base + i) == a[i]);
    assert(seg.node_count() == before);  // ★ 全区間なめてもノードは増えない
}

// -----------------------------------------------------------------------------
//  操作の組み合わせ (add→assign, assign→affine, affine→add, set→... など)
// -----------------------------------------------------------------------------
static void test_op_pairs(std::mt19937_64& rng) {
    enum Op { ADD, ASSIGN, AFFINE, PSET, PADD, OPCNT };
    for (int o1 = 0; o1 < OPCNT; o1++) {
        for (int o2 = 0; o2 < OPCNT; o2++) {
            for (int rep = 0; rep < 40; rep++) {
                int n = (int)rnd_range(rng, 1, 12);
                ll base = rnd_range(rng, -50, 50);
                std::vector<ll> a(n);
                DynamicLazySegTree seg;
                for (int i = 0; i < n; i++) {
                    a[i] = rnd_range(rng, -5, 5);
                    seg.set(base + i, a[i]);
                }
                int ops[2] = {o1, o2};
                for (int k = 0; k < 2; k++) {
                    int l = (int)rnd_range(rng, 0, n), r = (int)rnd_range(rng, 0, n);
                    if (l > r) std::swap(l, r);
                    ll x = rnd_range(rng, -6, 6);
                    if (ops[k] == ADD) {
                        for (int i = l; i < r; i++) a[i] += x;
                        seg.range_add(base + l, base + r, x);
                    } else if (ops[k] == ASSIGN) {
                        for (int i = l; i < r; i++) a[i] = x;
                        seg.range_assign(base + l, base + r, x);
                    } else if (ops[k] == AFFINE) {
                        ll b = rnd_range(rng, -3, 3);
                        for (int i = l; i < r; i++) a[i] = b * a[i] + x;
                        seg.range_affine(base + l, base + r, b, x);
                    } else if (ops[k] == PSET) {
                        int p = (int)rnd_range(rng, 0, n - 1);
                        a[p] = x;
                        seg.set(base + p, x);
                    } else {
                        int p = (int)rnd_range(rng, 0, n - 1);
                        a[p] += x;
                        seg.point_add(base + p, x);
                    }
                    for (int i = 0; i < n; i++) assert(seg.get(base + i) == a[i]);
                    for (int l2 = 0; l2 <= n; l2++) {
                        for (int r2 = l2; r2 <= n; r2++) {
                            ll s = 0, mn = DynamicLazySegTree::INF, mx = -DynamicLazySegTree::INF;
                            for (int i = l2; i < r2; i++) {
                                s += a[i];
                                mn = std::min(mn, a[i]);
                                mx = std::max(mx, a[i]);
                            }
                            DynamicLazySegTree::Agg q = seg.prod(base + l2, base + r2);
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

    DynamicLazySegTree seg;
    for (size_t i = 0; i < keys.size(); i++) seg.set(keys[i], (ll)i + 1);
    for (size_t i = 0; i < keys.size(); i++) assert(seg.get(keys[i]) == (ll)i + 1);

    // 半開区間なので [LMIN, LMAX) は LMAX 以外の全部
    assert(seg.range_sum(LMIN, LMAX) == 78);
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
    assert(seg.range_min(5, 5) == DynamicLazySegTree::INF);
    assert(seg.range_max(5, 5) == -DynamicLazySegTree::INF);
    assert(seg.range_sum(10, 3) == 0);

    // 境界をまたぐ更新
    seg.range_add(LMIN, LMAX, 1000);
    for (size_t i = 0; i < keys.size(); i++) {
        ll expect = (keys[i] == LMAX) ? 13 : (ll)i + 1 + 1000;
        assert(seg.get(keys[i]) == expect);
    }
    seg.range_affine(LMIN + 1, LMAX - 1, -1, 0);
    assert(seg.get(LMIN) == 1001);
    assert(seg.get(LMIN + 1) == -1002);
    assert(seg.get(LMAX - 2) == -1011);
    assert(seg.get(LMAX - 1) == 1012);
    assert(seg.get(LMAX) == 13);
    // 未設定座標 (値 0) も -1 倍されているので 0 のまま
    assert(seg.get(12345) == -1000);

    // 点操作は LLONG_MAX でも使える
    seg.point_add(LMAX, 7);
    assert(seg.get(LMAX) == 20);
}

// 巨大区間に対する更新が「未生成の座標」にも効いていることを確認
static void test_huge() {
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();
    DynamicLazySegTree seg;

    seg.range_assign(LMIN, LMAX, 7);
    assert(seg.get(0) == 7);
    assert(seg.get(LMIN) == 7);
    assert(seg.get(LMAX - 1) == 7);
    assert(seg.get(LMAX) == 0);  // 半開なので含まれない
    assert(seg.range_sum(0, 10) == 70);
    assert(seg.range_min(-5, 5) == 7 && seg.range_max(-5, 5) == 7);

    seg.range_add(-1000000000000000000LL, 1000000000000000000LL, 3);
    assert(seg.get(0) == 10);
    assert(seg.get(-999999999999999999LL) == 10);
    assert(seg.get(1000000000000000000LL) == 7);  // 半開なので含まれない
    assert(seg.range_sum(-5, 5) == 100);

    seg.range_affine(-100, 100, -2, 1);  // x -> -2x + 1
    assert(seg.get(0) == -19);
    assert(seg.range_min(-100, 100) == -19 && seg.range_max(-100, 100) == -19);
    assert(seg.range_min(-101, 101) == -19 && seg.range_max(-101, 101) == 10);
    assert(seg.range_sum(-100, 100) == -19 * 200);

    seg.set(LMAX, 12345);
    assert(seg.get(LMAX) == 12345);
    seg.point_add(LMAX, 5);
    assert(seg.get(LMAX) == 12350);
    seg.set(LMIN, -1);
    assert(seg.get(LMIN) == -1);
    assert(seg.range_min(LMIN, LMIN + 1) == -1);

    assert(seg.range_sum(5, 5) == 0);
    assert(seg.range_min(5, 5) == DynamicLazySegTree::INF);
    assert(seg.range_max(5, 5) == -DynamicLazySegTree::INF);
}

// ★ 大量の取得でノードが 1 つも増えないこと (全域ランダム座標)
static void test_query_allocates_nothing(std::mt19937_64& rng) {
    DynamicLazySegTree seg;
    for (int i = 0; i < 300; i++) {
        ll l = (ll)rng(), r = (ll)rng();
        if (l > r) std::swap(l, r);
        int ty = (int)(rng() % 4);
        if (ty == 0) seg.range_add(l, r, rnd_range(rng, -1000, 1000));
        else if (ty == 1) seg.range_assign(l, r, rnd_range(rng, -1000, 1000));
        else if (ty == 2) seg.range_affine(l, r, rnd_range(rng, -2, 2), rnd_range(rng, -100, 100));
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
        0,                  // 0 起点
        -20,                // 負→正 をまたぐ (符号ビットの境界)
        LMIN,               // 最小値ぴったり
        LMIN + 3,           //
        LMAX - 60,          // 最大値ぎりぎり
        -(1LL << 62),       //
        (1LL << 62),        //
        -(1LL << 40) - 7,   //
        123456789012345LL,  //
    };
    for (ll base : bases) test_window(rng, base, 40, 1200);

    for (int rep = 0; rep < 12; rep++) {
        ll base = (ll)rng();
        if (base > LMAX - 100) base = LMAX - 100;
        if (base < LMIN + 100) base = LMIN + 100;
        test_window(rng, base, 30, 800);
    }

    test_op_pairs(rng);
    test_boundary();
    test_huge();
    test_query_allocates_nothing(rng);
    std::cout << "[dynamic_lazy_segtree] random test: OK" << std::endl;
}

// -----------------------------------------------------------------------------
//  ノード数ベンチ (A: 更新中心 / B: 取得中心 / C: 取得のみ)
// -----------------------------------------------------------------------------
static void bench_nodes() {
    std::cout << "\n----- ノード数ベンチ (座標は 64bit 全域ランダム) -----\n";
    const int Q = 200000;
    const int QA = 50000;  // A はノードを大量に作るので控えめに
    {                      // A. 更新中心
        std::mt19937_64 rng(1);
        DynamicLazySegTree seg;
        for (int i = 0; i < QA; i++) {
            ll l = (ll)rng(), r = (ll)rng();
            if (l > r) std::swap(l, r);
            int ty = (int)(rng() % 4);
            if (ty == 0) seg.range_add(l, r, (ll)(rng() % 1000));
            else if (ty == 1) seg.range_assign(l, r, (ll)(rng() % 1000));
            else if (ty == 2) seg.range_affine(l, r, 1, (ll)(rng() % 10));
            else seg.set((ll)rng(), (ll)(rng() % 1000));
        }
        std::cout << "A 更新 " << QA << " 回           : node_count = " << seg.node_count()
                  << "  (1 更新あたり " << (double)seg.node_count() / QA << " ノード, 約 "
                  << (seg.node_count() * 64 / 1048576) << " MB)\n";
    }
    {  // B. 少数更新 + 大量取得
        std::mt19937_64 rng(2);
        DynamicLazySegTree seg;
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
    {  // C. 取得のみ (少数の区間更新で遅延を残した状態から)
        std::mt19937_64 rng(3);
        DynamicLazySegTree seg;
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
    std::cout << "\n===== ⑤ Dynamic Lazy Segment Tree 使用例デモ =====\n";
    const ll LMIN = std::numeric_limits<ll>::min();
    const ll LMAX = std::numeric_limits<ll>::max();

    // 【構築】サイズも座標範囲も指定しない。未設定の座標はすべて 0。
    DynamicLazySegTree seg;
    // 【reserve】事前確保しておくと速い (任意)
    seg.reserve(1 << 12);

    // 【set】1 点代入
    seg.set(0, 10);
    seg.set(1, 20);
    seg.set(1000000000000LL, 30);  // 10^12 のような遠い座標も OK

    // 【point_add】1 点加算
    seg.point_add(1, 5);

    // 【get】未設定の座標は 0。専用の 1 本道降下なのでノードは増えない
    std::cout << "get(1)      = " << seg.get(1) << "\n";
    std::cout << "get(999)    = " << seg.get(999) << "   (未設定 → 0)\n";

    // 【range_add】[0, 5) に +100
    seg.range_add(0, 5, 100);
    std::cout << "range_add(0,5,+100) 後 get(3) = " << seg.get(3) << "\n";

    // 【range_assign】[2, 4) を 7 で塗りつぶす
    seg.range_assign(2, 4, 7);

    // 【range_affine】[0, 5) に x -> 2x + 1
    seg.range_affine(0, 5, 2, 1);

    // 【range_affine で符号反転】x -> -x  (min と max が入れ替わるケース)
    seg.range_affine(3, 5, -1, 0);

    std::cout << "座標 0..5 の値     : ";
    for (ll i = 0; i <= 5; i++) std::cout << seg.get(i) << (i == 5 ? '\n' : ' ');

    // 【range_sum / range_min / range_max】半開区間 [l, r)
    std::cout << "range_sum(0, 5)   = " << seg.range_sum(0, 5) << "\n";
    std::cout << "range_min(0, 5)   = " << seg.range_min(0, 5) << "\n";
    std::cout << "range_max(0, 5)   = " << seg.range_max(0, 5) << "\n";

    // 【prod】sum/min/max を 1 回の走査でまとめて取る (3 回呼ぶより速い)
    DynamicLazySegTree::Agg q = seg.prod(0, 10);
    std::cout << "prod(0,10) = {sum=" << q.sum << ", min=" << q.mn << ", max=" << q.mx << "}\n";

    // 【巨大区間への更新】未生成の座標にもちゃんと効く (ノードはほとんど作られない)
    seg.range_assign(-1000000000000000000LL, 1000000000000000000LL, 4);
    std::cout << "巨大 range_assign 後 get(-777777) = " << seg.get(-777777) << "\n";
    std::cout << "                     range_sum(0, 100)  = " << seg.range_sum(0, 100) << "\n";

    // 【LLONG_MAX の扱い】区間は半開なので LLONG_MAX は含まれない。点操作なら OK。
    seg.set(LMAX, 999);
    std::cout << "get(LLONG_MAX)        = " << seg.get(LMAX) << "\n";
    seg.set(LMIN, -999);
    std::cout << "get(LLONG_MIN)        = " << seg.get(LMIN) << "\n";
    std::cout << "range_min(LLONG_MIN, 0)     = " << seg.range_min(LMIN, 0) << "\n";

    // 【all_prod】全座標の集約 (長さ 2^64 なので sum は溢れる点に注意)
    DynamicLazySegTree::Agg all = seg.all_prod();
    std::cout << "all_prod()                 = {min=" << all.mn << ", max=" << all.mx << "}\n";

    // 【空区間】l >= r なら sum=0, min=INF, max=-INF
    std::cout << "range_sum(7, 7)             = " << seg.range_sum(7, 7) << "   (空区間)\n";

    // 【node_count】取得ではノードが増えないことを確認
    size_t before = seg.node_count();
    ll acc = 0;
    for (int i = 0; i < 10000; i++) acc ^= seg.range_sum(-i, i + 1) ^ seg.get(i);
    (void)acc;
    std::cout << "node_count()                = " << before << " → " << seg.node_count()
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
