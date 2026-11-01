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
//  Segment Tree  (1点更新 / 区間モノイド積)
//  （作り直し版：元ファイル セグメント木/追加/segtree.cpp）
// -----------------------------------------------------------------------------
//  説明
//    - //  ① Segment Tree  (1点更新 / 区間モノイド積)   --- AtCoder 用テンプレート
//    - // -----------------------------------------------------------------------------
//    - //  値の型      : long long 固定
//    - //  区間の流儀  : 半開区間 [l, r)  / 0-indexed
//    - //  依存        : 標準ライブラリのみ (C++17)。bits/stdc++.h も使わないので MSVC でも通る。
//    - //
//    - //  計算量
//    - //    構築            O(n)
//    - //    set       O(log n)
//    - //    get       O(1)
//    - //    prod     O(log n)
//    - //    max_right/left  O(log n)
//    - //
//    - //  モノイドの切り替え
//    - //    下の「★ 使用するモノイドをここで選ぶ ★」ブロックで、
//    - //    有効にしたい 1 行のコメントを外し、他を // でコメントアウトするだけ。
//  ファイル構成
//    [テンプレート本体]  ← 「ここから」〜「ここまで」をコピーして Main_base.cpp に貼り付ける
//    それ以降            … テスト・使い方デモ（貼り付け不要）
// =============================================================================
// =============================================================================
//  モノイド定義
//    S      : 扱う値の型 (long long 固定)
//    e()    : 単位元
//    op(a,b): 二項演算 (結合則を満たすこと)
// =============================================================================

// 区間和
// ======================= ここからテンプレート本体（コピー範囲） =======================
struct MonoidSum {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return a + b; }
};

// 区間最小値
struct MonoidMin {
    using S = ll;
    static S e() { return std::numeric_limits<ll>::max(); }
    static S op(S a, S b) { return a < b ? a : b; }
};

// 区間最大値
struct MonoidMax {
    using S = ll;
    static S e() { return std::numeric_limits<ll>::min(); }
    static S op(S a, S b) { return a > b ? a : b; }
};

// 区間 GCD  (注: 負数を渡すと絶対値の gcd になる。LLONG_MIN は使用不可)
struct MonoidGcd {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return std::gcd(a, b); }
};

// 区間 LCM  (注: 非負整数を想定。0 は吸収元。オーバーフローに注意)
struct MonoidLcm {
    using S = ll;
    static S e() { return 1; }
    static S op(S a, S b) {
        if (a == 0 || b == 0) return 0;
        S g = std::gcd(a, b);
        return (a / g) * b;  // a*b/g だとオーバーフローしやすいので先に割る
    }
};

// 区間 XOR
struct MonoidXor {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return a ^ b; }
};

// 区間 AND
struct MonoidAnd {
    using S = ll;
    static S e() { return ~0LL; }  // 全ビット 1
    static S op(S a, S b) { return a & b; }
};

// 区間 OR
struct MonoidOr {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return a | b; }
};

// =============================================================================
//  Segment Tree 本体 (非再帰)
// =============================================================================
template <class M>
struct SegTree {
    using S = typename M::S;

    SegTree() : SegTree(0) {}
    explicit SegTree(int n) : SegTree(std::vector<S>(n, M::e())) {}
    explicit SegTree(const std::vector<S>& v) : _n((int)v.size()) {
        _log = 0;
        while ((1 << _log) < _n) ++_log;
        _size = 1 << _log;
        d.assign(2 * (size_t)_size, M::e());
        for (int i = 0; i < _n; i++) d[_size + i] = v[i];
        for (int i = _size - 1; i >= 1; i--) pull(i);
    }

    int size() const { return _n; }

    // p 番目を x に書き換える
    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += _size;
        d[p] = x;
        for (int i = 1; i <= _log; i++) pull(p >> i);
    }

    // p 番目の値を取得
    S get(int p) const {
        assert(0 <= p && p < _n);
        return d[p + _size];
    }

    // [l, r) のモノイド積
    S prod(int l, int r) const {
        assert(0 <= l && l <= r && r <= _n);
        S sml = M::e(), smr = M::e();
        l += _size;
        r += _size;
        while (l < r) {
            if (l & 1) sml = M::op(sml, d[l++]);
            if (r & 1) smr = M::op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }
        return M::op(sml, smr);
    }

    // 全体のモノイド積 (= prod(0, n) だが O(1))
    S all_prod() const { return d[1]; }

    // f(prod(l, r)) が true となる最大の r を返す。
    //   f(e()) == true が前提。f は「単調」(true が続いた後 false になる) であること。
    //   例: l から右へ和を伸ばして総和が X 以下に収まる最大の r
    template <class F>
    int max_right(int l, F f) const {
        assert(0 <= l && l <= _n);
        assert(f(M::e()));
        if (l == _n) return _n;
        l += _size;
        S sm = M::e();
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(M::op(sm, d[l]))) {
                while (l < _size) {
                    l = 2 * l;
                    if (f(M::op(sm, d[l]))) {
                        sm = M::op(sm, d[l]);
                        l++;
                    }
                }
                return l - _size;
            }
            sm = M::op(sm, d[l]);
            l++;
        } while ((l & -l) != l);
        return _n;
    }

    // f(prod(l, r)) が true となる最小の l を返す (= 左へ最大限伸ばす)。
    //   f(e()) == true が前提。ACL の min_left と同じもの。
    template <class F>
    int min_left(int r, F f) const {
        assert(0 <= r && r <= _n);
        assert(f(M::e()));
        if (r == 0) return 0;
        r += _size;
        S sm = M::e();
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(M::op(d[r], sm))) {
                while (r < _size) {
                    r = 2 * r + 1;
                    if (f(M::op(d[r], sm))) {
                        sm = M::op(d[r], sm);
                        r--;
                    }
                }
                return r + 1 - _size;
            }
            sm = M::op(d[r], sm);
        } while ((r & -r) != r);
        return 0;
    }

   private:
    int _n = 0, _size = 1, _log = 0;
    std::vector<S> d;
    void pull(int k) { d[k] = M::op(d[2 * k], d[2 * k + 1]); }
};

// =============================================================================
//  ★ 使用するモノイドをここで選ぶ ★  (1 つだけ有効にする)
// =============================================================================
using M = MonoidSum;
// using M = MonoidMin;
// using M = MonoidMax;
// using M = MonoidGcd;
// using M = MonoidLcm;
// using M = MonoidXor;
// using M = MonoidAnd;
// using M = MonoidOr;

using Seg = SegTree<M>;
// ======================= ここまでテンプレート本体 =======================

// =============================================================================
//  ここから下は「ランダムテスト」と「使用例デモ」。提出時は丸ごと消してよい。
// =============================================================================

namespace tester {

// 非可換モノイド (左右の順序が保たれているかの検証用): 「最後の非ゼロ値」
struct MonoidLastNonZero {
    using S = ll;
    static S e() { return 0; }
    static S op(S a, S b) { return b != 0 ? b : a; }
};

static ll rnd_range(std::mt19937_64& rng, ll lo, ll hi) {  // [lo, hi]
    // 幅は符号なしで計算（全範囲 [LLONG_MIN, LLONG_MAX] でも溢れない）
    unsigned long long w = (unsigned long long)hi - (unsigned long long)lo + 1;
    return (ll)((unsigned long long)lo + (w ? rng() % w : rng()));
}

// LCM テスト用: 720720 (= 2^4*3^2*5*7*11*13) の約数だけを使えば lcm は 720720 以下で収まる
static std::vector<ll> divisors_of_720720() {
    std::vector<ll> d;
    for (ll x = 1; x <= 720720; x++)
        if (720720 % x == 0) d.push_back(x);
    return d;
}

// 素朴な区間積
template <class Mo>
static typename Mo::S naive_prod(const std::vector<ll>& a, int l, int r) {
    typename Mo::S s = Mo::e();
    for (int i = l; i < r; i++) s = Mo::op(s, a[i]);
    return s;
}

// 1 つのモノイドについて、set / get / prod / all_prod を検証
template <class Mo>
static void test_one_monoid(const char* name, std::mt19937_64& rng, int n, int q, ll lo, ll hi,
                            const std::vector<ll>& pool = std::vector<ll>()) {
    auto gen = [&]() -> ll { return pool.empty() ? rnd_range(rng, lo, hi) : pool[rng() % pool.size()]; };
    std::vector<ll> a(n);
    for (int i = 0; i < n; i++) a[i] = gen();
    SegTree<Mo> seg(a);

    // 構築直後の全区間チェック
    for (int l = 0; l <= n; l++)
        for (int r = l; r <= n; r++)
            assert(seg.prod(l, r) == naive_prod<Mo>(a, l, r));

    for (int t = 0; t < q; t++) {
        int type = (int)(rng() % 4);
        if (type == 0 && n > 0) {  // set
            int p = (int)(rng() % n);
            ll x = gen();
            a[p] = x;
            seg.set(p, x);
        } else if (type == 1 && n > 0) {  // get
            int p = (int)(rng() % n);
            assert(seg.get(p) == a[p]);
        } else if (type == 2) {  // prod
            int l = (int)(rng() % (n + 1)), r = (int)(rng() % (n + 1));
            if (l > r) std::swap(l, r);
            assert(seg.prod(l, r) == naive_prod<Mo>(a, l, r));
        } else {  // all_prod
            assert(seg.all_prod() == naive_prod<Mo>(a, 0, n));
        }
    }
    (void)name;
}

// max_right / min_left を素朴解と突き合わせる
template <class Mo, class F>
static void check_binsearch(const SegTree<Mo>& seg, const std::vector<ll>& a, F f) {
    int n = (int)a.size();
    for (int l = 0; l <= n; l++) {
        int expect = l;
        typename Mo::S sm = Mo::e();
        while (expect < n) {
            typename Mo::S nx = Mo::op(sm, a[expect]);
            if (!f(nx)) break;
            sm = nx;
            expect++;
        }
        assert(seg.max_right(l, f) == expect);
    }
    for (int r = 0; r <= n; r++) {
        int expect = r;
        typename Mo::S sm = Mo::e();
        while (expect > 0) {
            typename Mo::S nx = Mo::op(a[expect - 1], sm);
            if (!f(nx)) break;
            sm = nx;
            expect--;
        }
        assert(seg.min_left(r, f) == expect);
    }
}

static void test_binsearch(std::mt19937_64& rng) {
    for (int iter = 0; iter < 1500; iter++) {
        int n = (int)(rng() % 40);
        // --- Sum: 「和が X 以下」 (非負値なので単調) ---
        {
            std::vector<ll> a(n);
            for (auto& x : a) x = rnd_range(rng, 0, 20);
            SegTree<MonoidSum> seg(a);
            ll X = rnd_range(rng, 0, 200);
            check_binsearch(seg, a, [X](ll v) { return v <= X; });
        }
        // --- Min: 「最小値が X 以上」 ---
        {
            std::vector<ll> a(n);
            for (auto& x : a) x = rnd_range(rng, -20, 20);
            SegTree<MonoidMin> seg(a);
            ll X = rnd_range(rng, -25, 25);
            check_binsearch(seg, a, [X](ll v) { return v >= X; });
        }
        // --- Max: 「最大値が X 以下」 ---
        {
            std::vector<ll> a(n);
            for (auto& x : a) x = rnd_range(rng, -20, 20);
            SegTree<MonoidMax> seg(a);
            ll X = rnd_range(rng, -25, 25);
            check_binsearch(seg, a, [X](ll v) { return v <= X; });
        }
    }
}

static void random_test() {
    std::mt19937_64 rng(20260916);
    const std::vector<ll> lcm_pool = divisors_of_720720();

    // 小さい n を全網羅的に、全モノイドで (複数シード分繰り返す)
    for (int rep = 0; rep < 5; rep++) {
        for (int n = 0; n <= 17; n++) {
            test_one_monoid<MonoidSum>("Sum", rng, n, 300, -1000000000LL, 1000000000LL);
            test_one_monoid<MonoidMin>("Min", rng, n, 300, -1000000000LL, 1000000000LL);
            test_one_monoid<MonoidMax>("Max", rng, n, 300, -1000000000LL, 1000000000LL);
            test_one_monoid<MonoidGcd>("Gcd", rng, n, 300, 0, 720720);
            test_one_monoid<MonoidLcm>("Lcm", rng, n, 300, 0, 0, lcm_pool);
            test_one_monoid<MonoidXor>("Xor", rng, n, 300, 0, (1LL << 40));
            test_one_monoid<MonoidAnd>("And", rng, n, 300, -(1LL << 40), (1LL << 40));
            test_one_monoid<MonoidOr>("Or", rng, n, 300, -(1LL << 40), (1LL << 40));
            test_one_monoid<MonoidLastNonZero>("NonComm", rng, n, 300, -3, 3);
        }
    }
    // やや大きい n でも (2 のべき乗ちょうど / ±1 も含める)
    for (int rep = 0; rep < 6; rep++) {
        static const int ns[] = {63, 64, 65, 127, 128, 255};
        int n = ns[rep];
        test_one_monoid<MonoidSum>("Sum", rng, n, 4000, -1000000000000LL, 1000000000000LL);
        test_one_monoid<MonoidMin>("Min", rng, n, 4000, -1000000000000LL, 1000000000000LL);
        test_one_monoid<MonoidXor>("Xor", rng, n, 4000, 0, (1LL << 62));
        test_one_monoid<MonoidAnd>("And", rng, n, 4000, -(1LL << 62), (1LL << 62));
    }
    // 二分探索
    test_binsearch(rng);

    // 境界値 (単位元まわり)
    {
        std::vector<ll> a = {std::numeric_limits<ll>::max(), std::numeric_limits<ll>::min(), 0};
        SegTree<MonoidMin> smin(a);
        assert(smin.prod(0, 3) == std::numeric_limits<ll>::min());
        SegTree<MonoidMax> smax(a);
        assert(smax.prod(0, 3) == std::numeric_limits<ll>::max());
        SegTree<MonoidMin> s1(a);
        assert(s1.prod(1, 1) == MonoidMin::e());  // 空区間は単位元
    }
    std::cout << "[segtree] random test: OK" << std::endl;
}

}  // namespace tester

// =============================================================================
//  使用例デモ : 全機能を 1 箇所で使う。実戦ではここを見れば思い出せる。
// =============================================================================
static void demo() {
    std::cout << "\n===== ① Segment Tree 使用例デモ =====\n";

    // 【構築1】初期値の配列から構築する  O(n)
    std::vector<ll> a = {5, 3, 8, 1, 9, 2, 7, 4};
    Seg seg(a);

    // 【構築2】要素数だけ指定する。全要素が単位元 e() で初期化される。
    Seg seg_empty(8);
    std::cout << "size()            = " << seg.size() << "\n";
    std::cout << "空木の all_prod() = " << seg_empty.all_prod() << "   (= 単位元 e())\n";

    // 【set】1 点更新: a[3] を 10 にする
    seg.set(3, 10);

    // 【get】1 点取得
    std::cout << "get(3)      = " << seg.get(3) << "\n";

    // 【prod】半開区間 [2, 6) のモノイド積（今は M = MonoidSum なので区間和）
    std::cout << "prod(2, 6) = " << seg.prod(2, 6) << "\n";

    // 【prod の空区間】l == r のときは単位元が返る
    std::cout << "prod(4, 4) = " << seg.prod(4, 4) << "   (空区間 → e())\n";

    // 【all_prod】全体の積。prod(0, n) と同じだが O(1)
    std::cout << "all_prod()       = " << seg.all_prod() << "\n";

    // 【max_right】左端 l を固定し、条件を満たす最大の r を返す（f(e()) は true 必須）
    //   例: 区間和が 20 以下に収まる最大の r
    {
        SegTree<MonoidSum> s(a);
        int r = s.max_right(1, [](ll v) { return v <= 20; });
        std::cout << "max_right(1, sum<=20) = " << r
                  << "   → [1, " << r << ") の和 = " << s.prod(1, r) << "\n";
    }

    // 【min_left】右端 r を固定し、条件を満たす最小の l を返す（左へ最大限のばす）
    {
        SegTree<MonoidSum> s(a);
        int l = s.min_left(7, [](ll v) { return v <= 20; });
        std::cout << "min_left(7, sum<=20)  = " << l
                  << "   → [" << l << ", 7) の和 = " << s.prod(l, 7) << "\n";
    }

    // 【モノイドの切り替え】使いたいものを using M = ...; で選ぶだけ
    {
        std::vector<ll> b = {12, 18, 24, 6, 30};
        std::cout << "--- 同じ配列 {12,18,24,6,30} の [0,5) を各モノイドで ---\n";
        std::cout << "  Sum = " << SegTree<MonoidSum>(b).prod(0, 5) << "\n";
        std::cout << "  Min = " << SegTree<MonoidMin>(b).prod(0, 5) << "\n";
        std::cout << "  Max = " << SegTree<MonoidMax>(b).prod(0, 5) << "\n";
        std::cout << "  Gcd = " << SegTree<MonoidGcd>(b).prod(0, 5) << "\n";
        std::cout << "  Lcm = " << SegTree<MonoidLcm>(b).prod(0, 5) << "\n";
        std::cout << "  Xor = " << SegTree<MonoidXor>(b).prod(0, 5) << "\n";
        std::cout << "  And = " << SegTree<MonoidAnd>(b).prod(0, 5) << "\n";
        std::cout << "  Or  = " << SegTree<MonoidOr>(b).prod(0, 5) << "\n";
    }
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    tester::random_test();
    demo();
    return 0;
}
