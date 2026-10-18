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
//  辺の削除もできる Union-Find（オフライン版 / Offline Dynamic Connectivity）  --- AtCoder 用テンプレート
// -----------------------------------------------------------------------------
//  「辺を足す」「辺を消す」「u と v は連結か？」などが混ざった操作列を、
//  先に全部登録しておき、最後に run() でまとめて答えを求める。
//
//  ★ どれを使うか
//    - 辺の削除があり、全部の入力を読んでから答えてよい   → このファイル（オフライン版。速い）
//    - 前の答えで次の入力が決まる（強制オンライン）など、その場で答える必要がある
//                                                          → dynamic_connectivity_online.cpp
//    - 削除しかない（最初に辺があって消えていくだけ）     → 逆順に見て「追加だけ」にし、普通の UnionFind
//
//  仕組み
//    各辺が存在する期間は「質問番号の区間 [l, r)」になる。これを時間軸のセグメント木の
//    O(log Q) 個のノードに載せ、根から DFS する。ノードに入るときにその辺を unite し、
//    出るときに巻き戻す（RollbackUnionFind）。葉に着いたときの状態 = その質問の時点のグラフ。
//
//  頂点番号   : 0-indexed（0 ~ n-1）
//  依存       : 標準ライブラリのみ（C++14 以上）。bits/stdc++.h 不要なので MSVC でも通る
//  文字コード : UTF-8 (BOM 付き)。Windows のコンソールでも文字化けしないよう main で自動設定
//
//  高速化のポイント
//    - 辺の追加と削除の対応付けは sort 1 回（ハッシュ表なし）
//    - セグメント木の各ノードの辺は CSR 形式（1 本の配列）に詰めて持つ（vector の vector を作らない）
//    - 質問のない時間帯（右端の余り）は DFS しない
//    - RollbackUnionFind は union by size で木の高さ log2(n) 以下（find は経路圧縮なしで O(log n)）
//
//  機能一覧（OfflineDynamicConnectivity）                             計算量
//    OfflineDynamicConnectivity dc(n)   頂点数 n・辺なしで作る           O(1)
//    dc.add_edge(u, v)                  辺を追加する（多重辺・自己ループ可）  O(1)
//    dc.remove_edge(u, v)               辺を 1 本削除する（その時点である辺）  O(1)
//    dc.query_same(u, v)   -> int       「今 u と v は連結か」を登録し、質問番号を返す  O(1)
//    dc.query_size(x)      -> int       「今 x の連結成分の頂点数は」を登録        O(1)
//    dc.query_group_count()-> int       「今の連結成分の個数は」を登録             O(1)
//    dc.query()            -> int       汎用の質問を登録（run(f) の f で自分で答える）  O(1)
//    dc.run()              -> vector<int>  まとめて処理。ans[質問番号] が答え（same は 1 / 0）
//    dc.run(f)             -> vector<int>  汎用の質問の時点ごとに f(質問番号, uf) を呼ぶ
//                                       run 全体で O(n + K log Q log n + Q log n)  K = 辺の操作数
//    dc.n() / dc.query_count()          頂点数 / 登録した質問の個数          O(1)
//
//  機能一覧（RollbackUnionFind … 単体でも使える「巻き戻せる Union-Find」）
//    RollbackUnionFind uf(n)            頂点数 n で作る                      O(n)
//    uf.unite(x, y)        -> bool      併合。新しく併合したら true          O(log n)
//    uf.same(x, y) / uf.find(x) / uf.size(x) / uf.group_count()           O(log n) / O(1)
//    uf.undo()                          直前の unite 1 回分を取り消す        O(1)
//    uf.snapshot()         -> int       今の状態の番号                       O(1)
//    uf.rollback(s)                     snapshot() が s を返した状態まで戻す  O(戻す回数)
//    uf.groups()                        全グループ                           O(n log n)
//    uf.init(n) / uf.n()                作り直す / 頂点数
//
//  ファイル構成
//    [テンプレート本体]  ← 提出時はここだけコピーすればよい
//    odc_test::random_test()  … 愚直解（毎回グラフから連結成分を求め直す）と全機能を突き合わせる
//    odc_test::benchmark()    … 速度計測
//    odc_demo::demo()         … 全機能の使い方デモ（日本語の説明付き）
// =============================================================================

// ======================= ここからテンプレート本体（コピー範囲） =======================

// -----------------------------------------------------------------------------
//  RollbackUnionFind: 巻き戻し（undo）できる Union-Find
//    経路圧縮をしない代わりに、unite を 1 回ずつ取り消せる。find は O(log n)。
// -----------------------------------------------------------------------------
struct RollbackUnionFind {
  public:
    RollbackUnionFind() : RollbackUnionFind(0) {}
    explicit RollbackUnionFind(int n) { init(n); }

    // 頂点数 n で作り直す（履歴も消える）  O(n)
    void init(int n) {
        n_ = n;
        cnt_ = n;
        par_.assign(n, -1);
        hist_.clear();
    }

    // 頂点数
    int n() const { return n_; }

    // x の根（代表元）  O(log n)
    int find(int x) const {
        while (par_[x] >= 0) x = par_[x];
        return x;
    }

    // 併合。新しく併合したら true。false の場合も「1 回の unite」として履歴に積む  O(log n)
    bool unite(int x, int y) {
        x = find(x);
        y = find(y);
        if (x == y) {
            hist_.emplace_back(-1, 0);
            return false;
        }
        if (par_[x] > par_[y]) std::swap(x, y);  // サイズの大きい方 x を根にする
        hist_.emplace_back(y, par_[y]);
        par_[x] += par_[y];
        par_[y] = x;
        --cnt_;
        return true;
    }

    // x と y が同じグループなら true  O(log n)
    bool same(int x, int y) const { return find(x) == find(y); }

    // x のグループの要素数  O(log n)
    int size(int x) const { return -par_[find(x)]; }

    // グループの個数  O(1)
    int group_count() const { return cnt_; }

    // 直前の unite を 1 回分取り消す（false だった unite も 1 回と数える）  O(1)
    void undo() {
        assert(!hist_.empty());
        const int y = hist_.back().first, py = hist_.back().second;
        hist_.pop_back();
        if (y < 0) return;
        const int x = par_[y];
        par_[x] -= py;
        par_[y] = py;
        ++cnt_;
    }

    // 今の状態を表す番号（= 取り消していない unite の回数）。rollback(s) でこの状態に戻れる  O(1)
    int snapshot() const { return (int)hist_.size(); }

    // snapshot() が s を返した状態まで戻す  O(戻す回数)
    void rollback(int s) {
        while ((int)hist_.size() > s) undo();
    }

    // 全グループ（最小要素の昇順・各グループ内も昇順）  O(n log n)
    std::vector<std::vector<int>> groups() const {
        std::vector<int> id(n_, -1);
        std::vector<std::vector<int>> res;
        for (int i = 0; i < n_; ++i) {
            const int r = find(i);
            if (id[r] < 0) {
                id[r] = (int)res.size();
                res.emplace_back();
            }
            res[id[r]].push_back(i);
        }
        return res;
    }

  private:
    int n_ = 0;
    int cnt_ = 0;
    std::vector<int> par_;                   // 根なら -(サイズ)、それ以外は親
    std::vector<std::pair<int, int>> hist_;  // (吸収された根, その元の par_ 値)。false の unite は (-1, 0)
};

// -----------------------------------------------------------------------------
//  OfflineDynamicConnectivity: 辺の追加・削除がある連結性の問題を「オフライン」で解く
//    操作と質問を先に全部登録してから run() でまとめて答えを求める。
//    仕組み: 各辺が存在する期間（質問番号の区間）を時間軸のセグメント木に載せ、
//            DFS しながら RollbackUnionFind で unite / 巻き戻しをする。
//    計算量: O(n + (辺の操作数) log Q log n + Q log n)   Q = 質問の個数
// -----------------------------------------------------------------------------
struct OfflineDynamicConnectivity {
  public:
    OfflineDynamicConnectivity() : OfflineDynamicConnectivity(0) {}
    explicit OfflineDynamicConnectivity(int n) : n_(n) {}

    // 頂点数
    int n() const { return n_; }
    // これまでに登録した質問の個数
    int query_count() const { return (int)qs_.size(); }

    // 辺 (u, v) を追加する（多重辺・自己ループも可）
    void add_edge(int u, int v) { push_event(u, v, true); }
    // 辺 (u, v) を 1 本削除する（その時点で存在する辺であること）
    void remove_edge(int u, int v) { push_event(u, v, false); }

    // 「この時点で u と v は連結か」を登録し、質問番号を返す（答えは 1 / 0）
    int query_same(int u, int v) { return push_query(kSame, u, v); }
    // 「この時点で x の連結成分の頂点数は」を登録し、質問番号を返す
    int query_size(int x) { return push_query(kSize, x, x); }
    // 「この時点での連結成分の個数は」を登録し、質問番号を返す
    int query_group_count() { return push_query(kCount, 0, 0); }
    // 汎用の質問を登録し、質問番号を返す。run(f) で f(質問番号, uf) が呼ばれる
    int query() { return push_query(kCustom, 0, 0); }

    // まとめて処理し、ans[質問番号] を返す（汎用の質問の欄は -1）
    std::vector<int> run() {
        return run([](int, RollbackUnionFind&) {});
    }

    // まとめて処理する。汎用の質問 query() の時点ごとに f(質問番号, uf) を呼ぶ。
    //   uf はその時点のグラフの連結性を表す RollbackUnionFind（same / size / find / group_count 等を使う。
    //   unite した場合は f の中で同じ回数 undo して元に戻すこと）。f は質問番号の小さい順に呼ばれる。
    template <class F>
    std::vector<int> run(F f) {
        const int Q = (int)qs_.size();
        std::vector<int> ans(Q, -1);
        if (Q == 0) return ans;

        // 1) 同じ辺の追加と削除を対応付け、辺が存在する質問番号の区間 [l, r) を求める
        std::vector<int> ord(ev_.size());
        for (int i = 0; i < (int)ord.size(); ++i) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](int a, int b) {
            if (ev_[a].u != ev_[b].u) return ev_[a].u < ev_[b].u;
            if (ev_[a].v != ev_[b].v) return ev_[a].v < ev_[b].v;
            return a < b;
        });
        std::vector<int> iu, iv, il, ir, alive;
        for (int i = 0; i < (int)ord.size();) {
            int j = i;
            const int u = ev_[ord[i]].u, v = ev_[ord[i]].v;
            alive.clear();
            for (; j < (int)ord.size() && ev_[ord[j]].u == u && ev_[ord[j]].v == v; ++j) {
                const Event& e = ev_[ord[j]];
                if (e.add) {
                    alive.push_back(e.t);
                } else {
                    assert(!alive.empty());  // 存在しない辺の削除
                    const int l = alive.back();
                    alive.pop_back();
                    if (l < e.t && u != v) iu.push_back(u), iv.push_back(v), il.push_back(l), ir.push_back(e.t);
                }
            }
            for (int l : alive)
                if (l < Q && u != v) iu.push_back(u), iv.push_back(v), il.push_back(l), ir.push_back(Q);
            i = j;
        }

        // 2) 区間を時間軸のセグメント木のノードに載せる（CSR 形式）
        size_ = 1;
        while (size_ < Q) size_ <<= 1;
        start_.assign(2 * size_ + 1, 0);
        const int m = (int)iu.size();
        for (int k = 0; k < m; ++k)
            for (int l = il[k] + size_, r = ir[k] + size_; l < r; l >>= 1, r >>= 1) {
                if (l & 1) ++start_[l++];
                if (r & 1) ++start_[--r];
            }
        for (int k = 0; k < 2 * size_; ++k) start_[k + 1] += start_[k];
        eu_.resize(start_[2 * size_]);
        ev2_.resize(start_[2 * size_]);
        for (int k = m - 1; k >= 0; --k)
            for (int l = il[k] + size_, r = ir[k] + size_; l < r; l >>= 1, r >>= 1) {
                if (l & 1) {
                    const int p = --start_[l++];
                    eu_[p] = iu[k], ev2_[p] = iv[k];
                }
                if (r & 1) {
                    const int p = --start_[--r];
                    eu_[p] = iu[k], ev2_[p] = iv[k];
                }
            }

        // 3) DFS: ノードに入るとき unite、出るとき巻き戻し。葉で質問に答える
        RollbackUnionFind uf(n_);
        dfs(1, 0, size_, Q, uf, f, ans);
        return ans;
    }

  private:
    enum { kSame = 0, kSize = 1, kCount = 2, kCustom = 3 };
    struct Event {
        int u, v, t;  // t = この操作より前に登録された質問の個数
        bool add;
    };
    struct Query {
        int type, u, v;
    };
    int n_ = 0;
    std::vector<Event> ev_;
    std::vector<Query> qs_;
    int size_ = 1;
    std::vector<int> start_, eu_, ev2_;

    void push_event(int u, int v, bool add) {
        assert(0 <= u && u < n_ && 0 <= v && v < n_);
        if (u > v) std::swap(u, v);
        ev_.push_back(Event{u, v, (int)qs_.size(), add});
    }
    int push_query(int type, int u, int v) {
        qs_.push_back(Query{type, u, v});
        return (int)qs_.size() - 1;
    }

    template <class F>
    void dfs(int k, int lo, int hi, int Q, RollbackUnionFind& uf, F& f, std::vector<int>& ans) {
        if (lo >= Q) return;  // 質問のない部分は見ない
        const int snap = uf.snapshot();
        for (int i = start_[k]; i < start_[k + 1]; ++i) uf.unite(eu_[i], ev2_[i]);
        if (hi - lo == 1) {
            const Query& q = qs_[lo];
            if (q.type == kSame) ans[lo] = uf.same(q.u, q.v) ? 1 : 0;
            else if (q.type == kSize) ans[lo] = uf.size(q.u);
            else if (q.type == kCount) ans[lo] = uf.group_count();
            else f(lo, uf);
        } else {
            const int mid = (lo + hi) >> 1;
            dfs(2 * k, lo, mid, Q, uf, f, ans);
            dfs(2 * k + 1, mid, hi, Q, uf, f, ans);
        }
        uf.rollback(snap);
    }
};
// ======================= ここまでテンプレート本体 =======================


// =============================================================================
//  ランダムテスト: 愚直解（今ある辺から毎回 BFS で連結成分を求め直す）と全機能を突き合わせる
// =============================================================================
namespace odc_test {

// NDEBUG（VS の Release ビルド）でも無効にならないチェック
#define ODC_CHECK(cond)                                                                        \
    do {                                                                                       \
        if (!(cond)) {                                                                         \
            std::cerr << "[OfflineDynamicConnectivity] TEST FAILED  line " << __LINE__ << ": " \
                      << #cond << "\n";                                                        \
            std::exit(1);                                                                      \
        }                                                                                      \
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

// 愚直解: 辺の一覧から BFS で連結成分番号を求める（番号は最小頂点の順に 0, 1, 2, ...）
static std::vector<int> components(int n, const std::vector<std::pair<int, int>>& es) {
    std::vector<std::vector<int>> g(n);
    for (const auto& e : es) g[e.first].push_back(e.second), g[e.second].push_back(e.first);
    std::vector<int> lab(n, -1), q;
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
    return lab;
}
static int count_groups(const std::vector<int>& lab) {
    int k = 0;
    for (int x : lab) k = std::max(k, x + 1);
    return k;
}
static int comp_size(const std::vector<int>& lab, int x) {
    int c = 0;
    for (int y : lab) c += (y == lab[x]);
    return c;
}
static std::vector<std::vector<int>> label_groups(const std::vector<int>& lab) {
    std::vector<std::vector<int>> g(count_groups(lab));
    for (int v = 0; v < (int)lab.size(); ++v) g[lab[v]].push_back(v);
    return g;
}

// uf の状態が成分番号 lab と一致するか（find / same / size / group_count / groups）
static void check_uf(const RollbackUnionFind& uf, const std::vector<int>& lab) {
    const int n = (int)lab.size();
    ODC_CHECK(uf.n() == n);
    ODC_CHECK(uf.group_count() == count_groups(lab));
    std::vector<int> root_of(n, -1), sz(n, 0);
    for (int v = 0; v < n; ++v) ++sz[lab[v]];
    for (int v = 0; v < n; ++v) {
        const int r = uf.find(v);
        ODC_CHECK(0 <= r && r < n && uf.find(r) == r && lab[r] == lab[v]);
        if (root_of[lab[v]] < 0) root_of[lab[v]] = r;
        ODC_CHECK(root_of[lab[v]] == r);
        ODC_CHECK(uf.size(v) == sz[lab[v]]);
    }
    ODC_CHECK(uf.groups() == label_groups(lab));
}

// ---------------- RollbackUnionFind 単体のテスト ----------------
static long long rollback_case(std::mt19937_64& rng, int n, int q) {
    RollbackUnionFind uf(n);
    std::vector<std::vector<int>> st(1, std::vector<int>(n));  // st[s] = unite を s 回したときの成分番号
    for (int i = 0; i < n; ++i) st[0][i] = i;
    auto normalize = [](const std::vector<int>& raw) {  // 代表の番号 → 最小頂点順の番号
        std::vector<int> id(raw.size(), -1), lab(raw.size());
        int k = 0;
        for (size_t v = 0; v < raw.size(); ++v) {
            if (id[raw[v]] < 0) id[raw[v]] = k++;
            lab[v] = id[raw[v]];
        }
        return lab;
    };
    for (int t = 0; t < q; ++t) {
        const int type = rnd(rng, 0, 99);
        if (n == 0) break;
        const int x = rnd(rng, 0, n - 1), y = rnd(rng, 0, n - 1);
        if (type < 40) {  // unite
            std::vector<int> cur = st.back();
            const bool expect = cur[x] != cur[y];
            const int from = cur[y], to = cur[x];
            for (int& l : cur)
                if (l == from) l = to;
            st.push_back(cur);
            const int rx = uf.find(x), ry = uf.find(y), sx = uf.size(x), sy = uf.size(y);
            ODC_CHECK(uf.unite(x, y) == expect);
            ODC_CHECK(uf.snapshot() == (int)st.size() - 1);
            if (expect) ODC_CHECK(uf.find(y) == (sx >= sy ? rx : ry));  // union by size（同じなら x 側が根）
        } else if (type < 55) {  // undo
            if (st.size() > 1) {
                st.pop_back();
                uf.undo();
            }
        } else if (type < 62) {  // rollback
            const int s = rnd(rng, 0, (int)st.size() - 1);
            st.resize(s + 1);
            uf.rollback(s);
            ODC_CHECK(uf.snapshot() == s);
        } else if (type < 72) {  // same
            ODC_CHECK(uf.same(x, y) == (st.back()[x] == st.back()[y]));
        } else if (type < 80) {  // size
            ODC_CHECK(uf.size(x) == comp_size(normalize(st.back()), x));
        } else if (type < 86) {  // find
            ODC_CHECK((uf.find(x) == uf.find(y)) == (st.back()[x] == st.back()[y]));
        } else if (type < 90) {  // group_count
            ODC_CHECK(uf.group_count() == count_groups(normalize(st.back())));
        } else if (type < 99) {  // 全体チェック
            check_uf(uf, normalize(st.back()));
        } else {  // init
            uf.init(n);
            st.assign(1, std::vector<int>(n));
            for (int i = 0; i < n; ++i) st[0][i] = i;
            ODC_CHECK(uf.snapshot() == 0);
        }
    }
    check_uf(uf, normalize(st.back()));
    return q;
}

// ---------------- OfflineDynamicConnectivity のテスト ----------------
//   操作をランダムに登録しながら、その時点の正解を愚直に求めておき、run() の結果と比べる
static long long offline_case(std::mt19937_64& rng, int n, int q) {
    OfflineDynamicConnectivity dc(n);
    std::vector<std::pair<int, int>> alive;       // 今ある辺（多重辺は複数入る）
    std::vector<int> expect;                      // 型つき質問の正解（汎用は -1）
    std::vector<std::vector<int>> custom_lab;     // 汎用質問の時点の成分番号（型つき質問は空）
    std::vector<int> lab = components(n, alive);  // 今の成分番号（辺が変わったら作り直す）
    bool dirty = false;
    const int big = std::max(1, n);
    for (int t = 0; t < q; ++t) {
        const int type = rnd(rng, 0, 99);
        if (type < 30) {  // 辺の追加（たまに自己ループ・既存の辺と同じ辺）
            int u = rnd(rng, 0, big - 1), v = rnd(rng, 0, big - 1);
            if (!alive.empty() && rnd(rng, 0, 5) == 0) {  // 多重辺
                const int k = rnd(rng, 0, (int)alive.size() - 1);
                u = alive[k].second, v = alive[k].first;
            }
            if (n == 0) continue;
            dc.add_edge(u, v);
            alive.emplace_back(u, v);
            dirty = true;
        } else if (type < 55) {  // 辺の削除（端点の順番もランダム）
            if (alive.empty()) continue;
            const int k = rnd(rng, 0, (int)alive.size() - 1);
            std::swap(alive[k], alive.back());
            const std::pair<int, int> e = alive.back();
            alive.pop_back();
            if (rng() & 1) dc.remove_edge(e.first, e.second);
            else dc.remove_edge(e.second, e.first);
            dirty = true;
        } else {  // 質問
            if (n == 0) {
                ODC_CHECK(dc.query_group_count() == (int)expect.size());
                expect.push_back(0);
                custom_lab.emplace_back();
                continue;
            }
            if (dirty) lab = components(n, alive), dirty = false;
            const int x = rnd(rng, 0, n - 1), y = rnd(rng, 0, n - 1);
            int id;
            if (type < 70) id = dc.query_same(x, y), expect.push_back(lab[x] == lab[y] ? 1 : 0);
            else if (type < 80) id = dc.query_size(x), expect.push_back(comp_size(lab, x));
            else if (type < 88) id = dc.query_group_count(), expect.push_back(count_groups(lab));
            else id = dc.query(), expect.push_back(-1);
            ODC_CHECK(id == (int)expect.size() - 1);
            custom_lab.push_back(expect.back() == -1 ? lab : std::vector<int>());
        }
    }
    ODC_CHECK(dc.n() == n);
    ODC_CHECK(dc.query_count() == (int)expect.size());

    // run(f): 汎用質問の時点ごとに uf を全体チェック。f の中で unite → undo しても壊れないことも確認
    int last = -1, called = 0;
    const std::vector<int> ans = dc.run([&](int id, RollbackUnionFind& uf) {
        ODC_CHECK(id > last);  // 質問番号の小さい順
        last = id;
        ++called;
        ODC_CHECK(!custom_lab[id].empty());
        check_uf(uf, custom_lab[id]);
        if (n >= 2) {
            const int a = rnd(rng, 0, n - 1), b = rnd(rng, 0, n - 1);
            uf.unite(a, b);
            uf.unite(b, (a + 1) % n);
            uf.undo();
            uf.undo();
        }
    });
    ODC_CHECK(ans == expect);
    int customs = 0;
    for (int x : expect) customs += (x == -1);
    ODC_CHECK(called == customs);
    ODC_CHECK(dc.run() == expect);  // 何回呼んでも同じ
    return q;
}

// 形の決まったケース
static void structured_cases(std::mt19937_64& rng) {
    {  // 質問なし・頂点なし
        OfflineDynamicConnectivity dc;
        ODC_CHECK(dc.run().empty() && dc.n() == 0);
        OfflineDynamicConnectivity dc2(3);
        dc2.add_edge(0, 1);
        ODC_CHECK(dc2.run().empty());
    }
    {  // 最後の質問より後の操作は影響しない / 追加して質問なしで削除
        OfflineDynamicConnectivity dc(4);
        dc.add_edge(0, 1);
        dc.remove_edge(1, 0);
        const int a = dc.query_same(0, 1);
        dc.add_edge(0, 1);
        const int b = dc.query_same(0, 1);
        dc.add_edge(2, 3);
        const std::vector<int> ans = dc.run();
        ODC_CHECK(ans[a] == 0 && ans[b] == 1 && ans.size() == 2);
    }
    {  // 多重辺と自己ループ
        OfflineDynamicConnectivity dc(3);
        dc.add_edge(0, 1), dc.add_edge(1, 0), dc.add_edge(0, 1), dc.add_edge(2, 2);
        const int a = dc.query_size(0);
        dc.remove_edge(0, 1), dc.remove_edge(0, 1);
        const int b = dc.query_same(1, 0);
        dc.remove_edge(1, 0);
        const int c = dc.query_same(0, 1);
        const int d = dc.query_group_count();
        dc.remove_edge(2, 2);
        const std::vector<int> ans = dc.run();
        ODC_CHECK(ans[a] == 2 && ans[b] == 1 && ans[c] == 0 && ans[d] == 3);
    }
    {  // 長い道: 全部つないで、真ん中を切ったりつないだり
        const int n = 1000;
        OfflineDynamicConnectivity dc(n);
        for (int i = 0; i + 1 < n; ++i) dc.add_edge(i, i + 1);
        std::vector<int> exp;
        dc.query_size(0), exp.push_back(n);
        for (int k = 0; k < 200; ++k) {
            const int m = rnd(rng, 0, n - 2);
            dc.remove_edge(m, m + 1);
            dc.query_size(0), exp.push_back(m + 1);
            dc.query_group_count(), exp.push_back(2);
            dc.add_edge(m + 1, m);
            dc.query_same(0, n - 1), exp.push_back(1);
        }
        ODC_CHECK(dc.run() == exp);
    }
    {  // RollbackUnionFind: 失敗した unite も undo 1 回分
        RollbackUnionFind uf(3);
        ODC_CHECK(uf.unite(0, 1) && !uf.unite(1, 0) && uf.snapshot() == 2);
        uf.undo();
        ODC_CHECK(uf.same(0, 1));
        uf.undo();
        ODC_CHECK(!uf.same(0, 1) && uf.group_count() == 3);
    }
}

static void random_test() {
    std::mt19937_64 rng(20260923);
    long long ops = 0;
    int cases = 0;
    for (int it = 0; it < 400 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 0, 10), q = rnd(rng, 0, 80);
        ops += rollback_case(rng, n, q);
    }
    for (int it = 0; it < 30 * TEST_SCALE; ++it, ++cases) {
        const int n = rnd(rng, 11, 80);
        ops += rollback_case(rng, n, 4 * n);
    }
    for (int it = 0; it < 500 * TEST_SCALE; ++it, ++cases) {  // 小さい n（衝突が多い）
        const int n = rnd(rng, 0, 9), q = rnd(rng, 0, 80);
        ops += offline_case(rng, n, q);
    }
    for (int it = 0; it < 40 * TEST_SCALE; ++it, ++cases) {  // 中くらい
        const int n = rnd(rng, 10, 120);
        ops += offline_case(rng, n, 6 * n);
    }
    for (int it = 0; it < TEST_SCALE; ++it, ++cases) {  // 大きめ
        const int n = rnd(rng, 500, 2000);
        ops += offline_case(rng, n, 3 * n);
    }
    structured_cases(rng);
    std::cout << "[OfflineDynamicConnectivity] random test: OK  (" << cases << " cases, " << ops
              << " ops + structured cases)\n";
}

// 速度計測: n 頂点に、辺の追加 40%・削除 20%・質問 same 40% をランダムに合計 q 個
static void benchmark() {
    if (TEST_SCALE == 1)
        std::cout << "[OfflineDynamicConnectivity] benchmark: (デバッグビルドなので遅い。速度は Release ビルドで確認する)\n";
    std::mt19937_64 rng(1);
    const int n = 100000, q = TEST_SCALE == 1 ? 50000 : 300000;
    OfflineDynamicConnectivity dc(n);
    std::vector<std::pair<int, int>> alive;
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < q; ++i) {
        const int type = rnd(rng, 0, 9);
        if (type < 4 || alive.empty()) {
            const int u = rnd(rng, 0, n - 1), v = rnd(rng, 0, n - 1);
            dc.add_edge(u, v);
            alive.emplace_back(u, v);
        } else if (type < 6) {
            const int k = rnd(rng, 0, (int)alive.size() - 1);
            std::swap(alive[k], alive.back());
            dc.remove_edge(alive.back().first, alive.back().second);
            alive.pop_back();
        } else {
            dc.query_same(rnd(rng, 0, n - 1), rnd(rng, 0, n - 1));
        }
    }
    const std::vector<int> ans = dc.run();
    const auto t1 = std::chrono::steady_clock::now();
    long long chk = 0;
    for (int x : ans) chk += x;
    std::cout << "[OfflineDynamicConnectivity] benchmark: n = " << n << ", " << q
              << " ops (add 40% / remove 20% / same 40%)  "
              << (long long)std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms  (checksum " << chk
              << ")\n";
}

}  // namespace odc_test


// =============================================================================
//  使い方デモ: 全 API を 1 つずつ「返り値・内容・計算量・注意点」つきで実行して見せる。
//    出力の見方:   実際に書くコード  ->  返り値    // 何を返しているか（日本語）
// =============================================================================
namespace odc_demo {

static std::string to_s(bool b) { return b ? "true" : "false"; }
static std::string to_s(int v) { return std::to_string(v); }
static std::string to_s(const std::string& s) { return s; }
static std::string to_s(const std::vector<int>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + std::to_string(v[i]);
    return s + "}";
}
static std::string to_s(const std::vector<std::vector<int>>& g) {
    std::string s = "{";
    for (size_t i = 0; i < g.size(); ++i) s += (i ? ", " : "") + to_s(g[i]);
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
// RollbackUnionFind の今のグループ分けを 1 行で表示
static void state(const RollbackUnionFind& uf) {
    std::string s;
    for (const auto& g : uf.groups()) {
        s += "{";
        for (size_t i = 0; i < g.size(); ++i) s += (i ? "," : "") + std::to_string(g[i]);
        s += "} ";
    }
    std::cout << "    [状態] " << s << " (グループ数 " << uf.group_count() << ", snapshot " << uf.snapshot() << ")\n";
}
template <class V>
static void show(const char* expr, const V& v, const std::string& comment) {  // コード -> 返り値 // 説明
    std::string line = "  " + pad(expr, 30) + " -> " + to_s(v);
    if (!comment.empty()) line = pad(line, 44) + "  // " + comment;
    std::cout << line << "\n";
}
// SHOW(式, "説明") : 式そのものと返り値と説明を表示
#define SHOW(expr, comment) odc_demo::show(#expr, (expr), comment)
// RUN(文) : 文を表示してから実行（宣言にも使える）
#define RUN(...)                   \
    odc_demo::stmt(#__VA_ARGS__); \
    __VA_ARGS__

static void cheat_sheet() {
    struct Row {
        const char* use;
        const char* ret;
        const char* what;
    };
    const Row rows[] = {
        {"OfflineDynamicConnectivity dc(n);", "", "頂点 0 ~ n-1、辺なしで作る"},
        {"dc.add_edge(u, v)", "", "辺を追加する（多重辺・自己ループ可）"},
        {"dc.remove_edge(u, v)", "", "辺を 1 本削除する（その時点である辺）"},
        {"dc.query_same(u, v)", "int", "「今 u と v は連結か」を登録 → 質問番号"},
        {"dc.query_size(x)", "int", "「今 x の成分の頂点数は」を登録 → 質問番号"},
        {"dc.query_group_count()", "int", "「今の連結成分の個数は」を登録 → 質問番号"},
        {"dc.query()", "int", "汎用の質問を登録 → 質問番号（run(f) の f で答える）"},
        {"dc.run()", "vector<int>", "まとめて処理。ans[質問番号] が答え（same は 1 / 0、汎用は -1）"},
        {"dc.run(f)", "vector<int>", "汎用の質問の時点ごとに f(質問番号, uf) を呼ぶ"},
        {"dc.n() / dc.query_count()", "int", "頂点数 / 登録した質問の個数"},
        {"RollbackUnionFind uf(n);", "", "巻き戻せる Union-Find（単体でも使える）"},
        {"uf.unite(x, y)", "bool", "併合。新しく併合したら true  O(log n)"},
        {"uf.same(x, y)", "bool", "同じグループなら true  O(log n)"},
        {"uf.find(x)", "int", "根（代表元）  O(log n)"},
        {"uf.size(x)", "int", "x のグループの要素数  O(log n)"},
        {"uf.group_count()", "int", "グループの個数  O(1)"},
        {"uf.undo()", "", "直前の unite 1 回分を取り消す  O(1)"},
        {"uf.snapshot()", "int", "今の状態の番号（= 取り消していない unite の回数）  O(1)"},
        {"uf.rollback(s)", "", "snapshot() が s だった状態まで戻す  O(戻す回数)"},
        {"uf.groups()", "vector<vector<int>>", "全グループ  O(n log n)"},
        {"uf.init(n) / uf.n()", "", "作り直す / 頂点数"},
    };
    std::cout << "\n■ API 一覧（頂点番号は 0 ~ n-1。詳しくは下の【1】~【20】）\n";
    std::cout << "  ★ 流れ: ① 操作と質問を「起きる順番」に登録 → ② run() を 1 回呼ぶ → ③ ans[質問番号] を見る\n";
    for (const Row& r : rows) {
        const std::string ret = (r.ret[0] ? std::string("-> ") + r.ret : std::string());
        std::cout << "  " << pad(r.use, 36) << pad(ret, 24) << r.what << "\n";
    }
}

static void demo() {
    banner("辺の削除もできる Union-Find（オフライン版）使い方デモ    見方:  コード  ->  返り値    // 意味");
    cheat_sheet();

    // ================================================================== OfflineDynamicConnectivity
    api("【1】 OfflineDynamicConnectivity dc(n);", "（コンストラクタ）",
        "頂点 0 ~ n-1、辺なしで作る。ここに操作と質問を「起きる順番」に登録していく", "O(1)");
    RUN(OfflineDynamicConnectivity dc(5));
    note("質問は「登録した時点のグラフ」について答えが出る。その後の追加・削除は影響しない");
    note("答えは run() を呼ぶまで分からない（オフライン）。その場で答えが必要ならオンライン版を使う");

    api("【2】 dc.add_edge(u, v)", "なし", "辺 (u, v) を追加する（多重辺・自己ループも可。u と v の順番はどちらでもよい）",
        "O(1)（登録するだけ）");
    RUN(dc.add_edge(0, 1));
    RUN(dc.add_edge(1, 2));
    note("[この時点のグラフ] 0-1-2   3   4");

    api("【3】 dc.query_same(u, v)", "int … 質問番号（登録した順に 0, 1, 2, ...）",
        "「この時点で u と v は連結か」を登録する。答え（1 / 0）は run() の後に ans[質問番号] で見る", "O(1)");
    SHOW(dc.query_same(0, 2), "質問 0: 今 0 と 2 は連結？");
    SHOW(dc.query_same(0, 3), "質問 1: 今 0 と 3 は連結？");

    api("【4】 dc.query_size(x)", "int … 質問番号", "「この時点で x の連結成分の頂点数は」を登録する", "O(1)");
    SHOW(dc.query_size(1), "質問 2: 今 1 の成分の頂点数は？");

    api("【5】 dc.query_group_count()", "int … 質問番号", "「この時点での連結成分の個数は」を登録する", "O(1)");
    SHOW(dc.query_group_count(), "質問 3: 今の連結成分の個数は？");

    api("【6】 dc.remove_edge(u, v)", "なし",
        "辺 (u, v) を 1 本削除する。その時点で存在する辺であること（u と v の順番はどちらでもよい）", "O(1)（登録するだけ）");
    note("※ 存在しない辺を削除すると run() の中で assert で止まる");
    RUN(dc.remove_edge(2, 1));
    note("[この時点のグラフ] 0-1   2   3   4");
    SHOW(dc.query_same(0, 2), "質問 4: 今 0 と 2 は連結？（1-2 を消した後）");
    RUN(dc.add_edge(3, 4));
    RUN(dc.add_edge(2, 3));
    note("[この時点のグラフ] 0-1   2-3-4");
    SHOW(dc.query_size(4), "質問 5: 今 4 の成分の頂点数は？");
    SHOW(dc.query_group_count(), "質問 6: 今の連結成分の個数は？");
    note("多重辺: 同じ辺が 2 本あるとき、1 本消してもまだつながっている");
    RUN(dc.add_edge(1, 0));
    RUN(dc.remove_edge(0, 1));
    SHOW(dc.query_same(0, 1), "質問 7: 今 0 と 1 は連結？（0-1 が 1 本残っている）");

    api("【7】 dc.query()", "int … 質問番号",
        "汎用の質問を登録する。run(f) のときに、この時点のグラフを表す uf が f(質問番号, uf) に渡される", "O(1)");
    SHOW(dc.query(), "質問 8: 汎用（【10】の run(f) で自分で答える）");

    api("【8】 dc.n() / dc.query_count()", "int", "頂点数 / これまでに登録した質問の個数（質問番号は 0 ~ query_count()-1）",
        "O(1)");
    SHOW(dc.n(), "頂点数");
    SHOW(dc.query_count(), "質問 0 ~ 8 の 9 個");

    api("【9】 dc.run()", "vector<int> … ans[質問番号] = 答え",
        "登録した全部をまとめて処理する。same は 1 / 0、size・group_count は個数。汎用の質問の欄は -1",
        "O(n + K log Q log n + Q log n)  K = 辺の操作数、Q = 質問の個数");
    note("※ run() は何度呼んでもよい（登録内容は消えない）。run() の後にさらに登録して呼び直してもよい");
    RUN(std::vector<int> ans = dc.run());
    SHOW(ans, "質問 0 ~ 8 の答え");
    SHOW(ans[0], "質問 0: 0-1-2 なので 0 と 2 は連結 → 1");
    SHOW(ans[1], "質問 1: 3 は孤立 → 0");
    SHOW(ans[2], "質問 2: {0,1,2} → 3");
    SHOW(ans[3], "質問 3: {0,1,2} {3} {4} → 3");
    SHOW(ans[4], "質問 4: 1-2 を消したので → 0");
    SHOW(ans[5], "質問 5: {2,3,4} → 3");
    SHOW(ans[6], "質問 6: {0,1} {2,3,4} → 2");
    SHOW(ans[7], "質問 7: 多重辺の 1 本が残っている → 1");
    SHOW(ans[8], "質問 8: 汎用の質問は -1");

    api("【10】 dc.run(f)", "vector<int>（run() と同じ）",
        "汎用の質問 query() の時点ごとに f(質問番号, uf) を呼ぶ。uf はその時点のグラフを表す RollbackUnionFind",
        "run() と同じ + f の計算量");
    note("f は [&](int id, RollbackUnionFind& uf) { ... } の形。質問番号の小さい順に呼ばれる");
    note("uf では same / find / size / group_count / groups が使える（unite したら f の中で同じ回数 undo して戻す）");
    RUN(std::vector<std::string> info(dc.query_count()));
    RUN(auto f = [&](int id, RollbackUnionFind& uf) { info[id] = "groups = " + std::to_string(uf.group_count()) + ", size(2) = " + std::to_string(uf.size(2)); });
    RUN(dc.run(f));
    SHOW(info[8], "質問 8 の時点の情報を自分で作った");

    // ================================================================== RollbackUnionFind
    {
    api("【11】 RollbackUnionFind uf(n);", "（コンストラクタ）",
        "巻き戻せる Union-Find。上の run(f) の uf の型で、単体でも使える", "O(n)");
    note("普通の UnionFind と違って経路圧縮をしないので find などは O(log n)。その代わり unite を取り消せる");
    RUN(RollbackUnionFind uf(5));
    state(uf);

    api("【12】 uf.unite(x, y)", "bool … 新しく併合したら true / 既に同じグループなら false",
        "x と y のグループを併合する。false の場合も「1 回の unite」として履歴に残る（undo で 1 回分）", "O(log n)");
    SHOW(uf.unite(0, 1), "");
    SHOW(uf.unite(2, 3), "");
    state(uf);
    SHOW(uf.unite(1, 0), "既に同じ → false（でも履歴には 1 回分残る）");
    state(uf);

    api("【13】 uf.same(x, y)", "bool", "x と y が同じグループなら true", "O(log n)");
    SHOW(uf.same(0, 1), "");
    SHOW(uf.same(1, 2), "");

    api("【14】 uf.find(x)", "int … x の根（代表元）", "同じグループなら同じ値。どの頂点が根になるかに意味はない", "O(log n)");
    SHOW(uf.find(1), "");
    SHOW(uf.find(3), "");

    api("【15】 uf.size(x)", "int", "x のグループの要素数", "O(log n)");
    SHOW(uf.size(0), "{0,1}");
    SHOW(uf.size(4), "{4}");

    api("【16】 uf.group_count()", "int", "グループの個数", "O(1)");
    SHOW(uf.group_count(), "{0,1} {2,3} {4}");

    api("【17】 uf.undo()", "なし", "直前の unite を 1 回分取り消す（false だった unite も 1 回と数える）", "O(1)");
    RUN(uf.undo());
    note("↑ 直前の unite(1, 0) は false だったので、見た目は何も変わらない");
    state(uf);
    RUN(uf.undo());
    note("↑ unite(2, 3) が取り消された");
    state(uf);

    api("【18】 uf.snapshot() / uf.rollback(s)", "snapshot: int（今の状態の番号） / rollback: なし",
        "snapshot() で今の番号を覚えておき、rollback(番号) でその状態まで一気に戻す", "snapshot O(1) / rollback O(戻す回数)");
    RUN(int s = uf.snapshot());
    SHOW(s, "今の状態の番号（= 取り消していない unite の回数）");
    RUN(uf.unite(1, 2));
    RUN(uf.unite(3, 4));
    RUN(uf.unite(0, 4));
    state(uf);
    RUN(uf.rollback(s));
    note("↑ 3 回の unite がまとめて取り消された");
    state(uf);

    api("【19】 uf.groups()", "vector<vector<int>>", "全グループ（最小要素の昇順・各グループ内も昇順）", "O(n log n)");
    SHOW(uf.groups(), "");

    api("【20】 uf.init(n) / uf.n()", "init: なし / n: int", "頂点数 n で作り直す（履歴も消える） / 頂点数", "O(n) / O(1)");
    RUN(uf.init(3));
    SHOW(uf.n(), "");
    state(uf);
    }

    // ================================================================== 典型パターン
    head("典型1: 「1 u v: 辺を追加」「2 u v: 辺を削除」「3 u v: 連結なら Yes」のクエリ問題（1-indexed）");
    note("提出用のコード:");
    code("int N, Q;");
    code("std::cin >> N >> Q;");
    code("OfflineDynamicConnectivity dc(N);");
    code("for (int i = 0; i < Q; ++i) {");
    code("    int t, u, v;");
    code("    std::cin >> t >> u >> v;");
    code("    --u, --v;");
    code("    if (t == 1) dc.add_edge(u, v);");
    code("    else if (t == 2) dc.remove_edge(u, v);");
    code("    else dc.query_same(u, v);");
    code("}");
    code("for (int a : dc.run()) std::cout << (a ? \"Yes\" : \"No\") << \"\\n\";");
    note("↓ 同じ処理を solve(入力の文字列) にして実行した結果（答えを空白区切りで返す）");
    auto solve = [](const std::string& input) {
        std::istringstream in(input);
        int N, Q;
        in >> N >> Q;
        OfflineDynamicConnectivity d(N);
        for (int i = 0; i < Q; ++i) {
            int t, u, v;
            in >> t >> u >> v;
            --u, --v;
            if (t == 1) d.add_edge(u, v);
            else if (t == 2) d.remove_edge(u, v);
            else d.query_same(u, v);
        }
        std::string out;
        for (int a : d.run()) out += (out.empty() ? "" : " ") + std::string(a ? "Yes" : "No");
        return out;
    };
    SHOW(solve("3 6\n1 1 2\n1 2 3\n3 1 3\n2 1 2\n3 1 3\n3 2 3\n"), "1-2-3 → Yes、1-2 を消す → No、2-3 は Yes");

    head("典型2: 汎用の質問で「その時点の最大の連結成分の大きさ」を求める");
    {
        RUN(OfflineDynamicConnectivity g(6));
        RUN(g.add_edge(0, 1));
        RUN(g.add_edge(1, 2));
        RUN(g.query());
        RUN(g.add_edge(3, 4));
        RUN(g.add_edge(4, 5));
        RUN(g.add_edge(2, 3));
        RUN(g.query());
        RUN(g.remove_edge(2, 3));
        RUN(g.remove_edge(0, 1));
        RUN(g.query());
        RUN(std::vector<int> best(g.query_count()));
        RUN(g.run([&](int id, RollbackUnionFind& uf) { for (int v = 0; v < uf.n(); ++v) best[id] = std::max(best[id], uf.size(v)); }));
        SHOW(best, "{0,1,2} → 3、全部つながって 6、{3,4,5} が最大で 3");
        note("※ ここでは各質問で全頂点を見ている（O(n)）。質問が多いときは別の工夫がいる");
    }

    head("典型3: RollbackUnionFind 単体 … 試しに併合して調べ、元に戻す");
    {
        RUN(RollbackUnionFind r(4));
        RUN(r.unite(0, 1));
        RUN(int s = r.snapshot());
        RUN(r.unite(1, 2));
        RUN(r.unite(2, 3));
        SHOW(r.group_count() == 1, "この 2 本を足すと全体がつながる");
        RUN(r.rollback(s));
        SHOW(r.group_count(), "元に戻った（{0,1} {2} {3}）");
    }

    head("メモ: 削除しかない問題（最初に辺があり、消えていくだけ）");
    note("クエリを逆順に見ると「追加だけ」になるので、普通の UnionFind で解ける（こちらの方が簡単で速い）");
}

}

void solv(){
    odc_test::random_test();  // ランダムテスト（失敗したらその行を表示して終了）
    odc_test::benchmark();    // 速度計測
    odc_demo::demo();         // 使い方デモ
    return;
}

int main(){
    BUFF;
    int T=1;
    //cin>>T;
    while(T--)solv();
    return 0;
}