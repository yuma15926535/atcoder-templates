#ifndef GG_TRIE
#define GG_TRIE
// ---- Trie（文字列の多重集合。文字列第二世代） ----
// Trie TR(A)：A は文字の種類（既定は英小文字）。子は 1 次元の配列（ノード数 × σ 個の int。10^6 文字・σ = 26 で約 100MB）
// 同じ文字列を何個でも入れられる（数を持つ）。消してもノードは残る（数が 0 になるだけ）。ノード 0 が根（空文字列）
struct Trie {
    explicit Trie(const Alphabet& A = Alphabet::lower()) : A_(A), sg_(A.size()) { new_node(-1, -1); }

    // ---- 追加・削除 ----
    template <class Seq>
    int insert(const Seq& s0, long long cnt = 1) {  // s を cnt 個入れる（0 個でもノードは作る）。s のノードを返す  O(|s|)
        if (cnt < 0) str_internal::fail("Trie::insert", "cnt ≥ 0 で使う（減らすときは erase）");
        const auto& s = str_internal::seq(s0);
        int v = 0;
        nd_[0].pass += cnt;
        for (std::size_t i = 0; i < s.size(); ++i) {
            const int c = A_.id(s[i]);
            int u = next_[(std::size_t)v * sg_ + c];
            if (u < 0) u = new_node(v, c), next_[(std::size_t)v * sg_ + c] = u;
            v = u;
            nd_[v].pass += cnt;
        }
        nd_[v].end += cnt;
        return v;
    }
    template <class Seq>
    bool erase(const Seq& s, long long cnt = 1) {  // s を cnt 個消す（入っている数が cnt より少なければ何もせず false）  O(|s|)
        if (cnt < 0) str_internal::fail("Trie::erase", "cnt ≥ 0 で使う");
        const int v = find(s);
        if (v < 0) return cnt == 0;
        if (nd_[v].end < cnt) return false;
        nd_[v].end -= cnt;
        for (int u = v; u >= 0; u = nd_[u].par) nd_[u].pass -= cnt;
        return true;
    }
    template <class Seq>
    long long erase_prefix(const Seq& p) {  // p で始まる文字列を全部消す。消した数を返す  O(|p| + 消えたノード数 × σ)
        const int v = find(p);
        if (v < 0 || nd_[v].pass == 0) return 0;
        const long long c = nd_[v].pass;
        for (int u = nd_[v].par; u >= 0; u = nd_[u].par) nd_[u].pass -= c;
        std::vector<int> st = {v};  // 数が残っている部分木だけを 0 にする
        while (!st.empty()) {
            const int u = st.back();
            st.pop_back();
            nd_[u].pass = nd_[u].end = 0;
            for (int k = 0; k < sg_; ++k) {
                const int w = next_[(std::size_t)u * sg_ + k];
                if (w >= 0 && nd_[w].pass > 0) st.push_back(w);
            }
        }
        return c;
    }
    void clear() { *this = Trie(A_); }

    // ---- 数える ----
    template <class Seq>
    long long count(const Seq& s) const {  // s の数  O(|s|)
        const int v = find(s);
        return v < 0 ? 0 : nd_[v].end;
    }
    template <class Seq>
    bool contains(const Seq& s) const { return count(s) > 0; }
    template <class Seq>
    long long count_prefix(const Seq& p) const {  // p で始まる文字列の数（p 自身も）  O(|p|)
        const int v = find(p);
        return v < 0 ? 0 : nd_[v].pass;
    }
    template <class Seq>
    long long count_prefixes_of(const Seq& s) const {  // s の接頭辞である文字列の数（空文字列・s 自身も）  O(|s|)
        long long r = 0;
        for (int v : path(s)) r += nd_[v].end;
        return r;
    }
    long long size() const { return nd_[0].pass; }  // 入っている文字列の数（同じものも数える）
    bool empty() const { return nd_[0].pass == 0; }

    // ---- 辞書順 ----
    std::vector<int> kth_ids(long long k) const {  // 辞書順 k 番目（0 から。同じものも数える）を番号の列で  O(|答え| σ)
        if (k < 0 || k >= size()) str_internal::fail("Trie::kth", "k = " + std::to_string(k) + " は 0〜size()-1 の外");
        std::vector<int> r;
        int v = 0;
        while (true) {
            if (k < nd_[v].end) return r;
            k -= nd_[v].end;
            for (int c = 0; c < sg_; ++c) {
                const int u = next_[(std::size_t)v * sg_ + c];
                if (u < 0) continue;
                if (k < nd_[u].pass) {
                    v = u, r.push_back(c);
                    break;
                }
                k -= nd_[u].pass;
            }
        }
    }
    std::string kth(long long k) const { return from_ids(kth_ids(k), A_); }  // 辞書順 k 番目の文字列
    template <class Seq>
    long long rank(const Seq& s0) const {  // s より辞書順で小さい文字列の数  O(|s| σ)
        const auto& s = str_internal::seq(s0);
        long long r = 0;
        int v = 0;
        for (std::size_t i = 0; i < s.size(); ++i) {
            r += nd_[v].end;  // s の真の接頭辞は s より小さい
            const int c = A_.id(s[i]);
            for (int k = 0; k < c; ++k) {
                const int u = next_[(std::size_t)v * sg_ + k];
                if (u >= 0) r += nd_[u].pass;
            }
            v = next_[(std::size_t)v * sg_ + c];
            if (v < 0) break;
        }
        return r;
    }
    std::vector<std::vector<int>> words_ids() const {  // 全部を辞書順に（番号の列。同じものは数だけ）
        std::vector<std::vector<int>> r;
        std::vector<int> cur;
        for (int v : preorder()) {
            if (nd_[v].pass == 0) continue;
            cur.resize(nd_[v].dep);
            if (v != 0) cur[nd_[v].dep - 1] = nd_[v].lab;
            for (long long t = 0; t < nd_[v].end; ++t) r.push_back(cur);
        }
        return r;
    }
    std::vector<std::string> words() const {  // 全部を辞書順に（文字列。同じものは数だけ）
        std::vector<std::string> r;
        for (const auto& w : words_ids()) r.push_back(from_ids(w, A_));
        return r;
    }

    // ---- 共通接頭辞（LCP） ----
    template <class Seq>
    int max_lcp(const Seq& s0, bool exclude_self = false) const {  // 入っている文字列との LCP の最大（exclude_self で s を 1 つ除いて。ABC287 E）  O(|s|)
        const auto& s = str_internal::seq(s0);
        if (exclude_self && count(s) == 0) str_internal::fail("Trie::max_lcp", "exclude_self = true は s が入っているときだけ");
        const long long sub = exclude_self ? 1 : 0;  // s の道の上のノードは数を 1 減らして見る
        int v = 0, d = 0;
        for (std::size_t i = 0; i < s.size(); ++i) {
            const int u = next_[(std::size_t)v * sg_ + A_.id(s[i])];
            if (u < 0 || nd_[u].pass - sub <= 0) break;
            v = u, ++d;
        }
        return d;
    }
    template <class Seq>
    long long sum_lcp_with(const Seq& s) const {  // 入っている全部の文字列 t との LCP(s, t) の和  O(|s|)
        long long r = 0;
        for (int v : path(s))
            if (v != 0) r = str_internal::checked_add(r, nd_[v].pass, "Trie::sum_lcp_with");
        return r;
    }
    long long sum_lcp_pairs() const {  // 入っている文字列の全部の組（i < j）の LCP の和（ABC353 E）  O(ノード数)
        long long r = 0;
        for (int v = 1; v < node_count(); ++v) {
            long long x = nd_[v].pass, y = x - 1, t;  // C(pass, 2)。先に偶数の方を 2 で割る
            if (x % 2 == 0) x /= 2;
            else y /= 2;
            if (__builtin_mul_overflow(x, y, &t)) str_internal::fail("Trie::sum_lcp_pairs", "値が long long に収まらない");
            r = str_internal::checked_add(r, t, "Trie::sum_lcp_pairs");
        }
        return r;
    }

    // ---- ノード ----
    int root() const { return 0; }
    int node_count() const { return (int)nd_.size(); }
    template <class Seq>
    int find(const Seq& s0) const {  // s のノード（ノードがなければ -1。数が 0 でもノードがあれば返す）  O(|s|)
        const auto& s = str_internal::seq(s0);
        int v = 0;
        for (std::size_t i = 0; i < s.size() && v >= 0; ++i) v = next_[(std::size_t)v * sg_ + A_.id(s[i])];
        return v;
    }
    template <class Seq>
    std::vector<int> path(const Seq& s0) const {  // 根から s をたどったノードの列（根を含む。途中でなくなればそこまで）  O(|s|)
        const auto& s = str_internal::seq(s0);
        std::vector<int> r = {0};
        int v = 0;
        for (std::size_t i = 0; i < s.size(); ++i) {
            v = next_[(std::size_t)v * sg_ + A_.id(s[i])];
            if (v < 0) break;
            r.push_back(v);
        }
        return r;
    }
    template <class C>
    int child(int v, C c) const {  // 子（文字 c。なければ -1）  O(1)
        check_node(v, "Trie::child");
        return next_[(std::size_t)v * sg_ + A_.id(c)];
    }
    int child_by_id(int v, int id) const {  // 子（文字の番号 id。なければ -1）  O(1)
        check_node(v, "Trie::child_by_id");
        if (id < 0 || id >= sg_) str_internal::fail("Trie::child_by_id", "番号は 0〜σ-1");
        return next_[(std::size_t)v * sg_ + id];
    }
    int parent(int v) const { return check_node(v, "Trie::parent"), nd_[v].par; }  // 親（根は -1）
    int depth(int v) const { return check_node(v, "Trie::depth"), nd_[v].dep; }    // 深さ = 文字列の長さ
    int label(int v) const { return check_node(v, "Trie::label"), nd_[v].lab; }    // 親からの文字の番号（根は -1）
    long long pass_count(int v) const { return check_node(v, "Trie::pass_count"), nd_[v].pass; }  // v を通る（v の文字列で始まる）文字列の数
    long long end_count(int v) const { return check_node(v, "Trie::end_count"), nd_[v].end; }     // v で終わる文字列の数
    std::string string_of(int v) const {  // v の文字列  O(深さ)
        check_node(v, "Trie::string_of");
        std::vector<int> ids;
        for (int u = v; u > 0; u = nd_[u].par) ids.push_back(nd_[u].lab);
        std::reverse(ids.begin(), ids.end());
        return from_ids(ids, A_);
    }
    std::vector<int> preorder() const {  // 全部のノードを辞書順の行きがけ順に（親が子より先。逆順にすれば子が先。数が 0 のノードも）  O(ノード数 × σ)
        std::vector<int> r, st = {0};
        r.reserve(nd_.size());
        while (!st.empty()) {
            const int v = st.back();
            st.pop_back();
            r.push_back(v);
            for (int c = sg_ - 1; c >= 0; --c) {
                const int u = next_[(std::size_t)v * sg_ + c];
                if (u >= 0) st.push_back(u);
            }
        }
        return r;
    }
    const Alphabet& alphabet() const { return A_; }

   private:
    int new_node(int p, int c) {
        const int v = (int)nd_.size();
        next_.insert(next_.end(), (std::size_t)sg_, -1);
        nd_.push_back({0, 0, p, p < 0 ? 0 : nd_[p].dep + 1, c});
        return v;
    }
    void check_node(int v, const char* func) const {
        if (v < 0 || v >= node_count()) str_internal::fail(func, "ノード " + std::to_string(v) + " は 0〜node_count()-1 の外");
    }
    Alphabet A_;
    int sg_;
    struct Node {
        long long pass, end;  // 通る数・終わる数
        int par, dep, lab;    // 親・深さ・親からの文字の番号
    };
    std::vector<int> next_;  // 子（ノード数 × σ）
    std::vector<Node> nd_;
};
#endif  // GG_TRIE
