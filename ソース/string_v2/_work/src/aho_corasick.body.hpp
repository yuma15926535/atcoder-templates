#ifndef GG_AHO_CORASICK
#define GG_AHO_CORASICK
// ---- Aho–Corasick（多数のパターンの照合。文字列第二世代） ----
// AhoCorasick ac(A)：add(p) でパターンを足し（番号 0, 1, … を返す。同じ文字列を 2 回足してもよい）、build() してから使う
// 状態 = Trie のノード（0 が根）。next(state, c) は O(1)（遷移の表は ノード数 × σ 個の int。2×10^8 を超えると止める）。build の後に add はできない
// 本文・パターンの文字は Alphabet の範囲で（範囲外なら止める）
struct AhoCorasick {
    explicit AhoCorasick(const Alphabet& A = Alphabet::lower()) : A_(A), sg_(A.size()) { new_node(-1, 0); }

    template <class Seq>
    int add(const Seq& p0) {  // パターン p を足して番号を返す（空のパターンは止める）  O(|p|)
        if (built_) str_internal::fail("AhoCorasick::add", "build() の後には足せない");
        const auto& p = str_internal::seq(p0);
        if (p.size() == 0) str_internal::fail("AhoCorasick::add", "空のパターンは使えない");
        int v = 0;
        for (std::size_t i = 0; i < p.size(); ++i) {
            const int c = A_.id(p[i]);
            int u = go_[(std::size_t)v * sg_ + c];
            if (u < 0) u = new_node(v, dep_[v] + 1), go_[(std::size_t)v * sg_ + c] = u;
            v = u;
        }
        node_of_.push_back(v), len_.push_back((int)p.size());
        return (int)node_of_.size() - 1;
    }
    void build() {  // 失敗のリンクと遷移の表を作る  O(ノード数 × σ)
        if (built_) return;
        built_ = true;
        const int n = size();
        link_.assign(n, 0), out_.assign(n, -1), cnt_.assign(n, 0), longest_.assign(n, 0), order_.clear();
        // ノードごとのパターンの番号（CSR）
        ids_start_.assign(n + 1, 0);
        for (int v : node_of_) ++ids_start_[v + 1], ++cnt_[v];
        for (int v = 0; v < n; ++v) ids_start_[v + 1] += ids_start_[v];
        ids_.assign(node_of_.size(), 0);
        std::vector<int> fill(ids_start_.begin(), ids_start_.end() - 1);
        for (int id = 0; id < (int)node_of_.size(); ++id) ids_[fill[node_of_[id]]++] = id;
        // 幅優先
        order_.reserve(n);
        order_.push_back(0);
        for (int c = 0; c < sg_; ++c) {
            int& u = go_[c];
            if (u < 0) u = 0;
            else link_[u] = 0, order_.push_back(u);
        }
        for (std::size_t h = 1; h < order_.size(); ++h) {
            const int v = order_[h];
            for (int c = 0; c < sg_; ++c) {
                int& u = go_[(std::size_t)v * sg_ + c];
                const int f = go_[(std::size_t)link_[v] * sg_ + c];
                if (u < 0) u = f;
                else link_[u] = f, order_.push_back(u);
            }
        }
        for (int v : order_) {  // 根から順に：リンクの先までの数・最長・出力のリンク
            if (v == 0) continue;
            const int f = link_[v];
            out_[v] = (ids_start_[f + 1] > ids_start_[f]) ? f : out_[f];
            longest_[v] = (ids_start_[v + 1] > ids_start_[v]) ? dep_[v] : longest_[f];
            cnt_[v] += cnt_[f];
        }
    }

    // ---- 状態（DP 用） ----
    int size() const { return (int)dep_.size(); }                // 状態の数
    int pattern_count() const { return (int)node_of_.size(); }   // パターンの数
    template <class C>
    int next(int state, C c) const {  // 状態 state で文字 c を読んだ次の状態  O(1)
        check_built("AhoCorasick::next"), check_state(state, "AhoCorasick::next");
        return go_[(std::size_t)state * sg_ + A_.id(c)];
    }
    int next_id(int state, int id) const {  // 文字の番号 id で  O(1)
        check_built("AhoCorasick::next_id"), check_state(state, "AhoCorasick::next_id");
        if (id < 0 || id >= sg_) str_internal::fail("AhoCorasick::next_id", "番号は 0〜σ-1");
        return go_[(std::size_t)state * sg_ + id];
    }
    long long match_count(int state) const { return check_built("AhoCorasick::match_count"), check_state(state, "AhoCorasick::match_count"), cnt_[state]; }  // ここで終わるパターンの数（同じ文字列も別に数える）
    int longest_match(int state) const { return check_built("AhoCorasick::longest_match"), check_state(state, "AhoCorasick::longest_match"), longest_[state]; }  // ここで終わる最長のパターンの長さ（なければ 0）
    std::vector<int> match_ids(int state) const {  // ここで終わるパターンの番号（長い順）  O(個数)
        check_built("AhoCorasick::match_ids"), check_state(state, "AhoCorasick::match_ids");
        std::vector<int> r;
        for (int v = (ids_start_[state + 1] > ids_start_[state]) ? state : out_[state]; v > 0; v = out_[v])
            for (int k = ids_start_[v]; k < ids_start_[v + 1]; ++k) r.push_back(ids_[k]);
        return r;
    }
    int depth(int state) const { return check_state(state, "AhoCorasick::depth"), dep_[state]; }                                   // 根からの深さ（状態の文字列の長さ）
    int link(int state) const { return check_built("AhoCorasick::link"), check_state(state, "AhoCorasick::link"), link_[state]; }  // 失敗のリンク（根は 0）
    int node_of(int id) const { return check_id(id, "AhoCorasick::node_of"), node_of_[id]; }  // パターン id が終わる状態
    int pattern_length(int id) const { return check_id(id, "AhoCorasick::pattern_length"), len_[id]; }
    const std::vector<int>& bfs_order() const { return check_built("AhoCorasick::bfs_order"), order_; }  // 根から幅優先の順（リンクの先が先に来る）

    // ---- 文字列の中の出現 ----
    template <class Seq>
    std::vector<long long> count_each(const Seq& t0) const {  // パターンごとの出現回数（ABC362 G）  O(|t| + ノード数)
        check_built("AhoCorasick::count_each");
        const auto& t = str_internal::seq(t0);
        std::vector<long long> vis(size(), 0);
        int v = 0;
        for (std::size_t i = 0; i < t.size(); ++i) v = go_[(std::size_t)v * sg_ + A_.id(t[i])], ++vis[v];
        for (int k = size() - 1; k > 0; --k) vis[link_[order_[k]]] += vis[order_[k]];
        std::vector<long long> r(pattern_count());
        for (int id = 0; id < pattern_count(); ++id) r[id] = vis[node_of_[id]];
        return r;
    }
    template <class Seq>
    long long count_total(const Seq& t0) const {  // すべてのパターンの出現回数の合計  O(|t|)
        check_built("AhoCorasick::count_total");
        const auto& t = str_internal::seq(t0);
        long long s = 0;
        int v = 0;
        for (std::size_t i = 0; i < t.size(); ++i) v = go_[(std::size_t)v * sg_ + A_.id(t[i])], s += cnt_[v];
        return s;
    }
    template <class Seq>
    std::vector<int> first_occurrence_each(const Seq& t0) const {  // パターンごとの最初の出現の開始位置（なければ -1）  O(|t| + ノード数)
        check_built("AhoCorasick::first_occurrence_each");
        const auto& t = str_internal::seq(t0);
        const int INFPOS = std::numeric_limits<int>::max();
        std::vector<int> first(size(), INFPOS);
        int v = 0;
        for (std::size_t i = 0; i < t.size(); ++i) {
            v = go_[(std::size_t)v * sg_ + A_.id(t[i])];
            if (first[v] == INFPOS) first[v] = (int)i;
        }
        for (int k = size() - 1; k > 0; --k) first[link_[order_[k]]] = std::min(first[link_[order_[k]]], first[order_[k]]);
        std::vector<int> r(pattern_count());
        for (int id = 0; id < pattern_count(); ++id) r[id] = first[node_of_[id]] == INFPOS ? -1 : first[node_of_[id]] - len_[id] + 1;
        return r;
    }
    template <class Seq, class F>
    void for_each_match(const Seq& t0, F f) const {  // すべての出現で f(開始位置, パターンの番号)（終わる位置の順）  O(|t| + 出現の数)
        check_built("AhoCorasick::for_each_match");
        const auto& t = str_internal::seq(t0);
        int v = 0;
        for (std::size_t i = 0; i < t.size(); ++i) {
            v = go_[(std::size_t)v * sg_ + A_.id(t[i])];
            for (int u = (ids_start_[v + 1] > ids_start_[v]) ? v : out_[v]; u > 0; u = out_[u])
                for (int k = ids_start_[u]; k < ids_start_[u + 1]; ++k) f((int)i - len_[ids_[k]] + 1, ids_[k]);
        }
    }
    template <class Seq>
    bool contains_any(const Seq& t0) const {  // どれかのパターンが現れるか  O(|t|)
        check_built("AhoCorasick::contains_any");
        const auto& t = str_internal::seq(t0);
        int v = 0;
        for (std::size_t i = 0; i < t.size(); ++i) {
            v = go_[(std::size_t)v * sg_ + A_.id(t[i])];
            if (cnt_[v] > 0) return true;
        }
        return false;
    }
    std::vector<int> contained_patterns(int id) const {  // パターン id の中に（部分文字列として）現れるほかのパターンの番号（昇順。同じ文字列のパターンも入る）  O(|p| + 個数 log 個数)
        check_built("AhoCorasick::contained_patterns");
        check_id(id, "AhoCorasick::contained_patterns");
        // id の接頭辞の状態それぞれから出力のリンクをたどって、そこで終わるパターンを集める（一度見た状態から先は集めたあとなので止める）
        if (mark_.size() != (std::size_t)size()) mark_.assign(size(), 0);
        ++stamp_;
        std::vector<int> r;
        for (int v : path_nodes(id)) {
            for (int u = (ids_start_[v + 1] > ids_start_[v]) ? v : out_[v]; u > 0 && mark_[u] != stamp_; u = out_[u]) {
                mark_[u] = stamp_;
                for (int k = ids_start_[u]; k < ids_start_[u + 1]; ++k)
                    if (ids_[k] != id) r.push_back(ids_[k]);
            }
        }
        std::sort(r.begin(), r.end());
        return r;
    }
    const Alphabet& alphabet() const { return A_; }

   private:
    int new_node(int parent, int d) {
        if ((long long)(dep_.size() + 1) * sg_ > 200000000LL) str_internal::fail("AhoCorasick", "ノード数 × σ が 2×10^8（遷移の表が約 800MB）を超える（σ が大きい数の列には向かない）");
        go_.insert(go_.end(), (std::size_t)sg_, -1);
        dep_.push_back(d), par_.push_back(parent);
        return (int)dep_.size() - 1;
    }
    std::vector<int> path_nodes(int id) const {  // パターン id の接頭辞の状態（長さ 1〜|p|）
        std::vector<int> r(len_[id]);
        int v = node_of_[id];
        for (int k = len_[id] - 1; k >= 0; --k) r[k] = v, v = par_[v];
        return r;
    }
    void check_built(const char* func) const {
        if (!built_) str_internal::fail(func, "先に build() を呼ぶ");
    }
    void check_state(int s, const char* func) const {
        if (s < 0 || s >= size()) str_internal::fail(func, "状態 " + std::to_string(s) + " は 0〜size()-1 の外");
    }
    void check_id(int id, const char* func) const {
        if (id < 0 || id >= pattern_count()) str_internal::fail(func, "パターンの番号 " + std::to_string(id) + " は 0〜pattern_count()-1 の外");
    }
    Alphabet A_;
    int sg_;
    bool built_ = false;
    std::vector<int> go_, dep_, par_, node_of_, len_, link_, out_, longest_, order_, ids_start_, ids_;
    std::vector<long long> cnt_;
    mutable std::vector<int> mark_;  // contained_patterns の「見た」印
    mutable int stamp_ = 0;
};
#endif  // GG_AHO_CORASICK
