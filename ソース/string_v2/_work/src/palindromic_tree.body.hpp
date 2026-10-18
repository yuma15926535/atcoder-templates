#ifndef GG_PALINDROMIC_TREE
#define GG_PALINDROMIC_TREE
// ---- 回文木（eertree。文字列第二世代） ----
// 異なる回文は N 個以下で、それぞれが 1 つのノード。ノード 0 は長さ -1 の根、ノード 1 は空文字列、2 以降が回文（長さ 1 以上）
// add(c) で 1 文字ずつ足せる（オンライン）。遷移は ノード数 × σ 個の int（N = 10^6・σ = 26 で最大約 100MB。2×10^8 を超えると止める。
// 値の種類が多い数の列には向かない）
struct PalindromicTree {
    explicit PalindromicTree(const Alphabet& A = Alphabet::lower()) : A_(A), sg_(A.size()) {
        new_node(-1, 0, 0), new_node(0, 0, 0);
        suf_.push_back(1);
    }
    template <class Seq>
    explicit PalindromicTree(const Seq& s, const Alphabet& A = Alphabet::lower()) : PalindromicTree(A) {
        const auto& t = str_internal::seq(s);
        const std::size_t m = t.size() + 2;  // ノードは N + 2 個以下。先に取って、伸ばすときの写しをなくす
        for (auto* v : {&len_, &link_, &par_, &first_end_, &diff_, &slink_, &suf_}) v->reserve(m);
        s_.reserve(t.size()), cnt_.reserve(m), go_.reserve(std::min<std::size_t>(m * (std::size_t)sg_, 200000000));
        for (std::size_t i = 0; i < t.size(); ++i) add(t[i]);
    }
    template <class C>
    int add(C ch) {  // 1 文字足し、今の文字列の最長の回文の接尾辞のノードを返す  O(1)（ならし。新しい回文ができたときは O(σ)）
        const int c = A_.id(ch), pos = (int)s_.size();
        s_.push_back(c);
        int cur = last_;
        while (!fits(cur, pos, c)) cur = link_[cur];
        int v = go_[(std::size_t)cur * sg_ + c];
        if (v < 0) {
            int w = link_[cur];
            while (cur != 0 && !fits(w, pos, c)) w = link_[w];
            const int lk = len_[cur] == -1 ? 1 : go_[(std::size_t)w * sg_ + c];  // 長さ 1 なら空文字列へ
            v = new_node(len_[cur] + 2, lk, pos);
            par_[v] = len_[cur] >= 1 ? cur : 1;
            go_[(std::size_t)cur * sg_ + c] = v;
        }
        last_ = v, ++cnt_[v], dirty_ = true;
        suf_.push_back(v);
        return v;
    }

    // ---- 回文の数え上げ ----
    int count_distinct() const { return size() - 2; }  // 異なる回文（空でない）の数
    long long count_total() const {                  // 回文である部分文字列（位置が違えば別）の数  O(ノード数)
        long long r = 0;
        for (long long x : occurrences()) r += x;
        return r;
    }
    const std::vector<long long>& occurrences() const {  // occ[v] = ノード v の回文の出現回数（0・1 は 0。次の add までそのまま使える）  O(ノード数)（足した後の初回）
        if (!dirty_) return occ_;
        occ_.assign(cnt_.begin(), cnt_.end());
        for (int v = size() - 1; v >= 2; --v) occ_[link_[v]] += occ_[v];  // 作った順の逆（リンクの先は先に作られている）
        occ_[0] = occ_[1] = 0;
        dirty_ = false;
        return occ_;
    }
    std::vector<StrSpan> palindromes() const {  // 異なる回文それぞれの最初の出現（ノード 2, 3, … の順）
        std::vector<StrSpan> r;
        for (int v = 2; v < size(); ++v) r.push_back({first_end_[v] - len_[v] + 1, first_end_[v] + 1});
        return r;
    }

    // ---- 回文への分割（N log N） ----
    int min_partition() const {  // 回文に分ける最小の個数  O(N log N)
        std::vector<int> from;
        return partition_dp(from).back();
    }
    std::vector<int> min_partition_prefix() const {  // r[i] = s[0, i) を回文に分ける最小の個数（長さ N + 1）  O(N log N)
        std::vector<int> from;
        return partition_dp(from);
    }
    std::vector<StrSpan> min_partition_spans() const {  // その分け方（左から）  O(N log N)
        std::vector<int> from;
        partition_dp(from);
        std::vector<StrSpan> r;
        for (int i = length(); i > 0; i = from[i]) r.push_back({from[i], i});
        std::reverse(r.begin(), r.end());
        return r;
    }
    template <class T>
    T count_partitions() const {  // 回文に分ける方法の数（T は ModInt998・long long（答えが溢れたら止める）など）  O(N log N)
        std::vector<char> over;
        auto dp = count_dp<T>(over);
        if (over.back()) str_internal::count_overflow("PalindromicTree::count_partitions");
        return dp.back();
    }
    template <class T>
    std::vector<T> count_partitions_prefix() const {  // r[i] = s[0, i) を回文に分ける方法の数（長さ N + 1）  O(N log N)
        std::vector<char> over;
        auto dp = count_dp<T>(over);
        for (char o : over)
            if (o) str_internal::count_overflow("PalindromicTree::count_partitions_prefix");
        return dp;
    }

    // ---- ノード ----
    int size() const { return (int)len_.size(); }  // ノードの数（根 2 つを含む）
    int length() const { return (int)s_.size(); }  // 足した文字の数
    int len(int v) const { return check(v, "PalindromicTree::len"), len_[v]; }
    int link(int v) const { return check(v, "PalindromicTree::link"), link_[v]; }      // 最長の真の回文接尾辞のノード（長さ 1 なら 1 = 空。根は 0）
    int parent(int v) const { return check(v, "PalindromicTree::parent"), par_[v]; }   // 両端を 1 文字ずつ除いた回文のノード（長さ 2 以下なら 1 = 空。根は -1）
    template <class C>
    int next(int v, C c) const {  // 両側に c を足した回文のノード（なければ -1）  O(1)
        check(v, "PalindromicTree::next");
        return go_[(std::size_t)v * sg_ + A_.id(c)];
    }
    int node_of_prefix(int i) const {  // s[0, i) の最長の回文の接尾辞のノード（i = 0 なら 1）  O(1)
        if (i < 0 || i > length()) str_internal::fail("PalindromicTree::node_of_prefix", "i は 0〜N");
        return suf_[i];
    }
    int first_end(int v) const { return check(v, "PalindromicTree::first_end"), first_end_[v]; }  // 最初の出現の最後の文字の位置（根は -1）
    StrSpan span(int v) const {                                                                    // 最初の出現の位置
        check(v, "PalindromicTree::span");
        if (v < 2) return {0, 0};
        return {first_end_[v] - len_[v] + 1, first_end_[v] + 1};
    }
    const Alphabet& alphabet() const { return A_; }

   private:
    std::vector<int> partition_dp(std::vector<int>& from) const {  // 最小の分割（直列リンクでまとめて O(N log N)）
        const int n = (int)s_.size();
        const int INFV = std::numeric_limits<int>::max() / 2;
        std::vector<int> dp(n + 1, INFV);
        std::vector<std::pair<int, int>> ser(size(), {INFV, -1});  // (値, 最後の回文の始まり)
        from.assign(n + 1, -1);
        dp[0] = 0;
        for (int i = 1; i <= n; ++i) {
            for (int u = suf_[i]; len_[u] > 0; u = slink_[u]) {
                const int j = i - (len_[slink_[u]] + diff_[u]);
                ser[u] = {dp[j], j};
                if (diff_[u] == diff_[link_[u]] && ser[link_[u]].first < ser[u].first) ser[u] = ser[link_[u]];
                if (ser[u].first + 1 < dp[i]) dp[i] = ser[u].first + 1, from[i] = ser[u].second;
            }
        }
        return dp;
    }
    template <class T>
    std::vector<T> count_dp(std::vector<char>& over) const {  // 分割の数（over[i] = dp[i] が溢れたか）
        const int n = (int)s_.size();
        std::vector<T> dp(n + 1, T(0)), ser(size(), T(0));
        std::vector<char> sover(size(), 0);
        over.assign(n + 1, 0);
        dp[0] = T(1);
        for (int i = 1; i <= n; ++i) {
            bool o = false;
            for (int u = suf_[i]; len_[u] > 0; u = slink_[u]) {
                const int j = i - (len_[slink_[u]] + diff_[u]);
                ser[u] = dp[j], sover[u] = over[j];
                if (diff_[u] == diff_[link_[u]]) {
                    bool so = sover[u];
                    str_internal::add_count(ser[u], so, ser[link_[u]], sover[link_[u]]);
                    sover[u] = so;
                }
                str_internal::add_count(dp[i], o, ser[u], sover[u]);
            }
            over[i] = o;
        }
        return dp;
    }
    bool fits(int v, int pos, int c) const {  // ノード v の回文の前の文字が c か（s[pos - 1 - len] == c）
        const int k = pos - 1 - len_[v];
        return k >= 0 && s_[k] == c;
    }
    int new_node(int len, int lk, int fe) {
        if ((long long)(len_.size() + 1) * sg_ > 200000000LL) str_internal::fail("PalindromicTree", "ノード数 × σ が 2×10^8（遷移の表が約 800MB）を超える（σ が大きすぎる）");
        len_.push_back(len), link_.push_back(lk), par_.push_back(-1), cnt_.push_back(0), first_end_.push_back(len <= 0 ? -1 : fe);
        go_.insert(go_.end(), (std::size_t)sg_, -1);
        const int v = (int)len_.size() - 1;
        diff_.push_back(len <= 0 ? 0 : len - len_[lk]);
        slink_.push_back(len <= 0 ? 0 : (diff_[v] == diff_[lk] ? slink_[lk] : lk));
        return v;
    }
    void check(int v, const char* func) const {
        if (v < 0 || v >= size()) str_internal::fail(func, "ノード " + std::to_string(v) + " は 0〜size()-1 の外");
    }
    Alphabet A_;
    int sg_;
    int last_ = 1;
    std::vector<int> s_, go_, len_, link_, par_, first_end_, diff_, slink_, suf_;
    std::vector<long long> cnt_;
    mutable bool dirty_ = true;
    mutable std::vector<long long> occ_;
};
#endif  // GG_PALINDROMIC_TREE
