#ifndef GG_SUFFIX_AUTOMATON
#define GG_SUFFIX_AUTOMATON
// ---- 接尾辞オートマトン（文字列第二世代） ----
// SuffixAutomaton sam(s, A)：s のすべての部分文字列を受け付ける最小の DFA。状態の数は 2N 以下（0 が根 = 空文字列）
// 遷移は 状態数 × σ 個の int（N = 5×10^5・σ = 26 で約 100MB。状態数 × σ が 2×10^8 を超えると止める）。extend(c) で 1 文字ずつ、add_string(t) で複数の文字列もまとめられる
// 位置は足した文字の通し番号（1 つの文字列なら s の位置そのもの）
// 部分文字列の問い合わせ（contains・count_occurrences・first_occurrence・longest_common_substring）の t は種類の外の文字を含んでもよい（現れない文字として扱う）
// 出現回数（count_occurrences・endpos_count）と k 番目は、extend・add_string の後の最初の呼び出しで数え直す（O(状態数 × σ) まで）。
// 1 文字足すたびに回数を聞くと O(N^2) になる（1 文字ごとでも速いのは contains・first_occurrence・count_distinct_substrings・total_length_distinct）
struct SuffixAutomaton {
    explicit SuffixAutomaton(const Alphabet& A = Alphabet::lower()) : A_(A), sg_(A.size()) { new_state(0, -1, -1); }
    template <class Seq>
    explicit SuffixAutomaton(const Seq& s, const Alphabet& A = Alphabet::lower()) : SuffixAutomaton(A) {
        const auto& t = str_internal::seq(s);
        reserve_states(std::min<std::size_t>(2 * t.size() + 2, 200000000 / (std::size_t)sg_ + 1));  // 上限を超える分は new_state で止める
        strings_ = 1;
        for (std::size_t i = 0; i < t.size(); ++i) extend(t[i]);
    }
    template <class C>
    void extend(C c) {  // 1 文字足す（最後に足した文字列の続きとして）  O(σ)（ならし）
        if (strings_ == 0) strings_ = 1;
        extend_id(A_.id(c));
    }
    template <class Seq>
    void add_string(const Seq& s) {  // 別の文字列として足す（前の文字列とつながらない。このあとの extend はこの文字列の続き）  O(|s| σ)
        const auto& t = str_internal::seq(s);
        last_ = 0, ++strings_;
        for (std::size_t i = 0; i < t.size(); ++i) extend_id(A_.id(t[i]));
    }

    // ---- 部分文字列 ----
    template <class Seq>
    bool contains(const Seq& t) const { return walk(t) >= 0; }  // t が部分文字列か  O(|t|)
    template <class Seq>
    long long count_occurrences(const Seq& t) const {  // t の出現回数（重なりも。複数の文字列なら合計。空なら Σ(|s_i| + 1)）  O(|t|)（足した後の初回は数え直し）
        const int v = walk(t);
        if (v < 0) return 0;
        if (v == 0) return empty_count();  // 空文字列
        ensure_counts();
        return cnt_[v];
    }
    template <class Seq>
    long long first_occurrence(const Seq& t) const {  // t が最初に現れる開始位置（なければ -1。空なら 0）  O(|t|)
        const int v = walk(t);
        if (v < 0) return -1;
        const long long m = (long long)str_internal::seq(t).size();
        return v == 0 ? 0 : first_end_[v] - m + 1;
    }
    long long count_distinct_substrings() const { return distinct_; }  // 異なる部分文字列（空でない）の数  O(1)（足すたびに更新）
    long long total_length_distinct() const { return total_len_; }    // 異なる部分文字列の長さの合計  O(1)
    std::vector<int> kth_substring_ids(long long k, bool with_dup = false) const {  // 辞書順 k 番目（0 から）の部分文字列（番号の列）。with_dup なら位置が違えば別に数える  O(|答え| σ)（足した後の初回は O(状態数 × σ)）
        ensure_paths(with_dup);
        const std::vector<long long>& P = with_dup ? paths_dup_ : paths_;
        if (k < 0 || k >= P[0]) str_internal::fail("SuffixAutomaton::kth_substring", "k は 0〜（部分文字列の数 - 1）");
        std::vector<int> r;
        int v = 0;
        while (true) {
            int go = -1;
            for (int c = 0; c < sg_ && go < 0; ++c) {
                const int u = next_[(std::size_t)v * sg_ + c];
                if (u < 0) continue;
                const long long here = with_dup ? cnt_[u] : 1;  // c を足した文字列そのものの数
                if (k < here) {
                    r.push_back(c);
                    return r;
                }
                if (k < here + P[u]) k -= here, go = u, r.push_back(c);  // その先へ
                else k -= here + P[u];
            }
            v = go;
        }
    }
    std::string kth_substring(long long k, bool with_dup = false) const { return from_ids(kth_substring_ids(k, with_dup), A_); }  // 文字列で返す
    long long count_substrings_with_dup() const {  // 部分文字列（空でない。位置が違えば別）の数  O(状態数 × σ)（足した後の初回）
        ensure_paths(true);
        return paths_dup_[0];
    }
    template <class Seq>
    StrSpan longest_common_substring(const Seq& t0) const {  // t の部分文字列で、この文字列にも現れる最長のもの（t での位置。同じ長さなら左）  O(|t|)
        const auto& t = str_internal::seq(t0);
        int v = 0, l = 0;
        StrSpan best{0, 0};
        for (int i = 0; i < (int)t.size(); ++i) {
            const int c = A_.contains(t[i]) ? A_.id(t[i]) : -1;
            while (v != 0 && (c < 0 || next_[(std::size_t)v * sg_ + c] < 0)) v = link_[v], l = len_[v];
            if (c >= 0 && next_[(std::size_t)v * sg_ + c] >= 0) v = next_[(std::size_t)v * sg_ + c], ++l;
            if (l > best.len()) best = {i + 1 - l, i + 1};
        }
        return best;
    }

    // ---- 状態（DP 用） ----
    int size() const { return (int)len_.size(); }  // 状態の数
    int len(int v) const { return check(v, "SuffixAutomaton::len"), len_[v]; }      // 状態 v の最長の文字列の長さ
    int link(int v) const { return check(v, "SuffixAutomaton::link"), link_[v]; }   // 接尾辞リンク（根は -1）
    template <class C>
    int next(int v, C c) const {  // 文字 c の遷移（なければ -1。c が種類の外なら止める。番号なら next_id）  O(1)
        check(v, "SuffixAutomaton::next");
        return next_[(std::size_t)v * sg_ + A_.id(c)];
    }
    int next_id(int v, int id) const {
        check(v, "SuffixAutomaton::next_id");
        if (id < 0 || id >= sg_) str_internal::fail("SuffixAutomaton::next_id", "番号は 0〜σ-1");
        return next_[(std::size_t)v * sg_ + id];
    }
    long long endpos_count(int v) const {  // 状態 v の文字列が現れる回数（根は空文字列の数 Σ(|s_i| + 1)）  O(1)（足した後の初回は O(状態数)）
        check(v, "SuffixAutomaton::endpos_count");
        if (v == 0) return empty_count();
        ensure_counts();
        return cnt_[v];
    }
    long long first_end(int v) const { return check(v, "SuffixAutomaton::first_end"), first_end_[v]; }  // 状態 v の文字列が最初に現れる終わりの位置（根は -1）
    std::vector<int> order_by_len() const {  // 状態を len の小さい順に（リンクの先が先。逆にたどれば部分木の和）  O(状態数 + 文字の数)
        std::vector<int> c(chars_ + 2, 0), r(size());
        for (int v = 0; v < size(); ++v) ++c[len_[v]];
        for (std::size_t i = 1; i < c.size(); ++i) c[i] += c[i - 1];
        for (int v = size() - 1; v >= 0; --v) r[--c[len_[v]]] = v;
        return r;
    }
    int last() const { return last_; }  // 最後に足した文字列全体の状態（何も足していなければ 0）
    const Alphabet& alphabet() const { return A_; }

   private:
    int new_state(int len, int link, long long fe) {
        if ((long long)(len_.size() + 1) * sg_ > 200000000LL) str_internal::fail("SuffixAutomaton", "状態数 × σ が 2×10^8（遷移の表が約 800MB）を超える（σ が大きい数の列なら 03_suffix_array）");
        len_.push_back(len), link_.push_back(link), first_end_.push_back(fe), occ_.push_back(0);
        next_.insert(next_.end(), (std::size_t)sg_, -1);
        return (int)len_.size() - 1;
    }
    void reserve_states(std::size_t n) {
        len_.reserve(n), link_.reserve(n), first_end_.reserve(n), occ_.reserve(n), next_.reserve(n * (std::size_t)sg_);
    }
    int clone_of(int q, int len) {
        const int cl = new_state(len, link_[q], first_end_[q]);
        std::copy(next_.begin() + (std::size_t)q * sg_, next_.begin() + (std::size_t)(q + 1) * sg_, next_.begin() + (std::size_t)cl * sg_);
        link_[q] = cl;
        return cl;
    }
    void extend_id(int c) {
        dirty_ = true;
        const long long pos = chars_++;
        int p = last_;
        const int ex = next_[(std::size_t)p * sg_ + c];
        if (ex >= 0) {  // 複数の文字列：もう同じ遷移がある
            if (len_[ex] == len_[p] + 1) {
                last_ = ex;
            } else {
                const int cl = clone_of(ex, len_[p] + 1);
                for (; p >= 0 && next_[(std::size_t)p * sg_ + c] == ex; p = link_[p]) next_[(std::size_t)p * sg_ + c] = cl;
                last_ = cl;
            }
            ++occ_[last_];
            return;
        }
        const int cur = new_state(len_[last_] + 1, 0, pos);
        occ_[cur] = 1;
        for (; p >= 0 && next_[(std::size_t)p * sg_ + c] < 0; p = link_[p]) next_[(std::size_t)p * sg_ + c] = cur;
        if (p >= 0) {
            const int q = next_[(std::size_t)p * sg_ + c];
            if (len_[q] == len_[p] + 1) {
                link_[cur] = q;
            } else {
                const int cl = clone_of(q, len_[p] + 1);
                for (; p >= 0 && next_[(std::size_t)p * sg_ + c] == q; p = link_[p]) next_[(std::size_t)p * sg_ + c] = cl;
                link_[cur] = cl;
            }
        }
        last_ = cur;
        const long long a = len_[cur], b = len_[link_[cur]];  // 新しい部分文字列は長さ b+1〜a（複製では合計は変わらない）
        distinct_ += a - b, total_len_ += a * (a + 1) / 2 - b * (b + 1) / 2;
    }
    long long empty_count() const { return chars_ + std::max(strings_, 1); }
    template <class Seq>
    int walk(const Seq& t0) const {  // t をたどった状態（なければ -1）
        const auto& t = str_internal::seq(t0);
        int v = 0;
        for (std::size_t i = 0; i < t.size() && v >= 0; ++i) {
            if (!A_.contains(t[i])) return -1;
            v = next_[(std::size_t)v * sg_ + A_.id(t[i])];
        }
        return v;
    }
    void check(int v, const char* func) const {
        if (v < 0 || v >= size()) str_internal::fail(func, "状態 " + std::to_string(v) + " は 0〜size()-1 の外");
    }
    void ensure_counts() const {
        if (!dirty_ && !cnt_.empty()) return;
        dirty_ = false, paths_.clear(), paths_dup_.clear();
        cnt_.assign(occ_.begin(), occ_.end());
        const auto ord = order_by_len();
        for (int k = size() - 1; k > 0; --k) cnt_[link_[ord[k]]] += cnt_[ord[k]];
    }
    void ensure_paths(bool dup) const {
        ensure_counts();
        std::vector<long long>& P = dup ? paths_dup_ : paths_;
        if (!P.empty()) return;
        // P[v] = v から 1 文字以上進んでできる文字列の数（重みつき）
        P.assign(size(), 0);
        const auto ord = order_by_len();
        for (int k = size() - 1; k >= 0; --k) {
            const int v = ord[k];
            long long s = 0;
            for (int c = 0; c < sg_; ++c) {
                const int u = next_[(std::size_t)v * sg_ + c];
                if (u >= 0) s += (dup ? cnt_[u] : 1) + P[u];
            }
            P[v] = s;
        }
    }
    Alphabet A_;
    int sg_;
    int last_ = 0, strings_ = 0;
    long long chars_ = 0, distinct_ = 0, total_len_ = 0;
    std::vector<int> next_, len_, link_;
    std::vector<long long> first_end_, occ_;
    mutable bool dirty_ = true;
    mutable std::vector<long long> cnt_, paths_, paths_dup_;
};
#endif  // GG_SUFFIX_AUTOMATON
