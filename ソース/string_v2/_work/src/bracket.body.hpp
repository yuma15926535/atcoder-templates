#ifndef GG_BRACKET
#define GG_BRACKET
// ---- 括弧列（文字列第二世代）  判定・対応・深さ・足りない括弧・つなげる順・区間の判定（書き換えも）・? を決める方法の数・全列挙 ----
// 括弧の種類は文字列 pairs で渡す（既定 "()"。"()[]{}<>" なら 4 種類。偶数番目が開き、次が対応する閉じ）
// 括弧以外の文字は無視する（中身として扱う。"a(b)c" は正しい括弧列）
// 対応は「隣り合う開き・閉じ（同じ種類）を消すことを繰り返す」ときに一緒に消える組。種類の違う閉じ括弧があると、それより前の開き括弧は後ろと対応しない
namespace str_internal {
struct BracketKinds {
    std::array<int, 256> kind{};  // 0 = 括弧以外、k + 1 = 種類 k の開き、-(k + 1) = 種類 k の閉じ
    int count = 0;
    explicit BracketKinds(const std::string& pairs, const char* func) {
        if (pairs.empty() || pairs.size() % 2 != 0) fail(func, "pairs は \"()\" や \"()[]\" のように開き・閉じを並べた偶数の長さの文字列");
        for (std::size_t i = 0; i < pairs.size(); ++i) {
            int& k = kind[(unsigned char)pairs[i]];
            if (k != 0) fail(func, std::string("pairs に同じ文字が 2 回ある: '") + pairs[i] + "'");
            k = (i % 2 == 0) ? (int)(i / 2) + 1 : -(int)(i / 2) - 1;
        }
        count = (int)pairs.size() / 2;
    }
    template <class C>
    int of(C c) const {
        const long long v = code_of(c);
        return (v >= 0 && v < 256) ? kind[(std::size_t)v] : 0;
    }
};
inline void need_single(const std::string& pairs, const char* func) {
    if (pairs.size() != 2) fail(func, "1 種類の括弧（pairs の長さ 2）だけで使う");
}
// 対応する位置（なければ -1）と、各境目 i（0〜N）の「状態の番号」。s[l, r) が正しい括弧列 ⇔ id[l] == id[r]
// 状態の番号：開きを積むたびに新しい番号、対応する閉じで 1 つ前に戻る、対応しない閉じで新しい番号（それより前とは対応しない）
template <class Seq>
void bracket_scan(const Seq& s0, const std::string& pairs, const char* func, std::vector<int>* match, std::vector<int>* ids) {
    const auto& s = seq(s0);
    const BracketKinds K(pairs, func);
    const int n = (int)s.size();
    if (match) match->assign(n, -1);
    if (ids) ids->assign(n + 1, 0);
    std::vector<int> st_pos, st_id;  // 積んだ開きの位置と、積む前の状態の番号
    st_pos.reserve(n), st_id.reserve(n);
    int cur = 0, fresh = 1;
    for (int i = 0; i < n; ++i) {
        const int k = K.of(s[i]);
        if (k > 0) {
            st_pos.push_back(i), st_id.push_back(cur);
            cur = fresh++;
        } else if (k < 0) {
            if (!st_pos.empty() && K.of(s[st_pos.back()]) == -k) {
                if (match) (*match)[i] = st_pos.back(), (*match)[st_pos.back()] = i;
                cur = st_id.back();
                st_pos.pop_back(), st_id.pop_back();
            } else {
                st_pos.clear(), st_id.clear();
                cur = fresh++;
            }
        }
        if (ids) (*ids)[i + 1] = cur;
    }
}
}  // namespace str_internal

template <class Seq>
bool is_balanced(const Seq& s, const std::string& pairs = "()") {  // 正しい括弧列か（括弧以外は無視）  O(N)
    std::vector<int> ids;
    str_internal::bracket_scan(s, pairs, "is_balanced", nullptr, &ids);
    return ids.front() == ids.back();
}
template <class Seq>
std::vector<int> match_brackets(const Seq& s, const std::string& pairs = "()") {  // 対応する括弧の位置（対応しない括弧・括弧以外は -1）  O(N)
    std::vector<int> m;
    str_internal::bracket_scan(s, pairs, "match_brackets", &m, nullptr);
    return m;
}
template <class Seq>
std::vector<int> bracket_depths(const Seq& s, const std::string& pairs = "()") {  // r[i] = i を真に囲む（対応した）括弧の組の数  O(N)
    const auto m = match_brackets(s, pairs);
    const int n = (int)m.size();
    std::vector<int> d(n);
    int cur = 0;
    for (int i = 0; i < n; ++i) {
        if (m[i] > i) d[i] = cur++;
        else if (m[i] >= 0) d[i] = --cur;
        else d[i] = cur;
    }
    return d;
}
template <class Seq>
std::vector<int> bracket_parents(const Seq& s, const std::string& pairs = "()") {  // r[i] = i を囲む一番内側の組の開き括弧の位置（なければ -1。入れ子を木に）  O(N)
    const auto m = match_brackets(s, pairs);
    const int n = (int)m.size();
    std::vector<int> par(n, -1), st;
    for (int i = 0; i < n; ++i) {
        if (m[i] >= 0 && m[i] < i) st.pop_back();
        par[i] = st.empty() ? -1 : st.back();
        if (m[i] > i) st.push_back(i);
    }
    return par;
}
template <class Seq>
str_internal::seq_t<Seq> erase_matched_pairs(const Seq& s0, const std::string& pairs = "()") {  // 対応する括弧の組を中身ごと消した残り（ABC307 D）  O(N)
    const auto& s = str_internal::seq(s0);
    const auto m = match_brackets(s, pairs);
    str_internal::seq_t<Seq> r;
    for (int i = 0, n = (int)s.size(); i < n;) {
        if (m[i] > i) i = m[i] + 1;  // 組の中身ごと飛ばす
        else r.push_back(s[i++]);
    }
    return r;
}

// ---- 1 種類の括弧：対応しない括弧の数（連結も計算できる） ----
struct BracketBalance {
    long long need_open = 0;   // 前に足す必要のある開き括弧の数（対応しない閉じ括弧の数）
    long long need_close = 0;  // 後ろに足す必要のある閉じ括弧の数（対応しない開き括弧の数）
    bool balanced() const { return need_open == 0 && need_close == 0; }
    long long net() const { return need_close - need_open; }  // 開きの数 - 閉じの数
    friend BracketBalance operator+(const BracketBalance& a, const BracketBalance& b) {  // a の後ろに b をつなげたもの
        const long long m = std::min(a.need_close, b.need_open);
        return {a.need_open + b.need_open - m, a.need_close + b.need_close - m};
    }
    BracketBalance& operator+=(const BracketBalance& b) { return *this = *this + b; }
    friend bool operator==(const BracketBalance& a, const BracketBalance& b) { return a.need_open == b.need_open && a.need_close == b.need_close; }
    friend bool operator!=(const BracketBalance& a, const BracketBalance& b) { return !(a == b); }
    friend std::ostream& operator<<(std::ostream& os, const BracketBalance& b) { return os << "{" << b.need_open << ", " << b.need_close << "}"; }  // デバッグ用
};
template <class Seq>
BracketBalance balance_of(const Seq& s0, const std::string& pairs = "()") {  // 1 種類の括弧。括弧以外は無視  O(N)
    const auto& s = str_internal::seq(s0);
    str_internal::need_single(pairs, "balance_of");
    const str_internal::BracketKinds K(pairs, "balance_of");
    BracketBalance b;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const int k = K.of(s[i]);
        if (k > 0) ++b.need_close;
        else if (k < 0) {
            if (b.need_close > 0) --b.need_close;
            else ++b.need_open;
        }
    }
    return b;
}
inline BracketBalance balance_of(char c, const std::string& pairs = "()") { return balance_of(std::string(1, c), pairs); }  // 1 文字

struct BracketOrder {
    bool ok = false;         // 正しい括弧列にできるか
    std::vector<int> order;  // つなげる順（vs の添字。ok でなければ空）
};
template <class Seq>
BracketOrder bracket_concat_order(const std::vector<Seq>& vs, const std::string& pairs = "()") {  // 全部を並べ替えてつなげ、正しい括弧列にする順（ABC167 F）  O(Σ|vs[i]| + K log K)
    const int k = (int)vs.size();
    std::vector<BracketBalance> b(k);
    for (int i = 0; i < k; ++i) b[i] = balance_of(vs[i], pairs);
    std::vector<int> up, down;  // 開きが多い（net ≥ 0）もの / 閉じが多いもの
    for (int i = 0; i < k; ++i) (b[i].net() >= 0 ? up : down).push_back(i);
    std::sort(up.begin(), up.end(), [&](int x, int y) { return b[x].need_open != b[y].need_open ? b[x].need_open < b[y].need_open : x < y; });  // 必要な開きが少ない順
    std::sort(down.begin(), down.end(), [&](int x, int y) { return b[x].need_close != b[y].need_close ? b[x].need_close > b[y].need_close : x < y; });  // 残す開きが多い順
    BracketOrder r;
    BracketBalance cur;
    for (int i : up) cur += b[i], r.order.push_back(i);
    for (int i : down) cur += b[i], r.order.push_back(i);
    r.ok = cur.balanced();
    if (!r.ok) r.order.clear();
    return r;
}
template <class Seq>
long long min_insertions_to_balance(const Seq& s, const std::string& pairs = "()") {  // 正しい括弧列にするのに足す括弧の最小の数（1 種類）  O(N)
    const BracketBalance b = balance_of(s, pairs);
    return b.need_open + b.need_close;
}
template <class Seq>
str_internal::seq_t<Seq> fix_brackets(const Seq& s0, const std::string& pairs = "()") {  // 最小の数の括弧を足して正しくしたもの（前に開き、後ろに閉じを足す。1 種類）  O(N)
    const auto& s = str_internal::seq(s0);
    const BracketBalance b = balance_of(s, pairs);
    using E = typename str_internal::seq_t<Seq>::value_type;
    str_internal::seq_t<Seq> r;
    r.reserve(s.size() + (std::size_t)(b.need_open + b.need_close));
    for (long long i = 0; i < b.need_open; ++i) r.push_back((E)(unsigned char)pairs[0]);
    r.insert(r.end(), s.begin(), s.end());
    for (long long i = 0; i < b.need_close; ++i) r.push_back((E)(unsigned char)pairs[1]);
    return r;
}

template <class T, class Seq>
T count_balanced_completions(const Seq& s0, char wild = '?', const std::string& pairs = "()") {  // wild をそれぞれ開き・閉じに決めて正しい括弧列にする方法の数（1 種類。ABC312 D）  O(N^2)
    const auto& s = str_internal::seq(s0);
    const char* func = "count_balanced_completions";
    str_internal::need_single(pairs, func);
    const long long op = (unsigned char)pairs[0], cl = (unsigned char)pairs[1], w = (unsigned char)wild;
    if (w == op || w == cl) str_internal::fail(func, "wild が括弧と同じ文字");
    const int n = (int)s.size();
    std::vector<T> dp(n + 2, T(0)), nd(n + 2, T(0));  // dp[j] = 開き − 閉じ が j の数
    std::vector<char> ov(n + 2, 0), nov(n + 2, 0);
    dp[0] = T(1);
    for (int i = 0; i < n; ++i) {
        const long long c = str_internal::code_of(s[i]);
        if (c != op && c != cl && c != w) str_internal::fail(func, "開き・閉じ・wild 以外の文字がある");
        std::fill(nd.begin(), nd.end(), T(0)), std::fill(nov.begin(), nov.end(), 0);
        for (int j = 0; j <= i; ++j) {
            if (c != cl) {  // 開き
                bool o = nov[j + 1];
                str_internal::add_count(nd[j + 1], o, dp[j], ov[j]);
                nov[j + 1] = o;
            }
            if (c != op && j > 0) {  // 閉じ
                bool o = nov[j - 1];
                str_internal::add_count(nd[j - 1], o, dp[j], ov[j]);
                nov[j - 1] = o;
            }
        }
        std::swap(dp, nd), std::swap(ov, nov);
    }
    if (ov[0]) str_internal::count_overflow(func);
    return dp[0];
}
template <class Seq>
long long count_balanced_substrings(const Seq& s, const std::string& pairs = "()") {  // 正しい括弧列である空でない部分文字列 [l, r) の数（括弧以外だけの部分も数える）  O(N)
    std::vector<int> ids;
    str_internal::bracket_scan(s, pairs, "count_balanced_substrings", nullptr, &ids);
    std::vector<long long> c(ids.size() + 1, 0);
    long long r = 0;
    for (int id : ids) r += c[id]++;
    return r;
}
template <class Seq>
StrSpan longest_balanced_substring(const Seq& s, const std::string& pairs = "()") {  // 正しい括弧列である最長の部分文字列（同じ長さなら左。なければ空 [0, 0)）  O(N)
    std::vector<int> ids;
    str_internal::bracket_scan(s, pairs, "longest_balanced_substring", nullptr, &ids);
    std::vector<int> first(ids.size() + 1, -1);
    StrSpan best{0, 0};
    for (int i = 0; i < (int)ids.size(); ++i) {
        if (first[ids[i]] < 0) first[ids[i]] = i;
        else if (i - first[ids[i]] > best.len()) best = {first[ids[i]], i};
    }
    return best;
}

// ---- 区間 [l, r) が正しい括弧列か O(1)（書き換えなし。何種類でも） ----
struct StaticBracket {
    StaticBracket() = default;
    template <class Seq>
    explicit StaticBracket(const Seq& s, const std::string& pairs = "()") {
        str_internal::bracket_scan(s, pairs, "StaticBracket", &match_, &ids_);
    }
    int size() const { return (int)match_.size(); }
    bool is_balanced(int l, int r) const {  // s[l, r) が正しい括弧列か  O(1)
        str_internal::check_range("StaticBracket::is_balanced", l, r, size());
        return ids_[l] == ids_[r];
    }
    int partner(int i) const {  // 対応する括弧の位置（なければ -1）  O(1)
        if (i < 0 || i >= size()) str_internal::fail("StaticBracket::partner", "位置は 0〜N-1");
        return match_[i];
    }

   private:
    std::vector<int> match_, ids_ = std::vector<int>(1, 0);  // 既定は空の文字列
};

// ---- 1 文字の書き換えと区間の判定（1 種類の括弧。セグメント木） ----
struct BracketSegTree {
    BracketSegTree() = default;
    template <class Seq>
    explicit BracketSegTree(const Seq& s0, const std::string& pairs = "()") : pairs_(pairs) {
        const auto& s = str_internal::seq(s0);
        str_internal::need_single(pairs, "BracketSegTree");
        (void)str_internal::BracketKinds(pairs, "BracketSegTree");
        n_ = (int)s.size();
        sz_ = 1;
        while (sz_ < n_) sz_ <<= 1;
        t_.assign(2 * (std::size_t)sz_, BracketBalance{});
        ch_.resize(n_);
        for (int i = 0; i < n_; ++i) {
            const long long c = str_internal::code_of(s[i]);
            if (c < 0 || c > 255) str_internal::fail("BracketSegTree", "文字（値 0〜255）の列で使う（値 " + std::to_string(c) + "）");
            ch_[i] = (char)(unsigned char)c, t_[sz_ + i] = leaf(ch_[i]);
        }
        for (int v = sz_ - 1; v >= 1; --v) t_[v] = t_[2 * v] + t_[2 * v + 1];
    }
    int size() const { return n_; }
    char get(int i) const {  // i 文字目  O(1)
        check_pos(i, "BracketSegTree::get");
        return ch_[i];
    }
    void set(int i, char c) {  // i 文字目を c に（括弧以外の文字も可）  O(log N)
        check_pos(i, "BracketSegTree::set");
        ch_[i] = c;
        int v = sz_ + i;
        t_[v] = leaf(c);
        for (v >>= 1; v >= 1; v >>= 1) t_[v] = t_[2 * v] + t_[2 * v + 1];
    }
    void flip(int i) {  // 開き ⇔ 閉じ（括弧以外はそのまま）  O(log N)
        check_pos(i, "BracketSegTree::flip");
        if (ch_[i] == pairs_[0]) set(i, pairs_[1]);
        else if (ch_[i] == pairs_[1]) set(i, pairs_[0]);
    }
    BracketBalance balance(int l, int r) const {  // s[l, r) の対応しない括弧の数  O(log N)
        str_internal::check_range("BracketSegTree::balance", l, r, n_);
        BracketBalance L, R;
        for (int a = l + sz_, b = r + sz_; a < b; a >>= 1, b >>= 1) {
            if (a & 1) L = L + t_[a++];
            if (b & 1) R = t_[--b] + R;
        }
        return L + R;
    }
    bool is_balanced(int l, int r) const { return balance(l, r).balanced(); }  // s[l, r) が正しい括弧列か  O(log N)
    bool all_balanced() const { return t_[1].balanced(); }                       // 全体が正しい括弧列か  O(1)
    std::string to_string() const { return std::string(ch_.begin(), ch_.end()); }

   private:
    BracketBalance leaf(char c) const {
        if (c == pairs_[0]) return {0, 1};
        if (c == pairs_[1]) return {1, 0};
        return {0, 0};
    }
    void check_pos(int i, const char* func) const {
        if (i < 0 || i >= n_) str_internal::fail(func, "位置 " + std::to_string(i) + " は 0〜" + std::to_string(n_ - 1) + " の外");
    }
    std::string pairs_ = "()";
    int n_ = 0, sz_ = 1;
    std::vector<BracketBalance> t_ = std::vector<BracketBalance>(2);  // 既定は空の文字列
    std::vector<char> ch_;
};

// ---- 長さ n の正しい括弧列を辞書順（開き < 閉じ）に全部（典型90 002） ----
template <class F>
void for_each_balanced(int n, F f, const std::string& pairs = "()") {  // f(const string&) を辞書順に呼ぶ（n が奇数なら呼ばない）  O(個数 × n)
    str_internal::need_single(pairs, "for_each_balanced");
    if (n < 0 || n % 2 != 0) return;
    const char op = pairs[0], cl = pairs[1];
    std::string s(n / 2, op);
    s += std::string(n / 2, cl);
    std::vector<int> pre(n + 1);  // pre[i] = s[0, i) の（開き - 閉じ）
    while (true) {
        f((const std::string&)s);
        for (int i = 0; i < n; ++i) pre[i + 1] = pre[i] + (s[i] == op ? 1 : -1);
        int i = n - 1;
        while (i >= 0 && !(s[i] == op && pre[i] >= 1)) --i;  // 閉じに変えられる一番右の開き
        if (i < 0) return;
        s[i] = cl;
        const int bal = pre[i] - 1, rest = n - i - 1, opens = (rest - bal) / 2;
        for (int j = 0; j < rest; ++j) s[i + 1 + j] = j < opens ? op : cl;
    }
}
#endif  // GG_BRACKET
